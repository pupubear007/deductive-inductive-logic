"""The full trainable model: encoders -> fused embedding z -> three heads.

1. Gated chain head (deductive gates + defeasible step kernels, calibrated probability):
   level probabilities that satisfy the hard constraints by construction (mechanism.py).
2. Conditional outcome head: E[outcome | infection, inputs] for incidence, DSI, lesion length,
   and an unconditional regression for yield.
3. Generative heads: conditional diffusion over the outcome vector and over 2-D disease maps.
"""

from __future__ import annotations

from typing import Any

import torch
from torch import nn
from torch.nn import functional as F

from ..mechanism import ChainProbs, chain_probabilities
from .diffusion import MapDiffusion, VectorDiffusion
from .encoders import Fusion, ImageEncoder, SeriesEncoder, TabularEncoder

CHAIN_LOGITS = ("apothecia", "external", "flowering", "microclimate", "q_court", "q_inf")


def _logit(p: float) -> float:
    p = min(max(p, 1e-4), 1 - 1e-4)
    return float(torch.logit(torch.tensor(p)))


class SSRModel(nn.Module):
    def __init__(self, meta, cfg: dict[str, Any]):
        super().__init__()
        mc, dc = cfg["model"], cfg["diffusion"]
        dim = mc["dim"]
        self.meta, self.cfg = meta, cfg
        self.mods = list(meta.modalities)
        self.img = nn.ModuleDict({m: ImageEncoder(s["channels"], dim, mc["image_width"]) for m, s in meta.modalities.items()})
        self.series = SeriesEncoder(len(meta.series_vars), dim, mc["series_hidden"])
        self.tab = TabularEncoder(len(meta.num_cols), [len(meta.vocab[c]) for c in meta.cat_cols], dim)
        self.fuse = Fusion(2 + len(self.mods), dim, mc["modality_dropout"])
        self.chain_head = nn.Sequential(nn.Linear(dim, dim), nn.SiLU(), nn.Linear(dim, len(CHAIN_LOGITS)))
        self.cont_head = nn.Sequential(nn.Linear(dim, dim), nn.SiLU(), nn.Linear(dim, max(len(meta.cont_names), 1)))

        # Assay observation channels P(reading = 1 | level): sensitivity / specificity (Def. S10).
        se, sp, learn = [], [], []
        for a in meta.bin_names:
            spec = cfg["assays"][a]
            se.append(_logit(spec.get("se", 1.0)))
            sp.append(_logit(spec.get("sp", 1.0)))
            learn.append(bool(spec.get("learn", False)))
        self.se_logit = nn.Parameter(torch.tensor(se), requires_grad=any(learn))
        self.sp_logit = nn.Parameter(torch.tensor(sp), requires_grad=any(learn))
        self.register_buffer("learn_mask", torch.tensor(learn, dtype=torch.float32))

        nc = len(meta.cont_names)
        self.vdiff = VectorDiffusion(nc, dim, dc["hidden"], dc["steps"]) if dc.get("vector") and nc else None
        self.mdiff = MapDiffusion(meta.map_size, dim, dc["steps"]) if dc.get("map") else None

    # ------------------------------------------------------------------------------------------
    def embed(self, b: dict[str, Any]) -> torch.Tensor:
        parts = [self.series(b["series"], b["series_mask"]), self.tab(b["num"], b["cat"])]
        present = [(b["series_mask"].flatten(1).sum(1) > 0).float(), torch.ones_like(b["do_court"])]
        for m in self.mods:
            parts.append(self.img[m](b["images"][m]))
            present.append(b["present"][m])
        return self.fuse(parts, present)

    def forward(self, b: dict[str, Any]) -> dict[str, Any]:
        z = self.embed(b)
        lg = self.chain_head(z)
        logits = {n: lg[:, i] for i, n in enumerate(CHAIN_LOGITS)}
        chain = chain_probabilities(logits, b["do_court"])
        cont = self.cont_head(z)
        kinds = [self.cfg["assays"][a]["kind"] for a in self.meta.cont_names]
        cont_mean = torch.stack(
            [torch.sigmoid(cont[:, j]) if k == "fraction" else F.softplus(cont[:, j]) for j, k in enumerate(kinds)], -1
        ) if kinds else cont[:, :0]
        return {"z": z, "chain": chain, "cont_mean": cont_mean}

    def assay_channels(self) -> tuple[torch.Tensor, torch.Tensor]:
        """Sensitivity and specificity actually used (fixed ones are detached)."""
        se, sp = torch.sigmoid(self.se_logit), torch.sigmoid(self.sp_logit)
        m = self.learn_mask
        return m * se + (1 - m) * se.detach(), m * sp + (1 - m) * sp.detach()

    def p_positive(self, chain: ChainProbs) -> torch.Tensor:
        """P(assay reads positive) for every binary assay, shape [B, n_bin]."""
        se, sp = self.assay_channels()
        cols = []
        for j, a in enumerate(self.meta.bin_names):
            p = chain.get(self.cfg["assays"][a]["level"])
            cols.append(se[j] * p + (1 - sp[j]) * (1 - p))
        return torch.stack(cols, -1) if cols else chain.infection[:, None][:, :0]

    def expected_outcomes(self, out: dict[str, Any]) -> torch.Tensor:
        """E[reading] = P(level) * E[reading | level] for conditional outcomes."""
        cols = []
        for j, a in enumerate(self.meta.cont_names):
            spec = self.cfg["assays"][a]
            lvl = spec.get("conditional_on") or (spec.get("level") if spec["kind"] == "fraction" else None)
            m = out["cont_mean"][:, j]
            cols.append(out["chain"].get(lvl) * m if lvl else m)
        return torch.stack(cols, -1) if cols else out["cont_mean"]


def compute_losses(model: SSRModel, b: dict[str, Any], out: dict[str, Any], cfg: dict[str, Any],
                   diffusion: bool = False) -> dict[str, torch.Tensor]:
    """Stage 1 trains encoders + chain + outcome heads (diffusion=False). The diffusion heads are
    trained in stage 2 on the frozen embedding (see train.train_diffusion)."""
    w = cfg["loss_weights"]
    losses: dict[str, torch.Tensor] = {}

    # 1. Chain: likelihood of each binary reading through its observation channel.
    if b["m_bin"].numel() and b["m_bin"].sum() > 0:
        p = model.p_positive(out["chain"]).clamp(1e-6, 1 - 1e-6)
        bce = F.binary_cross_entropy(p, b["y_bin"], reduction="none")
        losses["chain"] = (bce * b["m_bin"]).sum() / b["m_bin"].sum()

    # 2. Conditional regression: fit E[y | level] on readings where the level is present (y > 0);
    #    unconditional outcomes (yield) on all readings.
    if b["m_cont"].numel() and b["m_cont"].sum() > 0:
        terms, weights = [], []
        for j, a in enumerate(model.meta.cont_names):
            spec = cfg["assays"][a]
            cond = spec.get("conditional_on") or spec["kind"] == "fraction"
            m = b["m_cont"][:, j] * ((b["y_cont"][:, j] > 0).float() if cond else 1.0)
            if m.sum() > 0:
                terms.append((F.smooth_l1_loss(out["cont_mean"][:, j], b["y_cont"][:, j], reduction="none", beta=0.1) * m).sum())
                weights.append(m.sum())
        if terms:
            losses["regression"] = torch.stack(terms).sum() / torch.stack(weights).sum()

    # 3. Generative heads.
    if diffusion and model.vdiff is not None and b["m_cont"].sum() > 0:
        losses["vector_diffusion"] = model.vdiff.loss(b["y_cont"].clamp(-0.5, 1.5), b["m_cont"], out["z"])
    if diffusion and model.mdiff is not None and b["m_map"].sum() > 0:
        losses["map_diffusion"] = model.mdiff.loss(b["y_map"], b["m_map"], out["z"])

    total = sum(w.get(k, 1.0) * v for k, v in losses.items())
    losses["total"] = total if isinstance(total, torch.Tensor) else out["z"].sum() * 0.0
    return losses
