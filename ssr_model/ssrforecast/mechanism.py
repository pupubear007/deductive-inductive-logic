"""The SSR mechanism: levels, gates, hard constraints and interventions.

This module is the bridge between the paper (Wang, "Deduction and Induction in One Diagram",
Sections 3-8) and the trainable model. Nothing here is learned; it fixes the *shape* that every
learned component must respect.

Levels (entailment order, Def. 6.3; a later stage entails the earlier one)::

    spore  <=  court  <=  infection  <=  incidence>0
    apothecia (local) entails spore (sufficient source), but spore does NOT entail local
    apothecia: ascospores may come from outside the field.

Gates (hard necessary conditions, mechanism condition M1)::

    court     |=  spore AND flowering        (petals present while spores present)
    infection |=  court AND microclimate     (permissive canopy microclimate)

Progressions (defeasible, mechanism condition M2) are the learned step kernels q_k:

    P(court | spore, flowering, x)        = q_court(x)
    P(infection | court, microclimate, x) = q_inf(x)

So, under conditional independence of the gate and the lower level given the inputs x,

    p_spore = 1 - (1 - p_apo) * (1 - p_ext)
    p_court = p_spore * p_flower * q_court
    p_inf   = p_court * p_micro  * q_inf

which makes p_inf <= p_court <= p_spore hold *by construction* (no world violating a hard
constraint ever receives probability), i.e. the model is supported on W_K (Definition S1).

Interventions (do-operator). A growth-chamber or greenhouse inoculation (cut petiole, mycelial
plug) sets the infection court by intervention: do(court = 1). It is not an observation of the
court level, and it bypasses the spore and flowering levels. This is why chamber data inform
q_inf and lesion growth but not the lower steps of the chain.
"""

from __future__ import annotations

from dataclasses import dataclass

import torch

# Binary levels the chain head produces, in causal order.
LEVELS: tuple[str, ...] = ("apothecia", "external", "spore", "flowering", "court", "microclimate", "infection")
# Which of these are gates (necessary side conditions, not levels of the chain).
GATES: tuple[str, ...] = ("flowering", "microclimate")
# Inoculation methods that intervene on the court level.
COURT_INTERVENTIONS: frozenset[str] = frozenset({"cut_petiole", "mycelial_plug", "colonized_petal", "wound"})


@dataclass
class ChainProbs:
    """Probabilities of every level and gate, shape [B] each."""

    apothecia: torch.Tensor
    external: torch.Tensor
    spore: torch.Tensor
    flowering: torch.Tensor
    court: torch.Tensor
    microclimate: torch.Tensor
    infection: torch.Tensor

    def get(self, name: str) -> torch.Tensor:
        return getattr(self, name)

    def stack(self) -> torch.Tensor:
        return torch.stack([self.get(n) for n in LEVELS], dim=-1)


def chain_probabilities(logits: dict[str, torch.Tensor], do_court: torch.Tensor | None = None) -> ChainProbs:
    """Combine raw logits into level probabilities that satisfy the hard constraints.

    ``logits`` must contain: ``apothecia``, ``external`` (source terms), ``flowering``,
    ``microclimate`` (gates) and ``q_court``, ``q_inf`` (defeasible step kernels).
    ``do_court`` is a {0,1} tensor; 1 means the court was set by inoculation.
    """
    s = torch.sigmoid
    p_apo = s(logits["apothecia"])
    p_ext = s(logits["external"])
    p_spore = 1.0 - (1.0 - p_apo) * (1.0 - p_ext)
    p_flower = s(logits["flowering"])
    p_court_nat = p_spore * p_flower * s(logits["q_court"])
    if do_court is not None:
        d = do_court.to(p_court_nat.dtype)
        p_court = d + (1.0 - d) * p_court_nat
    else:
        p_court = p_court_nat
    p_micro = s(logits["microclimate"])
    p_inf = p_court * p_micro * s(logits["q_inf"])
    return ChainProbs(p_apo, p_ext, p_spore, p_flower, p_court, p_micro, p_inf)


def violations(levels: dict[str, torch.Tensor], do_court: torch.Tensor | None = None) -> torch.Tensor:
    """Count hard-constraint violations in a batch of *binary* level assignments.

    Used to audit labels and generated samples. Returns a [B] integer tensor.
    Constraints: court -> spore and flowering (unless do(court)); infection -> court and microclimate;
    apothecia -> spore.
    """
    b = {k: v.bool() for k, v in levels.items()}
    nat = torch.ones_like(b["court"]) if do_court is None else ~do_court.bool()
    v = torch.zeros_like(b["court"], dtype=torch.long)
    v += (b["apothecia"] & ~b["spore"]).long()
    v += (nat & b["court"] & ~(b["spore"] & b["flowering"])).long()
    v += (b["infection"] & ~(b["court"] & b["microclimate"])).long()
    return v


def project_outcomes(y: torch.Tensor, names: list[str], eps: float = 1e-3) -> torch.Tensor:
    """Project generated continuous outcomes onto admissible worlds.

    * fractions are clipped to [0, 1], other outcomes to >= 0;
    * if incidence is (numerically) zero, severity and lesion length must be zero
      (no infection, no disease severity: definitional constraint).
    """
    y = y.clone()
    idx = {n: i for i, n in enumerate(names)}
    for n, i in idx.items():
        y[:, i] = y[:, i].clamp(0.0, 1.0) if n in ("incidence", "dsi") else y[:, i].clamp_min(0.0)
    if "incidence" in idx:
        no_inf = y[:, idx["incidence"]] < eps
        for dep in ("dsi", "lesion_length"):
            if dep in idx:
                y[no_inf, idx[dep]] = 0.0
    return y


def outcome_violations(y: torch.Tensor, names: list[str], eps: float = 1e-3) -> torch.Tensor:
    """Fraction of rows that violate the definitional outcome constraints (before projection)."""
    idx = {n: i for i, n in enumerate(names)}
    bad = torch.zeros(y.shape[0], dtype=torch.bool, device=y.device)
    for n, i in idx.items():
        lo_bad = y[:, i] < -eps
        hi_bad = (y[:, i] > 1 + eps) if n in ("incidence", "dsi") else torch.zeros_like(lo_bad)
        bad |= lo_bad | hi_bad
    if "incidence" in idx:
        no_inf = y[:, idx["incidence"]] < eps
        for dep in ("dsi", "lesion_length"):
            if dep in idx:
                bad |= no_inf & (y[:, idx[dep]] > eps)
    return bad.float().mean()
