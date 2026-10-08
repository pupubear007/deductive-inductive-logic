"""Encoders for the three kinds of input: images, environment series, static attributes."""

from __future__ import annotations

import torch
from torch import nn


class ImageEncoder(nn.Module):
    """Small CNN for one image modality (growth chamber, greenhouse, field, drone, satellite).

    Deliberately small so it trains on a laptop. For real drone / satellite data swap in a
    pretrained backbone (e.g. a torchvision ResNet) with the same (B, C, H, W) -> (B, dim) contract.
    """

    def __init__(self, in_ch: int, dim: int, width: int = 32):
        super().__init__()
        w = width
        self.net = nn.Sequential(
            nn.Conv2d(in_ch, w, 3, padding=1), nn.GroupNorm(8, w), nn.SiLU(),
            nn.Conv2d(w, w, 3, stride=2, padding=1), nn.GroupNorm(8, w), nn.SiLU(),
            nn.Conv2d(w, 2 * w, 3, stride=2, padding=1), nn.GroupNorm(8, 2 * w), nn.SiLU(),
            nn.Conv2d(2 * w, 4 * w, 3, stride=2, padding=1), nn.GroupNorm(8, 4 * w), nn.SiLU(),
            nn.AdaptiveAvgPool2d(1), nn.Flatten(), nn.Linear(4 * w, dim),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


class SeriesEncoder(nn.Module):
    """GRU over daily environment readings up to the issue day (inputs are already truncated at
    tau, so the encoder can only see the past). Missing readings are zeros with a mask channel."""

    def __init__(self, n_vars: int, dim: int, hidden: int = 64):
        super().__init__()
        self.n_vars = n_vars
        self.gru = nn.GRU(2 * max(n_vars, 1), hidden, batch_first=True)
        self.out = nn.Linear(hidden, dim)

    def forward(self, x: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
        if self.n_vars == 0:
            return self.out.weight.new_zeros(x.shape[0], self.out.out_features)
        h, _ = self.gru(torch.cat([x, mask], dim=-1))
        return self.out(h[:, -1])


class TabularEncoder(nn.Module):
    def __init__(self, n_num: int, vocab_sizes: list[int], dim: int, emb: int = 16):
        super().__init__()
        self.embs = nn.ModuleList([nn.Embedding(v + 1, emb) for v in vocab_sizes])  # index 0 = unknown
        self.mlp = nn.Sequential(nn.Linear(n_num + emb * len(vocab_sizes), dim), nn.SiLU(), nn.Linear(dim, dim))
        self.n_in = n_num + emb * len(vocab_sizes)

    def forward(self, num: torch.Tensor, cat: torch.Tensor) -> torch.Tensor:
        if self.n_in == 0:
            return num.new_zeros(num.shape[0], self.mlp[-1].out_features)
        parts = [num] + [e(cat[:, i]) for i, e in enumerate(self.embs)]
        return self.mlp(torch.cat(parts, dim=-1))


class Fusion(nn.Module):
    """Concatenate modality embeddings with presence flags; missing modalities contribute zeros.
    During training whole modalities are dropped at random so the model learns to forecast from
    whatever subset of data a field actually has."""

    def __init__(self, n_parts: int, dim: int, p_drop: float = 0.2):
        super().__init__()
        self.p_drop = p_drop
        self.mlp = nn.Sequential(nn.Linear(n_parts * dim + n_parts, 2 * dim), nn.SiLU(), nn.Linear(2 * dim, dim))

    def forward(self, parts: list[torch.Tensor], present: list[torch.Tensor]) -> torch.Tensor:
        pres = torch.stack(present, dim=-1)
        if self.training and self.p_drop > 0:
            keep = (torch.rand_like(pres) > self.p_drop).float()
            pres = pres * keep
        x = torch.cat([p * pres[:, i:i + 1] for i, p in enumerate(parts)] + [pres], dim=-1)
        return self.mlp(x)
