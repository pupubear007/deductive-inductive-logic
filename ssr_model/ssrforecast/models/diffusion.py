"""Conditional denoising diffusion (DDPM) components.

Two generative heads, both conditioned on the fused embedding z of a sample:

* ``VectorDiffusion`` generates the continuous outcome vector of a world
  (incidence, DSI, lesion length, yield, ...). Draws are samples from P(outcomes | inputs), i.e. a
  learned prior over worlds (Definition 5.6 / S12); they are projected onto admissible worlds
  with ``mechanism.project_outcomes`` before use.
* ``MapDiffusion`` generates a 2-D disease map (e.g. incidence per grid cell derived from a drone
  or satellite image). This is where spatial spread / dispersal structure is learned.

Training minimises the usual noise-prediction loss, masked to the targets a sample actually has.
"""

from __future__ import annotations

import math

import torch
from torch import nn
from torch.nn import functional as F


def timestep_embedding(t: torch.Tensor, dim: int) -> torch.Tensor:
    half = dim // 2
    freqs = torch.exp(-math.log(10000) * torch.arange(half, device=t.device) / half)
    a = t.float()[:, None] * freqs[None]
    return torch.cat([a.sin(), a.cos()], dim=-1)


class Schedule(nn.Module):
    def __init__(self, steps: int):
        super().__init__()
        betas = torch.linspace(1e-4, 0.02, steps) * (1000 / steps)
        betas = betas.clamp(max=0.999)
        alphas = 1 - betas
        abar = torch.cumprod(alphas, 0)
        self.steps = steps
        self.register_buffer("betas", betas)
        self.register_buffer("alphas", alphas)
        self.register_buffer("abar", abar)

    def q_sample(self, x0: torch.Tensor, t: torch.Tensor, noise: torch.Tensor) -> torch.Tensor:
        shape = (-1,) + (1,) * (x0.dim() - 1)
        ab = self.abar[t].view(shape)
        return ab.sqrt() * x0 + (1 - ab).sqrt() * noise

    @torch.no_grad()
    def p_sample_loop(self, eps_fn, shape, device) -> torch.Tensor:
        x = torch.randn(shape, device=device)
        view = (-1,) + (1,) * (len(shape) - 1)
        for i in reversed(range(self.steps)):
            t = torch.full((shape[0],), i, device=device, dtype=torch.long)
            eps = eps_fn(x, t)
            a, ab, b = self.alphas[i], self.abar[i], self.betas[i]
            mean = (x - b / (1 - ab).sqrt() * eps) / a.sqrt()
            x = mean + (b.sqrt() * torch.randn_like(x) if i > 0 else 0.0)
        return x


# ---------------------------------------------------------------------------------------------

class VectorDiffusion(nn.Module):
    """DDPM over a d-dimensional outcome vector (values scaled to roughly [0, 1])."""

    def __init__(self, d: int, cond_dim: int, hidden: int = 256, steps: int = 100):
        super().__init__()
        self.d = d
        self.sched = Schedule(steps)
        self.t_dim = 64
        self.net = nn.Sequential(
            nn.Linear(d + cond_dim + self.t_dim, hidden), nn.SiLU(),
            nn.Linear(hidden, hidden), nn.SiLU(),
            nn.Linear(hidden, d),
        )

    def eps(self, x: torch.Tensor, t: torch.Tensor, cond: torch.Tensor) -> torch.Tensor:
        return self.net(torch.cat([x, cond, timestep_embedding(t, self.t_dim)], dim=-1))

    def loss(self, y: torch.Tensor, mask: torch.Tensor, cond: torch.Tensor) -> torch.Tensor:
        """Noise-prediction MSE on the observed coordinates only. Data are mapped to [-1, 1]."""
        keep = mask.sum(-1) > 0
        if not keep.any():
            return y.new_zeros(())
        y, mask, cond = y[keep], mask[keep], cond[keep]
        x0 = (2 * y - 1) * mask  # unobserved coordinates pinned at 0 in data space
        t = torch.randint(0, self.sched.steps, (x0.shape[0],), device=x0.device)
        noise = torch.randn_like(x0)
        pred = self.eps(self.sched.q_sample(x0, t, noise), t, cond)
        return ((pred - noise) ** 2 * mask).sum() / mask.sum().clamp_min(1)

    @torch.no_grad()
    def sample(self, cond: torch.Tensor) -> torch.Tensor:
        x = self.sched.p_sample_loop(lambda x, t: self.eps(x, t, cond), (cond.shape[0], self.d), cond.device)
        return (x + 1) / 2


# ---------------------------------------------------------------------------------------------

class _Block(nn.Module):
    def __init__(self, cin: int, cout: int, emb: int):
        super().__init__()
        self.c1 = nn.Conv2d(cin, cout, 3, padding=1)
        self.c2 = nn.Conv2d(cout, cout, 3, padding=1)
        self.n1, self.n2 = nn.GroupNorm(8, cout), nn.GroupNorm(8, cout)
        self.film = nn.Linear(emb, 2 * cout)
        self.skip = nn.Conv2d(cin, cout, 1) if cin != cout else nn.Identity()

    def forward(self, x: torch.Tensor, e: torch.Tensor) -> torch.Tensor:
        h = F.silu(self.n1(self.c1(x)))
        g, b = self.film(e).chunk(2, dim=-1)
        h = h * (1 + g[:, :, None, None]) + b[:, :, None, None]
        h = F.silu(self.n2(self.c2(h)))
        return h + self.skip(x)


class SmallUNet(nn.Module):
    def __init__(self, ch: int, cond_dim: int, width: int = 32):
        super().__init__()
        w, emb = width, 128
        self.t_dim = 64
        self.emb = nn.Sequential(nn.Linear(self.t_dim + cond_dim, emb), nn.SiLU(), nn.Linear(emb, emb))
        self.inp = nn.Conv2d(ch, w, 3, padding=1)
        self.d1, self.d2 = _Block(w, w, emb), _Block(w, 2 * w, emb)
        self.mid = _Block(2 * w, 2 * w, emb)
        self.u2, self.u1 = _Block(4 * w, w, emb), _Block(2 * w, w, emb)
        self.out = nn.Conv2d(w, ch, 3, padding=1)

    def forward(self, x: torch.Tensor, t: torch.Tensor, cond: torch.Tensor) -> torch.Tensor:
        e = self.emb(torch.cat([timestep_embedding(t, self.t_dim), cond], dim=-1))
        h0 = self.inp(x)
        h1 = self.d1(h0, e)
        h2 = self.d2(F.avg_pool2d(h1, 2), e)
        m = self.mid(F.avg_pool2d(h2, 2), e)
        u2 = self.u2(torch.cat([F.interpolate(m, scale_factor=2), h2], 1), e)
        u1 = self.u1(torch.cat([F.interpolate(u2, scale_factor=2), h1], 1), e)
        return self.out(u1)


class MapDiffusion(nn.Module):
    """DDPM over 1 x S x S disease maps in [0, 1] (S divisible by 4)."""

    def __init__(self, size: int, cond_dim: int, steps: int = 100, width: int = 32):
        super().__init__()
        assert size % 4 == 0, "map_size must be divisible by 4"
        self.size = size
        self.sched = Schedule(steps)
        self.unet = SmallUNet(1, cond_dim, width)

    def loss(self, y: torch.Tensor, mask: torch.Tensor, cond: torch.Tensor) -> torch.Tensor:
        keep = mask > 0
        if not keep.any():
            return y.new_zeros(())
        x0 = 2 * y[keep] - 1
        c = cond[keep]
        t = torch.randint(0, self.sched.steps, (x0.shape[0],), device=x0.device)
        noise = torch.randn_like(x0)
        return F.mse_loss(self.unet(self.sched.q_sample(x0, t, noise), t, c), noise)

    @torch.no_grad()
    def sample(self, cond: torch.Tensor) -> torch.Tensor:
        x = self.sched.p_sample_loop(lambda x, t: self.unet(x, t, cond),
                                     (cond.shape[0], 1, self.size, self.size), cond.device)
        return ((x + 1) / 2).clamp(0, 1)
