"""The mechanism's hard constraints hold by construction (mechanism notes S1, S5, S7)."""

import torch

from ssrforecast.mechanism import chain_probabilities, outcome_violations, project_outcomes, violations
from ssrforecast.models.full import CHAIN_LOGITS


def random_logits(n=512, scale=4.0):
    g = torch.Generator().manual_seed(0)
    return {k: torch.randn(n, generator=g) * scale for k in CHAIN_LOGITS}


def test_chain_is_monotone():
    """p_inf <= p_court <= p_spore and p_inf <= p_micro: no mass on worlds violating K."""
    c = chain_probabilities(random_logits())
    assert torch.all(c.infection <= c.court + 1e-7)
    assert torch.all(c.court <= c.spore + 1e-7)
    assert torch.all(c.court <= c.flowering + 1e-7)
    assert torch.all(c.infection <= c.microclimate + 1e-7)
    assert torch.all(c.apothecia <= c.spore + 1e-7)


def test_gate_failure_blocks_infection():
    """Prop. 4.4 / S7: if the microclimate gate is closed, infection has probability ~0."""
    lg = random_logits()
    lg["microclimate"] = torch.full_like(lg["microclimate"], -30.0)
    assert torch.all(chain_probabilities(lg).infection < 1e-10)


def test_no_spores_no_court_without_intervention():
    lg = random_logits()
    lg["apothecia"] = torch.full_like(lg["apothecia"], -30.0)
    lg["external"] = torch.full_like(lg["external"], -30.0)
    c = chain_probabilities(lg)
    assert torch.all(c.court < 1e-10)
    # do(court = 1): inoculation bypasses the spore level
    c2 = chain_probabilities(lg, do_court=torch.ones(len(c.court)))
    assert torch.all(c2.court == 1.0)
    assert torch.all(c2.infection <= c2.microclimate + 1e-7)


def test_violation_counter():
    lv = {k: torch.tensor([1, 0]) for k in ("apothecia", "spore", "flowering", "court", "microclimate", "infection")}
    assert violations(lv).tolist() == [0, 0]
    lv["court"] = torch.tensor([0, 0])  # infection without court in world 0
    assert violations(lv).tolist() == [1, 0]


def test_projection_makes_outcomes_admissible():
    names = ["incidence", "dsi", "lesion_length", "yield"]
    y = torch.tensor([[0.0, 0.3, 5.0, 1.0], [1.4, -0.2, 0.1, -1.0], [0.2, 0.1, 0.3, 0.9]])
    assert outcome_violations(y, names) > 0
    assert outcome_violations(project_outcomes(y, names), names) == 0
