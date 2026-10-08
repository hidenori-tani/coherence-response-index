import math
from tmc.audit import observed_support_rate


def test_observed_support_rate():
    # GO claim genes that are actually repressed / all claim genes
    claim_genes = {"GeneA", "GeneB", "GeneC", "GeneD"}
    repressed = {"GeneA", "GeneB"}  # observed down in matched chemical
    assert observed_support_rate(claim_genes, repressed) == 0.5


def test_observed_support_rate_full_and_none():
    claim = {"A", "B"}
    assert observed_support_rate(claim, {"A", "B", "C"}) == 1.0
    assert observed_support_rate(claim, set()) == 0.0


def test_observed_support_rate_empty_claim_is_nan():
    assert math.isnan(observed_support_rate(set(), {"A"}))
