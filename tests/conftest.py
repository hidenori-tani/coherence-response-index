import pytest

@pytest.fixture
def tiny_log2fc():
    # gene -> log2FC for one chemical
    return {"GeneA": -2.0, "GeneB": -1.5, "GeneC": 0.1, "GeneD": 1.8, "GeneE": -0.2}

@pytest.fixture
def tiny_targets():
    # miRNA -> set of predicted target genes
    return {"mmu-miR-1": {"GeneA", "GeneB"}, "mmu-miR-2": {"GeneC", "GeneD"}}
