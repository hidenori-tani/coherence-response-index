import math
import pandas as pd
from tmc.io_expr import compute_log2fc, compute_log2fc_robust, chemical_log2fc

def test_log2fc_basic():
    treat = {"A": 4.0, "B": 1.0}
    ctrl  = {"A": 1.0, "B": 4.0}
    out = compute_log2fc(treat, ctrl)
    assert math.isclose(out["A"], 2.0)
    assert math.isclose(out["B"], -2.0)

def test_log2fc_zero_uses_min_nonzero():
    # Published rule (PLoS2017/TAP2020): a zero value is replaced by the smallest
    # nonzero value in its OWN vector. In real ~20k-feature data that floor is tiny,
    # so a feature that drops to zero under treatment yields a strong negative log2FC.
    # The fixture includes a small nonzero (B=0.5) so the floor is representative.
    treat = {"A": 8.0, "B": 0.5, "Z": 0.0}   # min nonzero in treatment = 0.5
    ctrl  = {"A": 2.0, "B": 0.5, "Z": 4.0}
    out = compute_log2fc(treat, ctrl)
    assert math.isclose(out["A"], 2.0)
    assert math.isclose(out["Z"], math.log2(0.5 / 4.0))  # = -3.0, strong downregulation
    assert out["Z"] < 0

def test_log2fc_robust_no_floor():
    # Robust pseudocount estimator: a transcript dropping to zero is NOT floored to
    # an arbitrarily large ratio (unlike compute_log2fc). log2((0+1)/(4+1)) = log2(0.2).
    # A tiny nonzero (B=0.001) makes the zero-replacement floor representatively extreme.
    treat = {"A": 8.0, "B": 0.001, "Z": 0.0}   # zero-replacement floor = 0.001
    ctrl  = {"A": 2.0, "B": 0.001, "Z": 4.0}
    out = compute_log2fc_robust(treat, ctrl)
    assert math.isclose(out["A"], math.log2(9.0 / 3.0))     # log2((8+1)/(2+1))
    assert math.isclose(out["Z"], math.log2(1.0 / 5.0))     # bounded (~-2.32), not floored
    # robust drop is far milder than the zero-replacement floor for the same data
    zr = compute_log2fc(treat, ctrl)
    assert zr["Z"] < -10 and out["Z"] > zr["Z"]             # floor ~-11.97 vs robust ~-2.32

def test_chemical_log2fc_robust_and_cpm():
    mat = pd.DataFrame({"r1": [8, 0], "r2": [8, 0], "c1": [2, 4], "c2": [2, 4]},
                       index=["A", "Z"])
    acc = [{"run_accession": "r1", "chemical": "benzene"},
           {"run_accession": "r2", "chemical": "benzene"},
           {"run_accession": "c1", "chemical": "dmso_control"},
           {"run_accession": "c2", "chemical": "dmso_control"}]
    out = chemical_log2fc(mat, acc, "benzene", estimator=compute_log2fc_robust)
    assert math.isclose(out["A"], math.log2(9.0 / 3.0))     # log2((8+1)/(2+1))
    # cpm=True normalises each column to sum 1e6 before averaging; runs without error
    # and (pseudocount negligible at CPM scale) gives ~log2(3) for the dominant feature.
    out_cpm = chemical_log2fc(mat, acc, "benzene", estimator=compute_log2fc_robust, cpm=True)
    assert out_cpm["A"] > 1.0 and math.isfinite(out_cpm["A"])

def test_chemical_log2fc_means_replicates():
    mat = pd.DataFrame({"r1": [4, 1], "r2": [4, 1], "c1": [1, 4], "c2": [1, 4]},
                       index=["A", "B"])
    acc = [{"run_accession": "r1", "chemical": "benzene"},
           {"run_accession": "r2", "chemical": "benzene"},
           {"run_accession": "c1", "chemical": "dmso_control"},
           {"run_accession": "c2", "chemical": "dmso_control"}]
    out = chemical_log2fc(mat, acc, "benzene")
    assert round(out["A"], 6) == 2.0 and round(out["B"], 6) == -2.0
