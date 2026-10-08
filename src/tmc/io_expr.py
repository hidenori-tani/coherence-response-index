import math
import pandas as pd

def _min_nonzero(values):
    nz = [v for v in values if v > 0]
    return min(nz) if nz else 1e-9

def compute_log2fc(treatment, control):
    """r(x,y)=log2(y/x); zeros replaced by smallest nonzero in their own vector
    (matches Tani2017/TAP2020 convention). Floor-prone at low counts; used only to
    REPRODUCE the originals' estimator for the floor-artifact demonstration (S3.1).
    The coherence analysis uses compute_log2fc_robust instead."""
    floor_t = _min_nonzero(treatment.values())
    floor_c = _min_nonzero(control.values())
    keys = set(treatment) | set(control)
    out = {}
    for k in keys:
        y = treatment.get(k, 0.0) or floor_t
        x = control.get(k, 0.0) or floor_c
        out[k] = math.log2(y / x)
    return out


def compute_log2fc_robust(treatment, control, pseudo=1.0):
    """Low-count-stable log2 fold change: log2((y+pseudo)/(x+pseudo)). No zero
    replacement, so a transcript dropping to 0 in one condition is not floored to
    an arbitrarily large ratio. Used for all coherence analyses."""
    keys = set(treatment) | set(control)
    out = {}
    for k in keys:
        y = treatment.get(k, 0.0)
        x = control.get(k, 0.0)
        out[k] = math.log2((y + pseudo) / (x + pseudo))
    return out

def load_matrix(path):
    return pd.read_csv(path, sep="\t", index_col=0)

def chemical_log2fc(matrix, accessions, chemical, control="dmso_control",
                    estimator=compute_log2fc, cpm=False):
    """Mean-of-replicates log2FC for one chemical vs control. accessions: list of
    dict rows (same layer) with run_accession/chemical.
    estimator: compute_log2fc (zero-replacement, default) or compute_log2fc_robust.
    cpm: if True, column-normalise counts to counts-per-million before averaging
    (used for the miRNA layer in the robust build so library size is accounted for)."""
    runs = lambda c: [r["run_accession"] for r in accessions if r["chemical"] == c]
    mat = matrix / matrix.sum() * 1e6 if cpm else matrix
    treat = mat[runs(chemical)].mean(axis=1).to_dict()
    ctrl = mat[runs(control)].mean(axis=1).to_dict()
    return estimator(treat, ctrl)
