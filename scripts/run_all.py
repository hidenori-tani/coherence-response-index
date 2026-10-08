"""End-to-end analysis (after quantification + build_log2fc).

Stages:
  coherence    -> results/coherence.tsv      (per chemical x responsive miRNA)
  specificity  -> results/specificity.tsv    (per miRNA: generic / benzene contrast)

Usage:
  python scripts/run_all.py coherence
  python scripts/run_all.py specificity
  python scripts/run_all.py all
"""
import sys
import numpy as np
import pandas as pd

from tmc.config import load_accessions, CONTROL
from tmc.io_expr import load_matrix
from tmc.targets import load_mirdb
from tmc.pipeline import run_coherence, summarize_specificity

COHERENCE_OUT = "results/coherence.tsv"
SPECIFICITY_OUT = "results/specificity.tsv"


# Expression filters: exclude floor-inflated low-expression features whose
# zero-replaced log2FC is unreliable (see docs/STATUS.md sanity findings).
MIN_CTRL_TPM = 1.0    # mRNA: mean control TPM
MIN_CTRL_CPM = 10.0   # miRNA: mean control CPM


def _ctrl_runs(layer):
    acc = [r for r in load_accessions("config/accessions.tsv") if r["layer"] == layer]
    return [r["run_accession"] for r in acc if r["chemical"] == CONTROL]


def _baseline_expr(expressed_genes=None):
    """Covariate for matched-null = mean control TPM per gene (log10-scaled)."""
    tpm = load_matrix("results/gene_tpm.tsv")
    base = tpm[_ctrl_runs("mrna")].mean(axis=1)
    if expressed_genes is not None:
        base = base[base.index.isin(expressed_genes)]
    return {g: float(np.log10(v + 1.0)) for g, v in base.items()}


def _expressed_genes():
    tpm = load_matrix("results/gene_tpm.tsv")
    m = tpm[_ctrl_runs("mrna")].mean(axis=1)
    return set(m[m >= MIN_CTRL_TPM].index)


def _expressed_mirnas():
    cnt = load_matrix("results/mirna_counts.tsv")
    cpm = cnt / cnt.sum() * 1e6
    m = cpm[_ctrl_runs("mirna")].mean(axis=1)
    return set(m[m >= MIN_CTRL_CPM].index)


def stage_coherence():
    # Robust pseudocount log2FC (not the floor-prone zero-replacement estimator,
    # which is reserved for the floor-artifact demonstration in scripts/floor_demo.py).
    mirna = load_matrix("results/mirna_log2fc_robust.tsv")
    mrna = load_matrix("results/mrna_log2fc_robust.tsv")
    # restrict to expressed features (floor-inflated low-expression excluded)
    exg, exm = _expressed_genes(), _expressed_mirnas()
    mrna = mrna[mrna.index.isin(exg)]
    mirna = mirna[mirna.index.isin(exm)]
    print(f"expressed: {len(mrna)} genes (>= {MIN_CTRL_TPM} TPM), "
          f"{len(mirna)} miRNAs (>= {MIN_CTRL_CPM} CPM)")
    mb = load_mirdb("config/mirdb_mmu.tsv", score_min=80)
    base = _baseline_expr(exg)
    # Primary analysis uses miRDB only (the database the source miRNA study used);
    # TargetScan is a separate robustness check (scripts/sensitivity.py). Passing an
    # empty TargetScan dict keeps target_mode="union" = miRDB targets only.
    df = run_coherence(mirna, mrna, {}, mb, base,
                       responsive_thresh=1.0, target_mode="union",
                       min_targets=5, n_per_target=20, seed=0)
    df.to_csv(COHERENCE_OUT, sep="\t", index=False)
    print(f"wrote {COHERENCE_OUT}: {len(df)} rows; "
          f"class counts = {df['class'].value_counts().to_dict()}")


def stage_specificity():
    df = pd.read_csv(COHERENCE_OUT, sep="\t")
    spec = summarize_specificity(df, k=5, benzene_key="benzene")
    spec.to_csv(SPECIFICITY_OUT, sep="\t", index=False)
    print(f"wrote {SPECIFICITY_OUT}: {len(spec)} miRNAs; "
          f"generic = {int(spec['generic'].sum())}")


if __name__ == "__main__":
    stage = sys.argv[1] if len(sys.argv) > 1 else "all"
    if stage in ("coherence", "all"):
        stage_coherence()
    if stage in ("specificity", "all"):
        stage_specificity()
