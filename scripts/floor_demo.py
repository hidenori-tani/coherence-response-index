"""Floor-artifact demonstration on a SINGLE quantification (answers the
apples-to-oranges concern: the inflation is the estimator, not the pipeline).

For benzene, we take the SAME kallisto gene-level TPM matrix and compute the
per-gene log2 fold change TWO ways:
  (a) zero-replacement log-ratio  r = log2(treat/ctrl), zeros -> min nonzero
      (the originals' estimator), and
  (b) robust pseudocount          log2((treat+1)/(ctrl+1)).
We then show that the extreme fold changes appear ONLY under (a) and that the
published benzene 'biomarker' anchors collapse to near-zero under (b), while the
data themselves are unchanged. Prints a small table for Results 3.1.
"""
import numpy as np
import pandas as pd

from tmc.config import load_accessions, CONTROL
from tmc.io_expr import load_matrix, compute_log2fc, compute_log2fc_robust

CHEM = "benzene"


def main():
    gt = load_matrix("results/gene_tpm.tsv")
    acc = [r for r in load_accessions("config/accessions.tsv") if r["layer"] == "mrna"]
    runs = lambda c: [r["run_accession"] for r in acc if r["chemical"] == c]
    treat = gt[runs(CHEM)].mean(axis=1).to_dict()
    ctrl = gt[runs(CONTROL)].mean(axis=1).to_dict()
    zr = compute_log2fc(treat, ctrl)
    rb = compute_log2fc_robust(treat, ctrl)

    z = pd.Series(zr); r = pd.Series(rb)
    print(f"=== {CHEM}: same kallisto TPM, two estimators ===")
    print(f"genes: {len(z)}")
    print(f"|log2FC| >= 5 (>=32-fold): zero-replacement = {(z.abs()>=5).sum()},"
          f"  robust = {(r.abs()>=5).sum()}")
    print(f"max |log2FC|: zero-replacement = {z.abs().max():.1f},"
          f"  robust = {r.abs().max():.2f}")
    print(f"SD of log2FC: zero-replacement = {z.std():.2f},  robust = {r.std():.2f}")

    # published benzene anchors mapped to symbols
    r2s = pd.read_csv("config/refseq2symbol.tsv", sep="\t")
    m = dict(zip(r2s.refseq, r2s.gene_symbol))
    p = pd.read_csv("config/published_top30.tsv", sep="\t")
    pb = p[(p.layer == "mrna") & (p.chemical == CHEM)].copy()
    pb["sym"] = pb.feature.map(m)
    print("\nPublished benzene anchors (same kallisto TPM):")
    print(f"{'refseq':>14} {'symbol':>14} {'pub_fold':>10} {'ctrlTPM':>8} {'benzTPM':>8}"
          f" {'zeroRepl':>9} {'robust':>8}")
    for _, row in pb.iterrows():
        s = row["sym"]
        if s in z.index:
            print(f"{row.feature:>14} {str(s):>14} {row.published_value:>10.4g}"
                  f" {ctrl.get(s,0):>8.1f} {treat.get(s,0):>8.1f}"
                  f" {z[s]:>+9.2f} {r[s]:>+8.2f}")


if __name__ == "__main__":
    main()
