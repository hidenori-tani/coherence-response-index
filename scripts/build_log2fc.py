"""Build per-chemical log2FC matrices from the recomputed count/abundance matrices.

Inputs (produced by scripts/quant_mirna.py and scripts/quant_mrna.py):
  results/mirna_counts.tsv   (rows=miRNA, cols=run)
  results/gene_tpm.tsv       (rows=gene,  cols=run)
Outputs:
  results/mirna_log2fc.tsv   (rows=miRNA, cols=9 chemicals)
  results/mrna_log2fc.tsv    (rows=gene,  cols=9 chemicals)
"""
import pandas as pd

from tmc.config import load_accessions, CHEMICALS
from tmc.io_expr import (load_matrix, chemical_log2fc,
                         compute_log2fc, compute_log2fc_robust)


def build(matrix_path, layer, out_path, acc_path="config/accessions.tsv",
          estimator=compute_log2fc, cpm=False):
    acc = [r for r in load_accessions(acc_path) if r["layer"] == layer]
    mat = load_matrix(matrix_path)
    df = pd.DataFrame({c: chemical_log2fc(mat, acc, c, estimator=estimator, cpm=cpm)
                       for c in CHEMICALS})
    df.index.name = "feature"
    df.to_csv(out_path, sep="\t")
    return df


if __name__ == "__main__":
    # Zero-replacement estimator (reproduces originals; used for floor-artifact demo)
    build("results/mirna_counts.tsv", "mirna", "results/mirna_log2fc.tsv")
    build("results/gene_tpm.tsv", "mrna", "results/mrna_log2fc.tsv")
    # Robust pseudocount estimator (used for the coherence analysis): mRNA on TPM,
    # miRNA on CPM (library-size normalised); log2((x+1)/(ctrl+1)) in both.
    build("results/mirna_counts.tsv", "mirna", "results/mirna_log2fc_robust.tsv",
          estimator=compute_log2fc_robust, cpm=True)
    build("results/gene_tpm.tsv", "mrna", "results/mrna_log2fc_robust.tsv",
          estimator=compute_log2fc_robust, cpm=False)
    print("wrote results/{mirna,mrna}_log2fc.tsv (zero-replacement) and "
          "results/{mirna,mrna}_log2fc_robust.tsv (robust pseudocount)")
