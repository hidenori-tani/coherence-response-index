"""Sanity check: do our recomputed log2FC values agree (in sign and rank) with the
published Top30/RT-qPCR anchors from PLoS2017 (mRNA) and TAP2020 (miRNA)?

Anchors are in config/published_top30.tsv. mRNA anchors are RefSeq IDs with fold
values; they are mapped to our gene-level matrix via an optional refseq->symbol
lookup (config/refseq2symbol.tsv, built in the data phase). miRNA anchors are
miRBase IDs with log2 values and map directly.

If sign concordance is poor, STOP and debug quantification before trusting any
coherence result.
"""
import math
import os
import pandas as pd
from scipy.stats import spearmanr

from tmc.io_expr import load_matrix


def _load_refseq2symbol(path="config/refseq2symbol.tsv"):
    if not os.path.exists(path):
        return {}
    m = pd.read_csv(path, sep="\t")
    return dict(zip(m.iloc[:, 0], m.iloc[:, 1]))


def _to_log2(value, value_type):
    return value if value_type == "log2" else math.log2(value)


def main():
    anchors = pd.read_csv("config/published_top30.tsv", sep="\t")
    mirna = load_matrix("results/mirna_log2fc.tsv")
    mrna = load_matrix("results/mrna_log2fc.tsv")
    r2s = _load_refseq2symbol()

    recomputed, published, rows = [], [], []
    for _, a in anchors.iterrows():
        chem, feat = a["chemical"], a["feature"]
        pub_log2 = _to_log2(float(a["published_value"]), a["value_type"])
        if a["layer"] == "mirna":
            key = feat if feat in mirna.index else None
            mat = mirna
        else:
            key = r2s.get(feat, feat)
            key = key if key in mrna.index else None
            mat = mrna
        if key is None or chem not in mat.columns:
            rows.append((a["layer"], feat, "UNMAPPED", "", ""))
            continue
        rec = float(mat.loc[key, chem])
        recomputed.append(rec)
        published.append(pub_log2)
        same_sign = (rec > 0) == (pub_log2 > 0)
        rows.append((a["layer"], feat, f"{rec:+.2f}", f"{pub_log2:+.2f}",
                     "OK" if same_sign else "SIGN-MISMATCH"))

    print(f"{'layer':6} {'feature':16} {'recomputed':>10} {'published':>10}  status")
    for layer, feat, rec, pub, status in rows:
        print(f"{layer:6} {feat:16} {rec:>10} {pub:>10}  {status}")

    mapped = [s for *_, s in rows if s in ("OK", "SIGN-MISMATCH")]
    if mapped:
        conc = sum(s == "OK" for s in mapped) / len(mapped)
        rho = spearmanr(recomputed, published).statistic if len(recomputed) > 2 else float("nan")
        print(f"\nsign concordance = {conc:.0%} over {len(mapped)} mapped anchors; "
              f"Spearman(recomputed, published log2) = {rho:.2f}")
        if conc < 0.9:
            print("WARNING: concordance < 90% -> debug quantification before proceeding.")
    else:
        print("\nNo anchors mapped (need config/refseq2symbol.tsv for mRNA).")


if __name__ == "__main__":
    main()
