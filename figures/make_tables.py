"""Table 1: per-chemical coherence summary (n responsive miRNA, % productive,
median AUC). Table 2: top productive miRNA->target responses by CRI.
Rendered from results/coherence.tsv into results/table1.tsv and results/table2.tsv."""
import pandas as pd

from _common import COHERENCE, require


def main():
    df = pd.read_csv(require(COHERENCE), sep="\t")
    # "tested" = every evaluated test (>=5 expressed targets); this includes
    # anti-coherent calls, which are evaluated tests with a significant shift in
    # the wrong direction. Only 'untested' (too few targets) is excluded.
    evaluated = ["productive", "passenger", "anti_coherent"]
    tested = df[df["class"].isin(evaluated)]

    rows = []
    for c, g in df.groupby("chemical"):
        gt = g[g["class"].isin(evaluated)]
        n_resp = len(g)
        n_tested = len(gt)
        n_prod = int((gt["class"] == "productive").sum())
        rows.append({
            "chemical": c,
            "n_responsive_miRNA": n_resp,
            "n_tested": n_tested,
            "n_productive": n_prod,
            "pct_productive": round(100 * n_prod / n_tested, 1) if n_tested else 0.0,
            "median_AUC": round(gt["auc"].median(), 3) if n_tested else float("nan"),
            "is_carcinogen": (c == "benzene"),
        })
    t1 = pd.DataFrame(rows).sort_values("chemical")
    t1.to_csv("results/table1.tsv", sep="\t", index=False)

    t2 = (tested[tested["class"] == "productive"]
          .sort_values("cri", ascending=False)
          .loc[:, ["chemical", "miRNA", "mirna_log2fc", "n_targets",
                   "median_diff", "auc", "fdr", "cri"]]
          .head(25)
          .round(3))
    t2.to_csv("results/table2.tsv", sep="\t", index=False)
    print(f"wrote results/table1.tsv ({len(t1)} chemicals) and "
          f"results/table2.tsv ({len(t2)} top productive responses)")


if __name__ == "__main__":
    main()
