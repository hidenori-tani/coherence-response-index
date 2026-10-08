"""Fig 1: target-vs-background repression AUC per chemical (do toxic-response
miRNAs' predicted targets shift? AUC>0.5 = targets more repressed than matched
background). Rendered from results/coherence.tsv."""
import matplotlib.pyplot as plt
import pandas as pd

from _common import COHERENCE, require, save


def main():
    df = pd.read_csv(require(COHERENCE), sep="\t")
    df = df[df["auc"].notna()]
    chems = sorted(df["chemical"].unique())
    data = [df.loc[df.chemical == c, "auc"].values for c in chems]

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.axhline(0.5, color="grey", ls="--", lw=1, label="no shift (AUC=0.5)")
    for i, (c, vals) in enumerate(zip(chems, data)):
        x = [i + 1] * len(vals)
        ax.scatter(x, vals, s=12, alpha=0.5,
                   color="crimson" if c == "benzene" else "steelblue")
    ax.boxplot(data, showfliers=False, widths=0.5)
    ax.set_xticks(range(1, len(chems) + 1))
    ax.set_xticklabels(chems, rotation=45, ha="right")
    ax.set_ylabel("target-vs-background repression AUC")
    ax.set_title("miRNA-target directional shift per chemical")
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    save(fig, "fig1_coherence_auc")


if __name__ == "__main__":
    main()
