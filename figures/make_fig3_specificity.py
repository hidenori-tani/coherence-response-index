"""Fig 3: (a) generic-vs-specific structure (how many chemicals each miRNA is
productive in); (b) top benzene-specific miRNAs by CRI contrast. Rendered from
results/specificity.tsv."""
import matplotlib.pyplot as plt
import pandas as pd

from _common import SPECIFICITY, require, save


def main():
    spec = pd.read_csv(require(SPECIFICITY), sep="\t")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))

    counts = spec["n_chemicals_productive"]
    active = counts[counts > 0]
    if len(active):
        bins = range(0, int(active.max()) + 2)
        ax1.hist(active, bins=bins, color="slateblue", align="left", rwidth=0.85)
    ax1.set_xlabel("# chemicals in which miRNA is productive")
    ax1.set_ylabel("# miRNAs")
    ax1.set_title("(a) generic vs chemical-specific")

    top = spec.sort_values("benzene_contrast", ascending=False).head(12)
    ax2.barh(range(len(top))[::-1], top["benzene_contrast"], color="crimson")
    ax2.set_yticks(range(len(top))[::-1])
    ax2.set_yticklabels(top["miRNA"], fontsize=7)
    ax2.set_xlabel("benzene CRI contrast (benzene - mean others)")
    ax2.set_title("(b) top benzene-specific miRNAs")

    fig.tight_layout()
    save(fig, "fig3_specificity")


if __name__ == "__main__":
    main()
