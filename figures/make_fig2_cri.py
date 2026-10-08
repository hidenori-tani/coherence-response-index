"""Fig 2: (a) CRI distribution coloured by class (the audit's core null: CRI
centred on 0, essentially no coherent responses); (b) productive fraction per
chemical with Wilson CIs. Rendered from results/coherence.tsv."""
import matplotlib.pyplot as plt
import pandas as pd

from _common import COHERENCE, require, save, wilson_ci

CLASSES = ["passenger", "anti_coherent", "productive"]
COLORS = {"passenger": "lightgrey", "anti_coherent": "indianred", "productive": "seagreen"}


def main():
    df = pd.read_csv(require(COHERENCE), sep="\t")
    tested = df[df["class"].isin(CLASSES)].copy()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))

    series = [tested.loc[tested["class"] == c, "cri"].dropna() for c in CLASSES]
    ax1.hist(series, bins=25, stacked=True,
             label=[f"{c} (n={len(s)})" for c, s in zip(CLASSES, series)],
             color=[COLORS[c] for c in CLASSES])
    ax1.axvline(0, color="k", lw=0.8)
    ax1.set_xlabel("Coherence Response Index (CRI)")
    ax1.set_ylabel("miRNA x chemical responses")
    ax1.set_title("(a) CRI distribution by class")
    ax1.legend(fontsize=8)

    # productive fraction per chemical (here ~0 everywhere = the audit finding)
    chems = sorted(tested["chemical"].unique())
    fracs, los, his = [], [], []
    for c in chems:
        sub = tested[tested.chemical == c]
        k = int((sub["class"] == "productive").sum())
        n = len(sub)
        f = k / n if n else 0.0
        lo, hi = wilson_ci(k, n)
        fracs.append(f)
        los.append(max(0.0, f - lo))
        his.append(max(0.0, hi - f))
    colors = ["crimson" if c == "benzene" else "steelblue" for c in chems]
    ax2.bar(range(len(chems)), fracs, yerr=[los, his], capsize=3, color=colors)
    ax2.set_xticks(range(len(chems)))
    ax2.set_xticklabels(chems, rotation=45, ha="right")
    ax2.set_ylabel("fraction productive (95% Wilson CI)")
    ax2.set_ylim(0, 1)
    ax2.set_title("(b) productive fraction per chemical")

    fig.tight_layout()
    save(fig, "fig2_cri")


if __name__ == "__main__":
    main()
