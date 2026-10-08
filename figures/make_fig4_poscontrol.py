"""Fig 4: positive control. The same CRI pipeline applied to an independent
miRNA-mimic transfection dataset (GSE269995, three neural cell types) detects
coherent target repression for the transfected miRNAs but not for a panel of
random non-transfected miRNAs treated identically. Rendered from
results/posctrl_coherence.tsv (no hardcoded numbers)."""
import matplotlib.pyplot as plt
import pandas as pd

from _common import require, save, wilson_ci

POSCTRL = "results/posctrl_coherence.tsv"


def main():
    df = pd.read_csv(require(POSCTRL), sep="\t")
    tested = df[df["class"].isin(["productive", "passenger", "anti_coherent"])].copy()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4))

    # Panel A: AUC by group (target-vs-background repression)
    groups = ["transfected", "random"]
    data = [tested.loc[tested.group == g, "auc"].values for g in groups]
    ax1.axhline(0.5, color="grey", ls="--", lw=1, label="no shift (AUC=0.5)")
    for i, (g, vals) in enumerate(zip(groups, data)):
        ax1.scatter([i + 1] * len(vals), vals, s=18, alpha=0.6,
                    color="seagreen" if g == "transfected" else "darkgrey")
    ax1.boxplot(data, showfliers=False, widths=0.5)
    ax1.set_xticks([1, 2])
    ax1.set_xticklabels(["transfected\nmiRNAs", "random\nmiRNAs"])
    ax1.set_ylabel("target-vs-background repression AUC")
    ax1.set_title("A  Target repression detected for transfected miRNAs")
    ax1.legend(loc="lower right", fontsize=8)

    # Panel B: productive fraction by group, with Wilson CI
    fracs, los, his, labels = [], [], [], []
    for g in groups:
        sub = tested[tested.group == g]
        k = int((sub["class"] == "productive").sum())
        n = len(sub)
        p = k / n if n else 0.0
        lo, hi = wilson_ci(k, n)
        fracs.append(p); los.append(max(0.0, p - lo)); his.append(max(0.0, hi - p))
        labels.append(f"{g}\n({k}/{n})")
    ax2.bar([0, 1], fracs, yerr=[los, his], capsize=5,
            color=["seagreen", "darkgrey"], width=0.6)
    ax2.set_xticks([0, 1]); ax2.set_xticklabels(labels)
    ax2.set_ylabel("fraction productive")
    ax2.set_ylim(0, 1.05)
    ax2.set_title("B  Productive fraction (CRI verdict)")

    fig.tight_layout()
    save(fig, "fig4_poscontrol")


if __name__ == "__main__":
    main()
