"""Figures 1-5 for the CBI revision (R1), drawn with research/scripts/figstyle.py.

The content of Figures 1-3 and of Figure 4A-B is the submitted content, redrawn to
the Elsevier artwork rules (TrueType Arial, >= 7.8 pt at print size, print widths 90 / 140 /
190 mm, nothing inside the Editorial Manager stamp band, colour-vision-safe blue/vermillion
with grey). Figure 4C is new: the Dicer1-knockout positive control (R1-1, R3-3, R3-7).
Figure 5 is recomputed: the submitted spike-in shifted every synthetic target set down; the
revision shifts it in each test's coherent direction (results/revision_r1/power_coherent.tsv).
Every number is read from results/*.tsv at run time.

  python3 figures/make_figs_R1.py        # from the project root; writes submission_R1/fig/
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "scripts"))
from _common import require, wilson_ci  # noqa: E402
from figstyle import (use_journal_style, savefig, figure_label, mm,  # noqa: E402
                      OKABE_ITO, WIDTHS_MM)

OUT = "submission_R1/fig"
BLUE, VERM = OKABE_ITO["blue"], OKABE_ITO["vermillion"]
GREY, LIGHT = "#8C8C8C", "#C8C8C8"
LABEL_TOP_MM = 9.5      # panel letters sit below the EM stamp band (y 1.49-8.55 mm)
CHEM = {"benzene": "Benzene", "chloroform": "Chloroform", "dehp": "DEHP",
        "p_cresol": "$p$-Cresol", "p_dcb": "$p$-Dichlorobenzene", "phenol": "Phenol",
        "pyrocatechol": "Pyrocatechol", "tce": "Trichloroethylene",
        "tnbp": "Tri-$n$-butyl phosphate"}
TESTED = ["productive", "passenger", "anti_coherent"]
LS = 1.6                # line spacing of two-line axis labels (lines must not touch)


def canvas(column, height_mm):
    use_journal_style("elsevier")
    plt.rcParams.update({"xtick.major.pad": 4.5, "ytick.major.pad": 4.5})
    w = WIDTHS_MM["elsevier"][column]
    return plt.figure(figsize=(mm(w), mm(height_mm))), w, height_mm


def axes_mm(fig, w, h, left, top, width, height):
    """Axes placed in millimetres from the top-left corner of the figure."""
    return fig.add_axes([left / w, 1 - (top + height) / h, width / w, height / h])


def jitter(n, seed, width=0.12):
    return np.random.default_rng(seed).uniform(-width, width, n)


def box_and_points(ax, groups, colour):
    for i, vals in enumerate(groups):
        ax.scatter(i + 1 + jitter(len(vals), i), vals, s=7, color=colour, alpha=0.55,
                   linewidths=0, zorder=2)
    ax.boxplot(groups, showfliers=False, widths=0.5, zorder=3,
               medianprops={"color": "black", "linewidth": 1.0},
               boxprops={"linewidth": 0.7}, whiskerprops={"linewidth": 0.7},
               capprops={"linewidth": 0.7})
    ax.axhline(0.5, color=GREY, ls="--", lw=0.8, zorder=1)


def fig1():
    df = pd.read_csv(require("results/coherence.tsv"), sep="\t")
    df = df[df["class"].isin(TESTED)]
    chems = sorted(df["chemical"].unique())
    fig, w, h = canvas("onehalf", 82)
    ax = axes_mm(fig, w, h, 16, 11, w - 20, 44)
    box_and_points(ax, [df.loc[df.chemical == c, "auc"].values for c in chems], BLUE)
    ax.set_xticks(range(1, len(chems) + 1))
    ax.set_xticklabels([CHEM[c] for c in chems], rotation=40, ha="right")
    ax.set_ylabel("Target-vs-background\nrepression AUC", linespacing=LS)
    ax.set_xlim(0.4, len(chems) + 0.6)
    return fig


def productive_fraction(ax, tested, key, order, labels, colour):
    fr, lo, hi = [], [], []
    for g in order:
        sub = tested[tested[key] == g]
        k, n = int((sub["class"] == "productive").sum()), len(sub)
        p = k / n if n else 0.0
        a, b = wilson_ci(k, n)
        fr.append(p); lo.append(max(0.0, p - a)); hi.append(max(0.0, b - p))
    x = np.arange(len(order))
    ax.bar(x, fr, yerr=[lo, hi], capsize=2.5, color=colour, width=0.6,
           error_kw={"elinewidth": 0.7, "capthick": 0.7,
                     "clip_on": False})  # caps at 0 sit on the axis line
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Fraction productive\n(95% Wilson interval)", linespacing=LS)


def fig2():
    df = pd.read_csv(require("results/coherence.tsv"), sep="\t")
    t = df[df["class"].isin(TESTED)].copy()
    fig, w, h = canvas("double", 88)
    a = axes_mm(fig, w, h, 17, 15, 62, 44)
    order = ["passenger", "anti_coherent", "productive"]
    names = {"passenger": "Passenger", "anti_coherent": "Anti-coherent",
             "productive": "Productive"}
    series = [t.loc[t["class"] == c, "cri"].dropna() for c in order]
    a.hist(series, bins=25, stacked=True, color=[LIGHT, VERM, BLUE],
           label=[f"{names[c]} (n = {len(s)})" for c, s in zip(order, series)])
    a.axvline(0, color="black", lw=0.7)
    a.set_xlabel("Coherence Response Index (CRI)")
    a.set_ylabel("miRNA × chemical tests")
    a.legend(loc="upper left", handlelength=1.0, borderaxespad=0.2)
    b = axes_mm(fig, w, h, 104, 15, w - 108, 44)
    chems = sorted(t["chemical"].unique())
    productive_fraction(b, t, "chemical", chems, [CHEM[c] for c in chems], BLUE)
    plt.setp(b.get_xticklabels(), rotation=40, ha="right")
    figure_label(fig, "A", w, h, left_mm=2, top_mm=LABEL_TOP_MM, fontsize=10)
    figure_label(fig, "B", w, h, left_mm=88, top_mm=LABEL_TOP_MM, fontsize=10)
    return fig


def fig3():
    spec = pd.read_csv(require("results/specificity.tsv"), sep="\t")
    fig, w, h = canvas("double", 72)
    a = axes_mm(fig, w, h, 17, 15, 55, 44)
    n_chem = 9
    counts = spec["n_chemicals_productive"].value_counts().reindex(range(n_chem + 1),
                                                                     fill_value=0)
    a.bar(counts.index, counts.values, color=BLUE, width=0.7)
    for x, v in counts.items():
        if v:
            a.text(x, v + 3.5, str(int(v)), ha="center", va="bottom")
    a.set_xticks(range(n_chem + 1))
    a.set_xlim(-0.6, n_chem + 0.6)
    a.set_ylim(0, counts.max() * 1.15)
    a.set_xlabel("Chemicals in which the miRNA is productive")
    a.set_ylabel("Evaluable miRNAs")
    b = axes_mm(fig, w, h, 122, 15, w - 126, 44)
    top = spec.sort_values("benzene_contrast", ascending=False).head(12)
    y = np.arange(len(top))[::-1]
    b.barh(y, top["benzene_contrast"], color=VERM, height=0.7)
    b.set_yticks(y)
    b.set_yticklabels([m.replace("mmu-", "") for m in top["miRNA"]])
    b.set_xlabel("Benzene CRI contrast\n(benzene − mean of other chemicals)", linespacing=LS)
    b.set_ylim(-0.6, len(top) - 0.4)
    figure_label(fig, "A", w, h, left_mm=2, top_mm=LABEL_TOP_MM, fontsize=10)
    figure_label(fig, "B", w, h, left_mm=88, top_mm=LABEL_TOP_MM, fontsize=10)
    return fig


def fig4():
    pc = pd.read_csv(require("results/posctrl_coherence.tsv"), sep="\t")
    t = pc[pc["class"].isin(TESTED)]
    fig, w, h = canvas("double", 150)
    # (A) mimic AUC by group
    a = axes_mm(fig, w, h, 20, 15, 56, 42)
    groups = ["transfected", "random"]
    box_and_points(a, [t.loc[t.group == g, "auc"].values for g in groups], BLUE)
    a.set_xticks([1, 2])
    a.set_xticklabels(["Transfected\nmiRNAs", "Random\nmiRNAs"])
    a.set_ylabel("Target-vs-background\nrepression AUC", linespacing=LS)
    # (B) productive fraction by group
    b = axes_mm(fig, w, h, 112, 15, 58, 42)
    labs = []
    for g in groups:
        sub = t[t.group == g]
        labs.append(f"{g.capitalize()}\n({int((sub['class'] == 'productive').sum())}/{len(sub)})")
    productive_fraction(b, t, "group", groups, labs, BLUE)
    # (C) Dicer1 knockout at endogenous abundance
    d = pd.read_csv(require("results/revision_r1/physctrl_coherence.tsv"), sep="\t")
    d = d[d["class"].isin(TESTED)].copy()
    d["x"] = np.log10(d["mesc_ctrl_cpm"] + 1.0)
    d["y"] = 0.5 - d["auc"]
    c = axes_mm(fig, w, h, 20, 83, 92, 50)
    c.axhline(0, color=GREY, ls="--", lw=0.8, zorder=1)
    handles = []
    for s, col, name in (("not expressed", GREY, "Not expressed"),
                         ("expressed", BLUE, "Expressed")):
        g = d[d["set"] == s]
        prod = g["class"] == "productive"
        c.scatter(g.loc[~prod, "x"], g.loc[~prod, "y"], s=9, facecolors="none",
                  edgecolors=col, linewidths=0.6, zorder=2)
        c.scatter(g.loc[prod, "x"], g.loc[prod, "y"], s=9, color=col, linewidths=0,
                  zorder=3)
        handles += [Line2D([], [], ls="", marker="o", ms=3.5, color=col,
                           label=f"{name}: productive ({int(prod.sum())}/{len(g)})"),
                    Line2D([], [], ls="", marker="o", ms=3.5, mfc="none", mec=col,
                           mew=0.6, label=f"{name}: not productive")]
    cl = d[d["miRNA"].str.match(r"mmu-miR-29[0-5][ab]?-3p$")]
    c.scatter(cl["x"], cl["y"], s=34, facecolors="none", edgecolors=VERM, linewidths=0.9,
              zorder=4)
    handles.append(Line2D([], [], ls="", marker="o", ms=6, mfc="none", mec=VERM, mew=0.9,
                          label=f"miR-290–295 cluster, 3p (n = {len(cl)})"))
    c.set_xlabel("Abundance in vehicle-control mESCs, log10(CPM + 1)")
    c.set_ylabel("De-repression of predicted\ntargets in Dicer1 KO (0.5 − AUC)", linespacing=LS)
    c.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1.0),
             borderaxespad=0, handletextpad=0.3)
    figure_label(fig, "A", w, h, left_mm=2, top_mm=LABEL_TOP_MM, fontsize=10)
    figure_label(fig, "B", w, h, left_mm=94, top_mm=LABEL_TOP_MM, fontsize=10)
    figure_label(fig, "C", w, h, left_mm=2, top_mm=76, fontsize=10)
    return fig


def fig5():
    df = pd.read_csv(require("results/revision_r1/power_coherent.tsv"), sep="\t").sort_values("delta")
    fig, w, h = canvas("single", 72)
    a = axes_mm(fig, w, h, 15, 11, 60, 45)
    y = df["pct_productive_seed0"]
    a.errorbar(df["delta"], y, yerr=[y - df["pct_productive_min_20_seeds"],
                                     df["pct_productive_max_20_seeds"] - y],
               fmt="o-", color=BLUE, ms=3.5, capsize=2, elinewidth=0.8)
    a.set_xlabel("Injected coherent shift δ (log2 units)")
    a.set_ylabel("Tests called productive (%)", color=BLUE)
    a.tick_params(axis="y", labelcolor=BLUE)
    a.set_ylim(-3, 103)
    r = a.twinx()
    r.spines["right"].set_visible(True)
    r.plot(df["delta"], df["mean_coherence_auc_seed0"], "s--", color=VERM, ms=3.5)
    r.axhline(0.5, color=GREY, ls=":", lw=0.8)
    r.set_ylabel("Mean coherence-signed AUC", color=VERM, labelpad=6)
    r.tick_params(axis="y", labelcolor=VERM)
    r.set_ylim(0.45, 1.0)
    return fig


ORDER = [("Fig1", fig1), ("Fig2", fig2), ("Fig3", fig3), ("Fig4", fig4), ("Fig5", fig5)]


def main():
    os.makedirs(OUT, exist_ok=True)
    for stem, make in ORDER:
        fig = make()
        for p in savefig(fig, os.path.join(OUT, stem), dpi=600):
            print("wrote", p)
        plt.close(fig)


if __name__ == "__main__":
    main()
