"""Fig 5: spike-in power analysis. A known coherent repression of size delta is
injected into the real target sets and the test re-run; the fraction called
productive (and the mean AUC) rises steeply with delta, while delta=0 returns ~0.
Shows the 0/225 real-data null is an absence of signal, not low power. Rendered
from results/power.tsv (no hardcoded numbers)."""
import matplotlib.pyplot as plt
import pandas as pd

from _common import require, save

POWER = "results/power.tsv"


def main():
    df = pd.read_csv(require(POWER), sep="\t").sort_values("delta")
    fig, ax1 = plt.subplots(figsize=(6, 4))
    ax1.plot(df["delta"], df["pct_productive"], "o-", color="seagreen", label="% productive")
    ax1.set_xlabel("injected coherent shift δ (log2 units)")
    ax1.set_ylabel("% of tests called productive", color="seagreen")
    ax1.tick_params(axis="y", labelcolor="seagreen")
    ax1.set_ylim(-3, 103)
    ax2 = ax1.twinx()
    ax2.plot(df["delta"], df["mean_auc"], "s--", color="steelblue", label="mean AUC")
    ax2.axhline(0.5, color="grey", ls=":", lw=1)
    ax2.set_ylabel("mean AUC", color="steelblue")
    ax2.tick_params(axis="y", labelcolor="steelblue")
    ax2.set_ylim(0.45, 1.0)
    ax1.set_title("Spike-in power: the test detects modest coherent shifts")
    fig.tight_layout()
    save(fig, "fig5_power")


if __name__ == "__main__":
    main()
