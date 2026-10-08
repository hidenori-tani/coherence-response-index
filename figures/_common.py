"""Shared figure helpers. All figures are rendered from results/*.tsv at run time
(no hardcoded numbers), and saved as both vector PDF and 300-dpi PNG."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

COHERENCE = "results/coherence.tsv"
SPECIFICITY = "results/specificity.tsv"


def require(path):
    if not os.path.exists(path):
        raise SystemExit(f"missing {path}; run scripts/run_all.py first")
    return path


def save(fig, stem):
    os.makedirs("figures/out", exist_ok=True)
    fig.savefig(f"figures/out/{stem}.pdf")
    fig.savefig(f"figures/out/{stem}.png", dpi=300)
    plt.close(fig)
    print(f"wrote figures/out/{stem}.pdf and .png")


def wilson_ci(k, n, z=1.96):
    """Wilson score interval for a proportion (for the productive-fraction bars)."""
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    denom = 1 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    half = z * ((p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5) / denom
    return (max(0.0, center - half), min(1.0, center + half))
