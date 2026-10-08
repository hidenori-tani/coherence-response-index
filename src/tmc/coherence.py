"""Core coherence test: do a miRNA's predicted targets shift in the expected
direction (vs a matched non-target background), and how coherent is the response?

Coherence Response Index (CRI), per (miRNA m, chemical c):
- L = miRNA log2FC (L>0 = up-regulated, L<0 = down-regulated).
- D = target_median_diff = median(target log2FC) - median(background log2FC).
      D<0 = targets more repressed than background.
- auc = P(target log2FC < background log2FC); auc>0.5 = targets more repressed.

An UP miRNA is coherent when its targets are repressed (D<0, auc>0.5); a DOWN
miRNA is coherent when its targets are de-repressed (D>0, auc<0.5). CRI is
positive for coherent responses, negative for anti-coherent ones, and its
magnitude scales with |L| and the target-vs-background separation |auc-0.5|.
"""
import statistics
from scipy.stats import mannwhitneyu


def directional_shift(target_log2fc, background_log2fc):
    """Compare target vs background log2FC. Returns median_diff (target-bg),
    AUC = P(target < background) for repression, and Mann-Whitney p (two-sided)."""
    t, b = list(target_log2fc), list(background_log2fc)
    if not t or not b:
        return {"median_diff": float("nan"), "auc": float("nan"),
                "p": float("nan"), "n": len(t)}
    # AUC that targets are MORE repressed (smaller) than background
    less = sum(1 for x in t for y in b if x < y)
    ties = sum(1 for x in t for y in b if x == y)
    auc = (less + 0.5 * ties) / (len(t) * len(b))
    u, p = mannwhitneyu(t, b, alternative="two-sided")
    return {"median_diff": statistics.median(t) - statistics.median(b),
            "auc": auc, "p": float(p), "n": len(t)}


def cri(mirna_log2fc, target_median_diff, auc):
    """Coherence Response Index. Sign = coherence direction, magnitude =
    |miRNA log2FC| * (2*|auc-0.5|). An up miRNA expects negative target_median_diff;
    a down miRNA expects positive. Coherent -> positive CRI; anti-coherent -> negative."""
    expected_target_sign = -1.0 if mirna_log2fc >= 0 else 1.0
    coherence_sign = 1.0 if (target_median_diff * expected_target_sign) > 0 else -1.0
    separation = abs(auc - 0.5) * 2.0  # in [0, 1]
    return coherence_sign * abs(mirna_log2fc) * separation


def classify_response(p, auc, mirna_log2fc, p_thresh=0.05, auc_margin=0.1):
    """A responsive miRNA is 'productive' when its targets are significantly and
    substantially shifted vs background; otherwise 'passenger'. (mirna_log2fc is
    part of the contract for callers that pre-filter responsive miRNAs upstream.)"""
    if p < p_thresh and abs(auc - 0.5) >= auc_margin:
        return "productive"
    return "passenger"
