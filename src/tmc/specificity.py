"""Generic-vs-chemical-specific structure of coherent (productive) miRNA responses,
and the benzene (carcinogen) contrast."""


def generic_fraction(by_chemical, k):
    """by_chemical: {chemical: set(productive miRNAs)}. Returns {miRNA: bool},
    True when the miRNA is productive in at least k chemicals (a 'generic' response)."""
    counts = {}
    for productive in by_chemical.values():
        for m in productive:
            counts[m] = counts.get(m, 0) + 1
    return {m: (c >= k) for m, c in counts.items()}


def benzene_contrast(cri_by_chemical, benzene_key="benzene"):
    """cri_by_chemical: {chemical: {miRNA: CRI}}. Returns {miRNA: CRI_benzene -
    mean(CRI in the other chemicals)} over every miRNA seen in any chemical
    (missing CRI treated as 0)."""
    mirnas = set()
    for d in cri_by_chemical.values():
        mirnas |= set(d)
    benzene = cri_by_chemical.get(benzene_key, {})
    others = [c for c in cri_by_chemical if c != benzene_key]
    out = {}
    for m in mirnas:
        b = benzene.get(m, 0.0)
        other_vals = [cri_by_chemical[c].get(m, 0.0) for c in others]
        mean_other = sum(other_vals) / len(other_vals) if other_vals else 0.0
        out[m] = b - mean_other
    return out
