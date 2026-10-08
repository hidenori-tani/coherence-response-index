import random, bisect

def matched_background(targets, covariate, n_per_target=20, seed=0, window=50):
    """For each target, sample n non-target genes with nearest covariate value
    (expression/length/GC): draws come from the +/-`window` genes nearest to the
    target in the covariate-sorted non-target list. Returns a set of background
    gene IDs."""
    rng = random.Random(seed)
    items = sorted((v, g) for g, v in covariate.items() if g not in targets)
    vals = [v for v, _ in items]
    bg = set()
    # iterate targets in a deterministic order so the background (and every
    # downstream statistic) is bit-reproducible given the seed, independent of
    # Python set hash order / PYTHONHASHSEED.
    for t in sorted(targets):
        if t not in covariate:
            continue
        i = bisect.bisect_left(vals, covariate[t])
        lo, hi = max(0, i - window), min(len(items), i + window)
        pool = [g for _, g in items[lo:hi] if g not in bg]
        rng.shuffle(pool)
        bg.update(pool[:n_per_target])
    return bg
