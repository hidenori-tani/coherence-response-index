"""Orchestration: turn the two log2FC matrices + predicted targets into the
coherence table and the specificity summary. Pure functions over DataFrames so
they can be unit-tested with synthetic data (no files, no network)."""
import pandas as pd
from statsmodels.stats.multitest import multipletests

from tmc.targets import target_set, harmonize_mirna_names
from tmc.nullmodel import matched_background
from tmc.coherence import directional_shift, cri, classify_response


def _harmonize_target_dict(d):
    """Collapse target-dict keys with harmonize_mirna_names, unioning sets."""
    out = {}
    for k, genes in d.items():
        hk = harmonize_mirna_names(k)
        out.setdefault(hk, set()).update(genes)
    return out


def run_coherence(mirna_l2fc, mrna_l2fc, targetscan, mirdb, baseline_expr,
                  responsive_thresh=1.0, target_mode="union", min_targets=5,
                  n_per_target=20, seed=0, fdr_thresh=0.05, auc_margin=0.1,
                  window=50):
    """For each chemical x responsive miRNA, test whether predicted targets shift
    against a matched non-target background. Returns a tidy DataFrame with one row
    per (chemical, miRNA): n_targets, median_diff, auc, p, fdr, cri, class.

    - mirna_l2fc / mrna_l2fc: DataFrames (index=feature, columns=chemicals).
    - targetscan / mirdb: dict miRNA -> set(gene) (any naming; harmonized here).
    - baseline_expr: dict gene -> covariate value for matched-background selection.
    """
    ts = _harmonize_target_dict(targetscan)
    mb = _harmonize_target_dict(mirdb)
    genes_present = set(mrna_l2fc.index)
    # covariate restricted to genes we can actually read a log2FC for
    covar = {g: v for g, v in baseline_expr.items() if g in genes_present}

    rows = []
    for chem in mirna_l2fc.columns:
        if chem not in mrna_l2fc.columns:
            continue
        mrna_col = mrna_l2fc[chem]
        for mirna, l in mirna_l2fc[chem].items():
            if pd.isna(l) or abs(l) < responsive_thresh:
                continue
            hk = harmonize_mirna_names(mirna)
            targets = (target_set(hk, ts, mb, mode=target_mode) & genes_present)
            n_targets = len(targets)
            row = {"chemical": chem, "miRNA": mirna, "mirna_log2fc": l,
                   "n_targets": n_targets, "median_diff": float("nan"),
                   "auc": float("nan"), "p": float("nan"), "cri": float("nan")}
            if n_targets >= min_targets:
                target_vals = [mrna_col[g] for g in targets]
                bg = matched_background(targets, covar, n_per_target=n_per_target,
                                        seed=seed, window=window)
                bg_vals = [mrna_col[g] for g in bg if g in mrna_col.index]
                shift = directional_shift(target_vals, bg_vals)
                row.update(median_diff=shift["median_diff"], auc=shift["auc"], p=shift["p"])
                row["cri"] = cri(l, shift["median_diff"], shift["auc"])
            rows.append(row)

    df = pd.DataFrame(rows)
    if df.empty:
        df["fdr"] = []
        df["class"] = []
        return df

    # BH-FDR across all tested rows (finite p)
    df["fdr"] = float("nan")
    tested = df["p"].notna()
    if tested.any():
        df.loc[tested, "fdr"] = multipletests(df.loc[tested, "p"], method="fdr_bh")[1]

    def _cls(r):
        if pd.isna(r["fdr"]):
            return "untested"
        base = classify_response(r["fdr"], r["auc"], r["mirna_log2fc"],
                                 p_thresh=fdr_thresh, auc_margin=auc_margin)
        if base != "productive":
            return "passenger"
        # significant shift: productive only if in the COHERENT direction (CRI>0),
        # otherwise anti-coherent (e.g. down miRNA whose targets are also repressed)
        return "productive" if r["cri"] > 0 else "anti_coherent"
    df["class"] = df.apply(_cls, axis=1)
    return df[["chemical", "miRNA", "mirna_log2fc", "n_targets",
               "median_diff", "auc", "p", "fdr", "cri", "class"]]


def summarize_specificity(coherence_df, k, benzene_key="benzene"):
    """From the coherence table, per miRNA: how many chemicals it is productive in,
    whether that is 'generic' (>=k chemicals), and the benzene CRI contrast."""
    from tmc.specificity import generic_fraction, benzene_contrast

    prod = coherence_df[coherence_df["class"] == "productive"]
    by_chem = {c: set(g["miRNA"]) for c, g in prod.groupby("chemical")}
    generic = generic_fraction(by_chem, k)

    cri_by_chem = {}
    for c, g in coherence_df.groupby("chemical"):
        cri_by_chem[c] = {row["miRNA"]: row["cri"] for _, row in g.iterrows()
                          if pd.notna(row["cri"])}
    contrast = benzene_contrast(cri_by_chem, benzene_key=benzene_key)

    n_prod = {}
    for m in set(generic) | set(contrast):
        n_prod[m] = sum(1 for s in by_chem.values() if m in s)

    rows = [{"miRNA": m, "n_chemicals_productive": n_prod.get(m, 0),
             "generic": bool(generic.get(m, False)),
             "benzene_contrast": contrast.get(m, float("nan"))}
            for m in sorted(set(generic) | set(contrast))]
    return pd.DataFrame(rows)
