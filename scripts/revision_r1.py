"""Analyses added for the CBI major revision (CHEMBIOINT-D-26-02865, reviewers R1/R3).

Every analysis reads only the gene-level result tables in results/ (the raw reads
are not needed) and writes to results/revision_r1/. The primary analysis is first
re-run in memory and must reproduce results/coherence.tsv exactly before any
variant is trusted.

  python scripts/revision_r1.py all        # everything below
  python scripts/revision_r1.py <name>     # one of: baseline window weight
                                            #   replicates mir200 offon targetscan

R1-2  replicates  : each replicate half (control rep k + chemical rep k) analysed alone
R1-3  mir200      : the single productive call (DEHP / miR-200a-5p) examined
R1-6  weight      : CRI without the |miRNA log2FC| weight
R1-7  window      : background window (+/-10..250) and draws per target (10/20/40)
R1-9  targetscan  : TargetScan with the miRBase-ID -> family map (was: exact family-name lookup)
R3-6  offon       : genuine off-to-on / on-to-off genes under the pseudocount estimator
R1-2  replicate_noise : replicate-level numbers of Results 3.1 / 3.6 (after 'replicates')
"""
import csv
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr, pearsonr
from statsmodels.stats.multitest import multipletests

from tmc.config import load_accessions, CHEMICALS, CONTROL
from tmc.io_expr import load_matrix, compute_log2fc_robust
from tmc.targets import load_mirdb
from tmc.pipeline import run_coherence

OUT = Path("results/revision_r1")
MIN_CTRL_TPM, MIN_CTRL_CPM = 1.0, 10.0
PRIMARY = dict(responsive_thresh=1.0, target_mode="union", min_targets=5,
               n_per_target=20, seed=0)


# ---------------------------------------------------------------- inputs
def runs(layer, chemical, replicate=None):
    acc = [r for r in load_accessions("config/accessions.tsv") if r["layer"] == layer]
    return [r["run_accession"] for r in acc if r["chemical"] == chemical
            and (replicate is None or str(r["replicate"]) == str(replicate))]


def matrices():
    tpm = load_matrix("results/gene_tpm.tsv")
    cnt = load_matrix("results/mirna_counts.tsv")
    cpm = cnt / cnt.sum() * 1e6
    return tpm, cpm


def expressed(tpm, cpm):
    g = tpm[runs("mrna", CONTROL)].mean(axis=1)
    m = cpm[runs("mirna", CONTROL)].mean(axis=1)
    return set(g[g >= MIN_CTRL_TPM].index), set(m[m >= MIN_CTRL_CPM].index), g


def baseline_covariate(g_ctrl_mean, exg):
    return {k: float(np.log10(v + 1.0)) for k, v in g_ctrl_mean.items() if k in exg}


def primary_inputs():
    tpm, cpm = matrices()
    exg, exm, g = expressed(tpm, cpm)
    mirna = load_matrix("results/mirna_log2fc_robust.tsv")
    mrna = load_matrix("results/mrna_log2fc_robust.tsv")
    return (mirna[mirna.index.isin(exm)], mrna[mrna.index.isin(exg)],
            load_mirdb("config/mirdb_mmu.tsv", score_min=80), baseline_covariate(g, exg))


def summary(df, label):
    t = df[df["class"] != "untested"]
    return {"variant": label, "tested": len(t),
            "productive": int((t["class"] == "productive").sum()),
            "anti_coherent": int((t["class"] == "anti_coherent").sum()),
            "auc_mean": round(float(t["auc"].mean()), 3),
            "auc_median": round(float(t["auc"].median()), 3),
            "auc_sd": round(float(t["auc"].std()), 3),
            "frac_auc_gt_half": round(float((t["auc"] > 0.5).mean()), 3)}


def write(df, name):
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / name, sep="\t", index=False)


# ---------------------------------------------------------------- analyses
def baseline():
    """The primary analysis, re-run: must equal results/coherence.tsv."""
    mi, mr, mb, base = primary_inputs()
    df = run_coherence(mi, mr, {}, mb, base, **PRIMARY)
    ref = pd.read_csv("results/coherence.tsv", sep="\t")
    same = (len(df) == len(ref)
            and (df["class"].values == ref["class"].values).all()
            and np.allclose(df["auc"].fillna(-1), ref["auc"].fillna(-1), atol=1e-12, rtol=0)
            and np.allclose(df["fdr"].fillna(-1), ref["fdr"].fillna(-1), atol=1e-12, rtol=0))
    if not same:
        raise SystemExit("ABORT: re-run does not reproduce results/coherence.tsv")
    print(f"[baseline] reproduces results/coherence.tsv: {summary(df, 'primary')}")
    return df


def window():
    mi, mr, mb, base = primary_inputs()
    rows = []
    for w in (10, 25, 50, 100, 250):
        df = run_coherence(mi, mr, {}, mb, base, **{**PRIMARY, "window": w})
        rows.append({**summary(df, f"window=+/-{w}, draws=20"), "window": w, "draws": 20})
    for n in (10, 40):
        df = run_coherence(mi, mr, {}, mb, base, **{**PRIMARY, "n_per_target": n})
        rows.append({**summary(df, f"window=+/-50, draws={n}"), "window": 50, "draws": n})
    for s in range(1, 11):
        df = run_coherence(mi, mr, {}, mb, base, **{**PRIMARY, "seed": s})
        rows.append({**summary(df, f"window=+/-50, draws=20, seed={s}"), "window": 50,
                     "draws": 20, "seed": s})
    out = pd.DataFrame(rows)
    write(out, "window_sensitivity.tsv")
    print("[window]\n" + out[["variant", "tested", "productive", "anti_coherent",
                              "auc_mean", "auc_sd"]].to_string(index=False))


def weight():
    """R1-6: the classification never reads |miRNA log2FC| (only FDR, |AUC-0.5| and
    the sign of the target shift), so dropping the weight cannot change a call. Here
    we measure how much the weight changes the *ranking* the CRI is used for."""
    df = pd.read_csv("results/coherence.tsv", sep="\t")
    t = df[df["class"] != "untested"].copy()
    t["cri_unweighted"] = np.sign(t["cri"]) * 2 * (t["auc"] - 0.5).abs()
    rho = spearmanr(t["cri"], t["cri_unweighted"]).statistic
    top_w = set(t.nlargest(12, "cri")["miRNA"] + "|" + t.nlargest(12, "cri")["chemical"])
    top_u = set(t.nlargest(12, "cri_unweighted")["miRNA"] + "|"
                + t.nlargest(12, "cri_unweighted")["chemical"])
    res = pd.DataFrame([{"n_tests": len(t), "spearman_weighted_vs_unweighted": round(rho, 3),
                         "top12_overlap": len(top_w & top_u),
                         "classification_changes": 0}])
    write(res, "weight.tsv")
    print(f"[weight] {res.to_dict('records')[0]}")


def _half_log2fc(mat, layer, rep):
    ctrl = mat[runs(layer, CONTROL, rep)].mean(axis=1).to_dict()
    return pd.DataFrame({c: compute_log2fc_robust(mat[runs(layer, c, rep)].mean(axis=1).to_dict(),
                                                  ctrl) for c in CHEMICALS})


def replicates():
    """R1-2: analyse each replicate half on its own (control rep k vs chemical rep k,
    both layers from the same RNA preparation), on the primary feature universe."""
    tpm, cpm = matrices()
    exg, exm, g = expressed(tpm, cpm)
    base = baseline_covariate(g, exg)
    mb = load_mirdb("config/mirdb_mmu.tsv", score_min=80)
    halves, rows = {}, []
    for rep in (1, 2):
        mi = _half_log2fc(cpm, "mirna", rep); mi = mi[mi.index.isin(exm)]
        mr = _half_log2fc(tpm, "mrna", rep); mr = mr[mr.index.isin(exg)]
        df = run_coherence(mi, mr, {}, mb, base, **PRIMARY)
        halves[rep] = (mi, mr, df)
        write(df.sort_values(["chemical", "miRNA"]), f"replicates_half{rep}_coherence.tsv")
        write(mi.reset_index().rename(columns={"index": "miRNA"}),
              f"replicates_half{rep}_mirna_log2fc.tsv")
        rows.append(summary(df, f"replicate half {rep} only"))
    # per-feature noise between halves (what n=2 means for single genes)
    for layer, i in (("mRNA", 1), ("miRNA", 0)):
        r = [pearsonr(halves[1][i][c], halves[2][i][c]).statistic for c in CHEMICALS]
        rows.append({"variant": f"{layer} log2FC, half1 vs half2 (Pearson r per chemical)",
                     "r_median": round(float(np.median(r)), 3),
                     "r_min": round(float(min(r)), 3), "r_max": round(float(max(r)), 3)})
    # does any per-test signal replicate across halves?
    a = halves[1][2]; b = halves[2][2]
    key = ["chemical", "miRNA"]
    m = a[a["class"] != "untested"].merge(b[b["class"] != "untested"], on=key,
                                          suffixes=("_1", "_2"))
    rho = spearmanr(m["auc_1"], m["auc_2"])
    both = int(((m["class_1"] == "productive") & (m["class_2"] == "productive")).sum())
    rows.append({"variant": "tests evaluated in both halves", "tested": len(m),
                 "productive_in_both": both,
                 "auc_spearman_half1_vs_half2": round(float(rho.statistic), 3),
                 "auc_spearman_p": float(rho.pvalue)})
    out = pd.DataFrame(rows)
    write(out, "replicates.tsv")
    write(m.sort_values(key), "replicates_paired_tests.tsv")   # row order is set-order
    print("[replicates]\n" + out.to_string(index=False))


def _repro_row(label, a, b, thresh=1.0):
    """Agreement of two independent log2FC vectors over the same features."""
    r = pearsonr(a, b).statistic
    sel = a.abs() >= thresh
    same_sign = float((np.sign(a[sel]) == np.sign(b[sel])).mean()) if sel.any() else np.nan
    both = float(((b[sel].abs() >= thresh) & (np.sign(a[sel]) == np.sign(b[sel]))).mean()) \
        if sel.any() else np.nan
    return {"set": label, "n_features": len(a), "pearson_r": round(float(r), 3),
            "n_called_in_A": int(sel.sum()),
            "same_sign_in_B_pct": round(100 * same_sign, 1),
            "called_again_in_B_pct": round(100 * both, 1)}


def reproducibility():
    """R1-2: how reproducible are the per-chemical fold changes themselves at n=2?
    Split A = control rep 1 vs chemical rep 1; split B = control rep 2 vs chemical rep 2
    (independent samples). 'Called' = robust |log2FC| >= 1 in split A. The same metric is
    applied to the positive-control transfection (mimic rep i vs control rep i), where a
    real perturbation is known to exist, as the comparator."""
    tpm, cpm = matrices()
    exg, exm, g = expressed(tpm, cpm)
    rows = []
    for layer, mat, feats in (("mRNA", tpm, exg), ("miRNA", cpm, exm)):
        A = _half_log2fc(mat, layer.lower(), 1); B = _half_log2fc(mat, layer.lower(), 2)
        A, B = A[A.index.isin(feats)], B[B.index.isin(feats)]
        for c in CHEMICALS:
            rows.append({"dataset": "toxicogenomic", "layer": layer, "condition": c,
                         **_repro_row(f"{layer} {c}", A[c], B[c])})
        if layer == "mRNA":  # by expression tertile (control mean TPM)
            lvl = g[g.index.isin(feats)]
            q = pd.qcut(lvl, 3, labels=["low", "mid", "high"])
            for t in ("low", "mid", "high"):
                idx = q[q == t].index
                rows.append({"dataset": "toxicogenomic", "layer": layer,
                             "condition": f"all chemicals pooled, {t} tertile",
                             **_repro_row(f"{layer} tertile {t}",
                                          A.loc[idx].stack(), B.loc[idx].stack())})
    # the calls the paper actually uses: responsive = |mean-of-replicates log2FC| >= 1
    mi_mean = load_matrix("results/mirna_log2fc_robust.tsv")
    mi_mean = mi_mean[mi_mean.index.isin(exm)]
    A = _half_log2fc(cpm, "mirna", 1).loc[mi_mean.index]
    B = _half_log2fc(cpm, "mirna", 2).loc[mi_mean.index]
    m, a, b = mi_mean.stack(), A.stack(), B.stack()
    sel = m.abs() >= 1.0
    agree = (np.sign(a[sel]) == np.sign(m[sel])) & (np.sign(b[sel]) == np.sign(m[sel]))
    agree_half = agree & (a[sel].abs() >= 0.5) & (b[sel].abs() >= 0.5)
    rows.append({"dataset": "toxicogenomic", "layer": "miRNA",
                 "condition": "responsive calls used in the paper (|mean log2FC| >= 1)",
                 "set": "primary responsive miRNA x chemical", "n_features": int(sel.sum()),
                 "both_halves_same_sign_pct": round(100 * float(agree.mean()), 1),
                 "both_halves_same_sign_and_abs_ge_0.5_pct": round(100 * float(agree_half.mean()), 1)})
    # comparator: positive-control transfection, mRNA, rep i vs rep i
    acc = list(csv.DictReader(open("config/posctrl_accessions.tsv"), delimiter="\t"))
    ptpm = load_matrix("results/posctrl_gene_tpm.tsv")
    for cell in sorted({r["cell"] for r in acc}):
        ctl = [r["run_accession"] for r in acc if r["cell"] == cell and r["condition"] == "control"]
        mim = [r["run_accession"] for r in acc if r["cell"] == cell and r["condition"] != "control"]
        base = ptpm[ctl].mean(axis=1)
        feats = base[base >= MIN_CTRL_TPM].index
        l2 = lambda i: pd.Series(compute_log2fc_robust(ptpm[mim[i]].to_dict(),
                                                       ptpm[ctl[i]].to_dict())).loc[feats]
        rows.append({"dataset": "positive control (mimic)", "layer": "mRNA", "condition": cell,
                     **_repro_row(f"posctrl {cell}", l2(0), l2(1))})
    out = pd.DataFrame(rows)
    write(out, "reproducibility.tsv")
    print("[reproducibility]\n" + out.drop(columns=["set"]).to_string(index=False))


def _half_text(chem, mir):
    """S2a: what happened to one test in each replicate half, stated explicitly."""
    out = []
    for k in (1, 2):
        l2 = pd.read_csv(OUT / f"replicates_half{k}_mirna_log2fc.tsv", sep="\t").set_index("miRNA").loc[mir, chem]
        c = pd.read_csv(OUT / f"replicates_half{k}_coherence.tsv", sep="\t")
        c = c[(c["chemical"] == chem) & (c["miRNA"] == mir)]
        if len(c) and c.iloc[0]["class"] != "untested":
            r = c.iloc[0]
            out.append(f"half {k}: {r['class']} (miRNA log2FC {l2:.2f}; AUC {r['auc']:.2f}, FDR {r['fdr']:.2f})")
        else:
            out.append(f"half {k}: not tested (miRNA log2FC {l2:.2f}, below the responsiveness threshold of 1)")
    return "; ".join(out)


def mir200(n_draws=2000):
    """R1-3: the single productive call. (i) the call itself; (ii) Storey pi0 over all
    225 p-values; (iii) an empirical per-test null: random expression-matched non-target
    gene sets of the same size in the same chemical; (iv) persistence across background
    seeds and replicate halves; (v) the other miR-200 family members in DEHP."""
    mi, mr, mb, base = primary_inputs()
    ref = pd.read_csv("results/coherence.tsv", sep="\t")
    t = ref[ref["class"] != "untested"]
    hit = ref[ref["class"] == "productive"]
    assert len(hit) == 1, hit
    h = hit.iloc[0]
    chem, mir = h["chemical"], h["miRNA"]
    lam = 0.5
    pi0 = min(1.0, float((t["p"] > lam).sum()) / ((1 - lam) * len(t)))
    # (iii) empirical null for this chemical and target-set size
    from tmc.targets import harmonize_mirna_names
    mbh = {}
    for k, v in mb.items():
        mbh.setdefault(harmonize_mirna_names(k), set()).update(v)
    targets = mbh[harmonize_mirna_names(mir)] & set(mr.index)
    col = mr[chem]
    covar = {k: v for k, v in base.items() if k in col.index}
    items = sorted((v, k) for k, v in covar.items() if k not in targets)
    vals = np.array([v for v, _ in items]); names = [k for _, k in items]
    rng = np.random.default_rng(0)
    tgt_sorted = sorted(targets)
    obs_auc, obs_p = float(h["auc"]), float(h["p"])
    from tmc.nullmodel import matched_background
    bg = matched_background(targets, covar, n_per_target=20, seed=0)
    bg_vals = col[list(bg)].values
    null_p, null_dev = [], []
    for _ in range(n_draws):
        # a random non-target gene set with the same expression profile as the targets
        pick = []
        for gname in tgt_sorted:
            i = int(np.searchsorted(vals, covar[gname]))
            lo, hi = max(0, i - 50), min(len(names), i + 50)
            pick.append(names[int(rng.integers(lo, hi))])
        pv = col[pick].values
        u, p = mannwhitneyu(pv, bg_vals, alternative="two-sided")
        null_p.append(p)
        null_dev.append(abs((1 - u / (len(pv) * len(bg_vals))) - 0.5))
    emp_p = (1 + sum(p <= obs_p for p in null_p)) / (n_draws + 1)
    emp_dev = (1 + sum(d >= abs(obs_auc - 0.5) for d in null_dev)) / (n_draws + 1)
    # (iv) background seeds
    seeds = []
    for s in range(0, 20):
        df = run_coherence(mi, mr, {}, mb, base, **{**PRIMARY, "seed": s})
        r = df[(df["chemical"] == chem) & (df["miRNA"] == mir)].iloc[0]
        seeds.append((r["class"], r["fdr"], r["auc"]))
    rep = pd.read_csv(OUT / "replicates_paired_tests.tsv", sep="\t") \
        if (OUT / "replicates_paired_tests.tsv").exists() else pd.DataFrame()
    rr = rep[(rep.get("chemical") == chem) & (rep.get("miRNA") == mir)] if len(rep) else rep
    # (v) miR-200 family in DEHP (the 5p arm is the minor arm of the miR-200a hairpin)
    tpm, cpm = matrices()
    fam = [x for x in cpm.index if any(s in x for s in
           ("miR-200a-", "miR-200b-", "miR-200c-", "miR-141-", "miR-429"))]
    ctrl = cpm.loc[fam, runs("mirna", CONTROL)].mean(axis=1)
    l2 = pd.read_csv("results/mirna_log2fc_robust.tsv", sep="\t", index_col=0).loc[fam, chem]
    famtab = pd.DataFrame({"miRNA": fam, "ctrl_CPM": ctrl.round(1).values,
                           f"{chem}_robust_log2FC": l2.round(2).values})
    fam_rows = ref[(ref["chemical"] == chem) & (ref["miRNA"].isin(fam))][
        ["miRNA", "mirna_log2fc", "n_targets", "auc", "fdr", "class"]]
    res = pd.DataFrame([{
        "chemical": chem, "miRNA": mir, "mirna_log2fc": round(float(h["mirna_log2fc"]), 3),
        "n_targets": int(h["n_targets"]), "auc": round(obs_auc, 3), "p": obs_p,
        "fdr": round(float(h["fdr"]), 4), "storey_pi0_all_tests": round(pi0, 3),
        "empirical_p_vs_random_matched_sets": round(emp_p, 4),
        "empirical_p_auc_deviation": round(emp_dev, 4), "n_draws": n_draws,
        "productive_in_seeds_0_19": sum(c == "productive" for c, _, _ in seeds),
        "fdr_range_seeds": f"{min(f for _, f, _ in seeds):.3f}-{max(f for _, f, _ in seeds):.3f}",
        "auc_range_seeds": f"{min(a for _, _, a in seeds):.3f}-{max(a for _, _, a in seeds):.3f}",
        "replicate_halves": _half_text(chem, mir),
    }])
    write(res, "mir200a5p_call.tsv")
    write(famtab, "mir200_family_dehp.tsv")
    write(fam_rows, "mir200_family_dehp_tests.tsv")
    print("[mir200]\n" + res.T.to_string(header=False))
    print(famtab.to_string(index=False))
    print(fam_rows.to_string(index=False))


def offon(off=0.1, on=1.0):
    """R3-6: genuine off-to-on / on-to-off events. 'off' = TPM < `off` in BOTH
    replicates; 'on' = TPM >= `on` in BOTH replicates (a presence/absence call that
    does not depend on any ratio)."""
    tpm, _ = matrices()
    c = tpm[runs("mrna", CONTROL)]
    rows = []
    for chem in CHEMICALS:
        t = tpm[runs("mrna", chem)]
        up = (c < off).all(axis=1) & (t >= on).all(axis=1)
        dn = (c >= on).all(axis=1) & (t < off).all(axis=1)
        l2 = pd.Series(compute_log2fc_robust(t.mean(axis=1).to_dict(), c.mean(axis=1).to_dict()))
        rows.append({"chemical": chem, "off_to_on": int(up.sum()), "on_to_off": int(dn.sum()),
                     "off_to_on_robust_log2fc_min": round(float(l2[up[up].index].min()), 2)
                     if up.any() else np.nan,
                     "on_to_off_robust_log2fc_max": round(float(l2[dn[dn].index].max()), 2)
                     if dn.any() else np.nan,
                     "off_to_on_genes": ",".join(sorted(up[up].index)),
                     "on_to_off_genes": ",".join(sorted(dn[dn].index))})
    out = pd.DataFrame(rows)
    write(out, "off_on.tsv")
    print("[offon]\n" + out.drop(columns=["off_to_on_genes", "on_to_off_genes"]).to_string(index=False))


def targetscan():
    """R1-9: TargetScan Mouse 8.0 targets looked up through the official miRBase-ID ->
    family map (config/targetscan_mmu_mature.tsv, built by scripts/build_targetscan.py).
    The original sensitivity run looked miRNAs up by exact family name, which only
    matches single-member families (0/20 evaluable)."""
    mi, mr, _, base = primary_inputs()
    d = {}
    with open("config/targetscan_mmu_mature.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            d.setdefault(r["miRNA"], set()).add(r["Gene"])
    df = run_coherence(mi, mr, d, {}, base, **PRIMARY)
    s = summary(df, "TargetScan Mouse 8.0 via miRBase-ID -> family map")
    # the old exact-family-name lookup, reproduced for the record
    old = {}
    with open("config/targetscan_mmu.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            old.setdefault(r["miRNA"], set()).add(r["Gene"])
    dfo = run_coherence(mi, mr, old, {}, base, **PRIMARY)
    so = summary(dfo, "TargetScan, exact family-name lookup (as submitted)")
    responsive = int(((mi.abs() >= 1.0)).sum().sum())
    out = pd.DataFrame([so, s])
    out["responsive_tests"] = responsive
    write(out, "targetscan.tsv")
    write(df, "targetscan_coherence.tsv")
    print("[targetscan]\n" + out.to_string(index=False))


def _coherence_with(mirna, mrna, mb, bgfun):
    """run_coherence with a caller-supplied background (bgfun(targets) -> set) and the
    same productive / anti-coherent rule (direction-aware)."""
    from tmc.targets import target_set, harmonize_mirna_names
    from tmc.coherence import directional_shift, cri
    from tmc.pipeline import _harmonize_target_dict
    mbh = _harmonize_target_dict(mb)
    genes = set(mrna.index)
    rows = []
    for chem in mirna.columns:
        col = mrna[chem]
        for m, l in mirna[chem].items():
            if pd.isna(l) or abs(l) < 1.0:
                continue
            tg = target_set(harmonize_mirna_names(m), {}, mbh, mode="union") & genes
            if len(tg) < 5:
                continue
            bg = bgfun(tg)
            s = directional_shift([col[g] for g in tg], [col[g] for g in bg if g in col.index])
            rows.append({"chemical": chem, "miRNA": m, "mirna_log2fc": l, "n_targets": len(tg),
                         "median_diff": s["median_diff"], "auc": s["auc"], "p": s["p"],
                         "cri": cri(l, s["median_diff"], s["auc"])})
    df = pd.DataFrame(rows)
    df["fdr"] = multipletests(df["p"], method="fdr_bh")[1]
    sig = (df["fdr"] < 0.05) & ((df["auc"] - 0.5).abs() >= 0.1)
    df["class"] = np.where(sig & (df["cri"] > 0), "productive",
                           np.where(sig, "anti_coherent", "passenger"))
    return df


def _pshape(df, label):
    """Excess of small p-values and whether it follows the miRNA's direction."""
    t = df[df["class"] != "untested"]
    small = t[t["p"] < 0.05]
    pi0 = min(1.0, float((t["p"] > 0.5).sum()) / (0.5 * len(t)))
    return {**summary(df, label), "p_lt_0.05": len(small),
            "p_lt_0.05_expected": round(0.05 * len(t), 1),
            "p_lt_0.05_targets_up": int((small["auc"] < 0.5).sum()),
            "p_lt_0.05_coherent_direction": int((small["cri"] > 0).sum()),
            "storey_pi0": round(pi0, 3)}


def length2d():
    """Background matched jointly on expression and transcript length (the
    submitted version's 2D sensitivity analysis, re-run with the direction-aware
    rule), for the toxicogenomic data and the Dicer1-KO physiological control. Also
    reports the shape of the p-value distribution (R1-3)."""
    from sensitivity import gene_length, bg_2d
    glen = gene_length()
    mi, mr, mb, base = primary_inputs()
    ref = pd.read_csv("results/coherence.tsv", sep="\t")
    rows = [_pshape(ref, "toxicogenomic, expression-matched (primary)")]
    d2 = _coherence_with(mi, mr, mb, lambda tg: bg_2d(tg, base, glen))
    rows.append(_pshape(d2, "toxicogenomic, expression+length matched"))
    write(d2, "length2d_coherence.tsv")
    import physctrl
    pm, pcov = physctrl.expression()
    pc = pd.read_csv(OUT / "physctrl_coherence.tsv", sep="\t")
    pmir = pd.DataFrame({"dicer1_ko": -5.0}, index=pc["miRNA"].unique())
    p2 = _coherence_with(pmir, pm, mb, lambda tg: bg_2d(tg, pcov, glen))
    p2["set"] = p2["miRNA"].map(dict(zip(pc["miRNA"], pc["set"])))
    for name, g in p2.groupby("set"):
        rows.append({**summary(g, f"Dicer1 KO, {name}, expression+length matched")})
    for name, g in pc.groupby("set"):
        rows.append({**summary(g, f"Dicer1 KO, {name}, expression-matched")})
    write(p2, "physctrl_length2d_coherence.tsv")
    out = pd.DataFrame(rows)
    write(out, "length2d.tsv")
    print("[length2d]\n" + out.to_string(index=False))


def sensdir():
    """The submitted robustness variants (pseudocount 0.5/5, miRDB >= 90, expression
    filters, up/down split), recounted with the paper's direction-aware rule through
    run_coherence. (scripts/sensitivity.py counted a significant shift in EITHER
    direction as 'productive'; since the primary analysis had no anti-coherent call the
    old counts were upper bounds, but the rule should be the same everywhere.)"""
    tpm, cpm = matrices()
    exg, exm, g = expressed(tpm, cpm)
    base = baseline_covariate(g, exg)
    mb = load_mirdb("config/mirdb_mmu.tsv", score_min=80)
    mb90 = load_mirdb("config/mirdb_mmu.tsv", score_min=90)

    def robust(mat, layer, ps):
        ctrl = mat[runs(layer, CONTROL)].mean(axis=1).to_dict()
        return pd.DataFrame({c: compute_log2fc_robust(mat[runs(layer, c)].mean(axis=1).to_dict(),
                                                      ctrl, pseudo=ps) for c in CHEMICALS})
    rows = []
    for ps in (0.5, 1.0, 5.0):
        mi = robust(cpm, "mirna", ps); mi = mi[mi.index.isin(exm)]
        mr = robust(tpm, "mrna", ps); mr = mr[mr.index.isin(exg)]
        df = run_coherence(mi, mr, {}, mb, base, **PRIMARY)
        rows.append(summary(df, f"pseudocount {ps}"))
    mi, mr, _, _ = primary_inputs()
    rows.append(summary(run_coherence(mi, mr, {}, mb90, base, **PRIMARY), "miRDB score >= 90"))
    ref = pd.read_csv("results/coherence.tsv", sep="\t")
    rows.append(summary(ref[ref["mirna_log2fc"] > 0], "up-regulated miRNAs only"))
    rows.append(summary(ref[ref["mirna_log2fc"] < 0], "down-regulated miRNAs only"))
    for tth, cth in ((0.5, 5.0), (2.0, 20.0)):
        eg = set(g[g >= tth].index)
        m = cpm[runs("mirna", CONTROL)].mean(axis=1); em = set(m[m >= cth].index)
        b2 = baseline_covariate(g, eg)
        mi = load_matrix("results/mirna_log2fc_robust.tsv"); mi = mi[mi.index.isin(em)]
        mr = load_matrix("results/mrna_log2fc_robust.tsv"); mr = mr[mr.index.isin(eg)]
        rows.append(summary(run_coherence(mi, mr, {}, mb, b2, **PRIMARY),
                            f"expression filter TPM >= {tth}, CPM >= {cth}"))
    out = pd.DataFrame(rows)
    write(out, "sensitivity_directional.tsv")
    print("[sensdir]\n" + out.to_string(index=False))


def replicate_noise():
    """R1-2: the replicate-level numbers quoted in Results 3.1 and 3.6 that are not already
    in replicates.tsv / reproducibility.tsv. Needs replicates_paired_tests.tsv (stage
    'replicates') to exist.
      abundance  : Pearson r of log2(x + 1) between replicate 1 and 2, all 10 conditions
      two-fold   : expressed miRNAs with robust |log2FC| >= 1 between two single samples
      paired     : per-test agreement between halves, direction and coherence-signed AUC"""
    tpm, cpm = matrices()
    exg, exm, g = expressed(tpm, cpm)
    rows = []
    for name, layer, mat in (("mRNA TPM", "mrna", tpm), ("miRNA CPM", "mirna", cpm)):
        r = [np.corrcoef(np.log2(mat[runs(layer, c, 1)[0]] + 1),
                         np.log2(mat[runs(layer, c, 2)[0]] + 1))[0, 1] for c in [CONTROL] + CHEMICALS]
        rows.append({"measure": f"{name}: Pearson r of log2(x+1), replicate 1 vs 2 (10 conditions)",
                     "n": len(r), "median": round(float(np.median(r)), 3),
                     "min": round(float(min(r)), 3), "max": round(float(max(r)), 3)})
    # the submitted "replicate r≈0.99" did not state its metric; it is reproduced only on linear abundance
    for name, layer, mat in (("mRNA TPM", "mrna", tpm), ("miRNA CPM", "mirna", cpm)):
        r = [np.corrcoef(mat[runs(layer, c, 1)[0]], mat[runs(layer, c, 2)[0]])[0, 1]
             for c in [CONTROL] + CHEMICALS]
        rows.append({"measure": f"linear scale (as submitted): {name} Pearson r, replicate 1 vs 2 (10 conditions)",
                     "n": len(r), "median": round(float(np.median(r)), 4),
                     "min": round(float(min(r)), 4), "max": round(float(max(r)), 4)})
    feats = list(exm)

    def n2fold(num, den):
        l2 = pd.Series(compute_log2fc_robust(cpm[num].to_dict(), cpm[den].to_dict())).loc[feats]
        return int((l2.abs() >= 1).sum())
    ctrl = n2fold(runs("mirna", CONTROL, 2)[0], runs("mirna", CONTROL, 1)[0])
    within = [n2fold(runs("mirna", c, 2)[0], runs("mirna", c, 1)[0]) for c in CHEMICALS]
    versus = [n2fold(runs("mirna", c, k)[0], runs("mirna", CONTROL, k)[0])
              for k in (1, 2) for c in CHEMICALS]
    for label, v in (("vehicle control, replicate 2 vs 1", [ctrl]),
                     ("each chemical, replicate 2 vs 1", within),
                     ("chemical replicate k vs control replicate k", versus)):
        rows.append({"measure": f"expressed miRNAs (of {len(feats)}) differing >= 2-fold: {label}",
                     "n": len(v), "median": float(np.median(v)), "min": min(v), "max": max(v)})
    m = pd.read_csv(OUT / "replicates_paired_tests.tsv", sep="\t")
    s1, s2 = np.sign(m["mirna_log2fc_1"]), np.sign(m["mirna_log2fc_2"])
    coh1, coh2 = (m["auc_1"] - 0.5) * s1, (m["auc_2"] - 0.5) * s2   # > 0 = coherent
    rho = spearmanr(coh1, coh2)
    same = int((s1 == s2).sum())
    rows += [{"measure": "tests evaluated in both halves", "n": len(m)},
             {"measure": "  miRNA changed in the same direction in both halves", "n": same,
              "median": round(100 * same / len(m), 1)},
             {"measure": "  coherence-signed AUC deviation, half 1 vs half 2 (Spearman rho; p in max)",
              "n": len(m), "median": round(float(rho.statistic), 3), "max": float(rho.pvalue)},
             {"measure": "  coherent in both halves", "n": int(((coh1 > 0) & (coh2 > 0)).sum())},
             {"measure": "  anti-coherent in both halves", "n": int(((coh1 < 0) & (coh2 < 0)).sum())}]
    out = pd.DataFrame(rows)
    write(out, "replicate_noise.tsv")
    print("[replicate_noise]\n" + out.to_string(index=False))


def zero_replacement():
    """Results 3.2: how many miRNA x chemical tests the zero-replacement estimator makes
    'responsive' (and evaluable), against the robust estimator, on the same expressed
    features and the same coherence pipeline."""
    tpm, cpm = matrices()
    exg, exm, g = expressed(tpm, cpm)
    base = baseline_covariate(g, exg)
    mb = load_mirdb("config/mirdb_mmu.tsv", score_min=80)
    rows = []
    for label, mfile, gfile in (
            ("zero-replacement log2FC (miRNA and mRNA)", "results/mirna_log2fc.tsv",
             "results/mrna_log2fc.tsv"),
            ("robust pseudocount log2FC (primary)", "results/mirna_log2fc_robust.tsv",
             "results/mrna_log2fc_robust.tsv")):
        mi = load_matrix(mfile); mi = mi[mi.index.isin(exm)]
        mr = load_matrix(gfile); mr = mr[mr.index.isin(exg)]
        df = run_coherence(mi, mr, {}, mb, base, **PRIMARY)
        t = df[df["class"] != "untested"]
        rows.append({"estimator": label, "responsive_tests": int((mi.abs() >= 1.0).sum().sum()),
                     **summary(df, label), "median_targets": float(t["n_targets"].median())})
    out = pd.DataFrame(rows).drop(columns=["variant"])
    write(out, "zero_replacement.tsv")
    print("[zero_replacement]\n" + out.to_string(index=False))


def physctrl_fisher():
    """Results 3.4: specificity of the Dicer1-knockout control. Productive fraction of
    expressed vs not-expressed (no shared seed) miRNAs, Fisher's exact test (two-sided),
    with the expression-matched and the expression-and-length-matched background; and the
    miR-290-295 cluster 3p members. Needs physctrl_coherence.tsv (scripts/physctrl.py) and
    physctrl_length2d_coherence.tsv (stage 'length2d')."""
    from scipy.stats import fisher_exact
    rows = []
    for label, name in (("expression-matched", "physctrl_coherence.tsv"),
                        ("expression+length matched", "physctrl_length2d_coherence.tsv")):
        d = pd.read_csv(OUT / name, sep="\t")
        d = d[d["class"] != "untested"]
        k = {s: int((d[d["set"] == s]["class"] == "productive").sum()) for s in ("expressed", "not expressed")}
        n = {s: int((d["set"] == s).sum()) for s in ("expressed", "not expressed")}
        res = fisher_exact([[k["expressed"], n["expressed"] - k["expressed"]],
                            [k["not expressed"], n["not expressed"] - k["not expressed"]]])
        cl = d[d["miRNA"].str.match(r"mmu-miR-29[0-5][ab]?-3p$")]
        rows.append({"background": label,
                     "expressed_productive": k["expressed"], "expressed_tested": n["expressed"],
                     "expressed_pct": round(100 * k["expressed"] / n["expressed"], 1),
                     "not_expressed_productive": k["not expressed"],
                     "not_expressed_tested": n["not expressed"],
                     "not_expressed_pct": round(100 * k["not expressed"] / n["not expressed"], 1),
                     "odds_ratio": round(float(res.statistic), 2), "fisher_p": float(res.pvalue),
                     "cluster_290_295_3p_tested": len(cl),
                     "cluster_290_295_3p_productive": int((cl["class"] == "productive").sum()),
                     "cluster_290_295_3p_auc_min": round(float(cl["auc"].min()), 3),
                     "cluster_290_295_3p_auc_max": round(float(cl["auc"].max()), 3)})
    out = pd.DataFrame(rows)
    write(out, "physctrl_fisher.tsv")
    print("[physctrl_fisher]\n" + out.to_string(index=False))


def mirdb():
    """Methods 2.3: provenance of config/mirdb_mmu.tsv. Rebuilds it from the public miRDB v6.0
    file (ref/miRDB.txt.gz; see ref/SOURCE_miRDB.txt) with the GENCODE vM36 RefSeq bridge and
    checks that the rebuild is byte-identical to the committed table."""
    import gzip
    import hashlib
    from build_targets import load_refseq2symbol
    r2s = load_refseq2symbol()
    best, total, mapped = {}, 0, 0
    with gzip.open("ref/miRDB.txt.gz", "rt") as fh:
        for line in fh:
            mir, refseq, score = line.rstrip("\n").split("\t")
            if not mir.startswith("mmu-"):
                continue
            total += 1
            sym = r2s.get(refseq.split(".")[0])
            if sym is None:
                continue
            mapped += 1
            key, sc = (mir, sym), float(score)
            if key not in best or sc > best[key]:
                best[key] = sc
    rebuilt = "".join(f"{m}\t{g}\t{sc:.2f}\n" for (m, g), sc in sorted(best.items()))
    same = rebuilt == Path("config/mirdb_mmu.tsv").read_text()
    out = pd.DataFrame([{"mmu_rows": total, "mapped_rows": mapped,
                         "mapped_pct": round(100 * mapped / total, 1), "pairs": len(best),
                         "identical_to_config": same,
                         "sha256_rebuilt": hashlib.sha256(rebuilt.encode()).hexdigest()}])
    write(out, "mirdb_provenance.tsv")
    print("[mirdb]\n" + out.to_string(index=False))


def margins():
    """Methods 2.4: the call under three |AUC - 0.5| margins, on the primary results."""
    c = pd.read_csv("results/coherence.tsv", sep="\t")
    t = c[c["class"] != "untested"]
    rows = []
    for m in (0.08, 0.10, 0.12):
        sig = (t["fdr"] < 0.05) & ((t["auc"] - 0.5).abs() >= m)
        rows.append({"margin": m, "tested": len(t), "productive": int((sig & (t["cri"] > 0)).sum()),
                     "anti_coherent": int((sig & (t["cri"] <= 0)).sum())})
    out = pd.DataFrame(rows)
    write(out, "margins.tsv")
    print("[margins]\n" + out.to_string(index=False))


def power_coherent(deltas=(0.0, 0.1, 0.2, 0.3, 0.5), seeds=range(20)):
    """Results 3.5 / Fig. 5: spike-in in the COHERENT direction of each test, classified exactly as the
    real data. The submitted results/power.tsv (scripts/power_check.py) shifted every synthetic target set
    DOWN and counted any significant shift (FDR < 0.05, |AUC - 0.5| >= 0.1) whatever its direction; for the
    152 tests of down-regulated miRNAs the coherent movement is UP. Here the shift is -delta for up-regulated
    and +delta for down-regulated miRNAs; a test is productive when FDR < 0.05 (BH over all tests at that
    delta), |AUC - 0.5| >= 0.1 and CRI > 0, anti-coherent when significant with CRI <= 0. Synthetic target
    values are drawn with replacement from the test's own expression-matched background (primary seed 0 for
    the background, as in the real analysis); the draw is repeated for 20 random seeds (seed 0 = primary)."""
    import random
    from tmc.targets import target_set, harmonize_mirna_names
    from tmc.coherence import directional_shift, cri
    from tmc.pipeline import _harmonize_target_dict
    from tmc.nullmodel import matched_background
    mi, mr, mb, base = primary_inputs()
    mb = _harmonize_target_dict(mb)
    exg = set(mr.index)
    tests = []
    for chem in mi.columns:
        col = mr[chem]
        for m, l in mi[chem].items():
            if pd.isna(l) or abs(l) < 1.0:
                continue
            tg = target_set(harmonize_mirna_names(m), {}, mb, mode="union") & exg
            if len(tg) < 5:
                continue
            bg = matched_background(tg, base)
            tests.append((float(l), len(tg), sorted(col[g] for g in bg if g in col.index)))
    rows, per_seed = [], []
    for d in deltas:
        for s in seeds:
            rng = random.Random(s)
            ps, aucs, cris = [], [], []
            for l, nt, bgv in tests:
                shift = -d if l > 0 else d
                tv = [rng.choice(bgv) + shift for _ in range(nt)]
                r = directional_shift(tv, bgv)
                ps.append(r["p"]); aucs.append(r["auc"]); cris.append(cri(l, r["median_diff"], r["auc"]))
            fdr = multipletests(ps, method="fdr_bh")[1]
            sig = (fdr < 0.05) & (np.abs(np.array(aucs) - 0.5) >= 0.1)
            cr = np.array(cris)
            # coherence-oriented AUC: > 0.5 when targets move in the coherent direction
            coh_auc = [a if l > 0 else 1 - a for a, (l, _, _) in zip(aucs, tests)]
            per_seed.append({"delta": d, "seed": s, "tests": len(tests), "productive": int((sig & (cr > 0)).sum()),
                             "anti_coherent": int((sig & (cr <= 0)).sum()), "mean_coherence_auc": float(np.mean(coh_auc))})
    ps_df = pd.DataFrame(per_seed)
    for d, g in ps_df.groupby("delta"):
        p0 = g[g.seed == 0].iloc[0]
        pct = 100 * g.productive / g.tests
        rows.append({"delta": d, "tests": int(p0.tests), "productive_seed0": int(p0.productive),
                     "pct_productive_seed0": round(100 * p0.productive / p0.tests, 1),
                     "anti_coherent_seed0": int(p0.anti_coherent),
                     "pct_productive_min_20_seeds": round(float(pct.min()), 1),
                     "pct_productive_max_20_seeds": round(float(pct.max()), 1),
                     "anti_coherent_max_20_seeds": int(g.anti_coherent.max()),
                     "mean_coherence_auc_seed0": round(float(p0.mean_coherence_auc), 3)})
    out = pd.DataFrame(rows)
    write(out, "power_coherent.tsv")
    print("[power_coherent]\n" + out.to_string(index=False))


def physctrl_family():
    """Results 3.4: the Dicer1 Fisher test with non-independent miRNAs collapsed. miRNAs sharing a seed
    (TargetScan Mouse 8.0 Seed+m8) or an identical predicted target set give correlated or identical tests,
    so the miRNA-level test counts some evidence more than once. Two collapses: (a) seed family, productive
    when at least half of its evaluable members are productive (miRNAs without a family entry stand alone);
    (b) identical miRDB target set (identical tests counted once)."""
    import physctrl
    from scipy.stats import fisher_exact
    seed = physctrl.seeds()
    rows = []
    for label, name in (("expression-matched", "physctrl_coherence.tsv"),
                        ("expression+length matched", "physctrl_length2d_coherence.tsv")):
        d = pd.read_csv(OUT / name, sep="\t")
        d = d[d["class"] != "untested"].copy()
        if "set" not in d:
            d["set"] = d["miRNA"].map(dict(zip(pd.read_csv(OUT / "physctrl_coherence.tsv", sep="\t")["miRNA"],
                                               pd.read_csv(OUT / "physctrl_coherence.tsv", sep="\t")["set"])))
        d["prod"] = d["class"] == "productive"
        d["family"] = [seed.get(m.lower(), m) for m in d["miRNA"]]
        d["tset"] = d["n_targets"].astype(str) + "|" + d["auc"].round(12).astype(str)
        for unit, key in (("seed family", "family"), ("identical target set", "tset")):
            g = d.groupby(["set", key])["prod"].mean().reset_index()
            g["prod"] = g["prod"] >= 0.5
            k = {s: int(g[g.set == s]["prod"].sum()) for s in ("expressed", "not expressed")}
            n = {s: int((g.set == s).sum()) for s in ("expressed", "not expressed")}
            res = fisher_exact([[k["expressed"], n["expressed"] - k["expressed"]],
                                [k["not expressed"], n["not expressed"] - k["not expressed"]]])
            rows.append({"background": label, "unit": unit, "expressed_productive": k["expressed"],
                         "expressed_units": n["expressed"], "not_expressed_productive": k["not expressed"],
                         "not_expressed_units": n["not expressed"], "odds_ratio": round(float(res.statistic), 2),
                         "fisher_p": float(res.pvalue)})
    out = pd.DataFrame(rows)
    write(out, "physctrl_family.tsv")
    print("[physctrl_family]\n" + out.to_string(index=False))


def multiplicity():
    """Methods 2.4: the primary call under Benjamini-Yekutieli, which controls the FDR under any dependence
    (Benjamini-Hochberg needs positive regression dependence, which shared controls make plausible but do
    not prove)."""
    c = pd.read_csv("results/coherence.tsv", sep="\t")
    t = c[c["class"] != "untested"].copy()
    rows = []
    for proc, meth in (("Benjamini-Hochberg (primary)", "fdr_bh"), ("Benjamini-Yekutieli", "fdr_by")):
        q = multipletests(t["p"], method=meth)[1]
        sig = (q < 0.05) & ((t["auc"] - 0.5).abs() >= 0.1)
        rows.append({"procedure": proc, "tested": len(t), "productive": int((sig & (t["cri"] > 0)).sum()),
                     "anti_coherent": int((sig & (t["cri"] <= 0)).sum()), "min_adjusted_p": round(float(q.min()), 4)})
    out = pd.DataFrame(rows)
    write(out, "multiplicity.tsv")
    print("[multiplicity]\n" + out.to_string(index=False))


def normalisation():
    """Results 3.1: per-chemical median log2FC over expressed features under library-size normalisation
    (miRNA CPM; mRNA TPM) and under median-of-ratios size factors (each sample divided by the median, over
    features detected in every sample, of its ratio to the feature's geometric mean across samples; applied
    to miRNA counts and to mRNA TPM). log2FC = log2((mean chemical + 1) / (mean control + 1))."""
    tpm, cpm = matrices()
    exg, exm, _ = expressed(tpm, cpm)
    cnt = load_matrix("results/mirna_counts.tsv")
    rows = []
    for layer, lib, raw, feats in (("miRNA", cpm, cnt, exm), ("mRNA", tpm, tpm, exg)):
        pos = raw[(raw > 0).all(axis=1)]
        geo = np.exp(np.log(pos).mean(axis=1))
        sf = pos.div(geo, axis=0).median(axis=0)
        mor = raw.div(sf, axis=1)
        for norm, mat in (("library size", lib), ("median of ratios", mor)):
            meds = []
            ctrl = mat[runs(layer.lower(), CONTROL)].mean(axis=1)
            for c in CHEMICALS:
                l2 = np.log2((mat[runs(layer.lower(), c)].mean(axis=1) + 1) / (ctrl + 1))
                meds.append(float(l2[l2.index.isin(feats)].median()))
            rows.append({"layer": layer, "normalisation": norm, "n_features": len(feats),
                         "median_log2fc_min": round(min(meds), 4), "median_log2fc_max": round(max(meds), 4)})
    out = pd.DataFrame(rows)
    write(out, "normalisation.tsv")
    print("[normalisation]\n" + out.to_string(index=False))


def gene_biotype():
    """GENCODE vM36 gene symbol -> protein_coding / lncRNA / other (a gene is protein-coding if any of its
    transcripts is), from the transcript FASTA headers (field 8)."""
    import gzip
    import collections
    bt = collections.defaultdict(set)
    with gzip.open("ref/transcripts.fa.gz", "rt") as fh:
        for line in fh:
            if line[0] == ">":
                f = line[1:].split("|")
                bt[f[5]].add(f[7])
    return {g: "protein_coding" if "protein_coding" in s else ("lncRNA" if "lncRNA" in s else "other")
            for g, s in bt.items()}


def biotype():
    """Every predicted target is protein-coding, but the expression-matched background was drawn from all
    expressed genes (lncRNAs included). Re-run the primary analysis and the Dicer1-KO control with the
    background restricted to (a) protein-coding genes and (b) genes that miRDB predicts as a target of any
    miRNA (score >= 80), and report the biotype composition of targets and backgrounds."""
    import collections
    from tmc.targets import target_set, harmonize_mirna_names
    from tmc.pipeline import _harmonize_target_dict
    from tmc.nullmodel import matched_background
    gbt = gene_biotype()
    mi, mr, mb, base = primary_inputs()
    mbh = _harmonize_target_dict(mb)
    universe = set().union(*mbh.values())
    exg = set(mr.index)
    comp = []
    tc, bc = collections.Counter(), collections.Counter()
    for chem in mi.columns:
        for m, l in mi[chem].items():
            if pd.isna(l) or abs(l) < 1.0:
                continue
            tg = target_set(harmonize_mirna_names(m), {}, mbh, mode="union") & exg
            if len(tg) < 5:
                continue
            tc.update(gbt.get(g, "unmapped") for g in tg)
            bc.update(gbt.get(g, "unmapped") for g in matched_background(tg, base))
    uc = collections.Counter(gbt.get(g, "unmapped") for g in exg)
    for label, cnt in (("expressed genes (universe)", uc), ("predicted targets, summed over tests", tc),
                       ("primary background draws, summed over tests", bc)):
        tot = sum(cnt.values())
        comp.append({"set": label, "n": tot, **{f"pct_{k}": round(100 * cnt.get(k, 0) / tot, 1)
                                               for k in ("protein_coding", "lncRNA", "other", "unmapped")}})
    from scipy.stats import binomtest, wilcoxon

    def direction(d, label):
        """summary + does the shift follow the miRNA? (coherence-oriented AUC, > 0.5 = coherent direction)"""
        t = d[d["class"] != "untested"]
        coh = np.where(t["mirna_log2fc"] > 0, t["auc"], 1 - t["auc"])
        p, a = int((t["class"] == "productive").sum()), int((t["class"] == "anti_coherent").sum())
        return {**summary(d, label), "coherence_auc_mean": round(float(coh.mean()), 3),
                "coherence_auc_wilcoxon_p": round(float(wilcoxon(coh - 0.5).pvalue), 2),
                "coherent_share_binomial_p": round(float(binomtest(p, p + a).pvalue), 2) if p + a else None}
    rows = [direction(pd.read_csv("results/coherence.tsv", sep="\t"), "toxicogenomic, all expressed genes (primary)"),
            direction(pd.read_csv(OUT / "length2d_coherence.tsv", sep="\t"), "toxicogenomic, expression+length matched")]
    restrict = {"protein-coding background": lambda g: gbt.get(g) == "protein_coding",
                "background = genes predicted for any miRNA": lambda g: g in universe}
    import physctrl
    pm, pcov = physctrl.expression()
    pc = pd.read_csv(OUT / "physctrl_coherence.tsv", sep="\t")
    pmir = pd.DataFrame({"dicer1_ko": -5.0}, index=pc["miRNA"].unique())
    from scipy.stats import fisher_exact
    for name, keep in restrict.items():
        cov = {g: v for g, v in base.items() if keep(g)}
        d = _coherence_with(mi, mr, mb, lambda tg: matched_background(tg, cov))
        rows.append(direction(d, f"toxicogenomic, {name}"))
        write(d, f"biotype_{name.split()[0].replace('-', '')}_coherence.tsv")
        pcv = {g: v for g, v in pcov.items() if keep(g)}
        p = _coherence_with(pmir, pm, mb, lambda tg: matched_background(tg, pcv))
        p["set"] = p["miRNA"].map(dict(zip(pc["miRNA"], pc["set"])))
        e, n = p[p.set == "expressed"], p[p.set == "not expressed"]
        a, b = int((e["class"] == "productive").sum()), int((n["class"] == "productive").sum())
        odds, pv = fisher_exact([[a, len(e) - a], [b, len(n) - b]])
        rows.append({"variant": f"Dicer1 KO, {name}", "tested": len(e) + len(n),
                     "productive": a + b, "anti_coherent": int((p["class"] == "anti_coherent").sum()),
                     "expressed_productive": a, "expressed_tested": len(e),
                     "not_expressed_productive": b, "not_expressed_tested": len(n),
                     "odds_ratio": round(float(odds), 2), "fisher_p": float(pv)})
    import random
    import posctrl_cri as P
    pmr = load_matrix("results/posctrl_mrna_log2fc.tsv")
    ptpm = load_matrix("results/posctrl_gene_tpm.tsv")
    pacc = pd.read_csv("config/posctrl_accessions.tsv", sep="\t")
    cm = ptpm[list(pacc.run_accession[pacc.condition == "control"])].mean(axis=1)
    pex = set(cm[cm >= 1.0].index)
    pmr = pmr[pmr.index.isin(pex)]
    pbase = {g: float(np.log10(v + 1.0)) for g, v in cm.items() if g in pex}
    rng = random.Random(0)
    panel = rng.sample([m for m in mb if m not in P.TRANSFECTED and len(mb[m]) >= 50], P.N_RANDOM)
    mmir = pd.DataFrame({c: {m: 5.0 for m in P.TRANSFECTED + panel} for c in P.CELLS})
    for name, keep in [("all expressed genes (as published)", lambda g: True)] + list(restrict.items()):
        cov = {g: v for g, v in pbase.items() if keep(g)}
        d = _coherence_with(mmir, pmr, mb, lambda tg: matched_background(tg, cov))
        tr, rd = d[d["miRNA"].isin(P.TRANSFECTED)], d[~d["miRNA"].isin(P.TRANSFECTED)]
        rows.append({"variant": f"mimic transfection, {name}", "tested": len(d),
                     "transfected_productive": int((tr["class"] == "productive").sum()), "transfected_tested": len(tr),
                     "transfected_auc_median": round(float(tr["auc"].median()), 3),
                     "random_productive": int((rd["class"] == "productive").sum()), "random_tested": len(rd)})
    write(pd.DataFrame(comp), "biotype_composition.tsv")
    out = pd.DataFrame(rows)
    write(out, "biotype.tsv")
    print("[biotype]\n" + pd.DataFrame(comp).to_string(index=False) + "\n" + out.to_string(index=False))


STAGES = {"baseline": baseline, "window": window, "weight": weight,
          "replicates": replicates, "reproducibility": reproducibility, "mir200": mir200, "offon": offon,
          "targetscan": targetscan, "length2d": length2d, "sensdir": sensdir,
          "replicate_noise": replicate_noise, "zero_replacement": zero_replacement,
          "physctrl_fisher": physctrl_fisher, "mirdb": mirdb, "margins": margins,
          "power_coherent": power_coherent, "biotype": biotype,
          "physctrl_family": physctrl_family, "multiplicity": multiplicity,
          "normalisation": normalisation}

if __name__ == "__main__":
    arg = sys.argv[1] if len(sys.argv) > 1 else "all"
    order = list(STAGES) if arg == "all" else [arg]
    for name in order:
        STAGES[name]()
