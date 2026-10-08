"""Robustness / sensitivity analyses for the coherence null (QA Stage 3, R2).
Addresses reviewer requests: (1) jointly matching the non-target background on
expression AND transcript length (proxy for 3'UTR length); (2) robustness to the
pseudocount; (3) robustness to the miRDB score threshold; (4) up- vs down-regulated
miRNA split. Prints a compact report; the headline number is the productive count
(expected 0 throughout if the null is robust)."""
import bisect
import statistics
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from tmc.config import load_accessions, CONTROL, CHEMICALS
from tmc.io_expr import load_matrix, chemical_log2fc, compute_log2fc_robust
from tmc.targets import load_mirdb, target_set, harmonize_mirna_names
from tmc.coherence import directional_shift, cri, classify_response

MIN_CTRL_TPM, MIN_CTRL_CPM = 1.0, 10.0


def gene_length():
    """Max transcript length per gene from the GENCODE fasta headers (field 7)."""
    import gzip
    glen = {}
    with gzip.open("ref/transcripts.fa.gz", "rt") as fh:
        for line in fh:
            if line[0] == ">":
                f = line[1:].split("|")
                if len(f) >= 7:
                    sym, L = f[5], int(f[6])
                    glen[sym] = max(glen.get(sym, 0), L)
    return glen


def ctrl_runs(layer):
    acc = [r for r in load_accessions("config/accessions.tsv") if r["layer"] == layer]
    return [r["run_accession"] for r in acc if r["chemical"] == CONTROL]


def expressed():
    tpm = load_matrix("results/gene_tpm.tsv")
    g = tpm[ctrl_runs("mrna")].mean(axis=1)
    exg = set(g[g >= MIN_CTRL_TPM].index)
    cnt = load_matrix("results/mirna_counts.tsv")
    cpm = cnt / cnt.sum() * 1e6
    m = cpm[ctrl_runs("mirna")].mean(axis=1)
    exm = set(m[m >= MIN_CTRL_CPM].index)
    base = {gn: float(np.log10(v + 1.0)) for gn, v in g.items() if gn in exg}
    return exg, exm, base


def bg_2d(targets, expr, length, n=20, win=80, seed=0):
    """Background matched on expression (nearest window) then length (nearest within)."""
    import random
    rng = random.Random(seed)
    items = sorted((v, g) for g, v in expr.items() if g not in targets)
    vals = [v for v, _ in items]
    bg = set()
    for t in sorted(targets):
        if t not in expr or t not in length:
            continue
        i = bisect.bisect_left(vals, expr[t])
        lo, hi = max(0, i - win), min(len(items), i + win)
        pool = [g for _, g in items[lo:hi] if g not in bg and g in length]
        pool.sort(key=lambda g: abs(length[g] - length[t]))  # nearest in length
        bg.update(pool[:n])
    return bg


def run(mirna, mrna, mb, base, length=None, score_min=80, pseudo_note="",
        n_per_target=20):
    exg = set(mrna.index)
    rows = []
    for chem in mirna.columns:
        col = mrna[chem]
        for m, l in mirna[chem].items():
            if pd.isna(l) or abs(l) < 1.0:
                continue
            hk = harmonize_mirna_names(m)
            tg = target_set(hk, {}, mb, mode="union") & exg
            if len(tg) < 5:
                continue
            tv = [col[g] for g in tg]
            if length is not None:
                bg = bg_2d(tg, base, length)
            else:
                from tmc.nullmodel import matched_background
                bg = matched_background(tg, base, n_per_target=n_per_target)
            bv = [col[g] for g in bg if g in col.index]
            s = directional_shift(tv, bv)
            rows.append({"chemical": chem, "miRNA": m, "l": l, "auc": s["auc"],
                         "p": s["p"], "n": len(tg)})
    df = pd.DataFrame(rows)
    # BH across tests, then classify
    df = df.sort_values("p").reset_index(drop=True)
    mtot = len(df)
    df["fdr"] = (df["p"] * mtot / (df.index + 1)).clip(upper=1.0)[::-1].cummin()[::-1]
    df["prod"] = (df["fdr"] < 0.05) & ((df["auc"] - 0.5).abs() >= 0.1)
    print(f"  {pseudo_note}: tested={len(df)}  productive={int(df['prod'].sum())}  "
          f"AUC mean={df['auc'].mean():.3f}")
    return df


def main():
    from tmc.pipeline import _harmonize_target_dict
    exg, exm, base = expressed()
    mb = _harmonize_target_dict(load_mirdb("config/mirdb_mmu.tsv", score_min=80))
    mb90 = _harmonize_target_dict(load_mirdb("config/mirdb_mmu.tsv", score_min=90))
    mrna = load_matrix("results/mrna_log2fc_robust.tsv"); mrna = mrna[mrna.index.isin(exg)]
    mirna = load_matrix("results/mirna_log2fc_robust.tsv"); mirna = mirna[mirna.index.isin(exm)]
    glen = gene_length()

    print("=== Sensitivity analyses (productive count expected near 0 = FDR noise) ===")
    print("[baseline] expression-matched background:")
    base_df = run(mirna, mrna, mb, base, pseudo_note="expr-matched")

    # 1) length imbalance + 2D matched background
    # target vs expression-matched-bg transcript length
    from tmc.nullmodel import matched_background
    tlens, blens = [], []
    for chem in mirna.columns:
        for m, l in mirna[chem].items():
            if pd.isna(l) or abs(l) < 1.0:
                continue
            tg = target_set(harmonize_mirna_names(m), {}, mb, mode="union") & exg
            if len(tg) < 5:
                continue
            bg = matched_background(tg, base)
            tlens += [glen[g] for g in tg if g in glen]
            blens += [glen[g] for g in bg if g in glen]
    u, p = mannwhitneyu(tlens, blens, alternative="two-sided")
    print(f"\n[length covariate] target median tx length={statistics.median(tlens):.0f} nt vs "
          f"expr-matched bg={statistics.median(blens):.0f} nt (MWU p={p:.2e})")
    print("[1] expression+length 2D-matched background:")
    run(mirna, mrna, mb, base, length=glen, pseudo_note="2D expr+length")

    # 2) pseudocount robustness
    print("\n[2] pseudocount robustness (rebuild robust log2FC):")
    for ps in (0.5, 5.0):
        mr = rebuild_robust(ps, exg, exm)
        run(mr["mirna"], mr["mrna"], mb, base, pseudo_note=f"pseudo={ps}")

    # 3) miRDB threshold
    print("\n[3] miRDB score >= 90:")
    run(mirna, mrna, mb90, base, pseudo_note="miRDB>=90")

    # 4) up vs down split (baseline)
    up = base_df[base_df.l > 0]; dn = base_df[base_df.l < 0]
    print(f"\n[4] up-miRNAs: n={len(up)} productive={int(up['prod'].sum())} AUC mean={up['auc'].mean():.3f}"
          f" | down-miRNAs: n={len(dn)} productive={int(dn['prod'].sum())} AUC mean={dn['auc'].mean():.3f}")

    # 5) alternative target source: TargetScan Mouse
    from tmc.targets import load_targetscan
    from tmc.pipeline import run_coherence
    ts = load_targetscan("config/targetscan_mmu.tsv")
    df = run_coherence(mirna, mrna, ts, {}, base, responsive_thresh=1.0,
                       target_mode="union", min_targets=5, n_per_target=20, seed=0)
    t = df[df["class"] != "untested"]
    print(f"\n[5] TargetScan Mouse targets: tested={len(t)} "
          f"productive={int((t['class']=='productive').sum())} AUC mean={t.auc.mean():.3f}")

    # 6) expression-filter robustness
    print("\n[6] expression-filter robustness:")
    tpm = load_matrix("results/gene_tpm.tsv"); g = tpm[ctrl_runs("mrna")].mean(axis=1)
    cnt = load_matrix("results/mirna_counts.tsv"); cpm = cnt / cnt.sum() * 1e6
    mc = cpm[ctrl_runs("mirna")].mean(axis=1)
    for tpmth, cpmth in [(0.5, 5.0), (1.0, 10.0), (2.0, 20.0)]:
        exg2 = set(g[g >= tpmth].index); exm2 = set(mc[mc >= cpmth].index)
        b2 = {gn: float(np.log10(v + 1.0)) for gn, v in g.items() if gn in exg2}
        mi = load_matrix("results/mirna_log2fc_robust.tsv"); mi = mi[mi.index.isin(exm2)]
        mr = load_matrix("results/mrna_log2fc_robust.tsv"); mr = mr[mr.index.isin(exg2)]
        df = run_coherence(mi, mr, {}, mb, b2, responsive_thresh=1.0,
                           target_mode="union", min_targets=5, n_per_target=20, seed=0)
        t = df[df["class"] != "untested"]
        print(f"  TPM>={tpmth} CPM>={cpmth}: tested={len(t)} "
              f"productive={int((t['class']=='productive').sum())} AUC mean={t.auc.mean():.3f}")


def rebuild_robust(pseudo, exg, exm):
    acc_m = [r for r in load_accessions("config/accessions.tsv") if r["layer"] == "mirna"]
    acc_g = [r for r in load_accessions("config/accessions.tsv") if r["layer"] == "mrna"]
    cnt = load_matrix("results/mirna_counts.tsv")
    tpm = load_matrix("results/gene_tpm.tsv")
    est = lambda t, c: compute_log2fc_robust(t, c, pseudo=pseudo)
    mi = pd.DataFrame({c: chemical_log2fc(cnt, acc_m, c, estimator=est, cpm=True) for c in CHEMICALS})
    mr = pd.DataFrame({c: chemical_log2fc(tpm, acc_g, c, estimator=est, cpm=False) for c in CHEMICALS})
    return {"mirna": mi[mi.index.isin(exm)], "mrna": mr[mr.index.isin(exg)]}


if __name__ == "__main__":
    main()
