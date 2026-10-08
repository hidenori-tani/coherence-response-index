"""Power / spike-in analysis (QA Stage 3, R2): show the coherence test is not
merely under-powered. For every real evaluated (miRNA, chemical) test we keep the
real target-set size and the real expression-matched background, but replace the
target log2FCs with values drawn from the background and shifted DOWN by a known
coherent amount delta (simulated repression). We then re-classify all tests under
the same BH-FDR rule and report the fraction called productive at each delta. A
test that detects modest shifts but returns 0 on the real data is sensitive, not
blind.
"""
import random
import numpy as np
import pandas as pd

from tmc.config import load_accessions, CONTROL
from tmc.io_expr import load_matrix
from tmc.targets import load_mirdb, target_set, harmonize_mirna_names
from tmc.coherence import directional_shift
from tmc.nullmodel import matched_background
from tmc.pipeline import _harmonize_target_dict

MIN_CTRL_TPM, MIN_CTRL_CPM = 1.0, 10.0
DELTAS = [0.0, 0.1, 0.2, 0.3, 0.5]
RNG = random.Random(0)


def setup():
    acc = [r for r in load_accessions("config/accessions.tsv")]
    cr = lambda L: [r["run_accession"] for r in acc if r["layer"] == L and r["chemical"] == CONTROL]
    tpm = load_matrix("results/gene_tpm.tsv")
    g = tpm[cr("mrna")].mean(axis=1)
    exg = set(g[g >= MIN_CTRL_TPM].index)
    base = {gn: float(np.log10(v + 1.0)) for gn, v in g.items() if gn in exg}
    cnt = load_matrix("results/mirna_counts.tsv")
    cpm = cnt / cnt.sum() * 1e6
    m = cpm[cr("mirna")].mean(axis=1)
    exm = set(m[m >= MIN_CTRL_CPM].index)
    mb = _harmonize_target_dict(load_mirdb("config/mirdb_mmu.tsv", score_min=80))
    mirna = load_matrix("results/mirna_log2fc_robust.tsv"); mirna = mirna[mirna.index.isin(exm)]
    mrna = load_matrix("results/mrna_log2fc_robust.tsv"); mrna = mrna[mrna.index.isin(exg)]
    return mirna, mrna, base, mb, exg


def bh_productive(pvals, aucs):
    df = pd.DataFrame({"p": pvals, "auc": aucs}).sort_values("p").reset_index(drop=True)
    m = len(df)
    df["fdr"] = (df["p"] * m / (df.index + 1)).clip(upper=1.0)[::-1].cummin()[::-1]
    return int(((df["fdr"] < 0.05) & ((df["auc"] - 0.5).abs() >= 0.1)).sum()), m


def main():
    mirna, mrna, base, mb, exg = setup()
    # collect the real tests (target set + matched background pool of log2FCs)
    tests = []
    for chem in mirna.columns:
        col = mrna[chem]
        for mi, l in mirna[chem].items():
            if pd.isna(l) or abs(l) < 1.0:
                continue
            tg = target_set(harmonize_mirna_names(mi), {}, mb, mode="union") & exg
            if len(tg) < 5:
                continue
            bg = matched_background(tg, base)
            # sort so sampling is deterministic given the RNG seed (set order varies)
            bgvals = sorted(col[g] for g in bg if g in col.index)
            tests.append((len(tg), bgvals))
    print(f"power/spike-in over {len(tests)} real tests (real sizes + backgrounds)")
    print(f"{'delta(log2)':>11} {'productive':>10} {'%':>6} {'mean AUC':>9}")
    out = []
    for d in DELTAS:
        ps, aucs = [], []
        for nt, bgvals in tests:
            if len(bgvals) < 5:
                continue
            # synthetic 'targets' = sample from background, shifted DOWN by delta
            tv = [RNG.choice(bgvals) - d for _ in range(nt)]
            s = directional_shift(tv, bgvals)
            ps.append(s["p"]); aucs.append(s["auc"])
        nprod, m = bh_productive(ps, aucs)
        out.append({"delta": d, "n_tests": m, "n_productive": nprod,
                    "pct_productive": round(100 * nprod / m, 1), "mean_auc": round(float(np.mean(aucs)), 3)})
        print(f"{d:>11.1f} {nprod:>10} {100*nprod/m:>5.0f}% {np.mean(aucs):>9.3f}")
    pd.DataFrame(out).to_csv("results/power.tsv", sep="\t", index=False)
    print("wrote results/power.tsv")


if __name__ == "__main__":
    main()
