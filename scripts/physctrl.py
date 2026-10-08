"""Physiological positive control (CBI revision; R1-1, R3-3, R3-7).

Dicer1-knockout vs wild-type mouse embryonic stem cells (GEO GSE139348; five KO and
four WT clones, poly(A) RNA-seq, gene-level counts deposited by the authors). Loss of
Dicer1 removes the cell's own miRNAs at their endogenous abundance, so the predicted
targets of miRNAs that are genuinely active in mESCs should be de-repressed. This is a
positive control at physiological miRNA:target stoichiometry, in the same cell type as
the toxicogenomic data.

The CRI procedure is applied unchanged (miRDB >= 80, expression-matched background,
Mann-Whitney AUC, BH over all tests, the same productive rule). Every miRNA is entered
as down-regulated in the KO (log2FC = -5; only the sign enters the call). Two sets:
  expressed      : the 441 miRNAs expressed in our own vehicle-control mESCs (>= 10 CPM),
                   ranked by that abundance
  not expressed  : miRNAs < 1 CPM in our control mESCs whose seed (nt 2-8, TargetScan
                   miR_Family_Info) is not shared with any miRNA >= 1 CPM -- their loss
                   cannot matter, so their targets should not move (specificity control)
Gene-level expression = counts-per-million; log2FC = log2((KO mean + 1)/(WT mean + 1));
expressed genes = WT mean >= 1 CPM; matching covariate = log10(WT mean CPM + 1).

  python scripts/physctrl.py
"""
import csv
import gzip

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from tmc.config import load_accessions, CONTROL
from tmc.io_expr import load_matrix
from tmc.targets import load_mirdb
from tmc.pipeline import run_coherence
from revision_r1 import OUT, write

COUNTS = "data/physctrl_GSE139348/GSE139348_mESCs_dicer_RNAseq_counts_DEseq.csv.gz"
FAMILY = "data/targetscan_mmu80/miR_Family_Info.txt"


def expression():
    with gzip.open(COUNTS, "rt") as fh:
        df = pd.read_csv(fh, index_col=0)
    ko = [c for c in df.columns if "_KO" in c]
    wt = [c for c in df.columns if "_WT" in c]
    assert len(ko) == 5 and len(wt) == 4, (ko, wt)
    cnt = df[ko + wt].astype(float)
    cnt = cnt[~cnt.index.duplicated()]
    cpm = cnt / cnt.sum() * 1e6
    wt_m, ko_m = cpm[wt].mean(axis=1), cpm[ko].mean(axis=1)
    keep = wt_m[wt_m >= 1.0].index
    l2 = np.log2((ko_m + 1.0) / (wt_m + 1.0)).loc[keep]
    covar = {g: float(np.log10(v + 1.0)) for g, v in wt_m.loc[keep].items()}
    return pd.DataFrame({"dicer1_ko": l2}), covar


def mesc_mirna_cpm():
    acc = [r for r in load_accessions("config/accessions.tsv") if r["layer"] == "mirna"]
    ctrl = [r["run_accession"] for r in acc if r["chemical"] == CONTROL]
    cnt = load_matrix("results/mirna_counts.tsv")
    return (cnt / cnt.sum() * 1e6)[ctrl].mean(axis=1)


def seeds():
    s = {}
    with open(FAMILY) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["Species ID"] == "10090":
                s[r["MiRBase ID"].lower()] = r["Seed+m8"]
    return s


def main():
    mrna, covar = expression()
    cpm = mesc_mirna_cpm()
    mb = load_mirdb("config/mirdb_mmu.tsv", score_min=80)
    mb_names = {k.lower() for k in mb}
    sd = seeds()
    expressed = [m for m in cpm.index if cpm[m] >= 10.0]
    seeds_present = {sd[m.lower()] for m in cpm.index if cpm[m] >= 1.0 and m.lower() in sd}
    absent = [m for m in cpm.index if cpm[m] < 1.0 and m.lower() in mb_names
              and m.lower() in sd and sd[m.lower()] not in seeds_present]
    mirna = pd.DataFrame({"dicer1_ko": -5.0}, index=expressed + absent)
    df = run_coherence(mirna, mrna, {}, mb, covar, responsive_thresh=1.0,
                       target_mode="union", min_targets=5, n_per_target=20, seed=0)
    df["set"] = np.where(df["miRNA"].isin(expressed), "expressed", "not expressed")
    df["mesc_ctrl_cpm"] = df["miRNA"].map(cpm).round(2)
    t = df[df["class"] != "untested"].copy()
    rows = []
    for name, g in t.groupby("set"):
        rows.append({"set": name, "tested": len(g),
                     "productive": int((g["class"] == "productive").sum()),
                     "anti_coherent": int((g["class"] == "anti_coherent").sum()),
                     "auc_median": round(float(g["auc"].median()), 3)})
    e = t[t["set"] == "expressed"].sort_values("mesc_ctrl_cpm", ascending=False)
    e["abundance_rank"] = range(1, len(e) + 1)
    for lo, hi in ((1, 10), (11, 50), (51, 100), (101, len(e))):
        g = e[(e["abundance_rank"] >= lo) & (e["abundance_rank"] <= hi)]
        rows.append({"set": f"expressed, abundance rank {lo}-{hi}", "tested": len(g),
                     "productive": int((g["class"] == "productive").sum()),
                     "anti_coherent": int((g["class"] == "anti_coherent").sum()),
                     "auc_median": round(float(g["auc"].median()), 3)})
    rho = spearmanr(np.log10(e["mesc_ctrl_cpm"] + 1), 0.5 - e["auc"])
    rows.append({"set": "expressed: Spearman(log abundance, de-repression 0.5-AUC)",
                 "tested": len(e), "rho": round(float(rho.statistic), 3),
                 "p": float(rho.pvalue)})
    out = pd.DataFrame(rows)
    write(df, "physctrl_coherence.tsv")
    write(out, "physctrl_summary.tsv")
    print(f"genes: {len(mrna)} expressed (WT >= 1 CPM); miRNAs: {len(expressed)} expressed, "
          f"{len(absent)} not expressed (seed not shared)")
    print(out.to_string(index=False))
    print(e.head(12)[["miRNA", "mesc_ctrl_cpm", "n_targets", "auc", "fdr", "class"]]
          .to_string(index=False))


if __name__ == "__main__":
    main()
