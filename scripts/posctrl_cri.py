"""Positive-control CRI demonstration.

The three transfected miRNAs are overexpressed by design, so their predicted
targets should be repressed (CRI productive). A panel of non-transfected miRNAs,
treated identically, should not be (passenger). Same run_coherence / miRDB / null
as the toxicogenomics analysis. Run after scripts/posctrl_quant.py.
"""
import random
import numpy as np
import pandas as pd

from tmc.io_expr import load_matrix
from tmc.targets import load_mirdb
from tmc.pipeline import run_coherence

CELLS = ["bv2", "astrocyte", "neuro2a"]
TRANSFECTED = ["mmu-miR-199a-3p", "mmu-miR-99b-3p", "mmu-miR-140-5p"]
MIN_CTRL_TPM = 1.0
N_RANDOM = 60
OUT = "results/posctrl_coherence.tsv"


def main():
    mrna = load_matrix("results/posctrl_mrna_log2fc.tsv")
    tpm = load_matrix("results/posctrl_gene_tpm.tsv")
    acc = pd.read_csv("config/posctrl_accessions.tsv", sep="\t")
    ctrl_runs = list(acc.run_accession[acc.condition == "control"])
    ctrl_mean = tpm[ctrl_runs].mean(axis=1)
    expressed = set(ctrl_mean[ctrl_mean >= MIN_CTRL_TPM].index)
    mrna = mrna[mrna.index.isin(expressed)]
    base = {g: float(np.log10(v + 1.0)) for g, v in ctrl_mean.items() if g in expressed}

    mb = load_mirdb("config/mirdb_mmu.tsv", score_min=80)
    for m in TRANSFECTED:
        assert m in mb, f"{m} not in miRDB mouse"

    rng = random.Random(0)
    pool = [m for m in mb if m not in TRANSFECTED and len(mb[m]) >= 50]
    random_panel = rng.sample(pool, N_RANDOM)

    # all candidate miRNAs treated as "up" (transfected = truly overexpressed;
    # random = counterfactual "what if this one were the perturbed miRNA")
    mirs = TRANSFECTED + random_panel
    mirna_l2fc = pd.DataFrame({c: {m: 5.0 for m in mirs} for c in CELLS})

    df = run_coherence(mirna_l2fc, mrna, {}, mb, base,
                       responsive_thresh=1.0, target_mode="union",
                       min_targets=5, n_per_target=20, seed=0)
    df["group"] = df["miRNA"].apply(lambda m: "transfected" if m in TRANSFECTED else "random")
    df.to_csv(OUT, sep="\t", index=False)

    print("=== POSITIVE CONTROL (CRI on miRNA-mimic transfection) ===")
    for grp in ("transfected", "random"):
        g = df[(df.group == grp) & df["class"].isin(["productive", "passenger", "anti_coherent"])]
        prod = (g["class"] == "productive").mean() if len(g) else float("nan")
        print(f"\n{grp}: tested={len(g)}  productive={int((g['class']=='productive').sum())} "
              f"({prod:.0%})  median AUC={g['auc'].median():.3f}")
    print("\nTransfected miRNAs per cell (AUC>0.5 = targets repressed):")
    t = df[df.group == "transfected"]
    print(t[["chemical", "miRNA", "n_targets", "auc", "fdr", "cri", "class"]]
          .rename(columns={"chemical": "cell"}).to_string(index=False))


if __name__ == "__main__":
    main()
