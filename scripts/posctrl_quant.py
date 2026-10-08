"""Positive control (GSE269995 / PRJNA1124575): miRNA-mimic transfection RNA-seq.
3 mouse neural cell types, each control (cel-miR-67) vs a pool of three mimics
(miR-199a-3p, miR-99b-3p, miR-140-5p), 3 reps. Single-end; quantified through the
SAME kallisto mouse pipeline as the toxicogenomics data. Streams per run.

Outputs:
  results/posctrl_gene_tpm.tsv     genes x 18 runs
  results/posctrl_mrna_log2fc.tsv  genes x 3 cells (mimic vs control)
"""
import csv
import subprocess
from pathlib import Path

import pandas as pd

from _data import download, assert_space, free_gb
from tmc.io_expr import compute_log2fc_robust

ACC = "config/posctrl_accessions.tsv"
INDEX = "ref/kallisto_mmu.idx"
TX2GENE = "ref/tx2gene.tsv"
CELLS = ["bv2", "astrocyte", "neuro2a"]


def rows():
    with open(ACC) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def quant_run(run, fq):
    outdir = Path("data/kallisto_pc") / run
    outdir.mkdir(parents=True, exist_ok=True)
    subprocess.run(["kallisto", "quant", "-i", INDEX, "-o", str(outdir),
                    "--single", "-l", "200", "-s", "20", "-t", "4", str(fq)],
                   check=True, capture_output=True)
    return outdir / "abundance.tsv"


def aggregate():
    t2g = pd.read_csv(TX2GENE, sep="\t")
    sym = dict(zip(t2g["tx_id"], t2g["gene_symbol"]))
    per = {}
    for d in sorted(Path("data/kallisto_pc").glob("*/abundance.tsv")):
        a = pd.read_csv(d, sep="\t")
        a["gene"] = a["target_id"].map(sym).fillna(a["target_id"])
        per[d.parent.name] = a.groupby("gene")["tpm"].sum()
    g = pd.DataFrame(per).fillna(0.0)
    g.index.name = "gene"
    g.to_csv("results/posctrl_gene_tpm.tsv", sep="\t")
    return g


def main():
    Path("data/fastq").mkdir(parents=True, exist_ok=True)
    Path("results").mkdir(exist_ok=True)
    rs = rows()
    for r in rs:
        run = r["run_accession"]
        if (Path("data/kallisto_pc") / run / "abundance.tsv").exists():
            print(f"{run}: cached"); continue
        assert_space()
        fq = download(r["fastq_ftp"], Path("data/fastq") / f"{run}.fastq.gz", r["fastq_md5"])
        quant_run(run, fq)
        fq.unlink(missing_ok=True)
        print(f"{run} ({r['cell']}/{r['condition']}): quantified; free {free_gb():.1f} GB")
    g = aggregate()
    # per-cell log2FC: mimic mean vs control mean
    rs_by = {(r["cell"], r["condition"]): [] for r in rs}
    for r in rs:
        rs_by[(r["cell"], r["condition"])].append(r["run_accession"])
    cols = {}
    for cell in CELLS:
        treat = g[rs_by[(cell, "mimic")]].mean(axis=1).to_dict()
        ctrl = g[rs_by[(cell, "control")]].mean(axis=1).to_dict()
        cols[cell] = compute_log2fc_robust(treat, ctrl)   # robust pseudocount (matches main analysis)
    df = pd.DataFrame(cols)
    df.index.name = "feature"
    df.to_csv("results/posctrl_mrna_log2fc.tsv", sep="\t")
    print(f"wrote results/posctrl_mrna_log2fc.tsv {df.shape}")


if __name__ == "__main__":
    main()
