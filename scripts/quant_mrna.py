"""mRNA/lncRNA paired-end RNA-seq -> gene-level TPM via kallisto (mouse).

STREAMING + disk guard (mRNA layer is ~48 GB total; never keep >1 run on disk):
per run -> check free space -> download _1/_2 (~2.7 GB) -> kallisto quant ->
keep small abundance.tsv -> delete fastqs -> next run. After all runs, aggregate
transcript TPM to gene level via ref/tx2gene.tsv.

Prereq: kallisto index ref/kallisto_mmu.idx + ref/tx2gene.tsv (scripts/fetch_refs.py).
"""
import subprocess
from pathlib import Path

import pandas as pd

from _data import layer_rows, download, assert_space, free_gb

INDEX = "ref/kallisto_mmu.idx"
TX2GENE = "ref/tx2gene.tsv"   # cols: tx_id, gene_symbol, biotype


def quant_run(run, fq1, fq2):
    outdir = Path("data/kallisto") / run
    outdir.mkdir(parents=True, exist_ok=True)
    subprocess.run(["kallisto", "quant", "-i", INDEX, "-o", str(outdir),
                    "-t", "4", str(fq1), str(fq2)],
                   check=True, capture_output=True)
    return outdir / "abundance.tsv"


def aggregate_to_gene():
    t2g = pd.read_csv(TX2GENE, sep="\t")
    sym = dict(zip(t2g["tx_id"], t2g["gene_symbol"]))
    bio = dict(zip(t2g["tx_id"], t2g["biotype"]))
    per_run = {}
    for d in sorted(Path("data/kallisto").glob("*/abundance.tsv")):
        run = d.parent.name
        a = pd.read_csv(d, sep="\t")
        a["gene"] = a["target_id"].map(sym).fillna(a["target_id"])
        per_run[run] = a.groupby("gene")["tpm"].sum()
    gene = pd.DataFrame(per_run).fillna(0.0)
    gene.index.name = "gene"
    gene.to_csv("results/gene_tpm.tsv", sep="\t")
    # gene -> biotype (first transcript's biotype per gene symbol)
    g2b = {}
    for tx, s in sym.items():
        g2b.setdefault(s, bio.get(tx, "NA"))
    pd.Series(g2b, name="biotype").to_csv("results/gene_biotype.tsv")
    print(f"wrote results/gene_tpm.tsv {gene.shape}")


def main():
    Path("data/fastq").mkdir(parents=True, exist_ok=True)
    Path("results").mkdir(exist_ok=True)
    for r in layer_rows("mrna"):
        run = r["run_accession"]
        if (Path("data/kallisto") / run / "abundance.tsv").exists():
            print(f"{run}: cached, skip")
            continue
        assert_space()
        urls = r["fastq_ftp"].split(";")
        md5s = r["fastq_md5"].split(";")
        fq1 = download(urls[0], Path("data/fastq") / f"{run}_1.fastq.gz", md5s[0])
        fq2 = download(urls[1], Path("data/fastq") / f"{run}_2.fastq.gz", md5s[1])
        quant_run(run, fq1, fq2)
        fq1.unlink(missing_ok=True)
        fq2.unlink(missing_ok=True)
        print(f"{run}: quantified; free now {free_gb():.1f} GB")
    aggregate_to_gene()


if __name__ == "__main__":
    main()
