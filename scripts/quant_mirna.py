"""miRNA small-RNA-seq -> mature miRNA counts (miRBase mouse) via cutadapt+bowtie2.

Streaming per run (each run is tiny, ~21 MB): download -> trim -> align -> count
-> delete fastq. Writes results/mirna_counts.tsv (rows=miRNA, cols=run).
Prereq: bowtie2 index at ref/mirbase_mmu_mature (built by scripts/fetch_refs.py).
"""
import collections
import subprocess
import sys
from pathlib import Path

import pandas as pd

from _data import layer_rows, download, assert_space

ADAPTER = "TGGAATTCTCGGGTGCCAAGG"        # Illumina TruSeq small-RNA 3' adapter
INDEX = "ref/mirbase_mmu_mature"


def quant_run(fastq, run):
    trimmed = f"data/trim/{run}.fq.gz"
    subprocess.run([sys.executable, "-m", "cutadapt", "-a", ADAPTER, "-m", "18",
                    "-M", "30", "-j", "4", "-o", trimmed, str(fastq)],
                   check=True, capture_output=True)
    sam = f"data/aln/{run}.sam"
    subprocess.run(["bowtie2", "-x", INDEX, "-U", trimmed, "--norc",
                    "-L", "18", "-N", "1", "--no-unal", "-p", "4", "-S", sam],
                   check=True, capture_output=True)
    counts = collections.Counter()
    with open(sam) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            f = line.split("\t")
            flag = int(f[1])
            if not (flag & 4) and f[2] != "*":   # mapped
                counts[f[2]] += 1
    Path(trimmed).unlink(missing_ok=True)
    Path(sam).unlink(missing_ok=True)
    return counts


def main():
    for d in ("data/fastq", "data/trim", "data/aln", "results"):
        Path(d).mkdir(parents=True, exist_ok=True)
    rows = layer_rows("mirna")
    mat, mirnas = {}, set()
    for r in rows:
        run = r["run_accession"]
        assert_space()
        fq = download(r["fastq_ftp"], Path("data/fastq") / f"{run}.fastq.gz", r["fastq_md5"])
        c = quant_run(fq, run)
        mat[run] = c
        mirnas |= set(c)
        fq.unlink(missing_ok=True)   # streaming: free immediately
        print(f"{run}: {sum(c.values())} mapped reads, {len(c)} miRNAs")
    runs = [r["run_accession"] for r in rows]
    df = pd.DataFrame({run: {m: mat[run].get(m, 0) for m in sorted(mirnas)} for run in runs})
    df.index.name = "miRNA"
    df.to_csv("results/mirna_counts.tsv", sep="\t")
    print(f"wrote results/mirna_counts.tsv {df.shape}")


if __name__ == "__main__":
    main()
