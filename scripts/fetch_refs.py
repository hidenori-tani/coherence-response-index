"""Fetch + build references (idempotent; files cached in ref/).

Builds:
  ref/kallisto_mmu.idx          kallisto index (GENCODE vM36 transcriptome)
  ref/tx2gene.tsv               tx_id(full GENCODE header) -> gene_symbol, biotype
  ref/mmu_mature.fa             miRBase mmu mature, U->T (DNA) for bowtie2
  ref/mirbase_mmu_mature.*      bowtie2 index of the above
  config/refseq2symbol.tsv      RefSeq RNA (no version) -> gene_symbol (committed; small)

Reference versions: GENCODE M36 (mouse), miRBase mature (current), via stable URLs.
"""
import gzip
import subprocess
import urllib.request
from pathlib import Path

GENCODE_TX = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/gencode.vM36.transcripts.fa.gz"
GENCODE_REFSEQ = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/gencode.vM36.metadata.RefSeq.gz"
MIRBASE_MATURE = "https://www.mirbase.org/download/mature.fa"

REF = Path("ref")


def _get(url, dest):
    dest = Path(dest)
    if not dest.exists() or dest.stat().st_size == 0:
        print(f"downloading {dest.name} ...")
        urllib.request.urlretrieve(url, dest)
    return dest


def _strip_ver(acc):
    return acc.split(".")[0]


def build_tx2gene(tx_fa):
    """tx_id = full GENCODE header token; gene_symbol = field 5; biotype = field 7.
    Also returns bare-ENSMUST -> symbol for the RefSeq bridge."""
    rows, enst2sym = [], {}
    with gzip.open(tx_fa, "rt") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            head = line[1:].strip()
            f = head.split("|")
            sym = f[5] if len(f) > 5 else head
            bio = f[7] if len(f) > 7 else "NA"
            rows.append((head, sym, bio))
            enst2sym[_strip_ver(f[0])] = sym
    with open(REF / "tx2gene.tsv", "w") as out:
        out.write("tx_id\tgene_symbol\tbiotype\n")
        for tx, sym, bio in rows:
            out.write(f"{tx}\t{sym}\t{bio}\n")
    print(f"tx2gene.tsv: {len(rows)} transcripts")
    return enst2sym


def build_refseq2symbol(meta_gz, enst2sym):
    n = 0
    with gzip.open(meta_gz, "rt") as fh, open("config/refseq2symbol.tsv", "w") as out:
        out.write("refseq\tgene_symbol\n")
        seen = set()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            enst, refseq = _strip_ver(p[0]), _strip_ver(p[1])
            sym = enst2sym.get(enst)
            if sym and refseq not in seen:
                out.write(f"{refseq}\t{sym}\n")
                seen.add(refseq)
                n += 1
    print(f"config/refseq2symbol.tsv: {n} RefSeq->symbol")


def build_mirbase_index(mature_fa):
    out_fa = REF / "mmu_mature.fa"
    n = 0
    with open(mature_fa) as fh, open(out_fa, "w") as out:
        name, seq, keep = None, [], False
        def flush():
            nonlocal n
            if keep and seq:
                out.write(f">{name}\n{''.join(seq).replace('U', 'T').replace('u', 't')}\n")
                n += 1
        for line in fh:
            if line.startswith(">"):
                flush()
                name = line[1:].split()[0]      # e.g. mmu-let-7g-5p
                seq, keep = [], name.startswith("mmu-")
            else:
                seq.append(line.strip())
        flush()
    print(f"mmu_mature.fa: {n} mature miRNAs")
    subprocess.run(["bowtie2-build", "--quiet", str(out_fa), str(REF / "mirbase_mmu_mature")],
                   check=True)
    print("bowtie2 index built: ref/mirbase_mmu_mature")


def main():
    REF.mkdir(exist_ok=True)
    Path("config").mkdir(exist_ok=True)
    tx = _get(GENCODE_TX, REF / "transcripts.fa.gz")
    meta = _get(GENCODE_REFSEQ, REF / "metaRefSeq.gz")
    mature = _get(MIRBASE_MATURE, REF / "mature.fa")

    if not (REF / "kallisto_mmu.idx").exists():
        print("building kallisto index (a few minutes) ...")
        subprocess.run(["kallisto", "index", "-i", str(REF / "kallisto_mmu.idx"), str(tx)],
                       check=True)
    enst2sym = build_tx2gene(tx)
    build_refseq2symbol(meta, enst2sym)
    build_mirbase_index(mature)
    print("references ready.")


if __name__ == "__main__":
    main()
