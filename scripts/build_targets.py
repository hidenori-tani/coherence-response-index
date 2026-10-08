"""Build predicted miRNA->gene target tables for the coherence test.

PRIMARY = miRDB (mouse). miRDB is exactly the database TAP2020 used for its
predicted-target GO interpretation, so auditing against it directly tests the
originals' own method. miRDB gives miRNA(miRBase name) -> RefSeq target + score;
RefSeq is mapped to gene symbol via config/refseq2symbol.tsv (GENCODE bridge).

Outputs:
  config/mirdb_mmu.tsv      miRNA<TAB>gene_symbol<TAB>score  (load_mirdb format)
  config/targetscan_mmu.tsv header-only placeholder (TargetScan = optional robustness, deferred)

Input: ref/miRDB.txt.gz (downloaded if missing).
"""
import gzip
import urllib.request
from pathlib import Path

MIRDB_URL = "https://mirdb.org/download/miRDB_v6.0_prediction_result.txt.gz"
MIRDB_GZ = "ref/miRDB.txt.gz"


def load_refseq2symbol(path="config/refseq2symbol.tsv"):
    m = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            r, s = line.rstrip("\n").split("\t")
            m[r] = s
    return m


def main():
    Path("ref").mkdir(exist_ok=True)
    Path("config").mkdir(exist_ok=True)
    if not Path(MIRDB_GZ).exists():
        print("downloading miRDB ...")
        urllib.request.urlretrieve(MIRDB_URL, MIRDB_GZ)

    r2s = load_refseq2symbol()
    best = {}          # (miRNA, symbol) -> max score
    total_mmu = mapped = 0
    with gzip.open(MIRDB_GZ, "rt") as fh:
        for line in fh:
            mir, refseq, score = line.rstrip("\n").split("\t")
            if not mir.startswith("mmu-"):
                continue
            total_mmu += 1
            sym = r2s.get(refseq.split(".")[0])
            if sym is None:
                continue
            mapped += 1
            key = (mir, sym)
            sc = float(score)
            if key not in best or sc > best[key]:
                best[key] = sc

    with open("config/mirdb_mmu.tsv", "w") as out:
        for (mir, sym), sc in sorted(best.items()):
            out.write(f"{mir}\t{sym}\t{sc:.2f}\n")
    print(f"config/mirdb_mmu.tsv: {len(best)} miRNA-gene pairs "
          f"({mapped}/{total_mmu} mmu rows mapped to symbol)")

    ts = Path("config/targetscan_mmu.tsv")
    if not ts.exists():
        ts.write_text("miRNA\tGene\tscore\n")   # placeholder; TargetScan robustness deferred
        print("config/targetscan_mmu.tsv: header-only placeholder (TargetScan deferred)")


if __name__ == "__main__":
    main()
