from pathlib import Path
from tmc.config import load_accessions

# CWD-independent path to the frozen map (tests may run from repo root or pkg dir)
ACCESSIONS = str(Path(__file__).resolve().parent.parent / "config" / "accessions.tsv")


def test_accessions_complete():
    acc = load_accessions(ACCESSIONS)
    assert {r["layer"] for r in acc} == {"mirna", "mrna"}
    for layer in ("mirna", "mrna"):
        rows = [r for r in acc if r["layer"] == layer]
        assert len(rows) == 20
        chems = {r["chemical"] for r in rows}
        assert "dmso_control" in chems
        assert len(chems) == 10  # 9 chemicals + control


def test_each_chemical_has_two_replicates():
    acc = load_accessions(ACCESSIONS)
    for layer in ("mirna", "mrna"):
        rows = [r for r in acc if r["layer"] == layer]
        from collections import Counter
        counts = Counter(r["chemical"] for r in rows)
        assert all(v == 2 for v in counts.values()), counts


def test_fastq_fields_present():
    acc = load_accessions(ACCESSIONS)
    for r in acc:
        assert r["fastq_ftp"].startswith("ftp.sra.ebi.ac.uk")
        assert r["fastq_md5"]
        # mRNA is paired-end -> two files; miRNA single-end -> one
        n = len(r["fastq_ftp"].split(";"))
        assert n == (2 if r["layer"] == "mrna" else 1)
