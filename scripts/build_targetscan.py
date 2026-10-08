"""Build the TargetScan Mouse 8.0 target tables from the official download files.

Inputs (https://www.targetscan.org/mmu_80/mmu_80_data_download/, fetched 2026-09-24;
sha256 in data/targetscan_mmu80/SHA256SUMS.txt):
  data/targetscan_mmu80/Predicted_Targets_Info.default_predictions.txt
  data/targetscan_mmu80/miR_Family_Info.txt
Outputs:
  config/targetscan_mmu.tsv         miR family -> gene (mouse, Species ID 10090).
                                    This is the table used in the submitted version
                                    (92,817 pairs, 218 families); rebuilding it here
                                    must reproduce it exactly.
  config/targetscan_mmu_mature.tsv  mature miRNA (miRBase ID) -> gene, through the
                                    official miRBase-ID -> family map. This is the
                                    lookup the data need: the data are named by mature
                                    miRNA (mmu-miR-146b-5p), TargetScan by family
                                    (miR-146-5p).
"""
import csv
from pathlib import Path

SRC = Path("data/targetscan_mmu80")
MOUSE = "10090"


def family_targets():
    pairs = set()
    with open(SRC / "Predicted_Targets_Info.default_predictions.txt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["Species ID"] == MOUSE:
                pairs.add((r["miR Family"], r["Gene Symbol"]))
    return pairs


def mature_to_family():
    m = {}
    with open(SRC / "miR_Family_Info.txt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["Species ID"] == MOUSE:
                m[r["MiRBase ID"]] = r["miR family"]
    return m


def main():
    pairs = family_targets()
    fam = {}
    for f, g in pairs:
        fam.setdefault(f, set()).add(g)
    rows = sorted(("mmu-" + f, g) for f, g in pairs)
    old = Path("config/targetscan_mmu.tsv")
    if old.exists():
        with open(old) as fh:
            prev = {(r["miRNA"], r["Gene"]) for r in csv.DictReader(fh, delimiter="\t")}
        if prev != set(rows):
            raise SystemExit("ABORT: rebuilt family table differs from config/targetscan_mmu.tsv")
        print(f"config/targetscan_mmu.tsv reproduced exactly: {len(rows)} pairs, {len(fam)} families")
    m2f = mature_to_family()
    out = sorted((mat, g) for mat, f in m2f.items() for g in fam.get(f, ()))
    with open("config/targetscan_mmu_mature.tsv", "w") as fh:
        fh.write("miRNA\tGene\tfamily\n")
        for mat, g in out:
            fh.write(f"{mat}\t{g}\t{m2f[mat]}\n")
    with_targets = sum(1 for f in m2f.values() if f in fam)
    print(f"config/targetscan_mmu_mature.tsv: {len(out)} pairs; {len(m2f)} mouse mature miRNAs, "
          f"{with_targets} of them in a family with default predictions")
    # positive control of the map: a known member resolves to its family
    assert m2f.get("mmu-miR-146b-5p") == "miR-146-5p"
    assert m2f.get("mmu-let-7a-5p") == "let-7-5p/miR-98-5p"


if __name__ == "__main__":
    main()
