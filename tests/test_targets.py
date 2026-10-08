from tmc.targets import load_targetscan, load_mirdb, target_set, harmonize_mirna_names

def test_load_targetscan(tmp_path):
    p = tmp_path/"ts.tsv"
    p.write_text("miRNA\tGene\tscore\nmmu-miR-1\tGeneA\t-0.5\nmmu-miR-1\tGeneB\t-0.2\n")
    d = load_targetscan(str(p))
    assert d["mmu-miR-1"] == {"GeneA", "GeneB"}

def test_load_mirdb_threshold(tmp_path):
    p = tmp_path/"mirdb.tsv"
    p.write_text("mmu-miR-1\tGeneA\t90\nmmu-miR-1\tGeneC\t50\n")
    d = load_mirdb(str(p), score_min=80)
    assert d["mmu-miR-1"] == {"GeneA"}  # GeneC filtered (score 50 < 80)

def test_target_set_union_and_intersection():
    ts = {"mmu-miR-1": {"GeneA", "GeneB"}}
    mb = {"mmu-miR-1": {"GeneB", "GeneC"}}
    assert target_set("mmu-miR-1", ts, mb, mode="union") == {"GeneA","GeneB","GeneC"}
    assert target_set("mmu-miR-1", ts, mb, mode="intersection") == {"GeneB"}

def test_harmonize_strips_version_and_lowercases_suffix():
    # miRBase ID forms must collapse so log2fc keys and target keys join
    assert harmonize_mirna_names("mmu-miR-877-5p") == "mmu-mir-877-5p"
    assert harmonize_mirna_names("mmu-miR-877-5p.1") == "mmu-mir-877-5p"
