import pandas as pd
from tmc.pipeline import run_coherence, summarize_specificity

GENES = [f"G{i}" for i in range(12)]
TARGETS = {f"G{i}" for i in range(6)}        # G0..G5 are targets of mmu-miR-up
NONTARGETS = [f"G{i}" for i in range(6, 12)]  # G6..G11


def _fixtures():
    mirna = pd.DataFrame(
        {"benzene": {"mmu-miR-up": 3.0, "mmu-miR-flat": 0.2},
         "phenol":  {"mmu-miR-up": 3.0, "mmu-miR-flat": 0.1}}
    )
    # benzene: targets strongly DOWN (coherent w/ up miRNA), non-targets ~flat
    benzene = {g: (-2.5 if g in TARGETS else v)
               for g, v in zip(GENES, [0]*6 + [0.1, -0.1, 0.0, 0.05, -0.05, 0.02])}
    # phenol: targets and non-targets both ~flat (no coherence -> passenger)
    phenol = dict(zip(GENES, [0.1, -0.1, 0.0, 0.2, -0.2, 0.05,
                              0.0, 0.15, -0.15, 0.1, -0.1, 0.05]))
    mrna = pd.DataFrame({"benzene": benzene, "phenol": phenol})
    targetscan = {"mmu-miR-up": set(TARGETS)}
    mirdb = {}
    baseline = {g: float(i) for i, g in enumerate(GENES)}
    return mirna, mrna, targetscan, mirdb, baseline


def test_run_coherence_productive_vs_passenger():
    mirna, mrna, ts, mb, base = _fixtures()
    df = run_coherence(mirna, mrna, ts, mb, base, responsive_thresh=1.0,
                       min_targets=5, n_per_target=20, seed=1)
    # mmu-miR-flat is below threshold -> excluded entirely
    assert set(df["miRNA"]) == {"mmu-miR-up"}
    bz = df[(df.chemical == "benzene") & (df.miRNA == "mmu-miR-up")].iloc[0]
    ph = df[(df.chemical == "phenol") & (df.miRNA == "mmu-miR-up")].iloc[0]
    assert bz["class"] == "productive"
    assert bz["cri"] > 0            # up miRNA + repressed targets -> coherent
    assert bz["median_diff"] < 0
    assert ph["class"] == "passenger"
    assert "fdr" in df.columns


def test_run_coherence_min_targets_gate():
    mirna, mrna, ts, mb, base = _fixtures()
    df = run_coherence(mirna, mrna, ts, mb, base, responsive_thresh=1.0,
                       min_targets=99, seed=1)  # nothing has 99 targets
    assert (df["class"] == "untested").all()
    assert df["p"].isna().all()


def test_anti_coherent_down_mirna_with_repressed_targets():
    # a DOWN miRNA whose targets are also repressed = significant but wrong
    # direction -> anti_coherent, NOT productive
    mirna = pd.DataFrame({"benzene": {"mmu-miR-dn": -2.5}})
    benzene = {g: (-2.5 if g in TARGETS else v)
               for g, v in zip(GENES, [0]*6 + [0.1, -0.1, 0.0, 0.05, -0.05, 0.02])}
    mrna = pd.DataFrame({"benzene": benzene})
    ts = {"mmu-miR-dn": set(TARGETS)}
    base = {g: float(i) for i, g in enumerate(GENES)}
    df = run_coherence(mirna, mrna, ts, {}, base, responsive_thresh=1.0,
                       min_targets=5, seed=1)
    row = df.iloc[0]
    assert row["cri"] < 0                 # down miRNA + repressed targets = anti-coherent
    assert row["class"] == "anti_coherent"


def test_summarize_specificity():
    mirna, mrna, ts, mb, base = _fixtures()
    df = run_coherence(mirna, mrna, ts, mb, base, responsive_thresh=1.0,
                       min_targets=5, seed=1)
    spec = summarize_specificity(df, k=2, benzene_key="benzene")
    row = spec[spec.miRNA == "mmu-miR-up"].iloc[0]
    assert row["n_chemicals_productive"] == 1      # productive only in benzene
    assert not row["generic"]                        # <2 chemicals (numpy bool)
    assert row["benzene_contrast"] > 0              # higher CRI in benzene
