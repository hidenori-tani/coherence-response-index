from tmc.coherence import directional_shift, cri, classify_response


def test_directional_shift_detects_repression():
    # targets clearly more negative than background -> AUC>0.5, p small, neg median diff
    target = [-3, -2.5, -2, -2.2, -1.8]
    background = [0.1, -0.2, 0.0, 0.3, -0.1, 0.2]
    s = directional_shift(target, background)
    assert s["median_diff"] < 0
    assert s["auc"] > 0.5
    assert s["p"] < 0.05


def test_directional_shift_null_when_same():
    target = [0.1, -0.1, 0.0]
    background = [0.1, -0.1, 0.0, 0.05, -0.05]
    s = directional_shift(target, background)
    assert s["p"] > 0.05


def test_directional_shift_empty_is_nan():
    s = directional_shift([], [1.0, 2.0])
    assert s["n"] == 0


def test_cri_sign_follows_coherence():
    # up miRNA (log2fc>0) with repressed targets -> productive -> positive CRI
    up = cri(mirna_log2fc=2.0, target_median_diff=-1.5, auc=0.8)
    down = cri(mirna_log2fc=2.0, target_median_diff=0.0, auc=0.5)
    assert up > down
    assert up > 0


def test_cri_anticoherent_is_negative():
    # up miRNA but targets go UP (de-repressed) -> anti-coherent -> negative CRI
    val = cri(mirna_log2fc=2.0, target_median_diff=1.2, auc=0.2)
    assert val < 0


def test_cri_down_mirna_coherent_when_targets_up():
    # down miRNA (log2fc<0) expects targets de-repressed (positive diff, auc<0.5)
    val = cri(mirna_log2fc=-2.0, target_median_diff=1.0, auc=0.2)
    assert val > 0


def test_classify_response():
    assert classify_response(p=0.001, auc=0.8, mirna_log2fc=2.0) == "productive"
    assert classify_response(p=0.6, auc=0.52, mirna_log2fc=2.0) == "passenger"
