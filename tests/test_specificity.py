from tmc.specificity import generic_fraction, benzene_contrast


def test_generic_fraction_counts_shared_responses():
    # miRNA productive in >=k chemicals = generic
    by_chem = {"benzene": {"mir-1", "mir-2"}, "phenol": {"mir-1"}, "tce": {"mir-1"}}
    g = generic_fraction(by_chem, k=3)
    assert g["mir-1"] is True
    assert g["mir-2"] is False


def test_benzene_contrast_flags_carcinogen_specific():
    cri_by_chem = {"benzene": {"mir-9": 3.0}, "phenol": {"mir-9": 0.1}, "tce": {"mir-9": 0.0}}
    res = benzene_contrast(cri_by_chem)
    assert res["mir-9"] > 0  # higher in benzene than mean of others


def test_benzene_contrast_negative_when_lower_in_benzene():
    cri_by_chem = {"benzene": {"mir-x": 0.0}, "phenol": {"mir-x": 2.0}, "tce": {"mir-x": 2.0}}
    res = benzene_contrast(cri_by_chem)
    assert res["mir-x"] < 0


def test_generic_fraction_threshold_boundary():
    by_chem = {"a": {"m"}, "b": {"m"}, "c": set()}
    assert generic_fraction(by_chem, k=2)["m"] is True
    assert generic_fraction(by_chem, k=3)["m"] is False
