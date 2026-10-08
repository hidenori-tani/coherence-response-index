from tmc.nullmodel import matched_background

def test_matched_background_size_and_excludes_targets():
    covar = {f"G{i}": float(i) for i in range(100)}  # 1-D covariate (e.g., expression)
    targets = {"G10", "G20", "G30"}
    bg = matched_background(targets, covar, n_per_target=5, seed=1)
    assert len(bg) == 15            # 5 per target
    assert targets.isdisjoint(bg)   # never sample a target

def test_matched_background_is_covariate_close():
    covar = {f"G{i}": float(i) for i in range(1000)}
    targets = {"G500"}
    bg = matched_background(targets, covar, n_per_target=10, seed=1)
    import statistics
    assert abs(statistics.mean(covar[g] for g in bg) - 500) < 50  # near matched bin


def test_matched_background_window_bounds_the_draws():
    # the +/-window pool: a narrow window keeps every draw within `window` ranks
    covar = {f"G{i}": float(i) for i in range(1000)}
    targets = {"G500"}
    bg = matched_background(targets, covar, n_per_target=10, seed=1, window=10)
    assert all(abs(covar[g] - 500) <= 10 for g in bg)
    # the default window is 50 and is what the published analysis used
    assert matched_background(targets, covar, n_per_target=10, seed=1) == \
        matched_background(targets, covar, n_per_target=10, seed=1, window=50)
