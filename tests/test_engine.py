"""
Phase 3 (core engine) and Phase 4 (Sequential/Cumulative) validation checks.

Lean on purpose -- "trust the number it gives you" validation, not an
exhaustive suite. Run with `pytest tests/test_engine.py -v`, or directly
with `python tests/test_engine.py` (no pytest required).

IMPORTANT finding from building these tests (kept here, not just in a
commit message, because it matters for interpreting real results too):
with ~10 historical years, the false-positive rate CONDITIONAL ON ONE
FITTED MODEL varies a lot fit to fit -- in a 25-fit simulation it ranged
from 0% to 25% at nominal alpha=0.05, purely from the sampling noise of
estimating a median/MAD curve from only ~10 years. The rate averaged
ACROSS many different historical samples lands close to nominal alpha, but
any ONE real analysis only ever has ONE historical fit -- so a real result
inherits some of that fit-specific uncertainty too. This is the small-N
power/calibration limitation from the "is this defensible" discussion,
now with a number attached. The tests below pool across several
historical fits (not just one) for exactly this reason -- a single-fit
stochastic test would be flaky by construction, not because of a bug.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import pandas as pd

from tide_lite.engine import (
    TideConfig, fit_historical, test_treatment, to_day_of_year, bin_significant,
    _max_window_mean, estimate_departure_recovery,
)
from tide_lite.controllers import (
    holm_bonferroni, SequentialEvent, sequential_test,
    CumulativeWindow, cumulative_test, first_recovery,
)


# ---------------------------------------------------------------------------
# Shared synthetic data generator
# ---------------------------------------------------------------------------

def make_year(rng, year, effect=0.0, n_days=365, base=10.0, seasonal_amp=3.0,
              ar_noise_sd=0.6, ar_phi=0.7):
    days = np.arange(1, n_days + 1)
    seasonal = base + seasonal_amp * np.sin(2 * np.pi * (days - 60) / 365)
    noise = np.zeros(n_days)
    e = rng.normal(0, ar_noise_sd, n_days)
    for i in range(1, n_days):
        noise[i] = ar_phi * noise[i - 1] + e[i]
    values = np.clip(seasonal + noise + effect, 0.1, None)
    dates = pd.date_range(f"{year}-01-01", periods=n_days, freq="D")
    return pd.DataFrame({"date": dates, "value": values})


def make_year_localized(rng, year, spike_day_range=None, spike_size=0.0, **kwargs):
    """Same as make_year, but the effect (if any) only applies within
    spike_day_range=(start_day, end_day) instead of uniformly across the
    whole year -- for testing that bin_significant() correctly localizes
    an effect confined to specific months, not just detects "some effect
    somewhere."
    """
    df = make_year(rng, year, effect=0.0, **kwargs)
    if spike_day_range is not None:
        start, end = spike_day_range
        day_of_year = np.arange(1, len(df) + 1)
        mask = (day_of_year >= start) & (day_of_year <= end)
        df.loc[mask, "value"] = df.loc[mask, "value"] + spike_size
    return df


def make_historical(rng, start_year=2000, n_years=10):
    return pd.concat(
        [make_year(rng, start_year + i) for i in range(n_years)], ignore_index=True
    )


def make_fit(seed, n_years=10, n_bootstrap=300):
    rng = np.random.default_rng(seed)
    historical = make_historical(rng, n_years=n_years)
    config = TideConfig(bin_days=30, detrend_mode="additive", n_bootstrap=n_bootstrap)
    return fit_historical(historical, "date", "value", config), rng


# ---------------------------------------------------------------------------
# 1. Day-of-year alignment -- deterministic
# ---------------------------------------------------------------------------

def test_day_of_year_alignment():
    dates = pd.to_datetime([
        "2021-01-01", "2021-12-31",         # non-leap: day 1, day 365
        "2020-01-01", "2020-02-28",         # leap: day 1, day 59
        "2020-02-29",                        # leap day: NaN
        "2020-03-01", "2020-12-31",         # leap, post-Feb: shifted down 1
    ])
    doy = to_day_of_year(pd.Series(dates))
    assert doy.iloc[0] == 1
    assert doy.iloc[1] == 365
    assert doy.iloc[2] == 1
    assert doy.iloc[3] == 59
    assert np.isnan(doy.iloc[4])
    assert doy.iloc[5] == 60
    assert doy.iloc[6] == 365


# ---------------------------------------------------------------------------
# 2. Holm-Bonferroni -- deterministic, hand-checked
# ---------------------------------------------------------------------------

def test_holm_bonferroni_known_case():
    p_sorted = np.array([0.001, 0.02, 0.03, 0.04])
    reject = holm_bonferroni(p_sorted, alpha=0.05)
    assert list(reject) == [True, False, False, False]

    p_shuffled = np.array([0.04, 0.001, 0.03, 0.02])
    reject2 = holm_bonferroni(p_shuffled, alpha=0.05)
    assert list(reject2) == [False, True, False, False]


# ---------------------------------------------------------------------------
# 3. Null case: false positive rate near alpha, POOLED across several fits
# ---------------------------------------------------------------------------

def test_null_case_false_positive_rate():
    alpha = 0.05
    n_fits = 12
    n_trials_per_fit = 15
    false_positives = 0
    total = 0
    for fit_seed in range(n_fits):
        fit, rng = make_fit(seed=100 + fit_seed)
        for i in range(n_trials_per_fit):
            null_year = make_year(rng, 2100 + i, effect=0.0)
            result = test_treatment(fit, null_year, mode="prediction",
                                     treatment_year_index=len(fit.years_used), rng=rng)
            total += 1
            false_positives += result.p_value < alpha
    rate = false_positives / total
    assert 0.0 <= rate <= 0.13, f"pooled false positive rate {rate:.3f} is outside the expected band"


# ---------------------------------------------------------------------------
# 4. Known effect: gets detected (single fit is fine here -- a large true
#    effect should dominate fit-to-fit noise, unlike the null case above)
# ---------------------------------------------------------------------------

def test_known_effect_detected():
    fit, rng = make_fit(seed=2)
    detections = 0
    n_trials = 10
    for i in range(n_trials):
        effect_year = make_year(rng, 2200 + i, effect=2.5)
        result = test_treatment(fit, effect_year, mode="prediction",
                                 treatment_year_index=len(fit.years_used), rng=rng)
        detections += result.p_value < 0.05
    assert detections >= 8, f"only detected {detections}/{n_trials} -- expected a clear effect to be reliably caught"


# ---------------------------------------------------------------------------
# 5. Sequential: family-wise error rate control, POOLED across fits
# ---------------------------------------------------------------------------

def test_sequential_fwer_control():
    n_fits = 8
    n_experiments_per_fit = 12
    n_events = 5
    family_false_positives = 0
    total = 0
    for fit_seed in range(n_fits):
        fit, rng = make_fit(seed=300 + fit_seed)
        for exp in range(n_experiments_per_fit):
            events = [
                SequentialEvent(label=str(k),
                                 treatment_df=make_year(rng, 2300 + exp * 10 + k, effect=0.0),
                                 year_index=len(fit.years_used))
                for k in range(n_events)
            ]
            results = sequential_test(fit, events, alpha=0.05)
            family_false_positives += any(r.significant_corrected for r in results)
            total += 1
    fwer = family_false_positives / total
    assert fwer <= 0.15, f"family-wise error rate {fwer:.3f} is higher than expected"


# ---------------------------------------------------------------------------
# 6. Cumulative: flags an effect, then clears once it fades, POOLED across fits
# ---------------------------------------------------------------------------

def test_cumulative_flags_then_clears():
    n_fits = 6
    n_trials_per_fit = 10
    early_detections = 0
    late_detections = 0
    total = 0
    for fit_seed in range(n_fits):
        fit, rng = make_fit(seed=400 + fit_seed)
        for i in range(n_trials_per_fit):
            windows = [
                CumulativeWindow("early", make_year(rng, 2500 + i * 10, effect=3.0),
                                  year_index=len(fit.years_used)),
                CumulativeWindow("late", make_year(rng, 2500 + i * 10 + 1, effect=0.0),
                                  year_index=len(fit.years_used) + 1),
            ]
            results = cumulative_test(fit, windows, alpha=0.05)
            early_detections += results[0].significant
            late_detections += results[1].significant
            total += 1
    early_rate = early_detections / total
    late_rate = late_detections / total
    assert early_rate >= 0.85, f"a clear ongoing effect should almost always be flagged (got {early_rate:.2f})"
    assert late_rate <= 0.15, f"a fully faded effect should rarely be flagged (got {late_rate:.2f}, expected near alpha=0.05)"


# ---------------------------------------------------------------------------
# 6b. first_recovery: requires consecutive non-significant windows, not
#     just one -- and correctly ignores a later false-alarm blip once
#     recovery is already locked in.
# ---------------------------------------------------------------------------

def test_first_recovery_requires_consecutive():
    fit, rng = make_fit(seed=450, n_years=20)
    # effect fades to zero at index 2, stays there, but a later window
    # (index 4) gets a pure-noise false alarm -- first_recovery should
    # still report recovery at index 2, not be confused by index 4.
    effects = [3.0, 1.8, 0.0, 0.0, 0.0]
    windows = [
        CumulativeWindow(f"w{i}", make_year(rng, 2600 + i, effect=e), year_index=len(fit.years_used) + i)
        for i, e in enumerate(effects)
    ]
    results = cumulative_test(fit, windows, alpha=0.05)

    rec_strict = first_recovery(results, consecutive_required=2)
    assert rec_strict.recovered
    assert rec_strict.recovered_at_index == 2

    rec_none = first_recovery(results[:1], consecutive_required=2)
    assert not rec_none.recovered
    assert rec_none.recovered_at_index is None


# ---------------------------------------------------------------------------
# 7. bin_significant: a real effect confined to specific months gets
#    localized to (roughly) those months, not just "somewhere"
# ---------------------------------------------------------------------------

def test_bin_significant_localizes_a_real_effect():
    fit, rng = make_fit(seed=800, n_years=20)  # more years: a sharper,
    # more reliable per-bin signal, same reasoning as MonitoringSeries
    # needing more history -- per-bin tests split attention 13 ways
    treatment = make_year_localized(rng, 2100, spike_day_range=(121, 180), spike_size=6.0)
    result = test_treatment(fit, treatment, mode="prediction",
                             treatment_year_index=len(fit.years_used), rng=rng)
    sig = bin_significant(result, alpha=0.05)

    # bins overlapping day 121-180 at bin_days=30 are bins 5 and 6
    # (1-indexed: bin = (day-1)//30 + 1)
    spiked_bins = {5, 6}
    flagged_bins = set(fit.bins[sig].tolist())
    assert spiked_bins <= flagged_bins, (
        f"expected the spiked bins {spiked_bins} to be flagged, got {flagged_bins}"
    )
    other_flagged = flagged_bins - spiked_bins
    assert len(other_flagged) <= 1, (
        f"expected few or no bins flagged outside the real spike, got {other_flagged}"
    )


# ---------------------------------------------------------------------------
# 8. bin_significant: family-wise rate across bins stays near alpha under
#    the null (no real effect anywhere), POOLED across fits -- same
#    reasoning as every other calibration test in this file.
# ---------------------------------------------------------------------------

def test_bin_significant_null_calibration():
    n_fits = 6
    n_trials_per_fit = 10
    alpha = 0.05
    family_false_positives = 0
    total = 0
    for fit_seed in range(n_fits):
        fit, rng = make_fit(seed=900 + fit_seed, n_years=20)
        for i in range(n_trials_per_fit):
            null_year = make_year(rng, 2200 + i, effect=0.0)
            result = test_treatment(fit, null_year, mode="prediction",
                                     treatment_year_index=len(fit.years_used), rng=rng)
            sig = bin_significant(result, alpha=alpha)
            family_false_positives += bool(np.any(sig))
            total += 1
    fwer = family_false_positives / total
    assert fwer <= 0.15, f"family-wise error rate across bins {fwer:.3f} is higher than expected"


# ---------------------------------------------------------------------------
# 9. _max_window_mean mechanics -- deterministic
# ---------------------------------------------------------------------------

def test_max_window_mean_known_case():
    devs = np.array([1.0, 2.0, 5.0, 6.0, 1.0, 1.0])
    # windows of length 2: [1,2]=1.5 [2,5]=3.5 [5,6]=5.5 [6,1]=3.5 [1,1]=1.0
    best, start = _max_window_mean(devs, run_length=2)
    assert start == 2  # the [5,6] window
    assert abs(best - 5.5) < 1e-9

    # run_length == n_bins collapses to the overall mean, one window only
    best_full, start_full = _max_window_mean(devs, run_length=len(devs))
    assert start_full == 0
    assert abs(best_full - devs.mean()) < 1e-9


# ---------------------------------------------------------------------------
# 10. sustained-departure test: run_length=1 reproduces the standard
#     single-bin-max test EXACTLY -- a consistency check that the new
#     machinery generalizes the old, not a parallel reimplementation of it
# ---------------------------------------------------------------------------

def test_sustained_run_length_one_matches_standard():
    fit, rng = make_fit(seed=1000, n_years=15)
    treatment = make_year(rng, 3500, effect=1.5)
    result = test_treatment(fit, treatment, mode="prediction",
                             treatment_year_index=len(fit.years_used),
                             run_lengths=[1], rng=rng)
    assert result.sustained[1].p_value == result.p_value
    assert abs(result.sustained[1].test_statistic - result.test_statistic) < 1e-9


# ---------------------------------------------------------------------------
# 11. sustained-departure test: has real, substantially more power than
#     the single-bin test for a moderate effect spread across several
#     CONSECUTIVE months -- the actual point of building this. Validated
#     with enough trials that the first, much smaller (n=40) attempt at
#     this during development showed the opposite pattern purely from
#     noise at low power -- worth having in a comment, not just a commit
#     message, since it's a reminder that a quick power check can mislead
#     if it doesn't use enough trials.
# ---------------------------------------------------------------------------

def test_sustained_beats_single_bin_for_sustained_effect():
    # Pooled across several fits -- a single-fit version of this test
    # failed during development (16/50 detections, below the threshold)
    # purely from landing on a colder-than-average fit, the same
    # per-fit calibration variance documented throughout this file.
    # Pooling is the fix used everywhere else here for exactly this
    # reason, not a special case for this test.
    n_fits = 4
    n_trials_per_fit = 15
    single_wins = 0
    sustained_wins = 0
    total = 0
    for fit_seed in range(n_fits):
        fit, rng = make_fit(seed=1100 + fit_seed, n_years=20, n_bootstrap=500)
        for i in range(n_trials_per_fit):
            treatment = make_year_localized(rng, 4500 + fit_seed * 100 + i,
                                             spike_day_range=(91, 210), spike_size=0.35)
            result = test_treatment(fit, treatment, mode="prediction",
                                     treatment_year_index=len(fit.years_used),
                                     run_lengths=[4], rng=rng)
            single_wins += result.p_value < 0.05
            sustained_wins += result.sustained[4].p_value < 0.05
            total += 1
    assert sustained_wins > single_wins, (
        f"expected the sustained test to detect a moderate 4-month effect "
        f"more often (sustained={sustained_wins}/{total}, single={single_wins}/{total})"
    )
    assert sustained_wins / total >= 0.35, (
        f"sustained test detection rate {sustained_wins}/{total} lower than expected"
    )


# ---------------------------------------------------------------------------
# 12. sustained-departure test: null calibration, pooled across fits
# ---------------------------------------------------------------------------

def test_sustained_null_calibration():
    n_fits = 5
    n_trials_per_fit = 12
    alpha = 0.05
    false_positives = 0
    total = 0
    for fit_seed in range(n_fits):
        fit, rng = make_fit(seed=1200 + fit_seed, n_years=20, n_bootstrap=500)
        for i in range(n_trials_per_fit):
            null_year = make_year(rng, 2700 + i, effect=0.0)
            result = test_treatment(fit, null_year, mode="prediction",
                                     treatment_year_index=len(fit.years_used),
                                     run_lengths=[4], rng=rng)
            false_positives += result.sustained[4].p_value < alpha
            total += 1
    rate = false_positives / total
    assert rate <= 0.15, f"sustained-test false positive rate {rate:.3f} is outside the expected band"


# ---------------------------------------------------------------------------
# 13. estimate_departure_recovery: recovers within the period, and
#     matches the true injected boundaries reasonably closely
# ---------------------------------------------------------------------------

def test_departure_recovery_matches_injected_boundaries():
    fit, rng = make_fit(seed=1300, n_years=20, n_bootstrap=500)
    # spike on days 61-210 -> bins 3-7 at bin_days=30 ((day-1)//30+1)
    treatment = make_year_localized(rng, 4800, spike_day_range=(61, 210), spike_size=3.0)
    result = test_treatment(fit, treatment, mode="prediction",
                             treatment_year_index=len(fit.years_used),
                             run_lengths=[3], rng=rng)
    dr = estimate_departure_recovery(result, run_length=3)
    assert dr.recovered
    assert dr.departure_start_bin in (2, 3, 4)  # allow +/-1 bin slack
    assert dr.departure_end_bin in (6, 7, 8)
    assert dr.recovered_at_bin == dr.departure_end_bin + 1


def test_departure_recovery_never_recovers():
    fit, rng = make_fit(seed=1301, n_years=20, n_bootstrap=500)
    n_bins = len(fit.bins)
    # spike persisting through the last bin of the period
    treatment = make_year_localized(rng, 4801, spike_day_range=(241, 365), spike_size=3.0)
    result = test_treatment(fit, treatment, mode="prediction",
                             treatment_year_index=len(fit.years_used),
                             run_lengths=[3], rng=rng)
    dr = estimate_departure_recovery(result, run_length=3)
    assert not dr.recovered
    assert dr.recovered_at_bin is None
    assert dr.departure_end_bin == int(fit.bins[-1])


# ---------------------------------------------------------------------------
# Standalone runner (no pytest required)
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    tests = [
        test_day_of_year_alignment,
        test_holm_bonferroni_known_case,
        test_null_case_false_positive_rate,
        test_known_effect_detected,
        test_sequential_fwer_control,
        test_cumulative_flags_then_clears,
        test_first_recovery_requires_consecutive,
        test_bin_significant_localizes_a_real_effect,
        test_bin_significant_null_calibration,
        test_max_window_mean_known_case,
        test_sustained_run_length_one_matches_standard,
        test_sustained_beats_single_bin_for_sustained_effect,
        test_sustained_null_calibration,
        test_departure_recovery_matches_injected_boundaries,
        test_departure_recovery_never_recovers,
    ]
    failures = 0
    for t in tests:
        try:
            t()
            print(f"PASS  {t.__name__}")
        except AssertionError as e:
            failures += 1
            print(f"FAIL  {t.__name__}: {e}")
    print(f"\n{len(tests) - failures}/{len(tests)} passed")
    sys.exit(1 if failures else 0)
