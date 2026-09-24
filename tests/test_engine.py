"""
Statistical validation of the method: does it flag normal years at about
the stated rate, detect real changes, and do its parts agree with each
other? Run with `pytest tests/test_engine.py -v`.

Stochastic checks are POOLED across several simulated histories, never
one: with ~10 historical years the false-alarm rate conditional on one
particular history varies a lot from history to history (roughly 0-15% at
a nominal 5%), even though it averages close to 5%. A single-history test
would be flaky by construction. The full-size calibration study behind the
USER_GUIDE tables is tests/calibration_study.py.
"""

import os
import sys
import warnings

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import pandas as pd
from scipy.signal import lfilter

from tide_lite.engine import (
    TideConfig, fit_historical, test_treatment, to_day_of_year, bin_significant,
    critical_value, significance_band, _max_window_mean, estimate_departure_recovery,
    n_bins_for, _window_scores,
)
from tide_lite.controllers import (
    holm_bonferroni, SequentialEvent, sequential_test,
    CumulativeWindow, cumulative_test, first_recovery,
)

warnings.simplefilter("ignore", UserWarning)


# ---------------------------------------------------------------------------
# Shared synthetic data: seasonal cycle + persistent day-to-day noise
# (AR(1)) + optional whole-year shift ("wet/dry year").
# ---------------------------------------------------------------------------

def make_year(rng, year, effect=0.0, n_days=365, base=10.0, seasonal_amp=3.0,
              ar_noise_sd=0.6, ar_phi=0.7, year_sd=0.0):
    days = np.arange(1, n_days + 1)
    seasonal = base + seasonal_amp * np.sin(2 * np.pi * (days - 60) / 365)
    noise = lfilter([1.0], [1.0, -ar_phi], rng.normal(0, ar_noise_sd, n_days))
    values = np.clip(seasonal + noise + rng.normal(0, year_sd) + effect, 0.1, None)
    dates = pd.date_range(f"{year}-01-01", periods=n_days, freq="D")
    return pd.DataFrame({"date": dates, "value": values})


def make_year_localized(rng, year, spike_day_range=None, spike_size=0.0, **kwargs):
    """make_year with the effect confined to days spike_day_range=(start, end)."""
    df = make_year(rng, year, **kwargs)
    if spike_day_range is not None:
        start, end = spike_day_range
        day = np.arange(1, len(df) + 1)
        df.loc[(day >= start) & (day <= end), "value"] += spike_size
    return df


def make_historical(rng, start_year=2000, n_years=10, **kwargs):
    return pd.concat([make_year(rng, start_year + i, **kwargs) for i in range(n_years)],
                     ignore_index=True)


def make_fit(seed, n_years=10, n_bootstrap=300, bin_days=30, **kwargs):
    rng = np.random.default_rng(seed)
    historical = make_historical(rng, n_years=n_years, **kwargs)
    config = TideConfig(bin_days=bin_days, detrend_mode="additive", n_bootstrap=n_bootstrap)
    return fit_historical(historical, "date", "value", config), rng


def next_year(fit, i=0):
    """A calendar year after the historical record (the default trend basis)."""
    return fit.years_used[-1] + 1 + i


# ---------------------------------------------------------------------------
# Deterministic pieces
# ---------------------------------------------------------------------------

def test_day_of_year_alignment():
    dates = pd.to_datetime(["2021-01-01", "2021-12-31", "2020-01-01", "2020-02-28",
                            "2020-02-29", "2020-03-01", "2020-12-31"])
    doy = to_day_of_year(pd.Series(dates))
    assert list(doy.iloc[[0, 1, 2, 3, 5, 6]]) == [1, 365, 1, 59, 60, 365]
    assert np.isnan(doy.iloc[4])


def test_short_leftover_bin_is_merged():
    assert n_bins_for("month") == 12
    assert n_bins_for(30) == 12      # 12 x 30 days + 5 left over -> merged
    assert n_bins_for(7) == 52       # 52 x 7 + 1 -> merged
    assert n_bins_for(73) == 5       # divides exactly
    assert n_bins_for(100) == 4      # 65 left over is >= half a bin -> kept


def test_holm_bonferroni_known_case():
    assert list(holm_bonferroni(np.array([0.001, 0.02, 0.03, 0.04]), 0.05)) == [True, False, False, False]
    assert list(holm_bonferroni(np.array([0.04, 0.001, 0.03, 0.02]), 0.05)) == [False, True, False, False]
    assert list(holm_bonferroni(np.array([0.01, 0.02, 0.025]), 0.05)) == [True, True, True]


def test_max_window_mean_known_case():
    devs = np.array([1.0, 2.0, 5.0, 6.0, 1.0, 1.0])
    best, start = _max_window_mean(devs, run_length=2)
    assert start == 2 and abs(best - 5.5) < 1e-9
    best_full, start_full = _max_window_mean(devs, run_length=len(devs))
    assert start_full == 0 and abs(best_full - devs.mean()) < 1e-9


# ---------------------------------------------------------------------------
# Calibration: normal years are flagged at about alpha (pooled)
# ---------------------------------------------------------------------------

def _pooled_false_alarm_rate(n_fits, n_trials, seed0, alpha=0.05, **kwargs):
    hits = total = 0
    for f in range(n_fits):
        fit, rng = make_fit(seed=seed0 + f, **kwargs)
        for i in range(n_trials):
            r = test_treatment(fit, make_year(rng, next_year(fit, i), year_sd=kwargs.get("year_sd", 0.0)),
                               rng=rng)
            hits += r.p_value <= alpha
            total += 1
    return hits / total


def test_null_case_false_positive_rate():
    rate = _pooled_false_alarm_rate(12, 15, seed0=100)
    assert rate <= 0.10, f"pooled false positive rate {rate:.3f} is too high"


def test_whole_year_shifts_do_not_inflate_false_alarms():
    """The regime the annual-level component exists for: strong wet/dry-year
    shifts with a 10-year record. The previous method flagged ~16% of normal
    years here at a nominal 5%."""
    rate = _pooled_false_alarm_rate(15, 12, seed0=150, year_sd=1.2, n_bootstrap=500)
    assert rate <= 0.10, f"false positive rate {rate:.3f} with whole-year shifts"


def test_known_effect_detected():
    fit, rng = make_fit(seed=2)
    hits = sum(test_treatment(fit, make_year(rng, next_year(fit, i), effect=2.5), rng=rng).p_value <= 0.05
               for i in range(10))
    assert hits >= 8, f"only detected {hits}/10"


# ---------------------------------------------------------------------------
# The p-value, the per-bin flags and the plotted band are one decision
# ---------------------------------------------------------------------------

def test_band_flags_and_p_value_agree_exactly():
    fit, rng = make_fit(seed=3, n_years=15, n_bootstrap=400)
    checked_sig = checked_not = 0
    for i, effect in enumerate([0.0, 0.8, 1.5, 3.0] * 3):
        for alt in ("two-sided", "greater", "less"):
            r = test_treatment(fit, make_year(rng, next_year(fit, i), effect=effect),
                               alternative=alt, rng=rng)
            for alpha in (0.01, 0.05, 0.1):
                flags = bin_significant(r, alpha)
                band = significance_band(fit, r, alpha)
                outside = (r.treatment_curve > band["upper"] * (1 + 1e-12)) | \
                          (r.treatment_curve < band["lower"] * (1 - 1e-12))
                outside &= ~r.missing_bins
                assert np.array_equal(flags, outside), "band and bin flags disagree"
                assert flags.any() == (r.p_value <= alpha * (1 + 1e-9)), "flags and p-value disagree"
                c = critical_value(r, alpha)
                assert (r.test_statistic > c) == (r.p_value <= alpha * (1 + 1e-9))
                checked_sig += flags.any()
                checked_not += not flags.any()
    assert checked_sig > 0 and checked_not > 0


def test_bin_significant_localizes_a_real_effect():
    fit, rng = make_fit(seed=800, n_years=20)
    treatment = make_year_localized(rng, next_year(fit), spike_day_range=(121, 180), spike_size=6.0)
    sig = bin_significant(test_treatment(fit, treatment, rng=rng), alpha=0.05)
    flagged = set(fit.bins[sig].tolist())
    assert {5, 6} <= flagged, f"expected bins 5 and 6 flagged, got {flagged}"
    assert len(flagged - {5, 6}) <= 1


def test_one_sided_test_ignores_the_other_direction():
    fit, rng = make_fit(seed=4, n_years=15, n_bootstrap=400)
    low = make_year(rng, next_year(fit), effect=-3.0)
    assert test_treatment(fit, low, alternative="two-sided", rng=1).p_value <= 0.05
    assert test_treatment(fit, low, alternative="less", rng=1).p_value <= 0.05
    assert test_treatment(fit, low, alternative="greater", rng=1).p_value > 0.5


# ---------------------------------------------------------------------------
# Sequential and Cumulative
# ---------------------------------------------------------------------------

def test_sequential_fwer_control():
    fams = total = 0
    for f in range(8):
        fit, rng = make_fit(seed=300 + f)
        for exp in range(12):
            events = [SequentialEvent(str(k), make_year(rng, 2300 + exp * 10 + k)) for k in range(5)]
            fams += any(r.significant_corrected for r in sequential_test(fit, events, rng=rng))
            total += 1
    assert fams / total <= 0.12, f"family-wise error rate {fams / total:.3f}"


def test_cumulative_flags_then_clears():
    early = late = total = 0
    for f in range(6):
        fit, rng = make_fit(seed=400 + f)
        for i in range(10):
            rows = cumulative_test(fit, [
                CumulativeWindow("early", make_year(rng, 2500 + i * 10, effect=3.0)),
                CumulativeWindow("late", make_year(rng, 2501 + i * 10, effect=0.0)),
            ], rng=rng)
            early += rows[0].significant
            late += rows[1].significant
            total += 1
    assert early / total >= 0.85
    assert late / total <= 0.12


def test_first_recovery_requires_consecutive():
    fit, rng = make_fit(seed=450, n_years=20)
    effects = [3.0, 2.5, 0.0, 0.0, 0.0]
    rows = cumulative_test(fit, [CumulativeWindow(f"w{i}", make_year(rng, 2600 + i, effect=e))
                                 for i, e in enumerate(effects)], rng=rng)
    assert rows[0].significant and rows[1].significant
    rec = first_recovery(rows, consecutive_required=2)
    assert rec.recovered and rec.recovered_at_index == 2
    assert not first_recovery(rows[:1], consecutive_required=2).recovered


# ---------------------------------------------------------------------------
# Sustained-departure test and departure extent
# ---------------------------------------------------------------------------

def test_sustained_run_length_one_matches_standard():
    fit, rng = make_fit(seed=1000, n_years=15)
    for alt in ("two-sided", "greater", "less"):
        r = test_treatment(fit, make_year(rng, next_year(fit), effect=1.5), run_lengths=[1],
                           alternative=alt, rng=rng)
        assert r.sustained[1].p_value == r.p_value
        assert abs(r.sustained[1].test_statistic - r.test_statistic) < 1e-9


def test_sustained_beats_single_bin_for_sustained_effect():
    single = sustained = total = 0
    for f in range(5):
        fit, rng = make_fit(seed=1100 + f, n_years=20, n_bootstrap=500)
        for i in range(15):
            t = make_year_localized(rng, next_year(fit, i), spike_day_range=(91, 210), spike_size=0.5)
            r = test_treatment(fit, t, run_lengths=[4], rng=rng)
            single += r.p_value <= 0.05
            sustained += r.sustained[4].p_value <= 0.05
            total += 1
    assert sustained > single, f"sustained={sustained}/{total}, single={single}/{total}"
    assert sustained / total >= 0.35


def test_sustained_requires_a_consistent_direction():
    """Two high bins then two low bins is not a sustained departure; the old
    statistic (mean of absolute scores) scored it like four high bins."""
    z = np.array([0.0, 0.0, 3.0, 3.0, -3.0, -3.0, 0.0, 0.0])
    two_sided = _window_scores(z, 4, "two-sided")[0]
    assert two_sided[2] == 0.0                     # the up-then-down window cancels
    assert two_sided.max() == 1.5                  # best window is half-in the rise
    assert _window_scores(z, 2, "greater")[0].max() == 3.0
    assert _window_scores(z, 2, "less")[0].max() == 3.0


def test_sustained_null_calibration():
    hits = total = 0
    for f in range(5):
        fit, rng = make_fit(seed=1200 + f, n_years=20, n_bootstrap=500)
        for i in range(12):
            r = test_treatment(fit, make_year(rng, next_year(fit, i)), run_lengths=[4], rng=rng)
            hits += r.sustained[4].p_value <= 0.05
            total += 1
    assert hits / total <= 0.12


def test_departure_recovery_matches_injected_boundaries():
    fit, rng = make_fit(seed=1300, n_years=20, n_bootstrap=500)
    t = make_year_localized(rng, next_year(fit), spike_day_range=(61, 210), spike_size=3.0)  # bins 3-7
    dr = estimate_departure_recovery(test_treatment(fit, t, run_lengths=[3], rng=rng), 3)
    assert dr.recovered
    assert dr.departure_start_bin in (2, 3, 4)
    assert dr.departure_end_bin in (6, 7, 8)
    assert dr.recovered_at_bin == dr.departure_end_bin + 1


def test_departure_recovery_never_recovers():
    fit, rng = make_fit(seed=1301, n_years=20, n_bootstrap=500)
    t = make_year_localized(rng, next_year(fit), spike_day_range=(241, 365), spike_size=3.0)
    dr = estimate_departure_recovery(test_treatment(fit, t, run_lengths=[3], rng=rng), 3)
    assert not dr.recovered and dr.recovered_at_bin is None
    assert dr.departure_end_bin == int(fit.bins[-1])


# ---------------------------------------------------------------------------
# Confidence mode
# ---------------------------------------------------------------------------

def test_confidence_mode_calibration_and_power():
    hits = det = total = 0
    for f in range(8):
        fit, rng = make_fit(seed=1400 + f, n_years=15, n_bootstrap=400)
        for i in range(8):
            base = next_year(fit, 3 * i)
            null = pd.concat([make_year(rng, base + j) for j in range(3)])
            shifted = pd.concat([make_year(rng, base + 30 + j, effect=1.5) for j in range(3)])
            hits += test_treatment(fit, null, mode="confidence", rng=rng).p_value <= 0.05
            det += test_treatment(fit, shifted, mode="confidence", rng=rng).p_value <= 0.05
            total += 1
    assert hits / total <= 0.12
    assert det / total >= 0.8


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
