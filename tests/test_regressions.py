"""
Regression tests: one per defect found in the code audit, each named for
the wrong answer it prevents rather than the function it touches.

Every test in here failed before the audit's fixes. They are separated
from test_engine.py / test_monitoring.py (which validate that the METHOD
behaves statistically) because these validate that the IMPLEMENTATION
does not silently produce a confident wrong number -- a different kind of
failure, and the more dangerous one: a miscalibrated p-value is visible
to a careful reader, a maximally-significant p-value invented out of
missing data is not.

Run with `pytest tests/test_regressions.py -v`, or directly with
`python tests/test_regressions.py` (no pytest required).
"""

import sys
import os
import warnings

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import pandas as pd

from tide_lite import (
    TideConfig, fit_historical, test_treatment, bin_significant, get_envelope,
    holm_bonferroni, sens_slope, estimate_departure_recovery,
    SequentialEvent, sequential_test, CumulativeWindow, cumulative_test,
    first_recovery, MonitoringSeries,
)
from tests.test_engine import make_year, make_fit


def _trending_year(rng, year, slope=0.5, start=2000, n_days=365):
    days = np.arange(1, n_days + 1)
    values = (10.0 + slope * (year - start)
              + 3.0 * np.sin(2 * np.pi * (days - 60) / 365)
              + rng.normal(0, 0.8, n_days))
    return pd.DataFrame({"date": pd.date_range(f"{year}-01-01", periods=n_days, freq="D"),
                         "value": values})


# ---------------------------------------------------------------------------
# 1. A treatment period with missing bins must not be reported as
#    maximally significant in exactly the bins it has no data for.
#
#    Before the fix: obs_bin_devs was NaN there, `null >= nan` is False for
#    every draw, so the bin scored 0 exceedances and came out at
#    1/(n_bootstrap+1) -- the SMALLEST p-value the bootstrap can produce.
#    A treatment year truncated in August was "significant" Sep-Dec.
# ---------------------------------------------------------------------------

def test_missing_treatment_bins_are_not_significant():
    fit, rng = make_fit(seed=7, n_years=12, n_bootstrap=500)
    truncated = make_year(rng, 3000, effect=0.0)
    # months 1-8 = days 1-243; bin = (day-1)//30+1, so bins 1-9 have data
    # and bins 10-13 are absent
    truncated = truncated[truncated["date"].dt.month <= 8]

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = test_treatment(fit, truncated, mode="prediction",
                                 treatment_year_index=12, rng=rng)
    assert any("no data in" in str(w.message) for w in caught), \
        "a treatment period with missing bins should say so"

    assert result.missing_bins.sum() == 4
    assert np.all(np.isnan(result.bin_p_values[result.missing_bins])), \
        "bins with no data must have NaN p-values, not the bootstrap minimum"
    assert np.all(~np.isnan(result.bin_p_values[~result.missing_bins]))
    assert not bin_significant(result, alpha=0.05)[result.missing_bins].any(), \
        "a bin with no data must never be flagged as significant"


def test_sustained_test_ignores_windows_containing_missing_bins():
    """Before the fix a window holding a NaN produced a NaN statistic, which
    np.argmax then selected as the maximum -- and no null draw can exceed
    NaN, so the sustained test returned the smallest possible p-value."""
    fit, rng = make_fit(seed=8, n_years=12, n_bootstrap=500)
    truncated = make_year(rng, 3001, effect=0.0)
    truncated = truncated[truncated["date"].dt.month <= 8]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = test_treatment(fit, truncated, mode="prediction",
                                 treatment_year_index=12, run_lengths=[3], rng=rng)
    sustained = result.sustained[3]
    assert not np.isnan(sustained.test_statistic), "statistic must be a real number"
    assert sustained.p_value > 1.5 / (result.n_reference + 1), \
        "a null treatment period must not land at the bootstrap's minimum p-value"
    # the reported window must lie entirely inside the bins that have data
    assert not result.missing_bins[[list(result.bins).index(b) for b in sustained.window_bins]].any()


def test_treatment_period_too_sparse_is_rejected():
    fit, rng = make_fit(seed=9, n_years=12, n_bootstrap=200)
    sparse = make_year(rng, 3002, effect=0.0)
    sparse = sparse[sparse["date"].dt.month <= 5]           # 7 of 13 bins absent (54%)
    try:
        test_treatment(fit, sparse, mode="prediction", treatment_year_index=12, rng=rng)
        assert False, "expected a ValueError for a treatment period past max_missing_frac"
    except ValueError as e:
        assert "missing" in str(e)


# ---------------------------------------------------------------------------
# 2. Trend slope must be per CALENDAR year, not per position in the list.
#    Records with gaps are normal (a year dropped by the completeness
#    filter, or never sampled); the positional axis silently rescaled the
#    slope, and every residual, MAD and p-value downstream of it.
# ---------------------------------------------------------------------------

def test_sens_slope_uses_the_time_axis_it_is_given():
    values = np.array([0.0, 1.0, 2.0, 13.0, 14.0])
    times = np.array([0.0, 1.0, 2.0, 13.0, 14.0])       # a gap between index 2 and 3
    assert abs(sens_slope(values, times=times) - 1.0) < 1e-9
    assert sens_slope(values) > 1.2, "positional axis should differ, and does"


def test_trend_slope_is_correct_on_a_gapped_record():
    rng = np.random.default_rng(0)
    years = list(range(2000, 2005)) + list(range(2010, 2015))   # five-year hole
    gapped = pd.concat([_trending_year(rng, y, slope=0.5) for y in years],
                       ignore_index=True)
    fit = fit_historical(gapped, "date", "value", TideConfig(n_bootstrap=100))
    assert fit.trend["applied"]
    assert 0.40 < fit.trend["slope_per_year"] < 0.60, (
        f"true drift is 0.5/year; got {fit.trend['slope_per_year']:.3f} "
        f"(the positional axis used to give ~0.87)"
    )


def test_year_index_for_converts_calendar_years():
    rng = np.random.default_rng(1)
    years = list(range(2000, 2005)) + list(range(2010, 2015))
    fit = fit_historical(pd.concat([_trending_year(rng, y) for y in years],
                                    ignore_index=True),
                         "date", "value", TideConfig(n_bootstrap=100))
    assert fit.first_year == 2000
    assert fit.year_index_for(2015) == 15
    assert fit.year_index_for(2000) == 0
    # len(years_used) would have said 10 for a period that is really 15 years out
    assert fit.year_index_for(2015) != len(fit.years_used)


def test_default_treatment_year_index_follows_the_last_historical_year():
    rng = np.random.default_rng(2)
    years = list(range(2000, 2005)) + list(range(2010, 2015))
    fit = fit_historical(pd.concat([_trending_year(rng, y) for y in years],
                                    ignore_index=True),
                         "date", "value", TideConfig(n_bootstrap=100))
    result = test_treatment(fit, _trending_year(rng, 2015), mode="prediction", rng=rng)
    assert result.treatment_year_index == 15, "one calendar year after 2014, not len(years_used)"


# ---------------------------------------------------------------------------
# 3. Misspelled options must raise instead of silently running a different
#    analysis.
# ---------------------------------------------------------------------------

def test_misspelled_options_are_rejected():
    fit, rng = make_fit(seed=10, n_years=10, n_bootstrap=100)
    for bad in ("predicton", "Prediction", "confidance"):
        try:
            test_treatment(fit, make_year(rng, 3100), mode=bad, treatment_year_index=10)
            assert False, f"mode={bad!r} should raise (it used to run 'confidence')"
        except ValueError:
            pass
    for kwargs in ({"detrend_mode": "log"}, {"detrend_mode": "Additive"}, {"agg": "avg"},
                   {"bin_days": 0}, {"n_bootstrap": 0}, {"max_missing_frac": 1.0},
                   {"mk_alpha": 0.0}, {"pool_window_radius": -1}):
        try:
            TideConfig(**kwargs)
            assert False, f"TideConfig({kwargs}) should raise"
        except ValueError:
            pass


def test_holm_bonferroni_rejects_nan_and_bad_alpha():
    for alpha in (0.0, 1.0, -0.1, 1.5):
        try:
            holm_bonferroni(np.array([0.01, 0.2]), alpha=alpha)
            assert False, f"alpha={alpha} should raise"
        except ValueError:
            pass
    try:
        holm_bonferroni(np.array([0.01, np.nan]), alpha=0.05)
        assert False, "NaN p-values should raise, not be sorted into the sequence"
    except ValueError:
        pass


def test_first_recovery_rejects_a_vacuous_rule():
    fit, rng = make_fit(seed=11, n_years=10, n_bootstrap=200)
    rows = cumulative_test(fit, [CumulativeWindow("w0", make_year(rng, 3200, effect=6.0),
                                                   year_index=10)], rng=rng)
    assert rows[0].significant, "sanity: a large effect should be flagged"
    try:
        first_recovery(rows, consecutive_required=0)
        assert False, "consecutive_required=0 used to declare recovery on a significant window"
    except ValueError:
        pass


# ---------------------------------------------------------------------------
# 4. A threshold below the bootstrap's resolution must be reported, not
#    silently returned as "nothing found".
# ---------------------------------------------------------------------------

def test_bin_significant_warns_when_no_bin_could_ever_be_flagged():
    fit, rng = make_fit(seed=12, n_years=10, n_bootstrap=50)   # min p = 1/51 = 0.0196
    result = test_treatment(fit, make_year(rng, 3300, effect=8.0), mode="prediction",
                             treatment_year_index=10, rng=rng)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        flagged = bin_significant(result, alpha=0.05)   # strictest threshold 0.05/13
    assert not flagged.any()
    assert any("No bin can be flagged" in str(w.message) for w in caught), (
        "an all-False result that is arithmetically forced must say so, not look "
        "like evidence the data was normal"
    )


def test_monitoring_rejects_a_horizon_it_could_never_flag():
    fit, _ = make_fit(seed=13, n_years=25, n_bootstrap=300)    # min p = 1/301 = 0.0033
    try:
        MonitoringSeries(fit, alpha_total=0.05, n_years_horizon=100)  # alpha/check 0.0005
        assert False, "a horizon whose per-check alpha is unreachable should raise"
    except ValueError as e:
        assert "cannot flag anything" in str(e)
    for bad in ({"alpha_total": 0.0}, {"alpha_total": 1.0}, {"n_years_horizon": 0}):
        try:
            MonitoringSeries(fit, **bad)
            assert False, f"MonitoringSeries({bad}) should raise"
        except ValueError:
            pass


# ---------------------------------------------------------------------------
# 5. Under a fitted trend the envelope, the effect size and the p-value must
#    all describe the same comparison. They used to disagree completely: a
#    perfectly normal year came back p=0.40 with effect_size=+9.06, sitting
#    above the plotted band in all 13 bins.
# ---------------------------------------------------------------------------

def test_envelope_and_effect_size_agree_with_the_p_value_under_a_trend():
    rng = np.random.default_rng(4)
    hist = pd.concat([_trending_year(rng, y, slope=0.6) for y in range(2000, 2015)],
                     ignore_index=True)
    fit = fit_historical(hist, "date", "value", TideConfig(n_bootstrap=400))
    assert fit.trend["applied"]

    normal_2015 = _trending_year(rng, 2015, slope=0.6)     # an ordinary continuation
    result = test_treatment(fit, normal_2015, mode="prediction",
                             treatment_year_index=15, rng=rng)
    assert result.p_value > 0.05, "sanity: continuing the trend is not an anomaly"

    env = get_envelope(fit, year_index=result.treatment_year_index, rng=rng)
    outside = int(np.sum((result.treatment_curve > env["upper"]) |
                         (result.treatment_curve < env["lower"])))
    assert outside <= 3, (
        f"a non-significant year sat outside the band in {outside}/13 bins -- the "
        f"envelope is not on the treatment year's trend basis"
    )
    assert abs(result.effect_size_trend_adjusted) < 1.0, \
        "the effect size the p-value is based on should be near zero"
    # effect_size is measured against the record's own typical level (mid-record),
    # so it is bounded by the drift over half the record -- not, as before, by the
    # drift since the very first historical year.
    assert abs(result.effect_size) < 0.6 * len(fit.years_used)


def test_envelope_is_unshifted_when_no_trend_was_applied():
    fit, rng = make_fit(seed=14, n_years=10, n_bootstrap=200)
    assert not fit.trend["applied"]
    a = get_envelope(fit, year_index=None, rng=np.random.default_rng(5))
    b = get_envelope(fit, year_index=99, rng=np.random.default_rng(5))
    assert np.allclose(a["median_curve"], b["median_curve"])
    assert a["trend_shift"] == 0.0 and b["trend_shift"] == 0.0


# ---------------------------------------------------------------------------
# 6. Reproducibility: a reported p-value must be reproducible from a seed.
# ---------------------------------------------------------------------------

def test_controllers_are_reproducible_from_a_seed():
    fit, rng = make_fit(seed=15, n_years=10, n_bootstrap=200)
    # A weak effect on purpose: a strong one pins the p-value at the
    # bootstrap floor 1/(n+1) for every seed, which cannot show seed use.
    year = make_year(rng, 3400, effect=0.3)

    events = [SequentialEvent("a", year, year_index=10)]
    assert (sequential_test(fit, events, rng=42)[0].result.p_value ==
            sequential_test(fit, events, rng=42)[0].result.p_value)

    windows = [CumulativeWindow("a", year, year_index=10)]
    assert (cumulative_test(fit, windows, rng=7)[0].result.p_value ==
            cumulative_test(fit, windows, rng=7)[0].result.p_value)

    # and the seed is actually used: different seeds give different draws
    spread = {sequential_test(fit, events, rng=s)[0].result.p_value for s in range(1, 6)}
    assert len(spread) > 1, f"every seed gave the same p-value: {spread}"


# ---------------------------------------------------------------------------
# 7. The between-year diagnostic must actually detect whole-year level
#    shifts -- it is what tells a user which calibration regime they are in.
# ---------------------------------------------------------------------------

def _year_with_level_shift(rng, year, year_sd):
    n = 365
    days = np.arange(1, n + 1)
    e = rng.normal(0, 0.6, n)
    noise = np.zeros(n)
    for i in range(1, n):
        noise[i] = 0.7 * noise[i - 1] + e[i]
    values = 10.0 + 3.0 * np.sin(2 * np.pi * (days - 60) / 365) + noise + rng.normal(0, year_sd)
    return pd.DataFrame({"date": pd.date_range(f"{year}-01-01", periods=n, freq="D"),
                         "value": values})


def test_between_year_var_frac_separates_the_calibration_regimes():
    def frac(year_sd):
        rng = np.random.default_rng(900)
        hist = pd.concat([_year_with_level_shift(rng, 2000 + i, year_sd) for i in range(15)],
                         ignore_index=True)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return fit_historical(hist, "date", "value",
                                   TideConfig(n_bootstrap=100)).between_year_var_frac

    clean, shifted = frac(0.0), frac(1.2)
    assert clean < 0.3, f"pure within-year noise should read low, got {clean:.2f}"
    assert shifted > 0.6, f"strong whole-year shifts should read high, got {shifted:.2f}"


def test_fit_warns_on_high_between_year_variance_with_a_short_record():
    rng = np.random.default_rng(901)
    hist = pd.concat([_year_with_level_shift(rng, 2000 + i, 1.2) for i in range(10)],
                     ignore_index=True)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        fit_historical(hist, "date", "value", TideConfig(n_bootstrap=100))
    assert any("whole-year level shifts" in str(w.message) for w in caught), (
        "a short record dominated by whole-year shifts is the anti-conservative "
        "regime and should say so"
    )


# ---------------------------------------------------------------------------
# 8. estimate_departure_recovery must not claim recovery it cannot see.
# ---------------------------------------------------------------------------

def test_departure_recovery_flags_recovery_into_missing_data():
    fit, rng = make_fit(seed=16, n_years=20, n_bootstrap=400)
    treatment = make_year(rng, 3500, effect=0.0)
    spike = (treatment["date"].dt.dayofyear >= 121) & (treatment["date"].dt.dayofyear <= 240)
    treatment.loc[spike, "value"] += 4.0
    treatment = treatment[treatment["date"].dt.month <= 9]   # bins 10-13 absent

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = test_treatment(fit, treatment, mode="prediction",
                                 treatment_year_index=20, run_lengths=[3], rng=rng)
        dr = estimate_departure_recovery(result, run_length=3)
    if dr.recovered and np.isnan(result.bin_p_values[
            list(result.bins).index(dr.recovered_at_bin)]):
        assert any("no treatment data" in str(w.message) for w in caught), (
            "'recovered' into a bin with no data must be called out"
        )


if __name__ == "__main__":
    _here = sys.modules[__name__]
    tests = [v for k, v in sorted(vars(_here).items())
             if k.startswith("test_") and callable(v)
             and getattr(v, "__module__", None) == __name__]
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
