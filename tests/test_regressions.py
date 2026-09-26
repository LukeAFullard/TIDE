"""
Regression tests: one per implementation defect found in audits, each named
for the wrong answer it prevents. Every one of these failed on the code as
it was before the fix. Separate from test_engine.py (does the METHOD behave
statistically?) because these check the IMPLEMENTATION never hands back a
confident wrong number. Run with `pytest tests/test_regressions.py -v`.
"""

import os
import sys
import warnings

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import pandas as pd
import pytest

import tide_lite.engine as E
from tide_lite import (
    TideConfig, fit_historical, test_treatment, bin_significant, get_envelope,
    holm_bonferroni, sens_slope, estimate_departure_recovery, significance_band,
    SequentialEvent, sequential_test, CumulativeWindow, cumulative_test,
    first_recovery, MonitoringSeries, sensitivity_grid, plot_sequential, summarize,
)
from tests.test_engine import make_year, make_fit, next_year

warnings.simplefilter("ignore", UserWarning)


def _trending_year(rng, year, slope=0.5, start=2000, n_days=365, sd=0.8):
    days = np.arange(1, n_days + 1)
    values = (10.0 + slope * (year - start) + 3.0 * np.sin(2 * np.pi * (days - 60) / 365)
              + rng.normal(0, sd, n_days))
    return pd.DataFrame({"date": pd.date_range(f"{year}-01-01", periods=n_days, freq="D"),
                         "value": values})


def _trending_fit(years, slope=0.5, seed=0, n_bootstrap=200):
    rng = np.random.default_rng(seed)
    hist = pd.concat([_trending_year(rng, y, slope) for y in years], ignore_index=True)
    return fit_historical(hist, "date", "value", TideConfig(n_bootstrap=n_bootstrap)), rng


# ---------------------------------------------------------------------------
# Missing treatment data
# ---------------------------------------------------------------------------

def test_missing_treatment_bins_are_not_significant():
    """Bins with no data once scored at the SMALLEST possible p-value
    (NaN >= null is False for every draw)."""
    fit, rng = make_fit(seed=7, n_years=12, n_bootstrap=500)
    truncated = make_year(rng, next_year(fit))
    truncated = truncated[truncated["date"].dt.month <= 8]   # days <= 243 of 30-day bins:
    with warnings.catch_warnings(record=True) as caught:     # bins 1-8 full, bin 9 has 3 days
        warnings.simplefilter("always")
        r = test_treatment(fit, truncated, rng=rng)
    assert any("No treatment data in bins" in str(w.message) for w in caught)
    assert r.missing_bins.sum() == 4                         # 9 (too thin) + 10-12 (empty)
    assert np.all(np.isnan(r.bin_p_values[r.missing_bins]))
    assert np.all(np.isnan(r.bin_p_adjusted[r.missing_bins]))
    assert not bin_significant(r)[r.missing_bins].any()


def test_null_uses_only_the_bins_the_treatment_has():
    """The warning said the p-value was based on the bins with data, but the
    null's maximum was still taken over ALL bins, a stricter bar than the
    observed statistic faced (a silently wrong, too-large p-value)."""
    fit, rng = make_fit(seed=17, n_years=12, n_bootstrap=400)
    t = make_year(rng, next_year(fit), effect=1.0)
    t = t[t["date"].dt.month <= 8]
    r = test_treatment(fit, t, rng=123)
    x0 = [float(next_year(fit) - fit.first_year)]
    null = E._null_scores(fit, 400, np.random.default_rng(123), x0)
    expected = np.abs(null[:, ~r.missing_bins]).max(axis=1)
    assert np.allclose(r.null_max_stats, expected)


def test_sustained_test_ignores_windows_containing_missing_bins():
    fit, rng = make_fit(seed=8, n_years=12, n_bootstrap=500)
    t = make_year(rng, next_year(fit))
    t = t[t["date"].dt.month <= 8]
    r = test_treatment(fit, t, run_lengths=[3], rng=rng)
    s = r.sustained[3]
    assert np.isfinite(s.test_statistic)
    assert s.p_value > 1.5 / (r.n_reference + 1)
    idx = [list(r.bins).index(b) for b in s.window_bins]
    assert not r.missing_bins[idx].any()


def test_treatment_period_too_sparse_is_rejected():
    fit, rng = make_fit(seed=9, n_years=12, n_bootstrap=200)
    sparse = make_year(rng, next_year(fit))
    sparse = sparse[sparse["date"].dt.month <= 4]            # 8 of 12 bins absent
    with pytest.raises(ValueError, match="missing"):
        test_treatment(fit, sparse, rng=rng)


# ---------------------------------------------------------------------------
# Historical data screening
# ---------------------------------------------------------------------------

def test_bins_with_too_little_history_are_dropped():
    """A bin with data in only a couple of historical years gave a median
    and MAD from 1-2 values (or NaN), which could score a treatment value
    there at the minimum p-value."""
    rng = np.random.default_rng(20)
    years = []
    for i in range(12):
        y = make_year(rng, 2000 + i)
        if i >= 2:
            y = y[y["date"].dt.month != 12]                    # December sampled in 2 years only
        years.append(y)
    fit = fit_historical(pd.concat(years), "date", "value", TideConfig(n_bootstrap=200))
    assert 12 not in fit.bins
    assert any(b == 12 for b, _ in fit.bins_dropped)
    r = test_treatment(fit, make_year(rng, 2012, effect=0.0), rng=1)
    assert len(r.bins) == 11


def test_non_numeric_values_raise():
    """'<0.5' style censored results were silently dropped, biasing the record up."""
    rng = np.random.default_rng(21)
    df = make_year(rng, 2000)
    df["value"] = df["value"].astype(object)
    df.loc[5, "value"] = "<0.5"
    with pytest.raises(ValueError, match="not numbers"):
        fit_historical(df, "date", "value")


def test_partial_bins_are_not_compared_with_full_ones():
    """A treatment period ending on 30 June put ONE day into the 30-day bin
    starting 30 June; a one-day 'median' is far noisier than the 30-day
    medians it was compared with, and normal half-years were flagged ~30% of
    the time at a nominal 5%."""
    fit, rng = make_fit(seed=32, n_years=15, n_bootstrap=300)
    half = make_year(rng, next_year(fit))
    half = half[half["date"].dt.month <= 6]                  # 30 June = day 181 = bin 7
    r = test_treatment(fit, half, rng=1)
    assert r.missing_bins[6] and not r.missing_bins[:6].any()


def test_calendar_month_bins_match_monthly_sampling():
    """Monthly samples on the 1st: 30-day bins leave one bin empty and put two
    samples in another; calendar-month bins (the default) do not."""
    dates = pd.date_range("2001-01-01", periods=12, freq="MS")
    df = pd.DataFrame({"date": dates, "value": np.arange(12.0)})
    by30 = E._bin_and_pivot(df, "date", "value", 30, "median")
    bymonth = E._bin_and_pivot(df, "date", "value", "month", "median")
    assert by30.isna().sum().sum() >= 1
    assert bymonth.notna().all().all()


# ---------------------------------------------------------------------------
# Trend handling
# ---------------------------------------------------------------------------

def test_sens_slope_uses_the_time_axis_it_is_given():
    values = np.array([0.0, 1.0, 2.0, 13.0, 14.0])
    assert abs(sens_slope(values, times=values) - 1.0) < 1e-9
    assert sens_slope(values) > 1.2


def test_trend_slope_is_correct_on_a_gapped_record():
    fit, _ = _trending_fit(list(range(2000, 2005)) + list(range(2010, 2015)))
    assert fit.trend["applied"]
    assert 0.40 < fit.trend["slope_per_year"] < 0.60


def test_missing_months_do_not_create_a_trend():
    """The trend test used each year's RAW median, which moves with which
    months a year has; early years missing summer looked like a step up."""
    rng = np.random.default_rng(22)
    years = []
    for i in range(14):
        y = make_year(rng, 2000 + i, ar_noise_sd=0.2)
        if i < 7:
            y = y[~y["date"].dt.month.isin([6, 7, 8, 9])]      # 4 of 12 months missing
        years.append(y)
    fit = fit_historical(pd.concat(years), "date", "value", TideConfig(n_bootstrap=100))
    assert not fit.trend["applied"], f"spurious trend {fit.trend}"


def test_default_year_index_comes_from_the_data_dates():
    """The default used to be 'one year after the record' whatever the data's
    date, so a 2020 period on a 2000-2014 trend was projected to 2015."""
    fit, rng = _trending_fit(range(2000, 2015))
    r = test_treatment(fit, _trending_year(rng, 2020), rng=1)
    assert r.treatment_year_index == 20
    assert r.p_value > 0.05
    assert abs(r.effect_size_trend_adjusted) < 1.0


def test_year_index_for_converts_calendar_years():
    fit, _ = _trending_fit(list(range(2000, 2005)) + list(range(2010, 2015)), seed=1)
    assert fit.year_index_for(2015) == 15 and fit.year_index_for(2000) == 0


def test_treatment_year_inside_the_history_is_rejected():
    fit, rng = make_fit(seed=23, n_years=10)
    with pytest.raises(ValueError, match="also in the historical record"):
        test_treatment(fit, make_year(rng, fit.years_used[3]), rng=1)


def test_envelope_and_effect_size_agree_with_the_p_value_under_a_trend():
    fit, rng = _trending_fit(range(2000, 2015), slope=0.6, seed=4, n_bootstrap=400)
    assert fit.trend["applied"]
    r = test_treatment(fit, _trending_year(rng, 2015, slope=0.6), rng=rng)
    assert r.p_value > 0.05
    band = significance_band(fit, r, 0.05)
    assert not np.any((r.treatment_curve > band["upper"]) | (r.treatment_curve < band["lower"]))
    env = get_envelope(fit, year_index=r.treatment_year_index, rng=rng)
    outside = np.sum((r.treatment_curve > env["upper"]) | (r.treatment_curve < env["lower"]))
    assert outside <= 4
    assert abs(r.effect_size_trend_adjusted) < 1.0
    assert abs(r.effect_size) < 0.6 * len(fit.years_used)


def test_envelope_is_unshifted_when_no_trend_was_applied():
    fit, _ = make_fit(seed=14, n_years=10, n_bootstrap=200)
    assert not fit.trend["applied"]
    a = get_envelope(fit, year_index=None, rng=np.random.default_rng(5))
    b = get_envelope(fit, year_index=99, rng=np.random.default_rng(5))
    assert np.allclose(a["median_curve"], b["median_curve"]) and b["trend_shift"] == 0.0
    with pytest.raises(ValueError):
        get_envelope(fit, method="bootstrapp")


# ---------------------------------------------------------------------------
# Confidence mode (several years averaged)
# ---------------------------------------------------------------------------

def test_confidence_mode_averages_in_log_space():
    """Years were averaged in raw units, then logged: log(mean) > mean(log),
    so a set of perfectly typical years looked systematically high."""
    rng = np.random.default_rng(24)
    hist = pd.concat([make_year(rng, 2000 + i) for i in range(12)], ignore_index=True)
    fit = fit_historical(hist, "date", "value", TideConfig(detrend_mode="log_additive", n_bootstrap=200))
    typical = np.exp(fit.median_curve)
    rows = []
    for year, factor in ((2013, np.e ** 0.5), (2014, np.e ** -0.5)):
        dates = pd.date_range(f"{year}-01-01", f"{year}-12-31", freq="D")
        rows.append(pd.DataFrame({"date": dates, "value": typical[dates.month - 1] * factor}))
    r = test_treatment(fit, pd.concat(rows), mode="confidence", rng=1)
    assert np.allclose(r.treatment_curve, typical, rtol=1e-9)
    assert np.allclose(r.bin_scores, 0.0, atol=1e-9)


def test_confidence_mode_projects_the_trend_to_each_year():
    fit, rng = _trending_fit(range(2000, 2015), seed=25)
    both = pd.concat([_trending_year(rng, 2015), _trending_year(rng, 2016)])
    r = test_treatment(fit, both, mode="confidence", rng=1)
    assert r.treatment_year_index == 15.5
    assert abs(r.effect_size_trend_adjusted) < 0.6
    with pytest.raises(ValueError):
        test_treatment(fit, both, mode="confidence", treatment_year_index=15)


# ---------------------------------------------------------------------------
# Options and thresholds
# ---------------------------------------------------------------------------

def test_misspelled_options_are_rejected():
    fit, rng = make_fit(seed=10, n_years=10, n_bootstrap=100)
    year = make_year(rng, next_year(fit))
    for bad in ("predicton", "Prediction", "confidance"):
        with pytest.raises(ValueError):
            test_treatment(fit, year, mode=bad)
    with pytest.raises(ValueError):
        test_treatment(fit, year, alternative="two_sided")
    for kwargs in ({"detrend_mode": "log"}, {"agg": "avg"}, {"bin_days": 0}, {"bin_days": "week"},
                   {"n_bootstrap": 0}, {"max_missing_frac": 1.0}, {"mk_alpha": 0.0},
                   {"pool_window_radius": -1}):
        with pytest.raises(ValueError):
            TideConfig(**kwargs)


def test_holm_bonferroni_rejects_nan_and_bad_alpha():
    for alpha in (0.0, 1.0, -0.1, 1.5):
        with pytest.raises(ValueError):
            holm_bonferroni(np.array([0.01, 0.2]), alpha=alpha)
    with pytest.raises(ValueError):
        holm_bonferroni(np.array([0.01, np.nan]))


def test_one_significance_rule_everywhere():
    """p exactly equal to alpha is significant in every code path (it was
    `<` in some and `<=` in others)."""
    assert E._le(100 / 2000, 0.05)
    assert list(holm_bonferroni(np.array([0.05]), 0.05)) == [True]


def test_first_recovery_rejects_a_vacuous_rule():
    fit, rng = make_fit(seed=11, n_years=10, n_bootstrap=200)
    rows = cumulative_test(fit, [CumulativeWindow("w0", make_year(rng, next_year(fit), effect=6.0))], rng=rng)
    assert rows[0].significant
    with pytest.raises(ValueError):
        first_recovery(rows, consecutive_required=0)


def test_bin_significant_warns_when_nothing_could_be_flagged():
    fit, rng = make_fit(seed=12, n_years=10, n_bootstrap=15)   # min p = 1/16 > 0.05
    r = test_treatment(fit, make_year(rng, next_year(fit), effect=8.0), rng=rng)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        assert not bin_significant(r, alpha=0.05).any()
    assert any("Nothing can be flagged" in str(w.message) for w in caught)


def test_monitoring_rejects_a_horizon_it_could_never_flag():
    fit, _ = make_fit(seed=13, n_years=25, n_bootstrap=300)
    with pytest.raises(ValueError, match="cannot flag anything"):
        MonitoringSeries(fit, alpha_total=0.05, n_years_horizon=100)
    for bad in ({"alpha_total": 0.0}, {"alpha_total": 1.0}, {"n_years_horizon": 0}):
        with pytest.raises(ValueError):
            MonitoringSeries(fit, **bad)


# ---------------------------------------------------------------------------
# Controllers, sensitivity grid, monitoring state
# ---------------------------------------------------------------------------

def test_sequential_rejects_overlapping_events():
    fit, rng = make_fit(seed=26, n_years=10, n_bootstrap=100)
    y = make_year(rng, next_year(fit))
    with pytest.raises(ValueError, match="overlap"):
        sequential_test(fit, [SequentialEvent("a", y), SequentialEvent("b", y.copy())])


def test_controllers_are_reproducible_from_a_seed():
    fit, rng = make_fit(seed=15, n_years=10, n_bootstrap=200)
    year = make_year(rng, next_year(fit), effect=0.3)
    ev = [SequentialEvent("a", year)]
    assert sequential_test(fit, ev, rng=42)[0].result.p_value == sequential_test(fit, ev, rng=42)[0].result.p_value
    w = [CumulativeWindow("a", year)]
    assert cumulative_test(fit, w, rng=7)[0].result.p_value == cumulative_test(fit, w, rng=7)[0].result.p_value
    assert len({sequential_test(fit, ev, rng=s)[0].result.p_value for s in range(1, 6)}) > 1


def test_sensitivity_grid_uses_common_random_numbers_and_survives_errors():
    """Each variant used to draw fresh random numbers, so rows differed by
    Monte Carlo noise even when the setting made no difference; and one
    impossible variant crashed the whole grid."""
    rng = np.random.default_rng(27)
    hist = pd.concat([make_year(rng, 2000 + i) for i in range(10)], ignore_index=True)
    hist.loc[hist["date"] < "2000-02-01", "value"] = 0.0      # a month of zeros: no logs
    treat = make_year(rng, 2010, effect=0.5)
    grid = sensitivity_grid(hist, treat, "date", "value", TideConfig(n_bootstrap=300),
                            {"agg": ["median"], "detrend_mode": ["log_additive"]}, rng=5)
    assert grid.loc[0, "p_value"] == grid.loc[1, "p_value"]   # same setting, same seed
    assert "strictly positive" in grid.loc[2, "error"]


def test_monitoring_refuses_a_different_fit_on_reload(tmp_path):
    fit, rng = make_fit(seed=28, n_years=20, n_bootstrap=2000)
    series = MonitoringSeries(fit, alpha_total=0.05, n_years_horizon=5)
    series.check(make_year(rng, next_year(fit)), rng=1)
    path = tmp_path / "state.json"
    series.save(str(path))
    assert len(MonitoringSeries.load(fit, str(path)).checks) == 1
    other, _ = make_fit(seed=29, n_years=20, n_bootstrap=2000)
    with pytest.raises(ValueError, match="different historical fit"):
        MonitoringSeries.load(other, str(path))
    with pytest.raises(ValueError, match="already been checked"):
        series.check(make_year(rng, next_year(fit)), rng=2)


# ---------------------------------------------------------------------------
# Diagnostics, departure extent, plotting, report
# ---------------------------------------------------------------------------

def test_between_year_var_frac_separates_the_regimes():
    def frac(year_sd):
        rng = np.random.default_rng(900)
        hist = pd.concat([make_year(rng, 2000 + i, year_sd=year_sd) for i in range(15)], ignore_index=True)
        return fit_historical(hist, "date", "value", TideConfig(n_bootstrap=100)).between_year_var_frac
    assert frac(0.0) < 0.3
    assert frac(1.2) > 0.6


def test_departure_recovery_flags_recovery_into_missing_data():
    fit, rng = make_fit(seed=16, n_years=20, n_bootstrap=400)
    t = make_year(rng, next_year(fit))
    day = t["date"].dt.dayofyear
    t.loc[(day >= 121) & (day <= 240), "value"] += 4.0
    t = t[t["date"].dt.month <= 9]
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        r = test_treatment(fit, t, run_lengths=[3], rng=rng)
        dr = estimate_departure_recovery(r, run_length=3)
    if dr.recovered and np.isnan(r.bin_p_values[list(r.bins).index(dr.recovered_at_bin)]):
        assert any("no treatment data" in str(w.message) for w in caught)


def test_plot_sequential_handles_a_missing_last_bin():
    import matplotlib
    matplotlib.use("Agg")
    fit, rng = make_fit(seed=30, n_years=10, n_bootstrap=100)
    t = make_year(rng, next_year(fit))
    t = t[t["date"].dt.month <= 10]
    ax = plot_sequential(fit, sequential_test(fit, [SequentialEvent("a", t)], rng=1))
    assert ax is not None


def test_summary_states_the_answer_and_how_to_reproduce_it():
    fit, rng = make_fit(seed=31, n_years=12, n_bootstrap=200)
    r = test_treatment(fit, make_year(rng, next_year(fit), effect=4.0), rng=42)
    text = summarize(fit, r)
    assert "Answer: YES" in text and "rng=42" in text and fit.fingerprint() in text


def test_summary_never_reports_the_floor_p_value_as_exact():
    """At the smallest attainable p-value the true p-value may be smaller;
    'p = 0.004975' read as exact, and 'NO' alone was easily read as 'complied'."""
    fit, rng = make_fit(seed=31, n_years=12, n_bootstrap=200)
    r = test_treatment(fit, make_year(rng, next_year(fit), effect=8.0), rng=42)
    assert r.p_value == 1 / 201
    assert "p <= 0.005" in summarize(fit, r)
    r0 = test_treatment(fit, make_year(rng, next_year(fit, 1)), rng=42)
    if r0.p_value > 0.05:
        assert "no departure from the site's normal range was detected" in summarize(fit, r0)


# ---------------------------------------------------------------------------
# Null distribution
# ---------------------------------------------------------------------------

def test_every_bin_is_represented_once_in_each_synthetic_year():
    """With blocks borrowed from +/-1 neighbouring bin (the old default), the
    first and last bins of the year appeared only ~0.6 times per synthetic
    year, so their own variability was under-represented in the null."""
    fit, _ = make_fit(seed=5, n_years=12, n_bootstrap=100, bin_days="month")
    nb = len(fit.bins)
    fit.within_year = np.tile(np.arange(nb, dtype=float), (fit.within_year.shape[0], 1))
    out = E._stitch_within_year(fit, 3000, np.random.default_rng(0))
    assert np.array_equal(out, np.tile(np.arange(nb, dtype=float), (3000, 1)))


# ---------------------------------------------------------------------------
# Sampling design
# ---------------------------------------------------------------------------

def test_sustained_window_never_spans_a_dropped_bin():
    """With June never sampled, a 3-bin window could be May-Jul-Aug and be
    reported as three consecutive months."""
    rng = np.random.default_rng(7)
    hist = pd.concat([make_year(rng, 2000 + i) for i in range(12)], ignore_index=True)
    hist = hist[hist["date"].dt.month != 6]
    fit = fit_historical(hist, "date", "value", TideConfig(n_bootstrap=200))
    t = make_year(rng, 2013)
    t.loc[t["date"].dt.month.isin([5, 7, 8]), "value"] += 6.0   # high May, Jul, Aug
    t = t[t["date"].dt.month != 6]
    r = test_treatment(fit, t, run_lengths=[3], rng=1)
    w = r.sustained[3].window_bins
    assert w[-1] - w[0] == 2, f"window {w} is not three consecutive months"


def test_projecting_a_trend_far_beyond_the_record_warns():
    """A removed trend projected 10 years past the record was used silently;
    when the trend was not real, false alarms reached 16-29% at a nominal 5%."""
    fit, rng = _trending_fit(range(2000, 2015))
    assert fit.trend["applied"]
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        test_treatment(fit, _trending_year(rng, 2016), rng=1)
        assert not any("projected" in str(x.message) for x in w)
        test_treatment(fit, _trending_year(rng, 2024), rng=1)
    assert any("projected 10 years beyond" in str(x.message) for x in w)


def test_denser_treatment_sampling_than_history_warns():
    """Daily data tested against a history of monthly grab samples was
    compared silently; its monthly medians vary far less than one grab."""
    rng = np.random.default_rng(8)
    hist = pd.concat([make_year(rng, 2000 + i) for i in range(12)], ignore_index=True)
    grab = hist.groupby([hist["date"].dt.year, hist["date"].dt.month]).head(1)
    fit = fit_historical(grab, "date", "value", TideConfig(n_bootstrap=100))
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        test_treatment(fit, make_year(rng, 2013), rng=1)
    assert any("more than twice the usual number" in str(x.message) for x in w)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
