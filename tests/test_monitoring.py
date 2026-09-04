"""
Tests for MonitoringSeries (the fixed-horizon, flat-split monitoring
mode). Same pooling rationale as tests/test_engine.py -- see that file's
docstring for why: per-fit calibration varies enough with ~10 historical
years that a single-fit stochastic test is flaky by construction.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import pandas as pd

from tide_lite import TideConfig, fit_historical, MonitoringSeries
from tests.test_engine import make_year, make_historical, make_fit


# ---------------------------------------------------------------------------
# 1. Deterministic mechanics: alpha split, check numbering, flagged logic
# ---------------------------------------------------------------------------

def test_alpha_split_and_numbering():
    fit, rng = make_fit(seed=500)
    series = MonitoringSeries(fit, alpha_total=0.05, n_years_horizon=10)
    assert series.alpha_per_check == 0.005
    assert series.remaining == 10

    r1 = series.check(make_year(rng, 3000), label="y1", year_index=len(fit.years_used), rng=rng)
    assert r1.check_number == 1
    assert r1.alpha_used == 0.005
    assert r1.flagged == (r1.p_value < 0.005)
    assert series.remaining == 9

    r2 = series.check(make_year(rng, 3001), label="y2", year_index=len(fit.years_used) + 1, rng=rng)
    assert r2.check_number == 2
    assert series.remaining == 8
    assert len(series.checks) == 2


# ---------------------------------------------------------------------------
# 2. Horizon exhaustion raises; renew() gives a fresh, empty budget
# ---------------------------------------------------------------------------

def test_horizon_exhaustion_and_renew():
    fit, rng = make_fit(seed=501)
    series = MonitoringSeries(fit, alpha_total=0.05, n_years_horizon=2)
    series.check(make_year(rng, 3100), year_index=len(fit.years_used), rng=rng)
    series.check(make_year(rng, 3101), year_index=len(fit.years_used) + 1, rng=rng)
    assert series.remaining == 0
    try:
        series.check(make_year(rng, 3102), year_index=len(fit.years_used) + 2, rng=rng)
        assert False, "expected a ValueError when the horizon is exhausted"
    except ValueError as e:
        assert "used up" in str(e)

    fresh = series.renew()
    assert fresh.remaining == 2
    assert len(fresh.checks) == 0
    assert len(series.checks) == 2  # original untouched


# ---------------------------------------------------------------------------
# 3. year_index required when a trend is active
# ---------------------------------------------------------------------------

def test_year_index_required_with_trend():
    rng = np.random.default_rng(502)
    trend_dfs = [pd.DataFrame({"date": pd.date_range(f"{y}-01-01", periods=365, freq="D"),
                                "value": 10 + 0.5 * (y - 2010) + rng.normal(0, 1, 365)})
                 for y in range(2010, 2020)]
    fit = fit_historical(pd.concat(trend_dfs), "date", "value", TideConfig())
    assert fit.trend["applied"]
    series = MonitoringSeries(fit)
    try:
        series.check(make_year(rng, 3200), rng=rng)  # no year_index
        assert False, "expected a ValueError"
    except ValueError as e:
        assert "year_index" in str(e)


# ---------------------------------------------------------------------------
# 4. Save/load round-trip
# ---------------------------------------------------------------------------

def test_save_load_round_trip(tmp_path_str="/tmp/tide_lite_monitoring_test.json"):
    fit, rng = make_fit(seed=503)
    series = MonitoringSeries(fit, alpha_total=0.05, n_years_horizon=5)
    series.check(make_year(rng, 3300), label="first", year_index=len(fit.years_used), rng=rng)
    series.check(make_year(rng, 3301), label="second", year_index=len(fit.years_used) + 1, rng=rng)
    series.save(tmp_path_str)

    reloaded = MonitoringSeries.load(fit, tmp_path_str)
    assert reloaded.alpha_total == series.alpha_total
    assert reloaded.n_years_horizon == series.n_years_horizon
    assert len(reloaded.checks) == 2
    assert reloaded.checks[0].label == "first"
    assert reloaded.checks[1].p_value == series.checks[1].p_value
    os.remove(tmp_path_str)


# ---------------------------------------------------------------------------
# 5. The actual guarantee: "ever falsely flag across the horizon" stays
#    reasonably near alpha_total when enough historical data is available.
#    Uses fresh, independent synthetic years (the same trusted method as
#    test_engine.py's calibration tests) -- NOT leave-one-out against the
#    real historical years themselves. A leave-one-out version was tried
#    during development and rejected: repeatedly holding out from a small
#    pool of real years (with replacement) produces heavy repetition that
#    understates the true risk, worst exactly where the risk is highest.
#    See the IMPORTANT note in monitoring.py's module docstring.
#
#    This mode needs more historical years than the others for its
#    guarantee to hold -- tested at n_years=20 here, not the n_years=10
#    used elsewhere in this project, for exactly that documented reason.
# ---------------------------------------------------------------------------

def test_horizon_wide_false_alarm_rate():
    n_fits = 5
    n_careers_per_fit = 12
    n_years_horizon = 10
    alpha_total = 0.05
    ever_flagged = 0
    total_careers = 0
    for fit_seed in range(n_fits):
        rng = np.random.default_rng(600 + fit_seed)
        historical = make_historical(rng, n_years=20)  # see comment above
        # n_bootstrap must comfortably resolve alpha_per_check (0.005 here).
        # At 300 draws the finest attainable p-value is 1/301 = 0.0033, only
        # 1.5x below the threshold being measured, so the "ever flagged" rate
        # would be governed by bootstrap granularity as much as by the data.
        # MonitoringSeries now warns about exactly this.
        config = TideConfig(bin_days=30, detrend_mode="additive", n_bootstrap=2000)
        fit = fit_historical(historical, "date", "value", config)
        for career in range(n_careers_per_fit):
            series = MonitoringSeries(fit, alpha_total=alpha_total, n_years_horizon=n_years_horizon)
            career_flagged = False
            for i in range(n_years_horizon):
                y = make_year(rng, 8000 + career * 20 + i, effect=0.0)
                r = series.check(y, year_index=len(fit.years_used), rng=rng)
                career_flagged = career_flagged or r.flagged
            ever_flagged += career_flagged
            total_careers += 1
    rate = ever_flagged / total_careers
    # generous band above alpha_total -- same reasoning as test_engine.py's
    # pooled calibration tests: small samples need room, not a tight bound
    assert rate <= 0.20, f"ever-flagged-falsely rate {rate:.3f} exceeds the {alpha_total} budget by more than expected noise allows, even at n_years=20"


if __name__ == "__main__":
    import time
    tests = [
        test_alpha_split_and_numbering,
        test_horizon_exhaustion_and_renew,
        test_year_index_required_with_trend,
        test_save_load_round_trip,
        test_horizon_wide_false_alarm_rate,
    ]
    failures = 0
    for t in tests:
        t0 = time.time()
        try:
            t()
            print(f"PASS  {t.__name__}  ({time.time()-t0:.1f}s)")
        except AssertionError as e:
            failures += 1
            print(f"FAIL  {t.__name__}: {e}")
    print(f"\n{len(tests) - failures}/{len(tests)} passed")
    sys.exit(1 if failures else 0)
