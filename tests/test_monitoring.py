"""
MonitoringSeries: mechanics (threshold split, numbering, horizon, state)
and the horizon-wide false-alarm rate, pooled across simulated histories
for the reason given at the top of test_engine.py.
"""

import os
import sys
import warnings

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import pandas as pd
import pytest

from tide_lite import TideConfig, fit_historical, MonitoringSeries
from tests.test_engine import make_year, make_historical, make_fit, next_year

warnings.simplefilter("ignore", UserWarning)


def test_alpha_split_and_numbering():
    fit, rng = make_fit(seed=500)
    series = MonitoringSeries(fit, alpha_total=0.05, n_years_horizon=10)
    assert series.alpha_per_check == 0.005 and series.remaining == 10
    r1 = series.check(make_year(rng, next_year(fit)), label="y1", rng=rng)
    assert r1.check_number == 1 and r1.alpha_used == 0.005
    assert r1.flagged == (r1.p_value <= 0.005)
    assert r1.calendar_year == next_year(fit)
    r2 = series.check(make_year(rng, next_year(fit, 1)), rng=rng)
    assert r2.check_number == 2 and r2.label == str(next_year(fit, 1))
    assert series.remaining == 8


def test_horizon_exhaustion_and_renew():
    fit, rng = make_fit(seed=501)
    series = MonitoringSeries(fit, alpha_total=0.05, n_years_horizon=2)
    series.check(make_year(rng, next_year(fit)), rng=rng)
    series.check(make_year(rng, next_year(fit, 1)), rng=rng)
    with pytest.raises(ValueError, match="used up"):
        series.check(make_year(rng, next_year(fit, 2)), rng=rng)
    fresh = series.renew()
    assert fresh.remaining == 2 and len(fresh.checks) == 0 and len(series.checks) == 2


def test_year_index_is_read_from_the_dates_under_a_trend():
    rng = np.random.default_rng(502)
    hist = pd.concat([pd.DataFrame({"date": pd.date_range(f"{y}-01-01", periods=365, freq="D"),
                                    "value": 10 + 0.5 * (y - 2010) + rng.normal(0, 1, 365)})
                      for y in range(2010, 2020)])
    fit = fit_historical(hist, "date", "value", TideConfig(n_bootstrap=2000))
    assert fit.trend["applied"]
    series = MonitoringSeries(fit, n_years_horizon=5)
    y2024 = pd.DataFrame({"date": pd.date_range("2024-01-01", periods=365, freq="D"),
                          "value": 10 + 0.5 * 14 + rng.normal(0, 1, 365)})
    check = series.check(y2024, rng=1)
    assert check.year_index == 14
    assert not check.flagged


def test_save_load_round_trip(tmp_path):
    fit, rng = make_fit(seed=503)
    series = MonitoringSeries(fit, alpha_total=0.05, n_years_horizon=5, alternative="greater")
    series.check(make_year(rng, next_year(fit)), label="first", rng=rng)
    series.check(make_year(rng, next_year(fit, 1)), label="second", rng=rng)
    path = str(tmp_path / "state.json")
    series.save(path)
    reloaded = MonitoringSeries.load(fit, path)
    assert reloaded.alpha_total == series.alpha_total
    assert reloaded.n_years_horizon == series.n_years_horizon
    assert reloaded.alternative == "greater"
    assert [c.label for c in reloaded.checks] == ["first", "second"]
    assert reloaded.checks[1].p_value == series.checks[1].p_value


def test_horizon_wide_false_alarm_rate():
    """Chance of EVER flagging a normal year over a 10-year horizon at a 5%
    budget, 20 historical years. Measured at 0.07-0.09 in
    tests/calibration_study.py (300 horizons); the bound here is loose
    because this pools only 80 horizons."""
    ever = total = 0
    for f in range(20):
        rng = np.random.default_rng(600 + f)
        fit = fit_historical(make_historical(rng, n_years=20), "date", "value",
                             TideConfig(n_bootstrap=2000))
        for h in range(4):
            series = MonitoringSeries(fit, alpha_total=0.05, n_years_horizon=10)
            start = next_year(fit, 10 * h)
            ever += any(series.check(make_year(rng, start + i), rng=rng).flagged for i in range(10))
            total += 1
    assert ever / total <= 0.20, f"ever-flagged rate {ever / total:.3f}"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
