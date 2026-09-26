"""
End-to-end tutorial on synthetic data. Read it top to bottom.

To use your own data, replace section 1 with something like
    historical = pd.read_csv("history.csv")      # columns: date, value
    treatment  = pd.read_csv("2021.csv")
Everything below only needs a date column and a value column.

Run:    python examples/run_example.py
Writes: examples/envelope_plot.png, sequential_plot.png,
        cumulative_plot.png, monitoring_plot.png, monitoring_state.json
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import pandas as pd
from scipy.signal import lfilter

from tide_lite import (
    TideConfig, fit_historical, test_treatment, bin_significant, summarize,
    estimate_departure_recovery,
    SequentialEvent, sequential_test,
    CumulativeWindow, cumulative_test, first_recovery,
    MonitoringSeries,
    estimate_power, sensitivity_grid,
    plot_envelope, plot_sequential, plot_cumulative, plot_monitoring,
)

OUT_DIR = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(0)


# ---------------------------------------------------------------------------
# 1. Data: one site, one variable, daily values. A seasonal cycle, day-to-day
#    noise that persists for a few days, and a random whole-year shift
#    (wet and dry years). `effect` adds a change we want the test to find.
# ---------------------------------------------------------------------------

def make_year(year, effect=0.0):
    days = np.arange(1, 366)
    seasonal = 10.0 + 3.0 * np.sin(2 * np.pi * (days - 60) / 365)
    noise = lfilter([1.0], [1.0, -0.7], rng.normal(0, 0.6, 365))
    values = seasonal + noise + rng.normal(0, 0.4) + effect
    return pd.DataFrame({"date": pd.date_range(f"{year}-01-01", periods=365, freq="D"),
                         "value": values})


historical = pd.concat([make_year(y) for y in range(2000, 2020)], ignore_index=True)
treatment_2021 = make_year(2021, effect=2.0)        # the year we want an answer for


# ---------------------------------------------------------------------------
# 2. Settings. Decide these BEFORE looking at the treatment data and record
#    them with the result. The defaults are sensible for most records.
# ---------------------------------------------------------------------------

config = TideConfig(
    bin_days="month",          # compare month with month (7 = weekly, 30 = 30-day)
    detrend_mode="additive",   # "log_additive" if changes are proportional (values > 0)
    n_bootstrap=4000,          # 4000 so MonitoringSeries (section 7) can resolve its threshold
)


# ---------------------------------------------------------------------------
# 3. Fit the historical record, then check what it did with your data.
# ---------------------------------------------------------------------------

fit = fit_historical(historical, "date", "value", config)
print(f"Historical years used: {fit.years_used[0]}-{fit.years_used[-1]} ({len(fit.years_used)})")
print(f"Years dropped: {fit.years_dropped or 'none'};  bins dropped: {fit.bins_dropped or 'none'}")
print(f"Trend removed: {fit.trend['applied']} (Mann-Kendall p={fit.trend['mk_p']:.3f})")
print(f"Share of variation that is whole-year shifts: {fit.between_year_var_frac:.0%}\n")


# ---------------------------------------------------------------------------
# 4. Standard test: was 2021 unusual? Always pass a seed (rng=) for a result
#    you will report, so it can be reproduced exactly.
# ---------------------------------------------------------------------------

result = test_treatment(fit, treatment_2021, run_lengths=[3], rng=42)
print(summarize(fit, result, alpha=0.05))
print()
print(f"Months flagged: {fit.bins[bin_significant(result, 0.05)].tolist()}")
print(f"Sustained 3-month departure: p={result.sustained[3].p_value:.4f}, "
      f"months {result.sustained[3].window_bins.tolist()}")
departure = estimate_departure_recovery(result, run_length=3)
print(f"Departure (descriptive): months {departure.departure_start_bin}-"
      f"{departure.departure_end_bin}, recovered={departure.recovered}\n")

ax = plot_envelope(fit, result, alpha=0.05, treatment_label="2021", departure=departure)
ax.figure.savefig(os.path.join(OUT_DIR, "envelope_plot.png"), dpi=150)


# ---------------------------------------------------------------------------
# 5. Sequential: several separate years at this site, corrected together.
# ---------------------------------------------------------------------------

events = [
    SequentialEvent("2020 (no change)", make_year(2020)),
    SequentialEvent("2021 (change)", treatment_2021),
    SequentialEvent("2022 (no change)", make_year(2022)),
]
seq = sequential_test(fit, events, alpha=0.05, rng=42)
print("Sequential (Holm-Bonferroni across the three years):")
for row in seq:
    print(f"  {row.label:18s} p={row.result.p_value:.4f}  flagged={row.significant_corrected}")
ax = plot_sequential(fit, seq)
ax.figure.savefig(os.path.join(OUT_DIR, "sequential_plot.png"), dpi=150)


# ---------------------------------------------------------------------------
# 6. Cumulative: follow one known event over the following years.
# ---------------------------------------------------------------------------

windows = [
    CumulativeWindow("2021", treatment_2021),
    CumulativeWindow("2022", make_year(2022, effect=1.0)),
    CumulativeWindow("2023", make_year(2023)),
    CumulativeWindow("2024", make_year(2024)),
]
cum = cumulative_test(fit, windows, alpha=0.05, rng=42)
print("\nCumulative (same baseline, no correction):")
for row in cum:
    print(f"  {row.label}: p={row.result.p_value:.4f}  "
          f"effect={row.result.effect_size_trend_adjusted:+.2f}  significant={row.significant}")
rec = first_recovery(cum, consecutive_required=2)
print(f"No longer distinguishable from normal from: {rec.recovered_at_label}")
ax = plot_cumulative(cum)
ax.figure.savefig(os.path.join(OUT_DIR, "cumulative_plot.png"), dpi=150)


# ---------------------------------------------------------------------------
# 7. MonitoringSeries: check each new year as it arrives, with a 5% budget
#    for the whole 10-year horizon (0.5% per check). Save between runs.
#    These are freshly simulated years (only 2022 changed), separate from
#    the 2021 data tested above.
# ---------------------------------------------------------------------------

series = MonitoringSeries(fit, alpha_total=0.05, n_years_horizon=10)
for year, effect in ((2020, 0.0), (2021, 0.0), (2022, 3.0)):
    series.check(make_year(year, effect), rng=year)
print(f"\nMonitoring ({series.alpha_per_check:.4f} per check, {series.remaining} checks left):")
for c in series.checks:
    print(f"  {c.label}: p={c.p_value:.4f}  flagged={c.flagged}")
series.save(os.path.join(OUT_DIR, "monitoring_state.json"))   # reload next year with .load
ax = plot_monitoring(series)
ax.figure.savefig(os.path.join(OUT_DIR, "monitoring_plot.png"), dpi=150)


# ---------------------------------------------------------------------------
# 8. Power: what size of change could this record detect? Effects are in the
#    data's units (log units if detrend_mode="log_additive"). The rate at
#    +0.0 is the false-alarm rate, about alpha (0.05) give or take a point
#    or two of random-draw noise.
# ---------------------------------------------------------------------------

power = estimate_power(fit, effect_sizes=[0.0, 0.5, 1.0, 1.5, 2.0], n_simulations=400, rng=1)
print("\nChance of detecting a whole-year shift of a given size (alpha=0.05):")
for effect, p in power.items():
    print(f"  {effect:+.1f}: {p:.2f}")


# ---------------------------------------------------------------------------
# 9. Sensitivity: does the 2021 answer hold under other reasonable settings?
# ---------------------------------------------------------------------------

grid = sensitivity_grid(historical, treatment_2021, "date", "value", config,
                        {"bin_days": [30, 7], "agg": ["mean"], "block_length_bins": [2, 4]},
                        rng=42)
print("\nSensitivity grid (same random seed in every row):")
print(grid[["variant", "p_value", "effect_size_trend_adjusted", "n_historical_years",
            "n_bins", "error"]].to_string(index=False))
