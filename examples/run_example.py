"""
End-to-end usage example (Phase 6 of the rebuild plan).

Generates synthetic seasonal, autocorrelated data standing in for real
water quality data, then walks through every piece: fit -> test -> plot ->
Sequential -> Cumulative -> MonitoringSeries -> power -> sensitivity. Read
this top to bottom as a tutorial; swap the synthetic data generator for a
real CSV load when you're ready (see the "swap in real data" note near
the top).

Run: python examples/run_example.py
Output: examples/envelope_plot.png, sequential_plot.png,
        cumulative_plot.png, monitoring_plot.png
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import pandas as pd

from tide_lite import (
    TideConfig, fit_historical, test_treatment, bin_significant,
    estimate_departure_recovery,
    SequentialEvent, sequential_test,
    CumulativeWindow, cumulative_test, first_recovery,
    MonitoringSeries,
    estimate_power, estimate_power_sequential, sensitivity_grid,
    plot_envelope, plot_sequential, plot_cumulative, plot_monitoring,
)

rng = np.random.default_rng(0)
OUT_DIR = os.path.dirname(__file__)


# ---------------------------------------------------------------------------
# 1. Data. Swap this whole section for:
#      df = pd.read_csv("your_data.csv")
#    as long as df has a date column and a value column, everything below
#    is unchanged -- fit_historical/test_treatment only need those two
#    column names (see date_col/value_col below).
# ---------------------------------------------------------------------------

def make_year(year, effect=0.0, n_days=365, base=10.0, seasonal_amp=3.0,
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


historical = pd.concat([make_year(y) for y in range(2010, 2020)], ignore_index=True)
treatment_2021 = make_year(2021, effect=2.5)  # the year we actually want an answer for

DATE_COL, VALUE_COL = "date", "value"


# ---------------------------------------------------------------------------
# 2. Phase 1 decisions, in one place (this is what your real METHODS.md
#    should record, with your own reasoning for each choice).
# ---------------------------------------------------------------------------

config = TideConfig(
    bin_days=30,               # ~monthly bins
    agg="median",
    detrend_mode="additive",   # use "log_additive" if your variable's
                                # spread scales with its level
    mk_alpha=0.10,
    max_missing_frac=0.5,
    n_bootstrap=2000,
)


# ---------------------------------------------------------------------------
# 3. Fit + test (Phase 2)
# ---------------------------------------------------------------------------

fit = fit_historical(historical, DATE_COL, VALUE_COL, config)
print(f"Historical years used: {fit.years_used}")
if fit.years_dropped:
    print(f"Years dropped (completeness filter): {fit.years_dropped}")
print(f"Trend check: Mann-Kendall p={fit.trend['mk_p']:.3f}, "
      f"applied={fit.trend['applied']}")
# How much of the leftover variation is whole-year level shifts ("this was a
# wet year") rather than within-year wiggle. The block bootstrap averages
# whole-year shifts away, so a high value here means the p-values run
# anti-conservative -- see "Calibration" in USER_GUIDE.md.
print(f"Between-year variance share: {fit.between_year_var_frac:.0%} "
      f"({'low -- well calibrated' if fit.between_year_var_frac < 0.3 else 'high -- see USER_GUIDE.md Calibration'})")
print(f"Smallest resolvable p-value at n_bootstrap={config.n_bootstrap}: "
      f"{config.min_attainable_p:.4f}\n")

result = test_treatment(fit, treatment_2021, mode="prediction",
                         treatment_year_index=fit.year_index_for(2021),
                         run_lengths=[3], rng=rng)
print(f"2021 vs. historical envelope: p={result.p_value:.4f}, "
      f"effect_size={result.effect_size:+.2f} (original units)\n")

# Per-bin (monthly) breakdown: the global p-value above is for the whole
# year; this asks which specific bins were driving it, corrected across
# bins so it doesn't reintroduce the multiple-comparisons problem the
# global test exists to avoid.
sig_bins = bin_significant(result, alpha=0.05)
print(f"Significant bins (Holm-Bonferroni corrected across bins): "
      f"{fit.bins[sig_bins].tolist()}")
print(f"Per-bin p-values (uncorrected): {np.round(result.bin_p_values, 4).tolist()}\n")

# Sustained-departure test (run_lengths=[3] above): was there a run of 3+
# consecutive months, together, that was unusual -- a different question
# from any single bin's own significance, and more robust to one noisy
# month. Then estimate roughly when it started and (if at all) recovered.
sustained_3 = result.sustained[3]
print(f"Sustained 3-month departure: p={sustained_3.p_value:.4f}, "
      f"core window bins={sustained_3.window_bins.tolist()}")
departure = estimate_departure_recovery(result, run_length=3)
print(f"Estimated departure: bins {departure.departure_start_bin}-{departure.departure_end_bin}, "
      f"recovered={departure.recovered}"
      + (f", by bin {departure.recovered_at_bin}\n" if departure.recovered else " (persisted through period end)\n"))


# ---------------------------------------------------------------------------
# 4. Plot: historical envelope + individual years + the treatment year
# ---------------------------------------------------------------------------

ax = plot_envelope(fit, result, treatment_label="2021", departure=departure)
out_path = os.path.join(OUT_DIR, "envelope_plot.png")
ax.figure.savefig(out_path, dpi=150)
print(f"Saved plot to {out_path}\n")


# ---------------------------------------------------------------------------
# 5. Sequential: two independent events tested together
# ---------------------------------------------------------------------------

events = [
    SequentialEvent(label="2020 (no effect)", treatment_df=make_year(2020, effect=0.0),
                     year_index=fit.year_index_for(2020)),
    SequentialEvent(label="2021 (real effect)", treatment_df=treatment_2021,
                     year_index=fit.year_index_for(2021)),
]
seq_results = sequential_test(fit, events, alpha=0.05, rng=rng)
# Note: "2020 (no effect)" flagging significant here on any given run isn't
# necessarily a bug -- at alpha=0.05, a true null crosses the threshold
# about 1 time in 20 purely by chance (see tests/test_engine.py for the
# calibration check across many replicate runs, where this averages out).
print("Sequential (Holm-Bonferroni corrected across both events):")
for row in seq_results:
    print(f"  {row.label}: p={row.result.p_value:.4f}  "
          f"raw_sig={row.significant_raw}  corrected_sig={row.significant_corrected}")

ax = plot_sequential(fit, seq_results)
out_path = os.path.join(OUT_DIR, "sequential_plot.png")
ax.figure.savefig(out_path, dpi=150)
print(f"Saved plot to {out_path}\n")


# ---------------------------------------------------------------------------
# 6. Cumulative: track the 2021 event's recovery over the following years
# ---------------------------------------------------------------------------

# year_index is the offset in CALENDAR years from the first historical
# year (2010 here), so 2021 -> 11. fit.year_index_for() does the arithmetic;
# prefer it over len(fit.years_used), which is only the same number when the
# historical record has no missing years.
windows = [
    CumulativeWindow("2021 (+0y)", treatment_2021, year_index=fit.year_index_for(2021)),
    CumulativeWindow("2022 (+1y)", make_year(2022, effect=1.0), year_index=fit.year_index_for(2022)),
    CumulativeWindow("2023 (+2y)", make_year(2023, effect=0.0), year_index=fit.year_index_for(2023)),
    CumulativeWindow("2024 (+3y)", make_year(2024, effect=0.0), year_index=fit.year_index_for(2024)),
]
# rng= makes the whole family reproducible: every p-value here comes from a
# Monte Carlo null, so without a seed a borderline window can change which
# side of alpha it lands on between runs -- and first_recovery() keys off
# exactly that.
cum_results = cumulative_test(fit, windows, alpha=0.05, rng=rng)
print("Cumulative (same baseline throughout, no correction):")
for row in cum_results:
    print(f"  {row.label}: p={row.result.p_value:.4f}  "
          f"effect_size={row.result.effect_size:+.2f}  significant={row.significant}")

# Formal recovery rule: require 2 consecutive non-significant windows
# before declaring recovered, not just one (a single non-significant
# check-in can be noise -- see the calibration notes in test_engine.py).
recovery = first_recovery(cum_results, consecutive_required=2)
print(f"Recovery: {recovery}\n")

ax = plot_cumulative(cum_results)
out_path = os.path.join(OUT_DIR, "cumulative_plot.png")
ax.figure.savefig(out_path, dpi=150)
print(f"Saved plot to {out_path}\n")


# ---------------------------------------------------------------------------
# 7. MonitoringSeries: ongoing yearly checks, no specific event, a fixed
#    lifetime budget split evenly across a committed horizon.
# ---------------------------------------------------------------------------

series = MonitoringSeries(fit, alpha_total=0.05, n_years_horizon=10)
series.check(make_year(2020, effect=0.0), label="2020", year_index=fit.year_index_for(2020), rng=rng)
series.check(make_year(2021, effect=0.0), label="2021", year_index=fit.year_index_for(2021), rng=rng)
series.check(make_year(2022, effect=3.0), label="2022", year_index=fit.year_index_for(2022), rng=rng)
print(f"MonitoringSeries (alpha_per_check={series.alpha_per_check:.4f}, "
      f"{series.remaining}/{series.n_years_horizon} checks remaining):")
for c in series.checks:
    print(f"  {c.label}: p={c.p_value:.4f}  flagged={c.flagged}")
series.save(os.path.join(OUT_DIR, "monitoring_state.json"))  # picks up here next year

ax = plot_monitoring(series)
out_path = os.path.join(OUT_DIR, "monitoring_plot.png")
ax.figure.savefig(out_path, dpi=150)
print(f"Saved plot to {out_path}\n")


# ---------------------------------------------------------------------------
# 8. Power analysis (Phase 5) -- keep n_simulations modest here for speed;
#    a few hundred is typically enough for a usable curve.
# ---------------------------------------------------------------------------

power_curve = estimate_power(fit, effect_sizes=[0, 1, 2, 3, 4], n_simulations=150, rng=rng)
print("Power curve (prediction mode):")
for effect, power in power_curve.items():
    print(f"  effect={effect:+.1f}: power={power:.2f}")
print()

seq_power = estimate_power_sequential(fit, effect_sizes=[2, 4], n_events=2,
                                       n_simulations=150, rng=rng)
print("Power under Sequential correction (n_events=2, conservative approximation):")
for effect, power in seq_power.items():
    print(f"  effect={effect:+.1f}: power={power:.2f}")
print()


# ---------------------------------------------------------------------------
# 9. Sensitivity grid (Phase 5) -- stress-test the 2021 result against
#    defensible alternative Phase 1 choices.
# ---------------------------------------------------------------------------

sens = sensitivity_grid(
    historical, treatment_2021, DATE_COL, VALUE_COL, config,
    variations={"bin_days": [7, 60], "detrend_mode": ["additive", "log_additive"]},
    treatment_year_index=fit.year_index_for(2021), rng=rng,
)
print("Sensitivity grid (is the 2021 result stable to alternative settings?):")
print(sens.to_string(index=False))
