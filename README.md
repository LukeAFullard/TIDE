# tide_lite

**Did this period's water quality look abnormal compared with your own
historical "normal" years — accounting for seasonality and
autocorrelation, without assuming a distribution?**

One site, one variable, one question. No control site required.

```python
from tide_lite import fit_historical, test_treatment

fit    = fit_historical(historical_df, date_col="date", value_col="value")
result = test_treatment(fit, treatment_df, rng=42)

print(result.p_value, result.effect_size)
```

## Four modes, one question each

| Mode | Question it answers |
|---|---|
| **Standard** (`test_treatment`) | Did this one period behave abnormally? |
| **Sequential** (`sequential_test`) | Across several separate events, which were real? (Holm-Bonferroni corrected.) |
| **Cumulative** (`cumulative_test`) | How is this one known event progressing — has it recovered? (Fixed baseline, uncorrected by design.) |
| **MonitoringSeries** (`MonitoringSeries`) | Is each new year still normal, checked year after year indefinitely, without false-alarm risk piling up? Stateful. Needs ~20+ historical years. |

Standard also answers **which months specifically** (`bin_significant`),
**was it a sustained multi-month shift** (`run_lengths=` →
`result.sustained`), and **when did it start and did it recover**
(`estimate_departure_recovery`).

Alongside those: `estimate_power` / `estimate_power_sequential` (can this
test even detect an effect you'd care about, given the history you
actually have?) and `sensitivity_grid` (is the result stable to
defensible alternative settings?).

## Setup

```
pip install -e .
```

Four dependencies: numpy, pandas, scipy, matplotlib. Mann-Kendall and
Sen's slope are implemented directly rather than pulling in another
package.

## Start here

```
python examples/run_example.py    # full walkthrough on synthetic data, writes four plots
python -m pytest tests/ -q        # 38 validation and regression checks
```

**New here? Read `QUICKSTART.md` (5 minutes), then `USER_GUIDE.md`** for
the complete picture: what each mode is for, how to read a result, how
much to trust it, and a troubleshooting entry for every error the code
raises.

`examples/run_example.py` is written as a tutorial — read it top to
bottom. Swap the synthetic data generator at the top for a real
`pd.read_csv(...)` and everything downstream is unchanged.

## Two things to know before you interpret a real result

**1. Calibration depends on your data, and the tool tells you which case
you're in.** Measured false-positive rate at a nominal alpha of 0.05:

| Historical years | Within-year noise only | With whole-year level shifts |
|---|---|---|
| 10 | 0.054 | **0.112 – 0.130** |
| 20 | 0.042 | **0.072 – 0.081** |
| 40 | 0.018 | — |

When year-to-year variation is just accumulated within-year wiggle, the
test is well calibrated. When your data also carries genuine whole-year
level shifts — a wet year sitting high all year, the norm in hydrology —
it runs at roughly twice its nominal false-positive rate on a 10-year
record. `fit_historical` reports `between_year_var_frac` so you know
which row applies, and warns when you're in the bad regime. Full
explanation, including why the obvious fix doesn't work, in
`USER_GUIDE.md` § 8 and the `tide_lite/engine.py` docstring.

**2. With ~10 historical years, a single fitted model's false-positive
rate varies a lot fit to fit** (0–25% at nominal alpha=0.05 in
simulation) even though it averages out correctly across many possible
historical samples. Your real analysis only ever has one historical fit.
This is exactly why `estimate_power` and `sensitivity_grid` exist — run
them on a real result before leaning on a single p-value.

## Layout

```
tide_lite/
  engine.py       -- fit_historical, test_treatment (core), plus the
                     calibration study behind the table above
  controllers.py  -- Sequential, Cumulative, Holm-Bonferroni
  power.py        -- power analysis, sensitivity harness
  monitoring.py   -- MonitoringSeries: ongoing yearly monitoring
  plotting.py     -- a dedicated plot function per mode
tests/
  test_engine.py      -- statistical validation of the method
  test_monitoring.py  -- MonitoringSeries validation
  test_regressions.py -- one test per implementation bug found in audit
examples/
  run_example.py  -- runnable tutorial, all four modes plotted
QUICKSTART.md     -- 5-minute start
USER_GUIDE.md     -- modes, interpretation, calibration, limitations, troubleshooting
```
