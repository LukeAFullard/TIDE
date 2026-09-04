# Quick Start

A real result in about 5 minutes.

## 1. Install

```
pip install -e .
```

from the repository root. Four dependencies (numpy, pandas, scipy,
matplotlib), nothing else to configure.

## 2. See it work

```
python examples/run_example.py
```

This generates synthetic data, runs all four modes, and saves four plots
into `examples/`. If it prints results and the plots appear, you're set.

## 3. Run it on your own data

You need a table with a date column and a value column. Historical
("normal") years and the period you're asking about can be two files or
two slices of one — either is fine.

```python
import pandas as pd
from tide_lite import fit_historical, test_treatment

historical = pd.read_csv("historical.csv")       # your "normal" years
treatment  = pd.read_csv("treatment_year.csv")   # the period you're asking about

fit    = fit_historical(historical, date_col="date", value_col="value")
result = test_treatment(fit, treatment, mode="prediction", rng=42)

print(f"p = {result.p_value:.4f}")
print(f"effect size = {result.effect_size:+.2f}")
```

A small p-value (conventionally < 0.05) means the period looked unusual
compared with your historical years. `effect_size` says how different, in
your data's own units.

**Pass `rng=` for anything you'll report.** Every p-value comes from a
Monte Carlo null, so without a seed two runs differ slightly and a
borderline result can flip sides.

## 4. Three checks before you trust that number

```python
print(fit.years_used, fit.years_dropped)   # did a year get filtered out?
print(fit.trend)                            # was a long-term drift removed?
print(fit.between_year_var_frac)            # which calibration regime are you in?
```

- **`years_dropped`** — years missing more than half their bins are
  dropped. Make sure that's what you wanted.
- **`trend["applied"]`** — if True, quote
  `result.effect_size_trend_adjusted`, which is what the p-value actually
  tested, and read the `year_index` note in `USER_GUIDE.md` § 6.
- **`between_year_var_frac`** — above ~0.5 with fewer than 20 historical
  years, the test runs at roughly twice its nominal false-positive rate.
  It's the single most important number for deciding how much weight a
  p-value near 0.05 can carry. `USER_GUIDE.md` § 8 explains it.

Then read the **Limitations** section of `USER_GUIDE.md`. Short version:
with fewer than ~15–20 historical years one p-value is noisier than it
looks, and this tool ships power and sensitivity checks for exactly that
reason.

## 5. Plot it

```python
from tide_lite import plot_envelope
ax = plot_envelope(fit, result, treatment_label="2024")
ax.figure.savefig("result.png", dpi=150)
```

Does the picture match the p-value? If a small p-value comes with a line
that barely leaves the band — or a large one with a line well outside —
stop and investigate before reporting.

## 6. Testing more than one period

- Several **separate events**, and you want to know which were real →
  `sequential_test`.
- **One known event** tracked over time, has it recovered →
  `cumulative_test` + `first_recovery`.
- **No event at all**, just checking each new year as it arrives →
  `MonitoringSeries` (needs ~20+ historical years).

See **"Which mode do I use?"** in `USER_GUIDE.md` for a decision guide,
or `examples/run_example.py` for working code for all four.
