# Quick Start

Get a real result in about 5 minutes.

## 1. Install

```
pip install numpy pandas scipy matplotlib
```

Nothing else — no extra packages.

## 2. See it work

```
python examples/run_example.py
```

This generates synthetic data, runs all three modes, and saves a plot to
`examples/envelope_plot.png`. If it prints results and the plot appears,
everything is working.

## 3. Run it on your own data

You need a table with a date column and a value column. Historical
("normal") years and the period you want to test can be two separate
files or two slices of one file — either is fine.

```python
import pandas as pd
from tide_lite import fit_historical, test_treatment

historical = pd.read_csv("historical.csv")   # your "normal" years
treatment = pd.read_csv("treatment_year.csv")  # the period you're asking about

fit = fit_historical(historical, date_col="date", value_col="value")
result = test_treatment(fit, treatment, mode="prediction")

print(f"p = {result.p_value:.4f}")
print(f"effect size = {result.effect_size:+.2f}")
```

A small p-value (conventionally < 0.05) means the treatment period looked
unusual compared to your historical years. `effect_size` tells you how
different, in your data's own units.

## 4. Before you trust that result

Read **`USER_GUIDE.md`** — specifically the **Limitations** section. The
short version: with fewer than ~15-20 historical years, one p-value on
its own is noisier than it looks, and this tool comes with a power/
sensitivity check built in for exactly that reason. Two minutes there will
save you from over-reading a single number.

## 5. Testing more than one period, or tracking one over time

That's `Sequential` and `Cumulative` — see **"Which mode do I use?"** in
`USER_GUIDE.md` for a two-question decision guide, or jump straight to
`examples/run_example.py` for working code for both.
