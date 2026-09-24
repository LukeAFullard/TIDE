# Quick start

A first result in about five minutes.

## 1. Install

From the repository folder:

```
pip install -e .
python examples/run_example.py
```

The example builds a synthetic site, runs every feature and saves plots in
`examples/`. If it prints results, you are set up.

## 2. Your data

Two tables, each with a **date** column and a **value** column (CSV,
Excel, database query: anything pandas can read):

- **history**: the years you consider normal for this site (10 or more if
  possible; 20+ is better);
- **treatment**: the period you are asking about (one calendar year, or
  part of one).

Values must be numbers. Convert results such as `<0.5` (below detection)
to a number first, the same way for every year.

## 3. Run it

```python
import pandas as pd
from tide_lite import fit_historical, test_treatment, summarize, plot_envelope

history   = pd.read_csv("history.csv")
treatment = pd.read_csv("2024.csv")

fit    = fit_historical(history, date_col="date", value_col="value")
result = test_treatment(fit, treatment, rng=42)
print(summarize(fit, result))
```

Example output:

```
Question: were 2024 values higher or lower than the site's normal range (20 historical years, 2004-2023)?
Answer: YES at alpha=0.05 (p = 0.0015; at most alpha).
Bins outside the normal range: Mar, Sep, Oct, Dec.
Typical difference from the historical level: +1.637 (units of 'value'; median over 12 tested bins).
Notes: ...
Reproduce with: tide_lite 0.3.0; TideConfig(...); rng=42; fit fingerprint 5c4d7e515a4e4380.
```

- **p-value**: how often a normal year would look at least this unusual.
  At or below 0.05 is the usual bar for "unusual".
- **Bins**: the months that fell outside the normal range.
- **Typical difference**: how far above (+) or below (−) normal, in your units.
- **rng=42** fixes the random numbers so the result can be reproduced
  exactly. Use any number, but record it.

## 4. Look at it

```python
ax = plot_envelope(fit, result, treatment_label="2024")
ax.figure.savefig("2024.png", dpi=150)
```

The shaded band is where a normal year stays, all year, 95% of the time.
The red line leaving the band is exactly the same statement as p ≤ 0.05;
stars mark the months that left it.

## 5. Before you rely on it

1. Read the **Notes** in the summary. They list what the tool dropped and
   which cautions apply to your data.
2. If the answer is **NO**, run `estimate_power(fit, [1, 2, 4])` (sizes in
   your units) to see how big a change your record could have detected.
3. If the answer will be used in a decision, work through the checklist in
   `USER_GUIDE.md` section 7.

## 6. More than one period?

- Several separate years, "which were unusual?" → `sequential_test`
- Following one known event over the next years → `cumulative_test`
- Checking each new year as it arrives → `MonitoringSeries`

See `USER_GUIDE.md` section 4 to choose, and `examples/run_example.py` for
working code.
