# tide_lite user guide

For a first result, read `QUICKSTART.md`. For the full method, its
assumptions and measured error rates, read `METHODS.md`. This guide covers
everything in between.

1. [What it answers](#1-what-it-answers)
2. [Preparing your data](#2-preparing-your-data)
3. [Settings: decide before you look](#3-settings-decide-before-you-look)
4. [Which mode?](#4-which-mode)
5. [Running each mode](#5-running-each-mode)
6. [Reading the results](#6-reading-the-results)
7. [Checklist before you report](#7-checklist-before-you-report)
8. [Troubleshooting](#8-troubleshooting)
9. [Glossary](#9-glossary)

---

## 1. What it answers

**"Was this period unusual compared with this site's own normal years?"**
Normal allows for the seasons, a steady long-term trend, and ordinary
year-to-year variation (wet years, dry years).

In brief, the tool:

1. summarises each year month by month (by default);
2. learns from your history what a normal year looks like and how much
   normal years vary, month by month and as a whole;
3. builds thousands of realistic synthetic normal years from pieces of
   your own history;
4. counts how often those synthetic normal years look at least as unusual
   as your period. That share is the **p-value**.

`METHODS.md` section 2 lists every step.

**It does not tell you:**

- **why** a period was unusual. A significant result says the period
  differed from normal, not that a particular activity caused it;
- whether values **meet a limit or guideline**. It compares the site with
  its own past; compare with a limit directly for that;
- that **nothing changed**, when the result is not significant. A short or
  noisy record can miss real changes; `estimate_power` says how big a
  change you could have detected.

## 2. Preparing your data

**Format.** Two tables, each with a date column and a value column (any
names; you pass them in). One site and one variable per analysis.

- **History**: the years you consider normal. Exclude any year affected by
  the activity you are testing, and any year you are unsure about.
- **Treatment**: the period you are asking about. It must be within one
  calendar year (January–December) for a single test. A part year is fine
  as long as at least half of its months have data.

**How much history?**

| Years | What to expect |
|---|---|
| under 3 | Refused. |
| 3–9 | Valid but only large changes are detectable. You get a warning. |
| 10–19 | Workable. Run `estimate_power`. |
| 20+ | Good. Needed for `MonitoringSeries`. |

**Measurement frequency.** Anything from continuous sensors to one grab
sample a month. Monthly samples suit the default monthly bins. A month in
a given year counts only if it has at least 75% of that month's usual
number of measurements (`min_bin_coverage`), so a 5-day "June" is never
compared with 30-day ones. This applies to history and treatment alike.

**Non-detects** (e.g. `<0.5`). Text values are refused with an error.
Replace them with a number, the same rule for every year (e.g. the
detection limit, or half of it), and say so in your report. If more than
half the values in a month are identical non-detects, you get a warning:
the spread in that month cannot be measured and results there are
approximate.

**Missing data.** By default a year may miss up to half its months
(`max_missing_frac`). Months missing from the treatment period are left
out of the test and named in a warning; they are never treated as
evidence either way.

**Zero or negative values** rule out `detrend_mode="log_additive"`; use
the default `"additive"`.

## 3. Settings: decide before you look

Choose settings **before** looking at the treatment data, and report them.
Trying several and keeping the one that gives the answer you want makes
the p-value meaningless. `sensitivity_grid` exists to show that the answer
does not depend on the choice.

```python
from tide_lite import TideConfig
config = TideConfig(bin_days="month", detrend_mode="additive", n_bootstrap=2000)
fit = fit_historical(history, "date", "value", config)
```

| Setting | Default | Change it when |
|---|---|---|
| `bin_days` | `"month"` | You have frequent data and care about weekly detail (`7`). |
| `detrend_mode` | `"additive"` | Values are positive and changes are naturally proportional (concentrations that vary by factors): `"log_additive"`. Use `"none"` to never remove a trend (see below). |
| `n_bootstrap` | `2000` | Use 4000+ for `MonitoringSeries`, or for reporting p-values below 0.01. |
| `max_missing_frac` | `0.5` | Rarely. Raising it admits sparser years. |
| `min_bin_coverage` | `0.75` | Rarely. |
| `mk_alpha` | `0.10` | Rarely. The trend is removed if Mann-Kendall p is below this. |

**Should a trend be removed?** By default, if the history shows a steady
rise or fall, it is removed and projected to the tested year, so a year
that just continues the trend counts as normal. That suits "was this year
unusual given how the site was already changing?". If your question is
"has the site changed from its historical level?", a continuing
deterioration *is* the change you are looking for: use
`detrend_mode="none"`. Decide which question you are asking before you
look. The further the tested year is from the end of the history, the
more a removed trend matters; you get a warning beyond 5 years.

And for each test:

| Argument | Default | Meaning |
|---|---|---|
| `alternative` | `"two-sided"` | `"greater"` if only an increase matters (e.g. "did concentrations rise?"); `"less"` for decreases. More sensitive in that direction, blind to the other. Decide in advance. |
| alpha (in `summarize`, `bin_significant`, plots, controllers) | `0.05` | The false-alarm rate you accept. 0.05 or 0.10 are well calibrated; at 0.01 the real rate can be up to about 3% (`METHODS.md` section 4). |
| `rng` | none | **Always pass a number** (e.g. `rng=42`) for a result you report. |

## 4. Which mode?

| Your question | Mode | Corrects for testing several periods? |
|---|---|---|
| Was this one period unusual? | **Standard**: `test_treatment` | Not needed |
| I have several separate years; which were unusual? | **Sequential**: `sequential_test` | Yes (Holm) |
| A known event happened; is its effect still visible in each following year? | **Cumulative**: `cumulative_test` | No, by design |
| No specific event; check every new year as it comes, for years | **MonitoringSeries** | Yes (fixed budget over a horizon) |

**Sequential vs Cumulative.** Sequential treats each year as its own claim
("2019 was unusual", "2021 was unusual") and controls the chance of *any*
false claim. Cumulative describes the course of one event you already
know about, so each year is read as a point on a recovery curve, not as a
new discovery. If each year's result must stand on its own as evidence,
use Sequential.

**Several sites or variables.** Fit and test each separately, then
correct across them: `holm_bonferroni([p1, p2, p3], alpha=0.05)`.

## 5. Running each mode

All modes start from one fit:

```python
from tide_lite import fit_historical, TideConfig
fit = fit_historical(history, "date", "value", TideConfig())
```

The fit is deterministic (no randomness). The test reads each period's
year from its dates and projects any trend to that year automatically.

### Standard

```python
from tide_lite import test_treatment, summarize, bin_significant
result = test_treatment(fit, df_2024, rng=42)
print(summarize(fit, result))
print(fit.bins[bin_significant(result, alpha=0.05)])     # months flagged
```

To test the **average of several years** (e.g. "were 2021–2023 as a
whole unusual?"): `test_treatment(fit, df_2021_to_2023, mode="confidence", rng=42)`.
A month missing in any of those years is left out.

**Sustained departure.** "Was there a run of 3 months that was
consistently high (or low)?":

```python
result = test_treatment(fit, df_2024, run_lengths=[3], rng=42)
result.sustained[3].p_value, result.sustained[3].window_bins
```

This is more sensitive than the main test to a moderate shift lasting
several months, and less sensitive to a single-month spike. Decide which
you are testing in advance. `estimate_departure_recovery(result, 3)` then
describes roughly when the departure started and ended (descriptive, not
a test; it tends to over-reach by about a month).

### Sequential

```python
from tide_lite import SequentialEvent, sequential_test
events = [SequentialEvent("2019", df_2019), SequentialEvent("2021", df_2021)]
for row in sequential_test(fit, events, alpha=0.05, rng=42):
    print(row.label, row.result.p_value, row.significant_corrected)
```

Report `significant_corrected`. Events must be different calendar years;
overlap raises an error.

### Cumulative

```python
from tide_lite import CumulativeWindow, cumulative_test, first_recovery
windows = [CumulativeWindow(str(y), df[y]) for y in (2022, 2023, 2024, 2025)]
rows = cumulative_test(fit, windows, alpha=0.05, rng=42)
first_recovery(rows, consecutive_required=2)
```

`first_recovery` reports the first point after which the effect was not
detectable for 2 years in a row. That means "no longer distinguishable
from normal", not proven recovery. Quote `estimate_power` alongside it.

### MonitoringSeries

```python
from tide_lite import MonitoringSeries
fit = fit_historical(history, "date", "value", TideConfig(n_bootstrap=4000))
series = MonitoringSeries(fit, alpha_total=0.05, n_years_horizon=10)
series.check(df_2024, rng=2024)
series.save("monitoring.json")

# next year, same history and same settings:
fit = fit_historical(history, "date", "value", TideConfig(n_bootstrap=4000))
series = MonitoringSeries.load(fit, "monitoring.json")
series.check(df_2025, rng=2025)
```

`alpha_total` is the budget for the **whole horizon**: each check uses
`alpha_total / n_years_horizon` (0.005 here). Loading with a different fit
(other data or settings) raises an error, because changing the baseline
part-way breaks the budget. When the horizon is used up, call
`series.renew()`. Measured over 10 years the chance of ever flagging a
normal year was 5–10% for a 5% budget with 20–40 years of history: close
to, but not within, the budget.

### Power and sensitivity

```python
from tide_lite import estimate_power, sensitivity_grid
estimate_power(fit, effect_sizes=[0.5, 1, 2], rng=1)
# e.g. {0.5: 0.17, 1.0: 0.5, 2.0: 0.98}: share of changed years detected
sensitivity_grid(history, df_2024, "date", "value", TideConfig(),
                 {"bin_days": [30, 7], "detrend_mode": ["log_additive"]}, rng=42)
```

Effect sizes are whole-year shifts in your units (natural-log units with
`log_additive`: 0.1 ≈ +10%). Power at effect 0 is about alpha by
construction. In `sensitivity_grid` every row uses the same seed, so
differences come from the setting alone. Watch `n_historical_years` and
`n_bins`: a variant that drops data is not a like-for-like comparison.

## 6. Reading the results

### The fit

| Field | Check |
|---|---|
| `fit.years_used`, `fit.years_dropped` | Were any years left out that you expected in? |
| `fit.bins_dropped` | Were any months left out? |
| `fit.trend` | `applied`, `slope_per_year` (log units if `log_space`), `mk_p`. If a trend was removed, the test projects it to the tested year. |
| `fit.between_year_var_frac` | Share of normal variation that is whole-year shifts. Above about 0.5: a change affecting the whole year is hard to detect, because it looks like a wet or dry year. |
| `fit.fingerprint()` | Identifies the fitted baseline; record it. |

### The result

| Field | Meaning |
|---|---|
| `p_value` | How often a normal year looks at least this unusual. Smallest possible: 1/(n_bootstrap+1). |
| `effect_size_trend_adjusted` | Typical difference from what was expected for that year (median over months), your units. **Quote this one.** |
| `effect_size` | Same, against the record's overall typical level without the trend. Equal to the above when no trend was removed. |
| `bin_scores` | Each month's difference divided by that month's normal spread. |
| `bin_p_adjusted` | Per month, already adjusted for looking at all months: compare with alpha directly. `bin_significant()` does this. |
| `bin_p_values` | Per month, *not* adjusted. Descriptive only. |
| `missing_bins` | Months with no data; not tested. |
| `treatment_years`, `treatment_year_index` | Which year(s) were tested and where on the trend line. |
| `seed` | The seed passed as `rng`. |

The period is significant **if and only if** at least one month is
flagged, and the months flagged are exactly those outside the band in
`plot_envelope`.

### The plots

- `plot_envelope(fit, result)`: thin lines are your historical years; the
  shaded **significance band** is where a normal year stays in every month
  95% of the time (for alpha=0.05). Leaving it anywhere is the same
  statement as p ≤ 0.05.
- `plot_sequential(fit, rows)`: each event against the historical 10th–90th
  percentile, coloured by its corrected result.
- `plot_cumulative(rows)`: effect size per year, coloured by significance.
- `plot_monitoring(series)`: each check's p-value against its threshold.

## 7. Checklist before you report

1. Settings, alternative and alpha were chosen before looking at the
   treatment data.
2. `rng=` was set, and the `summarize` output (settings, seed,
   fingerprint, version) is kept with the result.
3. The historical years are all genuinely normal for the question, with no
   change of site, method, laboratory or units within them.
4. `years_dropped`, `bins_dropped` and any warnings were read and are
   acceptable.
5. If a trend was removed, it is plausible, removing it fits the question
   (section 3), and `effect_size_trend_adjusted` is the number quoted.
6. The plot agrees with the answer (it must, by construction; if it
   seems not to, check you are looking at the right alpha).
7. For a NO: `estimate_power` shows what size of change could have been
   detected. Report it.
8. `sensitivity_grid` gives the same answer under other reasonable
   settings.
9. Other explanations (weather, flow, upstream events, method changes)
   were considered before attributing a cause.
10. Several sites, variables or periods tested? The correction across them
    was applied.

## 8. Troubleshooting

**Errors**

| Message begins | Meaning and fix |
|---|---|
| `Column '...' not found` | Check `date_col` / `value_col` spelling. |
| `... value(s) in '...' are not numbers` | Text such as `<0.5`. Convert non-detects to numbers consistently (section 2). |
| `Only N historical years survived` | Fewer than 3 usable years. Add history, use coarser bins, or check dates. |
| `No bin has data in enough years` | Usually a date-parsing or column problem. |
| `detrend_mode='log_additive' needs strictly positive values` | Zeros or negatives present. Use `"additive"`. |
| `mode='prediction' tests one calendar year` | The treatment spans more than one calendar year. Split a period crossing 1 January (e.g. a water year) into its parts and test each (Sequential), or use `mode="confidence"` for the average of complete years. |
| `Treatment year(s) ... are also in the historical record` | Remove those years from the history and refit. |
| `The treatment data has no values in any of the fitted bins` | Column or date mismatch, or data outside the months the history covers. |
| `The treatment period is missing N of M bins` | Too sparse to compare with a whole year. Supply more data or use coarser bins. |
| `treatment_year_index cannot be overridden in confidence mode` | Remove the argument; each year's own date is used. |
| `No run of k consecutive bins is free of missing data` | Use a shorter `run_lengths` value. |
| `Events ... both contain YYYY data` | Sequential events overlap. Each must be a different year. |
| `This horizon cannot flag anything` | `alpha_total / n_years_horizon` is below the smallest possible p-value. Raise `n_bootstrap` as suggested. |
| `The monitoring horizon of N checks is used up` | Call `.renew()`. |
| `Year(s) ... have already been checked` | Each year is checked once per horizon. |
| `This monitoring state was created with a different historical fit` | Re-create the fit with exactly the same history and `TideConfig` as before. |
| `consecutive_required must be >= 1` | Use 2 or more. |
| `... must be one of ...` / `... must be between ...` | A misspelled or impossible setting. |

**Warnings** (the analysis ran; read them)

| Message begins | Meaning |
|---|---|
| `N year-bin cell(s) had fewer than 75%` / `N treatment bin(s) had fewer than 75%` | Those months had too few measurements and were treated as missing. |
| `Bins left out of the analysis` | Some months have too little history and are not tested. |
| `Only N historical years` | Valid, but only large changes are detectable. |
| `In bins [...] at least half the historical values are identical` | Usually non-detects. Results in those months are approximate. |
| `No treatment data in bins` | Those months are not tested. |
| `treatment_year_index=... but the data is dated ...` | You overrode the year; check this was intended. |
| `Nothing can be flagged at alpha=...` | `n_bootstrap` too small for that alpha. |
| `The bin after the departure ... has no treatment data` | "Recovered" only means the departure stopped being visible. |
| `N treatment bin(s) have more than twice the usual number of measurements` | The tested period was sampled far more often than the history (e.g. a sensor against monthly grab samples). Test a subset that matches the historical frequency. |
| `The historical trend is projected N years beyond the historical record` | The answer depends on a trend projected far forward. Also run with `detrend_mode="none"` and report both. |
| `MonitoringSeries is fitted on N historical years` / `alpha_per_check ... rests on fewer than 10` | Use 20+ years and a larger `n_bootstrap`. |

**Two runs give different p-values.** Pass `rng=` with a fixed number.

## 9. Glossary

- **Bin**: a slice of the year (a calendar month by default) that values
  are summarised into.
- **p-value**: how often a normal year would look at least as unusual as
  the tested period. Small means unusual.
- **alpha**: the p-value at or below which you call a result significant;
  also the false-alarm rate you accept (0.05 = 1 in 20 normal years
  flagged).
- **Significance band**: the range a normal year stays within in every
  month with probability 1 − alpha.
- **Effect size**: how far from normal, in your own units.
- **Power**: the chance of detecting a change of a given size.
- **Trend**: a steady long-term rise or fall, detected with Mann-Kendall
  and measured with Sen's slope.
- **Whole-year shift**: a year sitting high or low throughout (e.g. a wet
  year).
- **Median / MAD**: the middle value, and the typical distance from it.
  Both are robust to odd values.
- **Bootstrap**: building many realistic synthetic years from pieces of
  your own data instead of assuming a textbook distribution.
- **Holm correction**: adjusts for testing several periods so that the
  chance of any false alarm stays at alpha.
- **Non-detect**: a result below the detection limit, reported as `<DL`.
