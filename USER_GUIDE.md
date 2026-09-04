# TIDE lite — User Guide

A complete guide to using this tool, written to be useful whether you're
running the code yourself or reading someone else's results. If you just
want to get started, read `QUICKSTART.md` first (5 minutes) — this is the
comprehensive version.

## Contents

1. [What is this, and what question does it answer?](#1-what-is-this-and-what-question-does-it-answer)
2. [How to use this guide](#2-how-to-use-this-guide)
3. [Setup](#3-setup)
4. [Preparing your data](#4-preparing-your-data)
5. [Which mode do I use?](#5-which-mode-do-i-use)
6. [Running each mode](#6-running-each-mode)
7. [Reading the output: every number, in order](#7-reading-the-output-every-number-in-order)
8. [Calibration — how much can you trust the p-value?](#8-calibration--how-much-can-you-trust-the-p-value)
9. [Limitations — read this before trusting a result](#9-limitations--read-this-before-trusting-a-result)
10. [A checklist before you report a result](#10-a-checklist-before-you-report-a-result)
11. [Troubleshooting](#11-troubleshooting)
12. [Additional info](#12-additional-info)
13. [Glossary](#13-glossary)

---

## 1. What is this, and what question does it answer?

You have water quality data covering many "normal" years. Something
happened — a restoration project, a treatment plant upgrade, a policy
change, a pollution event — during or after a specific period. You want
to know: **did that period actually look different from what's normal,
once you account for the fact that water quality naturally moves with the
seasons and from year to year?**

This tool answers that one question with a statistical test built for
exactly this situation:

- It doesn't need a "control site" to compare against — many real
  monitoring situations only have one site.
- It respects seasonality automatically — it never compares a July value
  to a January baseline and calls that unusual.
- It doesn't assume your data follows a bell curve — real concentration
  data is often skewed, and this test doesn't care.

The output is a p-value (how surprising the observed period would be if
nothing had really changed), an effect size (how different, in your
data's own units), and a plot showing the period against the historical
range.

### What it does NOT answer

Being clear about this upfront saves more trouble than anything else in
this guide:

- **It does not prove causation.** A significant result says the period
  was unusual, not that your treatment made it unusual. See
  [Limitations](#9-limitations--read-this-before-trusting-a-result) #2.
- **It does not compare two sites.** One site, one variable, one
  analysis. There is no control-site design here.
- **It does not tell you whether a value is acceptable.** It compares
  against your own history, not against a guideline or standard. A site
  that has always been polluted reads as "normal" here.
- **It does not handle multiple variables jointly.** Run it once per
  variable; if you run it on ten variables and report the best p-value,
  you have re-created the multiple-comparisons problem the tool works
  hard to avoid internally.

## 2. How to use this guide

You don't need to read this top to bottom.

| If you want to… | Go to |
|---|---|
| Just run something | [Section 6](#6-running-each-mode), or `QUICKSTART.md` |
| Decide between the four modes | [Section 5](#5-which-mode-do-i-use) |
| Understand a number you got | [Section 7](#7-reading-the-output-every-number-in-order) |
| Know how much to trust it | [Section 8](#8-calibration--how-much-can-you-trust-the-p-value) — don't skip |
| Decide whether to report it | [Sections 9](#9-limitations--read-this-before-trusting-a-result) and [10](#10-a-checklist-before-you-report-a-result) |
| Fix an error message | [Section 11](#11-troubleshooting) |
| See the statistics | [Section 12](#12-additional-info) |
| Look up a term | [Section 13](#13-glossary) |

## 3. Setup

```
pip install -e .
```

from the repository root. That installs the package and its four
dependencies (numpy, pandas, scipy, matplotlib), so `import tide_lite`
works from any directory. No account, no API key, nothing else to
configure.

If you'd rather not install anything, the dependencies alone are enough
as long as you run from the repository root:

```
pip install numpy pandas scipy matplotlib
```

Check it works:

```
python examples/run_example.py     # full walkthrough, writes four plots
python -m pytest tests/ -q         # 38 checks, ~80 seconds
```

## 4. Preparing your data

You need two things, each as a table with a date column and a value
column (a pandas DataFrame, or anything you can load into one — CSV,
Excel via pandas, a database query result):

- **Historical data**: your "normal" years. More is better — see
  [Calibration](#8-calibration--how-much-can-you-trust-the-p-value).
  Ten years is a reasonable practical minimum; three is the hard floor
  the tool enforces.
- **Treatment data**: the period (or periods) you actually want an
  answer about.

Practical notes:

- **One variable at a time.** One column, one analysis.
- **Missing days are fine, up to a point.** The tool bins data (weekly
  or monthly, your choice). By default it drops historical years missing
  more than half their bins, and refuses a treatment period missing more
  than half of its bins. Individual missing bins in a treatment period
  are excluded from the test and reported — they are never scored.
- **You do not need to remove seasonality yourself.** That's the tool's
  job. Give it raw dated values.
- **Dates can be strings or datetimes.** Both work.
- **Gaps between years are fine.** A record of 2000–2004 and 2010–2014
  is handled correctly; the trend axis uses calendar years, not
  positions in a list.
- **Feb 29 is dropped** so every year has a 1–365 axis. At most one day
  every four years.
- If your variable can be legitimately zero or negative (a stage
  relative to a datum, an anomaly, redox potential), read the note on
  `detrend_mode` in [Section 6](#6-running-each-mode) before you start.

### Do you have enough history?

| Historical years | What's realistic |
|---|---|
| < 3 | The tool refuses to fit. |
| 3–9 | Runs, but a single p-value carries little weight. Exploratory only. |
| 10–15 | Workable for Standard / Sequential / Cumulative. Always pair with power and sensitivity checks. |
| 20+ | Comfortable, and the minimum for `MonitoringSeries`. |
| 40+ | The test becomes conservative — a significant result is strong. |

## 5. Which mode do I use?

Four modes, four different questions. Answer these in order:

**Are you asking about more than one distinct thing that happened?**

- **No, just one period, once** → **Standard**. *"Did this period behave
  abnormally?"*
- **Yes, several separate events, and you want to know which were real**
  → **Sequential**. *"Across these several things that happened, which
  were real effects and which could just be noise?"*
- **Yes, but it's really one known event you're following over time**
  → **Cumulative**. *"How is this one known event progressing? Has it
  recovered?"*
- **There's no event at all — you just want to check each new year as it
  arrives, indefinitely** → **MonitoringSeries**. *"Is this year tracking
  the same as history, checked year after year, without the false-alarm
  risk piling up the longer I keep watching?"*

| | Standard | Sequential | Cumulative | MonitoringSeries |
|---|---|---|---|---|
| Question | Is this one period unusual? | Which of these separate events were real? | Is this one event still showing an effect? | Is each new year, checked as it arrives, still normal? |
| How many things tested | One | Several, independent | One, tracked across follow-up windows | Open-ended, one at a time |
| Corrects for multiple testing? | N/A (one test) | Yes (Holm-Bonferroni) | No — deliberately, see below | Yes (fixed budget split across a committed horizon) |
| Needs all the data up front? | Yes | Yes | Yes | **No** — that's the point |
| Minimum historical years | ~10 | ~10 | ~10 | **~20** |
| Typical use | A single before/after check | Several sites, or several separate incidents | Post-restoration recovery tracking | Routine ongoing monitoring |

**Why doesn't Cumulative correct for multiple testing, if it runs several
tests too?** Sequential asks "is each of these a new discovery" — test
enough things carelessly and you'll rack up false alarms, so it corrects.
Cumulative asks "how is this **one** already-established event doing over
time" — each check-in isn't a fresh discovery claim, and correcting it
the same way would cost you the power to actually watch recovery happen.

**A borderline case:** if you have one underlying event but want each
year's check to stand alone as its own strict "did this happen" claim
rather than as a trajectory, use Sequential even though it's one event.
The deciding factor is the question you're asking, not how many things
caused it.

**MonitoringSeries vs. Sequential:** both correct for multiple testing,
but Sequential needs the whole batch up front (it ranks every p-value
together to set thresholds), so it can't be revealed one year at a time.
MonitoringSeries is built for exactly that: you don't need to know today
how many years you'll ultimately check.

## 6. Running each mode

All four start the same way: fit the historical data once.

```python
from tide_lite import TideConfig, fit_historical

config = TideConfig(
    bin_days=30,              # ~monthly bins (use 7 for weekly)
    detrend_mode="additive",  # "log_additive" if spread scales with level
    n_bootstrap=2000,         # Monte Carlo draws; see "choosing n_bootstrap"
)
fit = fit_historical(historical_df, date_col="date", value_col="value", config=config)

print(fit.years_used)                 # which years made it in
print(fit.years_dropped)              # and which didn't, with the reason
print(fit.trend)                      # was a long-term drift detected and removed?
print(fit.between_year_var_frac)      # see Section 8 — check this every time
```

Misspelled options raise immediately rather than silently running a
different analysis, so `detrend_mode="log"` is an error, not a surprise.

### Choosing `n_bootstrap`

The smallest p-value the test can ever produce is `1 / (n_bootstrap + 1)`
(`config.min_attainable_p`). At the default 2000 that's 0.0005, which is
ample for a Standard test at alpha=0.05. It starts to matter when you
push thresholds down:

- `bin_significant()` across 13 bins needs p as small as `alpha / 13`
  (0.0038 at alpha=0.05) — fine at 2000 draws.
- `MonitoringSeries` needs p as small as `alpha_total / n_years_horizon`
  (0.0025 for a 20-year horizon at alpha_total=0.05).

A useful rule: **`n_bootstrap` should be at least 10 / your smallest
threshold**, so a flag rests on a reasonable number of draws rather than
the two or three most extreme. For a 20-year monitoring horizon that
means `n_bootstrap=4000`.

Both functions warn when `n_bootstrap` is marginal for the threshold
being asked of it, and `MonitoringSeries` refuses outright when the
threshold is unreachable — rather than quietly returning "nothing found"
forever. Raising `n_bootstrap` costs run time and nothing else.

### Reproducibility — pass `rng=`

Every p-value here comes from a Monte Carlo null, so two runs on the same
data differ slightly. **For anything you intend to report, pass a seed**:

```python
result = test_treatment(fit, treatment_df, rng=42)
seq    = sequential_test(fit, events, rng=42)
cum    = cumulative_test(fit, windows, rng=42)
grid   = sensitivity_grid(..., rng=42)
```

Without it a borderline result can land on either side of alpha between
runs — and `first_recovery()` keys off exactly which windows cleared
alpha.

### `year_index` — the one parameter people get wrong

Whenever a trend was detected and removed (`fit.trend["applied"]` is
True), the tool has to know **when** your treatment period happened, so
it can extrapolate the historical trend to that point for a fair
comparison. That's `year_index`: the offset in **calendar years** from
your first historical year.

**Use `fit.year_index_for(...)` and don't do the arithmetic yourself:**

```python
result = test_treatment(fit, df_2024, treatment_year_index=fit.year_index_for(2024))
```

`len(fit.years_used)` is a common stand-in and is only the same number
when the record has no missing years. With history running 2000–2004 and
2010–2014, `len(fit.years_used)` is 10 but a 2015 treatment period is 15
years out.

Modes that test several periods (Sequential, Cumulative,
MonitoringSeries) **require** an explicit index when a trend is active —
they can't assume their periods are adjacent in time.

### Standard

```python
from tide_lite import test_treatment

result = test_treatment(fit, treatment_df, mode="prediction", rng=42)
print(result.p_value, result.effect_size)
```

`mode="prediction"` (the default) tests exactly one year and raises if
`treatment_df` spans more. `mode="confidence"` averages several years
into one curve and tests that average — use it only when the average
really is the thing you want to test.

### Sequential

```python
from tide_lite import SequentialEvent, sequential_test

events = [
    SequentialEvent(label="site A", treatment_df=site_a_df, year_index=fit.year_index_for(2021)),
    SequentialEvent(label="site B", treatment_df=site_b_df, year_index=fit.year_index_for(2021)),
]
results = sequential_test(fit, events, alpha=0.05, rng=42)
for r in results:
    print(r.label, r.result.p_value, "corrected significant:", r.significant_corrected)
```

Read `significant_corrected`, not `significant_raw` — the corrected
column is the whole point of using this mode.

**Your responsibility:** the date ranges behind the events must not
overlap. Overlapping data silently breaks the independence the correction
assumes, and the function can't check this because it only receives
already-separated per-event frames.

### Cumulative

```python
from tide_lite import CumulativeWindow, cumulative_test, first_recovery

windows = [
    CumulativeWindow("year 1 after", y1_df, year_index=fit.year_index_for(2022)),
    CumulativeWindow("year 2 after", y2_df, year_index=fit.year_index_for(2023)),
    CumulativeWindow("year 3 after", y3_df, year_index=fit.year_index_for(2024)),
]
results = cumulative_test(fit, windows, alpha=0.05, rng=42)

recovery = first_recovery(results, consecutive_required=2)
print(recovery.recovered, recovery.recovered_at_label)
```

`first_recovery` requires `consecutive_required` windows in a row to be
non-significant before declaring recovery — a lone clear check-in can be
noise. Use 2 or more for anything you'd report. A later false alarm after
recovery is locked in doesn't change the answer.

### MonitoringSeries

Unlike the other three, this one is stateful: create it once, then call
`.check()` each time a new year's data is ready — in practice a separate
script run, once a year.

```python
from tide_lite import MonitoringSeries

# n_bootstrap=4000: a 20-year horizon splits alpha down to 0.0025, and the
# rule above wants at least 10 / 0.0025 draws behind that threshold.
fit    = fit_historical(historical_df, "date", "value", TideConfig(n_bootstrap=4000))
series = MonitoringSeries(fit, alpha_total=0.05, n_years_horizon=20)

# this year:
check = series.check(this_years_df, label="2024",
                     year_index=fit.year_index_for(2024), rng=42)
print(check.p_value, check.flagged)
series.save("monitoring_state.json")

# next year, in a new script run:
fit = fit_historical(historical_df, date_col="date", value_col="value")  # re-fit
series = MonitoringSeries.load(fit, "monitoring_state.json")
check = series.check(next_years_df, label="2025", year_index=fit.year_index_for(2025))
```

`alpha_total` is your **total budget across the whole horizon**, not per
year — each check uses `alpha_total / n_years_horizon`, so a longer
horizon buys open-endedness at the cost of per-year sensitivity. Commit
to the shortest horizon you'll actually need. When the horizon
is used up, `.check()` raises; call `series.renew()` to start a fresh one
(deliberately not automatic — extending the guarantee is a real
decision).

Saved state holds the budget and the check log, not the fit itself:
re-run `fit_historical` on the same historical data each session.

**This mode needs ~20+ historical years.** It warns below that, and
[Limitations](#9-limitations--read-this-before-trusting-a-result) #6
explains why.

### Sustained departures and recovery within one period

Three tools for a question the whole-period p-value can't answer: was the
departure sustained (not just one noisy month), and did it come back?

**Sustained-departure test.** Asks "was there a run of several
consecutive bins that were *together* unusual" — more robust to one noisy
month, more sensitive to a persistent moderate shift no single month
flags alone:

```python
result = test_treatment(fit, treatment_df, run_lengths=[3], rng=42)
print(result.sustained[3].p_value, result.sustained[3].window_bins)
```

A genuine trade-off, not a strictly better test: a sharp one-month spike
is easier for the single-bin test to catch, since averaging over a window
dilutes it. Use both if you don't know which shape to expect.

Note the statistic averages the *absolute* deviation per bin, so a period
that runs high for two months and then low for two months scores like one
held high for four. Check the plot and `effect_size` for direction.

**Where did it start and end.**

```python
from tide_lite import estimate_departure_recovery

dep = estimate_departure_recovery(result, run_length=3)
print(dep.departure_start_bin, dep.departure_end_bin, dep.recovered, dep.recovered_at_bin)
```

This is **descriptive, not a hypothesis test**. The sustained test already
established significance for the core window; this traces outward from it
using a looser, uncorrected criterion, because "roughly how long did this
last" is a different question from "prove each additional month
independently". Expect it to over-reach by a bin or so.

### Plotting

```python
from tide_lite import plot_envelope, plot_sequential, plot_cumulative, plot_monitoring

plot_envelope(fit, result, treatment_label="2024")   # Standard
plot_sequential(fit, seq_results)                    # Sequential
plot_cumulative(cum_results)                         # Cumulative
plot_monitoring(series)                              # MonitoringSeries
```

- **`plot_envelope`** — historical band and individual years, treatment
  period highlighted, a star on each bin that's individually significant
  after correction, and (pass `departure=`) a shaded departure span with
  a recovery line. When a trend was removed, the envelope is
  automatically shifted onto the treatment period's own year so the
  picture matches the p-value.
- **`plot_sequential`** — every event on one envelope, coloured by
  whether it survived correction.
- **`plot_cumulative`** — trend-adjusted effect size across follow-ups,
  coloured by significance, dashed line at zero ("back to normal"). A
  declining, colour-changing line reads as a recovery curve.
- **`plot_monitoring`** — every check against the flat per-check bar it
  had to clear, log scale.

All four return a matplotlib `Axes` — call `.figure.savefig(...)`
yourself, or pass `ax=` to place one in a larger figure.

## 7. Reading the output: every number, in order

**`p_value`** — "if this period were really no different from normal, how
surprising would data like this be?" Small (conventionally < 0.05) means
surprising. It is **not** the probability that there's no effect, and
**not** a measure of how big the effect is.

**`effect_size`** — how different the period looked from the historical
typical level, in your data's own units, signed (positive = higher).
Measured against the historical record's own typical level, ignoring any
trend.

**`effect_size_trend_adjusted`** — the same, with the fitted historical
trend extrapolated to the period's own year and subtracted first. **This
is the one the p-value actually tested.** The two are identical when no
trend was detected. When they differ a lot, the difference *is* the trend
— quote the adjusted one alongside the p-value.

**`test_statistic`** — the largest standardized deviation in any single
bin. Useful for comparing periods; not interpretable on its own.

**`bin_p_values`** — per-bin, **uncorrected**. Do not compare these to
0.05 directly (see below). `NaN` marks a bin with no treatment data.

**`missing_bins`** — boolean per bin: no treatment data there. Those bins
take no part in the test.

**`n_reference`** — Monte Carlo draws behind the null. The smallest
possible p-value is `1 / (n_reference + 1)`; a p-value sitting exactly at
that floor means "at least this extreme", not a precise value.

**`treatment_year_index`** — the resolved calendar-year offset the period
was tested at. Worth checking when a trend is active.

### Which months, specifically — `bin_significant()`

```python
from tide_lite import bin_significant

sig = bin_significant(result, alpha=0.05)   # bool array, one per bin
print(fit.bins[sig])
```

This is deliberately **not** the same as comparing each bin's raw p-value
to 0.05. With 13 independent bins, that would flag something in about
**half** of all perfectly normal years (1 − 0.95¹³ ≈ 0.49) purely from
running 13 tests instead of one. `bin_significant()` applies
Holm-Bonferroni across the bins, so the false-alarm rate stays
controlled. Bins with no data are never flagged and don't make the other
bins' thresholds stricter.

It answers a genuinely different question from the global p-value — a
period can be significant overall with no single bin surviving this
stricter correction, and (less often) vice versa. Check both.

### The plot

The shaded band is where historical years typically fall; thin lines are
the actual historical years, so you can see real spread rather than a
summary. A bold treatment line clearly outside the band for a sustained
stretch is what a genuinely significant result looks like. If the p-value
is small but the line barely strays, look again — and check
`between_year_var_frac` (next section).

### Power curve and sensitivity grid

**`estimate_power`** — "given how much history I have, how reliably could
this test detect an effect of a given size?" A flat, low curve at sizes
you'd care about means a non-significant result may just reflect an
underpowered test rather than "no effect". Note the numbers are the
*optimistic* end: the simulated periods are drawn from the same pool that
builds the null, so real-world power is somewhat lower. Power at effect 0
comes back at ≈ alpha by construction — that's arithmetic, not evidence
of calibration.

**`sensitivity_grid`** — does the result hold up under small, defensible
changes to bin size or detrending? A result that flips from significant
to not under a minor reasonable change is much weaker than one stable
across the grid. Watch `n_historical_years` per row: a variant that
quietly drops years isn't a fair comparison.

## 8. Calibration — how much can you trust the p-value?

This section exists because the honest answer is "it depends on your
data, and the tool can now tell you which case you're in."

### The measurement

Feeding the pipeline fresh synthetic "normal" years and counting how
often it wrongly flags them, at a nominal alpha of 0.05 (900 trials per
cell):

| Historical years | Within-year noise only | With whole-year level shifts |
|---|---|---|
| 10 | 0.054 | **0.112 – 0.130** |
| 20 | 0.042 | **0.072 – 0.081** |
| 40 | 0.018 | — |

### What this means

When year-to-year variation is just accumulated within-year wiggle, the
test is well calibrated, and gets conservative on long records.

When your data **also** carries genuine whole-year level shifts — a wet
year sitting high from January to December, which is the norm in
hydrology — the test runs at roughly **twice** its nominal false-positive
rate with 10 historical years, and about 1.5× with 20.

### Why

The test builds its "what does a normal year look like" reference by
stitching each synthetic year from about four independently drawn
sub-year chunks of your real history. That handles within-year
autocorrelation well. But a whole-year offset present in every real year
gets averaged away in the stitching: measured directly, the spread of
annual means across synthetic years is about **half** that of the real
historical years once a year-level effect exists. The reference is too
narrow, so real years look more extreme than they should.

This is a real property of the method rather than a bug. Resampling whole
years instead is worse — with N historical years the smallest possible
p-value becomes 1/(N+1), which at N=10 is 0.09, above alpha, making
significance unreachable at any effect size. Adding a resampled
whole-year offset back onto each draw was prototyped and did **not** fix
it: with N historical years that offset distribution is itself capped at
the most extreme year ever observed, while a genuinely new year exceeds
all N of them about 2/(N+1) of the time.

### What to do about it

`fit_historical` reports the diagnostic that decides which row of the
table applies to you:

```python
print(fit.between_year_var_frac)    # 0.0 - 1.0
```

| Value | Reading |
|---|---|
| **< 0.3** | Whole-year shifts are minor. Top column of the table — well calibrated. |
| **0.3 – 0.5** | Mixed. Treat p-values within ~2× of alpha with care. |
| **> 0.5** | Whole-year shifts dominate. Bottom column. With fewer than 20 historical years you'll get a warning, and a p-value of 0.03 should be read as "suggestive", not "significant". |

If you're in the high regime: get more historical years if you can, lean
on `effect_size` and the plot rather than the p-value alone, require a
sustained departure rather than a single bin, and treat anything within
2× of alpha as inconclusive.

## 9. Limitations — read this before trusting a result

These aren't small print. Each one changes how much weight a result can
reasonably carry.

**1. With fewer than ~15–20 historical years, one p-value is noisier than
it looks.** Beyond the calibration issue above, the false-positive rate
*conditional on one fitted model* varies a lot fit to fit — in a 25-fit
simulation it ranged from 0% to 25% at nominal alpha=0.05, purely from
the sampling noise of estimating a median/MAD curve from ~10 years. It
averages out correctly across many possible historical samples, but your
real analysis only ever has the one. **What to do:** run `estimate_power`
and `sensitivity_grid` alongside any real result, especially near 0.05.

**2. "Unusual" is not "caused by the treatment."** Without a control
site, a significant result says the period was statistically unusual
relative to history — it doesn't rule out a coincident cause (an unusually
wet or dry period, an upstream event). **What to do:** check at least one
obvious alternative (was rainfall or flow also anomalous?) before
attributing it.

**3. The method assumes your historical years are genuinely comparable.**
A real slow change that isn't a roughly linear drift (land-use change, a
step change in method or lab) won't be fully captured by the built-in
trend check. **What to do:** look at the historical years overlaid — they
should read as variations on one process, not a visibly drifting or
splitting pattern.

**4. This is a new combination of established methods, not itself a
published, peer-reviewed method.** The pieces — block bootstrap,
permutation testing, Holm-Bonferroni — are well established. This
assembly of them isn't published. **What to do:** for internal research
and monitoring this is a non-issue; cite the underlying methods. If a
result needs to survive external legal or regulatory scrutiny, treat that
as a separate requirement this tool doesn't address.

**5. `detrend_mode="log_additive"` requires strictly positive data.** If
your variable can be zero or negative, use `"additive"` or offset the
data first. The tool raises a clear error rather than silently producing
NaNs.

**6. `MonitoringSeries` needs more historical years than the other three
modes.** Splitting an already-small budget across many checks pushes each
one deep into the tail of the bootstrap distribution, which is unreliable
when built from limited data. In simulation, a nominal 5% lifetime budget
over a 10-year horizon had a *true* false-alarm rate of roughly 17% with
10 historical years, 7% with 20, 3% with 30. **What to do:** use 20+
years; treat anything less as indicative, not a 5% guarantee. (A tool to
check this against your own data was attempted and deliberately left out
— see `tide_lite/monitoring.py`'s docstring for why the first design
would have given false reassurance exactly where a warning was needed.)

**7. The last bin of the year is short.** At `bin_days=30` the year
splits into twelve 30-day bins plus a 5-day remainder; at `bin_days=7`,
52 weeks plus one day. That final stub bin is noisier than the others.
It's handled (each bin is standardized by its own MAD) but worth knowing
if bin 13 looks odd.

**8. `mode="confidence"` averages bins across years using whatever data
is present.** If one of three treatment years is missing June, June's
value is a one-year average while the reference expects a three-year one.
Prefer complete years in confidence mode.

## 10. A checklist before you report a result

1. **Did you pass `rng=`?** If not, the number isn't reproducible.
2. **Check `fit.years_used` and `fit.years_dropped`.** Did a year you
   expected get silently filtered out?
3. **Check `fit.between_year_var_frac`** ([Section 8](#8-calibration--how-much-can-you-trust-the-p-value)).
   Which calibration regime are you in?
4. **Check `fit.trend`.** If a trend was applied, quote
   `effect_size_trend_adjusted`, and confirm your `year_index` values.
5. **Look at the plot.** Does the picture match the p-value? If not,
   something is off — investigate before reporting.
6. **Run `estimate_power`** at an effect size you'd care about. If power
   is low, a non-significant result means "couldn't tell", not "no
   effect".
7. **Run `sensitivity_grid`.** Does the conclusion survive reasonable
   alternative settings?
8. **Ask what else happened that period.** Rainfall, flow, upstream
   works, a change of lab or method.
9. **Report the effect size and the interval of your data, not just the
   p-value.**

## 11. Troubleshooting

**`Only N historical years survived the completeness filter (need >= 3)`**
— Fewer than 3 historical years had enough data (by default at least half
their bins filled). Provide more or more-complete history, or loosen
`max_missing_frac` if you accept noisier per-year estimates.

**`The treatment period is missing N of M bins (X%), above max_missing_frac`**
— Your treatment period is too sparse to compare against a whole-year
reference. Supply more complete data, use coarser bins (larger
`bin_days`), or raise `max_missing_frac` deliberately.

**`The treatment period has no data in any of the historical fit's N bins`**
— Usually a column-name or date-parsing mismatch, or a treatment period
covering a different part of the year than the history.

**`UserWarning: The treatment period has no data in N of M bins`** — Not
an error. Those bins are excluded and their `bin_p_values` are `NaN`; the
p-value uses the bins that do have data.

**`detrend_mode='log_additive' requires strictly positive values`** — Zero
or negative values present. Switch to `"additive"` or offset the data.

**`detrend_mode must be one of [...]` / `agg must be one of [...]` /
`mode must be one of [...]`** — A misspelled option. These used to fall
through silently to a different analysis; now they raise.

**`mode='prediction' expects exactly one year of treatment data, but
treatment_df spans N years`** — Restrict to one year, or use
`mode="confidence"` if you genuinely want the multi-year average.

**`every SequentialEvent needs an explicit year_index`** (or the
Cumulative / MonitoringSeries equivalent) — A trend was detected and at
least one period didn't say when it happened. Use
`fit.year_index_for(2024)`.

**`This horizon cannot flag anything`** — `alpha_total / n_years_horizon`
is below the smallest p-value your `n_bootstrap` can produce. Raise
`n_bootstrap` (the message tells you to what), shorten the horizon, or
raise `alpha_total`.

**`UserWarning: No bin can be flagged at alpha=...`** — Same arithmetic,
for `bin_significant()`. The all-False result is forced by resolution,
not by your data. Raise `n_bootstrap`.

**`UserWarning: N% of this record's residual variation is whole-year
level shifts`** — You're in the anti-conservative regime with a short
record. See [Section 8](#8-calibration--how-much-can-you-trust-the-p-value).

**`Monitoring horizon of N checks is used up`** — Call `.renew()` to start
a fresh horizon. Deliberately not automatic.

**`consecutive_required must be >= 1`** — `0` made the recovery rule
vacuous (an empty run trivially satisfies "all non-significant"), so it
declared recovery at the first window regardless.

**A result seems "too significant" given how the plot looks** — Check
`effect_size` against `effect_size_trend_adjusted`. If a trend was
applied, the p-value is computed on the adjusted basis.

**Two runs give different p-values** — Expected without `rng=`. Pass a
seed.

## 12. Additional info

### Other documents in this project

| File | What it's for |
|---|---|
| `README.md` | One-page overview and project layout. |
| `QUICKSTART.md` | Get a first result in ~5 minutes. |
| `USER_GUIDE.md` | This file — modes, interpretation, calibration, limitations. |
| `examples/run_example.py` | A runnable tutorial: every mode end to end on synthetic data. Read it top to bottom. |
| `tests/test_engine.py` | Statistical validation of the method (calibration, power, FWER). |
| `tests/test_monitoring.py` | `MonitoringSeries` validation. |
| `tests/test_regressions.py` | One test per implementation bug found in the code audit, named for the wrong answer it prevents. |
| `tide_lite/engine.py` | The core, with the reasoning for every design choice in comments — including the calibration section behind [Section 8](#8-calibration--how-much-can-you-trust-the-p-value). |

### The statistics, briefly

- **Block bootstrap** (Künsch 1989; Politis & Romano 1992/1994): instead
  of assuming a distribution, this resamples chunks of your own
  historical data to build a picture of plausible normal variation.
  Chunks are matched to the calendar position they came from (so a
  volatile month's pattern never gets tested against a stable month's
  spread) and pulled from a small neighbourhood of nearby positions for
  enough variety to resolve a precise p-value. See
  [Section 8](#8-calibration--how-much-can-you-trust-the-p-value) for
  what this handles and what it doesn't.
- **Median and MAD** (median absolute deviation, 1.4826×-corrected)
  instead of mean and standard deviation: robust to the skew common in
  concentration data, without assuming a bell curve.
- **Mann-Kendall test and Sen's slope** (Mann 1945; Kendall 1975; Sen
  1968): detect and estimate a long-term trend in the historical data
  before testing, so a real ongoing drift isn't mistaken for a treatment
  effect. The slope is computed on a calendar-year axis, so records with
  gaps are handled correctly. Note the Mann-Kendall normal approximation
  is weak for very short records (< ~10 years).
- **Holm-Bonferroni correction** (Holm 1979): controls the chance of a
  false alarm across several tests at once, without being as conservative
  as plain Bonferroni. Used across events in Sequential and across bins
  in `bin_significant()`.
- **Max-deviation test statistic:** the whole-period p-value is based on
  the single largest standardized bin deviation, which is what lets one
  number cover a whole year without a per-bin multiple-comparisons
  problem.

## 13. Glossary

- **p-value** — how surprising your result would be if there were truly
  no effect. Small = surprising = evidence of a real difference.
- **Effect size** — how big the difference was, in your data's own units
  (separate from how *sure* you are it's real).
- **Envelope** — the band of historical "normal" variation a treatment
  period is compared against.
- **Bin** — a fixed slice of the year (30 days by default) that values
  are aggregated into, so years can be compared like with like.
- **Seasonality** — the regular within-year pattern this tool accounts
  for automatically.
- **Autocorrelation** — the tendency for a value on one day to resemble
  nearby days; ignoring it makes tests overconfident.
- **Between-year variance share** (`between_year_var_frac`) — how much of
  the leftover variation is whole-year level shifts rather than
  within-year wiggle. Drives calibration; see
  [Section 8](#8-calibration--how-much-can-you-trust-the-p-value).
- **Bootstrap / resampling** — building a picture of plausible variation
  by repeatedly drawing from your own data, instead of assuming a
  textbook distribution.
- **Detrending** — removing a genuine long-term drift from the historical
  data before testing.
- **Family-wise error rate (FWER)** — the chance of at least one false
  alarm across a whole set of tests.
- **Exchangeability** — the assumption that your historical years are
  comparable draws from "the same normal".
- **Pointwise vs. corrected (bin-level)** — pointwise compares each
  month's raw p-value to 0.05 on its own, racking up false alarms;
  corrected (`bin_significant()`) controls for that.
- **Anti-conservative** — a test that flags more often than its stated
  alpha claims, i.e. more false alarms than advertised.
