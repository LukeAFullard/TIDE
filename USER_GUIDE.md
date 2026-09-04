# TIDE lite — User Guide

A complete guide to using this tool, written so it's useful whether you're
running the code yourself or reading someone else's results. If you just
want to get started immediately, see `QUICKSTART.md` instead — this guide
is the comprehensive version.

## Contents

1. [What is this?](#1-what-is-this)
2. [How to use this guide](#2-how-to-use-this-guide)
3. [Setup](#3-setup)
4. [Preparing your data](#4-preparing-your-data)
5. [Which mode do I use?](#5-which-mode-do-i-use)
6. [Running each mode](#6-running-each-mode)
7. [Interpreting your results](#7-interpreting-your-results)
8. [Limitations — read this before trusting a result](#8-limitations--read-this-before-trusting-a-result)
9. [Troubleshooting](#9-troubleshooting)
10. [Additional info](#10-additional-info)
11. [Glossary](#11-glossary)

---

## 1. What is this?

You have water quality data covering many "normal" years. Something
happened — a restoration project, a treatment plant upgrade, a policy
change, a pollution event — during or after a specific period. You want
to know: **did that period actually look different from what's normal,
once you account for the fact that water quality naturally moves with the
seasons and from year to year?**

This tool answers that question with a statistical test built for exactly
this situation:

- It doesn't need a "control site" to compare against — many real
  monitoring situations only have one site.
- It respects seasonality automatically — it doesn't compare a July value
  to a January baseline and call that unusual.
- It doesn't assume your data follows a bell curve — real concentration
  data is often skewed, and this test doesn't care.

The output is a p-value (how surprising the observed period would be if
nothing had really changed), an effect size (how different, in your
data's own units), and a plot showing the period against the historical
range.

## 2. How to use this guide

You don't need to read this top to bottom.

- **Just want to run something?** → [Section 6](#6-running-each-mode), or
  `QUICKSTART.md`.
- **Deciding between the three modes?** → [Section 5](#5-which-mode-do-i-use).
- **Have a result and want to know what it means?** → [Section 7](#7-interpreting-your-results).
- **Deciding whether to trust or report a result?** → [Section 8](#8-limitations--read-this-before-trusting-a-result) — don't skip this one.
- **Hit an error message?** → [Section 9](#9-troubleshooting).
- **Want the statistics behind it, or the project's own build history?** → [Section 10](#10-additional-info).
- **A term doesn't make sense?** → [Section 11](#11-glossary).

## 3. Setup

```
pip install numpy pandas scipy matplotlib
```

That's the complete list of dependencies. No account, no API key, nothing
else to configure.

## 4. Preparing your data

You need two things, each as a table with a date column and a value
column (a pandas DataFrame, or anything you can load into one — CSV,
Excel via pandas, a database query result):

- **Historical data**: your "normal" years. More is better — see
  [Limitations](#8-limitations--read-this-before-trusting-a-result) for
  why. Ten years is a reasonable practical minimum; three is the hard
  floor the tool enforces.
- **Treatment data**: the period (or periods) you actually want an answer
  about.

Practical notes:

- One variable at a time (one column, one analysis). If you're testing
  several variables, run this once per variable.
- Missing days are fine up to a point — the tool bins data (weekly or
  monthly, your choice) and tolerates gaps, dropping only historical
  years missing more than half their bins by default.
- You do **not** need to remove seasonality yourself — that's what the
  tool does. Just give it raw dated values.
- Dates can be strings or actual datetime values; both work.
- If your variable can be legitimately zero or negative (a stage relative
  to a datum, an anomaly measurement, redox potential), read the note on
  `detrend_mode` in [Section 6](#6-running-each-mode) before you start.

## 5. Which mode do I use?

Three modes, three different questions. Answer these in order:

**Are you asking about more than one distinct thing that happened?**

- **No, just one period, once** → use **Standard**. *"Did this period
  behave abnormally?"*
- **Yes, several separate events, and you want to know which (if any)
  were real** → use **Sequential**. *"Across these several things that
  happened, which were real effects and which could just be noise?"*
- **Yes, but it's really one known event you're following over time**
  (is it still showing an effect? has it recovered?) → use
  **Cumulative**. *"How is this one known event progressing?"*
- **There's no event at all — you just want to check each new year as it
  arrives, indefinitely** (routine monitoring, not testing a hypothesis)
  → use **MonitoringSeries**. *"Is this year tracking the same as
  history, checked year after year, without the false-alarm risk piling
  up the longer you keep watching?"*

| | Standard | Sequential | Cumulative | MonitoringSeries |
|---|---|---|---|---|
| Question | Is this one period unusual? | Which of these separate events were real? | Is this one event still showing an effect? | Is each new year, checked as it arrives, still normal? |
| How many things tested | One | Several, independent | One, tracked across several follow-up windows | Open-ended — one at a time, indefinitely |
| Corrects for multiple testing? | N/A (only one test) | Yes (Holm-Bonferroni) | No — deliberately, see below | Yes (a fixed budget split evenly across a committed horizon) |
| Typical use | A single before/after check | Comparing several sites, or several separate incidents | Post-restoration recovery tracking | Routine, ongoing monitoring with no specific "treatment" |

**Why doesn't Cumulative correct for multiple testing, if it's running
several tests too?** Sequential is answering "is each of these a new
discovery" — treat that carelessly and you'll rack up false alarms just
by testing enough things, so it corrects for that. Cumulative is
answering "how is this **one** known, already-established event doing
over time" — each check-in isn't a fresh discovery claim, and correcting
it the same way would cost you the power to actually see recovery happen.

**A borderline case:** if you have one underlying event but want each
year's check to stand alone as its own strict "did this happen"
claim rather than a trajectory, use Sequential even though it's a single
event — the deciding factor is the question you're asking, not how many
things caused it.

**MonitoringSeries vs. Sequential:** both correct for multiple testing,
but Sequential needs the whole batch of events up front (it ranks every
p-value together to set its thresholds) — it can't be revealed one year
at a time as new data arrives. MonitoringSeries is built for exactly
that: you don't need to know today how many years you'll ultimately
check.

## 6. Running each mode

All four start the same way: fit the historical data once.

```python
from tide_lite import TideConfig, fit_historical

config = TideConfig(
    bin_days=30,              # ~monthly bins (use 7 for weekly)
    detrend_mode="additive",  # use "log_additive" if your variable's
                               # spread scales with its level (common for
                               # concentrations) -- requires all-positive
                               # data, see Troubleshooting if it errors
)
fit = fit_historical(historical_df, date_col="date", value_col="value", config=config)
```

**A note on "normalization," if you know the term from elsewhere:** it
means at least three different things people sometimes conflate --
worth being precise. Rescaling relative to a defined reference/baseline
period isn't done here at all (deliberately left out — one fewer
place for a baseline window to accidentally overlap the period you're
testing, a real risk if you build that kind of feature). Removing a
long-term drift across your historical years IS done (`detrend_mode`)
and is a different operation, about time, not a reference level.
Log-transforming skewed data to stabilize its variance is ALSO done,
via `detrend_mode="log_additive"` — the one that actually matches what
"normalization" usually means statistically — bundled into the same
setting as detrending because in practice you want both done in the
same space, not because they're the same operation. And every
historical year always gets expressed as a deviation from the typical
seasonal shape, which you could loosely call normalizing too, but it's
the core mechanism the test runs on, not a setting you choose.

### Standard

```python
from tide_lite import test_treatment

result = test_treatment(fit, treatment_df, mode="prediction")
print(result.p_value, result.effect_size)
```

### Sequential

```python
from tide_lite import SequentialEvent, sequential_test

events = [
    SequentialEvent(label="site A", treatment_df=site_a_df, year_index=10),
    SequentialEvent(label="site B", treatment_df=site_b_df, year_index=10),
]
results = sequential_test(fit, events, alpha=0.05)
for r in results:
    print(r.label, r.result.p_value, "corrected significant:", r.significant_corrected)
```

`year_index` is required whenever your historical data has a detectable
trend (check `fit.trend["applied"]`) — Sequential events aren't assumed
to happen back-to-back, so it can't guess this for you. It's the number
of years after your first historical year that this event occurred (see
[Section 9](#9-troubleshooting) if you get an error about it).

### Cumulative

```python
from tide_lite import CumulativeWindow, cumulative_test

windows = [
    CumulativeWindow("year 1 after", year1_df, year_index=10),
    CumulativeWindow("year 2 after", year2_df, year_index=11),
    CumulativeWindow("year 3 after", year3_df, year_index=12),
]
results = cumulative_test(fit, windows, alpha=0.05)
for r in results:
    print(r.label, r.result.p_value, r.result.effect_size, r.significant)
```

Same `year_index` requirement as Sequential, if a trend was detected.

### MonitoringSeries

Unlike the other three, this one is stateful: create it once, then call
`.check()` each time a new year's data is ready — in practice, that's
usually a separate script run, once a year, not all in one sitting.

```python
from tide_lite import MonitoringSeries

series = MonitoringSeries(fit, alpha_total=0.05, n_years_horizon=20)

# this year:
result = series.check(this_years_df, label="2024")
print(result.p_value, result.flagged)
series.save("monitoring_state.json")   # so next year picks up where this left off

# next year, in a new script run:
fit = fit_historical(historical_df, date_col="date", value_col="value")  # re-fit, same historical data
series = MonitoringSeries.load(fit, "monitoring_state.json")
result = series.check(next_years_df, label="2025")
```

`alpha_total` is your total budget across the whole `n_years_horizon`-year
commitment (not per year) — each check uses `alpha_total / n_years_horizon`.
When the horizon is used up, call `series.renew()` to start a fresh one.

**Before relying on this mode, read the note on historical data size in
[Limitations](#8-limitations--read-this-before-trusting-a-result)** — it
needs more historical years than the other three modes to work as
intended.

### Sustained departures and recovery

Three related tools for a question the whole-year p-value can't answer on
its own: was the departure sustained (not just one noisy month), and did
it come back to normal?

**Sustained-departure test.** `bin_significant()` (above) asks "was any
single bin unusual." This asks a different question: "was there a run of
several consecutive bins, *together*, unusual" — more robust to one noisy
month, more sensitive to a persistent but moderate shift that no single
month is extreme enough to flag alone. Pass `run_lengths` to
`test_treatment`:

```python
result = test_treatment(fit, treatment_df, mode="prediction", run_lengths=[3])
print(result.sustained[3].p_value, result.sustained[3].window_bins)
```

This is a genuine trade-off, not a strictly better test — a sharp,
one-month spike is easier for the *single*-bin test to catch, since
averaging over a window dilutes it with its normal neighbors. Use both
if you're not sure which shape to expect.

**Where did it start and end — `estimate_departure_recovery()`.** Given
a sustained result, traces outward from the significant core window,
bin by bin, while each additional bin still looks unusual, and reports
where (if anywhere within the period) it recovered:

```python
from tide_lite import estimate_departure_recovery

dep = estimate_departure_recovery(result, run_length=3)
print(dep.departure_start_bin, dep.departure_end_bin, dep.recovered, dep.recovered_at_bin)
```

This is descriptive, not a fresh hypothesis test in its own right — the
sustained test already established significance for the core; this
traces how far it extends using a looser, uncorrected criterion, since
requiring full correction at every candidate boundary would be a
stricter, different question ("prove each additional month
independently") than "roughly how long did this last." `plot_envelope`
shades this span automatically when you pass `departure=dep`.

**Recovery across years — `first_recovery()`.** Cumulative already
tracks effect size across follow-up windows; this adds a formal rule for
when to call it "recovered," rather than eyeballing the plot:

```python
from tide_lite import first_recovery

recovery = first_recovery(cum_results, consecutive_required=2)
print(recovery.recovered, recovery.recovered_at_label)
```

Requires `consecutive_required` windows in a row to be non-significant,
not just one — a lone clear check-in can be noise (this project's own
testing found real fits can show a false-clear on a single check a
meaningful fraction of the time). A later false alarm after recovery is
already locked in doesn't change the answer.

### Plotting — every mode has a dedicated function

Each answers that mode's actual question visually, not just the same
generic picture reused four times:

```python
from tide_lite import plot_envelope, plot_sequential, plot_cumulative, plot_monitoring

plot_envelope(fit, result, treatment_label="2024")          # Standard
plot_sequential(fit, seq_results)                           # Sequential
plot_cumulative(cum_results)                                # Cumulative
plot_monitoring(series)                                     # MonitoringSeries
```

- `plot_envelope` — historical band and individual years, with the
  treatment period highlighted, a star on any specific bin
  (month, if `bin_days=30`) that's individually significant after
  correction, and (pass `departure=`) a shaded span for an estimated
  sustained departure with a line marking recovery — straightforward "is
  this line outside the shaded area, specifically where, and did it come
  back."
- `plot_sequential` — every tested event on the same envelope, colored
  by whether it survived correction — flagged events pop out visually.
- `plot_cumulative` — effect size across successive follow-ups, colored
  by significance, with a dashed line at zero ("back to normal"). A
  declining, color-changing line reads as a recovery curve at a glance.
- `plot_monitoring` — every check so far against the flat per-check bar
  it had to clear, on a log scale (thresholds are often small).

All four return a matplotlib `Axes` — call `.figure.savefig(...)` or
`.figure.show()` yourself, or pass `ax=` to place one inside a larger
figure. See `examples/run_example.py` for all four running end to end.

## 7. Interpreting your results

**p-value.** Roughly: "if the treatment period were really no different
from normal, how surprising would data like this be?" A small p-value
(conventionally < 0.05) means very surprising — probably not just natural
variation. It is **not** the probability that there's no effect, and it's
**not** a measure of how big the effect is (see effect size for that).

**effect_size.** How different the treatment period looked from the
historical typical, in your data's own units, with the sign telling you
direction (positive = higher than typical).

**effect_size_trend_adjusted.** The same thing, but with any ongoing
historical trend subtracted first. If your historical data has a real
long-term drift (checked automatically), the p-value is computed on this
trend-adjusted basis — so this number is the one that actually matches
what was tested. The two are equal whenever no trend was detected.

**Which months, specifically — `bin_significant()`.** `result.p_value`
is a single number for the *whole* period. To ask which specific months
were driving it:

```python
from tide_lite import bin_significant

sig = bin_significant(result, alpha=0.05)      # bool array, one per bin
print(fit.bins[sig])                            # which bins were flagged
print(result.bin_p_values)                       # raw, uncorrected, per bin
```

This is deliberately **not** the same as comparing each bin's raw p-value
to 0.05 directly — with 13 bins, that would flag something roughly a
third of the time in a perfectly normal year, purely from running 13
tests instead of one (the same problem the whole-year statistic exists
to avoid, just recreated one level down). `bin_significant()` applies
Holm-Bonferroni across the bins within the one period, the same
correction Sequential uses across events, so the false-alarm rate stays
controlled. It's a genuinely different question from the global
p-value — a whole year can be significant overall without any single
bin surviving this stricter, per-bin correction, and vice versa; check
both rather than assuming they always agree. `plot_envelope` marks
these bins with a star automatically (`mark_significant_bins=True` by
default).

**The plot.** The shaded band is where historical years typically fall;
individual thin lines are the actual historical years, so you can see the
real spread, not just a summary. The bold treatment line sitting clearly
outside the band, for a sustained stretch, is what a genuinely significant
result should look like — if the p-value is small but the line barely
strays from the band, that's worth a second look. Black stars mark
individually significant bins (see above).

**Power curve** (if you ran `estimate_power`). "Given how much historical
data you actually have, how reliably could this test detect an effect of
a given size?" A flat, low power curve at effect sizes you'd care about
means a non-significant result might just reflect an underpowered test,
not genuine "no effect."

**Sensitivity grid** (if you ran `sensitivity_grid`). Whether your result
holds up under small, defensible changes to settings like bin size or
detrending method. A result that flips from significant to not under a
minor, reasonable change is a much weaker finding than one that's stable
across the grid.

## 8. Limitations — read this before trusting a result

These aren't small print. Each one changes how much weight a result can
reasonably carry.

**1. With fewer than ~15-20 historical years, one p-value is noisier than
it looks.** In simulation, the same test procedure's false-positive rate
varied from close to 0% up to 25-28% *depending on which specific 10
years happened to be available as history* — even though it's correctly
calibrated on average across many possible historical samples. Any real
analysis only ever has the one set of historical years you actually have.
**What to do:** always run `estimate_power` and `sensitivity_grid`
alongside a real result, especially if the p-value is anywhere near 0.05.
More historical years directly reduces this problem.

**2. "Unusual" is not "proven caused by the treatment."** Without a
comparison/control site, a significant result tells you the treatment
period was statistically unusual relative to history — it doesn't rule
out some other coincident cause (an unusually wet or dry period, a
different upstream event). **What to do:** check at least one obvious
alternative explanation (was rainfall/flow also anomalous that period?)
before attributing a result to the treatment specifically.

**3. The method assumes your historical years are genuinely comparable.**
If there's a real slow change across your historical period (climate
trend, land-use change) that isn't a simple, roughly-linear drift, the
built-in trend check may not fully capture it. **What to do:** look at
the historical years overlaid (the plot in Section 6) — they should read
as variations on one process, not a visibly drifting or splitting pattern.

**4. This is a new combination of established methods, not itself a
published, peer-reviewed method.** The individual pieces — block
bootstrap resampling, permutation testing, Holm-Bonferroni correction —
are well-established statistical tools with decades of use. This specific
assembly of them for this purpose hasn't been published or used
elsewhere. **What to do:** for internal research and monitoring, this is
a non-issue — cite the underlying methods directly. If a result
specifically needs to survive external legal or regulatory scrutiny,
treat that as a separate, additional requirement this tool doesn't yet
address (see Section 10 on what's deliberately not built).

**5. `detrend_mode="log_additive"` requires strictly positive data.** If
your variable can be zero or negative, either use `detrend_mode="additive"`
or offset your data first — the tool will raise a clear error rather than
silently mishandling this, but it's worth knowing upfront.

**6. `MonitoringSeries` needs more historical years than the other three
modes — this isn't optional guidance, it's a measured requirement.**
Splitting an already-small budget across many checks pushes each one
into the tail of the historical bootstrap distribution, which is
unreliable when built from limited data. In simulation, a nominally 5%
lifetime budget (checked over a 10-year horizon) had a *true*
false-alarm rate of roughly 17% with 10 historical years, 7% with 20,
and 3% with 30. **What to do:** use at least 20 historical years before
relying on `MonitoringSeries`; treat anything built on fewer as
indicative, not a real 5% guarantee. (A tool to check this against your
own specific data was attempted and deliberately left out — see the
module docstring in `tide_lite/monitoring.py` for why a first design
would have given false reassurance in exactly the cases that need a
warning most.)

## 9. Troubleshooting

**`Only N historical years survived the completeness filter (need >= 3)`**
— Fewer than 3 of your historical years had enough data (by default, at
least half their bins filled). Either provide more/more-complete
historical data, or loosen `max_missing_frac` in `TideConfig` if you
understand the tradeoff (more data, but noisier per-year estimates).

**`detrend_mode='log_additive' requires strictly positive values`** — Your
data has zero or negative values. Switch to `detrend_mode="additive"`, or
fix/offset the data if you specifically need log-space detrending.

**`mode='prediction' expects exactly one year of treatment data, but
treatment_df spans N years`** — Your treatment data accidentally covers
more than one year. Either restrict it to one year, or use
`mode="confidence"` if you actually intend to test the average across
multiple years.

**`every SequentialEvent needs an explicit year_index`** (or the
Cumulative/MonitoringSeries equivalent) — Your historical data has a
detected trend, and at least one event/window/check didn't specify
`year_index`. Set it explicitly — it's the number of years after your
first historical year that this specific period occurred.

**`Monitoring horizon of N checks is used up`** — You've called
`.check()` more times than `n_years_horizon` allows on a
`MonitoringSeries`. Call `.renew()` to start a fresh horizon with a new
budget — deliberately not automatic, since extending the guarantee is a
real decision.

**A result seems "too significant" given how the plot looks** — Check
`effect_size` vs. `effect_size_trend_adjusted` (Section 7). If your
historical data has a real trend, the p-value is computed on the
trend-adjusted basis, which can look different from what the raw plot
suggests at a glance.

## 10. Additional info

**Other documents in this project**, and when to read them:

- `METHODS_template.md` — fill this in before running a real analysis. It
  turns every setting above into a documented decision.
- `tide_rebuild_plan.md` — the phase-by-phase plan this project was built
  from, including why block bootstrap was chosen over alternatives, and
  what's deliberately *not* built (audit trails, legal-defensibility
  documentation, a benchmark suite against other methods — all gated
  behind an actual need, not built speculatively).
- `ADVERSARIAL_REVIEW.md` — a self-review that found and fixed six issues
  after the initial build, including two genuine statistical bugs. Worth
  reading if you want to know exactly what's been checked.

**The statistics, briefly, for the technically curious:**

- **Block bootstrap** (Künsch 1989; Politis & Romano 1992/1994): instead
  of assuming a distribution, this resamples chunks of your own
  historical data to build up a picture of plausible "normal" variation.
  Chunks are matched to the calendar position they came from (so a
  volatile month's pattern never gets tested against a stable month's
  typical spread) and pulled from a small neighborhood of nearby
  positions for enough variety to resolve a precise p-value.
- **Median and MAD** (median absolute deviation) instead of mean and
  standard deviation: robust to the skew that's common in concentration
  data, without assuming a bell curve.
- **Mann-Kendall test and Sen's slope** (Mann 1945; Kendall 1975; Sen
  1968): detect and estimate a long-term trend in the historical data
  before testing, so a real ongoing drift doesn't get mistaken for a
  treatment effect (or vice versa).
- **Holm-Bonferroni correction** (Holm 1979): used in Sequential mode to
  control the chance of a false alarm across several tests at once,
  without being as overly conservative as the simpler Bonferroni method.

Full implementation detail, with the reasoning for every design choice
written as comments, is in `tide_lite/engine.py`, `controllers.py`, and
`power.py` — the code is meant to be read, not just run.

## 11. Glossary

- **p-value** — how surprising your result would be if there were truly
  no effect. Small = surprising = evidence of a real difference.
- **Effect size** — how big the difference actually was, in your data's
  own units (separate from how *sure* you are that it's real).
- **Envelope** — the band of historical "normal" variation a treatment
  period is compared against.
- **Seasonality** — the regular within-year pattern (e.g. warmer in
  summer) that this tool automatically accounts for.
- **Autocorrelation** — the tendency for a variable's value on one day to
  be related to its value on nearby days; ignoring this makes statistical
  tests overconfident.
- **Bootstrap / resampling** — building up a picture of plausible
  variation by repeatedly drawing from your own historical data, instead
  of assuming a textbook distribution.
- **Detrending** — removing a genuine long-term drift from the historical
  data before testing, so it isn't mistaken for a treatment effect.
- **Family-wise error rate (FWER)** — the chance of at least one false
  alarm across a whole set of tests, as opposed to just one.
- **Exchangeability** — the assumption that your historical years are all
  reasonably comparable draws from "the same normal," rather than
  fundamentally different from each other.
- **Pointwise vs. corrected (bin-level)** — pointwise means comparing
  each month's raw p-value to 0.05 on its own, which racks up false
  alarms just from testing many months; corrected (`bin_significant()`)
  controls for that, the same way Sequential does across events.
