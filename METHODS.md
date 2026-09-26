# tide_lite: method statement

This document describes everything tide_lite does, what it assumes, and
how well it performs. It is written for reviewers, for technical
appendices, and for anyone who has to explain a result to a regulator or a
court. For how to *use* the tool, see `USER_GUIDE.md`.

Version: tide_lite 0.4.0.

---

## 1. The question

> Was this period's data unusual compared with the same site's own
> normal years, once the seasonal cycle, any steady long-term trend, and
> ordinary year-to-year variation are allowed for?

**Null hypothesis (H0).** The tested year behaves like a randomly chosen
historical year: same seasonal pattern, same spread, same kind of
year-to-year variation, and (if a trend was detected) the same trend
continued to its date.

Note what the last part means: when a trend is removed, a year that simply
continues the historical trend counts as normal. If the question is
"has the site changed from its historical level?", including any steady
drift, set `detrend_mode="none"` before testing.

**Alternative.** Chosen before looking at the data. `"two-sided"`
(default): the year is unusually high *or* low in at least one bin.
`"greater"` / `"less"`: unusually high only / low only.

**Output.** A p-value for the whole period, the bins (months by default)
responsible, the typical size of the difference in the data's own units
(the median over tested bins of tested value minus expected value), and a
plot. The p-value is the chance that a normal year would look at least
this unusual.

What the test does **not** establish: the *cause* of a difference, or
whether values are acceptable against an external limit or guideline.

## 2. The procedure, step by step

Every step below is in `tide_lite/engine.py` (`fit_historical`,
`test_treatment`). Settings are in `TideConfig`; defaults in brackets.

| # | Step | What happens | Why |
|---|---|---|---|
| 1 | **Bin** | Each calendar year is cut into bins (calendar months, or fixed blocks of `bin_days` days; with fixed blocks, 29 February is dropped so every year has 365 days, and a short leftover at year end joins the last bin). Each bin is summarised by the median (`agg`) of its measurements. | Compare like with like: January with January. Medians resist odd single readings. |
| 2 | **Screen** | A bin in a given year counts as missing if it has fewer than 75% (`min_bin_coverage`) of that bin's usual number of measurements; this applies to the tested period too. Then, in order: bins with data in fewer than half the years are dropped, then years missing more than half their bins (`max_missing_frac`), then bins left with under 3 years. Every drop is reported. | A 3-day "month" is far noisier than a 30-day one and must not be compared with it. Too little data cannot define "normal". |
| 3 | **Transform** (optional) | `detrend_mode="log_additive"` works on natural logarithms. | Makes changes proportional (e.g. +20%) for skewed, positive data. |
| 4 | **Trend** | Each year's typical anomaly (median over bins of its departure from that bin's historical median) is tested with Mann-Kendall. If p < `mk_alpha` (0.10), the Sen slope per calendar year is removed from every year. Never done with `detrend_mode="none"`. | A steady long-term drift is not treated as an effect. Anomalies, not raw values, so a year missing its summer does not look like a trend. |
| 5 | **Baseline** | Per bin: the median across years (typical level) and 1.4826 × the median absolute deviation, MAD (typical spread). If a bin's spread is zero (half or more of its values identical, e.g. non-detects), the median spread of the other bins is used and a warning is given. | Robust to skew and outliers; the 1.4826 factor puts MAD on the scale of a standard deviation. |
| 6 | **Leave-one-year-out** | Each historical year is compared with the median and MAD of the *other* years. | The tested year is compared with a median and spread it did not help estimate; historical years must be too, or they look more ordinary than a new year. |
| 7 | **Split** | Each historical year's departures are split into an **annual level** (mean departure across bins, data units: "was it a high year overall?") and a **within-year pattern** (the rest, divided by the spread). | Whole-year shifts (wet/dry years) and month-to-month wobble behave differently and are resampled differently. |
| 8 | **Null distribution** | `n_bootstrap` (2000) synthetic normal years. Each = a within-year pattern stitched from blocks of consecutive bins, each block taken from a random historical year at exactly the bins it fills (`block_length_bins`, default a quarter of a year; seams at a random position) + one annual level drawn from a Student-t prediction distribution fitted to the historical annual levels, divided by each bin's spread. | Blocks keep the persistence between neighbouring months; drawing from many years gives far more distinct synthetic years than the N real ones. The Student-t step lets a new year's level be more extreme than any seen, as real new years sometimes are, widening for short records and for projecting a trend forward. |
| 9 | **Score** | Each bin of the tested year: score = (value − typical level at its date) / spread, the trend being projected to the tested year's own calendar year. | Same yardstick as the null. |
| 10 | **p-value** | Statistic = largest bin score (largest absolute score for two-sided), over the bins the tested period has data for; the synthetic years are scored over the same bins. p = (number of synthetic years with a statistic at least as large + 1) / (`n_bootstrap` + 1). | One number for the whole period with no inflation from testing many bins. The +1 terms mean the p-value is never 0 and random draws cannot make it too small (Phipson & Smyth 2010). |
| 11 | **Which bins** | A bin is flagged when its own score exceeds the same critical value (single-step max-T). | The period is significant **if and only if** at least one bin is flagged, and the plotted band is exactly that critical value: the plot, the flags and the p-value cannot disagree. |

The tested period is refused if it spans more than one calendar year
(except in `mode="confidence"`), overlaps the historical years, or is
missing more than `max_missing_frac` of the bins, the same completeness
bar historical years must clear. Missing bins are named in a warning and
never count as evidence either way. A warning is also given if the tested
period is sampled more than twice as often as the history.

**Optional extras on the same null draws:**
- *Sustained departure* (`run_lengths=[k]`): the largest average score over
  k consecutive bins, which must lean one way. It is more sensitive to a
  moderate multi-month shift, and less sensitive to a one-month spike.
  Windows never span a missing or dropped bin.
- *Several years at once* (`mode="confidence"`): years averaged bin by bin
  and compared with averages of the same number of synthetic years. The
  error in the historical medians is shared by all the averaged years, so
  it is not allowed to shrink with averaging.
- *Several periods*: Holm (1979) correction across events (`sequential_test`),
  or a fixed per-check threshold α/horizon (`MonitoringSeries`).

## 3. Assumptions

The error rates in section 4 hold when:

1. **The historical years are genuinely "normal"** for the question being
   asked: no treatment, no disturbance, and comparable with each other
   (same site, sampling method, laboratory, units; no step changes).
2. **Years are roughly independent** apart from any straight-line trend.
   Multi-year runs (a several-year drought) reduce the effective number of
   years and make the test flag too often.
3. **Any long-term change is roughly a straight line** over the record
   (in logs, if `log_additive`). Curved or stepped changes are not removed.
4. **Dependence within a year is short-range**, spanning a few bins at
   most, so it is carried by the resampled blocks.
5. **Whole-year shifts are not strongly lopsided or heavy-tailed.** This is
   the method's one distributional assumption, and it applies to the
   single annual level only. Section 4 shows the measured cost when it is
   violated.
6. **Measurements are consistent over time.** Non-detects are converted to
   numbers the same way every year, and the sampling frequency of the
   tested period is similar to the history's (enforced below 75%, warned
   above 200%).
7. **The analysis settings were fixed before looking at the tested
   period.** Trying several settings or directions and reporting the best
   invalidates the p-value; use `sensitivity_grid` to show stability.

## 4. Validation: measured error rates

How often does the test flag a year that is genuinely normal? Measured by
simulating a site's history from a known process, fitting it, and testing
brand-new years from the same process (so every flag is a false alarm);
3,000 tests per row, spread over 300 simulated histories. Reproduce with
`python tests/calibration_study.py` (about 10 minutes), which also prints
each figure's standard error.

The simulated site has a seasonal cycle, day-to-day noise that persists
for days, and a random whole-year shift. As a yardstick, call one **unit**
the typical year-to-year wobble of a single month's value caused by the
day-to-day noise (its standard deviation). "Moderate" whole-year shifts
have a standard deviation of about 1.35 units; "strong" ones about 3.2
units.

**Precision of these figures.** Each rate has a standard error of about
0.5 percentage points (0.8 for the three-year average), so rates within
about 1 point of 5% are consistent with exactly 5%. The figures are
averages over many histories. The rate for any *one* history varies more,
as for any method that estimates "normal" from a sample: for 9 in 10
simulated histories it was between 0% and 16% at a nominal 5% with 10
years of history, and between 0.5% and 12% with 20.

**False-alarm rate at a nominal 5% (α = 0.05):**

| Historical years | No whole-year shifts | Moderate | Strong |
|---|---|---|---|
| 5 | 5.2% | 3.1% | 2.6% |
| 10 | 4.0% | 3.0% | 4.3% |
| 20 | 5.1% | 5.0% | 5.1% |
| 40 | 5.0% | 5.4% | 5.0% |

| Other conditions (10 and 20 years, moderate shifts unless stated) | At α = 0.05 |
|---|---|
| Heavy-tailed whole-year shifts (strong) | 5.0%, 6.3% |
| Strongly skewed whole-year shifts (strong) | 6.6%, 6.2% |
| Real trend present, next year tested | 3.8%, 5.8% |
| Log-normal data, `log_additive` | 3.4%, 5.2% |
| One grab sample per month | 3.1%, 4.9% |
| Part year tested (Jan–Jun only) | 4.2%, 4.4% |
| Weekly bins (`bin_days=7`) | 3.9%, 4.2% |
| One-sided (`"greater"`) | 5.0%, 4.7% |
| Sustained test, 3 months | 4.1%, 4.8% |
| Three years averaged (`mode="confidence"`) | 4.2%, 5.8% |

**Summary.** With normally distributed whole-year shifts, the test flagged
**2.6%–5.8% of normal years at a nominal 5%**, **6.1%–11.0% at 10%** and
**0.2%–2.2% at 1%**. With strongly skewed or heavy-tailed whole-year
shifts (assumption 5), it flagged up to **6.6% at 5%** and **3.3% at 1%**.
The test is close to its stated rate at 5% and 10%. At 1% it can flag up
to about 3 times more often than stated, because the extreme tail of any
resampling method is limited by the number of historical years. So a
p-value at or below 0.01 is strong evidence, but should not be quoted as a
1% error rate; use α = 0.05 or 0.10 for formal decisions.

**MonitoringSeries** (chance of *ever* flagging a normal year over a
10-year horizon, 5% budget, so α = 0.005 per check; standard error about
1–1.5 points): 1.8%–12.3% with 10 historical years, 5.0%–9.8% with 20,
5.3%–8.8% with 40. The per-check threshold is deep in the tail, so expect
up to about twice the stated budget even with a long record.

**Earlier versions.** Version 0.3.0 flagged 3.4%–7.1% at a nominal 5% and
up to 4.2% at a nominal 1%; version 0.2.0 flagged 2.7%–16.1% and up to
12.6%. Results from earlier versions should be re-run (`CHANGELOG.md`).

**Power** (chance of flagging a year that really changed, α = 0.05,
units as above; 1,000 tests per cell):

| Historical years, whole-year shifts | Whole year +2 units | Whole year +4 units | May–Jun +4 units | May–Jun +8 units |
|---|---|---|---|---|
| 10, none | 41% | 89% | 46% | 91% |
| 10, moderate | 13% | 42% | 16% | 57% |
| 20, none | 64% | 99% | 73% | 99% |
| 20, moderate | 20% | 60% | 30% | 89% |

A change that lifts the whole year is hard to tell apart from a naturally
wet or dry year, especially with a short record. Use `estimate_power` on
your own fit for your own numbers.

## 5. Limitations

- **Unusual is not caused.** Without a control site the test cannot
  separate an intervention from a coincident cause (weather, flow, an
  upstream event). Check the obvious alternatives before attributing a
  result.
- **Not significant is not "no change".** It means no change was
  detected. Report power (`estimate_power`) alongside a non-significant
  result; with a short record, only large changes are detectable.
- **Short records.** Below about 10 years only large changes are
  detectable. The false-alarm rate stays at or below the stated rate
  because the Student-t step widens the null, but the method cannot create
  information the record does not hold.
- **Projecting a trend.** A removed trend is projected to the tested
  year. The further past the record, the more the answer depends on that
  projection. When the removed trend was not real, the false-alarm rate at
  a nominal 5% was about 8% one year past the record and 16–29% ten years
  past (warned beyond 5 years). A real trend that was too weak to detect
  and so not removed makes distant years look unusual, which is correct if
  the question is change from the historical level. Test years close to
  the record where possible, and show the result with `detrend_mode="none"`
  too.
- **Calendar years.** Periods are calendar years (January–December).
  A part year is fine; a water year spanning two calendar years must be
  tested as two parts or as separate events.
- **One site, one variable.** Several variables or sites are separate
  tests; correct across them (e.g. `holm_bonferroni`) rather than
  reporting the smallest p-value.
- **Monte Carlo p-values** vary slightly with the random seed; the
  smallest attainable value is 1/(`n_bootstrap` + 1). Always report the
  seed.

## 6. What is established and what is new

Every component is standard, published statistics:

| Component | Reference |
|---|---|
| Mann-Kendall trend test; Sen's slope | Mann (1945); Kendall (1975); Sen (1968); Hirsch, Slack & Smith (1982); Helsel et al. (2020) |
| Median and MAD | Rousseeuw & Croux (1993) |
| Block bootstrap for dependent data | Künsch (1989); Politis & Romano (1992) |
| Blocks matched to the same season | Politis (2001); Dudek, Leśkow, Paparoditis & Politis (2014) |
| Monte Carlo p-value (b+1)/(B+1) | Davison & Hinkley (1997); Phipson & Smyth (2010) |
| Max-T simultaneous test and band | Westfall & Young (1993) |
| Prediction distribution for a new observation | Hahn & Meeker (1991) |
| Holm step-down correction | Holm (1979) |

**Closest established practice.** Comparing a site's new data with its own
history is the idea behind the *intrawell* prediction limits and control
charts of the US EPA Unified Guidance for groundwater monitoring (USEPA
2009). tide_lite asks the same kind of question, but allows for a seasonal
cycle, a trend, whole-year shifts and dependence between months without
assuming the data follow a bell curve, and reports which months differed.

**What is new** is the combination, tailored to "was this year unusual for
this site": the leave-one-year-out scoring of historical years, the
separate treatment of the annual level with a Student-t prediction step,
and the single max-T decision that ties the p-value, the flagged months
and the plotted band together. To our knowledge this combination has not
been published as a single procedure, so it has not been peer-reviewed as
a whole. Instead, this repository provides a reproducible simulation
study with measured error rates under a range of conditions (section 4),
and a test suite (`tests/`) that checks each step.

## 7. Using a result as evidence

A reviewer or opposing expert can reasonably ask the questions below.
tide_lite is designed so that each has a documented answer.

| Question | Answer |
|---|---|
| Is the method testable, and what is its error rate? | Yes: section 4, reproducible with `tests/calibration_study.py`. |
| Are its parts published and peer-reviewed? | Yes, each component (section 6). The combination is not (section 6). |
| Is it standard practice? | The approach follows established intrawell comparison (section 6); the implementation is specific to this package. |
| Were the settings chosen after seeing the answer? | Record settings, alternative and alpha before the analysis (assumption 7). `sensitivity_grid` shows whether other reasonable settings change the answer. |
| Can the number be reproduced? | Yes, exactly: section 8. |
| Were the historical years really normal? | A matter of site knowledge, not statistics (assumption 1). Document why each year was included. |
| Does "significant" prove the activity caused it? | No (section 5). |

Admissibility rules differ between jurisdictions. Have the analysis
reviewed by a qualified statistician when it matters.

## 8. Reproducing a result

A result is fully reproducible from: the input data, the tide_lite
version, every `TideConfig` setting, the mode and alternative, and the
random seed (`rng=`). `summarize(fit, result)` prints all of these plus a
fingerprint of the fitted baseline, which changes if the data or any
setting changes. Two runs with the same inputs and seed give identical
numbers.

## References

- Davison, A.C. & Hinkley, D.V. (1997). *Bootstrap Methods and their Application*. Cambridge University Press.
- Dudek, A.E., Leśkow, J., Paparoditis, E. & Politis, D.N. (2014). A generalized block bootstrap for seasonal time series. *Journal of Time Series Analysis* 35(2), 89–114.
- Hahn, G.J. & Meeker, W.Q. (1991). *Statistical Intervals: A Guide for Practitioners*. Wiley.
- Helsel, D.R., Hirsch, R.M., Ryberg, K.R., Archfield, S.A. & Gilroy, E.J. (2020). *Statistical Methods in Water Resources*. USGS Techniques and Methods 4-A3.
- Hirsch, R.M., Slack, J.R. & Smith, R.A. (1982). Techniques of trend analysis for monthly water quality data. *Water Resources Research* 18(1), 107–121.
- Holm, S. (1979). A simple sequentially rejective multiple test procedure. *Scandinavian Journal of Statistics* 6(2), 65–70.
- Kendall, M.G. (1975). *Rank Correlation Methods*, 4th ed. Griffin.
- Künsch, H.R. (1989). The jackknife and the bootstrap for general stationary observations. *Annals of Statistics* 17(3), 1217–1241.
- Mann, H.B. (1945). Nonparametric tests against trend. *Econometrica* 13(3), 245–259.
- Phipson, B. & Smyth, G.K. (2010). Permutation p-values should never be zero. *Statistical Applications in Genetics and Molecular Biology* 9(1), Article 39.
- Politis, D.N. (2001). Resampling time series with seasonal components. *Proceedings of the 33rd Symposium on the Interface of Computing Science and Statistics*.
- Politis, D.N. & Romano, J.P. (1992). A circular block-resampling procedure for stationary data. In *Exploring the Limits of Bootstrap*, Wiley, 263–270.
- Rousseeuw, P.J. & Croux, C. (1993). Alternatives to the median absolute deviation. *Journal of the American Statistical Association* 88(424), 1273–1283.
- Sen, P.K. (1968). Estimates of the regression coefficient based on Kendall's tau. *Journal of the American Statistical Association* 63(324), 1379–1389.
- USEPA (2009). *Statistical Analysis of Groundwater Monitoring Data at RCRA Facilities: Unified Guidance*. EPA 530/R-09-007. US Environmental Protection Agency, Office of Resource Conservation and Recovery.
- Westfall, P.H. & Young, S.S. (1993). *Resampling-Based Multiple Testing*. Wiley.
