# tide_lite: method statement

This document is the complete description of what tide_lite does, what it
assumes, and how well it performs. It is written for reviewers, for
technical appendices, and for anyone who has to explain a result to a
regulator or a court. For how to *use* the tool, see `USER_GUIDE.md`.

Version: tide_lite 0.3.0.

---

## 1. The question

> Was this period's data unusual compared with the same site's own
> normal years, once the seasonal cycle, any steady long-term trend, and
> ordinary year-to-year variation are allowed for?

**Null hypothesis (H0).** The tested year behaves like a randomly chosen
historical year: same seasonal pattern, same spread, same kind of
year-to-year variation, and (if a trend was detected) the same trend
continued to its date.

**Alternative.** Chosen before looking at the data. `"two-sided"`
(default): the year is unusually high *or* low in at least one bin.
`"greater"` / `"less"`: unusually high only / low only.

**Output.** A p-value for the whole period, the bins (months by default)
responsible, the typical size of the difference in the data's own units,
and a plot. The p-value is the chance that a normal year would look at
least this unusual.

What the test does **not** establish: the *cause* of a difference, or
whether values are acceptable against an external limit or guideline.

## 2. The procedure, step by step

Every step below is in `tide_lite/engine.py` (`fit_historical`,
`test_treatment`). Settings are in `TideConfig`; defaults in brackets.

| # | Step | What happens | Why |
|---|---|---|---|
| 1 | **Bin** | Each calendar year is cut into bins (calendar months, or fixed blocks of `bin_days` days). Each bin is summarised by the median (`agg`) of its measurements. | Compare like with like: January with January. Medians resist odd single readings. |
| 2 | **Screen** | A bin in a given year counts as missing if it has fewer than 75% (`min_bin_coverage`) of that bin's usual number of measurements. Bins covered in fewer than half the years, then years missing more than half their bins (`max_missing_frac`), then bins with under 3 years of data, are dropped. All drops are reported. | A 3-day "month" is far noisier than a 30-day one and must not be compared with it. Too little data cannot define "normal". |
| 3 | **Transform** (optional) | `detrend_mode="log_additive"` works on natural logarithms. | Makes changes proportional (e.g. +20%) for skewed, positive data. |
| 4 | **Trend** | Each year's typical anomaly (median over bins of its departure from that bin's historical median) is tested with Mann-Kendall. If p < `mk_alpha` (0.10), the Sen slope per calendar year is removed from every year. | A steady long-term drift is not treated as an effect. Anomalies, not raw values, so a year missing its summer does not look like a trend. |
| 5 | **Baseline** | Per bin: the median across years (typical level) and 1.4826 × the median absolute deviation, MAD (typical spread). | Robust to skew and outliers; the 1.4826 factor puts MAD on the scale of a standard deviation. |
| 6 | **Leave-one-year-out** | Each historical year is compared with the median and MAD of the *other* years. | This is exactly how the tested year is compared with all of them, so historical years carry the same estimation error a genuinely new year does. |
| 7 | **Split** | Each historical year's departures are split into an **annual level** (mean departure across bins, data units: "was it a high year overall?") and a **within-year pattern** (the rest, divided by the spread). | Whole-year shifts (wet/dry years) and month-to-month wobble behave differently and are resampled differently. |
| 8 | **Null distribution** | `n_bootstrap` (2000) synthetic normal years. Each = a within-year pattern stitched from blocks of consecutive bins taken from historical years (`block_length_bins`, default a quarter of a year; each block from the same time of year ± `pool_window_radius` = 1 bin; seams at a random position) + one annual level drawn from a Student-t prediction distribution fitted to the historical annual levels. | Blocks keep the persistence between neighbouring months; drawing from many years gives far more distinct synthetic years than the N real ones. The Student-t step lets a new year's level be more extreme than any seen, as real new years sometimes are, widening appropriately for short records and for projecting a trend forward. |
| 9 | **Score** | Each bin of the tested year: score = (value − typical level at its date) / spread, the trend being projected to the tested year's own calendar year. | Same yardstick as the null. |
| 10 | **p-value** | Statistic = largest bin score (largest absolute score for two-sided). p = (number of synthetic years with a statistic at least as large + 1) / (`n_bootstrap` + 1). | One number for the whole period with no inflation from testing many bins. The +1 terms make the Monte Carlo test exact rather than optimistic. |
| 11 | **Which bins** | A bin is flagged when its own score exceeds the same critical value (single-step max-T). | The period is significant **if and only if** at least one bin is flagged, and the plotted band is exactly that critical value: the plot, the flags and the p-value cannot disagree. |

**Optional extras on the same null draws:**
- *Sustained departure* (`run_lengths=[k]`): the largest average score over
  k consecutive bins, which must lean one way. It is more sensitive to a
  moderate multi-month shift, and less sensitive to a one-month spike.
- *Several years at once* (`mode="confidence"`): years averaged bin by bin
  and compared with averages of the same number of synthetic years.
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
5. **Whole-year shifts are not extremely lopsided.** This is the method's
   one distributional assumption, and it applies to the single annual
   level only. Section 4 shows the measured cost when it is violated.
6. **Measurements are consistent over time.** Non-detects are converted to
   numbers the same way every year; sampling frequency is similar (the
   75% coverage rule enforces the worst mismatches).
7. **The analysis settings were fixed before looking at the tested
   period.** Trying several settings or directions and reporting the best
   invalidates the p-value; use `sensitivity_grid` to show stability.

## 4. Validation: measured error rates

How often does the test flag a year that is genuinely normal? Measured by
simulating a site's history from a known process, fitting it, and testing
brand-new years from the same process (so every flag is a false alarm);
1,500 tests per row, spread over 150 simulated histories. Reproduce with
`python tests/calibration_study.py` (about 5 minutes).

The simulated site has a seasonal cycle, day-to-day noise that persists
for days, and a random whole-year shift. As a yardstick, call one **unit**
the typical year-to-year wobble of a single month's value caused by the
day-to-day noise (its standard deviation). "Moderate" whole-year shifts
have a standard deviation of about 1.35 units; "strong" ones about 3.2
units.

**False-alarm rate at a nominal 5% (α = 0.05):**

| Historical years | No whole-year shifts | Moderate | Strong |
|---|---|---|---|
| 5 | 6.5% | 3.4% | 3.4% |
| 10 | 3.9% | 3.9% | 5.1% |
| 20 | 6.5% | 6.1% | 5.3% |
| 40 | 4.9% | 6.9% | 4.5% |

| Other conditions (10 and 20 years, moderate shifts) | At α = 0.05 |
|---|---|
| Heavy-tailed whole-year shifts (strong) | 5.2%, 6.6% |
| Strongly skewed whole-year shifts (strong) | 7.1%, 6.5% |
| Real trend present, next year tested | 4.9%, 7.0% |
| Log-normal data, `log_additive` | 3.8%, 6.0% |
| One grab sample per month | 4.0%, 5.3% |
| Part year tested (Jan–Jun only) | 5.3%, 4.9% |
| Weekly bins (`bin_days=7`) | 3.7%, 6.1% |
| One-sided (`"greater"`) | 6.0%, 5.1% |
| Sustained test, 3 months | 3.4%, 4.0% |
| Three years averaged (`mode="confidence"`) | 7.5%, 8.3% |

Across all standard single-year conditions: **3.4%–7.1% at a nominal 5%**,
**6.8%–12.9% at a nominal 10%**, and **0.3%–4.2% at a nominal 1%**. The
test is close to its stated rate at 5% and 10%. At 1% it can flag up to
about 4 times more often than stated, because the extreme tail of any
resampling method is limited by the number of historical years. So a
p-value at or below 0.01 is at least as strong as a well-calibrated 5%
result, but should not be quoted as a 1% result; do not use α below 0.05
for a formal decision unless you have 20+ years and accept that caveat.

**MonitoringSeries** (chance of *ever* flagging a normal year over a
10-year horizon, 5% budget, so α = 0.005 per check): 2.7%–15% with 10
historical years, 7.0%–8.7% with 20, 5.3%–9.7% with 40. The per-check
threshold is deep in the tail, so expect up to about twice the stated
budget even with a long record.

**For comparison, version 0.2.0** flagged 2.7%–16.1% of normal years at a
nominal 5% and up to 12.6% at a nominal 1%; strong whole-year shifts with a
10-year record were its worst case. Results from 0.2.0 should be re-run.

**Power** (chance of flagging a year that really changed, α = 0.05,
units as above):

| Historical years, whole-year shifts | Whole year +2 units | Whole year +4 units | Two months +4 units | Two months +8 units |
|---|---|---|---|---|
| 10, none | 43% | 89% | 51% | 93% |
| 10, moderate | 18% | 50% | 18% | 64% |
| 20, none | 64% | 99% | 75% | 100% |
| 20, moderate | 23% | 64% | 36% | 93% |

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
  detectable; the method stays valid because the Student-t step widens
  the null, but it cannot create information the record does not hold.
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

What is new is the combination, tailored to "was this year unusual for
this site": the leave-one-year-out scoring of historical years, the
separate treatment of the annual level with a Student-t prediction step,
and the single max-T decision that ties the p-value, the flagged months
and the plotted band together. To our knowledge this combination has not
been published as a single procedure. It has therefore not been
peer-reviewed as a whole. In its place this repository provides a
reproducible simulation study with measured error rates under a range of
conditions (section 4), and a test suite (`tests/`) that checks each step.

## 7. Reproducing a result

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
- Westfall, P.H. & Young, S.S. (1993). *Resampling-Based Multiple Testing*. Wiley.
