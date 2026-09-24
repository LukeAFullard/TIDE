# tide_lite

**Was this period unusual compared with the site's own normal years?**

You have a record of measurements from one site (water quality, air
quality, flows, anything measured through the year) and a period you want
to ask about: the year after a discharge started, a restoration, a
pollution event, or simply this year. tide_lite compares that period with
the site's own history, allowing for:

- **the seasons**: January is only ever compared with January;
- **a steady long-term trend**, if the record has one;
- **ordinary year-to-year variation**, such as wet and dry years;

and it does not assume the data follow a bell curve.

```python
from tide_lite import fit_historical, test_treatment, summarize, plot_envelope

fit    = fit_historical(history_df, date_col="date", value_col="value")
result = test_treatment(fit, this_year_df, rng=42)      # rng: seed, for reproducibility

print(summarize(fit, result))                           # plain-language answer
plot_envelope(fit, result).figure.savefig("result.png")
```

`summarize` answers in words: whether the period was outside the normal
range, which months, by how much, what caveats apply to your data, and the
exact settings needed to reproduce the result.

## Is this a compliance test?

It answers **"did this site change from its own normal?"** That is the
question behind conditions such as "no significant change from baseline"
or "no more than minor effects". It does **not** check values against a
fixed limit or guideline; compare with the limit directly for that. Two
cautions apply to every result:

- **Unusual is not the same as caused by an activity.** Check other
  explanations (weather, flow, upstream events, changes of method).
- **"Not significant" means "no change detected", not "no change".** A
  short or noisy record can only detect large changes; `estimate_power`
  tells you how large.

## Four ways to use it

| You want to know | Use |
|---|---|
| Was this one period unusual? | `test_treatment` |
| Which of several separate years were unusual? (corrected for testing several) | `sequential_test` |
| After a known event, is its effect still visible, year by year? | `cumulative_test` + `first_recovery` |
| Each new year as it arrives, with a fixed false-alarm budget over many years | `MonitoringSeries` (needs 20+ years of history) |

Also: which months (`bin_significant`), a sustained multi-month departure
(`run_lengths=[3]`), one-sided questions such as "higher than normal?"
(`alternative="greater"`), detectable effect size (`estimate_power`), and
stability under other reasonable settings (`sensitivity_grid`).

## How much can you trust it?

Measured by simulation (`tests/calibration_study.py`): with 5–40 years of
history, the test flagged **3.4%–7.1% of genuinely normal years at a
stated 5%**. That covers whole-year shifts, trends, skewed data, monthly
grab samples and part years. Details, assumptions and limitations are in
`METHODS.md`.

## Install and check

```
pip install -e .                  # numpy, pandas, scipy, matplotlib
python examples/run_example.py    # tutorial on synthetic data; writes example plots
python -m pytest tests/ -q        # 57 checks, about 30 seconds
```

## Documents

| File | For |
|---|---|
| `QUICKSTART.md` | A first result in 5 minutes |
| `USER_GUIDE.md` | Preparing data, choosing settings and modes, reading results, troubleshooting |
| `METHODS.md` | Every step of the method, its assumptions, measured error rates, references |
| `CHANGELOG.md` | What changed between versions (re-run results from 0.2.0) |
| `examples/run_example.py` | Every feature, end to end, on synthetic data |

## Layout

```
tide_lite/
  engine.py       fit_historical, test_treatment, bins, band, trend tests
  controllers.py  sequential_test, cumulative_test, first_recovery
  monitoring.py   MonitoringSeries
  power.py        estimate_power, sensitivity_grid
  plotting.py     one plot per mode
  report.py       summarize
tests/
  test_engine.py         statistical behaviour of the method
  test_regressions.py    one test per implementation defect found in audits
  test_monitoring.py     MonitoringSeries
  calibration_study.py   the simulation study behind METHODS.md section 4
```
