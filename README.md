# tide_lite

A from-scratch, lean rebuild: is a treatment period's water quality data
unusual relative to historical "normal" years, accounting for seasonality
and autocorrelation, without assuming a distribution?

Four modes, one question each:

- **Standard** (`test_treatment`) -- did this one period behave abnormally?
  Also answers "which months specifically" (`bin_significant`), "was it a
  sustained multi-month shift, not just one noisy month"
  (`run_lengths=` / `result.sustained`), and "when did it start and did
  it recover" (`estimate_departure_recovery`).
- **Sequential** (`sequential_test`) -- across several independent events,
  which were real? (Holm-Bonferroni corrected.)
- **Cumulative** (`cumulative_test`) -- how is this one known event
  progressing? (Same fixed baseline throughout, uncorrected.) A formal
  rule for declaring "recovered" -- requiring several consecutive
  non-significant windows, not just one -- is `first_recovery`.
- **MonitoringSeries** (`MonitoringSeries`) -- ongoing, open-ended yearly
  monitoring with no specific "treatment" -- check each new year as it
  arrives, indefinitely, without the false-alarm risk piling up the
  longer you keep watching. Stateful (save/load between runs), unlike
  the other three. **Needs more historical years than the other modes to
  be reliable** -- see `USER_GUIDE.md`.

Plus `estimate_power` / `estimate_power_sequential` (can this test even
detect an effect size you'd care about, given your real historical sample)
and `sensitivity_grid` (is a real result stable to defensible alternative
settings).

## Setup

```
pip install numpy pandas scipy matplotlib
```

No other dependencies -- Mann-Kendall and Sen's slope are implemented
directly rather than pulling in an extra package.

## Start here

```
python examples/run_example.py   # full walkthrough on synthetic data,
                                  # produces examples/envelope_plot.png
python tests/test_engine.py      # Phase 3/4 validation checks
```

**New to this project? Read `QUICKSTART.md` first (5 minutes), then
`USER_GUIDE.md`** for the complete picture — what each mode is for, how
to interpret a result, limitations worth knowing before you trust one,
and a troubleshooting section for every error message the code raises.

`examples/run_example.py` is written as a tutorial too — read it top to
bottom. Swap the synthetic data generator at the top for a real
`pd.read_csv(...)` and everything downstream is unchanged.


## A finding worth knowing before you interpret a real result

With roughly 10 historical years, the false-positive rate of a single
fitted model can vary a fair amount fit to fit (0-25% at nominal alpha=0.05
in simulation) even though it averages out correctly across many different
possible historical samples. Your real analysis only ever has one
historical fit. This is exactly why Phase 5 (power and sensitivity
analysis) exists -- run them on your real result before leaning on a
single p-value. See the docstring at the top of `tests/test_engine.py` for
the simulation this is based on.

## Layout

```
tide_lite/
  engine.py       -- Phase 2: fit_historical, test_treatment (core)
  controllers.py  -- Phase 4: Sequential, Cumulative, Holm-Bonferroni
  power.py        -- Phase 5: power analysis, sensitivity harness
  monitoring.py   -- MonitoringSeries: ongoing yearly monitoring (added after Phase 6)
  plotting.py     -- a dedicated plot function per mode (added on request)
tests/
  test_engine.py     -- Phase 3-4 validation checks
  test_monitoring.py -- MonitoringSeries validation checks
examples/
  run_example.py  -- Phase 6: full walkthrough, all four modes plotted
QUICKSTART.md         -- 5-minute start
USER_GUIDE.md          -- comprehensive guide: modes, interpretation, limitations, troubleshooting
```
