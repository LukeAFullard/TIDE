# Changelog

## 0.3.0

An independent audit re-measured the method's error rates instead of
relying on earlier documentation. It found that 0.2.0 flagged normal years
far more often than stated in common conditions. **Re-run any result
produced with 0.2.0.** Measured false-alarm rates at a stated 5% went from
2.7%–16.1% (0.2.0) to 3.4%–7.1% (0.3.0); at a stated 1%, from up to 12.6%
to up to 4.2%. Details are in `METHODS.md` section 4, and the study that
measures them is `tests/calibration_study.py`.

### Method changes (why the numbers changed)

- **Leave-one-year-out scoring.** Historical years were scored against a
  median and spread that they themselves helped estimate, which made them
  look more ordinary than a genuinely new year. Each is now scored against
  the other years only.
- **Whole-year shifts in the null.** Wet and dry years were averaged away
  when synthetic years were stitched together, so the null was too narrow
  whenever such shifts exist. A separate annual-level step (Student-t
  prediction) now puts them back.
- **One decision for the p-value, the flagged months and the plot.**
  Months are flagged with a single-step max-T rule against the same
  critical value as the overall p-value, replacing Holm across months. The
  two could previously disagree. `plot_envelope` now draws that critical
  value as the band, so the line leaving the band means p ≤ alpha.
- **Sustained test** now requires a consistent direction (the average
  score over the window), so two high months followed by two low ones no
  longer count as a sustained departure.

### Defects fixed (each has a regression test)

- Months missing from the treatment period were still included in the
  null's maximum, making p-values too large (warning said otherwise).
- The trend test used raw yearly medians, so a year missing its summer
  looked like a trend step. It now uses seasonal anomalies.
- The default trend year ignored the treatment data's own date (a 2020
  period on a 2000–2014 trend was projected to 2015). It is now read from
  the dates; `year_index` is no longer required anywhere.
- Confidence mode logged the average instead of averaging the logs (a set
  of typical years scored p=0.005), and projected the trend to one year for
  all averaged years.
- A partially covered bin (e.g. 1 day of a 30-day bin) was compared with
  full bins; normal half-years were flagged about 30% of the time. Bins
  with less than 75% of their usual measurements now count as missing.
- Bins with almost no history could produce the smallest possible p-value.
- Text values such as `<0.5` were silently dropped (biasing records
  upward); they now raise an error.
- A treatment year that was also in the history was compared with itself;
  now an error.
- Sequential events in the same year were not detected; now an error.
- `sensitivity_grid` gave each row different random numbers (so rows
  differed by chance) and stopped entirely on one failing variant.
- Monitoring state could be reloaded against a different fit, silently
  changing the baseline part-way through a horizon.
- `<` and `<=` were used inconsistently for significance; now `p <= alpha`
  everywhere.

### Other changes

- Default bins are calendar months (`bin_days="month"`), which match
  monthly sampling. With fixed-length bins, a short leftover at the end of
  the year is merged into the last bin (30-day bins give 12 bins, not 13).
- New: `alternative="greater"/"less"`, `summarize()`,
  `significance_band()`, `critical_value()`, `fit.fingerprint()`,
  `TideConfig.min_bin_coverage`, `n_bins_for()`.
- The warning about high between-year variance was removed (the method
  now accounts for it); `between_year_var_frac` remains as a diagnostic
  of how hard whole-year changes are to detect.
- Documentation rewritten; `METHODS.md` added.
