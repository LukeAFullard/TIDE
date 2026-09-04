# Adversarial Review — tide_lite

Self-review of the Phase 2–5 code, done by trying to break it rather than
confirm it works: adversarial inputs run against the actual code, plus
line-by-line reasoning about silent-failure modes. Six issues found, all
fixed and re-verified; the full test suite (6/6) and example script were
re-run clean after every fix, not just at the end.

## Findings

### 1. Block resampling mixed variance across seasons (most consequential)

**What:** the null distribution was built by drawing sub-year blocks from
the *entire* pool regardless of which calendar position they originally
came from. A naturally volatile month's anomaly pattern could land on a
naturally stable month's slot in a synthetic draw, and vice versa.
**Why it matters:** this distorts the null distribution's shape in a
position-dependent way — not necessarily false positives, but not a
faithful null either, and not guaranteed conservative.
**Fix:** position-matched pooling — a block filling calendar positions
`[s, s+L)` is only ever drawn from historical years' own data at that
position (or its immediate neighbors — see finding 3). Random starting
phase per draw keeps it genuinely circular rather than a fixed tiling.
**Verified:** reasoned through directly against `_draw_synthetic_residual`'s
source before touching code; full test suite re-run clean after the fix.

### 2. `log_additive` silently corrupted bins on non-positive data

**What:** `np.log()` of a zero or negative value produces `NaN` with a
warning, not an error — and `nanmedian` silently drops NaNs. A variable
that can legitimately be zero or negative (redox potential, a
stage/anomaly measurement, some indices) would silently lose data from
affected bins with no visible error under default warning settings.
**Why it matters:** this is the dangerous kind of bug — no crash, no
obvious symptom, just a quietly-wrong fit.
**Fix:** hard validation before the log transform, on both historical and
treatment data, raising a clear `ValueError` naming the count of bad
values.
**Verified:** adversarially constructed data with genuinely negative
aggregated bins (`min raw value: -5.3`, 1407 negative days) — reproduced
the silent corruption first (`median_curve` had no NaN despite negative
inputs, `RuntimeWarning: invalid value encountered in log`), then
confirmed the fix raises instead.

### 3. Position-matching alone left too little resampling diversity

**What:** found while re-running the test suite after fix #1 —
`test_sequential_fwer_control` started failing (18.75% family-wise rate
vs. an expected ~5%). Traced it to plain `test_treatment` calls on the
same fits, ruling out a Holm-Bonferroni bug — the *individual* test
calibration was running hot (up to 28% on some fits). Root cause: exact
position-matching left exactly N=10 candidate blocks per calendar
position (one per historical year), down from ~130 pooled globally
before fix #1's correction.
**Why it matters:** this is a direct, quantified trade-off between fix #1
(correctness) and calibration stability — worth knowing even though it's
now fixed, because it's the kind of interaction that's easy to miss when
fixes are reviewed one at a time instead of re-tested together.
**Fix:** pool each position with its immediate ±1 neighbors (still no
distant-season mixing), restoring the candidate pool to ~30 per position.
**Verified:** re-ran the exact 8 fit-seeds that had shown 14.6% pooled /
28% worst-case rates — dropped to 5.6% pooled / 10% worst-case after the
fix. Full suite re-run clean.

### 4. `mode="prediction"` silently used only the first year of multi-year input

**What:** if `treatment_df` accidentally contained more than one year's
data, prediction mode silently tested only the first (alphabetically/
numerically first `block_year`) and ignored the rest.
**Why it matters:** a plausible, undetectable wrong-scope result — e.g.
from an off-by-one date-slicing mistake on the caller's side.
**Fix:** raise a clear error if `mode="prediction"` and the data spans
more than one year, naming the years found and pointing to
`mode="confidence"` as the alternative if that's what was intended.
**Verified:** adversarially constructed a 2-year `treatment_df`, confirmed
the pre-fix code silently produced a "valid-looking" p-value using only
year one; confirmed the fix raises instead.

### 5. Sequential/Cumulative could silently misapply the trend correction

**What:** `year_index` (used to extrapolate the historical trend forward
for a fair comparison) defaults to "immediately after the historical
period" when not given. Sequential and Cumulative test multiple periods
that are *not* assumed adjacent in time — if a caller forgot to set
`year_index` per event/window, every one of them would silently get the
same (likely wrong) trend offset whenever a trend correction is active.
**Fix:** `sequential_test`/`cumulative_test` now raise a clear error
listing which events/windows are missing `year_index`, but only when the
fit actually has a trend correction applied (no cost when there's no
trend to misapply).
**Verified:** constructed a fit with a real injected trend
(`slope=0.5/year`, Mann-Kendall p≈8e-5), confirmed Sequential silently
ran before the fix and now raises.

### 6. Two smaller robustness gaps

- **`effect_size` vs. what the p-value tested:** `effect_size` was never
  trend-adjusted, while the p-value always is when a trend correction is
  active — correct on its own, but a reader could see a small raw effect
  size next to a very significant p-value and reasonably wonder why.
  Fixed by reporting both `effect_size` (raw) and
  `effect_size_trend_adjusted` (matches the test) on every result;
  they're equal when no trend was applied.
- **MAD=0 fallback used a hardcoded `1.0`:** the last-resort fallback for
  a bin where ≥half of historical years tie exactly used a fixed constant,
  a unit mismatch for arbitrarily-scaled variables (a variable in the
  thousands, or in the thousandths). Fixed to scale off the data's own
  level instead.

## Considered, not changed

- **`rng` defaults to fresh (unseeded) entropy everywhere.** Every
  function now also accepts a plain int as a convenience (`rng=42`), but
  the *default* is still non-reproducible. Left as-is: forcing a fixed
  default seed would make every unseeded call silently deterministic in a
  way that could mask real variability, which seems worse. Recommendation
  stands — pass an explicit seed for anything you intend to report.
- **Years that pass the completeness filter but contribute few/no blocks
  to the pool** (scattered rather than concentrated missingness). Not
  fixed — `min_blocks_per_position` on `TideFit` at least makes this
  visible rather than silent; a real fix (per-position completeness
  weighting) felt like more machinery than the current data warrants.
  Revisit if `min_blocks_per_position` ever comes out surprisingly low on
  real data.

## Re-verification after all fixes

- `tests/test_engine.py`: 6/6 passing (re-run after every individual fix,
  not just once at the end)
- `examples/run_example.py`: runs end-to-end, plot regenerates correctly
- Performance: unaffected — 0.065s per `test_treatment` call at the
  default `n_bootstrap=2000`, n_years=10 (was 0.054s pre-review; the
  position-matched, windowed pooling is doing more work but it's not
  where the cost is)
- Edge cases re-checked and still behave correctly: n_years=3 (the
  documented minimum), string dates, duplicate timestamps, empty/all-NaN
  input (already handled pre-review)
