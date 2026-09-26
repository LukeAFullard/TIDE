"""
tide_lite core engine.

Question: was a treatment period unusual compared with the same site's own
historical "normal" years, allowing for the seasonal cycle, a long-term
trend, and ordinary year-to-year variation, without assuming the data
follow a bell curve?

THE METHOD, EVERY STEP IN ORDER (METHODS.md section 2 is the
plain-language version; keep the two in step):

 1. Bin. Each calendar year is cut into seasonal bins (calendar months by
    default, or fixed day-of-year blocks) and each bin is summarised by the
    median (or mean) of the values in it.
 2. Screen. A bin in a given year counts as missing if it holds fewer than
    `min_bin_coverage` (75%) of that bin's usual number of measurements
    (applied to history and treatment alike). Then bins the record rarely
    covers are dropped, then years missing more than `max_missing_frac` of
    the remaining bins, then any bin left with fewer than 3 years of data.
    Everything dropped is reported.
 3. Transform (optional). detrend_mode="log_additive" works on natural
    logs, so changes are proportional rather than absolute.
 4. Trend. Each year's typical anomaly (median over bins of its departure
    from that bin's historical median) is tested for a monotonic trend with
    Mann-Kendall. If p < mk_alpha, the Sen slope (per calendar year) is
    removed from every year.
 5. Seasonal baseline. Per bin: the median across years (typical level) and
    1.4826 x the median absolute deviation (typical spread, "MAD").
 6. Leave-one-year-out departures. Each historical year is compared with
    the median and MAD of the OTHER years, exactly as the treatment year is
    compared with all of them, so these departures carry the same
    estimation noise a genuinely new year does.
 7. Each historical year's departures are split into an annual level (their
    mean across bins, in data units: "was this a high year overall?") and
    a within-year pattern (the rest, divided by the leave-one-out spread).
 8. Null distribution. Each synthetic normal year = a within-year pattern
    stitched together from blocks of consecutive bins taken from historical
    years (a block always comes from exactly the same bins of the year, so
    every month is represented once; the block seams fall at a random
    place) + one annual level drawn from a Student-t prediction
    distribution fitted to the historical annual levels, divided by each
    bin's spread. n_bootstrap synthetic years are drawn.
 9. The treatment period is scored against the step-5 baseline, with the
    step-4 trend extrapolated to its own calendar year.
10. Test statistic: the largest bin score (largest absolute score for a
    two-sided test). p-value = (number of synthetic years with a statistic
    at least as large + 1) / (n_bootstrap + 1).
11. Months: a bin is flagged when its own score exceeds the same critical
    value (single-step max-T, Westfall & Young 1993). So the period is
    significant if and only if at least one bin is flagged, and the plotted
    significance band is exactly that critical value.

WHY EACH NON-OBVIOUS STEP IS THERE (each was measured, not assumed; the
simulation study is in tests/calibration_study.py):

- Leave-one-year-out (step 6). Scoring historical years against a median
  and MAD that they themselves helped estimate makes them look more
  ordinary than a genuinely new year. With 10 historical years that alone
  made the test flag 2-3x too often at strict thresholds.
- Annual level (steps 7-8). Stitching blocks from different years averages
  away whole-year shifts (a wet year that sits high all year), so without
  step 7 the null is too narrow whenever such shifts exist, the norm in
  hydrology. Resampling the N observed annual levels instead cannot produce
  a level more extreme than any seen, although a new year exceeds all N
  about 2/(N+1) of the time. The Student-t prediction distribution is the
  textbook answer to "how far can the next value be from N past values"
  (it widens for small N and for extrapolating a trend). This is the one
  distributional assumption in the method, and it applies only to the
  single annual-level number per year.
- Sub-year blocks (step 8). Resampling whole years gives only N possible
  synthetic years, so the smallest attainable p-value would be 1/(N+1)
  (0.09 for N=10). Blocks of consecutive bins keep the within-year
  persistence (autocorrelation) that independent bins would destroy. A
  block is taken from the same bins it fills (pool_window_radius=0):
  borrowing from neighbouring bins could not reach past January or
  December, so those months were under-represented in the null and false
  alarms rose in every simulated condition.
- Max statistic + single-step max-T (steps 10-11). One number for the whole
  period, with no multiple-comparisons inflation across bins, and it makes
  the per-bin flags and the plot agree with the p-value by construction.
- Median/MAD (step 5). Robust to skew and outliers typical of concentration
  data; the 1.4826 factor only puts MAD on a standard-deviation scale.

Components and references: Mann (1945), Kendall (1975), Sen (1968) for the
trend; Kunsch (1989) and Politis & Romano (1992) for block bootstrap, with
blocks matched to the same season as in the seasonal block bootstraps of
Politis (2001) and Dudek, Leskow, Paparoditis & Politis (2014); Westfall &
Young (1993) for max-T; Davison & Hinkley (1997) and Phipson & Smyth
(2010) for the (b+1)/(B+1) Monte Carlo p-value; Hahn & Meeker (1991) for
prediction intervals; Rousseeuw & Croux (1993) for MAD; Holm (1979).
"""

from __future__ import annotations

from dataclasses import dataclass, asdict, field
import hashlib
import json
import warnings

import numpy as np
import pandas as pd
from scipy.stats import norm


# ---------------------------------------------------------------------------
# Seasonal bins
# ---------------------------------------------------------------------------

def to_day_of_year(dates: pd.Series) -> pd.Series:
    """Map dates to a 1-365 day-of-year axis, the same in leap and non-leap
    years. Feb 29 becomes NaN (dropped: at most one day every four years);
    in a leap year every date from Mar 1 shifts down by 1, so Dec 31 is
    always day 365. Only used for fixed-length (integer bin_days) bins;
    calendar-month bins need no alignment.
    """
    dates = pd.to_datetime(pd.Series(dates))
    doy = dates.dt.dayofyear.astype("float64").copy()
    is_leap = dates.dt.is_leap_year
    is_feb29 = is_leap & (dates.dt.month == 2) & (dates.dt.day == 29)
    shift = is_leap & (doy >= 61) & (~is_feb29)
    doy = doy.where(~shift, doy - 1)
    doy = doy.where(~is_feb29, np.nan)
    return doy


def n_bins_for(bin_days) -> int:
    """Number of bins per year. "month" -> 12. For an integer bin length,
    365 // bin_days full bins, plus one more only if the leftover days are
    at least half a bin; a shorter leftover is merged into the last bin
    (30-day bins -> 12 bins, the last 35 days long; 7-day bins -> 52).
    A 5-day or 1-day stub bin would be far noisier than the others.
    """
    if bin_days == "month":
        return 12
    n_full, rem = divmod(365, int(bin_days))
    if rem == 0 or (n_full >= 1 and rem < bin_days / 2):
        return max(n_full, 1)
    return n_full + 1


def assign_bins(day_of_year: pd.Series, bin_days: int) -> pd.Series:
    """Bin number (1-based) for each day of year, with a short leftover
    merged into the last bin (see n_bins_for)."""
    raw = ((day_of_year - 1) // bin_days).astype("Int64") + 1
    return raw.clip(upper=n_bins_for(bin_days))


# ---------------------------------------------------------------------------
# Trend detection (Mann-Kendall) and magnitude (Sen's slope)
# ---------------------------------------------------------------------------

def mann_kendall_test(values) -> tuple[float, float]:
    """Two-sided Mann-Kendall trend test (Mann 1945; Kendall 1975) with the
    tie-corrected variance and continuity correction. Values must be in
    time order. Returns (S statistic, p-value). NaNs dropped first."""
    x = np.asarray(values, dtype="float64")
    x = x[~np.isnan(x)]
    n = len(x)
    if n < 4:
        return 0.0, 1.0
    s = sum(np.sum(np.sign(x[i + 1:] - x[i])) for i in range(n - 1))
    _, counts = np.unique(x, return_counts=True)
    tie_term = np.sum(counts * (counts - 1) * (2 * counts + 5))
    var_s = (n * (n - 1) * (2 * n + 5) - tie_term) / 18.0
    if var_s <= 0:
        return float(s), 1.0
    if s > 0:
        z = (s - 1) / np.sqrt(var_s)
    elif s < 0:
        z = (s + 1) / np.sqrt(var_s)
    else:
        z = 0.0
    p = 2 * (1 - norm.cdf(abs(z)))
    return float(s), float(p)


def sens_slope(values, times=None) -> float:
    """Sen's slope (Sen 1968): median of all pairwise slopes, per unit of
    `times`. Pass real calendar times whenever observations are unevenly
    spaced (a gap in the record): the positional index silently rescales
    the slope (years 2000-04 + 2010-14 drifting 0.5/yr give ~0.87/yr
    positionally, 0.5/yr with real times)."""
    x = np.asarray(values, dtype="float64")
    n = len(x)
    t = np.arange(n, dtype="float64") if times is None else np.asarray(times, dtype="float64")
    if len(t) != n:
        raise ValueError(f"times has length {len(t)} but values has length {n}.")
    slopes = [
        (x[j] - x[i]) / (t[j] - t[i])
        for i in range(n - 1)
        for j in range(i + 1, n)
        if not (np.isnan(x[i]) or np.isnan(x[j]) or t[j] == t[i])
    ]
    return float(np.median(slopes)) if slopes else 0.0


def holm_bonferroni(p_values, alpha: float = 0.05) -> np.ndarray:
    """Holm (1979) step-down correction: controls the chance of ANY false
    alarm across len(p_values) tests at `alpha`, under any dependence.
    Sort p ascending; compare p(1) to alpha/n, p(2) to alpha/(n-1), ...;
    stop at the first that fails. Used across events in sequential_test.
    Returns a boolean array (reject) in the input order."""
    p_values = np.asarray(p_values, dtype="float64")
    if not 0 < alpha < 1:
        raise ValueError(f"alpha must be strictly between 0 and 1, got {alpha}.")
    if np.any(np.isnan(p_values)):
        raise ValueError(
            "holm_bonferroni received NaN p-values. Drop them first: a NaN is "
            "missing evidence, and sorting it into the sequence would change "
            "every other test's threshold."
        )
    n = len(p_values)
    order = np.argsort(p_values, kind="stable")
    reject_sorted = np.zeros(n, dtype=bool)
    for i, p in enumerate(p_values[order]):
        if _le(p, alpha / (n - i)):
            reject_sorted[i] = True
        else:
            break
    reject = np.zeros(n, dtype=bool)
    reject[order] = reject_sorted
    return reject


def _le(p, alpha) -> bool:
    """The one significance rule used everywhere: p <= alpha. A Monte Carlo
    p-value is (b+1)/(B+1), so `<=` is the rule that gives an exact level-
    alpha test; the tiny tolerance only absorbs floating-point rounding so
    that e.g. 100/2000 counts as <= 0.05 in every code path."""
    return bool(p <= alpha * (1 + 1e-9))


# ---------------------------------------------------------------------------
# rng convenience: accept None, an int seed, or an existing Generator.
# ---------------------------------------------------------------------------

def _resolve_rng(rng) -> np.random.Generator:
    if isinstance(rng, np.random.Generator):
        return rng
    return np.random.default_rng(rng)  # None -> fresh entropy; 42 -> seeded


# ---------------------------------------------------------------------------
# Config and result containers
# ---------------------------------------------------------------------------

_VALID_DETREND_MODES = {"additive", "log_additive", "none"}
_VALID_AGG = {"median", "mean"}
_VALID_MODES = {"prediction", "confidence"}
_VALID_ALTERNATIVES = {"two-sided", "greater", "less"}


@dataclass
class TideConfig:
    """Every analysis setting, with defaults. Record these with any result."""
    bin_days: "int | str" = "month"  # "month" = calendar months; or an integer
                                     # number of days (e.g. 7 = weekly)
    agg: str = "median"              # summary of the values within one bin
    detrend_mode: str = "additive"   # "additive": trend in original units;
                                     # "log_additive": logs first (proportional
                                     # changes; needs values > 0); "none": no
                                     # trend removal, no log
    mk_alpha: float = 0.10           # remove a trend if Mann-Kendall p < this
    max_missing_frac: float = 0.5    # a year may miss up to this share of bins
    min_bin_coverage: float = 0.75   # a bin in a given year counts only if it has at
                                     # least this share of that bin's usual number of
                                     # measurements (a 3-day "June" is not comparable
                                     # with a 30-day one)
    n_bootstrap: int = 2000          # synthetic normal years in the null
    block_length_bins: "int | None" = None  # None = auto: max(2, n_bins // 4)
    pool_window_radius: int = 0      # blocks may come from +/- this many bins
                                     # away from the slot they fill. Keep 0:
                                     # with 1, the first and last bins of the
                                     # year were under-represented in the null
                                     # (0.6x) and false alarms rose

    def __post_init__(self):
        """Reject misspelled or impossible settings loudly rather than
        silently running a different analysis."""
        if self.detrend_mode not in _VALID_DETREND_MODES:
            raise ValueError(
                f"detrend_mode must be one of {sorted(_VALID_DETREND_MODES)}, "
                f"got {self.detrend_mode!r}."
            )
        if self.agg not in _VALID_AGG:
            raise ValueError(f"agg must be one of {sorted(_VALID_AGG)}, got {self.agg!r}.")
        if self.bin_days != "month":
            if isinstance(self.bin_days, bool) or not isinstance(self.bin_days, (int, np.integer)):
                raise ValueError(f"bin_days must be 'month' or an integer, got {self.bin_days!r}.")
            if not 1 <= self.bin_days <= 365:
                raise ValueError(f"bin_days must be between 1 and 365, got {self.bin_days}.")
        if not 0 < self.mk_alpha < 1:
            raise ValueError(f"mk_alpha must be strictly between 0 and 1, got {self.mk_alpha}.")
        if not 0 <= self.max_missing_frac < 1:
            raise ValueError(f"max_missing_frac must be in [0, 1), got {self.max_missing_frac}.")
        if self.n_bootstrap < 1:
            raise ValueError(f"n_bootstrap must be >= 1, got {self.n_bootstrap}.")
        if self.block_length_bins is not None and self.block_length_bins < 1:
            raise ValueError(
                f"block_length_bins must be >= 1 or None (auto), got {self.block_length_bins}."
            )
        if not 0 <= self.min_bin_coverage <= 1:
            raise ValueError(f"min_bin_coverage must be in [0, 1], got {self.min_bin_coverage}.")
        if self.pool_window_radius < 0:
            raise ValueError(f"pool_window_radius must be >= 0, got {self.pool_window_radius}.")

    @property
    def min_attainable_p(self) -> float:
        """Smallest p-value n_bootstrap draws can produce: 1 / (n_bootstrap + 1)."""
        return 1.0 / (self.n_bootstrap + 1)


@dataclass
class TideFit:
    config: TideConfig
    bins: np.ndarray                 # bin labels kept (months 1-12 by default)
    median_curve: np.ndarray         # per-bin median, detrended, working units
    mad: np.ndarray                  # per-bin spread (1.4826 x MAD), working units
    residual_matrix: np.ndarray      # (years, bins) detrended value - median_curve
    historical_matrix: np.ndarray    # (years, bins) detrended values, working units
    loo_scores: np.ndarray           # (years, bins) leave-one-year-out scores (step 6)
    annual_levels: np.ndarray        # (years,) each year's mean leave-one-out
                                     # departure, working units (step 7)
    within_year: np.ndarray          # (years, bins) the rest of each year's departure,
                                     # divided by its leave-one-out spread
    level_mean: float                # Student-t prediction for the annual level:
    level_sd: float                  #   mean, SD and degrees of freedom of the
    level_df: int                    #   historical annual levels (working units)
    block_length_bins: int
    min_blocks_per_position: int     # fewest full-length blocks available for any slot
    years_used: list
    years_dropped: list              # [(year, reason)]
    bins_dropped: list               # [(bin, reason)]
    trend: dict
    date_col: str
    value_col: str
    median_curve_raw: np.ndarray     # per-bin median BEFORE detrending (the
                                     # record's own typical level; effect_size
                                     # is measured against this)
    between_year_var_frac: float     # share of score variance that is whole-year
                                     # level (0-1); high = wet/dry-year shifts
                                     # dominate, so uniform shifts are harder to detect
    year_offsets: np.ndarray = field(default=None, repr=False)
    typical_counts: np.ndarray = field(default=None, repr=False)  # usual measurements per bin

    @property
    def first_year(self) -> int:
        """Calendar year the trend axis is measured from."""
        return int(self.years_used[0])

    def year_index_for(self, calendar_year: int) -> int:
        """Calendar year -> offset from the first historical year, the
        `treatment_year_index` the trend is extrapolated to. Only needed to
        override the default, which is read from the treatment data's dates."""
        return int(calendar_year) - self.first_year

    def fingerprint(self) -> str:
        """Short hash of the settings, years and fitted baseline. Two fits
        with the same fingerprint give the same null for the same seed.
        Record it with a result; MonitoringSeries checks it on reload."""
        payload = {
            "config": {k: str(v) for k, v in asdict(self.config).items()},
            "years": [int(y) for y in self.years_used],
            "bins": [int(b) for b in self.bins],
            "median": np.round(self.median_curve, 9).tolist(),
            "mad": np.round(self.mad, 9).tolist(),
            "slope": round(float(self.trend["slope_per_year"]), 12),
            "cols": [self.date_col, self.value_col],
        }
        return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()[:16]


@dataclass
class SustainedResult:
    """Was there a run of `run_length` consecutive bins that together
    departed in one direction? Statistic: the largest window-average score
    (absolute value of the average for a two-sided test, so a window must
    lean consistently one way). Windows never wrap from the last bin to the
    first, and windows touching a bin with no data are skipped."""
    run_length: int
    p_value: float
    test_statistic: float
    window_start_bin: int
    window_bins: np.ndarray


@dataclass
class TideResult:
    p_value: float                   # whole-period p-value
    test_statistic: float            # largest bin score (see alternative)
    effect_size: float               # median over bins of (treatment - historical
                                     # typical level), original units, no trend
    effect_size_trend_adjusted: float  # same, against the trend-projected level
                                     # for the treatment year: the comparison
                                     # the p-value makes
    bin_scores: np.ndarray           # signed score per bin: (value - median)/spread
    bin_p_adjusted: np.ndarray       # per-bin p-value, max-T adjusted: compare
                                     # directly to alpha (bin_significant)
    bin_p_values: np.ndarray         # per-bin p-value, UNadjusted (descriptive)
    sustained: dict                  # {run_length: SustainedResult}
    mode: str
    alternative: str
    treatment_curve: np.ndarray      # binned treatment values, original units
    bins: np.ndarray
    n_reference: int                 # synthetic years in the null
    null_max_stats: np.ndarray = field(default=None, repr=False)
    treatment_year_index: float = 0.0  # year offset the trend was projected to
    treatment_years: list = None     # calendar year(s) tested
    missing_bins: np.ndarray = None  # True where the treatment has no data
    seed: "int | None" = None        # the int seed passed as rng, if any


# ---------------------------------------------------------------------------
# Binning and screening
# ---------------------------------------------------------------------------

def _bin_and_pivot(df, date_col, value_col, bin_days, agg, with_counts=False):
    """(calendar year x bin) table of per-bin summaries, all bins as columns,
    and optionally the matching table of how many measurements each holds."""
    for col in (date_col, value_col):
        if col not in df.columns:
            raise ValueError(f"Column {col!r} not found. Columns are: {list(df.columns)}.")
    d = df[[date_col, value_col]].copy()
    raw_vals = d[value_col]
    vals = pd.to_numeric(raw_vals, errors="coerce")
    bad = raw_vals.notna() & vals.isna()
    if bad.any():
        examples = raw_vals[bad].astype(str).unique()[:5].tolist()
        raise ValueError(
            f"{int(bad.sum())} value(s) in {value_col!r} are not numbers, e.g. "
            f"{examples}. Censored results such as '<0.5' must be converted to "
            f"numbers first, the same way for every year (for example the "
            f"detection limit, or half of it). Dropping them silently would bias "
            f"the record upwards."
        )
    d[value_col] = vals
    dates = pd.to_datetime(d[date_col])
    if bin_days == "month":
        d["bin"] = dates.dt.month
    else:
        d["bin"] = assign_bins(to_day_of_year(dates), int(bin_days))
    d["block_year"] = dates.dt.year
    d = d.dropna(subset=[value_col, "bin"])
    cols = range(1, n_bins_for(bin_days) + 1)
    pivot = d.pivot_table(index="block_year", columns="bin", values=value_col, aggfunc=agg)
    pivot = pivot.reindex(columns=cols)
    pivot.columns = pivot.columns.astype(int)
    if not with_counts:
        return pivot
    counts = d.pivot_table(index="block_year", columns="bin", values=value_col, aggfunc="count")
    counts = counts.reindex(index=pivot.index, columns=cols).fillna(0)
    counts.columns = counts.columns.astype(int)
    return pivot, counts


def _mask_thin_bins(pivot, counts, typical, min_coverage):
    """Blank year-bin cells holding fewer than min_coverage x the usual number
    of measurements for that bin. Returns (masked pivot, number blanked)."""
    thin = (counts > 0) & (counts < min_coverage * typical)
    return pivot.mask(thin), int(thin.to_numpy().sum())


def _check_positive_for_log(values, context: str):
    arr = np.asarray(values, dtype="float64")
    n_bad = int(np.nansum(arr <= 0))
    if n_bad > 0:
        raise ValueError(
            f"detrend_mode='log_additive' needs strictly positive values, but "
            f"{context} has {n_bad} value(s) <= 0. Use detrend_mode='additive', "
            f"or fix the data first."
        )


def _spread(dev: np.ndarray, fallback: "np.ndarray | None" = None) -> tuple[np.ndarray, np.ndarray]:
    """1.4826 x median absolute deviation per column. A zero spread (half or
    more of the values identical, typical of non-detects reported at a
    detection limit) is replaced by `fallback` if given, else by the median
    positive spread of the other columns. Returns (spread, was_zero)."""
    s = 1.4826 * np.nanmedian(np.abs(dev), axis=0)
    zero = ~(s > 0)
    if np.any(zero):
        if fallback is not None:
            s = np.where(zero, fallback, s)
        elif np.any(~zero):
            s = np.where(zero, np.nanmedian(s[~zero]), s)
        else:
            s = np.full_like(s, 1e-9)
    return s, zero


# ---------------------------------------------------------------------------
# fit_historical
# ---------------------------------------------------------------------------

def fit_historical(df: pd.DataFrame, date_col: str, value_col: str,
                   config: "TideConfig | None" = None) -> TideFit:
    """Steps 1-8 of the method (module docstring) on the historical
    ("normal") years: bin, screen, transform, detrend, baseline, and the
    ingredients of the null distribution. Deterministic: no randomness."""
    config = config or TideConfig()
    pivot, counts = _bin_and_pivot(df, date_col, value_col, config.bin_days, config.agg,
                                   with_counts=True)
    if len(pivot) == 0:
        raise ValueError("No usable rows in the historical data.")
    typical_counts = counts.where(counts > 0).median(axis=0)
    pivot, n_thin = _mask_thin_bins(pivot, counts, typical_counts, config.min_bin_coverage)
    if n_thin:
        warnings.warn(
            f"{n_thin} year-bin cell(s) had fewer than {config.min_bin_coverage:.0%} of the "
            f"usual number of measurements for that bin and were treated as missing.",
            UserWarning, stacklevel=2,
        )

    # Step 2: screening. Bins first (so a sampling programme that never
    # samples some months does not get every year thrown out), then years,
    # then any bin too thin to estimate a median and spread from.
    keep_share = 1.0 - config.max_missing_frac
    bins_dropped = []
    cover = pivot.notna().mean(axis=0)
    sparse = cover.index[cover < keep_share]
    bins_dropped += [(int(b), f"data in {cover[b]:.0%} of years") for b in sparse]
    pivot = pivot.drop(columns=sparse)
    if pivot.shape[1] == 0:
        raise ValueError("No bin has data in enough years. Check the date and value columns.")

    missing_frac = pivot.isna().mean(axis=1)
    keep = missing_frac <= config.max_missing_frac
    years_dropped = [(int(y), f"{missing_frac[y]:.0%} of bins missing")
                     for y in pivot.index[~keep]]
    pivot = pivot.loc[keep]
    if len(pivot) < 3:
        raise ValueError(
            f"Only {len(pivot)} historical years survived the completeness filter "
            f"(need >= 3, and >= 10 for useful results). Lower max_missing_frac, "
            f"use coarser bins, or check the data."
        )
    n_per_bin = pivot.notna().sum(axis=0)
    thin = n_per_bin.index[n_per_bin < 3]
    bins_dropped += [(int(b), f"only {n_per_bin[b]} year(s) of data") for b in thin]
    pivot = pivot.drop(columns=thin)
    if pivot.shape[1] == 0:
        raise ValueError("No bin has at least 3 years of historical data.")
    if bins_dropped:
        warnings.warn(
            f"Bins left out of the analysis (too little historical data): "
            f"{bins_dropped}. The test covers only the remaining bins.",
            UserWarning, stacklevel=2,
        )

    years_used = [int(y) for y in pivot.index]
    W = pivot.to_numpy(dtype="float64")
    bins = pivot.columns.to_numpy(dtype=int)
    N, nb = W.shape
    if N < 10:
        warnings.warn(
            f"Only {N} historical years. The test stays valid (it widens to allow "
            f"for the small sample) but can only detect large changes. Run "
            f"estimate_power before relying on a non-significant result.",
            UserWarning, stacklevel=2,
        )

    # Step 3: optional log transform.
    use_log = config.detrend_mode == "log_additive"
    if use_log:
        _check_positive_for_log(W, "the historical data")
        W = np.log(W)

    # Step 4: trend on deseasonalised annual anomalies, on a calendar-year axis.
    # (A year's raw median would move with WHICH months it has data for, so a
    # year missing its summer would look like a step in the trend.)
    year_offsets = np.asarray(years_used, dtype="float64") - float(years_used[0])
    median_curve_raw = np.nanmedian(W, axis=0)
    annual_anomaly = np.nanmedian(W - median_curve_raw, axis=1)
    _, mk_p = mann_kendall_test(annual_anomaly)
    slope, applied = 0.0, False
    if config.detrend_mode != "none" and mk_p < config.mk_alpha:
        slope = sens_slope(annual_anomaly, times=year_offsets)
        W = W - slope * year_offsets[:, None]
        applied = True

    # Step 5: seasonal baseline.
    median_curve = np.nanmedian(W, axis=0)
    residual_matrix = W - median_curve
    mad, zero_mad = _spread(residual_matrix)
    if np.any(zero_mad):
        warnings.warn(
            f"In bins {bins[zero_mad].tolist()} at least half the historical values "
            f"are identical (common with non-detects reported at a detection "
            f"limit), so their spread is zero. The median spread of the other "
            f"bins is used there instead; treat results in those bins as approximate.",
            UserWarning, stacklevel=2,
        )

    # Step 6: leave-one-year-out departures and spreads.
    loo_dev = np.full_like(W, np.nan)
    loo_spread = np.full_like(W, np.nan)
    for y in range(N):
        other = np.delete(W, y, axis=0)
        m_o = np.nanmedian(other, axis=0)
        loo_spread[y], _ = _spread(other - m_o, fallback=mad)
        loo_dev[y] = W[y] - m_o
    loo = loo_dev / loo_spread

    # Step 7: annual level (in working units, because a wet year shifts every
    # month by a similar AMOUNT, which is a larger score in a month with a
    # small spread) + within-year pattern (scored against the spread).
    annual_levels = np.nanmean(loo_dev, axis=1)
    within_year = (loo_dev - annual_levels[:, None]) / loo_spread
    ddof = 2 if applied else 1          # a fitted slope costs one more degree of freedom
    level_df = N - ddof
    level_sd = float(np.std(annual_levels, ddof=ddof))
    total_var = float(np.nanvar(loo_dev))
    between_year_var_frac = (
        float(np.clip(np.var(annual_levels) / total_var, 0.0, 1.0)) if total_var > 0 else 0.0
    )

    block_length_bins = min(config.block_length_bins or max(2, nb // 4), nb)
    min_blocks = min(
        len(_block_candidates(within_year, s, block_length_bins, config.pool_window_radius))
        for s in range(nb - block_length_bins + 1)
    )

    return TideFit(
        config=config, bins=bins, median_curve=median_curve, mad=mad,
        residual_matrix=residual_matrix, historical_matrix=W,
        loo_scores=loo, annual_levels=annual_levels, within_year=within_year,
        level_mean=float(np.mean(annual_levels)), level_sd=level_sd, level_df=int(level_df),
        block_length_bins=int(block_length_bins), min_blocks_per_position=int(min_blocks),
        years_used=years_used, years_dropped=years_dropped, bins_dropped=bins_dropped,
        trend={"mk_p": mk_p, "slope_per_year": slope, "applied": applied,
               "log_space": use_log, "first_year": int(years_used[0])},
        date_col=date_col, value_col=value_col, median_curve_raw=median_curve_raw,
        between_year_var_frac=between_year_var_frac, year_offsets=year_offsets,
        typical_counts=typical_counts.reindex(bins).to_numpy(dtype="float64"),
    )


# ---------------------------------------------------------------------------
# Null distribution (step 8)
# ---------------------------------------------------------------------------

def _block_candidates(U: np.ndarray, start: int, length: int, radius: int) -> np.ndarray:
    """All gap-free blocks of `length` consecutive bins, from any historical
    year, starting within +/- radius bins of `start` (never wrapping past
    either end of the year). Shape (n_candidates, length)."""
    n_years, nb = U.shape
    out = []
    for off in range(-radius, radius + 1):
        s = start + off
        if s < 0 or s + length > nb:
            continue
        blocks = U[:, s:s + length]
        out.append(blocks[~np.isnan(blocks).any(axis=1)])
    return np.concatenate(out) if out else np.empty((0, length))


def _stitch_within_year(fit: TideFit, n_draws: int, rng: np.random.Generator) -> np.ndarray:
    """(n_draws, n_bins) synthetic within-year patterns. Seams fall every
    block_length_bins bins starting from a random phase, so no bin is
    always at a seam; each block comes from the same time of year."""
    U = fit.within_year
    nb = U.shape[1]
    L = fit.block_length_bins
    radius = fit.config.pool_window_radius
    out = np.empty((n_draws, nb))
    phases = rng.integers(0, L, size=n_draws)
    cache = {}
    for phase in range(L):
        idx = np.flatnonzero(phases == phase)
        if idx.size == 0:
            continue
        seams = sorted({0, nb, *range(phase, nb, L)})
        for s, e in zip(seams[:-1], seams[1:]):
            key = (s, e - s)
            if key not in cache:
                cache[key] = _block_candidates(U, s, e - s, radius)
            cands = cache[key]
            if len(cands):
                out[idx, s:e] = cands[rng.integers(0, len(cands), size=idx.size)]
            else:
                # No gap-free block of this length here (patchy history):
                # fill bin by bin from the same neighbourhood instead.
                for p in range(s, e):
                    lo, hi = max(0, p - radius), min(nb, p + radius + 1)
                    vals = U[:, lo:hi].ravel()
                    vals = vals[~np.isnan(vals)]
                    out[idx, p] = vals[rng.integers(0, len(vals), size=idx.size)]
    return out


def _level_scale(fit: TideFit, year_indices) -> float:
    """Student-t prediction scale for the (average) annual level of
    len(year_indices) new years: sqrt(1/k + 1/N [+ trend extrapolation])."""
    k = len(year_indices)
    N = len(fit.years_used)
    var_factor = 1.0 / k + 1.0 / N
    if fit.trend["applied"]:
        x = fit.year_offsets
        sxx = float(np.sum((x - x.mean()) ** 2))
        var_factor += (float(np.mean(year_indices)) - x.mean()) ** 2 / sxx
    return float(np.sqrt(var_factor))


def _null_scores(fit: TideFit, n_draws: int, rng: np.random.Generator,
                 year_indices) -> np.ndarray:
    """(n_draws, n_bins) bin scores of synthetic normal periods: the
    average of len(year_indices) stitched within-year patterns plus one
    annual level from the Student-t prediction distribution."""
    k = len(year_indices)
    within = _stitch_within_year(fit, n_draws, rng)
    if k > 1:
        for _ in range(k - 1):
            within += _stitch_within_year(fit, n_draws, rng)
        within /= k
        # Each historical within-year pattern already includes the error of
        # the median it was compared with (leave-one-year-out), but that
        # error is the SAME for every year being averaged, so it must not
        # shrink by 1/k. Restore it using the large-sample variance of a
        # median, pi/2 x sigma^2 / N:  var(avg) = sigma^2 (1/k + v), v = pi/(2N).
        v = np.pi / (2 * len(fit.years_used))
        within *= np.sqrt((1 + k * v) / (1 + v))
    t = rng.standard_t(fit.level_df, size=n_draws)
    level = fit.level_mean + fit.level_sd * _level_scale(fit, year_indices) * t
    return within + level[:, None] / fit.mad


def _directional(z: np.ndarray, alternative: str) -> np.ndarray:
    if alternative == "two-sided":
        return np.abs(z)
    return z if alternative == "greater" else -z


def _window_scores(z: np.ndarray, run_length: int, alternative: str,
                   bins: "np.ndarray | None" = None) -> np.ndarray:
    """Directional score of every run of run_length consecutive bins (the
    mean score, then |.| for two-sided). Works on 1-D or (draws, bins);
    windows touching a NaN come back as -inf so they never win a max, and
    so do windows that jump over a bin dropped from the fit (given `bins`,
    the bin labels): those bins are not consecutive."""
    z = np.atleast_2d(z)
    n = z.shape[1]
    if run_length < 1:
        raise ValueError(f"run_length must be >= 1, got {run_length}.")
    if run_length > n:
        raise ValueError(f"run_length ({run_length}) exceeds the number of bins ({n}).")
    means = np.lib.stride_tricks.sliding_window_view(z, run_length, axis=1).mean(axis=-1)
    scores = _directional(means, alternative)
    scores = np.where(np.isnan(scores), -np.inf, scores)
    if bins is not None:
        b = np.asarray(bins)
        gap = (b[run_length - 1:] - b[:n - run_length + 1]) != run_length - 1
        scores[:, gap] = -np.inf
    return scores


def _max_window_mean(values: np.ndarray, run_length: int) -> tuple[float, int]:
    """Largest mean over runs of run_length consecutive values (no
    wrap-around; runs containing NaN excluded). Returns (value, start)."""
    scores = _window_scores(np.asarray(values, dtype="float64"), run_length, "greater")[0]
    if not np.any(np.isfinite(scores)):
        raise ValueError(
            f"No run of {run_length} consecutive bins is free of missing data in "
            f"the treatment period. Use a shorter run_length or more complete data."
        )
    best = int(np.argmax(scores))
    return float(scores[best]), best


# ---------------------------------------------------------------------------
# test_treatment
# ---------------------------------------------------------------------------

def test_treatment(fit: TideFit, treatment_df: pd.DataFrame,
                   treatment_year_index: "int | None" = None,
                   mode: str = "prediction",
                   run_lengths: "list[int] | None" = None,
                   alternative: str = "two-sided",
                   rng=None) -> TideResult:
    """Steps 9-11: test one treatment period against the historical fit.

    mode="prediction" (default): exactly one calendar year of data (a part
      year is fine up to max_missing_frac). Raises if it spans more.
    mode="confidence": several calendar years, averaged bin by bin, tested
      against the average of the same number of synthetic years. A bin
      missing in any of those years is left out.
    alternative: "two-sided" (default: unusually high OR low), "greater"
      (unusually high only) or "less". Choose before looking at the data.
    run_lengths: optional list, e.g. [3], to also test for a sustained
      departure of that many consecutive bins (result.sustained).
    treatment_year_index: where on the trend line the period sits. Default:
      read from the data's own dates. Only override if the dates are not
      the real ones; overriding is not allowed in confidence mode.
    rng: None (fresh randomness), an int seed, or a numpy Generator. Pass a
      seed for anything you report; the p-value is a Monte Carlo estimate.
    """
    config = fit.config
    if mode not in _VALID_MODES:
        raise ValueError(f"mode must be one of {sorted(_VALID_MODES)}, got {mode!r}.")
    if alternative not in _VALID_ALTERNATIVES:
        raise ValueError(
            f"alternative must be one of {sorted(_VALID_ALTERNATIVES)}, got {alternative!r}."
        )
    seed = int(rng) if isinstance(rng, (int, np.integer)) and not isinstance(rng, bool) else None
    rng = _resolve_rng(rng)
    nb = len(fit.bins)

    pivot, counts = _bin_and_pivot(treatment_df, fit.date_col, fit.value_col, config.bin_days,
                                   config.agg, with_counts=True)
    pivot, counts = pivot.reindex(columns=fit.bins), counts.reindex(columns=fit.bins).fillna(0)
    if fit.typical_counts is not None:
        pivot, n_thin = _mask_thin_bins(pivot, counts, fit.typical_counts, config.min_bin_coverage)
        if n_thin:
            warnings.warn(
                f"{n_thin} treatment bin(s) had fewer than {config.min_bin_coverage:.0%} of the "
                f"usual number of measurements and are treated as missing.",
                UserWarning, stacklevel=2,
            )
        n_dense = int((counts > 2 * fit.typical_counts).to_numpy().sum())
        if n_dense:
            warnings.warn(
                f"{n_dense} treatment bin(s) have more than twice the usual number of "
                f"measurements (e.g. continuous sensor data tested against a history of "
                f"grab samples). Their summaries vary less than the historical ones, so "
                f"the comparison is not like for like. Consider using only a subset of "
                f"the tested data that matches the historical sampling frequency.",
                UserWarning, stacklevel=2,
            )
    pivot = pivot.dropna(how="all")
    if len(pivot) == 0:
        raise ValueError(
            "The treatment data has no values in any of the fitted bins. Check "
            "the date/value columns and that it covers the same part of the year."
        )
    years = [int(y) for y in pivot.index]
    if mode == "prediction" and len(years) > 1:
        raise ValueError(
            f"mode='prediction' tests one calendar year, but the treatment data "
            f"spans {years}. A period crossing 1 January (e.g. a water year) must "
            f"be split into its calendar-year parts; to test the average of several "
            f"complete years, use mode='confidence'."
        )
    overlap = sorted(set(years) & set(fit.years_used))
    if overlap:
        raise ValueError(
            f"Treatment year(s) {overlap} are also in the historical record, so the "
            f"period would be compared partly with itself. Refit without them."
        )
    if treatment_year_index is not None:
        if mode == "confidence":
            raise ValueError(
                "treatment_year_index cannot be overridden in confidence mode; each "
                "year's own date is used."
            )
        from_dates = years[0] - fit.first_year
        if int(treatment_year_index) != from_dates:
            warnings.warn(
                f"treatment_year_index={treatment_year_index} but the data is dated "
                f"{years[0]} (index {from_dates}). Using the value you passed.",
                UserWarning, stacklevel=2,
            )
        year_indices = [float(treatment_year_index)]
    else:
        year_indices = [float(y - fit.first_year) for y in years]
    x0 = float(np.mean(year_indices))
    if fit.trend["applied"]:
        span = (fit.year_offsets[0], fit.year_offsets[-1])
        reach = max(max(x - span[1], span[0] - x) for x in year_indices)
        if reach > 5:
            warnings.warn(
                f"The historical trend is projected {reach:.0f} years beyond the historical "
                f"record. A trend estimated from the record becomes less reliable the further "
                f"it is projected: in simulation, when the removed trend was not real, the "
                f"false-alarm rate at a nominal 5% was about 8% one year beyond the record, "
                f"11-16% five years beyond and 16-29% ten years beyond. Check the result with "
                f"detrend_mode='none' (sensitivity_grid) and report both.",
                UserWarning, stacklevel=2,
            )

    raw = pivot.to_numpy(dtype="float64")
    use_log = fit.trend["log_space"]
    if use_log:
        _check_positive_for_log(raw, "the treatment data")
    work = np.log(raw) if use_log else raw.copy()
    slope = fit.trend["slope_per_year"] if fit.trend["applied"] else 0.0
    work = work - slope * np.asarray(year_indices)[:, None]
    missing = np.isnan(work).any(axis=0)
    curve = np.where(missing, np.nan, work.mean(axis=0) if len(work) else np.nan)

    n_missing = int(missing.sum())
    if n_missing == nb:
        raise ValueError("The treatment period has no complete bin to test.")
    if n_missing / nb > config.max_missing_frac:
        raise ValueError(
            f"The treatment period is missing {n_missing} of {nb} bins "
            f"({n_missing / nb:.0%}), above max_missing_frac={config.max_missing_frac:.0%}, "
            f"the same completeness bar historical years must clear. Supply more "
            f"complete data or use coarser bins."
        )
    if n_missing:
        warnings.warn(
            f"No treatment data in bins {fit.bins[missing].tolist()}"
            + (" (in at least one of the averaged years)" if mode == "confidence" else "")
            + f". They are left out of the test and of the null distribution; the "
              f"result covers the other {nb - n_missing} bins.",
            UserWarning, stacklevel=2,
        )

    # Scores and the null distribution, restricted to the same bins.
    z_obs = (curve - fit.median_curve) / fit.mad
    null_z = _null_scores(fit, config.n_bootstrap, rng, year_indices)
    null_z[:, missing] = np.nan
    obs_s = _directional(z_obs, alternative)
    null_s = _directional(null_z, alternative)
    B = config.n_bootstrap
    T = float(np.nanmax(obs_s))
    null_max = np.nanmax(null_s, axis=1)
    p_value = (np.sum(null_max >= T) + 1) / (B + 1)
    with np.errstate(invalid="ignore"):
        bin_p_adjusted = (np.sum(null_max[:, None] >= obs_s[None, :], axis=0) + 1) / (B + 1)
        bin_p_values = (np.sum(null_s >= obs_s[None, :], axis=0) + 1) / (B + 1)
    bin_p_adjusted = np.where(missing, np.nan, bin_p_adjusted)
    bin_p_values = np.where(missing, np.nan, bin_p_values)

    sustained = {}
    for k in (run_lengths or []):
        obs_w = _window_scores(z_obs, k, alternative, fit.bins)[0]
        if not np.any(np.isfinite(obs_w)):
            raise ValueError(
                f"No run of {k} consecutive bins is free of missing data in the "
                f"treatment period. Use a shorter run_length or more complete data."
            )
        start = int(np.argmax(obs_w))
        null_w = _window_scores(null_z, k, alternative, fit.bins).max(axis=1)
        sustained[k] = SustainedResult(
            run_length=k,
            p_value=float((np.sum(null_w >= obs_w[start]) + 1) / (B + 1)),
            test_statistic=float(obs_w[start]),
            window_start_bin=int(fit.bins[start]),
            window_bins=fit.bins[start:start + k],
        )

    # Report the tested curve back in original units, at the period's own
    # trend level (for one year this is exactly the binned data).
    shown = curve + slope * x0
    treatment_curve = np.exp(shown) if use_log else shown
    expected = fit.median_curve + slope * x0
    typical = fit.median_curve_raw
    if use_log:
        expected, typical = np.exp(expected), np.exp(typical)
    effect_size = float(np.nanmedian(treatment_curve - typical))
    effect_adj = float(np.nanmedian(treatment_curve - expected)) if fit.trend["applied"] else effect_size

    return TideResult(
        p_value=float(p_value), test_statistic=T, effect_size=effect_size,
        effect_size_trend_adjusted=effect_adj, bin_scores=z_obs,
        bin_p_adjusted=bin_p_adjusted, bin_p_values=bin_p_values,
        sustained=sustained, mode=mode, alternative=alternative,
        treatment_curve=treatment_curve, bins=fit.bins, n_reference=B,
        null_max_stats=null_max, treatment_year_index=x0,
        treatment_years=years, missing_bins=missing, seed=seed,
    )


test_treatment.__test__ = False   # stop pytest collecting it as a test


# ---------------------------------------------------------------------------
# Which bins, and the significance band
# ---------------------------------------------------------------------------

def critical_value(result: TideResult, alpha: float = 0.05) -> float:
    """The bin score a period must exceed to be significant at `alpha`: the
    period's p-value is <= alpha exactly when its largest (directional) bin
    score is above this number. Infinite when alpha is below the smallest
    p-value n_reference draws can produce."""
    if not 0 < alpha < 1:
        raise ValueError(f"alpha must be strictly between 0 and 1, got {alpha}.")
    B = result.n_reference
    k = int(np.floor(alpha * (B + 1) * (1 + 1e-9))) - 1
    if k < 0:
        return float("inf")
    return float(np.sort(result.null_max_stats)[::-1][k])


def bin_significant(result: TideResult, alpha: float = 0.05) -> np.ndarray:
    """Which bins were significantly unusual, with the chance of flagging
    ANY bin in a normal year held at alpha (single-step max-T adjustment).
    A bin is flagged exactly when its score is beyond critical_value(), so
    the period is significant if and only if at least one bin is flagged.
    Bins with no treatment data are never flagged. Boolean array aligned
    with result.bins."""
    p = np.asarray(result.bin_p_adjusted, dtype="float64")
    if 1.0 / (result.n_reference + 1) > alpha * (1 + 1e-9):
        warnings.warn(
            f"Nothing can be flagged at alpha={alpha}: {result.n_reference} bootstrap "
            f"draws cannot produce a p-value below {1 / (result.n_reference + 1):.4f}. "
            f"Raise n_bootstrap.",
            UserWarning, stacklevel=2,
        )
    return np.array([(not np.isnan(v)) and _le(v, alpha) for v in p], dtype=bool)


def significance_band(fit: TideFit, result: TideResult, alpha: float = 0.05) -> dict:
    """The band a normal period stays entirely inside with probability
    1 - alpha, in original units, on the tested period's trend level. The
    treatment curve leaves it in a bin exactly when that bin is flagged by
    bin_significant, and anywhere exactly when p_value <= alpha. One-sided
    tests have only one finite edge."""
    c = critical_value(result, alpha)
    slope = fit.trend["slope_per_year"] if fit.trend["applied"] else 0.0
    centre = fit.median_curve + slope * result.treatment_year_index
    half = c * fit.mad
    lower = centre - half if result.alternative != "greater" else np.full_like(centre, -np.inf)
    upper = centre + half if result.alternative != "less" else np.full_like(centre, np.inf)
    if fit.trend["log_space"]:
        centre, lower, upper = np.exp(centre), np.exp(lower), np.exp(upper)
    return {"bins": fit.bins, "lower": lower, "upper": upper, "centre": centre,
            "alpha": alpha, "critical_value": c}


# ---------------------------------------------------------------------------
# Departure extent within one period (descriptive)
# ---------------------------------------------------------------------------

@dataclass
class DepartureRecovery:
    run_length: int
    departure_start_bin: int        # first bin of the estimated departure
    departure_end_bin: int          # last bin still part of it
    recovered: bool                 # True if a measured, normal-looking bin followed
    recovered_at_bin: "int | None"  # first bin after the departure; None if it
                                    # lasted to the end of the period


def estimate_departure_recovery(result: TideResult, run_length: int,
                                extension_alpha: float = 0.1) -> DepartureRecovery:
    """DESCRIPTIVE, not a test. Starting from the core window found by the
    sustained test (result.sustained[run_length]), extend outwards one bin
    at a time while the neighbouring bin's UNadjusted p-value
    (bin_p_values) is below extension_alpha. Reports roughly when the
    departure started and ended, and whether a measured bin that looked
    normal followed it. Expect it to over-reach by about a bin."""
    if run_length not in result.sustained:
        raise ValueError(
            f"result.sustained has no entry for run_length={run_length}; call "
            f"test_treatment(..., run_lengths=[{run_length}]) first."
        )
    sw = result.sustained[run_length]
    n_bins = len(result.bins)
    core_start = {int(b): i for i, b in enumerate(result.bins)}[sw.window_start_bin]
    core_end = core_start + run_length - 1
    p = result.bin_p_values
    left = core_start
    while left - 1 >= 0 and p[left - 1] < extension_alpha:
        left -= 1
    right = core_end
    while right + 1 < n_bins and p[right + 1] < extension_alpha:
        right += 1
    recovered = right < n_bins - 1
    if recovered and np.isnan(p[right + 1]):
        warnings.warn(
            f"The bin after the departure (bin {int(result.bins[right + 1])}) has no "
            f"treatment data, so 'recovered' only means the departure stopped being "
            f"traceable, not that values were measured back at normal.",
            UserWarning, stacklevel=2,
        )
    return DepartureRecovery(
        run_length=run_length,
        departure_start_bin=int(result.bins[left]),
        departure_end_bin=int(result.bins[right]),
        recovered=recovered,
        recovered_at_bin=int(result.bins[right + 1]) if recovered else None,
    )


# ---------------------------------------------------------------------------
# Descriptive envelope (context for plots; the test's own band is
# significance_band)
# ---------------------------------------------------------------------------

def get_envelope(fit: TideFit, lower_q: float = 0.1, upper_q: float = 0.9,
                 method: str = "bootstrap", year_index: "float | None" = None,
                 rng=None) -> dict:
    """Pointwise percentile band of normal years, for context.

    method="bootstrap": percentiles of synthetic normal years from the same
      null as the test. method="empirical": quantiles of the N historical
      years themselves.
    year_index: trend level to draw it at (only matters when a trend was
      removed). None = one year after the last historical year.

    Each bin's band holds (upper_q - lower_q) of normal years in THAT bin;
    a normal year is expected to poke out of it in some bins. For the
    band that corresponds to the test, use significance_band.
    """
    if method not in ("bootstrap", "empirical"):
        raise ValueError(f"method must be 'bootstrap' or 'empirical', got {method!r}.")
    if not 0 <= lower_q < upper_q <= 1:
        raise ValueError("Need 0 <= lower_q < upper_q <= 1.")
    rng = _resolve_rng(rng)
    if year_index is None:
        year_index = float(fit.years_used[-1] - fit.first_year + 1)
    if method == "bootstrap":
        synth = fit.median_curve + _null_scores(fit, fit.config.n_bootstrap, rng,
                                                [year_index]) * fit.mad
        lower = np.percentile(synth, lower_q * 100, axis=0)
        upper = np.percentile(synth, upper_q * 100, axis=0)
    else:
        lower = np.nanquantile(fit.historical_matrix, lower_q, axis=0)
        upper = np.nanquantile(fit.historical_matrix, upper_q, axis=0)
    shift = fit.trend["slope_per_year"] * float(year_index) if fit.trend["applied"] else 0.0
    median_curve = fit.median_curve + shift
    year_curves = fit.historical_matrix + shift
    lower, upper = lower + shift, upper + shift
    if fit.trend["log_space"]:
        lower, upper = np.exp(lower), np.exp(upper)
        median_curve, year_curves = np.exp(median_curve), np.exp(year_curves)
    return {
        "bins": fit.bins, "lower": lower, "upper": upper,
        "median_curve": median_curve, "year_curves": year_curves,
        "years_used": fit.years_used, "n_years": len(fit.years_used),
        "method": method, "lower_q": lower_q, "upper_q": upper_q,
        "year_index": year_index, "trend_shift": float(shift),
    }
