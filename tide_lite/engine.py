"""
TIDE lite -- core engine (Phase 2 of the rebuild plan).

One job: given historical "normal" data and one treatment period, decide
whether the treatment period is unusual, accounting for seasonality and
WITHIN-year autocorrelation, without assuming a distribution.

Read the "Calibration and its limits" section below before trusting a
p-value near alpha: within-year dependence is handled, but whole-year
("this was a wet year") dependence is NOT, and that has a measured cost.

Method: historical years are aligned onto a common day-of-year axis and
binned (weekly/monthly). The typical seasonal shape is the bin-wise median
across historical years; each year's deviation from it is a residual
sequence. A circular block bootstrap (Politis & Romano 1992/1994; Kunsch
1989) resamples SUB-YEAR chunks of these residual sequences -- not whole
years -- and stitches them back together (wrapping at the year boundary,
since Dec is seasonally adjacent to Jan) into many synthetic "normal year"
residual sequences. Comparing the real treatment period's deviation
against this Monte Carlo null distribution gives the p-value.

Two corrections made after adversarial testing, not part of the original
plan -- both documented here because they matter for interpreting results,
not just as changelog entries:

1. Blocks are POSITION-MATCHED. A block used to fill calendar positions
   [s, s+L) is only ever drawn from historical years' OWN data at that
   exact position -- never from a different time of year. An earlier
   version pooled blocks globally regardless of origin, which could place
   naturally-volatile months' variance onto naturally-stable months'
   slots (or vice versa) in a synthetic draw, distorting the null
   distribution's shape in a position-dependent way. The starting phase is
   randomized per draw so coverage is still genuinely circular/moving,
   not a fixed tiling that always cuts seams at the same spots.

2. Resampling whole years (one block = one year) was tried first and
   rejected: with N historical years there are only N possible whole-year
   draws, so the smallest achievable p-value is 1/(N+1) -- with N=10 that
   floor is 0.09, *above* alpha=0.05, making significance mathematically
   unreachable regardless of effect size. Sub-year blocks fix this by
   combining many positions' worth of independently-drawn history into
   one synthetic year, which is why block_length_bins is deliberately
   shorter than a full year (see TideConfig).

CALIBRATION AND ITS LIMITS (measured by simulation, not derived on paper)

Feeding the pipeline fresh synthetic "normal" years and counting how often
it wrongly flags them, at a nominal alpha of 0.05, 900 trials per cell:

    historical years   within-year noise only   + whole-year level shifts
    ----------------   ----------------------   -------------------------
            10                  0.054                 0.112 - 0.130
            20                  0.042                 0.072 - 0.081
            40                  0.018                      --

So: when year-to-year variation is just accumulated within-year wiggle,
the test is well calibrated and gets conservative on long records. When
the data ALSO carries genuine whole-year level shifts -- a wet year
sitting high from January to December, which is the norm in hydrology --
the test runs at roughly twice its nominal false-positive rate with 10
historical years, and about 1.5x with 20.

The mechanism is design point 2 below, seen from the other side. A
synthetic "normal year" is stitched from ~4 independently drawn sub-year
chunks, so a whole-year offset present in every real year is averaged
away in the null: measured directly, the SD of the annual mean residual
across synthetic draws is about half that of the real historical years
once a year-level effect exists. The null is too narrow, so real years
look more extreme than they should.

This is a real property of the method, not a bug to patch out. Resampling
whole years instead is worse (design point 2). Adding a resampled
whole-year offset back onto each draw was prototyped and did NOT fix the
calibration -- with N historical years the offset distribution is itself
capped at the most extreme year ever observed, while a genuinely new year
exceeds all N of them about 2/(N+1) of the time. Short of a parametric
model of the year-level distribution -- which would give back the
distributional assumption this whole approach exists to avoid -- the
honest move is to measure the exposure and report it.

fit_historical therefore reports `between_year_var_frac` (the share of
residual variance carried by whole-year level shifts) and warns when it
is high on a short record. Low value -> the top table column applies.
High value -> the bottom one does; prefer more historical years and treat
p-values near alpha as indicative.

Other design choices, and why:
  - Median + MAD (median absolute deviation, 1.4826x-corrected to be
    comparable to a standard deviation) instead of mean + SD: robust to
    the skew that's typical of concentration data, and avoids a normality
    assumption.
  - No leave-one-out correction: each historical year's residual
    contributes to the median/MAD it's later measured against. With
    resampling pooling blocks across many years this influence is small
    and diluted; documented here rather than silently assumed.
  - detrend_mode="log_additive" requires strictly positive data and
    raises a clear error otherwise -- numpy's log() of a non-positive
    value silently produces NaN/a warning, and nanmedian silently drops
    NaNs, which without this check would quietly corrupt affected bins
    with no visible error under default warning settings.
"""

from __future__ import annotations

from dataclasses import dataclass
import warnings

import numpy as np
import pandas as pd
from scipy.stats import norm


# ---------------------------------------------------------------------------
# Day-of-year alignment (leap-year handling lives HERE, once, and nowhere
# else -- Phase 1 chose calendar day-of-year over a full Unix-timestamp
# rearchitecture, accepting ~5 days of drift over 20 years as documented).
# ---------------------------------------------------------------------------

def to_day_of_year(dates: pd.Series) -> pd.Series:
    """Map dates to a 1-365 day-of-year axis, same length in leap and
    non-leap years. Feb 29 becomes NaN (caller drops it -- at most one day
    lost every four years); every date from Mar 1 onward in a leap year
    shifts down by 1 so Dec 31 is always day 365.
    """
    dates = pd.to_datetime(dates)
    doy = dates.dt.dayofyear.astype("float64").copy()
    is_leap = dates.dt.is_leap_year
    is_feb29 = is_leap & (dates.dt.month == 2) & (dates.dt.day == 29)
    shift = is_leap & (doy >= 61) & (~is_feb29)
    doy = doy.where(~shift, doy - 1)
    doy = doy.where(~is_feb29, np.nan)
    return doy


def assign_bins(day_of_year: pd.Series, bin_days: int) -> pd.Series:
    """Group the 1-365 day-of-year axis into bins of `bin_days` days.
    bin_days=7 ~ weekly (53 bins); bin_days=30 ~ roughly monthly (13 bins).
    """
    return ((day_of_year - 1) // bin_days).astype("Int64") + 1


# ---------------------------------------------------------------------------
# Trend detection (Mann-Kendall) and magnitude (Sen's slope).
# ---------------------------------------------------------------------------

def mann_kendall_test(values) -> tuple[float, float]:
    """Two-sided Mann-Kendall trend test (Mann 1945; Kendall 1975).
    Returns (S statistic, p-value). NaNs dropped first."""
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
    """Sen's slope estimator (Sen 1968): median of all pairwise slopes.

    `times` is the x-axis the slope is expressed per unit of. It MUST be
    passed whenever the observations are not evenly spaced -- historical
    records routinely have gaps (a year dropped by the completeness
    filter, or simply never sampled), and using the positional index as
    a stand-in for the calendar silently rescales the slope. With
    years [2000..2004, 2010..2014] and a true drift of 0.5/year, the
    positional version returns ~0.87/year; passing the real calendar
    offsets returns ~0.5/year.
    """
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


def holm_bonferroni(p_values: np.ndarray, alpha: float = 0.05) -> np.ndarray:
    """Controls family-wise error at `alpha` across len(p_values) tests.
    Sort p ascending; compare p(1) to alpha/n, p(2) to alpha/(n-1), ...;
    stop rejecting at the first comparison that fails. Uniformly more
    powerful than naive Bonferroni (alpha/n applied to every test) while
    keeping the same FWER guarantee -- naive Bonferroni is the version
    that's easy to implement by accident.

    General-purpose: used across events in Sequential (controllers.py)
    and across bins within one period in bin_significant() below -- same
    correction, two different families of tests being corrected.
    """
    p_values = np.asarray(p_values, dtype="float64")
    if not 0 < alpha < 1:
        raise ValueError(f"alpha must be strictly between 0 and 1, got {alpha}.")
    if np.any(np.isnan(p_values)):
        raise ValueError(
            "holm_bonferroni received NaN p-values. Drop or mask them first -- "
            "a NaN is missing evidence, and silently sorting it into the "
            "sequence would change every other test's threshold."
        )
    n = len(p_values)
    order = np.argsort(p_values)
    sorted_p = p_values[order]
    reject_sorted = np.zeros(n, dtype=bool)
    for i, p in enumerate(sorted_p):
        threshold = alpha / (n - i)
        if p <= threshold:
            reject_sorted[i] = True
        else:
            break  # Holm's procedure stops at the first non-rejection
    reject = np.zeros(n, dtype=bool)
    reject[order] = reject_sorted
    return reject


# ---------------------------------------------------------------------------
# rng convenience: accept None, an int seed, or an existing Generator.
# ---------------------------------------------------------------------------

def _resolve_rng(rng) -> np.random.Generator:
    if isinstance(rng, np.random.Generator):
        return rng
    return np.random.default_rng(rng)  # rng=None -> fresh entropy; rng=42 -> seeded


# ---------------------------------------------------------------------------
# Config and result containers
# ---------------------------------------------------------------------------

_VALID_DETREND_MODES = {"additive", "log_additive", "none"}
_VALID_AGG = {"median", "mean"}
_VALID_MODES = {"prediction", "confidence"}


@dataclass
class TideConfig:
    """The Phase 1 decisions, in one place, with defaults documented."""
    bin_days: int = 30              # ~monthly bins; use 7 for weekly
    agg: str = "median"             # "median" or "mean" within each bin
    detrend_mode: str = "additive"  # "additive", "log_additive", or "none" --
                                     # "log_additive" also serves as this
                                     # tool's variance-stabilizing normalization
                                     # for skewed data; see USER_GUIDE.md's note
                                     # on normalization for why it's bundled here
                                     # rather than a separate setting, and what
                                     # this project deliberately does NOT do
                                     # (baseline/reference-level rescaling)
    mk_alpha: float = 0.10          # trend-test threshold; lenient on purpose
    max_missing_frac: float = 0.5   # drop historical years missing > this
    n_bootstrap: int = 2000         # Monte Carlo draws for the null distribution
    block_length_bins: int | None = None  # None = auto (n_bins // 4, min 2)
    pool_window_radius: int = 1     # pool each position with its +/- this-many
                                     # immediate neighbors (see _extract_block_pool)

    def __post_init__(self):
        """Reject misspelled options loudly. Every one of these used to fall
        through silently to a DIFFERENT analysis: detrend_mode="log" ran an
        additive detrend with no log transform at all, and agg="avg" only
        surfaced later as a pandas AttributeError from deep inside a groupby.
        """
        if self.detrend_mode not in _VALID_DETREND_MODES:
            raise ValueError(
                f"detrend_mode must be one of {sorted(_VALID_DETREND_MODES)}, "
                f"got {self.detrend_mode!r}."
            )
        if self.agg not in _VALID_AGG:
            raise ValueError(
                f"agg must be one of {sorted(_VALID_AGG)}, got {self.agg!r}."
            )
        if self.bin_days < 1 or self.bin_days > 365:
            raise ValueError(f"bin_days must be between 1 and 365, got {self.bin_days}.")
        if not 0 < self.mk_alpha < 1:
            raise ValueError(f"mk_alpha must be strictly between 0 and 1, got {self.mk_alpha}.")
        if not 0 <= self.max_missing_frac < 1:
            raise ValueError(
                f"max_missing_frac must be in [0, 1), got {self.max_missing_frac}."
            )
        if self.n_bootstrap < 1:
            raise ValueError(f"n_bootstrap must be >= 1, got {self.n_bootstrap}.")
        if self.block_length_bins is not None and self.block_length_bins < 1:
            raise ValueError(
                f"block_length_bins must be >= 1 or None (auto), got {self.block_length_bins}."
            )
        if self.pool_window_radius < 0:
            raise ValueError(
                f"pool_window_radius must be >= 0, got {self.pool_window_radius}."
            )

    @property
    def min_attainable_p(self) -> float:
        """The smallest p-value this many bootstrap draws can ever produce:
        (0 + 1) / (n_bootstrap + 1). Any threshold below this is
        unreachable -- the test can never flag anything, at any effect
        size. Guards elsewhere in the package compare against this.
        """
        return 1.0 / (self.n_bootstrap + 1)


@dataclass
class TideFit:
    config: TideConfig
    bins: np.ndarray
    median_curve: np.ndarray
    mad: np.ndarray
    residual_matrix: np.ndarray     # (n_years, n_bins) detrended deviations
    historical_matrix: np.ndarray   # (n_years, n_bins) detrended values
    block_pool: dict                # {start_bin: (n_years_at_that_bin, block_length_bins)}
    block_length_bins: int          # resolved value actually used
    min_blocks_per_position: int    # smallest pool size across all start positions --
                                     # low values mean some calendar positions are
                                     # resampled from very few real historical years
    years_used: list
    years_dropped: list
    trend: dict
    date_col: str
    value_col: str
    median_curve_raw: np.ndarray = None   # bin-wise median BEFORE detrending, in the
                                           # same space as median_curve (log if log mode).
                                           # effect_size is measured against this so the
                                           # "raw" number means "vs. the historical period's
                                           # own typical level", not "vs. the level at the
                                           # START of the record" -- see test_treatment.
    between_year_var_frac: float = 0.0     # share of residual variance carried by
                                           # WHOLE-YEAR level shifts (wet year / dry year)
                                           # rather than within-year wiggle. High values
                                           # mean the block bootstrap's null is too narrow
                                           # and p-values are anti-conservative -- see the
                                           # module docstring's calibration section.

    @property
    def first_year(self) -> int:
        """Calendar year the trend axis is measured from (years_used[0])."""
        return int(self.years_used[0])

    def year_index_for(self, calendar_year: int) -> int:
        """Convert a calendar year into the `treatment_year_index` that
        test_treatment / Sequential / Cumulative / MonitoringSeries expect:
        the offset in CALENDAR years from the first historical year.

            fit.year_index_for(2024)   # -> 24 if history starts in 2000

        Use this instead of len(fit.years_used) -- the two agree only when
        the historical record has no gaps.
        """
        return int(calendar_year) - self.first_year


@dataclass
class SustainedResult:
    """One run-length's sustained-departure test: was there a window of
    `run_length` (or more, since a longer real departure will still
    contain a k-long sustained sub-window) CONSECUTIVE bins that,
    together, were unusual -- a different, complementary question to any
    single bin's own significance (see bin_significant()). Windows here
    do NOT wrap across the bin axis (bin 13 and bin 1 are not treated as
    adjacent) -- a treatment period is a bounded window in time, not a
    repeating cycle to resample from, unlike the historical block pool.
    """
    run_length: int
    p_value: float
    test_statistic: float           # the observed window's mean standardized deviation
    window_start_bin: int           # first bin (by label) of the window that achieved it
    window_bins: np.ndarray         # all bin labels in that window


@dataclass
class TideResult:
    p_value: float
    test_statistic: float
    effect_size: float              # signed, original units, NOT trend-adjusted --
                                     # "how different did the observed data look"
    effect_size_trend_adjusted: float  # same, but with the extrapolated historical
                                     # trend subtracted first -- matches what the
                                     # p-value actually tested. Equal to effect_size
                                     # when no trend was applied.
    bin_p_values: np.ndarray        # per-bin, UNCORRECTED pointwise p-values --
                                     # see bin_significant() for the corrected version;
                                     # do not compare these to alpha directly
    sustained: dict                 # {run_length: SustainedResult}, empty unless
                                     # test_treatment was called with run_lengths=[...]
    mode: str                       # "prediction" or "confidence"
    treatment_curve: np.ndarray
    bins: np.ndarray
    n_reference: int                # Monte Carlo draws used for the null
    treatment_year_index: int = 0   # resolved calendar-year offset from the first
                                     # historical year that this period was tested at.
                                     # Plotting needs it to put the envelope on the
                                     # same trend basis as the treatment curve.
    missing_bins: np.ndarray = None  # bool per bin: no treatment data there. Those
                                     # bins have bin_p_values = NaN and take no part
                                     # in the test -- they are not evidence either way.


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _bin_and_pivot(df, date_col, value_col, bin_days, agg):
    d = df[[date_col, value_col]].copy()
    d["doy"] = to_day_of_year(d[date_col])
    d = d.dropna(subset=["doy"])
    d["bin"] = assign_bins(d["doy"], bin_days)
    d["block_year"] = pd.to_datetime(d[date_col]).dt.year
    return d.pivot_table(index="block_year", columns="bin", values=value_col, aggfunc=agg)


def _check_positive_for_log(values, context: str):
    arr = np.asarray(values, dtype="float64")
    n_bad = int(np.nansum(arr <= 0))
    if n_bad > 0:
        raise ValueError(
            f"detrend_mode='log_additive' requires strictly positive values, "
            f"but {context} has {n_bad} value(s) <= 0. Use detrend_mode='additive' "
            f"instead, or fix/offset the data before fitting -- np.log of a "
            f"non-positive value fails silently (NaN + a warning you likely "
            f"won't see, not an error), so this is checked explicitly."
        )


def _extract_block_pool(residual_matrix: np.ndarray, block_length_bins: int,
                         window_radius: int = 1) -> dict:
    """{start_bin: array of shape (n_candidates, block_length_bins)} for
    every start_bin 0..n_bins-1, wrapped circularly at the year boundary.

    Each position pools from ITSELF plus its +/- window_radius immediate
    neighboring positions, across every historical year -- not the exact
    position alone (with only N historical years, that leaves just N
    candidates per position, which in testing produced enough fit-to-fit
    calibration variance to matter -- see tests/test_engine.py), and not
    the whole pool globally regardless of position (which was the
    original version's bug: it could place a naturally volatile month's
    variance onto a naturally stable month's slot). A radius of 1 keeps
    the pool restricted to seasonally adjacent positions only.
    """
    n_years, n_bins = residual_matrix.shape
    L = block_length_bins
    pool = {}
    for start in range(n_bins):
        year_blocks = []
        for offset in range(-window_radius, window_radius + 1):
            s = (start + offset) % n_bins
            for y in range(n_years):
                row = residual_matrix[y]
                extended = np.concatenate([row, row[:L]])
                block = extended[s:s + L]
                if not np.any(np.isnan(block)):
                    year_blocks.append(block)
        if year_blocks:
            pool[start] = np.array(year_blocks)
    if not pool:
        raise ValueError(
            "No complete blocks survived pooling at any calendar position -- "
            "historical data has too many gaps for the chosen "
            "block_length_bins. Try a shorter block_length_bins or check for "
            "large missing stretches."
        )
    return pool


def _nearest_pool(block_pool: dict, pos: int, n_bins: int) -> np.ndarray:
    """Fallback for the rare case where NO historical year has a valid
    block at this exact position (e.g. every year has a gap there):
    search outward, position by position, for the nearest one that does.
    """
    for offset in range(1, n_bins):
        for alt in ((pos - offset) % n_bins, (pos + offset) % n_bins):
            if alt in block_pool:
                return block_pool[alt]
    raise ValueError("Block pool is empty at every position.")  # pool non-empty overall, so unreachable


def _draw_synthetic_residual(block_pool: dict, n_bins: int, block_length_bins: int,
                              rng: np.random.Generator) -> np.ndarray:
    """One circular block bootstrap draw. The starting phase is randomized
    so coverage is genuinely moving/circular rather than a fixed tiling
    that always cuts seams at the same spots; each block filling positions
    [s, s+L) is drawn only from historical years' own data at position s
    (position-matched -- see module docstring).
    """
    phi = int(rng.integers(0, n_bins))
    aligned = np.empty(n_bins)
    pos = phi
    filled = 0
    while filled < n_bins:
        candidates = block_pool.get(pos)
        if candidates is None:
            candidates = _nearest_pool(block_pool, pos, n_bins)
        block = candidates[rng.integers(0, len(candidates))]
        for v in block:
            aligned[pos % n_bins] = v
            pos += 1
            filled += 1
            if filled >= n_bins:
                break
    return aligned


def _max_standardized_deviation(curve, median_curve, mad) -> float:
    return float(np.nanmax(np.abs(curve - median_curve) / mad))


def _max_window_mean(bin_devs: np.ndarray, run_length: int) -> tuple[float, int]:
    """Max, over all LINEAR (non-wrapping) windows of run_length
    consecutive bins, of the mean deviation within that window. Returns
    (max_mean, start_index). Non-wrapping deliberately: a treatment
    period is a bounded stretch of time, not a repeating cycle -- unlike
    the historical block pool, there's no "next January" to wrap into
    within one period.
    """
    n = len(bin_devs)
    if run_length < 1:
        raise ValueError(f"run_length must be >= 1, got {run_length}.")
    if run_length > n:
        raise ValueError(f"run_length ({run_length}) exceeds the number of bins ({n}).")
    # Windows containing a missing bin are EXCLUDED rather than averaged
    # over what is present. A partial window is not comparable to the null,
    # which always averages run_length complete bins; and np.mean of a
    # window holding a NaN returns NaN, which np.argmax then selects as the
    # maximum -- previously handing the sustained test a NaN statistic that
    # no null draw could exceed, i.e. the smallest possible p-value for a
    # window that was really just missing data.
    window_means = np.array([
        bin_devs[s:s + run_length].mean() for s in range(n - run_length + 1)
    ])
    valid = ~np.isnan(window_means)
    if not np.any(valid):
        raise ValueError(
            f"No window of {run_length} consecutive bins is free of missing "
            f"data in the treatment period, so the sustained-departure test "
            f"has nothing comparable to the null to measure. Use a shorter "
            f"run_length, coarser bins (larger bin_days), or more complete "
            f"treatment data."
        )
    window_means = np.where(valid, window_means, -np.inf)
    best = int(np.argmax(window_means))
    return float(window_means[best]), best


# ---------------------------------------------------------------------------
# fit_historical
# ---------------------------------------------------------------------------

def fit_historical(df: pd.DataFrame, date_col: str, value_col: str,
                    config: "TideConfig | None" = None) -> TideFit:
    """Align, bin, filter, detrend, and summarize the historical
    ("business as usual") years, and build the block pool used to test
    against later.
    """
    config = config or TideConfig()
    pivot = _bin_and_pivot(df, date_col, value_col, config.bin_days, config.agg)

    missing_frac = pivot.isna().mean(axis=1)
    keep_mask = missing_frac <= config.max_missing_frac
    years_dropped = [
        (int(y), f"{missing_frac[y]:.0%} of bins missing")
        for y in pivot.index[~keep_mask]
    ]
    pivot = pivot.loc[keep_mask]
    if len(pivot) < 3:
        raise ValueError(
            f"Only {len(pivot)} historical years survived the completeness "
            f"filter (need >= 3). Lower max_missing_frac or check the data."
        )

    years_used = list(pivot.index.astype(int))
    matrix = pivot.to_numpy(dtype="float64")
    bins = pivot.columns.to_numpy()
    n_bins = len(bins)

    use_log = config.detrend_mode == "log_additive"
    if use_log:
        _check_positive_for_log(matrix, "the historical data")
    work = np.log(matrix) if use_log else matrix.copy()

    # The trend axis is CALENDAR years elapsed since the first historical
    # year, not the position of each year in the list. Those coincide only
    # for a gapless record; with any year missing (never sampled, or dropped
    # by the completeness filter above) the positional version silently
    # rescales the slope -- see sens_slope's docstring for a worked case.
    year_offsets = np.asarray(years_used, dtype="float64") - float(years_used[0])
    median_curve_raw = np.nanmedian(work, axis=0)

    yearly_summary = np.nanmedian(work, axis=1)
    mk_s, mk_p = mann_kendall_test(yearly_summary)
    slope = 0.0
    trend_applied = False
    if config.detrend_mode != "none" and mk_p < config.mk_alpha:
        slope = sens_slope(yearly_summary, times=year_offsets)
        work = work - (slope * year_offsets)[:, None]
        trend_applied = True

    median_curve = np.nanmedian(work, axis=0)
    residual_matrix = work - median_curve
    # 1.4826x is the standard consistency correction making MAD comparable
    # to a standard deviation for roughly-normal data (Rousseeuw & Croux
    # 1993). A uniform per-bin rescale, so it doesn't change any p-value
    # (cancels out of the null-vs-observed rank comparison) -- included so
    # mad and the test statistic are interpretable on a familiar scale.
    mad = 1.4826 * np.nanmedian(np.abs(residual_matrix), axis=0)
    # Last-resort fallback (a bin where >= half of years tie exactly, so
    # MAD=0) is scaled to the data's own level rather than a fixed 1.0,
    # which would be a unit mismatch for arbitrarily-scaled variables.
    fallback = np.nanmedian(mad[mad > 0]) if np.any(mad > 0) else \
        0.01 * max(np.nanmedian(np.abs(median_curve)), 1e-9)
    mad = np.where(mad == 0, fallback, mad)

    # How much of the leftover variation is a WHOLE-YEAR level shift (a wet
    # year sitting high all year) versus within-year wiggle. The bootstrap
    # stitches each synthetic year from several independently drawn chunks,
    # which averages whole-year shifts away -- so the larger this fraction,
    # the narrower the null is relative to what a real new year does, and the
    # more anti-conservative the p-value. Measured, not assumed: see the
    # calibration table in the module docstring.
    _year_levels = np.nanmean(residual_matrix, axis=1)
    _total_var = float(np.nanvar(residual_matrix))
    between_year_var_frac = (
        float(np.clip(np.nanvar(_year_levels) / _total_var, 0.0, 1.0))
        if _total_var > 0 else 0.0
    )
    if between_year_var_frac >= 0.5 and len(years_used) < 20:
        warnings.warn(
            f"{between_year_var_frac:.0%} of this record's residual variation is "
            f"whole-year level shifts, and only {len(years_used)} historical years "
            f"are available. In simulation this regime runs at roughly 2x the "
            f"nominal false-positive rate (11-13% at a nominal 5% with 10 years). "
            f"Treat p-values near alpha as indicative, prefer more historical "
            f"years, and read the 'Calibration' section of USER_GUIDE.md.",
            UserWarning, stacklevel=2,
        )

    block_length_bins = config.block_length_bins or max(2, n_bins // 4)
    block_length_bins = min(block_length_bins, n_bins)
    block_pool = _extract_block_pool(residual_matrix, block_length_bins,
                                      window_radius=config.pool_window_radius)
    min_blocks_per_position = min(len(v) for v in block_pool.values())

    return TideFit(
        config=config, bins=bins, median_curve=median_curve, mad=mad,
        residual_matrix=residual_matrix, historical_matrix=work,
        block_pool=block_pool, block_length_bins=block_length_bins,
        min_blocks_per_position=min_blocks_per_position,
        years_used=years_used, years_dropped=years_dropped,
        trend={"mk_p": mk_p, "slope_per_year": slope, "applied": trend_applied,
               "log_space": use_log, "first_year": int(years_used[0])},
        date_col=date_col, value_col=value_col,
        median_curve_raw=median_curve_raw,
        between_year_var_frac=between_year_var_frac,
    )


# ---------------------------------------------------------------------------
# test_treatment
# ---------------------------------------------------------------------------

def test_treatment(fit: TideFit, treatment_df: pd.DataFrame,
                    treatment_year_index: "int | None" = None,
                    mode: str = "prediction",
                    run_lengths: "list[int] | None" = None,
                    rng=None) -> TideResult:
    """Test one treatment period against the fitted historical envelope.

    treatment_year_index: position (0 = first historical year) used to
    extrapolate the historical trend forward for a fair comparison. If
    None, assumes immediately after the historical years used in the fit;
    Sequential/Cumulative require it explicitly whenever a trend
    correction is active, since their events aren't assumed adjacent.
    mode: "prediction" (default) tests exactly one treatment year -- raises
    if treatment_df spans more than one (silently testing only the first
    year of accidentally-multi-year data is worse than erroring). Use
    mode="confidence" to average several years/periods in treatment_df
    into one curve and test that average.
    run_lengths: optional list of consecutive-bin window lengths (e.g.
    [2, 3]) to ALSO test for a sustained departure lasting at least that
    many bins in a row -- a different, complementary question to the main
    (single most extreme bin) test: more power for a persistent, moderate
    shift that no single bin is extreme enough to flag alone; less power
    for a sharp, isolated one-bin spike, since averaging within the
    window dilutes it. Reuses the same bootstrap draws already being
    generated for the main test, so this costs almost nothing extra.
    Results land in result.sustained[run_length]. Default None: skip it.
    rng: None (fresh entropy), an int (seed), or a numpy Generator.
    """
    config = fit.config
    rng = _resolve_rng(rng)
    n_bins = len(fit.bins)

    if mode not in _VALID_MODES:
        raise ValueError(
            f"mode must be one of {sorted(_VALID_MODES)}, got {mode!r}. "
            f"(A misspelling used to fall through to the 'confidence' branch "
            f"and silently run a different analysis.)"
        )

    pivot = _bin_and_pivot(treatment_df, fit.date_col, fit.value_col,
                            config.bin_days, config.agg)
    if len(pivot) == 0:
        raise ValueError("No usable rows in treatment_df after alignment.")
    if mode == "prediction" and len(pivot) > 1:
        raise ValueError(
            f"mode='prediction' expects exactly one year of treatment data, "
            f"but treatment_df spans {len(pivot)} years ({list(pivot.index)}). "
            f"Restrict treatment_df to one year, or use mode='confidence' if "
            f"you intend to test the average across multiple years."
        )
    n_treatment_years = len(pivot) if mode == "confidence" else 1
    row = pivot.mean(axis=0) if mode == "confidence" else pivot.iloc[0]
    treatment_curve = row.reindex(fit.bins).to_numpy(dtype="float64")

    # Bins the historical fit knows about but the treatment period has no
    # data for. These carry NO information and must not be scored: an
    # earlier version compared NaN against every null draw, and since
    # `null >= nan` is False everywhere the bin scored 0 exceedances and
    # came out at the SMALLEST p-value the bootstrap can produce -- so a
    # treatment year truncated in August was reported as maximally
    # significant in exactly the months it had no data for.
    missing_bins = np.isnan(treatment_curve)
    n_missing = int(missing_bins.sum())
    if n_missing == n_bins:
        raise ValueError(
            "The treatment period has no data in any of the historical fit's "
            f"{n_bins} bins. Check that treatment_df covers the same part of "
            f"the year as the historical data and uses the same units/columns."
        )
    if n_missing / n_bins > config.max_missing_frac:
        raise ValueError(
            f"The treatment period is missing {n_missing} of {n_bins} bins "
            f"({n_missing / n_bins:.0%}), above max_missing_frac="
            f"{config.max_missing_frac:.0%} -- the same completeness bar "
            f"historical years have to clear. Testing a period this sparse "
            f"compares a few months against a whole-year null. Supply more "
            f"complete data, use coarser bins (larger bin_days), or raise "
            f"max_missing_frac deliberately if you accept the tradeoff."
        )
    if n_missing:
        warnings.warn(
            f"The treatment period has no data in {n_missing} of {n_bins} bins "
            f"({fit.bins[missing_bins].tolist()}). Those bins are excluded from "
            f"the test and their bin_p_values are NaN; the p-value is based on "
            f"the {n_bins - n_missing} bins that do have data.",
            UserWarning, stacklevel=2,
        )

    use_log = fit.trend["log_space"]
    if use_log:
        _check_positive_for_log(treatment_curve, "treatment_df")
    work_curve = np.log(treatment_curve) if use_log else treatment_curve

    # Resolved regardless of whether a trend was applied, so the result can
    # report the basis it was tested on and plotting can match it.
    if treatment_year_index is None:
        # One calendar year after the LAST historical year. len(years_used)
        # only equals that for a gapless record.
        treatment_year_index = int(fit.years_used[-1]) - int(fit.years_used[0]) + 1
    if fit.trend["applied"]:
        work_curve = work_curve - fit.trend["slope_per_year"] * treatment_year_index

    stat_obs = _max_standardized_deviation(work_curve, fit.median_curve, fit.mad)
    obs_bin_devs = np.abs(work_curve - fit.median_curve) / fit.mad  # per-bin, pre-max

    null_stats = np.empty(config.n_bootstrap)
    null_bin_devs = np.empty((config.n_bootstrap, n_bins))  # kept per-bin, not just the max,
    # so a per-bin (monthly) breakdown is available at no extra bootstrap cost -- see
    # bin_p_values on TideResult and bin_significant() below.
    for b in range(config.n_bootstrap):
        draws = [
            _draw_synthetic_residual(fit.block_pool, n_bins, fit.block_length_bins, rng)
            for _ in range(n_treatment_years)
        ]
        synthetic_residual = np.mean(draws, axis=0)
        bin_devs = np.abs(synthetic_residual) / fit.mad  # residual already centered on 0
        null_bin_devs[b] = bin_devs
        null_stats[b] = float(np.nanmax(bin_devs))
    n_ref = config.n_bootstrap
    p_value = (np.sum(null_stats >= stat_obs) + 1) / (n_ref + 1)
    # Pointwise, UNCORRECTED per-bin p-values -- same exact-permutation
    # formula as the global one, just without taking the max first. Do
    # NOT compare these to alpha directly across many bins (that's the
    # exact multiple-comparisons inflation the max-statistic exists to
    # avoid) -- use bin_significant() below, which corrects properly.
    bin_p_values = (np.sum(null_bin_devs >= obs_bin_devs[None, :], axis=0) + 1) / (n_ref + 1)
    # A bin with no treatment data gets NaN, not a p-value. See the
    # missing_bins comment above for what the old behaviour did instead.
    bin_p_values = np.where(missing_bins, np.nan, bin_p_values)

    # Sustained-departure test(s): a different question from bin_p_values
    # above -- "was there a run of run_length-or-more consecutive bins
    # that were TOGETHER unusual," not "was any single bin unusual."
    # Reuses obs_bin_devs / null_bin_devs already computed above; no
    # extra bootstrap draws needed.
    sustained = {}
    for k in (run_lengths or []):
        obs_window_stat, obs_start = _max_window_mean(obs_bin_devs, k)
        null_window_stats = np.array([
            _max_window_mean(null_bin_devs[b], k)[0] for b in range(n_ref)
        ])
        p_sustained = (np.sum(null_window_stats >= obs_window_stat) + 1) / (n_ref + 1)
        sustained[k] = SustainedResult(
            run_length=k, p_value=float(p_sustained), test_statistic=obs_window_stat,
            window_start_bin=int(fit.bins[obs_start]),
            window_bins=fit.bins[obs_start:obs_start + k],
        )

    # Raw effect size: how different the observed data looks, no trend
    # adjustment. Trend-adjusted: matches what the p-value actually
    # tested. The two differ only when a trend correction was applied --
    # report both rather than picking one, since silently only reporting
    # the raw figure can look inconsistent next to a trend-adjusted p-value
    # ("why is such a small effect so significant") and vice versa.
    # effect_size is measured against median_curve_raw -- the historical
    # bin-wise median BEFORE detrending, i.e. the record's own typical level.
    # It used to be measured against median_curve, which is the detrended
    # curve and therefore sits at the level of the FIRST historical year: on
    # a record with a real drift that made "raw" effect size grow with the
    # length of the record rather than describe the treatment period. On a
    # 15-year record drifting 0.6/year, a perfectly ordinary continuation of
    # the trend reported effect_size = +9.06 against a p-value of 0.40.
    # Fall back for a TideFit built before median_curve_raw existed (or
    # constructed by hand): with no trend applied the two are identical.
    baseline_raw = fit.median_curve if fit.median_curve_raw is None else fit.median_curve_raw
    if use_log:
        effect_size = float(np.nanmedian(treatment_curve - np.exp(baseline_raw)))
        effect_size_trend_adjusted = float(
            np.nanmedian(np.exp(work_curve) - np.exp(fit.median_curve))
        ) if fit.trend["applied"] else effect_size
    else:
        effect_size = float(np.nanmedian(treatment_curve - baseline_raw))
        effect_size_trend_adjusted = float(
            np.nanmedian(work_curve - fit.median_curve)
        ) if fit.trend["applied"] else effect_size

    return TideResult(
        p_value=float(p_value), test_statistic=stat_obs, effect_size=effect_size,
        effect_size_trend_adjusted=effect_size_trend_adjusted, bin_p_values=bin_p_values,
        sustained=sustained,
        mode=mode, treatment_curve=treatment_curve, bins=fit.bins, n_reference=n_ref,
        treatment_year_index=int(treatment_year_index), missing_bins=missing_bins,
    )


test_treatment.__test__ = False


def bin_significant(result: TideResult, alpha: float = 0.05) -> np.ndarray:
    """Which specific bins (months, if bin_days=30) were significantly
    unusual, Holm-Bonferroni corrected across the bins within this one
    period -- the answer to "the p-value is for the whole year, can we
    tell month by month?" without reintroducing the multiple-comparisons
    inflation that a single global test statistic exists to avoid in the
    first place (the same reasoning as testing 365 days pointwise against
    an envelope, just at whatever resolution bin_days was set to).

    This is a SEPARATE correction from result.p_value: the global p-value
    already controls for "did the year as a whole look unusual" using the
    max-deviation statistic; this additionally asks "which specific bins
    were driving that," corrected across just those bins. The two are
    related but different statistics (one about the single most extreme
    bin relative to a null built the same way across all bins; the other
    about each bin's own value relative to its own null, then corrected
    across bins) and can disagree -- check both rather than assuming a
    small global p-value means some specific bin will always survive this
    correction, or vice versa.

    Bins with no treatment data (NaN p-value) are never flagged and take no
    part in the correction -- they are missing evidence, not evidence of
    normality, and including them would also make every other bin's
    threshold stricter for no reason.

    Returns a boolean array aligned with result.bins / result.bin_p_values.
    """
    p = np.asarray(result.bin_p_values, dtype="float64")
    testable = ~np.isnan(p)
    n_testable = int(testable.sum())
    out = np.zeros(len(p), dtype=bool)
    if n_testable == 0:
        return out
    min_p = 1.0 / (result.n_reference + 1)
    strictest = alpha / n_testable
    if min_p > strictest:
        warnings.warn(
            f"No bin can be flagged at alpha={alpha}: the strictest "
            f"Holm threshold across {n_testable} bins is {strictest:.5f}, but "
            f"{result.n_reference} bootstrap draws can never produce a p-value "
            f"below {min_p:.5f}. This returns all-False for arithmetic reasons, "
            f"not because the data looked normal. Raise n_bootstrap to at least "
            f"{int(np.ceil(n_testable / alpha))} to make this test resolvable.",
            UserWarning, stacklevel=2,
        )
    out[testable] = holm_bonferroni(p[testable], alpha=alpha)
    return out


@dataclass
class DepartureRecovery:
    run_length: int
    departure_start_bin: int        # first bin (label) of the estimated departure
    departure_end_bin: int          # last bin (label) still part of it
    recovered: bool                 # True if the departure ended before the period did
    recovered_at_bin: "int | None"  # first bin back to normal; None if the departure
                                     # persisted through the last bin of the period


def estimate_departure_recovery(result: TideResult, run_length: int,
                                 extension_alpha: float = 0.1) -> DepartureRecovery:
    """Within ONE tested period: given a sustained departure was found
    (result.sustained[run_length], from test_treatment's run_lengths=
    [...]), estimate how far it actually extends and whether it recovered
    before the period ended.

    This is deliberately a DESCRIPTIVE boundary trace, not a hypothesis
    test in its own right -- result.sustained[run_length] already
    establishes formal significance for the core window; this walks
    outward from that core bin by bin, in each direction, extending
    while the individual bin still looks unusual (its own bin_p_values
    entry below extension_alpha -- a looser, UNCORRECTED threshold on
    purpose, since the core's significance was already established and
    this step is tracing an extent, not re-testing from scratch; using
    bin_significant()'s full correction at every candidate boundary bin
    would be a different, much more conservative question -- "roughly
    how long did this last" rather than "prove each additional bin was
    independently significant").

    Raises if result.sustained has no entry for run_length -- call
    test_treatment(..., run_lengths=[run_length, ...]) first.
    """
    if run_length not in result.sustained:
        raise ValueError(
            f"result.sustained has no entry for run_length={run_length}; call "
            f"test_treatment(..., run_lengths=[{run_length}]) (or including it "
            f"alongside other lengths) before calling this."
        )
    sw = result.sustained[run_length]
    n_bins = len(result.bins)
    bin_to_idx = {int(b): i for i, b in enumerate(result.bins)}
    core_start = bin_to_idx[sw.window_start_bin]
    core_end = core_start + run_length - 1  # inclusive, linear (no wraparound -- see
                                             # SustainedResult's docstring for why)

    left = core_start
    while left - 1 >= 0 and result.bin_p_values[left - 1] < extension_alpha:
        left -= 1
    right = core_end
    while right + 1 < n_bins and result.bin_p_values[right + 1] < extension_alpha:
        right += 1

    # "Recovered" means the next bin was measured and looked normal. A bin
    # with no data is not evidence of recovery, so say so rather than
    # reporting a recovery the data cannot support.
    recovered = right < n_bins - 1
    if recovered and np.isnan(result.bin_p_values[right + 1]):
        warnings.warn(
            f"The bin after the estimated departure (bin "
            f"{int(result.bins[right + 1])}) has no treatment data, so "
            f"'recovered' here means 'the departure stopped being traceable', "
            f"not 'the variable was measured back at normal'.",
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
# Envelope for plotting
# ---------------------------------------------------------------------------

def get_envelope(fit: TideFit, lower_q: float = 0.1, upper_q: float = 0.9,
                  method: str = "bootstrap", year_index: "int | None" = None,
                  rng=None) -> dict:
    """method="bootstrap" (default): percentiles across many Monte Carlo
    synthetic years, using the same resampling as test_treatment -- smooth,
    and backed by the same mechanism that produces the p-value.
    method="empirical": raw quantiles from just the N historical years --
    coarser, but literally just your own data with no resampling involved.
    Individual historical year curves are returned either way for overlay.

    year_index: which year's trend basis to express the envelope on, as a
    calendar-year offset from the first historical year (see
    fit.year_index_for). This matters ONLY when a trend correction was
    applied: the fit detrends every historical year back to the first
    year's level, so an un-shifted envelope sits at that first year's
    level while a treatment curve is plotted at its own year's level.
    Overlaying the two then shows a gap that is entirely the trend --
    on a 15-year record drifting 0.6/year, a perfectly normal 2015 (p=0.40)
    plotted above the band in all 13 bins. Pass the treatment period's
    year_index (plot_envelope does this automatically from the result) to
    put both on the same basis. None leaves the envelope on the first
    historical year's basis.
    """
    rng = _resolve_rng(rng)
    n_bins = len(fit.bins)
    if method == "bootstrap":
        synth = np.array([
            _draw_synthetic_residual(fit.block_pool, n_bins, fit.block_length_bins, rng)
            + fit.median_curve
            for _ in range(fit.config.n_bootstrap)
        ])
        lower = np.percentile(synth, lower_q * 100, axis=0)
        upper = np.percentile(synth, upper_q * 100, axis=0)
    else:
        lower = np.nanquantile(fit.historical_matrix, lower_q, axis=0)
        upper = np.nanquantile(fit.historical_matrix, upper_q, axis=0)

    median_curve = fit.median_curve
    year_curves = fit.historical_matrix
    shift = 0.0
    if fit.trend["applied"] and year_index is not None:
        shift = fit.trend["slope_per_year"] * float(year_index)
        lower, upper = lower + shift, upper + shift
        median_curve = median_curve + shift
        year_curves = year_curves + shift

    if fit.trend["log_space"]:
        lower, upper = np.exp(lower), np.exp(upper)
        year_curves = np.exp(year_curves)
        median_curve = np.exp(median_curve)

    return {
        "bins": fit.bins, "lower": lower, "upper": upper,
        "median_curve": median_curve, "year_curves": year_curves,
        "years_used": fit.years_used, "n_years": len(fit.years_used),
        "method": method, "lower_q": lower_q, "upper_q": upper_q,
        "year_index": year_index, "trend_shift": float(shift),
    }
