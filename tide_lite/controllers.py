"""
TIDE lite -- Sequential and Cumulative controllers (Phase 4).

Both are thin orchestration layers over engine.test_treatment: they call it
repeatedly and differ only in whether/how they correct for doing that.

    Standard    -- one call.  "Was this one period unusual?"
    Sequential  -- many independent calls, Holm-Bonferroni corrected.
                   "Which of these several events were real?"
    Cumulative  -- many calls against ONE fixed baseline, uncorrected.
                   "How is this one known event progressing?"

If you're tempted to add logic here beyond looping + (correct or don't),
it probably belongs in engine.py instead -- these two classes are meant to
stay thin.
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np
import pandas as pd

from .engine import (
    TideFit, TideResult, test_treatment, holm_bonferroni, _resolve_rng,
)


# ---------------------------------------------------------------------------
# Holm-Bonferroni (Holm 1979) now lives in engine.py -- it's used there too
# (bin_significant(), correcting across bins within one period) as well as
# here (correcting across events). Imported above, re-exported here
# unchanged so existing `from tide_lite.controllers import holm_bonferroni`
# still works.
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Sequential -- independent events, corrected together
# ---------------------------------------------------------------------------

@dataclass
class SequentialEvent:
    label: str
    treatment_df: pd.DataFrame
    year_index: int | None = None  # for trend extrapolation; see engine.test_treatment


@dataclass
class SequentialResultRow:
    label: str
    result: TideResult
    significant_raw: bool          # p < alpha, uncorrected
    significant_corrected: bool    # survives Holm-Bonferroni


def sequential_test(fit: TideFit, events: list[SequentialEvent],
                     alpha: float = 0.05, rng=None) -> list[SequentialResultRow]:
    """Test each event independently against the same historical fit, then
    apply Holm-Bonferroni across the whole set.

    Caller's responsibility (Phase 1 non-overlap rule): make sure the
    date ranges behind `events` don't overlap -- one event's data leaking
    into another's window silently breaks independence. This function
    does not check for overlap because it only receives already-aggregated
    per-event dataframes, not a shared raw timeline.

    rng: None (fresh entropy), an int (seed), or a numpy Generator. Pass a
    seed for anything you intend to report -- every p-value here comes from
    a Monte Carlo null, so two runs of the same data differ in the third
    decimal, and a borderline event can flip which side of alpha it lands
    on between runs. One Generator is threaded through every event so the
    whole family is reproducible together.
    """
    rng = _resolve_rng(rng)
    if fit.trend["applied"]:
        missing = [ev.label for ev in events if ev.year_index is None]
        if missing:
            raise ValueError(
                f"The historical fit has an applied trend correction, so "
                f"every SequentialEvent needs an explicit year_index -- "
                f"Sequential events are not assumed to be adjacent in time, "
                f"so they can't share a default. Missing year_index for: "
                f"{missing}"
            )
    raw = [
        (ev.label, test_treatment(fit, ev.treatment_df, mode="prediction",
                                   treatment_year_index=ev.year_index, rng=rng))
        for ev in events
    ]
    p_values = np.array([r.p_value for _, r in raw])
    corrected = holm_bonferroni(p_values, alpha)
    return [
        SequentialResultRow(
            label=label, result=result,
            significant_raw=bool(result.p_value < alpha),
            significant_corrected=bool(reject),
        )
        for (label, result), reject in zip(raw, corrected)
    ]


# ---------------------------------------------------------------------------
# Cumulative -- one event, tracked over time against a fixed baseline
# ---------------------------------------------------------------------------

@dataclass
class CumulativeWindow:
    label: str
    treatment_df: pd.DataFrame
    year_index: int | None = None


@dataclass
class CumulativeResultRow:
    label: str
    result: TideResult
    significant: bool  # plain alpha, deliberately uncorrected -- see module docstring


def cumulative_test(fit: TideFit, windows: list[CumulativeWindow],
                     alpha: float = 0.05, rng=None) -> list[CumulativeResultRow]:
    """Test each window against the SAME fit (same baseline, no re-fitting,
    no correction). fit is not mutated or refit between calls -- reusing
    one TideFit across every window is what makes this "cumulative" rather
    than a fresh Standard test each time.

    rng: None (fresh entropy), an int (seed), or a numpy Generator. Seed it
    for reportable results -- see sequential_test. It matters more here:
    first_recovery() keys off which windows cleared alpha, so an unseeded
    borderline window can move the declared recovery point between runs.
    """
    rng = _resolve_rng(rng)
    if fit.trend["applied"]:
        missing = [w.label for w in windows if w.year_index is None]
        if missing:
            raise ValueError(
                f"The historical fit has an applied trend correction, so "
                f"every CumulativeWindow needs an explicit year_index -- "
                f"without it every window would silently extrapolate the "
                f"trend to the same point. Missing year_index for: {missing}"
            )
    return [
        CumulativeResultRow(
            label=w.label,
            result=(r := test_treatment(fit, w.treatment_df, mode="prediction",
                                         treatment_year_index=w.year_index, rng=rng)),
            significant=bool(r.p_value < alpha),
        )
        for w in windows
    ]


@dataclass
class RecoveryResult:
    recovered: bool
    recovered_at_index: "int | None"   # index into the results list where the
                                        # qualifying non-significant run STARTS
    recovered_at_label: "str | None"


def first_recovery(results: list[CumulativeResultRow],
                    consecutive_required: int = 2) -> RecoveryResult:
    """Scans Cumulative results IN ORDER and finds the first point where
    `consecutive_required` windows in a row are all non-significant.

    Why not just the first single non-significant window: a lone
    non-significant check-in can easily be noise rather than genuine
    recovery -- this project's own testing found real fits can show a
    false-clear on any single check a meaningful fraction of the time
    (see tests/test_engine.py's calibration notes). Requiring several in
    a row before declaring "recovered" is a simple, direct guard against
    that, at the cost of taking longer to declare recovery than eyeballing
    the plot would.

    consecutive_required=1 recovers the naive "first non-significant
    window" behavior if you want it, but 2+ is recommended for anything
    you'd actually report.
    """
    if consecutive_required < 1:
        raise ValueError(
            f"consecutive_required must be >= 1, got {consecutive_required}. "
            f"With 0 the rule is vacuous -- an empty run trivially satisfies "
            f"'all non-significant', so it declared recovery at the first "
            f"window no matter how significant that window was."
        )
    n = len(results)
    for i in range(n - consecutive_required + 1):
        window = results[i:i + consecutive_required]
        if all(not r.significant for r in window):
            return RecoveryResult(recovered=True, recovered_at_index=i,
                                   recovered_at_label=results[i].label)
    return RecoveryResult(recovered=False, recovered_at_index=None, recovered_at_label=None)
