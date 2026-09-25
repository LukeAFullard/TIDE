"""
tide_lite -- Sequential and Cumulative: several periods against one fit.

Both call engine.test_treatment once per period. They differ only in
whether they correct for testing several periods:

    Standard    -- one call.  "Was this one period unusual?"
    Sequential  -- several separate periods, Holm-Bonferroni corrected.
                   "Which of these periods were unusual?"
    Cumulative  -- successive periods after one known event, uncorrected.
                   "Is the effect of this event still visible?"
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np
import pandas as pd

from .engine import TideFit, TideResult, test_treatment, holm_bonferroni, _resolve_rng, _le


# ---------------------------------------------------------------------------
# Sequential -- separate periods, corrected together
# ---------------------------------------------------------------------------

@dataclass
class SequentialEvent:
    label: str
    treatment_df: pd.DataFrame
    year_index: "int | None" = None  # override only; default is read from the dates


@dataclass
class SequentialResultRow:
    label: str
    result: TideResult
    significant_raw: bool          # p <= alpha on its own (not corrected)
    significant_corrected: bool    # survives Holm-Bonferroni across all events: report this


def _calendar_years(df: pd.DataFrame, date_col: str) -> set:
    return set(pd.to_datetime(df[date_col]).dt.year.dropna().astype(int))


def sequential_test(fit: TideFit, events: "list[SequentialEvent]", alpha: float = 0.05,
                    alternative: str = "two-sided", rng=None) -> "list[SequentialResultRow]":
    """Test each event (one calendar year each) against the same fit, then
    apply Holm-Bonferroni across all of them, so the chance of flagging ANY
    event in a set of normal years stays at or below alpha.

    Events must be different calendar years at the same site (raises if two
    share a year). For several SITES, fit each site separately and pass the
    p-values to holm_bonferroni().

    rng: seed or Generator, shared across events so the whole set is
    reproducible from one seed.
    """
    rng = _resolve_rng(rng)
    seen = {}
    for ev in events:
        for y in _calendar_years(ev.treatment_df, fit.date_col):
            if y in seen:
                raise ValueError(
                    f"Events {seen[y]!r} and {ev.label!r} both contain {y} data. "
                    f"Sequential events must not overlap in time; for several sites, "
                    f"fit each site on its own history."
                )
            seen[y] = ev.label
    raw = [
        (ev.label, test_treatment(fit, ev.treatment_df, mode="prediction",
                                  treatment_year_index=ev.year_index,
                                  alternative=alternative, rng=rng))
        for ev in events
    ]
    corrected = holm_bonferroni(np.array([r.p_value for _, r in raw]), alpha)
    return [
        SequentialResultRow(label=label, result=result,
                            significant_raw=_le(result.p_value, alpha),
                            significant_corrected=bool(reject))
        for (label, result), reject in zip(raw, corrected)
    ]


# ---------------------------------------------------------------------------
# Cumulative -- one event, followed over time against a fixed baseline
# ---------------------------------------------------------------------------

@dataclass
class CumulativeWindow:
    label: str
    treatment_df: pd.DataFrame
    year_index: "int | None" = None  # override only; default is read from the dates


@dataclass
class CumulativeResultRow:
    label: str
    result: TideResult
    significant: bool  # p <= alpha, deliberately uncorrected (see module docstring)


def cumulative_test(fit: TideFit, windows: "list[CumulativeWindow]", alpha: float = 0.05,
                    alternative: str = "two-sided", rng=None) -> "list[CumulativeResultRow]":
    """Test each follow-up window (one calendar year each, in time order)
    against the same fit, with no multiple-testing correction. Use it to
    describe how one already-established event evolves; it does not
    control the chance of a false alarm across windows (use sequential_test
    for that).

    rng: seed or Generator, shared across windows. Seed anything you
    report: first_recovery() depends on which windows cleared alpha.
    """
    rng = _resolve_rng(rng)
    rows = []
    for w in windows:
        r = test_treatment(fit, w.treatment_df, mode="prediction",
                           treatment_year_index=w.year_index,
                           alternative=alternative, rng=rng)
        rows.append(CumulativeResultRow(label=w.label, result=r,
                                        significant=_le(r.p_value, alpha)))
    return rows


@dataclass
class RecoveryResult:
    recovered: bool
    recovered_at_index: "int | None"   # index of the first window of the run
    recovered_at_label: "str | None"


def first_recovery(results: "list[CumulativeResultRow]",
                   consecutive_required: int = 2) -> RecoveryResult:
    """First point where `consecutive_required` windows in a row are all
    non-significant (a single clear window is easily chance).

    Read "recovered" as "no longer distinguishable from normal", not as
    proof of a return to normal: a test that cannot see the effect (low
    power, see estimate_power) also reports non-significant windows.
    """
    if consecutive_required < 1:
        raise ValueError(
            f"consecutive_required must be >= 1, got {consecutive_required}. With 0 "
            f"the rule is vacuous and declares recovery at the first window."
        )
    n = len(results)
    for i in range(n - consecutive_required + 1):
        if all(not r.significant for r in results[i:i + consecutive_required]):
            return RecoveryResult(recovered=True, recovered_at_index=i,
                                  recovered_at_label=results[i].label)
    return RecoveryResult(recovered=False, recovered_at_index=None, recovered_at_label=None)
