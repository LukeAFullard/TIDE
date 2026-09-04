"""
TIDE lite -- MonitoringSeries: ongoing, open-ended yearly monitoring.

Different in kind from Sequential/Cumulative, not just in name. Those two
take a fixed, fully-known list of events/windows and return all results
in one call -- a batch operation. This is stateful and incremental: you
create it once, then call .check() each time a new year's data is ready,
which in real use happens across separate script runs, months or years
apart. That's why this is a class with persistence (save/load), not a
function.

Method: fixed horizon, flat alpha split. Commit to a horizon of
n_years_horizon checks; each one gets alpha_total / n_years_horizon,
decided in advance. Unlike Holm-Bonferroni (Sequential), no future data
is needed to know today's threshold -- the flat split is valid the
instant each check happens. By a union bound over at most
n_years_horizon checks, the probability of EVER falsely flagging a truly
normal year across the whole horizon is at most alpha_total *in theory*,
regardless of dependence between checks -- but see the IMPORTANT note
below: in practice, with the historical sample sizes common in real
monitoring, this theoretical guarantee is optimistic. This is the classic
idea behind group sequential trial spending functions (Pocock 1977;
O'Brien-Fleming 1979), simplified to a flat, non-adaptive split -- the
"lite" choice among the three spending strategies discussed for this
feature: it needs no new machinery beyond a threshold calculation, and
gives a guarantee you can say out loud ("at most a 5% chance of ever
falsely flagging across the next 20 years of checks"). It trades away
some power in later years compared to an open-ended, discovery-adaptive
scheme (LORD/SAFFRON-style) -- add that only if a real need for
indefinite, un-plannable-horizon monitoring shows up; it's a materially
bigger piece of new statistical machinery, not a small extension of
this one.

When the horizon is used up, .check() raises rather than silently
continuing without a guarantee. Call .renew() to start a fresh horizon --
a deliberate, visible decision, not an automatic one, since extending a
monitoring program's guarantee is a real statistical choice.

IMPORTANT, found by testing this rather than just deriving it on paper:
this mode needs MORE historical years than Standard/Sequential/
Cumulative for its stated guarantee to actually hold. Splitting an
already-small alpha_total across n_years_horizon checks pushes each
individual check's threshold deep into the tail of the block-bootstrap
null distribution -- and bootstrap tail estimates are a well-known weak
point when the original sample (here, historical years) is small,
regardless of how many times you resample from it. In simulation
(fresh, independent synthetic years tested against the real pipeline --
the same method used throughout this project's own test suite), the
true "ever flagged falsely across the horizon" rate for a nominal 5%
budget (n_years_horizon=10) came out to roughly 17% with 10 historical
years, 7% with 20, and 3% with 30 -- i.e. the stated guarantee is
genuinely unreliable at the historical sample sizes the other three
modes tolerate fine.

A tool that lets you check this empirically against your own real
historical fit (rather than the generic numbers above) was attempted and
deliberately left out of this version: a leave-one-out design was built,
and it looked reasonable, but its result pattern (calibration getting
WORSE with MORE historical years -- backwards) didn't match anything
else found in this project's testing. Investigating traced it to a real
structural issue, not a simple bug: holding out one real year at a time
from a small historical pool and repeating that draw many times (with
replacement) to build a simulated horizon produces heavy repetition of
the same held-out year, which understates the true horizon-wide risk --
worst exactly where the risk is highest. Shipping a self-check tool that
gives false reassurance in the regime it most needs to warn about would
be worse than not having one. Revisit this only with a design that
doesn't share that flaw (e.g. genuinely fresh synthetic data per check,
not repeated real-year holdouts) -- until then, use the historical-year
recommendation above directly: 20+ years for reasonable confidence in
this mode's stated guarantee.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
import json

from .engine import TideFit, test_treatment


@dataclass
class MonitoringCheck:
    check_number: int
    label: str
    year_index: "int | None"
    alpha_used: float
    p_value: float
    effect_size: float
    effect_size_trend_adjusted: float
    flagged: bool


class MonitoringSeries:
    def __init__(self, fit: TideFit, alpha_total: float = 0.05,
                 n_years_horizon: int = 20):
        self.fit = fit
        self.alpha_total = alpha_total
        self.n_years_horizon = n_years_horizon
        self.checks: list[MonitoringCheck] = []

    @property
    def alpha_per_check(self) -> float:
        return self.alpha_total / self.n_years_horizon

    @property
    def remaining(self) -> int:
        return self.n_years_horizon - len(self.checks)

    def check(self, treatment_df, label: "str | None" = None,
              year_index: "int | None" = None, rng=None) -> MonitoringCheck:
        """Test one new year against the historical fit, using this
        horizon's flat alpha allocation. Raises if the horizon is used up
        (call .renew()) or if the fit has an applied trend correction and
        year_index wasn't given (checks aren't assumed adjacent in time,
        same reasoning as Sequential/Cumulative).
        """
        if self.remaining <= 0:
            raise ValueError(
                f"Monitoring horizon of {self.n_years_horizon} checks is "
                f"used up. Call .renew(...) to start a fresh horizon with "
                f"a new budget -- this is deliberately not automatic, "
                f"since it's a real statistical decision, not a formality."
            )
        if self.fit.trend["applied"] and year_index is None:
            raise ValueError(
                "The historical fit has an applied trend correction, so "
                "each check needs an explicit year_index -- monitoring "
                "checks aren't assumed to be adjacent in time relative to "
                "the historical fit, the same reasoning as Sequential and "
                "Cumulative."
            )
        check_number = len(self.checks) + 1
        result = test_treatment(self.fit, treatment_df, mode="prediction",
                                 treatment_year_index=year_index, rng=rng)
        record = MonitoringCheck(
            check_number=check_number,
            label=label or f"check {check_number}",
            year_index=year_index,
            alpha_used=self.alpha_per_check,
            p_value=result.p_value,
            effect_size=result.effect_size,
            effect_size_trend_adjusted=result.effect_size_trend_adjusted,
            flagged=bool(result.p_value < self.alpha_per_check),
        )
        self.checks.append(record)
        return record

    def renew(self, alpha_total: "float | None" = None,
              n_years_horizon: "int | None" = None) -> "MonitoringSeries":
        """Start a fresh monitoring horizon with a new budget. Returns a
        NEW MonitoringSeries with an empty check log -- renewal is a real
        statistical decision (a fresh guarantee starting now), so it's
        explicit and doesn't carry check history into the new budget's
        bookkeeping. Read/save .checks first (or keep a reference to this
        object) if you want the prior horizon's record.
        """
        return MonitoringSeries(
            fit=self.fit,
            alpha_total=alpha_total if alpha_total is not None else self.alpha_total,
            n_years_horizon=n_years_horizon if n_years_horizon is not None else self.n_years_horizon,
        )

    def to_dict(self) -> dict:
        """Only the lightweight state -- budget and check history, all
        JSON-safe scalars. Does NOT include the fit itself (a TideFit
        holds numpy arrays and is expected to be re-created fresh each
        session via fit_historical on your historical data, or persisted
        separately by your own method if you want to skip re-fitting).
        """
        return {
            "alpha_total": self.alpha_total,
            "n_years_horizon": self.n_years_horizon,
            "checks": [asdict(c) for c in self.checks],
        }

    @classmethod
    def from_dict(cls, fit: TideFit, d: dict) -> "MonitoringSeries":
        series = cls(fit=fit, alpha_total=d["alpha_total"],
                     n_years_horizon=d["n_years_horizon"])
        series.checks = [MonitoringCheck(**c) for c in d["checks"]]
        return series

    def save(self, path: str) -> None:
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2)

    @classmethod
    def load(cls, fit: TideFit, path: str) -> "MonitoringSeries":
        """fit must be supplied fresh (re-run fit_historical on your
        historical data) -- see to_dict()."""
        with open(path) as f:
            d = json.load(f)
        return cls.from_dict(fit, d)
