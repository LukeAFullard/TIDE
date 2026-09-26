"""
tide_lite -- MonitoringSeries: check each new year as it arrives.

Stateful: create it once, call .check() each time a new year of data is
ready (usually a separate script run, once a year), and .save()/.load() in
between.

Method: fixed horizon, flat split. Commit in advance to at most
n_years_horizon checks; each one is tested at alpha_total / n_years_horizon.
By the union (Bonferroni) bound, the chance of EVER flagging a normal year
across the whole horizon is at most alpha_total, if each check's p-value
is accurate at that small threshold. The split is fixed in advance, so no
future data is needed to know today's threshold.

The catch is the "if": each check is tested far out in the tail of the
null distribution (0.0025 for alpha_total=0.05 over 20 years), and a null
built from N historical years is least accurate there. Measured over a
10-check horizon with a 5% budget, the chance of ever flagging a normal
year was up to 12% with 10 years of history and 5-10% with 20-40 years
(METHODS.md section 4). Use 20+ years and n_bootstrap >= 4000, and treat
the budget as approximate.

When the horizon is used up .check() raises; .renew() starts a new horizon
(a deliberate decision, not an automatic one). The saved state records a
fingerprint of the fit, and loading it against a different fit raises:
changing the baseline part-way through would break the pre-committed
guarantee.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
import json
import warnings

import numpy as np
import pandas as pd

from .engine import TideFit, test_treatment, _le, _VALID_ALTERNATIVES


@dataclass
class MonitoringCheck:
    check_number: int
    label: str
    year_index: "float | None"
    alpha_used: float
    p_value: float
    effect_size: float
    effect_size_trend_adjusted: float
    flagged: bool
    calendar_year: "int | None" = None


class MonitoringSeries:
    def __init__(self, fit: TideFit, alpha_total: float = 0.05,
                 n_years_horizon: int = 20, alternative: str = "two-sided"):
        if not 0 < alpha_total < 1:
            raise ValueError(f"alpha_total must be strictly between 0 and 1, got {alpha_total}.")
        if n_years_horizon < 1:
            raise ValueError(f"n_years_horizon must be >= 1, got {n_years_horizon}.")
        if alternative not in _VALID_ALTERNATIVES:
            raise ValueError(f"alternative must be one of {sorted(_VALID_ALTERNATIVES)}.")
        self.fit = fit
        self.alpha_total = alpha_total
        self.n_years_horizon = n_years_horizon
        self.alternative = alternative
        self.checks: "list[MonitoringCheck]" = []

        min_p = 1.0 / (fit.config.n_bootstrap + 1)
        if not _le(min_p, self.alpha_per_check):
            raise ValueError(
                f"This horizon cannot flag anything: alpha_per_check="
                f"{self.alpha_per_check:.6f} is below the smallest p-value "
                f"{fit.config.n_bootstrap} bootstrap draws can produce ({min_p:.6f}). "
                f"Refit with TideConfig(n_bootstrap>={int(np.ceil(10 / self.alpha_per_check))}), "
                f"shorten n_years_horizon, or raise alpha_total."
            )
        if min_p > self.alpha_per_check / 10:
            warnings.warn(
                f"alpha_per_check={self.alpha_per_check:.6f} rests on fewer than 10 "
                f"of the {fit.config.n_bootstrap} bootstrap draws. Refit with "
                f"n_bootstrap >= {int(np.ceil(10 / self.alpha_per_check))} for a stable threshold.",
                UserWarning, stacklevel=2,
            )
        if len(fit.years_used) < 20:
            warnings.warn(
                f"MonitoringSeries is fitted on {len(fit.years_used)} historical years. "
                f"Its per-check threshold ({self.alpha_per_check:.4f}) is far into the "
                f"tail of the null, which a short record cannot pin down: in simulation "
                f"a 5% budget over 10 checks ran at up to 12% with 10 years of history "
                f"(up to about 10% with 20-40). Use 20+ years and treat the budget as "
                f"approximate (METHODS.md section 4).",
                UserWarning, stacklevel=2,
            )

    @property
    def alpha_per_check(self) -> float:
        return self.alpha_total / self.n_years_horizon

    @property
    def remaining(self) -> int:
        return self.n_years_horizon - len(self.checks)

    def check(self, treatment_df: pd.DataFrame, label: "str | None" = None,
              year_index: "int | None" = None, rng=None) -> MonitoringCheck:
        """Test one new calendar year at alpha_per_check and log it. Pass
        rng=<int> to make the logged p-value reproducible."""
        if self.remaining <= 0:
            raise ValueError(
                f"The monitoring horizon of {self.n_years_horizon} checks is used up. "
                f"Call .renew() to start a new horizon; this is deliberately not automatic."
            )
        years = sorted(set(pd.to_datetime(treatment_df[self.fit.date_col]).dt.year.dropna().astype(int)))
        done = {c.calendar_year for c in self.checks if c.calendar_year is not None}
        repeated = sorted(set(years) & done)
        if repeated:
            raise ValueError(f"Year(s) {repeated} have already been checked in this horizon.")
        result = test_treatment(self.fit, treatment_df, mode="prediction",
                                treatment_year_index=year_index,
                                alternative=self.alternative, rng=rng)
        number = len(self.checks) + 1
        record = MonitoringCheck(
            check_number=number,
            label=label or str(result.treatment_years[0]),
            year_index=result.treatment_year_index,
            alpha_used=self.alpha_per_check,
            p_value=result.p_value,
            effect_size=result.effect_size,
            effect_size_trend_adjusted=result.effect_size_trend_adjusted,
            flagged=_le(result.p_value, self.alpha_per_check),
            calendar_year=int(result.treatment_years[0]),
        )
        self.checks.append(record)
        return record

    def renew(self, alpha_total: "float | None" = None,
              n_years_horizon: "int | None" = None) -> "MonitoringSeries":
        """A NEW series with a fresh budget and an empty log. Save or keep
        this one first if you need the previous horizon's record."""
        return MonitoringSeries(
            fit=self.fit,
            alpha_total=self.alpha_total if alpha_total is None else alpha_total,
            n_years_horizon=self.n_years_horizon if n_years_horizon is None else n_years_horizon,
            alternative=self.alternative,
        )

    def to_dict(self) -> dict:
        """Budget, check log and the fit's fingerprint (not the fit itself:
        re-create it each session with fit_historical on the same data and
        settings)."""
        return {
            "alpha_total": self.alpha_total,
            "n_years_horizon": self.n_years_horizon,
            "alternative": self.alternative,
            "fit_fingerprint": self.fit.fingerprint(),
            "checks": [asdict(c) for c in self.checks],
        }

    @classmethod
    def from_dict(cls, fit: TideFit, d: dict) -> "MonitoringSeries":
        saved = d.get("fit_fingerprint")
        if saved is not None and saved != fit.fingerprint():
            raise ValueError(
                "This monitoring state was created with a different historical fit "
                "(different data, years or settings). Re-create the fit exactly as "
                "before (same historical data and TideConfig), or start a new series."
            )
        if saved is None:
            warnings.warn("Saved state has no fit fingerprint (older version); "
                          "cannot confirm the fit is unchanged.", UserWarning, stacklevel=2)
        series = cls(fit=fit, alpha_total=d["alpha_total"],
                     n_years_horizon=d["n_years_horizon"],
                     alternative=d.get("alternative", "two-sided"))
        series.checks = [MonitoringCheck(**c) for c in d["checks"]]
        return series

    def save(self, path: str) -> None:
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2)

    @classmethod
    def load(cls, fit: TideFit, path: str) -> "MonitoringSeries":
        """fit must be re-created with the same data and settings; see to_dict()."""
        with open(path) as f:
            return cls.from_dict(fit, json.load(f))
