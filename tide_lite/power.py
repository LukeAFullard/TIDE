"""
tide_lite -- power analysis and sensitivity check.

estimate_power: "could this test detect a change of a given size with the
history I have?" Answered by simulation: synthetic normal years from the
fit, with a known shift added, run through the same test.

sensitivity_grid: "does my result hold under other reasonable settings?"
Re-fits and re-tests with one setting changed at a time.
"""

from __future__ import annotations

from dataclasses import replace
import warnings

import numpy as np
import pandas as pd

from .engine import (
    TideFit, TideConfig, fit_historical, test_treatment,
    _null_scores, _directional, _resolve_rng, _le, _VALID_ALTERNATIVES,
)


def _next_year_index(fit: TideFit) -> float:
    return float(fit.years_used[-1] - fit.first_year + 1)


def simulate_synthetic_curve(fit: TideFit, effect_shift: float = 0.0,
                             n_blocks: int = 1, rng=None) -> np.ndarray:
    """One synthetic normal period (the average of n_blocks years) from the
    test's own null, plus a constant shift. Returned in the fit's working
    units (natural logs if detrend_mode="log_additive"), on the trend level
    of the first historical year; effect_shift is in the same units."""
    rng = _resolve_rng(rng)
    idx = [_next_year_index(fit)] * n_blocks
    z = _null_scores(fit, 1, rng, idx)[0]
    return fit.median_curve + z * fit.mad + effect_shift


def estimate_power(fit: TideFit, effect_sizes, alpha: float = 0.05,
                   n_simulations: int = 500, n_treatment_years: int = 1,
                   alternative: str = "two-sided", rng=None) -> dict:
    """Share of simulated periods flagged (p <= alpha) at each effect size.
    Returns {effect_size: detection_rate}.

    effect_sizes are constant shifts in the fit's working units: original
    units for detrend_mode "additive"/"none"; natural-log units for
    "log_additive" (0.1 is about +10%). Positive = higher.

    Simulated periods are synthetic normal years from the test's own null
    plus the shift, so the rate at effect 0 is about alpha by construction
    (that is not a calibration check; see METHODS.md section 4). Real
    changes are rarely perfectly uniform, so treat a size this calls
    marginal as not reliably detectable.
    """
    if alternative not in _VALID_ALTERNATIVES:
        raise ValueError(f"alternative must be one of {sorted(_VALID_ALTERNATIVES)}.")
    rng = _resolve_rng(rng)
    idx = [_next_year_index(fit)] * n_treatment_years
    B = fit.config.n_bootstrap
    null_max = _directional(_null_scores(fit, B, rng, idx), alternative).max(axis=1)
    power_curve = {}
    for effect in effect_sizes:
        z = _null_scores(fit, n_simulations, rng, idx) + effect / fit.mad
        T = _directional(z, alternative).max(axis=1)
        p = (np.sum(null_max[None, :] >= T[:, None], axis=1) + 1) / (B + 1)
        power_curve[float(effect)] = float(np.mean([_le(v, alpha) for v in p]))
    return power_curve


def estimate_power_sequential(fit: TideFit, effect_sizes, n_events: int,
                              alpha: float = 0.05, n_simulations: int = 500,
                              alternative: str = "two-sided", rng=None) -> dict:
    """Power for one event in a Sequential set of n_events, using alpha /
    n_events (the strictest threshold Holm ever applies). Holm is never
    less powerful than this, so the true power is at least this high."""
    return estimate_power(fit, effect_sizes, alpha=alpha / n_events,
                          n_simulations=n_simulations, n_treatment_years=1,
                          alternative=alternative, rng=rng)


def sensitivity_grid(historical_df: pd.DataFrame, treatment_df: pd.DataFrame,
                     date_col: str, value_col: str, base_config: TideConfig,
                     variations: dict, mode: str = "prediction",
                     alternative: str = "two-sided",
                     treatment_year_index: "int | None" = None,
                     rng=None) -> pd.DataFrame:
    """Re-fit and re-test with one setting changed at a time, e.g.
    variations={"bin_days": [30, 7], "detrend_mode": ["log_additive"]}.

    Every row uses the SAME random seed (rng, if it is an int), so
    differences between rows come from the setting, not from Monte Carlo
    noise, and the "base" row reproduces test_treatment(..., rng=rng). A variant that cannot run
    (e.g. logs of zero values) gets its error message in the "error"
    column instead of stopping the grid.

    Watch n_historical_years, n_bins and trend_applied: a variant can change
    the answer by dropping years or bins rather than through the setting.
    """
    if isinstance(rng, (int, np.integer)) and not isinstance(rng, bool):
        seed = int(rng)          # base row then matches test_treatment(..., rng=seed)
    else:
        seed = int(_resolve_rng(rng).integers(0, 2**32 - 1))
    rows = []

    def _run(cfg_kwargs: dict, tag: str):
        row = {"variant": tag}
        try:
            cfg = replace(base_config, **cfg_kwargs)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                fit = fit_historical(historical_df, date_col, value_col, cfg)
                result = test_treatment(fit, treatment_df, mode=mode, alternative=alternative,
                                        treatment_year_index=treatment_year_index, rng=seed)
            row.update({
                "p_value": result.p_value,
                "effect_size": result.effect_size,
                "effect_size_trend_adjusted": result.effect_size_trend_adjusted,
                "test_statistic": result.test_statistic,
                "n_historical_years": len(fit.years_used),
                "n_bins": len(fit.bins),
                "trend_applied": fit.trend["applied"],
                "between_year_var_frac": round(fit.between_year_var_frac, 3),
                "error": "",
            })
        except ValueError as e:
            row["error"] = str(e)
        rows.append(row)

    _run({}, "base")
    for field_name, alt_values in variations.items():
        for val in alt_values:
            _run({field_name: val}, f"{field_name}={val}")
    return pd.DataFrame(rows)
