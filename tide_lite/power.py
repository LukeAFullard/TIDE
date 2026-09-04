"""
TIDE lite -- power analysis and sensitivity harness (Phase 5).

Both tools sit on top of the validated engine (Phase 2-4) and simulate --
there is no closed-form power formula for this kind of resampling-based
test, so "can this test even detect an effect I'd care about" has to be
answered by generating synthetic data with a known injected effect and
checking the detection rate.
"""

from __future__ import annotations

from dataclasses import replace
import numpy as np
import pandas as pd

from .engine import (
    TideFit, TideConfig, fit_historical, test_treatment,
    _draw_synthetic_residual, _max_standardized_deviation, _resolve_rng,
)


def simulate_synthetic_curve(fit: TideFit, effect_shift: float = 0.0,
                              n_blocks: int = 1, rng=None) -> np.ndarray:
    """Draw a synthetic 'what a normal treatment period could look like',
    optionally with a constant effect injected, using the same
    position-matched circular block bootstrap as the engine's null
    distribution.

    n_blocks: how many independent draws to average (1 = single-year/
    prediction-style; >1 = confidence-style average of several years).
    effect_shift is in the same units as fit.median_curve: log-space units
    if the fit used detrend_mode="log_additive" (roughly a proportional
    change), otherwise original data units (an absolute change).
    rng: None (fresh entropy), an int (seed), or a numpy Generator.
    """
    rng = _resolve_rng(rng)
    n_bins = len(fit.bins)
    draws = [
        _draw_synthetic_residual(fit.block_pool, n_bins, fit.block_length_bins, rng)
        for _ in range(n_blocks)
    ]
    residual = np.mean(draws, axis=0)
    return fit.median_curve + residual + effect_shift


def estimate_power(fit: TideFit, effect_sizes, alpha: float = 0.05,
                    n_simulations: int = 500, n_treatment_years: int = 1,
                    rng=None) -> dict:
    """Monte Carlo power curve: detection rate at each candidate effect
    size. Returns {effect_size: detection_rate}.

    The null distribution doesn't depend on the injected effect, so it's
    built once (fit.config.n_bootstrap draws) and reused across every
    effect size and simulation -- this is the expensive part; keep
    n_simulations modest (a few hundred) unless you have time to spare.
    rng: None (fresh entropy), an int (seed), or a numpy Generator.
    """
    rng = _resolve_rng(rng)
    n_bins = len(fit.bins)

    null_stats = np.array([
        _max_standardized_deviation(
            np.mean([
                _draw_synthetic_residual(fit.block_pool, n_bins, fit.block_length_bins, rng)
                for _ in range(n_treatment_years)
            ], axis=0) + fit.median_curve,
            fit.median_curve, fit.mad,
        )
        for _ in range(fit.config.n_bootstrap)
    ])
    n_ref = len(null_stats)

    power_curve = {}
    for effect in effect_sizes:
        detections = 0
        for _ in range(n_simulations):
            curve = simulate_synthetic_curve(fit, effect_shift=effect,
                                              n_blocks=n_treatment_years, rng=rng)
            stat = _max_standardized_deviation(curve, fit.median_curve, fit.mad)
            p = (np.sum(null_stats >= stat) + 1) / (n_ref + 1)
            if p < alpha:
                detections += 1
        power_curve[float(effect)] = detections / n_simulations
    return power_curve


def estimate_power_sequential(fit: TideFit, effect_sizes, n_events: int,
                               alpha: float = 0.05, n_simulations: int = 500,
                               rng=None) -> dict:
    """Approximate power for one event within a Sequential set of
    n_events, using the Bonferroni bound alpha/n_events as the effective
    per-test alpha. That's the STRICTEST threshold Holm's step-down
    procedure ever applies, so this UNDERSTATES true Holm power (Holm is
    uniformly more powerful than plain Bonferroni) -- treat it as a
    conservative estimate. An exact version would simulate all n_events
    tests jointly; add that only if the conservative estimate isn't
    good enough on its own.
    """
    effective_alpha = alpha / n_events
    return estimate_power(fit, effect_sizes, alpha=effective_alpha,
                           n_simulations=n_simulations, n_treatment_years=1, rng=rng)


def sensitivity_grid(historical_df: pd.DataFrame, treatment_df: pd.DataFrame,
                      date_col: str, value_col: str, base_config: TideConfig,
                      variations: dict, mode: str = "prediction") -> pd.DataFrame:
    """Re-fit and re-test under each single-parameter variation (holding
    everything else at base_config), to check whether a REAL result is
    stable to defensible alternative choices. A Phase 6 tool -- run once
    you have a result to stress-test, not for choosing Phase 1 settings.

    variations: {TideConfig field name: [alternative values to try]}
    e.g. {"bin_days": [7, 30], "detrend_mode": ["additive", "log_additive"]}

    One row per variant: p_value, effect_size, test_statistic, and how
    many historical years survived that variant's completeness filter
    (watch this -- a variant that quietly drops years isn't a fair
    comparison).
    """
    rows = []

    def _run(cfg: TideConfig, tag: str):
        fit = fit_historical(historical_df, date_col, value_col, cfg)
        result = test_treatment(fit, treatment_df, mode=mode)
        rows.append({
            "variant": tag, "p_value": result.p_value,
            "effect_size": result.effect_size,
            "test_statistic": result.test_statistic,
            "n_historical_years": len(fit.years_used),
        })

    _run(base_config, "base")
    for field_name, alt_values in variations.items():
        for val in alt_values:
            cfg = replace(base_config, **{field_name: val})
            _run(cfg, f"{field_name}={val}")

    return pd.DataFrame(rows)
