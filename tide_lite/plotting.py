"""
tide_lite -- one plot per mode.

  plot_envelope    -- Standard: the period against the significance band.
  plot_sequential  -- Sequential: every event against the historical range.
  plot_cumulative  -- Cumulative: effect size over successive windows.
  plot_monitoring  -- MonitoringSeries: every check against its threshold.

Each takes an optional matplotlib `ax` and returns it; none calls plt.show().
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from .engine import (TideFit, TideResult, get_envelope, bin_significant,
                     significance_band, DepartureRecovery)

_MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
           "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]


def _legend_outside(ax, **kw):
    ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0, **kw)


def _new_ax(figsize=(11, 5)):
    _, ax = plt.subplots(figsize=figsize)
    return ax


def _fmt_p(p, n_reference):
    """A bootstrap p-value, never printed as 0: at the floor 1/(n+1) it
    means 'this small or smaller'."""
    floor = 1.0 / (n_reference + 1)
    if p <= floor * (1 + 1e-9):
        return f"p≤{floor:.2g}"
    return f"p={p:.3f}" if p >= 0.001 else f"p={p:.2g}"


def _bin_name(fit: TideFit, b) -> str:
    return _MONTHS[int(b) - 1] if fit.config.bin_days == "month" else str(int(b))


def _label_bins(ax, fit: TideFit):
    x = np.asarray(fit.bins)
    if fit.config.bin_days == "month":
        ax.set_xticks(x)
        ax.set_xticklabels([_MONTHS[b - 1] for b in x])
        ax.set_xlabel("month")
    else:
        ax.set_xlabel(f"bin ({fit.config.bin_days}-day)")


def plot_envelope(fit: TideFit, result: TideResult, alpha: float = 0.05,
                  treatment_label: str = "treatment period",
                  mark_significant_bins: bool = True,
                  departure: "DepartureRecovery | None" = None,
                  show_history: bool = True, ax=None):
    """Standard: the treatment curve against the significance band.

    The shaded band is significance_band(): a normal period stays entirely
    inside it with probability 1 - alpha. The treatment curve leaving the
    band is exactly the same statement as p <= alpha, and the starred bins
    are exactly the bins bin_significant() flags. Thin lines are the
    historical years (moved to the tested year's trend level if a trend was
    removed). departure= shades an estimate_departure_recovery() span.
    """
    ax = ax or _new_ax()
    band = significance_band(fit, result, alpha)
    x = np.asarray(fit.bins)
    if show_history:
        env = get_envelope(fit, method="empirical", year_index=result.treatment_year_index)
        for i, row in enumerate(env["year_curves"]):
            ax.plot(x, row, color="steelblue", alpha=0.25, linewidth=1,
                    label="historical years" if i == 0 else None)
    lo = np.where(np.isfinite(band["lower"]), band["lower"], np.nan)
    hi = np.where(np.isfinite(band["upper"]), band["upper"], np.nan)
    if result.alternative == "greater":
        lo = np.full_like(hi, np.nanmin(np.r_[hi, band["centre"], result.treatment_curve]))
    elif result.alternative == "less":
        hi = np.full_like(lo, np.nanmax(np.r_[lo, band["centre"], result.treatment_curve]))
    if np.isfinite(band["critical_value"]):
        ax.fill_between(x, lo, hi, color="steelblue", alpha=0.18,
                        label=f"significance band (alpha={alpha}, {result.alternative})")
    ax.plot(x, band["centre"], color="steelblue", linewidth=2, label="historical median")
    ax.plot(x, result.treatment_curve, color="firebrick", linewidth=2.5, marker="o", markersize=4,
            label=f"{treatment_label} ({_fmt_p(result.p_value, result.n_reference)})")
    if mark_significant_bins:
        sig = bin_significant(result, alpha=alpha)
        if np.any(sig):
            ax.scatter(x[sig], result.treatment_curve[sig], marker="*", s=220,
                       color="black", zorder=5, label="outside band (significant)")
    if departure is not None:
        ax.axvspan(departure.departure_start_bin - 0.5, departure.departure_end_bin + 0.5,
                   color="orange", alpha=0.15,
                   label=f"estimated departure ({_bin_name(fit, departure.departure_start_bin)}-"
                         f"{_bin_name(fit, departure.departure_end_bin)}, descriptive)")
        if departure.recovered:
            ax.axvline(departure.recovered_at_bin, color="darkgreen", linewidth=1.5, linestyle=":",
                       label=f"back to normal by {_bin_name(fit, departure.recovered_at_bin)}")
    _label_bins(ax, fit)
    ax.set_ylabel(fit.value_col)
    title = f"{treatment_label} vs. historical normal"
    if fit.trend["applied"]:
        title += f"\n(history moved to the tested year's trend level: {fit.trend['slope_per_year']:+.3g}/year"
        title += " in log units)" if fit.trend["log_space"] else ")"
    ax.set_title(title)
    _legend_outside(ax)
    ax.figure.tight_layout()
    return ax


def plot_sequential(fit: TideFit, results, envelope: "dict | None" = None, ax=None):
    """Sequential: every event on one picture, coloured by whether it
    survived Holm-Bonferroni correction. The band is the historical 10th-90th
    percentile per bin, for context only; judge events by their corrected
    result, not by distance from the band."""
    ax = ax or _new_ax()
    idx = {r.result.treatment_year_index for r in results}
    common = next(iter(idx)) if len(idx) == 1 else None
    env = envelope or get_envelope(fit, method="empirical",
                                   year_index=common if common is not None else None)
    x = np.asarray(env["bins"])
    ax.fill_between(x, env["lower"], env["upper"], color="steelblue", alpha=0.15,
                    label="historical 10th-90th percentile")
    ax.plot(x, env["median_curve"], color="steelblue", linewidth=2, label="historical median")
    styles = ["-", "--", "-.", ":"]
    for i, row in enumerate(results):
        flagged = row.significant_corrected
        ax.plot(x, row.result.treatment_curve, linewidth=2, alpha=0.85,
                color="firebrick" if flagged else "gray", linestyle=styles[i % len(styles)],
                label=f"{row.label}: {_fmt_p(row.result.p_value, row.result.n_reference)}, "
                      f"{'flagged' if flagged else 'not flagged'} (corrected)")
    _label_bins(ax, fit)
    ax.set_ylabel(fit.value_col)
    title = "Sequential: each event vs. historical normal"
    if fit.trend["applied"] and common is None:
        title += "\n(events in different years under a trend: compare p-values, not the band)"
    ax.set_title(title)
    _legend_outside(ax)
    ax.figure.tight_layout()
    return ax


def plot_cumulative(results, alpha: float = 0.05, ax=None):
    """Cumulative: effect size across follow-up windows, coloured by
    significance; the dashed line at 0 is the historical typical level.
    Uses the trend-adjusted effect size, the one the p-value tests."""
    ax = ax or _new_ax(figsize=(8, 4.5))
    labels = [r.label for r in results]
    adjusted = any(r.result.effect_size != r.result.effect_size_trend_adjusted for r in results)
    effects = [r.result.effect_size_trend_adjusted for r in results]
    colors = ["firebrick" if r.significant else "steelblue" for r in results]
    x = np.arange(len(results))
    ax.plot(x, effects, color="gray", linewidth=1, zorder=1)
    ax.scatter(x, effects, c=colors, s=70, zorder=2)
    ax.axhline(0, color="black", linewidth=0.8, linestyle="--")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_ylabel("effect size, trend-adjusted (original units)"
                  if adjusted else "effect size (original units)")
    ax.set_title("Cumulative: effect size over time")
    ax.legend(handles=[
        Line2D([0], [0], marker="o", color="w", markerfacecolor="firebrick", markersize=9,
               label=f"significant (p ≤ {alpha})"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="steelblue", markersize=9,
               label="not significant"),
        Line2D([0], [0], color="black", linewidth=0.8, linestyle="--", label="historical typical level"),
    ])
    ax.figure.tight_layout()
    return ax


def plot_monitoring(series, ax=None):
    """MonitoringSeries: each check's p-value against the per-check
    threshold (log scale)."""
    ax = ax or _new_ax(figsize=(8, 4.5))
    if not series.checks:
        ax.text(0.5, 0.5, "no checks yet", ha="center", va="center", transform=ax.transAxes)
        return ax
    x = [c.check_number for c in series.checks]
    ax.scatter(x, [c.p_value for c in series.checks], s=70, zorder=2,
               c=["firebrick" if c.flagged else "steelblue" for c in series.checks])
    ax.axhline(series.alpha_per_check, color="black", linewidth=0.8, linestyle="--")
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels([c.label for c in series.checks], rotation=20, ha="right")
    ax.set_xlabel(f"check (horizon {series.n_years_horizon}, {series.remaining} remaining)")
    ax.set_ylabel("p-value (log scale)")
    ax.set_title("Monitoring history")
    ax.legend(handles=[
        Line2D([0], [0], marker="o", color="w", markerfacecolor="firebrick", markersize=9, label="flagged"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="steelblue", markersize=9, label="not flagged"),
        Line2D([0], [0], color="black", linewidth=0.8, linestyle="--",
               label=f"threshold per check ({series.alpha_per_check:.4g})"),
    ])
    ax.figure.tight_layout()
    return ax
