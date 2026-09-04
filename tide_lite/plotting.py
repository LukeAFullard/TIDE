"""
TIDE lite -- plotting (added on request, alongside MonitoringSeries).

One function per mode, each answering that mode's actual question
visually, not just reusing the same envelope picture for everything:

  plot_envelope    -- Standard.    Is this one period unusual?
  plot_sequential  -- Sequential.  Which of these events stand out?
  plot_cumulative  -- Cumulative.  Is it recovering over time?
  plot_monitoring  -- MonitoringSeries. Every check so far, against its bar.

All four take an optional `ax` (a matplotlib Axes) and return it, so you
can compose them into a larger figure or just call .show()/.savefig()
yourself -- these functions never call plt.show().
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from .engine import TideFit, TideResult, get_envelope, bin_significant, DepartureRecovery


def _new_ax(figsize=(9, 5)):
    _, ax = plt.subplots(figsize=figsize)
    return ax


def plot_envelope(fit: TideFit, result: TideResult, envelope: "dict | None" = None,
                   treatment_label: str = "treatment period",
                   mark_significant_bins: bool = True, alpha: float = 0.05,
                   departure: "DepartureRecovery | None" = None, ax=None):
    """Standard: the treatment period's curve against the historical
    envelope. This is the plot from examples/run_example.py, pulled out
    into a reusable function.

    mark_significant_bins: if True (default), marks which specific bins
    were significantly unusual using bin_significant() (Holm-Bonferroni
    corrected across bins) -- the visual answer to "the p-value is for
    the whole year, can we tell month by month," not just a global
    number. Set False to skip (e.g. if you're calling bin_significant()
    yourself with a different alpha and want to avoid the double work).

    departure: an optional estimate_departure_recovery() result -- shades
    the estimated departure span and marks the recovery bin, if any.
    """
    ax = ax or _new_ax()
    env = envelope or get_envelope(fit)
    x = env["bins"]
    for row in env["year_curves"]:
        ax.plot(x, row, color="steelblue", alpha=0.25, linewidth=1)
    ax.fill_between(x, env["lower"], env["upper"], color="steelblue", alpha=0.15,
                     label="10th-90th percentile (bootstrap)")
    ax.plot(x, env["median_curve"], color="steelblue", linewidth=2, label="historical median")
    ax.plot(x, result.treatment_curve, color="firebrick", linewidth=2.5,
            label=f"{treatment_label} (p={result.p_value:.3f})")
    if mark_significant_bins:
        sig = bin_significant(result, alpha=alpha)
        if np.any(sig):
            ax.scatter(np.asarray(x)[sig], result.treatment_curve[sig],
                       marker="*", s=180, color="black", zorder=5,
                       label=f"significant bin (alpha={alpha}, corrected)")
    if departure is not None:
        ax.axvspan(departure.departure_start_bin, departure.departure_end_bin,
                   color="orange", alpha=0.15,
                   label=f"estimated departure (bins {departure.departure_start_bin}-{departure.departure_end_bin})")
        if departure.recovered:
            ax.axvline(departure.recovered_at_bin, color="darkgreen", linewidth=1.5,
                       linestyle=":", label=f"recovered by bin {departure.recovered_at_bin}")
    ax.set_xlabel("bin")
    ax.set_ylabel(fit.value_col)
    ax.set_title("Treatment period vs. historical envelope")
    ax.legend()
    return ax


def plot_sequential(fit: TideFit, results, envelope: "dict | None" = None, ax=None):
    """Sequential: every tested event overlaid on the same historical
    envelope, colored by whether it survived Holm-Bonferroni correction --
    the direct visual answer to "which of these stand out".

    results: the list returned by sequential_test().
    """
    ax = ax or _new_ax()
    env = envelope or get_envelope(fit)
    x = env["bins"]
    ax.fill_between(x, env["lower"], env["upper"], color="steelblue", alpha=0.15,
                     label="10th-90th percentile (bootstrap)")
    ax.plot(x, env["median_curve"], color="steelblue", linewidth=2, label="historical median")

    sig_labeled, nonsig_labeled = False, False
    for row in results:
        if row.significant_corrected:
            color, label = "firebrick", "flagged (corrected)"
            already = sig_labeled
            sig_labeled = True
        else:
            color, label = "gray", "not flagged"
            already = nonsig_labeled
            nonsig_labeled = True
        ax.plot(x, row.result.treatment_curve, color=color, linewidth=2, alpha=0.85,
                 label=None if already else label)
        ax.annotate(row.label, (x[-1], row.result.treatment_curve[-1]),
                    fontsize=8, color=color, xytext=(4, 0), textcoords="offset points")

    ax.set_xlabel("bin")
    ax.set_ylabel(fit.value_col)
    ax.set_title("Sequential: every event vs. historical envelope")
    ax.legend()
    return ax


def plot_cumulative(results, alpha: float = 0.05, ax=None):
    """Cumulative: effect size across successive follow-up windows -- the
    direct visual answer to "is it recovering". Markers colored by
    whether that window was still significant; a line at 0 marks "back
    to the historical median."

    results: the list returned by cumulative_test().
    """
    ax = ax or _new_ax(figsize=(8, 4.5))
    labels = [r.label for r in results]
    effects = [r.result.effect_size for r in results]
    colors = ["firebrick" if r.significant else "steelblue" for r in results]
    x = np.arange(len(results))

    ax.plot(x, effects, color="gray", linewidth=1, zorder=1)
    ax.scatter(x, effects, c=colors, s=70, zorder=2)
    ax.axhline(0, color="black", linewidth=0.8, linestyle="--", label="historical median (fully recovered)")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_ylabel("effect size (original units)")
    ax.set_title("Cumulative: effect size over time")

    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="firebrick", markersize=9, label="still significant"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="steelblue", markersize=9, label="not significant"),
        Line2D([0], [0], color="black", linewidth=0.8, linestyle="--", label="historical median"),
    ]
    ax.legend(handles=handles)
    ax.figure.tight_layout()
    return ax


def plot_monitoring(series, ax=None):
    """MonitoringSeries: every check so far, p-value against the flat
    per-check threshold it had to clear. Log scale on p-value since
    alpha_per_check is often small.
    """
    ax = ax or _new_ax(figsize=(8, 4.5))
    if not series.checks:
        ax.text(0.5, 0.5, "no checks yet", ha="center", va="center", transform=ax.transAxes)
        return ax
    x = [c.check_number for c in series.checks]
    p = [c.p_value for c in series.checks]
    colors = ["firebrick" if c.flagged else "steelblue" for c in series.checks]
    labels = [c.label for c in series.checks]

    ax.scatter(x, p, c=colors, s=70, zorder=2)
    ax.axhline(series.alpha_per_check, color="black", linewidth=0.8, linestyle="--",
               label=f"alpha per check ({series.alpha_per_check:.4f})")
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_xlabel(f"check (horizon: {series.n_years_horizon}, {series.remaining} remaining)")
    ax.set_ylabel("p-value (log scale)")
    ax.set_title("Monitoring history")

    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="firebrick", markersize=9, label="flagged"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="steelblue", markersize=9, label="not flagged"),
        Line2D([0], [0], color="black", linewidth=0.8, linestyle="--", label="alpha per check"),
    ]
    ax.legend(handles=handles)
    ax.figure.tight_layout()
    return ax
