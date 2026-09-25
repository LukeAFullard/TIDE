"""
tide_lite -- plain-language summary of one result, including everything
needed to reproduce it.
"""

from __future__ import annotations

from dataclasses import asdict

from .engine import TideFit, TideResult, bin_significant, _le

_MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
           "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]


def _bin_names(fit: TideFit, bins) -> str:
    bins = [int(b) for b in bins]
    if fit.config.bin_days == "month":
        return ", ".join(_MONTHS[b - 1] for b in bins)
    return ", ".join(f"bin {b}" for b in bins)


def summarize(fit: TideFit, result: TideResult, alpha: float = 0.05) -> str:
    """A short plain-language report of `result`: what was compared, the
    answer, which bins drove it, the size of the difference, the caveats
    that apply to this data, and the settings needed to reproduce it."""
    from . import __version__
    sig = _le(result.p_value, alpha)
    years = result.treatment_years or []
    direction = {"two-sided": "higher or lower than", "greater": "higher than",
                 "less": "lower than"}[result.alternative]
    first, last = fit.years_used[0], fit.years_used[-1]
    out = [
        f"Question: were {', '.join(map(str, years))} values {direction} the site's "
        f"normal range ({len(fit.years_used)} historical years, {first}-{last})?",
        f"Answer: {'YES' if sig else 'NO'} at alpha={alpha} "
        f"(p = {result.p_value:.4g}; {'at most' if sig else 'more than'} alpha).",
    ]
    flagged = bin_significant(result, alpha)
    if sig:
        out.append(f"Bins outside the normal range: {_bin_names(fit, fit.bins[flagged])}.")
    tested = ~result.missing_bins
    out.append(
        f"Typical difference from the historical level: {result.effect_size_trend_adjusted:+.4g} "
        f"(units of {fit.value_col!r}; median over {int(tested.sum())} tested bins"
        + (", after allowing for the historical trend)" if fit.trend["applied"] else ")")
        + "."
    )
    if result.missing_bins.any():
        out.append(f"Not tested (no data): {_bin_names(fit, fit.bins[result.missing_bins])}.")

    notes = []
    if not sig:
        notes.append("Not significant means no difference was detected, not that none "
                     "exists; check estimate_power for the size of change this record can detect.")
    if fit.trend["applied"]:
        unit = "log units" if fit.trend["log_space"] else f"{fit.value_col}"
        notes.append(f"A historical trend of {fit.trend['slope_per_year']:+.3g} {unit}/year was "
                     f"removed (Mann-Kendall p={fit.trend['mk_p']:.3f}) and projected to the tested year.")
    if len(fit.years_used) < 10:
        notes.append(f"Only {len(fit.years_used)} historical years: only large changes are detectable.")
    if fit.between_year_var_frac >= 0.5:
        notes.append(f"{fit.between_year_var_frac:.0%} of normal variation is whole-year shifts "
                     f"(wet/dry years), so a change affecting the whole year is hard to "
                     f"tell apart from an unusual year.")
    if fit.years_dropped:
        notes.append(f"Historical years left out: {fit.years_dropped}.")
    if fit.bins_dropped:
        notes.append(f"Bins left out: {fit.bins_dropped}.")
    notes.append("Unusual is not the same as caused by a particular activity; "
                 "check other explanations (flow, rainfall, method changes).")
    out += ["Notes:"] + [f"  - {n}" for n in notes]

    settings = ", ".join(f"{k}={v!r}" for k, v in asdict(fit.config).items())
    out.append(
        f"Reproduce with: tide_lite {__version__}; TideConfig({settings}); "
        f"mode={result.mode!r}, alternative={result.alternative!r}, "
        f"rng={result.seed if result.seed is not None else 'NOT SEEDED'}; "
        f"fit fingerprint {fit.fingerprint()}."
    )
    return "\n".join(out)
