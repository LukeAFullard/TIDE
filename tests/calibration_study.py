"""
Calibration study: how often does the test flag a year that is genuinely
normal? This is the evidence behind the error-rate tables in METHODS.md
section 4. Not part of the pytest suite (it takes several minutes); run:

    python tests/calibration_study.py            # full study, ~5 min on 4 cores
    python tests/calibration_study.py --quick    # fewer trials, ~1-2 min

How it works: simulate a site's history from a known process, fit it, then
test brand-new years drawn from the SAME process (so every flag is a false
alarm). Repeat over many simulated histories and count. A well-calibrated
test at alpha=0.05 flags about 5% of these years.

The simulated site: a seasonal cycle (sine, amplitude 3 around 10) plus
day-to-day noise that persists for days (AR(1), phi=0.7, sd=0.6), plus,
where stated, a whole-year shift drawn once per year ("wet/dry years") with
standard deviation year_sd. For scale: the typical year-to-year wobble of
one month's median from the day-to-day noise alone is WOBBLE = 0.37 (SD,
measured), so year_sd=0.5 (about 1.35 x WOBBLE) is a moderate and
year_sd=1.2 (about 3.2 x WOBBLE) a dominant whole-year effect.
"""

import argparse
import os
import sys
import time
import warnings
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.signal import lfilter

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from tide_lite import TideConfig, fit_historical, test_treatment, MonitoringSeries  # noqa: E402


WOBBLE = 0.37


def make_year(rng, year, year_sd=0.0, year_dist="normal", effect=0.0, trend=0.0,
              log=False, grab=False, spike=0.0):
    n = 365
    days = np.arange(1, n + 1)
    noise = lfilter([1.0], [1.0, -0.7], rng.normal(0, 0.6, n))
    if year_dist == "normal":
        a = rng.normal(0, year_sd)
    elif year_dist == "t3":          # heavy-tailed year effects (sd = year_sd)
        a = year_sd * rng.standard_t(3) / np.sqrt(3.0)
    else:                            # right-skewed year effects (centred lognormal, sd = year_sd)
        z = rng.lognormal(0, 0.8)
        a = year_sd * (z - np.exp(0.32)) / np.sqrt((np.exp(0.64) - 1) * np.exp(0.64))
    seasonal = 10.0 + 3.0 * np.sin(2 * np.pi * (days - 60) / 365)
    level = trend * (year - 2000) + spike * ((days >= 121) & (days <= 181))   # May-Jun
    if log:   # multiplicative process, strictly positive, skewed
        values = np.exp(np.log(seasonal) + 0.1 * (noise + a) + level + effect)
    else:
        values = seasonal + noise + a + level + effect
    df = pd.DataFrame({"date": pd.date_range(f"{year}-01-01", periods=n, freq="D"),
                       "value": values})
    if grab:  # one grab sample per calendar month on a random day
        month = df["date"].dt.month.to_numpy()
        pick = [rng.choice(np.flatnonzero(month == m)) for m in range(1, 13)]
        df = df.iloc[pick]
    return df


def _cell(args):
    (seed, n_years, k, kw, cfg_kw, test_kw, conf_k) = args
    warnings.simplefilter("ignore")
    rng = np.random.default_rng(seed)
    kw = dict(kw)
    next_only = kw.pop("next_year_only", False)   # test only the year right after the record
    treat_months = kw.pop("treat_months", 12)     # treatment covers Jan..this month only
    treat_kw = kw.pop("treat_kw", {})             # changes applied to the treatment year only
    hist = pd.concat([make_year(rng, 2000 + i, **kw) for i in range(n_years)], ignore_index=True)
    fit = fit_historical(hist, "date", "value", TideConfig(**cfg_kw))
    out = []
    for i in range(k):
        if conf_k:
            base = 2000 + n_years + i * conf_k
            df = pd.concat([make_year(rng, base + j, **kw) for j in range(conf_k)])
            r = test_treatment(fit, df, mode="confidence", rng=rng, **test_kw)
        else:
            year = 2000 + n_years + (0 if next_only else i)
            t = make_year(rng, year, **{**kw, **treat_kw})
            t = t[t["date"].dt.month <= treat_months]
            r = test_treatment(fit, t, rng=rng, **test_kw)
        out.append([r.p_value] + [r.sustained[s].p_value for s in sorted(r.sustained)])
    return out


def run_cell(pool, n_fits, k, n_years, kw=None, cfg_kw=None, test_kw=None, conf_k=0, seed0=0):
    args = [(seed0 * 100003 + s, n_years, k, kw or {}, cfg_kw or {}, test_kw or {}, conf_k)
            for s in range(n_fits)]
    return np.array([row for chunk in pool.map(_cell, args) for row in chunk])


def _monitor(args):
    seed, n_years, horizon, n_series, kw = args
    warnings.simplefilter("ignore")
    rng = np.random.default_rng(seed)
    hist = pd.concat([make_year(rng, 2000 + i, **kw) for i in range(n_years)], ignore_index=True)
    fit = fit_historical(hist, "date", "value", TideConfig(n_bootstrap=4000))
    flagged = []
    for s in range(n_series):
        series = MonitoringSeries(fit, alpha_total=0.05, n_years_horizon=horizon)
        start = 2000 + n_years + s * horizon
        hit = any(series.check(make_year(rng, start + j, **kw), rng=rng).flagged
                  for j in range(horizon))
        flagged.append(hit)
    return flagged


def fmt(p, alphas=(0.10, 0.05, 0.01)):
    return "  ".join(f"{np.mean(p <= a * (1 + 1e-9)):.3f}" for a in alphas)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--only", default="ABCDE", help="sections to run, e.g. --only D")
    args = ap.parse_args()
    n_fits, k = (60, 10) if args.quick else (150, 10)
    t0 = time.time()
    with Pool(min(4, os.cpu_count() or 1)) as pool:
        print(f"False-alarm rate (share of normal years flagged), {n_fits * k} trials per row.")
        print("Columns: nominal alpha = 0.10  0.05  0.01\n")

        print("A. Monthly bins (default), by record length and whole-year effect size")
        for n_years in ((5, 10, 20, 40) if "A" in args.only else ()):
            for ysd in (0.0, 0.5, 1.2):
                p = run_cell(pool, n_fits, k, n_years, {"year_sd": ysd}, seed0=n_years * 10 + int(ysd * 10))
                print(f"   {n_years:2d} years, year_sd={ysd:3.1f}:  {fmt(p[:, 0])}", flush=True)

        print("\nB. Non-normal whole-year effects (year_sd=1.2)")
        for dist in (("t3", "skew") if "B" in args.only else ()):
            for n_years in (10, 20):
                p = run_cell(pool, n_fits, k, n_years, {"year_sd": 1.2, "year_dist": dist},
                             seed0=500 + n_years + len(dist))
                print(f"   {n_years:2d} years, {dist:5s}:  {fmt(p[:, 0])}", flush=True)

        print("\nC. Other settings and data shapes (year_sd=0.5, 10 and 20 years)")
        cases = {
            "real trend 0.05/yr, next year": ({"year_sd": 0.5, "trend": 0.05, "next_year_only": True}, {}, {}, 0),
            "log_additive, skewed data": ({"year_sd": 0.5, "log": True}, {"detrend_mode": "log_additive"}, {}, 0),
            "monthly grab samples": ({"year_sd": 0.5, "grab": True}, {}, {}, 0),
            "part year (Jan-Jun only)": ({"year_sd": 0.5, "treat_months": 6}, {}, {}, 0),
            "weekly bins (bin_days=7)": ({"year_sd": 0.5}, {"bin_days": 7}, {}, 0),
            "one-sided (greater)": ({"year_sd": 0.5}, {}, {"alternative": "greater"}, 0),
            "confidence mode, 3 years": ({"year_sd": 0.5}, {}, {}, 3),
            "sustained test, run_length=3": ({"year_sd": 0.5}, {}, {"run_lengths": [3]}, 0),
        }
        for j, (name, (kw, cfg_kw, test_kw, conf_k)) in enumerate(
                cases.items() if "C" in args.only else []):
            for n_years in (10, 20):
                p = run_cell(pool, n_fits, k, n_years, kw, cfg_kw, test_kw, conf_k,
                             seed0=900 + 10 * j + n_years)
                col = -1 if "run_lengths" in test_kw else 0
                print(f"   {name:30s} {n_years:2d} years:  {fmt(p[:, col])}", flush=True)

        print("\nD. Power: share of CHANGED years flagged at alpha=0.05 (whole-year shift, or")
        print(f"   a May-Jun spike), in multiples of the typical monthly wobble ({WOBBLE})")
        for n_years in ((10, 20) if "D" in args.only else ()):
            for ysd in (0.0, 0.5):
                cells = []
                for name, tkw in (("shift 1x", {"effect": WOBBLE}), ("shift 2x", {"effect": 2 * WOBBLE}),
                                  ("shift 4x", {"effect": 4 * WOBBLE}), ("spike 4x", {"spike": 4 * WOBBLE}),
                                  ("spike 8x", {"spike": 8 * WOBBLE})):
                    p = run_cell(pool, n_fits // 3, k, n_years, {"year_sd": ysd, "treat_kw": tkw},
                                 seed0=3000 + n_years + int(10 * ysd))
                    cells.append(f"{name} {np.mean(p[:, 0] <= 0.05):.2f}")
                print(f"   {n_years:2d} years, year_sd={ysd:3.1f}:  " + "   ".join(cells), flush=True)

        print("\nE. MonitoringSeries: chance of EVER flagging a normal year over a 10-year")
        print("   horizon (alpha_total=0.05, so alpha 0.005 per check), n_bootstrap=4000")
        n_m = 40 if args.quick else 100
        for n_years in ((10, 20, 40) if "E" in args.only else ()):
            for ysd in (0.0, 0.5):
                res = pool.map(_monitor, [(7000 + s * 7 + n_years, n_years, 10, 3, {"year_sd": ysd})
                                          for s in range(n_m)])
                rate = np.mean([x for chunk in res for x in chunk])
                print(f"   {n_years:2d} years, year_sd={ysd:3.1f}:  {rate:.3f}  ({3 * n_m} horizons)",
                      flush=True)
    print(f"\n({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
