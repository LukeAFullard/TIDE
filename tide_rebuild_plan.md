# TIDE Rebuild Plan

**Revisions so far:**
- All three detector modes — Standard, Sequential, Cumulative — are core v1 scope (real experimental data requires all three).
- Power and sensitivity analysis are also core v1 scope (Phase 5) — concrete need confirmed.
- Audit mode/legal-defensibility and the benchmark suite vs. other methods remain gated (Phase 7) — no concrete need for those yet.
- Phases 0-6 have been implemented, tested against synthetic data, and delivered as working code (`tide_lite_project.zip`) — see **Files** below. Two real findings came out of building and testing it, folded into the design rather than left as caveats: see **Findings from building this** below.
- A 4th mode, MonitoringSeries (ongoing yearly monitoring with no fixed treatment), was added post-Phase-6 on request — see **Phase 8** below, including a real finding about how much historical data it needs and a self-check tool that was built, found flawed, and deliberately not shipped.

## Files (this rebuild)

| Phase | File | Status |
|---|---|---|
| 0-1 | `METHODS_template.md` | Template — fill in before using the rest |
| 2 | `tide_lite/engine.py` | Built, tested |
| 3 | `tests/test_engine.py` | Built, passing (6/6) |
| 4 | `tide_lite/controllers.py` | Built, tested |
| 5 | `tide_lite/power.py` | Built, smoke-tested |
| 6 | `examples/run_example.py` | Built, runs end-to-end, produces `envelope_plot.png` |
| 7 | — | Still not built (still gated, as above) |
| 8 | `tide_lite/monitoring.py`, `tests/test_monitoring.py` | Built, passing (5/5) — see Phase 8 |

`README.md` in the project root is the entry point — setup, how to run the example and tests, layout.

## Findings from building this

Two things emerged from actually testing the code against synthetic data, not just designing it on paper — both are the kind of thing that's much cheaper to catch now than after trusting a result on real data.

1. **The original block design was degenerate.** Resampling whole years as the bootstrap block (one full year = one block, as originally planned) means the null distribution has only N possible values — with N=10 historical years, the smallest achievable p-value is 1/11 ≈ 0.09, which is *above* alpha=0.05. That makes significance mathematically unreachable at any effect size, not just unlikely. Fixed by resampling sub-year (quarter-year, by default) blocks pooled across all historical years and start positions instead — this is the standard block-length tradeoff in the bootstrap literature (shorter blocks trade a little long-range correlation preservation for a null distribution rich enough to actually resolve significance). `engine.py`'s module docstring has the full explanation.
2. **Calibration is noisier per-fit than it looks on paper.** With ~10 historical years, the false-positive rate of one specific fitted model varies a lot fit to fit — 0% to 25% at nominal alpha=0.05 across a 25-fit simulation — even though it averages out correctly across many different possible historical samples. Any one real analysis only ever has one historical fit. This is the small-N power/calibration limitation from the "is this defensible" discussion earlier, now with a number attached — and it's the concrete reason Phase 5 (power + sensitivity analysis) matters in practice, not just in principle. Full simulation in the docstring at the top of `tests/test_engine.py`.

A follow-up adversarial review (`ADVERSARIAL_REVIEW.md`) found and fixed six more issues, including a second, related calibration-diversity problem introduced by fixing a subtler version of finding #1 above — see that file for the full trail.

## Ground rules

Read this before adding anything the plan below doesn't call for.

- Every feature must trace to a concrete, current need — not "would be good practice." If you can't name the need in one sentence, it goes in a `FUTURE_IDEAS.md` parking lot, not the codebase.
- One methods doc, one README, until Phase 6 is done. No new `.md` files before that.
- If you're building this with an AI coding agent again, give it an explicit scope boundary up front: *"Do not add new classes, config options, validation modes, or documentation files unless I explicitly ask. If you think something is missing, list it as a suggestion — don't implement it."*

## Phase 0 — Pin down the actual question(s)

**Deliverable:** an event calendar, plain language.

- List every treatment event: date/window, site, variable.
- For each event, name which question you're asking: "is this one different" (Standard), "is this one of several independently real effects" (Sequential), or "how is this one evolving over time" (Cumulative).
- One explicit answer: does this need to survive external (legal/regulatory) scrutiny, or is it research/monitoring only? This gates Phase 7.

**Gate:** every event in your real data is mapped to Standard, Sequential, or Cumulative (or more than one).

## Phase 1 — Statistical design, on paper, before code

**Deliverable:** `METHODS.md` — keep it tight.

Core, applies to all three:
- **Time resolution:** weekly or monthly bins, not daily.
- **Alignment:** calendar day-of-year (1–365), Feb 29 dropped or folded into day 365. Accept ~5 days of drift over 20 years as documented, not something to rearchitect around.
- **Exchangeability check:** plot all historical years overlaid; confirm by eye they read as one process.
- **Detrending:** additive vs. multiplicative, decided now. For concentration-type variables, check whether variance scales with level — if it does, work in log space.
- **Baseline (V0):** define the window once; confirm it cannot overlap any treatment period.
- **Test statistic and effect-size definition:** the one global statistic the permutation test uses (e.g. max standardized deviation across the year), and exactly how "effect size" is defined in your units. Write both down precisely — Phase 5's power analysis is built directly on this definition, so vagueness here becomes vagueness in every power number later.
- **Envelope mode:** prediction interval for single-year treatments.

Sequential-specific:
- **Correction procedure, spelled out exactly:** Holm-Bonferroni is a step-down procedure — sort p-values ascending, compare the smallest to α/n, the next smallest to α/(n−1), and so on — not naive Bonferroni (α/n applied to every test uniformly), which is what gets implemented by mistake and is needlessly conservative. Write down which exact set of tests are corrected together.
- **Non-overlap rule:** how irregularly-spaced treatment windows get sliced so one event's data never leaks into another's test window.

Cumulative-specific:
- **Baseline-reuse rule:** one `fit_historical()` call, reused unchanged across every tracking window. Each window's "time since event" is measured from the same original event — never re-anchored to a later window.

**Gate:** `METHODS.md` covers all three modes' specific rules, and every implementation decision below traces back to a line in it.

## Phase 2 — Core engine (build this once, correctly)

**Deliverable:** one class, doing exactly four things:

1. `fit_historical()` — align years, apply the Phase 1 detrend method, build the historical matrix.
2. Bootstrap — resample year-blocks (Monte Carlo for ≥6 historical years, exact/exhaustive enumeration below that), build the envelope.
3. `test_treatment()` — align the treatment window, compute the Phase 1 test statistic, compare against the bootstrap null, return a p-value and effect size.
4. One plot — treatment trajectory against the envelope.

Load-bearing for all three modes and for the power/sensitivity tooling in Phase 5 — the single highest-value phase to get right.

**Explicitly not here:** audit mode, normalization/envelope modes beyond what Phase 1 decided, a sprawling config object.

**Gate:** runs correctly on synthetic data.

## Phase 3 — Lean validation of the core engine

**Deliverable:** 5–10 scripted checks, not a full suite:

- Null case: no injected effect → doesn't flag at a rate above α across repeated synthetic draws.
- Known effect: inject an effect of a defined size → it gets detected.
- The Phase 1 exchangeability plot, formalized as a check on your real historical data.

**Gate:** you trust the core engine's number before it becomes the foundation for everything after this.

## Phase 4 — Sequential and Cumulative controllers

Both are thin orchestration layers over the Phase 2 engine.

- **SequentialTIDE:** loop over your Phase 0 event calendar, call the core engine once per event on non-overlapping data, collect p-values, apply Holm-Bonferroni exactly as specified in Phase 1, return per-event corrected results.
- **CumulativeTIDE:** call `fit_historical()` once, then call the core engine repeatedly across advancing windows against that one fixed baseline, no correction applied, return per-window results.

Validation specific to these two:
- Simulate several independent null datasets, run Sequential repeatedly, confirm the family-wise false positive rate lands near the nominal rate.
- Simulate an "effect then recovery" synthetic series, run Cumulative, confirm it flags the early windows and clears on the later ones.

**Gate:** both controllers pass their specific checks against synthetic data with known ground truth.

## Phase 5 — Power and sensitivity analysis

Both are simulation-based tools built on the validated Phase 2–4 machinery — no closed-form power formula exists for a bootstrap/permutation test, so this has to be estimated, not looked up.

**Power analysis:**
1. Using your real historical data, generate many synthetic "null" treatment years via the same block-resampling process used to build the envelope (i.e., what a no-effect year would plausibly look like).
2. For a grid of candidate effect sizes (in the units fixed in Phase 1), inject each into the synthetic treatment years.
3. Run the full pipeline on each and record detection (p < α); repeat several hundred times per effect size.
4. The detection rate at each effect size is your power curve. This turns the generic "small-N power is limited" caveat into an actual number: the smallest effect size you could reliably catch given your real historical sample.
5. Do this separately for Sequential (power is necessarily lower than Standard for the same effect size, because of the Holm-Bonferroni correction — worth quantifying how much you lose per additional event tested) and for Cumulative (power may differ by window position, e.g. a partial early window vs. a full later one).

**Sensitivity analysis:**
- Build a small harness that re-runs the Phase 6 real analyses under a handful of defensible alternative Phase 1 choices: aggregation resolution, detrending method (including log vs. additive), bootstrap iteration count, and baseline window definition if a reasonable alternative exists.
- Tabulate p-value, effect size, and detected-or-not across that grid.
- The point isn't to relitigate Phase 1's choices — it's to know, once you have a real result, whether that result is stable or whether it flips under a small, defensible change to an arbitrary setting. A stable result across the grid is a materially stronger finding than one that only holds under the exact settings you happened to pick.

**Gate:** power curves exist for all three modes as used in your design; the sensitivity harness runs against synthetic data before you point it at anything real.

## Phase 6 — Run all three on your real data, write it up

- Run Standard/Sequential/Cumulative on whichever parts of your dataset match the Phase 0 event calendar.
- Run each real result through the Phase 5 sensitivity harness.
- Report, per result: plot, p-value, effect size, the power curve context (was this test even capable of detecting an effect size you'd care about?), and sensitivity across the config grid.
- State plainly: unusual-relative-to-history vs. proof of causation (no control site), and whether you checked an obvious confound.

**This deliverable answers your original questions, across all three modes, with the robustness checks to back them up.**

## Phase 7 — Still gated: build only if you also have concrete need

- **Audit mode / legal-defensibility docs** — only if Phase 0's scrutiny answer was "yes."
- **A full alternative-method benchmark suite** (vs. ARIMA, t-tests, etc.) — only if this is being packaged for other people to depend on.

## Phase 8 — MonitoringSeries (added post-Phase-6, on request)

Ongoing, open-ended yearly monitoring — check each new year as it
arrives, indefinitely, with no specific "treatment" and no fixed batch
size known in advance. None of Phases 0-6 cover this: Sequential's
Holm-Bonferroni correction needs the whole batch of p-values ranked
together, so it can't be revealed one year at a time; Cumulative is
anchored to one already-known event, and here there isn't one.

**Method:** fixed horizon, flat alpha split — commit to `n_years_horizon`
checks, split `alpha_total` evenly across them. Chosen over two
alternatives discussed and rejected for this "lite" pass: an open-ended
shrinking-budget schedule (valid forever, but gets progressively
underpowered) and a discovery-adaptive schedule in the LORD/SAFFRON
family (stays powered, but is materially more new statistical machinery
and controls a different guarantee — FDR, not "ever flag falsely").

**Built:** `tide_lite/monitoring.py` — `MonitoringSeries`, stateful and
incremental (unlike the batch-style Sequential/Cumulative), with
save/load since real usage spans separate sessions across years.
`tests/test_monitoring.py` — 5/5 passing.

**Finding, the same shape as the two in the section above:** this mode
needs meaningfully more historical years than the other three for its
budget to hold. Splitting an already-small `alpha_total` across many
checks pushes each one deep into the tail of the historical bootstrap
distribution, where bootstrap tail estimates are a known weak point at
small sample sizes. Measured: a nominal 5% lifetime budget's *true*
false-alarm rate was ~17% at 10 historical years, ~7% at 20, ~3% at 30.

**Also worth recording:** a self-check tool (let a user empirically
verify this against their own data) was attempted and rejected, not
just deferred. A leave-one-out design was built and initially looked
right, but its result pattern — calibration getting *worse* with *more*
historical years — matched nothing else found in this project. Tracing
it (not just re-deriving the theory) found a real structural issue:
repeatedly holding out one real year from a small pool, with
replacement, produces heavy repetition that understates the true risk,
worst exactly where the risk is highest. A tool with that flaw would
give false reassurance in precisely the cases most needing a warning, so
it was removed rather than shipped caveated. Full account in
`tide_lite/monitoring.py`'s module docstring.
