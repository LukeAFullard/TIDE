# METHODS.md -- TIDE project (fill in before writing code)

Phase 0 and Phase 1 of the rebuild plan produce this document. Fill in
every `___`, delete this instruction line, and treat it as the source of
truth: every setting in your `TideConfig` should trace back to a line here.

## Phase 0 -- the question(s)

**Site(s) and variable(s):** ___

**Event calendar** -- every treatment event, and which question it needs:

| Event | Date/window | Question | Mode |
|---|---|---|---|
| ___ | ___ | is this one different? | Standard |
| ___ | ___ | is this one of several independently real effects? | Sequential |
| ___ | ___ | how is this one evolving over time? | Cumulative |

**Decision this informs:** ___

**External scrutiny?** (research/monitoring only, or must survive legal/
regulatory review) -- ___. This single answer decides whether Phase 7
(audit mode, legal-defensibility docs) is in scope at all.

## Phase 1 -- statistical design

Core (applies to every mode):

- **Time resolution:** ___ (weekly / monthly / other) -- why: ___
- **Historical period:** years ___ through ___ (___ years total)
- **Exchangeability check:** [ ] plotted all historical years overlaid,
  they read as one process. If not, note the concern here: ___
- **Detrending:** ___ (additive / log_additive / none — note:
  `log_additive` is also this tool's normalization for skewed data, not
  a separate setting; see USER_GUIDE.md's note on normalization if
  you're expecting a baseline/reference-level rescale, which this
  project deliberately doesn't do). If your variable
  is a concentration (nutrients, turbidity, etc.), did you check whether
  variance scales with level? ___
- **Baseline (V0) window:** ___ -- confirmed not to overlap any treatment
  period: [ ]
- **Test statistic:** max standardized deviation across bins (the
  `tide_lite` default) -- note here if you're using something else: ___
- **Effect size definition:** signed median deviation from the historical
  median curve, in original units (the `tide_lite` default) -- note here
  if you're using something else: ___
- **Envelope mode:** prediction (single treatment year/period) /
  confidence (averaging several years) -- ___

Sequential-specific (skip if not using this mode):

- **Correction:** Holm-Bonferroni step-down (the `tide_lite` default) --
  confirm alpha: ___
- **Which events are corrected together:** ___
- **Non-overlap rule:** how are irregularly-spaced event windows sliced so
  one event's data can't leak into another's test window? ___

Cumulative-specific (skip if not using this mode):

- **Which event** is being tracked: ___
- **Baseline-reuse confirmed:** one `fit_historical()` call, reused
  unchanged across every tracking window, each measured from the same
  original event: [ ]

## Sign-off

- [ ] Every implementation decision below traces back to a line above.
- [ ] No new `.md` files exist yet beyond this one and the README.
