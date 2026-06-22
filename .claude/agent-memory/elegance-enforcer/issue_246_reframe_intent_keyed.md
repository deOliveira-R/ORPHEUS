---
name: issue-246-reframe-intent-keyed
description: #246 PASS-WITH-NITS — _reframe keyed on is_moment_valued (typed origin) not trailing-length coincidence; is_multi_moment named predicate; the has_spatial_moment_axis Pattern-6 TRIM + the rank-discriminator twin.
metadata:
  type: project
---

# #246 — `_reframe` keys the moment-frame involution on INTENT (branch `refactor/sn-foundation-cleanup`, uncommitted at review)

PASS-WITH-NITS (2 follow-up, 0 blocking). The core fix is a genuine elegance GAIN: replaces the
coincidental `arr.shape[-1] != frame_signs.shape[0]` guard (S4 hazard: a d=2 non-moment array of
trailing length `2^d==4` mis-fires the sign involution) with `frame_signs is None or not is_moment_valued`
— the caller keys moment-ness on the array's TYPED ORIGIN. New named predicates `is_multi_moment`
(scheme) / `has_spatial_moment_axis` (field).

**Why:** counterweight review BEFORE commit; user resolves nits.
**How to apply:** if this carve re-surfaces or the nits get deferred to an issue, the rulings below stand.

## The 5 scrutiny rulings (verified, not asserted)
1. ⭐ **`has_spatial_moment_axis` (`transport/fields/_bases.py:255`) → TRIM (Pattern-6).** ZERO production
   consumers (grepped whole tree) — only the 3 Gate-3 tests + the implementer's OWN docstring says the
   inner walk MUST use the scheme predicate not this field query. Abstraction whose own docs explain why
   prod will never use it = textbook over-abstraction. The P4' invariant it was minted for
   (`bare-LD-field → no moment factor`, the construct-general regression guard) does NOT need it:
   `spatial_moments_per_axis` is a PRE-EXISTING property (`_bases.py:231`, untouched by #246) and the
   invariant IS `spatial_moments_per_axis == 1`. Retarget the 3 tests onto it. KEEP only if a real
   production method-boundary consumer is added THIS commit (Pattern-6 ≥2 instances).
2. **Source-site rank test `Q_cells.ndim > sigt_cells.ndim + 1` (`sweep_graph.py:929`) = SECOND SPELLING
   of `_moment_broadcast_sigma`'s `moment_valued.ndim > sig.ndim + 1` (`loss_representation.py:515`).**
   Reads cleanly inline BUT it is the 3rd member of the rank-discriminator family. ACCEPTABLE-FOR-NOW
   (non-divergent, same convention) → my Pattern-1/Pattern-2 twin tripwire. Do-now remedy: extract one
   named `_is_moment_valued_by_rank(arr, ref)` in the convention-owning module (loss_representation,
   beside `_moment_broadcast_sigma`); minimum = reciprocal cross-ref comments + tracked removal trigger
   (the 3rd DIVERGENT edit). Latent: a future moment source carrying 2 extra axes lands on one spelling.
3. ⭐ **`is_multi_moment` on Base-NOT-Protocol = ARCHITECTURALLY CORRECT (clean PASS, not a judgment call).**
   The codebase EXPLICITLY documents (`scheme.py:553-559`) the Protocol declares ONLY the universal
   `update`/`residual` contract so a minimal synthetic mock conforms on 2 members; scan-family
   capabilities (`cell_kernel_batch`/`residual_kernel_batch` — exactly what `is_multi_moment` gates) are
   Base-only opt-in. Putting it on the Protocol would be THE VIOLATION (breaks minimal-mock conformance).
   Implementer's stated reason matches the codebase's own rationale verbatim.
4. **Per-callsite `is_moment_valued` mapping reads HONESTLY (non-uniform, correct).** `_CellResidual` Q
   static `False` verified: `_MATVEC_ZERO_SOURCE = np.zeros((1,1,1))` always rank-3 flat. `_CellSolve` Q
   rank test (dual entry: prod lifts to moment / low-level `sweep` API flat). All other sites
   `is_multi_moment`. Single-sourced — `_OneDimScanWalk._sweep_direction` hoists once
   (`loss_representation.py:2297`); no call site re-spells `spatial_basis_per_axis > 1` as a live boolean.
5. **Docstrings EXEMPLARY (Cardinal Rule 3).** `_reframe` includes a `#246:` paragraph documenting WHY
   the old shape-coincidence guard was wrong (the d=2 `2^d==4` mis-fire) — design rationale not behaviour.

## Equivalence VERIFIED (the docstring's tri-spelling claim is provably exact)
`is_multi_moment ⟺ frame_signs is not None ⟺ moment_tail != ()`:
- `octant_moment_frame_signs` returns `None` iff `per_axis == 1` (`_ubld.py:146`).
- `is_moment = face_moment_tail(per_axis**ndim) != ()` (`loss_representation.py:2819`) reduces to
  `per_axis**ndim != 1` = `per_axis != 1` = `per_axis > 1` = `is_multi_moment` for ALL ndim≥1, per_axis≥1.
- So `is_moment` (the layout-derived branch flag at `:2819`) IS the 3rd spelling, but LOAD-BEARING in that
  scope (`moment_tail` also builds buffer shapes `:2830/:2842/:2847`) → NOT a collapse target, keep.

## STANDING TELLS (promote-worthy if they recur)
- **A named typed-query predicate minted alongside a test pin, with ZERO production consumer and a
  docstring that says prod must NOT use it** = Pattern-6 TRIM on sight. The test pin almost always already
  has a pre-existing lower-level property to assert against — grep for it before keeping the new predicate.
- **`X.ndim > ref.ndim + 1` ("moment-valued by rank") is becoming a recurring family** across
  loss_representation + sweep_graph. Co-locate the eventual primitive with `face_moment_tail` /
  `octant_moment_frame_signs` in the moment-layout leaf (`numerics/moment_layout.py` / `_ubld.py`) when a
  3rd DIVERGENT consumer lands — same home as the other moment-layout single-sources (#245).
