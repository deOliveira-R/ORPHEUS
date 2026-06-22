---
name: issue-247-legA-moment-source
description: Elegance review of #247 Leg A (moment-resolved external source for 2-D Cartesian LD) on refactor/sn-foundation-cleanup — PASS-WITH-NITS
metadata:
  type: project
---

# Elegance Review: #247 Leg A — moment-resolved external (slope-)SOURCE for 2-D Cartesian LD

Branch `refactor/sn-foundation-cleanup`, UNCOMMITTED. Scope: `orpheus/sn/solver.py`
(`_build_fixed_source_rhs` typed-union validation + `_lift_external_source_to_moments`
widening) and the `#247` block of `tests/sn/verification/mms/test_mms_ld_2d.py`
(lines ~486-1037).

## Summary Verdict
**PASS-WITH-NITS** — the carve is architecturally sound: a genuine typed-union
widening (flat OR moment-resolved ndarray rank), single-sourced lift, S4-safe
rank discrimination, no Pattern-6 production projection branch, exemplary
structural-not-converged-flux teeth. Three NITS, all do-now, all collapse-or-
clarify (none block the math). The headline adjudications the brief asked for:
the lift's re-check is JUSTIFIED self-defense (keep); the rank primitive
`is_moment_valued_by_rank` should be REUSED in the lift (NIT-1, the one real
single-source gap); the inline single-cell test projector is ACCEPTABLE.

## Findings

### PASS — Pattern 4 (illegal states): typed union discriminated by RANK not trailing-size
**Location**: `solver.py:1898-1922` (`_build_fixed_source_rhs`), `:1968-1986` (lift).
The union {flat `(N,ng,*spatial)`, moment-resolved `(N,ng,*spatial,2^d)`} is
discriminated by RANK, with the explicit rationale that a coincidental spatial
dim `== 2^d` could mis-fire a `shape[-1]` probe — this is the exact S4 hazard
`moment_layout.is_moment_valued_by_rank` was minted to dodge (#246), and the
diff cites it. The `n_cell_moments > 1` guard on `is_moment_resolved` makes a
moment-resolved input on DD/Step (`per_axis == 1`) structurally unrepresentable
(falls through to the flat-shape ValueError), which the test pins at `:1003-1010`.
Reads like the domain.

### PASS — Pattern 7 (normalise at definition site): no per-octant reframe duplication
**Location**: lift docstring `:1962-1965`.
The new moment-resolved arm threads the GLOBAL-frame slope rows through
unchanged and relies on the EXISTING per-octant `octant_moment_frame_signs` /
`sweep_graph._CellSolve` reframe to re-sign external slopes global→sweep
"EXACTLY as it does the scattering slopes — no new cell branch." This is the
correct Pattern-7 posture: the sign convention lives at the one sweep
definition site; the source path does NOT re-apply it. No new convention
habitat opened. Verified at the sweep call site (`sweep_graph.py:931`) that the
SAME `is_moment_valued_by_rank` gate re-signs moment-valued sources.

### PASS — Pattern 6 (defer abstraction): NO callable-projection branch in production
**Location**: `_build_fixed_source_rhs` / `_lift_external_source_to_moments` (whole).
Confirmed by read + grep: production accepts only a pre-projected ndarray; the
caller projects (`∫q·Pₖ` via leggauss). The test does its own projection
(`_project_scalar_to_tensor_legendre`, `:508`). No `Callable`/`leggauss`/
projection consumer was added to production — correct: there is no production
consumer yet, so a projection primitive would be premature (the lone
`moment_projection` hits at `solver.py:379+` are the unrelated Phase-5c
ANGULAR moment-OUTPUT path). PASS.

### PASS — adjudication: the lift's re-check is JUSTIFIED self-defense, KEEP (do not collapse)
**Location**: lift `:1980-1985` (`bulk_values.shape[-1] != n_cell_moments`) vs
`_build_fixed_source_rhs:1905-1914` (the negative pin).
The brief asks: redundant defense or justified? **Justified — keep both.**
Reasoning:
1. The two checks have DIFFERENT shapes and DIFFERENT contracts. The caller's
   negative pin (`:1905`) fires on right-rank/wrong-width and names the full
   moment vector for the public API error surface. The lift's check
   (`:1980`) is the helper guarding its OWN postcondition (it is about to
   thread the array through as a `2^d` moment vector) — it is the "parse,
   don't validate" boundary for the helper's contract, exactly the
   anti-pattern-8 inversion done RIGHT.
2. The docstring's "shared by the fixed-source and eigenvalue paths" claim is
   NOT a false single-source claim (the recurring stale-single-source tell):
   grep confirms ONE direct caller (`_build_fixed_source_rhs:1931`), but BOTH
   the fixed-source entry AND the eigenvalue entry (`solver.py:2157`) route
   THROUGH that one caller. The lift IS transitively shared; the docstring is
   accurate. Good.
3. The reachability gap is real but benign: today the only caller pre-checks,
   so the lift's `:1980` raise is unreachable from production. That is the
   CORRECT direction of redundancy for a helper that is a single source —
   collapsing it would make the helper unsafe to call from a future second
   caller (the very sharing the docstring promises). This is a true Pattern-2
   "single source defends its own contract," not a twin path. KEEP.

### NIT-1 (do-now, Pattern 2 / single source) — reuse `is_moment_valued_by_rank` for the rank test
**Location**: lift `solver.py:1973-1974`
(`flat_ndim = 2 + len(sn_mesh.spatial_shape)`; `if bulk_values.ndim == flat_ndim:`).
**Skill reference**: Pattern 2 (single source of truth) + anti-pattern #1 (two
spellings of one predicate).
**Problem**: "is this buffer moment-valued, judged by rank" is a project
PRIMITIVE — `orpheus.numerics.moment_layout.is_moment_valued_by_rank(array,
reference)`, `array.ndim > reference.ndim + 1`, single-sourced across the cell
(`sweep_graph.py:931`) and the matvec broadcast (`loss_representation.py:515`).
The lift RE-DERIVES the same predicate inline as a third spelling
(`bulk_values.ndim == flat_ndim`, i.e. the negation). The diff's OWN docstring
(`:1949`) even says it discriminates "via `is_moment_valued_by_rank` against the
per-ordinate-stripped flat rank" — but the code does NOT call it; it open-codes
the rank arithmetic. Docstring claims the primitive; code spells a twin. The two
spellings are provably coextensive TODAY (a `(ng,*spatial)` reference has
`ndim = 1 + len(spatial)`, so `NOT is_moment_valued_by_rank(bulk, ref)` ⟺
`bulk.ndim == 2 + len(spatial) == flat_ndim` — I verified the arithmetic), which
is exactly why this is a NIT not a VIOLATION.
**Bug-habitat argument**: the day the moment-layout convention shifts (a second
leading axis, a `2^{d-1}` face variant, the Leg B boundary trace that #251 adds
a DIFFERENT rank for), the primitive will be updated at its definition and at
the two registered consumers — and THIS open-coded third spelling will silently
diverge, exactly the divergence the primitive was minted to prevent (#246). A
source that mis-classifies flat-vs-moment lifts the slope rows onto the wrong
slot — a silent wrong-physics source, the worst class.
**Required change**: the brief raised the "is the reference array on hand?"
objection. It is the only real friction. The lift takes `(bulk_values,
sn_mesh)`; it has no cross-section array, only `sn_mesh`. Two clean options,
either acceptable:
  (a) **Preferred** — mint/keep the rank reference from the shape the lift
      already knows: `is_moment_valued_by_rank` needs only `reference.ndim`, and
      the lift already computes the flat rank (`flat_ndim`). Add a thin overload
      OR a sibling `is_moment_valued_by_flat_rank(array, flat_ndim: int) -> bool`
      in `moment_layout.py` (`array.ndim > flat_ndim`) and call it here AND have
      the existing `is_moment_valued_by_rank(array, reference)` delegate to it
      (`return is_moment_valued_by_flat_rank(array, reference.ndim + 1)`). One
      primitive, three consumers, zero open-coded rank arithmetic.
  (b) If (a) is judged scope-creep for Leg A: at MINIMUM make the inline test
      read as the primitive's negation with a cross-ref comment naming
      `is_moment_valued_by_rank` and a removal trigger (route through it when the
      shared overload lands) — anti-pattern #11. But (a) is the principled fix
      and is small; prefer it.

### NIT-2 (do-now, anti-pattern #11) — docstring names the primitive the code does not call
**Location**: lift docstring `solver.py:1949`
("discriminated by RANK ... via `is_moment_valued_by_rank` against the
per-ordinate-stripped flat rank").
**Skill reference**: anti-pattern #11 (a claim the code does not honor is a bug
habitat — a future maintainer trusts the docstring, greps for the call, finds
none, and is misled about the single source).
**Problem**: the docstring asserts the lift uses the shared primitive; the body
open-codes it (NIT-1). A comment that misstates an invariant is itself a habitat.
**Required change**: resolved automatically by applying NIT-1(a) (the code then
matches the docstring). If NIT-1 is deferred, soften the docstring to
"discriminated by RANK (the same predicate as `is_moment_valued_by_rank`, here
open-coded against the flat rank — TODO route through the shared overload, #247
follow-up)".

### NIT-3 (adjudication, ACCEPTABLE-as-is) — test projector: mesh-wide helper vs inline single-cell copy
**Location**: `test_mms_ld_2d.py:508` (`_project_scalar_to_tensor_legendre`,
mesh-wide) vs `:577-591` (inline single-cell projector in
`test_tensor_legendre_projection_matches_hand_polynomial`).
**Skill reference**: Pattern 2 / anti-pattern #1, weighed against L11 (structural
independence) and Pattern 6.
**Adjudication**: **ACCEPTABLE — keep the inline copy.** The two are NOT a twin
to collapse: the foundation test's inline projector is the L11 STRUCTURAL-
INDEPENDENCE reference (it pins the SHARED helper's normalization against
hand-laid polynomial algebra — `q̄ = a00 + a10·xc + ...`). If that foundation
test consumed the shared helper, it would be checking the helper against itself
(tautology) — the exact L11 trap the #247 spec is careful to avoid ("the
projector calls leggauss directly — NEVER `_lift_external_source_to_moments`",
`:501`). The inline copy is the independent oracle; the docstring justifies it
("one cell to keep the hand reference trivial", `:577-578`). This is the
keep-the-reference-impl exception done right (the shared helper has its own
consumers `:666/:691`; the inline is the pin). The ONE residual smell: the
single-cell projector is itself ~12 lines duplicating the helper's per-orientation
loop body. That is the price of independence and is bounded (two copies, one is
the oracle). No change required; recording the reasoning so a future de-dup sweep
does NOT collapse the oracle into the helper.

### PASS — the structural teeth read like the contract (master standard, test side)
**Location**: `test_ld_2d_external_slope_source_threaded_through_lift` (`:733`),
the M1-M4 mutation controls (`:869-969`), the negative pin (`:979`).
The teeth target the ACTUAL production change (the slope rows are no longer
zeroed) at O(1), not the sub-floor converged flux — `test_..._threaded_through_lift`
asserts `lifted == Qm` exactly AND the flat negative control (`flat_lift[...,1:]
== 0`) in the same test, so a regression to zeroing reddens at machine
precision. The `_CONSUMPTION_TOL = 1e-8` band is named with a rationale comment
(Pattern 3) explaining why it is ~5e7× above the inner tol yet below the spec-§0
value-band trap. The Mode-10 asymmetry (flat gate stays GREEN under the same
flip because it feeds an already-zero row) is pinned IN the mutation test
(`:902-915`) — exemplary. Leg B correctly SKIP-stubbed and routed to #251
(`:1022`), a tracked deferral (anti-pattern #11 EXCEPTION, done right).

## Architectural Opportunities
The `is_moment_valued_by_rank` primitive currently demands a `reference`
ARRAY purely to read `reference.ndim + 1`; three consumers
(`sweep_graph`, `loss_representation`, and now the open-coded lift) all really
want "is this rank above the flat rank." A `is_moment_valued_by_flat_rank(array,
flat_ndim: int)` core with the array-taking form delegating to it would let the
SOURCE path (which has no cross-section in scope, only a shape) join the single
source cleanly. Small, in-scope-adjacent; fold into NIT-1(a) rather than filing
a separate issue unless the main agent prefers to defer.

## Approval Conditions
1. **NIT-1** — route the lift's flat/moment rank test through the shared
   `is_moment_valued_by_rank` primitive (preferred: add the `flat_ndim` overload
   in `moment_layout.py` and have both forms delegate; the source path then
   joins the single source). If deferred, apply NIT-1(b) (cross-ref + removal
   trigger comment) at minimum.
2. **NIT-2** — docstring `:1949` must match the code (auto-resolved by NIT-1(a);
   else soften the "via `is_moment_valued_by_rank`" claim to "open-coded, TODO
   route through").
3. NIT-3 requires NO change (recorded as a do-not-collapse note for future
   de-dup sweeps).

The lift re-check (KEEP), the rank discrimination (PASS), the no-projection-in-
production (PASS), and the teeth (PASS) are all sound. With NIT-1/NIT-2 applied
the source path joins the moment-layout single source and the carve is clean to
commit.
