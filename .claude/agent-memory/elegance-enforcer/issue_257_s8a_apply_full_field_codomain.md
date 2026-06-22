---
name: issue-257-s8a-apply-full-field-codomain
description: S8a apply(FullField)->FullField re-point (cofree-comonad operator-output timeless) — PASS; the role-grid discriminator is iterate-builder-vs-base-arrow NOT apply-vs-solve; M[f].solve flips timeless correctly
metadata:
  type: project
---

# #257 S8a — apply(FullField) re-point + helper widening (the cofree-comonad operator-output split)

PASS (`feature/field-typed-operator-algebra`, baseline `90ea1cd`, UNCOMMITTED, behavioral-NEUTRAL).
Re-points every SN operator matvec LEAF (`StreamingOperator`/`InvertibleOperator`/`MultiplicationOperator`/
`SNBoundaryOperator`/`FissionOperator`+`ScatteringOperator` TimedFullField OUTPUT arms) from emitting the
history-bearing `TimedFullField` to the TIMELESS `FullField` — dropping `_history=()`/`history_depth=psi.history_depth`
stamping. `TimedFullField = Cofree(FullField, d)`; an operator output is a base arrow `FullField -> FullField`,
the comonad lives ONLY on the iteration DRIVER. The S8a-immediate consumer of the S4.5 `_recombine` infra.

## The load-bearing ruling: the role-grid discriminator is ITERATE-BUILDER vs BASE-ARROW, not apply-vs-solve

The prior S3b/S4.5 "operator-output vs solve-trace" discriminator REFINES here. `InvertibleOperator.solve`
(the WDD sweep) STAYS `TimedFullField -> TimedFullField` — VERIFIED right: it is the driver-facing resolvent
that BUILDS the iterate (output type matches input; `_solve_timed_full_field` body `operator.py:1154` keeps the
intentional `TimedFullField(_history=(), history_depth=rhs.history_depth)` — NOT a half-done seam). But
`MultiplicationOperator.solve` (`q/f -> AngularFlux`, the S3b engine inverse) FLIPS to `FullField -> FullField`
and that is ALSO right — VERIFIED ZERO graph callers (`mcp__nexus__callers` empty; only `engine.solve` internal):
it is an ALGEBRAIC-INVERSE base-arrow leaf, never the iteration resolvent. So `solve` is NOT uniformly timed;
the true axis is "does the iteration driver consume this as the iterate it carries forward?" (InvertibleOperator.solve = YES → timed)
vs "is this a composed algebraic action?" (M[f].solve = NO → timeless). STANDING TELL: when judging a `solve`
re-type, run `callers` — a zero-caller `.solve` is an algebraic leaf (flip timeless), a driver-facing resolvent
stays timed. Do NOT flag the M[f].solve/InvertibleOperator.solve asymmetry as inconsistency.

## Driver re-attach (design choice #1) — ORDER-DEPENDENT but the fold is order-correct; verify the fold

No explicit `advance` call: the SI step is `rhs = q_ext; for g in gains: rhs = rhs + g.apply(psi)`
(`iteration.py:525-527`). `q_ext` is TIMED (LEFT operand) + timeless gain (RIGHT) → `TimedFullField.__add__`
runs `self._check_partner(other)` (passes: timeless IS-A FullField) then `self._recombine(...)` (the TIMED
override `timed_full_field.py:252`) → timed result. Timed stays timed at every fold step BECAUSE q_ext is the
left seed. ⚠ This is order-dependent: `gain + q_ext` (timeless LEFT) would `_recombine` to TIMELESS and silently
drop history-depth. The SI fold is timed-first → safe. Krylov leg (`loss_minus_gains` `iteration.py:765`) ravels
to ndarray → timed/timeless label irrelevant, final rebuild via `solution_template._recombine`. So "let `__add__`/
`_recombine` reattach the comonad" IS the elegant move (NOT a latent gap) — the absence of `advance` is correct
because steady-state never advances (history_depth carried, `_history=()` always). STANDING TELL: when a re-point
relies on a dunder fold to reattach a type, the elegance verdict HINGES on operand order in the driver — read the
actual fold (`rhs = rhs + g`) and confirm the type-bearing operand is the accumulator/left seed.

## 4 design choices ALL PASS

- **#1 driver re-attach via `__add__`/`_recombine`, no `advance`**: PASS (see above; single-sourced, order-correct).
- **#2 `_require_typed_composite` widened `TimedFullField -> FullField`** (`operator.py:170`): PASS, textbook
  Pattern-4 — ONE guard `isinstance(field, FullField)` accepts timed+timeless (TimedFullField IS-A FullField),
  rejects bare ndarray; SI/Krylov drivers passing a timed iterate still satisfy it. SSOT of the input contract.
- **#3 InvertibleOperator.solve stays timed / apply goes timeless**: PASS, principled + documented (the resolvent
  builds the iterate covariantly under the `["FullField"]` class binding; the matvec leaves are base arrows).
- **#4 `full_field_space._rebuild` routes `replace(x,...,_history=())` -> `x._recombine(...)`** (`:197`): PASS —
  this is the METRIC analogue of the dunder recombine (unifies on the SAME `_recombine` hook, NOT a parallel path).
  The helper was ALREADY duck-typed on untyped `x`; routing through `._recombine` adds ZERO pyright errors and
  now correctly handles both carriers (timeless metric input post-S8a). Single-source win.

## Scope discipline CONFIRMED (brief axis #4)

`(L+C)−C` arithmetic in StreamingOperator.apply/apply_transpose UNTOUCHED (`out_bulk = lpc.bulk.values -
self.sigma_t[None]*psi.bulk.values`). Scattering/Fission `@apply.register def _(self, psi: TimedFullField)` INPUT
dispatch UNTOUCHED — ONLY the `-> "FullField"` output annotation + the `FullField(...)` construction changed
(diff confirms). No drift into S8b (`(L+C)−C` drop) or S8c (Scattering/Fission fibration).

## Consistency — every leaf re-pointed, no half-done seam

Grepped `_history=()`/`history_depth=psi|phi|q` across all 6 operator files: the ONLY survivor is
`operator.py:1154` inside `_solve_timed_full_field` (the InvertibleOperator.solve body = intentionally timed).
Test migrations STRENGTHEN: every old `assert isinstance(out, TimedFullField)` -> `isinstance(out, FullField)`
+ `not isinstance(out, TimedFullField)` (the old assert would PASS on a timed output; the new explicitly forbids).
Retirement = test rewiring, not deletion.

## C5 test (brief axis #5) — genuine intrinsic-property gate

`tests/sn/operators/test_apply_full_field_codomain.py` (foundation, 14✓+1xfail). NOT shallow: C5a probes the
DEFINING property (`type(out) is FullField and not isinstance(.,TimedFullField)`) across slab/sphere/cylinder
+ across input history_depth 0..4 (the timeless codomain is depth-independent) + Mode-11 aware (calls `.apply`
DIRECTLY because the matvec leaf has zero graph callers — a solve-only test routes around it). C5b pins the
driver-reattach end-to-end against the STRUCTURALLY-INDEPENDENT closed-form `k_inf=νΣf/Σa` (NOT old-vs-new ULP,
NOT a baseline L28 would force-mutate) for SI+Krylov × geometry, plus a direct `SourceIteration` build proving
the iterate stays timed. C5c keeps the `advance` type-guard. `-O`-firing throughout (pytest.fail/raises, no bare assert).

## Gates

C5 14✓/1xfail; sn/solve + transport 327✓; operators dir 7 reds ALL PRE-EXISTING baseline (3 SPHERE
`test_sphere_*_apply_bit_identical` = stale snapshots #232/#250, 100%-mismatch ~3.4 rel NOT a neutrality break —
the C5b sphere k_inf gate PASSES, structural proof the sphere apply path is value-correct; `ymin requires genuine
mu_y` in test_native_matvec/test_boundary_conditions = #214 ny=1 cluster; failing test bodies NOT in the S8a diff
hunks, test_boundary_conditions untouched). pyright per-file counts entirely pre-existing: scattering/fission
`register`/`NoReturn` singledispatchmethod limit ([[issue-256-apply-endomorphism-step3]]) + full_field_space
untyped-`x` metric-method `reportOptionalMemberAccess` at :230-255 (the `_rebuild` :197 change adds 0). 
multiplication_operator pyright-CLEAN.
