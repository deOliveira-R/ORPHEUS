---
name: issue-257-s7-scipy-single-source-closeout
description: "#257 S7 — consolidate the two inline scipy LinearOperator closures in KrylovAcceleration.solve into ONE ravel-aware adapter `_as_scipy_linop`, and retire the dead public `as_scipy_linop` (flat-only twin). Bit-identical plumbing refactor; 0 net-new pyright; net -1 type:ignore."
metadata:
  type: project
---

# #257 S7 — scipy single-source (bit-identical plumbing refactor)

Branch `feature/field-typed-operator-algebra`, HEAD `9ccfec7`, NOT
committed (main agent commits after review). Host env, `.venv/bin/python`.

**Why:** Cardinal Rule 2 (architecture / single source of truth). The
live ORPHEUS↔scipy Krylov boundary was TWO inline closures (`A_scipy`,
`M_scipy`) in `KrylovAcceleration.solve` (`iteration.py`), both wrapping
a carrier→carrier function with the ravel pair `_ravel`/`_unravel_like`
— twin paths sharing the ravel-wrapping structure. Separately,
`as_scipy_linop` (`operator.py`, public + exported + documented) was a
flat-only adapter with **ZERO production callers** (audited pre-task via
Nexus + grep; only 5 tests + 3 doc refs). The flat case is the
degenerate case of the ravel-aware adapter, so the live adapter
subsumes it. Layering (`iteration.py → operator.py`, one-way) forced the
single adapter to live in `iteration.py`.

## Deliverables (file:line summary)

1. **`orpheus/numerics/iteration.py`** — ADD module-level private
   `_as_scipy_linop(carrier_matvec, template, n: int) -> spla.LinearOperator`
   right after `_unravel_like` (now at line 210, before `_zeros_like`).
   Body: an inner `_flat_matvec(flat)` = `_ravel(carrier_matvec(_unravel_like(template, flat)))`,
   returns `spla.LinearOperator((n, n), matvec=_flat_matvec, dtype=float)`.
   `carrier_matvec` + `template` UNANNOTATED (mirrors `_unravel_like`'s
   loose `template` — accepts typed carriers AND bare ndarrays;
   annotating invites a pyright union headache). This is the SOLE site
   that constructs a scipy LinearOperator for the Krylov accelerator.

2. **`orpheus/numerics/iteration.py`** — `KrylovAcceleration.solve`
   (~line 770): REPLACE the two inline closures with a NAMED carrier-space
   matvec `loss_minus_gains(psi: V) -> V` (the honest within-group system
   matvec, `out = self.L.apply(psi); for g in self.gains: out = out - g.apply(psi)`
   — SAME calls, SAME left-to-right order, bit-identical) + two
   `_as_scipy_linop(...)` calls:
   - `A_scipy = _as_scipy_linop(loss_minus_gains, solution_template, n)`
   - `M_scipy = _as_scipy_linop(self._preconditioner, q_ext, n) if ... else None`
   The TYPED `loss_minus_gains(psi: V) -> V` form (V is in scope from
   `KrylovAcceleration(Generic[V])`) produced ZERO pyright friction — kept
   it (reads like the domain).

3. **`orpheus/numerics/operator.py`** — DELETE `as_scipy_linop` (the
   `# scipy interop` banner + the whole def, was lines 1780-1836, end of
   file), incl. its `# type: ignore[attr-defined]` (net -1 ignore).
   Remove `"as_scipy_linop"` from `__all__`. Remove the dedicated
   module-docstring paragraph (was lines 46-50) that named it. **REMOVE
   the now-orphaned `import scipy.sparse.linalg as spla`** (line 67) —
   audited: `spla` was used ONLY by `as_scipy_linop` (lines 1789 + 1835),
   zero other refs. All other deps (`_has`, `CAP_APPLY`,
   `CAP_APPLY_TRANSPOSE`, `MissingCapability`) heavily used elsewhere, kept.

4. **`orpheus/numerics/__init__.py`** — remove `as_scipy_linop,` from the
   `from .operator import (...)` block (was line 31) + `"as_scipy_linop",`
   from package `__all__` (was line 87).

5. **`tests/numerics/test_operator.py`** — remove `as_scipy_linop` from
   the import; DELETE the 5 orphaned tests + their `# scipy interop`
   banner (was lines 415-460): `test_as_scipy_linop_matvec_matches_apply`,
   `_rmatvec_when_transpose_capable`, `_no_rmatvec_when_not_transpose_capable`,
   `_rejects_missing_apply`, `_works_with_composite`. Fixtures
   (`matrix_full`, `matrix_apply_only`, `vector`, `NoApplyOperator`) NOT
   orphaned (heavily used by other tests; audited via grep) — kept.

6. **`docs/theory/operator_algebra.rst`** — remove the Key-Facts bullet
   dedicated to `as_scipy_linop` (was ~226); reword the prose paragraph
   (was ~391) to state the scipy boundary is now internal to
   `KrylovAcceleration` (single ravel-aware adapter, single source of
   truth) instead of pointing at the retired adapter.

7. **`docs/theory/discrete_ordinates.rst`** — repoint the illustrative
   `inverter = lambda q: gmres(as_scipy_linop(L), ...)` example (was
   ~11685) to `KrylovAcceleration(L, ...).solve(q)[0]` + a sentence that
   the ORPHEUS↔scipy boundary is internal to `KrylovAcceleration`.

## Gates (all green)

- **pyright**: baseline **2307 errors / 19 warnings** → after **2297
  errors / 19 warnings** (−10 errors, 0 net-new). The −10 (more than the
  brief's predicted "exactly 2307") is the deleted `as_scipy_linop` body
  + 2 of the deleted tests' pre-existing `#226`-noise errors disappearing
  — NOT masking (focused per-file pyright on the 4 touched files shows 22
  errors, all pre-existing `#226` protocol/operator-overload noise on test
  fixtures, none introduced). `# type: ignore` in operator.py **14 → 13**
  (the deleted `[attr-defined]`). ZERO net-new ignores.
- **Krylov-path gates** (138 tests): `test_fixed_source_g1`,
  `test_fixed_source_2d_equivalence`,
  **`test_krylov_curvilinear_precond_safety`** (the M_scipy
  preconditioner branch via the `q_ext` template — MUST pass, did),
  `test_krylov_restart_signature`, `test_affine_carve_bit_identity`,
  `test_keff_2d`, **`test_timed_full_field`** (the `from_flat(to_flat(x))==x`
  round-trip the adapter relies on, incl. the white-box helper-name
  import test), `test_operator` → **138 passed** (-O).
- **Broad regression subset** (-O, deselecting `test_heterogeneous_absolute_keff`):
  **7 failed / 2018 passed / 7 skipped / 5 xfailed**. The 7 reds are
  EXACTLY the documented baseline: #250 SPHERE ×5 (`test_vacuum_bulk_bit_identical_1d[0/1/2-SPH]`
  + `test_sphere_1g/2g_apply_bit_identical`) + #232 mu_y ×2
  (`test_2d_mesh_resolution` + `test_two_d_cartesian_loss_action_returns_result`,
  both "Face 'ymin' requires genuine mu_y cosines"). ZERO non-baseline reds.
- **Sphinx -W**: build succeeded (only pre-existing SyntaxWarnings in 2
  test files via the graph builder + pre-existing `ld-cartesian-1d`/`ld-slab`
  Nexus equation-node skip notices — neither is a fatal rST warning; no
  dangling `:func:` xref, which `-W` would have caught).

## Helper-name preservation (constraint)

`_is_ravellable`, `_ravel`, `_unravel_like`, `_zeros_like`, `_l2_norm`
ALL survived verbatim (white-box test `test_timed_full_field.py:515-524`
imports them by name; passed). New `_as_scipy_linop` slots between
`_unravel_like` and `_zeros_like`.

## Posture / no manifest owed

PURE plumbing single-source refactor (no new reference solver, no new
verifiable claim, no equation minted). NO algebra-of-record Branch-1/L1/
Sphinx-stub manifest owed (same posture as a numerics-primitive
consolidation). NO new ERR (no bug found — the retired symbol was dead,
not wrong). Bit-identity inherited by construction: `loss_minus_gains`
makes the SAME `L.apply` + per-gain `out - g.apply` calls in the SAME
order; the ravel-wrapping moved verbatim into `_as_scipy_linop`'s
`_flat_matvec`. NOT committed.

## Surprises

- pyright dropped −10 (not the brief's predicted "exactly 2307"). Genuine
  (deleted code carried pre-existing errors), not masking — verified by
  focused per-file pyright. Cleaner than predicted.
- `spla` import in operator.py was orphaned by the deletion (brief flagged
  to check) — removed. This is the one edit beyond the brief's explicit
  symbol list; it's the brief's STEP 3 "remove an import if it becomes
  genuinely unused" branch, executed.
- The typed `loss_minus_gains(psi: V) -> V` form (brief said "try typed
  first, fall back if friction") produced zero friction — kept typed.

## Archivist dispatch

A `DISPATCH_REQUEST` to archivist (`followup: false`) was emitted for a
consolidated rich-narrative pass on `docs/theory/operator_algebra.rst`
§"scipy interop" — the doc edits here were minimal mechanical
reference-cleanup tied to the retirement, NOT new narrative.

LESSON: when consolidating N inline scipy-LinearOperator closures into
one adapter while retiring a flat-only public twin, the flat twin's
`spla` import is almost certainly orphaned (the twin is the only
scipy-construction site in that module once the live boundary moves to
the caller) — grep `spla.` after the delete; the carrier-space matvec
extracted as a named `(psi: V) -> V` closure reads like the domain AND
takes the project's `Generic[V]` typevar without pyright friction (no
fallback to unannotated needed).
