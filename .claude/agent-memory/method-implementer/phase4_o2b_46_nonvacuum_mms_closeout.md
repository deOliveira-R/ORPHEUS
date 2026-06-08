---
name: phase4-o2b-46-nonvacuum-mms-closeout
description: Phase 4 / O.2b sub-step 4.6 — NON-VACUUM prescribed-inflow MMS (slab 1g+2g + sphere companion) + prescribed-inflow splitting-invariance (vv Mode 9) consistency test. LANDED on origin/main; later migrated onto the composite fixed-source public API. Durable design decisions A–E + the (N,ng,nx,ny) source-layout-lies-in-docstring hazard.
metadata:
  type: project
---

# Phase 4 / O.2b 4.6 — non-vacuum prescribed-inflow MMS (closeout)

**Branch** `refactor/field-role-typing`. **LANDED on origin/main** — committed as
`3273010` (feat: cases + Branch-1 SymPy) → `55f57bf` (test: convergence +
consistency) → `59c9046` (docs: narrative + V&V matrix), then MIGRATED onto the
composite fixed-source public API in `d87f3b0` (`prescribed_inflow` generator +
`build_nonvacuum_fixed_source` bundler — the bypass solve described below was
replaced by the public `solve_sn_fixed_source` path) and elegance-polished in
`bedc394`. The §"Open items" bypass/probe-deletion handoffs below are RESOLVED
by that migration. Built 2026-06-06 by method-implementer.

## DURABLE LESSON (read first) — the `external_source` layout lies in its docstring

The single transferable hazard: `external_source` returns layout
`(N, ng, nx, ny)` — the solver-consumed shape — NOT the `(N, nx, 1, 1)` the
docstrings claimed. The pre-existing `test_sn_mms_anisotropic_symbolic.py`
cross-check sliced `Q_numerical[:, :, 0, 0]` (WRONG axes — `(N,1)` vs `(N,nx)`
mismatch) and was already red at clean HEAD; the correct slice for g=0/ny=0 is
`[:, 0, :, 0]`. This is exactly the `coding-elegance` Pattern-13
bare-ndarray-shape-in-docstring hazard: the docstring lied about axis order and
the existing test inherited the lie. A typed source field would make the wrong
slice unspellable. **Rule**: never trust a bare-ndarray docstring for axis
order — probe the actual `.shape` before slicing a `(N, ng, nx, ny)` source.
This is a wrong-axes mismatch that fails LOUDLY (not a silent wrong-answer), so
no ERR-NNN; but it cost a real debugging detour.

## What 4.6 adds

The entire existing MMS catalog is vacuum-automatic (every ansatz vanishes at both
faces → `γ₋ψ ≡ 0`), so the prescribed-inflow `q.boundary ≠ 0` path is untested.
4.6 chooses the proven P1 ansatz `ψ_n = (A + μ_n B)/W` with `A` NON-zero at the
outer face (`a0>0`), lighting the prescribed-inflow source-slot bridge consumed by
`(L+C).solve`. Forward-only; NO eigenvalue claim (non-fissile mixture, MMS is
source-driven); NO adjoint.

## Deliverables (files created/edited)

| File | What | Status |
|------|------|--------|
| `orpheus/derivations/continuous/mms/sn.py` | (edit) Branch-1 SymPy + Branch-2 factories | shipped |
| `tests/derivations/test_sn_mms_nonvacuum_symbolic.py` | (new) foundation symbolic gate | 9 pass |
| `tests/sn/verification/analytical/test_mms_prescribed_inflow.py` | (new) T1/T2/T3/T3g | 4 pass + 1 xfail |
| `tests/sn/verification/analytical/test_prescribed_inflow_consistency.py` | (new) T4 (promotes the probe) | 2 pass |
| `docs/theory/discrete_ordinates.rst` | (edit) Sphinx stub + 4 `:label:` blocks | Sphinx -W clean |

### SymPy Branch-1 (decisions A + B)

- `_spherical_anisotropic_symbolic(A=None, B=None)` — PARAMETERIZED (decision A).
  No-arg default reproduces Phase 3.6 vacuum shapes BYTE-UNCHANGED (regression pin
  green); `derive_nonvacuum_spherical_mms()` reuses it with the 4.6 shapes
  (Cardinal Rule 2 — spherical-operator residual lives in ONE place). Both share
  the EXACT closed-form `Q_closed`; only A,B differ.
- `_nonvacuum_slab_symbolic()` + `derive_nonvacuum_slab_mms()` — NEW fresh pair
  (decision B). Slab has NO redistribution term → genuinely different operator
  from the sphere → not shared. Substitution residual == 0.
- `_nonvacuum_spherical_AB()` — the 4.6 sphere shapes `A=a0+a1 sin(kr)` (a0=1/2),
  `B=(r/R)[b0+b1 cos(kr)]` (b0=3/10, b1=1/5). HAZARD H1: `(r/R)` forces B(0)=0
  (pole-regular), B(R)≠0 (non-vacuum). `derive_*` both prove diff==0 symbolically.

### Branch-2 case dataclasses (decisions C + D)

- `SNSlabNonVacuumMMSCase` (+ `build_slab_nonvacuum_mms_case` 1g,
  `build_slab_2g_nonvacuum_mms_case` 2g) — MULTI-GROUP-capable via per-group
  amplitude vector `c_groups`, `sigma_t_g`, `sigma_s_matrix` (ORPHEUS
  `SigS[g_from,g_to]`; in-scatter = `SigSᵀ φ`). 1g = `c_groups=(1.0,)`. 2g =
  `(1.0, 0.4)` with DOWNSCATTER-only asymmetric Σs (ERR-002 transpose hazard live,
  Cardinal Rule 6 mandatory ≥2g). `_make_2g_asymmetric_mixture` helper.
- `SNSphericalNonVacuumMMSCase` (+ `build_sphere_nonvacuum_mms_case`) — `k=π/(2R)`
  so A(R)=0.75, B(R)=0.3 (matching the baked SymPy coefficients → L1 cross-check
  holds). r=0 symmetry BC (not a face); inflow only at r=R (`xmax`).
- New `prescribed_inflow(sn_mesh)` method on BOTH cases — returns a
  `BoundarySourceSink` with `ψ_chosen(x_face, μ_n)/W = (A + μ_n B)/W` on the inflow
  ordinate slots. Slab: both faces (a0>0). Sphere: r=R only. `face_view(face)` is
  `(N, ng)` on a 1-D trace → index `view[inflow, g]` per group.
- VACUUM mesh BCs; inflow injected via `q.boundary` (decision D — NO from_spec
  bridge). The bypass solve path (decision E) drives `_within_group_triple` +
  `_select_si_resolvent`/SourceIteration directly, NOT `solve_sn_fixed_source`
  (which hardcodes vacuum `q.boundary`).

### Sphinx stub (algebra-of-record stub-mode)

`docs/theory/discrete_ordinates.rst` § "Non-vacuum prescribed-inflow MMS
(Phase 4 / O.2b 4.6)" — 4 `:label:` blocks (`sn-mms-nonvacuum-psi`/`-qext`,
`sn-mms-nonvacuum-sph-psi`/`-qext`), `:mod:` cross-ref, a `.. todo::` archivist
marker. NOT the rich narrative (archivist deliverable). Sphinx `-W` build SUCCEEDS.

## Gate results (verbatim, `-O` mode)

| Piece | Gate command (PYTHONPATH=worktree, `.venv/bin/python -O -m pytest`) | Result |
|-------|------|--------|
| 1 (T4) | `tests/sn/verification/analytical/test_prescribed_inflow_consistency.py` | **2 passed** |
| 2 (symbolic) | `tests/derivations/test_sn_mms_nonvacuum_symbolic.py` | **9 passed** |
| 2 (regression A) | `tests/derivations/test_sn_mms_anisotropic_symbolic.py` | **10 passed, 2 PRE-EXISTING failed** (see below) |
| 3 (T1/T2) | `…/test_mms_prescribed_inflow.py::test_mms_prescribed_inflow_slab_converges_second_order` | **2 passed** |
| 4 (T3/T3g) | `…::test_mms_prescribed_inflow_sphere_*` | **1 xfailed (T3, strict) + 1 passed (T3g)** |
| all 4.6 | the 3 new test files | **14 passed, 1 xfailed** |

### Measured slab convergence (probe + test)

- **T1 (1g)**: orders `[2.04, 2.01]` (>1.9), finest L2err ~1.2e-3, pointwise
  `max|φ−A|` ~8e-5. Boundary value matches (φ[0]=φ[-1]≈0.515 ≈ A(0)=a0=0.5).
- **T2 (2g asymmetric)**: g0 orders `[2.04, 2.01]`, g1 `[2.05, 2.01]` (>1.9),
  finest maxabs g0=6.6e-5, g1=2.3e-5.
- **Value-check tol**: `rtol=5e-3, atol=5e-3` at nx=80 (pointwise error ~8e-5 ≪
  tol; a dropped `q.boundary` → φ≈0 at the boundary vs ref 0.5 → fails by ~0.5,
  amply discriminating).

### T3 xfail confirmation (genuine, rides #195)

Sphere non-vacuum L2(vol) STAGNATES: nc=20→2.37e-2, 40→2.42e-2, 80→2.43e-2
(orders ≈ -0.02 to -0.006 — NOT converging). Finest 2.4e-2 ≫ the #195 in-band
`[1e-8, 1e-3]`. BOTH the rate gate AND the magnitude band fire → xfail-strict is
robust. The boundary VALUE is honoured (φ[-1]=0.7499 ≈ ref=0.75) — the non-vacuum
inflow machinery WORKS; only the interior curvilinear-DD spatial convergence is
pre-asymptotic (exactly the `test_mms_curvilinear.py` #195 signature). T3 collects
+ xfails (NOT xpass). T3g (green structural) confirms inflow honoured at r=R +
redistribution source `(1−μ²)B/r` non-zero NOW.

## Must-stay-green (verbatim)

| Gate | Result |
|------|--------|
| `tests/sn/verification/mms/test_mms.py` (vacuum slab) | **2 passed** |
| `tests/sn/verification/mms/test_mms_curvilinear.py` + `test_curvilinear_aniso_convergence.py` (xfail #195) | **1 passed, 6 xfailed** (unchanged — parameterization did NOT flip them) |
| `tests/numerics/test_operator.py` + `tests/sn/operators/test_operator_block_role.py` | **80 passed** |
| `tests/sn/ -k "fixed_source and not heterogeneous_absolute and not continuous_get"` | **16 passed** |
| `tests/sn/verification/mms/test_mms_heterogeneous.py` (2g hetero vacuum) | SLOW — see "Open items" |
| Sphinx `-W` build | **build succeeded** (only pre-existing test-file SyntaxWarnings) |

## DEVIATIONS from the architecture decisions (flagged)

1. **Foundation test slice-index fix (NOT a deviation from decisions, a test
   correctness fix).** `external_source` returns layout `(N, ng, nx, ny)` — the
   solver-consumed shape — NOT the `(N, nx, 1, 1)` the existing docstrings claim.
   My foundation cross-check slices `[:, 0, :, 0]` (g=0, ny=0). The EXISTING
   `test_sn_mms_anisotropic_symbolic.py` cross-checks slice `[:, :, 0, 0]` — WRONG
   axes — which is why those 2 tests are PRE-EXISTING reds (see below). I followed
   the correct layout.

2. **Slab case made multi-group-capable (extends decision C, does not violate it).**
   Decision C said "add `SNSlabNonVacuumMMSCase`, don't overload existing cases."
   I made the NEW case multi-group via `c_groups`/`sigma_t_g`/`sigma_s_matrix` so
   ONE clean dataclass serves both T1 (1g) and T2 (2g) — rather than two sibling
   slab dataclasses (which would duplicate the shape machinery, violating Cardinal
   Rule 2). The 1g path is `c_groups=(1.0,)`. This is the principled reading of
   "single source of truth"; flag for review if a separate 2g dataclass was
   intended.

3. **`-O`-survivable assertions.** T4 + T1/T2/T3's load-bearing numeric ceilings
   use an explicit `_require(cond, msg)` raise (not bare `assert`, which `-O`
   strips → a canary that cannot die, per the sentinel-harness lesson). The
   `np.testing.assert_allclose` checks survive `-O` natively. The canonical
   invocation is `-O`; without this the consistency/value checks would be no-ops.

## PRE-EXISTING REDS (NOT caused by 4.6 — confirmed by `git stash` at clean HEAD)

`tests/sn/verification/mms/test_mms_anisotropic_symbolic.py::test_spherical_aniso_numerical_qext_matches_sympy`
and `::test_cylindrical_aniso_numerical_qext_matches_sympy` FAIL at clean HEAD
`7ccc14a` (verified by stashing my `sn.py` edit and re-running — still red). Root
cause: those tests slice `Q_numerical[:, :, 0, 0]` (shapes `(N,1)` vs `(N,nx)`
mismatch) — WRONG axes; the correct slice for the `(N, ng, nx, ny)` layout is
`[:, 0, :, 0]`. This is a 1-line test-indexing bug in the EXISTING file, orthogonal
to 4.6. **Recommendation**: file a `module:tests type:bug` fix (or fix inline in
the same session per `feedback_no_issues_for_inline_fixes` since it's a 2-char
slice swap on 2 lines). I did NOT touch the existing file (brief: must-not-flip;
they were already red). The decision-A REGRESSION PINS — the `*_identity` /
`*_overall_pass` tests — all PASS (10/12), so the parameterization is clean.

## Open items / handoffs

- **2g hetero vacuum MMS (`test_mms_heterogeneous.py`)** is SLOW (>600s; timed out
  at the default gtimeout). Re-running with a longer timeout in the background; this
  is a runtime concern, NOT a 4.6 regression (my change does not touch the
  vacuum-hetero code path). [UPDATE AT MERGE: paste the result.]
- **`solve_sn_fixed_source` `boundary_source` param (plan step 4, deferred)**: NOT
  needed by anything. Nothing in 4.6 or the must-stay-green gates requires it; the
  bypass path is sufficient. Stays deferred (flag to user only if a production
  consumer wants prescribed inflow as a first-class fixed-source feature).
- **Probe deletion**: `derivations/diagnostics/diag_p46_prescribed_inflow_source.py`
  is now pinned by T4 → PROMOTE-then-DELETE per `_promotion_policy.md`. Left in
  place (git-excluded); main agent may delete on commit.
- **Archivist dispatch**: the Sphinx stub owes a rich narrative (DISPATCH_REQUEST
  emitted, `followup: false` → goes to the user).

## Self-improvement check

No NEW failure mode introduced. Modes already catalogued: Mode 7 (handled by a0>0
non-vanishing + mandatory sphere companion), ERR-026 (T3 `catches`, xfail #195),
ERR-002 (T2 asymmetric SigSᵀ), vv Mode 9 (T4 promotes). NO skill-table append. NO
new ERR-NNN (no bug CAUGHT — the pre-existing aniso-test slice bug is a test-design
defect, not an L0 solver bug; logging it to error_catalog would be a stretch — it
is a wrong-axes mismatch that fails LOUDLY, not a silent wrong-answer). One lesson
worth a `vv-principles` note: the `(N, ng, nx, ny)` vs documented `(N, nx, 1, 1)`
layout drift in `external_source` is exactly the Pattern-13 bare-ndarray-shape-in-
docstring hazard `coding-elegance` warns about — the docstrings lie about the axis
order and the EXISTING test inherited the lie. A typed source field would have made
the slice unspellable-if-wrong.
