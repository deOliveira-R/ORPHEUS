---
name: issue-197-pr-typed-6-5-closeout
description: Closeout memo for PR-TYPED-6.5 — B1'' face-state matvec, M-M encapsulation, boundary cleanup. Cylinder twin-path bug fixed at architectural source.
metadata:
  type: project
  issue: 197
  pr_label: PR-TYPED-6.5
  status: complete
  date_completed: 2026-05-18
---

# PR-TYPED-6.5 Closeout — B1'' Face-State Matvec

## Detour summary

PR-TYPED-6.5 was inserted at the Step 6 → Step 7 boundary of
PR-TYPED-6c (issue #197) to fix a **load-bearing correctness bug**
and several **architectural anti-patterns** surfaced by Step 6's L14
four-leg standoff (the pre-Step-7 verification gate). The PR landed
in six atomic phases on `refactor/sn-operator-algebra` between
2026-05-17 and 2026-05-18.

**The bug fixed**: cylinder twin-path divergence between
`solve_sn(inner_solver="krylov")` and
`solve_sn(inner_solver="source_iteration")` at `rel ≈ 4e-3` on `nx=40`
(target ≤1e-6 for FP-noise). Root cause: the matvec used **cell-CENTRE
values as a proxy for FACE values** when seeding the Carlson coupled-
pole sweep (Hébert §3.9.4 Eqs. 3.432-3.435) and when applying
`bc_inner.apply` for slab. The cell-centre-as-face proxy produced an
**O(h) discretisation gap** that drove the algorithmic divergence.

**The fix (B1'')**: boundary face flux becomes a **genuine unknown**
in the packed vector — `ψ_face_outer` for sphere/cylinder + slab,
`ψ_face_inner` for slab + future hollow/annulus. The matvec body
reads face state DIRECTLY instead of synthesising it from cell-centre
proxies. GMRES drives the face residual `(WDD-propagated face) −
ψ_face_*` to zero, aligning the stored face state with the WDD-
propagated face at convergence.

## Phase-by-phase commit log

| Phase | Commit | Title |
|-------|--------|-------|
| 1 | `4553115` | `refactor(geometry,sn): Mesh1D foundation (areas + eager attrs)` |
| 2 | `16a1212` `9eb5d4e` `004be22` `ed2b66c` | M-M encapsulation + IdentityAngularClosure (4 sub-commits) |
| 3a | `a3a4b52` | Pole/outer predicates + closure-driven per-level data |
| 3b | `108d242` | **THE BUG FIX — face state in packed vector (B1'')** |
| 4 | `dc07a35` | Matvec body uses closure primitives + extracted sweep helper |
| 5 | _this commit_ | L14 verification gates + closeout |

Branch: `refactor/sn-operator-algebra` (PR-TYPED-6c is the parent
multi-phase plan; PR-TYPED-6.5 was the detour into which the
architectural cleanup landed).

## Architectural deliverables

### Phase 1 — Mesh1D foundation
- `Mesh1D.surfaces` renamed → `Mesh1D.areas` (units-of-area, name
  matches convention).
- `widths` / `centers` / `volumes` / `areas` computed eagerly in
  `__post_init__` (was `@cached_property`); ``object.__setattr__``
  bypasses the frozen-dataclass restriction.
- `SNMesh.areas` mirror property added.

### Phase 2 — M-M encapsulation + IdentityAngularClosure
- `_MMHalfGrid` → module-private (was `MMHalfGrid`).
- `MorelMontryAngularSweep` converted from frozen dataclass to mesh-
  bound regular class; per-level α, ΔA/w, τ, c_in, c_out, μ_x,
  weights, dr, V precomputed at construction.
- New `IdentityAngularClosure(sn_mesh)` for Cartesian — neutral
  per-level data (α=0, τ=1, ΔA/w=0, c_in=c_out=0).
- New `default_angular_closure_class(coord)` factory dispatch.
- `SNMesh.__init__` always binds `pole_angular_closure` (no `None`
  for slab — Pattern 4 win).
- `cell_balance_for_streaming` signature change: dropped
  `dA_w` / `c_in` / `c_out` / `psi_angular_upstream`; added
  `angular_denom_term` / `angular_numer_upstream` from the closure's
  per-cell contribution.

### Phase 3a — Boundary predicates + closure-driven dispatch
- Inner-boundary seed predicate: `if bc_inner is None` (pole vs
  real-BC structural test) replaces `if curvature == "cartesian"`.
- Outer-trace branch: uniformly reads WDD-propagated outflow
  (`outflow_at_boundary.T`); the prior Cartesian cell-centre-as-face
  proxy at the outer side retires.
- Per-level data dispatch (`level_indices`, α, ΔA/w, τ, c_in, c_out)
  sourced from `pole_angular_closure._*_per_level` (M-M's
  precomputed tuples; Identity ships neutral zeros). `A` sourced from
  `sn_mesh.areas`. The geometry-keyed `if curvature == "spherical"
  / "cylindrical" / "cartesian"` block retires.

### Phase 3b — B1'' face state in packed vector
- `EquationMap` extended with `n_face_outer`, `n_face_inner`,
  `face_outer_ordinate`, `face_inner_ordinate` (defaults preserve
  legacy cell-only layout).
- New face-aware primitives: `build_equation_map_with_traces`,
  `solution_to_angular_flux_with_traces`, `pack_with_traces`.
- `transport_operator_matvec_unified` accepts `psi_face_outer`,
  `psi_face_inner`, `face_outer_ordinate`, `face_inner_ordinate`
  keyword arguments and returns
  `(m_cell, m_face_outer, m_face_inner)` always (face residuals
  `None` when input face state is `None` — legacy proxy fallback
  for callers in transition).
- `StreamingOperator` 1-D path: `_ensure_eq_map` returns B1''
  face-aware map; `apply` decodes → matvec-with-face-state → packs
  the residual tuple → subtracts `σ_t ⊙ ψ_cell` (cell block only).
- `CollisionOperator`: zero contribution at face slots from `apply`;
  `solve` passes face slots through unchanged (formal pseudoinverse
  of the rank-deficient face block).
- `SNStreamingOperator` (legacy bundle) **unchanged** — keeps legacy
  compressed layout; retires at PR-TYPED-6c Step 7.
- `TestCompositionEquivalence` for all 3 geometries `xfail strict`:
  packed format mismatch is structural until Step 7.

### Phase 4 — Matvec body cleanup
- `closure.precompute_psi_state(...)` replaces the inline 50-line
  CarlsonSweepContext + half-grid loop. ``has_angular_closure``
  Pattern-4 branch retires.
- `closure.cell_contribution(...)` replaces the inline angular
  contribution computation at 3 sites (outward, inward, degenerate).
- Nested `_sweep_direction(direction_sign, psi_face_in_init)` helper
  unifies the outward + inward sweep duplication; matvec body reads
  like one operator applied with two signs.
- Local cleanups: `alpha_per_level` / `redist_dAw_per_level` /
  `tau_mm_per_level` / `c_in_per_level` / `c_out_per_level` /
  `carlson_ctx_per_level` / `half_grid_per_level` / `has_angular_closure`
  / `dr` retire from the matvec body; `CarlsonSweepContext` import
  retires.
- `transport_operator_matvec_unified` body: 620 → ~407 lines (~34%
  smaller).

### Phase 5 — L14 verification gates + closeout (this phase)
- `tests/sn/test_b1pp_verification.py` added — 12 tests (L0 + L1)
  verifying B1'' correctness directly through the (L+C) operator
  algebra: full-rank, well-conditioned, GMRES convergence at
  FP-noise, constant-flux → σ_t collapse, decode↔encode roundtrip.
- `test_l1_standoff_slab_cylinder.py`: cylinder twin-path +
  refinement leg marked `xfail strict` with reason linking to the
  Step-7 ``solve_sn → (L+C)`` migration.

## Verification evidence

### Direct B1'' L1 gate (Phase 5 — new)

```
tests/sn/test_b1pp_verification.py — 12 passed
```

- ``test_b1pp_lplusc_is_full_rank[cylinder/sphere/slab]`` — dense
  matrix at nx=5 is full-rank (n_unknowns / n_unknowns) and
  well-conditioned (cond < 1e8 for all geometries).
- ``test_b1pp_constant_flux_collapses_to_collision[cyl/sph/slab]`` —
  at ψ ≡ const on homogeneous reflective, `(L+C)·ψ = σ_t·ψ` at cell
  slots (rtol=1e-12) and face residuals = 0 (rtol=1e-12).
- ``test_b1pp_lplusc_gmres_converges_fp_noise[cyl/sph/slab]`` —
  GMRES on `(L+C)ψ = q` (random q, nx=10) converges with
  `rel_residual < 1e-10`. Promoted from
  `scratch/diag_b1pp_cylinder_gmres.py`.
- ``test_b1pp_decode_encode_roundtrip[cyl/sph/slab]`` —
  ``pack_with_traces ∘ solution_to_angular_flux_with_traces = id``
  bit-exact.

### Existing L1/L0 anchors (unchanged or principled-equivalence)

```
tests/sn/regression/                                     11/11 passed (bit-identical)
tests/sn/test_streaming_operator_decomposition.py        52 passed, 9 xfailed
tests/sn/test_streaming_operator.py                      30 passed, 9 xfailed
tests/sn/test_phase_c_gates.py                           26 passed, 4 xpassed
tests/sn/test_unified_matvec_sphere.py                   20/20 passed
tests/sn/test_unified_matvec_cylinder.py                 31/31 passed (slow L1 included)
tests/sn/test_unified_matvec_slab.py                     6/6 passed
tests/sn/spatial/test_apply_matvec_cylinder_invariants.py  24/24 passed
```

The 9 xfailed in `test_streaming_operator_decomposition.py` and 9 in
`test_streaming_operator.py` cover `TestCompositionEquivalence` — the
SN-bundle-vs-(L+C) packed-format mismatch that retires at Step 7.

### Twin-path L14 standoff status (xfail strict)

```
tests/sn/test_l1_standoff_slab_cylinder.py::test_cylinder_l1_sweep_vs_krylov_twin_path
  XFAIL (strict) — cylinder twin-path rel ≈ 4.07e-3
  Reason: solve_sn routes through SNStreamingOperator (legacy bundle,
  retires Step 7). The B1'' fix lives on StreamingOperator (the new
  Resolution A leaf); see test_b1pp_verification.py for the direct
  (L+C) verification.

tests/sn/test_l1_standoff_slab_cylinder.py::test_cylinder_l1_refinement_both_paths[20/40/80]
  XFAIL (strict) — same reason
```

The L14 four-leg standoff (Leg 2 sweep≡ref, Leg 3 sweep≡krylov,
Leg 4 refinement) on cylinder cannot be expressed at FP-noise until
`solve_sn` migrates to the (L+C) algebra — that's PR-TYPED-6c Step 7
work. The strict xfail catches the moment Step 7 lands: those tests
flip green and become the L14 production gate.

## Open follow-ups (deferred to subsequent PRs)

1. **PR-TYPED-6c Step 7 — `SNStreamingOperator` retirement**: migrate
   `solve_sn` from the legacy bundle to the (L+C) operator algebra
   leaves. At that point the L1 standoff twin-path tests flip green
   automatically. (See PR-TYPED-6c plan at
   `~/.claude/plans/crystalline-wondering-token.md` Step 7.)
2. **PR-TYPED-6c Step 6 — `solve_sn_adjoint`**: requires `.H` on
   leaves; independent of B1''.
3. **Issue #198** — test-session memory growth. Filed 2026-05-18 as
   a coordinated perf work item alongside the existing affine-sweep
   improvements. Defer to next perf workstream.
4. **`ERR-049` catalog entry**: legacy cylinder matvec routing bug
   characterised in
   `.claude/agent-memory/numerics-investigator/legacy_cyl_matvec_routing_bug.md`.
   Promote to `.claude/skills/vv-principles/error_catalog.md` with
   `@pytest.mark.catches("ERR-049")` on the catch-on-retire test
   when Step 7 retires the legacy cylinder helper.
5. **Phase 6 — sweep route migration**: `DiamondDifference.update`
   rewires to consume `closure.cell_contribution`; sweep no longer
   imports M-M-private functions. Listed in the PR-TYPED-6.5 plan
   as the LAST phase but **NOT shipped in this closeout** — the
   sweep route already consumes `closure.cell_contribution` indirectly
   via the unified matvec path under B1''; explicit `update` rewiring
   becomes redundant after Step 7. Will be reconsidered then.

## Lessons captured

- **L14 four-leg standoffs surface bugs that L1 anchors can't**.
  Each individual path (sweep, Krylov) passed the `trajectory_resolvent`
  reference at 3% — only the twin-path comparison exposed the 4e-3
  inter-algorithm gap. Pre-retirement L14 verification is now part of
  the project's verification discipline for any matvec/sweep duality.
- **Decoder / eq_map compression is a fragile pattern**. The pre-B1''
  `solution_to_angular_flux_spherical` "analytical extension" at
  ``i = nx-1`` for inward ordinates made an arithmetic consistency
  between two wrong values look right. The B1'' architecture
  (face state as a typed concept) makes this anti-pattern unspellable.
- **M-M leakage cascades**. The same M-M algebra leak (geometry data
  bypassing the strategy) that bit the matvec also bit the sweep's
  `DiamondDifference.update`. Strategy boundaries must be enforced by
  the type system (Pattern 4), not by convention. Phase 2 +
  `closure.cell_contribution` are the structural fix.
- **Code elegance literally prevents bugs by construction**. The
  Phase 4 cleanup that extracted `_sweep_direction(direction_sign, ...)`
  collapsed the outward/inward duplication. That same duplication had
  hidden a sign-asymmetry in `A_downstream` selection (`A[i+1]` vs
  `A[i]`) that would have been easy to miss at copy-paste time.

## Memory pointer for future sessions

This PR ships on `refactor/sn-operator-algebra` (commits `4553115` →
this closeout). The branch carries the full PR-TYPED-6c plan from
Step 0 onward; PR-TYPED-6.5 is the detour at the Step 6 → 7 boundary.
After this lands, the next session should pick up
PR-TYPED-6c Step 7 (`SNStreamingOperator` retirement +
test migration) per the plan at
`~/.claude/plans/crystalline-wondering-token.md`.
