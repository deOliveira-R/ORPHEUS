---
name: issue-196-sn-operator-architecture-audit
description: Architectural audit of the SN operator stack around Issue #196 manifestation #7 (SI-vs-Krylov per-cell drift). Diagnoses why the M-M angular recurrence has two structurally-mirrored implementations (sweep + apply-matvec) and proposes an operator-algebra unification target in which the SI sweep, the Krylov matvec, and the BCs are all compositions of a single set of `LinearOperator` primitives — so that twin-path defects become structurally impossible.
metadata:
  type: project
---

# Issue #196 — SN operator-algebra architectural audit

**Audit dispatched**: 2026-05-12. **Branch**:
`refactor/sn-operator-algebra`. **Trigger**: Issue #168 Phase F closed
the spherical/cylindrical sweep-path twin of a Carlson-seed bug
(`carlson_inward_sweep_from_source` factored, both paths fixed). The
twin-bug existing AT ALL is the architectural smoking gun: under a
correct operator algebra both paths would be applications of the same
operator, so fixing one fixes the other by construction.

## TL;DR — the architectural diagnosis

The operator-algebra layer (`orpheus.numerics.operator`) installs a
clean, mathematically-honest `LinearOperator` Protocol with full algebra
dunders (`+`, `-`, `*`, `@`, `&`, `**`, `.H`, `__call__`) and capability
sets propagating through compositions. The boundary-condition layer
(`orpheus.sn.boundary_realizer`) realises BCs into Wave-0 primitives
(`PermutationOperator`, `IncomingOrdinateMaskTensor`,
`AngularAverageOperator`, `IdentityOperator`, `ZeroOperator`,
`ScaledOperator`) that ARE first-class `LinearOperator`s.

That algebra is **almost completely unused by the load-bearing
transport paths.** The four functions

* `_sweep_1d_spherical` (sweep.py:397-595)
* `_sweep_1d_cylindrical` (sweep.py:602-780)
* `transport_operator_matvec_spherical` (operator.py:571-838)
* `transport_operator_matvec_cylindrical` (operator.py:851-1107)

are 850 lines of **procedural Python for-loops over ordinates and
cells**. They each re-implement, from scratch, the same five
primitives: WDD spatial closure, M-M angular closure, pole-face
boundary condition, outer-face boundary condition, source assembly.
The `SNStreamingOperator.apply()` dispatcher (operator.py:1287-1360)
is a 4-line dispatch shim — its name promises math, its body delivers
verbatim function-call delegation into the procedural matvec.

Manifestation #7 (Issue #196 — residual O(h) SI-vs-Krylov drift on
heterogeneous problems) is an **inevitable consequence** of this
structure: two procedural implementations of the same continuous
operator that build their numerical compositions independently, on
different vector layouts (packed eq-vector for matvec, structured
`(N,nx,ny,ng)` for sweep), with different BC entry points
(face-value `bc.apply(outflow_at_boundary)` in matvec; per-ordinate
`bc.apply(bc_outer)[n]` in sweep). They could not be bit-identical
even if every numerical decision were "correct" — the operator
formulation does not exist at a layer where a single decision can be
shared.

---

## Section 1 — Current state: call graphs

### 1.1 SI dispatch (spherical)

```
solve_sn_fixed_source(materials, mesh, quad, ext_src)              solver.py:924-1064
  │
  ├─ SNMesh(mesh, quad)                                            geometry.py
  ├─ if inner_solver is None and curvature in (sph, cyl):
  │     inner_solver = "krylov"      ◄── auto-flip per Phase D     solver.py:1030-1034
  ├─ SNSolver(materials, sn_mesh, inner_solver=...)                solver.py:151-262
  │     │
  │     ├─ ScatteringOperator.from_solver_data(...)                solver.py:240
  │     ├─ FissionOperator.from_solver_data(...)                   solver.py:253
  │     └─ self.L = SNStreamingOperator(sn_mesh, sig_t)            solver.py:260
  │
  └─ _solve_fixed_source_si(solver, sn_mesh, ext_src, ...)         solver.py:1067-1129
        │
        for n_inner in range(max_inner):
          Q = np.zeros_like(phi)
          solver._add_scattering_source(Q, phi)                    → scattering_op.add_iso_source
          solver._add_n2n_source(Q, phi)                           → scattering_op.add_n2n_source
          Q_aniso_p1 = solver._build_aniso_scattering(angular)     → scattering_op.build_aniso_source
          angular, phi = transport_sweep(Q, sig_t, sn_mesh, ...)   → sweep.py:112
            │
            └─ _curvilinear_sweep(...)                             sweep.py:202
                  │
                  └─ _sweep_1d_spherical(Q, sig_t, sn_mesh,        sweep.py:397-595
                                         psi_bc, Q_aniso)
                        │
                        │  PER-SWEEP-PREAMBLE:                                  ►► OPERATOR #1: BC apply
                        │    inflow_full = bc_outer_obj.apply(bc_outer)         sweep.py:508
                        │    bc_outer_value = inflow_full[most_inward_idx, :]   sweep.py:509
                        │
                        │  PER-SWEEP-PREAMBLE:                                  ►► OPERATOR #2: Carlson seed
                        │    phi_aux = carlson_inward_sweep_from_source(        sweep.py:518
                        │                Q_bar=0.5*Q_1d.T, sigma_t=...,
                        │                dr=..., bc_outer_value=...)            ← psi_half_angle_seed.py:363
                        │    psi_angle = phi_aux.T.copy()                       sweep.py:529
                        │
                        │  PER-ORDINATE LOOP (n = 0 .. N-1):
                        │    mu_n, w_n = mu[n], weights[n]
                        │    QV = QV_iso  [+ Q_aniso_1d[n] * V[:, None]]
                        │
                        │    PER-ORDINATE BC APPLY (inward):                    ►► OPERATOR #1 (AGAIN)
                        │      if mu_n < 0:
                        │        psi_in_full = bc_outer_obj.apply(bc_outer)     sweep.py:556
                        │        psi_spatial_in = psi_in_full[n]                sweep.py:557
                        │      else:
                        │        psi_spatial_in = np.zeros(ng)
                        │
                        │    PER-CELL LOOP (visit = topo order):
                        │      for visit in sn_mesh.iter_cell_visits(ord=n):
                        │        upstream = UpstreamState(spatial_upstream,
                        │                                  angular_upstream)
                        │        result = cell_update.update(                   ►► OPERATOR #3+#4:
                        │            visit, sig_t_1d[i], QV[i], upstream)          WDD + M-M closures
                        │                                                       diamond.py:528-624
                        │
                        │        psi = result.cell_average_flux
                        │        psi_angle[i] = result.outgoing_angular_state   ←   M-M closed half-flux
                        │        angular_flux[n, i, 0] = psi
                        │        scalar_flux[i] += w_n * psi                    ►► OPERATOR #5: ang-int
                        │        psi_spatial_in = result.outgoing_spatial_flux  ←   WDD-closed face
                        │
                        │    POST-SWEEP STORE (outward only):
                        │      if mu_n >= 0: bc_outer[n] = psi_spatial_in       sweep.py:593
```

### 1.2 SI dispatch (cylindrical)

Same shape as 1.1 with an extra outer level loop over
`quad.level_indices` (sweep.py:697); per-level Carlson seed
(sweep.py:705-714); per-level azimuthal recurrence; degenerate
`|η| < 1e-15` case is handled inside `cell_update.update` via
`outgoing_spatial_flux is None` (diamond.py:628-689).

### 1.3 Krylov dispatch (spherical)

```
solve_sn_fixed_source(..., inner_solver="krylov")                  solver.py:1132-1289
  │
  └─ _solve_fixed_source_krylov(solver, sn_mesh, ...)
        │
        │  PRE-LOOP:
        │    eq_map = build_equation_map_spherical(nx, quad, ng)   operator.py:465-488
        │    ext_packed_flat = ext_src ÷ sum_w, packed via eq_map  solver.py:1176-1181
        │    L_scipy = scipy.LinearOperator(matvec=solver.L.apply) solver.py:1186
        │    precond = solver._make_sweep_preconditioner(...)      solver.py:553-625
        │           ◄── this calls transport_sweep INSIDE GMRES,
        │               so the sweep code DOES run on the Krylov path
        │
        │  OUTER LOOP (scattering self-consistency):
        │    for n_outer in range(max_inner):
        │      rhs_iso = _build_rhs_spherical(...)                 solver.py:764-800
        │      rhs = rhs_iso + ext_packed_flat
        │      solution, info = gmres(L_scipy, rhs, ..., M=precond) solver.py:1232
        │             │
        │             └─ scipy iterates:
        │                  L_scipy.matvec(psi_k) = solver.L.apply(psi_k)
        │                                       ↓
        │                                       operator.py:1287-1360
        │                                       SNStreamingOperator.apply()
        │                                       │
        │                                       └─ transport_operator_matvec_spherical(
        │                                            psi, eq_map, quad, sig_t,
        │                                            nx, ng, face_areas, volumes,
        │                                            alpha_half, redist_dAw, tau_mm,
        │                                            sn_mesh=..., bc_outer=...,
        │                                            pole_angular_closure=...)
        │                                          operator.py:571-838
        │                                          │
        │                                          │ PREAMBLE:
        │                                          │   if bc_outer is None:
        │                                          │     bc_outer = realize(spec_law)         ►► OPERATOR #1: BC
        │                                          │   if pole_angular_closure is None:
        │                                          │     pole_angular_closure =
        │                                          │       LegacyTauSymmetricInterp.()
        │                                          │   fi = solution_to_angular_flux_sph(...) operator.py:491-568
        │                                          │
        │                                          │ CARLSON CONTEXT BUILD:
        │                                          │   outer_inflow_est =                     ►► OPERATOR #1 (AGAIN)
        │                                          │     bc_outer.apply(fi[:, :, -1, 0].T)    operator.py:713
        │                                          │   bc_outer_value = outer_inflow_est[
        │                                          │                       most_inward_idx, :]
        │                                          │   carlson_ctx = CarlsonSweepContext(
        │                                          │     sigma_t, dr, mu, w, bc_outer_value)
        │                                          │
        │                                          │ POLE CLOSURE / ANGULAR REDISTRIB:
        │                                          │   redist_full =                          ►► OPERATOR #3: pole closure
        │                                          │     pole_angular_closure(                 (M-M angular recurrence)
        │                                          │       fi[..., 0], alpha_half,
        │                                          │       redist_dAw, tau_mm, V,
        │                                          │       carlson_context=carlson_ctx)
        │                                          │
        │                                          │ OUTWARD SWEEP (μ > 0):
        │                                          │   psi_face_in = fi[:, outgoing, i0, 0]   ◄ pole-face anchor
        │                                          │   for visit in iter_cells_by_direction(+1):
        │                                          │     psi_cell = fi[:, outgoing, i, 0]
        │                                          │     psi_face_out = 2*psi_cell -          ►► OPERATOR #4: WDD
        │                                          │                    psi_face_in              (inlined here, NOT via cell_update)
        │                                          │     streaming = mu*(A[i+1]*psi_face_out
        │                                          │                  - A[i]*psi_face_in) / V[i]
        │                                          │     redistribution = redist_full[..., i]
        │                                          │     collision = sig_t[i] * psi_cell
        │                                          │     lhs[:, ks] = streaming + redistribution + collision
        │                                          │     psi_face_in = psi_face_out
        │                                          │   outflow_at_boundary[:, outgoing] = psi_face_out
        │                                          │
        │                                          │ BC APPLY (face-trace contract):          ►► OPERATOR #1 (FINAL CALL)
        │                                          │   inflow_full =                          operator.py:815
        │                                          │     bc_outer.apply(outflow_at_boundary.T)
        │                                          │
        │                                          │ INWARD SWEEP (μ < 0):
        │                                          │   psi_face_in = inflow_full[incoming, :].T
        │                                          │   for visit in iter_cells_by_direction(-1):
        │                                          │     [...analogous WDD recurrence, A[i+1]→psi_face_in,
        │                                          │     A[i]→psi_face_out...]
        │                                          │
        │                                          └─ return lhs.ravel(order='F')
        │
        │      fi = solution_to_angular_flux_spherical(solution, ...)
        │      phi = _scalar_flux_from_angular(fi, ...)
        │      [convergence check]
```

### 1.4 Krylov dispatch (cylindrical)

Same shape as 1.3 plus per-level loop (operator.py:990-1090) and a
degenerate-azimuthal Phase 3 block (operator.py:1092-1105). The
`pole_angular_closure` is called ONCE before the level loop with a
per-level list of `CarlsonSweepContext` (operator.py:976-985); the
WDD spatial recurrence is inlined per-level inside the matvec.

### 1.5 Krylov is doubly procedural

The Krylov path is NOT a clean separation between "operator" and
"inverse" — the preconditioner inside GMRES (`_make_sweep_preconditioner`,
solver.py:553-625) IS the sweep code (`transport_sweep`) called with
`Q_iso = 0` and `Q_aniso` decoded from the packed RHS. So the sweep
math runs both on the SI path AND inside GMRES as a left
preconditioner. Manifestation #7's residual drift is the gap between

* `transport_operator_matvec_spherical`'s **inlined** WDD recurrence
  (operator.py:792, 825), and
* `DiamondDifference._update_curvilinear`'s **factored** WDD recurrence
  via `cell_update.update(...)` (diamond.py:613-618).

These are two implementations of `ψ_face_out = 2·ψ_avg − ψ_face_in` ON
DIFFERENT VECTOR LAYOUTS, with different anchor conventions at the
pole face (matvec: cell-centre value `fi[:, outgoing_mask, i0, 0]`,
operator.py:784; sweep: implicitly `psi_spatial_in = zeros(ng)` for
outward + `cell_update`'s WDD then closing `psi_face_out = 2·ψ_avg −
ψ_face_in`).

---

## Section 2 — Duplication inventory

For every operator-algebra concept needed to apply L (streaming +
collision) and L⁻¹ (sweep), there are TWO independent implementations.

| #  | Operator concept                  | SI / sweep path                                                                | Krylov / apply-matvec path                                                   | Bit-identical on flat ψ? | Diverges on non-flat ψ?                                                                                                                                                |
|----|-----------------------------------|---------------------------------------------------------------------------------|------------------------------------------------------------------------------|---------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 1  | **WDD spatial closure**           | `DiamondDifference._update_curvilinear` (diamond.py:613-618) — via cell-balance form `ψ_avg = numer/denom`, then `ψ_face_out = 2·ψ_avg − ψ_face_in` | Inlined `psi_face_out = 2.0 * psi_cell - psi_face_in` (operator.py:792, 825). `psi_cell` here is `fi[:, outgoing, i, 0]` ← READS FROM INPUT VECTOR, NOT FROM A SOLVE | **NO** — sweep computes `ψ_cell` from balance equation, matvec reads `ψ_cell` from input vector | **YES** — this is the fundamental algebraic gap. Sweep applies L⁻¹·q; matvec applies L. They're algebraically inverses, not algebraically equal. Manifestation #7's residual drift comes from the layout mismatch (packed vs structured) — the same operator concept implemented on different storage |
| 2  | **M-M angular recurrence**        | Inside `cell_update.update` for sphere/cylinder: `psi_angle_out = (psi_avg − (1−τ)·psi_angle_in)/τ` (diamond.py:618) | Inside `pole_angular_closure(...)` → `MorelMontryAngularSweep` → DD-style recurrence `phi_{n+1/2} = (phi_n − (1−τ)·phi_{n−1/2})/τ` evaluated over (g, n, i) once before the cell loop | YES (post-Phase F)        | YES if seed differs. Phase F backported Carlson seed to BOTH paths via `carlson_inward_sweep_from_source` → ROOT CAUSE: TWO IMPLEMENTATIONS                              |
| 3  | **Carlson μ=−1 seed**             | `carlson_inward_sweep_from_source(0.5*Q_1d.T, sigma_t, dr, bc_outer_value)` called once per sweep, sweep.py:518 | `MorelMontryAngularSweep.psi_half_angle_seed = CarlsonInwardSweep()`. CarlsonInwardSweep.__call__ ALSO calls `carlson_inward_sweep_from_source` (psi_half_angle_seed.py:574-579) after folding `Q_bar = 0.5 · σ_t · φ_0` | YES post Phase F          | Same recurrence both sides via shared helper (#168 Phase F). But the **call sites are different**, the **inputs are different** (sweep: `0.5·Q_1d`; matvec: `0.5·σ_t·φ_0`), and the **algebraic invariant** (`σ_t·φ_0 = Q_1d` at the fixed point) is the only thing that keeps them consistent — sub-tle, fragile |
| 4  | **BC apply (outer face)**         | `bc_outer_obj.apply(bc_outer)` where `bc_outer` is the persistent `(N, ng)` outflow buffer, sweep.py:508, 556; called once for the Carlson context AND once **per inward ordinate** | `bc_outer.apply(outflow_at_boundary.T)` where `outflow_at_boundary` is the WDD-propagated `(N, ng)` face vector; called once globally at end of outward phase (operator.py:815); ALSO `bc_outer.apply(fi[:, :, -1, 0].T)` on cell-centred values for the Carlson context build (operator.py:713) | YES                       | YES on the **Carlson context build** — matvec uses `fi[:, :, -1, 0].T` (cell-centre values), sweep uses the persistent `bc_outer` outflow buffer. Different inputs → different `bc_outer_value` → different Carlson seed → different angular redistribution. This is structural, NOT a numerical accident |
| 5  | **Pole-face anchor (μ>0 inward sweep start)** | `psi_spatial_in = zeros(ng)` for `mu_n >= 0` (sweep.py:559) → A[0]·0 = 0 at r = 0 streaming term cancels | `psi_face_in = fi[:, outgoing_mask, i0, 0].copy()` (cell-centre value, operator.py:784) — analytical extension making flat-ψ recurrence preserve flat-flux invariant | YES (both give A[0]·anchor = 0 because A[0] = 0) | YES on non-flat ψ: matvec's anchor is the cell-centre value, sweep's anchor is zero. The streaming-term coefficient A[0]=0 absorbs the difference at the pole face, but only in the **streaming contribution** — the WDD recurrence at cell i=0 propagates the anchor downstream via `psi_face_out = 2·psi_cell − psi_face_in`, so a different anchor produces a different cell-1 face flux |
| 6  | **Source assembly**              | `Q_1d = Q[:, 0, :]` is fed direct (already includes scatter + ext); `QV_iso = Q_1d * V[:, None] * weight_norm` (sweep.py:476); per-ordinate `QV = QV_iso + Q_aniso_1d[n] * V[:, None]` | `rhs_iso = _build_rhs_spherical(zeros, phi, ...)` calls `_build_rhs_spherical` which manually builds `qS = sig_s[mid][0].T @ phi_cell / sum_w` + `qF + q2` per (ordinate, cell) (solver.py:764-800); separately, `ext_packed_flat = external_source[ord, ix, iy, :] / sum_w` (solver.py:1176-1181) | YES                       | The scatter source builds via two DIFFERENT code paths: `ScatteringOperator.add_iso_source` on the SI side; `_build_rhs_spherical`'s inline loop on the Krylov side. Both eventually compute `sig_s @ phi / sum_w`, but the cycle of normalization is different (SI: `Q_1d / sum_w` inside sweep; Krylov: `/sum_w` at build_rhs time, `× sum_w` undone in `_make_sweep_preconditioner`). Easy to drift                              |
| 7  | **Angular integration φ = Σ wₙ ψₙ** | Inline inside sweep `scalar_flux[i] += w_n * psi` (sweep.py:583, 763) | `_scalar_flux_from_angular(fi, quad, nx, ny, ng)` — loop over (iy, ix) doing `sf[ix,iy,:] = sum(fi[:,:,ix,iy] * weights[None,:], axis=1)` (solver.py:655-676) | YES                       | YES — both compute the same Σ to round-off. But it IS a duplication of the same conceptual operator (the "AngularIntegrationOperator" = `DiagonalOperator(weights) @ Σ` over ordinate axis) and the `numerics.operator.DiagonalOperator` (operator_algebra.rst §270-316) is exactly the type for this — used nowhere here                       |
| 8  | **Ordinate iteration order (DAG)** | `sn_mesh.iter_cell_visits(ordinate_idx=n)` — per-ordinate topological order, sweep.py:566 | `sn_mesh.iter_cells_by_direction(+1)` / `iter_cells_by_direction(-1)` — per-sign-class topological order, operator.py:752-753, 1009, 1063 | YES                       | NO — these are bit-identical when given a representative ordinate (foundation test pinned, `test_iter_cells_by_direction.py`). But they're TWO APIs on `SNMesh` for what is mathematically the same `Σ_{cell∈cells_in_dag_order(direction)} (cell-update on cell)` quantifier              |
| 9  | **Vector layout**                | Structured `(N, nx, ny, ng)` angular flux + `(nx, ny, ng)` scalar flux        | Packed 1-D `(n_unknowns,)` solution vector with `eq_map` selecting which `(ord, ix, iy)` are unknowns; F-order `solution.reshape(ng, n_eq, order='F')` | n/a — different shapes    | The packed layout EXISTS because `scipy.sparse.linalg.LinearOperator.matvec` wants a 1-D vector. The conversion is `solution_to_angular_flux*` (operator.py:491-568) which does its own analytical-extension fill at the outer cell for inward ordinates (operator.py:564-568); the sweep never does this conversion |
| 10 | **Half-angle face flux state**   | `psi_angle = phi_aux.T.copy()` — `(nx, ng)` state mutated per ordinate by `cell_update` (sweep.py:529, 580) | `redist_full = pole_angular_closure(...)` — precomputed `(ng, N, nx)` array of all angular redistribution contributions before the cell loop (operator.py:737-741) | YES (post-Phase F)         | The sweep ITERATES the M-M recurrence inside the cell loop; the matvec PRECOMPUTES the entire redistribution array via a separate strategy. Same math, two evaluation orders, two storage layouts. The matvec form composes cleanly with linearity (it's a single `LinearOperator`-style call); the sweep form is locked inside an imperative loop |

**Pattern**: every row of the table shows the same operator concept
implemented twice on different data structures. The duplication is
structural, not accidental — it reflects the fact that
`SNStreamingOperator.apply()` and `SNStreamingOperator.solve()`
historically were built for different consumers (BiCGSTAB vs source
iteration) and never unified. Wave H's Phase 0 algebra (`numerics/
operator.py`) is the foundation for the unification but the
unification itself has not happened on the load-bearing transport
paths.

---

## Section 3 — Operator-algebra surface audit

### 3.1 What dunders exist (and where)

| Module / class                                   | `__call__` | `__add__` | `__sub__` | `__mul__` | `__matmul__` | `__and__` (⊗) | `__neg__` | `__pow__` | `.H` (adjoint) | `apply` | `solve` | `apply_transpose` |
|--------------------------------------------------|:----------:|:---------:|:---------:|:---------:|:-------------:|:--------------:|:----------:|:----------:|:----------------:|:-------:|:-------:|:-----------------:|
| `LinearOperatorMixin` (numerics/operator.py)     | ✓          | ✓         | ✓         | ✓         | ✓             | ✓              | ✓          | ✓          | ✓                | (subclass) | (subclass) | (subclass) |
| `SNStreamingOperator` (sn/operator.py:1115)      | (mixin)    | (mixin)   | (mixin)   | (mixin)   | (mixin)       | (mixin)        | (mixin)    | (mixin)    | (mixin)          | ✓       | ✓ (via sweep) | ✓ (dense-probe) |
| `ScatteringOperator` (sn/scattering.py)          | (mixin)    | (mixin)   | (mixin)   | (mixin)   | (mixin)       | (mixin)        | (mixin)    | (mixin)    | (mixin)          | ✓       | —       | —                 |
| `FissionOperator` (sn/fission.py)                | (mixin)    | (mixin)   | (mixin)   | (mixin)   | (mixin)       | (mixin)        | (mixin)    | (mixin)    | (mixin)          | ✓       | —       | —                 |
| `PermutationOperator` (BC reflective realisation) | (mixin)   | (mixin)   | (mixin)   | (mixin)   | (mixin)       | (mixin)        | (mixin)    | (mixin)    | (mixin)          | ✓       | —       | ✓                 |
| `IncomingOrdinateMaskTensor` (BC vacuum)         | (mixin)    | (mixin)   | (mixin)   | (mixin)   | (mixin)       | (mixin)        | (mixin)    | (mixin)    | (mixin)          | ✓       | —       | ✓                 |
| `AngularAverageOperator` (BC white)              | (mixin)    | (mixin)   | (mixin)   | (mixin)   | (mixin)       | (mixin)        | (mixin)    | (mixin)    | (mixin)          | ✓       | —       | —                 |
| `IncomingSourceOperator` (BC prescribed inflow)  | (mixin)    | (mixin)   | (mixin)   | (mixin)   | (mixin)       | (mixin)        | (mixin)    | (mixin)    | (mixin)          | ✓       | —       | —                 |

The algebra surface is **broad and complete**. Every type that
participates can be summed, scaled, composed, tensored, adjointed,
called as `op(x)` for math-like notation. The capability set
propagation is enforced at composition time with `MissingCapability`
raised cleanly (numerics/operator.py:93-99).

### 3.2 What dunders are USED by the load-bearing paths

```
$ grep -rn '__matmul__\|@ self\.\|@ L\|@ S\|@ F\|self\.L @\|self\.S @\|self\.F @\|L_op @\|S_op @\|F_op @' orpheus/sn/*.py
```

Result: **zero hits**. Not a single `@` operator is used to compose
SN operators. `(L - S - F)·ψ` is documented in
`docs/theory/operator_algebra.rst:Eq.operator-fixed-source` as the
target form — and not used anywhere in the SN solver paths.

```
$ grep -rn '__call__\|\.L(\|\.S(\|\.F(\|sn_op(' orpheus/sn/*.py
```

The Mixin's `__call__` (numerics/operator.py:326-333) is a
universal alias for `apply`, intended to enable
`(L - S) @ psi` or `(L - S)(psi)` notation. It is never called from
the load-bearing paths — they all call `.apply(...)` explicitly, or
they bypass the algebra entirely and call the underlying
`transport_operator_matvec_*` / `transport_sweep` functions.

### 3.3 What's MISSING (operators that don't exist but should)

For the algebra to read like math at the SI/Krylov dispatch level,
the following operator types need to grow:

| Missing primitive | Mathematical role | Where current code re-implements it |
|---|---|---|
| **`CellUpdateOperator`** | `L_cell⁻¹` on a single cell — affine in `(spatial_upstream, angular_upstream, source)` | Built imperatively inside `DiamondDifference.update` (diamond.py:528-624) and AGAIN inside `transport_operator_matvec_spherical` via inlined `2·psi_cell − psi_face_in` (operator.py:792). The `CellUpdate` Protocol exists but is NOT a `LinearOperator` — it's a "thing with an `update(...)` method", not an algebra citizen |
| **`SweepOperator`** | `L⁻¹` on the full mesh — composition of per-cell `CellUpdateOperator` along DAG topological order, with BC operator inserted at the inflow face | `transport_sweep` (sweep.py:112-153) is a free function returning a tuple; it composes nothing in the algebra sense — it is `L⁻¹` only by convention |
| **`AngularRedistributionOperator`** | The M-M closed redistribution `(ΔA/w)·(α_{n+½}·φ_{n+½} − α_{n−½}·φ_{n−½})/V`. Linear in ψ via the Carlson-seeded M-M recurrence | `MorelMontryAngularSweep.__call__` (in pole_angular_closure.py) — NOT a LinearOperator subclass (it implements `__call__` but doesn't inherit `LinearOperatorMixin` or advertise capabilities) |
| **`StreamingOperator`** (curvilinear, sweep-direction-aware) | The pure `μ·(A∂/∂r)/V` action with WDD closure | Inlined per-sweep-direction inside `transport_operator_matvec_spherical` lines 791-803 / 824-836. Symmetric to inlined inside `cell_update` (diamond.py:601-616) |
| **`AngularIntegrationOperator`** (`Φ = ∫ ψ dΩ`) | `DiagonalOperator(weights=quad.weights, axis=ordinate_axis)` composed with a sum-along-ordinate reduction. Exactly the use-case `DiagonalOperator` was built for (operator_algebra.rst:270-316) | `_scalar_flux_from_angular` (solver.py:655-676) on the matvec side; inline `scalar_flux[i] += w_n * psi` on the sweep side |
| **`VectorLayoutAdapter`** (PackedToStructured / StructuredToPacked) | A `LinearOperator` that maps packed `(n_unknowns,)` ↔ structured `(N, nx, ny, ng)` via `eq_map`. Linear by construction. Carries the BC-fill-from-outflow logic | `solution_to_angular_flux_*` + the inverse code path in `_make_sweep_preconditioner` (solver.py:606-623). Two independent ad-hoc implementations |
| **`BoundaryRealizationOperator`** at the **face-trace level** vs. per-ordinate level — needs a uniform contract for "apply BC to the boundary outflow trace and emit the boundary inflow trace" | Already exists at the leaf (`PermutationOperator.apply`, `IncomingOrdinateMaskTensor.apply` etc. — single-arg, axis-aware) | BC.apply is called THREE different ways on the matvec path (operator.py:713, 815, 954) and TWO different ways on the sweep path (sweep.py:508, 556) — each call site does its own indexing dance to extract a specific ordinate's inflow value |

### 3.4 What dunders exist but are dormant

`SNStreamingOperator.solve` (operator.py:1364-1407) exists and
delegates to `transport_sweep`. It is **never called by Krylov**:
the preconditioner inside `_make_sweep_preconditioner` calls
`transport_sweep` DIRECTLY (solver.py:606-610), not through
`self.L.solve`. So the operator algebra's `solve` capability is
declared but bypassed by the load-bearing consumer.

`SNStreamingOperator.apply_transpose` (operator.py:1442-1472) is
implemented via dense-matrix probing (one-time O(n²) construction).
It works on small test problems but is not used by any production
code path. The algebra advertises `CAP_APPLY_TRANSPOSE` but the
adjoint flux solver hasn't been written.

`L - S - F` as an algebraic expression: the `__sub__` dunder gives
`OperatorSum(L, ScaledOperator(-1, S))`. Constructing
`(self.L - self.S - self.F)` would produce a valid
`LinearOperator` with `apply` capability — and it is exactly what
`Lψ = (1/k)·F·ψ + S·ψ + q` would consume at the iteration-driver
level. But the production paths build their per-iteration RHS
manually (solver.py:401-403, 1097-1099) rather than calling the
composed operator.

---

## Section 4 — Architectural diagnosis

**The operator formulation is inconsistent in TWO places, not one.**

### 4.1 Primary inconsistency: TWO implementations of the same continuous operator

`SNStreamingOperator.apply` and `SNStreamingOperator.solve` are
documented in the class docstring (operator.py:1150-1180) as
implementing **different closures** of the same operator —
"symmetric" (apply) vs "WDD asymmetric" (solve). The docstring
labels this "by design".

**This is the architectural smoking gun.** A linear operator has one
discrete representation. Saying `apply` and `solve` use "different
closures" is saying `L` and `L⁻¹` are different operators — which
cannot be true under operator algebra. What it actually means is

* `apply` represents L_sym (the symmetric-closure discretisation),
* `solve` represents L_WDD⁻¹ (the WDD-closure discretisation's
  inverse),

and **the sweep solves L_WDD ψ = q while Krylov-on-apply solves L_sym
ψ = q**. Two operators, both calling themselves "L", with two
different fixed points whose difference is Manifestation #7's residual
O(h) drift.

Post-Phase-F: the Carlson seed (formerly the only ERR-026-relevant
divergence) is now consistent across both paths — `apply` and `solve`
now produce angular redistribution from the SAME M-M recurrence with
the SAME seed. The remaining O(h) gap is the WDD-vs-symmetric spatial
closure difference: matvec inlines `psi_face_out = 2·psi_cell −
psi_face_in` where `psi_cell` comes from the input ψ; sweep computes
`psi_cell = numer/denom` from the balance equation, then closes
`psi_face_out = 2·psi_cell − psi_face_in`. These give the same
flat-flux fixed point but the per-iteration residual differs on
non-flat ψ; the Krylov iteration converges to the symmetric-closure
fixed point and the SI iteration to the WDD-closure fixed point.

### 4.2 Secondary inconsistency: load-bearing paths don't compose operators

Even if the apply and solve closures were unified, the SI/Krylov
dispatch wouldn't be expressed in operator algebra. The SI loop is a
hand-coded fixed-point iteration on `Q + S·φ + F·φ/k`; the Krylov
dispatch is a hand-coded GMRES call with manually-built RHS. Neither
reads `(L - S - F) @ psi = q` (the documented form in
`operator_algebra.rst:Eq.operator-fixed-source`).

The Wave E Round 1 iteration primitives (`SourceIteration`,
`PowerIteration` per the docs at `operator_algebra.rst:222-239`) are
declared in plans but the solver paths bypass them.

### 4.3 Tertiary inconsistency: BC is realised but not woven in

The boundary realiser produces clean Wave-0 `LinearOperator`
instances per the table in `boundary_realizer.py:14-44`. But the
matvec and sweep each call `bc.apply(...)` THREE-to-FIVE times per
iteration with different inputs (face vector, cell-centre vector,
per-ordinate slice), choosing which inflow ordinate's value to
extract by ad-hoc indexing (e.g. `outer_inflow_estimate[
most_inward_global_idx, :]`, operator.py:716-717). The user's exact
phrasing — "both the Krylov and the sweep are application of the
Operators including Boundary Operators" — is the principled target;
the current code uses BC as a stateless callable, not as a
first-class algebra primitive.

---

## Section 5 — Proposed unification

### 5.1 New / promoted operator types

| Operator                                     | Type                                | Capabilities                             | Role                                                                                  |
|----------------------------------------------|-------------------------------------|------------------------------------------|----------------------------------------------------------------------------------------|
| `SNCellOperator`                             | `LinearOperator`                    | `{apply, solve}`                         | The 8-line affine cell-balance op. `apply(psi_cell)` returns the residual `L_cell·ψ - q`. `solve(...)` returns `ψ_cell` from `(spatial_upstream, angular_upstream, source)` (current `DiamondDifference.update`) |
| `SNSweepOperator`                            | `LinearOperator`                    | `{apply, solve}`                         | The factored sweep. `apply(ψ)` performs ONE complete sweep iteration on a global flux vector (= what `transport_operator_matvec_spherical` currently does). `solve(q)` performs the source-iteration-style sweep (= what `transport_sweep` currently does). Both compose `SNCellOperator` over `iter_cells_by_direction(±1)` |
| `AngularRedistribution`                      | `LinearOperator`                    | `{apply}`                                | The `(ΔA/w)·(α·Δφ)/V` term as a linear op on `fi`. Wraps the current `MorelMontryAngularSweep.__call__` and `LinearOperatorMixin`. Composes the Carlson seed as an inner operator (also linear: `CarlsonInwardSweep` from psi_half_angle_seed.py) |
| `PoleFaceAnchor`                             | `LinearOperator`                    | `{apply}`                                | The `A[0]·ψ(0) = 0` symmetry condition as an explicit identity-on-the-pole-cell-centre operator. Linear by construction. Replaces the current ad-hoc `psi_face_in = zeros(ng)` / `fi[:, outgoing, i0, 0].copy()` decisions |
| `AngularIntegration`                         | `LinearOperator` (already exists via `DiagonalOperator @ ReductionOperator`) | `{apply, apply_transpose}` | `φ = ∫ ψ dΩ = Σ wₙ ψₙ`. Replaces `_scalar_flux_from_angular` and the inline `scalar_flux[i] += w_n * psi` |
| `PackedStructuredAdapter`                   | `LinearOperator`                    | `{apply, apply_transpose}`               | Pure layout adapter: packed `(n_unknowns,)` ↔ structured `(N, nx, ny, ng)` via `eq_map`. Linear because the BC fill is itself a `bc.apply()` composition. Currently inlined as `solution_to_angular_flux*` (matvec side) and the inverse code path in `_make_sweep_preconditioner` (Krylov side) |
| `SNBoundaryFaceTraceOperator`               | `LinearOperator`                    | `{apply}`                                | The user's principled target: the BC realised at the face-trace level, taking a `(N, ng)` outflow trace and returning a `(N, ng)` inflow trace. Already EXISTS — it's just `bc_outer` after realisation. The new contract is: **all** call sites consume it through the same `apply` shape, no per-ordinate slicing |

### 5.2 What `_sweep_1d_spherical` collapses to under the unified architecture

Current: 199 lines, 16 named local intermediates, 2 nested loops,
3 BC `apply` call sites.

Target:

```python
def _sweep_1d_spherical(Q, sig_t, sn_mesh, psi_bc, Q_aniso=None):
    """SI sweep = one solve of the sweep operator on the within-group source."""
    L_sweep = sn_mesh.sweep_operator(sig_t)  # = SNSweepOperator
    q_full = (Q + Q_aniso) * sn_mesh.volume_measure  # source assembly
    psi = L_sweep.solve(q_full, psi_bc=psi_bc)
    phi = sn_mesh.angular_integration(psi)
    return psi, phi
```

The for-loops over ordinates and cells now live INSIDE
`SNSweepOperator.solve`, which composes:

```python
class SNSweepOperator(LinearOperatorMixin):
    capabilities = frozenset({CAP_APPLY, CAP_SOLVE})

    def __init__(self, sn_mesh, sig_t):
        self.sn_mesh = sn_mesh
        self.cell_op = SNCellOperator(sig_t, sn_mesh.cell_update)
        self.bc_outer = sn_mesh.bc_right
        self.pole_anchor = PoleFaceAnchor(sn_mesh)
        self.redistribution = AngularRedistribution(sn_mesh)  # M-M closed

    def solve(self, q, psi_bc):
        # ψ_half = Carlson seed (currently linear, factored)
        psi_half = self.redistribution.seed(q, self.bc_outer)
        psi_face = self.pole_anchor(psi_half)
        # Sweep direction = +1 (outward), then BC apply, then -1 (inward)
        psi_face = self._sweep_outward(q, psi_half, psi_face)
        psi_face = self.bc_outer @ psi_face  # face-trace BC
        psi_face = self._sweep_inward(q, psi_half, psi_face)
        return self._assemble_global_flux(psi_face, psi_half)
```

Each `_sweep_outward` / `_sweep_inward` is `for visit in
iter_cells_by_direction(±1): psi_face = self.cell_op.solve(visit, q,
psi_face, psi_half)`. The math IS the M-M angular recurrence + WDD
spatial recurrence — exactly what `DiamondDifference._update_curvilinear`
already does. The sweep operator is just composition.

### 5.3 What `_solve_fixed_source_si` collapses to

```python
def _solve_fixed_source_si(solver, sn_mesh, external_source, ...):
    """Fixed-source SI = within-group iteration on (L - S - F⁻¹·F)·ψ = q."""
    L = solver.L   # SNStreamingOperator -- still curvilinear, post-Phase-F clean
    S = solver.S
    F = solver.F
    q_ext = external_source  # per-ordinate, source units already
    return SourceIteration(L, S, q_ext).solve(  # numerics.iteration primitive
        max_iter=max_inner, tol=inner_tol,
    )
```

Where `SourceIteration` is the primitive declared but not yet wired
in (`docs/theory/operator_algebra.rst:222-239`).

### 5.4 What `_solve_fixed_source_krylov` collapses to

```python
def _solve_fixed_source_krylov(solver, ...):
    """Fixed-source Krylov = GMRES on (L - S)·ψ = q_ext + F·φ/k."""
    op = (solver.L - solver.S) @ AngularToScalarAdapter(solver.quad)
    precond = solver.L  # the sweep is L⁻¹
    return PreconditionedGMRES(op, precond).solve(q_ext, ...)
```

Or equivalently, since outer-source-iteration is `SourceIteration`
wrapping a Krylov solve of `L·ψ = (S·φ + q_ext)`:

```python
return SourceIteration(L=solver.L, S=solver.S, q=q_ext,
                       inner_solver=PreconditionedGMRES(...)).solve()
```

### 5.5 How BC operator algebra plugs in

BC operators already ARE `LinearOperator` instances per the
realisation table. The unification is: every consumer (`SNSweepOperator`,
`SNStreamingOperator.apply`, `PackedStructuredAdapter`) accepts a
`SNBoundaryFaceTraceOperator` (= a realised BC) as a constructor
argument, NOT as an embedded call-site `bc.apply(...)`. The
operator-composition expression then carries the BC explicitly:

```
L_sweep = composition_of(
    PackedToStructured(eq_map),
    SNSweepOperator(sn_mesh, sig_t, bc_face_trace=sn_mesh.bc_right),
    StructuredToPacked(eq_map),
)
```

Read this as: pack → sweep (which internally consumes its
`bc_face_trace` once at the boundary edge) → unpack. The BC operator
is composed at the `__init__` site, not called per-ordinate inside
the sweep body.

### 5.6 Does manifestation #7 dissolve by construction?

**Partially.** The unification eliminates the *structural* cause of
twin-path drift (one operator, one composition, one BC entry point
per direction). What remains is the **numerical** difference between
the symmetric-closure and WDD-asymmetric-closure discretisations.
Under the proposed architecture this becomes a single, explicit
choice: which closure does `SNCellOperator` use? If both
`apply` and `solve` route through the same `SNCellOperator`
configured with the same closure (e.g. WDD), then Manifestation #7
**dissolves exactly** — bit-identical fixed points by construction,
not by careful synchronisation. If they're configured with
different closures (deliberately, e.g. for the Wave-E reconciliation
plan), then the Manifestation #7 drift is **principled** — it's the
documented gap between two discretisations of the continuous
operator, not a twin-bug.

Phase F's "viable options (a) sweep WDD-closure refinement / (b)
flip default to Krylov" become **operator-configuration decisions**
under the unified architecture, not code rewrites.

---

## Section 6 — Migration order

The campaign's Wave H acceptance criteria (LinearOperator,
DiscreteMeasure, ReducedStreamingOperator, CellUpdate are foundations,
not demolition targets) are preserved. The migration extends them.

### Step 1 (small, ~2-3 days). Make `CellUpdate` a `LinearOperator`.

`DiamondDifference` currently implements `update(visit, total_xs,
source, upstream_state) → CellResult`. Wrap it (or rebase it) onto
`LinearOperatorMixin` with the structured input
`(spatial_upstream, angular_upstream, source)` → output
`(cell_average, outgoing_spatial, outgoing_angular)`. Capabilities:
`{apply, solve}` where `apply(cell_average) → residual` is the cell
operator action and `solve(source, upstream) → cell_average` is the
existing update. No mathematical change; pure type-system promotion.

**Acceptance criterion**: hand-calc tests
(`tests/sn/spatial/test_diamond.py`) still pass at `np.array_equal`.
**Independently committable.** **Tests**: `pytest tests/sn/spatial/`.

### Step 2 (small, ~2 days). Make `MorelMontryAngularSweep` a `LinearOperator`.

The pole_angular_closure Protocol family becomes
`AngularRedistribution(LinearOperatorMixin)`. The `__call__(psi_level,
α, redist_dAw, τ, V, level_indices, carlson_context)` becomes
`apply(psi_level)` where the geometry coefficients and carlson_context
are bound at `__init__`. `CarlsonInwardSweep` similarly wraps as a
`LinearOperator` with `apply(psi_level, carlson_ctx)` →
half-angle seed.

**Acceptance criterion**: 4 closure-foundation tests
(`tests/sn/spatial/test_pole_angular_closure.py`) green.
**Tests**: `pytest tests/sn/spatial/`.

### Step 3 (small, ~3 days). Promote the SI sweep to a `SNSweepOperator`.

Refactor `_sweep_1d_spherical` and `_sweep_1d_cylindrical` into the
methods of `SNSweepOperator(LinearOperatorMixin)`. The free
function `transport_sweep` becomes a thin compatibility wrapper
that builds `SNSweepOperator(sn_mesh, sig_t).solve(Q, psi_bc,
Q_aniso)`. Capabilities: `{apply, solve}` where `apply(ψ)` is
"apply L to ψ via WDD closure" and `solve(q)` is "apply L⁻¹ to q
via the sweep".

**Acceptance criterion**: 11 regression snapshots
(`tests/sn/regression/snapshots/`) PASS at `assert_allclose
rtol=1e-12`. Bit-equivalence on flat-ψ is structural; numerical
equivalence on stored snapshots is via the wrapper.
**Tests**: `pytest tests/sn/regression/ tests/sn/spatial/`.

### Step 4 (medium, ~5 days). Replace the inlined matvec WDD recurrence with composition.

Refactor `transport_operator_matvec_spherical` and `_cylindrical` to
delegate the per-cell WDD recurrence to the **same**
`SNCellOperator.apply` that `SNSweepOperator` uses internally. The
matvec becomes:

```python
def transport_operator_matvec_spherical(...):
    L = build_streaming_operator_spherical(sn_mesh, sig_t, bc_outer, pole_closure)
    psi_structured = packed_to_structured.apply(solution)
    L_psi_structured = L.apply(psi_structured)
    return structured_to_packed.apply(L_psi_structured)
```

The WDD recurrence in operator.py:792 / 825 disappears — it now lives
ONCE inside `SNCellOperator`. Manifestation #7 dissolves at this step
IF the two paths converge to the same fixed point under the chosen
closure. Empirical validation: `tests/sn/spatial/test_sweep_vs_apply_consistency.py`
(57 tests) should now exhibit bit-identity, not just flat-flux
equivalence, on representative non-flat ψ.

**Acceptance criterion**: (a) `test_sweep_vs_apply_consistency.py`
all 57 GREEN at machine precision on non-flat ψ;
(b) `test_phase_e_trajectory_resolvent_flux_shape_crosscheck` xfail
flips to xpass and the strict marker can be removed.
**Tests**: full SN test suite + the Phase E sentinel.

### Step 5 (medium, ~5 days). Replace ad-hoc layout adapters with `PackedStructuredAdapter`.

The two layout-conversion code paths (`solution_to_angular_flux_*` in
operator.py; the inverse in `_make_sweep_preconditioner`) become
methods of a single `PackedStructuredAdapter(LinearOperator)`
parameterised on the `EquationMap`. The BC fill becomes the adapter's
`__init__`-bound `bc_*` operators applied uniformly.

**Acceptance criterion**: bit-identical regression. **Tests**:
`pytest tests/sn/regression/`.

### Step 6 (medium, ~5 days). Wire the iteration primitives (`SourceIteration`, `PreconditionedGMRES`).

`_solve_source_iteration`, `_solve_krylov`, `_solve_fixed_source_si`,
`_solve_fixed_source_krylov` collapse to the 3-5 line versions in
§5.2-5.4. The primitives live in `orpheus.numerics.iteration` (per
the operator_algebra.rst Phase 0 stub at line 226-230).

**Acceptance criterion**: full eigenvalue + fixed-source MMS suite
green. **Tests**: full `pytest tests/sn/`.

### Step 7 (medium, ~5 days). Choose and document the closure convention.

With the operator-algebra unification complete, the
symmetric-vs-WDD-asymmetric closure choice is now a parameter on
`SNCellOperator`. Document the trade-off explicitly in
`docs/theory/discrete_ordinates.rst`; choose ONE closure for the
default (likely WDD given the asymptotic-diffusion-limit
properties); make the other available via the strategy. Manifestation
#7 either dissolves entirely (one closure for both paths) or becomes
a principled documented choice.

**Acceptance criterion**: Phase E xfail removed; ERR-026 closeable to
fully CLOSED (along with Issue #195 pre-asymptotic magnitude
investigation).

### Total: 7 atomic steps, ~25 days estimated effort.

Each step preserves the regression contract; each step is
independently committable; each step extends the Wave H primitives
rather than demolishing them.

---

## Section 7 — Risks and decision points

### 7.1 Closure choice — WDD vs symmetric

**Decision needed at Step 4**. The current `SNStreamingOperator.apply`
uses symmetric-closure; the current `SNStreamingOperator.solve` uses
WDD. Choosing one as the default for `SNCellOperator` decides
Manifestation #7.

**Data to resolve**: comparison of L1 MMS convergence rates for both
closures across (slab, sphere, cylinder) at fixed nx; comparison of
the asymptotic-diffusion-limit accuracy (Bailey-Morel-Chang 2010 is
the canonical reference, available locally per
`scratch/literature/`).

**Provisional recommendation**: WDD asymmetric. Reasons: (a) it
yields the formal asymptotic-diffusion-limit accuracy that
Bailey-Morel-Chang prove; (b) it's the production sweep's existing
closure, so the 11 frozen snapshots are stable; (c) the symmetric
closure was historically only chosen for the BiCGSTAB FD operator's
boundary-face approximation, which is now replaced by face-trace BC
algebra. The Wave E "Krylov-on-apply with WDD closure" reconciliation
recommends WDD as the canonical choice.

### 7.2 Sweep direction state machine

The current SI sweep mutates `psi_bc["bc_sph"]` (sweep.py:469-471)
to carry per-iteration outflow into the next iteration's reflective
BC. Under the unified architecture this state belongs to either (a)
`SNSweepOperator` as instance state, OR (b) the user-supplied
`psi_bc` dict as before.

**Decision**: option (b) — keep state external, makes
`SNSweepOperator` stateless and re-usable. Required for
`_make_sweep_preconditioner` to remain stateless inside GMRES.

### 7.3 Pole-face anchor convention

The current matvec anchors `psi_face_in = fi[:, outgoing, i0, 0]`
(cell-centre value, operator.py:784); the current sweep anchors
`psi_spatial_in = zeros(ng)` (sweep.py:559). On flat ψ both work
(A[0]·anchor = 0). On non-flat ψ they differ — the matvec's choice
is the "analytical extension" Phase C added; the sweep's choice
matches the BC-spec "no flux flows out of the pole".

**Decision needed**: unify on the analytical-extension form (Phase C
chose this), but document the convention in `PoleFaceAnchor` so the
choice is explicit and tested at the type level. Alternative
positions for the convention live as named subclasses of
`PoleFaceAnchor`.

### 7.4 Backward-compatibility hooks

The transitional Wave H windows demand that legacy callers continue
to work. Strategy: keep the free functions (`transport_sweep`,
`transport_operator_matvec_*`) as thin wrappers around the new
operator types throughout the migration. Each migration step
deprecates one wrapper at a time; the final cleanup wave removes
them. The 11 frozen regression snapshots gate every step.

### 7.5 Where the architecture might still have multiple valid forms

(a) Should `SNStreamingOperator` survive as a top-level type, or
should it dissolve into `SNSweepOperator` (sweep-direction-symmetric)
plus `_curvilinear_streaming_op_apply` (a function returning a
composition expression `pack ∘ Σ_dir SweepSlice ∘ unpack`)? Both
forms are mathematically equivalent. The function form reads more
like math; the type form integrates more cleanly with the wider
operator algebra. **Decision deferred to Step 4 design.**

(b) Should the per-cell `CellUpdate` strategy become a
`LinearOperator` directly (option A) or be wrapped by an
adapter (option B)? Option A is cleaner; option B preserves the
existing `CellUpdate` Protocol untouched. **Recommendation**: option
A, with `CellUpdate` Protocol updated to declare it inherits from
`LinearOperator`.

---

## Pointers

- Plan: `.claude/plans/sn_reshape.md` (Wave H Phase A-F closeouts)
- Phase F closeout: `.claude/agent-memory/method-implementer/issue_168_phase_f_closeout.md`
- ERR-026 catalogue: `.claude/skills/vv-principles/error_catalog.md:1602+`
- Operator algebra: `docs/theory/operator_algebra.rst` (Phase 0 stub)
- BC realisation: `orpheus/sn/boundary_realizer.py`
- Wave H plan: `.claude/plans/sn_reshape.md:153-180` (Wave H rows)
- Issue #196 + Manifestation #7: GitHub
- Issue #168 + Phase F: GitHub
- Test consistency suite: `tests/sn/spatial/test_sweep_vs_apply_consistency.py`
- Phase E sentinel: `tests/sn/test_phase_c_crosscheck.py::test_phase_e_trajectory_resolvent_flux_shape_crosscheck`

## Linked memories

- [[issue-168-phase-f-closeout]] — the twin-bug that triggered this audit
- [[issue-168-phase-d-step3-closeout]] — Carlson seed on apply path
- [[c188-curvilinear-realizer-unify]] — BC-realiser unification across slab + curvilinear

## Headline summary

The operator algebra is correctly designed in `orpheus.numerics.operator`
and `orpheus.sn.boundary_realizer`; it is correctly wired into
`SNStreamingOperator` as a Protocol-conformant type. **But the
load-bearing transport paths (SI sweep, Krylov matvec, fixed-source
dispatch) do not COMPOSE with it.** They re-implement, in procedural
Python, the same five operator-algebra primitives (cell-update, sweep
composition, angular redistribution, angular integration, BC face-trace
application) twice, on different vector layouts, with different BC
entry points. Manifestation #7 is the load-bearing visible symptom of
this duplication. The unification target is "every per-cell, per-sweep,
and per-iteration step is the application of a named operator
composition" — and the unification is achievable in 7 atomic
independently-committable migration steps that extend the Wave H
primitives rather than demolishing them.
