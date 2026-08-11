# #340 N6b — what is in scope at the five `_warn_if_unconverged` sites

**Tree:** branch `refactor/operator-strategy-layers`, HEAD `0d38d331`.
`orpheus/` and `tests/` were CLEAN at open. ⚠ **Intra-dispatch drift:**
`orpheus/geometry/mesh.py` went clean → `+35/−5` while this ran (a `BC.vacuum` /
`BC.reflective` / `BC.white` `ClassVar` declaration for autodoc, #346 W1
follow-up). It touches nothing in this map — no solver, no `volumes`, no
`volume_measure` — so every finding below stands. Re-verified at close.

**Method.** Every claim about scope is MEASURED, not read: a spy wrapping
`orpheus.sn.solver._warn_if_unconverged` (a module-global, so patching the
module attribute takes at all five sites) dumps `sys._getframe(1).f_locals`
at the moment of the call, then the residual is actually computed from those
names. Probes in `/Users/rodrigo/.claude/jobs/c30e4f25/tmp/probe_{scope,compute,cost,edge,lag,windowed,ld}.py`.

**Tree-wide census:** exactly **five** production call sites (grep
`_warn_if_unconverged` over `orpheus/`), plus 8 direct-call test probes in
`tests/sn/solve/test_convergence_contract.py`.

---

## 0. Headline — three findings that change the N6b design

1. ⛔ **The plan's `[M]` "`solve_sn` discards the solver it builds" is FALSE.**
   `solver` is bound at `orpheus/sn/solver.py:2364` and is **live at the
   warning call** (`solver.compute_fission_source` at 2389, `solver.mat_xs` /
   `solver.scattering_op` at 2422). Measured: 38 locals at the call, `solver :
   orpheus.sn.solver.SNSolver` among them, on BOTH inner drivers. **Nothing has
   to be threaded onto `IterationHistory` for scope reasons.** The plan's
   premise was about a *post-hoc* computation from a returned `Solution`; at
   the call site the operators are all still local.

2. ✅ **All five sites can compute the residual — measured, all five.** Each
   one succeeded and produced a finite balance projection. See §1.

3. ⛔ **Two paths where the number is NOT available**, and they are different
   in kind:
   - **LD / moment-tailed schemes — NOT COMPUTABLE AT ALL.** `evaluate_residual`
     raises `ValueError: AngularResidual: values.shape (8, 2, 40, 2) does not
     match space.shape (8, 2, 40)`. This is the un-built residual-mint widening
     already named in `_certify_within_group_exit`'s docstring
     (`solver.py:613-622`). Reach is limited: `solve_sn` and `solve_sn_adjoint`
     have **no `scheme=` parameter** (measured) and `_as_sn_mesh` always builds a
     fresh `SNMesh`, so LD is structurally unreachable from both eigenvalue
     entries; `solve_sn_adjoint_fixed_source` **raises earlier** on LD
     (`AngularSourceSink: values.shape (8,2,40,2)…`), so its warning site is
     never reached. ⟹ **the hole is exactly the two `solve_sn_fixed_source`
     arms (3622, 3844) under `scheme=LinearDiscontinuous()`.**
   - **The windowed 2-D SI arm — computable, but NOT at line 3622 as ordered.**
     Measured: at 3622 `psi_typed.interior` is a `HarmonicMomentFlux` of shape
     `(1, 1, 2, 8, 8)`; `evaluate_residual` raises `IndexError: index 13 is out
     of bounds for axis 0 with size 1`. The full-angular reconstruction that
     fixes it (`angular_out`, `solver.py:3650`) happens **28 lines AFTER** the
     warning — and `'angular_out' in f_locals` measured `False` at the call.
     Moving the warning below the reconstruction makes it computable. This is
     an *ordering* defect, not a capability gap.

   Contrast: **`solve_sn`'s exit is windowing-immune.** Measured on a 2-D
   Lebedev(5) windowed eigenvalue solve, `final_psi_a` is a full-angular
   `TimedFullField`/`AngularFlux` `(14, 2, 8, 8)` and `evaluate_residual`
   succeeds — because the reconstruction sweep already ran (2490/2497).

---

## 1. The five-row table

All names verified live in `f_locals` at the call; `file:line` = the **binding**
line in `orpheus/sn/solver.py`.

### Row 1 — `solve_sn` @ `orpheus/sn/solver.py:2518` (`where="solve_sn"`)

| Q | answer |
|---|---|
| **1. live `SNSolver`?** | ✅ **YES.** `2364: solver = SNSolver(` — used at 2389/2390/2391/2422; identical on `inner_solver="krylov"` (measured, the exit block 2378–2530 is driver-independent). |
| **2. system / loss / grid?** | ⚠ **Half.** `2421: final_implicit = build_within_group_system(` … `.implicit_operator` — the `WithinGroupSystem` IS constructed at the exit but **only `.implicit_operator` (= `M = L+C`) is kept**; the record, and hence `.loss` / `_bare_loss_arm`, is discarded on the same expression. Binding it costs **zero extra work**. Rebuilding it costs `[M] 0.051–0.067 ms`. |
| **3. typed final iterate?** | ✅ `2497: final_psi_a = final_implicit.solve(final_rhs_a)` → `TimedFullField`; on a carrying mesh `2494: final_psi_a = _system_a_member(final_state)` **and** `2490: final_state = final_implicit.solve(…)` → `CoupledField` (measured on a sphere: `final_state: CoupledField`, `final_ray: RadialCharacteristicField`, `corner_state: RadialCharacteristicField` all in scope). ⚠ On a carrying mesh you MUST pass `final_state` + `system.loss`: measured, the bare System-A call is refused — `ValueError: evaluate_residual: this mesh carries starting-direction levels (R12a) …`. |
| **4. rhs?** | ⚠ **Not directly — and the obvious candidate is a trap.** `2420: q_final_per_ord = AngularSourceSink.from_isotropic(Q_final, sn_mesh)` wraps `Q_final`, which is the **TOTAL** source: `2389: Q_final = solver.compute_fission_source(scalar_flux, keff)` then `2390: solver._add_scattering_source(Q_final, scalar_flux)` and `2391: solver._add_n2n_source(Q_final, scalar_flux)` — both **mutate in place**. `[M] ‖Q_final − Fφ/k‖/‖Fφ/k‖ = 1.37e+01`. Using it against `A = L+C−S−B` would double-count scattering. **Reconstruct:** `solver.compute_fission_source(φ, keff)` (`[M] 0.021 ms`) + the `from_isotropic` + `TimedFullField` wrap (`[M] 0.026 ms`). Both `keff` (`2374`) and `scalar_flux` (`2374`) are local, and `AngularSourceSink` / `AngularBoundarySourceSink` are **already imported into the function's locals** (2385, 2451). On a carrying mesh the coupled rhs additionally needs `_radial_characteristic_fission_seed(q_scalar, sn_mesh)` + `_coupled_source_state(…)` — both module-level, measured working. |
| **5. `sn_mesh`?** | ✅ `2362: sn_mesh = _as_sn_mesh(mesh, quadrature, materials, mat_map=mat_map)` |

**Verdict: fully computable.** Measured end-to-end at this site, both rhs
variants: `defect = 7.73e-02`, `balance = 7.75e-02` with `φ` from the returned ψ
(the N5 headline convention); `7.49e-02` / `7.41e-02` with `φ = scalar_flux`.

### Row 2 — `solve_sn_adjoint` @ `orpheus/sn/solver.py:2794`

| Q | answer |
|---|---|
| **1. live `SNSolver`?** | ⛔ **NO.** This entry never builds one. Its stand-in is `2755: implicit_operator, gain, F_posed, template = _adjoint_posing_parts(` — which internally does `mat_xs = sn_mesh.material_xs_field()` (`2634`), i.e. the same XS source `SNSolver.__init__` uses at `1189`. Nothing is missing; the shape is just operators-not-solver. |
| **2. system / loss?** | ⚠ **Reconstructible in one expression, no rebuild.** `_adjoint_posing_parts` returns the *parts*, not the record — but `build_within_group_system` builds `A_AA = LC − S − B_a` (`coupled_system.py:509`) and returns `explicit_gains=(S, B_a)` (`532`), while the adjoint entry folds `gain = S + B_a` (`2638-2640`). So `implicit_operator − gain` **is** the loss by construction, and the daggered loss is `implicit_operator.H − gain.H` (`LinearOperator.__sub__` at `numerics/operator.py:779`). Measured working. |
| **3. typed final iterate?** | ✅ `2773: k_adj, keff_history, psi_star = (outcome.keff, …, outcome.flux_distribution)` → **`FullField`** (measured), not `TimedFullField`. `evaluate_residual` accepts it (measured OK) even though its annotation reads `TimedFullField \| CoupledField` (`solver.py:243`) — ⚠ the annotation is narrower than the behaviour; widen it or the site needs a cast. |
| **4. rhs?** | ✅ **Cheapest of all five.** `F_posed` is posed on the FULL-FIELD space (`2641-2643`), so `F_posed.H.apply(psi_star)` returns a `FullField` **directly** — no scalar detour, no `from_isotropic` wrap. `[M] 0.190 ms`. `k_adj` (`2773`) is local; the rhs is `F_posed.H.apply(psi_star) * (1/k_adj)`. |
| **5. `sn_mesh`?** | ✅ `2754: sn_mesh = _as_sn_mesh(mesh, quadrature, materials, mat_map=mat_map)` |

**Verdict: fully computable.** Measured: `defect = 9.19e-06`, `balance = 1.66e-06`.
⚠ *Interpretive* caveat (not a scope caveat): the balance projection's meaning —
"the component of the residual the functional `k = production/absorption` can
see" — was measured by N5 on the **forward** eigenvalue only. On the adjoint the
same arithmetic is a formally-valid rate defect but its overlap number has never
been measured.

### Row 3 — `solve_sn_adjoint_fixed_source` @ `orpheus/sn/solver.py:2940`

| Q | answer |
|---|---|
| **1. live `SNSolver`?** | ⛔ **NO** — same posing-parts shape as Row 2. |
| **2. system / loss?** | ⚠ Same as Row 2: `2888: implicit_operator, gain, _F, template = _adjoint_posing_parts(` ⟹ `implicit_operator.H − gain.H`. Carrying meshes are refused up front at `2880-2887`, so the coupled arm never arises here. |
| **3. typed final iterate?** | ✅ `2938: psi_star, record = si.solve(q_star, initial_guess=template)` → `FullField`. |
| **4. rhs?** | ✅ **Already there, nothing to reconstruct.** `2899: q_star = detector_response` (composite branch) or `2924: q_star = FullField(` (the angle-flat scalar lift). It is a genuine given — no lag, no reconstruction. |
| **5. `sn_mesh`?** | ✅ `2876: sn_mesh = _as_sn_mesh(` |

**Verdict: fully computable, zero reconstruction.** Measured at `max_inner=5`:
`defect = 5.80e-01`, `balance = 5.46e-01`.

### Row 4 — `_solve_fixed_source_si` @ `orpheus/sn/solver.py:3622` (`where="solve_sn_fixed_source"`)

| Q | answer |
|---|---|
| **1. live `SNSolver`?** | ✅ `3481: solver: SNSolver,` (parameter). |
| **2. system / loss?** | ✅ **Best-provisioned of all five.** `3557: system = build_within_group_system(` — the full record is bound, plus `3565: coupled = isinstance(system.implicit_operator, CoupledOperator)`. The loss is the certificate's own expression, verbatim: `system.loss if coupled else _bare_loss_arm(system)` (`3600`). Also `3560: si, base_implicit, gains, windowed = _within_group_si(`. |
| **3. typed final iterate?** | ⚠ `3592: psi_typed, record = si.solve(` and `3607: psi_full = _system_a_member(psi_typed)` — full-angular `TimedFullField` when `windowed=False`; **a `HarmonicMomentFlux`-bulk moment iterate when `windowed=True`** (measured `(1,1,2,8,8)`), on which `evaluate_residual` raises. The full-angular `angular_out` is bound at `3650`/`3656`, i.e. **after** this call (`'angular_out' in f_locals` measured `False`). |
| **4. rhs?** | ✅ `3483: q_ext_composite: "TimedFullField \| CoupledField",` (parameter) — the real external rhs, built once by `_build_fixed_source_rhs`. No lag, no reconstruction. |
| **5. `sn_mesh`?** | ✅ `3482: sn_mesh: SNMesh,` (parameter). |

**Verdict: computable when `windowed=False` (measured `defect = 8.66e-01`,
`balance = 8.14e-01` at `max_inner=5`); NOT computable at this line when
`windowed=True`, fixable by ordering; NOT computable at all under LD.**

### Row 5 — `_solve_fixed_source_krylov` @ `orpheus/sn/solver.py:3844` (`where="solve_sn_fixed_source"`)

| Q | answer |
|---|---|
| **1. live `SNSolver`?** | ✅ `3686: solver: SNSolver,` (parameter). |
| **2. system / loss?** | ✅ `3782: system = build_within_group_system(`, `3785: coupled = …`. Certificate expression verbatim at `3804`. |
| **3. typed final iterate?** | ✅ `3797: psi_typed, record = krylov.solve(`, `3811: psi_full = _system_a_member(psi_typed)`. **Always full-angular** — windowing is SI-only (`solver.py:3800-3802`). Also `3822: bulk = psi_full.interior` (`AngularFlux`) and `3829: phi = _average_moment_scalar(…)`. |
| **4. rhs?** | ✅ `3688: q_ext_composite:` (parameter). |
| **5. `sn_mesh`?** | ✅ `3687: sn_mesh: SNMesh,` (parameter). |

**Verdict: fully computable, zero reconstruction, no windowed caveat.**
Measured `defect = 2.81e-13`, `balance = 8.75e-14` at `max_inner=5`.

### Summary matrix

| # | site | `SNSolver` | loss reachable | typed ψ | rhs | `sn_mesh` | computable? |
|---|---|---|---|---|---|---|---|
| 1 | `solve_sn` 2518 | ✅ 2364 | ⚠ record dropped at 2421 (rebind = free) | ✅ 2494/2497 | ⚠ rebuild `Fφ/k` (0.05 ms) | ✅ 2362 | ✅ |
| 2 | `solve_sn_adjoint` 2794 | ⛔ none | ⚠ `impl.H − gain.H` (2755) | ✅ 2773 (`FullField`) | ✅ `F_posed.H.apply(ψ*)/k` | ✅ 2754 | ✅ |
| 3 | `…adjoint_fixed_source` 2940 | ⛔ none | ⚠ `impl.H − gain.H` (2888) | ✅ 2938 (`FullField`) | ✅ `q_star` 2899/2924 | ✅ 2876 | ✅ |
| 4 | `…fixed_source_si` 3622 | ✅ 3481 | ✅ `system` 3557 | ⚠ moment iterate if windowed | ✅ `q_ext_composite` 3483 | ✅ 3482 | ⚠ order / ⛔ LD |
| 5 | `…fixed_source_krylov` 3844 | ✅ 3686 | ✅ `system` 3782 | ✅ 3797 | ✅ `q_ext_composite` 3688 | ✅ 3687 | ⚠ ⛔ LD only |

⚠ Note that rows 1–2 and rows 3–5 answer **different equations**: rows 1–2 are
the OUTER eigenvalue residual `Aψ − Fφ(ψ)/k` (the equation N5 measured); rows
3–5 are the fixed-source residual `Aψ − q_ext`, where `q_ext` is a genuine
given and there is no lag at all. Both are honest exit residuals for their
entry; only rows 1–2 inherit N5's 4.64× overlap figure.

---

## 2. Certificate reuse — VERDICT: **NO. They are not the same equation.**

### 2a. `q_driver` at the two eigenvalue sites is the **LAGGED** fission source

The chain, quoted:

`orpheus/numerics/eigenvalue.py:404-405` — the outer loop:
```python
        fission_source = solver.compute_fission_source(flux_distribution, keff)
        flux_distribution = solver.solve_fixed_source(fission_source, flux_distribution)
```
`flux_distribution` here is the **previous** outer's iterate (stashed as
`flux_old` at 397) and `keff` the **previous** outer's eigenvalue
(`keff_old = keff` at 402, `keff` not updated until 432).

`orpheus/sn/solver.py:1728-1729` — the SI inner receives it as a parameter:
```python
    def _solve_source_iteration(
        self, fission_source: np.ndarray, flux_distribution: np.ndarray,
```

`orpheus/sn/solver.py:1805-1807` — wrapped verbatim, no recomputation:
```python
        q_ext_per_ord = AngularSourceSink.from_isotropic(
            fission_source, self.sn_mesh,
        )
```

`orpheus/sn/solver.py:1884-1894` — `q_driver` is that composite (plus, on a
carrying mesh, the ψ½ fold of **the same** `fission_source`):
```python
        q_driver = (
            _coupled_source_state(
                q_ext_composite,
                _radial_characteristic_fission_seed(
                    fission_source, self.sn_mesh,
                ),
                self.sn_mesh,
                context="SNSolver._solve_source_iteration",
            )
            if coupled else q_ext_composite
        )
```
The Krylov arm is byte-for-byte the same construction at `2070-2080`.

⟹ **LAGGED.** `q_driver = F φ^(n−1) / k^(n−1)`; the ψ it is paired with is
ψ^(n). Nothing recomputes it from the current ψ.

### 2b. Measured: how far apart the two equations are

Same fixture (§3), both defects evaluated on the **same** ψ at every outer:

| driver | outer | `record.converged` | defect vs `q_driver` (LAGGED) | defect vs `Fφ(ψ)/k` (OUTER) | `‖q_lag − q_now‖/‖q_now‖` |
|---|---|---|---|---|---|
| SI, `max_outer=3` | 1 | False | 1.361e-01 | 3.099e-01 | 2.583e-01 |
| | 2 | False | 1.063e-01 | 2.289e-01 | 1.908e-01 |
| | 3 | False | 2.037e-02 | 3.807e-02 | 3.174e-02 |
| SI, `max_outer=10` | 10 | False | 7.778e-08 | 1.409e-07 | 1.174e-07 |
| **Krylov, `max_outer=3`** | 1 | **True** | **9.457e-09** | **9.214e-03** | 9.214e-03 |
| | 2 | True | 6.838e-09 | 1.179e-04 | 1.179e-04 |
| | 3 | True | 9.396e-09 | 1.628e-06 | 1.628e-06 |

The Krylov row is decisive: the certificate's residual is **~10⁶× smaller** than
the outer residual at outer #1 and still **173×** smaller at the truncation
point. It measures *"did the inner solve the equation it was handed?"* — which,
by construction, a converged inner answers "yes" to. It carries essentially
**zero** information about the outer truncation the warning is about.

Two further structural reasons reuse fails, independent of the numbers:

- **It usually does not evaluate at all.** `_certify_within_group_exit` returns
  early on `not record.converged` (`solver.py:626-627`) — and on the SI arm with
  `max_inner=40, inner_tol=1e-8` the inner record measured `converged=False` at
  **every** outer, so the certificate computed nothing on the whole solve. Which
  branch you get depends on the inner driver, not on the outer's health.
- **It is not evaluated on the returned ψ.** `solve_sn` performs a **final
  full-angular reconstruction sweep** (`2490`/`2497`) after the loop; the object
  the user receives is `final_psi_a`, never the `psi_typed` the certificate saw.

### 2c. How many times each site is reached — **once per OUTER iteration**

Structurally: `_certify_within_group_exit` sits inside
`SNSolver._solve_source_iteration` (1902) / `_solve_krylov` (2089), which
`power_iteration` calls once per outer via `solver.solve_fixed_source`
(`eigenvalue.py:405` → `solver.py:1373-1376`). Measured hit counts:

| solve | hits |
|---|---|
| `solve_sn[source_iteration](max_outer=3)` | `{('_solve_source_iteration', 1902): 3}` |
| `solve_sn[source_iteration](max_outer=10)` | `{('_solve_source_iteration', 1902): 10}` |
| `solve_sn[krylov](max_outer=3)` | `{('_solve_krylov', 2089): 3}` |
| `_solve_fixed_source_si` (one fixed-source solve) | `{(…, 3599): 1}` |
| `_solve_fixed_source_krylov` | `{(…, 3803): 1}` |

⭐ **Two sites reach the certificate ZERO times, ever:** `solve_sn_adjoint` and
`solve_sn_adjoint_fixed_source` measured `{}`. The adjoint eigenvalue path runs
`KEigenvalue.solve → power_iteration(KEigenvalue) → SourceIteration.solve`
entirely inside `numerics/` and never enters SN solver code (the L-024 nesting
asymmetry); the adjoint fixed-source entry builds its own `SourceIteration`
inline at `2931-2934`. So on the two adjoint entries there is no certificate
residual to reuse even in principle.

⟹ **The projection needs its own evaluation at the top-level exit.** Do not
substitute the certificate's number; if you ever want to, the 4.64× overlap
must be re-measured for the lagged residual and quoted as a different figure.

---

## 3. Cost — one extra forward apply is **≈1 %** of a truncated eigenvalue solve

### The fixture, completely stated

- **Geometry:** 1-D slab, `StructuredGeometry(geometry="SLB")`, two regions —
  `Region(mat_id=2, outer_thickness_cm=1.0)` + `Region(mat_id=0,
  outer_thickness_cm=1.0)`; `bcs=(BC.reflective, BC.reflective)`.
- **Mesh:** `Mesh1D.from_geometry(..., region_meshes=(RegionMesh(n_cells=20),
  RegionMesh(n_cells=20)))` ⟹ **40 cells**, uniform `Δ = 0.05 cm`.
- **Materials:** `{0: get_mixture("B", "2g"), 2: get_mixture("A", "2g")}` from
  `orpheus.derivations.common.xs_library` ⟹ **ng = 2**.
- **Quadrature:** `Quadrature.gauss_legendre(8)` ⟹ **N = 8**.
- **Scheme:** default `DiamondDifference` (`spatial_basis_per_axis = 1`).
- **Solver:** `solve_sn(..., inner_solver="source_iteration")` (default),
  `inner_schedule="jacobi"` (default), `inner_tol=1e-8` (default),
  `keff_tol=1e-7`, `flux_tol=1e-6` (defaults).
- **Machine/env:** host `.venv`, Python 3.14, macOS (Darwin 25.4.0), serial.
  Timings are **medians of 20 reps after 3 warm-up reps**; solve times are
  medians of 5 reps after 1 warm-up. Measured 2026-08-10.

### Numbers

| budget | total inner iters | solve (median) | `evaluate_residual` | whole diagnostic | fraction |
|---|---|---|---|---|---|
| `max_outer=3, max_inner=40` | 120 | 92.58 ms | 2.221 ms | **2.371 ms** | **2.56 %** |
| `max_outer=3, max_inner=10` | 30 | 30.22 ms | 2.081 ms | **2.229 ms** | **7.38 %** |
| `max_outer=10, max_inner=40` | 400 | 302.07 ms | 2.052 ms | **2.184 ms** | **0.72 %** |

"Whole diagnostic" = `evaluate_residual` + `compute_fission_source`
(`[M] 0.021 ms`) + the `TimedFullField`/`from_isotropic` wrap (`[M] 0.022 ms`)
+ `build_within_group_system` (`[M] 0.058 ms`, and **free** if the record is
bound at 2421 instead of dropped) + the balance projection ×2
(`[M] 0.031 ms` for both).

**The honest scaling law, which is what makes the number reusable:** the
diagnostic costs **≈ 3 inner iterations** (2.2 ms ÷ 0.77 ms/inner at the 120-inner
point; 2.2 ms ÷ 1.0 ms/inner at 30). So

> fraction ≈ **3 / (total inner iterations)**

⟹ **order of magnitude: ~1 %.** It is ~10 % only on a pathologically truncated
solve (≲30 inner iterations total) and ~0.1 % on anything with ≳3000.

For contrast, on the same fixture: one bare `(L+C).solve` sweep = `[M] 0.176 ms`;
the prototype's from-scratch `SNSolver(sn_mesh)` construction = `[M] 2.591 ms`
(i.e. **the rebuild the prototype does costs more than the residual itself** —
and it is entirely avoidable at all five sites).

Balance values measured alongside, for reference:
`max_outer=3/max_inner=40` → `balance = 7.75e-02`, `R_g = [0.0377, 0.0510]`,
`Q_g = [0.8188, 0.0]`; `max_outer=3/max_inner=10` → `balance = 1.396`;
`max_outer=10/max_inner=40` (outer converged, all inners starved) →
`balance = 2.96e-07`.

Carrying-mesh cost (sphere, `gauss_legendre(8)`, 20 cells, 2 g, vacuum outer):
the coupled `evaluate_residual(system.loss, final_state, q_coupled)` =
`[M] 2.30 ms` — same order.

---

## 4. ⭐ What I would do

**Compute it at each site, pass it as one new keyword to
`_warn_if_unconverged`; do not store it on `IterationHistory` for scope
reasons — there are none.** The plan's blocker ("the projection cannot be
computed at the entry from public state") is true of a *post-hoc* computation
from a returned `Solution` and false at the call site: measured, all five
frames carry the mesh, the operators (or a one-expression reconstruction of
them), the typed iterate, and either the rhs or a `[M] 0.05 ms` reconstruction
of it. A field on `IterationHistory` may still be the right *product* decision
(so a caller can read the number without parsing a warning string), but it is
now a design choice, not a forced move — and that inverts the §N6b argument
that led to it.

The cheapest honest shape is **one SN-local helper** —
`_exit_balance_projection(loss, psi, q, sn_mesh) -> tuple[np.ndarray, np.ndarray] | None`
returning `(R_g, Q_g)` — called once per solve at each site with that site's
own triple, and `None` on the LD refusal. Per site, the threading cost:

- **3844 (Krylov fixed-source) and 3622 (SI fixed-source)** — *zero threading*.
  The triple is literally the certificate's own argument list, already on the
  line above. Start here.
- **2940 (adjoint fixed-source)** — *zero reconstruction*: `implicit_operator.H
  − gain.H`, `psi_star`, `q_star`.
- **2794 (adjoint eigenvalue)** — one operator apply: `F_posed.H.apply(psi_star)
  * (1/k_adj)`, `[M] 0.19 ms`.
- **2518 (`solve_sn`)** — the only one needing an edit above the call, and it is
  a *free* one: split `2421` into `final_system = build_within_group_system(…)`
  / `final_implicit = final_system.implicit_operator`. The system is already
  built there; only the record is thrown away. Then recompute the fission rhs
  from `psi.interior.integrate_angular()` (`[M] 0.05 ms`) — recompute rather
  than snapshot `Q_final` before line 2390, because recomputing from the
  **returned** ψ is exactly the convention `scratch/n5_outer_cert_lib.py`
  measured the 4.64× against, and a snapshot would silently give the
  `defect_pi` variant instead.

**Two things must be decided before the field is designed, not after:**

1. **⛔ The number is unavailable on `solve_sn_fixed_source(...,
   scheme=LinearDiscontinuous())` — both arms — and no site-local change fixes
   it.** `AngularResidual.from_balance` refuses the trailing `2^d` spatial-moment
   axis; that widening is the un-built residual carve already tracked in
   `_certify_within_group_exit`'s docstring (#310 deferred-out list). So the
   message WILL carry the number on some paths and not others. Decide the
   spelling now: I would have the helper return `None` and the message say
   nothing about balance on that path (a silent omission), **not** print
   "unavailable" — but the choice must be explicit, and whichever way it goes,
   the LD row needs an `xfail(strict=True)`-style pin so the day the mint widens,
   the omission surfaces instead of persisting.

2. **⚠ Site 3622 needs a one-block reorder, or it joins the LD hole for a
   reason that is purely accidental.** Move `history = IterationHistory(record=record)`
   + `_warn_if_unconverged(...)` (currently `3621-3624`) down to just above the
   `return Solution(...)` at `3675`, so `angular_out` (`3650`/`3656`) is bound.
   `history` is consumed at `3680` anyway, so the move is local; the only
   behavioural change is that a truncated *windowed* solve now hears the warning
   after one extra reconstruction sweep. Without the move, the windowed 2-D SI
   arm silently loses the number — and it is the one arm where a reader would
   least suspect it, because the *sibling* Krylov arm and the *eigenvalue* entry
   both work fine on the same 2-D mesh.

**Refuted candidates, with the structural reason** (so nobody re-attacks them):

- *Reuse `_certify_within_group_exit`'s residual.* ⛔ Different equation
  (lagged rhs — §2a), measured up to 10⁶× apart (§2b), evaluated once per outer
  rather than once per solve (§2c), skipped entirely whenever the inner did not
  claim convergence, never evaluated on the reconstructed ψ the user receives,
  and never reached at all on the two adjoint entries.
- *Reuse `q_final_per_ord` / `Q_final` at 2518 as the rhs.* ⛔ It is the TOTAL
  reconstruction source (fission + P0 scatter + (n,2n)); `[M] 13.7×` larger than
  the fission source. Against `A = L+C−S−B` it double-counts scattering.
- *Compute post-hoc from the returned `Solution` (the prototype's shape).* ⛔
  Requires a fresh `SNSolver(sn_mesh)` — `[M] 2.59 ms`, i.e. more than the
  residual it enables — and re-derives state that is already local. It was the
  right shape for a measurement harness and is the wrong one for production.
