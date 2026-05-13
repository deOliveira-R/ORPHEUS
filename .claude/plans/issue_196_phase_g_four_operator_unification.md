# Phase G — Four-operator unification (L, C, S, F) under one algebra

**Tracking issue**: Continuation of Issue [#196](https://github.com/deOliveira-R/ORPHEUS/issues/196) (ERR-026 manifestation #7 — residual SI-vs-Krylov WDD asymmetry). Phase F closed the structural pole-cell defect; Phase G closes the architectural duplication that produced the defect by routing **all** SN load-bearing paths through one operator algebra.
**Branch**: `refactor/sn-operator-algebra` (continue from Phase F tip `b0cc1b1`)
**Predecessor closeouts + audits**:
- `.claude/agent-memory/method-implementer/issue_168_phase_f_closeout.md` (Phase F tactical fix)
- `.claude/agent-memory/explorer/issue_196_sn_operator_architecture_audit.md` (815-line architectural audit; **7 proposed operators reduced to 4 + 1 by reconciliation**)
- `.claude/agent-memory/general-purpose/phase_g_four_operator_architecture_reconciliation.md` (reconciliation with `neutron_transport_grand_report_v3.md`; **canonical four-operator architecture**)

---

## Context (load into a fresh-context Claude)

### What Phase F left structurally unresolved

Phase F (commits `70bad29..b0cc1b1`, 2026-05-12) closed the headline ERR-026 manifestation #6 by backporting the Phase D Carlson coupled-pole seed to the SI/sweep path. Empirical wins were strong (sf[0]/sf[1] 0.522 → 0.778 STABLE; cv(ψ@i=0) 0.520 → 0.404; outer reflective sf[-1]/sf[-2] 0.887 → 0.997). But it patched a duplication rather than removing it — there are still TWO implementations of the M-M angular recurrence (one in `_sweep_1d_spherical`, one in `transport_operator_matvec_spherical`), and the user diagnosed the architectural smoking gun:

> "In principle, both the Krylov and the sweep are application of the Operators (including Boundary Operators) in a specific way/algorithm. Which means that if fixing this on the Krylov path didn't automatically fix the sweep, there is some inconsistency in the Operator formulation. It's not about auditing the sister path — the existence of a twin path itself is the smoking gun under the operator architecture."

The remaining manifestation #7 (residual O(h) SI-vs-Krylov drift) is the direct empirical signature of the duplication. Under correct operator algebra it CANNOT exist — they would literally execute the same code.

### The architectural target (user's verbatim mandate)

> "In the end we should have **L (streaming), C (collision), F (fission) and S (scattering)**. The end target should be that simple, and under the correct application, those operators look like math, and everything else using them looks like math as well (adjoint for example). This is why code elegance matters! It literally prevents bugs by construction."

Five user-facing operator types, period:
- **L** — `StreamingOperator` (`Ω · ∇`). Boundary operator B composed INSIDE L as constructor parameter (§16A).
- **C** — `CollisionOperator` (`Σ_t · ψ`). Diagonal, self-adjoint.
- **S** — `ScatteringOperator` (Legendre integral kernel `R · Λ · M`; harmonic factoring hidden behind `.factors`).
- **F** — `FissionOperator` (`F = χ ⊗ (νΣ_f)^T`).

Plus B (`BoundaryOperator`) — first-class `LinearOperator`, but composed inside L, never a top-level peer.

The user's elegance acceptance criterion (the adjoint solver):

```python
def solve_sn_adjoint(materials, mesh, quad, response, *, tol):
    V = DiscreteOrdinatesPhaseSpace.from_mesh(mesh, quad, materials.groups)
    L, C, S, F = build_sn_operators(V, materials)
    A_loss = L + C - S
    return GMRESSolver(tol=tol).solve(A_loss.H, response.as_source())
```

`.H` propagates through `OperatorSum` to leaves (`L.H` = reverse-direction sweep; `C.H = C`; `S.H` = transpose group transfer). **Zero new code beyond leaf `.H` implementations**. The adjoint writes itself.

### Resolved design questions (from reconciliation §5)

These are FIXED before migration starts — they shape the architecture, not just the migration sequence:

| # | Decision | Rationale |
|---|---|---|
| §5.1 | L exposes `L.solve_streaming(q)` (void / diagnostic) + `(L + C).solve(q)` via OperatorSum fusion (standard sweep). `L.solve_streaming` raises on non-void problems; users see `(L + C).solve` as the public sweep API. | User choice — exposes streaming-only for diagnostics while keeping the standard sweep as `(L + C).solve(q)`. |
| §5.2 | S hidden by default, `S.factors` property exposes `(R, Λ, M)` for users who want moment-level composition. | Keeps `A_loss = L + C - S` clean; enables PN / sensitivity work without polluting common API. |
| §5.3 | B composed INSIDE L's constructor: `L = StreamingOperator(V, sigma_t, boundary=B)`. Never a top-level summand. | User mandate — 4 user-facing operators, B is co-equal but not peer. |
| §5.4 | Angular integration = `M0 = M.restrict(ℓ=0)`; `phi = M0 @ psi`. | Most algebra-native; shares moment infrastructure with S = R·Λ·M factoring. |
| §5.5 | `CurvilinearSweepStrategy` dispatched by `V.mesh.geometry` inside L's `solve` / `apply` methods. Carlson seed + M-M recurrence become methods on the strategy. | L is one class — geometry dispatch is internal. Single closure config = manifestation #7 dissolves. |

---

## Goal

After Phase G ships, ALL of these hold simultaneously:

1. The four user-facing operator types L, C, S, F (+ composed B) are the ONLY public operator API. The audit memo's `SNCellOperator`, `SNSweepOperator`, `AngularRedistribution`, `PoleFaceAnchor`, `AngularIntegration`, `PackedStructuredAdapter`, `SNBoundaryFaceTraceOperator` are either dissolved (inside L), implementation details (Representation, Strategy), or composed expressions (`A_loss = L + C - S`).
2. `_sweep_1d_spherical`, `_sweep_1d_cylindrical` collapse to ≤ 10 lines that read like math (no procedural for-loops over ordinates / cells at the public level).
3. `transport_operator_matvec_spherical`, `_cylindrical` collapse to the operator-algebra expression form OR are retired entirely in favor of `L.apply` / `(L+C-S-F/k).apply`.
4. `_solve_fixed_source_si`, `_solve_fixed_source_krylov` collapse to ≤ 10 lines of `SourceIteration(L, S, F, q).solve(...)` / `PreconditionedGMRES(A_loss, precond=L).solve(...)` — no inline for-loops, no inline source assembly, no inline angle integration.
5. **Issue #196 manifestation #7 (SI-vs-Krylov O(h) WDD asymmetry) CLOSED by construction** — SI and Krylov route through the SAME `SNCellOperator` (private representation) with the SAME closure config. The Phase E flux-shape sentinel `test_phase_e_trajectory_resolvent_flux_shape_crosscheck` XPASSES; marker removed.
6. The adjoint flux solver `solve_sn_adjoint(materials, mesh, quad, response)` is implemented in ≤ 10 lines. The 6-line target above runs and matches a 2G slab reference at machine precision.
7. The 11 frozen regression snapshots pass at `assert_allclose(rtol=1e-12)` throughout the migration (bit-identity preserved unless principled per `vv-principles` §"Bit-identity vs principled-equivalence"; documented per snapshot if broken).
8. Sphinx narrative: `docs/theory/operator_algebra.rst` rewritten from Phase 0 stub into the full Wave H closure narrative; `discrete_ordinates.rst` Phase G section explains the unification.
9. ERR-026 manifestation table: #6 CLOSED (Phase F), #7 CLOSED (Phase G), #5 OPEN (Issue #195 — L1 magnitude pre-asymptotic, separate). ERR-026 narrows to single-manifestation OPEN.

---

## Step 1 — Promote `CellUpdate` + `MorelMontryAngularSweep` to `LinearOperator` (~3 days)

**Goal**: pure type-system promotion. The existing `DiamondDifference.update(visit, total_xs, source, upstream_state) → CellResult` becomes the `apply` / `solve` capability set on a `LinearOperator` subclass. No mathematical change; no behavioral change; bit-identical to current `np.array_equal`.

### Pre-step: test-architect dispatch

Before any code change, dispatch the **test-architect** to design the verification gates for Step 1:
- Bit-identity invariant: `SNCellOperator.solve(visit, ...) ≡ DiamondDifference.update(visit, ...)` to `np.array_equal` on every parameter combination already covered by `tests/sn/spatial/test_diamond.py`.
- Apply-vs-solve consistency: `SNCellOperator.apply(SNCellOperator.solve(q)) ≡ q` to `assert_allclose(rtol=1e-12)`.
- Capability surface: every capability declared on `LinearOperator` Protocol that `SNCellOperator` claims (`CAP_APPLY`, `CAP_SOLVE`, `CAP_TRANSPOSE` if applicable) must have a passing test.
- Curvilinear + Cartesian coverage.
- Verification gates DOC must call out: Phase D's twin-path bug shipped under apparently-green tests because the bit-identity invariant test only covered ONE of the twin paths.

### Implementation

1. **`SNCellOperator(LinearOperator)`** — wrap `DiamondDifference`. The `apply(cell_avg) → residual` method computes `L_cell · ψ_avg - q` (the cell's discretised operator residual). The `solve(spatial_upstream, angular_upstream, source) → cell_avg` method IS the current `DiamondDifference.update` body.
2. **`AngularRedistribution(LinearOperator)`** — wrap `MorelMontryAngularSweep.__call__` and the Carlson seed (`carlson_inward_sweep_from_source`) composition. The `apply(fi) → fi'` method computes the M-M angular recurrence transform. The Carlson seed becomes a private inner LinearOperator (also linear → composable).
3. **Both register with the existing `LinearOperator` Protocol** in `orpheus/numerics/operator.py`. Declare `capabilities = frozenset({CAP_APPLY, CAP_SOLVE})` per the type.

### Verification

- All tests in `tests/sn/spatial/test_diamond.py` and `tests/sn/spatial/test_psi_half_angle_seed.py` pass unchanged.
- New tests from the test-architect's gate design pass.
- 11 regression snapshots stay bit-identical to current.

### Deliverable

Commit pair: `feat(sn): SNCellOperator + AngularRedistribution LinearOperator promotion (Issue #196 Step 1)` + `test(sn): apply-vs-solve invariants for cell + angular ops (Issue #196 Step 1)`. Memo at `.claude/agent-memory/method-implementer/issue_196_phase_g_step1_closeout.md`.

---

## Step 2 — Unify L's closure: ONE WDD strategy for both `L.apply` and `L.solve` (~5 days)

**Goal**: this is the load-bearing step — where manifestation #7 dissolves. Both `_sweep_1d_spherical` (current SI sweep) AND `transport_operator_matvec_spherical` (current apply-matvec) internally compose the SAME `SNCellOperator` with the SAME closure config. They're allowed to remain as two free functions for this step, but the inline WDD recurrence + M-M recurrence + Carlson seed + BC dispatch are factored out and routed through the new operators.

### Pre-step: numerics-investigator dispatch

Dispatch numerics-investigator to characterise the residual O(h) drift between the current SI and Krylov paths at the per-cell level (the diagnostic that Phase F deferred to Issue #196). For each cell on `sphere_2g_3reg` n=40 at the converged fixed point, compute the per-cell, per-ordinate ψ from BOTH paths. The drift's location informs which primitive's closure must unify.

### Implementation

1. **Construct `StreamingOperator (L)`** in `orpheus/sn/operator.py` (or a new module). Constructor: `L = StreamingOperator(V_phase_space, boundary=B, sigma_t_source=None)`. The σ_t is NOT a constructor parameter on L — it lives on C.
2. **`L.solve_streaming(q)`** — streaming-only inverse via the existing curvilinear `iter_cells_by_direction` + `SNCellOperator.solve` over the DAG. Raises on problems where σ_t > 0 anywhere (the WDD doesn't converge without collision); useful for void / diagnostic problems.
3. **`L.apply(psi)`** — the directional derivative action. Routes through `iter_cells_by_direction` + `SNCellOperator.apply` cell-by-cell. The current `transport_operator_matvec_spherical` becomes a backward-compat alias.
4. **`(L + C).solve(q)`** — via OperatorSum fusion. The `(L + C)` summation expression detects the SN-context fusion and synthesizes a `FusedStreamingCollisionRepresentation` that runs the existing WDD sweep with `σ_t` from C. Internally, this IS today's `_sweep_1d_spherical` body, refactored to route through `SNCellOperator`.
5. **CurvilinearSweepStrategy** — dispatched by `V.mesh.geometry`. Holds the M-M recurrence + Carlson seed for sphere / cylinder.
6. **`C = CollisionOperator(V, sigma_t_field)`** — minimal, just `apply` and `solve` (trivial division).
7. **BC integration**: `L.boundary` is the realised B operator, applied ONCE at the boundary edge inside `(L + C).solve`. No per-ordinate `bc.apply` slicing inside the sweep body.

### Verification

- The Step 2 dispatch's per-cell diagnostic memo confirms the SI-vs-Krylov drift has REDUCED (or eliminated) at heterogeneous MR n=40 after unification.
- Issue #196 manifestation #7 closes empirically: `test_phase_e_trajectory_resolvent_flux_shape_crosscheck` xpasses → marker removed.
- Adjoint preview test: a slab `(L + C - S - F).H @ psi_dagger ≡ q_dagger` symmetric to the forward problem.

### Deliverable

Atomic commits per sub-feature. Memo at `.claude/agent-memory/method-implementer/issue_196_phase_g_step2_closeout.md` documenting manifestation #7 closure with empirical evidence + the new closure-configuration knobs.

---

## Step 3 — Wire `SourceIteration` + `PreconditionedGMRES` as algebra consumers (~3 days)

**Goal**: replace the 60-line `_solve_fixed_source_si` and 158-line `_solve_fixed_source_krylov` procedural dispatch functions with ≤ 10-line operator-algebra expressions. `SourceIteration` and `PreconditionedGMRES` already exist as primitives in `orpheus/numerics/iteration.py`; Step 3 just wires them.

### Implementation

1. **`SourceIteration(L, C, S, F, q_ext, *, max_iter, tol)`** — Richardson iteration: `ψ_{k+1} = (L+C).solve(q_ext + S @ ψ_k + (1/k_iter)·F @ ψ_k)`. Returns `(ψ, k_eff, residual_history)`. Consumes the four operators as constructor args.
2. **`PreconditionedGMRES(operator, precond)`** — wraps `scipy.sparse.linalg.gmres` with the operator's `apply` as `matvec` and the precond's `apply` as `M_matvec`. The operator expression `A_loss = L + C - S - F/k` is built outside.
3. **Update `_solve_fixed_source_si`**:
   ```python
   def _solve_fixed_source_si(...):
       L, C, S, F = build_sn_operators(V, materials, geometry)
       return SourceIteration(L, C, S, F=None, q_ext=external_source).solve(
           tol=inner_tol, max_iter=max_inner,
       )
   ```
4. **Update `_solve_fixed_source_krylov`**:
   ```python
   def _solve_fixed_source_krylov(...):
       L, C, S, F = build_sn_operators(V, materials, geometry)
       A_loss = L + C - S
       return PreconditionedGMRES(A_loss, precond=(L + C)).solve(
           q_ext=external_source, tol=inner_tol, max_iter=max_inner,
       )
   ```

### Verification

- All existing fixed-source tests pass unchanged (`tests/sn/regression/`, `tests/sn/test_phase_c_*`).
- The line count of `_solve_fixed_source_*` drops by ~85%.

### Deliverable

`refactor(sn): wire SourceIteration + PreconditionedGMRES as algebra consumers (Issue #196 Step 3)`.

---

## Step 4 — Promote BC to first-class via BoundaryRealizer (~3 days)

**Goal**: B is a `LinearOperator` composed inside L's constructor. The current per-ordinate `bc.apply(...)` slicing inside the sweep body is replaced by a face-trace-level `B @ psi_face_out → psi_face_in` operation done ONCE at the boundary edge.

### Implementation

1. **`BoundaryRealizer(V_trace_space, kind, params)`** — per-kind realiser for vacuum / reflective / albedo / symmetry. Returns a `LinearOperator` whose `apply(psi_outflow) → psi_inflow` is the BC trace law (§16A.4). Per the grand report §16B for spherical symmetry.
2. **`L = StreamingOperator(V, boundary=B)`** — B is a constructor arg; L's internal sweep applies B once per direction reversal.
3. **Retire `BoundaryFaceFlux` Protocol if not already done** (Phase C retired the Protocol but possibly left adapters; verify clean).
4. **All BC tests** at `tests/sn/test_snstreamingoperator.py` and `tests/sn/test_phase_c_*` continue to pass.

### Verification

- BC behavior bit-identical on all 11 regression snapshots.
- The per-ordinate `bc.apply` slicing inside `_sweep_1d_*` and `transport_operator_matvec_*` is removed.
- New test: `B.apply` linearity, idempotence (for symmetric BCs), surface measure consistency.

### Deliverable

`refactor(sn): promote BoundaryRealizer to first-class operator composed inside L (Issue #196 Step 4)`. Memo.

---

## Step 5 — Implement `.H` on each leaf → adjoint solver writes itself (~2 days)

**Goal**: the elegance acceptance criterion. Implement `.H` on `L`, `C`, `S`, `F` leaves. `.H` propagation through `OperatorSum` and `OperatorProduct` is already wired in `orpheus/numerics/operator.py` (Wave H Phase 0). `solve_sn_adjoint(materials, mesh, quad, response)` writes itself.

### Implementation

1. **`L.H`** = reverse-direction streaming. `iter_cells_by_direction(-1)` instead of `+1` (and vice versa). The sweep solves `L^T · ψ^† = source^†`.
2. **`C.H = C`** — collision is self-adjoint (diagonal, real).
3. **`S.H`** — transpose group transfer matrix and reverse moment-fold direction. Reuses S's existing harmonic infrastructure with `M.H` and `R.H`.
4. **`F.H`** — for `F = χ ⊗ (νΣ_f)^T`, `F.H = (νΣ_f) ⊗ χ^T` (swap rank-1 factors).
5. **`solve_sn_adjoint(materials, mesh, quad, response, *, tol)`**:
   ```python
   def solve_sn_adjoint(materials, mesh, quad, response, *, tol):
       V = DiscreteOrdinatesPhaseSpace.from_mesh(mesh, quad, materials.groups)
       L, C, S, F = build_sn_operators(V, materials)
       A_loss = L + C - S
       return GMRESSolver(tol=tol).solve(A_loss.H, response.as_source())
   ```
6. **`solve_sn_k_adjoint`** — k-eigenvalue adjoint via Arnoldi on `K.H = F.H @ A_loss.H.inverse()`.

### Verification

- 2G slab adjoint reference test (analytical k-eigenvalue adjoint for a homogeneous reflective slab). The adjoint flux should equal the forward flux at the same eigenvalue for a self-adjoint operator (`S = S.T` for L=0 isotropic).
- 2G slab response-weighted reaction rate test: `<response, ψ_dagger> = <forward_source, ψ_forward>` (reciprocity theorem).
- Sphere adjoint test (curvilinear): verify the adjoint flux solver runs and produces a physically-sensible result on a 2G heterogeneous reflective sphere.

### Deliverable

`feat(sn): implement .H on L/C/S/F leaves + solve_sn_adjoint (Issue #196 Step 5)`. New test module `tests/sn/test_adjoint.py`. Sphinx narrative additions.

---

## Cross-step deliverables

### Sphinx documentation campaign (parallel to Steps 2-5)

Dispatch **Archivist** in parallel starting at Step 2:

1. Rewrite `docs/theory/operator_algebra.rst` from current stub into the full Wave H closure narrative.
2. Cross-reference Phase F's narrative to point at Phase G's unification.
3. Add `docs/theory/sn_adjoint.rst` (or extend `discrete_ordinates.rst`) — the adjoint theory + how it falls out of `.H` propagation.
4. Wire `verifies(...)` decorators on the new adjoint tests against operator-algebra equation labels.

### Cross-solver leverage (out of scope for Phase G but the long-term pay-off)

Once L, C, S, F land in SN, the architecture pattern transfers to MoC, CP, diffusion (each with their own L, C, S, F realisations sharing the same Protocol). Per the user's `feedback_unify_after_two_instances.md` policy — don't unify across solvers until ≥2 instances exist; SN is instance 1, the next solver (likely MoC) is instance 2.

---

## Verification

### Acceptance criteria for Phase G closure

ALL of:
1. The 9 goals enumerated above.
2. The 11 frozen regression snapshots passed throughout (any bit-identity break documented per `vv-principles` §"Bit-identity vs principled-equivalence").
3. New foundation test `tests/sn/test_adjoint.py` ≥ 4 cases passing (slab self-adjoint reciprocity + sphere physical-sense + reaction-rate reciprocity).
4. Sphinx `-W` builds clean. `operator_algebra.rst` Phase 0 stub rewritten.
5. ERR-026 manifestation table: #7 CLOSED.
6. Issue [#196](https://github.com/deOliveira-R/ORPHEUS/issues/196) closes.
7. Phase E flux-shape sentinel xpasses; marker removed.

### Test sequence (end-to-end gate)

```bash
# Step 1 verification (~5 min): cell + angular op promotion
.venv/bin/python -m pytest tests/sn/spatial/ -q

# Step 2 verification (~10 min): manifestation #7 closure
.venv/bin/python -m pytest tests/sn/test_phase_c_*.py -v

# Step 3 verification (~10 min): SI / Krylov wiring
.venv/bin/python -m pytest tests/sn/regression/ -v

# Step 4 verification (~5 min): BC promotion
.venv/bin/python -m pytest tests/sn/test_snstreamingoperator.py tests/sn/test_bc_*.py -v

# Step 5 verification (~10 min): adjoint solver
.venv/bin/python -m pytest tests/sn/test_adjoint.py -v

# Phase G end-to-end gate (~25 min)
.venv/bin/python -m pytest tests/sn/ -v
```

---

## Sub-agent dispatch chain (proposed)

| Step | Pre-step | Implementation | Verification |
|---|---|---|---|
| 1 | test-architect (verification gates for type promotion) | method-implementer | qa |
| 2 | numerics-investigator (per-cell SI-vs-Krylov drift characterisation) | method-implementer | numerics-investigator |
| 3 | (no pre-step; pure wiring) | method-implementer | qa |
| 4 | test-architect (BC operator invariants) | method-implementer | qa |
| 5 | test-architect (adjoint reciprocity reference values + reciprocity theorem) | method-implementer | qa |
| Sphinx | (parallel) | Archivist | (final review by qa) |

---

## Risks + mitigations

- **Risk**: Step 2 closes manifestation #7 but discovers the residual closure-config drift is principled (e.g. WDD-asymmetric vs symmetric IS the correct choice, current is wrong). The "unification" might require choosing a closure config that breaks bit-identity on the homogeneous Cartesian snapshots.
  - **Mitigation**: the numerics-investigator pre-step characterises the drift; if it points at a principled closure refinement, document explicitly and accept the snapshot regen (per `vv-principles`). Worst case: Phase G splits into Phase G + Phase H.

- **Risk**: Step 5 `.H` implementations diverge subtly between L/C/S/F leaves (e.g. `L.H` reverses direction but the BC operator B is not its own adjoint for non-symmetric BCs).
  - **Mitigation**: test-architect pre-dispatch designs explicit reciprocity invariants (`<L @ ψ, ψ^†> = <ψ, L.H @ ψ^†>`) that catch any leaf `.H` bug. The L.H test must explicitly include curvilinear + non-symmetric BC cases.

- **Risk**: `(L + C).solve(q)` via OperatorSum fusion adds dispatch overhead vs the current direct sweep function call. Performance regression on the inner solve loop (called every SI iteration).
  - **Mitigation**: the FusedStreamingCollisionRepresentation cache is constructed ONCE per solver and reused per iteration. Benchmark on `sphere_2g_3reg` n=40: pre-Phase-G SI takes ~280 ms / iteration (per Step 2 numerics memo); post-Phase-G must stay within 2× that or refactor the fusion site.

- **Risk**: Step 3 SI / Krylov wiring exposes that `SourceIteration` / `PreconditionedGMRES` in `orpheus/numerics/iteration.py` are not yet capability-complete for the SN application (missing edge cases the procedural dispatch handled).
  - **Mitigation**: Step 3 dispatch includes a "capability gap audit" — list which procedural-dispatch features are not yet in the iteration primitives, extend the primitives, then wire.

- **Risk**: Step 4 BC promotion breaks the §16A.5 affine-BC contract for Phase C work. The current `bc.apply(...)` chain handles affine BCs with the `(R @ G) ψ + q_bc` form; the new `B @ psi_outflow` must too.
  - **Mitigation**: B's interface is `apply(psi_outflow) → psi_inflow` with the affine intercept carried INSIDE B (e.g. `B(psi) = R @ psi + q_intercept`). Existing affine-BC tests must pass unchanged.

---

## Out of scope (tracked separately)

- **Phase D + E + F Sphinx already in place** — Phase G ADDS a unification narrative; doesn't rewrite Phase D-F sections.
- **#170 / #171** F_N pillar — not blocked by Phase G; can land independently.
- **#195** L1 magnitude pre-asymptotic — separate ERR-026 manifestation #5, not closed by Phase G.
- **#193 / #194** Phase D follow-ups — pre-existing technical debt, not Phase G scope.
- **MoC / CP / diffusion operator-algebra extension** — out of scope per `feedback_unify_after_two_instances.md`; SN is instance 1.

---

## Final note for the post-compaction agent

When you pick this up:

1. **Re-read this plan first** before any code action.
2. **Re-load** in order:
   - `.claude/agent-memory/general-purpose/phase_g_four_operator_architecture_reconciliation.md` — the reconciliation memo, which has the canonical 4-operator architecture extracted from `neutron_transport_grand_report_v3.md`.
   - `.claude/agent-memory/explorer/issue_196_sn_operator_architecture_audit.md` — the audit, USEFUL for the duplication inventory (table form) but its 7-operator proposal is SUPERSEDED by the 4-operator reconciliation.
   - `.claude/agent-memory/method-implementer/issue_168_phase_f_closeout.md` — Phase F's tactical fix this builds on.
3. **Start with Step 1's test-architect dispatch** — the verification gates are non-negotiable before code change. Phase D's twin-path bug shipped under apparently-green tests; Step 1's tests must prevent that.
4. **At each step's empirical-decision point**, write a short closeout memo BEFORE proceeding. Each step is independently committable; do NOT batch.
5. **Bit-identity break on snapshots is expected at Step 2** (the closure unification regenerates the SI fixed point). Per `vv-principles` §"Bit-identity vs principled-equivalence": named intermediate + structurally-independent reference are required. The new SI fixed point IS the Krylov fixed point (the structurally-independent reference); regen and document.
6. **Phase E flux-shape sentinel re-enablement at Step 2 end** — the `xfail-strict` marker on `test_phase_e_trajectory_resolvent_flux_shape_crosscheck` must self-promote to xpass once Step 2 closes manifestation #7. If it doesn't, STOP and reinvestigate; do not improvise.
7. **`.H` correctness at Step 5 is the elegance acceptance criterion**, not just a feature. If `solve_sn_adjoint` is more than 6 lines of operator algebra, the architecture is wrong somewhere — re-scope.
8. **The user's standard**: "the implementation should read like math under the dunder methods implemented with operator expressiveness." Every step's deliverable must be defensible against this — if a sub-step produces procedural code that LOOKS like the old style, it's wrong.
