---
name: wave-t-t4-streaming-verification-spec
description: Verification SPEC for Wave T step T.4 — lifting StreamingOperator.apply into the tensor-network L = L_spatial + L_angular_redist decomposition. Pre-implementation test plan to drive turn-by-turn T.4 execution.
metadata:
  type: project
---

# Wave T T.4 — Streaming as `L_spatial + L_angular_redist` — verification spec

**Branch.** `refactor/sn-operator-algebra` (worktree
`moment-space-and-layering`). T.4 follows
`feedback_no_method_implementer_for_surgical_carves` — main agent
direct authorship with turn-by-turn user steering. This spec is the
test driver, not the implementation plan.

**Plan source.** `.claude/plans/wave_t_tensor_network.md` §6 step T.4
+ §1.2 (two views of the connection-coefficient operator) + §1.3
(expected complications) + §5 verification strategy + §7 risk row 1
(the highest-leverage curvilinear concern).

**Grand report.** §15 line 2003-2019 (V = X ⊗ Ω ⊗ G + native sum-of-
tensor-products form); §15.1 lines 2031-2045 (`L = Σ_axis (D_axis ⊗
Ω_axis ⊗ I_g)`); north-star line 5697 ("tensor products are
foundational").

**Bifurcation point (V&V pillars).** T.4 is a code-internal refactor
of an already-verified solver. The verification chain runs against
pre-T.4 numerics (regression / snapshot — see §6 for the
principled-equivalence three-criteria gate) backed by structurally-
independent ground truth at L1 — homogeneous `k_∞ = νΣ_f / Σ_a`
(closed-form pillar) for the eigenvalue claim AND P1 anisotropic MMS
(MMS pillar) for the flux-shape / convergence-order claim. **No
eigenvalue claim is asserted from MMS evidence.** Per
`vv-principles` §"The three pillars of verification": MMS is
source-driven and cannot prove eigenvalues; that is what L4-3 k_∞
homogeneous is for.

**T.3 precedent.** T.3 elected the SOTP/Route-B shape with default
nulp-relaxation gate per `vv-principles` §"Bit-identity vs
principled-equivalence". T.4 follows the same discipline AND extends
it: T.4's reductions reorder more aggressively than T.3 (the
spatial-axis sweep is sequential along nx; see §6 R3 for the
expanded drift-bound analysis).

---

## §1 Architectural questions (resolved before T.4b code is written)

These resolve the open Q1–Q5 in the dispatch brief AND the augmented
MA-Q1–MA-Q5 from the main agent. Each recommendation is the
test-architect's hard call; the user may override.

### Q1 — 2-D Cartesian fork — recommend **OPTION (c) HYBRID**

T.4 introduces the new `L_spatial + L_angular_redist` algebra and
routes the **1-D path** through it (slab, sphere, cylinder).
`_apply_2d_cartesian` stays untouched as the procedural path that
intercepts 2-D at `sn/operator.py:809-812`. T.4 adds a **structural-
pin test** (§4 A2D-1) that records the value of the procedural 2-D
output on a fixed fixture and asserts it remains bit-identical post-
T.4 — i.e. T.4 demonstrably did not touch the 2-D path's numerics.

**Why not (a) tensor-lift 2-D in scope.** The 2-D path has known
cell-centre-proxy semantics (`sn/operator.py:843-868`: the
`psi.boundary.face_view` is passive; the BC's `apply(outgoing)` writes
into bulk at the boundary cell). Reformulating `face_view` as a true
trace causes a 10% k_inf drift per the existing comment at L862-868.
That rewire is its OWN architectural payload — it requires the BC
realizers to gain a "proper composable algebra" (the comment names
this explicitly) — and bundling it into T.4 violates `feedback_unify_after_two_instances`
(unify after ≥2 working instances; we have one working 2-D path, not
two). The 2-D lift becomes a follow-up Wave T+ (or Wave O) once the
1-D tensor-network instances validate the abstraction's shape.

**Why not (b) leave 2-D entirely undocumented.** Option (c)'s pin
test is cheap insurance: it freezes the 2-D output as a regression
baseline so that ANY future drift (intentional or accidental) trips
the gate. Option (b) gives up that protection.

**MA-Q3 connection.** The user's `feedback_unify_after_two_instances`
applies directly here. T.4 has three working 1-D instances after the
lift (slab, sphere, cylinder); 2-D Cartesian is structurally
different (no curvilinear redistribution, no level_indices iteration,
no half-grid recurrence — the operator is `μ_x ∂_x + μ_y ∂_y`). The
2-D path's natural algebraic home is `L = L_x + L_y` with each
summand a `TensorProductOperator` per §15.1. That is genuinely
different architecture from the 1-D `L_spatial + L_angular_redist`
shape AND from the existing FD wavefront. The right call is to ship
the 1-D unification first, validate it across three geometries, then
return to 2-D with informed structure.

**Action.** Set the L4-4 test (2-D bit-identity) as part of T.4a's
snapshot capture; verify in T.4e regression sweep. NO new tensor-
network code touches `_apply_2d_cartesian` in T.4 scope.

### Q2 — Curvilinear angular redistribution factor structure — recommend **OPTION (iii) PROCEDURAL SUMMAND**

`R_polar` (M-M angular redistribution) is **NOT a clean tensor-
product factor.** The Morel-Montry half-grid (`pole_angular_closure.py:487+`)
is sequentially coupled along the angular axis via Hébert §3.9.4
Eqs. 3.432-3.435 — `α_{m+1/2}` recurs from `α_{m-1/2}` with
absorption coefficients that depend on the current cell's
cross-section AND the upstream half-angle ψ. This is a SWEEP
operator, not a diagonal-on-angular-axis operator. There is no clean
`R_polar = D_μ ⊗ A_polar` factorisation; option (ii) is aspirational.

Recommendation: ship `L_angular_redist` as a **procedural summand**
— a bespoke `LinearOperator` leaf (NOT a `TensorProductOperator`)
that wraps M-M's existing `cell_contribution` algebra. The leaf
advertises `{CAP_APPLY}` and has explicit `axes` metadata for code-
inspection. The summand participates in the operator-algebra at the
OUTER `L = L_spatial + L_angular_redist` level via `OperatorSum`,
NOT as a `SumOfTensorProductsOperator`.

**Why not (i) `AngularSweepOperator & I_x & I_group` (the 3-factor
wrap).** It looks like a tensor product, but the leaf factor is the
entire sweep recurrence — wrapping it with `I_x & I_group` adds zero
structural information; it is `assert_separable`-failing by axes
analysis (the leaf factor touches `(N, nx)` jointly via the
sequential recurrence). The wrap is cosmetic; it asserts a structure
that doesn't actually hold. **NEVER assert separability the structure
doesn't support** — flagged by `coding-elegance` Pattern 4 (make
illegal states unrepresentable — the converse: don't represent
states that aren't actually legal).

**Why not (ii) `R_polar = D_μ ⊗ A_polar`.** The angular axis is
sequentially coupled; there is no clean diagonal factor. Pushing this
factorisation would require introducing a new "lower-triangular
operator" primitive whose `apply` is a sweep — which IS what
option (iii) already does at the leaf level, just without claiming
to be a tensor product.

**Plan §7 risk 1 mitigation explicitly allows this fallback:**
"factor what's factorable, document non-factorisable pieces, ensure
the algebra still composes via OperatorSum. The architectural
commitment is ONE algebraic FORM, not that every factor is a clean
`TensorProductOperator`."

**MA-Q1 connection — the master condition.** T.3 (per-material per-ℓ
einsum coupling group + spatial) and T.4 curvilinear (sequential
angular redistribution coupling ordinates along the half-grid) both
hit the same wall. **The structural pattern is: SOTP requires
disjoint-axes per factor; "coupled physics" (per-material XS lookup,
per-ordinate recurrence, BC face-trace propagation) couples axes by
design.** The master condition for "this CAN be SOTP" is:

> **SOTP requires that EACH summand decomposes as a Cartesian product
> of independent per-axis operations** — `f(x_1, ..., x_d) = f_1(x_1)
> ⊗ ... ⊗ f_d(x_d)`. If the physics couples two axes inside a
> summand (e.g., per-material XS makes group depend on cell;
> per-half-grid α makes angular ordinates depend on cross-section
> AND upstream angular state), the summand is **NOT SOTP-able**.
> Falls back to `OperatorSum` over `LinearOperator` summands.

This is load-bearing for Wave O's typing commitments — BulkOperator
and FullOperator MUST accommodate `OperatorSum` summands that are
NOT `TensorProductOperator`. The grand-report aspiration "everything
is SOTP" is contradicted by the physics in two of the three Wave T
substeps (T.3, T.4-curvilinear); only T.1 (BC realizers) and T.2
(fission rank-1) cleanly support SOTP. Wave O should not assume
universal SOTP-ability.

**Action.** Define the bespoke summand type explicitly. Tentative
name: `AngularRedistributionOperator` (in
`orpheus/sn/operator.py` adjacent to `StreamingOperator`). It is a
`LinearOperator` leaf, NOT a `TensorProductOperator`. The L2-1 test
(§4) asserts `isinstance(L_angular_redist, AngularRedistributionOperator)`
for curvilinear AND `isinstance(L_angular_redist, ZeroOperator)` for
slab (per Q4).

### Q3 — The `L = M - C` subtraction location — recommend **OPTION (γ) APPLY-BOUNDARY**

`StreamingOperator.apply` retains the Resolution A subtractive
contract (the cell-collision σ_t·ψ subtraction lives at the apply
boundary, OUTSIDE `L_spatial` and `L_angular_redist`). Per the
existing docstring at `sn/operator.py:656-668`:

> `L.apply(ψ) := M(ψ; σ_t) − σ_t ⊙ ψ.bulk`

The decomposition exposes the math primitives at the **property
level**:

```python
class StreamingOperator:
    @cached_property
    def L_spatial(self) -> SumOfTensorProductsOperator: ...
    @cached_property
    def L_angular_redist(self) -> LinearOperator: ...

    def apply(self, psi):
        # M = L + C, where L = L_spatial + L_angular_redist,
        # so M(ψ; σ_t) = L_spatial(ψ) + L_angular_redist(ψ) + σ_t · ψ
        M = self.L_spatial.apply(psi) + self.L_angular_redist.apply(psi)
        # Resolution A subtractive contract: L = M - C
        return M_minus_sigma_t_psi(M, psi, self.sigma_t)
```

**Why not (α) inside L_spatial.** σ_t · ψ is per-cell scalar
(broadcast across `(N, ng, nx, ny)`); it has no natural home inside a
per-axis spatial tensor-product summand. Distributing it across the
spatial-axis summands ("each axis subtracts its share") requires a
fictive σ_t-per-axis decomposition — non-physical, non-elegant.

**Why not (β) inside the L_spatial+L_angular_redist sum.** Same
problem: σ_t · ψ is a CELL quantity, not a spatial OR angular
quantity. Burying it in a non-streaming summand fights the elegance
mandate.

**Why (γ) is right.** The math primitives `L_spatial` and
`L_angular_redist` produce the **un-subtracted** streaming-plus-
collision contributions — they ARE `M_spatial` and `M_angular_redist`
in the L-vs-M convention. The user-facing `L.apply` algebra is
`M_spatial(ψ) + M_angular_redist(ψ) − σ_t · ψ`. This mirrors T.3's
Pattern 7 producer-side normalisation (the `/sum_w` factor stays at
the apply boundary, OUTSIDE the kernel summands).

**Naming convention.** To avoid the L-vs-M confusion, **rename the
properties**: `M_spatial` and `M_angular_redist` (NOT `L_spatial`
and `L_angular_redist`). The grand-report names them `L_spatial`,
but that's the continuous L (without collision); in the discrete
matvec the spatial term comes pre-coupled with collision per the
cell-balance algebra (`cell_balance_for_streaming` at
`sn/spatial/cell_balance.py:120-198`). The denom carries
`streaming_denom_term + angular_denom_term + collision_denom_term` —
all three additive — and the natural property names mirror that.

**Action.** Properties named `M_spatial` and `M_angular_redist`. The
operator-algebra invariant test (§4 A-3) pins
`(L + C).apply(ψ) ≡ M_spatial(ψ) + M_angular_redist(ψ)` bit-exact
(modulo principled-equivalence ULP if applicable per §6 R3).

**MA-Q2 connection.** The cell-balance algebra (`denom = streaming_denom_term
+ angular_denom_term + collision_denom_term`; `numer_upstream = spatial_upstream_term
+ angular_numer_upstream`) decomposes into THREE additive terms; the
T.4 operator-algebra decomposition is TWO summands (`M_spatial`,
`M_angular_redist`) plus the apply-boundary subtraction. These are
**consistent, not in conflict**:

- `M_spatial(ψ)` produces the contribution of (`streaming_denom_term`
  + `collision_denom_term`) · ψ_cell − `spatial_upstream_term` —
  i.e., per-axis WDD streaming WITH the collision share that lives at
  the cell.
- `M_angular_redist(ψ)` produces the contribution of
  `angular_denom_term` · ψ_cell − `angular_numer_upstream` — the
  redistribution piece alone, geometry-blind by data.
- `M.apply(ψ) := M_spatial(ψ) + M_angular_redist(ψ)` then equals
  `(L + C).apply(ψ)` per the cell-balance algebra's additive
  structure; the apply boundary subtracts `σ_t · ψ` to deliver L.

T.4 is **NOT re-architecting the cell-balance contract** — it is
exposing the existing three-term decomposition at the operator
level. The cell-balance is still the single source of truth; T.4
groups its terms (streaming + collision → M_spatial; angular →
M_angular_redist). This passes Pattern 2 (single source of truth):
both M_spatial and M_angular_redist internally invoke
`cell_balance_for_streaming` and `pole_angular_closure.cell_contribution`
— NO duplication, just reorganisation of which summand drives which
call.

### Q4 — Slab degeneracy — recommend **STRUCTURAL-PRESERVATION FORM**

`L_angular_redist` is a `ZeroOperator` (or equivalent) for slab. The
property structure is preserved — slab and curvilinear both expose
`M_spatial` AND `M_angular_redist`; the slab's `M_angular_redist`
just happens to be the zero operator.

**Why.** Per `coding-elegance` Pattern 4 (illegal states
unrepresentable) and Pattern 7 (single producer site): the structural
shape `M = M_spatial + M_angular_redist` is **uniform across
geometries**. Slab is the degenerate case where the redistribution
contribution is identically zero (because `IdentityAngularClosure`
returns zeros from `cell_contribution`). Hiding the property for
slab would force every consumer to check geometry before accessing
the operator algebra — a layering violation.

`ZeroOperator` is cheap: its `apply(ψ)` allocates a zero output
matching the codomain shape (or returns a sentinel that
`OperatorSum.apply` short-circuits — minor optimisation, deferred).

**Action.** L2-3 (§4) asserts:
- For slab: `isinstance(op.M_angular_redist, ZeroOperator)` AND
  `op.M_angular_redist.apply(ψ).values == 0` bit-exact.
- For sphere/cylinder: `isinstance(op.M_angular_redist,
  AngularRedistributionOperator)`.

### Q5 — `InvertibleOperator` interaction — **CONFIRMED OUT OF SCOPE**

`StreamingOperator.__add__(other)` at `sn/operator.py:958-969`
dispatches to `InvertibleOperator(self, other)` when `other` is
`CollisionOperator`. The `InvertibleOperator.solve` path is OUT OF
SCOPE per plan §3.5 and §1.3 second bullet (it's the WDD procedural
sweep; not tensor-product-factorable).

**Contract.** T.4 MUST NOT touch:
- `StreamingOperator.__add__` (the L+C composition dispatcher).
- `InvertibleOperator.solve` / `_solve_timed_full_field` /
  `transport_sweep` (the procedural inverse).

The `__add__` behaviour stays identical because the dispatched
`InvertibleOperator(self, other)` still receives the same
`StreamingOperator` instance — only its `.apply` has been
re-implemented internally; the `__add__` API is unchanged.

**Action.** L4-5 / L4-6 (§4) explicitly test the `(L+C).apply(ψ)`
composition path AND the `(L+C).solve(q)` path on the same fixture,
verifying both are unchanged post-T.4. The `.solve` test is the
defensive guard against accidental modification of `.__add__`.

---

## §2 Substep gate sequence (T.4a → T.4f)

T.4 carries more architectural payload than T.3 (touches three
geometries + introduces a bespoke leaf type + reorganises the
StreamingOperator class layout). Each substep is a separate commit
gated on the listed tests passing.

### Substep T.4a — Capture pre-T.4 snapshots + perf baseline

**Action.** Before any T.4 code change, run a one-shot script at
`tools/wave_t/capture_pre_t4_snapshots.py` (or inline in the spec's
companion script) that:
1. Exercises `StreamingOperator.apply` on fixed (seed, mesh, mat)
   triples across slab/sphere/cylinder 1-D AND 2-D Cartesian.
2. Captures `(L+C).apply(ψ)` (the composite) at each fixture — used
   to verify the algebra-decomposition invariant in T.4e.
3. Captures `(L+C).solve(q)` at each fixture — pinning the
   InvertibleOperator path is unchanged (Q5 contract).
4. Captures wall-clock baseline for the 1-D slab `(L+C).apply` Krylov
   benchmark (≥1000 matvecs to amortise warmup). This is the L5-1
   reference walltime.

Snapshot file: `tests/sn/_fixtures/wave_t_t4/pre_t4_snapshots.npz`.
Perf baseline: `tests/sn/_fixtures/wave_t_t4/pre_t4_walltime.json`.

**Gate (T.4a).**
- Snapshot file exists; schema validated (one ndarray per dispatch
  arm × geometry × group count).
- Perf baseline file exists; walltime recorded with deterministic
  seed.

**Why first.** Per plan §5.1 — bit-identity (or principled-
equivalence) is the contract. Cannot demonstrate post-T.4 numerics
without pre-T.4 numbers in hand. ALSO: T.3g was deferred because no
pre-T.3 baseline was captured at T.3a (the spec called for it but
the implementer skipped it under time pressure). T.4f is enforceable
ONLY if T.4a captures the baseline.

### Substep T.4b — Slab `M_spatial` (no curvilinear touched)

**Action.** Introduce the `M_spatial` cached_property on
`StreamingOperator` for the **slab path ONLY**. The property
returns a `SumOfTensorProductsOperator` with ONE summand for slab
(`D_x & Ω_x & I_g` — the §15.1 1-D form). Wire `apply` to delegate
to `M_spatial.apply(ψ) − σ_t · ψ` for slab geometry; curvilinear
arms continue to call `transport_operator_matvec_unified` directly
as the fallback (parallel paths during the substep).

**Gate (T.4b).**
- L2-1: kernel structure for slab — `M_spatial` is a SOTP with one
  summand; the summand is a 3-factor `TensorProductOperator`.
- L2-3: slab `M_angular_redist` is a `ZeroOperator`.
- L4-1 slab: bit-identity (or principled-equivalence per §6 R3)
  against pre-T.4 snapshot.
- L4-3 k_∞ slab: closed-form match within existing tolerance.
- A-1 algebraic identity for slab: `assert_separable()` passes on
  `M_spatial`; on `M_angular_redist`, `apply(zeros) == zeros`.

**Why slab first.** Slab is the structurally simplest case (no
angular redistribution, no level_indices, no half-grid). If the
abstraction's shape fails for slab, it will fail catastrophically
for curvilinear. **NEVER ship the harder case before the easier case
that shares the abstraction** (Pattern 6 + general-first development
order per `coding-elegance` Pattern 15).

### Substep T.4c — Curvilinear `M_angular_redist` wiring

**Action.** Add the bespoke `AngularRedistributionOperator` leaf
type. Wire `M_angular_redist` for sphere and cylinder. Wire the
sphere `M_spatial` (one summand). Wire the cylinder `M_spatial` (one
summand). The leaf's `apply` body invokes
`pole_angular_closure.cell_contribution` per cell, mirroring the
existing `transport_operator_matvec_unified` body but extracted to
the leaf level.

**Gate (T.4c).**
- L2-2: each summand has the chosen factor count
  (`EXPECTED_FACTOR_COUNT=3`).
- L2-4: capabilities = `{CAP_APPLY}` for both M_spatial and
  M_angular_redist (composition inherits via OperatorSum
  capability intersection).
- L4-2 sphere: bit-identity (or principled-equivalence) against
  pre-T.4 snapshot.
- L4-2 cylinder: same.
- L4-3 k_∞ sphere and cylinder: closed-form match.
- A-2: `M.apply(ψ) ≡ M_spatial.apply(ψ) + M_angular_redist.apply(ψ)`
  bit-exact (Resolution A invariant; defines the algebra contract).

**Why curvilinear in one substep, not two.** Sphere and cylinder
share the same `MorelMontryAngularSweep` strategy (different
`level_indices`); the leaf's body is geometry-blind by data. Splitting
sphere/cylinder into two substeps doubles the commit count without
catching new bug classes. The L4-2 fixture covers both geometries.

### Substep T.4d — 2-D Cartesian regression pin (NO new code)

**Action.** No code change here — this is a verification step. The
2-D Cartesian path (`_apply_2d_cartesian`) was NOT touched by T.4b
or T.4c; this substep VERIFIES that.

**Gate (T.4d).**
- L4-4: 2-D Cartesian bit-identity against pre-T.4 snapshot
  (specular all sides AND vacuum all sides — both fixtures).

**Why this is a separate substep.** Per plan §6 T.3 paragraph 3
and T.3 spec §4 substep T.3e, the "defensive substep" pattern pins
contracts that the author MUST NOT modify while in the file. The
explicit gate catches a T.4 author who accidentally modernises
`_apply_2d_cartesian` while it's open.

### Substep T.4e — Composition invariants + L1 MMS gates

**Action.** Run the full L1 suite. Run the algebra-composition
invariant tests.

**Gate (T.4e).**
- A-3: `(L + C).apply(ψ)` value-equal to pre-T.4 snapshot for slab,
  sphere, cylinder, 2-D Cartesian. (This is the §1.3 Q3 contract
  — `(L+C).apply` is bit-identical to pre-T.4, because the
  subtraction at the apply boundary unwinds.)
- A-4: `(L + C).solve(q)` value-equal to pre-T.4 snapshot. Verifies
  Q5 contract — the `.solve` path was not touched.
- L1-1: `tests/sn/test_mms_aniso.py` 1-D slab P1 aniso O(h²) gate
  stays green.
- L1-2: `tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py`
  curvilinear gate stays at the same xfail set.
- L1-3: `tests/sn/test_mms.py` (P0 MMS) stays green.
- L1-4: `tests/sn/test_mms_heterogeneous.py` stays green.

**Why mandatory.** Per plan §5.2 line 162: "the L1 MMS gates are the
ground truth. They must stay green at every Wave T commit. If a
substep breaks an L1 gate, the substep is wrong (regardless of
bit-identity claims — bit-identity is necessary but not sufficient
when reductions reorder)."

### Substep T.4f — Performance regression gate

**Action.** Run the 1-D slab `(L+C).apply` Krylov benchmark.

**Gate (T.4f).**
- L5-1: ≤ 5% slowdown vs pre-T.4 walltime captured at T.4a.

**Why last.** Optimisation is a separate concern from correctness;
shipping a correct but slow T.4 is a regression but not a
disqualifier. If exceeded, investigate factor caching (the
`cached_property` is already in place) or SOTP factor fusion
(deferred per plan §4).

### Master gate: COMMIT only when ALL of {T.4a–T.4f} pass

Per `feedback_no_method_implementer_for_surgical_carves`: each
substep is a separate commit IFF the previous gate passes cleanly.

---

## §3 Snapshot capture (T.4a) — explicit fixture schema

The pre-T.4 snapshot script writes a single npz file with the
following entries. Deterministic seeds: `np.random.default_rng(42)`
for ψ; fixed mesh sizes; fixed material data from
`orpheus.derivations._xs_library`.

| Key | Shape | Description |
|-----|-------|-------------|
| `slab_1g_apply` | `(N, 1, nx, 1)` | `StreamingOperator.apply(ψ)` for slab 1G P0 vacuum BC, nx=40, S8 GL quadrature |
| `slab_1g_apply_boundary` | `(N, 1)` per face (flat per layout) | corresponding face residual for slab (xmin + xmax) |
| `slab_2g_apply` | `(N, 2, nx, 1)` | slab 2G P1 asymmetric SigS vacuum BC, nx=40, S8 GL |
| `slab_2g_apply_boundary` | per layout | corresponding face residual |
| `slab_2g_white_apply` | `(N, 2, nx, 1)` | slab 2G white-BC variant (catches BC convention drift) |
| `sphere_1g_apply` | `(N, 1, nx, 1)` | sphere 1G P0 vacuum at r=R, nx=20, S8 GL |
| `sphere_2g_apply` | `(N, 2, nx, 1)` | sphere 2G P1 asymmetric SigS vacuum at r=R |
| `sphere_2g_apply_boundary` | per layout | corresponding face residual at outer face |
| `cyl_1g_apply` | `(N, 1, nx, 1)` | cylinder 1G P0 vacuum at r=R, nx=20, LevelSym S6 |
| `cyl_2g_apply` | `(N, 2, nx, 1)` | cylinder 2G P1 asymmetric SigS vacuum |
| `cart2d_1g_specular` | `(N, 1, nx, ny)` | 2-D Cartesian 1G specular all sides, nx=ny=10, ProductGL S6 |
| `cart2d_1g_vacuum` | `(N, 1, nx, ny)` | 2-D Cartesian 1G vacuum all sides, same mesh/quad |
| `cart2d_2g_specular` | `(N, 2, nx, ny)` | 2-D Cartesian 2G P1 specular all sides |
| `LpC_slab_2g_apply` | `(N, 2, nx, 1)` | `(L+C).apply(ψ)` slab 2G — verifies algebra-decomposition |
| `LpC_sphere_2g_apply` | `(N, 2, nx, 1)` | same for sphere |
| `LpC_cyl_2g_apply` | `(N, 2, nx, 1)` | same for cylinder |
| `LpC_cart2d_2g_apply` | `(N, 2, nx, ny)` | same for 2-D Cartesian |
| `LpC_slab_2g_solve` | `(N, 2, nx, 1)` | `(L+C).solve(q)` slab 2G — verifies `.solve` path untouched |
| `LpC_sphere_2g_solve` | `(N, 2, nx, 1)` | same for sphere |
| `LpC_cyl_2g_solve` | `(N, 2, nx, 1)` | same for cylinder |
| `seed_psi_<arm>` | match shape | the input ψ values used for each capture (seed traceability) |
| `seed_q_<arm>` | match shape | the input q values for the `.solve` captures |

**Why 2G mandatory.** Per `vv-principles` §"1-group degeneracy
canonical statement": 1-group eigenvalue tests are degenerate
(`k = νΣ_f/Σ_a` is flux-shape independent). The 1G captures are kept
for **regression integrity** at the lowest dimension, but the
**load-bearing** verification rows are 2G with asymmetric SigS —
they exercise group coupling and catch ERR-002-class transpose
convention drift.

**Why asymmetric SigS.** Per `numerical-bug-signatures` Signature 3
(ERR-002 ScatteringMatrix transpose) — although T.4 doesn't touch
scattering, the M-M angular redistribution couples the cell-balance
arithmetic across all groups via `total_xs` and `sigma_t_gx`. A
2-group symmetric XS would mask group-axis bugs in M_angular_redist.
The fixture uses
`orpheus.derivations._xs_library.get_mixture("A", "2g")` material A
(asymmetric).

**Why white BC for slab.** The white BC angular average couples all
ordinates of opposite sign — a different convention path than
vacuum. T.3 didn't need this; T.4 introduces `M_spatial` whose
treatment of incoming face flux passes through the BC realizer. White
BC catches the regression where the realizer convention drifts.

**Why specular AND vacuum for 2-D.** Specular exercises the
permutation-operator BC path (the D-B+1 tensor-network instance);
vacuum exercises the inflow-mask BC. The 2-D path is currently
procedural FD; both BCs must remain unchanged.

**Perf baseline schema** (`pre_t4_walltime.json`):
```json
{
  "fixture": "slab_2g_p1_n_cells_40_S8GL_1000_matvecs",
  "walltime_seconds_median": <float>,
  "walltime_seconds_p95": <float>,
  "n_iterations": 1000,
  "python_version": "<sys.version>",
  "numpy_version": "<numpy.__version__>",
  "rng_seed": 42
}
```

---

## §4 New algebraic-identity tests (T.4-specific)

Per plan §5.3 and the architectural decisions in §1 above.

### 🔧 Post-T.4b implementation deviations (added 2026-05-31, T.5.1 amendment)

**The §4 test table below was authored against the SOTP form for `M_spatial`** (one summand per spatial axis, each a clean 3-factor `TensorProductOperator` `(D_axis & Ω_axis & I_g)`).  **The shipped form deviates** per the AskUserQuestion architectural follow-up that landed BEFORE T.4b code was written (the user surfaced the inconsistency between Q1's SOTP claim and Q2's bespoke-leaf treatment — same MA-Q1 master condition applies to both).

The user-endorsed honest decomposition (commit `c55b505` T.4b):

* `M_spatial` IS `_MSpatialOperatorSum(OperatorSum)` — a SUBCLASS of `OperatorSum`, NOT `SumOfTensorProductsOperator`.
* `M_spatial.a` and `M_spatial.b` are `_SpatialSweepDirection` BESPOKE LEAVES — NOT `TensorProductOperator` instances.  Each carries `direction_sign ∈ {+1, -1}` metadata.
* `M_spatial.apply` OVERRIDES the default `OperatorSum.apply` to run ONE bidirectional sweep with shared state (Design B — see plan §6 T.4 and the T.4b commit message for the rationale).

**Revised test rows (what actually shipped in `tests/sn/test_streaming_operator.py`)**:

| Spec row | Original claim (table below) | Shipped form (T.4b/T.4c) | Implementation |
|----------|------------------------------|--------------------------|----------------|
| **L2-1** | `op.M_spatial` is `SumOfTensorProductsOperator` with `len(summands) == n_spatial_axes` | `op.M_spatial` is `_MSpatialOperatorSum(OperatorSum)` with 2 `_SpatialSweepDirection` summands (forward + backward, regardless of spatial dimension count) | `TestT4bMSpatialStructure.test_M_spatial_is_operator_sum_with_two_per_direction_summands_slab` |
| **L2-2** | Each summand is a 3-factor `TensorProductOperator` `(D_axis & Ω_axis & I_g)` | Each summand is a `_SpatialSweepDirection` BESPOKE LEAF — NOT a `TensorProductOperator` wrap.  Test **DROPPED** (the claim was structurally false per MA-Q1) | (no test) |
| **L2-7** | Summand's 3 factors act on disjoint axes (D on spatial, Ω on angular, I on group) | Summands' `direction_sign` values are disjoint (one +1, one -1).  Direction-sign is the structural exposure; axes-disjoint claim doesn't apply to bespoke leaves | `TestT4bMSpatialStructure.test_M_spatial_summand_direction_signs_disjoint` |
| **L2-5** | `M_spatial.assert_separable()` passes (L1 invariant fires) | Test **DROPPED** — `assert_separable` is structurally inapplicable to `OperatorSum`-of-bespoke-leaves.  The narrower per-operator separability check fires only on T.1 BC realizers and T.2 fission TP | (no test) |
| **A-1** | Slab algebra-decomposition invariant at `rtol=1e-13, atol=1e-14` (Route A bit-id) | `np.testing.assert_array_equal` strict bit-identity for slab (M_spatial.apply IS bit-exact to the unified matvec for slab) | `TestT4bAlgebraDecompositionInvariantSlab.test_slab_apply_decomposition_invariant` |
| **A-2 / A-3** | Curvilinear algebra-decomposition invariant at `rtol=1e-13, atol=1e-14` | Curvilinear: `np.testing.assert_allclose(rtol=1e-14, atol=1e-14)` per the principled-equivalence ULP gate (M_spatial = `(L+C) − M_ang` introduces ~16·ULP drift) | `TestT4cAlgebraDecompositionInvariantCurvilinear.test_{sphere,cylinder}_LpC_equals_M_spatial_plus_M_angular_redist` |
| **A-9** | `M_spatial.apply(zeros) == zeros` | Test **DROPPED** as redundant with A-1 + A-8 algebra invariants (any failure here would surface via the additive decomposition tests at a stronger gate) | (no test) |
| **L4-4** | 2-D Cartesian bit-identity to pre-T.4 | Covered via A2D-1 defensive source-hash pin + the StreamingOperator.apply 2-D dispatch is untouched (verified inline at T.4b/T.4c snapshot re-capture; no pytest case shipped — overlaps with A2D-1) | `TestT4dApply2DCartesianSourceHashPin.test_apply_2d_cartesian_source_hash_unchanged` |
| **L4-5 / L4-6** | `(L+C).apply` / `(L+C).solve` bit-identical to pre-T.4 snapshot | Verified inline at snapshot re-capture (commits `cb18fdb`/`c55b505`/`90e7d4e` close-out checks all bit-identical).  No pytest case shipped — the inline check is the load-bearing gate for the algebra-decomposition contract; explicit test could be added in a future maintenance pass | (verified inline; no test row) |
| **L5-1** | `(L+C).apply` walltime ≤ 1.05× pre-T.4 baseline on 1000-iter slab P1 | Verified inline at T.4c close-out (median 1.040×, p95 0.998×; commit `90e7d4e`).  200 iterations instead of 1000 for fixture wallclock cost.  No pytest case shipped — perf gate would need a marker like `@pytest.mark.slow` to avoid CI noise; deferred | (verified inline; no test row) |

The **shipped test count** is 21 new tests in `tests/sn/test_streaming_operator.py` (vs the spec's projected 27 = 18 foundation + 7 L4 + 1 L5 + 1 A2D-1).  The 6 missing tests are either DROPPED (L2-2, L2-5, A-9 — structurally inapplicable to the OperatorSum/bespoke-leaf shape) or verified inline (L4-4, L4-5, L4-6, L5-1 — the snapshot re-capture proves bit-identity at close-out time; an explicit pytest case is redundant for these but could be added in a future maintenance pass for CI discipline).

The **architectural value** of the deviation is captured in the T.4 spec §1 "post-implementation deviations" subsection and the Wave T plan §6 T.4 + §9 amendments (T.5.1, also landing in this commit).

---

**Original §4 table** (the spec as authored 2026-05-30):



| #   | Test name | File | Class | Level | What it verifies | Reference | Pass criterion |
|-----|-----------|------|-------|-------|------------------|-----------|----------------|
| L2-1 | `test_M_spatial_is_sum_of_tensor_products` | `tests/sn/test_streaming_operator.py` | `TestTensorNetworkDecomposition` | foundation | `op.M_spatial` is a `SumOfTensorProductsOperator` with one summand per spatial axis (slab/sphere/cylinder: 1; 2-D Cartesian: out of scope per Q1) | Type-introspection | `isinstance(op.M_spatial, SumOfTensorProductsOperator)` AND `len(op.M_spatial.summands) == op.sn_mesh.n_spatial_axes` |
| L2-2 | `test_each_M_spatial_summand_is_3_factor_TP` | same | same | foundation | Each summand is a 3-factor `TensorProductOperator` per §15.1 `D_axis & Ω_axis & I_g` | Type-introspection | `len(summand.ops) == 3` AND each factor is a `LinearOperator` |
| L2-3 | `test_M_angular_redist_slab_is_zero` | same | same | foundation | Slab's `M_angular_redist` is `ZeroOperator` (or equivalent); `M_angular_redist.apply(ψ).values == 0` bit-exact | Algebraic ground truth (slab has no curvilinear redistribution) | `isinstance(op.M_angular_redist, ZeroOperator)` for slab; `np.testing.assert_array_equal(zeros)` for curvilinear at uniform-ψ probe |
| L2-4 | `test_M_angular_redist_curvilinear_is_bespoke_leaf` | same | same | foundation | Sphere and cylinder's `M_angular_redist` is the bespoke `AngularRedistributionOperator` leaf — NOT a TensorProductOperator (Q2 decision) | Type-introspection | `isinstance(op.M_angular_redist, AngularRedistributionOperator)` for sphere AND cylinder |
| L2-5 | `test_M_spatial_assert_separable_passes` | same | same | foundation | `M_spatial.assert_separable()` passes (the L1 invariant fires) | Built-in invariant | `M_spatial.assert_separable()` returns None (no raise) |
| L2-6 | `test_M_spatial_capabilities_apply_only` | same | same | foundation | `M_spatial.capabilities == {CAP_APPLY}` AND `M_angular_redist.capabilities == {CAP_APPLY}` | Capability set algebra | `frozenset({CAP_APPLY})` for both |
| L2-7 | `test_disjoint_axes_per_M_spatial_summand` | same | same | foundation | Each `M_spatial` summand's factors act on disjoint axes (D_axis on spatial; Ω_axis on angular; I_g on group). | Type-introspection on factor `axis` metadata | distinct `axis` values across the 3 factors per summand |
| A-1 | `test_apply_decomposition_invariant_slab` | same | `TestAlgebraDecomposition` | foundation | Slab: `op.apply(ψ).values + sigma_t · ψ.values == op.M_spatial.apply(ψ).values + op.M_angular_redist.apply(ψ).values`. **Defines the operator algebra contract.** | Algebraic identity (Resolution A) | `rtol=1e-13, atol=1e-14` (Route A bit-id if reductions don't reorder) |
| A-2 | `test_apply_decomposition_invariant_sphere` | same | same | foundation | Sphere variant of A-1; ALSO verifies that the bespoke `M_angular_redist` leaf produces the same contribution as the inline `pole_angular_closure.cell_contribution` summed over cells | Algebraic identity | same |
| A-3 | `test_apply_decomposition_invariant_cylinder` | same | same | foundation | Cylinder variant — exercises `level_indices` (multi-level recurrence) | Algebraic identity | same |
| A-4 | `test_LpC_apply_unchanged_post_t4_slab` | same | `TestLpCAlgebraInvariant` | foundation | `(L + C).apply(ψ)` value equals pre-T.4 snapshot bit-exact (the algebra-decomposition contract: subtracting σ_t·ψ at the apply boundary and re-adding at the C operator level returns the original M(ψ;σ_t)) | Pre-T.4 snapshot | `np.testing.assert_array_equal` (Route A) — see §6 R3 |
| A-5 | `test_LpC_apply_unchanged_post_t4_curvilinear` | same | same | foundation | Same as A-4 for sphere + cylinder | Pre-T.4 snapshot | same |
| A-6 | `test_LpC_solve_unchanged_post_t4` | same | same | foundation | `(L + C).solve(q)` value equals pre-T.4 snapshot bit-exact for all 3 geometries — the Q5 contract that `.solve` was not touched | Pre-T.4 snapshot + Q5 contract | `np.testing.assert_array_equal` |
| A-7 | `test_M_spatial_linearity_in_psi` | same | `TestAlgebraicIdentities` | foundation | `M_spatial.apply(α·ψ + β·φ) ≈ α·M_spatial.apply(ψ) + β·M_spatial.apply(φ)` at FP precision | Algebraic identity (structural) | `rtol=1e-12, atol=1e-13` |
| A-8 | `test_M_angular_redist_linearity` | same | same | foundation | Same for `M_angular_redist` (sphere + cylinder); identically zero for slab | Algebraic identity | same |
| A-9 | `test_M_spatial_zero_input_zero_output` | same | same | foundation | `M_spatial.apply(zeros).values == zeros` — linearity zero-guard | Structural (linearity) | `np.testing.assert_array_equal` |
| L4-1 | `test_apply_slab_bit_identical_to_pre_t4` | same | `TestPreT4RegressionSnapshot` | foundation | Slab 1G + 2G arms: per-dispatch-arm regression snapshot vs `pre_t4_snapshots.npz` | Pre-T.4 captured numerics | **Route A** (preserve reduction order): `np.testing.assert_array_equal`. **Route B** (sum-of-summands reorder): `assert_array_almost_equal_nulp(nulp=4*N_summands)` per principled-equivalence gate |
| L4-2 | `test_apply_sphere_bit_identical_to_pre_t4` | same | same | foundation | Sphere arm — both 1G and 2G | Pre-T.4 captured numerics | same |
| L4-3 | `test_apply_cylinder_bit_identical_to_pre_t4` | same | same | foundation | Cylinder arm — both 1G and 2G | Pre-T.4 captured numerics | same |
| L4-4 | `test_apply_2d_cartesian_bit_identical_to_pre_t4` | same | same | foundation | 2-D Cartesian (specular AND vacuum); strict bit-identity (Q1 decision: 2-D path untouched in T.4 scope) | Pre-T.4 captured numerics | `np.testing.assert_array_equal` (strict — no reductions reorder because the path is untouched) |
| L4-5 | `test_LpC_composition_apply_unchanged` | same | same | foundation | `(L + C).apply` composition is bit-identical pre/post-T.4 across all 3 geometries (Q3 + Q5: the subtraction unwinds; the composite output is the original M) | Pre-T.4 snapshot | `np.testing.assert_array_equal` |
| L4-6 | `test_InvertibleOperator_solve_unchanged` | `tests/sn/test_invertible_operator.py` (extend) | (existing class) | foundation | `(L + C).solve(q)` bit-identical pre/post-T.4 (Q5 contract) | Pre-T.4 snapshot | `np.testing.assert_array_equal` |
| L4-7 | `test_boundary_face_residuals_unchanged` | `tests/sn/test_streaming_operator.py` | `TestPreT4RegressionSnapshot` | foundation | `op.apply(ψ).boundary.face_view("xmax")` and `face_view("xmin")` (for slab) bit-identical pre/post-T.4. Verifies face-residual semantics (per plan §6 T.4 paragraph 1: "L is the ONLY operator that emits a non-zero face residual on its output .boundary") | Pre-T.4 snapshot | Route A: `np.testing.assert_array_equal`; Route B: nulp gate |
| A2D-1 | `test_apply_2d_cartesian_path_untouched` | same | same | foundation | Structural test: `_apply_2d_cartesian` source-code hash unchanged (or: `inspect.getsource(StreamingOperator._apply_2d_cartesian) == saved_signature`). Defensive pin against author drift. | Source-introspection snapshot | hash equality |
| L5-1 | `test_wave_t_t4_apply_overhead_below_5pct` | `tests/sn/performance/test_wave_t_overhead.py` | `TestT4PerformanceRegression` | (not vv-tagged) | Per plan §5.4. `(L+C).apply` walltime delta on 1-D slab P1 ≥1000 iterations. Pre-T.4 baseline from `pre_t4_walltime.json`. | Pre-T.4 wallclock baseline | `post_t4_walltime / pre_t4_baseline <= 1.05` |
| L1-1 | (delegated to existing) | `tests/sn/test_mms_aniso.py` | (module-level) | l1 | 1-D slab P1 aniso MMS O(h²) — must stay green | MMS pillar (structurally independent) | Existing xfail set unchanged |
| L1-2 | (delegated to existing) | `tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py` | (existing) | l1 | Curvilinear P1 aniso MMS — same xfail set | MMS pillar | Same xfail set |
| L1-3 | (delegated to existing) | `tests/sn/test_mms.py`, `tests/sn/test_mms_curvilinear.py`, `tests/sn/test_mms_heterogeneous.py`, `tests/sn/test_mms_2d.py` | (existing) | l1 | The full MMS suite (P0, multi-region, 2D) — must stay green | MMS pillar | All passing tests still pass |
| L1-4 | (delegated to existing) | `tests/sn/l1_analytical/test_kinf_homogeneous.py` | (existing) | l1 | `k_∞ = νΣ_f/Σ_a` closed-form on homogeneous reflective — slab AND curvilinear. **The structurally-independent closed-form reference for the principled-equivalence three-criteria gate.** | Closed-form pillar (homogeneous infinite medium) | Existing tolerance |
| L1-5 | (delegated to existing) | `tests/sn/test_heterogeneous_transport.py` | (existing) | l2 | Heterogeneous mesh-refined eigenvalue convergence — catches Mode #1 sign flip + Mode #4 wrong recursion (per `vv-principles` table) | L2 self-convergence | Existing tolerance |

### Test count summary

- **New foundation tests (algebraic-identity + structural)**: 18.
- **New regression-snapshot tests (L4)**: 7 (one per dispatch arm × group count).
- **New performance gate (L5)**: 1.
- **New defensive pin (A2D-1)**: 1.
- **Existing L1 gates that must stay green**: 5 test files.
- **Total new test files**: 0 (extend `test_streaming_operator.py` + extend `test_invertible_operator.py`).
- **Pre-T.4 snapshot fixtures required**: 2 files under `tests/sn/_fixtures/wave_t_t4/` (snapshots + walltime).

---

## §5 L1 MMS gates — explicit reference table

Per plan §5.2 and the test-architect's pillar-gate discipline.

| Test file | Geometry | Pillar | Claim layer | Status pre-T.4 | T.4 contract |
|-----------|----------|--------|-------------|----------------|--------------|
| `tests/sn/test_mms_aniso.py::test_sn_p1_aniso_mms_converges_second_order` | 1-D slab | MMS | Convergence-order + flux-shape | Green | Must stay green |
| `tests/sn/test_mms_aniso.py::test_sn_p1_aniso_mms_source_degrades_to_p0` | 1-D slab | MMS (degeneracy) | Source consistency | Green | Must stay green |
| `tests/sn/test_mms.py` (all) | 1-D slab P0 | MMS | Convergence-order | Green | Must stay green |
| `tests/sn/test_mms_curvilinear.py` (all) | 1-D sphere + cylinder P0 | MMS | Convergence-order | Green | Must stay green |
| `tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py` (all) | 1-D sphere + cylinder P1 | MMS | Convergence-order | 10 pre-existing DD-regression failures (per plan §6 line 295) | Same failure set; NO new failures |
| `tests/sn/test_mms_heterogeneous.py` (all) | Multi-region | MMS | Convergence-order + heterogeneity | Green | Must stay green |
| `tests/sn/test_mms_2d.py` (all) | 2-D Cartesian | MMS | Convergence-order | Green | Must stay green |
| `tests/sn/l1_analytical/test_kinf_homogeneous.py` (all) | Homogeneous reflective | Closed-form | Eigenvalue (k_∞) | Green | Must stay green |
| `tests/sn/l1_analytical/test_kinf_homogeneous_tolerance.py` | Homogeneous reflective | Closed-form | Eigenvalue tolerance | Green | Must stay green |
| `tests/sn/test_heterogeneous_transport.py` (all) | Heterogeneous | L2 self-convergence | Mesh refinement | Green | Must stay green |

**Why heterogeneous + mesh-refinement mandatory.** Per
`vv-principles` H2 ("Homogeneous eigenvalue is degenerate to
redistribution"): flat ψ makes every redistribution term identically
zero (`angular_denom_term`, `angular_numer_upstream` both zero). The
1-D curvilinear M_angular_redist is invisible to homogeneous tests.
ERR-006 (curvilinear sweep divergence) survived 20 passing tests
exactly because of this degeneracy.

**Why the L1 gate is GROUND TRUTH (not L4-1..L4-7).** The pre-T.4
snapshot regression tests (L4-1..L4-7) are bit-identity against the
OLD code. They cannot catch a bug that was ALREADY in the old code
and survived. The L1 MMS gate (against manufactured source —
structurally independent of any code path) AND the L1 k_∞ closed-form
gate (analytical reference) are the irreducible structural-independence
ground. If snapshot tests pass AND L1 fails, the snapshot itself was
wrong (both old and new code agree on a wrong value — a coincidence
this kind of test cannot detect).

---

## §6 Risk register — T.4-specific

Cross-referenced from plan §7 (the global Wave T register) and the
`numerical-bug-signatures` recognition catalogue. Risks ranked by
severity × catching-test coverage.

### R1 (HIGHEST) — Curvilinear angular redistribution lift breaks k_∞ agreement

**Mechanism.** The bespoke `AngularRedistributionOperator` leaf
re-extracts the M-M `cell_contribution` algebra from
`transport_operator_matvec_unified`'s nested `_sweep_direction` closure
(`sn/operator.py:442-502`). A subtle drift in which level_indices are
iterated, which face_in seed is used, or how the WDD `psi_face_in =
2.0 * psi_cell - psi_face_in` update interacts with the per-level
loop, breaks the recurrence. The bug fingerprint is **k_inf drift on
homogeneous reflective curvilinear** (sphere or cylinder) — a
`vv-principles` H2 violation (homogeneous degeneracy) AND a Signature-1
recurrence-coupled curvilinear bug (`numerical-bug-signatures`
Signature 1; ERR-006/ERR-007/ERR-026 family).

**Caught by.** L1-4 k_∞ homogeneous on sphere + cylinder at substep
T.4c. If k_∞ drifts even slightly, the regression snapshot tests
(L4-2, L4-3) ALSO fail (they share the converged eigenvector). The
TWO independent failures (closed-form + snapshot) together pinpoint
this as a real solver bug, not a reduction-reorder ULP artefact.

**Mitigation.**
- Gate the curvilinear lift behind L1-4 + L4-2/3 BEFORE merging T.4c.
- The bespoke leaf MUST call the same `pole_angular_closure.cell_contribution`
  primitive that `transport_operator_matvec_unified` calls — single
  source of truth, Pattern 2.
- Order T.4 substeps: slab (T.4b) before curvilinear (T.4c) — slab
  passes by construction (`M_angular_redist = ZeroOperator`) and
  validates the algebra-decomposition shape before introducing the
  bespoke leaf.

### R2 — 2-D Cartesian face_view-as-trace dependency surfaces

**Mechanism.** If the T.4 author misreads the Q1 scope decision and
attempts to lift `_apply_2d_cartesian` into `L_spatial_x +
L_spatial_y`, the cell-centre-proxy semantics break (per the
existing comment at `sn/operator.py:862-868`: 10% k_inf drift). The
2-D `psi.boundary.face_view` is currently passive; promoting it to a
real trace requires a BC realizer rewire that is OUT OF SCOPE for
T.4.

**Caught by.**
- L4-4: 2-D Cartesian bit-identity test fails immediately on the
  specular fixture (large drift, not ULP).
- A2D-1: source-code hash of `_apply_2d_cartesian` fails first
  (defensive pin).
- L1-3 (`test_mms_2d.py`) fails the manufactured-source convergence.

**Mitigation.** Q1 decision is binding: T.4 does NOT touch
`_apply_2d_cartesian`. The A2D-1 defensive pin catches accidental
modifications. The 2-D lift becomes its own follow-up wave.

### R3 — Reduction-order drift exceeds principled-equivalence ULP bound (HIGH for T.4)

**Mechanism.** This is the MA-Q5 concern. The slab WDD sweep
iterates sequentially over `nx` cells per ordinate, accumulating
`psi_face_in = 2.0 * psi_cell - psi_face_in` (WDD recurrence,
sequential). The reduction-tree depth at the matvec level is `O(nx)`
— for typical `nx=40`, drift bounded by `40 × 2.2e-16 ≈ 9e-15`. **For
typical eigenvalue Krylov iteration counts (~50), the cumulative drift
over the inner Krylov path is `~50 × 40 × ULP ≈ 4e-13`.** This is
ABOVE the typical `rtol=1e-13` regression contract for `apply` AND
near (but below) typical `rtol=1e-12` for inner Krylov solves.

**Resolution.** **Route A by default**: the spatial-axis sweep in
T.4b's `M_spatial.apply` MUST preserve the existing WDD recurrence
order verbatim. The lift is structural (extract from the closure
into a `TensorProductOperator` factor + summand) NOT a reduction-
order rewrite. With Route A, L4-1/L4-2/L4-3 pass strict
`np.testing.assert_array_equal`.

**Route B fallback** (per `vv-principles` §"Bit-identity vs
principled-equivalence"): IF a Route A is structurally infeasible
(e.g., the `TensorProductOperator.apply` fold inserts a `np.add` at
a different stage than the legacy fused einsum), apply the
three-criteria gate:

1. **Principled at every step**: each summand `D_axis & Ω_axis &
   I_g` IS the per-axis streaming contribution at the cell-balance
   denom-numer level — a named primitive (`cell_balance_for_streaming`
   already exposes the three named denom-numer terms). Not an
   unnamed intermediate.
2. **Structurally-independent reference**: L1-4 k_∞ homogeneous
   closed-form `k_∞ = νΣ_f/Σ_a`. The k_∞ value is independent of the
   spatial reduction order.
3. **FP-non-associativity, dimensionally explainable**: drift bounded
   by `iter_count × nx × ULP`. For 1000-matvec benchmark walltime
   measurement: drift bounded; for single-matvec snapshot tests:
   drift bounded by `nx × ULP ≈ 1e-14`.

**Pass criterion for Route B**: `assert_array_almost_equal_nulp(nulp =
nx_max × n_summands)` for snapshot tests, where `nx_max = 40` and
`n_summands = 1` (slab) or `3` (curvilinear: spatial summand +
angular summand + the apply-boundary subtraction reordering). For
`(L+C).solve` tests: `rtol = max(1e-12, iter_count × nx × ULP)` with
`iter_count` captured at solve time.

**Action.** L4-1..L4-7 tests parameterise the assertion gate by
`ROUTE`. Set `ROUTE=A` as default; flip to `ROUTE=B` only if T.4b's
gate fails strict equality and the three-criteria check passes.

### R4 — Q3 subtraction location drift breaks `(L+C).apply` bit-identity

**Mechanism.** If the T.4 author moves the `σ_t · ψ` subtraction
INTO `M_spatial.apply` (instead of at the apply boundary per the Q3
recommendation), the algebra `(L+C).apply = M_spatial + M_angular_redist
+ σ_t·ψ` no longer unwinds correctly. The `(L+C).apply` value drifts
because the subtraction-then-addition pair now happens at a different
position in the reduction tree.

**Caught by.**
- L4-5 (`test_LpC_composition_apply_unchanged`) fails at FP-non-
  associativity at minimum, or at structural level if the subtraction
  is now in the WRONG place algebraically.
- A-1/A-2/A-3 (`test_apply_decomposition_invariant_*`) defines the
  contract; if the subtraction lives inside `M_spatial`, the
  decomposition identity fails by `σ_t · ψ` (large drift, not ULP).

**Mitigation.** Q3 decision is binding: the subtraction lives at the
apply boundary, OUTSIDE `M_spatial` and `M_angular_redist`. The
properties produce `M_spatial(ψ)` and `M_angular_redist(ψ)` — the
un-subtracted contributions per `cell_balance_for_streaming`'s
additive structure.

### R5 — Q5 InvertibleOperator path accidentally touched

**Mechanism.** The T.4 author refactors `StreamingOperator` and
inadvertently changes `__add__` dispatch behaviour OR modifies
`InvertibleOperator.solve` (e.g., to route through the new
`M_spatial.apply`). The procedural WDD sweep at `transport_sweep`
breaks.

**Caught by.**
- L4-6 (`test_InvertibleOperator_solve_unchanged`) fails strict
  bit-identity against pre-T.4 snapshot.
- A-4/A-5/A-6 (`test_LpC_apply/solve_unchanged_post_t4_*`) cover
  both arms.

**Mitigation.** Q5 contract is binding. The L4-6 gate is enforced
at every substep that touches `operator.py` (i.e., all of T.4b–T.4e).

### R6 — Group-axis convention drift (ERR-002 family)

**Mechanism.** Per `numerical-bug-signatures` Signature 3, the
project's `Σ_t` and `SigS` conventions are group-leading (`(ng, nx,
ny)`). The bespoke `AngularRedistributionOperator` leaf MUST consume
sigma_t and total_xs at the right shape. A naive transpose drift in
the leaf's body (using `sigma_t[..., None]` instead of
`sigma_t[None, :, :, :]`, etc.) silently broadcasts wrong and
produces a 2G-asymmetric divergence.

**Caught by.**
- L4-2 / L4-3 fixtures use 2G ASYMMETRIC SigS (per §3 schema).
- L1-2 curvilinear MMS exercises 2G P1 anisotropic — the load-bearing
  group-coupling probe.

**Mitigation.** The bespoke leaf's body MUST call the SAME
`pole_angular_closure.cell_contribution` primitive that
`transport_operator_matvec_unified` calls — Pattern 2 single source
of truth. No re-implementation of the per-cell algebra in the leaf.

### R7 — Slab degeneracy hides M_angular_redist bugs (the silent failure)

**Mechanism.** If T.4c introduces a bug in `M_angular_redist` (e.g.,
wrong sign on `angular_denom_term`) that vanishes for ZeroOperator
(slab) but is active for curvilinear, the slab gate (T.4b) passes
green and the test-design loop falsely concludes the abstraction is
working. ERR-006 lived exactly here: homogeneous-reflective tests
passed because `angular_denom_term = (ΔA/w) · c_out` is structurally
zero at flat ψ.

**Caught by.**
- L1-2 curvilinear MMS aniso (NOT homogeneous) — the
  manufactured-source ψ is non-flat by construction.
- L1-5 heterogeneous mesh-refined (`test_heterogeneous_transport.py`)
  — flat ψ is forbidden by the heterogeneity.
- The `vv-principles` H2 cardinal rule: never accept homogeneous-only
  verification.

**Mitigation.** Mandatory mesh-refined heterogeneous + P1 anisotropic
test coverage at substep T.4c (NOT just T.4b's slab). The §3 snapshot
fixture explicitly includes 2G P1 asymmetric SigS at curvilinear
geometries.

### R8 — Performance regression > 5% (the perf-tax)

**Mechanism.** Per plan §7 risk 2: SOTP.apply is a Python-level fold
over summands; each summand is a 3-factor TP fold. For slab (1
summand × 3 factors): 3 Python calls per matvec vs ~1 numpy
operation pre-T.4. For curvilinear: same +
`AngularRedistributionOperator.apply` which iterates per cell.

**Caught by.** L5-1 perf gate at substep T.4f.

**Mitigation.**
- `cached_property` on `M_spatial` and `M_angular_redist` so
  construction happens once per StreamingOperator instance, not per
  apply.
- Factor caching at `TensorProductOperator._build` (already exists
  for the apply path; `solve` fold also caches per
  `numerics/operator.py:1058+`).
- If exceeded, investigate `(D_x & Ω_x & I_g)` factor fusion — the
  three factors might fuse at construction time into a single
  pre-cached apply closure. Deferred to a follow-up wave per plan
  §4.

---

## §7 Risk-to-test cross-map

For traceability — every risk above maps to at least one test row
in §4.

| Risk | Catching test(s) | Substep gate fires at |
|------|------------------|----------------------|
| R1 (k_∞ curvilinear) | L1-4 + L4-2/L4-3 (two independent failures) | T.4c |
| R2 (2-D face_view) | L4-4 + A2D-1 + L1-3 | T.4d |
| R3 (reduction drift) | L4-1..L4-7 with ROUTE flag | T.4b (slab first) |
| R4 (subtraction location) | L4-5 + A-1/A-2/A-3 | T.4b + T.4c |
| R5 (InvertibleOperator) | L4-6 + A-6 | All of T.4b–T.4e |
| R6 (group convention) | L4-2/L4-3 (2G asymmetric) + L1-2 | T.4c |
| R7 (slab hides bugs) | L1-2 curvilinear aniso MMS + L1-5 heterogeneous | T.4c (NOT T.4b) |
| R8 (perf regression) | L5-1 | T.4f |

---

## §8 Pickup checklist (main agent direct authorship)

Per `feedback_no_method_implementer_for_surgical_carves`: T.4 is the
main agent's surgical work with turn-by-turn user steering. NO
method-implementer dispatch.

### Files to touch (T.4b–T.4c surgery)

1. **`orpheus/sn/operator.py`** (the primary surgery target):
   - L644-969 `StreamingOperator` class — add `M_spatial`,
     `M_angular_redist` cached_properties; rewire `apply` for 1-D
     path to delegate.
   - **DO NOT touch** L830-954 `_apply_2d_cartesian` (Q1 contract).
   - **DO NOT touch** L958-969 `__add__` (Q5 contract).
   - **DO NOT touch** `InvertibleOperator.solve` /
     `_solve_timed_full_field` / `transport_sweep` (Q5 contract).
   - L203-642 `transport_operator_matvec_unified` — may STAY as the
     procedural body during the substep; eventually the body is
     called by `M_spatial.apply` + `M_angular_redist.apply` in
     two separate invocations. **Final state**: this function MAY
     become a thin shim that composes the two new property apply
     calls, OR it may be retired if every consumer has migrated.
     Retirement is out of T.4 scope; defer to T.5 cleanup.

2. **NEW (T.4c)**: bespoke leaf type
   `AngularRedistributionOperator` adjacent to `StreamingOperator`
   in `operator.py` (consistent with the file's existing leaf-class
   layout). Its `apply` body invokes
   `pole_angular_closure.cell_contribution` per cell — Pattern 2
   single source of truth.

3. **DO NOT touch**:
   - `orpheus/geometry/reduced_operator.py:151-511` (`StreamingTerms`
     + `ReducedStreamingOperator`) — geometric primitives stay.
   - `orpheus/sn/spatial/cell_balance.py:120-198`
     (`cell_balance_for_streaming`) — single source of truth, stays.
   - `orpheus/sn/spatial/pole_angular_closure.py:191-1200`
     (`PoleAngularClosure` Protocol + `MorelMontryAngularSweep` +
     `IdentityAngularClosure`) — strategy primitives stay.
   - `orpheus/numerics/operator.py:1058-1269` (`TensorProductOperator`
     + `SumOfTensorProductsOperator` + `assert_separable`) — L1
     primitives stay.

### Files to extend (new tests)

1. **`tests/sn/test_streaming_operator.py`** — extend with new
   classes:
   - `TestTensorNetworkDecomposition` (L2-1..L2-7)
   - `TestAlgebraDecomposition` (A-1..A-3)
   - `TestLpCAlgebraInvariant` (A-4..A-6)
   - `TestAlgebraicIdentities` (A-7..A-9)
   - `TestPreT4RegressionSnapshot` (L4-1..L4-5, L4-7, A2D-1)

2. **`tests/sn/test_invertible_operator.py`** — extend with L4-6.

3. **`tests/sn/performance/test_wave_t_overhead.py`** — extend (or
   create if no file) with `TestT4PerformanceRegression` (L5-1).

### Files to write (snapshot + perf fixtures)

1. **`tools/wave_t/capture_pre_t4_snapshots.py`** — one-shot script
   per §3 schema.
2. **`tests/sn/_fixtures/wave_t_t4/pre_t4_snapshots.npz`** —
   generated by step 1.
3. **`tests/sn/_fixtures/wave_t_t4/pre_t4_walltime.json`** —
   generated by step 1.

### Sequence (turn-by-turn)

1. **Main agent + user**: review this spec; confirm Q1-Q5 + MA-Q1-MA-Q5
   resolutions; flag any architectural objection BEFORE T.4a.
2. **Main agent**: write `capture_pre_t4_snapshots.py`; run once;
   commit fixtures (substep T.4a commit).
3. **Main agent direct authorship**: execute T.4b (slab `M_spatial`).
   Run L2-1..L2-7 (slab subset) + L4-1 (slab) + L1-1 + L1-4 (slab).
   Commit IFF all gate tests pass.
4. **Main agent direct authorship**: execute T.4c (curvilinear
   `M_angular_redist`). Run L2-1..L2-7 (curvilinear subset) +
   L4-2/L4-3 + L1-2 + L1-4 (curvilinear). Commit IFF all gate tests
   pass.
5. **Main agent**: execute T.4d (2-D regression pin only — no code).
   Run L4-4 + A2D-1. Commit IFF pass.
6. **Main agent**: execute T.4e (full L1 sweep + algebra invariants).
   Run A-1..A-9 + L4-5/L4-6/L4-7 + L1-1..L1-5. Commit IFF pass.
7. **Main agent**: execute T.4f (perf gate). Run L5-1. Commit IFF
   pass.
8. **Main agent**: dispatch **qa** sub-agent for post-T.4 review
   per `subagent-handoff-protocol`. The qa audit covers:
   - Pattern 7 producer-side normalisation integrity (the σ_t · ψ
     subtraction at the apply boundary).
   - Pattern 2 single source of truth (the bespoke leaf calls
     `pole_angular_closure.cell_contribution`; does not
     re-implement).
   - Pattern 4 illegal-states-unrepresentable (slab's
     `M_angular_redist = ZeroOperator` is structural, not a runtime
     branch).
   - ERR-002 transpose convention preservation in the lifted form.
   - No drift in `InvertibleOperator.solve` path.

---

## §9 Architectural concerns to flag BEFORE T.4 starts

### Concern 1 — Plan §6 T.4 says `R_polar(connection_coefficients) & I_x & I_group`, which is option (i) per Q2

The plan as written suggests the 3-factor wrap is the target form.
My Q2 recommendation (option iii) deviates: `AngularRedistributionOperator`
as a bespoke leaf, NOT wrapped as `(leaf & I_x & I_g)`. The wrap is
cosmetic and false-asserts separability the M-M recurrence does not
support.

**Recommendation to user**: confirm the plan §6 T.4 form is
ASPIRATIONAL (matches §15 grand-report aspiration) but NOT
load-bearing. The plan §7 risk 1 mitigation explicitly allows the
fallback: "ensure the algebra still composes via OperatorSum. The
architectural commitment is ONE algebraic FORM, not that every
factor is a clean TensorProductOperator." If the user agrees, the
plan §6 T.4 wording should be UPDATED to reflect the option (iii)
fallback as the production target.

### Concern 2 — The "M_spatial / M_angular_redist" naming deviates from plan

The plan §6 T.4 names them `L_spatial` and `L_angular_redist`. My Q3
recommendation renames to `M_spatial` and `M_angular_redist` because
they produce the UN-subtracted (L+C-share) contributions per the
cell-balance algebra. The discrete primitives are `M`-shaped, not
`L`-shaped; the apply boundary's subtraction recovers L.

**Recommendation to user**: confirm the rename. The continuous-vs-
discrete distinction is load-bearing for clarity; `L_spatial` is
ambiguous about whether it includes the cell-collision share.

### Concern 3 — `transport_operator_matvec_unified` deprecation timing

After T.4c, the function `transport_operator_matvec_unified` becomes
a thin orchestrator that the new property apply paths call. **Does
it survive T.4, or is it retired in T.5?** Per
`feedback_aggressive_retirement` and Cardinal Rule 2: superseded
code is noise. If `M_spatial.apply` + `M_angular_redist.apply` cover
the full matvec, the unified function is dead weight.

**Recommendation to user**: explicitly schedule the retirement of
`transport_operator_matvec_unified` in T.5 (or after T.4c if it
makes sense). The spec does NOT include the retirement in T.4 scope
to keep T.4 commits focused, but the test architecture supports the
retirement: every existing test that calls
`transport_operator_matvec_unified` directly (none in the production
path; only test code) would migrate to calling
`StreamingOperator.apply` instead.

### Concern 4 — The `assert_separable` invariant for `AngularRedistributionOperator`

The bespoke leaf is NOT a tensor product. Calling `assert_separable`
on `M = M_spatial + M_angular_redist` (an OperatorSum) requires that
`assert_separable` is defined for `OperatorSum` and either:
- (a) succeeds if ALL summands are separable, OR
- (b) raises if ANY summand is non-separable.

If (b), the slab `M = M_spatial + ZeroOperator` would pass (both are
separable trivially), but curvilinear `M = M_spatial +
AngularRedistributionOperator` would FAIL because the second summand
is NOT a tensor product. The L2-5 test as written passes
`M_spatial.assert_separable()` alone — fine for slab and curvilinear
both.

**Recommendation to user**: the test asserts `M_spatial.assert_separable()`
only, NOT the full `(L+C)` operator. The non-separable
`M_angular_redist` is acknowledged at the type level (L2-4) and the
operator-algebra accepts heterogeneous summands per plan §7 risk 1
mitigation.

---

## §10 Test-architect self-assessment

### Pillar gate (per `vv-principles` §1.5)

| Test row | Claim layer | Pillar | Structural independence? |
|----------|-------------|--------|--------------------------|
| L2-1..7 | Software invariant (type/structure) | n/a (foundation) | Pure type introspection |
| A-1..9 | Algebraic identity | Closed-form (linearity, additive decomposition) | Structural — pure math |
| L4-1..7 | Bit-identity (or principled-equivalence) to pre-T.4 | Regression-snapshot | Pre-T.4 numerics ARE the reference; pillar #2 (L1-4 k_∞ closed-form) backs the nulp gate |
| A2D-1 | Source-code invariance | Defensive structural | Source-hash identity |
| L5-1 | Performance (not vv-tagged) | n/a | n/a |
| L1-1..3 | Convergence + flux-shape | MMS | MMS is the pillar; manufactured source structurally independent of solver primitives |
| L1-4 | Eigenvalue | Closed-form | `k_∞ = νΣ_f/Σ_a` — homogeneous infinite medium, exact |
| L1-5 | Mesh refinement | L2 self-convergence + heterogeneity | Cross-mesh ratio |

**MMS does NOT prove eigenvalues** — confirmed in row L1-1..3
(convergence + flux-shape claims). L1-4 is the EIGENVALUE claim via
closed-form, the only pillar that supports it.

### Anti-pattern audit (per `vv-principles` §0)

1. NOT claiming verification on L4 alone — L1 + L4 stack, with L1-4
   closed-form k_∞ as the independent reference for the
   principled-equivalence Route B nulp gate.
2. NOT asserting `np.allclose` against another solver — comparing to
   pre-T.4 captured numerics (snapshot) + closed-form k_∞ + MMS
   exact source.
3. NOT accepting 1-group as verification — 2G test data mandated in
   §3 snapshot schema (all `*_2g_*` entries) AND L1-1/L1-2/L1-4 use
   ≥2G.
4. NOT homogeneous-only — L1-2 curvilinear aniso MMS provides
   non-flat ψ; L1-5 heterogeneous mesh-refined adds material
   heterogeneity; §3 fixture uses asymmetric SigS for group
   coupling.
5. NOT trusting convergence rate alone — L1-1/L1-2 verify the
   converged-to value via MMS exact solution; L1-4 verifies the
   converged-to value via closed-form k_∞.
6. NOT trusting an untraced reference — pre-T.4 snapshot is
   regression/control-on-self, L1 chain terminates in closed-form +
   MMS exact solution (both structurally independent).
7. NOT trusting "two derivations agree" — pre-T.4 vs post-T.4 is
   bit-identity at the IMPLEMENTATION level (Route A) or
   nulp-bounded (Route B) BACKED BY independent closed-form (L1-4
   pillar #2 evidence).
8. NOT particle balance as L0 — per-summand independence (A-1..3) +
   decomposition invariant replaces telescoping-sum traps.
9. NOT conflating validation with verification — T.4 is a pure
   refactor; no equation change; pure V&V.
10. NOT "reasonable numbers" — every algebraic identity verified
    structurally; every snapshot regression verified to ULP or
    near-ULP.

### Multi-group + heterogeneous + mesh-refined?

- **Multi-group**: §3 snapshot fixtures use 2G asymmetric SigS for
  all geometries (R6 mitigation).
- **Heterogeneous**: L1-5 (`test_heterogeneous_transport.py`) is
  the load-bearing heterogeneous gate; §3 fixtures also include
  multi-material in some captures.
- **Mesh-refined**: L1-1, L1-2, L1-3 are convergence ladders.

All three mandatory dimensions covered.

### New failure modes introduced — self-improvement directive

Per `subagent-handoff-protocol` and the test-architect self-improvement
trigger: T.4 introduces ONE new failure-mode pattern not yet in the
`vv-principles` §6 anti-pattern table OR the `numerical-bug-signatures`
recognition catalog:

**Pattern: structural-pin defensive test (A2D-1)** — asserting that
the source code of a method is unchanged via hash. This is a NEW
test class — not a failure mode per se, but a defensive technique
that pre-emptively catches author drift in the multi-substep workflow
of `feedback_no_method_implementer_for_surgical_carves`. It is the
test-architect's first use of source-introspection as a regression
gate.

**Recommendation**: do NOT add this as an anti-pattern entry; it is
a positive defensive pattern. Mention in the T.5 close-out memo if
it proves load-bearing.

---

## Pointers

- **Plan**: `.claude/plans/wave_t_tensor_network.md` §6 step T.4 +
  §1.2 (two views of the connection-coefficient operator) + §1.3
  (expected complications) + §5 verification strategy + §7 risk
  register row 1.
- **T.3 precedent spec**:
  `.claude/agent-memory/test-architect/wave_t_t3_scattering_verification_spec.md`
  — the document this spec is modelled on.
- **Grand report**: §15 lines 2003-2019 (V = X ⊗ Ω ⊗ G);
  §15.1 lines 2031-2045 (streaming as sum of tensor products); §35
  line 5675 (commandment); north-star line 5697.
- **Production code (primary surgery target)**:
  - `orpheus/sn/operator.py:644-969` (`StreamingOperator`).
  - `orpheus/sn/operator.py:203-642`
    (`transport_operator_matvec_unified` — geometry-agnostic 1-D
    body).
  - `orpheus/sn/operator.py:830-954` (`_apply_2d_cartesian` — DO NOT
    TOUCH per Q1).
- **Production code (READ but DO NOT MODIFY)**:
  - `orpheus/geometry/reduced_operator.py:151-511` (`StreamingTerms`
    + `ReducedStreamingOperator`).
  - `orpheus/sn/spatial/cell_balance.py:120-198`
    (`cell_balance_for_streaming` — single source of truth, Pattern
    2).
  - `orpheus/sn/spatial/pole_angular_closure.py:191-1200`
    (`PoleAngularClosure` + `MorelMontryAngularSweep` +
    `IdentityAngularClosure`).
- **L1 primitives**:
  - `orpheus/numerics/operator.py:1058-1177` (`TensorProductOperator`).
  - `orpheus/numerics/operator.py:1179-1269`
    (`SumOfTensorProductsOperator` + `assert_separable`).
  - `orpheus/numerics/space.py:301` (`TensorProductSpace`).
- **L1 MMS gates** (must stay green):
  - `tests/sn/test_mms_aniso.py` (1-D slab P1 aniso).
  - `tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py`
    (curvilinear P1, same xfail set).
  - `tests/sn/l1_analytical/test_kinf_homogeneous.py` (closed-form
    eigenvalue — the principled-equivalence pillar #2 reference).
  - `tests/sn/test_mms.py`, `test_mms_curvilinear.py`,
    `test_mms_heterogeneous.py`, `test_mms_2d.py`.
- **`vv-principles` §"Bit-identity vs principled-equivalence"**: the
  three-criteria gate for Route B (nulp-relaxation); every criterion
  answered in §1 Q3 + §6 R3.
- **`vv-principles` §"1-group degeneracy"**: justifies the 2G
  fixture mandate in §3.
- **`coding-elegance` Pattern 2** (single source of truth) — the
  bespoke leaf's body MUST call `pole_angular_closure.cell_contribution`,
  not re-implement it.
- **`coding-elegance` Pattern 4** (illegal states unrepresentable) —
  slab's `M_angular_redist = ZeroOperator` is structural at the type
  level, NOT a runtime branch.
- **`coding-elegance` Pattern 7** (producer-side normalisation) —
  the σ_t · ψ subtraction at the apply boundary (Q3 decision).
- **`numerical-bug-signatures` Signature 1**: curvilinear sweep
  divergence under refinement (ERR-006/ERR-007/ERR-026) — guarded
  by §6 R1 + §6 R7.
- **`numerical-bug-signatures` Signature 3**: ERR-002 SigS transpose
  convention — guarded by §6 R6 (2G asymmetric fixtures).
- **`feedback_unify_after_two_instances`** — Q1 decision rationale
  (3 working 1-D instances unify; 2-D is genuinely different).
- **`feedback_no_method_implementer_for_surgical_carves`** — T.4 is
  main-agent direct authorship.
- **`feedback_aggressive_retirement`** — §9 Concern 3 calls out the
  `transport_operator_matvec_unified` retirement question for T.5.
