---
name: issue-196-phase-g-replan-algebra
description: Issue #196 Phase G replan — exhaustive inventory of EXISTING four-operator algebra infrastructure (`LinearOperator` Protocol, composers, L/C/S/F leaves, SourceIteration/KEigenvalue) so the next plan does not reinvent any wheel.
metadata:
  type: project
---

# Issue #196 Phase G — operator-algebra infrastructure inventory

**Audit date**: 2026-05-13. **Branch**:
`refactor/sn-operator-algebra`. **Scope**: enumerate every existing
piece of the `(L + C - S - F/k) ψ = q` four-operator algebra so the
Phase G replan starts from "what's already there", not from a clean
sheet.

Use `(file:line)` citations to keep this verifiable. Reading order is
roughly bottom-up: protocol → composers → leaves → iteration →
narrative.

---

## TL;DR — one-line status per piece

| Component                            | Status            | Where                                                  |
| ------------------------------------ | ----------------- | ------------------------------------------------------ |
| `LinearOperator` Protocol            | exists fully      | `orpheus/numerics/operator.py:116-183`                 |
| `LinearOperatorMixin` (algebra dunders) | exists fully   | `orpheus/numerics/operator.py:202-407`                 |
| `OperatorSum`, `OperatorProduct`     | exists fully, `.H` propagates correctly | `operator.py:489-638` |
| `ScaledOperator`, `IdentityOperator`, `ZeroOperator` | exists fully | `operator.py:641-736`                          |
| `_AdjointOperator` (the `.H` wrapper)| exists fully      | `operator.py:415-486`                                  |
| `as_scipy_linop`                     | exists fully      | `operator.py:1318-1368`                                |
| `MissingCapability`, `IncompatibleOperatorComposition` | exists fully | `operator.py:93-112`                       |
| `L` — `SNStreamingOperator` (monolith: streaming+collision) | exists fully | `orpheus/sn/operator.py:1115-1472`     |
| `L` — `StreamingOperator` (pure-streaming, Phase G Step 2 first half) | exists partially | `orpheus/sn/streaming.py:251-408` |
| `C` — `CollisionOperator` (Σ_t alone) | exists fully     | `orpheus/sn/streaming.py:417-512`                      |
| `S` — `ScatteringOperator`           | exists fully (P0+Pℓ+n2n; `R·Λ·M` already composed) | `orpheus/sn/scattering.py:254-524` |
| `F` — `FissionOperator`              | exists fully (rank-1 `χ ⊗ νΣ_f`) | `orpheus/sn/fission.py:82-165`         |
| `SourceIteration(L, S, F, q_ext)`    | exists fully — already the Step 3 shape | `orpheus/numerics/iteration.py:164-352` |
| `KEigenvalue(L, S, F)`               | exists fully      | `orpheus/numerics/iteration.py:360-569`                |
| `power_iteration` (legacy)           | DEPRECATED on import | `orpheus/numerics/eigenvalue.py:120-156`            |
| `compute_keff` / `compute_group_production_rate` / `compute_group_absorption_rate` | exists fully | `orpheus/sn/solver.py:295-358` |
| Sphinx narrative `:ref:operator-algebra` | exists (Phase 0 stub, 796 lines) | `docs/theory/operator_algebra.rst`             |
| `HarmonicMomentProjection` (M), `HarmonicMomentReconstruction` (R), `LegendreMomentScattering` (Λ) | exists fully | `orpheus/numerics/projection.py:1-477`, `orpheus/sn/scattering.py:132-251` |
| `build_sn_operators(V, materials) → (L, C)` | exists, returns `(L, C)` only — S, F not yet wired in | `orpheus/sn/streaming.py:697-` |
| Builder for the full `(L, C, S, F)` 4-tuple from materials | DOES NOT EXIST yet | (only `(L, C)` returned at `streaming.py:697`) |
| `solve_within_group_sn_curvilinear(L, C, q)` (Picard/GMRES) | exists, only curvilinear, only `Q_aniso=None` for now | `orpheus/sn/streaming.py:520-688` |

---

## 1. The `LinearOperator` Protocol + algebra machinery

File: `orpheus/numerics/operator.py` (1368 lines, 24 May-2026).
Module docstring (`operator.py:1-51`) frames the algebra as the
`(L - S - F) ψ = q` form already, citing Trefethen & Bau §3.2 — the
mandate's "four-operator" form is one composition rewrite away.

### 1.1 Capability tag literals (`operator.py:88-90`)

```
CAP_APPLY:           "apply"
CAP_SOLVE:           "solve"
CAP_APPLY_TRANSPOSE: "apply_transpose"
```

These are bare strings (not an Enum) so user operators may advertise
method-specific tags without subclassing. Composers ignore tags they
don't recognise; only the three primitive ones propagate.

### 1.2 Exceptions

- `MissingCapability(TypeError)` (`operator.py:93-99`) — raised at
  composition time when a constituent lacks a needed capability.
  Composers veto incompatible compositions *at construction*, never
  at call time.
- `IncompatibleOperatorComposition(ValueError)`
  (`operator.py:102-112`) — domain/range mismatch in `+` or `@`.
  Skipped when either operand has `domain=None` / `range=None`
  (legacy backward-compat for ops predating Issue 9.6).

### 1.3 Protocol surface (`operator.py:115-183`)

```python
@runtime_checkable
class LinearOperator(Protocol):
    capabilities: frozenset[str]                # MANDATORY tag set
    @property
    def domain(self) -> Optional["FunctionSpace"]: ...   # may return None
    @property
    def range(self) -> Optional["FunctionSpace"]: ...    # may return None
    def apply(self, x: np.ndarray) -> np.ndarray: ...   # MANDATORY

    # Optional — declared in capabilities, not in the type:
    # def solve(self, b)
    # def apply_transpose(self, x)
```

Note: shape and dtype are deliberately NOT part of the protocol
(`operator.py:142-150`) — numpy duck-typing handles them at
`apply` call time. The Wave-A philosophy was "the contract should
be 'I can do *these* things' not 'I can do everything or fail'".

### 1.4 `LinearOperatorMixin` — every dunder it ships (`operator.py:202-407`)

Inheriting this gains the full algebra (no further effort):

| Dunder         | Returns                                                       | Line       |
| -------------- | ------------------------------------------------------------- | ---------- |
| `__add__`      | `OperatorSum(self, other)`                                    | 254-255    |
| `__radd__`     | `OperatorSum(other, self)`                                    | 257-258    |
| `__sub__`      | `OperatorSum(self, ScaledOperator(-1.0, other))`              | 260-261    |
| `__rsub__`     | `OperatorSum(other, ScaledOperator(-1.0, self))`              | 263-264    |
| `__mul__`      | `ScaledOperator(scalar, self)` (numeric scalar)               | 266-269    |
| `__rmul__`     | delegates to `__mul__`                                        | 271-272    |
| `__neg__`      | `ScaledOperator(-1.0, self)` — supports `-A`                  | 274-287    |
| `__truediv__`  | `ScaledOperator(1.0/scalar, self)` — supports `A / k`         | 289-304    |
| `__matmul__`   | `OperatorProduct(self, other)` — the `A @ B` form             | 306-307    |
| `__and__`/`__rand__` | `TensorProductOperator._build(...)` — `A & B` per Grand Report §6.3 | 309-324    |
| `__call__`     | alias for `apply` — supports `A(x)` math notation             | 326-333    |
| `__pow__`      | `A**n` for n ≥ 0; `A**0 = IdentityOperator()`                 | 335-360    |
| `adjoint()`    | returns `_AdjointOperator(self)`                              | 366-388    |
| `H` (property) | alias for `adjoint()` — supports `A.H` Hilbert-dagger         | 390-394    |
| `__repr__`     | uniform reporting with class, domain, range, capabilities     | 400-407    |

Default `domain`/`range` return `None` (`operator.py:242-248`) —
concrete operators may override.

### 1.5 `_AdjointOperator` — how `.H` propagates (`operator.py:415-486`)

When user calls `A.H` on any `LinearOperatorMixin`:

- Constructs `_AdjointOperator(A)` (`operator.py:432-442`).
- `_AdjointOperator.domain` returns `A.range`;
  `_AdjointOperator.range` returns `A.domain` — swapped (lines 445-451).
- Capability rule: the adjoint advertises `CAP_APPLY` **iff** the
  inner advertises `CAP_APPLY_TRANSPOSE` (line 436-437). `solve`
  does **not** propagate to the adjoint — the algebra of `A.H.solve`
  is deferred until a consumer demands it (line 440-442).
- `apply(y)` (line 453-479) performs the weight-aware Hilbert adjoint
  action `(A* y)_V = (1/w_V) · A.apply_transpose(w_W · y)` where
  `w_V`, `w_W` are the inner-product weights from `domain` and
  `range`. When both weight arrays are `None` (Euclidean inner
  product on both sides), this reduces to the bare
  `apply_transpose`.

Important: `_AdjointOperator.apply_transpose` raises
`NotImplementedError` (lines 481-486). If a consumer wants the
transpose of an adjoint, take the adjoint of the original's
transpose directly. This is fine for the four-operator algebra —
`(L+C-S-F/k).H = L.H + C.H - S.H - (1/k)·F.H`, which is a
composition of `.H` of leaves; `_AdjointOperator` never needs to
re-transpose itself.

### 1.6 `OperatorSum` (`operator.py:489-566`)

Constructor (`489-549`) does ALL of:

- `MissingCapability` if either operand lacks `CAP_APPLY` (513-522).
- `IncompatibleOperatorComposition` if domains or ranges
  mismatch (526-541) — skipped if either operand has `None`.
- Capability closure (545-549):
  - `apply` always (both summands have it by construction).
  - `apply_transpose` iff **both** have it
    (transposition distributes over sums: `(A+B)^T = A^T + B^T`).
  - `solve` does **NOT** propagate — there is no general algorithm
    for `(A+B)^{-1}` from `A^{-1}, B^{-1}` (Sherman-Morrison-Woodbury
    is the special case under low-rank structure).

`apply` (`561-562`):  `return self.a.apply(x) + self.b.apply(x)`.
`apply_transpose` (`564-565`): linearity.

### 1.7 `OperatorProduct` (`operator.py:568-638`)

Constructor (`587-617`):

- `MissingCapability` if either operand lacks `CAP_APPLY` (588-597).
- `IncompatibleOperatorComposition` if `A.domain != B.range`
  for `A @ B` (600-608) — skipped if either is `None`.
- Capability closure:
  - `apply` always.
  - `solve` iff **both** have `solve`, with order REVERSED:
    `(A B)^{-1} = B^{-1} A^{-1}` (613-614, 632-634).
  - `apply_transpose` iff **both** have it, with order REVERSED:
    `(A B)^T = B^T A^T` (615-616, 636-638).

`A @ B` returns shape: `domain = B.domain`, `range = A.range`.

### 1.8 `ScaledOperator` (`operator.py:641-694`)

- All inner capabilities preserved (unitary on the cap set).
- Zero scalar rejected (`657-665`): use `ZeroOperator` explicitly.
- `apply` accepts `*extra, **kwextra` (line 686) — load-bearing!
  This is how `ScaledOperator(1/k, F).apply(phi)` works even though
  `F.apply(phi)` has only one positional arg; and it tolerates the
  multi-arg `StreamingOperator.apply(psi, *, sigma_t=...)` signature
  (`streaming.py:310-314`).
- `solve` (`689-691`): `(αL)^{-1} q = (1/α) L^{-1} q` — the
  reciprocal applies on the way out.

### 1.9 `IdentityOperator` (`operator.py:697-716`)

`capabilities = frozenset({CAP_APPLY, CAP_SOLVE, CAP_APPLY_TRANSPOSE})`.
All three actions return `x` unchanged. Used by `__pow__(0)`.

### 1.10 `ZeroOperator` (`operator.py:719-736`)

`capabilities = frozenset({CAP_APPLY, CAP_APPLY_TRANSPOSE})` —
**not** `solve` (deliberate; line 723-727). Used by `KEigenvalue`'s
inner solve (`iteration.py:523`) to express the operator triple
`(L, S, 0)` for `(L-S) ψ = F·ψ/k`.

### 1.11 Other primitives in the module

- `PermutationOperator` (`operator.py:739`) — index permutation
  along a tagged axis, with `apply_transpose` = inverse permutation,
  `solve = apply_transpose` (for permutations).
- `IncomingOrdinateMaskTensor` (`operator.py:843`) — incoming-half
  mask, capabilities `{apply, apply_transpose}`.
- `PeriodicWrapOperator` (`operator.py:933`).
- `TensorProductOperator` (`operator.py:981`),
  `SumOfTensorProductsOperator` (`operator.py:1102`) — Grand
  Report §15.1 forms.
- `DiagonalOperator` (`operator.py:1195`) — `apply(x) = w * x`
  along a tagged axis; `solve` divides; full caps. Has
  `from_measure(measure, axis=0)` classmethod (line 1263).

### 1.12 `as_scipy_linop` adapter (`operator.py:1318-1368`)

```python
as_scipy_linop(op: LinearOperator, shape, dtype=float) → scipy LinearOperator
```

- Raises `MissingCapability` if `op` lacks `CAP_APPLY` (1355-1360).
- `matvec = op.apply` (line 1362).
- `rmatvec = op.apply_transpose` ONLY if op advertises
  `CAP_APPLY_TRANSPOSE` (`1363-1367`) — else `None`. The capability
  set is the single source of truth for whether `rmatvec` is exposed.

### 1.13 Tests for the algebra

- `tests/numerics/test_operator.py` (699 lines) — foundation-tagged
  invariants for every composer and every leaf. The protocol
  contracts are gated here.
- `tests/numerics/test_iteration.py` (462 lines) — see §3.4 below.
- `tests/numerics/test_diagonal_operator.py`,
  `test_tensor_product_operator.py`,
  `test_permutation_operator.py`,
  `test_periodic_wrap_operator.py`,
  `test_incoming_ordinate_mask_tensor.py`,
  `test_projection_operators.py`.

There is NO public `Representation` hierarchy in `operator.py`
(grand-report §32.8 forecast it; not landed). What exists is the
`(R, Π) = (Reconstruction, Projection)` pair in
`orpheus/numerics/projection.py` (§5 below).

---

## 2. Existing `ScatteringOperator` / `FissionOperator`

### 2.1 `S` — `ScatteringOperator` (`orpheus/sn/scattering.py`)

- Inherits `LinearOperatorMixin` (`scattering.py:254-255`).
- `capabilities = frozenset({CAP_APPLY})` (`scattering.py:317-319`).
  No `solve` (`scattering.py:21-29` — rank prevents inversion); no
  `apply_transpose` yet (not currently consumed; deferrable to a
  future wave).
- Constructed via factory `from_solver_data(...)`
  (`scattering.py:321-359`) — the canonical path the SN solver uses
  in `__init__` (`solver.py:240-252`).
- Data members (`scattering.py:303-319`): `n_ordinates`, `nx`, `ny`,
  `ng`, `scattering_order`, `sig_s` (per-material per-ℓ list of
  `(ng, ng)`), `sig2` (per-material n2n), `sig_s0` (P0 alias),
  `Y` (real spherical harmonics on quadrature), `weights`,
  `cells_by_mat`.
- **`apply(psi)` body** (`scattering.py:468-524`): the full operator
  surface — does P0 in-scatter + n2n + Pℓ Galerkin composition
  `R · Λ · M` on the angular flux. Returns per-ordinate
  `(N, nx, ny, ng)` source. **The "harmonic factoring `S = R · Λ · M`"
  the user asked about is ALREADY implemented inside
  `apply`** (lines 511-524).
- Convenience helpers used by the current SI sweep
  (`scattering.py:363-464`):
  - `add_iso_source(Q, phi)` — P0 in-place.
  - `add_n2n_source(Q, phi)` — n2n in-place.
  - `build_aniso_source(psi)` — Pℓ ≥ 1, returns the
    `(R Λ M ψ)` result as `(N, nx, ny, ng)` or `None`.

### 2.2 The `R · Λ · M` factoring already exists

In `scattering.py:447-464` (inside `build_aniso_source`):

```python
M   = HarmonicMomentProjection(weights=self.weights, Y=self.Y, L=L)
Lam = LegendreMomentScattering(sig_s=..., cells_by_mat=..., L=L, skip_l0=True)
R   = HarmonicMomentReconstruction.from_Y(self.Y)
return R.apply(Lam.apply(M.apply(angular_flux)))
```

All three (`R`, `Λ`, `M`) are `LinearOperator`s in their own right:

- `HarmonicMomentProjection` (M):
  `orpheus/numerics/projection.py` — a `GalerkinProjection`
  subclass. `capabilities = {apply, apply_transpose}` (line 127).
- `HarmonicMomentReconstruction` (R):
  `orpheus/numerics/projection.py`. Same surface.
- `LegendreMomentScattering` (Λ): in `scattering.py:132-251`,
  inherits `LinearOperatorMixin`, `capabilities = {CAP_APPLY}`
  (line 208-210). The per-ℓ block-diagonal `Σ_ℓ Pℓ ⊗ Σ_{s,ℓ}` per
  the reconciliation memo §2.5 / Grand Report §15.2.

So the "is any of that infrastructure already there?" question is
answered: **yes — both the operators and the composition; the
composition is just done as imperative `.apply(.apply(.apply(...)))`
calls rather than declared as `(R @ Lam @ M)`.** The composer is
on the bench, ready to be wired.

### 2.3 `F` — `FissionOperator` (`orpheus/sn/fission.py`)

- Inherits `LinearOperatorMixin` (`fission.py:82-83`).
- `capabilities = frozenset({CAP_APPLY})` (`fission.py:120-122`).
  No `solve` (rank-1 in energy, line 58-63); no `apply_transpose`
  yet (line 64).
- Factory `from_solver_data(chi, sig_p)` (`fission.py:124-135`).
- **`apply(phi)`** (`fission.py:137-165`) — returns `F·φ`, NOT
  `F·φ / k`. The `1/k` division is deliberately kept at the solver
  level (`fission.py:37-54`) — see Issue #169 / reconciliation
  memo. The user mandate `(L + C - S - F/k) ψ = q` is therefore
  implementable as `L + C - S - (1/k) * F` via the existing
  `ScaledOperator` composer; `F` itself is the linear operator
  `(F φ)_g = χ_g · Σ_{g'} νΣ_f^{g'} · φ_{g'}`.

### 2.4 Construction call sites

The SN solver builds `(S, F, L)` in its `__init__`
(`orpheus/sn/solver.py:235-262`):

```python
self.scattering_op = ScatteringOperator.from_solver_data(...)
self.fission_op    = FissionOperator.from_solver_data(chi=self.chi, sig_p=self.sig_p)
self.L = SNStreamingOperator(sn_mesh=sn_mesh, sig_t=self.sig_t)
self.S = self.scattering_op
self.F = self.fission_op
```

The four-operator algebra is **one step away**: the `(L+C)` split
is what `orpheus/sn/streaming.py` already did for the curvilinear
within-group sweep (see §4 below); the `SNSolver` itself still uses
the monolithic `L = SNStreamingOperator` (which carries Σ_t
internally).

### 2.5 The `S = R · Λ · M` Sphinx narrative

`docs/theory/operator_algebra.rst:362-363` already calls out
"`S_{SN} = R Λ M` becomes the operator-algebra…" and
`scattering.py:404-429` carries the full Galerkin derivation. The
algebra-of-record vocabulary is in place.

---

## 3. Existing iteration infrastructure

File: `orpheus/numerics/iteration.py` (569 lines).

### 3.1 `SourceIteration` (`iteration.py:164-352`)

**Signature** (`iteration.py:241-292`):

```python
SourceIteration(
    L: LinearOperator,        # MUST advertise CAP_APPLY
    S: LinearOperator,        # MUST advertise CAP_APPLY (use ZeroOperator if absent)
    F: LinearOperator,        # MUST advertise CAP_APPLY (use ZeroOperator if absent)
    *,
    inverter: Callable[[ndarray], ndarray] | None = None,
    max_iter: int = 1000,
    tol: float = 1e-8,
)
```

- **The Step 3 vision is already the shape.** `(L, S, F)` plus
  external source — exactly the user mandate, modulo `C` being
  rolled into `L` for now (the existing `L` is the monolith
  `SNStreamingOperator` carrying Σ_t).
- Capability checks happen at construction (`253-278`) — no stub
  failures mid-iteration.
- The `inverter` hook (`209-229`) is the load-bearing Wave-E design:
  when `None`, routes through `L.solve` (the WDD sweep — which is
  ERR-026-affected on curvilinear); when provided, the caller
  injects `lambda q: gmres(as_scipy_linop(L), q, M=...)` to use the
  symmetric closure with the sweep as a preconditioner. The same
  `SourceIteration` class drives both paths.
- `solve(q_ext, initial_guess=None)` (`294-352`):
  the iteration body is

  ```python
  rhs = self.F.apply(psi) + self.S.apply(psi) + q_ext
  psi = self._inverter(rhs)        # L^{-1} · rhs
  res = ||psi - psi_prev||_2 / max(||psi||, 1e-30)
  ```

  shape-agnostic; works equally on `(N,)` synthetic vectors and on
  `(nx, ny, ng)` SN scalar-flux arrays.
- Returns `(psi, residual_history)`.

### 3.2 `KEigenvalue` (`iteration.py:360-569`)

**Signature** (`iteration.py:428-471`):

```python
KEigenvalue(
    L, S, F,
    *,
    inverter=None,
    max_outer=500,
    keff_tol=1e-7, flux_tol=1e-6,
    max_inner=1000, inner_tol=1e-8,
    eigenvalue_method="power",                          # FEAST hook, only "power" today
    keff_estimator: Callable | None = None,             # default: _default_keff_estimator
)
```

- Capability validation done by trial-constructing the inner
  `SourceIteration` at `__init__` (line 477).
- `solve(initial_guess)` body (`479-569`):
  - inner solve: `inner = SourceIteration(L, S, ZeroOperator(), ...)`
    (line 523-531) — fission is the **external source** for the
    inner loop, not a within-group fixed-point term.
  - outer: `q_fission = F.apply(psi) / keff`, then
    `psi, _ = inner.solve(q_fission, initial_guess=psi_old)`.
  - keff update: `keff = self._keff_estimator(L, S, F, psi)`.
  - convergence requires both `dk < keff_tol` and
    `dphi < flux_tol`, with min 3 outer iterations
    (mirrors `SNSolver.converged` at `solver.py:360-370`).
- The `_default_keff_estimator` (`iteration.py:128-156`) is a
  generic Rayleigh quotient
  `k = sum(F·ψ) / (sum(L·ψ) - sum(S·ψ))`. The SN consumer should
  inject `sn_keff_estimator = lambda L,S,F,phi: solver.compute_keff(phi)`
  to preserve volume-weighted production/absorption math.

### 3.3 The deprecated `power_iteration` (`orpheus/numerics/eigenvalue.py`)

- 174 lines.
- `EigenvalueSolver` Protocol (`eigenvalue.py:44-117`) — the
  per-solver duck-typed contract (`initial_flux_distribution`,
  `compute_fission_source`, `solve_fixed_source`, `compute_keff`,
  `converged`).
- `power_iteration(solver, max_iter=500)` (`120-156`) — the legacy
  outer loop that `solve_sn` still calls (`solver.py:899`).
- **Fires `DeprecationWarning` at import** (`eigenvalue.py:166-173`)
  with the migration target named explicitly. The deprecation says
  it stays alive through the cross-solver migration sequence
  (CP, diffusion, MoC, homogeneous).

There is **no `PreconditionedGMRES`** wrapper class anywhere in
`orpheus.numerics`. The GMRES wiring is done inline at two sites:
`SNSolver._solve_krylov` (`solver.py:424-551`) and
`_solve_fixed_source_krylov` (`solver.py:1132-1289`). The natural
home for a `PreconditionedGMRES(L, q, M=...)` LinearOperator-aware
wrapper does not exist yet.

### 3.4 Tests (`tests/numerics/test_iteration.py`, 462 lines)

- `test_source_iteration_recovers_direct_solve` (L0, 4×4 dense
  matrices) — line 92.
- `test_source_iteration_with_fission_term` — line 131.
- `test_source_iteration_inverter_override` — line 158.
- `test_keigenvalue_recovers_dominant_eigenvalue` — line 193.
- Negative capability tests (`246-310`) — verify
  `MissingCapability` is raised at construction for every required
  capability.
- **`test_keigenvalue_matches_solve_sn_2g_slab`** (line 328) —
  L1 SN bridge: constructs `(SNStreamingOperator, ScatteringOperator,
  FissionOperator)` and runs `KEigenvalue` directly against the
  legacy `solve_sn` reference. **The test exists and passes
  TODAY** — proof the iteration primitives are ready to consume
  the SN operator triple. But it requires **adapter shims** on each
  of L/S/F to reconcile shape mismatches (line 402-440). Those
  shims are the "small bridge" that Step 3 of the replan needs to
  remove.

---

## 4. The current `_solve_fixed_source_si` and `_solve_fixed_source_krylov`

File: `orpheus/sn/solver.py` (1289 lines).

### 4.1 Dispatch tree

- `solve_sn(...)` (`solver.py:824-921`) — eigenvalue entry point.
  Always uses `inner_solver="source_iteration"` by default on ALL
  geometries (Cartesian + curvilinear). Calls legacy
  `power_iteration(solver)` (line 899). The inner solve is
  `SNSolver.solve_fixed_source` (line 281-293).
- `solve_sn_fixed_source(...)` (`solver.py:924-1064`) — fixed-source
  entry point. Auto-flips `inner_solver` based on curvature
  (line 1030-1034): curvilinear → `"krylov"`, Cartesian →
  `"source_iteration"`. Phase D rationale at `solver.py:1010-1029`.
- `SNSolver.solve_fixed_source(...)` (`solver.py:281-293`) routes
  internally to `_solve_source_iteration` or `_solve_krylov`.

### 4.2 `_solve_source_iteration` (`solver.py:374-420`, ~47 lines)

Approach: bit-identical preservation of the Wave-A-D inline loop.
Structure (conceptually a `SourceIteration` realisation with
`(L, S, ZeroOperator)` and `q_ext = fission_source`, but it does
NOT use the iteration primitive):

```python
for n_inner in range(max_inner):
    Q = fission_source.copy()
    solver._add_scattering_source(Q, phi)   # delegates to scattering_op.add_iso_source
    solver._add_n2n_source(Q, phi)          # delegates to scattering_op.add_n2n_source
    Q_aniso = solver._build_aniso_scattering(angular)  # delegates to scattering_op.build_aniso_source
    angular, phi = transport_sweep(Q, sig_t, sn_mesh, psi_bc, Q_aniso=Q_aniso)
    # convergence check on relative ||phi - phi_prev|| / ||phi||
```

This is the Step 3 target: should become

```python
SourceIteration(L, S, ZeroOperator()).solve(q_ext=fission_source)
```

once L/S adapter shims are resolved. Comment at
`solver.py:386-389` already flags this:

> "SourceIteration cannot be directly composed because the Pℓ
> anisotropic source requires per-ordinate angular-flux state
> across iterations — a future cleanup is to thread that state
> through ScatteringOperator so the primitive composes cleanly."

### 4.3 `_solve_krylov` (`solver.py:424-625`, ~200 lines)

GMRES on `self.L.apply` (the symmetric closure) with sweep
preconditioner. The packed-vector layout is built via
`EquationMap`. Highlights:

- `L_scipy = ScipyLinearOperator((n, n), matvec=self.L.apply,
  dtype=float)` (`solver.py:495-497`) — this is the manual
  equivalent of `as_scipy_linop(self.L, shape=(n,n))`. **The
  existing `as_scipy_linop` helper is not used here yet** —
  candidate for cleanup.
- Sweep wrapped as left preconditioner via
  `self._make_sweep_preconditioner` (`solver.py:553-625`). A
  stateless scipy `LinearOperator` that decodes the packed RHS
  → angular flux → runs `transport_sweep` → re-packs.
- Warm start, RHS construction (cartesian / spherical / cylindrical
  branches), GMRES dispatch with TypeError fallback for older
  scipy `tol` vs `rtol` (lines 513-526).

### 4.4 `_solve_fixed_source_si` (`solver.py:1067-1129`, ~63 lines)

Module-level helper called by `solve_sn_fixed_source` when
`inner_solver="source_iteration"`. Mirrors
`SNSolver._solve_source_iteration` — same `for n_inner in range:`
loop, but builds `Q` from `external_source` (`Q_aniso_total =
Q_aniso_p1 + external_source`) instead of `fission_source`. Bit
identical math to Wave A-D.

### 4.5 `_solve_fixed_source_krylov` (`solver.py:1132-1289`, ~158 lines)

Module-level helper. The pattern is the same as `_solve_krylov` but
runs an OUTER source iteration with a GMRES inner solve, because
scattering self-consistency requires re-building the RHS from the
current scalar flux each outer step. Reuses
`solver._make_sweep_preconditioner` (`solver.py:1189`). At
`solver.py:1228`: `rhs = rhs_iso + ext_packed_flat` — the external
per-ordinate source is added on top of the in-scatter RHS.

### 4.6 `SNStreamingOperator` — the existing monolithic `L` (`orpheus/sn/operator.py:1115-1472`)

This is the `L` that the current `(L, S, F)` triple uses.

- `capabilities = frozenset({CAP_APPLY, CAP_SOLVE, CAP_APPLY_TRANSPOSE})`
  (`operator.py:1244-1248`).
- **`L = Ω·∇ + Σ_t`** — carries Σ_t internally. Constructor takes
  `(sn_mesh, sig_t)` (lines 1241-1242). **This is what the four-
  operator algebra splits**: under the new mandate, `L` should be
  pure streaming (`Ω·∇` + curvilinear angular-redist) and `Σ_t`
  lives on a separate `C`.
- `apply(psi)` (`1287-1360`): dispatches by curvature to
  `transport_operator_matvec` (Cartesian),
  `transport_operator_matvec_spherical`,
  `transport_operator_matvec_cylindrical`. **The 850 lines of
  procedural matvec are the audit's smoking gun** — the prior
  audit memo (`.claude/agent-memory/explorer/issue_196_sn_operator_architecture_audit.md`)
  documents this in full.
- `solve(Q, psi_bc=None, Q_aniso=None)` (`1364-1407`): delegates to
  `transport_sweep`. Takes structured arrays, not the packed vector
  `apply` consumes — the shape-mismatch documented at
  `operator.py:1192-1212`. This is one of the "small bridges" the
  test shims expose.
- `apply_transpose(psi)` (`1442-1472`): assembled by probing
  `apply` with each unit basis vector once → `_dense_matrix`.
  O(n²) one-time cost; suitable for the small reciprocity test
  problems, NOT for production. Wave E flagged for an analytic
  matvec.

### 4.7 The Phase G Step 2 `(L, C)` split — already in tree

`orpheus/sn/streaming.py` (745 lines) — module docstring at
`streaming.py:1-90` explicitly calls itself "Issue #196 Phase G
Step 2".

**Status**: PARTIALLY landed. Maps the four-operator algebra's two
balance leaves:

- `StreamingOperator` (`streaming.py:251-408`) — **pure streaming
  (+ M-M angular redistribution for curvilinear)**.
  `capabilities = frozenset({CAP_APPLY})` (line 305). `apply(psi,
  *, sigma_t)` (line 310-393): currently still delegates to the
  monolithic `transport_operator_matvec_*` functions in `operator.py`
  (line 338-393). Step 2 comment: "_In Step 2 these are still the
  source of truth; Step 4 will promote them into the unified
  streaming module._" `solve_streaming` raises (line 395-408) —
  L alone is not invertible by the WDD sweep; only `(L+C)` is.
- `CollisionOperator` (`streaming.py:417-512`) — `C·ψ = Σ_t·ψ`,
  self-adjoint, fully invertible by per-element divide.
  `capabilities = {CAP_APPLY, CAP_SOLVE, CAP_APPLY_TRANSPOSE}`
  (line 446-448). Handles both packed and structured shapes
  (line 479-512).
- `DiscreteOrdinatesPhaseSpace` (`streaming.py:121-141`) — typed
  bundle for `(sn_mesh, quadrature, n_groups)`. Pattern 4 (illegal
  states unrepresentable).
- `solve_within_group_sn_curvilinear(L, C, q, ...)`
  (`streaming.py:520-688`) — the canonical curvilinear "sweep"
  via GMRES on `L.apply` (the `(L+C)` action — see comment at
  `streaming.py:631-636`: "In the Step 2 wiring, C is consumed
  inside the (L+C).solve fusion which handles the reshape … Step 4
  will refactor L.apply to pure streaming + angular-redist and
  route C separately").
- `build_sn_operators(V, materials, *, boundary=None)`
  (`streaming.py:697-`) — returns `(StreamingOperator,
  CollisionOperator)` ONLY. **Does NOT yet return `(L, C, S, F)`.**
  The Phase G "grand-report §29.2 canonical line" `L, C, S, F =
  build_sn_operators(V, materials)` is half built; S and F are
  still constructed independently by `SNSolver.__init__`
  (`solver.py:240-262`).

### 4.8 Existing wiring sites that consume Phase G Step 2

`orpheus/sn/sweep.py:266-292` invokes
`solve_within_group_sn_curvilinear` via `build_sn_operators` for
**curvilinear geometry only**. Cartesian still uses the classical
SI sweep. Tracking comment at `sweep.py:266-292`:

```python
L = StreamingOperator(V=V, boundary=sn_mesh.bc_right)
C = CollisionOperator(V=V, sigma_t=sig_t)
return solve_within_group_sn_curvilinear(L, C, ...)
```

This is the toehold the Phase G replan can pull on.

---

## 5. `compute_keff`, `compute_group_production_rate`, etc.

All in `orpheus/sn/solver.py:295-358`.

- `SNSolver.compute_group_production_rate(phi)` — `solver.py:295-331`.
  Returns `(ng,)` array. Components: `νΣ_f·φ` (volume-integrated
  via `mesh.volume_measure` — Issue 9.6 wiring) plus per-material
  n2n contribution (factor of 2 for two-neutron-out yield, loop on
  `_cells_by_mat`).
- `SNSolver.compute_group_absorption_rate(phi)` — `solver.py:333-345`.
  Returns `(ng,)` array. `Σ_a·φ` volume-integrated.
- `SNSolver.compute_keff(phi)` — `solver.py:347-358`. Composition:
  `k = compute_group_production_rate(phi).sum() /
       compute_group_absorption_rate(phi).sum()`. The
  named-intermediates per Issue #169.

`KEigenvalue` does NOT consume these directly — its default
`_default_keff_estimator` is the generic Rayleigh quotient
(`iteration.py:128-156`). The SN-aware bridge is supplied by the
caller via the `keff_estimator` kwarg (`iteration.py:441`). The
`test_keigenvalue_matches_solve_sn_2g_slab` test
(`test_iteration.py:444-445`) shows exactly this:

```python
def sn_keff_estimator(L_, S_, F_, phi):
    return solver.compute_keff(phi)

ke = KEigenvalue(..., keff_estimator=sn_keff_estimator)
```

CP / MoC / homogeneous also implement `compute_keff` on their
solver objects (Nexus shows `CPSolver.compute_keff`,
`MOCSolver.compute_keff`), so the named-intermediate refactor of
Issue #169 propagated cross-solver.

---

## 6. The Sphinx narrative on the algebra

### 6.1 `docs/theory/operator_algebra.rst` (796 lines, the Phase 0 stub)

Key sections and their `:label:` references the implementation
should pin `verifies(...)` decorators against:

- `.. _operator-algebra:` (line 1) — the top-level reference target.
- **Key Facts** (line 24-76) — the operator-algebra is grounded in
  Trefethen & Bau §3.2; full algebra dunders documented; the
  `MissingCapability`-at-composition-time invariant.
- `:label: operator-fixed-source` (line 38) — `(L - S - F)ψ = q`.
  vv-status: "documented" (line 35).
- `:label: operator-eigenvalue` (line 51) — `(L - S)ψ = (1/k)Fψ`.
  vv-status: "documented" (line 48).
- `:label: operator-apply` (line 91), `:label: operator-solve`
  (line 103), `:label: operator-apply-transpose` (line 115) — the
  three primitive actions.
- **Capability set semantics** (lines 132-160) — explicit
  rationale that S has no `solve` (rank), F has no `solve`
  (rank-1 in energy), Jacobi preconditioner has `apply` but no
  `solve`. The harmful-stub anti-pattern is named.
- **Composition algebra** (lines 163-216) — full closure-laws
  table for `OperatorSum`, `OperatorProduct`, `ScaledOperator`,
  `IdentityOperator`, `ZeroOperator`. The exact rules from §1.6–
  1.10 above are documented authoritatively here.
- **Cross-solver consumption (forward reference)** (lines 219-240)
  — explicitly cites Wave-H Issue 14 / Issue 15 as the
  cross-solver migration target. "Today the existing
  `build_transport_linear_operator` functions in
  `orpheus.sn.operator` continue to wrap their matvec in
  `scipy.sparse.linalg.LinearOperator` directly, untouched by
  this module" (line 232-235) — historical context for why
  `as_scipy_linop` is not yet consumed inside the solver.
- **Diagonal operator on a tagged axis** (lines 244-319;
  `:label: diagonal-operator-action`).
- **Tensor product algebra** (lines 320-562;
  `:label: tensor-product-action`,
  `:label: streaming-as-tensor-product-sum`,
  `:label: scattering-as-tensor-product-sum`,
  `:label: octant-direct-sum-tensor-product`).
- **Trace decompositions** (lines 660-796;
  `:label: trace-half-decomposition`,
  `:label: per-face-inflow-mask`,
  `:label: bc-rank-n-as-sum-of-products`).

The Sphinx page is **the algebra-of-record**. Every Phase G
implementation step should pin its tests to one of these labels
via `@pytest.mark.verifies("…")`.

### 6.2 `docs/theory/discrete_ordinates.rst` — SN-side algebra references

Relevant labels for Phase G:

- `:label: transport-cartesian` (185),
  `:label: transport-cartesian-2d` (202),
  `:label: transport-spherical` (220),
  `:label: transport-cylindrical` (238) — the continuous balance.
- `:label: multigroup` (261).
- `:label: hebert-3-432`, `:label: hebert-3-432-source`,
  `:label: hebert-3-434`, `:label: hebert-3-435` (3375-3452) —
  Hébert §3.4.3 numbering: the canonical algebra-of-record for
  the SN discrete balance.
- `:label: dd-cartesian-1d`, `:label: dd-cartesian-2d`,
  `:label: wdd-closure`, `:label: wdd-face`, `:label: mm-weights`,
  `:label: dd-solve` (lines 432-818) — closures for the WDD/MM
  paths.
- The branch's existing narrative around the operator-algebra
  capabilities is at line 3844 ("operator-algebra capabilities
  of ..."). Phase F closeout (line 3326, line 4197) already
  references this branch as `refactor/sn-operator-algebra`.

`docs/theory/operator_algebra.rst:362` directly identifies
`S_{SN} = R Λ M` and acknowledges
`scattering-as-tensor-product-sum` as the formal algebraic form.

---

## 7. Anti-recommendations — DO NOT REINVENT

Pin these to memory before writing the next sub-agent brief:

1. **Do NOT reimplement `OperatorSum`, `OperatorProduct`, or
   `ScaledOperator`.** Full algebra dunders ship at
   `orpheus/numerics/operator.py:489-694`. Composition propagates
   `.H` correctly through `_AdjointOperator`
   (`operator.py:415-486`). `(L + C - S - F/k).H` works today as
   `L.H + C.H - S.H - (1/k) * F.H` — the only thing missing is
   that S and F haven't advertised `CAP_APPLY_TRANSPOSE` yet.

2. **Do NOT add a new `inner_solver` switch.** The existing
   string switch at `orpheus/sn/solver.py:162-170` already gates
   `"source_iteration"` vs `"krylov"`. The auto-flip is at
   `solver.py:1030-1034`. Step 3 of the replan should consolidate
   onto `SourceIteration`-with-inverter, not add a third path.

3. **Do NOT write a new `SourceIteration` or `KEigenvalue`.**
   `orpheus/numerics/iteration.py:164-569` already takes
   `(L, S, F, q_ext)` and `(L, S, F)` respectively — exactly the
   Step-3 target shape. The test
   `test_keigenvalue_matches_solve_sn_2g_slab`
   (`tests/numerics/test_iteration.py:328`) is the L1 gate that
   PROVES the SN triple works with these primitives today.

4. **Do NOT reinvent `S = R · Λ · M`.** The composition exists
   at `orpheus/sn/scattering.py:447-464` (inside
   `build_aniso_source`). `HarmonicMomentProjection` (M) and
   `HarmonicMomentReconstruction` (R) are first-class
   `LinearOperator`s in `orpheus/numerics/projection.py`;
   `LegendreMomentScattering` (Λ) is in
   `orpheus/sn/scattering.py:132-251`. The harmonic factoring is
   on the bench — what's missing is exposing `S` itself as
   `R @ Lam @ M` declaratively (currently `.apply().apply().apply()`
   imperatively).

5. **Do NOT reinvent the `(L, C)` split.**
   `orpheus/sn/streaming.py:251-512` already defines
   `StreamingOperator` and `CollisionOperator` separately. Step 4
   of the streaming-module docstring (`streaming.py:36-41`)
   explicitly identifies the remaining work: "_Step 4 will refactor
   L.apply to pure streaming + angular-redist and route C
   separately_". The promotion is the next concrete step, not a
   from-scratch rebuild.

6. **Do NOT reinvent `as_scipy_linop`.** It's at
   `operator.py:1318-1368`. Use it inside `_solve_krylov` and
   `_solve_fixed_source_krylov` instead of the bare
   `ScipyLinearOperator((n, n), matvec=self.L.apply, ...)` lines
   at `solver.py:495-497` and `solver.py:1186-1188`. That swap is
   ONE LINE per call site — a free win.

7. **Do NOT reinvent the keff-estimator hook.** It's already a
   first-class kwarg on `KEigenvalue.__init__`
   (`iteration.py:441`). The SN-aware estimator is a 2-line
   wrapper around `SNSolver.compute_keff` per
   `test_iteration.py:444-445`. The `(production rate, absorption
   rate)` named intermediates are in `solver.py:295-358` (Issue
   #169 landed). Nothing new needs writing.

8. **Do NOT write a new `PreconditionedGMRES` LinearOperator
   wrapper before checking inline use.** The GMRES wiring at
   `solver.py:495-526` and `solver.py:1186-1242` is the only
   consumer. Wrap-as-operator only if a SECOND consumer needs it
   (per `.claude/agent-memory/.../feedback_unify_after_two_instances.md`
   — build standalone first, unify after ≥2 working instances).

9. **Do NOT add `apply_transpose` stubs to `S` and `F` that raise
   `NotImplementedError`.** That's the harmful-stub anti-pattern
   (`operator_algebra.rst:149-156`). When adjoint sensitivity
   becomes a real consumer, advertise `CAP_APPLY_TRANSPOSE` and
   implement the action. Until then, `S.H` and `F.H` *should*
   raise `MissingCapability` at composition time — that's the
   correct behaviour, not a gap.

10. **Do NOT split `_solve_source_iteration` and
    `_solve_fixed_source_si` into different code paths.** They
    are the SAME loop — within-group fixed-point on
    `(L - S - F) ψ = q_ext` (with `F = 0` for fixed-source). A
    single `SourceIteration(L, S, F, ...).solve(q_ext)` call
    site is sufficient. Same for `_solve_krylov` and
    `_solve_fixed_source_krylov` — both are GMRES on the same
    operator action with a slightly different RHS assembly.

11. **Do NOT use `power_iteration` for new code.** It fires a
    `DeprecationWarning` on import
    (`orpheus/numerics/eigenvalue.py:166-173`). The migration
    target is `KEigenvalue`; `power_iteration` stays alive only
    through the cross-solver migration sequence (CP, diffusion,
    MoC, homogeneous).

---

## 8. "Almost there" pieces — small bridges away from the algebra

1. **`SNStreamingOperator` carries Σ_t internally.** Splitting into
   `L = SNStreamingOperator(sn_mesh)` (pure Ω·∇) plus
   `C = CollisionOperator(V, sigma_t)` is half-done in
   `orpheus/sn/streaming.py`. The remaining work (per the
   Step 4 forward reference at `streaming.py:36-41`) is to make
   `StreamingOperator.apply` pure-streaming + angular-redist,
   and route Σ_t exclusively through `C`. The Sphinx narrative
   (`operator_algebra.rst`) already calls the four-operator
   form `(L + C - S - F/k)` even though the implementation has
   `L = streaming + collision`.

2. **`S.apply` and `F.apply` shape mismatch with `L.solve`.** The
   adapter shims in `test_keigenvalue_matches_solve_sn_2g_slab`
   (`test_iteration.py:402-440`) wrap each operator into a
   scalar-flux-only `(nx, ny, ng)` facade. The comment at
   `solver.py:386-389` and `operator.py:1192-1212` flags this.
   The fix is to thread the angular-flux state through
   `ScatteringOperator` so the Pℓ source can be a pure linear
   operator on `psi` (not need-an-extra-cached-angular-flux);
   `S = R @ Lam @ M` then composes cleanly without shims.

3. **`build_sn_operators` returns only `(L, C)`.** Extending it
   to return `(L, C, S, F)` per Grand Report §29.2 is a small
   refactor. Materials → S, F is already in
   `SNSolver.__init__:240-255`.

4. **Two GMRES call sites duplicate scipy wiring.**
   `_solve_krylov` (`solver.py:495-497`) and
   `_solve_fixed_source_krylov` (`solver.py:1186-1188`) both
   manually wrap `self.L.apply` in
   `ScipyLinearOperator((n, n), matvec=...)`. Should use
   `as_scipy_linop(self.L, shape=(n, n))`.

5. **`compute_keff` not wired into `KEigenvalue` yet.**
   `SNSolver.compute_keff` exists; `KEigenvalue.keff_estimator`
   exists; the bridge is a 2-line `keff_estimator=lambda
   L_,S_,F_,phi: solver.compute_keff(phi)` at the consumer site.
   The L1 test already does this (`test_iteration.py:444-445`).

6. **`_AdjointOperator.apply_transpose` raises
   `NotImplementedError`** (`operator.py:481-486`). Not a Phase G
   blocker — the four-operator algebra `(L+C-S-F/k).H` only needs
   `.H` of leaves, not `.H.H`. Document this constraint in the
   replan if any adjoint-of-adjoint consumer is anticipated.

---

## 9. Source files referenced (absolute paths)

- `/Users/rodrigo/git/nuclear/ORPHEUS/orpheus/numerics/operator.py`
  (1368 lines, the algebra core)
- `/Users/rodrigo/git/nuclear/ORPHEUS/orpheus/numerics/iteration.py`
  (569 lines, SourceIteration + KEigenvalue)
- `/Users/rodrigo/git/nuclear/ORPHEUS/orpheus/numerics/eigenvalue.py`
  (174 lines, deprecated power_iteration)
- `/Users/rodrigo/git/nuclear/ORPHEUS/orpheus/numerics/projection.py`
  (477 lines, R/Π hierarchy + HarmonicMomentProjection +
  HarmonicMomentReconstruction)
- `/Users/rodrigo/git/nuclear/ORPHEUS/orpheus/sn/scattering.py`
  (524 lines, ScatteringOperator + LegendreMomentScattering)
- `/Users/rodrigo/git/nuclear/ORPHEUS/orpheus/sn/fission.py`
  (165 lines, FissionOperator)
- `/Users/rodrigo/git/nuclear/ORPHEUS/orpheus/sn/operator.py`
  (1473 lines — includes the procedural matvec sweep math at
  lines 412-1107 and SNStreamingOperator at 1115-1472)
- `/Users/rodrigo/git/nuclear/ORPHEUS/orpheus/sn/streaming.py`
  (745 lines, Phase G Step 2 partial: StreamingOperator,
  CollisionOperator, solve_within_group_sn_curvilinear,
  build_sn_operators)
- `/Users/rodrigo/git/nuclear/ORPHEUS/orpheus/sn/solver.py`
  (1289 lines — SNSolver, _solve_source_iteration, _solve_krylov,
  solve_sn, solve_sn_fixed_source, _solve_fixed_source_si,
  _solve_fixed_source_krylov)
- `/Users/rodrigo/git/nuclear/ORPHEUS/docs/theory/operator_algebra.rst`
  (796 lines, the algebra-of-record narrative)
- `/Users/rodrigo/git/nuclear/ORPHEUS/docs/theory/discrete_ordinates.rst`
  (5100+ lines, SN labels listed in §6.2)
- `/Users/rodrigo/git/nuclear/ORPHEUS/tests/numerics/test_iteration.py`
  (462 lines)
- `/Users/rodrigo/git/nuclear/ORPHEUS/tests/numerics/test_operator.py`
  (699 lines)

---

## 10. Pointer to the existing architecture audit

Prior memo (read before drafting the Phase G replan):
`.claude/agent-memory/explorer/issue_196_sn_operator_architecture_audit.md`

That memo diagnoses Manifestation #7 (SI-vs-Krylov O(h) drift) as
two procedural implementations of the same five operator-algebra
primitives. This inventory complements it by mapping what the
algebra-of-record SHIPS today, so the Phase G replan plugs into
existing infrastructure rather than reinventing it.
