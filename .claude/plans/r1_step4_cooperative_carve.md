# R-1 Step 4 — Cooperative carve plan (post-revert)

**Branch**: `refactor/sn-operator-algebra`
**Status (2026-05-19, end of cooperative carve session 1)**: Steps A–F landed (eigenvalue path fully on typed AngularFlux + InvertibleOperator algebra); Steps G + H deferred to a follow-up session because the legacy `solve_sn_fixed_source` MMS path also reads through `SNStreamingOperator` / `EquationMap` / `_build_rhs_*` and requires its own carve before those symbols can retire.

## Session 1 closeout

| Step | Status | Commit |
|---|---|---|
| A — AngularFlux history container | landed | `fa5d762` |
| B — drop SI `inverter`, rename Krylov's to `preconditioner` | landed | `4ad9cd3` |
| C — `InvertibleOperator` + `__add__` dispatch | landed | `2c634ab` |
| D — carve `_solve_krylov` | landed | `8ff0881` |
| E — carve `_solve_source_iteration` | landed | `22c4786` |
| F — `test_collision_operator.py` typed migration + `CollisionOperator.solve/apply_transpose` typed overloads | landed | `e7d715e` |
| G — retire SNStreamingOperator + EquationMap + legacy decoders | **DEFERRED** | — |
| H — Sphinx narrative + closeout | **DEFERRED** | — |

**L1 anchor status (after session 1):** `tests/sn/l1_analytical/`: 35 passed + 5 xfailed. The 5 xfails are all `sphere-4eg-krylov` cases — unpreconditioned GMRES exceeds the 300-iter budget on 4-group curvilinear; pending [#200](https://github.com/deOliveira-R/ORPHEUS/issues/200)'s block-inverse face preconditioner.

**Step F follow-on test sites** (the plan listed ~70 across 10 files but the actual failing-pre-Step-F count was 39, all in `test_collision_operator.py`). Other listed files (`test_streaming_operator.py`, `test_phase_c_gates.py`, `test_b1pp_verification.py`, `test_unified_matvec_*.py`) test the LEGACY `SNStreamingOperator` directly and stay until Step G retires that class.

## Step G deferral rationale — scope analysis

The plan listed Step G as "~1100 LoC retirement". Closer inspection during session 1 showed the production-code coupling is bigger than a pure retirement:

* `SNStreamingOperator` is read by `solve_sn_fixed_source` (public API for MMS) via `solver.L.apply` and `solver._make_sweep_preconditioner`.
* `EquationMap` + `build_equation_map*` + `solution_to_angular_flux*` are read by `_solve_fixed_source_krylov` and `_solve_fixed_source_si` (both module-level helpers behind `solve_sn_fixed_source`).
* `_build_rhs_cartesian` / `_spherical` / `_cylindrical` are also read by the legacy fixed-source path.

Step G as written would orphan the MMS path. The principled completion requires Step G to ALSO carve `_solve_fixed_source_krylov` / `_solve_fixed_source_si` onto `KrylovAcceleration` / `SourceIteration` + typed `AngularFlux` — analogous to Steps D + E for the eigenvalue path. After that, the legacy symbols can retire together.

## Session 2 plan

A future session should:

1. **G1 — carve `_solve_fixed_source_krylov`** onto `KrylovAcceleration` consuming `AngularFlux`, with the same `preconditioner=lambda q: q` identity-precond + `S/sum_w` rescaling we ironed out for the eigenvalue path.  L1 anchor: `tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py` (MMS path).
2. **G2 — carve `_solve_fixed_source_si`** onto `SourceIteration` consuming `AngularFlux`, with the `*sum_w` bridge inside `InvertibleOperator.solve` carrying through (already landed in Step E).
3. **G3 — retire**: once nothing reads them,
   - `SNStreamingOperator` (in `orpheus/sn/operator.py`)
   - `EquationMap` + `build_equation_map` + `build_equation_map_spherical` + `build_equation_map_cylindrical`
   - `solution_to_angular_flux` + `_spherical` + `_cylindrical`
   - `_make_sweep_preconditioner` (SNSolver method)
   - `_build_rhs_cartesian` + `_spherical` + `_cylindrical` (free functions)
   - Tests: `tests/sn/test_snstreamingoperator.py` (delete; superseded by typed-operator tests)
4. **G4 — relocate the keepers**: `pack_with_traces`, `solution_to_angular_flux_with_traces`, `build_equation_map_with_traces` → private helpers of `AngularFlux.from_flat_with_traces` / `to_flat_with_traces`.
5. **H — Sphinx narrative + closeout** as the plan originally describes.

The deferred work continues from session-1 HEAD `e7d715e`.

---

**Original plan continues below (preserved for completeness — pre-session-1 plan as authored on 2026-05-19).**

---

**Approach**: cooperative carve in conversation with user; no `method-implementer` agent dispatch in this session

---

## Context (after compaction, read this first)

### What's already done in R-1 (cumulative through HEAD `9410659`)

| Commit | What it landed |
|---|---|
| `fb64049` | `KrylovAcceleration` sibling primitive in `orpheus/numerics/iteration.py` |
| `b667f8e` | R-1 Step 0 — pole iteration history evicted to `initial_guess: AngularFlux` parameter on `transport_sweep` (trace-space sweep) |
| `5cce245` | R-1 Step 1 — `BoundaryFlux` arithmetic dunders (immutable arithmetic) |
| `43bb8e0` | R-1 Step 2 — `AngularFlux` embeds `BoundaryFlux` via `boundary` field with auto-allocation; `to_flat_with_traces` / `from_flat_with_traces` scipy boundary methods |
| `7ce5d53` | Test fixes: DD regression curvilinear tolerance + `kinf_homogeneous_spectrum` axis bug |
| `c6ba262` | C1 inner-tol amplification: tight tolerances in `test_kinf_homogeneous` + L1 companion test |
| `eeac45f` | R-1 Step 3 — four-operator typed `apply(AngularFlux) → AngularFlux` overloads (L, C, S, F). L is the only operator emitting non-zero face residuals. |
| `2499859` | R-1 Step 4a — `KrylovAcceleration` + `SourceIteration` accept `AngularFlux` via duck-typed `to_flat_with_traces` protocol; `ZeroOperator.apply` uses `0.0 * x` to preserve typed flux |
| `578fb46` | Safe retirement: `EquationMap.unknowns_at_cell_for_mask` + `_cell_ord_to_k` + 5 dependent tests (zero production callers) |
| `9410659` | Revert of `86b1da5` (bad agent attempt) — agent ignored S split, added +279 LoC of accretion, used identity-on-face hack labeled as "linearity bug fix" |

### Architectural decisions made with the user

1. **The packed-1D path retires.** Iteration primitives consume typed AngularFlux end-to-end via the ravellable protocol (`to_flat_with_traces` at scipy boundary only).
2. **The S split (foldable / residual)** — postponed as a convergence-rate optimization. R-1 lands the basic L+C inversion; the σ_r foldable preconditioner is GitHub issue [#200](https://github.com/deOliveira-R/ORPHEUS/issues/200).
3. **Krylov has NO preconditioner in R-1.** The block-inverse face-block preconditioner is also issue [#200](https://github.com/deOliveira-R/ORPHEUS/issues/200). User direction: "I'm more interested in consolidating the foundational architecture for now."
4. **2-D Cartesian is deferred.** R-1 carves 1-D only (slab + sphere + cylinder). 2-D Cartesian dispatch raises `NotImplementedError` with a clear deferral message. Phase A redoes 2-D properly.
5. **AngularFlux gets a list-like history container** (Form B2). User confirmed this scaffolds future time-derivative work.
6. **Dimensional typing follow-up** filed as issue [#201](https://github.com/deOliveira-R/ORPHEUS/issues/201) — `AngularFlux` vs `AngularSource` vs `AngularResidual` as distinct types. Out of scope for R-1.
7. **`inverter` parameter retires.** SourceIteration uses `L.solve` directly. The redundant override goes away. KrylovAcceleration renames `inverter` → `preconditioner` (different concept).
8. **`InvertibleOperator`** (specialization of `OperatorSum`) is the type that carries the "this sum is sweep-invertible" algebraic identity. Both `L + C_t` AND `L + C_r` (where `C_r` is built with `σ_r = σ_t - Σ_s0(g→g)`) construct InvertibleOperator with `.solve = sweep`.

### Critical preserved findings from the bad agent

Despite the bad commit being reverted, **two findings from `scratch/r1_step4_*.py` diagnostics are real and load-bearing**:

* **Face-block linearity** (issue [#200](https://github.com/deOliveira-R/ORPHEUS/issues/200)) — encoding the sweep's outflow face value into the GMRES preconditioner's face block produces wrong answers at machine-zero residual because `M·A` becomes rank-deficient. Block-diagonal "identity on face" is the deliberate clean choice when preconditioning is enabled. R-1 sidesteps via "no preconditioner".
* **SI Carlson seed dependency** — SourceIteration's iteration step needs the previous-iterate angular flux as the Carlson seed for the curvilinear sweep. This is what motivates the AngularFlux history container.

---

## Plan — full implementation in 8 steps

### Step A — AngularFlux as list-like history container (Form B2)

**File**: `orpheus/sn/angular_flux.py`

**Approach**: keep the public surface (`values`, `boundary`, `mesh`) intact; add an internal `_history` tuple and `history_depth` constructor argument; add `__call__`, `stash`, `__lshift__`; arithmetic drops history.

```python
@dataclass(frozen=True)
class AngularFlux:
    values: np.ndarray
    mesh: "SNMesh"
    boundary: "BoundaryFlux | None" = None
    history_depth: int = 2            # max frames including current
    _history: tuple[tuple[np.ndarray, "BoundaryFlux"], ...] = ()

    def __post_init__(self):
        # Existing shape validation
        # Existing boundary auto-allocation
        # NEW: validate history_depth >= 1

    def __call__(self, lag: int = 0) -> "AngularFlux | None":
        """Return the frame at the given lag (0 = current, 1 = previous, ...).
        
        Returns a fresh single-frame AngularFlux (no history attached).
        Returns None if the lag is beyond the available history.
        """
        if lag == 0:
            return self
        idx = lag - 1
        if idx >= len(self._history):
            return None
        v, b = self._history[idx]
        return AngularFlux(v, self.mesh, b, history_depth=1)

    def stash(self, new: "AngularFlux") -> "AngularFlux":
        """Shift register: new becomes head; self becomes previous; oldest dropped.
        
        Self's history_depth governs the resulting depth.  Returns a fresh
        AngularFlux; self is unchanged (frozen dataclass).
        """
        self._validate_partner(new)
        # New history slot = (self.values, self.boundary) prepended;
        # truncate at history_depth - 1 (one less because the head IS the new one).
        new_history = ((self.values, self.boundary),) + self._history[: self.history_depth - 2]
        return AngularFlux(
            new.values, self.mesh, new.boundary,
            history_depth=self.history_depth,
            _history=new_history,
        )

    def __lshift__(self, new: "AngularFlux") -> "AngularFlux":
        return self.stash(new)

    def __len__(self) -> int:
        return 1 + len(self._history)   # current + history

    # Arithmetic dunders — operate on HEAD only; result has no history.
    def __add__(self, other):
        return AngularFlux(
            self.values + other.values,
            self.mesh,
            self.boundary + other.boundary,
            history_depth=self.history_depth,
            _history=(),
        )
    # ... similar __sub__, __mul__, __rmul__, __truediv__, __neg__
```

**Tests**: `tests/sn/test_angular_flux_with_boundary.py` extended with:
- `test_history_depth_default_is_2`
- `test_call_zero_returns_self`
- `test_call_lag_returns_previous`
- `test_call_beyond_history_returns_none`
- `test_stash_shift_register_semantics`
- `test_stash_truncates_at_history_depth`
- `test_arithmetic_drops_history`
- `test_lshift_alias_for_stash`
- `test_len_returns_total_frames`
- `test_kinetics_depth_5_holds_5_levels`

**Gotcha**: `_history` must be a tuple-of-tuples (immutable) for the frozen dataclass. The `(values, boundary)` pair gets stored verbatim — `boundary` shares state with the original BoundaryFlux (which is mutable). For frozen-arithmetic semantics, callers must not mutate stashed boundaries. Document this.

**Acceptance**: all existing AngularFlux tests pass unchanged (the new fields default to depth=2 and empty `_history` — pre-existing call sites are unaffected). New tests pass.

**Commit message**: `feat(sn): AngularFlux history container — list-like shift register for iterate history (R-1 Step 4 A)`

---

### Step B — Drop `inverter` from SourceIteration; rename Krylov's

**File**: `orpheus/numerics/iteration.py`

**Why**: per user's question — "Why do we have to create an .inverter when .solve is the action?" Answer: legacy redundancy. The correct architecture is `SourceIteration` uses `L.solve` directly (requires CAP_SOLVE on L); the caller controls behavior by constructing L appropriately. KrylovAcceleration's `preconditioner` parameter is a different concept (approximation to A⁻¹, optional).

**Changes**:

```python
# orpheus/numerics/iteration.py

class SourceIteration:
    def __init__(self, L, S, F):
        # OLD: also took inverter parameter
        if not _has(L, CAP_APPLY):
            raise MissingCapability(...)
        if not _has(S, CAP_APPLY):
            raise MissingCapability(...)
        if not _has(F, CAP_APPLY):
            raise MissingCapability(...)
        if not _has(L, CAP_SOLVE):
            raise MissingCapability(
                f"SourceIteration requires {CAP_SOLVE!r} on L; the iteration "
                f"step is ψ_(n+1) = L.solve(F·ψ_n + S·ψ_n + q_ext).  Pass an "
                f"L operator whose .solve method realises L⁻¹ for your problem."
            )
        self.L, self.S, self.F = L, S, F
        # ... tol, max_iter

    def solve(self, q_ext, initial_guess=None):
        if initial_guess is None:
            psi = _zeros_like(q_ext)
        else:
            psi = initial_guess
        residual_history = []
        for _ in range(self.max_iter):
            psi_prev = psi
            # Compose iteration RHS via operator algebra.
            rhs_raw = self.F.apply(psi) + self.S.apply(psi) + q_ext
            # R-1 Step 4 — attach psi as the .previous of rhs so L.solve
            # can read it for the Carlson seed (curvilinear sweeps need it).
            # The arithmetic above dropped history; we explicitly thread psi
            # into rhs as its previous frame so the sweep adapter sees it.
            rhs = _attach_previous(rhs_raw, psi)  # see helper below
            psi = self.L.solve(rhs)  # ← direct .solve call, no inverter
            # Convergence check ...

class KrylovAcceleration:
    def __init__(self, L, S, F, *, preconditioner: Callable | None = None, ...):
        # OLD parameter name: ``inverter``
        # NEW: ``preconditioner`` — clearer semantic (approximation to A⁻¹).
        # User direction for R-1: pass None (no preconditioner).
        # ...
        if preconditioner is not None:
            self._preconditioner = preconditioner
        elif _has(L, CAP_SOLVE):
            # By default, L.solve serves as preconditioner if L has CAP_SOLVE.
            # User can opt out by passing preconditioner=None explicitly.
            self._preconditioner = lambda q: self.L.solve(q)
        else:
            self._preconditioner = None
```

**Helper** (private to iteration.py):

```python
def _attach_previous(rhs, previous):
    """Attach ``previous`` as the lag-1 frame of ``rhs`` (typed-flux only).
    
    For bare ndarray rhs, returns rhs unchanged (no protocol; SI on bare
    arrays doesn't need a Carlson seed).
    """
    if not _is_ravellable(rhs):
        return rhs
    # rhs is an AngularFlux-like; create a new instance with previous as 
    # the lag-1 frame.  Use the same history_depth as rhs.
    return type(rhs)(
        rhs.values, rhs.mesh, rhs.boundary,
        history_depth=rhs.history_depth,
        _history=((previous.values, previous.boundary),),
    )
```

**Tests** to update: `tests/numerics/test_iteration.py` — 8 sites that pass `inverter=...`. Migrate by:
- Wrapping the inverter in a tiny operator class: `class _TestInverseOp: capabilities = {CAP_APPLY, CAP_SOLVE}; apply = lambda self, x: ...; solve = lambda self, q: inverter(q)`
- Pass this wrapped op as L
- Result: 0 `inverter=` references remaining in the codebase

**Tests for KrylovAcceleration**: rename `inverter` → `preconditioner` everywhere in `test_iteration.py` + `test_iteration_angular_flux.py`. Add a test that `KrylovAcceleration(L, S, F, preconditioner=None)` works (no preconditioning).

**Acceptance**: all 27 existing iteration tests pass with the renamed parameter; no `inverter` references remain.

**Commit message**: `refactor(numerics): retire inverter param — SI uses L.solve directly; Krylov renames to preconditioner (R-1 Step 4 B)`

---

### Step C — `InvertibleOperator` class + `StreamingOperator.__add__(CollisionOperator)` dispatch

**File**: `orpheus/sn/operator.py`

**Why**: `(L + C).solve` is the SN-specific algebraic identity (sweep IS the inverse algorithm). The class encodes this identity at the type level. Constructed via `L + C` automatically.

**The class**:

```python
class InvertibleOperator(OperatorSum):
    r"""Streaming + diagonal collision — sweep-invertible.
    
    Specialization of :class:`~orpheus.numerics.operator.OperatorSum` for
    the SN-specific algebraic identity:
    
    .. math::
    
        (L_{\rm streaming} + C_{\rm diagonal})^{-1} \approx \text{WDD sweep}
    
    The discrete-ordinates WDD sweep IS the inverse algorithm for this
    specific composite — that's the algebraic foundation of the entire SN
    method (Lewis & Miller §3, Adams & Larsen 2002).  No generic
    ``(A+B)^{-1}`` formula exists; this is a domain-specific identity.
    
    Constructed two ways:
    
    * ``L + C`` — when L is a :class:`StreamingOperator` and C is a
      :class:`CollisionOperator`, ``StreamingOperator.__add__`` dispatches
      to this class automatically.  The composite reads as math.
    * ``InvertibleOperator(L, C)`` — explicit construction; useful when
      composing variants such as ``InvertibleOperator(L, C_r)`` where
      ``C_r`` carries the foldable-scattering-modified diagonal
      :math:`\sigma_r = \sigma_t - \Sigma_{s,0}^{g \to g}` (Adams & Larsen
      2002 §III; tracked by issue #200).
    
    The class advertises both ``apply`` (inherited from OperatorSum,
    sums the operand actions) and ``solve`` (new — the sweep adapter).
    
    Parameters
    ----------
    streaming : StreamingOperator
        The streaming operator :math:`L = \Omega \cdot \nabla`.  Resolution A
        subtractive form: ``L.apply(\psi) = M(\psi; \sigma_t) - \sigma_t \odot \psi_{\rm cell}``.
    diagonal : CollisionOperator
        The diagonal collision operator :math:`C = \sigma \odot`.  Its
        ``.sigma`` attribute is the per-cell per-group coefficient used by
        the sweep (canonically :math:`\sigma_t`; can be :math:`\sigma_r`
        for the foldable variant).
    """
    
    capabilities = frozenset({CAP_APPLY, CAP_SOLVE})
    
    def __init__(self, streaming: "StreamingOperator", diagonal: "CollisionOperator"):
        # Validate mesh consistency
        if streaming.sn_mesh is not diagonal.sn_mesh:
            raise ValueError(
                "InvertibleOperator: streaming and diagonal operators must "
                "share the same SNMesh instance."
            )
        # Validate that the diagonal coefficient is positive-definite (a 
        # necessary condition for the sweep to be well-defined).
        if not np.all(diagonal.sigma > 0):
            raise ValueError(
                "InvertibleOperator: diagonal coefficient must be strictly "
                "positive everywhere for the sweep to be well-defined.  "
                "Got non-positive values in CollisionOperator.sigma."
            )
        super().__init__(streaming, diagonal)
    
    @property
    def streaming(self) -> "StreamingOperator":
        return self.a   # OperatorSum stores operands as .a and .b
    
    @property
    def diagonal(self) -> "CollisionOperator":
        return self.b
    
    @property
    def sn_mesh(self) -> "SNMesh":
        return self.streaming.sn_mesh
    
    @property
    def sigma(self) -> np.ndarray:
        """The diagonal coefficient (σ_t or σ_r) used by ``solve``."""
        return self.diagonal.sigma
    
    def solve(self, rhs: "AngularFlux") -> "AngularFlux":
        r"""Apply the SN sweep using ``self.sigma`` as the diagonal coefficient.
        
        The sweep inverts ``(L + C)`` by Wave-Diamond-Difference cell
        sweep in inflow-to-outflow order, with the angular closure
        handled by the mesh's :class:`PoleAngularClosure` strategy
        (Cartesian → identity; curvilinear → Morel-Montry).
        
        Parameters
        ----------
        rhs : AngularFlux
            The right-hand-side source.  Convention:
            
            * ``rhs.values``: per-ordinate source ``(N, ng, nx, ny)``.
            * ``rhs.boundary``: face source (typically zero for volumetric
              sources from S, F, q_ext).
            * ``rhs(1)`` (the previous frame, if present): the angular
              flux at the iterate that GENERATED this source; used as
              the Carlson seed for curvilinear sweeps.  ``None`` (or
              ``history_depth=1``) means cold start: the Carlson seed
              falls back to the in-iteration-source default.
        
        Returns
        -------
        AngularFlux
            The angular flux satisfying ``(L + C) ψ = rhs``.  Carries:
            
            * ``.values``: the swept cell-centre angular flux.
            * ``.boundary``: the sweep's outflow face state (the natural
              face flux at the converged BC).
            * ``.history_depth``: inherited from ``rhs.history_depth``.
            * No history attached (single-frame result).
        """
        from .angular_flux import AngularFlux
        from .sweep import transport_sweep
        from .sources import PerOrdinateSource
        
        # Extract Carlson seed from rhs(1) if present.
        carlson_seed = rhs(1)  # AngularFlux | None
        
        # The sweep API takes (iso_source, sigma_t, sn_mesh, boundary, 
        # aniso_source=..., initial_guess=...).  Treat all of rhs.values 
        # as the per-ordinate aniso_source; iso = zero.  This preserves
        # the per-ordinate variation (no integration to scalar).
        iso_source = self.sn_mesh.zeros_isotropic_source()
        aniso_source = PerOrdinateSource(rhs.values, self.sn_mesh)
        
        # Sweep operates with a write-through boundary buffer; the user-
        # facing AngularFlux is immutable, so we copy boundary state in
        # and read it out.  See "Boundary mutation contract" in the
        # docstring of transport_sweep.
        boundary_buf = self.sn_mesh.zeros_boundary_flux()
        # Copy rhs.boundary into the buffer (typically zero for SI rhs,
        # potentially non-zero for Krylov preconditioner inputs).
        if rhs.boundary is not None:
            _copy_boundary_flux(rhs.boundary, boundary_buf)
        
        angular, _scalar = transport_sweep(
            iso_source, self.sigma, self.sn_mesh, boundary_buf,
            aniso_source=aniso_source,
            initial_guess=carlson_seed,
        )
        
        return AngularFlux(
            angular, self.sn_mesh, boundary_buf,
            history_depth=rhs.history_depth,
        )
```

**Dispatch from `__add__`**:

```python
class StreamingOperator(LinearOperatorMixin):
    # ... existing fields and methods ...
    
    def __add__(self, other):
        """Compose with another LinearOperator.
        
        If ``other`` is a :class:`CollisionOperator`, returns a
        :class:`InvertibleOperator` carrying the sweep-invertible
        algebraic identity (CAP_SOLVE).  Otherwise falls through to
        generic :class:`OperatorSum`.
        """
        if isinstance(other, CollisionOperator):
            return InvertibleOperator(self, other)
        return OperatorSum(self, other)
    
    def __sub__(self, other):
        """Streaming - operator.  No special-case; delegate to OperatorSum + ScaledOperator(-1)."""
        return OperatorSum(self, ScaledOperator(-1.0, other))
```

**Note on `__radd__`**: also implement `CollisionOperator.__radd__(self, streaming)` to handle `C + L` (commutative) so user can write either order.

**The `_copy_boundary_flux` helper**:

```python
def _copy_boundary_flux(src: "BoundaryFlux", dst: "BoundaryFlux") -> None:
    """Copy non-None face buffers from src to dst (same mesh required)."""
    if src.xmin_face is not None and dst.xmin_face is not None:
        dst.xmin_face[...] = src.xmin_face
    if src.xmax_face is not None and dst.xmax_face is not None:
        dst.xmax_face[...] = src.xmax_face
    if src.xmin_xmax_buf is not None and dst.xmin_xmax_buf is not None:
        dst.xmin_xmax_buf[...] = src.xmin_xmax_buf
    if src.ymin_ymax_buf is not None and dst.ymin_ymax_buf is not None:
        dst.ymin_ymax_buf[...] = src.ymin_ymax_buf
```

**Tests** (new file: `tests/sn/test_invertible_operator.py`):
- `test_dispatch_via_add` — `L + C` returns InvertibleOperator
- `test_dispatch_via_radd` — `C + L` also returns InvertibleOperator (commutative)
- `test_dispatch_mismatched_mesh` — raises ValueError
- `test_dispatch_negative_sigma_rejected` — InvertibleOperator with negative σ raises
- `test_apply_equals_operatorsum` — `(L + C).apply(ψ)` ≡ `L.apply(ψ) + C.apply(ψ)` bit-exact
- `test_solve_inverts_apply` — `(L + C).solve((L + C).apply(ψ)) ≈ ψ` for various test problems (slab, sphere, cylinder)
- `test_solve_carlson_seed_from_history` — InvertibleOperator.solve reads `rhs(1)` for the Carlson seed
- `test_solve_cold_start_no_history` — InvertibleOperator.solve with no history falls back to in-iteration-source default (the trace-space sweep we landed in Step 0)
- `test_solve_with_sigma_r` — InvertibleOperator(L, CollisionOperator(σ_r)) computes the correct sweep for the foldable variant. Verify against an analytical 1-G case where σ_r = σ_t - σ_s0.

**Acceptance**: all new tests pass. Existing AngularFlux + iteration tests pass unchanged.

**Commit message**: `feat(sn): InvertibleOperator — typed (L+C) with sweep .solve via __add__ dispatch (R-1 Step 4 C)`

---

### Step D — Carve `_solve_krylov` (1-D only; 2-D NotImplementedError)

**File**: `orpheus/sn/solver.py`

**Goal**: replace the existing `_solve_fixed_source_krylov` (~200 LoC of inline scipy.gmres orchestration + packed-1D RHS builders + `_make_sweep_preconditioner`) with a thin wrapper around `KrylovAcceleration.solve(q_ext_AngularFlux)` on typed operators.

**In `SNSolver.__init__`**:

```python
# Existing: self.scattering_op, self.fission_op constructed
# Existing: self.L = SNStreamingOperator(...)  ← KEEP (other callers use it)

# NEW: typed leaf operators consumed by the carved paths.
self.L_leaf = StreamingOperator(sn_mesh, self.mat_xs.total_cross_section)
self.C_t = CollisionOperator(sn_mesh, self.mat_xs.total_cross_section)
# Composite — automatically returns InvertibleOperator via __add__ dispatch.
self.LC = self.L_leaf + self.C_t   # InvertibleOperator with CAP_SOLVE
```

**Replace `_solve_fixed_source_krylov`** (the inner method that's called per outer eigenvalue step):

```python
def _solve_krylov(self, fission_source: np.ndarray, flux_distribution: np.ndarray) -> np.ndarray:
    r"""Inner within-group solve via :class:`KrylovAcceleration` on typed AngularFlux.
    
    R-1 Step 4 — replaces the inline scipy.gmres orchestration + packed-1D
    RHS builders + identity-on-face preconditioner with the typed operator-
    algebra path:
    
    .. math::
    
        (L_{\rm leaf} + C_t - S - F) \psi = q_{\rm ext}
    
    where ``F = ZeroOperator`` for within-group inner solves.  GMRES runs
    unpreconditioned per R-1 user direction; the σ_r foldable preconditioner
    and the block-inverse face preconditioner are tracked under issue #200.
    
    1-D meshes (slab + sphere + cylinder) only.  2-D Cartesian raises
    :class:`NotImplementedError` — 2-D B1'' face layout work is deferred to
    Phase A.
    """
    if self.sn_mesh.reduced is None:
        raise NotImplementedError(
            "R-1 Step 4 — 2-D Cartesian carve deferred to Phase A.  The "
            "B1'' face block is 1-D-only; 2-D needs a separate 4-face "
            "layout (xmin, xmax, ymin, ymax).  Use a 1-D mesh for now."
        )
    
    from .angular_flux import AngularFlux
    
    N = self.quad.N
    nx, ny, ng = self.sn_mesh.nx, self.sn_mesh.ny, self.ng
    
    # ── Build the typed RHS ─────────────────────────────────────────
    # The fission source enters at the per-ordinate level — same scalar 
    # value broadcast over every ordinate (iso emission).  The convention 
    # is "per-ordinate source rate without /W normalisation" (same as 
    # ScatteringOperator's typed branch).
    q_per_ordinate = np.broadcast_to(
        fission_source[None, :, :, :], (N, ng, nx, ny),
    ).copy()
    q_ext_typed = AngularFlux(q_per_ordinate, self.sn_mesh)
    
    # ── Build the operator triple ───────────────────────────────────
    # InvertibleOperator carries .apply (sum) + .solve (sweep).  We do NOT
    # use .solve as a preconditioner per user direction (issue #200 tracks 
    # adding it); KrylovAcceleration runs unpreconditioned via 
    # preconditioner=None.
    krylov = KrylovAcceleration(
        L=self.LC,
        S=self.scattering_op,
        F=ZeroOperator(),
        preconditioner=None,
        tol=self.inner_tol,
        max_iter=self.max_inner,
        restart=min(50, N * ng * nx * ny),
    )
    
    # ── Warm start from previous iterate ───────────────────────────
    initial_guess = getattr(self, "_psi_typed", None)
    
    psi_typed, _residuals = krylov.solve(q_ext_typed, initial_guess=initial_guess)
    
    self._psi_typed = psi_typed   # cache for next outer step's warm start
    
    # Reduce angular → scalar for the eigenvalue solver's contract.
    return psi_typed.integrate_angular().values
```

**KEEP** (don't delete yet): `self.L` (SNStreamingOperator), `_make_sweep_preconditioner`, `_build_rhs_*`, the legacy `_solve_fixed_source_krylov` module-level function. They have other callers that retire in Step H.

**Acceptance gates**:
- `tests/sn/l1_analytical/test_kinf_homogeneous.py` passes (both inner_solver paths × all coords × all ng)
- `tests/sn/regression/test_dd_regression.py` passes
- `tests/sn/test_b1pp_verification.py` passes
- All 2-D Cartesian eigenvalue tests still pass (they route to `_solve_source_iteration` which is untouched in this step; if any route to `_solve_krylov`, they'll need to switch to `inner_solver="source_iteration"` or wait for Phase A)

**Risk**: GMRES without preconditioning may converge slower than the legacy path (~5-10× more iterations expected). Allow more `max_inner` budget; the L1 tests use tight tolerances and may need higher `max_inner`. **Tune via the failing-test feedback loop, not by pre-guessing.**

**Commit message**: `refactor(sn): _solve_krylov via KrylovAcceleration on typed AngularFlux (R-1 Step 4 D)`

---

### Step E — Carve `_solve_source_iteration` similarly

**File**: `orpheus/sn/solver.py`

**Goal**: replace the inline SI loop with `SourceIteration.solve(q_ext_AngularFlux)`.

```python
def _solve_source_iteration(self, fission_source, flux_distribution) -> np.ndarray:
    r"""Inner within-group solve via :class:`SourceIteration` on typed AngularFlux.
    
    R-1 Step 4 — the iteration ψ_(n+1) = L.solve(F·ψ_n + S·ψ_n + q_ext)
    where L = self.LC (InvertibleOperator) and L.solve IS the WDD sweep.
    
    The previous-iterate Carlson seed is threaded into the rhs's history 
    by the SourceIteration primitive — InvertibleOperator.solve reads 
    ``rhs(1)`` to populate the sweep's ``initial_guess``.
    """
    if self.sn_mesh.reduced is None:
        raise NotImplementedError(
            "R-1 Step 4 — 2-D Cartesian SI carve deferred to Phase A."
        )
    
    from .angular_flux import AngularFlux
    
    N = self.quad.N
    nx, ny, ng = self.sn_mesh.nx, self.sn_mesh.ny, self.ng
    
    q_per_ordinate = np.broadcast_to(
        fission_source[None, :, :, :], (N, ng, nx, ny),
    ).copy()
    q_ext_typed = AngularFlux(q_per_ordinate, self.sn_mesh)
    
    si = SourceIteration(
        L=self.LC,
        S=self.scattering_op,
        F=ZeroOperator(),
        tol=self.inner_tol,
        max_iter=self.max_inner,
    )
    
    initial_guess = getattr(self, "_psi_typed", None)
    psi_typed, _residuals = si.solve(q_ext_typed, initial_guess=initial_guess)
    
    self._psi_typed = psi_typed
    return psi_typed.integrate_angular().values
```

**Acceptance gates**: same L1 tests as Step D, with `inner_solver="source_iteration"`. The Carlson seed plumbing via `rhs(1)` is the load-bearing change here — verify with a curvilinear MG test (e.g. `test_kinf_homogeneous[sphere-2eg]`).

**Risk**: If Carlson seed plumbing is wrong, curvilinear SI regresses. Mitigation: add a focused test (`test_carlson_seed_threaded_via_history`) that constructs a curvilinear SourceIteration solve and verifies bit-identity against the legacy path on a small problem.

**Commit message**: `refactor(sn): _solve_source_iteration via SourceIteration on typed AngularFlux (R-1 Step 4 E)`

---

### Step F — Migrate packed-vector tests to typed AngularFlux

**Files**: ~10 test files, ~70 reference sites:
- `tests/sn/test_collision_operator.py` — 12 sites (currently fails because packed sizes mismatch; this step fixes)
- `tests/sn/test_streaming_operator_decomposition.py` — 4 sites
- `tests/sn/test_phase_c_gates.py` — 11 sites
- `tests/sn/test_b1pp_verification.py` — 5 sites
- `tests/numerics/test_iteration.py` — 1 site
- `tests/sn/test_snstreamingoperator.py` — 16 sites (will be deleted in Step H)
- `tests/sn/test_streaming_operator.py` — 8 sites
- `tests/sn/test_unified_matvec_*.py` — 7 sites total

**Pattern**:
```python
# OLD:
rng = np.random.default_rng(seed)
psi = rng.standard_normal(eq_map.n_unknowns)   # packed-1D
result = op.apply(psi)
# ... assertions on result ...

# NEW:
rng = np.random.default_rng(seed)
psi_typed = AngularFlux.from_flat_with_traces(
    rng.standard_normal(eq_map.n_unknowns), sn_mesh,
)
result_typed = op.apply(psi_typed)
result_flat = result_typed.to_flat_with_traces()
# ... same assertions on result_flat ...
```

For tests that exercise multiple operators, the typed surface is cleaner — direct AngularFlux throughout, only flatten at the test-assertion boundary.

**Acceptance**: zero `eq_map.n_unknowns`-sized packed vectors fed to `op.apply(np.ndarray)` remain. The 39 pre-existing CollisionOperator test failures resolve mechanically.

**Commit message**: `test(sn): migrate packed-vector test sites to typed AngularFlux (R-1 Step 4 F)`

---

### Step G — Retirement: kill SNStreamingOperator + EquationMap + legacy decoders

**Files**: `orpheus/sn/operator.py`, `orpheus/sn/solver.py`, ~10 test files

**Symbols to retire** (cumulative LoC ≈ 1100):

In `orpheus/sn/operator.py`:
- `SNStreamingOperator` class (bundled L+C with no B1'' face — sunset since Phase G)
- `EquationMap` class (~170 LoC including module-scope helpers)
- `build_equation_map` (2-D Cartesian + legacy 1-D)
- `build_equation_map_spherical` (legacy compressed sphere)
- `build_equation_map_cylindrical` (legacy compressed cylinder)
- `solution_to_angular_flux` (legacy decoder)
- `solution_to_angular_flux_spherical` (legacy compressed decoder)
- `solution_to_angular_flux_cylindrical` (legacy compressed decoder)

In `orpheus/sn/solver.py`:
- `_solve_fixed_source_krylov` (module-level MMS helper — migrate or retire)
- `_solve_fixed_source_si` (module-level MMS helper — migrate or retire)
- `_make_sweep_preconditioner` (method)
- `_build_rhs_cartesian`, `_build_rhs_spherical`, `_build_rhs_cylindrical` (free functions)
- `SNSolver.L` attribute (the SNStreamingOperator) — replace with `self.LC` (InvertibleOperator)

**Move** (not retire — relocate to private methods on AngularFlux):
- `pack_with_traces`, `solution_to_angular_flux_with_traces`, `build_equation_map_with_traces` → privatized inside `angular_flux.py` as the implementation of `to_flat_with_traces` / `from_flat_with_traces`.

**Tests to retire**:
- `tests/sn/test_snstreamingoperator.py` (delete entirely — operator is gone)
- Per-step gate: if any test imports a retired symbol, that test gets migrated to the typed equivalent OR retired if redundant.

**Acceptance**:
- `grep -rn "SNStreamingOperator\|EquationMap\|build_equation_map\|solution_to_angular_flux\|_make_sweep_preconditioner\|_build_rhs_" orpheus/ tests/` returns empty (or only in retired-history docstrings).
- All test suites pass.
- LoC reduction: ~1100 retired; ~200 moved into AngularFlux as private helpers. Net retirement: ~900 LoC.

**Commit message**: `refactor(sn): retire SNStreamingOperator + EquationMap + legacy decoders (R-1 Step 4 G)`

---

### Step H — Sphinx documentation + closeout

**Files**: `docs/theory/sn_*.rst`, `docs/development.rst`

**Required**:
- Update `docs/theory/sn_operator_algebra.rst` (or equivalent) to describe the typed-AngularFlux algebra end-to-end. Include:
  - The four operators (L, C, S, F) with typed apply signatures
  - `InvertibleOperator` and the algebraic identity it encodes
  - The ravellable protocol (`to_flat_with_traces` / `from_flat_with_traces`)
  - The AngularFlux history container and shift-register semantics
- Cross-references to issues #200 (preconditioner) and #201 (dimensional typing).
- Memory entries in `.claude/agent-memory/method-implementer/` for the carve approach.

**Acceptance**:
- `sphinx-build docs docs/_build/html` succeeds with no warnings.
- The new operator-algebra page reads as a textbook chapter (full derivations, gotchas, design rationale).
- ERR-NNN entries in `error_catalog.md` for any new bug classes caught during the carve.

**Commit message**: `docs(sn): operator-algebra narrative + closeout (R-1 Step 4 H)`

---

## Order & sequencing

Strict order:
1. **A** (AngularFlux history) — pure addition; tests unaffected
2. **B** (drop inverter) — refactor iteration primitives; test migration in-place
3. **C** (InvertibleOperator) — new class; no consumer changes yet
4. **D** (carve _solve_krylov) — first real architectural change; gate on L1 tests
5. **E** (carve _solve_source_iteration) — relies on Carlson-seed-via-history (A) + L.solve (C) + dropped inverter (B)
6. **F** (test migration) — packed-vector tests retire
7. **G** (retirement) — once no consumers remain, delete
8. **H** (docs) — narrate the architecture

Each step is a separate commit. Each gate is the existing test suite + new tests for that step.

## Verification gates — consolidated

```bash
# Per-step:
.venv/bin/python -m pytest tests/sn/test_angular_flux_with_boundary.py -q              # Step A
.venv/bin/python -m pytest tests/numerics/ -q                                          # Step B
.venv/bin/python -m pytest tests/sn/test_invertible_operator.py -q                     # Step C
.venv/bin/python -m pytest tests/sn/l1_analytical/test_kinf_homogeneous.py -q -k krylov  # Step D
.venv/bin/python -m pytest tests/sn/l1_analytical/test_kinf_homogeneous.py -q -k source  # Step E
.venv/bin/python -m pytest tests/sn/test_collision_operator.py tests/sn/test_streaming_operator.py -q  # Step F
.venv/bin/python -m pytest tests/sn/ -q                                                # Step G full sweep

# End-to-end:
.venv/bin/python -m pytest tests/ -q                                                   # Step H final
```

## Risks & rollback

| Step | Risk | Mitigation | Rollback |
|---|---|---|---|
| A | History container API churn | Default values match current contract; new fields are additive | Revert single commit |
| B | Test migration drift | Wrap inverters in tiny operator classes mechanically | Revert + replay |
| C | InvertibleOperator.solve semantics wrong | Bit-identity test against transport_sweep direct call | Revert single commit |
| D | GMRES convergence slow | Increase max_inner; if pathological, file issue and consider re-enabling identity-on-face preconditioner | Keep legacy path one more commit cycle |
| E | Carlson seed plumbing wrong on curvilinear | Targeted test reproducing the legacy SI path bit-identical | Revert single commit |
| F | Test migration boring but error-prone | Mechanical pattern; review diff systematically | Revert per-file as needed |
| G | Hidden caller of retired symbol | Grep before each delete; CI catches missed imports | Revert single delete; re-grep |
| H | Sphinx build breaks | Build locally before commit | Fix in place |

## Key files

| File | What changes |
|---|---|
| `orpheus/sn/angular_flux.py` | A: add history container + dunders |
| `orpheus/numerics/iteration.py` | B: retire inverter; rename for Krylov |
| `orpheus/sn/operator.py` | C: add InvertibleOperator + __add__ dispatch; G: retire SNStreamingOperator + EquationMap + legacy decoders |
| `orpheus/sn/solver.py` | D, E: carve _solve_krylov, _solve_source_iteration; G: retire _build_rhs_*, _make_sweep_preconditioner, module-level _solve_fixed_source_* |
| Test files (~10) | F: packed-vector migration |
| Sphinx docs | H: narrate |

## Cross-references

- Issue [#200](https://github.com/deOliveira-R/ORPHEUS/issues/200) — block-inverse Krylov preconditioner (deferred from R-1)
- Issue [#201](https://github.com/deOliveira-R/ORPHEUS/issues/201) — dimensional typing of AngularFlux / AngularSource / AngularResidual
- Closeout memo from bad agent attempt: `.claude/agent-memory/method-implementer/issue_197_r1_step4b_closeout.md`
- Diagnostic scripts (preserved): `scratch/r1_step4_*.py`

## After R-1 Step 4 completes — pending work

Tracked elsewhere; not in this plan:
- Phase R-2a/b — (n_mask, ng) flip + F-order packed-vec retirement (task #32)
- Phase R-4 — Memory growth investigation (Issue #198, task #34)
- Phase A — Algebra completion + retirement (Issue #199, task #36)
- Consolidated Phase BC — BC verification + cleanup (task #10)
- Consolidated Phase H — Adjoint on leaves + solve_sn_adjoint (task #11)
- Consolidated Phase D — Sphinx narrative (task #12)
- Consolidated Phase V — V&V & catalog hygiene (task #13)
