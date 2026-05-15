---
name: typed-field-contracts-phase-g
description: Canonical typed data-structure contract for the SN four-operator algebra (L + C − S − F/k)ψ = q. Frozen dataclass signatures for AngularFlux/ScalarFlux/Solution + relationship to Issue #197. Synthesis of the grand report v3 vocabulary with the c21c2ef leaf state.
metadata:
  type: project
---

# Typed Field Contracts for Phase G — `(L + C − S − F/k)ψ = q`

## Executive summary (5 lines)

1. **Canonical typed fields:** `AngularFlux(values=(N, nx, ny, ng))`, `ScalarFlux(values=(nx, ny, ng))`, `IsotropicSource`, `PerOrdinateSource`, `Solution(keff, angular_flux, scalar_flux, iter_history)` — all frozen dataclasses carrying a numpy array + SNMesh ref; arithmetic dunders delegate to numpy so block ops stay vectorised.
2. **Storage layout:** energy `g` is **always the trailing axis** (block over ng); ordinate `n` is the leading axis when present. Mesh axes (`nx, ny`) live between. One numpy op multiplies every field across every group; no per-group Python loops.
3. **Operator-algebra coupling:** every leaf becomes `apply(psi: AngularFlux) → AngularFlux`. F is `(internally) AngularFlux → ScalarFlux (project) → ScalarFlux (rank-1 ops) → AngularFlux (broadcast)`. `(L + C − S_foldable − F/k).apply(ψ)` distributes cleanly because the four operands agree on a typed domain; the existing `OperatorSum.apply: a.apply(x) + b.apply(x)` works unchanged.
4. **Issue #197 partial close:** `AngularFlux`, `ScalarFlux`, `IsotropicSource`, `PerOrdinateSource` land as **structural dataclasses** (more capable than bare NewType — they carry mesh ref + dunders); the remaining types from #197 Wave 1 (`CrossSection`, `ReactionRate`, `GroupRate`, `Volume`, `Position`, `KEff`, `EmissionSpectrum`, etc.) stay as bare `NewType` per the original plan.
5. **Migration impact on `c21c2ef`:** the three b.i v2 leaves keep their constructors (`StreamingOperator(sn_mesh, sigma_t)`, `CollisionOperator(sn_mesh, sigma)`, `ScatteringOperator.foldable_part()` API), but their `apply` rewires to `AngularFlux → AngularFlux`; the bodies switch from `transport_operator_matvec_*` (which hardcodes a discretisation) to `sn_mesh.cell_update.residual` folded through `sn_mesh.dag_walk(...)` — making the discretisation strategy-driven per the user's hard requirement. Tests in `tests/sn/test_streaming_operator{,_decomposition}.py` and `tests/sn/test_collision_operator.py` are mechanically rewritten (build typed input → assert typed output); the algebraic identity tests survive verbatim.

---

## §1 Inventory from the grand report

This section enumerates every solver-state data structure the grand report v3 names, with units, shape, and existing-codebase counterpart.

### 1.1 Field hierarchy (§5.5, lines 577–605)

| Report name | Physical meaning | Units | Shape | Existing counterpart | Notes |
|---|---|---|---|---|---|
| `AngularFlux` | ψ(r, Ω, g) — phase-space directional flux | 1/(cm²·s·sr) | `(N, nx, ny, ng)` | Bare `np.ndarray` returned by `transport_sweep` | Hard requirement: structural type |
| `ScalarFlux` | φ(r, g) = ∫ψ dΩ — angle-integrated flux | 1/(cm²·s) | `(nx, ny, ng)` | Bare `np.ndarray` returned by `transport_sweep` | Hard requirement: structural type |
| `Field` (base) | Generic abstract base — carries `space`, `data`, `representation` | — | — | None | Optional ancestor for §32.5 (line 5919) |
| `CoefficientField` | Coefficients in a basis (PN moments, Galerkin DoFs) | — | basis-dependent | Implicit in `HarmonicMomentProjection` output | Not needed for SN |
| `DiscreteField` | Per-cell discrete values | — | `(nx, ny, ng)` | The ScalarFlux above | Same data shape |
| `MeshField` | Field tagged with mesh (cell-centered, face-valued) | — | varies | `sn_mesh.volumes` etc. | Mesh side, not solution state |
| `TraceField` | Field on `Γ_-` or `Γ_+` | varies | `(N_inflow, ng)` per face | `psi_bc["bc_*"]` dict entries (line 282-289 sweep.py) | Boundary trace storage |
| `HarmonicMomentField` | φ_ℓm coefficients | 1/(cm²·s) | `(L+1, 2L+1, nx, ny, ng)` | `HarmonicMomentProjection.apply` output (orpheus/sn/scattering.py:212-251) | Internal to S |
| `GroupFlux` | Flux on one energy group | 1/(cm²·s) | `(N, nx, ny)` or `(nx, ny)` | Slice of AngularFlux/ScalarFlux | Not a separate type |
| `StochasticField` | UQ extension | varies | + `Xi` axis | None | Future |
| `TensorTrainField` | TT-compressed flux | varies | TT-cores | None | Future |
| `ParticleEnsemble` | MC weighted sample bank | — | `(n_particles,)` | None | MC, separate solver family |

### 1.2 Source/RHS vocabulary (§4, lines 41-77; §5.5, line 583)

| Report name | Meaning | Units | Shape | Existing counterpart |
|---|---|---|---|---|
| `q` (in `(L + C − S − F)ψ = q`) | External (driven) source | 1/(cm³·s·sr) | `(N, nx, ny, ng)` if anisotropic; `(nx, ny, ng)` if isotropic | Sweep arg `Q` (isotropic) + `Q_aniso` (per-ordinate) — see `transport_sweep(Q, sig_t, sn_mesh, psi_bc, Q_aniso)` at orpheus/sn/sweep.py:96 |
| `Source` (§28, line 5382) | Boundary or volumetric source | — | — | Implicit |
| Per-ordinate source `Q^aniso_n(r)` | RHS of the per-ordinate within-group system; rank ≥ 2 in angle | 1/(cm³·s·sr) | `(N, nx, ny, ng)` | `Q_aniso` in `transport_sweep`; output of `ScatteringOperator.build_aniso_source` (orpheus/sn/scattering.py:397-464) |
| Isotropic source | `(P0+(n,2n))` source per group, direction-independent | 1/(cm³·s·sr) (broadcast) | `(nx, ny, ng)` | `Q` (isotropic part), `Q_iso` in `ScatteringOperator.apply` (orpheus/sn/scattering.py:742) |
| `ResidualSource` (§25.1, line 4614) | `r = q − Aψ_D` for hybrid corrections | 1/(cm³·s·sr) | matches q | None — Phase G doesn't need this |
| `BoundarySource` (§16A.2, lines 2788) | Prescribed inflow at Γ_- | 1/(cm²·s·sr) | `(N_inflow,)` per face | Implicit in BC-applied face buffers |

### 1.3 Rates and tallies (§5.5, §5.6 lines 612-625; §27 lines 5144-5240)

| Report name | Meaning | Units | Shape | Existing counterpart |
|---|---|---|---|---|
| `ReactionRate` | σ·φ — collision/absorption/fission rate density | 1/(cm³·s) | `(nx, ny, ng)` | Computed inline in `compute_group_*_rate` (orpheus/sn/solver.py:334-385) |
| `GroupRate` | ∫σ·φ dV — volume-integrated per group | 1/s | `(ng,)` | Return of `compute_group_production_rate` / `compute_group_absorption_rate` |
| `TallyCosheaf` (§15A.6 line 2284) | Additive aggregation of scores | varies | — | None — MC concept |
| `CurrentCochain` (§15A.10) | Face-summed currents | n/cm²/s | `(N_faces, ng)` | None |
| `Functional` / `ReactionRateFunctional` (§5.7 line 633, §13 line 1892) | Map field → scalar response | — | — | `compute_keff` is a degenerate case |

### 1.4 Iteration state (§21.2 lines 4163-4188, §36.1 lines 6705-6743)

| Report name | Meaning | Units | Shape | Existing counterpart |
|---|---|---|---|---|
| `keff` (k) | Multiplication eigenvalue | dimensionless | scalar | `SNResult.keff: float` (orpheus/sn/solver.py:109) |
| `keff_history` | Outer iteration trajectory | dimensionless | `list[float]` | `SNResult.keff_history` (orpheus/sn/solver.py:110) |
| `Eigenpair` (§21.5 lines 4253-4257) | (value, right, left, residual_norm) | varies | — | Implicit — SNResult is a degenerate Eigenpair |
| `EigenSpectrum` (§21.5 lines 4259-4269) | Multiple Eigenpairs | — | — | None — Phase G is single-mode |
| `ResidualHistory` (Trefethen & Bau §3.2) | Per-iter Krylov residual norm | dimensionless | `list[float]` | None today; could land with Solution |
| `DominanceRatio` (§21.2 line 4188, §22) | `|k₁/k₀|` convergence quotient | dimensionless | scalar | None |
| `n_inner` / iter counts | Inner solver iteration count | int | scalar | `SNResult` lacks; `SNFixedSourceResult.n_inner` has it |
| `Diagnostic` (§39 lines 6989-7029) | Named conservation/positivity/etc. defects | varies | scalar | Distributed across vv-harness; no Solution-side carrier |

### 1.5 Solution-class proper (§29.2 lines 5503-5525, §31, §38)

The grand report does NOT use a single literal "Solution" class name in §29.2 — instead, the API sketch is `k, psi = PowerIteration().solve(problem)` with the result being a tuple. But §21.5 / §29.2 / §38 imply the canonical container:

| Report name | Meaning | Existing counterpart |
|---|---|---|
| `Eigenpair(value, right, left, residual_norm)` | The k-eigenvalue solution container | `SNResult` |
| (implicit) `(scalar_flux, angular_flux, mesh, iteration_history)` bundle | Fixed-source / eigenvalue result | `SNFixedSourceResult` (orpheus/sn/solver.py:82) + `SNResult` (orpheus/sn/solver.py:101) |

**Recommendation:** name the new class `Solution` (matching the user's verbatim ask). Place under `orpheus/sn/solution.py` so it can hold SN-specific iteration metadata, while leaving a `BaseSolution` Protocol surface in `orpheus/numerics/solution.py` for future MoC/CP analogues.

### 1.6 Operator vocabulary (§1 lines 25-110, §5.7 lines 636-672) — already implemented in ORPHEUS

| Report name | Code name | File:line |
|---|---|---|
| `L` — streaming + collision (legacy bundle) | `SNStreamingOperator` | orpheus/sn/operator.py:1120 |
| `L` — pure streaming (Resolution A) | `StreamingOperator` | orpheus/sn/operator.py:1551 |
| `C` — collision | `CollisionOperator` | orpheus/sn/operator.py:1748 |
| `S` — scattering | `ScatteringOperator` | orpheus/sn/scattering.py:255 |
| `S.foldable_part()` — P0 within-group self-scatter | `ScatteringOperator.foldable_part` | orpheus/sn/scattering.py:503 |
| `S.residual_part()` — cross-group/Pℓ/n2n | `ScatteringOperator.residual_part` | orpheus/sn/scattering.py:550 |
| `F` — fission | `FissionOperator` | orpheus/sn/fission.py:83 |
| `T` — time mass | (not yet) | — |
| `K` = `A_loss⁻¹ F` — multiplication operator | Implicit in `SNSolver.solve_eigenvalue` | — |
| Algebra: `OperatorSum`, `OperatorProduct`, `ScaledOperator` | Same | orpheus/numerics/operator.py:489, 568, 641 |
| Adapter to scipy | `as_scipy_linop` | orpheus/numerics/operator.py:1318 |

### 1.7 Boundary-trace vocabulary (§16A.1, §16A.5)

| Report name | Existing counterpart | File:line |
|---|---|---|
| `BoundaryTraceLaw` | Same | orpheus/geometry/boundary/*.py |
| `BoundaryRealizer` | `SNBoundaryRealizer` | orpheus/sn/boundary_realizer.py |
| `InflowTraceSpace`, `OutflowTraceSpace` | Same | orpheus/numerics/trace_space.py |
| `VacuumInflow` / `ZeroInflowTrace` | `VacuumInflow` | orpheus/geometry/boundary/vacuum.py |
| `ReflectiveBoundary` (physical) | Same | orpheus/geometry/boundary/reflective.py |
| `PermutationOperator`, `IncomingOrdinateMaskTensor` | Same | orpheus/numerics/operator.py:739, 843 |

### 1.8 Diagnostic and historical state (§32 lines 5721-5749, §39 lines 6993-7029)

| Report name | Meaning | Existing counterpart |
|---|---|---|
| `Axis(name, size, coordinate, measure)` | Labeled-axis primitive | None — bare numpy shapes |
| `AxisProduct` | Tuple of `Axis` | None |
| `DomainMismatchError`, `CodomainMismatchError` | Operator algebra type errors | `IncompatibleOperatorComposition` (orpheus/numerics/operator.py:102) |
| `ConservationDefect`, `PositivityViolation` | V&V residuals | Distributed across `tests/_harness/` |

**Recommendation:** the Phase G typed fields do NOT need `Axis`/`AxisProduct` yet — the report uses them aspirationally for §32, while the SN production stays on plain numpy. Defer.

---

## §2 Proposed canonical type hierarchy

### 2.1 Design tenets

Driven by the user's verbatim requirements:

1. **Frozen dataclasses + numpy array, not bare NewType.** A NewType is a label; a structural type carries operations (point-lookup, block-multiply, mesh ref).
2. **Block storage:** energy `g` always trails; `numpy` block ops just work — `flux_a + flux_b`, `2.0 * flux`, `sigma_t * flux.values` are one-step operations.
3. **Mesh ref always attached.** Every field knows where its values live, so plotting can ask `flux.at(position)` without threading the mesh through every call site.
4. **Pattern 4 (illegal-states-unrepresentable):** an `AngularFlux` paired with a mesh that has a different `(nx, ny)` is impossible — the constructor validates.
5. **`dataclasses.replace()` for mutation, not `replace_psi`.** Standard Python pattern.
6. **No `WithinGroupFlux` per-group dataclass** (anti-recommendation #4).
7. **Performance-neutral:** dunders delegate to numpy; no per-element overhead.

### 2.2 `AngularFlux` — phase-space directional flux

```python
# orpheus/sn/typed_fields.py  (NEW)

from __future__ import annotations
from dataclasses import dataclass, field, replace
from typing import TYPE_CHECKING
import numpy as np

if TYPE_CHECKING:
    from orpheus.sn.geometry import SNMesh


@dataclass(frozen=True, slots=True)
class AngularFlux:
    r"""Angular flux ψ(r, Ω, g) over an SN method space.

    Block-storage layout: leading axis is ordinate index n; trailing
    axis is energy group g. The interior (nx, ny) carries the spatial
    mesh. **Energy is trailing by design** — a per-ordinate per-cell
    block operator that mixes groups (e.g. `Σ_s,0[g→g'] · ψ`) acts as
    a matmul on the last axis with no reshape.

    Units: 1/(cm²·s·sr).

    Constructor invariants
    ----------------------
    * ``values.shape == (sn_mesh.quad.N, sn_mesh.nx, sn_mesh.ny, ng)``
      for some ``ng ≥ 1``. Validated at construction; mismatched shape
      raises ValueError (Pattern 4 — illegal states unrepresentable).
    * ``sn_mesh`` is held by reference, NOT copied. Multiple flux
      objects sharing the same mesh costs zero memory.

    Block-operation dunders
    -----------------------
    ``__add__`` / ``__sub__`` / ``__neg__`` — element-wise on
    ``values``; returned object has the same ``sn_mesh``. The two
    operands MUST have the same ``sn_mesh`` (identity check); a
    ``TypeError`` is raised otherwise — fluxes on different meshes
    are not algebraically combinable.

    ``__mul__`` / ``__rmul__`` — scalar or numpy-broadcastable array
    multiply. ``2.0 * flux`` returns ``AngularFlux``; ``sigma_t *
    flux`` (where ``sigma_t.shape == (nx, ny, ng)``) broadcasts
    correctly across the leading ordinate axis without manual
    reshape.

    ``__truediv__`` — scalar division for normalisation
    (``F.apply(psi) / k``).

    Conversion
    ----------
    :meth:`to_scalar` returns ``ScalarFlux``: ``φ_g(r) = Σ_n w_n
    ψ_n(r, g)``. Pure function — both operands stay untouched.

    Point lookup
    ------------
    :meth:`at` returns a ``(N, ng)`` value array at a given (x, y)
    cell index. The mesh's ``centers`` / ``edges`` are consulted via
    ``sn_mesh.mesh``. Future extension: spatial interpolation.

    Cross-references
    ----------------
    * Grand Report v3 §5.5 line 584 (``AngularFlux`` in the Field
      hierarchy); §1 Table line 32 (the symbol L acts on this type);
      §9 line 1185 (the SN cell/ordinate/group block layout); §29.2
      line 5524 (``k, psi = PowerIteration().solve(problem)`` — the
      ``psi`` is this type).
    * orpheus/sn/sweep.py:96 — ``transport_sweep`` returns the bare
      values today (``(N, nx, ny, ng)`` tuple of ``(angular_flux,
      scalar_flux)``).
    * orpheus/sn/scattering.py:695 — ``ScatteringOperator.apply``
      consumes ``(N, nx, ny, ng)`` today.

    Parameters
    ----------
    values : np.ndarray, shape (N, nx, ny, ng)
        The actual flux array.
    sn_mesh : SNMesh
        The augmented geometry whose quadrature / cells / boundary
        conditions this flux is defined on.
    """

    values: np.ndarray
    sn_mesh: "SNMesh"

    def __post_init__(self) -> None:
        expected = (
            self.sn_mesh.quad.N,
            self.sn_mesh.nx,
            self.sn_mesh.ny,
            self.values.shape[-1] if self.values.ndim == 4 else -1,
        )
        if self.values.shape[:3] != expected[:3] or self.values.ndim != 4:
            raise ValueError(
                f"AngularFlux.values.shape={self.values.shape} does "
                f"not match (N, nx, ny, ng) for the given sn_mesh "
                f"(N={self.sn_mesh.quad.N}, nx={self.sn_mesh.nx}, "
                f"ny={self.sn_mesh.ny})."
            )

    # ── Properties ────────────────────────────────────────────────
    @property
    def ng(self) -> int:
        return int(self.values.shape[-1])

    @property
    def shape(self) -> tuple[int, int, int, int]:
        return self.values.shape  # type: ignore[return-value]

    # ── Dunders (delegate to numpy; mesh-identity check) ──────────
    def _assert_same_mesh(self, other: "AngularFlux") -> None:
        if self.sn_mesh is not other.sn_mesh:
            raise TypeError(
                "AngularFlux arithmetic requires identical sn_mesh; "
                "cross-mesh combination is not algebraically defined."
            )

    def __add__(self, other: "AngularFlux") -> "AngularFlux":
        if not isinstance(other, AngularFlux):
            return NotImplemented
        self._assert_same_mesh(other)
        return replace(self, values=self.values + other.values)

    def __sub__(self, other: "AngularFlux") -> "AngularFlux":
        if not isinstance(other, AngularFlux):
            return NotImplemented
        self._assert_same_mesh(other)
        return replace(self, values=self.values - other.values)

    def __neg__(self) -> "AngularFlux":
        return replace(self, values=-self.values)

    def __mul__(self, other) -> "AngularFlux":
        # Scalar OR broadcastable ndarray.
        if isinstance(other, AngularFlux):
            return NotImplemented  # Hadamard products are not algebra
        return replace(self, values=self.values * other)

    def __rmul__(self, other) -> "AngularFlux":
        return self.__mul__(other)

    def __truediv__(self, scalar: float) -> "AngularFlux":
        return replace(self, values=self.values / scalar)

    # ── Conversion ────────────────────────────────────────────────
    def to_scalar(self) -> "ScalarFlux":
        r"""Reduce over ordinates: ``φ_g(r) = Σ_n w_n ψ_n(r, g)``.

        One numpy einsum; vectorised over (nx, ny, ng). Pure function.
        """
        phi = np.einsum(
            "n,nxyg->xyg", self.sn_mesh.quad.weights, self.values,
        )
        return ScalarFlux(values=phi, sn_mesh=self.sn_mesh)

    # ── Point lookup ──────────────────────────────────────────────
    def at(
        self, *, ix: int, iy: int = 0, n: int | None = None,
    ) -> np.ndarray:
        r"""Value at a given (cell, ordinate) slot.

        Returns ``(ng,)`` if ``n`` supplied; ``(N, ng)`` otherwise.
        """
        if n is None:
            return self.values[:, ix, iy, :]
        return self.values[n, ix, iy, :]

    @property
    def positions(self) -> np.ndarray:
        """Cell-centre positions, shape (nx,) for 1-D or (nx, ny) for 2-D.

        Forwards to ``sn_mesh.mesh.centers`` (orpheus/geometry/mesh.py:304).
        """
        return self.sn_mesh.mesh.centers  # type: ignore[union-attr]
```

### 2.3 `ScalarFlux` — angle-integrated flux

```python
@dataclass(frozen=True, slots=True)
class ScalarFlux:
    r"""Scalar flux φ(r, g) = ∫_{4π} ψ(r, Ω, g) dΩ over an SN method space.

    Block-storage layout: spatial axes lead, energy trails. The same
    block-multiply convention as :class:`AngularFlux`: ``sigma_t *
    phi.values`` (with ``sigma_t.shape == (nx, ny, ng)``) is one
    elementwise op.

    Units: 1/(cm²·s).

    Constructor invariants
    ----------------------
    ``values.shape == (sn_mesh.nx, sn_mesh.ny, ng)`` for some ``ng ≥
    1``. Validated at construction.

    Block-operation dunders
    -----------------------
    Same conventions as :class:`AngularFlux` — element-wise add/sub,
    scalar/array multiply, scalar divide; mesh-identity check on
    binary ops.

    Conversion
    ----------
    :meth:`broadcast_to_ordinates` returns ``AngularFlux``:
    ``ψ_n(r, g) = φ(r, g)`` for every n — the isotropic embedding
    used when an operator that wants to project a scalar response
    back into the angular phase space (e.g. ``F`` returning to an
    `AngularFlux` codomain in §3.2 below).

    Cross-references
    ----------------
    * Grand Report v3 §5.5 line 585 (``ScalarFlux``).
    * orpheus/sn/sweep.py:96 — second element of the
      ``transport_sweep`` return tuple.
    * orpheus/sn/fission.py:137 — ``FissionOperator.apply``
      consumes this shape today.

    Parameters
    ----------
    values : np.ndarray, shape (nx, ny, ng)
    sn_mesh : SNMesh
    """

    values: np.ndarray
    sn_mesh: "SNMesh"

    def __post_init__(self) -> None:
        expected = (self.sn_mesh.nx, self.sn_mesh.ny)
        if self.values.ndim != 3 or self.values.shape[:2] != expected:
            raise ValueError(
                f"ScalarFlux.values.shape={self.values.shape} does not "
                f"match (nx, ny, ng) for the given sn_mesh "
                f"(nx={self.sn_mesh.nx}, ny={self.sn_mesh.ny})."
            )

    @property
    def ng(self) -> int:
        return int(self.values.shape[-1])

    # Dunders identical to AngularFlux (omitted for brevity).
    # _assert_same_mesh, __add__, __sub__, __neg__, __mul__, __rmul__,
    # __truediv__ — same code with the type swap.

    def broadcast_to_ordinates(self) -> "AngularFlux":
        r"""Isotropic embedding: ``ψ_n(r, g) = φ(r, g)`` for all n.

        Returns an :class:`AngularFlux` whose values are the scalar
        flux broadcast across the ordinate axis. Pure function —
        the new instance owns its array (``np.broadcast_to`` + copy
        to make it writable; or store a view if the consumer
        promises read-only — Pattern 2 dual-view, future).
        """
        broadcast = np.broadcast_to(
            self.values[None, :, :, :],
            (self.sn_mesh.quad.N, *self.values.shape),
        ).copy()
        return AngularFlux(values=broadcast, sn_mesh=self.sn_mesh)

    def at(self, *, ix: int, iy: int = 0) -> np.ndarray:
        """Value at a given (cell) slot, shape ``(ng,)``."""
        return self.values[ix, iy, :]
```

### 2.4 `IsotropicSource` and `PerOrdinateSource` — RHS of `(L+C-S-F)ψ = q`

The current sweep splits the source pair `(Q, Q_aniso)` deliberately (orpheus/sn/sweep.py:97-101). The typed contract surfaces this split — `IsotropicSource` is direction-independent (one number per cell per group; broadcast across ordinates inside the sweep); `PerOrdinateSource` is the genuine per-ordinate Pℓ≥1 / boundary contribution. The pair lets the inner sweep avoid a wasteful broadcast.

```python
@dataclass(frozen=True, slots=True)
class IsotropicSource:
    r"""Direction-independent volumetric source ``q(r, g)``.

    Units: 1/(cm³·s·sr). Storage layout matches :class:`ScalarFlux`
    so a scalar-flux-derived source (the P0 + (n,2n) contribution
    in ``ScatteringOperator.add_iso_source`` /
    ``add_n2n_source`` at orpheus/sn/scattering.py:363-395) plugs
    directly in.

    Algebra: same dunders as :class:`ScalarFlux`. The naming makes
    the consumer's intent explicit (this is a SOURCE term, not a
    FLUX) — same shape, different mathematical role (Pattern 3:
    named intermediates with domain semantics).
    """

    values: np.ndarray
    sn_mesh: "SNMesh"
    # __post_init__, ng, dunders — same as ScalarFlux.


@dataclass(frozen=True, slots=True)
class PerOrdinateSource:
    r"""Per-ordinate anisotropic source ``q_n(r, g)`` — RHS of the
    Pℓ≥1 in-scatter sweep.

    Units: 1/(cm³·s·sr). Storage layout matches :class:`AngularFlux`.

    Output of ``ScatteringOperator.build_aniso_source`` (orpheus/sn/
    scattering.py:397-464); consumed by ``transport_sweep`` as
    ``Q_aniso`` (orpheus/sn/sweep.py:101, 271-274).
    """

    values: np.ndarray
    sn_mesh: "SNMesh"
    # Same dunders as AngularFlux.
```

These two types together form the typed RHS pair `q = (q_iso, q_aniso)` consumed by `(L+C-S_foldable).solve`.

### 2.5 `Solution` — k-eigenvalue / fixed-source result container

```python
@dataclass(frozen=True, slots=True)
class Solution:
    r"""Result of an SN transport calculation.

    Carries the k-eigenvalue (1.0 for a fixed-source problem),
    both flux representations (angular AND scalar — the user's
    "...we need an angular flux and a scalar flux..." requirement
    is honoured by carrying both, derived from each other by
    :meth:`AngularFlux.to_scalar`), and the iteration history
    (residuals + dominance ratios + per-outer keff for diagnostic
    plotting).

    Point lookup
    ------------
    Forwards to ``angular_flux.at(...)`` / ``scalar_flux.at(...)``.
    The convenience method ``flux_at(position)`` looks up the
    spatial bin via ``sn_mesh.mesh.edges`` (a future extension —
    today, callers query via ``sn_mesh.mesh.centers`` and pass
    integer indices).

    Parameters
    ----------
    keff : float
        The fundamental eigenvalue (1.0 for a fixed-source result).
    angular_flux : AngularFlux
        The phase-space directional flux.
    scalar_flux : ScalarFlux
        The angle-integrated flux. Held redundantly with
        ``angular_flux`` because both shapes are useful downstream
        (plots/tally/post-processing) AND because the algebraic
        identity ``scalar_flux == angular_flux.to_scalar()`` is
        an invariant the solver promises but consumers do not have
        to re-derive on every query.
    iter_history : IterationHistory | None
        Per-outer-iteration keff trajectory + per-inner residual
        trajectory. ``None`` for synthetic test results that don't
        track iteration state.
    elapsed_seconds : float
        Wall-clock time. Useful for performance regression tests.

    Cross-references
    ----------------
    * Grand Report v3 §21.5 lines 4252-4269 — the ``Eigenpair`` +
      ``EigenSpectrum`` containers (this is the SN k-eigenvalue
      specialisation of ``Eigenpair``).
    * orpheus/sn/solver.py:101-117 — the legacy ``SNResult`` carries
      the same fields as bare arrays; this is its typed evolution.
    """

    keff: float
    angular_flux: AngularFlux
    scalar_flux: ScalarFlux
    iter_history: "IterationHistory | None" = None
    elapsed_seconds: float = 0.0

    @property
    def sn_mesh(self) -> "SNMesh":
        """Shared mesh — both flux fields agree by construction."""
        return self.angular_flux.sn_mesh

    @property
    def positions(self) -> np.ndarray:
        """Cell-centre positions for plot consumers."""
        return self.sn_mesh.mesh.centers  # type: ignore[union-attr]


@dataclass(frozen=True, slots=True)
class IterationHistory:
    r"""Outer + inner iteration trajectories for diagnostic plots.

    Cross-references
    ----------------
    * Grand Report v3 §21.2 lines 4163-4188 (power iteration with
      dominance-ratio tracking).
    * orpheus/sn/solver.py:110 — today the bare ``list[float]``.
    """
    keff: list[float]
    outer_residuals: list[float]
    inner_residuals: list[list[float]]  # per-outer list of per-inner residuals
    dominance_ratios: list[float]
```

### 2.6 Construction conveniences

```python
@dataclass(frozen=True, slots=True)
class AngularFlux:
    # ... fields and methods above ...

    @classmethod
    def zeros(cls, sn_mesh: "SNMesh", ng: int) -> "AngularFlux":
        """Zero flux on the given mesh + group count."""
        return cls(
            values=np.zeros(
                (sn_mesh.quad.N, sn_mesh.nx, sn_mesh.ny, ng)
            ),
            sn_mesh=sn_mesh,
        )

    @classmethod
    def from_scalar(
        cls, scalar: "ScalarFlux",
    ) -> "AngularFlux":
        """Isotropic embedding — see ``ScalarFlux.broadcast_to_ordinates``."""
        return scalar.broadcast_to_ordinates()
```

`ScalarFlux.zeros` is analogous.

---

## §3 Operator-algebra coupling

### 3.1 The four leaves under the typed API

Each leaf's new signature:

| Leaf | New signature | Internal flow |
|---|---|---|
| `StreamingOperator.apply` | `(psi: AngularFlux) → AngularFlux` | Per ordinate `n`, iterate `sn_mesh.dag_walk(ordinate_idx=n)`, fold over `cell_update.residual(...)` (orpheus/sn/spatial/cell_update.py:518) to compute per-cell `(L_cell ψ̄ − q)` with `q=0`. Subtract `σ_t · ψ` packed-vector style (Resolution A's algebraic identity). Pack the per-cell residuals into `(N, nx, ny, ng)` shape. Wrap as `AngularFlux`. |
| `CollisionOperator.apply` | `(psi: AngularFlux) → AngularFlux` | `values = sigma[None, :, :, :] * psi.values` — one numpy broadcast multiply. Wrap as `AngularFlux`. |
| `ScatteringOperator.apply` | `(psi: AngularFlux) → AngularFlux` | Internally project to `(nx, ny, ng)` (line 738), apply P0 + (n,2n) (line 743-744), broadcast, add Pℓ≥1 (line 747-750). Return packed as `AngularFlux`. **No shape change** — this is what `apply` already does at line 695. |
| `FissionOperator.apply` | `(psi: AngularFlux) → AngularFlux` | Project to scalar (`psi.to_scalar()`); compute `chi[g] · Σ_g' νΣ_f(g') · φ(g')` as today (line 162-165); broadcast back across ordinates. The broadcast makes the algebra close — F's domain and range are both `AngularFlux`, so `(L + C − S − F/k).apply(ψ)` distributes. |

The strategy-driven discretisation is the crux of `StreamingOperator.apply`'s rewrite. Today its body calls `transport_operator_matvec_*` (orpheus/sn/operator.py:1700-1736), which hardcodes finite-difference upwind + τ-symmetric M-M closures regardless of what `sn_mesh.cell_update` is. The new body uses the strategy:

```python
# Pseudocode — new StreamingOperator.apply
def apply(self, psi: AngularFlux) -> AngularFlux:
    sn_mesh = self.sn_mesh
    cell_update = sn_mesh.cell_update          # strategy
    out_values = np.zeros_like(psi.values)

    for n in range(sn_mesh.quad.N):
        psi_in_face = ...                      # BC-applied at first cell
        for visit in sn_mesh.dag_walk(ordinate_idx=n):
            upstream = UpstreamState(
                spatial_upstream=psi_in_face,
                angular_upstream=...,           # for curvilinear
            )
            # apply direction = residual(ψ̄=current value, source=0)
            residual = cell_update.residual(
                cell_avg=psi.values[n, visit.cell_idx, 0, :],
                visit=visit,
                total_xs=self.sigma_t[visit.cell_idx, 0, :],
                source=np.zeros(psi.ng),       # streaming only — no source
                upstream_state=upstream,
            )
            out_values[n, visit.cell_idx, 0, :] = residual
            # propagate downstream face flux via cell_update.update
            ...                                # advance psi_in_face

    # Resolution A: subtract σ_t · ψ to get pure L
    # (cell_update.residual computes the FULL (L + C) action;
    #  the subtraction isolates L).
    out_values -= self.sigma_t[None, :, :, :] * psi.values

    return AngularFlux(values=out_values, sn_mesh=sn_mesh)
```

This is the user's hard requirement made literal: *"DiamondDifference is a selectable strategy. The operator should be formed by the strategy"* — `cell_update.residual` IS the strategy. The legacy matvec primitives become unused (and retire per the aggressive-retirement directive); the discretisation choice flows through `sn_mesh.cell_update` for both `apply` and `solve` (Pattern 2 dual-view, with `apply` ↔ `update` already an established round-trip invariant at cell_update.py:537-545).

### 3.2 `F`'s domain — `AngularFlux → AngularFlux`?

The §7 brief raises this question explicitly. Today (orpheus/sn/fission.py:137) F's apply takes `ScalarFlux`-shaped `phi` and returns `ScalarFlux`-shaped output. Under the typed contract, F **must** be `AngularFlux → AngularFlux` for the operator algebra to compose: `OperatorSum.apply(ψ)` distributes as `a.apply(ψ) + b.apply(ψ)`, so every operand needs the same domain.

Mathematically, F is rank-1 in energy + rank-0 in angle. The natural action is:

```
F(ψ)(n, r, g) = χ_g(r) · Σ_g' νΣ_f,g'(r) · φ_g'(r)
              = χ_g(r) · Σ_g' νΣ_f,g'(r) · [Σ_n' w_n' ψ_n'(r, g')]
```

so F internally:
1. Projects: `phi = psi.to_scalar()`
2. Applies rank-1 op: `q_per_cell_g = chi[..., g] * Σ_g' sig_p[..., g'] * phi[..., g']` — the existing fission_rate × spectrum computation (fission.py:163-165)
3. Broadcasts back: `q_aniso_n(r, g) = q_per_cell(r, g)` for every n

The k-eigenvalue iteration becomes:

```python
q_fission: AngularFlux = F.apply(psi) / k        # F as AngularFlux → AngularFlux
psi_next: AngularFlux = A_loss_inv.apply(q_fission)
```

This is consistent: every operand in `(L + C − S).solve` and in F sees `AngularFlux`. There's no ambiguity; F's internal projection/broadcast is an implementation detail of the operator, not a contract leak.

### 3.3 `OperatorSum.apply` composes cleanly

The existing `OperatorSum.apply` at orpheus/numerics/operator.py:561-562 is:

```python
def apply(self, x: np.ndarray) -> np.ndarray:
    return self.a.apply(x) + self.b.apply(x)
```

Under the typed contract this is unchanged. The `x` parameter is now `AngularFlux`; both operand applies return `AngularFlux`; the `+` in the body dispatches to `AngularFlux.__add__` (which delegates to numpy `+` on the underlying values with the mesh-identity check). One line — no special-casing needed.

`(L + C − S_foldable − F/k).apply(ψ)` then evaluates to:
```
((L.apply(ψ) + C.apply(ψ)) + (-1.0 * S_foldable).apply(ψ)) + (-1.0/k * F).apply(ψ)
= L(ψ) + C(ψ) − S_foldable(ψ) − F(ψ)/k     # all AngularFlux, summed elementwise
```

### 3.4 `(L + C − S_foldable).solve(b: AngularFlux) → AngularFlux`

The fusion hook lives on `OperatorSum`: when it sees its operands match the pattern `L + C − S` where `S.is_foldable_into_sigma_r()` returns True (orpheus/sn/scattering.py:605), it routes `solve` through the within-group sweep rather than letting `OperatorSum.solve` reject the request (today `OperatorSum` doesn't propagate `solve` — orpheus/numerics/operator.py:548).

Sketch of the fusion hook:

```python
# orpheus/sn/within_group_solve.py  (NEW — fusion hook entry)

def fuse_within_group_solve(
    op_sum: OperatorSum, b: AngularFlux,
) -> AngularFlux:
    r"""Specialised solve for ``L + C - S_foldable`` patterns.

    Detects the algebra ``L + C - S_foldable`` (where
    ``S.is_foldable_into_sigma_r() == True``) and routes through
    the within-group sweep at ``transport_sweep`` (orpheus/sn/
    sweep.py:96).

    Pattern 6 — different fold pattern: the algebraic recognition
    lives at this hook; the heavy lift stays in ``transport_sweep``.
    """
    L, C, neg_S = _decompose_LCS(op_sum)          # raises if not the pattern
    sigma_r = build_removal_xs(C.sigma, neg_S.foldable_sigma())
    # Repackage b as (Q, Q_aniso) — the sweep's current contract.
    Q_iso, Q_aniso = _project_to_sweep_pair(b)
    angular_values, scalar_values = transport_sweep(
        Q=Q_iso, sig_t=sigma_r,
        sn_mesh=L.sn_mesh,
        psi_bc=L.sn_mesh._psi_bc,
        Q_aniso=Q_aniso,
    )
    return AngularFlux(values=angular_values, sn_mesh=L.sn_mesh)
```

Whether the hook lives on `OperatorSum.__init__` (detecting the pattern at construction) or on `OperatorSum.solve` (detecting at call time) is a substep 3+4.b.ii design decision — orthogonal to this memo's type proposal.

### 3.5 The scipy adapter — single pack/unpack at the boundary

scipy's BiCGSTAB/GMRES require `np.ndarray` of shape `(n_unknowns,)`. The typed contract pushes the pack/unpack to a single adapter site (orpheus/numerics/operator.py:1318 today):

```python
def as_scipy_linop_typed(
    op: LinearOperator, *, sn_mesh: SNMesh, ng: int,
) -> spla.LinearOperator:
    r"""Wrap a typed AngularFlux-domain LinearOperator for scipy.

    Pack/unpack ``AngularFlux ↔ (n_unknowns,)`` at the boundary so
    scipy's matvec sees a flat ndarray while the operator's
    ``apply`` sees a typed AngularFlux.
    """
    N, nx, ny = sn_mesh.quad.N, sn_mesh.nx, sn_mesh.ny
    n_unknowns = N * nx * ny * ng

    def pack(psi: AngularFlux) -> np.ndarray:
        return psi.values.ravel()

    def unpack(vec: np.ndarray) -> AngularFlux:
        return AngularFlux(
            values=vec.reshape(N, nx, ny, ng),
            sn_mesh=sn_mesh,
        )

    def matvec(vec: np.ndarray) -> np.ndarray:
        return pack(op.apply(unpack(vec)))

    return spla.LinearOperator(
        (n_unknowns, n_unknowns), matvec=matvec, dtype=float,
    )
```

This **replaces** the legacy `EquationMap`-packed layout used by `SNStreamingOperator.apply` today (orpheus/sn/operator.py:1199-1218). The new layout is the natural `(N, nx, ny, ng).ravel()` — simpler, no eq_map dispatch, no `(ix, iy)`+`.T.ravel(order='F')` gymnastics that the c21c2ef Resolution A subtraction relies on at line 1741-1744. The legacy `EquationMap` retires with the legacy matvec primitives.

---

## §4 Relationship to Issue #197

### 4.1 What Issue #197's Wave 1 said (orpheus/numerics/quantities.py — proposed)

The issue's literal proposal (§"What lands"):

```python
Flux            = NewType("Flux", np.ndarray)            # (nx, ng)
AngularFlux     = NewType("AngularFlux", np.ndarray)     # (N, nx, ng)
CrossSection    = NewType("CrossSection", np.ndarray)
ReactionRate    = NewType("ReactionRate", np.ndarray)
GroupRate       = NewType("GroupRate", np.ndarray)
Volume          = NewType("Volume", np.ndarray)
Position        = NewType("Position", np.ndarray)
Lethargy        = NewType("Lethargy", np.ndarray)
EmissionSpectrum = NewType("EmissionSpectrum", np.ndarray)
KEff            = NewType("KEff", float)
OpticalThickness = NewType("OpticalThickness", np.ndarray)
Albedo          = NewType("Albedo", np.ndarray)
ScatteringRatio = NewType("ScatteringRatio", np.ndarray)
```

### 4.2 What this memo proposes — partial close of #197

| #197 type | Proposed action | Reason |
|---|---|---|
| `Flux` (a.k.a. `ScalarFlux`) | **Promote to structural** (this memo's `ScalarFlux`) | Has rich operations (point-lookup, broadcast-to-ordinates, mesh ref, arithmetic), all of which a bare NewType cannot carry. |
| `AngularFlux` | **Promote to structural** | Same reasoning. Also: the shape in #197 is `(N, nx, ng)` (no `ny`); the current 2-D-Cartesian path uses `(N, nx, ny, ng)` (orpheus/sn/scattering.py:695). The structural type encodes the actual shape. |
| (new) `IsotropicSource` | **Add as structural** | Pattern 3 — names the SOURCE role of a shape-`(nx, ny, ng)` array, distinct from `ScalarFlux`. |
| (new) `PerOrdinateSource` | **Add as structural** | Same reasoning. |
| `CrossSection` | **Keep as bare NewType per #197** | Just a label; no domain of operations beyond numpy arithmetic. The Mixture dataclass already carries cross-sections (orpheus/data/macro_xs/mixture.py); a structural type would duplicate. |
| `ReactionRate` | **Keep as bare NewType** | Same. |
| `GroupRate` | **Keep as bare NewType** | Same. |
| `Volume`, `Position`, `Lethargy` | **Keep as bare NewType** | Same — these are scalar/array labels, no behaviour. |
| `KEff` | **Keep as bare NewType** | Scalar; nothing to add beyond the float. |
| `EmissionSpectrum`, `Albedo`, `OpticalThickness`, `ScatteringRatio` | **Keep as bare NewType** | Same. |

### 4.3 Partial-close proposal text

> Issue #197 partial close (Phase G Step 3+4.b): the angular-flux and
> scalar-flux types (`AngularFlux`, `ScalarFlux`) land as **structural
> frozen dataclasses** in `orpheus/sn/typed_fields.py` (not bare
> NewType) because they have a non-trivial domain of operations
> (mesh-aware point lookup, ordinate↔scalar conversions, mesh-identity
> arithmetic). Two new typed source fields, `IsotropicSource` and
> `PerOrdinateSource`, surface the existing `(Q, Q_aniso)` source pair
> that `transport_sweep` already consumes. The remaining Wave 1 types
> from #197 (`CrossSection`, `ReactionRate`, `GroupRate`, `Volume`,
> `Position`, `KEff`, `Lethargy`, `EmissionSpectrum`, `Albedo`,
> `OpticalThickness`, `ScatteringRatio`) stay as bare `NewType`
> definitions in `orpheus/numerics/quantities.py` per #197's plan.

### 4.4 mypy/pyright compatibility

`@dataclass(frozen=True, slots=True)` types work seamlessly with mypy strict mode. The NewType pattern from #197 also continues to work. Mixing them creates no friction. Example signature:

```python
def compute_keff(
    flux: ScalarFlux,                  # structural type
    mesh: Mesh1D,                      # plain class
    materials: dict[int, Mixture],     # plain class
) -> KEff:                             # NewType("KEff", float) from #197
    ...
```

---

## §5 How this solves the current shape mismatch

### 5.1 The shape problem today

Recap from the brief:

| Operator | Today's apply input/output shape |
|---|---|
| `SNStreamingOperator.apply` (orpheus/sn/operator.py:1293) | packed `(n_unknowns,)` |
| `StreamingOperator.apply` (orpheus/sn/operator.py:1666) | packed `(n_unknowns,)` |
| `CollisionOperator.apply` (orpheus/sn/operator.py:1840) | packed `(n_unknowns,)` |
| `ScatteringOperator.apply` (orpheus/sn/scattering.py:695) | unpacked `(N, nx, ny, ng)` |
| `FissionOperator.apply` (orpheus/sn/fission.py:137) | scalar `(nx, ny, ng)` |

The `OperatorSum.apply` distribution `a.apply(x) + b.apply(x)` (orpheus/numerics/operator.py:561) passes the SAME `x` to both operands. With three different incompatible shape contracts, `(L + C − S_foldable).apply(packed_ψ)` blows up at distribution: `S.apply(packed_ψ)` reshape-fails or gives garbage.

### 5.2 Under the typed contract — every leaf agrees on `AngularFlux`

```python
psi: AngularFlux                          # canonical typed input
L_psi: AngularFlux = L.apply(psi)         # AngularFlux → AngularFlux
C_psi: AngularFlux = C.apply(psi)         # AngularFlux → AngularFlux
S_psi: AngularFlux = S.apply(psi)         # AngularFlux → AngularFlux  (current shape is already this — no change)
F_psi: AngularFlux = F.apply(psi)         # AngularFlux → AngularFlux  (internal scalar projection)

# OperatorSum.apply(psi) = L.apply(psi) + C.apply(psi) - ...  — distributes cleanly
A_psi: AngularFlux = (L + C - S_foldable).apply(psi)
A_psi: AngularFlux = (L + C - S - F/k).apply(psi)
```

### 5.3 Solve direction — fusion hook routes through the sweep

```python
b: AngularFlux                            # RHS
A_wg: OperatorSum = L + C - S_foldable
psi: AngularFlux = A_wg.solve(b)          # fusion hook detects pattern,
                                          # routes through transport_sweep,
                                          # repackages bare-ndarray
                                          # output as AngularFlux
```

### 5.4 Round-trip identity at machine precision

Because both `L.apply` and `L.solve` consume the SAME `sn_mesh.cell_update` strategy (`.residual` for apply; `.update` for solve), they are the **single source of truth** for the discrete operator. The strategy-layer round-trip invariant at cell_update.py:537-545 promises:

```
strategy.residual(strategy.update(...), ...) ≈ 0   at atol=1e-13
```

This propagates to the operator level:

```
A.apply(A.solve(b)) ≈ b                   at machine precision
```

— the ERR-026 closure-bias divergence between the historical FD apply and the WDD solve goes away by construction. The Pattern 6 differential-fold-shape is recognised explicitly (apply = fold over `residual`; solve = fold over `update`), but the bodies are the same algebra with the cell-balance equation solved vs evaluated.

---

## §6 Migration plan for `c21c2ef` (b.i v2)

The user's directive is explicit: do not revert `c21c2ef`. The leaves stay; their `apply` bodies refactor in place. This section is the method-implementer's brief.

### 6.1 Files to edit (all `orpheus/sn/`)

| File | Change |
|---|---|
| `typed_fields.py` (NEW) | Land `AngularFlux`, `ScalarFlux`, `IsotropicSource`, `PerOrdinateSource`, `Solution`, `IterationHistory` per §2. |
| `operator.py:1551-1744` (`StreamingOperator`) | Change `apply(psi: np.ndarray) → np.ndarray` to `apply(psi: AngularFlux) → AngularFlux`. Body switches from `transport_operator_matvec_*` (line 1700-1736) to a per-ordinate fold over `sn_mesh.dag_walk` calling `sn_mesh.cell_update.residual`. The Resolution A subtraction at line 1741-1744 becomes `out_values -= self.sigma_t[None, :, :, :] * psi.values`. |
| `operator.py:1748-1879` (`CollisionOperator`) | Change `apply(psi: np.ndarray) → np.ndarray` to `apply(psi: AngularFlux) → AngularFlux`. Body becomes one numpy broadcast (`self.sigma[None, :, :, :] * psi.values`); the `_sigma_at_unknowns` packed-vector gather (line 1827) retires. `solve` becomes `q.values / sigma[None, ...]`. `apply_transpose == apply` (self-adjoint) stays. |
| `operator.py:1120-1478` (`SNStreamingOperator`) | **Retire** per Cardinal Rule 2 + the aggressive-retirement directive. The bundled `L + C` legacy is subsumed by the leaf algebra. The legacy `EquationMap` (`operator.py` lines ~80-410 of the file) retires too — its only consumers were the matvec primitives. |
| `operator.py:414+, 573+, 853+` (the 3 `transport_operator_matvec_*` primitives) | **Retire** — strategy-driven discretisation supersedes them. |
| `scattering.py:695-751` (`ScatteringOperator.apply`) | Type the signature: `apply(psi: AngularFlux) → AngularFlux`. The body stays — current shape is already `(N, nx, ny, ng)`, so only the type wrapping changes. |
| `fission.py:137-165` (`FissionOperator.apply`) | Change to `apply(psi: AngularFlux) → AngularFlux`. Body adds a `psi.to_scalar()` at entry and a `.broadcast_to_ordinates()` at exit (or builds the broadcast values directly with `chi[None, :, :, :] * rate[None, :, :, None]`). Eigenvalue division `/k` stays at the solver level (fission.py module docstring lines 39-53). |
| `sweep.py:96-300` (`transport_sweep` family) | Two options: (a) keep the bare-ndarray inner interface (sweep is implementation, not interface); (b) wrap at the entry to consume an `IsotropicSource + PerOrdinateSource` pair. Recommendation: **option (a)** — the sweep is the "hot" implementation surface and its bare-numpy interior is performance-critical. Wrap at `OperatorSum.solve` (fusion hook) and at the public `SNSolver.solve_*` entry points. |
| `solver.py:101-117` (`SNResult`), `solver.py:82-99` (`SNFixedSourceResult`) | Optional in this Phase: rename to `Solution`, change field types from bare arrays to `AngularFlux`/`ScalarFlux`. Backward compatibility: keep `SNResult` as a deprecated alias for one release. |

### 6.2 Tests — what survives, what rewrites

| Test file | What stays | What rewrites |
|---|---|---|
| `tests/sn/test_streaming_operator.py` (orpheus repo line ranges from rg above) | `TestCapabilities` (line 136), `TestConstructor` (line 178), `TestLinearity` (line 304), `TestSumCapabilities` (line 338) — pure algebra tests, no shape contract | `TestApplyShape` (line 201), `TestCompositionEquivalence` (line 228) — must rewrite to construct `AngularFlux` inputs and assert `AngularFlux` outputs. The composition-equivalence test gains power: `(L+C).apply(psi)` now equals `M(psi)` *and* the result is an `AngularFlux`. |
| `tests/sn/test_collision_operator.py` | `TestCapabilities` (line 110), `TestSigmaInterpretation` (line 294), `TestSigmaLayout` (line 339), `TestApplySolveIdentity` (line 247), `TestApplyTranspose` (line 274) — algebra | `TestApply` (line 135), `TestSolve` (line 190) — typed inputs |
| `tests/sn/test_streaming_operator_decomposition.py` | `TestResolutionADecomposition` (line 163), `TestSubtractiveDefinition` (line 209), `TestResolutionADifferentFromPriorWrong` (line 247) — the Resolution A algebraic identity, which survives because the same algebra holds at the AngularFlux level (one broadcast-subtract instead of packed-vector subtract) | The eq_map construction helpers `_eq_map_for` (line 113), `_call_matvec` (line 125) — retire alongside the matvec primitives. |
| `tests/sn/test_scattering_operator.py` | Everything — the ScatteringOperator's shape contract is already `(N, nx, ny, ng)`. | Only the entry/exit wrapping: `op.apply(psi)` becomes `op.apply(AngularFlux(values=psi, sn_mesh=...))`. |

### 6.3 Suggested PR boundaries (method-implementer guidance)

| PR | Scope | Tests touched |
|---|---|---|
| PR-1: Typed fields landing | `orpheus/sn/typed_fields.py`; foundation tests pinning arithmetic + dunders + `at()` + `to_scalar()`/`broadcast_to_ordinates()` round-trips | `tests/sn/test_typed_fields.py` (NEW) |
| PR-2: Strategy-driven `StreamingOperator` | `operator.py:1551-1744` rewrite to use `dag_walk` + `cell_update.residual`. Keep `SNStreamingOperator` for now. | `tests/sn/test_streaming_operator{,_decomposition}.py` |
| PR-3: `CollisionOperator` typed migration | `operator.py:1748-1879` typed apply/solve | `tests/sn/test_collision_operator.py` |
| PR-4: `ScatteringOperator` and `FissionOperator` typed migration | `scattering.py:695`, `fission.py:137`; minor wrapping only | `tests/sn/test_scattering_operator.py`, `tests/sn/test_fission_operator.py` |
| PR-5: `OperatorSum.solve` fusion hook | The within-group `(L+C-S_foldable).solve` routing. Lives in a new `orpheus/sn/within_group_solve.py`. | New tests pinning `(L+C-S_foldable).apply(L+C-S_foldable.solve(b)) ≈ b` |
| PR-6: `SNSolver` rewires to typed leaves + `Solution` | `solver.py` consumes typed leaves; deprecation shim for `SNResult`. Legacy matvec + `SNStreamingOperator` retired. | All SN regression tests; `SNResult` deprecation messages |

PR-1 through PR-4 are independent; PR-5 depends on 1-4; PR-6 depends on 1-5.

### 6.4 Regression contract during migration

The 11 frozen regression snapshots at `tests/sn/regression/snapshots/` (referenced in orpheus/sn/geometry.py:138 and elsewhere) must stay green at `rtol=1e-12`. Two ways to keep them green:

1. **Algebraic equivalence:** the new `L.apply(ψ).values - sigma_t * ψ` equals the old packed-form output up to numpy broadcast vs Fortran-ravel reordering. Tests need a reshape, not a tolerance bump.
2. **Same `cell_update` strategy default (`DiamondDifference`):** the strategy choice is preserved (orpheus/sn/geometry.py:142-144), so the discretised operator is bit-identical between the legacy matvec and the new strategy-driven path on DD problems.

The Bailey-Morel curvilinear-ERR-026 evidence table cited in `SNStreamingOperator`'s docstring (operator.py:1156-1186) flagged an O(h) drift between apply (FD upwind) and solve (WDD asymmetric). Under the new contract, both share `cell_update.residual` / `.update`, so this drift vanishes — but the regression snapshots were generated under the old apply, so they may need re-baselining for curvilinear with a re-pinned tolerance. This is a Phase G Step 3+4 sub-substep, not Phase G's tail.

---

## §7 Open questions / risks

### 7.1 Performance — dataclass overhead

A frozen+slotted dataclass adds ~16 bytes of Python overhead per instance and one indirection per attribute access. The arithmetic dunders create a new instance per binary op. For a typical SN solve at `N=16, nx=160, ny=1, ng=2` (the regression-scale problem), an outer iteration touches the angular flux ~few times. Cost: ~few μs per outer; negligible vs the milliseconds-per-sweep cost.

Mitigation: the dunders return `replace(self, values=...)` which creates ONE new instance + reuses `sn_mesh`. No deep copy. The hot path (the sweep itself) stays bare-numpy inside `transport_sweep` and only wraps at the public boundary.

**Risk: low.** Verify with a benchmark in PR-1.

### 7.2 Memory — SNMesh reference cycles

`AngularFlux.sn_mesh` is a strong reference. Multiple flux objects can share the same mesh — no memory bloat (the mesh is the larger object, holding the precomputed caches at `sn_mesh._geom_cache` / `sn_mesh._coll_cache`, see orpheus/sn/solver.py:271-283).

No reference cycle: the mesh doesn't hold back-pointers to fluxes. Standard Python GC handles it.

**Risk: very low.**

### 7.3 JAX-compatibility — PyTree registration

JAX requires user types to register as PyTrees for `jit`/`vmap`. The registration is one decorator + two functions (`tree_flatten`, `tree_unflatten`) per dataclass. The trick is that `sn_mesh` is **static** (not traced) and `values` is **dynamic** (traced).

```python
import jax.tree_util as jtu

@jtu.register_pytree_node_class
@dataclass(frozen=True, slots=True)
class AngularFlux:
    values: np.ndarray            # dynamic — traced
    sn_mesh: "SNMesh"             # static — not traced

    def tree_flatten(self):
        return (self.values,), self.sn_mesh

    @classmethod
    def tree_unflatten(cls, aux_data, children):
        return cls(values=children[0], sn_mesh=aux_data)
```

This is mechanical and additive. It does NOT have to land in Phase G — defer to the JAX-migration window (#197 Option D's deferred work).

**Risk: low. Mitigation: PyTree registration is single-PR additive work, deferrable.**

### 7.4 F's `AngularFlux → AngularFlux` contract — over-allocation

F internally projects to scalar (cost-free numpy reduce) then broadcasts back to `(N, nx, ny, ng)`. The broadcast allocates `(N×nx×ny×ng)` floats for what is essentially `(nx×ny×ng)` data times N identical broadcasts.

For `N=16, nx=160, ng=2` the cost is ~80 KB per F apply. The k-eigenvalue iteration applies F once per outer (a few-to-few-dozen times per solve). Total: <few MB; negligible.

**Optimisation if needed:** F's apply could return a "lazy broadcast" `AngularFlux` whose `values` is a numpy stride-tricks broadcast view. Read-only fluxes accept this; arithmetic forces a copy. Defer to a future PR if benchmarks demand.

**Risk: low.**

### 7.5 The "points" for plotting

The user's verbatim: *"...the points where the fields were calculated (which it probably needs to get from SNMesh), so that a plot can read the solution and get point,value out of it."*

`SNMesh` exposes:
- `sn_mesh.mesh.centers` — cell-centre positions (orpheus/geometry/mesh.py:304)
- `sn_mesh.mesh.edges` — cell-face positions (orpheus/geometry/mesh.py:237)
- `sn_mesh.dx` — cell widths (orpheus/sn/geometry.py:197)
- `sn_mesh.quad.mu_x`, `sn_mesh.quad.weights` — ordinates

For a Solution at a given `(x, n, g)`, the plot needs:
- 1-D: `x → sn_mesh.mesh.centers[ix]` for some `ix`
- 2-D: `(x, y) → sn_mesh.mesh.centers` (via `Mesh2D.centers_x`, `centers_y`)

The point-lookup API `Solution.scalar_flux.at(ix=..., iy=...)` returns the value; the position comes from `Solution.positions`. A higher-level `Solution.flux_at(x)` that does the bin-search (`np.searchsorted(sn_mesh.mesh.edges, x)`) is a one-liner future extension — defer to whenever a plot consumer asks for it.

**Risk: very low — the point-lookup is already supported by SNMesh's existing API.**

### 7.6 Boundary trace types

The grand report §16A names `InflowTraceSpace` / `OutflowTraceSpace` (orpheus/numerics/trace_space.py exists). The angular flux on `Γ_-` is structurally distinct from the volumetric angular flux (different DoF count). For Phase G, the boundary-trace data lives in the `psi_bc: dict` carried by the sweep (orpheus/sn/sweep.py:286-308). Typing those dict entries is a separate refactor (the dict is keyed by string face names; bare-array values). Defer; not blocking the four-operator algebra.

**Risk: low — the boundary trace problem is orthogonal to Phase G Step 3+4.b.**

### 7.7 `Solution` vs `SNResult` rename — backward compatibility

`SNResult` is consumed by `examples/`, `tests/`, and possibly downstream notebooks. A clean rename breaks API. Recommendation:

```python
# orpheus/sn/solver.py
@dataclass
class SNResult:
    """Deprecated alias for Solution. Will be removed in v0.X."""
    # ... legacy fields ...

# orpheus/sn/typed_fields.py
@dataclass(frozen=True, slots=True)
class Solution:
    # ... new structural fields ...

    @classmethod
    def from_legacy_result(cls, result: SNResult) -> "Solution":
        ...
```

One release with a `DeprecationWarning`, then drop. Mechanical migration.

**Risk: low — standard Python deprecation pattern.**

---

## Next dispatch brief (for the method-implementer)

> **Goal:** land the typed contracts proposed in this memo and rewire the four-operator algebra leaves to consume them. Branch: `refactor/sn-operator-algebra` at tip `c21c2ef`.
>
> **Order of operations:**
> 1. Land `orpheus/sn/typed_fields.py` with `AngularFlux`, `ScalarFlux`, `IsotropicSource`, `PerOrdinateSource`, `Solution`, `IterationHistory` per §2 of this memo. New foundation test file `tests/sn/test_typed_fields.py` pinning arithmetic, dunders, mesh-identity validation, `to_scalar()`/`broadcast_to_ordinates()` round-trips, and `at()` lookups.
> 2. Refactor `StreamingOperator.apply` (orpheus/sn/operator.py:1666) per §6.1 — body switches from `transport_operator_matvec_*` to a fold over `sn_mesh.dag_walk` calling `sn_mesh.cell_update.residual`. Retire the matvec primitives + the `EquationMap` packed layout. Update `tests/sn/test_streaming_operator{,_decomposition}.py` for the typed API.
> 3. Refactor `CollisionOperator.apply`/`solve` (orpheus/sn/operator.py:1840, 1853) per §6.1 — one broadcast multiply + one broadcast divide. Update `tests/sn/test_collision_operator.py`.
> 4. Type `ScatteringOperator.apply` (orpheus/sn/scattering.py:695) and `FissionOperator.apply` (orpheus/sn/fission.py:137). Internal flow unchanged; only the type-wrapping at entry/exit. Update related tests.
> 5. Land `OperatorSum.solve` fusion hook per §3.4 (`orpheus/sn/within_group_solve.py` new module). Foundation test pinning `(L+C-S_foldable).apply((L+C-S_foldable).solve(b)) ≈ b` at machine precision.
> 6. Rewire `SNSolver` (orpheus/sn/solver.py:152+) to consume the typed leaves; return `Solution` from `solve_*` methods. Keep `SNResult` as a deprecated alias for one release.
> 7. Retire `SNStreamingOperator` (orpheus/sn/operator.py:1120-1478) — the bundled `L+C` is now `OperatorSum(L_leaf, C_leaf)`.
> 8. File the partial close of Issue #197 per §4.3.
>
> **Anti-recommendations (DO NOT do these):**
> * Do NOT modify the `c21c2ef` constructors (`StreamingOperator(sn_mesh, sigma_t)`, `CollisionOperator(sn_mesh, sigma)`) — the user explicitly said the previous leaves stay.
> * Do NOT introduce per-energy-group dataclasses (block storage is the user's hard requirement).
> * Do NOT propose Pint or any runtime units library (Issue #197 rejected this).
> * Do NOT register the typed fields as JAX PyTrees in this Phase — defer to the JAX migration window.
> * Do NOT pre-emptively type `BoundaryTraceLaw` or `psi_bc` dict entries — orthogonal refactor, separate dispatch.
