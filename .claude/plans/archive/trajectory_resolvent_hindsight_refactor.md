# Trajectory Resolvent — Hindsight Architectural Review

**Author**: method-implementer
**Date**: 2026-05-03
**Branch**: `main` (review only — no code change in this dispatch)
**Inputs read**: every file in
`orpheus/derivations/continuous/trajectory_resolvent/` (8 modules,
4992 LOC); the prior closeouts
`peierls_layout_reorganization.md`,
`cylinder_variant_alpha_phase2_unification.md`,
`numerics-investigator/peierls_greens_phase1_closeout.md` (referenced
by manifest only); the cross-domain-attacker memos
`variant_alpha_2surface_bie_frame.md`,
`variant_alpha_family_hindsight.md` (load-bearing — frame-level
analysis already done); `.claude/scratch/folder_naming_taxonomy.md`;
`docs/theory/trajectory_resolvent.rst` (Nexus-confirmed degree 199);
`orpheus/numerics/eigenvalue.py` (existing `power_iteration` driver
+ `EigenvalueSolver` Protocol); `orpheus/derivations/common/kernels.py`
(existing `chord_half_lengths` primitive).

This document is a **review + plan deliverable**. No code changes
land in this dispatch. The user's directive: *"after regression
cases are built we will start implementing the architectural
improvements identified."* Implementation comes later, gated by
the test-architect's cross-method regression net.

The cross-domain-attacker has already produced
`variant_alpha_family_hindsight.md` covering the *frame-level*
analysis (fiber bundle / G-structure as the top match). This memo
is the **implementer-perspective companion**: same target, but
ground level — concrete code-level proposals with cost estimates,
risk register, and the regression-net plan that gates each
refactor under Cardinal Rule 1.

---

## Phase 1 — Inventory

### Files and structural roles

| File                                  | LOC  | Role                                                | Public solvers                                                                         |
| ------------------------------------- | ---: | --------------------------------------------------- | -------------------------------------------------------------------------------------- |
| `__init__.py`                         |   30 | Package docstring                                   | —                                                                                      |
| `variant_alpha_core.py`               |  436 | **Shared core**: rank-1 + rank-2 closure primitives | `compute_resolvent_T`, `compute_resolvent_T_rank2`, `apply_variant_alpha_closure`, `apply_variant_alpha_closure_rank2` |
| `greens_function.py`                  | 1321 | Sphere (rank-1, 1G + MG + multi-region + fixed-source) | `solve_greens_function_sphere{,_mg,_mr,_mr_fixed_source}`                              |
| `greens_function_cylinder.py`         |  714 | Cylinder (rank-1, 1G + MG)                          | `solve_greens_function_cylinder{,_mg}`                                                 |
| `greens_function_slab.py`             |  343 | Slab symmetric (thin facade over slab_asymmetric)   | `solve_greens_function_slab{,_mg}`                                                     |
| `greens_function_slab_asymmetric.py`  |  613 | Slab asymmetric (rank-2, 1G + MG)                   | `solve_greens_function_slab_asymmetric{,_mg}`                                          |
| `greens_function_hollow_sphere.py`    |  735 | Hollow sphere (rank-1 + rank-2 b-partition)         | `solve_greens_function_hollow_sphere{,_mg}`                                            |
| `greens_function_annulus.py`          |  800 | Annulus (rank-1 + rank-2 b-partition, 3D angular)   | `solve_greens_function_annulus{,_mg}`                                                  |

### Power-iteration loop count

| File                                | `for it in range(max_iter)` loops |
| ----------------------------------- | --------------------------------: |
| `greens_function.py` (sphere)       | 4 (1G + MG + MR + MR-fixed-src)   |
| `greens_function_slab_asymmetric.py`| 2 (1G + MG)                       |
| `greens_function_cylinder.py`       | 2 (1G + MG)                       |
| `greens_function_hollow_sphere.py`  | 2 (1G + MG)                       |
| `greens_function_annulus.py`        | 2 (1G + MG)                       |
| `greens_function_slab.py`           | 0 (delegates)                     |
| **Total**                           | **12**                            |

All twelve loops are structurally identical:

1. compute scalar flux `phi(r)` (or `phi(x)`) from `psi`
2. build per-group source profile `q_g(r)/(4π) = inv_4pi · (Σ_s_in + νΣ_f/k · F)`
3. apply per-geometry operator (the only geometry-specific call) → `psi_new`
4. compute fission rates via `fission_rate(phi)` (geometry-specific volume integral)
5. Rayleigh quotient `k_new = k_eff · F_new / F_old`
6. normalise `psi_new / F_new`
7. convergence on `|Δk|/|k| < tol`

Steps 1, 2, 4, 5, 6, 7 are **byte-identical algebra** across every file.
Only step 3 (the operator) and the volume Jacobian inside step 4
(`r²` for sphere, `r` for cylinder, `1` for slab) carry geometry.

### Result dataclasses (currently 12 distinct types)

| Type                                | Geometry-specific axes              |
| ----------------------------------- | ----------------------------------- |
| `GreensFunctionResult`              | sphere `(n_r, n_mu)`                |
| `GreensFunctionMGResult`            | sphere MG `(G, n_r, n_mu)`          |
| `GreensFunctionMRResult`            | sphere MR `(G, n_r, n_mu)` + region |
| `GreensFunctionFixedSourceResult`   | sphere MR fixed-source              |
| `CylinderGreensResult`              | cylinder `(n_r, n_mu_axial, n_phi_az)` |
| `CylinderGreensMGResult`            | cylinder MG                         |
| `SlabGreensResult`                  | slab `(n_x, n_mu)`                  |
| `SlabGreensMGResult`                | slab MG                             |
| `SlabAsymmetricGreensResult`        | slab asym                           |
| `SlabAsymmetricGreensMGResult`      | slab asym MG                        |
| `HollowSphereGreensResult`          | hollow sphere                       |
| `HollowSphereGreensMGResult`        | hollow sphere MG                    |
| `AnnulusGreensResult`               | annulus                             |
| `AnnulusGreensMGResult`             | annulus MG                          |

12 distinct dataclasses for what is structurally the same object
(`k_eff, psi, phi, mesh-metadata, iterations, converged`).
Differences are *axis labels and shapes*, not *kinds of payload*.

### Per-geometry operator-helper functions (chord/trajectory primitives)

| File                                | Helpers                                                                              |
| ----------------------------------- | ------------------------------------------------------------------------------------ |
| `greens_function.py` (sphere)       | `_apply_operator{,_with_source_profile,_mr}`, `_region_at_radius`, `_trajectory_segments`, `_chord_segments` |
| `greens_function_cylinder.py`       | `_first_leg_2d_chord`, `_impact_parameter`, `_bounce_period_2d_chord`, `_apply_operator_cylinder`, `_scalar_flux_from_psi` |
| `greens_function_slab_asymmetric.py`| `_apply_operator_slab_asymmetric`, `_scalar_flux_from_psi`                           |
| `greens_function_hollow_sphere.py`  | `_apply_operator_hollow_sphere`, `_scalar_flux_from_psi`                             |
| `greens_function_annulus.py`        | `_apply_operator_annulus`, `_scalar_flux_from_psi`                                   |

Three distinct `_scalar_flux_from_psi` implementations exist
(slab/slab-asym, cyl/annulus, hollow-sphere). They are byte-identical
between slab and hollow_sphere (`2π · ∫ψdμ`) and share the
einsum pattern `(...rm,m,p->...r)` between cyl and annulus. Sphere's
scalar-flux extraction is inlined inside the power-iteration loop
(`2π · np.sum(psi · mu_weights, axis=...)`); it is NOT factored out.

### Boundary between "core" and "per-geometry"

- **In core**: 4 functions (`compute_resolvent_T`,
  `compute_resolvent_T_rank2`, `apply_variant_alpha_closure`,
  `apply_variant_alpha_closure_rank2`). All are rank-1/rank-2
  *closure* operations on already-computed `(F, B, τ_first, τ_period)`
  scalar/array tuples.
- **Per geometry**: chord arithmetic, trajectory parameterisation,
  scalar-flux extraction, fission-rate volume Jacobian, power-
  iteration scaffold, result dataclass, multi-region piecewise
  segmentation (sphere only).

The **boundary lives at the closure**. Everything *upstream* of the
closure (computing `F` and `B`) is per-geometry; the closure is
shared. This is the Phase-2 unification's verdict, and it has
held — but Phase-3 added 4 more geometries on top of the same
boundary, so the per-geometry overhead has now become 6× the
core's contribution. Time to push the boundary outward.

---

## Phase 2 — shared-concept gap analysis

### Class A — already unified in `variant_alpha_core.py`

- **Rank-1 boundary-to-boundary scattering resolvent**
  `T(α, τ_period) = 1/(1 − α·exp(−τ_period))` — `compute_resolvent_T`.
- **Rank-2 monodromy resolvent**
  `T = (I − S)^{-1}` with `S = antidiag(α_L·e^{−τ}, α_R·e^{−τ})` —
  `compute_resolvent_T_rank2`.
- **Rank-1 closure** `ψ_new = F + e^{−τ_first}·αB·T` —
  `apply_variant_alpha_closure`.
- **Rank-2 closure** with surface selection (left/right or
  inner/outer) — `apply_variant_alpha_closure_rank2`.
- **`alpha_per_period` 1-bounce-vs-2-bounce slab handling** —
  encoded as kwarg in the closure (Phase 3A pattern).

These four primitives compose to handle every geometry. Six
geometries × two topologies all consume the same four functions.
Ninety-four foundation tests pin them
(`tests/derivations/test_peierls_variant_alpha_core.py`).

### Class B — could/should be unified, but isn't (the prize)

Ranked by leverage × low-risk. Each candidate has a Branch-1 SymPy
anchor (the **algebra-of-record** — refactors don't touch it),
a Branch-2 Python sketch, a cost estimate, a risk assessment,
and a regression-net specification.

#### B1. **`ChordOracle` ABC + per-geometry implementations** (HIGH leverage, MEDIUM cost)

**Currently scattered across:**
- `greens_function.py` lines 184-188 (sphere first-leg + bounce-period)
- `greens_function.py` lines 717-794 (sphere multi-region segmentation)
- `greens_function_cylinder.py` lines 163-209 (`_first_leg_2d_chord`,
  `_impact_parameter`, `_bounce_period_2d_chord`)
- `greens_function_cylinder.py` lines 285-336 (3D-lifted chord
  arithmetic with `s_in_plane` axial projection)
- `greens_function_slab_asymmetric.py` lines 240-268 (`L_first`,
  `L_transit`, sign-of-μ branching for surface selection)
- `greens_function_hollow_sphere.py` lines 220-298 (b-partition,
  `L_step` shell-traversal, sign-of-μ for inner/outer first arrival)
- `greens_function_annulus.py` (3D-lifted hollow-sphere mirror)

The **mathematical concept is one**: given a starting interior
state `(r, ω)` and a domain Ω, a `ChordOracle` answers *"how long
is the backward characteristic to the first surface, what is the
bounce-period chord at the conserved impact parameter, and which
surface(s) does the trajectory touch?"* The 6 geometries differ
only in evaluating these primitives.

The b-partition logic in hollow_sphere/annulus is *structurally
identical* — `if b > R_in: rank-1 outer-only; else: rank-2
through-ray`. The slab-asymmetric `if mu > 0: surface=left, else
surface=right` branch is the same pattern in disguise (it's the
slab's b-partition reduced to a sign test because slab's "impact
parameter" degenerates).

**Unified abstraction sketch** (Branch-2 only — algebra-of-record
in `origins/specular/*.py` SymPy modules untouched):

```python
# orpheus/derivations/continuous/trajectory_resolvent/chord_oracle.py

from typing import Protocol, Literal
import numpy as np
from dataclasses import dataclass

@dataclass(frozen=True)
class ChordResult:
    """All chord primitives needed by the Variant α closure for one
    interior phase-space point. Geometry-agnostic."""
    L_first: float                # backward chord to first surface
    L_period: float | None        # bounce-period chord (rank-1) or None
    L_transit: float | None       # single-transit (rank-2) or None
    rank: Literal[1, 2]           # closure rank for this trajectory
    surface_first: str            # which surface absorbs the first arrival
                                  # ('left'/'right'/'inner'/'outer'/'closed')
    impact_parameter: float       # b — conserved across the bounce period
    # Trajectory parameterisations (callable for vectorised quadrature):
    r_first: callable             # s ∈ [0, L_first] → position-along-chord
    r_period: callable | None     # s ∈ [0, L_period] → position
    r_transit_LR: callable | None # rank-2 only
    r_transit_RL: callable | None # rank-2 only

class ChordOracle(Protocol):
    """Per-geometry chord/trajectory primitive."""
    def __call__(self, r: float, mu: float, *, phi_az: float = 0.0
                ) -> ChordResult: ...

# Concrete implementations (one per geometry):
def sphere_oracle(R: float, sigma_t: float) -> ChordOracle: ...
def cylinder_oracle(R: float, sigma_t: float) -> ChordOracle: ...
def slab_oracle(L: float, sigma_t: float) -> ChordOracle: ...
def slab_asymmetric_oracle(L, sigma_t) -> ChordOracle: ...  # = slab + surface routing
def hollow_sphere_oracle(R_in, R_out, sigma_t) -> ChordOracle: ...
def annulus_oracle(R_in, R_out, sigma_t) -> ChordOracle: ...
```

The operator becomes:

```python
def apply_operator(
    oracle: ChordOracle,
    source_interp: callable,
    grid: PhaseSpaceGrid,
    closure: ClosureKind,  # rank-1 or rank-2
    sigma_t: float,
    alpha: tuple[float, ...],   # (α,) for rank-1, (α_L, α_R) for rank-2
    *,
    n_traj_quad: int,
) -> np.ndarray:
    """Geometry-agnostic Variant α operator."""
    # 1. For each (r, ω), call oracle → ChordResult
    # 2. Quadrature on r_first → F (always)
    # 3. If rank-1: quadrature on r_period → B; closure rank-1
    # 4. If rank-2: quadratures on r_transit_LR, r_transit_RL → B_LR, B_RL;
    #    closure rank-2
    # 5. ψ_new[i, q] = closure(F, B_*, τ_*, alpha_*)
```

Multi-region sphere/cylinder is the next composition step:
`MultiRegionChordOracle(base_oracle, radii, sigma_t_per_region)`
wraps a base oracle and decomposes its `r_first`/`r_period`
parameterisations into segments. The current
`_trajectory_segments`/`_chord_segments` becomes a method on the
multi-region wrapper.

**Cost**: M (medium) — ~2-3 days.
- Day 1: oracle Protocol + sphere + slab implementations,
  port `solve_greens_function_sphere` and `solve_greens_function_slab`.
  Foundation test: each oracle's `ChordResult` matches the inlined
  primitives at ≥10 phase-space points.
- Day 2: cylinder + hollow_sphere + annulus + slab_asymmetric
  oracles. The cylinder's 3D `s_in_plane` axial projection is the
  trickiest abstraction — needs a ChordResult variant for "2D + 3D
  attenuation factor" (or absorb the projection into `r_first`/
  `r_period` callables, which is cleaner).
- Day 3: multi-region sphere oracle composition + bit-equality
  tests against current code.

**Risk**: MEDIUM. The 3D-lift in cylinder/annulus introduces a
parametrisation-axis distinction (2D arclength `s_2D` vs 3D
arclength `s_3D` connected by `ds_3D = ds_2D / sqrt(1-μ_axial²)`).
The abstraction must encode this without forcing slab/sphere/
hollow_sphere to carry a degenerate `s_in_plane = 1` factor (a
classic premature-abstraction smell).

**Failure-mode-induced bugs** (`numerical-bug-signatures`):
- **Variable swap (#2)** — confusing `L_2D` and `L_3D` would
  produce a bug invisible at α=1 (closed-shape eigenmode is
  μ-axial-symmetric so the projection cancels) but visible at
  intermediate α with anisotropic flux.
- **Sign flip (#1)** — first-arrival surface selection (`mu > 0`
  → inner vs outer in hollow_sphere) is exactly the kind of
  mechanically-explainable AI failure mode that survives 1G
  homogeneous closed-shape testing.

**Regression net** (Cardinal Rule 1 anchor):
1. **Branch-1 SymPy is the algebra-of-record** — `origins/specular/*.py`
   modules verify the algebraic identity that the chord primitives
   are *supposed* to satisfy. The refactor cannot change Branch-1.
   Every `derive_*()` foundation test must pass byte-equal.
2. **Bit-equality preservation** — ChordResult outputs at fixed
   `(r, μ, φ_az)` tuples must match the pre-refactor inlined
   formulas to IEEE-754 exact (not just `np.allclose`). 12-point
   stratified test grid per geometry (interior, near-surface,
   grazing, multi-region interior, multi-region interface, etc.).
3. **Cross-geometry V_α1 reproduction** — every geometry's V_α1
   identity (`solve_*` at α=1 gives k_inf) must reproduce to the
   same accuracy floor as today (1e-15 sphere/cylinder, 1e-9
   hollow_sphere, etc.).
4. **PS-1982 cross-check unaffected** — the only structurally-
   independent L1 reference (sphere vacuum). Must hold at ≤ 1e-4
   across 4 configurations.
5. **Garcia 2021 multi-region cross-check unaffected** — the
   load-bearing MR test for Issue #132.
6. **Cross-method reference comparison** — once test-architect
   ships the cross-method regression net, the chord-oracle refactor
   gains additional structurally-independent comparison against
   peierls_nystrom and singular_eigenfunction at the same
   geometries. The refactor MUST land *after* that net is in place.

#### B2. **Unified `power_iterate` driver** (HIGH leverage, LOW cost)

**Currently scattered across:** 12 power-iteration loops (catalogued
above). Each is 30-60 lines, byte-identical except for the
`apply_operator` call signature, the `fission_rate` Jacobian, and
the `_scalar_flux_from_psi` reduction.

**Existing infrastructure**: `orpheus/numerics/eigenvalue.py` already
defines an `EigenvalueSolver` Protocol + a `power_iteration` driver.
trajectory_resolvent does NOT use it. The driver was designed for
the SN/CP/diffusion solvers and pre-dates the Variant α work.

**Unified abstraction sketch** (consume the existing protocol with
a Variant-α-tailored adapter, OR write a new lightweight driver
inside trajectory_resolvent that exactly fits Variant α's needs —
recommendation: the latter, because the existing driver assumes
solver-internal `solve_fixed_source` semantics that don't map
cleanly onto "apply operator once" of Variant α):

```python
# orpheus/derivations/continuous/trajectory_resolvent/power_iteration.py

@dataclass(frozen=True)
class PowerIterationResult:
    k_eff: float
    psi: np.ndarray        # converged angular flux (any shape — caller
                           # owns the axis labels)
    phi: np.ndarray        # converged scalar flux (axis-0 = position)
    iterations: int
    converged: bool
    k_history: list[float] # iter → k_eff
    residual_history: list[float]  # iter → |Δk|/|k|

def power_iterate_variant_alpha(
    *,
    apply_operator: callable[[np.ndarray, float], np.ndarray],
        # (psi, k_eff) → psi_new — caller closes over geometry, XS, alpha
    extract_scalar_flux: callable[[np.ndarray], np.ndarray],
        # (psi) → phi(r) — geometry's angular reduction
    fission_rate: callable[[np.ndarray], float],
        # (phi) → ∫νΣ_f φ dV — geometry's volume Jacobian
    initial_psi: np.ndarray,
    initial_k: float,
    *,
    max_iter: int = 200,
    tol: float = 1e-10,
) -> PowerIterationResult:
    """Geometry-agnostic Variant α power iteration."""
    psi = initial_psi.copy()
    k_eff = initial_k
    k_history = []
    residual_history = []

    for it in range(max_iter):
        psi_new = apply_operator(psi, k_eff)

        phi_old = extract_scalar_flux(psi)
        phi_new = extract_scalar_flux(psi_new)

        F_old = fission_rate(phi_old)
        F_new = fission_rate(phi_new)
        if F_old < 1e-30:
            raise RuntimeError(f"Fission rate vanished at iter {it+1}")

        k_new = k_eff * F_new / F_old
        psi = psi_new / F_new

        rel_dk = abs(k_new - k_eff) / max(abs(k_eff), 1e-30)
        k_history.append(k_new)
        residual_history.append(rel_dk)
        k_eff = k_new

        if rel_dk < tol:
            return PowerIterationResult(k_eff, psi, phi_new / F_new,
                                        it+1, True, k_history,
                                        residual_history)

    return PowerIterationResult(k_eff, psi, phi_new / F_new, max_iter,
                                False, k_history, residual_history)
```

Each `solve_greens_function_*` becomes 30→8 lines:

```python
def solve_greens_function_sphere(R, sigma_t, sigma_s, nu_sigma_f, *,
                                 alpha=1.0, ...):
    grid = SphereGrid(n_r, n_mu, R)
    initial_psi = np.ones((grid.n_r, grid.n_mu)) if initial_psi is None else initial_psi
    return power_iterate_variant_alpha(
        apply_operator=lambda psi, k:
            apply_sphere_operator(psi, grid, sigma_t, sigma_s, nu_sigma_f, k, alpha),
        extract_scalar_flux=lambda psi: grid.scalar_flux(psi),
        fission_rate=lambda phi: nu_sigma_f * grid.volume_integral(phi),
        initial_psi=initial_psi,
        initial_k=nu_sigma_f / (sigma_t - sigma_s),
        max_iter=max_iter, tol=tol,
    )
```

The `~k_history` and `~residual_history` fields are NEW — they did
not exist in any current result dataclass. They emerge from the
unification because every iteration already computes them; only
the storage was geometry-specific noise that erased them.

**Cost**: S (small) — ~half a day.
- Write driver + 6 foundation tests.
- Port one geometry (sphere 1G is the canonical) and verify
  bit-equality.
- Port the other 11 loops (mechanical).

**Risk**: LOW. The driver is pure Python flow control; the
geometric content lives in the closure callables. Bit-equality is
guaranteed by IEEE-754 reproducibility (same operations in same
order). The ONLY footgun is the `extract_scalar_flux` and
`fission_rate` closures must capture the right grid axes — a 1-line
test per geometry catches this.

**Failure-mode-induced bugs**: minimal. The driver does no
geometry-specific arithmetic; it composes callables. The most
plausible error is a closure-capture bug (Python's late binding in
`for`-loop closures) — a `pytest.parametrize` over geometries with
explicit loops in the test file catches this in seconds.

**Regression net**:
1. Bit-equality at every test of every geometry — same as B1.
2. New foundation tests for the driver:
   - `test_power_iterate_converges_at_α=1` (every geometry; reduces
     to k_inf).
   - `test_power_iterate_returns_history` (k_history strictly
     monotone after iter 2 in homogeneous closed shapes).
   - `test_power_iterate_zero_fission_raises` (no silent NaN).
3. Iteration-count parity must match within ±1 (one extra check
   on history vs convergence detection ordering).

#### B3. **`GreensResult` standardised dataclass** (MEDIUM leverage, LOW cost)

**Currently scattered across:** 12 dataclasses (catalogued in Phase 1).
The `psi`, `phi`, `iterations`, `converged`, `k_eff` fields are universal.
Geometry differs only in mesh fields (`r_nodes` vs `x_nodes`; whether
`phi_az_nodes` exists; whether `region_at_node` exists).

**Unified abstraction sketch:**

```python
@dataclass(frozen=True)
class GreensResult:
    """Universal Variant α power-iteration result.

    Geometry/group/region structure is encoded in the `mesh` field
    (a frozen dataclass produced by the geometry's grid factory) and
    in the `psi`/`phi` array shapes. The user reads `result.psi.shape`
    or `result.mesh.kind` to dispatch.
    """
    k_eff: float
    psi: np.ndarray           # angular flux on the geometry's phase-space
    phi: np.ndarray           # scalar flux on the spatial mesh
    mesh: PhaseSpaceMesh      # tagged-union over {Sphere, Cylinder, Slab,
                              # SlabAsymmetric, HollowSphere, Annulus}
                              # × {SingleRegion, MultiRegion}
    n_groups: int             # 1 for 1G, G for MG; psi axis-0 is group
                              # iff n_groups > 1
    iterations: int
    converged: bool
    k_history: list[float]    # NEW (from B2)
    residual_history: list[float]  # NEW (from B2)
```

`PhaseSpaceMesh` is a tagged union — frozen dataclass per geometry-
topology pair, carrying the axis nodes and any geometry metadata
(`region_at_node` for MR, `R_in`/`R_out` for hollow, `phi_az_nodes`
for cylinder/annulus).

**Cost**: S (small) — ~4 hours.
- Define `PhaseSpaceMesh` hierarchy + `GreensResult`.
- Replace 12 returns. Existing call sites use `result.k_eff` /
  `result.psi` / `result.phi` — those keep working unchanged.
  Mesh-axis access changes from `result.r_nodes` to
  `result.mesh.r_nodes`; this requires touching tests but is
  mechanical (45 test files have ~80 such accesses).

**Risk**: LOW with one CAVEAT — the rename is API-breaking. If
external code (downstream notebooks, ad-hoc scripts not in the test
suite) reads `result.r_nodes` directly, those break. Mitigation:
keep the old attribute names as `@property` delegators with a
DeprecationWarning for one release cycle. After that, drop them.

**Failure-mode-induced bugs**: extremely low. The dataclass holds
data, doesn't compute. The most plausible error is a wrong axis
ordering (e.g., MG `psi_g` was always `(G, n_r, n_mu)` but a code
path could have reversed to `(n_r, n_mu, G)` in some callsite).
A `psi.shape == (n_groups, *mesh.shape)` invariant test catches
this in seconds.

**Regression net**: same as B2. Bit-equality at every test —
result *contents* are identical, only the *type signature* changes.

#### B4. **`extract_scalar_flux` is one operator parametrised by axes** (LOW leverage, LOW cost)

Currently 5 implementations (3 distinct, 2 byte-equal copies):

| Geometry          | Formula                                          |
| ----------------- | ------------------------------------------------ |
| sphere            | `2π · np.sum(psi · mu_weights, axis=-1)` (inlined) |
| slab_asymmetric   | `2π · np.einsum('...xm,m->...x', psi, mu_weights)` |
| hollow_sphere     | identical to slab_asymmetric                      |
| cylinder          | `np.einsum('...rmp,m,p->...r', psi, mu_weights, phi_az_weights)` |
| annulus           | identical to cylinder                             |

**Unified**:

```python
def scalar_flux_from_angular(
    psi: np.ndarray,
    angular_weights: tuple[np.ndarray, ...],
    azimuthal_factor: float = 2 * np.pi,  # 1 if azimuth is in `psi`'s axes
) -> np.ndarray:
    """Reduce angular flux to scalar by integrating each angular axis.

    angular_weights: one ndarray per angular axis, in the order they
    appear in psi's axes.
    azimuthal_factor: 2π for slab/sphere (μ only), 1.0 for cyl/annulus
    (μ_axial × φ_az already covers the full 4π solid angle).
    """
    # Reduce trailing-most angular axes one at a time
    out = psi
    for w in reversed(angular_weights):
        out = np.einsum('...i,i->...', out, w)
    return azimuthal_factor * out
```

**Cost**: trivially small. Roll into B2 if both land together; otherwise
~1 hour.

**Risk**: LOW. Single-line replacement. The `azimuthal_factor=1.0`
for cyl/annulus is a load-bearing convention — must be tested
explicitly.

**Failure-mode-induced bugs**: Failure mode #4 — quadrature-
dependent constant hardcoded (Signature 4 in
`numerical-bug-signatures`). The `2π` vs `1.0` factor is precisely
the kind of plausible-wrong constant that survives one
quadrature-family but breaks another. A streaming-equilibrium test
(`Q` uniform, `Σ_t` uniform, `α=0` open-domain → `φ ≠ Q/Σ_t`
because of leakage; check `4π · ∫ψdμd(...)` invariant per ordinate
rather than via the unified extractor) catches it.

**Regression net**: bit-equality + the streaming-equilibrium probe
above.

#### B5. **`fission_rate(phi)` volume-integral abstraction** (LOW-MEDIUM leverage, LOW cost)

Currently scattered as inner-defined `def fission_rate(phi):` at
12 sites. The Jacobian carries the geometry:
- sphere/hollow_sphere: `4π · ∫ φ(r) r² dr`
- cylinder/annulus: `2π · ∫ φ(r) r dr` per unit length
- slab/slab_asym: `∫ φ(x) dx`

**Unified** as a method on `PhaseSpaceMesh`:

```python
class PhaseSpaceMesh:
    def volume_integral(self, integrand_on_spatial_grid: np.ndarray) -> float:
        """Integrate an array shaped like `phi` over the geometry."""
        ...
```

Each `mesh.volume_integral` is a one-liner in its concrete subclass.

**Cost**: trivially small (S). Roll into B3 (`PhaseSpaceMesh` lives
there).

**Risk**: LOW. The fission-rate formulas are short and well-tested;
the danger is *forgetting* to multiply by `νΣ_f` in some callsite
(currently it's done manually by each `fission_rate` closure).
Solution: make the abstraction `mesh.volume_integral(phi)` return
the integral, and require callers to multiply by `νΣ_f` at the
callsite — the current pattern. This keeps the unification minimal
and transparent.

**Failure-mode-induced bugs**: Failure mode #3 (missing factor) —
the `4π` vs `2π` vs `1` Jacobian volume factor is exactly the kind
of plausible-wrong constant. The current tests catch it because the
heterogeneous k_eff cross-checks (Garcia 2021 sphere MR, slab/cyl
benchmarks) pin the absolute value. The unified abstraction is
trivially correct as long as the per-geometry concrete subclasses
are tested individually.

**Regression net**: every existing k_eff test serves as the
regression net. Bit-equality at every test.

#### B6. **`region_at_radius` + segmentation can become a `MultiRegionChordOracle` decorator** (MEDIUM leverage, MEDIUM cost)

Currently sphere-only at lines 717-794. Cylinder MR (Issue #129) and
slab MR (when shipped) will need the same logic. The pattern is
*decorator over a base oracle*: given a single-region oracle, wrap
it to decompose its `r_first`/`r_period` parameterisations into
segments at region boundaries, returning multi-segment quadrature
rules.

This is a Class B candidate that **lands as a corollary of B1**
(ChordOracle extraction). Once `ChordOracle` is a Protocol, the MR
sphere code becomes `MultiRegionChordOracle(SphereOracle(R, ...),
radii, sigma_t_per_region)`.

**Cost**: rolled into B1 (free with the abstraction).

**Risk**: LOW given B1 is sound. The MR segmentation logic is
already tested at sphere; the abstraction does not change the
segmentation algebra, only its callsite.

**Regression net**: Garcia 2021 multi-region sphere benchmark
(load-bearing). The current 21-test MR + Garcia 2021 suite
guards this.

### Class C — looks unifiable, actually shouldn't be (premature abstraction)

#### C1. **Phase-space dimensionality should NOT be parametrised**

Sphere is `(n_r, n_mu)`. Cylinder is `(n_r, n_mu_axial, n_phi_az)`.
Hollow sphere is `(n_r, n_mu)`. Annulus is `(n_r, n_mu_axial, n_phi_az)`.

It is tempting to write `psi[*radial_axes, *angular_axes]` with
parametrised number of angular axes. **Resist.** The angular structure
is fundamentally different: sphere reduces by full SO(3) symmetry to
a single `μ`; cylinder retains the SO(2)×R subgroup so needs
`(μ_axial, φ_az)`. A unified type with `n_angular_axes` parameter
hides a load-bearing physical distinction (the residual symmetry
group), and a future Christoffel-style symmetry refactor (cf the
cross-domain-attacker's "wait-for-instance-N+1" candidate) would
need to re-introduce it.

**Verdict**: keep the per-geometry phase-space shape explicit. The
`PhaseSpaceMesh` tagged union (B3) holds them as distinct types.

#### C2. **rank-1 vs rank-2 closure should NOT be merged into a "rank-N" abstraction**

Tempting to write a generic `compute_resolvent_T_rankN(alphas,
chord_lengths)` that handles 1 ≤ N ≤ ... arbitrary surface counts
via a generalised matrix `(I − S)^{-1}`. **Resist for now.**

Reasons:
1. Rank-3 has no current test case. Building the abstraction now
   means inventing the test against an analytical reference that
   doesn't exist.
2. The numerical conditioning of `(I − S)^{-1}` differs by rank.
   The rank-1 division-by-`(1 − α·exp(−τ))` is a 1-line analytic
   formula. Rank-2 is a 2×2 closed form. Rank-N requires
   `np.linalg.solve` and the conditioning depends on the
   eigenvalue spectrum of `S`. Treating these uniformly would
   slow rank-1 (which is the most-called) for a hypothetical
   rank-N case that doesn't exist.
3. The cross-domain-attacker's `variant_alpha_2surface_bie_frame.md`
   memo found rank-N has *non-monotone* numerical behaviour as a
   function of N — the abstraction would need to encode that
   non-monotonicity, and currently we don't know what the right
   shape is.

**Verdict**: keep rank-1 and rank-2 as separate, fully-explicit
functions. A rank-N abstraction lands AFTER a rank-3 use case
materialises and provides an analytical reference.

#### C3. **Sphere multi-region segmentation should NOT be unified with hollow-sphere/annulus b-partition**

Both cases segment a chord — but for fundamentally different
reasons:
- Sphere MR: chord crosses **interior region boundaries** (smooth
  σ_t change; chord is *one* trajectory).
- Hollow sphere through-ray: chord crosses an **inner-cavity
  surface** (vacuum → reflection at boundary; the chord's
  trajectory class changes — the cavity transit is "no source,
  no attenuation" while the shell segments are "with source and
  σ_t").

These are *not* the same problem. MR segments a smooth-σ_t
trajectory. b-partition handles a *trajectory class change* at the
inner surface. Conflating them would mean the b-partition logic is
buried inside a "segmentation" abstraction whose primary use-case
(MR sphere) doesn't need the trajectory-class change.

**Verdict**: B1 (`ChordOracle`) handles this correctly — sphere MR
is a *decorator over* SphereOracle (B6); hollow_sphere is a
*distinct oracle* with its own b-partition logic. Same Protocol,
different concrete implementations. The Protocol is unified; the
mathematical content is not forced together.

---

## Phase 3 — Elegant-architecture gap analysis

Beyond shared-concept unification, here are architectural patterns
NOT YET implemented that would simplify the family.

### A1. **Strategy pattern for boundary closure** (RECOMMEND DEFER)

Current: each geometry hardcodes "specular Variant α" closure. Other
BCs (vacuum, white, Marshak diffuse, isotropic-incident) are not yet
parametrised at the closure level.

Sketch:

```python
class BoundaryClosure(Protocol):
    def __call__(self, F, B_state, tau_state, *params) -> ndarray: ...

class VariantAlphaSpecular(BoundaryClosure): ...
class VacuumClosure(BoundaryClosure): ...        # ψ_surf = 0; degenerate
class WhiteClosure(BoundaryClosure): ...         # diffuse re-emission
class MarshakClosure(BoundaryClosure): ...
```

**Defer**: Variant α already absorbs vacuum (`α = 0`) and specular
(`α = 1`) and partial-reflection (`α ∈ (0, 1)`) into a single
closure. White and Marshak BCs DO require a different closure
structure (they integrate over outgoing μ rather than acting
trajectory-by-trajectory), but they're not currently in scope for
trajectory_resolvent — peierls_nystrom handles them in a different
discretisation. Premature to unify until we have a use case.

**Tripwire**: when a 2nd closure kind is added to trajectory_resolvent
(e.g., diffuse white BC alongside specular), introduce the Strategy.
Today: do not.

### A2. **Composable mesh layers: hollow = base + topology operator** (POWERFUL — ROLL INTO B1)

Hollow sphere = sphere + inner-vacuum exclusion. Annulus = cylinder +
hole. Slab asymmetric = slab + per-wall reflectivity.

This pattern is ALREADY captured by ChordOracle composition (B1):

- `SphereOracle(R)` is the base.
- `HollowSphereOracle(R_in, R_out) = AnnularBaseDecorator(SphereOracle(R_out))`
  could be the formal way to encode it.
- `MultiRegionChordOracle` is a separate decorator (B6).

**Verdict**: roll into B1. The decorator pattern naturally falls out
of the Protocol abstraction. Don't introduce a separate "topology
operator" abstraction — let composition do the work.

### A3. **Symmetry decomposition (group-theoretic)** (MAJOR — DEFER until 7th geometry)

Sphere has SO(3) symmetry. Cylinder has SO(2) × R. Slab has Z₂
(reflective) or trivial (vacuum). Hollow sphere has SO(3) on the
shell. Annulus has SO(2) × R on the shell.

A symmetry-aware solver would decompose the angular flux into
isotypic components (spherical harmonics for SO(3), Fourier modes
for SO(2)). For closed shapes the dominant eigenmode is in a
SINGLE isotypic component (the trivial representation for the
fundamental); the rest can be discarded or used for accelerated
convergence.

**Cost**: XL (extra-large). Group-equivariant operator
infrastructure is a multi-week investment. The cross-domain-attacker
correctly classified this as a "wait-for-instance-N+1" candidate.

**Tripwire**: when a 7th geometry forces the third rank-1 closed
shape (toroid? polar cap?), revisit this. Until then, the explicit
quadrature on `(r, μ, φ_az)` is honest, fast, and verifiable.

### A4. **Multi-group as tensor product** (MEDIUM leverage, low-medium cost — RECOMMEND CONDITIONALLY)

Currently MG is wired into each geometry's `solve_*_mg` as
"per-group operator + scattering matrix coupling". The pattern is:

```
for group g:
    psi_new[g] = apply_per_group_operator(psi[g], source_profile_g[g])
```

The `source_profile_g` carries the cross-group scattering and
fission. The per-group operator carries the *spatial-angular*
content. These factor cleanly:

```
MultiGroupOperator(per_group_op: PerGroupOperator,
                   coupling: GroupCouplingMatrix)
```

where `GroupCouplingMatrix` carries `Σ_s[g_from, g_to]` and
`χ ⊗ νΣ_f` (the fission outer product).

**Cost**: M (medium) — ~1 day. The factoring exists implicitly;
making it explicit is a tagged-union refactor.

**Risk**: LOW. The factoring is mathematically clean. The danger
is performance (one extra dispatch per group), but the per-group
spatial-angular operator dominates.

**Recommendation**: do this AFTER B2 (power_iterate driver). Once
the driver is geometry-agnostic, the MG ↔ 1G distinction also
becomes a callable parametrisation. Land it as the natural
generalisation.

**Failure-mode**: this is exactly the place where Failure mode #2
(scattering matrix transpose convention drift, ERR-002) bit the
SN solver. The MG abstraction MUST carry an explicit convention
test (`Σ_s[g_from, g_to]` shape, `q_g = Σ_g' Σ_s[g',g] φ_g'`)
and a 2G-asymmetric test that matches the ERR-002 fingerprint
test pattern in `tests/sn/test_quadrature.py`.

### A5. **Result-type discipline + naming consistency** (LOW leverage, but high-value style)

Today: `solve_greens_function_sphere`, `_sphere_mg`, `_sphere_mr`,
`_sphere_mr_fixed_source`, `_cylinder`, `_cylinder_mg`, ...

Proposed (after B2 + B3):

```python
def solve(geometry: Geometry, *, group_count: int = 1,
          n_regions: int = 1, fixed_source: bool = False,
          alpha=...) -> GreensResult: ...
```

A single facade. Each `geometry` is a tagged union value carrying
its own parameters (`Geometry.sphere(R=...)`, `Geometry.cylinder(R=...)`,
`Geometry.hollow_sphere(R_in=..., R_out=...)`, ...).

**Cost**: S (small) once B1/B2/B3 are in place.

**Risk**: API-breaking. Mitigation: keep the current
`solve_greens_function_sphere` etc. as one-line aliases for one
release cycle.

**Recommendation**: defer to a final API-cleanup pass after B1+B2+B3
have landed and bedded in. Don't rush this; the call-site sprawl
isn't expensive while we're refactoring.

---

## Phase 4 — Cross-pollination with adjacent methods

### P1. **`ChordOracle` belongs in `orpheus.derivations.common`**, not in trajectory_resolvent

The Nyström family (`peierls_nystrom`) ALREADY has unified-geometry
chord primitives — see the cross-domain-attacker memo
(`variant_alpha_family_hindsight.md` candidate 1). The chord oracle
is a method-agnostic geometric primitive.

**Recommendation**: when B1 lands, place `ChordOracle` in
`orpheus/derivations/common/chord_oracle.py`, not inside
`trajectory_resolvent/`. Both Nyström and Variant α consume it.
The MoC family (when it joins as a continuous reference) consumes
it. CP consumes it (per `common/kernels.py` `chord_half_lengths`,
which is already a pre-form of this).

The existing `chord_half_lengths` in `common/kernels.py` should be
absorbed into the new module (it's a 2D-cylinder-on-disk specialisation
of one ChordOracle method). This pollinates **across the entire
continuous-references package**, not just within trajectory_resolvent.

**Tripwire**: do this only after B1 is validated inside
trajectory_resolvent. Promoting an unstable abstraction to `common/`
creates an unnecessarily wide blast radius.

### P2. **`power_iterate` belongs in `orpheus.numerics`**, not in trajectory_resolvent

`orpheus/numerics/eigenvalue.py` already has a `power_iteration`
driver + `EigenvalueSolver` Protocol (lines 80-139). It is currently
unused by trajectory_resolvent.

**Recommendation**: when B2 lands, EITHER consume the existing
`numerics.power_iteration` (preferred — converges the project on
one driver) OR write a Variant-α-tailored driver and upgrade the
existing one to share infrastructure. The existing driver assumes
solver-internal `solve_fixed_source` semantics that don't quite
fit Variant α (which is "apply operator once, normalise"); a
small adapter or a generic-enough refactor of the protocol
resolves this.

**Tripwire**: do this only after B2 is validated inside
trajectory_resolvent. Then evaluate whether the SN/CP/diffusion
solvers' usage of `numerics.power_iteration` survives the unification
unchanged.

### P3. **Resolvent `T = (I − S)^{-1}` is a structural asset across pillars**

The fn_method's `peierls_atkinson_nystrom.py` already has
`power_iterate_dominant_eigenmode` that converges to the dominant
eigenvalue of a Peierls integral operator. The relationship between
that and the trajectory-resolvent multi-bounce closure is structural:
both are dominant-eigenmode extraction of integral operators with
exponential kernels.

**Cross-pollination opportunity**: a unified
`orpheus.numerics.dominant_eigenmode` driver that
- consumes either an `(N, N)` matrix (fn_method case) or a callable
  operator (trajectory_resolvent case),
- offers Arnoldi/GMRES as alternative eigensolvers (the
  cross-domain-attacker's "wait" candidate 4),
- returns the full `PowerIterationResult` with history.

**Cost**: M (medium) — ~1 day, after B2 lands in trajectory_resolvent.
**Tripwire**: do this when an Arnoldi-vs-power test case materialises
(thick-optical-depth multi-region where power iteration is slow).

### P4. **`chord_half_lengths` in `common/kernels.py`** has been the proto-ChordOracle

```python
def chord_half_lengths(radii: np.ndarray, y_pts: np.ndarray) -> np.ndarray:
    """Half-chord lengths sqrt(R² − y²) for impact parameter y."""
    return np.sqrt(np.maximum(radii[:, None]**2 - y_pts[None, :]**2, 0.0))
```

This is the "circular cross-section, conserved impact parameter
y" specialisation of one ChordOracle method. The trajectory_resolvent
sphere's `_apply_operator_with_source_profile` line 188
(`L_p = 2.0 * R * mu_surf`) and cylinder's `_bounce_period_2d_chord`
line 200-209 are both sphere/disk-specialisations of the same
formula `L = 2·sqrt(R² − b²)`.

**Recommendation**: when B1 lands and `ChordOracle` is moved to
`common/`, absorb `chord_half_lengths` into the new module as the
disk-cross-section instance. Keep the old name as an alias for one
release cycle.

---

## Implementation order (the actionable plan)

The user's directive: *"after regression cases are built we will
start implementing the architectural improvements identified."*
This means the test-architect's cross-method regression net lands
FIRST. Each refactor below is gated on the regression net it needs.

### Phase R0 — Wait for the regression net (test-architect dispatch in flight)

The test-architect is building the cross-method regression framework
in parallel with this review. **Do not begin any implementation
below until that net is green at HEAD.** The refactors below are
shipped only when their L1 cross-checks (peierls_nystrom comparison,
singular_eigenfunction comparison, fn_method comparison at
geometries where both methods cover the case) are alive in CI.

### Phase R1 — Land B2 first (`power_iterate` driver) — LOW RISK, HIGH LEVERAGE

1. Write `orpheus/derivations/continuous/trajectory_resolvent/power_iteration.py`
   with `power_iterate_variant_alpha` + `PowerIterationResult`.
2. Write 6 foundation tests at
   `tests/derivations/test_trajectory_resolvent_power_iterate.py`.
3. Port sphere 1G as proof-of-concept (one commit).
4. Run full suite — acceptance: 178/178 + new 6 = 184/184. Bit-equality
   at every test.
5. Port the remaining 11 power-iteration loops (one commit per
   geometry pair).
6. Acceptance after each commit: full suite passes, k_eff bit-equal,
   iteration count within ±1.

**Why first**: lowest risk. No mathematical content moves; only
control flow. The driver is small, the tests are mechanical, and
bit-equality is provably preserved by IEEE-754 reproducibility.
Builds the team's confidence before the harder B1 refactor.

### Phase R2 — Land B3 (`GreensResult` + `PhaseSpaceMesh`) — LOW RISK

1. Define `PhaseSpaceMesh` hierarchy + `GreensResult` (one commit).
2. Replace 12 returns + 12 dataclasses (one commit per geometry).
3. Add `@property` shims for old attribute names (one commit).
4. Acceptance: bit-equality at every test, no `result.r_nodes`-style
   regression.

**Why second**: another low-risk control-flow refactor. Now the
result types are uniform, which makes B1 much easier.

### Phase R3 — Land B1 (`ChordOracle` Protocol + concrete implementations) — MEDIUM RISK, HIGHEST LEVERAGE

1. Define `ChordOracle` Protocol + `ChordResult` dataclass at
   `orpheus/derivations/continuous/trajectory_resolvent/chord_oracle.py`.
2. Implement `SphereOracle` + 12-point bit-equality test at
   `tests/derivations/test_chord_oracle_sphere.py`. **Do NOT move any
   code** in this commit — just the new oracle, tested for parity.
3. Port `_apply_operator_with_source_profile` (sphere) to consume
   `SphereOracle`. Bit-equal verification.
4. Repeat for slab → cylinder → hollow_sphere → annulus →
   slab_asymmetric.
5. Implement `MultiRegionChordOracle` decorator. Port sphere MR.
6. Acceptance after each step: full suite + cross-method regression
   net + bit-equality at 12-point test grid per geometry.

**Why third**: highest leverage but most risk. The cross-method
regression net (P1-P4 cross-pollination references) is now alive and
serves as the structurally-independent guard against any silent
algebra drift introduced by the Protocol abstraction. The
chord-oracle is the load-bearing geometric primitive — getting it
right closes the family.

### Phase R4 — Promote `ChordOracle` to `common/` (P1)

After R3 has bedded for at least one merge cycle and is stable,
move `ChordOracle` to `orpheus/derivations/common/chord_oracle.py`.
Absorb `chord_half_lengths`. Update peierls_nystrom to consume the
shared oracle (this is a separate dispatch — let the
cross-domain-attacker or method-implementer plan it as
peierls_nystrom-side work).

### Phase R5 — Optional polish: B4, B5, A4, A5 — LOW PRIORITY

`extract_scalar_flux` consolidation (B4), `fission_rate`
abstraction via `mesh.volume_integral` (B5), MG tensor-product
factoring (A4), unified `solve(geometry, ...)` facade (A5). These
all fall out as natural cleanups after R1-R3 are done. Land them
opportunistically when touching the affected files.

---

## Risk register — Cardinal Rule 1 anchors

### Regression coverage requirement

**Every refactor commit MUST run, at minimum:**

1. Full trajectory_resolvent test suite (all 178+ tests).
2. The cross-method regression net (test-architect deliverable —
   not yet shipped; gates start of R1-R3).
3. Bit-equality probe: `git diff` of stored k_eff, ψ-norm, φ-norm
   values at a fixed random seed across the suite. Any non-zero
   diff blocks the commit.

### Bit-equal preservation discipline

The Phase-2 closeout (`cylinder_variant_alpha_phase2_unification.md`)
preserved bit-equality across the sphere/cylinder unification by
keeping the inlined formulas in algebraically-canonical order and
mirroring it in the shared closure. **The same discipline applies
to all of B1-B6 + R1-R3.** Concretely:

- The `apply_operator` closure inside `power_iterate_variant_alpha`
  must call the same arithmetic operations in the same order as the
  current inlined loops.
- The `ChordOracle` callables must produce the same trajectory
  parameter values (in IEEE-754) as the current inlined formulas.
- Any deviation from bit-equality must be flagged and either
  justified by a documented accuracy improvement (with new tests)
  or rolled back.

### Branch-1 vs Branch-2 separation (algebra-of-record)

All refactors are **Branch-2** (Python production code) only. The
SymPy modules in `origins/specular/*.py` are the canonical algebra-
of-record and **MUST NOT be touched** by any of B1-B6 or P1-P4.
Every `derive_*()` foundation test must pass byte-equal before and
after every refactor commit.

If a refactor reveals a Branch-1 SymPy issue (e.g., the current
algebra-of-record is incomplete for a corner case the abstraction
exposes), STOP, log the finding, and dispatch back to the user.
Do NOT silently extend Branch-1 to accommodate Branch-2's
abstraction.

### Failure-mode-induced risks (numerical-bug-signatures cross-reference)

| Refactor | Most-plausible failure mode               | Catching test                                                                           |
| -------- | ----------------------------------------- | --------------------------------------------------------------------------------------- |
| B1       | Variable swap (`L_2D` vs `L_3D`)          | Per-ordinate bit-equality at 12-point grid + cylinder MG asymmetric-XS test            |
| B1       | Sign flip (first-arrival surface)         | Hollow sphere R_in → 0 ↔ solid sphere bit-equality + slab μ→−μ symmetry                |
| B1       | Missing factor (`s_in_plane` projection)  | Cylinder grazing-ray stability (already exists in current suite)                       |
| B2       | Closure-capture (Python late-binding)     | `pytest.parametrize` over geometries with explicit fixture                              |
| B3       | Wrong axis ordering (psi shape)           | `psi.shape == (n_groups, *mesh.shape)` invariant test                                  |
| B4       | Quadrature-dependent constant (`2π` vs 1) | Streaming-equilibrium probe per geometry (matches numerical-bug-signatures Signature 4) |
| B5       | Missing `νΣ_f` factor                     | Garcia 2021 absolute k_eff (already exists)                                            |
| A4       | Scattering matrix transpose (ERR-002)     | 2G-asymmetric XS heterogeneous k_eff (already exists in MG suite)                       |

### Reference contamination guard

The cross-method regression net (test-architect deliverable) must
include at least one **structurally-independent** L1 reference per
geometry. Verbatim from the prior closeouts:

- **Sphere**: PS-1982 vacuum closed-form (existing); F_N method
  bare-critical (`fn_method_sphere_fn_proper`); Sood Ua-1-0-SP
  registry (1e-5).
- **Cylinder**: Sood Ua-1-0-CY F_N benchmark (8.5e-6, existing);
  `singular_eigenfunction/cylinder` Mitsis-WM 1973 reference
  (≤ 3e-7, existing).
- **Slab**: `fn_method/slab` F_N collocation (existing); Sood
  Ua-1-0-SL.
- **Hollow sphere / annulus**: R_in → 0 limit recovers solid
  sphere/cylinder (existing); independently-implemented MoC or
  ray-tracing as the structurally-independent comparison
  (NOT YET SHIPPED — flag for test-architect).

If a refactor regresses any of these L1 references, it is reverted.
A bit-equal numerical agreement with another trajectory_resolvent
solver is L4 (cross-implementation), NOT correctness evidence.

### Cross-method comparison gate

**The refactors must NOT begin until the cross-method regression
framework (test-architect deliverable) is alive.** Specifically,
each geometry must have at least one peer-method comparison test
that exercises a structurally-independent reference at the same
configuration. Without that gate, a bug introduced by the refactor
that compensates an existing bug would be invisible.

---

## Summary table — the 12 refactor candidates

| #  | Name                                       | Class | Leverage | Cost | Risk    | Lands after | Foundation tests + L1 net protect against ... |
| -- | ------------------------------------------ | ----- | -------- | ---- | ------- | ----------- | --------------------------------------------- |
| B1 | `ChordOracle` Protocol                     | B     | HIGH     | M    | MEDIUM  | R0+B2+B3    | Variable swap, sign flip, missing factor      |
| B2 | `power_iterate_variant_alpha` driver       | B     | HIGH     | S    | LOW     | R0          | Closure capture                               |
| B3 | `GreensResult` + `PhaseSpaceMesh`          | B     | MED      | S    | LOW     | R0+B2       | Wrong axis ordering                           |
| B4 | `extract_scalar_flux` parametric           | B     | LOW      | S    | LOW     | B2 or B3    | Quadrature-constant (Signature 4)             |
| B5 | `mesh.volume_integral` for fission rate    | B     | LOW-MED  | S    | LOW     | B3          | Missing factor                                |
| B6 | `MultiRegionChordOracle` decorator         | B     | MED      | M    | LOW (corollary of B1) | R3 | MR segmentation regression               |
| C1 | n_angular_axes parametric                  | C     | (avoid)  | —    | (high)  | NEVER       | —                                             |
| C2 | rank-N closure abstraction                 | C     | (avoid)  | —    | (high)  | NEVER       | —                                             |
| C3 | sphere MR + b-partition unified            | C     | (avoid)  | —    | (med)   | NEVER       | —                                             |
| A1 | Boundary-closure Strategy                  | A     | LOW      | M    | LOW     | DEFER       | (no use case yet)                             |
| A2 | Composable mesh layers                     | A     | (rolls into B1) | — | — | B1            | —                                             |
| A3 | Symmetry decomposition (group-theoretic)   | A     | HIGH (eventually) | XL | HIGH | DEFER (7th geometry tripwire) | —                                |
| A4 | MG as tensor product                       | A     | MED      | M    | LOW     | B2          | Scattering convention drift (ERR-002)         |
| A5 | Unified `solve(geometry, ...)` facade      | A     | LOW      | S    | LOW (API-break) | After B1+B2+B3 | (style, not correctness)              |
| P1 | Promote `ChordOracle` to `common/`         | P     | HIGH (cross-method) | S | LOW (after B1 stable) | After R3 stable | Cross-family alignment             |
| P2 | Consume `numerics.power_iteration`         | P     | MED      | S    | LOW     | After B2 stable | Single-driver consolidation             |
| P3 | Arnoldi/GMRES eigensolver                  | P     | LOW (until thick-τ use case) | M | MED | DEFER | (no use case yet)                          |
| P4 | Absorb `chord_half_lengths` into ChordOracle | P  | LOW      | S    | LOW     | R4          | Naming/cross-method                           |

**Top-3 immediate-action items** (ordered by ratio of leverage to risk):

1. **B2** — `power_iterate_variant_alpha` driver. Eliminates 12
   structurally-identical loops. Bit-equality preserved by IEEE-754.
2. **B3** — `GreensResult` + `PhaseSpaceMesh`. Collapses 12
   dataclasses to 1+6. Mechanical refactor.
3. **B1** — `ChordOracle` Protocol. Highest mathematical leverage —
   makes the family extensible (7th, 8th geometry are now trivial)
   and cross-pollinates with peierls_nystrom + MoC + CP.

These three together are the "Variant α 2.0 architecture" the
cross-domain-attacker memo flagged. Land them in order R1 → R2 → R3
after the cross-method regression net is alive.

---

## Pointers

- **Existing infrastructure**:
  - `orpheus/numerics/eigenvalue.py` (`EigenvalueSolver` Protocol +
    `power_iteration` driver — currently unused by trajectory_resolvent)
  - `orpheus/derivations/common/kernels.py:52`
    (`chord_half_lengths` — proto-ChordOracle)
  - `orpheus/derivations/common/eigenvalue.py`
    (`kinf_homogeneous`, `kinf_and_spectrum_homogeneous` — used as
    initial-k estimator in MG/MR solvers)
- **Cross-domain-attacker memos**:
  - `.claude/agent-memory/cross-domain-attacker/variant_alpha_family_hindsight.md`
    — frame-level analysis (fiber bundle / G-structure)
  - `.claude/agent-memory/cross-domain-attacker/variant_alpha_2surface_bie_frame.md`
    — rank-2 BIE frame validation (Phase 3B foundation)
- **Prior phase closeouts** (worked-example evidence):
  - `cylinder_variant_alpha_phase1.md` (sphere → cylinder unification)
  - `cylinder_variant_alpha_phase2_unification.md` (variant_alpha_core
    extraction; the bit-equality discipline this plan inherits)
  - `slab_variant_alpha_phase3a.md` (slab; `alpha_per_period` kwarg
    pattern for 2-bounce-per-period geometries)
  - `slab_asymmetric_variant_alpha_phase3b.md` (rank-2 closure;
    surface routing; ERR-034 / ERR-035 catches)
  - `err035_phase3a_delegation_fix.md` (Phase-3A becoming a thin
    facade over Phase-3B — the slab-symmetric-as-special-case
    pattern that will recur in B1's hollow_sphere/annulus)
  - `hollow_sphere_variant_alpha_phase3c1.md` and
    `annulus_variant_alpha_phase3c2.md` (rank-2 + b-partition;
    byte-equal-shared `variant_alpha_core` validation)
- **Theory page**: `docs/theory/trajectory_resolvent.rst` (Nexus
  degree 199 — load-bearing for the family)
- **Folder taxonomy**: `.claude/scratch/folder_naming_taxonomy.md`
  (the `peierls_greens_function` → `trajectory_resolvent` rename
  rationale)

---

## Final notes

This plan is the **implementer-perspective companion** to the
cross-domain-attacker's `variant_alpha_family_hindsight.md`. The
two memos converge on the same three top-priority refactors
(B1-B2-B3 ↔ candidates 1-2-3 of the frame analysis). The
cross-domain-attacker's memo names the *frame* (fiber bundle / G-
structure with BaseAtlas + AngularFiber + ChordOracle); this memo
names the *concrete code-level moves* (Protocol + Mesh + driver)
and gates them on a regression net.

The two memos should be read together. Their convergence is
evidence that the architectural diagnosis is solid.

The user's stated discipline — *"build cylinder first, test it,
then unify and see if the tests still hold"* — applies here
recursively: ship the cross-method regression net first; then
unify B1-B2-B3 only after that net is alive and the tests still
hold. Bit-equality is the floor, not the ceiling: any refactor
that improves accuracy is welcome, but only with documented
explicit justification.
