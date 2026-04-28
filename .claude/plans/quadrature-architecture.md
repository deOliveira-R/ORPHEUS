# Plan: Unified Quadrature Architecture for ORPHEUS

**Author**: Claude Opus 4.7, 2026-04-28
**Supersedes**: `.claude/plans/visibility-cone-substitution-rollout.md` (visibility-cone is now plumbing, not policy — it slots into this contract).
**Scope**: Replace the ~10 ad-hoc quadrature/grid builders and ~20 indexed-loop consumers with a single contract for 1-D quadrature, plus 2 geometry-aware recipes that consolidate the chord and observer-angular patterns.

## 1. Why

The vis-cone rollout (Phase 1A–1D, commits `7f7971f`, `bc2f10b`) surfaced an architectural smell that predates it. The codebase has, today:

- **Two redundant composite-GL implementations** (`peierls_geometry.composite_gl_r`, `cp/solver._composite_gauss_legendre`).
- **Two near-duplicate panel builders** for the visibility-cone substitution (`_build_visibility_cone_mu_grid_sphere`, `_build_visibility_cone_alpha_grid_cylinder`) — same algorithm, different parametrization.
- **20+ indexed-loop consumers** (every `compute_*_mode` and the four `compute_T_specular_*` / `compute_P_ss_*` callers) using `for k in range(n_quad): arr[k]` instead of `np.dot(wts, f_at_pts)`.
- **No first-class "observer-centred angular sweep" quadrature** — the kink-aware variant is in `build_volume_kernel` only; the 20 P_esc/G_bc/per-face mode primitives all reinvent a dumber version.
- **Float vs mpmath split** with the same physical integral implemented twice (`Ki_n` / `ki_n_float`, `build_volume_kernel` / `build_volume_kernel_adaptive`) without a precision dial in the contract.
- **Composition is ad-hoc**: panels concatenate by `np.concatenate` and `panel_bounds` is a parallel list passed by convention.

The full survey is in the explorer report (in-context). The bottom line: quadrature is foundational infrastructure, used everywhere, but it has no coherent abstraction — every consumer reinvents nodes/weights/integrate. Before continuing the visibility-cone rollout we have to fix the substrate.

## 2. Goals

1. **Single contract.** A `Quadrature1D` value object that exposes `(pts, wts)` as ndarrays plus an `integrate(f)` method. Every quadrature kind (plain GL, vis-cone GL, composite, mpmath-adaptive) returns the same type.
2. **First-class composition.** Panel concatenation is an operator (`q1 | q2`) on the contract, not a list-append idiom in each caller. Panel structure is preserved as a typed attribute (`panel_bounds`).
3. **First-class substitution.** Vis-cone, tanh, Gauss-Jacobi, etc. are *constructors* of `Quadrature1D` — not buried inside per-caller helpers.
4. **Pythonic consumption.** `q.integrate(f)` for a callable; `q.integrate_array(f_at_pts)` for a precomputed array. Either erases the `for k in range(n_quad)` smell. The contract should make indexed access *unnatural*, not just possible.
5. **Precision dial.** `dps` is a property of the constructor, not the consumer. Float-only (`float64`) is the degenerate `dps=53`-bit case of the same API.
6. **Co-located with chord helpers.** Lives next to `chord_half_lengths` in `_kernels.py` (or a sibling `_quadrature.py`). Never imports from `peierls_geometry.py` — only the other way.
7. **Geometry-aware recipes.** Two reusable factories sit one level above the primitive contract: `chord_quadrature(radii, R, …)` for surface-to-surface impact-parameter integrals, and `observer_angular_quadrature(r_obs, ω_low, ω_high, radii, …)` for observer-centred angular sweeps with tangent-angle subdivision. These are where the visibility-cone subdivision logic finally lives — *once*, not twice.

## 3. The contract — `Quadrature1D`

```python
# orpheus/derivations/_quadrature.py
from dataclasses import dataclass, field
from typing import Callable, Sequence
import numpy as np

@dataclass(frozen=True)
class Quadrature1D:
    """A 1-D quadrature rule on a (possibly composite) interval.

    The contract: ``pts`` and ``wts`` are 1-D float64 ndarrays of equal
    length such that ``∫_a^b f(y) dy ≈ Σ_i wts[i] f(pts[i])`` whenever
    f is in the family the rule was built for. Compositions of rules
    concatenate the arrays and the panel-bounds tuple, preserving the
    "an integral is a sum of panel integrals" semantics.
    """
    pts: np.ndarray
    wts: np.ndarray
    interval: tuple[float, float]
    # Tuple of (a_k, b_k) panel boundaries; for a single-panel rule
    # this is just (interval,).  Subdivisions push more entries.
    panel_bounds: tuple[tuple[float, float], ...] = ()

    def __post_init__(self) -> None:
        if self.pts.shape != self.wts.shape:
            raise ValueError(...)
        if self.pts.ndim != 1:
            raise ValueError(...)
        # If panel_bounds empty, fill with the single interval.
        if not self.panel_bounds:
            object.__setattr__(self, "panel_bounds", (self.interval,))

    def __len__(self) -> int:
        return self.pts.size

    def __iter__(self):
        # Iterate as (point, weight) pairs — never indexed.
        return zip(self.pts.tolist(), self.wts.tolist())

    def integrate(self, f: Callable[[np.ndarray], np.ndarray]) -> float:
        """Evaluate f at the nodes (vectorised) and return the weighted sum."""
        return float(np.dot(self.wts, f(self.pts)))

    def integrate_array(self, f_at_pts: np.ndarray) -> float:
        """Sum f_at_pts pre-evaluated at ``self.pts`` against the weights."""
        return float(np.dot(self.wts, f_at_pts))

    def __or__(self, other: "Quadrature1D") -> "Quadrature1D":
        """Concatenate panels.  ``q1 | q2`` requires q1.interval[1] == q2.interval[0]."""
        if not np.isclose(self.interval[1], other.interval[0]):
            raise ValueError(
                f"Cannot concatenate {self.interval} | {other.interval}: "
                f"intervals do not abut."
            )
        return Quadrature1D(
            pts=np.concatenate([self.pts, other.pts]),
            wts=np.concatenate([self.wts, other.wts]),
            interval=(self.interval[0], other.interval[1]),
            panel_bounds=self.panel_bounds + other.panel_bounds,
        )
```

That's the entire contract. ~30 LoC. The constraints come from `__post_init__`, the consumer ergonomics from `integrate` and `integrate_array`, and the composition from `__or__`.

## 4. Constructors

Each constructor is a free function that returns a `Quadrature1D`. Same return type means the consumer doesn't care which it called.

```python
def gauss_legendre(a: float, b: float, n: int) -> Quadrature1D:
    """Plain Gauss-Legendre on [a, b].  The degenerate case of every
    other rule below."""

def gauss_legendre_visibility_cone(
    a: float, b: float, n: int, *, singular_endpoint: str = "lower",
) -> Quadrature1D:
    """GL on [a, b] absorbing a √-endpoint vanishing.  Phase 1A primitive,
    re-expressed against the contract.  ``singular_endpoint`` is unchanged."""

def composite_gauss_legendre(
    breakpoints: Sequence[float], n_per_panel: int,
) -> Quadrature1D:
    """Plain GL on each subpanel of the breakpoint list, concatenated.
    Subsumes ``peierls_geometry.composite_gl_r`` and
    ``cp/solver._composite_gauss_legendre``."""

def gauss_laguerre(n: int, *, scale: float = 1.0) -> Quadrature1D:
    """GL on [0, ∞) with weight e^{-y/scale} folded in.  For exponentially
    decaying integrands.  (Currently used only in diagnostics; promoting
    keeps it consistent.)"""
```

All constructors accept an optional `dps: int = 53` (or equivalent) so the same call site can be promoted to mpmath precision without changing the consumer. Internally the high-precision path computes nodes via `mpmath.gauss_quadrature` and casts the returned arrays to `float64` — exactly what `gl_float` does today.

### Geometry-aware recipes (Layer 2)

These compose `gauss_legendre_visibility_cone` and `composite_gauss_legendre` into rules that match the recurring ORPHEUS patterns:

```python
def chord_quadrature(
    radii: np.ndarray, n_per_panel: int, *,
    split_first_panel: bool = True,
) -> Quadrature1D:
    """Impact-parameter quadrature on h ∈ [0, R] for surface-to-surface
    chord integrals on concentric annular geometries.  Subdivides at the
    interior shell radii, applies vis-cone-upper per panel (lower-endpoint
    smooth), and optionally splits [0, r_1] at r_1/2 to dodge the y_min=0
    degeneracy of the upper variant.  Both sphere and cylinder T_specular
    / P_ss reduce to a single chord_quadrature call after the
    µ → h or sin α → h substitution; the vis-cone subdivision lives here,
    not in two near-duplicate per-coordinate helpers."""

def observer_angular_quadrature(
    r_obs: float, ω_low: float, ω_high: float, radii: np.ndarray,
    n_per_panel: int,
) -> Quadrature1D:
    """ω-quadrature on [ω_low, ω_high] for an observer-centred ray sweep.
    Subdivides at the tangent angles ω_k = arcsin(r_k / r_obs) for shells
    visible from r_obs, plain GL per panel.  This is the kink-aware
    quadrature that today lives only in build_volume_kernel; promoting it
    to a primitive lets all 20 P_esc / G_bc / per-face mode call-sites
    inherit it for free."""
```

These are the only two recipes that capture all of P1+P2+P3 (chord) and P5 (observer angular) from the survey. P4 (slab geometric immunity) is plain `gauss_legendre`. P6 (volume kernel) tensor-products them. P7+P8 (Bickley / log-singular) keep `mpmath.quad` because they're inner kernels — the contract doesn't try to subsume mpmath-adaptive yet (see § 6 Q5 below).

## 5. Module layout

```
orpheus/derivations/
    _kernels.py              # chord_half_lengths, e_n, ki_n*, _shifted_legendre_*
    _quadrature.py           # NEW: Quadrature1D + constructors
    _quadrature_recipes.py   # NEW: chord_quadrature, observer_angular_quadrature
    peierls_geometry.py      # consumers; imports from above, exports nothing back
    cp_geometry.py           # consumers
    cp/solver.py             # consumers
    sn/quadrature.py         # SN angular ordinates — already a self-contained
                             # quadrature; align with Quadrature1D in Q4
```

The Sphinx home stays in §22 (Coordinate transformations in Nyström quadrature). The current §22.7 (visibility-cone substitution) is preserved as the *math* of the substitution; the *contract* gets a new §22.0 ("Quadrature contract") right after the principle subsection §22.1, before the catalogue of substitutions.

## 6. Migration phases

Each phase is its own commit. No phase ships without all consumers it touches passing the existing test suite. The vis-cone rollout (`visibility-cone-substitution-rollout.md`) is folded into Q2/Q3.

### Q1 — Ship the contract (1 session, ~3 hours)

1. Write `_quadrature.py` with `Quadrature1D` and the four primitive constructors (`gauss_legendre`, `gauss_legendre_visibility_cone`, `composite_gauss_legendre`, `gauss_laguerre`). The vis-cone constructor simply re-wraps the Phase 1A `gauss_legendre_visibility_cone` body in a `Quadrature1D` return.
2. Write `_quadrature_recipes.py` with `chord_quadrature` and `observer_angular_quadrature`.
3. L0 tests for each constructor in `tests/derivations/test_quadrature.py`:
   - Identity-on-polynomials (GL exact for poly degree ≤ 2n−1).
   - `Quadrature1D | Quadrature1D` concatenation correctness.
   - `integrate(callable)` and `integrate_array(values)` agreement.
   - Vis-cone constructor: re-uses the existing six tests in `test_kernels.py`, adapted to the new return type.
   - `chord_quadrature`: integrates known smooth integrands on `[0, R]` and verifies the sum partitions over `panel_bounds`.
   - `observer_angular_quadrature`: integrates a smooth angular function and verifies the tangent-angle breakpoints land at `arcsin(r_k/r_obs)`.
4. **Unwind commit `bc2f10b`** (Phase 1B). It introduced the per-coordinate helpers that this contract obsoletes. Phase 1A (`7f7971f`) stays — `gauss_legendre_visibility_cone` is the right primitive, just needs to return `Quadrature1D`.
5. Sphinx §22.0 + an updated §22.7 referencing the contract. Build clean.

**Acceptance**: contract shipped, L0 tests pass, Sphinx clean, branch back to a state functionally equivalent to `4dc03cf` plus the new contract.

### Q2 — Migrate the four chord call-sites (1 session, ~2 hours)

`compute_T_specular_sphere`, `compute_T_specular_cylinder_3d`, `compute_P_ss_sphere`, `compute_P_ss_cylinder` all become ~15-line functions that:

1. Build the impact-parameter quadrature: `q = chord_quadrature(radii, n_per_panel)`.
2. Build τ via `chord_half_lengths`: `tau = sig_t @ (2 * chord_half_lengths(radii, q.pts))`.
3. Evaluate the angular factor in h-space (`µ = √(1−h²/R²)` or `cos α = √(1−h²/R²)`).
4. Return `(prefactor) * q.integrate_array(integrand_at_h)` (or einsum for the rank-N modes).

The `T_00 = P_ss` algebraic identity becomes *trivial* (same `q`, same τ, identical integrand). `_build_visibility_cone_mu_grid_sphere` and `_build_visibility_cone_alpha_grid_cylinder` are deleted. ~250 LoC removed; ~60 LoC added.

**Acceptance**: 24/24 specular tests pass, white_hebert tests pass, T_00 = P_ss verified at 1e-16 for single AND multi-region sphere/cylinder.

### Q3 — Migrate observer-centred angular sweeps (1–2 sessions, ~4 hours)

The 20 `compute_P_esc_*` / `compute_G_bc_*` / `compute_*_mode_marshak` primitives all currently do:

```python
omega_pts, omega_wts = gl_float(n_angular, omega_low, omega_high, dps)
# … then for k in range(n_angular): scalar work …
```

Replace with:

```python
q = observer_angular_quadrature(r_obs=r_i, ω_low=omega_low, ω_high=omega_high,
                                radii=radii, n_per_panel=n_angular)
# vectorised body, integrand assembled as ndarray, result is q.integrate_array(integrand)
```

Per-mode kernel assembly (the `Σ c_m^k µ^k Ki_{2+k}(τ)` Knyazev pattern in cylinder) is hoisted into a `knyazev_mode_kernel(coefs, mu_2d, tau, max_kk)` helper that returns ndarrays. The 20 functions collapse to ~5 well-named helpers behind `compute_P_esc_mode` / `compute_G_bc_mode`.

This is the largest single LoC reduction (~500 lines removed, ~150 added).

**Acceptance**: rank-N specular and white_hebert tests pass; per-face Mark closure tests pass; the indexed-loop count drops to <5.

### Q4 — Subsume the two redundant composite-GL implementations (0.5 session)

`peierls_geometry.composite_gl_r` and `cp/solver._composite_gauss_legendre` both become thin wrappers (or delete-and-import) over `composite_gauss_legendre`. `panel_bounds` returned by the contract replaces the bespoke tuple format.

**Acceptance**: CP tests + the Phase-4.2 sphere/cylinder reference tests pass without behavior change. `_composite_gauss_legendre` is gone.

### Q5 — Wrap the mpmath-adaptive K_vol path (1 session, conditional)

`build_volume_kernel_adaptive` uses `mpmath.quad` with breakpoint hints. We expose this as a constructor:

```python
def adaptive_mpmath(
    a: float, b: float, *, breakpoints: Sequence[float] = (), dps: int = 30,
) -> Quadrature1D | "AdaptiveQuadrature1D":
```

Adaptive quadratures don't have a fixed `(pts, wts)`, so they live in a sibling type `AdaptiveQuadrature1D` that exposes only `.integrate(f)` (delegates to `mpmath.quad`). Both share a `Protocol` so consumers don't care.

**Decision point**: profile first. If `mpmath.quad` already converges spectrally (likely; that's the survey's finding), Q5 is a wrapping exercise (no convergence gain, just contract uniformity). If profiling shows it doesn't, Q5 needs more thought. Not blocking Q1–Q4.

### Q6 — Slab E₁-Nyström log-singular weights (deferred)

`peierls_slab._build_kernel_matrix` uses `mpmath.quad([p_a, x_eval, p_b])` to handle a log singularity with breakpoints. Wrap as a constructor when needed. Low priority.

## 7. Testing strategy

- **L0 contract tests** in `tests/derivations/test_quadrature.py` (NEW): polynomial-exactness for GL, composition correctness for `|`, integrate-vs-integrate_array agreement, panel_bounds invariants.
- **L0 substitution tests** preserve the six Phase 1A vis-cone tests, adapted to consume `Quadrature1D`.
- **L0 recipe tests** for `chord_quadrature` (smooth integrand on `[0, R]`, partition consistency over panel_bounds, sphere/cyl T_00 = P_ss algebraic identity at 1e-16) and `observer_angular_quadrature` (tangent-angle breakpoints, smooth integrand on partial intervals).
- **Transitive coverage** comes from the existing 24 specular tests + ~50 Peierls reference tests + CP convergence tests. We don't add transitive tests; we verify they keep passing.
- **Indexed-loop regression**: a one-line grep test (`pytest -k loops`) that fails if `for [a-z] in range(n_quad)` reappears in the migrated functions. Cheap, prevents regression.

## 8. Sphinx documentation

- New §22.0 "Quadrature contract" (~150 lines): the `Quadrature1D` value object, the composition operator, the consumer ergonomics, the precision dial.
- §22.7 (visibility-cone) keeps the math, switches its API examples from `(y_pts, y_wts) = gauss_legendre_visibility_cone(...)` to `q = gauss_legendre_visibility_cone(...)`. Same content, ~30 lines updated.
- New §22.8 "Geometry-aware recipes": `chord_quadrature` and `observer_angular_quadrature`. The §22.7 forward-references that are currently to "Phase 1B/1C/1D rollout sites" become cleaner because the rollout *is* the recipes.

## 9. Risks

| Risk | Likelihood | Mitigation |
|------|-----------|------------|
| Q2/Q3 break a tight-tolerance test | Medium | Existing tightest pin is 1e-10/1e-12 (slab T identity, white_hebert convergence); h-space rules still hit machine precision so identities are preserved. We verify rank-1 algebraic identities at 1e-14 explicitly in Q1 L0 tests before migration. |
| `Quadrature1D` immutability inconvenient | Low | Frozen dataclass is the right call (concatenation is non-destructive). If someone needs a mutable workspace, `Quadrature1D.with_replace(...)` (or just construct a new one) covers it. |
| `panel_bounds` API not used | Low | If unused after Q4, drop it — it's an attribute, not a method, low-cost to keep or remove. |
| mpmath path doesn't fit the contract (Q5) | Medium | Use a Protocol so adaptive and static quadratures share consumer-side ergonomics without sharing implementation. Don't force `(pts, wts)` on adaptive. |
| Indexed-loop "regression" test catches false positives | Low | Restrict the grep test to specific files and patterns. |

## 10. Estimated scope

| Phase | LoC added | LoC removed | Net |
|-------|-----------|-------------|-----|
| Q1    | ~250 (contract + tests + Sphinx) | 0 | +250 |
| Q2    | ~80 | ~250 | −170 |
| Q3    | ~150 | ~500 | −350 |
| Q4    | ~10 | ~120 | −110 |
| Q5    | ~80 | ~50 | +30 |
| **Total** | **~570** | **~920** | **−350** |

Plus a documentation reorg (§22.0 + §22.8) of ~250 lines. Net code reduction with vastly better factoring.

## 11. Non-goals

- 2-D / tensor-product quadratures. The volume kernel is currently the only consumer; tensor-product over two `Quadrature1D` instances inside the volume-kernel call is fine.
- Replacing mpmath as the high-precision engine. Quadrature1D *delegates* to mpmath when `dps > 53`.
- SN angular ordinates (`orpheus.sn.quadrature.GaussLegendre1D`) — that's a separate quadrature family (transport-sweep weights) with different invariants. Q4 only checks alignment, doesn't merge.
- Slab E₁-Nyström log-singular weights. Deferred to Q6.

## 12. Branch strategy

Current branch: `feature/visibility-cone-quadrature` (commits `7f7971f`, `bc2f10b`).

Action:

1. Revert `bc2f10b` (Phase 1B) — its per-coordinate helpers are obsolete under this plan.
2. Keep `7f7971f` (Phase 1A) — `gauss_legendre_visibility_cone` is the right primitive.
3. Rename branch to `feature/quadrature-architecture`.
4. Q1 lands as a new commit on this branch.
5. Q2–Q4 land as separate commits.
6. PR title: `feat(derivations): unified Quadrature1D contract and geometry-aware recipes`.

## 13. Open questions for the user

1. Should `Quadrature1D` be a frozen dataclass (current proposal) or a regular class with `__slots__`? Frozen is more functional/safer; non-frozen lets a future "adaptive workspace" mutate.
2. Should `composite_gauss_legendre` accept a `Sequence[float]` (breakpoints) or `Sequence[tuple[float, int]]` (breakpoints with per-panel `n` overrides)? The latter is more flexible for non-uniform panel quadrature.
3. Q5's `AdaptiveQuadrature1D` Protocol — necessary now or defer? My read: defer until a non-mpmath adaptive consumer appears.
4. Should we also align `orpheus.sn.quadrature.GaussLegendre1D` with `Quadrature1D` in this work, or is that a separate concern? My read: separate — SN ordinates carry transport-sweep semantics (forward/backward halves) that don't fit a generic 1-D quadrature contract.
