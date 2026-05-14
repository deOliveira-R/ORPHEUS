---
name: issue-196-phase-g-step2-5c-native-expression
description: Seven-frame attack on the Step 2.5c candidate "precomputed coefficient tensor + Blelloch scan". CONFIRMED: the candidate IS the right native expression, but the cache axis-factoring is wrong. Native immutability has THREE strata (geometry-only / cross-section-only / source-mutable), not one. Frames 1 (banded solve) REJECTED, 3 (lax.scan API) ADOPT-VOCABULARY, 4 (cyclic reduction) DEFER, 5 (operator-splitting view L⊕C) CONFIRMS the stratification, 6 (einsum / cumprod-cumsum decomposition) is what the body literally does, 7 (generating-function) DEFERRED. Frame 2 (banded triangular solve) refuted on the same grounds as the prior memo's Frame 4 (matrix-formation regression). Promotion candidate Smell #16 retired in favor of new Smell #16 "cache shape that mixes immutability strata". Cardinal action: factor cache as `SweepCoefficientCache` with three sub-dataclasses by mutation cadence.
metadata:
  type: project
---

# Phase G Step 2.5c — native expression of the precomputed cache

Branch: `refactor/sn-operator-algebra`. Date: 2026-05-14.

The Step 2.5b primitive (`ordinate_scan` + `affine_coefficients`)
is verified at rtol=1e-14 (pair-monoid associativity) and ships
2× faster than Step 2.5. The remaining bottleneck is the per-cell
`np.fromiter` work in `affine_coefficients` and the per-cell
`iter_cell_visits` generator. Step 2.5c proposes to PRECOMPUTE
the coefficients once at problem setup, then apply many sources
through pure tensor ops.

The user's question: is the candidate
"`SweepCoefficientCache` + three tensor ops in `slab_sweep`"
the right native expression?

---

## STRUCTURAL FEATURES

1. The per-cell `affine_coefficients(visits, sig_t_chain, QV_chain, psi_angle_chain)`
   has THREE input strata with DIFFERENT mutation cadences:
   - `visits` (geometry data: `abs_mu, A_inner, A_outer, A_down, dA_w, alpha_in,
     alpha_out, tau, V`) — **immutable across the entire problem lifetime**
     (set at SNMesh construction, never re-read).
   - `sig_t_chain` — **immutable within a fixed-source / k-eigenvalue outer iteration**
     (re-evaluated only on material change / depletion step / thermal feedback).
   - `QV_chain` + `psi_angle_chain` — **mutable every source iteration**
     (each Picard step re-evaluates `Q = Σs·φ + Σf·φ`).
2. The DD coefficient algebra splits along strata:
   - `denom = 2|μ|·A_down + dA·c_out + Σ_t·V` — depends on geometry + Σ_t.
   - `a = 2|μ|·A_total / denom − 1` — depends on geometry + Σ_t.
   - `b_iso = 2·source / denom` — depends on geometry + Σ_t + source.
   - `b_ang = 2·dA·c_in·ψ_a_in / denom` — depends on geometry + Σ_t + previous-iter angular state.
3. `inverse_denom = 1/denom` is computed twice per ordinate (once for `a`, once for `b`).
   This is the user's "literally no need to be computing stuff that is immutable" trigger.
4. The per-ordinate chain reordering `[chain_cell_idx]` flips arrays for inward
   sweeps. The chain-order is a **geometry-only invariant** per ordinate;
   currently recomputed every sweep.
5. The two cumulatives `cumprod(a)` and `cumsum(b/cumprod_a)` are SCALAR scans
   over `(nx,)` per `(ordinate, group)` pair. `(N, ng)` is the **batch axis**:
   independent, embarrassingly batchable.
6. The angular-state chain `ψ_a_out[i] = (ψ_avg[i] − (1−τ)·ψ_a_in[i])/τ` is
   ALSO an affine recurrence — but in the **angular index within a μ-level**,
   not the spatial index. Different chain axis, same monoid.
7. The candidate `SweepCoefficientCache.inverse_denom: (N, nx, ng)` collapses
   strata 1 and 2 into one tensor. `(N, nx)` is geometry-only;
   `(N, nx, ng)` is geometry+Σ_t. Mixing them ties the cache lifetime
   to Σ_t changes when most of the work is geometry-only.
8. The slab algebra `a = (g_s − Σ_t·V)/(g_s + Σ_t·V)` with `g_s = 2|μ|·A_total`
   is **separable in (geometry, Σ_t)** by the structure `f(g_s, x) = (g_s − x)/(g_s + x)`
   — a Möbius transform in `x = Σ_t·V` parameterised by `g_s`. Same for curvilinear
   `denom`.
9. The four return arrays (`a, b, inverse_denom-implicit, chain_perm`) have natural
   shapes `(N, nx, ng), (N, nx, ng), (N, nx, ng), (N, nx)` — but `chain_perm` is
   geometry-only and `(inverse_denom, a)` for slab is `(N, nx) × (1, 1, ng)`-broadcastable
   if `Σ_t` is energy-only (constant in space).

---

## ELEGANCE DETECTOR HITS

- **Smell #5** (the "compute twice what was computed once"): `denom` is built
  inside `affine_coefficients` every sweep call from per-cell scalars that
  haven't changed since SNMesh construction. The user named this directly:
  "no need to be computing stuff that is immutable once the problem is set up."
- **Smell #8** (nested loops with hidden algebraic structure): `for n in range(N)`
  in the sweep body wraps a per-ordinate scan. The `(N, nx, ng)` cache shape
  exposes that the ordinate axis IS just a batch axis — there is no per-ordinate
  computation that is not a numpy elementwise op.
- **Smell #13 — NEW CANDIDATE** (mixed-cadence cache shape): the proposed
  `SweepCoefficientCache.a_attenuation: (N, nx, ng)` collapses three immutability
  strata into one shape. The cache lifetime is determined by its most-mutable
  member; here geometry-only quantities ride the same lifetime as the
  Σ_t-dependent ones. Promotion candidate to skill Part C.

---

## FRAME CANDIDATES

### Frame 1 — Banded triangular solve (scipy.linalg.solve_banded / LAPACK ?gtsv)

**Trigger**: Feature 2 (`(L+C)` is bidiagonal lower-triangular per ordinate)
+ the user's question "does LAPACK give better performance".

**Reformulation**: For each `(ordinate, group)`, build the explicit
`nx×nx` bidiagonal matrix
```
[ denom_0,           0,           0,   ...  ]
[ −2|μ|·A_in_1,  denom_1,         0,   ...  ]
[      0,        −2|μ|·A_in_2, denom_2, ...  ]
```
and call `scipy.linalg.solve_banded((1, 0), AB, q)`. Or assemble the block-diagonal
`(N·nx, N·nx)` and one LAPACK call.

**Elegance payoff** (theoretical):
- *Algorithmic-advantage*: cuSolver `batched_gtsv` on GPU exists.

**REJECTED** — same grounds as the prior memo's Frame 4, sharpened:
- **Matrix-formation regression**: forming `(N·nx, N·nx)` banded storage is O(N·nx·ng)
  numpy ops to populate vs the current 3 numpy ops for `ordinate_scan`. The current
  candidate's bottleneck is NOT the scan itself (1% of time per the brief);
  it is per-cell `np.fromiter`. LAPACK does not help with `np.fromiter`.
- **Lost separation of cache strata**: the matrix encodes (geometry × Σ_t)
  jointly. The Möbius-separable structure (Feature 8) is destroyed. No
  pre-factorisation amortisation is possible because `(L+C)` changes with
  Σ_t between outer iterations — and within a single Picard fixed-point iteration
  the SAME `(L+C)` IS solved against many right-hand-sides (one per iteration).
  Step 2.5c's cache (one `inverse_denom`, many `b`) IS the amortisation. LAPACK
  pre-factor + back-substitute is the same idea with extra ceremony.
- **Bit-identity loss**: `STRSV` does forward substitution in a LAPACK-defined
  order that differs from `cumprod + cumsum` at ULP. The 11 regression snapshots
  would re-baseline for zero gain.
- **No second consumer**: GPU is not in scope; LAPACK on CPU buys nothing
  over numpy elementwise for `nx ≤ 1000`. (Confirmed by the [Tridiagonal matrix algorithm article](https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm) — Thomas algorithm IS forward substitution on a bidiagonal-form factor; the cumprod/cumsum decomposition is the elementwise unroll.)

**Frame-match score**: 5/5 (the structural identification is exact).
**Elegance score**: 1/5 (named primitive, but adds matrix formation as overhead).
**Performance**: SLOWER than candidate (matrix formation dominates).
**Recommendation**: REJECT. Keep `ordinate_scan` as the numpy scan.

**Structural attack on current**: none. The candidate already avoids matrix
formation; that is the right call.

---

### Frame 2 — Tensor network / MPO

**Trigger**: per-cell transition matrix `[[a, b], [0, 1]]` is a rank-1 MPO bond.

**REJECTED for DD/Step 2.5c** (same as prior memo's Frame 5):
- Bond dimension is 1. MPS framing degenerates to chain-scan.
- Reopens at rank-N closures (LD: bond 2; characteristic with multiple moments: higher).
  Not Step 2.5c's problem.

**Frame-match score**: 2/5 (technically correct, structurally degenerate).
**Elegance score**: 0/5 (overkill for bond-1 case).
**Performance**: equivalent to candidate.
**Recommendation**: REJECT for current scope. Already in Memory as deferral.

---

### Frame 3 — JAX-style typed associative scan (`lax.scan(f, init, xs)`)

**Trigger**: Feature 5 (ordinate axis is a pure batch axis; the per-cell
operation is a functional fold).

**Reformulation** — vocabulary, not implementation:
```python
# The math notation:  ψ_{i+1} = f(ψ_i, params_i), with params_i = (a_i, b_i)
# JAX:                 ψ = lax.scan(step, ψ_0, xs=(a, b))
# Our numpy:           ψ = ordinate_scan(a, b, ψ_0)   # identical signature
```
The candidate `slab_sweep(cache, source, bc_inflow) → ψ` IS the
`lax.scan(step, init, xs)` API, signature-by-signature. `cache` plays
the role of `xs` (precomputed per-step parameters); `bc_inflow` is `init`;
the returned `ψ` is `lax.scan`'s output trajectory.

**Elegance payoff**:
- *Expressive*: the API name `ordinate_scan(a, b, ψ_0)` reads as the math
  notation `ψ = scan(step, ψ_0, xs=(a,b))`. The user's elegance test passes
  literally. NO API CHANGE needed.
- *Structure-exposing*: confirms that the signature is the right shape;
  the cache `xs = (a, b)` (geometry-derived) is naturally separated
  from `init = ψ_0` (boundary) and from `source` (mutable).
- *Algorithmic-advantage*: zero NOW; available LATER if JAX is adopted as a
  backend (no API redesign needed).

**Frame-match score**: 5/5 (the API IS `lax.scan`).
**Elegance score**: 4/5 (validates the design vocabulary).
**Performance**: equivalent now; latent GPU path for later.
**Recommendation**: ADOPT VOCABULARY. Confirms the candidate's signature shape.

**First test**: rewrite `ordinate_scan`'s docstring to declare it as the numpy
backend for `lax.scan`-style API with a one-paragraph signature-equivalence note.
No code change. The test is "a JAX user reading our scan finds the API familiar".

**Structural attack on current**: the candidate's
`slab_sweep(cache, source, bc_inflow)` should be renamed `apply_sweep` or
`ordinate_apply(cache, source, ψ_0)` — `slab_` in the name is a geometry
leak (sphere/cylinder use the SAME signature once `cache` is built). The
right name is geometry-blind: it's a CACHE-APPLY operation.

---

### Frame 4 — Cyclic reduction / Stone's algorithm

**Trigger**: alternative parallel algorithm for bidiagonal systems.

**REJECTED for Step 2.5c**:
- nx=160 chain gives log₂(160) ≈ 7 reduction levels. Each level is
  a numpy op of size nx/2^k. Total work ~2·nx (Brent-optimal); same
  as cumprod/cumsum.
- Stone's algorithm has well-known *numerical-stability* issues for
  unstable recurrences (`|a| > 1`). Our DD recurrence has `|a| < 1`
  in well-resolved regime so this is not blocking — but cumprod is
  equally well-conditioned in the same regime.
- ZERO win on CPU sequential; non-trivial reorganisation cost.

**Frame-match score**: 3/5 (correct alt-algorithm, no payoff).
**Elegance score**: 1/5 (replaces one named primitive with another).
**Performance**: equivalent on CPU; latent for GPU (subsumed by Frame 3).
**Recommendation**: REJECT. (Brent's blocked-scan equivalence test already
covers the parallel-form theorem.)

---

### Frame 5 — Operator-splitting view (L and C separated)

**Trigger**: Feature 1 (THREE strata: geometry, Σ_t, source). The split
`(L+C)` has `L` = streaming + curvature (geometry only) and `C` = `Σ_t·V·I`
(collision; the Σ_t-dependent piece).

**Reformulation**: factor the cache as

```python
@dataclass(frozen=True, slots=True)
class GeometryCoefficients:
    """Stratum 1: geometry-only. Lifetime = SNMesh lifetime.

    Precomputed at SNMesh construction; NEVER reread once built.
    """
    g_streaming: np.ndarray   # (N, nx)   = 2|μ|·A_down + dA·c_out  (slab→2|μ|; sphere/cyl→curvature)
    g_attenuator: np.ndarray  # (N, nx)   = 2|μ|·(A_inner + A_outer)
    g_angular: np.ndarray     # (N, nx)   = dA·c_in           (slab→0; only curvilinear nonzero)
    g_volume: np.ndarray      # (N, nx)   = V[i] per cell
    chain_perm: np.ndarray    # (N, nx)   int — chain order per ordinate (sign-resolved)


@dataclass(frozen=True, slots=True)
class CollisionCache:
    """Stratum 2: geometry × Σ_t. Lifetime = constant Σ_t epoch.

    Rebuilt only on material / cross-section change (depletion, thermal
    feedback step). Stable across all source iterations within an outer.
    """
    inverse_denom: np.ndarray   # (N, nx, ng) = 1 / (g_streaming + Σ_t · g_volume)
    a: np.ndarray               # (N, nx, ng) = g_attenuator · inverse_denom − 1


def sweep_step(geom: GeometryCoefficients,
               coll: CollisionCache,
               source: np.ndarray,        # (nx, ng)
               psi_angle: np.ndarray | None,  # (N, nx, ng) or None for slab
               bc_inflow: np.ndarray) -> np.ndarray:
    """Stratum 3: source iteration. Apply the cache to a source.

    Pure tensor ops — no Python loop, no fromiter, no visit list.
    """
    # b[n, i, g] = 2·(source[i, g] + g_angular[n, i]·psi_angle[n, i, g]) · inverse_denom[n, i, g]
    ang_term = 0.0 if psi_angle is None else geom.g_angular[..., None] * psi_angle
    b = 2.0 * (source[None, :, :] + ang_term) * coll.inverse_denom
    # Vectorised scan over (nx,) for every (N, ng) pair.
    return ordinate_scan(coll.a, b, bc_inflow[..., None, :])
```

**Elegance payoff** (multi-criterion hit):
- *Structure-exposing*: the operator-split `(L + C)^{-1}` IS the cache;
  `L` lives in `GeometryCoefficients` (the dimensionful, geometry-pinned
  factors) and `C` lives in `CollisionCache` (the Σ_t-dependent
  combination). The cache reads as the math expression of the operator.
- *Expressive*: `sweep_step(geom, coll, source, ψ_a, ψ_in)` reads as
  `(L+C)^{-1} q`. One line per stratum.
- *Structurally-simpler*: collapses `affine_coefficients` from 1 method
  with 8 `np.fromiter` calls into ONE construction-time builder
  (`geom = GeometryCoefficients.from_mesh(sn_mesh)`) plus ONE
  Σ_t-epoch builder (`coll = CollisionCache.from_geom(geom, sig_t)`).
  Source-iteration body has zero geometry plumbing.
- *Algorithmic-advantage*: source iteration body becomes
  `O(N · nx · ng)` flops with ~5 numpy ops total. The per-cell
  `np.fromiter` (the actual current bottleneck per brief) is GONE
  because the per-cell scalars live in `geom`, packed at SNMesh time.

**Frame-match score**: 5/5 (this IS the operator split applied to the cache shape).
**Elegance score**: 5/5 (every criterion hits — structural, expressive, simpler, faster).
**Performance verdict**: FASTER than the candidate's flat
`SweepCoefficientCache: (N, nx, ng)`. The win is not from algorithmic
change but from PROPERLY-PLACED IMMUTABILITY: the per-cell `np.fromiter`
walks (40+ per ordinate in `affine_coefficients`) happen ONCE at SNMesh
construction, not every sweep.
**Recommendation**: **ADOPT**. This is the right native expression.

**First test**: `test_collision_cache_invariance_under_source_iteration` —
under a converged Picard loop, assert `coll.inverse_denom` and `coll.a` are
NEVER rebuilt (instrument via a counter). Pass condition: source-iteration
calls only `sweep_step`; cache rebuilds happen only at `update_sig_t`.

**Structural attack on current**: `SweepCoefficientCache.a_attenuation: (N, nx, ng)`
collapses Σ_t-dependence into the same tensor as the geometry factors.
This pretends `a` is geometry-only when it is actually geometry × Σ_t.
The cache cannot then be reused across depletion steps without full
rebuild; with the split (`GeometryCoefficients` separate), depletion
only rebuilds `CollisionCache`, half the work.

**Precedent**: matches the pattern in CP / MOC where chord-track geometry
is cached separately from optical-depth coefficients (which depend on Σ_t).
Cross-method pollination confirmed below.

---

### Frame 6 — `einsum` / closed-form cumulative decomposition

**Trigger**: Feature 5 (the cumprod-cumsum decomposition IS the literal
math); the user's elegance test ("reads like math notation").

**Reformulation**: examine what's "the math". The Blelloch closed form

`ψ_n = (∏_{i=0}^{n−1} a_i) · (ψ_0 + Σ_{k=0}^{n−1} b_k / ∏_{j=0}^{k} a_j)`

can be written compactly via the **resolvent / Green's function** form

`ψ_n = G_{n,0} · ψ_0 + Σ_{k=0}^{n−1} G_{n, k+1} · b_k`

with `G_{n, m} = ∏_{i=m}^{n−1} a_i` the discrete Green's function of
the chain. This IS what `ordinate_scan` computes: `cumprod_a` IS the
row of G against ψ_0; `cumsum(b / cumprod_a)` is the discrete convolution
of the source against the Green's function (Toeplitz-with-multiplicative
structure).

**Elegance payoff**:
- *Structure-exposing*: names the cumprod/cumsum decomposition for what
  it IS — the discrete Green's function representation. The "two scalar
  scans" are the diagonal of the resolvent matrix and the convolution
  of the source against it.
- *Structurally-simpler*: `ordinate_scan(a, b, ψ_0)` is the right
  abstraction; documenting the Green's-function reading in the docstring
  is enough.

**Frame-match score**: 5/5 (exact algebraic identity).
**Elegance score**: 3/5 (vocabulary win; no code change).
**Performance**: identical.
**Recommendation**: ADOPT VOCABULARY. Sharpen the `ordinate_scan`
docstring to note the Green's-function reading. Cross-link from the
Sphinx theory page.

**First test**: read-only; assert `ordinate_scan(a, 0·b, ψ_0)` equals
`G[:, 0] · ψ_0` (the homogeneous-source case is one Green's-function row).
Already covered by `test_ordinate_scan_zero_source`.

**Structural attack on current**: the docstring talks about
"Blelloch §1.5 first-order linear recurrence". That is the discrete-numerics
provenance. The PHYSICS provenance is "discrete Green's function of the
streaming-plus-collision chain". Two complementary readings; the
physics one is the user's elegance criterion.

---

### Frame 7 — Generating-function / shift-operator algebra

**Trigger**: Feature 1 (the recurrence is the discrete form of an ODE);
Karp-Miller-Winograd 1967 operational calculus.

**Reformulation**: define the generating function
`Ψ(z) = Σ_n ψ_n z^n`. The recurrence `ψ_{n+1} = a_n ψ_n + b_n` in
operator form is `(1 − a·E^{−1}) Ψ = b` where `E^{−1}` is the
shift-down operator. Inverting: `Ψ = (1 − a·E^{−1})^{−1} b`. This is
the formal Neumann series in the shift operator.

**REJECTED for current scope**:
- The generating-function machinery shines when `a` is *constant* (then
  the recurrence becomes a convolution with `a^n` kernel — Z-transform).
  Our `a_i` is spatially varying; the closed form degenerates to the
  Blelloch form (Frame 6). No additional insight.
- For *constant* Σ_t per region, the generating function COULD give
  closed forms for the homogeneous-region Green's function. That is
  exactly the `exp(−Σ_t · Δx / |μ|)`-style MOC ray-trace formula.
  Not a new tool — it IS the MOC path.

**Frame-match score**: 4/5 (correct generalisation, no payoff).
**Elegance score**: 1/5 (more abstraction, less code).
**Performance**: equivalent / N/A.
**Recommendation**: DEFER. Reopen if multi-region homogeneous-coefficient
fast path is implemented (then Z-transform gives closed-form
region-traversal coefficients — see cross-method pollination below).

---

## CROSS-METHOD POLLINATION

Current method: SN sweep, Step 2.5c precomputed-cache phase.

**Borrowings**:

- **From MOC (per-ray track geometry vs optical-depth coefficients)**:
  MOC long-tracking precomputes track length and segment-cell intersection
  geometry ONCE at problem setup. Optical depths `τ_seg = Σ_t · ℓ_seg`
  rebuild on Σ_t change. Source-iteration applies `exp(−τ_seg)` weights.
  Three strata: geometry / Σ_t / source — IDENTICAL to the SN split
  proposed in Frame 5. Trigger: same `(L+C)` operator with the same
  immutability cadence. **This validates Frame 5's three-stratum cache**
  via existing precedent in an adjacent solver. First test: structural
  consistency — pattern-match `MOC.TrackGeometry` vs
  `SN.GeometryCoefficients`; the API shapes should rhyme.

- **From CP (response-matrix cache)**: CP precomputes region-to-region
  first-flight probabilities once at problem setup; they're Σ_t-dependent
  but source-independent. CP source-iteration is one matrix-vector product
  per Picard step. Same stratification. The SN cache analog: source-iteration
  is one `sweep_step` per Picard step against a pre-built `(geom, coll)` pair.

- **From eigenvalue iteration**: the outer power iteration / Krylov consumes
  the SAME `(geom, coll)` across many fission-source updates. Cache lifetime
  spans the entire eigenvalue solve, not the inner SI. Three-stratum cache
  prevents redundant rebuilds.

- **From RES (resonance self-shielding) / depletion**: Σ_t changes per
  burn step / per resonance update. Geometry doesn't. `CollisionCache`
  rebuilds on these events; `GeometryCoefficients` survives. Mixed-stratum
  cache would force full rebuild on every depletion micro-step.

---

## UNEXPLORED

- **Z-transform / Wiener-Hopf** — same as prior memo; no half-space convolution
  trigger.
- **Sparse Cholesky / Krylov preconditioning** — applies at outer linear-solver
  layer, not sweep primitive.
- **Cache-oblivious algorithms (Frigo-Strumpen)** — recurrence is already
  cache-friendly (sequential access in `nx`); no benefit.
- **Symbolic / sympy-derived coefficient builder** — adds dependency for
  zero algorithmic win; the closed forms are already in `cell_balance.py`.
- **Block triangularisation** — `(L+C)` IS already block-triangular;
  exploited by the per-ordinate decomposition.
- **Padé / continued-fraction expansion of the per-cell update** — applies
  to higher-order closures (LD, EC); the DD update IS its own continued
  fraction at order 1.

---

## THE STRUCTURAL VERDICT

The candidate `SweepCoefficientCache + slab_sweep` IS the right
direction. Three sharpenings:

### 1. The cache must be factored by immutability stratum (THE main point)

The candidate's flat `SweepCoefficientCache: (N, nx, ng)` violates Smell
#13. Replace with TWO dataclasses:

- **`GeometryCoefficients`** — built once at SNMesh construction:
  `g_streaming, g_attenuator, g_angular, g_volume, chain_perm`.
  Shape `(N, nx)` (no `ng`). Geometry-only invariants.
- **`CollisionCache`** — built when Σ_t changes:
  `inverse_denom, a`. Shape `(N, nx, ng)`. Cached for the Σ_t epoch.

Source iteration consumes `(GeometryCoefficients, CollisionCache, source,
ψ_angle, ψ_in)`; the cache rebuilds at the cadence each piece actually
demands.

### 2. The physics-units / dimensional shape

| Quantity            | Shape    | Stratum                | Reads as                                                 |
| ------------------- | -------- | ---------------------- | -------------------------------------------------------- |
| `g_streaming`       | `(N, nx)`  | Geometry only          | Streaming flux at the downstream face per ψ-unit         |
| `g_attenuator`      | `(N, nx)`  | Geometry only          | Spatial-upstream multiplier `2|μ|·(A_in + A_out)`        |
| `g_angular`         | `(N, nx)`  | Geometry only          | Angular-redistribution multiplier `dA·c_in`              |
| `g_volume`          | `(N, nx)`  | Geometry only          | Cell volume per ordinate                                 |
| `chain_perm`        | `(N, nx)`  | Geometry only (int)    | Sweep-direction cell ordering                            |
| `inverse_denom`     | `(N, nx, ng)`| Geometry × Σ_t        | Resolvent diagonal — `1 / (L+C)_{ii}` per group          |
| `a`                 | `(N, nx, ng)`| Geometry × Σ_t        | Per-cell affine multiplier; `\|a\| < 1` in well-resolved regime |

Each row reads as a named physics quantity, not a "coefficient". The user's
elegance test passes: the cache fields ARE reactor-physics invariants.

### 3. The math notation

The slab sweep reads:

```
                                                geometry         × Σ_t      × source
                              ┌──────────────┐ ┌───────────┐ ┌──────────────┐
ψ_{out, n, ·, g}  =  ordinate_scan(a_{n, ·, g}, b_{n, ·, g}, ψ_{in, n, g})
                                                                ↑
                              with b = 2·(source + g_ang · ψ_a) · inverse_denom
```

For sphere/cylinder, the SAME line with `g_ang · ψ_a` non-zero. For slab,
`g_ang = 0` so `b = 2·source · inverse_denom`. ONE line, ONE expression,
geometry-blind by data.

This IS the math notation `(L+C)^{-1} q` decomposed by the discrete
Green's function: `(L+C)^{-1} = scan(a, b=q·inverse_denom)`.

The user's elegance command is satisfied: the code reads as the math.

---

## What "literally no need to be computing" gets eliminated

| Currently computed per sweep             | Stratum demands         | Action                                          |
| ---------------------------------------- | ----------------------- | ----------------------------------------------- |
| `np.fromiter(abs_mu)` × 8 fields × N     | Geometry only           | Hoist to `GeometryCoefficients` (built ONCE).   |
| `list(sn_mesh.iter_cell_visits(...))`    | Geometry only           | Replaced by `chain_perm` indexing.              |
| `chain_cell_idx = np.fromiter(...)`      | Geometry only           | Replaced by `chain_perm[n]`.                    |
| `denom = ... per ordinate`               | Geometry × Σ_t          | Hoist to `CollisionCache` (built per Σ_t epoch).|
| `inverse_denom = 1/denom`                | Geometry × Σ_t          | Cache field (single division per cell, ONCE).   |
| `a = ... per ordinate`                   | Geometry × Σ_t          | Cache field.                                    |
| `b = 2·(source + ang_contrib)·inv_denom` | Source iteration        | Computed each sweep (correctly mutable).        |
| `ordinate_scan(a, b, ψ_in)`              | Source iteration        | The scan (1% of brief's profile — keep).        |

Brief's stated "80% is per-cell fromiter + visit iter + per-cell streaming_terms"
moves to `GeometryCoefficients.from_mesh(sn_mesh)` which runs ONCE at SNMesh
construction. Source-iteration cost collapses to `b`-build + scan.

---

## Frame scores summary

| Frame                              | Match | Elegance | Perf vs candidate | Recommendation         |
| ---------------------------------- | :---: | :------: | :---------------: | ---------------------- |
| 1. Banded triangular solve         | 5     | 1        | slower            | REJECT                 |
| 2. Tensor network / MPO            | 2     | 0        | equiv             | REJECT (rank-1 degen)  |
| 3. JAX `lax.scan` API              | 5     | 4        | equiv             | ADOPT VOCABULARY       |
| 4. Cyclic reduction                | 3     | 1        | equiv             | REJECT                 |
| 5. Operator-splitting cache strata | 5     | 5        | FASTER            | **ADOPT — cardinal**   |
| 6. Green's-function / cumprod      | 5     | 3        | identical         | ADOPT VOCABULARY       |
| 7. Generating-function             | 4     | 1        | equiv             | DEFER                  |

---

## SELF-CORRECTION

None. The reclassification held. Two rejections (Frames 1, 2, 4) bound by
explicit conditions; Frame 5 hits all four elegance criteria so adopted
unambiguously. The user's "physics units" question is answered by the
table above, not by hedging.

---

## Memory updates (apply next skill revision)

1. **Promote Smell #16 candidate** — "Cache shape that mixes immutability
   strata": fires when a precomputed cache has fields with different
   mutation cadences packed into the same tensor / dataclass. Catches
   the Step 2.5c first-draft `SweepCoefficientCache: (N, nx, ng)`
   that mixed geometry-only and Σ_t-dependent fields.

2. **Sharpen A.7 row "First-order affine recurrence (Blelloch §1.5)"**
   to add the **three-stratum cache pattern**: when the per-cell affine
   coefficients themselves split by immutability cadence (geometry /
   Σ_t / source), factor the cache to match. Cross-method confirmed
   (MOC track-vs-optical-depth, CP region-pair Σ_t cache).

3. **Cross-link** this memo to `[[issue-196-phase-g-chain-dag-scan-frame-attack]]`
   (which validated the scan primitive) — Step 2.5b validated the operator;
   Step 2.5c validates the cache structure consuming it.

4. **Pollination precedent**: MOC long-tracking three-stratum cache (track
   geometry / Σ_t optical depth / source). Add to scripts/ as a validated
   precedent file (`scripts/validated_three_stratum_cache.md`) on next
   skill revision; the pattern is broader than SN.

Sources:
- [Tridiagonal matrix algorithm (Wikipedia)](https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm)
- [Sparse Direct Solvers — Techniques of High-Performance Computing](https://tbetcke.github.io/hpc_lecture_notes/sparse_direct_solvers.html)
- [scipy.linalg.solve docs](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.solve.html)
