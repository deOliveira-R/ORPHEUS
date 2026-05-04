---
name: Variant α 6-geometry family hindsight review
description: Hindsight elegance review of the 6-geometry × 2-orbit-space-class Variant α Green's-function family on variant_alpha_core. Identified 3 ready-now refactors and 2 wait-for-next-instance abstractions. Top frame: fiber bundle (BaseAtlas, AngularFiber, ChordOracle).
type: project
---

Reviewed 2026-05-02 on branch `feature/peierls-greens-cylinder` HEAD `8ae1d37`. 178/178 tests pass. Sphinx clean. The 6-geometry × 2-orbit-space-class family (the "2 classes" are one-surface-compact and two-surface — see Sphinx §`orbit-space-m-g-classification`) routes through `variant_alpha_core` (rank-1 + rank-2 closure primitives).

**Why:** User asked the meta-question — with hindsight, is THIS the most elegant formulation, or is there a foreign frame / shared structure / coordinate trick that would simplify the family further?

**How to apply:** When the next geometry is requested (3-region nested? toroidal?), revisit this memo BEFORE extending. The ready-now refactors should ship first; the wait-for items have tripwires.

## Ranked elegance-improvement candidates

### 1. (DO NOW) ChordOracle + BaseAtlas extraction shared with Nyström family

The b-partition logic (4-case sgn(μ)×sgn(r vs R_in)) is mathematically identical in slab-asym, hollow_sphere, and annulus. The chord formulas are 6 instances → 3 by axis class (planar / spherical / cylindrical). The Nyström family (peierls_nystrom) ALREADY validated unified-geometry — same chord primitives apply to BOTH families.

- Rule-of-three: PASSED. Three rank-2 + Nyström precedent = four data points.
- Cost: ~1 day refactor; risk low (178 tests guard).
- Refactor BEFORE shipping a 7th geometry.

### 2. (DO NOW) Unified `power_iterate` driver + single `GreensResult` dataclass

6 outer loops are structurally identical; only `apply_operator` and mesh structure vary. 6 nearly-identical Greens result dataclasses (k_eff, phi, iterations, converged, phi_g, iter_history).

- Rule-of-three: PASSED at 6 instances. Overdue.
- Cost: minimal abstraction overhead; iteration is already structurally identical.
- Phase-space dimension becomes a parameter, not a per-module type.

### 3. (WAIT for instance N+1) Full bundle/G-structure abstraction

`(BaseAtlas, AngularFiber, SymmetryGroup)`: sphere = (CompactRadial, S², SO(3)); hollow_sphere = (AnnularRadial, S², SO(3)); etc. Conserved invariants become Casimirs of momentum maps J: T*∂Ω → g*. G-equivariance gives automatic isotypic decomposition (spherical harmonics for sphere, Fourier modes for cylinder).

- Rule-of-three: NOT YET. Rank-1 has 2 instances; the third hasn't materialized.
- Cost: significant — equivariant-operator infrastructure.
- Tripwire: build when 3rd rank-1 geometry (toroidal? polar cap?) is requested.

### 4. (WAIT) Arnoldi/GMRES replacing power iteration

T = (I − S)^{-1} is a Neumann-convergent resolvent. Power iteration's convergence rate is the spectral gap of S. Arnoldi/eigs gives spectral information for free.

- Trigger present (resolvent + eigenvalue problem) but payoff is real only for thick optical depth.
- Tripwire: a use case where current power iteration is slow.

### 5. (WAIT, flag in algebra-of-record) Single SymPy V_α1 oracle parameterized by `chord(r, ω)`

Branch-1 V_α1 proofs are structurally parallel: surface fixed-point → constant ψ → operator on constant → ω₀. Could mechanize across geometries.

- Strengthens algebra-of-record (same canonical proof structure runs on every chord oracle), does NOT violate it.
- Tripwire: 7th geometry forcing Branch-1 revisit.

## Top frame match: Fiber bundle / G-structure

Trigger: phase-space pairs (sphere, hollow_sphere) and (cylinder, annulus) share the same angular fiber over the same radial base; only the base SUPPORT differs (full disk vs annular).

Reformulation: `(BaseAtlas, AngularFiber, ChordOracle)`:
- BaseAtlas: `Interval([0, R])` (compact) or `Interval([R_in, R_out])` (annular). Topology = π₀ of base.
- AngularFiber: SO(3)-orbit space modulo geometry symmetry — S² for sphere/hollow_sphere, S¹×[-1,1] for cylinder/annulus, [-1,1] for slabs.
- ChordOracle: function `chord(r, ω, base) → (length, exit_surface)` — the ONLY geometry-specific primitive.

The b-partition IS the pullback of the inner-boundary characteristic function through the chord map. Currently coded as 4-case if/else; under the bundle frame it's a single integration measure restricted to {chord crosses R_in}.

## Rejected frames (UNEXPLORED with reason)

- **Differential geometry / Christoffel**: chord algebra is linear-segment Euclidean within each geometry; no curvature term to redistribute. Would only become load-bearing if medium were inhomogeneous along chord.
- **Hilbert-Schmidt separable kernel**: phase-5 memo already promoted this for µ-resolved kernel; not load-bearing for boundary-resolvent family.
- **Category theory**: no compositional-structure trigger; geometries don't compose under current API.
- **QMC / Koksma-Hlawka**: no high-dimensional integration trigger; chord integrals are 1D-2D.

## Cross-method pollination

- ChordOracle extraction (candidate 1) directly pollinates Nyström family — highest-leverage refactor because it consolidates ACROSS solver families.
- Eigensolver abstraction (candidate 4) pollinates SN method's k-eigenvalue outer loop and diffusion-eigenvalue solver.
- Unified driver (candidate 2) is solver-family-internal — no cross-method pollination.

## Memory-discipline notes for next session

- The "rank-1 vs rank-2 closure" naming in `variant_alpha_core` is a SYMPTOM of the fiber-bundle frame not yet being adopted. Under the bundle frame, both are the same closure restricted to different connected components of the momentum-map level set.
- The Nyström family's unified-geometry validation is a precedent that the Variant α family has not yet harvested. This is the lowest-effort highest-payoff opportunity in the family.
