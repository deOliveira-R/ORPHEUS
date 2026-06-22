---
name: issue-196-phase-g-chain-dag-scan-frame-attack
description: Six-frame structural attack on the chain_dag_scan + affine_coefficients primitive (Phase G Step 2.5b/2.6). CONFIRMED: associative-scan / 2x2 matrix monoid is the canonical frame; chain_dag_scan IS the affine-recurrence-via-scan identity restricted to sequential scan. CONFIRMED: per-cell update is the monoid element (a, b). REFUTED: triangular-solve, MPS, full JAX integration. RECOMMENDED: rename to ordinate_scan with explicit monoid op; keep numpy cumprod path for sequential; carve a sequential→associative seam for 2D.
metadata:
  type: project
---

# Issue #196 Phase G — chain_dag_scan frame attack

Branch: `refactor/sn-operator-algebra`. Date: 2026-05-14.

The candidate primitive `chain_dag_scan(a, b, ψ₀) = cumprod_a · (ψ₀ + cumsum(b/cumprod_a))` is the closed-form solution of the affine recurrence `ψ[i+1] = a[i]·ψ[i] + b[i]`. The question: is this the right primitive, or is there a structurally-deeper one that subsumes it AND extends to 2D wavefronts / GPU?

---

## STRUCTURAL FEATURES

1. Per-cell update is **scalar affine in chain state**: `ψ_out = a·ψ_in + b` (both slab and curvilinear, after factoring the M-M angular state into `b`).
2. The sweep is **forward substitution on (L+C)** under DAG ordering induced by direction sign (already documented in sweep.py:15-22).
3. The chain-DAG (1D) is a **totally ordered sequence** of cell visits per ordinate.
4. The 2D wavefront DAG has **partial order**: anti-diagonal levels are mutually independent given previous-level outputs; cells within a level have no chain dependence on each other.
5. The recurrence operator is **associative in pair form**: `(α₁,β₁) ⊕ (α₂,β₂) = (α₁·α₂, α₂·β₁+β₂)` — this IS the 2×2 lower-triangular matrix multiplication `[[a,0],[b,1]]` (Blelloch 1990 §1.5).
6. Curvilinear adds a **per-cell exogenous input**: the M-M angular state `ψ_angle[i]` is read from a buffer populated by previous-ordinate sweeps within the same μ-level. From the *current ordinate's* perspective `ψ_angle[i]` is known data; the recurrence stays scalar-affine.
7. The 2D update is **two-channel affine**: `ψ_avg = (Q + sx·ψ_x_in + sy·ψ_y_in)/(Σt+sx+sy)` — two upstream chains converge at one cell. Outputs are then scattered to two downstream chains. This is NOT scalar affine; it's a tree-fan-in pattern over a 2D grid.
8. Slab's pre-Step-2.5 `_sweep_1d_cumprod` was the special case `b[i]=2·q/denom`, `a[i]=(2|μ|/Δ−Σt)/(2|μ|/Δ+Σt)` — verified by reading the historical body via the closeout memo §"Before — 11 concepts".
9. Existing precedent in numpy: `np.cumprod`/`np.cumsum` IS associative scan with monoid `(·, +)` respectively. JAX `lax.associative_scan` generalises to any user monoid.

## ELEGANCE DETECTOR HITS

- **Smell #1** (repeated near-identical code across geometry variants): the cumprod fast path existed only for slab; the curvilinear bodies couldn't share it because the M-M angular state was perceived as breaking the affine recurrence. The frame attack reveals it does NOT break it — the angular state is per-cell exogenous data, not part of the chain state.
- **Smell #6** ("iterative" without fixed-point structure analyzed): the per-cell fold body's loop carries a *scalar* state across cells. The fact that this is associative-scannable (and therefore parallel-scannable) is not yet named in the codebase.
- **Smell #8** (nested loops with algebraic structure): `for n: for visit:` with a scalar accumulator threading the inner loop is the textbook scan pattern. The Pythonic outer-loop over `range(N)` ordinates is the architectural anti-pattern Phase G is fighting.

## FRAME CANDIDATES

### Frame 1 — Associative-scan / pair-monoid (Blelloch 1990)

**Trigger**: Feature 1 (scalar affine recurrence) + Feature 5 (pair-monoid associativity).

**Reformulation**: Define the per-cell update as the monoid element `m[i] = (a[i], b[i])` with composition `(α₁,β₁) ⊕ (α₂,β₂) = (α₁·α₂, α₂·β₁+β₂)`. The sweep result is `ψ[i+1] = (Π m[0..i]) ⊕ (1,ψ₀)` where Π is the scan over `⊕`. The numpy implementation IS `chain_dag_scan` (cumprod/cumsum path). The JAX implementation is `lax.associative_scan(monoid_op, m)`.

The closed-form derivation: write `m[i]` as the 2×2 affine-projective matrix
```
M[i] = [[a[i], b[i]],
        [0,    1   ]]
```
applied to the augmented vector `(ψ, 1)`. Matrix multiplication composes them. The scan is the prefix-product of these matrices. The two cumulatives (cumprod_a, cumsum(b/cumprod_a)) are the two non-trivial entries of the prefix-product matrix; the numpy form decomposes the matrix scan into two scalar scans because the 2×2 has special lower-triangular structure.

**Elegance payoff**:
- **Structure-exposing**: names the per-cell update as an *element of a monoid*, not a "coefficient pair". The composition rule `⊕` IS the algebra of forward substitution on (L+C). The cumprod/cumsum decomposition becomes a *theorem* (2×2 lower-triangular factors), not an implementation trick.
- **Algorithmic-advantage**: any associative scan admits a `O(N)` sequential implementation (current numpy cumprod/cumsum) AND a `O(N log N)` parallel implementation (Blelloch up-sweep + down-sweep, Brent's theorem). Same primitive, two backends — sequential CPU now, parallel GPU later if profile demands.
- **Structurally-simpler**: subsumes slab's pre-Step-2.5 cumprod (special case) AND the curvilinear per-cell fold (same monoid, longer `b` expression). One primitive, two call sites.

**First test**: implement `ordinate_scan(a, b, ψ₀)` as the numpy cumprod path; add `tests/sn/test_ordinate_scan.py::test_monoid_associativity` that randomises `a, b`, splits the array, scans both halves, composes with `⊕`, and checks bit-equality (modulo ULP) against the unsplit scan. Pass condition: associativity holds to FP-non-associativity bound. This test discriminates the frame from the "just a numpy trick" framing because monoid associativity is what *justifies* the cumprod/cumsum decomposition.

**Structural attack on current**: the candidate `chain_dag_scan` signature exposes `(a, b, ψ₀)` as three numpy arrays. This is the implementation, not the math. The math is "a sequence of monoid elements, scanned". The current signature hides that `(a[i], b[i])` is one *element* and that the scan operator `⊕` is the math entity. Naming `chain_dag_scan` over `ordinate_scan` is also a leaky abstraction: "DAG" is true (single chain IS a DAG) but it suggests the primitive will extend to general DAGs. It will NOT (see Frame 7).

**Precedent**: Blelloch 1990 §1.5 "First-order linear recurrences" — exact match. The recurrence `x[i] = a[i]·x[i−1] + b[i]` is the canonical example.

---

### Frame 2 — Segmented associative scan (Sengupta-Harris-Owens 2007)

**Trigger**: Feature 4 (2D wavefront has partial order — independent anti-diagonal levels) + Feature 7 (two-chain fan-in is NOT a chain monoid).

**Reformulation**: For 2D wavefronts, frame the sweep as: outer loop is a sequential scan over levels (cells on level k depend only on level k−1 outputs); inner is a *parallel reduce* within a level. Segmented-scan vocabulary (Sengupta) gives one primitive that handles both: scan with segment boundaries at level transitions, where within a segment the operation is "level update" (the cell receives two upstream face fluxes from level k−1, produces two outgoing face fluxes for level k+1).

**Elegance payoff**:
- **Expressive**: names what 2D actually does — *fold over levels*, not fold over cells. The current `update_batch(SweepCellSlice)` shape already implements this; the segmented-scan framing makes the structure legible.
- **Structurally-simpler**: ONE vocabulary (segmented scan) covers 1D (one segment, scalar monoid) and 2D (many segments, batched parallel reduce within each). The number of distinct fold concepts in the codebase drops from 2 to 1.

**Refuted as a unification path** for current code, however: the *operator* in 2D's per-level inner step is NOT the 1D monoid `⊕`. It's a parallel batched application of the 2D balance equation `(Q + sx·ψ_x_in + sy·ψ_y_in)/(Σt+sx+sy)` with scatter writes into the outgoing face buffers. The segmented-scan frame *names* this structure but does not collapse 1D and 2D into one numpy call. Forcing them into one call pessimises 1D (per the explorer memo §Q2).

**First test**: extend the L1 wavefront regression to expose the level boundaries as data (`level_boundaries: list[int]`) and verify that the per-level computation is invariant under permutation of cells within a level. Pass condition: scrambling cells within a level changes scatter-write ORDER but not the converged answer (modulo FP non-associativity). This proves the within-level operation is genuinely a parallel reduce, not a hidden serial scan.

**Structural attack on current**: `_sweep_2d_wavefront` calls `sweep_graph.apply(...)` which iterates levels but the level-fold structure is buried inside `SweepDependencyGraph.apply`. Hoisting it explicitly — `for level in levels: cell_update.update_batch(level_slice)` — would make the segmented-scan structure visible at the sweep-skeleton level.

**Precedent**: Sengupta-Harris-Zhang-Owens 2007 §3 segmented scan.

---

### Frame 3 — JAX-style typed associative scan (`lax.associative_scan`)

**Trigger**: Feature 5 (monoid associativity) + observation that the existing numpy cumprod/cumsum path IS a special-cased associative-scan implementation.

**Reformulation**: Replace `chain_dag_scan` with `associative_scan(monoid_op, elements, axis=0)` where `monoid_op` is the user-supplied `⊕`. For numpy backend, recognise the lower-triangular monoid and dispatch to cumprod/cumsum; for JAX backend, get the parallel scan for free.

**Elegance payoff**:
- **Algorithmic-advantage**: same source code runs CPU-sequential (numpy) and GPU-parallel (JAX) by choosing backend. Identical semantics by associativity (modulo FP non-associativity, already accepted under `vv-principles`).
- **Structure-exposing**: the `monoid_op` parameter is the algebra. Reading the call site reads as math.

**Refuted as immediate path**:
- **Pattern 6 violation**: introducing JAX requires (1) jax as dependency, (2) backend-dispatch infrastructure, (3) test matrix across backends. Phase G has NO concrete second consumer for this primitive yet (the cylinder/sphere/slab unification needs only sequential scan; the 2D path uses a DIFFERENT operator that this monoid does not capture). Two concrete future use cases are required; we have one (1D unification) plus an analogy to 2D that the analysis above refutes.
- The PURE-numpy implementation of `associative_scan` for arbitrary monoids does not match the cumprod/cumsum O(N) numpy ops speed — a pure-Python loop over the monoid op runs at slab's current per-cell speed (the regression we're trying to fix). For numpy backend, the cumprod/cumsum special case IS the right implementation.

**Cost**: jax dependency, backend-dispatch shim, test matrix. Bookkeeping cost dominates the (currently zero) algorithmic benefit.

**First test**: deferred. Reopen if (a) GPU performance becomes a profiled bottleneck, OR (b) a second consumer of the abstract `associative_scan` shows up.

**Structural attack on current**: none for this phase — the candidate `chain_dag_scan` is *already* the right primitive at numpy granularity. JAX is the next layer.

**Precedent**: `jax.lax.associative_scan` docs; jackd's RWKV cumulative-sums reformulation.

---

### Frame 4 — Triangular-solve primitive (`scipy.linalg.solve_triangular`)

**Trigger**: Feature 2 (forward substitution on block-triangular (L+C)).

**Reformulation**: Build the explicit `(L+C)` matrix per ordinate (sparse bidiagonal); call `scipy.linalg.solve_triangular` or `scipy.sparse.linalg.spsolve_triangular`.

**Elegance payoff (theoretical)**: a standard library primitive that GPU-vendor implementations (cuSolver) accelerate for free.

**Refuted**:
- **Algorithmic mis-match**: SN sweeps NEVER form the matrix. The per-cell coefficient builder runs in O(nx) and the back-substitution runs in O(nx); matrix formation is O(nx) but with materially higher constants (allocating COO/CSR, packing, calling LAPACK). For nx ~ 80 this is a regression, not a gain.
- **Loss of structure**: the algebra `cell_balance_terms` exposes the curvilinear coefficients with named intermediates (`numer_upstream`, `denom`). The matrix form hides that the coefficients ARE the streaming-term components.
- **Bit-identity**: `scipy.linalg.solve_triangular` uses LAPACK `STRSV`, which is bit-identity-different from the elementwise recurrence. Would require regression-snapshot regen for no win.
- **Multi-group blocking**: each ordinate × group pair is one independent bidiagonal solve. Batching across groups buys nothing in numpy.

**First test**: deferred. Closed.

**Structural attack on current**: none — the candidate primitive is *correctly NOT* doing this.

---

### Frame 5 — Tensor networks / MPS

**Trigger**: Feature 6 (per-cell state includes M-M angular state); higher-order spatial schemes (LD, EC) would extend the per-cell state dimension.

**Reformulation**: Treat the sequence of per-cell transition matrices as a Matrix Product Operator (MPO); the sweep is MPO application to the initial-state MPS.

**Refuted for current scope**:
- The per-cell state is **scalar in the chain direction**. The M-M angular state is *per-cell exogenous*, not chained. The bond dimension of the would-be MPO is 1 (trivially). MPS framing degenerates to the chain-scan we already have.
- The non-trivial bond dimension would emerge ONLY if we later treat the angular DAG and spatial DAG as a single 2D MPS over (cells × ordinates within a level). This IS plausible for higher-order angular closures, but it's not what 2.5/2.5b/2.6 is solving.

**First test**: deferred. Reopen if (a) higher-rank angular closures (rank > 2) become production AND (b) memory-vs-rank tradeoffs require compression.

**Structural attack on current**: none for this phase.

**Precedent**: trajectory_resolvent memo (related Variant α work has the MPO/MPS match for an open kernel; here the open structure is degenerate).

---

### Frame 6 — Recurrence transformation / blocked scan (Karp-Miller-Winograd 1967, Leiserson-Saxe 1991)

**Trigger**: Feature 1 + the question "is `chain_dag_scan` the principled output of recurrence transformation, or just one valid form".

**Reformulation**: The textbook blocked scan: split the chain into B blocks of size N/B, scan each block sequentially (sequential within = good cache behaviour), reduce across block boundaries (parallel), then forward-correct each block with the inter-block prefix. For numpy this is *equivalent* to cumprod/cumsum at FP-non-associativity. For GPU this is the Blelloch up-sweep/down-sweep tree.

**Confirmation**: `chain_dag_scan` IS the canonical output of recurrence transformation for the *sequential* case. The blocked form is the canonical output for the *parallel* case. They are equivalent under associativity. No deeper transformation reduces work below O(N).

**Elegance payoff (confirmation, not new structure)**:
- Confirms Frame 1's `ordinate_scan` is correctly named and at the right level of abstraction.
- Provides the future-extension path (block-parallel within a single ordinate) without changing the primitive's signature.

**First test**: same as Frame 1's monoid-associativity test — that test is the precondition for any blocked-scan implementation.

**Structural attack on current**: confirms the candidate; no separate attack.

**Precedent**: Blelloch 1990 §4 work-efficient up-sweep/down-sweep; Snyder 1986; Leiserson-Saxe retiming as the operations-research framing.

---

## UNEXPLORED

- **de Rham cohomology / FEEC** — no compatible-element trigger; SN is purely operator algebra at this layer.
- **Symplectic geometry** — no Hamiltonian / phase-space flow trigger; SN sweep is fixed-point iteration, not characteristic integration.
- **Bloch / crystallographic groups** — no lattice-periodicity trigger at the sweep level.
- **Wiener-Hopf** — no half-space convolution trigger; the per-cell update is local.
- **Hilbert-Schmidt / separable kernels** (A.7) — already validated at Phase 5 in a *different* context (continuous-µ multi-bounce); not applicable to the cell-fold layer.
- **Multiplication-operator spectral theorem** (A.3) — same as above; this is the matvec/eigenvalue layer, not the sweep layer.
- **Graph theory / DAG scheduling** (A.7) — confirms the structure but adds no new primitive beyond Frame 1; the *scheduling* axis is what enables Frame 2 (2D segmented scan), already attacked.
- **Krylov subspace theory** — applies at the outer iteration (SI → GMRES), not the inner sweep primitive.

---

## CROSS-METHOD POLLINATION

Current method: SN sweep (forward substitution on block-triangular).

**Borrowings**:

- **From MOC (chain-fold over rays)**: MOC's ray sweep is the *exact same* affine scan along the chord: `ψ_out = exp(−τ)·ψ_in + (1−exp(−τ))·Q/Σt` is the affine update with `a = exp(−τ)` and `b = (1−exp(−τ))·Q/Σt`. Naming `ordinate_scan` exposes that SN-1D and MOC are the SAME primitive at different chord-parameter granularity (cells vs ray segments). Trigger: both are affine chain folds. First test: implement MOC's per-ray sweep as a call to `ordinate_scan` with `a, b` built from `(τ, Σt, Q)`; verify the L0 fixed-source slab MOC benchmark matches to discretization error.

- **From CP (Neumann series = source iteration)**: A.7 Fredholm theory tells us SI is the Neumann series of the transport operator's compact resolvent. This is the *outer* loop; the per-sweep scan is the per-iterate transport application. Naming `ordinate_scan` as the inner primitive and `SourceIteration` as the outer Neumann-series driver (Phase G Step 4) makes the bridge legible.

- **From parallel sweep literature (Adams-Larsen-Brown KBA, Pautz)**: 2D/3D wavefront sweeps are already known in the SN parallel-computing literature; Frame 2's segmented-scan framing matches the KBA "diagonal pipelining" decomposition. The blocked-scan extension (Frame 6) maps onto the spatial-domain-decomposition layer that ORPHEUS does not currently implement.

---

## RECOMMENDED PATH FORWARD — locked design

Frame name: **Associative-scan with pair-monoid (Blelloch §1.5), sequential numpy backend**.

Key primitives:

```python
# orpheus/sn/spatial/scan.py — new

def ordinate_scan(
    a: np.ndarray,           # (nx, ng) — per-cell multiplier
    b: np.ndarray,           # (nx, ng) — per-cell additive
    psi_0: np.ndarray,       # (ng,)   — initial face flux
) -> np.ndarray:             # (nx, ng) — outgoing face flux at each cell
    """Sequential associative scan of the affine recurrence

        psi[i+1] = a[i] * psi[i] + b[i]

    Monoid: (a1, b1) ⊕ (a2, b2) = (a1·a2, a2·b1 + b2),
    represented as the 2×2 lower-triangular affine matrix
    [[a, 0], [b, 1]] acting on (psi, 1).  Decomposes into
    cumprod_a · (psi_0 + cumsum(b / cumprod_a)) — two scalar
    scans because the 2×2 has lower-triangular structure.

    Blelloch 1990 §1.5 (first-order linear recurrences).
    """
    cumprod_a = np.cumprod(a, axis=0)
    return cumprod_a * (psi_0 + np.cumsum(b / cumprod_a, axis=0))


class CellUpdate(Protocol):
    """As today, with one ADDED method."""

    def affine_coefficients(
        self,
        visits: list[CellVisit],   # the per-ordinate fold's visit sequence
        total_xs: np.ndarray,      # (nx, ng)
        source: np.ndarray,        # (nx, ng)
        angular_state: np.ndarray | None,   # (nx, ng) or None for slab
    ) -> tuple[np.ndarray, np.ndarray]:
        """Vectorised per-cell affine coefficients (a, b).

        The (a, b) pair is dual to update(visit, ...): if the
        per-cell update returns ``psi_out = update(visit, ...,
        psi_in)``, then ``a[i]·psi_in + b[i] == psi_out`` cell-wise.
        For slab DD with neutral curvature this is the classical
        cumprod fast path.  For curvilinear DD with angular state
        as exogenous data, this still produces (a, b) — the angular
        contribution is absorbed into b[i].
        """
        ...
```

Sweep skeleton (per-ordinate body) becomes:

```python
a, b = cell_update.affine_coefficients(visits, sig_t, QV, angular_state=psi_angle)
psi_face_out_per_cell = ordinate_scan(a, b, psi_face_in)
# Cell averages and outgoing spatial face fluxes are derivable in
# one vectorised pass from psi_face_in / psi_face_out_per_cell.
```

For curvilinear, the M-M angular state thread is computed in a separate vectorised pass *after* the spatial scan (the angular closure `ψ_a_out = (ψ_avg − (1−τ)·ψ_a_in)/τ` is itself an affine recurrence in the angular index — second `ordinate_scan` call). Two scans per ordinate per level: one spatial, one angular.

Why this is the right primitive:

1. **It collapses concepts**: slab cumprod fast path + curvilinear per-cell fold collapse to ONE primitive with ONE call shape.
2. **It is named at the math level**: `ordinate_scan` reads as "scan the affine recurrence along the ordinate's cell chain". The `monoid op` framing is in the docstring; the call site reads as math.
3. **It restores slab performance**: cumprod/cumsum are ~3 numpy ops per ordinate, exactly the pre-Step-2.5 cost.
4. **It extends to MOC** (cross-method pollination above): same primitive, different `(a, b)` builder.
5. **It does NOT pretend to handle 2D**: 2D's two-channel fan-in is genuinely a different operator. The naming `ordinate_scan` (not `chain_dag_scan`) does NOT lie about scope.
6. **Pattern 5 respected**: ONE primitive, one builder method on `CellUpdate`, no Union types, no dispatch.
7. **Pattern 6 respected**: no JAX, no parallel-scan infrastructure, no MPS — those are reopened only when concrete second consumers appear.

What this does NOT do:

- 2D wavefront stays as-is (Frame 2 named the structure but didn't unify the operator).
- No GPU. No JAX. No parallel scan. All available *under associativity* if profile demands them later.

---

## SELF-CORRECTION

None. The reclassification held throughout: no hedging on trigger matches, no balanced views, no closing register. The two refutations (Frames 4, 5) and one deferral (Frame 3) are bound by specific conditions (no second consumer; matrix-form regression; rank-1 MPS degeneracy) rather than "may not apply here".

---

## Memory updates (apply next skill revision)

1. **Promote Frame 1 to A.7 trigger table** as new row "First-order affine recurrence (Blelloch §1.5)":
   - Trigger: per-cell update of form `ψ[i+1] = a[i]·ψ[i] + b[i]` with per-cell exogenous `(a, b)`.
   - Levers: 2×2 lower-triangular monoid, cumprod/cumsum decomposition, blocked-scan parallel extension, JAX `associative_scan` backend.
   - Payoff: SN-1D slab fast path + curvilinear unification under one primitive (this attack); future MOC pollination; sequential-now / parallel-later seam.

2. **Promote Frame 2 to A.7 trigger table** as new row "Segmented scan (Sengupta-Harris-Owens 2007)":
   - Trigger: outer fold over independent batches/levels + inner fan-in/fan-out (not chain monoid).
   - Payoff: names what 2D wavefront IS; does NOT unify with 1D's chain monoid (REFUTED unification, KEPT vocabulary).

3. **Promotion candidate Smell #16**: "Per-cell scalar accumulator threaded through an inner loop with no associativity argument named." Specifically catches Issue #196 Phase G Step 2.5 where the unified fold body inadvertently retired the cumprod fast path because the monoid structure of the affine recurrence was not named. The smell fires when an inner `for` loop with `psi_in = result.outgoing_spatial_flux` lacks a docstring or skill anchor naming the scan structure.

4. **Cross-link**: this memo links to `[[issue-196-phase-g-step2-5-closeout]]` (which documents the cumprod regression) and `[[issue-196-phase-g-step2-5-further-collapse]]` (which arrived at the same Q1 unification path without naming the monoid).
