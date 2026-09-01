# ORPHEUS Symmetry & Quotient Layer: Theory, Placement, and Action Plan

**Scope.** Encoding group actions on axes, quotient axes, and the descent of measures, bases,
quadratures, kernels, and operators. Written against the state recorded in
`orpheus_operator_machinery_report.md` (v2), `moment_space_and_layering_plan.md`, and the
axis/space taxonomy (five continuous space kinds, retract lattice, `WithTrace` per-axis
descriptor map).

Date: 2026-08-31

---

# Part 0 — The decision, stated first

**Build the catalog and the descent machinery. Refuse the general orbit-space engine.**

Computing an orbit space from scratch is mechanical (invariant ring → syzygy ideal →
Procesi–Schwarz PSD condition), but the groups that ever occur in transport number about a
dozen. A Gröbner-based engine would be abstraction without a consumer, and its failure mode is
debugging elimination orderings instead of transport. Each catalog entry is *derived* by the
procedure once, then hard-coded with the derivation in the docstring and a symbolic
regression test.

**The layer is nearly free given the kernel/operator split.** This is the load-bearing
architectural claim of this document:

> Symmetry reduction is a **second binding of the same kernel**. `Λ` is representation-free;
> `.on(V)` supplies `(M, R)` from `V`'s frame. A quotient space `V/H` is just another `V`.
> Nothing in the kernel layer changes, and no operator is ever assembled-then-projected.

Corollary, and the reason for the placement ruling in Part III: the quotient is applied at
**Realization (stage 2)**, as a property of the space being bound to. Post-hoc projection of an
assembled operator pays the full cost and defeats the purpose.

**What justifies the layer** (three named consumers, in dependency order):

1. Curvilinear S_N — the quotient's connection term *is* the angular redistribution operator.
2. The compatibility gates (Part II) — a runtime check on quadrature/symmetry mismatch that
   catches a real error class currently invisible.
3. Isotypic decomposition — required by higher harmonics (Riesz contour) and reactor noise
   regardless; symmetry is where it comes from.

Absent all three, do not build this.

---

# Part I — The theoretical spine

## I.1 The derivation procedure (catalog generator, run by hand, once per entry)

For `G` compact acting orthogonally on `ℝⁿ`, `X ⊆ ℝⁿ` a `G`-stable real algebraic set:

1. Minimal generators `p₁,…,p_k` of the invariant ring `ℝ[x]^G` (finite by Hilbert–Weyl).
2. `p : ℝⁿ → ℝᵏ` is proper and separates orbits ⟹ `ℝⁿ/G ≅ p(ℝⁿ)`, closed semialgebraic.
3. Syzygy ideal `I = ker(ℝ[y] → ℝ[x], yᵢ ↦ pᵢ)` by elimination.
4. **Procesi–Schwarz** (Invent. Math. 81, 1985): `p(ℝⁿ) = { y ∈ V(I) : P(y) ⪰ 0 }` with
   `P_ij = ⟨∇pᵢ, ∇pⱼ⟩` re-expressed in `y`. `V(I)` gives equalities; `P ⪰ 0` gives inequalities.
5. Intersect with the ideal of `X`.

**Worked entry — `S²/SO(2)`.** Invariants `p₁ = z`, `p₂ = x²+y²`; empty syzygy ideal;

```
P = [[1, 0], [0, 4p₂]],    det P = 4p₂
```

so `ℝ³/SO(2) = {p₂ ≥ 0}`. Adjoining `p₁² + p₂ = 1`:

```
S²/SO(2) = { μ ∈ ℝ : 1 − μ² ≥ 0 } = [−1, 1],    μ = Ω̂ · ẑ
```

Three facts fall out that are load-bearing downstream:

- **The action is not free.** The poles have full stabilizer. `[−1,1]` is a
  manifold-with-boundary / orbifold, not a quotient manifold; `μ = ±1` is the **singular
  stratum**. Any design assuming quotients are smooth submersions is wrong exactly there.
- **`det P = 4(1−μ²)` is the squared orbit radius**, and the orbit-space boundary is where it
  vanishes. This is the same polynomial as the coefficient in the curvilinear angular
  redistribution term `(1/r)∂_μ[(1−μ²)ψ]`, and the same locus as the characteristic angular
  boundary where the trace measure vanishes. See §I.4 for the status of that identification.
- **`(SO(3), SO(2))` is a Gelfand pair**, and `S² = SO(3)/SO(2)`, so the object above is really
  the double coset space `SO(2)\SO(3)/SO(2)`. Bi-invariant functions on it are the zonal
  spherical functions = Legendre polynomials. **Funk–Hecke is Schur's lemma on a Gelfand
  pair.** `RΛM` with `Λ` diagonal is a theorem of harmonic analysis on the quotient, not a
  chosen factorization.

## I.2 Encode the contravariant side, not the point set

`C(X/G) ≅ C(X)^G` and `L²(X/G, π_*μ) ≅ L²(X,μ)^G`. ORPHEUS never touches orbit spaces; it
touches function spaces over them. Reasons to build the invariant side, in order of force:

1. Bases and quadratures live on function spaces, not orbit spaces.
2. The invariant side is finitely presented and computable; the topological side is
   well-behaved only because `G` is compact and buys nothing further.
3. **It is already in the taxonomy.** For `G` compact acting unitarily, the Haar average

   ```
   ℛf = ∫_G f(g⁻¹ · ) dg
   ```

   satisfies `ℛ² = ℛ` and `ℛ* = ℛ`: the orthogonal projection onto the invariant subspace. It
   factors as `ℛ = π* ∘ 𝔄` with `𝔄 π* = id`. **`(π*, 𝔄)` is a retract pair.** Quotient-by-
   compact-group is not a new relational kind; it is an instance of the retract lattice with
   the Haar average as retraction. No new relational type is introduced by this plan.

**Consequence — operator descent is retract conjugation**, exactly as blocks already are:

```
A_H = 𝔄 ∘ A ∘ π*
```

and `A_H` equals the restriction of `A` to the invariant subspace **iff `[A, U_g] = 0 ∀ g ∈ H`**.
That commutation is a per-operator condition, checked at Declaration, and it is where physics
enters: geometry, material field, and BC must each be `H`-equivariant.

**Quotienting is an axis-level operation, not a space-kind-level one.** `S²/SO(2)` with `dμ/2`
is still a bulk data `L²` space over a different axis. The five space kinds are untouched. The
one genuinely new piece of axis data is the **stratification descriptor** (§I.5).

## I.3 `ℛ* = ℛ` is another instance of "adjoint = dual iff tight"

`ℛ² = ℛ` is algebraic and always holds. **Self-adjointness does not.** `ℛ* = ℛ` requires the
Gram to be `H`-invariant, i.e. the group must act by isometries of the *discrete* inner product,
not merely of the continuum measure. On a non-tight frame, the adjoint and the dual diverge and
`ℛ` is a projection only for the dual pairing — the Reynolds "projection" becomes an oblique
projection in the adjoint pairing.

This is your existing cardinal principle recurring on a new object. It means the tightness gate
(Phase 1.4 of the operator machinery plan) and the symmetry gate are **the same gate evaluated
on a different operator**, and the implementation should share the check.

**Tightness descends for free.** If `{fᵢ}` is a tight frame for `V` with bound `A` and `P` is an
orthogonal projection, `{Pfᵢ}` is a tight frame for `PV` with the same bound. Taking `P = ℛ`:
a tight frame on the full angular space restricts to a tight frame on the invariant subspace.
No new tightness work is required for the quotient — the `B/A` ratio is inherited.

## I.4 Descent of the four kernel kinds

| Term | Kernel type | Descends iff | Notes |
|---|---|---|---|
| `S` | integral, `σ_s(Ω̂·Ω̂′)` | always (bi-invariant under all of `SO(3)`) | descent *is* Funk–Hecke; `Λ` diagonal in `P_ℓ` |
| `F` | rank-1, isotropic emission | material field `H`-invariant | `\|χ/4π⟩⟨νΣ_f\|`; trivially bi-invariant in angle |
| `C` | multiplier `Σ_t(r)` | `Σ_t ∘ g = Σ_t ∀g ∈ H` | pure material-symmetry condition |
| `L` | differential | mesh and geometry `H`-equivariant | **acquires a connection term** — see below |

Only `L` is nontrivial. Reducing a differential operator along a non-free group action produces
a term supported on the quotient's stratification. **Claim to verify (not to assert):** for the
diagonal `SO(3)` action on 1-D spherical phase space, the reduced streaming operator is
`∂_r + (1/r) ∂_μ[(1−μ²) · ]`, and the `(1−μ²)` is `det P` from §I.1. The arithmetic matches
exactly. Whether that is an identity of the phase-space reduction or a coincidence of the
1-D spherical case specifically is **Task Q0.3** and blocks Phase Q4's documentation, not its
code.

If it holds, the payoff is structural: the curvilinear angular redistribution term stops being
a geometry-specific hand-derived addition and becomes the connection of a quotient, derived
uniformly, with the `μ = ±1` degeneracy explained rather than special-cased.

## I.5 The singular stratum and `WithTrace`

`μ = ±1` is the fixed-point set of the `SO(2)` action. There, the orbit collapses to a point
and the pushforward trace measure vanishes. This is exactly the characteristic angular boundary
already flagged as the forcing consumer for the per-axis `WithTrace` descriptor map.

**The two are the same object.** The quotient axis must carry a stratification descriptor whose
degenerate-measure entries are precisely the fixed-point strata of the acting group. This
argues the per-axis descriptor map over any bulk/trace/full decomposition a second time and from
an independent direction: the degeneracy is *per-axis and per-point*, computed from the group,
not a global property of the space.

## I.6 Isotypic decomposition and the eigenvalue payoff

If `[A, U_g] = 0`, Schur block-diagonalizes `A` over the irreps of `H`. The quotient space is the
trivial isotypic component.

**Theorem (high confidence).** The `k`-eigenvalue fundamental mode lies in the trivial isotypic
component. Proof sketch: Krein–Rutman gives a unique non-negative dominant eigenvector `ψ`; `ℛψ`
is non-negative, non-zero, and an eigenvector with the same eigenvalue; uniqueness forces
`ℛψ = ψ`. **So `k` may be computed entirely on the quotient**, and this is a theorem, not a
modelling assumption. The same argument applies to the dominant `α` mode when it exists and the
prompt semigroup is positive (moderate confidence — the positivity of the prompt semigroup is
standard, the uniqueness argument needs the discrete `α` to be the Perron root).

**Higher harmonics live in the nontrivial isotypic blocks.** The Riesz-projection work
(Phase 4.1) gains an exact labelling: azimuthal modes are indexed by `SO(2)` irrep. This is the
strongest argument for building the isotypic layer at all — it is required work that symmetry
supplies for free.

**Ground field interacts.** Over `ℝ` the nontrivial `SO(2)` irreps are 2-dimensional rotation
blocks; over `ℂ` they are 1-dimensional (`e^{imφ}`). The isotypic decomposition is clean only
over `ℂ`. This is an independent argument for the `𝔽 ∈ {ℝ, ℂ}` space-level tag, and it means the
noise application and the isotypic layer should land together, not separately.

---

# Part II — The three well-posedness gates

These are distinct conditions on distinct objects and must not be fused into one flag.

| Gate | Statement | Object | Stage | Failure mode when violated |
|---|---|---|---|---|
| **G1 measure** | `H ≤ Isom(axis metric)` and `μ_axis` is `H`-invariant | axis | 2 | `π_*` ill-defined; `ℛ` not self-adjoint (§I.3) |
| **G2 discretization** | `H ≤ Sym(Q)` and the retained basis set is `H`-stable | frame / quadrature | 2 | discretization and quotienting do not commute; residue = **ray effects** |
| **G3 commutation** | `[A, U_g] = 0 ∀ g ∈ H`, per operator term | operator | 3 | `𝔄 A π*` is not the restriction; silently wrong answer |

**G2 is the one that catches real errors.** Level-symmetric `LQ_N` sets carry octahedral
symmetry. `C₆` about a hexagonal axis is **not** a subgroup of `O_h`. A 1/6-core hexagonal
reduction with a standard `LQ_N` set is therefore not well-posed, the commuting square fails,
and the residue is ray effects aligned with the symmetry mismatch. This is a live configuration
class at INL and, as far as I know, no production code checks it. If this layer ships nothing
else, ship G2.

**G3's group is computed, not asserted.** The admissible group is

```
H_max = Sym(geometry) ∩ Sym(material field) ∩ Sym(BC)
```

an intersection over the catalog, evaluated at setup. The user does not declare "this is a
quarter-core problem"; the code computes the largest admissible `H` and reports it beside the
acyclicity certificate. A user-supplied `H` is validated against `H_max`, never trusted.

---

# Part III — Placement in the pipeline and the layer stack

## Pipeline (7-stage, from the operator machinery report)

```
0  DATA          cross sections, geometry            → Sym(materials), Sym(geometry) computed here
1  KERNEL        Λ, representation-free              → UNCHANGED. Kernels do not know about H.
2  REALIZATION   kernel + frame + measure → arrow    → QUOTIENT LIVES HERE.
                 axes minted, partitions minted         OrbitAxis minted; G1, G2 evaluated;
                                                        retract pair (π*, 𝔄) minted;
                                                        stratification descriptor set
3  DECLARATION   A_loss = L + C − S − B              → G3 evaluated per term
4  POSING        Pencil(A, Π)                        → isotypic block selection (trivial block
                                                        for k; labelled blocks for harmonics)
5  STRATEGY      splitting recipe                    → recipe is computed on the quotient digraph
6  INSTANTIATION recipe applied to A − σΠ            → unchanged
7  SOLUTION      iterate / invert / normalize        → unchanged
```

**Ruling: there is no "symmetry reduction" step.** Reduction is not an operation applied to an
assembled operator. It is a choice of which space to bind to, made at stage 2. Any design that
introduces a `reduce(A, H)` function has committed the stage-inversion error the operator
machinery report catalogues (later-stage machinery bound into an earlier stage, or here its
mirror: a stage-2 choice deferred to stage 3+).

## Layer stack (import direction preserved)

```
L1  numerics/symmetry/     group catalog, Haar/Reynolds, orbit charts, pushforward measures,
                           retract-pair minting. Knows nothing of transport.
    numerics/basis/        H-stable sub-basis selection (Y_ℓm → Y_ℓ0 under SO(2))
L2  transport/symmetry.py  Sym(geometry), Sym(materials), Sym(BC), H_max intersection
L3  sn/, moc/, cp/         method-specific G2 checks (quadrature symmetry groups)
```

`numerics/symmetry` must not import from `transport/`. The determination "this reactor problem
has symmetry `H`" is L2. The catalog entry "`SO(2)` acting on `S²` yields `[−1,1]` with `dμ/2`"
is L1.

## Naming (each name earned by a distinguishing invariant)

| Name | Invariant that earns it | Alternative rejected |
|---|---|---|
| `IsometryGroup` | acts by isometries of the axis metric — the precise condition for G1 | `SymmetryGroup` (states no invariant); `CompactIsometryGroup` (compactness is an assertion on the class, not a distinction between siblings) |
| `OrbitAxis` | its index set is literally the orbit space | `QuotientAxis` (names the construction, not the object) |
| `ReynoldsProjection` | idempotent + self-adjoint-under-invariant-Gram = orthogonal projection; Reynolds names the Haar average | `Symmetrizer` (verb, no invariant) |
| `IsotypicDecomposition` | Schur block structure indexed by irrep | `SymmetryBlocks` (collides with the existing partition/block vocabulary) |
| *(no new type)* | `(π*, 𝔄)` is a **`Retract` instance**, not a new relation | `QuotientMap` — would duplicate the retract lattice |

**Compact groups only.** Periodic lattices (`ℤⁿ`) are isometries but non-compact; Haar averaging
over infinite orbits is not available. They enter the catalog as **finite cyclic groups on the
discretized index set**, which is what the mesh actually has. Document this restriction in the
`IsometryGroup` docstring; it is a real limitation, not an oversight.

---

# Part IV — Action plan, by implementation order

## Phase Q0 — Doc and derivation only. No code. Prevents rework.

| # | Item | Target | Gate |
|---|---|---|---|
| Q0.1 | **Catalog derivations.** Run §I.1 for each entry: `SO(2)/S²`, `ℤ₂` antipodal, `C_n` and `D_n` about an axis, the `O_h` subgroup lattice relevant to octant symmetry, `SO(3)` (diagonal action, 1-D spherical), `SO(2)×ℝ_z` (1-D cylindrical). Record generators, syzygies, `P`-matrix, resulting chart, pushforward measure. | `docs/theory/foundations/symmetry.rst` (new) | each entry has an explicit `det P` and named singular stratum |
| Q0.2 | **Write the Gelfand-pair reading of Funk–Hecke.** `RΛM` with diagonal `Λ` is spherical-transform-on-a-double-coset-space; the `Λ` entries are the spherical transforms. Supersedes any framing of harmonic projection as an optimization. | `frame.rst` / `scattering.py` docstring | Funk–Hecke cited as a consequence of Schur, not as a special function identity |
| Q0.3 | **Verify or kill the connection-term identification** (§I.4): is the curvilinear `(1/r)∂_μ[(1−μ²)ψ]` the connection of the `SO(3)` phase-space quotient, with `(1−μ²) = det P`? Derive the reduction directly; do not argue from the coincidence. | `derivations/` script + `symmetry.rst` | either a derivation, or an explicit "coincidental in 1-D spherical" ruling — **no third outcome ships** |
| Q0.4 | **Record the three gates (Part II)** as a posing-table precondition block, with the `C₆`/`LQ_N` incompatibility as the worked negative example. | `operator_algebra.rst` posing table | G1/G2/G3 stated as distinct conditions on distinct objects |
| Q0.5 | **Record the placement ruling**: no `reduce(A, H)`; quotient is a stage-2 space choice. Pre-commits against the fork an agent is primed to implement. | `symmetry.rst` + agent-memory file | ruling written before any code lands |
| Q0.6 | **Record `ℛ* = ℛ ⟺ H`-invariant Gram** as an instance of adjoint-vs-dual (§I.3), and tightness descent through `ℛ` as a free inheritance. | `operator_adjoint.rst` (the D3 naming-discipline block) | the quotient appears as a fifth entry in the adjoint-coincidence table |

## Phase Q1 + Q2 — Core machinery and its first consumer. **Land together; Q1 does not ship alone.**

| # | Item | Depends on | Gate |
|---|---|---|---|
| Q1.1 | `numerics/symmetry/groups.py`: `IsometryGroup` catalog from Q0.1. Each entry carries: action on the axis index set, Haar measure (or counting measure for finite), fixed-point strata, and a symbolic test reproducing its Q0.1 derivation. | Q0.1 | every catalog entry's chart and measure reproduced by its own regression test |
| Q1.2 | `numerics/symmetry/reynolds.py`: `ReynoldsProjection` as a first-class operator with fixed domain/codomain, constructed via `.on(V)` like every other bound arrow. Factored as `π* ∘ 𝔄`. | Q1.1; frame hierarchy (exists) | `ℛ² = ℛ` and `𝔄π* = id` to machine precision (T0) |
| Q1.3 | `OrbitAxis` minting: `Quotient(axis, H) → (OrbitAxis, Retract)`. The retract pair is an **instance of the existing retract machinery** — no new relational type. Stratification descriptor populated from `H`'s fixed-point set. | Q1.1, `WithTrace` descriptor map | retract lattice tests pass with the new instance; no new `Relation` subclass added |
| Q1.4 | `H`-stable sub-basis selection in `numerics/basis/`: `SphericalHarmonicBasis` → trivial-isotypic sub-basis `{Y_ℓ0} ≅ {P_ℓ}` under `SO(2)`. | Q1.1, `SphericalHarmonicBasis` (in flight, ERR-039 branch) | sub-basis is `H`-stable by construction; tightness `B/A` inherited, not recomputed |
| **Q2.1** | **Retrofit the 1-D slab angular axis as `S²/SO(2)`.** The existing hard-coded `μ`-only path becomes the zero-configuration instance of the general machinery. No new physics. | Q1.1–Q1.4 | **all existing 1-D results bit-stable**; the hard-coded path deleted, not shadowed |
| Q2.2 | **Normalization audit, forced by Q2.1.** The pushforward `dΩ/4π ↦ dμ/2` fixes the factor of `2π` exactly once. Prediction: the retrofit surfaces at least one inconsistent normalization convention between the quadrature weights, the moment definitions, and the scattering `Λ`. | Q2.1 | every angular normalization traced to a single pushforward statement; any discrepancy found is recorded as an ERR, not silently patched |
| Q2.3 | **Realization/quotient commuting-square test** (T1): `ℛ ∘ bind(K).on(V) == bind(K).on(V/H)` for `S`, `F`, `C`. This is the central correctness assertion of the layer and it is cheap. | Q1.2, Q2.1 | passes for all three kernel kinds; a deliberately `H`-asymmetric `Σ_t` field REDs it (negative control) |

## Phase Q3 — The gates. The highest-value-per-line item in this plan.

| # | Item | Depends on | Gate |
|---|---|---|---|
| Q3.1 | Quadrature symmetry groups: annotate each quadrature set in the tree with its symmetry group (`LQ_N` → `O_h`, Gauss–Legendre product → `SO(2)`-compatible by construction, Lebedev → its stated point group). | Q1.1 | every quadrature set carries a group; unannotated sets refuse to participate in a quotient |
| Q3.2 | **G2 check**: `H ≤ Sym(Q)` and retained basis `H`-stable, evaluated at realization. | Q3.1, Q1.4 | **`C₆` hexagonal + `LQ_N` REDs with a message naming the missing subgroup relation** (the flagship negative control) |
| Q3.3 | `transport/symmetry.py`: `H_max = Sym(geometry) ∩ Sym(materials) ∩ Sym(BC)`, computed at setup, reported beside the acyclicity certificate and the `τ` mesh diagnostic. User-supplied `H` validated against it. | Q1.1 | a problem with an asymmetric material patch reports the reduced `H_max`, not the geometric one |
| Q3.4 | **G3 check** at Declaration: per-term commutation `[A_term, U_g] = 0`. Cheap randomized test on a few `g`, plus exact structural checks where available (multiplier: field invariance; integral: bi-invariance of `Λ`). | Q3.3, Q2.3 | each of `L, C, S, B, F` reports its own descent verdict; the sum descends iff all terms do |
| Q3.5 | **Ray-effect residue metric**: when G2 fails but the user forces the reduction, quantify the non-commutation of the square. Pairs with the angular-entropy diagnostic already in the research tier. | Q3.2 | forced-`C₆` case produces a measurable, reportable residue rather than a silent wrong answer |

## Phase Q4 — Curvilinear S_N. The forcing consumer.

| # | Item | Depends on | Gate |
|---|---|---|---|
| Q4.1 | 1-D spherical angular axis as the diagonal `SO(3)` quotient; `OrbitAxis` carries the `μ = ±1` singular stratum with degenerate trace measure. | Q2.1, Q1.3, Q0.3 | the `WithTrace` descriptor at `μ = ±1` is degenerate-measure *by derivation from `H`*, not by hand |
| Q4.2 | The angular redistribution term as the quotient connection, **iff Q0.3 verified**. If Q0.3 comes back negative, this item degrades to "hand-derived term, singular stratum still supplied by the quotient" — the stratification payoff survives either way. | Q0.3, Q4.1 | α-recursion's degenerate starting direction lands on the singular stratum without a special case |
| Q4.3 | 1-D cylindrical as the `SO(2) × ℝ_z` quotient; chart, measure, and stratification from Q0.1. | Q4.1 | cylindrical and spherical share one code path differing only in catalog entry |
| Q4.4 | Curvilinear sweep admissibility on the quotient: the causal order must be compatible with the stratification (the degenerate direction is the recursion seed). | Q4.1, `sweep_acyclicity` | acyclicity certificate emitted for curvilinear geometries |

## Phase Q5 — Spatial quotients. Independent of Q4; can interleave.

| # | Item | Depends on | Gate |
|---|---|---|---|
| Q5.1 | **Reflective BC as chart boundary, not boundary operator.** A reflective face is the statement "the domain is a fundamental domain of `H`." On the quotient the identification is structural and there is no `B` term for that face. | Q3.3 | quarter-core case constructs with `B` absent on the reflective faces; answer bit-consistent with the full-core solve averaged |
| Q5.2 | **Monodromy rank under quotient.** Prediction: quotienting reduces closing-edge rank by `\|H\|` where the action is free on the closing edge set, less where it is not. **Measure, do not assume** — fixed strata break freeness and I do not have a clean theorem. | Q5.1, `sweep_acyclicity` SCC-rank pass (item 2.5) | rank table extended with a measured quotient column; `reflective\|reflective` 1-D slab under midplane `ℤ₂` reported empirically |
| Q5.3 | 2-D quarter symmetry (W+S reflective, E+N vacuum) derived as a quotient rather than a hard-coded acyclicity gate case; subsumes existing item 2.8. | Q5.1 | the 4-orthant DAG is produced by the quotient, one pass, zero boundary iterations |
| Q5.4 | Hexagonal `C₆`/`D₆` 1/6- and 1/12-core reduction, **gated on Q3.2 passing** — i.e. requires a hexagonally-compatible quadrature. Surfacing that requirement is the point. | Q3.2, Q5.1 | either a compatible product quadrature is supplied, or the configuration is refused with the reason |

## Phase Q6 — Isotypic decomposition. Depends on the pencil (existing Phase 3) and `𝔽 = ℂ`.

| # | Item | Depends on |
|---|---|---|
| Q6.1 | `IsotypicDecomposition(A, H)`: Schur blocks over irreps. Trivial block = the quotient space already built. | Q1.2, existing 3.1–3.2 |
| Q6.2 | `k` restricted to the trivial block by the §I.6 theorem, with the theorem cited in the posing table as a precondition-free reduction. | Q6.1 |
| Q6.3 | Higher harmonics labelled by irrep; Riesz contour per block. Collapses part of existing item 4.1. | Q6.1, existing 4.1 |
| Q6.4 | Reactor noise: azimuthal mode index = `SO(2)` irrep; requires `𝔽 = ℂ` for 1-dimensional irreps (§I.6). Lands **with** the ground-field tag, not before. | Q6.1, `𝔽` tag, existing 4.2 |

## Explicitly refused

| Item | Reason |
|---|---|
| General Gröbner/Procesi–Schwarz orbit-space engine | ~12 catalog entries ever; engine failure modes are worse than the cost it saves |
| `reduce(A, H)` post-assembly projection | stage inversion; pays full assembly cost; defeats the purpose |
| Non-compact groups as first-class | no Haar average; lattices enter as finite cyclic on the discrete index set |
| A new relational type for quotients | it is a `Retract`; the lattice already has the structure |
| A sixth space kind | quotienting is axis-level; the five kinds are untouched |
| Symmetry in Monte Carlo (tally symmetrization) | different mechanism entirely — research tier, not this layer |

## V&V battery, by tier

**T0 (algebraic, machine precision):** `ℛ² = ℛ`; `𝔄π* = id`; `ℛ* = ℛ` under the invariant Gram
(and *fails* under a deliberately non-invariant one — negative control tying to §I.3).

**T1 (structural):** realization/quotient commuting square (Q2.3); `π_*Q = Q_quot` on the
invariant basis with order of accuracy preserved; sub-basis `H`-stability; `B/A` tightness
inherited rather than recomputed.

**T2 (solve-level):** full-space solve averaged == quotient solve, to iteration tolerance
(Q5.1); 1-D bit-stability under the Q2.1 retrofit; monodromy rank measured (Q5.2).

**T3 (physics):** `k` from quotient == `k` from full space; forced-G2-violation ray-effect
residue measured and reported (Q3.5); curvilinear verification case against a known solution.

---

# Confidence ledger

**High:**
- The `S²/SO(2) = [−1,1]` derivation and `det P = 4(1−μ²)` as the squared orbit radius
  (computed above, not recalled)
- Reynolds operator as retraction; `(π*, 𝔄)` as a retract pair instance; no new relational type
- `ℛ* = ℛ` requiring an `H`-invariant Gram, as an instance of adjoint-vs-dual
- Tightness descent through an orthogonal projection (standard frame fact)
- Gelfand pair `(SO(3), SO(2))`; Funk–Hecke as Schur on zonal spherical functions
- `k` fundamental mode lies in the trivial isotypic component (Krein–Rutman uniqueness argument)
- G2 as a genuine, uncaught error class; `C₆ ⊄ O_h`
- The singular stratum `μ = ±1` coinciding with the vanishing trace measure
- Quotienting is axis-level, not space-kind-level
- Placement at stage 2; symmetry as a second binding of the same kernel

**Moderate — verify before the dependent item ships:**
- The connection-term identification (§I.4, blocks Q4.2's *documentation*, not Q4.1's code) —
  the arithmetic matches exactly but I have not derived the reduction
- Dominant `α` mode in the trivial isotypic block (needs the Perron-root argument for the
  prompt semigroup specifically)
- Monodromy rank dividing by `\|H\|` (Q5.2) — stated as a measurement, deliberately not a theorem
- That the Q2.2 normalization audit will find something. Strong prior from the `4π`/`2`/`1`
  convention collision, but it is a prediction

**I do not know:**
- Whether a hexagonally-compatible high-order angular quadrature with the right symmetry exists
  in usable form for Q5.4, or whether that item terminates in "refuse the configuration"
- How the `H_max` intersection behaves on unstructured meshes where `Sym(geometry)` is only
  approximately realized by the discretization — the continuum group may not be a symmetry of
  the mesh, and I have not thought through the tolerance question
