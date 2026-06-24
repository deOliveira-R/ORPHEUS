---
name: projection-discipline-hierarchy-frames
description: Structural-detection verdict on the Galerkin/Petrov-Galerkin projection-discipline hierarchy in numerics/projection.py, for condensation/homogenization as Petrov-Galerkin on the existing Frame(basis,measure). VOCABULARY Π=analysis operator (Christensen)/projection/restriction, R∘Π=the projector (the idempotent — "projector" is RESERVED for it). HIERARCHY Galerkin & PG are SIBLINGS not Galerkin(PG) (Liskov fails — Galerkin STRENGTHENS Π*=R), but the deeper verdict is NEITHER earns a class. LATTICE one object = Saad's (K,L) trial/test pair = ORPHEUS (basis,measure); collocation=Dirac measure, LSQ=A*-pushed measure (LSQ needs A the frame lacks). VERDICT discipline = a PROPERTY (measure==basis.canonical/self-dual measure), the two ABCs are type-theatrics → retire, matches frame.py CAP_SOLVE iso/non-iso ruling + prior verdict 3.2. Read on main (clean, file:line ground truth — L-005 N/A).
metadata:
  type: project
---

# Projection-discipline hierarchy — `numerics/projection.py` (condensation/homogenization scope)

DESIGN VERDICT feeding the Petrov-Galerkin-condensation prototype. Read on
`main`, clean tree — `projection.py`/`frame.py`/`measure.py`/`basis/base.py`
file:line is GROUND TRUTH (no Nexus/worktree staleness; L-005 N/A). Extends
`projection_reconstruction_frame_pair.md` verdict 3.2 (Galerkin/PG = canonical/
oblique dual → ON THE FRAME, DEFER) one level UP: the ABC hierarchy itself.

## The five answers (durable)

1. **VOCABULARY — "projector" is reserved for the idempotent `R∘Π`.**
   - Π: V→W (fine→coarse, V≠W) is the **analysis operator** `T` (Christensen,
     *Frames & Riesz Bases* 2ed §1.1/§5.1 — `T:H→ℓ²` analysis, `T*` synthesis,
     `S=T*T` frame operator). It is NOT idempotent (different spaces) → calling
     it a "projector" is a category error. In FEM/approximation-theory idiom it
     is the **projection** or **restriction** `P` (Brenner–Scott §0.5/§3.4) —
     the "interpolation/projection operator" `Π_h`. The op `A_h = P A R`
     (`projection.py:8`) names it `P`, fine.
   - `R∘Π: V→V` IS the idempotent — THE **projector** (oblique in general,
     orthogonal in the Galerkin case). Saad *Iterative Methods* 2ed §1.12.4 +
     ch.5: an oblique projector onto `K` orthogonal to `L`; orthogonal iff `L=K`.
     For the SH frame `R∘Π` = the **band-limited projector** (`frame.py:35`,
     `projection.py:167`).
   - NAME RULING for ORPHEUS: keep `AnalysisOperator` for Π (Christensen-exact,
     already the abstract name `projection.py:114`). The idempotent `R∘Π` has NO
     type and needs none (it is a COMPOSITE, formed by the operator algebra) —
     do NOT mint `Projector`. If a name is wanted in prose, "band-limited
     projector" / "oblique projector" per the metric.

2. **HIERARCHY — Galerkin & PG are SIBLINGS, and the Liskov test KILLS the
   `Galerkin(PetrovGalerkin)` reading even though Galerkin ⊂ PG set-theoretically.**
   - Mathematically Galerkin = the constraint L=K specialization of PG (L≠K). So
     "Galerkin is-a PG" is TRUE as a set inclusion. But a type hierarchy encodes
     SUBSTITUTABILITY, not set inclusion. Galerkin ADDS a guarantee PG lacks:
     `Π* = R` up to the diagonal metric (`projection.py:169-176`, the
     `Π*=g_C·S0` adjoint-pair identity). A subtype that STRENGTHENS the
     postcondition is Liskov-LEGAL **only** for code that reads the base
     contract; but generic PG code is written to handle `Π*≠R`
     (`projection.py:203-204`) — it carries a SEPARATE adjoint path. Feed it a
     Galerkin instance and that path is dead weight / a redundant solve. The
     is-a is vacuous: nobody writes "PG-generic" code that a Galerkin instance
     usefully refines; they write code that NEEDS one or the other. → siblings,
     under the common `AnalysisOperator`. Current code (`:149`, `:190` both
     `(AnalysisOperator, ABC)`) is structurally RIGHT on sibling-vs-subclass.
   - BUT the deeper verdict (3+4) makes the question moot: neither earns a class.

3. **LATTICE of weighted-residual projections — one generalizing object, and it
   IS the (basis, measure) pair (= Saad's (K, L)).**
   - Weighted-residual unification: every method picks a **trial space** (where
     the solution lives, basis `K`) and a **test/weight functional family**
     (against which the residual is set to zero, `L`). Saad ch.5 names this the
     `(K, L)` projection pair verbatim; Galerkin = `L=K`, Petrov-Galerkin =
     `L≠K` (oblique). This maps 1:1 onto ORPHEUS: `K` = `Frame.basis` (trial,
     `frame.py:79`), `L` = `Frame.measure` (the test functionals = the Dirac-
     atom weights, `measure.py:124`). The frame docstring ALREADY says this
     ("basis fixes codomain/trial; measure fixes domain/test", `frame.py:14-20`).
   - The other weighted-residual members and where they sit:
     - **L²-orthogonal projection** = Galerkin with `L=K` under the L² measure.
       Special case, already the SH frame.
     - **Collocation** = test functionals are **Dirac deltas at nodes** →
       `L`=point-evaluation measure. EXPRESSIBLE as a `DiscreteMeasure` with
       unit weights at the collocation points (`measure.py` Dirac-atom def
       `:7-13` is literally this). FITS the (basis, measure) pair.
     - **Least-squares / minimal-residual** = test space `L = A*K` (test
       functions are the operator-adjoint of the trial functions; Saad §5.2.2,
       the "L=AK" minimal-residual class). This needs the **operator A** to
       build `L` from `K`. The (basis, measure) pair does NOT own A → LSQ is the
       ONE lattice member the bare Frame cannot express; it needs a
       `(basis, A·measure)` i.e. an A-pushforward of the measure
       (`measure.pushforward` `:538` is the right primitive, but the push is by
       A, supplied at solve time, not a frame-construction datum). RECORD this as
       the frame's expressivity boundary.
   - So: ONE generalizing object (the (K,L)=(basis,measure) pair) covers
     Galerkin, PG, collocation, L²-projection. Least-squares is adjacent — same
     pair shape but `L` is operator-dependent.

4. **DISCIPLINE = a PROPERTY of the (basis, measure) relation, NOT a type. The
   two ABCs are type-theatrics — RETIRE them.**
   - DECIDABLE CRITERION: the projection is **Galerkin iff
     `measure == basis.canonical_measure`** (the basis's self-dual / Plancherel
     measure — the measure under which the basis is its own canonical dual,
     equivalently the measure making the discrete Gram `T*T` a scalar multiple of
     identity on the band-limited subspace). Petrov-Galerkin iff
     `measure ≠ basis.canonical_measure` (a flux-weighted / spectrum-weighted /
     volume-weighted measure). This is a PROPERTY computable from the
     (basis, measure) pair the frame ALREADY holds — no class needed.
   - Type-vs-property gate (lessons L-004, coding-standards §"Type vs property"):
     mint a type iff ≥2 non-iso realizations + an APPLIED non-identity morphism
     that distinguishes them. Here: the morphism that would distinguish Galerkin
     from PG is `id` (both are the SAME `_FrameAnalysis.apply`, `frame.py:142`;
     they differ only in WHICH measure was bound, a DATUM not a code path). No
     consumer dispatches on `isinstance(op, GalerkinProjection)` — grep-checkable
     (the ABCs have ZERO concrete subclasses, `projection.py:206-210`; PG ships
     none). → both ABCs are CEREMONY: a conversion seam with no morphism reading
     it. RETIRE `GalerkinProjection` + `PetrovGalerkinProjection`; keep
     `AnalysisOperator`. The discipline becomes `frame.is_galerkin` (derived) or
     stays implicit in the bound measure.
   - This is CONSISTENT with `frame.py:27-35` which ALREADY made the SAME ruling
     for the SIBLING axis (iso vs non-iso = a `CAP_SOLVE` capability, "not a
     reason for a second mechanism"). Galerkin-vs-PG is the metric-of-the-measure
     analogue of iso-vs-non-iso: both are PROPERTIES of the given frame, both
     were wrongly heading toward a class. My prior verdict 3.2 deferred the MINT;
     this verdict says the existing ABCs are the premature mint to UNDO.
   - Smell #16 shape-2 (one concept in two homes): the discipline lives in BOTH
     a candidate class hierarchy AND the bound measure's relation to the basis's
     canonical measure. Collapse onto the measure (the primary object).

5. **CROSS-POLLINATION — the structural surprise.**
   - PN moment lift (`projection.py:35`), stochastic Galerkin / polynomial chaos
     (`:40`), MC adjoint moments (`:42`), SN sensitivity (adjoint-weighted inner
     products), energy condensation (`:36`), spatial homogenization (`:38`) are
     ALL one `Frame(basis, measure)` — the (K,L) pair with different `K` (SH /
     PN moments / PC polynomials / response functions) and different `L` (the
     measure: quadrature / spectrum-weighted / flux-volume-weighted / adjoint-
     weighted). Already asserted in `projection.py:23-43`; this analysis CONFIRMS
     it structurally.
   - **SURPRISE**: condensation/homogenization are the FIRST genuinely-oblique
     (L≠K) production instance — every prior Frame consumer (SH, PN) is Galerkin
     (L=K, canonical measure). So the thing that earns the "Petrov-Galerkin"
     distinction is NOT a subclass — it is the **non-canonical measure** (the
     within-group spectrum φ_g for condensation; the flux-volume weight for
     homogenization). The right place to encode "this is Petrov-Galerkin" is a
     spectrum-weighted / flux-weighted `DiscreteMeasure`, bound to the same
     coarse-energy / coarse-region basis. The PG-ness is a MEASURE, full stop.
     This is ALSO the unify-after-two-instances TRIGGER from prior verdict §10:
     condensation (§17) + homogenization (§18) are TWO oblique instances → they
     justify the oblique path, and they justify it AS a measure-construction
     (`spectrum_weighted_measure`, `flux_volume_measure`), not a class.
   - SN sensitivity is the adjoint-weighted sibling: `L` = the adjoint flux
     measure. Same shape; the "importance weighting" IS the test measure.

## Refuted frames (durable UNEXPLORED for this problem class)

- **Lattice/order theory (the literal "lattice of projections")** — the user's
  word "lattice" is the von Neumann lattice of projections on a Hilbert space
  (orthocomplemented, meet=∩ of ranges). It is REAL but LOW-SIGNAL here: the
  weighted-residual methods do NOT form a meet-semilattice under a natural
  order; they are POINTS in the (K,L) PARAMETER space, not elements of a
  projection lattice. The generalizing structure is the (K,L) product, not a
  lattice meet/join. Name the (K,L) pair; list lattice-theory UNEXPLORED.
- **Category theory / adjunction (Π ⊣ R)** — the analysis/synthesis duality is
  fully captured by the frame `T`/`T*` pair with the named diagonal metrics
  (`g_C`, `(2ℓ+1)`); no discriminating test the frame frame lacks. Vocabulary,
  not a lever. (L-001 pattern; same refutation as prior verdict.)
- **Tensor networks / MPO** — no bond-dimension knob; the basis tables are dense
  rank-full. (L-001.)
- **Spectral theory** — the Galerkin frame operator `S=T*T=4π·I` is a single-
  eigenvalue spectrum; `g_C` the only per-ℓ refinement, already named. A PG
  frame operator `T*_L T_K` is non-self-adjoint but its spectrum is not the
  lever — the (K,L) pair + the non-canonical measure is. No non-degenerate
  spectral payoff.
- **Differential geometry** — no curvature / `1/r` redistribution term; pure
  algebraic contraction on fixed spaces. (L-001/L-008.)

## First tests (each DISCRIMINATES — L-002)

- VERDICT-4 (discipline is a property, derivable from the measure): build a
  Frame with the SH basis and (a) its canonical quadrature measure, (b) a
  spectrum-weighted measure; assert `frame_a.is_galerkin` and NOT
  `frame_b.is_galerkin`, computed PURELY from `measure == basis.canonical_measure`
  — NO class tag read. A "discipline lives in the class" impl CANNOT produce this
  (both frames use the same `_FrameAnalysis` class → no class to read). The
  derivation must come from the (basis,measure) relation, which RED-s any
  class-based impl that can't tell the two frames apart.
- VERDICT-2 (siblings, not subclass): a NEGATIVE/substitutability test — generic
  PG code that builds the separate `Π*≠R` adjoint path, fed a Galerkin frame,
  must NOT silently produce a wrong/doubled adjoint. Equivalently: assert that
  for a Galerkin frame `Π* == R` up to `g_C` (`frame.py` `.H` gives `g_C·S0`),
  while a PG frame has `Π* != R` — a single `AnalysisOperator` + the bound
  measure decides it, no subclass. A `Galerkin(PetrovGalerkin)` hierarchy where
  Galerkin inherits the PG adjoint path FAILS (it carries the oblique path it
  should have specialized away).
- VERDICT-3 (least-squares is the expressivity boundary): assert that
  collocation (Dirac measure) and L²-projection (canonical measure) both
  construct as `Frame(basis, measure)` with NO extra argument, but a
  least-squares projection CANNOT be built without supplying `A` — the test that
  a `Frame.from_least_squares(basis, measure)` factory is impossible without an
  operator argument is the discriminator (a frame that claims to build LSQ from
  (basis, measure) alone is wrong — it silently drops the `L=A*K` push).

## Elegance assessment

- RETIRE `GalerkinProjection` + `PetrovGalerkinProjection` (and do NOT mint
  `Projector`): structurally-simpler (deletes 2 zero-subclass ABCs + the
  discipline-as-class seam), structure-exposing (the discipline becomes the
  THEOREM `measure==basis.canonical_measure`, computable, not a declared tag).
- The (K,L)=(basis,measure) recognition is expressive (one object covers
  Galerkin/PG/collocation/L²) + structure-exposing (condensation's PG-ness IS a
  non-canonical measure, not a method).
- Net: the prototype should build condensation/homogenization as a
  **non-canonical (spectrum/flux-weighted) DiscreteMeasure bound to a coarse
  basis via the existing Frame** — zero new operator classes. The two existing
  discipline ABCs are the premature mint to undo on the cleanup pass before the
  capability lands (coding-standards §"Clean before extending").
