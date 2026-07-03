---
name: frame-projection-coarsening-shape
description: Elegance rulings on the homogenize/condense coarsening machinery (frame projection campaign) — the axis-yields-views law, verb-layer-not-frame-layer symmetry, and the latent project() Gram misfire
metadata:
  type: project
---

The "frame projection machinery" campaign's two coarsening ops — spatial
homogenization (P3, merged) and energy condensation (P5, draft branch
`feature/sn-energy-condensation`) — reviewed 2026-06-27. Six durable rulings;
reuse these when ANY new coarsening/projection axis (CP/MoC homogenization,
GEC rank>0 #275) lands.

**Why:** the user explicitly re-opened the SHAPE of both (not bless-the-draft);
the abstraction was ~80% built but mis-located. The reshape is main-agent-steered.

**How to apply:** when a coarsening verb or a new frame consumer appears, check
it against these before crediting it.

1. **AXIS-YIELDS-ITS-VIEWS is the keystone law.** `Mesh` carries BOTH
   `volume_measure` (→DiscreteMeasure, source role) AND `indicator_basis`
   (→IndicatorBasis, target role) — and it is the UNARY axis object that does.
   The frame's two ingredients (a measure + a basis) come FROM the axis. Any new
   axis (EnergyGrid) MUST mirror this. The P5 draft VIOLATED it: the energy
   views were stranded on `GroupCondensation` (the BINARY fine+coarse object),
   so energy got-frame-from-pairing while space got-frame-from-axis — a live
   structural twin. Fix: `EnergyGrid.as_measure()`/`.as_basis()` (or match the
   spatial spelling `indicator_basis` for grep-symmetry). Retire/delegate the
   `GroupCondensation` views. See [[feedback-naming-consistency-greppable]].

2. **The collapse VERB lives on the XS CONTAINER, taking frames as ingredients
   — NOT on the frame, axis, or a free function.** `MaterialXSField.project_through(
   sigma_frame, emission_frame)` is the correct merged-P3 pattern: the field owns
   the channel→weighting taxonomy and routes Σ→sigma, χ→emission. P5's
   `Mixture.condense` re-enumerated the SAME XS channel list
   {SigC,SigL,SigF,SigP,SigT,SigS,Sig2,chi} that project_through enumerates = a
   twin-of-the-taxonomy. **RESOLUTION (P5.5 impl review 2026-06-27 — my original
   "dissolve `Mixture.condense` into project_through" was WRONG, layer-illegal):**
   `Mixture` is data, `project_through` is transport, and data↛transport — so a
   data-layer verb CANNOT route through the transport verb (the "single home"
   I named inverted the dependency arrow). The implementer instead kept BOTH
   verbs and extracted a shared **data-layer** assembler `Mixture.from_dense_channels(
   *, SigC=…,…, eg)` that BOTH `condense` (data) and `project_through` (transport→
   data return) call. This is SOUND and is the GENERAL pattern for a cross-layer
   coarsening twin: **when two verbs at DIFFERENT layers assemble the same
   container, extract the assembler at the CONTAINER's own (lowest) layer; do NOT
   route the lower verb through the higher one.** The cut is placed at the
   *mechanics* seam (csr-wrap + eg-thread → single-sourced; verified only two
   `csr_matrix(...)` sites in the whole XS surface, both inside from_dense_channels)
   and leaves the *weighting* taxonomy per-verb (energy marginalize-vs-average ≠
   space flux-vs-production frame.project) — correct, that asymmetry IS the domain
   content (different-axis-of-variation → no twin). RESIDUAL (downgraded to CONCERN-
   accept, NOT violation): the channel *enumeration* still appears as kwargs at 3
   sites, BUT from_dense_channels has an EXPLICIT keyword signature (not **kwargs),
   so a new channel is a LOUD pyright failure at every site, never a silent
   divergence — fails VIOLATION-leg-2. Leave it; the static-guarantee trade-off
   beats a dynamic `fields(Mixture)` iteration. Return-type asymmetry
   (homogenize→mesh-coupled MaterialMesh; condense→portable dict[int,Mixture]) is
   principled — two methods on one container is correct.

3. **`HomogenizationFrame`/`CondensationFrame` as verb-carrying PG subclasses is a
   FALSE SYMMETRY — mint NEITHER.** The asymmetry is structural and the code's own
   "one frame = one test weighting" law forces it: homogenization is intrinsically
   TWO frames (flux-weighted Σ + production-weighted χ), condensation is ONE frame
   + a bare `chi @ table` sum (χ over energy birth-groups is identity-preserving,
   needs no averaging). A HomogenizationFrame would bundle 2 weightings into 1
   frame = violates the law. Concept-count test: minting both ADDS two types
   without removing any (count goes up = wrong abstraction). The TRUE symmetry is
   at the VERB layer (`<XS-container>.project_through(frames)`), both consuming PG
   frames they do NOT subclass. `HarmonicFrame(GalerkinFrame)` is the ONLY
   legit frame subclass — it adds CARRIER-TYPED verbs forced by layering
   (numerics can't import transport carriers), not a coarsening verb. The χ
   frame-vs-table-sum asymmetry is the residual/displacement-style "different axis
   of variation → different mechanism → no twin" — document it deliberate in BOTH
   docstrings so nobody "unifies" it.

4. **`FrameBase.project` is the diagonal/PoU SPECIAL CASE wearing the general name
   — TYPE the precondition, do NOT build the dense path. [CLOSED in P5.5 impl
   2026-06-27 — exactly as ruled: `GramStructure` enum {DIAGONAL,PARTITION_OF_UNITY,
   DENSE} on `Basis` (base DEFAULTS to DENSE=safe-refusal), `FrameBase.gram` raises
   `MissingCapability` on DENSE, ~9 lines, NO dense (MR)⁻¹M solve built. The #20-
   family prose-doing-the-type-system's-job tell is gone.]** `project` = unconditional `gram.apply_inverse_metric(analysis(field))`;
   `gram` = unconditional row-sum probe `analysis(reconstruction(ones))`. Correct
   ONLY for disjoint-support (diagonal Gram) OR partition-of-unity trials — and
   the docstrings ADMIT this in prose while the code never checks it (a #20-family
   tell: prose doing the type system's job). A tapered/GEC-rank>0 basis (#275)
   bound into the frame calls project, gets the WRONG G silently, fails no assert.
   Latent because no violating basis EXISTS yet (gated by absence — #275 removes
   the gate). Fix per Pattern 4 ∩ Pattern 6: declare Gram-structure
   (DIAGONAL/PARTITION_OF_UNITY/DENSE) on the BASIS, have project/gram RAISE
   MissingCapability on a DENSE-Gram basis — ~5 lines, NO dense solve built.

5. **`Mixture.eg` is bare-numpy across the boundary (anti-pattern #13, LIVE).**
   `eg: np.ndarray | None`; `Solution.condense` does `EnergyGrid(material.eg)` —
   re-wrapping at the consumer. A `(N+1,)` descending boundary array type-checks
   identical to any float array. `Mixture.energy_grid` property is the min fix but
   a half-measure if `eg` stays bare (the EnergyGrid.__post_init__ descending/≥2
   invariant never runs on stored data). Honest fix: `eg: EnergyGrid | None`,
   parsed once at compute_macro_xs ingest. Larger blast radius — separable
   correctness carve, flag to user.

6. **`OverlapBasis` for SPATIAL homogenization = abstraction-ahead-of-evidence,
   DEFER.** Energy NEEDS fractional overlap (421→WIMS-69 genuinely non-nested);
   space does NOT — `homogenize` CONTRACTUALLY requires coarse cells = contiguous
   unions of fine cells (validates boundary alignment), so one-hot IndicatorBasis
   is exactly right and OverlapBasis buys nothing. Lifting it now abstracts over
   the DIFFERENCE (the nesting assumption) with ONE witness (energy) = the
   unify-after-two-by-abstracting-the-difference trap. The seam is already right
   (the basis is the only thing that changes; project_through is table-agnostic) —
   the lift is FREE when a real non-aligned-coarse-mesh consumer arrives (CP/MoC
   assembly homogenization). Leave a discoverable pointer in IndicatorBasis.

Minor: the `within_group is None` optional-arg dance is spelled 3× (EnergyGrid.
condense_to, GroupCondensation.from_grids, the dataclass default_factory) — a
Pattern-2 micro-twin (NIT, all produce identical objects). The default_factory
already supplies the default; condense_to/from_grids are ceremony. Keep at most
`from_mixture` and only if it does real work (deriving source grid from Mixture)
— a one-liner once EnergyGrid has the dual views (ruling 1).

**The THREE NICE-TO-HAVES (enumerated 2026-06-27 on committed HEAD `364f7b8` —
my original PASS closeout promised "three nice-to-haves" but only detailed the
CONCERN, so they were never listed for the main agent). All strictly optional,
none gates merge; ranked, #1 highest-leverage. Settled shape NOT re-litigated.**
1. **`OverlapBasis` carries a redundant `edges_per_axis`** (`overlap_basis.py:43-74`,
   inherited from `IndicatorBasis`, used "only for n_cells + space" per its own
   docstring) — but a fractional table's coarse count IS `overlap_table.shape[1]`,
   so it is a 2nd spelling of the coarse-group count (Pattern-2 micro-twin). Cost
   lands at the ONE ctor site `energy_grid.py:246` as the awkward
   `coarse.as_basis().edges_per_axis` round-trip (builds a throwaway IndicatorBasis
   to extract synthetic index edges it doesn't want). Fix: override `n_cells` →
   `int(self.overlap_table.shape[1])` on OverlapBasis + drop the `as_basis()` reach.
   Pays Pattern-2 + Pattern-3 together. (An `OverlapBasis` IS its table; its coarse
   dim is the column count.)
2. **`condense` vector-path lacks the morphism comment its siblings have**
   (`mixture.py:333-353`): the matrix collapse has `# sink g_to MARGINALIZED…source
   g_from AVERAGED…` (:338) and χ has its marginalize rationale (:350-351), but the
   five vector `average(self.SigX)` calls (:343-347) have no inline "AVERAGE/rate
   arm" anchor. Symmetry-in-math-→-symmetry-in-comment + Pattern-7 (these comments
   ARE the marginalize-vs-average crosswalk). Fix: one inline comment above the five
   vector calls mirroring the matrix one. No code change.
3. **`_overlap_table` clamp held in an undomain'd local `ce`** (`energy_grid.py:262-265`):
   the clamp ("coarse edges clamped to the fine span") is the load-bearing capture
   guarantee but `ce` reads identical to raw `coarse.edges`; a future vectorizing
   edit grabs raw edges, drops the clamp, reintroduces the mass leak. Fix: rename
   `ce` → `clamped_edges`. Pure local rename. Lowest delta of the three.
REJECTED as out-of-scope (recorded so not re-chased): (a) `indicator_basis.py:19-27`
still says homogenization is "Galerkin/L²(φV)-orthogonal" — REAL drift vs the PG
ruling + `frame.py:47-56`, but PRE-EXISTING P3 code about SPACE (homogenization),
belongs to #268 doc sweep NOT this energy pass; (b) `wims.py` has no `__all__`
(others do) — module hygiene, not a reads-like-domain improvement; its surface is
re-exported via the package `__init__`. Skip both unless user opens a hygiene sweep.
