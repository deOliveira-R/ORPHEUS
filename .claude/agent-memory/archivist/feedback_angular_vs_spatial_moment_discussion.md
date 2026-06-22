---
name: angular-vs-spatial-moment-discussion-and-spatial-moment-space-stubs
description: #240 D5b-S3 in discrete_ordinates.rst — author a NEW headline discussion contrasting the two orthogonal "moment" axes (angular HarmonicMomentField φ_ℓ^m vs spatial SpatialMomentSpace ψ̂) + expand the two committed S3 stubs (scattering-lift + SpatialMomentSpace). USER-requested headline ("a discussion of moments is a nice addition to the sphinx"). Sibling of affine-in-σ / scanmarch but on discrete_ordinates.rst (NOT loss_representations.rst) and PARTIAL (S3-A in-flight, not a close-out).
metadata:
  type: feedback
---

#240 D5b-S3 docs task on `docs/theory/discrete_ordinates.rst` (NOT the
loss_representations.rst page the affine/scanmarch siblings live on). THREE
deliverables: (1) a NEW headline H3 subsection `two-moment-axes` contrasting
the two orthogonal "moment" notions (the USER explicitly asked for this in a
design discussion); (2) expand the `ld-ubld-scattering-moment-lift` stub
(Σ_s⊗I scattering lift, PARTIAL/in-flight); (3) expand the
`spatial-moment-space` stub (the new first-class SpatialMomentSpace factor).

**Why this one is distinct from the affine/scanmarch siblings:**

- **The headline is a NEW from-scratch CONCEPT subsection, not a stub
  expansion.** The user said a discussion of moments (angular, spatial) "is a
  nice addition to the sphinx." The deliverable is the connecting tissue: the
  naming-collision warning + a contrast `.. list-table::` (resolves/basis/
  truncation-knob/count/typed-space/role/set-by) + the `φ_ℓ^m,p` tensor-product
  picture. CORE distinction to nail: angular = REPLACEMENT representation (hold
  ψ OR φ_ℓ^m, bridged by M/R); spatial = ADDITIONAL axis (rides on either).
  Lead with a `.. warning::` on the collision; place the new subsection BEFORE
  the scattering-lift stub (it is the tissue both stubs reference back to via
  `:ref:two-moment-axes`).
- **It is PARTIAL, not a close-out.** S3-A is IN-FLIGHT: the Σ_s⊗I lift LANDED
  but the φ̂ ITERATE that fills the slope rows is OWED (blocked on the
  typed-field widening → which IS the SpatialMomentSpace, the sibling stub).
  So the scattering-lift section needs a STATUS admonition (`.. admonition::
  :class: caution`) saying the lift scatters a slope source that is still ZERO
  in production (no field carries φ̂ yet) → the FP does NOT change until the
  iterate lands. Be honest about wired-vs-owed (the user's explicit "mark what
  is wired vs owed" instruction).

**The physics-completion-vs-iteration-only distinction (load-bearing, the
vv-principles Mode-9 nuance):**

- S2 solves `(L+C−S_flat)ψ=Q`, `S_flat=Σ_s⊗e₀e₀ᵀ` (scatter ONLY the cell
  average — slope rows zero) = O(h²) but diffusion-limit-INCONSISTENT. S3
  solves `(L+C−S_full)ψ=Q`, `S_full=Σ_s⊗I_spatial` (scatter every spatial
  moment) = diffusion-limit-CONSISTENT. DIFFERENT operators → DIFFERENT spectra
  → DIFFERENT fixed points. The converged flux CHANGES — that is the point (the
  thick-diffusion xfail flips xfail→pass BECAUSE the limit becomes correct).
- ⚠ This is NOT a Mode-9 FP-invariance change. Verifying S3 against the S2 FP
  would be the Mode-9 MIS-application (asserting FP-invariance of a change that
  legitimately moves the FP). The genuine Mode-9 invariant for S3 = SI-with-
  lagged-moment ≡ Krylov on the SAME `(L+C−S_full)`. Minted a derivation
  eq-label `ld-ubld-s2-s3-operators` for the S2/S3 operator contrast.
- The diffusion-limit WHY (a `.. note::` + `^^^^` subsection): the bilinear
  basis choice {1,x,y,xy} over simplex {1,x,y} AND the Σ_s⊗I lift are TWO
  HALVES OF THE SAME asymptotic argument (Adams 2001 simplex-LD discrete-
  diffusion-limit-invalid-on-quads; the xy cross-moment carries the limit; the
  lift makes sure the slopes/cross-moment receive scattering). Frame them as
  paired, not separate.

**Literature = PLAIN-TEXT inline, NOT [Key]_.** Adams2001 / LMM1987 / BLA1992 /
MRM2016 are NOT defined as `.. [Key]` bib entries ANYWHERE in docs/ (grep-
confirmed) — loss_representations.rst cites them plain-text ("Adams (2001)",
"Maginot, Ragusa & Morel (2016)"). So I wrote "Adams (2001)", "Border, Lewis &
Adams 1992", "Larsen, Morel & Miller 1987" inline. Using `[Adams2001]_` would
fire an undefined-citation warning (citations DO warn even without nitpicky).
(Decided NOT to mint new `.. [Key]` bib entries — keep consistent with the
existing plain-text-everywhere convention for these four refs.)

**Code-xref reality (matches the affine/scanmarch sibling memos):** this page's
`:meth:`/`:class:` render PLAIN-TEXT (page not member-autodoc'd; adding
automodule would dup-label the rich `:label:` docstrings). BUT the RST SOURCE
must name the OWNING class (Cardinal Rule 1, even invisible in output). CAUGHT
MY OWN BUG: I wrote `MaterialCrossSectionField` (4 sites) but the class is
`MaterialXSField` — `grep "^class " material_xs_field.py` is the gate, NOT the
build (`-W` does NOT catch a wrong-but-plain-text code-xref). ALWAYS grep the
real class name before citing. Other cited symbols verified: MomentProjection/
HarmonicMomentReconstruction (numerics.projection), SphericalHarmonicSpace,
HarmonicMomentField.from_mesh_and_L, AngularField/ScalarField (_bases),
D1ClosedForm.kernel_rhs (_ubld), find_factor (numerics.space),
spatial_moment_tail/average_moment_index (spatial_moment_space).

**Eq-label policy (matches siblings):** preserved the two committed
verifies-orphan labels (`ld-ubld-scattering-moment-lift`,
`spatial-moment-space-size`) UNCHANGED. Added 7 NEW derivation/representational
eq-labels (`two-moment-angular`, `-spatial`, `-tensor-product`,
`-carrier-space`, `ld-ubld-s2-s3-operators`, `spatial-moment-kronecker-order`,
`spatial-moment-append-policy`) — ALL untagged (no vv-status), they auto-appear
in matrix.rst's "Orphan equations" list right beside the affine/scanmarch
siblings (`loss-rep-affine` etc.). matrix.rst is AUTO-GENERATED → the orphan
bump is expected, NEVER hand-edit, NOT a build warning. Rationale for untagged:
they are algebra OF a within-narrative chain (the size law, the Kronecker
order, the S2/S3 contrast), not independent solver claims — same call as the
established siblings.

**`:eq:` vs `:ref:` discipline:** `harmonic-moment-projection` is an EQ-label
(`.. math:: :label:` in operator_algebra.rst) → cite with `:eq:`, NOT `:ref:`
(I initially wrote `:ref:harmonic-moment-projection` — would render plain-text
cross-doc but is WRONG; fixed to `:eq:`). `two-moment-axes` /
`spatial-moment-space` are SECTION anchors (`.. _name:`) → `:ref:`. Grep-gate
BOTH namespaces separately before building.

**Section-marker ladder = FILE-LOCAL `=/-/~/^`** (4 levels). The two stubs are
H3 `~~~~` (parent = the LD-UBLD H2 `----`). So I COULD and DID use H4 `^^^^`
sub-subsections inside each expanded stub — this DIFFERS from the scanmarch
sibling on loss_representations.rst which had NO `^` in its file (single-H3
only). Always scan the FILE's ladder first. Underline length = `len(title)`
code points (em-dash —, Σ, ⊗, § are 1 code point each) — verified all 14
title/underline pairs in the edited region with a python `len()` check.

**Build gate:** `-E -W --keep-going` baseline = EXIT 1 / 1 ERROR (the standing
`Mesh1D.from_geometry :paramref:` — sphinx-paramlinks not installed, out of
scope). Post-edit IDENTICAL = 1 (count-unchanged = pass). The 4
"has no matching equation node" lines are INFO (verifies-section-anchor
registrations), NOT warnings. Verified: 0 `.. todo::` remain in my region,
new section/eq anchors render in HTML, 7 new labels in matrix orphan list.

**Source-reading order that worked:** the SpatialMomentSpace MODULE DOCSTRING
(`spatial_moment_space.py`) was 80% of the headline-discussion prose seed (it
already carries the angular-vs-spatial contrast TABLE + the Kronecker order +
the orthogonal-axes framing). Then the two closeout memos
(`issue_240_d5b_s3_a0_..._closeout.md` for the typed-factor/Pattern-4/#207
story + the mutation-verified gate-has-teeth; `issue_240_d5b_s3_a_inc_c_...md`
for the physics-completion/BLOCKER/byte-id-rank-2-exact). Then the crosswalk
plan (`issue_240_d5b_s3_crosswalk.md`) for the FP-resolution + Mode-9 nuance +
the literature list. Then the CODE (`material_xs_field.py` einsum subscripts,
`_bases.py` _compose_spatial_moments, `_ubld.py` AVERAGE_MOMENT/face_moment_tail,
`space.py` find_factor) to ground every cited symbol. Closeouts/docstrings ARE
the prose seed; the code is the symbol-name + signature ground truth.

**Quality self-assessment (1-5):**

| Dimension | Score | Note |
|-----------|:---:|------|
| Derivation depth | 5 | full angular-vs-spatial contrast w/ table + φ_ℓ^m,p tensor-product picture; S2/S3 operator algebra; diffusion-limit WHY (slope rows load-bearing, bilinear+lift two-halves-one-argument); Kronecker order d=1/2/3; append-iff>1 ()-not-(1,) byte-id |
| Cross-references | 4 | every :eq:/:ref: resolves (grep-gated both namespaces); :meth:/:class: render plain-text by page convention (not a regression); fixed my own MaterialXSField class-name bug |
| Numerical evidence | 4 | byte-id-rank-2-exact negative control; mutation-verified gate-has-teeth (auto-read scheme → LD reds); foundation-test inventory; the thick-diffusion VALUE anchor named — but no convergence TABLE (S3-A in-flight, the value evidence is owed not landed) |
| Failed approaches | 5 | the bare-int-axis rejected (Pattern-4 table); the crosswalk's false "already broadcasts" assumption (.. warning::); the Mode-9 mis-application (verify-against-S2-FP would be wrong); why default-OFF even for LD |
| Code traceability | 5 | every symbol named to owning class + grep-verified; the 3 einsum lift subscripts table; the single-source delegation chain (spatial_moment_tail→face_moment_tail→_compose) |
| Derivation source | 5 | module docstring (the algebra-of-record for a typed-space half is the space module + the closeout memos, not a SymPy module) cross-checked against code @ 96dfc96 |

Weakest = numerical evidence (4) and cross-refs (4): numerical because S3-A is
in-flight (the diffusion-limit VALUE table is owed — I named the anchor + the
xfail-flip but the data does not exist yet, correctly); cross-refs because of
the page-wide plain-text code-xref convention (surfacing the whole transport/
sn package is its own arch-docs task — do NOT half-surface).
