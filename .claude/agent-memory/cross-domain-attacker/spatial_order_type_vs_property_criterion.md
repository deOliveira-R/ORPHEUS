---
name: spatial-order-type-vs-property-criterion
description: VALIDATED criterion for property-vs-type — a representation earns a distinct TYPE iff ≥2 NON-canonically-iso bases coexist connected by a MODELED, APPLIED change-of-basis morphism (carries truncation error, has adjoint, in the operator algebra). Applied to ORPHEUS spatial-vs-angular moment asymmetry on branch feature/field-typed-operator-algebra. Verdict: spatial-order stays a PROPERTY (one tensor-Legendre basis, only morphism=id); angular-order is correctly TWO types (ordinate↔harmonic, M/R the modeled Vandermonde). Defer-with-trigger.
metadata:
  type: project
---

# Property-vs-type criterion (the sharpened "dual must coexist") — branch `feature/field-typed-operator-algebra`

DESIGN MEMO answering: should expanded SPATIAL order (within-cell DG degree)
be a TYPE symmetric to `HarmonicMomentField`, given angular order IS a type?
User's developing thesis CONFIRMED + SHARPENED, not refuted.

## The criterion (the durable kernel of this attack)

> A representation earns a distinct **type** iff there exist **≥2 bases that
> are NOT canonically isomorphic** (the iso depends on a quadrature/node
> choice), connected by a **change-of-basis operator that is itself MODELED
> and APPLIED** — it carries truncation error, has an adjoint, and
> participates in the operator algebra.

Three clauses; ALL must hold. This is the sharp form of "a dual must
coexist and not mix." Decidable by grep: count the within-X representations
and the applied non-id morphisms between them.

- **Angular order PASSES all three:** `AngularFlux` (collocation, N ordinates
  on S²) ↔ `HarmonicMomentField` (real-SH modal, (L+1)(2L+1)). Non-canonically
  iso (depends on the quadrature `Y_ℓ^m(Ω_n)`). M (`MomentProjection`) / R
  (`HarmonicMomentReconstruction`) are APPLIED morphisms with truncation
  content + adjoints, in `numerics/projection.py`. → TWO TYPES, load-bearing.
  `flux+moments` type-rejected (Field Layer-1 gate). Correct as-is.
- **Spatial order FAILS clause 1 (today):** ONE within-cell basis
  (tensor-Legendre DG, `SpatialMomentSpace`, `per_axis**ndim`). Only
  change-of-basis = identity (and `truncate`/inclusion WITHIN the tower,
  returns the same family). → PROPERTY, sufficient. `spatial_moments:int`
  field + `SpatialMomentSpace` factor composed onto BOTH angular types via
  `BulkField._compose_spatial_moments` (`_bases.py:157-202`). Correct as-is.

## Code ground truth (branch-verified, NOT Nexus — stale)

- `transport/fields/harmonic_moment_field.py:41-52` — "why distinct from
  AngularFlux": moment space vs per-ordinate space, M/R the legit route.
- `numerics/spaces/spatial_moment_space.py:21-28` — the two axes are
  ORTHOGONAL; spatial rides on EITHER angular rep. `:111-127` —
  `spatial_moment_tail` returns `()` not `(1,)` at DD (byte-identity).
- `transport/fields/_bases.py:460-472` — `MomentField` ABC is a THIN family
  marker: "a second moment representation would lift its shared machinery
  here (feedback_unify_after_two_instances)". THE PROJECT ALREADY LEGISLATED
  the user's thesis at the type-hierarchy level.
- archivist `feedback_angular_vs_spatial_moment_discussion.md:24` — the
  load-bearing distinction: angular = REPLACEMENT rep (hold ψ OR φ, M/R
  bridge); spatial = ADDITIONAL axis (rides on either). A unified
  `MomentField[Basis]` would ERASE exactly this → premature.

## Cross-method consumer verdict (does ANY ORPHEUS method supply the dual?)

NO current/planned method supplies the spatial nodal↔modal dual:
- **MoC** = flat-source ONLY (`moc/core.py:5,27`). No linear-source/region-moment.
- **Diffusion** = finite-difference POINT-VALUE (`diffusion/solver.py:3-4`).
  NOT nodal-expansion. grep `nodal|NEM|SANM` = ZERO source matches.
- **Nodal diffusion (NEM/SANM/ANM)** = THE canonical spatial dual (transverse-
  Legendre moments ↔ face partial currents, bridged by coupling coefficients =
  a modeled morphism). The STRONGEST for-case — but NOT planned (no roadmap
  entry, grep-confirmed zero).
- **LS-MoC / FET-MC** = modal-only spatial expansions → argue for SHARING the
  `SpatialMomentSpace` at `transport/` layer (it's already there), NOT a type.
- **Higher-order SN (QD)** = `per_axis=2→3` WIDENING in the SAME tower → the
  clearest evidence the axis is a parameter/fiber, the OPPOSITE of a new type.

## Foreign-frame convergence (where is spatial order a TYPE elsewhere?)

EVERYWHERE order is a PARAMETER except where a nodal dual coexists:
- hp-FEM degree p, p-multigrid, hierarchical-Legendre (Karniadakis-Sherwin),
  FET — order = PARAMETER (single-tower hierarchical `P_p ⊂ P_{p+1}`,
  prolong=inclusion, restrict=adjoint — within ONE rep).
- nodal-DG (Hesthaven-Warburton), nodal diffusion — TYPE, because nodal
  point-values coexist with modal coeffs, bridged by the APPLIED Vandermonde /
  coupling-coefficient morphism. This is the user's "dual must coexist" verbatim.

## VERDICT

AGAINST the type NOW; defer with trigger. Against is stronger: the type's sole
structural payoff (forbid-mixing a dual) has NO referent today — one spatial
rep, so the modeled morphism would be `id` (theatrics by the branch's own
"if type-hinting doesn't help it's theatrics" standard). Pose/diagnose
utilities split: slope-overshoot/under-resolution diagnostics are REAL but
attach to the existing space + a `BulkField` method (no type needed);
p-restriction/error-estimator are REAL-IF-BUILT but are typed OPERATORS within
one family (like `truncate`), not a new field type. Unification to
`MomentField[Basis]` = mathematically true (both are L²-Galerkin projection)
but premature — the right shared abstraction is the SPACE layer
(`FunctionSpace`/`TensorProductSpace`/`find_factor`, where SphericalHarmonicSpace
& SpatialMomentSpace are ALREADY siblings), not the FIELD layer.

## TRIGGER CONDITION that flips it (Pattern-6 defer-until-the-dual-exists)

1. (DECISIVE) A nodal/point-value within-cell OR face-current spatial rep
   enters production AND a MODELED, APPLIED nodal↔modal morphism is written
   between it and `SpatialMomentSpace` (Vandermonde for nodal-DG SN; coupling
   coefficients for nodal diffusion). Then lift `SpatialMomentField` + its
   nodal dual into the `MomentField` ABC, mirroring the M/R pair. Nodal
   diffusion is the most likely arrival.
2. A second method carries a within-cell modal spatial expansion AND it's
   shared at `transport/` — lifts the SPACE/ABC machinery (warrants the field
   type only WITH trigger 1's dual).
3. p-multigrid/p-adaptivity needs modeled `prolong`/`restrict` — flips toward
   typed OPERATORS, not a field type (the morphisms are canonical within one
   Legendre tower).

## Elegance criteria hit by the criterion frame

Structure-exposing (the property-vs-type call becomes a theorem about
coexisting-rep count) + Structurally-simpler (replaces an unbounded taste call
with a grep-decidable count). The criterion is a PROMOTION CANDIDATE to skill
Part C as an elegance smell: "expanded-order modeled as a TYPE when only one
basis exists (the change-of-basis is identity)" = type-theatrics; mirror of the
existing #16 shapes. Promote after a second independent sighting.
