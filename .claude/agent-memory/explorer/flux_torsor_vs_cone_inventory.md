---
name: flux-torsor-vs-cone-inventory
description: Durable shape of the #208 flux-torsor architecture and the cone-side footholds at HEAD 2026-08-19 — who enforces, who consumes, what superposition actually is in the tree
metadata:
  type: project
---

# Flux torsor vs cone — the durable shape (measured 2026-08-19, HEAD 7aae9bf1)

Full memo: `scratch/omr_v2_grounding/D_flux_ontology.md` (per-campaign; may be
consumed/deleted). This file keeps only what a future flux-ontology exploration
would re-derive at full cost. **The ontology ruling is the user's and was NOT made
as of this writing** — if the cone campaign has since landed, most of this
describes the RETIRED side; check `git log -- orpheus/transport/fields/_flux_role.py`
first.

**Why:** the plan orpheus-operator-machinery-report-v2 (§I.7, correction #15,
item 0.7) proposes overturning the torsor doctrine for a positive-cone ontology,
and scoped the fix as prose-sweeping; measurement showed it is an implemented
architecture, and the campaign cannot be scoped without this map.

**How to apply:** any brief touching flux typing, displacement types, FluxRole,
convergence diagnostics, or positivity/cone work starts from these facts instead
of re-measuring them.

## The torsor side (incumbent)

- **Enforcement is ONE mixin**: `FluxRole` (`orpheus/transport/fields/_flux_role.py`)
  — `⊖` mints the Rep-keyed sibling displacement, `⊕ displacement` is the torsor
  action, `flux + flux` TypeError, `affine_combination` (Σλ=1) the only blend.
  7 flux leaves ↔ 7 displacement leaves, paired structurally via
  `Displacement.sibling_of(Rep)` registry. Containers (`Composite`/`FullField`,
  `TimedFullField`, `RadialCharacteristicField`, `CoupledField`) deliberately do
  NOT pre-check member types — member dunders are the single source of truth.
- **Exactly 2 production consumers of the displacement TYPE**:
  (1) `numerics/iteration.py` `SourceIteration` — duck-typed `_DisplacementLeaf`
  Protocol (numerics must not import transport); feeds `contraction_ratios` +
  `last_displacement` DIAGNOSTICS ONLY — **the stop rides flat norms of
  source-typed rhs differences, never the displacement**;
  (2) `sn/acceleration/dsa.py` `DSACorrection.apply` — isinstance-guards its
  input and MINTS displacement-typed returns so `ψ ⊕ Δψ` is well-formed.
- **`affine_combination` has ZERO production callers** (definition + error
  message + tests only).
- **Superposition in the tree is SOURCE-level, always**: `rhs = q_ext + Σ g·ψ`
  (sources are a closed vector space); no site anywhere adds two flux-typed
  objects; solve-level linearity gates spell additivity torsor-legally via
  base-point independence (`B(ψ₁+σ)−B(ψ₂+σ) == B(ψ₁)−B(ψ₂)`, #331's P5 fallback).
  Krylov's V-space algebra runs on FLAT ndarrays (role erased at scipy boundary,
  rebuilt flux-typed from templates).
- **The incumbent's own concessions**: scalar scaling (ψ/k) deliberately legal,
  zero fluxes constructed freely, DSA docstring says "the swept vector IS the
  displacement from zero".

## The cone side (challenger footholds at HEAD)

- `is_positivity_preserving` on the scheme Protocol/ABC
  (`transport/spatial/scheme.py`): DD `False` (+ behavioral witness
  `tests/sn/sweep/core/test_diamond.py::TestPositivityFailure`), LD `False`,
  **0 production readers**; "gates negative-flux diagnostics" is an
  unimplemented docstring claim. **No Step / step-characteristics realization
  exists** — `class Step(... key="step")` is a docstring example only ⟹ any gate
  "flag True for step/SC" has no witness (plan-authoring §6c).
- Production cone refusal: `sn/boundary/realizer.py` refuses `ZeroFluxBoundary`
  (𝒜=−1) — "a negative inflow is outside that cone (ψ ≥ 0 ⟹ J± ≥ 0)".
- The COEFFICIENT family already states cone-as-predicate doctrine
  (`cross_section_field.py`: "nonnegativity is the cone, a property — not a type
  invariant") with its own battery (`TestCrossSectionConeAlgebra`).
- `power_iteration` already ray-normalizes (unit production rate, `flux / p`).
- Negative DD fluxes today: no fixup, no warning, nothing reads the flag.

## Cross-references

- **#331** is the sharpest statement of the tension: the three leaves of
  `A=(L+C)−S−B` disagree on whether their domain includes V (`L` accepts a
  displacement, `S`/`B` refuse). Under the cone ontology the disagreement
  dissolves (one tensor type); under the torsor it needs the A-arm/V-arm ruling.
- ⚠ `tests/sn/sweep/core/test_phase_c_gates.py` is issue **#168**'s Phase C —
  name collision with the plan's cone-battery "Phase C".
- Unaccelerated SI is the DEFAULT (`acceleration=None`); `power_iteration` keeps
  extensible NAMED criterion trajectories per outer — the CW-bracket attachment
  point, no new history object needed.
- Doc surface: `field_algebra.rst` (602 lines, §"Why affine" is the overturned
  argument, 4 `affine-*` equation labels all vv-status:documented) + 11 ref
  sites in `operator_algebra.rst` + [[affine-operator-split-convention]].
