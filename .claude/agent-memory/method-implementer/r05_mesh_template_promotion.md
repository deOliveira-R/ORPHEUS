---
name: R0.5 MeshTemplate promotion closeout
description: 2026-05-03 R0.5 refactor — promoted MeshTemplate to derivations/common/geometry_template, refactored cross_method adapters to consume it directly, added alpha_from_bc mapper.
type: project
---

# R0.5 — MeshTemplate promotion + cross_method refactor (2026-05-03)

**Branch**: `refactor/r05-mesh-template-promotion` (4 commits, off `main`).
**Plan**: `.claude/plans/trajectory_resolvent_hindsight_refactor.md` +
`.claude/scratch/geometry_handling_unification_audit.md` §Q5.

## Why

The geometry-handling audit identified `MeshTemplate` as a
half-complete unification: it was already method-agnostic in
content but lived inside `sood_registry/la13511.py`, where the
trajectory_resolvent and fn_method adapter layers couldn't
naturally consume it without a circular feel. The audit's Q5
rollout proposed a 5-step plan; this work landed steps 1-4, with
step 5 (multi-region MeshTemplate) deferred per
`feedback_unify_after_two_instances.md` until a 2nd instance
arrives.

## What landed (4 atomic commits)

1. `refactor(common): promote MeshTemplate to derivations/common/geometry_template`
   — Class moves to a new module; sood_registry's la13511 re-exports
   for backward compatibility. No call-site changes. `field` import
   dropped from la13511 (no longer used post-extraction).
2. `refactor(tests): inline materials/mesh_template fields on CrossMethodCase`
   — New optional fields + `__post_init__` validation. Closed-sphere
   k_inf case migrated from notes-based parsing to inline materials +
   `MeshTemplate(geometry="sphere", bc_left/right=BC.reflective)`.
   `_closed_sphere_params` helper deleted.
3. `refactor(tests): trajectory_resolvent adapters read mesh_template directly`
   — Drops `_sphere_R_truth_mfp` and `_slab_a_truth_mfp` (the truth →
   cm re-derivation via σ_t). Adapters now read
   `mesh_template.critical_dimension_cm` (sphere) and
   `mesh_template.domain_extent_cm` (slab) directly via the new
   `_mesh_template_for` helper. The cross-method agreement tests
   shadow `mesh_template` (not `truth_value`) when feeding F_N's
   predicted thickness to trajectory_resolvent.
4. `refactor(tests): add alpha_from_bc mapper, drop hardcoded BC↔α constants`
   — `alpha_from_bc(bc: BC) -> float`: vacuum→0, reflective→1,
   partial→albedo, others→`NotImplementedError`. Three trajectory_
   resolvent adapter classes lose their `alpha: float` field;
   `case.mesh_template.bc_right` is now the source of truth.

## Validation

* All 84 cross_method tests pass identically across all 4 commits.
* `sphinx-build -W` clean across all 4 commits.
* `tests/derivations/test_sood_registry_compatibility.py` (110 tests)
  + `test_fn_la13511_kinf.py` + `test_fn_la13511_slab.py` +
  `test_fn_la13511_sphere.py` all clean — no regression in upstream
  consumers of MeshTemplate.

## Bit-equality preservation note

The slab `L_full_cm` value drifted by ~6e-7 (from
`2*0.93772556/0.32640 = 5.7458674...` to `2*2.872934 = 5.745868`)
because the registry's published `critical_dimension_cm` is rounded
to 7 digits while `truth_value/sigma_t` carries full precision. This
drift is well below the cross-method test tolerance of 5e-5 on
k_eff; all 84 tests pass identically with the new code path. The
audit's recommendation to read mesh_template directly was the right
call — eliminates the unit-conversion redundancy and trusts the
registry's stored cm value.

## Validation rules on CrossMethodCase

`__post_init__` enforces:
* `materials` set ⇒ `mesh_template` MUST be set (no partial inline).
* `materials` set ⇒ `registry_case` MUST be None (one XS source).
* `mesh_template` alone alongside a `registry_case` IS allowed (the
  "override" path: registry XS + inline geometry shadow). Used by
  cross-method agreement tests.

## Coordination

Worked alongside R1 (parallel `power_iterate_variant_alpha` driver
in `trajectory_resolvent/`). Zero file overlap: R0.5 touched only
`orpheus/derivations/common/`, `orpheus/derivations/continuous/sood_registry/`,
and `tests/cross_method/`. R1 touched only `trajectory_resolvent/`
and added `tests/derivations/test_trajectory_resolvent_power_iterate.py`.

## What did NOT land (Step 5 — deferred)

The reflected-slab cases (4 cases in `REFLECTED_SLAB_CASES`) still
use `case.notes` + `_parse_notes_kv` + `_reflected_slab_params` for
their `(c_core, c_reflector, reflector_half_thickness_mfp)` triple.
The right shape for these is a multi-region `MeshTemplate` carrying
per-region `Mixture` plus the radial breakpoints. Per
`feedback_unify_after_two_instances.md`, multi-region MeshTemplate
lands when a 2nd consumer (e.g. a multi-region trajectory_resolvent
sphere adapter) arrives.

## Files touched

- `orpheus/derivations/common/__init__.py` (new export)
- `orpheus/derivations/common/geometry_template.py` (NEW)
- `orpheus/derivations/continuous/sood_registry/la13511.py`
  (MeshTemplate moved out, re-exported)
- `tests/cross_method/protocol.py` (added inline fields + validation)
- `tests/cross_method/cases.py` (closed-sphere case migrated)
- `tests/cross_method/adapters.py` (helpers refactored, alpha_from_bc added)
- `tests/cross_method/test_eigenvalue.py` (cross-method agreement
  tests use `_shadow_with_thickness_mfp` to mutate mesh_template)
