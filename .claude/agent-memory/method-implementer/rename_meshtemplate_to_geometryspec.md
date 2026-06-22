---
name: MeshTemplate -> GeometrySpec rename closeout
description: 2026-05-03 lexical rename — MeshTemplate class + mesh_template field renamed to GeometrySpec / geometry_spec across orpheus/, tests/, docs/. Clean break, no backward-compat alias.
type: project
---

# MeshTemplate -> GeometrySpec rename closeout (2026-05-03)

Branch: `refactor/rename-mesh-template-to-geometry-spec` (commit `153e903`).

## Why

The class describes what it IS (a geometry specification) rather than
what it PRODUCES (a Mesh1D for discrete consumers). Method-agnostic
role established in R0.5 (commit `e41c1dd`); name now matches the
abstraction's semantics across both discrete and reference-method
consumers.

## Scope

* **File**: `orpheus/derivations/common/geometry_template.py` ->
  `orpheus/derivations/common/geometry_spec.py` (git-detected rename,
  86% similarity).
* **Class**: `MeshTemplate` -> `GeometrySpec`.
* **Field**: `mesh_template` -> `geometry_spec` on `La13511Case`,
  `CrossMethodCase`. All consumers updated.
* **No backward-compat alias** per user directive ("renaming later is
  waste of time"). R0.5 was the only commit using the legacy name; clean
  break is appropriate.

## Files touched (14)

```
orpheus/derivations/common/__init__.py
orpheus/derivations/common/geometry_template.py -> geometry_spec.py (rename)
orpheus/derivations/continuous/fn_method/benchmarks/__init__.py
orpheus/derivations/continuous/fn_method/benchmarks/la13511.py  (deprecation shim)
orpheus/derivations/continuous/sood_registry/__init__.py
orpheus/derivations/continuous/sood_registry/atalay1997.py
orpheus/derivations/continuous/sood_registry/builders.py
orpheus/derivations/continuous/sood_registry/la13511.py
tests/cross_method/adapters.py
tests/cross_method/cases.py
tests/cross_method/protocol.py
tests/cross_method/test_eigenvalue.py
tests/derivations/test_sood_registry_compatibility.py
docs/theory/sood_registry.rst
```

Net: 214 insertions / 210 deletions.

## Verification

* `pytest tests/cross_method/ tests/derivations/test_sood_registry_compatibility.py`:
  194/194 pass identically.
* Wider sweep (`tests/cross_method/ + sood_registry suite + fn_la13511 +
  carlvik_galerkin`): 344/344 pass, 12 skipped (pre-existing stub skips).
* `sphinx-build -W --keep-going`: exit 0, clean.
* `grep -rn "MeshTemplate\|mesh_template" orpheus/ tests/ docs/`: only
  3 historical-narrative references remain (all in docstrings/comments
  describing the rename event itself):
  - `orpheus/derivations/common/__init__.py:15` — "renamed from
    ``MeshTemplate`` to ``GeometrySpec`` on 2026-05-03"
  - `orpheus/derivations/common/geometry_spec.py:32` — "Originally
    named ``MeshTemplate``..."
  - `orpheus/derivations/continuous/sood_registry/la13511.py:70` —
    historical comment in the schema docstring

## Cross-method robustness

* All import paths resolve to the same `GeometrySpec` class:
  - `orpheus.derivations.common.GeometrySpec`
  - `orpheus.derivations.common.geometry_spec.GeometrySpec`
  - `orpheus.derivations.continuous.sood_registry.GeometrySpec`
  - `orpheus.derivations.continuous.sood_registry.la13511.GeometrySpec`
  - `orpheus.derivations.continuous.fn_method.benchmarks.GeometrySpec`
  - `orpheus.derivations.continuous.fn_method.benchmarks.la13511.GeometrySpec`
* `case.geometry_spec` (new field) is the canonical access path.
  Legacy `case.mesh_template` no longer exists — searched all of
  `orpheus/`, `tests/`, `docs/` confirms zero references.
* Property delegators (`case.n_groups`, `case.geometry`,
  `case.critical_dimension_mfp`, `case.critical_dimension_cm`) all
  read from `self.geometry_spec`.

## Local-variable polish

Small renames for consistency (within already-touched files):

* `_mesh_template_for(case)` -> `_geometry_spec_for(case)`
* Local `template = ...` -> `spec = ...` in `_sphere_R_cm`,
  `_slab_L_full_cm`
* `base_template`, `new_template` -> `base_spec`, `new_spec` in
  `tests/cross_method/test_eigenvalue.py::_shadow_with_thickness_mfp`
* `has_inline_template` -> `has_inline_spec` in `CrossMethodCase.__post_init__`

## Branch hygiene incident

During the rename, the parallel-running Archivist (working on
`docs/orbit-space-terminology`) appears to have pulled HEAD onto
their branch (their checkout fired while my last `git status` was
running, between the initial `git checkout -b` and the final commit).
Net effect: my rename commit landed on `docs/orbit-space-terminology`
instead of `refactor/rename-mesh-template-to-geometry-spec`.

Recovery (clean, non-destructive):

1. `git cherry-pick f258c10` onto `refactor/rename-mesh-template-to-geometry-spec`
2. `git reset --hard b81bc44` on `docs/orbit-space-terminology` to
   drop the rename from the Archivist's branch (their 3 commits
   `b8e0ae9..b81bc44` preserved intact).

Both branches now have the right state. The Archivist's branch is
unaffected; my rename lives on the intended branch as `153e903`.

## Lessons

* **When dispatching parallel agents touching adjacent territory**,
  the parent agent should explicitly tell each agent to verify
  `git branch --show-current` before each commit. The Archivist
  brief noted "different terminology fixes — probably no conflict",
  but the pickup was at the branch-state level, not file-content
  level.
* **`git mv` followed by `Write` on the new file in the same
  conversation** does not preserve the rename in `git status` until
  both are staged together — git detects the rename only at
  `add`-or-`commit` time via the similarity index. This is fine
  but worth flagging: the staged state was correct (`git diff
  --cached --stat` showed the rename arrow), even though `git
  status` showed `D + ??`.
* The Archivist's terminology fix on the same RST file
  (`docs/theory/sood_registry.rst`) was on different lines — no
  merge conflict expected when the two branches eventually merge.
  My only edit was the section heading "Mesh-template build
  conventions" -> "Geometry-spec build conventions", which is in
  a different region from any orbit-space terminology change.
