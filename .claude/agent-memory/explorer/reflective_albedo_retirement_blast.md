---
name: reflective-albedo-retirement-blast
description: Blast-set shape for retiring ReflectiveBoundary.albedo and refusing deck laws in LawScaled/LawSum (tree 0a5a23fa, 2026-10-01)
metadata:
  type: project
---

Census at `main` `0a5a23fa` (2026-10-01), full table in `scratch/boundary_ontology/reflective_cleanup_census.md`.

Durable shapes (re-verify before repeating):
- Two algebras read identically: `0.3 * law` builds a descriptor `LawScaled`; `0.3 * realize(law)` builds a Wave-0 `ScaledOperator`. The snapshot generator's "mixed" case is the second, so it is NOT a deck-composition site. Census the descriptor algebra with a runtime spy wrapping `LawScaled.__init__`/`LawSum.__init__` (frozen dataclasses: wrap `__init__`, a late `__post_init__` is never called).
- `SNProblem.reflective_axes` is an `isinstance(law, ReflectiveBoundary)` door: an RB(α<1) pair counts as a mirror loop, the albedo spelling does not. 1-D k is bitwise-blind to it (gauge needs two axes, #344).
- Retired-body transcriptions in tests (`_retired_diffusion_albedo`, `_pre_b32_face_action`) read `law.albedo` on RB rows: a field removal breaks them by AttributeError, not by assertion.
- Polymorphic mints of RB with albedo: `BoundaryTraceLaw.create("reflective", albedo=...)` and `cls(albedo=0.42) if "albedo" in fields` over the registry.
- The base `response_kernel` returns `None`: RB must keep an override returning `ScalarResponse(1.0)`.

**Why:** the reflective cleanup (W3) removes the albedo; these are the members a name grep misses.
**How to apply:** for any further deck/response carve, start from these doors and mints; see [[boundary-law-method-matrix]].
