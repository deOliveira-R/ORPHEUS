---
name: Cross-method test protocol design idioms
description: Pattern for shared cross-solver regression infrastructure consuming an existing registry case schema rather than reinventing it
type: feedback
---

When designing cross-method regression infrastructure that consumes
an existing case registry (e.g. `La13511Case`), reuse rather than
reinvent. Build a thin `CrossMethodCase` wrapper that holds:

- a pointer back to the registry case (`registry_case: La13511Case
  | None`)
- per-solver tolerances `tolerances: dict[str, float]` (adapter
  name → absolute tolerance)
- a `pillar` tag (closed-form / MMS / semi-analytical / ancillary)
- a `claim_layer` tag (eigenvalue / flux-shape / convergence-order)
- a `truth_source` string (literature citation, surfaced in error
  messages)
- a `notes` field that doubles as a key=value parameter store for
  cases without a registry entry (NM 1980 reflected slab, closed-
  sphere k_inf fixtures)

**Why:** the canonical mistake here is to write a "richer"
parallel schema that duplicates the registry case fields. That
loads the migration cost (when wave3 lifts La13511Case to
PaperCase, the parallel schema has to migrate too) and creates
a second source of truth that drifts. Reuse + thin wrapper is the
clean factoring.

**How to apply:**

1. Identify the registry the cross-method infra consumes. For
   transport: `sood_registry.la13511`. The registry case carries
   `materials: dict[int, Mixture]` + `mesh_template` + `truth`.
   Don't recreate that schema.

2. Adapter dataclass shape: frozen dataclass with `name`,
   `method`, `geometry` fields + a `solve(case) -> ScalarResult`
   method. Each adapter handles unit conversions (mfp ↔ cm,
   half-thickness ↔ full slab, mfp ↔ cm via σ_t division)
   internally so test code stays declarative.

3. Cross-method agreement tolerance MUST be `max(tol_a, tol_b)`
   — tighter is reference contamination per `vv-principles` §6
   AI failure mode #6. Implement as an `agreement_tolerance(case,
   a, b)` helper.

4. V&V tagging convention from existing shipped cross-method tests
   (`test_fn_la13511_slab_xverif.py`,
   `test_fn_la13511_sphere_xverif.py`): tag both truth gates and
   pairwise agreement gates as **L1**, NOT L4. The conceptual
   level per `vv-principles` is L4 ("two solvers agreeing produces
   zero correctness info"), but the codebase convention is L1
   because each method's individual L1 truth-backing makes the
   pairwise agreement L1-strength evidence when structural
   independence is genuine.

5. Foundation gates on the protocol metadata: every case has at
   least one tolerance, every named adapter is registered, every
   pillar ∈ {closed-form, MMS, semi-analytical} (rejecting
   "ancillary" for truth values), case_ids are disjoint. Plus a
   coverage-matrix diagnostic that prints the (case × adapter)
   matrix on every run.

6. **Don't aim for L4 marker registration.** The harness only
   knows L0..L3 + foundation. Use L1 for both truth gates and
   pairwise agreement gates — the existing convention. If L4
   becomes useful as a marker later, register it then.

7. **Honest one-sided coverage**: when only one method natively
   supports a case (reflected slab — fn_method only;
   multi-region fixed-source — trajectory_resolvent only),
   declare the case with a single-adapter tolerances mapping and
   document it as one-sided. Don't pretend coverage that doesn't
   exist.

8. **Multi-group ≥2G coverage gap is real.** Bare-critical
   slab/sphere is inherently 1G (neither fn_method nor
   trajectory_resolvent ships native multi-group critical-
   dimension solve). Document the gap. The closed-sphere α=1 V_α1
   identity is the cleanest path to ≥2G cross-method coverage:
   k_eff = k_inf at any group structure with α=1, no leakage,
   rank-1 isotropic eigenmode. Extending the closed-sphere
   adapter to 2G is the natural next step.

9. **Truth values must trace to primary citations.** When
   transcribing from memory (e.g. NM Table 1 reflector cases),
   verify against the literature memo (`literature-researcher/`)
   before writing the case — I made up two truth values from
   faulty recall and the tests caught it.

10. **Test count target**: ~80-100 collected tests is reasonable
    for a multi-method, multi-geometry regression net. Most are
    parametrised over case lists; the test-row count is much
    smaller (~10 def's) but each parametrize blows it out.
