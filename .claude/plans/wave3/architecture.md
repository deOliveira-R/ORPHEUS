# Wave 3 Architecture — Per-Paper Registries + Unified Meta-Registry

Plan returned 2026-05-03 by Plan agent. See main session for context
(post-compaction continuation of Wave 1/2 work).

---

## Executive summary

The Sood registry is well-factored: a tight 5-symbol public surface
(`Case`, `Truth`, `MeshTemplate`, `LA13511_CASES`, plus `build_*` /
`mixture_to_fn_arrays` / cache). Wave 2 has already smuggled non-Sood
papers (Atalay 1997) into `sood_registry/`, which is the violation
the user is asking us to undo. The plan below carves the existing
schema out into a paper-agnostic core, leaves Sood with only LA-13511
inside it, hosts each paper in its own sub-registry, and introduces a
thin meta-registry as an aggregator — without changing solver call
sites until ≥2 sub-registries are live (per
`feedback_unify_after_two_instances.md`).

---

## A. Folder layout

**Decision: meta-registry lives at `orpheus/derivations/registry/`**.

Layout:

```
orpheus/derivations/
├── common/
├── continuous/
│   ├── sood_registry/                 # Sood-only after Wave 3.0 cleanup
│   │   ├── __init__.py
│   │   ├── la13511.py                 # CASES + truth values for LA-13511
│   │   ├── builders.py
│   │   ├── extractors.py
│   │   └── cache.py
│   ├── kll_registry/                  # NEW Wave 3.1
│   ├── dahl_sjostrand_registry/       # NEW Wave 3.2
│   ├── atalay_registry/               # MOVED from sood_registry/atalay1997.py
│   ├── bis_registry/                  # NEW (Burkart-Ishiguro-Siewert 1976)
│   ├── nm_registry/                   # NEW (Neshat-Maiorino 1980)
│   ├── burkart_1971_registry/
│   ├── westfall_metcalf_registry/
│   ├── mitsis_registry/
│   └── ...
└── registry/                          # NEW: shared schema + meta-registry
    ├── __init__.py
    ├── schema.py
    ├── meta.py
    ├── canonical.py
    ├── cache.py
    └── solvers.py
```

## B. Shared schema

Promote to `derivations/registry/schema.py`:
- `MeshTemplate` (already paper-agnostic)
- `PaperCase` (renamed from `La13511Case`) with `paper_id`,
  `paper_table`, `primary_reference`, `notes`
- `Truth` generalised with `critical_thickness_mfp`,
  `critical_radius_mfp`, `extrapolated_endpoint_mfp`, `eigenvalue_c`
- `ParameterSpace` (new): `c`, `geometry`, `n_groups`,
  `scattering_order`, `bc_kind`, `multiplying`
- `SolverHandle` protocol

Defer: F_N legacy convenience properties (`case.sigma_t`) — keep on
`La13511Case` only; reassess after Wave 3.2.

## C. Meta-registry public API

`derivations/registry/meta.py`:
- `get(case_id)`
- `all_cases()`
- `by_paper(paper_id)`
- `by_parameter(c=..., geometry=..., ...)`
- `canonical_truth(case_id)`
- `agreement_matrix()` — cross-paper comparison
- `register_solver(paper_id, handle)`, `solver_for(paper_id)`

Auto-discovery via `pkgutil.walk_packages` (mirrors existing
`reference_values._build_continuous_registry`).

## D. Canonical-truth ordering

1. Most accurate published value (sig figs)
2. Earliest independent computation
3. Method most structurally different from consuming solver
4. Encoded as `canonical_priority: int` per case (data, not code)

`agreement_matrix()` returns (case_id × paper) table without
committing one as primary at the case-row level.

## E. Cache strategy

Per-paper cache, shared infra, namespaced by `paper_id`. Carve
`SoodResultCache` (515 LoC) → `RegistryResultCache` parameterised
by paper. Sood cache directory stays canonical for back-compat.

## F. Pytest / Sphinx integration

Two-tier fixtures:
- Per-registry (existing pattern: `LA13511_CASES`, future `KLL_CASES`)
- Meta-fixture in `tests/_harness/registry.py`:
  `meta_cases(paper=..., geometry=..., ...)` for cross-cutting tests

V&V matrix generator (`tools/verification/generate_matrix.py`) gains:
- "Per-paper registry coverage" section
- "Cross-paper agreement" section (renders
  `MetaRegistry.agreement_matrix()`)

## G. Order of operations

### Wave 3.0 — Cleanup (ZERO schema changes)
Move Atalay out of `sood_registry/`. Forwarding shim. Update Sphinx
split.

### Wave 3.1 — KLL 1974 sub-registry
Slab Table I + sphere Table V. Reuse existing `La13511Case` schema.
New tests `test_kll_registry_compatibility.py`.

### Wave 3.2 — Dahl-Sjöstrand sub-registry, then unify
DS Tables I + II. After 3.1+3.2 green: lift schema to
`registry/schema.py`, `La13511Case = PaperCase` alias, build
`MetaRegistry`.

### Wave 3.3+ — Remaining sub-registries
BIS, NM, Burkart 1971, WM, Mitsis as pure-data PRs (~1 day each).

### Wave 3.N — Cross-paper benchmarks (deferred)
Sphinx-rendered agreement matrix; cross-paper transcription smoke
tests.

## H. Test strategy (foundation tier, no solvers)

1. Schema invariants (every case has `materials`, `MeshTemplate`,
   `Truth`)
2. Cross-registry uniqueness of `case_id`
3. Auto-discovery integrity (parametrised count check per paper)
4. Canonical-truth resolution determinism
5. Solver-handle registration coverage
6. Cache namespacing (filesystem isolation)
7. Compatibility smoke (every case → `Mixture` + `Mesh1D`)

## I. Sphinx outline

`docs/theory/registry.rst` — top-level meta-registry page
`docs/theory/<paper>_registry.rst` — per-paper short pages (50-80
lines each)
`docs/verification/agreement.rst` — auto-generated cross-paper
agreement table

## J. Risk register

1. Atalay extraction breakage → back-compat shim, deprecation only
   after Wave 3.3
2. Schema lift accuracy regression → CI gate via verification matrix
3. Cache directory churn → keep old paths canonical for Sood
4. Sphinx build doubling → stub pages short, archivist defers

---

## Critical files

- `orpheus/derivations/continuous/sood_registry/la13511.py` (schema source)
- `orpheus/derivations/continuous/sood_registry/__init__.py` (back-compat surface)
- `orpheus/derivations/continuous/sood_registry/cache.py` (cache core)
- `orpheus/derivations/reference_values.py` (auto-discovery pattern)
- `tests/_harness/verify.py` (`vv_cases` entry point)
- `tools/verification/generate_matrix.py` (matrix renderer)

---

## Status

Plan approved by Plan agent. Implementation deferred until
literature-researcher (Atkinson) and numerics-investigator
(convention drift) results are in — those may surface schema
constraints not currently anticipated (e.g., if convention-drift fix
requires schema fields for anisotropy convention tag, the schema
lift in 3.2 absorbs that).
