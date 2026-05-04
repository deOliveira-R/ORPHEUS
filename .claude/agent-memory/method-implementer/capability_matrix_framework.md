---
name: Capability matrix framework — meta-generator + auto-discovery
description: 2026-05-03 (HEAD `045afec`) atomic 3-commit landing of the capability-matrix infrastructure. Renames peierls→peierls_nystrom, adds fn_method matrix, replaces per-method generators with a single auto-discovery meta-generator + foundation-tagged regression tests.
type: project
---

# Capability Matrix Framework — 3-commit landing

`feature/peierls-greens-cylinder` 2026-05-03 commits `e2f2e82 → 045afec`.
Plan: `.claude/plans/capability_matrix_framework.md`.

## What landed (in order)

1. **`e2f2e82`** — Rename Peierls capability matrix generator to
   peierls_nystrom (filename + module name + Sphinx hook + theory-page
   xref). The list-table title preserves "Peierls" as the math label;
   the filename uses `peierls_nystrom` as the package label. 6 files,
   27 insertions / 27 deletions.

2. **`f46e050`** — Add capability matrix to fn_method: new
   `orpheus/derivations/continuous/fn_method/cases.py` exposing
   `capability_rows()` (9 rows covering 1G/2G/mG k_inf, slab F_N
   bare-critical, sphere F_N bare-critical, NM 1980 reflected slab,
   cylinder stub) + interim per-method generator
   `tools/verification/generate_fn_method_matrix.py` + Sphinx hook +
   "Capabilities at a glance" section in fn_method.rst. 5 files,
   500 insertions.

3. **`045afec`** — Meta-generator
   `tools/verification/generate_capability_matrices.py` auto-discovers
   `cases.py:capability_rows()` across `orpheus.derivations.continuous`
   via `pkgutil.iter_modules`. Schema contract:
   - **Required** keys: `name`, `geometry`, `n_groups`, `n_regions`,
     `bc`, `status`.
   - **Optional** auto-detected: `r0_over_R`, `closure`, `accuracy`,
     `scattering_order`, `multiplying`, `topology_class`. A column
     appears in the rendered table iff at least one row carries the
     key.
   `--check` mode compares regenerated output to on-disk; exits 1 on
   drift. Two per-method generators deleted; conf.py hook collapsed
   from 2 → 1. 3 foundation-tagged tests at
   `tests/derivations/test_capability_matrices.py` (sync check,
   determinism, minimum-discovery floor). 8 files, 637 insertions /
   415 deletions.

## Schema-contract trade-off

The plan called for `bc` + `status` as required keys, but
`peierls_nystrom/cases.py` shipped without them. Resolution:
**augment the registry**, not weaken the contract. Added `bc` (e.g.
"white (rank-2 per-face)", "vacuum + F.4 cavity closure") + `status`
("shipped (registry-anchored)", "shipped (k_eff gate pending —
Issue #104)") to all 13 peierls_nystrom rows. Net effect: the
peierls_nystrom matrix gained 2 columns vs the pre-rename render —
but the data was always knowable, just not surfaced.

The plan's "regenerate identically (no behavioural change)" line is
satisfied at the **row-set** level (same 13 references, same
left-half columns) but NOT at the byte level (new BC + Status
columns). The discrepancy is justified by the schema benefit.

## Backtick formatting

The original peierls per-method generator wrapped names with
double-backticks (`` ``peierls_slab_2eg_2rg`` ``) inside its renderer.
The meta-generator's `_format_cell` is field-agnostic (passes through
verbatim), so the backticks now live in the row data itself
(``peierls_nystrom/cases.py:capability_rows()`` returns
`"name": "``peierls_slab_2eg_2rg``"`). This pushes formatting choice
to the source-of-truth function, where each package controls its
own name presentation. fn_method names are plain prose (no
backticks); peierls_nystrom names are package-internal symbols
(backticked).

## Key insight — the silo problem

The fn_method package shipped a far richer set of references than
peierls_nystrom (k_inf cases, bare-critical slab + sphere, reflected
slab) but had **zero discoverability surface**. The reason is exactly
what the plan called out: the per-method generator pattern (a
`tools/verification/generate_<x>_matrix.py` file + a `conf.py` hook +
a `cases.py:capability_rows()` function) was undocumented anywhere
that crossed the FN development path. The agent who built fn_method
never saw it.

The meta-generator inverts the burden: now adding a capability matrix
to a new method is a 2-step zero-tools change:

1. Add `capability_rows()` to the package's `cases.py`.
2. Add `.. include:: _<package>_capability_matrix.inc.rst` to the
   theory page.

The Sphinx build does the rest. Future packages
(singular_eigenfunction, galerkin_spectral, trajectory_resolvent,
flat_source_cp) can be backfilled in ~30 min each per the plan's
implementation order.

## V&V pillar interaction

These are foundation tests — the meta-generator is documentation
infrastructure, not a numerical solver. The tests verify a
**software invariant** (cases.py is the source of truth, the matrix
on disk matches it, generation is deterministic, both packages are
discovered) — they do NOT carry a `vv_level` or `verifies(...)` mark
because there is no equation to verify. The
`@pytest.mark.foundation` tag is correct per the V&V architecture.

## Architecture notes for future agents

- The meta-generator is in `tools/`, not in `orpheus/derivations/`.
  Build-tool ↔ runtime category separation per Plan agent's Wave-3
  namespace reservation.
- Wipe-then-write semantics: at the start of every run, every
  `docs/theory/_*_capability_matrix.inc.rst` is removed. A package
  whose `cases.capability_rows()` is deleted will see its matrix
  disappear automatically. `--check` mode also flags stale files
  (on-disk matrices with no matching package).
- Discovery uses `pkgutil.iter_modules(orpheus.derivations.continuous)`
  with `pkg.ispkg` filter — modules without `__init__.py` are
  skipped. Imports use `importlib.import_module(f"{pkg.name}.cases")`
  with ImportError catch — a package without `cases.py` is silently
  skipped.

## Files

- `tools/verification/generate_capability_matrices.py` — meta-generator (390 LOC)
- `orpheus/derivations/continuous/fn_method/cases.py` — fn_method registry (220 LOC)
- `orpheus/derivations/continuous/peierls_nystrom/cases.py` — augmented (bc + status)
- `tests/derivations/test_capability_matrices.py` — 3 foundation tests
- `docs/theory/fn_method.rst` — added "Capabilities at a glance" section
- `docs/theory/_fn_method_capability_matrix.inc.rst` — auto-generated
- `docs/theory/_peierls_nystrom_capability_matrix.inc.rst` — auto-generated
- `docs/conf.py` — single `_regenerate_capability_matrices` hook

## Outstanding work (per plan)

- Step 4: `algebra-of-record` skill update — main agent handles in
  parallel (NOT in scope for these 3 commits).
- Steps 5-8: backfill matrices for `singular_eigenfunction`,
  `galerkin_spectral`, `trajectory_resolvent`, `flat_source_cp` —
  deferred to future sessions; estimated ~30 min each per the plan.
