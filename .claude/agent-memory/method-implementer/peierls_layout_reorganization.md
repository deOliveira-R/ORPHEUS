---
name: Peierls layout reorganization
description: Split orpheus.derivations.continuous.peierls into peierls_nystrom + peierls_greens_function on 2026-05-02 (feature/peierls-greens-cylinder). Mechanical refactor with zero functional change.
type: project
---

# Peierls layout reorganization — 2026-05-02

Branch: `feature/peierls-greens-cylinder`
Pre-refactor head: `53efba3 docs(agent-memory): Phase-2 cylinder Variant α unification closeout`
Post-refactor commit: `3dfd76e refactor(layout): split peierls/ into peierls_nystrom/ + peierls_greens_function/`

## What

Split the single `orpheus.derivations.continuous.peierls` package into
two named-by-discretization siblings:

- **`peierls_nystrom/`** — scalar-flux integral equation discretized
  via Nyström collocation. Hosts `geometry.py` (god-module),
  `cases.py`, `cylinder.py`/`sphere.py`/`slab.py` facades,
  `ps1982_reference.py`, `reference.py`, and the Nyström-side
  specular origins (`r_matrix`, `slab`, `continuous_mu`).
- **`peierls_greens_function/`** — angle-resolved Variant α Green's
  function on the shared rank-1 resolvent. Hosts `greens_function.py`
  (sphere prod), `greens_function_cylinder.py` (cylinder prod),
  `variant_alpha_core.py` (shared resolvent + closure), and the
  Variant α SymPy origins under `origins/specular/`.

The previously-empty `greens_function/` placeholder was promoted into
this new home (it had been a placeholder waiting precisely for this
content).

## Why

Single-namespace ambiguity: both Nyström and Variant α share Peierls
ancestry but discretize *different* things — Nyström discretizes the
scalar-flux integral equation; Variant α discretizes the angle-
resolved Green's function and reduces analytically via rank-1.
Co-locating them under `peierls/` blurred the structural intent and
would have made Phase-3 (slab Variant α, 2-surface BIE frame) a
moving target. The split now lets Phase-3 land cleanly into
`peierls_greens_function/` without further reorganization.

## How (single-commit form)

The package directory rename breaks all 100+ `from orpheus.derivations.continuous.peierls...` imports until the substitution lands. A 2-commit split (moves first, imports second) would leave the intermediate commit in a non-building state, breaking the
"main is always green" invariant. Single atomic commit is the pragmatic choice.

### Mechanics

1. `git mv` for ALL renames and moves (history preserved):
   - `peierls/` → `peierls_nystrom/`
   - `greens_function/` → `peierls_greens_function/`
   - `peierls/{greens_function,greens_function_cylinder,variant_alpha_core}.py` →
     `peierls_greens_function/...`
   - `peierls/origins/specular/{greens_function,greens_function_cylinder}.py` →
     `peierls_greens_function/origins/specular/...`

2. New `__init__.py` files for the new `peierls_greens_function/origins/` and
   `peierls_greens_function/origins/specular/` (re-exports the 7
   `derive_*` Variant α functions for sphere + cylinder).

3. Trimmed `peierls_nystrom/origins/specular/__init__.py` to drop the
   moved Variant α re-exports (kept only Nyström-side: r_matrix, slab,
   continuous_mu).

4. Updated package docstrings on `peierls_nystrom/__init__.py` and
   `peierls_greens_function/__init__.py` to cross-reference each other
   (the Peierls ancestry is shared; the discretization is the split).

5. Import substitution (Python script `/tmp/refactor_imports.py`) in
   3 ordered passes to avoid swallowing the Variant α paths under the
   catch-all rule:
   - Pass 1: longest Variant α paths → `peierls_greens_function`.
   - Pass 2: in the 2 greens-symbolic test files only, redirect
     `peierls.origins.specular` → `peierls_greens_function.origins.specular`
     (those tests import the V_α derive functions).
   - Pass 3 (catch-all, regex with negative lookahead for `_`):
     remaining `orpheus.derivations.continuous.peierls.X` →
     `orpheus.derivations.continuous.peierls_nystrom.X`.

## Files moved (counts by category)

- 13 modules renamed via directory rename (peierls/ → peierls_nystrom/).
- 5 modules relocated via explicit git mv (Variant α prod + SymPy +
  shared resolvent → peierls_greens_function/).
- 1 placeholder package init removed; replaced by the new package
  init under the renamed directory.
- 2 brand-new `__init__.py` files created (`peierls_greens_function/origins/`
  and `peierls_greens_function/origins/specular/`).

## Files modified for imports

67 files modified by the import-substitution script:

- 4 in `orpheus/derivations/common/` (docstring cross-refs).
- 7 in `orpheus/derivations/continuous/` (cross-refs from sibling
  packages; package init docstrings; the two greens-prod modules and
  one origins module needed `peierls.variant_alpha_core` →
  `peierls_greens_function.variant_alpha_core`).
- 5 in `orpheus/derivations/continuous/peierls_nystrom/` (internal
  cross-refs).
- 39 in `tests/` (peierls + cp, every importer caught by the regex).
- 1 in `tools/verification/` (`generate_peierls_matrix.py`).
- 6 in `docs/theory/` (every `:mod:`/`:func:` cross-ref).
- `docs/verification/matrix.rst` auto-regenerated by Sphinx.

Representative diff hunk (a typical test-file change):

```diff
-from orpheus.derivations.continuous.peierls.greens_function import (
+from orpheus.derivations.continuous.peierls_greens_function.greens_function import (
     solve_greens_function_sphere,
 )
-from orpheus.derivations.continuous.peierls.geometry import (
+from orpheus.derivations.continuous.peierls_nystrom.geometry import (
     solve_peierls_1g,
 )
```

Total: 82 files changed in the commit, 558 insertions / 493 deletions
(most of the line delta is the path-string lengthening from
`peierls` to `peierls_nystrom` / `peierls_greens_function`, with the
docstring updates and verification-matrix regeneration accounting for
the rest).

## Test counts (zero regression)

| Suite                                | Pre   | Post          |
| ------------------------------------ | ----- | ------------- |
| Primary 9-file Variant α gate        | 65/65 | 65/65         |
| `tests/cp/test_*peierls*.py`         | 24/24 | 24/24         |
| Broader `tests/derivations/test_peierls_*.py` | 625 collected, exit 0 | 625 collected, exit 0 |
| `sphinx-build -W`                    | clean | clean         |

Accuracy floors preserved (the primary suite re-runs at the same
tolerances the Phase-1 / Phase-2 closeouts certified: 1e-15 at α=1,
~8e-8 self-consistency at α=0 vacuum, ≤1e-4 PS-1982 cross-verification
across 4 configurations).

Sphinx-side: no NEW warnings beyond the 4 pre-existing
`pytest.mark.verifies('sn-mms-*')` info-level skips that were already
present at baseline. The verification matrix's auto-regeneration
discovered new entries for `test_peierls_variant_alpha_core` and the
`peierls-greens-cylinder-*` equation labels (these were collectable
before the refactor but the matrix had been stale; rebuild brought
them in — a side benefit of the refactor running Sphinx).

## Subtle breakages encountered

None. The only nontrivial subtlety was the substitution-order
discipline (Variant α paths must be substituted FIRST so they don't
get eaten by the `peierls.X → peierls_nystrom.X` catch-all). The
Python script's 3-pass design handles this; the regex used for the
catch-all has a negative lookahead `(?![_a-zA-Z0-9])` so
`peierls_nystrom` and `peierls_greens_function` are NOT re-touched.

The two greens-symbolic test files (`test_peierls_greens_function_symbolic.py`
and `test_peierls_greens_function_cylinder_symbolic.py`) deserved
explicit handling because they imported from `peierls.origins.specular`
expecting the V_α `derive_*` symbols, which had moved out from under
that path. Those imports are now redirected to the new
`peierls_greens_function.origins.specular`.

No back-compat shims needed. No circular imports surfaced.

## Verdict on Phase-3 fitness

The new layout disambiguates Nyström vs Variant α cleanly enough to
support Phase-3 (slab Variant α, the 2-surface BIE frame from
`cross-domain-attacker/variant_alpha_2surface_bie_frame.md`)
without further reorganization.

Phase-3 will land:
- Branch-2 production code at
  `orpheus/derivations/continuous/peierls_greens_function/greens_function_slab.py`.
- Branch-1 SymPy at
  `orpheus/derivations/continuous/peierls_greens_function/origins/specular/greens_function_slab.py`.
- Re-exports added to
  `peierls_greens_function/origins/specular/__init__.py`.
- Tests at `tests/derivations/test_peierls_greens_function_slab_symbolic.py` /
  `..._slab_solver.py` / `..._slab_xverif.py`.
- Sphinx narrative at `docs/theory/peierls_greens.rst` (extended).

The shared `variant_alpha_core.py` rank-1 resolvent generalizes to
rank-2 cleanly (per the BIE frame attack); rank-1 sphere/cylinder
will keep using their existing entry points, rank-2 slab/annulus/
hollow-sphere will mount on a new `variant_alpha_core_rank2.py`
sibling. The Nyström side is unaffected by Phase-3 — it stays in
`peierls_nystrom/`.

## Pointers

- Commit: `3dfd76e refactor(layout): split peierls/ into peierls_nystrom/ + peierls_greens_function/`
- Acceptance gate evidence: log of 65/65 pass at same accuracy floors
  preserved in the commit body.
- Cross-domain-attacker memory:
  `.claude/agent-memory/cross-domain-attacker/variant_alpha_2surface_bie_frame.md`
  (rank-2 BIE frame, ready for Phase-3 deployment).
- Related closeouts:
  `cylinder_variant_alpha_phase1.md`, `cylinder_variant_alpha_phase2_unification.md`.
