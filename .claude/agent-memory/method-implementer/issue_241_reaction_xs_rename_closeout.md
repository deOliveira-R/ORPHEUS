---
name: issue-241-reaction-xs-rename-closeout
description: #241 model-agnostic PARAMETER rename across the DiscretizationScheme interface — sigt_cells/sig_t→reaction_xs, cell_average_weight→face_blend_weight. Pure rename, bit-identical, source-hash sentinel re-hashed.
metadata:
  type: project
---

# #241 — model-agnostic DiscretizationScheme parameter rename (closeout)

Branch `refactor/sn-foundation-cleanup` (2026-06-18; NOT committed — main
agent commits). Pure identifier rename, ZERO numerical change, bit-identical.
Authoritative verdict: `.claude/agent-memory/cross-domain-attacker/discretization_scheme_naming_signal.md`.

## What was renamed (the SCOPE decision — load-bearing)

The brief's AUDIT clause was the binding scope: "ZERO `sigt_cells` /
`cell_average_weight` in `orpheus/`; `sig_t` survives ONLY as caller-side
variables". `total_xs` was NOT in the zero-clause and NOT on the brief's
surface, so it was LEFT (it is the per-cell scalar `update`/`residual` form,
not the model-agnostic coefficient layer).

1. **`sigt_cells` → `reaction_xs`** (ALL occurrences — it is a scheme-specific
   param name, every occurrence IS the reaction-rate concept). Renamed in:
   `scheme.py` (`cell_kernel_batch`/`residual_kernel_batch` base sigs +
   docstrings), `diamond.py` (`cell_kernel_batch`/`residual_kernel_batch` +
   the private `_cartesian_streaming_diagonal` helper they share),
   `linear_discontinuous.py` (`cell_kernel_batch`/`residual_kernel_batch` +
   private `_ubld_system`), `sweep_graph.py` (`_CellSolve.cell` /
   `_CellResidual.cell` params + bodies + the `is_moment_valued_by_rank`
   comment), `moment_layout.py` (one docstring example reference).
2. **`sig_t` PARAMETER (scan-coefficient/closure methods only) → `reaction_xs`**:
   `affine_scan_coefficients` (scheme/diamond/LD), `cartesian_scan_coefficients`
   (scheme/diamond), `moment_scan_closure` (LD — NOT in the brief's "4 methods"
   list but carries `sig_t` as a scheme kwarg, so renamed). Param sig + body +
   docstring.
3. **`cell_average_weight` → `face_blend_weight`**: the `CollisionCache` FIELD
   (`sweep_cache.py:375`) + its 2 consumers (`loss_representation.py` 2948/3199
   `coll.cell_average_weight`→`coll.face_blend_weight`) + the local unpack var
   in `from_geometry`. The math symbol `w` is KEPT (only the long handle moved).

## The sig_t-PARAM-vs-VARIABLE disambiguation (the brief's CRITICAL constraint)

`sig_t` is the CODEBASE-WIDE Σ_t variable. ONLY the `sig_t` *keyword at the
scheme-method call boundary* changed; caller-side `sig_t` VARIABLES kept their
name:
- **RENAMED (keyword only, value expr untouched):** at the scheme scan/kernel
  call sites — `sweep_cache.py:475` `reaction_xs=sig_t_chain`;
  `loss_representation.py` 1565 (`cartesian_scan_coefficients`
  `reaction_xs=sig_t[None,:,:,j]`), 1686/2354 (`residual_kernel_batch`
  `reaction_xs=sig_t[...]`/`reaction_xs=sigma_gx[...]`), 3024
  (`moment_scan_closure` `reaction_xs=sig_t_chain`); `sweep_graph.py` 762/814
  (`level_op.cell(reaction_xs=sig_t[...])`).
- **LEFT (caller-side Σ_t holders):** `walk_full`/`walk_windowed`'s own `sig_t`
  PARAM (`sweep_graph.py` 731/776); `CollisionCache.from_geometry`'s `sig_t`
  param + its `sig_t_chain` local (`sweep_cache.py` 407/460); every
  `sig_t=` keyword in `loss_representation.py` (768/1057/1131/1311/1380/3367)
  that feeds `_ApplyOperands`/`_SolveOperands` dataclasses or
  `walk_full`/`walk_windowed`; ALL `sig_t=` across `tests/derivations`,
  `orpheus/derivations`, `data/macro_xs`, CellXS/Mixture constructors — these
  are NOT scheme-method calls.

## KEPT unchanged (memo's explicit KEEP list — confirmed NOT touched)

`streaming` / `s_axes` (the one genuine SN/CFD frame conflict — diffusion has
no advective μ-coefficient); `a` / `a_attenuation`; `inverse_denom`; the
reconstruction ops `source_emission` / `cell_average` /
`outgoing_face_from_average`; the `DiscretizationScheme` class name (deferred
high-churn hub rename); and `total_xs` (the per-cell scalar `update`/`residual`
form, intentionally out of scope per the audit clause). Verified by grep: none
of these appear in any of my diffs.

## Docstrings updated to past-tense "DONE"

`scheme.py` class docstring §"Generic advection–reaction interface": the
"reaction-rate is an EXPLICIT argument (`sig_t`/`total_xs`/`sigt_cells`)" line
now reads "`reaction_xs` on the coefficient + kernel methods; `total_xs` on the
per-cell update/residual"; the "DEFERRED to #241" note rewritten to "DONE
(#241)" — phrased WITHOUT the literal old tokens so a future grep audit of
`orpheus/` stays clean of `sigt_cells`/`cell_average_weight`.

## THE ONE SURPRISE — the source-hash sentinel (re-hashed, documented)

`tests/sn/sweep/core/test_cell_kernel_batch.py::TestKernelSourceOfRecord`
sha256-pins the SOURCE TEXT of `DiamondDifference.cell_kernel_batch` /
`residual_kernel_batch` as "the FP reduction tree of record". The rename
changed the body's source TEXT (param `sigt_cells`→`reaction_xs` + the
`_cartesian_streaming_diagonal(reaction_xs,...)` call), so the sentinel
RED-flagged (2 failures). This is the guard working as designed — an
identifier rename touches no operand/operation/fold-order, so the FP tree is
bit-for-bit preserved. ACTION (per the test's own instruction + vv-principles
bit-identity discipline): updated both EXPECTED hashes to the new source text
+ added a re-hash comment documenting #241 (rename-only, FP tree PRESERVED).
New hashes: cell_kernel_batch
`0e3195534ec7fa4257c5ad6ae294b6d531fb9282997ef0e348153c85e5cef03c`;
residual_kernel_batch
`96f5a7d8fff90f153fb9f83a63f0a05592e7d84be1e6197df72a318c9adc5d7b`.

LESSON (worth a coding-elegance note): a source-text sentinel (sha256 of
`inspect.getsource`) is bit-identity-load-bearing but NOT rename-transparent —
ANY identifier rename inside the pinned body reddens it even when the FP tree
is provably untouched. When a rename crosses such a sentinel, the re-hash is
MANDATORY and must cite the proof (the DD strict regression gate stayed
green at the SAME within-tol ULP values), not just assert it.

## Test-file keyword updates (gate sets include them)

The renamed methods are keyword-only (`*`), so a stale `sigt_cells=`/`sig_t=`
keyword raises TypeError (loud). Updated:
`test_cell_kernel_batch.py` (sed `sigt_cells`→`reaction_xs`, 26 occ — both
keyword + local var, test-local),
`test_sweep_graph_nd_admission.py` (2 kwargs; local `sigt` Σ_t scalar KEPT),
`test_ld_ubld_primitive.py` (`cell_kernel_batch reaction_xs=` +
`affine_scan_coefficients reaction_xs=`),
`test_linear_discontinuous.py` (sed `sigt_cells=`→`reaction_xs=` 7 kwargs;
local `sigt` Σ_t var KEPT + `affine_scan_coefficients reaction_xs=`),
`test_scheme_reaction_rate_contract.py` (2 `affine_scan_coefficients
reaction_xs=`; the `_sig_t()` Σ_t-field HELPER name KEPT — caller-side).

## GATES (all as expected)

- **`tests/sn/spatial` + `tests/sn/sweep/core` + `tests/sn/sweep/cartesian_2d`
  + LD MMS (2d+slab)**: 596 passed (the 2 sentinel-hash reds were the rename
  artefact; FIXED by the re-hash → re-ran the 2 sentinels = 2 passed; touched-
  test re-run 98 passed).
- **`tests/sn/operators`**: 505 passed, 7 FAILED = the documented baseline
  (5 sphere snapshot #250 — 3 `test_vacuum_bulk_bit_identical_1d[*-SPH]` +
  2 `test_sphere_{1g,2g}_apply_bit_identical`; 2 `mu_y` #232 — `Face 'ymin'`).
  ZERO new failures.
- **DD strict regression** (`tests/sn/regression -W
  error::...DriftWarning`): 13 passed, 13 within-tol DriftWarnings (harness
  reports within-tol drifts as warnings, NOT escalated despite `-W error`) —
  IDENTICAL pre-existing baseline ULP values (e.g. `2d_2g_p1_aniso` scalar_flux
  6920 ULP / 9.81e-13, the known #240-D5a value). NO snapshot regenerated.
  THIS is the bit-identity proof: the FP tree is preserved.

## Deliverable-manifest note

This is a PURE RENAME refactor (issue #241), NOT a new reference-solver build —
the algebra-of-record deliverable manifest (Branch-1 SymPy module, foundation
test gate, L1 cross-check, Sphinx stub, archivist dispatch) does NOT apply.
No new math, no new claim, no new equation label. The scheme docstrings carry
the model-agnostic vocabulary rationale already (advection–reaction frame,
cited at scheme.py §"Generic advection–reaction interface"); no Sphinx stub
needed. Did NOT commit (main agent commits per brief). Did NOT run the full
Sphinx build (main agent batches it). NO git checkout/restore/stash/reset (L28).
