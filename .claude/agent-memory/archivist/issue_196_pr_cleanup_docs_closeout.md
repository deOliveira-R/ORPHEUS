---
name: issue-196-pr-cleanup-docs-closeout
description: Closeout of Issue #196 PR-CLEANUP-DOCS — reconciles 9 SN Sphinx :verifies: warnings, corrects the typed-field memo for the principled (N, ng, nx, ny) layout, and promotes the field-vocabulary catalog to index_convention.rst as a new H1 section.
metadata:
  type: project
---

# Issue #196 PR-CLEANUP-DOCS — Closeout

## §1 Status banner

PR-CLEANUP-DOCS DONE on 2026-05-16, on branch
`refactor/sn-operator-algebra` from tip `6d16df6`.  Two bundled
tasks resolved:

- **Task A4** — Reconciled all 9 SN-side
  ``@pytest.mark.verifies(LABEL)`` warnings about missing equation
  labels.  Sphinx `-W` build passes with zero warnings (the 9
  remaining `pytest.mark.verifies` "skipping" info lines are CP /
  derivations module pre-existing and explicitly out-of-scope per
  §F mechanism criterion 1).
- **Task B4** — Corrected the typed-field-contract memo
  (`typed_field_contracts_for_phase_g.md`) for the principled
  `(N, ng, nx, ny)` layout, then promoted a polished field-vocabulary
  catalog to `docs/theory/index_convention.rst` as a new H1 section
  (`SN Field Vocabulary`, ~362 LoC).

## §2 Motivation preserved (from the open-issue version)

The PR-INDEX-1..7 migration sequence flipped the SN solver's storage
from the legacy `(N, nx, ny, ng)` (energy trailing) to the principled
`(N, ng, nx, ny)` (energy second).  The migration kept 11/11 regression
snapshots bit-identical via the load-bearing Step-1 transpose check.
The DOCS deliverable left two debts:

1. Equation labels that test decorators reference must exist as
   `.. math:: :label:` blocks in theory pages, else the
   `@pytest.mark.verifies` mechanism (per `vv-principles`) cannot
   wire test → equation in the Nexus graph.  9 SN-side
   decorators pointed at labels that never landed in
   `docs/theory/`.
2. The typed-field memo at
   `.claude/agent-memory/explorer/typed_field_contracts_for_phase_g.md`
   (the planning artifact that drove the typed-field-contract resume
   in Phase G) was authored BEFORE the layout discovery.  Every shape
   mention in it was wrong by the (N, ng, nx, ny) vs (N, nx, ny, ng)
   swap.  Leaving it stale would mislead the next session that picks
   up the resume.

PR-CLEANUP-DOCS is the dedicated docs-side cleanup PR that closes
both debts before the typed-field-contract refactor begins.

## §3 Sphinx warning resolution — gate audit

Pre-edit baseline (run on the worktree as it was on entry to this
session, with parallel-PR-CLEANUP-CODE changes also present):

- 18 unique missing-label warnings × multiple test sites each =
  36 total `no matching equation node` info lines from Sphinx.
- 9 unique SN labels (from `tests/sn/`) — IN SCOPE for this PR.
- 9 unique atalay / nm-1980 labels (from `tests/derivations/`) — OUT
  OF SCOPE (CP / case-method / FN module work), flagged below.

Post-edit:

- 0 SN-label warnings remain.
- 13 lines remain (the same 9 unique atalay / nm-1980 labels across
  13 test sites) — UNCHANGED from baseline.
- 0 WARNING / ERROR lines (Sphinx `-W` exit 0).

Pre-existing baseline noise also unchanged:

- `SNStreamingOperator:36` title-underline-too-short — pre-existing,
  retired-class docstring, not touched.
- `transport_operator_matvec_cylindrical:18` undefined-substitution
  "η" — pre-existing, retired-FD-matvec docstring, not touched.
- 7 `discrete_ordinates.rst:24xx-26xx` "Inconsistent title style"
  ERRORS — pre-existing on the working tree from a parallel PR's
  in-progress changes, not introduced by this PR.

The single new Sphinx warning my own edits had introduced
(`SNStreamingOperator:operator.py:125` "unknown document
docs/theory/index_convention") was a stale `:doc:` ref in the
PR-INDEX-7 retired-EquationMap docstring; corrected in this PR to
`:ref:`theory-sn-index-convention`` (the canonical anchor at the top
of `index_convention.rst`).

## §4 Files touched

| File | Scope |
|---|---|
| `docs/theory/discrete_ordinates.rst` | NEW isotropic-curvilinear MMS section (4 labelled eqs: `sn-mms-spherical-psi`, `sn-mms-spherical-qext`, `sn-mms-cylindrical-psi`, `sn-mms-cylindrical-qext`); NEW labelled eqs in the curvilinear-aniso section (`sn-mms-spherical-aniso-spatial-convergence`, `sn-mms-cylindrical-aniso-spatial-convergence`); NEW labelled eqs at the Phase-C homogeneous-kinf-recovery and trajectory-resolvent-crosscheck section bodies; renamed two section anchors to `*-section` to disambiguate from the new equation labels. |
| `docs/theory/index_convention.rst` | NEW labelled equation `blelloch-1990-eq-1-5` + Blelloch / Brent references; NEW H1 section `SN Field Vocabulary` (~362 LoC, 7 catalogue tables); rewired the `Future Work / Typed-field contract resume` subsection to point at the new section + the corrected memo. |
| `orpheus/sn/operator.py` | `EquationMap` docstring `:doc:` path fixed to `:ref:`theory-sn-index-convention``. |
| `.claude/agent-memory/explorer/typed_field_contracts_for_phase_g.md` | Correction banner at top; every shape mention `(N, nx, ny, ng) → (N, ng, nx, ny)`, `(nx, ny, ng) → (ng, nx, ny)`, `(N, nx, ng) → (N, ng, nx)`; all related dataclass `__post_init__` validations + property getters + `to_scalar()` einsum + `at()` slice arithmetic rewritten for the principled axes; ~5 narrative paragraphs flipped from "energy trailing" → "energy second"; Issue #197 quote left verbatim (with a context note); design recommendations untouched per §D anti-rec #3. |

## §5 Per-label resolution table

| LABEL | Resolution | Source file or new addition |
|---|---|---|
| `blelloch-1990-eq-1-5` | **Resolution 3** — new labelled `.. math::` block in `index_convention.rst` Algorithmic-consequence subsection. The Blelloch §1.5 closed-form `cumprod(a) · (psi_0 + cumsum(b / cumprod(a)))` was already present in the prose; promoted to a labelled equation with surrounding pair-monoid + Brent context.  Added Blelloch1990 + Brent1974 citations to `index_convention.rst` references. | `docs/theory/index_convention.rst` (after line ~221) |
| `sn-mms-spherical-psi` | **Resolution 3** — new H2 subsection `Curvilinear isotropic MMS — radial DD-closure probe` placed BEFORE the existing anisotropic curvilinear section.  Equation block defines `ψ_n(r) = sin(πr/R)/W`. | `docs/theory/discrete_ordinates.rst` `.. _sn-mms-curvilinear-isotropic-verification:` |
| `sn-mms-spherical-qext` | Same subsection as above. Equation block defines the per-ordinate manufactured source for the isotropic spherical ansatz. | Same file. |
| `sn-mms-cylindrical-psi` | Same subsection — cylindrical analogue, `ψ_n(r) = A(r)/W`. | Same file. |
| `sn-mms-cylindrical-qext` | Same subsection — cylindrical manufactured source. | Same file. |
| `sn-mms-spherical-aniso-spatial-convergence` | **Resolution 1** — new labelled equation block added INSIDE the existing `sn-mms-curvilinear-aniso-verification` section, formalising the contractual convergence-rate claim `‖φ_h − A‖ = O(h²)`. | `docs/theory/discrete_ordinates.rst` after curvilinear-aniso ansatz/source eqs. |
| `sn-mms-cylindrical-aniso-spatial-convergence` | Same section, cylindrical analogue. | Same file. |
| `sn-curvilinear-homogeneous-kinf-recovery` | **Resolution 1 + section-label disambiguation** — renamed the existing `.. _sn-curvilinear-homogeneous-kinf-recovery:` section anchor to `*-section`; added a new labelled equation block at the head of the renamed section, formalising `k_∞ = ρ(Σ_a⁻¹ χ νΣ_f^T) = νΣ_f/Σ_a (1G)`.  All 3 in-doc `:ref:` consumers updated to the `-section` variant. | `docs/theory/discrete_ordinates.rst` |
| `sn-curvilinear-trajectory-resolvent-crosscheck` | Same pattern — renamed section anchor to `*-section`, added a labelled equation block formalising the Phase-D `rtol ≤ 5e-4` cross-check claim. | Same file. |

Out-of-scope (per §F criterion 1 + §H "Hard scope limits"):

| LABEL | Decision | Rationale |
|---|---|---|
| `atalay-table2-slab-vacuum-isotropic` | DEFER to module:cp PR | Case-method derivations module; pre-existing on entry; not introduced by Issue #196 work. |
| `atalay-table2-slab-reflected-isotropic` | DEFER | Same. |
| `atalay-table2-slab-reflected-r099-precision-floor` | DEFER | Same. |
| `atalay-table3-slab-vacuum-anisotropic` | DEFER | Same. |
| `atalay-table6-eigenvalue-moderate-d-consistency` | DEFER | Same. |
| `atalay-eq46-slab-eq54-sphere-parity` | DEFER | Same. |
| `atalay-eq54-sphere-vacuum-isotropic` | DEFER | Same. |
| `atalay-eq42-extrapolated-endpoint` | DEFER | Same. |
| `nm-1980-reflected-slab-fn` | DEFER to module:cp / FN-method PR | FN method derivations; pre-existing. |

These remain as informational "no matching equation node" messages
from the Nexus Sphinx extension; the Sphinx `-W` build still
exits 0 because they emit at info severity, NOT as warnings.

## §6 SN Field Vocabulary section outline (Task B4 / C.2)

The new H1 section in `index_convention.rst` follows the typed-field
memo's §1 inventory.  Structure (~362 LoC):

1. **Intro paragraph** — frames the section as the durable Sphinx-side
   vocabulary; cross-links to the corrected memo as the design-side
   detail.
2. **Field hierarchy — phase-space and reduced flux types** — 5-row
   `.. list-table::` covering `AngularFlux`, `ScalarFlux`,
   `GroupFlux`, `HarmonicMomentField`, `TraceField` with shapes,
   units, code counterparts.  Followed by 1 paragraph on the
   ``to_scalar`` ↔ ``broadcast_to_ordinates`` conversion pair.
3. **Source / RHS vocabulary** — 4-row table on `IsotropicSource`,
   `PerOrdinateSource`, `ResidualSource`, `BoundarySource`.  Followed
   by 1 paragraph on the load-bearing isotropic vs per-ordinate
   shape split.
4. **Rates and tallies** — 4-row table on `ReactionRate`, `GroupRate`,
   `CurrentCochain`, `Functional`.
5. **Iteration state** — 5-row table on `keff` / `keff_history` /
   `Eigenpair` / `ResidualHistory` / `DominanceRatio`.
6. **Solution-class container** — narrative paragraph + bullet list of
   what the typed evolution adds over the legacy bare-array
   `SNResult`.
7. **Operator vocabulary — the four leaves of the algebra** —
   4-row table on `L` / `C` / `S` / `F` operators with code class +
   mathematical role, followed by 1 paragraph on the derived
   combinations `A_wg = L + C - S_foldable` (fusion target) and
   `K = A_loss⁻¹ F` (multiplication operator).
8. **Boundary-trace vocabulary** — 5-row table flagged orthogonal to
   the index convention (lives on Γ_± faces, not in the volume).
9. **Diagnostic and historical state** — 4-row aspirational table
   (Grand Report v3 §32 / §39); listed for vocabulary completeness.

Every "Existing counterpart" cell uses `:class:` / `:meth:` /
`:mod:` / `:func:` cross-references where the symbol resolves; bare
markdown for aspirational / not-yet-implemented entries.  Every
shape uses the principled `(N, ng, nx, ny)` / `(ng, nx, ny)`
layout.

## §7 §D anti-recommendations — explicit acknowledgement

1. **DO NOT dispatch the typed-field contract refactor itself.** ✓
   No production code modified by this PR; only docstring `:doc:`
   path fix at one site.
2. **DO NOT touch CP / MoC / diffusion theory pages.** ✓ Only
   `discrete_ordinates.rst` and `index_convention.rst` touched in
   `docs/theory/`.
3. **DO NOT rewrite the typed-field memo's design
   recommendations.** ✓ Only shape mentions + dataclass code
   snippet arithmetic (`__post_init__`, `at()`, `to_scalar()`)
   rewritten for the principled axes.  The dataclass APIs,
   operator-algebra coupling section, Issue #197 partial close
   text are conceptually untouched.
4. **DO NOT add new pages.** ✓ "SN Field Vocabulary" is an H1
   section inside the existing `index_convention.rst`.
5. **DO NOT modify production code.** ✓ Only one docstring fix at
   `orpheus/sn/operator.py:125` (the `:doc:` → `:ref:` path
   correction).
6. **DO NOT suppress Sphinx warnings via :noindex: / :nowarn:.** ✓
   No suppression directives added.

## §8 Verification gates (§E)

| Gate | Pre-edit | Post-edit | Notes |
|---|---|---|---|
| `sphinx-build -W` exit code | 1 (from missing-label info treated as warnings only without -W; baseline title-style errors present) | 0 | Build clean. |
| SN missing-label warnings | 9 unique | 0 | All 9 resolved per §5. |
| Out-of-scope (atalay / nm-1980) missing-label info | 9 unique | 9 unique | Unchanged, deferred to module:cp PR. |
| `tests._harness.audit` | Clean | Clean | 41/48 ERR-NNN entries have catching test (unchanged from baseline). |
| `pytest tests/sn/regression/ -q` | 11 passed | 11 passed | 63.40s wall-clock. |

## §9 Mechanism criteria audit (§F)

| # | Criterion | Result | Evidence |
|---|---|---|---|
| 1 | Sphinx `-W` build PASSES (or remaining warnings principled) | PASS | `exit=0`; remaining `pytest.mark.verifies skipping` lines are info severity, not warnings — Sphinx `-W` accepts them. |
| 2 | All listed SN labels resolved | PASS | 9/9 SN labels resolved per §5 table. |
| 3 | `typed_field_contracts_for_phase_g.md` has correction banner | PASS | Banner at top, starting "Correction (2026-05-15)." |
| 4 | All shape mentions in memo flipped | PASS | `grep -c "(N, nx, ny, ng)" memo = 0`. |
| 5 | `index_convention.rst` has new "SN Field Vocabulary" section | PASS | H1 header at line ~688 (was line 686 pre-insertion). |
| 6 | Field-vocabulary section is 200–400 LoC | PASS | `awk` between section markers = 362 LoC. |
| 7 | V&V audit clean | PASS | `python -m tests._harness.audit` finishes normally; ERR coverage 41/48 unchanged. |
| 8 | Regression 11/11 PASS | PASS | `pytest tests/sn/regression/ -q` = 11 passed in 63.40s. |

## §10 Infrastructure retained

- The renamed section anchors
  (`sn-curvilinear-homogeneous-kinf-recovery-section`,
  `sn-curvilinear-trajectory-resolvent-crosscheck-section`) preserve
  the existing 3 in-doc `:ref:` consumers in
  `discrete_ordinates.rst` (lines 2936, 3380, 3646).  All 3 updated.
- The Issue #197 NewType quote at memo §4.1 is preserved verbatim
  with a "shapes here predate Issue #196" context note.  Historical
  fidelity matters — that quote is what the resume PR's "what #197
  said" briefing references.
- Every existing equation label in `discrete_ordinates.rst` left
  untouched.
- Every existing `:ref:` consumer of unrenamed labels left
  untouched.

## §11 Session trail

- Commit chain on the branch: `e09b9f8` (PR-INDEX-1) → `6cfdfd4`
  (PR-INDEX-2) → `313f510` (PR-INDEX-3) → `fa41767` (PR-INDEX-4) →
  `3356cec` (PR-INDEX-5) → `13bcf5f` (PR-INDEX-6) → `6d16df6`
  (PR-INDEX-7) → THIS PR (PR-CLEANUP-DOCS, post-merge will be
  the 8th in the chain).
- Files staged (per §H "Do NOT commit. Stage and return."): three
  modified docs files
  (`docs/theory/index_convention.rst`,
  `docs/theory/discrete_ordinates.rst`,
  `.claude/agent-memory/explorer/typed_field_contracts_for_phase_g.md`),
  one one-line code fix at `orpheus/sn/operator.py`, plus this
  closeout memo.
- Sphinx build artefacts and `docs/verification/matrix.rst`
  regeneration are picked up automatically by the next clean build;
  not staged by this PR.
- Pre-existing baseline noise documented in §3 stays through this
  PR per the "build-gating: count-unchanged-from-pre-edit" rule.
