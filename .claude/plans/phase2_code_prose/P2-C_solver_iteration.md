# P2-C — `orpheus/sn/solver.py` + `orpheus/numerics/iteration.py` code-prose rebalance

Phase 2 of #231, batch P2-C. Branch `docs/sn-doc-architecture`. Pilot
calibration: `P2-A_scattering.md` (lesson L-033). Written INCREMENTALLY.

Pre-edit (measured, `tokenize` + `ast`):

| file | lines | docstring | comment | prose% | code-tokens |
|---|---|---|---|---|---|
| `orpheus/sn/solver.py` | 3137 | 1164 | 621 | 56% | 6743 |
| `orpheus/numerics/iteration.py` | 1334 | 715 | 183 | 67% | 2168 |

## Headline finding (consistent with the pilot)

**ZERO MOVED. The operator-algebra + SN book is exhaustively complete.**
Every teaching/history concept these two DRIVER files carry has a book
home — verified by grepping the landing chapters:

- k-estimator R7 / #291 leakage / net-removal / balance identity →
  `methods/sn/solver.rst` §§ (`:eq:sn-keff-update`, `:eq:sn-leakage-functional`).
- ERR-026 / ERR-058 / Carlson closure-seed history → `curvilinear_numerics.rst`,
  `verification.rst`, `history.rst`.
- RULING P1 / reflective `B_a+B_b` coupling → `boundary_conditions.rst`,
  `slab_multigroup.rst`, `history.rst`.
- #282 lag-death / Mode-12 / convergence certificate → `solver.rst`,
  `loss_representation.rst`, `curvilinear_numerics.rst`.
- SI ≡ Krylov twin / `k_inf=1.875` → `verification.rst` (1182–3494).

**These are CONTRACT files, not teaching files.** Unlike the pilot's
scattering.py (73% derivation-prose, cut −36%), solver.py/iteration.py
docstrings are mostly driver/interface CONTRACT: the keff estimators,
convergence math, the #282 threat model, the ERR-053 restart trap, the
reflective-coupling SSOT. The cut surface is therefore the standalone
**retirement tombstones** + the **SPLIT-docstring HISTORY tails** +
pure-narration inline comments — NOT the CONTRACT interfaces. Expect a
much smaller cut ratio than the pilot; that is the correct outcome for a
driver file (reported honestly, like the pilot reported ZERO MOVED).

## Structural cut-safety facts

- **`orpheus.numerics.iteration` is NOT automodule'd** (grep clean) → its
  docstrings are invisible to the Sphinx build; no link breaks on a cut
  (same as the pilot's scattering.py). NO sphinx gate needed for it.
- **`orpheus.sn.solver` IS automodule'd** at `docs/api/discrete_ordinates.rst:25`
  (`:members: :undoc-members: :show-inheritance: :noindex:`) → its
  docstrings ARE rendered. **Sphinx `-E` gate REQUIRED for solver.py**
  (cutting must not leave malformed RST / drop-into-a-warning). The
  `:noindex:` means solver.py's own symbols are not xref targets, but
  refs FROM its docstrings to other (registered) targets still resolve.
- **NO `:label:` / `vv-status` / `verifies`/`catches` marker inside either
  `.py`** (grep clean) → cutting docstring math orphans no `:eq:` target
  and no test marker points at these docstrings (targets are RST labels).
- Pre-edit `-E` baseline build: **0 warnings / 0 errors / 0 critical**
  (cleaner than AGENT.md's expected 1). Acceptance gate: post-edit `-E`
  set unchanged (still 0).

## Batch-special results

- **SPECIAL 1 (row-8 dual-A bridge) — SATISFIED, no fix needed.**
  `notation.rst` crosswalk row 8 promises the bridge is "stated at both
  ends (the module heads of `orpheus.numerics.iteration` and
  `orpheus.numerics.operator`, and the solver page)". `iteration.py`'s
  module head (lines 1–30) FULLY states it: line 7/12 pose
  `(A − Σᵢ gᵢ)ψ = q_ext` / `= (1/k)Fψ`; lines 15–16 declare A the
  INVERTIBLE resolvent operand; lines 18–25 give the SN binding
  (`A = L + C`, streaming `L=Ω·∇` + collision `C=Σ_t·`, gains `(S, B)`,
  the honest `L+C−S−B`), the matrix-method shorter-gain-tuple case, and
  "fission `F` is never a gain … the outer loop scales it by `1/k`". The
  spelling matches the row-8 record verbatim in substance. This module
  head is **CONTRACT — untouched** by the rebalance.
- **SPECIAL 2 (dated `(A−S−F)ψ=q_ext` posings) — CONFORMANT, both files.**
  Grep for `A[-−]S[-−]F` / `L[-−]S[-−]F` / `(A - S - F)` / `(L - S - F)`
  across BOTH files: **zero hits**. The catalog's old deferred-list line
  numbers (~1/7/32/861) all now use the variadic `(A − Σᵢ gᵢ)` record.
  Nothing to harmonize.

## MOVED candidates

**NONE.** (See headline — the book is complete.)

---

## Adjudication table — `orpheus/numerics/iteration.py`

One row per docstring ≥10 ln / comment ≥5 ln (live line ranges).

| lines | symbol | verdict | destination / record | content id |
|---|---|---|---|---|
| 1–138 | MODULE head | **CONTRACT (untouched)** | notation.rst §notation-crosswalk row 8 (oracle) | dual-A bridge (SPECIAL 1); posing + SN binding + Carlson-seed + forward refs |
| 165–191 | `SupportsSeededApply` | CONTRACT (kept) | — | seeded-apply static contract; #285 structural conformance |
| 196–203 | comment: KeffEstimator injection retired | **HISTORY → trim** | KEigenvalue.Notes (1077–1086, same file) owns the constraint | #259 retirement narration on the `Preconditioner` alias |
| 206–221 | `_SeededExactApply` | CONTRACT (kept) | — | two-kinds inverse split; accept-and-drop seed (a constraint) |
| 231–273 | `_wrap_delegate_member`/`_seeded_inverse` | CONTRACT (kept) | — | the ONE inverse→driver adaptation; #285 structural dispatch |
| 276–292 | comment: Ravellable protocol | CONTRACT (kept) | — | L1↛L2 decoupling constraint; to_flat/from_flat duck-type |
| 295–433 | ravel helpers + `_DisplacementLeaf`/`_flux_displacement_leaf` | CONTRACT (kept) | — | protocol faces, diagnostics, systems[0] recurse convention |
| 436–441 | comment: module-level estimator fns retired | **HISTORY → cut** | git log + KEigenvalue.Notes | pure retirement archaeology (#259 fold) |
| 449–565 | `SourceIteration` class | CONTRACT (kept) | — | the ρ-honest residual contract; #282 blindness; exact-M identity |
| 567–599 | `SourceIteration.__init__` | CONTRACT (kept) | — | eager apply-guards; pre-inverted step operator |
| 601–708 | `SourceIteration.solve` + body | CONTRACT (kept) | — | ravellable protocol; the stop logic; seed threading |
| 716–819 | `KrylovAcceleration` class | CONTRACT (kept, 1 trim) | Adams&Larsen 2002 (kept) | GMRES contract; preconditioner≠inverse-step (constraint). Trim only the dated `(R-1 Step B, 2026-05-19)` provenance |
| 821–1003 | `KrylovAcceleration.__init__`/`solve` + body | CONTRACT (kept) | — | GMRES setup; ERR-053 non-convergence surface; exact-breakdown invariant |
| 1011–1098 | `KEigenvalue` class | CONTRACT (kept) | — | k-posing realization; the Notes owns the #259 estimator constraint |
| 1100–1334 | `KEigenvalue.__init__` … `.solve` | CONTRACT (kept) | — | construct-time A-invertible; the estimators; boundary methods |

**iteration.py verdict: near-total KEEP.** Cuts = the two pure-narration
retirement comments (196–203 cut, 436–441 cut) + one dated-provenance
trim (KrylovAcceleration §, kept the category-mistake CONSTRAINT). This
is the correct light touch for a dense CONTRACT driver-primitive file.

---

## Adjudication table — `orpheus/sn/solver.py`

Pre-edit line ranges. CONTRACT interfaces (`_apply_default_bcs`,
`_as_sn_mesh`, `_system_a_residual`, `evaluate_residual` core,
`_certify_within_group_exit`, `_within_group_krylov`,
`_select_si_resolvent`, `_within_group_si`, all `compute_*` estimators,
`_boundary_leakage_rate`, `_face_area_of`, `_reflect_outflow_into_inflow`,
`solve_sn`/`solve_sn_fixed_source` bodies, `_build_fixed_source_rhs`,
`_lift_external_source_to_moments`, the delegators) are KEPT and not
tabulated individually. Only adjudicated CUTS/TRIMS below.

| pre-edit lines | block | verdict | destination / record | content id |
|---|---|---|---|---|
| 166–170 | comment: PR-TYPED-5 SNResult retired | **HISTORY → cut** | git; Solution.is_eigenvalue/is_fixed_source own the discriminator | the `from .solution import` line kept |
| 174–175 | comment: `_zero_within_group_fission` retired | HISTORY → cut (constraint folded) | — | within-group-fission-is-`q_ext` kept, folded into the trimmed block below |
| 178–189 | comment: `_within_group_triple`/`_lagged_gains` retired | **HISTORY → trim** | git | kept the SSOT pointer (build_within_group_system → WithinGroupSystem record; composition contract on its docstring); cut what the retired helpers WERE |
| 504–512 | comment: `_GaussSeidelResolvent` tombstone | **HISTORY → cut** | live design on `_select_si_resolvent` (the G-S arm) | reified-ScheduledInvertibleOperator design already there |
| 515–522 | comment: `_MomentWindowedResolvent` tombstone | **HISTORY → cut** | live design on `_maybe_window` | `P @ A.inverse()` design already there |
| 949–969 | `__init__` PR-TYPED-0/1/2 comments | **CONTRACT (trim)** | — | kept the two single-source-of-truth constraints (mat_xs; ng raises); cut tags + the seven-attrs-collapse + thin-props-retired narration |
| 1068–1078 | comment: PR-TYPED-2 read-through shims retired | **HISTORY → cut** | mat_xs.* accessors self-doc on MaterialXSField | accessor-list archaeology |
| 1189, 1221 | inline "PR-TYPED-2: consumes mat_xs (no shim)" | HISTORY → cut (2×) | — | redundant with the code (reads `self.mat_xs.*`); shape constraints above kept |
| 268–274 | `evaluate_residual` box-7 tail | **CONTRACT (trim)** | verification/loss_representation (Mode-12 kept in core) | cut "#208 box-7 consumer / mint previously unconsumed"; kept diagnostic role + `#2` DSA live pointer + not-in-convergence-path |
| 1466–1534 | `_solve_source_iteration` docstring Scope+Verified | **CONTRACT (trim) + HISTORY tail cut** | verification.rst (k_inf values + ERR-026/058 history) | kept twin-of-`_solve_krylov` + WDD sweep + seed threading + test pointer; cut R-1 Step E date, Wave-O date, the B1'' retirement ¶, k_inf literals |
| 1704–1740 | `_solve_krylov` docstring Scope | **CONTRACT (trim)** | — | cut R-1 Step D date + D-H.2-C4d/S6.3/S6.9 codenames + B1''; kept 2-D-supported + 4-face + matvec walk |
| 1867–1873 | comment: `_make_sweep_preconditioner` tombstone | **HISTORY → cut** | git | zero-callers archaeology |
| 1875–1882 | comment: source-helper delegators | **CONTRACT (trim) + FIX** | — | kept the delegation-surface constraint; trimmed Wave-D tag; **FIXED stale path** `orpheus/sn/scattering.py` → the class ref (live path is transport/operators/scattering.py) |
| 1934–1957 | comment: P1.7 `_build_rhs_*` tombstone (24 ln) | **HISTORY → cut** | git + the named pinning test | packed-1D RHS-builder archaeology |
| 2136–2139 | `solve_sn` Returns PR-TYPED-5 | HISTORY → trim | — | cut "legacy SNResult retired"; kept the Solution-type description |
| 2626–2640 | `solve_sn_fixed_source` inner_solver "History:" | **HISTORY tail → trim** | verification/history (Notes already points there) | cut the dated Phase-D→ERR-058 evolution; kept current-behavior + ERR-058/053 refs + Mode-9 cross-check |
| 2675–2697 | body "Issue #168 status" (23 ln) | **HISTORY → cut to 1-line** | git; the docstring | Phase A/B/C/D status archaeology → 1-line current-behavior constraint |
| 2741–2770 | `_solve_fixed_source_si` header + "loop RETIRED" | **CONTRACT (trim)** | — | cut Phase-1 date + "previous hand-rolled loop RETIRED, rebuilt via _add_*"; kept the coupling-gain mechanism + the `_reflect_outflow_into_inflow` helper-survival constraint |
| 2953–3014 | `_solve_fixed_source_krylov` header + RETIRED list + Scope | **HISTORY → cut + CONTRACT (trim)** | git (dated retirement list w/ commit hashes) | cut G1 date + the 6-item RETIRED-machinery list + the B1'' Scope ¶; kept twin + verified `q/Σ_t` + test pointer |

### solver.py verdict counts (adjudicated cut/trim blocks)

- **HISTORY → cut (standalone tombstones + status blocks)**: 8
  (PR-TYPED-5 SNResult, `_zero_within_group_fission`, `_GaussSeidelResolvent`,
  `_MomentWindowedResolvent`, PR-TYPED-2 shims, `_make_sweep_preconditioner`,
  P1.7 `_build_rhs_*`, Issue #168 status) + 2 inline PR-TYPED-2 tags.
- **CONTRACT (trim: keep the interface + constraints, cut embedded
  HISTORY/date tail)**: 9 (the `__init__` PR-TYPED, `evaluate_residual`
  tail, both inner-solver Scopes, both fixed-source-driver headers/tails,
  the inner_solver History, the source-helper comment, the Returns tail).
- **1 FIX (Cardinal Rule 1)**: stale path `orpheus/sn/scattering.py` →
  the class ref (the module lives at `transport/operators/scattering.py`).
- **MOVED**: 0. **TWIN-to-new-pointer**: 0 (the book already carries
  every teaching concept; the CONTRACT docstrings already carry their
  cross-refs, so cuts are subtractive, not pointer-adding).

**Deliberately NOT cut (scope ceiling):** the 3 remaining
`PR-TYPED-5: build typed Solution at the boundary` inline comments
(solve_sn 2165 / SI 2724 / Krylov 2945) — constraint comments (the
single-boundary Solution build, the `.mesh.*` handle exposure, the
`TimedFullField` MODULE-LEVEL-import shadow trap) carrying terse
campaign-tag prefixes. Per the pilot discipline (cut standalone
TWIN/HISTORY blocks + SPLIT-tails; do NOT churn every embedded
campaign parenthetical in constraint prose — negligible line gain,
non-zero constraint risk), left in place.

---

## Before / after metrics (measured, `tokenize` + `ast`)

| metric | solver.py before | after | Δ | iteration.py before | after | Δ |
|---|---|---|---|---|---|---|
| total lines | 3137 | 2965 | −172 | 1334 | 1318 | −16 |
| **docstring lines** | 1164 | 1103 | **−61** | 715 | 713 | −2 |
| **comment lines** | 621 | 520 | **−101** | 183 | 171 | −12 |
| prose share | 56% | 54% | −2pt | 67% | 67% | 0 |
| **code tokens** | **6743** | **6743** | **0** | **2168** | **2168** | **0** |

`git diff --stat`: 93 insertions, 281 deletions, 2 files (net −188 ln).

The cut is comment-heavy (−101 solver.py comments) because the dominant
HISTORY surface on these driver files is standalone retirement tombstones
+ campaign-status blocks (`#`-comments), not docstring derivation. The
docstring cut (−61) is the SPLIT-tail HISTORY net of CONTRACT rephrasing.
Appropriately small vs the pilot's −36% docstring — solver.py is a
CONTRACT file, not a teaching file.

## Gates (all GREEN)

1. `import orpheus.sn.solver, orpheus.numerics.iteration` → **OK**.
2. `pytest --collect-only -q` → **6652** (unchanged).
3. **Code-token invariance vs HEAD** (non-string/non-comment tokens):
   solver.py **6743 == 6743**, iteration.py **2168 == 2168** — IDENTICAL.
   Proves ZERO behavioral/code change.
4. **Sphinx `-E` build** (REQUIRED — solver.py IS automodule'd at
   `api/discrete_ordinates.rst:25`): EXIT 0, **0 WARNING / 0 ERROR /
   0 CRITICAL — UNCHANGED from the pre-edit `-E` baseline (0)**. The
   automodule'd solver.py docstrings render clean; no malformed RST.
   (iteration.py is not automodule'd → docstring cuts invisible to the
   build anyway.)
5. **Pointer resolution**: no NEW `:ref:`/`:eq:`/`:func:`/`:class:`/`:meth:`
   refs invented — all kept refs are pre-existing; new pointers are
   plain-text prose ("the SN verification theory page"). The clean `-E`
   build confirms every warning-class ref (`:ref:`/`:doc:`/`:eq:`/`[Key]_`)
   still resolves; `:func:`/`:class:`/`:meth:` render plain-text per the
   `:noindex:`-automodule convention (L-002).

## The hardest judgment calls (calibrate future driver-file batches)

1. **A driver file's cut surface is COMMENTS (tombstones + status
   blocks), not docstrings.** The pilot's teaching-file heuristic
   (docstring −36%) does NOT transfer: solver.py docstrings are the
   estimator/convergence/threat-model CONTRACT a modifier needs (kept
   near-whole), while the HISTORY lives in standalone `#`-comment
   retirement tombstones + campaign-status blocks. Result: comment
   −101 (−16%) dwarfs docstring −61 (−5%). For a driver file, hunt the
   `#`-comment tombstones first, not the docstrings.

2. **The classifier's SPLIT/MOVED verdicts on the `compute_*`
   estimators are WRONG — they are CONTRACT whole.** compute_keff /
   compute_production_rate / _boundary_leakage_rate carry the R7 role
   split, #291 leakage, balance identity, scale-bridge, #251 face-moment
   guard — every one a wrong-simplification guard (L-033 test c) a
   modifier of the eigenvalue path needs at point of use. The book ALSO
   carries the teaching (solver.rst §sn-keff-update / §sn-leakage-functional),
   but that makes the docstring CONTRACT-that-is-also-TWIN → the CONTRACT
   stays inline (the pilot's rule: even when the book carries the
   teaching, the method's own contract stays). Zero cut here.

3. **SPECIAL 1 needed VERIFICATION, not a fix.** The row-8 dual-A bridge
   was ALREADY fully stated in iteration.py's module head (posing +
   A=invertible-resolvent + SN binding A=L+C gains (S,B) → L+C−S−B +
   fission-never-a-gain). The brief said "if MISSING or drifted, FIX" —
   it was neither. The disciplined move was to READ notation.rst row 8,
   READ the module head, confirm the match, and REPORT satisfied — not
   to touch a correct CONTRACT block.

4. **A rebalance pass surfaces a Cardinal-Rule-1 staleness bug — FIX
   it.** The source-helper comment carried a stale path
   `orpheus/sn/scattering.py`; the class lives at
   `transport/operators/scattering.py` (solver.py:71). A raw path string
   in a comment never warns (`-W` blind), so only reading-while-trimming
   catches it. Fixed to a class ref (Cardinal Rule 1 supreme), folded
   into the tag-trim, and flagged in the return.

5. **The category-mistake CONSTRAINT survives a rename-narration cut.**
   iteration.py's KrylovAcceleration "Previously named ``inverter`` …
   (R-1 Step B, 2026-05-19)" is rename HISTORY, but the DISTINCTION it
   guards (a GMRES left-preconditioner ≈ the FULL inverse ≠
   SourceIteration's exact inverse STEP) is CONTRACT — a modifier who
   "tidies" `preconditioner`→`inverter` would conflate them. Cut the
   dated rename story, KEPT the distinction (folded tersely into the
   definition ¶). L-033 test (c) in action.

## Reported discrepancies (not defects introduced)

- **STALE PATH FIXED** (Cardinal Rule 1): `orpheus/sn/scattering.py` in
  the source-helper comment → the class ref. The live module is
  `orpheus.transport.operators.scattering`. Pre-existing; fixed in-pass.
- The 3 `PR-TYPED-5: build typed Solution` inline comments retain their
  campaign-tag prefix by scope-discipline choice (constraint comments,
  see "scope ceiling" above) — reported, not silently trimmed.
