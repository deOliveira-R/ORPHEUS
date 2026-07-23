# P2-B — `orpheus/sn/loss_representation/` code-prose rebalance

Phase 2 of #231, batch B (second after the P2-A `scattering.py` pilot).
Branch `docs/sn-doc-architecture`. Two files:

1. `orpheus/sn/loss_representation/__init__.py` (4371 ln)
2. `orpheus/sn/loss_representation/sweep_graph.py` (1069 ln)

## Headline finding (calibrates the remaining operator/machinery batches)

**ZERO MOVED — the book is exhaustively complete (confirmed by concept grep).**
Every teaching concept in both files has a verified TWIN home in the theory
corpus (one-walk, carrying-mesh, Wave-O/active-trace, coefficient-model,
Carlson/pole-seed, Fork-B2, PR-INDEX, two-frames, removal-form, Blelloch,
unified-moment-matvec, degenerate/pure-z — all grep-confirmed in
`loss_representation.rst` / `cartesian_multid.rst` / `curvilinear_numerics.rst` /
`curvilinear_one_group.rst` / `wavefront_cochain.rst` / `operator_tensor_network.rst`).

**But unlike the P2-A operator-file pilot (docstring −36 %), this batch cuts
MUCH LESS.** `loss_representation` is not an operator file whose teaching is
duplicated — it is the SWEEP MACHINERY (the actual numerical realization of the
`(L+C)` traversal). Its docstrings/comments are overwhelmingly CONTRACT-grade:
precise shape conventions, the Resolution-A double-subtract invariant, the L21
sweep≡matvec twin, the O.4b active-trace boundary block, the moment-frame
involution (the ERR-061 diffusion-limit root cause), the pole-seed physics,
aliasing-safety guards. These state THIS method's local contract *referencing* a
book concept — they are not the concept's teaching, so they STAY. A `_run_transpose`
comment is the algebra-of-record for the hand-transposed reverse scan; cutting it
would be a Cardinal-Rule-1 hazard.

The genuinely-duplicated prose that IS cut:
- module-head **carve history** (HISTORY → pointer) + **governing-principle essay**
  (TWIN → compressed rule + pointer) + the "historically scattered" motivation;
- **campaign-relocation provenance** in method docstrings ("relocated verbatim off
  `_MSpatialOperatorSum`", the "#206 Phase B/C", "#196 PR-INDEX arc", "renamed from
  `_sweep_2d_wavefront`") — the CONTRACT stays, the provenance sentence goes;
- `sweep_graph.py` module-head **tensorial-framing** teaching (TWIN → the constraint
  "no Python loop over ordinates/cells/groups; ordinate axis internal to the kernel"
  stays, the tensor-notation essay goes);
- a SELECTIVE set of the **repeated moment-axis re-teachings** (~10× "(#240 D5b-S3)…
  the iterate is the SSOT for the width / φ̂ travels between sweeps") → trimmed to the
  terse shape fact + byte-identical negative-control + one pointer.

## automodule DIVERGENCE from the pilot (build gate is LIVE)

`loss_representation` IS `automodule`'d — `docs/api/discrete_ordinates.rst:82`
(`:members: :undoc-members: :show-inheritance: :noindex:`). Unlike the pilot's
un-automodule'd `scattering.py`, the docstrings ARE rendered, so a malformed
docstring after an edit breaks the `-E -W` build. The Sphinx gate is REQUIRED here.
(`:noindex:` means no xref *targets* are registered — cross-refs TO these symbols
render plain-text elsewhere — but the docstrings still render, per L-002.)

Both files carry **ZERO `.. math:: :label:` and ZERO `vv-status`** (grep-confirmed),
so cutting orphans no `:eq:` target and no `verifies(...)`/`catches(...)` marker
(same safety as the pilot on that axis).

## Pointer convention

Literal greppable `docs/theory/<part>/<file>.rst §<label>` (ratified for this phase,
per brief). Every cited label grep-verified to resolve. NB judgment call: because
this file IS automodule'd (unlike the pilot), a `:ref:`<label>`` role would render
as a working hyperlink and match the file's native idiom — flagged for the main
agent; I used the ratified literal form per the brief.

Anchors cited (all resolve):
`loss_representation.rst`: `loss-rep-history`, `loss-rep-selection`,
`loss-rep-removal-form-matvec`, `loss-rep-four`, `loss-rep-fork-b2`,
`loss-rep-orientation-two-frames`, `loss-rep-one-walk-one-instance`,
`loss-rep-scanmarch-coefficient-model`.
`cartesian_multid.rst`: `ld-ubld-unified-moment-matvec`.
`wavefront_cochain.rst` (module framing), `operator_tensor_network.rst`.

---

## Adjudication table — `__init__.py` (pre-edit line ranges vs the 4371-line file)

Only the blocks that were **CUT / TRIMMED** are itemised; everything else is
KEEP-as-CONTRACT (see the "why so little cut" note). Edit-IDs E1–E14.

| pre-edit lines | block | verdict | action / destination |
|---|---|---|---|
| 15–23 | MODULE — "historically scattered" motivation | HISTORY/TWIN | E1a: compress; the 3-spellings detail → §loss-rep-history |
| 54–75 | MODULE — "governing principle" 3-layer essay | TWIN | E1b: compress to the rule + capability-vs-policy contract + `§loss-rep-selection` |
| 97–116 | MODULE — "Carve history (S0–S6.9)" bullet arc | HISTORY | E1c: cut bullets → one-line pointer `§loss-rep-history` |
| 24–34, 36–52, 77–95, 117–124 | MODULE — posing / hierarchy / selection-SSOT / See-also | **CONTRACT (KEEP)** | front-door orientation + the three-consumers SSOT |
| 364–383 | `streaming_action` protocol docstring | CONTRACT+TWIN | E2: keep single-sourcing contract; σ-affine/Carlson derivation → `§loss-rep-removal-form-matvec` |
| 1252–1254 | `_DAGWavefront` "in S1 … S3 … S4 widens" | HISTORY | E12: cut staging arc; keep "naturally d-general" + current scope |
| 1307–1314 | `MovingFrontierWindow` "SELECTABLE PEER" | CONTRACT+TWIN | E13: keep intentional-keep decision + coverage gate; measured 1.2–1.8× → `§loss-rep-fork-b2` |
| 1810–1853 | `ScanMarch` measured-basis block | CONTRACT+TWIN | E14: keep default-status + 1-D-selection + Mode-9 gate; sweep/matvec basis → `§loss-rep-fork-b2` |
| 2317–2326 | `_OneDimScanWalk` "#206 Phase B relocation" | HISTORY | E5: cut relocation provenance; keep frame-ownership + orientation contract |
| 2510–2524 | `sweep` "Bit-identity contract" | CONTRACT+HISTORY | E7: keep Pattern-2 dual-view + pinning-test pointer; cut `affine_coefficients` provenance |
| 2586–2587 | `_apply_walk` "#206 Phase C relocated verbatim off `_MSpatialOperatorSum`" | HISTORY | E3: cut provenance; keep L21 twin + coefficient-model arms |
| 3249–3260 | `_run` "#196 PR-INDEX-1..5" | CONTRACT+HISTORY | E6: keep the `(N,ng,nx,ny=1)` layout + cache/geom shapes + `:ref:`; cut PR-INDEX arc |
| 3308–3315 | `_run` "Spatial-moment width" comment | CONTRACT (COMMENT) | E8: add canonical pointer `§ld-ubld-unified-moment-matvec`; trim cross-ref chatter + stale "OWED-2" |
| 2939–2941 | `loss_action_transpose` "#206 Phase C relocated verbatim off `_MSpatialOperatorSum..._transpose`" | HISTORY | E4: cut provenance; keep reverse-substitution contract |
| 4218–4226 | `_sweep_scheduled` "Phase 5b storage-B" para | TWIN | E9: keep oracle-retention contract + `§loss-rep-four`; cut the O(...) storage teaching |
| 4318–4322 | `_sweep_jacobi` "renamed from `_sweep_2d_wavefront`" | HISTORY | E10: cut rename history; keep BARE-sweep + all-3-representations contract |
| **~2600–4160 bulk** (`_apply_walk`/`_run`/`_run_transpose` bodies) | the hand-transposed adjoint + Wave-O + carrying-mesh + pole-seed comments | **CONTRACT (KEEP)** | algebra-of-record for the reverse scan + O.4b active-trace + moment-frame involution — cutting = Cardinal-Rule-1 hazard |

## Adjudication table — `sweep_graph.py` (pre-edit line ranges vs the 1069-line file)

| pre-edit lines | block | verdict | action / destination |
|---|---|---|---|
| 28–43 | MODULE — "Tensorial framing" tensor-notation essay | CONTRACT+TWIN | E11: keep the "no Python loop over ordinates/cells/groups" vectorisation invariant; tensor-notation → `§loss-rep-four` |
| 1–16 | MODULE — §15A.2 primitive intro + derived-object | **CONTRACT (KEEP)** | module's raison d'être + "mesh-time not per-call" constraint |
| 18–26 | MODULE — Cardinal-Rule-2 MoC framing | **CONTRACT (KEEP)** | "no shared `SweepGraph` Protocol — do NOT add one" |
| 45–71 | MODULE — References + See-also future-direction | **CONTRACT (KEEP)** | literature + the `assert_*` invariant set the tests pin |
| 89–121, 128–174, 181–366, 374–1063 | `_reframe` / `OctantLabel` / frontier plan / DAG / cell ops | **CONTRACT (KEEP)** | moment-frame involution (ERR-061 root cause), cochain algebra, DAG topology, byte-identical negative control — all constraint |

---

## Verdict counts

`__init__.py` (14 CUT/TRIM edits E1–E14):
- **HISTORY (cut → pointer)**: 6 (module carve-history; `_DAGWavefront` S1/S3/S4;
  `_OneDimScanWalk`/`_apply_walk`/`loss_action_transpose` relocation provenance;
  `_sweep_jacobi` rename).
- **TWIN (cut → pointer)**: 4 (governing-principle essay; `streaming_action`
  σ-affine derivation; `_sweep_scheduled` storage-B; ScanMarch/MFW measured basis).
- **CONTRACT-TRIM (kept, provenance/chatter trimmed)**: 4 (`sweep` bit-identity;
  `_run` PR-INDEX layout; `_run` moment-width comment; the module-head motivation).
- **CONTRACT (KEEP, untouched)**: the overwhelming majority — every method contract,
  the `_run`/`_run_transpose`/`_apply_walk` bodies (the reverse-scan algebra-of-record),
  the O.4b active-trace, the moment-frame involution, all guards.
- **MOVED**: **0**.

`sweep_graph.py` (1 edit E11):
- **TWIN (cut → pointer)**: 1 (module tensorial-framing tensor-notation → `§loss-rep-four`).
- **CONTRACT (KEEP)**: everything else — the file is pure sweep-DAG contract.
- **MOVED**: **0**.

## Before/after (measured, `ast` docstring + `tokenize` comment)

| file | metric | HEAD | after | Δ |
|---|---|---|---|---|
| `__init__.py` | total | 4371 | 4328 | −43 |
| | docstring | 1385 | 1342 | −43 (−3.1 %) |
| | comment | 947 | 947 | 0 |
| `sweep_graph.py` | total | 1069 | 1064 | −5 |
| | docstring | 433 | 428 | −5 |
| | comment | 180 | 180 | 0 |

`git diff --stat`: 94 insertions, 142 deletions, 2 files (net −48). **Docstring
−48 (−2.6 %) — deliberately ≈14× smaller than the pilot's −36 %** (see headline).

## Gates (all pass)

1. `import orpheus.sn.loss_representation, …sweep_graph` → **OK**.
2. `pytest --collect-only -q` → **6652** (baseline 6652, unchanged).
3. **Token invariance** (non-string/non-comment stream vs HEAD): `__init__.py`
   **12061 == 12061**, `sweep_graph.py` **2937 == 2937** — IDENTICAL count AND
   sequence. Zero code change proven.
4. **Cited §label resolution**: all 6 distinct labels resolve
   (`loss-rep-{selection,history,removal-form-matvec,four,fork-b2}`,
   `ld-ubld-unified-moment-matvec`); the kept `:ref:theory-sn-index-convention`
   resolves (`conventions/indexing_and_layout.rst`).
5. **NO theory-page edits** — `git diff --name-only -- docs/` empty (MOVED reported, not executed).
6. **Sphinx `-E -W` build** (REQUIRED — the file IS automodule'd): pre-edit baseline
   = build succeeded / 0 W/E/C / exit 0; post-edit = **build succeeded / 0 W/E/C /
   exit 0** — W/E/C set unchanged.

## Hardest judgment calls (calibrate the machinery-file batches)

1. **A sweep-MACHINERY package is NOT the pilot's operator file — expect a ~14×
   smaller cut, and that is CORRECT, not under-delivery.** The pilot cut −36 % because
   an operator file's *teaching* is 100 % TWIN in the operator-algebra book. This file's
   prose STATES the local contract that *references* a book concept ("returns the FULL
   loss (L+C)ψ, Resolution A" is not teaching Resolution A — it's THIS method's return
   contract). The "would a modifier who never leaves the file do the wrong thing without
   this line?" test keeps it. ZERO MOVED still holds (grep-confirmed), but the TWIN-cut
   surface is the module-head essays + campaign provenance, not the method bodies.
2. **The `_run_transpose` / `_apply_walk` comment bodies are the algebra-of-record for
   the hand-transposed reverse scan — KEEP even though they read like narration.** The
   cotangent routing, the seed-fold transpose, the degenerate-diagonal adjoint, the O.4b
   active-trace block: a modifier editing the adjoint MUST have these; cutting them is a
   Cardinal-Rule-1 hazard (the "3 constraint-bearing blocks misgraded HISTORY" the brief
   warned of, at scale). The classify catalog marked most of these COMMENT-cut [L]; every
   one re-adjudicated to CONTRACT.
3. **The moment-frame involution comments (`_CellSolve.cell`, `_CellResidual.cell`,
   `_reframe`, the ~10 "trailing 2^d axis; DD byte-identical" annotations) are the
   ERR-061 diffusion-limit root cause — CONSTRAINT, not narration.** I added the canonical
   convention pointer (`§ld-ubld-unified-moment-matvec`) at ONE primary site (`_run`) and
   trimmed only cross-ref chatter + a stale "OWED-2" token; the shape fact + the
   byte-identical negative control + the single-source-via-`face_moment_tail` STAY at every
   site. Trimming these to save 1–2 lines each would strip a real invariant.
4. **Measured performance numbers duplicated 3–4× → single-source to the canonical theory
   home, but keep ONE inline.** The Fork-B2 basis (0.57–0.84×, etc.) is TWIN with
   `§loss-rep-fork-b2`; I pointed the two *class* docstrings there but LEFT the headline
   number in `default_for` (the selection factory — the point-of-decision). DRY (Cardinal
   Rule 2) without stripping evidence from the code.
5. **automodule + `:noindex:` makes the Sphinx gate LIVE (pilot divergence).** Unlike the
   un-automodule'd `scattering.py`, this package renders under
   `discrete_ordinates.rst:82` (`:members: :undoc-members: … :noindex:`), so a malformed
   docstring breaks the build — I ran the `-E -W` gate both sides. `:noindex:` means no
   xref *targets* (cross-refs TO these symbols render plain-text elsewhere), but the
   docstrings still render. Both files carry 0 `.. math:: :label:` / 0 `vv-status`, so
   cutting orphaned no `:eq:`/`verifies` target.

## MOVED candidates

**None.** Concept grep confirmed every teaching concept has a book home (see headline).

## Discrepancies vs the brief

- **The classify catalog's file-split header is misplaced.** Its `## sweep_graph.py`
  header sits above content whose line ranges (2194–4348) are still `__init__.py` (the
  second half); the real `sweep_graph.py` content is the un-headed `0–72 MODULE …` block
  at the very end. I used the LIVE files as ground truth (block boundaries/content-ids
  were otherwise accurate, as the brief said).
- **Cut magnitude ≪ pilot.** The brief framed P2-B as "sibling operator batches = TWIN-
  CUTTING". `loss_representation` is a machinery PACKAGE, not an operator surface, so the
  honest outcome is −2.6 % docstring (concentrated on module-head essays + campaign
  provenance), not the pilot's −36 %. Flagged as judgment call #1 for the reviewer.
- **Pointer form in an automodule'd file (judgment call).** I used the ratified literal
  `docs/theory/<part>/<file>.rst §<label>` form per the brief. Because this file IS
  automodule'd (unlike the pilot), a `:ref:`<label>`` role would render as a working
  hyperlink and match the file's native `:ref:`/`:doc:`/`:meth:` idiom — the main agent
  may prefer to convert my literal pointers to `:ref:` roles. Surfaced, not executed.
