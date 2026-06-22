---
name: issue-247-legA-mode10-slope-source-closeout
description: "#247 Leg A docs close-out on discrete_ordinates.rst — DIRECT SEQUEL to the D6 honest-scope note that DEFERRED this (the slope-SOURCE external half of the LM-1989 trap). Replaced the `.. todo:: Archivist` inside the existing ld-cartesian-2d-slope-source honest-scope note with a rich `^^^^` Leg A subsection (the FIRST vv Mode-10 close-out where the converged value is genuinely sub-floor → the WHOLE narrative is about WHY the teeth are NOT the converged flux). New H5 `'''` marker bit me (ladder collision). Branch refactor/sn-foundation-cleanup."
metadata:
  type: feedback
---

#247 Leg A docs close-out on `docs/theory/discrete_ordinates.rst`. SEQUEL to
[[d5b-d6-multi-d-ld-campaign-closeout]]: D6 kept the `ld-cartesian-2d` section's
honest-scope `.. note::` VERBATIM (slope-SOURCE half deferred, vv Mode-10, #247).
THIS task IS that deferral closing — the method-implementer landed Leg A (BULK
external moment-source consumption) and dispatched me to expand the
`.. todo:: Archivist expansion needed (Leg A, #247)` that sat at the BOTTOM of the
honest-scope note. The note's three bullets (Leg A verified / scattering M4 /
Leg B #251) were already written by the method-implementer — I REPLACED only the
`.. todo::` with the rich subsection, and rewrote the todo into a `:ref:` forward-
pointer ("the full Leg A narrative is the subsection :ref:`ld-cartesian-2d-legA`
immediately below").

## What made this one distinct (the vv Mode-10 narrative is the WHOLE point)

This is the FIRST close-out where the **converged value is genuinely sub-floor
sensitive to the bug class** — so unlike ERR-061 (D6, where the matvec-self-
consistency was the gate and the VALUE was the missing teeth), here BOTH the
order leg AND a value-band are blind, and the entire pedagogical content is WHY
and HOW the teeth come from elsewhere. The narrative spine = TWO O(1) structural
teeth, NOT the converged flux:
1. the production LIFT threads the projection through at MACHINE PRECISION
   (production-change proof, catches a regression to zeroing — the EXACT #247 bug);
2. a CONSUMED slope-row sign flip moves the converged flux ≫ solver tol
   (consumption proof, catches sign-blindness), PAIRED with the FLAT no-op leg
   that pins the Mode-10 ASYMMETRY (the old scalar gate is correctly blind).
The convergence-ORDER leg is NECESSARY (slope consistency) but explicitly NOT the
sign teeth. Recorded as a current-invariant `.. note::` lesson (per L10: design
rationale, not failed-experiment narrative → belongs in Sphinx). The closeout
memo + gate-spec §0 ARE the prose seed (they carry the probed O(h²)-ride table
8.18e-3→1.18e-4 both ways at orders ~2.0, the 1.4-1.5× constant ratio, the
6e-5 xy smallest-signal). NO new ERR (proactive-gap close, not a caught bug — the
lift correctly zeroed an unverified-but-honest default q̂=0).

## Subsection shape (5 H5 `'''` under one `^^^^` H4 "Leg A")

`.. _ld-cartesian-2d-legA:` + `^^^^` "Leg A — the external slope-moment source Q̂
(#247)" then five `'''` H5:
1. **The tensor-Legendre projection convention (the CRUX)** — the BARE per-volume
   Legendre coeff (q̄=avg slot0; q̂_a=⟨q,P₁⟩/⟨P₁,P₁⟩, NO θ/h/V — kernel mass
   M=diag(h,θh) adds them); d=2 Kronecker `[ψ̄,ψ̂_y,ψ̂_x,ψ̂_xy]` (x=slot2);
   GLOBAL frame (prod reframes per-octant via octant_moment_frame_signs); the
   hand-bilinear closed-form coeff block (2 new eq-labels); L11 leggauss-only;
   the slot-0=cell-AVERAGE-not-CENTRE caveat (differs O(h²) from legacy flat
   producer → cross-check vs fine-quad average NOT the centre producer).
2. **The typed-union bulk widening** — 3-row contract table (flat / moment-resolved
   LD-only / reject); WHY RANK-not-trailing-size; WHY DD/Step rejects (Pattern 4);
   the negative pin (ValueError naming 2^d); the 3-arm lift table; Pattern-6 no-
   callable-projection (would tautology the gate).
3. **The Mode-10 structural-teeth design** — the sub-floor table (correct vs
   flipped, both O(h²)); order-leg-blind + value-band-fragile; the two O(1) teeth;
   the necessary order leg + positive improves-on-flat leg; the lesson `.. note::`.
4. **The mutation-control table (M1–M4)** — M1-M3=NEW external Q̂ consumption,
   M4=EXISTING scattering Σ_s·φ̂ never sign-blind; each "NEW gate RED / FLAT scalar
   gate GREEN" (the asymmetry); -O-safe (Mode 8) finally-reverted; no ERR.
5. **Sources and gates** — prod sites + the 8 gate names + DriftWarning bit-id.

## THE TRAP THAT BIT ME (new lesson): H5 marker-ladder COLLISION

I used `"""` (double-quote) for the new H5 underlines. Build ERROR'd
"Inconsistent title style: skip from level 5 to 7" at TWO PRE-EXISTING sections
(lines 5865/5933) that used `'''` (apostrophe). Mechanism: RST marker levels are
FILE-LOCAL and assigned by FIRST-APPEARANCE ORDER. The file's H5 was `'''`,
appearing only ~line 5865. My `"""` at line 3267 appeared EARLIER → claimed
level 5 → the later `'''` jumped to level 7 (skipping 6) → ERROR at the OLD
sections, NOT mine. FIX = use the marker the file ALREADY uses at that depth
(`'''`), not a fresh one. LESSON (sharpens the AGENT.md "marker ladder is
file-local" rule): before introducing a NEW deeper level, GREP THE FILE for an
existing marker at that depth (`grep -nE '^([punct]){10,}$'` + tally first chars)
and REUSE it. Introducing a never-before-seen marker EARLIER in the file than an
existing same-depth marker SHIFTS the existing one down a level → skip-error at a
section you did not touch. The error line points at the OLD section, not yours —
do not be misled into thinking the OLD section is broken. The D6/scanmarch memos
say "scan the file's ladder before picking a level" for the 4 known levels
(=/-/~/^); this extends it: when you ADD a 5th level, match the file's existing
5th marker. (This file's 5th = `'''`.)

## Eq-label / cross-ref policy (matches all #240 siblings)

- MINTED 2 NEW derivation eq-labels (`ld-cartesian-2d-projection-coeff`,
  `ld-cartesian-2d-bilinear-coeffs`) — UNTAGGED (no vv-status); they auto-appear
  in matrix.rst orphan list (verified lines 664-665). Rationale: representational/
  derivation identities (the projection normalization + the hand-bilinear closed
  form), NOT solver claims — same call as the established ld-ubld orphan siblings.
- PRESERVED `ld-cartesian-2d` (11 verifies tests, grep-confirmed BEFORE touching)
  + the section anchors `ld-cartesian-2d-slope-source` (the honest-scope note) +
  the existing `ld-ubld-*` labels. MINTED 1 section anchor `ld-cartesian-2d-legA`
  (`.. _name:`) → cited via `:ref:` from the honest-scope note's forward-pointer.
- CITED (don't redefine) the D6 labels: `ld-ubld-mass-weights`,
  `ld-ubld-d1-reduction`, `ld-ubld-octant-moment-frame-signs`,
  `spatial-moment-kronecker-order`, `ld-ubld-scattering-moment-lift` (`:eq:`),
  `ld-ubld-sweep-global-frame` (`:ref:`), and `octant-sign-predicate` (`:eq:`,
  cross-doc to discrete_measures.rst — resolved 1 link). The frame involution +
  global-vs-sweep machinery was ALREADY documented by D6 → I REUSE it ("the
  external Q̂ rides the SAME global→sweep machinery") rather than re-deriving.
- matrix.rst AUTO-REGEN (NOT hand-edit); orphan bump is expected, NOT a warning.

## Code-xref + lit reality (matches siblings)

- This page's `:func:`/`:meth:`/`:class:` render PLAIN-TEXT (not member-autodoc'd;
  automodule would dup-label the rich `:label:` docstrings). RST SOURCE still
  names the OWNING module (Cardinal Rule 1). GREP-GATE every cited symbol (-W is
  BLIND to a wrong-but-plain-text code-xref). Verified @ branch:
  `solver._build_fixed_source_rhs` / `_lift_external_source_to_moments`;
  `_ubld.octant_moment_frame_signs`;
  `numerics.moment_layout.is_moment_valued_by_rank` (true home — solver imports it
  via _ubld re-export, but cite the canonical module);
  `sweep_graph._CellSolve.cell` (line 894, .cell method confirmed);
  `scattering.ScatteringOperator._assemble_per_ordinate_source`;
  `linear_discontinuous.LinearDiscontinuous`.
- LM-1989 / Larsen-Morel-Miller NOT a `.. [Key]` bib entry anywhere → wrote
  "Larsen, Morel & Miller (1989)" PLAIN-TEXT inline (grep-confirmed; page bib keys
  are [BaileyMorelChang2010]/[MorelMontry1984]/[AdamsLarsen2002] etc., NOT LM1989).

## Build gate

`-E -b html` baseline = EXIT 0 / "build succeeded, 1 warning" (the standing
`mesh.py Mesh1D.from_geometry :paramref:` ERROR — sphinx-paramlinks not installed,
out of scope). FIRST post-edit build = 1 warning + 3 ERRORS (the 2 H5-collision
skips + the paramref); after the `"""`→`'''` fix, IDENTICAL to baseline = 1
warning, EXIT 0. V&V audit EXIT 0; `ld-cartesian-2d` keeps 11 tests; the 2 new
labels in orphan list; 0 `.. todo::` remain in the Leg A region. All new anchors
render in HTML; all `:eq:`/`:ref:` resolve as links (grep-gated both namespaces).
MAIN checkout (NOT a worktree) → baseline=1 paramref (NOT 11).

## Quality self-assessment (1-5)

| Dimension | Score | Note |
|-----------|:---:|------|
| Derivation depth | 5 | projection normalization derived from the mass-matrix constraint (M=diag(h,θh) ⟹ q̂ is BARE coeff); hand-bilinear closed form (4 coeffs); the sub-floor O(h²)-ride table; the two-O(1)-teeth design; the Mode-10 lesson note |
| Cross-references | 4 | every :eq:/:ref: grep-gated both namespaces + render as links (incl cross-doc octant-sign-predicate); :func:/:meth:/:class: plain-text by page convention (do NOT half-surface) |
| Numerical evidence | 5 | the correct-vs-flipped sub-floor table (8.18e-3→1.18e-4, orders ~2.0 both); _CONSUMPTION_TOL band reasoning (5e7× over FP, 6000× under §0 trap); M1-M4 probed magnitudes (3e-3/1e-2/6e-5/2.6e-3); improves-on-flat 3.4e-3<5.9e-3 |
| Failed approaches | 5 | the WHOLE Mode-10 narrative (order-leg blind, value-band fragile — the two REJECTED gate designs); why the lift-zeroing was honest-not-wrong (no ERR); the slot0-centre-vs-average cross-check trap |
| Code traceability | 5 | every symbol grep-verified to owning module @ branch; the production-change two named sites; the 8 gate names + the foundation sub-gates |
| Derivation source | 5 | the test-architect gate-spec §0/§1 (the probed sub-floor evidence + the projection CRUX) + the method-implementer closeout (the teeth design + mutation table) = the algebra-of-record; cross-checked vs the projector code + solver code |

Weakest = cross-refs (4): the page-wide plain-text code-xref convention (same as
every #240 sibling — surfacing the whole sn package is its own arch-docs task).
