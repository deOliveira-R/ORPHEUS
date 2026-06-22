---
name: issue-251-legB-boundary-trace-closeout
description: "#251 Leg B docs close-out on discrete_ordinates.rst — the BOUNDARY twin of #247 Leg A (the bulk slope source). Expanded the `.. todo:: Archivist — expand the Leg B narrative (#251)` inside the ld-cartesian-2d-slope-source honest-scope note into a rich `^^^^` Leg B subsection (8 `'''` H5) mirroring Leg A. The SHARPEST vv Mode-10 yet: improves-on-flat is UNACHIEVABLE (boundary slope sub-floor for ANY value claim), so the O(1)-isolating-companion HALF of the Mode-10 recipe is genuinely unavailable → structural teeth ONLY. Branch refactor/sn-foundation-cleanup. NEW: marker reuse worked first try (learned from Leg A's collision); production docstrings were the dominant prose seed."
metadata:
  type: feedback
---

#251 Leg B docs close-out on `docs/theory/discrete_ordinates.rst`. DIRECT SEQUEL
to [[issue-247-legA-mode10-slope-source-closeout]]: Leg A closed the EXTERNAL
half of the slope-SOURCE trap (bulk `_lift_external_source_to_moments`); Leg B
closes the BOUNDARY half (the transverse face-slope carried on `mesh.trace`).
The method-implementer landed #251 and dispatched me to expand the
`.. todo:: Archivist — expand the Leg B narrative (#251)` that sat at the BOTTOM
of the `ld-cartesian-2d-slope-source` honest-scope note (after the
`ld-cartesian-2d-legA` subsection). I REPLACED only the `.. todo::` block, and
rewrote the note's forward-pointer prose into a `:ref:ld-cartesian-2d-legB`
breadcrumb that lists the subsection's contents (the same idiom Leg A used for
its own `:ref:`).

## Subsection shape (8 `'''` H5 under one `^^^^` H4 "Leg B")

`.. _ld-cartesian-2d-legB:` + `^^^^` "Leg B — the boundary transverse-face-slope
(#251)" then EIGHT `'''` H5 (Leg A had five; Leg B's production change spans more
sites so the shape is wider):
1. **The `boundary_face_layout` moment-tail storage lever (the CRUX)** — the
   ONE-site storage change (no carrier existed, unlike Leg A's free ride on the
   bulk iterate); the 3 properties (TraceSpace ZERO-changes "moment-ready by
   accident" / DD byte-identical via `face_moment_tail(1)==()` / 2-D-and-higher
   by construction since `d-1` exponent vanishes for a 1-D point face); minted
   eq-label `ld-cartesian-2d-face-slot-shape`.
2. **The transverse face-moment normalization** — the apples-to-apples trap
   transposed to the transverse axis (cochain's `assemble_inflow_axis` transverse
   mass = Kron of `mass_1d=diag(h_t,θh_t)` adds the weighting → trace carries BARE
   coeff); minted `-face-projection-coeff` + `-face-bilinear-coeffs`; slot-0
   centre-vs-average subtlety (SHARPER than Leg A — today's scalar trace is the
   CENTRE, projection slot-0 is the AVERAGE, differ O(h²) → gate compares slot-1
   ONLY).
3. **The rank-discriminated inflow lift + the 4 DROP sites** — `_inflow_to_moments`
   3-arm table (DD identity / scalar widen-slopes-zero / moment pass-through) via
   the Leg-A primitive `is_moment_valued_by_flat_rank`; FFW oracle
   `_octant_face_cochain` mirror; the 4 outflow capture-collapse DROP sites table.
4. **The producer — both ranks through one slot** — `prescribed_inflow` 2-arm
   table + the LOAD-BEARING `.. note::` "audit EXISTING scalar producers when you
   widen a slot" (the 8-gates-reddened under-scope lesson).
5. **The sharper Mode-10: structural teeth only, no improves-on-flat leg** — the
   sub-floor A-err table (real WORSE than flat, flipped BETTER — improves? NO);
   the 2 O(1) structural teeth (threading + consumption) with probed tables; the
   current-invariant `.. note::` (the companion-UNAVAILABLE branch of Mode-10).
6. **The negative pin, bit-identity guard, reflective follow-up** — width-reject
   pin; 3-way bit-id guard; the `.. note:: Reflective-BC transverse-slope sign
   (#252)` (storage-correct PermutationOperator passthrough / SIGN unverified /
   tracked in #252).
7. (folded into 6 — reflective note)
8. **Sources and gates** — 3 prod files + 6 gate names, -O-safe.

## WHAT MADE THIS DISTINCT (the Mode-10 sharpening is the WHOLE point AGAIN)

This is the THIRD Mode-10 close-out (after #240 D5b-S4 and #247 Leg A) and the
FIRST where the converged value is sub-floor for ANY value claim, not just the
sign. Leg A's "improves-on-flat" positive leg (probed 3.4e-3<5.9e-3) is
UNACHIEVABLE for Leg B (probed: real slope WORSE, flipped BETTER — a boundary
slope is too localized + O(h)-small to register against the bulk O(h²) floor).
So the canonical Mode-10 recipe's "add a companion gate that isolates the term so
its error is O(1)" half FAILS — there is no regime where a boundary-trace slope
is the dominant forcing. The current-invariant `.. note::` records this as the
companion-UNAVAILABLE branch: structural teeth (producer-threading at machine
precision + consumed-flip ≫ tol, paired with the scalar no-op asymmetry) are the
COMPLETE resolution. This generalizes the Mode-10 row — feed it back to
vv-principles (the test-architect's §13 + method-implementer's §1 both
recommended the one-line addition).

## SOURCE-READING ORDER (production docstrings were the DOMINANT seed this time)

closeout (method-implementer) → gate-spec §0 (the probed sub-floor evidence) →
explorer map (the file:line lever choices) → THEN the PRODUCTION CODE itself.
⭐ NEW vs Leg A: the method-implementer wrote RICH production docstrings on
`boundary_face_layout` (geometry.py:1004 — the whole §"Spatial-moment tail" block),
`_inflow_to_moments` (the 3-arm prose), `_octant_face_cochain`, and
`prescribed_inflow` (the 2-arm comment). These docstrings ARE the algebra-of-record
prose seed — narrate them, don't re-derive. The transverse-mass normalization came
straight from reading `assemble_inflow_axis` (_ubld.py:270 — `transverse_mass` =
Kron of `mass_1d`, applied BEFORE the active-axis trace) + the test helper
`_face_transverse_legendre` docstring (the bare-coeff projector). The probed number
TABLES (A-err sub-floor, consumption |Δφ|/|φ|, linearity in k) came VERBATIM from
the gate-spec §0b/§0c.

## THE MARKER-LADDER TRAP DID NOT BITE (Leg A's lesson applied)

Leg A's memory warned: the file's H5 is `'''` (apostrophe), NOT `"""`. I used
`'''` from the start → ZERO marker-collision error. CONFIRMED the lesson works:
before introducing the H5 underlines, I trusted the Leg-A memory's "this file's
5th marker = `'''`" and reused it. (The H4 is `^^^^`, matching the file ladder.)
The ONLY underline cleanup needed: 4 underlines were 1 cp too LONG (over-runs).
RST tolerates over-runs (Leg A's own underline is 64cp for a 62cp title — 2
over), but I normalized to exact `len(title)` code-points via a python rewrite
(em-dash = 1 cp / 3 bytes — measured with `len()` not `wc -c`, the title "Leg B —
…" is 49 cp). NOT a build error either way; normalized for cleanliness scoped to
MY lines only.

## EQ-LABEL / CROSS-REF POLICY (matches all #240/#247 siblings)

- MINTED 3 NEW derivation eq-labels — `ld-cartesian-2d-face-slot-shape`,
  `ld-cartesian-2d-face-projection-coeff`, `ld-cartesian-2d-face-bilinear-coeffs`
  — UNTAGGED (no vv-status); auto-appear in matrix.rst orphan list (verified
  lines 665-667). Rationale: representational/derivation identities (the slot-shape
  storage formula + the bare-coeff projection normalization + the hand-bilinear
  closed form), NOT solver claims — exactly the Leg-A orphan-sibling call.
- PRESERVED `ld-cartesian-2d` (11 verifies tests, grep-confirmed in matrix line
  463 BEFORE + AFTER) + minted 1 section anchor `ld-cartesian-2d-legB` (`.. _name:`)
  → cited via `:ref:` from the honest-scope note forward-pointer AND from inside
  the Leg B subsection's opening (a back-ref to Leg A).
- CITED (don't redefine): `spatial-moment-kronecker-order` (`:eq:`, the cell-tail
  exponent contrast d vs d-1), `ld-cartesian-2d-projection-coeff` +
  `ld-cartesian-2d-bilinear-coeffs` (Leg A's labels — Leg B's face coeffs are the
  1-D-transverse FACTOR of Leg A's tensor coeffs, framed as such). All resolve as
  HTML links (grep-gated `href="#equation-..."`).
- matrix.rst AUTO-REGEN (NOT hand-edit); orphan bump expected, NOT a warning.

## CODE-XREF + LIT REALITY (matches siblings)

- This page's `:func:`/`:meth:`/`:class:`/`:attr:` render PLAIN-TEXT (not
  member-autodoc'd; automodule would dup-label the rich `:label:` docstrings).
  RST SOURCE still names the OWNING module/class (Cardinal Rule 1). GREP-GATE every
  cited symbol (-W BLIND to a wrong-but-plain-text code-xref). Verified @ branch:
  `geometry.SNMesh.boundary_face_layout`; `moment_layout.face_moment_tail` +
  `is_moment_valued_by_flat_rank`; `_LossRepresentation._inflow_to_moments` (base,
  loss_representation.py:297/360) + `._n_face_moments` + `._spatial_moment_tail`;
  `FullFieldWavefront._octant_face_cochain` (staticmethod, 1213) — the 4 DROP
  sites live on MovingFrontierWindow(1007)+FullFieldWavefront(1180);
  `boundary_source_sink.BoundarySourceSink.prescribed_inflow`;
  `_ubld.assemble_inflow_axis` + `mass_1d`; `solver._lift_external_source_to_moments`
  (the Leg-A twin); `boundary_operator.SNBoundaryOperator._reflect_trace`;
  `numerics.operator.PermutationOperator` (961).
- NO new `.. [Key]` bib added — Leg B is a production-CONSUMPTION widening, not a
  new reference. LM-1989 already plain-text inline from Leg A; reused the
  established page bib keys implicitly (no new lit citation needed).
- #252 (the reflective-sign follow-up) is NEW this session — wrote it as plain
  inline "#252" prose (GitHub issue refs are plain text on this page, like #247/#251).

## BUILD GATE

`-E -b html` MAIN checkout (NOT a worktree) → baseline = EXIT 0 / "build
succeeded, 1 warning" (the standing `mesh.py Mesh1D.from_geometry :paramref:`
ERROR — sphinx-paramlinks not installed, out of scope). Post-edit build =
IDENTICAL: EXIT 0, 1 warning, SAME paramref ERROR, NO new WARNING/ERROR/CRITICAL.
V&V matrix regenerated: `ld-cartesian-2d` keeps 11 tests (line 463); the 3 new
labels in orphan list (665-667). 0 `.. todo::` remain in the Leg B region. All
new anchors render in HTML (`id="ld-cartesian-2d-legb"` + the 3 `equation-*`);
both forward-pointers (`:ref:ld-cartesian-2d-legA`/`-legB`) resolve as links.

## QUALITY SELF-ASSESSMENT (1-5)

| Dimension | Score | Note |
|-----------|:---:|------|
| Derivation depth | 5 | the storage-lever formula (eq-label) with the d-1 vs d exponent contrast derived (face is codim-1); the transverse-mass normalization derived from `assemble_inflow_axis` (BARE coeff because the Kron-mass adds h_t/θ); the hand-bilinear face coeffs as the 1-D-transverse FACTOR of Leg A's tensor coeffs; the sub-floor A-err table + the 2-O(1)-teeth design + the consumption |Δφ|/|φ| + linearity tables |
| Cross-references | 4 | every :eq:/:ref: grep-gated both namespaces + render as HTML links; :func:/:meth:/:class:/:attr: plain-text by page convention (do NOT half-surface the sn package — its own arch-docs task) |
| Numerical evidence | 5 | sub-floor A-err table (1.015e-2 real-WORSE / 1.009e-2 flipped-BETTER); consumption table (4.10e-3 near-bdy / 3.27e-3 global at nc16, halves at nc32); linearity {2.4,4.8,9.5}e-3 in k={.05,.1,.2}; the slot shapes (24,2,6,2) probed; ~5.6 orders over _CONSUMPTION_TOL |
| Failed approaches | 5 | the WHOLE sharper-Mode-10 narrative (improves-on-flat UNACHIEVABLE — the companion-unavailable branch); the producer under-scope (8 gates reddened by the unconditional reject → the audit-existing-scalar-producers lesson note); the slot-0 centre-vs-average trap; reflective-sign UNVERIFIED (#252, not claimed verified) |
| Code traceability | 5 | every symbol grep-verified to owning module/class @ branch; the 3 prod files; the 4 DROP sites; the 6 gate names; the production docstrings narrated (algebra-of-record) |
| Derivation source | 5 | method-implementer closeout (teeth design + producer-carve lesson) + gate-spec §0 (probed sub-floor/consumption/linearity) + explorer map (lever choices) + the PRODUCTION DOCSTRINGS themselves (the dominant seed this time) + read `assemble_inflow_axis`/`_face_transverse_legendre` for the normalization |

Weakest = cross-refs (4): the page-wide plain-text code-xref convention (same as
every #240/#247 sibling — surfacing the whole sn package is its own arch task).

## LESSON (sharpens the AGENT.md Mode-10 + producer-rank rules)

1. The companion-UNAVAILABLE branch of vv Mode-10: when a consumed term has NO
   regime where it is the dominant forcing (a boundary-trace slope is intrinsically
   a localized O(h) perturbation below the bulk O(h²) floor), the
   "O(1)-isolating companion gate" half of the recipe is genuinely unavailable
   → structural teeth ONLY (producer-threading at machine precision + consumed-flip
   ≫ tol + scalar no-op asymmetry) are the COMPLETE close. This is a one-line
   addition to the vv Mode-10 row (NO new mode, NO ERR) — propose it to
   vv-principles. Both the test-architect §13 and method-implementer §1 converged
   on this recommendation.
2. When a docs task narrates a production-CONSUMPTION carve (a slot/trace widening),
   the method-implementer's PRODUCTION DOCSTRINGS are often the richest prose seed
   — they were written WITH the algebra-of-record discipline (the 3-arm
   `_inflow_to_moments`, the §"Spatial-moment tail" block on `boundary_face_layout`).
   Read them BEFORE re-deriving; narrate + cross-link, don't compete (the
   algebra-of-record archivist rule applies to production docstrings, not only
   SymPy modules).
3. Marker-ladder reuse (Leg A's hard-won lesson) WORKS prophylactically: trust the
   prior memory's "this file's 5th marker = `'''`" and reuse it from the start →
   zero collision. The over-run normalization (len(title) code-points, em-dash=1cp)
   is cosmetic, not a build error — scope it to YOUR lines.
