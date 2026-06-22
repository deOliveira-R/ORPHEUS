---
name: issue-257-s9-coherent-promise-criterion-seam
description: "#257 S9 docs — the LD 2-D Cartesian BOUNDARY coherent-promise + the #263 property-vs-type criterion seam. FOURTH consecutive vv Mode-10 companion-unavailable close-out (D5b-S4 → #247 Leg A → #251 Leg B → S9), and the FIRST to add a cross-page FIELD-TYPE-ALGEBRA seam (a NEW = H1 criterion section on operator_algebra.rst). Two deliverables, two pages: a 3rd ^^^^ H4 peer 'S9' subsection on discrete_ordinates.rst (after Leg A/Leg B) + the field-type-vs-property-criterion = H1 section on operator_algebra.rst, bidirectionally :ref:-linked. NO new eq-label (REUSED ld-cartesian-2d); 2 NEW section anchors. Branch feature/field-typed-operator-algebra."
metadata:
  type: feedback
---

#257 S9 docs on TWO pages. DIRECT SEQUEL to [[issue-251-legB-boundary-trace-closeout]]
(and its parent [[issue-247-legA-mode10-slope-source-closeout]]): #251 Leg B
landed the boundary CONSUMPTION path (moment-resolved trace + `_inflow_to_moments`
pass-through end to end); the MMS *producer* still dropped the slope (built a
scalar). S9 closes that producer-blindness (the case `prescribed_inflow` now
EMITS the moment slot) AND adds the #263 field-type-vs-property criterion seam.
The brief was tightly scoped (S9 + #263 seam ONLY — NOT the broader owed
S3b/S5/S6/S8b/S8c consolidated pass), and explicitly told me to verify my
sections render `-W`-clean (build to /tmp; main agent runs the full build).

## THE THREE NARRATIVE POINTS (all three landed, the user's framing)

1. **The COHERENT PROMISE is TRUE and ALREADY DELIVERED — by the AVERAGE moment
   alone.** The load-bearing subtlety: NOT delivered by the transverse
   boundary-slope moment. The prescribed inflow is exact at the face cells; the
   bulk LD closure carries it inward at O(h²); the boundary cell integrates the
   inflow against its own linear basis and the cell-AVERAGE moment is what that
   integral needs to O(h²). So there was NEVER an O(h) deficiency to "recover" —
   reframe S9 from "recover 2nd order at the boundary" → "make the producer
   honest + LOCK the no-asterisk promise". EVIDENCE = `test_first_cell_row_already_second_order`
   (first-cell-row i=0 orders ~2.0 for BOTH flat and mom).
2. **The transverse slope is an inflow-REPRESENTATION refinement (O(h)→O(h²) on
   the face trace) that is genuinely consumed (flip moves flux ~4.1e-3 ≫ tol) +
   threaded at machine precision — but SUB-FLOOR for the converged flux.** So NO
   value/order gate keyed on the slope (would falsely RED a correct term);
   STRUCTURAL teeth only. This is the FOURTH vv Mode-10 companion-unavailable
   instance (D5b-S4 → #247 Leg A → #251 Leg B → S9). NO new vv mode, NO ERR
   (proactive-gap producer-blindness close, not a caught bug). The optical-depth
   + amplitude sweeps (`test_optical_sweep_slope_never_beats_floor`,
   `test_amplified_boundary_slope_still_subfloor`) prove the sub-floor wall is
   FUNDAMENTAL to a boundary-trace moment (codim-1, measure-zero in the limit),
   not regime-specific (amplify 20× → MONOTONICALLY WORSE, never better).
3. **The #263 CRITERION SEAM — a moment earns a TYPE iff a non-canonical DUAL
   coexists.** Angular order PASSES (AngularFlux↔HarmonicMomentField, the
   integrating quadrature makes them non-canonically iso, M/R the applied bridge
   → TWO types). The transverse SPATIAL moment FAILS (one tensor-Legendre tower,
   only morphism = identity) → it is a PROPERTY (the flat face-buffer moment tail
   via boundary_face_layout), exactly as the bulk carries spatial_moments as a
   property. A BoundaryMomentField would be a vacuous naming leaf. First-class
   SpatialMomentField DEFERRED to the collocation trigger (nodal-DG / Lagrange-FEM
   — where a nodal dual + Vandermonde morphism arrive), tracked #263 (OPEN).

## WHAT MADE THIS DISTINCT vs the THREE prior Mode-10 close-outs

- **First Mode-10 close-out to span TWO pages + add a cross-page field-type-algebra
  seam.** #240/#247/#251 were all single-page (discrete_ordinates.rst). S9's #263
  criterion belongs on the FIELD-TYPE-ALGEBRA home (operator_algebra.rst, which
  already discusses HarmonicMomentField as a distinct TYPE + the tensor-product
  algebra), NOT buried in the LD page. I put the criterion as a NEW `=` H1 section
  `field-type-vs-property-criterion` right after the tensor-product family (before
  the bc-tensor-primitives `=` section), and a BRIEF specialised statement +
  forward-:ref: from the LD page's S9 subsection. Bidirectional cross-doc :ref:
  (LD → operator_algebra criterion; operator_algebra → LD coherent-promise) —
  BOTH resolve as HTML hrefs (grep-gated `href="...html#anchor"`).
- **The brief offered a CHOICE for the #263 seam** ("focused subsection in the
  field-type page IF one exists, OR a cross-referenced note from the LD page").
  `operator_algebra.rst` EXISTS and already owns HarmonicMomentField/tensor
  algebra → the criterion's durable home is THERE (future readers reasoning about
  field types find it where the type algebra lives). I did BOTH halves: full
  criterion section on operator_algebra + brief+forward-ref on the LD page. The
  cross-domain-attacker's `spatial_order_type_vs_property_criterion.md` memo IS
  the algebra-of-record prose seed (the 3-clause criterion verbatim).
- **The PRODUCTION DOCSTRINGS + the gate-module docstring were the DOMINANT prose
  seed** (the Leg B lesson held again): the case module HONEST-SCOPE block
  (sn.py:1132-1191, already flipped to the now-CARRIED state by the
  method-implementer), the `prescribed_inflow` + `_project_inflow_to_face_moments`
  docstrings (the bare-Legendre normalization + L11 + collapse trigger), and the
  `test_ld_2d_boundary_promise.py` MODULE docstring (the durable verdict + the
  mechanism record). Narrate them, don't re-derive.

## SECTION SHAPE (3rd ^^^^ H4 peer of Leg A/Leg B; 6 ''' H5)

`.. _ld-cartesian-2d-coherent-promise:` + `^^^^` "S9 — the coherent boundary
promise and the property-vs-type seam (#257)" then SIX `'''` H5:
1. The coherent promise: LD is 2nd-order at the boundary, no asterisk (WHAT
   delivers it = the average moment; reframe of the motivation).
2. First-cell-row evidence (the decisive scaling argument; list-table of orders ~2.0
   flat AND mom; the optical-depth + amplitude sweep summary).
3. Producer honesty: the MMS case emits the moment slot (the face_moment_count gate;
   DD/Step byte-identical; the leggauss L11 projector; the bare-coeff math block
   reusing Eq ld-cartesian-2d-face-projection-coeff; GATE C single-source; the
   second-order test-side-baseline interaction `.. note::`).
4. The fourth Mode-10 instance — why NO value gate (the companion-unavailable
   branch; the Mode-11 sentinel; a `.. warning::` "do NOT write 'S9 recovers 2nd
   order at the boundary'" — the L8/L10 brain-discipline prophylactic).
5. The property-vs-type seam (#263): boundary moment is a PROPERTY (brief criterion
   + forward :ref: to operator_algebra's field-type-vs-property-criterion section).
6. Sources and gates (1 prod file; 4 promise/verdict gates + 2 producer-stamp gates;
   -O-safe; NO ERR).

operator_algebra.rst: `.. _field-type-vs-property-criterion:` + `=` H1 + 3 `-` H2
(the criterion blockquote / angular-PASSES-spatial-FAILS / the defer-with-trigger
#263 decision).

## EQ-LABEL / CROSS-REF POLICY (matches all #240/#247/#251 siblings)

- **NO new eq-label** — S9's verdict is STRUCTURAL narrative; the relevant
  projection labels already exist from Leg A/B. REUSED `ld-cartesian-2d` (the
  promise gate `@verifies` it; grep-confirmed 11→still its count, S9 mints none).
- **2 NEW section anchors** (`.. _name:` → `:ref:`): `ld-cartesian-2d-coherent-promise`
  (the S9 subsection) + `field-type-vs-property-criterion` (the operator_algebra
  criterion section). BOTH grep-confirmed ABSENT before minting; both render as
  HTML ids. NOT eq-labels (no vv-status — section anchors, different namespace).
- **CITED (don't redefine):** `ld-cartesian-2d` (`:eq:`), `ld-cartesian-2d-face-projection-coeff`
  (`:eq:`, Leg B's bare-coeff normalization REUSED by the producer), `ld-cartesian-2d-legA`/`-legB`
  (`:ref:`), `harmonic-moment-projection` (`:eq:`, intra operator_algebra — line 5912,
  AFTER my insert at 1437; backward+forward intra-doc :eq: both resolve in Sphinx).
- matrix.rst AUTO-REGEN — but S9 mints NO eq-label so NO orphan bump (the prior
  3 close-outs each minted derivation eq-labels; S9 is the first that adds none).
  Brief said do NOT regenerate the matrix (main agent's batched build does it).

## CODE-XREF + LIT REALITY (matches siblings)

- discrete_ordinates.rst `:func:`/`:meth:`/`:class:`/`:attr:` render PLAIN-TEXT
  (not member-autodoc'd; automodule would dup-label the rich `:label:` docstrings).
  GREP-GATED every cited symbol @ branch: `SN2DCartesianLDStressMMSCase` (sn.py:1344)
  + `.prescribed_inflow` (1557) + `._project_inflow_to_face_moments` (1644);
  `moment_layout.face_moment_count` (64) + `AVERAGE_MOMENT` (61); `SNMesh.boundary_face_layout`;
  the 5 gate functions (test_ld_2d_boundary_promise.py:170/231/280/316 +
  test_mms_ld_2d.py:1152). operator_algebra.rst cites `AngularFlux` (angular_flux.py:63),
  `HarmonicMomentField` (harmonic_moment_field.py:113), `SpatialMomentSpace`
  (spatial_moment_space.py:136), `MomentField` ABC (_bases.py:461) — all grep-verified.
- `:doc:/theory/reference_solvers` for the algebra-of-record pillar — doc EXISTS,
  resolves (2 hrefs). LITERATURE: none new (Hesthaven-Warburton nodal-DG + NEM/SANM/ANM
  named PLAIN-TEXT inline, NOT `.. [Key]` bib entries — same convention as LM-1989
  on this page). #257/#263/#251/#247/#252 all plain inline "#NNN" prose (GitHub
  refs are plain text on these pages).

## BUILD GATE

MAIN checkout (NOT a worktree) → baseline = EXIT 0 / "build succeeded, 1 warning"
(the standing `mesh.py Mesh1D.from_geometry :paramref:` ERROR — sphinx-paramlinks
not installed, out of scope). FIRST `-E` post-edit build = 1 NEW WARNING ("Title
underline too short" at the property-vs-type H5 — my underline was 64 cp for a
65 cp title). FIX = normalized 3 of my underlines to exact len(title) code-points
(the H5 under-run = the actual warning; the H4 + Sources-and-gates over-runs were
tolerated but I normalized for cleanliness, scoped to MY lines — left the
PRE-EXISTING `Cartesian 2D` 45-vs-43 over-run untouched). Second `-E` build =
IDENTICAL to baseline (EXIT 0, 1 warning = the paramref, NO new WARNING/ERROR/
CRITICAL). All anchors render as HTML ids; all cross-doc + intra-doc :ref: + :eq:
resolve as hrefs (grep-gated); 0 `.. todo::`; admonitions/list-table/math render.

## QUALITY SELF-ASSESSMENT (1-5)

| Dimension | Score | Note |
|-----------|:---:|------|
| Derivation depth | 5 | the WHY-the-average-delivers-it physics (cell integrates inflow against its own linear basis → cell-AVG is O(h²)-adequate); the bare-coeff producer projection reusing Leg B's normalization; the first-cell-row scaling argument as the mechanism oracle; the 3-clause criterion (non-canonical dual + applied morphism) specialised to the boundary moment |
| Cross-references | 5 | first close-out to add a BIDIRECTIONAL cross-DOC :ref: (LD↔operator_algebra) — both resolve; all :eq: (ld-cartesian-2d, face-projection-coeff, harmonic-moment-projection) + :ref: (legA/legB/criterion/coherent-promise) grep-gated as hrefs; :func:/:meth:/:class: plain-text by page convention |
| Numerical evidence | 5 | first-cell-row orders ~2.0 (flat 1.993/2.004/2.001, mom 1.998/2.005/2.003); the optical-depth sweep (Σ_tL 0.1→2.0 × c{0,0.5}); the 20× amplification +17%-worse; the consumption 4.1e-3≫tol; the sub-floor <20% band |
| Failed approaches | 5 | the WHOLE Mode-10 companion-unavailable narrative + the `.. warning::` prophylactic against "S9 recovers 2nd order at the boundary"; the producer-baseline-inheritance trap `.. note::`; the BoundaryMomentField-as-vacuous-naming-leaf rejection; the defer-with-trigger #263 |
| Code traceability | 5 | every symbol grep-verified to owning module/class @ branch; the 1 prod file; the 6 gate names; the production+gate docstrings narrated (algebra-of-record) |
| Derivation source | 5 | method-implementer closeout + test-architect gate-spec + numerics-investigator optical verdict + cross-domain-attacker criterion memo + the PRODUCTION DOCSTRINGS + the gate-module docstring (the dominant seeds) |

## LESSON (sharpens the Mode-10 + field-type-seam doc rules)

1. The FIELD-TYPE-vs-PROPERTY criterion seam belongs on the field-type-algebra
   HOME page (operator_algebra.rst — where HarmonicMomentField + tensor algebra
   live), NOT buried in the consuming solver page. Put the FULL criterion section
   there (a `=` H1 with the 3-clause blockquote + angular-PASSES/spatial-FAILS +
   defer-with-trigger), and a BRIEF specialised statement + forward-:ref: from the
   solver page. Bidirectional cross-doc :ref: — grep-gate BOTH hrefs in the
   rendered HTML (cross-doc dangling renders plain-text, -W blind). The
   cross-domain-attacker's criterion memo IS the verbatim prose seed (the
   3-clause "non-canonical dual must coexist").
2. A vv Mode-10 close-out CAN mint ZERO eq-labels when its content is purely the
   STRUCTURAL/verdict narrative (S9 is the first in the D5b-S4→#247→#251→S9 chain
   to do so — the projection labels it needs already exist from Leg A/B). REUSE
   the verifies-target label, mint only SECTION anchors. → no matrix.rst orphan
   bump (the prior 3 each minted derivation eq-labels).
3. The L8/L10 brain-discipline prophylactic: when the WHOLE point is "the obvious
   reading is WRONG" (here: 'S9 recovers 2nd order at the boundary' — false, the
   average already delivers it), add an explicit `.. warning::` naming the
   forbidden sentence + the correct framing. A future session quoting the page
   for V&V reasoning will hit the warning before the misreading. Mirrors the
   vv-principles "NEVER write 'MMS verifies the eigenvalue'" anti-pattern style.
4. The producer-honesty SECOND-ORDER interaction (a `.. note::` worth keeping):
   when a refactor makes a production path honest about a quantity a test was
   controlling externally, the test's controlled toggle MUST be built TEST-SIDE
   so it does not inherit the production change it is meant to be orthogonal to.
   (General doc lesson distilled from the method-implementer's `_solve_with_boundary_slope`
   re-baseline fix.)
