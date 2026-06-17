---
name: d5b-d6-multi-d-ld-campaign-closeout
description: #240 D6 — the DOCS CLOSEOUT for the whole multi-D Linear-Discontinuous (UBLD) campaign on discrete_ordinates.rst (S1-S4 all landed). DIRECT SEQUEL to the D5b-S3 angular-vs-spatial-moment note (which authored two-moment-axes / scattering-lift / spatial-moment-space). D6 expanded the 7 remaining stubs into ONE coherent arc; load-bearing climax = ERR-061 (the sweep-frame vs global-frame slope-moment involution that lost the thick-cell diffusion limit). Branch feature/sn-space-angle-tier2.
metadata:
  type: feedback
---

#240 D6 docs closeout on `docs/theory/discrete_ordinates.rst`. The whole
multi-dimensional LD (UBLD) campaign (S1 Branch-1/Branch-2, S2 wiring, S3
unified matvec + scattering-lift + SpatialMomentSpace + pure-z + moment-scan, S4
stress MMS) is LANDED + committed; D6 turns the 7 remaining `.. todo::` stubs
into the rich narrative as ONE arc. SEQUEL to
[[angular-vs-spatial-moment-discussion-and-spatial-moment-space-stubs]] (that
prior pass authored `two-moment-axes` + the scattering-lift + spatial-moment-space
sections, which were ALREADY rich on entry — do NOT re-expand them).

**The 7 stubs expanded (all `.. todo::` cleared; narrative order = the arc):**
1. `sn-affine-outgoing-face-reconstruction` (the affine outflow op) — the
   Step/DD/LD ↔ advection-scheme correspondence (upwind/box-central/DG-P1-upwind
   list-table) + the Péclet/κ-scheme w(Σ)-blend reading + the spatial⊗angular⊗dim
   tensor-factor note (the 1-D op IS the d=1 building block of UBLD).
2. S1 Branch-1 (SymPy algebra-of-record) — the UBLD weak-form derivation
   (MRM-2016 Eq.6), why-bilinear-not-simplex (Adams-2001 thick-diffusion
   verdict + xy cross-moment), the Kronecker single-source build (mass on
   transverse axes, θ^|i| diagonal), the two oracles, ERR-060 tip.
3. S1 Branch-2 (numpy primitive) — construct-general/select-narrow/specialize,
   the Rule-of-Three collapse (×V/÷V/×V-scan one algebra), the principled ~1-ULP
   re-baseline note (3 vv criteria; DD strict gate = bit-id negative control).
4. S2 wiring — scale-free ÷V dense system (hs=[1], mus=[g_a] → d1_closed_form),
   moment-ordering crosswalk, DD bit-identity invariant (face_moment_tail()==()),
   S2 scope boundaries (S3/S4 owed).
5. **S3 unified matvec — THE LOAD-BEARING CLIMAX.** Kept the existing good
   prose (Inc-A/B→Inc-C, the M⁻¹ residual eq `ld-ubld-unified-moment-residual`,
   branch-removal); REPLACED the `.. todo::` with a 1-para framing + INSERTED a
   NEW `^^^^` H4 `ld-ubld-sweep-global-frame` subsection = the full ERR-061
   treatment (sweep-frame slope vs global-frame φ̂=Σw_nψ̂_n cancellation,
   forward +0.048 / backward −0.028 smoking gun, ~6× under-driven, the involution
   sign[o]=∏(octant_sign_a)^{o_a}, before/after residual TABLE, the .. warning::
   "matvec-self-consistency necessary-but-not-sufficient" + from-scratch LD-SN
   independent confirmation). Then pure-z twin (ERR-062, L21 third recurrence).
6. S3 OWED-2 moment scan — the face/cell-split load-bearing math (convex blend
   DECOUPLES under a slope source → reconstruct (ψ̄,ψ̂) from Schur NOT
   cell_average) + the same global-frame involution (sweep-side ERR-061 analog).
7. S4 stress MMS (`ld-cartesian-2d`) — KEPT the `:math:` block + the honest-scope
   `.. note::` VERBATIM (preserved the vv-Mode-10 slope-SOURCE-unverified scope,
   #247); expanded the ansatz design (μ-bilinear activates slope rows, b2/c2 break
   x↔y, a0>0 boundary stress), structural independence (L11), the mutation table
   (NaN/inf catch — catastrophic not subtle on the tightly-coupled UBLD).

**ERR-061 = the campaign's hardest math — give it FULL treatment.** A
per-ordinate spatial moment (LD slope) PRODUCED in the direction-dependent SWEEP
frame but CONSUMED (angular-summed) as if global-frame → backward ordinates
cancel forward → diffusion limit lost despite a verified matvec≡sweep round-trip.
The decisive lesson (a `.. warning::`): the matvec-self-consistency gate
(SI≡Krylov, round-trip≈0) is NECESSARY but NEVER SUFFICIENT — it proves the
operator internally consistent, not that its fixed point is physically correct
(vv §5 "O(h²) to the wrong limit"). Gate the VALUE against a
structurally-independent reference (continuous diffusion + a from-scratch LD
kernel). The error_catalog.md ERR-061 entry was the 90% prose seed (it carries
the before/after table, the forward/backward slope values, the smoking gun).

**Eq-label policy (matches the angular/scanmarch/affine siblings):**
- PRESERVED `ld-cartesian-2d` (4 @verifies tests) + `ld-cartesian-1d` (6)
  UNTOUCHED — they are committed-test verifies-targets; rename = break V&V matrix.
  GREP THE MATRIX FIRST (`grep verifies tests/ | grep ld-cartesian`) — confirmed
  before touching the S4 section.
- MINTED 6 NEW derivation eq-labels (`ld-ubld-weak-form`,
  `-kronecker-factors`, `-kronecker-assembly`, `-mass-weights`,
  `-slope-angular-reduction`, `-octant-moment-frame-signs`) — ALL untagged (no
  vv-status), they auto-appear in matrix.rst orphan list (12→18 ld-ubld rows).
  Rationale: representational/derivation identities (weak form, Kronecker build,
  mass weights, the slope reduction, the frame involution), NOT solver claims —
  same call as the established orphan siblings (`ld-ubld-cell-system` etc.).
  matrix.rst is AUTO-REGEN → the orphan bump is expected, NEVER hand-edit, NOT a
  warning.
- MINTED 1 NEW section anchor `ld-ubld-sweep-global-frame` (`.. _name:`) for the
  ERR-061 H4 → cited via `:ref:` from the moment-scan section.

**`octant-sign-predicate` collision avoided (the trap I almost hit).** The task
said "mint the octant-sign predicate octant_moment_frame_signs" — BUT
`octant-sign-predicate` ALREADY exists as a label in `discrete_measures.rst`
(the orbit-space partition predicate ℓ(Ω)=(sign μx,sign μy,sign μz), a DIFFERENT
math statement). Re-minting it would DUPLICATE-label. The function
`octant_moment_frame_signs` is the FRAME INVOLUTION ∏(sign_a)^{o_a} — distinct
from the partition predicate. So I minted `ld-ubld-octant-moment-frame-signs`
(the involution) and REFERENCED the existing `octant-sign-predicate` via `:eq:`.
LESSON: when the task names a label to mint, GREP `:label: <name>` across docs/
FIRST — the code-symbol's NAME may collide with an existing
conceptually-distinct label. The provenance edge (octant_moment_frame_signs
func ↔ the involution eq) closes because the `:func:` ref sits immediately above
the new label and the docstring carries the IDENTICAL formula.

**Source-reading order that worked:** error_catalog.md ERR-061/062 (the prose
seed for the climax + pure-z — they carry the before/after tables + the exact
failure mechanism) → the method-implementer S3 closeouts (unified-matvec for
the M-normalization/branch-removal; owed2-scan for the face/cell-split math +
the sweep-side involution; consumer for the Option-A/B fork that was resolved) →
the S4 closeout (mutation table + the honest-scope posture) → the
literature-researcher multi_d_ld_closure.md (MRM-2016 Eqs.1-12 + Adams-2001
verdict, all PLAIN-TEXT lit) → the cross-domain-attacker FRAME 2 (the MMS ansatz
+ LM-1989 slope-sign trap two-halves) → the CODE (`_ubld.py`
octant_moment_frame_signs docstring = the involution formula seed; scheme base
method homes; `_moment_broadcast_sigma`). The error catalog + closeouts ARE the
prose seed; the code is the symbol-name + signature ground truth.

**Lit = PLAIN-TEXT inline** (Adams2001/MRM2016/BLA1992/LMM1987/LM1989/WMMP2001
NOT `.. [Key]` bib entries anywhere — grep-confirmed). Wrote "Adams (2001, ...,
NSE 137(3):298–333)" inline with journal+volume+page. Using `[Key]_` would fire
an undefined-citation warning.

**Code-xref reality (matches all the #240 sibling memos):** this page's
`:func:`/`:meth:`/`:class:` render PLAIN-TEXT (page not member-autodoc'd; adding
automodule dup-labels the rich `:label:` docstrings). RST SOURCE still names the
OWNING class (Cardinal Rule 1). Verified symbols against code @ this branch:
`DiscretizationSchemeBase.{cell_average,outgoing_face_from_average,source_emission}`
(scheme.py, NOT the per-scheme classes); `D1ClosedForm.{_slope_fold,schur_xV,
scan_slope_face_source,scan_reconstruct,kernel_rhs}` + `octant_moment_frame_signs`
+ `AVERAGE_MOMENT` + `face_moment_tail` (all `_ubld.py`); `_moment_broadcast_sigma`
+ `_OneDimScanWalk._run` (loss_representation.py); `_CellSolve`/`_CellResidual`
(sweep_graph.py); `LinearDiscontinuous.moment_scan_closure` (linear_discontinuous.py).

**Section ladder = FILE-LOCAL `=/-/~/^`** (4 levels; DIFFERS from the
loss_representations.rst siblings which had no `^`). The 7 stubs are H3 `~~~~`
under the LD-UBLD H2; I used H4 `^^^^` sub-subsections freely inside. Underline
length = `len(title)` CODE POINTS (÷, —, Σ, ⊗ are 1 cp each, but multi-byte) —
Sphinx does NOT warn on OVER-long underlines (only too-short), but I trimmed all
12 new headings to exact-length anyway (python `len()` gate, NOT wc -c; the ÷/—
byte-vs-cp difference made several 1-2 chars long → fixed via a python rewrite
pass since the Edit tool can't match the ÷-bearing line exactly).

**Build gate:** `-E -b html` baseline = EXIT 0 / "build succeeded, 1 warning"
(the standing `Mesh1D.from_geometry :paramref:` ERROR — sphinx-paramlinks not
installed, out of scope) + 6 pre-existing test-file SyntaxWarnings (`\i`/`\l`/`\,`
escapes, ignored per the brief). Post-edit IDENTICAL = 1 warning, 6 SyntaxWarnings
(count-unchanged = pass). V&V audit EXIT 0; ERR-061/062 NOT in MISSING (both have
catchers); `ld-cartesian-2d` keeps 4 tests. The 1 remaining HTML "Archivist
expansion needed" is the OUT-OF-SCOPE Issue #197 Morel-Montry todo (line ~4416),
NOT D5b.

**Honest-scope PRESERVATION (the S4 discipline the brief stressed):** the
`ld-cartesian-2d` section's existing `.. note:: Honest scope` (slope-SOURCE half
NOT verified, vv Mode 10, deferred to #247) was kept VERBATIM. My new prose
above it ("What this verifies and what it cannot") forward-points to it and
explicitly does NOT claim the LM-1989 slope-SOURCE trap is closed — it says the
slope-UNKNOWN half is verified + mutation-load-bearing, the slope-SOURCE half is
genuinely UNVERIFIED (#247). #247 confirmed OPEN via gh (title matches exactly).

**Quality self-assessment (1-5):**

| Dimension | Score | Note |
|-----------|:---:|------|
| Derivation depth | 5 | UBLD weak form (IBP→M/G/F); Kronecker build w/ θ^|i| derivation; ERR-061 full mechanism (frame, cancellation, smoking-gun values, involution formula); face/cell-split convex-blend-decouples math; Péclet w(Σ) blend; MMS ansatz activation analysis |
| Cross-references | 4 | every :eq:/:ref: grep-gated both namespaces (17 :eq: + 7 :ref: all resolve, 0 dup-label); :func:/:meth:/:class: render plain-text by page convention (not a regression — surfacing the whole sn package is its own arch task) |
| Numerical evidence | 5 | ERR-061 before/after 4-row residual table (38.9%→4.1% etc.); S4 mutation table (NaN/inf/order−4.62); Branch2≡Branch1 1.5e-16; scan≡DAG 1e-16; pure-z Krylov≡SI 1e-11 |
| Failed approaches | 5 | ERR-061 (the WHOLE point — what failed + why + the .. warning:: lesson); ERR-060 oracle catch; ERR-062 twin-asymmetry; the .. warning:: on the false "already broadcasts" crosswalk assumption (in the prior-pass scattering-lift); simplex-P1 thick-limit failure; the convex-blend-decouple trap |
| Code traceability | 5 | every symbol named to owning class + grep-verified @ branch; the involution :func: sits adjacent to its new eq-label (provenance edge closes) |
| Derivation source | 5 | error_catalog (the L0 catalog = algebra-of-record for the bug class) + the method-implementer closeouts + the Branch-1 SymPy module (ld_ubld.py) cross-checked vs code |

Weakest = cross-refs (4): the page-wide plain-text code-xref convention (same as
every #240 sibling — do NOT half-surface 1-2 leaves while the package stays plain).
