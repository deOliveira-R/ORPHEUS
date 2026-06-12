---
name: type-retirement-but-concept-survives docs
description: Recording a typed-field RETIREMENT in Sphinx when the carrier dies but the CONCEPT lives on in two native realizations (S6.4(f) WavefrontFlux/InteriorFaceSpace, #222) — retitle to the concept not the type, prominent succession note, de-role every dead module-path role to a literal, past-tense the type, repoint present-tense claims to the live realizations, KEEP the derivation/grids/tables as history
type: feedback
---

The sequel to [[feedback_sweep_walk_collapse_relayering]] (S6.4(e)): the
walk re-layering that collapsed the four-method product ALSO dissolved a
typed interior-face field's two verbs into the shared frame, so the type
(`WavefrontFlux` + its space `InteriorFaceSpace` + its 25-test foundation
suite) was git-rm'd at S6.4(f). The docs task is a TYPE-RETIREMENT sweep,
distinct from a normal de-role: the underlying **concept** (the interior
1-cochain `C¹_int`, the `C¹ = C¹_int ⊕ C¹_∂` biproduct, the `ι_*`/`ι*`
trace algebra) remains VALID THEORY and survives in two native code
realizations. The user's framing is worth honoring verbatim in the doc:
"a useful concept (at least while it lived)".

**The decisive distinction from a plain "retired-feature" sweep**: do NOT
delete the rich content. The derivation, the storage×role×locus grid, the
two-loci biproduct table, the ι-operator table are HISTORY-VALUABLE — they
explain WHY the cochain frame is the right one, which a future session
still needs even though the carrier is gone. Keep them; flip present→past
tense where they speak of the type; repoint the LIVE present-tense claims
to the realizations.

**The succession shape (the reusable moves)**:
- **Retitle to the CONCEPT, not the type.** `The interior face-flux
  cochain — :class:`WavefrontFlux`` → `The interior face-flux cochain —
  :math:`C^1_{\rm int}``. KEEP the `.. _<anchor>:` (it is cross-referenced
  from other files — grep `^\.\. _<anchor>:` across docs/+orpheus/ first;
  here api/numerics.rst + transport/fields/__init__.py both cite it).
- **A prominent succession `.. note::` immediately under the anchor**, BEFORE
  the preserved derivation. Four labelled paragraphs: (1) **Succession**
  banner (lived #205-Phase-5→S6.4, modules+tests DELETED, math REMAINS); (2)
  **Why it retired** (the verbs moved to the shared frame; the whole-N
  mesh-bound storage had no place in the per-octant walk — an octant
  transient cannot be a persistent mesh-bound field); (3) **Where it lives
  now** — name BOTH realizations with their facet: the rolling front
  (`_MovingFrontier`, per-level seed/shed — the "front itself", arguably
  TRUER) and the full-cochain oracle history (`_octant_face_cochain` /
  `_edge_outflow` — the "complete history", the fuller-view oracle); (4) the
  user's "useful concept while it lived" honoring sentence.
- **Match the tone of the ALREADY-LANDED succession notes.** The brief
  named them (transport/fields/__init__.py, docs/api/numerics.rst, the
  `_FrontierPlan` docstring, the oracle test docstring). READ them first —
  they fix the canonical naming (`_MovingFrontier` seed/shed; in-edge ι_*;
  `_edge_outflow` = ι*) and the literal-vs-role convention (they use plain
  double-backtick literals for the retired type AND the surviving private
  `_`-prefixed realizations — match that exactly).

**De-role discipline (the load-bearing correctness gate)**: an unresolvable
`:class:`/`:meth:` to a DELETED module renders as plain text WITHOUT a `-W`
warning (Cardinal-Rule-1 staleness, lesson L1) — so the build gate canNOT
catch a missed one. The ONLY gate is the explicit grep
`grep -n "wavefront_flux\.WavefrontFlux\|interior_face_space\.InteriorFaceSpace\|fields\.wavefront_flux\|spaces\.interior_face_space"`
returning EMPTY post-edit (match on the MODULE-PATH fragment, since bare
`WavefrontFlux` literals are legitimately kept as history). Every dead
`:class:`X`` / `:meth:`X.m`` → ``X`` / ``X.m`` plain literal.

**Realizations are plain literals, NOT cross-refs.** `_MovingFrontier` is in
`sweep_graph` (NOT automodule'd anywhere → no target). `FullFieldWavefront`
is in `loss_representation` (automodule'd but `:noindex:` + no
`:private-members:` → `_`-prefixed members are NOT global targets). So
ALL realization mentions are double-backtick literals — which also matches
the landed notes. Don't try to surface them with new automodule (lesson:
inconsistent half-surfacing).

**Repoint, don't just retire**: tables/notes that made a PRESENT-tense claim
about the type ("the interior face fluxes are now the typed cochain
`WavefrontFlux`") get the claim repointed to the realization ("...are the
interior 1-cochain `C^1_int` — carried since S6.4(f) by `_MovingFrontier`").
The storage-locus grid's interior-face cell becomes
"`_MovingFrontier` front / `_octant_face_cochain` history (#205 `WavefrontFlux`,
retired S6.4(f))". The ι-operator table's "Method" column becomes
"Realization verb" (`_MovingFrontier.seed`/`.shed`; `_octant_face_cochain`
in-edge / `_edge_outflow`).

**Forward-looking sections that have since landed**: the honest-scope warning
said "the transient WavefrontFlux ... is the TARGET of Phase 5b's window" and
the 1-D-scan note said the nd_foundation collapse was DEFERRED. Both LANDED.
Flip "is the target of"→"was the target of ... that window IS the
`_MovingFrontier`" and "deferred to a future session"→"that collapse has
since landed (`CumprodScan` + `_DAGWavefront` share the DAG)". A retirement
sweep is the right moment to also reconcile the now-stale future tense.

**Build gate (this worktree)**: forced `-E` cold build surfaces a WIDE
baseline that is ALL ENV-artifact + pre-existing (8 verification.rst include
CRITICALs + 8 `_generated/*.rst` InputErrors + 2 `.h5` FileNotFound + 2
homogeneous.rst plotting WARNINGs + 1 `loss_representation.py:92` title-
underline WARNING [a docstring in a NON-edited module] + 1 mesh.py paramref
ERROR = 12 warnings). NONE touch the edited files or target symbols. Gate =
the post-edit severity SET is byte-identical to pre-edit AND
`grep -iE "WARNING:|ERROR:|CRITICAL:" log | grep -iE "wavefront|interior_face|operator_algebra|discrete_ordinates"`
is EMPTY, EXIT=0, both files appear only as `reading sources`/`writing
output` progress lines.
