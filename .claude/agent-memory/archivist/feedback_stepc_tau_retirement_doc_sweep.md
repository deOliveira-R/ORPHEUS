---
name: stepc-tau-retirement-doc-sweep
description: "#236 Phase 2 Step C — doc-only τ-reference sweep across discrete_ordinates.rst + structured_geometry.rst flipping geometry-side-τ-as-home references to closure-owned after the geometry producer was deleted. 5th/terminal sibling of the A/B1/B2/B3 carve chain; a pure stale-ref sweep + future→landed flips, NOT a stub expansion."
metadata:
  type: feedback
---

#236 Phase 2 **Step C** doc-only sweep — the terminal step of the τ
carve. Sibling of [[b3-live-tau-fold-expansion]] but a DIFFERENT task
TYPE: B3 ELEVATED a stub to rich narrative; Step C is a **stale-reference
sweep across two pages** correcting "geometry owns τ" → "closure owns τ"
after the geometry-side producer was DELETED (migrate-then-delete). No
new claim, no stub; flip references + add ONE close-out note.

**Source order (algebra-of-record for a retirement):** method-implementer
CLOSEOUT MEMO (`issue_236_phase2_stepc_tau_retirement_closeout.md`, the
migrate-then-delete narrative + the oracle re-point + mutation-RED proof)
+ the L20 3-surface AUDIT (`issue_236_stepC_dependency_audit.md`, names
the exact surgical excisions + the §6 Sphinx line list) + READ the
post-deletion CODE (`reduced_operator.py` field defs gone,
`pole_angular_closure.morel_montry_tau_per_level` is the sole producer,
`geometry.py:1330-1332` stamp, `geometry.py:1550-1555` confirming the
Wave-E shim retirement). The closeout's "OWED follow-ups: Sphinx doc-only
update" + the §6 archivist hand-off line ARE the brief.

**The triage that mattered — most τ references were NOT stale:** a τ-rich
page has ~30 `τ` hits; only the ones asserting τ STILL LIVES on the
geometry side present-tense are stale. Classify EACH:
- **STALE (fix):** factory-output lists naming `tau_mm`/`tau_mm_per_level`
  (Key Facts §2 + the cyl/sphere factory API bullets); "The code stores
  these as `SNMesh.tau_mm`" (DOUBLY stale — attribute gone AND was never
  SNMesh, audit had it on ReducedStreamingOperator); the catcher-prose
  saying the oracle reference is `st.tau_mm` (re-pointed to
  `contamination.morel_montry_weights`); "Step C CAN delete" future tense.
- **HISTORICAL/correct (KEEP):** B3 narrative "read off the closure
  rather than off the geometry-owned `StreamingTerms.tau_mm`", "no live
  reader of", "rather than from `visit.streaming_terms.tau_mm`", "the
  former `st.tau_mm` read carried" — all use "no longer / rather than /
  former", which Step C only STRENGTHENS. Do NOT touch.
- **PHYSICS-math (KEEP):** the `mm-weights`/`morel-montry-clamp` eq blocks
  (literature definitions, geometry-agnostic), the τ∈[0,1] / clamp / dome
  bullets. The DEFINITION stands regardless of where it's evaluated — say
  so when noting the producer moved.

**The catcher-prose reconciliation was the load-bearing fix.** The B3
close-out section's "production-stamp catcher" bullet said the catcher's
independent reference IS the geometry-produced `st.tau_mm`. Step C
DELETED that field and RE-POINTED the catcher onto the structurally-
independent `contamination.morel_montry_weights` (a different code path
to the same BMC-2010-Eq.43 weight, cylinder clamp replicated in the test
surrogate). The doc MUST flip this or it describes a vv-L11 reference
that no longer exists — and the new reference IS the structural-
independence story. Cite `morel_montry_weights` with `:func:` (already
used elsewhere on the page lines 999/1480 → follows page convention; it
renders plain-text since contamination isn't automodule'd, NOT a
regression).

**The close-out note shape (NEW H4 inside the B3 H3 section):** mint
`.. _sn-tau-step-c-closeout:` + a `^^^^` H4 "Step C close-out — the
geometry-side τ producer is retired". 3-numbered-step arc mirroring the
migrate-then-delete: (1) oracles re-pointed FIRST (green while geometry-τ
present = migration faithful), (2) producers excised SURGICALLY (τ blocks
interleaved with live mu_start/mu_edge/sin_theta — whole-function delete
WRONG; α-dome/redist/face-areas STAY = genuinely geometric), (3) deletion
proven inert (DriftWarning-escalated gates 0-fail + count-delta reconciles
EXACTLY to retired tests + mutation-RED re-verified the catcher isn't a
tautology). Plus a `.. note::` that the **T5 legacy `__call__`-arg
`tau_mm`** on the unbound `MorelMontryAngularSweep(sn_mesh=None)` path is
a SEPARATE surface that SURVIVES Step C (tracked #248) — pre-empts "didn't
you miss this τ?".

**Deleted-attribute role-ref hygiene.** The "can delete" paragraph had
`:attr:` refs to `ReducedStreamingOperator.tau_mm` — that attribute is
now GONE, so the role would render plain-text (no -W warning, but a
Cardinal-Rule-1 staleness bug pointing at a deleted symbol). CONVERT to
plain `` `` `` literals. Grep-gate post-edit: `:attr:` `:func:` to any
DELETED `*.tau_mm` / `*.tau_mm_per_level` must return EMPTY. (`CellVisit.tau`
SURVIVES → its `:attr:` refs stay valid.)

**structured_geometry.rst (geometry-LAYER page) — the deprecation-shim
double-staleness.** The "SNMesh as router (post-Round-1.1)" §'s present-
tense "legacy attribute names (...`tau_mm`...`tau_mm_per_level`...) survive
as @property accessors routing to self.reduced" is DOUBLY stale: (a) the
curvature shims were ALREADY retired in Wave E Round 2 (#164,
`geometry.py:1550-1555`), AND (b) Step C removed `tau_mm`/`_per_level` as
FIELDS so they can't route to `self.reduced` at all. Preserve the Wave-D
historical WHY (flip the paragraph to PAST tense "survived transitionally")
+ add a `.. note::` retirement tombstone naming BOTH retirements (#164
shim removal + #236 field removal) and stating "neither name routes
through self.reduced today; the τ entries are retained only to record the
Wave-D-era shim set". Do NOT delete the historical attribute list.

**Build gate:** `-E -W --keep-going` MAIN-checkout baseline = **1** (the
standing `Mesh1D.from_geometry :paramref:` ERROR). Post-edit = 1, EXIT=1
(the 1 warning-as-error). Count-unchanged = PASS. The 10 "no matching
equation ... ld-cartesian-1d/ld-slab — skipping" are pre-existing INFO
`verifies(...)` section-anchor registrations, NOT warnings, NOT mine.
Nexus graph.db rebuilt by the same `sphinx-build` (`_build/html/_nexus/`).

**Title-underline:** new H4 "Step C close-out — the geometry-side τ
producer is retired" = 58 code points (em-dash 1 cp / 3 bytes, τ 1 cp /
2 bytes); underline `^^^^` ≥ 58. Size with `len()` in python, NEVER `wc -c`.

**Quality self-assessment (1-5):**

| Dimension | Score | Note |
|-----------|:---:|------|
| Derivation depth | 4 | the WHY-surgical (τ interleaved with live mu_start/mu_edge → can't whole-delete) + WHY-migrate-then-delete (oracle-faithful-before-delete); not a fresh derivation (a retirement) so 4 not 5 |
| Cross-references | 5 | every :ref:/:eq:/:func:/:attr:/:class:/:mod: grep-gated; deleted-attr refs converted to literals; new anchor + #248/#164 issue links resolve |
| Numerical evidence | 5 | DriftWarning-0-fail + exact count-delta reconciliation + mutation-RED (×1.1 + c_in↔c_out) against the independent oracle |
| Failed approaches | 5 | whole-function-delete-is-WRONG (interleaved), oracle-tautology-risk (why re-point to contamination.py not closure), T5-survives-don't-conflate |
| Code traceability | 5 | every excised block cited (mu_edge loop / eta_edge loop / slab synthetic), the surviving geometry enumerated, the 3 catchers + their re-point named |
| Derivation source | 5 | closeout memo + L20 audit §6 + post-deletion code reads (the algebra-of-record for a retirement IS the deletion + the migrated oracle) |

Weakest = derivation depth (4): a field-retirement sweep has no NEW
algebra to derive — the value is the surgical-excision WHY + the
migrate-then-delete discipline, both fully captured.
