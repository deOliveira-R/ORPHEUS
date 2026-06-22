---
name: s6-9-fork-b2-default-flip
description: S6.9 Fork-B2 default flip (#222) — registry-reorder MovingFrontierWindow→ScanMarch as multi-D Cartesian default; PASS-with-conditions (stale-default doc cluster).
metadata:
  type: project
---

S6.9 Fork-B2 DEFAULT FLIP (#222, 2026-06-11) — NO retirement (user decision:
"multiple selectable methods ARE the point"). Multi-D Cartesian production default
flips MovingFrontierWindow → ScanMarch via REGISTRY REORDER only. Reviewed PASS-
with-conditions.

**Why:** the mechanism is the cleanest possible — `default_for` walks
`LOSS_REPRESENTATIONS` first-`supports`-match; reorder `(CumprodScan, ScanMarch,
MovingFrontierWindow, FullFieldWavefront)` makes ScanMarch the first multi-D match.
`__post_init__` construction guard (`type(self).supports(self.mesh)` → raise
`IncompatibleRepresentation`) makes illegal pairings unrepresentable. ZERO production
`MovingFrontierWindow(...)` constructions outside the registry; `default_for` is the
SOLE selection entry (no parallel hardcoded default map). The window stays reachable
ONLY by explicit construction (peer), unreachable via `default_for`.

**How to apply (verdicts that recur on a registry-reorder default flip):**

1. PASS the flip mechanism iff: (a) selection is first-match over a registry,
   (b) `__post_init__` guard rejects illegal pairings, (c) grep
   `<OldDefault>(` in orpheus/ shows ONLY the class def + registry (no fresh
   production constructions), (d) `default_for`/registry is the SOLE selection door.

2. The inverted-forcing test helper guard `isinstance(rep, ScanMarch) and not
   mesh.is_1d`: the `not is_1d` is REDUNDANT-BUT-LOAD-BEARING-DEFENSIVE. Redundant
   because `real_default(1-D)`→CumprodScan never yields ScanMarch at 1-D; load-
   bearing because it stops the helper constructing an illegal
   `MovingFrontierWindow(1-D)` (window.supports = `is_cartesian and ndim==2`) that
   would raise instead of producing a test result, IF the registry invariant ever
   broke. PASS — documents the asymmetry (ScanMarch.supports = `is_1d OR
   is_cartesian` ⊋ window.supports).

3. The golden REGENERATION (4× 2-D sha256) is the G5 regenerate-in-commit textbook:
   history block names date + cause (schedule change → FP-association) + class
   (principled-equiv NOT numerics, vv §bit-id cross-ref) + output-identity evidence
   (G4.a/G4.b Mode-9 FP-invariance gates + G2.c nulp oracle) + blast-radius pin
   (1-D slab hashes UNTOUCHED, CumprodScan unchanged at 1-D). PASS.

**The standing CONCERN — stale-default DOC cluster the flip didn't sweep up.**
loss_action is representation-INVARIANT (the S6.5 one-instance proof), so NO test
passes on wrong numbers — purely documentation. But ~9 test docstrings + ≥5 prod
`orpheus/` docstrings still assert "2-D Cartesian production drives
`MovingFrontierWindow.loss_action`" / "2-D Cartesian → MovingFrontierWindow", now
FALSE (2-D prod drives ScanMarch.loss_action). The sharpest: `operator.py:1403`
(`apply` docstring) is a DIRECT TWIN of the line the diff FIXED at
`operator.py:1538` (cached_property) — same "1-D→CumprodScan; 2-D Cartesian→
MovingFrontierWindow" selection fact, one fixed, the twin survived IN THE SAME FILE.
Pattern-7 remedy: the selection fact should live in ONE docstring (cached_property /
`default_for`) and the others `:meth:`-cross-reference it, not restate it. DISCRIMINATOR
for which window refs are stale vs fine: a ref that names what 2-D *production* drives
= STALE; a ref to an explicit window-vs-oracle equivalence test (e.g.
`test_2d_full_field_oracle.py` `MovingFrontierWindow.loss_action ≡
FullFieldWavefront.loss_action`, built by direct construction) = STILL CORRECT (the
window is a live peer). See [[sweep_strategy_carve]], [[s6_5_one_representation_instance]].
