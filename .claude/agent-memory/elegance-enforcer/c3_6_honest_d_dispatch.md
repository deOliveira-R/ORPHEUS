---
name: c3-6-honest-d-dispatch
description: C3.6 honest-d-dispatch cut (#220) verdict — the N-D campaign TAIL. _sweep_jacobi twin-kill + ScanMarch.supports narrow + streaming(axis) tuple. PASS-with-conditions (one stale-doc line).
metadata:
  type: project
---

C3.6 "honest d-dispatch" cut (branch `worktree-sn-nd-layout`, #220, the N-D
campaign TAIL after #222 closed). Continuation of [[nd-generic-carve-tuple-axis]]
(C3.0–C3.3) + [[sweep-strategy-carve]] (S1–S6.9). Three sub-steps:
`4729660` (a) d-generic schedule, `c1a72c6` (b) per-axis streaming tuple,
uncommitted (c) honest ScanMarch.supports narrow + `_sweep_2d_wavefront`→
`_sweep_jacobi` rename + twin-kill.

**VERDICT: PASS-with-conditions.** ONE CONCERN (a single stale-doc line),
inline-fixable, NON-blocking, no architecture impact. The two committed
sub-steps are clean as-is.

**Durable rulings (reinforce — RIGHT):**

1. **The `_sweep_jacobi` twin-path kill is the headline win, correctly
   single-sourced at the POLICY.** Pre-C3.6.c, all 3 multi-D `sweep` doors
   (`MovingFrontierWindow:849`/`FullFieldWavefront:1060`/`ScanMarch:1248`)
   independently inlined `_sweep_scheduled(..., schedule=SweepSchedule.jacobi(
   mesh), reflect=None, ...)` — a real 3-instance Pattern-2/7 spelling of one
   "Jacobi-schedule, no-reflect" policy. Now ONE place (`_sweep_jacobi:2308`)
   names that policy; each door supplies only its `interior` kernel +
   `moment_projection`. VERIFIED grep: only OTHER `_sweep_scheduled` prod caller
   is the G-S/SI resolvent (`solver.py:366`, genuinely different schedule).
   This is the named collapse destination my S5.1a/b/S6.4(a/b) memory tracked.
   `interior=` REQUIRED-kwarg (no default after `*`) preserves the S6.5
   Pattern-4 win (defaulted kernel = construction door outside `default_for`).

2. **`ScanMarch.supports` narrow `is_1d or (is_cartesian and ndim==2)` =
   "select narrow," the HONEST INVERSE of the S3 FullFieldWavefront widen.**
   The row-march kernels unpack d=2 today, so advertising d≥3 = Pattern-4
   over-claim lie (shape-crash vs clean `default_for` fall-through to the
   d-generic FullFieldWavefront spine). d-generalize-now correctly rejected:
   no constructible d=3 mesh, ScanMarch mesh-bound w/ `__post_init__` re-check
   → d-generic march = untestable dead code until C5. Pinned at predicate
   level (`TestD3SupportsMatrix`, duck-typed stand-ins — right level when the
   mesh is unconstructible). Widen trigger documented ("WITH the kernel
   generalization, never before") = anti-pattern-11 satisfied.

3. **`AXIS_NAMES` SSOT consolidation (a) closes a d=3 ERR-056 EXISTENCE hole
   by construction.** Retired `_AXIS_NAMES` literal (walk) + `_OUT_FACE` x/y
   hand-list dict (schedule) → one `AXIS_NAMES` (`sweep_graph.py:96`). The
   hand-list would never shed `zmax` at d=3 → reflective z-face gets no G-S
   group → silently wrong FP (shed-EXISTENCE sibling of ERR-056's shed-ORDER).
   `_octant_sweep` SOLE-projection-site + walk dropped its 2nd `signs[:ndim]`
   re-truncation (mis-sized label now fails LOUD at face zips, not silently
   re-projected). Illegal-states-unrepresentable applied to a convention table.

4. **streaming(axis) tuple (b) = bit-identical BY CONSTRUCTION.** `mu_x` IS a
   property view of `axis_cosines(0)` → same operands/expression. G-b3 EXACT
   `assert_array_equal` (NOT allclose) guards view-not-copy. Phantom `ny=1`
   `streaming_y` gone; `streaming(1)` at d=1 = honest IndexError. Accessor
   guard stays SEMANTIC `is_cartesian` (AttributeError curvilinear), not a
   substrate-None proxy (project ruling honored).

5. **`_FACE_OF` mirror-not-import balance is CORRECT (brief's Q).** Production
   derives faces (`AXIS_NAMES` f-string + min/max); test oracle hand-lists
   (`_FACE_OF` literal dict, `test_sweep_schedule.py:83`). Independent → swap
   observable (vv L11). nd test IMPORTS `_expected_outgoing` rather than
   re-hand-listing → keeps the test ORACLE single-sourced (a 2nd literal copy
   = Pattern-2 twin of the test, drift w/ no prod-vs-test signal). Mirror at
   the prod/test boundary; single-source WITHIN each side. Right structure.

6. **Duck-typed SimpleNamespace fixtures = ACCEPTABLE, not a typed-fixture gap
   (brief's Q).** Schedule layer mesh-LIGHT (`_octant_sweep`/`_outgoing_faces`
   pure on labels; `gauss_seidel` reads only ndim/octants/faces/bc_*). Typed
   d=3 fixture IMPOSSIBLE (no Mesh3D — the point) + would be Pattern-6 spec.
   Fixture docstring scopes itself (vv Mode 7: STRUCTURE not d=3 VALUE).

7. **`ScanMarch.sweep` moment ValueError still coherent post-narrow (brief's
   Q).** Fires in `if is_1d` branch ("2-D Cartesian only"); the only OTHER
   admitted mesh is `ndim==2` exactly where moment works. Message names the
   live contract truthfully. No rot.

8. **Docstring trim Cardinal-Rule-3-faithful.** The ~60-line `_sweep_2d_
   wavefront` history (L7-trap, Wave-2 per-octant-batch, PR-INDEX-5 bit-id)
   survives in `discrete_ordinates.rst` (28 anchors; L7-trap §`:1824`
   repointed to `_sweep_jacobi`). Trimmed docstring kept caller-contract
   facts only, delegated deep history to the Sphinx brain. No history lost.

**THE ONE CONCERN (condition before commit):** `loss_representations.rst:760`
G4.b test-description STALE: "deferred until a **3-D quadrature** exists" — the
EXACT phrasing the author corrected in TWO siblings (`sweep_graph.py:223`
`_FrontierPlan` + `loss_representations.rst:392`). The cut's whole premise +
the d=3 duck pins + the test-architect memo = quadrature is NOT the blocker;
only `Mesh3D` is (the quad is already 3-cosine, all 8 sign-octants). Misstated-
load-bearing-dependency = planning-error habitat (maintainer builds a redundant
quad / defers unblocked Mesh3D work). Same stale-residue class as S6.4(b) 3
CONCERNs + S6.9 stale-default cluster. Fix: name `Mesh3D`, not "3-D quadrature."

**Recurring-tell confirmation for the campaign:** rename/narrowing sub-steps
leave descriptive PROSE naming the superseded premise (here: which dependency
blocks d=3). The author corrected 2 of 3 sibling spellings in-file; the 3rd
(a test-description line in a different doc section) survived — SAME shape as
S6.9's `:1403` apply-docstring twin surviving the `:1538` fix. DISCRIMINATOR
for this carve: a doc naming the d=3 BLOCKER → "3-D quadrature" is STALE,
"Mesh3D" is correct. Grep both spellings on every C3/C5 sub-step:
`grep -rn "3-D quadrature\|3-axis" docs/theory/`.

Carve scorecard (this campaign): S1 PASS-nits, S2 PASS, S3 PASS, S4 PASS-nits,
S5.1a PASS-nit, S5.1b PASS, S6.3 PASS-nit, S6.4(a) PASS, S6.4(b) PASS-3-CONCERN,
S6.4(c/d/e/f) PASS, S6.5 PASS, S5.2 PASS, S6.9 PASS-cond, **C3.6 PASS-cond**.
NEXT after C3.6 commit = #220 close, then C5 (Mesh3D) is where the deferred
ScanMarch d≥3 row-march + the d=3 value tests land.
