---
name: w2-pole-cell-characterization-gate
description: W2 (#233) permanent L∞/per-cell-rate CHARACTERIZATION gate for the curvilinear DD central-cell (r→0) spatial closure — pins the WONTFIX first-order-at-pole limitation in test form WITHOUT calcifying it. GREEN under -O 2026-06-13. Zero solver code.
metadata:
  type: project
---

W2 of the curvilinear-aniso program (branch `fix/curvilinear-aniso-pole-and-clamp`,
W1=`b2d8a6d`). File: `tests/sn/verification/mms/test_curvilinear_pole_cell_characterization.py`
(4 tests, `@l1`, GREEN under `.venv/bin/python -O`). DISTINCT from W1 (angular
τ-clamp) and #229 (angular half-angle-thread interpolation floor) — this is the
**spatial** central-cell closure `dd-curvilinear-scalar` (Hébert §3.9.4 / Stacey
§9.9 plain diamond; A(0)=0 ⟹ diamond over-predicts pole outer face +50% mesh-indep;
balance itself inconsistent at A_in=0). WONTFIX for DD; genuine fix = LD #6 / nodal
#158/#233.

⭐ CHARACTERIZATION-GATE DESIGN (the reusable part — pin TRUE + protect floor,
do NOT calcify the limitation):
- **GUARANTEE tests** (carry `verifies(dd-curvilinear-scalar, transport-{spherical,
  cylindrical}, sn-mms-*)`): global volume-weighted L2 `√Σ V·diff²` is O(h²)
  (`orders > 1.9`). SPHERE asserted under BOTH refs — midpoint AND Hébert-3.430
  **shell-volume-average** (the principled ref: the curvilinear DD unknown IS the
  shell avg `(4π/V)∫r²φ dr`, so comparing to it compares the unknown to its
  definition; agreement-on-order across two structurally-different refs proves the
  L2 order is REAL not a midpoint artifact). Shell-avg built from
  `scipy.integrate.quad` (trusted-lib, structurally indep of solver). CYLINDER
  global L2 O(h²) (midpoint).
- **CHARACTERIZATION tests** (NO `verifies` — pin a LIMITATION not a correctness
  claim): pole L∞ order **bounded BELOW only** (`> 0.8`, "≥ first-order, does not
  regress") + pole IS the L∞-dominant cell (`fraction > 0.99`) + interior clean
  O(h²) (`> 1.8`). ⭐⭐ NO UPPER BOUND on the pole order → a future LD/nodal fix that
  lifts pole to O(h²) keeps the gate GREEN (2.0 > 0.8). This is vv anti-pattern
  #5/#17: a characterization gate pins what is TRUE + the regression floor WITHOUT
  blocking a legitimate improvement.

MEASURED (reproduced exactly from diag_14/diag_30, sphere n_ord=16 / cyl n_mu4 n_phi8,
ladder [40,80,160,320], in the docstrings per vv discipline):
- sphere L2 midpt orders 2.01×3; L2 shell-avg 2.00×3; L∞(pole) 0.91/0.95/0.97;
  interior(max r/R>0.1) 1.84/1.92/1.96; pole-fraction 1.00 every mesh.
- cyl pole-vs-midpt 1.94/1.97/1.98 (ACCIDENTALLY O(h²) — r dr weight puts
  vol-centroid where diamond ½A(h)≈midpoint); pole-vs-volavg 0.99/0.99/1.00 (the
  REAL O(h) defect — SAME diamond inconsistency as sphere, MASKED by midpoint).
  Cyl test pins BOTH halves (`ord_mid>1.8` AND `ord_va>0.8` AND `median(ord_va)<1.5`)
  so the "cyl pole is clean O(h²)" misreading can't creep back.

⭐ THE HEADLINE the docstring must carry: L∞(pole)=O(h) COEXISTS with L2(global)=O(h²)
because the pole is ONE cell of V~h³ → √V~h^1.5 dilution → ~h^2.5 contribution to L2
→ subdominant → INVISIBLE to the production gate + keff. This is WHY #233 needed an
L∞/per-cell probe to surface.

Mode-8 VERIFIED LIVE this session: pytest rewrites asserts in collected `test_*`
modules so bare `assert` FIRES under `-O` (the PytestConfigWarning only covers
NON-test-module asserts). Probe: flipped `pole_orders > 0.8`→`> 1.5` → test FAILED
under `-O` exposing the real `[0.907,0.947,0.970]`. All assertions kept in collected
test fns; the shell-avg helper computes refs only (no asserts).

DECISIONS: marker `@l1` (MMS convergence-order = math claim); ladder [40,80,160,320]
(decisive AND ~0.3s sphere / 2.8s cyl → NO `@slow`); `verifies` on GUARANTEE tests
only (label `dd-curvilinear-scalar` EXISTS at discrete_ordinates.rst:1292 = the exact
cell-update); NO `@catches` yet (pole-cell ERR-NNN minted in W5 — docstring notes W5
adds the link). FINDING SOUND: diag_14/diag_30 reproduced bit-for-bit; cylinder claim
CONFIRMED (not papered over). Extends [[w1-sphere-clamp-removal-verification]] +
[[curvilinear-aniso-229-9-verification]].
