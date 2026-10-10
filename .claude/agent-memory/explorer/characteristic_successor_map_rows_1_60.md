---
name: characteristic-successor-map-rows-1-60
description: P1 (e1b) successor census for old trajectory-resolvent rows 1-60 against the characteristic gates; which laws have no successor (2026-10-10)
metadata:
  type: project
---

Census of 2026-10-10 (branch test/characteristic-successors, HEAD f66c1c45 plus in-flight edits),
detail in `scratch/characteristic_architecture/p1_step_e/ta_e1b/map_A.tsv`. 15 COVERED, 20 PARTIAL,
25 NONE.

**Why:** step (e) deletes 22 old test files; the old rows' contracts must land somewhere first.

**How to apply:** before re-measuring, re-grep; the test-architect was adding gates during the census.
- D5 (`test_characteristic_system.py`, closed body k = k_inf + flat + group ratio) is the one broad
  successor: sphere, slab, hollow sphere, solid cylinder (slow, URRb only). No annulus, no 1G, no 4G.
- No D8 (ordering k(vac) < k(alpha) < k_inf): 0 of 266 new test functions compare eigenvalues
  across albedos (AST scan; positive controls system.py:527, :627). No D9 k-ladder except D11a
  (1G, degree axis, monotonicity only). No D10 mode-shape, no WM-72 / Ua-1-0-CY, no method-of-images
  eigen row (only the slab reading vs the E_2 image series), no Branch-1 SymPy bridge (B7).
- C9 cylinder IS covered (assembly WC1 T_w = P_ss, slow): the dependency audit's "C9 no successor"
  is wrong for the cylinder.
- Trap: a successor imported constants from a file being deleted (reading.py:60 from the old Garcia
  file); the test-architect re-pointed it to `_garcia2021.py` mid-census. [[trajectory-resolvent-retirement-audit]]
