---
name: ordinate-scan-pole-cell-resonance
description: ERR-054 / Issue #209 — Blelloch cumprod-divide form of ordinate_scan returns NaN when any chain a equals 0. Cylindrical pole-cell algebraic resonance at (thick=2, n=20, LS-S8, mix-A 2G) lands a exactly on 0.
metadata:
  type: project
---

# ERR-054 ordinate_scan Blelloch NaN at cylindrical pole-cell resonance

2026-05-29. Bug: `solve_sn(inner_solver='source_iteration')` on
1-D CYL reflective, thickness=2.0, n_cells=20, LS-S8, mixture A 2G
returns `keff = NaN`. Same problem via `inner_solver='krylov'`
returns analytical k_inf = 1.875 in 0.03 s.

**Why:** Lookup in `vv-principles` skill, decision matrix. The root
cause IS NOT a sweep-vs-matvec inconsistency (Phase F class); it IS
an algorithmic-stability defect in `orpheus/sn/spatial/scan.py:138`.
The Blelloch closed-form

  `cumprod_a * (psi_0 + cumsum(b / cumprod_a))`

divides by cumprod_a. The pole cell has `A_down = 0` (inner radial
face at r=0 has zero area). The cache builds `a = 2|μ|·A_total /
(2|μ|·A_down + dA_w·c_out + Σ_t·V) − 1`. At pole cell this reduces
to `a = 2|μ|·A_total / (dA_w·c_out + Σ_t·V) − 1`, which is exactly
zero ⇔ `2|μ|·A_total = dA_w·c_out + Σ_t·V`. At `μ_x = -1/√20`
(smallest |μ| in LS-S8), `dr = 0.1`, `Σ_t = 1.0` (group 1 of A),
this identity holds bit-exactly at `2π/5`. Cumprod collapses, divide
gives NaN, NaN propagates through cumsum and out. Math is fine —
the explicit `ψ[i+1] = a·ψ + b` at a=0 gives `ψ = b`, finite. Bug
is in the algorithm, not the recurrence.

**How to apply:**

1. **When seeing NaN keff at sharp single-point (mesh, quad)
   configuration:** the bug is almost certainly an algebraic
   cancellation in a closed-form algorithm, NOT a spectral-radius
   divergence. Drop FP warnings as errors (`np.seterr(divide='raise',
   invalid='raise')`) — the traceback names the exact file:line of
   the FP source in one shot. Bypasses 6 of the 8 standard probe
   cascade steps.

2. **Diagnostic cascade short-circuits when traceback points at the
   bug:** the standard 8-step cascade was a degenerate 2-step here
   (step 1 characterize → step 5 root-cause). Document that the
   `np.seterr` trick is the highest-leverage first action for any
   NaN-suspected solver bug.

3. **Existing test_ordinate_scan suite has a coverage gap:** the
   small-attenuation test uses `a ∈ [0.05, 0.2]` (positive,
   bounded). No positive contract test exists for `a ∋ 0`. When
   replacing an algorithm with a documented regime constraint,
   ALWAYS write the contract test at the regime boundary first.
   Anti-pattern #10 in vv-principles ("docstring caveat without
   enforcement") generalises.

4. **L1 standoff slowness attribution:** the briefing claimed the
   L1 standoff test is slowed by this bug. Investigation showed the
   standoff uses thickness=2.0, n_cells=40 (dr=0.05), which does
   NOT land on the resonance. If there is genuine L1 standoff
   slowness, it belongs to a different bug class — recommend a
   separate timing trace.

5. **Fix is principled but architectural:** replace `ordinate_scan`
   Blelloch closed form with a pair-monoid prefix scan. The
   existing `test_pair_monoid_associativity` already proves the
   algebra. The fix is one function in `scan.py`; no other code
   change. Fix is orthogonal to the in-flight D-H.2-C5 angular_flux
   retirement; can land independently.

**Related:**

- [[krylov_restart_truncation_bug]] — ERR-053, similar SI/Krylov
  divergence but different root cause (Krylov restart=50 silently
  truncated). Both are SI/Krylov disagreement signatures localised
  via the structural-independence cross-check.

- Phase F (ERR-026 manifestation #6 + #7) — twin-path bug (sweep
  vs matvec algebraically diverged). Different from ERR-054: that
  was an algebraic-DIFFERENCE bug between the two paths; THIS is
  an algebraic-STABILITY bug in only ONE path (SI). Krylov uses a
  different algorithm with the same operator — gives correct answer
  on the same problem.

**Filed artifacts:**

- GitHub Issue #209 (module:sn, type:bug, level:L0).
- ERR-054 entry appended to
  `.claude/skills/vv-principles/error_catalog.md`.
- Regression catcher:
  `tests/sn/test_si_cyl_20cell_nan_regression.py` (4 tests; 2 fail
  pre-fix pinning the bug, 1 passes anti-control, 1 slow
  correctness check).
- Diagnostic scripts:
  `derivations/diagnostics/diag_si_cyl_20cell_nan_step1_characterize.py`
  (6 tests pinning the sharp-resonance fingerprint),
  `derivations/diagnostics/diag_si_cyl_20cell_nan_step5_root_cause.py`
  (4 tests pinning cache-level `a=0` algebraic identity,
  ordinate_scan NaN injection, explicit-loop finiteness, Krylov
  bypass invariant).

**Methodology lesson (for future numerics-investigator sessions):**

For NaN-class solver bugs, the highest-leverage first move is to
replace the standard 8-step probe cascade's step 1 with:

```python
import numpy as np
np.seterr(divide='raise', invalid='raise')   # promote FP warnings
# ... run failing case ...
```

The Python traceback then names the exact file:line where the FP
arithmetic produced NaN. From there, the next step is to read the
suspect formula and trace which input drove the singularity. For
ERR-054 this collapsed an estimated 4-hour 8-step cascade to a
30-minute 2-step direct identification.

NOT applicable when the bug is divergence (a non-singular but
divergent iteration) — in that case the standard mesh-refinement +
spectral-radius cascade is required.

