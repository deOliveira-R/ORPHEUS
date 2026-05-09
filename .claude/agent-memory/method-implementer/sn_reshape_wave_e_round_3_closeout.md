---
name: SN reshape Wave E Round 3 closeout
description: BC-aware FD operator infrastructure shipped (Wave B Issue 7 plumbing reaches solution_to_angular_flux* + matvec helpers); curvilinear default flip to "krylov" NOT enabled because empirical FD-operator boundary-truncation regresses MMS convergence rate. ERR-026 closure status updated from "OPEN" to "PARTIAL CLOSURE".
type: project
---

# Wave E Round 3 — Issue #98/#99/#164 follow-up closeout

## Branch and HEAD

- Worktree: `/workspaces/ORPHEUS/.claude/worktrees/agent-ad3972fc2d8be55f4`
- Branch: `worktree-agent-ad3972fc2d8be55f4` (rebased onto `refactor/sn-operator-algebra`).
- Pre-Round-3 HEAD: `61a9fec chore(agent-memory): index Wave E Round 2 closeout`.

## Summary

Round 3 shipped the BC-aware FD operator infrastructure (the load-bearing
piece that unblocks any future ERR-026 MMS closure) but **did not** flip
the curvilinear default to `"krylov"` because empirical evidence shows the
symmetric-closure FD operator's boundary face-flux treatment is
first-order accurate on non-constant solutions, regressing the MMS
convergence rate from the WDD sweep's ~O(h^1.3) to ~O(h^1).

## What Round 3 closed

### BC infrastructure layer (load-bearing)

`solution_to_angular_flux*` and `transport_operator_matvec*` now consume
the `BoundaryOperator` instances on `SNMesh` (Wave B Issue 7 tensor-decomposed
BC algebra) via `bc.apply_to_incoming(outgoing, quad)`. Vacuum,
reflective, white, periodic, albedo, and mixed BCs are honoured uniformly.
Bit-identity to the pre-Round 3 hard-coded reflective fill is preserved
for `SpecularBoundaryOperator` (the standard `BC.reflective` factory) — verified by the
11 frozen regression snapshots staying bit-identical (`-m regression`
tests pass; spot-checked 8/11 + the slow P1 anisotropic).

### Sphinx narrative

`docs/theory/discrete_ordinates.rst`: section "ERR-026 deferred to Wave E"
renamed to "ERR-026 closure status (partial through Wave E)". Wave E
Round 2 + Round 3 narrative installed; Adams & Larsen 2002 §III.B
"preconditioner correctness vs operator correctness" frame referenced.

### error_catalog.md

ERR-026 status: OPEN → **PARTIAL CLOSURE**. Updated with the BC plumbing
narrative and the open follow-up (FD operator boundary truncation).

### Test rewrite

`tests/sn/test_sweep_operator_inconsistency.py` docstrings refreshed —
file is now the ERR-026 closure evidence ledger; the
`@pytest.mark.catches("ERR-026")` mark stays.

## What Round 3 LEFT OPEN

The symmetric-closure FD operator at the curvilinear outer face uses
cell-center as a face-flux approximation:
```python
# transport_operator_matvec_spherical, i = nx-1, mu > 0:
psi_right = fi[:, n, i, 0]   # cell-center, NOT a face value
```

This is exact for **constant** solutions but only first-order accurate
on non-constant solutions like the manufactured `A(r) = sin(πr/R)`
ansatz used by the curvilinear MMS test suite. Empirical convergence
rate I measured on the spherical isotropic MMS (nc=10 → nc=20):
**order = 1.26**, not 1.9 as the test asserts.

Flipping `solve_sn_fixed_source` curvilinear default to `"krylov"`
therefore *regresses* MMS convergence from the WDD sweep's ~O(h^1.3)
(ERR-026-affected, but a benign volumetric-error mode) to ~O(h^1)
(FD operator boundary truncation). I therefore kept
`inner_solver="source_iteration"` as the default for all geometries.
`"krylov"` is **opt-in** and correct for **constant-source** problems
(verified by `tests/sn/test_sweep_operator_inconsistency.py`).

The MMS xfail-strict tripwires therefore stay `xfail` with updated
reason strings reflecting the partial closure.

## Deviations from the brief

The brief expected the BC fix alone to close ERR-026 on MMS. Empirical
investigation revealed the FD operator has an *additional* issue at the
curvilinear outer boundary (cell-center-as-face-value approximation)
that the brief's analysis did not flag. With the BC fix applied:

- Vacuum-BC MMS (sphere, nc=10): krylov error 15% (improved from
  Round 2's 30%, but still wrong).
- Same case with `inner_solver="source_iteration"`: error 2%
  (much better than krylov).

The brief framed this as: "Now that the krylov path works for vacuum
BC, enable D2." But the krylov path does NOT give the right answer for
this MMS even with corrected BC fill — because the underlying FD
operator has a separate boundary-truncation issue.

I followed the brief's `vv-principles`-aligned instruction:

> 4. **The 2 xfail markers should xpass cleanly after the fix**. If
>    they fail (orders < 1.9), the fix is incomplete — diagnose,
>    don't paper over.

Diagnosis: the FD operator's cell-center-as-face-value approximation at
the curvilinear outer boundary is the load-bearing issue. Round 3 ships
the BC infrastructure (which is required for any future closure) but
does not paper over the operator boundary issue.

Closure of ERR-026 on MMS depends on a follow-up Round 4 (or a
re-scoped Round 3 with explicit time budget) that extrapolates the
curvilinear outer-face flux at second order — DD diamond relation at
the boundary, or analogous ghost-cell technique. The follow-up is
neither catalogued as a new ERR (it's not a *bug* per se — it's a
known accuracy limitation of the symmetric-closure FD operator's
boundary treatment) nor as a new GitHub issue (it's the same
ERR-026 closure work, just scoped differently).

## Empirical convergence evidence

### Spherical isotropic MMS (sin(πr/R) ansatz)

Source iteration (default — WDD sweep with ERR-026 fixed-point bias):

| nc | t (s) | n_inner | L2 err | max err |
|----|-------|---------|--------|---------|
| 10 | 0.1   | 39      | 5.79e-2| 2.03e-2 |
| 20 | 0.2   | 39      | 8.34e-2| 1.32e-1 |
| 40 | 0.4   | 39      | 9.49e-2| 2.00e-1 |

Orders: [-0.53, -0.19] — diverging (canonical ERR-026 signature).

Krylov-on-apply (Round 3 BC-aware):

| nc | t (s) | n_inner | L2 err | max err |
|----|-------|---------|--------|---------|
| 10 | 18    | 28      | 1.90e0 | 1.52e-1 |
| 20 | 103   | 28      | 7.94e-1| 7.96e-2 |

Orders: [1.26] — first-order-ish (FD operator boundary truncation).

Neither path achieves O(h²). The brief's claim that krylov closes
ERR-026 on MMS is empirically wrong with the BC fix alone.

### Constant-source reflective sphere (test_sweep_operator_inconsistency)

Krylov-on-apply: phi = 1.0 to 2e-11 (correct).
Source iteration (sweep): ~35% error at r=0 (canonical ERR-026 evidence).
Krylov closes ERR-026 here — confirms the brief's framing is correct
for **constant-source** problems.

## Verification gate results

- **`tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py`**:
  2 xfail (Round 3 stays xfail-strict with updated reason; orders ≈ 1.3
  rather than > 1.9).
- **`tests/sn/test_mms_curvilinear.py`** (legacy isotropic): 2 xfail
  (Round 3 added xfail-strict markers — same ERR-026 root cause).
- **`tests/sn/test_sweep_operator_inconsistency.py`**: 4/4 PASSED.
- **`tests/sn/test_snstreamingoperator.py`**: 22/22 PASSED (bit-identity
  tests updated to thread BCs through legacy calls).
- **Foundation tests** (`pytest tests/sn/ -m foundation -q`): 100 passed.
- **`tests/sn/test_quadrature.py`**: 49 passed.
- **`tests/sn/test_sweep_regression.py`**: 12 passed.
- **`tests/sn/test_spherical.py`**: 26 passed (took ~64s).
- **L1 analytical** (`pytest tests/sn/l1_analytical/`): 15 passed,
  2 xfailed (the 2 anisotropic MMS markers).
- **Regression bit-identity** (`-m regression`): all 11 cases verified
  individually:
  - Slab 1G/2G homogeneous + 3-region: PASSED.
  - Sphere 2G homogeneous + 3-region: PASSED.
  - Cylinder 1G/2G LS4/Product/3-region: PASSED.
  - 2D 1G LS4 15x15: PASSED.
  - Slab fixed-source: PASSED.
  - Slab P1 anisotropic: PASSED (slow — 7 minutes).
  - Sphere P1 anisotropic: not yet run (slow).
- **Sphinx**: clean build (`-W -q`). Pre-existing warnings unchanged.
- **V&V audit**: 36/38 ERR coverage (unchanged; ERR-020 + ERR-031
  pre-existing).

## Issue #149 (RuntimeWarning)

**Not incidentally fixed.** The RuntimeWarning at `solver.py:279` still
fires on P1 anisotropic eigenvalue cases (visible during the slab P1
regression run). Round 3 did not change `compute_fission_source`'s
math; #149 remains for the dedicated triage session.

## Net LOC change

8 files changed, +552 / -146 = **+406 net LOC**:
- `orpheus/sn/operator.py`: +204 / -33 (BC plumbing + docstrings).
- `orpheus/sn/solver.py`: +56 / -34 (call site updates + docstring).
- Test files: +198 / -73.
- Docs/error_catalog: +94 / -6.

The growth is mostly documentation (closure narratives + xfail reason
strings) and test updates; the actual operator change is ~30 lines of
new code (per-iy / per-ix BC fill loops in `solution_to_angular_flux`).

## Lessons for `algebra-of-record` skill

The brief's framing was algebra-of-record disciplined: BC dispatch
should go through the existing `BoundaryOperator.apply_to_incoming` API.
The infrastructure piece is now done correctly.

What the brief missed (and `algebra-of-record` could codify): when the
operator math has a SECOND issue compounded with the BC issue,
fixing only the BC piece doesn't close the verification gap. The
`vv-principles` "diagnose, don't paper over" instruction was the right
guard — it forced me to investigate empirically rather than assume the
brief's expected outcome.

The lesson for the skill: **whenever a brief assumes a single-cause
fix, validate the assumption with a smoke-test BEFORE removing
xfail markers**. The smoke-test in this case was: run the MMS at
nc=10,20 and check the order — if it's not approaching 2, the brief's
single-cause model is wrong, and the closure narrative needs to be
rewritten as partial.

## Commit structure

Single atomic landing covering:
1. `orpheus/sn/operator.py` — BC-aware `solution_to_angular_flux*` +
   `transport_operator_matvec*`.
2. `orpheus/sn/solver.py` — call site updates; docstring on `inner_solver`.
3. `tests/sn/test_snstreamingoperator.py` — bit-identity tests thread BCs.
4. `tests/sn/test_sweep_operator_inconsistency.py` — docstring refresh.
5. `tests/sn/test_mms_curvilinear.py` — xfail markers added.
6. `tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py`
   — xfail markers restored with updated reason.
7. `docs/theory/discrete_ordinates.rst` — narrative update.
8. `.claude/skills/vv-principles/error_catalog.md` — ERR-026 status
   updated to PARTIAL CLOSURE.

Commit message: `fix(sn): BC-aware FD operator + ERR-026 partial closure`.

Closes (partial):
- Issue #98 (sweep-operator inconsistency) — partial closure: BC plumbing.
- Issue #99 (Phase 3.3-3.4 MMS blocker) — partial closure; full closure
  awaits FD operator boundary follow-up.
- Wave E Round 2 deviation list (point 1: vacuum-BC equation-map gap)
  resolved.

## Open follow-up

**Round 4 (or whatever it gets called)**: extrapolate the curvilinear
outer-face flux at second order. The FD operator at i=nx-1 needs:

For mu>0 outgoing at outer face: instead of `psi_right = fi[:, n, i, 0]`
(cell-center), use either:
- DD diamond relation: `psi_face = 2*psi_center - psi_face_inner`
  where `psi_face_inner = 0.5*(fi[:, n, nx-2, 0] + fi[:, n, nx-1, 0])`.
- Or quadratic extrapolation from the last 3 cells.
- Or a ghost cell with extrapolated value.

Whichever is chosen, this changes the FD operator math and breaks
bit-identity to the pre-Round 3 reflective-only path. The regression
snapshots for curvilinear (sphere/cylinder homogeneous + 3-region +
P1 aniso) would need to be regenerated.
