---
name: issue-240-phase2-step-d1-closeout
description: #240 Phase 2 Step D1 — added the generic affine-scheme outflow reconstruction op `outgoing_face_from_average` (inverse of `cell_average`) and routed the 5 inlined DD/LD direct-reconstruction sites through it (single-source the realization-leak duplication, Cardinal Rule 2). DD sites BYTE-IDENTICAL; LD sites a principled ~1-ULP re-baseline (no snapshot needed re-baselining — LD has no bit-identity snapshots, only MMS/two-paths exactness). One source-of-record SHA gate updated (deliberate edit, byte-id preserved). NEVER all tests/sn (#212 hang).
metadata:
  type: project
---

# #240 Phase 2 Step D1 closeout

**Branch** `feature/sn-space-angle-tier2` (NOT committed — user commits after
elegance + qa review). Host env, canonical `python -O`. HEAD at start `4937c3a`
(Step B records). Plan `.claude/plans/issue_240_phase2_step_d_homing.md` §D1.
Test-architect spec `.claude/agent-memory/test-architect/issue_240_phase2_step_b_sigma_leak.md`.

## What was built

ONE generic free function in `orpheus/sn/spatial/affine_closure.py`:

```python
def outgoing_face_from_average(psi_bar, face_in, w):
    return (psi_bar - (1.0 - w) * face_in) / w
```

The INVERSE of `cell_average` (`ψ̄ = (1−w)·ψ_in + w·ψ_out`). It single-sources
the per-scheme inlined outflow reconstruction (Cardinal Rule 2 — kills the
realization-leak duplication): DD's `2ψ̄ − ψ_in` is the `w=½` case; LD's
`(1+k)ψ̄ − k·ψ_in` is the `w=1/(1+k)` case. KEPT as a free function for now
(D2 homes it onto the scheme Base alongside `cell_average`/`source_emission`).

## SymPy-verified identities (`.venv/bin/python` + sympy)

All four returned 0 (exact):

1. Round-trip inverse: `outgoing_face_from_average(cell_average(in,out,w),in,w) − out = 0`.
2. DD `w=½`: `(ψ̄ − 0.5·in)/0.5 − (2ψ̄ − in) = 0`.
3. LD: `(ψ̄ − (1−w)·in)/w − [ψ̄ + (g/θ)(ψ̄−in)/D₂] = 0` with `w = 1/(1+k)`, `k = (g/θ)/D₂`.
4. Both LD forms equal `(1+k)ψ̄ − k·in`.

Numerical (numpy): DD `w=½` `array_equal` over 5000 random values (byte-identical —
`÷0.5` exact `×2`, round-to-nearest commutes with exact doubling). LD generic-vs-inline
is NOT byte-identical (different reduction tree). The synthetic worst-case ULP looks
large (1e4) but that is a relative-spacing artifact near `ψ_out ≈ 0` (catastrophic
cancellation; absolute Δ ~2e-19). In the real LD scan/kernel `ψ̄` and `ψ_in` are
correlated (smooth flux) so the absolute drift is FP-noise — confirmed by the LD MMS
O(h²) exactness + group3≡group2 two-paths gate staying green.

## Sites routed (files + lines)

| File | Line(s) (pre-edit) | Scheme | w | Bit-identity |
|------|---------------------|--------|---|--------------|
| `orpheus/sn/spatial/diamond.py` | `update` ~169 | DD | ½ | BYTE-IDENTICAL |
| `orpheus/sn/spatial/diamond.py` | `cell_kernel_batch` ~346 | DD | ½ | BYTE-IDENTICAL |
| `orpheus/sn/spatial/diamond.py` | `residual_kernel_batch` ~384 | DD | ½ | BYTE-IDENTICAL |
| `orpheus/sn/spatial/linear_discontinuous.py` | `cell_kernel_batch` ~459 | LD | 1/(1+k) | ~1-ULP principled |
| `orpheus/sn/spatial/linear_discontinuous.py` | `residual_kernel_batch` ~484 | LD | 1/(1+k) | ~1-ULP principled |
| `orpheus/sn/spatial/scan.py` | `_scanmarch_row` ~344 | DD (2-D scan) | ½ | BYTE-IDENTICAL |
| `orpheus/sn/loss_representation.py` | `_loss_action_interior` out_y ~1452 | DD (2-D scan) | ½ | BYTE-IDENTICAL |

Imports added: `from .affine_closure import outgoing_face_from_average` in
`diamond.py`, `linear_discontinuous.py`, `scan.py`; extended the existing
`affine_closure` import in `loss_representation.py`. No circular import (chain
`loss_representation → scan → affine_closure → numpy`).

## Sites DEFERRED to D5 (per brief — scan-recurrence, not direct reconstruction)

- `loss_representation.py:1431` (now ~1435) `_x_scan_faces(alpha_reflect, 2.0*psi_bar_row, …)`
  — the `2ψ̄` is the apply-direction β SCAN SOURCE for the pure-reflection x-scan
  recurrence, NOT a `ψ_out = reconstruct(ψ̄, in)` site. Folds onto the 2-D
  coefficient model in D5 (#239).
- The `alpha`/`beta` scan-recurrence terms in `_sweep_interior`/`_loss_action_interior`
  (the `α = 2 s_x/D − 1`, `β = 2(Q + s_y ψ_y_in)/D` solve coefficients) — affine-SCAN
  recurrence, D5 territory.

## Test-of-record update (deliberate, byte-id preserved)

`tests/sn/sweep/core/test_cell_kernel_batch.py::TestKernelSourceOfRecord` is a
`inspect.getsource` sha256 pin on the two DD kernel BODIES (the FP-reduction-tree
of record). The bodies' SOURCE TEXT changed (routed through the op) so the SHAs
changed — the gate explicitly instructs "update EXPECTED in THIS commit if
deliberate". Updated both:
- `cell_kernel_batch`: `645e063a…`
- `residual_kernel_batch`: `23b9787b…`
The numerical byte-identity is preserved (`w=½` op ≡ `2ψ̄ − in` bit-for-bit) —
all 503+ numerical regression snapshots stayed green with NO DriftWarning escalated.

## Gate results

- **Strict DriftWarning gate** (CORRECT path `tests.sn.regression._regression_assert.DriftWarning`):
  `tests/sn/sweep/core tests/sn/solve -W error::…DriftWarning` → **505 passed / 1 skipped /
  4 xfailed** (baseline restored after the SHA update). DD bit-identity confirmed.
- **Full route-around set** (`tests/sn/operators spatial sweep/core sweep/cartesian_2d solve`
  with the 7-red `-k` filter): **1083 passed / 6 skipped / 7 deselected / 5 xfailed**, exit 0.
  Baseline was 1080; +3 = the new `test_affine_closure.py`. The ONE DriftWarning shown
  (`vacuum_bulk_SLB_seed1`, 1 ULP within tol nulp=5) is PRE-EXISTING — confirmed by
  `git stash` (present on clean HEAD, on BOTH seeds). NOT new drift from D1.
- **LD scheme + MMS**: `test_linear_discontinuous.py` + `test_mms_ld_slab.py` →
  24 passed / 1 xfailed. LD reconstruction reduction-tree change preserves O(h²)
  linear exactness + group3≡group2 two-paths agreement.
- **New unit test** `tests/sn/spatial/test_affine_closure.py` (3 `@pytest.mark.foundation`):
  exact-inverse round-trip (FP tol, w∈{½,1,0.3,array,random}); DD `w=½` byte-identity
  (`array_equal` via `pytest.fail` — Mode-8-clean under `-O`); LD `w=1/(1+k)` algebraic
  equality (FP tol). 3 passed.
- **V&V audit**: `python -O -m tests._harness.audit` exit 0.
- **Imports**: all 4 touched modules + `affine_closure` import clean.

## No snapshot re-baselined

LD has NO bit-identity regression snapshots in the touched dirs — it is verified by
MMS O(h²) exactness + the group3≡group2 two-paths gate (both value-tolerance, not
byte-identity). So the LD ~1-ULP reduction-tree change tripped NOTHING. DD stayed
byte-identical → no DriftWarning. Only the deliberate source-of-record SHA gate
needed updating.

## Lessons / skill notes

- The `algebra-of-record` discipline applied to a CODE single-source: the SymPy
  round-trip (`outgoing_face_from_average ∘ cell_average = id`) is the foundation
  proof that the inverse op is correct; the per-scheme `w` instantiation (½ vs
  1/(1+k)) is the bifurcation. Verified BEFORE touching code (the load-bearing
  procedural rule).
- vv-principles "Bit-identity vs principled-equivalence" §: DD passes all 3 criteria
  TRIVIALLY (byte-identical via power-of-2 commute), LD passes via (1) named generic
  op intermediate `w`, (2) verified against the inlined Schur form AND the SymPy
  identity, (3) drift = reduction-tree FP-noise. No new ERR-NNN (no wrong value
  ever shipped; principled-equivalence carve like Step A/B).
- Mode-8 catch in the NEW test: my first draft used bare `assert np.array_equal(…)`
  for the DD byte-identity — INERT under `-O`. Converted to `pytest.fail`. The
  `np.testing.assert_allclose` calls fire under `-O` (function calls).
- The source-of-record SHA gate is a GENUINE source-hash exception (the kernels ARE
  the FP reduction tree of record) — a body edit that LOOKS like a refactor must
  update the SHA + re-verify byte-identity anchors in the SAME commit. The gate
  caught my edit exactly as designed.

## Next (D2..D6)

D2 homes `cell_average`/`source_emission`/`outgoing_face_from_average` onto the
DiscretizationScheme Base. D5 (#239) lifts the 2-D ScanMarch onto the coefficient
model (consumes the deferred scan-recurrence sites above). D6 archivist docs +
new Spatial × Angular tensor-product campaign issue.
