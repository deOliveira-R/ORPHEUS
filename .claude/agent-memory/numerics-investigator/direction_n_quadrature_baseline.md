---
name: Direction N F.4 quadrature baseline (Issue #123)
description: Empirical F.4 structural-floor baseline on the canonical 6-point grid. All six points UNRESOLVED at the 120 s/run devcontainer budget; RICH vs RICH+panels trajectories pinned into tests/cp/test_peierls_rank_n_protocol.py as the L17/L19 evidence. 4/6 points confirm L17 sign-flip or magnitude-growth under one-panel refinement.
type: project
---

# Direction N — F.4 structural-floor baseline (Issue #123, 2026-04-22)

## STATUS: BASELINE COMPLETE / ALL POINTS UNRESOLVED AT 120 S BUDGET

The L19-compliant protocol helper (`assert_rank_n_structural_win`) and the
per-point F.4 sign-stability scanner both ship in:

- `derivations/diagnostics/diag_f4_structural_floor_baseline.py`
- `tests/cp/test_peierls_rank_n_protocol.py` (14 unit tests of helper + 12
  parametrized F.4 baseline tests; 8 of the 14 unit tests are fast, the
  12 baseline tests are all `@pytest.mark.slow`)

The 6-point scan reveals that **F.4's convergence below 0.005 % is out of
reach of the 120 s/run budget on the devcontainer**. RICH (4, 8, 64) takes
44–67 s per point; RICH+panels (5, 8, 64) takes 78–119 s; RICH+pp (5, 10, 64)
and ULTRA (5, 10, 96) both exceed 120 s at every point. The pinned data is
the RICH vs RICH+panels pair.

## Per-point trajectory (2026-04-22 devcontainer scan)

| sigma_t*R | rho | RICH signed | RICH+panels signed | trajectory       | status |
|-----------|-----|-------------|---------------------|------------------|--------|
| 5.0       | 0.3 | +0.057825 % | +0.032331 %         | monotone_down    | unresolved, > 0.005 % |
| 5.0       | 0.5 | +0.297309 % | +0.322013 %         | magnitude_grew   | unresolved, quadrature regime |
| 10.0      | 0.3 | +0.003285 % | −0.008221 %         | **sign_flip**    | unresolved, L17 canonical |
| 10.0      | 0.5 | +0.016683 % | +0.007990 %         | monotone_down    | unresolved, close to threshold |
| 20.0      | 0.3 | −0.005982 % | −0.007481 %         | magnitude_grew   | unresolved, marginal growth |
| 20.0      | 0.5 | +0.005430 % | −0.003942 %         | **sign_flip**    | unresolved, L17 canonical |

Wall times per run, F.4 RICH / RICH+panels respectively:
44/96, 44/78, 48/100, 47/81, 67/119, 58/93 s. No run exceeded the 120 s
budget; the scan's total wall was 2323 s (the subprocess timeouts dominate —
4 quads × 6 points × up to 120 s = 2880 s theoretical max; 8 of 12 slow-quad
cells hit the ceiling).

**4/6 points confirm L17 directly**: two sign flips reproduce the L17
signatures reported in the research log (10.0, 0.3 and 20.0, 0.5), and
two magnitude growths (5.0, 0.5 and 20.0, 0.3) show F.4's RICH residual
is dominated by quadrature noise at those points. Only (5.0, 0.3) and
(10.0, 0.5) are on a clean trailing-edge convergence but still above the
0.005 % resolution threshold.

## Hardware cost to resolve each point (extrapolation)

At the anchor (10.0, 0.3) the historical record (L17, pre-kill probe from
Direction C) puts F.4 ULTRA = 0.0010 % — a 3× improvement over RICH. On
this devcontainer ULTRA on the anchor was previously timed at ~150–200 s
(extrapolated from the 100 s RICH+panels run at this point). A 180 s/run
budget would pin the anchor. A 300 s/run budget would pin all 6 points.

For the sign-flip points, nothing below ULTRA is decisive — RICH+panels
revealed the flip but didn't go past it. RICH+pp would probably reveal
whether RICH+panels is already on the trailing edge of the flip or in a
second crossing zone; ULTRA closes the question.

## Helper signature (in `tests/cp/test_peierls_rank_n_protocol.py`)

```python
def assert_rank_n_structural_win(
    closure_fn: Callable[[float, float, dict], float],  # returns k_eff
    f4_fn:      Callable[[float, float, dict], float],  # returns k_eff
    point: tuple[float, float],
    quads: Sequence[tuple[str, dict]],  # >= 2 (label, quad_dict)
) -> StabilityReport
```

Raises `AssertionError` if ANY of five L19 conditions fail:
1. `len(quads) >= 2`
2. `|closure| < |F.4|` at every quadrature (STRICT beat, not tie)
3. closure signed err does not flip sign across quads
4. closure `|err|` is monotone non-increasing across quads
5. F.4 `|err|` itself is monotone non-increasing AND F.4 signed err is
   sign-stable (otherwise the reference is in F.4's own quadrature-noise
   regime, L17 signature — unverifiable at these quads)

## Is the baseline crisp enough to anchor future rank-N comparisons?

**Not on this hardware.** At the 120 s/run budget the pinned per-point
data captures only RICH and RICH+panels. That is exactly the regime L17
identified as NOISE-LIMITED for F.4 at sigma_t * R >= 10. A future closure
claim cannot pass `assert_rank_n_structural_win` at the 6-point grid
with only these two quadratures — condition (5) (F.4 sign-stable + monotone)
fails at 4/6 points.

The protocol is crisp (all 14 helper unit tests pass, synthetic failure
modes are all caught). The baseline is partially crisp: it PROVES L17
on the devcontainer at RICH vs RICH+panels resolution (the sign flips
reproduce exactly), but it CANNOT YET supply a sub-0.005 % reference
quadrature for any of the 6 points. Any serious rank-N comparison must
either be run on faster hardware or use a >= 180 s/run budget.

## Files shipped

- `derivations/diagnostics/diag_f4_structural_floor_baseline.py` — scanner
- `tests/cp/test_peierls_rank_n_protocol.py` — helper + unit tests + pin
- This memo: `.claude/agent-memory/numerics-investigator/direction_n_quadrature_baseline.md`

## Follow-ups (not priorities now)

1. Re-run the scan with a 300 s/run budget on any faster machine; fill in
   the `ref`/`prev` fields in `F4_REFERENCE_BASELINE`; flip the per-point
   `status` from the unresolved reason string to `"resolved"`. The
   `test_f4_is_sign_stable_at_its_reference_quadrature` parametrize will
   then run instead of skip.
2. Add the `:label: peierls-rank-n-stability` block to
   `docs/theory/peierls_unified.rst` next to `peierls-rank-n-bc-closure`
   (marker-coverage gap flagged by the current `@pytest.mark.verifies`).
3. Issue #121 (PCA falsification) is ready to close referencing this memo.
4. Issue #123 itself stays open until (1) lands; the test file is already
   the L19 enforcement point.
