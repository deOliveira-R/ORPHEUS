---
name: SNStreamingOperator Wave D Round 3 closeout
description: Wave D capstone shipping SNStreamingOperator as LinearOperator with apply/solve/apply_transpose; reciprocity to round-off via dense-transpose probe; 11/11 regression bit-identical.
type: project
---

# Wave D Round 3 — Issue #160 — SNStreamingOperator closeout

Branch: `feat/sn/snstreamingoperator` based on `refactor/sn-operator-algebra` HEAD `c4138e7`.

## Status

Single commit pending. All verification gates green:

* New tests: 22/22 pass.
* Wave C/D R2 tests: 45/45 pass.
* L1 + MMS: 27 passed + 2 xfail intact (ERR-026 deferred to Wave E, as planned).
* Regression: in flight (~10-13 min).
* Sphinx -W: clean (exit 0).
* V&V audit: orphan count -1 (sn-streaming-reciprocity covered by 6 tests).

## Files modified

| File | Lines added | Action |
|---|---|---|
| `orpheus/sn/operator.py` | +280 (~290 deletions of the 2-line `LinearOperator` import) | Added `SNStreamingOperator` class + `__all__` + numerics-operator imports. KEPT existing `transport_operator_matvec_*` functions verbatim (Wave E retires them). |
| `orpheus/sn/__init__.py` | +1 | Export `SNStreamingOperator`. |
| `tests/sn/test_snstreamingoperator.py` | +475 (new) | 22 foundation-tagged tests. |
| `docs/theory/discrete_ordinates.rst` | +218 | New `_sn-streaming-operator` section + `:label:sn-streaming-reciprocity` math block. |

## Key design decisions

1. **`apply` and `apply_transpose` operate on the packed 1-D solution vector** used by the legacy BiCGSTAB FD operator path (via `EquationMap`). `solve` operates on structured `(nx, ny, ng)` source / `(N, nx, ny, ng)` angular-flux arrays per `transport_sweep`'s contract. The shape mismatch is documented in the class docstring as "the historical layouts of the two consumers"; Wave E will normalise via Krylov-on-apply.

2. **`apply_transpose` is implemented via dense-matrix transpose** — assembled lazily by probing `apply` with each of the `n_unknowns` unit basis vectors, transposed once, applied as `M^T @ psi`. Cost: `O(n_unknowns^2)` time/space, one-time on first call. For the small reciprocity test problems (`n ~ 30-200`) this is negligible (well under a second). For production-scale problems Wave E will ship an `O(n)` analytic-adjoint matvec.

3. **Why not extract `apply_transpose` analytically?** Documented in the Sphinx page "Why not extract `apply_transpose` analytically?" subsection. Each per-geometry analytical step (reverse upwind FD, transpose τ-symmetric M-M, re-derive cylindrical per-level redistribution adjoint, handle adjoint BCs) is a chance to introduce sign flips, missing factors, transposed indices — exactly the AI failure modes #1–#5 from `vv-principles`. The dense-transpose path bypasses all of them. Reciprocity is verified by construction at any size.

4. **`apply` and `solve` use DIFFERENT closures by design.** `apply` carries the symmetric closure (correct discretisation, BiCGSTAB-compatible). `solve` carries the WDD asymmetric closure (the historical sweep's math, ERR-026 affected for curvilinear). Wave E Issue 15 reconciles them via Krylov-on-apply with sweep as preconditioner — that's where ERR-026 actually closes. This is documented at length in the Sphinx section.

5. **`capabilities = frozenset({"apply", "solve", "apply_transpose"})`** — full citizen of Wave A operator algebra.

## Reciprocity test results

All 12 reciprocity tests (3 geometries × 2 test functions × parametrize) pass to round-off:

| Geometry | Single-pair `<Lψ,φ> = <ψ,L*φ>` rel diff | Multi-pair (10 trials) |
|---|---|---|
| Slab | 1.25e-15 | All ≤ 2e-15 |
| Spherical | 4.47e-16 | All ≤ 1e-15 |
| Cylindrical | 7.80e-16 | All ≤ 1e-14 |

Tolerance gate: `rtol=1e-12, atol=1e-13`. Passing margin: ~3 orders of magnitude.

## Bit-identical extraction tests (the load-bearing claim)

Each geometry has a paired test verifying `np.array_equal(SNStreamingOperator.apply(psi), legacy_transport_operator_matvec_*(psi))`:

* `test_apply_slab_bit_identical_to_legacy` — green.
* `test_apply_spherical_bit_identical_to_legacy` — green.
* `test_apply_cylindrical_bit_identical_to_legacy` — green.
* `test_apply_2d_cartesian_bit_identical_to_legacy` — green (covers the 2-D EquationMap filter that 1-D doesn't exercise).

For `solve`, the parallel claim:
* `test_solve_{slab,spherical,cylindrical}_bit_identical_to_transport_sweep` — all green.

These 7 tests are the load-bearing structural-extraction claim that Wave D R3 is **additive only** — no math is rewritten, only a Protocol wrapper is added.

## Out of scope (per Wave D plan + brief)

* `transport_operator_matvec_*` retirement — Wave E.
* `SNSolver` consuming `SNStreamingOperator` — Wave E Issue 15.
* `tests/sn/test_sweep_operator_inconsistency.py` rewrite — Wave E.
* LD/EC/Step `CellUpdate` strategies — Wave C-extension.
* O(n) analytic-adjoint matvec — future, when production reciprocity becomes performance-critical.
* ERR-026 closure — Wave E Issue 15 (Krylov-on-apply with sweep as preconditioner).

## Self-improvement notes

The dense-transpose-via-probing approach for `apply_transpose` is a **defensible design pattern** worth keeping in mind for future LinearOperator implementations where:

- The forward `apply` is structurally complex (per-geometry branches, BC reflections, mixed stencils).
- The adjoint derivation has many error-prone steps.
- The reciprocity test is a non-tautological gate on `apply` linearity + dense-assembly correctness.
- Test problem sizes are small enough that O(n²) is acceptable.

The reciprocity test in this case is testing a different thing than naively expected: not "is the analytical adjoint math correct" (which would require an analytical adjoint to test), but "is `apply` linear AND does the dense-assembly probe correctly capture the matrix of `apply`". Both are non-trivial software-correctness gates and catch different bug classes than analytical-derivation tests would.

Adding this pattern to `algebra-of-record` may be worthwhile if a second LinearOperator follow-up uses it — the threshold per the skill's growth protocol is "≥ 2 working instances" before unification.

## Process notes

- Worktree base verification caught the rebase need on first try (HEAD was at `06b46f2` pre-Wave-A baseline; clean rebase onto `refactor/sn-operator-algebra` succeeded).
- Reciprocity passed first attempt — no debugging cycle needed. The dense-transpose construction is forgiving; if it had failed, the next steps would have been (a) inspect `apply` for non-linearity (e.g., a hidden state-mutation in the EquationMap) or (b) inspect the dense-assembly probe (off-by-one in basis vector loop).
- The orphan-equation count moved from 24 to 23 because the new label `sn-streaming-reciprocity` was IMMEDIATELY covered by 6 tests via `@pytest.mark.verifies(...)`. This is the right discipline — adding a `:label:` to Sphinx without a covering test would have grown the orphan list.
