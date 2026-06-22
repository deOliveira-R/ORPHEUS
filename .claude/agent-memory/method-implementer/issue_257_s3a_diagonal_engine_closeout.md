---
name: issue-257-s3a-diagonal-engine-closeout
description: #257 S3a — generalize the dead numerics DiagonalOperator into the shared N-D diagonal-multiply broadcast ENGINE (expand_dims(coeff, broadcast_axes)*carrier); 1-D (weights,axis) is the rank-1 special case, bit-identical; the σ_t leading-axis case reduces to sigma[None]*carrier (the S3b bit-identity hinge)
metadata:
  type: project
---

# #257 S3a — DiagonalOperator → N-D broadcast engine (closeout)

Branch `feature/field-typed-operator-algebra`, HEAD `1ce727a`. NOT committed
(main agent commits after review). This is the **numerics-layer** half of the
`MultiplicationOperator` fold; S3b (transport `MultiplicationOperator(CrossSectionField)`
built ON this engine) is NOT mine.

## Scope (per the test-architect spec §1 broadcast oracle + §4b migration + §5 crosswalk)
Generalize the DEAD `DiagonalOperator` (`orpheus/numerics/operator.py`, ZERO prod
callers, 2 test callers) from "1-D weights on ONE axis" → "N-D coefficient on a
SUB-PRODUCT of the carrier's axes, broadcast over the complement". The law-suite
(§2) + resolvent (§3) are S3b — NOT touched.

## Final API
```python
DiagonalOperator(coefficient: np.ndarray,
                 broadcast_axes: tuple[int, ...] | None = None,
                 *, axis: int = 0)
```
Two construction modes, unified through ONE broadcast helper `_broadcast(ndim)`
(Pattern 2, single source of truth):

- **1-D special case** (`broadcast_axes is None`): `coefficient` MUST be 1-D;
  `axis` selects the single carrier axis it occupies; broadcast is
  RANK-AGNOSTIC (carrier rank read at apply-time → the SAME operator acts on
  1-D / 2-D / N-D carriers). `_broadcast` returns the historical
  `(1,…,N,…,1)` reshape — **byte-identical to the old `_reshape`**. This is the
  path the `TensorProductOperator` Kronecker compositions ride.
- **broadcast_axes mode** (explicit tuple): `coefficient` is N-D;
  `_broadcast` returns `np.expand_dims(coefficient, broadcast_axes)`. The
  coefficient's rank MUST equal `carrier.ndim - len(broadcast_axes)`, its axes
  map IN ORDER onto the carrier axes NOT in `broadcast_axes`.

Body (the load-bearing 2 lines):
```python
# _broadcast(ndim)
if self.broadcast_axes is None:
    shape = [1]*ndim; shape[self.axis] = -1
    return self.coefficient.reshape(shape)          # 1-D legacy, bit-id
return np.expand_dims(self.coefficient, self.broadcast_axes)  # N-D engine
# apply: return self._broadcast(x.ndim) * x_arr     (after _check_shape)
```
The **σ_t case S3b needs**: `DiagonalOperator(sigma, broadcast_axes=(0,))` with
`sigma:(ng,*spatial)` on a `(N,ng,*spatial)` carrier reduces EXACTLY to
`sigma[None]*carrier` (verified 0 ULP) — THE bit-identity hinge for S3b's
`CollisionOperator → MultiplicationOperator` promotion.

PRESERVED: `from_measure(measure, axis=0)` (1-D → engine); self-adjoint
(`apply_transpose==apply`, real coeff); `CAP_SOLVE iff np.all(coeff!=0.0)`
(eager, honest — Pattern 4; the legacy bare-σ collision path had NO σ=0 gate →
silent IEEE NaN, this is the chance to add it). The name `DiagonalOperator`
KEPT.

## Design decisions (elegance)
- `axis` is **keyword-only** (`*, axis`) so the positional slot is the
  N-D-first `(coefficient, broadcast_axes)`. The 2 test files + `from_measure`
  ALL call `axis=` by keyword → migration clean.
- `self.axis` stored as a plain `int` ALWAYS (default 0, harmless/ignored in
  broadcast mode) — NOT `int | None`. This was the pyright fix: my first draft
  set `self.axis = None` in broadcast mode → `shape[self.axis]` flagged
  `int | None` not indexable. Keeping it int is both clean AND correct (the
  attribute is consulted only when `broadcast_axes is None`).
- `weights` is a **@property** (not a stored attr) returning `self.coefficient`
  in 1-D mode, **raising AttributeError in broadcast mode** (Pattern 4 — reading
  a 1-D-named handle off a multi-axis coefficient is an illegal state; do NOT
  silently return an N-D array under the `weights` name). This kept the property
  cleanly typed `np.ndarray` → killed the `NDArray | None` errors the
  `np.outer(D.weights, …)` Kronecker tests would otherwise inherit.
- `_check_shape` rank-discriminates: 1-D mode checks `x.shape[axis]==N`; N-D
  mode checks `coefficient.ndim == x.ndim - len(broadcast_axes)` with a message
  naming the expected rank.

## Tests (spec §1 + §4b)
File `tests/numerics/test_diagonal_operator.py` migrated:
- ALL legacy 1-D tests (`TestDiagonalApplyShape/SelfAdjoint/CapabilitiesAndSolve/
  Composition/FromMeasure`) pass UNCHANGED against the generalized engine (the
  1-D path is bit-identical) — only the `solve` ValueError-match string updated
  (`"non-zero weights"`→`"non-zero coefficient"`, message reworded for the
  coefficient vocabulary).
- `TestDiagonalRejectsBadInput` migrated to the NEW contract: 2-D coeff WITHOUT
  `broadcast_axes` rejected (`"1-D special case"`), coefficient-rank-must-match-
  complement-rank rejected (`"rank-1 coefficient"` on a wrong-rank carrier),
  distinct-broadcast-axes rejected. (Replaced the old `test_2d_weights_rejected`
  which asserted "must be 1-D".)
- NEW `@pytest.mark.foundation TestDiagonalBroadcastOracle` (spec §1):
  (a) `test_leading_axis_equals_sigma_None_times_carrier` — the σ_t form on a
  carrier `(N=6, ng=2, nx=5, ny=3)` with **nx≠ny** (the vv-mode-#2 variable-swap
  discriminator) `≡ sigma[None]*carrier`, `assert_array_equal` 0 ULP;
  (b) `test_one_d_axis_mode_equals_legacy_reshape` — 1-D mode bit-identical to an
  INDEPENDENT hand-built `_reshape` across 4 axis/rank combos;
  (c) `test_non_leading_broadcast_axis` — engine is NOT hardcoded to axis=0
  (`(N,nx)` coeff on the MIDDLE axis of `(N,ng,nx)` ≡ `coeff[:,None,:]*carrier`).
  Docstring DOCUMENTS that 0-ULP bit-identity here is **inherited verification,
  not coincidence** (reduction depth 1, pure broadcast-multiply, NO sum) so a
  reviewer must NOT silently loosen to allclose/nulp — a non-bit-id result would
  be acceptable ONLY from a more-principled construction, which a pure broadcast
  cannot produce.

File `tests/numerics/test_tensor_product_operator.py`: **byte-UNCHANGED** (git
diff empty). All `DiagonalOperator(w, axis=k) & …` Kronecker compositions are
rank-1 instances of the engine; `.weights`/`.axis`/`T.ops` reads all preserved.

## Verification (paste-back, lessons-L12)
- `.venv/bin/python -O -m pytest -q tests/numerics/test_diagonal_operator.py
  tests/numerics/test_tensor_product_operator.py -p no:cacheprovider`:
  **`49 passed, 1 warning in 0.30s`** (the warning = the standard `-O` "assert
  statements not executed" config-warning; every test uses `np.testing.*` /
  `pytest.raises` → Mode-8 SAFE, NO inert bare-assert).
- Full `tests/numerics/`: **657 passed** (export + `__init__` import path clean).
- `npx pyright orpheus/numerics/operator.py tests/numerics/test_diagonal_operator.py
  tests/numerics/test_tensor_product_operator.py`: **22 errors, 0 warnings** —
  ALL 22 PRE-EXISTING #226 rooting noise, ZERO net-new, NO `# type: ignore`
  added. Proof: my first draft showed 35 (the +13 = `self.axis`/`self.weights`
  None-union + `bcast` possibly-unbound I created); after the int-axis +
  weights-property fix → 22. The 22 break down as: operator.py 631/744/875
  (`block_role` attr-assign covariance on `_AdjointOperator`/`OperatorSum`/
  `ScaledOperator`), 1835 (`spla.LinearOperator` scipy-stub overload in
  `as_scipy_linop`) — all OUTSIDE the DiagonalOperator class; test files
  (`test_tensor` byte-unchanged → 100% pre-existing) the `+`/`@`/`&` operator-
  overload `reportOperatorIssue` + `np.outer`/`assert_array_equal` arg-type
  diagnostics, the standard pyright-doesn't-see-the-LinearOperatorMixin-dunders
  noise.

## MUTATION-verified (vv mode #2 catch)
Out-of-band probe: correct engine `== sigma[None]*psi` True; a WRONG
`broadcast_axes=(3,)` (broadcast over a spatial axis instead of the leading
ordinate) RAISES ValueError on the `(6,2,5,3)` nx≠ny carrier — the rank/shape
contract makes the transposed-axis mistake non-silent (an equal-nx-ny carrier
would MASK a leading↔spatial swap; nx≠ny is what exposes it).

## Crosswalk (spec §5, Pattern 7 — the ONE bridge row)
| layer | input | internal | output |
|---|---|---|---|
| legacy collision | `σ:(ng,*sp)` + `ψ:(N,ng,*sp)` | `σ[None]*ψ` leading-axis | flux-shaped |
| **engine (this)** | `coeff` + `broadcast_axes` | `expand_dims(coeff,bcast)*carrier` | bare ndarray |
| bridge (S3b owns) | — | transport `MultiplicationOperator` DELEGATES broadcast to engine on `ψ.bulk.values`, `broadcast_axes` DERIVED from the carrier's ordinate-axis position | applies typed codomain (AngularSourceSink) |
The §5 "broadcast_axes-must-be-DERIVED" point: S3b must compute `broadcast_axes`
from the carrier layout (the ordinate axis position), NOT hardcode `(0,)` — the
engine is layout-agnostic, the transport wrapper supplies the axis.

## NO algebra-of-record manifest
Pure numerics-primitive generalization (a broadcast widening), NOT a new
reference solver — no Branch-1 SymPy / L1 cross-check / Sphinx theory stub /
archivist dispatch owed. The "reference" for the broadcast oracle is
`np.expand_dims`/`sigma[None]` itself (trusted-library, structurally
independent of the engine's `_broadcast`). NO new ERR (proactive generalization,
no bug caught; next free ERR-063 reserved per #251 note).

## Deliverables
- `orpheus/numerics/operator.py` — generalized `DiagonalOperator` (the engine).
- `tests/numerics/test_diagonal_operator.py` — migrated + foundation oracle.
- `tests/numerics/test_tensor_product_operator.py` — UNCHANGED (rank-1 survives).

## LESSON (coding-elegance / pyright-clean generalization)
When generalizing a primitive from "1-D-on-one-axis" to "N-D-on-a-sub-product",
the lazy rank-agnostic 1-D mode (carrier rank read at apply-time) and the
explicit-axes N-D mode are genuinely DIFFERENT broadcast policies — unify them
through ONE `_broadcast(ndim)` helper (Pattern 2) but do NOT try to collapse the
1-D `axis` into a degenerate `broadcast_axes` (the 1-D mode must stay
rank-agnostic to preserve byte-identity for unknown-rank carriers). Two pyright
traps when adding a second construction mode: (1) do NOT store a
mode-discriminated attr as `T | None` if it's indexed in the OTHER mode — keep
it the always-valid type (here `axis: int`, default 0, ignored when unused);
(2) a legacy attribute that is meaningful in only ONE mode is best a `@property`
that RAISES in the other mode (Pattern 4 illegal-state) rather than a stored
`ndarray | None` field — the property stays cleanly typed and the N-D-coeff
consumers can't read a 1-D-named handle by accident.
