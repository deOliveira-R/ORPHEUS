---
name: issue-257-s3a-diagonal-engine
description: PASS-WITH-NITS verdict on #257 S3a (DiagonalOperator generalized from 1-D-weights-on-one-axis into the N-D diagonal-multiply broadcast ENGINE; the numerics layer the transport MultiplicationOperator(σ_t) will delegate to in S3b) on branch feature/field-typed-operator-algebra.
metadata:
  type: project
---

# #257 S3a — DiagonalOperator → N-D broadcast engine (numerics layer of M[σ_t])

PASS-WITH-NITS, `feature/field-typed-operator-algebra`, HEAD `1ce727a`, pre-commit.
DEAD-CODE reshape (ZERO production `DiagonalOperator(...)` construction sites —
grep-verified; all construction is in tests). 1-D path byte-identical;
`from_measure`/self-adjoint/`CAP_SOLVE`-iff-nonzero survive; ZERO new `# type: ignore`
(region 1436-1657 pyright-clean; the 2 operator.py errors @875 block_role-ClassVar +
@1835 scipy-LinearOperator-call are PRE-EXISTING, outside region). 49/49 green
(test_diagonal + test_tensor_product, the latter byte-unchanged).

## The four judgment calls

1. **Dual-mode API (`broadcast_axes=None` 1-D vs explicit N-D) = KEEP, NOT the
   boolean-flag anti-pattern (#3).** The discriminator: `broadcast_axes` is not a
   bool toggling two behaviours of ONE operation — it is the DATA the N-D broadcast
   needs (which carrier axes the coeff is constant over). `None` is genuine absence
   (no axes named → the rank-agnostic 1-D reshape, which CANNOT name its complement
   because the carrier rank is unknown at construction). The two modes are not two
   operators in a trenchcoat; they are one operation (`expand_dims(c, bcast)*x`)
   whose 1-D case is forced to defer axis-placement to apply-time. The "always
   broadcast_axes + from_axis constructor" alternative was WEIGHED and is WORSE: it
   would force every 1-D caller to know the carrier rank at construction (defeating
   the rank-agnostic property the SN angular weight matrix relies on — same `W` acts
   on 1-D/2-D/N-D). LESSON: `X | None` where None=structural-absence-of-data is NOT
   the flag anti-pattern; the flag anti-pattern is `bool` selecting between two
   complete behaviours. NIT: the lingering `*, axis` is dead in broadcast mode —
   acceptable (well-typed default-0, documented "ignored"), but see NIT-2.

2. **`weights` @property RAISING in broadcast mode = KEEP (Pattern 4, textbook).**
   Mirrors the S1 Spectrumon-simplex move: an N-D coefficient under a 1-D name
   (`weights`) is the illegal state; returning it would be a units/rank lie. Raising
   AttributeError (not returning `.coefficient`) is correct — `weights` is RETAINED
   only as the back-compat 1-D alias (the historical name), and there are ZERO
   production readers of it on a DiagonalOperator instance (grep-verified: all
   `.weights` prod hits are `quad.weights`/`measure.weights`/`emit.weights`/
   `context.weights`, NONE a DiagonalOperator). So `weights` is a pure tests +
   ergonomics alias. Could be retired entirely (no prod consumer), but keeping it as
   a raising-in-illegal-mode alias is the conservative dead-code-reshape call and
   reads honestly.

3. **`_broadcast`/`_check_shape` helpers = KEEP, NOT over-factored.** `_broadcast`
   is the genuine SSOT — called by BOTH `apply` and `solve` (the `*` and the `/`),
   single-sources the two-mode reshape. `_check_shape` likewise single-sources the
   validation across apply+solve (the OLD code DUPLICATED the axis-size check inline
   in both apply and solve — the diff REMOVES that twin → Pattern-2 win). For a
   2-method 2-mode operator this is exactly right-sized, not ceremony.

4. **Naming `DiagonalOperator` for an N-D-coefficient multiply = KEEP.** Still
   diagonal in the full product basis (pointwise multiply = diagonal matrix); rank
   of the coefficient is just how many axes the diagonal varies over. The
   numerics/transport split reads RIGHT: `DiagonalOperator` (numerics) = the
   basis-agnostic broadcast engine; S3b `MultiplicationOperator(σ_t)` (transport) =
   the physics-named promotion that DELEGATES to it. Two names = two layers, not two
   names for one concept (the Cardinal-2 trap) — the transport name carries the
   physical meaning (collision/multiplier algebra) the numerics engine must not know.

## The one real nit (DO-NOW, Cardinal-Rule-3 doc accuracy)

⭐ **NIT-1 (the docstring asserts a coupling that does NOT exist).** Class docstring
:1506-1510 + the `weights` property docstring :1586-1587 both claim "the
tensor-product Kronecker compositions (`TensorProductOperator`) ... read it
[`weights`/`axis`]". FALSE — VERIFIED: `TensorProductOperator.apply/.solve` route
PURELY through `op.apply(out)`/`op.solve(out)` (operator.py:1312-1340); they read
NEITHER `.weights` NOR `.axis`. The `.axis` is consumed INTERNALLY by
`DiagonalOperator._broadcast`, never by the composition. This is exactly the
[[issue-236-phase2-stepc-tau-retirement]] STANDING-TELL class: a rationale-comment
asserting a load-bearing coupling the code does not depend on = a bug habitat (a
future maintainer "preserves `weights` for TensorProductOperator", or conversely
deletes the `axis` attr believing TP reads it, and the false claim misdirects).
Fix: reword to "the 1-D `weights`/`axis` attributes are retained as the back-compat
alias for `from_measure` ergonomics and existing 1-D call sites" — drop the
TensorProductOperator coupling claim entirely (TP only needs `.apply`).

## Minor (file/note, no patch)
- The `numerics/vector.py:46` note calls DiagonalOperator a "flat axis-primitive
  that acts on one tagged axis" — now outdated (acts on N-D coeff too). NOT in this
  diff; tiny. Sweep when S3b lands (the note will want the M[σ_t] mention anyway).
- `_check_shape` 1-D branch validates only `x.shape[axis]==coeff.shape[0]`; the
  broadcast branch validates only RANK (not the per-axis sizes of the complement —
  numpy's `*` will catch a size mismatch at the multiply with its own error). Fine
  for dead code; if S3b surfaces a confusing numpy broadcast error, add a
  complement-shape check. Not now.

## Standing tells reinforced
- `X | None` param where None = absence-of-the-data-the-other-mode-needs is NOT the
  boolean-flag anti-pattern; distinguish from `bool` selecting two complete behaviours.
- Property-raises-in-illegal-mode (Pattern 4) is the right retirement-conservative
  move for a legacy alias name with zero prod readers — keep it raising, don't
  silently widen it to return the N-D array.
- DOCSTRING-asserts-a-coupling-the-code-lacks: same habitat as the rationale-comment-
  asserts-false-ordering tell — grep the named collaborator (here TensorProductOperator)
  to verify the claimed read actually happens before trusting the doc.
