---
name: issue-257-s3b-multiplication-operator
description: S3b — CollisionOperator IS MultiplicationOperator (C=M[σ_t]); PASS-WITH-NITS, dead-import + engine-rebuild nits
metadata:
  type: project
---

# #257 S3b — `C = M[σ_t]` operator-as-promotion (MultiplicationOperator)

**PASS-WITH-NITS** (`feature/field-typed-operator-algebra`, HEAD `c1da42d`, pre-commit,
uncommitted working tree). The §5.7 promotion: `CollisionOperator(MultiplicationOperator)`
is the field promoted to a diagonal operator. NEW `orpheus/transport/multiplication_operator.py`
(`MultiplicationOperator(LinearOperatorMixin["TimedFullField"])`, `@dataclass(eq=False)`);
`CollisionOperator` retyped to a THIN subclass (apply/solve/apply_transpose bodies DELETED,
now inherited). 11 new tests ✓; 63 resolvent-leg gates (`test_invertible_operator` +
`test_kinf_homogeneous`) stay ✓ → carve bit-identical.

## What landed (the consumer step of S3a)
S3a built the N-D `DiagonalOperator(c, broadcast_axes=(0,))` engine; S3b is the transport-layer
consumer that delegates `f·ψ`/`q/f` to it and owns ONLY the typed codomain. `apply`→
`AngularSourceSink` (collision RATE density; `flux × 1/cm = rate`), `solve`→`AngularFlux`,
`apply_transpose==apply` (real diagonal self-adjoint, `M.H==M`). 3 solver.py sites rewired
`total_cross_section`→`total_cross_section_field` (the S2 zero-copy typed lens → bit-id).

## 4 axes ALL PASS
- **P2 single-source**: multiply lives ONCE in `DiagonalOperator.apply/solve` (operator.py:1639-1656).
  Deleting CollisionOperator's old `sigma[None]*psi` body removed REAL duplication (both σ_t and σ_r
  collision now route the same engine). MRO-verified: `CollisionOperator` inherits apply/solve/
  apply_transpose verbatim; subclass adds only `__init__`+`.sigma`+`__add__`. No twin survives.
- **P1 reads-like-math**: `C=M[σ_t]`, codomain-asymmetry reads like units flow; docstring law-suite
  (M[1]=I/M[0]=0/lin/homo/self-adj/spec=ess-range) faithfully mirrored test-for-test. Homomorphism
  test correctly DROPS to bare engine on raw product (σ·σ' is cm⁻², a numerics fact NOT a XS field).
- **P4 illegal-states**: `__post_init__` builds an engine, inherits `engine.capabilities` WHOLESALE
  → "CAP_SOLVE iff min|σ|>0" single-sourced at `DiagonalOperator.__init__:1563`, NOT a copied
  `np.all(σ!=0)`. Behavioural STRENGTHENING over legacy (advertised CAP_SOLVE unconditionally → silent
  NaN). Role-grid (note #4) satisfied: operator output→source/sink, solve trace→flux.
- **LAYER CONTRACT (not a wart)**: `MultiplicationOperator` (L2) has NO `__add__` → falls to
  `LinearOperatorMixin.__add__`→`OperatorSum` (cannot reach L3 `InvertibleOperator`). Only
  `CollisionOperator.__add__` (L3) dispatches `L+C→InvertibleOperator`. numerics↛transport↛sn
  layering respected; WDD-sweep dispatch stays L3. MRO-confirmed.

## 2 NITS (both do-now, neither blocking architecture)
- ⭐ **NIT-1 DEAD IMPORTS** (P5): `CAP_APPLY`/`CAP_APPLY_TRANSPOSE`/`CAP_SOLVE` imported at :101-103
  but appear in body ONLY in docstrings/comments (:140/:144/:146/:163/:218) — NEVER a runtime symbol
  (`capabilities` inherited at :171). pyright config doesn't flag unused-import → CLI silent but dead.
  Bug habitat = the IMPORT LIST asserts a coupling the code lacks (inverse of the S3a docstring-asserts-
  coupling tell) → a DRY-sweep reader infers this op assembles its own caps = the re-derivation P4
  avoided. Fix: drop the 3 from the import block, keep BlockRole/DiagonalOperator/LinearOperatorMixin
  (docstring `:data:` refs stay valid as prose).
- ⭐ **NIT-2 `_engine()` REBUILD-PER-CALL** (P2/P7): the SAME
  `DiagonalOperator(coefficient.values, broadcast_axes=(0,))` built in 3 places — `__post_init__:168`
  (discarded after reading `.capabilities`) + `_engine():183` on EVERY apply + EVERY solve.
  NOT free: `DiagonalOperator.__init__` runs `np.asarray(dtype=float)` + an O(size) `np.all(coeff!=0)`
  spectrum scan PER matvec — recomputing a capability already frozen in `self.capabilities`. RECOMMENDATION
  = STORE-ONCE: `engine: DiagonalOperator = field(init=False, repr=False)` set in `__post_init__`,
  apply/solve read `self.engine`, DELETE `_engine()`. Bug habitat = 3 spellings of `broadcast_axes=(0,)`;
  a future re-home (metric-weighted projection / N-D N>1) lands in `_engine()` not `__post_init__` →
  caps gate and action DIVERGE. Frozen-lens-over-cache shape (engine derived-once from already-cached
  `.values`, stateless → can't go stale) — distinguish from S2's "don't cache the FIELD" (that was a
  2nd mutable slot; this is a derived-once stateless lens). STANDING RULING: build-once-and-store BEATS
  build-as-method whenever the built object is stateless-derived-from-an-immutable-source AND a capability
  is computed from it (the cap freeze + the action must read the SAME object).

## Minor observations (NOT findings)
- `.sigma`→`coefficient.values` single-sourced; mesh-identity invariant holds (`__init__` stores SAME
  `sn_mesh` before `super().__init__`); ndarray↔CrossSectionField via the SAME `from_mesh` factory →
  structurally identical coefficients on both paths.
- Intra-expression asymmetry at the 3 solver.py sites: `StreamingOperator` still takes raw
  `total_cross_section` while `CollisionOperator` takes the field — HONEST (streaming not yet
  field-promoted), not drift. Retires when a future S-stage promotes StreamingOperator (track there).

STANDING TELLS reinforced: (1) a dead IMPORT (not just a docstring) can assert a false coupling —
grep the body for each imported name as a runtime symbol; (2) build-once-store > build-per-call when
the built object freezes a capability AND is stateless-derived-from-immutable.
