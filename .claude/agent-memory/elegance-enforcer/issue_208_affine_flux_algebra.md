---
name: issue-208-affine-flux-algebra
description: #208 affine flux algebra (FluxRole torsor + Displacement + typed equation-residual) durable facts — the residual/displacement asymmetry principle, the affine/torsor design, the sha256 bit-identity golden tooling, and the loss-matvec fold twin. All LANDED.
metadata:
  type: project
---

Durable facts from the #208 affine-typed SN field algebra (Pieces 1–3, all LANDED:
`8c2f355`→`9ea15d3`→`dfdd564`; skill promotion `9504146` closed #208/#201). The blow-by-blow
piece diaries are dropped; these are the reusable verdicts.

## The design (the affine/torsor frame in code)
Flux leaves get `FluxRole` mixed in BEFORE the storage base
(`AngularFlux(FluxRole, AngularField)`) → MRO precedence over `Field.__add__/__sub__`.
`__sub__` → `_mint_displacement` (dataclass field-copy swapping `values`, lands in the affine
difference space V); `__add__` is the torsor action `flux + displacement → flux` AND the gate
`flux + flux → TypeError` (illegal-states-unrepresentable: a state is not addable to a state;
scalar-mul KEPT, `__sub__` KEPT — the SI iterate-delta + source-add are legitimate).
`affine_combination(Σλ=1)`. Displacements get a `Displacement` mixin (marker + diagnostics) and
INHERIT vector-space `+`/·/`−` from Field (V is a genuine vector space). The state/displacement/
residual triple is the torsor frame promoted into the `coding-elegance` skill (commit `9504146`,
"state-typed-increment anti-pattern").

## Durable principle — the residual/displacement ASYMMETRY is principled (promoted AGENT.md #6)
Residuals have NO shared mixin (`_residual.py` does not exist) — pure thin leaves; their
distinguishing behaviour is their CONSTRUCTION (`from_balance`, a class-transition factory →
lives on the engine `Field._from_balance`). Displacements DO get a mixin — their distinguishing
behaviour is their METHODS (`contraction_ratio`/`true_error_estimate`/`where_largest`). Different
axis of variation → different mechanism → NO twin path. Likewise `_mint_displacement` correctly
does NOT reuse `_from_balance`: the latter rebuilds via `cls.from_mesh` (no `L` arg, wrong space
`"angular_residual"`); the displacement field-copy shares the flux's space + carries `L`. When you
see "X has a mixin, sibling Y doesn't," check construction-engine vs methods-mixin before flagging.

## Box-7 — `evaluate_residual` is the from_balance consumer (1-consumer-but-earned)
`evaluate_residual` (solver.py) reads like `r = (L+C−S−B)ψ − q`: composes the loss op at the call
site, routes the defect through `AngularResidual.from_balance`/`BoundaryResidual.from_balance`.
The residual math `lhs.values − rhs.values` is SINGLE-SOURCED in `Field._from_balance`; both leaves
delegate. It has NO production caller (test-only; `as_dsa_source` / DSA #2 deferred) — but it is
NOT speculation: it is the named-composition that exposes the FULL equation residual as a typed
object, distinct from the per-call matvec defect B.5.2 already retyped. The from_balance mint
existed unconsumed (role-grid completion); this is the production-reachable site that completes the
type. Earned per "1 consumer fine when benefit established" (DSA is the genuine future consumer).

## The ONE Pattern-2 nit — the loss-matvec FOLD appears twice (acceptable-for-now)
`evaluate_residual` folds `(L+C−S−B)·ψ` via `OperatorSum.apply` (clean, off-path); the Krylov
driver (`iteration.py:750`) inlines `out = L.apply(psi); for g: out -= g.apply(psi)` (justified
"no intermediate OperatorSum allocation"). SAME algebra, TWO fold shapes, single-sourced at the
LEAF (`L/S/B.apply` each defined once). Same shape as the O.4 twin-matvec "verified identical,
acceptable-for-now" judgment — a perf-motivated scipy-allocation specialization, not a divergence
risk while leaves are the single source. Flag if a THIRD fold appears. (This fold-twin is the same
family as [[wave-o-operator-algebra]] phased-delivery and is promoted to AGENT.md #2.)

## Tooling worth reusing — the sha256 bit-identity GOLDEN
`test_affine_carve_bit_identity.py` freezes sha256 of raw float64 bytes at the pre-carve commit,
covering ALL THREE displacement leaf types (MomentDisplacement via windowed 2-D SI;
AngularDisplacement via 1-D slab SI; the Krylov flat-residual path). Sharper than a tolerance-gated
DD snapshot suite (those pre-drift ~6920 ULP) — this is the RIGHT tool for a zero-numerical-change
type-only carve. Recommend it whenever a carve claims bit-identity.

## Latent seam (resolved) worth remembering as a pattern
Piece 1's `psi − psi_prev` in SI already minted a Displacement; bit-identity held only because the
displacement reached `_l2_norm` as a `TimedFullField` composite (HAS `to_flat` → reads `.values`),
NEVER the bare-leaf `np.asarray` fallback (which yields a 0-d OBJECT array that breaks the norm).
The duck-type guard `_flux_displacement_leaf` (mirrors `_is_ravellable`, avoids numerics↛transport
L1↛L2 import) reaches `.bulk` for a composite, the field for a direct leaf, None for bare ndarray.
Lesson: a typed wrapper that adds a class to an existing arithmetic path can silently change the
`np.asarray`/duck-type fallthrough — verify the new type still hits the same accessor (`to_flat`),
not the object-array fallback.

Related: [[wave-o-operator-algebra]] (the role grid + the loss-op the residual evaluates),
[[field-role-and-block-role-attrs]] (the per-leaf mixin/attr mechanism).
