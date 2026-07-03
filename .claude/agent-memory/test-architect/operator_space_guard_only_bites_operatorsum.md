---
name: operator-space-guard-only-bites-operatorsum
description: The OperatorSum/Product domain/codomain space guard is INVISIBLE to the SI/Krylov matvec (realised as L.apply−Σgᵢ.apply, no OperatorSum); it bites only on actually-composed OperatorSums. Verify space activation on those, not on the converging path.
metadata:
  type: project
---

When verifying an operator-algebra carve that activates the
`OperatorSum`/`OperatorProduct` runtime composability guard
(`operator.py` `IncompatibleOperatorComposition` on domain/codomain
mismatch), do NOT assume the within-group convergence path is
protected by it.

**Fact (ORPHEUS SN, re-validated 2026-06-25 at HEAD `0a02cf7`, branch
`refactor/operator-inverse-algebra`, post W-A/W-B/W-C):**
`KrylovAcceleration` / `SourceIteration` realise the within-group matvec
as `L.apply(ψ) − Σ gᵢ.apply(ψ)` per call — *"no intermediate
OperatorSum allocation"* (`iteration.py:600-601`; gains looped at
`:526-527` via the carrier `rhs = rhs + g.apply(psi)`, NOT operator
algebra). The drivers validate ONLY `capabilities` at construction,
never `domain`/`codomain`. So the space guard NEVER fires on
SI/Krylov/eigenvalue. **F NEVER enters any OperatorSum in production** —
it is always `F.apply(ψ)` (production estimator `iteration.py:290/319`;
external source `q_ext = F.apply(ψ)/k`). The guard's only
production-reachable `OperatorSum` is `InvertibleOperator(L, C)` (the
`L + C` build inside `_within_group_triple`, `solver.py:225`). The
SECOND site from the old memo — `evaluate_residual`'s `loss_op =
(LC − S − B)` — is **NOT production-reachable at HEAD**:
`evaluate_residual` (`solver.py:231`) has ZERO production callers (the
"#208 box-7 consumer", currently consumed ONLY by
`test_typed_residual_evaluation.py:74` + MMS tests). The composition
`(LC − S − B)` it documents is built in the TEST, not production.

**Post-W-C carrier state (drifted from the old memo):** at HEAD,
**L, LC, B already report `FullFieldSpace('sn_full_field', shape=(…,))`**
(the IDENTICAL `SNMesh.full_field_space` cached object); **C and S
report `None`** (the W-D targets along with F). B carrying a real space
means the `B`-arm of `(LC − S − B)` is ALREADY guarded today (B vs LC,
both correct); only the `S`-arm is skipped (S=None). `_within_group_triple`
returns `(L+C, S, B)` as a TUPLE at `solver.py:224-228`.

**Why it matters (the inverted assumption):** the naive read is "the
matvec composes the operators, so giving S/F real spaces protects the
matvec." FALSE — the matvec calls `.apply()` separately and subtracts
in the carrier (`TimedFullField`/`SourceSink`) `__sub__` algebra; the
guard is bypassed. Activating S/F spaces changes behaviour ONLY at the
two OperatorSum sites.

**How to apply:** when a carve gives previously-`None`-space operators
real spaces, target the §3.1-style activation gate at the
ACTUALLY-COMPOSED `OperatorSum` (`LC - S - B`), not the converging
path. The blocker condition to pin: the new space's `(name, shape)`
must equal the sibling composite's — `FunctionSpace.__eq__` is by
`(name, shape)`, NOT object identity (`space.py:157-163`), and
`SNMesh.full_field_space` is a `@cached_property` so same-mesh
operators share the identical object. All of L+C, S, F, B converge to
`('sn_full_field', (N·ng·…,))` — a probe confirmed `LC - S - B`
composes both today (S=None, guard skipped) AND when S carries
`full_field_space` (guard active). The blocker only manifests if the
carve mints a NEW space name for S/F. See [[lessons]] L4 (prove every
gate's teeth) — mutate `type(S).domain` in-process (monkeypatch, never
git checkout) to a wrong-named space and confirm the OperatorSum reds.
