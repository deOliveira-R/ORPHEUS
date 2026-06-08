---
name: wave-o-operator-algebra
description: Wave O #208 SN operator-algebra durable facts — the role grid, the G-adjoint capability-lattice (why C/S/F=None is illegal-states-unrepresentable), and the phased-carve twin-delivery pattern. Merges the 5 O.4/B.5.2/R5 review notes (all LANDED).
metadata:
  type: project
---

Durable facts from the Wave O (#208) SN operator-algebra carve series — all LANDED on
`refactor/field-role-typing` (commits `89b2f62`..`7ccc14a` for O.2b/R5; the −B and B.5.2
carves earlier). The transient "O.2-deferred-gap" diaries are dropped (git is authoritative
that they landed); these are the reusable verdicts.

## The role grid (load-bearing review axis for any typed-field retype)
```
            .apply (= Aψ)        .solve (trace)   from_balance (defect)
  bulk      AngularSourceSink    AngularFlux      AngularResidual
  boundary  BoundarySourceSink   BoundaryFlux     BoundaryResidual
```
Discriminator for any `*.zeros_on` / `zeros_for_mesh` flip: **operator OUTPUT vs solve-TRACE**.
Operator outputs (matvec, source builds, codomain zeros) flip to source/sink; solve traces,
converged finals, persistent partner-flux, and cold-start iterates STAY flux. A residual exists
ONLY via `from_balance`. (B.5.2 reviewed every `*.zeros_on` site against this and all held; the
13-flip / 7-stay split was correct.) This REVERSED an earlier plan to type matvec output as a
residual — guard against that two-hat cross-class throw recurring. Promoted to AGENT.md #4.

## The G-adjoint capability lattice (an illegal-states-unrepresentable WIN worth remembering)
R5 gave SN composites the metric-correct Hilbert adjoint `A† = G⁻¹AᵀG` via `FullFieldSpace`
(direct sum `V_bulk ⊕ V_trace`, dispatching `apply_metric`/`inner_product` per block — pure
Pattern-2 composition, no new metric arithmetic). L and B report `full_field_space`; **C/S/F
left at mixin-default `None`, and this is SOUND by a sharper reason than domain-plumbing:**
S (Scattering) and F (Fission) advertise only `{CAP_APPLY}` — no `CAP_APPLY_TRANSPOSE`. So
`(L+C−S−F−B)` as an `OperatorSum` DROPS `CAP_APPLY_TRANSPOSE` (sum needs BOTH operands have it),
and `_AdjointOperator.apply` raises `MissingCapability` BEFORE reading any domain. The feared
"`(C−S).H` silently goes Euclidean" CANNOT happen — it raises. The only reachable adjoint
composites are transpose-closed (`(L+C)`, `(L+C−B)`, `(L+B)`), every one contains L which carries
the composite metric. **Illegal-states-unrepresentable holds via the CAPABILITY LATTICE, not the
domain plumbing.** A shared `SNOperatorMixin` threading the mesh into C/S/F would be premature
abstraction AND would muddy the bulk-only types. This is the canonical example for: when a
feared-unsafe path is actually UNREACHABLE, find the capability/type fact that makes it
unrepresentable, don't add defensive plumbing.

## Phased-carve twin-DELIVERY pattern (the dominant recurring smell — promoted to AGENT.md #1)
The −B extraction (Krylov-commit / SI-commit / O.4b-2D) single-sourced the `SNBoundaryOperator`
but left parallel *delivery* routes (driver `S+B` fold vs `_reflect_outflow_into_inflow` helper;
`OperatorSum.apply` fold vs Krylov-inline fold). NOT a math twin (one operator) — a plumbing
twin. Verdict CONCERN-not-VIOLATION because both routes consume the one B; hazard is a future
transform landing in one route only. O.4b was the FIRST in the series to collapse a twin at the
MATH level (routed the 2-D matvec through the same `DiamondDifference`/`SweepDependencyGraph`
wavefront the sweep uses, retiring the bespoke `_apply_2d_cartesian` FD stencil that missed
k_inf ~10%). The standing collapse destination is "honest composition: drivers take whole
`L+C−S−F−B`."

## Specific recurring tells caught in this series (all promoted to AGENT.md #7)
- False ordering-rationale comment ("S summed FIRST so the domain check skips" — the check skips
  whenever EITHER domain is None, symmetrically; ordering irrelevant).
- Keystone deletion leaves `bc_outer`/`bc_inner` dangling, surviving only as a curvature proxy
  (`bc_left if cartesian else None`) — test curvature directly.
- Mesh-identity invariant documented-not-enforced in `SNBoundaryOperator._apply_faces` (reads
  `psi.bulk.mesh` for buffers but `self.sn_mesh.trace` for indices) → now ENFORCED with a raise
  (good; the fix pattern is "assert the identity the docstring claims").
- `reduced is not None` vs `is_1d (== ny==1)` two-spelling partition: dispatcher+guards sharing
  ONE predicate is a STRENGTH; the smell is a second spelling on one side.

Related: [[field-role-and-block-role-attrs]] (block_role + the role-grid attrs),
[[issue-208-affine-flux-algebra]] (the affine `from_balance` consumer that completes the grid).
