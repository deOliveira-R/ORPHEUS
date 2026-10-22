---
name: s7-axis-identity-design-census
description: S7 design census at dd17f678 — axis mint homes, 4 generated-eq space subclasses, alias-hidden Axis mints, 4 MomentHead doors
metadata:
  type: project
---

Facts at dd17f678 (re-verify before repeating):
- Every angular NODAL axis is minted through `Axis.of_quadrature` (5 prod callers). Nodal mints with no atoms in scope: EnergyAxis from_grid/synthetic, material_mesh rank-d spatial, homogeneous `Axis("spatial",(1,))`.
- `material_mesh.py` imports `Axis as SpaceFactorAxis`: a name-keyed AST census misses 2 mints without alias resolution.
- Four FunctionSpace subclasses (AngularTraceSpace, AngularFaceTraceSpace, RadialCharacteristicInterior/BoundarySpace) are re-decorated `@dataclass(frozen=True)`, so their eq/hash are GENERATED and override the base's delegation. They compare (name, shape) with an exact-class check and ignore axes. Detect this at runtime: `__eq__.__code__.co_filename == "<string>"`.
- There are 4 `isinstance(.., MomentHead)` doors; the 4th is `truncated_head` itself.
- In a worktree dispatch, the write-scope hook admits memory only under the MAIN checkout's `.claude/agent-memory/explorer/`, and refuses a subagent Write of report.md. Put the tables in the handback.

**Why:** S7 (numerics untangling, #631) moves the MomentHead surface onto the head axis and makes nodes and support axis content.
**How to apply:** start any S7 follow-up census from these; [[s6bc-census-facts]] is adjacent.
