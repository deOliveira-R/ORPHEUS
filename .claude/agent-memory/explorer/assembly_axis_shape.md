---
name: assembly-axis-shape
description: Assembly axis at 2375992d (pre-S3) — 12 assemble owners, 4 composite laws not 3, one production route consumer (MatrixInverseOperator), spy recipe
metadata:
  type: project
---

Shape of the assembly axis at `2375992d` (S3 census; tables in
`scratch/numerics_untangling/packages/S3_census.md` while it lives).

- Owners of `assemble`: the composites Sum/Product/FusedAtom (operator.py) AND `CoupledOperator`
  (`block_array`, a 4th n-ary law the plan row omits); unary wrappers `FlattenedOperator` (returns the
  inner's carrier by identity, `is`-pinned) and `BulkLift` (embed); leaves Leakage, DiffusionBoundary,
  IsotropicTransfer, MultiplicationOperator; the carrier itself (`return self`, `is`-pinned).
- Leaves are not pure matrix builders: DiffusionBoundary probes each law via `as_matrix`;
  IsotropicTransfer extracts via ng `apply` calls.
- The `as_matrix` → assemble route has ONE production consumer: `MatrixInverseOperator.__init__`
  (diffusion resolvent, homogeneous solver, `CoupledOperator.inverse`). SN produces 0 routed calls.
- Concept-order gate counts TYPE_CHECKING and lazy imports; same-tier edges are violations too.

**Why:** S3 (ruling 24) moves the laws downstream and makes leaves return raw matrices.
**How to apply:** a `-p` plugin wrapping `LinearOperator.as_matrix` / `_as_matrix_by_probing` / each
class's own `assemble` with caller frames is the census; the 5 gate suites serial `-m "not slow"` take
~34 min (xdist refuses: random-array param ids). Re-verify lines; operator.py moves every wave.
Related: [[operator-composite-algebra-seam]].
