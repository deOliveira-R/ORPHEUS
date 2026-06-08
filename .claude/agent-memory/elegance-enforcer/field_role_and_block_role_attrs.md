---
name: field-role-and-block-role-attrs
description: ORPHEUS typed-field convention — role-determined constants live as per-leaf class attrs; the dataclass-field trap dictates ClassVar (non-dataclass base) vs plain-attr (frozen-dataclass leaf). Covers field UNITS/_SPACE_NAME and operator block_role.
metadata:
  type: project
---

The durable mechanism behind two ORPHEUS typed-class conventions. (The recurring
PRINCIPLE — role-not-family + the dataclass trap — is promoted into AGENT.md §
Institutional knowledge #5; this note holds the concrete instances for verification.)

## Field-role constants (`orpheus/transport/fields|sources|residuals`)
Constants determined by a field's ROLE (not its storage family) are declared on each
concrete role leaf, NOT on the storage-family bases (`AngularField`/`ScalarField`/
`BoundaryField` in `_bases.py`). Instances: `_SPACE_NAME: ClassVar[str]` (bare annotation
on the bases, value per leaf; `from_mesh` reads `cls._SPACE_NAME`); `UNITS: ClassVar[Unit]`
(4 named dimensional-class constants in `orpheus/numerics/units.py`; the 10 leaves map onto
them, units coarser than classes).

**Why per-leaf is CORRECT, not duplication:** roles cut ACROSS families — `BoundarySource`
(a source) carries flux units because its trace is all-flux. A constant on `AngularField`/
`ScalarField` would force the wrong value on a sibling leaf. Per-leaf is enumeration of a
genuine many-to-fewer map, SSOT held in the named constants the leaves reference. Inline
construction of the value across leaves (vs referencing the shared constant) IS a duplication
VIOLATION — that is the line to police.

The arithmetic gate is class identity (`Field._check_partner`), NOT units — `UNITS` is
diagnostic metadata only. Do NOT flag "same units on different classes" as illegal-states;
that is intended (`test_field_units.py::TestSameUnitsDifferentRole`).

## Operator block_role (`orpheus/numerics/operator.py`)
`BlockRole` enum {BULK, FULL, BOUNDARY} = SSOT for the 2×2 bulk/boundary block structure.
`block_role: ClassVar[Optional[BlockRole]] = None` on `LinearOperatorMixin` (not a dataclass →
ClassVar there is documentary). **Frozen-dataclass leaves tag with a PLAIN unannotated attr**
`block_role = BlockRole.X` — deliberately NOT `ClassVar[...]`, because under
`from __future__ import annotations` a stringized ClassVar slips past dataclass field-detection
and becomes a FIELD. Same robustness reason as the field-role pattern, OPPOSITE mechanism
(here the leaf is a frozen dataclass so the annotation is DROPPED, not added). Do NOT "fix" the
missing ClassVar.

`_BlockRoleMeta.__instancecheck__` makes `isinstance(L, FullOperator)` read
`op.block_role is cls._role`. Value-based (not `@runtime_checkable Protocol`) is CORRECT: every
operator carries the attr so a structural Protocol would match all three roles; classification
must be instance-level (generic Identity/Scaled/Permutation play a boundary role).

**LANDED (was a latent O.2 gap, now history — git authoritative):** `_join_block_roles`
(operator.py:158/725) derives composite role from operands; `_AdjointOperator`/`ScaledOperator`
forward via `getattr(inner,'block_role',None)` (612/856); `InvertibleOperator`'s hand-set FULL
stamp retired in favour of derivation (`7ccc14a`). The per-site `_as_boundary` stamp in
`SNBoundaryRealizer.realize` correctly STAYS as the leaf-tagging mechanism (its leaves are
genuine TP/Scaled primitives); only composition gains derived role. So when reviewing a NEW
operator: a leaf hand-stamps, a composite must DERIVE — never hand-stamp a composite as a
parallel special case.

Related: [[issue-208-affine-flux-algebra]] (the role grid these attrs live under).
