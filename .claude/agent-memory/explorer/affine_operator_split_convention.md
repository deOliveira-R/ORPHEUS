---
name: affine-operator-split-convention
description: ORPHEUS has NO affine-operator type and this is a ruling, not a gap — an affine map in a linear operator slot broke GMRES's Arnoldi relation (measured |B(0)| = q). Every affine law is SPLIT into a linear operator plus a source on the rhs. Governs any future "should I mint an affine operator?" question.
metadata:
  type: reference
---

Durable architectural convention, verified 2026-08-13 on
`refactor/operator-strategy-layers` @ `bbeb5000`. Line numbers drift; the
RULING does not.

## There is no `AffineOperator`, and that is deliberate

A tree-wide grep for `class Affine` / `AffineOperator` / `AffineMap` /
`affine_part` / `linear_part` over `orpheus/` returns **zero** operator types.
The only `linear_part` is `RigidMotion.linear_part` in
`orpheus/geometry/transformation.py` — an E(d) group element, unrelated.
"Affine" in `orpheus/transport/` names the **field-role torsor** algebra
(`flux + flux → TypeError`), not an operator concept.

## The pattern every affine law follows: SPLIT, never one object

`γ₋ψ = L γ₊ψ + q` is realized as **two tiers**:

- the **linear factor `L = R·G`** is an ordinary `LinearOperator`, stamped
  `BlockRole.BOUNDARY`. For prescribed-inflow and vacuum alike, `L` is the
  ZERO morphism (`ZeroOperator`, both spaces bound) — a zero morphism is an
  ordinary linear operator.
- the **inhomogeneous `q`** travels a separate typed SOURCE channel
  (`AngularBoundarySourceSink` → `q.boundary`) and is added by the solve.

Doc of record: `docs/theory/foundations/boundary_conditions.rst`
`.. _bc-affine-source-channel:` (the `list-table` "Every SN-realizable law as
a linear factor plus a source" enumerates all six laws' `L` and `q`).

## ⭐ WHY — the measured failure that produced the ruling

Before campaign phase P3 (2026-08-05) one law realized BOTH halves into a
single object declared linear. Recorded verbatim in `BoundaryOperator`'s
docstring in `orpheus/numerics/operator.py`:

> measured `|B(0)| = q` and `|B(2x) − 2B(x)| = q`, and on the Krylov path a
> raised `ConvergenceCertificateError`, because an affine map breaks the
> Arnoldi relation GMRES's residual depends on.

Two further consequences worth carrying:

- **`q` was delivered TWICE** — once by the affine operator inside `B`, once
  by the source channel. The affine object was double-counting.
- **The `block_role` stamp is NOT the fence.**
  `SNBoundaryOperator._face_laws` collects every face's law with no
  `block_role` filter, so an unstamped affine leaf reached the block anyway.
  Do not reach for "leave it unstamped" as a containment strategy; it was
  tried and it does not contain anything.

## How to apply

When a campaign proposes realizing some march / recurrence / update as a
"first-class **affine** operator with a linear part and a residual", the
answer is NOT a new affine base class. It is:

1. name and realize the **linear part** as a `LinearOperator`
   (`orpheus/numerics/operator.py`, the `Protocol[Domain, Codomain]`);
2. route the **inhomogeneous part** as a typed SOURCE on the rhs;
3. if the inhomogeneous part is a DEFECT rather than a datum, the tree's
   existing vocabulary for it is the `*Residual` role leaf
   (e.g. `AngularBoundaryResidual` = "the defect of the affine boundary law").

Anything else re-opens a bug the tree already paid for. The worked precedent
to copy is the four-block ψ½ split in
`orpheus/sn/operators/radial_characteristic.py`: `A_BB` forward + resolvent,
`A_AB` the coupling — each a plain `LinearOperator`, the data carried by typed
composite fields, no affine object anywhere.
