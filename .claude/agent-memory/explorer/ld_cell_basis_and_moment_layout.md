---
name: ld-cell-basis-and-moment-layout
description: The LD tensor-Legendre cell basis (P_k(1)=1, Gram diag(1,1/3) per axis, axis-0-outer kron) vs numerics' angular LegendreBasis; moment_layout member classes; Declared moment-axis hazards (2026-10-25, HEAD 2258223d)
metadata:
  type: project
---

Snapshot from 2026-10-25 (ruling 19 of the layer-DAG campaign). Detail is in `scratch/layer_dag/moment_basis/context.md`.

- The scheme's cell basis per axis is `{1, ξ}` on `ξ ∈ [-1, 1]`, normalised at the endpoint (`P_k(1) = 1`). Its per-volume Gram is `diag(1, θ)`, and θ = 1/3 is EXACTLY the Legendre Gram (a basis fact written as `LinearDiscontinuous.theta`). The tensor is a kron with axis 0 outermost. The face trace sums the `o_a` blocks (`B(+1) = [1, 1]`). The octant sign map is the `(-1)^k` parity.
- `numerics.basis.LegendreBasis` is ANGULAR. It lives on `S^2/O(2)_a`, its Gram is `4π/(2ℓ+1)`, and its `LegendreLayout` is a DegreeLayout, so `moment_head_axis` would accept it as an angular head. Do not reuse it for the cell basis. Only `legendre_table` is reusable.
- `moment_layout` members fall into five classes: basis (`AVERAGE_MOMENT`, `cell_moment_count`), face trace (`face_moment_count`), axis label, rank predicates, and the "append iff > 1" tail policy, which is space composition and not a basis fact.
- No gate pins the LD axis as `Declared`: `spaces.rst` cites G5_6a, which does not test it. `FromBasis.kind` is a closed match over the layouts.
- Issue #423 is this work, written in stale vocabulary (`generator`, `mass_matrix`).

Related: [[moment-metric-slots-and-cross-axis-forms]], [[s7-axis-identity-design-census]].
