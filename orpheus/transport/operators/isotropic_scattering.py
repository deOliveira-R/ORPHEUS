r"""The model-independent isotropic energy operators on the scalar flux.

Every transport model — discrete ordinates (SN), collision probability (CP),
method of characteristics (MoC), diffusion, the homogeneous 0-D solver, Monte
Carlo — builds the **isotropic** (:math:`\ell=0`) in-scatter source by the SAME
per-cell energy operation on the **scalar flux** :math:`\phi`:

.. math::

    Q^{\rm iso}_g(\vec r) \;=\; \sum_{g'} \Sigma_{s,0}(g'\to g)\,\phi_{g'}(\vec r)
        \;+\; 2\,\sum_{g'} \Sigma_{2n}(g'\to g)\,\phi_{g'}(\vec r) .

What differs per model is ONLY the *angular wrapper* around it — how :math:`\phi`
is obtained from the model's flux (SN: :math:`\phi=\sum_n w_n\psi_n`; diffusion /
homogeneous: :math:`\phi` is the native unknown; CP / MoC: per region / per FSR)
and how the scalar source is embedded back (SN: :math:`/W` isotropic broadcast;
diffusion / homogeneous: fold :math:`-\Sigma_{s,0}^{T}` into the loss matrix
:math:`A`; CP / MoC: through the transport kernel). The **energy operation is
model-independent** (cross-domain-attacker + explorer, 2026-06-28: the same
``einsum("fg,f->g")`` is re-implemented 6× across the solver families — a
Cardinal-Rule-2 single-concept-many-places situation).

This module owns that shared energy operation as two
:class:`~orpheus.numerics.operator.LinearOperator`\ s on the scalar flux:

* :class:`IsotropicScattering` — :math:`\Sigma_{s,0}` (P0 in-scatter);
* :class:`IsotropicN2N` — :math:`2\,\Sigma_{2n}` (the (n,2n) doubling, a DISTINCT
  *multiplication* channel — physics-faithful, and it also feeds the keff
  *production* numerator, not the loss — so it stays its own operator, summed with
  :class:`IsotropicScattering` for the in-scatter source).

Both are the **scalar (:math:`\ell=0`) realization** of the moment-space operators
:class:`~orpheus.transport.operators.scattering.LegendreMomentScattering` (at
:math:`\ell=0`) and
:class:`~orpheus.transport.operators.scattering.N2NMomentOperator`: they route
through the SAME per-material kernel-field verbs
(:meth:`~orpheus.transport.material_field.ScatteringMaterialField.add_p0_source` /
:meth:`~orpheus.transport.material_field.N2NMaterialField.add_emission` + the
``…_transpose`` twins — CS4c step 3, the O-6 landing), so the cross-section
DATA and the per-material dispatch live ONCE (``coding-elegance`` Pattern 2). The harmonic-moment
:attr:`~orpheus.transport.operators.scattering.ScatteringOperator.full_scatter_kernel`
is the permanent verification oracle for this scalar form.

Layout & carriers
=================

Bare ``np.ndarray`` in / out — the model-portable contract (CP / MoC / diffusion
feed raw scalar-flux arrays; SN passes ``phi.values`` and re-wraps). The action is
**spatial-moment-axis-agnostic** (the trailing ``…`` rides an LD ``2^d`` spectator
axis through, #240 D5b-S3), so a φ̂-carrying scalar flux scatters every spatial
moment by the same :math:`\Sigma_s`.

Capabilities
============

apply + a working ``apply_transpose``. No inverse — the per-cell group-transfer
matrix is generally singular (a thermal group with no up-scatter source has a zero
row); it is *applied*, never inverted. The transpose IS advertised (campaign #276):
:math:`\Sigma_{s,0}^{T}` / :math:`2\Sigma_{2n}^{T}` (the group-axis flip) is the
group-asymmetric factor of the adjoint isotropic scattering source.

References
----------

* explorer (2026-06-28): cross-model grounding — the 6× / 5× duplication inventory.
* cross-domain-attacker (2026-06-28):
  ``iso_source_frame_conjugation_unification.md`` — relocate the shared energy
  operator; do NOT mint a ``ConstantBasis`` iso-frame (it forks ``R∘M``).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

from orpheus.numerics.operator import (
    BlockRole,
    LinearOperator,
)
from orpheus.transport.material_field import (
    N2NMaterialField,
    ScatteringMaterialField,
)
from orpheus.transport.operators.bound_operator import BoundOperator
# Runtime import for the composite-arm isinstance parse (mirrors
# fission.py / multiplication_operator.py): ``FullField`` is a leaf in
# the transport dependency graph (it imports no operators), so this
# module-level import is cycle-free.
from orpheus.transport.full_field import FullField

if TYPE_CHECKING:
    from orpheus.numerics.assembled_operator import SparseAssembledOperator
    from orpheus.numerics.space import FunctionSpace
    from orpheus.transport.mesh.material_xs_field import MaterialXSField


__all__ = ["IsotropicScattering", "IsotropicN2N"]


def _values_of(phi: "np.ndarray | object") -> np.ndarray:
    """Read the bare ``(ng, *spatial[, 2^d])`` array off a flux carrier or ndarray."""
    return np.asarray(getattr(phi, "values", phi))


def _scalar_composite_source(op: "LinearOperator", psi: FullField) -> FullField:
    r"""The shared scalar-composite arm of the two iso energy operators (#290 P4).

    Parses the composite's bulk as the SCALAR family (the iso operators
    are scalar-flux operators — an angular composite routes through
    :class:`~orpheus.transport.operators.scattering.ScatteringOperator`,
    never here), runs the operator's own bare-ndarray kernel on the bulk
    values (single source of truth — the SAME per-cell energy matmul),
    and wraps the result as the closed scalar source composite:
    :class:`~orpheus.transport.source_sinks.scalar_source_sink.ScalarSourceSink`
    bulk ⊕ implicit-zero
    :class:`~orpheus.transport.source_sinks.scalar_boundary_source_sink.ScalarBoundarySourceSink`
    boundary (an energy-transfer operator has no face-trace action —
    the in-scatter term is a CELL quantity). Shared by
    :class:`IsotropicScattering` and :class:`IsotropicN2N` (Pattern 2:
    one arm body, two kernels through ``op.apply``).
    """
    from orpheus.transport.fields._bases import ScalarField
    from orpheus.transport.source_sinks import (
        ScalarBoundarySourceSink,
        ScalarSourceSink,
    )

    bulk = psi.interior
    if not isinstance(bulk, ScalarField):
        raise TypeError(
            f"{type(op).__name__} composite apply: scalar-family bulk "
            f"required (the iso energy operators act on the scalar flux; "
            f"an angular composite's in-scatter routes through "
            f"ScatteringOperator); got {type(bulk).__name__}."
        )
    # CS4b S4 — the space route: output blocks ride the operand's blocks.
    return FullField(
        interior=ScalarSourceSink(values=op.apply(bulk.values), space=bulk.space),
        boundary=ScalarBoundarySourceSink.zeros(psi.boundary.space),
    )


def _iso_is_assemblable(op: "IsotropicScattering | IsotropicN2N") -> bool:
    r"""Shared assembly-axis predicate of the two iso energy operators.

    Emittable iff the composite flat layout is known AND the declared
    bulk is the SCALAR family's ``(ng, *cells)`` energy-first tensor —
    the same contract the bare-ndarray kernel arm assumes (a
    block-bearing :class:`FullFieldSpace` whose bulk leading axis is
    the group axis). A scalar-space binding (the diffusion / SN
    internal consumers — ends mandatory since CS4c step 3) honestly
    refuses: no composite flat layout to emit into.
    """
    from orpheus.numerics.spaces.full_field_space import FullFieldSpace

    space = op.domain
    return (
        isinstance(space, FullFieldSpace)
        and space.interior_space is not None
        and len(space.interior_space.shape) >= 1
        and int(space.interior_space.shape[0]) == int(op.data_ng)
    )


def _assemble_iso_energy_operator(
    op: "IsotropicScattering | IsotropicN2N",
) -> "SparseAssembledOperator":
    r"""Shared assembly arm of the two iso energy operators (Pattern 2 —
    one body, two kernels through ``op.apply``, the
    :func:`_scalar_composite_source` sibling).

    **Group-impulse extraction — one source with the production
    kernel.** The per-cell energy blocks are read THROUGH the
    operator's own bare-ndarray ``apply`` (the einsum kernels
    ``apply_p0_in_scatter`` / ``apply_n2n``): impulse ``g'`` is the
    field ``e_{g'} ⊗ 1_cells``, whose image column IS every cell's
    block column ``M_cell[:, g']`` exactly (unit inputs make the
    kernel's products coefficient reads — no summation error). ng
    kernel calls emit the whole cell-block-diagonal bulk. Deliberately
    NOT ``dense_per_material()`` — that is the storage-side
    transpose-convention ORACLE (vv L11) and must stay
    realization-independent of production consumption.

    The emission is bulk-only (an energy-transfer operator has no
    face-trace action — the same zero-boundary fact the composite arm
    encodes), scattered into the composite flat layout ``[bulk C-ravel
    | trace]``.
    """
    from scipy import sparse

    from orpheus.numerics.assembled_operator import SparseAssembledOperator
    from orpheus.numerics.operator import MissingAssembly
    from orpheus.numerics.spaces.full_field_space import FullFieldSpace

    space = op.domain
    if not _iso_is_assemblable(op):
        raise MissingAssembly(
            f"{type(op).__name__}.assemble requires a block-bearing "
            f"FullFieldSpace with the scalar (ng, *cells) energy-first "
            f"bulk; this instance is space-anonymous or its bulk leading "
            f"axis is not the group axis."
        )
    assert isinstance(space, FullFieldSpace)  # narrowed by the predicate
    assert space.interior_space is not None
    interior_shape = tuple(space.interior_space.shape)
    ng = interior_shape[0]
    n_cells = int(np.prod(interior_shape[1:])) if len(interior_shape) > 1 else 1
    n_interior = ng * n_cells
    n_total = int(space.shape[0])

    bulk_rows = np.arange(n_interior)
    cell_of_row = bulk_rows % n_cells
    rows: list[np.ndarray] = []
    cols: list[np.ndarray] = []
    vals: list[np.ndarray] = []
    for g_from in range(ng):
        impulse = np.zeros(interior_shape)
        impulse[g_from] = 1.0
        image = np.asarray(op.apply(impulse)).ravel()   # (n_interior,) g-major
        nonzero = image != 0.0
        rows.append(bulk_rows[nonzero])
        cols.append((g_from * n_cells + cell_of_row)[nonzero])
        vals.append(image[nonzero])

    matrix = sparse.coo_array(
        (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
        shape=(n_total, n_total),
    )
    return SparseAssembledOperator(matrix, domain=space, codomain=space)


@dataclass(eq=False)
class IsotropicScattering(BoundOperator):
    r"""The P0 isotropic in-scatter energy operator :math:`\Sigma_{s,0}` on the scalar flux.

    Per cell of material :math:`m`, :math:`(\Sigma_{s,0}^{T}\phi)_g =
    \sum_{g'}\Sigma_{s,0}^m(g'\to g)\,\phi_{g'}` — the per-material group-transfer
    matmul, vectorised over cells. The model-independent half of the isotropic
    in-scatter source (the other half is :class:`IsotropicN2N`); both route through
    :class:`~orpheus.transport.mesh.material_xs_field.MaterialXSField` so the
    per-material dispatch lives once.

    Parameters
    ----------
    scattering : ScatteringMaterialField
        The scattering channel's kernel field — the per-material Legendre
        stacks over the mesh layout (only the ``p0`` head is consumed;
        the tier-2 mint truncates to order 0 so the instance retains
        exactly what it reads).
    domain, codomain : FunctionSpace
        The two mandatory ends (kw-only, write-once —
        :class:`~orpheus.transport.operators.bound_operator.BoundOperator`,
        CS4c step 3): the scalar-flux space, both — the endomorphism
        sugar lives on :meth:`from_material_xs`. The OperatorSum
        composition guard and the per-END ng-conformity admission both
        read them; the pre-CS4c ``space=None`` anonymity is retired.
    """

    scattering: "ScatteringMaterialField"
    # A BULK energy operator (the scalar flux is the bulk block); no boundary action.
    # Class-level constant (unannotated so the dataclass does not treat it as a field).
    block_role = BlockRole.BULK

    def __post_init__(self) -> None:
        # CS4a K2 → CS4c per-END form: refuse ends whose EnergyAxis
        # contradicts the data's ng (reach + declared inertness:
        # _energy_conformity docstring).
        self._assert_energy_extent_both_ends(
            self.scattering.ng, operator="IsotropicScattering",
        )

    @classmethod
    def from_material_xs(
        cls, mat_xs: "MaterialXSField", *, space: "FunctionSpace",
    ) -> "IsotropicScattering":
        r"""Tier-2 extract-and-mint: the P0 truncation of the facade's
        scattering channel, endomorphic on one ``space=``."""
        return cls(
            ScatteringMaterialField.from_material_xs(mat_xs).truncated(0),
            domain=space,
            codomain=space,
        )

    @property
    def data_ng(self) -> int:
        """The bound data's group count (the assembly helpers' read)."""
        return self.scattering.ng

    @property
    def is_adjointable(self) -> bool:
        # apply_transpose realises Σ_{s,0}χ (the group-flip A^T); caps ⊇
        # apply_transpose. is_invertible inherits base False.
        return True

    def apply(self, phi: "np.ndarray | FullField | object") -> "np.ndarray | FullField":
        r""":math:`\Sigma_{s,0}^{T}\phi` — the per-cell P0 in-scatter source.

        Bare ndarray / flux-carrier in → bare ndarray out (the
        model-portable contract). A scalar :class:`FullField` composite
        in → the closed scalar source composite out (#290 P4 — see
        :func:`_scalar_composite_source`; the kernel is the SAME bare
        arm either way).
        """
        if isinstance(phi, FullField):
            return _scalar_composite_source(self, phi)
        arr = _values_of(phi)
        out = np.zeros_like(arr)
        self.scattering.add_p0_source(out, arr)
        return out

    def apply_transpose(self, chi: "np.ndarray | object") -> np.ndarray:
        r""":math:`\Sigma_{s,0}\chi` — the group-flip transpose (the bare Euclidean :math:`A^{T}`).

        Bare-ndarray surface only: the composite transpose arm lands
        with its first consumer (the adjoint diffusion chain, #281) —
        refuse a composite loudly rather than mis-reading it.
        """
        if isinstance(chi, FullField):
            raise TypeError(
                f"{type(self).__name__}.apply_transpose: composite "
                f"FullField transpose is not yet wired (lands with the "
                f"adjoint diffusion consumer, #281); pass the bulk "
                f"values."
            )
        arr = _values_of(chi)
        out = np.zeros_like(arr)
        self.scattering.add_p0_source_transpose(out, arr)
        return out

    # ── The assembly mode (stencil-assembly 2b) ────────────────────────

    @property
    def is_assemblable(self) -> bool:
        # Shared predicate — see :func:`_iso_is_assemblable`.
        return _iso_is_assemblable(self)

    def assemble(self) -> "SparseAssembledOperator":
        r""":math:`[\Sigma_{s,0}^{T}]` — cell-block-diagonal bulk emission
        through the production kernel (:func:`_assemble_iso_energy_operator`)."""
        return _assemble_iso_energy_operator(self)

    def dense_per_material(self) -> dict[int, np.ndarray]:
        r"""Per-material operator matrix :math:`\Sigma_{s,0}^{T}` (``[g_to, g_from]``)
        — the STORAGE-SIDE oracle view, not a production consumption mode.

        Returns ``{mid: M}`` with ``M @ φ_cell == apply(φ)_cell``, i.e.
        ``M = sig_s0.T`` (``sig_s0`` is stored ``[g_from, g_to]``), read
        DIRECTLY off the stored cross sections — structurally independent
        of the :meth:`~orpheus.transport.mesh.material_xs_field.MaterialXSField.apply_p0_in_scatter`
        einsum that realizes :meth:`apply`.  That independence is this
        method's job: the verification gates use it as the
        transpose-convention oracle for ``apply``/``apply_transpose``
        (`vv-principles` L11 — the two sides of an oracle pair must not
        share a realization).  Production materialization — the LHS-fold
        ``A = C − K_iso`` the homogeneous solver densifies — goes through
        the operator's own
        :meth:`~orpheus.numerics.operator.LinearOperator.as_matrix`
        apply-to-basis instead (taxonomy §12 step 5, which retired this
        docstring's earlier claim that the fold consumers read THIS
        method — they never did).  Each entry is a fresh copy.
        """
        return {
            mid: np.ascontiguousarray(kernel.p0.T)
            for mid, kernel in self.scattering.per_material.items()
        }


@dataclass(eq=False)
class IsotropicN2N(BoundOperator):
    r"""The (n,2n) isotropic energy operator :math:`2\,\Sigma_{2n}` on the scalar flux.

    Per cell, :math:`(2\Sigma_{2n}^{T}\phi)_g = 2\sum_{g'}\Sigma_{2n}^m(g'\to g)\,
    \phi_{g'}`. A DISTINCT *multiplication* channel (each event emits two neutrons),
    NOT scattering — kept its own operator (physics-faithful; it also feeds the keff
    *production* numerator) and summed with :class:`IsotropicScattering` for the
    isotropic in-scatter source. Routes through
    :meth:`~orpheus.transport.mesh.material_xs_field.MaterialXSField.apply_n2n`
    (Pattern 2).

    Parameters
    ----------
    n2n : N2NMaterialField
        The :math:`(n,2n)` channel's kernel field (per-material reaction
        matrices over the mesh layout; the multiplicity enters from
        :attr:`~orpheus.transport.kernels.N2NKernel.multiplicity` inside
        the field's verbs).
    domain, codomain : FunctionSpace
        The two mandatory ends — see :class:`IsotropicScattering`.
    """

    n2n: "N2NMaterialField"
    block_role = BlockRole.BULK

    def __post_init__(self) -> None:
        # CS4a K2 → CS4c per-END form (one shared guard, per end).
        self._assert_energy_extent_both_ends(
            self.n2n.ng, operator="IsotropicN2N",
        )

    @classmethod
    def from_material_xs(
        cls, mat_xs: "MaterialXSField", *, space: "FunctionSpace",
    ) -> "IsotropicN2N":
        r"""Tier-2 extract-and-mint: the facade's :math:`(n,2n)` channel,
        endomorphic on one ``space=``."""
        return cls(
            N2NMaterialField.from_material_xs(mat_xs),
            domain=space,
            codomain=space,
        )

    @property
    def data_ng(self) -> int:
        """The bound data's group count (the assembly helpers' read)."""
        return self.n2n.ng

    @property
    def is_adjointable(self) -> bool:
        # apply_transpose realises 2Σ_{2n}χ (the group-flip A^T); caps ⊇
        # apply_transpose. is_invertible inherits base False.
        return True

    def apply(self, phi: "np.ndarray | FullField | object") -> "np.ndarray | FullField":
        r""":math:`2\Sigma_{2n}^{T}\phi` — the per-cell (n,2n) source.

        Bare ndarray / flux-carrier in → bare ndarray out; a scalar
        :class:`FullField` composite in → the closed scalar source
        composite out (#290 P4 — the shared
        :func:`_scalar_composite_source` arm over this kernel).
        """
        if isinstance(phi, FullField):
            return _scalar_composite_source(self, phi)
        arr = _values_of(phi)
        out = np.zeros_like(arr)
        self.n2n.add_emission(out, arr)
        return out

    def apply_transpose(self, chi: "np.ndarray | object") -> np.ndarray:
        r""":math:`2\Sigma_{2n}\chi` — the group-flip transpose.

        Bare-ndarray surface only (composite transpose lands with #281 —
        see :meth:`IsotropicScattering.apply_transpose`).
        """
        if isinstance(chi, FullField):
            raise TypeError(
                f"{type(self).__name__}.apply_transpose: composite "
                f"FullField transpose is not yet wired (lands with the "
                f"adjoint diffusion consumer, #281); pass the bulk "
                f"values."
            )
        arr = _values_of(chi)
        out = np.zeros_like(arr)
        self.n2n.add_emission_transpose(out, arr)
        return out

    # ── The assembly mode (stencil-assembly 2b) ────────────────────────

    @property
    def is_assemblable(self) -> bool:
        # Shared predicate — see :func:`_iso_is_assemblable`.
        return _iso_is_assemblable(self)

    def assemble(self) -> "SparseAssembledOperator":
        r""":math:`[2\,\Sigma_{2n}^{T}]` — cell-block-diagonal bulk emission
        through the production kernel (:func:`_assemble_iso_energy_operator`)."""
        return _assemble_iso_energy_operator(self)

    def dense_per_material(self) -> dict[int, np.ndarray]:
        r"""Per-material operator matrix :math:`2\Sigma_{2n}^{T}` — ``M @ φ == apply(φ)``.

        The storage-side oracle view, read directly off the stored
        :math:`\Sigma_{2n}` — see
        :meth:`IsotropicScattering.dense_per_material` for the oracle-pair
        rationale (production materialization goes through ``as_matrix``).
        """
        return {
            mid: kernel.emission_matrix()
            for mid, kernel in self.n2n.per_material.items()
        }
