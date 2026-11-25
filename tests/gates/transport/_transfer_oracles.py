"""Test-side oracles of the transfer operators.

Three free functions that were methods of production classes with no production caller (retired from
``orpheus.transport.operators`` 2026-11-20); the gates that read them read them here.

* :func:`dense_per_material` is the STORAGE-SIDE oracle view of an isotropic transfer: the per-material operator
  matrix :math:`y\\,\\Sigma_{c,0}^{T}` (``[g_to, g_from]``), read DIRECTLY off the stored kernels through
  ``TransferKernel.emission_matrix``, structurally independent of the ``add_p0_source`` einsum that realizes ``apply``
  (``vv-principles`` L11: the two sides of an oracle pair must not share a realization). Each entry is a fresh copy.
* :func:`on_basis` and :func:`on_frame` are the frame-less and the frame-dressed tier-2 mints of
  :class:`~orpheus.transport.operators.transfer.LegendreMomentTransfer`. Production binds the moment factor through
  ``TransferOperator._moment_transfer`` on the frame's Parseval-dressed ``basis_space``.
"""

from __future__ import annotations

import numpy as np

from orpheus.numerics import frame_faces
from orpheus.numerics.basis.base import TruncatedBasis
from orpheus.numerics.frame import FrameBase
from orpheus.numerics.space import FunctionSpace
from orpheus.transport.material_field import TransferMaterialField
from orpheus.transport.operators.isotropic_transfer import IsotropicTransfer
from orpheus.transport.operators.transfer import LegendreMomentTransfer
from orpheus.numerics.factor import Factor, Kind, VarianceChart
from orpheus.numerics.space import TensorSpace
from orpheus.numerics.axis import Axis
from tests.gates.numerics._declared_support import whole_on_mod_kernel


def dense_per_material(isotropic_transfer: IsotropicTransfer) -> dict[int, np.ndarray]:
    r"""Per-material operator matrix :math:`y\,\Sigma_{c,0}^{T}` (``[g_to, g_from]``): ``{mid: M}`` with ``M @ phi_cell == apply(phi)_cell``."""
    return {
        mid: kernel.emission_matrix()
        for mid, kernel in isotropic_transfer.transfer.per_material.items()
    }


def on_basis(
    transfer: TransferMaterialField,
    basis: TruncatedBasis,
    *,
    skip_l0: bool = True,
) -> LegendreMomentTransfer:
    r"""Tier-2 mint: bring a channel's field to the basis's order and bind the endomorphic ends on the BASIS's own
    coefficient space (``FunctionSpace.of_basis(basis)`` supplying both). The basis names the FAMILY, so the ends are
    that basis's continuum-Gram space, never re-minted from ``L``."""
    ends = TensorSpace.of_axes(whole_on_mod_kernel(Axis.of_basis(basis), Kind.MEASURE, VarianceChart.DENSITY))
    return LegendreMomentTransfer(
        transfer.at_order(basis.L),
        skip_l0=skip_l0,
        domain=ends,
        codomain=ends,
    )


def on_frame(
    transfer: TransferMaterialField,
    frame: FrameBase,
    *,
    skip_l0: bool = True,
) -> LegendreMomentTransfer:
    r"""Tier-2 mint on the FRAME's Parseval-dressed coefficient space, the ends production binds, so the factor composes
    with the frame's faces (``conjugate(frame, Λ)``) by construction. The order is the frame's bound family's."""
    basis = frame.basis
    if not isinstance(basis, TruncatedBasis):
        raise TypeError(
            f"on_frame: the frame's trial {type(basis).__name__} carries no truncation order."
        )
    ends = frame_faces.basis_space(frame)
    return LegendreMomentTransfer(
        transfer.at_order(basis.L),
        skip_l0=skip_l0,
        domain=ends,
        codomain=ends,
    )
