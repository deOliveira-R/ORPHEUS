r"""Scalar (angle-integrated) residual field :math:`r(\vec r, g)`.

The L2 typed wrapper for the **angle-integrated transport-balance
residual** — the scalar reduction of
:class:`~orpheus.transport.residuals.angular_residual.AngularResidual`,

.. math::

    r_g(\vec r) \;=\; \sum_n w_n \,
        \bigl[(L + C - S - F)\,\psi - q\bigr]_{n,g}(\vec r),

the per-cell-per-group defect of the balance after the angular
quadrature contraction. Same storage shape as
:class:`~orpheus.transport.fields.scalar_flux.ScalarFlux` and
:class:`~orpheus.transport.source_sinks.scalar_source_sink.ScalarSourceSink`
(``(ng, nx, ny)``), but a structurally distinct physical quantity: the
residual of an equation, not a flux and not an external source.

Why a distinct class
====================

:class:`ScalarResidual` and
:class:`~orpheus.transport.source_sinks.scalar_source_sink.ScalarSourceSink` share a
dimensional signature (both rate densities,
:math:`1/(\mathrm{cm^3 \cdot s \cdot eV})`) — yet they are DIFFERENT
classes. Field's Layer-1 class-identity gate
(:meth:`~orpheus.numerics.field.Field._check_partner`) therefore
rejects ``scalar_residual - scalar_source`` even though units align:
same units give permission to add, not meaning. The balance combining
a residual and a source goes through an explicit named composition
(planned ``IterationResidual.from_balance``, B.5 / Issue #201), never a
bare cross-class ``-``. See the module docstring of
:mod:`orpheus.transport.residuals.angular_residual` for the full
"dimensional sin" rationale.

Units (B.4 — declared as the ``UNITS`` class constant)
======================================================

Scalar reaction-rate-density defect:
:math:`[1/(\mathrm{cm^3 \cdot s})]`
(:data:`~orpheus.numerics.units.SCALAR_RATE_UNITS`) — the ``ScalarFlux``
units times an inverse length; eV-free per the binned-energy convention.
Under View-G the units are NOT a space property; they ARE the
role-leaf's ``UNITS`` class constant (class identity *is* units
identity). See :mod:`orpheus.numerics.units`.

References
==========

* ``.claude/plans/field_role_typing_view_g.md`` — Phase B, step B.3
  (new role leaves); B.5 (named composition + dimensional-sin rewire).
* Grand Report v3 §5.5 (Field hierarchy), §32.5 (Field primitive spec).
* ``coding-elegance`` Pattern 4 (illegal states unrepresentable).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

from orpheus.numerics.units import SCALAR_RATE_UNITS, Unit
from orpheus.transport.fields._bases import ScalarField


__all__ = ["ScalarResidual"]


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class ScalarResidual(ScalarField):
    r"""L2 angle-integrated transport-balance residual on ``(ng, nx, ny)``.

    Parameters
    ----------
    values : NDArray
        Field values of shape ``(ng, nx, ny)`` (the principled layout —
        Issue #196 PR-INDEX-7).
    space : FunctionSpace
        The function space. Must have ``shape == (mesh.ng, mesh.nx,
        mesh.ny)``. Use :meth:`from_mesh` to derive automatically.
    mesh : SNMesh
        The SN phase-space carrier.

    Notes
    -----
    A thin role leaf: all storage, validation, algebra, and the
    ``from_mesh`` / ``from_ndarray`` factories are inherited from
    :class:`~orpheus.transport.fields._bases.ScalarField` /
    :class:`~orpheus.transport.fields._bases.BulkField`. The leaf
    carries no residual-specific behaviour beyond its class identity,
    which keeps residual arithmetic from silently mixing with
    :class:`~orpheus.transport.source_sinks.scalar_source_sink.ScalarSourceSink`
    (same shape AND same units) or
    :class:`~orpheus.transport.fields.scalar_flux.ScalarFlux` (same
    shape, different units). Same-class arithmetic is closed.

    Like :class:`AngularResidual`, it exposes NO bespoke factory — a
    residual is *produced* by an operator balance, not constructed from
    thin air. Build via the inherited :meth:`from_mesh`.
    """

    #: The :class:`FunctionSpace` name for this leaf — distinct from
    #: ``"scalar_flux"`` / ``"scalar_source_sink"`` so the space identity
    #: tracks the role. All mesh/shape/algebra/factory machinery is
    #: inherited from :class:`ScalarField` / :class:`BulkField`.
    _SPACE_NAME: ClassVar[str] = "scalar_residual"

    #: Dimensional identity (View-G, B.4): scalar rate density
    #: ``1/(cm³·s)`` — :data:`~orpheus.numerics.units.SCALAR_RATE_UNITS`,
    #: shared with ``ScalarSourceSink`` (same units, different role → the gate
    #: is class identity, NOT units). Metadata, not the gate.
    UNITS: ClassVar[Unit] = SCALAR_RATE_UNITS
