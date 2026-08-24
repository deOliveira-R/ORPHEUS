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
(:meth:`ScalarResidual.from_balance`, B.5 / Issue #201), never a
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
        The function space — the carrier's cached ``mesh.bulk_space``
        (CS4b S5).

    Notes
    -----
    A thin role leaf: all storage, validation, and algebra are
    inherited from
    :class:`~orpheus.transport.fields._bases.ScalarField` /
    :class:`~orpheus.transport.fields._bases.BulkField`. The leaf
    carries no residual-specific behaviour beyond its class identity,
    which keeps residual arithmetic from silently mixing with
    :class:`~orpheus.transport.source_sinks.scalar_source_sink.ScalarSourceSink`
    (same shape AND same units) or
    :class:`~orpheus.transport.fields.scalar_flux.ScalarFlux` (same
    shape, different units). Same-class arithmetic is closed.

    Like :class:`AngularResidual`, its one bespoke factory is the named
    composition :meth:`from_balance` — a residual is *produced* by an
    operator balance, not minted from thin air. Tests and direct
    producers may also build it space-primary on ``mesh.bulk_space``.
    """

    #: Dimensional identity (View-G, B.4): scalar rate density
    #: ``1/(cm³·s)`` — :data:`~orpheus.numerics.units.SCALAR_RATE_UNITS`,
    #: shared with ``ScalarSourceSink`` (same units, different role → the gate
    #: is class identity, NOT units). Metadata, not the gate.
    UNITS: ClassVar[Unit] = SCALAR_RATE_UNITS

    @classmethod
    def from_balance(
        cls, lhs: ScalarField, rhs: ScalarField,
    ) -> "ScalarResidual":
        r"""Construct the angle-integrated transport-balance residual.

        The residual leaf's bespoke factory: the scalar (angle-integrated)
        transport-balance defect ``r = lhs − rhs``, formed from two
        same-class
        :class:`~orpheus.transport.source_sinks.scalar_source_sink.ScalarSourceSink`
        operands. See :meth:`~orpheus.numerics.field.Field._from_balance` for
        the three guards (same-class operands, ``sr``-exact units, same space
        + mesh) and why the result lands on the ``"scalar_residual"`` space.
        """
        return cls._from_balance(lhs, rhs)
