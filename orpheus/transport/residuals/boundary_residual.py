r"""Boundary-trace residual field :math:`r_\Gamma(\vec r_\Gamma, \hat\Omega, g)`.

The L2 typed wrapper for the **defect of the boundary balance** — the
residual of the affine boundary law

.. math::

    r_\Gamma \;=\; \gamma_-\psi \;-\; \bigl(R\,G\,\gamma_+\psi + q\bigr),

stored as a field on the unified
:class:`~orpheus.numerics.spaces.trace_space.TraceSpace` (flat
``(layout.total_size,)`` buffer). It is the *residual* role leaf of
:class:`~orpheus.transport.fields._bases.BoundaryField`, sibling to
:class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux` (flux)
and :class:`~orpheus.transport.source_sinks.boundary_source_sink.BoundarySourceSink`
(source).

Consumer: the named balance, NOT an operator output (B.5.2)
===========================================================

The boundary residual is the *defect* of the affine boundary law,
formed **only by an explicit balance** — never straight off an operator
output. B.5.2 (#208) settled this with the governing principle: *a
residual arises only when an operator output* ``Aψ`` *is compared
against something else (the source* ``b`` *) and the difference is
taken*; the operator output ITSELF is ``Aψ`` — a **source/sink**, not a
residual. Accordingly, every SN operator's ``.apply`` boundary
(``StreamingOperator`` non-zero; ``Collision`` / ``Scattering`` /
``Fission`` / :class:`SNBoundaryOperator` zero-or-coupling) is the
*source/sink* role leaf
:class:`~orpheus.transport.source_sinks.boundary_source_sink.BoundarySourceSink`
(mirroring the bulk's ``AngularSourceSink``). The GMRES flat residual
``b − Aψ`` is formed internally on the raveled vector (via
:meth:`TimedFullField.to_flat`) and is never typed as a field.

This leaf's consumer is therefore the **named composition**
:meth:`BoundaryResidual.from_balance` (minted B.5.1) — the boundary
counterpart of the (likewise not-yet-consumer-driven) bulk
:class:`~orpheus.transport.residuals.angular_residual.AngularResidual`.
The honest driver that will write ``BoundaryResidual.from_balance(
lhs=ψ.inflow, rhs=B·ψ.outflow + q)`` at the solver level is Wave O step
**O.2** (the ``L+C−S−F−B`` loss-operator driver). Until then the leaf is
minted, units-tagged, and ready, with no end-to-end consumer — exactly
the role-grid-completion status its bulk sibling holds.

The completed boundary role grid (B.5.2) mirrors the bulk exactly::

    .apply        →  BoundarySourceSink   (Aψ — a source/sink)
    .solve        →  BoundaryFlux         (the swept solution trace)
    from_balance  →  BoundaryResidual     (the defect — O.2 honest driver)

Units (B.4 — declared as the ``UNITS`` class constant)
======================================================

The boundary balance is a flux-matching equation, so its defect is
flux-typed: :math:`[1/(\mathrm{cm^2 \cdot s \cdot sr})]`
(:data:`~orpheus.numerics.units.ANGULAR_FLUX_UNITS`) — the SAME as
``BoundaryFlux`` / ``BoundarySourceSink``, and NOT the volumetric
``cm⁻³`` of the *bulk* residual (the trace is all-flux). eV-free per
the binned-energy convention. Same units, different role — the gate is
class identity. See :mod:`orpheus.numerics.units`.

References
==========

* ``.claude/plans/field_role_typing_view_g.md`` — Phase B step B.3
  (mint) and B.5 (the matvec operator-output retype + ``from_balance``).
* Grand Report v3 §16A.1 (affine boundary form), §5.5 (Field hierarchy).
* ``coding-elegance`` Pattern 4 (illegal states unrepresentable).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

from orpheus.numerics.units import ANGULAR_FLUX_UNITS, Unit
from orpheus.transport.fields._bases import BoundaryField


__all__ = ["BoundaryResidual"]


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class BoundaryResidual(BoundaryField):
    r"""L2 boundary-balance residual — the *residual* role leaf of
    :class:`~orpheus.transport.fields._bases.BoundaryField`.

    Parameters
    ----------
    values : NDArray
        Flat backing buffer, shape ``(space.layout.total_size,)``.
    space : TraceSpace
        The unified boundary
        :class:`~orpheus.numerics.spaces.trace_space.TraceSpace`
        (canonically ``mesh.trace``), carrying the per-geometry
        :class:`~orpheus.numerics.face_layout.FaceLayout`.
    mesh : SNMesh
        The SN phase-space carrier (cross-mesh-arithmetic guard).

    Notes
    -----
    A thin role leaf: all storage, validation, algebra, per-face access,
    and the ``zeros_on`` / ``from_face_arrays`` factories are
    inherited from :class:`BoundaryField`. The leaf carries no
    residual-specific behaviour beyond its class identity. All three
    boundary leaves share the SAME ``TraceSpace`` (``mesh.trace``), so it
    is the **class** gate (not the space gate) that rejects
    ``BoundaryResidual ± BoundaryFlux`` / ``BoundaryResidual ±
    BoundarySourceSink``. Same-class arithmetic is closed.

    See the module docstring for the already-computed matvec consumer and
    why the wiring (retype of the operator-output boundary) is the
    boundary half of the B.5 carve.
    """

    #: Dimensional identity (View-G, B.4): the trace is all-flux, so
    #: ``1/(cm²·s·sr)`` — :data:`~orpheus.numerics.units.ANGULAR_FLUX_UNITS`,
    #: shared with ``BoundaryFlux`` / ``BoundarySourceSink`` (same units,
    #: different role → class gate). Metadata, not the gate.
    UNITS: ClassVar[Unit] = ANGULAR_FLUX_UNITS

    @classmethod
    def from_balance(
        cls, lhs: BoundaryField, rhs: BoundaryField,
    ) -> "BoundaryResidual":
        r"""Construct the boundary-balance residual.

        The residual leaf's bespoke factory: the defect of the affine
        boundary law

        .. math::

            r_\Gamma \;=\; \gamma_-\psi \;-\; \bigl(R\,G\,\gamma_+\psi + q\bigr),

        formed as ``lhs − rhs`` from two same-class boundary operands (both
        flux-typed on the trace,
        :data:`~orpheus.numerics.units.ANGULAR_FLUX_UNITS`). See
        :meth:`~orpheus.numerics.field.Field._from_balance` for the three
        guards (same-class operands, ``sr``-exact units, same space + mesh);
        the result lands on the shared ``mesh.trace`` space.

        The factory is minted in B.5.1 so the leaf is ready; the matvec
        wiring that *feeds* it (the operator-output face-defect retype) is
        deferred to B.5.2 / #208 — see the module docstring.
        """
        return cls._from_balance(lhs, rhs)
