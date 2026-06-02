r"""Boundary-trace source field :math:`q(\vec r_\Gamma, \hat\Omega, g)`.

The L2 typed wrapper for the **prescribed inflow source** materialised
as a stored field on the boundary trace — the :math:`q` term of the
affine boundary law

.. math::

    \gamma_-\psi \;=\; R\,G\,\gamma_+\psi \;+\; q ,

realised as concrete per-face values packed into the unified
:class:`~orpheus.numerics.spaces.trace_space.TraceSpace` flat layout
(``(layout.total_size,)``). It is the *source* role leaf of
:class:`~orpheus.transport.fields._bases.BoundaryField`, sibling to
:class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux` (flux)
and :class:`~orpheus.transport.residuals.boundary_residual.BoundaryResidual`
(residual).

Recipe → snapshot: relationship to :class:`InflowSourceSpec`
============================================================

There are TWO distinct objects for the inflow :math:`q`, related as
**recipe → snapshot** (NOT duplicates):

* :class:`~orpheus.geometry.boundary._source.InflowSourceSpec` (the
  geometry-layer **generator**, formerly named ``BoundarySource``):
  a lazy, per-face, ``runtime_checkable`` Protocol — ``evaluate(shape)
  -> ndarray`` — that produces inflow values on demand from a bare
  shape, with no trace/mesh handle (deliberately decoupled from the
  trace; the impls :class:`NoSource` / :class:`ConstantInflowSource`
  derive their output from construction-time data only). It is the
  *recipe* the affine boundary law / sweep consumes inline per face
  (:meth:`IncomingSourceOperator.apply`).
* :class:`BoundarySourceSink` (THIS class, the L2 transport **field**):
  the *eager, whole-boundary, mesh-bound, role-typed snapshot* — the
  materialised :math:`q` as a stored
  :class:`~orpheus.numerics.field.Field` on ``mesh.trace``.

Two renames produced the current name. First (B.3) the geometry
Protocol ``BoundarySource`` was renamed to ``InflowSourceSpec``,
freeing ``BoundarySource`` for this field leaf and completing the
``{Angular,Scalar,Boundary}×{Flux,Source,Residual}`` role grid
(issues #205 / #201). Then (B.5) the whole **source** column was
renamed ``Source`` → ``SourceSink`` so the name carries the signed
rate-density semantics — a *source* (production) and a *sink* (an
operator-loss output such as :math:`L\psi`) are the same quantity with
opposite sign, and the role-leaf type holds both. Hence
``BoundarySourceSink``.

Consumer status (role-grid completion, not yet consumer-driven)
===============================================================

As of this commit there is **no end-to-end consumer** of a stored
boundary-source field: every MMS case in the suite vanishes on the
incoming trace :math:`\Gamma_-` by construction (vacuum :math:`q\equiv
0`), and :func:`solve_sn_fixed_source` has no per-face inflow argument.
The prescribed-inflow stack (``PrescribedInflow`` →
``IncomingSourceOperator`` → :class:`InflowSourceSpec`) is fully built
but exercised only by isolated BC-realizer unit tests. Minting this
leaf is therefore a deliberate **role-grid completion** (the
``{Boundary}×{Source}`` cell reserved at
:class:`~orpheus.transport.fields._bases.BoundaryField`), so the
vocabulary is coherent and the type is ready the moment a
beam / incident-flux problem arrives.

The **recipe → snapshot bridge** ``BoundarySourceSink.from_spec(spec,
mesh)`` (materialise an :class:`InflowSourceSpec` onto the trace by
looping ``spec.evaluate(face_shape)`` per face and packing the flat
layout) is intentionally NOT added yet — per
``feedback_unify_after_two_instances`` it waits for the first real
consumer that both declares a non-trivial ``InflowSourceSpec`` AND
drives a sweep that consumes a typed boundary-source field (rather than
the current inline ``evaluate(shape)`` call). Until then the field is
built directly from known per-face arrays via the inherited
:meth:`~orpheus.transport.fields._bases.BoundaryField.from_face_arrays`.

Units (B.4 — declared as the ``UNITS`` class constant)
======================================================

The boundary trace is **all-flux**: the prescribed inflow ``q`` is a
flux added to :math:`\gamma_-\psi`, so it carries the angular-flux units
:math:`[1/(\mathrm{cm^2 \cdot s \cdot sr})]`
(:data:`~orpheus.numerics.units.ANGULAR_FLUX_UNITS`) — the SAME as
``BoundaryFlux`` and ``BoundaryResidual``. So on the trace, unlike the
bulk, a *source* does NOT pick up the volumetric ``cm⁻³``. eV-free per
the binned-energy convention. Same units, different role — the gate is
class identity. See :mod:`orpheus.numerics.units`.

References
==========

* ``.claude/plans/field_role_typing_view_g.md`` — Phase B step B.3.
* Grand Report v3 §16A.1 (affine boundary form), §5.5 (Field hierarchy).
* ``coding-elegance`` Pattern 2 (single source of truth — the
  ``InflowSourceSpec`` rename removes the same-name-two-concepts smell),
  Pattern 4 (illegal states unrepresentable — cross-role boundary
  arithmetic raises by class identity).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

from orpheus.numerics.units import ANGULAR_FLUX_UNITS, Unit
from orpheus.transport.fields._bases import BoundaryField


__all__ = ["BoundarySourceSink"]


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class BoundarySourceSink(BoundaryField):
    r"""L2 boundary-trace inflow source — the *source* role leaf of
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
    A thin role leaf: all storage, validation, algebra, per-face access
    (:meth:`face_view`), and the
    :meth:`~orpheus.transport.fields._bases.BoundaryField.zeros_on`
    / :meth:`~orpheus.transport.fields._bases.BoundaryField.from_face_arrays`
    factories are inherited from :class:`BoundaryField`. The leaf carries
    no source-specific behaviour beyond its class identity — which is
    exactly what Field's Layer-1 gate uses to keep boundary source, flux,
    and residual arithmetic from silently mixing. Note that all three
    boundary leaves share the SAME ``TraceSpace`` (``mesh.trace``); the
    space comparison would pass, so it is the **class** gate (not the
    space gate) that rejects ``BoundarySourceSink ± BoundaryFlux`` /
    ``BoundarySourceSink ± BoundaryResidual``. Same-class arithmetic is
    closed.

    See the module docstring for the recipe→snapshot relationship to the
    geometry-layer :class:`InflowSourceSpec` generator (formerly named
    ``BoundarySource``) and the deferred ``from_spec`` bridge.
    """

    #: Dimensional identity (View-G, B.4): the boundary is all-flux, so
    #: ``1/(cm²·s·sr)`` — :data:`~orpheus.numerics.units.ANGULAR_FLUX_UNITS`,
    #: shared with ``BoundaryFlux`` / ``BoundaryResidual`` (same units,
    #: different role → class gate). Metadata, not the gate.
    UNITS: ClassVar[Unit] = ANGULAR_FLUX_UNITS
