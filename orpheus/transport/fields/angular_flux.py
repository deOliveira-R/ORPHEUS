r"""Pure-Field :class:`AngularFlux` — bulk angular flux on (N, ng, nx, ny).

L2 typed wrapper for the angular flux :math:`\psi(\vec r, \hat\Omega, g)`
on the discrete-ordinates phase space, inheriting
:class:`~orpheus.numerics.field.Field`.

Migration status (Depth B step D-H.1)
======================================

Pre-D-H ``orpheus.sn.angular_flux.AngularFlux`` conflated THREE
concerns on a single class:

1. Bulk volumetric flux values (the legitimate Field role).
2. Boundary trace state (``self.boundary: BoundaryFlux``).
3. Iteration / time history (``self._history`` shift register).

Post-D-H this class is a **pure Field**:

* Bulk flux values + space + mesh (the Field role).
* No ``boundary`` attribute (boundary state lives on
  :class:`~orpheus.transport.timed_full_field.TimedFullField`).
* No ``_history`` attribute (history lives on TimedFullField).
* No ``from_flat_with_traces`` / ``to_flat_with_traces`` (those move
  to TimedFullField, where the trace-data target is structurally
  correct).

Symmetric with :class:`~orpheus.transport.fields.scalar_flux.ScalarFlux`
in all respects — both are pure typed Fields on their respective phase
spaces.

References
==========

* Depth B plan §3.3 (L2 field types — pure AngularFlux) and §3.8
  (TimedFullField container).
* Grand Report v3 §5.5 (Field hierarchy), §32.5 (Field primitive
  spec).
* ``coding-elegance`` Pattern 4 (illegal states unrepresentable) —
  bulk flux without coupled boundary state is now structurally
  representable; the composite state is the explicit named TimedFullField.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, ClassVar

import numpy as np

from orpheus.numerics.units import ANGULAR_FLUX_UNITS, Unit
from orpheus.transport.displacements.angular_displacement import AngularDisplacement
from orpheus.transport.fields._bases import AngularField
from orpheus.transport.fields._flux_role import FluxRole

if TYPE_CHECKING:
    from orpheus.sn.geometry import SNMesh


__all__ = ["AngularFlux"]


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class AngularFlux(FluxRole, AngularField):
    r"""L2 typed angular flux: pure :class:`Field` on (N, ng, nx, ny).

    Parameters
    ----------
    values : NDArray
        Field values of shape ``(N, ng, nx, ny)`` in the principled
        layout (Issue #196 PR-INDEX-5/7).
    space : FunctionSpace
        The function space. Must have ``shape == (mesh.quad.N, mesh.ng,
        mesh.nx, mesh.ny)``. Use :meth:`from_mesh` to derive
        automatically.
    mesh : SNMesh
        The SN phase-space carrier.

    Notes
    -----
    Algebra is inherited from :class:`~orpheus.numerics.field.Field`
    (dunders ``+``, ``-``, unary ``-``, scalar ``*``, scalar ``/``,
    plus diagnostics ``linf``, ``l2``, ``inner_product``, ``copy``).
    Field's Layer 1 class-identity gate rejects cross-class arithmetic
    with :class:`~orpheus.transport.fields.scalar_flux.ScalarFlux`
    (different shape, different physical kind), with
    :class:`~orpheus.transport.source_sinks.angular_source_sink.AngularSourceSink`
    (same shape, different physical kind — source vs flux), and with
    any other Field subclass. Same-class arithmetic is closed.

    Cross-class same-units operations (e.g., AngularResidual −
    AngularSourceSink per future Issue #201) require an explicit named
    composition, not direct ``-``.

    No coupled boundary or history (pure Field). The composite
    iteration state (flux + boundary + history) is
    :class:`~orpheus.transport.timed_full_field.TimedFullField`.
    """

    #: The :class:`FunctionSpace` name for this leaf (preserves the
    #: pre-B.1 space identity). All mesh/shape/algebra/factory machinery
    #: is inherited from :class:`AngularField` / :class:`BulkField`.
    _SPACE_NAME: ClassVar[str] = "angular_flux"

    #: The affine difference-space sibling minted by ``ψ ⊖ ψ`` (#208) — the
    #: tangent vector :class:`AngularDisplacement` (same space + metric, distinct
    #: role). See :class:`~orpheus.transport.fields._flux_role.FluxRole`.
    _DISPLACEMENT_CLS: ClassVar[type] = AngularDisplacement

    #: Dimensional identity (View-G, B.4): areal per-solid-angle flux
    #: density ``1/(cm²·s·sr)``. Metadata, NOT the arithmetic gate (class
    #: identity is) — see :mod:`orpheus.numerics.units`.
    UNITS: ClassVar[Unit] = ANGULAR_FLUX_UNITS

    # ── Reductions ───────────────────────────────────────────────────

    def integrate_angular(self) -> "object":
        r"""Reduce to the scalar flux :math:`\phi = \sum_n w_n \psi_n`.

        Returns a :class:`~orpheus.transport.fields.scalar_flux.ScalarFlux`
        of shape ``(ng, nx, ny)`` carrying the per-cell-per-group
        weighted angular sum. The L² metric on
        :class:`~orpheus.transport.fields.scalar_flux.ScalarFlux` is
        the cell-volume-weighted norm; the angular reduction here is
        a contraction over the leading ``N`` axis with the quadrature
        weights ``w_n``.
        """
        from orpheus.transport.fields.scalar_flux import ScalarFlux
        weights = self.mesh.quad.weights
        # values shape (N, ng, nx, ny); contract over leading N axis
        # with weights (N,) → (ng, nx, ny).
        scalar_values = np.einsum("n,ng...->g...", weights, self.values)
        return ScalarFlux.from_mesh(scalar_values, self.mesh)
