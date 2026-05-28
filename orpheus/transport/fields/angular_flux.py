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
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from orpheus.numerics.field import Field
from orpheus.numerics.space import FunctionSpace

if TYPE_CHECKING:
    from orpheus.sn.geometry import SNMesh


__all__ = ["AngularFlux"]


@dataclass(frozen=True, eq=False, kw_only=True)
class AngularFlux(Field):
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
    :class:`~orpheus.transport.sources.per_ordinate_source.PerOrdinateSource`
    (same shape, different physical kind — source vs flux), and with
    any other Field subclass. Same-class arithmetic is closed.

    Cross-class same-units operations (e.g., AngularResidual −
    AngularSource per future Issue #201) require an explicit named
    composition, not direct ``-``.

    No coupled boundary or history (pure Field). The composite
    iteration state (flux + boundary + history) is
    :class:`~orpheus.transport.timed_full_field.TimedFullField`.
    """

    mesh: "SNMesh"

    # ── Construction validation ──────────────────────────────────────

    def __post_init__(self) -> None:
        super().__post_init__()
        expected = (self.mesh.quad.N, self.mesh.ng, self.mesh.nx, self.mesh.ny)
        if self.space.shape != expected:
            raise ValueError(
                f"AngularFlux: space.shape {self.space.shape!r} does "
                f"not match (mesh.quad.N, mesh.ng, mesh.nx, mesh.ny) "
                f"= {expected!r}"
            )

    # ── Algebra extensions (over Field) ──────────────────────────────

    def _check_partner(self, other: object) -> None:
        super()._check_partner(other)
        if self.mesh is not other.mesh:  # type: ignore[attr-defined]
            raise ValueError(
                "AngularFlux arithmetic across distinct SNMesh "
                "instances is forbidden — the field is mesh-bound."
            )

    # ── Construction factories ───────────────────────────────────────

    @classmethod
    def from_mesh(
        cls, values: NDArray, mesh: "SNMesh",
    ) -> "AngularFlux":
        r"""Construct from raw values + mesh, deriving the space."""
        space = FunctionSpace(
            name="angular_flux",
            shape=(mesh.quad.N, mesh.ng, mesh.nx, mesh.ny),
        )
        return cls(values=values, space=space, mesh=mesh)

    @classmethod
    def from_ndarray(
        cls, arr: NDArray, mesh: "SNMesh",
    ) -> "AngularFlux":
        r"""Test-ergonomic alias for :meth:`from_mesh`."""
        return cls.from_mesh(arr, mesh)

    @classmethod
    def zeros_for_sn_mesh(cls, mesh: "SNMesh") -> "AngularFlux":
        r"""Construct an all-zero :class:`AngularFlux` sized to ``mesh``."""
        N = mesh.quad.N
        return cls.from_mesh(
            np.zeros((N, mesh.ng, mesh.nx, mesh.ny)), mesh,
        )

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
        scalar_values = np.einsum("n,ngij->gij", weights, self.values)
        return ScalarFlux.from_mesh(scalar_values, self.mesh)

    # ── Metadata read-throughs ───────────────────────────────────────

    @property
    def N(self) -> int:  # noqa: N802 — matches AngularQuadrature.N
        r"""Number of angular ordinates."""
        return self.mesh.quad.N

    @property
    def ng(self) -> int:
        r"""Number of energy groups."""
        return self.mesh.ng

    @property
    def nx(self) -> int:
        r"""Number of cells along x."""
        return self.mesh.nx

    @property
    def ny(self) -> int:
        r"""Number of cells along y."""
        return self.mesh.ny
