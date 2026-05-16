r"""Angular flux field on an SN phase space.

Issue #197 PR-TYPED-2 — the typed wrapper for
:math:`\psi(\vec r, \hat\Omega_n, g)` sampled on the SN phase space
(quadrature × spatial grid × energy groups).

Frozen dataclass (``coding-elegance`` Pattern 4 — illegal states
unrepresentable): the underlying ``(N, ng, nx, ny)`` ndarray is wrapped
with its :class:`SNMesh` reference so shape mismatches surface at
construction time.

Reads as the math via dunder arithmetic
(``coding-elegance`` Pattern 1 — match the algebra of the domain) and
a typed angular reduction (``coding-elegance`` Pattern 3 — named
intermediates):

* ``psi_a + psi_b`` returns a new :class:`AngularFlux`.
* ``alpha * psi`` returns a new :class:`AngularFlux`.
* ``psi.integrate_angular()`` reduces along the ordinate axis,
  returning a :class:`ScalarFlux` (:math:`\phi = \sum_n w_n \psi_n`).
* ``psi.at_ordinate(n)`` returns the per-ordinate ``(ng, nx, ny)``
  slice view.

Units: :math:`[1/(\rm cm^2\,s\,sr\,eV)]` per energy group bin.  The
per-bin energy density is absorbed into the cross-section convention.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from .geometry import SNMesh
    from .scalar_flux import ScalarFlux


__all__ = ["AngularFlux"]


@dataclass(frozen=True)
class AngularFlux:
    r"""Angular flux field :math:`\psi(\vec r, \hat\Omega_n, g)`.

    Parameters
    ----------
    values : np.ndarray
        Field values of shape ``(N, ng, nx, ny)`` in the principled
        layout (Issue #196 PR-INDEX-5/7).
    mesh : SNMesh
        The SN phase-space carrier — validates shape and supplies the
        canonical ``(N, ng, nx, ny)`` quadruple.
    """

    values: np.ndarray
    mesh: "SNMesh"

    def __post_init__(self) -> None:
        N = self.mesh.quad.N
        expected = (N, self.mesh.ng, self.mesh.nx, self.mesh.ny)
        if self.values.shape != expected:
            raise ValueError(
                f"AngularFlux expects shape (N, ng, nx, ny) = {expected}; "
                f"got {self.values.shape}"
            )

    # ── Reduction → ScalarFlux ────────────────────────────────────────

    def integrate_angular(self) -> "ScalarFlux":
        r"""Reduce along the ordinate axis: :math:`\phi = \sum_n w_n \psi_n`.

        Returns a :class:`~orpheus.sn.scalar_flux.ScalarFlux` of shape
        ``(ng, nx, ny)``.  The per-ordinate quadrature weights live on
        the mesh's :class:`AngularQuadrature`; the contraction is the
        canonical angular average that produces the scalar flux that
        the within-group scattering source consumes.
        """
        from .scalar_flux import ScalarFlux
        w = self.mesh.quad.weights
        return ScalarFlux(
            np.einsum("n,ngxy->gxy", w, self.values), self.mesh,
        )

    # ── Dunder algebra (within-type) ──────────────────────────────────

    def __add__(self, other: "AngularFlux") -> "AngularFlux":
        self._validate_partner(other)
        return AngularFlux(self.values + other.values, self.mesh)

    def __sub__(self, other: "AngularFlux") -> "AngularFlux":
        self._validate_partner(other)
        return AngularFlux(self.values - other.values, self.mesh)

    def __mul__(self, scalar: float) -> "AngularFlux":
        return AngularFlux(self.values * float(scalar), self.mesh)

    __rmul__ = __mul__

    def __truediv__(self, scalar: float) -> "AngularFlux":
        return AngularFlux(self.values / float(scalar), self.mesh)

    def __neg__(self) -> "AngularFlux":
        return AngularFlux(-self.values, self.mesh)

    def _validate_partner(self, other: "AngularFlux") -> None:
        if not isinstance(other, AngularFlux):
            raise TypeError(
                "AngularFlux arithmetic requires an AngularFlux partner; "
                f"got {type(other).__name__}"
            )
        if other.mesh is not self.mesh:
            raise ValueError(
                "AngularFlux arithmetic across distinct SNMesh instances "
                "is forbidden — the field is mesh-bound."
            )

    # ── Selectors ─────────────────────────────────────────────────────

    def at_ordinate(self, n: int) -> np.ndarray:
        """Return the per-ordinate slice ``values[n]``, shape ``(ng, nx, ny)``."""
        return self.values[n]

    # ── Diagnostics ───────────────────────────────────────────────────

    def linf(self) -> float:
        r"""Return :math:`\lVert\psi\rVert_\infty` over all entries."""
        return float(np.abs(self.values).max())

    def copy(self) -> "AngularFlux":
        """Return a deep copy carrying an owned ndarray."""
        return AngularFlux(self.values.copy(), self.mesh)

    # ── Metadata read-throughs ────────────────────────────────────────

    @property
    def N(self) -> int:  # noqa: N802 — matches AngularQuadrature.N
        """Number of angular ordinates."""
        return self.mesh.quad.N

    @property
    def ng(self) -> int:
        """Energy group count."""
        return self.mesh.ng

    @property
    def nx(self) -> int:
        """Spatial extent in x."""
        return self.mesh.nx

    @property
    def ny(self) -> int:
        """Spatial extent in y."""
        return self.mesh.ny
