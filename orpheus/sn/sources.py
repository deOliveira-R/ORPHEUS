r"""Source-density field types for SN transport.

Issue #197 PR-TYPED-3 — typed wrappers for the SN external + secondary
source-density fields:

* :class:`IsotropicSource` —
  :math:`Q(\vec r, g)`, shape ``(ng, nx, ny)``, the isotropic
  volumetric source density.  Aggregates the per-group P0 in-scatter +
  (n,2n) + fission contributions that every ordinate sees identically
  (the ``1/W`` factor is applied later by the sweep).
* :class:`PerOrdinateSource` —
  :math:`Q^{\rm aniso}(\vec r, \hat\Omega_n, g)`, shape
  ``(N, ng, nx, ny)``, the per-ordinate anisotropic source density.
  Carries the :math:`P_\ell \ge 1` Galerkin reconstruction plus the
  MMS external source for verification problems.

Both classes are frozen dataclasses; both validate shape against the
:class:`SNMesh` at construction time
(``coding-elegance`` Pattern 4 — illegal states unrepresentable).

The load-bearing pattern in this module
========================================

The cross-type :meth:`IsotropicSource.__add__` accepting a
:class:`PerOrdinateSource` operand handles the broadcast that used to
live as a procedural ``np.broadcast_to(Q_iso[None, :, :, :],
psi.shape).copy()`` + ``Q += Q_aniso`` pair inside
:meth:`ScatteringOperator.apply`.  Under this PR the call site reads as
math (``coding-elegance`` Pattern 1 — match the algebra of the
domain)::

    combined: PerOrdinateSource = iso_source + aniso_source

The :class:`IsotropicSource.__add__` resolves either operand type and
returns the more specific :class:`PerOrdinateSource` when the partner
is per-ordinate.  Within-type addition stays within type.

Why distinct types from :class:`ScalarFlux` / :class:`AngularFlux`
==================================================================

The source types have **the same storage shape** as
:class:`~orpheus.sn.scalar_flux.ScalarFlux` /
:class:`~orpheus.sn.angular_flux.AngularFlux` respectively but are
**different physical quantities**.  Cross-type addition between a flux
and a source is undefined (dimensional mismatch: flux carries
:math:`[1/(\mathrm{cm}^2\,\mathrm s\,\mathrm{sr}\,\mathrm{eV})]`,
source carries the SAME units only because the transport equation
folds the per-bin energy density into the cross-section convention,
but the algebraic role is distinct: flux is the unknown of the
operator equation; source is the right-hand side).  Keeping the types
distinct enforces this at the dunder-algebra layer — ``phi + Q`` does
not type-check; ``phi * sigma`` (a future reaction-rate wiring) does.

Units
=====

Both types carry units
:math:`[1/(\mathrm{cm}^3\,\mathrm s\,\mathrm{sr}\,\mathrm{eV})]` when
expressed as an angular density (the transport equation's
:math:`Q^{\rm aniso}` term consumes it in this form).  The transport
equation eats both types under the same convention; the
:class:`IsotropicSource` differs in dimensionality (rank-3 over phase
space, no ordinate axis) but not in units.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from .geometry import SNMesh


__all__ = ["IsotropicSource", "PerOrdinateSource"]


@dataclass(frozen=True)
class IsotropicSource:
    r"""Isotropic external + secondary source density :math:`Q(\vec r, g)`.

    Parameters
    ----------
    values : np.ndarray
        Field values of shape ``(ng, nx, ny)`` in the principled layout
        (Issue #196 PR-INDEX-7).
    mesh : SNMesh
        The SN phase-space carrier — validates shape and supplies the
        canonical ``(ng, nx, ny)`` triple.

    Notes
    -----
    Same storage shape as :class:`~orpheus.sn.scalar_flux.ScalarFlux`,
    different physical quantity.  Cross-type addition with
    :class:`~orpheus.sn.scalar_flux.ScalarFlux` is deliberately
    undefined (dimensional mismatch).

    The cross-type :meth:`__add__` accepting a
    :class:`PerOrdinateSource` partner is the load-bearing dunder that
    replaces the procedural ``np.broadcast_to + copy + add`` pattern
    that used to live inside :meth:`ScatteringOperator.apply`.
    """

    values: np.ndarray
    mesh: "SNMesh"

    def __post_init__(self) -> None:
        expected = (self.mesh.ng, self.mesh.nx, self.mesh.ny)
        if self.values.shape != expected:
            raise ValueError(
                f"IsotropicSource expects shape (ng, nx, ny) = "
                f"{expected}; got {self.values.shape}"
            )

    # ── Dunder algebra ────────────────────────────────────────────────

    def __add__(
        self, other: "IsotropicSource | PerOrdinateSource",
    ) -> "IsotropicSource | PerOrdinateSource":
        r"""Add an :class:`IsotropicSource` or :class:`PerOrdinateSource`.

        * ``iso_a + iso_b`` → :class:`IsotropicSource` (within type).
        * ``iso + aniso`` → :class:`PerOrdinateSource` (broadcast across
          the N axis, then add).  The broadcast happens via
          ``values[None, :, :, :] + aniso.values`` — numpy returns a
          freshly-allocated array, no ``.copy()`` needed.

        The cross-type case is the load-bearing dunder for PR-TYPED-3:
        it dissolves the procedural
        ``np.broadcast_to(Q_iso[None, ...], psi.shape).copy(); Q +=
        Q_aniso`` pattern into one line that reads as math.
        """
        if isinstance(other, PerOrdinateSource):
            self._validate_mesh(other)
            return PerOrdinateSource(
                self.values[None, :, :, :] + other.values, self.mesh,
            )
        if isinstance(other, IsotropicSource):
            self._validate_mesh(other)
            return IsotropicSource(self.values + other.values, self.mesh)
        return NotImplemented

    def __radd__(
        self, other: "PerOrdinateSource",
    ) -> "PerOrdinateSource":
        """Reflected add when ``other`` is a :class:`PerOrdinateSource`.

        Python tries ``PerOrdinateSource.__add__(IsotropicSource)``
        first; that path delegates back to this :class:`IsotropicSource`
        via ``return other + self`` (commutative by construction).  If
        Python falls through to ``IsotropicSource.__radd__`` we honour
        the same algebra.
        """
        if isinstance(other, PerOrdinateSource):
            return self + other
        return NotImplemented

    def __sub__(self, other: "IsotropicSource") -> "IsotropicSource":
        if not isinstance(other, IsotropicSource):
            return NotImplemented
        self._validate_mesh(other)
        return IsotropicSource(self.values - other.values, self.mesh)

    def __mul__(self, scalar: float) -> "IsotropicSource":
        return IsotropicSource(self.values * float(scalar), self.mesh)

    __rmul__ = __mul__

    def __truediv__(self, scalar: float) -> "IsotropicSource":
        return IsotropicSource(self.values / float(scalar), self.mesh)

    def __neg__(self) -> "IsotropicSource":
        return IsotropicSource(-self.values, self.mesh)

    def _validate_mesh(
        self, other: "IsotropicSource | PerOrdinateSource",
    ) -> None:
        if other.mesh is not self.mesh:
            raise ValueError(
                "Source arithmetic across distinct SNMesh instances "
                "is forbidden — the field is mesh-bound."
            )

    # ── Conversion ────────────────────────────────────────────────────

    def as_per_ordinate(self) -> "PerOrdinateSource":
        r"""Broadcast this isotropic source to a per-ordinate source.

        Returns a :class:`PerOrdinateSource` whose every ordinate slice
        equals ``self.values``.  Uses ``np.broadcast_to(...).copy()`` so
        the resulting array is writable and owns its data.

        Used when a consumer needs to merge the isotropic source into a
        per-ordinate buffer without invoking the cross-type
        :meth:`__add__` dunder (e.g., when the partner is a zero
        :class:`PerOrdinateSource`, the dunder works but is equivalent
        to a plain ``as_per_ordinate``).
        """
        N = self.mesh.quad.N
        target_shape = (N, *self.values.shape)
        return PerOrdinateSource(
            np.broadcast_to(self.values[None, :, :, :], target_shape).copy(),
            self.mesh,
        )

    # ── Diagnostics ───────────────────────────────────────────────────

    def linf(self) -> float:
        """Return :math:`\\lVert Q\\rVert_\\infty` over all entries."""
        return float(np.abs(self.values).max())

    def copy(self) -> "IsotropicSource":
        """Return a deep copy carrying an owned ndarray."""
        return IsotropicSource(self.values.copy(), self.mesh)

    # ── Metadata read-throughs ────────────────────────────────────────

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


@dataclass(frozen=True)
class PerOrdinateSource:
    r"""Per-ordinate source density :math:`Q^{\rm aniso}(\vec r, \hat\Omega_n, g)`.

    Parameters
    ----------
    values : np.ndarray
        Field values of shape ``(N, ng, nx, ny)`` in the principled
        layout (Issue #196 PR-INDEX-5/7).
    mesh : SNMesh
        The SN phase-space carrier — validates shape and supplies the
        canonical ``(N, ng, nx, ny)`` quadruple.

    Notes
    -----
    Same storage shape as :class:`~orpheus.sn.angular_flux.AngularFlux`,
    different physical quantity.  Cross-type addition with
    :class:`~orpheus.sn.angular_flux.AngularFlux` is deliberately
    undefined.

    Carries the :math:`P_\ell \ge 1` Galerkin reconstruction
    contribution and any MMS-style per-ordinate external source.
    """

    values: np.ndarray
    mesh: "SNMesh"

    def __post_init__(self) -> None:
        N = self.mesh.quad.N
        expected = (N, self.mesh.ng, self.mesh.nx, self.mesh.ny)
        if self.values.shape != expected:
            raise ValueError(
                f"PerOrdinateSource expects shape (N, ng, nx, ny) = "
                f"{expected}; got {self.values.shape}"
            )

    # ── Dunder algebra ────────────────────────────────────────────────

    def __add__(
        self, other: "IsotropicSource | PerOrdinateSource",
    ) -> "PerOrdinateSource":
        r"""Add a partner source — always returns :class:`PerOrdinateSource`.

        ``iso + aniso`` and ``aniso + iso`` are both routed through this
        method.  For the isotropic partner we delegate to
        :meth:`IsotropicSource.__add__` so the broadcast logic lives at
        ONE site (Pattern 2 — single source of truth).
        """
        if isinstance(other, IsotropicSource):
            # Delegate to IsotropicSource.__add__ — commutative by
            # construction.  Avoids duplicating the broadcast logic.
            return other + self  # type: ignore[return-value]
        if isinstance(other, PerOrdinateSource):
            self._validate_mesh(other)
            return PerOrdinateSource(self.values + other.values, self.mesh)
        return NotImplemented

    def __sub__(
        self, other: "PerOrdinateSource",
    ) -> "PerOrdinateSource":
        if not isinstance(other, PerOrdinateSource):
            return NotImplemented
        self._validate_mesh(other)
        return PerOrdinateSource(self.values - other.values, self.mesh)

    def __mul__(self, scalar: float) -> "PerOrdinateSource":
        return PerOrdinateSource(self.values * float(scalar), self.mesh)

    __rmul__ = __mul__

    def __truediv__(self, scalar: float) -> "PerOrdinateSource":
        return PerOrdinateSource(self.values / float(scalar), self.mesh)

    def __neg__(self) -> "PerOrdinateSource":
        return PerOrdinateSource(-self.values, self.mesh)

    def _validate_mesh(self, other: "PerOrdinateSource") -> None:
        if other.mesh is not self.mesh:
            raise ValueError(
                "PerOrdinateSource arithmetic across distinct SNMesh "
                "instances is forbidden — the field is mesh-bound."
            )

    # ── Selectors ─────────────────────────────────────────────────────

    def at_ordinate(self, n: int) -> np.ndarray:
        """Return the per-ordinate slice ``values[n]``, shape ``(ng, nx, ny)``."""
        return self.values[n]

    # ── Diagnostics ───────────────────────────────────────────────────

    def linf(self) -> float:
        r"""Return :math:`\lVert Q^{\rm aniso}\rVert_\infty`."""
        return float(np.abs(self.values).max())

    def copy(self) -> "PerOrdinateSource":
        """Return a deep copy carrying an owned ndarray."""
        return PerOrdinateSource(self.values.copy(), self.mesh)

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
