r"""Scalar flux field on an SN phase space.

Issue #197 PR-TYPED-2 — the typed wrapper for the scalar flux
:math:`\phi(\vec r, g) = \int_{4\pi} \psi(\vec r, \hat\Omega, g)\,d\hat\Omega`,
or in discrete form
:math:`\phi_g(\vec r) = \sum_n w_n \psi_{n,g}(\vec r)`.

Frozen dataclass (``coding-elegance`` Pattern 4 — illegal states
unrepresentable): the underlying ``(ng, nx, ny)`` ndarray is wrapped
with its :class:`SNMesh` reference so shape mismatches surface at
construction time, not at the consumer of a downstream einsum.

Reads as the math via dunder arithmetic
(``coding-elegance`` Pattern 1 — match the algebra of the domain):

* ``phi_a + phi_b`` returns a new :class:`ScalarFlux`.
* ``alpha * phi`` returns a new :class:`ScalarFlux` for scalar ``alpha``.
* ``phi.at_group(g)`` returns the per-group ``(nx, ny)`` slice.

In this PR :class:`ScalarFlux` ALSO doubles as the container for the
isotropic source ``Q`` (same shape, same mesh) — the dedicated
:class:`IsotropicSource` type is introduced in PR-TYPED-3.  Operations
that mutate ``Q`` in place (``add_iso_source`` / ``add_n2n_source``)
continue to accept and modify ``np.ndarray`` directly, leaving the
:class:`ScalarFlux` immutability contract intact at the user-visible
boundary.

Units: :math:`[1/(\rm cm^2\,s\,eV)]` per energy group bin.  The per-bin
energy density is absorbed into the cross-section convention.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from .geometry import SNMesh


__all__ = ["ScalarFlux"]


@dataclass(frozen=True)
class ScalarFlux:
    r"""Scalar flux field :math:`\phi(\vec r, g)`.

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
    :class:`ScalarFlux` is **frozen** — the ndarray must not be mutated
    externally.  When in doubt about ownership, call
    :meth:`copy` and operate on the fresh instance.

    Per-group selectors (:meth:`at_group`) return ``np.ndarray`` VIEWS
    into ``values``; downstream callers must not mutate them.
    """

    values: np.ndarray
    mesh: "SNMesh"

    def __post_init__(self) -> None:
        expected = (self.mesh.ng, self.mesh.nx, self.mesh.ny)
        if self.values.shape != expected:
            raise ValueError(
                f"ScalarFlux expects shape (ng, nx, ny) = {expected}; "
                f"got {self.values.shape}"
            )

    # ── Dunder algebra (within-type) ──────────────────────────────────

    def __add__(self, other: "ScalarFlux") -> "ScalarFlux":
        self._validate_partner(other)
        return ScalarFlux(self.values + other.values, self.mesh)

    def __sub__(self, other: "ScalarFlux") -> "ScalarFlux":
        self._validate_partner(other)
        return ScalarFlux(self.values - other.values, self.mesh)

    def __mul__(self, scalar: float) -> "ScalarFlux":
        return ScalarFlux(self.values * float(scalar), self.mesh)

    __rmul__ = __mul__

    def __truediv__(self, scalar: float) -> "ScalarFlux":
        return ScalarFlux(self.values / float(scalar), self.mesh)

    def __neg__(self) -> "ScalarFlux":
        return ScalarFlux(-self.values, self.mesh)

    def _validate_partner(self, other: "ScalarFlux") -> None:
        if not isinstance(other, ScalarFlux):
            raise TypeError(
                "ScalarFlux arithmetic requires a ScalarFlux partner; "
                f"got {type(other).__name__}"
            )
        if other.mesh is not self.mesh:
            raise ValueError(
                "ScalarFlux arithmetic across distinct SNMesh instances "
                "is forbidden — the field is mesh-bound."
            )

    # ── Selectors ─────────────────────────────────────────────────────

    def at_group(self, g: int) -> np.ndarray:
        """Return the per-group slice ``values[g]``, shape ``(nx, ny)``."""
        return self.values[g]

    # ── Diagnostics ───────────────────────────────────────────────────

    def linf(self) -> float:
        """Return :math:`\\lVert\\phi\\rVert_\\infty` over all entries."""
        return float(np.abs(self.values).max())

    def copy(self) -> "ScalarFlux":
        """Return a deep copy carrying an owned ndarray."""
        return ScalarFlux(self.values.copy(), self.mesh)

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
