r"""Per-ordinate source field :math:`Q^{\rm aniso}(\vec r, \hat\Omega_n, g)`.

The L2 typed wrapper for an anisotropic / per-ordinate source. Same
storage shape as :class:`~orpheus.sn.angular_flux.AngularFlux` (``(N,
ng, nx, ny)``), but a structurally distinct physical quantity (source
density, not flux density). Carries the :math:`P_\ell \ge 1` Galerkin
reconstruction contribution AND any MMS-style per-ordinate external
source.

Migration status (Depth B step D-F)
====================================

Moved from ``orpheus.sn.sources.PerOrdinateSource`` to here. Inherits
:class:`~orpheus.numerics.field.Field`. The cross-class arithmetic
with :class:`~orpheus.transport.sources.isotropic_source.IsotropicSource`
is **PRESERVED** — the refined Issue #207 principle (recorded in the
conversation 2026-05-27) permits cross-class dunders when there is a
canonical subspace-containment relation. :class:`IsotropicSource`
lives in the subspace of :class:`PerOrdinateSource` where every
ordinate carries the same value; the canonical injection
``iso → 1 ⊗ iso`` (broadcast across the Ω axis) is applied inside
the dunder. The result type is the LARGER (containing) space:

.. code-block:: python

    Q_total = aniso + iso          # → PerOrdinateSource (commutative)
    Q_total = iso + aniso          # → PerOrdinateSource (same result)

Two named factories complement the dunder for Pattern 7
producer-side normalisation (``/sum_w`` baked in vs not):

* :meth:`from_isotropic(values, mesh)` — broadcast + ``/sum_w`` (the
  Pattern 7 entry point).
* :meth:`IsotropicSource.as_per_ordinate` — broadcast WITHOUT
  ``/sum_w`` (when caller has already done it).

Units (informational, not enforced at dunder level)
====================================================

Per-ordinate reaction-rate density:
:math:`[1/(\mathrm{cm^3 \cdot s \cdot sr \cdot eV})]`, differing
from :class:`IsotropicSource`'s units by a ``1/sr`` factor. The
cross-class dunder applies the structural injection (broadcast)
without ``/(4π sr)`` normalisation; the Pattern 7 ``/sum_w``
discipline is the caller's concern (apply before the dunder) or
goes through :meth:`from_isotropic` (named factory bakes it in).
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from orpheus.numerics.field import Field
from orpheus.numerics.space import FunctionSpace

if TYPE_CHECKING:
    from orpheus.sn.geometry import SNMesh


__all__ = ["PerOrdinateSource"]


@dataclass(frozen=True, eq=False, kw_only=True)
class PerOrdinateSource(Field):
    r"""Per-ordinate source field :math:`Q^{\rm aniso}(\vec r, \hat\Omega_n, g)`.

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
    Algebra is inherited from :class:`~orpheus.numerics.field.Field`.
    Same-class arithmetic is closed. Cross-class arithmetic — including
    with :class:`~orpheus.transport.fields.angular_flux.AngularFlux`
    (same shape, different physical kind) — is rejected by Field's
    Layer 1 class-identity gate.
    """

    mesh: "SNMesh"

    # ── Construction validation ──────────────────────────────────────

    def __post_init__(self) -> None:
        super().__post_init__()
        N = self.mesh.quad.N
        expected = (N, self.mesh.ng, self.mesh.nx, self.mesh.ny)
        if self.space.shape != expected:
            raise ValueError(
                f"PerOrdinateSource: space.shape {self.space.shape!r} "
                f"does not match (mesh.quad.N, mesh.ng, mesh.nx, mesh.ny) "
                f"= {expected!r}"
            )

    # ── Algebra extensions (over Field) ──────────────────────────────

    def _check_partner(self, other: object) -> None:
        super()._check_partner(other)
        if self.mesh is not other.mesh:  # type: ignore[attr-defined]
            raise ValueError(
                "PerOrdinateSource arithmetic across distinct SNMesh "
                "instances is forbidden — the field is mesh-bound."
            )

    def __add__(self, other):
        r"""Add a partner source — dispatches by partner type.

        * :class:`PerOrdinateSource` partner → within-class add via
          Field's inherited algebra.
        * :class:`IsotropicSource` partner → canonical subspace-
          containment injection: broadcast the iso operand across the
          Ω axis (the structural map ``iso → 1 ⊗ iso``) and add into
          this per-ordinate operand's space. Returns
          :class:`PerOrdinateSource`. Commutative — delegates to
          :meth:`IsotropicSource.__add__` for the broadcast-and-add
          (single source of truth for the iso injection logic).

        See module docstring for the principled justification (refined
        Issue #207: cross-class dunders permitted via canonical
        subspace containment).
        """
        from orpheus.transport.sources.isotropic_source import IsotropicSource
        if isinstance(other, IsotropicSource):
            # Delegate to IsotropicSource.__add__ — Pattern 2 (single
            # source of truth for the broadcast-and-add).
            return other + self
        # Within-class (and any other-type → NotImplemented) via Field.
        return super().__add__(other)

    # ── Construction factories ───────────────────────────────────────

    @classmethod
    def from_mesh(
        cls, values: NDArray, mesh: "SNMesh",
    ) -> "PerOrdinateSource":
        r"""Construct from raw values + mesh, deriving the space."""
        space = FunctionSpace(
            name="per_ordinate_source",
            shape=(mesh.quad.N, mesh.ng, mesh.nx, mesh.ny),
        )
        return cls(values=values, space=space, mesh=mesh)

    @classmethod
    def from_ndarray(
        cls, arr: NDArray, mesh: "SNMesh",
    ) -> "PerOrdinateSource":
        r"""Test-ergonomic alias for :meth:`from_mesh`."""
        return cls.from_mesh(arr, mesh)

    @classmethod
    def from_isotropic(
        cls, iso_values: NDArray, mesh: "SNMesh",
    ) -> "PerOrdinateSource":
        r"""Project an iso scalar source :math:`Q(\vec r, g)` to per-ordinate.

        Producer-side Pattern 7 normalisation: applies the ``/sum_w``
        projection and broadcasts across the :math:`N` ordinates.
        Returns a per-ordinate source whose every ordinate slice equals
        ``iso_values / sum_w``.

        This is the canonical entry for **external, user-supplied** iso
        scalar sources (fixed-source problems) that must enter the
        per-ordinate SN transport equation
        :math:`(\Omega\cdot\nabla + \sigma_t)\psi_n = Q/W + q_{\rm aniso,n}`.
        It is also the canonical migration path from the retired
        cross-class :meth:`IsotropicSource.__add__(PerOrdinateSource)`
        dunder (see module docstring).

        Parameters
        ----------
        iso_values : NDArray
            Iso scalar source, shape ``(ng, nx, ny)``.
        mesh : SNMesh
            Phase-space carrier.

        Returns
        -------
        PerOrdinateSource
            Per-ordinate source with ``Q/sum_w`` broadcast across N.

        See also
        --------
        :meth:`IsotropicSource.as_per_ordinate` — broadcast WITHOUT
            the ``/sum_w`` normalisation (use when caller has already
            divided by sum_w).
        """
        expected = (mesh.ng, mesh.nx, mesh.ny)
        if iso_values.shape != expected:
            raise ValueError(
                f"PerOrdinateSource.from_isotropic expects iso shape "
                f"(ng, nx, ny) = {expected}; got {iso_values.shape}"
            )
        sum_w = float(mesh.quad.weights.sum())
        N = mesh.quad.N
        per_ord_values = np.broadcast_to(
            (iso_values / sum_w)[None, :, :, :],
            (N, mesh.ng, mesh.nx, mesh.ny),
        ).copy()
        return cls.from_mesh(per_ord_values, mesh)

    # ── Selectors ────────────────────────────────────────────────────

    def at_ordinate(self, n: int) -> NDArray:
        r"""Return the per-ordinate slice ``values[n]``, shape ``(ng, nx, ny)``."""
        return self.values[n]

    # ── Metadata read-throughs ───────────────────────────────────────

    @property
    def N(self) -> int:  # noqa: N802 — matches AngularQuadrature.N
        r"""Number of angular ordinates."""
        return self.mesh.quad.N

    @property
    def ng(self) -> int:
        return self.mesh.ng

    @property
    def nx(self) -> int:
        return self.mesh.nx

    @property
    def ny(self) -> int:
        return self.mesh.ny
