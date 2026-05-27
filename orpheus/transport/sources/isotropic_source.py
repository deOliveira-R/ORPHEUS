r"""Isotropic source field :math:`Q(\vec r, g)`.

The L2 typed wrapper for an isotropic external + secondary source
density. Same storage shape as
:class:`~orpheus.transport.fields.scalar_flux.ScalarFlux` (``(ng, nx,
ny)``), but a structurally distinct physical quantity (units of
reaction-rate density / source density, not flux density).

Migration status (Depth B step D-F)
====================================

Moved from ``orpheus.sn.sources.IsotropicSource`` to here. Inherits
:class:`~orpheus.numerics.field.Field`. The cross-class arithmetic
with :class:`~orpheus.transport.sources.per_ordinate_source.PerOrdinateSource`
is **PRESERVED** — the refined Issue #207 principle (recorded in the
conversation 2026-05-27) permits cross-class dunders when there is a
canonical subspace-containment relation between the operands' spaces.
:class:`IsotropicSource` lives in the subspace of
:class:`PerOrdinateSource` where every ordinate carries the same
value; the canonical injection ``iso → 1 ⊗ iso`` (broadcast across
the Ω axis) is the structural map the dunder applies during the
operation.

The math reads cleanly as

.. math::

    Q_{\text{total}} = Q_{\text{iso}} + Q_{\text{aniso}}

with the dunder using ``isinstance`` dispatch (or
:func:`functools.singledispatchmethod` — the formal version) to
recognize the containment direction:

.. code-block:: python

    # Within-class: returns IsotropicSource (Field's inherited algebra).
    Q_total = iso_a + iso_b

    # Cross-class with containment: returns the LARGER type
    # (PerOrdinateSource), broadcast inside the dunder.
    Q_total = iso + aniso          # → PerOrdinateSource
    Q_total = aniso + iso          # commutative; same result

Units (informational, not enforced at dunder level)
====================================================

Reaction-rate density per cell per group:
:math:`[1/(\mathrm{cm^3 \cdot s \cdot eV})]` for the isotropic form,
vs :math:`[1/(\mathrm{cm^3 \cdot s \cdot sr \cdot eV})]` for the
per-ordinate form. The dimensionality differs by ``1/sr``; physically
the iso → per-ordinate broadcast also requires a ``/(4π sr)``
normalisation. The dunder does NOT apply that conversion (the
broadcast is the pure structural injection); callers that need
Pattern 7 producer-side normalisation apply ``/sum_w`` explicitly
BEFORE the dunder (e.g. ``(iso / sum_w) + aniso``) — or use the
named factory :meth:`PerOrdinateSource.from_isotropic` which bakes
``/sum_w`` in. The choice between caller-applied ``/sum_w`` and the
named factory is a Pattern 7 discipline concern, not a dunder concern.
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
    from orpheus.transport.sources.per_ordinate_source import PerOrdinateSource


__all__ = ["IsotropicSource"]


@dataclass(frozen=True, eq=False, kw_only=True)
class IsotropicSource(Field):
    r"""Isotropic source field :math:`Q(\vec r, g)`.

    Parameters
    ----------
    values : NDArray
        Field values of shape ``(ng, nx, ny)`` in the principled layout
        (Issue #196 PR-INDEX-7).
    space : FunctionSpace
        The function space. Must have ``shape == (mesh.ng, mesh.nx,
        mesh.ny)``. Use :meth:`from_mesh` to derive automatically.
    mesh : SNMesh
        The SN phase-space carrier.

    Notes
    -----
    Algebra is inherited from :class:`~orpheus.numerics.field.Field`
    (dunders ``+``, ``-``, unary ``-``, scalar ``*``, scalar ``/``,
    plus diagnostics ``linf``, ``l2``, ``inner_product``, ``copy``).
    Field's Layer 1 class-identity gate rejects cross-class arithmetic
    with :class:`~orpheus.transport.fields.scalar_flux.ScalarFlux`
    (same shape, different physical kind), with
    :class:`PerOrdinateSource` (different shape AND different physical
    kind), and with any other Field subclass. Same-class arithmetic is
    closed.

    The cross-class ``iso + aniso`` dunder of the pre-D-F class is
    RETIRED — see module docstring for the named-composition migration.
    """

    mesh: "SNMesh"

    # ── Construction validation ──────────────────────────────────────

    def __post_init__(self) -> None:
        super().__post_init__()
        expected = (self.mesh.ng, self.mesh.nx, self.mesh.ny)
        if self.space.shape != expected:
            raise ValueError(
                f"IsotropicSource: space.shape {self.space.shape!r} does "
                f"not match (mesh.ng, mesh.nx, mesh.ny) = {expected!r}"
            )

    # ── Algebra extensions (over Field) ──────────────────────────────

    def _check_partner(self, other: object) -> None:
        super()._check_partner(other)
        if self.mesh is not other.mesh:  # type: ignore[attr-defined]
            raise ValueError(
                "IsotropicSource arithmetic across distinct SNMesh "
                "instances is forbidden — the field is mesh-bound."
            )

    def __add__(self, other):
        r"""Add a partner source — dispatches by partner type.

        * :class:`IsotropicSource` partner → within-class add via
          Field's inherited algebra.
        * :class:`PerOrdinateSource` partner → canonical subspace-
          containment injection: broadcast ``self.values`` across the
          Ω axis (the structural map ``iso → 1 ⊗ iso``) and add into
          the per-ordinate operand's space. Returns
          :class:`PerOrdinateSource`. Commutative — see
          :meth:`PerOrdinateSource.__add__`.

        See module docstring for the principled justification (refined
        Issue #207: cross-class dunders permitted via canonical
        subspace containment).
        """
        from orpheus.transport.sources.per_ordinate_source import (
            PerOrdinateSource,
        )
        if isinstance(other, PerOrdinateSource):
            if self.mesh is not other.mesh:
                raise ValueError(
                    "IsotropicSource + PerOrdinateSource across "
                    "distinct SNMesh instances is forbidden — both "
                    "fields are mesh-bound."
                )
            # Canonical injection: broadcast iso → per-ordinate, add.
            return replace(
                other, values=self.values[None, :, :, :] + other.values,
            )
        # Within-class (and any other-type → NotImplemented) via Field.
        return super().__add__(other)

    def __radd__(self, other):
        r"""Reflected add — Python falls here when
        ``PerOrdinateSource.__add__(IsotropicSource)`` returns
        NotImplemented. Delegates to :meth:`__add__` (commutative).
        """
        from orpheus.transport.sources.per_ordinate_source import (
            PerOrdinateSource,
        )
        if isinstance(other, PerOrdinateSource):
            return self.__add__(other)
        return NotImplemented

    # ── Construction factories ───────────────────────────────────────

    @classmethod
    def from_mesh(
        cls, values: NDArray, mesh: "SNMesh",
    ) -> "IsotropicSource":
        r"""Construct from raw values + mesh, deriving the space."""
        space = FunctionSpace(
            name="isotropic_source",
            shape=(mesh.ng, mesh.nx, mesh.ny),
        )
        return cls(values=values, space=space, mesh=mesh)

    @classmethod
    def from_ndarray(
        cls, arr: NDArray, mesh: "SNMesh",
    ) -> "IsotropicSource":
        r"""Test-ergonomic alias for :meth:`from_mesh`."""
        return cls.from_mesh(arr, mesh)

    # ── Conversion (named cross-class composition) ───────────────────

    def as_per_ordinate(self) -> "PerOrdinateSource":
        r"""Broadcast this isotropic source to a per-ordinate source.

        Returns a :class:`PerOrdinateSource` whose every ordinate slice
        equals ``self.values``. Uses ``np.broadcast_to(...).copy()`` so
        the resulting array is writable and owns its data.

        Use this when the producer-side normalisation (``/sum_w``,
        Pattern 7) has ALREADY been applied to ``self.values`` and the
        consumer just needs the per-ordinate broadcast. When the
        normalisation has NOT been applied yet, use
        :meth:`PerOrdinateSource.from_isotropic` instead.
        """
        from orpheus.transport.sources.per_ordinate_source import (
            PerOrdinateSource,
        )
        N = self.mesh.quad.N
        target_shape = (N, *self.values.shape)
        per_ord_values = np.broadcast_to(
            self.values[None, :, :, :], target_shape,
        ).copy()
        return PerOrdinateSource.from_mesh(per_ord_values, self.mesh)

    # ── Metadata read-throughs ───────────────────────────────────────

    @property
    def ng(self) -> int:
        return self.mesh.ng

    @property
    def nx(self) -> int:
        return self.mesh.nx

    @property
    def ny(self) -> int:
        return self.mesh.ny
