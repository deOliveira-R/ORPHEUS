r"""Scalar flux field on a discretized (group × spatial) phase space.

The L2 typed wrapper for

.. math::

    \phi(\vec r, g) = \int_{4\pi} \psi(\vec r, \hat\Omega, g) \, d\hat\Omega

(continuous form) or, in the discrete-ordinates form consumed by SN
solvers,

.. math::

    \phi_g(\vec r) = \sum_n w_n \, \psi_{n,g}(\vec r)

where :math:`w_n` is the angular quadrature weight on ordinate
:math:`\hat\Omega_n`.

Migration status (Depth B, step D-D)
====================================

This class moved from ``orpheus.sn.scalar_flux`` to
``orpheus.transport.fields.scalar_flux`` and now inherits from
:class:`~orpheus.numerics.field.Field` rather than carrying a
hand-coded dunder skeleton. The migration:

* Drops the six per-class hand-coded dunders (Cardinal Rule 2 —
  single source of truth; the algebra is now inherited from
  :class:`Field` via :func:`dataclasses.replace`).
* Adds the ``space: FunctionSpace`` field (mandatory, per the
  Field ABC contract).
* Keeps ``mesh: SNMesh`` as an additive field, but now annotated
  under ``TYPE_CHECKING`` because :class:`ScalarFlux` lives at L2
  while :class:`SNMesh` lives at L3 (``orpheus.sn.geometry``).
  Runtime-wise the ``mesh`` field is duck-typed; the import-linter's
  TYPE_CHECKING exemption keeps the layer contract clean.
* Preserves the existing strict semantics: arithmetic across two
  :class:`ScalarFlux` instances with DIFFERENT :class:`SNMesh`
  identities is forbidden, even when the meshes have matching
  shapes. This is enforced via a :meth:`_check_partner` override
  that adds the mesh-identity check on top of Field's class-and-
  space gate. See the docstring of :meth:`_check_partner` for the
  full enforcement story.
* Introduces :meth:`from_mesh` and :meth:`from_ndarray` classmethods
  for ergonomic 2-arg construction (the new dataclass is
  ``kw_only=True`` per Depth B plan §8 risk #1; positional
  construction is no longer available on the raw constructor).

Why "SN-bound for now"
======================

Conceptually :class:`ScalarFlux` is method-agnostic — a scalar flux
distribution on a discretized phase space is a concept shared by
SN, CP, MoC, and diffusion. But today the ONLY consumer is the SN
solver chain; CP / MoC / diffusion still carry their own field
representations. Per ``feedback_unify_after_two_instances``, the
abstract :class:`TransportMesh` Protocol is DEFERRED until a second
consumer arrives. Until then ``mesh: SNMesh`` is the most honest
annotation: ScalarFlux is SN-bound by usage, even though its
algebra is method-agnostic.

Units (informational, not yet enforced)
=======================================

:math:`[1/(\mathrm{cm^2 \cdot s \cdot eV})]` per energy group bin.
The per-bin energy density is absorbed into the cross-section
convention. Under View-G (issues #205 / #207) this label is NOT a
space property — it becomes the role-leaf's ``UNITS`` class constant
(Phase B of the field-vocabulary plan), and the operator-side unit-
gain check gates composition at operator-construction time (#208).

References
----------

* Lewis, E.E. & Miller, W.F. (1993). *Computational Methods of
  Neutron Transport*. American Nuclear Society. §1.2 — scalar /
  angular flux definitions.
* Depth B plan §3.3 (L2 field type spec), §6 step D-D
  (migration plan), §8 risk #1 (kw_only mitigation).
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


__all__ = ["ScalarFlux"]


@dataclass(frozen=True, eq=False, kw_only=True)
class ScalarFlux(Field):
    r"""Scalar flux field :math:`\phi(\vec r, g)`.

    Parameters
    ----------
    values : NDArray
        Field values of shape ``(ng, nx, ny)`` (the principled
        layout — Issue #196 PR-INDEX-7).
    space : FunctionSpace
        The function space this flux lives on. Must satisfy
        ``space.shape == (mesh.ng, mesh.nx, mesh.ny)``. Construction
        via :meth:`from_mesh` is the canonical path; direct kw-only
        construction is for callers that already hold a constructed
        space.
    mesh : SNMesh
        The SN discretisation handle. Carries the per-cell volumes,
        coordinate system, and angular quadrature reference that
        downstream operators read.

    Notes
    -----
    Algebra is inherited from :class:`~orpheus.numerics.field.Field`
    (dunders ``+``, ``-``, unary ``-``, scalar ``*``, scalar ``/``,
    plus diagnostics ``linf``, ``l2``, ``inner_product``, ``copy``).
    The :meth:`_check_partner` override adds the SN-specific
    mesh-identity check on top of Field's class-and-space gate.

    Per-group selectors (:meth:`at_group`) return ``np.ndarray``
    VIEWS into ``values`` — downstream callers must not mutate them.

    The ``mesh`` annotation is a ``TYPE_CHECKING``-guarded forward
    reference because :class:`SNMesh` lives at L3 and this class
    lives at L2. At runtime ``mesh`` is duck-typed.
    """

    mesh: "SNMesh"

    # ── Construction validation ──────────────────────────────────────

    def __post_init__(self) -> None:
        super().__post_init__()  # Field validates shape against space.
        expected = (self.mesh.ng, self.mesh.nx, self.mesh.ny)
        if self.space.shape != expected:
            raise ValueError(
                f"ScalarFlux: space.shape {self.space.shape!r} does not "
                f"match mesh phase-space shape (mesh.ng, mesh.nx, mesh.ny) "
                f"= {expected!r}"
            )

    # ── Algebra extensions (over Field) ──────────────────────────────

    def _check_partner(self, other: object) -> None:
        r"""Reject cross-mesh arithmetic on top of Field's class/space
        gate.

        Field's :meth:`~Field._check_partner` rejects on class
        identity (Layer 1) and space equality. This override adds
        the SN-specific mesh-identity check: two :class:`ScalarFlux`
        instances built on DIFFERENT :class:`SNMesh` instances are
        non-additive, even when the meshes have matching shapes.

        Rationale: an SNMesh carries per-cell geometric metadata
        (volumes, edge coordinates, BC tags) that two same-shape
        meshes may disagree on. Silently mixing fluxes across such
        meshes produces a physically meaningless result. The
        existing pre-Depth-B contract was strict on this; the
        migration preserves it.
        """
        super()._check_partner(other)
        if self.mesh is not other.mesh:  # type: ignore[attr-defined]
            raise ValueError(
                "ScalarFlux arithmetic across distinct SNMesh instances "
                "is forbidden — the field is mesh-bound."
            )

    # ── Construction factories ───────────────────────────────────────

    @classmethod
    def from_mesh(
        cls, values: NDArray, mesh: "SNMesh",
    ) -> "ScalarFlux":
        r"""Construct a :class:`ScalarFlux` from raw values + mesh,
        deriving the :class:`FunctionSpace` automatically.

        This is the ergonomic 2-arg constructor — the canonical
        migration path from the pre-Depth-B positional
        ``ScalarFlux(values, mesh)`` API. The space is built with
        ``name='scalar_flux'`` and ``shape=(mesh.ng, mesh.nx, mesh.ny)``.

        Parameters
        ----------
        values : NDArray
            Field values of shape ``(mesh.ng, mesh.nx, mesh.ny)``.
        mesh : SNMesh
            The SN phase-space carrier.

        Returns
        -------
        ScalarFlux
            Typed scalar-flux instance.
        """
        space = FunctionSpace(
            name="scalar_flux",
            shape=(mesh.ng, mesh.nx, mesh.ny),
        )
        return cls(values=values, space=space, mesh=mesh)

    @classmethod
    def from_ndarray(
        cls, arr: NDArray, mesh: "SNMesh",
    ) -> "ScalarFlux":
        r"""Test-ergonomic alias for :meth:`from_mesh`.

        Per Depth B plan §3.7, every typed field exposes
        ``from_ndarray(arr, mesh)`` as the migration path from the
        retired ``apply(np.ndarray)`` singledispatch handlers (D-I).
        For :class:`ScalarFlux` this is just an alias of
        :meth:`from_mesh`; subclasses with more complex construction
        (e.g. :class:`AngularFlux` with its boundary partner) carry
        a distinct ``from_ndarray`` body.
        """
        return cls.from_mesh(arr, mesh)

    # ── Selectors ────────────────────────────────────────────────────

    def at_group(self, g: int) -> NDArray:
        r"""Return the per-group slice ``values[g]``, shape ``(nx, ny)``."""
        return self.values[g]

    # ── Metadata read-throughs ───────────────────────────────────────

    @property
    def ng(self) -> int:
        r"""Energy group count (delegated to ``mesh.ng``)."""
        return self.mesh.ng

    @property
    def nx(self) -> int:
        r"""Spatial extent in x (delegated to ``mesh.nx``)."""
        return self.mesh.nx

    @property
    def ny(self) -> int:
        r"""Spatial extent in y (delegated to ``mesh.ny``)."""
        return self.mesh.ny
