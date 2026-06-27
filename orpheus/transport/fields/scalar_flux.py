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
  while :class:`SNMesh` lives at L3 (``orpheus.sn.mesh.augmented_mesh``).
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

Units (B.4 — declared as the ``UNITS`` class constant)
======================================================

:math:`[1/(\mathrm{cm^2 \cdot s})]` — areal angle-integrated flux,
:data:`~orpheus.numerics.units.SCALAR_FLUX_UNITS`. **eV-free**: a stored
flux is always integrated over an energy *bin* (a multigroup group, or a
Monte-Carlo tally bin), so :math:`\phi_g = \int_{E_g}\phi(E)\,dE` is
group-integrated by construction — the ``eV`` cancels. Continuous energy
lives in the cross-section data / collision kernel, not in this field
(so a CE-MC tally and an MG-deterministic solve share this signature).
Under View-G (issues #205 / #207) units are NOT a space property; they
are the role-leaf's ``UNITS`` constant, and the operator-side unit-gain
check gates composition at operator-construction time (#208). See
:mod:`orpheus.numerics.units` for the full convention.

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
from typing import TYPE_CHECKING, ClassVar

from numpy.typing import NDArray

from orpheus.numerics.units import SCALAR_FLUX_UNITS, Unit
from orpheus.transport.fields._bases import ScalarField
from orpheus.transport.fields._flux_role import FluxRole

if TYPE_CHECKING:
    from orpheus.sn.mesh.augmented_mesh import SNMesh


__all__ = ["ScalarFlux"]


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class ScalarFlux(FluxRole, ScalarField):
    r"""Scalar flux field :math:`\phi(\vec r, g)`.

    Parameters
    ----------
    values : NDArray
        Field values of shape ``(ng, nx, ny)`` (the principled
        layout — Issue #196 PR-INDEX-7).
    space : FunctionSpace
        The function space this flux lives on. Must satisfy
        ``space.shape == (mesh.ng, *mesh.spatial_shape)``. Construction
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

    #: The :class:`FunctionSpace` name for this leaf (preserves the
    #: pre-B.1 space identity). All mesh/shape/algebra/factory machinery
    #: is inherited from :class:`ScalarField` / :class:`BulkField`.
    _SPACE_NAME: ClassVar[str] = "scalar_flux"

    #: Dimensional identity (View-G, B.4): areal angle-integrated flux
    #: ``1/(cm²·s)`` (eV-free — see module docstring). Metadata, not the
    #: arithmetic gate. See :mod:`orpheus.numerics.units`.
    UNITS: ClassVar[Unit] = SCALAR_FLUX_UNITS

    # ── Selectors ────────────────────────────────────────────────────

    def at_group(self, g: int) -> NDArray:
        r"""Return the per-group slice ``values[g]``, shape ``(nx, ny)``."""
        return self.values[g]
