r"""Real-spherical-harmonic moment **source/sink** field :math:`q_\ell^m(\vec r, g)`.

The source/sink-role sibling of
:class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`:
a moment-space rate density rather than a flux. It is the missing leaf that
completes the ``(angular ⊗ moment) × (flux ⊗ source)`` carrier grid::

                M : AngularFlux ───────────────▶ HarmonicMomentFlux
   (project, role-preserving, Galerkin Π)                │
                                                          │ Λ : scatter, role-CHANGING,
                                                          │     axis-preserving  (Σ_{s,ℓ})
                                                          ▼
        AngularSourceSink ◀──────────────── HarmonicMomentSourceSink
                R : (reconstruct, role-preserving, Galerkin R)

The anisotropic in-scatter source :math:`S_{\ell\ge 1} = (1/W)\,R\,\Lambda\,M\,\psi`
threads through this type: :math:`\Lambda` (Legendre moment scattering
:math:`\Sigma_{s,\ell}`) maps a :class:`HarmonicMomentFlux` to a
:class:`HarmonicMomentSourceSink` — the **sole role-changing edge** of the
grid — and :math:`R` reconstructs it to an
:class:`~orpheus.transport.source_sinks.angular_source_sink.AngularSourceSink`,
role preserved. Before this leaf the flux→source role change had nowhere to
live and leaked to the scattering consumer as a raw ``np.ndarray``.

Why bare :class:`MomentField` (no :class:`~orpheus.transport.fields._flux_role.FluxRole`)
=========================================================================================

A source/sink is a **rate density** — a signed term in the balance
:math:`L\psi + C\psi - S\psi - F\psi = q` — and rate densities **add
vectorially** (gains accumulate; a source plus a sink is their net rate). So
its additive algebra is the plain vector-space algebra a generic
:class:`~orpheus.numerics.field.Field` inherits: ``source + source → source``
is CLOSED, with no affine/torsor gate and no displacement mint. This is
exactly the bare-vs-``FluxRole`` split that distinguishes
:class:`~orpheus.transport.source_sinks.angular_source_sink.AngularSourceSink`
(bare) from :class:`~orpheus.transport.fields.angular_flux.AngularFlux`
(``FluxRole``), and
:class:`~orpheus.transport.source_sinks.scalar_source_sink.ScalarSourceSink`
(bare) from :class:`~orpheus.transport.fields.scalar_flux.ScalarFlux`
(``FluxRole``). The flux state is an affine point; the source/sink is a
vector — different role, different algebra, same storage family.

All construction (the ``L`` / ``spatial_moments`` fields, the
``(L+1, 2L+1, ng, *spatial[, …])`` shape, the ``SphericalHarmonicSpace(L) ⊗
CellGroup`` :class:`~orpheus.numerics.space.TensorProductSpace`, the
:meth:`~orpheus.transport.fields._bases.MomentField.from_mesh_and_L` /
:meth:`~orpheus.transport.fields._bases.MomentField.zeros_for_mesh_and_L`
factories, the ``L``-match ``_check_partner``) is inherited from
:class:`~orpheus.transport.fields._bases.MomentField` — the shared
moment-space machinery lifted there when this second moment leaf arrived
(``feedback_unify_after_two_instances``; the "clean before extending" pass).
This leaf adds only its **role identity**: the rate-density :attr:`UNITS` and
a distinct cell-group space name.

Units (rate density, not flux)
==============================

:math:`[1/(\mathrm{cm^3 \cdot s})]` —
:data:`~orpheus.numerics.units.SCALAR_RATE_UNITS`, the angle-integrated rate
shared with :class:`~orpheus.transport.source_sinks.scalar_source_sink.ScalarSourceSink`
(the ``ℓ=0`` moment of a moment source IS a scalar source). This is the
``1/sr`` difference from the per-ordinate
:class:`~orpheus.transport.source_sinks.angular_source_sink.AngularSourceSink`
(:data:`~orpheus.numerics.units.ANGULAR_RATE_UNITS`): a moment is
angle-integrated, the quadrature weights carry the ``sr`` that cancels.
Same units, different class — the gate is class identity, not units.

Scope
=====

Constructed-and-algebra-verified; the production wiring of the
:math:`\Lambda : \text{flux} \to \text{source}` edge onto this type lands with
the :class:`~orpheus.transport.frames.harmonic_frame.HarmonicFrame` typed-face
migration (Frame campaign P4). Cross-class arithmetic with
:class:`HarmonicMomentFlux` (same shape, different role) is forbidden by
Field's Layer-1 class-identity gate — the legitimate route is
:math:`\Lambda` (scattering), the role-changing edge.

References
----------

* Lewis, E.E. & Miller, W.F. (1993). *Computational Methods of Neutron
  Transport*. ANS. §3.5 — spherical-harmonic moments; §3.2 — the in-scatter
  source.
* ``.claude/plans/frame_projection_machinery.md`` P4 — the carrier-grid
  completion.
* ``coding-elegance`` Pattern 2 (single source of truth — the lift),
  Pattern 4 (illegal states unrepresentable — the role/class gate).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

from orpheus.numerics.units import SCALAR_RATE_UNITS, Unit
from orpheus.transport.fields._bases import MomentField


__all__ = ["HarmonicMomentSourceSink"]


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class HarmonicMomentSourceSink(MomentField):
    r"""Real-spherical-harmonic moment source/sink :math:`q_\ell^m(\vec r, g)`.

    Parameters
    ----------
    values : NDArray
        Moment coefficients of shape ``(L+1, 2L+1, ng, *spatial[, …])`` (the
        same layout as :class:`HarmonicMomentFlux`; the leading two axes
        encode the truncation order ``L``).
    space : FunctionSpace
        The function space — canonically a
        :class:`~orpheus.numerics.space.TensorProductSpace`
        ``SphericalHarmonicSpace(L) ⊗ CellGroup``. Construction via
        :meth:`~orpheus.transport.fields._bases.MomentField.from_mesh_and_L`
        is the canonical path.
    mesh : SNMesh
        The SN phase-space carrier.
    L : int
        Maximum harmonic order retained.

    Notes
    -----
    Algebra is the plain vector-space algebra inherited from
    :class:`~orpheus.numerics.field.Field` (``source + source → source`` is
    closed — a rate density, NOT a flux state, so no
    :class:`~orpheus.transport.fields._flux_role.FluxRole` affine gate).
    Construction, the ``L``-match ``_check_partner``, and the shape contract
    are inherited from :class:`~orpheus.transport.fields._bases.MomentField`.
    Cross-class arithmetic with :class:`HarmonicMomentFlux` (same shape,
    different physical kind) is rejected by Field's Layer-1 class-identity
    gate — see the module docstring.
    """

    #: Distinct cell-group factor name (the source/sink space identity,
    #: keeping it non-``==`` to the flux leaf's ``"cell_group"`` — the
    #: class-identity gate is the arithmetic guard, this is for clarity).
    _CELL_GROUP_NAME: ClassVar[str] = "cell_group_source_sink"

    #: Dimensional identity: angle-integrated rate density ``1/(cm³·s)`` —
    #: :data:`~orpheus.numerics.units.SCALAR_RATE_UNITS`, shared with
    #: :class:`~orpheus.transport.source_sinks.scalar_source_sink.ScalarSourceSink`
    #: (the ``ℓ=0`` moment source IS a scalar source). Metadata, not the
    #: gate (class identity is).
    UNITS: ClassVar[Unit] = SCALAR_RATE_UNITS
