r"""The §5.6 reaction-rate functional :math:`\langle \Sigma_x, \cdot\rangle`.

The reaction rate is the contraction of a flux against a reaction cross
section over the energy-group axis — the most fundamental scalar diagnostic
of reactor analysis:

.. math::

    r_x(\vec r) \;=\; \langle \Sigma_x, \phi\rangle
                \;=\; \sum_{g'} \Sigma_{x,g'}(\vec r)\,\phi_{g'}(\vec r) .

It is a :class:`~orpheus.numerics.functional.Functional` (field → scalar-field,
fiberwise over space) — a **co-vector** :math:`\langle \Sigma_x|`, the dual
companion of a flux :class:`~orpheus.numerics.vector.Vector`. The two named
instances are the **production rate** (``Σx = νΣf``) and the **absorption
rate** (``Σx = Σa``); the Rayleigh-quotient eigenvalue is their ratio
``k = ⟨νΣf,φ⟩ / ⟨Σa,φ⟩``.

Where this sits in the operator algebra
=======================================

:class:`ReactionRateFunctional` is the **row-factor of the fission rank-1
dyad**. Fission is the single-mode (rank-1 / ℓ=0) reaction kernel

.. math::

    F \;=\; |\chi\rangle\langle \nu\Sigma_f|
      \;=\; \texttt{outer}(\chi,\; \mathrm{ReactionRateFunctional}(\nu\Sigma_f)) ,

so :meth:`evaluate` (the contraction ``⟨νΣf, φ⟩``) **is** the production-rate
contraction the fission matvec performs — not a parallel description of it.
This is the rank-1 degenerate of the multi-mode scattering kernel
``R∘Λ∘M`` (the :class:`~orpheus.transport.frames.harmonic_frame.HarmonicFrame`
manages the analogous *stack* of co-vectors); see
:mod:`orpheus.transport.operators.integral_kernel_operator`.

Specialisation of :class:`InnerProductFunctional`
=================================================

:class:`ReactionRateFunctional` **specialises** the generic numerics co-vector
:class:`~orpheus.numerics.functional.InnerProductFunctional` (``⟨w, ·⟩``) by
carrying the weight as a domain-typed
:class:`~orpheus.transport.fields.cross_section_field.CrossSectionField` — which
brings ``.mesh``, the ``1/cm`` units, and the group/spatial shape — exactly as
:class:`~orpheus.transport.frames.harmonic_frame.HarmonicFrame` specialises
:class:`~orpheus.numerics.frame.GalerkinFrame`. The contraction is inherited
unchanged (the group axis, ``axis=0``); only the typed construction surface is
added.

References
----------

* Grand Report v3 §5.6 — the suffix law (Operator / Kernel / Functional);
  the reaction rate as the canonical Functional.
* Lewis, E.E. & Miller, W.F. (1993). *Computational Methods of Neutron
  Transport*. ANS. §1.1 — reaction rates and the fission source.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from typing import TYPE_CHECKING

from orpheus.numerics.functional import InnerProductFunctional
from orpheus.transport.fields.cross_section_field import CrossSectionField

if TYPE_CHECKING:
    from orpheus.numerics.space import Covector
    from orpheus.transport.fields.scalar_flux import ScalarFlux

__all__ = ["ReactionRateFunctional", "IntegratedReactionRate"]


@dataclass(frozen=True, eq=False, init=False, repr=False)
class ReactionRateFunctional(InnerProductFunctional):
    r"""The reaction-rate co-vector :math:`\langle \Sigma_x, \phi\rangle`.

    A domain-typed :class:`~orpheus.numerics.functional.InnerProductFunctional`
    whose weight is a reaction cross section
    :math:`\Sigma_x` carried as a
    :class:`~orpheus.transport.fields.cross_section_field.CrossSectionField`.
    The group axis (``axis=0``) is contracted; the spatial axes survive as the
    per-cell reaction-rate **density** (no volume measure folded in — the
    volume-integrated rate is a separate spatial-measure concern).

    Build with the production cross section for the fission emission density
    (``ReactionRateFunctional(CrossSectionField.fission_production_of(mat_xs))``) or with the
    absorption cross section for the absorption rate. The weight is a
    snapshot of ``cross_section.values`` at construction: build the
    functional again after a depletion or feedback update. (The fission
    operator's own production rate is a separate, cached gather; the two
    sources of per-cell νΣ_f are #648.)

    Parameters
    ----------
    cross_section : CrossSectionField
        The reaction cross section :math:`\Sigma_x` (``1/cm``), shape
        ``(ng, *spatial)``. Production = :math:`\nu\Sigma_f`; absorption =
        :math:`\Sigma_a`.
    """

    cross_section: CrossSectionField

    def __init__(self, cross_section: CrossSectionField) -> None:
        # Set the inherited (frozen) base fields directly — the GalerkinFrame
        # specialisation pattern. The weight IS the cross section's value array
        # (axis 0 = groups); ``cross_section`` is retained for the domain handle
        # (``.mesh``, units).
        object.__setattr__(self, "cross_section", cross_section)
        object.__setattr__(self, "weight", cross_section.values)
        object.__setattr__(self, "axis", 0)

    def __repr__(self) -> str:
        return f"ReactionRateFunctional({self.cross_section!r})"


@dataclass(frozen=True, eq=False, init=False, repr=False)
class IntegratedReactionRate:
    r"""The volume-integrated reaction rate :math:`\langle \Sigma_x, \phi\rangle_V`, the co-vector :math:`\flat_V \Sigma_x`.

    A :class:`~orpheus.numerics.functional.Functional` returning a **scalar**
    (not a scalar-field): the total reaction rate over the whole spatial domain,

    .. math::

        R_x(\phi) \;=\; \int_V \sum_{g'} \Sigma_{x,g'}(\vec r)\,\phi_{g'}(\vec r)\,dV
                  \;=\; \langle \Sigma_x, \phi\rangle_V
                  \;=\; \sum_i c_i\,\phi_i ,
        \qquad c = \flat_V \Sigma_x = G\,\Sigma_x .

    It IS the co-vector of its cross section on the cross section's space
    :math:`V`: :attr:`covector` is ``V.flat(Σx)`` (a
    :class:`~orpheus.numerics.space.Covector`, minted once), and the value on
    a flux is the canonical pairing. The space's metric is the product of
    the energy axis's counting measure (a theorem of the multigroup
    convention, :class:`~orpheus.numerics.axis.EnergyAxis`) and the spatial
    axis's cell volumes, so :math:`G\,\Sigma_x = \Sigma_x\,dV` and the
    pairing is :math:`\sum_g \sum_i \Sigma_{x,gi}\,\phi_{gi}\,V_i`. The
    space answers the measure (CS4b S3); this class reads no axis.

    It equals the volume integral of the per-cell reaction-rate density: the
    fibrewise co-vector :class:`ReactionRateFunctional` contracts the group
    axis and leaves a density over the cells, and integrating that density
    over the cells gives :math:`R_x(\phi)`. The two are one integral in two
    association orders, which is why they agree to round-off and not bit
    for bit.

    Until the space-membership unit U2 (2026-10-30) the body contracted the
    group axis through the per-cell density :class:`ReactionRateFunctional`
    and then summed against the last axis's weights, which assumed the
    spatial measure sat on the last axis. The two spellings re-associate one
    product and reduce in another order: `[M]` the spec's drift battery, 144
    of 202 calls bit-equal and the rest within 4.2e-16 relative
    (``scratch/membership_covector/spec/spec.md`` §U2.3).

    This is the canonical typed object for the **k-eigenvalue numerator and
    denominator** — ``k = R_{νΣf}(φ) / R_{Σa}(φ)`` (plus per-method terms such as
    the SN ``(n,2n)`` production or the diffusion leakage, which are explicit
    additive corrections, NOT reaction rates). The unweighted call is the
    **φ†=1 degenerate** of the adjoint-weighted bilinear
    :math:`\langle \phi^\dagger, M[\Sigma_x]\,\phi\rangle = \int_V \sum_g
    \phi^\dagger_g \Sigma_{x,g} \phi_g\,dV` — which is **LIVE** (P6, #281):
    pass ``adjoint=`` (the importance :math:`\varphi^*`, a
    :class:`~orpheus.transport.fields.scalar_flux.ScalarFlux`) to
    :meth:`evaluate` for the bilinear itself, the vector-channel worth
    numerator of the eigenvalue-consistent collapse (theorem T1 of the
    algebra of record, :mod:`orpheus.derivations.common.homogenization`).

    Parameters
    ----------
    cross_section : CrossSectionField
        The reaction cross section :math:`\Sigma_x`, a member of the space
        whose metric supplies the volume integral. Production =
        :math:`\nu\Sigma_f`; absorption = :math:`\Sigma_a`.
    """

    cross_section: CrossSectionField

    def __init__(self, cross_section: CrossSectionField) -> None:
        object.__setattr__(self, "cross_section", cross_section)

    @cached_property
    def covector(self) -> Covector:
        r"""The co-vector :math:`\flat_V \Sigma_x = G\,\Sigma_x` on :math:`V^{*}`, minted once per rate."""
        return self.cross_section.space.flat(self.cross_section)

    def evaluate(self, phi, *, adjoint: "ScalarFlux | None" = None) -> float:
        r"""Return :math:`\langle \Sigma_x, \phi\rangle_V` — or the bilinear.

        With ``adjoint=None`` (the default): the volume-integrated reaction
        rate, the pairing of :attr:`covector` with ``phi``. With ``adjoint=``
        an importance field :math:`\varphi^*`, the **adjoint-weighted
        bilinear**

        .. math::

            \big\langle \varphi^*\!\odot\Sigma_x,\; \phi\big\rangle_V
            \;=\; \int_V \sum_g \varphi^*_g\,\Sigma_{x,g}\,\phi_g\,dV ,

        the vector-channel worth numerator of the eigenvalue-consistent
        collapse (T1, :mod:`orpheus.derivations.common.homogenization`).
        The importance is folded into the DENSITY (it pairs with the test
        side, never the carrier) and the result is the flat of the folded
        density. The degenerate law is exact: a flat :math:`\varphi^* = 1`
        reproduces the unweighted call bit-for-bit (``1.0 * x == x`` in IEEE
        floats, and both route through one mint).

        Parameters
        ----------
        phi : Field or NDArray
            The flux, a member of the cross section's space: a typed member
            is compared by content on every run; a bare array must have the
            space's exact shape, checked under ``__debug__``.
        adjoint : ScalarFlux, optional
            The importance :math:`\varphi^*` (e.g.
            ``adjoint_solution.importance``), on the SAME space as the cross
            section (content equality, the σ↔geometry pairing).

        Raises
        ------
        ValueError
            If ``phi`` or ``adjoint`` is not a member of the cross section's
            space.
        """
        if adjoint is None:
            return self.covector.evaluate(phi)
        # The importance folds into the density (it pairs with the test
        # side, never the carrier); the space's own membership verb admits it.
        space = self.cross_section.space
        return space.flat(space.coordinates_of(adjoint) * self.cross_section.values).evaluate(phi)

    def __repr__(self) -> str:
        return f"IntegratedReactionRate({self.cross_section!r})"
