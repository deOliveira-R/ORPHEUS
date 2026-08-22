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
from typing import TYPE_CHECKING

import numpy as np

from orpheus.numerics.functional import InnerProductFunctional
from orpheus.transport.fields.cross_section_field import CrossSectionField

if TYPE_CHECKING:
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
    (``ReactionRateFunctional(mat_xs.fission_production_field)``) or with the
    absorption cross section for the absorption rate. The fission operator
    builds it fresh per access (read-through to a depleting / feedback-updated
    :class:`MaterialXSField`), so a snapshot of ``cross_section.values`` as the
    base weight is faithful.

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
    r"""The volume-integrated reaction rate :math:`\int_V \langle \Sigma_x, \phi\rangle\,dV`.

    A :class:`~orpheus.numerics.functional.Functional` returning a **scalar**
    (not a scalar-field): the total reaction rate over the whole spatial domain,

    .. math::

        R_x(\phi) \;=\; \int_V \sum_{g'} \Sigma_{x,g'}(\vec r)\,\phi_{g'}(\vec r)\,dV
                  \;=\; \sum_{\text{cells}} V_{\text{cell}}\,
                        \langle \Sigma_x, \phi\rangle(\text{cell}) .

    It is the **volume integral of the per-cell density**
    :class:`ReactionRateFunctional` — so the group contraction has a *single
    source* (the density functional) and the spatial integral reuses the mesh's
    canonical ``volume_measure`` (the same measure
    :meth:`~orpheus.sn.solver.SNSolver.compute_group_production_rate` uses): no
    independent re-derivation of either reduction.

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
        The reaction cross section :math:`\Sigma_x` (carries ``.mesh``, whose
        ``volume_measure`` supplies the spatial integral). Production =
        :math:`\nu\Sigma_f`; absorption = :math:`\Sigma_a`.
    """

    density: ReactionRateFunctional

    def __init__(self, cross_section: CrossSectionField) -> None:
        # The per-cell group-contracted density is the single source of the
        # ``Σx·φ`` contraction; this functional adds only the volume integral.
        object.__setattr__(self, "density", ReactionRateFunctional(cross_section))

    @property
    def cross_section(self) -> CrossSectionField:
        return self.density.cross_section

    def evaluate(self, phi, *, adjoint: "ScalarFlux | None" = None) -> float:
        r"""Return :math:`\int_V \langle \Sigma_x, \phi\rangle\,dV` — or the bilinear.

        With ``adjoint=None`` (the default): the volume-integrated reaction
        rate, exactly as always.  With ``adjoint=`` an importance field
        :math:`\varphi^*`, the **adjoint-weighted bilinear**

        .. math::

            \int_V \big\langle \varphi^*\!\odot\Sigma_x,\; \phi\big\rangle\,dV
            \;=\; \int_V \sum_g \varphi^*_g\,\Sigma_{x,g}\,\phi_g\,dV ,

        the vector-channel worth numerator of the eigenvalue-consistent
        collapse (T1, :mod:`orpheus.derivations.common.homogenization`).
        The degenerate law is exact: a flat :math:`\varphi^* = 1` reproduces
        the unweighted call bit-for-bit (``1.0 * x == x`` in IEEE floats and
        the contraction body is shared).

        The importance is folded into the WEIGHT (:math:`\varphi^*` pairs
        with the test side, never the carrier), and the folded contraction
        routes through the same generic
        :class:`~orpheus.numerics.functional.InnerProductFunctional` body
        the density functional inherits — one contraction source, no
        parallel reduction.

        Parameters
        ----------
        phi : Field or NDArray
            The carrier flux, shape ``(ng, *spatial)``.
        adjoint : ScalarFlux, optional
            The importance :math:`\varphi^*` (e.g.
            ``adjoint_solution.importance``).  Must live on the SAME mesh
            object as the cross section (identity, not shape-equality —
            the σ↔geometry pairing tier).

        Raises
        ------
        ValueError
            If ``adjoint`` lives on a different mesh object than the
            cross section.
        """
        mesh = self.cross_section.mesh
        if adjoint is None:
            per_cell = np.asarray(self.density.evaluate(phi))  # (1, *spatial)
        else:
            if adjoint.space != self.cross_section.space:
                raise ValueError(
                    "IntegratedReactionRate.evaluate: the adjoint weight and "
                    "the cross section must agree in space content "
                    "(space-content invariant) — the σ↔geometry pairing; got "
                    f"{adjoint.space!r} vs {self.cross_section.space!r}."
                )
            folded = InnerProductFunctional(
                np.asarray(adjoint.values) * self.cross_section.values, axis=0,
            )
            per_cell = np.asarray(folded.evaluate(phi))       # (1, *spatial)
        # ``volume_measure`` consumes a flat ``(N_cells, k)`` view and contracts
        # the cell axis → ``(k,)``; here k = 1 (the density's collapsed group
        # axis), so the sum is the scalar volume integral.
        flat = np.moveaxis(per_cell, 0, -1).reshape(-1, 1)
        return float(mesh.volume_measure(flat).sum())

    def __repr__(self) -> str:
        return f"IntegratedReactionRate({self.cross_section!r})"
