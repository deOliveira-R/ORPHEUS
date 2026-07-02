r"""The abstract analysis / reconstruction operator **roles**.

Per Grand Report v3 §5.7 and §32.7, every Galerkin-style discretisation
factors as a **(reconstruction, analysis)** pair :math:`(R, M)` of linear
operators. The fine→coarse analysis :math:`M` produces coefficients in a
smaller space; the coarse→fine reconstruction :math:`R` synthesises a fine-space
field from those coefficients. The discretised operator is :math:`A_h = M\,A\,R`.

This module carries the two abstract **roles** only:

* :class:`AnalysisOperator` (:math:`M : V \to W`, the fine→coarse / measured side);
* :class:`ReconstructionOperator` (:math:`R : W \to V`, the coarse→fine synthesis).

The concrete realisation is the discrete frame: a
:class:`~orpheus.numerics.frame.FrameBase` binds a
:class:`~orpheus.numerics.basis.Basis` to a
:class:`~orpheus.numerics.measure.DiscreteMeasure` and emits the
:attr:`~orpheus.numerics.frame.FrameBase.analysis` face (an
:class:`AnalysisOperator`) and the
:attr:`~orpheus.numerics.frame.FrameBase.reconstruction` face (a
:class:`ReconstructionOperator`). Build the angular one with
``quadrature.angular_frame(L)``.

The discipline lives in the frame TYPE, not in an operator marker
=================================================================

The **discipline** — whether the *test* space equals the *trial* space (Galerkin)
or differs from it (Petrov-Galerkin) — is NOT a property of the analysis
operator. It is a kind of *frame*, carried by the frame TYPE
(:class:`~orpheus.numerics.frame.GalerkinFrame` vs
:class:`~orpheus.numerics.frame.PetrovGalerkinFrame`; GitHub #268). The earlier
``GalerkinProjection`` / ``PetrovGalerkinProjection`` marker ABCs that lived here
were retired with that move — they declared a discipline on the *role* when the
discipline belongs to the *frame*. A consumer that needs the discipline reads the
frame's type (or its strengthened promise, e.g. ``Π* = R`` on a
:class:`~orpheus.numerics.frame.GalerkinFrame`), not a mixin on the operator.

Cross-method consumers
----------------------

The analysis / reconstruction roles are **infrastructure**, not SN-only. Every
method that lifts an angular / energy / spatial axis between fine and coarse
representations builds on them, always through a frame:

* **SN scattering** — :math:`Y^* W` Galerkin analysis of the angular flux onto real
  spherical harmonics; the per-ordinate Pℓ source is :math:`R\,\Lambda\,M\,\psi`.
* **Spatial homogenisation** — flux-weighted, eigenvalue-consistent adjoint-weighted
  region reduction (a :class:`~orpheus.numerics.frame.PetrovGalerkinFrame`).
* **Energy condensation** — within-group spectrum-weighted coarse-group reduction
  (also Petrov-Galerkin).
* **PN solver** (future, §10) — moment-space analysis / lift.
* **Stochastic Galerkin** (future, §15A.7) — polynomial-chaos analysis on the
  random-input axis.

References
----------

* Galerkin / Petrov-Galerkin general framework: Brenner, S. C. and Scott, L. R.
  (2008). *The Mathematical Theory of Finite Element Methods*, 3rd ed. Springer. §3.4.
* Spherical-harmonic moment projection in transport: Bell, G. I. and Glasstone, S.
  (1970). *Nuclear Reactor Theory*. Van Nostrand Reinhold. §1.6.
* Energy condensation as a Petrov-Galerkin projection: Hébert, A. (2009).
  *Applied Reactor Physics*. Polytechnique. §6.2.
"""

from __future__ import annotations

from abc import ABC, abstractmethod

from orpheus.numerics.operator import (
    Codomain,
    Domain,
    LinearOperator,
)


__all__ = [
    "AnalysisOperator",
    "ReconstructionOperator",
]


# ─────────────────────────────────────────────────────────────────────
# Abstract roles
# ─────────────────────────────────────────────────────────────────────


class AnalysisOperator(LinearOperator[Domain, Codomain], ABC):
    r"""Abstract :math:`M : V \to W` — the fine→coarse (measured) side of a frame.

    Generic over the two carrier types ``[Domain, Codomain]`` because the role
    IS a map between distinct spaces (the measure space :math:`V` → the
    coefficient space :math:`W`) — the two-parameter shape is the role's
    definition, not a speculative generalisation. The numerics frame face
    (:class:`~orpheus.numerics.frame._FrameAnalysis`) inherits this BARE: it
    acts on raw ``ndarray`` — a carrier-agnostic numerics primitive that carries
    no ``Vector``-typed static carrier (same as the bare operator primitives,
    :class:`~orpheus.numerics.operator.DiagonalOperator` et al.) — so the
    :math:`V`/:math:`W` distinction lives in the runtime
    :attr:`domain`/:attr:`codomain` :class:`FunctionSpace`\\ s. A transport-layer
    realisation instead instantiates distinct typed carriers
    (e.g. ``AnalysisOperator[AngularFlux, HarmonicMomentFlux]``).

    The ABC carries no implementation — the concrete realisation is the
    :attr:`~orpheus.numerics.frame.FrameBase.analysis` face of a frame, which
    provides :meth:`apply` (and :meth:`apply_transpose`, so its Hilbert adjoint
    ``.H`` falls out of the metric-aware ``_AdjointOperator``). The discipline
    (Galerkin vs Petrov-Galerkin) is the frame's TYPE, not a property of this role.

    Adjointable by default (a working ``apply_transpose`` is part of the
    role's contract). A realisation that cannot cheaply apply the
    transpose MAY override :attr:`is_adjointable` to ``False`` — but
    then it forfeits the free Hilbert adjoint.

    Notes
    -----

    No domain / codomain validation here — the ABC stays abstract enough to host
    both disciplines. The concrete frame face carries its own shape contract.
    """

    @property
    def is_adjointable(self) -> bool:
        # The analysis role carries a working apply_transpose, so concrete
        # frame faces get the metric-aware ``.H`` for free. A realisation
        # that forfeits the free adjoint MUST override this to False.
        # (ReconstructionOperator is apply-only at the role level → its
        # is_adjointable inherits the base False; the frame FACE overrides.)
        return True

    @abstractmethod
    def apply(self, x: Domain, /) -> Codomain:
        ...


class ReconstructionOperator(LinearOperator[Domain, Codomain], ABC):
    r"""Abstract :math:`R : W \to V`. Sibling of :class:`AnalysisOperator`.

    Generic over ``[Domain, Codomain]`` for the same reason as
    :class:`AnalysisOperator` — reconstruction maps the coefficient space back
    into the fine space, an honestly two-typed role. The numerics frame face
    inherits it BARE (raw ``ndarray``); a transport realisation carries distinct
    typed flux/moment carriers.

    A reconstruction operator lifts coefficients from a coarse space :math:`W` back
    into the fine space :math:`V`. The canonical pair :math:`(R, M)` exposes both
    primitives: :class:`AnalysisOperator` carries :math:`M`,
    :class:`ReconstructionOperator` carries :math:`R`.

    The concrete realisation — the
    :attr:`~orpheus.numerics.frame.FrameBase.reconstruction` face of a frame —
    provides :meth:`apply` (and :meth:`apply_transpose`, so ``R.H`` falls out for
    free). The role is apply-only; the frame face adds the transpose and
    overrides :attr:`is_adjointable`.

    Per Grand Report v3 §5.7 — the Operator hierarchy. Pair this with
    :class:`AnalysisOperator` when shipping a discretisation: the :math:`(R, M)` pair
    IS the discretisation, modulo a metric correction on the operators' codomain
    spaces.
    """

    @abstractmethod
    def apply(self, coefficients: Domain, /) -> Codomain:
        ...
