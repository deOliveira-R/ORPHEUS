r"""Galerkin and Petrov-Galerkin projection / reconstruction primitives.

Per Grand Report v3 §5.7 and §32.7, every Galerkin-style discretisation
factors as a **(reconstruction, projection)** pair :math:`(R, P)` of
linear operators. The fine-to-coarse projection :math:`P` produces
coefficients in a smaller space; the coarse-to-fine reconstruction
:math:`R` synthesises a fine-space field from those coefficients. The
discretised operator is :math:`A_h = P \, A \, R`.

Two disciplines

* **Galerkin** — test space equals trial space; the canonical case
  for :math:`L^2`-orthogonal moment projection (e.g., real spherical
  harmonics of the SN angular flux). Invariant: :math:`P R = I` on
  the coefficient space, and :math:`R P` is the band-limited
  projector on the fine space.
* **Petrov-Galerkin** — test space differs from trial space; the
  canonical case for cross-section energy condensation (within-group
  spectrum weighting on the test side; identity on the trial side)
  and for spatial homogenisation (flux-volume weighting). Sibling of
  :class:`GalerkinProjection`; built when energy condensation lands.

Cross-method consumers
----------------------

The Galerkin / Petrov-Galerkin pair is **infrastructure**, not SN-only.
Every method that lifts an angular / energy / spatial axis between
fine and coarse representations builds on these primitives:

* **SN scattering** (this work, Wave 1) — :math:`Y^* W` Galerkin
  projection of the angular flux onto real spherical harmonics; the
  per-ordinate Pℓ source is :math:`R \Lambda M \psi`. Realised by the
  :class:`~orpheus.numerics.frame.Frame`'s analysis (:math:`M`) and
  reconstruction (:math:`R`) faces (``quadrature.angular_frame(L)``).
* **PN solver** (future, §10) — moment-space projection / lift.
* **Energy condensation** (future, §17) — Petrov-Galerkin coarse-energy
  cross-section reduction with within-group spectrum weighting.
* **Spatial homogenisation** (future, §18) — Petrov-Galerkin flux-volume
  weighted region reduction.
* **Stochastic Galerkin** (future, §15A.7) — polynomial-chaos basis
  projection on the random-input axis.
* **MC adjoint moments** (future) — variance reduction via response
  moments built against :math:`Y_\ell^m`.

Roles and disciplines (abstract vocabulary)
--------------------------------------------

This module carries the **abstract roles** and **disciplines** only; the
concrete spherical-harmonic realisation lives in
:class:`~orpheus.numerics.frame.Frame`.

* The two roles are :class:`AnalysisOperator` (:math:`\Pi : V \to W`, the
  fine→coarse analysis side) and :class:`ReconstructionOperator`
  (:math:`R : W \to V`, the coarse→fine synthesis side).
* The discipline — whether the **test** space equals the **trial** space —
  is declared by the analysis side inheriting from
  :class:`GalerkinProjection` (test = trial) or
  :class:`PetrovGalerkinProjection` (test ≠ trial).

The concrete realisation is the discrete frame: a
:class:`~orpheus.numerics.frame.Frame` binds a
:class:`~orpheus.numerics.basis.Basis` to a
:class:`~orpheus.numerics.measure.DiscreteMeasure` and emits the
:attr:`~orpheus.numerics.frame.Frame.analysis` face (the :math:`Y^* W`
Galerkin projection :math:`M`) and the
:attr:`~orpheus.numerics.frame.Frame.reconstruction` face (the
addition-theorem synthesis :math:`R = (2\ell+1)\,S_0`). Build one with
``quadrature.angular_frame(L)``; the §10 PN moment-space lift and the future
Petrov-Galerkin energy/spatial frames are the SAME ``Frame`` abstraction with
a different basis and measure. The "harmonic"-ness of the moments lives in the
frame's basis (:class:`~orpheus.numerics.basis.SphericalHarmonicBasis`) and its
coefficient space (:class:`~orpheus.numerics.spaces.SphericalHarmonicSpace`),
not in an operator class name.

References
----------

* Galerkin / Petrov-Galerkin general framework: Brenner, S. C. and
  Scott, L. R. (2008). *The Mathematical Theory of Finite Element
  Methods*, 3rd ed. Springer. §3.4.
* Spherical-harmonic moment projection in transport: Bell, G. I.
  and Glasstone, S. (1970). *Nuclear Reactor Theory*. Van Nostrand
  Reinhold. §1.6.
* Energy condensation as a Galerkin projection: Hébert, A. (2009).
  *Applied Reactor Physics*. Polytechnique. §6.2.
"""

from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np

from orpheus.numerics.operator import (
    CAP_APPLY,
    CAP_APPLY_TRANSPOSE,
    LinearOperatorMixin,
)


__all__ = [
    "AnalysisOperator",
    "GalerkinProjection",
    "PetrovGalerkinProjection",
    "ReconstructionOperator",
]


# ─────────────────────────────────────────────────────────────────────
# Abstract bases
# ─────────────────────────────────────────────────────────────────────


class AnalysisOperator(LinearOperatorMixin, ABC):
    r"""Abstract :math:`\Pi : V \to W`. Most general projection ABC.

    The ABC carries no implementation — concrete subclasses provide
    :meth:`apply`. The discipline (Galerkin vs Petrov-Galerkin) is
    declared by inheriting from :class:`GalerkinProjection` or
    :class:`PetrovGalerkinProjection`, which add the corresponding
    invariants.

    Attributes
    ----------
    capabilities : frozenset[str]
        ``{"apply", "apply_transpose"}`` by default. Concrete subclasses
        SHOULD provide :meth:`apply_transpose` (the paired
        :class:`Reconstruction`-shaped operator). If a subclass cannot
        cheaply apply the transpose, it MAY override ``capabilities``
        to drop ``CAP_APPLY_TRANSPOSE`` — but then the
        Galerkin invariants below cannot be checked symbolically.

    Notes
    -----

    No domain / codomain validation here — the ABC stays abstract enough
    to host both Galerkin and Petrov-Galerkin disciplines. The concrete
    realisation — the :attr:`~orpheus.numerics.frame.Frame.analysis` face of
    a :class:`~orpheus.numerics.frame.Frame` — carries its own shape contract.
    """

    capabilities: frozenset[str] = frozenset({CAP_APPLY, CAP_APPLY_TRANSPOSE})

    @abstractmethod
    def apply(self, x: np.ndarray, /) -> np.ndarray:
        ...


class GalerkinProjection(AnalysisOperator, ABC):
    r"""Galerkin discipline: test space equals trial space.

    Invariants the concrete realisation MUST satisfy (pinned, for the SH
    :class:`~orpheus.numerics.frame.Frame`, by
    ``tests/numerics/test_spherical_harmonic_space.py``; other concrete
    frames assert the same identities through their own test files):

    * :math:`\Pi \, R = (\sum_n w_n / |V|) \cdot I` on the band-limited
      coefficient space, where :math:`R` is the paired reconstruction
      and the proportionality constant arises from the basis's
      Gram-matrix normalisation. For the no-prefactor real-SH basis on
      :math:`S^2` (the
      :class:`~orpheus.numerics.spaces.SphericalHarmonicSpace`
      instance), this reduces to :math:`\Pi R = 4\pi \cdot I` —
      pinned by
      ``test_pi_R_is_4pi_identity_on_band_limited`` and the legacy
      ``test_pi_R_is_identity_on_band_limited``.
    * :math:`R \, \Pi` is the band-limited projector on :math:`V` —
      idempotent on :math:`V`.
    * :math:`\Pi^* = g_C \cdot S_0` where :math:`g_C` is the
      Gram-matrix diagonal on the coefficient space. This is the
      **adjoint-pair** identity that distinguishes Galerkin from
      Petrov-Galerkin: a Galerkin projection's Hilbert adjoint is the
      metric-weighted naked synthesis, NOT the addition-theorem
      reconstruction :math:`R = (2\ell+1) \cdot S_0` (the two differ
      by the diagonal :math:`g_C \cdot (2\ell+1) = 4\pi`). Pinned by
      ``test_H_equals_g_C_times_S0``.

    .. note::

       The earlier ``assert_galerkin_idempotency`` method that
       hand-checked ``Π R c = c`` (without the :math:`4\pi` factor)
       was retired in P1.6 of the moment-space + layering plan
       (ERR-039 follow-up). The method asserted the wrong invariant
       under the no-prefactor SH convention; the genuine
       :math:`4\pi \cdot I` identity is now pinned by the test files
       above.
    """


class PetrovGalerkinProjection(AnalysisOperator, ABC):
    r"""Petrov-Galerkin discipline: test space ≠ trial space.

    Sibling of :class:`GalerkinProjection`. The canonical use cases:

    * **Energy condensation** (§17) — within-group spectrum weighting
      on the test side, identity (or volume) on the trial side.
    * **Spatial homogenisation** (§18) — flux-volume weighting on the
      test side, identity on the trial side.

    For the most general Petrov-Galerkin pair :math:`(R, \Pi)` the
    invariant is :math:`\Pi \, R = I` on the coefficient space (same
    as Galerkin, by definition of "projection"). The DIFFERENCE is
    that :math:`\Pi^* \ne R` — the test space's inner product makes
    the adjoint distinct from the reconstruction.

    No concrete subclasses ship in this commit (the work that needs
    them — energy condensation, spatial homogenisation — is on
    follow-up branches). The ABC is included here so the symmetry of
    the two disciplines is visible to future readers, and the type
    signature is available to consumers.
    """


class ReconstructionOperator(LinearOperatorMixin, ABC):
    r"""Abstract :math:`R : W \to V`. Sibling of :class:`AnalysisOperator`.

    A reconstruction operator lifts coefficients from a coarse space
    :math:`W` back into the fine space :math:`V`. The canonical
    Galerkin / Petrov-Galerkin pair :math:`(R, \Pi)` exposes both
    primitives: :class:`AnalysisOperator` carries :math:`\Pi`,
    :class:`ReconstructionOperator` carries :math:`R`.

    The concrete realisation — the
    :attr:`~orpheus.numerics.frame.Frame.reconstruction` face of a
    :class:`~orpheus.numerics.frame.Frame` — provides :meth:`apply`. The
    default capability set advertises :pydata:`CAP_APPLY`; subclasses MAY add
    :pydata:`CAP_APPLY_TRANSPOSE` if they can cheaply expose the
    representation transpose.

    Per Grand Report v3 §5.7 — the Operator hierarchy. Pair this with
    :class:`AnalysisOperator` (also §5.7) when shipping a Galerkin
    discretisation: the (R, Π) pair IS the discretisation, modulo a
    metric correction that lives on the operators' codomain spaces.
    """

    capabilities: frozenset[str] = frozenset({CAP_APPLY})

    @abstractmethod
    def apply(self, coefficients: np.ndarray, /) -> np.ndarray:
        ...
