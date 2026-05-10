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
  per-ordinate Pℓ source is :math:`R \Lambda M \psi`.
* **PN solver** (future, §10) — moment-space projection / lift.
* **Energy condensation** (future, §17) — Petrov-Galerkin coarse-energy
  cross-section reduction with within-group spectrum weighting.
* **Spatial homogenisation** (future, §18) — Petrov-Galerkin flux-volume
  weighted region reduction.
* **Stochastic Galerkin** (future, §15A.7) — polynomial-chaos basis
  projection on the random-input axis.
* **MC adjoint moments** (future) — variance reduction via response
  moments built against :math:`Y_\ell^m`.

Naming hierarchy
----------------

The naming is **deliberately three-level** so the type itself signals
which discipline a concrete projection follows. A reader of
``class HarmonicMomentProjection(GalerkinProjection)`` immediately
knows (a) it's a projection, (b) it follows Galerkin discipline (test
space = trial space), (c) the basis is real spherical harmonics — no
docstring needed for those facts.

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
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from orpheus.numerics.operator import (
    CAP_APPLY,
    CAP_APPLY_TRANSPOSE,
    LinearOperator,
    LinearOperatorMixin,
)

if TYPE_CHECKING:
    from orpheus.numerics.measure import DiscreteMeasure


__all__ = [
    "ProjectionOperator",
    "GalerkinProjection",
    "PetrovGalerkinProjection",
    "HarmonicMomentProjection",
    "HarmonicMomentReconstruction",
]


# ─────────────────────────────────────────────────────────────────────
# Abstract bases
# ─────────────────────────────────────────────────────────────────────


class ProjectionOperator(LinearOperatorMixin, ABC):
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

    No domain / range validation here — the ABC stays abstract enough
    to host both Galerkin and Petrov-Galerkin subclasses. The concrete
    classes :class:`HarmonicMomentProjection`, etc. carry their own
    shape contracts in their docstrings.
    """

    capabilities: frozenset[str] = frozenset({CAP_APPLY, CAP_APPLY_TRANSPOSE})

    @abstractmethod
    def apply(self, x: np.ndarray) -> np.ndarray:
        ...


class GalerkinProjection(ProjectionOperator, ABC):
    r"""Galerkin discipline: test space equals trial space.

    Invariants concrete subclasses MUST satisfy:

    * :math:`\Pi \, R = I` on the coefficient space, where :math:`R`
      is the paired reconstruction operator. This is the
      **idempotency-on-coefficients** invariant.
    * :math:`R \, \Pi` is the band-limited projector on :math:`V` —
      idempotent on :math:`V`.
    * :math:`\Pi^* = R` under the :math:`V`-inner-product (the test
      space's inner product). This is the **adjoint-pair** identity
      that distinguishes Galerkin from Petrov-Galerkin: a Galerkin
      projection's adjoint IS its reconstruction.

    The :meth:`assert_galerkin_idempotency` method materialises the
    first invariant for testing purposes.
    """

    def assert_galerkin_idempotency(
        self,
        reconstruction: "LinearOperator",
        sample_coefficients: np.ndarray,
        atol: float = 1e-12,
        rtol: float = 1e-12,
    ) -> None:
        r"""Assert :math:`\Pi \, R \, c = c` for sample coefficients.

        Apply the reconstruction then the projection; compare against
        the identity on the input coefficients. Used by tests of
        concrete Galerkin projections to validate the
        idempotency-on-coefficients invariant.

        Parameters
        ----------
        reconstruction : LinearOperator
            The paired reconstruction operator :math:`R`.
        sample_coefficients : np.ndarray
            Sample input :math:`c` in the coefficient space; the
            invariant is checked as :math:`\Pi(R(c)) = c`.
        atol, rtol : float
            Absolute / relative tolerance for the comparison.

        Raises
        ------
        AssertionError
            If :math:`\Pi(R(c))` does not match :math:`c` within
            tolerance — signals a violation of Galerkin idempotency.
        """
        out = self.apply(reconstruction.apply(sample_coefficients))
        np.testing.assert_allclose(
            out, sample_coefficients, atol=atol, rtol=rtol,
            err_msg=(
                "Galerkin idempotency Π R c = c violated; check the "
                "(test, trial) basis pair and the inner-product "
                "weighting in Π."
            ),
        )


class PetrovGalerkinProjection(ProjectionOperator, ABC):
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


# ─────────────────────────────────────────────────────────────────────
# Concrete: harmonic-moment projection / reconstruction (real SH)
# ─────────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class HarmonicMomentProjection(GalerkinProjection):
    r"""Discrete Galerkin projection on real spherical harmonics.

    The single-axis primitive

    .. math::
       :label: harmonic-moment-projection

       \phi_\ell^m
       \;=\; \sum_n w_n \, Y_\ell^m(\hat\Omega_n) \, \psi_n,

    mapping an angular field :math:`\psi : (N, \ldots) \to`
    moment field :math:`\phi : (L+1, 2L+1, \ldots)`. The remaining
    axes (cell, group, ...) broadcast through unchanged.

    Galerkin discipline: the basis :math:`\{Y_\ell^m\}` is used for
    BOTH the trial space (where :math:`\psi` lives via the
    addition-theorem reconstruction) and the test space (where the
    inner product is taken). The paired reconstruction is
    :class:`HarmonicMomentReconstruction`.

    Mathematical content
    --------------------

    Under the no-:math:`4\pi/(2\ell+1)`-prefactor convention used by
    :func:`orpheus.numerics.spherical_harmonics.evaluate_real_sh`, the
    addition theorem reads :math:`\sum_m Y_\ell^m(\Omega) Y_\ell^m(\Omega')
    = P_\ell(\Omega \cdot \Omega')`. The Galerkin idempotency
    :math:`\Pi \, R \, c = c` holds on the band-limited subspace
    :math:`\mathrm{span}\{Y_\ell^m : \ell \le L\}` provided the
    quadrature integrates :math:`Y_\ell^m Y_{\ell'}^{m'}` exactly
    for :math:`\ell + \ell' \le 2 L` (the Lebedev / level-symmetric
    rules' ``degree_of_exactness`` field tracks this).

    Cross-method use
    ----------------

    * **SN aniso scattering** (this work, Wave 1) —
      :math:`Q^{\ell\ge 1}_n = R \Lambda M \psi` builds the per-ordinate
      Pℓ source in :class:`orpheus.sn.scattering.ScatteringOperator`.
    * **PN solver** (future, §10) — same projection, just on the PN
      moment-space side of the streaming-coupling.

    Parameters
    ----------
    weights : np.ndarray, shape (N,)
        Quadrature weights :math:`w_n` from the angular cubature.
        Typically ``DiscreteMeasure.weights`` for a Lebedev /
        level-symmetric rule.
    Y : np.ndarray, shape (N, L+1, 2L+1)
        Spherical-harmonic table from
        :func:`orpheus.numerics.spherical_harmonics.evaluate_real_sh`.
    L : int
        Maximum harmonic order. Must match ``Y.shape[1] - 1``.

    Notes
    -----

    Implementation: ONE :func:`numpy.einsum` on the angular axis. No
    Python loop over :math:`n`, :math:`\ell`, or :math:`m` — the
    iteration is internal to the einsum. This is the core "metric
    knows its iterative structure" property the Wave 0 plan is
    organised around.
    """

    weights: np.ndarray              # (N,)
    Y: np.ndarray                    # (N, L+1, 2L+1)
    L: int

    def __post_init__(self) -> None:
        # Validate the array shapes match the declared L. Frozen
        # dataclass — uses object.__setattr__ trick if we wanted
        # derived fields, but here we only validate.
        if self.Y.shape != (
            self.weights.shape[0], self.L + 1, 2 * self.L + 1,
        ):
            raise ValueError(
                f"HarmonicMomentProjection: Y shape "
                f"{self.Y.shape} inconsistent with N={self.weights.shape[0]}, "
                f"L={self.L}; expected ({self.weights.shape[0]}, "
                f"{self.L + 1}, {2 * self.L + 1})."
            )

    @classmethod
    def from_measure(
        cls,
        measure: "DiscreteMeasure",
        L: int,
    ) -> "HarmonicMomentProjection":
        r"""Construct from a :class:`DiscreteMeasure` over :math:`S^2`.

        Convenience constructor that pulls the weights from the
        measure and builds the :math:`Y` table by calling
        :func:`evaluate_real_sh` on the measure's nodes.

        The measure's ``nodes`` MUST be a ``(N, 3)`` array of
        direction cosines :math:`(\mu_x, \mu_y, \mu_z)` per ordinate.
        """
        from orpheus.numerics.spherical_harmonics import evaluate_real_sh
        nodes = measure.nodes
        if nodes.ndim != 2 or nodes.shape[1] != 3:
            raise ValueError(
                f"HarmonicMomentProjection.from_measure expects "
                f"nodes of shape (N, 3); got {nodes.shape}."
            )
        Y = evaluate_real_sh(L, nodes[:, 0], nodes[:, 1], nodes[:, 2])
        return cls(weights=measure.weights, Y=Y, L=L)

    def apply(self, psi: np.ndarray) -> np.ndarray:
        r"""Apply :math:`\Pi \psi` along the leading angular axis.

        :math:`\psi` shape: ``(N, ...)`` — leading axis is the ordinate
        axis, trailing axes broadcast.
        Output shape: ``(L+1, 2L+1, ...)`` — leading two axes are the
        :math:`\ell, m` moment axes.
        """
        return np.einsum("n,nlm,n...->lm...", self.weights, self.Y, psi)

    def apply_transpose(self, x: np.ndarray) -> np.ndarray:
        r"""Apply :math:`\Pi^*` (the W-weighted adjoint).

        Under the :math:`L^2(S^2, dW)` inner product on the trial
        space, :math:`\Pi^* = R` (the addition-theorem reconstruction
        with the :math:`(2\ell+1)` factor). This method builds a
        transient :class:`HarmonicMomentReconstruction` and delegates
        to its :meth:`apply` — the Galerkin invariant
        :math:`\Pi^* = R` becomes the implementation.
        """
        two_l_plus_one = np.arange(2 * self.L + 1, step=2) + 1.0
        # arange(0, 2L+1, step=2) = [0, 2, 4, ..., 2L] — that's NOT (2l+1).
        # Build (2l+1) for l = 0, 1, ..., L:
        two_l_plus_one = (2.0 * np.arange(self.L + 1) + 1.0)
        R = HarmonicMomentReconstruction(Y=self.Y, two_l_plus_one=two_l_plus_one)
        return R.apply(x)


@dataclass(frozen=True)
class HarmonicMomentReconstruction(LinearOperatorMixin):
    r"""Adjoint reconstruction :math:`R` paired with
    :class:`HarmonicMomentProjection`.

    Maps moment-space coefficients to a per-ordinate angular field
    via the addition-theorem reconstruction with the
    :math:`(2\ell+1)` factor:

    .. math::
       :label: harmonic-moment-reconstruction

       q_n
       \;=\; \sum_\ell (2\ell+1) \sum_m Y_\ell^m(\hat\Omega_n)
              \, \phi_\ell^m,

    consuming a moment field :math:`\phi : (L+1, 2L+1, \ldots) \to`
    angular field :math:`q : (N, \ldots)`. The remaining axes
    broadcast through unchanged.

    Adjoint of :class:`HarmonicMomentProjection` under the
    W-weighted :math:`L^2(S^2)` inner product — the Galerkin
    discipline's :math:`\Pi^* = R` identity. This relationship is the
    structural content the no-:math:`4\pi/(2\ell+1)`-prefactor
    convention is designed to make literal.

    Parameters
    ----------
    Y : np.ndarray, shape (N, L+1, 2L+1)
        Spherical-harmonic table from
        :func:`orpheus.numerics.spherical_harmonics.evaluate_real_sh`
        — same table consumed by the projection.
    two_l_plus_one : np.ndarray, shape (L+1,)
        The addition-theorem factor :math:`2\ell+1` precomputed for
        :math:`\ell = 0, 1, \ldots, L`.

    Notes
    -----

    The capability set advertises :pydata:`CAP_APPLY` and
    :pydata:`CAP_APPLY_TRANSPOSE` (the latter is the
    :class:`HarmonicMomentProjection`'s :meth:`apply` modulo the
    :math:`W` weight — but to compute :math:`R^*` in isolation we
    would need the weights, which this dataclass does not carry. The
    canonical adjoint pairing flows through
    :class:`HarmonicMomentProjection` instead).

    Implementation: ONE :func:`numpy.einsum`. No Python loops over
    :math:`n`, :math:`\ell`, or :math:`m`.
    """

    Y: np.ndarray                    # (N, L+1, 2L+1)
    two_l_plus_one: np.ndarray       # (L+1,)
    capabilities: frozenset = frozenset({CAP_APPLY})

    def __post_init__(self) -> None:
        n_l = self.two_l_plus_one.shape[0]
        if self.Y.shape[1] != n_l:
            raise ValueError(
                f"HarmonicMomentReconstruction: Y has L+1 = {self.Y.shape[1]} "
                f"but two_l_plus_one has length {n_l}; mismatch."
            )

    @classmethod
    def from_Y(cls, Y: np.ndarray) -> "HarmonicMomentReconstruction":
        """Construct from a :math:`Y` table; compute :math:`(2\\ell+1)`."""
        L = Y.shape[1] - 1
        two_l_plus_one = (2.0 * np.arange(L + 1) + 1.0)
        return cls(Y=Y, two_l_plus_one=two_l_plus_one)

    def apply(self, moments: np.ndarray) -> np.ndarray:
        r"""Apply :math:`R \phi` along the leading harmonic-coefficient axes.

        ``moments`` shape: ``(L+1, 2L+1, ...)`` — leading two axes are
        the :math:`\ell, m` moment axes.
        Output shape: ``(N, ...)`` — leading axis is the ordinate axis.
        """
        return np.einsum(
            "nlm,l,lm...->n...", self.Y, self.two_l_plus_one, moments,
        )
