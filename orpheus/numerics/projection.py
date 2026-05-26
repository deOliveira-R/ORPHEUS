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
``class MomentProjection(GalerkinProjection)`` immediately knows
(a) it's a projection, (b) it follows Galerkin discipline (test space
= trial space). The "harmonic"-ness of the moments is carried by the
typed codomain :class:`~orpheus.numerics.spaces.SphericalHarmonicSpace`
exposed via :attr:`MomentProjection.codomain`, rather than baked into
the operator class name — so the §10 PN moment-space projection can
land as a sibling under the same name with a different codomain.

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
from functools import cached_property
from typing import TYPE_CHECKING

import numpy as np

from orpheus.numerics.operator import (
    CAP_APPLY,
    CAP_APPLY_TRANSPOSE,
    LinearOperator,
    LinearOperatorMixin,
)
from orpheus.numerics.space import FunctionSpace

if TYPE_CHECKING:
    from orpheus.numerics.measure import DiscreteMeasure
    from orpheus.numerics.spaces.spherical_harmonic_space import (
        SphericalHarmonicSpace,
    )


__all__ = [
    "ProjectionOperator",
    "GalerkinProjection",
    "PetrovGalerkinProjection",
    "ReconstructionOperator",
    "MomentProjection",
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
    classes :class:`MomentProjection`, etc. carry their own
    shape contracts in their docstrings.
    """

    capabilities: frozenset[str] = frozenset({CAP_APPLY, CAP_APPLY_TRANSPOSE})

    @abstractmethod
    def apply(self, x: np.ndarray) -> np.ndarray:
        ...


class GalerkinProjection(ProjectionOperator, ABC):
    r"""Galerkin discipline: test space equals trial space.

    Invariants concrete subclasses MUST satisfy (pinned by
    ``tests/numerics/test_spherical_harmonic_space.py`` for the SH
    instance; other concrete pairs assert the same identities through
    their own test files):

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


class ReconstructionOperator(LinearOperatorMixin, ABC):
    r"""Abstract :math:`R : W \to V`. Sibling of :class:`ProjectionOperator`.

    A reconstruction operator lifts coefficients from a coarse space
    :math:`W` back into the fine space :math:`V`. The canonical
    Galerkin / Petrov-Galerkin pair :math:`(R, \Pi)` exposes both
    primitives: :class:`ProjectionOperator` carries :math:`\Pi`,
    :class:`ReconstructionOperator` carries :math:`R`.

    Concrete subclasses (e.g. :class:`HarmonicMomentReconstruction`)
    provide :meth:`apply`. The default capability set advertises
    :pydata:`CAP_APPLY`; subclasses MAY add :pydata:`CAP_APPLY_TRANSPOSE`
    if they can cheaply expose the representation transpose.

    Per Grand Report v3 §5.7 — the Operator hierarchy. Pair this with
    :class:`ProjectionOperator` (also §5.7) when shipping a Galerkin
    discretisation: the (R, Π) pair IS the discretisation, modulo a
    metric correction that lives on the operators' codomain spaces.
    """

    capabilities: frozenset[str] = frozenset({CAP_APPLY})

    @abstractmethod
    def apply(self, coefficients: np.ndarray) -> np.ndarray:
        ...


# ─────────────────────────────────────────────────────────────────────
# Concrete: harmonic-moment projection / reconstruction (real SH)
# ─────────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class MomentProjection(GalerkinProjection):
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

    .. note::

       The class name dropped the ``Harmonic`` prefix in P1.3 of the
       moment-space + layering plan; the typed codomain
       :class:`~orpheus.numerics.spaces.SphericalHarmonicSpace`
       (exposed via :attr:`codomain`) is now the carrier of
       "harmonic-ness" — the operator class itself is the generic
       :class:`MomentProjection`. The §10 PN moment-space projection
       (still future work) will be a sibling under the same name.

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
                f"MomentProjection: Y shape "
                f"{self.Y.shape} inconsistent with N={self.weights.shape[0]}, "
                f"L={self.L}; expected ({self.weights.shape[0]}, "
                f"{self.L + 1}, {2 * self.L + 1})."
            )

    @classmethod
    def from_measure(
        cls,
        measure: "DiscreteMeasure",
        L: int,
        *,
        Y: np.ndarray | None = None,
    ) -> "MomentProjection":
        r"""Construct from a :class:`DiscreteMeasure` over :math:`S^2`.

        The canonical factory.  Pulls the weights from the measure and
        either builds the spherical-harmonic table :math:`Y` from
        :class:`~orpheus.numerics.basis.SphericalHarmonicBasis` or
        accepts a precomputed ``Y`` (the optional kwarg) when the
        caller already holds it — e.g.,
        :class:`~orpheus.sn.scattering.ScatteringOperator.build_aniso_source`
        builds :math:`M` and :math:`R` from the same evaluator output.

        The measure's ``nodes`` MUST be a ``(N, 3)`` array of
        direction cosines :math:`(\mu_x, \mu_y, \mu_z)` per ordinate.

        Parameters
        ----------
        measure : DiscreteMeasure
            Angular quadrature on :math:`S^2`. ``measure.nodes`` is the
            ``(N, 3)`` direction-cosine array; ``measure.weights`` is
            the :math:`(N,)` quadrature weights.
        L : int
            Maximum harmonic order.
        Y : np.ndarray, shape ``(N, L+1, 2L+1)``, optional
            Precomputed spherical-harmonic table. If ``None`` (default),
            ``SphericalHarmonicBasis(L=L).evaluate(measure.nodes)`` is
            called internally.
        """
        nodes = measure.nodes
        if nodes.ndim != 2 or nodes.shape[1] != 3:
            raise ValueError(
                f"MomentProjection.from_measure expects "
                f"nodes of shape (N, 3); got {nodes.shape}."
            )
        if Y is None:
            from orpheus.numerics.basis.spherical_harmonic_basis import (
                SphericalHarmonicBasis,
            )
            Y = SphericalHarmonicBasis(L=L).evaluate(nodes)
        return cls(weights=measure.weights, Y=Y, L=L)

    @cached_property
    def codomain(self) -> "SphericalHarmonicSpace":
        r"""The :class:`~orpheus.numerics.spaces.SphericalHarmonicSpace`
        this projection targets — moment space of order :math:`L`.

        Built on access from :meth:`SphericalHarmonicSpace.from_L` so the
        metric carried by :attr:`FunctionSpace.inner_product_weights` is
        the canonical SH Gram :math:`g_C = \mathrm{diag}(4\pi/(2\ell+1))`
        padded to the moment-field storage layout ``(L+1, 2L+1)``. This
        is the metric that the generic ``A.H`` machinery uses (along with
        the domain's quadrature weights ``W``) to compute the W-weighted
        Hilbert adjoint :math:`\Pi^* = g_C \cdot S_0` correctly.

        ``range`` (the pre-P3.5 framework name read by
        :class:`~orpheus.numerics.operator._AdjointOperator`) aliases
        this property; both return the same :class:`SphericalHarmonicSpace`
        instance per call.
        """
        from orpheus.numerics.spaces.spherical_harmonic_space import (
            SphericalHarmonicSpace,
        )
        return SphericalHarmonicSpace.from_L(self.L)

    @cached_property
    def range(self) -> "SphericalHarmonicSpace":
        r"""Alias of :attr:`codomain` — the pre-P3.5 framework name.

        :class:`~orpheus.numerics.operator._AdjointOperator` reads
        ``range`` to find the codomain's inner-product weights for the
        Hilbert-adjoint calculation. P3.5 renames the framework attribute
        ``range`` → ``codomain``; until then both are defined here
        with identical return values (cached so the Krylov inner loop
        doesn't reallocate the :class:`SphericalHarmonicSpace` per
        matvec).
        """
        return self.codomain

    @cached_property
    def domain(self) -> FunctionSpace:
        r"""The angular-ordinate :class:`FunctionSpace` this projection consumes.

        Shape ``(N,)`` with the quadrature weights ``w_n`` as
        :attr:`FunctionSpace.inner_product_weights`. The 1-D shape is
        layout-agnostic; the
        :func:`~orpheus.numerics.operator._broadcast_for_leading_axes`
        helper pads the weights with trailing 1s to match higher-rank
        data tensors at adjoint-application time.

        This is the angular metric :math:`W = \mathrm{diag}(w_n)` that
        the generic ``A.H`` machinery reads to compute the W-weighted
        Hilbert adjoint correctly.
        """
        return FunctionSpace(
            name="angular_ordinate",
            shape=(self.weights.shape[0],),
            inner_product_weights=self.weights,
        )

    def apply(self, psi: np.ndarray) -> np.ndarray:
        r"""Apply :math:`\Pi \psi` along the leading angular axis.

        :math:`\psi` shape: ``(N, ...)`` — leading axis is the ordinate
        axis, trailing axes broadcast.
        Output shape: ``(L+1, 2L+1, ...)`` — leading two axes are the
        :math:`\ell, m` moment axes.
        """
        return np.einsum("n,nlm,n...->lm...", self.weights, self.Y, psi)

    def apply_transpose(self, c: np.ndarray) -> np.ndarray:
        r"""Apply the representation transpose :math:`\Pi^\top`.

        Since :math:`\Pi_{(\ell m), n} = w_n \, Y_\ell^m(\hat\Omega_n)`,
        the matrix transpose is :math:`\Pi^\top_{n, (\ell m)} = w_n \,
        Y_\ell^m(\hat\Omega_n)`, so

        .. math::
           :label: moment-projection-transpose

           (\Pi^\top c)_n
           \;=\; w_n \sum_{\ell, m} Y_\ell^m(\hat\Omega_n)\, c_\ell^m,

        the naked synthesis :math:`S_0(c)` post-multiplied by the
        quadrature weight :math:`w_n` on each ordinate. This is the
        **representation transpose** — the operation needed by the
        generic :class:`~orpheus.numerics.operator._AdjointOperator`
        machinery to compute the W-weighted Hilbert adjoint
        :math:`\Pi^* = g_C \cdot S_0` via

        .. math::

           (\Pi^* c)_n \;=\; \frac{1}{w_n}\, \Pi^\top(g_C \cdot c)_n
                       \;=\; \sum_{\ell, m} g_C^{\ell}\,
                              Y_\ell^m(\hat\Omega_n)\, c_\ell^m,

        with :math:`g_C^\ell = 4\pi/(2\ell+1)` carried on
        :attr:`codomain`'s ``inner_product_weights``.

        Numerical relationship (pinned by
        ``tests/numerics/test_spherical_harmonic_space.py``):

        .. math::

            \Pi \, R \;=\; 4\pi \cdot I \quad
              \text{(addition-theorem identity on band-limited inputs)}

            \Pi \, \Pi^* \;=\; \frac{4\pi}{2\ell+1} \cdot I_\ell \quad
              \text{(adjoint composition, diagonal in } \ell\text{)}.

        ERR-039 in one sentence: pre-P1.4 this method returned the bare
        :math:`S_0(c)` and was labeled the W-weighted Hilbert adjoint —
        but the true representation transpose carries :math:`w_n`, and
        the Hilbert adjoint additionally carries :math:`g_C`. The two
        are now separately typed: ``.apply_transpose`` is the
        representation transpose; ``.H`` is the Hilbert adjoint
        constructed by the metric-aware adjoint machinery.
        """
        return np.einsum("n,nlm,lm...->n...", self.weights, self.Y, c)


@dataclass(frozen=True)
class HarmonicMomentReconstruction(ReconstructionOperator):
    r"""Addition-theorem reconstruction :math:`R` paired with
    :class:`MomentProjection`.

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

    .. note::

       :math:`R` is **not** the strict Hilbert adjoint of
       :class:`MomentProjection` — see ERR-039. The addition-theorem
       reconstruction carries the :math:`(2\ell+1)` factor that the
       :math:`P_\ell` scattering reconstruction needs (Bell & Glasstone
       1970, §1.6); the strict W-weighted Hilbert adjoint
       :math:`\Pi^* = g_C \cdot S_0` lives on
       :attr:`MomentProjection.H` (Phase 1 P1.4) where it is computed
       generically by the metric-aware ``_AdjointOperator`` machinery
       using the codomain's :class:`~orpheus.numerics.spaces.SphericalHarmonicSpace`
       Gram. The two differ by exactly the
       :class:`~orpheus.numerics.spaces.SphericalHarmonicSpace`
       Gram :math:`g_C = \mathrm{diag}(4\pi/(2\ell+1))`:
       :math:`R = (2\ell+1) S_0 = 4\pi \cdot g_C^{-1} \cdot S_0`
       whereas :math:`\Pi^* = g_C \cdot S_0`. Both are useful
       primitives; the identity :math:`\Pi R = 4\pi I` is pinned by
       ``tests/numerics/test_spherical_harmonic_space.py``.

    Parameters
    ----------
    Y : np.ndarray, shape (N, L+1, 2L+1)
        Spherical-harmonic table from
        :meth:`~orpheus.numerics.basis.SphericalHarmonicBasis.evaluate`
        — same table consumed by the projection.
    two_l_plus_one : np.ndarray, shape (L+1,)
        The addition-theorem factor :math:`2\ell+1` precomputed for
        :math:`\ell = 0, 1, \ldots, L`. Use
        :meth:`from_spherical_harmonic_space` to source this from
        :attr:`SphericalHarmonicSpace.addition_theorem_factor` (the
        canonical single-source path); :meth:`from_Y` is the
        back-compat shim that internally constructs the space.

    Notes
    -----

    The capability set advertises :pydata:`CAP_APPLY` only. The
    canonical adjoint pairing flows through :class:`MomentProjection`
    via its ``.H`` property (which carries the W and :math:`g_C`
    metrics needed for the W-weighted Hilbert adjoint).

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
        r"""Construct from a :math:`Y` table.  The canonical factory.

        Sources the addition-theorem factor :math:`(2\ell+1)` from
        :attr:`SphericalHarmonicSpace.addition_theorem_factor` (which
        in turn delegates to
        :attr:`SphericalHarmonicBasis.addition_theorem_factor`).  The
        SH convention's :math:`(2\ell+1)` literal lives in exactly one
        place — :class:`~orpheus.numerics.basis.SphericalHarmonicBasis`
        — and the typed-space derivation here keeps the chain explicit.

        Parameters
        ----------
        Y : np.ndarray, shape ``(N, L+1, 2L+1)``
            Spherical-harmonic table built by
            :meth:`~orpheus.numerics.basis.SphericalHarmonicBasis.evaluate`
            on a quadrature's direction cosines.  :math:`L` is derived
            as ``Y.shape[1] - 1``.
        """
        from orpheus.numerics.spaces.spherical_harmonic_space import (
            SphericalHarmonicSpace,
        )
        L = Y.shape[1] - 1
        space = SphericalHarmonicSpace.from_L(L)
        return cls(Y=Y, two_l_plus_one=space.addition_theorem_factor)

    def apply(self, moments: np.ndarray) -> np.ndarray:
        r"""Apply :math:`R \phi` along the leading harmonic-coefficient axes.

        ``moments`` shape: ``(L+1, 2L+1, ...)`` — leading two axes are
        the :math:`\ell, m` moment axes.
        Output shape: ``(N, ...)`` — leading axis is the ordinate axis.
        """
        return np.einsum(
            "nlm,l,lm...->n...", self.Y, self.two_l_plus_one, moments,
        )
