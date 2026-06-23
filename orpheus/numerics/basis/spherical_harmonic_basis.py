r"""Real spherical-harmonic basis on :math:`S^2` truncated at degree :math:`L`.

Naming: for :math:`Y_\ell^m`, :math:`\ell` is the **degree** and :math:`m` the
**order** — so the truncation parameter :math:`L` is the maximum *degree*, not
the order. (The transport idiom "P_L scattering order" is a separate concept —
there "order" qualifies the Legendre expansion, not the SH azimuthal index.)

The canonical home of the :math:`Y_\ell^m` evaluator and of the Gram matrix
:math:`g_C = \mathrm{diag}(4\pi/(2\ell+1))` of ERR-039 fame. Pre-frame the
evaluator lived as a free ``evaluate_real_sh`` function and the Gram literal
:math:`(2\ell+1)` was carried as a raw array on the since-retired
harmonic-reconstruction operator. Both move here under a single typed class —
the SH convention, the addition-theorem factor, and the discrete Gram now have
one home (the frame's reconstruction face reads the factor live from
:attr:`addition_theorem_factor`).

Convention
==========

The project uses the **no-:math:`4\pi/(2\ell+1)`-prefactor** normalisation
of the real spherical harmonics, in which the addition theorem reads

.. math::
   :label: real-sh-addition-theorem

   \sum_{m=-\ell}^{\ell} Y_\ell^m(\hat\Omega)\,Y_\ell^m(\hat\Omega')
   \;=\; P_\ell(\hat\Omega \cdot \hat\Omega'),

so the :math:`P_\ell`-scattering reconstruction takes the form

.. math::

   q_n \;=\; \sum_{\ell=0}^{L} (2\ell+1) \sum_m Y_\ell^m(\hat\Omega_n)\,
            \phi^{\ell m},

and the discrete inner product against an exact quadrature satisfies

.. math::
   :label: sh-mass-matrix-diagonal

   \sum_n w_n \, Y_\ell^m(\hat\Omega_n) \, Y_{\ell'}^{m'}(\hat\Omega_n)
   \;=\; \frac{4\pi}{2\ell+1} \, \delta_{\ell\ell'} \delta_{m m'}.

The polar axis is :math:`\mu_x` (so :math:`\cos\theta = \mu_x`,
:math:`\sin\theta = \sqrt{1-\mu_x^2}`); azimuth is measured in the
:math:`(\mu_y,\mu_z)` plane:

.. math::

   \cos\phi \;=\; \frac{\mu_y}{\sin\theta}, \qquad
   \sin\phi \;=\; \frac{\mu_z}{\sin\theta}.

The :math:`\ell\le 1` branch is hard-coded to bit-identical values for
the legacy :math:`P_0/P_1` regression tests:

.. math::

   Y_0^0 = 1, \quad Y_1^{-1} = \mu_z, \quad Y_1^0 = \mu_x,
   \quad Y_1^{+1} = \mu_y.

For :math:`\ell \ge 2` the formula uses :func:`scipy.special.lpmv` with
the Condon–Shortley phase :math:`(-1)^m` removed and norm
:math:`\sqrt{2(\ell-m)!/(\ell+m)!}` for :math:`m \ne 0`.

Cross-method consumers
======================

This module is generic infrastructure: every method that integrates an
angular field against the spherical-harmonic basis consumes it.

* **SN scattering** (:mod:`orpheus.sn.scattering`) — the
  :math:`Y^* W` projection that builds the per-ordinate :math:`P_\ell`
  source.
* **PN solver** (future, §10 of the architecture report) —
  spherical-harmonic moment basis is the native space.
* **MC adjoint moments** — variance reduction with response moments
  built against :math:`Y_\ell^m`.
* **Energy-condensation diagnostics** — when within-group anisotropy
  needs an angular characterisation.

References
----------

* Bell, G. I. and Glasstone, S. (1970). *Nuclear Reactor Theory*.
  Van Nostrand Reinhold. §1.6 (real spherical harmonics in transport).
* Lewis, E. E. and Miller, W. F. Jr. (1993). *Computational Methods
  of Neutron Transport*. ANS. §4.7 (:math:`P_\ell` scattering Galerkin
  reconstruction with the :math:`(2\ell+1)` factor).
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from math import factorial, sqrt
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray
from scipy.special import lpmv

from orpheus.numerics.basis.base import Basis

if TYPE_CHECKING:
    from orpheus.numerics.measure import DiscreteMeasure
    from orpheus.numerics.spaces.spherical_harmonic_space import (
        SphericalHarmonicSpace,
    )


__all__ = ["SphericalHarmonicBasis"]


# ─────────────────────────────────────────────────────────────────────
# Basis class
# ─────────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class SphericalHarmonicBasis(Basis):
    r"""Real spherical harmonics on :math:`S^2`, truncated at degree :math:`L`.

    The first concrete :class:`~orpheus.numerics.basis.base.Basis` — the
    synthesis (trial) side of the spherical-harmonic frame. Implements the
    three fundamental basis operations (:meth:`evaluate`, :meth:`synthesize`,
    :meth:`mass_matrix`) and carries the SH convention (the no-prefactor
    normalisation, the addition-theorem factor :math:`2\ell+1`, the continuous
    Gram diagonal :math:`4\pi/(2\ell+1)`).

    Parameters
    ----------
    L : int
        Maximum harmonic degree retained. ``L == 0`` returns the :math:`P_0`
        table only. Negative ``L`` is rejected.

    Notes
    -----
    Frozen dataclass — equality and hashing are by ``L`` alone (the basis
    IS its truncation degree; two ``SphericalHarmonicBasis(L=2)`` instances
    are the same basis).

    The basis carries the convention (the no-prefactor normalisation
    documented in the module docstring) and the addition-theorem factor
    :math:`2\ell+1`. The continuous (theoretical) mass-matrix diagonal
    :math:`4\pi/(2\ell+1)` is the :meth:`metric_per_ell` property; the
    DISCRETE mass matrix against a quadrature is computed by
    :meth:`discrete_mass_matrix`. The two agree to within the
    quadrature's degree of exactness.
    """

    L: int

    def __post_init__(self) -> None:
        if self.L < 0:
            raise ValueError(
                f"SphericalHarmonicBasis: L must be non-negative, got L={self.L}"
            )

    # ── Convention-bearing properties ────────────────────────────────

    @cached_property
    def addition_theorem_factor(self) -> NDArray:
        r"""The :math:`(2\ell+1)` array, shape ``(L+1,)``.

        Used by the addition-theorem reconstruction
        :math:`R = (2\ell+1) Y` and equal to
        :math:`4\pi \cdot g_C^{-1}` where :math:`g_C` is the SH Gram matrix.
        """
        return 2.0 * np.arange(self.L + 1) + 1.0

    @cached_property
    def metric_per_ell(self) -> NDArray:
        r"""The Gram-matrix diagonal :math:`4\pi/(2\ell+1)` per :math:`\ell`, shape ``(L+1,)``.

        This is the THEORETICAL (continuous-:math:`L^2`-on-:math:`S^2`)
        metric under the project's no-prefactor SH convention. The
        DISCRETE counterpart against a quadrature is the diagonal of
        :meth:`discrete_mass_matrix`; the two agree iff the quadrature
        is exact for :math:`Y_\ell^m Y_{\ell'}^{m'}` of degree
        :math:`\ell + \ell' \le 2L`.
        """
        return 4.0 * np.pi / self.addition_theorem_factor

    # ── Tabulation ────────────────────────────────────────────────────

    def evaluate(self, directions: NDArray) -> NDArray:
        r"""Tabulate :math:`Y_\ell^m(\hat\Omega_n)` at the given direction set.

        Parameters
        ----------
        directions : NDArray, shape ``(N, 3)``
            Direction cosines :math:`(\mu_x, \mu_y, \mu_z)` per ordinate.

        Returns
        -------
        NDArray, shape ``(N, L+1, 2L+1)``
            ``Y[n, l, l+m]`` is :math:`Y_\ell^m(\hat\Omega_n)` under the
            no-prefactor convention; entries outside :math:`|m| \le \ell`
            are zero.
        """
        directions = np.asarray(directions)
        if directions.ndim != 2 or directions.shape[1] != 3:
            raise ValueError(
                f"SphericalHarmonicBasis.evaluate expects directions of shape "
                f"(N, 3); got {directions.shape}."
            )
        return _evaluate_real_sh(
            self.L, directions[:, 0], directions[:, 1], directions[:, 2],
        )

    def evaluate_from_components(
        self,
        mu_x: NDArray,
        mu_y: NDArray,
        mu_z: NDArray,
    ) -> NDArray:
        r"""Tabulate from separated direction-cosine arrays.

        Equivalent to :meth:`evaluate` but accepts three :math:`(N,)`
        arrays instead of one :math:`(N, 3)` array. Provided for
        :class:`~orpheus.numerics.quadrature.Quadrature` and other
        per-component consumers (the SN quadratures historically expose
        ``mu_x`` / ``mu_y`` / ``mu_z`` as separate attributes).

        Parameters
        ----------
        mu_x, mu_y, mu_z : NDArray, shape ``(N,)``
            Direction-cosine arrays.

        Returns
        -------
        NDArray, shape ``(N, L+1, 2L+1)``
            Same layout as :meth:`evaluate`.
        """
        return _evaluate_real_sh(self.L, mu_x, mu_y, mu_z)

    # ── Mass matrix ───────────────────────────────────────────────────

    def mass_matrix(self, measure: "DiscreteMeasure") -> NDArray:
        r"""Discrete Gram matrix :math:`\sum_n w_n Y_\ell^m Y_{\ell'}^{m'}` over a quadrature.

        Computes the :math:`(L+1, 2L+1, L+1, 2L+1)` 4-tensor

        .. math::

            g[\ell, m, \ell', m']
            \;=\; \sum_n w_n \, Y_\ell^m(\hat\Omega_n) \, Y_{\ell'}^{m'}(\hat\Omega_n).

        For an exact quadrature of degree :math:`\ge 2L` this equals

        .. math::

            \mathrm{diag}\!\left(\frac{4\pi}{2\ell+1}\right)
            \delta_{\ell\ell'} \delta_{m m'}

        per :eq:`sh-mass-matrix-diagonal`, agreeing with
        :attr:`metric_per_ell` along the diagonal and vanishing
        off-diagonal. The match-to-:attr:`metric_per_ell` is the
        ERR-039 test gate pinned by
        ``tests/numerics/test_spherical_harmonic_space.py``.

        Parameters
        ----------
        measure : DiscreteMeasure
            Angular quadrature on :math:`S^2`. ``measure.nodes`` MUST be
            a ``(N, 3)`` array of direction cosines; ``measure.weights``
            the :math:`(N,)` quadrature weights.

        Returns
        -------
        NDArray, shape ``(L+1, 2L+1, L+1, 2L+1)``
            Full 4-tensor Gram matrix. Off-diagonal entries that vanish
            under an exact quadrature are bit-zero only at exact arithmetic;
            small FP residuals are expected and used by tests as a
            quadrature-exactness diagnostic.
        """
        Y = self.evaluate(measure.nodes)
        return np.einsum("n,nlm,nLM->lmLM", measure.weights, Y, Y)

    # ── Table contractions (the Frame caches the table, delegates here) ───

    def synthesize(self, coefficients: NDArray, table: NDArray) -> NDArray:
        r"""Naked synthesis :math:`S_0(c)_n = \sum_{\ell, m} Y_\ell^m(\hat\Omega_n)\, c_\ell^m`.

        The bare reconstruction with NO :math:`(2\ell+1)` factor and NO
        :math:`w_n` weight — the pure synthesis :math:`S_0` (the frame-theory
        synthesis operator :math:`T^*`). The three weighted contractions
        (:meth:`analyze`, :meth:`analyze_transpose`, :meth:`reconstruct`) are each
        this kernel post-multiplied by ONE diagonal, but kept as separate fused
        einsums for 0-ULP bit-identity (FP non-associativity).

        Parameters
        ----------
        coefficients : NDArray, shape ``(L+1, 2L+1, ...)``
            Moment-space input; entries outside :math:`|m| \le \ell` are zero by
            construction (dotted with zero ``table`` entries).
        table : NDArray, shape ``(N, L+1, 2L+1)``
            The :math:`Y_\ell^m(\hat\Omega_n)` table from :meth:`evaluate`.

        Returns
        -------
        NDArray, shape ``(N, ...)`` — :math:`S_0(c)` per ordinate.
        """
        return np.einsum("nlm,lm...->n...", table, coefficients)

    def analyze(
        self, values: NDArray, table: NDArray, weights: NDArray,
    ) -> NDArray:
        r"""Analysis :math:`(M\psi)_\ell^m = \sum_n w_n\, Y_\ell^m(\hat\Omega_n)\,\psi_n`.

        The W-weighted Galerkin projection of a per-ordinate field onto harmonic
        moments — the frame's *analysis operator* :math:`T`. ONE fused
        :func:`numpy.einsum` over the ordinate axis; trailing axes broadcast.
        """
        return np.einsum("n,nlm,n...->lm...", weights, table, values)

    def analyze_transpose(
        self, coefficients: NDArray, table: NDArray, weights: NDArray,
    ) -> NDArray:
        r"""Representation transpose :math:`(M^\top c)_n = w_n \sum_{\ell m} Y_\ell^m(\hat\Omega_n)\, c_\ell^m`.

        The matrix transpose of :meth:`analyze` (:math:`= w_n \cdot S_0`) — NOT the
        Hilbert adjoint. The metric-aware ``_AdjointOperator`` combines it with the
        domain/codomain Gram to give :math:`M^* = g_C \cdot S_0`, so the Frame's
        analysis face gets ``.H`` for free.
        """
        return np.einsum("n,nlm,lm...->n...", weights, table, coefficients)

    def reconstruct(self, coefficients: NDArray, table: NDArray) -> NDArray:
        r"""Reconstruction :math:`(R\phi)_n = \sum_\ell (2\ell+1) \sum_m Y_\ell^m(\hat\Omega_n)\, \phi_\ell^m`.

        The addition-theorem (canonical-dual) synthesis — :math:`S_0` weighted by
        the dual factor :math:`2\ell+1 = 4\pi\, g_C^{-1}`, read live from
        :attr:`addition_theorem_factor` (no stored copy). Measure-free.
        """
        return np.einsum(
            "nlm,l,lm...->n...", table, self.addition_theorem_factor, coefficients,
        )

    def reconstruct_transpose(self, values: NDArray, table: NDArray) -> NDArray:
        r"""Representation transpose :math:`(R^\top v)_\ell^m = (2\ell+1) \sum_n Y_\ell^m(\hat\Omega_n)\, v_n`.

        The matrix transpose of :meth:`reconstruct` (:math:`= 2\ell+1` times the
        naked analysis :math:`S_0^\top`) — NOT the Hilbert adjoint, and **measure-free**:
        no :math:`w_n` is baked in (symmetric with :meth:`reconstruct`, asymmetric
        with :meth:`analyze_transpose`, whose forward bakes the weights in). The
        metric-aware ``_AdjointOperator`` combines it with the codomain (measure
        :math:`w_n`) and domain (basis :math:`g_C`) Gram to give the W-weighted
        Hilbert adjoint :math:`(R^* v)_\ell^m = \frac{(2\ell+1)^2}{4\pi} \sum_n w_n
        Y_\ell^m(\hat\Omega_n)\, v_n`, so the Frame's reconstruction face gets ``.H``
        for free.
        """
        return np.einsum(
            "nlm,l,n...->lm...", table, self.addition_theorem_factor, values,
        )

    # ── Coefficient space ─────────────────────────────────────────────────

    @property
    def space(self) -> "SphericalHarmonicSpace":
        r"""The :class:`SphericalHarmonicSpace` of degree :math:`L` this basis spans.

        Its Gram :math:`g_C = \mathrm{diag}(4\pi/(2\ell+1))` (carried as
        ``inner_product_weights``) is the codomain metric the
        :class:`~orpheus.numerics.frame.Frame` and the Hilbert-adjoint machinery
        read. Lazy import: ``SphericalHarmonicSpace`` imports this basis, so a
        top-level import would cycle.
        """
        from orpheus.numerics.spaces.spherical_harmonic_space import (
            SphericalHarmonicSpace,
        )
        return SphericalHarmonicSpace.from_L(self.L)


# ─────────────────────────────────────────────────────────────────────
# Algorithm body (private — preserved bit-identical from
# legacy spherical_harmonics.evaluate_real_sh)
# ─────────────────────────────────────────────────────────────────────


def _evaluate_real_sh(
    L: int,
    mu_x: NDArray,
    mu_y: NDArray,
    mu_z: NDArray,
) -> NDArray:
    r"""Implementation body of :meth:`SphericalHarmonicBasis.evaluate`.

    Kept as a free function so the back-compat shim at
    :mod:`orpheus.numerics.spherical_harmonics` can re-export this
    without instantiating a class. Algorithm preserved bit-identical to
    the legacy ``evaluate_real_sh`` so the snapshots at
    ``tests/sn/regression/snapshots/`` continue to pass.

    See the module docstring for the convention and citations.
    """
    mu_x = np.asarray(mu_x)
    mu_y = np.asarray(mu_y)
    mu_z = np.asarray(mu_z)
    N = len(mu_x)
    if L < 0:
        return np.zeros((N, 0, 0))
    Y = np.zeros((N, L + 1, 2 * L + 1))

    Y[:, 0, 0] = 1.0
    if L == 0:
        return Y

    Y[:, 1, 0] = mu_z   # m = -1
    Y[:, 1, 1] = mu_x   # m =  0
    Y[:, 1, 2] = mu_y   # m = +1
    if L == 1:
        return Y

    cos_theta = mu_x
    sin_theta = np.sqrt(np.maximum(1.0 - cos_theta * cos_theta, 0.0))
    on_axis = sin_theta < 1e-15
    safe_st = np.where(on_axis, 1.0, sin_theta)
    cos_phi = np.where(on_axis, 1.0, mu_y / safe_st)
    sin_phi = np.where(on_axis, 0.0, mu_z / safe_st)
    phi = np.arctan2(sin_phi, cos_phi)

    for l in range(2, L + 1):
        Y[:, l, l] = lpmv(0, l, cos_theta)  # m = 0: P_l(μ_x)
        for m in range(1, l + 1):
            P_lm = lpmv(m, l, cos_theta)
            sign = (-1.0) ** m   # remove Condon–Shortley phase
            norm = sqrt(2.0 * factorial(l - m) / factorial(l + m))
            cos_mphi = np.cos(m * phi)
            sin_mphi = np.sin(m * phi)
            Y[:, l, l + m] = sign * norm * P_lm * cos_mphi   # m > 0
            Y[:, l, l - m] = sign * norm * P_lm * sin_mphi   # m < 0
    return Y
