r"""Bare-critical multiplying slab eigenvalue spectrum via Dahl-Sjostrand 1979.

Solves the Galerkin spectral expansion of Carlvik's integral equation
:cite:`DahlSjostrand1979` Eq. (1), in the slab case (q=0, even-Legendre
basis :math:`P_0, P_2, \ldots, P_{2(N-1)}`), via the block-matrix
linearization Eq. (4). Returns the **full eigenvalue spectrum**
(typically 2N values, including complex-conjugate pairs at high
:math:`\bar\mu`), with the dominant real eigenvalue being the
critical secondaries-ratio :math:`c_{\rm crit}`.

Reproduces Dahl-Sjostrand 1979 Table II to 7 significant figures
across :math:`\bar\mu \in \{0, 0.10, 0.15, 0.20, 0.25, 0.30\}` and
:math:`d \in \{0.2, 1, 2, 5, 8, 20\}` mfp.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..core.galerkin_matrix import solve_eq4_eigenproblem


@dataclass(frozen=True)
class CarlvikGalerkinSlabResult:
    r"""Result of :func:`solve_galerkin_spectral_slab`.

    Attributes
    ----------
    c_critical : float
        Fundamental (dominant) real eigenvalue of Eq. (4); the
        secondaries-per-collision ratio at criticality
        :math:`c_{\rm crit}(\bar\mu, d)`.
    eigenvalue_spectrum : np.ndarray
        Length-:math:`2N` complex array of all eigenvalues
        :math:`c_j`, sorted by ascending :math:`\Re c_j`. Complex
        eigenvalues appear as conjugate pairs (only at high
        :math:`\bar\mu \ge 0.15` for higher modes per Dahl-Sjostrand
        Figs. 3-6).
    eigenvectors : np.ndarray
        :math:`(N, 2N)` matrix of upper-half eigenvectors :math:`F_j`
        (one per column). Each is an :math:`N`-dim vector of Legendre
        coefficients :math:`F_n` (where the basis is :math:`P_0,
        P_2, \ldots, P_{2(N-1)}`). The lower-half eigenvector
        :math:`K_j = c_j F_j` is redundant (Dahl-Sjostrand p. 119).
    n_modes : int
        Galerkin expansion order :math:`N`.
    a : float
        Half-thickness in mfp.
    d : float
        Total slab thickness in mfp = :math:`2a`.
    mu_bar : float
        Linearly-anisotropic mean cosine.
    n_quad : int
        Outer Gauss-Legendre quadrature order used.
    indices : np.ndarray
        :math:`P_n` basis indices (slab: :math:`[0, 2, 4, \ldots,
        2(N-1)]`).
    """

    c_critical: float
    eigenvalue_spectrum: np.ndarray
    eigenvectors: np.ndarray
    n_modes: int
    a: float
    d: float
    mu_bar: float
    n_quad: int
    indices: np.ndarray


def solve_galerkin_spectral_slab(
    c: float,
    d: float,
    mu_bar: float = 0.0,
    *,
    n_modes: int = 9,
    n_quad: int = 128,
    return_full_spectrum: bool = True,
) -> CarlvikGalerkinSlabResult:
    r"""Anisotropic-scattering slab eigenvalue spectrum.

    .. note::

       Per Dahl-Sjostrand 1979, the eigenvalue formulation is
       **material-agnostic**: the eigenvalues :math:`c_j(\bar\mu, d)`
       depend only on geometry parameters :math:`(\bar\mu, d)`, not
       on cross sections. The ``c`` argument here is the user's
       **claimed** value; the function returns the *actual* critical
       :math:`c_{\rm crit}(d, \bar\mu)` from the eigenvalue equation
       and the user can compare to determine criticality.

       For a critical slab (:math:`k_{\rm eff} = 1`), :math:`c_{\rm crit}`
       equals :math:`c = (\Sigma_s + \nu \Sigma_f)/\Sigma_t`. For a
       sub/supercritical slab, the implied :math:`k_{\rm eff} =
       c_{\rm material} / c_{\rm crit}(\bar\mu, d)`.

    Parameters
    ----------
    c : float
        User's claimed mean number of secondaries per collision.
        Stored on the result for diagnostic purposes; does NOT
        affect the eigenvalue computation.
    d : float
        Total slab thickness in mfp = :math:`2a`. Must be > 0.
    mu_bar : float, default 0.0
        Linearly-anisotropic mean scattering cosine
        :math:`\bar\mu = \int \mu f(\mu) d\mu`. Range:
        :math:`[-1, +1]`; physical isotropic = 0; forward-peaked
        = positive; back-peaked = negative.
    n_modes : int, default 9
        Galerkin expansion order :math:`N`. Dahl-Sjostrand
        recommends 9 for saturation against the published 7-sig-fig
        tables.
    n_quad : int, default 128
        Gauss-Legendre quadrature order for matrix-element
        evaluation. Default 128 reproduces Dahl-Sjostrand Table II
        to 7 sig figs across all tabulated :math:`(\bar\mu, d)`.
    return_full_spectrum : bool, default True
        If ``True``, return the full :math:`2N`-eigenvalue spectrum
        with corresponding eigenvectors. If ``False``, return only
        :math:`c_{\rm critical}` and the dominant eigenvector
        (sets the other fields to length-1 stubs for memory
        efficiency).

    Returns
    -------
    CarlvikGalerkinSlabResult
    """
    if d <= 0.0:
        raise ValueError(f"d must be positive, got d={d}")
    if n_modes < 1:
        raise ValueError(f"n_modes must be >= 1, got {n_modes}")
    if abs(mu_bar) > 1.0:
        raise ValueError(f"mu_bar must be in [-1, +1], got {mu_bar}")

    a = d / 2.0
    indices = np.arange(0, 2 * n_modes, 2, dtype=int)

    c_vals, F_eigvecs, _ = solve_eq4_eigenproblem(
        a=a,
        n_modes=n_modes,
        mu_bar=mu_bar,
        geometry="slab",
        n_quad=n_quad,
    )

    # Dominant eigenvalue: smallest real positive c. Per
    # Dahl-Sjostrand p. 124: "the eigenvalue of the fundamental
    # mode is real over the whole range of μ̄ values." So we filter
    # for real eigenvalues, and take the smallest among them.
    real_mask = np.abs(c_vals.imag) < 1e-9
    real_c_vals = c_vals[real_mask].real
    real_c_positive = real_c_vals[real_c_vals > 0.0]
    if len(real_c_positive) == 0:
        raise RuntimeError(
            f"No positive real eigenvalue found at d={d}, mu_bar={mu_bar}, "
            f"n_modes={n_modes}. Spectrum: {c_vals}"
        )
    c_critical = float(np.min(real_c_positive))

    if return_full_spectrum:
        return CarlvikGalerkinSlabResult(
            c_critical=c_critical,
            eigenvalue_spectrum=c_vals,
            eigenvectors=F_eigvecs,
            n_modes=n_modes,
            a=a,
            d=d,
            mu_bar=mu_bar,
            n_quad=n_quad,
            indices=indices,
        )
    else:
        # Find the dominant-real-eigenvalue index in the sorted spectrum.
        idx = int(np.argmin(np.where(real_mask & (c_vals.real > 0), c_vals.real, np.inf)))
        return CarlvikGalerkinSlabResult(
            c_critical=c_critical,
            eigenvalue_spectrum=np.array([c_vals[idx]]),
            eigenvectors=F_eigvecs[:, idx:idx+1],
            n_modes=n_modes,
            a=a,
            d=d,
            mu_bar=mu_bar,
            n_quad=n_quad,
            indices=indices,
        )
