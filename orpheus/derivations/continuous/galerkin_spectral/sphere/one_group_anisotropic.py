r"""Bare-critical multiplying sphere eigenvalue spectrum via Dahl-Sjostrand 1979.

Solves the Galerkin spectral expansion of Carlvik's integral equation
:cite:`DahlSjostrand1979` Eq. (1) in the sphere case (q=1, odd-Legendre
basis :math:`P_1, P_3, \ldots, P_{2N-1}` for the reduced flux
:math:`r \phi(r)`), via Eq. (4) block-matrix linearization.

Reproduces Dahl-Sjostrand 1979 Table I to 7 significant figures across
:math:`\bar\mu \in \{0, 0.10, 0.15, 0.20, 0.25, 0.30\}` and
:math:`d \in \{0.2, 1, 2, 5, 8, 20\}` mfp (where :math:`d = 2 R` is
the diameter).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..core.galerkin_matrix import solve_eq4_eigenproblem


@dataclass(frozen=True)
class CarlvikGalerkinSphereResult:
    r"""Result of :func:`solve_galerkin_spectral_sphere`.

    Attributes
    ----------
    c_critical : float
        Fundamental real eigenvalue (smallest positive real
        :math:`c_j`). For a critical sphere :math:`c_{\rm crit}` =
        the secondaries-per-collision ratio at criticality.
    eigenvalue_spectrum : np.ndarray
        Length-:math:`2N` complex array of all eigenvalues, sorted
        by ascending :math:`\Re c`.
    eigenvectors : np.ndarray
        :math:`(N, 2N)` upper-half eigenvectors :math:`F_j` (one
        per column). The basis is :math:`P_1, P_3, \ldots, P_{2N-1}`,
        so :math:`F_j` is the expansion of :math:`r \phi(r)`.
    n_modes : int
        Number of Legendre terms used.
    R : float
        Sphere radius in mfp = :math:`a`. (Note: Dahl-Sjostrand uses
        :math:`a` for radius and :math:`d = 2 a` for "diameter".)
    d : float
        Sphere diameter in mfp = :math:`2 R`.
    mu_bar : float
        Linearly-anisotropic mean cosine.
    n_quad : int
        Quadrature order.
    indices : np.ndarray
        :math:`P_n` basis indices (sphere: :math:`[1, 3, 5, \ldots,
        2N-1]`).
    """

    c_critical: float
    eigenvalue_spectrum: np.ndarray
    eigenvectors: np.ndarray
    n_modes: int
    R: float
    d: float
    mu_bar: float
    n_quad: int
    indices: np.ndarray


def solve_galerkin_spectral_sphere(
    c: float,
    d: float,
    mu_bar: float = 0.0,
    *,
    n_modes: int = 9,
    n_quad: int = 128,
    return_full_spectrum: bool = True,
) -> CarlvikGalerkinSphereResult:
    r"""Anisotropic-scattering sphere eigenvalue spectrum.

    See :func:`...slab.solve_galerkin_spectral_slab` for the parameter
    semantics — they are identical between slab and sphere; only the
    Legendre basis (odd vs even) and the boundary-chord parameter
    :math:`q` differ.

    Parameters
    ----------
    c : float
        Diagnostic only.
    d : float
        Sphere diameter in mfp = :math:`2 R`. Must be > 0.
    mu_bar : float, default 0.0
        Mean scattering cosine.
    n_modes : int, default 9
        Galerkin expansion order.
    n_quad : int, default 128
        Quadrature order.
    return_full_spectrum : bool, default True
        If False, return only the fundamental.

    Returns
    -------
    CarlvikGalerkinSphereResult
    """
    if d <= 0.0:
        raise ValueError(f"d must be positive, got d={d}")
    if n_modes < 1:
        raise ValueError(f"n_modes must be >= 1, got {n_modes}")
    if abs(mu_bar) > 1.0:
        raise ValueError(f"mu_bar must be in [-1, +1], got {mu_bar}")

    R = d / 2.0  # radius
    indices = np.arange(1, 2 * n_modes + 1, 2, dtype=int)

    c_vals, F_eigvecs, _ = solve_eq4_eigenproblem(
        a=R,
        n_modes=n_modes,
        mu_bar=mu_bar,
        geometry="sphere",
        n_quad=n_quad,
    )

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
        return CarlvikGalerkinSphereResult(
            c_critical=c_critical,
            eigenvalue_spectrum=c_vals,
            eigenvectors=F_eigvecs,
            n_modes=n_modes,
            R=R,
            d=d,
            mu_bar=mu_bar,
            n_quad=n_quad,
            indices=indices,
        )
    else:
        idx = int(np.argmin(np.where(real_mask & (c_vals.real > 0), c_vals.real, np.inf)))
        return CarlvikGalerkinSphereResult(
            c_critical=c_critical,
            eigenvalue_spectrum=np.array([c_vals[idx]]),
            eigenvectors=F_eigvecs[:, idx:idx+1],
            n_modes=n_modes,
            R=R,
            d=d,
            mu_bar=mu_bar,
            n_quad=n_quad,
            indices=indices,
        )
