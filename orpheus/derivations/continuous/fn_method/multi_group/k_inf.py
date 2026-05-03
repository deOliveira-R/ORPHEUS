r"""Branch-2 numpy implementations of the LA-13511 :math:`k_\infty` formulae.

These are direct translations of the closed forms verified
symbolically in :mod:`..origins.k_inf_derivations` — no eigenvalue
iteration, no matrix inversion beyond the formula's own algebra.
That structural simplicity is what makes them a credible
*structurally-independent* cross-check against
:func:`orpheus.derivations.common.eigenvalue.kinf_homogeneous`, which
solves the same problem via :func:`numpy.linalg.eig` of
:math:`A^{-1}F`.

ORPHEUS group convention is used throughout (``g=0`` fast → ``g=N-1``
slow). The 2G entry points take ORPHEUS-ordered cross sections; the
Sood-side renaming (``g=2`` ↔ ``g=0``, ``g=1`` ↔ ``g=1``) happens
inside the function so the user never deals with the reverse
convention.
"""
from __future__ import annotations

import numpy as np


def compute_kinf_1g(
    sigma_t: float,
    sigma_s: float,
    nu_sigma_f: float,
) -> float:
    r"""1G infinite-medium :math:`k_\infty` via Sood Eq 19.

    .. math::

        k_\infty = \frac{\nu\Sigma_f}{\Sigma_t - \Sigma_s}

    Parameters
    ----------
    sigma_t : float
        Total cross section.
    sigma_s : float
        Scattering cross section (total, all moments — for isotropic
        this is the only moment).
    nu_sigma_f : float
        Production cross section :math:`\nu\Sigma_f`.

    Returns
    -------
    float
        :math:`k_\infty`. Raises if :math:`\Sigma_t \le \Sigma_s`
        (subcritical or pure-absorber configuration).
    """
    Sigma_a = sigma_t - sigma_s
    if Sigma_a <= 0.0:
        raise ValueError(
            f"sigma_t ({sigma_t}) must exceed sigma_s ({sigma_s}); "
            f"absorption Sigma_a = {Sigma_a} ≤ 0."
        )
    return float(nu_sigma_f / Sigma_a)


def compute_kinf_2g_no_upscatter(
    sigma_t: np.ndarray,
    sigma_s: np.ndarray,
    nu_sigma_f: np.ndarray,
    chi: np.ndarray,
) -> float:
    r"""2G infinite-medium :math:`k_\infty` via Sood Eq 29 (no upscatter).

    Sood's printed Eq 29 uses **g=2 fast, g=1 slow** convention. This
    function takes ORPHEUS-ordered inputs (g=0 fast, g=1 slow) and
    relabels internally before applying:

    .. math::

        k_\infty = \frac{\chi_1\,\nu_1\Sigma_{1f}}{\Sigma_1^{\rm rem}}
                 + \chi_2\left[
                     \frac{\nu_1\Sigma_{1f}\,\Sigma_{12s}}
                          {\Sigma_1^{\rm rem}\,\Sigma_2^{\rm rem}}
                     + \frac{\nu_2\Sigma_{2f}}{\Sigma_2^{\rm rem}}
                   \right]

    where :math:`\Sigma_g^{\rm rem} = \Sigma_g - \Sigma_{ggs}` and the
    1↔fast, 2↔slow indices in the formula are SOOD-side. The check
    that the no-upscatter assumption holds is done at function entry.

    Parameters
    ----------
    sigma_t : (2,) array
        ORPHEUS-ordered total cross sections, ``[fast, slow]``.
    sigma_s : (2, 2) array
        ORPHEUS-ordered scattering matrix with ``[from, to]``
        convention. ``sigma_s[1, 0]`` MUST be 0 for the no-upscatter
        formula to apply.
    nu_sigma_f : (2,) array
        ORPHEUS-ordered production cross sections.
    chi : (2,) array
        ORPHEUS-ordered fission spectrum.

    Returns
    -------
    float
        :math:`k_\infty`.
    """
    sigma_t = np.asarray(sigma_t, dtype=float)
    sigma_s = np.asarray(sigma_s, dtype=float)
    nu_sigma_f = np.asarray(nu_sigma_f, dtype=float)
    chi = np.asarray(chi, dtype=float)
    if sigma_t.shape != (2,) or sigma_s.shape != (2, 2):
        raise ValueError(
            f"Expected 2-group inputs; got sigma_t.shape={sigma_t.shape}, "
            f"sigma_s.shape={sigma_s.shape}."
        )

    # ORPHEUS convention: sigma_s[from, to]. sigma_s[1, 0] is from
    # slow to fast = upscatter. Must be zero for Eq 29.
    if sigma_s[1, 0] > 0.0:
        raise ValueError(
            f"compute_kinf_2g_no_upscatter requires no upscatter "
            f"(sigma_s[1, 0] = 0); got {sigma_s[1, 0]}. Use "
            f"compute_kinf_2g_general for cases with upscatter."
        )

    # Map to Sood notation: 2 = fast = ORPHEUS index 0, 1 = slow = ORPHEUS index 1.
    Sigma_2 = sigma_t[0]    # fast total
    Sigma_1 = sigma_t[1]    # slow total
    Sigma_22s = sigma_s[0, 0]  # from fast to fast = self
    Sigma_11s = sigma_s[1, 1]  # from slow to slow = self
    Sigma_12s = sigma_s[0, 1]  # from fast to slow (Sood: from group 2 to group 1)
    nuSf_2 = nu_sigma_f[0]  # fast
    nuSf_1 = nu_sigma_f[1]  # slow
    chi_2 = chi[0]
    chi_1 = chi[1]

    Sigma_1rem = Sigma_1 - Sigma_11s
    Sigma_2rem = Sigma_2 - Sigma_22s

    if Sigma_1rem <= 0.0 or Sigma_2rem <= 0.0:
        raise ValueError(
            f"Sigma_g^rem must be positive; got Sigma_1rem={Sigma_1rem}, "
            f"Sigma_2rem={Sigma_2rem}."
        )

    k = (
        chi_1 * nuSf_1 / Sigma_1rem
        + chi_2 * (
            nuSf_1 * Sigma_12s / (Sigma_1rem * Sigma_2rem)
            + nuSf_2 / Sigma_2rem
        )
    )
    return float(k)


def compute_kinf_2g_general(
    sigma_t: np.ndarray,
    sigma_s: np.ndarray,
    nu_sigma_f: np.ndarray,
    chi: np.ndarray,
) -> float:
    r"""2G infinite-medium :math:`k_\infty` via the corrected Sood Eq 28.

    Sood Eq 28 as printed contains a typo (the :math:`\chi_g` numerator
    has the wrong :math:`\Sigma_g^{\rm rem}` factor; see
    :func:`..origins.k_inf_derivations.derive_kinf_2g_general_from_matrix`
    for the SymPy proof of this). This function uses the **corrected
    Eq 28** as derived from :math:`\det(M(k)) = 0` of Sood Eq 25:

    .. math::

        k_\infty = \frac{\chi_1\,(\nu_2\Sigma_{2f}\,\Sigma_{21s} + \Sigma_2^{\rm rem}\,\nu_1\Sigma_{1f})
                       + \chi_2\,(\nu_1\Sigma_{1f}\,\Sigma_{12s} + \Sigma_1^{\rm rem}\,\nu_2\Sigma_{2f})}
                      {\Sigma_1^{\rm rem}\,\Sigma_2^{\rm rem} - \Sigma_{12s}\,\Sigma_{21s}}

    For the no-upscatter case (:math:`\Sigma_{21s} = 0`) this should
    agree bit-for-bit with :func:`compute_kinf_2g_no_upscatter`; this
    is verified in the test suite.

    Parameters and return: same as :func:`compute_kinf_2g_no_upscatter`.
    """
    sigma_t = np.asarray(sigma_t, dtype=float)
    sigma_s = np.asarray(sigma_s, dtype=float)
    nu_sigma_f = np.asarray(nu_sigma_f, dtype=float)
    chi = np.asarray(chi, dtype=float)
    if sigma_t.shape != (2,) or sigma_s.shape != (2, 2):
        raise ValueError(
            f"Expected 2-group inputs; got sigma_t.shape={sigma_t.shape}, "
            f"sigma_s.shape={sigma_s.shape}."
        )

    # Sood notation mapping:
    Sigma_2 = sigma_t[0]
    Sigma_1 = sigma_t[1]
    Sigma_22s = sigma_s[0, 0]
    Sigma_11s = sigma_s[1, 1]
    Sigma_12s = sigma_s[0, 1]   # ORPHEUS [fast, slow] = Sood Σ_12s (from 2 to 1)
    Sigma_21s = sigma_s[1, 0]   # ORPHEUS [slow, fast] = Sood Σ_21s (from 1 to 2)
    nuSf_2 = nu_sigma_f[0]
    nuSf_1 = nu_sigma_f[1]
    chi_2 = chi[0]
    chi_1 = chi[1]

    Sigma_1rem = Sigma_1 - Sigma_11s
    Sigma_2rem = Sigma_2 - Sigma_22s

    if Sigma_1rem <= 0.0 or Sigma_2rem <= 0.0:
        raise ValueError(
            f"Sigma_g^rem must be positive; got Sigma_1rem={Sigma_1rem}, "
            f"Sigma_2rem={Sigma_2rem}."
        )

    denom = Sigma_1rem * Sigma_2rem - Sigma_12s * Sigma_21s
    if denom == 0.0:
        raise ValueError(
            "Denominator Sigma_1rem·Sigma_2rem - Sigma_12s·Sigma_21s = 0; "
            "the 2x2 system is singular."
        )

    numer = (
        chi_1 * (nuSf_2 * Sigma_21s + Sigma_2rem * nuSf_1)
        + chi_2 * (nuSf_1 * Sigma_12s + Sigma_1rem * nuSf_2)
    )
    return float(numer / denom)


def compute_flux_ratio_2g_no_upscatter(
    sigma_t: np.ndarray,
    sigma_s: np.ndarray,
    nu_sigma_f: np.ndarray,
    chi: np.ndarray,
    *,
    return_in_orpheus_order: bool = True,
) -> float:
    r"""2G flux ratio :math:`\phi_2/\phi_1` via Sood Eq 32.

    Sood publishes :math:`\phi_2/\phi_1` = (fast)/(slow). This function
    can return that ratio directly (Sood-style) or the inverse
    (ORPHEUS-style: :math:`\phi_{\rm slow}/\phi_{\rm fast}`).

    .. math::

        \frac{\phi_{\rm fast}}{\phi_{\rm slow}} =
            \frac{\Sigma_1^{\rm rem} - \nu_1\Sigma_{1f}/k_\infty}
                 {\nu_2\Sigma_{2f}/k_\infty - \Sigma_2^{\rm rem} + \Sigma_{12s}}

    where 1=slow, 2=fast in Sood's convention.

    Parameters
    ----------
    sigma_t, sigma_s, nu_sigma_f, chi : same as the k_inf functions.
    return_in_orpheus_order : bool
        If True (default), return :math:`\phi_{\rm slow}/\phi_{\rm fast}`
        — i.e. ``phi[g=1]/phi[g=0]`` consistent with how the
        kinf_and_spectrum_homogeneous solver returns the spectrum.
        If False, return Sood's published ratio :math:`\phi_2/\phi_1`
        = fast/slow.

    Returns
    -------
    float
        Flux ratio. Computed by:

        1. Solving Eq 29 (no upscatter) for :math:`k_\infty`.
        2. Substituting :math:`k_\infty` into Eq 32 for fast/slow.
        3. Inverting if ORPHEUS-order is requested.
    """
    k_inf = compute_kinf_2g_no_upscatter(sigma_t, sigma_s, nu_sigma_f, chi)

    Sigma_2 = sigma_t[0]
    Sigma_1 = sigma_t[1]
    Sigma_22s = sigma_s[0, 0]
    Sigma_11s = sigma_s[1, 1]
    Sigma_12s = sigma_s[0, 1]
    nuSf_2 = nu_sigma_f[0]
    nuSf_1 = nu_sigma_f[1]
    Sigma_1rem = Sigma_1 - Sigma_11s
    Sigma_2rem = Sigma_2 - Sigma_22s

    # Sood Eq 32: phi_fast/phi_slow = phi_2/phi_1
    numer = Sigma_1rem - nuSf_1 / k_inf
    denom = nuSf_2 / k_inf - Sigma_2rem + Sigma_12s
    fast_over_slow = float(numer / denom)

    if return_in_orpheus_order:
        return 1.0 / fast_over_slow
    return fast_over_slow


def compute_kinf_mg(
    sigma_t: np.ndarray,
    sigma_s: np.ndarray,
    nu_sigma_f: np.ndarray,
    chi: np.ndarray,
) -> float:
    r"""General multi-group :math:`k_\infty` via Sood Eq 76.

    .. math::

        k_\infty = \overline{\nu\Sigma_f}\,
                   (\overline{\overline{\Sigma_t}} -
                    \overline{\overline{\Sigma_s}})^{-1}\,
                   \bar\chi

    is the **single matrix inversion** form. Note: Sood Eq 76 uses
    the Sood scattering convention :math:`(\Sigma_s)_{gh}` = scattering
    FROM h TO g. ORPHEUS stores ``sigma_s[g, h]`` = scattering FROM
    g TO h (transpose). Therefore in this function we use
    ``sigma_s.T`` to match Sood's notation.

    Parameters
    ----------
    sigma_t : (G,) array
        Total cross sections in ORPHEUS group order.
    sigma_s : (G, G) array
        Scattering matrix in ORPHEUS ``[from, to]`` convention.
    nu_sigma_f : (G,) array
        Production cross sections.
    chi : (G,) array
        Fission spectrum (sums to 1).

    Returns
    -------
    float
        :math:`k_\infty`.
    """
    sigma_t = np.asarray(sigma_t, dtype=float)
    sigma_s = np.asarray(sigma_s, dtype=float)
    nu_sigma_f = np.asarray(nu_sigma_f, dtype=float)
    chi = np.asarray(chi, dtype=float)
    G = sigma_t.shape[0]
    if sigma_s.shape != (G, G):
        raise ValueError(
            f"sigma_s.shape={sigma_s.shape} inconsistent with sigma_t.shape={sigma_t.shape}."
        )

    # Sood Eq 76: A = Sigma_t - Sigma_s_sood, where Sigma_s_sood[g,h]
    # is FROM h TO g. ORPHEUS sigma_s is FROM g TO h, so transpose.
    A = np.diag(sigma_t) - sigma_s.T
    A_inv_chi = np.linalg.solve(A, chi)
    return float(nu_sigma_f @ A_inv_chi)
