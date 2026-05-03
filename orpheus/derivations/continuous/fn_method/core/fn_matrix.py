r"""F_N collocation matrix assembler — slab and sphere.

The F_N method's load-bearing object is the :math:`(N+1) \times (N+1)`
collocation matrix :math:`M(R)`. For 1G bare-critical problems the
matrix entry has the unified form

.. math::

   M_{\beta,\alpha} = B_\alpha(\xi_\beta) + s\,e^{-2R/\xi_\beta}\,
   A_\alpha(\xi_\beta), \qquad
   \beta, \alpha = 0, 1, \ldots, N

where

* :math:`s \in \{+1, -1\}` is the **geometry sign**: :math:`s = +1`
  for the slab (Siewert-Benoist 1979 Eq. 4 BC, symmetric flux),
  :math:`s = -1` for the sphere (Siewert-Thomas 1986 Eq. 46 BC,
  anti-symmetric flux);
* :math:`B_\alpha`, :math:`A_\alpha` are the half-range moment
  integrals defined in :mod:`.moments`;
* :math:`\xi_\beta` are the Grandjean-Siewert collocation points
  defined in :func:`.moments.collocation_points`;
* :math:`R` is the critical dimension (slab half-thickness or
  sphere radius, in mean-free paths).

The structural insight (Siewert-Thomas 1986 p. 268, verbatim quote):

    "the critical sphere problem requires only that Eq. (4) be
    changed to read [Eq. 46], and that we interpret a as the
    critical radius, we have incorporated the relevant minus sign
    in our developed equations"

— means slab and sphere F_N differ by **exactly one sign on the
attenuation block**. All other machinery (moment recursion,
collocation grid, dispersion-root finder, X-function definition
since X depends only on :math:`\Lambda(z) = 1 - (cz/2)\,
\mathrm{atanh}(1/z)` which is a medium property) is shared.

This module factors that shared machinery out so :mod:`..slab.one_group`
and :mod:`..sphere.one_group` reuse the same assembler with different
``geometry_sign``.

References
----------

* Siewert & Benoist 1979, *Nucl. Sci. Eng.* **69**, 156, Eq. 49.
* Grandjean & Siewert 1979, *Nucl. Sci. Eng.* **69**, 161, Eq. 25.
* Siewert & Thomas 1986, *Nucl. Sci. Eng.* **94**, 264, Eq. 46.
"""
from __future__ import annotations

import cmath
from typing import Literal

import numpy as np


GeometrySign = Literal[+1, -1]


def _moment_arrays_complex(
    N: int, xi: complex, c: float
) -> tuple[np.ndarray, np.ndarray]:
    r"""Compute :math:`(B_\alpha(\xi))` and :math:`(A_\alpha(\xi))` for
    :math:`\alpha = 0, \ldots, N` at a complex :math:`\xi`.

    Special-cases :math:`\xi = 0` to avoid division-by-zero in the
    :math:`\xi \log(1 + 1/\xi)` term. The :math:`\xi = 0` limits are:

    * :math:`B_0(0) = 2/c - 1` (the :math:`\xi \log` term vanishes).
    * :math:`B_\alpha(0) = -1/(\alpha+1)` for :math:`\alpha \ge 1`.
    * :math:`A_0(0) = 1`.
    * :math:`A_\alpha(0) = 1/(\alpha+1)` for :math:`\alpha \ge 1`.
    """
    Bvec = np.empty(N + 1, dtype=complex)
    Avec = np.empty(N + 1, dtype=complex)
    if xi == 0.0:
        Bvec[0] = (2.0 / c) - 1.0
        Avec[0] = 1.0
    else:
        log_term = xi * cmath.log(1.0 + 1.0 / xi)
        Bvec[0] = (2.0 / c) - 1.0 - log_term
        Avec[0] = 1.0 - log_term
    for n in range(1, N + 1):
        Bvec[n] = xi * Bvec[n - 1] - 1.0 / (n + 1)
        Avec[n] = -xi * Avec[n - 1] + 1.0 / (n + 1)
    return Bvec, Avec


def assemble_fn_matrix(
    c: float,
    R_mfp: float,
    geometry_sign: int,
    n_modes: int,
    xis: np.ndarray,
) -> np.ndarray:
    r"""Assemble the F_N matrix :math:`M(R)` for slab or sphere.

    Row :math:`\beta` (collocation point :math:`\xi_\beta`) has entries

    .. math::

       M_{\beta,\alpha} = B_\alpha(\xi_\beta) + s\,e^{-2R/\xi_\beta}\,
       A_\alpha(\xi_\beta), \qquad \alpha = 0, \ldots, N

    where :math:`s \in \{+1, -1\}` is ``geometry_sign``.

    The :math:`\xi = 0` row uses :math:`e^{-\infty} = 0` (since
    :math:`R > 0`), leaving only the :math:`B_\alpha(0)` part — and
    therefore the geometry_sign does not affect the :math:`\xi = 0`
    row. This is consistent with both slab and sphere
    formulations: the :math:`\xi = 0` collocation contributes only
    the moment-integral term in both cases.

    Parameters
    ----------
    c : float
        Mean number of secondaries per collision, :math:`c = (\Sigma_s
        + \nu\Sigma_f)/\Sigma_t`. Determines the :math:`B_0` seed.
    R_mfp : float
        Critical dimension in mean-free paths (slab half-thickness or
        sphere radius, depending on ``geometry_sign``).
    geometry_sign : int
        :math:`+1` for slab (Siewert-Benoist Eq. 4 symmetric BC),
        :math:`-1` for sphere (Siewert-Thomas Eq. 46 anti-symmetric BC).
    n_modes : int
        F_N order :math:`N`. The matrix is :math:`(N+1) \times (N+1)`.
    xis : np.ndarray
        Length-:math:`(N+1)` complex array of collocation points.
        Typically: :math:`\xi_0 = \nu_0` (Case discrete eigenvalue,
        possibly imaginary for :math:`c > 1`), :math:`\xi_1 = 0`,
        :math:`\xi_2 = 1`, plus :math:`N - 2` interior points.

    Returns
    -------
    np.ndarray
        Complex :math:`(N+1) \times (N+1)` matrix. The critical
        condition is :math:`\det M = 0`.

    Notes
    -----
    Slab and sphere share **everything** except ``geometry_sign``:

    * Moment recursion (:func:`._moment_arrays_complex`).
    * Collocation grid (passed via ``xis``).
    * Dispersion-root finder (the supercritical case for both
      gives :math:`\nu_0 = i u_0` with the same :math:`u_0`).
    * Bisection / Newton on :math:`R` such that :math:`\det M = 0`.

    This is the load-bearing structural fact validated by the SymPy
    derivations in :mod:`..origins.fn_sphere_derivations`.
    """
    if geometry_sign not in (+1, -1):
        raise ValueError(
            f"geometry_sign must be +1 (slab) or -1 (sphere), got {geometry_sign}"
        )
    if n_modes < 0:
        raise ValueError(f"n_modes must be ≥ 0, got {n_modes}")
    if len(xis) != n_modes + 1:
        raise ValueError(
            f"xis must have length n_modes+1 = {n_modes+1}, got {len(xis)}"
        )
    s = float(geometry_sign)
    M = np.empty((n_modes + 1, n_modes + 1), dtype=complex)
    for beta, xi in enumerate(xis):
        Bvec, Avec = _moment_arrays_complex(n_modes, xi, c)
        if xi == 0.0:
            exp_factor = 0.0 + 0.0j
        else:
            exp_factor = cmath.exp(-2.0 * R_mfp / xi)
        M[beta, :] = Bvec + s * exp_factor * Avec
    return M


def fn_determinant(
    c: float,
    R_mfp: float,
    geometry_sign: int,
    n_modes: int,
    xis: np.ndarray,
) -> complex:
    r"""Return :math:`\det M(R)` as a complex number.

    Convenience wrapper around :func:`assemble_fn_matrix` +
    :func:`numpy.linalg.det`. The critical condition is
    :math:`\det M(R) = 0`; the implementation typically uses a
    real-part bisection on this value.
    """
    M = assemble_fn_matrix(c, R_mfp, geometry_sign, n_modes, xis)
    return complex(np.linalg.det(M))
