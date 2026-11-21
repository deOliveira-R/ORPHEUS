r"""Test-side dense eigen oracle of the generalized eigenproblem :math:`A\,\varphi = \tfrac1k F\,\varphi` (moved from
``orpheus.numerics.eigenvalue``, dead-B rows 49, 130 and 131).

Production solves the posed problem by :func:`~orpheus.numerics.eigenvalue.power_iteration` (iterative, any operator that can be
applied) and the homogeneous solver reads :math:`k_\infty` off the rank-one fission dyad; neither forms a dense resolvent. The dense
engine lives here, as the oracle the gates compare production against: ``direct_eigenvalue`` forms the resolvent
:math:`A^{-1}F` with :func:`numpy.linalg.solve` and ``dominant_eigenpair`` extracts the Perron--Frobenius pair from it, both in
ONE module, so ``direct_eigenvalue`` reaches ``dominant_eigenpair`` through a module-global lookup (a spy planted on this module's
``dominant_eigenpair`` is a one-home proof: it moves ``direct_eigenvalue``).
"""

from __future__ import annotations

import numpy as np


def _sign_normalised(v: np.ndarray) -> np.ndarray:
    r"""Sign-normalise an eigenvector to the physical, non-negative convention.

    A criticality eigenvector is a flux distribution, defined up to scale and in particular up to SIGN. The convention is
    ``v.sum() >= 0`` (the physical, non-negative spectrum); this helper is its single source here, shared by
    :func:`dominant_eigenpair` and through it :func:`direct_eigenvalue`.
    """
    return v if v.sum() >= 0.0 else -v


def dominant_eigenpair(
    M: np.ndarray,
    *,
    imag_tol: float = 1e-9,
) -> tuple[float, np.ndarray]:
    r"""Dominant (Perron--Frobenius) eigenpair of a materialized resolvent.

    Given the dense resolvent :math:`M = A^{-1}F`, however it was formed, return the dominant eigenpair with the criticality
    contract enforced: by Krein--Rutman / Perron--Frobenius the fundamental mode of a well-posed criticality resolvent is the
    unique non-negative eigenvector and its eigenvalue is **real and positive**, so a complex dominant eigenvalue is rejected as a
    malformed problem (fail loud, never return a non-eigenvalue).

    Parameters
    ----------
    M : np.ndarray
        The materialized resolvent (square, ``(n, n)``, real) whose dominant eigenpair is sought.
    imag_tol : float, optional
        Tolerance on the dominant eigenvalue's imaginary part before it is rejected as complex. Default ``1e-9``.

    Returns
    -------
    k : float
        The dominant eigenvalue :math:`\lambda_{\max}(M)` (largest real part; real by the Perron--Frobenius contract).
    phi : np.ndarray
        The corresponding right eigenvector (real-dtype, ``(n,)``), sign-normalised so ``phi.sum() >= 0``; its absolute scale is
        arbitrary.

    Raises
    ------
    ValueError
        If ``M`` is not a square 2-D matrix, or the dominant eigenvalue is complex (``|Im λ| > imag_tol``).
    """
    M = np.asarray(M, dtype=float)
    if M.ndim != 2 or M.shape[0] != M.shape[1]:
        raise ValueError(
            f"dominant_eigenpair: M must be a square 2-D matrix; got shape {M.shape}."
        )
    eigvals, eigvecs = np.linalg.eig(M)

    # The dominant (Perron–Frobenius) mode: the largest real part.
    dominant = int(np.argmax(np.real(eigvals)))
    lam = eigvals[dominant]
    if abs(float(np.imag(lam))) > imag_tol:
        raise ValueError(
            f"dominant_eigenpair: the dominant eigenvalue is complex "
            f"({lam:.6g}); the resolvent A⁻¹F of a well-posed criticality "
            f"problem has a real dominant eigenvalue (Perron–Frobenius). "
            f"A complex dominant signals a malformed problem."
        )
    k = float(np.real(lam))
    phi = _sign_normalised(np.real(eigvecs[:, dominant]))
    return k, phi


def direct_eigenvalue(
    A: np.ndarray,
    F: np.ndarray,
    *,
    imag_tol: float = 1e-9,
) -> tuple[float, np.ndarray]:
    r"""Exact dominant eigenpair of the generalized eigenproblem :math:`A\,\varphi = \tfrac{1}{k}\,F\,\varphi`, via the dense
    resolvent.

    The **direct (non-iterative) sibling of** :func:`~orpheus.numerics.eigenvalue.power_iteration`: it forms the dense resolvent
    :math:`A^{-1}F` and returns the EXACT dominant eigenpair in one LAPACK call,

    .. math::

        k \;=\; \lambda_{\max}\!\bigl(A^{-1}F\bigr), \qquad A^{-1}F\,\varphi \;=\; k\,\varphi .

    It is the right realization for **small, densifiable** operators (0-D infinite-medium spectra, few-group / few-region
    problems), where the iterative engine would only approximate an answer the dense solve gives exactly. Both engines solve the
    SAME posed problem :math:`(A_{\rm loss}, M) = (A, F)`; they differ only in exact-dense vs iterative.

    Parameters
    ----------
    A : np.ndarray
        The loss matrix :math:`A` (square, ``(n, n)``). Must be invertible: a singular loss matrix is a malformed problem and
        :func:`numpy.linalg.solve` raises :class:`numpy.linalg.LinAlgError`.
    F : np.ndarray
        The production matrix :math:`F` (same shape as ``A``).
    imag_tol : float, optional
        Tolerance on the dominant eigenvalue's imaginary part before it is rejected as complex. Default ``1e-9``.

    Returns
    -------
    k : float
        The dominant eigenvalue :math:`\lambda_{\max}(A^{-1}F)` (real).
    phi : np.ndarray
        The corresponding right eigenvector (real-dtype, ``(n,)``), sign-normalised so ``phi.sum() >= 0``.

    Raises
    ------
    ValueError
        If ``A`` / ``F`` are not square matrices of equal shape, or the dominant eigenvalue is complex (raised by
        :func:`dominant_eigenpair`, the extraction home).
    numpy.linalg.LinAlgError
        If ``A`` is singular (propagated from :func:`numpy.linalg.solve`).
    """
    A = np.asarray(A, dtype=float)
    F = np.asarray(F, dtype=float)
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError(
            f"direct_eigenvalue: A must be a square 2-D matrix; got shape {A.shape}."
        )
    if F.shape != A.shape:
        raise ValueError(
            f"direct_eigenvalue: A and F must have the same shape; "
            f"got A {A.shape}, F {F.shape}."
        )

    # The dense resolvent A⁻¹F.  np.linalg.solve raises LinAlgError on a singular A, left to fail loud.  The
    # Perron–Frobenius extraction lives in ONE home, ``dominant_eigenpair`` of this module (a global-name lookup).
    return dominant_eigenpair(np.linalg.solve(A, F), imag_tol=imag_tol)
