r"""Eq. (3) and Eq. (4) matrix assembly + eigenvalue extraction.

This module is the bridge between the matrix-element evaluator
(:mod:`.carlvik_recurrences`) and the geometry-specific Branch-2
solvers in :mod:`...slab` and :mod:`...sphere`.

Eq. (3) of Dahl-Sjostrand 1979:

.. math::

   \big[A - 3\bar\mu(c-1) B\big] F = \frac{1}{c d} F .

Eq. (4) block-matrix linearization:

.. math::

   \begin{pmatrix} G & H \\ I & 0 \end{pmatrix}
   \begin{pmatrix} F \\ K \end{pmatrix}
   = \frac{1}{c}
   \begin{pmatrix} F \\ K \end{pmatrix} ,

with :math:`G = d(A + 3\bar\mu B)`, :math:`H = -3\bar\mu d B`,
:math:`K = c F`. Solving Eq. (4) as a standard non-symmetric
eigenproblem gives **all 2N eigenvalues** :math:`1/c_j` in one call.
"""
from __future__ import annotations

import numpy as np
import scipy.linalg

from .carlvik_recurrences import compute_A_matrix, compute_B_matrix


def assemble_eq4_block_matrix(
    a: float,
    n_modes: int,
    mu_bar: float,
    *,
    geometry: str,
    n_quad: int = 64,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    r"""Assemble the Eq. (4) block matrix and component matrices.

    Parameters
    ----------
    a : float
        Half-thickness (slab) or radius (sphere) in mean free paths.
    n_modes : int
        Number of Legendre modes :math:`N` in the expansion. The
        slab uses :math:`P_0, P_2, \ldots, P_{2(N-1)}`; the sphere
        uses :math:`P_1, P_3, \ldots, P_{2N-1}`.
    mu_bar : float
        Linearly-anisotropic mean cosine :math:`\bar\mu` (dimensionless,
        in :math:`[-1, +1]`). :math:`\bar\mu = 0` recovers the
        isotropic case.
    geometry : str
        ``"slab"`` (q=0, even-Legendre basis) or ``"sphere"`` (q=1,
        odd-Legendre basis).
    n_quad : int, default 64
        Gauss-Legendre quadrature order for matrix-element evaluation.

    Returns
    -------
    block : np.ndarray
        :math:`2N \times 2N` block matrix from Eq. (4):
        :math:`\begin{pmatrix} G & H \\ I & 0 \end{pmatrix}`.
    A : np.ndarray
        :math:`N \times N` matrix :math:`A_{m,n}(a)`.
    B : np.ndarray
        :math:`N \times N` matrix :math:`B_{m,n}(a)`.
    indices : np.ndarray
        Length-:math:`N` array of Legendre indices used.
    """
    if geometry == "slab":
        q = 0
        indices = np.arange(0, 2 * n_modes, 2, dtype=int)
    elif geometry == "sphere":
        q = 1
        indices = np.arange(1, 2 * n_modes + 1, 2, dtype=int)
    else:
        raise ValueError(f"geometry must be 'slab' or 'sphere', got {geometry!r}")

    if n_modes < 1:
        raise ValueError(f"n_modes must be >= 1, got {n_modes}")

    # The "natural" convention for the matrix elements (Carlvik 1968,
    # half-prefactor (2n+1)/2 with the Galerkin LHS = 2 F_m) yields
    # the eigenvalue equation A^my F = (4/(cd)) F. Dahl-Sjostrand
    # Eq. (3) prints the equation as A F = (1/(cd)) F, with their A
    # absorbing the factor of 4 (i.e. A^DS = A^my / 4). We follow
    # the printed Dahl-Sjostrand convention so eigenvalues come out
    # equal to the published table values.
    A_my = compute_A_matrix(a, indices, n_quad=n_quad)
    B_my = compute_B_matrix(a, indices, q=q, n_quad=n_quad)
    A = A_my / 4.0
    B = B_my / 4.0

    d = 2.0 * a
    G = d * (A + 3.0 * mu_bar * B)
    H = -3.0 * mu_bar * d * B
    I = np.eye(n_modes)
    Z = np.zeros((n_modes, n_modes))

    block = np.block([[G, H], [I, Z]])
    return block, A, B, indices


def solve_eq4_eigenproblem(
    a: float,
    n_modes: int,
    mu_bar: float,
    *,
    geometry: str,
    n_quad: int = 64,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    r"""Solve the Eq. (4) standard non-symmetric eigenproblem.

    Returns the **full eigenvalue spectrum** of the :math:`2N \times 2N`
    block matrix, with eigenvalues :math:`1/c_j` and corresponding
    eigenvectors. Sorted by :math:`\Re c_j` (real parts).

    Parameters
    ----------
    a : float
        Half-thickness or radius in mfp.
    n_modes : int
        Legendre expansion order.
    mu_bar : float
        Mean scattering cosine.
    geometry : str
        ``"slab"`` or ``"sphere"``.
    n_quad : int, default 64
        Quadrature order for matrix-element evaluation.

    Returns
    -------
    c_values : np.ndarray
        Length-:math:`2N` complex array of eigenvalues
        :math:`c_j = 1/\lambda_j`, sorted by ascending :math:`\Re c_j`.
    F_eigvecs : np.ndarray
        :math:`(N, 2N)` array of upper-half eigenvectors
        :math:`F_j` (one per column). The auxiliary
        :math:`K_j = c_j F_j` is the lower-half eigenvector,
        omitted from the return per Dahl-Sjostrand p. 119: "the
        useful length of the eigenvectors is the same as for
        isotropic scattering".
    block : np.ndarray
        The :math:`2N \times 2N` block matrix that was solved (returned
        for diagnostic purposes, e.g., residual checking).
    """
    block, _, _, _ = assemble_eq4_block_matrix(
        a, n_modes, mu_bar, geometry=geometry, n_quad=n_quad
    )

    # Solve the standard non-symmetric eigenproblem M v = λ v,
    # where λ = 1/c. Take c = 1/λ.
    eigvals, eigvecs = scipy.linalg.eig(block)
    # eigvals[j] = 1/c_j → c_j = 1/eigvals[j].
    # Filter out near-zero eigenvalues (would correspond to c → ∞).
    eps = 1e-300
    nonzero_mask = np.abs(eigvals) > eps
    c_vals = np.full_like(eigvals, np.inf, dtype=complex)
    c_vals[nonzero_mask] = 1.0 / eigvals[nonzero_mask]

    # Sort by real part of c.
    order = np.argsort(c_vals.real)
    c_sorted = c_vals[order]
    # Upper-half eigenvectors: indices 0..n_modes-1.
    F_eigvecs = eigvecs[:n_modes, order]

    return c_sorted, F_eigvecs, block
