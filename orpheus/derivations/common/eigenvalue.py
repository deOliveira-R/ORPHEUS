"""Shared eigenvalue solvers for analytical verification cases.

Three function families cover all derivation eigenvalue computations:

* :func:`kinf_homogeneous` / :func:`kinf_and_spectrum_homogeneous` —
  infinite medium, any group count (k∞ and the forward flux spectrum).
* :func:`kinf_and_adjoint_spectrum_homogeneous` — the same 0-D problem's
  **adjoint** (left-eigenvector) spectrum (campaign #276 A4 / P1.4).
* :func:`kinf_from_cp` — multi-region CP transport with P_inf matrices.

All support an optional ``sig_2`` / ``sig_2_mats`` parameter for
(n,2n) reactions.  When omitted, (n,2n) is treated as zero.
"""

from __future__ import annotations

import numpy as np


def _dense_transfer_matrix(m) -> np.ndarray:
    r"""Densify a group-transfer matrix: ndarray-like OR scipy-sparse.

    ``Mixture.SigS`` stores per-order transfer matrices SPARSE, and the
    pre-A4 reference accepted them (accidentally, through the
    ``ndarray − sparse.T → np.matrix`` overload chain).  A bare
    ``np.asarray(sparse, dtype=float)`` does NOT densify — it raises
    ``ValueError: setting an array element with a sequence`` — so the
    input tolerance is restored EXPLICITLY here (duck-typed on
    ``.toarray``), producing a plain ndarray rather than the old
    accidental ``np.matrix``.  (Caught by the A4 phase-end full tree:
    every A4 battery densified via ``.todense()`` first and was blind
    to the narrowed input domain.)
    """
    if hasattr(m, "toarray"):
        return np.asarray(m.toarray(), dtype=float)
    return np.asarray(m, dtype=float)


def _infinite_medium_matrices(
    sig_t: np.ndarray,
    sig_s: np.ndarray,
    nu_sig_f: np.ndarray,
    chi: np.ndarray,
    sig_2: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    r"""The 0-D loss/production pair :math:`(\mathbf{A}, \mathbf{F})`.

    .. math::

        \mathbf{A} = \text{diag}(\Sigma_t) - (\Sigma_s + 2\Sigma_2)^T,
        \qquad
        \mathbf{F} = \chi \otimes (\nu\Sigma_f)

    The single assembly site shared by the forward spectrum
    (:func:`kinf_and_spectrum_homogeneous`) and the adjoint spectrum
    (:func:`kinf_and_adjoint_spectrum_homogeneous`) — the two references
    MUST agree on the operator pair or their k's would not be comparable.
    """
    sig_s_eff = _dense_transfer_matrix(sig_s)
    if sig_2 is not None:
        sig_s_eff = sig_s_eff + 2.0 * _dense_transfer_matrix(sig_2)
    A = np.diag(np.asarray(sig_t, dtype=float)) - sig_s_eff.T
    F = np.outer(chi, nu_sig_f)
    return A, F


def kinf_homogeneous(
    sig_t: np.ndarray,
    sig_s: np.ndarray,
    nu_sig_f: np.ndarray,
    chi: np.ndarray,
    sig_2: np.ndarray | None = None,
) -> float:
    r"""Infinite-medium eigenvalue for any number of energy groups.

    Solves :math:`k = \lambda_{\max}(\mathbf{A}^{-1}\mathbf{F})` where

    .. math::

        \mathbf{A} = \text{diag}(\Sigma_t) - (\Sigma_s + 2\Sigma_2)^T

        \mathbf{F} = \chi \otimes (\nu\Sigma_f)

    For 1-group, this reduces to :math:`k = \nu\Sigma_f / \Sigma_a`.

    Parameters
    ----------
    sig_t : (ng,) total XS.
    sig_s : (ng, ng) P0 scattering matrix, ``[from, to]`` convention.
    nu_sig_f : (ng,) production XS (nu * sigma_f).
    chi : (ng,) fission spectrum.
    sig_2 : (ng, ng) or None — (n,2n) transfer matrix.
    """
    k, _ = kinf_and_spectrum_homogeneous(sig_t, sig_s, nu_sig_f, chi, sig_2)
    return k


def kinf_and_spectrum_homogeneous(
    sig_t: np.ndarray,
    sig_s: np.ndarray,
    nu_sig_f: np.ndarray,
    chi: np.ndarray,
    sig_2: np.ndarray | None = None,
) -> tuple[float, np.ndarray]:
    r"""Infinite-medium :math:`k_\infty` **and** the dominant flux spectrum.

    Phase-1.1 extension of :func:`kinf_homogeneous` that also returns
    the right eigenvector — the multigroup flux spectrum — of
    :math:`\mathbf{A}^{-1}\mathbf{F}`. This is the continuous
    reference for any test that wants to verify *flux shape* on a
    reflective (infinite-medium) problem, not just :math:`k_{\text{eff}}`.

    The spectrum is sign-normalised so every component is
    non-negative (physical flux) and normalised to unit
    :math:`\ell^{2}` norm. Callers that need a different
    normalisation (e.g. :math:`\phi_1 = 1` or volume-integrated
    total flux) should rescale the returned vector.

    Parameters
    ----------
    sig_t, sig_s, nu_sig_f, chi, sig_2 : see :func:`kinf_homogeneous`.

    Returns
    -------
    k : float
        Dominant eigenvalue, identical to :func:`kinf_homogeneous`.
    phi_spectrum : (ng,) ndarray
        Right eigenvector, :math:`\ell^{2}`-normalised, non-negative.
    """
    A, F = _infinite_medium_matrices(sig_t, sig_s, nu_sig_f, chi, sig_2)
    M = np.linalg.solve(A, F)

    eigvals, eigvecs = np.linalg.eig(M)
    real_vals = np.real(eigvals)
    dominant = int(np.argmax(real_vals))
    k = float(real_vals[dominant])
    phi = np.real(eigvecs[:, dominant])

    # Sign-normalise so the spectrum is non-negative (physical)
    if phi.sum() < 0:
        phi = -phi
    # Safety: if numerical noise leaves tiny negatives, clamp to zero.
    # For a well-posed downscatter problem the physical flux is
    # strictly non-negative; any component < -1e-12 indicates a real
    # numerical problem and is propagated as-is so the caller notices.
    phi = np.where(np.abs(phi) < 1e-14, 0.0, phi)

    norm = float(np.linalg.norm(phi))
    if norm == 0:
        raise ValueError(
            "Dominant eigenvector of A^{-1}F is identically zero — "
            "cross-section inputs are likely degenerate."
        )
    return k, phi / norm


def kinf_and_adjoint_spectrum_homogeneous(
    sig_t: np.ndarray,
    sig_s: np.ndarray,
    nu_sig_f: np.ndarray,
    chi: np.ndarray,
    sig_2: np.ndarray | None = None,
) -> tuple[float, np.ndarray]:
    r"""Infinite-medium :math:`k_\infty` **and** the dominant ADJOINT spectrum.

    The adjoint flux solves the DAGGERED eigenproblem — the defining law —

    .. math::

        \mathbf{A}^T\,\varphi^{*} \;=\; \tfrac{1}{k}\,\mathbf{F}^T\,\varphi^{*}
        \quad\Longleftrightarrow\quad
        (\mathbf{A}^T)^{-1}\mathbf{F}^T\,\varphi^{*} = k\,\varphi^{*} ,

    i.e. :math:`\varphi^*` is the dominant right eigenvector of
    :math:`(\mathbf{A}^T)^{-1}\mathbf{F}^T` — equivalently the **left**
    eigenvector of the REVERSED product :math:`\mathbf{F}\mathbf{A}^{-1}`,
    and **NOT** of :math:`\mathbf{M} = \mathbf{A}^{-1}\mathbf{F}`.

    .. warning::

        **The factor-order trap (Mode 12's transpose-side sibling,
        caught live at #276 A4).**  :math:`\mathbf{M}^T =
        \mathbf{F}^T\mathbf{A}^{-T}` is SIMILAR to
        :math:`(\mathbf{A}^T)^{-1}\mathbf{F}^T` (conjugation by
        :math:`\mathbf{A}^T` — the same algebra as the #226 step-5b
        ``A·(A⁻¹F)·A⁻¹ = FA⁻¹`` finding), so every k-level functional is
        blind to the swap — but the eigenVECTORS differ, and for the
        rank-1 :math:`\mathbf{F} = \chi\otimes\nu\Sigma_f` the wrong
        product's dominant eigenvector degenerates to EXACTLY
        :math:`\widehat{\nu\Sigma_f}` (since :math:`\mathbf{F}^T x
        \propto \nu\Sigma_f` for every :math:`x`) — a reference carrying
        ZERO :math:`\mathbf{A}`-physics.  The first spelling of this
        function used :math:`\text{eig}(\mathbf{M}^T)` (following the P6
        spec's P1.4 text) and its pin test encoded the same wrong law
        self-consistently; the structurally-independent SN daggered solve
        disagreed on first contact and exposed both.

    :math:`k` is IDENTICAL to the forward problem's; only the vector
    differs.  This is the closed-form flux-shape reference for the SN
    adjoint solve (campaign #276 A4, gate P1.4): an :math:`F^\dagger`
    role-swap error (the missing χ↔νΣf swap) moves the computed adjoint
    spectrum O(1) while every k-level functional stays exactly equal —
    the Mode-12 spectral invisibility this vector-level reference exists
    to break.

    :math:`(\mathbf{A}, \mathbf{F})` come from the SAME assembly as the
    forward spectrum (:func:`_infinite_medium_matrices`), so the two
    references pin the same operator pair.  Eigen-extraction rides
    :func:`~orpheus.numerics.eigenvalue.dominant_eigenpair` — the shared
    Perron–Frobenius primitive (structural independence from the SN
    solver lives in the INPUT assembly: dense XS matrices here vs the
    operator-algebra composition there, per the #276
    ``direct_eigenvalue`` ruling).  The returned vector carries the same
    convention as the forward spectrum (:math:`\ell^2`-normalised,
    non-negative sum), so the two are directly comparable.

    Returns
    -------
    k : float
        Dominant eigenvalue — identical to :func:`kinf_homogeneous`.
    phi_star_spectrum : (ng,) ndarray
        Left eigenvector (adjoint spectrum), :math:`\ell^{2}`-normalised.
    """
    # Local import keeps this Branch-1 reference module numpy-only at
    # import time; the shared extraction primitive is pulled on use.
    from orpheus.numerics.eigenvalue import dominant_eigenpair

    A, F = _infinite_medium_matrices(sig_t, sig_s, nu_sig_f, chi, sig_2)
    # The DAGGERED resolvent (Aᵀ)⁻¹Fᵀ — NOT (A⁻¹F)ᵀ; see the warning above.
    k, phi_star = dominant_eigenpair(np.linalg.solve(A.T, F.T))
    # dominant_eigenpair contracts sign (sum >= 0) but leaves scale
    # arbitrary — normalise explicitly, exactly as the forward sibling
    # does, so the two spectra are comparable by construction.
    return k, phi_star / float(np.linalg.norm(phi_star))


def kinf_from_cp(
    P_inf_g: np.ndarray,
    sig_t_all: np.ndarray,
    V_arr: np.ndarray,
    sig_s_mats: list[np.ndarray],
    nu_sig_f_mats: list[np.ndarray],
    chi_mats: list[np.ndarray],
    sig_2_mats: list[np.ndarray] | None = None,
) -> float:
    r"""Multi-region CP eigenvalue via dense eigensolver.

    Builds the generalized eigenvalue problem
    :math:`\mathbf{A}\phi = \frac{1}{k}\mathbf{B}\phi` and returns
    :math:`k = \lambda_{\max}(\mathbf{A}^{-1}\mathbf{B})`.

    The matrices are:

    .. math::

        A_{ig,jg'} = \delta_{ij}\delta_{gg'}\,\Sigma_{t,ig} V_i
                    - P_{ji,g}\, V_j\, (\Sigma_{s,j} + 2\Sigma_{2,j})_{g' \to g}

        B_{ig,jg'} = P_{ji,g}\, V_j\, \chi_{j,g}\, (\nu\Sigma_f)_{j,g'}

    Parameters
    ----------
    P_inf_g : (N_reg, N_reg, ng) — CP matrices per group.
    sig_t_all : (N_reg, ng) — total XS per region and group.
    V_arr : (N_reg,) — region volumes.
    sig_s_mats : list of (ng, ng) — scattering matrices per region.
    nu_sig_f_mats : list of (ng,) — production XS per region.
    chi_mats : list of (ng,) — fission spectrum per region.
    sig_2_mats : list of (ng, ng) or None — (n,2n) matrices per region.
    """
    N_reg = P_inf_g.shape[0]
    ng = P_inf_g.shape[2]
    dim = N_reg * ng

    A_mat = np.zeros((dim, dim))
    B_mat = np.zeros((dim, dim))

    for i_reg in range(N_reg):
        for g in range(ng):
            row = i_reg * ng + g
            A_mat[row, row] = sig_t_all[i_reg, g] * V_arr[i_reg]

            for j_reg in range(N_reg):
                for gp in range(ng):
                    col = j_reg * ng + gp
                    Pji = P_inf_g[j_reg, i_reg, g]

                    scat = sig_s_mats[j_reg][gp, g]
                    if sig_2_mats is not None:
                        scat = scat + 2.0 * sig_2_mats[j_reg][gp, g]

                    A_mat[row, col] -= Pji * V_arr[j_reg] * scat
                    B_mat[row, col] += (
                        Pji * V_arr[j_reg]
                        * chi_mats[j_reg][g]
                        * nu_sig_f_mats[j_reg][gp]
                    )

    M = np.linalg.solve(A_mat, B_mat)
    return float(np.max(np.real(np.linalg.eigvals(M))))
