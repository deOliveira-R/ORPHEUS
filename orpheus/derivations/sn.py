r"""SN (Discrete Ordinates) eigenvalue reference solutions.

Three problem families:

1. **Homogeneous slab** (:func:`_derive_sn_homogeneous`) — analytical
   :math:`k_\infty` from the transport equation with flat flux and
   reflective BCs. Pure T1, closed form.

2. **Heterogeneous spatial operator** (MMS) — verified in
   :mod:`orpheus.derivations.sn_mms` with a manufactured solution on
   smooth :math:`\Sigma(x)` (Phase 2.1a).

3. **Heterogeneous eigenvalue** (this module, Phase 2.1b) — Case
   singular-eigenfunction expansion of the discrete-:math:`S_N` slope
   matrix on a two-region reflective slab, giving a semi-analytical
   :math:`k_\text{eff}` reference that stands on its own mathematical
   merits. The reference uses the **same Gauss–Legendre quadrature
   order** as the solver under test, so it is the exact discrete-SN
   eigenvalue — not the continuous-angle limit — and should agree with
   ``solve_sn`` to machine precision once spatial discretisation
   converges. See :func:`derive_sn_heterogeneous_continuous` and the
   "Case singular-eigenfunction method" section of
   :doc:`/theory/discrete_ordinates`.

The Case-method implementation follows the same real-basis
mode-decomposition pattern as :mod:`orpheus.derivations.diffusion`:

- Per region, diagonalise the slope matrix
  :math:`\mathbf S(k) = \mathbf M^{-1}(\mathbf K(k) - \Sigma_t\mathbf I)`
  via :func:`numpy.linalg.eig`, where :math:`\mathbf M = \mathrm{diag}(\mu_n)`.
- Each real eigenvalue :math:`\lambda` gives one real basis mode
  :math:`\exp(\lambda(x - x_\text{anchor}))\,\mathbf v`, anchored at the
  nearer region edge (right edge for :math:`\lambda \ge 0`, left edge
  for :math:`\lambda < 0`) to keep the mode bounded on the region.
- Each complex-conjugate pair :math:`(\lambda, \lambda^*)` gives two
  real basis modes via the standard
  :math:`\exp(\Re\lambda\,\xi)\,(\cos(\Im\lambda\,\xi)\,\mathbf v_R -
  \sin(\Im\lambda\,\xi)\,\mathbf v_I)` combination, with the same
  anchored argument :math:`\xi = x - x_\text{anchor}` to bound the
  exponential envelope.
- Assemble the :math:`2N \times 2N` real matching matrix
  :math:`\mathbf C(k)` enforcing reflective BCs at both outer edges
  and full angular continuity at the material interface.
- Solve :math:`\det(\mathbf C(k)) = 0` by coarse scan + brentq refinement,
  and **physically validate** every candidate root by reconstructing the
  BC residuals from the null vector. This rejects spurious sign changes
  caused by the non-continuous eigenvalue ordering of
  :func:`numpy.linalg.eig` across :math:`k` crossings of the
  per-region spectrum.

The Phase 2.1b target problem is 1-group, two-region (fuel + moderator)
at reflective BCs, :math:`S_8` quadrature — the configuration for which
the ERR-025 diamond-difference bug in ``_sweep_1d_cumprod`` was
diagnosed and fixed. This module produces that reference; the consumer
test is :func:`tests.sn.test_heterogeneous_transport.test_sn_2region_reflective_case_eigenvalue`.

Multigroup extension is a mechanical generalisation: the slope matrix
becomes block-diagonal in ordinate with per-group off-diagonal blocks
from scattering and fission. The current implementation is 1-group only;
multigroup is a Phase 2.1b+ follow-up.

See :doc:`/theory/discrete_ordinates` "Case singular-eigenfunction
method" section for the full mathematical treatment with equation
labels, and :doc:`/verification/reference_solutions` for the campaign
philosophy.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from ._eigenvalue import kinf_homogeneous
from ._reference import (
    ContinuousReferenceSolution,
    ProblemSpec,
    Provenance,
)
from ._types import VerificationCase
from ._xs_library import get_mixture, get_xs


# ═══════════════════════════════════════════════════════════════════════
# Legacy homogeneous cases (T1 analytical)
# ═══════════════════════════════════════════════════════════════════════

def _derive_sn_homogeneous(ng_key: str) -> VerificationCase:
    """Derive SN eigenvalue for a homogeneous slab from the transport equation."""
    xs = get_xs("A", ng_key)
    ng = len(xs["sig_t"])

    k_val = kinf_homogeneous(
        sig_t=xs["sig_t"], sig_s=xs["sig_s"],
        nu_sig_f=xs["nu"] * xs["sig_f"], chi=xs["chi"],
    )

    if ng == 1:
        latex = (
            r"From the 1D S\ :sub:`N` equation with "
            r":math:`\partial\psi_m/\partial x = 0` (homogeneous, reflective BCs):"
            "\n\n"
            r".. math::" "\n"
            rf"   k = \nu\Sigma_f / \Sigma_a = {k_val:.6f}"
        )
    else:
        latex = (
            r"Multi-group S\ :sub:`N` homogeneous: flat flux reduces to "
            r":math:`k = \lambda_{\max}(\mathbf{A}^{-1}\mathbf{F})`."
            "\n\n"
            r".. math::" "\n"
            rf"   k_\infty = {k_val:.10f}"
        )

    labels: list[str] = ["transport-cartesian", "reflective-bc"]
    if ng == 1:
        labels.append("one-group-kinf")
    else:
        labels += ["matrix-eigenvalue", "mg-balance", "multigroup"]

    return VerificationCase(
        name=f"sn_slab_{ng}eg_1rg",
        k_inf=k_val,
        method="sn",
        geometry="slab",
        n_groups=ng,
        n_regions=1,
        materials={0: get_mixture("A", ng_key)},
        geom_params={},
        latex=latex,
        description=f"SN 1D reflective slab, {ng}G homogeneous — from transport equation",
        tolerance="< 1e-8",
        vv_level="L1",
        equation_labels=tuple(labels),
    )


def all_cases() -> list[VerificationCase]:
    """Return analytical SN cases (homogeneous only).

    Heterogeneous SN verification is now covered by the
    Phase-2.1a MMS continuous reference in
    :mod:`orpheus.derivations.sn_mms` (spatial operator) and
    the Phase-2.1b Case singular-eigenfunction reference in
    :func:`continuous_cases` below (eigenvalue). This module no
    longer produces heterogeneous ``VerificationCase`` objects.
    """
    return [_derive_sn_homogeneous(ng_key) for ng_key in ["1g", "2g", "4g"]]


# ═══════════════════════════════════════════════════════════════════════
# Phase 2.1b — Case singular-eigenfunction method (heterogeneous slab)
# ═══════════════════════════════════════════════════════════════════════
#
# The following functions implement a real-basis mode decomposition of
# the discrete-SN slope matrix, producing a semi-analytical eigenvalue
# reference for two-region reflective slabs. The pattern mirrors the
# Phase 1.2 diffusion transcendental reference
# (:func:`orpheus.derivations.diffusion.derive_2rg_continuous`); all
# four diffusion dead ends documented in
# ``docs/theory/diffusion_1d.rst`` "Investigation history" apply
# verbatim here and are handled the same way.


def _build_sn_slope_matrix_1g(
    sig_t: float,
    sig_s: float,
    nu_sig_f: float,
    k: float,
    mu: np.ndarray,
    weights: np.ndarray,
) -> np.ndarray:
    r"""Build the 1-group discrete-SN slope matrix :math:`\mathbf S(k)`.

    From the 1D slab transport equation
    :math:`\mu_n\,\partial\psi_n/\partial x + \Sigma_t\,\psi_n =
    (1/W)\,(\Sigma_s + \nu\Sigma_f/k)\,\phi`, with
    :math:`\phi = \sum_m w_m\,\psi_m` and :math:`W = \sum_m w_m`, the
    per-ordinate ODE is

    .. math::

        \mu_n\,\frac{d\psi_n}{dx}
          = -\Sigma_t\,\psi_n + \frac{c_\text{eff}(k)}{W}\,\sum_m w_m\,\psi_m

    where :math:`c_\text{eff}(k) = \Sigma_s + \nu\Sigma_f/k`. Stacking
    the angular flux into :math:`\mathbf y \in \mathbb R^N` with
    components :math:`y_n = \psi_n`, the system is
    :math:`d\mathbf y/dx = \mathbf S(k)\,\mathbf y` with

    .. math::

        \mathbf S(k)[n, m]
          = \frac{1}{\mu_n}\left(
               -\Sigma_t\,\delta_{nm}
               + \frac{c_\text{eff}(k)}{W}\,w_m
            \right).

    Note the **row-scaling** :math:`1/\mu_n` — the slope matrix is
    generally non-symmetric even for symmetric GL quadrature. Its
    spectrum dictates the per-region spatial mode structure.
    """
    N = len(mu)
    W = float(weights.sum())
    c_eff = sig_s + nu_sig_f / k

    M_inv = 1.0 / mu                           # (N,)
    K = (c_eff / W) * weights[None, :]         # (1, N) → broadcast to (N, N)
    K_minus_T = K - sig_t * np.eye(N)          # (N, N)
    S = M_inv[:, None] * K_minus_T             # row-scale by 1/mu
    return S


def _region_spatial_modes_sn_1g(
    sig_t: float,
    sig_s: float,
    nu_sig_f: float,
    k: float,
    mu: np.ndarray,
    weights: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    r"""Diagonalise the 1-group SN slope matrix :math:`\mathbf S(k)`.

    Returns ``(eigvals, eigvecs)`` from :func:`numpy.linalg.eig` on
    :math:`\mathbf S(k)`. The eigenvalues are generally complex because
    :math:`\mathbf S(k)` is non-symmetric; for :math:`c_\text{eff}(k) < 1`
    (subcritical region) they come in :math:`\pm` real pairs, for
    :math:`c_\text{eff}(k) > 1` (supercritical region, e.g. pure fuel at
    ``k`` below ``k_inf_fuel``) some pairs become complex-conjugate.

    :func:`numpy.linalg.eig` ordering is not a continuous function of
    ``k`` — at crossings the eigenvalue labels permute discontinuously,
    which drives the spurious-root pathology that ``_physical_validation``
    filters out downstream (see Phase 1.2 dead end #3 in the diffusion
    investigation history).
    """
    S = _build_sn_slope_matrix_1g(sig_t, sig_s, nu_sig_f, k, mu, weights)
    return np.linalg.eig(S)


def _pair_complex_eigenvalues(
    lams: np.ndarray,
    tol: float = 1e-10,
) -> list[tuple[str, list[int]]]:
    r"""Group eigenvalues into real entries and complex-conjugate pairs.

    Returns a list of ``("real", [i])`` and ``("complex", [i, j])``
    tuples that together partition the ``len(lams)`` eigenvalue indices.
    Complex-conjugate pairs are detected by matching
    :math:`\Re\lambda_i \approx \Re\lambda_j` and
    :math:`\Im\lambda_i \approx -\Im\lambda_j`, both with tolerance
    ``tol``. Real eigenvalues are those with ``|Im λ| < tol``.

    A complete pairing is required; the function raises
    :class:`RuntimeError` if any complex eigenvalue cannot be matched
    with its conjugate, which would indicate a numerical artefact in
    :func:`numpy.linalg.eig` rather than a physical spectrum.
    """
    N = len(lams)
    used = [False] * N
    groups: list[tuple[str, list[int]]] = []
    for i in range(N):
        if used[i]:
            continue
        if abs(lams[i].imag) < tol:
            groups.append(("real", [i]))
            used[i] = True
        else:
            for j in range(i + 1, N):
                if (
                    not used[j]
                    and abs(lams[j].imag) > tol
                    and abs(lams[j].real - lams[i].real) < tol
                    and abs(lams[j].imag + lams[i].imag) < tol
                ):
                    groups.append(("complex", [i, j]))
                    used[i] = True
                    used[j] = True
                    break
            if not used[i]:
                raise RuntimeError(
                    f"Eigenvalue lams[{i}] = {lams[i]} has no conjugate "
                    "partner in the spectrum; numpy.linalg.eig returned "
                    "an incomplete complex pairing."
                )
    return groups


def _region_basis_at_x(
    sig_t: float,
    sig_s: float,
    nu_sig_f: float,
    k: float,
    L: float,
    x_rel: float,
    mu: np.ndarray,
    weights: np.ndarray,
) -> np.ndarray:
    r"""Evaluate the :math:`N` real bounded basis modes at ``x_rel``.

    ``x_rel`` is the position relative to the region's left edge, so
    :math:`0 \le \text{x\_rel} \le L`. Returns an ``(N, N)`` matrix
    ``psi[:, j]`` whose columns are the :math:`N`-vector angular-flux
    values of the :math:`j`-th real basis mode at the requested point.

    The real basis is built from the per-region eigendecomposition of
    :math:`\mathbf S(k)`:

    - **Real eigenvalue** :math:`\lambda` — one mode,
      :math:`\exp(\lambda\,\xi)\,\mathbf v` with anchored argument
      :math:`\xi = x - L` if :math:`\lambda \ge 0` (anchor at right edge,
      mode bounded by 1 at :math:`x = L`) and :math:`\xi = x` if
      :math:`\lambda < 0` (anchor at left edge).
    - **Complex-conjugate pair** :math:`(\lambda, \lambda^*)` — two
      real modes,

      .. math::

          \psi^{(c)}(x) &= \exp(\Re\lambda\,\xi)\,
            (\cos(\Im\lambda\,\xi)\,\mathbf v_R
             - \sin(\Im\lambda\,\xi)\,\mathbf v_I) \\
          \psi^{(s)}(x) &= \exp(\Re\lambda\,\xi)\,
            (\sin(\Im\lambda\,\xi)\,\mathbf v_R
             + \cos(\Im\lambda\,\xi)\,\mathbf v_I)

      with :math:`\mathbf v = \mathbf v_R + i\mathbf v_I` the complex
      eigenvector and :math:`\xi` anchored the same way as real modes,
      via :math:`\Re\lambda`'s sign.

    Every mode is bounded by :math:`|\mathbf v|` on the region, so the
    matching matrix assembled from these values has :math:`\mathcal O(1)`
    entries and :math:`\det\mathbf C(k)` is computed without
    catastrophic cancellation — exactly the Phase 1.2 fix for the
    diffusion analogue.
    """
    lams, V = _region_spatial_modes_sn_1g(
        sig_t, sig_s, nu_sig_f, k, mu, weights,
    )
    groups = _pair_complex_eigenvalues(lams)

    N = len(mu)
    psi = np.zeros((N, N))
    col = 0
    for kind, idxs in groups:
        if kind == "real":
            i = idxs[0]
            lam = float(lams[i].real)
            v = V[:, i].real
            if lam >= 0:
                xi = x_rel - L
            else:
                xi = x_rel
            psi[:, col] = np.exp(lam * xi) * v
            col += 1
        else:
            i = idxs[0]
            lam = lams[i]
            v = V[:, i]
            a = lam.real
            b = lam.imag
            if a >= 0:
                xi = x_rel - L
            else:
                xi = x_rel
            env = np.exp(a * xi)
            c_arg = np.cos(b * xi)
            s_arg = np.sin(b * xi)
            psi[:, col]     = env * (c_arg * v.real - s_arg * v.imag)
            psi[:, col + 1] = env * (s_arg * v.real + c_arg * v.imag)
            col += 2
    return psi


def _assemble_matching_matrix_sn_1g(
    xs_A: dict,
    xs_B: dict,
    H_A: float,
    H_B: float,
    k: float,
    mu: np.ndarray,
    weights: np.ndarray,
) -> np.ndarray:
    r"""Assemble the :math:`2N \times 2N` matching matrix for a 2-region
    reflective slab.

    The eigenvalue problem has :math:`2N` unknowns — :math:`N` mode
    coefficients per region — and :math:`2N` linear constraints:

    1. **Reflective BC at** :math:`x = 0`:
       :math:`\psi^A_n(0) = \psi^A_{N-1-n}(0)` for each positive
       ordinate :math:`n \in [0, N/2)`. Gauss–Legendre ordinates are
       sorted in ascending :math:`\mu`, so the reflection-partner
       index of ordinate :math:`n` is :math:`N-1-n`. :math:`N/2`
       equations.
    2. **Interface continuity at** :math:`x = H_A`:
       :math:`\psi^A_n(H_A) = \psi^B_n(H_A)` for **all** :math:`n`.
       :math:`N` equations.
    3. **Reflective BC at** :math:`x = H_A + H_B`:
       :math:`\psi^B_n(L) = \psi^B_{N-1-n}(L)` for each positive
       ordinate. :math:`N/2` equations.

    Total: :math:`N + N = 2N` equations in :math:`2N` unknowns. The
    coefficient vector is :math:`[\mathbf c^A;\,\mathbf c^B]`, and
    the eigenvalue condition is :math:`\det\mathbf C(k) = 0`.

    Notes
    -----
    This mirrors :func:`orpheus.derivations.diffusion._assemble_matching_matrix`
    in structure, with three differences:

    - The BC is **reflective**, not vacuum (angular-flux pairing
      instead of flux nulling).
    - The operator is **first-order** in :math:`x`, so only the angular
      flux :math:`\psi` is matched — not the derivative. Interface
      continuity of :math:`\psi` already implies continuity of
      :math:`\mu_n\,\partial_x\psi_n` via the transport equation, since
      :math:`\Sigma_t` changes discontinuously across the interface
      but the source term
      :math:`(c_\text{eff}/W)\sum_m w_m\psi_m` is continuous in
      :math:`\psi`.
    - The dimension is :math:`2N`, not :math:`2\cdot 2N` as in the
      second-order diffusion case.
    """
    N = len(mu)
    n_half = N // 2
    ref_idx = np.arange(N - 1, -1, -1)  # GL: ordinate n ↔ N-1-n

    psi_A_left = _region_basis_at_x(
        xs_A["sig_t"], xs_A["sig_s"], xs_A["nu_sig_f"],
        k, H_A, 0.0, mu, weights,
    )
    psi_A_right = _region_basis_at_x(
        xs_A["sig_t"], xs_A["sig_s"], xs_A["nu_sig_f"],
        k, H_A, H_A, mu, weights,
    )
    psi_B_left = _region_basis_at_x(
        xs_B["sig_t"], xs_B["sig_s"], xs_B["nu_sig_f"],
        k, H_B, 0.0, mu, weights,
    )
    psi_B_right = _region_basis_at_x(
        xs_B["sig_t"], xs_B["sig_s"], xs_B["nu_sig_f"],
        k, H_B, H_B, mu, weights,
    )

    C = np.zeros((2 * N, 2 * N))

    # Rows 0 .. n_half-1: reflective BC at x=0
    #   ψ_A(0, μ_{N-1-n}) - ψ_A(0, μ_n) = 0  for positive ordinates n
    # i.e. the positive-direction flux equals the reflected negative-
    # direction flux at the left edge. Implementing in GL index
    # convention: the positive half of ordinates is [n_half, N-1],
    # its reflection partners (negative half) are [n_half-1, ..., 0].
    # The constraint is psi_A_left[n, :] − psi_A_left[ref_idx[n], :] = 0
    # applied for n in the positive half.
    for i in range(n_half):
        n_pos = n_half + i
        n_neg = ref_idx[n_pos]  # = n_half - 1 - i
        C[i, 0:N] = psi_A_left[n_pos, :] - psi_A_left[n_neg, :]

    # Rows n_half .. n_half+N-1: interface continuity at x=H_A, all ordinates
    #   ψ_A(H_A, μ_n) − ψ_B(0, μ_n) = 0
    for n in range(N):
        C[n_half + n, 0:N] = psi_A_right[n, :]
        C[n_half + n, N:2 * N] = -psi_B_left[n, :]

    # Rows n_half + N .. 2N-1: reflective BC at x=L
    for i in range(n_half):
        row = n_half + N + i
        n_pos = n_half + i
        n_neg = ref_idx[n_pos]
        C[row, N:2 * N] = psi_B_right[n_pos, :] - psi_B_right[n_neg, :]

    return C


def _extract_real_null_vector_sn(
    C: np.ndarray, tol: float = 1e-3,
) -> np.ndarray:
    r"""Return the real 1-D null vector of ``C`` via SVD.

    ``C`` is the :math:`2N \times 2N` real matching matrix. At the
    true eigenvalue ``k``, its smallest singular value collapses to
    :math:`\mathcal O(\epsilon_{\text{cond}})`, and the corresponding
    right singular vector gives the coefficient vector in the
    :math:`2N`-dimensional real mode basis.
    """
    U, s, Vh = np.linalg.svd(C)
    if s[-1] > tol * s[0]:
        raise RuntimeError(
            f"Matching matrix is not singular: "
            f"s[-1]/s[0] = {s[-1] / s[0]:.2e} > {tol}"
        )
    return Vh[-1, :]


def _evaluate_flux_at(
    xs_A: dict,
    xs_B: dict,
    H_A: float,
    H_B: float,
    k: float,
    x: float,
    c_A: np.ndarray,
    c_B: np.ndarray,
    mu: np.ndarray,
    weights: np.ndarray,
) -> tuple[np.ndarray, float]:
    r"""Evaluate the angular flux and the scalar flux at a point ``x``.

    Returns ``(psi, phi)`` where ``psi`` is the length-:math:`N` angular
    flux vector and ``phi`` is the weighted-sum scalar flux. The caller
    provides the region coefficient vectors ``c_A``, ``c_B`` extracted
    from the null space of the matching matrix.
    """
    if x <= H_A:
        basis = _region_basis_at_x(
            xs_A["sig_t"], xs_A["sig_s"], xs_A["nu_sig_f"],
            k, H_A, x, mu, weights,
        )
        psi = basis @ c_A
    else:
        basis = _region_basis_at_x(
            xs_B["sig_t"], xs_B["sig_s"], xs_B["nu_sig_f"],
            k, H_B, x - H_A, mu, weights,
        )
        psi = basis @ c_B
    phi = float(weights @ psi)
    return psi, phi


def _physical_validation_sn(
    xs_A: dict,
    xs_B: dict,
    H_A: float,
    H_B: float,
    k: float,
    mu: np.ndarray,
    weights: np.ndarray,
    tol: float = 1e-6,
) -> bool:
    r"""Return True if ``k`` is a physical eigenvalue of the 2-region
    reflective SN slab.

    Rebuilds :math:`\mathbf C(k)`, extracts the null vector via SVD,
    and explicitly reconstructs:

    1. Reflective-BC residual at :math:`x = 0` — every positive
       ordinate's angular flux must equal its reflection-partner
       negative ordinate's angular flux.
    2. Interface continuity at :math:`x = H_A` — the region-A and
       region-B angular flux vectors must agree to machine precision.
    3. Reflective-BC residual at :math:`x = H_A + H_B` — same as (1).

    Rejects spurious sign changes caused by the non-continuous
    eigenvalue ordering of :func:`numpy.linalg.eig` across critical
    ``k`` values. This is the SN analogue of the Phase 1.2 diffusion
    ``_physical_validation`` — same pattern, different operator.
    """
    N = len(mu)
    n_half = N // 2
    ref_idx = np.arange(N - 1, -1, -1)

    try:
        C = _assemble_matching_matrix_sn_1g(
            xs_A, xs_B, H_A, H_B, k, mu, weights,
        )
        nv = _extract_real_null_vector_sn(C, tol=1e-3)
    except RuntimeError:
        return False

    c_A = nv[:N]
    c_B = nv[N:]

    L = H_A + H_B
    psi_0, phi_0 = _evaluate_flux_at(
        xs_A, xs_B, H_A, H_B, k, 0.0, c_A, c_B, mu, weights,
    )
    psi_if_A, _ = _evaluate_flux_at(
        xs_A, xs_B, H_A, H_B, k, H_A, c_A, c_B, mu, weights,
    )
    psi_if_B, _ = _evaluate_flux_at(
        xs_A, xs_B, H_A, H_B, k, H_A + 1e-12, c_A, c_B, mu, weights,
    )
    psi_L, _ = _evaluate_flux_at(
        xs_A, xs_B, H_A, H_B, k, L, c_A, c_B, mu, weights,
    )
    # Scale by the peak angular flux magnitude so the tolerance is
    # dimensionless and not sensitive to the null-vector normalisation.
    scale = max(
        float(np.max(np.abs(psi_if_A))),
        float(np.max(np.abs(psi_0))),
        1e-30,
    )

    # Reflective-BC residuals: for each positive ordinate, the
    # positive and reflection-partner-negative values must match.
    refl_left = 0.0
    refl_right = 0.0
    for i in range(n_half):
        n_pos = n_half + i
        n_neg = ref_idx[n_pos]
        refl_left += abs(psi_0[n_pos] - psi_0[n_neg])
        refl_right += abs(psi_L[n_pos] - psi_L[n_neg])
    refl_left /= scale
    refl_right /= scale
    interface = float(np.linalg.norm(psi_if_A - psi_if_B)) / scale

    return bool(
        refl_left < tol
        and interface < tol
        and refl_right < tol
    )


def _solve_2region_reflective_sn_eigenvalue(
    xs_A: dict,
    xs_B: dict,
    H_A: float,
    H_B: float,
    mu: np.ndarray,
    weights: np.ndarray,
    k_low: float = 0.3,
    k_high: float = 2.0,
    n_scan: int = 400,
    brentq_xtol: float = 1e-14,
) -> float:
    r"""Find the fundamental :math:`k_\text{eff}` of a 2-region
    reflective SN slab.

    Builds :math:`\mathbf C(k)` on a coarse ``k`` scan, brackets every
    sign change in :math:`\det\mathbf C(k)` with :func:`scipy.optimize.brentq`,
    and physically validates each candidate via
    :func:`_physical_validation_sn`. The fundamental is the largest
    validated root.
    """
    def det_c(k: float) -> float:
        C = _assemble_matching_matrix_sn_1g(
            xs_A, xs_B, H_A, H_B, k, mu, weights,
        )
        return float(np.linalg.det(C))

    ks = np.linspace(k_low, k_high, n_scan)
    dets = np.array([det_c(k) for k in ks])

    candidates: list[float] = []
    for i in range(len(ks) - 1):
        if dets[i] * dets[i + 1] < 0:
            try:
                root = brentq(det_c, ks[i], ks[i + 1], xtol=brentq_xtol)
                candidates.append(float(root))
            except Exception:
                pass

    validated = [
        k for k in candidates
        if _physical_validation_sn(
            xs_A, xs_B, H_A, H_B, k, mu, weights,
        )
    ]

    if not validated:
        raise RuntimeError(
            f"No validated SN eigenvalue in [{k_low}, {k_high}] for "
            f"(Σt_A, Σs_A, νΣf_A) = ({xs_A['sig_t']}, {xs_A['sig_s']}, "
            f"{xs_A['nu_sig_f']}), (Σt_B, Σs_B, νΣf_B) = ({xs_B['sig_t']}, "
            f"{xs_B['sig_s']}, {xs_B['nu_sig_f']}), H_A={H_A}, H_B={H_B}. "
            f"Found {len(candidates)} sign-change candidates but none "
            "passed physical validation. Widen the scan range or check "
            "the cross sections and quadrature."
        )
    return max(validated)


# ═══════════════════════════════════════════════════════════════════════
# Public API: continuous reference for 1-group 2-region reflective slab
# ═══════════════════════════════════════════════════════════════════════

def _xs_1g_dict(mixture) -> dict:
    """Flatten a 1-group :class:`Mixture` to a plain dict of floats.

    Convenience adapter so the Case-method primitives can work with the
    same numeric inputs regardless of whether they were pulled from
    :func:`orpheus.derivations._xs_library.get_mixture` or built
    manually for testing. Fails loudly if the mixture is not 1-group.
    """
    if mixture.ng != 1:
        raise ValueError(
            f"_xs_1g_dict requires a 1-group mixture, got ng={mixture.ng}. "
            "Multigroup Case-method reference is a Phase 2.1b+ follow-up."
        )
    sig_s = float(mixture.SigS[0][0, 0])
    return dict(
        sig_t=float(mixture.SigT[0]),
        sig_s=sig_s,
        nu_sig_f=float(mixture.SigP[0]),
    )


def derive_sn_heterogeneous_continuous(
    mat_A_name: str = "A",
    mat_B_name: str = "B",
    H_A: float = 0.5,
    H_B: float = 0.5,
    n_ordinates: int = 8,
) -> ContinuousReferenceSolution:
    r"""Phase 2.1b Case singular-eigenfunction reference for a 2-region
    reflective slab.

    Produces a :class:`~orpheus.derivations.ContinuousReferenceSolution`
    with ``operator_form="differential-sn"``,
    ``is_eigenvalue=True``, and a callable ``phi(x, g)`` from the
    back-substituted real-basis mode expansion. The reference is the
    exact discrete-:math:`S_N` eigenvalue at Gauss–Legendre order
    ``n_ordinates``, and should be consumed only by SN tests at the
    same quadrature order — the Case method is **not** the
    continuous-angle limit.

    Parameters
    ----------
    mat_A_name, mat_B_name : str
        Names in the :mod:`orpheus.derivations._xs_library` library.
        Default ``"A"`` (fuel) and ``"B"`` (moderator) reproduce the
        Phase 2.1b diagnostic configuration used to catch ERR-025.
    H_A, H_B : float
        Region thicknesses in cm. Default ``0.5 + 0.5`` matches the
        ERR-025 diagnostic problem.
    n_ordinates : int
        Gauss–Legendre quadrature order :math:`N` (number of ordinates).
        Default ``8`` matches ``solve_sn``'s default. This MUST match
        the quadrature order of any consumer test.

    Notes
    -----
    See :doc:`/theory/discrete_ordinates` "Case singular-eigenfunction
    method" section for the full derivation. The investigation
    history section of that page records the dead ends encountered
    during the Phase 2.1b prototype work — in particular the
    disagreement between an earlier Case implementation and
    ``solve_sn`` that eventually diagnosed ERR-025.
    """
    if n_ordinates % 2 != 0:
        raise ValueError(
            f"n_ordinates must be even for Gauss-Legendre, got {n_ordinates}"
        )
    mat_A = get_mixture(mat_A_name, "1g")
    mat_B = get_mixture(mat_B_name, "1g")
    xs_A = _xs_1g_dict(mat_A)
    xs_B = _xs_1g_dict(mat_B)

    mu, weights = np.polynomial.legendre.leggauss(n_ordinates)
    N = n_ordinates

    k_val = _solve_2region_reflective_sn_eigenvalue(
        xs_A, xs_B, H_A, H_B, mu, weights,
    )

    # Extract the 2N-d coefficient vector once at the fundamental root.
    C = _assemble_matching_matrix_sn_1g(
        xs_A, xs_B, H_A, H_B, k_val, mu, weights,
    )
    null_vec = _extract_real_null_vector_sn(C)
    c_A = null_vec[:N]
    c_B = null_vec[N:]

    # Sign-normalise so phi is non-negative at the midpoint, and
    # amplitude-normalise so max|phi| on a dense sweep is exactly 1.
    _dense_x = np.linspace(0.0, H_A + H_B, 2001)
    _dense_phi = np.zeros(_dense_x.size)
    for i, xi in enumerate(_dense_x):
        _, _dense_phi[i] = _evaluate_flux_at(
            xs_A, xs_B, H_A, H_B, k_val, float(xi),
            c_A, c_B, mu, weights,
        )
    _peak_idx = int(np.argmax(np.abs(_dense_phi)))
    if _dense_phi[_peak_idx] < 0:
        c_A = -c_A
        c_B = -c_B
        _dense_phi = -_dense_phi
    _max_abs = float(np.max(np.abs(_dense_phi)))
    if _max_abs == 0:
        raise RuntimeError(
            "Back-substituted flux is identically zero — "
            "null-vector extraction failed for the SN 2-region reference."
        )
    c_A = c_A / _max_abs
    c_B = c_B / _max_abs

    def phi(x: np.ndarray, g: int = 0) -> np.ndarray:
        if g != 0:
            raise IndexError(
                f"1-group reference has only g=0, got g={g}"
            )
        x_arr = np.asarray(x, dtype=float)
        orig_shape = x_arr.shape
        x_flat = x_arr.ravel()
        out = np.zeros(x_flat.size)
        for i, xi in enumerate(x_flat):
            _, out[i] = _evaluate_flux_at(
                xs_A, xs_B, H_A, H_B, k_val, float(xi),
                c_A, c_B, mu, weights,
            )
        return out.reshape(orig_shape)

    return ContinuousReferenceSolution(
        name=f"sn_slab_1eg_2rg_S{N}",
        problem=ProblemSpec(
            materials={0: mat_A, 1: mat_B},
            geometry_type="slab",
            geometry_params={
                "fuel_height": H_A,
                "refl_height": H_B,
                "length": H_A + H_B,
                "n_ordinates": N,
            },
            boundary_conditions={"left": "reflective", "right": "reflective"},
            external_source=None,
            is_eigenvalue=True,
            n_groups=1,
        ),
        operator_form="differential-sn",
        phi=phi,
        provenance=Provenance(
            citation=(
                "Case & Zweifel 1967, *Linear Transport Theory*, Ch. 6 "
                "(singular-eigenfunction expansion for multi-region "
                "slabs); Gülderen & Türeci 2023, Cumhuriyet Sci. J. "
                "44(4), 775–784, DOI 10.17776/csj.1240161 (canonical "
                "reflective-BC matching matrix derivation); Garis & "
                "Sjöstrand 1990, Ann. Nucl. Energy 17, DOI "
                "10.1016/0306-4549(90)90029-d (one-speed two-medium "
                "slab lattice, mathematically equivalent to reflective "
                "2-region slab by unit-cell symmetry)."
            ),
            derivation_notes=(
                "Real-basis mode decomposition of the discrete-S_N slope "
                "matrix S(k) = M^{-1}(K(k) - Σ_t I) for a 1-group 2-region "
                "reflective slab. Per region, diagonalise S(k) via "
                "numpy.linalg.eig; for each real eigenvalue keep one "
                "bounded exponential mode anchored at the nearer edge "
                "(right edge for λ ≥ 0, left edge for λ < 0); for each "
                "complex-conjugate pair keep two real modes built from "
                "the cos/sin/Re-Im combination times an anchored "
                "exponential envelope. This bounded real basis is the "
                "Phase 1.2 diffusion fix ported verbatim to the "
                "first-order SN operator.\n\n"
                "The 2N×2N matching matrix C(k) enforces N/2 reflective "
                "BC equations at each outer edge (ψ_A(0, +μ_n) = "
                "ψ_A(0, -μ_n) and likewise at x = L) plus N interface "
                "continuity equations at x = H_A (ψ_A(H_A, μ_n) = "
                "ψ_B(H_A, μ_n) for all n). Gauss-Legendre ordinate pairing "
                "via the N-1-n reflection index. det(C(k)) = 0 is solved "
                "by coarse scan plus scipy.optimize.brentq to xtol=1e-14.\n\n"
                "Spurious sign changes from np.linalg.eig's non-continuous "
                "eigenvalue ordering across k crossings are filtered by "
                "physical validation: every candidate root rebuilds the "
                "null vector via SVD of C(k), then explicitly reconstructs "
                "the reflective-BC residuals at x=0 and x=L and the "
                "interface continuity residual at x=H_A. Only candidates "
                "with all three residuals below 1e-6 (relative to the "
                "peak angular flux) are accepted. The fundamental is the "
                "largest validated root.\n\n"
                "Back-substitution gives the continuous scalar flux phi(x) "
                "via psi(x) = basis(x) @ c_region, phi(x) = Σ w_n psi_n(x). "
                "All basis modes are bounded by O(1) on their region, so "
                "evaluation is stable to machine precision.\n\n"
                "This is the exact discrete-S_N eigenvalue at the "
                "specified quadrature order, NOT the continuous-angle "
                "limit. Consumer tests must use the same GL order. "
                "The reference does NOT consume or cross-check against "
                "solve_sn — it is an independent mathematical construction "
                "from the transport equation, standing on its own merits "
                "in line with Cardinal Rule 2 of the verification "
                "campaign."
            ),
            sympy_expression=None,
            precision_digits=None,
        ),
        k_eff=k_val,
        psi=None,
        equation_labels=(
            "transport-cartesian",
            "reflective-bc",
            "one-group-kinf",
            "dd-cartesian-1d",
            "dd-recurrence",
            "sn-case-per-ordinate",
            "sn-case-slope-matrix",
            "sn-case-spatial-modes",
            "sn-case-real-basis",
            "sn-case-matching-matrix",
            "sn-case-physical-validation",
            "sn-case-back-substitution",
        ),
        vv_level="L1",
        description=(
            f"1-group Case singular-eigenfunction reference, 2-region "
            f"reflective slab ({mat_A_name}+{mat_B_name}, H_A={H_A} + "
            f"H_B={H_B} cm), S{N} Gauss-Legendre. Phase 2.1b — "
            f"k_eff = {k_val:.10f}."
        ),
        tolerance="< 1e-8 at matching quadrature order",
    )


def continuous_cases() -> list[ContinuousReferenceSolution]:
    """Return all SN Phase-2.1b continuous reference solutions.

    Currently ships the 1-group 2-region reflective slab at :math:`S_8`
    (the Phase 2.1b diagnostic configuration that caught ERR-025).
    Additional ordinate orders and multigroup extensions can be added
    here as follow-ups.
    """
    return [derive_sn_heterogeneous_continuous()]
