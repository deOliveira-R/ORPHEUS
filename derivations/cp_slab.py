"""SymPy derivations for slab collision probability eigenvalues.

Derives analytical k_inf for 2-region slab geometry using the E₃
exponential integral kernel. The CP matrix structure is expressed
symbolically; E₃ values are computed numerically.
"""

from __future__ import annotations

import numpy as np
import sympy as sp
from scipy.sparse import csr_matrix

from data.macro_xs.mixture import Mixture

from ._kernels import e3
from ._types import VerificationCase


def _make_mixture(
    sig_t: np.ndarray,
    sig_c: np.ndarray,
    sig_f: np.ndarray,
    nu: np.ndarray,
    chi: np.ndarray,
    sig_s: np.ndarray,
) -> Mixture:
    """Build a Mixture from N-group arrays."""
    ng = len(sig_t)
    eg = np.logspace(7, -3, ng + 1)
    return Mixture(
        SigC=sig_c.copy(), SigL=np.zeros(ng),
        SigF=sig_f.copy(), SigP=(nu * sig_f).copy(),
        SigT=sig_t.copy(), SigS=[csr_matrix(sig_s)],
        Sig2=csr_matrix((ng, ng)), chi=chi.copy(), eg=eg.copy(),
    )


# ═══════════════════════════════════════════════════════════════════════
# Slab CP matrix from E₃ kernel (N regions, N_g groups)
# ═══════════════════════════════════════════════════════════════════════

def _slab_cp_matrix(
    sig_t_all: np.ndarray,
    t_arr: np.ndarray,
) -> np.ndarray:
    """Compute the infinite-lattice CP matrix for a slab.

    Uses the E₃ second-difference formula with white boundary condition.

    Parameters
    ----------
    sig_t_all : (N_reg, ng) — total XS per region and group
    t_arr : (N_reg,) — region thicknesses

    Returns
    -------
    P_inf : (N_reg, N_reg, ng) — collision probability matrix
    """
    N_reg = len(t_arr)
    ng = sig_t_all.shape[1]
    P_inf_g = np.zeros((N_reg, N_reg, ng))

    for g in range(ng):
        sig_t_g = sig_t_all[:, g]
        tau = sig_t_g * t_arr  # optical thicknesses

        # Boundary positions in optical coordinates
        bnd_pos = np.zeros(N_reg + 1)
        for i in range(N_reg):
            bnd_pos[i + 1] = bnd_pos[i] + tau[i]

        # Reduced collision probability via E₃ second differences
        rcp = np.zeros((N_reg, N_reg))
        for i in range(N_reg):
            # Self-collision (diagonal)
            rcp[i, i] += 0.5 * sig_t_g[i] * (
                2 * t_arr[i] - (2.0 / sig_t_g[i]) * (0.5 - e3(tau[i]))
            )

            for j in range(N_reg):
                tau_i, tau_j = tau[i], tau[j]

                # Direct path (different regions only)
                if j > i:
                    gap_d = max(bnd_pos[j] - bnd_pos[i + 1], 0.0)
                elif j < i:
                    gap_d = max(bnd_pos[i] - bnd_pos[j + 1], 0.0)
                else:
                    gap_d = None

                dd = 0.0
                if gap_d is not None:
                    dd = (e3(gap_d) - e3(gap_d + tau_i)
                          - e3(gap_d + tau_j) + e3(gap_d + tau_i + tau_j))

                # Reflected path (always present)
                gap_c = bnd_pos[i] + bnd_pos[j]
                dc = (e3(gap_c) - e3(gap_c + tau_i)
                      - e3(gap_c + tau_j) + e3(gap_c + tau_i + tau_j))

                rcp[i, j] += 0.5 * (dd + dc)

        # Normalise to get P_cell
        P_cell = np.zeros((N_reg, N_reg))
        for i in range(N_reg):
            P_cell[i, :] = rcp[i, :] / (sig_t_g[i] * t_arr[i])

        # White-BC correction: P_inf = P_cell + P_out * P_in^T / (1 - P_inout)
        P_out = np.maximum(1.0 - P_cell.sum(axis=1), 0.0)
        P_in = sig_t_g * t_arr * P_out
        P_inout = max(1.0 - P_in.sum(), 0.0)
        P_inf_g[:, :, g] = P_cell + np.outer(P_out, P_in) / (1.0 - P_inout)

    return P_inf_g


def _kinf_from_cp(
    P_inf_g: np.ndarray,
    sig_t_all: np.ndarray,
    V_arr: np.ndarray,
    sig_s_mats: list[np.ndarray],
    nu_sig_f_mats: list[np.ndarray],
    chi_mats: list[np.ndarray],
) -> float:
    """Solve the CP eigenvalue problem for k_inf.

    The generalised eigenvalue problem is:

    .. math::
        (\\text{diag}(\\Sigma_t V) - P^T \\Sigma_s V) \\phi
        = \\frac{1}{k} P^T \\chi \\nu\\Sigma_f V \\phi
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
                    A_mat[row, col] -= Pji * V_arr[j_reg] * sig_s_mats[j_reg][gp, g]
                    B_mat[row, col] += (
                        Pji * V_arr[j_reg]
                        * chi_mats[j_reg][g]
                        * nu_sig_f_mats[j_reg][gp]
                    )

    M = np.linalg.solve(A_mat, B_mat)
    return float(np.max(np.real(np.linalg.eigvals(M))))


def _symbolic_latex_cp_eigenvalue(ng: int, n_reg: int) -> str:
    """Generate LaTeX for the CP eigenvalue problem."""
    dim = n_reg * ng
    return (
        r"The slab CP eigenvalue problem with white boundary conditions:"
        "\n\n"
        r".. math::" "\n"
        r"   \left(\text{diag}(\Sigma_{t,i} V_i) "
        r"- \mathbf{P}^T \text{diag}(V_j \Sigma_{s,j})\right) \phi "
        r"= \frac{1}{k} \mathbf{P}^T \text{diag}(V_j \chi_j \nu\Sigma_{f,j}) \phi"
        "\n\n"
        rf"where :math:`\mathbf{{P}} \in \mathbb{{R}}^{{{n_reg} \times {n_reg}}}` "
        r"is the collision probability matrix per energy group, computed from "
        r"the :math:`E_3` exponential integral kernel with the second-difference "
        r"formula."
        "\n\n"
        rf"The full eigenvalue problem has dimension {dim} "
        rf"({n_reg} regions :math:`\times` {ng} groups)."
    )


# ═══════════════════════════════════════════════════════════════════════
# Cross-section sets for verification cases
# ═══════════════════════════════════════════════════════════════════════

# Region A: fissile
_XS_A = dict(
    sig_t_1g=np.array([1.0]), sig_c_1g=np.array([0.2]),
    sig_f_1g=np.array([0.3]), nu_1g=np.array([2.5]),
    chi_1g=np.array([1.0]), sig_s_1g=np.array([[0.5]]),

    sig_t_2g=np.array([0.50, 1.00]), sig_c_2g=np.array([0.01, 0.02]),
    sig_f_2g=np.array([0.01, 0.08]), nu_2g=np.array([2.50, 2.50]),
    chi_2g=np.array([1.00, 0.00]),
    sig_s_2g=np.array([[0.38, 0.10], [0.00, 0.90]]),
)

# Region B: non-fissile scatterer
_XS_B = dict(
    sig_t_1g=np.array([2.0]), sig_c_1g=np.array([0.1]),
    sig_f_1g=np.array([0.0]), nu_1g=np.array([0.0]),
    chi_1g=np.array([1.0]), sig_s_1g=np.array([[1.9]]),

    sig_t_2g=np.array([0.60, 2.00]), sig_c_2g=np.array([0.02, 0.05]),
    sig_f_2g=np.array([0.00, 0.00]), nu_2g=np.array([0.00, 0.00]),
    chi_2g=np.array([1.00, 0.00]),
    sig_s_2g=np.array([[0.40, 0.18], [0.00, 1.95]]),
)


def _build_case(
    suffix: str,
    xs_a: dict,
    xs_b: dict,
    t_a: float = 0.5,
    t_b: float = 0.5,
) -> VerificationCase:
    """Build a 2-region slab verification case for given group count."""
    mix_a = _make_mixture(
        xs_a[f"sig_t_{suffix}"], xs_a[f"sig_c_{suffix}"],
        xs_a[f"sig_f_{suffix}"], xs_a[f"nu_{suffix}"],
        xs_a[f"chi_{suffix}"], xs_a[f"sig_s_{suffix}"],
    )
    mix_b = _make_mixture(
        xs_b[f"sig_t_{suffix}"], xs_b[f"sig_c_{suffix}"],
        xs_b[f"sig_f_{suffix}"], xs_b[f"nu_{suffix}"],
        xs_b[f"chi_{suffix}"], xs_b[f"sig_s_{suffix}"],
    )

    ng = len(xs_a[f"sig_t_{suffix}"])
    sig_t_all = np.vstack([xs_a[f"sig_t_{suffix}"], xs_b[f"sig_t_{suffix}"]])
    t_arr = np.array([t_a, t_b])

    P_inf_g = _slab_cp_matrix(sig_t_all, t_arr)

    k_inf = _kinf_from_cp(
        P_inf_g=P_inf_g,
        sig_t_all=sig_t_all,
        V_arr=t_arr,
        sig_s_mats=[xs_a[f"sig_s_{suffix}"], xs_b[f"sig_s_{suffix}"]],
        nu_sig_f_mats=[
            xs_a[f"nu_{suffix}"] * xs_a[f"sig_f_{suffix}"],
            xs_b[f"nu_{suffix}"] * xs_b[f"sig_f_{suffix}"],
        ],
        chi_mats=[xs_a[f"chi_{suffix}"], xs_b[f"chi_{suffix}"]],
    )

    geom_params = dict(
        n_fuel=1, n_clad=0, n_cool=1,
        fuel_half=t_a, clad_thick=0.0, cool_thick=t_b,
    )

    name = f"cp_slab_{ng}eg_2rg"
    latex = _symbolic_latex_cp_eigenvalue(ng, 2) + (
        "\n\n"
        rf"For the {ng}-group, 2-region slab (t_A = {t_a}, t_B = {t_b}):"
        "\n\n"
        r".. math::" "\n"
        rf"   k_\infty = {k_inf:.10f}"
    )

    return VerificationCase(
        name=name,
        k_inf=k_inf,
        method="cp",
        geometry="slab",
        n_groups=ng,
        n_regions=2,
        materials={2: mix_a, 0: mix_b},
        geom_params=geom_params,
        latex=latex,
        description=f"{ng}-group 2-region slab CP (E₃ kernel, white BC)",
        tolerance="< 1e-6",
    )


def derive_1g_slab() -> VerificationCase:
    """1-group, 2-region slab CP verification case."""
    return _build_case("1g", _XS_A, _XS_B)


def derive_2g_slab() -> VerificationCase:
    """2-group, 2-region slab CP verification case."""
    return _build_case("2g", _XS_A, _XS_B)


def all_cases() -> list[VerificationCase]:
    """Return all slab CP verification cases."""
    return [derive_1g_slab(), derive_2g_slab()]
