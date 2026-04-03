"""SymPy derivations for cylindrical collision probability eigenvalues.

Derives analytical k_inf for 2-region Wigner-Seitz cylindrical geometry
using the Ki₃/Ki₄ Bickley-Naylor kernel. The CP matrix structure is
expressed symbolically; Ki₄ values are computed numerically.
"""

from __future__ import annotations

import numpy as np
import sympy as sp
from scipy.sparse import csr_matrix

from data.macro_xs.mixture import Mixture

from ._kernels import bickley_tables
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


def _chord_half_lengths(radii: np.ndarray, y_pts: np.ndarray) -> np.ndarray:
    """Half-chord lengths l_k(y) for each annular region.

    Parameters
    ----------
    radii : (N,) outer radius of each annular region
    y_pts : (n_y,) perpendicular distance quadrature points

    Returns
    -------
    chords : (N, n_y) half-chord length per region per y-point
    """
    N = len(radii)
    chords = np.zeros((N, len(y_pts)))
    r_inner = np.zeros(N)
    r_inner[1:] = radii[:-1]
    y2 = y_pts**2

    for k in range(N):
        r_out, r_in = radii[k], r_inner[k]
        outer = np.sqrt(np.maximum(r_out**2 - y2, 0.0))
        if r_in > 0:
            inner = np.sqrt(np.maximum(r_in**2 - y2, 0.0))
            mask = y_pts < r_in
            chords[k, mask] = outer[mask] - inner[mask]
        mask_p = (y_pts >= r_in) & (y_pts < r_out)
        chords[k, mask_p] = outer[mask_p]

    return chords


def _cylinder_cp_matrix(
    sig_t_all: np.ndarray,
    radii: np.ndarray,
    volumes: np.ndarray,
    r_cell: float,
    n_quad_y: int = 64,
) -> np.ndarray:
    """Compute the infinite-lattice CP matrix for a cylindrical cell.

    Uses the Ki₄ second-difference formula with y-quadrature and
    white boundary condition.

    Parameters
    ----------
    sig_t_all : (N_reg, ng)
    radii : (N_reg,) outer radii of annular regions
    volumes : (N_reg,) annular areas (pi*(r_k^2 - r_{k-1}^2))
    r_cell : equivalent cell radius
    n_quad_y : number of GL quadrature points per segment

    Returns
    -------
    P_inf : (N_reg, N_reg, ng)
    """
    N_reg = len(radii)
    ng = sig_t_all.shape[1]

    # Build Ki₄ lookup table
    tables = bickley_tables()

    # Build y-quadrature with breakpoints at each radius
    gl_pts, gl_wts = np.polynomial.legendre.leggauss(n_quad_y)
    breakpoints = np.concatenate(([0.0], radii))
    y_all, w_all = [], []
    for seg in range(len(breakpoints) - 1):
        a, b = breakpoints[seg], breakpoints[seg + 1]
        y_all.append(0.5 * (b - a) * gl_pts + 0.5 * (b + a))
        w_all.append(0.5 * (b - a) * gl_wts)
    y_pts = np.concatenate(y_all)
    y_wts = np.concatenate(w_all)

    chords = _chord_half_lengths(radii, y_pts)
    n_y = len(y_pts)

    ki4_0 = tables.ki4_vec(np.zeros(n_y))

    P_inf_g = np.empty((N_reg, N_reg, ng))

    for g in range(ng):
        sig_t_g = sig_t_all[:, g]
        tau = sig_t_g[:, None] * chords  # (N_reg, n_y)

        # Boundary positions in optical coordinates
        bnd_pos = np.zeros((N_reg + 1, n_y))
        for k in range(N_reg):
            bnd_pos[k + 1, :] = bnd_pos[k, :] + tau[k, :]

        rcp = np.zeros((N_reg, N_reg))

        for i in range(N_reg):
            tau_i = tau[i, :]
            sti = sig_t_g[i]
            if sti == 0:
                continue

            # Self-same collision
            self_same = 2.0 * chords[i, :] - (2.0 / sti) * (
                ki4_0 - tables.ki4_vec(tau_i)
            )
            rcp[i, i] += 2.0 * sti * np.dot(y_wts, self_same)

            for j in range(N_reg):
                tau_j = tau[j, :]

                # Direct path
                if j > i:
                    gap_d = np.maximum(bnd_pos[j, :] - bnd_pos[i + 1, :], 0.0)
                elif j < i:
                    gap_d = np.maximum(bnd_pos[i, :] - bnd_pos[j + 1, :], 0.0)
                else:
                    gap_d = None

                if gap_d is not None:
                    dd = (tables.ki4_vec(gap_d)
                          - tables.ki4_vec(gap_d + tau_i)
                          - tables.ki4_vec(gap_d + tau_j)
                          + tables.ki4_vec(gap_d + tau_i + tau_j))
                else:
                    dd = np.zeros(n_y)

                # Reflected path
                gap_c = bnd_pos[i, :] + bnd_pos[j, :]
                dc = (tables.ki4_vec(gap_c)
                      - tables.ki4_vec(gap_c + tau_i)
                      - tables.ki4_vec(gap_c + tau_j)
                      + tables.ki4_vec(gap_c + tau_i + tau_j))

                rcp[i, j] += 2.0 * np.dot(y_wts, dd + dc)

        # Normalise to P_cell
        P_cell = np.zeros((N_reg, N_reg))
        for i in range(N_reg):
            if sig_t_g[i] * volumes[i] > 0:
                P_cell[i, :] = rcp[i, :] / (sig_t_g[i] * volumes[i])

        # White-BC correction
        P_out = np.maximum(1.0 - P_cell.sum(axis=1), 0.0)
        S_cell = 2.0 * np.pi * r_cell
        P_in = sig_t_g * volumes * P_out / S_cell
        P_inout = max(1.0 - P_in.sum(), 0.0)

        P_inf = P_cell.copy()
        if P_inout < 1.0:
            P_inf += np.outer(P_out, P_in) / (1.0 - P_inout)

        P_inf_g[:, :, g] = P_inf

    return P_inf_g


def _kinf_from_cp(
    P_inf_g: np.ndarray,
    sig_t_all: np.ndarray,
    V_arr: np.ndarray,
    sig_s_mats: list[np.ndarray],
    nu_sig_f_mats: list[np.ndarray],
    chi_mats: list[np.ndarray],
) -> float:
    """Solve the CP eigenvalue problem for k_inf."""
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


# ═══════════════════════════════════════════════════════════════════════
# Cross-section sets (same as cp_slab — abstract regions A and B)
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


def _build_case(suffix: str, xs_a: dict, xs_b: dict) -> VerificationCase:
    """Build a 2-region cylindrical verification case."""
    r_fuel, r_cell = 0.5, 1.0

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

    # Cylindrical geometry
    radii = np.array([r_fuel, r_cell])
    volumes = np.array([
        np.pi * r_fuel**2,
        np.pi * (r_cell**2 - r_fuel**2),
    ])

    P_inf_g = _cylinder_cp_matrix(sig_t_all, radii, volumes, r_cell)

    k_inf = _kinf_from_cp(
        P_inf_g=P_inf_g,
        sig_t_all=sig_t_all,
        V_arr=volumes,
        sig_s_mats=[xs_a[f"sig_s_{suffix}"], xs_b[f"sig_s_{suffix}"]],
        nu_sig_f_mats=[
            xs_a[f"nu_{suffix}"] * xs_a[f"sig_f_{suffix}"],
            xs_b[f"nu_{suffix}"] * xs_b[f"sig_f_{suffix}"],
        ],
        chi_mats=[xs_a[f"chi_{suffix}"], xs_b[f"chi_{suffix}"]],
    )

    # Geometry params matching solver convention
    pitch = r_cell * np.sqrt(np.pi)
    geom_params = dict(
        n_fuel=1, n_clad=0, n_cool=1,
        r_fuel=r_fuel, r_clad=r_fuel,
        pitch=pitch,
    )

    name = f"cp_cyl1D_{ng}eg_2rg"
    latex = (
        r"The cylindrical CP eigenvalue problem with white boundary conditions "
        r"uses the Ki₄ Bickley-Naylor kernel:"
        "\n\n"
        r".. math::" "\n"
        r"   \text{Ki}_3(x) = \int_0^{\pi/2} "
        r"\exp\!\left(-\frac{x}{\sin\theta}\right) \sin\theta\, d\theta"
        "\n\n"
        r".. math::" "\n"
        r"   \text{Ki}_4(x) = \int_x^\infty \text{Ki}_3(t)\, dt"
        "\n\n"
        r"The reduced collision probability uses the Ki₄ second-difference "
        r"formula integrated over chord perpendicular distance :math:`y`."
        "\n\n"
        rf"For the {ng}-group, 2-region cylinder "
        rf"(:math:`r_\text{{fuel}} = {r_fuel}`, "
        rf":math:`r_\text{{cell}} = {r_cell}`):"
        "\n\n"
        r".. math::" "\n"
        rf"   k_\infty = {k_inf:.10f}"
    )

    return VerificationCase(
        name=name,
        k_inf=k_inf,
        method="cp",
        geometry="cyl1D",
        n_groups=ng,
        n_regions=2,
        materials={2: mix_a, 0: mix_b},
        geom_params=geom_params,
        latex=latex,
        description=f"{ng}-group 2-region cylindrical CP (Ki₄ kernel, white BC)",
        tolerance="< 1e-5",
    )


def derive_1g_cylinder() -> VerificationCase:
    """1-group, 2-region Wigner-Seitz cylinder verification case."""
    return _build_case("1g", _XS_A, _XS_B)


def derive_2g_cylinder() -> VerificationCase:
    """2-group, 2-region Wigner-Seitz cylinder verification case."""
    return _build_case("2g", _XS_A, _XS_B)


def all_cases() -> list[VerificationCase]:
    """Return all cylindrical CP verification cases."""
    return [derive_1g_cylinder(), derive_2g_cylinder()]
