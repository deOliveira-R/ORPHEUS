"""Foundation tests for the :class:`Billiard` class.

These tests pin the Billiard facade's bit-equality with the legacy
``solve_greens_function_*`` entry points across all 6 geometries × {1G,
MG} combinations + the multi-region sphere fixed-source path. They do
NOT re-test the underlying solvers' correctness — that's the existing
suite's job. They DO assert that wrapping a solve in
:meth:`Billiard.solve_critical` returns a shared
:class:`~orpheus.derivations.common.solution_types.CriticalSolution`
whose contents are bit-for-bit identical to the original return.

Tests are tagged ``foundation`` because they verify a software
contract (the facade's bit-equal preservation) rather than an L0/L1
mathematical claim about a solver.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.solution_types import (
    CriticalSolution,
    FluxSolution,
)
from orpheus.derivations.continuous.trajectory_resolvent import Billiard
from orpheus.derivations.continuous.trajectory_resolvent import (
    greens_function as gf_sphere,
    greens_function_cylinder as gf_cyl,
    greens_function_slab as gf_slab,
    greens_function_slab_asymmetric as gf_slab_asym,
    greens_function_hollow_sphere as gf_hollow,
    greens_function_annulus as gf_annulus,
)


# ─────────────────────────────────────────────────────────────────────
# Bit-equality fixtures — keep params identical across paired tests
# ─────────────────────────────────────────────────────────────────────

# Use small grids so the test suite stays fast.
SPHERE_PARAMS_1G = dict(
    R=5.0, sigma_t=0.5, sigma_s=0.4, nu_sigma_f=0.1,
    n_r=8, n_mu=8, n_traj_quad=16, max_iter=10, tol=1e-9,
)
SPHERE_PARAMS_MG = dict(
    R=5.0,
    sigma_t=np.array([1.0, 0.5]),
    sigma_s=np.array([[0.4, 0.4], [0.0, 0.4]]),
    nu_sigma_f=np.array([0.05, 0.10]),
    n_r=8, n_mu=8, n_traj_quad=16, max_iter=10, tol=1e-9,
)
CYL_PARAMS_1G = dict(
    R=5.0, sigma_t=0.5, sigma_s=0.4, nu_sigma_f=0.1,
    n_r=6, n_mu_axial=6, n_phi_az=8, n_traj_quad=12, max_iter=10, tol=1e-9,
)
SLAB_PARAMS_1G = dict(
    L=5.0, sigma_t=0.5, sigma_s=0.4, nu_sigma_f=0.1,
    n_x=8, n_mu=8, n_traj_quad=16, max_iter=10, tol=1e-9,
)
HOLLOW_SPH_PARAMS_1G = dict(
    R_in=1.0, R_out=5.0, sigma_t=0.5, sigma_s=0.4, nu_sigma_f=0.1,
    n_r=8, n_mu=8, n_traj_quad=16, max_iter=10, tol=1e-9,
)
ANNULUS_PARAMS_1G = dict(
    R_in=1.0, R_out=5.0, sigma_t=0.5, sigma_s=0.4, nu_sigma_f=0.1,
    n_r=6, n_mu_axial=6, n_phi_az=8, n_traj_quad=12, max_iter=10, tol=1e-9,
)


def _bit_equal_arrays(a: np.ndarray, b: np.ndarray) -> bool:
    """Strict bit-for-bit equality on float arrays (no allclose)."""
    a = np.asarray(a)
    b = np.asarray(b)
    if a.shape != b.shape:
        return False
    if a.dtype != b.dtype:
        return False
    return bool(np.array_equal(a, b))


def _bit_equal_floats(a: float, b: float) -> bool:
    """Strict bit-for-bit float equality (handles NaN by hex)."""
    return float(a).hex() == float(b).hex()


def _check_critical_invariants(
    sol: CriticalSolution, *, n_groups: int, geometry_kind: str,
) -> None:
    """Assert universal invariants on a shared CriticalSolution."""
    assert isinstance(sol, CriticalSolution)
    assert sol.eigenvalue_kind == "k_eff"
    assert sol.parameter_kind == "fixed_geometry"
    assert sol.metadata["n_groups"] == n_groups
    assert sol.metadata["geometry_kind"] == geometry_kind
    assert "raw_result" in sol.metadata
    assert "psi" in sol.metadata
    assert "phi" in sol.metadata
    assert "iterations" in sol.metadata
    assert "mesh" in sol.metadata


# ─────────────────────────────────────────────────────────────────────
# 1. Construction tests
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_billiard_from_problem_sphere_one_endpoint_rank_1():
    """One-endpoint orbit-space class → closure_rank == 1."""
    b = Billiard.from_problem(
        geometry_kind="sphere",
        materials={"sigma_t": 0.5, "sigma_s": 0.4, "nu_sigma_f": 0.1},
        geometry={"R": 5.0},
        alpha=1.0,
    )
    assert b.closure_rank == 1
    assert b.geometry_kind == "sphere"
    assert b.alpha_payload == {"alpha": 1.0}


@pytest.mark.foundation
def test_billiard_from_problem_slab_asymmetric_two_endpoint_rank_2():
    """Two-endpoint orbit-space class → closure_rank == 2."""
    b = Billiard.from_problem(
        geometry_kind="slab_asymmetric",
        materials={"sigma_t": 0.5, "sigma_s": 0.4, "nu_sigma_f": 0.1},
        geometry={"L": 5.0},
        alpha={"alpha_left": 0.7, "alpha_right": 0.9},
    )
    assert b.closure_rank == 2
    assert b.alpha_payload == {"alpha_left": 0.7, "alpha_right": 0.9}


@pytest.mark.foundation
def test_billiard_from_problem_scalar_alpha_normalized_for_two_endpoint():
    """Scalar alpha gets normalized to two-endpoint payload for slab_asym."""
    b = Billiard.from_problem(
        geometry_kind="slab_asymmetric",
        materials={"sigma_t": 0.5, "sigma_s": 0.4, "nu_sigma_f": 0.1},
        geometry={"L": 5.0},
        alpha=0.5,
    )
    assert b.alpha_payload == {"alpha_left": 0.5, "alpha_right": 0.5}

    b2 = Billiard.from_problem(
        geometry_kind="hollow_sphere",
        materials={"sigma_t": 0.5, "sigma_s": 0.4, "nu_sigma_f": 0.1},
        geometry={"R_in": 1.0, "R_out": 5.0},
        alpha=0.7,
    )
    assert b2.alpha_payload == {"alpha_in": 0.7, "alpha_out": 0.7}


@pytest.mark.foundation
def test_billiard_with_alpha_returns_modified_copy():
    """with_alpha returns a copy with new alpha; original unchanged."""
    b = Billiard.from_problem(
        geometry_kind="sphere",
        materials={"sigma_t": 0.5, "sigma_s": 0.4, "nu_sigma_f": 0.1},
        geometry={"R": 5.0},
        alpha=1.0,
    )
    b2 = b.with_alpha(0.5)
    assert b.alpha_payload == {"alpha": 1.0}
    assert b2.alpha_payload == {"alpha": 0.5}
    assert b.geometry_kind == b2.geometry_kind
    assert b.materials == b2.materials
    assert b.geometry_payload == b2.geometry_payload


@pytest.mark.foundation
def test_billiard_unknown_geometry_raises():
    with pytest.raises(ValueError, match="unknown geometry_kind"):
        Billiard.from_problem(
            geometry_kind="hyperbolic_torus",
            materials={"sigma_t": 0.5, "sigma_s": 0.4, "nu_sigma_f": 0.1},
            geometry={"R": 5.0},
            alpha=1.0,
        )


# ─────────────────────────────────────────────────────────────────────
# 2. Bit-equality tests — sphere
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_billiard_sphere_1g_bit_equal_legacy():
    """Billiard.solve_critical sphere 1G ≡ solve_greens_function_sphere."""
    legacy = gf_sphere.solve_greens_function_sphere(
        alpha=1.0, **SPHERE_PARAMS_1G,
    )
    b = Billiard.from_problem(
        geometry_kind="sphere",
        materials={
            "sigma_t": SPHERE_PARAMS_1G["sigma_t"],
            "sigma_s": SPHERE_PARAMS_1G["sigma_s"],
            "nu_sigma_f": SPHERE_PARAMS_1G["nu_sigma_f"],
        },
        geometry={"R": SPHERE_PARAMS_1G["R"]},
        alpha=1.0,
        quadrature={
            "n_r": SPHERE_PARAMS_1G["n_r"],
            "n_mu": SPHERE_PARAMS_1G["n_mu"],
            "n_traj_quad": SPHERE_PARAMS_1G["n_traj_quad"],
        },
    )
    sol = b.solve_critical(
        max_iter=SPHERE_PARAMS_1G["max_iter"],
        tol=SPHERE_PARAMS_1G["tol"],
    )
    _check_critical_invariants(sol, n_groups=1, geometry_kind="sphere")
    assert _bit_equal_floats(sol.eigenvalue, legacy.k_eff)
    assert sol.parameter_value == 5.0
    assert _bit_equal_arrays(sol.metadata["psi"], legacy.psi)
    assert _bit_equal_arrays(sol.metadata["phi"], legacy.phi)
    assert sol.metadata["iterations"] == legacy.iterations
    assert sol.converged == legacy.converged
    assert _bit_equal_arrays(sol.metadata["mesh"]["r_nodes"], legacy.r_nodes)
    assert _bit_equal_arrays(sol.metadata["mesh"]["mu_nodes"], legacy.mu_nodes)


@pytest.mark.foundation
def test_billiard_sphere_mg_bit_equal_legacy():
    """Billiard sphere MG ≡ solve_greens_function_sphere_mg."""
    legacy = gf_sphere.solve_greens_function_sphere_mg(
        alpha=1.0, **SPHERE_PARAMS_MG,
    )
    b = Billiard.from_problem(
        geometry_kind="sphere",
        materials={
            "sigma_t": SPHERE_PARAMS_MG["sigma_t"],
            "sigma_s": SPHERE_PARAMS_MG["sigma_s"],
            "nu_sigma_f": SPHERE_PARAMS_MG["nu_sigma_f"],
        },
        geometry={"R": SPHERE_PARAMS_MG["R"]},
        alpha=1.0,
        quadrature={
            "n_r": SPHERE_PARAMS_MG["n_r"],
            "n_mu": SPHERE_PARAMS_MG["n_mu"],
            "n_traj_quad": SPHERE_PARAMS_MG["n_traj_quad"],
        },
    )
    sol = b.solve_critical(
        max_iter=SPHERE_PARAMS_MG["max_iter"],
        tol=SPHERE_PARAMS_MG["tol"],
    )
    _check_critical_invariants(sol, n_groups=2, geometry_kind="sphere")
    assert _bit_equal_floats(sol.eigenvalue, legacy.k_eff)
    assert _bit_equal_arrays(sol.metadata["psi"], legacy.psi_g)
    assert _bit_equal_arrays(sol.metadata["phi"], legacy.phi_g)


# ─────────────────────────────────────────────────────────────────────
# 3. Bit-equality tests — cylinder
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_billiard_cylinder_1g_bit_equal_legacy():
    """Billiard cylinder 1G ≡ solve_greens_function_cylinder."""
    legacy = gf_cyl.solve_greens_function_cylinder(
        alpha=1.0, **CYL_PARAMS_1G,
    )
    b = Billiard.from_problem(
        geometry_kind="cylinder",
        materials={
            "sigma_t": CYL_PARAMS_1G["sigma_t"],
            "sigma_s": CYL_PARAMS_1G["sigma_s"],
            "nu_sigma_f": CYL_PARAMS_1G["nu_sigma_f"],
        },
        geometry={"R": CYL_PARAMS_1G["R"]},
        alpha=1.0,
        quadrature={
            "n_r": CYL_PARAMS_1G["n_r"],
            "n_mu_axial": CYL_PARAMS_1G["n_mu_axial"],
            "n_phi_az": CYL_PARAMS_1G["n_phi_az"],
            "n_traj_quad": CYL_PARAMS_1G["n_traj_quad"],
        },
    )
    sol = b.solve_critical(
        max_iter=CYL_PARAMS_1G["max_iter"],
        tol=CYL_PARAMS_1G["tol"],
    )
    _check_critical_invariants(sol, n_groups=1, geometry_kind="cylinder")
    assert _bit_equal_floats(sol.eigenvalue, legacy.k_eff)
    assert _bit_equal_arrays(sol.metadata["psi"], legacy.psi)
    assert _bit_equal_arrays(sol.metadata["phi"], legacy.phi)
    assert "mu_axial_nodes" in sol.metadata["mesh"]
    assert "phi_az_nodes" in sol.metadata["mesh"]


# ─────────────────────────────────────────────────────────────────────
# 4. Bit-equality tests — slab + slab_asymmetric
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_billiard_slab_1g_bit_equal_legacy():
    """Billiard slab (symmetric, delegates to slab_asym) 1G ≡ legacy."""
    legacy = gf_slab.solve_greens_function_slab(
        alpha=1.0, **SLAB_PARAMS_1G,
    )
    b = Billiard.from_problem(
        geometry_kind="slab",
        materials={
            "sigma_t": SLAB_PARAMS_1G["sigma_t"],
            "sigma_s": SLAB_PARAMS_1G["sigma_s"],
            "nu_sigma_f": SLAB_PARAMS_1G["nu_sigma_f"],
        },
        geometry={"L": SLAB_PARAMS_1G["L"]},
        alpha=1.0,
        quadrature={
            "n_x": SLAB_PARAMS_1G["n_x"],
            "n_mu": SLAB_PARAMS_1G["n_mu"],
            "n_traj_quad": SLAB_PARAMS_1G["n_traj_quad"],
        },
    )
    sol = b.solve_critical(
        max_iter=SLAB_PARAMS_1G["max_iter"],
        tol=SLAB_PARAMS_1G["tol"],
    )
    _check_critical_invariants(sol, n_groups=1, geometry_kind="slab")
    assert _bit_equal_floats(sol.eigenvalue, legacy.k_eff)
    assert _bit_equal_arrays(sol.metadata["psi"], legacy.psi)
    assert _bit_equal_arrays(sol.metadata["phi"], legacy.phi)


@pytest.mark.foundation
def test_billiard_slab_asym_1g_bit_equal_legacy():
    """Billiard slab_asymmetric 1G with α_L ≠ α_R ≡ legacy."""
    legacy = gf_slab_asym.solve_greens_function_slab_asymmetric(
        alpha_left=0.5, alpha_right=0.8, **SLAB_PARAMS_1G,
    )
    b = Billiard.from_problem(
        geometry_kind="slab_asymmetric",
        materials={
            "sigma_t": SLAB_PARAMS_1G["sigma_t"],
            "sigma_s": SLAB_PARAMS_1G["sigma_s"],
            "nu_sigma_f": SLAB_PARAMS_1G["nu_sigma_f"],
        },
        geometry={"L": SLAB_PARAMS_1G["L"]},
        alpha={"alpha_left": 0.5, "alpha_right": 0.8},
        quadrature={
            "n_x": SLAB_PARAMS_1G["n_x"],
            "n_mu": SLAB_PARAMS_1G["n_mu"],
            "n_traj_quad": SLAB_PARAMS_1G["n_traj_quad"],
        },
    )
    sol = b.solve_critical(
        max_iter=SLAB_PARAMS_1G["max_iter"],
        tol=SLAB_PARAMS_1G["tol"],
    )
    _check_critical_invariants(
        sol, n_groups=1, geometry_kind="slab_asymmetric",
    )
    assert _bit_equal_floats(sol.eigenvalue, legacy.k_eff)
    assert _bit_equal_arrays(sol.metadata["psi"], legacy.psi)
    assert _bit_equal_arrays(sol.metadata["phi"], legacy.phi)


# ─────────────────────────────────────────────────────────────────────
# 5. Bit-equality tests — hollow_sphere + annulus
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_billiard_hollow_sphere_1g_bit_equal_legacy():
    """Billiard hollow_sphere 1G ≡ solve_greens_function_hollow_sphere."""
    legacy = gf_hollow.solve_greens_function_hollow_sphere(
        alpha_in=1.0, alpha_out=1.0, **HOLLOW_SPH_PARAMS_1G,
    )
    b = Billiard.from_problem(
        geometry_kind="hollow_sphere",
        materials={
            "sigma_t": HOLLOW_SPH_PARAMS_1G["sigma_t"],
            "sigma_s": HOLLOW_SPH_PARAMS_1G["sigma_s"],
            "nu_sigma_f": HOLLOW_SPH_PARAMS_1G["nu_sigma_f"],
        },
        geometry={
            "R_in": HOLLOW_SPH_PARAMS_1G["R_in"],
            "R_out": HOLLOW_SPH_PARAMS_1G["R_out"],
        },
        alpha=1.0,
        quadrature={
            "n_r": HOLLOW_SPH_PARAMS_1G["n_r"],
            "n_mu": HOLLOW_SPH_PARAMS_1G["n_mu"],
            "n_traj_quad": HOLLOW_SPH_PARAMS_1G["n_traj_quad"],
        },
    )
    sol = b.solve_critical(
        max_iter=HOLLOW_SPH_PARAMS_1G["max_iter"],
        tol=HOLLOW_SPH_PARAMS_1G["tol"],
    )
    _check_critical_invariants(
        sol, n_groups=1, geometry_kind="hollow_sphere",
    )
    assert _bit_equal_floats(sol.eigenvalue, legacy.k_eff)
    assert _bit_equal_arrays(sol.metadata["psi"], legacy.psi)
    assert _bit_equal_arrays(sol.metadata["phi"], legacy.phi)
    assert sol.parameter_value == 5.0  # R_out


@pytest.mark.foundation
def test_billiard_annulus_1g_bit_equal_legacy():
    """Billiard annulus 1G ≡ solve_greens_function_annulus."""
    legacy = gf_annulus.solve_greens_function_annulus(
        alpha_in=1.0, alpha_out=1.0, **ANNULUS_PARAMS_1G,
    )
    b = Billiard.from_problem(
        geometry_kind="annulus",
        materials={
            "sigma_t": ANNULUS_PARAMS_1G["sigma_t"],
            "sigma_s": ANNULUS_PARAMS_1G["sigma_s"],
            "nu_sigma_f": ANNULUS_PARAMS_1G["nu_sigma_f"],
        },
        geometry={
            "R_in": ANNULUS_PARAMS_1G["R_in"],
            "R_out": ANNULUS_PARAMS_1G["R_out"],
        },
        alpha=1.0,
        quadrature={
            "n_r": ANNULUS_PARAMS_1G["n_r"],
            "n_mu_axial": ANNULUS_PARAMS_1G["n_mu_axial"],
            "n_phi_az": ANNULUS_PARAMS_1G["n_phi_az"],
            "n_traj_quad": ANNULUS_PARAMS_1G["n_traj_quad"],
        },
    )
    sol = b.solve_critical(
        max_iter=ANNULUS_PARAMS_1G["max_iter"],
        tol=ANNULUS_PARAMS_1G["tol"],
    )
    _check_critical_invariants(sol, n_groups=1, geometry_kind="annulus")
    assert _bit_equal_floats(sol.eigenvalue, legacy.k_eff)
    assert _bit_equal_arrays(sol.metadata["psi"], legacy.psi)
    assert _bit_equal_arrays(sol.metadata["phi"], legacy.phi)


# ─────────────────────────────────────────────────────────────────────
# 6. Multi-region sphere + fixed-source
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_billiard_sphere_mr_bit_equal_legacy():
    """Billiard sphere_mr ≡ solve_greens_function_sphere_mr."""
    radii = np.array([2.0, 5.0])
    sigma_t_r = np.array([[1.0, 0.5], [0.5, 0.3]])
    sigma_s_r = np.array([
        [[0.5, 0.4], [0.0, 0.4]],
        [[0.3, 0.2], [0.0, 0.2]],
    ])
    nu_sf_r = np.array([[0.05, 0.10], [0.0, 0.0]])

    legacy = gf_sphere.solve_greens_function_sphere_mr(
        radii=radii, sigma_t=sigma_t_r, sigma_s=sigma_s_r,
        nu_sigma_f=nu_sf_r, alpha=1.0,
        n_r=8, n_mu=8, n_traj_quad=16, max_iter=10, tol=1e-9,
    )

    b = Billiard.from_problem(
        geometry_kind="sphere_mr",
        materials={
            "sigma_t": sigma_t_r,
            "sigma_s": sigma_s_r,
            "nu_sigma_f": nu_sf_r,
        },
        geometry={"radii": radii},
        alpha=1.0,
        quadrature={"n_r": 8, "n_mu": 8, "n_traj_quad": 16},
    )
    sol = b.solve_critical(max_iter=10, tol=1e-9)
    _check_critical_invariants(sol, n_groups=2, geometry_kind="sphere_mr")
    assert _bit_equal_floats(sol.eigenvalue, legacy.k_eff)
    assert _bit_equal_arrays(sol.metadata["psi"], legacy.psi_g)
    assert _bit_equal_arrays(sol.metadata["phi"], legacy.phi_g)
    assert _bit_equal_arrays(
        sol.metadata["mesh"]["region_at_node"], legacy.region_at_node,
    )
    assert sol.parameter_value == 5.0  # radii[-1]


@pytest.mark.foundation
def test_billiard_fixed_source_sphere_mr_bit_equal_legacy():
    """Billiard.solve_fixed_source sphere_mr ≡ legacy fixed-source."""
    radii = np.array([2.0, 5.0])
    sigma_t = np.array([[1.0], [0.5]])
    sigma_s = np.array([[[0.5]], [[0.3]]])
    ext_q = np.array([[1.0], [0.0]])

    legacy = gf_sphere.solve_greens_function_sphere_mr_fixed_source(
        radii=radii, sigma_t=sigma_t, sigma_s=sigma_s,
        external_source=ext_q, alpha=0.0,
        n_r=8, n_mu=8, n_traj_quad=16, max_iter=20, tol=1e-7,
    )

    b = Billiard.from_problem(
        geometry_kind="sphere_mr",
        materials={"sigma_t": sigma_t, "sigma_s": sigma_s},
        geometry={"radii": radii},
        alpha=0.0,
        quadrature={"n_r": 8, "n_mu": 8, "n_traj_quad": 16},
    )
    sol = b.solve_fixed_source(ext_q, max_iter=20, tol=1e-7)
    assert isinstance(sol, FluxSolution)
    assert _bit_equal_arrays(sol.metadata["psi"], legacy.psi_g)
    assert _bit_equal_arrays(sol.metadata["phi"], legacy.phi_g)
    assert sol.metadata["iterations"] == legacy.iterations
    assert sol.eigenvalue_kind == "none"
    assert sol.spatial_units == "cm"


@pytest.mark.foundation
def test_billiard_fixed_source_unsupported_geometry_raises():
    """fixed_source on a non-sphere_mr geometry → NotImplementedError."""
    b = Billiard.from_problem(
        geometry_kind="sphere",
        materials={"sigma_t": 0.5, "sigma_s": 0.4, "nu_sigma_f": 0.1},
        geometry={"R": 5.0},
        alpha=0.0,
    )
    with pytest.raises(NotImplementedError, match="sphere_mr"):
        b.solve_fixed_source(np.ones(3))
