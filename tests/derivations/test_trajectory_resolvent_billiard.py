"""Foundation tests for the :class:`Billiard` class.

These tests pin the Billiard facade's bit-equality with the underlying
``solve_greens_function_*`` entry points across the geometry families
that :meth:`Billiard.from_problem` can construct today (sphere /
cylinder / slab / slab_asymmetric, in 1G and MG variants). They do
NOT re-test the underlying solvers' correctness — that's the existing
suite's job. They DO assert that wrapping a solve in
:meth:`Billiard.solve_critical` returns a shared
:class:`~orpheus.derivations.common.solution_types.CriticalSolution`
whose contents are bit-for-bit identical to the original return.

The dispatcher additionally supports ``hollow_sphere``, ``annulus``,
``sphere_mr`` for callers that construct ``Billiard`` instances
directly, but :meth:`from_problem` cannot construct those today (the
:class:`~orpheus.derivations.common.geometry_spec.GeometrySpec` schema
does not yet carry inner radii or per-region zone descriptions).
The bit-equality coverage for those geometries returns once Step 3
of the input-cleanup track lands the multi-region GeometrySpec
extension; for now those tests are removed (they cannot be expressed
through the production-protocol factory).

Tests are tagged ``foundation`` because they verify a software
contract (the facade's bit-equal preservation) rather than an L0/L1
mathematical claim about a solver.
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.derivations.common.geometry_spec import GeometrySpec
from orpheus.derivations.common.solution_types import CriticalSolution
from orpheus.derivations.continuous.trajectory_resolvent import Billiard
from orpheus.derivations.continuous.trajectory_resolvent import (
    greens_function as gf_sphere,
    greens_function_cylinder as gf_cyl,
    greens_function_slab as gf_slab,
    greens_function_slab_asymmetric as gf_slab_asym,
)
from orpheus.geometry.mesh import BC


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


# ─────────────────────────────────────────────────────────────────────
# Inline Mixture / GeometrySpec helpers
# ─────────────────────────────────────────────────────────────────────


def _mixture_from_xs(
    sigma_t: float | np.ndarray,
    sigma_s: float | np.ndarray,
    nu_sigma_f: float | np.ndarray,
    chi: np.ndarray | None = None,
) -> Mixture:
    """Build a minimal :class:`Mixture` from raw XS scalars / arrays.

    Used by these tests to feed the production-protocol factory
    surface from the same raw cross sections the legacy
    ``solve_greens_function_*`` entry points consume directly. The
    factory then re-extracts those XS through
    :func:`_mixture_to_solver_xs_payload` — bit-equality of the
    resulting solver call is the invariant under test.
    """
    sig_t_arr = np.atleast_1d(np.asarray(sigma_t, dtype=float))
    sig_s_arr = np.atleast_2d(np.asarray(sigma_s, dtype=float))
    nu_sf_arr = np.atleast_1d(np.asarray(nu_sigma_f, dtype=float))
    if sig_s_arr.shape != (sig_t_arr.size, sig_t_arr.size):
        sig_s_arr = sig_s_arr.reshape(sig_t_arr.size, sig_t_arr.size)
    if chi is None:
        chi_arr = np.zeros(sig_t_arr.size, dtype=float)
        chi_arr[0] = 1.0
    else:
        chi_arr = np.atleast_1d(np.asarray(chi, dtype=float))
    ng = sig_t_arr.size
    return Mixture(
        SigC=np.zeros(ng),
        SigL=np.zeros(ng),
        SigF=np.zeros(ng),
        SigP=nu_sf_arr.copy(),
        SigT=sig_t_arr.copy(),
        SigS=[csr_matrix(sig_s_arr.copy())],
        Sig2=csr_matrix((ng, ng)),
        chi=chi_arr,
        eg=np.logspace(7, -3, ng + 1),
    )


def _sphere_spec(R_cm: float, n_groups: int = 1) -> GeometrySpec:
    """Closed homogeneous sphere spec at radius :math:`R_{\\rm cm}`."""
    return GeometrySpec(
        geometry="sphere",
        critical_dimension_mfp=None,
        critical_dimension_cm=float(R_cm),
        n_groups=n_groups,
        bc_left=BC.reflective,
        bc_right=BC.reflective,
    )


def _cylinder_spec(R_cm: float, n_groups: int = 1) -> GeometrySpec:
    """Closed homogeneous cylinder spec at radius :math:`R_{\\rm cm}`."""
    return GeometrySpec(
        geometry="cylinder",
        critical_dimension_mfp=None,
        critical_dimension_cm=float(R_cm),
        n_groups=n_groups,
        bc_left=BC.reflective,
        bc_right=BC.reflective,
    )


def _slab_spec(L_cm: float, n_groups: int = 1) -> GeometrySpec:
    """Symmetric slab spec at full thickness :math:`L_{\\rm cm}`.

    GeometrySpec convention: ``critical_dimension_cm`` is the slab
    half-thickness; :attr:`GeometrySpec.domain_extent_cm` then
    returns ``2 * critical_dimension_cm``, which matches the ``L``
    argument of :func:`solve_greens_function_slab`.
    """
    return GeometrySpec(
        geometry="slab",
        critical_dimension_mfp=None,
        critical_dimension_cm=0.5 * float(L_cm),
        n_groups=n_groups,
        bc_left=BC.vacuum,
        bc_right=BC.vacuum,
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
        materials={0: _mixture_from_xs(0.5, 0.4, 0.1)},
        geometry_spec=_sphere_spec(5.0),
        alpha=1.0,
    )
    assert b.closure_rank == 1
    assert b.geometry_kind == "sphere"
    assert b.alpha_payload == {"alpha": 1.0}


@pytest.mark.foundation
def test_billiard_from_problem_slab_asymmetric_two_endpoint_rank_2():
    """Two-endpoint orbit-space class → closure_rank == 2.

    The ``slab_asymmetric`` family is selected automatically when the
    factory's ``alpha`` argument is a dict carrying
    ``alpha_left`` / ``alpha_right`` keys on a slab GeometrySpec.
    """
    b = Billiard.from_problem(
        materials={0: _mixture_from_xs(0.5, 0.4, 0.1)},
        geometry_spec=_slab_spec(5.0),
        alpha={"alpha_left": 0.7, "alpha_right": 0.9},
    )
    assert b.closure_rank == 2
    assert b.geometry_kind == "slab_asymmetric"
    assert b.alpha_payload == {"alpha_left": 0.7, "alpha_right": 0.9}


@pytest.mark.foundation
def test_billiard_from_problem_scalar_alpha_on_slab_stays_symmetric():
    """A scalar alpha on a slab GeometrySpec stays in the symmetric
    ``slab`` family — never silently promoted to ``slab_asymmetric``.

    The asymmetric family is selected ONLY by the alpha-dict shape
    (presence of ``alpha_left`` / ``alpha_right`` keys). A scalar
    alpha lands in the one-endpoint payload ``{"alpha": float}``
    with ``closure_rank=1``, regardless of geometry. Pin this so a
    future refactor cannot accidentally re-introduce the dead
    "scalar → broadcast to slab_asymmetric" branch.
    """
    b = Billiard.from_problem(
        materials={0: _mixture_from_xs(0.5, 0.4, 0.1)},
        geometry_spec=_slab_spec(5.0),
        alpha=0.5,
    )
    assert b.geometry_kind == "slab"
    assert b.closure_rank == 1
    assert b.alpha_payload == {"alpha": 0.5}


@pytest.mark.foundation
def test_billiard_with_alpha_returns_modified_copy():
    """with_alpha returns a copy with new alpha; original unchanged."""
    b = Billiard.from_problem(
        materials={0: _mixture_from_xs(0.5, 0.4, 0.1)},
        geometry_spec=_sphere_spec(5.0),
        alpha=1.0,
    )
    b2 = b.with_alpha(0.5)
    assert b.alpha_payload == {"alpha": 1.0}
    assert b2.alpha_payload == {"alpha": 0.5}
    assert b.geometry_kind == b2.geometry_kind
    # The materials/geometry_spec pair is preserved; only alpha differs.
    assert b.materials is b2.materials or b.materials == b2.materials
    assert b.geometry_spec == b2.geometry_spec
    # The synthesised solver-facing payload is bit-equal too.
    assert b.xs_payload == b2.xs_payload
    assert b.geometry_payload == b2.geometry_payload


@pytest.mark.foundation
def test_billiard_unsupported_geometry_spec_raises():
    """Infinite-medium / ISLC GeometrySpecs are rejected.

    The factory supports ``slab`` / ``sphere`` / ``cylinder``;
    constructing on ``infinite`` (no boundary, no chord algebra) or
    ``ISLC`` (still on the roadmap) raises ``ValueError`` so the user
    learns immediately rather than crashing inside the dispatcher.
    """
    inf_spec = GeometrySpec(
        geometry="infinite",
        critical_dimension_mfp=None,
        critical_dimension_cm=None,
        n_groups=1,
    )
    with pytest.raises(
        ValueError, match="cannot construct on geometry_spec"
    ):
        Billiard.from_problem(
            materials={0: _mixture_from_xs(0.5, 0.4, 0.1)},
            geometry_spec=inf_spec,
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
        materials={0: _mixture_from_xs(
            SPHERE_PARAMS_1G["sigma_t"],
            SPHERE_PARAMS_1G["sigma_s"],
            SPHERE_PARAMS_1G["nu_sigma_f"],
        )},
        geometry_spec=_sphere_spec(SPHERE_PARAMS_1G["R"]),
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
        materials={0: _mixture_from_xs(
            SPHERE_PARAMS_MG["sigma_t"],
            SPHERE_PARAMS_MG["sigma_s"],
            SPHERE_PARAMS_MG["nu_sigma_f"],
        )},
        geometry_spec=_sphere_spec(SPHERE_PARAMS_MG["R"], n_groups=2),
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
        materials={0: _mixture_from_xs(
            CYL_PARAMS_1G["sigma_t"],
            CYL_PARAMS_1G["sigma_s"],
            CYL_PARAMS_1G["nu_sigma_f"],
        )},
        geometry_spec=_cylinder_spec(CYL_PARAMS_1G["R"]),
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
        materials={0: _mixture_from_xs(
            SLAB_PARAMS_1G["sigma_t"],
            SLAB_PARAMS_1G["sigma_s"],
            SLAB_PARAMS_1G["nu_sigma_f"],
        )},
        geometry_spec=_slab_spec(SLAB_PARAMS_1G["L"]),
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
        materials={0: _mixture_from_xs(
            SLAB_PARAMS_1G["sigma_t"],
            SLAB_PARAMS_1G["sigma_s"],
            SLAB_PARAMS_1G["nu_sigma_f"],
        )},
        geometry_spec=_slab_spec(SLAB_PARAMS_1G["L"]),
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
# 5. Fixed-source — unsupported-geometry rejection
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_billiard_fixed_source_unsupported_geometry_raises():
    """fixed_source on a non-sphere_mr geometry → NotImplementedError.

    The dispatcher currently routes fixed-source solves only for the
    multi-region sphere family (Garcia 2021 benchmarks). Constructing
    a homogeneous-sphere :class:`Billiard` and calling
    :meth:`solve_fixed_source` on it must raise ``NotImplementedError``
    with a message that names ``sphere_mr`` so the user learns the
    supported shape.
    """
    b = Billiard.from_problem(
        materials={0: _mixture_from_xs(0.5, 0.4, 0.1)},
        geometry_spec=_sphere_spec(5.0),
        alpha=0.0,
    )
    with pytest.raises(NotImplementedError, match="sphere_mr"):
        b.solve_fixed_source(np.ones(3))


# ─────────────────────────────────────────────────────────────────────
# Note on removed coverage
# ─────────────────────────────────────────────────────────────────────
#
# The hollow_sphere / annulus / sphere_mr bit-equality tests have been
# removed because :meth:`Billiard.from_problem` cannot construct on
# those geometries — :class:`GeometrySpec` does not yet carry the
# inner-radius / per-region-zone fields these families require, and
# the factory rejects unsupported specs. Those tests will return when
# the multi-region :class:`GeometrySpec` extension lands (Step 3 of
# the input-cleanup track). The dispatcher continues to support them
# for callers that build :class:`Billiard` instances directly via
# the dataclass constructor.
