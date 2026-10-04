"""Foundation tests for the :class:`Billiard` class (Phase D).

These tests pin the Billiard facade's bit-equality with the underlying
``solve_greens_function_*`` entry points across the geometry families
that direct :class:`Billiard` construction supports today (sphere /
cylinder / slab / slab_asymmetric, in 1G and MG variants). They do
NOT re-test the underlying solvers' correctness — that's the existing
suite's job. They DO assert that wrapping a solve in
:meth:`Billiard.solve_critical` returns a shared
:class:`~orpheus.derivations.common.solution_types.CriticalSolution`
whose contents are bit-for-bit identical to the original return.

Phase D consumes :class:`~orpheus.geometry.structured_geometry.StructuredGeometry`
directly; its albedos are the geometry's declared laws, and a slab whose
two faces declare different albedos selects the asymmetric-slab branch.

Tests are tagged ``foundation`` because they verify a software
contract (the facade's bit-equal preservation) rather than an L0/L1
mathematical claim about a solver.
"""
from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.derivations.common.solution_types import CriticalSolution
from orpheus.derivations.continuous.trajectory_resolvent import Billiard
from orpheus.derivations.continuous.trajectory_resolvent import (
    greens_function as gf_sphere,
    greens_function_cylinder as gf_cyl,
    greens_function_slab as gf_slab,
    greens_function_slab_asymmetric as gf_slab_asym,
)
from orpheus.geometry import BC, CoordSystem, StructuredGeometry


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
# Inline Mixture / StructuredGeometry helpers
# ─────────────────────────────────────────────────────────────────────


def _mixture_from_xs(
    sigma_t: float | np.ndarray,
    sigma_s: float | np.ndarray,
    nu_sigma_f: float | np.ndarray,
    chi: np.ndarray | None = None,
) -> Mixture:
    """Build a minimal :class:`Mixture` from raw XS scalars / arrays."""
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
    # Synthetic XS for trajectory-resolvent billiard test (Phase E):
    # no physical energy grid.
    #
    # This is a PRODUCING (multiplying) medium: ``SigP = νΣ_f > 0`` drives a
    # fission source through ``chi`` (the billiard solver reads ``SigP`` and
    # ``chi``, never ``SigF``). The χ guard keys on PRODUCTION (``SigP > 0``),
    # so ``is_producing`` is True and its χ — the default simplex
    # ``[1, 0, ...]`` — is correctly required to be a probability simplex.
    # ``SigF`` is the fission cross-section, a distinct quantity the billiard
    # path never reads; it stays zero (no SigF stand-in needed).
    return Mixture(
        SigC=np.zeros(ng),
        SigL=np.zeros(ng),
        SigF=np.zeros(ng),
        SigP=nu_sf_arr.copy(),
        SigT=sig_t_arr.copy(),
        SigS=[csr_matrix(sig_s_arr.copy())],
        Sig2=[csr_matrix((ng, ng))],
        chi=chi_arr,
    )


def _sphere_geom(R_cm: float) -> StructuredGeometry:
    """Closed homogeneous sphere geometry at radius :math:`R_{\\rm cm}`."""
    return StructuredGeometry(
        coord=CoordSystem.SPHERICAL,
        breakpoints=(0.0, float(R_cm)),
        mat_ids=(0,),
        boundaries=(BC.reflective,),
    )


def _cylinder_geom(R_cm: float) -> StructuredGeometry:
    """Closed homogeneous cylinder geometry at radius :math:`R_{\\rm cm}`."""
    return StructuredGeometry(
        coord=CoordSystem.CYLINDRICAL,
        breakpoints=(0.0, float(R_cm)),
        mat_ids=(0,),
        boundaries=(BC.reflective,),
    )


def _slab_geom(L_cm: float) -> StructuredGeometry:
    """Symmetric slab geometry at FULL thickness :math:`L_{\\rm cm}`.

    StructuredGeometry convention: ``domain_extent_cm`` is the full
    slab width, which matches the ``L`` argument of
    :func:`solve_greens_function_slab`.
    """
    return StructuredGeometry(
        coord=CoordSystem.CARTESIAN,
        breakpoints=(0.0, float(L_cm)),
        mat_ids=(0,),
        boundaries=(BC.vacuum, BC.vacuum),
    )


def _law(albedo: float) -> BC:
    """The law declaring a specular albedo: vacuum, mirror, or partial."""
    if albedo == 0.0:
        return BC.vacuum
    if albedo == 1.0:
        return BC.reflective
    return BC("partial", {"albedo": albedo})


def _laws(geometry: StructuredGeometry, alpha) -> StructuredGeometry:
    """``geometry`` with the laws that declare ``alpha``: a scalar on every
    boundary point, or ``{"alpha_left", "alpha_right"}`` on a slab's faces.
    Billiard reads its albedos from these laws (its ``alpha`` parameter
    retired in P1 step 2b)."""
    if isinstance(alpha, dict):
        laws = (_law(alpha["alpha_left"]), _law(alpha["alpha_right"]))
    else:
        laws = tuple(_law(alpha) for _ in geometry.boundaries)
    return replace(geometry, boundaries=laws)


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
def test_billiard_sphere_one_endpoint_rank_1():
    """One-endpoint orbit-space class → closure_rank == 1."""
    b = Billiard(
        materials={0: _mixture_from_xs(0.5, 0.4, 0.1)},
        geometry=_laws(_sphere_geom(5.0), 1.0),
    )
    assert b.closure_rank == 1
    assert b.geometry_kind == "sphere"
    assert b.alpha_payload == {"alpha": 1.0}


@pytest.mark.foundation
def test_billiard_slab_asymmetric_two_endpoint_rank_2():
    """Two-endpoint orbit-space class → closure_rank == 2.

    The ``slab_asymmetric`` family is selected automatically when
    ``alpha`` is a dict carrying ``alpha_left`` / ``alpha_right``
    keys on a SLB :class:`StructuredGeometry`.
    """
    b = Billiard(
        materials={0: _mixture_from_xs(0.5, 0.4, 0.1)},
        geometry=_laws(_slab_geom(5.0), {"alpha_left": 0.7, "alpha_right": 0.9}),
    )
    assert b.closure_rank == 2
    assert b.geometry_kind == "slab_asymmetric"
    assert b.alpha_payload == {"alpha_left": 0.7, "alpha_right": 0.9}


@pytest.mark.foundation
def test_billiard_scalar_alpha_on_slab_stays_symmetric():
    """A scalar alpha on a slab :class:`StructuredGeometry` stays in the
    symmetric ``slab`` family — never silently promoted to
    ``slab_asymmetric``.
    """
    b = Billiard(
        materials={0: _mixture_from_xs(0.5, 0.4, 0.1)},
        geometry=_laws(_slab_geom(5.0), 0.5),
    )
    assert b.geometry_kind == "slab"
    assert b.closure_rank == 1
    assert b.alpha_payload == {"alpha": 0.5}


@pytest.mark.foundation
def test_billiard_reads_its_albedo_from_the_declared_law():
    """RE-POSED from ``test_billiard_with_alpha_returns_modified_copy``: the
    ``alpha`` parameter and ``with_alpha`` retired (P1 step 2b), so a
    different albedo is a different declared law on the geometry, and it
    moves only the albedo payload."""
    materials = {0: _mixture_from_xs(0.5, 0.4, 0.1)}
    b = Billiard(materials=materials, geometry=_laws(_sphere_geom(5.0), 1.0))
    b2 = Billiard(materials=materials, geometry=_laws(_sphere_geom(5.0), 0.5))
    assert b.alpha_payload == {"alpha": 1.0}
    assert b2.alpha_payload == {"alpha": 0.5}
    assert b.geometry_kind == b2.geometry_kind
    assert b.xs_payload == b2.xs_payload
    assert b.geometry_payload == b2.geometry_payload


@pytest.mark.foundation
def test_billiard_dispatches_every_coordinate_system():
    """A homogeneous body in every ``CoordSystem`` member routes to an arm:
    the dispatch covers the closed enum (P1 step 2's re-pose of the old
    tag-refusal gate, re-posed again when the kind table became one match
    over the body shape in step 2b)."""
    materials = {0: _mixture_from_xs(0.5, 0.4, 0.1)}
    geometries = {
        CoordSystem.CARTESIAN: _slab_geom(5.0),
        CoordSystem.SPHERICAL: _sphere_geom(5.0),
        CoordSystem.CYLINDRICAL: _cylinder_geom(5.0),
    }
    assert set(geometries) == set(CoordSystem)
    kinds = {
        coord: Billiard(materials=materials, geometry=g).geometry_kind
        for coord, g in geometries.items()
    }
    assert kinds == {
        CoordSystem.CARTESIAN: "slab",
        CoordSystem.SPHERICAL: "sphere",
        CoordSystem.CYLINDRICAL: "cylinder",
    }


# ─────────────────────────────────────────────────────────────────────
# 2. Bit-equality tests — sphere
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_billiard_sphere_1g_bit_equal_legacy():
    """Billiard.solve_critical sphere 1G ≡ solve_greens_function_sphere."""
    legacy = gf_sphere.solve_greens_function_sphere(
        alpha=1.0, **SPHERE_PARAMS_1G,
    )
    b = Billiard(
        materials={0: _mixture_from_xs(
            SPHERE_PARAMS_1G["sigma_t"],
            SPHERE_PARAMS_1G["sigma_s"],
            SPHERE_PARAMS_1G["nu_sigma_f"],
        )},
        geometry=_laws(_sphere_geom(SPHERE_PARAMS_1G["R"]), 1.0),
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
    b = Billiard(
        materials={0: _mixture_from_xs(
            SPHERE_PARAMS_MG["sigma_t"],
            SPHERE_PARAMS_MG["sigma_s"],
            SPHERE_PARAMS_MG["nu_sigma_f"],
        )},
        geometry=_laws(_sphere_geom(SPHERE_PARAMS_MG["R"]), 1.0),
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
    b = Billiard(
        materials={0: _mixture_from_xs(
            CYL_PARAMS_1G["sigma_t"],
            CYL_PARAMS_1G["sigma_s"],
            CYL_PARAMS_1G["nu_sigma_f"],
        )},
        geometry=_laws(_cylinder_geom(CYL_PARAMS_1G["R"]), 1.0),
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
    """Billiard slab (symmetric) 1G ≡ legacy."""
    legacy = gf_slab.solve_greens_function_slab(
        alpha=1.0, **SLAB_PARAMS_1G,
    )
    b = Billiard(
        materials={0: _mixture_from_xs(
            SLAB_PARAMS_1G["sigma_t"],
            SLAB_PARAMS_1G["sigma_s"],
            SLAB_PARAMS_1G["nu_sigma_f"],
        )},
        geometry=_laws(_slab_geom(SLAB_PARAMS_1G["L"]), 1.0),
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
    b = Billiard(
        materials={0: _mixture_from_xs(
            SLAB_PARAMS_1G["sigma_t"],
            SLAB_PARAMS_1G["sigma_s"],
            SLAB_PARAMS_1G["nu_sigma_f"],
        )},
        geometry=_laws(_slab_geom(SLAB_PARAMS_1G["L"]), {"alpha_left": 0.5, "alpha_right": 0.8}),
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
    """fixed_source on a non-sphere_mr geometry → NotImplementedError."""
    b = Billiard(
        materials={0: _mixture_from_xs(0.5, 0.4, 0.1)},
        geometry=_laws(_sphere_geom(5.0), 0.0),
    )
    with pytest.raises(NotImplementedError, match="sphere_mr"):
        b.solve_fixed_source(np.ones(3))


# ─────────────────────────────────────────────────────────────────────
# #405 P2 step 7b.2.3 — the billiard refuses a body material it would read only in part (gates R7b2.3.<n>)
# ─────────────────────────────────────────────────────────────────────


def _partial_read_mixtures():
    """Isotropic fuel A, the library's moderator B (a P1 moment, mean cosine 0.6), and fuel A with an (n,2n) matrix."""
    from orpheus.derivations.common.xs_library import get_mixture, get_xs, make_mixture

    xs = get_xs("A", "2g")
    common = dict(sig_t=xs["sig_t"], sig_c=xs["sig_c"], sig_f=xs["sig_f"], nu=xs["nu"], chi=xs["chi"], sig_s=xs["sig_s"])
    isotropic = make_mixture(**common)
    n2n = make_mixture(**common, sig_2=np.array([[0.0, 0.01], [0.0, 0.0]]))
    return isotropic, get_mixture("B", "2g"), n2n


def _layered_sphere(outer=BC.reflective) -> StructuredGeometry:
    return StructuredGeometry(coord=CoordSystem.SPHERICAL, breakpoints=(0.0, 0.5, 2.0), mat_ids=(0, 1), boundaries=(outer,))


@pytest.mark.foundation
def test_r7b2_3_1_a_body_material_read_in_part_is_refused() -> None:
    """Every trajectory-resolvent solver reads ``SigS[0]`` alone and no ``Sig2``: a body material with a P1 moment
    or an (n,2n) reaction would be solved as another problem, silently. Refused at the root, keyed; the
    isotropic control is admitted (the activation)."""
    isotropic, anisotropic, n2n = _partial_read_mixtures()
    assert len(anisotropic.SigS) == 2 and anisotropic.SigS[1].count_nonzero() > 0  # the activation
    assert n2n.Sig2[0].count_nonzero() > 0
    Billiard(geometry=_layered_sphere(), materials={0: isotropic, 1: isotropic})
    with pytest.raises(ValueError, match="material 1 scatters anisotropically"):
        Billiard(geometry=_layered_sphere(), materials={0: isotropic, 1: anisotropic})
    with pytest.raises(ValueError, match=r"material 1 carries an \(n,2n\) reaction"):
        Billiard(geometry=_layered_sphere(), materials={0: isotropic, 1: n2n})


@pytest.mark.foundation
def test_r7b2_3_2_a_material_the_body_does_not_hold_is_never_read() -> None:
    """A spectator material (id 7, not in the geometry) carrying a P1 moment and an (n,2n) matrix is admitted: the
    refusal reads only the body's materials (the leak principle)."""
    isotropic, anisotropic, n2n = _partial_read_mixtures()
    billiard = Billiard(geometry=_layered_sphere(), materials={0: isotropic, 1: isotropic, 7: anisotropic, 8: n2n})
    assert billiard.geometry_kind == "sphere_mr"


@pytest.mark.foundation
def test_r7b2_3_3_the_law_door_keeps_its_message() -> None:
    """The refusal runs after the route: a white outer law on an anisotropic body is refused with the LAW's message,
    so the doors' messages do not depend on which defect the input also has."""
    _, anisotropic, _ = _partial_read_mixtures()
    with pytest.raises(NotImplementedError, match="another angular shape"):
        Billiard(geometry=_layered_sphere(BC("white")), materials={0: anisotropic, 1: anisotropic})
