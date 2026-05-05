r"""Compatibility smoke tests for the new sood_registry.

These tests are the **bridge gate** between the method-agnostic
:mod:`orpheus.derivations.continuous.sood_registry` package and:

* Semi-analytical reference solvers (F_N, transfer-matrix k_inf) —
  consume Mixture → numpy arrays via
  :func:`extractors.mixture_to_fn_arrays`.
* Production discrete solvers (CP, SN, MOC) — consume
  ``materials: dict[int, Mixture]`` + ``mesh: Mesh1D`` directly.

The gates are intentionally minimal — they verify that the case
objects are *consumable* by both flavors of solver, NOT that the
answers are correct (the existing ``test_fn_la13511_*`` files cover
correctness; production-solver cross-checks at full Sood precision
land in Phase B).

Test classification:

* :func:`test_*_consumable_by_*` — ``foundation``-tagged structural
  invariants of the registry / extractor interface.
* :func:`test_PUa_1_0_IN_kinf_consumable_by_kinf_homogeneous` —
  ``foundation`` smoke that the existing transfer-matrix k_inf
  reference solver can consume the new registry directly and
  reproduces the published value.
* :func:`test_*_smoke_via_solve_cp` /
  :func:`test_*_smoke_via_solve_sn` — ``l1`` smoke that production
  solvers accept the registry's mesh + materials. These are slow
  (full power iteration) and run with ``@pytest.mark.slow`` so the
  default ``foundation`` collection skips them.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.derivations.common.eigenvalue import kinf_homogeneous
from orpheus.derivations.continuous.sood_registry import (
    LA13511_CASES,
    PU_2_0_IN,
    PUA_1_0_IN,
    UA_1_0_SL_STUB,
    UA_1_0_SP_STUB,
    La13511Case,
    build_materials,
    build_mesh,
)
from orpheus.geometry.coord import CoordSystem
from orpheus.geometry.mesh import BC, Mesh1D


def _xs_from_mixture(mix: Mixture):
    """Inline replacement for the retired ``mixture_to_fn_arrays``.

    Returns ``(sigma_t, sigma_s_p0, nu_sigma_f, chi)`` from the
    production-protocol Mixture surface. Phase D retired the shared
    extractor; tests that historically used it now read the same
    fields directly off the Mixture.
    """
    sigma_t = np.asarray(mix.SigT, dtype=float)
    sigma_s = mix.SigS[0].toarray().astype(float)
    nu_sigma_f = np.asarray(mix.SigP, dtype=float)
    chi = np.asarray(mix.chi, dtype=float)
    return sigma_t, sigma_s, nu_sigma_f, chi


# ═══════════════════════════════════════════════════════════════════
# Foundation: registry shape invariants
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.parametrize("case_id", list(LA13511_CASES.keys()))
def test_case_materials_is_dict_of_mixture(case_id: str) -> None:
    """Every registered case has materials: dict[int, Mixture]."""
    case = LA13511_CASES[case_id]
    assert isinstance(case.materials, dict), (
        f"{case_id}: materials must be a dict, got {type(case.materials)}"
    )
    assert len(case.materials) >= 1, f"{case_id}: at least one material required"
    for mat_id, mixture in case.materials.items():
        assert isinstance(mat_id, int), (
            f"{case_id}: material key must be int, got {type(mat_id)}"
        )
        assert isinstance(mixture, Mixture), (
            f"{case_id}: materials[{mat_id}] must be Mixture, got {type(mixture)}"
        )


@pytest.mark.foundation
@pytest.mark.parametrize("case_id", list(LA13511_CASES.keys()))
def test_case_geometry_kind_well_formed(case_id: str) -> None:
    """Every registered case carries a valid geometry_kind tag."""
    case = LA13511_CASES[case_id]
    assert case.geometry_kind in {"slab", "sphere", "cylinder", "infinite"}


@pytest.mark.foundation
@pytest.mark.parametrize(
    ("case_id", "expected_geom"),
    [
        ("PUa-1-0-IN", "infinite"),
        ("PU-2-0-IN", "infinite"),
        ("Ua-1-0-SL", "slab"),
        ("Ua-1-0-CY", "cylinder"),
        ("Ua-1-0-SP", "sphere"),
    ],
)
def test_case_geometry_matches_expected(case_id: str, expected_geom: str) -> None:
    """Geometry kind matches the case-id suffix."""
    assert LA13511_CASES[case_id].geometry_kind == expected_geom


# ═══════════════════════════════════════════════════════════════════
# Foundation: extractor round-trip — published values preserved exactly
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_PUa_1_0_IN_extractor_matches_published_xs() -> None:
    r"""Extractor returns Sood's published XS for PUa-1-0-IN exactly.

    Sood Table 2: :math:`\Sigma_t = 0.32640`, :math:`\Sigma_s = 0.225216`,
    :math:`\nu \Sigma_f = 3.24 \cdot 0.0816 = 0.264384`, :math:`\chi = 1.0`.
    """
    sigma_t, sigma_s, nu_sigma_f, chi = _xs_from_mixture(
        PUA_1_0_IN.materials[0]
    )
    np.testing.assert_allclose(sigma_t, [0.32640], atol=0)
    np.testing.assert_allclose(sigma_s, [[0.225216]], atol=0)
    np.testing.assert_allclose(nu_sigma_f, [3.24 * 0.0816], atol=0)
    np.testing.assert_allclose(chi, [1.0], atol=0)


@pytest.mark.foundation
def test_PU_2_0_IN_extractor_matches_published_xs() -> None:
    r"""Extractor returns Sood's published 2G XS for PU-2-0-IN exactly.

    Verifies the ORPHEUS-order conversion (Sood g=2 fast → ORPHEUS g=0
    fast) is preserved end-to-end through the Mixture round-trip.
    """
    sigma_t, sigma_s, nu_sigma_f, chi = _xs_from_mixture(
        PU_2_0_IN.materials[0]
    )
    # ORPHEUS convention: g=0 fast, g=1 slow.
    np.testing.assert_allclose(sigma_t, [0.2208, 0.3360], atol=0)
    np.testing.assert_allclose(
        sigma_s,
        [[0.0792, 0.0432], [0.0, 0.23616]],
        atol=0,
    )
    np.testing.assert_allclose(
        nu_sigma_f, [3.10 * 0.0936, 2.93 * 0.08544], atol=0
    )
    np.testing.assert_allclose(chi, [0.575, 0.425], atol=0)


# ═══════════════════════════════════════════════════════════════════
# Foundation: mesh template build conventions
# ═══════════════════════════════════════════════════════════════════


def _critical_dimension_cm(case: La13511Case) -> float:
    """Derive the critical dimension in cm from the case.

    ``cm = mfp / Σ_t``, with mfp from :attr:`La13511Truth.critical_dimension_mfp`
    and Σ_t from the primary mixture's group-0 SigT.
    """
    return float(case.truth.critical_dimension_mfp) / float(case.materials[0].SigT[0])


@pytest.mark.foundation
def test_slab_mesh_builds_full_symmetric_domain() -> None:
    r"""Slab mesh covers the FULL symmetric slab :math:`[0, 2a]`.

    F_N convention: ``truth.critical_dimension_mfp`` corresponds to the
    half-thickness :math:`a`; in cm this is
    :math:`a_{\rm cm} = a_{\rm mfp} / \Sigma_t`. The constructed mesh
    covers ``[0, 2a_{cm}]`` with vacuum BCs at both ends — matching
    the F_N vs Variant α slab xverif test convention.
    """
    case = UA_1_0_SL_STUB
    mesh = build_mesh(case, n_cells=16)
    assert mesh.coord == CoordSystem.CARTESIAN
    expected_total = 2.0 * _critical_dimension_cm(case)
    assert mesh.total_width == pytest.approx(expected_total, abs=1e-12)
    assert mesh.bc_left == BC.vacuum
    assert mesh.bc_right == BC.vacuum
    assert mesh.N == 16


@pytest.mark.foundation
def test_sphere_mesh_builds_radial_domain_with_reflective_centre() -> None:
    """Sphere mesh covers ``[0, R]``; centreline reflective is implicit.

    :meth:`Mesh1D.from_geometry` for SPH sets ``bc_left=None`` (the
    centreline at the coordinate origin is implicit reflective and is
    interpreted by each solver's augmented mesh, e.g. ``SNMesh``
    defaults ``None`` → ``BC.reflective``).
    """
    case = UA_1_0_SP_STUB
    mesh = build_mesh(case, n_cells=24)
    assert mesh.coord == CoordSystem.SPHERICAL
    assert mesh.total_width == pytest.approx(_critical_dimension_cm(case), abs=1e-12)
    assert mesh.bc_left is None  # centreline implicit reflective
    assert mesh.bc_right == BC.vacuum


@pytest.mark.foundation
def test_cylinder_mesh_builds_radial_domain_with_reflective_axis() -> None:
    """Cylinder mesh covers ``[0, R]``; axis reflective is implicit."""
    case = LA13511_CASES["Ua-1-0-CY"]
    mesh = build_mesh(case, n_cells=20)
    assert mesh.coord == CoordSystem.CYLINDRICAL
    assert mesh.total_width == pytest.approx(_critical_dimension_cm(case), abs=1e-12)
    assert mesh.bc_left is None  # axis implicit reflective
    assert mesh.bc_right == BC.vacuum


@pytest.mark.foundation
def test_infinite_geometry_kind_raises_on_build() -> None:
    """build_mesh on infinite-medium case raises ValueError."""
    for case in (PUA_1_0_IN, PU_2_0_IN):
        with pytest.raises(ValueError, match="infinite"):
            build_mesh(case)


# ═══════════════════════════════════════════════════════════════════
# Foundation: registry-side production-protocol round-trip via
# transfer-matrix k_inf reference solver
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.parametrize("case", [PUA_1_0_IN, PU_2_0_IN], ids=lambda c: c.case_id)
def test_kinf_homogeneous_consumes_registry_case_directly(
    case: La13511Case,
) -> None:
    r"""The existing transfer-matrix
    :func:`orpheus.derivations.common.eigenvalue.kinf_homogeneous`
    consumes the registry's Mixture round-tripped to numpy arrays and
    reproduces Sood's published :math:`k_\infty` to ≤ 1e-5.

    This is a lightweight production-protocol smoke: if the Mixture
    round-trip is wrong, this test catches it without booting the
    full CP / SN solvers.
    """
    sigma_t, sigma_s, nu_sigma_f, chi = _xs_from_mixture(
        case.materials[0]
    )
    k = kinf_homogeneous(sigma_t, sigma_s, nu_sigma_f, chi)
    assert k == pytest.approx(case.truth.k_eff_or_kinf, abs=1e-5), (
        f"{case.case_id}: k_inf via registry+kinf_homogeneous = {k}, "
        f"truth = {case.truth.k_eff_or_kinf}"
    )


# ═══════════════════════════════════════════════════════════════════
# L1: production-solver structural-bridge smoke gates
# ═══════════════════════════════════════════════════════════════════
#
# These tests confirm the registry's case objects can be passed to
# production solvers (solve_cp, solve_sn) WITHOUT raising. They do
# NOT check answer correctness — production solvers were not
# originally tuned for the bare-critical Sood configurations and may
# need refinement / feature additions before they reach Sood
# precision. Phase B will land the full cross-checks.


@pytest.mark.l1
@pytest.mark.slow
def test_solve_cp_accepts_registry_sphere_case() -> None:
    r"""Bridge smoke: ``solve_cp`` accepts the registry's
    ``Ua-1-0-SP`` materials + mesh without raising.

    Asserts only that the call returns a result object with a finite
    ``keff`` field — this is a structural-bridge gate, NOT a
    correctness gate. Phase B will add the full cross-check.
    """
    import math

    from orpheus.cp.solver import CPParams, solve_cp

    case = UA_1_0_SP_STUB
    materials = build_materials(case)
    mesh = build_mesh(case, n_cells=32)
    params = CPParams(max_outer=50, keff_tol=1e-4, flux_tol=1e-3)

    result = solve_cp(materials, mesh, params)
    assert math.isfinite(result.keff), (
        f"solve_cp(Ua-1-0-SP) returned non-finite keff: {result.keff}"
    )
    assert result.keff > 0, (
        f"solve_cp(Ua-1-0-SP) returned non-positive keff: {result.keff}"
    )


@pytest.mark.l1
@pytest.mark.slow
def test_solve_sn_accepts_registry_slab_case() -> None:
    r"""Bridge smoke: ``solve_sn`` accepts the registry's ``Ua-1-0-SL``
    materials + mesh + a Gauss-Legendre quadrature without raising.

    Same structural-bridge gate as above. The mesh BCs are vacuum at
    both ends; SN consumes them via :class:`SNMesh` boundary
    handling.
    """
    import math

    from orpheus.sn.quadrature import GaussLegendre1D
    from orpheus.sn.solver import solve_sn

    case = UA_1_0_SL_STUB
    materials = build_materials(case)
    mesh = build_mesh(case, n_cells=32)
    quadrature = GaussLegendre1D.create(n_ordinates=16)

    result = solve_sn(
        materials, mesh, quadrature,
        max_outer=50, keff_tol=1e-4, flux_tol=1e-3,
    )
    assert math.isfinite(result.keff), (
        f"solve_sn(Ua-1-0-SL) returned non-finite keff: {result.keff}"
    )
    assert result.keff > 0, (
        f"solve_sn(Ua-1-0-SL) returned non-positive keff: {result.keff}"
    )
