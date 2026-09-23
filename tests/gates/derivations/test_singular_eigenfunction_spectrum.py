r"""Foundation tests for the Case singular-eigenfunction spectrum (Phase D).

These are software-invariant tests on the
:class:`~orpheus.derivations.continuous.singular_eigenfunction.spectrum.Spectrum`
facade — bit-equality preservation against the function-level API,
StructuredGeometry input acceptance, and result-type contract
preservation.

Phase D
-------

Spectrum now consumes a :class:`StructuredGeometry` directly via
its frozen ``__init__`` (no GeometrySpec, no from_problem factory).
Cylinder + linear-anisotropy validation moved to ``__post_init__``.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.solution_types import (
    CriticalSolution,
    FluxSolution,
)
from orpheus.derivations.common.xs_library import make_mixture
from orpheus.derivations.continuous.singular_eigenfunction import Spectrum
from orpheus.derivations.continuous.singular_eigenfunction.cylinder.one_group import (
    solve_singular_eigenfunction_cylinder_bare_critical,
)
from orpheus.derivations.continuous.singular_eigenfunction.slab.one_group import (
    solve_case_method_slab_critical,
)
from orpheus.derivations.continuous.singular_eigenfunction.sphere.one_group import (
    solve_case_method_sphere_critical,
)
from orpheus.geometry.mesh import BC
from orpheus.geometry.structured_geometry import (
    Region,
    StructuredGeometry,
)


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------


def _make_1g_mixture(
    sigma_t: float, sigma_s: float, nu_sigma_f: float
):
    r"""Build an isotropic 1G ``Mixture`` for the c-only Spectrum tests."""
    return make_mixture(
        sig_t=np.array([sigma_t]),
        sig_c=np.array([0.0]),
        sig_f=np.array([nu_sigma_f / 2.0]),
        nu=np.array([2.0]),
        chi=np.array([1.0]),
        sig_s=np.array([[sigma_s]]),
    )


def _make_1g_anisotropic_mixture(
    sigma_t: float, sigma_s: float, nu_sigma_f: float, f1: float
):
    r"""Build a linearly-anisotropic 1G ``Mixture``."""
    return make_mixture(
        sig_t=np.array([sigma_t]),
        sig_c=np.array([0.0]),
        sig_f=np.array([nu_sigma_f / 2.0]),
        nu=np.array([2.0]),
        chi=np.array([1.0]),
        sig_s=np.array([[sigma_s]]),
        sig_s1=np.array([[f1 * sigma_s]]),
    )


def _slab_geometry(R_albedo: float = 0.0) -> StructuredGeometry:
    r"""Slab geometry with reflection coefficient on the right (outer) BC."""
    if R_albedo == 0.0:
        outer_bc: BC = BC.vacuum
    else:
        outer_bc = BC("partial", {"albedo": R_albedo})
    return StructuredGeometry(
        geometry="SLB",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0 * 0.93772556),),
        bcs=(BC.vacuum, outer_bc),
    )


def _sphere_geometry(R_albedo: float = 0.0) -> StructuredGeometry:
    r"""Sphere geometry with reflection coefficient on the outer BC."""
    if R_albedo == 0.0:
        outer_bc: BC = BC.vacuum
    else:
        outer_bc = BC("partial", {"albedo": R_albedo})
    return StructuredGeometry(
        geometry="SPH",
        regions=(Region(mat_id=0, outer_thickness_cm=2.4248249802),),
        bcs=(outer_bc,),
    )


def _cylinder_geometry() -> StructuredGeometry:
    r"""Bare cylinder geometry (Sood ``Ua-1-0-CY``)."""
    return StructuredGeometry(
        geometry="CYL",
        regions=(Region(mat_id=0, outer_thickness_cm=1.72500292),),
        bcs=(BC.vacuum,),
    )


# ----------------------------------------------------------------------
# Construction and validation
# ----------------------------------------------------------------------


@pytest.mark.foundation
def test_spectrum_constructs_with_structured_geometry() -> None:
    r"""``Spectrum`` accepts ``StructuredGeometry`` + ``dict[int, Mixture]``."""
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    geom = _slab_geometry()
    spec = Spectrum(geometry=geom, materials={0: mix})
    assert isinstance(spec, Spectrum)
    assert spec.geometry is geom
    assert spec.materials[0] is mix
    assert spec.n_modes == 8  # default


@pytest.mark.foundation
def test_spectrum_rejects_anisotropic_cylinder() -> None:
    r"""Cylinder + linear anisotropy is out of pillar — Westfall-Metcalf
    1972 covers isotropic only. Validation runs in __post_init__."""
    mix = _make_1g_anisotropic_mixture(1.0, 0.7, 0.6, f1=0.20)
    geom = _cylinder_geometry()
    with pytest.raises(NotImplementedError, match="cylinder \\+ linear"):
        Spectrum(geometry=geom, materials={0: mix})


@pytest.mark.foundation
def test_spectrum_c_property_isotropic_matches_textbook() -> None:
    r"""``Spectrum.c`` returns :math:`(\Sigma_s + \nu\Sigma_f)/\Sigma_t`."""
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    spec = Spectrum(geometry=_slab_geometry(), materials={0: mix})
    expected = (0.7 + 0.6) / 1.0
    assert spec.c == expected
    assert spec.n_groups == 1
    assert spec.f1 == 0.0


@pytest.mark.foundation
def test_spectrum_f1_property_extracts_from_sigs1() -> None:
    r"""``Spectrum.f1`` returns :math:`\Sigma_{s,1}/\Sigma_{s,0}`."""
    mix = _make_1g_anisotropic_mixture(1.0, 0.7, 0.6, f1=0.25)
    spec = Spectrum(geometry=_slab_geometry(), materials={0: mix})
    assert abs(spec.f1 - 0.25) < 1e-15


# ----------------------------------------------------------------------
# Bit-equality preservation: solve_critical
# ----------------------------------------------------------------------


@pytest.mark.foundation
def test_cylinder_bit_equality_with_function_api() -> None:
    r"""Bit-equality: ``Spectrum.solve_critical()`` returns the same
    ``r_c_mfp`` as direct
    ``solve_singular_eigenfunction_cylinder_bare_critical(c=spec.c,
    n_grid=spec.n_modes)``."""
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    spec = Spectrum(
        geometry=_cylinder_geometry(), materials={0: mix}, n_modes=24
    )
    c_extracted = spec.c

    sol = spec.solve_critical()
    res = solve_singular_eigenfunction_cylinder_bare_critical(
        c=c_extracted, n_grid=24
    )

    assert sol.parameter_value.hex() == res.r_c_mfp.hex(), (
        f"Bit-equality broken: class={sol.parameter_value!r} hex="
        f"{sol.parameter_value.hex()}; func={res.r_c_mfp!r} hex="
        f"{res.r_c_mfp.hex()}"
    )
    assert sol.parameter_kind == "radius_mfp"
    assert sol.eigenvalue_kind == "k_eff"
    assert sol.eigenvalue == 1.0
    assert sol.metadata["raw_result"] is res or (
        sol.metadata["raw_result"].r_c_mfp == res.r_c_mfp
    )


@pytest.mark.foundation
def test_slab_bit_equality_with_function_api() -> None:
    r"""Bit-equality: ``Spectrum.solve_critical()`` for slab returns the
    same ``d_critical_mfp`` as direct
    ``solve_case_method_slab_critical``."""
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    spec = Spectrum(
        geometry=_slab_geometry(), materials={0: mix}, n_modes=8
    )
    c_extracted = spec.c

    sol = spec.solve_critical()
    res = solve_case_method_slab_critical(
        c=c_extracted, R=0.0, f1=0.0, maxdegree=8
    )
    assert sol.parameter_value.hex() == res.d_critical_mfp.hex()
    assert sol.parameter_kind == "half_thickness_mfp"
    assert sol.eigenvalue_kind == "k_eff"


@pytest.mark.foundation
def test_sphere_bit_equality_with_function_api() -> None:
    r"""Bit-equality: ``Spectrum.solve_critical()`` for sphere returns the
    same ``R_critical_mfp`` as direct
    ``solve_case_method_sphere_critical``."""
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    spec = Spectrum(
        geometry=_sphere_geometry(), materials={0: mix}, n_modes=8
    )
    c_extracted = spec.c

    sol = spec.solve_critical()
    res = solve_case_method_sphere_critical(
        c=c_extracted, R_refl=0.0, f1=0.0
    )
    assert sol.parameter_value.hex() == res.R_critical_mfp.hex()
    assert sol.parameter_kind == "radius_mfp"
    assert sol.eigenvalue_kind == "k_eff"


@pytest.mark.foundation
def test_slab_anisotropic_bit_equality_with_function_api() -> None:
    r"""Bit-equality with linear anisotropy: Spectrum extracts ``f1``
    from the Mixture and forwards it to the underlying solver."""
    mix = _make_1g_anisotropic_mixture(1.0, 0.7, 0.6, f1=0.20)
    spec = Spectrum(
        geometry=_slab_geometry(), materials={0: mix}, n_modes=8
    )

    f1_extracted = spec.f1
    sol = spec.solve_critical()
    res = solve_case_method_slab_critical(
        c=spec.c, R=0.0, f1=f1_extracted, maxdegree=8
    )
    assert sol.parameter_value.hex() == res.d_critical_mfp.hex()
    assert sol.metadata["f1"] == f1_extracted


@pytest.mark.foundation
def test_slab_partial_reflection_bc_routing() -> None:
    r"""``BC("partial", {"albedo": R})`` flows through to the underlying
    solver's ``R`` parameter, returning the reflected-slab critical
    half-thickness."""
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    geom = _slab_geometry(R_albedo=0.50)
    spec = Spectrum(geometry=geom, materials={0: mix}, n_modes=8)

    sol = spec.solve_critical()
    res = solve_case_method_slab_critical(
        c=spec.c, R=0.50, f1=0.0, maxdegree=8
    )
    assert sol.parameter_value.hex() == res.d_critical_mfp.hex()
    assert sol.metadata["R"] == 0.50


# ----------------------------------------------------------------------
# Cross-method shared-result-type contract
# ----------------------------------------------------------------------


@pytest.mark.foundation
def test_critical_solution_is_correct_type() -> None:
    r"""Every ``solve_critical()`` call returns the SHARED
    :class:`CriticalSolution` type."""
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    spec = Spectrum(
        geometry=_cylinder_geometry(), materials={0: mix}, n_modes=24
    )
    sol = spec.solve_critical()
    assert isinstance(sol, CriticalSolution)
    assert sol.eigenvalue_kind == "k_eff"
    assert sol.eigenvalue == 1.0
    assert sol.parameter_kind == "radius_mfp"
    assert sol.parameter_value > 0.0
    assert sol.converged is True
    assert sol.metadata["method"] == (
        "solve_singular_eigenfunction_cylinder_bare_critical"
    )
    assert "raw_result" in sol.metadata
    assert "u_0" in sol.metadata


@pytest.mark.foundation
def test_solve_fixed_source_returns_flux_solution_for_cylinder() -> None:
    r"""Cylinder ``solve_fixed_source`` returns the Westfall-Metcalf
    bare-cylinder dominant Case eigenfunction wrapped as a
    :class:`FluxSolution`."""
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    spec = Spectrum(
        geometry=_cylinder_geometry(), materials={0: mix}, n_modes=24
    )
    flux = spec.solve_fixed_source(n_eval=32)
    assert isinstance(flux, FluxSolution)
    assert flux.spatial_units == "mfp"
    assert flux.eigenvalue == 1.0
    assert flux.eigenvalue_kind == "k_eff"
    assert flux.scalar_flux.shape == (32,)
    assert abs(flux.scalar_flux[0] - 1.0) < 1e-15
    assert np.all(np.diff(flux.scalar_flux) <= 0.0)


@pytest.mark.foundation
def test_solve_fixed_source_rejects_slab() -> None:
    r"""Slab flux reconstruction is owned by the F_N pillar."""
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    spec = Spectrum(
        geometry=_slab_geometry(), materials={0: mix}, n_modes=8
    )
    with pytest.raises(NotImplementedError, match="F_N pillar"):
        spec.solve_fixed_source()


@pytest.mark.foundation
def test_solve_fixed_source_rejects_sphere() -> None:
    r"""Sphere flux reconstruction is owned by the F_N pillar (KLL)."""
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    spec = Spectrum(
        geometry=_sphere_geometry(), materials={0: mix}, n_modes=8
    )
    with pytest.raises(NotImplementedError, match="F_N pillar"):
        spec.solve_fixed_source()


# ----------------------------------------------------------------------
# Validation guards
# ----------------------------------------------------------------------


@pytest.mark.foundation
def test_spectrum_rejects_n_modes_too_small() -> None:
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    with pytest.raises(ValueError, match="n_modes"):
        Spectrum(
            geometry=_slab_geometry(), materials={0: mix}, n_modes=1
        )


@pytest.mark.foundation
def test_spectrum_rejects_missing_mat_id() -> None:
    """Region's mat_id must be present in the materials dict."""
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    geom = StructuredGeometry(
        geometry="SLB",
        regions=(Region(mat_id=42, outer_thickness_cm=1.88),),
        bcs=(BC.vacuum, BC.vacuum),
    )
    with pytest.raises(ValueError, match="mat_id=42"):
        Spectrum(geometry=geom, materials={0: mix})


@pytest.mark.foundation
def test_spectrum_rejects_subcritical_medium() -> None:
    r"""``c <= 1`` (sub-multiplying medium) has no real bare-critical
    configuration — Spectrum rejects it explicitly."""
    mix = _make_1g_mixture(1.0, 0.7, 0.2)  # c = 0.9 < 1
    spec = Spectrum(
        geometry=_cylinder_geometry(), materials={0: mix}, n_modes=24
    )
    with pytest.raises(ValueError, match="c > 1"):
        spec.solve_critical()
