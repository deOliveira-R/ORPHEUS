r"""Foundation tests for the Case singular-eigenfunction spectrum (math-heart class).

These are software-invariant tests on the
:class:`~orpheus.derivations.continuous.singular_eigenfunction.spectrum.Spectrum`
facade — bit-equality preservation against the function-level API,
production-protocol input acceptance, and result-type contract
preservation. They are NOT verification claims about the
singular-eigenfunction method's accuracy: those are owned by the
existing ``test_singular_eigenfunction_cylinder.py``,
``test_case_method_slab.py``, ``test_case_method_sphere.py`` test
files (which exercise the function-level API and are unchanged by
this dispatch).

The tests pin three invariants:

1. **Bit-equality with the function-level API** — the class produces
   IDENTICAL float results to direct calls of
   ``solve_case_method_slab_critical`` /
   ``solve_case_method_sphere_critical`` /
   ``solve_singular_eigenfunction_cylinder_bare_critical`` when given
   the same numerical inputs. Verified via ``float.hex(...)``
   exact-bit comparison.

2. **Production-protocol input acceptance** — the class consumes the
   same ``(materials: dict[int, Mixture], GeometrySpec)`` pair as
   production CP/SN/MOC solvers and the parallel ``MomentSpace`` /
   ``Billiard`` math-heart classes. Verified by constructing from
   :func:`orpheus.derivations.common.xs_library.make_mixture` +
   :class:`GeometrySpec`.

3. **Cross-method shared-result-type contract** — the class returns
   :class:`~orpheus.derivations.common.solution_types.CriticalSolution`
   /
   :class:`~orpheus.derivations.common.solution_types.FluxSolution`
   from :mod:`orpheus.derivations.common.solution_types`, with the
   correct ``parameter_kind`` / ``eigenvalue_kind`` / ``spatial_units``
   tags so cross-method consumers can compare results without
   method-specific knowledge.

The third sibling
=================

These tests pin :class:`Spectrum` as the **3rd concrete instance of
the math-heart pattern** — alongside ``MomentSpace`` (F_N method,
Galerkin moment projection) and ``Billiard`` (trajectory_resolvent,
Birkhoff transfer-operator resolvent). The "≥3 instances" threshold
empirically validates the unifying ``TransportSolver`` Protocol
designed in the parallel agent's spec — Spectrum's structural
conformance (factory + properties + Protocol-shaped solve methods)
is part of the cross-method contract verification.

Tests are foundation-tagged (``@pytest.mark.foundation``) because they
verify *structural* invariants on the class facade, not L0/L1/L2
claims about a solver's accuracy.
"""
from __future__ import annotations

import warnings

import numpy as np
import pytest

from orpheus.derivations.common.geometry_spec import GeometrySpec
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
    r"""Build a linearly-anisotropic 1G ``Mixture`` for Atalay slab/sphere tests.

    Convention: :math:`\Sigma_{s,1} = f_1\,\Sigma_{s,0}` so the
    extracted ``f1 = SigS[1] / SigS[0]`` returns the requested value.
    """
    return make_mixture(
        sig_t=np.array([sigma_t]),
        sig_c=np.array([0.0]),
        sig_f=np.array([nu_sigma_f / 2.0]),
        nu=np.array([2.0]),
        chi=np.array([1.0]),
        sig_s=np.array([[sigma_s]]),
        sig_s1=np.array([[f1 * sigma_s]]),
    )


def _slab_geometry(R_albedo: float = 0.0) -> GeometrySpec:
    r"""Slab geometry with reflection coefficient ``R_albedo`` on bc_right."""
    if R_albedo == 0.0:
        bc_right: BC = BC.vacuum
    else:
        bc_right = BC("partial", {"albedo": R_albedo})
    return GeometrySpec(
        geometry="slab",
        critical_dimension_mfp=0.93772556,
        critical_dimension_cm=0.93772556,
        n_groups=1,
        bc_left=BC.vacuum,
        bc_right=bc_right,
    )


def _sphere_geometry(R_albedo: float = 0.0) -> GeometrySpec:
    r"""Sphere geometry with reflection coefficient ``R_albedo`` on bc_right."""
    if R_albedo == 0.0:
        bc_right: BC = BC.vacuum
    else:
        bc_right = BC("partial", {"albedo": R_albedo})
    return GeometrySpec(
        geometry="sphere",
        critical_dimension_mfp=2.4248249802,
        critical_dimension_cm=2.4248249802,
        n_groups=1,
        bc_left=BC.reflective,
        bc_right=bc_right,
    )


def _cylinder_geometry() -> GeometrySpec:
    r"""Bare cylinder geometry (Sood ``Ua-1-0-CY``)."""
    return GeometrySpec(
        geometry="cylinder",
        critical_dimension_mfp=1.72500292,
        critical_dimension_cm=1.72500292,
        n_groups=1,
        bc_left=BC.reflective,
        bc_right=BC.vacuum,
    )


# ----------------------------------------------------------------------
# Construction and validation
# ----------------------------------------------------------------------


@pytest.mark.foundation
def test_spectrum_constructs_with_production_protocol_inputs() -> None:
    r"""``Spectrum.from_problem`` accepts ``dict[int, Mixture] +
    GeometrySpec`` (the production-protocol input shape — same as
    ``MomentSpace`` and ``Billiard``)."""
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    geom = _slab_geometry()
    spec = Spectrum.from_problem(materials={0: mix}, geometry=geom)
    assert isinstance(spec, Spectrum)
    assert spec.geometry is geom
    assert spec.materials[0] is mix
    assert spec.n_modes == 8  # default


@pytest.mark.foundation
def test_spectrum_rejects_infinite_geometry() -> None:
    r"""Infinite-medium k_inf is out of pillar — singular-eigenfunction
    criticality requires a finite spatial domain. ``MomentSpace`` is
    the right tool for k_inf calls."""
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    inf_geom = GeometrySpec(
        geometry="infinite",
        critical_dimension_mfp=None,
        critical_dimension_cm=None,
        n_groups=1,
    )
    with pytest.raises(ValueError, match="out of pillar"):
        Spectrum.from_problem(materials={0: mix}, geometry=inf_geom)


@pytest.mark.foundation
def test_spectrum_rejects_anisotropic_cylinder() -> None:
    r"""Cylinder + linear anisotropy is out of pillar — Westfall-Metcalf
    1972 covers isotropic only."""
    mix = _make_1g_anisotropic_mixture(1.0, 0.7, 0.6, f1=0.20)
    geom = _cylinder_geometry()
    with pytest.raises(NotImplementedError, match="cylinder \\+ linear"):
        Spectrum.from_problem(materials={0: mix}, geometry=geom)


@pytest.mark.foundation
def test_spectrum_protocol_method_name_is_canonical() -> None:
    r"""``Spectrum.method_name == "singular_eigenfunction"`` — used by
    cross-method adapters to dispatch on method type without an
    isinstance ladder."""
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    spec = Spectrum.from_problem(
        materials={0: mix}, geometry=_slab_geometry()
    )
    assert spec.method_name == "singular_eigenfunction"


@pytest.mark.foundation
def test_spectrum_protocol_geometry_spec_alias() -> None:
    r"""``Spectrum.geometry_spec`` is an alias for the ``geometry`` field
    — Protocol-conforming surface (the Protocol contract names this
    property ``geometry_spec``; the dataclass field is ``geometry`` to
    match :class:`MomentSpace`)."""
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    geom = _slab_geometry()
    spec = Spectrum.from_problem(materials={0: mix}, geometry=geom)
    assert spec.geometry_spec is spec.geometry
    assert spec.geometry_spec is geom


@pytest.mark.foundation
def test_spectrum_c_property_isotropic_matches_textbook() -> None:
    r"""``Spectrum.c`` returns :math:`(\Sigma_s + \nu\Sigma_f)/\Sigma_t`
    for isotropic 1G mixtures."""
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    spec = Spectrum.from_problem(
        materials={0: mix}, geometry=_slab_geometry()
    )
    expected = (0.7 + 0.6) / 1.0
    assert spec.c == expected
    assert spec.n_groups == 1
    assert spec.f1 == 0.0


@pytest.mark.foundation
def test_spectrum_f1_property_extracts_from_sigs1() -> None:
    r"""``Spectrum.f1`` returns :math:`\Sigma_{s,1}/\Sigma_{s,0}` from
    a linearly-anisotropic mixture's ``SigS[1] / SigS[0]``."""
    mix = _make_1g_anisotropic_mixture(1.0, 0.7, 0.6, f1=0.25)
    spec = Spectrum.from_problem(
        materials={0: mix}, geometry=_slab_geometry()
    )
    # Floating-point exact: 0.25 * 0.7 / 0.7 = 0.25
    assert abs(spec.f1 - 0.25) < 1e-15


# ----------------------------------------------------------------------
# Bit-equality preservation: solve_critical
# ----------------------------------------------------------------------


@pytest.mark.foundation
def test_cylinder_bit_equality_with_function_api() -> None:
    r"""Bit-equality: ``Spectrum.solve_critical()`` returns the same
    ``r_c_mfp`` as direct
    ``solve_singular_eigenfunction_cylinder_bare_critical(c=spec.c,
    n_grid=spec.n_modes)``.

    Verified by ``float.hex`` exact-bit comparison.
    """
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    spec = Spectrum.from_problem(
        materials={0: mix}, geometry=_cylinder_geometry(), n_modes=24
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
    ``solve_case_method_slab_critical(c=spec.c, R=0.0, f1=0.0,
    maxdegree=spec.n_modes)``."""
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    spec = Spectrum.from_problem(
        materials={0: mix}, geometry=_slab_geometry(), n_modes=8
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
    ``solve_case_method_sphere_critical(c=spec.c, R_refl=0.0, f1=0.0)``."""
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    spec = Spectrum.from_problem(
        materials={0: mix}, geometry=_sphere_geometry(), n_modes=8
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
    from the Mixture and forwards it to the underlying solver.

    Sood-equivalent c=1.30 with :math:`f_1 = 0.20` (within Atalay's
    Eq 5 validity bound :math:`c \le 1 + 1/(3 f_1) = 2.667`)."""
    mix = _make_1g_anisotropic_mixture(1.0, 0.7, 0.6, f1=0.20)
    spec = Spectrum.from_problem(
        materials={0: mix}, geometry=_slab_geometry(), n_modes=8
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
    half-thickness.  Atalay 1997 Table 2 R=0.50 column for c=1.30
    isotropic."""
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    geom = _slab_geometry(R_albedo=0.50)
    spec = Spectrum.from_problem(materials={0: mix}, geometry=geom, n_modes=8)

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
    :class:`CriticalSolution` type from
    :mod:`orpheus.derivations.common.solution_types` — making the
    output substitutable across MomentSpace / Billiard / Spectrum."""
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    spec = Spectrum.from_problem(
        materials={0: mix}, geometry=_cylinder_geometry(), n_modes=24
    )
    sol = spec.solve_critical()
    assert isinstance(sol, CriticalSolution)
    assert sol.eigenvalue_kind == "k_eff"
    assert sol.eigenvalue == 1.0
    assert sol.parameter_kind == "radius_mfp"
    assert sol.parameter_value > 0.0
    assert sol.converged is True
    # Method-specific metadata must be present.
    assert sol.metadata["method"] == (
        "solve_singular_eigenfunction_cylinder_bare_critical"
    )
    assert "raw_result" in sol.metadata
    assert "u_0" in sol.metadata


@pytest.mark.foundation
def test_solve_fixed_source_returns_flux_solution_for_cylinder() -> None:
    r"""Cylinder ``solve_fixed_source`` returns the Westfall-Metcalf
    bare-cylinder dominant Case eigenfunction
    :math:`\rho(r) = J_0(r/u_0)` wrapped as a :class:`FluxSolution`."""
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    spec = Spectrum.from_problem(
        materials={0: mix}, geometry=_cylinder_geometry(), n_modes=24
    )
    flux = spec.solve_fixed_source(n_eval=32)
    assert isinstance(flux, FluxSolution)
    assert flux.spatial_units == "mfp"
    assert flux.eigenvalue == 1.0
    assert flux.eigenvalue_kind == "k_eff"
    assert flux.scalar_flux.shape == (32,)
    # rho(0) = J_0(0) = 1.
    assert abs(flux.scalar_flux[0] - 1.0) < 1e-15
    # Monotone-decreasing on (0, R_c).
    assert np.all(np.diff(flux.scalar_flux) <= 0.0)


@pytest.mark.foundation
def test_solve_fixed_source_rejects_slab() -> None:
    r"""Slab flux reconstruction is owned by the F_N pillar — using
    Spectrum for it would violate structural independence."""
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    spec = Spectrum.from_problem(
        materials={0: mix}, geometry=_slab_geometry(), n_modes=8
    )
    with pytest.raises(NotImplementedError, match="F_N pillar"):
        spec.solve_fixed_source()


@pytest.mark.foundation
def test_solve_fixed_source_rejects_sphere() -> None:
    r"""Sphere flux reconstruction is owned by the F_N pillar (KLL
    Fredholm iteration)."""
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    spec = Spectrum.from_problem(
        materials={0: mix}, geometry=_sphere_geometry(), n_modes=8
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
        Spectrum.from_problem(
            materials={0: mix}, geometry=_slab_geometry(), n_modes=1
        )


@pytest.mark.foundation
def test_spectrum_rejects_missing_mat_id() -> None:
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    geom = GeometrySpec(
        geometry="slab",
        critical_dimension_mfp=0.94,
        critical_dimension_cm=0.94,
        n_groups=1,
        mat_id=42,  # not present in materials dict
        bc_left=BC.vacuum,
        bc_right=BC.vacuum,
    )
    with pytest.raises(ValueError, match="mat_id=42"):
        Spectrum.from_problem(materials={0: mix}, geometry=geom)


@pytest.mark.foundation
def test_spectrum_rejects_subcritical_medium() -> None:
    r"""``c <= 1`` (sub-multiplying medium) has no real bare-critical
    configuration — Spectrum rejects it explicitly."""
    mix = _make_1g_mixture(1.0, 0.7, 0.2)  # c = 0.9 < 1
    spec = Spectrum.from_problem(
        materials={0: mix}, geometry=_cylinder_geometry(), n_modes=24
    )
    with pytest.raises(ValueError, match="c > 1"):
        spec.solve_critical()
