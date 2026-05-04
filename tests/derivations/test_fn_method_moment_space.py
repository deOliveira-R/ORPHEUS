r"""Foundation tests for the F_N moment space (math-heart class).

These are software-invariant tests on the :class:`MomentSpace` facade
— bit-equality preservation against the function-level API,
production-protocol input acceptance, and result-type contract
preservation. They are NOT verification claims about the F_N method's
accuracy: those are owned by ``test_fn_la13511_*.py`` (which exercise
the function-level API and are unchanged by this dispatch).

The tests pin three invariants:

1. **Bit-equality with the function-level API** — the class produces
   IDENTICAL float results to direct calls of
   ``solve_fn_slab_bare_critical`` / ``solve_fn_sphere_bare_critical``
   / ``compute_kinf_*`` when given the same numerical inputs.
   Verified via ``float.hex(...)`` exact-bit comparison.

2. **Production-protocol input acceptance** — the class consumes the
   same ``(materials: dict[int, Mixture], GeometrySpec)`` pair as
   production CP/SN/MOC solvers, so a single problem definition
   serves both. Verified by constructing from
   :func:`orpheus.derivations.common.xs_library.make_mixture` +
   :class:`GeometrySpec`.

3. **Cross-method shared-result-type contract** — the class returns
   :class:`CriticalSolution` / :class:`FluxSolution` from
   :mod:`orpheus.derivations.common.solution_types`, with the
   correct ``parameter_kind`` / ``eigenvalue_kind`` / ``spatial_units``
   tags so cross-method consumers can compare results without method-
   specific knowledge.

Tests are foundation-tagged because they verify *structural* invariants
on the class facade, not L0/L1/L2 claims about a solver's accuracy —
the latter live in the LA-13511 test files that exercise the
underlying functions.
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
from orpheus.derivations.continuous.fn_method.moment_space import MomentSpace
from orpheus.derivations.continuous.fn_method.multi_group.k_inf import (
    compute_kinf_1g,
    compute_kinf_2g_general,
    compute_kinf_mg,
)
from orpheus.derivations.continuous.fn_method.slab import (
    solve_fn_slab_bare_critical,
)
from orpheus.derivations.continuous.fn_method.sphere import (
    solve_fn_sphere_bare_critical,
)
from orpheus.geometry.mesh import BC


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------


def _make_1g_mixture(sigma_t: float, sigma_s: float, nu_sigma_f: float):
    r"""Build a 1G ``Mixture`` with the cross sections needed for c-only F_N."""
    return make_mixture(
        sig_t=np.array([sigma_t]),
        sig_c=np.array([0.0]),
        sig_f=np.array([nu_sigma_f / 2.0]),  # nu=2 → nu*sig_f = nu_sigma_f
        nu=np.array([2.0]),
        chi=np.array([1.0]),
        sig_s=np.array([[sigma_s]]),
    )


def _ua10sl_geometry() -> GeometrySpec:
    """Sood Ua-1-0-SL slab geometry."""
    return GeometrySpec(
        geometry="slab",
        critical_dimension_mfp=0.93772556,
        critical_dimension_cm=0.93772556,
        n_groups=1,
        bc_left=BC.vacuum,
        bc_right=BC.vacuum,
    )


def _ua10sp_geometry() -> GeometrySpec:
    """Sood Ua-1-0-SP sphere geometry."""
    return GeometrySpec(
        geometry="sphere",
        critical_dimension_mfp=2.4248249802,
        critical_dimension_cm=2.4248249802,
        n_groups=1,
        bc_left=BC.reflective,
        bc_right=BC.vacuum,
    )


def _infinite_geometry() -> GeometrySpec:
    """Infinite-medium k_inf-only geometry."""
    return GeometrySpec(
        geometry="infinite",
        critical_dimension_mfp=None,
        critical_dimension_cm=None,
        n_groups=1,
    )


# ----------------------------------------------------------------------
# Construction and validation
# ----------------------------------------------------------------------


@pytest.mark.foundation
def test_moment_space_constructs_with_production_protocol_inputs() -> None:
    r"""``MomentSpace.from_problem`` accepts ``dict[int, Mixture] +
    GeometrySpec`` (the production-protocol input shape)."""
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    geom = _ua10sl_geometry()
    ms = MomentSpace.from_problem(materials={0: mix}, geometry=geom)
    assert isinstance(ms, MomentSpace)
    assert ms.geometry is geom
    assert ms.materials[0] is mix
    assert ms.fn_order == 9  # default
    assert ms.flux_reconstruction == "atkinson_nystrom"  # default


@pytest.mark.foundation
def test_moment_space_rejects_cylinder_geometry() -> None:
    r"""Cylinder is out of pillar (Westfall-Metcalf 1972); the class
    rejects it at construction time with an explicit error message."""
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    cyl = GeometrySpec(
        geometry="cylinder",
        critical_dimension_mfp=1.72500292,
        critical_dimension_cm=1.72500292,
        n_groups=1,
        bc_left=BC.reflective,
        bc_right=BC.vacuum,
    )
    with pytest.raises(ValueError, match="out of pillar"):
        MomentSpace.from_problem(materials={0: mix}, geometry=cyl)


@pytest.mark.foundation
def test_moment_space_rejects_unknown_flux_strategy() -> None:
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    geom = _ua10sl_geometry()
    with pytest.raises(ValueError, match="flux_reconstruction"):
        MomentSpace.from_problem(
            materials={0: mix},
            geometry=geom,
            flux_reconstruction="bogus",  # type: ignore[arg-type]
        )


@pytest.mark.foundation
def test_moment_space_c_property_matches_textbook_formula() -> None:
    r"""``MomentSpace.c`` returns
    :math:`(\Sigma_s + \nu\Sigma_f)/\Sigma_t` for 1G mixtures."""
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    ms = MomentSpace.from_problem(materials={0: mix}, geometry=_ua10sl_geometry())
    expected = (0.7 + 0.6) / 1.0
    assert ms.c == expected
    assert ms.n_groups == 1


# ----------------------------------------------------------------------
# Bit-equality preservation: solve_critical
# ----------------------------------------------------------------------


@pytest.mark.foundation
def test_slab_bit_equality_with_function_api() -> None:
    r"""Bit-equality: ``MomentSpace.solve_critical()`` returns the same
    ``a_critical_mfp`` as direct ``solve_fn_slab_bare_critical(c=ms.c)``.

    Verified by ``float.hex`` exact-bit comparison — the floats must
    be identical at the bit level, not merely close.
    """
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    ms = MomentSpace.from_problem(
        materials={0: mix}, geometry=_ua10sl_geometry(), fn_order=10
    )
    c_val = ms.c

    sol = ms.solve_critical()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = solve_fn_slab_bare_critical(c=c_val, n_modes=10)

    # Bit-equal float comparison via hex rep.
    assert sol.parameter_value.hex() == res.a_critical_mfp.hex(), (
        f"Bit-equality broken: class={sol.parameter_value!r} hex="
        f"{sol.parameter_value.hex()}; func={res.a_critical_mfp!r} "
        f"hex={res.a_critical_mfp.hex()}"
    )
    assert sol.parameter_kind == "half_thickness_mfp"
    assert sol.eigenvalue_kind == "k_eff"
    assert sol.eigenvalue == 1.0
    assert sol.metadata["raw_result"] is res or (
        # Also accept structural equality (the class re-runs the solve;
        # the SlabFNResult field-by-field equality is the load-bearing check).
        sol.metadata["raw_result"].a_critical_mfp == res.a_critical_mfp
    )


@pytest.mark.foundation
def test_sphere_bit_equality_with_function_api() -> None:
    r"""Bit-equality: ``MomentSpace.solve_critical()`` returns the same
    ``R_critical_mfp`` as direct ``solve_fn_sphere_bare_critical(c=ms.c)``."""
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    ms = MomentSpace.from_problem(
        materials={0: mix}, geometry=_ua10sp_geometry(), fn_order=10
    )
    c_val = ms.c

    sol = ms.solve_critical()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = solve_fn_sphere_bare_critical(c=c_val, n_modes=10)

    assert sol.parameter_value.hex() == res.R_critical_mfp.hex(), (
        f"Bit-equality broken: class={sol.parameter_value!r}; "
        f"func={res.R_critical_mfp!r}"
    )
    assert sol.parameter_kind == "radius_mfp"
    assert sol.eigenvalue_kind == "k_eff"
    assert sol.eigenvalue == 1.0


@pytest.mark.foundation
def test_kinf_1g_bit_equality_with_function_api() -> None:
    r"""1G k_inf path: ``MomentSpace.solve_critical()`` returns the same
    :math:`k_\infty` as direct ``compute_kinf_1g``."""
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    ms = MomentSpace.from_problem(materials={0: mix}, geometry=_infinite_geometry())

    sol = ms.solve_critical()
    expected = compute_kinf_1g(1.0, 0.7, 0.6)
    assert sol.eigenvalue.hex() == expected.hex()
    assert sol.eigenvalue_kind == "k_inf"
    assert sol.parameter_kind == "k_inf_only"
    assert sol.parameter_value == 0.0
    # solve_kinf shorthand returns the same value.
    assert ms.solve_kinf() == sol.eigenvalue


@pytest.mark.foundation
def test_kinf_2g_general_bit_equality() -> None:
    r"""2G k_inf with upscatter: ``MomentSpace.solve_critical()`` agrees
    bit-for-bit with ``compute_kinf_2g_general``.

    Uses Sood ``URRb-2-0-IN`` parameters (with-upscatter case).
    """
    sigma_t = np.array([0.65696, 1.52150])
    sigma_s = np.array([
        [0.62568, 0.029227],   # fast → fast, fast → slow
        [0.000767, 1.42710],   # slow → fast (upscatter), slow → slow
    ])
    nu_sigma_f = np.array([0.005305, 0.0976])
    chi = np.array([1.0, 0.0])
    mix = make_mixture(
        sig_t=sigma_t, sig_c=np.zeros(2),
        sig_f=nu_sigma_f / 2.0, nu=np.array([2.0, 2.0]),
        chi=chi, sig_s=sigma_s,
    )
    geom_inf = GeometrySpec(
        geometry="infinite", critical_dimension_mfp=None,
        critical_dimension_cm=None, n_groups=2,
    )
    ms = MomentSpace.from_problem(materials={0: mix}, geometry=geom_inf)
    sol = ms.solve_critical()
    expected = compute_kinf_2g_general(sigma_t, sigma_s, nu_sigma_f, chi)
    assert sol.eigenvalue.hex() == expected.hex()
    assert sol.eigenvalue_kind == "k_inf"
    assert sol.metadata["n_groups"] == 2
    assert sol.metadata["method"] == "compute_kinf_2g_general"


# ----------------------------------------------------------------------
# Cross-method shared-result-type contract
# ----------------------------------------------------------------------


@pytest.mark.foundation
def test_critical_solution_carries_correct_metadata_for_slab() -> None:
    """The :class:`CriticalSolution` returned by slab F_N populates the
    metadata fields that cross-method adapters consume."""
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    ms = MomentSpace.from_problem(
        materials={0: mix}, geometry=_ua10sl_geometry(), fn_order=8
    )
    sol = ms.solve_critical()
    assert isinstance(sol, CriticalSolution)
    assert "n_modes" in sol.metadata
    assert sol.metadata["n_modes"] == 8
    assert "determinant_residual" in sol.metadata
    assert "c" in sol.metadata
    assert "raw_result" in sol.metadata


@pytest.mark.foundation
def test_critical_solution_carries_correct_metadata_for_sphere() -> None:
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    ms = MomentSpace.from_problem(
        materials={0: mix}, geometry=_ua10sp_geometry(), fn_order=8
    )
    sol = ms.solve_critical()
    assert isinstance(sol, CriticalSolution)
    assert sol.metadata["n_modes"] == 8
    assert "determinant_residual" in sol.metadata
    assert "raw_result" in sol.metadata


# ----------------------------------------------------------------------
# reconstruct_flux: contract and end-to-end smoke
# ----------------------------------------------------------------------


@pytest.mark.foundation
def test_reconstruct_flux_none_strategy_raises() -> None:
    r"""Constructing with ``flux_reconstruction='none'`` makes
    :meth:`reconstruct_flux` raise — :meth:`solve_critical` is the
    only callable in that mode."""
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    ms = MomentSpace.from_problem(
        materials={0: mix}, geometry=_ua10sl_geometry(),
        flux_reconstruction="none",
    )
    # solve_critical still works.
    ms.solve_critical()
    # reconstruct_flux raises.
    with pytest.raises(ValueError, match="flux_reconstruction='none'"):
        ms.reconstruct_flux()


@pytest.mark.foundation
def test_reconstruct_flux_slab_returns_flux_solution() -> None:
    r"""End-to-end smoke: slab Atkinson reconstruction returns a
    :class:`FluxSolution` with the right shape and tags."""
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    ms = MomentSpace.from_problem(
        materials={0: mix}, geometry=_ua10sl_geometry(), fn_order=8,
    )
    flux = ms.reconstruct_flux(n_panels=64)
    assert isinstance(flux, FluxSolution)
    assert flux.spatial_nodes.shape == (64,)
    assert flux.scalar_flux.shape == (64,)
    assert flux.angular_flux is None
    assert flux.eigenvalue == 1.0
    assert flux.eigenvalue_kind == "k_eff"
    assert flux.spatial_units == "mfp"
    assert flux.metadata["geometry"] == "slab"
    assert flux.metadata["flux_reconstruction"] == "atkinson_nystrom"
    # Scalar flux must be positive everywhere on the eigenmode.
    assert np.all(flux.scalar_flux > 0)


@pytest.mark.foundation
def test_reconstruct_flux_sphere_returns_flux_solution() -> None:
    """End-to-end smoke: sphere KLL reconstruction returns a FluxSolution."""
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    ms = MomentSpace.from_problem(
        materials={0: mix}, geometry=_ua10sp_geometry(), fn_order=8,
    )
    flux = ms.reconstruct_flux(n_panels=32)
    assert isinstance(flux, FluxSolution)
    assert flux.spatial_nodes.shape == (32,)
    assert flux.scalar_flux.shape == (32,)
    assert flux.eigenvalue_kind == "k_eff"
    assert flux.spatial_units == "mfp"
    assert flux.metadata["geometry"] == "sphere"
    assert flux.metadata["flux_reconstruction"] == "kll_fredholm"
    # Flux is positive everywhere; centerline value is the maximum (sphere).
    assert np.all(flux.scalar_flux > 0)
    # Sphere flux is monotonically decreasing from r=0 outward (the
    # natural shape of the bare-critical eigenmode).
    assert flux.scalar_flux[0] >= flux.scalar_flux[-1]


# ----------------------------------------------------------------------
# Error handling
# ----------------------------------------------------------------------


@pytest.mark.foundation
def test_subcritical_slab_raises() -> None:
    r"""F_N bare-critical solve requires :math:`c > 1`; subcritical
    inputs are rejected at the boundary, not silently."""
    # c = (sigma_s + nu_sigma_f)/sigma_t = (0.5+0.4)/1.0 = 0.9 < 1
    mix = _make_1g_mixture(1.0, 0.5, 0.4)
    ms = MomentSpace.from_problem(
        materials={0: mix}, geometry=_ua10sl_geometry()
    )
    with pytest.raises(ValueError, match="multiplying medium"):
        ms.solve_critical()
