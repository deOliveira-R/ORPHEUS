r"""Foundation tests for the F_N moment space (math-heart class, Phase D).

These are software-invariant tests on the :class:`MomentSpace` facade
— bit-equality preservation against the function-level API,
StructuredGeometry input acceptance, and result-type contract
preservation. They are NOT verification claims about the F_N method's
accuracy: those are owned by ``test_fn_la13511_*.py`` (which exercise
the function-level API and are unchanged by this dispatch).

Phase D consumes :class:`StructuredGeometry` directly via
``MomentSpace(geometry=..., materials=...)``; infinite-medium
:math:`k_\infty` is a separate static method
``MomentSpace.solve_kinf(mixture)`` (no geometry needed).
"""
from __future__ import annotations

import warnings

import numpy as np
import pytest

from orpheus.derivations.common.solution_types import (
    CriticalSolution,
    FluxSolution,
)
from orpheus.derivations.common.xs_library import make_mixture
from orpheus.derivations.continuous.fn_method.moment_space import MomentSpace
from orpheus.derivations.continuous.fn_method.multi_group.k_inf import (
    compute_kinf_1g,
    compute_kinf_2g_general,
)
from orpheus.derivations.continuous.fn_method.slab import (
    solve_fn_slab_bare_critical,
)
from orpheus.derivations.continuous.fn_method.sphere import (
    solve_fn_sphere_bare_critical,
)
from orpheus.geometry.mesh import BC
from orpheus.geometry.structured_geometry import (
    Region,
    StructuredGeometry,
)


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


def _ua10sl_geometry() -> StructuredGeometry:
    """Sood Ua-1-0-SL slab geometry (full slab width = 2 × half-thickness)."""
    return StructuredGeometry(
        geometry="SLB",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0 * 0.93772556),),
        bcs=(BC.vacuum, BC.vacuum),
    )


def _ua10sp_geometry() -> StructuredGeometry:
    """Sood Ua-1-0-SP sphere geometry."""
    return StructuredGeometry(
        geometry="SPH",
        regions=(Region(mat_id=0, outer_thickness_cm=2.4248249802),),
        bcs=(BC.vacuum,),
    )


def _cylinder_geometry() -> StructuredGeometry:
    """Sood-style cylinder geometry (out-of-pillar for F_N)."""
    return StructuredGeometry(
        geometry="CYL",
        regions=(Region(mat_id=0, outer_thickness_cm=1.72500292),),
        bcs=(BC.vacuum,),
    )


# ----------------------------------------------------------------------
# Construction and validation
# ----------------------------------------------------------------------


@pytest.mark.foundation
def test_moment_space_constructs_with_structured_geometry() -> None:
    r"""``MomentSpace`` accepts ``StructuredGeometry`` + ``dict[int, Mixture]``."""
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    geom = _ua10sl_geometry()
    ms = MomentSpace(geometry=geom, materials={0: mix})
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
    cyl = _cylinder_geometry()
    with pytest.raises(ValueError, match="out of pillar"):
        MomentSpace(geometry=cyl, materials={0: mix})


@pytest.mark.foundation
def test_moment_space_rejects_unknown_flux_strategy() -> None:
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    geom = _ua10sl_geometry()
    with pytest.raises(ValueError, match="flux_reconstruction"):
        MomentSpace(
            geometry=geom,
            materials={0: mix},
            flux_reconstruction="bogus",  # type: ignore[arg-type]
        )


@pytest.mark.foundation
def test_moment_space_c_property_matches_textbook_formula() -> None:
    r"""``MomentSpace.c`` returns
    :math:`(\Sigma_s + \nu\Sigma_f)/\Sigma_t` for 1G mixtures."""
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    ms = MomentSpace(geometry=_ua10sl_geometry(), materials={0: mix})
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
    """
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    ms = MomentSpace(
        geometry=_ua10sl_geometry(), materials={0: mix}, fn_order=10
    )
    c_val = ms.c

    sol = ms.solve_critical()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = solve_fn_slab_bare_critical(c=c_val, n_modes=10)

    assert sol.parameter_value.hex() == res.a_critical_mfp.hex(), (
        f"Bit-equality broken: class={sol.parameter_value!r} hex="
        f"{sol.parameter_value.hex()}; func={res.a_critical_mfp!r} "
        f"hex={res.a_critical_mfp.hex()}"
    )
    assert sol.parameter_kind == "half_thickness_mfp"
    assert sol.eigenvalue_kind == "k_eff"
    assert sol.eigenvalue == 1.0
    assert sol.metadata["raw_result"] is res or (
        sol.metadata["raw_result"].a_critical_mfp == res.a_critical_mfp
    )


@pytest.mark.foundation
def test_sphere_bit_equality_with_function_api() -> None:
    r"""Bit-equality: ``MomentSpace.solve_critical()`` returns the same
    ``R_critical_mfp`` as direct ``solve_fn_sphere_bare_critical(c=ms.c)``."""
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    ms = MomentSpace(
        geometry=_ua10sp_geometry(), materials={0: mix}, fn_order=10
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
    r"""1G k_inf path: ``MomentSpace.solve_kinf(mix)`` (no-geometry
    static method) returns the same :math:`k_\infty` as direct
    ``compute_kinf_1g``."""
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    sol = MomentSpace.solve_kinf(mix)
    expected = compute_kinf_1g(1.0, 0.7, 0.6)
    assert sol.eigenvalue.hex() == expected.hex()
    assert sol.eigenvalue_kind == "k_inf"
    assert sol.parameter_kind == "k_inf_only"
    assert sol.parameter_value == 0.0


@pytest.mark.foundation
def test_kinf_2g_general_bit_equality() -> None:
    r"""2G k_inf with upscatter: ``MomentSpace.solve_kinf(mix)`` agrees
    bit-for-bit with ``compute_kinf_2g_general``.

    Uses Sood ``URRb-2-0-IN`` parameters (with-upscatter case).
    """
    sigma_t = np.array([0.65696, 1.52150])
    sigma_s = np.array([
        [0.62568, 0.029227],
        [0.000767, 1.42710],
    ])
    nu_sigma_f = np.array([0.005305, 0.0976])
    chi = np.array([1.0, 0.0])
    mix = make_mixture(
        sig_t=sigma_t, sig_c=np.zeros(2),
        sig_f=nu_sigma_f / 2.0, nu=np.array([2.0, 2.0]),
        chi=chi, sig_s=sigma_s,
    )
    sol = MomentSpace.solve_kinf(mix)
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
    ms = MomentSpace(
        geometry=_ua10sl_geometry(), materials={0: mix}, fn_order=8
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
    ms = MomentSpace(
        geometry=_ua10sp_geometry(), materials={0: mix}, fn_order=8
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
    ms = MomentSpace(
        geometry=_ua10sl_geometry(), materials={0: mix},
        flux_reconstruction="none",
    )
    ms.solve_critical()
    with pytest.raises(ValueError, match="flux_reconstruction='none'"):
        ms.reconstruct_flux()


@pytest.mark.foundation
def test_reconstruct_flux_slab_returns_flux_solution() -> None:
    r"""End-to-end smoke: slab Atkinson reconstruction returns a
    :class:`FluxSolution` with the right shape and tags."""
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    ms = MomentSpace(
        geometry=_ua10sl_geometry(), materials={0: mix}, fn_order=8,
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
    assert np.all(flux.scalar_flux > 0)


@pytest.mark.foundation
def test_reconstruct_flux_sphere_returns_flux_solution() -> None:
    """End-to-end smoke: sphere KLL reconstruction returns a FluxSolution."""
    mix = _make_1g_mixture(1.0, 0.7, 0.6)
    ms = MomentSpace(
        geometry=_ua10sp_geometry(), materials={0: mix}, fn_order=8,
    )
    flux = ms.reconstruct_flux(n_panels=32)
    assert isinstance(flux, FluxSolution)
    assert flux.spatial_nodes.shape == (32,)
    assert flux.scalar_flux.shape == (32,)
    assert flux.eigenvalue_kind == "k_eff"
    assert flux.spatial_units == "mfp"
    assert flux.metadata["geometry"] == "sphere"
    assert flux.metadata["flux_reconstruction"] == "kll_fredholm"
    assert np.all(flux.scalar_flux > 0)
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
    ms = MomentSpace(
        geometry=_ua10sl_geometry(), materials={0: mix}
    )
    with pytest.raises(ValueError, match="multiplying medium"):
        ms.solve_critical()
