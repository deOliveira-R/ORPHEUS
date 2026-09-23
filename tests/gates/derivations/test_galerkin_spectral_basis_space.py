r"""Foundation tests for the Galerkin spectral basis space (Phase D).

These are software-invariant tests on the :class:`BasisSpace` facade
— bit-equality preservation against the function-level API,
StructuredGeometry input acceptance, and result-type contract
preservation.

Phase D
-------

BasisSpace now consumes a :class:`StructuredGeometry` directly via
its frozen ``__init__``. The ``d`` parameter (Galerkin spectral
spatial dimension in mfp) MUST be passed explicitly to
``solve_critical`` — the pre-Phase-D fallback that read from
``GeometrySpec.critical_dimension_mfp`` was retired.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.solution_types import CriticalSolution
from orpheus.derivations.common.xs_library import make_mixture
from orpheus.derivations.continuous.galerkin_spectral import BasisSpace
from orpheus.derivations.continuous.galerkin_spectral.slab import (
    solve_galerkin_spectral_slab,
)
from orpheus.derivations.continuous.galerkin_spectral.sphere import (
    solve_galerkin_spectral_sphere,
)
from orpheus.geometry.mesh import BC
from orpheus.geometry.structured_geometry import (
    Region,
    StructuredGeometry,
)


# ════════════════════════════════════════════════════════════════════════
# Helpers
# ════════════════════════════════════════════════════════════════════════


def _make_1g_isotropic_mixture(
    sigma_t: float, sigma_s: float, nu_sigma_f: float
):
    r"""Build a 1G isotropic-scattering ``Mixture``."""
    return make_mixture(
        sig_t=np.array([sigma_t]),
        sig_c=np.array([0.0]),
        sig_f=np.array([nu_sigma_f / 2.0]),
        nu=np.array([2.0]),
        chi=np.array([1.0]),
        sig_s=np.array([[sigma_s]]),
    )


def _make_1g_p1_mixture(
    sigma_t: float, sigma_s: float, nu_sigma_f: float, mu_bar: float
):
    r"""Build a 1G mixture with linearly anisotropic (P_1) scattering."""
    return make_mixture(
        sig_t=np.array([sigma_t]),
        sig_c=np.array([0.0]),
        sig_f=np.array([nu_sigma_f / 2.0]),
        nu=np.array([2.0]),
        chi=np.array([1.0]),
        sig_s=np.array([[sigma_s]]),
        sig_s1=np.array([[sigma_s * mu_bar]]),
    )


def _slab_geometry(d_mfp: float = 2.0) -> StructuredGeometry:
    """Bare-critical slab geometry. The full slab width = d_mfp (in mfp).

    The geometry's role for BasisSpace is structural (tag + materials
    lookup) — the ``d`` parameter is passed explicitly to
    :meth:`solve_critical`. We use the same value for the
    StructuredGeometry's full extent so the geometry is dimensionally
    coherent, but BasisSpace does not consult it.
    """
    return StructuredGeometry(
        geometry="SLB",
        regions=(Region(mat_id=0, outer_thickness_cm=float(d_mfp)),),
        bcs=(BC.vacuum, BC.vacuum),
    )


def _sphere_geometry(d_mfp: float = 2.0) -> StructuredGeometry:
    """Bare-critical sphere geometry."""
    R = d_mfp / 2.0
    return StructuredGeometry(
        geometry="SPH",
        regions=(Region(mat_id=0, outer_thickness_cm=float(R)),),
        bcs=(BC.vacuum,),
    )


def _hex(x: float) -> str:
    """Exact bit-pattern of a float for IEEE-754 equality tests."""
    return float(x).hex()


# ════════════════════════════════════════════════════════════════════════
# Construction and validation
# ════════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_basis_space_constructs_with_structured_geometry() -> None:
    r"""``BasisSpace`` accepts ``StructuredGeometry`` + ``dict[int, Mixture]``."""
    mix = _make_1g_isotropic_mixture(1.0, 0.7, 0.6)
    geom = _slab_geometry()
    bs = BasisSpace(geometry=geom, materials={0: mix})
    assert isinstance(bs, BasisSpace)
    assert bs.geometry is geom
    assert bs.materials[0] is mix
    assert bs.basis_order == 9  # default
    assert bs.n_quad == 128  # default


@pytest.mark.foundation
def test_basis_space_rejects_cylinder_geometry() -> None:
    r"""Cylinder is out of pillar (Westfall-Metcalf 1972)."""
    mix = _make_1g_isotropic_mixture(1.0, 0.7, 0.6)
    cyl = StructuredGeometry(
        geometry="CYL",
        regions=(Region(mat_id=0, outer_thickness_cm=1.72500292),),
        bcs=(BC.vacuum,),
    )
    with pytest.raises(ValueError, match="singular_eigenfunction"):
        BasisSpace(geometry=cyl, materials={0: mix})


@pytest.mark.foundation
def test_basis_space_rejects_invalid_basis_order() -> None:
    mix = _make_1g_isotropic_mixture(1.0, 0.7, 0.6)
    geom = _slab_geometry()
    with pytest.raises(ValueError, match="basis_order"):
        BasisSpace(geometry=geom, materials={0: mix}, basis_order=0)


@pytest.mark.foundation
def test_basis_space_rejects_invalid_n_quad() -> None:
    mix = _make_1g_isotropic_mixture(1.0, 0.7, 0.6)
    geom = _slab_geometry()
    with pytest.raises(ValueError, match="n_quad"):
        BasisSpace(geometry=geom, materials={0: mix}, n_quad=0)


@pytest.mark.foundation
def test_basis_space_rejects_missing_mat_id() -> None:
    mix = _make_1g_isotropic_mixture(1.0, 0.7, 0.6)
    geom = StructuredGeometry(
        geometry="SLB",
        regions=(Region(mat_id=99, outer_thickness_cm=1.0),),
        bcs=(BC.vacuum, BC.vacuum),
    )
    with pytest.raises(ValueError, match="mat_id=99"):
        BasisSpace(geometry=geom, materials={0: mix})


@pytest.mark.foundation
def test_basis_space_rejects_multi_group_mixture() -> None:
    r"""The shipping Galerkin spectral path is 1G only."""
    mix_2g = make_mixture(
        sig_t=np.array([1.0, 1.0]),
        sig_c=np.array([0.1, 0.1]),
        sig_f=np.array([0.2, 0.2]),
        nu=np.array([2.0, 2.0]),
        chi=np.array([1.0, 0.0]),
        sig_s=np.array([[0.5, 0.0], [0.0, 0.5]]),
    )
    geom = _slab_geometry()
    with pytest.raises(ValueError, match="1G"):
        BasisSpace(geometry=geom, materials={0: mix_2g})


@pytest.mark.foundation
def test_basis_space_is_frozen() -> None:
    r"""The dataclass is frozen — instances are immutable."""
    mix = _make_1g_isotropic_mixture(1.0, 0.7, 0.6)
    bs = BasisSpace(geometry=_slab_geometry(), materials={0: mix})
    with pytest.raises((AttributeError, Exception)):
        bs.basis_order = 5  # type: ignore[misc]


# ════════════════════════════════════════════════════════════════════════
# Derived primary parameters
# ════════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_basis_space_c_property_isotropic() -> None:
    r"""``BasisSpace.c`` returns
    :math:`(\Sigma_s + \nu\Sigma_f) / \Sigma_t` for 1G mixtures."""
    mix = _make_1g_isotropic_mixture(1.0, 0.7, 0.6)
    bs = BasisSpace(geometry=_slab_geometry(), materials={0: mix})
    expected = (0.7 + 0.6) / 1.0
    assert bs.c == pytest.approx(expected, rel=1e-15)


@pytest.mark.foundation
def test_basis_space_mu_bar_property_isotropic_returns_zero() -> None:
    r"""For mixtures without a P_1 moment, ``mu_bar`` returns 0."""
    mix = _make_1g_isotropic_mixture(1.0, 0.7, 0.6)
    bs = BasisSpace(geometry=_slab_geometry(), materials={0: mix})
    assert bs.mu_bar == 0.0


@pytest.mark.foundation
def test_basis_space_mu_bar_property_p1_anisotropic() -> None:
    r"""For mixtures with P_1, ``mu_bar = SigS_p1 / SigS_p0``."""
    mix = _make_1g_p1_mixture(1.0, 0.7, 0.6, mu_bar=0.30)
    bs = BasisSpace(geometry=_slab_geometry(), materials={0: mix})
    assert bs.mu_bar == pytest.approx(0.30, rel=1e-15)


# ════════════════════════════════════════════════════════════════════════
# Bit-equality with function-level API
# ════════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_basis_space_solve_critical_slab_isotropic_bit_equal() -> None:
    r"""``BasisSpace.solve_critical(d=...)`` produces IDENTICAL c_critical
    to a direct ``solve_galerkin_spectral_slab`` call."""
    sigma_t, sigma_s, nu_sigma_f = 1.0, 0.7, 0.6
    d = 2.0
    mix = _make_1g_isotropic_mixture(sigma_t, sigma_s, nu_sigma_f)
    bs = BasisSpace(
        geometry=_slab_geometry(d_mfp=d),
        materials={0: mix},
        basis_order=6,
        n_quad=64,
    )

    res_direct = solve_galerkin_spectral_slab(
        c=bs.c, d=d, mu_bar=0.0,
        n_modes=6, n_quad=64,
    )
    sol = bs.solve_critical(d=d)

    assert _hex(sol.eigenvalue) == _hex(res_direct.c_critical)
    assert np.array_equal(
        sol.metadata["eigenvalue_spectrum"],
        res_direct.eigenvalue_spectrum,
    )


@pytest.mark.foundation
def test_basis_space_solve_critical_slab_anisotropic_bit_equal() -> None:
    r"""Slab P_1-anisotropic case — bit-equal c_critical at
    non-trivial :math:`\bar\mu = 0.20`."""
    sigma_t, sigma_s, nu_sigma_f = 1.0, 0.7, 0.6
    mu_bar = 0.20
    d = 2.0
    mix = _make_1g_p1_mixture(sigma_t, sigma_s, nu_sigma_f, mu_bar=mu_bar)
    bs = BasisSpace(
        geometry=_slab_geometry(d_mfp=d),
        materials={0: mix},
        basis_order=6,
        n_quad=64,
    )

    res_direct = solve_galerkin_spectral_slab(
        c=bs.c, d=d, mu_bar=bs.mu_bar,
        n_modes=6, n_quad=64,
    )
    sol = bs.solve_critical(d=d)
    assert _hex(sol.eigenvalue) == _hex(res_direct.c_critical)


@pytest.mark.foundation
def test_basis_space_solve_critical_sphere_isotropic_bit_equal() -> None:
    r"""Sphere isotropic — bit-equal c_critical."""
    sigma_t, sigma_s, nu_sigma_f = 1.0, 0.7, 0.6
    d = 4.0
    mix = _make_1g_isotropic_mixture(sigma_t, sigma_s, nu_sigma_f)
    bs = BasisSpace(
        geometry=_sphere_geometry(d_mfp=d),
        materials={0: mix},
        basis_order=6,
        n_quad=64,
    )

    res_direct = solve_galerkin_spectral_sphere(
        c=bs.c, d=d, mu_bar=0.0,
        n_modes=6, n_quad=64,
    )
    sol = bs.solve_critical(d=d)
    assert _hex(sol.eigenvalue) == _hex(res_direct.c_critical)


@pytest.mark.foundation
def test_basis_space_solve_critical_sphere_anisotropic_bit_equal() -> None:
    r"""Sphere P_1-anisotropic — bit-equal at :math:`\bar\mu = 0.20`."""
    sigma_t, sigma_s, nu_sigma_f = 1.0, 0.7, 0.6
    mu_bar = 0.20
    d = 4.0
    mix = _make_1g_p1_mixture(sigma_t, sigma_s, nu_sigma_f, mu_bar=mu_bar)
    bs = BasisSpace(
        geometry=_sphere_geometry(d_mfp=d),
        materials={0: mix},
        basis_order=6,
        n_quad=64,
    )

    res_direct = solve_galerkin_spectral_sphere(
        c=bs.c, d=d, mu_bar=bs.mu_bar,
        n_modes=6, n_quad=64,
    )
    sol = bs.solve_critical(d=d)
    assert _hex(sol.eigenvalue) == _hex(res_direct.c_critical)


@pytest.mark.foundation
def test_basis_space_solve_critical_requires_explicit_d() -> None:
    r"""Phase D: ``d`` MUST be passed; no GeometrySpec-based fallback."""
    sigma_t, sigma_s, nu_sigma_f = 1.0, 0.7, 0.6
    mix = _make_1g_isotropic_mixture(sigma_t, sigma_s, nu_sigma_f)
    bs = BasisSpace(
        geometry=_slab_geometry(d_mfp=1.0),
        materials={0: mix},
        basis_order=6,
        n_quad=64,
    )
    with pytest.raises(ValueError, match="explicit"):
        bs.solve_critical()


# ════════════════════════════════════════════════════════════════════════
# CriticalSolution result type contract
# ════════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_basis_space_solve_critical_returns_critical_solution_with_correct_tags() -> None:
    r"""The return type and field tags conform to the cross-method
    :class:`CriticalSolution` contract."""
    mix = _make_1g_isotropic_mixture(1.0, 0.7, 0.6)
    bs = BasisSpace(
        geometry=_slab_geometry(d_mfp=2.0),
        materials={0: mix},
        basis_order=6,
        n_quad=64,
    )
    sol = bs.solve_critical(d=2.0)

    assert isinstance(sol, CriticalSolution)
    assert sol.eigenvalue_kind == "c_critical"
    assert sol.parameter_kind == "domain_extent_cm"
    assert sol.eigenvalue > 0.0
    spectrum = sol.metadata["eigenvalue_spectrum"]
    assert spectrum.shape == (2 * bs.basis_order,)
    real_mask = np.abs(spectrum.imag) < 1e-9
    real_positive = spectrum[real_mask].real
    real_positive = real_positive[real_positive > 0.0]
    assert sol.eigenvalue == pytest.approx(np.min(real_positive))


@pytest.mark.foundation
def test_basis_space_solve_critical_metadata_carries_diagnostics() -> None:
    r"""Method-specific diagnostics are reachable via ``metadata``."""
    mix = _make_1g_p1_mixture(1.0, 0.7, 0.6, mu_bar=0.20)
    bs = BasisSpace(
        geometry=_slab_geometry(d_mfp=2.0),
        materials={0: mix},
        basis_order=6,
        n_quad=64,
    )
    sol = bs.solve_critical(d=2.0)
    md = sol.metadata
    assert md["geometry"] == "slab"
    assert md["d_mfp"] == 2.0
    assert md["mu_bar"] == pytest.approx(0.20, rel=1e-14)
    assert md["basis_order"] == 6
    assert md["n_quad"] == 64
    assert md["c"] == pytest.approx(1.3, rel=1e-15)
    assert md["method"] == "solve_galerkin_spectral_slab"
    assert md["raw_result"].c_critical == sol.eigenvalue


# ════════════════════════════════════════════════════════════════════════
# solve_full_spectrum — distinguishing feature
# ════════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_basis_space_solve_full_spectrum_returns_2N_eigenvalues() -> None:
    r"""The full-spectrum entry returns 2N eigenvalues."""
    mix = _make_1g_isotropic_mixture(1.0, 0.7, 0.6)
    bs = BasisSpace(
        geometry=_slab_geometry(d_mfp=2.0),
        materials={0: mix},
        basis_order=5,
        n_quad=64,
    )
    spec = bs.solve_full_spectrum(d=2.0)
    assert spec["eigenvalue_spectrum"].shape == (2 * 5,)
    assert spec["eigenvectors"].shape == (5, 2 * 5)
    assert spec["c_critical"] > 0.0
