r"""Foundation tests for the Galerkin spectral basis space (math-heart class).

These are software-invariant tests on the :class:`BasisSpace` facade
— bit-equality preservation against the function-level API,
production-protocol input acceptance, and result-type contract
preservation. They are NOT verification claims about the Galerkin
spectral method's accuracy: those are owned by
``test_carlvik_galerkin_*.py`` (which exercise the function-level API
and are unchanged by this dispatch).

The tests pin three invariants:

1. **Bit-equality with the function-level API** — the class produces
   IDENTICAL float results to direct calls of
   ``solve_galerkin_spectral_slab`` / ``solve_galerkin_spectral_sphere``
   when given the same numerical inputs. Verified via ``float.hex(...)``
   exact-bit comparison on the reported ``c_critical`` and via
   ``np.array_equal`` on the eigenvalue spectrum.

2. **Production-protocol input acceptance** — the class consumes the
   same ``(materials: dict[int, Mixture], GeometrySpec)`` pair as
   production CP/SN/MOC solvers, so a single problem definition serves
   both. Verified by constructing from
   :func:`orpheus.derivations.common.xs_library.make_mixture` +
   :class:`GeometrySpec`.

3. **TransportSolver Protocol conformance** — :class:`BasisSpace`
   exposes ``materials``, ``geometry_spec``, ``method_name``, and
   ``solve_critical()`` as required by
   :class:`~orpheus.derivations.common.solver_protocol.TransportSolver`.

Tests are foundation-tagged because they verify *structural* invariants
on the class facade, not L0/L1/L2 claims about the solver's accuracy —
the latter live in ``test_carlvik_galerkin_slab.py``,
``test_carlvik_galerkin_sphere.py``, ``test_carlvik_galerkin_xverif_fn.py``.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.geometry_spec import GeometrySpec
from orpheus.derivations.common.solution_types import CriticalSolution
from orpheus.derivations.common.xs_library import make_mixture
from orpheus.derivations.continuous.galerkin_spectral import BasisSpace

# The TransportSolver Protocol may land on main via a parallel feature
# branch (see ``.claude/plans/transport_solver_protocol_spec.md``). The
# Protocol-conformance tests below skip gracefully when the module is
# absent; the rest of the foundation suite (bit-equality, type
# contracts) is unaffected.
try:
    from orpheus.derivations.common.solver_protocol import (  # type: ignore[import-not-found]
        KNOWN_TRANSPORT_SOLVERS,
        TransportSolver,
    )

    _HAS_PROTOCOL = True
except ImportError:  # pragma: no cover (parallel-branch dependency)
    _HAS_PROTOCOL = False
from orpheus.derivations.continuous.galerkin_spectral.slab import (
    solve_galerkin_spectral_slab,
)
from orpheus.derivations.continuous.galerkin_spectral.sphere import (
    solve_galerkin_spectral_sphere,
)
from orpheus.geometry.mesh import BC


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
    r"""Build a 1G mixture with linearly anisotropic (P_1) scattering.

    The mixture's :math:`P_1` scattering moment is set to
    :math:`\Sigma_{s,1} = \bar\mu \cdot \Sigma_{s,0}`, the
    scattering-only convention.
    """
    return make_mixture(
        sig_t=np.array([sigma_t]),
        sig_c=np.array([0.0]),
        sig_f=np.array([nu_sigma_f / 2.0]),
        nu=np.array([2.0]),
        chi=np.array([1.0]),
        sig_s=np.array([[sigma_s]]),
        sig_s1=np.array([[sigma_s * mu_bar]]),
    )


def _slab_geometry(d_mfp: float = 2.0) -> GeometrySpec:
    """Bare-critical slab geometry. ``d = 2 * critical_dimension_mfp``."""
    a = d_mfp / 2.0
    return GeometrySpec(
        geometry="slab",
        critical_dimension_mfp=a,
        critical_dimension_cm=a,
        n_groups=1,
        bc_left=BC.vacuum,
        bc_right=BC.vacuum,
    )


def _sphere_geometry(d_mfp: float = 2.0) -> GeometrySpec:
    """Bare-critical sphere geometry. ``d = 2 * critical_dimension_mfp``."""
    R = d_mfp / 2.0
    return GeometrySpec(
        geometry="sphere",
        critical_dimension_mfp=R,
        critical_dimension_cm=R,
        n_groups=1,
        bc_left=BC.reflective,
        bc_right=BC.vacuum,
    )


def _hex(x: float) -> str:
    """Exact bit-pattern of a float for IEEE-754 equality tests."""
    return float(x).hex()


# ════════════════════════════════════════════════════════════════════════
# Construction and validation
# ════════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_basis_space_constructs_with_production_protocol_inputs() -> None:
    r"""``BasisSpace.from_problem`` accepts ``dict[int, Mixture] +
    GeometrySpec`` (the production-protocol input shape)."""
    mix = _make_1g_isotropic_mixture(1.0, 0.7, 0.6)
    geom = _slab_geometry()
    bs = BasisSpace.from_problem(materials={0: mix}, geometry=geom)
    assert isinstance(bs, BasisSpace)
    assert bs.geometry_spec is geom
    assert bs.materials[0] is mix
    assert bs.basis_order == 9  # default
    assert bs.n_quad == 128  # default


@pytest.mark.foundation
def test_basis_space_rejects_cylinder_geometry() -> None:
    r"""Cylinder is out of pillar (Westfall-Metcalf 1972); the class
    rejects it at construction time with an explicit error message
    naming the alternative pillar."""
    mix = _make_1g_isotropic_mixture(1.0, 0.7, 0.6)
    cyl = GeometrySpec(
        geometry="cylinder",
        critical_dimension_mfp=1.72500292,
        critical_dimension_cm=1.72500292,
        n_groups=1,
        bc_left=BC.reflective,
        bc_right=BC.vacuum,
    )
    with pytest.raises(ValueError, match="singular_eigenfunction"):
        BasisSpace.from_problem(materials={0: mix}, geometry=cyl)


@pytest.mark.foundation
def test_basis_space_rejects_infinite_geometry() -> None:
    r"""Infinite medium has no spatial Galerkin expansion; reject."""
    mix = _make_1g_isotropic_mixture(1.0, 0.7, 0.6)
    inf_geom = GeometrySpec(
        geometry="infinite",
        critical_dimension_mfp=None,
        critical_dimension_cm=None,
        n_groups=1,
    )
    with pytest.raises(ValueError, match="kinf_homogeneous"):
        BasisSpace.from_problem(materials={0: mix}, geometry=inf_geom)


@pytest.mark.foundation
def test_basis_space_rejects_invalid_basis_order() -> None:
    mix = _make_1g_isotropic_mixture(1.0, 0.7, 0.6)
    geom = _slab_geometry()
    with pytest.raises(ValueError, match="basis_order"):
        BasisSpace.from_problem(materials={0: mix}, geometry=geom, basis_order=0)


@pytest.mark.foundation
def test_basis_space_rejects_invalid_n_quad() -> None:
    mix = _make_1g_isotropic_mixture(1.0, 0.7, 0.6)
    geom = _slab_geometry()
    with pytest.raises(ValueError, match="n_quad"):
        BasisSpace.from_problem(materials={0: mix}, geometry=geom, n_quad=0)


@pytest.mark.foundation
def test_basis_space_rejects_missing_mat_id() -> None:
    mix = _make_1g_isotropic_mixture(1.0, 0.7, 0.6)
    geom = GeometrySpec(
        geometry="slab",
        critical_dimension_mfp=1.0,
        critical_dimension_cm=1.0,
        n_groups=1,
        mat_id=99,  # not present in materials
    )
    with pytest.raises(ValueError, match="mat_id=99"):
        BasisSpace.from_problem(materials={0: mix}, geometry=geom)


@pytest.mark.foundation
def test_basis_space_rejects_multi_group_mixture() -> None:
    r"""The shipping Galerkin spectral path is 1G only; multi-group is
    deferred follow-up."""
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
        BasisSpace.from_problem(materials={0: mix_2g}, geometry=geom)


@pytest.mark.foundation
def test_basis_space_is_frozen() -> None:
    r"""The dataclass is frozen — instances are immutable."""
    mix = _make_1g_isotropic_mixture(1.0, 0.7, 0.6)
    bs = BasisSpace.from_problem(materials={0: mix}, geometry=_slab_geometry())
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
    bs = BasisSpace.from_problem(materials={0: mix}, geometry=_slab_geometry())
    expected = (0.7 + 0.6) / 1.0
    assert bs.c == pytest.approx(expected, rel=1e-15)


@pytest.mark.foundation
def test_basis_space_mu_bar_property_isotropic_returns_zero() -> None:
    r"""For mixtures without a P_1 moment, ``mu_bar`` returns 0."""
    mix = _make_1g_isotropic_mixture(1.0, 0.7, 0.6)
    bs = BasisSpace.from_problem(materials={0: mix}, geometry=_slab_geometry())
    assert bs.mu_bar == 0.0


@pytest.mark.foundation
def test_basis_space_mu_bar_property_p1_anisotropic() -> None:
    r"""For mixtures with P_1, ``mu_bar = SigS_p1 / SigS_p0``."""
    mix = _make_1g_p1_mixture(1.0, 0.7, 0.6, mu_bar=0.30)
    bs = BasisSpace.from_problem(materials={0: mix}, geometry=_slab_geometry())
    assert bs.mu_bar == pytest.approx(0.30, rel=1e-15)


# ════════════════════════════════════════════════════════════════════════
# Bit-equality with function-level API
# ════════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_basis_space_solve_critical_slab_isotropic_bit_equal() -> None:
    r"""``BasisSpace.solve_critical()`` produces IDENTICAL c_critical to
    a direct ``solve_galerkin_spectral_slab(c=..., d=..., mu_bar=0.0,
    n_modes=basis_order, n_quad=n_quad)`` call (verified via
    ``float.hex(...)`` exact-bit comparison)."""
    sigma_t, sigma_s, nu_sigma_f = 1.0, 0.7, 0.6
    d = 2.0
    mix = _make_1g_isotropic_mixture(sigma_t, sigma_s, nu_sigma_f)
    geom = _slab_geometry(d_mfp=d)
    bs = BasisSpace.from_problem(
        materials={0: mix}, geometry=geom, basis_order=6, n_quad=64
    )

    # Direct function call.
    res_direct = solve_galerkin_spectral_slab(
        c=bs.c, d=d, mu_bar=0.0,
        n_modes=6, n_quad=64,
    )

    # Class facade.
    sol = bs.solve_critical()

    # Bit-equal c_critical.
    assert _hex(sol.eigenvalue) == _hex(res_direct.c_critical)
    # Bit-equal full eigenvalue spectrum.
    assert np.array_equal(
        sol.metadata["eigenvalue_spectrum"],
        res_direct.eigenvalue_spectrum,
    )


@pytest.mark.foundation
def test_basis_space_solve_critical_slab_anisotropic_bit_equal() -> None:
    r"""Slab P_1-anisotropic case — bit-equal c_critical from the class
    facade vs direct function call at non-trivial :math:`\bar\mu = 0.20`.

    Note: we compare against the SAME ``mu_bar`` derived from the
    mixture (``bs.mu_bar``) — passing the constructor's input
    ``mu_bar`` directly to the function-level API would round-trip
    ``Sigma_s_p1 / Sigma_s_p0`` and lose 1 ULP. The class facade
    is bit-exact w.r.t. ``bs.mu_bar``, which is the right contract.
    """
    sigma_t, sigma_s, nu_sigma_f = 1.0, 0.7, 0.6
    mu_bar = 0.20
    d = 2.0
    mix = _make_1g_p1_mixture(sigma_t, sigma_s, nu_sigma_f, mu_bar=mu_bar)
    geom = _slab_geometry(d_mfp=d)
    bs = BasisSpace.from_problem(
        materials={0: mix}, geometry=geom, basis_order=6, n_quad=64
    )

    res_direct = solve_galerkin_spectral_slab(
        c=bs.c, d=d, mu_bar=bs.mu_bar,
        n_modes=6, n_quad=64,
    )
    sol = bs.solve_critical()
    assert _hex(sol.eigenvalue) == _hex(res_direct.c_critical)


@pytest.mark.foundation
def test_basis_space_solve_critical_sphere_isotropic_bit_equal() -> None:
    r"""Sphere isotropic — bit-equal c_critical from the class facade vs
    direct function call."""
    sigma_t, sigma_s, nu_sigma_f = 1.0, 0.7, 0.6
    d = 4.0
    mix = _make_1g_isotropic_mixture(sigma_t, sigma_s, nu_sigma_f)
    geom = _sphere_geometry(d_mfp=d)
    bs = BasisSpace.from_problem(
        materials={0: mix}, geometry=geom, basis_order=6, n_quad=64
    )

    res_direct = solve_galerkin_spectral_sphere(
        c=bs.c, d=d, mu_bar=0.0,
        n_modes=6, n_quad=64,
    )
    sol = bs.solve_critical()
    assert _hex(sol.eigenvalue) == _hex(res_direct.c_critical)


@pytest.mark.foundation
def test_basis_space_solve_critical_sphere_anisotropic_bit_equal() -> None:
    r"""Sphere P_1-anisotropic case — bit-equal c_critical at
    non-trivial :math:`\bar\mu = 0.20`. Pass ``bs.mu_bar`` (not the
    constructor input) for the same reason as the slab case."""
    sigma_t, sigma_s, nu_sigma_f = 1.0, 0.7, 0.6
    mu_bar = 0.20
    d = 4.0
    mix = _make_1g_p1_mixture(sigma_t, sigma_s, nu_sigma_f, mu_bar=mu_bar)
    geom = _sphere_geometry(d_mfp=d)
    bs = BasisSpace.from_problem(
        materials={0: mix}, geometry=geom, basis_order=6, n_quad=64
    )

    res_direct = solve_galerkin_spectral_sphere(
        c=bs.c, d=d, mu_bar=bs.mu_bar,
        n_modes=6, n_quad=64,
    )
    sol = bs.solve_critical()
    assert _hex(sol.eigenvalue) == _hex(res_direct.c_critical)


@pytest.mark.foundation
def test_basis_space_solve_critical_explicit_d_overrides_geometry_spec() -> None:
    r"""Passing ``d=`` to ``solve_critical`` overrides the value derived
    from ``GeometrySpec.critical_dimension_mfp``."""
    sigma_t, sigma_s, nu_sigma_f = 1.0, 0.7, 0.6
    mix = _make_1g_isotropic_mixture(sigma_t, sigma_s, nu_sigma_f)
    # Set GeometrySpec to one dimension; pass a different d to solve_critical.
    geom = _slab_geometry(d_mfp=1.0)  # critical_dimension_mfp = 0.5
    bs = BasisSpace.from_problem(
        materials={0: mix}, geometry=geom, basis_order=6, n_quad=64
    )

    sol_default = bs.solve_critical()  # uses d = 2 * 0.5 = 1.0
    sol_override = bs.solve_critical(d=2.0)  # explicit d
    # The two configurations are different — eigenvalues should differ.
    assert sol_default.eigenvalue != sol_override.eigenvalue
    # The override should match a direct call.
    res_direct = solve_galerkin_spectral_slab(
        c=bs.c, d=2.0, mu_bar=0.0, n_modes=6, n_quad=64
    )
    assert _hex(sol_override.eigenvalue) == _hex(res_direct.c_critical)


# ════════════════════════════════════════════════════════════════════════
# CriticalSolution result type contract
# ════════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_basis_space_solve_critical_returns_critical_solution_with_correct_tags() -> None:
    r"""The return type and field tags conform to the cross-method
    :class:`CriticalSolution` contract."""
    mix = _make_1g_isotropic_mixture(1.0, 0.7, 0.6)
    bs = BasisSpace.from_problem(
        materials={0: mix},
        geometry=_slab_geometry(d_mfp=2.0),
        basis_order=6,
        n_quad=64,
    )
    sol = bs.solve_critical()

    assert isinstance(sol, CriticalSolution)
    # Galerkin spectral reports c_critical, NOT k_eff. The tag
    # MUST disambiguate so cross-method comparators don't compare to
    # k_eff = 1 by accident.
    assert sol.eigenvalue_kind == "c_critical"
    assert sol.parameter_kind == "domain_extent_cm"
    # Eigenvalue is positive real (fundamental).
    assert sol.eigenvalue > 0.0
    # Spectrum is 2 * basis_order long (full Eq. 4 block matrix).
    spectrum = sol.metadata["eigenvalue_spectrum"]
    assert spectrum.shape == (2 * bs.basis_order,)
    # The reported c_critical equals the smallest positive real eigenvalue.
    real_mask = np.abs(spectrum.imag) < 1e-9
    real_positive = spectrum[real_mask].real
    real_positive = real_positive[real_positive > 0.0]
    assert sol.eigenvalue == pytest.approx(np.min(real_positive))


@pytest.mark.foundation
def test_basis_space_solve_critical_metadata_carries_diagnostics() -> None:
    r"""Method-specific diagnostics are reachable via ``metadata``."""
    mix = _make_1g_p1_mixture(1.0, 0.7, 0.6, mu_bar=0.20)
    bs = BasisSpace.from_problem(
        materials={0: mix},
        geometry=_slab_geometry(d_mfp=2.0),
        basis_order=6,
        n_quad=64,
    )
    sol = bs.solve_critical()
    md = sol.metadata
    assert md["geometry"] == "slab"
    assert md["d_mfp"] == 2.0
    # mu_bar round-trips Sigma_s_p1 / Sigma_s_p0 → loses ≤ 1 ULP relative
    # to the user-input mu_bar; rel-1e-15 tolerance suffices.
    assert md["mu_bar"] == pytest.approx(0.20, rel=1e-14)
    assert md["basis_order"] == 6
    assert md["n_quad"] == 64
    assert md["c"] == pytest.approx(1.3, rel=1e-15)
    assert md["method"] == "solve_galerkin_spectral_slab"
    # raw_result preserves bit-equal access to the legacy result type.
    assert md["raw_result"].c_critical == sol.eigenvalue


# ════════════════════════════════════════════════════════════════════════
# solve_full_spectrum — distinguishing feature
# ════════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_basis_space_solve_full_spectrum_returns_2N_eigenvalues() -> None:
    r"""The full-spectrum entry returns 2N eigenvalues — Galerkin
    spectral's distinguishing feature relative to other math-heart classes."""
    mix = _make_1g_isotropic_mixture(1.0, 0.7, 0.6)
    bs = BasisSpace.from_problem(
        materials={0: mix},
        geometry=_slab_geometry(d_mfp=2.0),
        basis_order=5,
        n_quad=64,
    )
    spec = bs.solve_full_spectrum()
    assert spec["eigenvalue_spectrum"].shape == (2 * 5,)
    assert spec["eigenvectors"].shape == (5, 2 * 5)
    # The reported c_critical is the smallest positive real.
    assert spec["c_critical"] > 0.0


# ════════════════════════════════════════════════════════════════════════
# TransportSolver Protocol conformance
# ════════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_basis_space_exposes_protocol_shaped_attributes() -> None:
    r""":class:`BasisSpace` exposes the ``materials``, ``geometry_spec``,
    ``method_name`` attributes required by the
    :class:`TransportSolver` Protocol — independent of whether the
    Protocol module is currently importable.

    The class is designed Protocol-shaped: when the parallel
    ``feat/transport-solver-protocol`` branch lands on main,
    :class:`BasisSpace` automatically conforms to the structural
    contract via :func:`runtime_checkable` Protocol.
    """
    mix = _make_1g_isotropic_mixture(1.0, 0.7, 0.6)
    bs = BasisSpace.from_problem(
        materials={0: mix}, geometry=_slab_geometry()
    )
    # The three required Protocol attributes.
    assert isinstance(bs.materials, dict)
    assert isinstance(bs.geometry_spec, GeometrySpec)
    assert bs.method_name == "galerkin_spectral"
    # The Protocol method.
    sol = bs.solve_critical()
    assert isinstance(sol, CriticalSolution)


@pytest.mark.foundation
@pytest.mark.skipif(
    not _HAS_PROTOCOL,
    reason="TransportSolver Protocol module not yet on this branch's HEAD; "
    "lands via parallel feat/transport-solver-protocol branch",
)
def test_basis_space_satisfies_transport_solver_protocol() -> None:
    r""":class:`BasisSpace` is structurally compatible with the
    :class:`TransportSolver` Protocol via :func:`runtime_checkable`."""
    mix = _make_1g_isotropic_mixture(1.0, 0.7, 0.6)
    bs = BasisSpace.from_problem(
        materials={0: mix}, geometry=_slab_geometry()
    )
    # Protocol membership via runtime_checkable.
    assert isinstance(bs, TransportSolver)
    # Each Protocol field is exposed.
    assert isinstance(bs.materials, dict)
    assert isinstance(bs.geometry_spec, GeometrySpec)
    assert bs.method_name == "galerkin_spectral"
    sol = bs.solve_critical()
    assert isinstance(sol, CriticalSolution)


@pytest.mark.foundation
@pytest.mark.skipif(
    not _HAS_PROTOCOL,
    reason="TransportSolver Protocol module not yet on this branch's HEAD; "
    "the registry tuple is defined alongside the Protocol on the "
    "parallel feat/transport-solver-protocol branch",
)
def test_galerkin_spectral_in_known_transport_solvers_registry() -> None:
    r"""The pillar tag ``"galerkin_spectral"`` is registered in
    :data:`KNOWN_TRANSPORT_SOLVERS`. The registry is the foundation
    gate that catches silent removal of ``method_name`` from a
    math-heart class."""
    assert "galerkin_spectral" in KNOWN_TRANSPORT_SOLVERS


@pytest.mark.foundation
@pytest.mark.skipif(
    not _HAS_PROTOCOL,
    reason="TransportSolver Protocol module not yet on this branch's HEAD",
)
def test_basis_space_method_name_matches_known_pillar_tag() -> None:
    r"""The ``method_name`` field on every instance equals the registered
    pillar tag — no per-instance drift."""
    mix = _make_1g_isotropic_mixture(1.0, 0.7, 0.6)
    bs = BasisSpace.from_problem(
        materials={0: mix}, geometry=_slab_geometry()
    )
    assert bs.method_name in KNOWN_TRANSPORT_SOLVERS
    assert bs.method_name == "galerkin_spectral"
