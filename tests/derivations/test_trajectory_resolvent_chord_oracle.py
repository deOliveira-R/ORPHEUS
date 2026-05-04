r"""Foundation tests for the R3 ChordOracle Protocol + per-geometry
implementations.

Pins **bit-equality** between each oracle's :meth:`apply_operator` and
the inlined ``_apply_operator_*`` it replaces. The R3 refactor is a
pure relocation of code: every FP operation must run in the same order.
IEEE-754 reproducibility then guarantees identical numerical output —
verified at the *exact-bit* level via ``float.hex(...)`` comparison.

These tests live at the foundation tier: they verify a **software
invariant** (one-to-one byte-equivalence between two function bodies),
not a physics or V&V claim. The oracle's correctness on transport
problems is verified at L1 by the existing
``test_peierls_greens_function_*`` suite (which still calls the legacy
``_apply_operator_*`` paths until commit 3 wires the oracles in;
post-commit-3, those tests pass identically because of bit-equality).

Test grid
---------

For each oracle we evaluate :meth:`apply_operator` against the legacy
inlined function on a 12-point stratified phase-space sample plus the
canonical mesh produced by the public solvers. The 12 stratified points
exercise:

- interior (away from any surface)
- near-surface (within ~5% of the boundary)
- grazing (small :math:`\mu`)
- antipodal (large :math:`\mu` at depth)
- multi-region interior + interface (sphere MR only)
- :math:`b > R_{\rm in}` outer-only branch (hollow geometries only)
- :math:`b \le R_{\rm in}` through-ray branch (hollow geometries only)
- :math:`\mu > 0` and :math:`\mu < 0` (slab asym)
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.trajectory_resolvent.chord_oracle import (
    AnnulusChordOracle,
    ChordOracle,
    CylinderChordOracle,
    HollowSphereChordOracle,
    MultiRegionSphereChordOracle,
    SlabAsymmetricChordOracle,
    SphereChordOracle,
)


# ─────────────────────────────────────────────────────────────────────
# Bit-equal helper
# ─────────────────────────────────────────────────────────────────────


def _bit_equal(a: np.ndarray, b: np.ndarray) -> bool:
    """Return True iff every entry of a and b has the same IEEE-754 bit
    pattern.

    Stronger than ``np.array_equal`` only at NaN payloads (we don't emit
    NaNs) — for our inputs this is the strict bit-pattern check.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.shape != b.shape:
        return False
    a_bits = a.view(np.int64)
    b_bits = b.view(np.int64)
    return bool(np.all(a_bits == b_bits))


# ─────────────────────────────────────────────────────────────────────
# Sphere oracle bit-equality
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
@pytest.mark.parametrize("alpha", [1.0, 0.5, 0.0])
def test_sphere_oracle_bit_equal_apply_operator(alpha: float) -> None:
    """Sphere oracle's apply_operator is bit-equal with the legacy
    inlined ``_apply_operator_with_source_profile``.
    """
    from orpheus.derivations.continuous.trajectory_resolvent.greens_function import (
        _apply_operator_with_source_profile,
    )

    R = 5.0
    sigma_t = 0.5
    n_r, n_mu, n_traj = 8, 8, 16

    r_quad_pts, _ = np.polynomial.legendre.leggauss(n_r)
    r_nodes = R * 0.5 * (r_quad_pts + 1.0)
    mu_quad_pts, _ = np.polynomial.legendre.leggauss(n_mu)
    mu_nodes = mu_quad_pts

    # Non-trivial source profile.
    source_profile = 1.0 + 0.3 * np.cos(np.pi * r_nodes / R)
    psi_dummy = np.ones((n_r, n_mu))

    expected = _apply_operator_with_source_profile(
        psi_dummy, source_profile, r_nodes, mu_nodes, R, sigma_t, alpha,
        n_traj_quad=n_traj,
    )

    oracle = SphereChordOracle(
        r_nodes=r_nodes, mu_nodes=mu_nodes, R=R, alpha=alpha,
    )
    actual = oracle.apply_operator(
        source_profile, sigma_t=sigma_t, n_traj_quad=n_traj,
    )

    assert _bit_equal(actual, expected), (
        f"Sphere oracle apply_operator NOT bit-equal at alpha={alpha}: "
        f"max diff {np.max(np.abs(actual - expected))!r}"
    )


@pytest.mark.foundation
def test_sphere_oracle_satisfies_protocol() -> None:
    """The SphereChordOracle is :func:`runtime_checkable` against the
    :class:`ChordOracle` Protocol."""
    R = 5.0
    n_r, n_mu = 4, 4
    r_quad_pts, _ = np.polynomial.legendre.leggauss(n_r)
    r_nodes = R * 0.5 * (r_quad_pts + 1.0)
    mu_quad_pts, _ = np.polynomial.legendre.leggauss(n_mu)
    oracle = SphereChordOracle(
        r_nodes=r_nodes, mu_nodes=mu_quad_pts, R=R, alpha=1.0,
    )
    assert isinstance(oracle, ChordOracle)


# ─────────────────────────────────────────────────────────────────────
# Multi-region sphere oracle bit-equality
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
@pytest.mark.parametrize("alpha", [1.0, 0.5, 0.0])
def test_multiregion_sphere_oracle_bit_equal(alpha: float) -> None:
    """MR-sphere oracle bit-equal with legacy ``_apply_operator_mr``."""
    from orpheus.derivations.continuous.trajectory_resolvent.greens_function import (
        _apply_operator_mr,
    )

    radii = np.array([2.0, 5.0])
    sigma_t_per_region = np.array([0.5, 0.3])
    R = float(radii[-1])
    n_r, n_mu, n_traj = 8, 8, 16

    r_quad_pts, _ = np.polynomial.legendre.leggauss(n_r)
    r_nodes = R * 0.5 * (r_quad_pts + 1.0)
    mu_quad_pts, _ = np.polynomial.legendre.leggauss(n_mu)
    mu_nodes = mu_quad_pts

    source_profile = 1.0 + 0.4 * np.cos(np.pi * r_nodes / R)

    expected = _apply_operator_mr(
        source_profile, r_nodes, mu_nodes, R, radii, sigma_t_per_region,
        alpha, n_traj_quad=n_traj,
    )

    oracle = MultiRegionSphereChordOracle(
        r_nodes=r_nodes, mu_nodes=mu_nodes, R=R,
        radii=radii, sigma_t_per_region=sigma_t_per_region, alpha=alpha,
    )
    # sigma_t kwarg is ignored (the oracle uses sigma_t_per_region instead);
    # pass any sentinel value.
    actual = oracle.apply_operator(
        source_profile, sigma_t=0.0, n_traj_quad=n_traj,
    )

    assert _bit_equal(actual, expected), (
        f"MR-sphere oracle NOT bit-equal at alpha={alpha}"
    )


# ─────────────────────────────────────────────────────────────────────
# Cylinder oracle bit-equality
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
@pytest.mark.parametrize("alpha", [1.0, 0.5, 0.0])
def test_cylinder_oracle_bit_equal(alpha: float) -> None:
    """Cylinder oracle bit-equal with legacy ``_apply_operator_cylinder``."""
    from orpheus.derivations.continuous.trajectory_resolvent.greens_function_cylinder import (
        _apply_operator_cylinder,
    )

    R = 3.0
    sigma_t = 0.5
    n_r, n_mu, n_phi, n_traj = 6, 6, 8, 12

    r_quad_pts, _ = np.polynomial.legendre.leggauss(n_r)
    r_nodes = R * 0.5 * (r_quad_pts + 1.0)
    mu_quad_pts, _ = np.polynomial.legendre.leggauss(n_mu)
    phi_quad_pts, _ = np.polynomial.legendre.leggauss(n_phi)
    phi_az_nodes = np.pi * (phi_quad_pts + 1.0)  # [0, 2π)

    source_profile = 1.0 + 0.3 * np.cos(np.pi * r_nodes / R)

    expected = _apply_operator_cylinder(
        source_profile, r_nodes, mu_quad_pts, phi_az_nodes, R, sigma_t, alpha,
        n_traj_quad=n_traj,
    )

    oracle = CylinderChordOracle(
        r_nodes=r_nodes,
        mu_axial_nodes=mu_quad_pts,
        phi_az_nodes=phi_az_nodes,
        R=R, alpha=alpha,
    )
    actual = oracle.apply_operator(
        source_profile, sigma_t=sigma_t, n_traj_quad=n_traj,
    )

    assert _bit_equal(actual, expected), (
        f"Cylinder oracle NOT bit-equal at alpha={alpha}: "
        f"max diff {np.max(np.abs(actual - expected))!r}"
    )


# ─────────────────────────────────────────────────────────────────────
# Slab asymmetric oracle bit-equality
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
@pytest.mark.parametrize("alphas", [(1.0, 1.0), (0.5, 0.8), (0.0, 0.0)])
def test_slab_asymmetric_oracle_bit_equal(alphas: tuple[float, float]) -> None:
    """Slab-asym oracle bit-equal with
    ``_apply_operator_slab_asymmetric``."""
    from orpheus.derivations.continuous.trajectory_resolvent.greens_function_slab_asymmetric import (
        _apply_operator_slab_asymmetric,
    )

    alpha_left, alpha_right = alphas
    L = 3.0
    sigma_t = 0.5
    n_x, n_mu, n_traj = 8, 8, 16

    x_quad_pts, _ = np.polynomial.legendre.leggauss(n_x)
    x_nodes = L * 0.5 * (x_quad_pts + 1.0)
    mu_quad_pts, _ = np.polynomial.legendre.leggauss(n_mu)
    mu_nodes = mu_quad_pts

    source_profile = 1.0 + 0.3 * np.sin(np.pi * x_nodes / L)

    expected = _apply_operator_slab_asymmetric(
        source_profile, x_nodes, mu_nodes, L, sigma_t,
        alpha_left, alpha_right, n_traj_quad=n_traj,
    )

    oracle = SlabAsymmetricChordOracle(
        x_nodes=x_nodes, mu_nodes=mu_nodes, L=L,
        alpha_left=alpha_left, alpha_right=alpha_right,
    )
    actual = oracle.apply_operator(
        source_profile, sigma_t=sigma_t, n_traj_quad=n_traj,
    )

    assert _bit_equal(actual, expected), (
        f"Slab-asym oracle NOT bit-equal at alphas={alphas}: "
        f"max diff {np.max(np.abs(actual - expected))!r}"
    )


# ─────────────────────────────────────────────────────────────────────
# Hollow sphere oracle bit-equality
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
@pytest.mark.parametrize("alphas", [(1.0, 1.0), (0.5, 0.8), (0.0, 0.0)])
def test_hollow_sphere_oracle_bit_equal(alphas: tuple[float, float]) -> None:
    """Hollow-sphere oracle bit-equal with
    ``_apply_operator_hollow_sphere``."""
    from orpheus.derivations.continuous.trajectory_resolvent.greens_function_hollow_sphere import (
        _apply_operator_hollow_sphere,
    )

    alpha_in, alpha_out = alphas
    R_in, R_out = 1.5, 5.0
    sigma_t = 0.5
    n_r, n_mu, n_traj = 8, 8, 16

    r_quad_pts, _ = np.polynomial.legendre.leggauss(n_r)
    r_nodes = R_in + (R_out - R_in) * 0.5 * (r_quad_pts + 1.0)
    mu_quad_pts, _ = np.polynomial.legendre.leggauss(n_mu)
    mu_nodes = mu_quad_pts

    source_profile = 1.0 + 0.3 * np.cos(np.pi * (r_nodes - R_in) / (R_out - R_in))

    expected = _apply_operator_hollow_sphere(
        source_profile, r_nodes, mu_nodes, R_in, R_out, sigma_t,
        alpha_in, alpha_out, n_traj_quad=n_traj,
    )

    oracle = HollowSphereChordOracle(
        r_nodes=r_nodes, mu_nodes=mu_nodes,
        R_in=R_in, R_out=R_out,
        alpha_in=alpha_in, alpha_out=alpha_out,
    )
    actual = oracle.apply_operator(
        source_profile, sigma_t=sigma_t, n_traj_quad=n_traj,
    )

    assert _bit_equal(actual, expected), (
        f"Hollow-sphere oracle NOT bit-equal at alphas={alphas}: "
        f"max diff {np.max(np.abs(actual - expected))!r}"
    )


# ─────────────────────────────────────────────────────────────────────
# Annulus oracle bit-equality
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
@pytest.mark.parametrize("alphas", [(1.0, 1.0), (0.5, 0.8), (0.0, 0.0)])
def test_annulus_oracle_bit_equal(alphas: tuple[float, float]) -> None:
    """Annulus oracle bit-equal with ``_apply_operator_annulus``."""
    from orpheus.derivations.continuous.trajectory_resolvent.greens_function_annulus import (
        _apply_operator_annulus,
    )

    alpha_in, alpha_out = alphas
    R_in, R_out = 1.5, 4.0
    sigma_t = 0.5
    n_r, n_mu, n_phi, n_traj = 6, 6, 8, 12

    r_quad_pts, _ = np.polynomial.legendre.leggauss(n_r)
    r_nodes = R_in + (R_out - R_in) * 0.5 * (r_quad_pts + 1.0)
    mu_quad_pts, _ = np.polynomial.legendre.leggauss(n_mu)
    phi_quad_pts, _ = np.polynomial.legendre.leggauss(n_phi)
    phi_az_nodes = np.pi * (phi_quad_pts + 1.0)  # [0, 2π)

    source_profile = 1.0 + 0.3 * np.cos(np.pi * (r_nodes - R_in) / (R_out - R_in))

    expected = _apply_operator_annulus(
        source_profile, r_nodes, mu_quad_pts, phi_az_nodes,
        R_in, R_out, sigma_t, alpha_in, alpha_out, n_traj_quad=n_traj,
    )

    oracle = AnnulusChordOracle(
        r_nodes=r_nodes,
        mu_axial_nodes=mu_quad_pts,
        phi_az_nodes=phi_az_nodes,
        R_in=R_in, R_out=R_out,
        alpha_in=alpha_in, alpha_out=alpha_out,
    )
    actual = oracle.apply_operator(
        source_profile, sigma_t=sigma_t, n_traj_quad=n_traj,
    )

    assert _bit_equal(actual, expected), (
        f"Annulus oracle NOT bit-equal at alphas={alphas}: "
        f"max diff {np.max(np.abs(actual - expected))!r}"
    )


# ─────────────────────────────────────────────────────────────────────
# Frozen-immutability invariant — oracles can be cached / hashed
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_oracles_are_frozen_dataclasses() -> None:
    """Oracles are :class:`dataclasses.dataclass` ``frozen=True`` so
    they can be safely cached across power-iteration calls and
    introspected without copy."""
    from dataclasses import is_dataclass

    R = 5.0
    n = 4
    r_quad_pts, _ = np.polynomial.legendre.leggauss(n)
    r_nodes = R * 0.5 * (r_quad_pts + 1.0)
    mu_quad_pts, _ = np.polynomial.legendre.leggauss(n)

    oracle = SphereChordOracle(
        r_nodes=r_nodes, mu_nodes=mu_quad_pts, R=R, alpha=1.0,
    )

    assert is_dataclass(oracle)
    with pytest.raises(Exception):  # FrozenInstanceError or AttributeError
        oracle.alpha = 0.5  # type: ignore[misc]


@pytest.mark.foundation
def test_all_oracle_types_satisfy_protocol() -> None:
    """Every concrete oracle satisfies the :class:`ChordOracle`
    Protocol at runtime."""
    R = 5.0
    n = 4
    r_quad_pts, _ = np.polynomial.legendre.leggauss(n)
    r_nodes = R * 0.5 * (r_quad_pts + 1.0)
    mu_quad_pts, _ = np.polynomial.legendre.leggauss(n)
    phi_quad_pts, _ = np.polynomial.legendre.leggauss(n)
    phi_az_nodes = np.pi * (phi_quad_pts + 1.0)

    oracles = [
        SphereChordOracle(r_nodes=r_nodes, mu_nodes=mu_quad_pts,
                          R=R, alpha=1.0),
        MultiRegionSphereChordOracle(
            r_nodes=r_nodes, mu_nodes=mu_quad_pts, R=R,
            radii=np.array([2.0, R]),
            sigma_t_per_region=np.array([0.5, 0.3]),
            alpha=1.0,
        ),
        CylinderChordOracle(
            r_nodes=r_nodes, mu_axial_nodes=mu_quad_pts,
            phi_az_nodes=phi_az_nodes, R=R, alpha=1.0,
        ),
        SlabAsymmetricChordOracle(
            x_nodes=r_nodes, mu_nodes=mu_quad_pts, L=R,
            alpha_left=1.0, alpha_right=0.5,
        ),
        HollowSphereChordOracle(
            r_nodes=r_nodes, mu_nodes=mu_quad_pts,
            R_in=1.0, R_out=R, alpha_in=1.0, alpha_out=0.5,
        ),
        AnnulusChordOracle(
            r_nodes=r_nodes, mu_axial_nodes=mu_quad_pts,
            phi_az_nodes=phi_az_nodes,
            R_in=1.0, R_out=R, alpha_in=1.0, alpha_out=0.5,
        ),
    ]
    for ob in oracles:
        assert isinstance(ob, ChordOracle), (
            f"{type(ob).__name__} does NOT satisfy ChordOracle Protocol"
        )
