"""L0 parity tests for the unified flat-source CP builder.

Phase B.2a gate: the unified :func:`~cp_geometry.build_cp_matrix`
must reproduce the legacy ``_slab_cp_matrix``, ``_cylinder_cp_matrix``,
and ``_sphere_cp_matrix`` outputs **bit-identically** — all three
kernel paths use the same underlying primitives on both sides:

- Slab: scipy ``expn(3, x)`` via :func:`_kernels.e3_vec`.
- Cylinder: legacy :class:`_kernels.BickleyTables` (retires in
  Phase B.4 / :issue:`94` when the precision upgrade to
  ``ki_n_mp(3, ·, 30)`` lands).
- Sphere: ``np.exp(-tau)``.

These tests pin the B.2a commit before B.2b swaps the three facades
to delegate to the unified builder. If this file passes in Commit 1
and the facades' L1 eigenvalue tests pass in Commit 2, the refactor
is safe."""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations import cp_cylinder as _legacy_cyl
from orpheus.derivations import cp_slab as _legacy_slab
from orpheus.derivations import cp_sphere as _legacy_sph
from orpheus.derivations.cp_geometry import (
    CYLINDER_1D,
    SLAB,
    SPHERE_1D,
    build_cp_matrix,
)
from orpheus.derivations._xs_library import LAYOUTS, get_xs


# ═══════════════════════════════════════════════════════════════════════
# Helpers — materialise the same test inputs each facade uses
# ═══════════════════════════════════════════════════════════════════════

def _slab_inputs(ng_key: str, n_regions: int):
    layout = LAYOUTS[n_regions]
    t_arr = np.array(_legacy_slab._THICKNESSES[n_regions], dtype=float)
    sig_t_all = np.vstack([get_xs(r, ng_key)["sig_t"] for r in layout])
    return sig_t_all, t_arr


def _cyl_inputs(ng_key: str, n_regions: int):
    layout = LAYOUTS[n_regions]
    radii = np.array(_legacy_cyl._RADII[n_regions], dtype=float)
    r_inner = np.zeros(n_regions)
    r_inner[1:] = radii[:-1]
    volumes = np.pi * (radii ** 2 - r_inner ** 2)
    sig_t_all = np.vstack([get_xs(r, ng_key)["sig_t"] for r in layout])
    return sig_t_all, radii, volumes, float(radii[-1])


def _sph_inputs(ng_key: str, n_regions: int):
    layout = LAYOUTS[n_regions]
    radii = np.array(_legacy_sph._RADII[n_regions], dtype=float)
    r_inner = np.zeros(n_regions)
    r_inner[1:] = radii[:-1]
    volumes = (4.0 / 3.0) * np.pi * (radii ** 3 - r_inner ** 3)
    sig_t_all = np.vstack([get_xs(r, ng_key)["sig_t"] for r in layout])
    return sig_t_all, radii, volumes, float(radii[-1])


# ═══════════════════════════════════════════════════════════════════════
# Slab: bit-identity (same scipy kernel on both sides)
# ═══════════════════════════════════════════════════════════════════════

@pytest.mark.l0
@pytest.mark.parametrize("ng_key", ["1g", "2g", "4g"])
@pytest.mark.parametrize("n_regions", [1, 2, 4])
def test_unified_slab_matches_legacy(ng_key, n_regions):
    """Unified slab P_inf matches legacy ``_slab_cp_matrix`` to machine
    precision (same scipy ``expn(3, x)`` on both sides)."""
    sig_t_all, t_arr = _slab_inputs(ng_key, n_regions)
    R_cell = float(t_arr.sum())

    legacy = _legacy_slab._slab_cp_matrix(sig_t_all, t_arr)
    unified = build_cp_matrix(
        SLAB, sig_t_all, t_arr, t_arr, R_cell, n_quad_y=64,
    )

    np.testing.assert_allclose(
        unified, legacy, rtol=1e-12, atol=1e-14,
        err_msg=f"slab {ng_key} {n_regions}rg unified vs legacy drift",
    )


# ═══════════════════════════════════════════════════════════════════════
# Cylinder: bit-identity (same BickleyTables on both sides in Phase B.2)
# ═══════════════════════════════════════════════════════════════════════

@pytest.mark.l0
@pytest.mark.parametrize("ng_key", ["1g", "2g", "4g"])
@pytest.mark.parametrize("n_regions", [1, 2, 4])
def test_unified_cylinder_matches_legacy(ng_key, n_regions):
    """Unified cylinder P_inf matches legacy ``_cylinder_cp_matrix``
    bit-identically — same ``BickleyTables.Ki3_vec`` on both sides in
    Phase B.2. Phase B.4 will introduce the Ki_3 precision upgrade."""
    sig_t_all, radii, volumes, R_cell = _cyl_inputs(ng_key, n_regions)

    legacy = _legacy_cyl._cylinder_cp_matrix(
        sig_t_all, radii, volumes, R_cell, n_quad_y=64,
    )
    unified = build_cp_matrix(
        CYLINDER_1D, sig_t_all, radii, volumes, R_cell, n_quad_y=64,
    )

    np.testing.assert_allclose(
        unified, legacy, rtol=1e-12, atol=1e-14,
        err_msg=f"cylinder {ng_key} {n_regions}rg unified vs legacy drift",
    )


# ═══════════════════════════════════════════════════════════════════════
# Sphere: bit-identity (same np.exp on both sides)
# ═══════════════════════════════════════════════════════════════════════

@pytest.mark.l0
@pytest.mark.parametrize("ng_key", ["1g", "2g", "4g"])
@pytest.mark.parametrize("n_regions", [1, 2, 4])
def test_unified_sphere_matches_legacy(ng_key, n_regions):
    """Unified sphere P_inf matches legacy ``_sphere_cp_matrix`` to
    machine precision (same ``np.exp(-tau)`` on both sides)."""
    sig_t_all, radii, volumes, R_cell = _sph_inputs(ng_key, n_regions)

    legacy = _legacy_sph._sphere_cp_matrix(
        sig_t_all, radii, volumes, R_cell, n_quad_y=64,
    )
    unified = build_cp_matrix(
        SPHERE_1D, sig_t_all, radii, volumes, R_cell, n_quad_y=64,
    )

    np.testing.assert_allclose(
        unified, legacy, rtol=1e-12, atol=1e-14,
        err_msg=f"sphere {ng_key} {n_regions}rg unified vs legacy drift",
    )


# ═══════════════════════════════════════════════════════════════════════
# Δ² operator geometry-invariance sanity check
# ═══════════════════════════════════════════════════════════════════════

@pytest.mark.l0
@pytest.mark.verifies("cp-second-difference-operator")
class TestSecondDifferenceOperator:
    """The Δ² operator is geometry-invariant: it depends on the kernel
    callable and three τ arguments only. Same inputs → same output
    regardless of which geometry 'owns' the kernel."""

    def test_identity_holds_for_all_three_kernels(self):
        from orpheus.derivations.cp_geometry import _second_difference

        gap = np.array([0.5, 1.0, 2.0])
        tau_i = np.array([0.3, 0.7, 1.1])
        tau_j = np.array([0.2, 0.4, 0.9])

        for geom in (SLAB, CYLINDER_1D, SPHERE_1D):
            direct = (geom.kernel_F3(gap)
                      - geom.kernel_F3(gap + tau_i)
                      - geom.kernel_F3(gap + tau_j)
                      + geom.kernel_F3(gap + tau_i + tau_j))
            via_op = _second_difference(
                geom.kernel_F3, gap, tau_i, tau_j,
            )
            np.testing.assert_allclose(
                via_op, direct, rtol=1e-14, atol=1e-14,
                err_msg=f"Δ² operator mismatch for {geom.kind}",
            )

    def test_degenerate_zero_tau_gives_zero(self):
        from orpheus.derivations.cp_geometry import _second_difference

        gap = np.array([0.5, 1.0, 2.0])
        zero = np.zeros_like(gap)
        for geom in (SLAB, CYLINDER_1D, SPHERE_1D):
            val = _second_difference(geom.kernel_F3, gap, zero, zero)
            np.testing.assert_allclose(val, 0.0, atol=1e-14)


# ═══════════════════════════════════════════════════════════════════════
# Phase B.4 marker: Ki_3 precision upgrade is deferred
# ═══════════════════════════════════════════════════════════════════════

@pytest.mark.l0
def test_cylinder_kernel_is_bickley_in_phase_b2():
    """Phase B.2 pins the cylinder kernel to the legacy BickleyTables
    for bit-identity with the pre-refactor builder. Phase B.4 will
    swap to ``ki_n_mp(3, ·, 30)`` and retire :class:`BickleyTables`
    (:issue:`94`). When that commit lands, this test flips."""
    from orpheus.derivations.cp_geometry import _ki3_legacy
    from orpheus.derivations._kernels import bickley_tables

    tau = np.array([0.0, 0.5, 1.0, 2.5, 5.0, 10.0])
    np.testing.assert_array_equal(
        _ki3_legacy(tau), bickley_tables().Ki3_vec(tau),
    )
