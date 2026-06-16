"""Hash-equality + property tests for the reduced streaming operator.

The load-bearing contract: arrays produced by the geometry-layer
factories (:func:`slab_streaming`, :func:`spherical_streaming`,
:func:`cylindrical_streaming`) MUST be **bit-identical** to those
produced by the legacy in-line setup methods on
:class:`~orpheus.sn.geometry.SNMesh`.

These tests are tagged ``@pytest.mark.foundation`` — they verify
software invariants (hash-equality of array bit patterns, property
correctness) rather than equation-level math claims.  The
equation-level math is verified transitively via the 11 SN regression
snapshots that gate the entire reshape campaign.

Wave B Round 1 — Issue 6 (.claude/plans/sn_reshape.md, GH #155).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import (
    BC,
    CoordSystem,
    Mesh1D,
    ReducedStreamingOperator,
    StreamingTerms,
    cylindrical_streaming,
    slab_streaming,
    spherical_streaming,
)
from orpheus.sn.geometry import SNMesh
from orpheus.numerics.quadrature import Quadrature
from tests.sn._test_helpers import placeholder_materials


# ═══════════════════════════════════════════════════════════════════════
# Mesh fixtures
# ═══════════════════════════════════════════════════════════════════════

def _slab_mesh() -> Mesh1D:
    """Cartesian slab — 5-cell uniform mesh on [0, 1]."""
    return Mesh1D(
        edges=np.linspace(0.0, 1.0, 6),
        mat_ids=np.zeros(5, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )


def _spherical_mesh() -> Mesh1D:
    """Spherical 5-cell mesh on [0, 1]."""
    return Mesh1D(
        edges=np.linspace(0.0, 1.0, 6),
        mat_ids=np.zeros(5, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )


def _cylindrical_mesh() -> Mesh1D:
    """Cylindrical 5-cell mesh on [0, 1]."""
    return Mesh1D(
        edges=np.linspace(0.0, 1.0, 6),
        mat_ids=np.zeros(5, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )


# ═══════════════════════════════════════════════════════════════════════
# Hash-equality tests (the primary gate)
# ═══════════════════════════════════════════════════════════════════════

class TestHashEqualitySpherical:
    """Spherical: every precomputed array bit-identical to SNMesh."""

    @pytest.fixture
    def pair(self):
        mesh = _spherical_mesh()
        quad = Quadrature.gauss_legendre(8)
        sn_mesh = SNMesh(mesh, quad, placeholder_materials())
        reduced = spherical_streaming(mesh, quad)
        return sn_mesh, reduced

    @pytest.mark.foundation
    def test_face_areas_bit_identical(self, pair):
        sn_mesh, reduced = pair
        assert reduced.face_areas is not None
        assert np.array_equal(reduced.face_areas, sn_mesh.face_areas)

    @pytest.mark.foundation
    def test_delta_A_bit_identical(self, pair):
        sn_mesh, reduced = pair
        assert reduced.delta_A is not None
        assert np.array_equal(reduced.delta_A, sn_mesh.delta_A)

    @pytest.mark.foundation
    def test_alpha_half_bit_identical(self, pair):
        sn_mesh, reduced = pair
        assert reduced.alpha_half is not None
        assert np.array_equal(reduced.alpha_half, sn_mesh.reduced.alpha_half)

    @pytest.mark.foundation
    def test_redist_dAw_bit_identical(self, pair):
        sn_mesh, reduced = pair
        assert reduced.redist_dAw is not None
        assert np.array_equal(reduced.redist_dAw, sn_mesh.reduced.redist_dAw)

    @pytest.mark.foundation
    def test_tau_mm_bit_identical(self, pair):
        sn_mesh, reduced = pair
        assert reduced.tau_mm is not None
        assert np.array_equal(reduced.tau_mm, sn_mesh.reduced.tau_mm)

    @pytest.mark.foundation
    @pytest.mark.parametrize("N", [4, 8, 16, 32])
    def test_alpha_dome_bit_identical_across_N(self, N):
        """Hash equality holds for every quadrature order."""
        mesh = _spherical_mesh()
        quad = Quadrature.gauss_legendre(N)
        sn_mesh = SNMesh(mesh, quad, placeholder_materials())
        reduced = spherical_streaming(mesh, quad)
        assert np.array_equal(reduced.alpha_half, sn_mesh.reduced.alpha_half)
        assert np.array_equal(reduced.tau_mm, sn_mesh.reduced.tau_mm)
        assert np.array_equal(reduced.redist_dAw, sn_mesh.reduced.redist_dAw)


class TestHashEqualityCylindrical:
    """Cylindrical: per-level arrays bit-identical to SNMesh."""

    @pytest.fixture
    def pair(self):
        mesh = _cylindrical_mesh()
        quad = Quadrature.product(n_mu=2, n_phi=4)
        sn_mesh = SNMesh(mesh, quad, placeholder_materials())
        reduced = cylindrical_streaming(mesh, quad)
        return sn_mesh, reduced

    @pytest.mark.foundation
    def test_face_areas_bit_identical(self, pair):
        sn_mesh, reduced = pair
        assert np.array_equal(reduced.face_areas, sn_mesh.face_areas)

    @pytest.mark.foundation
    def test_delta_A_bit_identical(self, pair):
        sn_mesh, reduced = pair
        assert np.array_equal(reduced.delta_A, sn_mesh.delta_A)

    @pytest.mark.foundation
    def test_alpha_per_level_bit_identical(self, pair):
        sn_mesh, reduced = pair
        assert reduced.alpha_per_level is not None
        assert len(reduced.alpha_per_level) == len(sn_mesh.reduced.alpha_per_level)
        for lvl, (rdc, snm) in enumerate(
            zip(reduced.alpha_per_level, sn_mesh.reduced.alpha_per_level)
        ):
            assert np.array_equal(rdc, snm), f"level {lvl} mismatch"

    @pytest.mark.foundation
    def test_redist_dAw_per_level_bit_identical(self, pair):
        sn_mesh, reduced = pair
        assert reduced.redist_dAw_per_level is not None
        for lvl, (rdc, snm) in enumerate(
            zip(reduced.redist_dAw_per_level, sn_mesh.reduced.redist_dAw_per_level)
        ):
            assert np.array_equal(rdc, snm), f"level {lvl} mismatch"

    @pytest.mark.foundation
    def test_tau_mm_per_level_bit_identical(self, pair):
        sn_mesh, reduced = pair
        assert reduced.tau_mm_per_level is not None
        for lvl, (rdc, snm) in enumerate(
            zip(reduced.tau_mm_per_level, sn_mesh.reduced.tau_mm_per_level)
        ):
            assert np.array_equal(rdc, snm), f"level {lvl} mismatch"

    @pytest.mark.foundation
    @pytest.mark.parametrize("n_mu,n_phi", [(2, 4), (4, 4), (4, 8)])
    def test_full_per_level_hash_equality(self, n_mu, n_phi):
        """Per-level hash equality across multiple quadrature shapes."""
        mesh = _cylindrical_mesh()
        quad = Quadrature.product(n_mu=n_mu, n_phi=n_phi)
        sn_mesh = SNMesh(mesh, quad, placeholder_materials())
        reduced = cylindrical_streaming(mesh, quad)
        for rdc, snm in zip(
            reduced.alpha_per_level, sn_mesh.reduced.alpha_per_level
        ):
            assert np.array_equal(rdc, snm)
        for rdc, snm in zip(
            reduced.redist_dAw_per_level, sn_mesh.reduced.redist_dAw_per_level
        ):
            assert np.array_equal(rdc, snm)
        for rdc, snm in zip(
            reduced.tau_mm_per_level, sn_mesh.reduced.tau_mm_per_level
        ):
            assert np.array_equal(rdc, snm)


# ═══════════════════════════════════════════════════════════════════════
# Property tests
# ═══════════════════════════════════════════════════════════════════════

class TestProperties:
    """ReducedStreamingOperator advertises the right metadata."""

    @pytest.mark.foundation
    def test_slab_no_upstream_no_axis(self):
        mesh = _slab_mesh()
        quad = Quadrature.gauss_legendre(8)
        op = slab_streaming(mesh, quad)
        assert op.requires_upstream_angular_state is False
        assert op.angular_marching_axis is None
        assert op.coord is CoordSystem.CARTESIAN

    @pytest.mark.foundation
    def test_slab_curvature_arrays_are_none(self):
        mesh = _slab_mesh()
        quad = Quadrature.gauss_legendre(8)
        op = slab_streaming(mesh, quad)
        assert op.face_areas is None
        assert op.delta_A is None
        assert op.alpha_half is None
        assert op.redist_dAw is None
        assert op.tau_mm is None
        assert op.alpha_per_level is None
        assert op.redist_dAw_per_level is None
        assert op.tau_mm_per_level is None

    @pytest.mark.foundation
    def test_sphere_requires_upstream(self):
        mesh = _spherical_mesh()
        quad = Quadrature.gauss_legendre(8)
        op = spherical_streaming(mesh, quad)
        assert op.requires_upstream_angular_state is True
        assert op.angular_marching_axis == "mu"
        assert op.coord is CoordSystem.SPHERICAL

    @pytest.mark.foundation
    def test_cylinder_requires_upstream(self):
        mesh = _cylindrical_mesh()
        quad = Quadrature.product(n_mu=2, n_phi=4)
        op = cylindrical_streaming(mesh, quad)
        assert op.requires_upstream_angular_state is True
        assert op.angular_marching_axis == "mu"
        assert op.coord is CoordSystem.CYLINDRICAL


# ═══════════════════════════════════════════════════════════════════════
# Per-direction extraction tests
# ═══════════════════════════════════════════════════════════════════════

class TestStreamingTermsExtraction:
    """streaming_terms() returns the right per-cell-per-direction subset."""

    @pytest.mark.foundation
    def test_slab_streaming_terms_neutral_curvature(self):
        """Slab carries neutral curvature values (Issue #196 Step 2.5).

        Pre-Step-2.5 the curvature fields were ``None`` and the
        cell-update strategies branched on ``alpha_in is None``.
        Step 2.5 retires that branch by populating neutral values:
        ``face_area_inner = face_area_outer = 1.0``,
        ``delta_A_over_w = 0.0``, ``alpha_in = alpha_out = 0.0``,
        ``tau_mm = 1.0`` (synthetic neutral element of the M-M
        closure — slab has no half-angles).  These values make the
        unified cell-balance helper consume the same packet for
        slab and curvilinear without geometry dispatch.
        """
        mesh = _slab_mesh()
        quad = Quadrature.gauss_legendre(8)
        op = slab_streaming(mesh, quad)
        st = op.streaming_terms(cell_idx=2, direction_idx=3)
        assert isinstance(st, StreamingTerms)
        assert st.chord_length == float(mesh.widths[2])
        assert st.mu == float(quad.mu_x[3])
        # Neutral curvature: slab carries the values that make the
        # unified cell-balance algebra collapse to the slab form.
        assert st.face_area_inner == 1.0
        assert st.face_area_outer == 1.0
        assert st.delta_A_over_w == 0.0
        assert st.alpha_in == 0.0
        assert st.alpha_out == 0.0
        assert st.tau_mm == 1.0

    @pytest.mark.foundation
    def test_sphere_streaming_terms_match_arrays(self):
        mesh = _spherical_mesh()
        quad = Quadrature.gauss_legendre(8)
        op = spherical_streaming(mesh, quad)
        i, n = 1, 5
        st = op.streaming_terms(cell_idx=i, direction_idx=n)
        assert st.chord_length == float(mesh.widths[i])
        # Geometric labels: inner = closer to r=0 (A[i]),
        #                   outer = farther (A[i+1]).
        assert st.face_area_inner == float(op.face_areas[i])
        assert st.face_area_outer == float(op.face_areas[i + 1])
        assert st.delta_A_over_w == float(op.redist_dAw[i, n])
        assert st.alpha_in == float(op.alpha_half[n])
        assert st.alpha_out == float(op.alpha_half[n + 1])
        assert st.tau_mm == float(op.tau_mm[n])
        # Sphere also exposes signed mu (global ordinate index == n).
        assert st.mu == float(quad.mu_x[n])

    @pytest.mark.foundation
    def test_cylinder_streaming_terms_match_per_level(self):
        mesh = _cylindrical_mesh()
        quad = Quadrature.product(n_mu=4, n_phi=4)
        op = cylindrical_streaming(mesh, quad)
        i, level, m = 2, 1, 0
        st = op.streaming_terms(
            cell_idx=i, direction_idx=m, mu_level_idx=level,
        )
        assert st.chord_length == float(mesh.widths[i])
        # Geometric labels: inner / outer relative to r=0,
        # independent of sweep direction.
        assert st.face_area_inner == float(op.face_areas[i])
        assert st.face_area_outer == float(op.face_areas[i + 1])
        assert st.delta_A_over_w == float(
            op.redist_dAw_per_level[level][i, m]
        )
        assert st.alpha_in == float(op.alpha_per_level[level][m])
        assert st.alpha_out == float(op.alpha_per_level[level][m + 1])
        assert st.tau_mm == float(op.tau_mm_per_level[level][m])
        # Cylinder mu is signed eta from the GLOBAL ordinate index
        # (resolved via level_indices) — bug 2 fix anchor.
        global_n = int(quad.level_indices[level][m])
        assert st.mu == float(quad.mu_x[global_n])
        assert st.abs_mu == float(abs(quad.mu_x[global_n]))

    @pytest.mark.foundation
    def test_cylinder_streaming_terms_requires_level(self):
        mesh = _cylindrical_mesh()
        quad = Quadrature.product(n_mu=2, n_phi=4)
        op = cylindrical_streaming(mesh, quad)
        with pytest.raises(ValueError, match="mu_level_idx"):
            op.streaming_terms(cell_idx=0, direction_idx=0)


# ═══════════════════════════════════════════════════════════════════════
# Wave C Round 1: volume + abs_mu fields on StreamingTerms
# ═══════════════════════════════════════════════════════════════════════

class TestStreamingTermsVolumeAndAbsMu:
    """``volume`` and ``abs_mu`` are populated by all three factories.

    These are the two new fields added in Wave C Round 1 (Issue #157)
    so that downstream cell-update strategies receive a self-contained
    per-cell, per-direction packet without reaching back into
    ``SNMesh`` or the ``AngularQuadrature``.
    """

    @pytest.mark.foundation
    @pytest.mark.parametrize("cell_idx", [0, 4])
    @pytest.mark.parametrize("direction_idx", [0, 5])
    def test_slab_volume_matches_mesh(self, cell_idx, direction_idx):
        mesh = _slab_mesh()
        quad = Quadrature.gauss_legendre(8)
        op = slab_streaming(mesh, quad)
        st = op.streaming_terms(cell_idx=cell_idx, direction_idx=direction_idx)
        assert st.volume == float(mesh.volumes[cell_idx])
        assert st.volume is not None
        assert st.volume > 0.0

    @pytest.mark.foundation
    @pytest.mark.parametrize("cell_idx", [0, 4])
    @pytest.mark.parametrize("direction_idx", [0, 5])
    def test_slab_abs_mu_matches_quadrature(self, cell_idx, direction_idx):
        mesh = _slab_mesh()
        quad = Quadrature.gauss_legendre(8)
        op = slab_streaming(mesh, quad)
        st = op.streaming_terms(cell_idx=cell_idx, direction_idx=direction_idx)
        assert st.abs_mu == float(abs(quad.mu_x[direction_idx]))
        assert st.abs_mu is not None
        assert st.abs_mu >= 0.0

    @pytest.mark.foundation
    @pytest.mark.parametrize("N", [4, 8, 16, 32])
    def test_spherical_volume_and_abs_mu_match(self, N):
        mesh = _spherical_mesh()
        quad = Quadrature.gauss_legendre(N)
        op = spherical_streaming(mesh, quad)
        for cell_idx in (0, mesh.N - 1):
            for direction_idx in (0, N - 1):
                st = op.streaming_terms(
                    cell_idx=cell_idx, direction_idx=direction_idx,
                )
                assert st.volume == float(mesh.volumes[cell_idx])
                assert st.volume > 0.0
                assert st.abs_mu == float(abs(quad.mu_x[direction_idx]))

    @pytest.mark.foundation
    @pytest.mark.parametrize("n_mu,n_phi", [(2, 4), (4, 4), (4, 8)])
    def test_cylindrical_volume_and_abs_mu_match(self, n_mu, n_phi):
        mesh = _cylindrical_mesh()
        quad = Quadrature.product(n_mu=n_mu, n_phi=n_phi)
        op = cylindrical_streaming(mesh, quad)
        # Pick a cell-idx pair and use level 0's first ordinate.
        level = 0
        m = 0
        absolute_idx = quad.level_indices[level][m]
        for cell_idx in (0, mesh.N - 1):
            st = op.streaming_terms(
                cell_idx=cell_idx, direction_idx=m, mu_level_idx=level,
            )
            assert st.volume == float(mesh.volumes[cell_idx])
            assert st.volume > 0.0
            assert st.abs_mu == float(abs(quad.mu_x[absolute_idx]))

    @pytest.mark.foundation
    def test_existing_curvilinear_fields_unchanged(self):
        """Anchor: extending StreamingTerms didn't drift the existing fields.

        Re-checks that ``alpha_in``, ``alpha_out``, ``delta_A_over_w``,
        ``tau_mm`` still match the underlying arrays after the
        ``volume`` / ``abs_mu`` extension.  Defense against accidental
        reordering or drift.
        """
        mesh = _spherical_mesh()
        quad = Quadrature.gauss_legendre(8)
        op = spherical_streaming(mesh, quad)
        i, n = 1, 5
        st = op.streaming_terms(cell_idx=i, direction_idx=n)
        assert st.alpha_in == float(op.alpha_half[n])
        assert st.alpha_out == float(op.alpha_half[n + 1])
        assert st.delta_A_over_w == float(op.redist_dAw[i, n])
        assert st.tau_mm == float(op.tau_mm[n])


# ═══════════════════════════════════════════════════════════════════════
# Geometric labels: face_area_inner/outer are direction-independent
# ═══════════════════════════════════════════════════════════════════════

class TestStreamingTermsGeometricLabels:
    """``face_area_inner`` / ``face_area_outer`` are **purely geometric**.

    The two face-area fields encode position relative to :math:`r=0`
    (inner = closer, outer = farther) — they do NOT depend on the
    sweep's marching direction.  The same cell yields the same
    inner / outer values regardless of which ordinate is queried.
    Sweep-direction resolution lives in the SN module
    (:class:`orpheus.sn.spatial.scheme.CellVisit`).
    """

    @pytest.mark.foundation
    def test_sphere_faces_invariant_under_direction(self):
        """Sphere: inner/outer faces are the same for ±μ ordinates."""
        mesh = _spherical_mesh()
        quad = Quadrature.gauss_legendre(8)  # μ ordered ascending
        op = spherical_streaming(mesh, quad)
        i = 2
        # n_neg has μ < 0 (inward); n_pos has μ > 0 (outward).
        n_neg, n_pos = 0, quad.N - 1
        st_neg = op.streaming_terms(cell_idx=i, direction_idx=n_neg)
        st_pos = op.streaming_terms(cell_idx=i, direction_idx=n_pos)
        # Geometric labels — same for both directions.
        assert st_neg.face_area_inner == st_pos.face_area_inner
        assert st_neg.face_area_outer == st_pos.face_area_outer
        # Anchor: inner is A[i], outer is A[i+1] regardless of μ.
        assert st_neg.face_area_inner == float(op.face_areas[i])
        assert st_neg.face_area_outer == float(op.face_areas[i + 1])
        # Signed mu is direction-dependent — the discriminator.
        assert st_neg.mu < 0
        assert st_pos.mu > 0

    @pytest.mark.foundation
    def test_cylinder_faces_invariant_under_direction(self):
        """Cylinder: inner/outer faces are direction-independent."""
        mesh = _cylindrical_mesh()
        quad = Quadrature.product(n_mu=4, n_phi=4)
        op = cylindrical_streaming(mesh, quad)
        i, level = 2, 1
        # Find a within-level pair with η < 0 and η > 0 by scanning.
        level_idx = quad.level_indices[level]
        m_neg = next(
            j for j in range(len(level_idx))
            if quad.mu_x[level_idx[j]] < 0
        )
        m_pos = next(
            j for j in range(len(level_idx))
            if quad.mu_x[level_idx[j]] > 0
        )
        st_neg = op.streaming_terms(
            cell_idx=i, direction_idx=m_neg, mu_level_idx=level,
        )
        st_pos = op.streaming_terms(
            cell_idx=i, direction_idx=m_pos, mu_level_idx=level,
        )
        # Geometric labels — same.
        assert st_neg.face_area_inner == st_pos.face_area_inner
        assert st_neg.face_area_outer == st_pos.face_area_outer
        assert st_neg.face_area_inner == float(op.face_areas[i])
        assert st_neg.face_area_outer == float(op.face_areas[i + 1])
        # Signed eta is direction-dependent.
        assert st_neg.mu < 0
        assert st_pos.mu > 0


# ═══════════════════════════════════════════════════════════════════════
# Cylindrical abs_mu uses the GLOBAL ordinate index (Bug 2 regression)
# ═══════════════════════════════════════════════════════════════════════

class TestCylindricalAbsMuUsesGlobalOrdinate:
    """``abs_mu`` for cylindrical reads :math:`|\\eta|` at the GLOBAL ordinate.

    Pre-fix (Wave D), the cylindrical
    :meth:`ReducedStreamingOperator.streaming_terms` computed
    ``abs_mu = abs(quad.mu_x[direction_idx])`` where ``direction_idx``
    is the **within-level** azimuthal index — pulling the wrong
    global ordinate's :math:`|\\eta|`.  The fix resolves the global
    ordinate via ``level_indices[mu_level_idx][direction_idx]``.

    The Wave D Round 2 sweep had a workaround
    ``_streaming_terms_with_abs_mu`` that overrode the wrong value;
    this test pins the in-place fix so the workaround can be deleted.
    """

    @pytest.mark.foundation
    def test_cylindrical_abs_mu_per_level(self):
        """For every (level, m_local), abs_mu == |mu_x[level_idx[m_local]]|."""
        mesh = _cylindrical_mesh()
        quad = Quadrature.product(n_mu=4, n_phi=4)
        op = cylindrical_streaming(mesh, quad)
        for level, level_idx in enumerate(quad.level_indices):
            for m_local in range(len(level_idx)):
                global_n = int(level_idx[m_local])
                st = op.streaming_terms(
                    cell_idx=0,
                    direction_idx=m_local,
                    mu_level_idx=level,
                )
                # The fix: read from the GLOBAL ordinate index.
                expected = float(abs(quad.mu_x[global_n]))
                assert st.abs_mu == expected, (
                    f"level={level} m_local={m_local} "
                    f"global_n={global_n}: "
                    f"abs_mu={st.abs_mu} != {expected}"
                )

    @pytest.mark.foundation
    def test_cylindrical_signed_mu_per_level(self):
        """Signed ``mu`` (= η) also reads from the GLOBAL ordinate."""
        mesh = _cylindrical_mesh()
        quad = Quadrature.product(n_mu=4, n_phi=4)
        op = cylindrical_streaming(mesh, quad)
        for level, level_idx in enumerate(quad.level_indices):
            for m_local in range(len(level_idx)):
                global_n = int(level_idx[m_local])
                st = op.streaming_terms(
                    cell_idx=0,
                    direction_idx=m_local,
                    mu_level_idx=level,
                )
                assert st.mu == float(quad.mu_x[global_n])


# ═══════════════════════════════════════════════════════════════════════
# Misuse / coordinate-mismatch guards
# ═══════════════════════════════════════════════════════════════════════

class TestGuards:
    """Misuse raises ``ValueError`` at the factory entry point."""

    @pytest.mark.foundation
    def test_slab_factory_rejects_spherical_mesh(self):
        with pytest.raises(ValueError, match="CARTESIAN"):
            slab_streaming(_spherical_mesh(), Quadrature.gauss_legendre(4))

    @pytest.mark.foundation
    def test_spherical_factory_rejects_slab_mesh(self):
        with pytest.raises(ValueError, match="SPHERICAL"):
            spherical_streaming(_slab_mesh(), Quadrature.gauss_legendre(4))

    @pytest.mark.foundation
    def test_cylindrical_factory_rejects_sphere_mesh(self):
        with pytest.raises(ValueError, match="CYLINDRICAL"):
            cylindrical_streaming(
                _spherical_mesh(), Quadrature.product(2, 4),
            )

    @pytest.mark.foundation
    def test_cylindrical_factory_requires_level_indices(self):
        mesh = _cylindrical_mesh()
        quad = Quadrature.gauss_legendre(8)  # no level_indices
        with pytest.raises(ValueError, match="level structure"):
            cylindrical_streaming(mesh, quad)
