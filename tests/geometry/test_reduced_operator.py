"""Factory-binding + packet-plumbing tests for the reduced streaming operator.

What this file pins
-------------------
The **wiring**: that :class:`~orpheus.sn.mesh.augmented_mesh.SNMesh`
builds ``self.reduced`` by calling the geometry-layer factories
(:func:`slab_streaming`, :func:`spherical_streaming`,
:func:`cylindrical_streaming`) with the caller's mesh and quadrature —
so every array a curvilinear sweep reads is the one *this* module
produced, and not a divergent second copy that a future refactor
re-introduces.  Plus the per-(cell, direction) ``streaming_terms``
plumbing: which array lands in which packet slot, and that the
geometric ``inner``/``outer`` labels stay direction-independent.

What this file is blind to — MEASURED, 2026-08-03
-------------------------------------------------
The connection-coefficient **math**.  The file's original claim — the
factories are bit-identical to the legacy in-line
``SNMesh._setup_spherical`` / ``._setup_cylindrical`` — was a genuine
two-implementation gate *while both were live*.  Those methods were
retired at the Wave-B carve (#159), and ``SNMesh.__init__`` now calls
these very factories, so the comparison silently became
``factory(mesh, quad) == SNMesh(mesh, quad).reduced``: the SAME
producer on both sides.  Measured by replacing every array the
factories emit with deterministic garbage (in-process rebind of all
12 module bindings, so both call sites receive the SAME wrong values,
exactly as a real sign flip in ``reduced_operator.py`` would):

* **all 47 tests in this file stay green**;
* 29 gates in ``tests/sn/primitives/test_quadrature.py`` +
  ``tests/sn/sweep/curvilinear/test_alpha_closed_form.py`` redden in
  the same run (the positive control, ``vv-principles`` anti-#17);
* garbaging only the ``SNMesh``-side binding reddens exactly the 15
  factory-binding cases below — which is the claim they honestly carry.

``face_areas`` is not merely equal but the *same object* on both sides
(both are the shared ``Mesh1D.areas`` cached property), so that leg is
``array_equal(x, x)`` for any face-area math whatsoever.

Where the math IS pinned — measured by the same mutation
--------------------------------------------------------
Every array below has at least one **structurally-independent**
(closed-form) catcher; the regression snapshots corroborate but are
nowhere the sole evidence.

* ``delta_A`` — ``tests/sn/primitives/test_quadrature.py::
  TestL0TermVerification::test_delta_A_magnitude`` (closed form
  :math:`4\\pi\\,\\Delta(r^2)` / :math:`2\\pi\\,\\Delta r`).
  ⛔ This entry used to add "sole catcher — the snapshots are blind to
  it, which is correct: ``delta_A`` has no production consumer."  That
  was true until 2026-08-26 and is now FALSE: retiring the fused
  ``redist_dAw`` cache made ``delta_A`` the SPATIAL factor that both
  ``streaming_terms`` and the angular closure read, so every
  curvilinear snapshot now rides on it.
* ``angular.alpha_per_level`` (was ``alpha_half`` / ``alpha_per_level``)
  — ``…::test_per_ordinate_flat_flux_consistency`` on both arms, the L0
  per-ordinate flat-flux identity (``catches("ERR-006", "ERR-007")``);
  ``tests/sn/sweep/curvilinear/test_alpha_closed_form.py`` (the
  Dirichlet-kernel closed form); plus both geometries' snapshots.
* ``redist_dAw`` / ``redist_dAw_per_level`` — **RETIRED 2026-08-26** as
  a fused product (``ΔA ⊗ 1/w``) that neither of its two consumers
  owned.  Its former catcher — ``tests/sn/sweep/curvilinear/
  test_streaming_equilibrium_curvilinear.py``, the L0 closed-form
  :math:`\\varphi = Q/(\\Sigma_t(1-c))` gate — still covers the
  QUANTITY, which is now formed at each consumer from ``delta_A`` and
  the measure's weights.  (The historical note that the L0 flat-flux
  identity "RECOMPUTES ``dA / w`` rather than reading the production
  array" is why that gate needed no migration: it was already computing
  the product rather than reading the cache.)
* ``face_areas`` — ``tests/geometry/test_geometry.py`` pins the
  producer ``compute_areas_1d`` against the closed form; the
  snapshots pin that the factory forwards it.

These tests are tagged ``@pytest.mark.foundation`` — software
invariants (wiring, packet plumbing), never equation-level math claims.

Wave B Round 1 — Issue 6 (.claude/plans/sn_reshape.md, GH #155);
re-scoped 2026-08-03 after the demoted claim was measured.
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
    alpha_dome,
    cylindrical_streaming,
    slab_streaming,
    spherical_streaming,
)
from orpheus.geometry.reduced_operator import (
    _ALPHA_CLOSURE_ATOL,
    _assert_alpha_dome_closes,
)
from orpheus.sn.mesh.augmented_mesh import SNMesh
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
# Factory-binding tests — SNMesh routes to THIS module's factories
# ═══════════════════════════════════════════════════════════════════════

class TestSNMeshBindsSphericalFactory:
    """Sphere: ``SNMesh.reduced`` is what :func:`spherical_streaming` built.

    The claim is **routing**, not math: re-running the factory on the
    same ``(mesh, quadrature)`` reproduces every array ``SNMesh`` holds,
    so ``SNMesh`` cannot have grown a private second implementation.  A
    change to the connection-coefficient math moves BOTH sides and is
    invisible here by construction — see the module docstring for the
    measurement and for where the math is actually pinned.
    """

    @pytest.fixture
    def pair(self):
        mesh = _spherical_mesh()
        quad = Quadrature.gauss_legendre(8)
        sn_mesh = SNMesh(mesh, quad, placeholder_materials())
        reduced = spherical_streaming(mesh, quad)
        return sn_mesh, reduced

    @pytest.mark.foundation
    def test_face_areas_read_through_is_the_factory_value(self, pair):
        """The deprecated ``SNMesh.face_areas`` lands on the factory array.

        Note both sides are the SAME ``Mesh1D.areas`` object, so this
        leg cannot see a change in how face areas are computed — only a
        read-through that stops forwarding to ``self.reduced``.  The
        stronger ``is``-identity form of that claim lives at
        ``tests/sn/primitives/test_snmesh_consumes_reduced.py``.
        """
        sn_mesh, reduced = pair
        assert reduced.face_areas is not None
        assert np.array_equal(reduced.face_areas, sn_mesh.face_areas)

    @pytest.mark.foundation
    def test_delta_A_read_through_is_the_factory_value(self, pair):
        """The deprecated ``SNMesh.delta_A`` lands on the factory array."""
        sn_mesh, reduced = pair
        assert reduced.delta_A is not None
        assert np.array_equal(reduced.delta_A, sn_mesh.delta_A)

    @pytest.mark.foundation
    def test_angular_factor_is_the_factory_value(self, pair):
        """The ANGULAR factor on ``SNMesh`` came from this module's producer.

        Successor of ``test_alpha_half_is_the_factory_value`` and
        ``test_redist_dAw_is_the_factory_value`` (2026-08-26 un-weld).  The
        routing claim is unchanged; what moved is the object it routes.
        ``redist_dAw`` had no successor because it was never one quantity —
        it was the fused product ``ΔA ⊗ 1/w``, and each of its two factors
        is routed by a test of its own (``delta_A`` above; the weights ride
        the quadrature both sides receive)."""
        sn_mesh, reduced = pair
        assert np.array_equal(
            reduced.angular.alpha_per_level[0],
            sn_mesh.reduced.angular.alpha_per_level[0],
        )
        assert (reduced.angular.mu_start_per_level
                == sn_mesh.reduced.angular.mu_start_per_level)

    # Issue #236 Step C: the geometry-side tau_mm producer was retired (the
    # M-M angular weight is now closure-owned).  ``test_tau_mm_bit_identical``
    # (geometry-vs-geometry τ bit-identity) is subsumed by the closure
    # producer-equivalence floor at
    # ``tests/sn/sweep/curvilinear/test_tau_producer_equivalence.py`` (closure
    # τ vs the structurally-independent contamination.morel_montry_weights).

    @pytest.mark.foundation
    @pytest.mark.parametrize("N", [4, 8, 16, 32])
    def test_factory_binding_holds_across_quadrature_order(self, N):
        """The routing claim is not an artefact of one quadrature order."""
        mesh = _spherical_mesh()
        quad = Quadrature.gauss_legendre(N)
        sn_mesh = SNMesh(mesh, quad, placeholder_materials())
        reduced = spherical_streaming(mesh, quad)
        assert np.array_equal(
            reduced.angular.alpha_per_level[0],
            sn_mesh.reduced.angular.alpha_per_level[0],
        )
        assert np.array_equal(reduced.delta_A, sn_mesh.reduced.delta_A)


class TestSNMeshBindsCylindricalFactory:
    """Cylinder: ``SNMesh.reduced`` is what :func:`cylindrical_streaming` built.

    Routing claim only — see :class:`TestSNMeshBindsSphericalFactory`
    and the module docstring.

    **Why ``folded_product`` and not ``product``** (Q5.6): a cylindrical
    :class:`SNMesh` admits only a rule whose every μ-level is CARRYING
    (``assert_carrying_quadrature``, the R12a march-start predicate).  An
    unquotiented ``product`` rule puts an ordinate at ξ = 0 on the most
    inward node of level 0, so the seed slot is a rank-duplicate of ψ₀ and
    admission REFUSES it.  The σ_y quotient is the fix the guard itself
    names.  The routing claim under test is indifferent to which rule is
    used — both sides receive the same ``quad`` — so this is a fixture
    repair, not a change of claim.  Admissibility of all three shapes below
    is pinned independently by
    ``tests/sn/mesh/test_cylindrical_quadrature_admission.py``.
    """

    @pytest.fixture
    def pair(self):
        mesh = _cylindrical_mesh()
        quad = Quadrature.folded_product(n_mu=2, n_phi=4)
        sn_mesh = SNMesh(mesh, quad, placeholder_materials())
        reduced = cylindrical_streaming(mesh, quad)
        return sn_mesh, reduced

    @pytest.mark.foundation
    def test_face_areas_read_through_is_the_factory_value(self, pair):
        """The deprecated ``SNMesh.face_areas`` lands on the factory array."""
        sn_mesh, reduced = pair
        assert np.array_equal(reduced.face_areas, sn_mesh.face_areas)

    @pytest.mark.foundation
    def test_delta_A_read_through_is_the_factory_value(self, pair):
        """The deprecated ``SNMesh.delta_A`` lands on the factory array."""
        sn_mesh, reduced = pair
        assert np.array_equal(reduced.delta_A, sn_mesh.delta_A)

    @pytest.mark.foundation
    def test_angular_factor_is_the_factory_value(self, pair):
        """The per-level ANGULAR factor came from this module's producer.

        Successor of the per-level α and ΔA/w routing tests (2026-08-26
        un-weld) — see the spherical twin for why ``redist_dAw_per_level``
        has no successor of its own."""
        sn_mesh, reduced = pair
        assert (reduced.angular.n_levels
                == sn_mesh.reduced.angular.n_levels)
        for lvl, (rdc, snm) in enumerate(zip(
            reduced.angular.alpha_per_level,
            sn_mesh.reduced.angular.alpha_per_level,
        )):
            assert np.array_equal(rdc, snm), f"level {lvl} mismatch"
        assert (reduced.angular.mu_start_per_level
                == sn_mesh.reduced.angular.mu_start_per_level)

    # Issue #236 Step C: ``test_tau_mm_per_level_bit_identical`` retired —
    # the geometry-side cylinder τ producer was deleted; the closure τ
    # producer-equivalence floor (``test_tau_producer_equivalence.py``) is
    # the structurally-independent successor.

    @pytest.mark.foundation
    @pytest.mark.parametrize("n_mu,n_phi", [(2, 4), (4, 4), (4, 8)])
    def test_factory_binding_holds_across_quadrature_shape(self, n_mu, n_phi):
        """The routing claim is not an artefact of one quadrature shape."""
        mesh = _cylindrical_mesh()
        quad = Quadrature.folded_product(n_mu=n_mu, n_phi=n_phi)
        sn_mesh = SNMesh(mesh, quad, placeholder_materials())
        reduced = cylindrical_streaming(mesh, quad)
        for rdc, snm in zip(
            reduced.angular.alpha_per_level,
            sn_mesh.reduced.angular.alpha_per_level,
        ):
            assert np.array_equal(rdc, snm)
        assert np.array_equal(reduced.delta_A, sn_mesh.reduced.delta_A)


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
    def test_slab_carries_the_neutral_angular_element(self):
        """Slab has no curvature — and since the 2026-08-26 un-weld it
        SPELLS that, rather than leaving the angular fields ``None``.

        ⭐ This is a STRICTLY STRONGER claim than the one this test made
        before (``alpha_half is None`` etc.).  ``None`` says only "the
        sphere/cylinder arm did not populate this"; the neutral element
        says the dome IS identically zero and the starting direction IS
        the diameter ray — i.e. it pins the VALUES a slab must carry so
        that the angular closure's ``c_in = c_out = 0`` fall out of the
        general body instead of being special-cased.  The per-coordinate
        ``Optional`` union died with it (Pattern 4)."""
        mesh = _slab_mesh()
        quad = Quadrature.gauss_legendre(8)
        op = slab_streaming(mesh, quad)
        # The SPATIAL chart is still absent for a slab (no curvature).
        assert op.face_areas is None
        assert op.delta_A is None
        # The ANGULAR factor is present and NEUTRAL.
        assert op.angular.n_levels == 1
        assert np.array_equal(
            op.angular.alpha_per_level[0], np.zeros(quad.N + 1)
        )
        assert op.angular.mu_start_per_level == (-1.0,)

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
        cell-update strategies branched on ``delta_A_over_w is None``.
        Step 2.5 retires that branch by populating neutral values:
        ``face_area_inner = face_area_outer = 1.0``,
        ``delta_A_over_w = 0.0``.  These values make the unified
        cell-balance helper consume the same packet for slab and
        curvilinear without geometry dispatch.  (Issue #236 Step C
        retired the M-M ``alpha_in`` / ``alpha_out`` / ``tau_mm``
        packing on ``StreamingTerms`` — the angular weight is now
        closure-owned, stamped on ``CellVisit``.)
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
        assert st.delta_A_over_w == float(op.delta_A[i] / quad.weights[n])
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
        # Cylinder mu is signed eta from the GLOBAL ordinate index
        # (resolved via level_indices) — bug 2 fix anchor.
        global_n = int(quad.level_indices[level][m])
        # ΔA/w is formed from its two factors (the fused ``redist_dAw``
        # cache retired 2026-08-26).  The claim class is unchanged: this
        # pins the INDEXING — cell i against GLOBAL ordinate n, not the
        # within-level m, and not transposed.
        assert st.delta_A_over_w == float(
            op.delta_A[i] / quad.weights[global_n]
        )
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
    ``SNMesh`` or the ``Quadrature``.
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

        Re-checks that ``delta_A_over_w`` still matches the underlying
        array after the ``volume`` / ``abs_mu`` extension.  Defense
        against accidental reordering or drift.  (Issue #236 Step C
        retired the M-M ``alpha_in`` / ``alpha_out`` / ``tau_mm``
        packing on ``StreamingTerms`` — the angular weight is now
        closure-owned, stamped on ``CellVisit``.)
        """
        mesh = _spherical_mesh()
        quad = Quadrature.gauss_legendre(8)
        op = spherical_streaming(mesh, quad)
        i, n = 1, 5
        st = op.streaming_terms(cell_idx=i, direction_idx=n)
        assert st.delta_A_over_w == float(op.delta_A[i] / quad.weights[n])


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
    (:class:`orpheus.transport.spatial.scheme.CellVisit`).
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


# ═══════════════════════════════════════════════════════════════════════
# The α-dome: ONE recursion, and a contract that survives ``-O``
# ═══════════════════════════════════════════════════════════════════════


class TestAlphaDomeClosureContract:
    r"""The admission contract :math:`\alpha_{M+1/2} = 0`, on both arms.

    ⛔ **What these rows exist to prevent, measured 2026-08-12.** The
    recursion :math:`\alpha_{m+1/2} = \alpha_{m-1/2} - w_m\mu_m` is
    strictly ONE-SIDED — :math:`\alpha_{1/2} = 0` is an axiom, the far-end
    closure is *not*; it is a consequence of the measure's antisymmetry
    and therefore a real contract on every admitted quadrature. Before
    this batch it was enforced:

    * on the SPHERE by a **bare** ``assert abs(alpha[N]) < 1e-12``;
    * on the CYLINDER **not at all**.

    And ``.claude/rules/vv-testing.md`` makes ``python -O -m pytest``
    canonical, which strips every ``assert``. `[M]` a measure closing at
    ``alpha[N] = +0.2`` was REFUSED under plain ``python`` and **ACCEPTED**
    under ``python -O`` — so the one check that existed did not run in the
    suite that matters.

    ⭐ **These rows are therefore only meaningful because they use
    ``pytest.raises``, not ``assert``.** A future contributor who
    "simplifies" the guard back to an ``assert`` reddens them under the
    canonical runner; that is the whole point of the class.
    """

    #: Every curvilinear rule the tree ships through these factories.
    _SPHERE_N = (4, 8, 16, 32, 64, 128)
    _CYL_RULES = ((4, 8), (4, 16), (4, 32), (8, 16), (6, 12), (2, 6))

    # ---- positive leg: the shipped rules are admitted --------------------
    # (`vv-principles` #11 — a negative-only guard test validates the
    # RAISING behaviour but never the claim that the contract is right.)

    @pytest.mark.foundation
    @pytest.mark.parametrize("N", _SPHERE_N)
    def test_every_shipped_gauss_legendre_dome_closes(self, N):
        """`[M]` residual 5.6e-17 (N=4) … 8.2e-17 (N=64) — ≈1 ULP of the
        dome peak, and NOT drifting with N."""
        quad = Quadrature.gauss_legendre(N)
        alpha = alpha_dome(quad.mu_x, quad.weights)
        _assert_alpha_dome_closes(alpha, coord=CoordSystem.SPHERICAL)
        assert abs(alpha[-1]) < 1e-15, (
            f"GL{N} dome residual {alpha[-1]!r} is far above the ~1e-16 "
            f"floor these rules have always delivered — the guard's 1e-12 "
            f"admission band would still pass it, so this row is the "
            f"tighter statement about the PRODUCERS (`vv-principles` #16: "
            f"the type permits 1e-12; Gauss-Legendre achieves 1e-16)."
        )

    @pytest.mark.foundation
    @pytest.mark.parametrize("n_mu,n_phi", _CYL_RULES)
    def test_every_shipped_folded_product_dome_closes_on_every_level(
        self, n_mu, n_phi,
    ):
        """`[M]` worst residual over the shipped rules is 2.78e-16 at
        ``folded_product(4, 32)``; ``(2, 6)`` closes at exactly 0.0."""
        quad = Quadrature.folded_product(n_mu=n_mu, n_phi=n_phi)
        for p, idx in enumerate(quad.level_indices):
            alpha = alpha_dome(quad.mu_x[idx], quad.weights[idx])
            _assert_alpha_dome_closes(
                alpha, coord=CoordSystem.CYLINDRICAL, level=p,
            )
            assert abs(alpha[-1]) < 1e-15, (
                f"folded_product({n_mu}, {n_phi}) level {p} residual "
                f"{alpha[-1]!r} is far above the ~1e-16 floor"
            )

    # ---- negative leg: an inadmissible dome is REFUSED, under -O ---------

    @pytest.mark.foundation
    @pytest.mark.parametrize("coord,level", [
        (CoordSystem.SPHERICAL, None),
        (CoordSystem.CYLINDRICAL, 3),
    ])
    def test_a_dome_that_does_not_close_is_refused(self, coord, level):
        """The row that the retired ``assert`` could not carry: this must
        raise under ``python -O`` as well as plain ``python``."""
        bad = np.zeros(5)
        bad[-1] = 0.2
        with pytest.raises(ValueError, match="does not close"):
            _assert_alpha_dome_closes(bad, coord=coord, level=level)

    @pytest.mark.foundation
    def test_the_cylinder_message_names_the_offending_level(self):
        """A folded rule closes on every level or none, but a
        level-symmetric rule can lose antisymmetry on ONE level — a
        whole-measure check would only say "somewhere"."""
        bad = np.zeros(5)
        bad[-1] = 0.2
        with pytest.raises(ValueError, match="on level 2"):
            _assert_alpha_dome_closes(
                bad, coord=CoordSystem.CYLINDRICAL, level=2,
            )

    @pytest.mark.foundation
    def test_the_guard_admits_exactly_up_to_its_own_stated_tolerance(self):
        """`vv-principles` #16 — the gate quotes the guard's OWN threshold
        rather than a tighter one the guard does not promise."""
        ok, bad = np.zeros(3), np.zeros(3)
        ok[-1] = 0.5 * _ALPHA_CLOSURE_ATOL
        bad[-1] = 2.0 * _ALPHA_CLOSURE_ATOL
        _assert_alpha_dome_closes(ok, coord=CoordSystem.SPHERICAL)  # no raise
        with pytest.raises(ValueError, match="does not close"):
            _assert_alpha_dome_closes(bad, coord=CoordSystem.SPHERICAL)

    # ---- the single-source claim ----------------------------------------

    @pytest.mark.foundation
    def test_the_derivations_alpha_dome_IS_the_production_one(self):
        r"""Cardinal Rule 2 — the recursion had THREE spellings until
        2026-08-12 (sphere factory, cylinder factory, and the derivations
        analysis module), which is why the contract could be enforced on
        one arm only.

        ⚠ This row is **not** a two-implementation agreement gate and must
        never be described as one: the derivations name now delegates, so
        no input can make the two disagree (``coding-standards``' "single-
        sourcing a duplicate demotes every gate that compared its copies"
        — the demotion is CORRECT here and the fix stays). What it pins is
        that the delegation is still in place, i.e. that a fourth copy has
        not been re-introduced under the old name.
        """
        from orpheus.derivations.discrete.sn.angular_differencing import (
            alpha_dome as derivations_alpha_dome,
        )

        quad = Quadrature.gauss_legendre(8)
        np.testing.assert_array_equal(
            derivations_alpha_dome(quad.mu_x, quad.weights),
            alpha_dome(quad.mu_x, quad.weights),
        )

    @pytest.mark.foundation
    def test_the_derivations_recursion_stays_UNGUARDED(self):
        r"""The contract lives at the production ADMISSION point, not
        inside the recursion — deliberately.

        ``derivations.discrete.sn.angular_differencing``'s P0/P4 predicate
        ladder exists precisely to characterise measures whose dome does
        NOT close. A guard welded into the shared recursion would make
        that analysis unspellable, so the split into
        :func:`alpha_dome` (pure) + :func:`_assert_alpha_dome_closes`
        (admission) is load-bearing, not stylistic.
        """
        from orpheus.derivations.discrete.sn.angular_differencing import (
            alpha_dome as derivations_alpha_dome,
        )

        # A deliberately non-antisymmetric measure: the recursion must
        # still RETURN it (open dome and all), not raise.
        mu = np.array([-0.8, -0.3, 0.3, 0.8])
        w = np.array([0.50, 0.25, 0.25, 0.25])
        alpha = derivations_alpha_dome(mu, w)
        assert alpha[-1] == pytest.approx(0.2), (
            "the analysis recursion must report the OPEN dome, not refuse it"
        )
        # …and the production admission point must refuse that same dome.
        with pytest.raises(ValueError, match="does not close"):
            _assert_alpha_dome_closes(alpha, coord=CoordSystem.SPHERICAL)
