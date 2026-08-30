"""Factory-binding + packet-plumbing tests for the reduced streaming operator.
What this file pins
-------------------
The **wiring**: that :class:`~orpheus.sn.mesh.augmented_mesh.SNMesh`
builds ``self.reduced`` by calling this module's chart factories
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

Moved here 2026-08-28 (un-weld arc P4.4)
----------------------------------------
This file was ``tests/geometry/test_reduced_operator.py`` until the
system under test moved from ``orpheus/geometry/reduced_operator.py`` to
:mod:`orpheus.sn.mesh.reduced_operator`.  ⚠ The "all 47 tests in this
file stay green" measurement above was taken 2026-08-03 on the file as it
then stood; the α-dome closure contract (7 tests) stayed behind in
``tests/geometry/test_reduced_operator.py`` with :func:`~orpheus.sn.angular.redistribution.alpha_dome`,
so this file now carries 29 of those cases.  The mutation has **not** been
re-run since — the blindness finding is unchanged in kind, but the count
is historical.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.transport.spatial.scheme import StreamingTerms
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.mesh.reduced_operator import (
    ReducedStreamingOperator,
    cylindrical_streaming,
    slab_streaming,
    spherical_streaming,
)
from orpheus.sn.angular.redistribution import angular_redistribution
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
        assert np.array_equal(reduced.face_areas, sn_mesh.reduced.face_areas)

    @pytest.mark.foundation
    def test_delta_A_read_through_is_the_factory_value(self, pair):
        """The deprecated ``SNMesh.delta_A`` lands on the factory array."""
        sn_mesh, reduced = pair
        assert reduced.delta_A is not None
        assert np.array_equal(reduced.delta_A, sn_mesh.reduced.delta_A)

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
        snm_reduced = sn_mesh.reduced
        assert snm_reduced is not None  # 1-D mesh => minted by the ctor (narrowing)
        reduced = spherical_streaming(mesh, quad)
        assert np.array_equal(
            reduced.angular.alpha_per_level[0],
            snm_reduced.angular.alpha_per_level[0],
        )
        assert np.array_equal(reduced.delta_A, snm_reduced.delta_A)


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
        assert np.array_equal(reduced.face_areas, sn_mesh.reduced.face_areas)

    @pytest.mark.foundation
    def test_delta_A_read_through_is_the_factory_value(self, pair):
        """The deprecated ``SNMesh.delta_A`` lands on the factory array."""
        sn_mesh, reduced = pair
        assert np.array_equal(reduced.delta_A, sn_mesh.reduced.delta_A)

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
        snm_reduced = sn_mesh.reduced
        assert snm_reduced is not None  # 1-D mesh => minted by the ctor (narrowing)
        reduced = cylindrical_streaming(mesh, quad)
        for rdc, snm in zip(
            reduced.angular.alpha_per_level,
            snm_reduced.angular.alpha_per_level,
        ):
            assert np.array_equal(rdc, snm)
        assert np.array_equal(reduced.delta_A, snm_reduced.delta_A)


# ═══════════════════════════════════════════════════════════════════════
# Property tests
# ═══════════════════════════════════════════════════════════════════════

class TestProperties:
    """ReducedStreamingOperator advertises the right metadata.

    ⭐ **What the three chart tests became, and why (P4.1a, 2026-08-27).**
    Each used to read ``assert op.coord is CoordSystem.X`` -- a stored
    literal read straight back, i.e. a tautology waiting for the field to
    be recognised as a copy.  It was: every factory *validates*
    ``mesh.coord`` against the literal it then stored, so
    ``op.coord is op.mesh.coord`` held **by construction**, and the field
    retired.

    That retirement PROMOTED the three validation guards from redundant to
    load-bearing -- they are now the sole reason ``op.mesh.coord`` is the
    operator's chart -- and ``[M]`` 2026-08-27 they had **zero witnesses**
    anywhere in the tree (`grep "requires .* mesh"` returned only the three
    production lines).  That is the ``coding-standards`` MIRROR clause: a
    retirement silently strengthening a claim raises no alarm, because the
    suite only gets greener.  So each test below is now the guard's
    witness, in ``vv-principles`` #11 form -- a **positive** leg (the right
    chart is accepted, and the operator reports it) and a **negative** leg
    per wrong chart (each other chart is refused, by name).
    """

    @pytest.mark.foundation
    def test_slab_is_posed_on_the_cartesian_chart(self):
        """Until 2026-08-26 this test also asserted two flags that were
        each exactly ``coord is not CARTESIAN`` -- one line below this
        very assertion.  They had no production reader and are retired;
        the fact they pinned is pinned here."""
        quad = Quadrature.gauss_legendre(8)
        op = slab_streaming(_slab_mesh(), quad)
        assert op.mesh.coord is CoordSystem.CARTESIAN

        for wrong in (_spherical_mesh(), _cylindrical_mesh()):
            with pytest.raises(ValueError, match="requires CARTESIAN mesh"):
                slab_streaming(wrong, quad)

    @pytest.mark.foundation
    def test_slab_carries_the_neutral_element_on_BOTH_factors(self):
        """Slab has no curvature — and it SPELLS that on both factors now,
        rather than leaving either one ``None``.

        ⭐ This is a STRICTLY STRONGER claim than the one this test made
        before (``alpha_half is None`` etc.).  ``None`` says only "the
        sphere/cylinder arm did not populate this"; the neutral element
        says the dome IS identically zero and the starting direction IS
        the diameter ray — i.e. it pins the VALUES a slab must carry so
        that the angular closure's ``c_in = c_out = 0`` fall out of the
        general body instead of being special-cased.  The per-coordinate
        ``Optional`` union died with it (Pattern 4).

        ⛔ Until P4.1b (2026-08-27) the SPATIAL half was still exempt: this
        test asserted ``op.face_areas is None`` / ``op.delta_A is None``
        under the comment *"The SPATIAL chart is still absent for a slab"*,
        two lines below the paragraph above celebrating the same upgrade on
        the ANGULAR half.  It was not absent — it was underived.
        ``compute_areas_1d`` returns a real unit cross-section on CARTESIAN,
        so the slab's face areas are ``ones(nx+1)`` and its connection
        integral is ``zeros(nx)``: **a slab has no AREA CHANGE, which is
        what "no curvature" means.**  Spelling it is what let
        ``streaming_terms``' three arms collapse to one body — the slab arm
        was the sphere body with those two values written out by hand.
        """
        mesh = _slab_mesh()
        quad = Quadrature.gauss_legendre(8)
        op = slab_streaming(mesh, quad)
        # The SPATIAL factor is present and NEUTRAL: unit cross-section,
        # zero area change.  These are the values the collapsed body reads
        # where the retired CARTESIAN arm carried the literals 1.0 and 0.0.
        assert np.array_equal(op.face_areas, np.ones(mesh.N + 1))
        assert np.array_equal(op.delta_A, np.zeros(mesh.N))
        # The ANGULAR factor is present and NEUTRAL.
        assert op.angular.n_levels == 1
        assert np.array_equal(
            op.angular.alpha_per_level[0], np.zeros(quad.N + 1)
        )
        assert op.angular.mu_start_per_level == (-1.0,)

    @pytest.mark.foundation
    def test_sphere_is_posed_on_the_spherical_chart(self):
        quad = Quadrature.gauss_legendre(8)
        op = spherical_streaming(_spherical_mesh(), quad)
        assert op.mesh.coord is CoordSystem.SPHERICAL

        for wrong in (_slab_mesh(), _cylindrical_mesh()):
            with pytest.raises(ValueError, match="requires SPHERICAL mesh"):
                spherical_streaming(wrong, quad)

    @pytest.mark.foundation
    def test_cylinder_is_posed_on_the_cylindrical_chart(self):
        quad = Quadrature.product(n_mu=2, n_phi=4)
        op = cylindrical_streaming(_cylindrical_mesh(), quad)
        assert op.mesh.coord is CoordSystem.CYLINDRICAL

        for wrong in (_slab_mesh(), _spherical_mesh()):
            with pytest.raises(ValueError, match="requires CYLINDRICAL mesh"):
                cylindrical_streaming(wrong, quad)


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
        # (P4.7: the packet's ``chord_length`` / ``mu`` / ``delta_A_over_w``
        # retired — zero production readers; the ΔA/w fusion is
        # closure-owned and cache-interned, its zero pinned at the factor.)
        # Neutral curvature: slab carries the values that make the
        # unified cell-balance algebra collapse to the slab form.
        #
        # ⭐ These three assertions were PROMOTED by P4.1b (2026-08-27)
        # without a character of the test body changing.  Until then the
        # CARTESIAN arm of ``streaming_terms`` *returned the literals*
        # ``1.0`` / ``1.0`` / ``0.0`` and this restated them — a tautology.
        # That arm is retired; the values now come from ``mesh.areas``
        # through the one shared body, so these lines test that the face-area
        # path reaches the packet.  A wrong ``compute_areas_1d`` on CARTESIAN
        # reddens here (the producer itself is pinned independently, against
        # a hand-written literal, at ``tests/geometry/test_geometry.py``
        # ``TestSurfaces1D::test_cartesian``).
        assert st.face_area_inner == 1.0
        assert st.face_area_outer == 1.0
        assert float(op.delta_A[2]) == 0.0

    @pytest.mark.foundation
    def test_sphere_streaming_terms_match_arrays(self):
        mesh = _spherical_mesh()
        quad = Quadrature.gauss_legendre(8)
        op = spherical_streaming(mesh, quad)
        i, n = 1, 5
        st = op.streaming_terms(cell_idx=i, direction_idx=n)
        # Geometric labels: inner = closer to r=0 (A[i]),
        #                   outer = farther (A[i+1]).  (P4.7: the packet's
        # chord/mu/ΔA-w retired; abs_mu pins the ordinate resolution.)
        assert st.face_area_inner == float(op.face_areas[i])
        assert st.face_area_outer == float(op.face_areas[i + 1])
        assert st.abs_mu == float(abs(quad.mu_x[n]))

    @pytest.mark.foundation
    def test_cylinder_streaming_terms_match_per_level(self):
        mesh = _cylindrical_mesh()
        quad = Quadrature.product(n_mu=4, n_phi=4)
        op = cylindrical_streaming(mesh, quad)
        i, level, m = 2, 1, 0
        st = op.streaming_terms(
            cell_idx=i, direction_idx=m, mu_level_idx=level,
        )
        # Geometric labels: inner / outer relative to r=0,
        # independent of sweep direction.
        assert st.face_area_inner == float(op.face_areas[i])
        assert st.face_area_outer == float(op.face_areas[i + 1])
        # The GLOBAL-ordinate resolution (bug 2 fix anchor) now rides
        # ``abs_mu`` — same ``level_indices`` path the retired signed-mu
        # and ΔA/w packings used (P4.7; the cache's own gates pin the
        # interned ΔA/w row against the factors at this resolution).
        global_n = int(quad.level_indices[level][m])
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

# ═══════════════════════════════════════════════════════════════════════
# Geometric labels: face_area_inner/outer are direction-independent
# ═══════════════════════════════════════════════════════════════════════

class TestStreamingTermsGeometricLabels:
    """``face_area_inner`` / ``face_area_outer`` are **purely geometric**.

    The two face-area fields encode position relative to :math:`r=0`
    (inner = closer, outer = farther) — they do NOT depend on the
    sweep's marching direction.  The same cell yields the same
    inner / outer values regardless of which ordinate is queried.
    Sweep-direction resolution is stamped by the sweep producer onto
    :class:`orpheus.transport.spatial.scheme.CellVisit`.
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
        # Direction-dependence lives on the quadrature (P4.7 — the
        # signed packet field retired).
        assert float(quad.mu_x[n_neg]) < 0
        assert float(quad.mu_x[n_pos]) > 0

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
        # Signed eta is direction-dependent (off the quadrature, P4.7).
        assert float(quad.mu_x[level_idx[m_neg]]) < 0
        assert float(quad.mu_x[level_idx[m_pos]]) > 0


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

    # (P4.7: ``test_cylindrical_signed_mu_per_level`` retired with the
    # packet's signed ``mu`` — the GLOBAL-ordinate resolution claim it
    # pinned lives on in the ``abs_mu`` sibling above, which walks the
    # same ``level_indices`` path.)


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
# P4-remainder (2026-08-29): the producer binds the angular AXIS
# ═══════════════════════════════════════════════════════════════════════

def _scale_decoy(quad):
    """The measured keystone decoy: nodes x 0.9, weights preserved.

    [M] the ONLY decoy admissible at BOTH tiers — but the refusals are
    TWO contracts, not one (archivist re-measure 2026-08-29, correcting
    this docstring's first wording): the α-dome guard (Σ w·µ = 0,
    `redistribution.py`) refuses only the ROLL (it breaks the pairing to
    ±0.3); negated/reversed nodes PASS the dome (±5.6e-17) and die one
    tier later at the closure's P3 τ ∈ [0,1] guard (`closure.py`).
    Scaling preserves both. Weight-preserving ⟹ the decoy axis is
    IDENTITY-EQUAL to the true one (the space cannot see it) while every
    nonzero |cosine| moves.
    """
    import dataclasses as _dc

    dm = quad.measure
    decoy_m = _dc.replace(dm, nodes=(np.asarray(dm.nodes) * 0.9).copy())
    return Quadrature(
        measure=decoy_m,
        level_structure=quad.level_structure,
        folded_by=quad.folded_by,
    )


def _level_roll_decoy(quad):
    """K2's isolator: SAME measure, level list rolled by one.

    Moves only what the ``level_indices`` read resolves — ``mu_x`` is
    untouched — so it separates the ``:517`` read from the ``:528`` read
    (vv #17's per-arm discipline; the scale decoy moves both together).
    """
    import dataclasses as _dc

    ls = quad.level_structure
    assert ls is not None
    rolled = _dc.replace(
        ls,
        level_indices=tuple(
            ls.level_indices[(k - 1) % len(ls.level_indices)]
            for k in range(len(ls.level_indices))
        ),
    )
    return Quadrature(
        measure=quad.measure, level_structure=rolled, folded_by=quad.folded_by
    )


class TestP4RemTheProducerBindsTheAxis:
    """The binding's own gates: the two mints agree, the courier is dead,
    the generator-less axis refuses by name, and — the KEYSTONE — the
    packet reads go THROUGH the axis (route gates; every value comparison
    over the re-point is X == X because the bound generator IS the
    courier's old payload, three-way ``is``-identical)."""

    @pytest.mark.foundation
    @pytest.mark.parametrize(
        "coord,quad_factory",
        [
            (CoordSystem.CARTESIAN, lambda: Quadrature.gauss_legendre(4)),
            (CoordSystem.SPHERICAL, lambda: Quadrature.gauss_legendre(4)),
            (CoordSystem.CYLINDRICAL, lambda: Quadrature.folded_product(4, 6)),
        ],
    )
    def test_the_two_mints_agree_on_the_1d_arm(self, coord, quad_factory):
        """Gate 4.3 (1-D arm) — the factory's axis IS the space's axis.

        ``==`` deliberately, never ``is`` (the mint is fresh per call);
        plus the generator identity and the field-not-property row. This
        is the ONLY gate that reds on a wrong LABEL at either mint site
        (the label defaults at the generator so a twin spelling is
        unspellable, and this row is its witness).
        """
        from tests.sn._test_helpers import placeholder_materials

        quad = quad_factory()
        mesh = Mesh1D(
            edges=np.linspace(0.0, 1.0, 5),
            mat_ids=np.zeros(4, dtype=int),
            coord=coord,
        )
        from orpheus.sn.mesh.augmented_mesh import SNMesh

        sn = SNMesh(mesh, quad, placeholder_materials())
        r = sn.reduced
        assert r is not None
        assert r.angular_axis == sn.angular_bulk_space.axis("angular")
        assert r.angular_axis.generator is sn.quad
        assert r.angular_axis is r.angular_axis  # field, not per-access mint

    @pytest.mark.foundation
    def test_the_two_mints_agree_on_the_d2_cartesian_arm(self):
        """Gate 4.3 (d ≥ 2 arm) — the RSO-less branch's own mint agrees.

        [M] the branch M1's corpus cannot redden (verification plan §2.1):
        ``reduced is None`` there, so the hub mints the closure's axis
        itself — the second mint site, one label typo away from the
        first until the generator-side default killed the spelling. The
        arm's premise (``reduced is None``) is asserted as its own row.
        """
        from tests.sn._test_helpers import placeholder_materials
        from orpheus.sn.mesh.augmented_mesh import SNMesh
        from orpheus.transport.mesh.axis import AxisMesh

        axes = tuple(
            AxisMesh(edges=np.linspace(0.0, ext, n + 1))
            for ext, n in zip((1.0, 2.0), (2, 3))
        )
        sn = SNMesh.from_axes(
            axes, Quadrature.level_symmetric(sn_order=4),
            placeholder_materials(ng=2),
        )
        assert sn.reduced is None  # the arm's own premise
        # The closures consume the axis at construction and do not store
        # it (a stored-but-unread field is dead weight, not provenance),
        # so the arm's mint is witnessed by: the two mint EXPRESSIONS
        # agreeing (the label-twin witness — both spell ``quad.axis()``
        # with the label defaulted at the generator), and the narrow's
        # observable consequence (the closure's ordinate count followed
        # the axis's generator).
        assert sn.quad.axis() == sn.angular_bulk_space.axis("angular")
        assert sn.quad.axis().generator is sn.quad
        from orpheus.sn.angular.closure import IdentityAngularClosure

        closure = sn.angular_closure
        assert isinstance(closure, IdentityAngularClosure)
        assert closure.tau_per_ordinate.size == sn.quad.N

    @pytest.mark.foundation
    def test_the_courier_is_dead_by_field_set(self):
        """Gate 4.2 — a STRUCTURAL row, not a value row.

        ``dataclasses.fields`` (never ``hasattr`` — a defaulted field or
        a ``getattr`` fallback would still answer), and EQUALITY of the
        set so a re-addition reds too. This is the row that keeps a
        partial re-point unspellable: with no courier on the angular
        factor, a consumer that still wants the quadrature has exactly
        one place to get it — through the axis.
        """
        import dataclasses

        from orpheus.sn.angular.redistribution import AngularRedistribution

        assert tuple(
            f.name for f in dataclasses.fields(AngularRedistribution)
        ) == ("coord", "alpha_per_level", "mu_start_per_level")

    @pytest.mark.foundation
    def test_a_generator_less_axis_refuses_naming_streaming_terms(self):
        """G5.2 (producer row) — the refusal names BOTH parties.

        Two-fragment match: the axis label AND the consumer's name — a
        generic message keeps a wrong reason true (L31). Fragment
        disjointness (G5.4): the neighbouring space refusal spells
        'axis lookup:', which this message deliberately does not.
        """
        from orpheus.numerics.axis import Axis, BasisKind

        quad = Quadrature.gauss_legendre(4)
        bare = Axis(
            "angular", (quad.N,),
            weights=np.asarray(quad.weights, float), kind=BasisKind.NODAL,
        )
        op = ReducedStreamingOperator(
            mesh=_slab_mesh(),
            angular=angular_redistribution(quad, CoordSystem.CARTESIAN),
            angular_axis=bare,
        )
        with pytest.raises(ValueError, match="angular"):
            op.streaming_terms(0, 0)
        with pytest.raises(ValueError, match="streaming_terms"):
            op.streaming_terms(0, 0)
        with pytest.raises(ValueError, match="minted through"):
            op.streaming_terms(0, 0)

    @pytest.mark.foundation
    @pytest.mark.parametrize(
        "coord,quad_factory,moved_min",
        [
            (CoordSystem.CARTESIAN, lambda: Quadrature.gauss_legendre(4), 4),
            (CoordSystem.SPHERICAL, lambda: Quadrature.gauss_legendre(4), 4),
            (CoordSystem.CYLINDRICAL, lambda: Quadrature.folded_product(4, 6), 8),
        ],
    )
    def test_K1_the_packet_reads_THROUGH_the_axis(
        self, coord, quad_factory, moved_min
    ):
        """KEYSTONE K1 — the route gate, because every value gate is X == X.

        [M] the bound generator IS the courier's old payload (three-way
        ``is``-identity), so only a decoy that the space cannot see can
        prove the route. Four legs (verification plan §5.3): the anti-dud
        control, the route (a decoy generator MOVES |µ|), the decoy's
        invisibility to space identity (vv #19 — the instrument's own
        precondition), and the AR-was-not-reached documentation leg.
        Red-before record: pre-carve, the §5.2 simulation (decoy on a
        true-dome AR) is the honest analogue — [M] slab 4/4, sphere 4/4,
        cylinder 8/12 packets moved; this row pins those same counts.
        """
        quad = quad_factory()
        mesh = Mesh1D(
            edges=np.linspace(0.0, 1.0, 5),
            mat_ids=np.zeros(4, dtype=int),
            coord=coord,
        )
        ar_true = angular_redistribution(quad, coord)
        decoy = _scale_decoy(quad)

        op_true = ReducedStreamingOperator(
            mesh=mesh, angular=ar_true, angular_axis=quad.axis()
        )
        op_decoy = ReducedStreamingOperator(
            mesh=mesh, angular=ar_true, angular_axis=decoy.axis()
        )
        # leg 3: the decoy is INVISIBLE to space identity
        assert quad.axis() == decoy.axis()
        assert hash(quad.axis()) == hash(decoy.axis())
        # leg 4: the AR was NOT reached (nothing on it to reach — the
        # discriminating power comes from the axis alone)
        assert op_decoy.angular is ar_true

        def _packets(op):
            out = []
            if coord is CoordSystem.CYLINDRICAL:
                ls = quad.level_structure
                assert ls is not None
                for p_idx, level in enumerate(ls.level_indices):
                    for d_idx in range(len(level)):
                        out.append(
                            op.streaming_terms(1, d_idx, mu_level_idx=p_idx).abs_mu
                        )
            else:
                for d_idx in range(quad.N):
                    out.append(op.streaming_terms(1, d_idx).abs_mu)
            return np.asarray(out)

        base = _packets(op_true)
        # leg 1: the anti-dud control — the true axis reproduces the
        # factory-built operator's packets exactly
        factory = {
            CoordSystem.CARTESIAN: slab_streaming,
            CoordSystem.SPHERICAL: spherical_streaming,
            CoordSystem.CYLINDRICAL: cylindrical_streaming,
        }[coord]
        np.testing.assert_array_equal(base, _packets(factory(mesh, quad)))
        # leg 2: the ROUTE — the decoy moves the answer
        moved = int(np.sum(~np.isclose(_packets(op_decoy), base)))
        if moved < moved_min:
            pytest.fail(
                f"{coord.name}: only {moved} packets moved under the decoy "
                f"(expected >= {moved_min}) — the producer did not read "
                f"through the axis"
            )

    @pytest.mark.foundation
    def test_K2_the_cylinder_index_read_is_a_separate_route(self):
        """KEYSTONE K2 — the level-roll decoy isolates the ``level_indices``
        read from the ``mu_x`` read (same measure, rolled level list —
        ``mu_x`` untouched). The scale decoy moves both reads together,
        so without this row a partial re-point is indistinguishable from
        a complete one.

        [M] 2026-08-29, this gate's own fixture: **4 of 12** packets move
        on fp(4,6) — the level |µ| sequence is a PALINDROME ([A,B,B,A]),
        so roll-1 fixes half the levels, and each level's degenerate
        0-|µ| member never moves (the L43e palindrome family, partially
        biting the roll as it fully bites the reversal). 4 > 0 with
        ``mu_x`` untouched is exactly the isolation this row needs.
        """
        quad = Quadrature.folded_product(4, 6)
        mesh = Mesh1D(
            edges=np.linspace(0.0, 1.0, 5),
            mat_ids=np.zeros(4, dtype=int),
            coord=CoordSystem.CYLINDRICAL,
        )
        ar_true = angular_redistribution(quad, CoordSystem.CYLINDRICAL)
        decoy = _level_roll_decoy(quad)
        assert quad.axis() == decoy.axis()  # same measure ⟹ invisible

        op_true = ReducedStreamingOperator(
            mesh=mesh, angular=ar_true, angular_axis=quad.axis()
        )
        op_decoy = ReducedStreamingOperator(
            mesh=mesh, angular=ar_true, angular_axis=decoy.axis()
        )
        ls = quad.level_structure
        assert ls is not None

        def _packets(op):
            return np.asarray([
                op.streaming_terms(1, d_idx, mu_level_idx=p_idx).abs_mu
                for p_idx, level in enumerate(ls.level_indices)
                for d_idx in range(len(level))
            ])

        moved = int(np.sum(~np.isclose(_packets(op_decoy), _packets(op_true))))
        if moved < 4:
            pytest.fail(
                f"only {moved} packets moved under the level-roll decoy "
                f"(expected >= 4, the measured floor on this palindromic "
                f"rule) — the level_indices read did not go through the "
                f"axis"
            )
