r"""#291/#259 — the SN k-estimator cross-engine gate + its mutation teeth.

THE CLAIM UNDER GATE: the reported k is the eigenvalue of the problem
``solve_fixed_source`` poses (only fission scaled by 1/k; scattering and
the (n,2n) emission plain gains) — spelled
``k = R_νΣf / (R_Σa + L − E_2n)`` in :meth:`SNSolver.compute_keff`.

GROUND TRUTH: the converged-fixed-point **map ratio**
``k* = P(Mφ*) / P(φ*)`` with ``M: φ ↦ solve_fixed_source(Fφ)`` — at the
dominant eigenvector ANY positive linear functional gives the
eigenvalue. Structurally independent of the estimator's spelling: it
consumes only the iteration map and a linear functional, never the
absorption/leakage/n2n decomposition under test. (The SN analogue of
the diffusion ``TestEngineCrossGate`` — no vacuum absolute-k anchor
existed in tests/sn before this gate; every prior k assertion consumed
the estimator itself.)

Characterized pre-fix defects this gate REDs on (2026-07-03, pre-fix
@ ``d1daaac``): vacuum slab reported 1.8377 vs true 0.9816 (**+87%**);
het vacuum sphere +22.5%; reflective Σ₂≠0 reported 1.9286 vs true
2.6128 (**−26%**, zero leakage — the R7 convention defect). Post-fix
agreement ≤ 6e-10 at these tolerances.

TEETH (vv discipline — a gate's bite is proven, permanently): the
leakage-drop mutation REDs the vacuum leg and leaves the reflective leg
**bitwise** green (that asymmetry IS the #291 catcher — reflective
faces are a structural zero, never a computed one); the leakage
sign-flip REDs the vacuum leg; the old-convention (n,2n)-in-numerator
mutation REDs the Σ₂≠0 leg; the d=3 transverse-enumeration mutants
(reversed / transposed) red the 3-D vacuum leg, whose face-measure
ARRAY is additionally object-pinned against the boundary layer's
``volumes / Δ_axis`` (Mode-12 discipline). All mutations are
in-process monkeypatches (never a git-checkout revert —
process-discipline rule).

``foundation`` — a solver-invariant claim (estimator ≡ posed problem),
pinned by an exact self-reference, not a literature value.
"""
from __future__ import annotations

import copy

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import (
    BC, CoordSystem, Mesh1D, Region, RegionMesh, StructuredGeometry,
)
from orpheus.numerics.eigenvalue import power_iteration
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.solver import SNSolver, _as_sn_mesh
from orpheus.transport.mesh.axis import AxisMesh
from orpheus.transport.reaction_rate_functional import IntegratedReactionRate

pytestmark = pytest.mark.foundation

_COORD_TAG = {"CARTESIAN": "SLB", "CYLINDRICAL": "CYL", "SPHERICAL": "SPH"}


def _mesh(regions, bc, coord):
    """regions = [(mat_id, thickness_cm, n_cells), ...]"""
    tag = _COORD_TAG[coord.name]
    bcs = (bc, bc) if tag == "SLB" else (bc,)
    geom = StructuredGeometry(
        geometry=tag,
        regions=tuple(
            Region(mat_id=m, outer_thickness_cm=t) for m, t, _ in regions
        ),
        bcs=bcs,
    )
    return Mesh1D.from_geometry(
        geom,
        region_meshes=tuple(RegionMesh(n_cells=n) for *_, n in regions),
    )


def _solve(materials, mesh, scattering_order=0):
    """Converged eigenpair via the production driver, solver retained."""
    sn_mesh = _as_sn_mesh(mesh, Quadrature.gauss_legendre(8), materials)
    solver = SNSolver(
        sn_mesh,
        scattering_order=scattering_order,
        keff_tol=1e-9, flux_tol=1e-8, max_inner=2000, inner_tol=1e-11,
    )
    keff, _history, phi = power_iteration(solver, max_iter=2000)
    return solver, keff, phi


def _map_ratio_kstar(solver, phi):
    """One application of the power-iteration map at the fixed point.

    ``k* = P(Mφ*)/P(φ*)`` — P is the (linear) total production
    functional; at the converged eigenvector the ratio is the dominant
    eigenvalue of the posed problem regardless of the estimator's
    decomposition.
    """
    q = solver.compute_fission_source(phi, 1.0)
    phi_next = solver.solve_fixed_source(q, phi.copy())
    return solver.compute_production_rate(phi_next) / \
        solver.compute_production_rate(phi)


def _mixture_with_n2n():
    """Library fuel + an asymmetric injected Sig2 (the #269 fixture idiom).

    Every xs_library mixture has Sig2 = 0, so the (n,2n) convention leg
    would be vacuously green on library data (vv Mode 10); the asymmetry
    also catches a g_from↔g_to swap.
    """
    from scipy.sparse import csr_matrix

    fuel = copy.deepcopy(get_mixture("A", "2g"))
    sig2 = np.array([[0.0, 0.03], [0.01, 0.0]])
    fuel.Sig2 = csr_matrix(sig2)
    fuel.SigT = np.asarray(fuel.SigT) + sig2.sum(axis=1)
    return fuel


# The d=3 vacuum box for the transverse-area-product arm (#291 3-D wire):
# distinct per-axis cell counts AND distinct, strongly non-uniform width
# patterns (x increasing, y decreasing, z mixed — Mode-2 asymmetry), so a
# transposed or reversed transverse enumeration moves the leakage sum O(1)
# instead of cancelling on a symmetric coincidence.
_D3_EDGES = (
    np.array([0.0, 0.3, 0.8, 1.6, 3.0]),   # x: 4 cells, widths ↑
    np.array([0.0, 1.2, 1.9, 2.4]),        # y: 3 cells, widths ↓
    np.array([0.0, 0.9, 1.4]),             # z: 2 cells
)


def _d3_vacuum_sn_mesh():
    vac = BC("vacuum")
    axes = tuple(
        AxisMesh(edges=e, bc_low=vac, bc_high=vac) for e in _D3_EDGES
    )
    return SNMesh.from_axes(
        axes, Quadrature.level_symmetric(sn_order=4),
        {0: get_mixture("A", "2g")},
    )


def _d3_vacuum_solve():
    """Converged d=3 all-vacuum eigenpair, solver retained."""
    solver = SNSolver(
        _d3_vacuum_sn_mesh(),
        keff_tol=1e-9, flux_tol=1e-8, max_inner=2000, inner_tol=1e-11,
    )
    keff, _history, phi = power_iteration(solver, max_iter=2000)
    return solver, keff, phi


# ═══════════════════════════════════════════════════════════════════════
# The cross-engine gate: reported k ≡ the posed problem's eigenvalue.
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.verifies("sn-keff-update")
class TestReportedKeffIsThePosedEigenvalue:
    """reported k ≡ map-ratio k* across the BC × Σ₂ defect matrix."""

    RTOL = 5e-8  # measured post-fix agreement ≤ 6e-10; ~80× headroom

    @pytest.mark.verifies("sn-leakage-functional")
    @pytest.mark.catches("ERR-064")
    def test_vacuum_slab(self):
        """[#291] Vacuum-bounded slab — the leakage term is load-bearing.

        Pre-fix this leg reads 1.8377 vs true 0.9816 (+87%): the
        absorption-only denominator reports a near-critical leaky slab
        as strongly supercritical.
        """
        solver, keff, phi = _solve(
            {0: get_mixture("A", "2g")},
            _mesh([(0, 8.0, 40)], BC.vacuum, CoordSystem.CARTESIAN),
        )
        k_star = _map_ratio_kstar(solver, phi)
        np.testing.assert_allclose(keff, k_star, rtol=self.RTOL)

    @pytest.mark.verifies("sn-leakage-functional")
    @pytest.mark.catches("ERR-064")
    def test_vacuum_sphere_curvilinear_face_area(self):
        """[#291] Vacuum sphere — pins the 4πR² face-area convention.

        The leakage functional integrates the trace's net current over
        ``MaterialMesh.areas`` — the SAME geometric convention the cell
        volumes integrate under. A missing/extra 4π (vv Mode 3) moves
        the denominator's leakage term ~an order of magnitude and REDs
        this leg.
        """
        solver, keff, phi = _solve(
            {0: get_mixture("A", "2g")},
            _mesh([(0, 4.0, 20)], BC.vacuum, CoordSystem.SPHERICAL),
        )
        k_star = _map_ratio_kstar(solver, phi)
        np.testing.assert_allclose(keff, k_star, rtol=self.RTOL)

    @pytest.mark.verifies("sn-leakage-functional")
    def test_vacuum_3d_box_transverse_area_product(self):
        """[#291 d=3 arm] Vacuum box — the transverse-area outer product.

        The d=3 leakage term integrates each face's net current over
        the outer product of the two transverse axes' edge widths in
        ascending-axis order (the ``face_shape`` enumeration). k* (map
        ratio) never consumes the face areas, so on the Mode-2
        asymmetric box any transposed / reversed / wrong-axis
        enumeration moves the reported k off k* (teeth proven in
        ``TestGateTeeth``). This is also the not-slow WIRING pin for
        the d=3 vacuum eigenvalue class (the ERR-069 discipline —
        before the wire, this path raised ``NotImplementedError`` and
        only a slow-tier gate could see it).
        """
        solver, keff, phi = _d3_vacuum_solve()
        k_star = _map_ratio_kstar(solver, phi)
        np.testing.assert_allclose(keff, k_star, rtol=self.RTOL)

    def test_reflective_lattice_bitwise_degenerate(self):
        """Reflective ⟹ the estimator IS the lattice functional, bitwise.

        Reflective faces are a STRUCTURAL zero in the leakage loop
        (skipped, never accumulated) and Σ₂ = 0 makes the emission an
        exact 0.0 — so the reported k must equal the historical
        ``production/absorption`` arithmetic EXACTLY (float ``==``, not
        allclose): ``P/(A + 0.0 − 0.0)`` is the same division. This is
        the pin that keeps every reflective anchor and regression
        snapshot in the suite unmoved by #291.
        """
        solver, keff, phi = _solve(
            {0: get_mixture("A", "2g")},
            _mesh([(0, 8.0, 40)], BC.reflective, CoordSystem.CARTESIAN),
        )
        production = IntegratedReactionRate(
            solver.mat_xs.fission_production_field
        ).evaluate(phi)
        absorption = IntegratedReactionRate(
            solver.mat_xs.absorption_cross_section_field
        ).evaluate(phi)
        assert keff == production / absorption, (
            f"reflective Σ₂=0 reported k={keff!r} must be BITWISE the "
            f"lattice functional {production / absorption!r} — the "
            f"leakage/n2n terms must be structural zeros, not computed "
            f"±noise."
        )
        # And the lattice functional is the posed eigenvalue here too.
        np.testing.assert_allclose(
            keff, _map_ratio_kstar(solver, phi), rtol=self.RTOL,
        )

    @pytest.mark.catches("ERR-065")
    def test_reflective_n2n_convention(self):
        """[R7] Σ₂≠0, zero leakage — the (n,2n) placement is load-bearing.

        Pre-fix this leg reads 1.9286 vs true 2.6128 (−26%): the
        (νΣf+2Σ₂)/Σa spelling equals the posed problem's eigenvalue only
        when Σ₂=0 or k=1. The unified estimator puts the emission on the
        net-removal side and lands on k* for any Σ₂.
        """
        solver, keff, phi = _solve(
            {0: _mixture_with_n2n()},
            _mesh([(0, 8.0, 40)], BC.reflective, CoordSystem.CARTESIAN),
        )
        k_star = _map_ratio_kstar(solver, phi)
        np.testing.assert_allclose(keff, k_star, rtol=self.RTOL)


# ═══════════════════════════════════════════════════════════════════════
# The face-measure OBJECT pin (vv Mode 12: pin the object, not only the
# functional through which the k gate consumes it).
# ═══════════════════════════════════════════════════════════════════════


class TestFaceAreaObjectPin:
    """``_face_area_of`` ≡ boundary-layer ``volumes / Δ_axis``, all faces.

    The k gate above sees the face measure only through the leakage
    functional; this pin fixes the ARRAY itself against an
    independently-assembled oracle — the mesh's own cell volumes
    divided by the face axis's boundary width. The shaped ``volumes``
    array carries the mesh's ascending-axis cell enumeration, which is
    exactly the ``face_shape`` enumeration the trace layout uses, so
    agreement here localizes any ordering defect to the exact face
    array instead of a k residual.
    """

    def test_d3_face_measure_matches_volume_ratio_on_all_six_faces(self):
        sn = _d3_vacuum_sn_mesh()
        solver = SNSolver(sn)
        volumes = np.asarray(sn.volumes, dtype=float).reshape(
            sn.spatial_shape
        )
        widths = [np.diff(e) for e in _D3_EDGES]
        for axis, prefix in enumerate("xyz"):
            for face, layer in ((f"{prefix}min", 0), (f"{prefix}max", -1)):
                area = np.asarray(solver._face_area_of(face))
                # Trace-layout authority: the slot is (N, ng, *face_spatial).
                slot = sn.angular_trace.layout.faces[face]
                np.testing.assert_equal(
                    area.shape, slot.shape[2:],
                    err_msg=f"{face}: face measure vs trace-slot spatial "
                            f"shape",
                )
                oracle = (
                    np.take(volumes, layer, axis=axis) / widths[axis][layer]
                )
                np.testing.assert_allclose(
                    area, oracle, rtol=1e-13,
                    err_msg=f"{face}: face measure vs boundary-layer "
                            f"volumes/Δ oracle",
                )

    def test_d2_face_measure_unchanged_by_the_d3_wire(self):
        """The d=2 arm is the single-axis degenerate of the same body.

        ``reduce(outer, [w])`` IS ``w`` — assert the 2-D face measure
        still equals the transverse width vector bitwise, so the #291
        2-D estimator arithmetic is untouched by the generalization.
        """
        vac = BC("vacuum")
        axes = (
            AxisMesh(edges=_D3_EDGES[0], bc_low=vac, bc_high=vac),
            AxisMesh(edges=_D3_EDGES[1], bc_low=vac, bc_high=vac),
        )
        sn = SNMesh.from_axes(
            axes, Quadrature.level_symmetric(sn_order=4),
            {0: get_mixture("A", "2g")},
        )
        solver = SNSolver(sn)
        np.testing.assert_array_equal(
            solver._face_area_of("xmin"), np.diff(_D3_EDGES[1]),
        )
        np.testing.assert_array_equal(
            solver._face_area_of("ymax"), np.diff(_D3_EDGES[0]),
        )


# ═══════════════════════════════════════════════════════════════════════
# Mutation teeth — the gate's bite, proven in-process and permanently.
# ═══════════════════════════════════════════════════════════════════════


class TestGateTeeth:
    """Seeded estimator defects MUST red the gate exactly where claimed."""

    def test_leakage_drop_reds_vacuum_stays_bitwise_on_reflective(
        self, monkeypatch,
    ):
        """The #291 catcher asymmetry: drop L ⟹ vacuum REDs, reflective exact.

        With ``_boundary_leakage_rate`` muted to 0.0 the vacuum leg must
        re-develop an O(L/A) ≈ 87% bias against k* — and the reflective
        leg must stay BITWISE identical to the unmutated solve, proving
        the reflective arithmetic never touched the leakage path at all.
        """
        # Reflective reference BEFORE mutation (clean process state).
        _, k_refl_clean, _ = _solve(
            {0: get_mixture("A", "2g")},
            _mesh([(0, 8.0, 40)], BC.reflective, CoordSystem.CARTESIAN),
        )

        with monkeypatch.context() as m:
            m.setattr(
                SNSolver, "_boundary_leakage_rate",
                lambda self, fission_production: 0.0,
            )
            solver, k_vac_mut, phi = _solve(
                {0: get_mixture("A", "2g")},
                _mesh([(0, 8.0, 40)], BC.vacuum, CoordSystem.CARTESIAN),
            )
            k_star = _map_ratio_kstar(solver, phi)
            rel_gap = abs(k_vac_mut / k_star - 1.0)
            assert rel_gap > 0.05, (
                f"leakage-drop mutation must re-open the #291 bias on the "
                f"vacuum leg (expected ~0.87 relative), got {rel_gap:.3e} "
                f"— the gate has no teeth against the omission it exists "
                f"to catch."
            )

            _, k_refl_mut, _ = _solve(
                {0: get_mixture("A", "2g")},
                _mesh([(0, 8.0, 40)], BC.reflective, CoordSystem.CARTESIAN),
            )
            assert k_refl_mut == k_refl_clean, (
                f"reflective reported k moved under the leakage-drop "
                f"mutation ({k_refl_mut!r} vs {k_refl_clean!r}) — the "
                f"reflective path is supposed to be a STRUCTURAL zero "
                f"that never evaluates the leakage functional."
            )

    def test_leakage_sign_flip_reds_vacuum(self, monkeypatch):
        """A sign-flipped L (A−L denominator) must red the vacuum leg.

        Detection is accepted in EITHER mode — the flipped denominator
        ``A−L ≈ 0.07`` drives a huge-k / negative-k transient whose
        fission source flips the iterate negative, so the observed
        outcome is the scale-bridge guard RAISING on the degenerate
        state (a crash-RED, the loudest catch); if a variant mutant
        survived to convergence instead, the reported k must disagree
        with k* O(1). The only forbidden outcome is a mutant that
        solves AND matches k*.
        """
        true_rate = SNSolver._boundary_leakage_rate

        with monkeypatch.context() as m:
            m.setattr(
                SNSolver, "_boundary_leakage_rate",
                lambda self, p: -true_rate(self, p),
            )
            try:
                solver, k_mut, phi = _solve(
                    {0: get_mixture("A", "2g")},
                    _mesh([(0, 8.0, 40)], BC.vacuum, CoordSystem.CARTESIAN),
                )
            except RuntimeError:
                return  # the guard caught the mutant mid-iteration — RED
            k_star = _map_ratio_kstar(solver, phi)
            assert abs(k_mut / k_star - 1.0) > 0.05, (
                "a sign-flipped leakage term must move the vacuum "
                "reported k O(1) — the gate failed to catch A−L."
            )

    def test_old_convention_reds_n2n_leg(self, monkeypatch):
        """The retired (νΣf+2Σ₂)/Σa spelling must red the Σ₂≠0 leg (R7)."""

        def old_convention_keff(self, flux_distribution):
            production = self.compute_production_rate(flux_distribution)
            absorption = IntegratedReactionRate(
                self.mat_xs.absorption_cross_section_field
            ).evaluate(flux_distribution)
            return production / absorption

        with monkeypatch.context() as m:
            m.setattr(SNSolver, "compute_keff", old_convention_keff)
            solver, k_mut, phi = _solve(
                {0: _mixture_with_n2n()},
                _mesh([(0, 8.0, 40)], BC.reflective, CoordSystem.CARTESIAN),
            )
            k_star = _map_ratio_kstar(solver, phi)
            assert abs(k_mut / k_star - 1.0) > 0.05, (
                "the old n2n-in-numerator convention must disagree with "
                "the posed eigenvalue on a Σ₂≠0 problem (measured −26% "
                "pre-fix) — the Σ₂≠0 leg failed to catch it."
            )

    def test_transverse_ordering_mutants_red_the_3d_vacuum_leg(
        self, monkeypatch,
    ):
        """d=3 face-measure enumeration mutants must be caught.

        Two mutants of ``_face_area_of`` on the Mode-2 asymmetric box:

        * REVERSED — the first transverse axis enumerated descending
          (``area[::-1]``): shape-silent, value-wrong on the
          non-uniform widths; the reported k must leave the k* band.
        * TRANSPOSED — the transverse axes swapped (``area.T``): with
          distinct per-axis counts the broadcast against
          ``(ng, *face_spatial)`` cannot align, so the accepted
          detection is a loud crash-RED; a square-face survivor would
          have to miss k*. The only forbidden outcome is a mutant that
          solves AND matches k*.
        """
        true_area = SNSolver._face_area_of

        with monkeypatch.context() as m:
            m.setattr(
                SNSolver, "_face_area_of",
                lambda self, face: true_area(self, face)[::-1],
            )
            solver, k_mut, phi = _d3_vacuum_solve()
            k_star = _map_ratio_kstar(solver, phi)
            gap = abs(k_mut / k_star - 1.0)
            assert gap > 100 * TestReportedKeffIsThePosedEigenvalue.RTOL, (
                f"the reversed transverse enumeration must move the d=3 "
                f"reported k off the posed eigenvalue (got relative gap "
                f"{gap:.3e}) — the 3-D vacuum leg has no teeth against "
                f"face-cell ordering defects."
            )

        with monkeypatch.context() as m:
            m.setattr(
                SNSolver, "_face_area_of",
                lambda self, face: true_area(self, face).T,
            )
            try:
                solver, k_mut, phi = _d3_vacuum_solve()
            except ValueError:
                return  # broadcast refuses the swapped enumeration — RED
            k_star = _map_ratio_kstar(solver, phi)
            assert abs(k_mut / k_star - 1.0) > \
                100 * TestReportedKeffIsThePosedEigenvalue.RTOL, (
                    "a transposed transverse enumeration survived the "
                    "broadcast AND matched k* — the d=3 vacuum leg has "
                    "no teeth against the axis swap."
                )
