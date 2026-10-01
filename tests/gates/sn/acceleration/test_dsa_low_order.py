r"""The 3b object gates — the production low-order build (D7/D8 + ties).

Three foundation clusters:

1. **The production tie**: :class:`DSALowOrderSystem` (the SN-side
   production build) equals the derivation-side reference builder
   ``orpheus.derivations.discrete.sn.dsa.build_consistent_dd_system``
   entry-for-entry, per group, on the heterogeneous non-uniform slab —
   both BC variants. The derivation is the algebra of record (its rows
   are theorems); a drift in either spelling reds here instead of
   forking silently.
2. **Admission teeth**: the loud seams (non-DD scheme, unsupported
   walls, the Σw = 2 convention boundary, the D-positivity guard)
   actually refuse.
3. **D7/D8 — the R/P object laws**: R conservation (hand-posed
   ``⟨1, R r⟩ = ⟨1, r⟩`` with explicit :math:`w_n, V_i`), the R
   identifications (``integrate_angular`` ≡ the ℓ=0 frame row ≡ the
   displacement tangent map, 0-ULP), and the exact round-trip
   ``R ∘ P = I`` with the frame's :math:`Y_0^0 = 1` table backing the
   normalized-injection spelling of P.

Levels: foundation (object identities; no physics reference).
"""
from __future__ import annotations

import re

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.derivations.discrete.sn import dsa as dsa_reference
from orpheus.geometry import BC, StructuredGeometry
from orpheus.geometry.boundary import (
    AlbedoBoundary,
    PrescribedInflow,
    ReflectiveBoundary,
    SpecularReturn,
    VacuumInflow,
    WhiteBoundary,
)
from orpheus.mesh import CellEdges, CellsByCount, Mesher
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.acceleration import DSACorrection, DSALowOrderSystem
from orpheus.sn.problem import SNProblem
from orpheus.transport.fields.angular_flux import AngularFlux

pytestmark = pytest.mark.foundation


def _four_cell_mesh(left: BC, right: BC):
    """The non-uniform 4-cell slab: materials 0 | 1 | 0 on [0, 0.5, 3, 5],
    the middle region split at 1.5."""
    return Mesher(StructuredGeometry.slab(
        (0.0, 0.5, 3.0, 5.0), (0, 1, 0), left=left, right=right,
    )).partition((
        CellEdges(np.array([0.0, 0.5])),
        CellEdges(np.array([0.5, 1.5, 3.0])),
        CellEdges(np.array([3.0, 5.0])),
    )).mesh


def _four_cell_problem(left, right) -> SNProblem:
    """Heterogeneous, non-uniform 4-cell slab, S4, 2 groups (the 3a tie
    fixture — mixtures carry real P1 data, so the (23c) D is exercised
    beyond the bare-P0 coincidence), with any declared laws, typed or
    tagged."""
    return SNProblem(
        _four_cell_mesh(left, right), Quadrature.gauss_legendre(n_ordinates=4),
        {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")},
    )


def _slab(left: str = "vacuum", right: str = "vacuum") -> SNProblem:
    """The tie fixture with the laws named by their ``BC`` tags."""
    return _four_cell_problem(BC(left), BC(right))


def _reference_inputs(problem: SNProblem):
    """The reference builder's per-group inputs, gathered the same way
    the production build gathers them (data path shared; the FORMULAS
    are what the tie discriminates)."""
    h = np.diff(np.asarray(getattr(problem.mesh, "edges"), float))
    xs = problem.mat_xs
    sigma_t = np.asarray(xs.total_cross_section_field.values, float)
    ng = sigma_t.shape[0]
    mat_ids = np.asarray(problem.mat_map, int).ravel()
    fold = xs.foldable_sigma()
    sigma_s0 = np.stack([fold[int(m)] for m in mat_ids], axis=1)
    residual = xs.residual_sig_s()
    s1 = {
        mid: (
            np.asarray(np.diag(mats[1]), float)
            if len(mats) > 1
            else np.zeros(ng)
        )
        for mid, mats in residual.items()
    }
    sigma_s1 = np.stack([s1[int(m)] for m in mat_ids], axis=1)
    mu = np.asarray(problem.quad.mu_x, float)
    w = np.asarray(problem.quad.weights, float)
    return h, sigma_t, sigma_s0, sigma_s1, mu, w


class TestProductionTie:
    """Weld: the production build ≡ the derivation reference builder."""

    @pytest.mark.parametrize(
        "bc", [("vacuum", "vacuum"), ("reflective", "vacuum")]
    )
    @pytest.mark.verifies("sn-dsa-consistent-low-order")
    def test_low_order_matches_reference_builder(self, bc):
        """At ``scattering_order=1`` so the mixtures' real P1 rows
        exercise the (23c) transport-corrected D, not the bare-P0
        coincidence."""
        problem = _slab(*bc)
        system = DSALowOrderSystem.from_problem(problem.with_scattering_order(1))
        h, sigma_t, sigma_s0, sigma_s1, mu, w = _reference_inputs(problem)
        for g in range(sigma_t.shape[0]):
            a_ref, g_ref = dsa_reference.build_consistent_dd_system(
                h, sigma_t[g], sigma_s0[g], sigma_s1[g], mu, w, bc=bc
            )
            np.testing.assert_allclose(
                system.a_low[g], a_ref, rtol=0, atol=1e-15,
                err_msg=f"A_low group {g} must equal the proven reference",
            )
            np.testing.assert_allclose(
                system.g_map[g], g_ref, rtol=0, atol=1e-15,
                err_msg=f"G group {g} must equal the proven reference",
            )

    def test_solve_correction_realizes_the_reference_solve(self):
        """f0 = A⁻¹(G d) and the (28a) cell average, against a direct
        dense solve of the reference system."""
        problem = _slab()
        system = DSALowOrderSystem.from_problem(problem.with_scattering_order(1))
        h, sigma_t, sigma_s0, sigma_s1, mu, w = _reference_inputs(problem)
        rng = np.random.default_rng(7)
        d0 = rng.standard_normal((sigma_t.shape[0], h.shape[0]))
        f0 = system.solve_correction(d0)
        for g in range(sigma_t.shape[0]):
            a_ref, g_ref = dsa_reference.build_consistent_dd_system(
                h, sigma_t[g], sigma_s0[g], sigma_s1[g], mu, w
            )
            d = np.concatenate([d0[g], np.zeros_like(d0[g])])
            f_ref = np.linalg.solve(a_ref, g_ref @ d)
            np.testing.assert_allclose(f0[g], f_ref, rtol=1e-12, atol=1e-14)
            np.testing.assert_allclose(
                system.cell_update(f0)[g],
                0.5 * (f_ref[:-1] + f_ref[1:]),
                rtol=1e-12, atol=1e-14,
            )


class TestAdmissionTeeth:
    """The loud seams refuse — each guard actually bites."""

    @pytest.mark.parametrize(
        "unadmitted",
        [WhiteBoundary(), AlbedoBoundary(albedo=0.3), PrescribedInflow()],
        ids=["white", "albedo", "prescribed"],
    )
    def test_unsupported_boundary_refused(self, unadmitted):
        """The albedo/white seam guard, on a structural stub carrying the
        admission surface (geometry, scheme, per-face laws). The stub is for
        the bare ``AlbedoBoundary(0.3)`` row: the SN realizer refuses an
        albedo with no re-emission closure before this guard is reached.
        The other two are not pre-refused: a typed law reaches
        :class:`SNProblem` whatever its tag registry lists (since
        ``985497b5``, 2026-08-05), so this guard is the live door for a
        white or prescribed-inflow face (``[M]`` 2026-09-30: both build an
        ``SNProblem`` and are refused here), and for the partial specular
        reflector of ERR-094 (:meth:`test_a_partial_specular_reflector_is_refused`).

        The stub's faces carry REAL laws. Until campaign phase B2 they were
        ``SimpleNamespace(kind="white")`` tag surrogates, which could only
        exercise a string comparison; the guard now reads the laws' affine
        factors, so a surrogate would be testing a fiction. ``PrescribedInflow``
        is in the set deliberately: its response IS zero, so the factor test
        alone would ADMIT it and silently build a Marshak row that drops ``q``
        — it is refused by family, and this leg is what holds that.
        """
        from types import SimpleNamespace

        real = _slab()
        stub = SimpleNamespace(
            is_cartesian=True,
            ndim=1,
            mesh=real.mesh,
            scheme=SimpleNamespace(key="diamond_difference"),
            bc={
                "xmin": SimpleNamespace(law=unadmitted),
                "xmax": SimpleNamespace(law=VacuumInflow()),
            },
        )
        with pytest.raises(NotImplementedError, match="Marshak-albedo"):
            DSALowOrderSystem.from_problem(stub)  # type: ignore[arg-type]

    _PARTIAL = [
        (face, spelling)
        for face in ("xmin", "xmax")
        for spelling in ("reflective", "albedo-specular")
    ]

    @pytest.mark.catches("ERR-094")
    @pytest.mark.parametrize(
        "face, spelling", _PARTIAL, ids=[f"{f}-{s}" for f, s in _PARTIAL],
    )
    def test_a_partial_specular_reflector_is_refused(self, face, spelling):
        r"""A specular law of amplitude 0.7 permutes ordinates like the mirror,
        but its net current is :math:`(1 - \alpha) J^+`, which the mirror's
        low-order row (39), :math:`f_1 = 0`, does not state. It has no proven
        row and is refused, on a real :class:`SNProblem`, on either face.

        First red, ``[M]`` 2026-09-30 with the amplitude condition removed
        (the defect, ERR-094): all four rows build a system. What that system
        did, on a 1-group slab of 40 cells, Gauss-Legendre 8, with
        ``ReflectiveBoundary("x", 0.7)`` on both faces and source iteration
        with DSA: at :math:`c = 0.9`, :math:`\sigma_t h = 1` it converged to the
        right fixed point (rate 0.24); at :math:`c = 0.99`,
        :math:`\sigma_t h = 5` it diverged (residual ``inf`` after 4000
        iterations, flux 8e151 times the plain-SI answer); at
        :math:`\sigma_t h = 20` it raised on a NaN. The mirror's row on a
        partial reflector is the inconsistent low-order system of the
        diffusive regime.
        """
        law = (
            ReflectiveBoundary("x", 0.7) if spelling == "reflective"
            else AlbedoBoundary(0.7, SpecularReturn("x"))
        )
        laws = {"xmin": VacuumInflow(), "xmax": VacuumInflow(), face: law}
        problem = _four_cell_problem(laws["xmin"], laws["xmax"])
        with pytest.raises(
            NotImplementedError,
            match=rf"{re.escape(repr(law))}.*Marshak-albedo",
        ):
            DSALowOrderSystem.from_problem(problem)

    @pytest.mark.rests_on(
        "tests/gates/sn/acceleration/test_dsa_low_order.py::TestAdmissionTeeth::"
        "test_a_partial_specular_reflector_is_refused",
    )
    @pytest.mark.catches("ERR-094")
    def test_the_dsa_entry_refuses_a_partial_reflector(self):
        """The public route: ``solve_sn_fixed_source(..., acceleration="dsa")``
        on a partial reflector refuses before a sweep runs, with the same
        message (the admission is the operator build the entry calls)."""
        from orpheus.sn.solver import solve_sn_fixed_source

        mesh = Mesher(StructuredGeometry.slab(
            (0.0, 5.0), (0,),
            left=ReflectiveBoundary("x", 0.7), right=ReflectiveBoundary("x", 0.7),
        )).partition(CellsByCount.uniform_width(5)).mesh
        with pytest.raises(NotImplementedError, match="Marshak-albedo"):
            solve_sn_fixed_source(
                {0: get_mixture("A", "2g")}, mesh,
                Quadrature.gauss_legendre(n_ordinates=4),
                np.ones((4, 2, 5)), acceleration="dsa",
            )

    _EDGES = [
        ("reflective-0", ReflectiveBoundary("x", 0.0), VacuumInflow()),
        ("albedo-specular-0", AlbedoBoundary(0.0, SpecularReturn("x")), VacuumInflow()),
        ("albedo-specular-1", AlbedoBoundary(1.0, SpecularReturn("x")),
         ReflectiveBoundary("x")),
    ]

    @pytest.mark.rests_on(
        "tests/gates/sn/acceleration/test_dsa_low_order.py::TestProductionTie::"
        "test_low_order_matches_reference_builder",
    )
    @pytest.mark.parametrize(
        "edge, law, foundation", _EDGES, ids=[e[0] for e in _EDGES],
    )
    def test_the_edges_build_the_foundation_rows(self, edge, law, foundation):
        r"""The positive legs of the refusal above: a specular law returning
        nothing builds the vacuum (Marshak) system, and one returning
        everything, spelled as an albedo, builds the mirror's system, both
        bitwise (``[M]`` 2026-09-30: ``a_low`` and ``g_map`` equal to the last
        bit). The foundation laws' systems are tied to the derivation's
        reference builder by the row this rests on, and differ from each
        other (asserted), so equality picks the right one. A guard reading
        the pairing on the geometry tier only (where the albedo spelling
        carries the identity) refuses ``albedo-specular-1`` and reds it
        (``[M]`` 2026-09-30, ``scratch/boundary_ontology/battery_err094.md``).
        """
        built = DSALowOrderSystem.from_problem(_four_cell_problem(law, law))
        expected = DSALowOrderSystem.from_problem(
            _four_cell_problem(foundation, foundation),
        )
        other = DSALowOrderSystem.from_problem(_four_cell_problem(
            *((VacuumInflow(),) * 2 if foundation == ReflectiveBoundary("x")
              else (ReflectiveBoundary("x"),) * 2),
        ))
        if np.array_equal(expected.a_low, other.a_low):
            pytest.fail("the vacuum and mirror systems coincide: the edge rows "
                        "cannot tell which foundation they reproduce")
        np.testing.assert_array_equal(
            built.a_low, expected.a_low,
            err_msg=f"{edge}: the low-order matrix is not the foundation's",
        )
        np.testing.assert_array_equal(
            built.g_map, expected.g_map,
            err_msg=f"{edge}: the low-order source map is not the foundation's",
        )

    def test_non_dd_scheme_refused(self):
        from orpheus.transport.spatial.linear_discontinuous import (
            LinearDiscontinuous,
        )

        mesh1d = _four_cell_mesh(BC("vacuum"), BC("vacuum"))
        problem = SNProblem(
            mesh1d,
            Quadrature.gauss_legendre(n_ordinates=4),
            {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")},
            scheme=LinearDiscontinuous(),
        )
        with pytest.raises(NotImplementedError, match="diamond"):
            DSALowOrderSystem.from_problem(problem)

    def test_quadrature_convention_guard_fires(self):
        problem = _slab()
        h, sigma_t, sigma_s0, sigma_s1, mu, w = _reference_inputs(problem)
        with pytest.raises(ValueError, match="Σw = 2"):
            DSALowOrderSystem._build(
                h, sigma_t, sigma_s0, sigma_s1, mu, w / 2.0,
                (VacuumInflow(), VacuumInflow()),
            )

    def test_d_positivity_guard_fires(self):
        problem = _slab()
        h, sigma_t, sigma_s0, _sigma_s1, mu, w = _reference_inputs(problem)
        bad_s1 = np.full_like(sigma_t, 2.0) * sigma_t  # σ_s1 > σ_t
        with pytest.raises(ValueError, match="positive"):
            DSALowOrderSystem._build(
                h, sigma_t, sigma_s0, bad_s1, mu, w,
                (VacuumInflow(), VacuumInflow()),
            )


class TestRestrictionProlongation:
    """D7/D8 — the R/P object laws on the DSA carrier."""

    @pytest.fixture()
    def psi(self):
        problem = _slab()
        rng = np.random.default_rng(11)
        values = rng.standard_normal((4, 2, 4))
        # CS4b S4: the field no longer carries the mesh — the tests that
        # need carrier data receive the pair.
        return problem, AngularFlux(values=values, space=problem.angular_bulk_space)

    @pytest.mark.verifies("sn-dsa-restriction")
    def test_d7_restriction_conserves_particles(self, psi):
        r"""⟨1, R r⟩ = ⟨1, r⟩ — hand-posed with explicit w_n and V_i
        (structurally independent of the einsum body)."""
        problem, psi = psi
        w = np.asarray(problem.quad.weights, float)
        volumes = np.diff(np.asarray(problem.mesh.edges, float))
        reduced = psi.integrate_angular().values  # (ng, nx)
        lhs = 0.0
        rhs = 0.0
        for g in range(reduced.shape[0]):
            for i in range(reduced.shape[1]):
                lhs += volumes[i] * reduced[g, i]
                for n in range(w.shape[0]):
                    rhs += volumes[i] * w[n] * psi.values[n, g, i]
        np.testing.assert_allclose(lhs, rhs, rtol=1e-14, atol=0)

    @pytest.mark.verifies("sn-dsa-restriction")
    def test_d8_restriction_is_the_frame_moment_row(self, psi):
        r"""``integrate_angular`` ≡ the ℓ=0 analysis row of
        ``Quadrature.angular_frame(0)`` (Y⁰₀ = 1 under the no-prefactor
        SH convention ⟹ the row IS the weight vector) — 0-ULP."""
        problem, psi = psi
        frame = problem.quad.angular_frame(0)
        table = np.asarray(frame.table)  # (N, 1, 1); Y00 ≡ 1
        np.testing.assert_array_equal(table.ravel(), np.ones(4))
        w = np.asarray(problem.quad.weights, float)
        frame_row = np.einsum("n,ng...->g...", w * table.ravel(), psi.values)
        np.testing.assert_array_equal(
            psi.integrate_angular().values, frame_row,
        )

    @pytest.mark.parametrize("n_ordinates", [4, 8])
    def test_d8_prolongation_round_trip_is_identity(self, n_ordinates):
        r"""The P0 injection IS the angular section (#520), and
        R ∘ P = I to the N-term product-sum's rounding.

        The section's divisor is the frame's Gram entry, so DSA follows
        the one convention instead of re-deriving it. GL8 is the witness:
        its Gram entry is ``1.9999999999999998`` against
        ``w.sum() == 2.0``, so the former hand-rolled ``x / w.sum()``
        injection differs from the section there (the equality below is
        red under that spelling), while GL4 cannot tell the two apart.
        R ∘ P is nulp-tier, not bit-exact, whichever divisor is used
        (``[M]`` 2026-09-28: max rel 2.2e-16 at GL8 for both; the same
        re-association ``test_g61_retraction_of_section_is_the_identity``
        records)."""
        mesh1d = _four_cell_mesh(BC("vacuum"), BC("vacuum"))
        problem = SNProblem(
            mesh1d,
            Quadrature.gauss_legendre(n_ordinates=n_ordinates),
            {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")},
        )
        corr = DSACorrection.from_problem(problem)
        section = problem.angular_bulk_space.section("angular")
        np.testing.assert_array_equal(corr._sum_w, section.total_weight)
        phi = np.random.default_rng(0).uniform(0.5, 2.0, (2, 4))
        injection = corr._section.apply(phi)
        np.testing.assert_array_equal(injection, section.apply(phi))
        injected = AngularFlux(values=injection, space=problem.angular_bulk_space)
        # An N-term product-sum re-associates to at most ~N ULP.
        np.testing.assert_allclose(
            injected.integrate_angular().values, phi,
            rtol=n_ordinates * np.finfo(float).eps, atol=0,
        )


class TestApplyAdmission:
    """CS3 §4.5 — the rewritten input guard of ``DSACorrection.apply`` has
    teeth (net-new: [M] no ``pytest.raises`` targeted it before the carve).

    Since the cone carve the SI sweep increment and the Krylov swept vector
    are ONE type (``AngularFlux``), so the guard admits exactly that and
    refuses a moment-windowed carrier by name.
    """

    def _increment(self, problem, seed=5):
        from orpheus.transport.full_field import FullField
        from orpheus.transport.fields.angular_boundary_flux import (
            AngularBoundaryFlux,
        )

        rng = np.random.default_rng(seed)
        return FullField(
            interior=AngularFlux(values=rng.standard_normal(
                    (problem.quad.N, problem.ng, *problem.spatial_shape)
                ), space=problem.angular_bulk_space),
            boundary=AngularBoundaryFlux.zeros(problem.angular_trace),
        )

    def test_moment_windowed_interior_refuses(self):
        """NEGATIVE — a moment-windowed carrier is outside the arm-1
        admission; the message names it (pin the SHORT fragment only)."""
        from orpheus.transport.full_field import FullField
        from orpheus.transport.fields.angular_boundary_flux import (
            AngularBoundaryFlux,
        )
        from orpheus.transport.fields.harmonic_moment_flux import (
            HarmonicMomentFlux,
        )

        problem = _slab()
        corrector = DSACorrection.from_problem(problem)
        L = 1
        # the angular head is READ off the frame (#429): the slab's is FLAT.
        shape = (
            *problem.quad.angular_frame(L).basis.space.shape,
            problem.ng,
            *problem.spatial_shape,
        )
        windowed = FullField(
            interior=HarmonicMomentFlux.from_problem_and_L(
                np.ones(shape), problem, L
            ),
            boundary=AngularBoundaryFlux.zeros(problem.angular_trace),
        )
        with pytest.raises(TypeError, match="moment-windowed"):
            corrector.apply(windowed)

    def test_flux_increment_admitted_and_flux_typed(self):
        """POSITIVE (vv #11 pairing) — a full-angular increment is admitted
        and the correction comes back FLUX-typed on both blocks (the CS3
        re-typing: the update is the plain vector add ψ + Δψ_corr).

        The boundary block is NONZERO here and on vacuum alike ([M] probe
        2026-08-19: ‖trace‖ = 3.665 vacuum / 5.985 reflective) — the trace
        arm always writes the wall-edge f₀ solutions; "inert on vacuum"
        means UNREAD downstream, which this unit test structurally cannot
        see (the consumption claim lives in the acceleration-level gates).
        """
        from orpheus.transport.fields.angular_boundary_flux import (
            AngularBoundaryFlux,
        )

        for bcs in [("vacuum", "vacuum"), ("reflective", "vacuum")]:
            problem = _slab(*bcs)
            corrector = DSACorrection.from_problem(problem)
            out = corrector.apply(self._increment(problem))
            if type(out.interior) is not AngularFlux:
                pytest.fail(
                    f"{bcs}: correction interior is "
                    f"{type(out.interior).__name__}, not AngularFlux"
                )
            if type(out.boundary) is not AngularBoundaryFlux:
                pytest.fail(
                    f"{bcs}: correction boundary is "
                    f"{type(out.boundary).__name__}, not AngularBoundaryFlux"
                )
            if not np.linalg.norm(out.boundary.values) > 0.0:
                pytest.fail(f"{bcs}: the trace arm wrote a zero boundary block")
