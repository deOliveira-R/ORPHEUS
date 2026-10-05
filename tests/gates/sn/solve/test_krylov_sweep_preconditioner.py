r"""The within-group Krylov solve preconditioned by the sweep (#200): its acceptance gates.

Since #200, :func:`orpheus.sn.solver._within_group_krylov` hands GMRES the
sweep, the operator ``LC.inverse()`` (the exact inverse of the within-group
``M`` of the splitting, block back-substitution on a carrying mesh), as its
left preconditioner when no DSA corrector is posed; before, it handed it the
identity. The preconditioner is the public ``KrylovAcceleration.preconditioner``
operator, which the caller states (it has no default since 2026-10-04). Three
claims, one per acceptance item of the issue:

1. **The preconditioner is the full-space inverse, and linear on the full typed
   state.** GMRES feeds the preconditioner RESIDUAL vectors, which populate
   every block: bulk, inflow and outflow trace rows, and on a carrying mesh the
   ψ½ System-B blocks. The issue records the failure this guards: a sweep whose
   outflow ignored the residual's boundary block converged to a wrong answer
   at machine-zero residual (9953 iterations, 0.57–0.86 against 1.875), because
   ``M·A`` was rank-deficient. That map is LINEAR (a projection followed by a
   sweep), so linearity alone is blind to it (designed-green, ``vv-principles``
   mode 12); the row that sees it is the round trip of a boundary-only
   residual, ``A(M q_b) = q_b``. The rows here gate the ROUTE (the production
   preconditioner is, bit for bit, the inverse whose full-space identity
   ``tests/gates/sn/operators/test_sweep_inverse_identity.py`` proves), the
   linearity the issue names, and the boundary round trip.
2. **The fixed point is unchanged.** A preconditioner changes GMRES's
   trajectory, never its solution: the preconditioned and the unpreconditioned
   solves agree on k and on the scalar flux to ten times the inner tolerance,
   on a slab and a cylinder whose faces are reflective (a live boundary block).
   The independent references are elsewhere and rest on the same route:
   ``test_krylov_curvilinear_precond_safety.py`` (k_inf, closed form, every
   coordinate system) and ``test_l1_standoff_slab_cylinder.py`` (the Case
   reference).
3. **The rate.** The largest inner iteration count of a solve stays flat as the
   mesh refines, where the identity's grows with the number of unknowns.
4. **One preconditioner for both inner postures under DSA.** With a
   consistent-DSA corrector :math:`\mathcal{C}`, the Krylov preconditioner is
   :math:`P = (I + \mathcal{C})(L+C)^{-1}`, and a source-iteration step is
   Richardson iteration preconditioned by the same :math:`P`:
   :math:`\psi_{n+1} = \psi_n + P\,(q - (M - N)\,\psi_n)`. Source iteration
   realises it through the increment, ``ψ½ + C(ψ½ − ψ_n)``, so the two postures
   spell :math:`P` twice; the row holding them equal is the witness that keeps
   the spellings from drifting (the user's ruling, 2026-10-04: the full merge
   waits for a second corrector).

Every row here is ``python -O``-safe (asserts in a collected module, or
``np.testing``). Evidence: ``scratch/reference_architecture/p3/krylov200/``.
"""
from __future__ import annotations

import dataclasses
import functools
from typing import Any

import numpy as np
import pytest

import orpheus.sn.solver as sn_solver
from orpheus.geometry import CoordSystem
from orpheus.numerics.coupled_system import CoupledField
from orpheus.numerics.operator import IdentityOperator
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.splitting import Splitting, resolve_schedule
from tests.gates.sn.operators._full_space_states import (
    GEOMS,
    MESHES,
    random_composite,
    random_source_composite,
    random_state,
    system_a,
)
from tests.gates.sn.verification.analytical._aba_reference import aba_materials, aba_uniform_width_mesh
from tests.gates.sn.verification.analytical._case_slab_reference import case_slab_mesh

pytestmark = [pytest.mark.l1]

_HERE = "tests/gates/sn/solve/test_krylov_sweep_preconditioner.py"
_SWEEP_IDENTITY = (
    "tests/gates/sn/operators/test_sweep_inverse_identity.py::TestSweepInverseIdentity::"
    "test_forward_of_inverse_is_identity_on_a_random_composite"
)
#: Slab vacuum, slab reflective, 2-D Cartesian box, folded cylinder, GL sphere: every trace row live (the
#: cylinder's tangential slots aside). The 2-D row was added on qa's #200 review (2026-10-05); ``[M]`` its three
#: rows redden under the identity, boundary-drop and affine preconditioners
#: (``scratch/reference_architecture/p3/gates200b/arm4_*.log``).
_GEOMS = GEOMS


# ── the production preconditioner, read off the production driver ─────────


def _production_krylov(geom: str) -> tuple[Any, Any, Any]:
    """``(problem, M, Krylov)`` for ``geom``: the within-group splitting as the solver builds it, and the
    ``KrylovAcceleration`` the production driver ``_within_group_krylov`` builds on it (no corrector: the
    posture #200 changed)."""
    problem = MESHES[geom]()
    system = problem.system
    splitting = Splitting.from_schedule(system, resolve_schedule(problem, "jacobi"))
    krylov = sn_solver._within_group_krylov(
        splitting.implicit, *splitting.explicit, n_dof=1, max_iter=1, tol=1e-12,
    )
    return problem, splitting.implicit, krylov


def _boundary_only(q: Any) -> Any:
    """``q`` with its System-A bulk and every System-B block zeroed: a residual living on the trace alone."""
    a = system_a(q)
    a_b = dataclasses.replace(a, interior=a.interior * 0.0)
    if not isinstance(q, CoupledField):
        return a_b
    system_b: tuple[Any, ...] = tuple(q.systems[1:])
    return CoupledField(systems=(a_b, *(member * 0.0 for member in system_b)))


def _flat(x: Any) -> np.ndarray:
    return np.asarray(x.to_flat())


# ── 1. the preconditioner is the full-space inverse, linear on the full state ──


@pytest.mark.parametrize("geom", _GEOMS)
@pytest.mark.rests_on(*(f"{_SWEEP_IDENTITY}[{g}]" for g in _GEOMS))
def test_p200_1_the_preconditioner_is_the_full_space_sweep_inverse(geom: str) -> None:
    """ROUTE: the production preconditioner IS ``M.inverse()`` of the splitting's implicit part, bit for bit, on a
    random residual with every block populated. That inverse's identity ``M ∘ M⁻¹ = I`` on the full composite
    (outflow rows and ψ½ blocks included, ERR-071 / ERR-078) is the supporting gate's; this row says GMRES
    receives that object. Reddens on the identity (the pre-#200 production) and on any other preconditioner."""
    _problem, implicit, krylov = _production_krylov(geom)
    q = random_state(_problem, implicit, seed=29)
    np.testing.assert_array_equal(
        _flat(krylov.preconditioner.apply(q)), _flat(implicit.inverse().apply(q)),
        err_msg=f"{geom}: the Krylov preconditioner is not the sweep inverse of the splitting's implicit part",
    )


#: The linearity band, relative to the largest preconditioned value: the sweep is a linear recurrence, so the
#: defect is rounding. ``[M]`` 2026-10-04 (``krylov200/probes/p1_measure.log``, ``p2_blocks.log``): 1.2e-16 to
#: 3.4e-16 over the four meshes and the three splits; the band is about thirty times that.
_LINEARITY_RTOL = 1e-14


@pytest.mark.parametrize("geom", _GEOMS)
@pytest.mark.rests_on(*(f"{_HERE}::test_p200_1_the_preconditioner_is_the_full_space_sweep_inverse[{g}]" for g in _GEOMS))
def test_p200_1_the_preconditioner_is_linear_on_the_full_typed_state(geom: str) -> None:
    """M(q1 + q2) = M(q1) + M(q2) and M(αq) = αM(q) at every bulk, trace (inflow AND outflow) and ψ½ slot, for
    full random residuals and for residuals living on the boundary block alone (the block the issue's failure
    ignored). Reddens on an AFFINE preconditioner (a sweep carrying a seed between calls, the stateful-inverse
    class of ERR-050): additivity defects 6.9e-4 to 9.4e-4 (``[M]`` ``krylov200/battery/affine.log``). The
    additivity law is blind to a rank-deficient LINEAR map, the issue's own failure: that map reddens this row only
    through its non-vacuity check (the boundary-only residual is annihilated, scale 0), and the next row is the one
    that states it."""
    problem, implicit, krylov = _production_krylov(geom)
    precondition = krylov.preconditioner.apply
    full = [random_state(problem, implicit, seed=s) for s in (31, 37)]
    pairs = {"full": full, "boundary-only": [_boundary_only(q) for q in full]}
    for label, (q1, q2) in pairs.items():
        m1, m2, m12 = _flat(precondition(q1)), _flat(precondition(q2)), _flat(precondition(q1 + q2))
        scale = max(np.max(np.abs(m1)), np.max(np.abs(m2)))
        assert scale > 0.1, f"{geom}/{label}: the preconditioned residual vanished (scale {scale:.2e})"
        defect = np.max(np.abs(m12 - m1 - m2)) / scale
        assert defect <= _LINEARITY_RTOL, f"{geom}/{label}: additivity defect {defect:.2e} > {_LINEARITY_RTOL:.0e}"
        alpha = -2.5
        homogeneity = np.max(np.abs(_flat(precondition(q1 * alpha)) - alpha * m1)) / (abs(alpha) * scale)
        assert homogeneity <= _LINEARITY_RTOL, f"{geom}/{label}: homogeneity defect {homogeneity:.2e}"


@pytest.mark.parametrize("geom", _GEOMS)
def test_p200_1_a_boundary_only_residual_round_trips(geom: str) -> None:
    """THE ISSUE'S FAILURE: a residual on the trace alone (random on every inflow and outflow slot, zero bulk and
    ψ½) is preconditioned and mapped back by the forward, ``A(M q_b) = q_b``, on the bulk, the live trace rows and
    the ψ½ blocks. A preconditioner whose output ignores the boundary block maps ``q_b`` to zero and fails here;
    so does the identity (``A q_b ≠ q_b``). The tangential (μ_r = 0) trace slots of the folded cylinder are
    structural zero rows of the forward and are excluded, as in the supporting identity gate."""
    problem, implicit, krylov = _production_krylov(geom)
    q_b = _boundary_only(random_state(problem, implicit, seed=41))
    back = implicit.apply(krylov.preconditioner.apply(q_b))
    q_a, back_a = system_a(q_b), system_a(back)
    np.testing.assert_allclose(
        np.asarray(back_a.interior.values), 0.0, atol=1e-12,
        err_msg=f"{geom}: a boundary-only residual must map back to a zero bulk",
    )
    if isinstance(q_b, CoupledField):
        np.testing.assert_allclose(
            _flat(back.systems[1]), 0.0, atol=1e-12, err_msg=f"{geom}: and to a zero ψ½ state",
        )
    trace = problem.angular_trace
    n_live = 0
    for face in q_a.boundary.layout.faces:
        live = np.union1d(trace.inflow_indices_for_face(face), trace.outflow_indices_for_face(face))
        n_live += live.size
        np.testing.assert_allclose(
            np.asarray(back_a.boundary.face_view(face))[live], np.asarray(q_a.boundary.face_view(face))[live],
            rtol=1e-12, atol=1e-12,
            err_msg=f"{geom}/{face}: A(M q_b) must return the boundary residual on the live trace rows",
        )
    assert n_live > 0, f"{geom}: no live trace rows (vacuous)"


# ── 2 and 3. the solves: the fixed point, and the rate ──────────────────────

_INNER_TOL = 1e-10
_KEFF_TOL = 1e-12
#: Preconditioned against unpreconditioned: ten times the inner tolerance, for k and for the gauged scalar flux.
_FIXED_POINT_RTOL = 10.0 * _INNER_TOL


def _problem(body: str, n: int):
    """``(materials, mesh, quadrature)``: the Case two-region slab (reflective both faces) with ``n`` cells per
    region, or the A|B|A cylinder (reflective outer face, 2 groups) with ``n`` cells."""
    if body == "slab":
        mesh, materials, n_ordinates = case_slab_mesh(n)
        return materials, mesh, Quadrature.gauss_legendre(n_ordinates)
    return dict(aba_materials()), aba_uniform_width_mesh(CoordSystem.CYLINDRICAL, n), Quadrature.folded_product(n_mu=4, n_phi=8)


def _identity_driver(*args, **kwargs):
    """``_within_group_krylov`` with its preconditioner replaced by the identity: the unpreconditioned solve."""
    krylov = _honest_driver(*args, **kwargs)
    krylov.preconditioner = IdentityOperator()
    return krylov


_honest_driver = sn_solver._within_group_krylov


def _krylov_solve(body: str, n: int):
    """The Krylov eigenvalue solve of :func:`_problem` at the gates' tolerances, through whatever driver
    ``sn_solver._within_group_krylov`` currently names (the fixed-point row swaps in the identity)."""
    materials, mesh, quadrature = _problem(body, n)
    return sn_solver.solve_sn(materials, mesh, quadrature, inner_solver="krylov", max_outer=500, max_inner=500,
                              keff_tol=_KEFF_TOL, inner_tol=_INNER_TOL)


#: The production solves, cached across the fixed-point and the rate rows.
_production_solve = functools.cache(_krylov_solve)


def _gauged_flux(solution) -> np.ndarray:
    flux = np.asarray(solution.scalar_flux.values, dtype=float)
    return flux / flux.sum()


@pytest.mark.parametrize("body", ["slab", "cylinder"])
@pytest.mark.rests_on(f"{_HERE}::test_p200_1_the_preconditioner_is_the_full_space_sweep_inverse[slab_reflective]",
                      f"{_HERE}::test_p200_1_the_preconditioner_is_the_full_space_sweep_inverse[cyl_folded]")
def test_p200_2_the_fixed_point_is_the_unpreconditioned_one(body: str, monkeypatch: pytest.MonkeyPatch) -> None:
    """The preconditioned and the unpreconditioned Krylov solves (the latter: the identity put back in the
    production driver) agree on k and on the gauged scalar flux to :data:`_FIXED_POINT_RTOL`, at 10 cells per
    region (slab) or 10 cells (cylinder), reflective faces. ``[M]`` 2026-10-04 (``krylov200/probes/
    p4_fixed_point.log``): |Δk|/k 8.5e-12 and 7.2e-12, gauged flux 9.2e-11 and 4.2e-11. Reddens on a
    preconditioner that ignores the boundary block: GMRES then reports a running residual of 0 to 7e-34 (the
    issue's machine-zero wrong answer) and production's convergence-claim check refuses it, honest residual 0.23
    and 0.52 (``krylov200/battery/boundary_drop.log``). Against the identity itself the row is green by
    construction: the rate row is what tells the two apart."""
    preconditioned = _production_solve(body, 10)
    monkeypatch.setattr(sn_solver, "_within_group_krylov", _identity_driver)
    unpreconditioned = _krylov_solve(body, 10)
    k_p, k_u = float(preconditioned.outcome.keff), float(unpreconditioned.outcome.keff)
    assert abs(k_p - k_u) <= _FIXED_POINT_RTOL * k_u, f"{body}: k {k_p!r} against {k_u!r}"
    phi_p, phi_u = _gauged_flux(preconditioned), _gauged_flux(unpreconditioned)
    gap = np.max(np.abs(phi_p - phi_u)) / np.max(np.abs(phi_u))
    assert gap <= _FIXED_POINT_RTOL, f"{body}: gauged scalar flux differs by {gap:.2e}"


#: The ladder of the rate row, and the growth it allows between its ends.
_RATE_LADDER = (10, 20, 40)
#: ``[M]`` 2026-10-04 (``krylov200/probes/p1_measure.log``), the largest inner iteration count of a solve at
#: 10, 20, 40 cells: preconditioned, slab 15, 15, 15 and cylinder 23, 24, 22; the identity, slab 176, 336, 656
#: and cylinder 238, 490, 950 (a growth of 3.7 and 4.0 over the ladder, as the unknowns).
_RATE_GROWTH = 1.25


@pytest.mark.parametrize("body", ["slab", "cylinder"])
@pytest.mark.rests_on(*(f"{_HERE}::test_p200_2_the_fixed_point_is_the_unpreconditioned_one[{b}]" for b in ("slab", "cylinder")))
def test_p200_3_inner_iterations_stay_flat_under_refinement(body: str) -> None:
    """The largest inner (GMRES) iteration count of a solve, at the finest mesh of :data:`_RATE_LADDER`, is at
    most :data:`_RATE_GROWTH` times that at the coarsest: the sweep preconditions the streaming, so the count
    follows the scattering, not the mesh. Reddens on the identity (growth 3.7 slab, 4.0 cylinder)."""
    largest = {}
    for n in _RATE_LADDER:
        record = _production_solve(body, n).record
        counts = [child.n_iterations for child in record.children]
        assert counts, f"{body} n={n}: the record carries no inner solves"
        largest[n] = max(counts)
    coarse, fine = largest[_RATE_LADDER[0]], largest[_RATE_LADDER[-1]]
    assert fine <= _RATE_GROWTH * coarse, (
        f"{body}: the largest inner iteration count grows {fine / coarse:.2f}x over cells {_RATE_LADDER} "
        f"({largest}), above {_RATE_GROWTH}"
    )


# ── 4. under DSA, source iteration is Richardson on the Krylov preconditioner ──

#: ``[M]`` 2026-10-04 (``scratch/reference_architecture/p3/operator_contract/probe_witness.py``): 6.2e-16 and
#: 4.3e-16 relative on the two slabs; with the corrector dropped from ``P``, 0.42 and 0.39.
_RICHARDSON_RTOL = 1e-13


@pytest.mark.parametrize("geom", ["slab_vacuum", "slab_reflective"])
@pytest.mark.rests_on(f"{_SWEEP_IDENTITY}[slab_vacuum]", f"{_SWEEP_IDENTITY}[slab_reflective]")
def test_p200_4_a_dsa_source_iteration_step_is_richardson_on_the_krylov_preconditioner(geom: str) -> None:
    """One production source-iteration step with the consistent-DSA corrector, from a random flux ``ψ_n`` and a
    random source ``q`` (every block populated), equals ``ψ_n + P (q − (M − N) ψ_n)`` with ``P`` the preconditioner
    the production Krylov builder hands GMRES under the same corrector. The identity behind it is
    ``M⁻¹(q + Nψ_n) − ψ_n = M⁻¹ r_n``, so the row rests on the exact inverse (the supporting gate). Reddens when
    either posture drops or alters the corrector: ``P`` without it misses by 0.42 and 0.39. Slab only: consistent
    DSA is admitted on the 1-D Cartesian slab alone (``DSALowOrderSystem.from_problem``)."""
    from orpheus.sn.acceleration.dsa import DSACorrection

    problem = MESHES[geom]()
    splitting = Splitting.from_schedule(problem.system, resolve_schedule(problem, "jacobi"))
    implicit, gains = splitting.implicit, splitting.explicit
    corrector = DSACorrection.from_problem(problem)
    q, psi_n = random_source_composite(problem, seed=7), random_composite(problem, seed=11)

    source_iteration, *_ = sn_solver._within_group_si(splitting, problem, max_iter=1, tol=0.0, corrector=corrector)
    step, record = source_iteration.solve(q, initial_guess=psi_n)
    assert record.iterations_run == 1, f"{geom}: expected one corrected step, ran {record.iterations_run}"

    preconditioner = sn_solver._within_group_krylov(
        implicit, *gains, n_dof=1, max_iter=1, tol=1e-12, corrector=corrector,
    ).preconditioner
    residual = q - implicit.apply(psi_n)
    for gain in gains:
        residual = residual + gain.apply(psi_n)
    richardson = psi_n + preconditioner.apply(residual)

    gap = np.max(np.abs(_flat(step) - _flat(richardson))) / np.max(np.abs(_flat(step)))
    assert gap <= _RICHARDSON_RTOL, (
        f"{geom}: the DSA source-iteration step differs from Richardson on the Krylov preconditioner by {gap:.2e}"
    )


# ── 5. the inner record reads the loop's own test under the sweep ────────────

#: qa's thin slab (#200 review, item 3), shrunk to the smallest mesh that still shows the defect: two 10 cm
#: regions of a 2-group mixture with upscatter, scaled to σ_t = 0.02 and 0.04 /cm (the second region at half
#: density), vacuum faces, 2 cells per region, Gauss-Legendre S8. The sweep then amplifies the source,
#: ‖M b‖/‖b‖ ≈ 26 (``[M]`` 2026-10-05, ``scratch/reference_architecture/p3/gates200b/probe_l1_size.log``: of
#: five sizes the one where the raw reading misreads the most inner solves, 4 of 8, in 0.03 s).
_THIN_SCALE = 0.02
#: The activation the row needs: under the raw reading the record's residual is ‖M b‖/‖b‖ times the loop's, so
#: an amplification below a few would leave the raw reading inside the tolerance and the row blind.
_THIN_MIN_AMPLIFICATION = 10.0


def _thin_mixture(scale: float):
    from orpheus.derivations.common.xs_library import make_mixture

    sig_t = np.array([1.0, 2.0]) * scale
    sig_s = np.array([[0.5, 0.4], [0.05, 1.5]]) * scale
    sig_f = np.array([0.02, 0.15]) * scale
    return make_mixture(
        sig_t=sig_t, sig_c=sig_t - sig_s.sum(axis=1) - sig_f, sig_f=sig_f,
        nu=np.array([2.5, 2.5]), chi=np.array([1.0, 0.0]), sig_s=sig_s,
    )


@pytest.mark.rests_on(
    "tests/gates/numerics/test_iteration.py::test_krylov_record_verdict_agrees_with_scipy_acceptance[accepted_s16]",
    "tests/gates/numerics/test_iteration.py::test_krylov_record_is_invariant_under_preconditioner_scaling",
    f"{_HERE}::test_p200_1_the_preconditioner_is_the_full_space_sweep_inverse[slab_vacuum]",
)
def test_p200_5_no_inner_record_reads_unconverged_when_scipy_accepted(monkeypatch: pytest.MonkeyPatch) -> None:
    """On a thin slab, where the sweep amplifies the residual about 26 times, every inner (GMRES) record reads
    converged when scipy accepted the solve (``info = 0``), no "inner(gmres) hit" budget warning fires, and the
    solve reads fully converged. The record's criterion is ``‖M r‖/‖M b‖``, the quantity scipy's inner loop stops
    on (the user's ruling, 2026-10-04). Reddens on the raw callback reading ``‖M r‖/‖b‖``: ``[M]`` 4 of 8 inner
    records then read not-converged, the budget warning fires, ``fully_converged`` is False
    (``scratch/reference_architecture/p3/gates200b/``). The opposite error, a record claiming a convergence the
    true residual does not support, is the end-of-solve claim check's; on 72 slab configurations it never arose
    beyond that check's safety factor (README, item 3)."""
    import warnings

    import orpheus.numerics.iteration as iteration_mod
    from orpheus.geometry import BC, StructuredGeometry
    from orpheus.mesh import CellsByCount, Mesher

    calls: list[tuple[int, float]] = []
    honest = iteration_mod.spla.gmres

    def spy(A: Any, b: np.ndarray, *args: Any, M: Any = None, **kwargs: Any) -> tuple[np.ndarray, int]:
        x, info = honest(A, b, *args, M=M, **kwargs)
        calls.append((int(info), float(np.linalg.norm(M.matvec(b)) / np.linalg.norm(b))))
        return x, info

    monkeypatch.setattr(iteration_mod.spla, "gmres", spy)
    geometry = StructuredGeometry.slab((0.0, 10.0, 20.0), (0, 1), left=BC.vacuum, right=BC.vacuum)
    mesh = Mesher(geometry).partition(CellsByCount.uniform_width(2)).mesh
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        solution = sn_solver.solve_sn(
            {0: _thin_mixture(_THIN_SCALE), 1: _thin_mixture(_THIN_SCALE / 2)}, mesh, Quadrature.gauss_legendre(8),
            inner_solver="krylov", keff_tol=1e-10, inner_tol=1e-8, max_outer=300,
        )
    inner = [child for child in solution.record.children if child.label == "inner(gmres)"]
    assert inner and len(inner) == len(calls), f"{len(inner)} inner(gmres) records against {len(calls)} gmres calls"
    amplification = min(ratio for _info, ratio in calls)
    assert amplification > _THIN_MIN_AMPLIFICATION, (
        f"the sweep no longer amplifies the residual (‖M b‖/‖b‖ ≥ {amplification:.1f}): the row is blind"
    )
    misread = [i for i, ((info, _r), child) in enumerate(zip(calls, inner)) if info == 0 and not child.converged]
    assert not misread, (
        f"inner solves {misread} of {len(inner)}: scipy accepted (info = 0) but the record reads not converged "
        f"(last pr_residual {[inner[i].criteria[0].last for i in misread]}, tol {inner[0].criteria[0].tolerance:.0e})"
    )
    budget = [str(w.message) for w in caught if "inner(gmres) hit" in str(w.message)]
    assert not budget, f"a false inner budget warning fired: {budget[0][:160]}"
    assert solution.record.fully_converged, "the solve does not read fully converged"
