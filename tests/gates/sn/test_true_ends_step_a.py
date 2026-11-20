r"""#707 Step A, the S_N rows: the data side V* = A ⊗ (E ⊗ X)*, the operator ends, the refusals, the adjoint, the gauge.

Spec: ``scratch/true_ends/stepA/gate_spec.md`` (row ids in each docstring). Written BEFORE the carve: every row
marked RED below is red today for the reason Step A fixes, and the rows marked GREEN are the positive legs that must
stay green (the composition that is legal, the state side that does not move).

The problems are the chart-law gate's (``tests/gates/cross_method/_source_chart_cases.py``): a non-uniform
fuel | moderator slab (DD, P1), a two-region sphere (DD), a non-uniform 2-D box (DD, P1), and their LD siblings; all
2-group, heterogeneous, with vacuum and mirror faces, so no row sits on a fixture where the measure is 1.

The side of a space is read through its metric, never through a class or a name (R8, the ruling of 2026-11-12): a
data-side bulk space carries the angular weight ``w_n`` on its angular factor (primal: angle is a function variable)
and the inverse spatial form ``1 / (V θ^{|o|})`` on energy ⊗ space ⊗ in-cell moments (dual). Whatever object the
carve spells it with (a ``TensorProductSpace`` of the angular factor and the dual of the rest, or a new type), that is
the metric it must carry, so the row is spelling-independent.
"""

from __future__ import annotations

import warnings
from functools import cache

import numpy as np
import pytest

from orpheus.numerics.coordinate_space import NotInCoordinateSpace
from orpheus.numerics.operator import IncompatibleOperatorComposition, MissingAdjoint
from tests.gates._true_ends import H_REFUSAL_FRAGMENT, pairing_of
from tests.gates.cross_method import _source_chart_cases as cc

pytestmark = [pytest.mark.foundation, pytest.mark.filterwarnings("ignore::Warning")]


_PROBLEMS = ("slab_dd_p1", "sphere_dd", "2d_dd_p1", "slab_ld_p1", "2d_ld")
_FLUX_TO_SOURCE = ("StreamingCollisionOperator", "ScatteringOperator", "N2NOperator", "SNBoundaryOperator",
                   "FissionOperator")


@cache
def _problem(name: str):
    return cc.problems()[name]()


def _atom(name: str, atom: str):
    found = cc.atoms(_problem(name))
    if atom not in found:
        pytest.skip(f"{atom} is not on {name}'s pencil")  # FissionOperator on the sphere (System B carries it)
    return found[atom]


def _psi(name: str):
    return cc.state(_problem(name), np.random.default_rng(cc.SEED))


def _cases() -> list[tuple[str, str]]:
    return [(p, a) for p in _PROBLEMS for a in _FLUX_TO_SOURCE]


def _expected_data_metric(problem, shape: tuple[int, ...]) -> np.ndarray:
    """``w_n ⊗ 1_g ⊗ 1/(V θ^{|o|})``, from raw quadrature weights and mesh volumes (R8)."""
    w = np.asarray(problem.quad.weights, dtype=float).reshape((-1,) + (1,) * (len(shape) - 1))
    return w / cc.sn_bulk_measure(problem, shape)


def _expected_state_metric(problem, shape: tuple[int, ...]) -> np.ndarray:
    w = np.asarray(problem.quad.weights, dtype=float).reshape((-1,) + (1,) * (len(shape) - 1))
    return w * cc.sn_bulk_measure(problem, shape)


# ── SA-1: the side identity of a source (RED) and of a flux (GREEN control) ───────────────────────────────────────


@pytest.mark.parametrize(("name", "atom"), _cases())
@pytest.mark.xfail(strict=True, reason="#707 Step A: red by design until the true-ends carve (ontology chart W1 to W4); the flip-edit removes this marker")
def test_sa1_codomain_is_the_data_side(name: str, atom: str) -> None:
    """SA-1 RED. A flux → source atom's codomain is not its domain, and its bulk block carries the data-side metric.

    First red today: ``codomain == domain`` (census verdict c: one cached ``FullFieldSpace`` at both ends).
    """
    op = _atom(name, atom)
    assert op.codomain != op.domain, f"{atom} on {name}: codomain == domain (the role lives in the class, not the space)"
    interior = op.codomain.interior_space
    got = np.asarray(interior.metric_diagonal(), dtype=float)
    want = _expected_data_metric(_problem(name), got.shape)
    np.testing.assert_allclose(got, want, rtol=4 * np.finfo(float).eps, atol=0,
                               err_msg=f"{atom} on {name}: the codomain's bulk metric is not w_n ⊗ 1/(V θ^|o|)")


@pytest.mark.parametrize("name", _PROBLEMS)
def test_sa1_state_side_is_unchanged(name: str) -> None:
    """SA-1 GREEN control: the state side V keeps ``w_n ⊗ V θ^{|o|}`` (R6: no flux number moves)."""
    interior = _problem(name).full_field_space.interior_space
    got = np.asarray(interior.metric_diagonal(), dtype=float)
    np.testing.assert_allclose(got, _expected_state_metric(_problem(name), got.shape), rtol=4 * np.finfo(float).eps)


@pytest.mark.parametrize(("name", "atom"), _cases())
@pytest.mark.xfail(strict=True, reason="#707 Step A: red by design until the true-ends carve (ontology chart W1 to W4); the flip-edit removes this marker")
def test_sa2_the_member_lives_in_the_codomain(name: str, atom: str) -> None:
    """SA-2 RED. The atom's output is a member of the data side (its space is the codomain, not the flux space).

    First red today: the output's interior space is the flux's interior space.
    """
    op = _atom(name, atom)
    try:
        out = op.apply(_psi(name))
    except ValueError as exc:  # FissionOperator on LD raises today (spec CL-LDF), outside Step A
        pytest.xfail(f"pre-existing LD fission defect: {exc}")
    assert out.interior.space != _problem(name).full_field_space.interior_space, (
        f"{atom} on {name}: the source member lives in the flux space")
    assert out.interior.space == op.codomain.interior_space


# ── SA-3: the inverse runs V* → V, and refuses a flux (RED), admits a source (GREEN) ─────────────────────────────


@pytest.mark.parametrize("name", _PROBLEMS)
@pytest.mark.xfail(strict=True, reason="#707 Step A: red by design until the true-ends carve (ontology chart W1 to W4); the flip-edit removes this marker")
def test_sa3_the_sweep_ends_are_reversed(name: str) -> None:
    """SA-3 RED. ``(L + C)⁻¹ : V* → V``: its domain is the loss's codomain and differs from its codomain."""
    lc = _atom(name, "StreamingCollisionOperator")
    sweep = lc.inverse()
    assert sweep.domain == lc.codomain and sweep.codomain == lc.domain
    assert sweep.domain != sweep.codomain, f"{name}: the sweep declares one space at both ends"


@pytest.mark.parametrize("name", _PROBLEMS)
@pytest.mark.xfail(strict=True, reason="#707 Step A: red by design until the true-ends carve (ontology chart W1 to W4); the flip-edit removes this marker")
def test_sa3_the_sweep_refuses_a_flux(name: str) -> None:
    """SA-3 RED. A flux member handed to the sweep is refused at the membership leg.

    First red today: accepted (census: 19 448 production calls through ``FlattenedOperator``). The refusal type is the
    algebra's membership refusal; its message fragment is pinned when the carve writes it (spec NEEDS).
    """
    sweep = _atom(name, "StreamingCollisionOperator").inverse()
    with pytest.raises((NotInCoordinateSpace, IncompatibleOperatorComposition)):
        sweep.apply(_psi(name))


@pytest.mark.parametrize("name", ("slab_dd_p1", "2d_dd_p1"))
def test_sa3_the_sweep_inverts_the_loss(name: str) -> None:
    """SA-3 GREEN. ``(L + C)⁻¹ (L + C) ψ = ψ`` on the bulk: a source member is admitted and the inverse is exact."""
    lc = _atom(name, "StreamingCollisionOperator")
    psi = _psi(name)
    back = lc.inverse().apply(lc.apply(psi))
    scale = float(np.max(np.abs(psi.interior.values)))
    err = float(np.max(np.abs(back.interior.values - psi.interior.values)))
    assert err <= 1e-10 * scale, f"{name}: sweep ∘ loss != id on the bulk: {err:.3e}"


# ── SA-4: composition reads the sides ─────────────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("name", _PROBLEMS)
def test_sa4_loss_terms_compose(name: str) -> None:
    """SA-4 GREEN. ``L + C − S`` composes: the summands share both ends, and the sum keeps them."""
    lc, s = _atom(name, "StreamingCollisionOperator"), _atom(name, "ScatteringOperator")
    a = lc - s
    assert a.domain == lc.domain and a.codomain == lc.codomain


@pytest.mark.parametrize("name", _PROBLEMS)
@pytest.mark.xfail(strict=True, reason="#707 Step A: red by design until the true-ends carve (ontology chart W1 to W4); the flip-edit removes this marker")
def test_sa4_a_loss_does_not_compose_with_a_loss(name: str) -> None:
    """SA-4 RED. ``(L + C) ∘ (L + C)`` is ill-typed (V → V* then V → V*) and refuses.

    First red today: it composes, because both ends are one space. (``L + I`` would be the other spelling of the
    side mismatch, but ``IdentityOperator`` declares no ends and composes with anything: spec SA-4b, no first red.)
    """
    lc = _atom(name, "StreamingCollisionOperator")
    with pytest.raises(IncompatibleOperatorComposition):
        _ = lc @ lc


@pytest.mark.parametrize("name", _PROBLEMS)
def test_sa4_sweep_after_loss_composes(name: str) -> None:
    """SA-4 GREEN. ``(L + C)⁻¹ ∘ S`` (V → V* → V) composes: the source iteration's operator."""
    lc, s = _atom(name, "StreamingCollisionOperator"), _atom(name, "ScatteringOperator")
    k = lc.inverse() @ s
    assert k.domain == s.domain and k.codomain == lc.domain


# ── SA-5 / SA-6: the Hilbert adjoint refuses a V → V* operator; the pairing adjoint is the transport adjoint ─────


@pytest.mark.parametrize(("name", "atom"), _cases())
@pytest.mark.xfail(strict=True, reason="#707 Step A: red by design until the true-ends carve (ontology chart W1 to W4); the flip-edit removes this marker")
def test_sa5_hilbert_adjoint_refuses_across_sides(name: str, atom: str) -> None:
    """SA-5 RED. ``.H`` (♯ ∘ dual ∘ ♭ with both ends' metrics) is not the transport adjoint of V → V*, and refuses,
    naming its replacement, ``DualPairing.adjoint`` (Pattern 4 corollary (3): the replacement refuses the old idiom).

    First red today: ``.H`` constructs (15 of 18 production ``.H`` sites take it, spec §2.4).
    """
    op = _atom(name, atom)
    with pytest.raises(MissingAdjoint, match=H_REFUSAL_FRAGMENT):
        op.H


def _dense(op, name: str, *, transpose: bool = False) -> np.ndarray:
    """``op`` (or its transpose) applied to every unit member of the full field, as a dense matrix on (bulk ⊕ trace)."""
    problem = _problem(name)
    psi = cc.state(problem, np.random.default_rng(0))
    nb, nt = psi.interior.values.size, psi.boundary.values.size
    cols = []
    for j in range(nb + nt):
        e = np.zeros(nb + nt)
        e[j] = 1.0
        x = cc.state_from(problem, e[:nb].reshape(psi.interior.values.shape), e[nb:])
        y = op.apply_transpose(x) if transpose else op.apply(x)
        cols.append(np.concatenate([np.ravel(y.interior.values), np.ravel(y.boundary.values)]))
    return np.stack(cols, axis=1)


def _pairing_weights(name: str) -> np.ndarray:
    """The V* × V pairing on (bulk ⊕ trace): ``w_n`` on the bulk (angular metric ⊗ canonical spatial), and on the trace
    today's trace metric ``|Ω·n| w_n`` (the trace keeps its chart in Step A; Step B re-poses it)."""
    problem = _problem(name)
    psi = cc.state(problem, np.random.default_rng(0))
    shape = psi.interior.values.shape
    w = np.broadcast_to(np.asarray(problem.quad.weights, float).reshape((-1,) + (1,) * (len(shape) - 1)), shape)
    trace = np.asarray(problem.full_field_space.trace_space.metric_diagonal(), dtype=float).ravel()
    return np.concatenate([w.ravel(), trace])


@pytest.mark.parametrize("atom", ("StreamingCollisionOperator", "ScatteringOperator", "FissionOperator"))
def test_sa6r_hilbert_adjoint_today_record(atom: str) -> None:
    """SA-6r RECORD, GREEN today and designed to redden at the carve (``.H`` refuses, SA-5). Runs SA-6's apparatus
    (``_dense``, the live-slot mask) against today's ``.H = G⁻¹AᵀG`` with the FULL Hilbert weights (bulk
    ``w_n V``, trace ``|Ω·n| w_n``), so the strict-xfail row's machinery has executed before its verb exists."""
    name = "slab_dd_p1"
    op = _atom(name, atom)
    A = _dense(op, name)
    problem = _problem(name)
    psi = cc.state(problem, np.random.default_rng(0))
    shape = psi.interior.values.shape
    G = np.concatenate([_expected_state_metric(problem, shape).ravel(),
                        np.asarray(problem.full_field_space.trace_space.metric_diagonal(), float).ravel()])
    live = G != 0.0
    want = (A.T * G[None, :]) / np.where(live, G, 1.0)[:, None]
    got = _dense(op.H, name)
    scale = float(np.max(np.abs(want[np.ix_(live, live)])))
    err = float(np.max(np.abs(got[np.ix_(live, live)] - want[np.ix_(live, live)])))
    assert err <= 256 * np.finfo(float).eps * scale, f"{atom}: today's .H != G⁻¹AᵀG: {err / scale:.3e}"


@pytest.mark.rests_on("tests/gates/sn/test_true_ends_step_a.py::test_sa6r_hilbert_adjoint_today_record", "tests/gates/sn/test_true_ends_step_a.py::test_sa5_hilbert_adjoint_refuses_across_sides")
@pytest.mark.xfail(strict=True, reason="#707 Step A: DualPairing (which owns the transport adjoint) is not built")
@pytest.mark.parametrize("atom", ("StreamingCollisionOperator", "ScatteringOperator", "FissionOperator"))
def test_sa6_pairing_adjoint_is_angular_legs_only(atom: str) -> None:
    """SA-6 (strict xfail until ``DualPairing`` lands). ``pairing.adjoint(A) = P⁻¹ Aᵀ P`` with ``P`` the pairing weights (``w_n`` on the bulk:
    the ANGULAR legs only; the spatial legs retired with R6), on the live slots; and ``(A†)† = A`` bit for bit.

    The flip-edit removes the marker; the body's values are the carve's (Mode 8(10) check: it is not a copy of SA-5).
    Fixture: the non-uniform slab, P1 (scattering couples ordinates, so ``P`` does not commute with ``Aᵀ``).
    """
    name = "slab_dd_p1"
    op = _atom(name, atom)
    pairing = pairing_of(op.domain, op.codomain)
    adj = pairing.adjoint(op)
    A = _dense(op, name)
    P = _pairing_weights(name)
    live = P != 0.0
    want = (A.T * P[None, :]) / np.where(live, P, 1.0)[:, None]
    got = _dense(adj, name)
    scale = float(np.max(np.abs(want[np.ix_(live, live)])))
    err = float(np.max(np.abs(got[np.ix_(live, live)] - want[np.ix_(live, live)])))
    assert err <= 256 * np.finfo(float).eps * scale, f"{atom}: A† != P⁻¹AᵀP: {err / scale:.3e}"
    again = pairing.adjoint(adj)
    assert np.array_equal(_dense(again, name), A), f"{atom}: (A†)† != A"


# ── SA-7 / SA-8: the adjoint eigen gauge and the Rayleigh estimator are pairings ──────────────────────────────────


def _transport_pairing_with_unit_fission(name: str, psi_star: np.ndarray) -> float:
    r"""``⟨F†ψ*, 1⟩ = ⟨ψ*, F 1⟩ = Σ_n w_n Σ_{g,i} ψ*_{n,g,i} · V_i χ_g(i) Σ_{g'} νΣf_{g'}(i)``, from raw data.

    ``F 1`` per steradian: the isotropic emission of the unit flux (φ = Σ w_n = W) is ``χ νΣf φ / W = χ νΣf``.
    """
    problem = _problem(name)
    mats = cc._mixtures()
    mat_of_cell = np.asarray(cc._slab_mesh().mat_ids, dtype=int)  # per cell, read from the mesh
    chi = np.stack([np.asarray(mats[m].chi, float) for m in mat_of_cell], axis=1)  # (ng, nx)
    nusf = np.stack([np.asarray(mats[m].SigP, float) for m in mat_of_cell], axis=1).sum(axis=0)  # (nx,)
    V = np.asarray(problem.volumes, float)
    f_unit = chi * nusf[None, :] * V[None, :]  # (ng, nx)
    w = np.asarray(problem.quad.weights, float)
    return float(np.einsum("n,ngi,gi->", w, psi_star, f_unit))


@cache
def _adjoint_eigen():
    from orpheus.sn import solve_sn_adjoint

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mesh = cc._slab_mesh()
        return solve_sn_adjoint(cc._mixtures(), mesh, _problem("slab_dd_p1").quad, scattering_order=1,
                                keff_tol=1e-11, flux_tol=1e-10, max_inner=2000, inner_tol=1e-12,
                                boundary_conditions=cc._conditions(mesh))


@pytest.mark.rests_on("tests/gates/cross_method/test_true_ends_acceptance.py::test_state_is_unmoved", "tests/gates/cross_method/test_source_chart_law.py::test_atom_output_is_integrated")
@pytest.mark.xfail(strict=True, reason="#707 Step A: red by design until the true-ends carve (ontology chart W1 to W4); the flip-edit removes this marker")
def test_sa7_adjoint_gauge_is_the_pairing() -> None:
    """SA-7 RED (census C8, ``iteration.py:1590``; ruled 2026-11-12). The adjoint importance is gauged by the
    transport pairing, ``⟨F†ψ*, 1⟩ = 1``, not by the coordinate sum of ``F†ψ*``; the importance moves by one declared
    global factor (the acceptance gauge event).

    Oracle (hand-built, from raw χ, νΣf, V, w_n and the returned ψ*): ``⟨ψ*, F 1⟩``. Today it reads 0.05 while the
    coordinate gauge reads 1.
    """
    sol = _adjoint_eigen()
    psi_star = np.asarray(sol.angular_flux.interior.values, dtype=float)
    value = _transport_pairing_with_unit_fission("slab_dd_p1", psi_star)
    assert abs(value - 1.0) <= 1e-9, f"⟨F†ψ*, 1⟩ = {value!r}, the gauge value is 1"


@pytest.mark.xfail(strict=True, reason="#707 Step A: red by design until the true-ends carve (ontology chart W1 to W4); the flip-edit removes this marker")
def test_sa7_the_pairing_evaluates_the_gauge() -> None:
    """SA-7p RED. The production leg through ``DualPairing``: ``pairing(ψ*, F·1)`` equals the hand-built oracle
    ``⟨ψ*, F 1⟩`` (to 1e-13 relative). First red today: no ``DualPairing``."""
    problem = _problem("slab_dd_p1")
    F = _atom("slab_dd_p1", "FissionOperator")
    pairing = pairing_of(F.domain, F.codomain)
    sol = _adjoint_eigen()
    psi_star = cc.state_from(problem, np.asarray(sol.angular_flux.interior.values, float),
                             np.asarray(sol.angular_flux.boundary.values, float))
    unit = cc.state_from(problem, np.ones(psi_star.interior.values.shape), np.zeros(psi_star.boundary.values.shape))
    want = _transport_pairing_with_unit_fission("slab_dd_p1", np.asarray(psi_star.interior.values, float))
    got = float(pairing(psi_star, F.apply(unit)))
    assert abs(got - want) <= 1e-13 * abs(want), (got, want)


def test_sa7_adjoint_k_is_forward_k() -> None:
    """SA-7 GREEN (REFERENCE across the two daggered solves): ``k† = k`` to 10 × keff_tol on the P1 slab."""
    from orpheus.sn import solve_sn

    mesh = cc._slab_mesh()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fwd = solve_sn(cc._mixtures(), mesh, _problem("slab_dd_p1").quad, scattering_order=1, keff_tol=1e-11,
                       flux_tol=1e-10, max_inner=2000, inner_tol=1e-12, boundary_conditions=cc._conditions(mesh))
    assert abs(float(fwd.outcome.keff) - float(_adjoint_eigen().outcome.keff)) <= 1e-10


def _anchored_slab():
    name = "slab_dd_p1"
    a = cc_anchor()
    problem = _problem(name)
    psi = cc.state_from(problem, a[f"{name}/state/bulk"], a[f"{name}/state/trace"])
    return name, a, problem, psi


def _loss(a: dict, name: str, block: str) -> np.ndarray:
    """``A = (L+C) − S − N2n − B`` from the anchored density outputs (the pencil's signs)."""
    g = lambda atom: a[f"{name}/{atom}/{block}"]  # noqa: E731
    return g("StreamingCollisionOperator") - g("ScatteringOperator") - g("N2NOperator") - g("SNBoundaryOperator")


@pytest.mark.rests_on("tests/gates/cross_method/test_source_chart_law.py::test_atom_output_is_integrated")
@pytest.mark.xfail(strict=True, reason="#707 Step A: red by design until the true-ends carve (ontology chart W1 to W4); the flip-edit removes this marker")
def test_sa8_rayleigh_bulk_part_is_the_pairing() -> None:
    """SA-8 RED (census C9, ``posing.py:163``; ruled 2026-11-12: the Rayleigh estimator pairs over the FULL space,
    bulk ⊕ trace). Step A gates the BULK part: ``DualPairing`` evaluated on (1, y) with y's trace zeroed, for
    ``y = Fψ`` and ``y = Aψ`` at a non-eigen state, equals ``Σ_n w_n Σ_i m (y)_density`` from the anchor.

    Oracle: the anchored density outputs × m × w_n (hand-built). Production leg: ``DualPairing`` on today's
    operators. First red today: no ``DualPairing``. (Today's ``rayleigh(w=1)`` reads 0.012553693834414816, the
    coordinate sum over bulk and trace exactly; the bulk pairing ratio is 0.02111190814237901.)
    """
    name, a, problem, psi = _anchored_slab()
    shape = a[f"{name}/state/bulk"].shape
    m = cc.sn_bulk_measure(problem, shape)
    w = np.asarray(problem.quad.weights, float).reshape((-1,) + (1,) * (len(shape) - 1))
    lc = _atom(name, "StreamingCollisionOperator")
    pairing = pairing_of(lc.domain, lc.codomain)
    unit = cc.state_from(problem, np.ones(shape), np.zeros(psi.boundary.values.shape))
    loss = lc - _atom(name, "ScatteringOperator") - _atom(name, "N2NOperator") - _atom(name, "SNBoundaryOperator")
    for label, op, density in (("F", _atom(name, "FissionOperator"), a[f"{name}/FissionOperator/bulk"]),
                               ("A", loss, _loss(a, name, "bulk"))):
        y = op.apply(psi)
        y.boundary.values[...] = 0.0
        want = float(np.sum(w * m * density))
        got = float(pairing(unit, y))
        assert abs(got - want) <= 1e-12 * abs(want), f"⟨1, ({label}ψ)_bulk⟩: {got!r} != {want!r}"


@pytest.mark.xfail(strict=True, reason="#707 Step B: the trace term of the Rayleigh pairing needs the current-posed "
                                       "trace (R4, R9); in Step A the trace is still the density chart")
@pytest.mark.rests_on(f"{__name__.replace('.', '/')}.py::test_sa8_rayleigh_bulk_part_is_the_pairing")
def test_sa8b_rayleigh_pairs_the_full_space() -> None:
    """SA-8b (strict xfail, pending Step B). ``rayleigh(ψ, w=1) = ⟨1, Fψ⟩ / ⟨1, Aψ⟩`` over bulk ⊕ trace, the trace
    weighted by its current measure ``|Ω·n| w_n A_f`` (R9).

    Oracle: the anchor; bulk ``w_n m``, trace ``|Ω·n| w_n A_f`` applied to today's density trace values. Blindness
    declared: on the slab ``A_f = 1`` on both faces, so the face area is not seen here (a curvilinear or 2-D row owes it).
    """
    from orpheus.numerics.coupled_system import CoupledField

    name, a, problem, psi = _anchored_slab()
    shape = a[f"{name}/state/bulk"].shape
    m = cc.sn_bulk_measure(problem, shape)
    w = np.asarray(problem.quad.weights, float).reshape((-1,) + (1,) * (len(shape) - 1))
    current = np.asarray(problem.full_field_space.trace_space.metric_diagonal(), float).ravel()  # |Ω·n| w_n, A_f = 1
    num = np.sum(w * m * a[f"{name}/FissionOperator/bulk"]) + np.sum(current * a[f"{name}/FissionOperator/trace"])
    den = np.sum(w * m * _loss(a, name, "bulk")) + np.sum(current * _loss(a, name, "trace"))
    got = problem.eigen_posing.rayleigh(CoupledField(systems=(psi,)), w=1.0)
    assert abs(got - num / den) <= 1e-12 * abs(num / den), (got, num / den)


# ── LD-F: the fission operator on an LD problem (fixed IN Step A, ruled 2026-11-12) ──────────────────────────────


@pytest.mark.parametrize("name", ("slab_ld_p1", "2d_ld"))
@pytest.mark.xfail(strict=True, reason="#707 Step A: red by design until the true-ends carve (ontology chart W1 to W4); the flip-edit removes this marker")
def test_ldf_fission_acts_on_an_ld_problem(name: str) -> None:
    """LD-F RED. ``F`` on an LD problem returns the integrated, mass-weighted isotropic emission:
    ``(Fψ)[n, g, i, o] = V_i θ^{|o|} χ_g(i) Σ_{g'} νΣf_{g'}(i) φ_{g', i, o} / W``, with ``φ_o = Σ_n w_n ψ_{n, o}`` and
    ``W = Σ_n w_n`` (per steradian, R8; mass-weighted moments, Q-D).

    Oracle: raw χ, νΣf, V, θ, w_n. First red today: ``ValueError`` "FissionMaterialField gather: the layout indexes
    1 spatial axes but the requested shape (10, 2) has 2" (``transport/material_field.py:483``): the in-cell moment
    axis is read as a spatial axis, so F, the LD eigen question and the chart-law LD fission rows are all unreachable.
    """
    problem = _problem(name)
    F = cc.atoms(problem)["FissionOperator"]
    psi = _psi(name)
    try:
        out = np.asarray(F.apply(psi).interior.values, float)
    except ValueError as exc:
        pytest.fail(f"F raises on the LD problem {name}: {exc}")
    mats = cc._mixtures()
    mesh_mats = np.asarray(cc._slab_mesh().mat_ids if name == "slab_ld_p1" else cc._mesh_2d().mat_map)  # per cell
    w = np.asarray(problem.quad.weights, float)
    vals = np.asarray(psi.interior.values, float)  # (N, ng, *spatial, M)
    phi = np.einsum("n,n...->...", w, vals)  # (ng, *spatial, M)
    chi = np.stack([np.asarray(mats[int(k)].chi, float) for k in mesh_mats.ravel()], axis=-1)  # (ng, cells)
    nusf = np.stack([np.asarray(mats[int(k)].SigP, float) for k in mesh_mats.ravel()], axis=-1)  # (ng, cells)
    spatial = vals.shape[2:-1]
    chi = chi.reshape((chi.shape[0],) + spatial)
    nusf = nusf.reshape((nusf.shape[0],) + spatial)
    production = np.einsum("g...,g...m->...m", nusf, phi)  # (*spatial, M)
    emission = chi[..., None] * production[None] / w.sum()  # (ng, *spatial, M)
    want = cc.sn_bulk_measure(problem, out.shape) * emission[None]
    scale = float(np.max(np.abs(want)))
    assert float(np.max(np.abs(out - want))) <= 256 * np.finfo(float).eps * scale


@pytest.mark.xfail(strict=True, reason="#707 Step A: red by design until the true-ends carve (ontology chart W1 to W4); the flip-edit removes this marker")
def test_ldf_the_ld_eigen_question_is_posable() -> None:
    """LD-F RED. The LD eigen question reads F: ``eigen_posing.rayleigh`` on an LD problem returns a finite number.
    (There is no public LD eigen entry: ``solve_sn`` takes no ``scheme``; this is the posing it would solve.) First red
    today: the same ``ValueError`` from ``material_field.py:483``."""
    from orpheus.numerics.coupled_system import CoupledField

    problem = _problem("slab_ld_p1")
    try:
        k = problem.eigen_posing.rayleigh(CoupledField(systems=(_psi("slab_ld_p1"),)), w=1.0)
    except ValueError as exc:
        pytest.fail(f"the LD eigen question raises: {exc}")
    assert np.isfinite(k)


@cache
def cc_anchor() -> dict[str, np.ndarray]:
    with np.load(cc.ANCHOR) as data:
        return {k: data[k] for k in data.files}
