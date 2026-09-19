r"""The OUTCOMES' defining laws — the kind is the posing's TYPE, fused with its answer and gauge.

Consumers campaign step 3, U1 (2026-09-17; ``consumers_step3_design.md`` §3.2, the
test-architect's rows 1.1–1.5, 1.13, 1.14). Shape (B) at the Solution tier: two
kind-typed outcomes and ZERO Optionals, ONE fused member so every illegal pairing
is UNCONSTRUCTIBLE (Pattern 4) — an eigen question beside a kernel gauge, a source
question beside a λ, a state beside a question of the other kind.

The ``Evidence`` sum is asserted CLOSED (pyright's ``assert_never`` on an exhaustive
match); its members' reachability through the real entries is U2's row, not this
file's (the certificate EVALUATOR lands with the SN mint).

``CoupledField.space`` (row 1.14) closes an inert ends check: ``[M]`` until this
step every SN source was lifted to a one-system coupled state that had NO
``.space``, so ``SourcePosing.__post_init__``'s "the source lives on the
operator's codomain" law refused nothing. Now it refuses a cross-hub source.
"""

from __future__ import annotations

from typing import assert_never

import dataclasses

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry.coord import CoordSystem
from orpheus.geometry.mesh import BC, Mesh1D
from orpheus.homogeneous.solver import HomogeneousProblem
from orpheus.numerics.gauge import ScaleGauge
from orpheus.numerics.outcome import (
    Certified,
    EigenOutcome,
    Evidence,
    ExitCertificate,
    Measured,
    NotApplicable,
    NotYet,
    SourceOutcome,
)
from orpheus.numerics.posing import ALPHA_MAP, EigenPosing, SourcePosing
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solver import _as_problem, _build_fixed_source_rhs

pytestmark = pytest.mark.foundation


def _require(condition: object, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def _zero_d(fam: str = "A", ng: str = "2g"):
    """The 0-D problem, its exact k∞ and eigenvector as the posed (ng, 1) column."""
    problem = HomogeneousProblem(get_mixture(fam, ng))
    A = np.asarray(problem.loss.as_matrix(), dtype=float)
    F = np.asarray(problem.production.as_matrix(), dtype=float)
    w, V = np.linalg.eig(np.linalg.solve(A, F))
    i = int(np.argmax(np.real(w)))
    phi = np.real(V[:, i])
    phi = phi if phi.sum() >= 0 else -phi
    return problem, float(np.real(w[i])), phi.reshape(-1, 1)


def _eigen_outcome(fam: str = "A", ng: str = "2g") -> tuple[HomogeneousProblem, EigenOutcome]:
    problem, k, phi = _zero_d(fam, ng)
    gauge = ScaleGauge(problem.production_rate.evaluate, 100.0)
    return problem, EigenOutcome(problem.eigen_posing, gauge.apply(phi), k, (k,), gauge)


_EPS = np.finfo(float).eps


# ── the kind is the posing's type ────────────────────────────────────────


class TestLawTheKindIsThePosingsType:
    def test_law_an_eigen_outcome_constructs_and_reads_lam_as_keff(self) -> None:
        problem, out = _eigen_outcome()
        _require(isinstance(out.posing, EigenPosing), "the question is the eigen one")
        _require(out.keff == out.lam, "keff IS lam under the k map")
        _require(out.adjoint_posing.pencil.lhs.H is out.posing.pencil.lhs, "the nullary adjoint daggers the same pencil (A.H.H is A by identity)")

    def test_negative_an_eigen_outcome_REFUSES_a_source_question(self) -> None:
        problem, k, phi = _zero_d()
        gauge = ScaleGauge(problem.production_rate.evaluate, 100.0)
        with pytest.raises(TypeError, match="must be an EigenPosing"):
            EigenOutcome(SourcePosing(problem.pencil.lhs, phi), phi, k, (k,), gauge)  # type: ignore[arg-type]  # the refusal pyright also sees

    def test_negative_a_source_outcome_REFUSES_an_eigen_question(self) -> None:
        problem, k, phi = _zero_d()
        with pytest.raises(TypeError, match="must be a SourcePosing"):
            SourceOutcome(problem.eigen_posing, phi, _TrivialKernelGauge())  # type: ignore[arg-type]

    def test_law_a_source_outcome_has_no_eigen_verbs(self) -> None:
        problem, k, phi = _zero_d()
        A = np.asarray(problem.pencil.lhs.as_matrix(), dtype=float)
        q = np.ones((A.shape[0], 1))
        psi = np.linalg.solve(A, q)
        out = SourceOutcome(SourcePosing(problem.pencil.lhs, q), psi, _TrivialKernelGauge())
        _require(not hasattr(out, "keff") and not hasattr(out, "lam") and not hasattr(out, "adjoint_posing"), "a coset has no λ and no nullary adjoint — the attributes do not exist")
        _require(float(np.max(np.abs(np.asarray(out.residual())))) <= 1e-12, "its residual closes at the solution")

    def test_law_the_trajectory_is_co_indexed_with_lam(self) -> None:
        problem, k, phi = _zero_d()
        gauge = ScaleGauge(problem.production_rate.evaluate, 100.0)
        EigenOutcome(problem.eigen_posing, phi, k, (0.9 * k, k), gauge)
        with pytest.raises(ValueError, match="last entry"):
            EigenOutcome(problem.eigen_posing, phi, k, (k, 0.9 * k), gauge)

    def test_law_keff_is_refused_under_a_non_k_map(self) -> None:
        problem, k, phi = _zero_d()
        gauge = ScaleGauge(problem.production_rate.evaluate, 100.0)
        alpha = EigenOutcome(EigenPosing(problem.pencil, ALPHA_MAP), phi, -0.3, (-0.3,), gauge)
        _require(alpha.lam == -0.3, "lam is the physical eigenvalue under whichever map")
        with pytest.raises(ValueError, match="not 'k'"):
            alpha.keff


class TestLawTheOutcomeAgreesWithItsQuestion:
    @pytest.mark.parametrize("ng", ["1g", "2g", "4g"])
    def test_law_rayleigh_agrees_with_lam_on_the_exact_pair(self, ng: str) -> None:
        """``[M]`` 2g/4g gap 2.22e-16, 1g exactly 0 — a dense ``eig`` vs ONE division of two pairings."""
        problem, out = _eigen_outcome("A", ng)
        _require(abs(out.rayleigh() - out.lam) <= 8 * _EPS * abs(out.lam), f"rayleigh − lam = {out.rayleigh() - out.lam:.3e}")
        _require(abs(out.balance()) <= 64 * _EPS * float(np.max(np.abs(out.state))), "β(ψ, λ) = 0 at the exact pair")

    def test_negative_a_perturbed_lam_is_seen(self) -> None:
        problem, k, phi = _zero_d()
        gauge = ScaleGauge(problem.production_rate.evaluate, 100.0)
        wrong = EigenOutcome(problem.eigen_posing, phi, k * (1 + 1e-9), (k * (1 + 1e-9),), gauge)
        _require(abs(wrong.rayleigh() - wrong.lam) > 1e-10 * abs(k), "a 1e-9 perturbation of λ is visible against the posing's own quotient")

    def test_law_the_dagger_is_an_involution_on_the_posing(self) -> None:
        """No negative leg: the involution is a theorem of the adjoint wrapper
        (``[M]`` ``A.H.H is A`` by identity, ``operator.py:1674-1678``)."""
        problem, out = _eigen_outcome()
        twice = out.adjoint_posing.H()
        _require(twice.pencil.lhs is out.posing.pencil.lhs and twice.pencil.rhs is out.posing.pencil.rhs, "H∘H = id on the pencil's ends")


# ── the evidence sum ─────────────────────────────────────────────────────


def _describe(evidence: Evidence) -> str:
    match evidence:
        case Measured(value=v):
            return f"measured {v}"
        case Certified(bound=b, by=by):
            return f"certified ≤ {b} by {by}"
        case NotApplicable(reason=r):
            return f"n/a: {r}"
        case NotYet(issue=i, reason=r):
            return f"not yet (#{i}): {r}"
        case _:
            assert_never(evidence)


class TestLawTheEvidenceSum:
    def test_law_every_member_constructs_and_the_sum_is_closed(self) -> None:
        members = (Measured(0.31), Certified(1e-7, "within-group certificate"), NotApplicable("zero source"), NotYet(354, "the carrying eigen exit's coupled rhs"))
        for m in members:
            _require(_describe(m), f"{m!r} is described by the exhaustive match")
        cert = ExitCertificate(balance=members[0], gauge=members[2], rayleigh_gap=members[0], admissibility=members[3])
        _require(cert.balance == Measured(0.31), "members are values")
        _require(not hasattr(cert, "value_or_none"), "no accessor re-imports the None leak")


# ── CoupledField.space closes the ends check ─────────────────────────────


def _slab(length: float = 4.0):
    materials = {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")}
    mesh = Mesh1D(edges=np.linspace(0.0, length, 9), mat_ids=np.array([0] * 4 + [1] * 4, dtype=int), bc_left=BC("reflective"), bc_right=BC("vacuum"))
    return _as_problem(mesh, Quadrature.gauss_legendre(8), materials)


def _carrying_sphere():
    mesh = Mesh1D(edges=np.linspace(0.0, 3.0, 7), mat_ids=np.zeros(6, dtype=int), coord=CoordSystem.SPHERICAL, bc_right=BC("vacuum"))
    return _as_problem(mesh, Quadrature.gauss_legendre(8), {0: get_mixture("A", "2g")})


class TestLawTheCoupledStateKnowsItsSpace:
    @pytest.mark.parametrize("hub_of, nx", [(_slab, 8), (_carrying_sphere, 6)], ids=["slab (1 system)", "carrying sphere (2 systems)"])
    def test_law_the_source_lives_on_the_pencils_domain(self, hub_of, nx: int) -> None:
        hub = hub_of()
        q = _build_fixed_source_rhs(np.ones((hub.quad.N, 2, nx)), hub)
        posing = hub.source_posing(q)
        _require(posing.source.space == hub.pencil.at(1.0).domain, "the lifted source's DERIVED space equals the pencil's domain (content identity)")
        _require(posing.source.space == hub.system.space, "…and the builder-minted coupled space")
        _require(hub.system.space.zeros().space == hub.system.space, "the space's own zero element lives on it")

    def test_negative_a_cross_hub_source_is_REFUSED_at_the_posing(self) -> None:
        """The ends check is LIVE now: a source built on another Problem's space
        cannot be posed on this one (until step 3 the coupled state had no
        ``.space`` and this raised nothing)."""
        here, there = _slab(4.0), _slab(6.0)
        foreign = _build_fixed_source_rhs(np.ones((there.quad.N, 2, 8)), there)
        with pytest.raises(ValueError, match="pose q on the operator's ends"):
            here.source_posing(foreign)


class _TrivialKernelGauge:
    """A zero-block section: the identity, dimension 0 (the honest 'no freedom' value)."""

    def gauge(self, trace):
        return trace

    @property
    def dimension(self) -> int:
        return 0


# ── the dominance ratio is a reading of the outcome's trajectory ──────────
#
# Migrated from the retired SN view's rows (step 3 U6): the late-iteration
# ratio |k_n − k_{n−1}| / |k_{n−1}| is a reading of the TRAJECTORY the eigen
# outcome carries, so it lives on the outcome, not on a diagnostics view.


def _with_trajectory(trajectory: tuple[float, ...]) -> EigenOutcome:
    _, out = _eigen_outcome()
    return dataclasses.replace(out, lam=trajectory[-1], trajectory=trajectory)


class TestDominanceRatioIsTheOutcomes:
    def test_a_single_entry_has_no_ratio(self) -> None:
        assert _with_trajectory((1.0,)).dominance_ratio() is None

    def test_three_iterates_read_the_last_step(self) -> None:
        ratio = _with_trajectory((1.0, 1.05, 1.1)).dominance_ratio()
        assert ratio == pytest.approx(abs(1.1 - 1.05) / 1.05)

    def test_a_zero_previous_iterate_has_no_ratio(self) -> None:
        assert _with_trajectory((0.0, 1.0)).dominance_ratio() is None
