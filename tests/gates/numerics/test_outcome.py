r"""The OUTCOMES' defining laws — the kind is the posing's TYPE, fused with its answer and gauge.

Consumers campaign step 3, U1 (2026-09-17; ``consumers_step3_design.md`` §3.2,
rows 1.1–1.5, 1.13, 1.14); re-based layer DAG Phase 2
(``scratch/layer_dag/refiling/specs_2r/numerics_laws.md``) onto
:mod:`tests.gates.numerics._energy_pencil`. Shape (B) at the Solution tier: two
kind-typed outcomes and ZERO Optionals, ONE fused member so every illegal pairing
is UNCONSTRUCTIBLE (Pattern 4).

The ``Evidence`` sum is asserted CLOSED (pyright's ``assert_never`` on an exhaustive
match).

Row 1.14, the ends check of ``SourcePosing``: the LAW (a source on another space
of the same shape is refused; one on the operator's codomain is admitted) is
posed here on numerics spaces; the SN half — the SN-lifted coupled source
carries a ``.space`` equal to the hub's — moved to
``tests/gates/sn/architecture/test_posing.py``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import assert_never

import dataclasses

import numpy as np
import pytest

from tests.gates.numerics._energy_pencil import FISSILE, POINT, EnergyPencil

from orpheus.numerics.axis import EnergyAxis
from orpheus.numerics.field import Field
from orpheus.numerics.gauge import ScaleGauge
from orpheus.numerics.evidence import Asserted, Evidence, Measured, NotApplicable, NotYet
from orpheus.numerics.fredholm import Solve
from orpheus.numerics.outcome import EigenOutcome, ExitReport, SourceOutcome
from orpheus.numerics.posing import ALPHA_MAP, EigenPosing, SourcePosing
from orpheus.numerics.space import FunctionSpace

pytestmark = pytest.mark.foundation


def _require(condition: object, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def _zero_d(ng: int = 2) -> tuple[EnergyPencil, float, np.ndarray]:
    """The numerics pencil, its exact k and eigenvector."""
    problem = EnergyPencil.of(FISSILE[ng])
    k, phi = problem.exact_pair
    return problem, k, phi


def _eigen_outcome(ng: int = 2) -> tuple[EnergyPencil, EigenOutcome]:
    problem, k, phi = _zero_d(ng)
    gauge = ScaleGauge(problem.production_rate, 100.0)
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
        gauge = ScaleGauge(problem.production_rate, 100.0)
        with pytest.raises(TypeError, match="must be an EigenPosing"):
            EigenOutcome(SourcePosing(problem.pencil.lhs, phi), phi, k, (k,), gauge)  # type: ignore[arg-type]  # the refusal pyright also sees

    def test_negative_a_source_outcome_REFUSES_an_eigen_question(self) -> None:
        problem, k, phi = _zero_d()
        with pytest.raises(TypeError, match="must be a SourcePosing"):
            SourceOutcome(problem.eigen_posing, phi)  # type: ignore[arg-type]

    def test_law_a_source_outcome_has_no_eigen_verbs(self) -> None:
        problem, k, phi = _zero_d()
        A = np.asarray(problem.pencil.lhs.as_matrix(), dtype=float)
        q = np.ones(problem.space.shape)
        psi = np.linalg.solve(A, q.ravel()).reshape(problem.space.shape)
        out = SourceOutcome(SourcePosing.of_pencil(problem.pencil, 0.0, q), psi)
        _require(out.gauge.dimension == 0, "the gauge is READ off the posing's kernel (the fixture's pencil predicts it trivial)")
        _require(not hasattr(out, "keff") and not hasattr(out, "lam") and not hasattr(out, "adjoint_posing"), "a coset has no λ and no nullary adjoint — the attributes do not exist")
        _require(float(np.max(np.abs(np.asarray(out.residual())))) <= 1e-12, "its residual closes at the solution")

    def test_law_the_trajectory_is_co_indexed_with_lam(self) -> None:
        problem, k, phi = _zero_d()
        gauge = ScaleGauge(problem.production_rate, 100.0)
        EigenOutcome(problem.eigen_posing, phi, k, (0.9 * k, k), gauge)
        with pytest.raises(ValueError, match="last entry"):
            EigenOutcome(problem.eigen_posing, phi, k, (k, 0.9 * k), gauge)

    def test_law_keff_is_refused_under_a_non_k_map(self) -> None:
        problem, k, phi = _zero_d()
        gauge = ScaleGauge(problem.production_rate, 100.0)
        alpha = EigenOutcome(EigenPosing(problem.pencil, ALPHA_MAP), phi, -0.3, (-0.3,), gauge)
        _require(alpha.lam == -0.3, "lam is the physical eigenvalue under whichever map")
        with pytest.raises(ValueError, match="not 'k'"):
            alpha.keff


class TestLawTheOutcomeAgreesWithItsQuestion:
    @pytest.mark.parametrize("ng", [1, 2, 4])
    @pytest.mark.rests_on(
        "tests/gates/numerics/test_posing.py::TestLawTheEigenPosing::test_law_the_rayleigh_at_the_solution_is_the_eigenvalue",
    )
    def test_law_rayleigh_agrees_with_lam_on_the_exact_pair(self, ng: int) -> None:
        """``[M]`` 2g/4g gap 2.22e-16, 1g exactly 0 — a dense ``eig`` vs ONE division of two pairings."""
        problem, out = _eigen_outcome(ng)
        _require(abs(out.rayleigh() - out.lam) <= 8 * _EPS * abs(out.lam), f"rayleigh − lam = {out.rayleigh() - out.lam:.3e}")
        _require(abs(out.balance()) <= 64 * _EPS * float(np.max(np.abs(out.state))), "β(ψ, λ) = 0 at the exact pair")

    @pytest.mark.rests_on(
        "tests/gates/numerics/test_posing.py::TestLawTheEigenPosing::test_law_the_rayleigh_at_the_solution_is_the_eigenvalue",
    )
    def test_negative_a_perturbed_lam_is_seen(self) -> None:
        problem, k, phi = _zero_d()
        gauge = ScaleGauge(problem.production_rate, 100.0)
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
        case Asserted(bound=b, by=by):
            return f"certified ≤ {b} by {by}"
        case NotApplicable(reason=r):
            return f"n/a: {r}"
        case NotYet(issue=i, reason=r):
            return f"not yet (#{i}): {r}"
        case _:
            assert_never(evidence)


class TestLawTheEvidenceSum:
    def test_law_every_member_constructs_and_the_sum_is_closed(self) -> None:
        members = (Measured(0.31), Asserted(1e-7, "within-group claim check"), NotApplicable("zero source"), NotYet(354, "the carrying eigen exit's coupled rhs"))
        for m in members:
            _require(_describe(m), f"{m!r} is described by the exhaustive match")
        cert = ExitReport(balance=members[0], gauge=members[2], rayleigh_gap=members[0], admissibility=members[3],
                          solvability=Solve(because="a fixture", conditioning=members[2]))
        _require(cert.balance == Measured(0.31), "members are values")
        _require(not hasattr(cert, "value_or_none"), "no accessor re-imports the None leak")


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


# ── the ends check of the source posing (row 1.14's law) ─────────────────


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class _Source(Field):
    """A concrete field: values on a space (the source's ``.space`` is what the ends check reads)."""


def _partitioned_twin(space: FunctionSpace) -> FunctionSpace:
    """A space of EQUAL shape and different content (a real group partition)."""
    ng = space.shape[0]
    return FunctionSpace.of_axes(EnergyAxis.of_edges(np.geomspace(2.0e7, 1.0e-5, ng + 1)), POINT)


class TestLawTheSourceLivesOnTheOperatorsEnds:
    def test_law_a_source_on_the_codomain_is_ADMITTED(self) -> None:
        problem, _, _ = _zero_d()
        member = problem.pencil.at(1.0)
        _require(member.domain == problem.space, "the member's DERIVED domain is the pencil's space (the sum's agreement law)")
        q = _Source(values=np.ones(problem.space.shape), space=problem.space)
        posing = SourcePosing(member, q)
        _require(posing.source.space == member.codomain, "the posed source lives on the member's codomain")

    def test_negative_a_foreign_source_is_REFUSED_at_the_posing(self) -> None:
        """Equal SHAPE, different space: a shape-keyed check would admit it."""
        problem, _, _ = _zero_d()
        foreign = _partitioned_twin(problem.space)
        _require(foreign.shape == problem.space.shape and foreign != problem.space, "the sharp witness: equal shapes, different spaces")
        with pytest.raises(ValueError, match="pose q on the operator's ends"):
            SourcePosing(problem.pencil.at(1.0), _Source(values=np.ones(foreign.shape), space=foreign))
