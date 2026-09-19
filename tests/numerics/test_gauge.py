r"""The GAUGES' defining laws — a section of a torsor's quotient, realized twice.

Consumers campaign step 3, U1 (2026-09-17; ``consumers_step3_design.md`` §3.2, the
test-architect's rows 1.6–1.12). A solve returns a REPRESENTATIVE of a solution
SET carrying a free transitive group action — the eigen ray under
:math:`(\mathbb{R}_+, \times)`, the source coset under :math:`(\ker A, +)` — and
the gauge is the SECTION that picked it. Each realization ships a test of its
defining laws (``feedback_test_intrinsic_properties``):

* it LANDS on its target (``n(apply(ψ)) ≈ t``; ``Π(gauge(ψ)) ≈ 0``);
* it is IDEMPOTENT — ``allclose``, deliberately NOT bit-equal for the scale gauge:
  the second rescale multiplies by :math:`1/(1 \pm \varepsilon)` (a first draft of
  the frame pass claimed 0 ULP; refuted by arithmetic);
* it is Γ-INVARIANT — the group element for the kernel gauge is drawn from
  :math:`\Pi`'s OWN RANGE on the TRACE space (the kernel is a trace object —
  ``loss_kernel_gauge.py``'s class docstring — so an interior perturbation is
  not a group element and must not be used as one; the frame pass's "interior
  kernel modes survive" leg was refuted on exactly this point);
* the kernel gauge is RESIDUAL-NEUTRAL: :math:`A(\psi - \Pi\psi) = A\psi` because
  :math:`\Pi\psi \in \ker A` (the production claim at ``_exit_gauge_trace``).

And the NEGATIVE law that motivated the type: two gauges of one ray under two
functionals are DISTINGUISHABLE values — until this step the tree applied four
"production rate" functionals and recorded none of them (``[M]``
``eigenvalue.py:437-447``, ``homogeneous/solver.py:436``, ``eigenvalue.py:506``).

Fixtures: the 0-D ``HomogeneousProblem`` pencil (≥ 2G — ``[M]`` the reaction-rate
functional silently accepts an ``(ng,)`` vector and returns the wrong number
where the posed ``(ng, 1)`` column reads the right one; at 1G they coincide) and
the entry ledger's gauge-singular all-reflective 2-D box (``dim ker A = 12``).
"""

from __future__ import annotations

from dataclasses import replace
from typing import cast

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture, make_mixture
from orpheus.geometry import BC, Mesh2D
from orpheus.homogeneous.solver import HomogeneousProblem
from orpheus.numerics.coupled_system import CoupledField
from orpheus.numerics.gauge import KernelGauge, ScaleGauge
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solver import _as_problem
from orpheus.transport.full_field import FullField

pytestmark = pytest.mark.foundation


def _require(condition: object, message: str) -> None:
    if not condition:
        raise AssertionError(message)


# ── fixtures ─────────────────────────────────────────────────────────────


def _zero_d(fam: str = "A", ng: str = "2g"):
    """The 0-D problem and its exact gauged eigenvector as the posed (ng, 1) column."""
    problem = HomogeneousProblem(get_mixture(fam, ng))
    A = np.asarray(problem.loss.as_matrix(), dtype=float)
    F = np.asarray(problem.production.as_matrix(), dtype=float)
    w, V = np.linalg.eig(np.linalg.solve(A, F))
    i = int(np.argmax(np.real(w)))
    phi = np.real(V[:, i])
    phi = phi if phi.sum() >= 0 else -phi
    return problem, phi.reshape(-1, 1)


_R = BC("reflective")
_QUAD = Quadrature.level_symmetric(sn_order=4)


def _absorber(ng: int = 2):
    sig_t = np.linspace(0.8, 1.6, ng)
    return make_mixture(
        sig_t=sig_t, sig_c=sig_t.copy(), sig_f=np.zeros(ng),
        nu=np.zeros(ng), chi=np.zeros(ng), sig_s=np.zeros((ng, ng)),
    )


def _gauge_singular_hub():
    """The ledger's all-reflective (3, 4) box — ``[M]`` ``dim ker A = 12`` (every
    mesh there is singular; parity only decides whether a SYMMETRIC source
    excites it, which is irrelevant to a law stated on a random trace)."""
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 1.0, 4), edges_y=np.linspace(0.0, 2.0, 5),
        mat_map=np.zeros((3, 4), dtype=int),
        bc_xmin=_R, bc_xmax=_R, bc_ymin=_R, bc_ymax=_R,
    )
    return _as_problem(mesh, _QUAD, {0: _absorber()})


def _random_state(hub, rng):
    """A random 1-system coupled state on the hub — interior AND trace populated."""
    zero = hub.system.space.zeros()
    return CoupledField.from_flat(rng.random(zero.to_flat().size), zero)


def _member(state: CoupledField) -> FullField:
    """System A's composite (the member contract is structural; the fixture's is a FullField)."""
    return cast(FullField, state.systems[0])


def _trace_values(state: CoupledField) -> np.ndarray:
    return np.asarray(_member(state).boundary.values, dtype=float)


def _with_trace(state: CoupledField, trace_values: np.ndarray) -> CoupledField:
    member = _member(state)
    return CoupledField(systems=(replace(member, boundary=replace(member.boundary, values=trace_values)),))


# ── ScaleGauge ───────────────────────────────────────────────────────────


class TestLawTheScaleGauge:
    @pytest.mark.parametrize("ng", ["2g", "4g"])
    def test_law_the_section_lands_on_its_target(self, ng: str) -> None:
        problem, phi = _zero_d("A", ng)
        gauge = ScaleGauge(problem.production_rate.evaluate, 100.0)
        gauged = gauge.apply(3.7 * phi)
        _require(abs(gauge.functional(gauged) - 100.0) <= 1e-12 * 100.0, "n(apply(ψ)) = t to 1e-12 (one multiplication, no reduction)")
        _require(gauged.shape == phi.shape, "the section preserves the posed (ng, 1) shape")

    def test_law_idempotence_to_allclose_not_to_the_bit(self) -> None:
        problem, phi = _zero_d()
        gauge = ScaleGauge(problem.production_rate.evaluate, 100.0)
        once = gauge.apply(phi)
        twice = gauge.apply(once)
        _require(np.allclose(twice, once, rtol=1e-14, atol=0.0), "apply∘apply = apply to rtol 1e-14 — the second rescale is ×1/(1±ε)")

    @pytest.mark.parametrize("c", [0.1, 0.5, 2.0, 10.0])
    def test_law_gamma_invariance_under_the_scale_group(self, c: float) -> None:
        problem, phi = _zero_d()
        gauge = ScaleGauge(problem.production_rate.evaluate, 100.0)
        _require(np.allclose(gauge.apply(c * phi), gauge.apply(phi), rtol=1e-12, atol=0.0), "the section is invariant under ψ ↦ cψ, c ∈ ℝ₊")

    def test_law_the_sign_leg_is_one_case_of_the_section(self) -> None:
        """The dense engines' ℤ/2 fix (``_sign_normalised``) is subsumed: a
        NEGATIVE functional value gives a negative displacement and ``apply``
        lands on the positive representative."""
        problem, phi = _zero_d()
        gauge = ScaleGauge(problem.production_rate.evaluate, 100.0)
        _require(gauge.displacement(-phi) < 0.0, "a negative reading gives a negative displacement")
        _require(np.allclose(gauge.apply(-phi), gauge.apply(phi), rtol=1e-12, atol=0.0), "the sign is fixed by the section, not by a second convention")

    def test_negative_a_zero_functional_reading_has_no_section(self) -> None:
        problem, phi = _zero_d()
        gauge = ScaleGauge(lambda state: 0.0, 100.0)
        with pytest.raises(ValueError, match="no section"):
            gauge.displacement(phi)

    def test_negative_two_sections_under_two_functionals_are_distinguishable(self) -> None:
        """The law that motivated the type — the tree applied four functionals
        under one name and recorded none. Two gauges of one ray under two
        functionals are DIFFERENT values (the TYPE leg) and read DIFFERENT
        numbers on the same state (the VALUE leg)."""
        problem, phi = _zero_d()
        by_production = ScaleGauge(problem.production_rate.evaluate, 1.0)
        by_absorption = ScaleGauge(problem.absorption_rate.evaluate, 1.0)
        _require(by_production != by_absorption, "the functional is part of the value — the two sections are distinguishable")
        _require(by_production == ScaleGauge(problem.production_rate.evaluate, 1.0), "…and equal to a section built from the same functional and target")
        p, a = by_production.functional(phi), by_absorption.functional(phi)
        _require(abs(p - a) > 1e-2 * abs(p), f"the two functionals read different numbers on one state ({p:.6g} vs {a:.6g})")
        _require(not np.allclose(by_production.apply(phi), by_absorption.apply(phi)), "so the two representatives differ by a positive scalar the recorded gauge names")


# ── KernelGauge (the Protocol, realized by the SN loss-kernel gauge) ─────


class TestLawTheKernelGauge:
    def test_law_the_sn_gauge_satisfies_the_protocol_structurally(self) -> None:
        hub = _gauge_singular_hub()
        _require(isinstance(hub.loss_kernel_gauge, KernelGauge), "LossKernelGauge conforms to KernelGauge without inheriting (numerics never imports sn)")
        _require(hub.loss_kernel_gauge.dimension > 0, "non-vacuity: the fixture's kernel is non-trivial")

    def test_law_idempotence_of_the_projector_and_the_section(self) -> None:
        hub = _gauge_singular_hub()
        gauge = hub.loss_kernel_gauge
        rng = np.random.default_rng(3)
        trace = rng.random(_trace_values(hub.system.space.zeros()).shape)
        pi = gauge.apply(trace)
        _require(float(np.linalg.norm(pi)) > 0.0, "non-vacuity: a random trace has a kernel component")
        _require(np.allclose(gauge.apply(pi), pi, rtol=1e-12, atol=1e-14), "Π² = Π")
        _require(np.allclose(gauge.gauge(gauge.gauge(trace)), gauge.gauge(trace), rtol=1e-12, atol=1e-14), "(I − Π)² = I − Π")

    def test_law_gamma_invariance_with_a_group_element_from_the_kernel(self) -> None:
        """The group element is drawn from Π's OWN RANGE on the TRACE space."""
        hub = _gauge_singular_hub()
        gauge = hub.loss_kernel_gauge
        rng = np.random.default_rng(5)
        shape = _trace_values(hub.system.space.zeros()).shape
        trace = rng.random(shape)
        v = gauge.apply(rng.random(shape))
        _require(float(np.linalg.norm(v)) > 0.0, "non-vacuity: the group element is non-zero")
        _require(np.allclose(gauge.gauge(trace + v), gauge.gauge(trace), rtol=1e-12, atol=1e-14), "the section is invariant under ψ ↦ ψ + v, v ∈ ker A")

    def test_law_residual_neutrality_of_the_section(self) -> None:
        """A(ψ − Πψ) = Aψ because Πψ ∈ ker A — the full residual VECTOR, not a
        mirror-even projection of it."""
        hub = _gauge_singular_hub()
        gauge = hub.loss_kernel_gauge
        rng = np.random.default_rng(7)
        state = _random_state(hub, rng)
        trace = _trace_values(state)
        gauged = _with_trace(state, gauge.gauge(trace))
        _require(float(np.linalg.norm(trace - _trace_values(gauged))) > 0.0, "non-vacuity: the section moved the trace")
        before = hub.system.loss.apply(state).to_flat()
        after = hub.system.loss.apply(gauged).to_flat()
        scale = float(np.max(np.abs(before)))
        _require(float(np.max(np.abs(before - after))) <= 1e-12 * scale, "the loss residual is unchanged by the section to 1e-12 relative")

    def test_negative_a_projector_off_the_kernel_is_not_residual_neutral(self) -> None:
        """The law has teeth: moving the trace by a NON-kernel vector moves the residual."""
        hub = _gauge_singular_hub()
        rng = np.random.default_rng(11)
        state = _random_state(hub, rng)
        trace = _trace_values(state)
        moved = _with_trace(state, trace + rng.random(trace.shape))
        before = hub.system.loss.apply(state).to_flat()
        after = hub.system.loss.apply(moved).to_flat()
        _require(float(np.max(np.abs(before - after))) > 1e-6 * float(np.max(np.abs(before))), "a non-kernel move IS visible in the residual")
