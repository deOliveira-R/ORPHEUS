r"""The SCALE gauge's defining laws — a section of the :math:`(\mathbb{R}_+, \times)` quotient.

Consumers campaign step 3, U1 (2026-09-17; ``consumers_step3_design.md`` §3.2,
rows 1.6–1.9); re-based layer DAG Phase 2
(``scratch/layer_dag/refiling/specs_2r/numerics_laws.md``) onto
:mod:`tests.gates.numerics._energy_pencil`. A solve returns a REPRESENTATIVE of
the eigen ray, and the gauge is the SECTION that picked it:

* it LANDS on its target (``n(apply(ψ)) ≈ t``);
* it is IDEMPOTENT — ``allclose``, deliberately NOT bit-equal: the second
  rescale multiplies by :math:`1/(1 \pm \varepsilon)`;
* it is Γ-INVARIANT under :math:`\psi \mapsto c\psi`, :math:`c > 0`;
* the sign leg is one case of the section.

And the NEGATIVE law that motivated the type: two gauges of one ray under two
functionals are DISTINGUISHABLE values.

The KERNEL gauge's laws are numerics' own since the space-membership campaign's U3:
:class:`~orpheus.numerics.kernel_projection_gauge.KernelProjectionGauge` realises
:class:`~orpheus.numerics.gauge.KernelGauge` and is gated on a numerics fixture with a
non-uniform metric in ``tests/gates/numerics/test_kernel_projection_gauge.py``; the SN
hub's instance (SN supplies the generators) is gated in
``tests/gates/sn/operators/test_loss_kernel_gauge.py``. The five kernel-gauge rows this
file carried on the SN hub were retired as dominated by those two files (layer DAG
Phase 2 refresh, ``specs_2r/numerics_laws.md`` §1.1, with the arms that red each).

The rate functionals are co-vectors (``space.flat(w)``, :ref:`spaces-covector`); the
gauge reads their ``evaluate``.
"""

from __future__ import annotations


import numpy as np
import pytest

from tests.gates.numerics._energy_pencil import FISSILE, EnergyPencil

from orpheus.numerics.gauge import ScaleGauge

pytestmark = pytest.mark.foundation


def _require(condition: object, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def _ray(ng: int = 2) -> tuple[EnergyPencil, np.ndarray]:
    """The numerics pencil and its exact eigenvector (a representative of the ray)."""
    p = EnergyPencil.of(FISSILE[ng])
    return p, p.exact_pair[1]


class TestLawTheScaleGauge:
    @pytest.mark.parametrize("ng", [2, 4])
    def test_law_the_section_lands_on_its_target(self, ng: int) -> None:
        p, phi = _ray(ng)
        gauge = ScaleGauge(p.production_rate, 100.0)
        gauged = gauge.apply(3.7 * phi)
        _require(abs(gauge.functional(gauged) - 100.0) <= 1e-12 * 100.0, "n(apply(ψ)) = t to 1e-12 (one multiplication, no reduction)")
        _require(gauged.shape == phi.shape, "the section preserves the carrier's shape")

    def test_law_idempotence_to_allclose_not_to_the_bit(self) -> None:
        p, phi = _ray()
        gauge = ScaleGauge(p.production_rate, 100.0)
        once = gauge.apply(phi)
        twice = gauge.apply(once)
        _require(np.allclose(twice, once, rtol=1e-14, atol=0.0), "apply∘apply = apply to rtol 1e-14 — the second rescale is ×1/(1±ε)")

    @pytest.mark.parametrize("c", [0.1, 0.5, 2.0, 10.0])
    def test_law_gamma_invariance_under_the_scale_group(self, c: float) -> None:
        p, phi = _ray()
        gauge = ScaleGauge(p.production_rate, 100.0)
        _require(np.allclose(gauge.apply(c * phi), gauge.apply(phi), rtol=1e-12, atol=0.0), "the section is invariant under ψ ↦ cψ, c ∈ ℝ₊")

    def test_law_the_sign_leg_is_one_case_of_the_section(self) -> None:
        """A NEGATIVE functional value gives a negative displacement and ``apply``
        lands on the positive representative (the dense engines' ℤ/2 fix subsumed)."""
        p, phi = _ray()
        gauge = ScaleGauge(p.production_rate, 100.0)
        _require(gauge.displacement(-phi) < 0.0, "a negative reading gives a negative displacement")
        _require(np.allclose(gauge.apply(-phi), gauge.apply(phi), rtol=1e-12, atol=0.0), "the sign is fixed by the section, not by a second convention")

    def test_negative_a_zero_functional_reading_has_no_section(self) -> None:
        _, phi = _ray()
        gauge = ScaleGauge(lambda state: 0.0, 100.0)
        with pytest.raises(ValueError, match="no section"):
            gauge.displacement(phi)

    def test_negative_two_sections_under_two_functionals_are_distinguishable(self) -> None:
        """Two gauges of one ray under two functionals are DIFFERENT values (the
        TYPE leg) and read DIFFERENT numbers on the same state (the VALUE leg)."""
        p, phi = _ray()
        by_production = ScaleGauge(p.production_rate, 1.0)
        by_absorption = ScaleGauge(p.absorption_rate, 1.0)
        _require(by_production != by_absorption, "the functional is part of the value — the two sections are distinguishable")
        _require(by_production == ScaleGauge(p.production_rate, 1.0), "…and equal to a section built from the same functional and target")
        prod, absn = by_production.functional(phi), by_absorption.functional(phi)
        _require(abs(prod - absn) > 1e-2 * abs(prod), f"the two functionals read different numbers on one state ({prod:.6g} vs {absn:.6g})")
        _require(not np.allclose(by_production.apply(phi), by_absorption.apply(phi)), "so the two representatives differ by a positive scalar the recorded gauge names")
