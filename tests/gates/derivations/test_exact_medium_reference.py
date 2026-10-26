r"""The exact infinite medium as a reference: its factory and its derivation (#405 P2 step 6, R6.6; step 7b.1, R6.8).

Specified by the test-architect (2026-10-03, ``.claude/plans/reference_p2_spec.md``
§1.6 and §1.7b.1). These rows lived in ``tests/gates/reference/test_reference_solution.py``
until the layer-DAG re-filing (Phase 1b, 2026-10-26, ``.claude/plans/layer_dag.md``):
each states a claim about :mod:`orpheus.derivations.common.exact_homogeneous`
(the factory ``exact_infinite_medium_reference`` and its
``ExactInfiniteMediumDerivation``), which the reference layer cannot name. What
the reference layer claims of ANY derivation (a read is the derivation's
evaluation; an uncertified value reads as itself) is gated there on the
test-side ``TableDerivation`` (R6.3, R6.5).

* R6.6: the factory answers the fundamental k question on the infinite medium
  only, and refuses every other question; its reference reads the exact k, the
  exact group fluxes in the declared gauge, ratios and weighted sums exactly;
* R6.8: the derivation evaluates a symbolic weight it cannot certify
  ``Uncertified`` (carrying ``Uncertifiable.approximation``), certifies every
  rational reading, and lets every other refusal of ``Exact`` propagate.
"""

from __future__ import annotations

import math
from fractions import Fraction
from typing import Any

import numpy as np

import pytest

from orpheus.derivations.common import exact_homogeneous
from tests.gates._content_identity_helpers import require
from tests.gates.reference import _step5 as s5
from tests.gates.reference import _step6 as s6
from tests.gates.specification import _fixtures as sf

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/derivations/test_exact_medium_reference.py"


# ═════════════════════════════════════════════════════════════════════════════
# R6.6 — the first instance: the exact infinite medium
# ═════════════════════════════════════════════════════════════════════════════


_MIXTURES = (("fuel-2g-upscatter-n2n", sf.fuel), ("fission-only-2g", sf.fission_only), ("fuel-4g", lambda: sf.fuel(4)))


def _medium(mixture: Any) -> Any:
    from orpheus.data.cells import CellCoefficient, Channel
    from orpheus.numerics.question import Eigen
    from orpheus.specification import InfiniteMediumSpecification

    # The fission gauge is DECLARED (the user's ruling of 2026-10-08: the default production counts the (n,2n)
    # emission too): these rows judge the exact flux of ``ExactInfiniteMedium``, whose own gauge is fission alone at the density 100.
    fission = CellCoefficient.every(Channel.FISSION_EMISSION)
    return InfiniteMediumSpecification(0, mixture, Eigen(fission, gauge=fission))


@pytest.mark.parametrize("make", [m for _, m in _MIXTURES], ids=[n for n, _ in _MIXTURES])
def test_r6_6_the_exact_infinite_medium_reads_the_exact_answer(make: Any) -> None:
    """REFERENCE (exact rationals). For each mixture (2 and 4 groups; anti-#3: a
    1-group k is flux-shape blind), the factory's reference is certified
    ``Valid``; its ``Eigenvalue`` reading contains the exact k_inf of
    ``exact_infinite_medium_of`` (read in the k chart, NEEDS 1) with a bound at
    most one ulp; every group's indicator flux integral contains the exact
    group flux under the module's gauge; and the ratio of the first two groups
    contains the exact ratio, gauge-free; an UNCLAIMED weight (g + 1)/4 reads
    the exact weighted sum (the derivation establishes every flux integral of a
    regionwise-constant weight exactly)."""
    mixture = make()
    exact = exact_homogeneous.exact_infinite_medium_of(mixture)
    ref = exact_homogeneous.exact_infinite_medium_reference(_medium(mixture))
    require(ref.certificate is not None and s5.state_kind(ref.certificate) == "Valid", f"{ref.certificate!r}")
    k = ref.read(s5.eigenvalue())
    if not isinstance(k, s5.ENCLOSURE.Enclosure):
        raise AssertionError(f"k reads {k!r}, not an enclosure")
    require(s6.contains(k, exact.k_inf) and k.bound <= math.ulp(k.value), f"k reading {k!r} against {float(exact.k_inf)!r}")
    groups = len(exact.flux)
    for g in range(groups):
        phi = ref.read(s6.group_flux(g, groups=groups))
        require(s6.contains(phi, exact.flux[g]), f"group {g}: {phi!r} against {float(exact.flux[g])!r}")
    ratio = ref.read(s5.OBSERVABLE.Ratio(s6.group_flux(0, groups=groups), s6.group_flux(1, groups=groups)))
    require(s6.contains(ratio, exact.flux[0] / exact.flux[1]), f"ratio {ratio!r}")
    weights = [(g + 1) / 4 for g in range(groups)]  # dyadic, so the exact weighted sum is the reference
    mixed = ref.read(s5.flux((tuple(weights),)))
    require(s6.contains(mixed, sum((Fraction(w) * f for w, f in zip(weights, exact.flux)), Fraction(0))), f"weighted {mixed!r}")


def test_r6_6_the_factory_refuses_a_question_it_does_not_answer() -> None:
    """The exact family answers the k question on the infinite medium only: a
    geometry specification, and a fixed-source question, are refused."""
    with pytest.raises((TypeError, ValueError)):
        exact_homogeneous.exact_infinite_medium_reference(s5.eigen_slab())
    from orpheus.numerics.mesh_free_function import RegionwiseConstant
    from orpheus.numerics.question import FixedSource
    from orpheus.specification import InfiniteMediumSpecification
    import numpy as np

    source = InfiniteMediumSpecification(0, sf.fuel(), FixedSource(RegionwiseConstant(np.ones((1, 2)))))
    with pytest.raises((TypeError, ValueError)):
        exact_homogeneous.exact_infinite_medium_reference(source)


def test_r6_6_the_factory_refuses_another_eigen_direction() -> None:
    """An eigen question along another direction (every SCATTERING emission,
    the classical c-eigenvalue) is an Eigen question, so the claimed
    eigenvalue's admission does not refuse it: only the factory's own
    k-question check does, naming the k-eigenvalue (the fixed-source row above
    is refused by the claim's admission first, ``[M]`` battery arm E3)."""
    from orpheus.data.cells import CellCoefficient, Channel
    from orpheus.numerics.question import Eigen
    from orpheus.specification import InfiniteMediumSpecification

    c_question = InfiniteMediumSpecification(0, sf.fuel(), Eigen(CellCoefficient.every(Channel.SCATTERING_EMISSION)))
    with pytest.raises(ValueError, match="k-eigenvalue"):
        exact_homogeneous.exact_infinite_medium_reference(c_question)



# ═════════════════════════════════════════════════════════════════════════════
# The step-6 review round (qa and the elegance review, 2026-10-03)
# ═════════════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("variant", ["nearest-mode", "offset-point"])
def test_r6_6_the_factory_answers_only_the_fundamental_at_the_physical_point(variant: str) -> None:
    """The exact family answers ONE question, the fundamental k at the physical
    point; another mode or an offset point is another question, refused (it
    read the unperturbed k, certified Valid, before the review: qa F2,
    elegance 1)."""
    from orpheus.data.cells import CellCoefficient, Channel
    from orpheus.numerics.question import Eigen, Nearest
    from orpheus.specification import InfiniteMediumSpecification

    k = CellCoefficient.every(Channel.FISSION_EMISSION)
    question = (Eigen(k, mode=Nearest(0.3)) if variant == "nearest-mode"
                else Eigen(k, point={CellCoefficient.every(Channel.SCATTERING_EMISSION): 0.5}))
    with pytest.raises(ValueError, match="k-eigenvalue"):
        exact_homogeneous.exact_infinite_medium_reference(InfiniteMediumSpecification(0, sf.fuel(), question))


def _default_medium(mixture: Any) -> Any:
    """The infinite medium's k question with NO declared gauge: the specification writes in the default production."""
    from orpheus.data.cells import CellCoefficient, Channel
    from orpheus.numerics.question import Eigen
    from orpheus.specification import InfiniteMediumSpecification

    return InfiniteMediumSpecification(0, mixture, Eigen(CellCoefficient.every(Channel.FISSION_EMISSION)))


@pytest.mark.rests_on(f"{_HERE}::test_r6_6_the_exact_infinite_medium_reads_the_exact_answer")
def test_r6_6_the_default_gauge_counts_the_n2n_emission() -> None:
    """The user's ruling of 2026-10-08: an undeclared gauge is the default production, what fission and the
    (n,2n) reaction emit. On the fuel (fission, upscatter and a non-zero (n,2n) transfer), the reference
    with no declared gauge reads the weight nu Sigma_f + 2 sum_g' Sigma_2,g->g' (written here from the
    mixture's arrays, the literal 2 the multiplicity) at the density 100, exactly (the reading contains 100);
    the declared fission gauge reads nu Sigma_f at 100 and the production weight above 100 (the (n,2n)
    emission is not zero, the activation). First reds: the reference reading its flux in the fission gauge
    whatever the question declares; the (n,2n) emission counted once per reaction."""
    from orpheus.numerics.mesh_free_function import RegionwiseConstant
    from orpheus.numerics.observable import FluxIntegral

    mixture = sf.fuel()
    nu_sigma_f = np.asarray(mixture.SigP, dtype=float)
    n2n = 2.0 * mixture.Sig2[0].toarray().sum(axis=1)
    production = FluxIntegral(RegionwiseConstant((nu_sigma_f + n2n)[None, :]))
    fission = FluxIntegral(RegionwiseConstant(nu_sigma_f[None, :]))
    require(n2n.min() >= 0.0 and n2n.max() > 0.0, f"activation: the (n,2n) emission {n2n}")
    default = exact_homogeneous.exact_infinite_medium_reference(_default_medium(mixture))
    declared = exact_homogeneous.exact_infinite_medium_reference(_medium(mixture))
    require(s6.contains(default.read(production), Fraction(100)), f"default gauge: {default.read(production)!r}")
    require(s6.contains(declared.read(fission), Fraction(100)), f"fission gauge: {declared.read(fission)!r}")
    require(declared.read(production).value > 100.0 + 1e-6, f"activation: {declared.read(production)!r}")


def test_r6_6_a_non_dyadic_weight_reads_exactly() -> None:
    """Every R6.6 weight above is dyadic, so a float32 rounding of the weight
    left them green (qa F4); 0.1 is not dyadic: the reading must contain the
    exact sum of the DOUBLE 0.1 times the exact flux."""
    mixture = sf.fuel()
    exact = exact_homogeneous.exact_infinite_medium_of(mixture)
    ref = exact_homogeneous.exact_infinite_medium_reference(_medium(mixture))
    weights = (0.1, 0.3)
    reading = ref.read(s5.flux((weights,)))
    truth = sum((Fraction(w) * f for w, f in zip(weights, exact.flux)), Fraction(0))
    require(s6.contains(reading, truth), f"{reading!r} against {float(truth)!r}")


@pytest.mark.parametrize("digits", [15, 3])
def test_r6_6_a_symbolic_float_weight_reads_its_exact_binary_value(digits: int) -> None:
    """A SymPy Float multiplies at its own precision, so the product was rounded
    before Exact certified it: ``Float('0.1', 3)`` read 2.1e-5 relative off with
    a bound of 1e-59 (qa F1). The weight's Float is replaced by the exact binary
    value it holds; the reading must contain that exact product."""
    import sympy

    from orpheus.numerics.mesh_free_function import Symbolic

    mixture = sf.fuel()
    exact = exact_homogeneous.exact_infinite_medium_of(mixture)
    ref = exact_homogeneous.exact_infinite_medium_reference(_medium(mixture))
    weight = sympy.Float("0.1", digits)
    held_rational = sympy.Rational(weight)
    held = Fraction(str(held_rational))
    reading = ref.read(s5.OBSERVABLE.FluxIntegral(Symbolic.of(weight, 0)))
    require(s6.contains(reading, held * exact.flux[0]), f"{reading!r} against {float(held * exact.flux[0])!r}")


def test_r6_6_admission_and_evaluation_share_one_definition_of_constant() -> None:
    """``sin(φ)² + cos(φ)²`` is admitted on the infinite medium (it does not
    depend on φ, :meth:`Symbolic.depends_on`), so it must READ: as the weight 1
    (the elegance review, finding 2: the raw expression was refused by Exact)."""
    import sympy

    from orpheus.numerics.mesh_free_function import Symbolic

    mixture = sf.fuel()
    exact = exact_homogeneous.exact_infinite_medium_of(mixture)
    ref = exact_homogeneous.exact_infinite_medium_reference(_medium(mixture))
    one = sympy.sin(Symbolic.phi) ** 2 + sympy.cos(Symbolic.phi) ** 2
    reading = ref.read(s5.OBSERVABLE.FluxIntegral(Symbolic.of(one, one)))
    require(s6.contains(reading, sum(exact.flux, Fraction(0))), f"{reading!r}")


@pytest.mark.parametrize("table", [np.ones((2, 2)), np.ones((1, 3))], ids=["two-regions", "one-group-too-many"])
def test_r6_review_the_derivation_refuses_a_weight_that_does_not_fit(table: Any) -> None:
    """The exact derivation read ``values[0]`` and ``zip`` truncated, so a
    2-region weight or one group too many returned an Exact value while
    production refused it (the elegance review of step 7a, finding 1). Both
    readers now share ``values_without_position``, which refuses it."""
    from orpheus.numerics.mesh_free_function import RegionwiseConstant

    mixture = sf.fuel()
    ref = exact_homogeneous.exact_infinite_medium_reference(_medium(mixture))
    with pytest.raises(ValueError):
        ref.derivation.evaluate(s5.OBSERVABLE.FluxIntegral(RegionwiseConstant(table)))



# ═════════════════════════════════════════════════════════════════════════════
# R6.8 — the exact family reads uncertified where Exact cannot certify (step 7b.1)
# ═════════════════════════════════════════════════════════════════════════════


def _machin_weight() -> Any:
    from orpheus.numerics.mesh_free_function import Symbolic

    zero = s6.machin_zero()
    return s5.OBSERVABLE.FluxIntegral(Symbolic.of(zero, zero))


def test_r6_8_the_exact_family_reads_an_uncertifiable_value_uncertified() -> None:
    """A symbolic weight whose value SymPy cannot certify (an exactly-zero
    Machin combination) is admitted on the infinite medium (it is constant),
    so it must READ: the exact family evaluates it ``Uncertified`` rather than
    refusing (step 7b.1's API item 3), and ``read`` returns that value, equal
    to what ``evaluate`` returned. Its value is the uncertified float
    evaluation of an exact 0, so it lies within 1e-100 of 0 (``[M]``
    1.94e-121). The certificate stays ``Valid`` (it claims k and the group
    fluxes, all certified). Activation: ``Exact`` itself refuses the value."""
    import sympy

    mixture = sf.fuel()
    ref = exact_homogeneous.exact_infinite_medium_reference(_medium(mixture))
    require(ref.certificate is not None and s5.state_kind(ref.certificate) == "Valid", f"{ref.certificate!r}")
    with pytest.raises(ValueError, match="certified digits"):
        s5.CERTIFICATE.Exact(sympy.srepr(s6.machin_zero()), "activation")
    observable = _machin_weight()
    reading = ref.read(observable)
    require(type(reading) is s6.uncertified_class(), f"read {reading!r}, not an Uncertified")
    require(reading == ref.derivation.evaluate(observable), f"read {reading!r}")
    require(abs(reading.value) <= 1e-100, f"the uncertified value of an exact 0 reads {reading.value!r}")


def test_r6_8_the_exact_family_still_certifies_every_rational_reading() -> None:
    """The control leg: k and a tabulated weight still evaluate ``Exact`` (the
    fallback is reached only where ``Exact`` refuses), so their readings are
    enclosures, and a ratio of the uncertified weight over group 0 reads
    uncertified."""
    mixture = sf.fuel()
    ref = exact_homogeneous.exact_infinite_medium_reference(_medium(mixture))
    for observable in (s5.eigenvalue(), s6.group_flux(0), s5.flux(((0.1, 0.3),))):
        evaluation = ref.derivation.evaluate(observable)
        require(type(evaluation) is s5.CERTIFICATE.Exact, f"{observable!r} evaluates {evaluation!r}")
        require(type(ref.read(observable)) is s5.ENCLOSURE.Enclosure, f"{observable!r} does not read enclosed")
    ratio = ref.read(s5.OBSERVABLE.Ratio(_machin_weight(), s6.group_flux(0)))
    require(type(ratio) is s6.uncertified_class(), f"the ratio reads {ratio!r}")


def test_r6_8_the_exact_familys_point_arm_refuses() -> None:
    """``read`` admits first, so a point value never reaches the derivation
    (R6.2); called directly, the derivation's unreachable ``PointValue`` arm
    raises ``ValueError`` naming the missing position (it raised
    ``NotCertified``, retired by step 7b.1), never an ``Uncertified``."""
    mixture = sf.fuel()
    ref = exact_homogeneous.exact_infinite_medium_reference(_medium(mixture))
    with pytest.raises(ValueError, match="position"):
        ref.derivation.evaluate(s5.point(0.5, 0))



# ── the review round (qa F2 and the narrowed catch, 2026-10-03) ─────────────


def _offset_weight(groups: tuple[bool, bool], exponent: int = 130) -> Any:
    """The flux integral whose weight is the Machin zero plus 10^-exponent on the flagged groups, 0 elsewhere."""
    import sympy

    from orpheus.numerics.mesh_free_function import Symbolic

    w = s6.machin_zero() + sympy.Integer(10) ** -exponent
    return s5.OBSERVABLE.FluxIntegral(Symbolic.of(*(w if on else 0 for on in groups)))


def _offset_truth(groups: tuple[bool, bool], exponent: int = 130) -> Fraction:
    exact = exact_homogeneous.exact_infinite_medium_of(sf.fuel())
    return Fraction(1, 10**exponent) * sum((f for on, f in zip(groups, exact.flux) if on), Fraction(0))


def test_r6_8_an_uncertified_value_is_the_value_of_its_expression() -> None:
    """qa F2, closing the declared-blind arm X2 (the fallback read as 0.0): the
    weight ``z + 10^-130`` on group 0 (z the Machin zero) is uncertifiable
    (strict ``evalf`` cancels; ``Exact`` raises ``Uncertifiable``, the
    activation) and NONZERO, so its uncertified reading must carry the value
    of the expression: within 1e-15 relative of the exact ``10^-130 φ_0``
    (``[M]`` 7.9e-17)."""
    import sympy

    groups = (True, False)
    ref = exact_homogeneous.exact_infinite_medium_reference(_medium(sf.fuel()))
    observable = _offset_weight(groups)
    with pytest.raises(s5.CERTIFICATE.Uncertifiable):
        s5.CERTIFICATE.Exact(sympy.srepr((s6.machin_zero() + sympy.Integer(10) ** -130) * 2), "activation")
    reading = ref.read(observable)
    require(type(reading) is s6.uncertified_class(), f"read {reading!r}, not an Uncertified")
    truth = _offset_truth(groups)
    relative = abs(Fraction(reading.value) - truth) / truth
    require(relative <= Fraction(1, 10**15), f"the uncertified value {reading.value!r} is {float(relative):.3g} relative off")


def _machin_inverse() -> Any:
    return 1 / s6.machin_zero()


def _acos_two() -> Any:
    import sympy

    return sympy.acos(2)


_DEFECTS = (("one-over-the-machin-zero", _machin_inverse), ("acos-of-two", _acos_two))


@pytest.mark.parametrize("make", [m for _, m in _DEFECTS], ids=[n for n, _ in _DEFECTS])
def test_r6_8_a_defective_expression_is_refused_never_read_uncertified(make: Any) -> None:
    """The review round narrowed the fallback to ``Uncertifiable`` alone (the
    one refusal of ``Exact`` meaning "no bound can be derived"). Every other
    refusal is a defect of the expression and propagates: ``1/z`` (not
    provably finite; it read ``Uncertified(7e131)`` before) and ``acos(2)``
    (not real) raise a ``ValueError`` that is NOT ``Uncertifiable``, and
    nothing is read."""
    from orpheus.numerics.mesh_free_function import Symbolic

    e = make()
    ref = exact_homogeneous.exact_infinite_medium_reference(_medium(sf.fuel()))
    with pytest.raises(ValueError, match="finite real constant") as refusal:
        ref.read(s5.OBSERVABLE.FluxIntegral(Symbolic.of(e, e)))
    require(not isinstance(refusal.value, s5.CERTIFICATE.Uncertifiable), f"refused as uncertifiable: {refusal.value}")
