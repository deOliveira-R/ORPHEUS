r"""Gates of the typed realized boundary operators ``orpheus/sn/boundary/realized_operators.py`` (wave 2, slice S2b).

Spec: ``.claude/plans/wave2_boundary_carve.md`` section 2.6 (the types), section 3 (laws L-A, L-E) and section 6.2 (S2b
gates 14 to 18 and 20). PLACEMENT: the spec names ``tests/gates/sn/operators/``; this file lives in ``tests/gates/sn/architecture/``
because the ``operators`` directory's ``conftest.py`` imports ``tests/gates/sn/_test_helpers.py``, which still imports the
retired ``BC``: a gate there has a conftest abort for its base red, which hides its own reason and aborts every file of the
same pytest invocation. ``architecture/`` has no conftest, so each gate fails at the line it names. The module under test
does not exist at the base tip: every row's first red is the missing module (the
accessor ``ro()`` of the helper imports it at first use), except the mass census, which reads only numerics and is green
today (it is the PRECONDITION of the divisor rule, run first as the spec asks).

What the types are
==================

Realization is a functor on the response network: the empty network goes to ``ZeroMap``, ``s * SpecularChain`` to
``ScaledOrdinatePermutation`` (``s P_pi``), ``d * DiffuseChain`` to ``RankOneReemission`` (``d (E o Rt)``), a sum to
``RealizedSum``. Each fact the consumers read (``return_fraction``, ``permutes_ordinates``, ``domain_face``,
``is_vacuum``) is a PROPERTY derived from the constructor's parameters ``(alpha, ..., face, trace)``, and the operator
``inner`` is derived too, from the spaces read off the trace ladder: no fact is a field beside the operator, so a lying
``alpha`` is unspellable.

The ladder (``rests_on`` edges are declared on the tests)::

    test_the_mass_census_...                    precondition (numerics only; green at the base tip)
    test_the_constructors_take_exactly_...      foundation: the signatures
    test_the_typed_operators_answer_...         foundation: the facts of each type
        ^
    test_return_fraction_is_the_measured_...    edge: the fact equals the measured J-/J+ of the built operator
        ^
    test_reciprocity_...                        interior: L-A (i) on S16
    test_the_left_eigenvector_rider_...         interior: L-A (iv), with an independent witness
    test_the_algebra_object_delegates           composition: apply, transpose, spaces, identity
    test_the_declared_and_the_realized_composition_meet   composition: the chains' tiers meet the spaces

Mutations for the reviewers (each on a COPY of ``realized_operators.py`` built by transforming ``inspect.getsource``;
the verdict is per arm, and the measured red sets are in ``scratch/wave2_s2b/battery/``):

* M-free-norm: ``RankOneReemission`` takes a free ``norm`` argument (the lying-parameter shape DRAFT 2 left open):
  reddens ``test_the_constructors_take_exactly_their_parameters``.
* M-alpha-squared: the builders scale by ``alpha**2``: reddens the ratio gate and the rider (the one-vector ratio is
  the weaker witness, the array identity of the rider is the stronger).
* M-outflow-divisor: ``RankOneReemission`` divides by the OUTFLOW mass: reddens only on the asymmetric-trace rows
  (``M_+ != M_-``), the first red of law L-A (iv): on every shipped quadrature ``M_+ == M_-`` bitwise (the census), so
  no symmetric row can see it.
* M-sum-fsum: ``RealizedSum.return_fraction`` is a left-to-right addition of floats (``reduce(add, ...)``; the built-in
  ``sum`` is compensated since Python 3.12): reddens the bit-exact row on the weights ``(0.1, 0.2, 0.3)``.
* M-reciprocity: a Lambertian weighted by ``c |mu|`` (``R ~ mu'^2``): reddens the reciprocity gate.
* M-stored-alpha: ``return_fraction`` becomes a stored field set by hand: reddens ``test_a_fact_cannot_be_set_by_hand``.
"""

from __future__ import annotations

import dataclasses
import functools
import inspect
import math
import operator

import numpy as np
import pytest

from tests.gates.sn.architecture import _realized_boundary_helpers as H

_THIS = "tests/gates/sn/architecture/test_realized_boundary.py"

#: Slab traces of the four quadratures of the B3.4c measurement, per face.
_QUADS = ("gl8", "product(2,4)", "ls6", "lebedev17")
_FACES = ("xmin", "xmax")


def _build(kind: str, quad_name: str, face: str, *, alpha: float = 0.7):
    """One typed operator of ``kind`` on ``face`` of the slab trace of ``quad_name``; returns ``(op, quad, trace)``."""
    quad = H.QUADRATURES[quad_name]()
    trace = H.trace_of(quad)
    R = H.ro()
    if kind == "zero":
        op = R.ZeroMap(face=face, trace=trace)
    elif kind == "specular":
        op = H.scaled_ordinate_permutation(alpha, H.mirror_pairing(quad, face), face, face, trace)
    elif kind == "translation":
        partner = "xmax" if face == "xmin" else "xmin"
        op = H.scaled_ordinate_permutation(1.0, H.identity_pairing(quad), face, partner, trace)
    elif kind == "diffuse":
        op = R.RankOneReemission(alpha=alpha, face=face, trace=trace)
    elif kind == "sum":
        a = H.scaled_ordinate_permutation(0.3, H.mirror_pairing(quad, face), face, face, trace)
        b = R.RankOneReemission(alpha=0.5, face=face, trace=trace)
        op = R.RealizedSum(terms=(a, b))
    else:  # pragma: no cover
        raise ValueError(kind)
    return op, quad, trace


# ═════════════════════════════════════════════════════════════════════════════
# The precondition: the divisor rule (section 2.6, "The divisor rule"; R-11)
# ═════════════════════════════════════════════════════════════════════════════

#: The registry's population, listed HERE (a gate over a list owes its input count, X2). Gauss-Legendre, level-symmetric,
#: Lebedev and product rules, the four families of ``quadrature_registry``.
_CENSUS = (
    [(f"gl{n}", lambda n=n: H.quadrature_class().gauss_legendre(n_ordinates=n)) for n in (2, 4, 8, 16, 32, 64)]
    + [(f"ls{s}", lambda s=s: H.quadrature_class().level_symmetric(sn_order=s)) for s in (2, 4, 6, 8, 10, 12, 14, 16)]
    + [(f"lebedev{o}", lambda o=o: H.quadrature_class().lebedev(o)) for o in (3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31)]
    + [(f"product({a},{b})", lambda a=a, b=b: H.quadrature_class().product(n_mu=a, n_phi=b)) for a, b in ((2, 4), (4, 4), (4, 8), (8, 8), (8, 16), (16, 16))]
)

#: Pairs ``(quadrature, face)`` whose inflow and outflow masses differ in a bit. Section 8 row R-11 lists them by
#: quadrature BEFORE the commit; empty today (``[M]`` 70 of 70 pairs equal, 2026-10-23).
DECLARED_MASS_EXCEPTIONS: frozenset[tuple[str, str]] = frozenset()


@pytest.mark.foundation
def test_the_mass_census_inflow_equals_outflow_on_every_registry_rule() -> None:
    r"""``M_+ == M_-`` bitwise (the plain ``np.sum`` of ``w |Omega.n|`` over each half-trace) for every rule of the
    registry and both faces of a slab, or the pair is listed in :data:`DECLARED_MASS_EXCEPTIONS` (section 8, R-11).

    This is the divisor rule's precondition: the Lambertian image divides by the INFLOW mass, which equals the shipped
    OUTFLOW divisor exactly where these two sums are bit-equal, so the carry is bit-identical there. Population: 35
    rules times 2 faces, printed in the failure message. The control is the asymmetric trace of law L-A (iv): a rule
    whose outflow weights are scaled reads UNEQUAL, so a census that compared nothing could not pass it.

    Red input: none today (the census is green: a precondition, not a defect). It reddens on a registry rule whose
    half-trace masses differ, which is exactly the case section 8 row R-11 must then list.
    """
    unequal = []
    for name, build in _CENSUS:
        quad = build()
        trace = H.trace_of(quad)
        for face in _FACES:
            m_out, m_in = H.trace_weights(trace, face)
            if float(np.sum(m_out)) != float(np.sum(m_in)):
                unequal.append((name, face))
    population = len(_CENSUS) * len(_FACES)
    assert population == 70, f"the census population changed: {population} pairs"
    assert frozenset(unequal) == DECLARED_MASS_EXCEPTIONS, f"{len(unequal)} of {population} (quadrature, face) pairs have M_+ != M_-: {unequal}; list them in R-11"
    proxy, trace = H.asymmetric_trace()
    m_out, m_in = H.trace_weights(trace, "xmax")
    assert float(np.sum(m_out)) != float(np.sum(m_in)), "the control: the scaled trace must read M_+ != M_-"


# ═════════════════════════════════════════════════════════════════════════════
# Gate 14: the facts are properties of the types
# ═════════════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_the_constructors_take_exactly_their_parameters() -> None:
    r"""The parameters ARE the type: ``ZeroMap(face, trace)``, ``ScaledOrdinatePermutation(alpha, permutation, face,
    partner, trace)``, ``RankOneReemission(alpha, face, trace)``, ``RealizedSum(terms)``. No ``inner`` argument, and no
    ``norm``, ``contraction_weights`` or ``inflow_extent`` (all three are READ from the trace's half-spaces), so the
    lying-parameter shape DRAFT 2 left open is unspellable.

    Red input: the module does not exist at the base tip. Mutation M-free-norm: a ``RankOneReemission`` that takes a
    free ``norm`` reddens it.
    """
    R = H.ro()
    expected = {
        "ZeroMap": {"face", "trace"},
        "ScaledOrdinatePermutation": {"alpha", "permutation", "face", "partner", "trace"},
        "RankOneReemission": {"alpha", "face", "trace"},
        "RealizedSum": {"terms"},
    }
    for name, parameters in expected.items():
        got = set(inspect.signature(getattr(R, name)).parameters)
        assert got == parameters, f"{name} takes {sorted(got)}, the spec says {sorted(parameters)}"
    for forbidden in ("inner", "norm", "contraction_weights", "inflow_extent"):
        for name in expected:
            assert forbidden not in inspect.signature(getattr(R, name)).parameters, f"{name} takes a {forbidden!r} argument"


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_constructors_take_exactly_their_parameters")
@pytest.mark.parametrize("quad_name", ("gl8", "lebedev17"))
@pytest.mark.parametrize("face", _FACES)
def test_the_typed_operators_answer_from_their_type(quad_name: str, face: str) -> None:
    r"""Each type answers the five questions the consumers ask of a realized face, from its own parameters:

    * ``ZeroMap``: ``is_vacuum`` True, ``return_fraction`` 0.0, ``permutes_ordinates`` False, ``domain_face == face``;
    * ``ScaledOrdinatePermutation`` of a mirror: ``permutes_ordinates`` True, ``return_fraction == alpha``,
      ``domain_face == partner``; of the identity with a partner (the translation wrap): ``permutes_ordinates`` False;
    * ``RankOneReemission``: ``permutes_ordinates`` False, ``return_fraction == alpha``, ``domain_face == face``;
    * ``RealizedSum``: ``return_fraction`` is the ``math.fsum`` of its terms' (bit-exact on a sum whose left-to-right
      float sum differs), ``permutes_ordinates`` False, ``is_vacuum`` False.

    Red input: the module does not exist at the base tip. Mutation: ``return_fraction`` of the sum as a plain ``sum``
    (M-sum-fsum, rows below); a mirror reporting ``permutes_ordinates`` False.
    """
    zero, _, _ = _build("zero", quad_name, face)
    assert zero.is_vacuum is True and zero.return_fraction == 0.0 and zero.permutes_ordinates is False and zero.domain_face == face

    mirror, _, _ = _build("specular", quad_name, face, alpha=0.7)
    assert mirror.is_vacuum is False and mirror.permutes_ordinates is True
    assert mirror.return_fraction == 0.7 and mirror.domain_face == face

    wrap, _, _ = _build("translation", quad_name, face)
    partner = "xmax" if face == "xmin" else "xmin"
    assert wrap.permutes_ordinates is False and wrap.return_fraction == 1.0 and wrap.domain_face == partner

    diffuse, _, _ = _build("diffuse", quad_name, face, alpha=0.3)
    assert diffuse.permutes_ordinates is False and diffuse.return_fraction == 0.3 and diffuse.domain_face == face and diffuse.is_vacuum is False

    total, quad, trace = _build("sum", quad_name, face)
    assert total.permutes_ordinates is False and total.is_vacuum is False and total.domain_face == face
    assert total.return_fraction == math.fsum([0.3, 0.5])

    # the bit-exact row: three weights whose left-to-right float sum is NOT their correctly rounded sum
    R = H.ro()
    parts = tuple(R.RankOneReemission(alpha=a, face=face, trace=trace) for a in (0.1, 0.2, 0.3))
    # (``sum`` of floats is compensated since Python 3.12, so the plain left-to-right addition is spelled by ``reduce``)
    assert math.fsum([0.1, 0.2, 0.3]) != functools.reduce(operator.add, [0.1, 0.2, 0.3]), "the fixture must separate fsum from left-to-right addition"
    assert R.RealizedSum(terms=parts).return_fraction == math.fsum([0.1, 0.2, 0.3])


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_constructors_take_exactly_their_parameters")
def test_alpha_is_a_scale_in_the_half_open_unit_interval_and_the_sum_is_parsed() -> None:
    r"""``alpha`` lies in ``(0, 1]``: ``0`` is the ``ZeroMap`` (refused here), ``> 1`` is super-Markov, NaN and negative
    are not scales. ``RealizedSum`` is parsed: fewer than two terms, a ``ZeroMap`` term (a zero term is the empty
    network's image and is never summed) and terms on different faces are refused at construction. The controls
    construct (``alpha`` 1.0 and ``1e-12``; two terms on one face), so a type that refused everything cannot pass.

    Red input: the module does not exist at the base tip.
    """
    R = H.ro()
    quad = H.QUADRATURES["gl8"]()
    trace = H.trace_of(quad)
    pairing = H.mirror_pairing(quad, "xmax")
    for bad in (0.0, -0.1, 1.0000000000000002, 2.0, float("nan"), float("inf")):
        with pytest.raises(ValueError):
            H.scaled_ordinate_permutation(bad, pairing, "xmax", "xmax", trace)
        with pytest.raises(ValueError):
            R.RankOneReemission(alpha=bad, face="xmax", trace=trace)
    for good in (1.0, 1e-12):
        H.scaled_ordinate_permutation(good, pairing, "xmax", "xmax", trace)
        R.RankOneReemission(alpha=good, face="xmax", trace=trace)
    a = R.RankOneReemission(alpha=0.2, face="xmax", trace=trace)
    b = R.RankOneReemission(alpha=0.3, face="xmax", trace=trace)
    c = R.RankOneReemission(alpha=0.3, face="xmin", trace=trace)
    R.RealizedSum(terms=(a, b))
    with pytest.raises(ValueError):
        R.RealizedSum(terms=(a,))
    with pytest.raises(ValueError):
        R.RealizedSum(terms=(a, R.ZeroMap(face="xmax", trace=trace)))
    with pytest.raises(ValueError):
        R.RealizedSum(terms=(a, c))


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_typed_operators_answer_from_their_type")
def test_a_fact_cannot_be_set_by_hand() -> None:
    r"""The facts are derived, never stored: assigning ``alpha`` or a fact raises, and ``dataclasses.replace(op,
    alpha=...)`` REBUILDS the operator, so the rebuilt ``return_fraction`` and the rebuilt matrix both carry the new
    scale (the matrix is the array witness that ``inner`` is derived from ``alpha``, not stored beside it).

    Red input: the module does not exist at the base tip. Mutation M-stored-alpha: a ``return_fraction`` stored as a
    settable field.
    """
    R = H.ro()
    op, quad, trace = _build("specular", "gl8", "xmax", alpha=0.7)
    for name, value in (("alpha", 0.2), ("return_fraction", 0.2), ("permutes_ordinates", False), ("domain_face", "xmin"), ("inner", None)):
        with pytest.raises(AttributeError):
            setattr(op, name, value)
    n = H.n_outflow(trace, "xmax")
    before = H.dense(op, n)
    rebuilt = dataclasses.replace(op, alpha=0.3)
    after = H.dense(rebuilt, n)
    assert rebuilt.return_fraction == 0.3 and after.max() == 0.3 and before.max() == 0.7
    assert isinstance(rebuilt, R.ScaledOrdinatePermutation)
    diffuse, _, _ = _build("diffuse", "gl8", "xmax", alpha=0.7)
    scaled = dataclasses.replace(diffuse, alpha=0.35)
    assert scaled.return_fraction == 0.35
    np.testing.assert_allclose(H.dense(scaled, n), 0.5 * H.dense(diffuse, n), rtol=4 * np.finfo(float).eps, atol=0.0)


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_typed_operators_answer_from_their_type")
@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason=(
        "rider (b) of section 2.6, UNRATIFIED (NEEDS of the S2b gates): ordinate_pairing returns a bare Permutation, so an "
        "uncertified quadrature.ordinate_permutation(motion) cannot be told from a certified one. Lands with a frozen "
        "OrdinatePairing carrier (and a translation arm in ordinate_pairing); the strict xfail then flips and this marker goes."
    ),
)
def test_an_uncertified_permutation_cannot_be_realized() -> None:
    r"""Pattern 4 on the rider: the permutation a ``ScaledOrdinatePermutation`` realizes must be one ``ordinate_pairing``
    certified. The raw ``quadrature.ordinate_permutation(motion)`` of the same rule (what the realizer's own builder used
    to take) is REFUSED at construction (``TypeError`` or ``ValueError``), the certified one is accepted (the positive leg).

    NOT YET LANDED: today both are a bare ``Permutation``, indistinguishable, so the second leg does not raise. The
    paired RECORD row below is green today and reddens at the carve that lands the carrier.
    """
    from orpheus.geometry.deck import mirror_motion

    quad = H.QUADRATURES["gl8"]()
    trace = H.trace_of(quad)
    raw = quad.ordinate_permutation(mirror_motion("x", 1))
    H.scaled_ordinate_permutation(1.0, H.mirror_pairing(quad, "xmax"), "xmax", "xmax", trace)
    refused = False
    try:
        H.scaled_ordinate_permutation(1.0, raw, "xmax", "xmax", trace)
    except (TypeError, ValueError):
        refused = True
    assert refused, "an uncertified permutation was realized"


@pytest.mark.foundation
def test_record_ordinate_pairing_returns_a_bare_permutation() -> None:
    r"""RECORD (what the code does today, designed to redden at the carve): ``ordinate_pairing`` returns the bare
    ``Permutation`` of ``Quadrature.ordinate_permutation``, equal to the raw matcher result. The day it returns a carrier
    (rider (b) ratified) this row reddens and ``test_an_uncertified_permutation_cannot_be_realized`` flips.

    Green today (numerics only); a RECORD row is not coverage of the rider.
    """
    from orpheus.geometry.deck import mirror_motion
    from orpheus.geometry.transformation import Permutation

    quad = H.QUADRATURES["gl8"]()
    certified = H.mirror_pairing(quad, "xmax")
    assert type(certified) is Permutation
    np.testing.assert_array_equal(certified.indices, quad.ordinate_permutation(mirror_motion("x", 1)).indices)


# ═════════════════════════════════════════════════════════════════════════════
# Gate 15: return_fraction is the measured partial-current ratio (ERR-094's collapse trigger)
# ═════════════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.catches("ERR-094")
@pytest.mark.rests_on(f"{_THIS}::test_the_typed_operators_answer_from_their_type", f"{_THIS}::test_the_mass_census_inflow_equals_outflow_on_every_registry_rule")
@pytest.mark.parametrize("kind", ("zero", "specular", "translation", "diffuse", "sum"))
@pytest.mark.parametrize("quad_name", _QUADS)
@pytest.mark.parametrize("face", _FACES)
def test_return_fraction_is_the_measured_partial_current_ratio(kind: str, quad_name: str, face: str) -> None:
    r"""For each type built over a synthetic positive outflow, ``J^- / J^+`` of the built operator's image equals
    ``return_fraction`` to ``1e-12``, with ``J^+ = sum(w|mu| psi)`` over ``Gamma_+`` and ``J^-`` the same over
    ``Gamma_-``, both in the ladder's metric (read from the half-spaces).

    This is the law-L-A witness that collapses ERR-094 structurally: the day an angle-dependent reflectivity ships, its
    realizer computes ``return_fraction`` on the actual trace and the type carries it. A ``ZeroMap`` returns nothing
    (``J^- = 0``). It is the WEAKER witness of the pair (one vector): the array identity of the rider gate is the
    stronger. Red input: the module does not exist at the base tip. Mutation M-alpha-squared: a builder scaling by
    ``alpha**2`` reads ``0.49`` against ``0.7``.
    """
    op, quad, trace = _build(kind, quad_name, face, alpha=0.7)
    n_out = H.n_outflow(trace, face)
    psi = H.positive_outflow(n_out)
    if kind == "translation":
        domain_face = op.domain_face
        m_out, _ = H.trace_weights(trace, domain_face)
        _, m_in = H.trace_weights(trace, face)
        psi = H.positive_outflow(H.n_outflow(trace, domain_face))
    else:
        m_out, m_in = H.trace_weights(trace, face)
    j_plus = float(np.sum(m_out * psi))
    j_minus = float(np.sum(m_in * np.asarray(op.apply(psi))))
    assert abs(j_minus / j_plus - op.return_fraction) <= 1e-12, f"{kind} on {quad_name}/{face}: J-/J+ = {j_minus / j_plus!r}, return_fraction = {op.return_fraction!r}"


# ═════════════════════════════════════════════════════════════════════════════
# Gate 16: reciprocity (law L-A (i))
# ═════════════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.rests_on(f"{_THIS}::test_return_fraction_is_the_measured_partial_current_ratio")
@pytest.mark.parametrize("face", _FACES)
def test_reciprocity(face: str) -> None:
    r"""Law L-A (i): :math:`M_+^{-1} B^{\mathsf T} M_- = P_{\rm rev}\, B\, P_{\rm rev}` for the specular and the
    Lambertian realized operators on an S16 slab face (``M_\pm`` the half-space metrics read from the ladder, ``P_rev``
    the :math:`\mu \to -\mu` reversal built by matching direction cosines, never from the realizer). ``array_equal`` for
    the specular image (its entries are ``alpha`` and the paired weights are bit-equal), ``1e-15`` for the Lambertian.

    The control makes the gate able to fail: the same check applied to the non-reciprocal matrix
    ``R(mu' -> mu) ~ mu'^2`` (the Lambertian built with ``c |mu|``) reads a defect orders above the bound.

    Red input: the module does not exist at the base tip; mutation M-reciprocity.
    """
    quad = H.QUADRATURES["gl16"]()
    trace = H.trace_of(quad)
    n = H.n_outflow(trace, face)
    spec, _, _ = _build("specular", "gl16", face, alpha=0.7)
    diff, _, _ = _build("diffuse", "gl16", face, alpha=0.3)
    assert H.reciprocity_defect(H.dense(spec, n), quad, trace, face) == 0.0
    assert H.reciprocity_defect(H.dense(diff, n), quad, trace, face) <= 1e-15

    m_out, m_in = H.trace_weights(trace, face)
    mu_out = np.abs(np.asarray(quad.mean_axis_cosine(0), dtype=float))[np.asarray(trace.outflow_indices_for_face(face))]
    wrong = np.ones((m_in.size, 1)) * (m_out * mu_out)[None, :] * 0.3 / float(np.sum(m_in))
    assert H.reciprocity_defect(wrong, quad, trace, face) > 1e-3, "the control: a non-reciprocal response must read a defect"


# ═════════════════════════════════════════════════════════════════════════════
# Gate 17: the algebra object delegates; equality and hash are over the parameters
# ═════════════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_typed_operators_answer_from_their_type")
@pytest.mark.parametrize("kind", ("zero", "specular", "translation", "diffuse", "sum"))
def test_the_algebra_object_delegates(kind: str) -> None:
    r"""``apply``, ``apply_transpose``, ``domain``, ``codomain`` and ``block_role`` of each realized type equal those of
    the operator it builds (the shim contract of ``bound_compat.py:109-250``), and that operator carries
    ``BlockRole.BOUNDARY``. Equality and hash are over the PARAMETERS (``type``, ``face``, ``partner``, ``alpha``, the
    permutation's index bytes), not over a stored operator:

    * two types with equal ``alpha`` and different ``partner`` compare unequal;
    * a type rebuilt from equal parameters (a new trace of the same quadrature) is equal, with an equal hash;
    * a type compared with its builder's bare operator is unequal, in both directions;
    * the old shim's string ``__eq__`` arm (#553) is NOT carried: ``op == "reflective"`` and ``op == "vacuum"`` are False.

    Red input: the module does not exist at the base tip. Mutation: a ``__eq__`` that keeps the string arm.
    """
    from orpheus.numerics.operator import BlockRole

    op, quad, trace = _build(kind, "gl8", "xmax")
    n_out = H.n_outflow(trace, "xmax") if kind != "translation" else H.n_outflow(trace, op.domain_face)
    x = H.positive_outflow(n_out)
    np.testing.assert_array_equal(op.apply(x), op.inner.apply(x))
    y = np.asarray(op.apply(x))
    np.testing.assert_array_equal(op.apply_transpose(y), op.inner.apply_transpose(y))
    assert op.domain == op.inner.domain and op.codomain == op.inner.codomain
    assert op.block_role == op.inner.block_role == BlockRole.BOUNDARY

    twin, _, _ = _build(kind, "gl8", "xmax")
    assert twin == op and hash(twin) == hash(op), "equal parameters must compare and hash equal"
    assert op != op.inner and op.inner != op, "a type is not its builder's bare operator"
    for text in ("reflective", "vacuum", "white", "specular", "albedo"):
        assert (op == text) is False, f"the retired string arm answers {text!r} (#553)"
    assert len({op, twin}) == 1


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_algebra_object_delegates")
def test_equality_reads_the_partner_the_alpha_and_the_permutation() -> None:
    r"""The parameter list of the equality is each of its members: ``partner``, ``alpha``, the permutation's index bytes
    and the face each separate two values that agree on the others.

    Red input: the module does not exist at the base tip. Mutation: an equality that drops ``partner``.
    """
    quad = H.QUADRATURES["gl8"]()
    trace = H.trace_of(quad)
    ident = H.identity_pairing(quad)
    mirror = H.mirror_pairing(quad, "xmin")
    def make(**kw):
        values = {"alpha": 1.0, "permutation": ident, "face": "xmin", "partner": "xmax", **kw}
        return H.scaled_ordinate_permutation(values["alpha"], values["permutation"], values["face"], values["partner"], trace)

    base = make()
    assert base == make()
    assert base != make(partner="xmin")
    assert base != make(alpha=0.5)
    assert base != make(permutation=mirror)
    assert base != make(face="xmax", partner="xmin")


# ═════════════════════════════════════════════════════════════════════════════
# Gate 18: the left-eigenvector rider (law L-A (iv))
# ═════════════════════════════════════════════════════════════════════════════


def _split_pair_witness(trace, face: str, alpha: float, psi: np.ndarray) -> np.ndarray:
    r""":math:`\alpha\,(E \circ R_t)\psi` built from ``AxisRetractionOperator`` and ``AxisSectionOperator`` over a
    TEST-LOCAL nodal axis carrying the measure :math:`w|\mu|` (the trace half-spaces carry no axes, so the numerics pair
    cannot be read off them: section 2.6, decision 57). A trivial one-group energy axis supplies the second axis a
    marginal needs. The divisors are ``math.fsum`` here and a plain ``np.sum`` in the realized operator, which differ by
    a rounding on lebedev 17 and GL 16 (``[M]``)."""
    from orpheus.numerics.axis import Axis, EnergyAxis
    from orpheus.numerics.axis_marginal import AxisRetractionOperator, AxisSectionOperator
    from orpheus.numerics.index_set import IndexSet
    from orpheus.numerics.measure import DiscreteMeasure
    from orpheus.numerics.space import FunctionSpace

    def space_of(weights: np.ndarray) -> FunctionSpace:
        points = IndexSet(label="angular", n=weights.size)
        measure = DiscreteMeasure(nodes=np.arange(weights.size, dtype=float), weights=weights, support=points)
        return FunctionSpace.of_axes(Axis.of_measure(measure, "angular"), EnergyAxis.synthetic(1))

    m_out, m_in = H.trace_weights(trace, face)
    retraction = AxisRetractionOperator.of(space_of(m_out), "angular")
    section = AxisSectionOperator.of(space_of(m_in), "angular")
    return alpha * np.asarray(section.apply(retraction.apply(psi[:, None])))[:, 0]


@pytest.mark.l1
@pytest.mark.rests_on(f"{_THIS}::test_return_fraction_is_the_measured_partial_current_ratio")
@pytest.mark.parametrize("quad_name", _QUADS)
@pytest.mark.parametrize("face", _FACES)
def test_the_left_eigenvector_rider(quad_name: str, face: str) -> None:
    r"""Law L-A (iv): ``op.apply_transpose(w|mu| on Gamma_-) == rho * (w|mu| on Gamma_+)`` with ``rho`` the return
    fraction, for the specular image (``array_equal``) and the Lambertian image (within the split pair's residual,
    ``ceil(N/3)`` ulp per entry, ``N`` the number of outflow ordinates; ``[M]`` 0 ulp on all eight rows, the bound is
    the contract of ``axis_marginal.py``). The vector ``w|mu|`` on ``Gamma_-`` is the ladder's metric: the left
    eigenvector of the response, the property that makes the return fraction the eigenvalue.

    The INDEPENDENT witness of the Lambertian is built without the realizer's divisor: the numerics pair over a
    test-local axis with the measure ``w|mu|`` (``_split_pair_witness``), within ``ceil(N/3) + 1`` ulp (the pair's own
    residual plus one rounding for the plain-versus-``fsum`` divisor; ``[M]`` at most 1).

    Red input: the module does not exist at the base tip. Mutations: M-alpha-squared (the specular array identity
    reads ``alpha**2 m`` against ``alpha m``); the divisor mutation M-outflow-divisor is NOT visible here (the masses
    are bit-equal on every shipped rule) and reddens the asymmetric-trace row below.
    """
    spec, quad, trace = _build("specular", quad_name, face, alpha=0.7)
    diff, _, _ = _build("diffuse", quad_name, face, alpha=0.3)
    m_out, m_in = H.trace_weights(trace, face)
    np.testing.assert_array_equal(spec.apply_transpose(m_in), 0.7 * m_out)
    bound = math.ceil(m_out.size / 3)
    got = np.asarray(diff.apply_transpose(m_in))
    ulps = np.max(np.abs(got - 0.3 * m_out) / np.spacing(np.abs(0.3 * m_out)))
    assert ulps <= bound, f"rider of the Lambertian: {ulps} ulp > {bound}"

    psi = H.positive_outflow(m_out.size)
    witness = _split_pair_witness(trace, face, 0.3, psi)
    image = np.asarray(diff.apply(psi))
    witness_ulps = np.max(np.abs(image - witness) / np.spacing(np.abs(witness)))
    assert witness_ulps <= bound + 1, f"the realized Lambertian differs from alpha (E o Rt) of the numerics pair by {witness_ulps} ulp > {bound + 1}"


@pytest.mark.l1
@pytest.mark.rests_on(f"{_THIS}::test_the_left_eigenvector_rider")
def test_the_rider_on_a_trace_whose_masses_differ() -> None:
    r"""The first red of law L-A (iv): a delegating quadrature proxy whose OUTFLOW-side weights of ``xmax`` are scaled by
    1.1 builds a trace with ``M_+ = 0.5563404348779897 != M_- = 0.5057640317072634`` (``[M]`` ``_d3_trace_axes.py``;
    the unperturbed trace has ``M_+ == M_-``). On it the Lambertian built to the INFLOW-mass rule still has the
    eigenvalue its ``alpha`` declares, ``apply_transpose(w|mu| on Gamma_-) == alpha (w|mu| on Gamma_+)`` to
    ``ceil(N/3)`` ulp and ``J^-/J^+ = alpha``; an image dividing by the OUTFLOW mass has ``rho = alpha M_-/M_+``.

    Red input: the module does not exist at the base tip. Mutation M-outflow-divisor: a copy of ``RankOneReemission``
    that divides by the OUTFLOW mass reads ``rho = 0.3 * 0.909`` against ``0.3``, on this row only.
    """
    R = H.ro()
    _, trace = H.asymmetric_trace()
    m_out, m_in = H.trace_weights(trace, "xmax")
    assert float(np.sum(m_out)) != float(np.sum(m_in)), "the fixture must carry M_+ != M_-"
    op = R.RankOneReemission(alpha=0.3, face="xmax", trace=trace)
    got = np.asarray(op.apply_transpose(m_in))
    ulps = np.max(np.abs(got - 0.3 * m_out) / np.spacing(np.abs(0.3 * m_out)))
    assert ulps <= math.ceil(m_out.size / 3), f"rho moved: {ulps} ulp"
    psi = H.positive_outflow(m_out.size)
    j_plus = float(np.sum(m_out * psi))
    j_minus = float(np.sum(m_in * np.asarray(op.apply(psi))))
    assert abs(j_minus / j_plus - 0.3) <= 1e-12


# ═════════════════════════════════════════════════════════════════════════════
# Gate 20: the declared chain and the realized image meet
# ═════════════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_typed_operators_answer_from_their_type")
@pytest.mark.parametrize("face", _FACES)
def test_the_declared_and_the_realized_composition_meet(face: str) -> None:
    r"""For each chain of the closed sum, ``chain.factors[0].domain`` is the tier of the outflow space the typed image
    derives from ``(face, trace)`` and ``chain.factors[-1].codomain`` the tier of its inflow space; for the diffuse
    chain the factor boundary ``factors[0].codomain`` is the tier of the ladder's ``current`` space. A CHECK of the
    declared tiers against the spaces, never a dispatch on ``role``.

    Red input: the module does not exist at the base tip. Mutation: a copy of ``DiffuseChain`` ending on ``"current"``
    (``factors = (CurrentRetraction,)``) reads ``"current" != "inflow"``.
    """
    from typing import get_args

    B = H.bc()
    R = H.ro()
    quad = H.QUADRATURES["gl8"]()
    trace = H.trace_of(quad)
    images = {
        B.SpecularChain: H.scaled_ordinate_permutation(1.0, H.mirror_pairing(quad, face), face, face, trace),
        B.DiffuseChain: R.RankOneReemission(alpha=1.0, face=face, trace=trace),
    }
    assert set(images) == set(get_args(B.Chain)), "a chain of the closed sum has no image in this gate"
    outflow, inflow, current = trace.outflow_space(face), trace.inflow_space(face), trace.current_space(face)
    for chain, image in images.items():
        assert chain.factors[0].domain == outflow.role == image.inner.domain.role, f"{chain.__name__}: the first factor does not start on the outflow tier"
        assert chain.factors[-1].codomain == inflow.role == image.inner.codomain.role, f"{chain.__name__}: the last factor does not end on the inflow tier"
    diffuse = B.DiffuseChain
    assert diffuse.factors[0].codomain == diffuse.factors[1].domain == "current"
    assert current.name == f"angular_trace[{face}:current]", "the current space is named for its tier"
