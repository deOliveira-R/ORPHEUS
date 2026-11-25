r"""G04b, G04j, G04k, G04n: the Klein action of a factor lifts to the tensor space, slot by slot (N04, spec 7.5; laws L2, L4, L5 of 2.2).

G04b. On a fresh space ``V.dual().factor(l)`` is ``V.factor(l).dual()`` (flipped kind, chart and subquotient) and ``V.dual().axis(l) == V.axis(l)`` for every label, on a
whole space, on a partial dual and on a quotient (the zero-weight slot of ``Diag(0.3, -1.1, 0.9, 0)`` taken modulo its kernel).
G04j. Over a whole space, a ``Coord`` quotient of dimension 2 of 3 and a ``Spanned`` quotient, ``X.dual()`` and ``X.riesz()`` keep ``shape`` and each factor's ``stored_shape``.
G04k. L2: ``(S (x) T).dual() == S.dual() (x) T.dual()`` slot-wise and is not ``T.dual() (x) S.dual()``; L4: ``dual(on=L).dual(on=M) == dual(on=L triangle M)`` over all 64
label pairs of a three-slot space, ``dual(on=all) == dual()``, ``on=()`` is the identity (same object), ``KeyError`` on an unknown label; L5: ``riesz(on={"energy"})`` is
refused with the message "riesz is not defined: ". Each law runs on a positive and on a signed angular diagonal form.
G04n. Two spaces over equal axes whose variances differ in one slot are unequal and hash unequal; the one-slot ``riesz`` and ``dual`` images of a space differ from it and from each other.

First red at the base tip ``cf924aaa``: import-only (``TensorSpace`` is absent; the tip's ``V.dual()`` is a ``DualSpace`` with ``axes is None``, so ``.axis("angular")`` raises
``TypeError``). Mutations (in-process on the unit's tree, described): ``_make_dual`` returns the unflipped factors reddens G04b ("not flipped"); ``Factor.dual`` applied to
``axis.restricted(live)`` reddens G04j on the quotient and not on a whole space (the quotient is the input with teeth); permuting the slots in ``dual`` reddens L2; applying ``on=L``
twice reddens L4 only through ``dual(on=all) == dual()``; ``TensorSpace.identity_key`` dropping the variance reddens G04n. Exact on any host: identity relations.
"""

from __future__ import annotations

import itertools

import pytest

from tests.gates.numerics.w1 import _factor_probe as FP
from tests.gates.numerics.w1 import _tensor_probe as P

pytestmark = pytest.mark.foundation

LABELS = ("angular", "energy", "spatial")


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def _spaces() -> dict[str, object]:
    return {
        "positive": P.triple_space(),
        "signed": P.triple_space((1.3, -0.7, 2.1, 0.9)),
        "quotient": P.signed_pair_space(),
    }


def _quotients() -> dict[str, object]:
    from orpheus.numerics import factor as F
    from orpheus.numerics.declaration import Coord

    TS = P.tensor_space()
    coord = P.density_factor(P.diag_axis((1.0, 2.0, 3.0), "c"), sk=F.Subquotient(Coord(3, (0, 1, 2)), Coord(3, (2,))))
    axis, kernel, _g = FP.dense_axis(4, 3, True)
    spanned = F.Factor(axis, F.Variance(F.Kind.MEASURE, F.VarianceChart.DENSITY), F.Subquotient(Coord(4, (0, 1, 2, 3)), kernel))
    return {"whole": P.triple_space(), "coord 2 of 3": TS.of_axes(coord), "spanned": TS.of_axes(spanned)}


def test_dual_keeps_the_factors() -> None:
    for name, V in _spaces().items():
        labels = [f.axis.label for f in V.factors]
        D = V.dual()
        for l in labels:
            _require(D.factor(l) == V.factor(l).dual(), f"[{name}] dual().factor({l}) is not factor({l}).dual(): kind, chart or subquotient not flipped")
            _require(D.axis(l) == V.axis(l), f"[{name}] dual().axis({l}) differs from axis({l})")
        first = labels[0]
        part = V.dual(on={first})
        _require(part.factor(first) == V.factor(first).dual(), f"[{name}] the partial dual did not flip {first}")
        for l in labels[1:]:
            _require(part.factor(l) == V.factor(l), f"[{name}] the partial dual on {first} moved {l}")


def test_stored_shape_closed_on_quotients() -> None:
    for name, X in _quotients().items():
        for verb, image in (("dual", X.dual()), ("riesz", X.riesz())):
            _require(image.shape == X.shape, f"[{name}] {verb} moved shape {X.shape} -> {image.shape}")
            for f, g in zip(X.factors, image.factors):
                _require(g.stored_shape == f.stored_shape, f"[{name}] {verb} moved the stored shape of {f.axis.label}: {f.stored_shape} -> {g.stored_shape}")


def test_laws_l2_l4_l5() -> None:
    TS = P.tensor_space()
    for name in ("positive", "signed"):
        V = P.triple_space() if name == "positive" else P.triple_space((1.3, -0.7, 2.1, 0.9))
        a, e, s = V.factors
        # L2: slot-wise, order kept
        _require(V.dual() == TS.of_axes(a.dual(), e.dual(), s.dual()), f"[{name}] L2: the dual is not the slot-wise dual")
        _require(V.dual() != TS.of_axes(s.dual(), e.dual(), a.dual()), f"[{name}] L2: the dual equals the reversed slot-wise dual")
        # L4 over the 64 label pairs
        subsets = [set(c) for r in range(4) for c in itertools.combinations(LABELS, r)]
        _require(len(subsets) * len(subsets) == 64, "the label-pair population is not 64")
        for L in subsets:
            for M in subsets:
                _require(V.dual(on=L).dual(on=M) == V.dual(on=L ^ M), f"[{name}] L4: d_{sorted(L)} d_{sorted(M)} != d_{sorted(L ^ M)}")
        _require(V.dual(on=set(LABELS)) == V.dual(), f"[{name}] L4: dual(on=all) != dual()")
        _require(V.dual(on=()) is V, f"[{name}] L4: dual(on=()) is not the same object")
        with pytest.raises(KeyError):
            V.dual(on={"nope"})
        with pytest.raises(KeyError):
            V.riesz(on={"nope"})
        # L5: riesz is not defined on a formless slot
        with pytest.raises(Exception, match="riesz is not defined: "):
            V.riesz(on={"energy"})
        _require(V.riesz(on={"angular"}) == TS.of_axes(a.riesz(), e, s), f"[{name}] riesz(on=angular) is not the one-slot riesz")


def test_variance_is_in_the_identity() -> None:
    from orpheus.numerics import factor as F

    TS = P.tensor_space()
    a, e, s = P.triple_factors()
    other = F.Factor(a.axis, F.Variance(F.Kind.MEASURE, F.VarianceChart.INTEGRATED), a.sk)
    X, Y = TS.of_axes(a, e, s), TS.of_axes(other, e, s)
    _require(X != Y, "two spaces over equal axes with different angular variance are equal")
    _require(hash(X) != hash(Y), "two spaces whose variances differ in one slot hash equal")
    _require(X.riesz(on={"angular"}) == Y, "the one-slot riesz image is not the INTEGRATED-chart space")
    _require(X.riesz(on={"angular"}) != X.dual(on={"angular"}), "the one-slot riesz and dual images are equal")
    _require(X.riesz(on={"angular"}) != X and X.dual(on={"angular"}) != X, "a one-slot image equals its source")
