r"""G03a and G03k: the Klein action on a factor, exhaustively on a four-coordinate axis, and its stored shape (N03, spec 7.4).

G03a. The axis ``Diag([1, 2, 3, 0])`` has the kernel ``{3}``. Its 81 pairs ``(S, K)`` of coordinate subsets with ``K ⊆ S`` times the four
values of ``Variance`` are 324 states: 216 construct and 108 are refused (constructibility, P2-1: a density factor needs ``S ∩ ker G ⊆ K``,
an integrated factor needs its ``riesz`` image to be constructible). On the 216 the three laws hold (``riesz.riesz = id``, ``dual.dual = id``,
``dual.riesz = riesz.dual``). The laws are involutions of the variance and of ``sk.dual()`` and are true on all 324, so they carry no weight
alone: the COUNT and the FORM-LEVEL witness carry the constructibility rule (``f.riesz().form.apply(f.form.apply(psi))`` equals ``psi`` on
``S minus K``). On a formless axis (``EnergyAxis``) the 81 pairs times the two values of ``FormlessVariance`` are 162 states and all 162
construct (the 162 states the retired guard refused are not values of either type). A dense arm: a dense ``G`` of rank ``n - 1`` with a
``Spanned`` kernel: ``V/Spanned`` constructs, its ``riesz`` image ``Ann(Spanned)/Coord(n, ())`` constructs, ``riesz.riesz = id``, and a whole
``sk`` on the singular dense form is refused with the constructibility message.

G03k. Over the same 216 states, ``dual`` and ``riesz`` keep the stored shape (the axis's shape) and ``dim = dim S - dim K``.

First red at ``ff305e0f``: import-only (``orpheus.numerics.factor`` does not exist), read in each body at call time. Mutations: (i) the narrower
rule ``K == S ∩ ker G`` reddens the count (64 constructed, 260 refused); (ii) no constructibility check reddens the count (324 / 0) and the
form-level witness, while the three laws alone stay green; (iii) refusing the whole-space intersection ``Coord(n, all) ∩ d`` reddens the dense
arm. Exact on any host: counts of states and a ``1e-12`` comparison of a computation with a computation.
"""

from __future__ import annotations

import numpy as np
import pytest

from tests.gates.numerics.w1 import _factor_probe as P

pytestmark = pytest.mark.foundation


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def _enumerate() -> tuple[list, int]:
    from orpheus.numerics import factor as F

    ax = P.diag_axis([1.0, 2.0, 3.0, 0.0])
    built: list = []
    refused = 0
    for S, K in P.pairs(4):
        sk = P.sk_of(4, S, K)
        for v in P.variances():
            try:
                built.append(F.Factor(ax, v, sk))
            except ValueError:
                refused += 1
    return built, refused


def test_exhaustive_klein_on_a_four_coordinate_axis() -> None:
    from orpheus.numerics.axis import EnergyAxis
    from orpheus.numerics import factor as F

    built, refused = _enumerate()
    _require(len(P.pairs(4)) * 4 == 324, "the enumeration is not 81 pairs times 4 variances")
    _require((len(built), refused) == (216, 108), f"{len(built)} constructed and {refused} refused; 216 and 108 are the constructibility rule")
    for f in built:
        _require(f.riesz().riesz() == f, f"riesz.riesz != id on {f.sk.S.indices}/{f.sk.K.indices} {f.variance}")
        _require(f.dual().dual() == f, f"dual.dual != id on {f.sk.S.indices}/{f.sk.K.indices} {f.variance}")
        _require(f.dual().riesz() == f.riesz().dual(), f"dual.riesz != riesz.dual on {f.sk.S.indices}/{f.sk.K.indices} {f.variance}")
    # the form-level witness: the round trip is the identity on S minus K, on every constructed state of positive dimension
    ax = P.diag_axis([1.0, 2.0, 3.0, 0.0])
    psi = np.arange(1.0, 5.0)
    checked = 0
    for S, K in P.pairs(4):
        sk = P.sk_of(4, S, K)
        if sk.dim == 0:
            continue
        try:
            f = P.density(ax, sk)
        except ValueError:
            continue
        live = [i for i in S if i not in K]
        back = f.riesz().form.apply(f.form.apply(psi))
        err = float(np.max(np.abs(np.asarray(back)[live] - psi[live])))
        _require(err <= 1e-12, f"round trip fails on S={S} K={K}: |back - psi| on S minus K = {err:.3e}")
        checked += 1
    _require(checked > 20, f"the form-level witness read {checked} states")
    # the formless axis: 162 of 162
    e = EnergyAxis.synthetic(4)
    n_formless = 0
    for S, K in P.pairs(4):
        for v in F.FormlessVariance:
            F.Factor(e, v, P.sk_of(4, S, K))
            n_formless += 1
    _require(n_formless == 162, f"{n_formless} formless states construct, expected 162 of 162")


def test_the_dense_arm() -> None:
    from orpheus.numerics import factor as F
    from orpheus.numerics.declaration import Coord

    ax, kernel, _ = P.dense_axis()
    n = 4
    d = P.density(ax, F.Subquotient(Coord(n, tuple(range(n))), kernel))
    image = d.riesz()
    _require(image.sk == F.Subquotient(kernel.ann(), Coord(n, ())), f"the riesz image is {image.sk}, not Ann(Spanned)/Coord(n, ())")
    _require(image.riesz() == d, "riesz.riesz != id on the dense arm")
    try:
        P.density(ax, F.Subquotient.whole(n))
    except ValueError as e:
        _require("quotient" in str(e), f"the whole sk on the singular dense form is refused without naming the quotient: {e}")
    else:
        pytest.fail("a whole sk on a singular dense form was constructed")


@pytest.mark.rests_on("tests/gates/numerics/w1/test_factor_klein.py::test_exhaustive_klein_on_a_four_coordinate_axis")
def test_stored_shape_is_closed_under_both_involutions() -> None:
    built, _ = _enumerate()
    _require(len(built) == 216, f"{len(built)} constructed states, expected 216")
    for f in built:
        shape = f.axis.shape
        _require(f.stored_shape == shape, f"stored_shape {f.stored_shape} is not the axis shape {shape}")
        _require(f.dual().stored_shape == f.stored_shape, "dual moved the stored shape")
        _require(f.riesz().stored_shape == f.stored_shape, "riesz moved the stored shape")
        _require(f.dim == f.sk.S.dim - f.sk.K.dim, f"dim {f.dim} is not dim S - dim K")
