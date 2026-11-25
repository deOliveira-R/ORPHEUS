r"""G03f: the form a subquotient induces (N03, spec 7.4).

``Ann(Ann(d)) == d``; ``(in ∪ tang)/tang != Sub(in)`` and their duals differ; an undecided ``K ⊆ S`` raises. On a full-rank SPD dense ``G`` (``n = 4``):
the form of the density factor on ``S`` whole, ``K = Coord(n, (0,))`` is the Schur complement of ``G`` by index 0 to ``1e-12`` (two sections differ by
``e_0 r``: the coordinate restriction differs from it by O(1)), with kernel ``Coord(n, (0,))``; the form of the integrated factor on ``V/Coord(n, (0,))``
(the ``riesz`` image of the density factor on ``Sub(S)``, ``S = (1, 2, 3)``) is ``(G_SS)^-1`` to ``1e-12`` and not ``(G^-1)_SS``, and it is the density form
with its ``face`` flipped. Over EVERY pair ``(S, K)`` of ``n = 4`` and ``n = 6`` with a positive-dimensional quotient the integrated face has exactly
zero rows and columns on ``K ∪ (I minus S)`` (the leg that carries the inverse-route mutation: the full-matrix SVD route leaks on 4 of the 65 pairs
of ``n = 4``). The exact ``Coord`` admission: ``DenseMetric(diag([1, 2]), kernel=Coord(2, (1,)), face=METRIC)`` is refused ("kernel is not declared by
a spanned subspace") and ``DenseMetric(diag([1, 0]), kernel=Coord(2, (1,)), face=METRIC)`` builds. The dense ``Spanned`` arm: both forms are built and
``form(image) * form(V/Spanned)`` is the identity on the quotient to ``1e-12``. The zero-dimensional state on a dense axis is refused AT CONSTRUCTION
("the induced dense form is the zero form, which DenseMetric does not spell"); a deficient dense axis with a proper ``Coord`` ``S`` is refused with
"subspace intersection is undecided", and the raise in ``declaration.py::intersect`` carries a ``SCOPE-BOUNDARY[guard]`` tag the ledger's predicate
accepts. A constructed factor always reads: for every factor that constructs, over the full-rank, the ``Spanned``-kernel and the non-empty
``Coord``-kernel dense axes, ``.form`` of it, of its ``riesz`` and of its ``dual`` (each a construction that may refuse) raises nothing.

Exact on any host (a tolerance of ``1e-12`` between computations, exact zeros by construction of the live-block embedding, message fragments).
First red at ``ff305e0f``: import-only (``orpheus.numerics.factor``; ``Face`` and the ``face`` keyword do not exist). Mutations: replace the coset minimum by
the coordinate section (reddens by 3.09), use ``(G^-1)_SS`` for the integrated form (by 2.0e-01), the full-matrix SVD route for the inverse face (the
all-pairs leg), the exact ``Coord`` test dropped, the image form taken as the stated matrix with no flip (off the identity by 207).
"""

from __future__ import annotations

import numpy as np
import pytest

from tests.gates.numerics.w1 import _factor_probe as P

pytestmark = pytest.mark.foundation


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def _spd(n: int, seed: int = 5) -> np.ndarray:
    b = np.random.default_rng(seed).standard_normal((n, n))
    return b @ b.T + np.eye(n)


def _dense_axis_of(g: np.ndarray) -> object:
    from orpheus.numerics.axis import Axis, BasisKind, Declared
    from orpheus.numerics.declaration import Coord
    from orpheus.numerics.metric import DenseMetric, Face

    n = g.shape[0]
    return Axis("h", (n,), origin=Declared(BasisKind.NODAL, DenseMetric(g, kernel=Coord(n, ()), face=Face.METRIC)))


def test_subquotient_normal_form_and_induced_form() -> None:
    from orpheus.numerics import factor as F
    from orpheus.numerics.declaration import Ann, Coord, Spanned
    from orpheus.numerics.metric import DenseMetric, Face
    from orpheus.values.content import content_digest

    # normal forms
    d = Spanned(4, ("k", 1), 1)
    _require(d.ann().ann() == d, "Ann(Ann(d)) != d")
    inn, tang = Coord(6, (0, 1, 2)), Coord(6, (3, 4))
    union = Coord(6, (0, 1, 2, 3, 4))
    a = F.Subquotient(union, tang)
    b = F.Subquotient(inn, Coord(6, ()))
    _require(a != b, "(in ∪ tang)/tang equals Sub(in)")
    _require(a.dual() != b.dual(), "the duals of (in ∪ tang)/tang and Sub(in) are equal")
    # the form on a full-rank dense axis
    n, g = 4, _spd(4)
    ax = _dense_axis_of(g)
    whole = Coord(n, tuple(range(n)))
    f = P.density(ax, F.Subquotient(whole, Coord(n, (0,))))
    schur = g[1:, 1:] - np.outer(g[1:, 0], g[0, 1:]) / g[0, 0]
    got = np.asarray(f.form.matrix)[1:, 1:]
    _require(float(np.max(np.abs(got - schur))) <= 1e-12, f"the coset minimum differs from the Schur complement by {np.max(np.abs(got - schur)):.2e}")
    _require(float(np.max(np.abs(g[1:, 1:] - schur))) > 1e-1, "the contrast (coordinate restriction against the Schur complement) is inert")
    _require(f.form.kernel == Coord(n, (0,)), f"the kernel is {f.form.kernel}, not Coord(4, (0,))")
    # the integrated end: (G_SS)^-1 and not (G^-1)_SS; the density form with its face flipped
    S = (1, 2, 3)
    dens = P.density(ax, F.Subquotient(Coord(n, S), Coord(n, ())))
    img = dens.riesz()
    live = list(S)
    want = np.linalg.inv(g[np.ix_(live, live)])
    have = np.asarray(img.form.matrix)[np.ix_(live, live)]
    _require(float(np.max(np.abs(have - want))) <= 1e-12, f"the integrated form is off (G_SS)^-1 by {np.max(np.abs(have - want)):.2e}")
    _require(float(np.max(np.abs(np.linalg.inv(g)[np.ix_(live, live)] - want))) > 1e-1, "the (G^-1)_SS contrast is inert")
    flipped = dens.form.flipped()
    _require(content_digest(img.form) == content_digest(flipped) and flipped.face is Face.GRAM, "the integrated form is not the density form with its face flipped")
    # undecided K ⊆ S raises
    try:
        F.Subquotient(Coord(4, (0, 1)), Spanned(4, ("j", 2), 1))
    except ValueError:
        pass
    else:
        pytest.fail("a Subquotient with K = Spanned inside a Coord S was constructed (undecided inclusion must raise)")
    _require(not isinstance(d.ann(), Coord), "the annihilator of a Spanned is an Ann")
    _require(isinstance(d.ann(), Ann), "the annihilator of a Spanned is not an Ann")
    # the exact Coord admission
    try:
        DenseMetric(np.diag([1.0, 2.0]), kernel=Coord(2, (1,)), face=Face.METRIC)
    except ValueError as e:
        _require("kernel is not declared by a spanned subspace" in str(e), f"the refusal has another message: {e}")
    else:
        pytest.fail("a Coord kernel on a matrix with a nonzero row was admitted")
    m = DenseMetric(np.diag([1.0, 0.0]), kernel=Coord(2, (1,)), face=Face.METRIC)
    _require(m.kernel == Coord(2, (1,)), "the exact Coord kernel did not build")


@pytest.mark.parametrize("n", [4, 6])
@pytest.mark.rests_on("tests/gates/numerics/w1/test_factor_induced_form.py::test_subquotient_normal_form_and_induced_form")
def test_the_integrated_face_has_exact_zeros_on_the_kernel_and_outside_s(n: int) -> None:
    g = _spd(n, seed=7)
    ax = _dense_axis_of(g)
    checked = 0
    for S, K in P.pairs(n):
        sk = P.sk_of(n, S, K)
        if sk.dim == 0:
            continue
        f = P.density(ax, sk)
        face = np.asarray(f.form.inverse_matrix)
        dead = sorted(set(K) | (set(range(n)) - set(S)))
        if dead:
            _require(not bool(np.any(face[dead, :])) and not bool(np.any(face[:, dead])), f"n={n} S={S} K={K}: the integrated face has a nonzero entry on K ∪ (I minus S) = {dead}")
        checked += 1
    _require(checked == {4: 65, 6: 665}[n], f"{checked} pairs of positive dimension at n = {n}")


def test_dense_spanned_arm_and_the_declared_edges() -> None:
    from orpheus.numerics import factor as F
    from orpheus.numerics import declaration as D
    from orpheus.numerics.declaration import Coord

    ax, kernel, g = P.dense_axis()
    n = 4
    f = P.density(ax, F.Subquotient(Coord(n, tuple(range(n))), kernel))
    img = f.riesz()
    prod = np.asarray(img.form.matrix) @ np.asarray(f.form.matrix)
    w, v = np.linalg.eigh(g)
    ker = v[:, :1]
    q = np.eye(n) - ker @ ker.T  # projector on the representatives of the quotient, (ker G)-perp
    _require(float(np.max(np.abs(q @ (prod - np.eye(n)) @ q))) <= 1e-12, "form(image) * form(V/Spanned) is not the identity on the quotient")
    # the zero-dimensional state is refused at construction
    for S, K in (((0, 1), (0, 1)), ((0, 1, 2, 3), (0, 1, 2, 3))):
        try:
            P.density(ax, P.sk_of(n, S, K))
        except ValueError as e:
            _require(
                "the induced dense form is the zero form, which DenseMetric does not spell" in str(e) or "subspace intersection is undecided" in str(e) or "quotient" in str(e),
                f"S = K = {S}: refused with an unexpected message: {e}",
            )
        else:
            pytest.fail(f"the zero-dimensional state S = K = {S} on a dense axis was constructed")
    # the proper Coord S on the deficient axis: undecided meet
    try:
        P.density(ax, P.sk_of(n, (0, 1), ()))
    except ValueError as e:
        _require("subspace intersection is undecided" in str(e), f"a proper Coord S on a deficient dense axis: another refusal: {e}")
    else:
        pytest.fail("a proper Coord S on a deficient dense axis was constructed")
    # the zero-dimensional refusal on the whole space of a FULL-RANK dense axis carries the zero-form message
    full = _dense_axis_of(_spd(4))
    try:
        P.density(full, P.sk_of(4, (0, 1, 2, 3), (0, 1, 2, 3)))
    except ValueError as e:
        _require("the induced dense form is the zero form, which DenseMetric does not spell" in str(e), f"the zero-dimensional dense refusal reads: {e}")
    else:
        pytest.fail("the zero-dimensional state on a full-rank dense axis was constructed")
    # the intersect raise carries its declared-edge tag
    import inspect
    import re

    from tests.gates import test_elegance_debt_is_tagged as ledger

    src = inspect.getsource(D)
    lines = src.splitlines()
    sites = [i for i, ln in enumerate(lines) if "subspace intersection is undecided" in ln and "raise" in ln]
    _require(bool(sites), "the raise of intersect was not found in numerics/declaration.py")
    for i in sites:
        window = lines[max(0, i - 8): i + 4]
        _require(any(re.search(r"SCOPE-BOUNDARY\[guard\]", ln) for ln in window), f"the raise at declaration.py:{i + 1} carries no SCOPE-BOUNDARY[guard] tag")
    _require(not ledger._boundary_violations_in(src, label="declaration.py"), "declaration.py has a boundary tag missing machinery, ruling or revisit")


def _constructed(make):
    """The factor ``make()`` builds, or ``None`` when construction refuses it (a ``ValueError`` or ``TypeError``); a ``LinAlgError`` is a ``ValueError``
    and is NOT a refusal, so it is re-raised."""
    try:
        return make()
    except np.linalg.LinAlgError:
        raise
    except (ValueError, TypeError):
        return None


def test_a_constructed_factor_always_reads_its_form() -> None:
    from orpheus.numerics import factor as F
    from orpheus.numerics.declaration import Coord

    ax, kernel, _ = P.dense_axis()
    full = _dense_axis_of(_spd(4))
    coord_ax, _, _ = P.dense_coord_axis()
    n = 4
    pairs = P.pairs(n)
    spanned_extra = [F.Subquotient(Coord(n, tuple(range(n))), kernel), F.Subquotient(kernel.ann(), Coord(n, ()))]
    populations = {
        "spanned-kernel": [(P.sk_of(n, S, K), ax) for S, K in pairs] + [(sk, ax) for sk in spanned_extra],
        "full-rank": [(P.sk_of(n, S, K), full) for S, K in pairs],
        "coord-kernel": [(P.sk_of(n, S, K), coord_ax) for S, K in pairs],
    }
    for name, states in populations.items():
        built = 0
        for sk, axis in states:
            for v in P.variances():
                f = _constructed(lambda: F.Factor(axis, v, sk))
                if f is None:
                    continue  # refused at construction: the only place a state may be refused
                built += 1
                for label, other in (("itself", f), ("riesz", _constructed(f.riesz)), ("dual", _constructed(f.dual))):
                    if other is None:
                        continue  # the construction of the partner refused: that is its own construction refusal
                    try:
                        _ = other.form
                    except Exception as e:  # noqa: BLE001  any exception at a read of a constructed factor is the failure, whatever its type
                        pytest.fail(f"{name} ({sk.S}, {sk.K}, {v}) {label}: a constructed factor's .form raised {type(e).__name__}: {e}")
        _require(built > 0, f"no factor constructed on the {name} dense axis")
