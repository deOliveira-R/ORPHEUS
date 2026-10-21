r"""S6a — the dual space carries G⁻¹ (rulings 26, 75; issues #617, #625).

The numerics rows; the rows on the S\ :sub:`N` carriers (C1–C3, the mint count
K2) are ``tests/gates/sn/operators/test_dual_metric_sn.py``. Row ids are those of
the verification spec (``scratch/numerics_untangling/packages/S6a_gate_spec.md``,
2026-10-21). The objects: :meth:`~orpheus.numerics.metric.HilbertMetric.dual`
returns the face-swap :class:`~orpheus.numerics.metric.DualMetric`;
:meth:`~orpheus.numerics.space.FunctionSpace.dual` is memoised per space, and a
composite's dual is the composite of its members' duals; :func:`~orpheus.numerics.operator.posed`
reads only the ends and the metric-free flag.

The laws, as the rows state them:

* a metric's dual is its two faces swapped, bit for bit, and ``dual().dual()``
  is the metric (M1–M4);
* ``V*`` reads its metric as the primal's dual; ``⟨f, g⟩_{V*} = ⟨♯f, ♯g⟩_V``;
  ``V* ≠ V`` in both orders; ``V.dual()`` is memoised and ``V.dual().dual() is V``
  (S1–S5, K1);
* a composite's dual is the composite of its members' duals (C1–C3);
* the Riesz legs admit dual spaces, and ``♭.H`` is ``♯`` by reciprocity (L1–L3);
* ``posed`` reads only the ends; no class defines ``is_adjoint_posed``; the
  dual arrow's adjoint is the base door, and the dagger–dual square commutes
  bitwise (P1, P2, D1–D3);
* the DualSpace mint count of a solve does not grow with its iteration count (K2);
* a spaceless pointwise operator acts on coordinates; a bound composite holding
  one on a dense metric has the right adjoint (Q1–Q3, #625);
* the composition guard refuses ``R @ R.dual()`` (G1, G2);
* the honesty gate's three DualSpace rows open (H1, H2).

``-O``-safe: every assertion is ``_require`` (a raise) or ``np.testing``.
"""

from __future__ import annotations

import ast
import inspect
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest

from orpheus.numerics import frame_faces
from orpheus.numerics import operator as op_mod
from orpheus.numerics.axis import Axis, BasisKind
from orpheus.numerics.elementary_operators import DiagonalOperator
from orpheus.numerics.metric import DenseMetric, DiagonalMetric, DualMetric, FactoredMetric
from orpheus.numerics.operator import (
    AdjointOperator,
    IdentityOperator,
    IncompatibleOperatorComposition,
    LinearOperator,
    RieszLowerOperator,
    RieszRaiseOperator,
    ZeroOperator,
    posed,
)
from orpheus.numerics.quadrature import Quadrature
from dataclasses import replace

from orpheus.numerics.space import DualSpace, FunctionSpace, TensorProductSpace

pytestmark = [
    pytest.mark.foundation,
    pytest.mark.rests_on(
        "tests/gates/numerics/test_riesz_legs.py::test_lowering_realizes_the_metric_pairing",
        "tests/gates/numerics/test_riesz_legs.py::test_round_trip_is_the_range_projector_on_a_singular_metric",
        "tests/gates/numerics/test_dense_metric.py::TestPairingSpelling::"
        "test_the_pairing_keeps_the_legacy_reduction_tree_bit_exactly",
    ),
]

_EPS = float(np.finfo(float).eps)
_OPERATOR_PY = Path(inspect.getfile(op_mod))


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


# ═══════════════════════════════════════════════════════════════════════
# Fixtures — every metric class, a singular member, exact-binary values
# ═══════════════════════════════════════════════════════════════════════

_W_POS = np.array([0.5, 2.0, 4.0])
_W_SING = np.array([2.0, 0.0, 3.0, 0.0])
#: test_dense_metric.py's exact-binary SPD literal.
_G = np.array([[2.0, 0.5, 0.0], [0.5, 1.0, 0.25], [0.0, 0.25, 4.0]])
#: a rank-2 PSD Gram (one dead direction), the Parseval-head shape.
_G_SING = np.array([[1.0, 1.0, 0.0], [1.0, 1.0, 0.0], [0.0, 0.0, 2.0]])


def _metrics():
    return {
        "diag_pos": (DiagonalMetric(_W_POS), (3,)),
        "diag_sing": (DiagonalMetric(_W_SING), (4,)),
        "dense": (DenseMetric(_G), (3,)),
        "dense_inverse_of_singular": (DenseMetric.inverse_of(_G_SING), (3,)),
        "factored_dense_diag": (
            FactoredMetric((((3,), DenseMetric(_G)), ((2,), DiagonalMetric(np.array([2.0, 0.5]))))),
            (3, 2),
        ),
        "factored_euclid_block": (
            FactoredMetric((((2,), None), ((3,), DiagonalMetric(_W_POS)))),
            (2, 3),
        ),
    }


def _leaf_spaces():
    w = np.array([2.0, 4.0])
    return {
        "named_weights": FunctionSpace("s6a_V", (3,), inner_product_weights=_W_POS),
        "named_singular": FunctionSpace("s6a_Vs", (4,), inner_product_weights=_W_SING),
        "dense": FunctionSpace("s6a_VD", (3,), metric=DenseMetric(_G)),
        "of_axes": FunctionSpace.of_axes(Axis("a", (2,), weights=w, kind=BasisKind.NODAL)),
        "basis_space_gl8_L2": frame_faces.basis_space(
            frame_faces.angular_frame(Quadrature.gauss_legendre(8), 2)
        ),
        "basis_space_product48_L4": frame_faces.basis_space(
            frame_faces.angular_frame(Quadrature.product(4, 8), 4)
        ),
        "tensor_dense_diag": FunctionSpace("s6a_T1", (3,), metric=DenseMetric(_G))
        * FunctionSpace("s6a_T2", (2,), inner_product_weights=np.array([2.0, 0.5])),
    }


def _sample(shape, seed):
    return np.random.default_rng(seed).standard_normal(shape)


# ═══════════════════════════════════════════════════════════════════════
# M — the metric's dual: a face swap, an involution
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("name", list(_metrics()))
def test_M1_the_dual_metric_is_the_face_swap_bitwise(name):
    """``m.dual().apply ≡ m.apply_inverse`` and ``m.dual().apply_inverse ≡ m.apply``,
    ``array_equal`` on an element with a trailing axis (leading-aligned
    broadcast). First red: a dual realised as ``DiagonalMetric(1/w)`` (x·(1/w) ≠
    x/w in bits; ∞ on a zero weight) or ``DenseMetric(pinv(G))`` (a second
    pseudo-inversion); a dual returning the primal (G for G⁻¹, O(1))."""
    m, shape = _metrics()[name]
    x = _sample(shape + (2,), 7)
    d = m.dual()
    np.testing.assert_array_equal(d.apply(x), m.apply_inverse(x))
    np.testing.assert_array_equal(d.apply_inverse(x), m.apply(x))
    # activation: on a non-Euclidean metric the two faces differ
    _require(not np.array_equal(d.apply(x), m.apply(x)), f"[{name}] fixture metric is its own inverse: the row is blind")


@pytest.mark.parametrize("name", ["diag_pos", "diag_sing", "dense"])
def test_M1b_the_positioned_block_form_swaps_too(name):
    """``dual.apply_block(inverse=b) ≡ primal.apply_block(inverse=not b)`` — the
    form a ``FactoredMetric`` composes (an axis-built space's metric). First
    red: a dual whose ``apply_block`` ignores the swap."""
    m, shape = _metrics()[name]
    x = _sample((2,) + shape + (2,), 8)
    for inverse in (False, True):
        np.testing.assert_array_equal(
            m.dual().apply_block(x, start=1, block_shape=shape, inverse=inverse),
            m.apply_block(x, start=1, block_shape=shape, inverse=not inverse),
        )


@pytest.mark.parametrize("name", list(_metrics()))
def test_M2_the_dual_is_an_involution_by_identity(name):
    """``m.dual().dual() is m`` and ``m.dual() is m.dual()`` (ruling Q3, 2026-10-21).
    First red: an uncached dual (a fresh object per call)."""
    m, _ = _metrics()[name]
    _require(m.dual() is m.dual(), f"[{name}] m.dual() is not cached")
    _require(m.dual().dual() is m, f"[{name}] m.dual().dual() is not m")


def test_M3_a_dense_dual_computes_no_new_inverse(monkeypatch):
    """The dense dual reads the two stored faces: it is the one face-swap
    wrapper (:class:`DualMetric`) over the dense metric, no ``pinv`` runs when
    it is built or applied, and its faces are the stored matrices bit for bit
    (census §1). First red: a dual that re-inverts (``DenseMetric(m.inverse_matrix)``
    derives a second pseudo-inverse: ``pinv`` runs and the inverse face moves)."""
    import orpheus.numerics.metric as metric_mod

    m = DenseMetric(_G)
    inv = m.inverse_matrix if m.inverse_matrix is not None else pytest.fail("primal has no inverse face")
    calls = []
    real = np.linalg.pinv
    monkeypatch.setattr(metric_mod.np.linalg, "pinv", lambda *a, **k: calls.append(1) or real(*a, **k))
    d = m.dual()
    x = _sample((3, 2), 5)
    fwd, back = d.apply(x), d.apply_inverse(x)
    _require(calls == [], f"the dense dual ran pinv {len(calls)} time(s)")
    _require(type(d) is DualMetric and isinstance(d, DualMetric) and d.primal is m,
             f"the dense dual is a {type(d).__name__}, not the face-swap wrapper over the metric")
    np.testing.assert_array_equal(fwd, inv @ x)
    np.testing.assert_array_equal(back, m.matrix @ x)


@pytest.mark.parametrize("name", list(_metrics()))
def test_M4_the_dual_pairing_has_one_spelling(name):
    """``m.dual().pairing(x, y) == Σ (G⁺x)⊙y`` bitwise: the dual inherits the one
    pairing spelling. First red: a dual overriding ``pairing`` (a matmul spelling
    moves bits, `[M]` up to 1792 ULP at n=15, test_dense_metric A10)."""
    m, shape = _metrics()[name]
    x, y = _sample(shape, 9), _sample(shape, 10)
    _require(m.dual().pairing(x, y) == float(np.sum(m.apply_inverse(x) * y)),
             f"[{name}] the dual pairing is not Σ(G⁺x)⊙y")


# ═══════════════════════════════════════════════════════════════════════
# S — V* carries G⁻¹ and is a distinct space
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("name", list(_leaf_spaces()))
def test_S1_the_dual_space_reads_the_inverse_metric_bitwise(name):
    """``V*.apply_metric ≡ V.apply_inverse_metric`` and the reverse, ``array_equal``.
    First red TODAY (#617): ``V*`` carries G (``[M]`` the old C2 row read the
    primal's 23.25 on the dense dual)."""
    V = _leaf_spaces()[name]
    x = _sample(V.shape, 11)
    Vd = V.dual()
    np.testing.assert_array_equal(Vd.apply_metric(x), V.apply_inverse_metric(x))
    np.testing.assert_array_equal(Vd.apply_inverse_metric(x), V.apply_metric(x))
    _require(not np.array_equal(Vd.apply_metric(x), V.apply_metric(x)),
             f"[{name}] V* and V apply the same metric (activation)")


def test_S2_the_dual_pairing_is_the_pairing_of_the_raises_exact_binary():
    """``⟨f, g⟩_{V*} = ⟨♯f, ♯g⟩_V``. Exact-binary diagonal case, hand value:
    w = (½, 2, 4), f = (1, 1, 1), g = (2, 4, 8): Σ fg/w = 4 + 2 + 2 = 8.
    First red TODAY: the dual pairs in G and reads ½·2 + 2·4 + 4·8 = 41."""
    V = FunctionSpace("s6a_S2", (3,), inner_product_weights=_W_POS)
    f, g = np.ones(3), np.array([2.0, 4.0, 8.0])
    _require(V.dual().inner_product(f, g) == 8.0, f"⟨f,g⟩_V* = {V.dual().inner_product(f, g)!r}, hand value 8.0")
    raise_ = RieszRaiseOperator(V)
    _require(V.inner_product(raise_.apply(f), raise_.apply(g)) == 8.0, "⟨♯f,♯g⟩_V ≠ 8.0")


def test_S2b_dense_dual_pairing_against_the_exact_rational_inverse():
    """On the dense SPD literal the reference is G⁻¹ in exact rationals (a
    Gauss–Jordan over ``Fraction``, structurally independent of ``pinv``):
    ``⟨f, g⟩_{V*} = fᵀG⁻¹g`` to 8 ε·κ relative; and ``⟨♯f, ♯g⟩_V`` agrees with it."""
    n = 3
    A = [[Fraction(v).limit_denominator(64) for v in row] + [Fraction(int(i == j)) for j in range(n)]
         for i, row in enumerate(_G.tolist())]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[p] = A[p], A[c]
        piv = A[c][c]
        A[c] = [v / piv for v in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                fac = A[r][c]
                A[r] = [a - fac * b for a, b in zip(A[r], A[c])]
    Ginv = np.array([[float(v) for v in row[n:]] for row in A])
    V = FunctionSpace("s6a_S2b", (3,), metric=DenseMetric(_G))
    f, g = _sample((3,), 12), _sample((3,), 13)
    exact = float(f @ Ginv @ g)
    kappa = float(np.linalg.cond(_G))
    tol = 8 * _EPS * kappa * max(1.0, abs(exact))
    got = V.dual().inner_product(f, g)
    _require(abs(got - exact) <= tol, f"⟨f,g⟩_V* = {got!r} vs fᵀG⁻¹g = {exact!r} (tol {tol:.1e})")
    raised = V.inner_product(RieszRaiseOperator(V).apply(f), RieszRaiseOperator(V).apply(g))
    _require(abs(raised - exact) <= tol, f"⟨♯f,♯g⟩_V = {raised!r} vs {exact!r}")


@pytest.mark.parametrize("name", list(_leaf_spaces()))
def test_S2c_the_dual_pairing_is_the_pairing_of_the_raises(name):
    """The same law on every leaf family, singular members included
    (G⁺GG⁺ = G⁺): relative 64 ε (a few reductions of ≤ 45 terms; dense heads
    carry their conditioning, read at landing)."""
    V = _leaf_spaces()[name]
    f, g = _sample(V.shape, 14), _sample(V.shape, 15)
    lhs = V.dual().inner_product(f, g)
    R = RieszRaiseOperator(V)
    rhs = V.inner_product(R.apply(f), R.apply(g))
    _require(abs(lhs - rhs) <= 64 * _EPS * max(1.0, abs(lhs), abs(rhs)) * 1e2,
             f"[{name}] ⟨f,g⟩_V* = {lhs!r} vs ⟨♯f,♯g⟩_V = {rhs!r}")


def test_S3_the_primal_is_the_duals_content():
    """``V*`` copies no metric source: its own ``inner_product_weights``, ``axes``
    and ``metric`` are ``None``; ``V*.primal is V``; and the ROUTE: a metric
    change on the primal (an in-place weight write, the battery idiom
    ``_resolved_metric`` is per-call for) moves ``V*``. First red TODAY: the
    dual's copied weights (``test_dual_preserves_weights``'s pin)."""
    w = np.array([1.0, 2.0, 4.0])
    V = FunctionSpace("s6a_S3", (3,), inner_product_weights=w.copy())
    Vd = V.dual()
    _require(isinstance(Vd, DualSpace) and Vd.primal is V, "V*.primal is not V")
    for f in ("inner_product_weights", "axes", "metric"):
        _require(getattr(Vd, f) is None, f"V* carries its own {f} (a copy of the primal's content)")
    x = np.ones(3)
    before = Vd.apply_metric(x).copy()
    object.__setattr__(V, "inner_product_weights", 2.0 * w)
    after = Vd.apply_metric(x)
    _require(not np.array_equal(before, after), "V* did not follow its primal's metric (a copied metric)")
    np.testing.assert_array_equal(after, x / (2.0 * w))


def _distinctness_families():
    fam = dict(_leaf_spaces())
    return fam


@pytest.mark.parametrize("name", list(_distinctness_families()))
def test_S4_V_star_is_not_V_in_both_orders(name):
    """``V* ≠ V`` and ``V ≠ V*`` (``==`` and ``!=``, both orders), the hash is
    ``hash(("dual", V))``, a set/dict holding V does not hold V*, and two mints of
    the dual are equal. First red TODAY (ruling 75): the axis-built families
    (``of_axes``, ``basis_space_*``, ``tensor``) read ``V.dual() == V`` with equal
    hashes (census §3(b))."""
    V = _distinctness_families()[name]
    Vd = V.dual()
    _require((Vd == V) is False and (V == Vd) is False, f"[{name}] V* == V")
    _require((Vd != V) is True and (V != Vd) is True, f"[{name}] V* != V is not True")
    if isinstance(Vd, DualSpace):  # a leaf dual; a product's dual is the product of duals (R1)
        _require(hash(Vd) == hash(("dual", V)), f"[{name}] hash(V*) is not hash(('dual', V))")
    _require(Vd not in {V} and {V: 1}.get(Vd) is None, f"[{name}] a container holding V holds V*")
    _require(DualSpace.of(V) == Vd and hash(DualSpace.of(V)) == hash(Vd), f"[{name}] two mints of V* differ")


def test_S5_duals_are_equal_iff_their_primals_are():
    a = FunctionSpace("s6a_same", (2,), inner_product_weights=np.array([2.0, 3.0]))
    b = FunctionSpace("s6a_same", (2,), inner_product_weights=np.array([5.0, 7.0]))
    c = FunctionSpace("s6a_other", (2,))
    _require(a == b and DualSpace.of(a) == DualSpace.of(b), "equal primals, unequal duals")
    _require(a != c and DualSpace.of(a) != DualSpace.of(c), "unequal primals, equal duals")


# ═══════════════════════════════════════════════════════════════════════
# K — the dual is memoised per space
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("name", list(_leaf_spaces()))
def test_K1_the_dual_is_memoised_per_space(name):
    """``V.dual() is V.dual()`` and ``V.dual().dual() is V`` (ruling 2026-10-21:
    memoised per space). First red TODAY: ``V.dual()`` mints a fresh DualSpace."""
    V = _leaf_spaces()[name]
    _require(V.dual() is V.dual(), f"[{name}] V.dual() is not memoised")
    _require(V.dual().dual() is V, f"[{name}] V.dual().dual() is not V")


def test_K1b_a_dual_arrows_ends_are_the_memoised_duals():
    V = FunctionSpace("s6a_K1b_V", (3,), inner_product_weights=_W_POS)
    W = FunctionSpace("s6a_K1b_W", (2,), inner_product_weights=np.array([3.0, 0.25]))
    A = _Rect(V, W)
    D = A.dual()
    _require(D.domain is W.dual() and D.codomain is V.dual(), "the dual arrow's ends are fresh duals")


# ═══════════════════════════════════════════════════════════════════════
# L — the Riesz legs on V*
# ═══════════════════════════════════════════════════════════════════════


def test_L1_the_legs_admit_a_dual_space():
    """``♭_{V*} : V* → V`` applies G⁻¹; ``♯_{V*} : V → V*`` applies G. INVERTS
    ``test_riesz_legs.py::test_the_legs_refuse_a_dual_space`` (today
    ``TypeError: PRIMAL space only``)."""
    V = FunctionSpace("s6a_L1", (3,), inner_product_weights=_W_POS)
    Vd = V.dual()
    low, high = RieszLowerOperator(Vd), RieszRaiseOperator(Vd)
    _require(low.domain == Vd and low.codomain == V, "♭_{V*} ends are not V* → V")
    _require(high.domain == V and high.codomain == Vd, "♯_{V*} ends are not V → V*")
    x = np.array([0.7, -1.2, 0.9])
    np.testing.assert_array_equal(low.apply(x), V.apply_inverse_metric(x))
    np.testing.assert_array_equal(high.apply(x), V.apply_metric(x))


def test_L2_the_double_lowering_is_the_identity():
    """The G² trap's witness, inverted: ``♭_{V*}(♭_V x) = G⁻¹Gx = x`` exactly on
    w = (½, 2, 4) (power-of-two weights: every product and quotient exact).
    `[M]` with ``V*`` carrying G this read ``G²x = [0.25, 4, 16]`` for x = 1."""
    V = FunctionSpace("s6a_L2", (3,), inner_product_weights=_W_POS)
    x = np.array([1.0, 1.0, 1.0])
    np.testing.assert_array_equal(RieszLowerOperator(V.dual()).apply(RieszLowerOperator(V).apply(x)), x)
    Vs = FunctionSpace("s6a_L2s", (4,), inner_product_weights=_W_SING)
    xs = np.array([1.0, 2.0, 3.0, 4.0])
    np.testing.assert_array_equal(
        RieszLowerOperator(Vs.dual()).apply(RieszLowerOperator(Vs).apply(xs)), np.where(_W_SING != 0, xs, 0.0)
    )


@pytest.mark.parametrize("weights", [_W_POS, _W_SING, np.array([0.3, 1.7, 2.9])], ids=["pos", "singular", "nonbinary"])
def test_L3_lowering_adjoint_is_raising_by_reciprocity(weights):
    """``♭_V.H`` is POSED (an ``AdjointOperator``, no override) and equals ``♯_V``:
    reciprocity ``⟨♭x, f⟩_{V*} = ⟨x, ♭.H f⟩_V`` and ``‖♭.H f − ♯f‖ ≤ 8 ε‖♯f‖``
    (G⁺GG⁺ vs G⁺: three roundings); the mirror ``♯.H`` vs ``♭``. First red TODAY:
    ``MissingAdjoint`` (DualSpace end). Mutation: a leg override ``adjoint`` returning
    the other leg (the row's AST leg)."""
    V = FunctionSpace("s6a_L3", (weights.size,), inner_product_weights=weights)
    low, high = RieszLowerOperator(V), RieszRaiseOperator(V)
    for cls in (RieszLowerOperator, RieszRaiseOperator):
        _require("adjoint" not in vars(cls) and "H" not in vars(cls), f"{cls.__name__} overrides its adjoint")
    x, f = _sample(weights.shape, 16), _sample(weights.shape, 17)
    lH = low.H
    _require(isinstance(lH, AdjointOperator), f"♭.H is a {type(lH).__name__}")
    lhs, rhs = V.dual().inner_product(low.apply(x), f), V.inner_product(x, lH.apply(f))
    _require(abs(lhs - rhs) <= 16 * _EPS * max(1.0, abs(lhs)), f"♭ reciprocity {lhs!r} vs {rhs!r}")
    np.testing.assert_allclose(lH.apply(f), high.apply(f), rtol=8 * _EPS, atol=0.0)
    np.testing.assert_allclose(high.H.apply(x), low.apply(x), rtol=8 * _EPS, atol=0.0)


# ═══════════════════════════════════════════════════════════════════════
# P, D — posed reads the ends; the dual arrow's adjoint is the base door
# ═══════════════════════════════════════════════════════════════════════


class _Rect(LinearOperator):
    """A bound ``A : V → W`` with a coded transpose (test_riesz_legs' fixture)."""

    _M = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])

    def __init__(self, V, W, M=None):
        self._domain, self._codomain = V, W
        if M is not None:
            self._M = np.asarray(M, float)

    @property
    def domain(self):
        return self._domain

    @property
    def codomain(self):
        return self._codomain

    def apply(self, x):
        return self._M @ x

    def apply_transpose(self, y):
        return self._M.T @ y

    @property
    def is_adjointable(self) -> bool:
        return True


class _Unbound(_Rect):
    @property
    def domain(self):
        return None


def test_P1_posed_reads_the_ends_only():
    """A DualSpace end poses; an unbound end does not (negative leg, vv #11).
    AST: no ``DualSpace`` name in ``posed``/``_ends_give_riesz_legs``/the legs'
    ``__init__``/``AdjointOperator.__init__``, and no ``isinstance(…, DualSpace)``
    anywhere in operator.py (TODAY 4: the leg refusals, the clause, the message
    branch; census §3). First red TODAY: ``posed(♭_V)`` is False."""
    V = FunctionSpace("s6a_P1V", (3,), inner_product_weights=_W_POS)
    W = FunctionSpace("s6a_P1W", (2,), inner_product_weights=np.array([3.0, 0.25]))
    _require(posed(RieszLowerOperator(V)) is True, "posed(♭_V) is not True")
    _require(posed(_Rect(V, W).dual()) is True, "posed(A.dual()) is not True")
    _require(posed(_Unbound(V, W)) is False, "an unbound operator reads posed (negative leg)")
    tree = ast.parse(_OPERATOR_PY.read_text())
    hits = [n.lineno for n in ast.walk(tree)
            if isinstance(n, ast.Call) and getattr(n.func, "id", None) == "isinstance"
            and any(isinstance(a, ast.Name) and a.id == "DualSpace" for a in ast.walk(n))]
    n_calls = sum(1 for n in ast.walk(tree) if isinstance(n, ast.Call) and getattr(n.func, "id", None) == "isinstance")
    _require(n_calls > 0, "AST control: no isinstance call found in operator.py (broken filter)")
    _require(hits == [], f"operator.py still discriminates on DualSpace at lines {hits}")


def test_P2_no_class_defines_is_adjoint_posed():
    """Ruling Q2 (2026-10-21): the attribute retires from the base AND from
    ``_DualOperator``; ``posed`` IS ``_ends_give_riesz_legs``. Walk: every
    ``LinearOperator`` subclass importable from ``orpheus`` (count printed in the
    message; control: ``AdjointOperator`` is in the walk). First red TODAY: 2
    definitions (``LinearOperator``, ``_DualOperator``)."""
    import importlib
    import pkgutil

    import orpheus

    for m in pkgutil.walk_packages(orpheus.__path__, "orpheus."):
        try:
            importlib.import_module(m.name)
        except Exception:  # noqa: BLE001 — an unimportable module is out of the walk, counted below
            continue
    seen, stack = set(), [LinearOperator]
    while stack:
        c = stack.pop()
        if c not in seen:
            seen.add(c)
            stack.extend(c.__subclasses__())
    _require(AdjointOperator in seen and len(seen) > 50, f"walk control failed ({len(seen)} classes)")
    defining = sorted(c.__qualname__ for c in seen if "is_adjoint_posed" in vars(c))
    _require(defining == [], f"{len(defining)} of {len(seen)} classes define is_adjoint_posed: {defining}")
    fn = ast.parse(inspect.getsource(op_mod.posed)).body[0]
    if not isinstance(fn, ast.FunctionDef):
        pytest.fail("AST control: posed did not parse as a function")
    body = fn.body[1:] if (fn.body and isinstance(fn.body[0], ast.Expr)
                           and isinstance(fn.body[0].value, ast.Constant)) else fn.body
    reads = [n for stmt in body for n in ast.walk(stmt)
             if (isinstance(n, ast.Attribute) and n.attr == "is_adjoint_posed")
             or (isinstance(n, ast.Constant) and n.value == "is_adjoint_posed")]
    _require(reads == [], "posed's code still reads is_adjoint_posed (attribute or getattr string)")


def test_D1_the_dual_arrow_has_no_adjoint_of_its_own():
    """``_DualOperator`` defines neither ``adjoint`` nor ``is_adjoint_posed``; its
    ``.H`` is the base door ``AdjointOperator(A.dual())``. First red TODAY: the
    commutation override returns ``_DualOperator(A.H)``."""
    _require("adjoint" not in vars(op_mod._DualOperator), "_DualOperator.adjoint survives")
    _require("is_adjoint_posed" not in vars(op_mod._DualOperator), "_DualOperator.is_adjoint_posed survives")
    _require("dual" not in vars(AdjointOperator), "AdjointOperator.dual survives (a copy of the base dual)")
    _require("dual" in vars(LinearOperator), "AST control: the base dual is absent")
    V = FunctionSpace("s6a_D1V", (3,), inner_product_weights=_W_POS)
    W = FunctionSpace("s6a_D1W", (2,), inner_product_weights=np.array([3.0, 0.25]))
    A = _Rect(V, W)
    DH = A.dual().H
    if type(DH) is not AdjointOperator or not isinstance(DH, AdjointOperator):
        pytest.fail(f"A.dual().H is a {type(DH).__name__}, not the base door")
    _require(DH.inner == A.dual(), "A.dual().H's inner is not A.dual()")
    _require(DH.H is DH.inner, "A.dual().H.H is not A.dual()")
    AH = A.H
    HD = AH.dual()
    _require(type(HD) is op_mod._DualOperator and getattr(HD, "inner", None) is AH,
             f"A.H.dual() is a {type(HD).__name__}, not the base dual arrow over A.H")


def _square_fixtures():
    V = FunctionSpace("s6a_sqV", (3,), inner_product_weights=_W_POS)
    W = FunctionSpace("s6a_sqW", (2,), inner_product_weights=np.array([3.0, 0.25]))
    VD = FunctionSpace("s6a_sqVD", (3,), metric=DenseMetric(_G))
    Vs = FunctionSpace("s6a_sqVs", (4,), inner_product_weights=_W_SING)
    rng = np.random.default_rng(23)
    return {
        "diag_rect": (_Rect(V, W), 3),
        "dense_endo": (_Rect(VD, VD, rng.standard_normal((3, 3))), 3),
        "singular_to_dense": (_Rect(Vs, VD, rng.standard_normal((3, 4))), 4),
    }


@pytest.mark.parametrize("name", list(_square_fixtures()))
def test_D2_the_dagger_dual_square_commutes_bitwise(name):
    """``(A.dual()).H ≡ (A.H).dual()``, ``array_equal`` — a LAW after S6a (two
    routes: the base sandwich over the dual metrics, and the transpose of the
    primal adjoint; same three delegations, face-swapped). Today green by the
    override (one route; DEMOTED); after the carve PROMOTED. First red after the
    carve: a dual metric that is not the face swap (bits), or ``V*`` carrying G
    (O(1)): ``[M]`` arm ``euclid_dual`` 0.43 on the SN slab (C3)."""
    A, n = _square_fixtures()[name]
    probe = _sample((n,), 24)
    np.testing.assert_array_equal(A.dual().H.apply(probe), A.H.dual().apply(probe))


# ═══════════════════════════════════════════════════════════════════════
# Q — #625: a spaceless pointwise operator acts on coordinates
# ═══════════════════════════════════════════════════════════════════════

_VD = FunctionSpace("s6a_Q_VD", (3,), metric=DenseMetric(_G))


def _wip(space, x, y):
    return float(space.inner_product(x, y))


@pytest.mark.parametrize("shape", ["D@A", "A@D", "2D@A+A", "D+A"])
def test_Q1_a_bound_composite_with_a_pointwise_term_on_a_dense_metric(shape):
    """``⟨Sx, y⟩_G = ⟨x, S.H y⟩_G`` for words and sums holding a spaceless
    ``D`` on the dense metric ``VD`` (64 ε relative); and the naive route through
    the term's own ``.H`` (``D.H is D``) misses (loaded fixture: ``[M]`` 1.26e-3
    on ``D + A``, #625). ``D+A`` is test_linear_combination c25 (REUSED: kept
    here only as the family's control row; drop at landing if c25 is cited)."""
    rng = np.random.default_rng(51)
    A = _Rect(_VD, _VD, rng.standard_normal((3, 3)))
    D = DiagonalOperator(np.array([1.0, 2.0, 3.0]))
    S = {"D@A": lambda: D @ A, "A@D": lambda: A @ D, "2D@A+A": lambda: 2.0 * (D @ A) + A, "D+A": lambda: D + A}[shape]()
    x, y = rng.standard_normal(3), rng.standard_normal(3)
    lhs, rhs = _wip(_VD, S.apply(x), y), _wip(_VD, x, S.H.apply(y))
    _require(abs(lhs - rhs) <= 64 * _EPS * max(1.0, abs(lhs)), f"[{shape}] reciprocity {lhs!r} vs {rhs!r}")
    naive = {"D@A": lambda: A.H @ D, "A@D": lambda: D @ A.H, "2D@A+A": lambda: 2.0 * (A.H @ D) + A.H,
             "D+A": lambda: D + A.H}[shape]()
    miss = abs(_wip(_VD, x, naive.apply(y)) - lhs) / max(1.0, abs(lhs))
    _require(miss > 1e-6, f"[{shape}] the naive per-term adjoint also passes ({miss:.1e}): the fixture is blind")


@pytest.mark.parametrize("cls", [IdentityOperator, ZeroOperator])
def test_Q2_a_constant_multiplier_is_its_own_adjoint_on_any_metric(cls):
    """``cI`` commutes with every metric, so ``(cI).H is cI`` (a
    :class:`ConstantMultiplier`, review F5/F6, #625): spaceless, metric-free,
    algebra-closed, and Euclidean reciprocity holds. First red: the metric-free
    law moved off the constants (``.H`` refuses as unbound)."""
    D = cls()
    _require(isinstance(D, op_mod.ConstantMultiplier), f"{cls.__name__} is not a ConstantMultiplier")
    _require(D.domain is None and D.codomain is None, f"{cls.__name__} has a bound end")
    _require(D.is_metric_free_adjoint is True and D.H is D, f"{cls.__name__}.H is not itself")
    x, y = _sample((3,), 61), _sample((3,), 62)
    _require(float(D.apply(x) @ y) == float(x @ D.apply(y)), "Euclidean reciprocity")


def test_Q2b_a_field_multiplier_has_no_spaceless_adjoint():
    """A diagonal FIELD multiplier commutes only with diagonal metrics (#625: on
    the dense literal ``D.H is D`` missed by 1.4e-3), so a spaceless ``D`` is not
    metric-free and its ``.H`` refuses as unbound; its Euclidean transpose
    ``D.dual()`` still runs and acts as ``D``; bound inside a composite, its
    adjoint is taken at the composite's ends (Q1). First red: the pre-review
    ``PointwiseOperator.adjoint`` returning self."""
    d = np.array([1.0, 2.0, 3.0])
    D = DiagonalOperator(d)
    _require(not isinstance(D, op_mod.ConstantMultiplier), "DiagonalOperator is a ConstantMultiplier")
    _require(D.is_metric_free_adjoint is False, "a field multiplier claims a metric-free adjoint")
    with pytest.raises(op_mod.MissingAdjoint, match="UNBOUND"):
        _ = D.H
    x = _sample((3,), 63)
    np.testing.assert_array_equal(D.dual().apply(x), D.apply(x))
    from orpheus.transport.operators.multiplication_operator import MultiplicationOperator

    flag = vars(MultiplicationOperator).get("is_metric_free_adjoint")
    if isinstance(flag, property) and flag.fget is not None:
        _require(flag.fget(None) is False, "MultiplicationOperator still declares a metric-free adjoint")


def test_Q3_no_docstring_claims_a_field_multiplier_commutes_with_every_metric():
    """X3, a prose gate (weak by kind; Q2/Q2b are the theorems): operator.py
    states neither "every shipped (diagonal) inner product" nor "commutes with
    every diagonal metric, so" (the pre-#625 claims), with a positive control
    on the phrase the review wrote ("commute with every metric")."""
    text = " ".join(_OPERATOR_PY.read_text().split())
    _require("commute with every metric" in text or "commutes with every metric" in text,
             "control: the constant-multiplier law is not stated")
    for phrase in ("every shipped (diagonal) inner product", "commutes with every diagonal metric, so"):
        _require(phrase not in text, f"operator.py still says {phrase!r}")


# ═══════════════════════════════════════════════════════════════════════
# G — the composition guard sees covariance
# ═══════════════════════════════════════════════════════════════════════


def test_G1_the_guard_refuses_R_at_R_dual():
    """``R : basis → measure``, ``R.dual() : measure* → basis*``: ``R @ R.dual()``
    joins basis and basis*, refused; positive leg: ``R.dual() @ ♭_measure``
    (measure → measure* → basis*) composes and applies. First red TODAY: the
    axis-built basis space equals its dual, so the guard composes (census §3(b))."""
    fr = frame_faces.angular_frame(Quadrature.gauss_legendre(8), 2)
    R = frame_faces.FrameReconstruction.of(fr)
    with pytest.raises(IncompatibleOperatorComposition, match="domain"):
        _ = R @ R.dual()
    m = frame_faces.measure_space(fr)
    w = R.dual() @ RieszLowerOperator(m)
    _require(np.shape(w.apply(np.ones(m.shape))) == tuple(R.domain.shape), "the positive leg does not apply")


def test_G2_a_sum_refuses_to_join_V_and_V_star():
    """``A + A.dual()`` for an endomorphism of an axis-built V: the ends V and V*
    disagree, ``_agreed_space`` refuses. First red TODAY: they compare equal."""
    V = frame_faces.basis_space(frame_faces.angular_frame(Quadrature.gauss_legendre(8), 2))
    n = int(np.prod(V.shape))
    A = _Rect(V, V, np.eye(n))
    with pytest.raises(IncompatibleOperatorComposition):
        _ = A + A.dual()


# ═══════════════════════════════════════════════════════════════════════
# H — the honesty gate's three DualSpace rows open
# ═══════════════════════════════════════════════════════════════════════

_H_ROWS = ["riesz_lower_primal", "riesz_raise_primal", "riesz_lower_primal_dual"]


@pytest.mark.parametrize("row", _H_ROWS)
def test_H1_the_dual_space_rows_open(row):
    """The iff gate (``test_each_transpose_door_opens_iff_predicted``) is green
    whether these rows open or stay refused (it predicts from ``posed``): this
    row is the one that reds if the clause survives. ``posed`` True, ``.H``
    applies, reciprocity in the operator's OWN end metrics (no override pairing),
    1e-12. First red TODAY: ``MissingAdjoint`` (DualSpace end)."""
    import tests.gates.numerics.test_adjointability_is_honest as Hn

    built = Hn._REGISTRY[row]()
    rng = np.random.default_rng(1729)
    x = Hn._like(built.x, rng)
    y = built.codomain_sample(rng)
    _require(posed(built.op) is True, f"[{row}] posed is False")
    image = built.op.H.apply(y)
    lhs, rhs = built.inner_w(built.act(x), y), built.inner_v(x, image)
    _require(Hn._rel(lhs, rhs) < 1e-12, f"[{row}] reciprocity {lhs!r} vs {rhs!r}")


def test_H2_the_dual_of_lowering_commutes():
    """``(♭.dual()).H ≡ (♭.H).dual()`` bitwise (``riesz_lower_primal_dual``'s value)."""
    V = FunctionSpace("s6a_H2", (3,), inner_product_weights=_W_POS)
    low = RieszLowerOperator(V)
    f = np.array([1.0, -2.0, 0.5])
    np.testing.assert_array_equal(low.dual().H.apply(f), low.H.dual().apply(f))


# ═══════════════════════════════════════════════════════════════════════
# R — review rows (2026-10-21): the dual of a product, the memo, the refusals
# ═══════════════════════════════════════════════════════════════════════

_RA = np.array([2.0, 4.0])
_RB = np.array([1.0, 2.0, 4.0])


def _r_factors(kind):
    if kind == "named":
        return (FunctionSpace("s6a_Ra", (2,), inner_product_weights=_RA),
                FunctionSpace("s6a_Rb", (3,), inner_product_weights=_RB))
    return (FunctionSpace.of_axes(Axis("ra", (2,), weights=_RA, kind=BasisKind.NODAL)),
            FunctionSpace.of_axes(Axis("rb", (3,), weights=_RB, kind=BasisKind.NODAL)))


@pytest.mark.parametrize("kind", ["named", "axes"])
def test_R1_a_product_with_a_dual_factor_carries_that_factors_inverse_form(kind):
    """``(V*⊗W).apply_metric(1) = outer(1/w_V, w_W)`` and ``(V⊗W*)`` the
    mirror, exactly (power-of-two weights); and ``(V⊗W)* == V*⊗W*``, memoised,
    with the same metric bit for bit. `[M]` before the review the dual factor
    read Euclidean in the product builder (the raw-field read)."""
    V, W = _r_factors(kind)
    one = np.ones((2, 3))
    np.testing.assert_array_equal((V.dual() * W).apply_metric(one), np.outer(1.0 / _RA, _RB))
    np.testing.assert_array_equal((V * W.dual()).apply_metric(one), np.outer(_RA, 1.0 / _RB))
    P = V * W
    _require(P.dual() == V.dual() * W.dual() and P.dual() is P.dual(), "(V⊗W)* is not V*⊗W*, memoised")
    _require(isinstance(P.dual(), TensorProductSpace), f"(V⊗W)* is a {type(P.dual()).__name__}")
    _require(P.dual() != P and P.dual().dual() is P, "(V⊗W)* == V⊗W, or its dual is not V⊗W")
    x = _sample((2, 3), 71)
    np.testing.assert_array_equal(P.dual().apply_metric(x), P.apply_inverse_metric(x))


def test_R2_dual_space_of_is_the_memoised_dual():
    """``DualSpace.of(V) is V.dual()`` (qa F2), and ``of(of(V)) is V``."""
    V = FunctionSpace("s6a_R2", (3,), inner_product_weights=_W_POS)
    P = _r_factors("named")[0] * _r_factors("named")[1]
    for S in (V, P):
        _require(DualSpace.of(S) is S.dual(), f"DualSpace.of({S!r}) is not its memoised dual")
        _require(DualSpace.of(DualSpace.of(S)) is S, "of(of(V)) is not V")


def test_R3_the_memo_survives_replace():
    """A ``replace``d dual is a fresh object whose dual is still the primal
    (``DualSpace._make_dual`` returns the primal, it never wraps a dual); a
    ``replace``d primal mints a dual equal to the original's. First red: a dual
    whose ``_make_dual`` wraps it again (V** a DualSpace over V*)."""
    V = FunctionSpace("s6a_R3", (3,), inner_product_weights=_W_POS)
    Vd2 = replace(V.dual())
    _require(Vd2 is not V.dual() and Vd2 == V.dual(), "replace(V*) is not a fresh equal dual")
    _require(Vd2.dual() is V, f"replace(V*).dual() is a {type(Vd2.dual()).__name__}, not V")
    V2 = replace(V)
    _require(V2.dual() == V.dual() and V2.dual().dual() is V2, "replace(V).dual() is wrong")


def test_R4_a_dual_space_refuses_content_of_its_own():
    """``DualSpace.__post_init__`` (review F4): a shape unlike its primal's, or any
    metric source of its own (weights, metric, axes), is refused; positive leg:
    the honest construction builds. First red: the refusal removed."""
    V = FunctionSpace("s6a_R4", (3,), inner_product_weights=_W_POS)
    _require(DualSpace(name="s6a_R4*", shape=(3,), primal=V) == V.dual(), "the honest construction failed")
    with pytest.raises(ValueError, match="is not its primal"):
        DualSpace(name="s6a_R4*", shape=(4,), primal=V)
    with pytest.raises(ValueError, match="metric source of its own"):
        DualSpace(name="s6a_R4*", shape=(3,), primal=V, inner_product_weights=_W_POS)
    with pytest.raises(ValueError, match="metric source of its own"):
        DualSpace(name="s6a_R4*", shape=(3,), primal=V, metric=DiagonalMetric(_W_POS))
    with pytest.raises(ValueError, match="metric source of its own"):
        DualSpace(name="s6a_R4*", shape=(2,), primal=_r_factors("axes")[0],
                  axes=_r_factors("axes")[0].axes)
