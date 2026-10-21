r"""S2 — a sum is a formal linear combination, a product a flattened word.

The verification gates of package S2 of the numerics untangling
(``.claude/plans/numerics_untangling.md``, rulings 23 and 50-53, and the
orchestrator's answers Q1-Q5 of 2026-10-20). The spec that carries each
gate's first red, its mutation and its re-baseline partner is
``scratch/numerics_untangling/packages/S2_gate_spec.md``.

The objects (the design ``S2_design.md`` section 1):

* a COMBINATION is ``{word -> coefficient}``: unordered for equality,
  insertion-ordered for the floating-point fold of ``apply``; terms collect by
  a KEY that compares atoms by object identity and never calls an atom's
  ``__eq__`` (D2); the formal arrows (``AdjointOperator``, the dual arrow,
  ``InverseOperator``) are keyed by (kind, key of the inner), plus the bound
  ends when the arrow is taken at spaces other than its inner's (Q5);
* a WORD is an atom or a flattened product ``(f1, ..., fn)``, ``n >= 2``; a
  scalar inside a word moves to the coefficient; identity factors stay (D7);
* every zero is the empty combination, ``ZeroMorphism(V, W)``, and is
  absorbed by ``+`` while contributing its spaces (D6, Q3);
* the adjoint distributes through the DUAL, the Riesz legs at the
  combination's ends: ``(Σ cᵢwᵢ).H = ♯_V ∘ (Σ cᵢ wᵢᵀ) ∘ ♭_W`` (D5, Q1);
* a sum of two or more terms has no inverse of its own; the inverse of
  ``A = M − N`` is ``GreenOperator(M, N)`` (D3, ruling 50);
* ``+`` never fuses: the fused composites are atoms built by their owner, each
  carrying ``expansion()`` (D4, ruling 52).

Claim kind: every row is a THEOREM (a law true for every admissible input),
except the rows marked RECORD (a float-arithmetic fact recorded as such).
The references are hand-composed from the atoms' own ``apply`` (no composite
on the reference side) or dense ``numpy`` matrices; nothing on the reference
side is built by the algebra under test.

V&V level L0 (term verification: each law of the type against a hand-built
reference, `vv-testing` "a type that embodies a mathematical concept ships
the test of its defining laws"). Each THEOREM row ``verifies`` the labelled
law it asserts (``operator_algebra.rst`` §composition-algebra,
``operator_adjoint.rst`` ``composite-adjoint-at-the-ends``,
``operator_inverse_family.rst`` ``green-split-inverse``), mapped by the
row's own assertion. Rows that carry no ``verifies``, deliberately: c05 (a
RECORD of the float limit of the scalar law: a limitation carries no
``verifies``), c19 (the composition action, no label), c31 (the role fold,
no label) and c32 (the retired idioms refuse: a software invariant).

The ADAPTER block below is the one place that names the representation's
accessors (``OperatorSum.terms``, ``OperatorProduct.factors``,
``expansion()``): if the landed spelling differs, edit it there.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from orpheus.numerics.metric import DenseMetric
from orpheus.numerics.operator import (
    BlockRole,
    IdentityOperator,
    InverseOperator,
    LinearOperator,
    MissingAdjoint,
    NotInvertible,
    OperatorProduct,
    OperatorSum,
    ZeroMorphism,
    ZeroOperator,
)
from orpheus.numerics.elementary_operators import DiagonalOperator
from orpheus.numerics.space import FunctionSpace

pytestmark = pytest.mark.l0

_HERE = "tests/gates/numerics/test_linear_combination.py"


# ───────────────────────────────────────────────────────────────────────
# ADAPTER — the representation's accessors (one place to re-spell)
# ───────────────────────────────────────────────────────────────────────


def _terms(op) -> list[tuple[float, object]]:
    """The (coefficient, word) pairs of ``op`` in insertion order; an atom or
    a word is the one-term ``[(1.0, op)]``."""
    if isinstance(op, OperatorSum):
        return [(float(c), w) for c, w in op.terms]
    return [(1.0, op)]


def _factors(word) -> tuple:
    """The factors of a word; an atom is its own one-factor word."""
    if isinstance(word, OperatorProduct):
        return tuple(word.factors)
    return (word,)


# ───────────────────────────────────────────────────────────────────────
# Fixtures: bound dense-matrix doubles on non-Euclidean spaces
# ───────────────────────────────────────────────────────────────────────


def _space(name: str, weights) -> FunctionSpace:
    w = np.asarray(weights, dtype=float)
    return FunctionSpace(name, (w.size,), inner_product_weights=w)


#: Distinct, non-uniform diagonal metrics on spaces of distinct dimension, so
#: a mis-placed Riesz leg or an unreversed word is a value error or a shape
#: error, never a coincidence.
V1 = _space("s2_v1", [1.0, 2.0, 0.5])
V2 = _space("s2_v2", [4.0, 0.25, 1.0, 2.0])
V3 = _space("s2_v3", [0.5, 1.0, 2.0, 4.0, 8.0])
V = V2

#: The exact-binary dense metric of ``test_dense_metric.py`` (every entry a
#: multiple of 1/4): the Q1 witness space, where a pointwise term's own
#: ``.H is self`` is not its adjoint.
_G = np.array([[2.0, 0.5, 0.0], [0.5, 1.0, 0.25], [0.0, 0.25, 4.0]])
VD = FunctionSpace("s2_dense3", (3,), metric=DenseMetric(_G))


class Mat(LinearOperator):
    """A bound dense-matrix operator ``domain -> codomain`` (the test double:
    its predicates are its constructor flags)."""

    def __init__(self, M, domain, codomain, *, invertible=False, adjointable=True,
                 block_role=None):
        self.M = np.asarray(M, dtype=float)
        self._domain, self._codomain = domain, codomain
        self._inv, self._adj = invertible, adjointable
        self.block_role = block_role
        self.system_role = None

    @property
    def domain(self):
        return self._domain

    @property
    def codomain(self):
        return self._codomain

    @property
    def is_invertible(self) -> bool:
        return self._inv

    @property
    def is_adjointable(self) -> bool:
        return self._adj

    def apply(self, x):
        return self.M @ x

    def apply_transpose(self, y):
        return self.M.T @ y

    def solve(self, b):
        return np.linalg.solve(self.M, b)

    def inverse(self) -> LinearOperator:
        if not self._inv:
            raise NotInvertible("Mat(invertible=False) has no inverse")
        return InverseOperator(self)


class _LiarAtom(Mat):
    """An atom whose ``__eq__`` says every atom is equal and whose hash is
    constant: collection by ``__eq__`` would merge two of them."""

    def __eq__(self, other):
        return True

    def __hash__(self):
        return 0


class _UnhashableAtom(Mat):
    """The ``StreamingOperator`` shape: ``__eq__`` raises, ``__hash__`` is
    ``None`` (a non-frozen dataclass)."""

    def __eq__(self, other):
        raise AssertionError("an atom's __eq__ was called by the algebra")

    __hash__ = None  # type: ignore[assignment]


class _SeedRecordingInverse(LinearOperator):
    """An inverse-family member recording the seed its driver hands it."""

    def __init__(self, forward):
        self.forward = forward
        self.seeds: list = []

    @property
    def domain(self):
        return self.forward.codomain

    @property
    def codomain(self):
        return self.forward.domain

    def apply(self, b, /, *, initial_guess=None):
        self.seeds.append(initial_guess)
        return np.linalg.solve(self.forward.M, b)


class _SeededMat(Mat):
    def __init__(self, M, space):
        super().__init__(M, space, space, invertible=True)
        self.recorder = _SeedRecordingInverse(self)

    def inverse(self) -> LinearOperator:
        return self.recorder


def _rng(seed: int = 20261020) -> np.random.Generator:
    return np.random.default_rng(seed)


def _mat(n_out, n_in, seed, dom, cod, **kw) -> Mat:
    r = _rng(seed)
    M = r.standard_normal((n_out, n_in))
    if n_out == n_in:
        M += n_in * np.eye(n_in)
    return Mat(M, dom, cod, **kw)


@pytest.fixture
def ABC():
    """Three distinct bound endomorphisms of V (each its own object)."""
    return (_mat(4, 4, 1, V, V, invertible=True), _mat(4, 4, 2, V, V, invertible=True),
            _mat(4, 4, 3, V, V, invertible=True))


def _wip(space: FunctionSpace, x, y) -> float:
    return float(space.inner_product(x, y))


# ═══════════════════════════════════════════════════════════════════════
# 1. The abelian-group laws and scalar distributivity, as equality
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.verifies("combination-normal-form")
def test_c01_sum_commutes(ABC):
    """``A + B == B + A``, with equal hashes, and the terms are the multiset
    {A: 1, B: 1}. First red today: ``==`` is object identity on two distinct
    sums. Mutation: equality reads the ORDERED term tuple (arm EQ_ORDERED)."""
    A, B, _ = ABC
    s1, s2 = A + B, B + A
    assert s1 == s2
    assert hash(s1) == hash(s2)
    assert {id(w): c for c, w in _terms(s1)} == {id(A): 1.0, id(B): 1.0}


@pytest.mark.verifies("combination-normal-form")
@pytest.mark.rests_on(_HERE + "::test_c01_sum_commutes")
def test_c02_sum_associates_by_flattening(ABC):
    """``(A + B) + C == A + (B + C)``: three terms, no term a sum. First red
    today: two distinct binary trees. Mutation: a sum operand stays one term
    (arm NO_FLATTEN_SUM)."""
    A, B, C = ABC
    left, right = (A + B) + C, A + (B + C)
    assert left == right
    assert len(_terms(left)) == 3 and len(_terms(right)) == 3
    assert not any(isinstance(w, OperatorSum) for _, w in _terms(right))


@pytest.mark.verifies("combination-scalar-laws")
@pytest.mark.rests_on(_HERE + "::test_c02_sum_associates_by_flattening")
def test_c03_scalar_distributes_over_a_sum(ABC):
    """``α(A + B) == αA + αB`` for a non-dyadic α (one multiply per side, so
    the coefficients agree bitwise). First red today: a ``ScaledOperator``
    over a sum. Mutation: scaling keeps the sum as one word (NO_DISTRIBUTE)."""
    A, B, _ = ABC
    alpha = 0.3
    assert alpha * (A + B) == alpha * A + alpha * B
    assert sorted(c for c, _ in _terms(alpha * (A + B))) == [alpha, alpha]


@pytest.mark.verifies("combination-scalar-laws")
@pytest.mark.rests_on(_HERE + "::test_c02_sum_associates_by_flattening")
def test_c04_coefficients_collect_over_repeated_atoms(ABC):
    """``(α + β)A == αA + βA`` and associativity over a repeated atom, with
    DYADIC coefficients (exact float sums). First red today: no collection.
    Mutation: no collection (NO_COLLECT)."""
    A, B, _ = ABC
    assert (0.5 + 0.25) * A == 0.5 * A + 0.25 * A
    assert _terms(0.5 * A + 0.25 * A) == [(0.75, A)]
    assert (0.5 * A + 0.25 * A) + 0.125 * A == 0.5 * A + (0.25 * A + 0.125 * A)
    assert A + A == 2.0 * A
    assert (A + B) + A == 2.0 * A + B


@pytest.mark.rests_on(_HERE + "::test_c04_coefficients_collect_over_repeated_atoms")
def test_c05_record_non_dyadic_collection_is_coefficient_rounded(ABC):
    """RECORD (Q2): equality is equality of the REPRESENTATION, so it is
    exact only where the coefficient arithmetic is. ``(0.1 + 0.2) + 0.3`` and
    ``0.1 + (0.2 + 0.3)`` differ in the last bit, so the two collections are
    unequal; and ``0.1A + 0.2A − 0.3A`` leaves the term
    ``5.551115123125783e-17·A``, not the zero. Float coefficient arithmetic,
    by ruling (no ``Fraction``, no tolerance in ``__eq__``). Red if either
    is introduced; red today (no collection)."""
    A, _, _ = ABC
    assert (0.1 * A + 0.2 * A) + 0.3 * A != 0.1 * A + (0.2 * A + 0.3 * A)
    residue = 0.1 * A + 0.2 * A - 0.3 * A
    assert not isinstance(residue, ZeroMorphism)
    assert _terms(residue) == [((0.1 + 0.2) - 0.3, A)]
    assert (0.1 + 0.2) - 0.3 == 5.551115123125783e-17


# ═══════════════════════════════════════════════════════════════════════
# 2. Collection by identity, never by an atom's __eq__ (D2)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.verifies("combination-normal-form", "combination-scalar-laws")
def test_c06_two_equal_valued_atoms_stay_two_terms():
    """Two distinct atoms whose ``__eq__`` says equal (and whose hashes
    coincide) are two terms: the combination is free on OBJECTS. First red
    today: no ``terms``. Mutation: key terms by the atom (KEY_EQ)."""
    M = np.eye(4) * 2.0
    a1, a2 = _LiarAtom(M, V, V), _LiarAtom(M, V, V)
    s = a1 + a2
    assert len(_terms(s)) == 2
    assert {id(w) for _, w in _terms(s)} == {id(a1), id(a2)}
    assert _terms(a1 + a1) == [(2.0, a1)]


@pytest.mark.verifies("combination-normal-form", "combination-scalar-laws")
@pytest.mark.rests_on(_HERE + "::test_c06_two_equal_valued_atoms_stay_two_terms")
def test_c07_an_unhashable_atom_combines_compares_and_hashes(ABC):
    """An atom whose ``__eq__`` raises and whose ``__hash__`` is ``None``
    (``StreamingOperator``'s shape) sums, compares, hashes and dedups in a
    set, and the algebra never calls its ``__eq__``. First red today: ``==``
    on two sums is identity (False). Mutation: KEY_EQ (TypeError: unhashable)."""
    A, _, _ = ABC
    R = _UnhashableAtom(np.eye(4), V, V)
    s1, s2 = R + A, A + R
    assert s1 == s2
    assert len({s1, s2}) == 1
    assert _terms(R + R) == [(2.0, R)]


@pytest.mark.verifies("combination-zero-and-unit")
@pytest.mark.rests_on(_HERE + "::test_c06_two_equal_valued_atoms_stay_two_terms", _HERE + "::test_c04_coefficients_collect_over_repeated_atoms")
def test_c08_cancellation_returns_the_surviving_atom(ABC):
    """``(A + B) − A`` IS ``B`` (the one-term collapse returns the object
    itself when its ends are the combination's), and ``L + C − L`` is ``C``.
    First red today: a three-node tree. Mutation: no collapse (NO_COLLAPSE)."""
    A, B, _ = ABC
    assert (A + B) - A is B
    assert A + B - A is B


# ═══════════════════════════════════════════════════════════════════════
# 3. Zero and the one-term collapse (D6, Q3)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.verifies("combination-zero-and-unit")
@pytest.mark.rests_on(_HERE + "::test_c08_cancellation_returns_the_surviving_atom")
def test_c09_zero_is_the_bound_zero_morphism():
    """``0·A``, ``A − A`` and ``(A + B) − (A + B)`` are ``ZeroMorphism(A.domain,
    A.codomain)`` on a rectangular ``A : V1 → V3``, and the zero applies to
    the codomain's zeros. First red today: ``ScaledOperator`` refuses a zero
    scalar (``ValueError``, "zero scalar"). Mutation: a zero coefficient kept
    as a term (ZERO_KEEP)."""
    A = _mat(5, 3, 4, V1, V3)
    B = _mat(5, 3, 5, V1, V3)
    for z in (0 * A, 0.0 * A, A - A, (A + B) - (A + B)):
        assert isinstance(z, ZeroMorphism)
        assert z.domain == V1 and z.codomain == V3
        np.testing.assert_array_equal(z.apply(np.ones(3)), np.zeros(5))


@pytest.mark.verifies("combination-zero-and-unit")
@pytest.mark.rests_on(_HERE + "::test_c09_zero_is_the_bound_zero_morphism")
def test_c10_an_unbound_zero_refuses_naming_the_term():
    """``0·U`` for an unbound ``U`` cannot mint its ``ZeroMorphism`` and
    refuses, naming the term's type (the message is the gate). First red
    today: the ``ScaledOperator`` refusal names no operand type."""
    U = Mat(np.eye(3), None, None)
    with pytest.raises(ValueError, match="Mat"):
        _ = 0.0 * U
    with pytest.raises(ValueError, match="ZeroMorphism"):
        _ = U - U


@pytest.mark.verifies("combination-zero-and-unit")
@pytest.mark.rests_on(_HERE + "::test_c09_zero_is_the_bound_zero_morphism", _HERE + "::test_c12_one_term_collapse")
def test_c11_every_zero_is_absorbed_and_contributes_its_spaces(ABC):
    """Q3: ``ZeroMorphism(V, W)``, the pointwise ``ZeroOperator()`` and the
    scalar 0 of ``sum(...)`` are absorbed; an absorbed zero still contributes
    its ends, so an UNBOUND atom plus a bound zero is a one-term combination
    carrying the zero's spaces, never the atom itself. First red today: ``A +
    ZeroMorphism`` is a two-node sum. Mutations: a zero kept as a term
    (ZERO_ATOM_KEPT); collapse ignoring the ends (COLLAPSE_ALWAYS)."""
    A, B, _ = ABC
    assert A + ZeroMorphism(domain=V, codomain=V) is A
    assert ZeroMorphism(domain=V, codomain=V) + A is A
    assert A + ZeroOperator() is A
    assert sum([A]) is A
    assert sum([A, B]) == A + B
    U = Mat(np.eye(4), None, None)
    bound = U + ZeroMorphism(domain=V, codomain=V)
    assert bound is not U
    assert bound.domain == V and bound.codomain == V
    assert _terms(bound) == [(1.0, U)]
    x = _rng(7).standard_normal(4)
    np.testing.assert_array_equal(bound.apply(x), U.apply(x))


@pytest.mark.verifies("combination-zero-and-unit")
@pytest.mark.rests_on(_HERE + "::test_c04_coefficients_collect_over_repeated_atoms")
def test_c12_one_term_collapse(ABC):
    """``1·A is A``, ``−(−A) is A``, ``A / 1 is A``, ``(2A)·0.5 is A``. First
    red today: each mints a ``ScaledOperator``. Mutation: NO_COLLAPSE."""
    A, _, _ = ABC
    assert 1 * A is A
    assert 1.0 * A is A
    assert -(-A) is A
    assert A / 1.0 is A
    assert (2.0 * A) * 0.5 is A


@pytest.mark.verifies("combination-zero-and-unit")
@pytest.mark.rests_on(_HERE + "::test_c09_zero_is_the_bound_zero_morphism")
def test_c38_the_zero_of_pointwise_terms_is_the_pointwise_zero(ABC):
    """Q8 (2026-10-20): a zero combination whose cancelled terms are all
    pointwise is ``ZeroOperator()`` (spaceless by law); one with a bound term
    is ``ZeroMorphism`` on that term's spaces (pointwise terms contribute
    none); only an unbound NON-pointwise term refuses (c10). First red today:
    ``ScaledOperator(0, I)`` refuses. Mutation: POINTWISE_ZERO_REFUSES."""
    A, _, _ = ABC
    I = IdentityOperator()
    D = DiagonalOperator(np.array([1.0, 2.0, 3.0, 4.0]))
    for z in (0.0 * I, I - I, D - D):
        assert isinstance(z, ZeroOperator)
    z = (D + A) - (D + A)
    assert isinstance(z, ZeroMorphism) and z.domain == V and z.codomain == V


# ═══════════════════════════════════════════════════════════════════════
# 4. apply: a left fold in insertion order, bit for bit
# ═══════════════════════════════════════════════════════════════════════


def _fold_fixture():
    """Four endomorphisms of V whose outputs span 16 decades, so the left
    and right folds round differently (asserted as the activation leg)."""
    big = Mat(np.diag([1e16, 3e16, 7e15, 1e16]), V, V)
    s1 = Mat(np.diag([1.0, 3.0, 1.0, 5.0]), V, V)
    s2 = Mat(np.diag([1.5, 0.5, 2.5, 1.0]), V, V)
    s3 = Mat(_rng(11).standard_normal((4, 4)), V, V)
    x = np.array([1.0, 1.0, 1.0, 1.0]) + _rng(12).standard_normal(4) * 1e-3
    return big, s1, s2, s3, x


@pytest.mark.verifies("combination-left-fold")
@pytest.mark.rests_on(_HERE + "::test_c02_sum_associates_by_flattening")
def test_c13_apply_is_the_left_fold_of_todays_left_nested_spelling():
    """``((A − S) − N) − B`` applies as ``((Ax − Sx) − Nx) − Bx`` bit for bit
    (today's binary arithmetic; the SN loss ``A_AA = LC − S − N2N − B_a`` has
    this shape), and a scaled term as ``Ax + (−σ)(Fx)``. Green today and
    after: the fold reproduces today's rounding. Activation leg: the right
    fold differs on this fixture. Mutations: FOLD_RIGHT, FOLD_REVERSED."""
    A, S, N, B, x = _fold_fixture()
    expected = ((A.apply(x) - S.apply(x)) - N.apply(x)) - B.apply(x)
    right = A.apply(x) + ((-S.apply(x)) + ((-N.apply(x)) + (-B.apply(x))))
    assert not np.array_equal(expected, right), "activation: the folds must differ here"
    np.testing.assert_array_equal((A - S - N - B).apply(x), expected)
    sigma = 0.3
    np.testing.assert_array_equal((A - sigma * S).apply(x), A.apply(x) + (-sigma) * S.apply(x))


@pytest.mark.verifies("combination-left-fold")
@pytest.mark.rests_on(_HERE + "::test_c13_apply_is_the_left_fold_of_todays_left_nested_spelling")
def test_c39_subtraction_is_bit_identical_to_adding_the_negation():
    """The fold subtracts for a −1 coefficient (no negated temporary, the
    allocation ceiling of ``test_p3_peak_allocation_per_dof``), and ``a − b``
    is bit for bit today's ``a + (−1)·b``: negation is exact in IEEE 754, and
    this row is that claim's witness on the fold fixture. Green today.
    Mutation: SUB_SIGN (the −1 term added)."""
    A, S, N, B, x = _fold_fixture()
    np.testing.assert_array_equal((A - S).apply(x), A.apply(x) + (-1.0) * S.apply(x))
    np.testing.assert_array_equal((A - S - N - B).apply(x),
                                  ((A.apply(x) + (-1.0) * S.apply(x)) + (-1.0) * N.apply(x)) + (-1.0) * B.apply(x))


@pytest.mark.verifies("combination-left-fold")
@pytest.mark.rests_on(_HERE + "::test_c13_apply_is_the_left_fold_of_todays_left_nested_spelling")
def test_c14_apply_does_not_alias_its_operand(ABC):
    """``(I + A).apply(x)`` leaves ``x`` untouched and returns a new array:
    the fold must not accumulate in place into a term's output (the identity
    returns ``x`` itself). Green today. Mutation: in-place ``+=`` (FOLD_INPLACE)."""
    A, _, _ = ABC
    x = _rng(13).standard_normal(4)
    x0 = x.copy()
    y = (IdentityOperator() + A).apply(x)
    np.testing.assert_array_equal(x, x0)
    assert y is not x
    np.testing.assert_array_equal(y, x0 + A.apply(x0))


@pytest.mark.verifies("combination-dual-distributes")
@pytest.mark.rests_on(_HERE + "::test_c13_apply_is_the_left_fold_of_todays_left_nested_spelling")
def test_c15_transpose_is_the_fold_of_the_terms_transposes(ABC):
    """``(A − 2B)ᵀy = Aᵀy + (−2)(Bᵀy)`` bit for bit, and the combination is
    adjointable iff every term is. Green today (survives)."""
    A, B, _ = ABC
    y = _rng(14).standard_normal(4)
    np.testing.assert_array_equal(
        (A - 2.0 * B).apply_transpose(y), A.apply_transpose(y) + (-2.0) * B.apply_transpose(y)
    )
    lacks = Mat(np.eye(4), V, V, adjointable=False)
    assert (A + lacks).is_adjointable is False
    with pytest.raises(MissingAdjoint):
        (A + lacks).apply_transpose(y)


# ═══════════════════════════════════════════════════════════════════════
# 5. The flattened word (D7)
# ═══════════════════════════════════════════════════════════════════════


@pytest.fixture
def chain():
    """``A : V2 → V3``, ``B : V1 → V2``, ``C : V3 → V1`` (distinct dims)."""
    return (_mat(5, 4, 21, V2, V3), _mat(4, 3, 22, V1, V2), _mat(3, 5, 23, V3, V1))


@pytest.mark.verifies("word-normal-form")
def test_c16_product_associates_by_flattening(chain):
    """``(AB)C == A(BC)``, factors ``(A, B, C)``, no factor a product. First
    red today: two binary trees; no ``factors``. Mutation: NO_FLATTEN_WORD."""
    A, B, C = chain
    left, right = (A @ B) @ C, A @ (B @ C)
    assert left == right
    assert _factors(left) == (A, B, C) == _factors(right)


@pytest.mark.verifies("word-normal-form")
@pytest.mark.rests_on(_HERE + "::test_c16_product_associates_by_flattening", _HERE + "::test_c12_one_term_collapse")
def test_c17_a_scalar_leaves_the_word(chain):
    """``A(2B) == 2(AB) == (2A)B``: one term, coefficient 2, word ``(A, B)``.
    First red today. Mutation: the scaled factor stays in the word
    (SCALAR_STAYS)."""
    A, B, _ = chain
    s = A @ (2.0 * B)
    assert s == 2.0 * (A @ B) == (2.0 * A) @ B
    [(c, w)] = _terms(s)
    assert c == 2.0 and _factors(w) == (A, B)


@pytest.mark.verifies("word-normal-form")
@pytest.mark.rests_on(_HERE + "::test_c16_product_associates_by_flattening")
def test_c18_identity_factors_are_kept_and_pointwise_ends_resolve(chain):
    """``A @ I`` is the word ``(A, I)``, not ``A`` (D7: the unit law is not
    claimed); a pointwise end factor takes its neighbour's space. First red
    today: no ``factors``. Mutation: DROP_IDENTITY."""
    A, B, _ = chain
    I = IdentityOperator()
    w = A @ I
    assert w is not A and _factors(w) == (A, I)
    D = DiagonalOperator(np.array([1.0, 2.0, 3.0, 4.0, 5.0]))
    D1 = DiagonalOperator(np.array([1.0, 2.0, 3.0]))
    assert (D @ A).codomain == A.codomain and (D @ A).domain == A.domain
    assert (A @ B @ D1).domain == B.domain
    assert (D @ A @ B @ D1).codomain == A.codomain


@pytest.mark.rests_on(_HERE + "::test_c16_product_associates_by_flattening")
def test_c19_word_applies_right_to_left(chain):
    """``(ABC)x = A(B(Cx))`` bit for bit. Green today (survives). Mutation:
    left-to-right (WORD_LTR, a shape error on these distinct dims)."""
    A, B, C = chain
    x = _rng(24).standard_normal(5)
    np.testing.assert_array_equal((A @ B @ C).apply(x), A.apply(B.apply(C.apply(x))))


@pytest.mark.verifies("word-reversal-laws")
@pytest.mark.rests_on(_HERE + "::test_c19_word_applies_right_to_left")
def test_c20_word_inverse_and_transpose_reverse(ABC):
    """``(ABC)⁻¹ = C⁻¹B⁻¹A⁻¹`` and ``(ABC)ᵀ = CᵀBᵀAᵀ``, bit for bit against
    the hand-spelled sequence. Green today (survives)."""
    A, B, C = ABC
    y = _rng(25).standard_normal(4)
    np.testing.assert_array_equal(
        (A @ B @ C).inverse().apply(y), C.solve(B.solve(A.solve(y)))
    )
    np.testing.assert_array_equal(
        (A @ B @ C).apply_transpose(y),
        C.apply_transpose(B.apply_transpose(A.apply_transpose(y))),
    )


# ═══════════════════════════════════════════════════════════════════════
# 6. The adjoint distributes through the dual (D5, Q1, Q5)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.verifies("combination-dual-distributes", "composite-adjoint-at-the-ends")
@pytest.mark.rests_on(_HERE + "::test_c23_arrow_keys_separate_their_inners", _HERE + "::test_c24_sum_adjoint_reciprocity_on_weighted_spaces", _HERE + "::test_c03_scalar_distributes_over_a_sum")
def test_c21_sum_dual_distributes_and_its_adjoint_is_taken_at_the_ends():
    """Ruling 53's precise form (Q1/Q6, 2026-10-20): the DUAL distributes,
    ``(2A − 0.5B)ᵀ == 2Aᵀ − 0.5Bᵀ`` as representations; the ADJOINT is ``♯ ∘
    dual ∘ ♭`` taken once at the combination's ends, equal across separate
    builds, and bit for bit the hand composition ``♯_{V2}(2Aᵀ(♭_{V3}y) +
    (−0.5)Bᵀ(♭_{V3}y))`` (today's whole-sum route). The theorem ``S.H = 2A.H −
    0.5B.H`` holds here (bound, non-degenerate) and is gated in value only.
    First red today: ``S.dual()`` and ``S.H`` are not comparable across builds.
    Mutations: DUAL_UNDISTRIBUTED, ARROW_ID, LEGS_EUCLIDEAN, ADJ_NAIVE (bits)."""
    from orpheus.numerics.operator import RieszLowerOperator, RieszRaiseOperator

    A = _mat(5, 4, 31, V2, V3)
    B = _mat(5, 4, 32, V2, V3)
    S = 2.0 * A - 0.5 * B
    assert S.dual() == 2.0 * A.dual() - 0.5 * B.dual()
    assert S.H == (2.0 * A - 0.5 * B).H
    assert S.H.domain == V3 and S.H.codomain == V2
    y = _rng(33).standard_normal(5)
    z = RieszLowerOperator(V3).apply(y)
    expected = RieszRaiseOperator(V2).apply(2.0 * A.apply_transpose(z) + (-0.5) * B.apply_transpose(z))
    np.testing.assert_array_equal(S.H.apply(y), expected)
    theorem = (2.0 * A.H - 0.5 * B.H).apply(y)
    np.testing.assert_allclose(S.H.apply(y), theorem, rtol=64 * np.finfo(float).eps, atol=0.0)


@pytest.mark.verifies("word-reversal-laws", "composite-adjoint-at-the-ends")
@pytest.mark.rests_on(_HERE + "::test_c16_product_associates_by_flattening", _HERE + "::test_c23_arrow_keys_separate_their_inners")
def test_c22_word_dual_reverses_and_its_adjoint_is_taken_at_the_ends(chain):
    """``(AB)ᵀ == BᵀAᵀ`` as representations; ``(AB).H`` is one arrow at the
    word's ends, equal across builds, bit for bit ``♯_{V1}(Bᵀ(Aᵀ(♭_{V3}y)))``;
    the theorem ``(AB).H = B.H A.H`` holds on this non-degenerate middle and is
    gated in value only (c36 is the degenerate middle). First red today.
    Mutations: WORD_ADJ_UNREVERSED, ARROW_ID, LEGS_EUCLIDEAN, WORD_ADJ_NAIVE (bits)."""
    from orpheus.numerics.operator import RieszLowerOperator, RieszRaiseOperator

    A, B, _ = chain
    assert (A @ B).dual() == B.dual() @ A.dual()
    assert (A @ B).H == (A @ B).H
    y = _rng(34).standard_normal(5)
    expected = RieszRaiseOperator(V1).apply(B.apply_transpose(A.apply_transpose(RieszLowerOperator(V3).apply(y))))
    np.testing.assert_array_equal((A @ B).H.apply(y), expected)
    np.testing.assert_allclose((A @ B).H.apply(y), (B.H @ A.H).apply(y), rtol=64 * np.finfo(float).eps, atol=0.0)


@pytest.mark.verifies("combination-scalar-laws")
def test_c23_arrow_keys_separate_their_inners(ABC):
    """Q5's witness: a formal arrow is keyed by (kind, its inner's key), so two
    separately built ``A.H`` are equal and collect, while ``A.H != B.H`` and
    ``A.H + B.H`` keeps two terms; the same for ``.inverse()`` and ``.dual()``.
    First red today: arrows compare by identity. Mutations: ARROW_ID (the
    first legs), ARROW_KIND_ONLY (the second)."""
    A, B, _ = ABC
    assert A.H == A.H and A.H != B.H
    assert _terms(A.H + A.H)[0][0] == 2.0 and len(_terms(A.H + A.H)) == 1
    assert len(_terms(A.H + B.H)) == 2
    assert A.dual() == A.dual() and A.dual() != B.dual()
    assert len(_terms(A.inverse() + A.inverse())) == 1
    assert len(_terms(A.inverse() + B.inverse())) == 2


@pytest.mark.verifies("composite-adjoint-at-the-ends")
def test_c24_sum_adjoint_reciprocity_on_weighted_spaces():
    """``⟨Sx, y⟩_{V3} = ⟨x, S.H y⟩_{V2}`` for ``S = 2A − 0.5B : V2 → V3``,
    distinct non-uniform weights. Green today and after; the tolerance is
    the reduction depth (a few ulp, relative). Mutation: the Euclidean
    transpose in place of the adjoint (LEGS_EUCLIDEAN)."""
    A = _mat(5, 4, 41, V2, V3)
    B = _mat(5, 4, 42, V2, V3)
    S = 2.0 * A - 0.5 * B
    r = _rng(43)
    x, y = r.standard_normal(4), r.standard_normal(5)
    lhs, rhs = _wip(V3, S.apply(x), y), _wip(V2, x, S.H.apply(y))
    assert abs(lhs - rhs) <= 64 * np.finfo(float).eps * max(1.0, abs(lhs))


@pytest.mark.verifies("composite-adjoint-at-the-ends")
@pytest.mark.rests_on(_HERE + "::test_c24_sum_adjoint_reciprocity_on_weighted_spaces")
def test_c25_spaceless_term_adjoint_is_taken_at_the_combinations_spaces():
    """Q1's witness. ``D`` is spaceless: it acts on coordinates, where
    ``D.H is D`` is a theorem (#625, narrowed at S6a). Inside ``D + A`` on
    ``VD`` (a dense, non-diagonal metric) the adjoint is taken at the
    combination's spaces, never through the term's own ``.H``. `[M]` 2026-10-20
    (probe ``S2_gates/probes/p1_pointwise_dense.py``): today's whole-sum
    route holds to 2.1e-16 relative, the naive ``D.H + A.H`` misses by
    1.26e-3. Green today; red under ADJ_NAIVE."""
    A = Mat(_rng(51).standard_normal((3, 3)), VD, VD)
    D = DiagonalOperator(np.array([1.0, 2.0, 3.0]))
    S = D + A
    r = _rng(52)
    x, y = r.standard_normal(3), r.standard_normal(3)
    lhs, rhs = _wip(VD, S.apply(x), y), _wip(VD, x, S.H.apply(y))
    assert abs(lhs - rhs) <= 64 * np.finfo(float).eps * max(1.0, abs(lhs))


@pytest.mark.verifies("composite-adjoint-at-the-ends")
@pytest.mark.rests_on(_HERE + "::test_c21_sum_dual_distributes_and_its_adjoint_is_taken_at_the_ends", _HERE + "::test_c25_spaceless_term_adjoint_is_taken_at_the_combinations_spaces")
def test_c26_involution_equality_on_a_sum_identity_on_a_leaf(ABC):
    """``S.H.H == S`` on a sum (a new object, equal), including a spaceless
    term on the dense metric; ``A.H.H is A`` on a leaf. Green today (today
    the sum's is identity); the row's tooth is an adjoint-at arrow whose own
    adjoint re-wraps instead of returning its term (ADJAT_REWRAP)."""
    A, B, _ = ABC
    S = 2.0 * A - B
    assert S.H.H == S
    assert A.H.H is A
    Ad = Mat(_rng(53).standard_normal((3, 3)), VD, VD)
    D = DiagonalOperator(np.array([1.0, 2.0, 3.0]))
    assert (D + Ad).H.H == D + Ad


@pytest.mark.verifies("composite-adjoint-at-the-ends")
@pytest.mark.rests_on(_HERE + "::test_c24_sum_adjoint_reciprocity_on_weighted_spaces")
def test_c36_word_adjoint_reciprocity_through_a_singular_middle():
    """Q6's witness. ``⟨ABx, y⟩_W = ⟨x, (AB).H y⟩_V`` when the MIDDLE space U
    has a singular metric (a zero weight: the SN partial-current trace's
    tangential slots until S6). The word's adjoint must carry the Riesz legs
    at its ENDS only: ``B.H @ A.H`` inserts ``♭_U ♯_U = G_U G_U⁺``, a projector.
    `[M]` 2026-10-20 (probe ``S2_gates/probes/p2_word_adjoint_singular_middle.py``):
    the whole-word route holds exactly, the distributed one misses by 6.4e-2
    relative. Green today; red under the literal distribution (the
    prototype's default until Q6 is ruled)."""
    U = _space("s2_u_singular", [1.0, 0.0, 2.0])
    Vs = _space("s2_v_s", [1.0, 2.0, 0.5])
    Ws = _space("s2_w_s", [4.0, 1.0, 0.25])
    r = _rng(71)
    A = Mat(r.standard_normal((3, 3)), U, Ws)
    B = Mat(r.standard_normal((3, 3)), Vs, U)
    x, y = r.standard_normal(3), r.standard_normal(3)
    lhs, rhs = _wip(Ws, (A @ B).apply(x), y), _wip(Vs, x, (A @ B).H.apply(y))
    assert abs(lhs - rhs) <= 64 * np.finfo(float).eps * max(1.0, abs(lhs))


@pytest.mark.verifies("composite-adjoint-at-the-ends")
@pytest.mark.rests_on(_HERE + "::test_c22_word_dual_reverses_and_its_adjoint_is_taken_at_the_ends", _HERE + "::test_c13_apply_is_the_left_fold_of_todays_left_nested_spelling")
def test_c37_the_adjoint_pencil_member_keeps_todays_bits():
    """The production adjoint pencil member ``pencil.H.at(2.5)`` on a carrying
    fissile sphere (its fission ``rhs`` is the WORD ``CoupledOperator @
    SystemRestrictionOperator``) is bit for bit today's route:
    ``lhs.H y + (−2.5)·♯_V(rhsᵀ(♭_W y))``. `[M]` 2026-10-20: the literal
    word-adjoint distribution moved it by 1.8e-19 relative (bit probe under the
    prototype); ruling 53's precise form keeps it. Green today. Mutation:
    WORD_ADJ_NAIVE."""
    from orpheus.numerics.operator import RieszLowerOperator, RieszRaiseOperator
    from tests.gates.sn.operators.test_psi_half_coupling import _fissile_sphere

    from orpheus.numerics.operator import adjointable

    problem = _fissile_sphere()
    pencil = problem.pencil
    rhs = pencil.rhs
    if rhs.domain is None or rhs.codomain is None or not adjointable(rhs):
        pytest.fail("premise: the carrying sphere's fission word is bound and transposes")
    zero = problem.system.space.zeros()
    n = np.asarray(zero.to_flat()).size
    for seed in range(16):  # the moved bits are draw-dependent (anti-#31): sweep
        y = type(zero).from_flat(_rng(seed).random(n), zero)
        rhs_adjoint_y = RieszRaiseOperator(rhs.domain).apply(
            rhs.apply_transpose(RieszLowerOperator(rhs.codomain).apply(y))
        )
        expected = pencil.lhs.H.apply(y) + (-2.5) * rhs_adjoint_y
        np.testing.assert_array_equal(
            np.asarray(pencil.H.at(2.5).apply(y).to_flat()), np.asarray(expected.to_flat()),
            err_msg=f"seed {seed}",
        )


def _split():
    """``M`` diagonally dominant and invertible, ``N = S + F`` small, on V."""
    r = _rng(61)
    M = Mat(np.diag([4.0, 5.0, 6.0, 7.0]) + 0.1 * r.standard_normal((4, 4)), V, V,
            invertible=True)
    S = Mat(0.3 * r.standard_normal((4, 4)), V, V)
    F = Mat(0.2 * r.standard_normal((4, 4)), V, V)
    return M, S, F


@pytest.mark.verifies("combination-inverse")
@pytest.mark.rests_on(_HERE + "::test_c12_one_term_collapse")
@pytest.mark.parametrize("order", ["implicit_first", "implicit_last"])
def test_c27_a_sum_of_two_terms_has_no_inverse_of_its_own(order):
    """In either spelling order a two-term sum is not invertible and its
    ``inverse()`` refuses, naming the split. First red today: ``M − S`` is
    "invertible" by its head and returns a Green. Mutation: the head contract
    (HEAD)."""
    M, S, _ = _split()
    A = M - S if order == "implicit_first" else -S + M
    assert A.is_invertible is False
    assert isinstance(A, OperatorSum)  # the static narrow: the sum's own inverse()
    with pytest.raises(NotInvertible, match=r"GreenOperator\(M, N\)"):
        A.inverse()


@pytest.mark.verifies("green-split-inverse")
@pytest.mark.rests_on(_HERE + "::test_c27_a_sum_of_two_terms_has_no_inverse_of_its_own")
def test_c28_green_is_built_from_an_explicit_split():
    """``GreenOperator(M, N)`` returns ``(M − N)⁻¹q`` (dense reference,
    SAFETY × tol), its forward is ``M − N``, and it refuses a non-invertible
    implicit part naming its type. First red today: the constructor takes a
    sum. Mutation: the forward spelled ``M + N`` (GREEN_SIGN)."""
    from orpheus.numerics.green_operator import GreenOperator

    M, S, F = _split()
    N = S + F
    tol = 1e-12
    G = GreenOperator(M, N, tol=tol)
    q = _rng(62).standard_normal(4)
    ref = np.linalg.solve(M.M - S.M - F.M, q)
    np.testing.assert_allclose(G.apply(q), ref, rtol=10 * tol * np.linalg.cond(M.M - S.M - F.M))
    assert G.inverse() == M - N
    with pytest.raises(NotInvertible, match="Mat"):
        GreenOperator(S, M)


@pytest.mark.verifies("combination-inverse")
@pytest.mark.rests_on(_HERE + "::test_c23_arrow_keys_separate_their_inners", _HERE + "::test_c12_one_term_collapse")
def test_c29_one_term_inverse_is_algebra_closed(ABC):
    """``(2A)⁻¹ = ½A⁻¹``: equal as combinations (the inverse arrow keyed by
    its inner) and bit for bit in apply. First red today: two distinct
    ``ScaledOperator`` objects. Mutation: ARROW_ID."""
    A, _, _ = ABC
    y = _rng(63).standard_normal(4)
    inv = (2.0 * A).inverse()
    assert inv == 0.5 * A.inverse()
    np.testing.assert_array_equal(inv.apply(y), 0.5 * A.inverse().apply(y))


@pytest.mark.verifies("combination-inverse")
@pytest.mark.rests_on(_HERE + "::test_c29_one_term_inverse_is_algebra_closed")
def test_c30_a_scaled_inverse_scales_its_seed():
    """Q4: ``(αA)⁻¹`` seeds its inner solve with ``α·x₀``, because the inner
    answer ``A⁻¹b`` is α times the outer one. A latent defect: today the seed
    is forwarded unscaled (first red, ``[x₀]`` recorded). Mutation: SEED_RAW."""
    A = _SeededMat(np.diag([2.0, 3.0, 4.0, 5.0]), V)
    x0 = np.array([1.0, -2.0, 0.5, 4.0])
    scaled = 4.0 * A
    assert isinstance(scaled, OperatorSum)  # the one-term combination (static narrow)
    # The static face of ``OperatorSum.inverse()`` is ``LinearOperator``, whose
    # ``apply`` declares no seed; the one-term combination forwards it at
    # runtime (the S2 spec §1 finding 6). Typed as ``Any`` here, honestly.
    inverse: Any = scaled.inverse()
    inverse.apply(np.ones(4), initial_guess=x0)
    np.testing.assert_array_equal(A.recorder.seeds[-1], 4.0 * x0)


# ═══════════════════════════════════════════════════════════════════════
# 8. Roles, and the retired idioms refuse
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.rests_on(_HERE + "::test_c01_sum_commutes")
def test_c31_roles_fold_over_the_terms_whatever_the_order():
    """``block_role`` of ``L + C`` with roles FULL and BULK is FULL in both
    orders and with a third term. Green today for the two-term case; red
    under ROLE_HEAD (role from the first term)."""
    L = Mat(np.eye(4), V, V, block_role=BlockRole.FULL)
    C = Mat(2 * np.eye(4), V, V, block_role=BlockRole.BULK)
    X = Mat(3 * np.eye(4), V, V, block_role=BlockRole.BULK)
    assert (L + C).block_role is BlockRole.FULL
    assert (C + L).block_role is BlockRole.FULL
    assert (C + X - L).block_role is BlockRole.FULL


def test_c32_the_binary_idioms_are_gone(ABC):
    """The combination has no ``.a``/``.b``; ``ScaledOperator`` is not in the
    operator module; ``2A`` is an ``OperatorSum``. Pattern 4's corollary:
    the replacement refuses the old idiom. First red today."""
    import orpheus.numerics.operator as opm

    A, B, _ = ABC
    assert not hasattr(A + B, "a") and not hasattr(A + B, "b")
    assert not hasattr(A @ B, "a")
    assert not hasattr(opm, "ScaledOperator")
    assert type(2.0 * A) is OperatorSum


# ═══════════════════════════════════════════════════════════════════════
# 9. `+` never fuses; the fused composites are atoms (D4, ruling 52)
# ═══════════════════════════════════════════════════════════════════════


def _sn_legs():
    """``(L, C, S, B_lower, sn)`` on the heterogeneous 2G vacuum slab of the
    SN Green gates (its config notes say why: het σ_t, asymmetric 2G P0)."""
    from orpheus.sn.operators.boundary import SNBoundaryOperator
    from orpheus.sn.operators.streaming import StreamingOperator
    from orpheus.transport.operators.multiplication_operator import MultiplicationOperator
    from orpheus.transport.operators.scattering import ScatteringOperator
    from tests.gates.sn.operators.test_green_operator_sn import _het_scattering_slab

    sn = _het_scattering_slab()
    L = StreamingOperator.pose(sn)
    C = MultiplicationOperator.from_mesh(sn.mat_xs.total_cross_section, sn)
    S = ScatteringOperator.from_solver_data(
        mat_xs=sn.mat_xs, scattering_order=0, space=sn.full_field_space,
    )
    return L, C, S, SNBoundaryOperator(sn), sn


@pytest.mark.verifies("fused-atom-expansion")
@pytest.mark.rests_on(_HERE + "::test_c01_sum_commutes", _HERE + "::test_c27_a_sum_of_two_terms_has_no_inverse_of_its_own")
def test_c33_streaming_plus_collision_is_a_plain_order_free_sum():
    """``L + C`` and ``C + L`` are the same plain combination; neither is the
    fused atom, and neither has an inverse of its own. First red today:
    ``L + C`` fuses and ``C + L`` does not. Mutation: restore the
    ``StreamingOperator.__add__`` fused branch."""
    from orpheus.sn.operators.streaming import StreamingCollisionOperator

    L, C, _, _, _ = _sn_legs()
    assert type(L + C) is OperatorSum and type(C + L) is OperatorSum
    assert L + C == C + L
    assert not isinstance(L + C, StreamingCollisionOperator)
    assert (L + C).is_invertible is False


@pytest.mark.verifies("fused-atom-expansion")
@pytest.mark.rests_on(_HERE + "::test_c33_streaming_plus_collision_is_a_plain_order_free_sum")
def test_c34_the_fused_atom_is_built_explicitly_and_expands_to_its_algebra():
    """``StreamingCollisionOperator(L, C)`` is an atom: ``expansion() == L + C``
    (and ``== C + L``), its legs are the named fields, its inverse is the
    sweep, and in a further sum it is ONE term. First red today: no
    ``expansion``. Mutation: the atom's ``expansion`` built from ``C + L``
    alone is not a mutation (equal by D1); the tooth is a fused atom left an
    ``OperatorSum`` subclass (the ``sco − S`` row then flattens it)."""
    from orpheus.sn.operators.streaming import StreamingCollisionOperator
    from orpheus.sn.operators.sweep_operator import SweepOperator

    L, C, S, _, _ = _sn_legs()
    sco = StreamingCollisionOperator(L, C)
    assert not isinstance(sco, OperatorSum)
    assert sco.expansion() == L + C == C + L
    assert sco.streaming is L and sco.diagonal is C
    assert isinstance(sco.inverse(), SweepOperator)
    assert {id(w) for _, w in _terms(sco - S)} == {id(sco), id(S)}


@pytest.mark.verifies("fused-atom-expansion")
@pytest.mark.rests_on(_HERE + "::test_c34_the_fused_atom_is_built_explicitly_and_expands_to_its_algebra")
def test_c35_scheduled_and_windowed_atoms_expand_to_their_algebra():
    """``ScheduledInvertibleOperator(sco, B_lower).expansion() == L + C −
    B_lower`` (flattened over the original atoms) and ``WindowedSweep(p,
    sweep).expansion() == p @ sweep`` (a plain word), on the production
    Gauss-Seidel windowed config; ``p @ sweep`` itself no longer fuses, so the
    windowed atom is built explicitly (by the shared helper, as its owner,
    the solver's ``AngularMomentWindow``, does). First red today: no ``expansion``; ``@`` fuses."""
    from orpheus.sn.operators.scheduled_invertible import ScheduledInvertibleOperator
    from orpheus.sn.operators.windowing import BulkAnalysisOperator, WindowedSweep
    from tests.gates.sn.solve.test_2d_anisotropic_windowing import (
        _windowed_product_and_oracle_operands,
    )

    _, base, product, _, _ = _windowed_product_and_oracle_operands(inner_schedule="gauss_seidel")
    assert isinstance(base, ScheduledInvertibleOperator)
    sco = base.invertible
    assert base.expansion() == sco.streaming + sco.diagonal - base.lower
    # The helper builds the windowed atom explicitly (its owner's spelling).
    assert isinstance(product, WindowedSweep)
    assert not isinstance(product, OperatorProduct)
    p, sweep = product.p, product.sweep
    assert isinstance(p, BulkAnalysisOperator)
    assert type(p @ sweep) is OperatorProduct  # ``@`` never fuses
    assert product.expansion() == p @ sweep


@pytest.mark.rests_on(_HERE + "::test_c35_scheduled_and_windowed_atoms_expand_to_their_algebra", _HERE + "::test_c13_apply_is_the_left_fold_of_todays_left_nested_spelling")
@pytest.mark.verifies("fused-atom-expansion")
def test_c40_the_scheduled_atom_applies_its_expansion():
    r"""``ScheduledInvertibleOperator.apply`` (the sweep atom's own ``M(σ)ψ``
    minus ``B_lower ψ``) equals its expansion ``L + C − B_lower`` applied as a
    plain combination, to ``8ε`` relative. The two are independent spellings
    of one operator: the atom reads σ through the loss representation's walk,
    the expansion through the collision multiplier. `[M]` 2026-10-20 (the S2
    elegance review): 1.2e-16 relative on the production Gauss-Seidel config.
    Mutation: drop ``− B_lower`` from the atom's ``apply`` and this reds."""
    from orpheus.sn.operators.scheduled_invertible import ScheduledInvertibleOperator
    from tests.gates.sn.solve.test_2d_anisotropic_windowing import (
        _windowed_product_and_oracle_operands,
    )

    _, base, _, rhs, _ = _windowed_product_and_oracle_operands(inner_schedule="gauss_seidel")
    assert isinstance(base, ScheduledInvertibleOperator)
    psi = base.solve(rhs)  # a state with live boundary rows, so B_lower acts
    atom = np.asarray(base.apply(psi).interior.values)
    expanded = np.asarray(base.expansion().apply(psi).interior.values)
    boundary_atom = np.asarray(base.apply(psi).boundary.values)
    boundary_expanded = np.asarray(base.expansion().apply(psi).boundary.values)
    scale = max(float(np.max(np.abs(atom))), float(np.max(np.abs(boundary_atom))))
    assert scale > 0.0
    assert float(np.max(np.abs(atom - expanded))) <= 8.0 * np.finfo(float).eps * scale
    assert float(np.max(np.abs(boundary_atom - boundary_expanded))) <= 8.0 * np.finfo(float).eps * scale
    assert float(np.max(np.abs(np.asarray(base.lower.apply(psi).boundary.values)))) > 0.0, (
        "B_lower does not act on this state, so the row cannot see its term"
    )


@pytest.mark.rests_on(_HERE + "::test_c17_a_scalar_leaves_the_word", _HERE + "::test_c11_every_zero_is_absorbed_and_contributes_its_spaces")
@pytest.mark.verifies("word-normal-form")
def test_c41_a_one_term_combination_keeps_its_ends_inside_a_word():
    """A one-term combination gives its coefficient to the word only when its
    word carries the combination's ends. ``I + Z_W`` is the one-term ``1·I``
    bound to ``W`` by the absorbed zero; composed after ``A`` on ``V`` it must
    be refused like ``A @ Z_W``, never unwrapped to the spaceless ``I`` that
    passes. Mutation: unwrap every one-term combination, and the composition
    is accepted with domain ``V`` (the S2 elegance review's probe)."""
    from orpheus.numerics.operator import IncompatibleOperatorComposition

    a = Mat(np.eye(3), V1, V1)
    w_zero = ZeroMorphism(domain=VD, codomain=VD)
    bound_identity = IdentityOperator() + w_zero
    assert isinstance(bound_identity, OperatorSum) and bound_identity.domain == VD
    with pytest.raises(IncompatibleOperatorComposition, match="OperatorProduct requires each factor's domain"):
        a @ w_zero
    with pytest.raises(IncompatibleOperatorComposition, match="OperatorProduct requires each factor's domain"):
        a @ bound_identity


@pytest.mark.rests_on(_HERE + "::test_c31_roles_fold_over_the_terms_whatever_the_order", _HERE + "::test_c09_zero_is_the_bound_zero_morphism")
def test_c42_a_role_lives_on_the_leaves_and_survives_every_rebuild():
    r"""A block role stamped on a LEAF survives every combination the algebra
    rebuilds from it: ``−A``, ``αA``, ``0·A`` (the zero inherits the cancelled
    term's role), ``A − A``, and a mix of two leaves of the same role; a mix of
    two roles joins to FULL. The user's ruling of 2026-10-20 (the S2 review's
    F1): a role is a fact about an atom, and composites derive it by the join.
    Before it, a role stamped on a COMPOSITE was lost whenever the algebra
    rebuilt the composite from its terms (``-op`` of a realized boundary law
    read None at α = 0.5). Mutation: mint the zero without the role join, and
    the ``0·A`` and ``A − A`` rows red."""
    a = Mat(np.eye(4), V2, V2, block_role=BlockRole.BOUNDARY)
    b = Mat(2.0 * np.eye(4), V2, V2, block_role=BlockRole.BOUNDARY)
    bulk = Mat(3.0 * np.eye(4), V2, V2, block_role=BlockRole.BULK)
    for rebuilt in (-a, 0.5 * a, 0.0 * a, a - a, 0.3 * a + 0.7 * b, -(0.3 * a + 0.7 * b)):
        assert rebuilt.block_role is BlockRole.BOUNDARY, (rebuilt, rebuilt.block_role)
    assert (a + bulk).block_role is BlockRole.FULL


@pytest.mark.verifies("combination-left-fold")
@pytest.mark.rests_on(_HERE + "::test_c13_apply_is_the_left_fold_of_todays_left_nested_spelling", _HERE + "::test_c14_apply_does_not_alias_its_operand")
@pytest.mark.parametrize(
    "first, later",
    [
        (np.ones(3), np.ones((2, 3))),  # a later term broadcasts wider
        (np.ones(3, np.float32), np.full(3, 0.1)),  # a later term promotes the dtype
        (np.ones(3, int), np.full(3, 0.5)),  # an integer accumulator meets a float
    ],
    ids=["wider_broadcast", "float32_meets_float64", "int_meets_float"],
)
def test_c43_the_fold_value_never_depends_on_accumulating_in_place(first, later):
    """The fold accumulates in place only where that IS the out-of-place sum,
    so ``A + B + C`` applies, shape and dtype included, exactly as
    ``((Ax + Bx) + Cx)`` spelled out of place (unit coefficients, so no
    multiply promotes the accumulator before the in-place step). Found by the S2 qa review
    (F2): NumPy's ``+=`` keeps the accumulator's shape and dtype, so it raised
    on a wider later term and silently kept float32. No production path
    reaches these carriers today (`[M]` qa's fold census: 108 float64 arrays of
    matching shape, 251 typed carriers). Mutation: accumulate in place
    unconditionally, and all three cases red."""

    class Const(LinearOperator):
        def __init__(self, value):
            self.value = value

        @property
        def domain(self):
            return None

        @property
        def codomain(self):
            return None

        def apply(self, x):
            return self.value

    terms = (Const(first), Const(first), Const(later))
    combination = terms[0] + terms[1] + terms[2]
    reference = (first + first) + later
    result = combination.apply(np.zeros(3))
    assert result.dtype == reference.dtype and result.shape == reference.shape
    np.testing.assert_array_equal(result, reference)
