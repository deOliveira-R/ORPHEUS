"""MOVED from tests/gates/numerics/test_linear_combination.py (layer DAG Phase 2 refresh, specs_2r §5). Bodies unchanged."""

from __future__ import annotations
import numpy as np
import pytest
from orpheus.numerics.metric import DenseMetric
from orpheus.numerics.operator import InverseOperator, LinearOperator, NotInvertible, OperatorProduct, OperatorSum
from orpheus.numerics.space import FunctionSpace


pytestmark = pytest.mark.l0


_HERE = "tests/gates/numerics/test_linear_combination.py"   # the numerics laws these rows rest on
_SN = "tests/gates/sn/operators/test_linear_combination_sn.py"   # this file


def _terms(op) -> list[tuple[float, object]]:
    """The (coefficient, word) pairs of ``op`` in insertion order; an atom or
    a word is the one-term ``[(1.0, op)]``."""
    if isinstance(op, OperatorSum):
        return [(float(c), w) for c, w in op.terms]
    return [(1.0, op)]




def _space(name: str, weights) -> FunctionSpace:
    w = np.asarray(weights, dtype=float)
    return FunctionSpace(name, (w.size,), inner_product_weights=w)


V1 = _space("s2_v1", [1.0, 2.0, 0.5])


V2 = _space("s2_v2", [4.0, 0.25, 1.0, 2.0])


V3 = _space("s2_v3", [0.5, 1.0, 2.0, 4.0, 8.0])


V = V2


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






@pytest.fixture
def chain():
    """``A : V2 → V3``, ``B : V1 → V2``, ``C : V3 → V1`` (distinct dims)."""
    return (_mat(5, 4, 21, V2, V3), _mat(4, 3, 22, V1, V2), _mat(3, 5, 23, V3, V1))




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
@pytest.mark.rests_on(_SN + "::test_c33_streaming_plus_collision_is_a_plain_order_free_sum")
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
@pytest.mark.rests_on(_SN + "::test_c34_the_fused_atom_is_built_explicitly_and_expands_to_its_algebra")
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


@pytest.mark.rests_on(_SN + "::test_c35_scheduled_and_windowed_atoms_expand_to_their_algebra", _HERE + "::test_c13_apply_is_the_left_fold_of_todays_left_nested_spelling")
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


