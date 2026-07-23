r"""Core algebra of the assembly mode — the sparse carrier + the composer laws.

Stencil-assembly campaign 2b, commit 1 (the numerics home): pins the
:class:`~orpheus.numerics.assembled_operator.SparseAssembledOperator`
carrier contract, the three-layer assembly surface on
:class:`~orpheus.numerics.operator.LinearOperator` (predicate /
``SupportsAssembly`` / ``assemblable``), the composer homomorphism laws
(Sum → ``+``, Product → ``@``, Scaled → scalar ``*``), the R2
``as_matrix`` delegation WITH the retained probing pathway, and the
L16 negative test (assembly is structure-only — it must succeed with
every ``apply`` disabled).

These are STRUCTURAL algebra pins (exact float comparisons are legal
here: a CSR matvec against a basis vector extracts a stored coefficient
with no summation, and the composer laws are realized by scipy's own
sparse arithmetic on identical operands). The PHYSICS equivalence gates
(nulp / rtol, het + asymmetric-SigS fixtures) live with the emitters:
``tests/diffusion`` (the loss family) and ``tests/transport/spatial``
(DD / LD walk assembly).
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy import sparse

from orpheus.numerics.assembled_operator import SparseAssembledOperator
from orpheus.numerics.flat_operator import FlattenedOperator
from orpheus.numerics.operator import (
    IdentityOperator,
    MatrixTooLarge,
    MissingAssembly,
    OperatorProduct,
    OperatorSum,
    ScaledOperator,
    assemblable,
)

pytestmark = [pytest.mark.foundation]


# ── Fixtures: two small, deliberately ASYMMETRIC assembled leaves ──────
#
# Asymmetry (A != A^T) keeps a row/col-transposed carrier observable
# (anti-Mode-12); values are exact dyadic floats so sparse arithmetic
# reproduces dense arithmetic bit-for-bit.


def _leaf_a() -> SparseAssembledOperator:
    dense = np.array(
        [
            [2.0, 0.0, -1.0],
            [0.0, 4.0, 0.0],
            [0.5, 0.0, 8.0],
        ]
    )
    return SparseAssembledOperator(sparse.csr_array(dense))


def _leaf_b() -> SparseAssembledOperator:
    dense = np.array(
        [
            [1.0, -2.0, 0.0],
            [0.0, 0.25, 16.0],
            [0.0, 0.0, -4.0],
        ]
    )
    return SparseAssembledOperator(sparse.csr_array(dense))


# ── The carrier contract ───────────────────────────────────────────────


class TestSparseAssembledOperator:
    def test_apply_is_the_csr_matvec(self) -> None:
        op = _leaf_a()
        x = np.array([1.0, 2.0, 4.0])
        np.testing.assert_array_equal(op.apply(x), op.as_matrix() @ x)

    def test_apply_transpose_is_the_exact_transpose(self) -> None:
        op = _leaf_a()
        x = np.array([-1.0, 8.0, 0.5])
        np.testing.assert_array_equal(
            op.apply_transpose(x), op.as_matrix().T @ x
        )
        assert op.is_adjointable

    def test_coo_duplicates_are_summed(self) -> None:
        # The FEM scatter semantics: repeated (row, col) entries SUM on
        # CSR canonicalization — an emitter may scatter per-face /
        # per-cell contributions freely.
        rows = np.array([0, 0, 1])
        cols = np.array([0, 0, 1])
        vals = np.array([1.5, 2.5, -1.0])
        coo = sparse.coo_array((vals, (rows, cols)), shape=(2, 2))
        op = SparseAssembledOperator(coo)
        np.testing.assert_array_equal(
            op.as_matrix(), np.array([[4.0, 0.0], [0.0, -1.0]])
        )

    def test_assemble_is_idempotent(self) -> None:
        op = _leaf_a()
        assert op.is_assemblable
        assert op.assemble() is op  # the functor's fixed point

    def test_shape_and_input_guards(self) -> None:
        op = _leaf_a()
        assert op.shape == (3, 3)
        with pytest.raises(ValueError, match="input size 2"):
            op.apply(np.ones(2))
        with pytest.raises(ValueError, match="input size 2"):
            op.apply_transpose(np.ones(2))

    def test_as_matrix_contract(self) -> None:
        op = _leaf_a()
        # Consistent explicit basis_shape passes; a contradicting one
        # refuses; the size gate refuses densification.
        np.testing.assert_array_equal(
            op.as_matrix(basis_shape=(3,)), op.as_matrix()
        )
        with pytest.raises(ValueError, match="contradicts"):
            op.as_matrix(basis_shape=(4,))
        with pytest.raises(MatrixTooLarge):
            op.as_matrix(max_dimension=2)

    def test_not_invertible_by_default(self) -> None:
        # The matrix realization reads values, not structure: inversion
        # is an explicit consumer choice (MatrixInverseOperator), never
        # an automatic capability of the carrier.
        assert not _leaf_a().is_invertible


# ── The three-layer surface + composer homomorphism laws ──────────────


class TestAssemblyAxis:
    def test_base_default_is_false(self) -> None:
        assert not IdentityOperator().is_assemblable
        assert not assemblable(IdentityOperator())

    @pytest.mark.verifies("matrix-functor-homomorphism")
    def test_sum_law(self) -> None:
        a, b = _leaf_a(), _leaf_b()
        assembled = OperatorSum(a, b).assemble()
        np.testing.assert_array_equal(
            assembled.as_matrix(), a.as_matrix() + b.as_matrix()
        )

    def test_difference_via_the_dunder(self) -> None:
        # ``A - B`` routes through Sum(A, Scaled(-1, B)) — the composed
        # law [A] + (-1)[B], byte-equal to dense subtraction on exact
        # dyadic values.
        a, b = _leaf_a(), _leaf_b()
        assembled = (a - b).assemble()
        np.testing.assert_array_equal(
            assembled.as_matrix(), a.as_matrix() - b.as_matrix()
        )

    @pytest.mark.verifies("matrix-functor-homomorphism")
    def test_product_law(self) -> None:
        a, b = _leaf_a(), _leaf_b()
        assembled = OperatorProduct(a, b).assemble()
        np.testing.assert_array_equal(
            assembled.as_matrix(), a.as_matrix() @ b.as_matrix()
        )

    @pytest.mark.verifies("matrix-functor-homomorphism")
    def test_scaled_law(self) -> None:
        a = _leaf_a()
        assembled = ScaledOperator(-0.5, a).assemble()
        np.testing.assert_array_equal(
            assembled.as_matrix(), -0.5 * a.as_matrix()
        )

    def test_composite_predicates_recurse(self) -> None:
        a, b = _leaf_a(), _leaf_b()
        assert OperatorSum(a, b).is_assemblable
        assert OperatorProduct(a, b).is_assemblable
        assert ScaledOperator(2.0, a).is_assemblable
        # One non-emitting leg poisons the composite (the intersection
        # law, like the adjoint axis).
        assert not OperatorSum(a, IdentityOperator()).is_assemblable
        assert not OperatorProduct(IdentityOperator(), b).is_assemblable
        assert not ScaledOperator(2.0, IdentityOperator()).is_assemblable

    def test_missing_assembly_refusals(self) -> None:
        a = _leaf_a()
        with pytest.raises(MissingAssembly, match="both summands"):
            OperatorSum(a, IdentityOperator()).assemble()
        with pytest.raises(MissingAssembly, match="both factors"):
            OperatorProduct(IdentityOperator(), a).assemble()
        with pytest.raises(MissingAssembly, match="assemblable operand"):
            ScaledOperator(3.0, IdentityOperator()).assemble()

    def test_assembly_is_structure_only(self, monkeypatch) -> None:
        # L16's negative test: with every leaf ``apply`` disabled, the
        # composite still assembles — proving the emission never routes
        # through probing/apply (the whole point of the third mode).
        a, b = _leaf_a(), _leaf_b()

        def _boom(x, /):
            raise AssertionError("assembly must not call apply")

        monkeypatch.setattr(a, "apply", _boom)
        monkeypatch.setattr(b, "apply", _boom)
        assembled = OperatorSum(a, b).assemble()
        assert assembled.shape == (3, 3)


# ── R2: as_matrix delegation + the retained probing pathway ───────────


class TestAsMatrixDelegation:
    def test_delegation_matches_probing_on_the_carrier(self) -> None:
        # The core-level probed≡assembled pin: for the carrier itself a
        # probe column apply(e_j) extracts stored coefficients with no
        # summation, so the two pathways are EXACTLY equal. (The
        # physics-emitter pins at nulp live with the emitters.)
        a, b = _leaf_a(), _leaf_b()
        composed = OperatorSum(a, ScaledOperator(2.0, b))
        via_delegation = composed.as_matrix(basis_shape=(3,))
        via_probing = composed._as_matrix_by_probing((3,))
        np.testing.assert_array_equal(via_delegation, via_probing)

    def test_probing_remains_the_fallback(self) -> None:
        # A non-assemblable operator still materializes by probing —
        # the delegation is a no-op for the existing operator fleet.
        ident = IdentityOperator()
        np.testing.assert_array_equal(
            ident.as_matrix(basis_shape=(3,)), np.eye(3)
        )

    def test_delegation_honors_the_size_gate(self) -> None:
        a = _leaf_a()
        with pytest.raises(MatrixTooLarge):
            OperatorSum(a, _leaf_b()).as_matrix(
                basis_shape=(3,), max_dimension=2
            )


# ── FlattenedOperator passthrough ──────────────────────────────────────


class _FlatDuck:
    """Minimal FlatVectorLike: fixes the flat size the wrapper freezes."""

    def __init__(self, n: int) -> None:
        self._n = n

    def to_flat(self) -> np.ndarray:
        return np.zeros(self._n)

    @classmethod
    def from_flat(cls, flat, template):  # pragma: no cover - not exercised
        raise NotImplementedError


class TestFlattenedPassthrough:
    def test_passthrough_returns_the_inner_assembly(self) -> None:
        inner = _leaf_a()
        wrapped = FlattenedOperator(inner, _FlatDuck(3))
        assert wrapped.is_assemblable
        assert wrapped.assemble() is inner.assemble()

    def test_dimension_drift_refuses(self) -> None:
        inner = _leaf_a()
        wrapped = FlattenedOperator(inner, _FlatDuck(4))
        with pytest.raises(ValueError, match="drifted"):
            wrapped.assemble()

    def test_non_assemblable_inner_refuses(self) -> None:
        wrapped = FlattenedOperator(IdentityOperator(), _FlatDuck(3))
        assert not wrapped.is_assemblable
        with pytest.raises(MissingAssembly):
            wrapped.assemble()
