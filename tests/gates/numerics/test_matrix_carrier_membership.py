r"""A matrix carrier bound to spaces acts on MEMBERS (#683, option 1).

The membership law of U1b (``1c1d8d1c``): an element of a space takes exactly
the space's shape. A carrier holding a matrix :math:`M` realises a map
:math:`V \to W` in the spaces' C-order flat charts, so

* ``apply`` parses its operand at :math:`V`'s action verb
  (:meth:`~orpheus.numerics.space.FunctionSpace.leading_coordinates_of`, a
  trailing payload :math:`P` riding through as :math:`M \otimes I_P`) and
  returns an array of shape ``W.shape + payload``;
* ``apply_transpose`` is the same definition with the ends swapped;
* the matrix is admitted at construction only when it is
  ``(W.n_coordinates, V.n_coordinates)``; an unbound end is
  :math:`\mathbb R^{n}`, whose members are 1-D.

The independent reference on every value row is numpy on the chosen dense
matrix (``M @ x.ravel()`` reshaped by the test, ``np.linalg.solve``), and the
metric adjoint's reference is :math:`G_V^{-1} M^{\mathsf T} G_W` built from the
axis weights the test chose, never read back from the space.

**Activation.** The non-uniform weights sit on the multi-point axes (a rank-1
point axis has a scalar weight, which commutes with everything: Mode 12, space
side). Two rectangular fixtures: EQUAL ``n_coordinates`` with different shapes
(``(2, 1)`` vs ``(2,)``, the defect a size check cannot see, only the returned
shape) and DIFFERENT ``n_coordinates`` (``(3, 2)`` to ``(4,)``, so a swap of
ends cannot even reshape). Every matrix is asymmetric, so a transpose dropped
or doubled moves values.

**First red** (the existing-tree input): every row below is red at ``HEAD``
``56b54f76`` (``apply`` returns the flat CSR matvec and admits any matrix);
the in-process arms of the battery under ``scratch/683/battery/`` redden them
on the carved tree. The spec is ``scratch/683/gate_spec.md``.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy import sparse

from orpheus.numerics.assembled_operator import SparseAssembledOperator
from orpheus.numerics.axis import Axis, BasisKind, Declared, EnergyAxis
from orpheus.numerics.coordinate_space import NotInCoordinateSpace
from orpheus.numerics.kernel_projection_gauge import KernelProjectionGauge, KernelSummand
from orpheus.numerics.matrix_inverse_operator import MatrixInverseOperator
from orpheus.numerics.operator import LinearOperator
from orpheus.numerics.space import FunctionSpace

pytestmark = [pytest.mark.foundation]

_THIS = "tests/gates/numerics/test_matrix_carrier_membership.py"
_CSR_ROWS = (
    "tests/gates/numerics/test_assembled_operator.py::TestSparseAssembledOperator::test_apply_is_the_csr_matvec",
    "tests/gates/numerics/test_assembled_operator.py::TestSparseAssembledOperator::test_apply_transpose_is_the_exact_transpose",
)
#: The adjoint reference is a dense chain of three small products; the pairing
#: identity is two independently accumulated sums of six terms. ~10 ulp each.
_RTOL = 1e-14


# ── Fixtures: spaces built from weights the test chose ─────────────────────

#: The quotient point, the homogeneous pose's spelling (a one-cell counting axis).
_POINT = Axis("spatial", (1,), origin=Declared(BasisKind.NODAL))
_W_GROUP = np.array([1.0, 3.0])
_W_A = np.array([1.0, 3.0, 0.5])
_W_B = np.array([2.0, 0.25])
_W_C = np.array([1.0, 4.0, 0.5, 2.0])


def _axis(name: str, weights: np.ndarray) -> Axis:
    return Axis(name, (weights.shape[0],), weights, origin=Declared(BasisKind.NODAL))


def energy_point() -> FunctionSpace:
    """Energy ⊗ point under the counting metric, shape ``(2, 1)`` (#683's own fixture)."""
    return FunctionSpace.of_axes(EnergyAxis.synthetic(2), _POINT)


def weighted_group_point() -> FunctionSpace:
    """A weighted group axis ⊗ point, shape ``(2, 1)``, Gram ``diag(1, 3)``."""
    return FunctionSpace.of_axes(_axis("weighted-group", _W_GROUP), _POINT)


def weighted_group() -> FunctionSpace:
    """The same weighted group axis alone, shape ``(2,)``: equal ``n_coordinates``, different shape."""
    return FunctionSpace.of_axes(_axis("weighted-group", _W_GROUP))


def v32() -> FunctionSpace:
    """Shape ``(3, 2)``, Gram ``outer(_W_A, _W_B)`` raveled."""
    return FunctionSpace.of_axes(_axis("a", _W_A), _axis("b", _W_B))


def w4() -> FunctionSpace:
    """Shape ``(4,)``, Gram ``diag(_W_C)``."""
    return FunctionSpace.of_axes(_axis("c", _W_C))


_M22 = np.array([[2.0, -1.0], [0.5, 4.0]])
_M46 = np.array([
    [2.0, 0.0, -1.0, 0.5, 0.0, 3.0],
    [0.0, 4.0, 0.0, 0.0, -2.0, 0.0],
    [0.5, 0.0, 8.0, 1.0, 0.0, 0.0],
    [0.0, 0.25, 0.0, 0.0, 1.5, -4.0],
])


def _sao(m: np.ndarray, domain: FunctionSpace | None, codomain: FunctionSpace | None) -> SparseAssembledOperator:
    return SparseAssembledOperator(sparse.csr_array(m), domain=domain, codomain=codomain)


# ── Rung 1: the forward arm returns a member of the codomain ───────────────


class TestTheForwardArmReturnsAMember:
    @pytest.mark.rests_on(*_CSR_ROWS)
    def test_energy_point_member_in_member_out(self) -> None:
        """#683's case: bound to Energy ⊗ point (2, 1), a (2, 1) member maps to
        a (2, 1) member, valued ``M @ x.ravel()``. HEAD returns (2,)."""
        V = energy_point()
        op = _sao(_M22, V, V)
        x = np.array([[0.75], [-1.5]])
        y = op.apply(x)
        assert y.shape == V.shape, f"apply returned {y.shape}, not the codomain's {V.shape}"
        np.testing.assert_array_equal(y, (_M22 @ x.ravel()).reshape(V.shape))
        V.coordinates_of(y)  # the exact-shape parse: a member, or NotInCoordinateSpace

    @pytest.mark.rests_on(f"{_THIS}::TestTheForwardArmReturnsAMember::test_energy_point_member_in_member_out")
    def test_the_operators_own_adjoint_accepts_its_output(self) -> None:
        """``op.H.apply(op.apply(x))`` composes: the output is a member of the
        codomain, which ``.H`` parses. HEAD: ``.H`` refuses the (2,) output."""
        V = energy_point()
        op = _sao(_M22, V, V)
        x = np.array([[0.75], [-1.5]])
        z = op.H.apply(op.apply(x))
        assert np.shape(z) == V.shape
        np.testing.assert_array_equal(z, (_M22.T @ _M22 @ x.ravel()).reshape(V.shape))

    @pytest.mark.rests_on(*_CSR_ROWS)
    def test_rectangular_ends_differently_shaped(self) -> None:
        """V (3, 2) to W (4,): apply takes (3, 2) to (4,), apply_transpose takes
        (4,) to (3, 2), and each refuses the other end's member. A swap of ends
        on either arm reddens (it cannot even reshape: 6 ≠ 4)."""
        V, W = v32(), w4()
        op = _sao(_M46, V, W)
        x = np.linspace(0.5, 3.0, 6).reshape(V.shape)
        y = np.array([1.0, -2.0, 0.5, 3.0])
        fwd, bwd = op.apply(x), op.apply_transpose(y)
        assert fwd.shape == W.shape and bwd.shape == V.shape, (fwd.shape, bwd.shape)
        np.testing.assert_array_equal(fwd, _M46 @ x.ravel())
        np.testing.assert_array_equal(bwd, (_M46.T @ y).reshape(V.shape))
        with pytest.raises(NotInCoordinateSpace, match=r"leading shape \(3, 2\)"):
            op.apply(y)
        with pytest.raises(NotInCoordinateSpace, match=r"leading shape \(4,\)"):
            op.apply_transpose(x)

    @pytest.mark.rests_on(*_CSR_ROWS)
    def test_equal_dimension_ends_keep_their_own_shapes(self) -> None:
        """V (2, 1) to W (2,), the same n_coordinates: the returned SHAPE is the
        only witness of which end each arm lands on (a size check is blind)."""
        V, W = weighted_group_point(), weighted_group()
        op = _sao(_M22, V, W)
        x = np.array([[0.75], [-1.5]])
        y = np.array([2.0, -0.5])
        assert op.apply(x).shape == (2,)
        assert op.apply_transpose(y).shape == (2, 1)
        with pytest.raises(NotInCoordinateSpace):
            op.apply(y)
        # A (2, 1) array handed to the (2,)-end is read by the action verb as a
        # member of W with a one-column payload, never re-shaped into W.
        assert op.apply_transpose(x).shape == (2, 1, 1)


# ── Rung 2: the metric adjoint law on the bound ends ───────────────────────


class TestTheMetricAdjoint:
    @pytest.mark.rests_on(
        f"{_THIS}::TestTheForwardArmReturnsAMember::test_rectangular_ends_differently_shaped",
        f"{_THIS}::TestTheForwardArmReturnsAMember::test_equal_dimension_ends_keep_their_own_shapes",
    )
    @pytest.mark.parametrize("ends", ["v32_to_w4", "group_point_to_group"])
    def test_pairing_identity_and_the_dense_reference(self, ends: str) -> None:
        r""":math:`\langle Ax, y\rangle_W = \langle x, A^{\mathsf H}y\rangle_V` with
        non-uniform Grams, and :math:`A^{\mathsf H}y = G_V^{-1}M^{\mathsf T}G_W y`
        against the dense reference from the CHOSEN weights."""
        V, W, M, g_v, g_w = (
            (v32(), w4(), _M46, np.outer(_W_A, _W_B).ravel(), _W_C) if ends == "v32_to_w4"
            else (weighted_group_point(), weighted_group(), _M22, _W_GROUP, _W_GROUP)
        )
        op = _sao(M, V, W)
        rng = np.random.default_rng(683)
        x = rng.standard_normal(V.shape)
        y = rng.standard_normal(W.shape)
        h = op.H.apply(y)
        assert np.shape(h) == V.shape
        reference = ((M.T @ (g_w * y.ravel())) / g_v).reshape(V.shape)
        np.testing.assert_allclose(h, reference, rtol=_RTOL, atol=0.0)
        lhs = W.inner_product(op.apply(x), y)
        rhs = V.inner_product(x, h)
        np.testing.assert_allclose(lhs, rhs, rtol=_RTOL, atol=0.0)


# ── Rung 3: a trailing payload rides through as M ⊗ I_P ────────────────────


class TestThePayloadRidesThrough:
    @pytest.mark.rests_on(f"{_THIS}::TestTheForwardArmReturnsAMember::test_rectangular_ends_differently_shaped")
    def test_both_arms_act_column_by_column(self) -> None:
        """``V.shape + (p,)`` maps to ``W.shape + (p,)``; each payload column is
        the arm applied to that column alone (the per-column loop is the reference)."""
        V, W = v32(), w4()
        op = _sao(_M46, V, W)
        p = 3
        X = np.arange(1.0, 1.0 + 6 * p).reshape(V.shape + (p,))
        Y = np.arange(-4.0, -4.0 + 4 * p).reshape(W.shape + (p,))
        fwd, bwd = op.apply(X), op.apply_transpose(Y)
        assert fwd.shape == W.shape + (p,) and bwd.shape == V.shape + (p,), (fwd.shape, bwd.shape)
        for k in range(p):
            np.testing.assert_array_equal(fwd[..., k], op.apply(X[..., k]))
            np.testing.assert_array_equal(bwd[..., k], op.apply_transpose(Y[..., k]))


# ── Rung 4: admission at construction ──────────────────────────────────────


class TestAdmission:
    def test_a_matching_matrix_is_admitted(self) -> None:
        """The positive leg: (W.n_coordinates, V.n_coordinates) = (4, 6) constructs."""
        assert _sao(_M46, v32(), w4()).shape == (4, 6)

    @pytest.mark.parametrize(
        "case, matrix, domain, codomain, needs",
        [
            ("column mismatch", np.ones((4, 5)), "v32", "w4", r"4×6"),
            ("row mismatch", np.ones((3, 6)), "v32", "w4", r"4×6"),
            ("the issue comment's 2×2 on a (4,) domain", np.eye(2), "r4", "r2", r"2×4"),
        ],
    )
    def test_a_mismatched_matrix_is_refused_naming_its_owner(self, case, matrix, domain, codomain, needs) -> None:
        spaces = {"v32": v32, "w4": w4, "r4": lambda: FunctionSpace("r4", (4,)), "r2": lambda: FunctionSpace("r2", (2,))}
        with pytest.raises(ValueError, match=rf"SparseAssembledOperator\(.*needs {needs}"):
            _sao(matrix, spaces[domain](), spaces[codomain]())


# ── Rung 5: the unbound chart is ℝⁿ ────────────────────────────────────────


class TestTheUnboundChart:
    @pytest.mark.rests_on(*_CSR_ROWS)
    def test_a_one_d_operand_and_its_payload(self) -> None:
        op = _sao(_M46, None, None)
        x = np.linspace(-1.0, 1.0, 6)
        assert op.apply(x).shape == (4,)
        X = np.stack([x, 2.0 * x], axis=-1)
        np.testing.assert_array_equal(op.apply(X), np.stack([op.apply(x), op.apply(2.0 * x)], axis=-1))
        assert op.apply_transpose(np.ones((4, 2))).shape == (6, 2)

    def test_an_operand_that_only_ravels_to_the_length_is_refused(self) -> None:
        """HEAD ravels a (3, 2) operand into a 6-column matrix; ℝ⁶'s members are 1-D."""
        op = _sao(_M46, None, None)
        with pytest.raises(NotInCoordinateSpace):
            op.apply(np.ones((3, 2)))
        with pytest.raises(NotInCoordinateSpace):
            op.apply_transpose(np.ones((2, 2)))


# ── Rung 6: the matrix inverse acts between the inner's ends ───────────────


class TestTheMatrixInverseArms:
    """``(A⁻¹)``: apply maps inner.codomain to inner.domain, apply_transpose the
    reverse. HEAD reshapes both arms to the DOMAIN's basis shape, so on (2, 1)
    vs (2,) its transpose arm returns (2, 1) (#683's MIO defect)."""

    @staticmethod
    def _inverse() -> tuple[MatrixInverseOperator, FunctionSpace, FunctionSpace]:
        V, W = weighted_group_point(), weighted_group()
        return MatrixInverseOperator(_sao(_M22, V, W)), V, W

    @pytest.mark.rests_on(f"{_THIS}::TestTheForwardArmReturnsAMember::test_equal_dimension_ends_keep_their_own_shapes")
    def test_apply_maps_the_codomain_to_the_domain(self) -> None:
        inv, V, W = self._inverse()
        b = np.array([2.0, -0.5])
        out = inv.apply(b)
        assert np.shape(out) == V.shape
        np.testing.assert_allclose(out, np.linalg.solve(_M22, b).reshape(V.shape), rtol=_RTOL, atol=0.0)
        # A (2, 1) operand is a member of W = (2,) with a one-column payload.
        assert np.shape(inv.apply(np.array([[2.0], [-0.5]]))) == V.shape + (1,)

    @pytest.mark.rests_on(f"{_THIS}::TestTheForwardArmReturnsAMember::test_equal_dimension_ends_keep_their_own_shapes")
    def test_apply_transpose_maps_the_domain_to_the_codomain(self) -> None:
        inv, V, W = self._inverse()
        c = np.array([[0.75], [-1.5]])
        out = inv.apply_transpose(c)
        assert np.shape(out) == W.shape, f"apply_transpose returned {np.shape(out)}, not inner.codomain's {W.shape}"
        np.testing.assert_allclose(out, np.linalg.solve(_M22.T, c.ravel()), rtol=_RTOL, atol=0.0)
        with pytest.raises(NotInCoordinateSpace):
            inv.apply_transpose(np.array([0.75, -1.5]))

    @pytest.mark.rests_on(
        f"{_THIS}::TestTheMatrixInverseArms::test_apply_maps_the_codomain_to_the_domain",
        f"{_THIS}::TestTheMatrixInverseArms::test_apply_transpose_maps_the_domain_to_the_codomain",
    )
    def test_payload_and_the_metric_adjoint_of_the_inverse(self) -> None:
        r"""Both arms carry a payload column by column, and
        :math:`\langle A^{-1}b, c\rangle_V = \langle b, (A^{-1})^{\mathsf H}c\rangle_W`."""
        inv, V, W = self._inverse()
        B = np.array([[2.0, 1.0, -3.0], [-0.5, 4.0, 0.25]])
        out = inv.apply(B)
        assert np.shape(out) == V.shape + (3,)
        for k in range(3):
            np.testing.assert_allclose(out[..., k], inv.apply(B[:, k]), rtol=0.0, atol=0.0)
        b, c = B[:, 0], np.array([[0.75], [-1.5]])
        np.testing.assert_allclose(
            V.inner_product(inv.apply(b), c), W.inner_product(b, inv.H.apply(c)), rtol=_RTOL, atol=0.0,
        )


# ── Rung 7: a carrier is posed on both ends or on neither ──────────────────


class _HalfPosed(LinearOperator):
    """A dense action posed on ONE end: the inner a half-posed MatrixInverseOperator would invert."""

    def __init__(self, matrix: np.ndarray, domain: FunctionSpace | None, codomain: FunctionSpace | None) -> None:
        self.matrix, self._ends = np.asarray(matrix, dtype=float), (domain, codomain)

    @property
    def domain(self) -> FunctionSpace | None:
        return self._ends[0]

    @property
    def codomain(self) -> FunctionSpace | None:
        return self._ends[1]

    def apply(self, x, /):
        return self.matrix @ np.asarray(x, dtype=float).ravel()


class TestHalfPosedIsRefused:
    """``admit_matrix`` parses the ends: exactly one ``None`` is refused, on both
    matrix carriers, whichever end is missing (the adjoint is posed the same way)."""

    @pytest.mark.parametrize("missing", ["codomain", "domain"])
    def test_the_sparse_carrier(self, missing: str) -> None:
        V = weighted_group()
        ends = {"domain": None, "codomain": V} if missing == "domain" else {"domain": V, "codomain": None}
        with pytest.raises(ValueError, match=r"SparseAssembledOperator\(.*posed on both ends or on neither"):
            _sao(_M22, **ends)

    @pytest.mark.parametrize("missing", ["codomain", "domain"])
    def test_the_matrix_inverse(self, missing: str) -> None:
        V = weighted_group()
        inner = _HalfPosed(_M22, None, V) if missing == "domain" else _HalfPosed(_M22, V, None)
        with pytest.raises(ValueError, match=r"MatrixInverseOperator.*posed on both ends or on neither"):
            MatrixInverseOperator(inner, basis_shape=(2,))


# ── Rung 8: an unposed inverse acts on its basis carrier ───────────────────


class TestTheUnposedInverseChart:
    @pytest.mark.rests_on(f"{_THIS}::TestTheMatrixInverseArms::test_apply_maps_the_codomain_to_the_domain")
    def test_basis_shape_in_and_out_on_both_arms(self) -> None:
        """Inner unposed, ``basis_shape=(2, 3)``: the chart is ℝ^(2,3) on both
        ends, so a ``(2, 3)`` operand returns ``(2, 3)`` on each arm, valued by
        ``np.linalg.solve``, and a flat ``(6,)`` operand is refused."""
        rng = np.random.default_rng(6830)
        M = rng.standard_normal((6, 6)) + 6.0 * np.eye(6)
        inv = MatrixInverseOperator(_HalfPosed(M, None, None), basis_shape=(2, 3))
        b = rng.standard_normal((2, 3))
        fwd, bwd = inv.apply(b), inv.apply_transpose(b)
        assert np.shape(fwd) == (2, 3) and np.shape(bwd) == (2, 3), (np.shape(fwd), np.shape(bwd))
        np.testing.assert_allclose(fwd, np.linalg.solve(M, b.ravel()).reshape(2, 3), rtol=1e-13, atol=1e-15)
        np.testing.assert_allclose(bwd, np.linalg.solve(M.T, b.ravel()).reshape(2, 3), rtol=1e-13, atol=1e-15)
        with pytest.raises(NotInCoordinateSpace):
            inv.apply(b.ravel())
        with pytest.raises(NotInCoordinateSpace):
            inv.apply_transpose(b.ravel())


# ── Rung 9: the gauge's displacement is a scalar verb: an exact parse ──────


class TestTheDisplacementParsesExactly:
    def test_a_payload_operand_is_refused_by_the_parse(self) -> None:
        """``displacement`` has no payload: ``(2, 3, 4)`` on a ``(2, 3)`` space is
        refused by ``coordinates_of`` ("exact shape"), never by a numpy broadcast."""
        plane = FunctionSpace.of_axes(_axis("p", np.array([1.0, 2.0])), _axis("q", np.array([1.0, 1.0, 2.0])))
        g = KernelProjectionGauge.direct_sum(plane, (KernelSummand(np.array([1, 2, 4]), np.array([[1.0], [-1.0], [2.0]]), "plane"),))
        x = np.random.default_rng(29).standard_normal((2, 3))
        assert 0.0 < g.displacement(x) <= 1.0  # the positive leg: a member is read
        with pytest.raises(NotInCoordinateSpace, match="exact shape"):
            g.displacement(np.stack([x, x, x, x], axis=-1))


# ── Rung 10: a rank-2 payload is folded, never batched by the action ───────


class TestARankTwoPayloadIsFolded:
    """A payload ``(p, q)`` acts as ``M ⊗ I_{pq}``: the action sees the folded
    block ``(n, p·q)``. Handed ``(n, p, q)`` unfolded, scipy's ``lu_solve``
    batches the leading axes and answers wrongly (qa, 2026-11-02: max relative
    error 0.875), and a CSR matmul refuses the rank-3 operand."""

    @pytest.mark.rests_on(f"{_THIS}::TestThePayloadRidesThrough::test_both_arms_act_column_by_column")
    def test_both_carriers_both_arms(self) -> None:
        V, W = weighted_group_point(), weighted_group()
        p, q = 2, 3
        rng = np.random.default_rng(6831)
        X = rng.standard_normal(V.shape + (p, q))
        Y = rng.standard_normal(W.shape + (p, q))
        op = _sao(_M22, V, W)
        fwd, bwd = op.apply(X), op.apply_transpose(Y)
        np.testing.assert_array_equal(fwd, (_M22 @ X.reshape(2, p * q)).reshape(W.shape + (p, q)))
        np.testing.assert_array_equal(bwd, (_M22.T @ Y.reshape(2, p * q)).reshape(V.shape + (p, q)))
        inv = MatrixInverseOperator(op)
        sol, solT = inv.apply(Y), inv.apply_transpose(X)
        np.testing.assert_allclose(sol, np.linalg.solve(_M22, Y.reshape(2, p * q)).reshape(V.shape + (p, q)), rtol=_RTOL, atol=1e-15)
        np.testing.assert_allclose(solT, np.linalg.solve(_M22.T, X.reshape(2, p * q)).reshape(W.shape + (p, q)), rtol=_RTOL, atol=1e-15)
