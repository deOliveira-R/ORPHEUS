r"""S1 gates: the marginalisation lift with the axis owning its unit (ruling 36).

The gate spec is ``scratch/numerics_untangling/packages/S1_gate_spec.md``
(test-architect, 2026-10-19); the design is rulings 36 to 40 of
``.claude/plans/numerics_untangling.md``. The lift is
:mod:`orpheus.numerics.axis_marginal`; the unit is
:attr:`orpheus.numerics.axis.Axis.unit`.

The lift over axis k of V = V_k (x) V_rest, with G = G_k (x) G_rest:
u = the coordinates of the constant function 1 in the axis's representation,
m = G_k u (the covector of integration), R = m^T (x) id, pi* = R^H = u (x) id,
E = (u / <m,u>) (x) id.

Families (one class each): G1 the unit is the constant; G2 integration; G3
reciprocity; G4 the section laws; G5 bit identity with the pre-S1 nodal pair;
G6 the refusals. ⛔ G3 and G4 hold for EVERY nonzero u (the algebra is
self-consistent whatever the unit), so they are blind to a wrong unit by
construction: only G1 and G2 see it (battery arm ``unit_e0``).

Every assertion is a function call (``python -O`` safe).
"""
from __future__ import annotations

import math

import numpy as np
import numpy.testing as npt
import pytest

from orpheus.numerics import frame_faces
from orpheus.numerics.axis import Declared, Modal
from orpheus.numerics.axis_marginal import AxisPullbackOperator, AxisRetractionOperator, AxisSectionOperator
from tests.gates.numerics._axis_marginal_lift_support import (
    COEFFICIENT_HEAD_ROWS,
    FRAME_HEAD_ROWS,
    MODAL_ROWS,
    MOMENT_SPACE_ROWS,
    NODAL_ROWS,
    legendre_table,
    real_sh_table,
    resolve,
    sphere_rule,
)

pytestmark = pytest.mark.foundation

_EPS = np.finfo(float).eps
_HERE = "tests/gates/numerics/test_axis_marginal_lift.py::"
_G1C = _HERE + "TestG1UnitIsTheConstant::test_g1c_a_moment_head_unit_is_the_moments_of_one"
_G1B = _HERE + "TestG1UnitIsTheConstant::test_g1b_a_coefficient_head_unit_is_e0_and_evaluates_to_one"
_G5 = _HERE + "TestG5NodalBitIdentity::test_g5_the_nodal_arrows_are_the_pre_s1_pair"
_G3_PRODUCT = _HERE + "TestG3Reciprocity::test_g3_a_non_product_metric_cannot_be_built_on_axes"
_ALL_ROWS = NODAL_ROWS + MODAL_ROWS + [("folded(4,4) frame L4", "spatial"), ("xy LS4 moment L2", "spatial")]


def _ids(rows):
    return [f"{n}|{lab or 'head'}".replace(" ", "_") for n, lab in rows]


def _rng(seed):
    return np.random.default_rng(seed)


def _axis_block(V, k):
    """G_k from the RAW axis or its positioned form (never from the operator)."""
    ax = V.axes[k]
    nk = int(np.prod(ax.shape))
    form = V._overlay_forms()[k]
    if form is not None:
        return np.asarray(form.matrix, float)
    if ax.weights is None:
        return np.eye(nk)
    return np.diag(np.asarray(ax.weights, float).ravel())


def _dims(V, k):
    start = sum(len(a.shape) for a in V.axes[:k])
    return tuple(range(start, start + len(V.axes[k].shape)))


def _front(x, dims):
    moved = np.moveaxis(np.asarray(x), dims, range(len(dims)))
    return moved.reshape(-1, *moved.shape[len(dims):])


# ── the independent side: a head's own frame data read as INPUT, its basis re-evaluated HERE ──


def _independent_table(frame, L):
    """(N, n_slots) table of the frame's basis at its measure's nodes, by numpy/scipy,
    σ-even restricted on a folded rule (parity derived HERE from the evaluator)."""
    from orpheus.numerics.basis.legendre_basis import LegendreBasis

    nodes = np.asarray(frame.measure.nodes)
    if isinstance(frame.basis, LegendreBasis):
        mu = frame.basis.domain.quotient_map(nodes) if nodes.ndim == 2 and nodes.shape[1] == 3 else nodes.reshape(-1)
        return legendre_table(L, mu)
    table = real_sh_table(L, nodes).reshape(nodes.shape[0], -1)
    mirror = getattr(frame.basis, "mirror_axis", None)
    if mirror is not None:
        r = _rng(5).standard_normal((40, 3))
        r /= np.linalg.norm(r, axis=1)[:, None]
        rs = r.copy()
        rs[:, mirror] *= -1.0
        even = np.all(np.isclose(real_sh_table(L, r), real_sh_table(L, rs), atol=1e-12), axis=0).ravel()
        table = table * even[None, :]
    return table


def _head_frame(V):
    from orpheus.numerics.axis import FromFrame

    match V.axes[0].origin:
        case FromFrame(frame):
            return frame
        case origin:
            pytest.fail(f"precondition: {V.axes[0].label!r} is not frame-generated (origin {type(origin).__name__})")


def _head_L(V):
    return V.axes[0].shape[0] - 1


# ═══ G1 — the unit is the constant ═══════════════════════════════════════════


class TestG1UnitIsTheConstant:
    @pytest.mark.parametrize("name,label", NODAL_ROWS, ids=_ids(NODAL_ROWS))
    def test_g1a_a_nodal_unit_is_ones(self, name, label):
        """THEOREM: on point values the constant 1 is the all-ones vector. ``array_equal``.
        First red (before S1, at 0457f113): AttributeError (``Axis`` has no ``unit``). Arm ``unit_pinv``
        (u := G⁺m) reds the zero-weight row only."""
        V, label = resolve(name, label)
        ax = V.axis(label)
        npt.assert_array_equal(np.asarray(ax.unit), np.ones(ax.shape))

    def test_g1d_the_ld_moment_unit_is_e0_and_its_marginal_is_the_cell_average(self):
        """THEOREM: the LD scheme's within-cell moments are Legendre coefficients on the
        reference cell (P_0 = 1 at slot 0), so the unit of the spatial-moment axis is e_0,
        and with its mass diagonal [1, 1/3] the covector is e_0 too: the marginal over the
        axis is the average slot, the cell average. The second leg reads slot 0 of the
        array directly, never through the lift. Added after the S1 qa review (F1): the arm
        ``ld_unit_ones`` (u := ones) reddened nothing in 524 tests."""
        from orpheus.numerics.moment_layout import SPATIAL_MOMENT_AXIS_LABEL

        V, label = resolve("slab GL8 LD angular_trial", SPATIAL_MOMENT_AXIS_LABEL)
        ax = V.axis(label)
        e0 = np.zeros(ax.shape)
        e0[0] = 1.0
        npt.assert_array_equal(np.asarray(ax.unit), e0)
        assert V.axes is not None  # narrowing only
        k = V._axis_index(label)
        dim = sum(len(a.shape) for a in V.axes[:k])
        x = _rng(17).standard_normal(V.shape)
        npt.assert_array_equal(AxisRetractionOperator.of(V, label).apply(x), np.take(x, 0, axis=dim))

    @pytest.mark.parametrize("name,label", COEFFICIENT_HEAD_ROWS, ids=_ids(COEFFICIENT_HEAD_ROWS))
    def test_g1b_a_coefficient_head_unit_is_e0_and_evaluates_to_one(self, name, label):
        """THEOREM (no-prefactor convention: P_0 = Y_00 = 1): u = e_0 exactly, and the
        function it represents, evaluated by numpy/scipy at 64 random points, is 1.
        First red (before S1, at 0457f113): AttributeError. Arm ``sh_orthonormal`` reds."""
        V, label = resolve(name, label)
        ax = V.axis(label)
        u = np.asarray(ax.unit)
        e0 = np.zeros(ax.shape)
        e0.flat[0] = 1.0
        npt.assert_array_equal(u, e0)
        r = _rng(11).standard_normal((64, 3))
        r /= np.linalg.norm(r, axis=1)[:, None]
        L = ax.shape[0] - 1
        if len(ax.shape) == 1:  # Legendre about the spent axis
            spent = {"x": 0, "y": 1, "z": 2}[ax.spent_axis]
            f = legendre_table(L, r[:, spent]) @ u
        else:
            f = real_sh_table(L, r).reshape(64, -1) @ u.ravel()
        npt.assert_allclose(f, 1.0, rtol=0, atol=64 * _EPS)

    @pytest.mark.parametrize("name,label", MOMENT_SPACE_ROWS + FRAME_HEAD_ROWS, ids=_ids(MOMENT_SPACE_ROWS + FRAME_HEAD_ROWS))
    def test_g1c_a_moment_head_unit_is_the_moments_of_one(self, name, label):
        """REFERENCE: a frame-dressed head carries COVARIANT moments φ_i = Σ_n w_n Y_i(Ω_n) ψ_n
        (its metric is the inverse discrete Gram), so the coordinates of the constant are the
        moments of 1, computed HERE with an independent evaluator and ``math.fsum``; and the
        function they represent, Y·(G_k u), is 1 at the axis's nodes.
        [M] 2026-10-19 on the model: ≤ 5.7e-16 / ≤ 4.4e-16; with u = e_0 (the design's
        literal prescription) the second leg misses by 0.50 (Legendre) and 0.92-0.95 (harmonic).
        First red (before S1, at 0457f113): AttributeError. Arm ``unit_e0`` reds every row."""
        V, label = resolve(name, label)
        frame, L = _head_frame(V), _head_L(V)
        u = np.asarray(V.axes[0].unit).ravel()
        w = np.asarray(frame.measure.weights, float)
        Y = _independent_table(frame, L)
        truth = np.array([math.fsum(w * Y[:, i]) for i in range(Y.shape[1])])
        scale = max(1.0, float(np.max(np.abs(truth))))
        npt.assert_allclose(u, truth, rtol=0, atol=w.size * _EPS * scale)
        Gk = _axis_block(V, 0)
        npt.assert_allclose(Y @ (Gk @ u), 1.0, rtol=0, atol=64 * _EPS)


# ═══ G2 — integration ═════════════════════════════════════════════════════════


@pytest.mark.verifies("spaces-axis-marginal-integral")
class TestG2Integration:
    @pytest.mark.rests_on(_G1B)
    @pytest.mark.parametrize("name,label", COEFFICIENT_HEAD_ROWS, ids=_ids(COEFFICIENT_HEAD_ROWS))
    def test_g2a_a_coefficient_head_integrates_the_represented_function(self, name, label):
        """REFERENCE: R(c ⊗ y) = (∫_{S²} f_c dΩ)·y, the integral by a 40 × 80 product rule
        built HERE and f_c evaluated by numpy/scipy (never the axis's own measure).
        [M] model: rel ≤ 1.5e-15. First red (before S1, at 0457f113): TypeError (MODAL refused)."""
        V, label = resolve(name, label)
        ax = V.axis(label)
        c = _rng(21).standard_normal(ax.shape)
        y = _rng(22).standard_normal(V.without_axis(label).shape)
        dirs, wq = sphere_rule()
        L = ax.shape[0] - 1
        if len(ax.shape) == 1:
            f = legendre_table(L, dirs[:, {"x": 0, "y": 1, "z": 2}[ax.spent_axis]]) @ c
        else:
            f = real_sh_table(L, dirs).reshape(dirs.shape[0], -1) @ c.ravel()
        integral = math.fsum(wq * f)
        got = np.asarray(AxisRetractionOperator.of(V, label).apply(np.multiply.outer(c, y)))
        npt.assert_allclose(got, integral * y, rtol=1e-13, atol=1e-13 * abs(integral) * np.max(np.abs(y)))

    @pytest.mark.rests_on(_G1C, _G5)
    @pytest.mark.parametrize("name,label", MOMENT_SPACE_ROWS, ids=_ids(MOMENT_SPACE_ROWS))
    def test_g2b_the_moment_marginal_is_the_ordinate_marginal(self, name, label):
        """REFERENCE (a commuting square): integrating the moments of ψ over the harmonic
        head equals integrating ψ over the ordinates, R_head(Mψ) = R_angular(ψ), with Mψ
        computed HERE (independent evaluator, σ-even restriction). The right side rests on
        G5 (the nodal retraction pinned to the hand einsum). [M] model: max rel ≤ 2.2e-16;
        arm ``unit_e0`` reds by 0.50 (GL) / 0.92 (harmonic). First red (before S1, at 0457f113): TypeError."""
        from tests.gates.numerics._axis_marginal_lift_support import problems

        V, label = resolve(name, label)
        problem = problems()[name.rsplit(" moment", 1)[0]]
        A = problem.angular_bulk_space
        frame, L = _head_frame(V), _head_L(V)
        psi = _rng(31).standard_normal(A.shape)
        Y = _independent_table(frame, L)
        w = np.asarray(frame.measure.weights, float)
        moments = np.einsum("n,ni,n...->i...", w, Y, psi).reshape(V.shape)
        lhs = np.asarray(AxisRetractionOperator.of(V, label).apply(moments))
        rhs = np.asarray(AxisRetractionOperator.of(A, "angular").apply(psi))
        npt.assert_allclose(lhs, rhs, rtol=0, atol=w.size * Y.shape[1] * _EPS * float(np.max(np.abs(rhs))))

    @pytest.mark.rests_on(_G1B, _G1C)
    @pytest.mark.parametrize("name,label", MODAL_ROWS, ids=_ids(MODAL_ROWS))
    def test_g2c_the_mass_of_a_modal_head(self, name, label):
        """THEOREM: ⟨m, u⟩ is the mass of the measure the head carries: 4π on a continuum
        head (the sphere, for Legendre too: its measure is the sphere's pushed to μ), Σw on
        a frame-dressed head. [M] model: exact to 1 ULP. Arm ``unit_e0`` reds (0.5, 1/4π)."""
        V, label = resolve(name, label)
        ax = V.axis(label)
        E = AxisSectionOperator.of(V, label)
        from orpheus.numerics.axis import FromFrame

        origin = ax.origin
        mass = math.fsum(np.asarray(origin.frame.measure.weights, float)) if isinstance(origin, FromFrame) else 4.0 * math.pi
        npt.assert_allclose(E.mass, mass, rtol=8 * _EPS)


# ═══ G3 — reciprocity under the product metric ════════════════════════════════


def _product_residual(V, label, x):
    k = V._axis_index(label)
    dims = _dims(V, k)
    Gk = _axis_block(V, k)
    M = V.without_axis(label)
    full = np.asarray(V.apply_metric(x))
    fx = _front(x, dims)
    contracted = np.einsum("ij,j...->i...", Gk, fx)
    rebuilt = np.stack([np.asarray(M.apply_metric(s)) for s in contracted])
    rebuilt = np.moveaxis(rebuilt.reshape(*V.axes[k].shape, *M.shape), range(len(dims)), dims)
    return float(np.max(np.abs(full - rebuilt)) / max(np.max(np.abs(full)), 1e-300))


class TestG3Reciprocity:
    @pytest.mark.rests_on(_G3_PRODUCT)
    @pytest.mark.parametrize("name,label", _ALL_ROWS, ids=_ids(_ALL_ROWS))
    def test_g3_the_pullback_is_the_hilbert_adjoint(self, name, label):
        """THEOREM: G = G_k ⊗ G_rest (leg 1, from the raw axes and forms), so
        ⟨Rψ, φ⟩ = ⟨ψ, π*φ⟩ (leg 2, residual normalised by ‖Gψ‖‖π*φ‖ ≤ 16ε; [M] model
        ≤ 0.06ε). ⛔ Blind to the unit (arm ``unit_e0`` green, declared). Arms
        ``pullback_from_m``, ``adjoint_is_section``, ``pullback_scaled``, ``metric_inverse`` red.
        First red (before S1, at 0457f113): TypeError (modal rows), ValueError "two measures" (D2)."""
        V, label = resolve(name, label)
        x = _rng(41).standard_normal(V.shape)
        if _product_residual(V, label, x) > 4 * _EPS:
            pytest.fail(f"{name}: the metric is not the product over {label!r}; the closed form is not licensed")
        R = AxisRetractionOperator.of(V, label)
        psi, phi = _rng(42).standard_normal(V.shape), _rng(43).standard_normal(R.codomain.shape)
        P_phi = np.asarray(R.H.apply(phi))
        lhs = R.codomain.inner_product(np.asarray(R.apply(psi)), phi)
        rhs = V.inner_product(psi, P_phi)
        scale = np.linalg.norm(np.asarray(V.apply_metric(psi))) * np.linalg.norm(P_phi)
        if abs(lhs - rhs) > 16 * _EPS * scale:
            pytest.fail(f"{name}: reciprocity residual {abs(lhs - rhs) / scale:.2e} (normalised)")

    @pytest.mark.parametrize("name,label", MODAL_ROWS, ids=_ids(MODAL_ROWS))
    def test_g3_the_closed_form_agrees_with_the_sandwich_on_modal_heads(self, name, label):
        """X4 witness on the modal rows: the closed-form pullback u ⊗ id agrees with the
        generic sandwich ♯ ∘ Rᵀ ∘ ♭ (``AdjointOperator(R)``, built here; production does
        not route there). On a singular head (a frame head of deficient rank) the sandwich
        applies G⁺ and still lands on u ⊗ φ, because the moment unit u = B e₀ lies in
        range(G⁺). Bound 64·ε relative to the largest entry; [M] qa 2026-10-19: worst
        3.7e-16 over the 19 modal rows. The spec's §5 widening of S6.4, added after the
        S1 qa review (F4)."""
        from orpheus.numerics.operator import AdjointOperator

        V, label = resolve(name, label)
        R = AxisRetractionOperator.of(V, label)
        sandwich = AdjointOperator(R)  # pyright: ignore[reportArgumentType]
        phi = _rng(44).standard_normal(R.codomain.shape)
        closed, generic = np.asarray(R.H.apply(phi)), np.asarray(sandwich.apply(phi))
        npt.assert_allclose(closed, generic, rtol=0, atol=64 * _EPS * np.max(np.abs(closed)))

    def test_g3_a_non_product_metric_cannot_be_built_on_axes(self):
        """The product condition is a CONSTRUCTION invariant: an axis-built space with a
        full-shape metric object is refused, so no non-product space reaches the lift.
        [M] 2026-10-19: 0 of 51 axis-built census spaces fail the product (p1_census.py);
        the non-product spaces of the tree (FullFieldSpace, 6 of 6 problems) are axes-less
        and refused by the by-label resolver."""
        from orpheus.numerics.axis import Axis, BasisKind
        from orpheus.numerics.metric import DenseMetric
        from orpheus.numerics.space import FunctionSpace

        a = Axis("angular", (2,), weights=np.array([0.5, 1.5]), origin=Declared(BasisKind.NODAL))
        b = Axis("spatial", (3,), origin=Declared(BasisKind.NODAL))
        with pytest.raises(ValueError, match="DERIVED from its axes"):
            FunctionSpace("np", (2, 3), axes=(a, b), metric=DenseMetric(np.eye(6) + 0.1))


# ═══ G4 — the section laws ════════════════════════════════════════════════════


class TestG4Section:
    @pytest.mark.parametrize("name,label", _ALL_ROWS, ids=_ids(_ALL_ROWS))
    def test_g4_retraction_of_section_is_the_identity(self, name, label):
        """THEOREM: R∘E = id at reduction depth N_k = |axis| (nulp ≤ max(4, N_k); [M] worst N/3
        ULP over 120 shipped rules, p6), and E∘R idempotent. ⛔ Blind to the unit (declared).
        Arms ``section_no_divide``, ``section_euclid_mass`` red."""
        V, label = resolve(name, label)
        R, E = AxisRetractionOperator.of(V, label), AxisSectionOperator.of(V, label)
        nk = int(np.prod(V.axis(label).shape))
        phi = _rng(51).standard_normal(R.codomain.shape)
        npt.assert_array_almost_equal_nulp(np.asarray(R.apply(E.apply(phi))), phi, nulp=max(4, nk))
        x = _rng(52).standard_normal(V.shape)
        p = np.asarray(E.apply(R.apply(x)))
        npt.assert_allclose(np.asarray(E.apply(R.apply(p))), p, rtol=0, atol=max(4, nk) * _EPS * float(np.max(np.abs(p))))


# ═══ G5 — bit identity with the pre-S1 nodal pair ═════════════════════════════

_G5_ROWS = [r for r in NODAL_ROWS if r[0] not in ("D2",)]  # D2 had no pair before S1


class TestG5NodalBitIdentity:
    @pytest.mark.parametrize("name,label", _G5_ROWS, ids=_ids(_G5_ROWS))
    def test_g5_the_nodal_arrows_are_the_pre_s1_pair(self, name, label):
        """RECORD at bit level: R = the hand einsum over the raw weights; π* = the broadcast;
        E = broadcast(φ / E.mass) (divide, then scatter). ``array_equal`` each.
        Green before S1 too (the arrows existed); S1 must keep them. Arms ``R_uses_u`` (positive control),
        ``section_scale_first``, ``pullback_from_m``, ``unit_pinv`` (zero-weight row) red."""
        V, label = resolve(name, label)
        k = V._axis_index(label)
        dims = _dims(V, k)
        ax = V.axes[k]
        w = np.ones(int(np.prod(ax.shape))) if ax.weights is None else np.asarray(ax.weights, float).ravel()
        R, E = AxisRetractionOperator.of(V, label), AxisSectionOperator.of(V, label)
        x = _rng(61).standard_normal(V.shape)
        npt.assert_array_equal(np.asarray(R.apply(x)), np.einsum("n,n...->...", w, _front(x, dims)))
        phi = _rng(62).standard_normal(R.codomain.shape) * 3.7
        bcast = np.broadcast_to(np.expand_dims(phi, dims), V.shape)
        npt.assert_array_equal(np.asarray(R.H.apply(phi)), bcast)
        scaled = np.broadcast_to(np.expand_dims(phi / E.mass, dims), V.shape)
        npt.assert_array_equal(np.asarray(E.apply(phi)), scaled)


def _rule_population():
    import itertools

    from collections.abc import Callable

    from orpheus.numerics.quadrature import Quadrature

    rules: dict[str, Callable[[], Quadrature]] = {}
    rules |= {f"GL{n}": (lambda n=n: Quadrature.gauss_legendre(n_ordinates=n)) for n in range(2, 65)}
    rules |= {f"LS{s}": (lambda s=s: Quadrature.level_symmetric(sn_order=s)) for s in (2, 4, 6, 8, 10, 12, 14, 16)}
    for a, b in itertools.product((2, 4, 6, 8), (4, 8, 12, 16)):
        rules[f"product({a},{b})"] = lambda a=a, b=b: Quadrature.product(a, b)
        rules[f"folded({a},{b})"] = lambda a=a, b=b: Quadrature.folded_product(a, b)
    rules |= {f"lebedev{o}": (lambda o=o: Quadrature.lebedev(order=o))
              for o in (3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 35, 41)}
    return rules


_RULES = _rule_population()


@pytest.mark.parametrize("rule", list(_RULES), ids=list(_RULES))
def test_g5_the_nodal_mass_is_correctly_rounded(rule):
    """THEOREM (recommended contract, spec §G5): the section's divisor ⟨m, u⟩ = Σw is the
    correctly rounded sum of the stored weights, ``math.fsum``, ``array_equal``. First red
    before S1 (0457f113): 68 of 120 rules (that divisor is the retired frame's 3-operand Gram einsum,
    up to 39 ULP off at lebedev29). Arms ``mass_gram_einsum`` (68 red), ``mass_numpy_sum``
    (29 red, `[M]` qa 2026-10-19)."""
    from orpheus.numerics.axis import Axis, BasisKind
    from orpheus.numerics.space import FunctionSpace

    w = np.asarray(_RULES[rule]().weights, float)
    V = FunctionSpace.of_axes(Axis("angular", w.shape, weights=w, origin=Declared(BasisKind.NODAL)), Axis("x", (3,), origin=Declared(BasisKind.NODAL)))
    npt.assert_array_equal(AxisSectionOperator.of(V, "angular").mass, math.fsum(w))


# ═══ G6 — the refusals keep their messages ═══════════════════════════════════

class TestG6Refusals:
    def _cells(self):
        from orpheus.numerics.axis import Axis, BasisKind

        return Axis("spatial", (5,), weights=np.array([0.2, 0.3, 0.4, 0.7, 1.4]), origin=Declared(BasisKind.NODAL))

    def test_the_energy_axis_refuses_toward_condensation(self):
        from orpheus.numerics.axis import EnergyAxis
        from orpheus.numerics.space import FunctionSpace

        V = FunctionSpace.of_axes(EnergyAxis.synthetic(2), self._cells())
        with pytest.raises(TypeError, match="condensation"):
            AxisRetractionOperator.of(V, "energy")
        with pytest.raises(TypeError, match="condensation"):
            AxisSectionOperator.of(V, "energy")

    def test_a_lone_energy_axis_names_condensation_not_the_scalar(self):
        """Wiring order: the energy clause precedes the single-axis clause."""
        from orpheus.numerics.axis import EnergyAxis
        from orpheus.numerics.space import FunctionSpace

        with pytest.raises(TypeError, match="condensation") as info:
            AxisRetractionOperator.of(FunctionSpace.of_axes(EnergyAxis.synthetic(2)), "energy")
        if "bare scalar" in str(info.value):
            pytest.fail("the energy refusal must not carry the single-axis reason")

    @pytest.mark.parametrize("which", ["nodal", "moment head"])
    def test_a_single_axis_space_refuses_as_a_bare_scalar(self, which):
        """The moment-head leg is the discrimination row: before S1 it refused as MODAL (first
        red), after S1 as a bare scalar."""
        from orpheus.numerics.axis import Axis, BasisKind
        from orpheus.numerics.quadrature import Quadrature
        from orpheus.numerics.space import FunctionSpace

        V = (FunctionSpace.of_axes(Axis("angular", (4,), weights=np.array([0.3, 0.7, 0.5, 0.5]), origin=Declared(BasisKind.NODAL)))
             if which == "nodal" else frame_faces.basis_space(frame_faces.angular_frame(Quadrature.folded_product(4, 4), 2)))
        assert V.axes is not None  # narrowing only: both members are axis-built
        with pytest.raises(ValueError, match="bare scalar"):
            AxisRetractionOperator.of(V, V.axes[0].label)

    def test_an_axes_less_space_refuses(self):
        from orpheus.numerics.space import FunctionSpace

        with pytest.raises(TypeError, match="not axis-built"):
            AxisRetractionOperator.of(FunctionSpace(name="legacy", shape=(4, 2)), "angular")

    def test_an_unknown_label_refuses_naming_the_axes(self):
        V, _ = resolve("synthetic nodal", "angular")
        with pytest.raises(ValueError, match="names 0 axes"):
            AxisRetractionOperator.of(V, "azimuthal")

    @pytest.mark.parametrize("which", ["nodal signed", "modal"])
    def test_zero_mass_refuses_the_section_only(self, which):
        """⟨m,u⟩ = 0: no section; the retraction stays legal. The modal leg is new (a head whose
        metric vanishes on the constant); first red before S1 (0457f113): TypeError (MODAL)."""
        from orpheus.numerics.axis import Axis, BasisKind, FromBasis, LegendreAxis
        from orpheus.numerics.basis.legendre_basis import LegendreBasis
        from orpheus.numerics.space import FunctionSpace

        if which == "nodal signed":
            V, label = resolve("zero-mass signed", "angular")
        else:
            head = LegendreAxis("legendre", (3,), weights=np.array([0.0, 1.0, 1.0]),
                                origin=FromBasis(LegendreBasis(L=2, axis="z"), np.array([1.0, 0.0, 0.0])),
                                spent_axis="z")
            V, label = FunctionSpace.of_axes(head, Axis("spatial", (3,), origin=Declared(BasisKind.NODAL))), "legendre"
        AxisRetractionOperator.of(V, label)
        with pytest.raises(ValueError, match="zero total weight"):
            AxisSectionOperator.of(V, label)

    @pytest.mark.parametrize("leg", ["a string kind", "a zero unit", "a unit over another shape"])
    def test_a_modal_kind_carries_its_unit_or_is_refused(self, leg):
        """Rulings 40 and 41: the kind is the closed sum ``BasisKind.NODAL |
        Modal(unit)``, so a modal factor without its unit cannot be spelled.
        What is left to refuse at construction, one leg each so one red cannot
        mask another: a kind that is neither arm (a string, the pre-S1 test
        spelling the S1 review caught), a unit that cannot be the constant's
        coordinates (zero), and a unit over another shape, here through
        ``replace`` on a nodal axis."""
        import dataclasses

        from orpheus.numerics.axis import Axis, BasisKind, Modal

        if leg == "a string kind":
            with pytest.raises(TypeError, match="kind must be BasisKind.NODAL or Modal"):
                Declared("nodal")  # pyright: ignore[reportArgumentType]
        elif leg == "a zero unit":
            with pytest.raises(ValueError, match="finite and nonzero"):
                Modal(np.zeros(2))
        else:
            nodal = Axis("moment", (2,), weights=np.array([1.0, 3.0]), origin=Declared(BasisKind.NODAL))
            with pytest.raises(ValueError, match="must live over exactly the axis shape"):
                dataclasses.replace(nodal, origin=Declared(Modal(np.eye(3)[0])))

    def test_wrong_shapes_refuse_both_directions(self):
        V, _ = resolve("synthetic nodal", "angular")
        R, E = AxisRetractionOperator.of(V, "angular"), AxisSectionOperator.of(V, "angular")
        with pytest.raises(ValueError, match="full space"):
            R.apply(np.zeros((2, 5)))
        with pytest.raises(ValueError, match="marginal space"):
            E.apply(np.zeros(V.shape))

    def test_a_second_pullback_cannot_be_minted(self):
        V, _ = resolve("synthetic nodal", "angular")
        with pytest.raises(TypeError, match="minted by its retraction, once"):
            AxisPullbackOperator(AxisRetractionOperator.of(V, "angular"))

    def test_the_constructors_memoise_one_marginal(self):
        """R.of and E.of return the SAME objects on repeat, and share ONE marginal instance;
        the consumers ``phi.space is R.codomain`` (test_angular_flux, test_angular_lift) rest
        on it."""
        V, _ = resolve("slab GL8 angular_bulk", "angular")
        R1, R2 = AxisRetractionOperator.of(V, "angular"), AxisRetractionOperator.of(V, "angular")
        E1, E2 = AxisSectionOperator.of(V, "angular"), AxisSectionOperator.of(V, "angular")
        if R1 is not R2 or E1 is not E2:
            pytest.fail("the lift must be memoised per space per axis label")
        if R1.codomain is not E1.domain or R1.H.H is not R1:
            pytest.fail("one marginal instance and R.H.H is R")
