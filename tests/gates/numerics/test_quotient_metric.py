r"""S6b — a degenerate metric is a quotient (rulings 74 and 80): its rank and its live projector.

A metric owns its degeneracy: :meth:`~orpheus.numerics.metric.HilbertMetric.rank`
on each realisation, :attr:`~orpheus.numerics.space.FunctionSpace.rank` and
:meth:`~orpheus.numerics.space.FunctionSpace.project_live` (``P = G⁺G``, the
G-orthogonal projector onto range G) on a space. The adjoint is posed on
``V / ker G``, so its laws hold modulo the kernel and are gated in the
G-seminorm, never component-wise on the kernel. The rows on the S\ :sub:`N`
carriers (the trace ranks, ruling 80's kernel-preservation property, the sweep's
adjoint laws) are ``tests/gates/sn/operators/test_quotient_adjoint_sn.py``.

Hand counts are independent of the metric code: a diagonal metric's rank is a
count of nonzero weights written in the test; a Parseval head's rank is the rank
of the spherical harmonics sampled at the rule's nodes (``scipy.special.sph_harm_y``,
no ORPHEUS frame code), restricted to the harmonics even about the folded axis
for a folded rule.

``-O``-safe: every assertion is ``_require`` (a raise) or ``np.testing``.
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.special as sps

from orpheus.numerics import frame_faces
from orpheus.numerics.axis import Declared, Axis, BasisKind
from orpheus.numerics.axis_marginal import AxisRetractionOperator
from orpheus.numerics.metric import DenseMetric, DiagonalMetric, FactoredMetric
from orpheus.numerics.operator import AdjointOperator, RieszLowerOperator, RieszRaiseOperator, adjointable
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.space import FunctionSpace

pytestmark = [
    pytest.mark.foundation,
    pytest.mark.rests_on(
        "tests/gates/numerics/test_dual_metric.py::test_S1_the_dual_space_reads_the_inverse_metric_bitwise",
        "tests/gates/numerics/test_riesz_legs.py::test_round_trip_is_the_range_projector_on_a_singular_metric",
    ),
]

_EPS = float(np.finfo(float).eps)
_W_SING = np.array([2.0, 0.0, 3.0, 0.0])
#: rank-2 PSD Gram: one dead direction, (1, -1, 0)/√2.
_G_SING = np.array([[1.0, 1.0, 0.0], [1.0, 1.0, 0.0], [0.0, 0.0, 2.0]])
_G = np.array([[2.0, 0.5, 0.0], [0.5, 1.0, 0.25], [0.0, 0.25, 4.0]])


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def _rng(seed):
    return np.random.default_rng(seed)


# ═══════════════════════════════════════════════════════════════════════
# B1 — the rank of each metric realisation, against hand counts
# ═══════════════════════════════════════════════════════════════════════


def _metric_rows():
    dense_sing = DenseMetric(_G_SING)
    diag_sing = DiagonalMetric(_W_SING)
    diag_rank3 = DiagonalMetric(np.array([2.0, 0.0, 3.0, 1.0]))
    return [
        ("diag_singular", diag_sing, (4,), 2),
        ("diag_singular_trailing", diag_sing, (4, 3), 6),  # leading-aligned: 2 live rows × 3
        ("diag_positive", DiagonalMetric(np.array([0.5, 2.0, 4.0])), (3,), 3),
        ("dense_spd", DenseMetric(_G), (3,), 3),
        ("dense_rank2", dense_sing, (3,), 2),
        ("dense_rank2_trailing", dense_sing, (3, 2), 4),
        ("dense_inverse_of_rank2", DenseMetric.inverse_of(_G_SING), (3,), 2),
        # rank-3 diagonal block, so a product (6) and a sum (5) of block ranks differ
        # (`[M]` 2026-10-21: with two rank-2 blocks 2·2 = 2 + 2 hid a summed rank).
        ("factored_dense_x_diag", FactoredMetric((((3,), dense_sing), ((4,), diag_rank3))), (3, 4), 6),
        ("factored_euclid_x_diag", FactoredMetric((((2,), None), ((4,), diag_rank3))), (2, 4), 6),
        # a trailing axis beyond the factored blocks multiplies the rank (review (b):
        # `[M]` 2026-10-21 the pre-fix rank ignored it)
        ("factored_trailing", FactoredMetric((((3,), dense_sing), ((4,), diag_rank3))), (3, 4, 5), 30),
        ("dual_of_rank2", dense_sing.dual(), (3,), 2),
        ("dual_of_diag_singular", diag_sing.dual(), (4,), 2),
    ]


@pytest.mark.parametrize("row", _metric_rows(), ids=lambda r: r[0])
def test_B1_metric_rank_is_the_hand_count(row):
    """``m.rank(shape)`` equals the count written in the row: nonzero weights;
    the dense matrix's rank times the trailing multiplicity; a product of block
    ranks (a Euclidean block counting full); a dual's rank its primal's. First
    red: a rank counting the element size (blind to the kernel)."""
    name, m, shape, hand = row
    _require(m.rank(shape) == hand, f"[{name}] rank {m.rank(shape)} != hand count {hand}")


_Q5 = np.linalg.qr(np.random.default_rng(0).standard_normal((5, 5)))[0]


def _gram(eigs):
    G = _Q5 @ np.diag(eigs) @ _Q5.T
    return (G + G.T) / 2


def _inverted_directions(m, n):
    r"""The number of directions the metric's INVERSE face inverts, read through
    the two verbs (``apply_inverse ∘ apply`` on the unit vectors), never through
    ``rank``: the eigenvalues of that projector are 0 or 1, so the count is the
    number above ½."""
    P = np.array([m.apply_inverse(m.apply(np.eye(n)[i])) for i in range(n)]).T
    return int(np.sum(np.linalg.eigvalsh((P + P.T) / 2) > 0.5))


@pytest.mark.parametrize("rel", [1e-8, 1e-6, 1e-4, 1e-14],
                         ids=["1e-8", "1e-6", "1e-4", "below_cutoff_1e-14"])
def test_B1b_dense_rank_counts_the_directions_the_inverse_face_inverts(rel):
    """A dense metric with one eigenvalue at ``rel`` of the scale: its rank is
    the number of directions its inverse face inverts (5 above the 1e-12 cutoff,
    4 below), and the rank of ``G⁺G``. The three mid-band rows sit between the
    cutoff and 1e-3, where a rank with its own cutoff (``matrix_rank`` at rtol
    1e-6 or 1e-3) disagrees with the inverse face (qa, S6b review: neither was
    pinned). First red: a rank computed with a cutoff other than the inverse
    face's."""
    m = DenseMetric(_gram([1.0, 0.5, 0.3, 0.2, rel]))
    hand = 5 if rel > 1e-12 else 4
    _require(_inverted_directions(m, 5) == hand, f"fixture: the inverse face inverts {_inverted_directions(m, 5)}")
    _require(m.rank((5,)) == hand, f"rank {m.rank((5,))} vs the {hand} inverted directions (eigenvalue {rel:g})")
    _require(m.dual().rank((5,)) == hand, "the dual's rank differs")


def test_B1c_dense_rank_agrees_with_the_inverse_face_across_the_cutoff():
    """201 metrics with one eigenvalue at f · 1e-12, f in [0.9, 1.1] (the tie qa
    found): rank equals the inverted-direction count on every one, and both
    sides of the cutoff occur (activation). First red: ``matrix_rank`` with
    numpy's own tolerance, which split from ``pinv`` at the tie."""
    from orpheus.numerics.metric import _DENSE_METRIC_RCOND as RC

    counts, bad = set(), []
    for f in np.linspace(0.9, 1.1, 201):
        m = DenseMetric(_gram([1.0, 0.5, 0.3, 0.2, f * RC]))
        inv = _inverted_directions(m, 5)
        counts.add(inv)
        if m.rank((5,)) != inv:
            bad.append((round(float(f), 4), m.rank((5,)), inv))
    _require(bad == [], f"rank disagrees with the inverse face at {len(bad)} of 201: {bad[:5]}")
    _require(counts == {4, 5}, f"activation: the sweep did not straddle the cutoff ({counts})")


def test_B1d_inverse_of_a_gram_with_a_near_cutoff_eigenvalue():
    """``DenseMetric.inverse_of(g)`` with g's eigenvalue at 5e-13 (below the
    cutoff, where #634's Penrose refusal does not bite: `[M]` it bites from
    2e-12): rank 4, the inverted-direction count, and the dual's rank, agree."""
    m = DenseMetric.inverse_of(_gram([1.0, 0.5, 0.3, 0.2, 5e-13]))
    inv = _inverted_directions(m, 5)
    _require(inv == 4 and m.rank((5,)) == 4 and m.dual().rank((5,)) == 4,
             f"rank {m.rank((5,))}, inverted {inv}, dual {m.dual().rank((5,))}")


# ═══════════════════════════════════════════════════════════════════════
# B2 — the Parseval heads: dead slots and a rank-deficient Gram
# ═══════════════════════════════════════════════════════════════════════


def _sh_rank(q: Quadrature, L: int, even_about_fold: bool) -> int:
    r"""Rank of the complex spherical harmonics up to degree ``L`` sampled at the
    rule's nodes (rank is invariant under the unitary real/complex change of basis
    within each degree and under a rotation of the polar axis). For a folded rule
    the polar axis is the folded one and only the harmonics even under that mirror
    (``l + m`` even) are kept."""
    c = {a: np.asarray(getattr(q, a), float) for a in ("mu_x", "mu_y", "mu_z")}
    if even_about_fold:
        polar = next(a for a in c if c[a].min() >= 0.0)
    else:
        polar = "mu_z"
    rest = [a for a in c if a != polar]
    th = np.arccos(np.clip(c[polar], -1.0, 1.0))
    ph = np.arctan2(c[rest[1]], c[rest[0]])
    cols = [sps.sph_harm_y(l, m, th, ph) for l in range(L + 1) for m in range(-l, l + 1)
            if not even_about_fold or (l + m) % 2 == 0]
    s = np.linalg.svd(np.array(cols).T, compute_uv=False)
    return int(np.sum(s > 1e-10 * s[0]))


_HEADS = {
    #: name: (rule, L, folded, hand rank), `[M]` 2026-10-21 by _sh_rank
    "product48_L4": (lambda: Quadrature.product(4, 8), 4, False, 23),
    "folded26_L2": (lambda: Quadrature.folded_product(n_mu=2, n_phi=6), 2, True, 5),
    "folded48_L3": (lambda: Quadrature.folded_product(n_mu=4, n_phi=8), 3, True, 10),
}


@pytest.mark.parametrize("name", list(_HEADS))
def test_B2_a_parseval_heads_rank_is_the_sampled_harmonic_rank(name):
    """``basis_space(angular_frame(rule, L)).rank`` equals the rank of the sampled
    harmonics: product(4,8) L4, 23 of 45 layout slots (the census's 22 dead: 20
    layout zeros of the (L+1)×(2L+1) head and 2 aliased directions); the folded
    rules' mirror-even heads, 5 of 15 and 10 of 28. The hand count is pinned in
    the table AND recomputed (the reference is checked against itself first).
    Its dual has the same rank. First red: a rank reading the layout size."""
    rule, L, folded, hand = _HEADS[name]
    q = rule()
    _require(_sh_rank(q, L, folded) == hand, f"[{name}] the reference itself moved: {_sh_rank(q, L, folded)} vs {hand}")
    B = frame_faces.basis_space(frame_faces.angular_frame(q, L))
    size = int(np.prod(B.shape))
    _require(B.rank == hand and hand < size, f"[{name}] rank {B.rank} vs hand {hand} (size {size})")
    _require(B.dual().rank == B.rank, f"[{name}] the dual's rank {B.dual().rank} differs")


# ═══════════════════════════════════════════════════════════════════════
# B3 — project_live: an orthogonal projector that kills exactly ker G
# ═══════════════════════════════════════════════════════════════════════


def _spaces():
    return {
        "diag_singular": FunctionSpace("s6b_ds", (4,), inner_product_weights=_W_SING),
        "dense_rank2": FunctionSpace("s6b_dr", (3,), metric=DenseMetric(_G_SING)),
        "dense_spd": FunctionSpace("s6b_dp", (3,), metric=DenseMetric(_G)),
        "diag_positive": FunctionSpace("s6b_dpos", (3,), inner_product_weights=np.array([0.3, 1.7, 2.9])),
        "head_product48_L4": frame_faces.basis_space(frame_faces.angular_frame(Quadrature.product(4, 8), 4)),
        "head_folded26_L2": frame_faces.basis_space(
            frame_faces.angular_frame(Quadrature.folded_product(n_mu=2, n_phi=6), 2)),
    }


def _P_matrix(V):
    n = int(np.prod(V.shape))
    return np.array([np.ravel(V.project_live(np.eye(n)[i].reshape(V.shape))) for i in range(n)]).T


@pytest.mark.parametrize("name", list(_spaces()))
def test_B3_project_live_is_the_g_orthogonal_projector_onto_range_g(name):
    """``P = G⁺G``: idempotent (``P² = P``), G-self-adjoint
    (``⟨Px, y⟩_G = ⟨x, Py⟩_G``), of rank ``V.rank``, and it KILLS the kernel:
    ``G (x − Px) = 0`` (the residual is invisible to the metric). On a
    non-degenerate space it is the identity (to rounding: ``x·w/w`` is not
    bitwise ``x``). Tolerances: 16 ε relative, scaled by the stored Gram's
    conditioning on the dense heads. First red: ``P`` the identity (rank =
    size; arm ``project_identity``). A swapped leg order ``P = GG⁺`` cannot
    redden this row: for a symmetric G it is the same projector (qa F6); B4,
    the bitwise route through the legs, is that swap's only catcher."""
    V = _spaces()[name]
    Pm = _P_matrix(V)
    n = Pm.shape[0]
    tol = 1e3 * _EPS if "head" in name or "dense" in name else 16 * _EPS
    # (Not a catcher of a leg swap: for a symmetric G, GG⁺ = G⁺G, so ``P = GG⁺``
    # is the same projector; the swap is inside this row's stabiliser and only
    # B4's bitwise route sees it, qa F6.)
    np.testing.assert_allclose(Pm @ Pm, Pm, atol=tol, rtol=0)
    _require(int(np.linalg.matrix_rank(Pm, tol=1e-8)) == V.rank, f"[{name}] rank(P) != V.rank")
    x, y = _rng(1).standard_normal(V.shape), _rng(2).standard_normal(V.shape)
    lhs, rhs = V.inner_product(V.project_live(x), y), V.inner_product(x, V.project_live(y))
    _require(abs(lhs - rhs) <= tol * max(1.0, abs(lhs)), f"[{name}] P is not G-self-adjoint: {lhs!r} vs {rhs!r}")
    k = x - V.project_live(x)
    _require(np.abs(V.apply_metric(k)).max() <= tol * max(1.0, np.abs(V.apply_metric(x)).max()),
             f"[{name}] the residual x − Px is not in ker G")
    if V.rank == n:
        np.testing.assert_allclose(V.project_live(x), x, rtol=4 * _EPS, atol=0)


def test_B3b_project_live_zeroes_exactly_the_hand_kernel():
    """Diagonal: the dead coordinates are zeroed and nothing else moves beyond
    rounding; dense rank-2: the hand kernel vector (1, −1, 0)/√2 maps to 0 and the
    live vector (1, 1, 0) to itself."""
    V = FunctionSpace("s6b_k", (4,), inner_product_weights=_W_SING)
    x = np.array([1.0, 2.0, 3.0, 4.0])
    np.testing.assert_array_equal(V.project_live(x), np.where(_W_SING != 0, x, 0.0))
    D = FunctionSpace("s6b_kd", (3,), metric=DenseMetric(_G_SING))
    ker, live = np.array([1.0, -1.0, 0.0]) / np.sqrt(2.0), np.array([1.0, 1.0, 0.0])
    np.testing.assert_allclose(D.project_live(ker), 0.0, atol=8 * _EPS)
    np.testing.assert_allclose(D.project_live(live), live, rtol=8 * _EPS, atol=8 * _EPS)


@pytest.mark.parametrize("name", ["diag_singular", "dense_rank2", "head_product48_L4"])
def test_B4_the_riesz_round_trip_is_the_live_projector(name):
    """``♯ ∘ ♭ = P``, bitwise: the round trip runs the same two delegations as
    ``project_live`` (one definition of P). First red: a ``project_live`` spelled
    as its own arithmetic that drifts from the legs."""
    V = _spaces()[name]
    x = _rng(3).standard_normal(V.shape)
    np.testing.assert_array_equal(RieszRaiseOperator(V).apply(RieszLowerOperator(V).apply(x)), V.project_live(x))


# ═══════════════════════════════════════════════════════════════════════
# B5 — the retraction's closed-form adjoint on a zero weight: modulo the kernel
# ═══════════════════════════════════════════════════════════════════════


def test_B5_the_retraction_closed_form_agrees_with_the_sandwich_modulo_the_kernel():
    """On ``V = angular ⊗ spatial`` with ONE zero weight on the collapsed
    angular axis, ``R.H`` (the closed-form pullback) and the base sandwich
    ``AdjointOperator(R)`` differ by 0.500 in the full norm (the census) and by
    ~0 in the G-seminorm: `[M]` 2026-10-21 seminorm 0.0, the whole difference on
    the dead angular slot. Both legs are asserted (the full-norm gap is the
    activation: a fixture without the dead slot cannot see the quotient), plus
    reciprocity of each against the forward. First red: a closed form that
    writes into the LIVE slots (a mis-weighted pullback)."""
    spa = Axis("spatial", (5,), weights=np.array([1.0, 2.0, 0.5, 0.25, 3.0]), origin=Declared(BasisKind.NODAL))
    ang = Axis("angular", (4,), weights=np.array([0.3, 0.0, 0.5, 0.3]), origin=Declared(BasisKind.NODAL))
    V = FunctionSpace.of_axes(ang, spa)
    R = AxisRetractionOperator.of(V, "angular")
    phi = _rng(1).standard_normal(R.codomain.shape)
    if not adjointable(R):
        pytest.fail("fixture: the retraction is not adjointable")
    closed, sandwich = R.H.apply(phi), AdjointOperator(R).apply(phi)
    d = closed - sandwich
    full = float(np.linalg.norm(d) / np.linalg.norm(closed))
    semi = float(np.sqrt(max(V.inner_product(d, d), 0.0)) / np.sqrt(V.inner_product(closed, closed)))
    _require(full > 0.1, f"activation: the closed form and the sandwich agree in the full norm ({full:.2e})")
    _require(semi <= 16 * _EPS, f"the closed form differs from the sandwich in the G-seminorm: {semi:.2e}")
    x = _rng(2).standard_normal(V.shape)
    U = R.codomain if R.codomain is not None else pytest.fail("fixture: the retraction has no codomain")
    lhs = U.inner_product(R.apply(x), phi)
    for tag, img in (("closed", closed), ("sandwich", sandwich)):
        rhs = V.inner_product(x, img)
        _require(abs(lhs - rhs) <= 16 * _EPS * max(1.0, abs(lhs)), f"{tag} reciprocity {lhs!r} vs {rhs!r}")
