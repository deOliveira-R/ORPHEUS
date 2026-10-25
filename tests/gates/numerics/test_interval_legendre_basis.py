r"""The interval Legendre basis and the tensor product of bases — their defining laws (S1, #423).

The spatial-moment carve (``.claude/plans/moment_basis.md``) makes the LD cell
basis a numerics object: :class:`IntervalLegendreBasis` (the Legendre
polynomials on :math:`\xi \in [-1, 1]`, :math:`P_k(1) = 1`, Gram
:math:`1/(2k+1)` against :math:`\mathrm d\xi/2`) and :class:`TensorProductBasis`
(factor 0 outermost, the Kronecker order). This module is the obligation of
``vv-testing`` "a type that embodies a mathematical concept ships the test of
its defining laws".

Claim layer and kind
--------------------
Every row is **L0, foundation, THEOREM** against a ground the code does not
share: SymPy exact integration and SymPy's own ``legendre`` (closed form), the
``numpy.polynomial.legendre`` Clenshaw evaluator (a different algorithm from
the Bonnet recurrence production calls), ``fractions.Fraction`` for the
correctly rounded Gram, and the hand-typed LD tables of the UBLD assembler
(``_FOUT_1D``, ``_FIN_TRACE``) where the face trace must reproduce them. No row
reads ``moment_layout`` or the scheme: the scheme-side rows are
``tests/gates/transport/spatial/test_cell_basis_head.py``.

Bit-identity rows (``array_equal``) are bit-exact for a stated reason in each
docstring (a correctly rounded division, an exact sign flip, an exact product
by one); every other row carries its measured tolerance.
"""

from __future__ import annotations

import pathlib
from fractions import Fraction
from functools import reduce
from itertools import product

import numpy as np
import pytest
import sympy as sp

from orpheus.geometry.manifold import Interval, Product
from orpheus.numerics.basis import IntervalLegendreBasis, TensorProductBasis, spatial_moment_basis
from orpheus.numerics.coefficient_layout import PolynomialLayout, TensorLayout

pytestmark = [pytest.mark.foundation]

_THIS = "tests/gates/numerics/test_interval_legendre_basis.py"
_XI = sp.Symbol("xi", real=True)
_COUNTS = (1, 2, 3, 4, 5)        # n = 1 (DD/Step), 2 (LD), 3..5 (the family's laws are n-general)
_DIMS = (1, 2, 3)
_SEEDS = range(10)


def _exact_gram(n: int) -> sp.Matrix:
    """The per-volume Gram by exact integration against dxi/2, SymPy's own P_k."""
    return sp.Matrix(n, n, lambda k, l: sp.integrate(sp.legendre(k, _XI) * sp.legendre(l, _XI), (_XI, -1, 1)) / 2)


# ── The interval basis ─────────────────────────────────────────────────────


@pytest.mark.parametrize("n", _COUNTS)
def test_the_gram_is_the_exact_integral_correctly_rounded(n: int) -> None:
    r"""``coefficient_weights[k] == fl(1/(2k+1))``, and the exact Gram is diagonal.

    Ground: SymPy integrates :math:`\int P_k P_l\,\mathrm d\xi/2` exactly. The
    weights are ``array_equal`` to the exact diagonal rounded ONCE, because
    ``1.0 / (2k+1)`` is one IEEE division (correctly rounded) and
    ``float(Fraction)`` is correctly rounded too. The off-diagonals are exactly
    zero, which is the ``GramStructure.DIAGONAL`` claim.
    Tooth: a Gram of ``2/(2k+1)`` (the un-halved measure) or ``1/(2k-1)``.
    """
    exact = _exact_gram(n)
    assert all(exact[k, l] == 0 for k in range(n) for l in range(n) if k != l), "the exact Gram must be diagonal"
    expected = np.array([float(Fraction(int(sp.numer(exact[k, k])), int(sp.denom(exact[k, k])))) for k in range(n)])
    np.testing.assert_array_equal(IntervalLegendreBasis(n).coefficient_weights, expected)


@pytest.mark.parametrize("n", _COUNTS)
def test_the_endpoint_values_are_plus_minus_one_to_the_degree(n: int) -> None:
    r""":math:`P_k(\pm 1) = (\pm 1)^k`, bitwise, against SymPy's closed form.

    Bit-exact: at :math:`\mu = \pm 1` every Bonnet step is
    :math:`((2k-1)\cdot(\pm1)\cdot P - (k-1)P')/k`, integers over integers.
    Tooth: an endpoint table read at the wrong end, or a parity sign dropped.
    """
    b = IntervalLegendreBasis(n)
    for side in (1, -1):
        exact = np.array([float(sp.legendre(k, side)) for k in range(n)])
        np.testing.assert_array_equal(b.endpoint_values(side), exact)
        np.testing.assert_array_equal(b.evaluate(np.array([float(side)]))[0], exact)


def test_an_endpoint_other_than_plus_or_minus_one_is_refused() -> None:
    with pytest.raises(ValueError, match="side is"):
        IntervalLegendreBasis(2).endpoint_values(0)


@pytest.mark.parametrize("n", _COUNTS)
def test_the_parity_is_the_mirrors_action(n: int) -> None:
    r""":math:`P_k(-\xi) = (-1)^k P_k(\xi)`: the stored parity, the symbolic law, and the table.

    Three legs: the stored ``parity`` equals :math:`(-1)^k`; SymPy proves the law
    symbolically; the evaluated table at :math:`-\xi` equals the table at
    :math:`\xi` times ``parity``, ``array_equal`` over 10 seeds (a sign flip of
    the input commutes exactly with every Bonnet operation).
    Tooth: ``parity = ones`` (a reflection that flips nothing).
    """
    b = IntervalLegendreBasis(n)
    np.testing.assert_array_equal(b.parity, (-1.0) ** np.arange(n))
    for k in range(n):
        assert sp.expand(sp.legendre(k, -_XI) - (-1) ** k * sp.legendre(k, _XI)) == 0
    for seed in _SEEDS:
        x = np.random.default_rng(seed).uniform(-1.0, 1.0, 257)
        np.testing.assert_array_equal(b.evaluate(-x), b.evaluate(x) * b.parity)


@pytest.mark.parametrize("n", _COUNTS)
def test_evaluation_agrees_with_clenshaw(n: int) -> None:
    r"""The Bonnet table equals ``numpy.polynomial.legendre.legval`` (Clenshaw) within ``16 * eps``.

    Two algorithms, no shared code. ``[M]`` 2026-10-25, 10 seeds x 257 uniform
    points: max ``|Bonnet - legval|`` is 0 for n <= 3, 3.5 eps at n = 4, 5.0 eps
    at n = 5 and 6.5 eps at n = 6 (the zero at n <= 3 is a draw, not a theorem,
    so it is not pinned). ``16 * eps`` leaves a factor 3.2 over n = 5. Tooth: a
    wrong recurrence coefficient moves the table O(1).
    """
    b = IntervalLegendreBasis(n)
    for seed in _SEEDS:
        x = np.random.default_rng(seed).uniform(-1.0, 1.0, 257)
        clenshaw = np.stack([np.polynomial.legendre.legval(x, np.eye(n)[k]) for k in range(n)], axis=1)
        np.testing.assert_allclose(b.evaluate(x), clenshaw, rtol=0.0, atol=16 * np.finfo(float).eps)


@pytest.mark.parametrize("n", _COUNTS)
def test_reconstruction_inverts_analysis_on_the_span(n: int) -> None:
    r"""``reconstruct(analyze(f)) == f`` at the nodes for :math:`f = \sum_k c_k P_k`, and ``analyze(f) = G c``.

    Analysis returns the COVARIANT moments :math:`\varphi_k = \int P_k f\,\mathrm d\xi/2
    = c_k/(2k+1)` (the Gram times the coefficients); reconstruction scales by
    the dual factor :math:`2k+1` and synthesises values. On a Gauss-Legendre
    rule of ``n`` nodes, weights halved (total mass 1), every product of degree
    :math:`\le 2n-2` is exact, so both legs hold to rounding. Tooth: a dual
    factor ``2k - 1`` or ``2(2k+1)``; a Gram not halved.
    """
    b = IntervalLegendreBasis(n)
    x, w = np.polynomial.legendre.leggauss(n)
    table = b.evaluate(x)
    c = np.random.default_rng(n).standard_normal(n)
    f = b.synthesize(c, table)
    moments = b.analyze(f, table, w / 2.0)
    np.testing.assert_allclose(moments, c * b.coefficient_weights, rtol=0.0, atol=1e-14)
    np.testing.assert_allclose(b.reconstruct(moments, table), f, rtol=0.0, atol=1e-14)


def test_the_count_must_be_a_positive_integer() -> None:
    for bad in (0, -1, 1.5):
        with pytest.raises(ValueError, match="positive integer"):
            IntervalLegendreBasis(bad)  # type: ignore[arg-type]


def test_a_point_off_the_reference_interval_is_refused() -> None:
    with pytest.raises(ValueError, match="reference interval"):
        IntervalLegendreBasis(2).evaluate(np.array([1.5]))


def test_the_interval_layout_is_not_an_angular_moment_head() -> None:
    r"""The refuted reuse (the plan, "The object"): a spatial axis must not pass for an angular head.

    ``LegendreLayout(1)`` and the LD 1-D layout share the shape ``(2,)``; the
    head door (:func:`~orpheus.numerics.spaces.moment_head.head_layout`) must
    refuse the spatial one. Tooth: the interval basis given a ``LegendreLayout``.
    """
    from orpheus.numerics.axis import Axis
    from orpheus.numerics.coefficient_layout import DegreeLayout
    from orpheus.numerics.space import FunctionSpace
    from orpheus.numerics.spaces.moment_head import head_layout

    for basis in (IntervalLegendreBasis(2), spatial_moment_basis(2, 1)):
        assert not isinstance(basis.layout, DegreeLayout)
        with pytest.raises(TypeError, match="is not a moment head"):
            head_layout(FunctionSpace.of_axes(Axis.of_basis(basis)))


def test_the_domain_is_the_reference_interval() -> None:
    assert IntervalLegendreBasis(2).domain == Interval(-1.0, 1.0)
    assert IntervalLegendreBasis(2).layout == PolynomialLayout(2)


# ── The tensor product ─────────────────────────────────────────────────────


@pytest.mark.parametrize("n,d", list(product((1, 2, 3), _DIMS)))
def test_the_dimension_is_the_count_to_the_dimension(n: int, d: int) -> None:
    r"""``layout.shape == (n**d,)`` and the table is ``(N, n**d)``; the domain is the d-fold product."""
    b = spatial_moment_basis(n, d)
    assert b.layout.shape == (n ** d,)
    assert b.layout.factor_shape == (n,) * d
    pts = np.random.default_rng(d).uniform(-1.0, 1.0, (7, d))
    assert b.evaluate(pts).shape == (7, n ** d)
    assert b.domain == reduce(Product, [Interval(-1.0, 1.0)] * d)


def test_the_kronecker_order_is_the_recorded_one() -> None:
    r"""d=2 is ``[ψ̄, ψ̂_y, ψ̂_x, ψ̂_xy]`` (``:label: spatial-moment-kronecker-order``): axis 0 outermost.

    At the point :math:`(\xi_x, \xi_y) = (a, b)` the table row is exactly
    :math:`[1, b, a, ab]` (a product by one is exact). The d=3 row is the
    hand-written :math:`[1, c, b, bc, a, ac, ab, abc]`, with each product
    spelled in the order ``np.kron`` forms it (the leading factor first). The
    doc leg reads the label's bracketed order on the theory page, so the page
    and the code cannot drift apart (a RECORD of the page, not a theorem).
    Tooth: factor order reversed (axis 0 innermost), which swaps slots 1 and 2.
    """
    a, b_, c = 0.3, -0.7, 0.45
    np.testing.assert_array_equal(spatial_moment_basis(2, 2).evaluate(np.array([[a, b_]]))[0], [1.0, b_, a, a * b_])
    np.testing.assert_array_equal(
        spatial_moment_basis(2, 3).evaluate(np.array([[a, b_, c]]))[0],
        [1.0, c, b_, b_ * c, a, a * c, a * b_, (a * b_) * c],
    )
    layout = spatial_moment_basis(2, 2).layout
    assert [layout.multi_index(k) for k in range(4)] == [(0, 0), (0, 1), (1, 0), (1, 1)]
    import orpheus

    page = pathlib.Path(orpheus.__file__).resolve().parents[1] / "docs/theory/methods/sn/cartesian_multid.rst"
    text = page.read_text(encoding="utf-8")
    at = text.index(":label: spatial-moment-kronecker-order")
    window = text[at: at + 600]
    for line in (r"d{=}2:&\quad [\,\bar\psi,\ \hat\psi_y,\ \hat\psi_x,\ \hat\psi_{xy}\,]",
                 r"d{=}3:&\quad [\,\bar\psi,\ \hat\psi_z,\ \hat\psi_y,\ \hat\psi_{yz},"):
        assert line in window, f"the theory page no longer records {line!r} under its label"


@pytest.mark.parametrize("n,d", list(product((1, 2, 3), _DIMS)))
def test_the_tensor_gram_is_the_exact_product_integral(n: int, d: int) -> None:
    r"""The tensor Gram equals :math:`\prod_a 1/(2k_a+1)` within ``d - 1`` roundings, and is the kron of the factors.

    Ground: ``Fraction`` products, rounded once. The kron of correctly rounded
    factors carries at most ``d - 1`` further roundings, so the per-entry band
    is ``(2d - 1) * eps/2`` relative (each of ``d`` factors rounded once, ``d - 1``
    products rounded once); at ``n = 2`` the bound is met with ZERO
    ulp for ``d <= 3`` (``[M]`` 2026-10-25, ``probe_gram_bits.out``: every
    closed-form spelling bit-identical), and the second leg pins that. At
    ``n = 3`` the two spellings differ in 1 of 9 entries at d=2 and 7 of 27 at
    d=3 (``probe_interval_laws.out``), so the contract there is the band and the
    ORDER (kron of floats, factor 0 outermost), not the once-rounded value.
    Tooth: a Gram built in the reverse factor order is invisible here (the
    factors are equal); ``test_the_kronecker_order_is_the_recorded_one`` owns order.
    """
    b = spatial_moment_basis(n, d)
    idx = np.array(np.unravel_index(np.arange(n ** d), (n,) * d)).T
    exact = np.array([float(reduce(lambda x, y: x * y, (Fraction(1, 2 * int(k) + 1) for k in o), Fraction(1))) for o in idx])
    gram = b.coefficient_weights
    np.testing.assert_allclose(gram, exact, rtol=(2 * d - 1) * np.finfo(float).eps / 2, atol=0.0)
    np.testing.assert_array_equal(gram, reduce(np.kron, [f.coefficient_weights for f in b.factors], np.array([1.0])))
    if n == 2:
        np.testing.assert_array_equal(gram, exact)


@pytest.mark.parametrize("d", _DIMS)
def test_the_constant_mode_is_slot_zero(d: int) -> None:
    r"""The unit (the coefficients of 1) is :math:`e_0`; column 0 of the table is 1; synthesis of :math:`e_0` is 1.

    This is ``AVERAGE_MOMENT = 0`` as a basis fact. Bit-exact: :math:`P_0 = 1`
    and a product of ones. Tooth: the constant mode moved to the last slot.
    """
    b = spatial_moment_basis(2, d)
    e0 = np.zeros(2 ** d)
    e0[0] = 1.0
    np.testing.assert_array_equal(b.layout.unit, e0)
    pts = np.random.default_rng(d).uniform(-1.0, 1.0, (33, d))
    table = b.evaluate(pts)
    np.testing.assert_array_equal(table[:, 0], np.ones(33))
    np.testing.assert_array_equal(b.synthesize(e0, table), np.ones(33))


@pytest.mark.parametrize("d", (2, 3))
@pytest.mark.parametrize("side", (1, -1))
def test_the_face_trace_is_one_factor_evaluated_at_its_end(d: int, side: int) -> None:
    r"""The trace of an expansion on the face :math:`\xi_a = \pm1` is the face basis applied to the contracted coefficients.

    For random coefficients :math:`c`: evaluating the cell expansion at points
    with :math:`\xi_a = \pm 1` equals ``face(a)`` evaluated at the transverse
    coordinates, applied to :math:`c` contracted on factor ``a`` with
    ``endpoint_values(±1)``. Face dimension ``2**(d-1)`` (``face_moment_count``).
    Band: the two sides sum the same products in different orders, ``[R]``
    ``2**d`` terms of magnitude ``|c|``, so ``1e-14`` absolute on unit draws.
    Tooth: the endpoint pairing taken at the other end, or the face keeping
    factor ``a`` instead of dropping it.
    """
    b = spatial_moment_basis(2, d)
    rng = np.random.default_rng(10 * d + side)
    c = rng.standard_normal(2 ** d)
    transverse = rng.uniform(-1.0, 1.0, (17, d - 1))
    for a in range(d):
        face = b.face(a)
        assert face.layout.shape == (2 ** (d - 1),)
        pts = np.insert(transverse, a, float(side), axis=1)
        cell_trace = b.evaluate(pts) @ c
        contracted = np.tensordot(b.factors[a].endpoint_values(side), c.reshape((2,) * d), axes=([0], [a])).reshape(-1)
        np.testing.assert_allclose(cell_trace, face.evaluate(transverse) @ contracted, rtol=0.0, atol=1e-14)


def test_the_face_drops_exactly_the_normal_factor() -> None:
    r"""``face(a)`` is the product of the OTHER factors, in their order, on distinguishable factors.

    On the LD cell basis every factor is the same, so which factor ``face(a)``
    drops is unobservable there (``[M]`` 2026-10-25, ``BASIS_ARM=face_keeps_wrong``
    reddened 0 rows of the face-trace law: a Mode-12 space-side stabiliser).
    This row uses factor counts 2, 3, 4 so the dropped factor is named by the result.
    """
    b = TensorProductBasis((IntervalLegendreBasis(2), IntervalLegendreBasis(3), IntervalLegendreBasis(4)), "probe")
    assert b.face(0).factors == (IntervalLegendreBasis(3), IntervalLegendreBasis(4))
    assert b.face(1).factors == (IntervalLegendreBasis(2), IntervalLegendreBasis(4))
    assert b.face(2).factors == (IntervalLegendreBasis(2), IntervalLegendreBasis(3))


def test_the_face_dimension_is_the_product_of_the_other_counts() -> None:
    r"""``face_dimension(a)`` is the product of the other factors' counts, and 1 on a one-factor product.

    Total where ``face`` is not: the face of a one-factor product is a point,
    whose function space is the constants (the empty product). On factors of
    distinguishable counts the dropped factor is named by the result.
    """
    b = TensorProductBasis((IntervalLegendreBasis(2), IntervalLegendreBasis(3), IntervalLegendreBasis(4)), "probe")
    assert [b.face_dimension(a) for a in range(3)] == [12, 8, 6]
    assert [b.face_dimension(a) for a in range(3)] == [b.face(a).layout.shape[0] for a in range(3)]
    for d in _DIMS:
        assert spatial_moment_basis(2, d).face_dimension(0) == 2 ** (d - 1)
    assert spatial_moment_basis(2, 1).face_dimension(0) == 1


def test_the_face_of_a_one_factor_product_is_refused() -> None:
    with pytest.raises(ValueError, match="is a point"):
        spatial_moment_basis(2, 1).face(0)


def test_the_endpoint_values_reproduce_the_ubld_face_tables() -> None:
    r"""``_FOUT_1D = outer(B(+1), B(+1))`` and ``_FIN_TRACE = B(-1)`` are the interval basis at its ends.

    Bit-exact: both sides are the integers ±1. Since S3 the UBLD module spells
    its tables from the basis, so the reference is the hand-typed tables they
    replaced (Larsen and Morel 1989, the LD face traces), kept here.
    """
    from orpheus.transport.spatial import _ubld

    np.testing.assert_array_equal(_ubld._FOUT_1D, np.array([[1.0, 1.0], [1.0, 1.0]]))
    np.testing.assert_array_equal(_ubld._FIN_TRACE, np.array([1.0, -1.0]))


def _hand_octant_signs(octant: tuple[int, ...]) -> np.ndarray:
    r"""The sign of slot ``[o_0, ..., o_{d-1}]`` is :math:`\prod_a \mathrm{octant}[a]^{o_a}`, decoded by hand, axis 0 outermost."""
    d = len(octant)
    signs = np.ones(2 ** d)
    for flat in range(2 ** d):
        for a in range(d):
            if (flat >> (d - 1 - a)) & 1 and octant[a] < 0:
                signs[flat] = -signs[flat]
    return signs


@pytest.mark.parametrize("d", _DIMS)
def test_the_reflection_signs_are_the_mirrors_action_and_the_octant_frame_map(d: int) -> None:
    r"""``reflection_signs(mirrored)`` is the action of :math:`\xi_a \mapsto -\xi_a`, and equals ``octant_moment_frame_signs``.

    Leg 1 (THEOREM): ``evaluate(M x) == evaluate(x) * signs`` bitwise for every
    mirror subset (sign flips are exact). Leg 2 (REFERENCE):
    ``octant_moment_frame_signs(octant, 2)``, which since S3 reads
    ``reflection_signs``, equals the hand Kronecker decode it replaced
    (:func:`_hand_octant_signs`, axis 0 outermost) for all :math:`2^d` octants.
    Tooth: the decode with axis 0 innermost (d >= 2), ``SITE_ARM=kron_decode_144``.
    """
    from orpheus.transport.spatial._ubld import octant_moment_frame_signs

    b = spatial_moment_basis(2, d)
    x = np.random.default_rng(d).uniform(-1.0, 1.0, (19, d))
    for octant in product((1, -1), repeat=d):
        mirrored = tuple(s < 0 for s in octant)
        signs = b.reflection_signs(mirrored)
        np.testing.assert_array_equal(b.evaluate(x * np.array(octant, dtype=float)), b.evaluate(x) * signs)
        np.testing.assert_array_equal(octant_moment_frame_signs(octant, 2), _hand_octant_signs(octant))


def test_the_tensor_layout_is_the_kronecker_of_its_factors() -> None:
    r"""``TensorLayout`` content: two products with the same shape and different factor counts are different layouts."""
    assert spatial_moment_basis(2, 2).layout == TensorLayout((PolynomialLayout(2),) * 2, "spatial_moment")
    assert spatial_moment_basis(4, 1).layout != spatial_moment_basis(2, 2).layout  # both (4,)
