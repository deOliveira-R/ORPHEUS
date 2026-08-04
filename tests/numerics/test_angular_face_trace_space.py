r"""The Γ ladder — :class:`AngularFaceTraceSpace` and its defining laws.

G6.1 of the ``geometric_transformation_consolidation`` campaign (tracked
tree-wide as **#330**). These gates cover the type MINTED here; the gates that
cover BINDING operators to it land with G6.3.

Two design facts drive every gate below, and both are measured, not assumed:

1. ⭐ **The face is part of the space's IDENTITY.**
   :meth:`FunctionSpace.__eq__` is ``(name, shape)`` and
   ``inner_product_weights`` is ``compare=False``, so the metric offers NO
   secondary discrimination. Every shipped quadrature gives
   ``|Γ₊(xmin)| == |Γ₊(xmax)|`` over *different* ordinate indices — so a name
   omitting the face would make a space compare EQUAL to its opposite face's
   and let a cross-face composition type-check while being wrong.

2. ⚠ **``gauss_legendre(4)`` has ZERO tangential ordinates.** It is the
   degenerate fixture: a gate written only on it cannot see the middle tier at
   all (``vv-principles`` Mode 7 — the fixture nulls the term it was meant to
   exercise). Every partition-touching gate here is parameterised over
   ``product(4,4)`` and ``lebedev(17)`` as well, and
   :func:`test_the_degenerate_fixture_is_degenerate` pins WHY, so a future
   reader cannot quietly drop the expensive fixtures.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.face_layout import FaceLayout
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.spaces import AngularFaceTraceSpace, AngularTraceSpace

pytestmark = pytest.mark.foundation


# Lazily constructed: building a quadrature at COLLECTION time turns a
# construction failure into a collection error that takes the whole file down
# — the G5 lesson (a mutation battery reported "no output" for exactly this).
_QUADRATURES = {
    "gauss_legendre(4)": lambda: Quadrature.gauss_legendre(n_ordinates=4),
    "product(4,4)": lambda: Quadrature.product(n_mu=4, n_phi=4),
    "lebedev(17)": lambda: Quadrature.lebedev(order=17),
}

#: The fixtures that can SEE the tangential tier. ``gauss_legendre(4)`` is
#: excluded by measurement, not by taste — see the module docstring.
_TANGENTIAL_BEARING = ["product(4,4)", "lebedev(17)"]

_ALL = list(_QUADRATURES)


def _trace(name: str, ng: int = 2) -> AngularTraceSpace:
    """A slab trace space with the PRODUCTION per-face shape ``(N, ng)``.

    Deliberately ``ng > 1``: the metric is a 1-D vector along the leading
    (ordinate) axis and is re-expanded across the trailing group axis on
    application, so a single-group fixture would leave that broadcast
    untested — and the broadcast is exactly what makes the space's shape
    ownership work. Mirrors ``SNMesh.boundary_face_layout``'s slab arm.
    """
    quad = _QUADRATURES[name]()
    layout = FaceLayout.from_named_shapes(
        [("xmin", (int(quad.N), ng)), ("xmax", (int(quad.N), ng))]
    )
    return AngularTraceSpace.from_quadrature_and_layout(quad, layout)


# ─────────────────────────────────────────────────────────────────────
# Identity — the reason the face is in the name
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("quad_name", _ALL)
def test_same_role_on_different_faces_is_a_different_space(quad_name):
    """``Γ₊(xmin) != Γ₊(xmax)`` — the collision this type exists to refuse.

    The NEGATIVE leg is the point: these two have EQUAL shape on every shipped
    quadrature, so shape alone cannot separate them. If they compared equal, a
    boundary law realized at ``xmin`` would type-check when composed against
    ``xmax``'s trace — silently wrong, and invisible to every numerical gate.
    """
    trace = _trace(quad_name)
    xmin, xmax = trace.outflow_space("xmin"), trace.outflow_space("xmax")

    # The positive leg: the collision is REAL, not hypothetical.
    assert xmin.shape == xmax.shape, (
        "fixture no longer exercises the collision — the two faces have "
        "different shapes, so this gate would pass for the wrong reason"
    )
    # The negative leg: identity separates them anyway.
    assert xmin != xmax
    assert hash(xmin) != hash(xmax)


@pytest.mark.parametrize("quad_name", _ALL)
def test_opposite_roles_on_one_face_are_different_spaces(quad_name):
    """``Γ₊(f) != Γ₋(f)`` — the role is load-bearing too.

    Same argument one axis over: the two half-traces have equal shape on every
    shipped quadrature (measured 2/2, 4/4, 49/49) over DISJOINT index sets.
    """
    trace = _trace(quad_name)
    out, inn = trace.outflow_space("xmin"), trace.inflow_space("xmin")

    assert out.shape == inn.shape, "fixture no longer exercises the collision"
    assert out != inn
    assert set(out.ordinate_indices.tolist()).isdisjoint(
        inn.ordinate_indices.tolist()
    )


@pytest.mark.parametrize("quad_name", _ALL)
def test_a_tier_equals_itself_across_calls(quad_name):
    """The positive leg of identity — and that the ladder is CACHED.

    Without this, the two refusals above are satisfiable by a type that simply
    never compares equal to anything, which would be useless.
    """
    trace = _trace(quad_name)
    first, second = trace.outflow_space("xmin"), trace.outflow_space("xmin")

    assert first == second
    assert hash(first) == hash(second)
    # Cached, not merely equal: a half-trace space is per-(face, quadrature),
    # and rebuilding one per call would move allocation into the sweep path.
    assert first is second


# ─────────────────────────────────────────────────────────────────────
# Structure — the three tiers do NOT partition
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("quad_name", _ALL)
def test_the_halves_are_disjoint_subsets_of_the_full_tier(quad_name):
    """``Γ₊ ⊔ Γ₋ ⊆ Γ(f)``, always — containment and disjointness."""
    trace = _trace(quad_name)
    full = set(trace.face_space("xmin").ordinate_indices.tolist())
    out = set(trace.outflow_space("xmin").ordinate_indices.tolist())
    inn = set(trace.inflow_space("xmin").ordinate_indices.tolist())

    assert out.isdisjoint(inn)
    assert out | inn <= full


@pytest.mark.parametrize("quad_name", _TANGENTIAL_BEARING)
def test_the_full_tier_is_strictly_larger_than_the_two_halves(quad_name):
    """⭐ ``Γ(f) ⊋ Γ₊ ⊔ Γ₋`` — why the middle tier needs to exist.

    If the full slot were the direct sum of the halves, it would be
    derivable and would not earn a tier. It is not: the tangential ordinates
    (``|Ω·n| <= TANGENTIAL_EPS``) belong to neither half. This is the code-side
    statement of the THREE-way face partition, and it is the reason
    ``"not inflow"`` must never be spelled ``"outflow"``.
    """
    trace = _trace(quad_name)
    full = trace.face_space("xmin")
    out = trace.outflow_space("xmin")
    inn = trace.inflow_space("xmin")

    n_tangential = full.shape[0] - out.shape[0] - inn.shape[0]
    assert n_tangential > 0, (
        f"{quad_name} has no tangential ordinates, so it cannot exercise the "
        "middle tier — it does not belong in _TANGENTIAL_BEARING"
    )


def test_the_degenerate_fixture_is_degenerate():
    """⚠ Pins WHY ``gauss_legendre(4)`` is excluded from the tier gates.

    This gate exists so the exclusion above is a measured fact in the suite
    rather than a comment someone later "simplifies" away by re-adding the
    cheap fixture to ``_TANGENTIAL_BEARING`` and getting a green run that
    proves nothing.
    """
    trace = _trace("gauss_legendre(4)")
    full = trace.face_space("xmin")
    out = trace.outflow_space("xmin")
    inn = trace.inflow_space("xmin")

    assert full.shape[0] - out.shape[0] - inn.shape[0] == 0
    assert "gauss_legendre(4)" not in _TANGENTIAL_BEARING


@pytest.mark.parametrize("quad_name", _ALL)
@pytest.mark.parametrize("role", ["full", "outflow", "inflow"])
def test_indices_are_sorted_and_unique(quad_name, role):
    """The index set is a SET, presented canonically.

    ``TraceRestrictionOperator`` validates sortedness and uniqueness for
    hand-built index sets; the ladder must satisfy those guards by
    construction, not by luck.
    """
    trace = _trace(quad_name)
    space = {
        "full": trace.face_space,
        "outflow": trace.outflow_space,
        "inflow": trace.inflow_space,
    }[role]("xmin")
    idx = space.ordinate_indices

    assert np.array_equal(idx, np.sort(idx))
    assert len(set(idx.tolist())) == idx.size
    assert space.shape[0] == idx.size


# ─────────────────────────────────────────────────────────────────────
# The metric — a trace pairing is PHYSICAL, never Euclidean
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("quad_name", _ALL)
@pytest.mark.parametrize("role", ["outflow", "inflow"])
def test_the_half_trace_metric_is_strictly_positive(quad_name, role):
    """On a HALF, ``G = |Ω·n|·w > 0`` — the tangential rows are excluded.

    This is what makes the halves a genuine inner-product space (rather than
    merely semi-definite), and it is the precondition for the round-trip law
    below to be an identity rather than a projection.
    """
    trace = _trace(quad_name)
    space = (trace.outflow_space if role == "outflow" else trace.inflow_space)("xmin")
    weights = np.asarray(space.inner_product_weights)

    assert weights is not None
    assert weights.shape == (space.shape[0],)
    assert (weights > 0.0).all()


@pytest.mark.parametrize("quad_name", _ALL)
def test_the_metric_is_not_euclidean(quad_name):
    """⭐ The metric MATTERS — dropping it is not a harmless simplification.

    If ``G`` were constant, ``A† = G⁻¹AᵀG`` would collapse to ``Aᵀ`` for every
    operator and the whole binding exercise would be decoration. Measured
    max/min: 1.35 (gauss_legendre 4), 3.47 (product 4,4), 5.6 (lebedev 17).
    """
    trace = _trace(quad_name)
    weights = np.asarray(trace.inflow_space("xmin").inner_product_weights)

    assert not np.allclose(weights, weights[0]), (
        "the trace metric is constant on this fixture, so it cannot "
        "distinguish the Hilbert adjoint from the Euclidean transpose"
    )


@pytest.mark.parametrize("quad_name", _ALL)
@pytest.mark.parametrize("role", ["outflow", "inflow"])
def test_the_metric_round_trip_is_the_identity_on_a_half(quad_name, role):
    r"""``G⁺(G x) == x`` where ``G > 0`` — the adjoint's building block.

    :class:`_AdjointOperator` computes ``A† = G_V⁻¹AᵀG_W``, so this round-trip
    is the law that makes ``A†† == A``. Asserted at 2 nulp, not
    ``array_equal``: ``(g·x)/g`` is a genuine floating-point round trip and is
    NOT bit-identical. The tolerance is derived from the round-trip depth (one
    multiply, one divide), not from folklore — do not relax it to ``allclose``.
    """
    trace = _trace(quad_name)
    space = (trace.outflow_space if role == "outflow" else trace.inflow_space)("xmin")
    rng = np.random.default_rng(20260804)
    x = rng.standard_normal(space.shape)

    np.testing.assert_array_almost_equal_nulp(
        space.apply_inverse_metric(space.apply_metric(x)), x, nulp=2,
    )


@pytest.mark.parametrize("quad_name", _TANGENTIAL_BEARING)
def test_the_full_tier_round_trip_is_a_projection_not_the_identity(quad_name):
    r"""On ``Γ(f)`` the round trip is ``G⁺G = P``, the tangential projector.

    The NEGATIVE counterpart of the gate above, and it can only be written on
    a tangential-bearing fixture. ``|Ω·n| = 0`` on the tangential rows, so the
    Moore-Penrose pseudo-inverse correctly annihilates them rather than
    dividing by zero — the components carry zero ``⟨·,·⟩_G`` weight, so this
    is exact for the Hilbert adjoint, not an approximation.
    """
    trace = _trace(quad_name)
    space = trace.face_space("xmin")
    weights = np.asarray(space.inner_product_weights)
    tangential = weights == 0.0
    assert tangential.any(), "fixture cannot exercise the null space"

    rng = np.random.default_rng(20260804)
    x = rng.standard_normal(space.shape)
    round_tripped = space.apply_inverse_metric(space.apply_metric(x))

    assert (round_tripped[tangential] == 0.0).all()
    np.testing.assert_array_almost_equal_nulp(
        round_tripped[~tangential], x[~tangential], nulp=2,
    )


# ─────────────────────────────────────────────────────────────────────
# The inner product's defining laws (symmetry / bilinearity / definiteness)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("quad_name", _ALL)
def test_the_inner_product_obeys_its_defining_laws(quad_name):
    """A math-bearing type ships a test of the laws it claims, not just usage.

    Symmetry, linearity in the first argument, positive-definiteness, and the
    induced-norm relation — on ``Γ₋``, where the metric is strictly positive so
    definiteness is a genuine claim rather than semi-definiteness.
    """
    trace = _trace(quad_name)
    space = trace.inflow_space("xmin")
    rng = np.random.default_rng(20260804)
    x, y, z = (rng.standard_normal(space.shape) for _ in range(3))
    alpha = 2.5

    assert space.inner_product(x, y) == pytest.approx(space.inner_product(y, x))
    assert space.inner_product(alpha * x + y, z) == pytest.approx(
        alpha * space.inner_product(x, z) + space.inner_product(y, z)
    )
    assert space.inner_product(x, x) > 0.0
    assert space.inner_product(np.zeros(space.shape), x) == pytest.approx(0.0)
    assert space.norm(x) == pytest.approx(np.sqrt(space.inner_product(x, x)))


@pytest.mark.parametrize("quad_name", _ALL)
def test_the_weighted_pairing_differs_from_the_euclidean_one(quad_name):
    """The metric is not merely present — it CHANGES the answer.

    A metric that existed but happened to agree with the Euclidean pairing
    would leave every adjoint gate unable to tell "applied" from "dropped".
    This is the activation check for the whole G6 argument, and it is written
    to be reddenable: it fails if the two ever coincide.
    """
    trace = _trace(quad_name)
    space = trace.inflow_space("xmin")
    rng = np.random.default_rng(20260804)
    x, y = rng.standard_normal(space.shape), rng.standard_normal(space.shape)

    assert not np.isclose(
        space.inner_product(x, y), float(np.sum(x * y)), rtol=1e-9,
    )


# ─────────────────────────────────────────────────────────────────────
# Construction discipline
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("quad_name", _ALL)
def test_an_unknown_face_raises_naming_the_available_faces(quad_name):
    """The ladder inherits the layout's own error, not a KeyError."""
    trace = _trace(quad_name)
    for accessor in (trace.face_space, trace.outflow_space, trace.inflow_space):
        with pytest.raises(ValueError, match="Unknown face"):
            accessor("not_a_face")


@pytest.mark.parametrize("quad_name", _ALL)
def test_the_name_encodes_face_and_role(quad_name):
    """The spelling is the contract — ``angular_trace[<face>:<role>]``.

    Pinned because identity rides entirely on this string: a refactor that
    "tidied" the name into something face-blind would re-open the collision
    the whole type exists to close, and every other gate here would still pass
    if the two faces' names merely differed in some other way.
    """
    trace = _trace(quad_name)
    for face in ("xmin", "xmax"):
        assert trace.face_space(face).name == f"angular_trace[{face}:full]"
        assert trace.outflow_space(face).name == f"angular_trace[{face}:outflow]"
        assert trace.inflow_space(face).name == f"angular_trace[{face}:inflow]"


@pytest.mark.parametrize("quad_name", _ALL)
def test_a_tier_carries_its_face_and_role_as_data(quad_name):
    """``face`` / ``role`` are readable, so a diagnostic need not parse the name."""
    trace = _trace(quad_name)
    space = trace.outflow_space("xmax")

    assert space.face == "xmax"
    assert space.role == "outflow"
    assert isinstance(space, AngularFaceTraceSpace)
