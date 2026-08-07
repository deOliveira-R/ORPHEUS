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
from orpheus.numerics.operator import (
    IncompatibleOperatorComposition,
    TraceRestrictionOperator,
)
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.space import FunctionSpace
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

#: The fixtures that can build a MULTI-DIMENSIONAL layout at all —
#: ``gauss_legendre`` is a 1-D slab rule with no ``mu_y``, so a ``y``-face
#: raises rather than producing a degenerate answer.
#:
#: ⚠ Deliberately a SEPARATE list from ``_TANGENTIAL_BEARING`` even though the
#: members coincide today. They coincide for two unrelated reasons —
#: dimensionality here, tangential-ordinate count there — and a single list
#: whose values correlate with something other than its name is a field doing
#: two jobs; the day a 1-D rule grows tangential ordinates (or a 3-D rule loses
#: them) the two would need to diverge, and a merged list would silently pick
#: the wrong one.
_MULTIDIM = ["product(4,4)", "lebedev(17)"]

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
    # ⭐ The excess is not merely nonzero — it is EXACTLY the metric's null
    # space, so this gate and the degeneracy gate measure one fact from two
    # sides. See test_the_metric_null_space_IS_the_tangential_set.
    metric = np.asarray(full.inner_product_weights)
    assert int((metric == 0.0).sum()) == n_tangential


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

    ⭐ It is also the END-metric half of the factored-adjoint theorem's
    precondition: :math:`R^{*} = G_+^{-1}R^{\\mathsf T}G_-` reaches an end
    metric through its INVERSE, so ``.H`` on a half-trace is well-posed only
    because this gate holds — what reads as hygiene is what lets a half-trace
    sit at the end of an adjoint-bearing chain at all. The intermediate-tier
    half of the same precondition is
    :func:`test_the_current_space_is_ADMISSIBLE_as_a_chain_intermediate`; the
    theorem itself, including what breaks at a degenerate metric, is gated in
    ``tests/numerics/test_factored_adjoint_identity.py``.
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


# ─────────────────────────────────────────────────────────────────────
# G6.3 step 2 — S(f), the angle-INTEGRATED tier
# ─────────────────────────────────────────────────────────────────────


def _trace_2d(quad_name: str, ng: int = 2) -> AngularTraceSpace:
    """A 2-D Cartesian trace: x-faces carry ``ny``, y-faces carry ``nx``.

    Needed because ``S(f)`` must collapse the ordinate axis and leave EVERY
    trailing axis intact; a 1-D fixture (whose only trailing axis is the group)
    cannot distinguish "collapse axis 0" from "keep exactly one trailing axis".
    """
    quad = _QUADRATURES[quad_name]()
    n = int(quad.N)
    return AngularTraceSpace.from_quadrature_and_layout(
        quad,
        FaceLayout.from_named_shapes(
            [("xmin", (n, ng, 5)), ("xmax", (n, ng, 5)),
             ("ymin", (n, ng, 7)), ("ymax", (n, ng, 7))]
        ),
    )


@pytest.mark.parametrize("quad_name", _ALL)
def test_the_current_space_collapses_ONLY_the_ordinate_axis(quad_name):
    r"""``S(f)`` is ``Γ(f)`` with the ordinate axis integrated out.

    The tier the ladder was missing, and it is not a subset of ordinates — it is
    an axis COLLAPSE. The group axis must survive; its 2-D companion below adds
    the codim-1 spatial axis.
    """
    trace = _trace(quad_name, ng=2)
    for face in ("xmin", "xmax"):
        current = trace.current_space(face)
        assert current.shape == (2,)
        assert current.shape == trace.face_space(face).shape[1:]


@pytest.mark.parametrize("quad_name", _MULTIDIM)
def test_the_current_space_keeps_the_codim1_SPATIAL_axis_too(quad_name):
    r"""EVERY trailing axis survives, not just the group one.

    A 1-D fixture has exactly one trailing axis, so it cannot distinguish
    "collapse axis 0" from "keep exactly one trailing axis" — two different rules
    that agree there and diverge in 2-D, where x-faces carry ``ny`` and y-faces
    carry ``nx``. This is also the gate that would catch a hardcoded
    ``shape[:2]``-style slice.
    """
    trace = _trace_2d(quad_name)
    assert trace.current_space("xmin").shape == (2, 5)
    assert trace.current_space("ymin").shape == (2, 7)
    for face in ("xmin", "ymin"):
        assert (
            trace.current_space(face).shape
            == trace.face_space(face).shape[1:]
        )


@pytest.mark.parametrize("quad_name", _ALL)
def test_the_current_space_metric_is_EXPLICITLY_unit_not_None(quad_name):
    r"""⭐ ``None`` would encode two different states; ones encodes the intent.

    The angular measure :math:`|\Omega\cdot\hat n|\,w` is CONSUMED by the
    contraction that lands here, so nothing remains to weight — and the ladder
    carries no face area (``AngularTraceSpace``'s own metric is
    :math:`|\Omega\cdot\hat n|\,w` alone; the area-weighted boundary integral
    belongs to ``ScalarTraceSpace``). So unit is the *derived* answer, not a
    default.

    Gated as ``is not None`` because ``None`` means BOTH "deliberately
    Euclidean" and "nobody bound a metric", and no gate can separate one value
    from itself — the reason a generic space was rejected for this campaign.
    """
    trace = _trace(quad_name, ng=3)
    weights = trace.current_space("xmin").inner_product_weights

    assert weights is not None, "unit must be RECORDED, not left to the default"
    assert np.all(np.asarray(weights) == 1.0)


@pytest.mark.parametrize("quad_name", _ALL)
def test_the_metric_null_space_IS_the_tangential_set(quad_name):
    r"""⭐⭐ ``ker(G_{Γ(f)}) == Γ(f) ∖ (Γ₊ ⊔ Γ₋)`` — ONE theorem, not two facts.

    This module gates the non-partition
    (:func:`test_the_full_tier_is_strictly_larger_than_the_two_halves`) and the
    metric degeneracy
    (:func:`test_the_current_space_is_ADMISSIBLE_as_a_chain_intermediate`)
    separately, which reads as two independent coincidences. **They are the same
    fact**, and this gate is the identity that says so: the rows the halves
    exclude are *exactly* the rows the metric annihilates, because both are
    ``|Ω·n| ≤ TANGENTIAL_EPS``.

    Two consequences follow, and neither is derivable from either half alone:

    * **The full tier IS the direct sum of its halves — in the quotient.**
      ``Γ(f)/ker(G) ≅ Γ₊(f) ⊕ Γ₋(f)``. As a *Hilbert* space (the only category
      the adjoint cares about) the decomposition holds; it is only the storage
      array that carries the extra rows.
    * **Which is precisely why ``Γ(f)`` can never be a chain intermediate while
      the halves can.** The factored-adjoint theorem needs a non-degenerate
      intermediate metric, and the degeneracy IS the excess.

    Gated on every quadrature, including ``gauss_legendre(4)`` where both sides
    are empty — a vacuous-but-true instance keeps the identity honest about the
    degenerate case instead of quietly excluding it.
    """
    trace = _trace(quad_name, ng=2)
    for face in ("xmin", "xmax"):
        n_ordinates = int(trace.omega_dot_n.shape[1])
        halves = set(trace.outflow_indices_for_face(face).tolist()) | set(
            trace.inflow_indices_for_face(face).tolist()
        )
        excluded = set(range(n_ordinates)) - halves

        metric = np.asarray(trace.face_space(face).inner_product_weights)
        annihilated = set(np.flatnonzero(metric == 0.0).tolist())

        assert annihilated == excluded, (
            "the rows the halves exclude and the rows the metric annihilates "
            "have diverged; the non-partition and the degeneracy are no longer "
            "the same fact, and every argument resting on that identity needs "
            "re-deriving"
        )
        # ... and the halves are exactly the metric's support.
        assert set(np.flatnonzero(metric > 0.0).tolist()) == halves


@pytest.mark.parametrize("quad_name", _TANGENTIAL_BEARING)
def test_the_current_space_is_ADMISSIBLE_as_a_chain_intermediate(quad_name):
    r"""⭐ Non-degenerate — the factored-adjoint theorem's one requirement.

    ``R* = C*B* = G₊⁻¹RᵀG₋`` holds because the intermediate metric cancels, and
    that cancellation needs the intermediate metric to be **invertible**; at
    zero it breaks outright (see
    ``tests/numerics/test_factored_adjoint_identity.py``).

    The contrast is the point, and it is why this is parametrised over the
    tangential-bearing quadratures: ``Γ(f)``'s metric IS degenerate — exactly
    zero on tangential ordinates — so the full tier could **never** serve as a
    chain intermediate, while ``S(f)`` and the two halves can. This gate is what
    connects step 2 to the theorem rather than leaving the precondition as a
    docstring claim.

    ⭐ **Not an independent coincidence.** The degeneracy asserted here and the
    non-partition asserted in
    :func:`test_the_full_tier_is_strictly_larger_than_the_two_halves` are the
    SAME fact; :func:`test_the_metric_null_space_IS_the_tangential_set` is the
    identity that proves it. Read the three together — this one measures the
    metric, that one the index sets, and the identity says they describe one
    set of rows.
    """
    trace = _trace(quad_name, ng=2)
    current = np.asarray(trace.current_space("xmin").inner_product_weights)
    full_tier = np.asarray(trace.face_space("xmin").inner_product_weights)

    assert (current > 0.0).all(), "S(f) must be admissible as an intermediate"
    assert (full_tier == 0.0).any(), (
        "the full tier no longer has a degenerate metric, so this gate has "
        "lost the contrast it exists to draw"
    )


@pytest.mark.parametrize("quad_name", _ALL)
def test_the_current_space_has_its_own_identity(quad_name):
    """Face-encoded and distinct from every angular tier, like the rest.

    ``S(xmin)`` and ``S(xmax)`` share a shape on a symmetric mesh, so the same
    collision argument applies here as for the half-traces.
    """
    trace = _trace(quad_name, ng=2)
    xmin, xmax = trace.current_space("xmin"), trace.current_space("xmax")

    assert xmin.shape == xmax.shape, "fixture no longer exercises the collision"
    assert xmin != xmax
    assert xmin is trace.current_space("xmin")          # cached
    assert xmin != trace.face_space("xmin")
    assert xmin != trace.outflow_space("xmin")
    assert xmin != trace.inflow_space("xmin")


# ─────────────────────────────────────────────────────────────────────
# G6.3 step 1 — γ± BOUND to the ladder, and what the binding buys
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("quad_name", _ALL)
def test_the_trace_maps_are_bound_to_the_ladder(quad_name):
    r"""``γ₊ : Γ(f) → Γ₊(f)`` and ``γ₋ : Γ(f) → Γ₋(f)`` — in the TYPE system.

    B3.2 made this true in the guards and the prose; until G6.3 it was absent
    from the types, so ``LinearOperator.domain``'s own docstring applied — *"when
    either operand has ``None`` … the composability check is SKIPPED"*.
    """
    trace = _trace(quad_name)
    for face in ("xmin", "xmax"):
        assert trace.outflow_restriction(face).domain is trace.face_space(face)
        assert trace.outflow_restriction(face).codomain is trace.outflow_space(face)
        assert trace.inflow_restriction(face).domain is trace.face_space(face)
        assert trace.inflow_restriction(face).codomain is trace.inflow_space(face)


@pytest.mark.parametrize("quad_name", _ALL)
def test_the_bound_trace_map_satisfies_the_weighted_adjoint_law(quad_name):
    r"""``⟨γx, y⟩_{Γ₊} == ⟨x, γ*y⟩_{Γ(f)}`` — the point of binding.

    Unbound, ``.H`` silently degrades to the Euclidean transpose and this
    identity is asserted against nothing.
    """
    trace = _trace(quad_name)
    gamma = trace.outflow_restriction("xmin")
    domain, codomain = gamma.domain, gamma.codomain
    assert domain is not None and codomain is not None, "step 1 binds these"

    rng = np.random.default_rng(20260804)
    x = rng.standard_normal(domain.shape)
    y = rng.standard_normal(codomain.shape)

    assert codomain.inner_product(gamma.apply(x), y) == pytest.approx(
        domain.inner_product(x, gamma.H.apply(y)), rel=1e-13,
    )


@pytest.mark.parametrize("quad_name", _ALL)
def test_a_restrictions_hilbert_adjoint_IS_its_transpose_by_construction(quad_name):
    r"""⭐ For ``γ±`` the metric CANCELS — a theorem, not an oversight.

    ``Γ₊(f)``'s metric *is* ``Γ(f)``'s restricted to the same indices, so in
    ``G_dom⁻¹γᵀG_cod`` the weights cancel on the selected rows and vanish
    elsewhere. Hence binding ``γ±`` changes no numbers, and its value is the
    type-level refusal below, not a numerical correction.

    This gate exists so that fact is *recorded* rather than rediscovered: a
    reader who measures ``.H ≈ γᵀ`` and concludes the binding is inert would be
    drawing the wrong lesson. Its companion is
    :func:`test_a_restriction_bound_to_a_FOREIGN_metric_does_shift_its_adjoint`,
    which shows the metric is genuinely applied — without that negative leg this
    gate cannot distinguish "cancels" from "never ran"
    (``vv-principles`` anti-pattern #11).

    ⚠ **Asserted at 2 nulp, not bit-identity, and the bound is derived rather
    than fitted.** The cancellation is exact in exact arithmetic; in floating
    point it runs through ``(g·y)/g``, one multiply and one divide, so ≤2 ulp.
    `[M]` measured **1** nulp on ``product(4,4)`` and ``lebedev(17)``, and
    genuinely **0** on ``gauss_legendre(4)``. Do NOT relax this to ``allclose``
    — the tolerance is the claim, and a real metric error is O(1) relative, as
    the negative leg shows.
    """
    trace = _trace(quad_name)
    gamma = trace.outflow_restriction("xmin")
    domain, codomain = gamma.domain, gamma.codomain
    assert domain is not None and codomain is not None, "step 1 binds these"

    rng = np.random.default_rng(20260804)
    y = rng.standard_normal(codomain.shape)

    domain_metric = np.asarray(domain.inner_product_weights)
    codomain_metric = np.asarray(codomain.inner_product_weights)
    assert np.array_equal(codomain_metric, domain_metric[gamma.indices]), (
        "the codomain metric is no longer the restriction of the domain's, so "
        "the cancellation argument does not apply"
    )
    np.testing.assert_array_almost_equal_nulp(
        gamma.H.apply(y), gamma.apply_transpose(y), nulp=2,
    )


@pytest.mark.parametrize("quad_name", _ALL)
def test_a_restriction_bound_to_a_FOREIGN_metric_does_shift_its_adjoint(quad_name):
    """The negative leg — the metric is applied, it merely cancels.

    Rebinding the codomain to a scaled metric (no longer the restriction of the
    domain's) must move ``.H`` away from the bare transpose. If this passed
    trivially, the gate above would be pinning "the metric is ignored".

    The separation is not marginal: a ×3 metric shifts the adjoint by O(1)
    RELATIVE, against the ≤2 ulp the genuine cancellation leaves. Asserted as a
    relative floor so the two legs can never be confused for one another.
    """
    trace = _trace(quad_name)
    gamma = trace.outflow_restriction("xmin")
    codomain = gamma.codomain
    assert codomain is not None, "step 1 binds this"

    foreign = FunctionSpace(
        name="foreign", shape=codomain.shape,
        inner_product_weights=np.asarray(codomain.inner_product_weights) * 3.0,
    )
    rebound = TraceRestrictionOperator(
        gamma.indices, n_total=gamma.n_total, axis=0,
        domain=gamma.domain, codomain=foreign,
    )
    rng = np.random.default_rng(20260804)
    y = rng.standard_normal(codomain.shape)

    hilbert = rebound.H.apply(y)
    euclidean = rebound.apply_transpose(y)
    relative = np.abs(hilbert - euclidean).max() / np.abs(euclidean).max()
    assert relative > 0.1, (
        f"a x3 metric shifted the adjoint by only {relative:.2e} relative — "
        "indistinguishable from the round-trip noise the positive leg allows"
    )


@pytest.mark.parametrize("quad_name", _ALL)
def test_composing_across_faces_is_REFUSED_naming_both_spaces(quad_name):
    r"""⭐ The payoff: a cross-face composition raises instead of type-checking.

    ``|Γ₊(xmin)| == |Γ₊(xmax)|`` on every shipped quadrature, so **shape can
    never catch this** — the face had to be in the space's identity for the
    check to have anything to compare. This is the error class the whole binding
    exists to refuse, and before G6.3 the check was skipped outright.
    """
    trace = _trace(quad_name)
    with pytest.raises(IncompatibleOperatorComposition, match="angular_trace"):
        _ = trace.outflow_restriction("xmin") @ trace.outflow_restriction("xmax")


@pytest.mark.parametrize("quad_name", _ALL)
def test_composing_the_wrong_HALF_is_REFUSED(quad_name):
    r"""``γ₊ @ γ₋`` — ``Γ(f) != Γ₋(f)``, so the chain does not connect.

    The sibling of the cross-face gate, one axis over: it is the ROLE rather
    than the face that fails to match.
    """
    trace = _trace(quad_name)
    with pytest.raises(IncompatibleOperatorComposition, match="angular_trace"):
        _ = trace.outflow_restriction("xmin") @ trace.inflow_restriction("xmin")


@pytest.mark.parametrize("quad_name", _ALL)
def test_a_space_contradicting_its_own_length_is_REFUSED_at_construction(quad_name):
    """The lengths and the spaces are redundant; a disagreement must not pass.

    ``n_total`` / ``n_restricted`` and the bound spaces describe the same fact
    until G6.5 retires the former. A mis-binding is invisible at apply-time
    (the arrays still broadcast), so it is refused at construction.
    """
    trace = _trace(quad_name)
    gamma = trace.outflow_restriction("xmin")

    with pytest.raises(ValueError, match="along axis 0"):
        TraceRestrictionOperator(
            gamma.indices, n_total=gamma.n_total, axis=0,
            codomain=gamma.domain,          # the FULL tier, not the half
        )


@pytest.mark.parametrize("quad_name", _ALL)
def test_a_restriction_still_constructs_UNBOUND(quad_name):
    """Binding is optional at this step — the tree-wide mandate is #330.

    Pins the transitional state honestly: the realizer's own producer
    (``sn/boundary/realizer.py``) is still unbound at step 1, so a gate
    asserting binding is mandatory would be false today.
    """
    trace = _trace(quad_name)
    bare = TraceRestrictionOperator(
        trace.outflow_indices_for_face("xmin"),
        n_total=int(trace.omega_dot_n.shape[1]),
    )
    assert bare.domain is None
    assert bare.codomain is None


@pytest.mark.parametrize("quad_name", _ALL)
def test_a_tier_carries_its_face_and_role_as_data(quad_name):
    """``face`` / ``role`` are readable, so a diagnostic need not parse the name."""
    trace = _trace(quad_name)
    space = trace.outflow_space("xmax")

    assert space.face == "xmax"
    assert space.role == "outflow"
    assert isinstance(space, AngularFaceTraceSpace)


# ─────────────────────────────────────────────────────────────────────
# The ordinate DIRECTIONS each tier carries (campaign P6 step 1)
#
# A boundary source is a function of direction, q(Ω). A tier handed only a
# SHAPE can express nothing angular — which is why the prescribed-inflow spec
# had to smuggle μ through its constructor. These rows pin the field that
# removes the smuggling, and the property that makes it usable: row order.
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("quad_name", _ALL)
@pytest.mark.parametrize("role", ["full", "outflow", "inflow"])
def test_a_tier_carries_its_ordinate_directions_in_ROW_order(quad_name, role):
    r"""``directions[i]`` is the direction of ordinate ``ordinate_indices[i]``.

    ⭐ **Row order is the whole contract.** A consumer reads row ``i`` of a
    vector in this space and needs to know which direction it belongs to; if
    ``directions`` were merely the right SET, every angularly-varying source
    would be silently permuted and no shape check could see it. Asserted
    against the quadrature's own nodes indexed by the tier's own index set —
    two independently-stored arrays that must agree element-wise.
    """
    quad = _QUADRATURES[quad_name]()
    trace = _trace(quad_name)
    space = {
        "full": trace.face_space, "outflow": trace.outflow_space,
        "inflow": trace.inflow_space,
    }[role]("xmin")

    np.testing.assert_array_equal(
        np.asarray(space.directions),
        np.asarray(quad.nodes)[np.asarray(space.ordinate_indices)],
        err_msg=(
            f"[{quad_name}/{role}] the tier's directions are not its own rows' "
            f"directions — a source reading direction i would be applied to a "
            f"different ordinate's row"
        ),
    )


@pytest.mark.parametrize("quad_name", _ALL)
def test_the_directions_are_dimension_generic(quad_name):
    r"""``(m,)`` in 1-D and ``(m, 3)`` in 3-D, by the SAME leading-axis index.

    The restriction indexes the leading axis only, so it carries the ordinate
    axis and leaves the component axis alone. Without this row a 1-D-only
    fixture set would let a ``directions[indices]`` regress to something that
    happens to work on a flat array — `[M]` 1-D rules give ``(m,)`` and
    ``product(4,4)`` / ``lebedev(17)`` give ``(m, 3)``.
    """
    quad = _QUADRATURES[quad_name]()
    space = _trace(quad_name).inflow_space("xmin")
    m = int(np.size(space.ordinate_indices))
    expected = (m,) if int(quad.dim) == 1 else (m, 3)

    assert np.asarray(space.directions).shape == expected, (
        f"[{quad_name}] dim={quad.dim} tier directions have shape "
        f"{np.asarray(space.directions).shape}, expected {expected}"
    )


@pytest.mark.parametrize("quad_name", _ALL)
def test_the_two_halves_directions_partition_the_full_tiers(quad_name):
    r"""Γ₊ ⊔ Γ₋ directions are disjoint subsets of Γ(f)'s — no row is shared.

    The directions inherit the tiers' own disjointness (already gated on the
    INDEX sets), so this row would be redundant if the two were built from one
    source. They are: ``directions`` and ``ordinate_indices`` are restricted by
    the same expression. The row exists so that a future edit which sourced the
    directions differently — from a re-derived mask, say — reddens instead of
    silently admitting a direction into both halves.
    """
    trace = _trace(quad_name)
    plus = np.atleast_2d(np.asarray(trace.outflow_space("xmin").directions).T).T
    minus = np.atleast_2d(np.asarray(trace.inflow_space("xmin").directions).T).T
    full = np.atleast_2d(np.asarray(trace.face_space("xmin").directions).T).T

    rows = lambda a: {tuple(np.round(r, 15)) for r in a}  # noqa: E731
    assert rows(plus).isdisjoint(rows(minus)), (
        f"[{quad_name}] a direction appears in BOTH half-traces"
    )
    assert rows(plus) | rows(minus) <= rows(full), (
        f"[{quad_name}] a half-trace carries a direction the full tier does not"
    )


@pytest.mark.parametrize("quad_name", _ALL)
def test_the_inflow_directions_all_point_INWARD(quad_name):
    r"""⭐ The physics the index set encodes, asserted on the directions.

    :math:`\Gamma_-(f) = \{\Omega : \Omega\cdot\hat n_f < 0\}`. Every other row
    here checks bookkeeping (the directions match the indices); this one checks
    that the bookkeeping means what its name says — computed from the
    directions and the face normal, NOT read back from ``omega_dot_n``, which
    is the array the index set was built from and would make the check
    circular.

    ``xmin``'s outward normal is :math:`-\hat x`, so inflow is
    :math:`\mu_x > 0`. `[M]` on ``gauss_legendre(8)``: rows ``[4 5 6 7]``,
    ``μ ∈ [+0.183, +0.960]``.
    """
    space = _trace(quad_name).inflow_space("xmin")
    directions = np.asarray(space.directions)
    mu_x = directions if directions.ndim == 1 else directions[:, 0]

    assert np.all(mu_x > 0.0), (
        f"[{quad_name}] Γ₋(xmin) carries a direction with μ_x = "
        f"{float(mu_x.min()):.6e} ≤ 0 — that ordinate leaves the domain "
        f"through xmin, so it belongs to Γ₊"
    )
    # PC⁻: the opposite face must be the mirror claim, or the row above would
    # pass for a space that put every ordinate in every inflow set.
    mu_x_max = np.asarray(_trace(quad_name).inflow_space("xmax").directions)
    mu_x_max = mu_x_max if mu_x_max.ndim == 1 else mu_x_max[:, 0]
    assert np.all(mu_x_max < 0.0), (
        f"[{quad_name}] Γ₋(xmax) carries a direction with μ_x ≥ 0"
    )
