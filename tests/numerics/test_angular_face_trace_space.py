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
