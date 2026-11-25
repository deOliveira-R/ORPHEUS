r"""The Riesz legs — the Hilbert adjoint's factors as first-class arrows.

CS4c step 1 (R2 ruling, 2026-08-30) re-expressed ``AdjointOperator`` as the
named composition

.. math::

   A^{*} \;=\; \sharp_V \circ A^{\mathsf T} \circ \flat_W
   \qquad (\text{``RieszRaiseOperator(domain)} \circ \text{A.dual()} \circ
   \text{RieszLowerOperator(codomain)''})

with the legs delegating to the space's single-sourced metric verbs. This
module is the legs' DIRECT witness set (the composition's end-to-end
witnesses predate it: ``test_factored_adjoint_identity``,
``test_inverse_adjoint_coherence``, the monomorphic-leaves G1.4 grid):

* **G-H1** — the round trip ``♯ ∘ ♭`` is :math:`P_{\mathrm{range}(G)}`,
  NOT the identity: identity exactly on a strictly-positive metric (the
  positivity precondition is asserted so the row cannot silently stop
  meaning what it says), the null-space-zeroing projector on a singular
  one. `[M]` 2026-08-30 (pre-carve round §1.5): every monomorphic-leaves
  geometry carries 0 tangential trace slots — blind to this distinction —
  while a legal 2-D ``Quadrature.product(4,4)`` mesh has 32/64 tangential
  slots and round-trip trace defect ``2.871``. The zero-weight fixture
  here realizes the same law minimally. **First RED**: implement ``♯`` as
  ``1/G`` instead of Moore–Penrose — the singular row reads inf/nan.
* **G-H2** — the defining Riesz property
  :math:`\langle \flat x, y\rangle_{\mathrm{Euclid}} = \langle x,
  y\rangle_G`, RHS hand-spelled from the raw weight array. **First RED**:
  a leg that drops the weights reads O(1).
* **G-H3** — the composite USES the legs (vv Mode 11): counters
  monkeypatched over the leg classes' ``apply`` must move under one
  ``A.H.apply``. Without it, a composite silently calling ``apply_metric``
  directly is indistinguishable from one routing through the legs.
* **G-H4** — arrow bookkeeping and the dagger laws as object identities:
  ``A.H.H is A``, ``A.dual().dual() is A``, the dual arrow's ends are the
  swapped duals, the dagger–dual square commutes bitwise, and the legs
  ADMIT a ``DualSpace``, which carries the inverse form G⁻¹ (S6a, ruling
  26): ``lower_{V*}(lower_V(x)) = x``. Until S6a the legs refused a dual
  space, because the dual then carried G and the double lowering read
  ``G²x = [0.25, 4, 16]`` for ``w = [.5, 2, 4]`` (the G² trap).
* **G-A1 (theorem legs)** — ``apply`` bit-identical to the pre-leg inline
  formula, and ``apply_transpose`` bit-identical to the hand-built
  :math:`(A^{*})^{\mathsf T} = G_W\,A\,G_V^{+}` (#375's closure). Any
  rtol here is a red flag: the reduction tree is preserved by
  construction, so only ``np.array_equal`` is honest.

Ordering note: the negative leg of G-A1 (legs composed in the WRONG order)
is measured at ``4.226e-01 = |1 - c|`` on a flat metric and O(1) generally
— asserted as a genuine inequality so the identity rows cannot be
vacuously satisfied by a metric that cancels.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.operator import (
    LinearOperator,
    MissingAdjoint,
    RieszLowerOperator,
    RieszRaiseOperator,
    AdjointOperator,
)
from orpheus.numerics.face_layout import FaceLayout
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.space import DualSpace, FunctionSpace
from orpheus.numerics.spaces.angular_trace_space import AngularTraceSpace
from orpheus.numerics.factor import Kind, Variance, VarianceChart
from orpheus.numerics.space import TensorSpace
from orpheus.numerics.metric import DiagonalMetric
from tests.gates.numerics._declared_support import declared_mod_kernel


# ── fixtures ──────────────────────────────────────────────────────────────

_W_POSITIVE = np.array([0.5, 2.0, 4.0])          # strictly positive metric
_W_SINGULAR = np.array([2.0, 0.0, 3.0, 0.0])     # two null slots (the
# tangential-trace shape: |Ω·n| = 0 slots carry zero metric weight)


def _v() -> FunctionSpace:
    return TensorSpace.declared("V_pos", (3,), Variance(Kind.FUNCTION, VarianceChart.DENSITY), form=DiagonalMetric(_W_POSITIVE))


def _vs() -> FunctionSpace:
    return declared_mod_kernel("V_sing", (4,), Variance(Kind.FUNCTION, VarianceChart.DENSITY), form=DiagonalMetric(_W_SINGULAR))


def _w() -> FunctionSpace:
    return TensorSpace.declared("W_pos", (2,), Variance(Kind.FUNCTION, VarianceChart.DENSITY), form=DiagonalMetric(np.array([3.0, 0.25])))


class _Rect(LinearOperator):
    """A bound non-endomorphism ``A : V → W`` with a working transpose —
    the minimal fixture that exercises BOTH ends' metrics distinctly."""

    _M = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])

    def __init__(self) -> None:
        self._domain = _v()
        self._codomain = _w()

    @property
    def domain(self) -> FunctionSpace:
        return self._domain

    @property
    def codomain(self) -> FunctionSpace:
        return self._codomain

    def apply(self, x):
        return self._M @ x

    def apply_transpose(self, y):
        return self._M.T @ y

    @property
    def is_adjointable(self) -> bool:
        return True


class _NoTranspose(LinearOperator):
    """Adjointability negative control: no working transpose."""

    @property
    def domain(self):
        return None

    @property
    def codomain(self):
        return None

    def apply(self, x):
        return x


# ── G-H1: the round trip is P_range(G) ────────────────────────────────────


def test_round_trip_is_identity_on_a_strictly_positive_metric():
    """``♯(♭ x) == x`` exactly — WITH the positivity precondition asserted."""
    V = _v()
    if not np.all(_W_POSITIVE > 0):
        pytest.fail(
            "fixture precondition broken: this row's claim (round trip = "
            "identity) holds ONLY for a strictly positive metric"
        )
    x = np.array([0.2, 1.1, -0.4])
    rt = RieszRaiseOperator(V).apply(RieszLowerOperator(V).apply(x))
    np.testing.assert_array_equal(rt, x)


def test_round_trip_is_the_range_projector_on_a_singular_metric():
    """``♯ ∘ ♭ = P_range(G)``: identity on the support, ZERO on the null
    slots — the Moore–Penrose law. A ``1/G`` implementation reads inf here
    while every strictly-positive fixture stays green: that contrast is
    this row's value (`[M]` the monomorphic-leaves corpus is blind — all
    four geometries carry 0 tangential slots)."""
    Vs = _vs()
    x = np.array([1.0, 2.0, 3.0, 4.0])
    rt = RieszRaiseOperator(Vs).apply(RieszLowerOperator(Vs).apply(x))
    expected = np.where(_W_SINGULAR != 0.0, x, 0.0)
    assert np.all(np.isfinite(rt)), (
        f"round trip produced non-finite values {rt} — ♯ is 1/G, not the "
        f"Moore–Penrose pseudo-inverse the metric doctrine mandates"
    )
    np.testing.assert_array_equal(rt, expected)


# ── G-H2: the defining Riesz property ─────────────────────────────────────


@pytest.mark.l0
@pytest.mark.verifies("spaces-riesz-lower-raise")
def test_lowering_realizes_the_metric_pairing():
    """``⟨♭x, y⟩_Euclid == ⟨x, y⟩_G`` with the RHS hand-spelled from the
    raw weight array (the same reduction tree — bit-identity is honest)."""
    V = _v()
    x = np.array([0.7, -1.2, 0.9])
    y = np.array([2.0, 0.3, -0.8])
    lhs = float(np.sum(RieszLowerOperator(V).apply(x) * y))
    rhs = float(np.sum((_W_POSITIVE * x) * y))
    assert lhs == rhs, f"Riesz pairing broken: {lhs} != {rhs}"


# ── G-H3: the composite USES the legs (Mode 11) ───────────────────────────


def test_the_adjoint_composite_routes_through_the_legs(monkeypatch):
    """Counters over the leg classes' ``apply`` must MOVE under one
    ``A.H.apply`` — else a composite silently bypassing the legs (calling
    ``apply_metric`` directly) is indistinguishable from the real one, and
    the per-leg mutation battery upstream is mutating dead seams."""
    calls = {"lower": 0, "raise": 0}
    lower_apply = RieszLowerOperator.apply
    raise_apply = RieszRaiseOperator.apply

    def counting_lower(self, x):
        calls["lower"] += 1
        return lower_apply(self, x)

    def counting_raise(self, x):
        calls["raise"] += 1
        return raise_apply(self, x)

    monkeypatch.setattr(RieszLowerOperator, "apply", counting_lower)
    monkeypatch.setattr(RieszRaiseOperator, "apply", counting_raise)
    A = _Rect()
    A.H.apply(np.array([1.0, -1.0]))
    assert calls["lower"] == 1 and calls["raise"] == 1, (
        f"A.H.apply did not route through the legs (counters {calls}) — "
        f"the composition is bypassing its own factors"
    )


# ── G-H4: arrow bookkeeping + the dagger laws as object identities ────────


def test_the_dual_arrow_carries_the_swapped_duals():
    """For ``A : V → W``: ``A.dual() : W* → V*`` — both ends are
    ``DualSpace``s of the SWAPPED originals (the exchange realized in
    code, per the R2 ruling's articulation)."""
    A = _Rect()
    D = A.dual()
    assert D.domain == A.codomain.dual() and D.domain.dual() == A.codomain
    assert D.codomain == A.domain.dual() and D.codomain.dual() == A.domain


def test_dagger_and_dual_involutions_are_object_identities():
    """``A.H.H is A`` and ``A.dual().dual() is A`` — no arithmetic, no
    double wrapper. `[M]` pre-CS4c both were unreachable (#375:
    ``A.H.is_adjointable`` read False and the transpose raised)."""
    A = _Rect()
    assert A.H.H is A
    assert A.dual().dual() is A
    assert A.H.is_adjointable is True


def test_the_dagger_dual_square_commutes_bitwise():
    """``(A.dual()).H == (A.H).dual()`` — a LAW since S6a, with two routes:
    the left side is the base sandwich ``♯_{W*} ∘ A ∘ ♭_{V*}`` over the dual
    metrics, the right side the transpose of the primal adjoint
    ``♭_W ∘ A ∘ ♯_V``. Bit-identity holds because a dual metric is its
    primal with the faces swapped, so the two routes run the same three
    delegations in the same order. (Until S6a ``_DualOperator.adjoint``
    returned the right side, and this row compared one route with itself.)"""
    A = _Rect()
    probe = np.array([0.9, -1.2, 0.3])
    np.testing.assert_array_equal(
        A.dual().H.apply(probe), A.H.dual().apply(probe)
    )


def test_the_legs_admit_a_dual_space():
    """The G² trap is gone (S6a, ruling 26): ``V*`` carries ``G⁻¹``, so
    ``♭_{V*}`` applies ``G⁻¹`` and the double lowering is the identity,
    exactly, on the power-of-two weights ``w = [.5, 2, 4]`` (`[M]` with
    ``V*`` carrying ``G`` it read ``G²x = [0.25, 4, 16]``). Until S6a both
    legs refused a ``DualSpace`` (``TypeError: PRIMAL space only``)."""
    V = _v()
    dual = V.dual()
    low_d, high_d = RieszLowerOperator(dual), RieszRaiseOperator(dual)
    assert low_d.domain == dual and low_d.codomain == V
    assert high_d.domain == V and high_d.codomain == dual
    x = np.ones(3)
    np.testing.assert_array_equal(low_d.apply(RieszLowerOperator(V).apply(x)), x)
    # positive leg (vv #11): the primal constructs both, ends correct.
    low, high = RieszLowerOperator(V), RieszRaiseOperator(V)
    assert low.domain == V and low.codomain == V.dual()
    assert high.codomain == V and high.domain == V.dual()


def test_dual_refuses_a_transposeless_operator():
    """``.dual()`` is gated exactly like ``.H``: no working transpose, no
    dual arrow — eagerly, in the broken-stub-refusing house style."""
    with pytest.raises(MissingAdjoint, match="dual"):
        _NoTranspose().dual()


# ── G-A1: bit-identity of the re-expression + the transpose theorem ───────


@pytest.mark.l0
@pytest.mark.verifies("spaces-adjoint-riesz-composition")
def test_adjoint_apply_is_bit_identical_to_the_inline_formula():
    """The leg composition reproduces ``G_V⁺ ⊙ Aᵀ(G_W ⊙ y)`` to the BIT —
    same call order, same delegation targets. ``array_equal``, not
    allclose: an rtol here would mask a moved reduction tree."""
    A = _Rect()
    y = np.array([1.3, -0.7])
    inline = A.domain.apply_inverse_metric(
        A.apply_transpose(A.codomain.apply_metric(y))
    )
    np.testing.assert_array_equal(A.H.apply(y), inline)


@pytest.mark.l0
@pytest.mark.verifies("spaces-adjoint-riesz-composition")
def test_adjoint_transpose_is_the_leg_theorem():
    """``(A*)ᵀ = G_W A G_V⁺`` (#375's four-line composition, landed) —
    hand-built RHS, bit-identical; and the WRONG-ORDER composition
    genuinely differs (the identity rows are not vacuous)."""
    A = _Rect()
    x = np.array([0.2, 1.1, -0.4])
    hand = A.codomain.apply_metric(A.apply(A.domain.apply_inverse_metric(x)))
    H = A.H
    assert isinstance(H, AdjointOperator)  # type-narrowing only (-O strips)
    np.testing.assert_array_equal(H.apply_transpose(x), hand)
    # negative leg: swapping ♭/♯ roles moves the value O(1) on these
    # non-flat metrics — the composition order is load-bearing.
    wrong = A.codomain.apply_inverse_metric(A.apply(A.domain.apply_metric(x)))
    assert not np.allclose(hand, wrong), (
        "the negative control collapsed: the fixture's metrics cancel and "
        "these identity rows are vacuous on it — change the weights"
    )

# ── Trace-space rows: the raising leg on the SN boundary-trace metric ─────
#
# The inverse metric G^+ : V* -> V of a space is the raising leg (the arrow's two ends are G-H4's, above). Here its laws on the
# angular trace space, whose metric G = |Omega . n| w_n is EXACTLY zero on tangential ordinates, so the pseudo-inverse masking
# is the normal case. Layer DAG Phase 2: the trace space is built in numerics (``AngularTraceSpace.from_quadrature_and_layout``),
# not through an ``SNProblem`` (`[M]` the same space, see ``_trace_space``).

#: The four faces of a 3 × 3 two-group 2-D box, each slot ``(N, ng, 3)`` — the
#: layout the SN hub derives for that box (``SNProblem.boundary_face_layout``).
_FACES = ("xmin", "xmax", "ymin", "ymax")


def _trace_space(quadrature: Quadrature) -> AngularTraceSpace:
    """The boundary-trace space of a 3 × 3, two-group 2-D box, built from the
    quadrature and the face layout alone (numerics; no mesh, no hub).

    ``[M]`` 2026-10-30 (``scratch/layer_dag/refiling/specs_2/proto/probe_trace.py``):
    ``==`` the space the SN hub mints for the all-reflective 3 × 3 ``Mesh2D``,
    ``inner_product_weights`` ``array_equal``, for both quadratures below.
    """
    n = quadrature.n_ordinates
    layout = FaceLayout.from_named_shapes([(face, (n, 2, 3)) for face in _FACES])
    return AngularTraceSpace.from_quadrature_and_layout(quadrature, layout)


#: ``level_symmetric`` has NO tangential ordinate, ``product`` has many —
#: the pair is the point (a nondegenerate and a degenerate metric).
_QUADRATURES = {
    "level_symmetric(4)": Quadrature.level_symmetric(4),
    "product(4,4)": Quadrature.product(n_mu=4, n_phi=4),
}


@pytest.mark.foundation
@pytest.mark.parametrize("name", sorted(_QUADRATURES))
def test_the_raising_leg_is_self_adjoint_on_a_trace_space(name):
    """A real diagonal weight is its own transpose: ``apply_transpose(f)`` is bit-equal to ``apply(f)`` and
    ``<apply(f), g> = <f, apply(g)>`` over ``V*`` under the canonical (Euclidean) pairing."""
    space = _trace_space(_QUADRATURES[name])
    operator = RieszRaiseOperator(space)
    rng = np.random.default_rng(0)
    shape = np.asarray(space.inner_product_weights).shape
    f, g = rng.standard_normal(shape), rng.standard_normal(shape)
    assert operator.is_adjointable
    np.testing.assert_array_equal(operator.apply(f), operator.apply_transpose(f))
    np.testing.assert_allclose(np.sum(operator.apply(f) * g), np.sum(f * operator.apply(g)), rtol=1e-12)


@pytest.mark.foundation
def test_the_raising_leg_has_no_solve_and_is_spelled_by_ABSENCE():
    """``is_invertible`` is ``False`` AND there is no ``inverse`` method.

    On a degenerate metric this is a Moore–Penrose pseudo-inverse, so
    ``G⁺G ≠ I`` and "invert me" is not a question it can answer.  The
    house spelling for that is the ABSENCE of the method
    (``TraceRestrictionOperator``), never a raising stub — so misuse is
    a static error rather than a runtime one.
    """
    operator = RieszRaiseOperator(_trace_space(_QUADRATURES["product(4,4)"]))
    assert operator.is_invertible is False
    assert not hasattr(operator, "inverse")


@pytest.mark.foundation
def test_the_raising_leg_masks_the_trace_metrics_null_space():
    """`[M]` 192 of 384 ``product(4,4)`` trace rows carry ``G == 0`` exactly.

    On those rows the pseudo-inverse returns 0 — it does NOT divide.
    The ``level_symmetric`` control has zero such rows, which is exactly
    why it cannot witness this and the pair must both be exercised.
    """
    degenerate = _trace_space(_QUADRATURES["product(4,4)"])
    weights = np.asarray(degenerate.inner_product_weights)
    null_rows = np.abs(weights.ravel()) == 0.0
    assert null_rows.any(), "fixture no longer has a degenerate metric"

    ones = np.ones(weights.shape)
    out = RieszRaiseOperator(degenerate).apply(ones).ravel()
    assert np.all(out[null_rows] == 0.0)
    assert np.all(np.isfinite(out)), "a masked row must not divide by zero"

    control = np.asarray(
        _trace_space(_QUADRATURES["level_symmetric(4)"]).inner_product_weights
    )
    assert not np.any(np.abs(control.ravel()) == 0.0), (
        "level_symmetric is the NONDEGENERATE control; if it grows a "
        "tangential ordinate the pair stops discriminating."
    )


@pytest.mark.foundation
@pytest.mark.parametrize("name", sorted(_QUADRATURES))
def test_the_raising_leg_inverts_the_trace_metric_off_the_null_space(name):
    """``G (G⁺ x) == x`` wherever the metric is nonzero (``x`` annihilates ``ker G``) — the actual contract."""
    space = _trace_space(_QUADRATURES[name])
    weights = np.asarray(space.inner_product_weights)
    live = np.abs(weights.ravel()) > 0.0
    x = np.ones(weights.shape)
    round_trip = space.apply_metric(RieszRaiseOperator(space).apply(x))
    np.testing.assert_allclose(round_trip.ravel()[live], x.ravel()[live])
