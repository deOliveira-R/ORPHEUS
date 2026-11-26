r"""Pre-carve anchors for CS4c step 6 item 6.2c — the angular moment HEAD becomes axis-built.

This docstring is written in the tense of 2026-09-07, before the carve: its present tense describes ``main`` @ ``79d2944a``. The carve landed, and the
tensor space of factors (N04) has since replaced the axes-less arm, the ``axes`` field and ``inner_product_weights`` on a head; the rows below are read as the
record of what the heads were at that date.

**What this file is.** Item 6.2c makes the angular head of every moment space
an AXIS-BUILT space: one rank-``r`` :class:`~orpheus.numerics.axis.Axis`
(``(L+1, 2L+1)`` for the rectangular harmonics, ``(L+1,)`` for the flat
Legendre family), ``Modal(unit)``, carrying the head's own measure, with
``inner_product_weights = None`` (the "never both" construction guard).  It is
the last axes-less factor on every production tensor product.

Every row below is measured on the **unmodified tree** (`main` @ ``79d2944a``)
and is GREEN today.  Three of them are RECORDS of a state the carve re-poses —
they are marked ``⛔ RE-POSED BY 6.2c`` in their own docstrings, so a session
reading a red here knows whether it is a regression or the carve arriving.

**Why the head can be simulated without a production edit.**  The head classes
already inherit ``axes`` from :class:`~orpheus.numerics.space.FunctionSpace`,
and ``__post_init__`` accepts a single axis whose shape IS the head's shape.
So ``replace(head, inner_product_weights=None, axes=(Axis(...),))`` builds
exactly the object 6.2c will mint — with the SAME class, name, shape and ``L``
— and the whole equivalence is measurable before any code moves.

**The measured headline `[M]` 2026-09-07 (`scratch/_step6_2c/p12_axis_head_feasible.py`).**

===========================================================  ==========================
claim                                                        measured
===========================================================  ==========================
head ``apply_metric``  axis route vs dense-slot route        **bit-identical, 0.0**
head ``inner_product``                                       **exactly equal**
product ``apply_metric`` axis-threaded vs ``FactoredMetric``  **bit-identical, 0 ULP**
two axis-built mints of one head                             ``==`` and hash-equal
axis-built head ``==`` its name-built twin                   **False** (the identity flip)
``head * bulk`` with an axis-built bulk                       takes the AXIS arm: 3 axes,
                                                             ``metric is None``
``factors`` / ``factors[0]`` / ``MomentHead`` narrowing        survive
``has_coordinate_cone``                                       ``None`` → **False**
===========================================================  ==========================

⟹ 6.2c's metric unification is a **bit-identical re-expression** on every
diagonal-Gram row, not a ``nulp``-banded one (contrast item 6.2a, whose dense →
factored move measured 2 ULP).  That is the acceptance tier this file pins.

**Activation evidence.**  Each gate names, in its own docstring, the mutation
that reddens it.  The two structural RECORDS (``truncated`` and the
one-metric-source guard) are red-BEFORE facts for 6.2c and are asserted here in
their TODAY form so the carve's own red is attributable.

Foundation mark: software invariants (identity, metric arithmetic, layout
surface); no physics claim rides here.
"""

from __future__ import annotations

from orpheus.numerics.basis.spherical_harmonic_basis import SphericalHarmonicBasis
from orpheus.numerics.basis.legendre_basis import LegendreBasis

from collections.abc import Callable
from dataclasses import replace

import numpy as np
import pytest

from orpheus.numerics.axis import Declared, Axis, BasisKind
from orpheus.numerics.metric import DenseMetric, Face, FactoredMetric, NoForm
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.space import FunctionSpace, TensorSpace
from orpheus.numerics.coefficient_layout import HarmonicLayout, LegendreLayout
from orpheus.numerics.spaces.moment_head import head_layout, moment_head_axis
from orpheus.numerics import frame_faces
from orpheus.numerics.axis import Modal
from orpheus.numerics.factor import Factor, Kind, VarianceChart
from orpheus.numerics.factor import Variance
from orpheus.numerics.metric import DiagonalMetric
from orpheus.numerics.declaration import Coord
from tests.gates.numerics._declared_support import declared_mod_kernel
from tests.gates.numerics._declared_support import whole_on_mod_kernel

pytestmark = pytest.mark.foundation

_HARMONIC_AXIS_LABEL = "harmonic"


# ── the head's two spellings: the basis-built production one (since 6.2c-ii), the hand weights that are the oracle of the arithmetic rows,
#    and the name-built partner of the identity rows ──


#: A head is a plain axis-built space since S7 (ruling 81); its layout is read
#: through :func:`head_layout`, not off a family class.
Head = FunctionSpace


def _literal_weights(label: str) -> np.ndarray:
    """The head's Parseval weights written from the definition, NOT read from the space under test: :math:`4\\pi/(2l+1)` on the live entries of
    degree :math:`l` (the harmonic layout is ``(L+1, 2L+1)`` with the first :math:`2l+1` entries of row :math:`l` live and the rest zero, the
    Legendre layout is ``(L+1,)``)."""
    L = int(label.split("_")[1][1])
    if label.startswith("sh_"):
        w = np.zeros((L + 1, 2 * L + 1))
        for l in range(L + 1):
            w[l, : 2 * l + 1] = 4.0 * np.pi / (2 * l + 1)
        return w
    return np.array([4.0 * np.pi / (2 * l + 1) for l in range(L + 1)])


def _hand_apply(w: np.ndarray, x: np.ndarray, *, inverse: bool) -> np.ndarray:
    """A numpy oracle of the head's metric (or its pseudo-inverse): ``w`` broadcast on the leading axes of ``x``, and the masked division (exactly
    zero where the weight is zero) written by hand. No form object and no space verb is called."""
    tail = (1,) * (x.ndim - w.ndim)
    wb = w.reshape(w.shape + tail)
    if inverse:
        return np.divide(x, wb, out=np.zeros_like(x), where=wb != 0.0)
    return wb * x


def _dense_slot_twin[T: Head](head: T) -> T:
    """The name-built partner of ``head``: a declared nodal space of the same name and shape whose diagonal form is the weights ``head`` applies to
    the unit (the head's own metric read back, so it is a PARTNER for the identity rows and is never an oracle for an arithmetic row: the arithmetic
    rows use :func:`_literal_weights` and :func:`_hand_apply`). The retired spelling stored the measure in a dense ``inner_product_weights`` slot; that
    slot is gone, and this partner differs from the head only in carrying no basis origin.
    """
    assert len(head.axes) == 1, "a production head is single-axis"
    form = head.axes[0].form
    # a diagonal head's metric is what the space applies to the unit; a dense or formless head has no weight array
    w = None if isinstance(form, (NoForm, DenseMetric)) else np.asarray(head.apply_metric(np.ones(head.shape)))
    if w is None:
        return TensorSpace.declared(head.name, head.shape, Variance(Kind.FUNCTION, VarianceChart.DENSITY))
    return declared_mod_kernel(head.name, head.shape, Variance(Kind.FUNCTION, VarianceChart.DENSITY), form=DiagonalMetric(w))


def _axis_built[T: Head](head: T) -> T:
    """The production head IS axis-built since item 6.2c-ii — identity, kept
    so the rows below read as they were written (the premise landed)."""
    assert len(head.axes) >= 1
    return head


def _bulk(ng: int = 2, spatial: tuple[int, ...] = (4, 3)) -> FunctionSpace:
    """An axis-built cell group, as every SN carrier's ``bulk_space`` is."""
    return TensorSpace.of_axes(Factor.whole_on(Axis("energy", (ng,), origin=Declared(BasisKind.NODAL)), Kind.MEASURE, VarianceChart.DENSITY), Factor.whole_on(Axis("spatial", spatial, origin=Declared(BasisKind.NODAL, DiagonalMetric(np.full(spatial, 0.5)))), Kind.MEASURE, VarianceChart.DENSITY))


_HEADS: dict[str, Callable[[], Head]] = {
    "sh_L0": lambda: TensorSpace.of_axes(whole_on_mod_kernel(Axis.of_basis(SphericalHarmonicBasis(L=0)), Kind.MEASURE, VarianceChart.DENSITY)),
    "sh_L1": lambda: TensorSpace.of_axes(whole_on_mod_kernel(Axis.of_basis(SphericalHarmonicBasis(L=1)), Kind.MEASURE, VarianceChart.DENSITY)),
    "sh_L2": lambda: TensorSpace.of_axes(whole_on_mod_kernel(Axis.of_basis(SphericalHarmonicBasis(L=2)), Kind.MEASURE, VarianceChart.DENSITY)),
    "sh_L3": lambda: TensorSpace.of_axes(whole_on_mod_kernel(Axis.of_basis(SphericalHarmonicBasis(L=3)), Kind.MEASURE, VarianceChart.DENSITY)),
    "legendre_L0": lambda: TensorSpace.of_axes(whole_on_mod_kernel(Axis.of_basis(LegendreBasis(L=0, axis="x")), Kind.MEASURE, VarianceChart.DENSITY)),
    "legendre_L1": lambda: TensorSpace.of_axes(whole_on_mod_kernel(Axis.of_basis(LegendreBasis(L=1, axis="x")), Kind.MEASURE, VarianceChart.DENSITY)),
    "legendre_L2": lambda: TensorSpace.of_axes(whole_on_mod_kernel(Axis.of_basis(LegendreBasis(L=2, axis="x")), Kind.MEASURE, VarianceChart.DENSITY)),
    "legendre_L3_z": lambda: TensorSpace.of_axes(whole_on_mod_kernel(Axis.of_basis(LegendreBasis(L=3, axis="z")), Kind.MEASURE, VarianceChart.DENSITY)),
}


# ═══════════════════════════════════════════════════════════════════════
# G1 — the metric arithmetic is a BIT-IDENTICAL re-expression
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("label", sorted(_HEADS))
def test_the_head_metric_is_bit_identical_on_the_axis_route(label: str) -> None:
    r"""The head's measure, carried by ONE axis, applies exactly as the hand numpy product of its stated weights.

    The head's metric is the walk of its one factor's form (``DiagonalMetric.apply_block``), which reshapes the axis's weights to the
    head's own leading block and multiplies: the SAME single elementwise product a numpy broadcast of the stated weights performs
    (the oracle here, written from the definition :math:`4\\pi/(2l+1)`, with the masked reciprocal on the inverse leg), so
    there is no re-association and ``array_equal`` is the honest assertion
    (``vv-principles`` §bit-identity criterion 3: zero reductions reordered).

    Reddens under: the battery's ``A5`` arm (the per-axis metric skips an
    axis) — `[M]` 2026-09-07 it reddens all 8 rows of this gate, 17 rows of
    this file and 43 of the 4501-row battery scope.
    """
    axis_head = _HEADS[label]()
    w = _literal_weights(label)
    assert axis_head.axes is not None
    assert w.shape == axis_head.shape, "the hand weights span the head's own layout"

    rng = np.random.default_rng(20260907)
    x = rng.standard_normal(axis_head.shape + (2, 4))
    np.testing.assert_array_equal(
        np.asarray(axis_head.apply_metric(x)), _hand_apply(w, x, inverse=False),
        err_msg="the axis route must be bit-identical to the hand product of the stated weights",
    )
    np.testing.assert_array_equal(
        np.asarray(axis_head.apply_inverse_metric(x)), _hand_apply(w, x, inverse=True),
        err_msg="the pseudo-inverse must agree bit-for-bit with the hand masked reciprocal (dead slots -> 0)",
    )
    # The (2, 4) tail is a payload: the pairing over it is posed on the explicit tensor product, and is pinned
    # against the hand metric contracted by hand (comparing two tensor products with each
    # other was blind to a shared metric defect: the U1 qa review, F3).
    payload = TensorSpace.of_axes(Factor.whole_on(Axis("payload", (2, 4), origin=Declared(BasisKind.NODAL)), Kind.MEASURE, VarianceChart.DENSITY))
    by_hand = float(np.sum(_hand_apply(w, x, inverse=False) * x))
    assert (axis_head * payload).inner_product(x, x) == by_hand


@pytest.mark.parametrize("label", sorted(_HEADS))
def test_the_product_metric_is_bit_identical_when_the_head_gains_axes(label: str) -> None:
    r"""``head * bulk`` applies the hand product of the stated weights, bit for bit.

    The product's metric is the Kronecker walk of the factors' forms in slot order (the one arm since N04). Before it, an axes-less
    head took the factored arm of the retired ``from_factors`` and an axis-built head the axis arm; the two applied the same per-block
    diagonals in the same order, so the values were bit-identical, which is what made 6.2c's metric unification a re-expression rather
    than a re-baseline. The oracle is the hand product of the literal weights (:func:`_literal_weights`, :func:`_hand_apply`).

    Reddens under: the battery's ``A5`` arm (a per-axis measure dropped) —
    `[M]` all 8 rows of this gate.
    """
    axis_head = _HEADS[label]()
    bulk = _bulk()
    axis_product = axis_head * bulk

    assert len(axis_product.axes) == 3
    assert axis_product.factors == axis_head.factors + bulk.factors, "the product is not the concatenation of the slots' factors"
    assert axis_product.shape == axis_head.shape + bulk.shape

    # The oracle is the hand product of the stated weights: the head's (4 pi/(2l+1)) on its own leading block, then the bulk's cell measure
    # (0.5 on every spatial node; the energy slot has no form), one broadcast product per formed slot in slot order.
    w = _literal_weights(label)
    rng = np.random.default_rng(5)
    x = rng.standard_normal(axis_product.shape)
    hand = _hand_apply(w, x, inverse=False)
    hand = hand * 0.5
    np.testing.assert_array_equal(
        np.asarray(axis_product.apply_metric(x)), hand,
        err_msg="the axis-threaded product must reproduce the hand product of the stated weights bit-for-bit",
    )
    assert axis_product.inner_product(x, x) == float(np.sum(hand * x))


def test_no_reachable_array_on_the_axis_product_is_state_sized() -> None:
    r"""The axis arm stores per-AXIS measures, never their outer product.

    The structural half of the memory claim (``L59c``): the separation leg
    below is a synthetic ratio, this one is exact and free — with three axes of
    shapes ``(3, 5)``, ``(2,)``, ``(4, 3)`` the product's index set has 360
    entries and NO array the space can reach has that size.
    """
    head = TensorSpace.of_axes(whole_on_mod_kernel(Axis.of_basis(SphericalHarmonicBasis(L=2)), Kind.MEASURE, VarianceChart.DENSITY))
    product = head * _bulk()
    assert product.axes is not None
    n = int(np.prod(product.shape))
    assert n == 360
    reachable = [
        ax.form.weights for ax in product.axes
        if isinstance(ax.form, DiagonalMetric) and not bool((ax.form.weights == 1.0).all())
    ]
    assert reachable, "the fixture must carry at least one non-counting axis"
    assert all(w.size < n for w in reachable)
    assert sum(w.nbytes for w in reachable) < 8 * n


def test_the_axis_route_never_densifies_a_large_product() -> None:
    r"""SEPARATION, sized on a synthetic pair so the leg is fixture-honest.

    ``L59c``: a small production fixture gives only an 8× ratio, which is too
    weak to be a keystone.  Two 2000-point axes: per-axis storage is 2 × 2000
    weights (32 000 B), the outer product would be 4 000 000 entries (32 MB).
    """
    w = np.linspace(0.5, 1.5, 2000)
    a = Axis("a", (2000,), origin=Declared(BasisKind.NODAL, DiagonalMetric(w)))
    b = Axis("b", (2000,), origin=Declared(BasisKind.NODAL, DiagonalMetric(w)))
    space = TensorSpace.of_axes(Factor.whole_on(a, Kind.MEASURE, VarianceChart.DENSITY), Factor.whole_on(b, Kind.MEASURE, VarianceChart.DENSITY))
    assert len(space.factors) == 2 and all(isinstance(f.form, DiagonalMetric) for f in space.factors), "the product does not carry one diagonal form per slot"
    assert isinstance(a.form, DiagonalMetric) and isinstance(b.form, DiagonalMetric)
    stored = a.form.weights.nbytes + b.form.weights.nbytes
    dense = 8 * int(np.prod(space.shape))
    assert dense // stored >= 100, f"separation {dense}/{stored} is too weak to be a keystone"


# ═══════════════════════════════════════════════════════════════════════
# G4 — the MomentHead surface survives the axis-building
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("label", sorted(_HEADS))
def test_the_moment_head_surface_survives(label: str) -> None:
    r"""An axis-built head answers the moment-head reads through the one door,
    with its basis's layout (S7, ruling 81: the layout is a value the basis
    carries and the head axis reads it).

    The production readers — ``_moment_head_of`` (``transfer.py``),
    ``FissionMomentOperator._head``, the field's ``head_layout``,
    ``_block_contraction`` (keyed on ``len(head.shape)``) and the per-degree
    contraction (``head.degree_block``) — read :func:`head_layout`. The retired
    name-built twin carries no axis, so the door refuses it.

    Reddens under: minting the head as a ``Declared`` axis (no layout).
    """
    axis_head = _HEADS[label]()
    layout = head_layout(axis_head)
    assert layout == moment_head_axis(axis_head).origin.basis.layout
    assert layout.shape == axis_head.shape
    with pytest.raises(TypeError, match="is not a moment head"):
        head_layout(_dense_slot_twin(axis_head))


@pytest.mark.parametrize("label", sorted(_HEADS))
def test_the_product_still_exposes_its_head_as_factor_zero(label: str) -> None:
    r"""``factors[0]`` survives the axis-threading arm.

    ``TensorSpace`` keeps the head as its first slot, so the tree query the
    moment carriers rely on (``head_layout(space.sub(<first label>))``) reads
    the head the product was built from. (``find_factor`` retired at S7 with
    the head classes it was the typed bridge to.)
    """
    axis_head = _axis_built(_HEADS[label]())
    product = axis_head * _bulk()
    assert isinstance(product, TensorSpace)
    head = product.sub(product.labels[0])
    assert head == axis_head
    assert head_layout(head) == head_layout(axis_head)


# ═══════════════════════════════════════════════════════════════════════
# The identity flip on the head — the FORK's mechanism
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("label", sorted(_HEADS))
def test_two_axis_built_mints_of_one_head_are_one_space(label: str) -> None:
    """Structural identity: same axes ⟹ same space, and hash-consistent."""
    a, b = _HEADS[label](), _HEADS[label]()
    assert a is not b
    assert a == b and hash(a) == hash(b)
    assert len({a, b}) == 1


@pytest.mark.parametrize("label", sorted(_HEADS))
def test_an_axis_built_head_is_never_equal_to_its_name_built_twin(label: str) -> None:
    r"""The identity flip (6.1): an axis-built space is not a name-built one.

    This is why the head's §6b set is the whole COMPARISON population and not
    just the construction sites — every partner has to move in the same commit
    (``plan-authoring`` §6b).
    """
    axis_head = _HEADS[label]()
    head = _dense_slot_twin(axis_head)
    assert axis_head.shape == head.shape
    assert axis_head != head
    assert head != axis_head


@pytest.mark.parametrize("L", [0, 1, 2])
def test_the_metric_enters_the_identity_once_the_head_is_axis_built(L: int) -> None:
    r"""⭐ **THE FORK, in one assertion — RULED and LANDED (R-6.2c-1, item
    6.2c-ii).**  Two heads that differ ONLY in their measure — the basis's
    continuum head and the frame's Parseval-dressed head — are UNEQUAL now
    that the head is axis-built (``Axis.__eq__`` compares weights bytes);
    their RETIRED dense-slot twins were ``==`` until 2026-10-30 (the
    metric-blind ``(name, shape)`` identity of a name-built space), which is
    the seam the tree carried until 6.2c-ii and the reason the fork had to be
    ruled before the head could be axis-built. Since U1 (space membership,
    2026-10-30) an axes-less space's identity carries its metric too, so the
    dense-slot twins separate as well and the seam is closed on both
    spellings. The tree binds the dressed one.
    """
    frame = frame_faces.angular_frame(Quadrature.lebedev(11), L)
    continuum = TensorSpace.of_axes(whole_on_mod_kernel(Axis.of_basis(frame.basis), Kind.MEASURE, VarianceChart.DENSITY))
    dressed = frame_faces.basis_space(frame)
    assert head_layout(continuum) == head_layout(dressed) == HarmonicLayout(L)
    assert continuum.axes is not None and dressed.axes is not None
    assert not np.array_equal(
        np.asarray(continuum.axes[0].form.weights), np.asarray(dressed.axes[0].form.weights),
    )
    assert continuum != dressed, "axis-built: the measure IS the identity — the two heads separate"
    assert _dense_slot_twin(continuum) != _dense_slot_twin(dressed), (
        "the name-built spelling carries its metric since U1: the metric-blind seam is closed on both spellings"
    )


# ═══════════════════════════════════════════════════════════════════════
# The two RED-BEFORE structural records 6.2c must repair
# ═══════════════════════════════════════════════════════════════════════


def test_a_dense_gram_lives_positioned_on_the_axis_built_heads_own_object() -> None:
    r"""RE-POSED at CS4c step 6 item 6.2c-i (ruling R-6.2c-2, 2026-09-08) —
    the construction guard an axis-built head used to collide with.

    Until 6.2c-i ``FunctionSpace.__post_init__`` refused ``axes`` + ``metric``
    together ("one metric source only"), and on a DENSE-Gram frame the
    Parseval-dressed head is exactly ``inner_product_weights=None`` + a
    :class:`DenseMetric` — no diagonal measure to put on an axis, so the
    axis-built head was UNCONSTRUCTIBLE for that row (`[M]` 17 of 75 shipped
    (rule, L) rows are DENSE, 2 of them at ``L ≤ 2``). Now the space's
    metric is a DERIVED object over its axes and the dense form is a block
    POSITIONED on it: admitted, and APPLIED.

    The control the old record carried survives with its teeth: dropping
    the object for a counting axis silently changes the value.
    """
    G = np.eye(3) + 0.25
    dense = DenseMetric(DenseMetric.symmetrised(G), kernel=Coord(3, ()), face=Face.GRAM)
    axis = Axis(_HARMONIC_AXIS_LABEL, (3,), origin=Declared(Modal(np.eye(3)[0]), dense))
    head = TensorSpace.of_axes(Factor.whole_on(axis, Kind.MEASURE, VarianceChart.DENSITY))
    x = np.array([1.0, 2.0, 3.0])
    # The block applies through DenseMetric.apply_block (a reshape-and-matmul on the block) while the standalone object applies
    # through DenseMetric.apply (the leading matmul): the SAME form, two associations — [M] 1 ULP apart on this probe, so the
    # honest tier is a nulp band, not bit-equality (vv-principles §bit-identity criterion 3).
    np.testing.assert_array_almost_equal_nulp(
        head.apply_metric(x), np.asarray(dense.apply(x)), nulp=4,
    )

    counting = TensorSpace.of_axes(Factor.whole_on(Axis(_HARMONIC_AXIS_LABEL, (3,), origin=Declared(Modal(np.eye(3)[0]))), Kind.MEASURE, VarianceChart.DENSITY))
    np.testing.assert_array_equal(counting.apply_metric(x), x)
    assert not np.allclose(np.asarray(dense.apply(x)), x), (
        "the control must MOVE, or dropping the form would be inert"
    )


def test_the_cone_answer_flips_from_unanswerable_to_refusing() -> None:
    r"""⛔ RE-POSED BY 6.2c — ``has_coordinate_cone`` goes ``None`` → ``False``.

    A MODAL factor makes the per-component sign test meaningless, so
    :meth:`~orpheus.numerics.field.Field.cone_violations` REFUSES rather than
    answering.  `[M]` 2026-09-07 the moment path has **zero production callers**
    of ``cone_violations`` (47 sites tree-wide, all definitions/prose in
    ``orpheus/`` and 35 test calls, none on a moment field), so the flip is
    predicted inert — battery arm ``C1`` measures it rather than asserting it.
    """
    axis_head = TensorSpace.of_axes(whole_on_mod_kernel(Axis.of_basis(SphericalHarmonicBasis(L=1)), Kind.MEASURE, VarianceChart.DENSITY))
    head = _dense_slot_twin(axis_head)            # the name-built twin is a declared nodal space: it answered None while it was a leaf
    assert head.has_coordinate_cone is True
    assert axis_head.has_coordinate_cone is False
    assert (axis_head * _bulk()).has_coordinate_cone is False
    assert (head * _bulk()).has_coordinate_cone is True


# ═══════════════════════════════════════════════════════════════════════
# The production heads, not only the hand-built ones
# ═══════════════════════════════════════════════════════════════════════


_RULES = {
    "gauss_legendre(8)": lambda: Quadrature.gauss_legendre(8),
    "level_symmetric(4)": lambda: Quadrature.level_symmetric(4),
    "lebedev(11)": lambda: Quadrature.lebedev(11),
    "product(4,6)": lambda: Quadrature.product(4, 6),
    "folded_product(4,8)": lambda: Quadrature.folded_product(4, 8),
}


@pytest.mark.parametrize("label", sorted(_RULES))
@pytest.mark.parametrize("L", [0, 1, 2])
def test_the_production_head_of_every_family_is_axis_buildable(label: str, L: int) -> None:
    r"""The head a shipped rule's frame actually binds — built through the
    production chain, not hand-minted (``vv`` #28: build the operand).

    Both heads are covered: the **continuum** head (``TensorSpace.of_basis(frame.basis, variance)``,
    the basis its generator) and the frame's **dressed** head
    (``basis_space(frame)``, the frame its generator) — each axis-built, each
    bit-identical to its own retired dense-slot twin on the metric verbs; on
    a DENSE-Gram row the dressed head carries no axis measure and its dense
    form is POSITIONED on the space's own metric object (hazard H-1, ruling
    R-6.2c-2), so the twin comparison is made on the continuum head only there.
    """
    from orpheus.numerics.frame import GramStructure

    # The numerics frame itself: transport's HarmonicFrame.from_galerkin is a behaviour-identical upgrade that SHARES
    # the frame (tests/gates/transport/frames/test_harmonic_frame.py), so the head it binds is this one ([M] 15 of 15
    # (rule, L) cases equal on both arms, layer DAG Phase 2 refresh, specs_2r/proto/probe_harmonic.py).
    frame = frame_faces.angular_frame(_RULES[label](), L)
    rng = np.random.default_rng(11)
    for head, (arm, generator) in ((TensorSpace.of_axes(whole_on_mod_kernel(Axis.of_basis(frame.basis), Kind.MEASURE, VarianceChart.DENSITY)), ("basis", frame.basis)), (frame_faces.basis_space(frame), ("frame", frame))):
        assert head_layout(head) in (HarmonicLayout(L), LegendreLayout(L)), (
            "every shipped rule binds one of the two head families"
        )
        assert len(head.factors) == 1 and isinstance(head.factors[0].form, (DenseMetric, DiagonalMetric)), "a head is one factor with a form"
        assert getattr(head.axes[0].origin, arm) is generator
        assert head_layout(head).L == L
        if isinstance(head.axes[0].form, DenseMetric):
            assert frame.discrete_gram_structure is GramStructure.DENSE and head is frame_faces.basis_space(frame)
            assert isinstance(head.factors[0].form, DenseMetric)
            continue
        twin = _dense_slot_twin(head)
        x = rng.standard_normal(head.shape + (2, 3))
        np.testing.assert_array_equal(
            np.asarray(head.apply_metric(x)), np.asarray(twin.apply_metric(x)),
        )
