r"""An axis-built space's metric is a DERIVED object over its axes — and a
dense Gram is a block POSITIONED on it (CS4c step 6 item 6.2c-i, ruling
R-6.2c-2, 2026-09-08).

**The doctrine, stated once.** The axes ARE the metric source of an
axis-built space (a measure is diagonal by nature); the space REALIZES that
source as ONE :class:`~orpheus.numerics.metric.FactoredMetric` — a
:class:`~orpheus.numerics.metric.DiagonalMetric` on each weighted axis's
block, ``None`` on a counting axis — and applies it through
:meth:`DiagonalMetric.apply_block`, the explicit leading-1s / block /
trailing-1s reshape-and-multiply.  A form no diagonal measure can spell
(a moment head whose discrete Gram is DENSE — `[M]` 17 of 75 shipped
(rule, L) rows) is a :class:`~orpheus.numerics.metric.DenseMetric`
POSITIONED on that axis's block of the space's own object — an OVERLAY
merged into the derived entries, never a replacement of them: the
construction guard admits ``axes`` + a metric object iff the object is a
``FactoredMetric`` with exactly one entry per axis, in order, a form only
on an axis that carries no measure, and never a diagonal one (a diagonal
measure is the axis's own to carry). Every other axis keeps supplying its
own block, and a product of axis-built factors carries its factors'
overlays concatenated beside its axes.

**The pair this file is (ruling O-2's gate pair):**

* the dense block IS APPLIED — negative control: the counting-axis route
  with no object, `[M]` the two disagree by :math:`O(1)` (the 6.2c
  verification round measured ``max|Δ| = 14.6`` on a dense-Gram head and
  68.7 % / 1784 % on the two ``L ≤ 2`` dense rows — a VALUE bug, which is
  why a merely relaxed guard was refused: the old per-axis short-circuit
  read the axes FIRST and ignored the object);
* with NO dense block the derived object is **bit-identical** to the
  retired inline per-axis loop — reimplemented HERE from its own
  reshape-and-multiply spelling as the oracle, never routed through the
  production helper (``vv`` #22: an oracle that shares the SUT's body is a
  tautology), on synthetic axes (a rank-2 head-shaped axis among them)
  and on the four ledger carriers' trial spaces.

Foundation mark: software invariants of the metric realization, no physics
claim.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.axis import Declared, Axis, BasisKind, EnergyAxis
from orpheus.numerics.metric import DenseMetric, Face, DiagonalMetric, FactoredMetric, NoForm
from orpheus.numerics.space import FunctionSpace
from orpheus.numerics.axis import Modal
from orpheus.numerics.quadrature import Quadrature
from tests.gates.numerics._phase_space_twins import (
    MONOMORPHIC_EDGES,
    MONOMORPHIC_RECTANGLE,
    spatial_axis_1d,
    spatial_axis_2d,
    twin,
)
from orpheus.numerics.factor import Factor, Kind, VarianceChart
from orpheus.numerics.space import TensorSpace
from orpheus.numerics.declaration import Coord
from tests.gates.numerics._declared_support import whole_on_mod_kernel

pytestmark = pytest.mark.foundation

#: The ledger carriers' trial spaces, as numerics twins of the four monomorphic meshes (layer DAG Phase 4; until then
#: read off the SN problems, which are all Diamond Difference, so their trial space is the angular bulk), plus the LD
#: slab, whose trial space carries the moment axis the docstring below names (no row had it before Phase 4).
_GEOMETRIES = {
    "slab": lambda: twin(Quadrature.gauss_legendre(n_ordinates=4), spatial_axis_1d("slab", MONOMORPHIC_EDGES)).angular_bulk,
    "sphere": lambda: twin(Quadrature.gauss_legendre(n_ordinates=4), spatial_axis_1d("sphere", MONOMORPHIC_EDGES)).angular_bulk,
    "cylinder": lambda: twin(Quadrature.folded_product(n_mu=4, n_phi=8), spatial_axis_1d("cylinder", MONOMORPHIC_EDGES)).angular_bulk,
    "cart2d": lambda: twin(Quadrature.level_symmetric(sn_order=4), spatial_axis_2d(*MONOMORPHIC_RECTANGLE)).angular_bulk,
    "slab_ld": lambda: twin(Quadrature.gauss_legendre(n_ordinates=4), spatial_axis_1d("slab", MONOMORPHIC_EDGES)).ld_trial(),
}


# ── the ORACLE: the retired inline loop, spelled here ────────────────────

def _inline_per_axis(axes, x: np.ndarray, *, inverse: bool) -> np.ndarray:
    """The reshape-and-multiply the retired ``_apply_axes_weights`` performed
    (leading 1s for the preceding axes' ranks, the axis shape, trailing 1s
    for the rest; the Moore–Penrose pseudo-inverse per axis) — an oracle
    independent of the production realization."""
    out = np.asarray(x)
    ndim = out.ndim
    start = 0
    for ax in axes:
        rank = len(ax.shape)
        w = ax.form.weights if isinstance(ax.form, DiagonalMetric) else None
        if w is not None:
            wb = np.asarray(w).reshape((1,) * start + tuple(ax.shape) + (1,) * (ndim - start - rank))
            if inverse:
                nonzero = wb != 0.0
                out = np.where(nonzero, out / np.where(nonzero, wb, 1.0), 0.0)
            else:
                out = out * wb
        start += rank
    return out


def _synthetic_axes() -> tuple[Axis, ...]:
    """A rank-2 head-shaped MODAL axis (zero-padded weights, like a
    spherical-harmonic head at L = 1), a counting energy axis, a weighted
    spatial axis — every shape the realization must position."""
    head_w = np.array([[4.0 * np.pi, 0.0, 0.0], [4.0 * np.pi / 3.0] * 3])
    return (
        Axis(label="harmonic_probe", shape=(2, 3), origin=Declared(Modal(np.eye(1, 6)[0].reshape(2, 3)), DiagonalMetric(head_w))),
        EnergyAxis(label="energy", shape=(2,)),
        Axis(label="spatial", shape=(4,), origin=Declared(BasisKind.NODAL, DiagonalMetric(np.array([0.5, 1.0, 1.5, 2.0])))),
    )


# ═════════════════════════════════════════════════════════════════════════
# Leg 1 — with NO dense block, the derived object is bit-identical to the
#         retired inline loop (synthetic axes, then the four carriers)
# ═════════════════════════════════════════════════════════════════════════

#: The counting payload axis a pairing over the trailing element axis is posed on.
_COUNTING_PAYLOAD = TensorSpace.of_axes(Factor.whole_on(Axis("payload", (5,), origin=Declared(BasisKind.NODAL)), Kind.MEASURE, VarianceChart.DENSITY))


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_the_derived_object_is_bit_identical_to_the_inline_loop_on_synthetic_axes(seed):
    """``apply_metric`` / ``apply_inverse_metric`` / ``inner_product`` on an
    axis-built space ≡ the inline per-axis oracle, ``array_equal`` — including
    a rank-2 axis with zero-padded weights (the pseudo-inverse's kernel) and
    an extra trailing element axis (the leading-aligned convention)."""
    axes = _synthetic_axes()
    space = TensorSpace.of_axes(*(whole_on_mod_kernel(_a, Kind.MEASURE, VarianceChart.DENSITY) for _a in axes))
    if not isinstance(space._resolved_metric, FactoredMetric):
        pytest.fail("the axis-built space did not DERIVE a FactoredMetric from its axes")
    rng = np.random.default_rng(seed)
    for shape in (space.shape, space.shape + (5,)):
        x = rng.standard_normal(shape)
        y = rng.standard_normal(shape)
        if not np.array_equal(space.apply_metric(x), _inline_per_axis(axes, x, inverse=False)):
            pytest.fail(f"apply_metric differs from the inline loop on shape {shape}")
        if not np.array_equal(space.apply_inverse_metric(x), _inline_per_axis(axes, x, inverse=True)):
            pytest.fail(f"apply_inverse_metric differs from the inline loop on shape {shape}")
        want = float(np.sum(_inline_per_axis(axes, x, inverse=False) * y))
        # A trailing payload has no inner product of its own: the pairing over it is posed on the
        # explicit tensor product with a counting payload axis (membership is the space's verb).
        pairing_space = space if shape == space.shape else space * _COUNTING_PAYLOAD
        if pairing_space.inner_product(x, y) != want:
            pytest.fail(f"inner_product differs from the inline loop's reduction on shape {shape}")


def test_an_all_counting_space_applies_the_identity():
    """Every axis counting or formless ⟹ every form is a no-op (a ``NoForm`` slot, or a counting form applied as a product with one) and the
    verbs are the identity."""
    space = TensorSpace.of_axes(Factor.whole_on(EnergyAxis(label="energy", shape=(3,)), Kind.MEASURE, VarianceChart.DENSITY), Factor.whole_on(Axis(label="spatial", shape=(2,), origin=Declared(BasisKind.NODAL)), Kind.MEASURE, VarianceChart.DENSITY))
    for f in space.factors:
        counting = isinstance(f.form, DiagonalMetric) and bool((f.form.weights == 1.0).all())
        if not (isinstance(f.form, NoForm) or counting):
            pytest.fail(f"the form of {f.axis.label!r} is neither NoForm nor counting: {f.form!r}")
    x = np.arange(6.0).reshape(3, 2)
    assert np.array_equal(space.apply_metric(x), x)
    assert space.inner_product(x, x) == float(np.sum(x * x))


@pytest.mark.parametrize("geometry", list(_GEOMETRIES), ids=list(_GEOMETRIES))
def test_the_carriers_trial_space_reads_bit_identically(geometry):
    """On the four ledger carriers' ``angular_trial_space`` (the interior of
    every composite: ``V_cell × w_n`` per axis, a moment axis on LD) the
    derived object reproduces the inline loop bit-for-bit, both directions
    and the pairing."""
    space = _GEOMETRIES[geometry]()
    axes = space.axes
    assert axes is not None
    rng = np.random.default_rng(7)
    x = rng.standard_normal(space.shape)
    y = rng.standard_normal(space.shape)
    if not np.array_equal(space.apply_metric(x), _inline_per_axis(axes, x, inverse=False)):
        pytest.fail(f"[{geometry}] apply_metric differs from the inline loop")
    if not np.array_equal(space.apply_inverse_metric(x), _inline_per_axis(axes, x, inverse=True)):
        pytest.fail(f"[{geometry}] apply_inverse_metric differs from the inline loop")
    if space.inner_product(x, y) != float(np.sum(_inline_per_axis(axes, x, inverse=False) * y)):
        pytest.fail(f"[{geometry}] inner_product differs from the inline loop")


# ═════════════════════════════════════════════════════════════════════════
# Leg 2 — a POSITIONED dense block is admitted and APPLIED
# ═════════════════════════════════════════════════════════════════════════

def _dense_gram_on_the_head() -> tuple[tuple[Axis, ...], np.ndarray]:
    """Axes whose head block carries NO diagonal measure (a dense Gram has
    none) plus a symmetric positive-definite 6×6 form for that block."""
    rng = np.random.default_rng(11)
    a = rng.standard_normal((6, 6))
    gram = a @ a.T + 6.0 * np.eye(6)
    axes = (
        Axis(label="harmonic_probe", shape=(2, 3), origin=Declared(Modal(np.eye(1, 6)[0].reshape(2, 3)))),
        Axis(label="spatial", shape=(4,), origin=Declared(BasisKind.NODAL, DiagonalMetric(np.array([0.5, 1.0, 1.5, 2.0])))),
    )
    return axes, gram


def _dense_head_axis(axes, gram):
    """The head axis of ``axes`` re-minted with the dense form ``gram`` as its form (a coefficient axis, Modal)."""
    return Axis(label=axes[0].label, shape=axes[0].shape, origin=Declared(axes[0].kind, DenseMetric(gram, kernel=Coord(6, ()), face=Face.METRIC)))


def test_a_positioned_dense_block_is_admitted_and_applied():
    """The guard admits ``axes`` + a FactoredMetric positioned over them
    (the dense form on the measure-less head block, ``None`` on the
    spatial block), and ``apply_metric`` APPLIES the dense block on the
    head's flattened block (oracle: an explicit matmul on the reshaped
    array) while the spatial axis's OWN measure — never restated in the
    object — still applies on its block: the overlay is merged into the
    derived entries, not substituted for them.

    NEGATIVE CONTROL: the same axes with NO object (the counting-axis route
    the old short-circuit would have taken) differ at O(1)."""
    axes, gram = _dense_gram_on_the_head()
    spatial_w = axes[1].form.weights
    assert spatial_w is not None
    space = TensorSpace.of_axes(
        Factor.whole_on(_dense_head_axis(axes, gram), Kind.MEASURE, VarianceChart.DENSITY),
        Factor.whole_on(axes[1], Kind.MEASURE, VarianceChart.DENSITY),
    )

    x = np.random.default_rng(3).standard_normal((2, 3, 4))
    got = space.apply_metric(x)
    # Oracle: G on the flattened head block, then the spatial diagonal.
    want = (gram @ x.reshape(6, 4)).reshape(2, 3, 4) * spatial_w.reshape(1, 1, 4)
    np.testing.assert_allclose(got, want, rtol=1e-13, atol=0.0)
    inv = space.apply_inverse_metric(got)
    np.testing.assert_allclose(inv, x, rtol=1e-11, atol=1e-13)

    counting_route = TensorSpace.of_axes(*(Factor.whole_on(_a, Kind.MEASURE, VarianceChart.DENSITY) for _a in axes)).apply_metric(x)
    separation = float(np.abs(got - counting_route).max()) / float(np.abs(got).max())
    if separation < 1e-2:
        pytest.fail(
            f"CONTROL FAILED: the positioned dense block moved apply_metric by only "
            f"{separation:.3e} relative against the counting-axis route"
        )


def test_a_product_of_axis_built_factors_carries_the_head_form_beside_its_axes():
    """``head * bulk`` with a dense-Gram head (the 6.2c-ii occupant: the
    moment space of a DENSE-Gram frame): the product is axis-built (axes
    concatenated), carries the head's positioned form beside them, and
    applies it — the dense block on the head, the spatial measure on its
    own block, identity on the counting energy block. Oracle: the explicit
    matmul on the reshaped array."""
    axes, gram = _dense_gram_on_the_head()
    spatial_w = axes[1].form.weights
    assert spatial_w is not None
    head_axis = axes[0]
    head = TensorSpace.of_axes(Factor.whole_on(_dense_head_axis(axes, gram), Kind.MEASURE, VarianceChart.DENSITY))
    bulk = TensorSpace.of_axes(Factor.whole_on(EnergyAxis(label="energy", shape=(2,)), Kind.MEASURE, VarianceChart.DENSITY), Factor.whole_on(axes[1], Kind.MEASURE, VarianceChart.DENSITY))
    product = head * bulk
    assert len(product.axes) == 3
    if not isinstance(product.factors[0].form, DenseMetric):
        pytest.fail("the product dropped the head's dense form")

    y = np.random.default_rng(5).standard_normal((2, 3, 2, 4))
    got = product.apply_metric(y)
    want = (gram @ y.reshape(6, 8)).reshape(2, 3, 2, 4) * spatial_w.reshape(1, 1, 1, 4)
    np.testing.assert_allclose(got, want, rtol=1e-13, atol=0.0)
    # And the product reads the same pairing its factors do — the head's
    # form and the bulk's measures, one spelling.
    x_head = y[:, :, 0, 0]
    np.testing.assert_allclose(
        head.inner_product(x_head, x_head), float(x_head.ravel() @ gram @ x_head.ravel()),
        rtol=1e-13,
    )
