r"""The dense-metric family's own laws (campaign 1, P7 — group A of the
battery of record, ``scratch/p7_verification_plan.md`` §2).

L0 software invariants of a numerics primitive with no theory-page
``:label:`` — ``foundation``-marked, no ``verifies()`` markers.

⭐ The references are **hand-derived arithmetic in exact binary fractions**
(halves and quarters), so the value gates use ``==`` / ``array_equal``, not
tolerances. This is the ⛔⛔-mandated *independent* pin: reciprocity cannot
adjudicate a metric (``A† ≡ G⁻¹AᵀG`` is an identity for EVERY invertible
``G``), so the metric's correctness evidence is literals with no ORPHEUS
solver, frame, or quadrature in the chain (plus the wrong-metric
discriminator in ``test_frame.py``, which lands with the S3 dressing).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.metric import (
    _DENSE_METRIC_SYMMETRY_RTOL,
    DenseMetric,
    DiagonalMetric,
)
from orpheus.numerics.quadrature.directional import Quadrature
from orpheus.numerics.space import FunctionSpace

pytestmark = pytest.mark.foundation


def _require(condition: bool, message: str) -> None:
    """A ``-O``-firing assertion (NOT a bare ``assert``)."""
    if not condition:
        pytest.fail(message)


#: Exact-binary-fraction fixture: every entry a multiple of 1/4, so G@x and
#: yᵀGx are exact in IEEE double and the gates below may use ``==``.
_G_LITERAL = np.array(
    [
        [2.0, 0.5, 0.0],
        [0.5, 1.0, 0.25],
        [0.0, 0.25, 4.0],
    ]
)
_X = np.array([1.0, 2.0, 3.0])
_Y = np.array([0.5, -1.0, 2.0])


class TestDenseMetricLaws:
    def test_dense_metric_pairing_matches_hand_derived_literals(self):
        """A1 — G@x and ⟨x,y⟩ against hand-derived exact literals.

        By hand: G@x = (2+1+0, 1/2+2+3/4, 0+1/2+12) = (3, 13/4, 25/2) and
        yᵀGx = 3/2 − 13/4 + 25 = 93/4 = 23.25 — all exact binary fractions,
        so the comparison is ``==``, not a tolerance.
        """
        metric = DenseMetric(_G_LITERAL)
        _require(
            bool(np.array_equal(metric.apply(_X), [3.0, 3.25, 12.5])),
            f"G@x != hand literal: {metric.apply(_X)!r}",
        )
        space = FunctionSpace("dense3", (3,), metric=metric)
        got = space.inner_product(_X, _Y)
        _require(
            got == 23.25,
            f"<x,y> = {got!r} != the hand-derived 23.25",
        )

    def test_dense_metric_matrix_and_pairing_are_symmetric(self):
        """A2 — symmetry of the INSTALLED matrix, on the flagship Gram.

        The matrix leg is the load-bearing one: [M] 2026-08-30,
        ``np.linalg.pinv`` WITHOUT ``hermitian=True`` gives
        ``max|M−Mᵀ| = 4.74e-14`` on this Gram while the pairing asymmetry
        stays at ``1.5e-15`` — a pairing-only gate is BLIND to the
        spelling choice (battery arm M3's catcher is THIS matrix leg).
        """
        frame = Quadrature.gauss_legendre(8).angular_frame(2)
        metric = DenseMetric.inverse_of(frame.discrete_gram)
        asym = float(np.max(np.abs(metric.matrix - metric.matrix.T)))
        _require(
            asym <= 1e-15,
            f"installed matrix asymmetry {asym:.3e} > 1e-15",
        )
        rng = np.random.default_rng(20260830)
        x = rng.standard_normal(metric.dim)
        y = rng.standard_normal(metric.dim)
        space = FunctionSpace("slab_par", (metric.dim,), metric=metric)
        _require(
            space.inner_product(x, y)
            == pytest.approx(space.inner_product(y, x), rel=1e-14),
            "pairing is not symmetric",
        )

    def test_dense_metric_refuses_an_asymmetric_matrix(self):
        """A3 — the symmetry guard's negative leg (an asymmetric form is
        not an inner product; the producer symmetrizes, the type refuses).
        The threshold quoted is the type's own
        (``_DENSE_METRIC_SYMMETRY_RTOL``, vv-principles #16)."""
        _require(
            _DENSE_METRIC_SYMMETRY_RTOL == 1e-12,
            "the gate's docstring quotes a stale threshold",
        )
        with pytest.raises(ValueError, match="requires a symmetric matrix"):
            DenseMetric(np.array([[1.0, 2.0], [0.0, 1.0]]))

    def test_dense_metric_refuses_a_matrix_that_does_not_span_the_space_shape(self):
        """A4 — the admission check at space construction: a (3,3) matrix
        cannot serve a (4,)-shaped space (modelled on the axes' own shape
        refusal, so the vocabulary stays greppable)."""
        with pytest.raises(ValueError, match="cannot serve a space"):
            FunctionSpace("wrong", (4,), metric=DenseMetric(_G_LITERAL))

    def test_dense_inverse_metric_is_the_range_projector_not_the_identity(self):
        """A5 — the Moore–Penrose doctrine, hand-derived and exact.

        ``S = [[1,1,0],[1,1,0],[0,0,2]]`` has rank 2; ``S⁺S`` is the
        range projector ``[[.5,.5,0],[.5,.5,0],[0,0,1]]``. ON range the
        round trip is the identity; OFF range (the kernel direction
        ``(1,−1,0)``) it annihilates — ``max|Px − x| = 1`` exactly. The
        off-range leg is what makes the on-range leg informative
        (vv-principles #19): the separation is ~4.5e15×.
        """
        s = DenseMetric(np.array([[1.0, 1.0, 0.0], [1.0, 1.0, 0.0], [0.0, 0.0, 2.0]]))
        on_range = np.array([1.0, 1.0, 0.0])
        off_range = np.array([1.0, -1.0, 0.0])
        round_on = s.apply_inverse(s.apply(on_range))
        round_off = s.apply_inverse(s.apply(off_range))
        _require(
            float(np.max(np.abs(round_on - on_range))) <= 1e-15,
            f"on-range round trip broke: {round_on!r}",
        )
        _require(
            float(np.max(np.abs(round_off - off_range)))
            == pytest.approx(1.0, abs=1e-15),
            f"off-range component was not annihilated: {round_off!r}",
        )

    def test_dense_metric_of_a_diagonal_matrix_reduces_to_the_hadamard_apply(self):
        """A6 — ``DenseMetric(diag(w)).apply`` equals the Hadamard apply
        bit-for-bit. [M] seed sweep at design time: 0 of 2000 draws differ
        at n = 3 / 15 / 225 — ``array_equal`` is a property of the
        arithmetic (a matmul row with one nonzero is that one product),
        not of a lucky draw."""
        rng = np.random.default_rng(7)
        for n in (3, 15, 225):
            w = rng.standard_normal(n) ** 2 + 0.1
            x = rng.standard_normal(n)
            dense = DenseMetric(np.diag(w)).apply(x)
            hadamard = DiagonalMetric(w).apply(x)
            _require(
                bool(np.array_equal(dense, hadamard)),
                f"n={n}: diag-matrix apply differs from Hadamard",
            )

    def test_dense_metric_of_a_diagonal_matrix_reduces_to_the_hadamard_pairing(self):
        """A7 — the same reduction on the PAIRING face, at rtol=1e-12 and
        deliberately NOT ``array_equal``: [M] design-time sweep,
        ``sum((w·x)·y)`` vs ``y@(diag(w)@x)`` differ on 1360 of 2000 draws
        at n=15, worst 1792 ULP (rel 2.0e-13) — the two reduction trees
        are different, so a nulp pin taken from one green draw would be
        seed-fragile. This is also why the space's pairing is spelled
        ``sum(metric.apply(x) * y)`` and never densified (arm M13)."""
        rng = np.random.default_rng(11)
        for n in (3, 15, 225):
            w = rng.standard_normal(n) ** 2 + 0.1
            x = rng.standard_normal(n)
            y = rng.standard_normal(n)
            dense = FunctionSpace(
                f"d{n}", (n,), metric=DenseMetric(np.diag(w))
            ).inner_product(x, y)
            hadamard = FunctionSpace(
                f"h{n}", (n,), inner_product_weights=w
            ).inner_product(x, y)
            _require(
                dense == pytest.approx(hadamard, rel=1e-12),
                f"n={n}: {dense!r} vs {hadamard!r}",
            )

    @pytest.mark.parametrize("shape", [(4,), (3, 4), (2, 3, 2)])
    def test_dense_metric_apply_equals_an_independent_einsum_reference(self, shape):
        """A8 — the flattening/reshape convention pinned against a second
        spelling: an explicit C-order ravel + einsum. The metric's block
        is row-major over the space shape — the same order as
        ``FrameBase.discrete_gram``'s ``table.reshape(N, K)``."""
        rng = np.random.default_rng(sum(shape))
        k = int(np.prod(shape))
        a = rng.standard_normal((k, k))
        g = a.T @ a + np.eye(k)
        x = rng.standard_normal(shape)
        got = DenseMetric(g).apply(x)
        ref = np.einsum("jk,k->j", g, x.ravel(order="C")).reshape(shape)
        _require(
            bool(np.allclose(got, ref, rtol=1e-14, atol=0.0)),
            f"shape {shape}: apply disagrees with the einsum reference",
        )

    def test_a_dense_metric_space_keeps_metric_blind_equality_and_hash(self):
        """A9 — Ruling 2's identity clause: the metric slot is
        ``compare=False`` metadata, so two spaces differing only in metric
        are equal, hash-equal, and interchangeable as dict keys. ⚠ The
        constraint is structural, not taste: [M] a COMPARED metric field
        makes the dataclass ``__eq__`` of subclasses return an ndarray
        (truthiness-ambiguous) and ``hash()`` raise — the CS5
        generator-slot lesson recurring (battery arm M12).
        """
        bare = FunctionSpace("same", (3,))
        dressed = FunctionSpace("same", (3,), metric=DenseMetric(_G_LITERAL))
        _require((bare == dressed) is True, "metric moved space equality")
        _require(hash(bare) == hash(dressed), "metric moved the space hash")
        _require({bare: 1}[dressed] == 1, "metric broke dict-key interchange")


class TestDenseMetricPropagation:
    """Group C (partial — C1/C2): the sites that used to DROP a dense
    metric silently. C3 (the frame's gram probe) lives beside the frame
    gates; C4/C5 land with the S3 dressing."""

    def test_tensor_product_with_a_dense_metric_factor_does_not_go_silently_euclidean(self):
        """C1 — the product carries G ⊗ diag(w), lazily positioned.

        Hand-derived on separable probes: with x = xv ⊗ xw, y = yv ⊗ yw,
        ⟨x,y⟩_{G⊗w} = ⟨xv,yv⟩_G · ⟨xw,yw⟩_w = 23.25 · (2·1·1 + 0.5·2·0.5)
        = 23.25 · 2.5 = 58.125 — exact binary fractions, compared with
        ``==``. [M] the pre-repair behaviour dropped the dense factor
        entirely (Euclidean on that block); battery arm M8 is the revert.
        """
        v = FunctionSpace("V", (3,), metric=DenseMetric(_G_LITERAL))
        w = FunctionSpace(
            "W", (2,), inner_product_weights=np.array([2.0, 0.5])
        )
        product = v * w
        x = np.multiply.outer(_X, np.array([1.0, 2.0]))
        y = np.multiply.outer(_Y, np.array([1.0, 0.5]))
        got = product.inner_product(x, y)
        _require(
            got == 58.125,
            f"tensor-product pairing {got!r} != the hand-derived 58.125 "
            f"(the dense factor was dropped or mis-positioned)",
        )

    def test_the_dual_of_a_dense_metric_space_carries_the_same_pairing(self):
        """C2 — L²-Riesz: the dual carries the SAME metric as the primal.
        [M] the pre-repair ``DualSpace.of`` read the plain Euclidean
        pairing on a dense-metric primal (4.5 where the primal reads
        23.3-class values); battery arm M7 is the revert."""
        v = FunctionSpace("V", (3,), metric=DenseMetric(_G_LITERAL))
        dual = v.dual()
        _require(
            dual.inner_product(_X, _Y) == v.inner_product(_X, _Y) == 23.25,
            "the dual's pairing must equal the primal's",
        )


class TestPairingSpelling:
    def test_the_pairing_keeps_the_legacy_reduction_tree_bit_exactly(self):
        """A10 (battery arm M13's catcher — [M] the battery found every
        tolerance gate BLIND to a densified pairing, deviation ~2e-13).

        The diagonal pairing is bit-identical to the legacy
        ``np.sum(w * x * y)`` spelling — the guarantee that made the P7
        reroute invisible to every pinned number in the tree, promoted
        from a docstring claim to an ``==`` witness. [M] seed 1, n=15:
        the densified ``y @ (diag(w) @ x)`` spelling differs on this
        draw (−3.382498612232062 vs −3.3824986122320624), so this gate
        reds under the forbidden spelling while every rtol gate stays
        green.
        """
        rng = np.random.default_rng(1)
        n = 15
        w = rng.standard_normal(n) ** 2 + 0.1
        x = rng.standard_normal(n)
        y = rng.standard_normal(n)
        space = FunctionSpace("w15", (n,), inner_product_weights=w)
        legacy = float(np.sum(w * x * y))
        _require(
            space.inner_product(x, y) == legacy,
            "the pairing left the legacy reduction tree",
        )
