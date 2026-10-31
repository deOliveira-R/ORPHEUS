r"""The kernel gauge's laws on a numerics fixture: :class:`~orpheus.numerics.kernel_projection_gauge.KernelProjectionGauge`.

The projector :math:`\Pi = \sum_b \gamma_b^{\mathsf T} R_b G_b^{-1} M_b \gamma_b` onto a kernel given by spanning
summands, and its section :math:`\sigma(\psi) = \psi - \Pi\psi`, the minimum-:math:`G`-norm member of the coset
(space-membership campaign U3, plan ``space_membership_covector_gauge.md``; spec
``scratch/membership_covector/spec/spec.md`` §U3.2 and §R3). SN's loss kernel is one consumer, gated on its production
instances in ``tests/gates/sn/operators/test_loss_kernel_gauge.py``; these rows are the laws' numerics witness.

The fixture's metric is NON-uniform (max/min = 16, premise row): under a uniform :math:`G` the Euclidean projector is
the :math:`G`-orthogonal one and the Euclidean transpose is :math:`\Pi` itself, so the self-adjointness, minimum-norm
and transpose rows would be blind (a Mode-12 stabiliser of the space).

Teeth: every row names the arm that reddens it. `[M]` 2026-10-31, the in-process battery
``scratch/membership_covector/u3/battery/kpg_arms.py`` (a ``-p`` plugin rebinding its target at ``pytest_configure``),
one arm per mutation, under ``python -O -m pytest``.
"""
from __future__ import annotations

import inspect

import numpy as np
import pytest

from orpheus.numerics.axis import Axis, BasisKind, Declared
from orpheus.numerics.gauge import KernelGauge
from orpheus.numerics.kernel_projection_gauge import (
    SUMMAND_ORTHOGONALITY_TOL,
    KernelProjectionGauge,
    KernelSummand,
)
from orpheus.numerics.space import FunctionSpace

pytestmark = [
    pytest.mark.l0,
    pytest.mark.rests_on(
        "tests/gates/numerics/test_coordinate_space.py::TestThePrimitive",
        "tests/gates/numerics/test_coordinate_space.py::TestTheSpaceVerbsParseOnEveryRun",
    ),
]

_EPS = np.finfo(float).eps
_N = 9
_G = np.array([0.5, 1.0, 2.0, 0.25, 3.0, 1.5, 0.75, 4.0, 1.25])
SPACE = FunctionSpace.of_axes(Axis("dof", (_N,), weights=_G, origin=Declared(BasisKind.NODAL)))

#: Coordinate-aligned kernel K = {2, 5, 7}: A = A_RR (+) 0_K, so A e_t == 0 is bitwise (L90: gate the object).
_K = np.array([2, 5, 7])
_R = np.setdiff1d(np.arange(_N), _K)
_A = np.zeros((_N, _N))
_A[np.ix_(_R, _R)] = np.random.default_rng(20261030).standard_normal((_R.size, _R.size)) + 4.0 * np.eye(_R.size)

#: Two summands with NON-orthonormal spanning columns on disjoint supports (a legal input: numerics orthonormalises).
_ROTATED = ((np.array([0, 1, 2]), np.array([[1.0], [-2.0], [0.5]])), (np.array([5, 6]), np.array([[0.3], [1.0]])))

#: Over-complete: three columns spanning a rank-2 subspace. `[M]` sqrt(G)-weighted s/s0 = [1, 0.598, 2.9e-17]
#: (``scratch/membership_covector/u3/rank_gap.py``), so the rank cut ``SPAN_RANK_RTOL`` decides dimension 2, not 3.
_V1, _V2 = np.array([1.0, -2.0, 0.5, 0.0]), np.array([0.0, 0.7, 1.0, -1.3])
_OVERCOMPLETE = ((np.array([0, 1, 2, 3]), np.column_stack([_V1, _V2, 0.3 * _V1 - 1.7 * _V2])),)

#: Overlapping supports, G-orthogonal exactly in floats: on the overlap {2, 3}, G2*0.5*1 + G3*1*(-4) = 1 - 1 = 0.
_ORTHOGONAL_OVERLAP = ((np.array([0, 1, 2, 3]), np.array([[1.0], [-2.0], [0.5], [1.0]])),
                       (np.array([2, 3, 4]), np.array([[1.0], [-4.0], [0.7]])))
#: Overlapping, not orthogonal: both summands contain e_1.
_OVERLAP_NOT_ORTHOGONAL = ((np.array([0, 1]), np.eye(2)), (np.array([1, 2]), np.eye(2)))
_OVERLAP_REFUSAL = "not mutually G-orthogonal"


def _summands(pairs):
    return tuple(KernelSummand(support=idx, generators=cols, label=f"kpg_fixture_{k}") for k, (idx, cols) in enumerate(pairs))


def _build(space, pairs):
    return KernelProjectionGauge.direct_sum(space, _summands(pairs))


def _embed(idx, col):
    v = np.zeros(_N)
    v[idx] = col
    return v


def _ip(x, y):
    """The G inner product by hand: the oracle never asks the SUT's space (X4)."""
    return float(np.sum(_G * x * y))


def _x(seed):
    return np.random.default_rng(seed).standard_normal(_N)


@pytest.fixture
def aligned():
    return _build(SPACE, ((_K, np.eye(_K.size)),))


@pytest.fixture
def rotated():
    return _build(SPACE, _ROTATED)


@pytest.fixture
def overcomplete():
    return _build(SPACE, _OVERCOMPLETE)


def test_premise_the_metric_is_not_uniform():
    assert _G.max() / _G.min() > 10.0


def test_premise_the_axis_built_space_and_its_legacy_twin_apply_one_metric():
    """Localises a red of the twin row: green here and red there means the gauge reads the metric's SPELLING."""
    twin = FunctionSpace("twin", (_N,), inner_product_weights=_G)
    x = _x(1)
    assert np.array_equal(SPACE.apply_metric(x), twin.apply_metric(x))


# ── rung 1: the projector laws ─────────────────────────────────────────────────────────────────────────

class TestTheProjectorLaws:
    @pytest.mark.parametrize("which", ["aligned", "rotated", "overcomplete"])
    def test_idempotence_of_the_projector_and_of_the_section(self, which, request):
        """Arm ``leak`` (Pi x + 1e-3 x) reddens all three."""
        g = request.getfixturevalue(which)
        x = _x(2)
        atol = 16 * _EPS * np.abs(x).max()
        np.testing.assert_allclose(g.apply(g.apply(x)), g.apply(x), rtol=0, atol=atol)
        np.testing.assert_allclose(g.gauge(g.gauge(x)), g.gauge(x), rtol=0, atol=atol)

    @pytest.mark.parametrize("which", ["aligned", "rotated", "overcomplete"])
    def test_G_self_adjointness(self, which, request):
        """Arm ``euclid`` (the metric read as ones) reddens 'rotated' and 'overcomplete'; 'aligned' is inside that
        arm's stabiliser (coordinate columns are G-orthogonal under every diagonal G): the control leg."""
        g = request.getfixturevalue(which)
        x, y = _x(3), _x(4)
        assert _ip(g.apply(x), y) == pytest.approx(_ip(x, g.apply(y)), rel=1e-13)

    @pytest.mark.parametrize("which", ["rotated", "overcomplete"])
    def test_the_euclidean_transpose_is_not_the_projector(self, which, request):
        r""":math:`\langle \Pi x, y\rangle = \langle x, \Pi^{\mathsf T} y\rangle` in the EUCLIDEAN pairing. Arm
        ``transpose_is_apply`` (``apply_transpose`` returns ``apply``, SN's predecessor) reddens both: on a
        non-uniform G, :math:`\Pi^{\mathsf T} = G\Pi G^{-1} \ne \Pi` (asserted below, so the row is not blind)."""
        g = request.getfixturevalue(which)
        x, y = _x(5), _x(6)
        assert float(g.apply(x) @ y) == pytest.approx(float(x @ g.apply_transpose(y)), rel=1e-13)
        assert np.abs(g.apply_transpose(y) - g.apply(y)).max() > 1e-2

    @pytest.mark.parametrize("pairs, dim", [(_ROTATED, 2), (_OVERCOMPLETE, 2)], ids=["rotated", "overcomplete"])
    def test_the_range_is_the_span_of_the_given_columns_and_the_dimension_is_its_rank(self, pairs, dim):
        """Every given column (the over-complete one included) is fixed by Pi, and dimension = rank, not the column
        count. Arm ``no_rank_cut`` (the cut at 0) reddens 'overcomplete' (dimension 3); arm ``drop_summand``
        reddens 'rotated'."""
        g = _build(SPACE, pairs)
        for idx, cols in pairs:
            for c in cols.T:
                v = _embed(idx, c)
                np.testing.assert_allclose(g.apply(v), v, rtol=0, atol=16 * _EPS * np.abs(v).max())
        assert g.dimension == dim

    def test_the_section_is_G_orthogonal_to_the_kernel(self, rotated):
        """Arm ``gauge_sign`` (sigma = x + Pi x) reddens it."""
        x = _x(7)
        for idx, cols in _ROTATED:
            v = _embed(idx, cols[:, 0])
            assert abs(_ip(rotated.gauge(x), v)) < 1e-14 * np.sqrt(_ip(x, x) * _ip(v, v))

    def test_the_section_is_the_minimum_G_norm_member_of_the_coset(self, rotated):
        """No shift along a kernel vector shortens the section: the best shift ``t* = -<s, v>_G / <v, v>_G`` is the
        hand-computed minimiser, so a section off the minimum is shortened by it. Arm ``gauge_sign`` reddens it.
        (A fixed set of shifts does not: `[M]` the first draft's t in {-1, 0.3, 2} stayed green under that arm.)"""
        x = _x(8)
        s = rotated.gauge(x)
        for idx, cols in _ROTATED:
            v = _embed(idx, cols[:, 0])
            t_star = -_ip(s, v) / _ip(v, v)
            shifted = s + t_star * v
            assert _ip(s, s) <= _ip(shifted, shifted) * (1.0 + 64 * _EPS)
            for t in (-1.0, 0.3, 2.0):
                assert _ip(s, s) < _ip(s + t * v, s + t * v)


# ── rung 1: the summand basis that earns the DIAGONAL declaration (relocated from the SN gauge suite at U3) ──

class TestTheSummandBasis:
    @pytest.mark.parametrize("pairs", [_ROTATED, _OVERCOMPLETE, _ORTHOGONAL_OVERLAP], ids=["rotated", "overcomplete", "overlap"])
    def test_each_summand_gram_is_the_identity_which_is_what_earns_DIAGONAL(self, pairs):
        r""":math:`\Phi^{\mathsf T} G \Phi = I` per summand, on generators that are NOT G-orthogonal ('overcomplete':
        three columns, two modes). The frame's ``CrossGramInverse`` reads a DIAGONAL Gram by a row-sum probe, which is
        the true diagonal only if the Gram IS diagonal; arm ``euclid`` (the SVD unweighted) reddens 'rotated' and
        'overcomplete'. The SN production instances (d=3, 43 % off-diagonal raw generators) stay in
        ``test_loss_kernel_gauge.py``."""
        from orpheus.numerics.basis import GramStructure
        g = _build(SPACE, pairs)
        saw_multimode = False
        for s in g.summands:
            table = s.basis.table
            assert s.basis.gram_structure is GramStructure.DIAGONAL
            gram = table.T @ (_G[s.gather.indices][:, None] * table)
            saw_multimode |= table.shape[1] > 1
            np.testing.assert_allclose(gram, np.eye(table.shape[1]), rtol=0, atol=64 * _EPS)
        assert saw_multimode or pairs is not _OVERCOMPLETE

    def test_the_frame_gram_probe_returns_ones(self, overcomplete):
        """Closes the loop between the declaration and the machinery that consumes it: the value
        ``CrossGramInverse.of(frame)`` computes on the summand's frame."""
        from orpheus.numerics import frame_faces
        from orpheus.numerics.frame import GalerkinFrame
        from orpheus.numerics.measure import DiscreteMeasure
        for s in overcomplete.summands:
            idx = s.gather.indices
            frame = GalerkinFrame(s.basis, DiscreteMeasure(nodes=idx.astype(float), weights=_G[idx], support=s.basis.domain))
            diagonal = np.asarray(frame_faces.CrossGramInverse.of(frame).diagonal)
            np.testing.assert_allclose(diagonal, np.ones_like(diagonal), rtol=0, atol=64 * _EPS)

    def test_a_computed_unit_is_the_projection_of_the_constant_and_a_declared_one_is_zero(self):
        """The basis always computes ``Phi^T G 1``; a summand DECLARED orthogonal to the constants is checked and then
        stores the exact zero. The declared span on {2, 3}: G2*1 + G3*(-8) = 2 - 2 = 0 exactly in floats."""
        idx, cols = _ROTATED[0]
        computed = KernelProjectionGauge.direct_sum(SPACE, (KernelSummand(idx, cols, "u"),)).summands[0].basis
        np.testing.assert_allclose(computed.layout.unit, computed.table.T @ _G[idx], rtol=0, atol=8 * _EPS)
        assert np.abs(computed.layout.unit).max() > 0.1
        declared = KernelProjectionGauge.direct_sum(SPACE, (KernelSummand(
            np.array([2, 3]), np.array([[1.0], [-8.0]]), "d", orthogonal_to_constants=True),)).summands[0].basis
        assert declared.layout.unit.shape == (1,) and not np.any(declared.layout.unit)

    def test_a_false_declaration_of_orthogonality_to_the_constants_is_refused(self):
        """Arm ``no_declaration_check`` reddens it (the false zero would be stored)."""
        idx, cols = _ROTATED[0]
        with pytest.raises(ValueError, match="declared orthogonal to the constants"):
            KernelProjectionGauge.direct_sum(SPACE, (KernelSummand(idx, cols, "f", orthogonal_to_constants=True),))

    def test_a_zero_metric_support_row_is_refused(self):
        zero = FunctionSpace("zero_row", (3,), inner_product_weights=np.array([1.0, 0.0, 2.0]))
        with pytest.raises(ValueError, match="carry a vanishing metric"):
            KernelProjectionGauge.direct_sum(zero, (KernelSummand(np.array([0, 1]), np.eye(2), "z"),))


# ── rung 2: residual neutrality against an operator whose kernel is the span ───────────────────────────

@pytest.mark.rests_on("tests/gates/numerics/test_kernel_projection_gauge.py::TestTheProjectorLaws")
class TestResidualNeutrality:
    def test_the_kernel_columns_are_annihilated_bitwise(self):
        for t in _K:
            e = np.zeros(_N)
            e[t] = 1.0
            assert np.array_equal(_A @ e, np.zeros(_N))

    def test_the_section_cannot_move_the_residual(self, aligned):
        """Arm ``leak`` reddens it (the leaked 1e-3 x is not in ker A)."""
        x = _x(9)
        np.testing.assert_allclose(_A @ aligned.gauge(x), _A @ x, rtol=0, atol=64 * _EPS * np.abs(_A @ x).max())

    def test_negative_a_projector_off_the_kernel_is_not_residual_neutral(self):
        """Fixture-loading control, green under every arm: the same object on a span NOT in ker A moves A x."""
        off = _build(SPACE, ((_R[:2], np.eye(2)),))
        x = _x(10)
        assert np.abs(_A @ off.gauge(x) - _A @ x).max() > 1e-2


# ── rung 3: the direct sum ─────────────────────────────────────────────────────────────────────────────

@pytest.mark.rests_on("tests/gates/numerics/test_kernel_projection_gauge.py::TestTheProjectorLaws")
class TestTheDirectSum:
    def test_overlapping_G_orthogonal_summands_are_admitted_and_sum_to_a_projector(self):
        """R3: the SN diamond-difference modes and the physical constant mode overlap on the trace. Positive leg."""
        g = _build(SPACE, _ORTHOGONAL_OVERLAP)
        x, y = _x(11), _x(12)
        assert g.dimension == 2
        np.testing.assert_allclose(g.apply(g.apply(x)), g.apply(x), rtol=0, atol=16 * _EPS * np.abs(x).max())
        assert _ip(g.apply(x), y) == pytest.approx(_ip(x, g.apply(y)), rel=1e-13)

    def test_overlapping_summands_that_are_NOT_G_orthogonal_are_refused(self):
        """Negative leg. Arm ``no_orthogonality_check`` reddens this row and not the positive leg above, which is
        what tells the check from the law."""
        with pytest.raises(ValueError, match=_OVERLAP_REFUSAL):
            _build(SPACE, _OVERLAP_NOT_ORTHOGONAL)

    @pytest.mark.parametrize("cosine, admitted", [(1e-17, True), (1e-6, False)], ids=["R3-union-reading", "1e-6"])
    def test_the_orthogonality_tolerance_separates_round_off_from_a_real_overlap(self, cosine, admitted):
        """`[M]` ``scratch/membership_covector/u3/union_orthogonality.py``: in the R3 union case (the DD (3,4) box at
        Sigma_a = 0) and three siblings, every DD mode's G-cosine with the physical constant mode is <= 1e-17. Arm
        ``no_orthogonality_check`` reddens the 1e-6 leg."""
        assert 1e-17 < SUMMAND_ORTHOGONALITY_TOL < 1e-6
        v = _embed(_ORTHOGONAL_OVERLAP[0][0], _ORTHOGONAL_OVERLAP[0][1][:, 0])
        w = _embed(_ORTHOGONAL_OVERLAP[1][0], _ORTHOGONAL_OVERLAP[1][1][:, 0])
        tilted = w + cosine * np.sqrt(_ip(w, w) / _ip(v, v)) * v
        support = np.arange(5)
        pairs = (_ORTHOGONAL_OVERLAP[0], (support, tilted[support, None]))
        if admitted:
            assert _build(SPACE, pairs).dimension == 2
        else:
            with pytest.raises(ValueError, match=_OVERLAP_REFUSAL):
                _build(SPACE, pairs)

    @pytest.mark.parametrize("cut", [1e-14, 1e-12, 1e-10, 1e-8])
    def test_the_verdicts_do_not_depend_on_the_cut_anywhere_in_the_claimed_range(self, cut, monkeypatch):
        """Gates the production comment's claim that the verdict is the same for any cut in [1e-14, 1e-8]: the
        honest overlap (exactly orthogonal) and a 1e-17 tilt are admitted, a 1e-6 tilt refused, at every cut. The SN
        half of the claim (the DD modes against the constant, the R3 union) is the SN row
        ``test_the_kernel_modes_are_orthogonal_to_the_constants`` (<= 64 eps)."""
        import orpheus.numerics.kernel_projection_gauge as kpg
        monkeypatch.setattr(kpg, "SUMMAND_ORTHOGONALITY_TOL", cut)
        v = _embed(_ORTHOGONAL_OVERLAP[0][0], _ORTHOGONAL_OVERLAP[0][1][:, 0])
        w = _embed(_ORTHOGONAL_OVERLAP[1][0], _ORTHOGONAL_OVERLAP[1][1][:, 0])
        support = np.arange(5)
        assert _build(SPACE, _ORTHOGONAL_OVERLAP).dimension == 2
        for cosine, admitted in ((1e-17, True), (1e-6, False)):
            tilted = w + cosine * np.sqrt(_ip(w, w) / _ip(v, v)) * v
            pairs = (_ORTHOGONAL_OVERLAP[0], (support, tilted[support, None]))
            if admitted:
                assert _build(SPACE, pairs).dimension == 2
            else:
                with pytest.raises(ValueError, match=_OVERLAP_REFUSAL):
                    _build(SPACE, pairs)

    def test_summand_labels_must_be_distinct(self):
        """Each label names its support's manifold; two equal labels would make two coefficient spaces compare equal.
        Arm ``no_label_check`` reddens it."""
        twins = (KernelSummand(np.array([0, 1]), np.eye(2), "same"), KernelSummand(np.array([5, 6]), np.eye(2), "same"))
        with pytest.raises(ValueError, match="labels must be distinct"):
            KernelProjectionGauge.direct_sum(SPACE, twins)

    def test_a_support_outside_the_space_is_refused(self):
        with pytest.raises(ValueError, match="outside the space's"):
            KernelProjectionGauge.direct_sum(SPACE, (KernelSummand(np.array([7, 9]), np.eye(2), "o"),))

    def test_an_empty_direct_sum_is_the_identity_section_not_None(self):
        g = _build(SPACE, ())
        x = _x(13)
        assert g.dimension == 0 and np.array_equal(g.gauge(x), x)

    def test_a_rank_zero_summand_contributes_nothing(self):
        g = _build(SPACE, ((np.array([0, 1]), np.zeros((2, 1))), _ROTATED[1]))
        assert g.dimension == 1 and len(g.summands) == 1

    def test_no_inverse_and_it_is_spelled_by_absence(self, rotated):
        """Relocated from the SN file (``test_no_inverse_and_it_is_spelled_by_ABSENCE``): a property of the type."""
        assert rotated.is_invertible is False and not hasattr(rotated, "inverse")
        assert rotated.is_adjointable is True and rotated.domain is rotated.codomain


# ── rung 3: the Protocol, the metric and the parse ─────────────────────────────────────────────────────

class TestTheProtocolAndTheMetric:
    def test_it_conforms_to_the_protocol_and_neither_names_a_trace(self):
        """``runtime_checkable`` compares NAMES only, so the parameter is gated by signature. First red: the Protocol's
        ``gauge(self, trace)`` before U3."""
        assert isinstance(_build(SPACE, ()), KernelGauge)
        for cls in (KernelGauge, KernelProjectionGauge):
            assert "trace" not in inspect.signature(cls.gauge).parameters, cls

    def test_the_axis_built_space_and_its_legacy_twin_give_bit_identical_gauges(self):
        """The metric is read through the space's metric leg. Arm ``weights_spelling`` (SN's predecessor read
        ``np.asarray(space.inner_product_weights, dtype=float)``, which is ``array(None, dtype=object)`` on an
        axis-built space) reddens it."""
        twin = FunctionSpace("twin", (_N,), inner_product_weights=_G)
        x = _x(14)
        for pairs in (_ROTATED, _OVERCOMPLETE):
            assert np.array_equal(_build(SPACE, pairs).gauge(x), _build(twin, pairs).gauge(x))

    def test_a_non_diagonal_metric_is_a_declared_scope_boundary(self):
        from orpheus.numerics.metric import DenseMetric
        dense = FunctionSpace("dense", (2,), metric=DenseMetric(np.array([[2.0, 0.5], [0.5, 1.0]])))
        with pytest.raises(ValueError, match="the metric is not diagonal"):
            KernelProjectionGauge.direct_sum(dense, ())

    def test_a_space_that_is_not_one_dimensional_is_refused(self):
        """A summand's support is flat indices. Arm ``no_rank_check`` (the 1-D check removed) reddens it."""
        plane = FunctionSpace.of_axes(
            Axis("p", (2,), weights=np.array([1.0, 2.0]), origin=Declared(BasisKind.NODAL)),
            Axis("q", (3,), weights=np.array([1.0, 1.0, 2.0]), origin=Declared(BasisKind.NODAL)))
        with pytest.raises(ValueError, match="must be one-dimensional"):
            KernelProjectionGauge.direct_sum(plane, ())

    def test_a_gauge_on_a_dual_space_is_admitted_and_is_a_projector(self):
        """The user's ruling (a), 2026-10-31: a dual space's metric G^+ is diagonal when G is, so the gauge is
        admitted (before the ruling it was refused as non-diagonal). Its laws hold in V*'s inner product, G^+ = 1/G.
        Arm ``dual_refused`` (the dual metric's diagonal is ``None``, the state before the ruling) reddens it."""
        dual = SPACE.dual()
        g = KernelProjectionGauge.direct_sum(dual, _summands(_ROTATED))
        ginv = 1.0 / _G
        x, y = _x(15), _x(16)
        np.testing.assert_allclose(g.apply(g.apply(x)), g.apply(x), rtol=0, atol=16 * _EPS * np.abs(x).max())
        assert float(np.sum(ginv * g.apply(x) * y)) == pytest.approx(float(np.sum(ginv * x * g.apply(y))), rel=1e-13)
        assert float(np.sum(_G * g.apply(x) * y)) != pytest.approx(float(np.sum(_G * x * g.apply(y))), rel=1e-6)

    def test_the_algebra_adjoint_is_the_projector_itself(self, rotated):
        """Pi is spelled in the operator algebra (the user's ruling (c)), so its adjoint follows by law: a G-orthogonal
        projector is G-self-adjoint, Pi.H == Pi. Arm ``transpose_is_apply`` reddens it (then Pi.H = G^-1 Pi G)."""
        x = _x(17)
        np.testing.assert_allclose(rotated.H.apply(x), rotated.apply(x), rtol=0, atol=64 * _EPS * np.abs(x).max())

    @pytest.mark.rests_on("tests/gates/numerics/test_coordinate_space.py::TestThePrimitive")
    @pytest.mark.parametrize("verb", ["apply", "apply_transpose", "gauge"])
    def test_it_parses_its_operand_at_its_domain(self, verb, rotated):
        """An operand whose LEADING shape is not the domain's is refused on every run, -O included (an action admits
        only a trailing payload, A ⊗ I_P). Arm ``no_parse`` (the verbs index the bare array) reddens all three: the
        fancy-indexed gather of an (N-1,) array either raises IndexError or returns."""
        with pytest.raises(ValueError, match="leading shape"):
            getattr(rotated, verb)(np.ones(_N - 1))


# ── the summand parses its fields ──────────────────────────────────────────────────────────────────────

class TestTheSummandParsesItsFields:
    """``KernelSummand.__post_init__``: a 1-D integer support, sorted, unique and non-negative; generators
    ``(support.size, k)``. Arm ``no_summand_parse`` (the parse made a no-op) reddens every row."""

    @pytest.mark.parametrize("support, fragment", [
        (np.array([0.0, 1.0]), "integer flat indices"),
        (np.array([[0, 1]]), "integer flat indices"),
        (np.array([1, 0]), "sorted, unique"),
        (np.array([0, 0]), "sorted, unique"),
        (np.array([-1, 0]), "sorted, unique"),
    ], ids=["float", "2-D", "unsorted", "duplicate", "negative"])
    def test_a_malformed_support_is_refused(self, support, fragment):
        with pytest.raises(ValueError, match=fragment):
            KernelSummand(support, np.eye(2), "bad")

    @pytest.mark.parametrize("generators", [np.ones(2), np.ones((3, 1))], ids=["1-D", "wrong rows"])
    def test_malformed_generators_are_refused(self, generators):
        with pytest.raises(ValueError, match="generators must be"):
            KernelSummand(np.array([0, 1]), generators, "bad")


# ── the metric's diagonal verb (the user's ruling (a), 2026-10-31) ────────────────────────────────────

class TestTheMetricDiagonalVerb:
    r"""``HilbertMetric.diagonal(shape)``: :math:`G\mathbf 1` iff the metric is diagonal, else ``None``. The oracle is
    the hand-written weights (never the SUT's ``apply``)."""

    def test_a_diagonal_metric_returns_its_weights(self):
        from orpheus.numerics.metric import DiagonalMetric
        np.testing.assert_array_equal(DiagonalMetric(_G).diagonal((_N,)), _G)

    def test_a_dense_metric_is_not_known_to_be_diagonal(self):
        """Arm ``dense_says_diagonal`` reddens it, and the scope-boundary row."""
        from orpheus.numerics.metric import DenseMetric
        assert DenseMetric(np.diag(_G)).diagonal((_N,)) is None

    def test_a_factored_metric_is_diagonal_iff_every_block_is(self):
        """Arm ``factored_ignores_blocks`` reddens the Dense-block leg."""
        from orpheus.numerics.metric import DenseMetric, DiagonalMetric, FactoredMetric
        w = np.array([1.0, 3.0])
        diag = FactoredMetric((((2,), DiagonalMetric(w)), ((3,), None)))
        np.testing.assert_array_equal(diag.diagonal((2, 3)), np.repeat(w[:, None], 3, axis=1))
        assert FactoredMetric((((2,), DenseMetric(np.eye(2))), ((3,), None))).diagonal((2, 3)) is None

    def test_a_dual_metric_is_diagonal_iff_its_primal_is(self):
        """Arm ``dual_is_primal`` reddens the diagonal leg."""
        from orpheus.numerics.metric import DenseMetric, DiagonalMetric, DualMetric
        np.testing.assert_allclose(DualMetric(DiagonalMetric(_G)).diagonal((_N,)), 1.0 / _G, rtol=1e-15, atol=0)
        assert DualMetric(DenseMetric(np.diag(_G))).diagonal((_N,)) is None

    def test_a_euclidean_space_reads_ones(self):
        flat = FunctionSpace("euclid", (4,))
        np.testing.assert_array_equal(flat.metric_diagonal(), np.ones(4))


class TestTheCompositeMetricDiagonal:
    """A composite's diagonal is its blocks' in the flat layout, or ``None`` if any block is not diagonal; the oracle is
    each block's own ``apply_metric(ones)``. Arms ``full_field_reversed`` and ``coupled_reversed`` (the blocks
    concatenated in the wrong order) redden their rows."""

    _A = FunctionSpace.of_axes(Axis("bulk", (3,), weights=np.array([1.0, 2.0, 3.0]), origin=Declared(BasisKind.NODAL)))
    _B = FunctionSpace.of_axes(Axis("trace", (2,), weights=np.array([0.5, 4.0]), origin=Declared(BasisKind.NODAL)))

    def _expected(self):
        return np.concatenate([self._A.apply_metric(np.ones(3)), self._B.apply_metric(np.ones(2))])

    def test_a_full_field_space_concatenates_interior_then_trace(self):
        from orpheus.numerics.spaces.full_field_space import FullFieldSpace
        np.testing.assert_array_equal(FullFieldSpace.from_blocks(self._A, self._B).metric_diagonal(), self._expected())

    def test_a_coupled_space_concatenates_its_systems_in_order(self):
        from orpheus.numerics.coupled_system import CoupledSpace
        np.testing.assert_array_equal(CoupledSpace.from_systems((self._A, self._B)).metric_diagonal(), self._expected())

    def test_a_composite_with_a_dense_block_is_not_diagonal(self):
        from orpheus.numerics.metric import DenseMetric
        from orpheus.numerics.coupled_system import CoupledSpace
        from orpheus.numerics.spaces.full_field_space import FullFieldSpace
        dense = FunctionSpace("dense", (2,), metric=DenseMetric(np.array([[2.0, 0.5], [0.5, 1.0]])))
        assert FullFieldSpace.from_blocks(self._A, dense).metric_diagonal() is None
        assert CoupledSpace.from_systems((self._A, dense)).metric_diagonal() is None

    @pytest.mark.xfail(strict=True, raises=AttributeError, reason=(
        "#666: a composite space's metric legs take its composite carrier, not its flat coordinates, so the summand's "
        "gather.H (sharp . dual . flat) cannot act on the flat chart the gauge acts on. `[M]` 2026-10-31: FullFieldSpace "
        "raises \"'numpy.ndarray' object has no attribute 'interior'\", CoupledSpace \"... 'n_systems'\". The fix lands in "
        "U0 (the first gauge posed on diffusion's full field, R3's union kernel). An XPASS means it landed: drop the mark."))
    @pytest.mark.parametrize("composite", ["full_field", "coupled"])
    def test_a_gauge_on_a_composite_space_spans_its_blocks(self, composite):
        """R3's union kernel lives on bulk (+) trace: a summand crossing the block boundary is admitted and is a
        G-orthogonal projector in the composite's metric."""
        from orpheus.numerics.coupled_system import CoupledSpace
        from orpheus.numerics.spaces.full_field_space import FullFieldSpace
        space = (FullFieldSpace.from_blocks(self._A, self._B) if composite == "full_field"
                 else CoupledSpace.from_systems((self._A, self._B)))
        g = KernelProjectionGauge.direct_sum(space, (KernelSummand(np.array([1, 2, 3]), np.array([[1.0], [-1.0], [2.0]]), "x"),))
        w = self._expected()
        x, y = np.random.default_rng(18).standard_normal((2, 5))
        np.testing.assert_allclose(g.apply(g.apply(x)), g.apply(x), rtol=0, atol=16 * _EPS * np.abs(x).max())
        assert float(np.sum(w * g.apply(x) * y)) == pytest.approx(float(np.sum(w * x * g.apply(y))), rel=1e-13)

def test_the_composite_xfail_is_strict():
    marks = [m for m in TestTheCompositeMetricDiagonal.test_a_gauge_on_a_composite_space_spans_its_blocks.pytestmark
             if m.name == "xfail"]
    assert len(marks) == 1 and marks[0].kwargs.get("strict") is True
