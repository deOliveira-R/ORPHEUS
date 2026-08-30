r"""Intrinsic-property tests for the frame faces + the discipline-type hierarchy.

The angular :class:`~orpheus.numerics.frame.GalerkinFrame` binds a :class:`Basis` to a
:class:`DiscreteMeasure` and exposes the analysis / reconstruction faces. These tests pin:

* the **adjoint-for-free** — BOTH ``frame.analysis.H`` and ``frame.reconstruction.H``
  fall out of the frame's swapped spaces with no bespoke code (each pinned against an
  INDEPENDENT closed-form einsum: :math:`M^* = S_0 \circ G^{-1} = R/W` for analysis,
  :math:`R^* = W\,M` for reconstruction — the F-0 Parseval metric);
* the **F-0 Parseval metric** (``frame_square_recarve.md``) — the codomain metric is
  the INVERSE discrete trial Gram, so analysis is an isometry onto its image, the
  frame square closes with one scalar, and the slab's non-diagonal Gram is refused
  VISIBLY (the ``TestParseval``-prefixed suite at the end of this module);
* the symmetric **space pairing** (basis → ``basis_space``; measure → ``measure_space``);
* the faces **compose through ``OperatorProduct`` with real spaces** (no ``cast``) —
  the enabler for the Phase-C cast retirement;
* the structural Galerkin invariant :math:`\Pi R = 4\pi I` routed through the frame.

The full SH-space law suite (:math:`\Pi R = 4\pi I`, :math:`\Pi^* = g_C S_0`,
:math:`R = (2\ell+1) S_0`) lives in ``test_spherical_harmonic_space.py``,
constructed on the same frame faces.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.basis import (
    GramStructure,
    IndicatorBasis,
    OverlapBasis,
    SphericalHarmonicBasis,
    WeightedIndicatorBasis,
)
from dataclasses import replace

from orpheus.numerics.frame import FrameBase, GalerkinFrame, PetrovGalerkinFrame
from orpheus.numerics.measure import DiscreteMeasure
from orpheus.numerics.metric import _DENSE_METRIC_RCOND, DenseMetric
from orpheus.numerics.operator import NotInvertible
from orpheus.numerics.quadrature import Quadrature, lebedev_sphere
from orpheus.numerics.spaces import SphericalHarmonicSpace
from tests._harness.predicates import (
    STRUCTURAL_ABSENT,
    assert_inverse_adjoint_contract,
)


@pytest.fixture
def sh_frame():
    """An exact spherical-harmonic frame: Lebedev(13) ⋈ SH(L=3)."""
    measure = lebedev_sphere(13)
    L = 3
    basis = SphericalHarmonicBasis(L=L)
    return GalerkinFrame(basis, measure), L


@pytest.mark.foundation
def test_frame_faces_two_axis_contract(sh_frame):
    r"""The carve keystone (#226, P4 keystone v2) for the numerics frame faces.

    Both faces carry a working ``apply_transpose`` — so ``.H`` falls out of
    the metric-aware ``AdjointOperator`` — hence ``is_adjointable`` is True
    (and the eager ``.H`` returns the wrapper); neither face is invertible,
    STRUCTURALLY (a projection face declares no ``inverse()`` — misuse is a
    static error, Design C). The reconstruction face's ``is_adjointable`` is
    the OVERRIDE that lifts it above the bare
    :class:`~orpheus.numerics.projection.ReconstructionOperator` default
    (``is_adjointable`` False); the analysis face inherits ``True`` from
    :class:`~orpheus.numerics.projection.AnalysisOperator`.
    """
    frame, _ = sh_frame
    for face in (frame.analysis, frame.reconstruction):
        assert_inverse_adjoint_contract(
            face,
            invertible=False,
            adjointable=True,
            inverse_contract=STRUCTURAL_ABSENT,
        )


def _band_limited(rng, L, *trailing):
    """A random ``(L+1, 2L+1, *trailing)`` moment array with |m|>ℓ slots zeroed."""
    c = rng.standard_normal((L + 1, 2 * L + 1, *trailing))
    for ell in range(L + 1):
        c[ell, 2 * ell + 1 :] = 0.0
    return c


# ── the adjoint-for-free ──────────────────────────────────────────────────

@pytest.mark.foundation
def test_analysis_hilbert_adjoint_falls_out_of_the_frame_spaces(sh_frame):
    r"""``frame.analysis.H`` is the PHYSICAL Hilbert adjoint :math:`M^* = S_0 \circ G^{-1}`.

    No bespoke adjoint code — the frame's swapped ``(measure_space, basis_space)``
    metrics feed the generic ``AdjointOperator``, and the F-0 Parseval dressing
    (the codomain metric is the INVERSE discrete Gram — see
    :attr:`FrameBase.basis_space`) makes the sandwich the physical adjoint. Pinned
    against an INDEPENDENT reference: the direct :math:`S_0(G^{-1}c)` einsum with
    the closed-form inverse SH Gram :math:`(2\ell+1)/4\pi` (exact for the
    degree-exact Lebedev rule; NOT the frame's own contraction). Equivalently
    :math:`M^* = R/W` — the closure pinned family-wide by
    ``test_parseval_frame_square_closes`` below.

    ⛔ Pre-F-0 this gate pinned :math:`g_C \cdot S_0` with the CONTINUUM Gram
    :math:`g_C = 4\pi/(2\ell+1)` as the codomain metric — the WRONG side for
    the carried covariant moments (`[M]` Parseval ratio 118.7 vs 1.000;
    ``scratch/probe_f1_parseval.py``, 2026-08-24). The loaded-not-blind negative
    leg is ``test_parseval_reds_under_the_pre_repair_continuum_metric``.
    """
    frame, L = sh_frame
    rng = np.random.default_rng(14)
    c = _band_limited(rng, L, 4, 2)
    g_inv = (2.0 * np.arange(L + 1) + 1.0) / (4.0 * np.pi)  # closed-form G⁻¹ diag
    expected = np.einsum("nlm,l,lm...->n...", frame.table, g_inv, c)
    np.testing.assert_allclose(
        frame.analysis.H.apply(c), expected, rtol=1e-12, atol=1e-14,
    )


@pytest.mark.foundation
def test_reconstruction_hilbert_adjoint_falls_out_of_the_frame_spaces(sh_frame):
    r"""``frame.reconstruction.H`` is the PHYSICAL Hilbert adjoint :math:`R^* = W\,M`.

    Symmetric with the analysis face: the F-0-dressed domain metric
    (:math:`G^{-1}`) enters the sandwich through its pseudo-inverse :math:`G`,
    giving :math:`(R^* v)_\ell^m = d_\ell\,G_\ell \sum_n w_n
    Y_\ell^m(\hat\Omega_n)\, v_n` — and the SH identity :math:`d_\ell G_\ell =
    (2\ell+1)\cdot 4\pi/(2\ell+1) = 4\pi = W` collapses the per-:math:`\ell`
    factor to the ONE scalar :math:`W`. Pinned against that INDEPENDENT
    closed-form einsum (NOT the frame's own contraction).
    ``R : \text{basis} \to \text{measure}``, so ``R.H`` maps nodal values →
    coefficients. (Pre-F-0 this pinned :math:`(2\ell+1)^2/4\pi\cdot Y^{\mathsf T}W`
    — the continuum-metric sandwich.)
    """
    frame, L = sh_frame
    rng = np.random.default_rng(17)
    n = frame.measure.weights.shape[0]
    v = rng.standard_normal((n, 4, 2))
    expected = 4.0 * np.pi * np.einsum(
        "n,nlm,n...->lm...", frame.measure.weights, frame.table, v,
    )
    np.testing.assert_allclose(
        frame.reconstruction.H.apply(v), expected, rtol=1e-12, atol=1e-14,
    )


# ── the space pairing ─────────────────────────────────────────────────────

@pytest.mark.foundation
def test_basis_space_is_the_spherical_harmonic_space(sh_frame):
    frame, L = sh_frame
    # basis.space rebuilds per call (cheap); the Frame caches it. Equal by
    # (name, shape), and the Frame's cached instance is shared across the faces.
    assert frame.basis_space == frame.basis.space
    assert frame.basis_space == SphericalHarmonicSpace.from_L(L)
    # the analysis codomain / reconstruction domain are that same space
    assert frame.analysis.codomain is frame.basis_space
    assert frame.reconstruction.domain is frame.basis_space


@pytest.mark.foundation
def test_measure_space_carries_the_quadrature_weights_as_its_metric(sh_frame):
    frame, _ = sh_frame
    ms = frame.measure_space
    assert ms.shape == (frame.measure.weights.shape[0],)
    np.testing.assert_array_equal(ms.inner_product_weights, frame.measure.weights)
    # analysis domain / reconstruction codomain are that same space
    assert frame.analysis.domain is ms
    assert frame.reconstruction.codomain is ms


# ── composition through OperatorProduct (the cast-retirement enabler) ─────

@pytest.mark.foundation
def test_faces_compose_through_operator_product_with_real_spaces(sh_frame):
    """``reconstruction @ analysis`` builds an ``OperatorProduct`` (no ``cast``).

    Both faces carry real ``domain``/``codomain``, so the composition's space-
    compatibility check is live (not skipped on ``None``) and passes — which is what
    let Phase C drop the ``cast(LinearOperator, …)`` workarounds in the scattering
    kernel (now retired).
    """
    frame, _ = sh_frame
    rng = np.random.default_rng(15)
    psi = rng.standard_normal((frame.measure.weights.shape[0], 4, 2))
    product = frame.reconstruction @ frame.analysis
    expected = frame.reconstruction.apply(frame.analysis.apply(psi))
    np.testing.assert_array_equal(product.apply(psi), expected)


@pytest.mark.foundation
def test_pi_R_is_4pi_identity_through_the_frame(sh_frame):
    """``analysis ∘ reconstruction = 4π·I`` on band-limited coefficients (via the frame).

    Inherited from the structural ``test_spherical_harmonic_space`` invariant by the
    faces' bit-identity; pinned here on the frame-routed path.
    """
    frame, L = sh_frame
    rng = np.random.default_rng(16)
    c = _band_limited(rng, L)
    out = frame.analysis.apply(frame.reconstruction.apply(c))
    np.testing.assert_allclose(out, 4.0 * np.pi * c, rtol=1e-10, atol=1e-12)


# ── caching + capability surface ──────────────────────────────────────────

@pytest.mark.foundation
def test_table_and_faces_are_cached(sh_frame):
    frame, _ = sh_frame
    assert frame.table is frame.table
    assert frame.analysis is frame.analysis
    assert frame.reconstruction is frame.reconstruction


# (The former ``test_face_capabilities`` API-smoke — caps == {apply,
# apply_transpose} on both faces — was retired with the frozenset at carve
# P4; its surviving claim, both faces adjointable-not-invertible, is pinned
# in FULL by ``test_frame_faces_two_axis_contract`` above.)


# ── the discipline-type hierarchy (P1) ─────────────────────────────────────

@pytest.mark.foundation
def test_galerkin_is_a_petrov_galerkin_is_a_frame_base(sh_frame):
    """``GalerkinFrame ⊂ PetrovGalerkinFrame ⊂ FrameBase`` — discipline is the TYPE.

    Liskov-correct: a Galerkin frame IS-A Petrov-Galerkin frame (with ``test is
    trial``), strengthening — never weakening — the base promises.
    """
    frame, _ = sh_frame
    assert isinstance(frame, GalerkinFrame)
    assert isinstance(frame, PetrovGalerkinFrame)
    assert isinstance(frame, FrameBase)


@pytest.mark.foundation
def test_galerkin_frame_test_is_the_trial_basis(sh_frame):
    """The Galerkin specialisation fixes ``test = trial`` and reuses the trial caches.

    Reusing :attr:`table`/:attr:`basis_space` (not re-evaluating) is what keeps the
    Galerkin analysis 0-ULP-identical to the single-discipline frame this hierarchy
    replaced, and preserves the analysis-codomain ``is`` identity.
    """
    frame, _ = sh_frame
    assert frame.test is frame.basis
    assert frame.test_table is frame.table
    assert frame.test_space is frame.basis_space
    assert frame.analysis.codomain is frame.basis_space


@pytest.mark.foundation
def test_galerkin_frame_takes_no_test_basis():
    """``GalerkinFrame`` binds test=trial; its constructor takes no ``test_basis``.

    The ``test ≠ trial`` freedom is exactly what a :class:`PetrovGalerkinFrame` is for,
    so a distinct test basis is forbidden on a :class:`GalerkinFrame` by the constructor
    SIGNATURE itself (a ``TypeError`` on the extra argument), not a runtime guard.
    """
    measure = lebedev_sphere(13)
    # ``*args`` so the arity violation is a runtime TypeError, not a static type error.
    args = [SphericalHarmonicBasis(L=3), measure, SphericalHarmonicBasis(L=2)]
    with pytest.raises(TypeError):
        GalerkinFrame(*args)


@pytest.mark.foundation
def test_petrov_galerkin_degenerate_equals_galerkin_bit_identically(sh_frame):
    """A ``PetrovGalerkinFrame`` with ``test_basis = trial`` is the Galerkin degenerate.

    Passing the trial basis itself as the test basis resolves test→trial, so the GENERAL
    Petrov-Galerkin analysis (which reads the TEST table) must reduce BIT-IDENTICALLY to
    the Galerkin analysis. This pins the PG analysis machinery in the degenerate case
    here; the genuine ``test ≠ trial`` instance (flux-weighted homogenisation) lands
    with its consumer in a later phase.
    """
    galerkin, L = sh_frame
    pg = PetrovGalerkinFrame(galerkin.basis, galerkin.measure, galerkin.basis)
    assert pg.test is pg.basis
    rng = np.random.default_rng(24)
    psi = rng.standard_normal((galerkin.measure.weights.shape[0], 4, 2))
    np.testing.assert_array_equal(
        pg.analysis.apply(psi), galerkin.analysis.apply(psi),
    )
    c = _band_limited(rng, L)
    np.testing.assert_array_equal(
        pg.reconstruction.apply(c), galerkin.reconstruction.apply(c),
    )


# ── the project verb: G⁻¹ M (homogenise / condense) — P3 ────────────────────

def _indicator_frame(edges, centres, weights, test_weight=None):
    """A small hand-built indicator frame on an explicit measure.

    ``test_weight=None`` → :class:`GalerkinFrame` (test=trial=plain indicator);
    an array → :class:`PetrovGalerkinFrame` with a flux-weighted indicator test.
    """
    trial = IndicatorBasis((np.asarray(edges, dtype=float),))
    measure = DiscreteMeasure(
        nodes=np.asarray(centres, dtype=float),
        weights=np.asarray(weights, dtype=float),
        support="spatial_R1",
    )
    if test_weight is None:
        return GalerkinFrame(trial, measure)
    return PetrovGalerkinFrame(
        trial, measure, WeightedIndicatorBasis(trial, np.asarray(test_weight, float)),
    )


@pytest.mark.foundation
def test_project_is_gram_inverse_times_analysis():
    r"""``frame.project(f) = G⁻¹ M f`` on a hand-known diagonal-Gram frame.

    Galerkin indicator frame (test=trial), 4 fine nodes / non-uniform volumes, 3
    coarse cells the LAST of which is EMPTY (no fine node). The diagonal Gram is the
    region volume :math:`m_R = \sum_{i\in R} V_i`; ``project`` is the volume-weighted
    average :math:`(\sum_R V_i f_i)/m_R`. The empty region (``m_R = 0``) must project
    to ``0`` (the Moore–Penrose pseudo-inverse), NOT ``nan``/``inf`` — the verb-level
    pin of the zero-flux-region law.
    """
    centres = [0.5, 1.5, 2.5, 3.5]
    V = [1.0, 1.0, 2.0, 1.0]
    f = np.array([10.0, 20.0, 30.0, 40.0])
    frame = _indicator_frame([0.0, 2.0, 4.0, 5.0], centres, V)  # R2=[4,5] empty
    out = frame.project(f)
    expected = np.array([
        (1.0 * 10 + 1.0 * 20) / (1.0 + 1.0),    # R0: nodes 0,1
        (2.0 * 30 + 1.0 * 40) / (2.0 + 1.0),    # R1: nodes 2,3
        0.0,                                     # R2: empty → 0 (pseudo-inverse)
    ])
    np.testing.assert_allclose(out, expected, rtol=1e-14)
    assert np.isfinite(out).all(), "empty region produced nan/inf, not 0"


@pytest.mark.foundation
def test_petrov_galerkin_project_is_cross_gram_extraction():
    r"""``PetrovGalerkinFrame.project`` extracts coefficients against the CROSS Gram.

    A genuine ``test ≠ trial`` frame (test = a flux-weighted indicator ``w·1_R``,
    trial = ``1_R``): the diagonal cross Gram is :math:`G_R = \langle\chi_R,
    \mathbf 1_R\rangle_W = \sum_{i\in R} w_i V_i` and :math:`(M f)_R = \sum_{i\in R}
    w_i V_i f_i`, so ``project = M f / G``. Pinned against the independent hand
    arithmetic (NOT a re-call of the production einsum).
    """
    centres = [0.5, 1.5, 2.5, 3.5]
    V = np.array([0.5, 1.5, 1.0, 1.0])
    w = np.array([2.0, 3.0, 5.0, 7.0])
    f = np.array([10.0, 20.0, 30.0, 40.0])
    frame = _indicator_frame([0.0, 2.0, 4.0], centres, V, test_weight=w)  # 2 cells
    regions = [[0, 1], [2, 3]]
    expected = np.array([
        sum(w[i] * V[i] * f[i] for i in sel) / sum(w[i] * V[i] for i in sel)
        for sel in regions
    ])
    np.testing.assert_allclose(frame.project(f), expected, rtol=1e-13)


@pytest.mark.foundation
def test_petrov_galerkin_degenerate_project_equals_galerkin_project(sh_frame):
    r"""``PetrovGalerkinFrame(b, m, test=b).project ≡ GalerkinFrame(b, m).project``.

    The ``project``-verb analogue of the face-level degenerate test: when ``test is
    trial`` the general PG ``project`` (which reads the TEST Gram) must reduce to the
    SAME numpy chain BIT-IDENTICALLY (``array_equal``, the 0-ULP discipline).
    """
    galerkin, _ = sh_frame
    pg = PetrovGalerkinFrame(galerkin.basis, galerkin.measure, galerkin.basis)
    rng = np.random.default_rng(31)
    psi = rng.standard_normal((galerkin.measure.weights.shape[0], 4, 2))
    np.testing.assert_array_equal(pg.project(psi), galerkin.project(psi))


@pytest.mark.foundation
def test_petrov_galerkin_project_differs_from_galerkin_when_test_neq_trial():
    r"""The PG type is LOAD-BEARING: ``test ≠ trial`` gives a DIFFERENT answer.

    The same geometry projected (a) flux-weighted (the PG test ``w·1_R``) and (b)
    plain Galerkin (test=trial=``1_R``, the volume average) gives materially distinct
    coefficients — the type carries real information, it is not ceremony. Both match
    their respective independent hand references; the discrimination is asserted to
    have actually fired (no silent same-answer pass).
    """
    centres = [0.5, 1.5, 2.5, 3.5]
    V = np.array([0.5, 1.5, 1.0, 1.0])
    w = np.array([2.0, 3.0, 5.0, 7.0])
    f = np.array([10.0, 20.0, 30.0, 40.0])
    edges = [0.0, 2.0, 4.0]
    regions = [[0, 1], [2, 3]]

    pg = _indicator_frame(edges, centres, V, test_weight=w).project(f)
    galerkin = _indicator_frame(edges, centres, V).project(f)

    pg_ref = np.array([
        sum(w[i] * V[i] * f[i] for i in sel) / sum(w[i] * V[i] for i in sel)
        for sel in regions
    ])
    gal_ref = np.array([
        sum(V[i] * f[i] for i in sel) / sum(V[i] for i in sel) for sel in regions
    ])
    np.testing.assert_allclose(pg, pg_ref, rtol=1e-13)
    np.testing.assert_allclose(galerkin, gal_ref, rtol=1e-13)
    assert not np.allclose(pg, galerkin, rtol=1e-6), (
        "PG (flux-weighted) and Galerkin (volume) projections coincided — the test "
        "weight does not discriminate; the PG type would be ceremony here"
    )


# ── Gram-structure: the projection-validity declaration (P5.5a) ───────────


@pytest.mark.foundation
def test_basis_gram_structure_declarations():
    r"""Each trial basis declares the Gram structure that makes ``project`` valid.

    Intrinsic-property pin: ``IndicatorBasis`` (disjoint cells) and
    ``SphericalHarmonicBasis`` (orthogonal) are DIAGONAL; ``OverlapBasis`` (fractional
    rows summing to 1) is PARTITION_OF_UNITY — it MUST override the DIAGONAL it inherits
    from ``IndicatorBasis`` (a straddling row shares ≥2 columns ⟹ non-diagonal Gram).
    The base ``Basis`` default (here via the test-only ``WeightedIndicatorBasis``) is
    DENSE — the safe refusal a new basis inherits until it consciously declares.
    """
    ib = IndicatorBasis((np.array([0.0, 1.0, 2.0]),))
    ob = OverlapBasis(
        edges_per_axis=(np.array([-0.5, 0.5, 1.5]),),
        overlap_table=np.array([[1.0, 0.0], [0.5, 0.5], [0.0, 1.0]]),
    )
    assert ib.gram_structure is GramStructure.DIAGONAL
    assert SphericalHarmonicBasis(L=1).gram_structure is GramStructure.DIAGONAL
    assert ob.gram_structure is GramStructure.PARTITION_OF_UNITY
    # OverlapBasis IS-A IndicatorBasis but MUST NOT inherit its DIAGONAL claim.
    assert ob.gram_structure is not ib.gram_structure
    # The base default (a test-only basis never used as a trial) is the safe DENSE.
    assert WeightedIndicatorBasis(ib, np.ones(2)).gram_structure is GramStructure.DENSE


@pytest.mark.foundation
def test_project_refuses_dense_gram_trial():
    r"""``project`` / ``.gram`` REFUSE a DENSE-Gram trial — illegal state unrepresentable.

    The row-sum probe is wrong for a trial that is neither disjoint nor a partition of
    unity; rather than return a silently-wrong coarsening, the frame raises
    :class:`NotInvertible` (the dense ``(MR)⁻¹M`` solve is unbuilt — #275). Mutation
    gate: a trial declaring ``GramStructure.DENSE`` reddens BOTH ``.gram`` and
    ``.project``, while the otherwise-identical DIAGONAL trial succeeds — proving the
    refusal keys on the declaration, not on some unrelated failure.
    """

    class _DenseTrial(IndicatorBasis):
        @property
        def gram_structure(self) -> GramStructure:
            return GramStructure.DENSE

    edges = np.array([0.0, 1.0, 2.0])
    measure = DiscreteMeasure(
        nodes=np.array([0.5, 1.5]), weights=np.ones(2), support="spatial_R1",
    )
    dense = PetrovGalerkinFrame(
        _DenseTrial((edges,)), measure, WeightedIndicatorBasis(IndicatorBasis((edges,)), np.ones(2)),
    )
    with pytest.raises(NotInvertible, match="DENSE"):
        _ = dense.gram
    with pytest.raises(NotInvertible, match="DENSE"):
        dense.project(np.array([3.0, 5.0]))
    # Control: the SAME geometry with the honest DIAGONAL trial projects fine.
    ok = PetrovGalerkinFrame(
        IndicatorBasis((edges,)), measure,
        WeightedIndicatorBasis(IndicatorBasis((edges,)), np.ones(2)),
    )
    np.testing.assert_allclose(ok.project(np.array([3.0, 5.0])), [3.0, 5.0])


# ── the F-0 Parseval metric — the metric truth (frame_square_recarve F-0) ──
#
# THE THEOREM (exact, unconditional — algebra, not quadrature exactness): for a
# band-limited field ψ = S₀c the analysis output is φ = Mψ = Gc IDENTICALLY,
# with G the discrete TRIAL Gram of the pairing (basis ⊗ measure). So the
# codomain inner product under which analysis is an isometry onto its image is
# the INVERSE discrete Gram:
#
#     ‖φ‖²_{G⁻¹} = cᵀG c = ‖ψ‖²_W                     (Parseval)
#
# and with that metric both faces' .H are the PHYSICAL Hilbert adjoints. Two
# consequences, split by precondition:
#
#   * Parseval needs only a DIAGONAL discrete Gram — any values (LS4 at L=2 is
#     dressed with its true discrete inverse and closes exactly, whatever its
#     relation to the continuum Gram);
#   * the SH closure M* = R/W, R* = W·M additionally needs the per-ℓ identity
#     d_ℓ·G_ℓ = W (degree-exactness; [M] 2026-08-24 every shipped sphere
#     family measures exact to ~1e-15, incl. LS4/LS8 at L=2).
#
# [M] scratch/probe_f1_parseval.py + probe_f1_parseval_slab.py (2026-08-24):
# the pre-F-0 stored metric (continuum 4π/(2ℓ+1)) was the WRONG side — ratio
# 118.7 vs 1.000 on that probe's seed (the ratio is a moment-energy-weighted
# average of the per-ℓ factors, so it is draw-dependent; the draw-independent
# statement is the (4π/(2ℓ+1))² per-ℓ factor); closure ≤1e-15 once dressed;
# the slab GL live Gram has off-diagonals at 0.93 of the Cauchy–Schwarz
# scale, so NO diagonal Parseval metric exists there (the DENSE arm below).

def _overlap_frame() -> GalerkinFrame:
    """The PoU overlap frame — measured-DENSE with an INVERTIBLE Gram
    ([M] ``[[1.25, .25], [.25, 1.25]]``, cond 1.50): the non-angular
    member of the DENSE population."""
    ob = OverlapBasis(
        edges_per_axis=(np.array([-0.5, 0.5, 1.5]),),
        overlap_table=np.array([[1.0, 0.0], [0.5, 0.5], [0.0, 1.0]]),
    )
    measure = DiscreteMeasure(
        nodes=np.array([0.0, 0.5, 1.0]), weights=np.ones(3), support="spatial_R1",
    )
    return GalerkinFrame(ob, measure)


#: The six sphere families whose discrete Gram measures DIAGONAL — the
#: population of the DIAGONAL-arm dressing gate and of the scalar
#: frame-square collapse (a sphere-family property; see D3).
_DIAGONAL_FRAME_CASES = [
    pytest.param(lambda: Quadrature.level_symmetric(4).angular_frame(1), id="LS4-L1"),
    pytest.param(lambda: Quadrature.level_symmetric(4).angular_frame(2), id="LS4-L2"),
    pytest.param(lambda: Quadrature.level_symmetric(8).angular_frame(2), id="LS8-L2"),
    pytest.param(lambda: Quadrature.product(8, 8).angular_frame(2), id="product8x8-L2"),
    pytest.param(
        lambda: Quadrature.folded_product(8, 8).angular_frame(2), id="folded8x8-L2",
    ),
    pytest.param(
        lambda: GalerkinFrame(SphericalHarmonicBasis(L=2), lebedev_sphere(13)),
        id="lebedev13-L2",
    ),
]

#: Every Parseval-capable frame — the diagonal six PLUS the slab, whose
#: dressing is the matrix pseudo-inverse (P7). Until P7 the slab param
#: carried a skip mark (*"NON-DIAGONAL discrete Gram … no diagonal
#: Parseval metric exists; the matrix-metric home is the CS4c Riesz-leg
#: machinery"* — [M] 2026-08-23): the isometry now RUNS there through
#: the DenseMetric dressing, and the mark retired with the refusal.
_PARSEVAL_FRAME_CASES = _DIAGONAL_FRAME_CASES + [
    pytest.param(
        lambda: Quadrature.gauss_legendre(8).angular_frame(2), id="slab-GL8-L2",
    ),
]

#: The measured-DENSE population (D1) — four MECHANISMS, never one
#: family (vv-principles #13): a slab measure, a coarse product at L=2,
#: a coarse level-symmetric at L=3, and a non-angular partition-of-unity
#: basis.
_DENSE_FRAME_CASES = [
    pytest.param(
        lambda: Quadrature.gauss_legendre(8).angular_frame(2), id="slab-GL8-L2",
    ),
    pytest.param(lambda: Quadrature.product(4, 4).angular_frame(2), id="product4x4-L2"),
    pytest.param(lambda: Quadrature.level_symmetric(4).angular_frame(3), id="LS4-L3"),
    pytest.param(_overlap_frame, id="overlap-R1"),
]


@pytest.mark.foundation
@pytest.mark.parametrize("make_frame", _DIAGONAL_FRAME_CASES)
def test_parseval_dressing_installed_on_diagonal_frames(make_frame):
    r"""The verdict is DIAGONAL and ``basis_space`` carries the inverse discrete Gram.

    The dressed metric equals ``1/G_kk`` on live slots and EXACTLY ``0.0`` on
    dead ones (layout padding; the folded frame's σ-odd columns) — the dead-slot
    zeros are what make the Moore–Penrose inverse-metric path exact. The
    Galerkin override keeps the analysis codomain the SAME dressed object.
    The DENSE population's sibling is
    ``test_dense_frames_are_dressed_with_the_pseudo_inverse_gram`` (D1).
    """
    frame = make_frame()
    assert frame.discrete_gram_structure is GramStructure.DIAGONAL
    diag = np.diagonal(frame.discrete_gram).reshape(frame.basis.space.shape)
    live = diag > 0.0
    metric = frame.basis_space.inner_product_weights
    assert metric is not None
    np.testing.assert_allclose(metric[live], 1.0 / diag[live], rtol=1e-15)
    np.testing.assert_array_equal(metric[~live], 0.0)
    assert frame.test_space is frame.basis_space


@pytest.mark.l1
@pytest.mark.catches("ERR-039")
@pytest.mark.parametrize("make_frame", _PARSEVAL_FRAME_CASES)
def test_parseval_analysis_is_an_isometry_onto_its_image(make_frame):
    r"""``‖Mψ‖_{basis_space} = ‖ψ‖_W`` for band-limited ψ — Parseval, rtol 1e-12.

    The coefficient draw is deliberately UNMASKED (garbage in dead slots): a
    dead table column annihilates its coefficient in ψ = S₀c AND zeroes both
    its moment and its metric slot, so the identity must hold regardless —
    pinning the Moore–Penrose dead-slot handling for free.
    """
    frame = make_frame()
    rng = np.random.default_rng(1234)
    c = rng.standard_normal(frame.basis_space.shape)
    psi = frame.basis.synthesize(c, frame.table)
    phi = frame.analysis.apply(psi)
    np.testing.assert_allclose(
        frame.basis_space.inner_product(phi, phi),
        frame.measure_space.inner_product(psi, psi),
        rtol=1e-12,
    )


@pytest.mark.l1
@pytest.mark.catches("ERR-039")
@pytest.mark.verifies("hilbert-adjoint-equals-metric-times-S0")
@pytest.mark.parametrize("make_frame", _DIAGONAL_FRAME_CASES)
def test_parseval_frame_square_closes(make_frame):
    r"""``M.H = R/W`` and ``R.H = W·M`` — the frame square closes with ONE scalar.

    The SH-specific collapse of the general adjoints (:math:`M^* = S_0\circ
    G^{-1}`, :math:`R^* = d\,G\cdot(Y^{\mathsf T}W)`): per ℓ,
    :math:`d_\ell G_\ell = (2\ell+1)\cdot 4\pi/(2\ell+1) = 4\pi = W`, so the
    whole per-ℓ dressing collapses to the single scalar :math:`W = \sum_n w_n`
    — which IS the shipped scattering kernel's :math:`1/W` prefactor.
    `[M]` closure 5.6e-17 on the probe; every shipped sphere family measures
    degree-exact to ~1e-15.

    ⛔ The slab param is deliberately ABSENT, and D3
    (``test_diagonal_gram_suffices_for_the_collapse_and_dense_does_not_decide_it``)
    states why as a claim rather than a silent removal: `[M]` 2026-08-30,
    under the CORRECT dense Parseval metric the slab's collapse still
    fails (rel 2.65 on this file's seed) — its live ℓ=2 Gram diagonal
    ``[0.4, 0.8, 0.8]`` is not a per-ℓ scalar, so no :math:`G_\ell`
    exists there at ANY metric. ⚠ The honest quantifier (archivist,
    200-seed census): DIAGONAL is SUFFICIENT for the collapse — this
    gate's population — while DENSE does not decide it either way
    (``folded_product(4,6)`` L=3 is DENSE and satisfies it; the sphere
    rule ``product(4,4)`` L=2 breaks it).
    """
    frame = make_frame()
    W = float(frame.measure.weights.sum())
    rng = np.random.default_rng(4321)
    y = rng.standard_normal(frame.basis_space.shape)
    v = rng.standard_normal(frame.measure.weights.shape)
    np.testing.assert_allclose(
        frame.analysis.H.apply(y), frame.reconstruction.apply(y) / W,
        rtol=1e-12, atol=1e-13,
    )
    np.testing.assert_allclose(
        frame.reconstruction.H.apply(v), W * frame.analysis.apply(v),
        rtol=1e-12, atol=1e-13,
    )


@pytest.mark.foundation
def test_parseval_reds_under_the_pre_repair_continuum_metric():
    r"""The §6c witness: the Parseval isometry gate is LOADED on the metric, not blind.

    In-process pre-repair mutation (process-discipline: monkeypatch, never a git
    checkout): pre-seed the frame's cached-property slots with the UNDRESSED
    continuum-metric space (``SphericalHarmonicSpace.from_L`` — exactly what the
    pre-F-0 ``basis_space`` returned), and Parseval FAILS by the measured
    margin: the ratio is a weighted average of the per-ℓ factors
    :math:`(4\pi/(2\ell+1))^2 \ge 17.5` at L=1 (`[M]` 118.7 on the probe's
    seed). A BLIND gate would read 1.0 here (vv-principles #19: only the
    wrong-structure reading discriminates loaded from blind).
    """
    frame = Quadrature.level_symmetric(4).angular_frame(1)
    undressed = SphericalHarmonicSpace.from_L(1)
    # cached_property stores in the instance __dict__ — pre-seeding it IS the
    # pre-repair frame; no production code is touched and nothing needs undoing.
    frame.__dict__["basis_space"] = undressed
    frame.__dict__["test_space"] = undressed
    rng = np.random.default_rng(1234)
    c = rng.standard_normal(undressed.shape)
    psi = frame.basis.synthesize(c, frame.table)
    phi = frame.analysis.apply(psi)
    ratio = frame.basis_space.inner_product(phi, phi) / (
        frame.measure_space.inner_product(psi, psi)
    )
    assert ratio > 10.0, (
        f"pre-repair continuum metric read Parseval ratio {ratio:.3g} ≈ 1 — "
        f"the isometry gate would be BLIND to the wrong-side metric"
    )


@pytest.mark.foundation
def test_slab_frame_is_dressed_with_the_matrix_parseval_metric():
    r"""The slab GL frame measures DENSE and carries the MATRIX Parseval metric.

    The declared/measured separation stays (two questions, two
    properties): the SH basis DECLARES DIAGONAL (continuum-orthogonal);
    the measured verdict on THIS measure is DENSE (`[M]` 2026-08-23,
    discovery record ``scratch/probe_f1_parseval_slab.py``: total weight
    2 not 4π, live slots [1, 1, 3] per degree, live off-diagonals at
    0.93 of the Cauchy–Schwarz scale — NO diagonal candidate satisfies
    Parseval). Until P7 (2026-08-30) that verdict REFUSED the dressing —
    this gate then pinned the refusal (``basis_space`` array-equal to the
    undressed continuum metric), with the matrix home recorded as the
    CS4c Riesz-leg debt. P7 installed it: the space carries a
    :class:`~orpheus.numerics.metric.DenseMetric` whose matrix is the
    Moore–Penrose pseudo-inverse of the (symmetrized) measured Gram, the
    legacy weights slot reads ``None`` (it describes the DIAGONAL source
    only — a dense-metric space is NOT Euclidean), and the isometry gate
    runs on this frame (the retired skip). This flip IS the phase's §6c
    witness: the pre-P7 body was designed-red the moment the dressing
    landed.
    """
    frame = Quadrature.gauss_legendre(8).angular_frame(2)
    assert frame.basis.gram_structure is GramStructure.DIAGONAL       # declared
    assert frame.discrete_gram_structure is GramStructure.DENSE       # measured
    metric = frame.basis_space.metric
    assert isinstance(metric, DenseMetric)
    g = frame.discrete_gram
    np.testing.assert_allclose(
        metric.matrix,
        np.linalg.pinv(
            (g + g.T) / 2.0, hermitian=True, rcond=_DENSE_METRIC_RCOND
        ),
        rtol=1e-12, atol=1e-15,
    )
    assert frame.basis_space.inner_product_weights is None


@pytest.mark.foundation
def test_indicator_frame_parseval_metric_is_the_inverse_region_mass():
    r"""The SAME theorem on the indicator frame: the Parseval metric is ``1/m_R``.

    :math:`G_{RR} = \sum_{i\in R} w_i = m_R` (the region mass), so the dressed
    metric is ``1/m_R`` on occupied regions and EXACTLY ``0.0`` on the empty one
    (the dead-slot arm — matching ``project``'s Moore–Penrose convention).
    Parseval on a region-wise-constant (band-limited) field:
    :math:`\|Mf\|^2_{1/m} = \sum_R m_R \bar f_R^2 = \|f\|^2_V` exactly.
    """
    frame = _indicator_frame(
        [0.0, 2.0, 4.0, 5.0], [0.5, 1.5, 2.5, 3.5], [1.0, 1.0, 2.0, 1.0],
    )  # R2 = [4, 5] is empty
    assert frame.discrete_gram_structure is GramStructure.DIAGONAL
    np.testing.assert_allclose(
        frame.basis_space.inner_product_weights, [0.5, 1.0 / 3.0, 0.0],
    )
    f = np.array([10.0, 10.0, 30.0, 30.0])  # constant per region — band-limited
    phi = frame.analysis.apply(f)
    np.testing.assert_allclose(
        frame.basis_space.inner_product(phi, phi),
        frame.measure_space.inner_product(f, f),
        rtol=1e-13,
    )


@pytest.mark.foundation
def test_overlap_frame_measures_dense_while_declaring_partition_of_unity():
    r"""The declared/measured Gram facts are INDEPENDENT — the overlap witness.

    :class:`OverlapBasis` DECLARES PARTITION_OF_UNITY (the cross-Gram row-sum
    probe is valid for ``project`` — :math:`R\mathbf 1 = \mathbf 1`), while its
    TRIAL Gram MEASURES DENSE (a straddling row makes two columns share
    support). Since P7 the DENSE verdict means the Parseval dressing is
    INSTALLED (a matrix metric on ``basis_space``) while ``project``
    keeps working through the row-sum probe — which never inherits the
    dressing (the stripped ``gram``, pinned by the C3 gate below).
    """
    frame = _overlap_frame()
    assert frame.basis.gram_structure is GramStructure.PARTITION_OF_UNITY  # declared
    assert frame.discrete_gram_structure is GramStructure.DENSE    # measured
    assert isinstance(frame.basis_space.metric, DenseMetric)       # dressed (P7)
    assert frame.basis_space.inner_product_weights is None         # ≠ Euclidean
    np.testing.assert_allclose(
        frame.project(np.array([2.0, 4.0, 6.0])), [8.0 / 3.0, 16.0 / 3.0],
    )


@pytest.mark.foundation
def test_the_gram_row_sum_probe_survives_a_dense_dressed_test_space():
    r"""C3 (P7 S2, battery arm M14): ``gram``/``project`` are CROSS-Gram
    machinery and must never inherit the test space's Parseval dressing.

    The pre-P7 spelling ``replace(self.test_space,
    inner_product_weights=diagonal)`` carried a dense-dressed test
    space's metric OBJECT into the probe — [M] 2026-08-30 (pre-flight,
    ``scratch/p7/preflight.log``): ``frame.project([2,4,6])`` read
    ``[7.0, 11.0]`` against the true ``[8/3, 16/3]`` (rel 1.625), a
    silent VALUE error, with no guard involved. The repaired spelling
    strips the object while installing the row-sum diagonal.

    The dressed state is simulated by pre-seeding the ``basis_space``
    cache (the same idiom as the pre-repair-metric red gate) — so this
    witness has teeth NOW, before the S3 installer dresses anything.
    """
    from dataclasses import replace as _replace

    from orpheus.numerics.metric import DenseMetric

    ob = OverlapBasis(
        edges_per_axis=(np.array([-0.5, 0.5, 1.5]),),
        overlap_table=np.array([[1.0, 0.0], [0.5, 0.5], [0.0, 1.0]]),
    )
    measure = DiscreteMeasure(
        nodes=np.array([0.0, 0.5, 1.0]), weights=np.ones(3), support="spatial_R1",
    )
    frame = GalerkinFrame(ob, measure)
    dressed = _replace(
        frame.basis.space,
        metric=DenseMetric.inverse_of(frame.discrete_gram),
    )
    vars(frame)["basis_space"] = dressed  # the cached_property pre-seed idiom
    assert frame.test_space is dressed  # the Galerkin identity holds on the seed
    probe = frame.gram
    assert probe.metric is None, "the probe must not carry the dressing"
    np.testing.assert_array_equal(
        frame.project(np.array([2.0, 4.0, 6.0])), [8.0 / 3.0, 16.0 / 3.0],
    )


@pytest.mark.foundation
@pytest.mark.parametrize("make_frame", _DENSE_FRAME_CASES)
def test_dense_frames_are_dressed_with_the_pseudo_inverse_gram(make_frame):
    r"""D1 — the DENSE arm's counterpart of the diagonal dressing gate.

    Four mechanisms, one law: verdict DENSE ⟹ ``basis_space`` carries a
    :class:`DenseMetric` whose matrix is the Moore–Penrose pseudo-inverse
    of the measured (symmetrized) Gram at the module's pinned ``rcond``,
    and the Galerkin identity ``test_space is basis_space`` survives the
    dressing.
    """
    frame = make_frame()
    assert frame.discrete_gram_structure is GramStructure.DENSE
    metric = frame.basis_space.metric
    assert isinstance(metric, DenseMetric)
    g = frame.discrete_gram
    expected = np.linalg.pinv(
        (g + g.T) / 2.0, hermitian=True, rcond=_DENSE_METRIC_RCOND
    )
    np.testing.assert_allclose(metric.matrix, expected, rtol=1e-12, atol=1e-15)
    assert frame.test_space is frame.basis_space


@pytest.mark.l1
@pytest.mark.catches("ERR-039")
def test_no_diagonal_metric_can_satisfy_parseval_on_a_dense_frame():
    r"""D2 — THE WRONG-METRIC DISCRIMINATOR (vv-principles #19's loaded reading).

    Three readings of the SAME band-limited ψ on the slab GL8/L=2 frame
    (`[M]` 2026-08-30, THIS seed — the floors are draw-dependent, the
    dense ≈1 is a theorem): undressed continuum **25.53**, the best
    diagonal candidate ``1/diag(G)`` **1.806**, the dense pseudo-inverse
    **0.999999999999999**. The middle reading is the phase's whole
    justification: on this frame a diagonal metric is not merely
    unavailable, it is PROVABLY insufficient. And this family is the only
    correctness evidence the metric has — reciprocity holds to 1e-16 for
    EVERY invertible G (#409) and can never adjudicate one. ⚠ Slab-only
    by measurement: `[M]` the same diagonal candidate reads 1.066 on
    product(4,4) and 0.996 on LS4-L3 — the separation is frame-dependent,
    so only the slab's is gated.
    """
    frame = Quadrature.gauss_legendre(8).angular_frame(2)
    rng = np.random.default_rng(1234)
    c = rng.standard_normal(frame.basis.space.shape)
    psi = frame.basis.synthesize(c, frame.table)
    norm_w = frame.measure_space.inner_product(psi, psi)
    phi = frame.analysis.apply(psi)

    dense = frame.basis_space.inner_product(phi, phi) / norm_w
    diag = np.diagonal(frame.discrete_gram).reshape(frame.basis.space.shape)
    live = diag > 0.0
    diagonal_candidate = replace(
        frame.basis.space,
        inner_product_weights=np.where(
            live, 1.0 / np.where(live, diag, 1.0), 0.0
        ),
    )
    diagonal = diagonal_candidate.inner_product(phi, phi) / norm_w
    continuum = frame.basis.space.inner_product(phi, phi) / norm_w
    assert dense == pytest.approx(1.0, rel=1e-12)
    assert diagonal > 1.5, f"diagonal candidate read {diagonal:.4f}"
    assert continuum > 10.0, f"continuum metric read {continuum:.4f}"


@pytest.mark.foundation
def test_diagonal_gram_suffices_for_the_collapse_and_dense_does_not_decide_it():
    r"""D3 — why the frame-square gate keeps NO slab param.

    The positive replacement for a silent param removal: `[M]` 2026-08-30
    under the CORRECT dense Parseval metric, the isometry holds (≈1,
    rtol 1e-12) while ``M.H`` vs ``R/W`` reads rel **2.65** (this seed) —
    because the slab's live ℓ=2 Gram diagonal is ``[0.4, 0.8, 0.8]``,
    not a single per-ℓ scalar, so no :math:`G_\ell` exists and the
    collapse :math:`d_\ell G_\ell = W` is unspellable there at ANY
    metric.

    ⛔ REFRAMED 2026-08-30 (archivist refutation of this gate's first
    name, "…is a sphere-family property"): `[M]` ``product(4,4)`` L=2 IS
    a sphere rule and BREAKS the collapse (rel 3.1e-3–0.333 over 200
    seeds), while ``folded_product(4,6)`` L=3 measures DENSE and
    SATISFIES it to 2.8e-15 — its only live off-diagonal block is
    rank-1 (det −8.7e-17), so :math:`Y(G^{+} − \mathrm{diag}(d)/W) = 0`
    holds anyway. The decidable statement is the current name:
    a DIAGONAL verdict is SUFFICIENT for the collapse; a DENSE verdict
    does not decide it either way. The slab is a DENSE member where it
    fails, which is all this gate pins.
    """
    frame = Quadrature.gauss_legendre(8).angular_frame(2)
    rng = np.random.default_rng(1234)
    c = rng.standard_normal(frame.basis.space.shape)
    psi = frame.basis.synthesize(c, frame.table)
    phi = frame.analysis.apply(psi)
    assert frame.basis_space.inner_product(phi, phi) == pytest.approx(
        frame.measure_space.inner_product(psi, psi), rel=1e-12
    )
    live_l2 = np.diagonal(frame.discrete_gram).reshape(frame.basis.space.shape)[2]
    np.testing.assert_allclose(live_l2[live_l2 > 0.0], [0.4, 0.8, 0.8], rtol=1e-12)
    w_total = float(frame.measure.weights.sum())
    y = rng.standard_normal(frame.basis_space.shape)
    lhs = frame.analysis.H.apply(y)
    rhs = frame.reconstruction.apply(y) / w_total
    rel = float(np.max(np.abs(lhs - rhs)) / np.max(np.abs(rhs)))
    assert rel > 0.5, (
        f"the sphere collapse unexpectedly held on the slab (rel {rel:.3g})"
    )


@pytest.mark.foundation
def test_the_dense_dressing_reds_under_the_diagonal_and_the_pre_repair_metrics():
    r"""D4 — the DENSE arm's loadedness witness, modelled on the diagonal
    arm's (pre-seed the cached space with a wrong metric; the isometry
    must FAIL). Two wrong metrics, two floors, both `[M]` on this seed:
    the pre-repair continuum reads 25.53 (floor 10) and the best diagonal
    candidate reads 1.806 (floor 1.5). The dressed reading passes at
    1e-12 (D2) — the family is loaded, not blind.
    """
    for label, floor in (("continuum", 10.0), ("diagonal", 1.5)):
        frame = Quadrature.gauss_legendre(8).angular_frame(2)
        if label == "continuum":
            wrong = SphericalHarmonicSpace.from_L(2)
        else:
            diag = np.diagonal(frame.discrete_gram).reshape(
                frame.basis.space.shape
            )
            live = diag > 0.0
            wrong = replace(
                frame.basis.space,
                inner_product_weights=np.where(
                    live, 1.0 / np.where(live, diag, 1.0), 0.0
                ),
            )
        vars(frame)["basis_space"] = wrong
        vars(frame)["test_space"] = wrong
        rng = np.random.default_rng(1234)
        c = rng.standard_normal(wrong.shape)
        psi = frame.basis.synthesize(c, frame.table)
        phi = frame.analysis.apply(psi)
        ratio = frame.basis_space.inner_product(phi, phi) / (
            frame.measure_space.inner_product(psi, psi)
        )
        assert ratio > floor, (
            f"{label} metric read Parseval {ratio:.3g} — the gate would be "
            f"blind to it"
        )


@pytest.mark.l1
@pytest.mark.catches("ERR-039")
def test_the_dressing_lands_parseval_on_the_production_anisotropic_frame():
    r"""D5 — the tier-of-observability gate for the production adjoint move.

    ``product(4,4).angular_frame(2)`` is a production-reachable
    ScatteringOperator configuration (``scattering.py`` builds
    ``quadrature.angular_frame(scattering_order)``) and its Gram measures
    DENSE. `[M]` 2026-08-30 (design pre-flight): dressing it moves
    ``frame.analysis.H`` by ``max|Δ| = 8.246`` — **rel 0.8995** on that
    probe's draw; ⚠ draw-banded 0.879–0.986 over 200 seeds, with the
    draw-free operator-level Frobenius relative at **0.980–0.985**
    across three DENSE frames (archivist census) — the recorded F-0
    limitation repaired (the undressed ``.H`` was the stored-metric
    sandwich, NOT the physical Hilbert adjoint), and NOTHING else in
    the 4371-test pre-flight scope observed it (plan-authoring §8,
    measured). This gate is where the change is visible: Parseval holds
    post-dressing (`[M]` 1.000000000000 this seed) and the pre-repair
    continuum metric reads 65.66.
    """
    frame = Quadrature.product(4, 4).angular_frame(2)
    assert frame.discrete_gram_structure is GramStructure.DENSE
    rng = np.random.default_rng(1234)
    c = rng.standard_normal(frame.basis.space.shape)
    psi = frame.basis.synthesize(c, frame.table)
    phi = frame.analysis.apply(psi)
    norm_w = frame.measure_space.inner_product(psi, psi)
    assert frame.basis_space.inner_product(phi, phi) / norm_w == pytest.approx(
        1.0, rel=1e-12
    )
    pre_repair = frame.basis.space.inner_product(phi, phi) / norm_w
    assert pre_repair > 10.0, f"pre-repair metric read {pre_repair:.4f}"
