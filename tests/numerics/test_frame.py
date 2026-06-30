r"""Intrinsic-property tests for the frame faces + the discipline-type hierarchy.

The angular :class:`~orpheus.numerics.frame.GalerkinFrame` binds a :class:`Basis` to a
:class:`DiscreteMeasure` and exposes the analysis / reconstruction faces. These tests pin:

* the **adjoint-for-free** — BOTH ``frame.analysis.H`` and ``frame.reconstruction.H``
  fall out of the frame's swapped spaces with no bespoke code (each pinned against an
  INDEPENDENT closed-form einsum: :math:`g_C \cdot S_0` for analysis,
  :math:`\frac{(2\ell+1)^2}{4\pi} M` for reconstruction);
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
from orpheus.numerics.frame import FrameBase, GalerkinFrame, PetrovGalerkinFrame
from orpheus.numerics.measure import DiscreteMeasure
from orpheus.numerics.operator import (
    CAP_APPLY,
    CAP_APPLY_TRANSPOSE,
    MissingCapability,
)
from orpheus.numerics.quadrature import lebedev_sphere
from orpheus.numerics.spaces import SphericalHarmonicSpace
from tests._harness.predicates import assert_capability_faithful


@pytest.fixture
def sh_frame():
    """An exact spherical-harmonic frame: Lebedev(13) ⋈ SH(L=3)."""
    measure = lebedev_sphere(13)
    L = 3
    basis = SphericalHarmonicBasis(L=L)
    return GalerkinFrame(basis, measure), L


@pytest.mark.foundation
def test_frame_faces_predicates_faithful_to_caps(sh_frame):
    r"""The carve keystone (#226) for the numerics frame faces.

    Both faces advertise ``apply_transpose`` — so ``.H`` falls out of the
    metric-aware ``_AdjointOperator`` — hence ``is_adjointable`` is True and
    mirrors ``capabilities`` EXACTLY; neither face is invertible. The
    reconstruction face's ``is_adjointable`` is the OVERRIDE that lifts it above
    the bare :class:`~orpheus.numerics.projection.ReconstructionOperator`
    default (caps ``{apply}`` → ``is_adjointable`` False); the analysis face
    inherits ``True`` from :class:`~orpheus.numerics.projection.AnalysisOperator`.
    """
    frame, _ = sh_frame
    for face in (frame.analysis, frame.reconstruction):
        assert_capability_faithful(face)
        assert face.is_adjointable is True
        assert face.is_invertible is False


def _band_limited(rng, L, *trailing):
    """A random ``(L+1, 2L+1, *trailing)`` moment array with |m|>ℓ slots zeroed."""
    c = rng.standard_normal((L + 1, 2 * L + 1, *trailing))
    for ell in range(L + 1):
        c[ell, 2 * ell + 1 :] = 0.0
    return c


# ── the adjoint-for-free ──────────────────────────────────────────────────

@pytest.mark.foundation
def test_analysis_hilbert_adjoint_falls_out_of_the_frame_spaces(sh_frame):
    r"""``frame.analysis.H`` is the W-weighted Hilbert adjoint :math:`g_C \cdot S_0`.

    No bespoke adjoint code — the frame's swapped ``(measure_space, basis_space)``
    metrics feed the generic ``_AdjointOperator``. Pinned against an INDEPENDENT
    reference: the direct :math:`g_C \cdot S_0` einsum with the closed-form SH Gram
    diagonal :math:`g_C^\ell = 4\pi/(2\ell+1)` (NOT the frame's own contraction —
    a structurally distinct construction of the same adjoint).
    """
    frame, L = sh_frame
    rng = np.random.default_rng(14)
    c = _band_limited(rng, L, 4, 2)
    g_C = 4.0 * np.pi / (2.0 * np.arange(L + 1) + 1.0)  # closed-form SH Gram diag
    expected = np.einsum("nlm,l,lm...->n...", frame.table, g_C, c)
    np.testing.assert_allclose(
        frame.analysis.H.apply(c), expected, rtol=1e-12, atol=1e-14,
    )


@pytest.mark.foundation
def test_reconstruction_hilbert_adjoint_falls_out_of_the_frame_spaces(sh_frame):
    r"""``frame.reconstruction.H`` is the W-weighted Hilbert adjoint of :math:`R`.

    Symmetric with the analysis face: no bespoke adjoint code — the frame's swapped
    ``(basis_space, measure_space)`` metrics feed the generic ``_AdjointOperator``,
    giving :math:`(R^* v)_\ell^m = \frac{(2\ell+1)^2}{4\pi} \sum_n w_n
    Y_\ell^m(\hat\Omega_n)\, v_n`. Pinned against that INDEPENDENT closed-form einsum
    (NOT the frame's own contraction). ``R : \text{basis} \to \text{measure}``, so
    ``R.H`` maps nodal values → coefficients.
    """
    frame, L = sh_frame
    rng = np.random.default_rng(17)
    n = frame.measure.weights.shape[0]
    v = rng.standard_normal((n, 4, 2))
    factor = (2.0 * np.arange(L + 1) + 1.0) ** 2 / (4.0 * np.pi)  # (2ℓ+1)²/4π
    expected = np.einsum(
        "n,nlm,l,n...->lm...", frame.measure.weights, frame.table, factor, v,
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
    lets Phase C drop the ``cast(LinearOperator, …)`` workarounds in the scattering
    kernel.
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


@pytest.mark.foundation
def test_face_capabilities(sh_frame):
    frame, _ = sh_frame
    assert frame.analysis.capabilities == frozenset({CAP_APPLY, CAP_APPLY_TRANSPOSE})
    # Phase D: reconstruction now advertises apply_transpose too (→ R.H falls out),
    # symmetric with the analysis face.
    assert frame.reconstruction.capabilities == frozenset(
        {CAP_APPLY, CAP_APPLY_TRANSPOSE}
    )


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
    :class:`MissingCapability` (the dense ``(MR)⁻¹M`` solve is unbuilt — #275). Mutation
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
    with pytest.raises(MissingCapability, match="DENSE"):
        _ = dense.gram
    with pytest.raises(MissingCapability, match="DENSE"):
        dense.project(np.array([3.0, 5.0]))
    # Control: the SAME geometry with the honest DIAGONAL trial projects fine.
    ok = PetrovGalerkinFrame(
        IndicatorBasis((edges,)), measure,
        WeightedIndicatorBasis(IndicatorBasis((edges,)), np.ones(2)),
    )
    np.testing.assert_allclose(ok.project(np.array([3.0, 5.0])), [3.0, 5.0])
