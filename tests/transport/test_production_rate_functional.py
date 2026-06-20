r"""#257 S5 — Spec B: ``ProductionRateFunctional.evaluate`` correctness.

The functional computes the per-cell fission emission DENSITY (group-axis
contraction):

.. math::

    p(\vec r) \;=\; \sum_{g'} \nu\Sigma_{f,g'}(\vec r)\,\phi_{g'}(\vec r)

This is EXACTLY the anonymous ``inner`` einsum currently inside
``RankOneOperator.apply`` (``operator.py``):
``inner = (self.right * x).sum(axis=self.axis, keepdims=True)`` with
``right = νΣf``, ``axis = 0`` (group). In S5 the fission operator is NOT
rewired (that is S6) — the functional is created and cross-checked here,
ready for S6 to compose ``F = M_χ ∘ ProductionRateFunctional ∘ M_νΣf``.

vv claim layer (1.5 gate): a FLUX-SHAPE / value claim (``evaluate`` value
on a fixed input), NOT an eigenvalue claim. The correctness reference is
the STRUCTURALLY-INDEPENDENT hand-derived double-loop
:func:`hand_derived_production_density` — an explicit per-cell Python sum
that shares NO reduction primitive with the production code and NONE of
the ORPHEUS operator algebra. The bit-identity-vs-``RankOneOperator``
leg is the EQUIVALENCE check (de-risks the S6 composition) and is clearly
demarcated as such — it is NOT the correctness reference (L11: testing
only against the operator the SUT will replace is cross-implementation
agreement, not correctness).

Failure-mode exposure:

* Mode-2 (variable swap νΣf↔φ / wrong-axis contraction): the inputs are
  ASYMMETRIC ≥2G HETEROGENEOUS (distinct per group AND per cell, distinct
  ramp patterns) so a swap or a transposed contraction axis yields a
  DIFFERENT value than the hand loop. A flat or symmetric input would
  null this (H1/H2).
* Mode-3 (missing/extra factor — a measure fold): the density is
  group-collapsed but carries NO volume measure; the shape/units pin
  asserts the SUT did not accidentally multiply by ``V`` (that is a
  SEPARATE concern — the production-RATE integral lives in the
  estimators, Spec C).

vv Mode-8: structural assertions route through ``require``; value
assertions through ``np.testing.*`` — both fire under ``python -O``.

``foundation`` — L0 value verification against a hand-derived reference,
no theory-page ``:label:``.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.operator import RankOneOperator

from tests.transport._functional_helpers import (
    asymmetric_nu_sigma_f,
    asymmetric_phi,
    build_production_rate_functional,
    cartesian_2d_mesh,
    cross_section_field,
    hand_derived_production_density,
    require,
    slab_mesh,
    squeeze_density,
)

pytestmark = pytest.mark.foundation


# ───────────────────────────────────────────────────────────────────────
# Helper: build SUT + the matching phi field on a chosen mesh.
# ───────────────────────────────────────────────────────────────────────


def _setup(sn, ng):
    nu_sf = asymmetric_nu_sigma_f(ng=ng, spatial_shape=sn.spatial_shape)
    phi = asymmetric_phi(ng=ng, spatial_shape=sn.spatial_shape)
    func = build_production_rate_functional(cross_section_field(nu_sf, sn))
    return func, nu_sf, phi


def _evaluate(func, phi: np.ndarray) -> np.ndarray:
    """Call the SUT's ``evaluate`` with the raw ``(ng, *spatial)`` flux array.

    The SUT contract is ``evaluate(x: V) -> float | V``. The exact
    accepted carrier (raw ndarray vs a typed ScalarFlux/AngularFlux) is
    the method-implementer's latitude; the brief's value description is
    over the raw ``(ng, *spatial)`` array, so the spec passes the raw
    array and reads ``.values`` off the result if the SUT returned a
    typed field. Squeezing happens in the caller.
    """
    out = func.evaluate(phi)
    return np.asarray(getattr(out, "values", out))


# ═══════════════════════════════════════════════════════════════════════
# B.1 — CORRECTNESS reference: the hand-derived explicit double-loop.
# ═══════════════════════════════════════════════════════════════════════


class TestEvaluateAgainstHandDerivedReference:
    """The load-bearing correctness leg — structurally-independent ref."""

    def test_evaluate_matches_hand_loop_2g_heterogeneous_2d(self):
        """evaluate(φ) == explicit Python Σ_g νΣf_g·φ_g (2G, het, nx≠ny).

        The hand loop is a different STRUCTURAL path (explicit scalar
        accumulation, no numpy reduction). Agreement to high precision
        proves the SUT contracts the RIGHT axis with the RIGHT operand
        and applies NO spurious factor. A νΣf↔φ swap or a wrong-axis
        contraction disagrees with this loop on the asymmetric inputs.
        """
        sn = cartesian_2d_mesh(nx=5, ny=3, ng=2)
        func, nu_sf, phi = _setup(sn, ng=2)

        got = squeeze_density(_evaluate(func, phi))
        ref = squeeze_density(hand_derived_production_density(nu_sf, phi))

        require(
            got.shape == ref.shape,
            f"evaluate shape {got.shape} != hand-derived shape {ref.shape}.",
        )
        # Tolerance: the SUT may reduce with numpy (pairwise) and the loop
        # accumulates left-to-right. The reduction depth is ng=2 (a single
        # add per cell) so the FP drift is ≤ a few ULP. nulp=8 is generous
        # headroom; the values are O(1)–O(10).
        np.testing.assert_array_almost_equal_nulp(got, ref, nulp=8)

    def test_evaluate_matches_hand_loop_4g_heterogeneous_slab(self):
        """4G heterogeneous slab — more groups deepen the contraction.

        ng=4 makes the group-axis contraction non-trivial (4 terms per
        cell) so a partial-sum / off-by-one-group error surfaces that a
        2G case (single add) cannot see.
        """
        sn = slab_mesh(nx=6, ng=4)
        func, nu_sf, phi = _setup(sn, ng=4)

        got = squeeze_density(_evaluate(func, phi))
        ref = squeeze_density(hand_derived_production_density(nu_sf, phi))
        np.testing.assert_array_almost_equal_nulp(got, ref, nulp=16)

    def test_wrong_axis_contraction_disagrees_with_sut(self):
        """MUTATION (mode #2): the wrong contraction axis is detectably wrong.

        The genuine mode-2 hazard for this functional is contracting the
        SPATIAL axis instead of the GROUP axis (the pointwise product
        νΣf·φ is commutative, so a literal νΣf↔φ swap is value-invariant —
        the real trap is the axis). On the nx≠ny mesh a spatial-axis
        contraction produces a DIFFERENT SHAPE than the SUT's
        group-collapsed output, AND, where shapes happen to align, a
        DIFFERENT VALUE. This row proves the test bench discriminates the
        axis: the SUT's output matches the group-axis loop and does NOT
        match a spatial-axis contraction.
        """
        sn = cartesian_2d_mesh(nx=5, ny=3, ng=2)
        func, nu_sf, phi = _setup(sn, ng=2)
        got = squeeze_density(_evaluate(func, phi))

        # Group-axis (correct) vs an x-axis contraction (mode-2 wrong).
        group_axis = squeeze_density(
            hand_derived_production_density(nu_sf, phi)
        )
        wrong_axis = (nu_sf * phi).sum(axis=1)  # sum over x, not group

        np.testing.assert_array_almost_equal_nulp(got, group_axis, nulp=8)
        require(
            got.shape != wrong_axis.shape,
            f"On the nx≠ny mesh the group-collapsed output {got.shape} must "
            f"NOT share the spatial-axis-contracted shape {wrong_axis.shape} "
            f"— the axis discriminator a square mesh would hide.",
        )


# ═══════════════════════════════════════════════════════════════════════
# B.2 — EQUIVALENCE check (NOT correctness): bit-identity vs the live
#       RankOneOperator.inner extraction. De-risks the S6 composition.
# ═══════════════════════════════════════════════════════════════════════


class TestEvaluateEquivalentToRankOneInner:
    """Cross-implementation agreement with the operator S6 will compose.

    DEMARCATED as the EQUIVALENCE check, not the correctness reference.
    The correctness reference is the hand loop (B.1). This leg exists
    only to guarantee that when S6 wires
    ``F = M_χ ∘ ProductionRateFunctional ∘ M_νΣf``, the functional
    reproduces the EXACT bytes the legacy ``inner`` produced — so the
    S6 rewire is bit-identical by construction.
    """

    def test_evaluate_bit_identical_to_rank_one_inner_2d(self):
        """evaluate(φ) is 0-ULP identical to ``(νΣf*φ).sum(axis=0, keepdims=True)``.

        The production extraction is the ``inner`` line of
        ``RankOneOperator.apply``: ``(self.right * x).sum(axis=0,
        keepdims=True)`` with ``right = νΣf``. The functional MUST emit
        the same reduction tree (same numpy primitive, same axis, same
        order) so the S6 composition inherits bit-identity. Expect 0 ULP
        (``assert_array_equal``); if a future re-association appears,
        relax per vv-principles' 3 criteria with a documented justification.
        """
        sn = cartesian_2d_mesh(nx=5, ny=3, ng=2)
        func, nu_sf, phi = _setup(sn, ng=2)

        # The EXACT live extraction: build the RankOneOperator the way
        # FissionOperator.kernel does (chi is irrelevant to ``inner``;
        # we read the inner contraction directly via .right/.apply
        # structure). The ``inner`` is the (right * x).sum(axis,keepdims).
        chi = np.ones_like(nu_sf)
        rank_one = RankOneOperator(chi, nu_sf, axis=0)
        legacy_inner = (rank_one.right * phi).sum(
            axis=rank_one.axis, keepdims=True
        )

        got = _evaluate(func, phi)
        # LITERAL-RANK byte-identity (NOT squeezed): the SUT commits to
        # ``keepdims=True`` → ``(1, *spatial)``, reproducing the legacy
        # ``inner`` EXACTLY including rank. Pinning the un-squeezed arrays
        # makes a future ``keepdims`` change redden — the S6 composition
        # ``M_χ ∘ this`` broadcasts χ against this leading-1 axis, so the
        # rank is load-bearing for the bit-identity de-risk, not incidental.
        # (B.1's hand loop is the squeeze-agnostic VALUE reference; THIS leg
        # pins the literal reduction tree the docstring claims.)
        np.testing.assert_array_equal(got, legacy_inner)


# ═══════════════════════════════════════════════════════════════════════
# B.3 — Mode-3 guard: NO volume measure folded into the density.
# ═══════════════════════════════════════════════════════════════════════


class TestNoVolumeMeasureFolded:
    """The density is group-collapsed but carries NO volume weight.

    A later accidental measure-fold (``p·V``) would redden this. The
    volume INTEGRAL is a separate concern (the production-RATE estimator,
    Spec C) — the functional must stay the per-cell density.
    """

    def test_density_shape_is_group_collapsed_only(self):
        """Output collapses the group axis ONLY; spatial axes are preserved.

        Shape contract: ``(ng, *spatial) -> (1, *spatial)`` or
        ``(*spatial,)``. A measure-fold that sums over space would change
        the spatial shape (a scalar or a reduced grid) — caught here.
        """
        sn = cartesian_2d_mesh(nx=5, ny=3, ng=2)
        func, nu_sf, phi = _setup(sn, ng=2)
        got = squeeze_density(_evaluate(func, phi))
        require(
            got.shape == sn.spatial_shape,
            f"Production density must preserve the spatial shape "
            f"{sn.spatial_shape} (group axis collapsed only); got {got.shape}. "
            f"A spatial reduction here signals a spurious volume-measure fold.",
        )

    def test_density_unweighted_by_cell_volume(self):
        """The per-cell value equals Σ_g νΣf_g·φ_g — NOT Σ_g νΣf_g·φ_g·V.

        Direct value check: the SUT's per-cell number must match the
        UNWEIGHTED hand loop. If the SUT folded the cell volume, the
        values would differ by the per-cell ``V`` (non-uniform on the
        nx≠ny mesh: Δx = 2/5 ≠ Δy = 1/3, so V varies cell-to-cell and the
        discrepancy is detectable). Asserting equality to the unweighted
        loop is the measure-fold tripwire.
        """
        sn = cartesian_2d_mesh(nx=5, ny=3, ng=2)
        func, nu_sf, phi = _setup(sn, ng=2)
        got = squeeze_density(_evaluate(func, phi))
        unweighted = squeeze_density(
            hand_derived_production_density(nu_sf, phi)
        )
        np.testing.assert_array_almost_equal_nulp(got, unweighted, nulp=8)
