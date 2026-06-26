r"""``ReactionRateFunctional`` — ``⟨Σx, φ⟩`` correctness + the closed-form k∞ oracle.

The generalisation (and successor) of ``test_production_rate_functional.py``:
``ReactionRateFunctional(Σx)`` is the reaction-rate co-vector, with the
production (``Σx = νΣf``) and absorption (``Σx = Σa``) instances. Two legs:

* **A — hand-loop value correctness** (heterogeneous, asymmetric, ≥2G): the
  per-cell density ``Σ_g Σx_g·φ_g`` matches an explicit Python double-loop
  (``hand_derived_production_density``) — a structurally-independent reference
  sharing no reduction primitive with the SUT. Generalises the retired
  ProductionRateFunctional B.1.

* **B — the closed-form k∞ per-term oracle (the headline upgrade).** This
  REPLACES the retired procedural twin (``χ·evaluate ≡ apply``, the same numpy
  primitive both sides → vv-principles L11 weak) with a structurally-independent
  ground: the infinite-medium decomposition ``k∞ = λ_max(A⁻¹F)`` (an
  ``np.linalg.eig`` of the transfer matrix, NOT the rank-1 path). The
  **production AND absorption functionals are pinned INDEPENDENTLY** against
  ``⟨νΣf,φ*⟩`` / ``⟨Σa,φ*⟩`` at the converged spectrum φ*, so a shared-factor
  error that the ratio ``k = PROD/ABSN`` would mask (a mis-scaled accessor, a
  spurious volume fold on BOTH) is caught. Mixture A 2g + 4g; **4g is mandatory**
  — A 2g's spectrum is coincidentally flat (``[0.707, 0.707]``), so its
  flat-flux ratio equals k∞ and is flux-shape-blind; A 4g's φ* is genuinely
  non-flat and carries the flux-shape teeth.

* **C — no volume measure** (Mode-3): the functional returns the per-cell
  density, NOT ``∫·dV`` (the volume-integrated rate is a separate spatial-measure
  concern — the Phase-4 ``compute_keff`` fold).

vv claim-layer (1.5 gate): Leg B pairs an **eigenvalue** claim with the
**closed-form** pillar (``λ_max(A⁻¹F)``) — legal; all chains terminate in
``np.linalg.eig``, structurally independent of the SUT. ``foundation`` — L0
value correctness + software invariants, no theory ``:label:``.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.eigenvalue import kinf_and_spectrum_homogeneous
from orpheus.derivations.common.xs_library import get_mixture
from orpheus.numerics.functional import Functional, InnerProductFunctional
from orpheus.numerics.operator import LinearOperator
from orpheus.transport.reaction_rate_functional import ReactionRateFunctional

from tests.transport._functional_helpers import (
    asymmetric_nu_sigma_f,
    asymmetric_phi,
    cartesian_2d_mesh,
    cross_section_field,
    hand_derived_production_density,
    require,
    slab_mesh,
    squeeze_density,
)

pytestmark = pytest.mark.foundation


def _hand_dot(weight: np.ndarray, phi: np.ndarray) -> float:
    """``⟨weight, phi⟩`` over groups by an EXPLICIT Python loop (NOT ``@`` / ``.sum``).

    A ``@`` would re-share numpy's reduction with the functional's ``.sum`` —
    the procedural-twin risk recurs. The loop is the structurally-independent
    reference (mirrors ``hand_derived_production_density``).
    """
    return float(sum(weight[g] * phi[g] for g in range(len(weight))))


def _homogeneous_setup(ng_key: str):
    """Mixture A → ``(νΣf, Σa, φ*, k∞)`` via the closed-form ``λ_max(A⁻¹F)``."""
    mix = get_mixture("A", ng_key)
    sigp = np.asarray(mix.SigP)
    siga = np.asarray(mix.absorption_xs)
    sigt = np.asarray(mix.SigT)
    sigs0 = np.asarray(mix.SigS[0].todense())
    chi = np.asarray(mix.chi)
    k, phi = kinf_and_spectrum_homogeneous(sigt, sigs0, sigp, chi)
    return sigp, siga, phi, k


# ═══════════════════════════════════════════════════════════════════════
# A — CORRECTNESS reference: the hand-derived explicit double-loop.
# ═══════════════════════════════════════════════════════════════════════


class TestEvaluateAgainstHandDerivedReference:
    """Heterogeneous, asymmetric ≥2G — the structurally-independent value leg."""

    def test_evaluate_matches_hand_loop_2g_heterogeneous_2d(self):
        sn = cartesian_2d_mesh(nx=5, ny=3, ng=2)
        sigx = asymmetric_nu_sigma_f(ng=2, spatial_shape=sn.spatial_shape)
        phi = asymmetric_phi(ng=2, spatial_shape=sn.spatial_shape)
        func = ReactionRateFunctional(cross_section_field(sigx, sn))

        got = squeeze_density(func.evaluate(phi))
        ref = squeeze_density(hand_derived_production_density(sigx, phi))
        require(got.shape == ref.shape, f"evaluate shape {got.shape} != ref {ref.shape}.")
        np.testing.assert_array_almost_equal_nulp(got, ref, nulp=8)

    def test_evaluate_matches_hand_loop_4g_heterogeneous_slab(self):
        sn = slab_mesh(nx=6, ng=4)
        sigx = asymmetric_nu_sigma_f(ng=4, spatial_shape=sn.spatial_shape)
        phi = asymmetric_phi(ng=4, spatial_shape=sn.spatial_shape)
        func = ReactionRateFunctional(cross_section_field(sigx, sn))

        got = squeeze_density(func.evaluate(phi))
        ref = squeeze_density(hand_derived_production_density(sigx, phi))
        np.testing.assert_array_almost_equal_nulp(got, ref, nulp=16)


# ═══════════════════════════════════════════════════════════════════════
# B — the closed-form k∞ per-term oracle (the structural-independence UPGRADE).
# ═══════════════════════════════════════════════════════════════════════


class TestClosedFormKinfPerTermOracle:
    """Production AND absorption pinned INDEPENDENTLY vs ``λ_max(A⁻¹F)``."""

    @pytest.mark.parametrize(
        "ng_key, prod_snapshot, absn_snapshot, kinf_snapshot",
        [
            ("2g", 0.15909903, 0.08485281, 1.8750000000),
            ("4g", 0.19041882, 0.12799012, 1.4877619048),
        ],
    )
    def test_production_and_absorption_pinned_independently(
        self, ng_key, prod_snapshot, absn_snapshot, kinf_snapshot
    ):
        sigp, siga, phi, k = _homogeneous_setup(ng_key)
        ng = len(sigp)
        sn = slab_mesh(nx=1, ng=ng)
        phi_field = phi[:, None]  # flat in space, (ng, 1)

        prod = squeeze_density(
            ReactionRateFunctional(cross_section_field(sigp[:, None], sn)).evaluate(phi_field)
        ).item()
        absn = squeeze_density(
            ReactionRateFunctional(cross_section_field(siga[:, None], sn)).evaluate(phi_field)
        ).item()

        # PRIMARY: each functional pinned INDEPENDENTLY vs the Python-loop dot
        # (the structurally-independent reference; catches a shared-factor error
        # the ratio masks).
        np.testing.assert_allclose(prod, _hand_dot(sigp, phi), rtol=1e-13)
        np.testing.assert_allclose(absn, _hand_dot(siga, phi), rtol=1e-13)
        # COMPOSED: the Rayleigh quotient reproduces the closed-form eigenvalue.
        np.testing.assert_allclose(prod / absn, k, rtol=1e-10)
        # REGRESSION LOCK: end-to-end (mixture → eig → dot) snapshot. Guards a
        # drift in φ* / the closed-form that the live legs above would track.
        np.testing.assert_allclose([prod, absn, k], [prod_snapshot, absn_snapshot, kinf_snapshot], rtol=1e-7)

    def test_4g_spectrum_is_non_flat_so_4g_carries_flux_shape_teeth(self):
        """The config-blindness guard: A 4g's flat-flux ratio ≠ k∞.

        A flat φ (or A 2g's coincidentally flat spectrum) makes ``Σνσf/ΣΣa``
        equal k∞ — flux-shape-blind. The 4g φ* is genuinely non-flat, so its
        flat-flux ratio DIFFERS from k∞; that gap is what makes the 4g per-term
        pin able to catch a group-weight error.
        """
        sigp_2g, siga_2g, _, k2 = _homogeneous_setup("2g")
        sigp_4g, siga_4g, _, k4 = _homogeneous_setup("4g")
        flat_2g = float(sigp_2g.sum() / siga_2g.sum())
        flat_4g = float(sigp_4g.sum() / siga_4g.sum())
        require(
            np.isclose(flat_2g, k2, rtol=1e-3),
            f"A 2g spectrum is expected coincidentally flat (flat ratio {flat_2g:.4f} "
            f"≈ k∞ {k2:.4f}) — documents WHY 2g alone is flux-shape-blind.",
        )
        require(
            not np.isclose(flat_4g, k4, rtol=1e-3),
            f"A 4g flat-flux ratio {flat_4g:.4f} must DIFFER from k∞ {k4:.4f} "
            f"— else the mandatory 4g case would be flux-shape-blind too.",
        )


# ═══════════════════════════════════════════════════════════════════════
# C — Mode-3 guard: NO volume measure folded into the density.
# ═══════════════════════════════════════════════════════════════════════


class TestNoVolumeMeasureFolded:
    """The density is group-collapsed but carries NO volume weight."""

    def test_density_unweighted_by_cell_volume(self):
        sn = cartesian_2d_mesh(nx=5, ny=3, ng=2)  # Δx=2/5 ≠ Δy=1/3 → V varies per cell
        sigx = asymmetric_nu_sigma_f(ng=2, spatial_shape=sn.spatial_shape)
        phi = asymmetric_phi(ng=2, spatial_shape=sn.spatial_shape)
        got = squeeze_density(ReactionRateFunctional(cross_section_field(sigx, sn)).evaluate(phi))
        unweighted = squeeze_density(hand_derived_production_density(sigx, phi))
        require(
            got.shape == sn.spatial_shape,
            f"density must preserve the spatial shape {sn.spatial_shape} (group axis "
            f"collapsed only); got {got.shape} — a spatial reduction signals a volume fold.",
        )
        np.testing.assert_array_almost_equal_nulp(got, unweighted, nulp=8)


# ═══════════════════════════════════════════════════════════════════════
# Category membership — a Functional that specialises InnerProductFunctional.
# ═══════════════════════════════════════════════════════════════════════


class TestCategoryMembership:
    def test_reaction_rate_is_functional_specialising_inner_product(self):
        sn = slab_mesh(nx=1, ng=2)
        f = ReactionRateFunctional(cross_section_field(np.array([[0.025], [0.20]]), sn))
        require(isinstance(f, Functional), "must satisfy the Functional Protocol.")
        require(
            isinstance(f, InnerProductFunctional),
            "must specialise InnerProductFunctional (the HarmonicFrame(GalerkinFrame) pattern).",
        )
        require(
            not isinstance(f, LinearOperator),
            "a Functional must NOT satisfy LinearOperator (no apply/capabilities).",
        )
