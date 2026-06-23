r"""Intrinsic-property tests for :class:`orpheus.numerics.frame.Frame`.

The :class:`Frame` binds a :class:`Basis` to a :class:`DiscreteMeasure` and exposes
the analysis / reconstruction faces. These tests pin:

* **bit-identity** of the faces to the operators they retire
  (:class:`MomentProjection` / :class:`HarmonicMomentReconstruction`) — the carve is
  a re-homing, not a re-derivation, so ``np.array_equal`` is the contract (the
  scattering kernel is pinned at 0 ULP);
* the **adjoint-for-free** — ``frame.analysis.H`` falls out of the frame's swapped
  spaces with no bespoke code;
* the symmetric **space pairing** (basis → ``basis_space``; measure → ``measure_space``);
* the faces **compose through ``OperatorProduct`` with real spaces** (no ``cast``) —
  the enabler for the Phase-C cast retirement.

The structural Galerkin invariant :math:`\Pi R = 4\pi I` itself lives in
``test_spherical_harmonic_space.py`` (and is inherited here by bit-identity).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.basis import SphericalHarmonicBasis
from orpheus.numerics.frame import Frame
from orpheus.numerics.operator import CAP_APPLY, CAP_APPLY_TRANSPOSE
from orpheus.numerics.projection import (
    HarmonicMomentReconstruction,
    MomentProjection,
)
from orpheus.numerics.quadrature import lebedev_sphere
from orpheus.numerics.spaces import SphericalHarmonicSpace


@pytest.fixture
def sh_frame():
    """An exact spherical-harmonic frame: Lebedev(13) ⋈ SH(L=3)."""
    measure = lebedev_sphere(13)
    L = 3
    basis = SphericalHarmonicBasis(L=L)
    return Frame(basis, measure), L


def _band_limited(rng, L, *trailing):
    """A random ``(L+1, 2L+1, *trailing)`` moment array with |m|>ℓ slots zeroed."""
    c = rng.standard_normal((L + 1, 2 * L + 1, *trailing))
    for ell in range(L + 1):
        c[ell, 2 * ell + 1 :] = 0.0
    return c


# ── bit-identity of the faces to the operators they retire ────────────────

@pytest.mark.foundation
def test_analysis_face_is_bit_identical_to_moment_projection(sh_frame):
    frame, L = sh_frame
    rng = np.random.default_rng(11)
    psi = rng.standard_normal((frame.measure.weights.shape[0], 4, 2))  # (N, cell, group)
    M = MomentProjection(weights=frame.measure.weights, Y=frame.table, L=L)
    np.testing.assert_array_equal(frame.analysis.apply(psi), M.apply(psi))


@pytest.mark.foundation
def test_reconstruction_face_is_bit_identical_to_harmonic_reconstruction(sh_frame):
    frame, L = sh_frame
    rng = np.random.default_rng(12)
    c = _band_limited(rng, L, 4, 2)
    R = HarmonicMomentReconstruction.from_Y(frame.table)
    np.testing.assert_array_equal(frame.reconstruction.apply(c), R.apply(c))


@pytest.mark.foundation
def test_analysis_transpose_is_bit_identical_to_moment_projection(sh_frame):
    frame, L = sh_frame
    rng = np.random.default_rng(13)
    c = _band_limited(rng, L, 4, 2)
    M = MomentProjection(weights=frame.measure.weights, Y=frame.table, L=L)
    np.testing.assert_array_equal(
        frame.analysis.apply_transpose(c), M.apply_transpose(c),
    )


# ── the adjoint-for-free ──────────────────────────────────────────────────

@pytest.mark.foundation
def test_analysis_hilbert_adjoint_falls_out_of_the_frame_spaces(sh_frame):
    """``frame.analysis.H`` (= g_C·S₀) matches ``MomentProjection.H`` bit-identically.

    No bespoke adjoint code — the frame's swapped (measure_space, basis_space)
    metrics feed the generic ``_AdjointOperator`` exactly as the legacy operator's
    domain/codomain did.
    """
    frame, L = sh_frame
    rng = np.random.default_rng(14)
    c = _band_limited(rng, L, 4, 2)
    M = MomentProjection.from_measure(frame.measure, L, Y=frame.table)
    np.testing.assert_array_equal(
        frame.analysis.H.apply(c), M.H.apply(c),
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
    # reconstruction is apply-only in Phase B; apply_transpose (→ R.H) lands in Phase D
    assert frame.reconstruction.capabilities == frozenset({CAP_APPLY})
