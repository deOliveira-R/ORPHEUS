r"""Intrinsic gates for the model-independent isotropic energy operators (#276 / K_iso).

:class:`IsotropicScattering` (:math:`\Sigma_{s,0}`) and :class:`IsotropicN2N`
(:math:`2\Sigma_{2n}`) are the scalar-flux realization of the isotropic in-scatter
source shared across ALL transport models. These gates pin them in isolation
(synthetic per-material XS — the operators touch ONLY the per-material dispatch;
the real-mesh SN integration is the P2 forward / P3 adjoint bit-identity gates):

* **apply ≡ fast-path** — the OperatorSum re-expression of ``add_iso_source`` +
  ``add_n2n_source`` is byte-identical (same ``mat_xs`` verbs, same accumulation).
* **transpose** — Euclidean reciprocity ``⟨Kφ,χ⟩=⟨φ,Kᵀχ⟩`` + a structurally-
  independent dense ``Σ @ vec`` per-material reference (a wrong group axis reds).
* **as_dense ≡ apply** — the two consumption modes agree (the LHS-fold view).
* **predicates** — ``is_adjointable`` True; STRUCTURALLY non-invertible
  (``is_invertible`` False and no ``inverse`` method at all).

Asymmetric ``SigS`` + non-vacuous ``Sig2`` (≠ per material) so a group-axis flip or
a dropped n2n channel is detectable (L27 / #269). LD multi-moment (trailing 2^d)
included. vv Mode-8: ``np.testing.*`` / :func:`require` only (fire under ``-O``).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.transport.mesh.material_xs_field import MaterialXSField
from orpheus.transport.operators.isotropic_scattering import (
    IsotropicScattering,
    IsotropicN2N,
)

pytestmark = pytest.mark.foundation


def require(condition: bool, message: str) -> None:
    if not condition:
        pytest.fail(message)


# ASYMMETRIC P0 + non-vacuous, per-material-distinct (n,2n) (so a group-axis flip
# OR a dropped n2n channel is detectable — L27 / #269).
_SIGS0_A = np.array([[0.38, 0.10], [0.05, 0.90]])
_SIGS0_B = np.array([[0.55, 0.03], [0.12, 0.40]])
_SIG2_A = np.array([[0.00, 0.03], [0.01, 0.00]])
_SIG2_B = np.array([[0.02, 0.00], [0.00, 0.04]])


def _mat_xs(nx=6):
    """Synthetic two-material XS field over ``nx`` cells (ny=1)."""
    half = nx // 2
    ix = np.arange(nx)
    iy = np.zeros(nx, dtype=int)
    cells = {0: (ix[:half], iy[:half]), 1: (ix[half:], iy[half:])}
    return MaterialXSField._synthetic_for_tests(
        sig_s={0: [_SIGS0_A], 1: [_SIGS0_B]},
        sig2={0: _SIG2_A, 1: _SIG2_B},
        cells_by_mat=cells, ng=2, nx=nx, ny=1,
    )


def _phi(nx=6, spatial_moments=0, seed=1):
    rng = np.random.default_rng(seed)
    shape = (2, nx, 1) if spatial_moments == 0 else (2, nx, 1, spatial_moments)
    return rng.uniform(0.1, 1.0, size=shape)


# ═══════════════════════════════════════════════════════════════════════
# Gate #1 — apply ≡ legacy fast-path (the bit-identity-inheritance anchor).
# ═══════════════════════════════════════════════════════════════════════


class TestApplyEqualsFastPath:
    @pytest.mark.parametrize("sm", [0, 4], ids=["scalar", "LD-2^d=4"])
    def test_p0_apply_is_apply_p0_in_scatter(self, sm):
        mat = _mat_xs()
        phi = _phi(spatial_moments=sm)
        ref = np.zeros_like(phi)
        mat.apply_p0_in_scatter(ref, phi)
        np.testing.assert_array_equal(
            IsotropicScattering(mat).apply(phi), ref,
            err_msg="IsotropicScattering.apply must route through apply_p0_in_scatter (0-ULP).",
        )

    @pytest.mark.parametrize("sm", [0, 4], ids=["scalar", "LD-2^d=4"])
    def test_n2n_apply_is_apply_n2n(self, sm):
        mat = _mat_xs()
        phi = _phi(spatial_moments=sm)
        ref = np.zeros_like(phi)
        mat.apply_n2n(ref, phi)
        np.testing.assert_array_equal(
            IsotropicN2N(mat).apply(phi), ref,
            err_msg="IsotropicN2N.apply must route through apply_n2n (0-ULP).",
        )

    def test_sum_equals_inplace_iso_then_n2n(self):
        r"""``(IsoScat + IsoN2N).apply`` ≡ the in-place ``add_iso`` then ``add_n2n``
        accumulation (``0 + A`` is exact, so the OperatorSum reorders nothing)."""
        mat = _mat_xs()
        phi = _phi()
        combined = IsotropicScattering(mat).apply(phi) + IsotropicN2N(mat).apply(phi)
        inplace = np.zeros_like(phi)
        mat.apply_p0_in_scatter(inplace, phi)
        mat.apply_n2n(inplace, phi)
        np.testing.assert_array_equal(
            combined, inplace,
            err_msg="The OperatorSum order (P0 + n2n) must match the in-place accumulation.",
        )


# ═══════════════════════════════════════════════════════════════════════
# Gate #2 / #9 — the transpose: reciprocity + structurally-independent dense ref.
# ═══════════════════════════════════════════════════════════════════════


class TestTranspose:
    @pytest.mark.parametrize(
        "factory", [IsotropicScattering, IsotropicN2N], ids=["iso", "n2n"]
    )
    @pytest.mark.parametrize("sm", [0, 4], ids=["scalar", "LD-2^d=4"])
    def test_euclidean_reciprocity(self, factory, sm):
        r"""``⟨Kφ, χ⟩ = ⟨φ, Kᵀχ⟩`` (full Euclidean contraction — K_iso has no
        angular axis to telescope, so the full inner product is honest)."""
        mat = _mat_xs()
        op = factory(mat)
        phi, chi = _phi(spatial_moments=sm, seed=1), _phi(spatial_moments=sm, seed=2)
        lhs = float((op.apply(phi) * chi).sum())
        rhs = float((phi * op.apply_transpose(chi)).sum())
        np.testing.assert_allclose(
            lhs, rhs, rtol=1e-12,
            err_msg=f"{factory.__name__} reciprocity ⟨Kφ,χ⟩=⟨φ,Kᵀχ⟩ violated.",
        )

    def test_iso_transpose_matches_dense_per_material(self):
        r"""STRUCTURALLY-INDEPENDENT: ``IsotropicScattering.apply_transpose`` applies
        :math:`\Sigma_{s,0}` (un-transposed) — built here by an explicit per-(material,
        cell) ``sig @ vec`` loop. A wrong group axis in the production verb disagrees.
        """
        mat = _mat_xs()
        chi = _phi(seed=3)
        got = IsotropicScattering(mat).apply_transpose(chi)
        want = np.zeros_like(chi)
        for mid, (ix, iy), sig in [(0, mat.cells_by_material[0], _SIGS0_A),
                                   (1, mat.cells_by_material[1], _SIGS0_B)]:
            for c in range(len(ix)):
                # forward out_g = Σ_f sig[f,g] in_f  ⇒  transpose out_f = Σ_g sig[f,g] in_g
                want[:, ix[c], iy[c]] = sig @ chi[:, ix[c], iy[c]]
        np.testing.assert_allclose(
            got, want, rtol=1e-13, atol=0,
            err_msg="IsotropicScattering.apply_transpose group axis disagrees with the dense Σ_s0 @ vec.",
        )

    def test_n2n_distinct_and_factor_two(self):
        r"""(n,2n) is a DISTINCT channel: ``IsotropicN2N.apply`` = ``2·Σ_2nᵀ`` (NOT
        the scattering matrix), so it differs from ``IsotropicScattering.apply`` and
        carries the factor 2 (drop-factor-2 / fold-into-scatter would red)."""
        mat = _mat_xs()
        phi = _phi(seed=5)
        iso, n2n = IsotropicScattering(mat).apply(phi), IsotropicN2N(mat).apply(phi)
        require(not np.allclose(iso, n2n),
                "IsotropicN2N must be a distinct channel from IsotropicScattering.")
        # factor 2: n2n.apply == 2 × (Σ_2nᵀ @ φ) per material
        want = np.zeros_like(phi)
        for (ix, iy), sig in [(mat.cells_by_material[0], _SIG2_A),
                              (mat.cells_by_material[1], _SIG2_B)]:
            for c in range(len(ix)):
                want[:, ix[c], iy[c]] = 2.0 * (sig.T @ phi[:, ix[c], iy[c]])
        np.testing.assert_allclose(n2n, want, rtol=1e-13, atol=0,
                                   err_msg="IsotropicN2N must apply 2·Σ_2nᵀ.")


# ═══════════════════════════════════════════════════════════════════════
# Gate #6 — as_dense ≡ apply (the two consumption modes agree).
# ═══════════════════════════════════════════════════════════════════════


class TestDensePerMaterial:
    @pytest.mark.parametrize(
        "factory,xs", [(IsotropicScattering, {0: _SIGS0_A, 1: _SIGS0_B})],
        ids=["iso"],
    )
    def test_iso_dense_is_operator_matrix(self, factory, xs):
        mat = _mat_xs()
        dense = factory(mat).dense_per_material()
        for mid, sig in xs.items():
            np.testing.assert_array_equal(
                dense[mid], sig.T,
                err_msg="IsotropicScattering.dense_per_material must be Σ_s0ᵀ (M @ φ == apply).",
            )

    def test_n2n_dense_is_two_sig2_transpose(self):
        mat = _mat_xs()
        dense = IsotropicN2N(mat).dense_per_material()
        for mid, sig in {0: _SIG2_A, 1: _SIG2_B}.items():
            np.testing.assert_array_equal(dense[mid], 2.0 * sig.T)

    def test_dense_matvec_equals_apply(self):
        r"""``M @ φ_cell`` (the LHS-fold consumption mode) ≡ ``apply(φ)_cell``."""
        mat = _mat_xs()
        phi = _phi(seed=7)
        op = IsotropicScattering(mat)
        dense, got = op.dense_per_material(), op.apply(phi)
        mat_of = {0: mat.cells_by_material[0], 1: mat.cells_by_material[1]}
        for mid, (ix, iy) in mat_of.items():
            for c in range(len(ix)):
                np.testing.assert_allclose(
                    dense[mid] @ phi[:, ix[c], iy[c]], got[:, ix[c], iy[c]],
                    rtol=1e-13, atol=0,
                    err_msg="dense_per_material M @ φ must equal apply(φ).",
                )


# ═══════════════════════════════════════════════════════════════════════
# Gate #11 — per-axis structural predicates.
# ═══════════════════════════════════════════════════════════════════════


class TestPredicates:
    @pytest.mark.parametrize("factory", [IsotropicScattering, IsotropicN2N])
    def test_apply_and_transpose_not_solve(self, factory):
        op = factory(_mat_xs())
        require(callable(getattr(op, "apply", None)), "must expose apply.")
        require(op.is_adjointable,
                "must advertise the adjoint axis (campaign #276).")
        require(not op.is_invertible,
                "isotropic energy operator must NOT be invertible (singular group-transfer).")
        require(not hasattr(op, "inverse"),
                "isotropic energy operator is STRUCTURALLY non-invertible — "
                "it must not declare an inverse method at all.")
