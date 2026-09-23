r"""Pin the adjoint-weighted collapse taxonomy — campaign P6 (#281) B1.

:mod:`orpheus.derivations.common.homogenization` is the algebra of record
for the eigenvalue-consistent (adjoint-weighted) homogenization /
condensation rules: which per-channel collapse rule zeroes the first-order
XS-collapse worth (spatial axis, T1–T5) and which condensation convention
reproduces the fine :math:`k` exactly (energy axis, T6).  The production
``adjoint=`` arm (B2) implements exactly what these theorems prove, and the
C-gates cite them — so the theorems themselves must stay green under the
tree: this file runs each derivation as a test.

Every theorem verifies its identity in EXACT (rational/symbolic) arithmetic
inside the module and raises
:class:`~orpheus.derivations.common.homogenization.DerivationError` on any
failure — there is no tolerance anywhere in this file.  One test per
theorem so a regression localizes to the broken identity.

vv Mode-8 note: the module uses explicit raises (never bare ``assert``), so
the theorems fire under the canonical ``python -O`` invocation; the tests
below simply propagate.
"""
from __future__ import annotations

import pytest

from orpheus.derivations.common import homogenization as hom

pytestmark = pytest.mark.foundation


class TestAdjointWeightedCollapseTheorems:
    """The T0–T6 theorem battery (each raises ``DerivationError`` on failure)."""

    def test_t0_first_order_k_shift_formula(self):
        r"""δμ = ⟨x*,(δA − μδF)x⟩/⟨x*,Fx⟩ — exact ε-series check on a rational pencil."""
        hom.derive_first_order_k_shift()

    def test_t1_vector_channels_bilinear_pair_weight(self):
        r"""φ*⊙φ zeroes the vector-channel worth, uniquely; forward rule does not."""
        hom.derive_vector_channel_rule()

    def test_t1b_collision_channel_angular_pairing(self):
        r"""The exact Σt rule weights by ρ = Σ_n w ψ*ψ (unique); the scalar pair
        is its isotropic limit and fails on anisotropic shapes (user-ruled
        implemented — option 2 at the P6 open)."""
        hom.derive_angular_sigma_t_rule()

    def test_t2_matrix_channels_per_pair_weight(self):
        r"""Per-pair φ*_g·φ_g' zeroes every (g',g) worth; the source-product
        broadcast (what the per-(cell, group) plumbing would produce) does not."""
        hom.derive_matrix_channel_rule()

    def test_t3_fission_dyad_mixed_fold_factored_rule(self):
        r"""The mixed-fold νΣf rule zeroes the TOTAL fission worth for ANY simplex
        χ_R; the canonical χ_R is convex; flat-φ* degenerates to the forward pair."""
        hom.derive_fission_factored_rule()

    def test_t4_balance_tradeoff(self):
        r"""Worth-exactness ⟺ broken Σt balance (φ* non-flat); both restored flat."""
        hom.derive_balance_tradeoff()

    def test_t5_forward_rule_first_order_error(self):
        r"""The forward collapse worth is O(ε) in the adjoint tilt — the C2 basis."""
        hom.derive_forward_rule_first_order_error()

    def test_t6_energy_condensation_exactness(self):
        r"""Row-consistent bilinear condensation: k_C == k_fine EXACTLY at the true
        spectra; O(ε²) under spectrum error vs the forward rules' O(ε); mixed
        conventions are not exact."""
        hom.derive_energy_condensation_exactness()
