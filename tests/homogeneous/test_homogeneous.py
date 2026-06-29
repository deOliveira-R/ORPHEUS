"""Verify the infinite-medium eigenvalue solver against SymPy analytical solutions."""

import pytest

from orpheus.derivations import get
from orpheus.homogeneous.solver import solve_homogeneous_infinite

# File-level verifies marker: every test in this file exercises the
# homogeneous eigenvalue chain end-to-end by asserting k_inf matches
# the SymPy-derived analytical reference to 1e-12. That tolerance is
# tight enough to pin every step of the derivation — if any of the
# labelled equations below were implemented incorrectly the k mismatch
# would be far larger than 1e-12, so a passing test is equation-level
# (L1) verification for every link in the chain.
#
# Declared explicitly here (rather than inherited from
# VerificationCase) so the Nexus AST pass picks it up via decorator
# parsing and writes TESTS edges.
#
# The 2G labels (two-group-*) and the power-iteration step labels
# (fission-source, fixed-source-solve, keff-update) are all exercised
# by the homo_2eg / homo_4eg cases — the analytical k_inf is derived
# symbolically via exactly those equations, so a solver k that matches
# to 1e-12 implies every link in the chain is correct. The absorption-
# xs label is the derived property used inside keff-update.
pytestmark = [pytest.mark.l1, pytest.mark.verifies(
    "one-group-kinf",
    "inf-hom-balance",
    "matrix-eigenvalue",
    "removal-matrix",
    "fission-matrix",
    "mg-balance",
    # B.1 additions (issue #87): the full 2G analytical chain and the
    # power-iteration step labels, all verified end-to-end by the
    # homo_2eg and homo_4eg parametric cases.
    "two-group-A",
    "two-group-F",
    "two-group-Ainv",
    "two-group-M",
    "two-group-charpoly",
    "two-group-roots",
    "fission-source",
    "fixed-source-solve",
    "keff-update",
    "absorption-xs",
)]


@pytest.mark.parametrize("case_name", [
    "homo_1eg",
    "homo_2eg",
    "homo_4eg",
    # Non-trivial, asymmetric (n,2n): de-vacuums the n2n-in-A convention.
    # Every other case has Sig2=0, so the 2·Σ₂ᵀ loss term is never
    # exercised; here it moves k_inf ~0.6 (a drop/double-count reds).
    "homo_2eg_n2n",
])
def test_kinf_exact(case_name):
    """Eigenvalue must match analytical solution to machine precision.

    The refounded solver assembles A = C − K_iso from the transport
    operators (a different FP reduction tree than the oracle's fused
    ``(Σ_s + 2Σ_2)ᵀ``), so the tolerance is now principled-equivalence
    (FP-non-associativity, ~1 ULP), not bit-identity — still ≪ 1e-12.
    """
    case = get(case_name)
    mix = next(iter(case.materials.values()))
    result = solve_homogeneous_infinite(mix)
    assert abs(result.k_inf - case.k_inf) < 1e-12, (
        f"k_inf mismatch: solver={result.k_inf:.10f} "
        f"analytical={case.k_inf:.10f}"
    )


@pytest.mark.verifies("normalisation")
def test_post_solve_production_rate_is_100():
    """L1: post-convergence flux is normalised to 100 n/cm^3/s production.

    After :func:`solve_homogeneous_infinite` solves, the flux is rescaled
    so the **fission** production rate

    .. math::

       \\nu\\Sigma_f \\cdot \\boldsymbol{\\phi} = 100

    (see Eq. ``normalisation`` in docs/theory/homogeneous.rst).
    Production is :math:`\\nu\\Sigma_f` only — the (n,2n) reaction is a
    loss-side transfer folded into the loss matrix as
    :math:`2\\Sigma_2^T`, NOT a production channel (it does NOT enter
    the numerator).  The ``homo_2eg_n2n`` case (non-zero, asymmetric
    :math:`\\Sigma_2`) makes this non-vacuous: under the retired
    ``(\\Sigma_p + 2\\cdot\\text{colsum}(\\Sigma_2))`` formula the
    production would not equal 100.  The 1G case is a degenerate
    one-scalar normalisation a bug could accidentally satisfy, so it is
    excluded.
    """
    for case_name in ("homo_2eg", "homo_4eg", "homo_2eg_n2n"):
        case = get(case_name)
        mix = next(iter(case.materials.values()))
        result = solve_homogeneous_infinite(mix)

        # Production = νΣ_f @ φ  (fission only; n2n lives in A, not F).
        production = mix.SigP @ result.flux

        assert abs(production - 100.0) < 1e-9, (
            f"{case_name}: production rate = {production:.6e}, "
            f"expected 100.0 (normalisation constraint)"
        )
