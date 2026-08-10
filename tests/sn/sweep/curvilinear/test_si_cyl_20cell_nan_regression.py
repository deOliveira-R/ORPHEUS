"""Regression catcher: SI cylindrical 1-D pole-cell algebraic resonance
triggers ``ordinate_scan`` NaN.

Bug class: ``orpheus/sn/sweep/scan.py`` Blelloch-form
``cumprod_a * (psi_0 + cumsum(b / cumprod_a))`` produces NaN whenever
the per-cell attenuation chain ``a`` contains an exact zero.  At the
cylindrical pole cell ``A_down = 0`` (the inner radial face at r=0
has zero area); this makes ``a = 2|μ|·A_total / (dA_w·c_out + Σ_t·V)
− 1`` vanish at every (μ, dr, Σ_t) triple that satisfies the
algebraic identity ``2|μ|·A_total = dA_w·c_out + Σ_t·V``.  At the
then-canonical resonance point (μ_x = 1/√20 — an LS8 node of the
pre-#337 project-convention seed, since retired: no shipped rule
carries that node today; dr = 0.1, Σ_t = 1.0 of mixture A group 1)
``a`` was bit-EXACTLY zero, and the cumprod collapsed for the whole
tail of the chain.

REACHABILITY TODAY (Q5.6.3, measured 2026-08-08): on the folded
family ``SNMesh(CYLINDRICAL)`` admits, the pole-cell resonance is
UNREACHABLE at physical cross sections — the per-ordinate resonant
``Σ_t*`` solves ≤ 0 for every inward ordinate at every probed
(order, mesh).  The reachability tripwire lives in
``tests/sn/sweep/core/test_ordinate_scan_reset.py`` (it reds if the
6.4 absorber retirement changes ``c_out`` enough to bring the
resonance back); the ``a ∋ 0`` scan-form contract itself stays
pinned QUADRATURE-FREE by the last test in this file.  The solver
tests here keep their finite/k_inf pins on the admitted family.

CLOSED-FORM v EXPLICIT LOOP.  The explicit recurrence
``ψ[i+1] = a·ψ[i] + b[i]`` is well-defined at a=0; it gives ``ψ = b``
which is the limit ``Σ_t → ∞`` flux (fully-absorbed cell), a
physically meaningful solution.  The Blelloch *algorithm* trades
numerical stability for vectorisation.  When the chain contains an
exact zero, the algorithm fails by an algorithmic-stability
violation, not a math error.

WHY THIS WAS INVISIBLE TO EXISTING TESTS.

* ``tests/sn/sweep/core/test_ordinate_scan.py::test_ordinate_scan_small_attenuation``
  uses ``a ∈ [0.05, 0.2]`` (positive, bounded away from 0).
* ``tests/sn/sweep/core/test_ordinate_scan.py::test_ordinate_scan_zero_attenuation``
  tests ``a = 1`` everywhere (pure accumulation), the OPPOSITE
  regime.
* No test covers ``a ∋ 0`` (the pole-cell pathology).
* The integration tests on cylindrical SI all happen to use mesh
  refinements (n_cells = 10, 40, 80, 160) that do NOT land on the
  resonance.  n_cells = 20 with thickness = 2.0 is unique to the
  L1 standoff benchmark suite.

CATCHES: future ERR-NNN ``si-cyl-pole-cell-resonance``.

FIX SITE: ``orpheus/sn/sweep/scan.py``.  The Blelloch closed form
was replaced by a numerically-stable alternative — option (a) landed:
(a) the pair-monoid prefix scan ``(α, β) ⊕ (α', β') = (α·α', α'·β + β')``
which has no division (dispatched when the closed form is non-finite);
(rejected: (b) per-segment fall-back to the explicit recurrence at
cells where ``|a| < ε``; (c) a Brent blocked scan with bounded
condition number per block).  All three preserve the joint
vectorisation across (K, ng); none divide by cumprod_a.

The architectural fix point is the C5 retirement of the
``angular_flux.py`` legacy buffer; consider landing the scan fix
alongside.
"""
import numpy as np
import pytest

from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solver import solve_sn
from orpheus.derivations.common.xs_library import get_mixture

# SI cylindrical pole-cell NaN regression catcher: ``regression`` flags
# the frozen-behaviour drift gate; ``foundation`` gives it a V&V-level
# so the audit does not report it as an orphan. Both compose. (Was a
# V&V orphan before the taxonomy reorg forced a marker.)
pytestmark = [pytest.mark.regression, pytest.mark.foundation]


@pytest.fixture
def homog_cyl_2g_thick2_n20():
    """The bug report's spatial configuration (2 cm, 20 cells, mixture A
    2g, reflective).  The quadrature was LS8 in the report and rode the
    Q5.6.3 migration onto the admitted folded family; the resonance
    node itself was an artifact of the retired pre-#337 LS8 seed (see
    the module docstring), so the pins here are finite-ness and k_inf,
    not resonance-landing."""
    fuel = get_mixture('A', '2g')
    geom = StructuredGeometry(
        geometry='CYL',
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(BC.reflective,),
    )
    mesh = Mesh1D.from_geometry(
        geom, region_meshes=(RegionMesh(n_cells=20),),
    )
    quad = Quadrature.folded_product(n_mu=8, n_phi=16)
    return fuel, mesh, quad


@pytest.mark.filterwarnings(
    "ignore::orpheus.numerics.convergence.ConvergenceWarning"
)
def test_si_returns_finite_keff(homog_cyl_2g_thick2_n20):
    """SI must return a finite k_eff at the bug report's configuration.

    Pre-fix this test FAILED with keff=NaN — pinning the bug class.
    Post-fix (the division-free pair-monoid backend, dispatched when
    the closed form is non-finite) it passes with k_eff = k_inf =
    1.875.

    ⭐ **The truncation is the FIXTURE, and it is declared** (#340 R2/R4).
    ``max_inner=3`` is deliberate: the NaN this row pins appears in the
    FIRST inner iteration, so a converged solve would only make the test
    slower without making it sharper.  The solve therefore exits truncated
    and says so — which is correct, and is why the suppression above names
    the ONE category this row expects rather than silencing the module.

    ⛔ Until 2026-08-10 that suppression was a module-level, category-less
    ``warnings.filterwarnings('ignore')``, evaluated at IMPORT — so it
    silenced every warning of every category, and `[M]` the only thing it
    was actually hiding was this ``ConvergenceWarning``
    (``inner(source-iteration) hit max_inner=3 … set max_inner=1093``).
    A global silencer installed to hide one expected line is how the next
    UNexpected one goes unnoticed; `[M]` removing it left all 4 rows green.
    """
    fuel, mesh, quad = homog_cyl_2g_thick2_n20
    # NaN appears in the FIRST inner iteration — small caps suffice.
    res = solve_sn(
        materials={0: fuel}, mesh=mesh, quadrature=quad,
        inner_solver='source_iteration',
        keff_tol=1e-10, flux_tol=1e-9,
        max_outer=3, max_inner=3,
    )
    assert np.isfinite(res.keff), (
        f"SI returned non-finite k_eff: {res.keff}.  Bug class: "
        f"ordinate_scan Blelloch-form NaN at pole-cell algebraic "
        f"resonance (ERR-054; its denormal-underflow sibling is "
        f"ERR-057).  The scan-form contract and the reachability "
        f"tripwire live in "
        f"tests/sn/sweep/core/test_ordinate_scan_reset.py; the "
        f"division-free backend is orpheus/sn/sweep/scan.py."
    )


@pytest.mark.slow
def test_si_agrees_with_kinf_at_resonance(homog_cyl_2g_thick2_n20):
    """SI must converge to k_inf = νΣ_f/Σ_a = 1.875 for homogeneous
    reflective configurations (irrespective of mesh discretisation).

    Pre-fix: this test FAILS (NaN ≠ 1.875).
    Post-fix: this test passes within the eigenvalue tolerance.
    """
    fuel, mesh, quad = homog_cyl_2g_thick2_n20
    # Post-fix this test pins SI ↔ k_inf convergence on the bug
    # report's spatial configuration.  Pre-fix it FAILED with NaN
    # before convergence was checked.
    res = solve_sn(
        materials={0: fuel}, mesh=mesh, quadrature=quad,
        inner_solver='source_iteration',
        keff_tol=1e-10, flux_tol=1e-9,
        max_outer=200, max_inner=200,
    )
    assert np.isfinite(res.keff), f"keff is non-finite: {res.keff}"
    assert abs(res.keff - 1.875) < 1e-5, (
        f"SI k_eff = {res.keff}, expected k_inf = 1.875 "
        f"(homogeneous reflective ⇒ k = νΣ_f/Σ_a)"
    )


def test_krylov_unaffected(homog_cyl_2g_thick2_n20):
    """Krylov inner solver bypasses the buggy ``ordinate_scan`` path.

    This test exists to anchor the structural-independence claim:
    the bug is in the SI sweep code path, not in the operator
    construction.  Krylov consumes the ``transport_operator_matvec``
    path which does NOT call ``ordinate_scan``.

    If this test FAILS, the matvec path has been refactored to also
    consume ``ordinate_scan`` — at which point the bug would be visible
    in Krylov too, contradicting the diagnosis.
    """
    fuel, mesh, quad = homog_cyl_2g_thick2_n20
    res = solve_sn(
        materials={0: fuel}, mesh=mesh, quadrature=quad,
        inner_solver='krylov',
        keff_tol=1e-10, flux_tol=1e-9,
        max_outer=100, max_inner=100,
    )
    assert np.isfinite(res.keff)
    assert abs(res.keff - 1.875) < 1e-6


def test_ordinate_scan_at_a_zero_returns_finite_via_loop():
    """Pure numerical scan-form contract test.

    Independent of solver geometry: ``ordinate_scan`` must remain
    finite when the chain ``a`` contains exact zeros.  The
    mathematically-equivalent explicit recurrence does; the current
    Blelloch closed form does NOT.

    This is the "structural" test of the fix.  The fix (the
    division-free pair-monoid backend, dispatched when the closed
    form is non-finite) makes this assertion pass; the
    SI-cylinder-20cell test above passes by structural consequence.
    """
    from orpheus.sn.sweep.scan import ordinate_scan
    # Chain with an exact zero at the tail (mirrors the failing
    # cylindrical pole-cell config).
    a = np.array([0.7965626, 0.79170007, 0.78645589, 0.78077641,
                  0.77459667, 0.7678371, 0.76039895, 0.75215787,
                  0.74295459, 0.73258108, 0.72075922, 0.70710678,
                  0.69108049, 0.67187516, 0.64823152, 0.61803399,
                  0.57735027, 0.51763809, 0.41421356, 0.0])
    b = np.full_like(a, 0.5)
    psi_0 = np.array([0.3])

    out = ordinate_scan(a[:, None], b[:, None], psi_0)
    assert np.all(np.isfinite(out)), (
        f"ordinate_scan returned non-finite values when a contains 0: "
        f"out = {out.ravel()}.  Pre-fix: this is expected (the bug). "
        f"Post-fix: out[-1] should equal b[-1] = 0.5 by the explicit "
        f"recurrence."
    )
    # Post-fix assertion (after Blelloch is replaced with a stable form):
    assert abs(out.ravel()[-1] - 0.5) < 1e-12, (
        f"out[-1] = {out.ravel()[-1]}, expected b[-1] = 0.5"
    )
