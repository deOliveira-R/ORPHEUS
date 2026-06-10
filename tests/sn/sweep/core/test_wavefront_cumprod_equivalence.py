r"""wavefront(1-D) ≡ cumprod(1-D) — the spine-vs-optimization oracle (C3).

THIS FILE IS A TEST-ARCHITECT STUB — UNCOMMITTED, FOR MAIN-AGENT REVIEW.

Architecture (the EXPANDED C3 scope — settled, not re-litigated here)
====================================================================
The dimension-generic wavefront DAG walk
(``SweepDependencyGraph.from_cartesian(shape)`` + the d-generic
``DiamondDifference.cell_kernel_batch`` summing ``denom = σ_t +
Σ_axis s_axis``) is the SPINE: ONE code path for d = 1, 2, 3, and the
VERIFICATION ORACLE.

The existing 1-D cumprod parallel-prefix scan (``ordinate_scan`` /
``_sweep_1d_unified``, the Blelloch 1990 §1.5 closed form
``cumprod_a · (psi_0 + cumsum(b / cumprod_a))``) is reframed as the
d = 1 OPTIMIZATION — default-selected for d=1 production (zero perf
regression). This file proves the optimization computes the SAME 1-D
transport sweep as the spine.

Bit-identical or principled-equivalence?  → PRINCIPLED-EQUIVALENCE.
==================================================================
The cumprod scan and the wavefront forward-substitution walk compute
the SAME affine recurrence per ordinate

    psi_out[i] = 2·psi_bar − psi_in[i],
    psi_bar    = (Q + s·psi_in[i]) / (σ_t + s)

but in DIFFERENT floating-point ASSOCIATION ORDERS:

* cumprod path: ``cumprod_a · (psi_0 + cumsum(b / cumprod_a))`` —
  a parallel-prefix reduction TREE over all nx cells at once (the
  whole chain's products/sums are formed before the final
  multiply).
* wavefront spine: cell-by-cell forward substitution along the DAG
  levels — each cell's ``psi_in`` is the previous cell's ``psi_out``,
  evaluated sequentially.

Because IEEE-754 addition/multiplication is NOT associative, these two
real-arithmetic-identical computations disagree at ULP level. Per
``vv-principles`` §"Bit-identity vs principled-equivalence", this is a
principled-equivalence-at-nulp case — NOT a bug, an FP reduction-order
difference — and ALL THREE criteria hold:

  (1) Principled at every step. Both paths' intermediates are named
      reactor-physics quantities: the per-cell downstream face flux
      ``psi_out[i]`` and the cell average ``psi_bar``. Neither is an
      anonymous reduction.
  (2) Verified against a structurally-independent reference. The
      anchor is NOT old-vs-new ULP distance (necessary, never
      sufficient) — it is the ANALYTICAL closed form
      ``k_inf = λ_max(A⁻¹F) = 1.875`` for mixture A 2g (reflective
      infinite-medium eigenvalue), computed independently by
      ``kinf_and_spectrum_homogeneous`` (transfer-matrix), see
      ``test_cumprod_path_hits_analytical_kinf`` below. The cumprod
      path — the path the spine is validated against — is ITSELF
      anchored to closed-form ground truth, so "spine ≡ cumprod"
      transitively anchors the spine.
  (3) Drift is FP-non-associativity, dimensionally explainable: a
      single-step computation (one sweep, no outer iteration) is
      bounded by ``(reduction depth) × ULP`` ≈ ``nx × ε`` per cell.

The chosen nulp bound is derived in ``_NULP_BOUND`` below from the
reduction depth (the chain length nx) with a safety factor.

Cardinal Rule 6 + vv anti-patterns #3/#4 — what this test STRESSES
==================================================================
A 1-group flat-flux 1-D test proves NOTHING (H1 1-group degeneracy +
H2 homogeneous-flat degeneracy: a flat ψ nulls the redistribution and
makes the recurrence trivial). The equivalence test therefore uses:

* ≥2 GROUPS (mixture A 2g — asymmetric downscatter SigS[0,1]=0.1) so a
  group-coupling / x↔y-analogue convention drift would surface.
* NON-FLAT flux: a heterogeneous Σ_t profile (per-cell random, strictly
  positive) AND a non-uniform source AND a non-zero boundary inflow,
  so ``psi_in ≠ psi_out`` at every cell (the redistribution term is
  driven OUT of cancellation — the whole point of the recurrence).
* A heterogeneous-region slab for the converged anchor, so the
  per-region flux is non-flat (the eigenvalue is still k_inf because
  the medium is infinite/reflective, but the flux SHAPE is stressed).

``-O`` discipline (vv Mode 8): every assertion is a ``np.testing.*`` /
``pytest.fail`` FUNCTION CALL — fires under the canonical ``python -O``.
The per-sweep equivalence + the analytical anchor run under ``-O`` (no
bare assert). The perf measurement is ``@pytest.mark.slow`` + a LOOSE
final-gate tripwire only.

API ASSUMPTIONS (confirm spelling at review — only the adapters change)
=======================================================================
* ``transport_sweep(AngularSourceSink, sig_t (ng,nx), SNMesh,
  BoundaryFlux)`` — the cumprod path (CONFIRMED live; dispatches to
  ``_sweep_1d_unified`` for a reduced 1-D mesh). Returns
  ``(angular_flux (N,ng,nx,1), scalar_flux (ng,nx,1))``.
* ``SweepDependencyGraph.from_cartesian((nx,), label=<1-sign-tuple>)``
  — the C3 d-generic spine builder at d = 1 (the carve target; NOT yet
  landed). The wavefront-1D driver is isolated in
  ``_wavefront_1d_sweep`` (the ONLY adapter that changes if the final
  spine spelling differs). It builds the d=1 DAG, walks it with the
  d-generic ``cell_kernel_batch`` (single-axis ``s``), seeds the BC
  inflow, and returns the SAME ``(angular, scalar)`` layout as
  ``transport_sweep``.

Until the C3 spine lands, ``_wavefront_1d_sweep`` raises
``_SpineNotLanded`` and the equivalence test is marked
``xfail(strict=False)`` with a reason pointing at the unlocking API
(``feedback_vv_tagging`` idiom). The analytical-anchor test + the perf
measurement's cumprod arm run TODAY (they touch only the live cumprod
path), so the file is not inert.

``foundation`` — software invariants (spine ≡ optimization; cumprod ≡
analytical limit), no theory ``:label:``, no ``verifies(...)``.
"""
from __future__ import annotations

import time

import numpy as np
import pytest

from orpheus.derivations.common.eigenvalue import kinf_and_spectrum_homogeneous
from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.geometry import SNMesh
from orpheus.sn.solver import solve_sn
from orpheus.sn.spatial.diamond import DiamondDifference
from orpheus.sn.sweep import transport_sweep
from orpheus.sn.sweep_graph import SweepDependencyGraph
from orpheus.transport.fields.boundary_flux import BoundaryFlux
from orpheus.transport.source_sinks import AngularSourceSink

pytestmark = pytest.mark.foundation


# ─── nulp bound: reduction-depth × ULP × safety ─────────────────────
#
# Per vv-principles §"Bit-identity vs principled-equivalence" criterion
# (3): a single-step (one-sweep) computation's FP drift is bounded by
# (reduction depth) × ULP. The reduction depth of the cumprod chain is
# the cell count nx (the cumprod/cumsum span the whole chain); the
# wavefront does nx sequential forward substitutions. With the test's
# nx ≤ 16, a depth ≈ 16 plus a generous safety factor (×8, covering the
# 2-group coupling + the source/BC affine terms) gives:
#
#     _NULP_BOUND = 8 × 16 = 128
#
# This is a STARTING bound the carve confirms empirically (measure the
# actual ULP distance once the spine lands; tighten if the measured max
# is ≪ 128, loosen + RE-JUSTIFY here if it exceeds — never silently).
_NULP_BOUND = 128

_NX_EQUIV = 12           # chain length for the per-sweep equivalence
_LENGTH = 6.0            # slab thickness (cm)


class _SpineNotLanded(RuntimeError):
    """The C3 d-generic wavefront spine builder is not yet available."""


# ─── adapters (UPDATE ONLY THESE if the final C3 spelling differs) ───


def _slab_sn_mesh(nx: int, *, bc: str, ng_key: str = "2g") -> SNMesh:
    """Heterogeneous-capable slab SNMesh, mixture A, Gauss-Legendre S8."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, _LENGTH, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC(bc),
        bc_right=BC(bc),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    return SNMesh(mesh, quad, {0: get_mixture("A", ng_key)})


def _cumprod_1d_sweep(Q, sig_t, sn_mesh, boundary_flux):
    """Drive the cumprod optimization (CONFIRMED live API).

    ``transport_sweep`` dispatches a reduced 1-D mesh to
    ``_sweep_1d_unified`` → ``ordinate_scan`` (the Blelloch cumprod
    closed form). Returns ``(angular (N,ng,nx,1), scalar (ng,nx,1))``.
    """
    return transport_sweep(Q, sig_t, sn_mesh, boundary_flux)


def _wavefront_1d_sweep(Q_arr, sig_t, sn_mesh, bc_inflow):
    """Drive the d-generic wavefront SPINE at d = 1 (ASSUMED C3 API).

    This is the ONLY adapter that changes when the C3 spine lands. It
    must compute ONE forward-substitution sweep of ``(L+C)⁻¹q`` on the
    1-D chain via ``from_cartesian((nx,))`` + the d-generic
    ``cell_kernel_batch`` (single streaming axis), seeded with the same
    BC inflow as the cumprod path, and return the SAME
    ``(angular (N,ng,nx,1), scalar (ng,nx,1))`` layout.

    Wiring sketch the carve will realise (per-octant, per the d=2
    ``_sweep_2d_wavefront`` precedent):

        graph_plus  = SweepDependencyGraph.from_cartesian((nx,), label=(+1,))
        graph_minus = SweepDependencyGraph.from_cartesian((nx,), label=(-1,))
        for each octant (sign of μ_n):
            s = 2|μ_n| · A_face / V          # single-axis streaming coeff
            seed psi_in at the upwind domain face from ``bc_inflow``
            graph.apply(cell_update=DiamondDifference(), ..., s, Q, sig_t)
            accumulate scalar_flux += w_n · psi_bar
        return angular, scalar

    Raises ``_SpineNotLanded`` until ``from_cartesian`` exists at d = 1.
    """
    if not hasattr(SweepDependencyGraph, "from_cartesian"):
        raise _SpineNotLanded(
            "SweepDependencyGraph.from_cartesian (the C3 d-generic spine "
            "builder) is not yet landed. "
            "Once the carve lands, wire this adapter to the d=1 walk."
        )
    # ── C3-carve wiring goes here (delete the raise above) ──
    raise _SpineNotLanded(
        "from_cartesian present but the d=1 driver is not wired — "
        "complete _wavefront_1d_sweep per the docstring sketch."
    )


def _spine_landed() -> bool:
    """True once the d=1 wavefront spine is callable end-to-end."""
    try:
        _wavefront_1d_sweep(None, None, None, None)
    except _SpineNotLanded:
        return False
    except Exception:
        # Any OTHER exception means the API exists and the adapter is
        # mis-wired — surface it (don't silently skip).
        return True
    return True


# ─── NEW-1: wavefront(1-D) ≡ cumprod(1-D), principled @ nulp ─────────


@pytest.mark.xfail(
    not _spine_landed(),
    reason=(
        "C3 d-generic wavefront spine (SweepDependencyGraph.from_cartesian "
        "at d=1) not yet landed — unlocks when the spine builder + the d=1 "
        "driver in _wavefront_1d_sweep are wired. strict=False so it flips "
        "to xpass automatically on landing."
    ),
    strict=False,
)
@pytest.mark.parametrize("bc", ["vacuum", "reflective"])
def test_wavefront_1d_equals_cumprod_1d_per_sweep(bc):
    """NEW-1: the d=1 wavefront SPINE computes the SAME single 1-D sweep
    as the cumprod OPTIMIZATION — principled-equivalence at nulp.

    What this catches: a convention drift between the spine's d-generic
    single-axis closure and the cumprod recurrence (sign of the
    streaming coeff, the ``/W`` placement, the BC-inflow seed face, the
    group axis). ≥2G + heterogeneous Σ_t + non-uniform source + non-zero
    inflow ⟹ ``psi_in ≠ psi_out`` everywhere (redistribution driven out
    of cancellation; H1/H2 degeneracies avoided).

    Gate: ``assert_array_almost_equal_nulp`` (FP-association difference,
    NOT bit-identity — see module docstring). The ANALYTICAL anchor for
    the converged value lives in
    ``test_cumprod_path_hits_analytical_kinf`` (k_inf = 1.875); this row
    proves spine ≡ cumprod so the spine inherits that anchor.
    """
    nx = _NX_EQUIV
    sn_mesh = _slab_sn_mesh(nx, bc=bc)
    ng = sn_mesh.ng
    N = sn_mesh.quad.N
    rng = np.random.default_rng(20260609)

    # Heterogeneous, strictly-positive Σ_t (per cell, per group) — the
    # non-flat driver. Shape (ng, nx) is the 1-D cumprod-path layout.
    sig_t = rng.uniform(0.3, 3.0, size=(ng, nx))
    # Non-uniform per-ordinate source magnitude (already /W-projected by
    # using a non-trivial iso source through from_isotropic).
    iso = rng.uniform(0.2, 1.5, size=(ng, nx))
    Q = AngularSourceSink.from_isotropic(iso, sn_mesh)
    Q_arr = Q.values

    # Non-zero boundary inflow so the chain seed is non-trivial (a flat
    # zero-inflow vacuum sweep would weaken the BC-seed convention test).
    bf_cumprod = BoundaryFlux.zeros_on(sn_mesh)
    for face in bf_cumprod.layout.faces:
        fv = bf_cumprod.face_view(face)
        fv[...] = rng.uniform(0.0, 1.0, size=fv.shape)
    bc_inflow = BoundaryFlux.from_mesh(bf_cumprod.values.copy(), sn_mesh)

    ang_c, scal_c = _cumprod_1d_sweep(Q, sig_t, sn_mesh, bf_cumprod)
    ang_w, scal_w = _wavefront_1d_sweep(Q_arr, sig_t, sn_mesh, bc_inflow)

    np.testing.assert_array_almost_equal_nulp(
        ang_w, ang_c, nulp=_NULP_BOUND,
    )
    np.testing.assert_array_almost_equal_nulp(
        scal_w, scal_c, nulp=_NULP_BOUND,
    )


# ─── NEW-1 anchor: the cumprod path hits the analytical k_inf ────────


def test_cumprod_path_hits_analytical_kinf():
    """NEW-1 anchor: the cumprod 1-D path converges to the closed-form
    infinite-medium eigenvalue ``k_inf = λ_max(A⁻¹F) = 1.875`` for
    mixture A 2g (reflective slab → infinite medium).

    This is the STRUCTURALLY-INDEPENDENT ground (vv criterion 2) that
    anchors the whole oracle: old-vs-new ULP distance is necessary but
    never sufficient, so the path the spine is validated against must
    itself terminate in a closed-form reference. ``solve_sn`` on a
    reflective slab drives the cumprod sweep; the transfer-matrix
    ``kinf_and_spectrum_homogeneous`` is the independent analytical
    value (different structural angle: eigendecomposition of A⁻¹F, no
    sweep at all).

    ≥2G (H1): k = 1.875 ≠ a 1-group material ratio — the 2g downscatter
    structure is load-bearing in the transfer matrix.
    """
    mix = get_mixture("A", "2g")
    k_ref, _phi = kinf_and_spectrum_homogeneous(
        mix.SigT, mix.SigS[0], mix.SigP, mix.chi,
    )
    np.testing.assert_allclose(
        k_ref, 1.875, rtol=0, atol=1e-12,
        err_msg=f"transfer-matrix k_inf {k_ref} != 1.875 (mixture A 2g)",
    )

    # Reflective slab → infinite medium. Thick enough that the
    # eigenvalue is the bulk k_inf (BC leakage negligible).
    nx = 12
    mesh = Mesh1D(
        edges=np.linspace(0.0, 50.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("reflective"),
        bc_right=BC("reflective"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    sol = solve_sn(
        {0: mix}, mesh, quad, keff_tol=1e-9, flux_tol=1e-8,
    )
    np.testing.assert_allclose(
        sol.keff, k_ref, rtol=1e-6, atol=0,
        err_msg=(
            f"cumprod-path reflective slab k_eff {sol.keff} drifted from "
            f"analytical k_inf {k_ref} beyond rtol=1e-6"
        ),
    )


# ─── NEW-2: perf measurement — cumprod speedup over wavefront (d=1) ──


@pytest.mark.slow
@pytest.mark.xfail(
    not _spine_landed(),
    reason=(
        "C3 d-generic wavefront spine not yet landed — the perf "
        "measurement needs the d=1 wavefront arm to time against. "
        "strict=False so it flips to xpass on landing."
    ),
    strict=False,
)
def test_cumprod_faster_than_wavefront_d1_loose_tripwire():
    """NEW-2: MEASURE (record + LOOSE tripwire) the cumprod optimization's
    speedup over the naive d-generic wavefront walk at d = 1.

    This is the justification for keeping cumprod as the d=1 default. It
    is NOT a hard correctness gate — wall-clock is noisy (cf. the
    Phase-5b "Leg-2 loose discipline": a peak/time benchmark must stay
    loose, marked slow, final-gate only). The cumprod scan is
    O(N/log N)-work parallel-prefix (three numpy ops, no per-cell Python
    walk); the wavefront does nx sequential level-by-level forward
    substitutions. On a long-enough chain the ratio is large.

    The tripwire fires ONLY if cumprod is not even modestly faster
    (catching "the wavefront accidentally became the d=1 default and 1-D
    got slow"). The bound is deliberately VERY loose: cumprod must beat
    wavefront by at least 2× on a long chain. The printed ratio is the
    real deliverable for the perf record.
    """
    nx = 4096                                   # long chain — cumprod shines
    sn_mesh = _slab_sn_mesh(nx, bc="vacuum")
    ng = sn_mesh.ng
    rng = np.random.default_rng(99)
    sig_t = rng.uniform(0.3, 3.0, size=(ng, nx))
    iso = rng.uniform(0.2, 1.5, size=(ng, nx))
    Q = AngularSourceSink.from_isotropic(iso, sn_mesh)
    Q_arr = Q.values

    def _time(fn, *args, repeats=5):
        # Warm up (cache build on first call), then take the min over
        # repeats (min is the least noisy wall-clock estimator).
        fn(*args)
        best = float("inf")
        for _ in range(repeats):
            bf = BoundaryFlux.zeros_on(sn_mesh)
            t0 = time.perf_counter()
            fn(*((args[0], args[1], args[2], bf) if fn is _cumprod_1d_sweep
                 else (args[0], args[1], args[2], bf)))
            best = min(best, time.perf_counter() - t0)
        return best

    t_cumprod = _time(_cumprod_1d_sweep, Q, sig_t, sn_mesh)
    t_wavefront = _time(_wavefront_1d_sweep, Q_arr, sig_t, sn_mesh)

    ratio = t_wavefront / t_cumprod
    # Perf RECORD — surfaced via -s / captured output. Not an assertion.
    print(
        f"\n[NEW-2 perf] d=1 nx={nx}: cumprod={t_cumprod*1e3:.3f} ms, "
        f"wavefront={t_wavefront*1e3:.3f} ms, speedup={ratio:.1f}×"
    )

    # LOOSE tripwire only (final-gate). 2× is far below the expected
    # ratio; it fires only on a gross regression (wavefront ≈ cumprod ⟹
    # the optimization was bypassed).
    if ratio < 2.0:
        pytest.fail(
            f"cumprod is not meaningfully faster than the d=1 wavefront "
            f"(speedup {ratio:.2f}× < 2×) — the cumprod optimization may "
            f"have been bypassed for d=1 production."
        )


# Keep DiamondDifference import load-bearing (the spine's cell kernel).
def test_diamond_difference_importable():
    """Trivial guard: the spine's d-generic cell kernel is importable."""
    dd = DiamondDifference()
    if dd is None:
        pytest.fail("DiamondDifference() returned None")
