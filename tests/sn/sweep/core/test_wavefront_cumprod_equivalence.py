r"""CumprodScan(d=1) ≡ FullFieldWavefront(d=1) — the optimization-vs-spine oracle.

The LossRepresentation carve (S3) makes the d-generic wavefront DAG walk
(:class:`~orpheus.sn.loss_representation.FullFieldWavefront`) a *selectable*
strategy — the verification SPINE — and the 1-D Blelloch parallel-prefix scan
(:class:`~orpheus.sn.loss_representation.CumprodScan`) the d=1 production
OPTIMIZATION. This file proves the optimization computes the SAME 1-D transport
sweep as the spine, through the real strategy API (no hand-built adapters).

Bit-identical or principled-equivalence?  → PRINCIPLED-EQUIVALENCE.
==================================================================
The cumprod scan and the wavefront forward-substitution walk compute the SAME
affine recurrence per ordinate

    psi_out[i] = 2·psi_bar − psi_in[i],
    psi_bar    = (Q + s·psi_in[i]) / (σ_t + s)

but in DIFFERENT floating-point ASSOCIATION ORDERS:

* cumprod path: ``cumprod_a · (psi_0 + cumsum(b / cumprod_a))`` — a
  parallel-prefix reduction TREE over all nx cells at once.
* wavefront spine: cell-by-cell forward substitution along the DAG levels —
  each cell's ``psi_in`` is the previous cell's ``psi_out``, evaluated
  sequentially.

Because IEEE-754 addition/multiplication is NOT associative, these two
real-arithmetic-identical computations disagree at ULP level. Per
``vv-principles`` §"Bit-identity vs principled-equivalence", this is a
principled-equivalence-at-nulp case — NOT a bug — and ALL THREE criteria hold:

  (1) Principled at every step. Both paths' intermediates are named
      reactor-physics quantities (the per-cell downstream face flux and the
      cell average); neither is an anonymous reduction.
  (2) Verified against a structurally-independent reference. The anchor is the
      ANALYTICAL closed form ``k_inf = λ_max(A⁻¹F) = 1.875`` for mixture A 2g
      (reflective infinite-medium eigenvalue), computed independently by
      ``kinf_and_spectrum_homogeneous`` (transfer-matrix), see
      ``test_cumprod_path_hits_analytical_kinf`` below. The cumprod path — the
      path the spine is validated against — is ITSELF anchored to closed-form
      ground truth, so "spine ≡ cumprod" transitively anchors the spine.
  (3) Drift is FP-non-associativity, dimensionally explainable: a single-sweep
      computation is bounded by ``(reduction depth) × ULP`` ≈ ``nx × ε``.

Layout
======
Both strategies consume and return the genuine rank-1 ``(N, ng, nx)`` /
``(ng, nx)`` layout (spatial rank == ndim — the dimension-agnostic target the
C-phases landed; the legacy reduced-mesh phantom-y axis is gone). The source
``AngularSourceSink.from_isotropic`` projects to the SAME ``(N, ng, nx)`` both
strategies read, so they are directly comparable (no layout bridge).

Cardinal Rule 6 + vv anti-patterns #3/#4 / Mode 9 — what this STRESSES
=====================================================================
A 1-group flat-flux 1-D test proves NOTHING (H1 1-group + H2 flat
degeneracies). The equivalence test therefore uses ≥2 GROUPS (mixture A 2g,
asymmetric downscatter), a heterogeneous per-cell Σ_t, a non-uniform source,
AND a non-zero boundary inflow — so ``psi_in ≠ psi_out`` at every cell (the
redistribution term is driven OUT of cancellation), an anisotropic /
heterogeneous config (vv Mode 9: a splitting/optimization FP-invariance must be
verified off the degenerate isotropic-reflective box).

``-O`` discipline (vv Mode 8): every assertion is a ``np.testing.*`` function
call — fires under the canonical ``python -O``.

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
from orpheus.sn.loss_representation import CumprodScan, FullFieldWavefront
from orpheus.transport.fields.boundary_flux import BoundaryFlux
from orpheus.transport.source_sinks import AngularSourceSink

pytestmark = pytest.mark.foundation


# ─── nulp bound: reduction-depth × ULP × safety ─────────────────────
#
# Per vv-principles §"Bit-identity vs principled-equivalence" criterion (3): a
# single-sweep FP drift is bounded by (reduction depth) × ULP. The cumprod
# chain's reduction depth is the cell count nx; the wavefront does nx sequential
# forward substitutions. With nx ≤ 16, depth ≈ 16 plus a ×8 safety factor
# (covering the 2-group coupling + the source/BC affine terms) gives 128.
_NULP_BOUND = 128

_NX_EQUIV = 12           # chain length for the per-sweep equivalence
_LENGTH = 6.0            # slab thickness (cm)


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


def _seeded_inflow(sn_mesh: SNMesh, rng: np.random.Generator) -> BoundaryFlux:
    """A boundary flux with a random non-zero inflow trace on every face."""
    bf = BoundaryFlux.zeros_on(sn_mesh)
    for face in bf.layout.faces:
        fv = bf.face_view(face)
        fv[...] = rng.uniform(0.0, 1.0, size=fv.shape)
    return bf


# ─── NEW-1: CumprodScan(d=1) ≡ FullFieldWavefront(d=1), principled @ nulp ───


@pytest.mark.parametrize("bc", ["vacuum", "reflective"])
def test_cumprod_1d_equals_full_field_spine(bc):
    """The d=1 cumprod OPTIMIZATION computes the SAME single sweep as the
    d-generic full-field SPINE — principled-equivalence at nulp.

    Strategy-vs-strategy through the real API (no hand adapters): the production
    :class:`CumprodScan` vs the verification :class:`FullFieldWavefront`, both
    selected for the SAME 1-D mesh, fed the SAME per-ordinate source / Σ_t /
    boundary inflow. What this catches: a convention drift between the spine's
    d-generic single-axis closure (the DAG walk + whole-trace BC seed/absorb)
    and the cumprod recurrence (the chain seed) — sign of the streaming coeff,
    the BC-inflow seed face, the group axis.
    """
    nx = _NX_EQUIV
    sn_mesh = _slab_sn_mesh(nx, bc=bc)
    ng = sn_mesh.ng
    rng = np.random.default_rng(20260610)

    # Heterogeneous Σ_t (ng, nx) + non-uniform source — the non-flat drivers.
    sig_t = rng.uniform(0.3, 3.0, size=(ng, nx))
    iso = rng.uniform(0.2, 1.5, size=(ng, nx))
    Q = AngularSourceSink.from_isotropic(iso, sn_mesh)
    Q_arr = Q.values                     # (N, ng, nx) — both strategies, rank-1

    # Separate seeded inflow per strategy (each sweep mutates its own trace).
    bf_cumprod = _seeded_inflow(sn_mesh, rng)
    bf_spine = BoundaryFlux.from_mesh(bf_cumprod.values.copy(), sn_mesh)

    ang_c, scal_c = CumprodScan(sn_mesh).sweep(Q_arr, sig_t, bf_cumprod)
    ang_s, scal_s = FullFieldWavefront(sn_mesh).sweep(Q_arr, sig_t, bf_spine)

    np.testing.assert_array_almost_equal_nulp(ang_s, ang_c, nulp=_NULP_BOUND)
    np.testing.assert_array_almost_equal_nulp(scal_s, scal_c, nulp=_NULP_BOUND)


# ─── NEW-1 anchor: the cumprod path hits the analytical k_inf ────────


def test_cumprod_path_hits_analytical_kinf():
    """NEW-1 anchor: the cumprod 1-D path converges to the closed-form
    infinite-medium eigenvalue ``k_inf = λ_max(A⁻¹F) = 1.875`` for mixture A 2g
    (reflective slab → infinite medium).

    This is the STRUCTURALLY-INDEPENDENT ground (vv criterion 2) that anchors
    the whole oracle: old-vs-new ULP distance is necessary but never sufficient,
    so the path the spine is validated against must itself terminate in a
    closed-form reference. ``solve_sn`` on a reflective slab drives the cumprod
    sweep; the transfer-matrix ``kinf_and_spectrum_homogeneous`` is the
    independent analytical value (different structural angle: eigendecomposition
    of A⁻¹F, no sweep at all).

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

    # Reflective slab → infinite medium. Thick enough that the eigenvalue is the
    # bulk k_inf (BC leakage negligible).
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


# ─── NEW-2: perf measurement — cumprod speedup over the spine (d=1) ──


@pytest.mark.slow
def test_cumprod_faster_than_full_field_spine_d1():
    """NEW-2: MEASURE (record + LOOSE tripwire) the cumprod optimization's
    speedup over the d-generic full-field spine at d = 1 — the justification for
    keeping cumprod as the d=1 default.

    NOT a hard correctness gate — wall-clock is noisy (slow, final-gate only).
    The cumprod scan is O(N/log N)-work parallel-prefix (three numpy ops, no
    per-cell Python walk); the spine does nx sequential level-by-level forward
    substitutions. On a long chain the ratio is large. The tripwire fires ONLY
    if cumprod is not even modestly faster (catching "the spine accidentally
    became the d=1 default and 1-D got slow").
    """
    nx = 4096                                   # long chain — cumprod shines
    sn_mesh = _slab_sn_mesh(nx, bc="vacuum")
    ng = sn_mesh.ng
    rng = np.random.default_rng(99)
    sig_t = rng.uniform(0.3, 3.0, size=(ng, nx))
    iso = rng.uniform(0.2, 1.5, size=(ng, nx))
    Q = AngularSourceSink.from_isotropic(iso, sn_mesh)
    Q_arr = Q.values

    cumprod = CumprodScan(sn_mesh)
    spine = FullFieldWavefront(sn_mesh)

    def _time(strategy, repeats=5):
        bf = BoundaryFlux.zeros_on(sn_mesh)
        strategy.sweep(Q_arr, sig_t, bf)        # warm up (cache build)
        best = float("inf")
        for _ in range(repeats):
            bf = BoundaryFlux.zeros_on(sn_mesh)
            t0 = time.perf_counter()
            strategy.sweep(Q_arr, sig_t, bf)
            best = min(best, time.perf_counter() - t0)
        return best

    t_cumprod = _time(cumprod)
    t_spine = _time(spine)
    ratio = t_spine / t_cumprod
    print(
        f"\n[NEW-2 perf] d=1 nx={nx}: cumprod={t_cumprod*1e3:.3f} ms, "
        f"spine={t_spine*1e3:.3f} ms, speedup={ratio:.1f}×"
    )
    if ratio < 2.0:
        pytest.fail(
            f"cumprod is not meaningfully faster than the d=1 spine "
            f"(speedup {ratio:.2f}× < 2×) — the cumprod optimization may have "
            f"been bypassed for d=1 production."
        )
