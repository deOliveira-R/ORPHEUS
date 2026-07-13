"""Diagnostic: pin the root cause of the SI-CYL-20cell NaN to
``ordinate_scan`` Blelloch-form numerical breakdown.

Created by numerics-investigator on 2026-05-28.

Three independent assertions prove the chain:

1. ``a_attenuation[ord=72, g=1, chain_pos=19]`` is EXACTLY 0 in the
   collision cache for the failing configuration.  Verified
   algebraically: at the pole cell ``A_down = 0`` so the formula
   ``a = 2|μ|·A_total / (2|μ|·A_down + dA_w·c_out + Σ_t·V) − 1``
   reduces to ``2|μ|·A_total / (dA_w·c_out + Σ_t·V) − 1`` and lands
   on the algebraic identity ``2|μ|·A_total = dA_w·c_out + Σ_t·V``
   at this specific (μ=1/√20, dr=0.1, Σ_t=1.0) point.

2. ``np.cumprod(a, axis=0)`` collapses to 0 at chain_pos 19, and
   the subsequent ``b / cumprod_a`` produces NaN; the final
   ``cumprod_a * (psi_0 + cumsum(...))`` is ``0 · NaN = NaN``.

3. The mathematically-equivalent explicit ``psi[i+1] = a·psi[i] + b[i]``
   loop returns a FINITE value (``b[-1]``, since the last ``a`` is 0
   and the in-cell update reduces to ``psi = b``).

The closed-form algorithm is *correct in real arithmetic* but
numerically catastrophic in IEEE-754 when any chain entry equals 0.

If this test fails (i.e. ``a_attenuation`` is no longer exactly zero
at the resonance OR ``ordinate_scan`` no longer NaNs when a chain
entry is 0), the underlying behaviour has shifted — revisit the
diagnosis.  When the fix lands (a Blelloch implementation that does
not divide by ``cumprod_a``), update this file to pin the fix and
remove the NaN assertion.
"""
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pytest

from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.numerics.quadrature import Quadrature
from orpheus.derivations.common.xs_library import get_mixture
from orpheus.sn.geometry import SNMesh
from orpheus.sn.sweep.cache import GeometryCoefficients, CollisionCache
from orpheus.sn.sweep.scan import ordinate_scan


def _build_caches(thickness: float = 2.0, n_cells: int = 20, eg: str = '2g'):
    fuel = get_mixture('A', eg)
    geom_obj = StructuredGeometry(
        geometry='CYL',
        regions=(Region(mat_id=0, outer_thickness_cm=thickness),),
        bcs=(BC.reflective,),
    )
    mesh = Mesh1D.from_geometry(geom_obj, region_meshes=(RegionMesh(n_cells=n_cells),))
    quad = Quadrature.level_symmetric(sn_order=8)
    sn_mesh = SNMesh(mesh, quad, materials={0: fuel})
    geom = GeometryCoefficients.from_mesh_and_quad(sn_mesh)
    ng = fuel.SigT.shape[0]
    sig_t_1d = np.broadcast_to(fuel.SigT[:, None], (ng, n_cells)).copy()
    coll = CollisionCache.from_geometry(geom, sig_t_1d)
    return sn_mesh, geom, coll, fuel


def test_attenuation_exactly_zero_at_resonance():
    """``a_attenuation[72, 1, 19] == 0`` exactly at the resonance.

    The pole cell has ``A_down = 0`` (the inner radial face at r=0 has
    zero area).  The pole-cell ``a`` formula reduces to
    ``2|μ|·A_total / (dA_w·c_out + Σ_t·V) − 1``.  At ordinate 72
    (μ_x = -1/√20), in mixture A group 1 (Σ_t = 1), with dr = 0.1
    (thick=2 / n=20), the numerator and denominator are algebraically
    identical → a = 0.
    """
    _sn_mesh, _geom, coll, _ = _build_caches()
    a = coll.a_attenuation
    assert a.shape == (80, 2, 20), f"shape: {a.shape}"
    a_resonance = a[72, 1, 19]
    assert a_resonance == 0.0, (
        f"a_attenuation[72,1,19] = {a_resonance!r}, expected exactly 0.0 "
        f"(pole-cell algebraic resonance)"
    )


def test_ordinate_scan_returns_nan_when_a_contains_zero():
    """``ordinate_scan`` produces NaN whenever any chain entry of ``a`` is
    exactly 0 — because ``b/cumprod_a`` divides by 0.

    Reproduces the exact failing chain that triggers the user bug.
    """
    _sn_mesh, _geom, coll, _ = _build_caches()
    n_ord, n_grp = 72, 1
    a_chain = coll.a_attenuation[n_ord, n_grp]  # (20,)
    assert a_chain[-1] == 0.0, "Pre-condition: last entry must be 0"

    # Build a representative b and psi_0; the result will NaN regardless
    # of their values once cumprod_a hits 0.
    b_chain = np.full_like(a_chain, 0.5)
    psi_0 = np.array(0.3)

    out = ordinate_scan(a_chain[:, None], b_chain[:, None], np.atleast_1d(psi_0))
    out_flat = out.ravel()
    assert np.isnan(out_flat[-1]), (
        f"ordinate_scan last cell = {out_flat[-1]}, expected NaN. "
        f"If FINITE, the Blelloch closed-form has been fixed to avoid "
        f"the divide-by-cumprod pathology — invert this assertion."
    )


def test_explicit_loop_is_finite_when_a_contains_zero():
    """The mathematically-equivalent explicit recurrence is FINITE even
    when a chain entry is exactly 0.

    Closed loop ``psi[i+1] = a·psi[i] + b[i]`` with a[-1]=0 reduces to
    ``psi[-1] = b[-1]``, which is well-defined.
    """
    _sn_mesh, _geom, coll, _ = _build_caches()
    a_chain = coll.a_attenuation[72, 1]
    assert a_chain[-1] == 0.0
    b_chain = np.full_like(a_chain, 0.5)
    psi = 0.3
    for i in range(a_chain.size):
        psi = a_chain[i] * psi + b_chain[i]
    assert np.isfinite(psi)
    assert abs(psi - b_chain[-1]) < 1e-15, (
        f"psi at exit = {psi}, expected b[-1] = {b_chain[-1]}"
    )


def test_krylov_avoids_ordinate_scan_path():
    """Krylov ``apply``-path does NOT consume ``ordinate_scan``.

    The bug is invisible to the matvec path; only the SI ``solve`` path
    triggers it.  This test passes today only if the architectural
    invariant holds: ``ordinate_scan`` is exclusive to the SI sweep
    (``_sweep_1d_unified``).

    Verified by grep + manual review of ``orpheus/sn/operator.py``
    (apply paths) versus ``orpheus/sn/sweep.py`` (sweep paths).  This
    assertion is a structural invariant.  See
    ``orpheus/sn/sweep/scan.py`` and ``orpheus/sn/sweep.py``.
    """
    import orpheus.sn.sweep.scan as scan_mod
    import orpheus.sn.sweep as sweep_mod
    import orpheus.sn.operator as op_mod
    # ordinate_scan is exported from spatial/scan.py
    assert hasattr(scan_mod, 'ordinate_scan')
    # sweep.py imports ordinate_scan from spatial/scan.py
    assert 'ordinate_scan' in sweep_mod.__dict__, (
        "ordinate_scan must be imported in sweep.py (the failing path)"
    )
    # operator.py must NOT import ordinate_scan
    assert 'ordinate_scan' not in op_mod.__dict__, (
        "ordinate_scan must NOT be imported in operator.py "
        "(the Krylov-apply path).  If this assertion fires, the "
        "matvec path now also uses ordinate_scan — that path will "
        "ALSO hit the NaN bug."
    )
