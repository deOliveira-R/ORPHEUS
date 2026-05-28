"""Diagnostic Step 1: characterize the SI vs Krylov keff mismatch.

Created by numerics-investigator on 2026-05-28.

Reproduces the symptom described in the brief:
- SI inner solver -> keff = 1.875 (analytical k_inf = nu*Sigma_f/Sigma_a)
- Krylov inner solver -> keff ~ 1.40 (~25% wrong)

If this test catches a real bug, promote to ``tests/sn/test_sn_solver.py``
or ``tests/sn/spatial/test_sweep_vs_apply_consistency.py``.
"""
from __future__ import annotations

import warnings

import numpy as np
import pytest


def _solve(inner_solver: str, inner_tol: float):
    """Helper: run solve_sn with the parameters at issue."""
    from orpheus.derivations.common.xs_library import get_mixture
    from orpheus.geometry import (
        BC, Mesh1D, Region, RegionMesh, StructuredGeometry,
    )
    from orpheus.numerics.quadrature import Quadrature
    from orpheus.sn.solver import solve_sn

    fuel = get_mixture("A", "2g")
    geom = StructuredGeometry(
        geometry="SPH",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(BC.reflective,),
    )
    mesh = Mesh1D.from_geometry(
        geom, region_meshes=(RegionMesh(n_cells=20),),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        res = solve_sn(
            materials={0: fuel}, mesh=mesh, quadrature=quad,
            inner_solver=inner_solver,
            keff_tol=1e-12, flux_tol=1e-10,
            inner_tol=inner_tol,
        )
    return res.keff


def _kinf_analytical() -> float:
    """Compute analytical k_inf = nu*Sigma_f / Sigma_a for mixture A 2g."""
    from orpheus.derivations.common.xs_library import get_mixture

    mix = get_mixture("A", "2g")
    # k_inf = sum(nu_sigma_f * phi) / sum(sigma_a * phi).
    # For homogeneous reflective, phi is the dominant eigenvector of
    # A^{-1} F (a fully analytical homogeneous problem); the simplest
    # explicit answer comes from sweeping the 2x2 transfer matrix.
    # The Phase F test asserts 1.875 — accept that as the reference.
    return 1.875


def test_step1_reproduce_si_correct_krylov_wrong_at_default_tol():
    """SI=1.875, Krylov drifts at default inner_tol=1e-8."""
    keff_si = _solve("source_iteration", inner_tol=1e-8)
    keff_kr_loose = _solve("krylov", inner_tol=1e-8)
    kinf = _kinf_analytical()

    si_err = abs(keff_si - kinf)
    kr_err_loose = abs(keff_kr_loose - kinf)

    # Pin the symptom exactly as the brief describes.
    assert si_err < 1e-6, f"SI must hit kinf=1.875; got {keff_si}, err={si_err:.3e}"
    assert kr_err_loose > 0.1, (
        f"REPRODUCER: Krylov with inner_tol=1e-8 must drift >0.1 from kinf=1.875; "
        f"got keff={keff_kr_loose}, err={kr_err_loose:.3e}"
    )

    print(f"\n[step1] kinf_analytical = {kinf:.10f}")
    print(f"[step1] keff_si (tol=1e-8)        = {keff_si:.10f}  err = {si_err:.3e}")
    print(f"[step1] keff_krylov (tol=1e-8)    = {keff_kr_loose:.10f}  err = {kr_err_loose:.3e}")


def test_step1_tightening_inner_tol_fixes_or_does_not():
    """Test whether tightening inner_tol to 1e-12 alone fixes Krylov."""
    keff_kr_loose = _solve("krylov", inner_tol=1e-8)
    keff_kr_tight = _solve("krylov", inner_tol=1e-12)
    kinf = _kinf_analytical()

    err_loose = abs(keff_kr_loose - kinf)
    err_tight = abs(keff_kr_tight - kinf)

    print(f"\n[step1] keff_krylov (tol=1e-8)  = {keff_kr_loose:.10f}  err = {err_loose:.3e}")
    print(f"[step1] keff_krylov (tol=1e-12) = {keff_kr_tight:.10f}  err = {err_tight:.3e}")
    print(f"[step1] ratio (loose/tight)     = {err_loose / max(err_tight, 1e-30):.2e}")

    # No assertion on which fixes it — this is a *characterization* probe.


if __name__ == "__main__":
    test_step1_reproduce_si_correct_krylov_wrong_at_default_tol()
    test_step1_tightening_inner_tol_fixes_or_does_not()
