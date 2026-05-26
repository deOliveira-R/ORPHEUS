"""Generate the frozen-reference SN regression snapshots.

This script is the **Issue 16 redesign** (per
``docs/testing/sn_verification_matrix.rst``): instead of snapshotting
every existing SN test, it captures a **minimum frozen-reference set**
proving operator-algebra equivalence under DD across all geometry,
group-count, BC, and quadrature combinations the SN solver supports.

Run::

    python -m tests.sn.regression._generate_snapshots
    python -m tests.sn.regression._generate_snapshots --case slab_2g_homogeneous_dd_n20

Each snapshot writes ``snapshots/<name>.npz`` containing ``(keff,
scalar_flux, case_name, case_description, generator_commit)`` as the
authoritative pre-reshape reference. The companion test
``test_dd_regression.py`` re-runs each case and asserts bit-for-bit
agreement.

**The snapshots are NOT a verification reference.** They are a
regression gate detecting unintended numerical drift across refactors.
End-to-end correctness verification runs through L0/L1
analytical-reference tests (``l1_analytical/test_kinf_homogeneous.py``,
the curvilinear MMS suites, F_N / Case-method cross-verifications).

When this script is run on a different solver state than the snapshots
were generated against, the reproductions will mismatch — that is the
intended behaviour. The protocol for legitimate updates: (1) audit why
the new output is correct (with V&V evidence), (2) re-run the
generator, (3) commit both the new snapshot AND the audit narrative
in the same commit.
"""
from __future__ import annotations

import argparse
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import (
    BC,
    Mesh1D,
    Mesh2D,
    Region,
    RegionMesh,
    StructuredGeometry,
)
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solver import solve_sn, solve_sn_fixed_source


SNAPSHOT_DIR = Path(__file__).parent / "snapshots"

# ─── case configuration helpers ──────────────────────────────────────


def _slab_homogeneous(ng: str, n_cells: int) -> dict:
    fuel = get_mixture("A", ng)
    geom = StructuredGeometry(
        geometry="SLB",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(BC.reflective, BC.reflective),
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=n_cells),))
    return dict(
        materials={0: fuel}, mesh=mesh,
        quadrature=Quadrature.gauss_legendre(n_ordinates=8),
        scattering_order=0,
    )


def _slab_3region(ng: str, n_cells: int) -> dict:
    """Fuel | moderator | fuel, equal thicknesses summing to 2 cm."""
    fuel = get_mixture("A", ng)
    mod = get_mixture("B", ng)
    geom = StructuredGeometry(
        geometry="SLB",
        regions=(
            Region(mat_id=0, outer_thickness_cm=0.5),
            Region(mat_id=1, outer_thickness_cm=1.0),
            Region(mat_id=0, outer_thickness_cm=0.5),
        ),
        bcs=(BC.reflective, BC.reflective),
    )
    # equal subdivision across regions: total cells split per outer thickness
    n_per_region = (n_cells // 4, n_cells // 2, n_cells // 4)
    mesh = Mesh1D.from_geometry(
        geom,
        region_meshes=tuple(RegionMesh(n_cells=n) for n in n_per_region),
    )
    return dict(
        materials={0: fuel, 1: mod}, mesh=mesh,
        quadrature=Quadrature.gauss_legendre(n_ordinates=8),
        scattering_order=0,
    )


def _sphere_homogeneous(ng: str, n_cells: int) -> dict:
    fuel = get_mixture("A", ng)
    geom = StructuredGeometry(
        geometry="SPH",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(BC.reflective,),
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=n_cells),))
    return dict(
        materials={0: fuel}, mesh=mesh,
        quadrature=Quadrature.gauss_legendre(n_ordinates=8),
        scattering_order=0,
    )


def _sphere_3region(ng: str, n_cells: int) -> dict:
    fuel = get_mixture("A", ng)
    mod = get_mixture("B", ng)
    geom = StructuredGeometry(
        geometry="SPH",
        regions=(
            Region(mat_id=0, outer_thickness_cm=0.5),
            Region(mat_id=1, outer_thickness_cm=1.0),
            Region(mat_id=0, outer_thickness_cm=0.5),
        ),
        bcs=(BC.reflective,),
    )
    n_per_region = (n_cells // 4, n_cells // 2, n_cells // 4)
    mesh = Mesh1D.from_geometry(
        geom,
        region_meshes=tuple(RegionMesh(n_cells=n) for n in n_per_region),
    )
    return dict(
        materials={0: fuel, 1: mod}, mesh=mesh,
        quadrature=Quadrature.gauss_legendre(n_ordinates=8),
        scattering_order=0,
    )


def _cylinder_homogeneous(ng: str, n_cells: int, quad_kind: str) -> dict:
    fuel = get_mixture("A", ng)
    geom = StructuredGeometry(
        geometry="CYL",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(BC.reflective,),
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=n_cells),))
    if quad_kind == "LS4":
        quadrature = Quadrature.level_symmetric(sn_order=4)
    elif quad_kind == "product_2x4":
        quadrature = Quadrature.product(n_mu=2, n_phi=4)
    else:
        raise ValueError(f"unknown cylinder quadrature kind: {quad_kind}")
    return dict(
        materials={0: fuel}, mesh=mesh, quadrature=quadrature,
        scattering_order=0,
    )


def _cylinder_3region(ng: str, n_cells: int, quad_kind: str) -> dict:
    fuel = get_mixture("A", ng)
    mod = get_mixture("B", ng)
    geom = StructuredGeometry(
        geometry="CYL",
        regions=(
            Region(mat_id=0, outer_thickness_cm=0.5),
            Region(mat_id=1, outer_thickness_cm=1.0),
            Region(mat_id=0, outer_thickness_cm=0.5),
        ),
        bcs=(BC.reflective,),
    )
    n_per_region = (n_cells // 4, n_cells // 2, n_cells // 4)
    mesh = Mesh1D.from_geometry(
        geom,
        region_meshes=tuple(RegionMesh(n_cells=n) for n in n_per_region),
    )
    if quad_kind == "LS4":
        quadrature = Quadrature.level_symmetric(sn_order=4)
    elif quad_kind == "product_2x4":
        quadrature = Quadrature.product(n_mu=2, n_phi=4)
    else:
        raise ValueError(f"unknown cylinder quadrature kind: {quad_kind}")
    return dict(
        materials={0: fuel, 1: mod}, mesh=mesh, quadrature=quadrature,
        scattering_order=0,
    )


def _slab_p1_aniso(ng: str, n_cells: int) -> dict:
    """Slab with B mixture (μ_bar=0.6, strongly anisotropic, P1 data)."""
    mod = get_mixture("B", ng)
    geom = StructuredGeometry(
        geometry="SLB",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(BC.reflective, BC.reflective),
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=n_cells),))
    return dict(
        materials={0: mod}, mesh=mesh,
        quadrature=Quadrature.gauss_legendre(n_ordinates=8),
        scattering_order=1,
    )


def _sphere_p1_aniso(ng: str, n_cells: int) -> dict:
    """Sphere with B mixture; activates Pℓ Galerkin path on curvilinear sweep."""
    mod = get_mixture("B", ng)
    geom = StructuredGeometry(
        geometry="SPH",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(BC.reflective,),
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=n_cells),))
    return dict(
        materials={0: mod}, mesh=mesh,
        quadrature=Quadrature.gauss_legendre(n_ordinates=8),
        scattering_order=1,
    )


def _cartesian_2d(ng: str, n_per_side: int) -> dict:
    """2D Cartesian homogeneous, LS_4 (12 ordinates), single material.

    Pins the 2D wavefront sweep diagonal scheduling against any
    refactor that touches ``orpheus/sn/sweep.py::_sweep_2d_wavefront``.
    """
    fuel = get_mixture("A", ng)
    L = 2.0
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, L, n_per_side + 1),
        edges_y=np.linspace(0.0, L, n_per_side + 1),
        mat_map=np.zeros((n_per_side, n_per_side), dtype=int),
        bc_xmin=BC.reflective, bc_xmax=BC.reflective,
        bc_ymin=BC.reflective, bc_ymax=BC.reflective,
    )
    return dict(
        materials={0: fuel}, mesh=mesh,
        quadrature=Quadrature.level_symmetric(sn_order=4),
        scattering_order=0,
    )


def _slab_fixed_source(ng: str, n_cells: int) -> dict:
    """Slab vacuum-BC fixed-source with uniform isotropic external Q.

    No fission. The reference output is :math:`\\phi`-only (no k_eff);
    the ``kind="fixed_source"`` flag dispatches this case to
    ``solve_sn_fixed_source``.
    """
    fuel = get_mixture("A", ng)
    L = 2.0
    n_ord = 8
    n_groups = 1 if ng == "1g" else 2 if ng == "2g" else 4
    mesh = Mesh1D(
        edges=np.linspace(0.0, L, n_cells + 1),
        mat_ids=np.zeros(n_cells, dtype=int),
        bc_left=BC.vacuum,
        bc_right=BC.vacuum,
    )
    quadrature = Quadrature.gauss_legendre(n_ordinates=n_ord)
    # R-1 Step 4 A1 — ``external_source`` is **per-ordinate density**
    # (already projected via ``/sum_w`` at the caller boundary).
    # Iso scalar source magnitude 1.0 ⇒ per-ord density ``1.0/sum_w``.
    # The snapshot remains bit-identical because the scalar flux only
    # depends on the iso magnitude: ``Σ_n w_n · ψ_n = sum_w · (1/sum_w)/Σ_t``.
    # Issue #196 PR-INDEX-5: principled shape ``(N, ng, nx, ny)``.
    sum_w = float(quadrature.weights.sum())
    external_source = np.full(
        (n_ord, n_groups, n_cells, 1), 1.0 / sum_w,
    )
    return dict(
        materials={0: fuel}, mesh=mesh, quadrature=quadrature,
        scattering_order=0, kind="fixed_source",
        external_source=external_source,
    )


# ─── snapshot case registry ──────────────────────────────────────────


@dataclass(frozen=True)
class SnapshotCase:
    """One regression snapshot definition.

    The ``builder`` callable returns a kwargs dict consumed by both the
    generator and the regression test, so the two sides cannot drift
    apart in mesh / quadrature / material configuration.
    """
    name: str
    description: str
    builder: Callable[[], dict]


CASES: tuple[SnapshotCase, ...] = (
    SnapshotCase(
        "slab_2g_homogeneous_dd_n20",
        "DD slab 2G homogeneous, GL-8, n=20",
        lambda: _slab_homogeneous("2g", 20),
    ),
    SnapshotCase(
        "slab_2g_3reg_dd_n40",
        "DD slab fuel|moderator|fuel, 2G, GL-8, n=40",
        lambda: _slab_3region("2g", 40),
    ),
    SnapshotCase(
        "sphere_2g_homogeneous_dd_n20",
        "DD sphere 2G homogeneous, GL-8, n=20",
        lambda: _sphere_homogeneous("2g", 20),
    ),
    SnapshotCase(
        "sphere_2g_3reg_dd_n40",
        "DD sphere fuel|moderator|fuel, 2G, GL-8, n=40",
        lambda: _sphere_3region("2g", 40),
    ),
    SnapshotCase(
        "cyl_1g_homogeneous_LS4_dd_n20",
        "DD cylinder 1G homogeneous, LS_4 (12 ord), n=20",
        lambda: _cylinder_homogeneous("1g", 20, "LS4"),
    ),
    SnapshotCase(
        "cyl_1g_homogeneous_product_dd_n20",
        "DD cylinder 1G homogeneous, ProductQuadrature(2x4) (8 ord), n=20",
        lambda: _cylinder_homogeneous("1g", 20, "product_2x4"),
    ),
    SnapshotCase(
        "cyl_2g_3reg_LS4_dd_n40",
        "DD cylinder fuel|moderator|fuel, 2G, LS_4, n=40",
        lambda: _cylinder_3region("2g", 40, "LS4"),
    ),
    SnapshotCase(
        "slab_2g_p1_aniso_dd_n20",
        "DD slab 2G + B mixture P1 anisotropic scattering, GL-8, n=20",
        lambda: _slab_p1_aniso("2g", 20),
    ),
    SnapshotCase(
        "sphere_2g_p1_aniso_dd_n20",
        "DD sphere 2G + B mixture P1 anisotropic scattering, GL-8, n=20",
        lambda: _sphere_p1_aniso("2g", 20),
    ),
    SnapshotCase(
        "2d_1g_LS4_dd_15x15",
        "DD 2D Cartesian 1G homogeneous, LS_4 (12 ord), 15x15",
        lambda: _cartesian_2d("1g", 15),
    ),
    SnapshotCase(
        "slab_fixed_source_dd_n20",
        "DD slab fixed-source 1G, vacuum BCs, uniform isotropic Q, GL-8, n=20",
        lambda: _slab_fixed_source("1g", 20),
    ),
)

# Cases queued for follow-up snapshot work (remain unimplemented here):
#   - cyl_white_bc_dd_n20       (post-Issue 7 — BoundaryOperator tensor decomposition)


# ─── snapshot generation / loading ───────────────────────────────────


def _git_short_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short=12", "HEAD"],
            stderr=subprocess.DEVNULL,
        ).decode().strip()
    except Exception:
        return "unknown"


def run_case(cfg: dict):
    """Run the configured case (eigen or fixed-source) and return result.

    Shared between snapshot generation and the regression test so the
    two paths cannot drift in solver-tolerance configuration.
    """
    kind = cfg.get("kind", "eigen")
    if kind == "fixed_source":
        return solve_sn_fixed_source(
            cfg["materials"], cfg["mesh"], cfg["quadrature"],
            cfg["external_source"],
            scattering_order=cfg.get("scattering_order", 0),
            max_inner=500, inner_tol=1e-12,
        )
    if kind == "eigen":
        return solve_sn(
            cfg["materials"], cfg["mesh"], cfg["quadrature"],
            scattering_order=cfg.get("scattering_order", 0),
            max_outer=500, keff_tol=1e-12, flux_tol=1e-10,
            max_inner=300, inner_tol=1e-10,
        )
    raise ValueError(f"unknown case kind: {kind!r}")


def generate_one(case: SnapshotCase, *, sha: str | None = None) -> Path:
    """Run the case and write the snapshot .npz file."""
    cfg = case.builder()
    result = run_case(cfg)
    SNAPSHOT_DIR.mkdir(parents=True, exist_ok=True)
    out = SNAPSHOT_DIR / f"{case.name}.npz"

    payload: dict = dict(
        scalar_flux=np.asarray(result.scalar_flux.values, dtype=np.float64),
        case_name=np.array(case.name),
        case_description=np.array(case.description),
        case_kind=np.array(cfg.get("kind", "eigen")),
        generator_commit=np.array(sha or _git_short_sha()),
    )
    if cfg.get("kind", "eigen") == "eigen":
        payload["keff"] = np.float64(result.keff)

    np.savez_compressed(out, **payload)
    return out


def generate_all(case_names: list[str] | None = None) -> list[Path]:
    sha = _git_short_sha()
    written = []
    for case in CASES:
        if case_names and case.name not in case_names:
            continue
        path = generate_one(case, sha=sha)
        written.append(path)
        print(f"wrote  {path.relative_to(Path.cwd())}")
    return written


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate frozen-reference SN regression snapshots.",
    )
    parser.add_argument(
        "--case", action="append", default=None,
        help="Restrict to a specific case name (repeatable).",
    )
    parser.add_argument(
        "--list", action="store_true",
        help="List available cases and exit.",
    )
    args = parser.parse_args()
    if args.list:
        for case in CASES:
            print(f"  {case.name:50s}  {case.description}")
        return
    written = generate_all(case_names=args.case)
    print(f"generated {len(written)} snapshot(s) in {SNAPSHOT_DIR}")


if __name__ == "__main__":
    main()
