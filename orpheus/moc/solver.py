"""2D Method of Characteristics (MoC) neutron transport solver for a pin cell.

Solves the multi-group transport equation on a 2-D pin-cell geometry
(concentric annuli inside a square lattice cell) with reflective boundary
conditions.  Uses the flat-source MOC with a product angular quadrature
(azimuthal x Tabuchi-Yamamoto polar) and exact ray tracing through
circular boundaries.

Reference: Boyd et al. (2014) "The OpenMOC Method of Characteristics",
Ann. Nucl. Energy 68, 43-52.

.. seealso:: :ref:`theory-method-of-characteristics` — Key Facts, equations, gotchas.
"""

from __future__ import annotations

import time
from dataclasses import dataclass

import numpy as np

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.geometry import CoordSystem, Mesh1D
from orpheus.numerics.convergence import IterationRecord, warn_if_unconverged
from orpheus.numerics.eigenvalue import power_iteration

from .geometry import MOCMesh
from .quadrature import MOCQuadrature
from .core import MOCSolver


# ---------------------------------------------------------------------------
# Result dataclass
# ---------------------------------------------------------------------------

@dataclass
class MoCResult:
    """Results of a Method of Characteristics calculation."""

    keff: float
    keff_history: list[float]
    flux_per_material: dict[int, np.ndarray]  # mat_id -> (ng,) volume-averaged flux
    scalar_flux: np.ndarray   # (n_regions, ng) scalar flux per FSR
    moc_mesh: MOCMesh
    eg: np.ndarray | None     # (ng+1,) energy boundaries, or None for synthetic XS
    elapsed_seconds: float

    record: IterationRecord
    """The iteration tree this solve actually ran (#340 N4).

    Ask :attr:`~orpheus.numerics.convergence.IterationRecord.fully_converged`
    to learn whether the returned eigenvalue can be trusted — it answers for
    the outer AND every inner beneath it, which a single ``converged`` flag
    structurally cannot (#340 F1).  The children are the per-outer sweep
    loops, each declaring BOTH readings it is driven by: the scalar-flux
    increment and the boundary angular flux increment (see
    :meth:`MOCSolver.solve_fixed_source`).
    """

    @property
    def ng(self) -> int:
        """Number of energy groups (derived from the flux array)."""
        return self.scalar_flux.shape[1]

    @property
    def flux_fuel(self) -> np.ndarray:
        return self.flux_per_material.get(2, np.zeros(self.ng))

    @property
    def flux_clad(self) -> np.ndarray:
        return self.flux_per_material.get(1, np.zeros(self.ng))

    @property
    def flux_cool(self) -> np.ndarray:
        return self.flux_per_material.get(0, np.zeros(self.ng))


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def solve_moc(
    materials: dict[int, Mixture],
    mesh: Mesh1D | None = None,
    n_azi: int = 16,
    n_polar: int = 3,
    ray_spacing: float = 0.05,
    max_outer: int = 500,
    keff_tol: float = 1e-6,
    flux_tol: float = 1e-5,
    inner_tol: float = 1e-8,
    max_inner_sweeps: int = 200,
) -> MoCResult:
    """Run the 2D Method of Characteristics transport calculation.

    Parameters
    ----------
    materials : dict[int, Mixture]
        Macroscopic cross sections keyed by material ID.
    mesh : Mesh1D, optional
        Cylindrical 1-D Wigner-Seitz mesh. Defaults to the standard PWR
        pin cell built via
        :meth:`StructuredGeometry.wigner_seitz_pin_cell <orpheus.geometry.structured_geometry.StructuredGeometry.wigner_seitz_pin_cell>`
        →
        :meth:`Mesh1D.from_geometry <orpheus.geometry.mesh.Mesh1D.from_geometry>`.
    n_azi : int
        Number of azimuthal angles in [0, pi).
    n_polar : int
        Number of TY polar angles per half-space (1, 2, or 3).
    ray_spacing : float
        Perpendicular spacing between parallel rays (cm).
    max_outer : int
        Maximum number of outer (power) iterations.
    keff_tol, flux_tol : float
        Outer convergence tolerances (``|Δk|`` and relative ``‖Δφ‖``).
    inner_tol : float
        Inner convergence tolerance, applied to BOTH the scalar-flux and
        the boundary-angular-flux increments (#340 N4).  Mirrors CP's
        ``inner_tol``: three decades tighter than the outer, so the
        outer's stop test is not measuring inner noise.
    max_inner_sweeps : int
        Budget on transport sweeps per outer iteration — a CAP, not a
        schedule (⛔ was a fixed ``n_inner_sweeps=15`` count with no
        stopping test until 2026-08-11).  `[M]` a converged outer needs 1
        sweep; a cold boundary flux needs 80–110, so the default carries
        ~2x headroom and costs nothing once the loop can stop early.
    """
    t_start = time.perf_counter()

    if mesh is None:
        from orpheus.geometry import (
            Mesh1D as _M,
            RegionMesh as _RM,
            StructuredGeometry as _SG,
        )
        _geom = _SG.wigner_seitz_pin_cell()
        mesh = _M.from_geometry(_geom, region_meshes=(
            _RM(n_cells=10),  # fuel
            _RM(n_cells=3),   # clad
            _RM(n_cells=7),   # cool
        ))

    _any_mat = next(iter(materials.values()))
    eg = _any_mat.eg
    ng = _any_mat.ng

    # Build augmented geometry
    quad = MOCQuadrature.create(n_azi=n_azi, n_polar=n_polar)
    moc_mesh = MOCMesh(mesh, quad, ray_spacing=ray_spacing)

    print(f"  MOC: {moc_mesh.n_regions} regions, {len(moc_mesh.tracks)} tracks, "
          f"{n_azi} azi x {n_polar} polar, spacing={ray_spacing}")

    # Eigenvalue solve
    solver = MOCSolver(
        moc_mesh, materials,
        keff_tol=keff_tol, flux_tol=flux_tol,
        inner_tol=inner_tol, max_inner_sweeps=max_inner_sweeps,
    )
    outcome = power_iteration(solver, max_iter=max_outer, budget_name="max_outer")
    keff, keff_history, phi = (
        outcome.keff, outcome.keff_history, outcome.flux_distribution,
    )

    # Post-processing: volume-averaged spectra per material
    nr = moc_mesh.n_regions
    unique_mats = set(int(m) for m in moc_mesh.region_mat_ids)
    vol_per_mat = {m: 0.0 for m in unique_mats}
    flux_per_mat = {m: np.zeros(ng) for m in unique_mats}

    for i in range(nr):
        mat_id = int(moc_mesh.region_mat_ids[i])
        A_i = moc_mesh.region_areas[i]
        vol_per_mat[mat_id] += A_i
        flux_per_mat[mat_id] += phi[i, :] * A_i

    for m in unique_mats:
        if vol_per_mat[m] > 0:
            flux_per_mat[m] /= vol_per_mat[m]

    elapsed = time.perf_counter() - t_start
    print(f"  Elapsed: {elapsed:.1f}s")

    # ⭐ #340 N4.7 — see the note at the identical call in
    # :func:`~orpheus.cp.solver.solve_cp` for why this sits directly in the
    # entry and passes no ``balance_defect``.
    #
    # ⚠ MoC is the family where this warning has the most to say: its inner
    # is a *within-outer sweep* loop whose binding mode is the CYCLIC TRACK
    # boundary flux, not the volume moment (#340 N4 — the ‖Δφ‖-only criterion
    # was Mode-12 BLIND to it, `[M]` φ machine-zero at sweep 4 while ψ_b was
    # still at 1e-3).  A starved MoC inner therefore looks fine in every
    # volume diagnostic, which is exactly the class this warning exists for.
    # ``balance_defect=None`` explicitly — see the note at ``solve_cp``.
    warn_if_unconverged(outcome.record, where="solve_moc", balance_defect=None)

    return MoCResult(
        keff=keff_history[-1],
        keff_history=keff_history,
        flux_per_material=flux_per_mat,
        scalar_flux=phi,
        moc_mesh=moc_mesh,
        eg=eg,
        elapsed_seconds=elapsed,
        record=outcome.record,
    )
