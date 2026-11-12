r"""The #707 Step A acceptance problems, as DATA: what must not move when the data side V* is stored integrated.

Not a test module (no ``test_`` prefix). Two readers share this file and nothing else:

* ``_capture_true_ends_acceptance.py`` solves every case on the pre-carve tree and freezes the answers into
  ``true_ends_acceptance.npz`` (with the commit, the numpy version and the platform);
* ``test_true_ends_acceptance.py`` solves the SAME cases on the current tree and requires each frozen key bit-equal
  (``np.array_equal``), unless a declared event in ``true_ends_acceptance_rebaselines.json`` overrides it (read through
  :func:`expected`, the one reading of the fixture).

What is frozen is the STATE side only (R6: k and fluxes do not move): ``k`` (eigen and adjoint eigen), the scalar
flux, the bulk angular flux, the adjoint importance. Source VALUES are deliberately absent: R6 moves every one of them
by the spatial measure, and their migration is the chart-law gate (``tests/gates/sn/test_source_chart_law.py``).

Population (X2). 19 cases. SN 15: geometry slab 6, sphere 3, cylinder 1, 2-D Cartesian 5; scheme DD 12, LD 3 (the
registry holds two schemes, ``diamond_difference`` and ``linear_discontinuous``; there is no Step scheme in the tree,
`[M]` ``DiscretizationSchemeBase.registry`` 2026-11-12); direction forward 10, adjoint 5; question eigen 6, fixed
source 9; scattering order 0 on 11, P1 on 4. Diffusion 3 (slab and sphere eigen, slab fixed source). Homogeneous 1.
Every multi-group case is 2-group and heterogeneous except the homogeneous baseline and the 1-group vacuum slab
(leakage-dominated). LD has no eigen entry (``solve_sn`` takes no ``scheme``), so its rows are fixed-source only.

The cases already frozen elsewhere are not re-captured here (reuse > new): ``carry_head.npz`` holds the 1-D SN
eigen/fixed k and fluxes on slab, cylinder, sphere and hollow cylinder under every wall law, and the diffusion, CP,
MoC and MC k; ``tests/gates/sn/regression/snapshots`` holds the 14 DD regression cases. This module adds the gaps:
LD, the adjoint (eigen and fixed source), 2-D P1, diffusion fixed source, homogeneous. Several builders are the
regression module's own (``_generate_snapshots``), reused at smaller sizes so the capture stays under a minute.

Every setting an answer depends on is spelled here (tolerances, budgets, the inner solver), so a bit-equal value
is meaningful only because the solver is the same solver.
"""

from __future__ import annotations

import dataclasses
import json
from functools import cache
from pathlib import Path
from typing import Any, Callable

import numpy as np

FIXTURE = Path(__file__).with_name("true_ends_acceptance.npz")
REBASELINES = Path(__file__).with_name("true_ends_acceptance_rebaselines.json")

EIGEN: dict[str, Any] = dict(max_outer=500, keff_tol=1e-10, flux_tol=1e-9, max_inner=2000, inner_tol=1e-11)
FIXED: dict[str, Any] = dict(max_inner=3000, inner_tol=1e-12)


@dataclasses.dataclass(frozen=True)
class Case:
    """One frozen problem: ``solve()`` returns ``{field: ndarray}``, the keys stored as ``<name>/<field>``."""

    name: str
    solve: Callable[[], dict[str, np.ndarray]]
    #: The family a re-baseline law is stated for: ``state`` (bit-identity required), ``adjoint_eigen`` (the gauge
    #: restatement of census C8 rescales the importance; a declared ``gauge`` event is the only admitted move).
    family: str = "state"


# ── shared fixtures ──────────────────────────────────────────────────────────────────────────────────────────────


def _regression():
    from tests.gates.sn.regression import _generate_snapshots as g

    return g


def _fuel_mod():
    from orpheus.data.macro_xs.synthetic import get_mixture

    return get_mixture("A", "2g"), get_mixture("B", "2g")


def _slab_het_mesh(n: int, faces: tuple[Any, Any]):
    """Fuel | moderator slab, 2 regions of ``n`` cells each (``uniform_width`` counts per region), the two face roles given."""
    from orpheus.geometry import CoordSystem, StructuredGeometry
    from orpheus.mesh import CellsByCount, Mesher

    geom = StructuredGeometry(
        coord=CoordSystem.CARTESIAN, breakpoints=(0.0, 1.5, 3.0), mat_ids=(1, 0), face_roles=faces,
    )
    return Mesher(geom).partition(CellsByCount.uniform_width(n)).mesh


def _sn_fields(sol, *, adjoint: bool = False) -> dict[str, np.ndarray]:
    out: dict[str, np.ndarray] = {}
    outcome = sol.outcome
    if hasattr(outcome, "keff"):
        try:
            out["k"] = np.asarray(float(outcome.keff))
        except Exception:  # a fixed-source outcome refuses keff
            pass
    flux = sol.importance if adjoint else sol.scalar_flux
    out["scalar_flux"] = np.asarray(flux.values, dtype=float)
    out["angular_flux"] = np.asarray(sol.angular_flux.interior.values, dtype=float)
    return out


#: The one exterior tag the module's own meshes carry (the regression builders carry the regression table).
VACUUM_TAG = 1


def _conditions(mesh):
    from orpheus.transport.boundary_conditions import WallCondition
    from tests.gates._exterior import conditions_for

    return conditions_for(mesh, {VACUUM_TAG: WallCondition.vacuum()})


# ── SN forward eigen ─────────────────────────────────────────────────────────────────────────────────────────────


def _sn_eigen(builder: Callable[[], dict], *, adjoint: bool = False) -> Callable[[], dict[str, np.ndarray]]:
    def solve() -> dict[str, np.ndarray]:
        from orpheus.sn import solve_sn, solve_sn_adjoint

        g = _regression()
        cfg = builder()
        kw: dict[str, Any] = dict(scattering_order=cfg.get("scattering_order", 0), **EIGEN)
        if adjoint:
            sol = solve_sn_adjoint(cfg["materials"], cfg["mesh"], cfg["quadrature"], **kw,
                                   boundary_conditions=g.case_conditions(cfg))
        else:
            sol = solve_sn(cfg["materials"], cfg["mesh"], cfg["quadrature"], **kw,
                           inner_solver=cfg.get("inner_solver", "source_iteration"),
                           boundary_conditions=g.case_conditions(cfg))
        return _sn_fields(sol, adjoint=adjoint)

    return solve


def _sn_fixed(builder: Callable[[], dict], *, scheme: str | None = None) -> Callable[[], dict[str, np.ndarray]]:
    def solve() -> dict[str, np.ndarray]:
        from orpheus.sn import solve_sn_fixed_source
        from orpheus.transport.spatial.scheme import DiscretizationSchemeBase

        g = _regression()
        cfg = builder()
        kw: dict[str, Any] = dict(scattering_order=cfg.get("scattering_order", 0), **FIXED,
                  inner_solver=cfg.get("inner_solver", "source_iteration"))
        if scheme is not None:
            kw["scheme"] = DiscretizationSchemeBase.registry[scheme]()
        sol = solve_sn_fixed_source(cfg["materials"], cfg["mesh"], cfg["quadrature"], cfg["external_source"], **kw,
                                    boundary_conditions=g.case_conditions(cfg))
        return _sn_fields(sol)

    return solve


def _sn_adjoint_fixed(builder: Callable[[], dict], *, scheme: str | None = None) -> Callable[[], dict[str, np.ndarray]]:
    """The importance problem with a detector response in group 1 of the fuel half (a density, R10's input)."""

    def solve() -> dict[str, np.ndarray]:
        from orpheus.sn import solve_sn_adjoint_fixed_source
        from orpheus.transport.spatial.scheme import DiscretizationSchemeBase

        g = _regression()
        cfg = builder()
        q = np.asarray(cfg["external_source"])
        ng, spatial = q.shape[1], q.shape[2:]
        detector = np.zeros((ng, *spatial))
        detector[ng - 1, : max(1, spatial[0] // 2)] = 1.0
        kw: dict[str, Any] = dict(scattering_order=cfg.get("scattering_order", 0), **FIXED)
        if scheme is not None:
            kw["scheme"] = DiscretizationSchemeBase.registry[scheme]()
        sol = solve_sn_adjoint_fixed_source(cfg["materials"], cfg["mesh"], cfg["quadrature"], detector, **kw,
                                            boundary_conditions=g.case_conditions(cfg))
        return _sn_fields(sol, adjoint=True)

    return solve


def _slab_het_fixed_cfg(n: int = 10) -> dict:
    """2G fuel | moderator slab, vacuum | mirror, GL-8, a group-0 source in every cell (per-ordinate density)."""
    from orpheus.geometry.deck import MirrorDeck
    from orpheus.numerics.quadrature import Quadrature

    fuel, mod = _fuel_mod()
    mesh = _slab_het_mesh(n, (VACUUM_TAG, MirrorDeck()))
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    q = np.zeros((quad.n_ordinates, 2, 2 * n))
    q[:, 0] = 1.0 / float(quad.weights.sum())
    return dict(materials={1: fuel, 0: mod}, mesh=mesh, quadrature=quad, external_source=q, scattering_order=1,
                inner_solver="krylov")


def _2d_het_fixed_cfg() -> dict:
    """The regression module's 2-D P1 heterogeneous fixed-source case."""
    return _regression()._cartesian_2d_p1_aniso_het_si("2g")


def _sn_cases() -> list[Case]:
    g = _regression()
    return [
        # forward eigen: 1-D DD, curvilinear DD, 2-D DD
        Case("sn/slab_het_dd/eigen", _sn_eigen(lambda: g._slab_3region("2g", 20))),
        Case("sn/sphere_het_dd/eigen", _sn_eigen(lambda: g._sphere_3region("2g", 20))),
        Case("sn/cyl_het_dd/eigen", _sn_eigen(lambda: g._cylinder_3region("2g", 20, "folded_4x8"))),
        Case("sn/2d_het_dd/eigen", _sn_eigen(lambda: g._cartesian_2d_het_si("2g"))),
        # adjoint eigen (the gauge family)
        Case("sn/slab_het_dd/adjoint_eigen", _sn_eigen(lambda: g._slab_3region("2g", 20), adjoint=True),
             family="adjoint_eigen"),
        Case("sn/2d_het_dd/adjoint_eigen", _sn_eigen(lambda: g._cartesian_2d_het_si("2g"), adjoint=True),
             family="adjoint_eigen"),
        # forward fixed source, P1 anisotropic: slab, sphere, 2-D
        Case("sn/slab_p1_dd/fixed", _sn_fixed(lambda: g._slab_p1_aniso("2g", 20))),
        Case("sn/sphere_p1_dd/fixed", _sn_fixed(lambda: g._sphere_p1_aniso("2g", 20))),
        Case("sn/2d_p1_het_dd/fixed", _sn_fixed(_2d_het_fixed_cfg)),
        # leakage-dominated: 1G vacuum slab
        Case("sn/slab_vacuum_dd/fixed", _sn_fixed(lambda: g._slab_fixed_source("1g", 20))),
        # LD: 1-D heterogeneous P1, 2-D heterogeneous
        Case("sn/slab_het_ld/fixed", _sn_fixed(_slab_het_fixed_cfg, scheme="linear_discontinuous")),
        Case("sn/2d_het_ld/fixed", _sn_fixed(_2d_het_fixed_cfg, scheme="linear_discontinuous")),
        # adjoint fixed source (the detector door, R10): DD slab, LD slab, DD 2-D
        Case("sn/slab_het_dd/adjoint_fixed", _sn_adjoint_fixed(_slab_het_fixed_cfg)),
        Case("sn/slab_het_ld/adjoint_fixed", _sn_adjoint_fixed(_slab_het_fixed_cfg, scheme="linear_discontinuous")),
        Case("sn/2d_p1_het_dd/adjoint_fixed", _sn_adjoint_fixed(_2d_het_fixed_cfg)),
    ]


# ── diffusion and homogeneous ────────────────────────────────────────────────────────────────────────────────────


def _diffusion_mesh(coord: str):
    from orpheus.geometry import CoordSystem, StructuredGeometry
    from orpheus.mesh import CellsByCount, Mesher

    geom = StructuredGeometry(coord=getattr(CoordSystem, coord), breakpoints=(0.0, 1.5, 3.0), mat_ids=(1, 0),
                              face_roles=_diffusion_faces(coord))
    return Mesher(geom).partition(CellsByCount.uniform_width(6)).mesh


def _diffusion_faces(coord: str):
    """Slab: mirror at the left face, vacuum at the right. Solid sphere: the centre carries no role, vacuum outside."""
    from orpheus.geometry.deck import MirrorDeck

    return (MirrorDeck(), VACUUM_TAG) if coord == "CARTESIAN" else (VACUUM_TAG,)


def _diffusion_eigen(coord: str) -> Callable[[], dict[str, np.ndarray]]:
    def solve() -> dict[str, np.ndarray]:
        from orpheus.diffusion.solver import solve_diffusion_1d

        fuel, mod = _fuel_mod()
        mesh = _diffusion_mesh(coord)
        r = solve_diffusion_1d({1: fuel, 0: mod}, mesh, boundary_conditions=_conditions(mesh))
        return {"k": np.asarray(float(r.keff)), "scalar_flux": np.asarray(r.flux.interior.values, dtype=float)}

    return solve


def _diffusion_fixed() -> dict[str, np.ndarray]:
    from orpheus.diffusion.augmented_mesh import DiffusionMesh
    from orpheus.diffusion.solver import solve_diffusion_fixed_source
    from orpheus.transport.mesh.material_mesh import MaterialMesh

    fuel, mod = _fuel_mod()
    mesh = _diffusion_mesh("CARTESIAN")
    mm = MaterialMesh(mesh, {1: fuel, 0: mod})
    dm = DiffusionMesh.from_material_mesh(mm, boundary_conditions=_conditions(mesh))
    q = np.zeros((2, 12))
    q[0] = 1.0
    r = solve_diffusion_fixed_source(dm, q)
    return {"scalar_flux": np.asarray(r.flux.interior.values, dtype=float)}


def _homogeneous() -> dict[str, np.ndarray]:
    from orpheus.homogeneous.solver import solve_homogeneous_infinite

    fuel, _ = _fuel_mod()
    r = solve_homogeneous_infinite(fuel)
    return {"k": np.asarray(float(r.k_inf)), "scalar_flux": np.asarray(r.flux, dtype=float)}


def cases() -> list[Case]:
    return [
        *_sn_cases(),
        Case("diffusion/slab_het/eigen", _diffusion_eigen("CARTESIAN")),
        Case("diffusion/sphere_het/eigen", _diffusion_eigen("SPHERICAL")),
        Case("diffusion/slab_het/fixed", _diffusion_fixed),
        Case("homogeneous/fuel_2g/eigen", _homogeneous),
    ]


# ── the declared principled-equivalence bands (vv-principles, "Bit-identity vs principled-equivalence") ──────────
#
# Rows Step A is PREDICTED to move off bit-identity, each with its band and its reason. The three criteria:
# (1) principled: the moved arithmetic is the integrated chart itself (a named quantity, the source in a cell);
# (2) independent reference: each family keeps its own (2-D: the regression module's closed-form and balance
#     corroboration; LD: the slab MMS of ``ld_moments/memo.md`` and the derivations' Galerkin moment; diffusion: the
#     carry ``keff_diffusion/*`` rows and the continuous-reference gates);
# (3) the drift is FP re-association: iterative rows are bounded by 10 × the solver's own convergence tolerance (the
#     trajectory's rounding changes, the fixed point does not); a direct solve by its reduction depth × condition × eps.
# Every other row stays ``np.array_equal``. A move beyond its band is an algorithmic change, never a re-band.

#: ``case name prefix → (field → rtol, reason)``; rtol is relative to ``max|expected|``.
PRINCIPLED: dict[str, tuple[dict[str, float], str]] = {
    "sn/2d_": ({"k": 10 * EIGEN["keff_tol"], "scalar_flux": 10 * EIGEN["flux_tol"], "angular_flux": 10 * EIGEN["flux_tol"]},
               "multi-D DAG kernel (census C3): the per-volume cell balance gains ÷V on an integrated source, "
               "(QV)/V != Q in floating point; the inner and outer trajectories re-round. Fixed-source rows are "
               "bounded tighter by FIXED_BAND below."),
    "sn/slab_het_ld/": ({"scalar_flux": 0.0, "angular_flux": 0.0},
                        "LD: the 1-D closed form drops θ from its source half (Q-D) and the scan's ×V moves to the "
                        "door, re-associating θ·(Q̂·V); FIXED_BAND applies."),
    "diffusion/": ({"k": 0.0, "scalar_flux": 0.0},
                   "diffusion: the loss output stored integrated (×V retires the /V of operators.py:418) and k "
                   "restated as a pairing (solver.py:431): (x/V)·V != x; a direct LU solve, DIRECT_BAND applies."),
}
#: Fixed-source SN rows (inner_tol 1e-12): 10 × the inner tolerance.
FIXED_BAND = 10 * FIXED["inner_tol"]
#: Diffusion, direct: reduction depth (the 28-row flat loss) × a condition bound 1e3 (`[R]`, measure at landing) × eps.
DIRECT_BAND = 28 * 1e3 * float(np.finfo(float).eps)


def band(name: str, field: str) -> tuple[float, str] | None:
    """The declared band of ``<name>/<field>``, or ``None`` (bit-identity required)."""
    for prefix, (fields, reason) in PRINCIPLED.items():
        if name.startswith(prefix) and field in fields:
            rtol = fields[field]
            if prefix == "diffusion/":
                rtol = DIRECT_BAND
            elif name.endswith("/fixed") or name.endswith("/adjoint_fixed"):
                rtol = FIXED_BAND
            return rtol, reason
    return None


# ── the carried value of a key: the capture, unless a declared re-baseline overrides it ─────────────────────────


@cache
def head() -> Any:
    if not FIXTURE.exists():
        raise FileNotFoundError(f"{FIXTURE} is missing: run _capture_true_ends_acceptance.py on the pre-carve tree")
    return np.load(FIXTURE)


@cache
def rebaseline_events() -> tuple[dict[str, Any], ...]:
    if not REBASELINES.exists():
        return ()
    return tuple(json.loads(REBASELINES.read_text())["events"])


def expected(key: str) -> np.ndarray:
    """The newest declared re-baseline of ``key``, else the capture (the one reading of the fixture, X4)."""
    for event in reversed(rebaseline_events()):
        if key in event["entries"]:
            return np.asarray(event["entries"][key]["new"], dtype=float)
    return head()[key]
