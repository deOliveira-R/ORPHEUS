# pyright: basic, reportAttributeAccessIssue=false, reportCallIssue=false, reportArgumentType=false, reportMissingImports=false
r"""Capture the wave-2 carry fixture FROM THE HEAD TREE (``main`` ``5d04c296``). Run it once, from outside the tree.

Not a test module and NOT importable in the carve branch: it spells every problem in the typed-law vocabulary that
slice S2 retires (``orpheus.geometry.boundary``), so it only runs inside a detached worktree of ``5d04c296``::

    git worktree add --detach <tmp> 5d04c296
    cd /tmp && PYTHONPATH=<tmp> CARRY_HEAD_TREE=<tmp> python <this file> <output.npz>
    git worktree remove --force <tmp>

It prints ``orpheus.__file__`` first and REFUSES to run when that file is not inside ``CARRY_HEAD_TREE`` or when the
tree's commit is not :data:`HEAD_COMMIT` (a probe that silently imports the carve branch would freeze the answers of
the code the fixture is meant to be checked against, and the gate would then be green by construction: X4).

What is frozen (the sizes are in ``_carry_cases``): per (body, declaration) the eigenvalue ``keff`` and the scalar
flux of ``solve_sn``, the scalar flux of ``solve_sn_fixed_source`` under a bulk source, the slab's fixed-source flux
under a uniform inflow on the left face, the dense matrix of the realized operator of both faces of an S16 slab for
each slab declaration, and, for each of the seven ``bc_equivalence`` snapshots, the image the HEAD realizer gives of
the frozen probe (the snapshot files themselves are independent derivations and are NOT re-captured: their sha256 is
recorded so the gate can prove they were not re-baselined, section 8 row R-10).
"""

from __future__ import annotations

import hashlib
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

import importlib.util

#: ``_carry_cases`` is read BY PATH from beside this file: the HEAD tree has no such module, and importing it by
#: package name would need the carve branch on ``sys.path`` (which would shadow the HEAD tree's ``orpheus``).
_spec = importlib.util.spec_from_file_location("_carry_cases", Path(__file__).with_name("_carry_cases.py"))
assert _spec is not None and _spec.loader is not None
cases = importlib.util.module_from_spec(_spec)
sys.modules["_carry_cases"] = cases
_spec.loader.exec_module(cases)


def _law(declaration, face):
    from orpheus.geometry import BC
    from orpheus.geometry.boundary import AlbedoBoundary, IsotropicReturn, SpecularReturn

    kind = declaration[0]
    if kind == "vacuum":
        return BC.vacuum
    if kind == "mirror":
        return BC.reflective
    if kind == "specular":
        return AlbedoBoundary(declaration[1], SpecularReturn(axis="x"))
    if kind == "white":
        return AlbedoBoundary(declaration[1], IsotropicReturn(axis="x", outward_sign=-1 if face == "xmin" else +1))
    raise ValueError(declaration)


def _mesh(body, declaration, *, left=None, right=None):
    """The mesh of ``body`` with ``declaration`` on its outer wall(s); ``left``/``right`` override a slab's faces."""
    from orpheus.geometry import BC, CoordSystem, StructuredGeometry
    from orpheus.mesh import CellsByCount, Mesher

    coord = getattr(CoordSystem, body.coord)
    if body.coord == "CARTESIAN":
        boundaries = (left if left is not None else _law(declaration, "xmin"), right if right is not None else _law(declaration, "xmax"))
    elif body.hollow:
        boundaries = (BC.reflective, _law(declaration, "xmax"))
    else:
        boundaries = (_law(declaration, "xmax"),)
    geometry = StructuredGeometry(coord=coord, breakpoints=body.breakpoints, mat_ids=body.mat_ids, boundaries=boundaries)
    return Mesher(geometry).partition(tuple(CellsByCount.uniform_volume(n) for n in body.cells)).mesh


def _phi(solution):
    return np.asarray(solution.scalar_flux.values, dtype=float)


def _capture_solves(out):
    from orpheus.geometry.boundary import ConstantInflowSource, PrescribedInflow
    from orpheus.sn.problem import SNProblem
    from orpheus.sn.solver import solve_sn, solve_sn_fixed_source

    mats = cases.materials()
    for name, declaration in cases.declared_cases():
        body = cases.BODIES[name]
        quad = cases.quadrature(body)
        mesh = _mesh(body, declaration)
        eigen = solve_sn(mats, mesh, quad, **cases.EIGEN)
        out[cases.key(name, declaration, "keff")] = np.float64(eigen.outcome.keff)
        out[cases.key(name, declaration, "phi_eigen")] = _phi(eigen)
        sn = SNProblem(mesh, quad, mats)
        source = cases.bulk_source(quad.n_ordinates, sn.ng, tuple(sn.spatial_shape))
        fixed = solve_sn_fixed_source(mats, mesh, quad, source, **cases.FIXED)
        out[cases.key(name, declaration, "phi_fixed")] = _phi(fixed)
    # the uniform inflow on the left face of the slab (curvilinear inflow is refused, #653)
    body = cases.BODIES["slab"]
    quad = cases.quadrature(body)
    mesh = _mesh(body, None, left=PrescribedInflow(ConstantInflowSource(cases.INFLOW_VALUE)), right=_law(("vacuum",), "xmax"))
    sn = SNProblem(mesh, quad, mats)
    source = np.zeros((quad.n_ordinates, sn.ng, *sn.spatial_shape))
    fixed = solve_sn_fixed_source(mats, mesh, quad, source, **cases.FIXED)
    out[cases.inflow_key("phi_fixed")] = _phi(fixed)


def _other_law(declaration, chart):
    from orpheus.geometry import BC

    kind = declaration[0]
    return {"vacuum": BC.vacuum, "mirror": BC.reflective, "white": BC.white, "zero_flux": BC("zero_flux"),
            "albedo": BC("albedo", {"albedo": declaration[1] if len(declaration) > 1 else 0.0})}[kind]


def _chart_mesh(chart, declaration, *, mat_ids=(1, 0), cells=cases.OTHER_CELLS, left_mirror=False):
    from orpheus.geometry import CoordSystem, StructuredGeometry, BC
    from orpheus.mesh import CellsByCount, Mesher

    coord_name, breakpoints = cases.CHARTS[chart]
    law = _other_law(declaration, chart)
    if chart == "slab":
        boundaries = (BC.reflective if left_mirror else law, law)
    else:
        boundaries = (law,)
    geometry = StructuredGeometry(coord=getattr(CoordSystem, coord_name), breakpoints=breakpoints, mat_ids=mat_ids, boundaries=boundaries)
    return Mesher(geometry).partition(tuple(CellsByCount.uniform_volume(n) for n in cells)).mesh


def _capture_other_methods(out):
    import contextlib
    import io

    from orpheus.cp.solver import CPParams, solve_cp
    from orpheus.diffusion.solver import solve_diffusion_1d
    from orpheus.geometry import BC, CoordSystem, StructuredGeometry
    from orpheus.mc.solver import MCMesh, MCParams, solve_monte_carlo
    from orpheus.mesh import CellsByCount, Mesher
    from orpheus.moc.solver import solve_moc

    mats = cases.materials()
    quiet = io.StringIO()
    for chart in cases.CHARTS:
        for declaration in cases.DIFFUSION_DECLARATIONS:
            out[cases.other_key("diffusion", chart, declaration)] = np.float64(solve_diffusion_1d(mats, _chart_mesh(chart, declaration)).keff)
        for declaration in cases.CP_DECLARATIONS:
            with contextlib.redirect_stdout(quiet):
                result = solve_cp(mats, _chart_mesh(chart, declaration, left_mirror=True), CPParams())
            out[cases.other_key("cp", chart, declaration)] = np.float64(result.keff)
    cylinder = Mesher(StructuredGeometry(coord=CoordSystem.CYLINDRICAL, breakpoints=cases.CHARTS["cylinder"][1], mat_ids=(1, 0), boundaries=(BC.reflective,))).partition(
        tuple(CellsByCount.uniform_volume(n) for n in cases.MOC_CELLS)).mesh
    with contextlib.redirect_stdout(quiet):
        out[cases.other_key("moc", "cylinder", ("mirror",))] = np.float64(solve_moc(mats, cylinder, **cases.MOC_ARGS).keff)
    periodic = Mesher(StructuredGeometry(coord=CoordSystem.CARTESIAN, breakpoints=cases.CHARTS["slab"][1], mat_ids=(2, 0), boundaries=(BC("periodic"), BC("periodic")))).partition(
        tuple(CellsByCount.uniform_volume(n) for n in cases.OTHER_CELLS)).mesh
    mc_mats = {2: mats[1], 0: mats[0]}
    with contextlib.redirect_stdout(quiet):
        out[cases.other_key("mc", "slab", ("periodic",))] = np.float64(solve_monte_carlo(mc_mats, MCParams(**cases.MC_PARAMS, geometry=MCMesh(periodic, pitch=cases.MC_PARAMS["pitch"]))).keff)


def _capture_face_matrices(out):
    from orpheus.sn.problem import SNProblem

    body = cases.BODIES["slab"]
    quad = cases.quadrature(body, ordinates=cases.S16_ORDINATES)
    mats = cases.materials()
    for declaration in cases.DECLARATIONS["slab"]:
        sn = SNProblem(_mesh(body, declaration), quad, mats)
        trace = sn.angular_trace
        for face in ("xmin", "xmax"):
            n_out = len(trace.outflow_indices_for_face(face))
            op = sn.bc[face]
            out[cases.face_key(declaration, face)] = np.stack(
                [np.asarray(op.apply(np.eye(n_out)[j])) for j in range(n_out)], axis=1
            )


def _capture_snapshots(out, tree):
    from orpheus.sn.boundary.realizer import SNBoundaryRealizer
    from tests.gates.geometry._generate_bc_equivalence_snapshots import CASES
    from tests.gates.sn._test_helpers import face_method_space

    snapshots = Path(tree) / "tests/gates/geometry/snapshots"
    for case in CASES:
        path = snapshots / f"bc_equivalence_{case.case_id}.npz"
        snap = np.load(path)
        out[f"snapshot_sha256/{case.case_id}"] = np.array(hashlib.sha256(path.read_bytes()).hexdigest())
        space = face_method_space(case.build_quadrature(), face=case.face, faces=case.faces)
        realizer = SNBoundaryRealizer()
        composed = []

        def realize(law, _realizer=realizer, _space=space, _composed=composed):
            _composed.append(law)
            return _realizer.realize(law, _space)

        operator = case.compose(realize)
        probe = snap["psi_out"]
        if case.partner_face is not None:
            (law,) = composed
            if law.geometry_map.domain_face(case.face) == case.partner_face:
                probe = snap["psi_out_partner"]
        out[f"snapshot_image/{case.case_id}"] = np.asarray(operator.apply(probe))


def main(argv):
    import orpheus

    tree = Path(os.environ["CARRY_HEAD_TREE"]).resolve()
    here = Path(orpheus.__file__).resolve()
    print("orpheus.__file__ =", here)
    if tree not in here.parents:
        raise SystemExit(f"REFUSED: orpheus is imported from {here}, not from the HEAD tree {tree}")
    commit = subprocess.run(["git", "-C", str(tree), "rev-parse", "--short=8", "HEAD"], capture_output=True, text=True, check=True).stdout.strip()
    if commit != cases.HEAD_COMMIT:
        raise SystemExit(f"REFUSED: the tree is at {commit}, not {cases.HEAD_COMMIT}")
    out: dict[str, np.ndarray] = {}
    _capture_solves(out)
    _capture_face_matrices(out)
    _capture_snapshots(out, tree)
    _capture_other_methods(out)
    out["meta/commit"] = np.array(commit)
    out["meta/numpy"] = np.array(np.__version__)
    out["meta/platform"] = np.array(sys.platform)
    np.savez(argv[1], **out)
    print(f"froze {len(out)} arrays at {commit} into {argv[1]}")


if __name__ == "__main__":
    main(sys.argv)
