r"""The chart law of #707 Step A, as DATA: every data-side value is the spatial measure times today's density value.

Not a test module. Two readers:

* ``_capture_source_chart_anchor.py`` freezes, on the pre-carve tree, the DENSITY-chart output of every flux → source
  atom applied to one seeded state, and the stored value at each input door, into ``source_chart_anchor.npz``;
* ``test_source_chart_law.py`` recomputes the same quantities on the current tree and requires
  ``new == m ⊙ anchor``, where ``m`` is the spatial measure built HERE from raw mesh data (never from a production
  metric object).

The law (R6, R8, Q-D). The data side is V* = A ⊗ (E ⊗ X)*: integrated over the cell (and its in-cell trial basis),
per ordinate and per steradian, with no quadrature weight. So the measure ``m`` that maps today's density chart to the
integrated chart is

* SN bulk: ``m[n, g, i, o] = V_i · θ^{|o|}`` (θ = 1/3 per axis, LM-1989 Eq. 4.3b, typed here; ``|o|`` the number of
  slope indices of the tensor-Legendre moment ``o``); ``o`` is absent on a slopeless scheme. No ``w_n``.
* SN trace: ``m = 1`` IN STEP A (declared): the trace re-pose (the quotient by the tangential ordinates, the current
  view, R4/R9's ``|Ω·n|·w_n·A_f``) is Step B, and a Step-A carve that integrates the trace is ahead of its ruling.
* diffusion bulk: ``m[g, i] = V_i``; diffusion trace ``m = 1`` (method harmonisation, later).

Population (X2): 5 SN problems (non-uniform slab DD P1, the same with LD, a two-region sphere DD, a non-uniform 2-D
DD P1, the same 2-D with LD) × the atoms the pencil walk finds (``StreamingCollisionOperator``,
``ScatteringOperator``, ``N2NOperator``, ``SNBoundaryOperator``, ``FissionOperator``); 2 SN doors (the external-source
door flat and moment-resolved, the detector door); diffusion 1 problem × its loss and fission. ``N2NOperator`` reads a
zero ``Sig2`` on every library mixture, so its rows are VACUOUS (0 == 0 ⊙ m) and are kept only as a shape witness:
they are not counted as coverage. The radial-characteristic atoms (System B on the sphere) are not reached by the
full-field walk and are owed (spec row CL-RC).
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any, Callable

import numpy as np

ANCHOR = Path(__file__).with_name("source_chart_anchor.npz")
THETA = 1.0 / 3.0  # LM-1989 Eq. 4.3b, typed (never read from the scheme: the scheme and its metric share one basis)
SEED = 20261112

#: The atom classes a pencil walk records (by class name; a rename at the carve is a loud KeyError, not a pass).
ATOMS = ("StreamingCollisionOperator", "ScatteringOperator", "N2NOperator", "SNBoundaryOperator", "FissionOperator")
VACUOUS = ("N2NOperator",)


# ── the problems ─────────────────────────────────────────────────────────────────────────────────────────────────


def _mixtures():
    from orpheus.data.macro_xs.synthetic import get_mixture

    return {1: get_mixture("A", "2g"), 0: get_mixture("B", "2g")}


def _conditions(mesh):
    from orpheus.transport.boundary_conditions import WallCondition
    from tests.gates._exterior import conditions_for

    return conditions_for(mesh, {1: WallCondition.vacuum()})


def _slab_mesh():
    """Fuel 1.0 cm | moderator 2.0 cm, 5 cells each: widths 0.2 and 0.4 (V != 1, V non-uniform)."""
    from orpheus.geometry import CoordSystem, StructuredGeometry
    from orpheus.geometry.deck import MirrorDeck
    from orpheus.mesh import CellsByCount, Mesher

    geom = StructuredGeometry(coord=CoordSystem.CARTESIAN, breakpoints=(0.0, 1.0, 3.0), mat_ids=(1, 0),
                              face_roles=(1, MirrorDeck()))
    return Mesher(geom).partition(CellsByCount.uniform_width(5)).mesh


def _sphere_mesh():
    from orpheus.geometry import CoordSystem, StructuredGeometry
    from orpheus.mesh import CellsByCount, Mesher

    geom = StructuredGeometry(coord=CoordSystem.SPHERICAL, breakpoints=(0.0, 1.0, 2.0), mat_ids=(1, 0),
                              face_roles=(1,))
    return Mesher(geom).partition(CellsByCount.uniform_width(4)).mesh


def _mesh_2d():
    """2-D fuel | moderator, non-uniform x edges, vacuum on xmax, mirrors elsewhere."""
    from orpheus.geometry.deck import MirrorDeck
    from orpheus.mesh import Mesh2D

    edges_x = np.array([0.0, 0.3, 0.7, 1.2, 1.8, 2.5])
    edges_y = np.array([0.0, 0.4, 1.0, 1.5])
    mat = np.zeros((5, 3), dtype=int)
    mat[:2, :] = 1
    return Mesh2D(edges_x=edges_x, edges_y=edges_y, mat_map=mat,
                  face_roles={"xmin": MirrorDeck(), "xmax": 1, "ymin": MirrorDeck(), "ymax": MirrorDeck()})


def problems() -> dict[str, Callable[[], Any]]:
    from orpheus.numerics.quadrature import Quadrature
    from orpheus.sn.problem import SNProblem
    from orpheus.transport.spatial import LinearDiscontinuous

    def sn(mesh_fn, quad, order, ld=False):
        def build():
            mesh = mesh_fn()
            return SNProblem(mesh, quad(), _mixtures(), scheme=LinearDiscontinuous() if ld else None,
                             scattering_order=order, boundary_conditions=_conditions(mesh))
        return build

    gl8 = lambda: Quadrature.gauss_legendre(n_ordinates=8)  # noqa: E731
    ls4 = lambda: Quadrature.level_symmetric(sn_order=4)  # noqa: E731
    return {
        "slab_dd_p1": sn(_slab_mesh, gl8, 1),
        "slab_ld_p1": sn(_slab_mesh, gl8, 1, ld=True),
        "sphere_dd": sn(_sphere_mesh, gl8, 0),
        "2d_dd_p1": sn(_mesh_2d, ls4, 1),
        "2d_ld": sn(_mesh_2d, ls4, 0, ld=True),
    }


# ── the state, the atoms, the measure ────────────────────────────────────────────────────────────────────────────


def leaves(op: Any) -> list[Any]:
    """Every operator reachable from ``op`` through ``blocks`` / ``_terms`` / ``inner``, depth first, deduplicated."""
    out: list[Any] = []
    seen: set[int] = set()

    def walk(x: Any) -> None:
        if x is None or id(x) in seen:
            return
        if isinstance(x, (list, tuple)):
            for y in x:
                walk(y)
            return
        if not hasattr(x, "apply"):
            return
        seen.add(id(x))
        out.append(x)
        for attr in ("blocks", "_terms", "inner"):
            walk(getattr(x, attr, None))

    walk(op)
    return out


def atoms(problem: Any) -> dict[str, Any]:
    found: dict[str, Any] = {}
    for leaf in leaves(problem.pencil.lhs) + leaves(problem.pencil.rhs):
        name = type(leaf).__name__
        if name in ATOMS:
            found.setdefault(name, leaf)
    return found


def state(problem: Any, rng: np.random.Generator) -> Any:
    """A seeded flux composite on the problem's full field (bulk AND trace random, every term activated)."""
    from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
    from orpheus.transport.fields.angular_flux import AngularFlux
    from orpheus.transport.timed_full_field import TimedFullField

    psi = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, space=problem.full_field_space)
    psi.interior.values[...] = rng.standard_normal(psi.interior.values.shape)
    for face in psi.boundary.layout.faces:
        fv = psi.boundary.face_view(face)
        fv[...] = rng.standard_normal(fv.shape)
    return psi


def state_from(problem: Any, interior: np.ndarray, boundary: np.ndarray) -> Any:
    psi = state(problem, np.random.default_rng(0))
    psi.interior.values[...] = interior
    psi.boundary.values[...] = boundary
    return psi


def sn_bulk_measure(problem: Any, shape: tuple[int, ...]) -> np.ndarray:
    """``V_i · θ^{|o|}`` broadcast to the bulk ``shape`` ``(N, ng, *spatial[, 2^d])``; no ``w_n`` (R8)."""
    V = np.asarray(problem.volumes, dtype=float)
    m = V.reshape((1, 1) + V.shape)
    if len(shape) == 2 + V.ndim + 1:  # a trailing in-cell moment axis
        mm = np.array([1.0])
        for _ in range(V.ndim):
            mm = np.kron(mm, np.array([1.0, THETA]))
        if mm.size != shape[-1]:
            raise ValueError(f"moment axis {shape[-1]} != 2^d = {mm.size}")
        m = m[..., None] * mm
    return np.broadcast_to(m, shape)


def diffusion_problem() -> Any:
    from orpheus.diffusion.augmented_mesh import DiffusionMesh
    from orpheus.diffusion.solver import DiffusionSolver
    from orpheus.transport.mesh.material_mesh import MaterialMesh

    mesh = _slab_mesh()
    dm = DiffusionMesh.from_material_mesh(MaterialMesh(mesh, _mixtures()), boundary_conditions=_conditions(mesh))
    return DiffusionSolver(dm)


def diffusion_volumes(solver: Any) -> np.ndarray:
    return np.asarray(solver.mesh.volumes, dtype=float)


# ── the quantities the anchor freezes and the gate recomputes ───────────────────────────────────────────────────


def _bulk_and_trace(member: Any) -> tuple[np.ndarray, np.ndarray]:
    return np.array(member.interior.values, dtype=float), np.array(member.boundary.values, dtype=float)


def sn_atom_outputs(name: str, problem: Any, psi: Any) -> dict[str, np.ndarray]:
    """Each atom's output on ``psi``; an atom whose ``apply`` RAISES today is recorded under ``<name>/raised``.

    `[M]` 2026-11-12: ``FissionOperator.apply`` raises on both LD problems (``FissionMaterialField gather: the layout
    indexes 1 spatial axes but the requested shape (10, 2) has 2``): the in-cell moment axis is read as a spatial axis.
    A pre-existing defect, outside Step A, reported in the spec (row CL-LDF).
    """
    out: dict[str, np.ndarray] = {}
    raised: list[str] = []
    for atom_name, op in atoms(problem).items():
        try:
            result = op.apply(psi)
        except ValueError:
            raised.append(atom_name)
            continue
        bulk, trace = _bulk_and_trace(result)
        out[f"{name}/{atom_name}/bulk"] = bulk
        out[f"{name}/{atom_name}/trace"] = trace
    out[f"{name}/raised"] = np.array(sorted(raised), dtype=str)
    return out


def _first_system(source: Any) -> Any:
    """The full-field member of a (possibly coupled, carrying-mesh) source."""
    return source.systems[0] if hasattr(source, "systems") else source


def door_sources(problem: Any) -> dict[str, np.ndarray]:
    """The external-source door's stored bulk for a seeded flat density (and, on LD, a moment-resolved one)."""
    from orpheus.sn.solver import _build_fixed_source_rhs

    rng = np.random.default_rng(SEED + 1)
    N, ng, spatial = problem.quad.n_ordinates, problem.ng, tuple(problem.spatial_shape)
    flat = rng.uniform(0.5, 1.5, size=(N, ng, *spatial))
    out = {"flat_input": flat, "flat_stored": np.array(_first_system(_build_fixed_source_rhs(flat, problem)).interior.values)}
    if problem.scheme.spatial_basis_per_axis > 1:
        resolved = rng.uniform(-0.5, 1.5, size=(N, ng, *spatial, 2 ** len(spatial)))
        out["resolved_input"] = resolved
        out["resolved_stored"] = np.array(_first_system(_build_fixed_source_rhs(resolved, problem)).interior.values)
    return out


def detector_source(problem: Any, materials: Any, mesh: Any) -> tuple[np.ndarray, np.ndarray]:
    """The detector door: the stored adjoint source ``q*`` the posing receives for a seeded density ``Σ_d``.

    Captured by wrapping ``SourcePosing.of_pencil`` in ``orpheus.sn.solver``'s namespace for one solve.
    """
    import orpheus.sn.solver as solver_mod
    from orpheus.sn.solver import solve_sn_adjoint_fixed_source

    rng = np.random.default_rng(SEED + 2)
    sigma_d = rng.uniform(0.5, 1.5, size=(problem.ng, *problem.spatial_shape))
    seen: list[Any] = []
    original = solver_mod.SourcePosing.of_pencil

    def spy(pencil, sigma, source, *a, **k):
        seen.append(source)
        return original(pencil, sigma, source, *a, **k)

    solver_mod.SourcePosing.of_pencil = spy  # type: ignore[method-assign]
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            solve_sn_adjoint_fixed_source(materials, mesh, problem.quad, sigma_d, max_inner=3,
                                          inner_tol=1e-3, boundary_conditions=_conditions(mesh))
    except Exception:  # a truncated solve may refuse its exit; the source is already captured
        pass
    finally:
        solver_mod.SourcePosing.of_pencil = original  # type: ignore[method-assign]
    if not seen:
        raise RuntimeError("the detector door was not reached: SourcePosing.of_pencil never called")
    member = _first_system(seen[0])
    return sigma_d, np.array(member.interior.values, dtype=float)


def diffusion_outputs(solver: Any) -> dict[str, np.ndarray]:
    rng = np.random.default_rng(SEED + 3)
    phi = solver.template
    phi_vals = rng.uniform(0.5, 1.5, size=phi.to_flat().shape)
    member = type(phi).from_flat(phi_vals, phi)
    out = {"diffusion/state": phi_vals}
    for label, op in (("loss", solver.loss), ("fission", solver.fission)):
        y = op.apply(member)
        out[f"diffusion/{label}/bulk"] = np.array(y.interior.values, dtype=float)
        out[f"diffusion/{label}/trace"] = np.array(y.boundary.values, dtype=float)
    return out


def everything() -> dict[str, np.ndarray]:
    """Every anchored array, keyed; the seeded states are stored beside the outputs they produced."""
    payload: dict[str, np.ndarray] = {}
    for i, (name, build) in enumerate(problems().items()):
        problem = build()
        psi = state(problem, np.random.default_rng(SEED + 10 + i))
        payload[f"{name}/state/bulk"], payload[f"{name}/state/trace"] = _bulk_and_trace(psi)
        payload.update(sn_atom_outputs(name, problem, psi))
        for k, v in door_sources(problem).items():
            payload[f"{name}/door/{k}"] = v
    slab = problems()["slab_dd_p1"]()
    sd, stored = detector_source(slab, _mixtures(), _slab_mesh())
    payload["slab_dd_p1/detector/input"], payload["slab_dd_p1/detector/stored"] = sd, stored
    payload.update(diffusion_outputs(diffusion_problem()))
    return payload
