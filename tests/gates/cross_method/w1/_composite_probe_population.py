"""The population builder of the composite record (N05, spec 6.6 "Bit identity", ON05-11 and ON05-22).

Not a test module. Every function RETURNS data and none asserts (``python -O`` strips an ``assert`` outside a collected module).

A composite space is a direct sum: the full field (bulk and trace), the coupled space of the two radial systems, and System B's ray
instance. The builder reaches the 33 composites of the 15 harness fixtures through the problems' public getters
(``full_field_space``, ``system.space``, ``radial_characteristic_field_space`` of an SN problem, ``full_field_space`` of a diffusion
mesh) and names neither class that N05 retires, so that it runs unedited on the unchanged base and on the tree after the unit (a name
import of a retiring class would break it at the unit's commit). It calls only verbs both trees share: ``apply_metric``,
``apply_inverse_metric``, ``inner_product``, ``norm``, ``rank``, ``kernel_basis``, ``metric_diagonal``, ``n_coordinates``,
``member``, ``dual`` and the carriers' ``to_flat``.

Three answers of the sum change by ruling (``unit``, ``integral`` and ``has_coordinate_cone``, ON05-10 and ON05-25); the record holds
no key for them, and ``record_keys_are_admissible`` is the predicate the gate reads.

Keys are ``(fixture, quantity)`` pairs; the arrays are the flat C-order bytes of what each verb returns.
"""

from __future__ import annotations

import numpy as np

from tests.gates.cross_method.w1 import _fixtures as fx, _snapshot as snap

Key = tuple[str, str]

#: Verb names the record must never hold (answers that change by ruling, ON05-10, ON05-25).
FORBIDDEN_VERBS: tuple[str, ...] = ("unit", "integral", "has_coordinate_cone")

#: A second seed for the second operand of ``inner_product``.
_SECOND_SEED = 7


def _flat(array) -> np.ndarray:
    return np.ascontiguousarray(np.asarray(array)).reshape(-1)


def _payload(n: int) -> np.ndarray:
    """A two-column payload chart ``(n, 2)``: the probe vector and its reverse scaled by 3 (two different columns)."""
    v = snap._probe_vector(n)
    return np.ascontiguousarray(np.stack([v, 3.0 * v[::-1]], axis=1))


def _try(out: dict[str, np.ndarray], name: str, fn) -> None:
    """Record ``fn()`` under ``name``; a verb the tree refuses on this input leaves no key (the key-set comparison then reads it)."""
    try:
        out[name] = _flat(fn())
    except (TypeError, ValueError, AttributeError, NotImplementedError, RuntimeError, KeyError):
        return


def arrays_of(space, *, typed: bool = True) -> dict[str, np.ndarray]:
    """The verbs of ``space`` (primal and dual) in both charts, the member round trip, ``norm`` and the two-column payload chart.

    ``typed=False`` is for a space that holds no zero factory (the lifted full field): the typed-carrier rows are then not taken.
    """
    out: dict[str, np.ndarray] = {}
    n = int(space.n_coordinates)
    v = snap._probe_vector(n)
    x = space.member(v) if typed else None
    for tag, sp in (("primal", space), ("dual", space.dual())):
        for verb in ("apply_metric", "apply_inverse_metric"):
            out[f"{tag}.{verb}.flat"] = _flat(getattr(sp, verb)(v))
            if typed:
                out[f"{tag}.{verb}.carrier_to_flat"] = _flat(getattr(sp, verb)(x).to_flat())
            _try(out, f"{tag}.{verb}.payload", lambda sp=sp, verb=verb: getattr(sp, verb)(_payload(n)))
        out[f"{tag}.rank"] = np.array([sp.rank])
        out[f"{tag}.n"] = np.array([sp.n_coordinates])
        out[f"{tag}.kernel_basis"] = _flat(sp.kernel_basis())
        out[f"{tag}.kernel_basis.shape"] = np.array(np.shape(sp.kernel_basis()))
        diagonal = sp.metric_diagonal()
        if diagonal is not None:
            out[f"{tag}.metric_diagonal"] = _flat(diagonal)
    if typed:
        y = space.member(np.random.default_rng(_SECOND_SEED).standard_normal(n))
        out["primal.inner_product"] = np.array([space.inner_product(x, y)])
        out["dual.inner_product"] = np.array([space.dual().inner_product(x, y)])
        _try(out, "primal.norm", lambda: np.array([space.norm(x)]))
        out["primal.member_roundtrip"] = _flat(space.member(v).to_flat())
    return out


def _composites_of_sn(kind: str, ng: int) -> dict[str, object]:
    problem, _ = fx.sn_problem(kind, ng)
    name = snap.sn_fixture_name(kind, ng)
    found: dict[str, object] = {f"{name}.full_field_space": problem.full_field_space, f"{name}.system.space": problem.system.space}
    if problem.radial_characteristic_field_space is not None:
        found[f"{name}.radial_field_space"] = problem.radial_characteristic_field_space
    return found


def _lifted_full_field(problem):
    """The full field whose bulk is the harmonic-moment space, the codomain of ``BulkAnalysisOperator`` (``angular_lift.py`` re-binds on it)."""
    from orpheus.sn.operators.windowing import BulkAnalysisOperator

    scattering = problem.system.factors.scattering
    return BulkAnalysisOperator(scattering.flux_analysis, problem.full_field_space).codomain


def _timed_pair(problem):
    """``(space, timed)``: a timed full field over the SN full field with ``history_depth`` 2 and ONE frame in its history.

    The frame is the probe iterate advanced once: ``advance`` pushes the current frame into the history and installs the new one.
    """
    from orpheus.transport.timed_full_field import TimedFullField

    space = problem.full_field_space
    n = int(space.n_coordinates)
    first = space.member(snap._probe_vector(n))
    second = space.member(2.0 * snap._probe_vector(n)[::-1])
    timed = TimedFullField(interior=first.interior, boundary=first.boundary, _history=(), history_depth=2)
    return space, timed.advance(second.interior, second.boundary)


def _timed_arrays(problem) -> dict[str, np.ndarray]:
    space, timed = _timed_pair(problem)
    out: dict[str, np.ndarray] = {}
    n = int(space.n_coordinates)
    out["timed.history_length_before"] = np.array([timed.history_length])
    out["timed.history_depth_before"] = np.array([timed.history_depth])
    out["timed.to_flat"] = _flat(timed.to_flat())
    applied = space.apply_metric(timed)
    out["timed.apply_metric.to_flat"] = _flat(applied.to_flat())
    out["timed.apply_metric.history_length"] = np.array([applied.history_length])
    out["timed.apply_metric.history_depth"] = np.array([applied.history_depth])
    inverse = space.apply_inverse_metric(timed)
    out["timed.apply_inverse_metric.to_flat"] = _flat(inverse.to_flat())
    out["timed.apply_inverse_metric.history_length"] = np.array([inverse.history_length])
    out["timed.inner_product"] = np.array([space.inner_product(timed, timed)])
    out["timed.n"] = np.array([n])
    return out


def _kernel_arrays(kernel, prefix: str) -> dict[str, np.ndarray]:
    """The summands of a predicted kernel (``right`` and ``left``): the count, and each summand's support and generators as flat arrays."""
    out: dict[str, np.ndarray] = {}
    for side in ("right", "left"):
        summands = getattr(kernel, side)
        out[f"{prefix}.{side}.count"] = np.array([len(summands)])
        for i, s in enumerate(summands):
            out[f"{prefix}.{side}.{i}.support"] = _flat(s.support)
            out[f"{prefix}.{side}.{i}.generators"] = _flat(s.generators)
    return out


def _gauge_arrays(problem) -> dict[str, np.ndarray]:
    """The kernels the loss-kernel gauge predicts through its hand slices of the bulk size, where the problem can call them (the key set records which can)."""
    from orpheus.sn.operators.loss_kernel_gauge import predicted_discretization_kernel, predicted_physical_kernel

    out: dict[str, np.ndarray] = {}
    for name, fn in (("discretization", predicted_discretization_kernel), ("physical", predicted_physical_kernel)):
        try:
            kernel = fn(problem)
        except (TypeError, ValueError, AttributeError, NotImplementedError, RuntimeError, KeyError):
            continue
        out.update(_kernel_arrays(kernel, f"gauge.{name}"))
    return out


def gauge_problems() -> dict[str, object]:
    """Two SN problems on which the gauge's predictors answer non-trivially: the harness's 12 SN fixtures all predict an EMPTY kernel (``[M]`` 0 summands at the base tip), because each has a wall that is not
    all-reflective on a Cartesian mesh or is curvilinear, so these are built here through the public constructors, after ``tests/gates/sn/_singular_loss_box.py`` (``diamond_singular_box``,
    ``diamond_scattering_box``): the harness's 2-cell by 3-cell non-dyadic 2-D box with all four faces specular, ``level_symmetric(sn_order=4)`` (no tangential ordinate), 2 groups; the first with a pure absorber
    (the diamond-difference face modes), the second with a pure within-group scatterer with no capture (one absorption-free closed class per group: the physical modes too)."""
    from orpheus.data.macro_xs.synthetic import make_mixture
    from orpheus.mesh import Mesh2D
    from orpheus.numerics.quadrature import Quadrature
    from orpheus.sn.problem import SNProblem
    from orpheus.transport.boundary_conditions import WallCondition

    total = np.array([0.8, 1.6])
    zero = np.zeros(2)
    mixtures = {
        "absorber": make_mixture(sig_t=total, sig_c=total.copy(), sig_f=zero, nu=zero, chi=zero, sig_s=np.zeros((2, 2))),
        "scatterer": make_mixture(sig_t=total, sig_c=zero, sig_f=zero, nu=zero, chi=zero, sig_s=np.diag(total)),
    }
    out: dict[str, object] = {}
    for name, mixture in mixtures.items():
        nx, ny = fx.BOX_CELLS[0], 3
        mesh = Mesh2D(
            edges_x=np.linspace(0.0, fx.BOX_EXTENTS[0], nx + 1), edges_y=np.linspace(0.0, fx.BOX_EXTENTS[1], ny + 1), mat_map=np.zeros((nx, ny), dtype=int),
            face_roles={"xmin": 0, "xmax": 0, "ymin": 0, "ymax": 0})
        out[f"gauge_box_{name}"] = SNProblem(
            mesh, Quadrature.level_symmetric(sn_order=4), {0: mixture}, scheme=None, scattering_order=0,
            boundary_conditions=fx.conditions_for(mesh, {0: WallCondition.specular(fx.SPECULAR_ALBEDO)}))
    return out


def capture_composite_probes() -> dict[Key, np.ndarray]:
    """``{(fixture, quantity): flat array}`` over the 33 composites of the harness and the ON05-22 additions."""
    arrays: dict[Key, np.ndarray] = {}

    def take(prefix: str, space, *, typed: bool = True) -> None:
        for quantity, array in arrays_of(space, typed=typed).items():
            arrays[(prefix, quantity)] = array

    for kind in fx.SN_KINDS:
        for ng in fx.GROUPS:
            for label, space in _composites_of_sn(kind, ng).items():
                take(label, space)
            problem, _ = fx.sn_problem(kind, ng)
            name = snap.sn_fixture_name(kind, ng)
            for quantity, array in _gauge_arrays(problem).items():
                arrays[(f"{name}.loss_kernel_gauge", quantity)] = array
            if kind in ("slab", "box"):
                take(f"{name}.lifted_full_field", _lifted_full_field(problem), typed=False)
                for quantity, array in _timed_arrays(problem).items():
                    arrays[(f"{name}.timed_full_field", quantity)] = array
    for name, problem in gauge_problems().items():
        for quantity, array in _gauge_arrays(problem).items():
            arrays[(f"{name}.loss_kernel_gauge", quantity)] = array
    for ng in fx.GROUPS:
        mesh, solver, _ = fx.diffusion_problem(ng)
        name = snap.diffusion_fixture_name(ng)
        take(f"{name}.full_field_space", mesh.full_field_space)
        for quantity, array in _kernel_arrays(solver.predicted_kernel, "predicted_kernel").items():
            arrays[(f"{name}.solver", quantity)] = array
        arrays[(f"{name}.solver", "group_of_each_coordinate")] = _flat(solver._group_of_each_coordinate())
    return {k: np.ascontiguousarray(a) for k, a in arrays.items()}


def record_keys_are_admissible(keys) -> list[Key]:
    """The keys that name a verb whose answer changes by ruling (``unit``, ``integral``, ``has_coordinate_cone``); an admissible record returns ``[]``."""
    return sorted(k for k in keys if any(part in FORBIDDEN_VERBS for part in k[1].split(".")))


def write_record(directory) -> dict[Key, np.ndarray]:
    """Take the record: ``composite_probes.npz`` (keys ``fixture::quantity``), ``MANIFEST.txt`` and ``STAMP.txt`` under ``directory``.

    Run once, on the unchanged base, by the commit that adds the record; no later unit re-takes it.
    """
    from pathlib import Path

    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    arrays = capture_composite_probes()
    np.savez(directory / "composite_probes.npz", **{f"{f}::{n}": a for (f, n), a in arrays.items()})
    lines = [f"{f}\t{n}\t{tuple(a.shape)}\t{a.dtype}\t{snap.sha1(a)}" for (f, n), a in sorted(arrays.items())]
    (directory / "MANIFEST.txt").write_text("\n".join(lines) + "\n")
    (directory / "STAMP.txt").write_text("\t".join(snap.host_stamp()) + "\n")
    return arrays
