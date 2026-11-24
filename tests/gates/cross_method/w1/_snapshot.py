r"""The capture, the table of named getters, the compare helper and the snapshot store of the bit-identity harness.

Not a test module (no ``test_`` prefix). Every function here RETURNS data and none asserts: an ``assert`` in a module that a test
imports is removed by ``python -O`` (the canonical command), so a helper that asserted would turn a red gate green. The two test
modules next to this file hold every assertion.

Spaces are reached only through the tables of named getters below (attribute chains off a problem or a mesh) and read only
through the four verbs ``apply_metric``, ``apply_inverse_metric``, ``metric_diagonal`` and ``n_coordinates``. A later unit that
deletes or re-spells an accessor re-points the getter here and leaves the stored bytes alone. No space class is named. Three
classes are, each for a reason no accessor replaces and none on a retirement row of section 5 of
``.claude/plans/ontology_w1_spec.md``: ``FlattenedOperator`` densifies an operator (quantity group a), ``DenseMetric`` builds
the dense metric of a Gram (group c), and ``NotInCoordinateSpace`` is the refusal the function-space arm raises for a flat vector,
on which the verbs fall back to the space's own shape.

Taking the snapshots (a commit of its own, on the host named by the stamp)::

    python -m tests.gates.cross_method.w1._snapshot tests/gates/cross_method/w1/snapshots
"""

from __future__ import annotations

import hashlib
import os
import platform
import sys
from collections.abc import Callable, Mapping
from pathlib import Path
from typing import Any

import numpy as np
import scipy

from orpheus.numerics.coordinate_space import NotInCoordinateSpace
from orpheus.numerics.flat_operator import FlattenedOperator
from orpheus.numerics.metric import DenseMetric

from . import _fixtures as fx

#: The seed of the probe vector ``v``. Each space draws from a fresh generator, so ``v`` depends only on the seed and on
#: ``n_coordinates`` and never on the order of the getters.
SEED = 20261123

SNAPSHOT_DIR = Path(__file__).parent / "snapshots"

Key = tuple[str, str]  # (fixture, quantity name)
Families = list[tuple[str, str, Callable[[], np.ndarray]]]

# ─────────────────────────── the table of named getters ───────────────────────────


Row = tuple[str, Callable[[Any], Any], frozenset[str]]
TRACE, PROBED = "b", "d"  # quantity group b: ``apply_metric(ones)`` and ``metric_diagonal()``; group d: ``apply_metric(v)`` and ``apply_inverse_metric(v)``


def _trace_of_radial_field(problem: Any) -> Any:
    field = problem.radial_characteristic_field_space
    return None if field is None else field.trace_space  # the radial system exists on curvilinear problems only


#: One row per space of the SN problem: ``(name, getter(problem), groups)``. A getter returns ``None`` for a space the problem
#: does not have. A row is one value, so a group cannot name a space that has no getter.
SN_ROWS: tuple[Row, ...] = (
    ("bulk_space", lambda p: p.bulk_space, frozenset({PROBED})),
    ("angular_bulk_space", lambda p: p.angular_bulk_space, frozenset({PROBED})),
    ("angular_trial_space", lambda p: p.angular_trial_space, frozenset({PROBED})),
    ("moment_space(2)", lambda p: p.moment_space(2), frozenset({PROBED})),
    ("full_field_space", lambda p: p.full_field_space, frozenset({PROBED})),
    ("full_field_space.trace_space", lambda p: p.full_field_space.trace_space, frozenset({TRACE})),
    ("radial_characteristic_boundary_space", lambda p: p.radial_characteristic_boundary_space, frozenset({TRACE, PROBED})),
    ("radial_characteristic_interior_space", lambda p: p.radial_characteristic_interior_space, frozenset({PROBED})),
    ("radial_characteristic_field_space", lambda p: p.radial_characteristic_field_space, frozenset({PROBED})),
    ("radial_characteristic_field_space.trace_space", _trace_of_radial_field, frozenset({TRACE})),
)

#: The same for the diffusion slab, whose getters take the ``DiffusionMesh``.
DIFFUSION_ROWS: tuple[Row, ...] = (
    ("bulk_space", lambda m: m.bulk_space, frozenset({PROBED})),
    ("full_field_space", lambda m: m.full_field_space, frozenset({PROBED})),
    ("full_field_space.trace_space", lambda m: m.full_field_space.trace_space, frozenset({TRACE})),
    ("scalar_trace", lambda m: m.scalar_trace, frozenset({TRACE})),
)


# ─────────────────────────── the verbs ───────────────────────────


def _probe_vector(n: int) -> np.ndarray:
    return np.random.default_rng(SEED).standard_normal(n)


def _apply(space: Any, verb: str, x_flat: np.ndarray) -> np.ndarray:
    """Apply a metric verb to ``x_flat``, flat first and in ``space.shape`` when the flat form is refused.

    The function-space arm refuses a flat vector with ``NotInCoordinateSpace``; any other error propagates.
    """
    method = getattr(space, verb)
    try:
        return np.asarray(method(x_flat))
    except NotInCoordinateSpace:
        return np.asarray(method(x_flat.reshape(space.shape)))


def _space_arrays(problem: Any, rows: tuple[Row, ...]) -> dict[str, np.ndarray]:
    arrays: dict[str, np.ndarray] = {}
    for name, getter, groups in rows:
        space = getter(problem)
        if space is None:
            continue
        n = int(space.n_coordinates)
        if TRACE in groups:
            arrays[f"b.{name}.apply_metric_ones"] = _apply(space, "apply_metric", np.ones(n))
            diagonal = space.metric_diagonal()
            if diagonal is not None:
                arrays[f"b.{name}.metric_diagonal"] = np.asarray(diagonal)
        if PROBED in groups:
            v = _probe_vector(n)
            arrays[f"d.{name}.apply_metric"] = _apply(space, "apply_metric", v)
            arrays[f"d.{name}.apply_inverse_metric"] = _apply(space, "apply_inverse_metric", v)
    return arrays


def _dense(operators: Mapping[str, tuple[Any, Any]]) -> dict[str, np.ndarray]:
    """``{label: (operator, template)}`` to ``{"a.<label>": dense matrix}`` through ``FlattenedOperator.as_matrix``."""
    return {f"a.{label}": np.asarray(FlattenedOperator(op, template).as_matrix()) for label, (op, template) in operators.items()}


# ─────────────────────────── the capture ───────────────────────────


def sn_fixture_name(kind: str, ng: int) -> str:
    return f"sn_{kind}_ng{ng}"


def diffusion_fixture_name(ng: int) -> str:
    return f"diff_slab_ng{ng}"


def capture_problems() -> tuple[dict[Key, np.ndarray], list[tuple[str, str]]]:
    """The arrays of quantity groups a, b and d of the 15 fixtures and the recorded constructor call of each."""
    arrays: dict[Key, np.ndarray] = {}
    builders: list[tuple[str, str]] = []
    for kind in fx.SN_KINDS:
        for ng in fx.GROUPS:
            name = sn_fixture_name(kind, ng)
            problem, text = fx.sn_problem(kind, ng)
            system = problem.system
            produced = _dense({
                "A": (system.loss, system.space.zeros()),
                "F": (system.production, system.space.zeros()),
                # B lives on the full-field space, not on the coupled space of the sphere and the cylinder.
                "B": (system.factors.boundary, problem.full_field_space.zeros()),
            })
            produced.update(_space_arrays(problem, SN_ROWS))
            arrays.update({(name, key): value for key, value in produced.items()})
            builders.append((name, text + " | densified by FlattenedOperator(problem.system.loss | production, system.space.zeros()) "
                                          "and FlattenedOperator(problem.system.factors.boundary, problem.full_field_space.zeros())"))
    for ng in fx.GROUPS:
        name = diffusion_fixture_name(ng)
        mesh, solver, text = fx.diffusion_problem(ng)
        produced = _dense({
            "A": (solver.loss, solver.template),
            "F": (solver.fission, solver.template),
            "B": (solver.boundary, solver.template),
        })
        produced.update(_space_arrays(mesh, DIFFUSION_ROWS))
        arrays.update({(name, key): value for key, value in produced.items()})
        builders.append((name, text + " | densified by FlattenedOperator(solver.loss | fission | boundary, solver.template)"))
    return arrays, builders


def capture_frames(families: Families | None = None) -> tuple[dict[Key, np.ndarray], list[tuple[str, str]]]:
    """Quantity group c: for each frame family, the input Gram and the two faces of ``DenseMetric.inverse_of(gram)``.

    ``metric`` is ``.matrix`` (the pseudo-inverse of the Gram) and ``gram`` is ``.inverse_matrix`` (the symmetrised Gram).
    """
    arrays: dict[Key, np.ndarray] = {}
    builders: list[tuple[str, str]] = []
    for family, text, build in (fx.frame_families() if families is None else families):
        gram = np.asarray(build(), float)
        dense = DenseMetric.inverse_of(gram)
        arrays[("frames", f"c.{family}.input_gram")] = gram
        arrays[("frames", f"c.{family}.metric")] = np.asarray(dense.matrix)
        arrays[("frames", f"c.{family}.gram")] = np.asarray(dense.inverse_matrix)
        builders.append((f"frames.{family}", f"{text}; metric = DenseMetric.inverse_of(gram)"))
    return arrays, builders


def capture_with_builders(families: Families | None = None) -> tuple[dict[Key, np.ndarray], list[tuple[str, str]]]:
    """Rebuild every fixture and every family from scratch: ``({(fixture, name): array}, [(id, constructor call)])``."""
    arrays, builders = capture_problems()
    frame_arrays, frame_builders = capture_frames(families)
    arrays.update(frame_arrays)
    return arrays, builders + frame_builders


def capture(families: Families | None = None) -> dict[Key, np.ndarray]:
    """The 336 arrays on the live tree, keyed by ``(fixture, quantity name)``."""
    return capture_with_builders(families)[0]


# ─────────────────────────── compare and the stored form ───────────────────────────


def sha1(array: np.ndarray) -> str:
    return hashlib.sha1(np.ascontiguousarray(array).tobytes()).hexdigest()


def differing_keys(a: Mapping[Key, np.ndarray], b: Mapping[Key, np.ndarray]) -> list[Key]:
    """The sorted keys at which ``a`` and ``b`` differ: present in one only, or different in shape, dtype or ``.tobytes()``.

    Never asserts. Equality is by bytes, never by ``allclose``.
    """
    out = []
    for key in set(a) | set(b):
        if key not in a or key not in b:
            out.append(key)
            continue
        x, y = np.asarray(a[key]), np.asarray(b[key])
        if x.shape != y.shape or x.dtype != y.dtype or x.tobytes() != y.tobytes():
            out.append(key)
    return sorted(out)


def host_stamp() -> tuple[str, ...]:
    """``(system, platform, machine, numpy, scipy, blas)``: what decides whether stored bytes can be compared on this host."""
    blas = np.show_config(mode="dicts")["Build Dependencies"]["blas"]["name"]
    return (platform.system(), platform.platform(), platform.machine(), np.__version__, scipy.__version__, str(blas))


def manifest_lines(arrays: Mapping[Key, np.ndarray]) -> list[str]:
    """One tab-separated line per array: fixture, name, shape, dtype, sha1 of the bytes; sorted by key."""
    return [f"{f}\t{n}\t{tuple(arrays[(f, n)].shape)}\t{arrays[(f, n)].dtype}\t{sha1(arrays[(f, n)])}" for f, n in sorted(arrays)]


def read_manifest(path: Path) -> tuple[tuple[str, ...], list[str]]:
    """The stamp (the first line, split on tabs) and the array lines of ``MANIFEST.txt``."""
    lines = Path(path).read_text().splitlines()
    if not lines:
        raise ValueError(f"{path} is empty: it must start with the host stamp line")
    return tuple(lines[0].split("\t")), lines[1:]


def load_snapshots(directory: Path) -> dict[Key, np.ndarray]:
    """Every array of every ``*.npz`` file of ``directory``, keyed by ``(file stem, array name)``."""
    arrays: dict[Key, np.ndarray] = {}
    for path in sorted(Path(directory).glob("*.npz")):
        with np.load(path) as stored:
            for name in stored.files:
                arrays[(path.stem, name)] = stored[name]
    return arrays


def directory_bytes(directory: Path) -> int:
    """The total size in bytes of every file under ``directory``, subdirectories included."""
    return sum(p.stat().st_size for p in Path(directory).rglob("*") if p.is_file())


def manifest_keys(lines: list[str]) -> set[Key]:
    """The ``(fixture, name)`` of each array line of ``MANIFEST.txt``."""
    return {(f, n) for f, n, *_ in (line.split("\t") for line in lines)}


def read_builders(path: Path) -> list[tuple[str, str]]:
    """The ``(id, constructor call)`` lines of ``BUILDERS.txt``."""
    return [tuple(line.split("\t", 1)) for line in Path(path).read_text().splitlines()]  # type: ignore[misc]


def difference_report(a: Mapping[Key, np.ndarray], b: Mapping[Key, np.ndarray], limit: int = 10) -> list[tuple[Key, str, float]]:
    """``(key, reason, max|a-b| / max|b|)`` for the first ``limit`` keys of ``differing_keys(a, b)``. Never asserts.

    The reason is ``missing from a``, ``missing from b``, ``shape``, ``dtype`` or ``bytes``; the relative difference is ``nan``
    unless the two arrays have one shape and a finite common type.
    """
    rows = []
    for key in differing_keys(a, b)[:limit]:
        if key not in a or key not in b:
            rows.append((key, f"missing from {'a' if key not in a else 'b'}", float("nan")))
            continue
        x, y = np.asarray(a[key]), np.asarray(b[key])
        if x.shape != y.shape:
            rows.append((key, f"shape {x.shape} against {y.shape}", float("nan")))
        elif x.dtype != y.dtype:
            rows.append((key, f"dtype {x.dtype} against {y.dtype}", float("nan")))
        else:
            scale = float(np.max(np.abs(y))) if y.size else 0.0
            rel = float(np.max(np.abs(x.astype(float) - y.astype(float)))) / scale if scale else float("nan")
            rows.append((key, "bytes", rel))
    return rows


def write_snapshots(directory: Path) -> dict[Key, np.ndarray]:
    """Capture on this host and write the 16 ``.npz`` files, ``MANIFEST.txt`` and ``BUILDERS.txt`` into ``directory``."""
    directory = Path(directory)
    os.makedirs(directory, exist_ok=True)
    arrays, builders = capture_with_builders()
    for stem in sorted({f for f, _ in arrays}):
        np.savez_compressed(directory / f"{stem}.npz", **{n: arrays[(f, n)] for f, n in sorted(arrays) if f == stem})
    (directory / "MANIFEST.txt").write_text("\n".join(["\t".join(host_stamp()), *manifest_lines(arrays)]) + "\n")
    (directory / "BUILDERS.txt").write_text("".join(f"{a}\t{b}\n" for a, b in builders))
    return arrays


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: python -m tests.gates.cross_method.w1._snapshot OUTDIR")
    written = write_snapshots(Path(sys.argv[1]))
    print(f"wrote {len(written)} arrays to {sys.argv[1]}; stamp {host_stamp()}")
