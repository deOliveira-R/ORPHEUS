r"""Shared machinery of the gates on ``orpheus.transport.boundary_conditions`` (wave 2, slice S2a).

Not a test module (no ``test_`` prefix). Spec: ``.claude/plans/wave2_boundary_carve.md`` sections 2.4, 3, 6.2 (S2a).

What this module holds, and why each piece exists:

* :data:`bc`, the package read AT RUN TIME. Nothing here imports the package at collection: at the slice's base tip
  the package does not exist, and a collection error would hide which gate fails for which reason (the pattern of
  ``tests/gates/geometry/test_deck.py``). Every row reads ``bc.<name>`` at first use, so each red names its own line.
* :func:`load`, the isolation shim. ``orpheus/transport/__init__.py`` imports ``fields``, whose
  ``angular_boundary_source_sink`` imports ``TransportMethod`` from ``method.py``, which imports the retired
  ``orpheus.geometry.BC`` (slice S1 dropped the export). So, until slice S2b re-spells ``method.py``,
  ``import orpheus.transport.<anything>`` raises ``ImportError`` however correct the new package is, and the
  package could not be "green in isolation" (the S2a brief) without a way past that parent. ``load`` imports the
  package normally and, ONLY when that raises because the PARENT ``orpheus.transport`` cannot import, installs a bare
  namespace parent (a module with the real directory as its ``__path__``, no ``__init__`` run) and imports the
  package under it. After S2b the first ``import`` succeeds and the fallback never runs; the shim is then dead code
  and is deleted with the S2b amendment (``NEEDS`` of the S2a report). The shim is NOT a silent skip: a package that
  does not exist still raises ``ModuleNotFoundError`` from the second import, which is the base red of every row.
* the builders every gate file shares (the cone grid, the samplers, the test-local conditions, a real slab mesh, the
  scanners of the structure gates with their planted positive controls).

Nothing here depends on the package at import time: the test-local classes that subclass a production base are made
inside cached functions.
"""

from __future__ import annotations

import ast
import dataclasses
import functools
import importlib
import sys
import types
from pathlib import Path
from typing import Any

import numpy as np

from tests.gates._content_identity_helpers import require

import orpheus as _orpheus

#: The tree being tested, read from the imported ``orpheus`` (not from this file's path, which a symlinked or
#: copied test tree would resolve elsewhere).
ROOT = Path(_orpheus.__file__).parent.parent
PACKAGE = "orpheus.transport.boundary_conditions"
PACKAGE_DIR = Path(_orpheus.__file__).parent / "transport" / "boundary_conditions"

#: The ten modules of section 2.4's table, exactly (``__init__`` re-exports and nothing else).
MODULES = (
    "__init__", "admission", "binding", "condition", "declared", "displacement", "errors", "inflow", "network", "realization",
)


def load() -> types.ModuleType:
    """The package, imported; see the module docstring for the one fallback."""
    try:
        return importlib.import_module(PACKAGE)
    except ImportError as first:
        if PACKAGE in str(first) and isinstance(first, ModuleNotFoundError):
            raise  # the package itself does not exist: the honest base red
        parent = "orpheus.transport"
        if parent in sys.modules and getattr(sys.modules[parent], "__file__", None) is not None:
            raise  # the real parent imported; the failure is the package's own
        stub = types.ModuleType(parent)
        stub.__path__ = [str(PACKAGE_DIR.parent)]  # type: ignore[attr-defined]
        sys.modules[parent] = stub
        return importlib.import_module(PACKAGE)


class _Lazy:
    """The production package, read at first attribute access (see the module docstring)."""

    def __getattr__(self, attr: str) -> Any:
        return getattr(load(), attr)


bc: Any = _Lazy()


def node(path: str, name: str) -> str:
    """A pytest node id of a file under ``tests/gates/transport/`` (for ``rests_on``)."""
    return f"tests/gates/transport/{path}::{name}"


# ═════════════════════════════════════════════════════════════════════════════
# values the gates build
# ═════════════════════════════════════════════════════════════════════════════


@dataclasses.dataclass(frozen=True)
class ProfileSampler:
    """An ENCODABLE inflow sampler: a frozen ``eq=True`` dataclass with a float field (content by the encoder's
    frozen-dataclass rule), satisfying the ``InflowSampler`` protocol. Used where a roster leg needs a profile with
    content; the contentless case is :class:`OpaqueSampler`."""

    scale: float = 1.0

    def sample(self, space: Any) -> np.ndarray:
        return self.scale * np.ones(space.shape)


class OpaqueSampler:
    """A callable-like sampler with NO content: a plain object, keyed by identity (the three tree implementors of
    ``evaluate(space)`` are of this kind, section 2.4 "the frame is the type, the content is open")."""

    def __init__(self, value: float = 1.0) -> None:
        self.value = value
        self.calls: list[Any] = []

    def sample(self, space: Any) -> np.ndarray:
        self.calls.append(space)
        return self.value * np.ones(space.shape)


@functools.cache
def local() -> types.SimpleNamespace:
    """Test-local subclasses of ``BoundaryCondition`` (made lazily: they subclass a production base).

    * ``Probe(rho)``: a method-only condition with a chosen return fraction, face-local, no source;
    * ``ProbeChild``: a subclass of ``Probe`` (the exact-type lookup row);
    * ``SourcedProbe``: a method-only condition whose ``source`` is nonzero (the accessor-only reading of sources);
    * ``NegativeProbe``: ``return_fraction == -1.0`` (a ``ZeroFluxCondition`` stand-in: diffusion's type may not be
      imported here, and SN refuses it on the state-cone axis without importing it).
    """
    base = bc.BoundaryCondition

    @dataclasses.dataclass(frozen=True, eq=False)
    class Probe(base):  # type: ignore[misc, valid-type]
        rho: float = 0.5

        @property
        def return_fraction(self) -> float:
            return self.rho

        @property
        def source(self) -> Any:
            return bc.InflowTerms.empty()

        @property
        def face_local(self) -> bool:
            return True

    @dataclasses.dataclass(frozen=True, eq=False)
    class ProbeChild(Probe):
        pass

    @dataclasses.dataclass(frozen=True, eq=False)
    class NegativeProbe(base):  # type: ignore[misc, valid-type]
        @property
        def return_fraction(self) -> float:
            return -1.0

        @property
        def source(self) -> Any:
            return bc.InflowTerms.empty()

        @property
        def face_local(self) -> bool:
            return True

    @dataclasses.dataclass(frozen=True, eq=False)
    class SourcedProbe(base):  # type: ignore[misc, valid-type]
        """A method-only condition with a NONZERO source: ``Binding.sources`` must read it through the accessor of
        the base, never through an ``isinstance`` on ``WallCondition``."""

        amount: float = 3.0

        @property
        def return_fraction(self) -> float:
            return 0.0

        @property
        def source(self) -> Any:
            return bc.InflowTerms(self.amount, ())

        @property
        def face_local(self) -> bool:
            return True

    return types.SimpleNamespace(Probe=Probe, ProbeChild=ProbeChild, NegativeProbe=NegativeProbe, SourcedProbe=SourcedProbe)


def weights(**by_chain: float) -> Any:
    """``ResponseNetwork`` of ``specular=``/``diffuse=`` weights (the raw constructor: the cone parse is the
    condition's, the network parses only its own keys and reals)."""
    from orpheus.values.content import FrozenMapping

    table = {}
    if "specular" in by_chain:
        table[bc.SpecularChain()] = by_chain["specular"]
    if "diffuse" in by_chain:
        table[bc.DiffuseChain()] = by_chain["diffuse"]
    return bc.ResponseNetwork(FrozenMapping(table))


def cone_grid() -> list[tuple[float, float]]:
    """``(s, d)`` on and inside the cone ``s, d >= 0, s + d <= 1``: the dyadic lattice of step 1/8, plus the
    non-dyadic points the section names (0.3, 0.4), (0.9, 0.1), (0.7, 0.2)."""
    lattice = [(i / 8, j / 8) for i in range(9) for j in range(9 - i)]
    return lattice + [(0.3, 0.4), (0.9, 0.1), (0.7, 0.2)]


def wall(s: float, d: float) -> Any:
    """The condition at the cone point ``(s, d)`` built through ``from_amplitudes`` (the section of ``cone_point``)."""
    return bc.WallCondition.from_amplitudes(s, d)


def slab_mesh(roles: dict[str, Any], *, edges: tuple[float, ...] = (0.0, 1.0, 2.0), materials: tuple[int, ...] = (0,),
              region_ids: tuple[int, ...] | None = None) -> Any:
    """A real 1-D Cartesian ``Mesh1D`` (a ``BoundaryCarrier`` after slice S1): ``materials`` are the region
    materials, ``region_ids`` the region of each cell (default: the cells in order, one region per material)."""
    from orpheus.geometry import CoordSystem
    from orpheus.mesh import Mesh1D

    cart = CoordSystem.CARTESIAN
    arr = np.asarray(edges, dtype=float)
    ids = np.asarray(region_ids if region_ids is not None else range(len(edges) - 1), dtype=int)
    if region_ids is None and len(materials) == 1:
        ids = np.zeros(len(edges) - 1, dtype=int)
    return Mesh1D(coord=cart, edges=arr, volumes=cart.measure(arr), region_ids=ids, region_materials=materials, face_roles=roles)


def mesh2d(roles: dict[str, Any]) -> Any:
    """A real ``Mesh2D`` (Cartesian, 2 x 3 cells, one material)."""
    from orpheus.mesh import Mesh2D

    return Mesh2D(np.array([0.0, 1.0, 2.0]), np.array([0.0, 1.0, 2.0, 3.0]), np.zeros((2, 3), dtype=int), face_roles=roles)


@dataclasses.dataclass(frozen=True)
class Carrier:
    """A test-local ``BoundaryCarrier``: the object owning the faces and their colouring symmetry (decision 51)."""

    face_roles: dict[str, Any]
    symmetry: Any


def function_space(n: int = 3) -> Any:
    """A real ``FunctionSpace`` of shape ``(n,)`` (what ``sample_terms`` reads)."""
    from orpheus.numerics.axis import EnergyAxis
    from orpheus.numerics.space import FunctionSpace

    return FunctionSpace.of_axes(EnergyAxis.synthetic(n))


def no_orpheus_frame(error: BaseException) -> bool:
    """True iff no frame of the traceback is a file of ``orpheus/`` (the error came from Python's operator
    protocol, not from a production ``raise``)."""
    import orpheus

    package_dir = Path(orpheus.__file__).resolve().parent
    tb = error.__traceback__
    while tb is not None:
        if Path(tb.tb_frame.f_code.co_filename).resolve().is_relative_to(package_dir):
            return False
        tb = tb.tb_next
    return True


# ═════════════════════════════════════════════════════════════════════════════
# the scanners of the structure gates, each with the planted source that proves it can fire (X1)
# ═════════════════════════════════════════════════════════════════════════════


def imports_of(source: str, module_name: str = "orpheus.transport.boundary_conditions.x") -> list[str]:
    """Every module a source text imports, ABSOLUTE (a relative import is resolved against ``module_name``'s
    package), at ANY nesting (function-local, ``if TYPE_CHECKING:``): an ``ast.walk`` over the whole tree.

    ``import a.b`` names ``a.b``; ``from a import b`` names ``a`` AND ``a.b`` (``b`` may be a submodule), which is
    the over-approximation a forbidden-prefix test wants.
    """
    package = module_name.rsplit(".", 1)[0] if "." in module_name else module_name
    names: list[str] = []
    for item in ast.walk(ast.parse(source)):
        if isinstance(item, ast.Import):
            names += [alias.name for alias in item.names]
        elif isinstance(item, ast.ImportFrom):
            if item.level:
                base = package.split(".")
                base = base[: len(base) - (item.level - 1)] if item.level > 1 else base
                head = ".".join(base + ([item.module] if item.module else []))
            else:
                head = item.module or ""
            names.append(head)
            names += [f"{head}.{alias.name}" for alias in item.names if alias.name != "*"]
    return names


def violations_of_imports(source: str, forbidden: tuple[str, ...], module_name: str = "orpheus.transport.boundary_conditions.x") -> list[str]:
    """The imports of ``source`` that equal or lie under a forbidden module (``orpheus.sn`` forbids
    ``orpheus.sn.problem``, never ``orpheus.snake``)."""
    return sorted({
        name for name in imports_of(source, module_name)
        if any(name == f or name.startswith(f + ".") for f in forbidden)
    })


def matches_naming(source: str, class_names: frozenset[str]) -> tuple[int, list[str]]:
    """``(population, violations)`` of one source text: the ``match`` statements one of whose cases is a class
    pattern naming one of ``class_names`` (annotation-free: the subject of ``cone_point``'s match is an unannotated
    loop variable), and those among them whose last ``case`` is not a guard-free catch-all calling ``assert_never``.

    The same predicate as ``tests/gates/geometry/test_deck.py``'s scanner, keyed on the CLASSES a case names
    instead of the subject's annotation, because the subject of the chain match is not annotated.
    """
    population, violations = 0, []
    for fn in ast.walk(ast.parse(source)):
        if not isinstance(fn, ast.FunctionDef | ast.AsyncFunctionDef):
            continue
        for item in ast.walk(fn):
            if not isinstance(item, ast.Match):
                continue
            named = {
                sub.cls.id for case in item.cases for sub in ast.walk(case.pattern)
                if isinstance(sub, ast.MatchClass) and isinstance(sub.cls, ast.Name)
            }
            if not named & class_names:
                continue
            population += 1
            last = item.cases[-1]
            catch_all = isinstance(last.pattern, ast.MatchAs) and last.pattern.pattern is None and last.guard is None
            call = last.body[-1].value if isinstance(last.body[-1], ast.Expr) else None
            called = set()
            if isinstance(call, ast.Call):
                called = {call.func.id} if isinstance(call.func, ast.Name) else {call.func.attr} if isinstance(call.func, ast.Attribute) else set()
            if not (catch_all and "assert_never" in called):
                violations.append(f"{fn.name}:{item.lineno}")
    return population, violations


def chain_classes_in(source: str) -> set[str]:
    """The classes of a source text that carry a class-level ``factors`` attribute (what makes a class a chain)."""
    found = set()
    for item in ast.walk(ast.parse(source)):
        if isinstance(item, ast.ClassDef):
            for stmt in item.body:
                target = stmt.target if isinstance(stmt, ast.AnnAssign) else stmt.targets[0] if isinstance(stmt, ast.Assign) else None
                if isinstance(target, ast.Name) and target.id == "factors":
                    found.add(item.name)
    return found


def orientation_fields(cls: type) -> list[str]:
    """The names ``axis`` / ``outward_sign`` among a class's dataclass fields and class attributes."""
    names = {f.name for f in dataclasses.fields(cls)} if dataclasses.is_dataclass(cls) else set()
    names |= {n for k in cls.__mro__ for n in vars(k) if not n.startswith("__")}
    return sorted(names & {"axis", "outward_sign"})


def walk_subclasses(root: type) -> list[type]:
    """``root`` and every (transitive) subclass, once."""
    seen: list[type] = []
    stack = [root]
    while stack:
        cls = stack.pop()
        if cls not in seen:
            seen.append(cls)
            stack.extend(cls.__subclasses__())
    return seen


__all__ = [
    "Carrier", "MODULES", "OpaqueSampler", "PACKAGE", "PACKAGE_DIR", "ProfileSampler", "ROOT", "bc",
    "chain_classes_in", "cone_grid", "function_space", "imports_of", "load", "local", "matches_naming", "mesh2d",
    "no_orpheus_frame", "node", "orientation_fields", "require", "slab_mesh", "violations_of_imports",
    "walk_subclasses", "wall", "weights",
]
