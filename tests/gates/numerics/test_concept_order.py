"""The ``numerics`` layer imports in one acyclic concept order.

Each module of :mod:`orpheus.numerics` has a RANK: its place in the
concept order of the numerics-untangling plan (candidate E′,
``.claude/plans/numerics_untangling.md``). A concept is defined from
the concepts before it and never reaches one after it, so an import is
legal only from a lower rank into a higher one. All three kinds of
import count: top-level, lazy (inside a function body) and
``TYPE_CHECKING``. A lazy or typing import of a later concept is the
same dependence as a top-level one, only hidden from the interpreter.

``PENDING`` lists every edge that violates the order today, each with
the plan package that removes it. The list only shrinks: a package
lands with its entries deleted, and
:func:`test_pending_entries_still_violate` refuses an entry that no
longer violates, so a stale exemption cannot outlive its edge.
Adding an entry is a review blocker. A new module must be given a rank
(:func:`test_every_module_has_a_rank`).

The grain is the module directly under ``orpheus/numerics/``: the
packages ``basis``, ``quadrature`` and ``spaces`` are one node each,
and imports inside one node are not edges.
"""

from __future__ import annotations

import ast
import pathlib

import pytest

NUMERICS = pathlib.Path(__file__).resolve().parents[3] / "orpheus" / "numerics"

# The concept order, upstream first. Modules in one tier have no edge
# between them in the target order.
TIERS: tuple[tuple[str, ...], ...] = (
    # coordinate-space membership, the shape leg of the one parse: below every
    # concept, so that a measure, a basis and a metric block, which own their
    # domain shapes but cannot name a space, parse with it (U1b, 2026-10-30)
    ("coordinate_space",),
    # primitives, and the orthogonal-polynomial recurrence (upstream of both measures)
    # (``symmetry``, ``roots_of_unity`` and ``manifold``'s geometric half left for
    # geometry at layer DAG WP-B, 2026-10-22; ``index_set`` holds the finite index
    # sets, members of geometry's ``Manifold`` sum, and imports no numerics module;
    # ``tangential`` holds the one grazing tolerance both ``spaces`` and ``quadrature``
    # read, a leaf so that quadrature need not import spaces, layer DAG WF2.2)
    ("evidence", "face_layout", "gauge", "index_set", "instance_cache", "metric",
     "recurrence", "registry", "tangential", "units", "vector"),
    # the coefficient layouts are content values a basis carries (ruling 81); ``_import_time``
    # is the traced memo's rule for which import-time statements a run observed, and imports
    # nothing from ``orpheus`` (cone store B2, #689)
    ("_import_time", "_traced_memo_boot", "coefficient_layout", "enclosure", "mesh_free_function"),
    ("traced_memo",),
    ("basis", "exactness", "observable", "question"),
    ("measure",),
    # the frame's data = join(basis, measure); the Gauss rule. Invariance is no module:
    # S5 (ruling 42) made it the measure's composition of pushforward with its
    # G-space support's action, so it ranks with the measure.
    ("frame", "generating_measure"),
    ("quadrature",),
    ("axis",),
    ("space",),
    ("field",),
    ("functional",),
    # the algebra, with the adjoint, the dual and the Riesz operators (rulings 25, 26)
    ("operator",),
    ("assembled_operator", "axis_marginal", "elementary_operators", "flat_operator",
     "projection"),
    # the LU inverse densifies through the assembly functor (ruling 60)
    ("matrix_inverse_operator",),
    ("coupled_system", "frame_faces"),
    # the kernel gauge composes frame faces over spanning summands (U3, 2026-10-31)
    ("kernel_projection_gauge", "spaces"),
    # the Fredholm verdict decides on a predicted kernel's projector and annihilators (U0, 2026-10-31);
    # the pencil carries the predicted kernel and the posings derive their verdict from it (the user's
    # ruling, 2026-10-31), so both rank above it, and the exit report records it, so the outcomes rank above them.
    # The near-kernel witness measures what a predicted kernel claims and imports none of the verdict (U3b,
    # 2026-10-31): it reads the algebra, the assembly and the kernel gauge's metric door, so it ranks beside fredholm.
    ("fredholm", "near_kernel"),
    ("pencil",),
    ("posing",),
    ("outcome",),
    ("convergence",),
    ("eigenvalue",),
    ("iteration",),
    ("green_operator",),
)
RANK: dict[str, int] = {m: i for i, tier in enumerate(TIERS) for m in tier}

# The memo's boot split is deliberate: the boot module imports the memo
# lazily so the memo can boot itself (plan finding F7).
EXEMPT: frozenset[tuple[str, str]] = frozenset({("traced_memo", "_traced_memo_boot")})

# (imported module, importing module) -> the plan package that removes the edge.
PENDING: dict[tuple[str, str], str] = {}


def _node(path: pathlib.Path, root: pathlib.Path) -> str:
    return path.relative_to(root).parts[0].removesuffix(".py")


def _resolve(node: ast.ImportFrom, path: pathlib.Path, root: pathlib.Path) -> list[str]:
    """Dotted targets of one ``from ... import`` (each name may be a submodule)."""
    # the package a relative import is resolved against: the module's parent
    package = path.relative_to(root.parent.parent).with_suffix("").parts[:-1]
    if node.level:
        base = list(package[: len(package) - (node.level - 1)])
        module = ".".join(base + ([node.module] if node.module else []))
    else:
        module = node.module or ""
    return [f"{module}.{a.name}" for a in node.names] + [module]


def _is_type_checking_guard(test: ast.expr) -> bool:
    """``if TYPE_CHECKING:`` or ``if typing.TYPE_CHECKING:``, exactly (``if not TYPE_CHECKING:`` guards the runtime
    branch). The same predicate as ``tests/gates/test_layer_imports.py``'s; one home for both is #695."""
    return (isinstance(test, ast.Name) and test.id == "TYPE_CHECKING") or (
        isinstance(test, ast.Attribute) and test.attr == "TYPE_CHECKING")


#: The kinds of import edge: ``runtime`` a top-level import, ``lazy`` an import inside a def or a lambda,
#: ``typing`` an import under ``if TYPE_CHECKING:`` (erased at runtime wherever it sits, so it wins over lazy).
KINDS: tuple[str, ...] = ("runtime", "lazy", "typing")


def import_kinds(root: pathlib.Path = NUMERICS) -> dict[tuple[str, str], frozenset[str]]:
    """Every (imported, importer) pair between distinct top-level nodes of ``root``, with the kinds of import
    that make it (:data:`KINDS`). The one walk of this order's import graph: the gate reads its edges, and the
    cone store's placement (``tools/cone_store/placement.py``) reads its kinds to tell a test that composes two
    concepts from one that reaches a concept through an import the interpreter never runs."""
    nodes = {p.name.removesuffix(".py") for p in root.iterdir()
             if (p.suffix == ".py" and p.name != "__init__.py") or (p.is_dir() and (p / "__init__.py").exists())}
    prefix = f"{root.parent.name}.{root.name}."
    edges: dict[tuple[str, str], set[str]] = {}

    def visit(node: ast.AST, kind: str, path: pathlib.Path, importer: str) -> None:
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            targets = [a.name for a in node.names] if isinstance(node, ast.Import) else _resolve(node, path, root)
            for target in targets:
                if target.startswith(prefix):
                    head = target.removeprefix(prefix).split(".")[0]
                    if head in nodes and head != importer:
                        edges.setdefault((head, importer), set()).add(kind)
        if kind == "runtime" and isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            kind = "lazy"
        for child in ast.iter_child_nodes(node):
            if isinstance(child, ast.If) and _is_type_checking_guard(child.test):
                for statement in child.body:
                    visit(statement, "typing", path, importer)
                for statement in child.orelse:
                    visit(statement, kind, path, importer)
            else:
                visit(child, kind, path, importer)

    for path in root.rglob("*.py"):
        if path.parent == root and path.name == "__init__.py":
            continue
        visit(ast.parse(path.read_text(encoding="utf-8")), "runtime", path, _node(path, root))
    return {edge: frozenset(kinds) for edge, kinds in edges.items()}


def import_edges(root: pathlib.Path = NUMERICS) -> set[tuple[str, str]]:
    """Every (imported, importer) pair between distinct top-level nodes of ``root``, of any kind."""
    return set(import_kinds(root))


def violations(edges: set[tuple[str, str]]) -> set[tuple[str, str]]:
    return {(a, b) for a, b in edges if RANK[a] >= RANK[b]} - EXEMPT


def _modules() -> set[str]:
    return {p.name.removesuffix(".py") for p in NUMERICS.iterdir()
            if (p.suffix == ".py" and p.name != "__init__.py") or (p.is_dir() and (p / "__init__.py").exists())}


def test_every_module_has_a_rank():
    unranked = sorted(_modules() - RANK.keys())
    assert not unranked, f"numerics modules with no place in the concept order: {unranked}"


def test_the_concept_order_names_only_modules_that_exist():
    r"""The other half of :func:`test_every_module_has_a_rank`: no rank names a module that is gone.

    That test refuses an unranked module; this one refuses a ranked module
    that no longer exists, which would otherwise stay silent after a
    retirement (``moment_layout`` left tier 0 at the moment-basis carve, #423).
    """
    present = _modules()
    assert "frame" in present and "basis" in present, "the module walk lost known members (positive control)"
    stale = sorted(RANK.keys() - present)
    assert not stale, f"the concept order ranks modules that do not exist: {stale}"


def test_import_edges_follow_the_concept_order():
    unexplained = sorted(violations(import_edges()) - PENDING.keys())
    assert not unexplained, (
        "imports from a later concept into an earlier one (imported, importer): "
        f"{unexplained}. Re-home the dependence; never add it to PENDING."
    )


def test_pending_entries_still_violate():
    stale = sorted(PENDING.keys() - violations(import_edges()))
    assert not stale, f"PENDING entries whose edge no longer violates the order; delete them: {stale}"


def test_the_order_detects_a_backward_import(tmp_path):
    """Positive control: a planted lazy import of a later concept is found."""
    root = tmp_path / "orpheus" / "numerics"
    root.mkdir(parents=True)
    (root / "__init__.py").write_text("")
    (root / "measure.py").write_text(
        "def f():\n    from orpheus.numerics.axis import Axis\n    return Axis\n")
    (root / "axis.py").write_text(
        "from typing import TYPE_CHECKING\nif TYPE_CHECKING:\n"
        "    from .measure import DiscreteMeasure\n")
    edges = import_edges(root)
    assert edges == {("axis", "measure"), ("measure", "axis")}
    assert violations(edges) == {("axis", "measure")}


def test_each_import_edge_carries_its_kind(tmp_path):
    """Positive control of :func:`import_kinds`: one planted import of each kind is found with that kind, a
    ``TYPE_CHECKING`` import inside a def is ``typing`` (erased at runtime, wherever it sits), and an edge made by
    two kinds carries both."""
    root = tmp_path / "orpheus" / "numerics"
    root.mkdir(parents=True)
    (root / "__init__.py").write_text("")
    for name in ("a", "b", "c", "d"):
        (root / f"{name}.py").write_text("")
    (root / "e.py").write_text(
        "from typing import TYPE_CHECKING\n"
        "from orpheus.numerics import a\n"
        "def f():\n    from . import b\n    if TYPE_CHECKING:\n        from orpheus.numerics.d import X\n"
        "if TYPE_CHECKING:\n    from .c import Y\n"
        "g = lambda: __import__('x')\n"
        "def h():\n    import orpheus.numerics.a\n")
    assert import_kinds(root) == {
        ("a", "e"): frozenset({"runtime", "lazy"}), ("b", "e"): frozenset({"lazy"}),
        ("c", "e"): frozenset({"typing"}), ("d", "e"): frozenset({"typing"})}


def test_import_kinds_reads_the_branch_the_interpreter_runs(tmp_path):
    """Witnesses of :func:`import_kinds`' refinements (qa's F4 on B3, 2026-11-03), one planted import each: the
    ``else`` branch of ``if TYPE_CHECKING:`` runs, so its import is ``runtime``; ``if not TYPE_CHECKING:`` is not the
    guard (its body runs), so its import is ``runtime``; an import inside an ``async def`` is ``lazy``.

    First red ``[M]`` 2026-11-03, each a textual mutant of :func:`import_kinds` rebound by a ``-p`` plugin: the
    ``else`` branch visited as ``typing`` reddens on ``a``; the guard tested as a substring (``"TYPE_CHECKING" in
    ast.unparse(test)``) on ``b``; ``ast.AsyncFunctionDef`` dropped from the def tuple on ``c``."""
    root = tmp_path / "orpheus" / "numerics"
    root.mkdir(parents=True)
    (root / "__init__.py").write_text("")
    for name in ("a", "b", "c"):
        (root / f"{name}.py").write_text("")
    (root / "e.py").write_text(
        "from typing import TYPE_CHECKING\n"
        "if TYPE_CHECKING:\n    pass\nelse:\n    from orpheus.numerics import a\n"
        "if not TYPE_CHECKING:\n    from . import b\n"
        "async def f():\n    from .c import Z\n")
    assert import_kinds(root) == {
        ("a", "e"): frozenset({"runtime"}), ("b", "e"): frozenset({"runtime"}), ("c", "e"): frozenset({"lazy"})}
