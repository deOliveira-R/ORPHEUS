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
    # primitives, and the orthogonal-polynomial recurrence (upstream of both measures)
    ("face_layout", "gauge", "metric", "moment_layout", "recurrence", "registry",
     "roots_of_unity", "scalars", "units", "vector"),
    ("content", "symmetry"),
    ("_traced_memo_boot", "enclosure", "manifold", "mesh_free_function"),
    ("traced_memo",),
    ("basis", "exactness", "observable", "question"),
    ("measure",),
    # the frame's data = join(basis, measure); the Gauss rule; invariance = join(measure, G-space)
    ("frame", "generating_measure", "invariance"),
    ("quadrature",),
    ("axis",),
    ("space",),
    ("field",),
    ("functional",),
    # the algebra, with the adjoint, the dual and the Riesz operators (rulings 25, 26)
    ("operator",),
    ("assembled_operator", "axis_marginal", "elementary_operators", "flat_operator",
     "matrix_inverse_operator", "pencil", "projection"),
    ("coupled_system", "frame_faces", "posing"),
    ("outcome", "spaces"),
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
PENDING: dict[tuple[str, str], str] = {
    # S7: the coefficient space of a basis and the sampled space of a measure (ruling 35)
    ("space", "basis"): "S7",
    ("space", "measure"): "S7",
    ("spaces", "basis"): "S7",
    # M6: the measure stops reaching the generating measure (ruling 12)
    ("generating_measure", "measure"): "M6",
    # S1: the marginalisation lift is built on the operator side (ruling 10)
    ("axis", "frame"): "S1",
    ("operator", "frame"): "S1",
    ("operator", "space"): "S1",
    ("space", "frame"): "S1",
    # S3: assembly is a homomorphism downstream of the algebra (ruling 24)
    ("assembled_operator", "operator"): "S3",
    # S4: a sum's inverse is InverseOperator(sum), realised by the Strategy (ruling 20)
    ("green_operator", "operator"): "S4",
    # S5: invariance is a composition the measure does not own (ruling 13)
    ("invariance", "measure"): "S5",
}


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


def import_edges(root: pathlib.Path = NUMERICS) -> set[tuple[str, str]]:
    """Every (imported, importer) pair between distinct top-level nodes of ``root``."""
    nodes = {p.name.removesuffix(".py") for p in root.iterdir()
             if (p.suffix == ".py" and p.name != "__init__.py") or (p.is_dir() and (p / "__init__.py").exists())}
    prefix = f"{root.parent.name}.{root.name}."
    edges: set[tuple[str, str]] = set()
    for path in root.rglob("*.py"):
        if path.parent == root and path.name == "__init__.py":
            continue
        importer = _node(path, root)
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if isinstance(node, ast.Import):
                targets = [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom):
                targets = _resolve(node, path, root)
            else:
                continue
            for target in targets:
                if target.startswith(prefix):
                    head = target.removeprefix(prefix).split(".")[0]
                    if head in nodes and head != importer:
                        edges.add((head, importer))
    return edges


def violations(edges: set[tuple[str, str]]) -> set[tuple[str, str]]:
    return {(a, b) for a, b in edges if RANK[a] >= RANK[b]} - EXEMPT


def _modules() -> set[str]:
    return {p.name.removesuffix(".py") for p in NUMERICS.iterdir()
            if (p.suffix == ".py" and p.name != "__init__.py") or (p.is_dir() and (p / "__init__.py").exists())}


def test_every_module_has_a_rank():
    unranked = sorted(_modules() - RANK.keys())
    assert not unranked, f"numerics modules with no place in the concept order: {unranked}"


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
