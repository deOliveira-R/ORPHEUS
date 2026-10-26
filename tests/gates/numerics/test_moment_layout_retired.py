r"""``orpheus.numerics.moment_layout`` is retired, and its hand spellings with it (S4, #423).

Every row is a retirement witness (``retirement-audit`` A.2, D.18, E.19): it is
RED on the pre-carve tree, with the member it finds named in the failure, and
GREEN only when S3 and S4 have landed. Each AST census carries a positive
control on a synthetic source holding one member of the shape it hunts (X1:
a census that finds nothing in a clean tree and a census that cannot find
anything print the same zero), and reports its input count.

``[M]`` 2026-10-25 first reds at HEAD ``57f2d777`` (+ uncommitted S1): the
module exists; 14 production + 10 test ``ImportFrom`` statements name it; the
float root at ``fields/_bases.py:533``; the Kronecker decode at
``_ubld.py:144``; the θ literal ``1.0 / 3.0`` at ``linear_discontinuous.py:331``.
"""

from __future__ import annotations

import ast
import importlib.util
import pathlib

import pytest

import orpheus

pytestmark = [pytest.mark.foundation]

ROOT = pathlib.Path(orpheus.__file__).resolve().parents[1]
RETIRED = "orpheus.numerics.moment_layout"
#: The tail rule's home (ruling 4: numerics space composition, beside
#: ``FunctionSpace.of_axes``; named ``factor_tail`` by the user, 2026-10-25).
TAIL_RULE_HOME = ("orpheus.numerics.space", "factor_tail")


def _py_files(*roots: str) -> list[pathlib.Path]:
    out = [p for r in roots for p in (ROOT / r).rglob("*.py") if "__pycache__" not in p.parts]
    assert out, f"the census read no files under {roots}"
    return out


def _imports_of(tree: ast.AST, module: str) -> list[int]:
    lines = []
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module is not None:
            if node.module == module or (node.module == module.rsplit(".", 1)[0] and any(a.name == module.rsplit(".", 1)[1] for a in node.names)):
                lines.append(node.lineno)
        elif isinstance(node, ast.Import) and any(a.name == module for a in node.names):
            lines.append(node.lineno)
    return sorted(lines)


def test_the_module_is_gone() -> None:
    assert importlib.util.find_spec(RETIRED) is None, f"{RETIRED} still exists"


def test_no_file_imports_the_retired_module() -> None:
    r"""Top-level, late (in-function) and ``TYPE_CHECKING`` imports all count; ``orpheus/``, ``tests/``, ``tools/``."""
    control = ast.parse("def f():\n    from orpheus.numerics.moment_layout import AVERAGE_MOMENT\nfrom orpheus.numerics import moment_layout\n")
    assert _imports_of(control, RETIRED) == [2, 3], "the import census lost a known member (positive control)"
    files = _py_files("orpheus", "tests", "tools")
    hits = [f"{p.relative_to(ROOT)}:{line}" for p in files for line in _imports_of(ast.parse(p.read_text(encoding="utf-8")), RETIRED)]
    assert not hits, f"{len(hits)} import(s) of {RETIRED} over {len(files)} files: {hits}"


def _float_roots(tree: ast.AST) -> list[int]:
    r"""``x ** (1.0 / n)``: a per-axis width recovered by a float root."""
    return [n.lineno for n in ast.walk(tree)
            if isinstance(n, ast.BinOp) and isinstance(n.op, ast.Pow) and isinstance(n.right, ast.BinOp)
            and isinstance(n.right.op, ast.Div) and isinstance(n.right.left, ast.Constant) and n.right.left.value == 1.0]


def _kronecker_decodes(tree: ast.AST) -> list[int]:
    r"""``(flat // p ** e) % p``: a Kronecker index decoded by hand."""
    return [n.lineno for n in ast.walk(tree)
            if isinstance(n, ast.BinOp) and isinstance(n.op, ast.Mod) and isinstance(n.left, ast.BinOp)
            and isinstance(n.left.op, ast.FloorDiv) and isinstance(n.left.right, ast.BinOp) and isinstance(n.left.right.op, ast.Pow)]


def _theta_literals(tree: ast.AST) -> list[int]:
    r"""``1.0 / 3.0`` (or ``1 / 3``): the LD Gram typed as a constant."""
    return [n.lineno for n in ast.walk(tree)
            if isinstance(n, ast.BinOp) and isinstance(n.op, ast.Div) and isinstance(n.left, ast.Constant)
            and isinstance(n.right, ast.Constant) and n.left.value == 1 and n.right.value == 3]


_CENSUSES = {
    "float root": (_float_roots, "per_axis = round(tail[0] ** (1.0 / ndim))\n", ("orpheus/transport/fields/_bases.py",)),
    "Kronecker decode": (_kronecker_decodes, "o = (flat // (p ** (d - 1 - a))) % p\n", ("orpheus/transport/spatial/_ubld.py", "orpheus/sn")),
    "theta literal": (_theta_literals, "theta = 1.0 / 3.0\n", ("orpheus/transport/spatial/linear_discontinuous.py", "orpheus/transport/spatial/_ubld.py", "orpheus/transport/spatial/scheme.py")),
}


@pytest.mark.parametrize("name", list(_CENSUSES))
def test_no_hand_spelling_of_the_cell_basis_survives(name: str) -> None:
    r"""The four hand spellings of the plan's finding re-point to the basis; none survives in the files that held one.

    The literal ``2`` of the LD static methods has no reliable AST shape (a
    ``2`` is also an honest integer), so it is a grep OBLIGATION of the S4
    dependency audit, not a row here: say so in the audit.
    """
    census, witness, paths = _CENSUSES[name]
    assert census(ast.parse(witness)) == [1], f"the {name} census lost its known member (positive control)"
    files = [p for path in paths for p in ([ROOT / path] if path.endswith(".py") else _py_files(path))]
    assert all(p.exists() for p in files) and files
    hits = [f"{p.relative_to(ROOT)}:{line}" for p in files for line in census(ast.parse(p.read_text(encoding="utf-8")))]
    assert not hits, f"{len(hits)} {name}(s) survive over {len(files)} files: {hits}"


def test_the_tail_rule_lives_in_space_composition() -> None:
    r"""Ruling 4: the "append a trailing axis iff its width > 1" rule is numerics space composition.

    The value legs are the rule itself: width 1 appends nothing (DD/Step buffers byte-identical, #240
    D5b), width > 1 appends one axis.
    """
    module, name = TAIL_RULE_HOME
    rule = getattr(importlib.import_module(module), name)
    assert rule.__module__ == module
    assert rule(1) == () and rule(2) == (2,) and rule(4) == (4,)

