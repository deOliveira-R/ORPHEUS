"""Phase 3 import-inventory walker for ORPHEUS.

Walks every .py file under orpheus/, parses imports via ast, and emits a JSON
record of (importer_pkg, importee_pkg, file_relpath, line, symbol, is_type_checking).

Only edges whose importee resolves to a DIFFERENT orpheus top-level package
than the importer are recorded (intra-package edges are noise for the layer-
contract analysis).
"""
from __future__ import annotations

import ast
import json
import os
import sys
from pathlib import Path
from typing import Iterable

ORPHEUS_ROOT = Path(
    "/Users/rodrigo/git/nuclear/ORPHEUS/.claude/worktrees/moment-space-and-layering/orpheus"
)


def file_to_package(path: Path) -> str:
    """Return the top-level orpheus subpackage that owns this file.

    `orpheus/sn/solver.py` -> 'sn'
    `orpheus/__init__.py`  -> '__root__'
    `orpheus/plotting.py`  -> 'plotting'  (single-file top-level module)
    """
    rel = path.relative_to(ORPHEUS_ROOT)
    parts = rel.parts
    if len(parts) == 1:
        return "__root__" if parts[0] == "__init__.py" else parts[0].removesuffix(".py")
    return parts[0]


def importee_to_package(module: str | None) -> str | None:
    """Map a dotted module string to a top-level orpheus subpackage."""
    if module is None:
        return None
    if module == "orpheus":
        return "__root__"
    if not module.startswith("orpheus."):
        return None
    tail = module[len("orpheus.") :]
    return tail.split(".", 1)[0]


def resolve_relative_module(
    file_path: Path, level: int, module: str | None
) -> str | None:
    """Convert `from . import x`, `from ..foo import y` to an absolute dotted module."""
    rel = file_path.relative_to(ORPHEUS_ROOT)
    parts = rel.parts
    if parts[-1] == "__init__.py":
        pkg_parts = list(parts[:-1])
    else:
        pkg_parts = list(parts[:-1])
    pkg_parts = ["orpheus", *pkg_parts]
    if level > len(pkg_parts):
        return None
    base = pkg_parts[: len(pkg_parts) - (level - 1)]
    if module:
        absolute = ".".join([*base, module])
    else:
        absolute = ".".join(base)
    return absolute


class ImportVisitor(ast.NodeVisitor):
    def __init__(self, file_path: Path, importer_pkg: str):
        self.file_path = file_path
        self.importer_pkg = importer_pkg
        self.edges: list[dict] = []
        self.type_checking_stack: list[bool] = [False]

    @property
    def in_type_checking(self) -> bool:
        return any(self.type_checking_stack)

    def visit_If(self, node: ast.If) -> None:
        if self._is_type_checking_test(node.test):
            self.type_checking_stack.append(True)
            for stmt in node.body:
                self.visit(stmt)
            self.type_checking_stack.pop()
            for stmt in node.orelse:
                self.visit(stmt)
        else:
            self.generic_visit(node)

    @staticmethod
    def _is_type_checking_test(test: ast.expr) -> bool:
        if isinstance(test, ast.Name) and test.id == "TYPE_CHECKING":
            return True
        if (
            isinstance(test, ast.Attribute)
            and test.attr == "TYPE_CHECKING"
            and isinstance(test.value, ast.Name)
            and test.value.id in {"typing", "t"}
        ):
            return True
        return False

    def visit_Import(self, node: ast.Import) -> None:
        for alias in node.names:
            importee_pkg = importee_to_package(alias.name)
            if importee_pkg is None or importee_pkg == self.importer_pkg:
                continue
            self.edges.append(
                {
                    "importer_pkg": self.importer_pkg,
                    "importee_pkg": importee_pkg,
                    "importee_module": alias.name,
                    "symbol": alias.asname or alias.name,
                    "kind": "import",
                    "line": node.lineno,
                    "type_checking": self.in_type_checking,
                    "file": str(self.file_path.relative_to(ORPHEUS_ROOT.parent)),
                }
            )

    def visit_ImportFrom(self, node: ast.ImportFrom) -> None:
        if node.level and node.level > 0:
            module = resolve_relative_module(self.file_path, node.level, node.module)
        else:
            module = node.module
        importee_pkg = importee_to_package(module) if module else None
        if importee_pkg is None or importee_pkg == self.importer_pkg:
            return
        for alias in node.names:
            self.edges.append(
                {
                    "importer_pkg": self.importer_pkg,
                    "importee_pkg": importee_pkg,
                    "importee_module": module,
                    "symbol": alias.name,
                    "kind": "from",
                    "line": node.lineno,
                    "type_checking": self.in_type_checking,
                    "file": str(self.file_path.relative_to(ORPHEUS_ROOT.parent)),
                }
            )


def walk(root: Path) -> Iterable[dict]:
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if d != "__pycache__"]
        for fn in filenames:
            if not fn.endswith(".py"):
                continue
            path = Path(dirpath) / fn
            try:
                src = path.read_text(encoding="utf-8")
                tree = ast.parse(src, filename=str(path))
            except (SyntaxError, UnicodeDecodeError) as exc:
                print(f"# parse-error {path}: {exc}", file=sys.stderr)
                continue
            importer_pkg = file_to_package(path)
            visitor = ImportVisitor(path, importer_pkg)
            visitor.visit(tree)
            yield from visitor.edges


def main() -> None:
    edges = list(walk(ORPHEUS_ROOT))
    edges.sort(
        key=lambda e: (
            e["importer_pkg"],
            e["importee_pkg"],
            e["file"],
            e["line"],
        )
    )
    out_path = Path(__file__).with_name("phase3_import_edges.json")
    out_path.write_text(json.dumps(edges, indent=2))
    print(f"# wrote {len(edges)} edges -> {out_path}")

    summary: dict[str, dict[str, int]] = {}
    for e in edges:
        summary.setdefault(e["importer_pkg"], {}).setdefault(e["importee_pkg"], 0)
        summary[e["importer_pkg"]][e["importee_pkg"]] += 1
    print("\n# Per-package summary (importer -> {importee: count}):")
    for importer in sorted(summary):
        print(f"  {importer}:")
        for importee in sorted(summary[importer]):
            tc_count = sum(
                1
                for e in edges
                if e["importer_pkg"] == importer
                and e["importee_pkg"] == importee
                and e["type_checking"]
            )
            n = summary[importer][importee]
            tc_str = f" ({tc_count} TYPE_CHECKING)" if tc_count else ""
            print(f"    -> {importee}: {n}{tc_str}")


if __name__ == "__main__":
    main()
