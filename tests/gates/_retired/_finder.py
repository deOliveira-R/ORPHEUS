"""The one definition finder of the retirement engine, for every leg and for the error-catalogue ``Retired:`` leg (1.5)."""

from __future__ import annotations

import ast
from collections.abc import Iterator


def _target_names(target: ast.expr) -> Iterator[ast.Name]:
    if isinstance(target, ast.Name):
        yield target
    elif isinstance(target, (ast.Tuple, ast.List)):
        for element in target.elts:
            yield from _target_names(element)


def definitions(source: str) -> dict[str, int]:
    """Qualname to line of every definition of ``source``.

    A ``def``, a ``class``, a module- or class-level ``Assign`` or ``AnnAssign`` (each ``Name`` of a tuple or list target, C1),
    including those under a module- or class-level ``if``, ``try`` or ``with`` (C3), nested qualnames ``Outer.inner`` with no
    ``<locals>`` (C2), and ``self.x = ...`` anywhere in the body of a method, ``self`` being the method's first parameter (C5), which
    defines ``Class.x``. A variable a function merely assigns is not a definition. The first definition of a qualname wins.
    """
    found: dict[str, int] = {}

    def visit(body: list[ast.stmt], prefix: str, in_class: bool, defs_only: bool = False) -> None:
        for node in body:
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                found.setdefault(prefix + node.name, node.lineno)
                if in_class and node.args.args:
                    receiver = node.args.args[0].arg
                    for sub in ast.walk(node):
                        targets = sub.targets if isinstance(sub, ast.Assign) else [sub.target] if isinstance(sub, ast.AnnAssign) else []
                        for target in targets:
                            if isinstance(target, ast.Attribute) and isinstance(target.value, ast.Name) and target.value.id == receiver:
                                found.setdefault(prefix + target.attr, target.lineno)
                visit(node.body, f"{prefix}{node.name}.", False, True)
            elif isinstance(node, ast.ClassDef):
                found.setdefault(prefix + node.name, node.lineno)
                visit(node.body, f"{prefix}{node.name}.", True)
            elif defs_only:
                continue
            elif isinstance(node, ast.Assign):
                for target in node.targets:
                    for name in _target_names(target):
                        found.setdefault(prefix + name.id, name.lineno)
            elif isinstance(node, ast.AnnAssign):
                for name in _target_names(node.target):
                    found.setdefault(prefix + name.id, name.lineno)
            elif isinstance(node, (ast.If, ast.Try, ast.With)):
                blocks = [node.body, getattr(node, "orelse", []), getattr(node, "finalbody", [])]
                blocks += [handler.body for handler in getattr(node, "handlers", [])]
                for block in blocks:
                    visit(block, prefix, in_class)

    visit(ast.parse(source).body, "", False)
    return found
