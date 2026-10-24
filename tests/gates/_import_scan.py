r"""One AST scan of the modules a source imports, for the gates that fence a layer (wave 2, S5).

Not a test module. The gate that says "this file imports nothing from these packages" needs an instrument whose zero is evidence
(``instrument-doctrine`` X1): a ``grep`` for ``from orpheus.transport`` misses ``import orpheus.transport.x``, the
``from orpheus import transport`` spelling, an import inside a function body or under ``TYPE_CHECKING``, and a relative import that
reaches the package without naming it. :func:`imported_names` walks the whole tree with :mod:`ast` (so a function-body and a
``TYPE_CHECKING`` import are found) and returns every dotted name the source imports in EACH spelling it could be matched by: the
module, and, for ``from m import a``, the name ``m.a`` too (``a`` may be a submodule). A relative import is reported as
``"<relative>"`` unless the caller names the importing file's ``package`` (the dotted package a relative import is read from:
``orpheus.derivations.common`` for ``orpheus/derivations/common/x.py``), in which case it is RESOLVED to its absolute name, because a
relative import can reach a forbidden package (``from ...geometry.boundary import X``) without spelling it.

:func:`forbidden_imports` filters those names by dotted prefix (a prefix matches itself and its submodules, never a sibling that
shares leading letters: ``orpheus.transport`` does not match ``orpheus.transport_x``). Each consumer proves its fence with
:data:`CONTROLS`, one planted source per import spelling, each of which must be found; a clean source must not be (the negative
leg).
"""

from __future__ import annotations

import ast
from collections.abc import Iterable

__all__ = ["CONTROLS", "forbidden_imports", "imported_names"]


def imported_names(source: str, package: str | None = None) -> list[str]:
    """Every dotted name ``source`` imports, any spelling, any depth.

    A relative import is resolved against ``package`` when it is given, and is ``"<relative>"`` when it is not.
    """
    names: list[str] = []
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Import):
            names += [alias.name for alias in node.names]
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                if package is None:
                    names.append("<relative>")
                    continue
                parts = package.split(".")
                if node.level - 1 > len(parts) - 1:
                    names.append("<relative>")  # reaches above the top-level package: not resolvable, so not silent
                    continue
                anchor = ".".join(parts[: len(parts) - (node.level - 1)])
                base = f"{anchor}.{node.module}" if node.module else anchor
            else:
                base = node.module or ""
            names += [base] + [f"{base}.{alias.name}" for alias in node.names]
    return names


def forbidden_imports(source: str, prefixes: Iterable[str], package: str | None = None) -> list[str]:
    """The imported names of ``source`` that are, or lie under, one of ``prefixes``; an unresolved relative import is reported too."""
    prefixes = tuple(prefixes)
    return [
        name for name in imported_names(source, package)
        if name == "<relative>" or any(name == p or name.startswith(p + ".") for p in prefixes)
    ]


#: One planted source per import spelling, as ``(source, the prefix that must be found in it)``.
CONTROLS: tuple[tuple[str, str], ...] = (
    ("import orpheus.transport.core", "orpheus.transport"),
    ("from orpheus.transport import method", "orpheus.transport"),
    ("from orpheus import transport", "orpheus.transport"),
    ("from orpheus.geometry.boundary import _tag", "orpheus.geometry.boundary"),
    ("from orpheus.geometry import boundary", "orpheus.geometry.boundary"),
    ("def f():\n    from orpheus.transport.boundary_conditions import WallCondition\n", "orpheus.transport"),
    ("from typing import TYPE_CHECKING\nif TYPE_CHECKING:\n    from orpheus.geometry.boundary import BC\n", "orpheus.geometry.boundary"),
)
