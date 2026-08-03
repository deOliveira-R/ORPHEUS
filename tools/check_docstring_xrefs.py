#!/usr/bin/env python
r"""Resolve Sphinx Python-domain cross-references against the LIVE namespace.

Why this exists — no Sphinx-side gate can do it
------------------------------------------------

An unresolved Python-domain xref (``:func:`` / ``:class:`` / ``:meth:`` /
``:mod:`` / ``:attr:``) renders as **plain text** and draws no ``-W`` warning.
The usual next thought is "run a nitpicky build" — but ``-n`` cannot help
either, and the reason is structural rather than a config gap:

**Sphinx can only diagnose what it RENDERS.** A docstring in a module that is
never ``automodule``'d is never read, so its roles never become nodes and
nothing can warn about them. Measured 2026-08-03 on this tree: the doc source
carries **45** live ``automodule`` directives against **263** modules in
``orpheus/`` — so ~83 % of the package is invisible to Sphinx at every
severity, and ``tests/`` (which carries many xrefs) is never rendered at all.

The project's own knowledge graph inherits the same blind spot for the same
reason. Nexus classifies a reference as ``unresolved`` when it finds no target
NODE, and a target node exists only if Sphinx rendered the target — so live
symbols in un-rendered modules land in the ``unresolved`` bucket alongside
genuinely dead ones. Spot-checked: ``Quadrature.angular_frame``, ``.octants``
and ``.spherical_harmonics`` are all flagged ``unresolved`` and all EXIST.
That bucket therefore cannot be gated on.

This checker sidesteps rendering entirely: it reads docstrings with ``ast``
and resolves each target by **importing it**. Render coverage is irrelevant,
so it covers all 263 modules plus ``tests/``.

What it checks, and what it deliberately does not
-------------------------------------------------

Only **fully-qualified** targets rooted in a known top-level package are
decidable by import, and only those are reported. An unqualified ref
(``:meth:`Quadrature.product```) is resolved by Sphinx against the current
module's context, which this tool does not emulate — flagging those would
manufacture false positives, which is the fastest way to get a gate ignored.

A target is reported ONLY when the failure is a genuinely missing module or
attribute. If the module exists but raises on import, that is reported
separately as a tooling problem, never as a dead reference — the distinction
matters because the two need opposite fixes.

Usage
-----

    .venv/bin/python tools/check_docstring_xrefs.py            # orpheus + tests
    .venv/bin/python tools/check_docstring_xrefs.py orpheus    # narrower

Exit code is 1 when any dead reference is found, so it can gate CI.
"""

from __future__ import annotations

import argparse
import ast
import importlib
import pathlib
import re
from collections import defaultdict

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent

#: ``:role:`target``` — the Python-domain roles that render silently when dead.
ROLE_PATTERN = re.compile(r":(func|class|meth|mod|attr|exc|data|obj):`([^`]+)`")

#: Only refs rooted here are decidable by import. Extend as dependencies grow.
DECIDABLE_ROOTS = ("orpheus", "numpy", "scipy", "sympy", "pytest", "matplotlib")

#: Targets that are legitimately unresolvable — keep EMPTY if at all possible,
#: and always with a comment saying why. An allowlist is where a gate goes to die.
ALLOWLIST: frozenset[str] = frozenset()


def extract_target(raw: str) -> str | None:
    """Normalise a role body to a bare dotted path, or ``None`` if not one.

    Handles the two Sphinx spellings that wrap a target: ``~shortened.path``
    and ``Display Text <real.target>``.
    """
    raw = raw.strip()
    if "<" in raw and raw.endswith(">"):
        raw = raw[raw.index("<") + 1 : -1].strip()
    raw = raw.lstrip("~").lstrip(".")
    if not raw or any(c in raw for c in " \n\t()[]"):
        return None
    return raw


def resolve(dotted: str) -> tuple[bool, str | None]:
    """Resolve ``dotted`` by importing the longest importable prefix and
    ``getattr``-ing the remainder.

    Returns ``(resolved, problem)``. ``problem`` is ``None`` on success,
    ``"import-error:<Type>"`` when a module EXISTS but fails to import (a
    tooling problem, not a dead ref), or ``"missing"`` when the module or
    attribute genuinely is not there.
    """
    parts = dotted.split(".")
    for split in range(len(parts), 0, -1):
        module_name = ".".join(parts[:split])
        try:
            obj = importlib.import_module(module_name)
        except ModuleNotFoundError:
            continue  # try a shorter prefix — this segment may be an attribute
        except Exception as exc:  # exists, but importing it blew up
            return False, f"import-error:{type(exc).__name__}"
        for attribute in parts[split:]:
            try:
                obj = getattr(obj, attribute)
            except AttributeError:
                return False, "missing"
        return True, None
    return False, "missing"


def iter_docstrings(tree: ast.AST):
    """Yield ``(lineno, docstring)`` for every module/class/function docstring.

    Module docstrings report line 1 rather than ``ast.Module``'s absent
    ``lineno``, so every hit is clickable.
    """
    for node in ast.walk(tree):
        if isinstance(
            node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)
        ):
            doc = ast.get_docstring(node, clean=False)
            if doc:
                yield getattr(node, "lineno", 1), doc


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("roots", nargs="*", default=None,
                        help="directories to scan (default: orpheus tests)")
    parser.add_argument("--quiet", action="store_true",
                        help="summary only, no per-site listing")
    args = parser.parse_args()
    roots = args.roots or ["orpheus", "tests"]

    dead: dict[str, list[tuple[str, int, str]]] = defaultdict(list)
    broken_imports: dict[str, str] = {}
    cache: dict[str, tuple[bool, str | None]] = {}
    n_files = n_roles = n_decidable = 0

    for top in roots:
        for path in sorted((REPO_ROOT / top).rglob("*.py")):
            try:
                tree = ast.parse(path.read_text(encoding="utf-8"))
            except (SyntaxError, UnicodeDecodeError):
                continue
            n_files += 1
            for lineno, doc in iter_docstrings(tree):
                for role, raw in ROLE_PATTERN.findall(doc):
                    n_roles += 1
                    target = extract_target(raw)
                    if (
                        target is None
                        or not target.startswith(DECIDABLE_ROOTS)
                        or target in ALLOWLIST
                    ):
                        continue
                    n_decidable += 1
                    if target not in cache:
                        cache[target] = resolve(target)
                    resolved, problem = cache[target]
                    if resolved:
                        continue
                    if problem and problem.startswith("import-error"):
                        broken_imports[target] = problem
                        continue
                    rel = str(path.relative_to(REPO_ROOT))
                    dead[target].append((rel, lineno, role))

    n_sites = sum(len(v) for v in dead.values())
    print(f"files scanned          : {n_files}")
    print(f"xref roles found       : {n_roles}")
    print(f"fully-qualified (gated): {n_decidable}")
    print(f"DEAD TARGETS           : {len(dead)}  across {n_sites} sites")
    if broken_imports:
        print(f"unimportable (TOOLING, not dead refs): {len(broken_imports)}")
        for target, problem in sorted(broken_imports.items()):
            print(f"    {target} -> {problem}")

    if dead and not args.quiet:
        print()
        for target, sites in sorted(dead.items(), key=lambda kv: (-len(kv[1]), kv[0])):
            print(f"  {target}   [{len(sites)} site(s)]")
            for rel, lineno, role in sites:
                print(f"      {rel}:{lineno}  (:{role}:)")

    return 1 if dead else 0


if __name__ == "__main__":
    raise SystemExit(main())
