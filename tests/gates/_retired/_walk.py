"""The walk of the shared retirement engine (specification ``ontology_w0_spec.md`` section 1.5).

``walk(root)`` reads every ``.py`` of the five trees under ``root`` once, indexes it, and returns one verdict per row: the list of
``(file, line)`` findings (``line`` is ``0`` for a finding that has no line), empty when the row is green. The rows are the ``ROWS``
and ``KEPT`` of the ``rows_*.py`` modules of THIS package directory, never of ``root`` (decision E-D4: ``root`` is only the source
tree scanned, so a one-file mutation root still carries every row).

The import legs call the one resolver of ``tests.gates._import_scan`` and define no second one.
"""

from __future__ import annotations

import ast
import importlib.util
import sys
import re
import warnings
from collections import defaultdict
from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import TypeGuard, assert_never

from tests.gates._import_scan import forbidden_imports, imported_names
from tests.gates._retired import Kept, Moved, Retired, RetiredModule, RetiredName, Row
from tests.gates._retired._finder import definitions

HERE = Path(__file__).resolve().parent
SCAN_TREES = ("orpheus", "tests", "tools", "examples", "student_resources")
_SELF = "tests/gates/test_retirements.py"
_EXCLUDED = "tests/gates/_retired/"
_ROW_TYPES = (Retired, Moved, RetiredName, RetiredModule)
_REPO_PATH = re.compile(r"^(orpheus|tests|tools|examples|student_resources)/[A-Za-z0-9_/]+\.py$")
Verdicts = dict[tuple[type, object], list[tuple[str, int]]]
_STRING_ATTRIBUTE_ACCESSORS = ("getattr", "setattr", "delattr", "hasattr")  # read or write an attribute by a string: ``hasattr`` answers False silently once it is gone


class GateError(Exception):
    """A malformed rows module or an unparsable file: the gate is red, not silent."""


def load_rows_by_module(rows_dir: Path | None = None) -> dict[str, tuple[tuple[Row, ...], tuple[Kept, ...]]]:
    """The ``(ROWS, KEPT)`` of every ``rows_*.py`` of ``rows_dir`` (default: this package's directory), keyed by file name, in file-name order.

    The two are kept apart and typed: a ``Kept`` inside a ``ROWS`` or a row inside a ``KEPT`` is a gate error (1.5: ``Kept`` is not a
    member of ``ROWS``). The module of an entry is what a failure message names.
    """
    directory = HERE if rows_dir is None else Path(rows_dir)
    if not directory.is_dir():
        raise GateError(f"the rows directory is not a directory: {directory}")
    by_module: dict[str, tuple[tuple[Row, ...], tuple[Kept, ...]]] = {}
    for path in sorted(directory.glob("rows_*.py")):
        rows: list[Row] = []
        kept: list[Kept] = []
        name = f"_retired_rows_{path.stem}_{abs(hash(str(path)))}"
        spec = importlib.util.spec_from_file_location(name, path)
        if spec is None or spec.loader is None:
            raise GateError(f"cannot load the rows module {path}")
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        try:
            spec.loader.exec_module(module)
        finally:
            sys.modules.pop(name, None)
        for attribute in ("ROWS", "KEPT"):
            if not hasattr(module, attribute):
                raise GateError(f"{path.name} defines no {attribute}")
        for entry in module.ROWS:
            if not isinstance(entry, _ROW_TYPES):
                raise GateError(f"{path.name}: ROWS holds a non-row entry (a Kept belongs in KEPT): {entry!r}")
            rows.append(entry)
        for entry in module.KEPT:
            if not isinstance(entry, Kept):
                raise GateError(f"{path.name}: KEPT holds a non-Kept entry (a row belongs in ROWS): {entry!r}")
            kept.append(entry)
        by_module[path.name] = (tuple(rows), tuple(kept))
    return by_module


def load_rows(rows_dir: Path | None = None) -> tuple[tuple[Row, ...], tuple[Kept, ...]]:
    """The ``ROWS`` and the ``KEPT`` of every ``rows_*.py`` of ``rows_dir`` (default: this package's directory), in file-name order:
    the flat pair of :func:`load_rows_by_module`."""
    by_module = load_rows_by_module(rows_dir)
    return tuple(r for rows, _ in by_module.values() for r in rows), tuple(k for _, kept in by_module.values() for k in kept)


def validate(entries: Iterable[object]) -> list[Row | Kept]:
    """Parse ``entries`` into typed rows or refuse: a non-row, an empty field, a ``file`` or ``home`` that is not a repository
    ``.py`` path of the five trees, a row listed twice, a ``(file, qualname)`` named by two of ``Retired``, ``Moved`` and ``Kept``,
    a ``RetiredName`` that is the last component of a ``Kept`` qualname."""
    checked: list[Row | Kept] = []
    seen: set[tuple] = set()
    for entry in entries:
        if not isinstance(entry, (*_ROW_TYPES, Kept)):
            raise GateError(f"not a row type: {entry!r}")
        if any(not isinstance(value, str) or value == "" for value in entry):
            raise GateError(f"empty or non-string field in {entry!r}")
        for field in ("file", "home"):
            if field in entry._fields and not _REPO_PATH.match(getattr(entry, field)):
                raise GateError(f"{field} must be a repository path of the five trees ending in .py (no space, no ./, no absolute path): {entry!r}")
        # one retirement or keeping per (file, qualname) across Retired, Moved and Kept; one entry per fields otherwise
        key = ("definition", entry.file, entry.qualname) if isinstance(entry, (Retired, Moved, Kept)) else (type(entry), tuple(entry))
        if key in seen:
            raise GateError(f"duplicate entry (or one definition named by two of Retired, Moved and Kept): {entry!r}")
        seen.add(key)
        checked.append(entry)
    kept_last = {entry.qualname.split(".")[-1] for entry in checked if isinstance(entry, Kept)}
    for entry in checked:
        if isinstance(entry, RetiredName) and entry.name in kept_last:
            raise GateError(f"RetiredName {entry.name!r} is the last component of a Kept qualname")
    return checked


def _package_of(rel: str) -> str:
    """C4: the package of a file is the dotted path of its directory from the root, ``__init__.py`` included."""
    return ".".join(rel.split("/")[:-1])


def _module_of(rel: str) -> str:
    parts = rel[:-3].split("/")
    if parts[-1] == "__init__":
        parts = parts[:-1]
    return ".".join(parts)


def _ancestors(module: str) -> list[str]:
    parts = module.split(".")
    return [".".join(parts[:i]) for i in range(len(parts) - 1, 0, -1)]


def _dotted(expr: ast.expr) -> str | None:
    if isinstance(expr, ast.Name):
        return expr.id
    if isinstance(expr, ast.Attribute):
        head = _dotted(expr.value)
        return f"{head}.{expr.attr}" if head else None
    return None


def _is_str(node: ast.expr, value: str | None = None) -> TypeGuard[ast.Constant]:
    return isinstance(node, ast.Constant) and isinstance(node.value, str) and (value is None or node.value == value)


class _File:
    """One file of the scan, indexed in one pass over its tree."""

    def __init__(self, rel: str, tree: ast.Module) -> None:
        self.rel = rel
        self.tree = tree
        self.package = _package_of(rel)
        self.imports: list[tuple[int, str, list[str]]] = []  # (line, unparsed statement, resolved dotted names)
        self.alias: dict[str, set[str]] = defaultdict(set)  # C6, C15: a name bound by an import binds file-wide, to every target it is bound to
        self.plain_heads: set[str] = set()  # the head of `import a.b.m` (C19)
        self.identifiers: dict[str, set[int]] = defaultdict(set)  # RetiredName: name -> lines
        self.attributes: dict[str, list[tuple[ast.expr, int]]] = defaultdict(list)  # attr -> (value, line)
        self.all_members: dict[str, list[int]] = defaultdict(list)
        self.setters: dict[str, list[tuple[ast.expr, int]]] = defaultdict(list)  # NAME argument -> (first argument, line)
        self.setter_targets: list[tuple[str, int]] = []  # `<receiver>.setattr("dotted.target", v)`
        self.module_strings: list[tuple[str, int]] = []  # first argument of import_module / __import__
        self.strings: list[tuple[int, str]] = []
        for node in ast.walk(tree):
            self._index(node)
        self._string_imports: list[tuple[int, str, list[str]]] | None = None

    def _index(self, node: ast.AST) -> None:
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            statement = ast.unparse(node)
            names = imported_names(statement, self.package)
            self.imports.append((node.lineno, statement, names))
            if isinstance(node, ast.Import):
                for alias in node.names:
                    if alias.asname:
                        self.alias[alias.asname].add(alias.name)
                    else:
                        self.plain_heads.add(alias.name.split(".")[0])
            else:
                for part in filter(None, (node.module or "").split(".")):
                    self.identifiers[part].add(node.lineno)
                for alias in node.names:
                    if alias.name != "*":
                        self.alias[alias.asname or alias.name].add(f"{names[0]}.{alias.name}")
        elif isinstance(node, ast.alias):
            for part in node.name.split("."):
                self.identifiers[part].add(node.lineno)
            if node.asname:
                self.identifiers[node.asname].add(node.lineno)
        elif isinstance(node, ast.Name):
            self.identifiers[node.id].add(node.lineno)
        elif isinstance(node, ast.Attribute):
            self.identifiers[node.attr].add(node.lineno)
            self.attributes[node.attr].append((node.value, node.lineno))
        elif isinstance(node, (ast.keyword, ast.arg)):
            name = node.arg
            if name:
                self.identifiers[name].add(node.lineno)
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            self.identifiers[node.name].add(node.lineno)
        elif isinstance(node, ast.ExceptHandler):
            if node.name:
                self.identifiers[node.name].add(node.lineno)
        elif isinstance(node, (ast.Global, ast.Nonlocal)):
            for name in node.names:
                self.identifiers[name].add(node.lineno)
        elif isinstance(node, (ast.MatchAs, ast.MatchStar, ast.TypeVar, ast.ParamSpec, ast.TypeVarTuple)):
            if node.name:
                self.identifiers[node.name].add(node.lineno)
        elif isinstance(node, ast.MatchMapping):
            if node.rest:
                self.identifiers[node.rest].add(node.lineno)
        elif isinstance(node, (ast.Assign, ast.AnnAssign, ast.AugAssign)):
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            if any(isinstance(t, ast.Name) and t.id == "__all__" for t in targets) and isinstance(node.value, (ast.List, ast.Tuple)):
                for element in node.value.elts:
                    if _is_str(element):
                        self.all_members[element.value].append(element.lineno)
                        self.identifiers[element.value].add(element.lineno)
        elif isinstance(node, ast.Call):
            self._index_call(node)
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            self.strings.append((node.lineno, node.value))

    def _index_call(self, call: ast.Call) -> None:
        func = call.func
        name = func.id if isinstance(func, ast.Name) else func.attr if isinstance(func, ast.Attribute) else ""
        is_method_setter = isinstance(func, ast.Attribute) and func.attr == "setattr"  # C20: any receiver
        if (isinstance(func, ast.Name) and name in _STRING_ATTRIBUTE_ACCESSORS) or is_method_setter:
            if len(call.args) >= 2 and _is_str(call.args[1]):
                target = call.args[1].value
                self.setters[target].append((call.args[0], call.lineno))
                self.identifiers[target].add(call.lineno)
            if is_method_setter and call.args and _is_str(call.args[0]):
                self.setter_targets.append((call.args[0].value, call.lineno))
        if name in ("import_module", "__import__") and call.args and _is_str(call.args[0]):
            self.module_strings.append((call.args[0].value, call.lineno))

    def resolve(self, expr: ast.expr) -> frozenset[str]:
        """The dotted names an expression denotes when its head is bound by an import of this file (every binding of the head, C6), else none (C7)."""
        dotted = _dotted(expr)
        if dotted is None:
            return frozenset()
        head, _, rest = dotted.partition(".")
        suffix = "." + rest if rest else ""
        denoted = {target + suffix for target in self.alias.get(head, ())}
        if head in self.plain_heads:
            denoted.add(dotted)
        return frozenset(denoted)

    def string_imports(self) -> list[tuple[int, str, list[str]]]:
        """C16, C17: every string constant that parses as Python, read for its imports only (relative ones against this file's package)."""
        if self._string_imports is None:
            found: list[tuple[int, str, list[str]]] = []
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                for line, text in self.strings:
                    try:
                        parsed = ast.parse(text)
                    except (SyntaxError, ValueError, RecursionError):
                        continue
                    for node in ast.walk(parsed):
                        if isinstance(node, (ast.Import, ast.ImportFrom)):
                            statement = ast.unparse(node)
                            found.append((line, statement, imported_names(statement, self.package)))
            self._string_imports = found
        return self._string_imports


class _Tree:
    """The five trees of a root, parsed once, and the definitions of any file under the root."""

    def __init__(self, root: Path) -> None:
        self.root = root
        self.files: dict[str, _File] = {}
        for top in SCAN_TREES:
            for path in sorted((root / top).rglob("*.py")):
                rel = path.relative_to(root).as_posix()
                if "__pycache__" in path.parts or rel.startswith(_EXCLUDED) or rel == _SELF:
                    continue
                self.files[rel] = _File(rel, self._parse(rel, path))
        self._defs: dict[str, dict[str, int] | None] = {}

    @staticmethod
    def _parse(rel: str, path: Path) -> ast.Module:
        try:
            return ast.parse(path.read_text(encoding="utf-8"))
        except (SyntaxError, UnicodeDecodeError, ValueError) as error:
            raise GateError(f"unparsable file under the root: {rel}: {error}") from error

    def defs(self, rel: str) -> dict[str, int] | None:
        """The definitions of ``rel`` (any file under the root), or ``None`` when it is absent."""
        if rel not in self._defs:
            path = self.root / rel
            if not path.is_file():
                self._defs[rel] = None
            else:
                try:
                    self._defs[rel] = definitions(path.read_text(encoding="utf-8"))
                except (SyntaxError, UnicodeDecodeError, ValueError) as error:
                    raise GateError(f"unparsable file under the root: {rel}: {error}") from error
        return self._defs[rel]

    def _read(self, module: str) -> tuple[str, ast.Module] | None:
        """The package ``__init__`` or the ``.py`` of ``module`` under the root, as ``(rel, tree)``: the indexed tree when the scan holds it."""
        base = module.replace(".", "/")
        for rel in (f"{base}/__init__.py", f"{base}.py"):
            if rel in self.files:
                return rel, self.files[rel].tree
            if (self.root / rel).is_file():
                return rel, self._parse(rel, self.root / rel)
        return None

    def binding(self, module: str, name: str, target: str, depth: int = 0) -> str:
        """What ``module`` (a package or a module) binds ``name`` to: ``"row"`` (the module ``target``, through a chain of
        ``from X import name`` links into any module, relative or absolute), ``"other"`` (a definition, a renamed import, an import
        of another object or of an external module, or a submodule of that name) or ``"none"`` (nothing, or a link that leads to a
        module of the five trees that is gone: the leftover of a retirement counts; decision E-D1, C18). Only the module-level
        statements are read, and an ``as`` that repeats the name is no rename."""
        read = self._read(module)
        if read is None:
            return "none" if module.split(".")[0] in SCAN_TREES else "other"
        if depth > 8:
            return "none"
        rel, tree = read
        result = "none"
        for node in _module_level(tree.body):
            if isinstance(node, ast.ImportFrom):
                for alias in node.names:
                    if (alias.asname or alias.name) != name:
                        continue
                    base = imported_names(ast.unparse(node), _package_of(rel))[0]
                    if base == _package_of(rel) or alias.name != name:  # `from . import name`, or a rename: another object
                        result = "other"
                    else:
                        result = "row" if base == target else self.binding(base, name, target, depth + 1)
            elif isinstance(node, ast.Import):
                if any((alias.asname or alias.name.split(".")[0]) == name for alias in node.names):
                    result = "other"
        if result != "none":
            return result
        if name in (self.defs(rel) or {}):
            return "other"
        directory = self.root / module.replace(".", "/")
        if rel.endswith("__init__.py") and ((directory / f"{name}.py").exists() or (directory / name).is_dir()):
            return "other"
        return "none"

    def valid_ancestors(self, module: str, name: str) -> list[str]:
        """The ancestor packages whose import of ``name`` counts: those that do not bind ``name`` to a different object (E-D1)."""
        return [a for a in _ancestors(module) if self.binding(a, name, module) != "other"]


def _module_level(body: list[ast.stmt]) -> Iterator[ast.stmt]:
    """The statements of a module body in source order, descending into ``if``, ``try`` (and ``try``/``except*``) and ``with`` blocks and never into a ``def`` or ``class``."""
    for node in body:
        yield node
        if isinstance(node, (ast.If, ast.Try, ast.With, ast.TryStar)):
            blocks = [node.body, getattr(node, "orelse", []), getattr(node, "finalbody", []), *(h.body for h in getattr(node, "handlers", []))]
            for block in blocks:
                yield from _module_level(block)


def _leg_definition(row: Retired | Moved, tree: _Tree) -> set[tuple[str, int]]:
    out: set[tuple[str, int]] = set()
    defined = tree.defs(row.file)
    if defined is not None and row.qualname in defined:
        out.add((row.file, defined[row.qualname]))
    module = _module_of(row.file)
    dotted_qualname = "." in row.qualname
    last = row.qualname.split(".")[-1]
    owner = row.qualname.rsplit(".", 1)[0] if dotted_qualname else None
    ancestors = tree.valid_ancestors(module, row.qualname.split(".")[0])  # E-D1: the import, ``__all__``, attribute and setter legs alike
    readers = {module, *ancestors}  # E-D1: one rule for when an ancestor package counts, read by the attribute leg and the undotted setter leg
    owners = {f"{module}.{owner}", *(f"{a}.{owner}" for a in ancestors)} if dotted_qualname else readers  # what an attribute read or a string accessor must denote
    for f in tree.files.values():
        if not dotted_qualname:
            for line, _statement, names in f.imports:
                if f"{module}.{last}" in names or any(f"{a}.{last}" in names for a in ancestors):
                    out.add((f.rel, line))
            if last in f.all_members:
                here = _module_of(f.rel)
                if here == module or (f.rel.endswith("__init__.py") and here in ancestors):
                    out.update((f.rel, line) for line in f.all_members[last])
        for value, line in f.attributes.get(last, ()):  # ``<expr>.last``: the module or an ancestor, or for a dotted row the class (an instance read stays unread)
            if owners & f.resolve(value):
                out.add((f.rel, line))
        for first, line in f.setters.get(last, ()):
            if f.resolve(first) & owners:
                out.add((f.rel, line))
        out.update((f.rel, line) for target, line in f.setter_targets if target == f"{module}.{row.qualname}")
    if not dotted_qualname:
        for f in tree.files.values():
            for line, statement, names in f.string_imports():
                hit = f"{module}.{last}" in names or any(f"{a}.{last}" in names for a in ancestors)
                if not hit:  # C22: a relative fixture import whose module is a dotted suffix of the row's module
                    hit = _relative_suffix_hit(statement, last, module)
                if hit:
                    out.add((f.rel, line))
    if isinstance(row, Moved):
        home = tree.defs(row.home)
        if home is None or row.home_name not in home:
            out.add((row.home, 0))
    return out


def _relative_suffix_hit(statement: str, last: str, module: str) -> bool:
    if not statement.startswith("from ."):
        return False
    node = ast.parse(statement).body[0]
    if not (isinstance(node, ast.ImportFrom) and node.level and node.module):
        return False
    return any(alias.name == last for alias in node.names) and (module == node.module or module.endswith("." + node.module))


def _leg_name(row: RetiredName, tree: _Tree) -> set[tuple[str, int]]:
    return {(f.rel, line) for f in tree.files.values() for line in f.identifiers.get(row.name, ())}


def _leg_module(row: RetiredModule, tree: _Tree) -> set[tuple[str, int]]:
    out: set[tuple[str, int]] = set()
    for f in tree.files.values():
        for line, statement, _names in f.imports:
            if forbidden_imports(statement, [row.module], f.package):
                out.add((f.rel, line))
        out.update((f.rel, line) for value, line in f.module_strings if value == row.module or value.startswith(row.module + "."))
    base = row.module.replace(".", "/")  # E-D6: the absence leg: the module file, or any .py under the package directory
    if (tree.root / f"{base}.py").is_file():
        out.add((f"{base}.py", 0))
    directory = tree.root / base
    if directory.is_dir():
        out.update((p.relative_to(tree.root).as_posix(), 0) for p in directory.rglob("*.py") if "__pycache__" not in p.parts)
    return out


def _leg_kept(entry: Kept, tree: _Tree) -> set[tuple[str, int]]:
    defined = tree.defs(entry.file)
    return set() if defined is not None and entry.qualname in defined else {(entry.file, 0)}


def walk(root: Path, rows: Iterable[object] | None = None, rows_dir: Path | None = None) -> Verdicts:
    """The verdict of every row and ``Kept`` entry on ``root``: ``(type of the row, row)`` to its sorted ``(file, line)`` findings.

    ``rows`` (rows and ``Kept`` entries mixed) overrides the loader, which reads ``rows_dir`` (default: this package's directory).
    The key carries the type because ``RetiredName('x') == RetiredModule('x')`` as tuples (C12).
    """
    root = Path(root)
    if not root.is_dir():
        raise GateError(f"the root to walk is not a directory: {root}")
    if rows is None:
        loaded_rows, loaded_kept = load_rows(rows_dir)
        rows = [*loaded_rows, *loaded_kept]
    entries = validate(rows)
    tree = _Tree(root)
    verdicts: Verdicts = {}
    for entry in entries:
        match entry:
            case Retired() | Moved():
                findings = _leg_definition(entry, tree)
            case RetiredName():
                findings = _leg_name(entry, tree)
            case RetiredModule():
                findings = _leg_module(entry, tree)
            case Kept():
                findings = _leg_kept(entry, tree)
            case _:
                assert_never(entry)
        verdicts[(type(entry), entry)] = sorted(findings)
    return verdicts
