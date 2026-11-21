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
import warnings
from collections import defaultdict
from collections.abc import Iterable
from pathlib import Path

from tests.gates._import_scan import forbidden_imports, imported_names
from tests.gates._retired import Kept, Moved, Retired, RetiredModule, RetiredName
from tests.gates._retired._finder import definitions

HERE = Path(__file__).resolve().parent
SCAN_TREES = ("orpheus", "tests", "tools", "examples", "student_resources")
_SELF = "tests/gates/test_retirements.py"
_EXCLUDED = "tests/gates/_retired/"
_ENTRY_TYPES = (Retired, Moved, RetiredName, RetiredModule, Kept)
_ATTRIBUTE_SETTERS = ("getattr", "setattr", "delattr")


class GateError(Exception):
    """A malformed rows module or an unparsable file: the gate is red, not silent."""


def load_rows(rows_dir: Path | None = None) -> list:
    """The ``ROWS`` then the ``KEPT`` of every ``rows_*.py`` of ``rows_dir`` (default: this package's directory), in file-name order."""
    directory = HERE if rows_dir is None else Path(rows_dir)
    entries: list = []
    for path in sorted(directory.glob("rows_*.py")):
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
            entries += list(getattr(module, attribute))
    return entries


def validate(entries: list) -> None:
    """Refuse an invalid entry: a non-row, an empty field, a ``file`` or ``home`` not ending ``.py``, a duplicate, a collision."""
    seen: set[tuple] = set()
    for entry in entries:
        if not isinstance(entry, _ENTRY_TYPES):
            raise GateError(f"not a row type: {entry!r}")
        if any(not isinstance(value, str) or value == "" for value in entry):
            raise GateError(f"empty or non-string field in {entry!r}")
        for field in ("file", "home"):
            if field in entry._fields and not getattr(entry, field).endswith(".py"):
                raise GateError(f"{field} must end in .py: {entry!r}")
        key = (type(entry).__name__, tuple(entry))
        if key in seen:
            raise GateError(f"duplicate entry: {entry!r}")
        seen.add(key)
    kept_last = {entry.qualname.split(".")[-1] for entry in entries if isinstance(entry, Kept)}
    for entry in entries:
        if isinstance(entry, RetiredName) and entry.name in kept_last:
            raise GateError(f"RetiredName {entry.name!r} is the last component of a Kept qualname")


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


def _is_str(node: ast.expr, value: str | None = None) -> bool:
    return isinstance(node, ast.Constant) and isinstance(node.value, str) and (value is None or node.value == value)


class _File:
    """One file of the scan, indexed in one pass over its tree."""

    def __init__(self, rel: str, tree: ast.Module) -> None:
        self.rel = rel
        self.tree = tree
        self.package = _package_of(rel)
        self.imports: list[tuple[int, str, list[str]]] = []  # (line, unparsed statement, resolved dotted names)
        self.alias: dict[str, str] = {}  # C6, C15: a name bound by an import binds file-wide
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
                        self.alias[alias.asname] = alias.name
                    else:
                        self.plain_heads.add(alias.name.split(".")[0])
            else:
                for alias in node.names:
                    if alias.name != "*":
                        self.alias[alias.asname or alias.name] = f"{names[0]}.{alias.name}"
        elif isinstance(node, ast.alias):
            self.identifiers[node.name.split(".")[-1]].add(node.lineno)
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
        elif isinstance(node, ast.Assign):
            if any(isinstance(t, ast.Name) and t.id == "__all__" for t in node.targets) and isinstance(node.value, (ast.List, ast.Tuple)):
                for element in node.value.elts:
                    if _is_str(element):
                        self.all_members[element.value].append(element.lineno)  # type: ignore[attr-defined]
                        self.identifiers[element.value].add(element.lineno)  # type: ignore[attr-defined]
        elif isinstance(node, ast.Call):
            self._index_call(node)
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            self.strings.append((node.lineno, node.value))

    def _index_call(self, call: ast.Call) -> None:
        func = call.func
        name = func.id if isinstance(func, ast.Name) else func.attr if isinstance(func, ast.Attribute) else ""
        is_method_setter = isinstance(func, ast.Attribute) and func.attr == "setattr"  # C20: any receiver
        if (isinstance(func, ast.Name) and name in _ATTRIBUTE_SETTERS) or is_method_setter:
            if len(call.args) >= 2 and _is_str(call.args[1]):
                target = call.args[1].value  # type: ignore[attr-defined]
                self.setters[target].append((call.args[0], call.lineno))
                self.identifiers[target].add(call.lineno)
            if is_method_setter and call.args and _is_str(call.args[0]):
                self.setter_targets.append((call.args[0].value, call.lineno))  # type: ignore[attr-defined]
        if name in ("import_module", "__import__") and call.args and _is_str(call.args[0]):
            self.module_strings.append((call.args[0].value, call.lineno))  # type: ignore[attr-defined]

    def resolve(self, expr: ast.expr) -> str | None:
        """The dotted name an expression denotes when its head is bound by an import of this file, else ``None`` (C7)."""
        dotted = _dotted(expr)
        if dotted is None:
            return None
        head, _, rest = dotted.partition(".")
        if head in self.alias:
            return self.alias[head] + ("." + rest if rest else "")
        if head in self.plain_heads:
            return dotted
        return None

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
        self._inits: dict[str, tuple[ast.Module | None, Path]] = {}

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

    def _init(self, package: str) -> tuple[ast.Module | None, Path]:
        if package not in self._inits:
            directory = self.root / package.replace(".", "/")
            init = directory / "__init__.py"
            self._inits[package] = (self._parse(f"{package}/__init__.py", init) if init.is_file() else None, directory)
        return self._inits[package]

    def binding(self, package: str, name: str, target: str, depth: int = 0) -> str:
        """What ``package`` binds ``name`` to: ``"row"`` (the module ``target``, through relative or absolute import links),
        ``"other"`` (a definition, an import from elsewhere, or a submodule of that name) or ``"none"`` (decision E-D1, C18)."""
        init, directory = self._init(package)
        if init is None or depth > 8:
            return "none"
        for node in ast.walk(init):
            if isinstance(node, ast.ImportFrom):
                for alias in node.names:
                    if (alias.asname or alias.name) == name:
                        base = imported_names(ast.unparse(node), package)[0]
                        if alias.asname is None and alias.name == name:
                            if base == target:
                                return "row"
                            if target.startswith(base + "."):
                                return self.binding(base, name, target, depth + 1)
                        return "other"
            elif isinstance(node, ast.Import):
                if any((alias.asname or alias.name.split(".")[0]) == name for alias in node.names):
                    return "other"
        if name in definitions(ast.unparse(init)):
            return "other"
        if (directory / f"{name}.py").exists() or (directory / name).is_dir():
            return "other"
        return "none"

    def valid_ancestors(self, module: str, name: str) -> list[str]:
        """The ancestor packages whose import of ``name`` counts: those that do not bind ``name`` to a different object (E-D1)."""
        return [a for a in _ancestors(module) if self.binding(a, name, module) != "other"]


def _leg_definition(row: Retired | Moved, tree: _Tree) -> set[tuple[str, int]]:
    out: set[tuple[str, int]] = set()
    defined = tree.defs(row.file)
    if defined is not None and row.qualname in defined:
        out.add((row.file, defined[row.qualname]))
    module = _module_of(row.file)
    dotted_qualname = "." in row.qualname
    last = row.qualname.split(".")[-1]
    owner = row.qualname.rsplit(".", 1)[0] if dotted_qualname else None
    ancestors = [] if dotted_qualname else tree.valid_ancestors(module, last)
    for f in tree.files.values():
        if not dotted_qualname:
            for line, _statement, names in f.imports:
                if f"{module}.{last}" in names or any(f"{a}.{last}" in names for a in ancestors):
                    out.add((f.rel, line))
            for value, line in f.attributes.get(last, ()):
                if f.resolve(value) == module:
                    out.add((f.rel, line))
            if last in f.all_members:
                here = _module_of(f.rel)
                if here == module or (f.rel.endswith("__init__.py") and here in ancestors):
                    out.update((f.rel, line) for line in f.all_members[last])
        for first, line in f.setters.get(last, ()):
            resolved = f.resolve(first)
            if (not dotted_qualname and resolved == module) or (
                dotted_qualname and (resolved == f"{module}.{owner}" or any(resolved == f"{a}.{owner}" for a in _ancestors(module)))
            ):
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
    base = row.module.replace(".", "/")
    for candidate in (f"{base}.py", f"{base}/__init__.py"):  # E-D6: the absence leg
        if (tree.root / candidate).is_file():
            out.add((candidate, 0))
    return out


def _leg_kept(entry: Kept, tree: _Tree) -> set[tuple[str, int]]:
    defined = tree.defs(entry.file)
    return set() if defined is not None and entry.qualname in defined else {(entry.file, 0)}


def walk(root: Path, rows: Iterable | None = None, rows_dir: Path | None = None) -> dict[tuple[str, object], list[tuple[str, int]]]:
    """The verdict of every row and ``Kept`` entry on ``root``: ``(type name, row)`` to its sorted ``(file, line)`` findings.

    ``rows`` (rows and ``Kept`` entries mixed) overrides the loader, which reads ``rows_dir`` (default: this package's directory).
    The key carries the type name because ``RetiredName('x') == RetiredModule('x')`` as tuples (C12).
    """
    root = Path(root)
    entries = list(load_rows(rows_dir) if rows is None else rows)
    validate(entries)
    tree = _Tree(root)
    verdicts: dict[tuple[str, object], list[tuple[str, int]]] = {}
    for entry in entries:
        if isinstance(entry, (Retired, Moved)):
            findings = _leg_definition(entry, tree)
        elif isinstance(entry, RetiredName):
            findings = _leg_name(entry, tree)
        elif isinstance(entry, RetiredModule):
            findings = _leg_module(entry, tree)
        else:
            findings = _leg_kept(entry, tree)
        verdicts[(type(entry).__name__, entry)] = sorted(findings)
    return verdicts
