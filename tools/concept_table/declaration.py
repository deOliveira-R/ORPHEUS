"""The class declaration: one comment on the physical line that ends a class header. Reader and parser; imports nothing from orpheus.
    # concept: <id> | role: <role> | kind: <kind>[.<subkind>] [| welded: <id>, <id>, ...]
"""
from __future__ import annotations

import ast
import io
import re
import tokenize
from collections.abc import Iterator
from dataclasses import dataclass

ROLES = ("head", "arm", "weld", "part", "component", "interface", "twin", "annex")
KINDS = ("set", "measure", "space", "element", "morphism", "record", "apparatus")
SUBKINDS = {"element": ("member", "covector"), "morphism": ("operator",)}
_MARK = re.compile(r"#\s*concept:")
CONCEPT_ID = re.compile(r"(EXT:)?[a-z][a-z0-9-]*\Z")


class DeclarationError(ValueError):
    """A declaration that does not parse, or is not on the line that ends its class header."""


@dataclass(frozen=True)
class Declaration:
    concept: str
    role: str
    kind: str
    subkind: str | None  # None, "member", "covector" or "operator"
    welded: tuple[str, ...]  # () unless the role is "weld" or the row is a welded head; then welded[0] == concept

    def __post_init__(self) -> None:
        if self.role not in ROLES or self.kind not in KINDS:
            raise DeclarationError(f"role {self.role!r} or kind {self.kind!r} is not in the vocabulary")
        if self.subkind is not None and self.subkind not in SUBKINDS.get(self.kind, ()):
            raise DeclarationError(f"subkind {self.subkind!r} is not allowed after kind {self.kind!r}")
        if not CONCEPT_ID.match(self.concept):
            raise DeclarationError(f"concept id {self.concept!r} is not an id")
        if self.role == "weld" and not self.welded:
            raise DeclarationError("a weld names its concepts")
        if self.welded and self.role not in ("weld", "head"):
            raise DeclarationError("only a weld or a head carries a welded list")
        if self.welded and (self.welded[0] != self.concept or len(self.welded) < 2 or len(set(self.welded)) != len(self.welded)):
            raise DeclarationError("welded starts with the row's concept, holds at least two distinct ids")

    def render(self) -> str:
        kind = self.kind + (f".{self.subkind}" if self.subkind else "")
        text = f"# concept: {self.concept} | role: {self.role} | kind: {kind}"
        return text + (f" | welded: {', '.join(self.welded)}" if self.welded else "")


def parse_declaration(comment: str) -> Declaration:
    """Parse one comment token (the text from '#'). Text after a second '#' is another comment and is ignored."""
    body = comment.split("#")[1]
    keys: list[str] = []
    got: dict[str, str] = {}
    for field in body.split("|"):
        key, sep, value = field.strip().partition(":")  # an id may hold ':' ("EXT:citation"): only the first ':' splits
        if not sep or key not in ("concept", "role", "kind", "welded") or key in got:
            raise DeclarationError(f"field {field.strip()!r} is not 'concept:', 'role:', 'kind:' or 'welded:', or is repeated")
        keys.append(key)
        got[key] = value.strip()
    if keys not in (["concept", "role", "kind"], ["concept", "role", "kind", "welded"]):
        raise DeclarationError(f"fields are concept, role, kind and optionally welded, in that order; got {keys}")
    kind, _, sub = got["kind"].partition(".")
    welded = tuple(x.strip() for x in got["welded"].split(",")) if "welded" in got else ()
    return Declaration(got["concept"], got["role"], kind, sub or None, welded)


def header_end(tokens: list[tokenize.TokenInfo], node: ast.ClassDef) -> int:
    """The physical line holding the colon at bracket depth 0 that ends the header of `node` (from its `class` keyword)."""
    depth, started = 0, False
    for tok in tokens:
        if tok.start[0] < node.lineno:
            continue
        if not started:
            started = tok.type == tokenize.NAME and tok.string == "class" and tok.start[0] == node.lineno
            continue
        if tok.type == tokenize.OP:
            if tok.string in "([{":
                depth += 1
            elif tok.string in ")]}":
                depth -= 1
            elif tok.string == ":" and depth == 0:
                return tok.start[0]
    raise DeclarationError(f"no header colon for the class at line {node.lineno}")


def walk_classes(tree: ast.AST) -> Iterator[tuple[str, ast.ClassDef]]:
    """Every class of a module with its qualified name: ``Outer.Inner`` inside a class, ``outer.<locals>.Name`` inside a function.

    The one qualname walk: the declaration reader, the population, the fingerprint and the reference pass all take their
    classes from it, so no two of them can disagree on what a class is called.
    """
    def walk(n: ast.AST, prefix: str) -> Iterator[tuple[str, ast.ClassDef]]:
        for c in ast.iter_child_nodes(n):
            if isinstance(c, ast.ClassDef):
                yield prefix + c.name, c
                yield from walk(c, prefix + c.name + ".")
            elif isinstance(c, (ast.FunctionDef, ast.AsyncFunctionDef)):
                yield from walk(c, prefix + c.name + ".<locals>.")
            else:
                yield from walk(c, prefix)

    yield from walk(tree, "")


def read_declarations(source: str) -> tuple[dict[str, tuple[int, int, list[str]]], list[int]]:
    """({qualname: (class line, header-end line, [marker comments on the header-end line])}, lines of marker comments that are on no header-end line).
    A marker is a comment token matching `# concept:`. The qualname of a function-local class is `outer.<locals>.Name`."""
    tokens = list(tokenize.generate_tokens(io.StringIO(source).readline))
    markers = {t.start[0]: t.string for t in tokens if t.type == tokenize.COMMENT and _MARK.match(t.string)}
    out: dict[str, tuple[int, int, list[str]]] = {}
    used: set[int] = set()
    for q, c in walk_classes(ast.parse(source)):
        end = header_end(tokens, c)
        out[q] = (c.lineno, end, [markers[end]] if end in markers else [])
        used.add(end)
    return out, sorted(set(markers) - used)
