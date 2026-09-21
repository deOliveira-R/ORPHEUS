"""Plain-text citations by ID, resolved against the registries that define them.

A link with an anchor is checked by the build; ``X2``, ``Cardinal Rule 4``,
``Pattern 7``, ```vv-principles` #17``, ``mode 8``, ``ERR-026``, ``L28``,
``B.4`` (a retirement-audit item) and a plan-authoring tag such as ``VALIDATE-THE-FILTER`` are not, and
the review of 2026-09-21 counted 89 of them over 12 files before the
consolidation that made citing by ID the norm. Each registry below is parsed
from the page that defines its IDs; every citation in the pages under check
must resolve, so a renamed or retired definition reddens ``--check`` at every
site that still cites it.

Not checked, stated so the zero is read for what it is: a bare ``#N`` with no
page name before it in the same paragraph, outside the two pages that number
their own anti-patterns (it reads as an issue number); ``L1``–``L4`` (they are
also the V&V levels); a tag with fewer than two hyphens (indistinguishable
from a word in capitals); the evidence pages (historical text cites retired
numberings); anything inside a code span or a fence.
"""
from __future__ import annotations

import re
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from .source import DOCS_DEV, REPO_ROOT, rel

ERROR_CATALOG = REPO_ROOT / "docs" / "theory" / "verification" / "error_catalog.rst"
ANTI_PATTERN_PAGES = ("coding-elegance", "vv-principles")
_FENCE = re.compile(r"(?ms)^```.*?^```[ \t]*$")
_CODE = re.compile(r"`[^`\n]*`")
_ITEM = re.compile(r"(?m)^(\d+)\. \*\*")


def _text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


@dataclass(frozen=True)
class Registries:
    x: frozenset[str]                       # instrument-doctrine: X1..X4
    cardinal: frozenset[str]                # cardinal: 1..5
    pattern: frozenset[str]                 # coding-elegance: 1..7
    anti: dict[str, frozenset[str]]         # per page: its numbered anti-patterns
    mode: frozenset[str]                    # vv-principles: the 12 failure modes
    err: frozenset[str]                     # error_catalog.rst: ERR-NNN
    lesson: frozenset[str]                  # evidence/lessons.md: Lnn
    item: frozenset[str]                    # retirement-audit: A.1 .. G.25
    tag: frozenset[str]                     # plan-authoring and retirement-audit: the bold tags


def registries(root: Path = DOCS_DEV, catalog: Path = ERROR_CATALOG) -> Registries:
    rules, skills = root / "rules", root / "skills"
    ce, vv = _text(skills / "coding-elegance.md"), _text(skills / "vv-principles.md")
    cs = _text(skills / "retirement-audit.md")
    items: set[str] = set()
    section = ""
    for line in cs.splitlines():
        if m := re.match(r"^\*\*([A-G])\. ", line):
            section = m.group(1)
        elif (m := _ITEM.match(line)) and section:
            items.add(f"{section}.{m.group(1)}")
    modes = set(re.findall(r"(?m)^\| (\d+) \|", vv.split("## The 6 AI failure modes", 1)[1].split("## The 6 test-design", 1)[0]))
    modes |= set(_ITEM.findall(vv.split("## The 6 test-design", 1)[1].split("## Anti-patterns", 1)[0]))
    return Registries(
        x=frozenset(f"X{n}" for n in re.findall(r"(?m)^## X(\d+)\.", _text(rules / "instrument-doctrine.md"))),
        cardinal=frozenset(re.findall(r"(?m)^## (\d)\. ", _text(rules / "cardinal.md"))),
        pattern=frozenset(re.findall(r"(?m)^\*\*(\d) — ", ce)),
        anti={"coding-elegance": frozenset(re.findall(r"(?m)^(\d+)\. \*\*NEVER\*\*", ce)),
              "vv-principles": frozenset(re.findall(r"(?m)^(\d+)\. \*\*NEVER\*\*", vv))},
        mode=frozenset(modes),
        err=frozenset(re.findall(r"error-entry:: (ERR-\d+)", _text(catalog))) if catalog.exists() else frozenset(),
        lesson=frozenset(re.findall(r"(?m)^## (L\d+)\b", _text(root / "evidence" / "lessons.md"))),
        item=frozenset(items),
        tag=frozenset(re.findall(r"(?m)^- \*\*([A-Z][A-Z'\[\]-]+(?: [A-Z][A-Z'-]+)*)\*\*", _text(rules / "plan-authoring.md") + _text(skills / "retirement-audit.md"))),
    )


def pages_under_check(root: Path = DOCS_DEV) -> list[Path]:
    return sorted([*root.glob("rules/*.md"), *root.glob("skills/*.md"), *root.glob("agents/*.md"),
                   root / "onboarding.md", root / "workflows.md", root / "harness.md"])


@dataclass(frozen=True)
class Citation:
    kind: str
    id: str
    line: int


def citations(text: str, page_name: str) -> list[Citation]:
    """Every ID cited in ``text`` outside code, with the registry it names."""
    masked = _FENCE.sub(lambda m: "\n" * m.group(0).count("\n"), text)
    spans = [m.span() for m in _CODE.finditer(masked)]

    def in_code(pos: int) -> bool:
        return any(a <= pos < b for a, b in spans)

    def line_of(pos: int) -> int:
        return masked.count("\n", 0, pos) + 1

    found: list[Citation] = []
    for m in re.finditer(r"\bX(\d)\b", masked):
        if not in_code(m.start()):
            found.append(Citation("x", m.group(0), line_of(m.start())))
    for m in re.finditer(r"Cardinal Rules? (\d)((?:,? (?:and )?\d)*)", masked):
        for n in [m.group(1), *re.findall(r"\d", m.group(2))]:
            found.append(Citation("cardinal", n, line_of(m.start())))
    for m in re.finditer(r"\bPatterns? ?[- ]?(\d)\b((?:(?:,| and| ∩) (\d))*)", masked):
        for n in [m.group(1), *re.findall(r"\d", m.group(2))]:
            found.append(Citation("pattern", n, line_of(m.start())))
    for m in re.finditer(r"\b[Mm]ode[- ](\d{1,2})\b", masked):
        found.append(Citation("mode", m.group(1), line_of(m.start())))
    for m in re.finditer(r"\bERR-(\d{3})\b", masked):
        if not in_code(m.start()):
            found.append(Citation("err", m.group(0), line_of(m.start())))
    for m in re.finditer(r"\bL(\d{1,3})\b", masked):
        if int(m.group(1)) >= 5 and not in_code(m.start()):
            found.append(Citation("lesson", m.group(0), line_of(m.start())))
    for m in re.finditer(r"\b([A-G])\.(\d{1,2})\b", masked):
        if not in_code(m.start()):
            found.append(Citation("item", f"{m.group(1)}.{m.group(2)}", line_of(m.start())))
    for m in re.finditer(r"\b([A-Z][A-Z']*(?:-[A-Z'\[\]]+){2,})\b", masked):
        if not in_code(m.start()):
            found.append(Citation("tag", m.group(1), line_of(m.start())))
    # `#N`: the page's own anti-patterns, or those of the page named last in the same paragraph
    for para_m in re.finditer(r"(?s)[^\n](?:[^\n]|\n(?!\n))*", masked):
        para = para_m.group(0)
        qualifier = page_name if page_name in ANTI_PATTERN_PAGES else None
        for m in re.finditer(r"`(coding-elegance|vv-principles)`|(?<![\w#])#(\d{1,2})\b", para):
            if m.group(1):
                qualifier = m.group(1)
            elif not in_code(para_m.start() + m.start()):
                found.append(Citation(f"anti:{qualifier}", m.group(2), line_of(para_m.start() + m.start())))
    return found


def check(root: Path = DOCS_DEV, pages: Iterable[Path] | None = None, catalog: Path = ERROR_CATALOG) -> tuple[int, list[str]]:
    """(the number of citations resolved, the problems): a citation whose ID no registry defines."""
    reg = registries(root, catalog)
    resolved = 0
    problems: list[str] = []
    for path in (pages if pages is not None else pages_under_check(root)):
        for c in citations(_text(path), path.stem):
            where = f"{rel(path)}:{c.line}"
            if c.kind == "x":
                ok, defined = c.id in reg.x, "instrument-doctrine"
            elif c.kind == "cardinal":
                ok, defined = c.id in reg.cardinal, "cardinal"
            elif c.kind == "pattern":
                ok, defined = c.id in reg.pattern, "coding-elegance"
            elif c.kind == "mode":
                ok, defined = c.id in reg.mode, "vv-principles"
            elif c.kind == "err":
                ok, defined = c.id in reg.err, rel(catalog)
            elif c.kind == "lesson":
                ok, defined = c.id in reg.lesson, "evidence/lessons.md"
            elif c.kind == "item":
                ok, defined = c.id in reg.item, "retirement-audit"
            elif c.kind == "tag":
                ok, defined = c.id in reg.tag, "plan-authoring or retirement-audit"
            elif c.kind == "anti:None":
                problems.append(f"{where}: cites #{c.id} with no page named before it in the paragraph (an anti-pattern needs `vv-principles` or `coding-elegance` beside it; an issue is never a bare #N)")
                continue
            else:
                page = c.kind.split(":", 1)[1]
                ok, defined = c.id in reg.anti[page], page
            if ok:
                resolved += 1
            else:
                problems.append(f"{where}: cites {c.kind.split(':')[0]} {c.id}, which {defined} does not define")
    return resolved, problems
