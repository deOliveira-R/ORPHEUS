"""What a source page is.

A page under ``docs/development/`` that a harness renders declares itself in a
``harness:`` block of its YAML front matter, in block style::

    ---
    harness:
      kind: rule            # rule | skill | agent | index | onboarding
      budget_tokens: 1200
      paths: ["tests/**"]   # optional: a rule that loads only for these paths
    ---

Discovery is path-scoped, so an omission is a PROBLEM and never a silent
retirement: every ``.md`` under ``rules/``, ``skills/`` and ``agents/`` MUST
carry the block, and its ``kind`` MUST agree with the directory; a top-level
page may carry ``index`` or ``onboarding``; any other page (``evidence/``)
must not carry the block. The rest of the front matter is kept verbatim for a
skill (the Agent Skills standard's own keys, emitted as written) and is a
PROBLEM on any other kind. Two pages of one kind and name cannot exist: one
directory, one stem.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
DOCS_DEV = REPO_ROOT / "docs" / "development"
ERROR_INDEX_MARK = "<!-- harness: error_index -->"  # a skill's injection point for the error-catalogue index; the harness realises it


class Kind(StrEnum):
    RULE = "rule"
    SKILL = "skill"
    AGENT = "agent"
    INDEX = "index"
    ONBOARDING = "onboarding"


KIND_DIRS = {"rules": Kind.RULE, "skills": Kind.SKILL, "agents": Kind.AGENT}
TOP_LEVEL_KINDS = frozenset({Kind.INDEX, Kind.ONBOARDING})
BLOCK_KEYS = frozenset({"kind", "budget_tokens", "paths"})
_FRONT = re.compile(r"\A---\n(.*?)\n---\n", re.S)
_BLOCK = re.compile(r"(?m)^harness:[^\n]*\n(?:[ \t]+[^\n]*\n|[ \t]*\n(?=[ \t]+\S))*")  # the key line, its indented continuation, blank lines inside


@dataclass(frozen=True)
class Page:
    kind: Kind
    name: str
    path: Path
    rel: str            # the path as the stamp and every problem name it: relative to the repo root
    front_matter: str   # the page's own front matter, fences included, the harness block removed; "" when nothing remains
    body: str
    budget_tokens: int
    paths: tuple[str, ...] = ()


def rel(path: Path) -> str:
    return path.relative_to(REPO_ROOT).as_posix() if path.is_relative_to(REPO_ROOT) else path.as_posix()


def _load(text: str, where: str) -> tuple[dict | None, str | None]:
    """Front matter as a mapping, or the problem: not YAML, or not a mapping."""
    try:
        meta = yaml.safe_load(text) or {}
    except yaml.YAMLError as e:
        return None, f"{where}: front matter is not YAML ({e})"
    if not isinstance(meta, dict):
        return None, f"{where}: front matter is not a mapping"
    return meta, None


def read(path: Path, root: Path = DOCS_DEV) -> tuple[Page | None, list[str]]:
    """Parse one page. A page that is not generated returns ``(None, [])``; a
    malformed one ``(None, [problem])``."""
    parts = path.relative_to(root).parts
    first = parts[0] if parts else ""
    where = rel(path)
    text = path.read_text(encoding="utf-8")
    m = _FRONT.match(text)
    front_text = (m.group(1) + "\n") if m else ""
    meta: dict = {}
    if m:
        loaded, problem = _load(front_text, where)
        if loaded is None:
            return None, [problem or ""]
        meta = loaded
    if len(parts) > 2 and first in KIND_DIRS:
        return None, [f"{where}: a page under {first}/ is generated one level deep only; nested pages are not"]
    expected = KIND_DIRS.get(first) if len(parts) == 2 else None
    top_level = len(parts) == 1
    block = meta.get("harness")
    if block is None:
        if expected is not None:
            return None, [f"{where}: a page under {first}/ must declare a harness: block (kind, budget_tokens)"]
        return None, []
    if expected is None and not top_level:
        return None, [f"{where}: a harness: block is allowed only under rules/, skills/, agents/ or at the top level"]
    if not isinstance(block, dict):
        return None, [f"{where}: harness: must be a mapping"]
    if unknown := set(block) - BLOCK_KEYS:
        return None, [f"{where}: harness: has unknown keys {sorted(unknown)}"]
    try:
        kind = Kind(block.get("kind"))
    except ValueError:
        return None, [f"{where}: harness.kind {block.get('kind')!r} is not one of {[k.value for k in Kind]}"]
    if expected is not None and kind is not expected:
        return None, [f"{where}: harness.kind {kind.value} disagrees with its directory {first}/ ({expected.value})"]
    if top_level and kind not in TOP_LEVEL_KINDS:
        return None, [f"{where}: a top-level page may be {sorted(k.value for k in TOP_LEVEL_KINDS)}, not {kind.value}"]
    budget = block.get("budget_tokens")
    if not isinstance(budget, int) or isinstance(budget, bool):
        return None, [f"{where}: harness.budget_tokens must be an integer"]
    paths = block.get("paths", [])
    if not isinstance(paths, list) or not all(isinstance(p, str) for p in paths):
        return None, [f"{where}: harness.paths must be a list of strings"]
    remaining = _BLOCK.sub("", front_text)
    rest, problem = _load(remaining, where)
    if rest is None or rest != {k: v for k, v in meta.items() if k != "harness"}:
        return None, [f"{where}: the harness: block could not be removed textually (write it in block style, no anchors); {problem or 'the rest re-parses differently'}"]
    front_matter = f"---\n{remaining.rstrip()}\n---\n" if remaining.strip() else ""
    if front_matter and kind is not Kind.SKILL:
        return None, [f"{where}: front matter other than harness: is meaningful only on a skill"]
    if kind is Kind.SKILL and not {"name", "description"} <= set(rest):
        return None, [f"{where}: a skill source needs the Agent Skills front matter (name, description) beside its harness: block"]
    body = text[m.end():] if m else text
    return Page(kind, path.stem, path, where, front_matter, body, budget, tuple(paths)), []


def discover(root: Path = DOCS_DEV) -> tuple[list[Page], list[str]]:
    """Every generated page under ``root``, and every problem found on the way."""
    pages: list[Page] = []
    problems: list[str] = []
    for path in sorted(root.rglob("*.md")):
        page, found = read(path, root)
        problems += found
        if page is not None:
            pages.append(page)
    return pages, problems
