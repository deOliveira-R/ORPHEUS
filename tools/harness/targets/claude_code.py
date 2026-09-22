"""Claude Code: ``.claude/`` and the repo-root ``CLAUDE.md``.

* a **rule** is ``.claude/rules/<name>.md``, always-on unless it carries
  ``paths`` (then Claude Code loads it when a matching file is touched);
* a **skill** is ``.claude/skills/<name>/SKILL.md`` with its Agent-Skills front
  matter verbatim, and the error-index injection point becomes a ``!cat`` line
  that loads ``error_index.md`` at skill-load time (harness-specific, so it
  lives here and not in the docs);
* an **agent** is the definition block of ``.claude/agents/<name>/AGENT.md``,
  everything after the hand-maintained front matter (tools, model, memory),
  between markers: the role, the method and the return contract;
* the **index** is ``.claude/<name>.md``, not auto-loaded (the session-start
  batch reads it);
* the **onboarding** page is the on-boarding block of ``CLAUDE.md`` — the one
  target outside ``.claude/``, because the harness reads it at the root; its
  markers were placed by hand once (the file has no front matter to insert after).

Always-on, what every session and every inheriting dispatch pays: the rules
without ``paths`` and CLAUDE.md, generated or installed. An **installed** file
is one another tool wrote into ``.claude/`` and recorded in its manifest
(``nexus setup`` and ``nexus-install-manifest.json`` today); Claude Code loads
a rule by its location whoever wrote it, so an installed rule without
``paths`` is always-on and the cost line counts it, while its drift instrument
is the installer's own check, never this generator.
"""
from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import assert_never

from ..render import Block, Rendered, WholeFile, stamp
from ..source import ERROR_INDEX_MARK, REPO_ROOT, Kind, Page

ERROR_INDEX_LINE = (
    '!`cat "${CLAUDE_PROJECT_DIR:-.}/.claude/skills/vv-principles/error_index.md" '
    '2>/dev/null || echo "(error index unavailable — run: '
    '.venv/bin/python -m tools.verification.generate_error_index)"; exit 0`'
)
_ALL_KINDS = frozenset(Kind)
MANIFESTS = ("nexus-install-manifest.json",)   # under .claude/: {"files": {"<relative path>": {"version": …}}}
_PATHS_FRONT_MATTER = re.compile(r"\A---\n(?:(?!---\n).*\n)*?paths:")


@dataclass(frozen=True)
class ClaudeCode:
    repo_root: Path = REPO_ROOT

    @property
    def name(self) -> str:
        return "claude-code"

    @property
    def kinds(self) -> frozenset[Kind]:
        return _ALL_KINDS

    @property
    def harness_dir(self) -> Path:
        return self.repo_root / ".claude"

    @property
    def roots(self) -> tuple[Path, ...]:
        return (self.harness_dir, self.repo_root / "CLAUDE.md")

    def target(self, page: Page) -> Path:
        match page.kind:
            case Kind.RULE:
                return self.harness_dir / "rules" / f"{page.name}.md"
            case Kind.SKILL:
                return self.harness_dir / "skills" / page.name / "SKILL.md"
            case Kind.AGENT:
                return self.harness_dir / "agents" / page.name / "AGENT.md"
            case Kind.INDEX:
                return self.harness_dir / f"{page.name}.md"
            case Kind.ONBOARDING:
                return self.repo_root / "CLAUDE.md"
            case _:
                assert_never(page.kind)

    def render(self, page: Page, body: str) -> Rendered:
        match page.kind:
            case Kind.RULE:
                return WholeFile(_paths_front_matter(page.paths) + stamp(page.rel) + body.lstrip("\n"))
            case Kind.SKILL:
                return WholeFile(page.front_matter + stamp(page.rel)
                                 + body.lstrip("\n").replace(ERROR_INDEX_MARK, ERROR_INDEX_LINE))
            case Kind.INDEX:
                return WholeFile(stamp(page.rel) + body.lstrip("\n"))
            case Kind.AGENT:
                return Block("definition", body.strip("\n"))
            case Kind.ONBOARDING:
                return Block("on-boarding block", body.strip("\n"))
            case _:
                assert_never(page.kind)

    def always_on(self, page: Page) -> bool:
        match page.kind:
            case Kind.RULE:
                return not page.paths
            case Kind.ONBOARDING:
                return True
            case Kind.SKILL | Kind.AGENT | Kind.INDEX:
                return False
            case _:
                assert_never(page.kind)

    def installed(self) -> dict[Path, str]:
        """Every file an installer wrote under ``.claude/`` and recorded: path → ``<installer> <version>``."""
        found: dict[Path, str] = {}
        for name in MANIFESTS:
            manifest = self.harness_dir / name
            if not manifest.exists():
                continue
            files = json.loads(manifest.read_text(encoding="utf-8")).get("files") or {}
            installer = name.split("-", 1)[0]
            for key, entry in files.items():
                found[self.harness_dir / key] = f"{installer} {(entry or {}).get('version', '?')}"
        return found

    def installed_always_on(self) -> tuple[Path, ...]:
        rules = self.harness_dir / "rules"
        return tuple(sorted(
            p for p in self.installed()
            if p.parent == rules and p.suffix == ".md" and p.exists()
            and not _PATHS_FRONT_MATTER.match(p.read_text(encoding="utf-8"))
        ))


def _paths_front_matter(paths: tuple[str, ...]) -> str:
    """Claude Code's path-scoped rule: a ``paths:`` list in the rule's front matter."""
    if not paths:
        return ""
    return "---\npaths:\n" + "".join(f'  - "{p}"\n' for p in paths) + "---\n\n"
