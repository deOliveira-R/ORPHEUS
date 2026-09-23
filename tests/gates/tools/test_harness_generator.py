"""One red per problem arm of ``tools/harness`` (instrument-doctrine X1: an
instrument is evidence only if some realizable state changes its reading).
``tests/gates/test_harness_generated.py`` runs the real tree through ``--check``;
these run each pure function on a fault it must catch, the Claude Code
harness on every kind, the pipeline on both shapes, and the neutrality
invariant with its positive controls.
"""
from __future__ import annotations

import ast
import json
import pathlib

import pytest

from tools.harness import budget, paths, pipeline, source
from tools.harness.links import relink
from tools.harness.render import Block, WholeFile, block_text, markers, splice, stamp
from tools.harness.source import Kind, Page, discover
from tools.harness.targets.claude_code import ERROR_INDEX_LINE, ClaudeCode

pytestmark = pytest.mark.foundation

REPO_ROOT = pathlib.Path(__file__).resolve().parents[3]
PKG = REPO_ROOT / "tools" / "harness"


def write(path: pathlib.Path, text: str) -> pathlib.Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def front(kind: str, tokens: int = 100, extra: str = "") -> str:
    return f"---\nharness:\n  kind: {kind}\n  budget_tokens: {tokens}\n{extra}---\n\n"


SKILL_FRONT = "---\nname: s\ndescription: d\nharness:\n  kind: skill\n  budget_tokens: 400\nallowed-tools: Bash\n---\n\n# S\n"


def page(kind: Kind, name: str = "p", body: str = "# P\n\nbody\n", budget_tokens: int = 100,
         front_matter: str = "", paths: tuple[str, ...] = (), path: pathlib.Path | None = None) -> Page:
    path = path or REPO_ROOT / "docs" / "development" / f"{name}.md"
    return Page(kind, name, path, source.rel(path), front_matter, body, budget_tokens, paths)


# ---------------------------------------------------------------- budget

def test_budget_over_is_a_problem() -> None:
    assert "> budget" in (budget.check("x" * 3600, 999, "t") or "")


def test_budget_slack_over_is_a_problem() -> None:
    assert "SLACK_MAX" in (budget.check("x" * 360, 100 + budget.SLACK_MAX + 1, "t") or "")


def test_budget_within_slack_is_fine() -> None:
    assert budget.check("x" * 360, 100 + budget.SLACK_MAX, "t") is None


# ----------------------------------------------------------------- links

def test_relink_repoints_and_flags_missing_file_and_anchor(tmp_path: pathlib.Path) -> None:
    write(tmp_path / "docs" / "a.md", "# Title\n\ntext\n")
    src = write(tmp_path / "docs" / "sub" / "b.md", "[t](../a.md#title) [x](../a.md#nope) [m](../missing.md)\n")
    dst = tmp_path / "out" / "deep" / "b.md"
    problems: list[str] = []
    out = relink(src.read_text(), src, dst, {}, problems, "b.md")
    assert "](../../docs/a.md#title)" in out
    assert [p for p in problems if "missing heading" in p and "#nope" in p]
    assert [p for p in problems if "missing file" in p and "missing.md" in p]
    assert "](../missing.md)" in out, "a link that cannot be resolved is left as written"


def test_relink_writes_a_generated_source_to_its_copy(tmp_path: pathlib.Path) -> None:
    a = write(tmp_path / "docs" / "a.md", "# A\n")
    src = write(tmp_path / "docs" / "b.md", "[a](a.md)\n")
    dst, a_copy = tmp_path / "out" / "b.md", tmp_path / "out" / "a.md"
    out = relink(src.read_text(), src, dst, {a.resolve(): a_copy.resolve()}, [], "b.md")
    assert out == "[a](a.md)\n"


# ---------------------------------------------------------------- render

BLOCK = Block("definition", "# qa — role block\n\nrole")


def test_splice_replaces_between_markers_and_keeps_the_rest() -> None:
    begin, end = markers("definition", "s.md")
    current = f"---\nname: qa\n---\n\n{begin}\nOLD\n{end}\nkept tail\n"
    out, problem = splice(BLOCK, current, "s.md", "t")
    assert problem is None and "OLD" not in out
    assert out == f"---\nname: qa\n---\n\n{block_text(BLOCK, 's.md')}kept tail\n"


def test_splice_does_not_eat_a_character_after_the_end_marker() -> None:
    begin, end = markers("definition", "s.md")
    out, problem = splice(BLOCK, f"{begin}\nOLD\n{end}# next", "s.md", "t")
    assert problem is None and out.endswith("-->\n# next")


def test_splice_inserts_after_front_matter_on_first_run() -> None:
    out, problem = splice(BLOCK, "---\nname: qa\n---\ntail\n", "s.md", "t")
    assert problem is None and out.startswith("---\nname: qa\n---\n\n<!-- BEGIN GENERATED definition")


@pytest.mark.parametrize("current, fragment", [
    (f"---\nx: 1\n---\n{markers('definition', 's.md')[0]}\nbody\n", "BEGIN marker without END"),
    (f"---\nx: 1\n---\nbody\n{markers('definition', 's.md')[1]}\n", "END marker without BEGIN"),
    ("---\nx: 1\n---\n# qa — role block\n\nrole\n", "present without its markers"),
    ("no front matter\n", "no front matter to insert"),
])
def test_splice_faults_are_problems(current: str, fragment: str) -> None:
    out, problem = splice(BLOCK, current, "s.md", "t")
    assert out == current and problem and fragment in problem


def test_stamp_names_the_generator() -> None:
    assert stamp("docs/development/rules/x.md").startswith("<!-- GENERATED by tools/harness from docs/development/rules/x.md")


# ---------------------------------------------------------------- source

def test_discover_reads_every_kind_and_strips_only_the_block(tmp_path: pathlib.Path) -> None:
    root = tmp_path / "dev"
    write(root / "rules" / "r.md", front("rule", 300) + "# R\n")
    write(root / "skills" / "s.md", SKILL_FRONT)
    write(root / "agents" / "a.md", front("agent") + "# a — role block\n")
    write(root / "lessons.md", front("index") + "# L\n")
    write(root / "onboarding.md", front("onboarding", 1500, "  paths: []\n") + "# O\n")
    write(root / "evidence" / "e.md", "# evidence, no block\n")
    write(root / "harness.md", "# top-level page, not generated\n")
    pages, problems = discover(root)
    assert problems == []
    assert {p.name: p.kind for p in pages} == {"r": Kind.RULE, "s": Kind.SKILL, "a": Kind.AGENT, "lessons": Kind.INDEX, "onboarding": Kind.ONBOARDING}
    skill = next(p for p in pages if p.kind is Kind.SKILL)
    assert skill.front_matter == "---\nname: s\ndescription: d\nallowed-tools: Bash\n---\n"
    assert skill.body == "\n# S\n" and skill.budget_tokens == 400
    assert next(p for p in pages if p.kind is Kind.RULE).front_matter == ""


@pytest.mark.parametrize("style", [
    "---\nharness: {kind: rule, budget_tokens: 100}\n---\n# R\n",                    # one-line mapping
    "---\nharness:\n  kind: rule\n\n  budget_tokens: 100\n---\n# R\n",               # a blank line inside the block
])
def test_block_styles_removable_textually_are_accepted(tmp_path: pathlib.Path, style: str) -> None:
    write(tmp_path / "dev" / "rules" / "x.md", style)
    pages, problems = discover(tmp_path / "dev")
    assert problems == [] and pages[0].budget_tokens == 100 and pages[0].front_matter == ""


@pytest.mark.parametrize("rel, text, fragment", [
    ("rules/x.md", "# no block\n", "must declare a harness: block"),
    ("rules/x.md", front("skill") + "# R\n", "disagrees with its directory"),
    ("evidence/x.md", front("rule") + "# E\n", "allowed only under"),
    ("x.md", front("rule") + "# top\n", "a top-level page may be"),
    ("rules/x.md", front("law") + "# R\n", "is not one of"),
    ("rules/x.md", "---\nharness:\n  kind: rule\n---\n# R\n", "budget_tokens must be an integer"),
    ("rules/x.md", front("rule", 100, "  extra: 1\n") + "# R\n", "unknown keys"),
    ("rules/x.md", "---\nharness:\n  kind: rule\n  budget_tokens: 100\n  paths: tests\n---\n# R\n", "paths must be a list"),
    ("rules/x.md", "---\ntitle: t\nharness:\n  kind: rule\n  budget_tokens: 100\n---\n# R\n", "meaningful only on a skill"),
    ("rules/x.md", "---\nharness: [1, 2\n---\n# R\n", "not YAML"),
    ("rules/x.md", "---\nharness: {kind: rule,\nbudget_tokens: 100}\n---\n# R\n", "could not be removed textually"),
    ("rules/x.md", "---\nharness: &h\n  kind: rule\n  budget_tokens: 100\nother: *h\n---\n# R\n", "could not be removed textually"),
    ("skills/x.md", front("skill") + "# S\n", "needs the Agent Skills front matter"),
    ("rules/deep/x.md", front("rule") + "# R\n", "one level deep only"),
])
def test_source_faults_are_problems_never_exceptions(tmp_path: pathlib.Path, rel: str, text: str, fragment: str) -> None:
    root = tmp_path / "dev"
    write(root / rel, text)
    pages, problems = discover(root)
    assert pages == [] and len(problems) == 1 and fragment in problems[0]


# ------------------------------------------------------------ claude code

def test_claude_code_targets_every_kind_under_its_roots(tmp_path: pathlib.Path) -> None:
    cc = ClaudeCode(repo_root=tmp_path)
    targets = {kind: cc.target(page(kind, "n")) for kind in Kind}
    assert targets[Kind.RULE] == tmp_path / ".claude" / "rules" / "n.md"
    assert targets[Kind.SKILL] == tmp_path / ".claude" / "skills" / "n" / "SKILL.md"
    assert targets[Kind.AGENT] == tmp_path / ".claude" / "agents" / "n" / "AGENT.md"
    assert targets[Kind.INDEX] == tmp_path / ".claude" / "n.md"
    assert targets[Kind.ONBOARDING] == tmp_path / "CLAUDE.md"
    assert all(any(t == r or t.is_relative_to(r) for r in cc.roots) for t in targets.values())


def test_claude_code_render_shapes() -> None:
    cc = ClaudeCode()
    rule = cc.render(page(Kind.RULE, paths=("tests/**",)), "\n# R\n")
    assert isinstance(rule, WholeFile) and rule.text.startswith('---\npaths:\n  - "tests/**"\n---\n\n<!-- GENERATED by tools/harness')
    skill = cc.render(page(Kind.SKILL, front_matter="---\nname: s\n---\n"), f"# S\n{source.ERROR_INDEX_MARK}\n")
    assert isinstance(skill, WholeFile) and skill.text.startswith("---\nname: s\n---\n<!-- GENERATED") and ERROR_INDEX_LINE in skill.text
    assert cc.render(page(Kind.AGENT), "\n# a — role block\n\n") == Block("definition", "# a — role block")
    assert cc.render(page(Kind.ONBOARDING), "# O\n") == Block("on-boarding block", "# O")
    assert isinstance(cc.render(page(Kind.INDEX), "# L\n"), WholeFile)


def test_claude_code_always_on() -> None:
    cc = ClaudeCode()
    assert cc.always_on(page(Kind.RULE)) and cc.always_on(page(Kind.ONBOARDING))
    assert not cc.always_on(page(Kind.RULE, paths=("tests/**",)))
    assert not any(cc.always_on(page(k)) for k in (Kind.SKILL, Kind.AGENT, Kind.INDEX))


def test_the_real_vv_principles_skill_carries_the_injection_line() -> None:
    assert ERROR_INDEX_LINE in (REPO_ROOT / ".claude" / "skills" / "vv-principles" / "SKILL.md").read_text(encoding="utf-8")


# -------------------------------------------------------------- pipeline

class Stub:
    """A harness for one arm of the pipeline: rules only, every page to `out`."""
    name = "stub"
    kinds = frozenset({Kind.RULE})

    def __init__(self, root: pathlib.Path, out: pathlib.Path | None = None) -> None:
        self.roots = (root,)
        self.out = out or root / "one.md"

    def target(self, page: Page) -> pathlib.Path:
        return self.out

    def render(self, page: Page, body: str) -> WholeFile:
        return WholeFile(body)

    def always_on(self, page: Page) -> bool:
        return False

    def installed_always_on(self) -> tuple[pathlib.Path, ...]:
        return ()


def test_installed_rules_come_from_a_manifest_and_only_the_unscoped_ones_are_always_on(tmp_path: pathlib.Path) -> None:
    """X1: an installed rule without ``paths`` moves the always-on count; one with ``paths``,
    one the manifest never recorded, and a recorded file that is gone do not."""
    cc = ClaudeCode(repo_root=tmp_path)
    assert cc.installed() == {} and cc.installed_always_on() == (), "no manifest: nothing is installed"
    rules = tmp_path / ".claude" / "rules"
    write(rules / "routing.md", "# a routing rule\n")
    write(rules / "scoped.md", "---\npaths:\n  - \"tests/**\"\n---\n\n# scoped\n")
    write(rules / "hand.md", "# hand-written, recorded by nobody\n")
    write(tmp_path / ".claude" / "skills" / "s" / "SKILL.md", "# a skill\n")
    write(tmp_path / ".claude" / "nexus-install-manifest.json", json.dumps({"files": {
        "rules/routing.md": {"version": "0.18.0"}, "rules/scoped.md": {"version": "0.18.0"},
        "rules/gone.md": {"version": "0.18.0"}, "skills/s/SKILL.md": {"version": "0.18.0"}}}))
    assert cc.installed()[rules / "routing.md"] == "nexus 0.18.0"
    assert set(cc.installed()) == {rules / "routing.md", rules / "scoped.md", rules / "gone.md", tmp_path / ".claude" / "skills" / "s" / "SKILL.md"}
    assert cc.installed_always_on() == (rules / "routing.md",)


def test_the_real_tree_installs_exactly_the_nexus_routing_rule() -> None:
    cc = ClaudeCode()
    assert cc.installed_always_on() == (REPO_ROOT / ".claude" / "rules" / "nexus-tools.md",)
    assert cc.installed()[REPO_ROOT / ".claude" / "rules" / "nexus-tools.md"].startswith("nexus ")


def test_two_pages_on_one_target_is_a_problem_and_unrealised_kinds_are_skipped(tmp_path: pathlib.Path) -> None:
    root = tmp_path / "dev"
    a = write(root / "rules" / "a.md", front("rule") + "# A\n")
    b = write(root / "rules" / "b.md", front("rule") + "# B\n")
    write(root / "skills" / "s.md", SKILL_FRONT)
    pages, _ = discover(root)
    assert {p.kind for p in pages} == {Kind.RULE, Kind.SKILL}
    outputs, problems = pipeline.generate(Stub(tmp_path / "out"), pages)
    assert len(outputs) == 1 and len(problems) == 1 and "produced by both" in problems[0]
    assert source.rel(a) in problems[0] and source.rel(b) in problems[0]
    assert all(o.page.kind is Kind.RULE for o in outputs.values()), "the skill was never offered to a rules-only harness"


def test_a_target_outside_the_roots_is_a_problem(tmp_path: pathlib.Path) -> None:
    write(tmp_path / "dev" / "rules" / "a.md", front("rule") + "# A\n")
    pages, _ = discover(tmp_path / "dev")
    outputs, problems = pipeline.generate(Stub(tmp_path / "out", out=tmp_path / "elsewhere" / "a.md"), pages)
    assert outputs == {} and len(problems) == 1 and "outside the roots" in problems[0]


def test_generate_splices_a_block_and_surfaces_its_faults(tmp_path: pathlib.Path) -> None:
    cc = ClaudeCode(repo_root=tmp_path)
    root = tmp_path / "docs" / "development"
    write(root / "agents" / "a.md", front("agent") + "# a — role block\n\nrole\n")
    pages, _ = discover(root)
    dst = cc.target(pages[0])
    write(dst, "---\nname: a\ntools: [Read]\n---\n\n# Hand-written body\n")
    outputs, problems = pipeline.generate(cc, pages)
    assert problems == []
    text = outputs[dst].text
    assert text.startswith("---\nname: a\ntools: [Read]\n---\n\n<!-- BEGIN GENERATED definition") and text.endswith("# Hand-written body\n")
    assert outputs[dst].budgeted == block_text(Block("definition", "# a — role block\n\nrole"), pages[0].rel)
    write(dst, text.replace("<!-- END GENERATED definition -->\n", ""))
    _, problems = pipeline.generate(cc, pages)
    assert len(problems) == 1 and "BEGIN marker without END" in problems[0]


def test_drift_and_orphans(tmp_path: pathlib.Path) -> None:
    cc = ClaudeCode(repo_root=tmp_path)
    root = tmp_path / "docs" / "development"
    write(root / "rules" / "r.md", front("rule") + "# R\n")
    pages, _ = discover(root)
    outputs, problems = pipeline.generate(cc, pages)
    target = cc.target(pages[0])
    assert problems == [] and pipeline.drift(outputs) == [target]
    write(target, outputs[target].text)
    assert pipeline.drift(outputs) == []
    orphaned = [
        write(tmp_path / ".claude" / "rules" / "new.md", "<!-- GENERATED by tools/harness from x.md — DO NOT EDIT -->\n"),
        write(tmp_path / ".claude" / "agents" / "gone" / "AGENT.md", "---\nname: gone\n---\n<!-- BEGIN GENERATED definition — source: x; -->\n<!-- END GENERATED definition -->\n"),
        write(tmp_path / "CLAUDE.md", "<!-- BEGIN GENERATED on-boarding block — source: x; -->\n<!-- END GENERATED on-boarding block -->\n"),
    ]
    write(tmp_path / ".claude" / "skills" / "v" / "error_index.md", "<!-- GENERATED by tools/verification/generate_error_index.py — DO NOT EDIT. -->\n")
    write(tmp_path / ".claude" / "rules" / "hand.md", "# a hand-maintained rule\n")
    found = pipeline.orphans(cc, outputs)
    assert {f.split(": ")[0] for f in found} == {source.rel(p) for p in orphaned}


# ------------------------------------------------------------ neutrality

IMPLEMENTATIONS = {p.stem for p in (PKG / "targets").glob("*.py")} - {"base", "__init__"}
WIRING = {PKG / "__init__.py", PKG / "__main__.py", PKG / "targets" / "__init__.py"}  # may name and import harnesses
NEUTRAL = sorted(p for p in PKG.rglob("*.py") if p not in WIRING and p.stem not in IMPLEMENTATIONS)
HARNESS_WORDS = ("claude", "cursor", "codex")


def imports_of(path: pathlib.Path) -> set[str]:
    """Every imported dotted name: a module, and each ``from m import x`` as ``m.x``."""
    names: set[str] = set()
    for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
        if isinstance(node, ast.ImportFrom):
            names.add(node.module or "")
            names.update(f"{node.module or ''}.{alias.name}" for alias in node.names)
        elif isinstance(node, ast.Import):
            names.update(alias.name for alias in node.names)
    return names


def imports_an_implementation(path: pathlib.Path) -> bool:
    parts = {part for name in imports_of(path) for part in name.split(".")}
    return bool(parts & IMPLEMENTATIONS) or ("targets" in parts and "base" not in parts)


def names_a_harness(path: pathlib.Path) -> bool:
    text = path.read_text(encoding="utf-8").lower()
    return any(word in text for word in HARNESS_WORDS)


def test_neutral_modules_name_no_harness_and_import_no_implementation() -> None:
    assert len(NEUTRAL) >= 6 and IMPLEMENTATIONS == {"claude_code"}, (NEUTRAL, IMPLEMENTATIONS)
    for path in NEUTRAL:
        assert not names_a_harness(path), f"{path.name} names a harness"
        assert not imports_an_implementation(path), f"{path.name} imports an implementation"
    # positive controls: each detector fires where it must
    assert names_a_harness(PKG / "targets" / "claude_code.py")
    assert imports_an_implementation(PKG / "__main__.py"), "`from .targets import HARNESSES` is an implementation import"
    assert imports_an_implementation(PKG / "targets" / "__init__.py"), "`from .claude_code import ClaudeCode` is one"
    assert not imports_an_implementation(PKG / "pipeline.py"), "`from .targets.base import Harness` is the one allowed edge"


# ------------------------------------------------------------ brief rules

def test_brief_is_read_on_a_rule_and_refused_elsewhere(tmp_path: pathlib.Path) -> None:
    from tools.harness import brief as brief_mod

    root = tmp_path / "docs"
    write(root / "rules" / "b.md", front("rule", extra="  brief: >-\n    second sentence,\n    folded.\n") + "# B\n")
    write(root / "rules" / "a.md", front("rule", extra="  brief: one sentence.\n") + "# A\n")
    write(root / "rules" / "c.md", front("rule") + "# C\n")
    pages, problems = discover(root)
    assert problems == []
    assert {p.name: p.brief for p in pages} == {"a": "one sentence.", "b": "second sentence, folded.", "c": None}
    block = brief_mod.assemble(pages)
    assert block.label == brief_mod.LABEL
    assert block.text == "- `a`: one sentence.\n- `b`: second sentence, folded."  # page-name order; a rule without a brief contributes nothing
    write(root / "skills" / "s.md", SKILL_FRONT.replace("  budget_tokens: 400\n", "  budget_tokens: 400\n  brief: nope\n"))
    _, problems = discover(root)
    assert any("harness.brief is meaningful only on a rule" in p for p in problems)
    write(root / "rules" / "d.md", front("rule", extra="  brief: ''\n") + "# D\n")
    _, problems = discover(root)
    assert any("harness.brief must be a non-empty string" in p for p in problems)


def test_brief_block_is_spliced_and_a_hand_edit_inside_it_is_drift(tmp_path: pathlib.Path) -> None:
    from tools.harness import brief as brief_mod

    begin, end = markers(brief_mod.LABEL, brief_mod.SOURCE)
    page_path = write(tmp_path / "workflows.md", f"# W\n\ntext\n\n{begin}\nstale\n{end}\n\ntail\n")
    pages = [page(Kind.RULE, "r", path=tmp_path / "rules" / "r.md")]
    pages = [Page(p.kind, p.name, p.path, p.rel, p.front_matter, p.body, p.budget_tokens, p.paths, "the sentence.") for p in pages]
    text, problem = brief_mod.render(pages, page_path)
    assert problem is None
    assert f"{begin}\n- `r`: the sentence.\n{end}\n\ntail\n" in text and "stale" not in text and text.startswith("# W\n\ntext\n")
    page_path.write_text(text, encoding="utf-8")
    again, _ = brief_mod.render(pages, page_path)
    assert again == text  # idempotent: no drift once written
    page_path.write_text(text.replace("the sentence.", "an edit by hand."), encoding="utf-8")
    fixed, _ = brief_mod.render(pages, page_path)
    assert fixed == text  # the hand edit is drift: the generator restores the source's text


def test_the_real_workflows_page_carries_the_brief_block_of_every_rule_that_declares_one() -> None:
    from tools.harness import brief as brief_mod

    pages, problems = discover()
    assert problems == []
    declared = sorted(p.name for p in pages if p.kind is Kind.RULE and p.brief)
    assert declared, "no rule declares a brief: the block would be empty and the template would carry the copy by hand again"
    text = brief_mod.BRIEF_PAGE.read_text(encoding="utf-8")
    begin, end = markers(brief_mod.LABEL, brief_mod.SOURCE)
    inside = text.split(begin, 1)[1].split(end, 1)[0]
    assert [line.split("`")[1] for line in inside.strip().splitlines()] == declared


# ---------------------------------------------------------- citations by ID

def test_every_registry_reads_a_known_member_of_the_real_tree() -> None:
    from tools.harness import ids

    reg = ids.registries()
    assert "X2" in reg.x and "4" in reg.cardinal and "7" in reg.pattern and "8" in reg.mode
    assert "17" in reg.anti["vv-principles"] and "20" in reg.anti["coding-elegance"]
    assert "ERR-026" in reg.err and "L28" in reg.lesson and "B.4" in reg.item and "D.18" in reg.item and "E.18" not in reg.item
    assert "VALIDATE-THE-FILTER" in reg.tag and "A-LIST-IS-N-CENSUSES" in reg.tag


def test_citations_are_found_outside_code_and_qualified_in_their_paragraph() -> None:
    from tools.harness import ids

    text = ("X2 and Cardinal Rules 1, 4 and 5; Pattern 7 and Patterns 2 ∩ 4; mode 8; ERR-026 and `catches(\"ERR-999\")`;\n"
            "L28 but L1 and `L999`; B.4 and D.18; VALIDATE-THE-FILTER and A-LIST-IS-N-CENSUSES but ONE-HYPHEN.\n\n"
            "`vv-principles` #17, #11–#14 and `coding-elegance` #20.\n\n#3 alone.\n\n```\nX9 ERR-999 in a fence\n```\n")
    found = {(c.kind, c.id) for c in ids.citations(text, "some-page")}
    assert found == {("x", "X2"), ("cardinal", "1"), ("cardinal", "4"), ("cardinal", "5"), ("pattern", "7"),
                     ("pattern", "2"), ("pattern", "4"), ("mode", "8"), ("err", "ERR-026"), ("lesson", "L28"),
                     ("item", "B.4"), ("item", "D.18"), ("tag", "VALIDATE-THE-FILTER"), ("tag", "A-LIST-IS-N-CENSUSES"),
                     ("anti:vv-principles", "17"), ("anti:vv-principles", "11"), ("anti:vv-principles", "14"),
                     ("anti:coding-elegance", "20"), ("anti:None", "3")}
    own = {(c.kind, c.id) for c in ids.citations("#3 alone.", "vv-principles")}
    assert own == {("anti:vv-principles", "3")}


def test_a_dangling_citation_of_every_kind_is_a_problem(tmp_path: pathlib.Path) -> None:
    from tools.harness import ids

    seeded = write(tmp_path / "seeded.md", "X5. Cardinal Rule 6. Pattern 9. mode 13. ERR-999. L999. A.99. "
                                           "NOT-A-TAG-AT-ALL. `vv-principles` #99. `coding-elegance` #99.\n\n#3 alone.\n")
    control = write(tmp_path / "control.md", "X2, Cardinal Rule 4, Pattern 7, mode 8, ERR-026, L28, B.4, VALIDATE-THE-FILTER, `vv-principles` #17.\n")
    resolved, problems = ids.check(pages=[seeded, control])
    assert resolved == 9, problems
    dangling = sorted(p.split(": ", 1)[1].split(",")[0] for p in problems)
    assert dangling == sorted(["cites x X5", "cites cardinal 6", "cites pattern 9", "cites mode 13", "cites err ERR-999",
                               "cites lesson L999", "cites item A.99", "cites tag NOT-A-TAG-AT-ALL",
                               "cites anti 99", "cites anti 99", "cites #3 with no page named before it in the paragraph (an anti-pattern needs `vv-principles` or `coding-elegance` beside it; an issue is never a bare #N)"])


def test_the_real_tree_has_no_dangling_citation() -> None:
    from tools.harness import ids

    resolved, problems = ids.check()
    assert problems == [] and resolved > 200, (resolved, problems)


# ------------------------------------------------------------ citations by ID

def test_ids_check_reads_the_evidence_pages_with_the_three_decoder_rules(tmp_path: pathlib.Path) -> None:
    """X1 for the evidence-page check (2026-09-21): a dangling audit item reddens; a `task #N`,
    a `nexus#N`, a `mode-0` and a bare `#N` on a page whose subject is the anti-pattern list
    do not; the same bare `#N` on any other page does until a page is named before it."""
    from tools.harness import ids

    reg = ids.registries()
    item, anti = sorted(reg.item)[0], sorted(reg.anti["vv-principles"], key=int)[0]
    assert any(p.parent.name == "evidence" for p in ids.pages_under_check()), "the evidence pages are under check"
    on_anti_page = write(tmp_path / "vv-anti-patterns.md",
                         f"A case.\n\nAnti-pattern #{anti} fired; task #51 and nexus#85 are issues, so is Task #7; "
                         f"a mode-0 normalisation mismatch (ERR-030); item {item} applies; item G.99 does not.\n")
    elsewhere = write(tmp_path / "some-case.md", f"Bare #{anti} here, `vv-principles` #{anti} there.\n")
    resolved, problems = ids.check(pages=[on_anti_page, elsewhere])
    assert resolved == 4, problems                      # the anti-pattern, ERR-030, the live item, the qualified #N
    assert [p.split(":")[0].rsplit("/", 1)[-1] for p in problems] == ["vv-anti-patterns.md", "some-case.md"]
    assert "G.99" in problems[0] and "with no page named" in problems[1]
    assert not any(s in " ".join(problems) for s in ("#51", "#85", "#7", "mode 0", f"item {item}"))


# ------------------------------------------------------- code-span paths

def small_tree() -> "paths.Tree":
    files = frozenset({"orpheus/sn/solver.py", "orpheus/sn/__init__.py", "orpheus/__init__.py",
                       "tools/harness/ids.py", "tools/harness/__init__.py", "tools/__init__.py",
                       ".hidden/skills/v/scripts/case.md", "docs/theory/verification/error_catalog.rst"})
    dirs = frozenset("/".join(f.split("/")[:i]) for f in files for i in range(1, f.count("/") + 1))
    tails = frozenset("/".join(p.split("/")[i:]) for p in files | dirs for i in range(p.count("/") + 1))
    return paths.Tree(files, dirs, frozenset(f.rsplit("/", 1)[-1] for f in files),
                      frozenset(d for d in dirs if "/" not in d), tails)


@pytest.mark.parametrize("span, expected", [
    ("orpheus/sn/solver.py", True), ("orpheus/sn/solver.py::test_x", True), ("orpheus/sn/solver.py:120", True),
    ("sn/solver.py", True), ("scripts/case.md", True), ("orpheus/sn/", True),
    ("orpheus/sn/missing.py", False), ("docs/theory/discrete_ordinates.rst", False),
    ("tools.harness.ids", True), ("tools.harness.ids.check", True), ("tools._gone.check", False),
    ("error_catalog.rst", True), ("error_catalog.md", False),
    ("orpheus/sn/diag_<NN>.py", None), ("python -O -m pytest", None), ("np.testing.assert_allclose", None),
    ("SNMesh", None), ("a/b", None), ("tests/*.py", None),
])
def test_a_code_span_is_read_as_a_file_only_when_it_spells_one(span: str, expected: bool | None) -> None:
    assert small_tree().resolves(span) is expected


def test_a_dead_span_is_a_problem_at_its_file_line_and_fences_are_not_read(tmp_path: pathlib.Path) -> None:
    body = "text `orpheus/sn/solver.py`\n\n```\n`orpheus/sn/fenced_gone.py`\n```\n\nand `orpheus/sn/gone.py`\n"
    src = write(tmp_path / "p.md", front("rule") + body)
    resolved, problems = paths.check([page(Kind.RULE, body=body, path=src)], small_tree())
    assert resolved == 1 and len(problems) == 1, problems
    assert "gone.py" in problems[0] and "fenced" not in problems[0]
    assert f"p.md:{(front('rule') + body).count(chr(10))}:" in problems[0], "the line is the file's, front matter counted"


def test_a_path_git_ignores_is_local_by_design_not_dead(tmp_path: pathlib.Path) -> None:
    assert paths.ignored(["scratch/literature/x.pdf", "orpheus/sn/solver.py"]) == {"scratch/literature/x.pdf"}
    body = "the folder `scratch/literature/x.pdf`\n"
    src = write(tmp_path / "p.md", body)
    resolved, problems = paths.check([page(Kind.RULE, body=body, path=src)], paths.Tree.tracked())
    assert problems == [] and resolved == 0


def test_a_path_under_a_directory_pattern_is_ignored_even_where_the_directory_is_absent(tmp_path: pathlib.Path) -> None:
    """CI's checkout has no literature folder: a `dir/` pattern must still read the bare directory
    as ignored (the first push of the check read it dead there and ignored here)."""
    import subprocess
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    write(tmp_path / ".gitignore", "lit/\n")
    assert not (tmp_path / "lit").exists()
    assert paths.ignored(["lit", "lit/", "lit/x.pdf", "kept.md"], root=tmp_path) == {"lit", "lit/x.pdf"}

