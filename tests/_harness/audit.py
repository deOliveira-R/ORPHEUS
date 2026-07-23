"""V&V audit CLI — dump the level × module × equation grid.

Usage::

    python -m tests._harness.audit              # text report
    python -m tests._harness.audit --json       # machine-readable
    python -m tests._harness.audit --untagged   # list only unmarked tests
    python -m tests._harness.audit --gaps       # equations with no coverage

The tool runs ``pytest --collect-only`` under the hood so
:data:`tests._harness.registry.TEST_REGISTRY` is populated, then queries
the registry. No test code is executed.

Exit codes:
    0  clean run
    1  --strict was passed and the report has gaps or untagged tests
    2  collection failed, or a vv-status sentinel violation (unknown
       status / dead label / cross-file sentinel / malformed line —
       checked BEFORE collection, so bad V&V metadata fails fast)
"""

from __future__ import annotations

import argparse
import contextlib
import json
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, NamedTuple

import pytest

from tests._harness import registry
from tests._harness.registry import TestMetadata

# The canonical scan-exemption marker: a column-0 RST comment that
# excludes a file from BOTH the theory label/sentinel scan and the
# all-docs phantom census. The generator emits this exact string into
# the matrix header; the regex below is its (whitespace-tolerant)
# parse — the ONE spelling both scan functions share, so the two gates
# can never drift apart on what counts as exempt.
VV_AUDIT_SKIP_MARKER = ".. vv-audit: skip-file"
_VV_AUDIT_SKIP_RE = re.compile(
    r"^\.\.\s+vv-audit:\s*skip-file\s*$", re.MULTILINE
)


def _run_collection() -> int:
    """Invoke pytest in collect-only mode so TEST_REGISTRY gets populated.

    Pytest's own stdout (the flat list of 497 nodeids) is discarded so
    the audit report is the only thing the user sees. If pytest fails,
    its stderr is still routed to the real stream.
    """
    devnull = open(os.devnull, "w")  # noqa: SIM115 — must stay open for duration
    try:
        with contextlib.redirect_stdout(devnull):
            return pytest.main(
                [
                    "--collect-only",
                    "-q",
                    "--no-header",
                    "--disable-warnings",
                    "-p",
                    "no:cacheprovider",
                ]
            )
    finally:
        devnull.close()


def _module_of(file_path: str) -> str:
    """Return the pytest test-file's display name for the module grid.

    The nested tests/ layout (issue #77) groups files by solver
    module, so multiple files share a basename (e.g. both
    ``tests/cp/test_properties.py`` and ``tests/sn/test_properties.py``
    have stem ``test_properties``). Collapsing them to the bare stem
    hides the per-module breakdown, so the grid uses the
    ``parent/stem`` form — e.g. ``cp/test_properties`` — when the
    file lives below a subfolder of ``tests/``. Files at the
    ``tests/`` root (``test_pending_ports``,
    ``test_vv_harness_audit``, ``test_convergence``) keep their
    bare stem since there is no parent folder to disambiguate.
    """
    p = Path(file_path)
    parent = p.parent.name
    if parent and parent != "tests":
        return f"{parent}/{p.stem}"
    return p.stem


def _group_by_module_level(
    items: list[TestMetadata],
) -> dict[str, Counter]:
    out: dict[str, Counter] = defaultdict(Counter)
    for m in items:
        out[_module_of(m.file)][m.level or "unmarked"] += 1
    return out


def _equation_coverage(items: list[TestMetadata]) -> dict[str, list[str]]:
    coverage: dict[str, list[str]] = defaultdict(list)
    for m in items:
        for eq in m.equations:
            coverage[eq].append(m.nodeid)
    return coverage


def _caught_tags(items: list[TestMetadata]) -> dict[str, list[str]]:
    caught: dict[str, list[str]] = defaultdict(list)
    for m in items:
        for tag in m.catches:
            caught[tag].append(m.nodeid)
    return caught


def _phantom_verifies(
    coverage: dict[str, list[str]], doc_labels: set[str]
) -> dict[str, list[str]]:
    """Verifies-targets declared by tests with NO matching doc ``:label:``.

    The inverse of the orphan-equation gate (issue #224): a theory-page
    label rename/removal that is not migrated into the
    ``@pytest.mark.verifies`` strings silently drops the test from the
    V&V matrix — the Nexus graph build skips the edge with only a log
    line. This gate makes that drift loud at audit time.

    ``doc_labels`` is scanned from ALL of ``docs/`` (not just the
    theory tree) — the gate checks that a verifies-target *resolves*,
    wherever the label lives.
    """
    return {
        eq: tests for eq, tests in coverage.items() if eq not in doc_labels
    }


# ---------------------------------------------------------------------------
# Reporters
# ---------------------------------------------------------------------------


def _render_text(
    items: list[TestMetadata],
    *,
    theory_labels: set[str],
    documented_labels: set[str],
    all_doc_labels: set[str],
    err_tags: set[str],
    skipped_files: list[str],
) -> str:
    lines: list[str] = []
    total = len(items)
    level_totals = Counter(m.level or "unmarked" for m in items)
    source_totals = Counter(m.level_source for m in items)

    lines.append("=" * 72)
    lines.append("ORPHEUS V&V Test Audit")
    lines.append("=" * 72)
    lines.append(f"Total tests collected: {total}")
    lines.append("")
    lines.append("By V&V level:")
    # L0..L3 are the physics-verification ladder. foundation is the
    # orthogonal software-invariant bucket; it is reported here for
    # visibility but is not part of the L0..L3 progression. "unmarked"
    # is a gap that the --strict gate surfaces.
    for lvl in ("L0", "L1", "L2", "L3", "foundation", "unmarked"):
        count = level_totals.get(lvl, 0)
        pct = 100 * count / total if total else 0
        lines.append(f"  {lvl:11} {count:5}   ({pct:4.1f}%)")
    lines.append("")
    lines.append("By tagging source:")
    for src in (
        "explicit",
        "class-name",
        "func-name",
        "case",
        "unmarked",
    ):
        count = source_totals.get(src, 0)
        lines.append(f"  {src:12} {count:5}")
    lines.append("")

    # Module × level grid. ``FD`` column counts foundation-marker tests
    # (software invariants, orthogonal to the physics ladder).
    grid = _group_by_module_level(items)
    lines.append("Module × level grid:")
    header = (
        f"  {'module':<36} "
        f"{'L0':>4} {'L1':>4} {'L2':>4} {'L3':>4} {'FD':>4} {'??':>4}"
    )
    lines.append(header)
    lines.append("  " + "-" * (len(header) - 2))
    for module in sorted(grid):
        row = grid[module]
        lines.append(
            f"  {module:<36} "
            f"{row.get('L0', 0):>4} "
            f"{row.get('L1', 0):>4} "
            f"{row.get('L2', 0):>4} "
            f"{row.get('L3', 0):>4} "
            f"{row.get('foundation', 0):>4} "
            f"{row.get('unmarked', 0):>4}"
        )
    lines.append("")

    # Equation coverage
    coverage = _equation_coverage(items)
    lines.append("Equation coverage:")
    if not coverage:
        lines.append("  (no tests declare @pytest.mark.verifies yet)")
    else:
        for eq in sorted(coverage):
            lines.append(f"  {eq:40} {len(coverage[eq]):>3} test(s)")
    lines.append("")

    # Phantom verifies-targets (declared by tests, no matching doc
    # label anywhere under docs/) — the issue-#224 drift class.
    phantoms = _phantom_verifies(coverage, all_doc_labels)
    if phantoms:
        lines.append(
            f"PHANTOM verifies targets ({len(phantoms)} label(s) declared "
            "by tests with NO matching :label: anywhere under docs/ — "
            "these tests silently drop out of the V&V matrix):"
        )
        for eq in sorted(phantoms):
            lines.append(f"  {eq:40} {len(phantoms[eq]):>3} test(s)")
        lines.append("")

    # Orphan equations (declared in theory pages, never referenced by
    # any test) — excluding labels explicitly marked `.. vv-status: X
    # documented` as definitional or not-yet-implemented.
    testable_labels = theory_labels - documented_labels
    orphans = sorted(testable_labels - coverage.keys())
    if theory_labels:
        lines.append(
            f"Orphan equations ({len(orphans)} of {len(testable_labels)} "
            "testable theory labels have zero test coverage; "
            f"{len(documented_labels)} labels are :vv-status: documented "
            "and excluded from the orphan gate):"
        )
        for eq in orphans:
            lines.append(f"  {eq}")
        lines.append("")

    # Scan-exempt files — always visible, so the skip-file marker can
    # never become a silent exclusion channel.
    if skipped_files:
        lines.append(
            f"Theory files excluded from the label/sentinel scan by an "
            f"explicit '.. vv-audit: skip-file' marker "
            f"({len(skipped_files)}):"
        )
        for f in skipped_files:
            lines.append(f"  {f}")
        lines.append("")

    # ERR catalog cross-check
    caught = _caught_tags(items)
    if err_tags:
        missing = sorted(err_tags - caught.keys())
        lines.append(
            f"error_catalog.md ERR coverage "
            f"({len(err_tags) - len(missing)}/{len(err_tags)} entries have a "
            "catching test):"
        )
        for err in missing:
            lines.append(f"  MISSING {err}")
        lines.append("")

    return "\n".join(lines)


def _render_json(
    items: list[TestMetadata],
    *,
    theory_labels: set[str],
    documented_labels: set[str],
    all_doc_labels: set[str],
    err_tags: set[str],
    skipped_files: list[str],
) -> str:
    coverage = _equation_coverage(items)
    caught = _caught_tags(items)
    testable_labels = theory_labels - documented_labels
    payload: dict[str, Any] = {
        "phantom_verifies": {
            eq: tests
            for eq, tests in sorted(
                _phantom_verifies(coverage, all_doc_labels).items()
            )
        },
        "total": len(items),
        "by_level": dict(Counter(m.level or "unmarked" for m in items)),
        "by_source": dict(Counter(m.level_source for m in items)),
        "grid": {
            module: dict(counts)
            for module, counts in _group_by_module_level(items).items()
        },
        "equation_coverage": {eq: tests for eq, tests in coverage.items()},
        "orphan_equations": sorted(testable_labels - coverage.keys()),
        "documented_equations": sorted(documented_labels),
        "err_coverage": {err: caught.get(err, []) for err in sorted(err_tags)},
        "untagged": [m.nodeid for m in items if m.level is None],
        "skipped_theory_files": list(skipped_files),
    }
    return json.dumps(payload, indent=2, sort_keys=True)


# ---------------------------------------------------------------------------
# External inputs: theory equation labels + ERR catalog
# ---------------------------------------------------------------------------


class TheoryScan(NamedTuple):
    """Result of :func:`_scan_theory_equations`.

    A named carrier rather than a bare tuple because two fields
    (``violations`` and ``skipped``) share the type ``list[str]`` with
    opposite severities — a positional swap would silently convert a
    build-breaking sentinel violation into a report-only skip line.
    Named access makes the swap unspellable at the consumer.
    """

    all_labels: set[str]
    documented: set[str]
    violations: list[str]
    skipped: list[str]


def _scan_theory_equations(theory_dir: Path) -> TheoryScan:
    """Scan theory RST pages for equation labels and documented-only markers.

    Returns
    -------
    TheoryScan
        ``all_labels`` is every ``.. math:: :label: foo`` found under
        ``theory_dir`` (recursive). ``documented`` is the subset
        of those labels that carry the V&V-harness-specific sentinel

            .. vv-status: <label> documented

        in the **same file** as the ``:label:`` (a sentinel is
        point-of-use metadata; the same-file rule is enforced). This is
        a plain RST comment — the ``.. `` prefix followed by text that
        is not a known directive is silently stripped by Sphinx, so the
        sentinel has no effect on the rendered theory page. The audit
        tool uses it to exclude three kinds of label from the
        orphan-equation gate:

        1. **Pure definitional labels** — e.g. ``boltzmann``,
           ``transport-equation``, ``balance-general``. These name the
           governing equation or a mathematical identity that has no
           single "implementing function" to test against. They
           belong in the theory page for the narrative but cannot be
           paired with a verifying test.
        2. **Not-yet-implemented modules** — e.g. the TH / FB / RK
           equations whose Python ports do not exist yet. A real
           orphan (implemented but untested) is a V&V gap; a
           documented-but-not-implemented equation is a work-in-
           progress marker, not a gap.
        3. **Equations with a pending catching test** — when an
           author deliberately wants to defer writing a test and
           surface it as "acceptable gap" rather than "bug in audit",
           marking ``documented`` is the escape hatch. This should
           be rare and paired with a GitHub issue.

        The orphan gate (and ``--strict``) only fires for labels that
        are *not* in ``documented_labels``.

        ``documented`` is the ONLY status. ``tested`` / ``verified``
        are **derived** facts (from ``@pytest.mark.verifies``) — a
        hand-written coverage claim would be a second source of truth
        that can silently lie, so the sentinel vocabulary deliberately
        excludes them.

        ``violations`` is a list of ``path:lineno: message`` strings,
        one per malformed sentinel: an unknown status word, a sentinel
        whose label has no ``:label:`` in the same file (dead or
        misplaced), or a line that does not parse as
        ``.. vv-status: <label> <status>``. Violations are a hard
        audit error (exit 2) — invalid V&V metadata is the same
        failure class as a collection error, never a silent drop.
        (Sentinel-syntax is defined by ``sentinel_re`` below; the
        exemption marker by the module-level ``_VV_AUDIT_SKIP_RE``.)

        ``skipped`` (paths relative to ``theory_dir``) are files
        excluded from the scan by an explicit column-0
        ``.. vv-audit: skip-file`` comment. The scanner is line-based
        and cannot tell a literal-block *teaching example* of the
        ``:label:`` / ``.. vv-status:`` syntax from the real thing, so
        the two pages whose label mentions are not declarations — the
        harness architecture page (verbatim syntax examples) and the
        generated matrix page (prose about the census) — opt out with
        the marker. The exclusion is reported in every output mode,
        never silent; a real theory page must never carry it (that
        would hide genuine equations from the orphan gate).
    """
    if not theory_dir.is_dir():
        return TheoryScan(set(), set(), [], [])

    sentinel_re = re.compile(r"^\.\.\s+vv-status:(.*)$")

    # Pass 1 — collect per-file labels and sentinel lines, so the
    # same-file rule can be checked and a misplaced sentinel can be
    # distinguished from a dead one in the violation message.
    per_file_labels: dict[Path, set[str]] = {}
    per_file_sentinels: dict[Path, list[tuple[int, str]]] = {}
    skipped: list[str] = []
    for rst in theory_dir.rglob("*.rst"):
        try:
            text = rst.read_text(encoding="utf-8")
        except OSError:
            continue
        if _VV_AUDIT_SKIP_RE.search(text):
            skipped.append(str(rst.relative_to(theory_dir)))
            continue
        labels_here: set[str] = set()
        sentinels_here: list[tuple[int, str]] = []
        for lineno, line in enumerate(text.splitlines(), start=1):
            stripped = line.strip()
            if stripped.startswith(":label:"):
                labels_here.add(stripped.split(":", 2)[2].strip())
                continue
            m = sentinel_re.match(stripped)
            if m:
                sentinels_here.append((lineno, m.group(1).strip()))
        per_file_labels[rst] = labels_here
        per_file_sentinels[rst] = sentinels_here

    all_labels = set().union(*per_file_labels.values()) if per_file_labels else set()

    # Pass 2 — validate each sentinel against the single-status schema
    # and the same-file rule; collect the documented set from the valid
    # ones only.
    documented: set[str] = set()
    violations: list[str] = []
    for rst, sentinels in per_file_sentinels.items():
        labels_here = per_file_labels[rst]
        for lineno, rest in sentinels:
            parts = rest.split()
            if len(parts) != 2:
                violations.append(
                    f"{rst}:{lineno}: malformed vv-status line "
                    f"(expected '.. vv-status: <label> documented')"
                )
                continue
            label, status = parts
            if status != "documented":
                violations.append(
                    f"{rst}:{lineno}: unknown vv-status {status!r} for "
                    f"{label!r} — 'documented' is the only status; "
                    f"tested/verified are derived from "
                    f"@pytest.mark.verifies"
                )
                continue
            if label not in labels_here:
                where = (
                    "it exists in another file — move the sentinel there"
                    if label in all_labels
                    else "no such equation :label: anywhere under the "
                    "theory tree — dead sentinel"
                )
                violations.append(
                    f"{rst}:{lineno}: vv-status names {label!r} but no "
                    f"such :label: exists in this file ({where})"
                )
                continue
            documented.add(label)

    return TheoryScan(all_labels, documented, sorted(violations), sorted(skipped))


def _scan_all_doc_labels(docs_dir: Path) -> set[str]:
    """Every ``:label:`` under ``docs/`` (excluding ``_build``).

    The phantom-verifies gate compares against the FULL doc label set,
    not just the theory tree — a verifies-target must *resolve*,
    wherever the label lives. ``_build`` is excluded so stale build
    artifacts cannot mask a genuinely-removed label, and files carrying
    the ``.. vv-audit: skip-file`` marker are excluded so a syntax
    *example* on a teaching page can never mask a genuine phantom.
    """
    if not docs_dir.is_dir():
        return set()

    labels: set[str] = set()
    for rst in docs_dir.rglob("*.rst"):
        if "_build" in rst.parts:
            continue
        try:
            text = rst.read_text(encoding="utf-8")
        except OSError:
            continue
        if _VV_AUDIT_SKIP_RE.search(text):
            continue
        for line in text.splitlines():
            stripped = line.strip()
            if stripped.startswith(":label:"):
                labels.add(stripped.split(":", 2)[2].strip())
    return labels


def _scan_err_catalog(catalog: Path) -> set[str]:
    """Extract ``ERR-NNN`` IDs from the L0 error catalog markdown."""
    if not catalog.is_file():
        return set()
    text = catalog.read_text(encoding="utf-8")
    return set(re.findall(r"\bERR-\d{3}\b", text))


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m tests._harness.audit",
        description="ORPHEUS V&V test audit",
    )
    parser.add_argument("--json", action="store_true", help="JSON output")
    parser.add_argument(
        "--untagged",
        action="store_true",
        help="list only tests with no V&V level",
    )
    parser.add_argument(
        "--gaps",
        action="store_true",
        help="list orphan equations + missing ERR catchers",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="exit 1 if untagged tests or orphan equations are present",
    )
    parser.add_argument(
        "--theory-dir",
        type=Path,
        default=Path("docs/theory"),
        help="Sphinx theory directory (for orphan-equation scan)",
    )
    parser.add_argument(
        "--err-catalog",
        type=Path,
        default=Path(".claude/skills/vv-principles/error_catalog.md"),
        help="L0 error catalog markdown (for ERR coverage)",
    )
    args = parser.parse_args(argv)

    # Sentinel violations are a hard error BEFORE collection: invalid
    # V&V metadata is the same failure class as a collection error.
    # The Sphinx matrix hook surfaces the message as a build warning
    # (fatal under ``-W``), so a bad sentinel breaks the docs build.
    scan = _scan_theory_equations(args.theory_dir)
    if scan.violations:
        print(
            f"vv-status sentinel violations ({len(scan.violations)}):",
            file=sys.stderr,
        )
        for v in scan.violations:
            print(f"  {v}", file=sys.stderr)
        return 2

    rc = _run_collection()
    if rc != 0 and rc != 5:  # 5 == pytest "no tests"
        print(f"pytest collection failed (exit {rc})", file=sys.stderr)
        return 2

    items = sorted(registry.TEST_REGISTRY.values(), key=lambda m: m.nodeid)
    testable_labels = scan.all_labels - scan.documented
    all_doc_labels = _scan_all_doc_labels(args.theory_dir.parent)
    err_tags = _scan_err_catalog(args.err_catalog)

    if args.untagged:
        for m in items:
            if m.level is None:
                print(m.nodeid)
    elif args.gaps:
        coverage = _equation_coverage(items)
        caught = _caught_tags(items)
        orphans = sorted(testable_labels - coverage.keys())
        phantoms = sorted(_phantom_verifies(coverage, all_doc_labels))
        missing_err = sorted(err_tags - caught.keys())
        if orphans:
            print("# Orphan equations (no verifying tests)")
            for eq in orphans:
                print(eq)
        if phantoms:
            if orphans:
                print()
            print("# Phantom verifies targets (no matching doc :label:)")
            for eq in phantoms:
                print(eq)
        if missing_err:
            if orphans or phantoms:
                print()
            print("# ERR entries with no catching test")
            for err in missing_err:
                print(err)
    elif args.json:
        print(
            _render_json(
                items,
                theory_labels=scan.all_labels,
                documented_labels=scan.documented,
                all_doc_labels=all_doc_labels,
                err_tags=err_tags,
                skipped_files=scan.skipped,
            )
        )
    else:
        print(
            _render_text(
                items,
                theory_labels=scan.all_labels,
                documented_labels=scan.documented,
                all_doc_labels=all_doc_labels,
                err_tags=err_tags,
                skipped_files=scan.skipped,
            )
        )

    if args.strict:
        untagged = sum(1 for m in items if m.level is None)
        coverage = _equation_coverage(items)
        orphans = testable_labels - coverage.keys()
        phantoms = _phantom_verifies(coverage, all_doc_labels)
        if untagged or orphans or phantoms:
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
