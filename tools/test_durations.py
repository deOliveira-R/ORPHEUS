r"""Per-test durations of the whole gate tree, measured sharded on CI (#405).

Two subcommands, used by ``.github/workflows/test-durations.yml``:

``shard``
    Lists every tracked test file under ``tests/gates/`` and deals the files
    round-robin into ``--shards`` groups, so every tree spreads over every
    shard.  Prints the GitHub Actions matrix as JSON.

``aggregate``
    Reads the JUnit XML files the shards wrote, and writes one table: every
    test case with its wall time and outcome, and the runner configuration
    the times were measured on.  A timing is comparable only between two
    tables whose runner stamps match (ruling R2 of
    ``.claude/plans/test_runtime_405.md``); the stamp travels with every
    table for that reason.  Also writes a Markdown summary of the slowest
    tests for the run's summary page.

Timing is never a gate (ruling R6 of ``.claude/plans/vv_suite_layout.md``):
this measures the cost of running the gates, so that slow verification can
be removed at its source.
"""
from __future__ import annotations

import argparse
import json
import os
import pathlib
import platform
import subprocess
import sys
import xml.etree.ElementTree as ET
from dataclasses import asdict, dataclass

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
GATES = "tests/gates"


@dataclass(frozen=True)
class CaseTiming:
    """One JUnit test case: its pytest node id, its wall time, and its outcome."""

    node_id: str
    seconds: float
    outcome: str  # passed | failed | error | skipped | timeout


def gate_files() -> list[str]:
    """Every tracked ``test_*.py`` under ``tests/gates/``, sorted."""
    listed = subprocess.run(
        ["git", "ls-files", "-z", "--", f"{GATES}/**/test_*.py", f"{GATES}/test_*.py"],
        cwd=REPO_ROOT, capture_output=True, text=True, check=True,
    ).stdout
    return sorted(f for f in listed.split("\0") if f)


def deal(files: list[str], shards: int) -> list[list[str]]:
    """Deal ``files`` round-robin into ``shards`` non-empty groups."""
    groups: list[list[str]] = [[] for _ in range(min(shards, len(files)))]
    for i, f in enumerate(files):
        groups[i % len(groups)].append(f)
    return groups


def _node_id(case: ET.Element, root: pathlib.Path) -> str:
    """The pytest node id of a JUnit ``testcase``: ``path::Class::name``.

    pytest's JUnit writer puts the dotted module path (and any class) in
    ``classname`` and the test, with its parametrize id, in ``name``.  The
    module is the longest dotted prefix that names an existing ``.py`` file
    under ``root``; the rest of ``classname`` is the class nesting.
    """
    classname, name = case.get("classname", ""), case.get("name", "")
    parts = classname.split(".")
    module_end = next(
        (i for i in range(len(parts) - 1, -1, -1) if (root / ("/".join(parts[: i + 1]) + ".py")).is_file()),
        None,
    )
    if module_end is None:
        raise ValueError(f"no test module under {root} matches the JUnit classname {classname!r}")
    return "::".join(["/".join(parts[: module_end + 1]) + ".py", *parts[module_end + 1:], name])


def _outcome(case: ET.Element) -> str:
    for child in case:
        if child.tag in ("failure", "error"):
            text = (child.get("message", "") + (child.text or "")).lower()
            return "timeout" if "timeout" in text else ("failed" if child.tag == "failure" else "error")
        if child.tag == "skipped":
            return "skipped"
    return "passed"


def read_junit(path: pathlib.Path, root: pathlib.Path = REPO_ROOT) -> list[CaseTiming]:
    """Every test case in one JUnit XML file, with node ids resolved against ``root``.

    A case's time is pytest's JUnit ``time``: setup, call and teardown
    together, so a module- or session-scoped fixture's cost lands on the
    first test that requests it.
    """
    tree = ET.parse(path).getroot()
    return [
        CaseTiming(_node_id(case, root), float(case.get("time", "0") or 0), _outcome(case))
        for case in tree.iter("testcase")
    ]


def runner_stamp() -> dict[str, str]:
    """The configuration a timing is only comparable within."""
    import numpy
    import scipy
    import sympy

    return {
        "runner_os": os.environ.get("ImageOS", platform.system()),
        "runner_image_version": os.environ.get("ImageVersion", "local"),
        "cpu_count": str(os.cpu_count()),
        "machine": platform.machine(),
        "python": platform.python_version(),
        "numpy": numpy.__version__,
        "scipy": scipy.__version__,
        "sympy": sympy.__version__,
        "commit": os.environ.get("GITHUB_SHA", "local"),
    }


def tree_of(node_id: str) -> str:
    """The first directory under ``tests/gates/`` of a node id's file, or ``(root)``.

    Only the path part is read (before the first ``::``): a parametrize id
    may itself contain ``/``.
    """
    path = pathlib.PurePosixPath(node_id.split("::", 1)[0])
    parts = path.relative_to(GATES).parts
    return parts[0] if len(parts) > 1 else "(root)"


def summary_markdown(cases: list[CaseTiming], top: int) -> str:
    """The slowest ``top`` cases, and the time per tree, as Markdown."""
    by_tree: dict[str, float] = {}
    for c in cases:
        by_tree[tree_of(c.node_id)] = by_tree.get(tree_of(c.node_id), 0.0) + c.seconds
    lines = [f"## Test durations: {len(cases)} cases, {sum(c.seconds for c in cases) / 3600:.2f} h serial", ""]
    lines += ["| tree | hours |", "|---|---|"]
    lines += [f"| {t} | {s / 3600:.2f} |" for t, s in sorted(by_tree.items(), key=lambda kv: -kv[1])]
    lines += ["", f"### The {top} slowest", "", "| seconds | outcome | test |", "|---|---|---|"]
    for c in sorted(cases, key=lambda c: -c.seconds)[:top]:
        lines.append(f"| {c.seconds:.1f} | {c.outcome} | `{c.node_id}` |")
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Per-test durations of the gate tree, sharded on CI (#405).")
    sub = parser.add_subparsers(dest="command", required=True)
    shard = sub.add_parser("shard", help="print the shard matrix as JSON")
    shard.add_argument("--shards", type=int, required=True)
    aggregate = sub.add_parser("aggregate", help="merge JUnit XML into the duration table")
    aggregate.add_argument("junit_dir", type=pathlib.Path)
    aggregate.add_argument("--out", type=pathlib.Path, required=True)
    aggregate.add_argument("--summary", type=pathlib.Path, required=True)
    aggregate.add_argument("--top", type=int, default=50)
    args = parser.parse_args(argv)

    if args.command == "shard":
        groups = deal(gate_files(), args.shards)
        print(json.dumps({"include": [{"shard": i, "files": " ".join(g)} for i, g in enumerate(groups)]}))
        return 0

    files = sorted(args.junit_dir.rglob("*.xml"))
    if not files:
        print(f"no JUnit XML under {args.junit_dir}", file=sys.stderr)
        return 2
    cases = [case for f in files for case in read_junit(f)]
    table = {"runner": runner_stamp(), "junit_files": len(files), "cases": [asdict(c) for c in cases]}
    args.out.write_text(json.dumps(table, indent=1))
    args.summary.write_text(summary_markdown(cases, args.top))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
