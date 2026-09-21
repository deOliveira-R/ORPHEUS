"""``python -m tools.harness`` writes every harness's view of the docs;
``--check`` verifies only and exits 1 on a problem or on drift. The reader of
``--check`` is ``tests/test_harness_generated.py``; write mode is a
``_GENERATORS`` row in ``docs/conf.py``, so every Sphinx build regenerates."""
from __future__ import annotations

import argparse
import sys

from . import brief, budget, ids
from .pipeline import drift, generate, orphans
from .source import DOCS_DEV, discover, rel
from .targets import HARNESSES


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="write every harness's view of docs/development/; --check verifies only")
    ap.add_argument("--check", action="store_true", help="verify only; exit 1 on drift or any problem")
    args = ap.parse_args(argv)
    pages, source_problems = discover(DOCS_DEV)
    err = sys.stdout if args.check else sys.stderr
    for problem in source_problems:
        print(f"PROBLEM: {problem}", file=err)
    rc = 1 if source_problems else 0  # a malformed source blocks writing: the view would be partial
    for harness in HARNESSES:
        outputs, problems = generate(harness, pages)
        problems += orphans(harness, outputs)
        drifted = drift(outputs)
        always = [o for o in outputs.values() if harness.always_on(o.page)]
        installed = harness.installed_always_on()
        generated_cost = sum(budget.tokens(o.text) for o in always)
        installed_cost = sum(budget.tokens(p.read_text(encoding="utf-8")) for p in installed)
        cost = (f"always-on ≈{generated_cost + installed_cost} tokens: generated ≈{generated_cost} over {len(always)} files"
                + (f", installed ≈{installed_cost} over {len(installed)} ({', '.join(p.name for p in installed)})" if installed else ""))
        for problem in problems:
            print(f"PROBLEM: {problem}", file=err)
        if args.check:
            for path in drifted:
                print(f"DRIFT: {rel(path)} differs from its source (run the generator)")
            print(f"harness {harness.name}: {len(outputs)} targets, {len(problems)} problems, {len(drifted)} drifted; {cost}")
            rc |= int(bool(problems or drifted))
        elif problems or source_problems:
            rc = 1
        else:
            for path in drifted:
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(outputs[path].text, encoding="utf-8")
            print(f"harness {harness.name}: wrote {len(outputs)} targets ({len(drifted)} changed); {cost}")
    # Plain-text citations by ID resolve against their registries, whatever the harness.
    resolved, id_problems = ids.check()
    for problem in id_problems:
        print(f"PROBLEM: {problem}", file=err)
    rc |= int(bool(id_problems))
    if args.check:
        print(f"citations by ID: {resolved} resolved, {len(id_problems)} dangling")
    # The brief-rules block is a source-to-source derivation, the same for every harness.
    text, problem = brief.render(pages)
    if problem:
        print(f"PROBLEM: {problem}", file=err)
        rc = 1
    elif text != brief.BRIEF_PAGE.read_text(encoding="utf-8"):
        if args.check:
            print(f"DRIFT: {rel(brief.BRIEF_PAGE)} differs from the rules' harness.brief (run the generator)")
            rc = 1
        elif not source_problems:
            brief.BRIEF_PAGE.write_text(text, encoding="utf-8")
            print(f"brief rules: wrote the block on {rel(brief.BRIEF_PAGE)}")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
