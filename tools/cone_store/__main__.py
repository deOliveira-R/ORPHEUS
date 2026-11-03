"""``python -m tools.cone_store {capture,select,status,place} --repo R --store S [--package P] [--memo-root M] [ids]``.

``place`` takes ``--layer`` (the test-id prefix, ``tests/gates/numerics/``), ``--code-root`` (``orpheus/numerics``),
``--order FILE`` (the file holding the concept order: its ``TIERS`` and its import graph ``import_kinds``, read by
path) and ``--out`` (the JSON report; printed when absent). It measures the layer's declared ``rests_on`` edges as
pytest resolves them, collected in a fresh interpreter by the one resolver, ``tests/_harness/rests_on_manifest.py``."""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any


def main(argv: list[str] | None = None) -> int:
    started = time.time_ns()
    parser = argparse.ArgumentParser(prog="python -m tools.cone_store", description=__doc__)
    parser.add_argument("command", choices=("capture", "select", "status", "place"))
    parser.add_argument("ids", nargs="*", help="pytest node ids or paths (default: the project's testpaths)")
    parser.add_argument("--repo", default=".", help="the repository root (default: the working directory)")
    parser.add_argument("--store", default=None, help="the SQLite store (default: <repo>/.cache/cone_store/cones.sqlite)")
    parser.add_argument("--package", default="orpheus", help="the first-party package under the repository")
    parser.add_argument("--memo-root", default=None, help="the traced memo's store (default: the memo's own)")
    parser.add_argument("--route", choices=("auto", "fresh"), default="auto",
                        help="capture: auto forks where the memo's guard allows; fresh runs every test in its own interpreter")
    parser.add_argument("--why", action="store_true", help="select: print the reasons beside each id")
    parser.add_argument("--layer", default="tests/gates/numerics/", help="place: the test-id prefix of the layer")
    parser.add_argument("--code-root", default="orpheus/numerics", help="place: the layer's code, relative to the repository")
    parser.add_argument("--order", default="tests/gates/numerics/test_concept_order.py",
                        help="place: the file holding the concept order (TIERS) and its import graph (import_kinds)")
    parser.add_argument("--threshold", type=float, default=0.95, help="place: the containment a proposed support needs")
    parser.add_argument("--out", default=None, help="place: where to write the JSON report (default: print it)")
    options = parser.parse_intermixed_args(argv)
    repo = Path(options.repo).resolve()
    store = Path(options.store) if options.store else repo / ".cache" / "cone_store" / "cones.sqlite"
    import tools.cone_store as cone_store

    if options.command == "capture":
        print(cone_store.capture(repo, store, options.package, options.memo_root, options.ids, options.route, started))
    elif options.command == "select":
        selection = cone_store.select(repo, store, options.package, options.memo_root, options.ids)
        for test_id in selection.ids:
            print(f"{test_id}  # {'; '.join(selection.reasons[test_id])}" if options.why else test_id)
    elif options.command == "place":
        order = _order_module(repo, options.order)
        report = cone_store.place(store, options.layer, options.code_root, [list(tier) for tier in order.TIERS],
                                  order.import_kinds(repo / options.code_root), threshold=options.threshold,
                                  declared=_declared(repo, options.layer))
        text = json.dumps(report, indent=1, sort_keys=True)
        if options.out:
            Path(options.out).write_text(text + "\n")
            print(json.dumps(report["summary"] | {"straddlers": len(report["summary"]["straddlers"])}, sort_keys=True))
        else:
            print(text)
    else:
        print(cone_store.status(repo, store))
    return 0


def _declared(repo: Path, layer: str) -> list[tuple[str, str]]:
    """Every (test id, ``rests_on`` target) the layer's tests declare, from a collection in a fresh interpreter."""
    import os
    import subprocess
    import tempfile

    from tests._harness.rests_on_manifest import MANIFEST_ENV

    with tempfile.TemporaryDirectory() as scratch:
        manifest = Path(scratch) / "rests_on.json"
        run = subprocess.run([sys.executable, "-O", "-m", "pytest", "--collect-only", "-q", "-p", "no:cacheprovider",
                              "-p", "tests._harness.rests_on_manifest", layer], cwd=repo, capture_output=True,
                             text=True, env={**os.environ, MANIFEST_ENV: str(manifest)})
        if run.returncode != 0:
            raise SystemExit(f"cone_store place: collecting {layer} failed:\n{run.stdout[-2000:]}{run.stderr[-2000:]}")
        rows = json.loads(manifest.read_text(encoding="utf-8"))
    return [(row["nodeid"], target) for row in rows for target in row["rests_on"]]


def _order_module(repo: Path, file: str) -> Any:
    """The module holding a concept order, loaded by path (a test module is not a package member)."""
    import importlib.util

    loaded = importlib.util.spec_from_file_location("_cone_store_order", repo / file)
    if loaded is None or loaded.loader is None:
        raise SystemExit(f"cone_store place: cannot load the order from {file}")
    module = importlib.util.module_from_spec(loaded)
    loaded.loader.exec_module(module)
    return module


if __name__ == "__main__":
    sys.exit(main())
