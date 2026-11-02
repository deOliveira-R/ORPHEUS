"""``python -m tools.cone_store {capture,select,status} --repo R --store S [--package P] [--memo-root M] [ids]``."""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    started = time.time_ns()
    parser = argparse.ArgumentParser(prog="python -m tools.cone_store", description=__doc__)
    parser.add_argument("command", choices=("capture", "select", "status"))
    parser.add_argument("ids", nargs="*", help="pytest node ids or paths (default: the project's testpaths)")
    parser.add_argument("--repo", default=".", help="the repository root (default: the working directory)")
    parser.add_argument("--store", default=None, help="the SQLite store (default: <repo>/.cache/cone_store/cones.sqlite)")
    parser.add_argument("--package", default="orpheus", help="the first-party package under the repository")
    parser.add_argument("--memo-root", default=None, help="the traced memo's store (default: the memo's own)")
    parser.add_argument("--route", choices=("auto", "fresh"), default="auto",
                        help="capture: auto forks where the memo's guard allows; fresh runs every test in its own interpreter")
    parser.add_argument("--why", action="store_true", help="select: print the reasons beside each id")
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
    else:
        print(cone_store.status(repo, store))
    return 0


if __name__ == "__main__":
    sys.exit(main())
