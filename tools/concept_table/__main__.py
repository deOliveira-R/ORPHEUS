"""``python -m tools.concept_table report``: the counts of the concept table, the debt and the population.

The first line is the one a reader compares with the specification: the classes, the declared classes, the concepts and
the length of each of the nine pinned lists. The lines after it give the entries of the lists by wave, the references
the checks do not judge (those owned by module-level functions) and the files skipped.
"""

from __future__ import annotations

import collections
import sys
from pathlib import Path

from .debt import load_debt, wave_of
from .population import read_tree
from .refs import read_function_edges
from .table import load_table

ROOT = Path(__file__).resolve().parents[2]
ARCH = ROOT / "docs" / "architecture"


def report(root: Path = ROOT) -> list[str]:
    """The lines of the report on the tree at ``root``."""
    population = read_tree(root)
    concepts = load_table(root / "docs" / "architecture" / "concept_table.toml")
    debt = load_debt(root / "docs" / "architecture" / "concept_debt.toml")
    sources = {r.key.file: (root / r.key.file).read_text() for r in population.rows}
    lists = debt.lists()
    waves: collections.Counter[str] = collections.Counter()
    for entries in lists.values():
        for entry in entries:
            waves.update(sorted(wave_of(entry, concepts, population)))
    declared = sum(1 for r in population.rows if r.declaration is not None)
    first = " ".join([f"classes {len(population.rows)}", f"declared {declared}", f"concepts {len(concepts)}",
                      *(f"{name} {len(entries)}" for name, entries in lists.items())])
    by_wave = " ".join(f"{w} {n}" for w, n in sorted(waves.items()))
    return [first, f"entries by wave: {by_wave}",
            f"function-owned references, not judged: {len(read_function_edges(sources, population))}",
            f"files read {population.files_read}, files skipped under derivations/ {population.files_skipped_derivations}"]


def main(argv: list[str]) -> int:
    if argv != ["report"]:
        sys.stderr.write("usage: python -m tools.concept_table report\n")
        return 2
    sys.stdout.write("\n".join(report()) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
