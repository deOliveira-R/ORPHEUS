"""``python -m tools.concept_table report | seed``: the counts of the concept table, and the fingerprint seed rewritten.

The first line is the one a reader compares with the specification: the classes, the declared classes, the concepts and
the length of each of the nine pinned lists. The lines after it give the entries of the lists by wave, the references
the checks do not judge (those owned by module-level functions) and the files skipped.

``seed`` rewrites ``concept_seed.toml`` after a class on the twin or weld list is edited or a class joins a list: the
entry of every class on a list today takes that class's fingerprint, and every other entry is kept as it is. It never
removes an entry and never touches the entry of a class that is off the lists, so it cannot make a relabel pass.
"""

from __future__ import annotations

import collections
import sys
from pathlib import Path

from .checks import fingerprint, measure_twin, measure_weld
from .debt import SeedEntry, dumps_seed, load_debt, load_seed, wave_of
from .paths import DEBT_FILE, ROOT, SEED_FILE, TABLE_FILE
from .population import Population, read_sources
from .refs import read_function_edges
from .table import load_table


def report(root: Path = ROOT) -> list[str]:
    """The lines of the report on the tree at ``root``."""
    sources, skipped = read_sources(root)
    population = Population.from_sources(sources, skipped)
    concepts = load_table(root / TABLE_FILE)
    debt = load_debt(root / DEBT_FILE)
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


def reseed(root: Path = ROOT) -> str:
    """The seed file's text for the tree at ``root``: listed classes at today's fingerprint, every other entry kept unless its file left the tree."""
    sources, skipped = read_sources(root)
    population = Population.from_sources(sources, skipped)
    path = root / SEED_FILE
    header = "".join(line for line in path.read_text().splitlines(keepends=True) if line.startswith("#"))
    listed = {d.key for d in measure_twin(population) | measure_weld(population)}
    kept = [e for e in load_seed(path) if e.key not in listed and e.key.file in sources]
    return dumps_seed([*kept, *(SeedEntry(k, fingerprint(sources[k.file], k.qualname)) for k in listed)], header)


def main(argv: list[str]) -> int:
    if argv == ["report"]:
        sys.stdout.write("\n".join(report()) + "\n")
        return 0
    if argv == ["seed"]:
        (ROOT / SEED_FILE).write_text(reseed())
        return 0
    sys.stderr.write("usage: python -m tools.concept_table report | seed\n")
    return 2


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
