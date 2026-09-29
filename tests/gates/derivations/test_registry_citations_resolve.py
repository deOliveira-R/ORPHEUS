r"""The Sood and Atalay registries cite by :class:`~orpheus.data.Citation`, and every key resolves.

P1 step 2c of the reference-solution campaign (#405). The free-text
``Provenance`` record retired: a case's ``problem`` is the citation of the
source that defined the benchmark problem, and its truth's ``sources`` are the
citations of the published values. Whether a citation's key exists in
``docs/refs.bib`` cannot be a runtime check (``docs/`` does not ship), so it
is this gate.

What each row asserts, and what reddens it:

* the reader control: the ``refs.bib`` key reader finds a known key, reports
  an unknown one, and agrees with an independent BibTeX parser (``pybtex``) on
  the whole key set, so a green row (a) cannot come from a reader that finds
  nothing or everything;
* (a) every citation reachable from a case, found by walking the case's
  fields (not by naming them), names a key in ``docs/refs.bib``;
* (b) the problem's source: the 47 Sood cases cite ``SoodForsterParsons2003``
  at ``problem N`` with N Sood's problem number, and the Atalay cases cite
  ``Atalay1997`` at the table that prints them;
* (c) every case's truth credits at least one source, and one of them is the
  case's own paper at a place in it;
* (d) no citation names ``SoodLA13511_1999``: the registry cites the 2003
  edition (the user's ruling, 2026-09-29). The 1999 key IS in ``refs.bib``,
  so row (a) cannot see this;
* (e) Sood's problem numbers are distinct across the cases.

The expected problem numbers are a transcription of the problem column of
Sood, Forster and Parsons (2003), checked 47 of 47 against the registry's
pre-carve ``problem_number`` field on 2026-09-29; the Atalay table numbers
are the pre-carve ``Provenance.paper_table``. They are the carve's
value-preservation anchor: the carve moves each number into a locator
string, and a number moved onto the wrong case is red here.
"""

from __future__ import annotations

import ast
import dataclasses
import functools
import re
from collections import Counter
from collections.abc import Iterable, Iterator, Mapping
from pathlib import Path

import pytest

import orpheus.derivations.continuous.sood_registry as sood_registry
from orpheus.data import Citation

# The two collections, named in ONE place: P1 step 2c's step 3 renames the
# first to ``SOOD2003_CASES``, a one-line edit here.
from orpheus.derivations.continuous.sood_registry import (
    ATALAY_ALL_CASES as ATALAY_CASES,
    LA13511_CASES as SOOD_CASES,
)

pytestmark = pytest.mark.foundation

_SELF = "tests/gates/derivations/test_registry_citations_resolve.py"
_LAWS = "tests/gates/data/test_citation.py"

_REFS_BIB = Path(__file__).resolve().parents[3] / "docs" / "refs.bib"

SOOD_2003 = "SoodForsterParsons2003"
SOOD_1999 = "SoodLA13511_1999"
ATALAY_1997 = "Atalay1997"

#: Sood's problem number per case: the problem column of Sood, Forster and
#: Parsons (2003), whose problem numbers are the 1999 report's (1-75).
_SOOD_PROBLEM: dict[str, int] = {
    "PUa-1-0-IN": 1, "PU-2-0-IN": 44, "Ua-1-0-SL": 12, "Ua-1-0-CY": 13,
    "Ua-1-0-SP": 14, "PUb-1-0-IN": 5, "Ua-1-0-IN": 11, "Ub-1-0-IN": 15,
    "Uc-1-0-IN": 17, "Ud-1-0-IN": 19, "UD2O-1-0-IN": 21, "Ue-1-0-IN": 29,
    "PU-1-1-IN": 31, "UD2Oa-1-1-IN": 38, "UD2Ob-1-1-IN": 40, "UD2Oc-1-1-IN": 42,
    "U-2-0-IN": 47, "UAL-2-0-IN": 50, "URRa-2-0-IN": 53, "URRb-2-0-IN": 56,
    "URRc-2-0-IN": 57, "URRd-2-0-IN": 62, "UD2O-2-0-IN": 67, "URR-3-0-IN": 74,
    "URR-6-0-IN": 75, "PUa-1-0-SL": 2, "PUb-1-0-SL": 6, "UD2O-1-0-SL": 22,
    "PUb-1-0-SP": 8, "UD2O-1-0-SP": 24, "PUa-1-1-SL": 32, "PUb-1-1-SL": 34,
    "UD2Oa-1-1-SP": 39, "UD2Ob-1-1-SP": 41, "UD2Oc-1-1-SP": 43, "PUb-1-0-CY": 7,
    "UD2O-1-0-CY": 23, "PU-2-0-SL": 45, "PU-2-0-SP": 46, "U-2-0-SL": 48,
    "U-2-0-SP": 49, "UAL-2-0-SL": 51, "UAL-2-0-SP": 52, "URRa-2-0-SL": 54,
    "URRa-2-0-SP": 55, "UD2O-2-0-SL": 68, "UD2O-2-0-SP": 69,
}

#: The Atalay (1997) table printing each case. The f1 = 0 sphere retired in
#: step 2 as a twin of Sood's problem 14 (Atalay prints no sphere at f1 = 0).
_ATALAY_TABLE: dict[str, int] = {
    "atalay-1997-slab-c1.30-R0.00-f1_0.00": 2,
    "atalay-1997-slab-c1.30-R0.25-f1_0.00": 2,
    "atalay-1997-slab-c1.30-R0.50-f1_0.00": 2,
    "atalay-1997-slab-c1.30-R0.75-f1_0.00": 2,
    "atalay-1997-slab-c1.30-R0.00-f1_0.10": 3,
    "atalay-1997-slab-c1.30-R0.50-f1_0.10": 3,
}

_PROBLEM_LOCATOR = re.compile(r"problem (\d+)")
_TABLE_LOCATOR = re.compile(r"Table (\d+)\b")


def _require(condition: bool, message: str) -> None:
    """A ``-O``-surviving assertion (not a bare ``assert``)."""
    if not condition:
        pytest.fail(message)


# ── the instruments ────────────────────────────────────────────────────────


_ENTRY = re.compile(r"^@(\w+)\s*\{\s*([^,\s]+)\s*,", re.MULTILINE)
_NOT_AN_ENTRY = frozenset({"comment", "string", "preamble"})


def _bib_keys(text: str) -> frozenset[str]:
    """The entry keys of a BibTeX file, read from its ``@type{key,`` lines."""
    return frozenset(key for kind, key in _ENTRY.findall(text) if kind.lower() not in _NOT_AN_ENTRY)


@functools.cache
def _refs_bib_keys() -> frozenset[str]:
    return _bib_keys(_REFS_BIB.read_text(encoding="utf-8"))


def _unresolved(citations: Iterable[tuple[str, Citation]], keys: frozenset[str]) -> list[str]:
    """The ``path: citation`` of every citation whose key is not in ``keys``."""
    return [f"{path}: {citation!r}" for path, citation in citations if citation.bibkey not in keys]


def _cases(collection: object) -> tuple:
    """A collection's cases, whether it is a mapping by case id or a tuple."""
    if isinstance(collection, Mapping):
        return tuple(collection.values())
    if isinstance(collection, (tuple, list)):
        return tuple(collection)
    raise TypeError(f"not a case collection: {type(collection).__name__}")


def _citations(value: object, path: str = "case") -> Iterator[tuple[str, Citation]]:
    """Every :class:`Citation` reachable from ``value``, with its field path.

    Walks the registry's own dataclasses (the case, its truth), tuples, lists
    and mappings. It does not descend into other packages' objects (a
    ``Mixture`` holds no citation), so a new citation-bearing field on a
    registry type is found without editing this gate.
    """
    if isinstance(value, Citation):
        yield path, value
    elif (
        dataclasses.is_dataclass(value)
        and not isinstance(value, type)
        and type(value).__module__.startswith(sood_registry.__name__)
    ):
        for field in dataclasses.fields(value):
            yield from _citations(getattr(value, field.name), f"{path}.{field.name}")
    elif isinstance(value, (tuple, list)):
        for i, item in enumerate(value):
            yield from _citations(item, f"{path}[{i}]")
    elif isinstance(value, Mapping):
        for key, item in value.items():
            yield from _citations(item, f"{path}[{key!r}]")


def _duplicates(values: Iterable[object]) -> dict[object, int]:
    """Each value occurring more than once, with its count."""
    return {value: n for value, n in Counter(values).items() if n > 1}


_SOOD = _cases(SOOD_CASES)
_ATALAY = _cases(ATALAY_CASES)
_ALL = _SOOD + _ATALAY
_SOOD_IDS = [case.case_id for case in _SOOD]
_ATALAY_IDS = [case.case_id for case in _ATALAY]
_ALL_IDS = _SOOD_IDS + _ATALAY_IDS
_BY_ID = {case.case_id: case for case in _ALL}


# ── the controls: the instruments read positive ────────────────────────────


def test_refs_bib_reader_positive_control() -> None:
    """The key reader finds a known key, reports an unknown one, and agrees
    with ``pybtex`` (the parser ``sphinxcontrib.bibtex`` renders with) on the
    whole key set. The 1999 report's key is present, which is why row (d)
    exists beside row (a)."""
    from pybtex.database import parse_file

    _require(_REFS_BIB.is_file(), f"no bibliography at {_REFS_BIB}")
    keys = _refs_bib_keys()
    independent = frozenset(parse_file(str(_REFS_BIB), bib_format="bibtex").entries.keys())
    _require(len(keys) > 0, "the reader found no key")
    _require(
        keys == independent,
        f"the reader and pybtex disagree: reader only {sorted(keys - independent)}, "
        f"pybtex only {sorted(independent - keys)}",
    )
    _require(SOOD_2003 in keys and ATALAY_1997 in keys and SOOD_1999 in keys, "a known key is missing")
    probe = [("known", Citation(SOOD_2003, "problem 4")), ("unknown", Citation("NoSuchWork2099"))]
    reported = _unresolved(probe, keys)
    _require(
        reported == ["unknown: Citation(bibkey='NoSuchWork2099', locator=None)"],
        f"the resolution check reported {reported}",
    )


@pytest.mark.rests_on(
    f"{_LAWS}::test_legal_bibkeys_construct",
    f"{_LAWS}::test_refuses_illegal_bibkey",
    f"{_SELF}::test_refs_bib_reader_positive_control",
)
def test_every_refs_bib_key_is_a_legal_citation() -> None:
    """The ``Citation`` key grammar covers the bibliography: every key in
    ``docs/refs.bib`` constructs a citation. A key added outside the grammar
    is red here, before any registry tries to cite it."""
    refused = []
    for key in sorted(_refs_bib_keys()):
        try:
            Citation(key)
        except ValueError as error:
            refused.append(f"{key}: {error}")
    _require(not refused, f"{len(refused)} refs.bib keys refused by Citation: {refused}")


def test_the_collections_and_the_reference_tables_name_the_same_cases() -> None:
    """The input count: the reference tables and the collections hold the
    same cases, and case ids are distinct across both collections. A case
    added or retired owes a row in its table."""
    _require(len(_SOOD) == len(_SOOD_PROBLEM) == 47, f"{len(_SOOD)} Sood cases")
    _require(set(_SOOD_IDS) == set(_SOOD_PROBLEM), f"Sood ids differ: {set(_SOOD_IDS) ^ set(_SOOD_PROBLEM)}")
    _require(len(_ATALAY) == len(_ATALAY_TABLE) > 0, f"{len(_ATALAY)} Atalay cases")
    _require(
        set(_ATALAY_IDS) == set(_ATALAY_TABLE),
        f"Atalay ids differ: {set(_ATALAY_IDS) ^ set(_ATALAY_TABLE)}",
    )
    _require(len(set(_ALL_IDS)) == len(_ALL_IDS), f"duplicate case ids: {_duplicates(_ALL_IDS)}")
    _require(all(1 <= n <= 75 for n in _SOOD_PROBLEM.values()), "a problem number is outside 1-75")


@pytest.mark.rests_on(f"{_SELF}::test_the_collections_and_the_reference_tables_name_the_same_cases")
def test_citation_walker_reaches_the_problem_and_every_source() -> None:
    """The walker finds, on every case, the citations the schema names:
    ``case.problem`` and each element of ``case.truth.sources``."""
    missed = []
    for case in _ALL:
        found = set(_citations(case))
        named = {("case.problem", case.problem)} | {
            (f"case.truth.sources[{i}]", source) for i, source in enumerate(case.truth.sources)
        }
        if not named <= found:
            missed.append(f"{case.case_id}: {sorted(map(str, named - found))}")
    _require(not missed, f"the walker missed citations on {len(missed)} cases: {missed}")


# ── (a) every key resolves ─────────────────────────────────────────────────


@pytest.mark.rests_on(
    f"{_SELF}::test_refs_bib_reader_positive_control",
    f"{_SELF}::test_citation_walker_reaches_the_problem_and_every_source",
)
@pytest.mark.parametrize("case_id", _ALL_IDS)
def test_every_registry_citation_resolves_in_refs_bib(case_id: str) -> None:
    """Every citation reachable from the case names a key in ``docs/refs.bib``."""
    citations = list(_citations(_BY_ID[case_id]))
    _require(len(citations) >= 2, f"{case_id}: {len(citations)} citations found (a problem and a source at least)")
    missing = _unresolved(citations, _refs_bib_keys())
    _require(not missing, f"{case_id}: keys not in docs/refs.bib: {missing}")


# ── (b) the problem's source ───────────────────────────────────────────────


@pytest.mark.rests_on(f"{_SELF}::test_the_collections_and_the_reference_tables_name_the_same_cases")
@pytest.mark.parametrize("case_id", _SOOD_IDS)
def test_sood_case_cites_its_problem_in_the_2003_paper(case_id: str) -> None:
    """``problem == Citation("SoodForsterParsons2003", "problem N")``, N from
    Sood's problem column."""
    expected = Citation(SOOD_2003, f"problem {_SOOD_PROBLEM[case_id]}")
    problem = _BY_ID[case_id].problem
    _require(problem == expected, f"{case_id}: problem is {problem!r}, expected {expected!r}")


@pytest.mark.rests_on(f"{_SELF}::test_the_collections_and_the_reference_tables_name_the_same_cases")
@pytest.mark.parametrize("case_id", _ATALAY_IDS)
def test_atalay_case_cites_its_table_in_atalay_1997(case_id: str) -> None:
    """``problem`` cites ``Atalay1997`` at the table that prints the case
    (a locator opening ``Table T``)."""
    problem = _BY_ID[case_id].problem
    _require(problem.bibkey == ATALAY_1997, f"{case_id}: problem cites {problem.bibkey!r}")
    table = _TABLE_LOCATOR.match(problem.locator or "")
    if table is None:
        pytest.fail(f"{case_id}: locator {problem.locator!r} does not open 'Table T'")
    _require(
        int(table.group(1)) == _ATALAY_TABLE[case_id],
        f"{case_id}: cites Table {table.group(1)}, expected Table {_ATALAY_TABLE[case_id]}",
    )


@pytest.mark.rests_on(f"{_SELF}::test_the_collections_and_the_reference_tables_name_the_same_cases")
def test_atalay_problem_citations_are_distinct() -> None:
    """Two Atalay cases citing the same place would be one problem twice."""
    repeated = _duplicates(case.problem for case in _ATALAY)
    _require(not repeated, f"Atalay problem citations repeated: {repeated}")


# ── (c) the published values are credited ──────────────────────────────────


@pytest.mark.rests_on(f"{_SELF}::test_the_collections_and_the_reference_tables_name_the_same_cases")
@pytest.mark.parametrize("case_id", _ALL_IDS)
def test_truth_sources_are_nonempty_citations(case_id: str) -> None:
    """``truth.sources`` is a non-empty tuple of citations."""
    sources = _BY_ID[case_id].truth.sources
    _require(isinstance(sources, tuple), f"{case_id}: sources is a {type(sources).__name__}")
    _require(len(sources) > 0, f"{case_id}: no source credits the published values")
    bad = [s for s in sources if not isinstance(s, Citation)]
    _require(not bad, f"{case_id}: sources hold non-citations {bad!r}")


def test_truth_refuses_empty_sources() -> None:
    """The witness of the truth type's admission guard: a truth rebuilt with
    one source constructs, and with none is refused. The row above cannot
    see a removed guard (no registry case is built without a source); this
    row can."""
    truth = _SOOD[0].truth
    one = dataclasses.replace(truth, sources=(Citation(SOOD_2003, "p. 69"),))
    _require(one.sources == (Citation(SOOD_2003, "p. 69"),), f"the positive leg built {one.sources!r}")
    with pytest.raises(ValueError, match="sources is empty"):
        dataclasses.replace(truth, sources=())


@pytest.mark.rests_on(f"{_SELF}::test_truth_sources_are_nonempty_citations")
@pytest.mark.parametrize("case_id", _ALL_IDS)
def test_truth_sources_place_the_value_in_the_case_paper(case_id: str) -> None:
    """One source is the case's own paper at a place in it (Sood's table or
    page for a Sood case, Atalay's table for an Atalay case): the value is
    printed there, whatever primary source the paper credits."""
    paper = SOOD_2003 if case_id in _SOOD_PROBLEM else ATALAY_1997
    sources = _BY_ID[case_id].truth.sources
    placed = [s for s in sources if s.bibkey == paper and s.locator is not None]
    _require(bool(placed), f"{case_id}: no source cites {paper} at a place: {sources!r}")


# ── (d) the registry cites the 2003 edition ────────────────────────────────


@pytest.mark.rests_on(f"{_SELF}::test_citation_walker_reaches_the_problem_and_every_source")
@pytest.mark.parametrize("case_id", _ALL_IDS)
def test_no_registry_citation_names_the_1999_report(case_id: str) -> None:
    """No citation reachable from the case names ``SoodLA13511_1999``, and the
    case's notes do not cite the 1999 report by its number (the successor of
    the retired provenance gate's string check)."""
    case = _BY_ID[case_id]
    cited_1999 = [f"{path}: {c!r}" for path, c in _citations(case) if c.bibkey == SOOD_1999]
    _require(not cited_1999, f"{case_id}: cites the 1999 edition: {cited_1999}")
    _require("LA-13511" not in case.notes, f"{case_id}: notes cite the 1999 report: {case.notes!r}")


# ── (e) Sood's problem numbers are distinct ────────────────────────────────


@pytest.mark.rests_on(f"{_SELF}::test_the_collections_and_the_reference_tables_name_the_same_cases")
def test_sood_problem_numbers_are_distinct() -> None:
    """Each Sood problem is one case. Read from the locators themselves, not
    from the reference table."""
    _require(_duplicates([1, 2, 2, 3]) == {2: 2}, "the duplicate finder is broken")
    numbers = []
    for case in _SOOD:
        match = _PROBLEM_LOCATOR.fullmatch(case.problem.locator or "")
        if match is None:
            pytest.fail(f"{case.case_id}: locator {case.problem.locator!r} is not 'problem N'")
        numbers.append(int(match.group(1)))
    _require(len(numbers) == len(_SOOD), "a Sood case was not read")
    _require(not _duplicates(numbers), f"Sood problem numbers repeated: {_duplicates(numbers)}")


# ── the retired schema stays retired ───────────────────────────────────────


@pytest.mark.parametrize("case_id", _ALL_IDS)
def test_provenance_fields_are_retired(case_id: str) -> None:
    """``Provenance`` and ``problem_number`` retired into citations (step 2);
    the flat fields retired in step 2b. ``notes`` moved onto the case."""
    case = _BY_ID[case_id]
    back = [f for f in ("provenance", "problem_number", "sood_table", "primary_reference") if hasattr(case, f)]
    _require(not back, f"{case_id}: retired fields are back: {back}")
    _require(isinstance(case.notes, str), f"{case_id}: notes is a {type(case.notes).__name__}")
    _require(not hasattr(sood_registry, "Provenance"), "sood_registry still exports Provenance")


def test_the_case_schema_has_its_own_module() -> None:
    """The schema lives in ``sood_registry/case.py``, and ``atalay1997.py``
    imports it from there, not from the Sood cases' module (the weld step 2
    dissolves)."""
    case_module = f"{sood_registry.__name__}.case"
    for case in (_SOOD[0], _ATALAY[0]):
        _require(
            type(case).__module__ == case_module and type(case.truth).__module__ == case_module,
            f"{case.case_id}: schema from {type(case).__module__} / {type(case.truth).__module__}",
        )
    source = (Path(sood_registry.__file__).parent / "atalay1997.py").read_text(encoding="utf-8")
    relative = [node for node in ast.walk(ast.parse(source)) if isinstance(node, ast.ImportFrom) and node.level == 1]
    from_case = {alias.name for node in relative if node.module == "case" for alias in node.names}
    _require({"La13511Case", "La13511Truth"} <= from_case, f"atalay1997 imports {sorted(from_case)} from .case")
    from_sood = [node.module for node in relative if node.module in {"la13511", "sood2003"}]
    _require(not from_sood, f"atalay1997 still imports from the Sood module: {from_sood}")
