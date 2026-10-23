r"""The admission pin: one row per (method, value class, site class), pinned EXACTLY (wave 2, section 4.4).

Spec: ``.claude/plans/wave2_boundary_carve.md`` section 4 ("the admission table, as data"), 4.4 (the pin), 6.2 S2a.

This file is created by slice S2a with the MECHANISM's rows only: a test-local table whose rows copy, by hand, the SN
columns of section 4.2 (the values are typed here, never read from a production table, so a change to any admission
edits this file). Slices S2b, S3 and S4 each APPEND their method's rows to :data:`ROWS` and register the builder of
their real table in :data:`TABLES`; the same parametrised checker then pins them. A row is
``Row(method, value, site, expected)`` where ``value`` is a key of ``_admission_helpers.values()`` (the value classes
of 4.2: vacuum, specular of weight 1.0 and 0.5, white, both chains, a pure inflow, a response over an inflow), ``site``
is ``"CART"``, ``"RAD_OUTER"`` or ``"RAD_INNER"``, and ``expected`` is the admission's SUMMARY
(:func:`summary`): ``("R", faithful)``, ``("A", names, issues)`` or ``("F", axis, issue)`` (the reason text is prose
and is not pinned).

The claims (``vv-principles``), all REFERENCE rows (the expected values are the plan's table typed by hand):

* ``test_the_pin`` (per row): the table's ``admit`` equals the row.
* ``test_the_pin_covers_every_cell``: for each method the rows are exactly the product of the value classes and the
  sites, so a new value class or site owes its row (a pin with a hole pins nothing there).
* ``test_both_chains_is_the_meet``: the "both chains" row equals ``Admission.meet`` of the two single-chain rows (and
  a response over a source, the meet of the response and the source), at every site, so a change to either single
  row edits the mixed row (section 4.2's own rule).
* ``test_a_plumbing_refusal_names_a_cited_issue``: every plumbing row carries an issue from :data:`CITED_ISSUES`, the
  list for the orchestrator's ``gh issue view --json state`` read (section 4.4: ``Refused.__post_init__`` cannot check
  liveness; ``[M]`` 9 of 9 OPEN on 2026-10-22).
* ``test_the_pin_has_teeth``: the checker applied to a copy of the table with one admission drifted (the diffuse
  chain admitted on a radial outer face) reports exactly that cell (X1: a pin that cannot fail pins nothing).

First reds, at the base tip ``4ec95bc6``: ``ModuleNotFoundError`` (the first read of the package). Mutations: the
drifted table of ``test_the_pin_has_teeth``; a ``meet`` that returns the greater; a table whose zero row is consulted
for every wall (the vacuum and inflow rows are read by the functor row of the mechanism file).
"""

from __future__ import annotations

import dataclasses
import itertools
from collections.abc import Callable
from typing import Any

import pytest

from tests.gates._content_identity_helpers import require
from tests.gates.transport._admission_helpers import SITES, make_table, site, values
from tests.gates.transport._boundary_conditions_helpers import bc, node

pytestmark = pytest.mark.foundation

_MECHANISM = node("test_admission_mechanism.py", "test_the_lattice_of_admissions")

#: The issue numbers the pin cites, for the orchestrator's liveness read (section 4.4). Appended to as rows are.
CITED_ISSUES: frozenset[int] = frozenset({511, 513, 514, 546, 547, 653, 654, 655, 656})


@dataclasses.dataclass(frozen=True)
class Row:
    method: str
    value: str
    site: str
    expected: tuple


def R(faithful: bool = True) -> tuple:
    return ("R", faithful)


def F(axis: str, issue: int | None = None) -> tuple:
    return ("F", axis, issue)


def summary(admission: Any) -> tuple:
    """The pinned shape of an admission: tier and, for the lesser tiers, the names/issues or the axis and issue."""
    if isinstance(admission, bc.Realized):
        return ("R", admission.faithful)
    if isinstance(admission, bc.Approximated):
        return ("A", tuple(sorted(admission.names)), tuple(sorted(admission.issues)))
    require(isinstance(admission, bc.Refused), f"{admission!r} is not an admission")
    return ("F", admission.axis.name, admission.issue)


#: The mechanism's rows: the SN columns of section 4.2, typed by hand (CART, RAD_OUTER, RAD_INNER per value class).
_MECHANISM_COLUMNS: dict[str, tuple[tuple, tuple, tuple]] = {
    "vacuum":                        (R(), R(), F("PLUMBING", 511)),
    "specular 1.0":                  (R(), R(), R()),
    "specular 0.5":                  (R(), R(), F("PLUMBING", 511)),
    "white 0.5":                     (R(), F("PLUMBING", 653), F("PLUMBING", 511)),
    "both 0.3 0.5":                  (R(), F("PLUMBING", 653), F("PLUMBING", 511)),
    "inflow 1.0":                    (R(), F("PLUMBING", 653), F("PLUMBING", 511)),
    "specular 0.5 over inflow 1.0":  (R(), F("PLUMBING", 653), F("PLUMBING", 511)),
}

#: Slices S2b, S3 and S4 append their rows HERE (``ROWS += (Row(...), ...)``) and register their real table's builder in
#: ``TABLES``; the checkers below pin them without edit.
ROWS: tuple[Row, ...] = tuple(
    Row("mechanism", value, where, expected)
    for value, column in _MECHANISM_COLUMNS.items()
    for where, expected in zip(SITES, column)
)

TABLES: dict[str, Callable[[], Any]] = {"mechanism": lambda: make_table().table}


def _table(method: str) -> Any:
    return TABLES[method]()


def pin_disagreements(table: Any, rows: tuple[Row, ...]) -> list[str]:
    """The rows whose pinned summary differs from ``table.admit``'s (an empty list is a pinned table)."""
    declared = values()
    out = []
    for row in rows:
        got = summary(table.admit(declared[row.value], site(row.site)))
        if got != row.expected:
            out.append(f"{row.method} / {row.value} @ {row.site}: admit gives {got}, the pin says {row.expected}")
    return out


def _ids(rows: tuple[Row, ...]) -> list[str]:
    return [f"{r.method}-{r.value}-{r.site}".replace(" ", "_") for r in rows]


@pytest.mark.parametrize("row", ROWS, ids=_ids(ROWS))
@pytest.mark.rests_on(_MECHANISM)
def test_the_pin(row: Row) -> None:
    """The table's ``admit`` equals the row (summary: tier, names/issues or axis/issue)."""
    got = summary(_table(row.method).admit(values()[row.value], site(row.site)))
    require(got == row.expected, f"{row.method} / {row.value} @ {row.site}: admit gives {got}, the pin says {row.expected}")


@pytest.mark.rests_on(_MECHANISM)
def test_the_pin_covers_every_cell() -> None:
    """For each method the rows are exactly (value classes) x (sites): no cell unpinned, none pinned twice."""
    declared = set(values())
    require(len(declared) == 7 and set(SITES) == {"CART", "RAD_OUTER", "RAD_INNER"}, "the value classes or sites changed")
    for method in {row.method for row in ROWS}:
        cells = [(row.value, row.site) for row in ROWS if row.method == method]
        require(len(cells) == len(set(cells)), f"{method}: a cell is pinned twice")
        require(set(cells) == set(itertools.product(declared, SITES)), f"{method}: unpinned {sorted(set(itertools.product(declared, SITES)) - set(cells))}")


@pytest.mark.rests_on(_MECHANISM)
def test_both_chains_is_the_meet() -> None:
    """At every site, ``admit(both chains)`` is ``Admission.meet(admit(specular), admit(white))`` and ``admit(a
    response over a source)`` is ``meet(admit(response), admit(source))``: the mixed rows are derived, so a change to a
    single row edits them. The pieces are read from the same table, not from the pin (the pin rows are pinned
    separately by ``test_the_pin``).

    Reddens on: a table whose ``admit`` for a wall consults the chains in another order than the meet's, or a
    ``meet`` that returns the greater.
    """
    table, W = _table("mechanism"), bc.WallCondition
    for kind in SITES:
        where = site(kind)
        spec, white = table.admit(W.specular(0.3), where), table.admit(W.white(0.5), where)
        require(summary(table.admit(W.from_amplitudes(0.3, 0.5), where)) == summary(bc.Admission.meet(spec, white)), f"both chains @ {kind} is not the meet of the two single-chain rows")
        response = table.admit(W.specular(0.5), where)
        source = table.admit(bc.InflowTerms(1.0, ()), where)
        require(summary(table.admit(W.specular(0.5) + W.uniform(1.0), where)) == summary(bc.Admission.meet(response, source)), f"a response over a source @ {kind} is not the meet")
    require(summary(table.admit(W.from_amplitudes(0.3, 0.5), site("RAD_OUTER"))) == F("PLUMBING", 653), "the radial outer mixed row is the diffuse refusal")


@pytest.mark.rests_on(_MECHANISM)
def test_a_plumbing_refusal_names_a_cited_issue() -> None:
    """Every plumbing refusal a registered table RETURNS, over the whole grid of value classes and sites, carries an
    issue number of ``CITED_ISSUES``; so does every plumbing row of the pin; the list holds nine ints (the orchestrator
    reads each with ``gh``). The returned refusals are read from the tables, not from the rows, so a realizer that
    cites a number the pin does not list reds here.

    Reddens on: a realizer citing an issue outside ``CITED_ISSUES`` (e.g. the diffuse refusal's 653 changed to 999 in
    ``_admission_helpers.make_table``: this row and the pinned rows that read it); a plumbing row of the pin without a
    cited number.
    """
    plumbing = [row for row in ROWS if row.expected[0] == "F" and row.expected[1] == "PLUMBING"]
    require(len(plumbing) >= 10, f"the pin holds {len(plumbing)} plumbing rows")
    for row in plumbing:
        require(isinstance(row.expected[2], int) and row.expected[2] in CITED_ISSUES, f"{row}: a plumbing refusal cites {row.expected[2]!r}, not an issue of CITED_ISSUES")
    require(all(isinstance(n, int) for n in CITED_ISSUES) and len(CITED_ISSUES) == 9, "CITED_ISSUES is not the nine issues of section 4.4")

    declared, returned = values(), []
    for method in sorted({row.method for row in ROWS}):
        table = _table(method)
        for label, value in declared.items():
            for kind in SITES:
                admission = table.admit(value, site(kind))
                if isinstance(admission, bc.Refused) and admission.axis is bc.RefusalAxis.PLUMBING:
                    returned.append((method, label, kind, admission.issue))
    require(len(returned) >= 10, f"the tables returned {len(returned)} plumbing refusals")
    stray = [entry for entry in returned if entry[3] not in CITED_ISSUES]
    require(not stray, f"plumbing refusals citing an issue outside CITED_ISSUES: {stray}")


@pytest.mark.rests_on(_MECHANISM)
def test_the_pin_has_teeth() -> None:
    """The checker applied to the honest table reports nothing; applied to a copy with the diffuse chain ADMITTED on
    a radial outer face (a drift of one admission) it reports exactly the two cells that read that admission (the
    white row and the both-chains row at RAD_OUTER; the source rows do not read the diffuse unit and stay agreed).

    Reddens on: a checker that compares nothing (an empty report on a drifted table).
    """
    require(pin_disagreements(_table("mechanism"), ROWS) == [], f"the honest table disagrees with its pin: {pin_disagreements(_table('mechanism'), ROWS)}")

    drifted = bc.RealizerTable("drifted", state_cone=bc.StateCone.NONNEGATIVE_FLUX, unregistered_issue=654)
    R_, Refused = bc.Realized, bc.Refused

    class Spec(bc.ResponseRealizer, table=drifted, realizes=bc.SpecularChain):
        status = R_(True)

        def admit(self, unit, weight, site_, space=None):
            if site_.site_class is bc.SiteClass.RAD_INNER:
                return R_(True) if weight == 1.0 else Refused(bc.RefusalAxis.PLUMBING, "inner", 511)
            return self.status

        def realize(self, unit, weight, site_, space):
            raise NotImplementedError

    class Diff(bc.ResponseRealizer, table=drifted, realizes=bc.DiffuseChain):
        status = R_(True)

        def admit(self, unit, weight, site_, space=None):
            return Refused(bc.RefusalAxis.PLUMBING, "inner", 511) if site_.site_class is bc.SiteClass.RAD_INNER else self.status  # RAD_OUTER DRIFTED to admitted

        def realize(self, unit, weight, site_, space):
            raise NotImplementedError

    class Zero(bc.ResponseRealizer, table=drifted, realizes=bc.ZeroNetwork):
        status = R_(True)

        def admit(self, unit, weight, site_, space=None):
            return Refused(bc.RefusalAxis.PLUMBING, "inner", 511) if site_.site_class is bc.SiteClass.RAD_INNER else self.status

        def realize(self, unit, weight, site_, space):
            raise NotImplementedError

    class Src(bc.SourceRealizer, table=drifted):
        status = R_(True)

        def admit(self, terms, site_, space=None):
            match site_.site_class:
                case bc.SiteClass.CART:
                    return self.status
                case bc.SiteClass.RAD_OUTER:
                    return Refused(bc.RefusalAxis.PLUMBING, "corner", 653)
                case bc.SiteClass.RAD_INNER:
                    return Refused(bc.RefusalAxis.PLUMBING, "inner", 511)

    drifted.seal()
    report = pin_disagreements(drifted, tuple(row for row in ROWS if row.method == "mechanism"))
    cells = sorted(line.split(": ")[0] for line in report)
    require(cells == ["mechanism / both 0.3 0.5 @ RAD_OUTER", "mechanism / white 0.5 @ RAD_OUTER"], f"the drifted table's report is {cells}")
