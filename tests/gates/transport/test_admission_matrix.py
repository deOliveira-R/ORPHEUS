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

First reds, at the base tip ``4ec95bc6``: ``ModuleNotFoundError`` (the first read of the package). Slice S3 appends the diffusion rows
(``_DIFFUSION_COLUMNS``, 21 rows, plus the two derived-row checks); at its base tip ``0757a8d9`` each reads ``DIFFUSION_REALIZERS``,
whose module cannot be imported (``ImportError`` from the removed ``resolve_boundary_conditions``), so the 21 pin rows, the meet
row and the cited-issues row (which ranges over every registered method) are red until the diffusion table lands. Mutations: the
drifted table of ``test_the_pin_has_teeth``; a ``meet`` that returns the greater; a table whose zero row is consulted
for every wall (the vacuum and inflow rows are read by the functor row of the mechanism file).

Slice S4 appends the CP, MoC and MC rows (``_CP_COLUMNS``, 24 rows keyed on ``(site_class, face)`` through CP's slab split ``CART_XMIN`` /
``CART_XMAX``; ``_MOC_COLUMNS``, 8 rows at its one site; ``_MC_COLUMNS``, 16 rows), over the eight value classes of those columns (the seven
of ``_admission_helpers.values()`` plus ``white 1.0``, the weight at which CP's closure is the diffuse chain), each method pinned over its OWN
sites (``METHOD_SITES``; the mechanism, SN and diffusion keep the three of 4.2). At the base tip ``d4da4328`` each of the 48 rows reads
``CP_REALIZERS``, ``MOC_REALIZERS`` or ``MC_REALIZERS``, whose modules cannot be imported (the retired ``BC``), so the 48 pin rows, the
coverage row, the meet rows and the cited-issues row (which ranges over every registered method) are red until the three tables land.
The mixed cells of the MC column follow the MEET law (section 4.2's own rule: the mixed row is the meet of the single-chain rows), which at
MC's radial outer site picks the specular unit's ``Refused(TOPOLOGY, #547)`` over the diffuse unit's plumbing refusal (the lattice's meet
keeps the first refusal in chain order), where the section-4.2 table text reads ``Refused(PLUMBING, #654)``; see NEEDS 2 of the S4 gates report.
"""

from __future__ import annotations

import dataclasses
import itertools
from collections.abc import Callable
from typing import Any

import pytest

from tests.gates._content_identity_helpers import require
from tests.gates.transport import _wave2_s4_helpers as S4
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

# ── Slice S2b: the SN rows (section 4.2, the SN columns), typed by hand and pinned against ``SN_REALIZERS`` ──────────────

#: SN, CART / RAD_OUTER / RAD_INNER per value class. Typed here and nowhere else (the mechanism's copy above is the
#: test-local table's, not SN's): a change to an SN admission edits THIS literal. The hollow inner face admits a
#: specular chain of weight 1.0 only (the cavity the closure computes, #511); the radial outer face refuses the diffuse
#: chain and every source (the corner action, #653); the zero network is refused on the inner face (#511).
_SN_COLUMNS: dict[str, tuple[tuple, tuple, tuple]] = {
    "vacuum":                        (R(), R(), F("PLUMBING", 511)),
    "specular 1.0":                  (R(), R(), R()),
    "specular 0.5":                  (R(), R(), F("PLUMBING", 511)),
    "white 0.5":                     (R(), F("PLUMBING", 653), F("PLUMBING", 511)),
    "both 0.3 0.5":                  (R(), F("PLUMBING", 653), F("PLUMBING", 511)),
    "inflow 1.0":                    (R(), F("PLUMBING", 653), F("PLUMBING", 511)),
    "specular 0.5 over inflow 1.0":  (R(), F("PLUMBING", 653), F("PLUMBING", 511)),
}

ROWS += tuple(Row("sn", value, where, expected) for value, column in _SN_COLUMNS.items() for where, expected in zip(SITES, column))


# ── Slice S3: the diffusion rows (section 4.2, the diffusion column), typed by hand and pinned against ``DIFFUSION_REALIZERS`` ──

#: Diffusion, CART / RAD_OUTER / RAD_INNER per value class: the SAME cell at the Cartesian and the radial outer site (a per-face
#: scalar that reads nothing of the face). The hollow inner face has no trace slot, so the stencil computes zero net current there
#: whatever is declared: only ``specular 1.0`` is admitted at RAD_INNER and every other class is ``Refused(PLUMBING, #511)``, as in
#: SN's column (``[REVISED 2026-10-23]`` at the S3 integration: the section 4.2 draft admitted diffusion at every site, which dropped
#: the inner condition silently, ``[M]`` bit-identical k for 3 of 3 inner conditions). The empty network is tier 1 (the deciding line of the measurement is ``return response.amplitude``;
#: Marshak is the method's closure, not the vacuum law's admission); a chain is tier 2 (``Realized(False)``: diffusion annihilates
#: the shape) and so is a mixed wall (the meet of two tier-2 units, ``A = s + d``); a nonzero source, alone or over a network, is the
#: plumbing refusal of #656 (the affine channel). Typed here and nowhere else: a change to a diffusion admission edits THIS literal.
#: The method-only units (``ZeroFluxCondition``, the decks) are not value classes of the shared family; their rows are
#: ``tests/gates/diffusion/test_diffusion_realizer_table.py``'s.
_DIFFUSION_COLUMNS: dict[str, tuple[tuple, tuple, tuple]] = {
    "vacuum":                        (R(), R(), F("PLUMBING", 511)),
    "specular 1.0":                  (R(False), R(False), R(False)),
    "specular 0.5":                  (R(False), R(False), F("PLUMBING", 511)),
    "white 0.5":                     (R(False), R(False), F("PLUMBING", 511)),
    "both 0.3 0.5":                  (R(False), R(False), F("PLUMBING", 511)),
    "inflow 1.0":                    (F("PLUMBING", 656), F("PLUMBING", 656), F("PLUMBING", 511)),
    "specular 0.5 over inflow 1.0":  (F("PLUMBING", 656), F("PLUMBING", 656), F("PLUMBING", 511)),
}

ROWS += tuple(Row("diffusion", value, where, expected) for value, column in _DIFFUSION_COLUMNS.items() for where, expected in zip(SITES, column))


def _diffusion_table() -> Any:
    """``DIFFUSION_REALIZERS`` read at call time: at the base tip of slice S3 ``import orpheus.diffusion`` raises (its mesh module imports
    the resolver S2b removed) and the realizer module exports ``DiffusionBoundaryRealizer``."""
    from tests.gates.diffusion._wave2_helpers import table

    return table()


def _sn_table() -> Any:
    """``SN_REALIZERS`` read at call time: at the base tip of slice S2b the module still exports ``SNBoundaryRealizer``."""
    from tests.gates.sn.architecture._realized_boundary_helpers import sn_table

    return sn_table()




# ── Slice S4: the CP, MoC and MC rows (section 4.2, the CP, MoC and MC columns), typed by hand and pinned against ``CP_REALIZERS``,
# ``MOC_REALIZERS`` and ``MC_REALIZERS`` ──────────────────────────────────────────────────────────────────────────────────────────

#: The value classes of the S4 columns: the seven of ``values()`` and ``white 1.0``.
_S4_VALUES = ("vacuum", "specular 1.0", "specular 0.5", "white 0.5", "white 1.0", "both 0.3 0.5", "inflow 1.0", "specular 0.5 over inflow 1.0")

#: The sites each method is pinned over. CP's table is keyed on ``(site_class, face)``: a ``CART`` site is split by its face (``xmin`` is the
#: forced mirror), so the Cartesian column is two. MoC reads a solid cylinder only (its slab and hollow body are the guard's, #514).
METHOD_SITES: dict[str, tuple[str, ...]] = {
    "cp": ("CART_XMIN", "CART_XMAX", "RAD_OUTER"),
    "moc": ("RAD_OUTER",),
    "mc": ("CART", "RAD_OUTER"),
}

#: CP: CART ``xmin`` | CART ``xmax`` | RAD_OUTER. The left face is a forced mirror, so only a deck is admitted there (#513: every wall is refused);
#: elsewhere the empty network and the diffuse chain at weight 1.0 are realized (D4: ``k`` omits leakage, #546, carried by the strict xfails; the
#: albedo below one is refused, #654, today a silent drop) and everything else is the old registry's gap (#654).
_CP_COLUMNS: dict[str, tuple[tuple, ...]] = {
    "vacuum":                        (F("PLUMBING", 513), R(), R()),
    "specular 1.0":                  (F("PLUMBING", 513), F("PLUMBING", 654), F("PLUMBING", 654)),
    "specular 0.5":                  (F("PLUMBING", 513), F("PLUMBING", 654), F("PLUMBING", 654)),
    "white 0.5":                     (F("PLUMBING", 513), F("PLUMBING", 654), F("PLUMBING", 654)),
    "white 1.0":                     (F("PLUMBING", 513), R(), R()),
    "both 0.3 0.5":                  (F("PLUMBING", 513), F("PLUMBING", 654), F("PLUMBING", 654)),
    "inflow 1.0":                    (F("PLUMBING", 513), F("PLUMBING", 654), F("PLUMBING", 654)),
    "specular 0.5 over inflow 1.0":  (F("PLUMBING", 513), F("PLUMBING", 654), F("PLUMBING", 654)),
}

#: MoC, RAD_OUTER (a solid cylinder): the specular chain at weight 1.0 is a named DOMAIN approximation (the equal-area square cell, #547); every
#: other wall is the old registry's gap (#654), the diffuse chain included at every weight.
_MOC_COLUMNS: dict[str, tuple[tuple, ...]] = {
    "vacuum":                        (F("PLUMBING", 654),),
    "specular 1.0":                  (("A", (S4.MOC_DOMAIN_NAME,), (547,)),),
    "specular 0.5":                  (F("PLUMBING", 654),),
    "white 0.5":                     (F("PLUMBING", 654),),
    "white 1.0":                     (F("PLUMBING", 654),),
    "both 0.3 0.5":                  (F("PLUMBING", 654),),
    "inflow 1.0":                    (F("PLUMBING", 654),),
    "specular 0.5 over inflow 1.0":  (F("PLUMBING", 654),),
}

#: MC: CART | RAD_OUTER (a solid cylinder). The specular chain on the cylinder is the torus refusal, for every weight (#547, the user's ruling);
#: at a Cartesian site and for every other wall the old registry held nothing (#654). The mixed cells of the radial column are the MEET of the
#: single-chain cells (the specular unit's TOPOLOGY refusal comes first in chain order), see the module docstring.
_MC_COLUMNS: dict[str, tuple[tuple, ...]] = {
    "vacuum":                        (F("PLUMBING", 654), F("PLUMBING", 654)),
    "specular 1.0":                  (F("PLUMBING", 654), F("TOPOLOGY", 547)),
    "specular 0.5":                  (F("PLUMBING", 654), F("TOPOLOGY", 547)),
    "white 0.5":                     (F("PLUMBING", 654), F("PLUMBING", 654)),
    "white 1.0":                     (F("PLUMBING", 654), F("PLUMBING", 654)),
    "both 0.3 0.5":                  (F("PLUMBING", 654), F("TOPOLOGY", 547)),
    "inflow 1.0":                    (F("PLUMBING", 654), F("PLUMBING", 654)),
    "specular 0.5 over inflow 1.0":  (F("PLUMBING", 654), F("TOPOLOGY", 547)),
}

for _method, _columns in (("cp", _CP_COLUMNS), ("moc", _MOC_COLUMNS), ("mc", _MC_COLUMNS)):
    ROWS += tuple(Row(_method, value, where, expected) for value, column in _columns.items() for where, expected in zip(METHOD_SITES[_method], column))

METHOD_VALUES: dict[str, tuple[str, ...]] = {method: _S4_VALUES for method in METHOD_SITES}


def declared_values() -> dict[str, Any]:
    """The value classes of ``_admission_helpers.values()`` plus ``white 1.0`` (the S4 columns' extra class)."""
    return {**values(), "white 1.0": bc.WallCondition.white(1.0)}


def site_of(label: str) -> Any:
    """A site by label: the three of 4.2, or the Cartesian faces CP tells apart."""
    return site(label) if label in SITES else S4.site(label)


def method_values(method: str) -> tuple[str, ...]:
    return METHOD_VALUES.get(method, tuple(values()))


def method_sites(method: str) -> tuple[str, ...]:
    return METHOD_SITES.get(method, SITES)


def _cp_table() -> Any:
    return S4.cp_table()


def _moc_table() -> Any:
    return S4.moc_table()


def _mc_table() -> Any:
    return S4.mc_table()

TABLES: dict[str, Callable[[], Any]] = {"mechanism": lambda: make_table().table, "sn": _sn_table, "diffusion": _diffusion_table, "cp": _cp_table, "moc": _moc_table, "mc": _mc_table}


def _table(method: str) -> Any:
    return TABLES[method]()


def pin_disagreements(table: Any, rows: tuple[Row, ...]) -> list[str]:
    """The rows whose pinned summary differs from ``table.admit``'s (an empty list is a pinned table)."""
    declared = declared_values()
    out = []
    for row in rows:
        got = summary(table.admit(declared[row.value], site_of(row.site)))
        if got != row.expected:
            out.append(f"{row.method} / {row.value} @ {row.site}: admit gives {got}, the pin says {row.expected}")
    return out


def _ids(rows: tuple[Row, ...]) -> list[str]:
    return [f"{r.method}-{r.value}-{r.site}".replace(" ", "_") for r in rows]


@pytest.mark.parametrize("row", ROWS, ids=_ids(ROWS))
@pytest.mark.rests_on(_MECHANISM)
def test_the_pin(row: Row) -> None:
    """The table's ``admit`` equals the row (summary: tier, names/issues or axis/issue)."""
    got = summary(_table(row.method).admit(declared_values()[row.value], site_of(row.site)))
    require(got == row.expected, f"{row.method} / {row.value} @ {row.site}: admit gives {got}, the pin says {row.expected}")


@pytest.mark.rests_on(_MECHANISM)
def test_the_pin_covers_every_cell() -> None:
    """For each method the rows are exactly (its value classes) x (its sites): no cell unpinned, none pinned twice. The mechanism, SN and
    diffusion range over the seven value classes and the three sites of 4.2; CP, MoC and MC over the eight value classes of their columns
    and their own sites (``METHOD_SITES``)."""
    declared = set(values())
    require(len(declared) == 7 and set(SITES) == {"CART", "RAD_OUTER", "RAD_INNER"}, "the value classes or sites changed")
    require(set(declared_values()) == declared | {"white 1.0"}, "the S4 value classes changed")
    for method in {row.method for row in ROWS}:
        cells = [(row.value, row.site) for row in ROWS if row.method == method]
        expected = set(itertools.product(method_values(method), method_sites(method)))
        require(len(cells) == len(set(cells)), f"{method}: a cell is pinned twice")
        require(set(cells) == expected, f"{method}: unpinned {sorted(expected - set(cells))}, extra {sorted(set(cells) - expected)}")


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

    declared, returned = declared_values(), []
    for method in sorted({row.method for row in ROWS}):
        table = _table(method)
        for label in method_values(method):
            value = declared[label]
            for kind in method_sites(method):
                admission = table.admit(value, site_of(kind))
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
    mechanism_rows = tuple(row for row in ROWS if row.method == "mechanism")  # S3: the SN columns equalled the mechanism's, the diffusion column does not
    require(pin_disagreements(_table("mechanism"), mechanism_rows) == [], f"the honest table disagrees with its pin: {pin_disagreements(_table('mechanism'), mechanism_rows)}")

    drifted = bc.RealizerTable("drifted", state_cone=bc.StateCone.NONNEGATIVE_FLUX, unregistered_issue=654)
    R_, Refused = bc.Realized, bc.Refused

    class Spec(bc.ResponseRealizer, table=drifted, realizes=bc.SpecularChain):
        status = R_(True)

        def admit(self, unit, weight, site, space=None):
            if site.site_class is bc.SiteClass.RAD_INNER:
                return R_(True) if weight == 1.0 else Refused(bc.RefusalAxis.PLUMBING, "inner", 511)
            return self.status

        def realize(self, unit, weight, site, space):
            raise NotImplementedError

    class Diff(bc.ResponseRealizer, table=drifted, realizes=bc.DiffuseChain):
        status = R_(True)

        def admit(self, unit, weight, site, space=None):
            return Refused(bc.RefusalAxis.PLUMBING, "inner", 511) if site.site_class is bc.SiteClass.RAD_INNER else self.status  # RAD_OUTER DRIFTED to admitted

        def realize(self, unit, weight, site, space):
            raise NotImplementedError

    class Zero(bc.ResponseRealizer, table=drifted, realizes=bc.ZeroNetwork):
        status = R_(True)

        def admit(self, unit, weight, site, space=None):
            return Refused(bc.RefusalAxis.PLUMBING, "inner", 511) if site.site_class is bc.SiteClass.RAD_INNER else self.status

        def realize(self, unit, weight, site, space):
            raise NotImplementedError

    class Src(bc.SourceRealizer, table=drifted):
        status = R_(True)

        def admit(self, terms, site, space=None):
            match site.site_class:
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


# ── Slice S2b: what the SN table adds to the pin ─────────────────────────────────────────────────────────────────────────


@pytest.mark.rests_on(_MECHANISM)
def test_sn_both_chains_and_a_response_over_a_source_are_the_meet() -> None:
    """The SN mixed rows are DERIVED from the single-chain and source rows of the same table, at every site: a change
    to the diffuse row edits the "both" row and the "specular over inflow" row (section 4.2's own rule). The pieces
    are read from ``SN_REALIZERS``, not from the pin.

    Reddens on: an ``SNDiffuseChain.admit`` that admits a radial outer site while the wall's admission is cached from the
    specular row; a ``meet`` that returns the greater (the S2a gate of the mechanism shares the operator).
    """
    table, W = _table("sn"), bc.WallCondition
    for kind in SITES:
        where = site(kind)
        spec, white = table.admit(W.specular(0.3), where), table.admit(W.white(0.5), where)
        require(summary(table.admit(W.from_amplitudes(0.3, 0.5), where)) == summary(bc.Admission.meet(spec, white)), f"SN: both chains @ {kind} is not the meet of the two single-chain rows")
        response, source = table.admit(W.specular(0.5), where), table.admit(bc.InflowTerms(1.0, ()), where)
        require(summary(table.admit(W.specular(0.5) + W.uniform(1.0), where)) == summary(bc.Admission.meet(response, source)), f"SN: a response over a source @ {kind} is not the meet")


@pytest.mark.rests_on(_MECHANISM)
def test_sn_plumbing_refusals_cite_the_sn_issues_and_the_catalogue_is_an_upper_bound() -> None:
    """Every plumbing refusal the SN table returns over the whole grid cites #511 or #653 (the unregistered-unit door
    cites #659 and is read by ``test_the_unregistered_door_names_the_state_cone_...``); and the catalogue is an upper
    bound: SN's catalogue is ``Realized(True)`` for every unit, so no admitted cell may be anything but a faithful
    ``Realized``. The population is printed (X2): 7 value classes times 3 sites.

    Reddens on: an SN realizer whose ``admit`` cites an issue outside this list; a ``status`` below what ``admit``
    returns (the mutation of S2a gate 13 (f), applied to the real table).
    """
    table, declared = _table("sn"), values()
    seen, cited = 0, set()
    for label, value in declared.items():
        for kind in SITES:
            admission = table.admit(value, site(kind))
            seen += 1
            if isinstance(admission, bc.Refused) and admission.axis is bc.RefusalAxis.PLUMBING:
                cited.add(admission.issue)
            if isinstance(admission, bc.Realized):
                require(admission.faithful is True, f"{label} @ {kind}: SN realizes every admitted cell faithfully (its catalogue is Realized(True))")
    require(seen == 21, f"the SN grid is {seen} cells, not 7 x 3")
    require(cited == {511, 653}, f"the SN plumbing refusals cite {sorted(cited)}, not {{511, 653}}")


# ── Slice S3: what the diffusion table adds to the pin ─────────────────────────────────────────────────────────────────


@pytest.mark.rests_on(_MECHANISM)
def test_diffusion_both_chains_and_a_response_over_a_source_are_the_meet() -> None:
    """The diffusion mixed rows are DERIVED from the single-chain and source rows of the same table, at every site: a change to
    either chain's row edits the "both" row (``Realized(False)``: the meet of two tier-2 units) and the "specular over inflow" row
    (the meet with the plumbing refusal of #656). The pieces are read from ``DIFFUSION_REALIZERS``, not from the pin.

    Reddens on: a ``DiffusionSpecularChain.admit`` that returns tier 1 while the wall's admission is cached from the diffuse row; a
    ``meet`` that returns the greater (the S2a gate of the mechanism shares the operator). Red at the base tip: ``import
    orpheus.diffusion`` raises.
    """
    table, W = _table("diffusion"), bc.WallCondition
    for kind in SITES:
        where = site(kind)
        spec, white = table.admit(W.specular(0.3), where), table.admit(W.white(0.5), where)
        require(summary(table.admit(W.from_amplitudes(0.3, 0.5), where)) == summary(bc.Admission.meet(spec, white)), f"diffusion: both chains @ {kind} is not the meet of the two single-chain rows")
        response, source = table.admit(W.specular(0.5), where), table.admit(bc.InflowTerms(1.0, ()), where)
        require(summary(table.admit(W.specular(0.5) + W.uniform(1.0), where)) == summary(bc.Admission.meet(response, source)), f"diffusion: a response over a source @ {kind} is not the meet")


# ── Slice S4: what the CP, MoC and MC tables add to the pin ───────────────────────────────────────────────────────────────


@pytest.mark.parametrize("method", ["cp", "moc", "mc"])
@pytest.mark.rests_on(_MECHANISM)
def test_s4_both_chains_and_a_response_over_a_source_are_the_meet(method: str) -> None:
    """The mixed rows of CP, MoC and MC are DERIVED from the single-chain and source rows of the same table, at every site the method is
    pinned over: a change to either chain's row edits the "both" row, and a change to the diffuse or source row edits the "specular over
    inflow" row (section 4.2's own rule). The pieces are read from the table, not from the pin. The instance of the rule that is not the
    section text's reading, MC's radial outer site (the meet keeps the specular unit's TOPOLOGY refusal, which comes first in chain order),
    is the row the pin edits if the ruling goes the other way.

    Reddens on: a realizer pair whose wall admission is cached from one row; a ``meet`` that returns the greater. Red at the base tip:
    ``ImportError`` at the table accessor.
    """
    table, W = _table(method), bc.WallCondition
    for kind in method_sites(method):
        where = site_of(kind)
        spec, white = table.admit(W.specular(0.3), where), table.admit(W.white(0.5), where)
        require(summary(table.admit(W.from_amplitudes(0.3, 0.5), where)) == summary(bc.Admission.meet(spec, white)), f"{method}: both chains @ {kind} is not the meet of the two single-chain rows")
        response, source = table.admit(W.specular(0.5), where), table.admit(bc.InflowTerms(1.0, ()), where)
        require(summary(table.admit(W.specular(0.5) + W.uniform(1.0), where)) == summary(bc.Admission.meet(response, source)), f"{method}: a response over a source @ {kind} is not the meet")


@pytest.mark.parametrize("method", ["cp", "moc", "mc"])
@pytest.mark.rests_on(_MECHANISM)
def test_s4_tables_hold_only_what_the_old_registries_held(method: str) -> None:
    """CP, MoC and MC are NOT improved on their own terms (``CLAUDE.md``, "The direction of development"): each table realizes ONLY what
    its string registry realized (CP: ``white`` and ``vacuum``; MoC: ``reflective``; MC: ``periodic``) and states the rest as the
    plumbing refusal of #654 (or the theory axis of #547 for the torus). So the faithful cells of a table are few and enumerable: CP
    realizes exactly the empty network and ``white 1.0`` at ``xmax`` and ``RAD_OUTER`` (4 of its 24 cells), MoC realizes no cell
    faithfully and approximates one, MC realizes none of the shared value classes (its one realized unit is the translation deck). A new
    faithful cell is a change of the table's scope and owes an edit here.

    Reddens on: a table that widened its coverage (an admission that returns ``Realized`` for a class the registry did not hold).
    """
    realized = sorted((r.value, r.site) for r in ROWS if r.method == method and r.expected[0] == "R")
    expected = {
        "cp": [("vacuum", "CART_XMAX"), ("vacuum", "RAD_OUTER"), ("white 1.0", "CART_XMAX"), ("white 1.0", "RAD_OUTER")],
        "moc": [],
        "mc": [],
    }[method]
    require(realized == sorted(expected), f"{method}: the pin realizes {realized}")
    approximated = [(r.value, r.site) for r in ROWS if r.method == method and r.expected[0] == "A"]
    require(approximated == ([("specular 1.0", "RAD_OUTER")] if method == "moc" else []), f"{method}: approximated cells {approximated}")
