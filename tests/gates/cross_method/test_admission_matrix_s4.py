r"""The admission pin of the CP, MoC and MC tables (wave 2, slice S4; section 4.4): one row per (method, value class, site), pinned EXACTLY.

Spec: ``.claude/plans/wave2_boundary_carve.md`` sections 4.2 (the CP, MoC and MC columns) and 4.4 (the pin), 6.4 S4. The pin's machinery,
its mechanism rows and the SN and diffusion rows are ``tests/gates/transport/test_admission_matrix.py``'s; this file reuses its ``Row``,
``R``, ``F``, ``summary`` and ``CITED_ISSUES`` and holds the three method columns, which moved here on 2026-10-26 (layer DAG ruling 27: as
import statements, the tables' accessors were a transport gate importing three methods; the rows state method claims, so they live with
the cross-method gates). The ids map one to one (``test_the_pin[cp-vacuum-CART_XMIN]`` keeps its id under this path); the coverage and
cited-issue checks are the transport file's, restricted to these three methods.

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

import itertools
from typing import Any

import pytest

from tests.gates._content_identity_helpers import require
from tests.gates.cp import _cp_s4_helpers as CP
from tests.gates.mc import _mc_s4_helpers as MC
from tests.gates.moc import _moc_s4_helpers as MOC
from tests.gates.transport import _wave2_s4_helpers as S4
from tests.gates.transport._admission_helpers import SITES, site, values
from tests.gates.transport._boundary_conditions_helpers import bc, node
from tests.gates.transport.test_admission_matrix import CITED_ISSUES, F, R, Row, summary

pytestmark = pytest.mark.foundation

_MECHANISM = node("test_admission_mechanism.py", "test_the_lattice_of_admissions")


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
    "specular 1.0":                  (("A", (MOC.MOC_DOMAIN_NAME,), (547,)),),
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

ROWS: tuple[Row, ...] = ()
for _method, _columns in (("cp", _CP_COLUMNS), ("moc", _MOC_COLUMNS), ("mc", _MC_COLUMNS)):
    ROWS += tuple(Row(_method, value, where, expected) for value, column in _columns.items() for where, expected in zip(METHOD_SITES[_method], column))

METHOD_VALUES: dict[str, tuple[str, ...]] = {method: _S4_VALUES for method in METHOD_SITES}


def declared_values() -> dict[str, Any]:
    """The value classes of ``_admission_helpers.values()`` plus ``white 1.0`` (the S4 columns' extra class)."""
    return {**values(), "white 1.0": bc.WallCondition.white(1.0)}


def site_of(label: str) -> Any:
    """A site by label: the three of 4.2, or the Cartesian faces CP tells apart."""
    return site(label) if label in SITES else S4.site(label)


TABLES = {"cp": CP.cp_table, "moc": MOC.moc_table, "mc": MC.mc_table}


def _table(method: str) -> Any:
    return TABLES[method]()


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
    """For each of CP, MoC and MC the rows are exactly (the eight value classes of its column) x (its own sites, ``METHOD_SITES``): no
    cell unpinned, none pinned twice."""
    require(set(declared_values()) == set(values()) | {"white 1.0"}, "the S4 value classes changed")
    require({row.method for row in ROWS} == set(METHOD_SITES), f"the methods pinned here are {sorted({row.method for row in ROWS})}")
    for method in METHOD_SITES:
        cells = [(row.value, row.site) for row in ROWS if row.method == method]
        expected = set(itertools.product(METHOD_VALUES[method], METHOD_SITES[method]))
        require(len(cells) == len(set(cells)), f"{method}: a cell is pinned twice")
        require(set(cells) == expected, f"{method}: unpinned {sorted(expected - set(cells))}, extra {sorted(set(cells) - expected)}")


@pytest.mark.rests_on(_MECHANISM)
def test_a_plumbing_refusal_names_a_cited_issue() -> None:
    """Every plumbing refusal the CP, MoC and MC tables RETURN, over the whole grid of their value classes and sites, carries an issue
    number of ``CITED_ISSUES``; so does every plumbing row of this pin. The returned refusals are read from the tables, not from the
    rows, so a realizer that cites a number the pin does not list reds here.

    Reddens on: a CP, MoC or MC realizer citing an issue outside ``CITED_ISSUES``; a plumbing row of the pin without a cited number.
    """
    plumbing = [row for row in ROWS if row.expected[0] == "F" and row.expected[1] == "PLUMBING"]
    require(len(plumbing) >= 10, f"the pin holds {len(plumbing)} plumbing rows")
    for row in plumbing:
        require(isinstance(row.expected[2], int) and row.expected[2] in CITED_ISSUES, f"{row}: a plumbing refusal cites {row.expected[2]!r}, not an issue of CITED_ISSUES")
    declared, returned = declared_values(), []
    for method in sorted(METHOD_SITES):
        table = _table(method)
        for label in METHOD_VALUES[method]:
            for kind in METHOD_SITES[method]:
                admission = table.admit(declared[label], site_of(kind))
                if isinstance(admission, bc.Refused) and admission.axis is bc.RefusalAxis.PLUMBING:
                    returned.append((method, label, kind, admission.issue))
    require(len(returned) >= 10, f"the tables returned {len(returned)} plumbing refusals")
    stray = [entry for entry in returned if entry[3] not in CITED_ISSUES]
    require(not stray, f"plumbing refusals citing an issue outside CITED_ISSUES: {stray}")


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
    for kind in METHOD_SITES[method]:
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
