r"""The references import no production response, and the readers that did are gone (wave 2, S5; plan 2.9 and gate 6 of 6.5).

Specified by the test-architect (2026-10-23, ``.claude/plans/wave2_boundary_carve.md`` 6.5). A reference states its boundary in
its OWN closed values (``WallAmplitudes(specular, diffuse)``, a mirror and a wrap read from geometry's decks, plain data) and
imports no production response: the user's ruling "references closed, production versatile" (references are insulated from
production numerics), and the interface between references and production is plain data, crossed by ONE test-side loop
(:mod:`tests.gates._reference_pairing`). Two claims, three rows.

* **The fence** (``test_the_references_import_nothing_from_transport_or_the_retiring_boundary_package``). An AST scan of every
  module under ``orpheus/derivations`` (relative imports RESOLVED to their absolute name against the importing file's package,
  function-body and ``TYPE_CHECKING`` imports included) finds nothing under ``orpheus.transport`` or
  ``orpheus.geometry.boundary``. The insulation ledger ``tests/gates/test_layer_imports.py::
  test_closed_references_import_only_the_numerics_interface`` already forbids the transport package and every method package
  and stays green; this row's increment is ``orpheus.geometry.boundary``, which is NOT forbidden by that ledger (the geometry
  package is an allowed import of a reference) and which the carve retires: the two modules that still read the retiring
  law family are ``derivations/common/reference_body.py`` and ``derivations/continuous/characteristic/walls.py``, so the
  population this row retires is exactly those two. The instrument is :mod:`tests.gates._import_scan`, whose seven planted
  controls (one per import spelling) and whose negative leg (a clean source and a sibling that shares leading letters) the row
  runs FIRST, and whose population (files scanned, imported names read) the failure message states
  (``instrument-doctrine`` X1, X2).
* **The retired readers** (``test_the_retired_law_readers_are_gone``). ``TAG_REGISTRY``, ``_law_of``, ``_wall_of``, ``_parameters``
  (``walls.py``), ``specular_albedo`` (``reference_body.py``) and the ``BC.to_alpha`` method have no definition and no
  reference under ``orpheus/derivations`` (the census tree; ``to_alpha`` itself is still defined in ``orpheus/geometry/boundary/_tag.py``
  until WF2.4 deletes that package, so the claim is about the reference tree only) (identifier census over ``ast.Name``, ``ast.Attribute``, import aliases
  and function/class definitions, so a docstring mention is not a use). The census's control is a known live name of each file
  (``Walls``; ``specular_albedos``). A retirement that left a reader behind would keep the law family importable from a reference.
* **The carve is what moves it** (``test_the_reference_walls_module_names_no_production_type``): the ``walls.py`` module reads
  ``WallAmplitudes`` and geometry's decks (``MirrorDeck``, ``TranslationDeck`` from ``orpheus.geometry.deck`` or the geometry
  package) and no name of the retiring law family (``BoundaryTraceLaw``, ``AlbedoBoundary``, ``ReflectiveBoundary``,
  ``WhiteBoundary``, ``PeriodicBoundary``, ``VacuumInflow``, ``ZeroFluxBoundary``, ``PrescribedInflow``, ``BC``).

Declared red at the base tip (``de354496``), the input that reddens each (X1): the fence reads ``reference_body.py`` and
``walls.py`` importing ``orpheus.geometry.boundary`` (2 files of the 145 scanned); the retired-readers census finds the five names
defined in the two files; the walls module names seven of the nine law-family names (``PrescribedInflow`` and ``ZeroFluxBoundary`` it does not import). Named mutations the reviewers apply to the carved tree:
(1) re-add ``from orpheus.geometry.boundary import BC`` to any ``derivations`` module (even inside a function) reddens the fence;
(2) re-add ``from ...geometry.boundary import X`` (relative) reddens it too, which is why the scan resolves relative imports
instead of skipping them; (3) restore ``def _wall_of`` or ``TAG_REGISTRY`` in ``walls.py`` reddens the census; (4) re-add
``from orpheus.transport.boundary_conditions import WallCondition`` reddens the fence (the transport half; the insulation ledger
reddens beside it, and both being red is the expected reading: they are two instruments of one fence).

Classification of the surviving gates this change touches (``test-architect`` section 4, each owed in the commit body):
``test_characteristic_walls.py::test_the_registry_admits_exactly_the_kinds_this_file_twins``,
``test_a_tag_reads_the_same_walls_as_the_law_it_names``, ``test_the_white_tag_without_an_albedo_is_a_full_return``,
``test_the_registry_hands_the_white_law_the_walls_outward_sign`` and the ``unadmitted_tag_*``, ``*_with_a_parameter``,
``partial_without_albedo`` rows are DEAD (their subject, the string-tag registry, no longer exists; delete, never repair by
passing a new argument); the ``albedo_with_unstated_shape``, ``zero_flux`` and ``inflow_source`` rows are DEAD BY TYPE (the
specification has no field that can state them: ``tests/gates/specification/test_exterior_on_the_specification.py``
asserts the absence); ``wrap_not_paired*`` and ``wrap_on_*`` are MOVED to geometry (``test_curved_wraps_are_deck_errors.py``);
``amplitude_above_one`` and ``amplitude_negative`` are RE-POSED under R-8b (``ValueError`` from ``WallAmplitudes``, fragment kept);
``test_a_directly_built_walls_refuses_what_walls_of_refuses`` SURVIVES (``Walls.__post_init__`` keeps its radial-wrap and
unpaired-wrap checks for a directly built ``Walls``, so those guards are live, not DEAD) with its last line re-posed.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

from tests.gates import _import_scan
from tests.gates._content_identity_helpers import require

pytestmark = pytest.mark.foundation

_ROOT = Path(__file__).resolve().parents[3]
_DERIVATIONS = _ROOT / "orpheus" / "derivations"
_WALLS = _DERIVATIONS / "continuous" / "characteristic" / "walls.py"
_REFERENCE_BODY = _DERIVATIONS / "common" / "reference_body.py"
_FENCED = ("orpheus.transport", "orpheus.geometry.boundary")
_LAYER = "tests/gates/test_layer_imports.py::test_closed_references_import_only_the_numerics_interface"
_HERE = "tests/gates/derivations/test_references_import_no_production.py::"


def _package_of(path: Path) -> str:
    """The dotted package a relative import in ``path`` is read from."""
    rel = path.relative_to(_ROOT).with_suffix("")
    return ".".join(rel.parts[:-1])  # a module and an ``__init__`` alike: the directory is the package


def test_the_import_scan_finds_each_spelling_and_resolves_a_relative_import() -> None:
    """The instrument of the fence reads positive on each spelling, negative on clean sources (``instrument-doctrine`` X1)."""
    for source, prefix in _import_scan.CONTROLS:
        found = _import_scan.forbidden_imports(source, [prefix])
        require(bool(found), f"control: {source!r} not found under {prefix}")
    relative = _import_scan.forbidden_imports("from ...geometry.boundary import X", ["orpheus.geometry.boundary"], "orpheus.derivations.common")
    require(bool(relative), "control: a relative import reaching orpheus.geometry.boundary is not resolved")
    require(_import_scan.forbidden_imports("from . import y", list(_FENCED), "orpheus.derivations.common") == [], "control: a sibling import is flagged")
    require(_import_scan.forbidden_imports("import math\nfrom orpheus.transport_x import y", list(_FENCED)) == [], "control: a lookalike prefix is flagged")
    require(_package_of(_WALLS) == "orpheus.derivations.continuous.characteristic", f"control: package of walls.py is {_package_of(_WALLS)}")
    require(_package_of(_DERIVATIONS / "__init__.py") == "orpheus.derivations", "control: package of an __init__ is its own directory")


@pytest.mark.rests_on(_HERE + "test_the_import_scan_finds_each_spelling_and_resolves_a_relative_import", _LAYER)
def test_the_references_import_nothing_from_transport_or_the_retiring_boundary_package() -> None:
    """No module under ``orpheus/derivations`` imports ``orpheus.transport`` or ``orpheus.geometry.boundary`` in any spelling.

    Red at base: ``reference_body.py`` and ``walls.py`` import ``orpheus.geometry.boundary`` (2 of 145 files). The message lists
    each offender with the names it imports, and the scan's population."""
    files = sorted(_DERIVATIONS.rglob("*.py"))
    require(len(files) > 100, f"control: the scan found {len(files)} files under orpheus/derivations")
    require(_WALLS in files and _REFERENCE_BODY in files, "control: the two known readers are in the scanned population")
    offenders: dict[str, list[str]] = {}
    names_read = 0
    for path in files:
        source = path.read_text()
        names_read += len(_import_scan.imported_names(source, _package_of(path)))
        found = _import_scan.forbidden_imports(source, _FENCED, _package_of(path))
        if found:
            offenders[str(path.relative_to(_ROOT))] = sorted(set(found))
    require(names_read > 1000, f"control: the scan read only {names_read} imported names")
    require(
        not offenders,
        f"{len(offenders)} of {len(files)} files under orpheus/derivations import {list(_FENCED)} ({names_read} imported names read): "
        + "; ".join(f"{f}: {names[:4]}" for f, names in sorted(offenders.items())),
    )


def _defined_and_used(source: str) -> set[str]:
    """Every identifier a source defines (``def``, ``class``, assignment) or uses (``Name``, ``Attribute``, import alias)."""
    names: set[str] = set()
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            names.add(node.name)
        elif isinstance(node, ast.Name):
            names.add(node.id)
        elif isinstance(node, ast.Attribute):
            names.add(node.attr)
        elif isinstance(node, (ast.Import, ast.ImportFrom)):
            names |= {alias.asname or alias.name.split(".")[-1] for alias in node.names}
    return names


_RETIRED = ("TAG_REGISTRY", "_law_of", "_wall_of", "_parameters", "specular_albedo", "to_alpha")


@pytest.mark.rests_on(_HERE + "test_the_references_import_nothing_from_transport_or_the_retiring_boundary_package")
def test_the_retired_law_readers_are_gone() -> None:
    """None of the six retired names is defined or used anywhere under ``orpheus/derivations`` (an identifier census, docstrings excluded).

    Positive control: the census of ``walls.py`` finds ``Walls`` and that of ``reference_body.py`` finds ``specular_albedos``
    (the plural, a survivor with a changed signature: a substring match would confuse it with the retired singular, an exact
    identifier match does not). Red at base: ``TAG_REGISTRY``, ``_law_of``, ``_wall_of`` and ``_parameters`` are defined in
    ``walls.py``; ``specular_albedo`` in ``reference_body.py`` and ``to_alpha`` used by it."""
    require("Walls" in _defined_and_used(_WALLS.read_text()), "control: the census does not find Walls in walls.py")
    require("specular_albedos" in _defined_and_used(_REFERENCE_BODY.read_text()), "control: the census does not find specular_albedos")
    require("specular_albedo" not in _defined_and_used("def specular_albedos(): pass"), "control: the plural is read as the retired singular")
    found: dict[str, list[str]] = {}
    files = sorted(_DERIVATIONS.rglob("*.py"))
    for path in files:
        hits = [n for n in _RETIRED if n in _defined_and_used(path.read_text())]
        if hits:
            found[str(path.relative_to(_ROOT))] = hits
    require(not found, f"retired readers still present ({len(files)} files scanned): {found}")


_LAW_FAMILY = (
    "BoundaryTraceLaw", "AlbedoBoundary", "ReflectiveBoundary", "WhiteBoundary", "PeriodicBoundary",
    "VacuumInflow", "ZeroFluxBoundary", "PrescribedInflow", "BC",
)


@pytest.mark.rests_on(_HERE + "test_the_retired_law_readers_are_gone")
def test_the_reference_walls_module_names_no_production_type() -> None:
    """``walls.py`` reads ``WallAmplitudes``, ``TranslationDeck`` (the partner) and the one role reader ``amplitudes_at`` (which names both decks, in ``reference_body.py``), and names none of the retiring law family's nine names.

    Control: the same census over a planted source with ``BC`` finds it. Red at base: seven of the nine are named (``BoundaryTraceLaw``,
    ``AlbedoBoundary``, ``ReflectiveBoundary``, ``WhiteBoundary``, ``PeriodicBoundary``, ``VacuumInflow``, ``BC``). ``WallAmplitudes`` is named after the carve (``Wall`` composes it), so its presence is part of
    the row: a ``walls.py`` that dropped the law family and stated its amplitudes as bare floats would pass the first half
    only."""
    require("BC" in _defined_and_used("from orpheus.geometry.boundary import BC"), "control: the census misses a planted BC")
    names = _defined_and_used(_WALLS.read_text())
    named = [n for n in _LAW_FAMILY if n in names]
    require(not named, f"walls.py still names the retiring law family: {named}")
    require("WallAmplitudes" in names, "walls.py does not name WallAmplitudes: Wall must compose it (decision 54)")
    require({"TranslationDeck", "amplitudes_at"} <= names, "walls.py reads no deck of its own (TranslationDeck, for the partner) or no role reader (amplitudes_at)")
    # the role's amplitudes are read in ONE place (the amend of review finding 7): both decks are named there, and no law-family name
    shared = _defined_and_used(_REFERENCE_BODY.read_text())
    require({"MirrorDeck", "TranslationDeck"} <= shared, "reference_body.py reads no geometry deck (MirrorDeck, TranslationDeck)")
    require(not [n for n in _LAW_FAMILY if n in shared], f"reference_body.py names the retiring law family: {[n for n in _LAW_FAMILY if n in shared]}")
