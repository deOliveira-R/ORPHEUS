r"""Step 4 of #405 P3: the clients — the characteristic reference's reading and its solve; the exact medium is not one.

The step's first clients, the trajectory resolvent's two multi-region solves and its reading (rows M4.1b-M4.8), were
deleted with the family in step (e2) of the characteristic-reference campaign; their successors on the characteristic
reference are the rows below "the characteristic reading" (and the synthetic child-memo rows of
``test_traced_memo_process.py``).

Spec ``.claude/plans/reference_p3_spec.md`` §1.4, gates M4.1–M4.10. Client names are resolved through
``_traced_memo_api.CLIENTS``. The real-tree mutation witnesses (M4.4) run on a COPY of ``orpheus/``'s Python
files in ``tmp_path``, in a fresh interpreter whose ``sys.path`` starts at the copy; no tracked file is edited.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
import textwrap
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from . import _traced_memo_api as api

pytestmark = pytest.mark.foundation

REPO = Path(__file__).resolve().parents[3]


def _exact_specification():
    from orpheus.data.cells import CellCoefficient, Channel
    from orpheus.derivations.common.xs_library import get_mixture
    from orpheus.numerics.question import Eigen
    from orpheus.specification.specification import InfiniteMediumSpecification

    return InfiniteMediumSpecification(0, get_mixture("A", "2g"), Eigen(CellCoefficient.every(Channel.FISSION_EMISSION)))


def test_m4_1_the_clients_are_traced_memos():
    """M4.1 (a): every client of ``api.CLIENTS`` (the characteristic reading and its solve) is bound as a traced
    memo where every caller reaches it (the class attribute ``ReferenceSolution.read`` and ``answer`` call), and its
    function id is the one ``api.spelled_function_id`` spells."""
    memo_type = api.name("TracedMemo")
    for client in api.CLIENTS:
        assert isinstance(api.client(client), memo_type), f"{client} is not a traced memo"
        assert api.function_id(client) == api.spelled_function_id(client), client


# ── M4.4: the real-tree witnesses, on a copy ─────────────────────────────────────

def _neutral_edit(path: Path, anchor: str, method: str | None) -> None:
    """A behaviour-neutral AST edit: ``_ = None`` as the first statement of the def at ``anchor`` (or of the
    method ``method`` of the class at ``anchor``)."""
    text = path.read_text()
    if text.count(anchor) != 1:
        raise AssertionError(f"{anchor!r} occurs {text.count(anchor)} times in {path}")
    start = text.index(anchor)
    if method is not None:
        start = text.index(f"def {method}(", start)
    body = text.index(":\n", text.index(")", start)) + 2
    indent = len(text[body:]) - len(text[body:].lstrip(" "))
    if text[body:].lstrip(" ").startswith(('"""', "r\"\"\"")):  # step past the docstring
        close = text.index('"""', text.index('"""', body) + 3) + 4
        body = close
    new = text[:body] + " " * indent + "_ = None\n" + text[body:]
    path.write_text(new)


def _copy_tree(tmp_path: Path) -> Path:
    copy = tmp_path / "copy"
    for source in (REPO / "orpheus").rglob("*.py"):
        target = copy / source.relative_to(REPO)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
    return copy


# ── the readings' values ────────────────────────────────────────────────────────


def test_m4_9_the_exact_medium_is_read_in_the_asking_process(tmp_path, monkeypatch):
    """M4.9, re-posed by the user's ruling of 2026-10-04 (Q5: "It can execute if it's fast"): the exact infinite
    medium reads in 0.01 s, against 0.11 s warm and 6.3 s cold through a generating process, so it is not
    memoised. Its construction and every reading start no interpreter and write no entry, its ``evaluate`` is a
    plain method, and its refusal of a point value is raised in this process. Its certificate is unaffected
    either way: ``ReferenceSolution`` re-checks every claim against the evaluation at construction."""
    from orpheus.derivations.common.exact_homogeneous import ExactInfiniteMediumDerivation, exact_infinite_medium_reference
    from orpheus.numerics.observable import Eigenvalue, PointValue

    assert not isinstance(ExactInfiniteMediumDerivation.__dict__["evaluate"], api.name("TracedMemo"))
    spawns = api.SpawnCounter(monkeypatch)
    with api.cache_root(tmp_path):
        reference = exact_infinite_medium_reference(_exact_specification())
        reference.read(Eigenvalue())
        with pytest.raises(ValueError, match="no position"):
            reference.derivation.evaluate(PointValue(position=0.5, group=0))
    assert spawns.count == 0
    assert not api.entry_files(tmp_path)


def test_m4_10_the_default_store_is_gitignored():
    """M4.10: the default cache root is ``.cache/references`` under the checkout, and git ignores it (a
    ``git add`` of the tree can never commit an entry). The control: a tracked path is not ignored."""
    root = api.name("default_root")()
    assert root == REPO / ".cache" / "references", root
    probe = subprocess.run(["git", "check-ignore", "-q", str(root / "x" / "key.npz")], cwd=REPO)
    control = subprocess.run(["git", "check-ignore", "-q", str(REPO / "orpheus" / "__init__.py")], cwd=REPO)
    assert probe.returncode == 0 and control.returncode == 1


# ── the characteristic reading: the successors of M4.1b-M4.8 ─────────────────────
#
# P1 step (e1b) of the characteristic-reference campaign: step (e2) deleted the trajectory-resolvent family, and with
# it the three old clients rows M4.1b-M4.8 were posed on. ``CharacteristicDerivation.evaluate`` is the one memo client
# that survives; since #592 its solve is a second memo, ``CharacteristicDerivation.solve``, a child entry of every
# reading of the derivation (``ChildPin``: one, the solve), so every observable of one derivation shares one solve.
# The successors, row by row (``scratch/characteristic_architecture/p1_step_e/ta_e1b/README_consumers.md``):
#
# * M4.1b (the route, the activation leg): ``tests/gates/derivations/test_characteristic_reference.py::
#   test_c3_an_equal_derivation_reads_from_the_memo_and_a_new_resolution_does_not`` (a cold read starts an
#   interpreter, an equal derivation's read starts none and reads the same bits);
# * M4.2 (one child entry, two consumers) and M4.4's ``the-solver`` witness: the child-memo contract, on a synthetic
#   client: ``test_traced_memo_process.py::test_m3_12e_a_child_generated_by_its_parent_serves_a_direct_caller`` and
#   ``::test_m3_12b_a_mutation_only_a_child_traced_makes_the_parent_stale``;
# * M4.7 (served arrays read-only): ``test_traced_memo_process.py::test_m3_14_a_served_payload_is_fresh_read_only_
#   and_bit_identical`` (a characteristic reading serves a float, no array);
# * M4.3, M4.4's other three witnesses, M4.5 and M4.8: the rows below.

_CHARACTERISTIC_DRIVER_SETUP = '''
from orpheus.data.cells import CellCoefficient, Channel
from orpheus.data.materials import Materials
from orpheus.derivations.common.xs_library import get_xs, make_mixture
from orpheus.derivations.continuous.characteristic import CharacteristicDerivation, Resolution, TransportResolution
from orpheus.geometry import BC, StructuredGeometry
from orpheus.numerics.question import Eigen
from orpheus.specification.specification import GeometrySpecification


def tiny_sphere(radius=1.0):
    """A one-region sphere of mixture A (two groups, P0 scattering) under vacuum, at the smallest resolution."""
    xs = get_xs("A", "2g")
    mixture = make_mixture(sig_t=xs["sig_t"], sig_c=xs["sig_c"], sig_f=xs["sig_f"], nu=xs["nu"], chi=xs["chi"],
                           sig_s=xs["sig_s"])
    geometry = StructuredGeometry.sphere((0.0, radius), (0,), outer=BC.vacuum)
    question = Eigen(CellCoefficient.every(Channel.FISSION_EMISSION))
    specification = GeometrySpecification(Materials({0: mixture}), geometry, question)
    return CharacteristicDerivation(specification, Resolution(1, 0, 0.4, TransportResolution(4, 4, 4), 4))
'''

_tiny_namespace: dict[str, Any] = {}
exec(_CHARACTERISTIC_DRIVER_SETUP, _tiny_namespace)


def _tiny_characteristic_sphere(radius: float = 1.0) -> Any:
    """The characteristic sphere the rows below read (the same text the copied-tree driver execs)."""
    return _tiny_namespace["tiny_sphere"](radius)


def _characteristic_entries(root: Path) -> list[Path]:
    return [f for f in api.entry_files(root) if api.entry_function(f) == api.function_id("characteristic_reading")]


def _solve_entries(root: Path) -> list[Path]:
    return [f for f in api.entry_files(root) if api.entry_function(f) == api.spelled_function_id("characteristic_solve")]


#: What the solve runs and the reading, since #592, does not: the closure's least solution and the pencil's mode.
_SOLVE_ONLY = {"LinePeriod.inflow", "DensePencil.fundamental"}
#: What every reading of a derivation runs in its own process: the derivation's construction.
_CONSTRUCTION = {"CharacteristicDerivation.__post_init__", "Walls.of", "PanelBasis.of", "RegionCrossSections.of",
                 "GalerkinSystem.__post_init__"}


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_process.py::test_m3_9_arguments_cross_through_their_constructor',
                      'tests/gates/numerics/test_traced_memo_process.py::test_m3_12a_a_parent_entry_pins_its_child_by_reference',
                      'tests/gates/derivations/test_characteristic_reference.py::test_c3_an_equal_derivation_reads_from_the_memo_and_a_new_resolution_does_not')
def test_m4_3c_the_characteristic_reading_manifest_holds_its_construction_and_pins_one_solve(tmp_path):
    """M4.3, re-posed on the characteristic reading, and #592's first done-when: the reading entry's manifest holds
    the derivation's construction (``CharacteristicDerivation.__post_init__``, ``Walls.of``, ``PanelBasis.of``,
    ``RegionCrossSections.of``, ``GalerkinSystem.__post_init__``), which a pickled derivation would hide, and NOT the
    solve (``LinePeriod.inflow``, ``DensePencil.fundamental``: they ran in the solve's child); it carries exactly ONE
    ``ChildPin``, whose function id is ``CharacteristicDerivation.solve`` and whose key names the solve entry in the
    store; the solve entry's own manifest holds the solve and pins no child. The cylinder's polar rule, which a
    sphere never runs, is in neither.

    Until #592 this row asserted the reading pinned no child (one memo level). First red (``[M]`` 2026-10-10, the
    tree at ``2f56e69b``, ``scratch/characteristic_architecture/p1_step_e/ta_592/``): ``ChildPin`` is ``[]`` and the
    solve's functions are in the reading's manifest."""
    from orpheus.numerics.observable import Eigenvalue

    with api.cache_root(tmp_path):
        _tiny_characteristic_sphere().evaluate(Eigenvalue())
    (entry_file,) = _characteristic_entries(tmp_path)
    entry = api.read_entry(entry_file)
    functions = {q for _, q, *_ in api.manifest_rows(entry, "DefPin")}
    assert _CONSTRUCTION <= functions, sorted(_CONSTRUCTION - functions)
    assert not (_SOLVE_ONLY & functions), sorted(_SOLVE_ONLY & functions)
    assert "polar_rule" not in functions and "impact_per_polar" not in functions
    children = api.manifest_rows(entry, "ChildPin")
    assert [c[0] for c in children] == [api.spelled_function_id("characteristic_solve")], children
    (solve_file,) = _solve_entries(tmp_path)
    assert solve_file.stem == children[0][1], (solve_file.stem, children[0][1])
    solve = api.read_entry(solve_file)
    solve_functions = {q for _, q, *_ in api.manifest_rows(solve, "DefPin")}
    assert _SOLVE_ONLY <= solve_functions, sorted(_SOLVE_ONLY - solve_functions)
    assert "polar_rule" not in solve_functions
    assert api.manifest_rows(solve, "ChildPin") == []


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_clients.py::test_m4_3c_the_characteristic_reading_manifest_holds_its_construction_and_pins_one_solve')
def test_m4_11_two_observables_of_one_cold_derivation_share_one_solve(tmp_path):
    """#592's second done-when, read from the store: on a cold root, an eigenvalue and then a point value of one
    derivation write two reading entries whose ``ChildPin`` rows name the SAME solve key, the store holds exactly
    one solve entry, and that entry's file is not rewritten by the second reading (its inode and modification time,
    read after each reading, are unchanged: the store writes an entry by an atomic rename, so a regenerated solve
    is a new file).

    What is counted: solve entries in the store (1), distinct child keys over the reading entries (1), rewrites of
    the solve entry across the second reading (0). Not a spy: a monkeypatch reaches only this process, and the
    solve runs in the readings' generating interpreters. First reds (``[M]`` 2026-10-10,
    ``scratch/characteristic_architecture/p1_step_e/ta_592/``): on ``2f56e69b`` no solve entry and no child key
    (every reading solved in its own process); on a copied tree whose ``answer`` calls the solve with a fresh salt
    (mutant ``solve-salted``) the two readings name two child keys. The inode leg has no witness of its own: a
    regenerated solve under one key is not spellable without a store mutation, and the key leg reds first."""
    from orpheus.numerics.observable import Eigenvalue, PointValue

    with api.cache_root(tmp_path):
        _tiny_characteristic_sphere().evaluate(Eigenvalue())
        solves = _solve_entries(tmp_path)
        assert len(solves) == 1, solves
        first = solves[0].stat()
        _tiny_characteristic_sphere().evaluate(PointValue(position=0.4, group=1))
    readings = _characteristic_entries(tmp_path)
    assert len(readings) == 2, readings
    keys = {c[1] for f in readings for c in api.manifest_rows(api.read_entry(f), "ChildPin")}
    assert keys == {solves[0].stem}, keys
    assert _solve_entries(tmp_path) == solves
    again = solves[0].stat()
    assert (again.st_ino, again.st_mtime_ns) == (first.st_ino, first.st_mtime_ns), "the second reading regenerated the solve"


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_clients.py::test_m4_11_two_observables_of_one_cold_derivation_share_one_solve',
                      'tests/gates/numerics/test_traced_memo_process.py::test_m3_14_a_served_payload_is_fresh_read_only_and_bit_identical')
def test_m4_12_a_served_solve_is_read_only_and_its_readings_are_the_in_process_ones(tmp_path):
    """#592's third done-when: after a reading has written the solve entry, a fresh derivation's ``solve()`` in this
    process is a ``Hit``; its arrays (the flux and the emission) are read-only, and a write into one raises; they
    are the in-process (``bypass``) solve's bit for bit; and readings taken off the warm solve through the memo (a
    point value of group 1, a flux integral of the fission production, each a new reading entry whose generation
    hits the solve) equal the in-process readings bit for bit (``float.hex``). First reds (``[M]`` 2026-10-10,
    ``scratch/characteristic_architecture/p1_step_e/ta_592/``): on ``2f56e69b`` there is no ``solve`` (``KeyError``);
    mutant ``served-writeable`` (the store's decode leaving arrays writeable): the read-only leg; mutant
    ``solve-salted``: the ``Hit`` leg (the salted key is ``Absent``). The bit-identity legs carry no separate
    mutant here (``test_m4_5c`` is the readings' served-equals-fresh row)."""
    from orpheus.numerics.mesh_free_function import RegionwiseConstant
    from orpheus.numerics.observable import Eigenvalue, FluxIntegral, PointValue

    production = _tiny_characteristic_sphere().specification.materials[0]
    observables = (PointValue(position=0.4, group=1), FluxIntegral(RegionwiseConstant(np.asarray(production.SigP)[None, :])))
    with api.cache_root(tmp_path):
        _tiny_characteristic_sphere().evaluate(Eigenvalue())
        derivation = _tiny_characteristic_sphere()
        assert api.verdict_kind(type(derivation).__dict__["solve"].lookup(derivation)) == "Hit"
        served = derivation.solve()
        warm = [_tiny_characteristic_sphere().evaluate(o).value.hex() for o in observables]
    with api.bypass():
        fresh_answer = _tiny_characteristic_sphere().solve()
        fresh = [_tiny_characteristic_sphere().evaluate(o).value.hex() for o in observables]
    for field in ("flux", "emission"):
        array = getattr(served, field)
        assert not array.flags.writeable, field
        assert array.tobytes() == getattr(fresh_answer, field).tobytes(), field
    with pytest.raises(ValueError, match="read-only"):
        served.flux[0, 0] = 0.0
    assert served.k.hex() == fresh_answer.k.hex()
    assert warm == fresh, (warm, fresh)


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_clients.py::test_m4_11_two_observables_of_one_cold_derivation_share_one_solve',
                      'tests/gates/derivations/test_characteristic_reference.py::test_an_eigenvalue_of_a_source_answer_is_refused_before_any_solve')
def test_m4_13_a_refusal_before_the_solve_writes_no_solve_entry(tmp_path):
    """#592's fourth done-when: an eigenvalue of a source question, read through the memo, is refused with its
    message ("an eigenvalue is read from an eigen question") and the store holds no entry at all, a solve's or a
    reading's; the activation leg: the same derivation's flux integral then writes one solve entry and one reading
    entry, so the store can record a solve here. First red (``[M]`` 2026-10-10, ``2f56e69b``): the activation leg
    (no solve entry exists to be written). The main leg's own first red is the refusal placed after the solve
    (``self.answer`` read before the match), measured on a copied tree (the generating interpreter imports the
    tree, so no in-process mutation reaches it): mutant ``refusal-after-solve``, the solve entry is written before
    the refusal and the row reds on the empty-store leg."""
    from dataclasses import replace

    from orpheus.numerics.mesh_free_function import RegionwiseConstant
    from orpheus.numerics.observable import Eigenvalue, FluxIntegral
    from orpheus.numerics.question import FixedSource

    base = _tiny_characteristic_sphere()
    source = replace(base.specification, question=FixedSource(RegionwiseConstant(np.ones((1, 2)))))
    with api.cache_root(tmp_path):
        with pytest.raises(ValueError, match="an eigenvalue is read from an eigen question"):
            type(base)(source, base.resolution).evaluate(Eigenvalue())
        assert api.entry_files(tmp_path) == [], api.entry_files(tmp_path)
        type(base)(source, base.resolution).evaluate(FluxIntegral(RegionwiseConstant(np.ones((1, 2)))))
    assert len(_solve_entries(tmp_path)) == 1 and len(_characteristic_entries(tmp_path)) == 1, api.entry_files(tmp_path)


_CHARACTERISTIC_DRIVER = textwrap.dedent('''
    import sys
    copy, cache, action = sys.argv[1:4]
    sys.path[:0] = [copy]
    import orpheus
    assert orpheus.__file__.startswith(copy), orpheus.__file__
''') + _CHARACTERISTIC_DRIVER_SETUP + textwrap.dedent('''
    from orpheus.numerics.observable import Eigenvalue
    from orpheus.numerics.traced_memo import cache_root
    with cache_root(cache):
        d = tiny_sphere()
        if action == "generate":
            print("K", d.evaluate(Eigenvalue()).value.hex())
        print("VERDICT", type(type(d).evaluate.lookup(d, Eigenvalue())).__name__)
''')

#: (row, file under orpheus/, the def or statement to find, the edit, the reading's verdict after it)
CHARACTERISTIC_TREE_WITNESSES = [
    ("construction-helper", "derivations/continuous/characteristic/walls.py", "def _wall_of(", None, "Stale"),
    ("rule-of-another-chart", "derivations/continuous/characteristic/lines.py", "def polar_rule(", None, "Hit"),
    ("reading-module-constant", "derivations/continuous/characteristic/reference.py",
     'Lift = Callable[["sympy.Expr"], "sympy.Expr"]', "constant-before", "Stale"),
]


def _drive_characteristic(copy: Path, cache: Path, action: str) -> str:
    run = subprocess.run(
        [sys.executable, "-O", "-c", _CHARACTERISTIC_DRIVER, str(copy), str(cache), action],
        capture_output=True, text=True, cwd=str(copy.parent), env={k: v for k, v in os.environ.items() if k != "PYTHONPATH"},
    )
    lines = [line for line in run.stdout.splitlines() if line.startswith("VERDICT ")]
    if not lines:
        raise AssertionError(f"the driver failed:\n{run.stderr[-3000:]}")
    return lines[-1][len("VERDICT "):]


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_process.py::test_m3_2_an_edit_to_what_ran_regenerates_and_nothing_else_does',
                      'tests/gates/numerics/test_traced_memo_clients.py::test_m4_3c_the_characteristic_reading_manifest_holds_its_construction_and_pins_one_solve')
@pytest.mark.parametrize("row,file,anchor,method,reading", CHARACTERISTIC_TREE_WITNESSES,
                         ids=[w[0] for w in CHARACTERISTIC_TREE_WITNESSES])
def test_m4_4c_an_edit_to_the_copied_tree_misses_exactly_where_the_characteristic_reading_ran(tmp_path, row, file, anchor,
                                                                                            method, reading):
    """M4.4, re-posed on the characteristic reading, on a copy of ``orpheus/``: a neutral edit to a construction
    helper the reading ran (``_wall_of``) makes the reading stale; an edit to the cylinder's polar rule, which a
    sphere reading never runs, leaves it a hit; a new constant in the reading module's skeleton makes it stale.
    The copy's own driver reports ``Hit`` before the edit (the activation leg)."""
    copy = _copy_tree(tmp_path)
    cache = tmp_path / "cache"
    assert _drive_characteristic(copy, cache, "generate") == "Hit"
    path = copy / "orpheus" / file
    if method == "constant-before":
        text = path.read_text()
        if text.count(anchor) != 1:
            raise AssertionError(f"{anchor!r} occurs {text.count(anchor)} times in {path}")
        path.write_text(text.replace(anchor, "_NEUTRAL_CONSTANT = None\n" + anchor))
    else:
        _neutral_edit(path, anchor, method)
    assert _drive_characteristic(copy, cache, "lookup") == reading, row


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_validation.py::test_m2_6_floats_cross_the_payload_bit_for_bit',
                      'tests/gates/derivations/test_characteristic_reference.py::test_c3_an_equal_derivation_reads_from_the_memo_and_a_new_resolution_does_not')
def test_m4_5c_a_served_characteristic_reading_is_the_fresh_reading_bit_for_bit(tmp_path):
    """M4.5, re-posed: the memo's characteristic readings equal the in-process readings under ``bypass()`` bit for
    bit: k, a point value of group 1 and a flux integral of the fission production (``float.hex``). The second
    pass through the memo reads each from an entry (three entries written)."""
    from orpheus.numerics.mesh_free_function import RegionwiseConstant
    from orpheus.numerics.observable import Eigenvalue, FluxIntegral, PointValue

    production = _tiny_characteristic_sphere().specification.materials[0]
    weight = FluxIntegral(RegionwiseConstant(np.asarray(production.SigP)[None, :]))
    observables = (Eigenvalue(), PointValue(position=0.4, group=1), weight)
    with api.cache_root(tmp_path):
        generated = [_tiny_characteristic_sphere().evaluate(o).value.hex() for o in observables]
        served = [_tiny_characteristic_sphere().evaluate(o).value.hex() for o in observables]
    assert len(_characteristic_entries(tmp_path)) == 3
    with api.bypass():
        fresh = [_tiny_characteristic_sphere().evaluate(o).value.hex() for o in observables]
    assert served == generated == fresh


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_process.py::test_m3_11_a_raising_generator_writes_nothing_and_raises_its_own_error',
                      'tests/gates/derivations/test_characteristic_reference.py::test_c4_a_refusal_crosses_the_memo_with_its_type')
def test_m4_8c_a_characteristic_refusal_is_never_cached(tmp_path, monkeypatch):
    """M4.8, re-posed: the characteristic reference is a direct solve (no iteration to leave unconverged), so its
    refusal after the solve is the one M4.8's contract applies to: the flux of a mode nearest tau. Read through the
    memo it crosses back as ``NotImplementedError`` with its message and writes nothing, so a second reading starts
    an interpreter again; the mode's eigenvalue, read beside it, is cached. First red: a refusal cached as a value
    (the second read starts no interpreter)."""
    from dataclasses import replace

    from orpheus.numerics.observable import Eigenvalue, PointValue
    from orpheus.numerics.question import Nearest

    base = _tiny_characteristic_sphere()
    question = replace(base.specification.question, mode=Nearest(tau=0.05))
    derivation = type(base)(replace(base.specification, question=question), base.resolution)
    with api.cache_root(tmp_path):
        with pytest.raises(NotImplementedError, match="a higher mode has no flux scale"):
            derivation.evaluate(PointValue(position=0.4, group=0))
        spawns = api.SpawnCounter(monkeypatch)
        with pytest.raises(NotImplementedError, match="a higher mode has no flux scale"):
            type(base)(derivation.specification, base.resolution).evaluate(PointValue(position=0.4, group=0))
        refused_again = spawns.count
        derivation.evaluate(Eigenvalue())
    assert refused_again >= 1, f"the second refusal started {refused_again} interpreters: it was served from an entry"
    assert len(_characteristic_entries(tmp_path)) == 1


# ── #592 review: the (n,2n) gauge, whose key a hash seed could split ─────────────────────────────

_GAUGE_DRIVER = textwrap.dedent('''
    import os, sys
    repo, cache = sys.argv[1:3]
    sys.path.insert(0, repo)
    from orpheus.numerics import content
    if os.environ.get("ARM") == "iteration-order":
        content._Exactly.elements = lambda self, encoded: encoded     # the pre-fix spelling, in the parent only
    from tests.gates.derivations.test_characteristic_reference import _K, _TINY, _hetero_sphere
    from orpheus.derivations.continuous.characteristic import CharacteristicDerivation
    from orpheus.numerics.observable import Eigenvalue, PointValue
    from orpheus.numerics.question import Eigen
    from orpheus.numerics.traced_memo import cache_root
    make = lambda: CharacteristicDerivation(_hetero_sphere(Eigen(_K)), _TINY)
    print("ORDER", repr([c for c in make().specification.question.gauge.cells]))
    if sys.argv[3] == "run":
        with cache_root(cache):
            make().evaluate(Eigenvalue())
            make().evaluate(PointValue(position=0.4, group=1))
            d = make()
            print("LOOKUP", type(type(d).__dict__["solve"].lookup(d)).__name__)
            print("PARENT_KEY", type(d).__dict__["solve"].key(d))
            d.answer
''')


def _gauge_driver(tmp_path: Path, seed: str, action: str, arm: str = "none") -> dict[str, str]:
    script = tmp_path / "gauge_driver.py"
    script.write_text(_GAUGE_DRIVER)
    env = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"} | {"PYTHONHASHSEED": seed, "ARM": arm}
    run = subprocess.run([sys.executable, "-O", str(script), str(REPO), str(tmp_path / f"cache-{seed}-{arm}"), action],
                         capture_output=True, text=True, env=env, cwd=str(tmp_path), timeout=600)
    rows = dict(line.split(" ", 1) for line in run.stdout.splitlines() if " " in line)
    if "ORDER" not in rows or (action == "run" and "PARENT_KEY" not in rows):
        raise AssertionError(f"the driver failed under PYTHONHASHSEED={seed}:\n{run.stderr[-3000:]}")
    return rows


def _parent_seed_that_reorders_the_gauge(tmp_path: Path) -> str:
    """A hash seed under which the (n,2n) gauge's two cells iterate in another order than under seed 0, where the
    generations run: the precondition for the split to be observable (the instrument's positive control)."""
    generation = _gauge_driver(tmp_path, "0", "order")["ORDER"]
    for seed in ("1", "2", "3", "4", "5", "6", "7", "8"):
        if _gauge_driver(tmp_path, seed, "order")["ORDER"] != generation:
            return seed
    raise AssertionError("no seed in 1..8 reorders the gauge's cells: the leg would be blind")


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_clients.py::test_m4_12_a_served_solve_is_read_only_and_its_readings_are_the_in_process_ones',
                      'tests/gates/numerics/test_traced_memo_findings.py::test_s1_one_call_is_one_key_under_every_hash_seed')
def test_m4_11b_m4_12b_a_parent_at_another_hash_seed_shares_the_readings_solve(tmp_path):
    """The (n,2n)-gauge leg of M4.11 and M4.12 (qa, #592 review): ``_tiny_characteristic_sphere``'s gauge holds ONE
    cell, so no hash seed can reorder it (a space-side stabiliser); ``_hetero_sphere(Eigen(_K))`` (an (n,2n)
    material) holds two. The parent runs at a hash seed under which those two cells iterate differently than under
    seed 0, where every generation runs (found and asserted first, the positive control). Two readings write one
    solve entry; the parent's ``solve.lookup`` of the same derivation is a ``Hit``; its key is the one the readings
    pin; and the parent's ``.answer`` adds no second solve entry. First red (``[M]`` 2026-10-10): the parent's
    ``_Exactly.elements`` rebound to the identity (env ``ARM=iteration-order``): ``Absent``, two keys, two
    entries; qa measured the same split on the pre-fix tree (``probe5.py``)."""
    seed = _parent_seed_that_reorders_the_gauge(tmp_path)
    rows = _gauge_driver(tmp_path, seed, "run")
    root = tmp_path / f"cache-{seed}-none"
    solves = _solve_entries(root)
    readings = _characteristic_entries(root)
    assert rows["LOOKUP"] == "Hit", (seed, rows["LOOKUP"])
    assert len(solves) == 1, solves
    assert solves[0].stem == rows["PARENT_KEY"]
    keys = {c[1] for f in readings for c in api.manifest_rows(api.read_entry(f), "ChildPin")}
    assert keys == {rows["PARENT_KEY"]}, keys


def test_m4_11b_the_iteration_order_arm_splits_the_parent_from_the_readings(tmp_path):
    """The witness of the row above, kept as a row (``vv`` #17): with the parent's ``_Exactly.elements`` rebound to
    the identity, at the reordering seed the parent's lookup is ``Absent`` and its ``.answer`` writes a second solve
    entry. Green here means the row above can red."""
    seed = _parent_seed_that_reorders_the_gauge(tmp_path)
    rows = _gauge_driver(tmp_path, seed, "run", arm="iteration-order")
    assert rows["LOOKUP"] == "Absent", rows["LOOKUP"]
    assert len(_solve_entries(tmp_path / f"cache-{seed}-iteration-order")) == 2
