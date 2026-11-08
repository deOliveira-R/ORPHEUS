"""B5: the reads of a test's own run, resolved module-qualified (rulings 47 and 49), and the pins they mint.

Spec ``scratch/verification_tree/cone_store/b5/spec_b5.md`` §3.1-§3.3, rows R1-R17 (the resolution), P1-P8 (the pin)
and F1-F3 (the floor, and the selection cone left alone). The fixture is ``_READS`` in ``_cone_store_fixture.py``
(``build(root, reads=True)``), a layer ``pkg/rd`` with a test per shape of read; every expectation below is read off
those sources by hand.

The interface (the orchestrator, 2026-11-07): ``tools.cone_store.reads.Read(path, name, channel)``, ``path`` relative
to the repository, ``name`` ``None`` for a read of a module object (C8), ``channel`` one of ``direct``, ``decorator``,
``string``, ``via_apparatus``, ``reexport``, ``closure``, ``module``, or ``reads.Channel.UNRESOLVED`` (C11, ``name`` the dotted
spelling); ``ConeStore.own_reads(test_id)`` returns them sorted; the B5 statement pins are ordinary ``StatementsPin``
rows of the own cone (C9). Most fixture tests raise ``AssertionError`` by hand: the assertion rewriter adds bare loads
of each operand's module, which only R11 is about, so the rows below ignore channel ``module`` unless it is their
subject.
"""
from __future__ import annotations

import json
import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest

from orpheus.numerics import traced_memo as memo

from ._cone_store_fixture import build, capture, copy_captured, edit, environment, paths

pytestmark = pytest.mark.foundation

T = "tests/rd/test_rd.py::"
LOW, MID, REX, HIGH = "pkg/rd/low.py", "pkg/rd/mid.py", "pkg/rd/rex.py", "pkg/rd/high.py"
INIT, OTHER, HELPER = "pkg/rd/__init__.py", "pkg/other/o.py", "tests/rd/_helper.py"


@pytest.fixture(scope="module")
def captured(tmp_path_factory) -> Path:
    root = build(tmp_path_factory.mktemp("reads") / "tree", memo=False, reads=True)
    capture(root, "tests/rd")
    return root


def _store(root: Path):
    from tools.cone_store.store import ConeStore

    return ConeStore(paths(root)["store"])


@pytest.fixture(scope="module")
def reads(captured: Path) -> dict[str, set[tuple[str, str | None, str]]]:
    """Every stored test's ``own_reads`` as ``(path, name, channel)`` triples."""
    with _store(captured) as store:
        found = {test_id: {(r.path, r.name, str(r.channel)) for r in store.own_reads(test_id)} for test_id in store.tests()}
    assert len(found) >= 30, f"{len(found)} tests stored: the capture did not run the fixture"
    return found


def rows(reads: dict, name: str, *, module: bool = False) -> set[tuple[str, str | None, str]]:
    """The test's reads in the package ``pkg`` (closure rows included), channel ``module`` only when asked."""
    test_id = name if "::" in name else T + name
    assert test_id in reads, f"{test_id} not stored"
    return {r for r in reads[test_id] if r[0].startswith("pkg/") and (module or r[2] != "module")}


def node_rows(reads: dict, name: str) -> set[tuple[str, str | None, str]]:
    """As :func:`rows`, without the ``closure`` and ``reexport`` channels (those feed pins, never nodes)."""
    return {r for r in rows(reads, name) if r[2] not in ("closure", "reexport")}


# ── group R, the resolution ──


def test_r1_a_global_alias_reads_one_binding(reads: dict) -> None:
    """R1. ``low.TOL`` through a module-level ``from pkg.rd import low``: one row, in ``low``. First red: the
    as-written pairing of every ``co_names`` entry with every module the code reads."""
    assert rows(reads, "test_global") == {(LOW, "TOL", "direct")}


def test_r2_an_attribute_chain_resolves_to_its_last_module(reads: dict) -> None:
    """R2. ``pkg.rd.low.TOL`` after ``import pkg.rd.low``: the chain passes ``pkg`` and ``pkg.rd`` (submodules bound
    by the import system, no statement) and ends at ``low.TOL``. First red: the chain stopped after its first
    attribute (``pkg/rd/__init__.py`` or nothing)."""
    assert rows(reads, "test_chain") == {(LOW, "TOL", "direct")}


def test_r3_an_import_inside_the_body_binds_a_local(reads: dict) -> None:
    """R3 (C2). ``from pkg.rd import low as m`` in the body, then ``m.TOL``. First red: no code-local import table
    (the local ``m`` is bound by no module statement, and the read names nothing)."""
    assert rows(reads, "test_in_def") == {(LOW, "TOL", "direct")}


def test_r4_a_cell_read_resolves_through_the_outer_import(reads: dict) -> None:
    """R4 (C2). The inner function reads ``m`` as a cell, imported by the test's body. First red: ``LOAD_DEREF``
    resolved through the module's globals only."""
    assert rows(reads, "test_cell") == {(LOW, "TOL", "direct")}


def test_r5_a_reexport_is_followed_and_its_module_pinned(reads: dict) -> None:
    """R5 (C7). ``TOL`` imported from ``rex``, which imports it from ``low``: the read is ``low``'s, and ``rex``'s
    import statement is a ``reexport`` row (pinned, not a node). First red: stop at the first binding (``rex``)."""
    assert rows(reads, "test_reexport") == {(LOW, "TOL", "direct"), (REX, "TOL", "reexport")}


def test_r6_a_reexport_through_the_package_init(reads: dict) -> None:
    """R6 (C7). ``Rec`` from ``pkg.rd``, whose ``__init__`` imports it from ``low``. First red: stop at
    ``__init__``, or follow by name without the identity check."""
    assert node_rows(reads, "test_package_reexport") == {(LOW, "Rec", "direct")}
    assert (INIT, "Rec", "reexport") in rows(reads, "test_package_reexport")


def test_r7_a_fallback_import_stops_where_the_object_is_bound(reads: dict) -> None:
    """R7 (C7). ``low.SPEED`` is bound by ``from pkg.rd._fast import SPEED`` in a ``try`` and by ``SPEED = 0`` in
    its ``except``; ``pkg.rd._fast`` does not exist, so the read stops at ``low``. ``[M]`` 2026-11-07: dropping the
    identity check leaves this row green (the absent module stops the chain first); the identity check's catcher is
    R7b. First red: a resolver that follows an import to a module that does not exist (the read is dropped)."""
    assert node_rows(reads, "test_fallback") == {(LOW, "SPEED", "direct")}
    assert not any("_fast" in path for path, _, _ in rows(reads, "test_fallback"))


def test_r7b_an_import_rebound_by_a_later_statement_stops_where_the_object_is_bound(reads: dict) -> None:
    """R7b (C7; added 2026-11-07 because the battery's arm C7b, the identity check dropped, left R7 green: ``_fast``
    does not exist, so R7 stops at ``low`` without the check). ``low`` imports ``REBOUND`` from ``pkg.other.alt``
    (``99``) and then binds ``REBOUND = 5``: the import's module exists and binds another object, so only the identity
    check stops at ``low``. First red ``[M]``: the identity check replaced by ``True`` (the read resolves to
    ``pkg/other/alt.py``)."""
    assert node_rows(reads, "test_rebound") == {(LOW, "REBOUND", "direct")}
    assert not any(path == "pkg/other/alt.py" for path, _, _ in reads[T + "test_rebound"])


def test_r8_a_name_bound_by_a_def_and_an_assignment(captured: Path, reads: dict) -> None:
    """R8. ``low.f`` is bound by ``def f`` and by ``f = _wrap(f)``. The test runs it, so ``low`` holds a def pin, and
    ``own_reads`` holds ``low:f``; the projection over both statements is P2's and P4's. First red: the read of a name
    whose code ran dropped as already pinned."""
    assert (LOW, "f", "direct") in rows(reads, "test_def_and_assign")
    with _store(captured) as store:
        defs = {pin.path for pin in store.own_pins(T + "test_def_and_assign") if isinstance(pin, memo.DefPin)}
    assert LOW in defs, defs


def test_r9_a_dataclass_is_read_by_its_class_statement(reads: dict) -> None:
    """R9 (attack a12). ``low.Rec(1.0) == low.Rec(1.0)`` runs the generated ``__init__`` and ``__eq__``, code whose
    file is ``<string>``: no row comes from it. First red: generated code resolved as first-party."""
    assert node_rows(reads, "test_dataclass") == {(LOW, "Rec", "direct")}
    assert not any(path.startswith("<") for path, _, _ in reads[T + "test_dataclass"])


def test_r10_an_isinstance_against_a_protocol_names_the_protocol(reads: dict) -> None:
    """R10 (attack a13). The protocol check runs in ``typing``; the test's own load of ``low.Proto`` is the read.
    A control that a stdlib-dispatched read still names its class."""
    assert node_rows(reads, "test_protocol") == {(LOW, "Proto", "direct")}


def test_r11_a_constant_string_pairs_with_the_module_read(reads: dict) -> None:
    """R11 (C3). ``assert getattr(low, "TOL") > 0`` and ``assert hasattr(low, "ALIAS")``, under pytest's assertion
    rewriter, which moves each constant into a temporary so no load chain pairs it with ``low``. First red: C3
    dropped (the two names are not read)."""
    assert {r for r in rows(reads, "test_getattr") if r[2] == "string"} == {(LOW, "TOL", "string"),
                                                                           (LOW, "ALIAS", "string")}


def test_r12_a_name_tried_and_absent_reads_the_module(reads: dict) -> None:
    """R12 (C8). ``hasattr(low, "GONE")``: ``GONE`` is bound nowhere, so the read is of the module object, name
    ``None``; P6 shows its pin keeps ``GONE``. First red: a chain ending on a module dropped."""
    found = rows(reads, "test_absent", module=True)
    assert (LOW, None, "module") in found
    assert not any(name == "GONE" for _, name, _ in found)


def test_r13_a_module_dunder_reads_the_module(reads: dict) -> None:
    """R13 (C8). ``low.__file__`` is bound by the import system, not by a statement: a read of the module. First red:
    the dunder counted as an unresolved name."""
    assert (LOW, None, "module") in rows(reads, "test_dunder", module=True)
    assert not any(channel == "unresolved" for _, _, channel in reads[T + "test_dunder"])


def test_r14_a_star_import_is_followed_by_identity(reads: dict) -> None:
    """R14 (C11). ``from pkg.rd.low import *`` binds ``TOL`` with no statement naming it; followed by identity, the
    read is ``low:TOL`` and nothing is unresolved. First red: names no statement binds dropped or counted."""
    found = reads["tests/rd/test_star.py::test_star"]
    assert (LOW, "TOL", "direct") in found
    assert not any(channel == "unresolved" for _, _, channel in found)


def test_r15_a_class_attribute_read_through_an_instance_is_the_declared_boundary(reads: dict) -> None:
    """R15 (C11, one-sided). ``box.LIMIT``, ``box`` from a conftest fixture: the receiver is not a module, so
    ``low:LIMIT`` is NOT read (the declared boundary; this row asserts the limitation, it does not certify a
    capability). The fixture's own code reads ``low.Box``, so ``low`` is still reached. First red: a resolver that
    pairs an attribute with every class of the layer (it names ``low:LIMIT``)."""
    found = rows(reads, "test_instance_attr")
    assert (LOW, "Box", "direct") in found
    assert not any(name == "LIMIT" for _, name, _ in found)


def test_r16_a_test_helper_that_only_reexports_is_not_apparatus(reads: dict) -> None:
    """R16 (C7). ``BUILT`` imported from ``tests.rd._helper``, which imports it from ``mid``: the read is ``mid``'s;
    the helper appears at most as a ``reexport`` row; ``BUILT``'s builder (``low``) appears only as ``closure``.
    First red: a re-exporting test file read as apparatus."""
    test_id = T + "test_helper_reexport"
    assert node_rows(reads, test_id) == {(MID, "BUILT", "direct")}
    assert {channel for path, _, channel in reads[test_id] if path == HELPER} <= {"reexport"}
    assert {channel for path, _, channel in reads[test_id] if path == LOW} <= {"closure"}


def test_r17_a_parametrized_methods_decorator_is_read(reads: dict) -> None:
    """R17 (C12). ``TestParam.test_parametrized`` is parametrized over ``low.TOL``: the value reached the test as an
    argument, evaluated at collection. First red: the decorators of the top-level statement only (the class)."""
    ids = [t for t in reads if t.startswith(T + "TestParam::test_parametrized[")]
    assert len(ids) == 1, ids
    assert (LOW, "TOL", "decorator") in rows(reads, ids[0])


def test_r18_an_import_stored_and_loaded_by_one_superinstruction(reads: dict) -> None:
    """R18 (the elegance review's F8, p1 case ``h``; the orchestrator, 2026-11-07). ``import pkg as pk; value =
    pk.rd.high.h()`` on ONE line compiles to ``IMPORT_NAME pkg`` then ``STORE_FAST_LOAD_FAST (pk, pk)`` (a store and a
    load fused; ``[M]`` Python 3.14.3, ``dis``: on two lines, or with a ``LOAD_FAST_BORROW``, the two stay apart). The
    reader must take the store as the import's binding AND the load as the start of the chain ``pk.rd.high.h``.
    ``[M]`` 0 sites in the tree (the review's p2), so this row is the only witness. First red: the fused opcode
    unlisted (it falls to the walker's ``else``, and both the binding and the chain are dropped): RED on the build
    of 2026-11-07, by design, until the reader learns the opcode."""
    assert (HIGH, "h", "direct") in rows(reads, "test_superinstruction")
    # the fused store ends the import: ``value``, stored next, is not bound to ``pkg`` (a stray binding would read the
    # load of ``value`` as a read of the module ``pkg``)
    assert not any(path == "pkg/__init__.py" for path, _, _ in reads[T + "test_superinstruction"])


def test_r19_a_name_met_unbound_under_closure_is_not_the_tests_read(reads: dict) -> None:
    """R19 (the revision's ``record``, 2026-11-07). ``mid.FROM_LATE = low.LATE2``, and ``LATE2`` is set by
    ``globals()["LATE2"] = 2``, a store no statement's bindings name. ``test_closure_unbound`` reads ``mid.FROM_LATE``:
    following its provenance meets ``LATE2`` unbound, which joins ``low``'s pin (P9) and is NOT the test's unresolved
    read. First red: the unbound name recorded whatever the channel (an ``unresolved`` row ``pkg.rd.low.LATE2``)."""
    found = reads[T + "test_closure_unbound"]
    assert (MID, "FROM_LATE", "direct") in found
    assert not any(channel == "unresolved" for _, _, channel in found), sorted(found)


def test_r20_a_cell_is_read_through_its_own_enclosing_def(reads: dict) -> None:
    """R20 (QA's F3 on the build, 2026-11-07). Two helpers nested in the test each bind ``m`` by an import (``low``,
    then ``high``) and each hold an ``inner`` reading ``m``: each ``inner`` resolves through the cells of its own
    enclosing def. First red: cells keyed by file alone (both ``inner``s read one ``m``, and one of ``low:TOL`` and
    ``high:HIGH`` is lost)."""
    found = node_rows(reads, "test_two_cells")
    assert {(LOW, "TOL", "direct"), (HIGH, "HIGH", "direct")} <= found, sorted(found)
    assert not any(name in ("TOL", "HIGH") and path not in (LOW, HIGH) for path, name, _ in found)


def test_r21_an_attribute_after_a_fused_pair_of_loads(reads: dict) -> None:
    """R21 (the revision's ``_FUSED``). ``max(x, m.TOL)`` compiles to ``LOAD_FAST_BORROW_LOAD_FAST_BORROW (x, m)``
    then ``LOAD_ATTR TOL`` (``[M]`` Python 3.14.3, ``dis``): the attribute extends the SECOND half's chain. First red:
    the fused pair read as two bare loads (``m.TOL`` lost) or the attribute given to ``x``."""
    assert (LOW, "TOL", "direct") in rows(reads, "test_fast_pair")


# ── group P, the pin ──

_VALIDATE = r'''
import json, sys
from orpheus.numerics import traced_memo as memo
from tools.cone_store.store import ConeStore

checkout = memo._Checkout(memo.Store.current())
found = {}
with ConeStore(sys.argv[1]) as store:
    for test_id in json.loads(sys.argv[2]):
        found[test_id] = [[pin.path, pin.stale(checkout)] for pin in store.own_pins(test_id)
                          if type(pin).__name__ == "StatementsPin"]
print("B5-VALIDATE " + json.dumps(found))
'''


def _validate(root: Path, *names: str) -> dict[str, list[tuple[str, str | None]]]:
    """Each named test's own-cone ``StatementsPin``s judged on ``root``'s files now (a fresh process: the memo's
    reads are process-global), as ``(path, stale reason or None)``."""
    driver = root.parent / f"{root.name}.validate.py"
    driver.write_text(_VALIDATE)
    ids = [n if "::" in n else T + n for n in names]
    run = subprocess.run([sys.executable, "-O", str(driver), str(paths(root)["store"]), json.dumps(ids)], cwd=root,
                         env=environment(root), capture_output=True, text=True, timeout=600)
    line = next((l for l in run.stdout.splitlines() if l.startswith("B5-VALIDATE ")), None)
    if run.returncode != 0 or line is None:
        raise AssertionError(f"validate driver exited {run.returncode}:\n{run.stdout[-3000:]}\n{run.stderr[-3000:]}")
    return {k.split("::")[-1]: [tuple(row) for row in v] for k, v in json.loads(line[len("B5-VALIDATE "):]).items()}


def _stale(found: list[tuple[str, str | None]], path: str) -> bool:
    on = [reason for p, reason in found if p == path]
    assert on, f"no statement pin on {path}: {found}"
    return any(reason is not None for reason in on)


#: Tests of ``test_rd.py`` that run no ``pkg`` code: every statement pin on a ``pkg`` file in their own cone is B5's.
_READ_ONLY = ("test_global", "test_chain", "test_in_def", "test_cell", "test_reexport", "test_fallback",
              "test_protocol", "test_absent", "test_dunder", "test_helper_reexport", "test_case_table",
              "test_helper_table", "test_legendre_like", "test_two_reads")


def test_p1_every_b5_pin_validates_and_covers_every_file_read(captured: Path, reads: dict) -> None:
    """P1. On the unedited tree every statement pin of the read-only tests validates, and their ``pkg`` files are
    exactly the ``pkg`` files of their ``own_reads`` (every channel, closure and reexport included: those feed pins).
    A control for P2-P7; P3 is the row that tells a resolved projection from the run's whole closure."""
    found = _validate(captured, *_READ_ONLY)
    assert len(found) == len(_READ_ONLY)
    for name, pins in found.items():
        assert all(reason is None for _, reason in pins), (name, pins)
        pinned = {path for path, _ in pins if path.startswith("pkg/")}
        read = {path for path, _, channel in reads[T + name] if path.startswith("pkg/") and channel != "unresolved"}
        assert pinned == read, (name, sorted(pinned ^ read))


@pytest.fixture
def edited(captured: Path, tmp_path: Path):
    def make(relative: str, old: str, new: str) -> Path:
        root = copy_captured(captured, tmp_path / "tree")
        edit(root, relative, old, new)
        return root

    return make


def test_p2_an_edit_to_a_read_binding_stales_the_pin(edited) -> None:
    """P2. ``TOL = 1e-12`` to ``1e-11``: ``test_global``'s ``low`` pin is stale. First red: a pin projected on no
    names (effects only)."""
    assert _stale(_validate(edited(LOW, "TOL = 1e-12", "TOL = 1e-11"), "test_global")["test_global"], LOW)


def test_p3_an_edit_to_a_binding_not_read_leaves_the_pin_valid(edited) -> None:
    """P3. ``HIGH = 7`` to ``8``: ``test_global`` read nothing of ``high`` and its ``low`` pin is valid; it holds no
    ``high`` pin. Then ``ALIAS = TOL`` to ``ALIAS = 2.0`` (a ``low`` binding ``test_global`` did not read): still
    valid. First red: the pin projected on E′, the run's whole closure, which names ``ALIAS``."""
    found = _validate(edited(HIGH, "HIGH = 7", "HIGH = 8"), "test_global")["test_global"]
    assert not _stale(found, LOW) and not any(path == HIGH for path, _ in found)


def test_p3b_an_unread_binding_of_the_same_module(edited) -> None:
    """P3, second leg (see P3)."""
    assert not _stale(_validate(edited(LOW, "ALIAS = TOL", "ALIAS = 2.0"), "test_global")["test_global"], LOW)


def test_p4_a_new_statement_binding_a_read_name_stales_the_pin(edited) -> None:
    """P4. ``TOL = 2`` appended to ``low.py``: a new statement binds a resolved name. A control of the existing
    projection rule's reach."""
    root = edited(LOW, "    SPEED = 0\n", "    SPEED = 0\n\nTOL = 2\n")
    assert _stale(_validate(root, "test_global")["test_global"], LOW)


def test_p5_deleting_a_reexport_stales_the_reexporters_pin(edited) -> None:
    """P5 (C7). The re-export removed from ``rex.py``: ``test_reexport``'s ``rex`` pin is stale. First red: only the
    terminus pinned."""
    root = edited(REX, "from pkg.rd.low import TOL", "TOL_GONE = None")
    assert _stale(_validate(root, "test_reexport")["test_reexport"], REX)


def test_p6_binding_a_name_tried_and_absent_stales_the_module_pin(edited) -> None:
    """P6 (C8). ``GONE = 1`` appended to ``low.py``: ``test_absent`` read ``low`` and tried ``GONE``, so its pin is
    stale. First red: a module read stored with empty names."""
    root = edited(LOW, "    SPEED = 0\n", "    SPEED = 0\n\nGONE = 1\n")
    assert _stale(_validate(root, "test_absent")["test_absent"], LOW)


def test_p7_a_closure_read_is_pinned(edited) -> None:
    """P7 (C6). ``Rec``'s field annotated ``int``: ``test_case_table`` reached ``low:Rec`` through the apparatus
    ``CASES`` (a node read), and ``test_legendre_like`` reached it only as ``BUILT``'s builder (``closure``, a pin,
    not a node); both ``low`` pins are stale. First red: no pin for ``closure`` reads."""
    found = _validate(edited(LOW, "    value: float\n", "    value: int\n"), "test_case_table", "test_legendre_like")
    assert _stale(found["test_case_table"], LOW)
    assert _stale(found["test_legendre_like"], LOW)


def test_p9_a_name_met_unbound_in_provenance_joins_the_pin(edited) -> None:
    """P9 (R19's pin half). ``LATE2 = 3`` appended to ``low.py``: a statement now binds the name
    ``test_closure_unbound``'s provenance met unbound, so its ``low`` pin is stale. The test runs no ``low`` or ``mid``
    code (its def is ``high.h``), so its only ``low`` pin is B5's: ``[M]`` 2026-11-07, while the test ran ``mid.g``
    (hence ``low.f``) the memo's own ``low`` pin, projected on the run's closure (which names ``LATE2`` through
    ``FROM_LATE``'s statement), staled too and hid the B5 pin's names from this row. First red: the names tried or met
    unbound dropped from the pin."""
    root = edited(LOW, "    SPEED = 0\n", "    SPEED = 0\n\nLATE2 = 3\n")
    assert _stale(_validate(root, "test_closure_unbound")["test_closure_unbound"], LOW)


def test_p8_the_reads_round_trip_and_the_form_is_four(captured: Path, tmp_path: Path) -> None:
    """P8. ``FORMAT`` is 4; the ``own_reads`` rows read the same through two openings; a store marked with form 3 is
    emptied when opened (its reads too). First red ``[M]``: ``FORMAT`` left at 3. ``[M]`` 2026-11-07: leaving
    ``own_reads`` out of the tables a store of another form drops keeps this row green, because dropping ``tests``
    deletes its rows through the ``ON DELETE CASCADE`` foreign key; the third leg guards that the two mechanisms
    together empty the table, and the drop list matters only for a change of columns."""
    from tools.cone_store import store as store_module

    assert store_module.FORMAT == 4
    with _store(captured) as one:
        first = {t: one.own_reads(t) for t in one.tests()}
    with _store(captured) as two:
        assert {t: two.own_reads(t) for t in two.tests()} == first
    assert any(first.values())
    old = tmp_path / "form3.sqlite"
    old.write_bytes(paths(captured)["store"].read_bytes())
    with sqlite3.connect(old) as db:
        db.execute("UPDATE meta SET value = '3' WHERE key = 'format'")
    with store_module.ConeStore(old) as reopened:
        assert reopened.tests() == {} or all(not reopened.own_reads(t) for t in reopened.tests())
    with sqlite3.connect(old) as db:
        tables = {name for (name,) in db.execute("SELECT name FROM sqlite_master WHERE type = 'table'")}
        if "own_reads" in tables:
            assert db.execute("SELECT COUNT(*) FROM own_reads").fetchone()[0] == 0


# ── group F, the floor and the selection cone ──


def test_f1_a_pytest_hook_is_the_floor(reads: dict) -> None:
    """F1 (C4). ``tests/rd/conftest.py::pytest_runtest_setup`` runs in every test and reads ``pkg.other.o.OTHER``:
    no test but ``test_other`` holds a read in ``pkg/other``. First red: hooks resolved (every test gains an
    above-layer read)."""
    holding = sorted(t for t, found in reads.items() if any(path == OTHER for path, _, _ in found))
    print(f"F1: {len(reads)} tests read, {len(holding)} hold a pkg/other read")
    assert holding == [T + "test_other"], holding


def test_f2_the_floors_positive_control(reads: dict) -> None:
    """F2. ``test_other`` reads ``other_module.OTHER`` itself and holds the row. First red: a resolver that skips the
    code of test files (F1 would then pass for the wrong reason)."""
    assert (OTHER, "OTHER", "direct") in reads[T + "test_other"]


def test_f3_resolving_adds_nothing_to_the_selection_cone(captured: Path, reads: dict) -> None:
    """F3 (the orchestrator, 2026-11-07). The resolver runs in the child after the own recording stops, while the
    inherited selection recording may still be active. Of the tests that open no ``pkg`` file, none holds a
    ``DataPin`` or a ``ListingPin`` on a ``pkg`` path in either cone, although each resolved a ``pkg`` module. The
    population is printed and asserted non-empty, and the audit tests are the positive control (they do hold such
    ``DataPin``s). First red ``[M]``: the resolver run under the memo's ``recorded()`` before the inherited recording
    stops (its source reads become ``DataPin``s on modules the test never opened, an over-selection outside B4's
    bound). ``[M]`` 2026-11-07: removing ``cone.py``'s ``unrecorded()`` alone, or with the reorder, leaves the row
    green: the worker's ``pytest_runtestloop`` already runs under ``unrecorded()`` (``_quiet``), and the inherited
    recording is stopped before the own cone is read."""
    audits = {t for t in reads if t.startswith("tests/rd/test_audit.py::")}
    checked = 0
    with _store(captured) as store:
        for test_id, found in reads.items():
            data = [(type(p).__name__, p.path) for p in [*store.pins(test_id), *store.own_pins(test_id)]
                    if isinstance(p, (memo.DataPin, memo.ListingPin)) and str(p.path).startswith("pkg")]
            if test_id in audits:
                assert data, f"{test_id}: the audit's own read is missing (the control)"
            elif any(path.startswith("pkg/") for path, _, _ in found):
                assert not data, (test_id, data)
                checked += 1
    print(f"F3: {checked} resolving tests checked, {len(audits)} audit controls")
    assert checked >= 20 and len(audits) == 4
