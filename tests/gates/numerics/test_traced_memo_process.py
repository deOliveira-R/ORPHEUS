r"""Step 3 of #405 P3: the key, the process boundary and the store, end to end through the memo.

Spec ``.claude/plans/reference_p3_spec.md`` §1.3, gates M3.1–M3.15. Every generator is a synthetic one in
``tmp_path`` (``_traced_memo_synthetic``); every edit is to that copy; the cache root is the test's own. The
number of generations is read from the generator's own counter file, never from the memo (X3).
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import textwrap

import numpy as np
import pytest

from . import _traced_memo_api as api
from ._traced_memo_synthetic import Package

pytestmark = pytest.mark.foundation


@pytest.fixture
def package(tmp_path):
    pkg = Package(tmp_path)
    with api.cache_root(tmp_path / "cache"):
        yield pkg
    pkg.close()


def _entry_files(tmp_path):
    return api.entry_files(tmp_path / "cache")


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_validation.py::test_m2_1_validation_misses_exactly_when_something_that_ran_changed', 'tests/gates/numerics/test_traced_memo_validation.py::test_m2_6_floats_cross_the_payload_bit_for_bit', 'tests/gates/numerics/test_traced_memo_validation.py::test_m2_7_arrays_keep_dtype_shape_and_bits_and_load_fresh_and_read_only')
def test_m3_1_a_miss_generates_once_and_a_hit_never(package, tmp_path):
    """M3.1: the first call generates in one child (1 generation) and writes one entry; the second is a ``Hit``
    with 0 generations; both return the same float, bit for bit, as the function called directly here."""
    alpha = package.module("alpha")
    assert api.verdict_kind(alpha.generate.lookup(1.5)) == "Absent"
    first = alpha.generate(1.5)
    assert package.generations("generate") == 1 and len(_entry_files(tmp_path)) == 1
    assert api.verdict_kind(alpha.generate.lookup(1.5)) == "Hit"
    second = alpha.generate(1.5)
    assert package.generations("generate") == 1, "a hit ran the generator"
    direct = alpha.generate.__wrapped__(1.5)
    assert first.hex() == second.hex() == direct.hex()


#: (row, old, new, a miss?)
MEMO_WITNESSES = [
    ("traced-body", "return x * SCALE", "return x * SCALE + 1.0", True),
    ("constant-read", "SCALE = 3.0", "SCALE = 4.0", True),
    ("untraced-body", "return x - 1.0", "return x - 5.0", False),
    ("docstring", '"""Doc of a traced helper."""', '"""Changed doc."""', False),
]


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_validation.py::test_m2_1_validation_misses_exactly_when_something_that_ran_changed', 'tests/gates/numerics/test_traced_memo_process.py::test_m3_1_a_miss_generates_once_and_a_hit_never')
@pytest.mark.parametrize("row,old,new,miss", MEMO_WITNESSES, ids=[w[0] for w in MEMO_WITNESSES])
def test_m3_2_an_edit_to_what_ran_regenerates_and_nothing_else_does(package, row, old, new, miss):
    """M3.2, the X1 witnesses through the memo: after an edit to a traced def or to a constant it reads, the
    lookup is ``Stale`` and the next call generates again and returns the NEW answer; after an edit to an
    untraced def or a docstring, the lookup is ``Hit`` and the call generates nothing."""
    alpha = package.module("alpha")
    before = alpha.generate(1.5)
    package.edit("alpha", old, new)
    assert api.verdict_kind(alpha.generate.lookup(1.5)) == ("Stale" if miss else "Hit"), row
    after = alpha.generate(1.5)
    assert package.generations("generate") == (2 if miss else 1), row
    assert (after != before) is (row in ("traced-body", "constant-read")), (row, before, after)


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_process.py::test_m3_1_a_miss_generates_once_and_a_hit_never')
def test_m3_3_the_signature_binds_one_call_to_one_key(package):
    """M3.3: an omitted default, the default spelled explicitly, a positional argument spelled by keyword, and
    keyword-only defaults are ONE key and ONE generation (the census: without the binding, ``Billiard``'s
    omission of ``None`` settings and a direct caller's explicit default are two entries); another value is
    another key."""
    alpha = package.module("alpha")
    spellings = [((4,), {}), ((4, 2.0), {}), ((), {"n": 4, "scale": 2.0, "offset": 0.0}), ((4,), {"offset": 0.0})]
    keys = {alpha.solve.key(*a, **k) for a, k in spellings}
    assert len(keys) == 1, keys
    for a, k in spellings:
        alpha.solve(*a, **k)
    assert package.generations("solve") == 1
    assert alpha.solve.key(4, 2.5) not in keys and alpha.solve.key(5) not in keys


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_process.py::test_m3_3_the_signature_binds_one_call_to_one_key')
def test_m3_4_the_key_separates_what_content_identifies_and_the_function_tells_apart(package):
    """M3.4: ``content_digest`` identifies a float32 array with its float64 twin and ``8`` with ``8.0``
    (``[M]`` 2026-10-04 on the main tree, ``probes/probe_key_hazards.log``), and a function can tell each pair
    apart. The key does not: float32 and float64 arguments are two entries with two different answers, and
    ``solve(8.0)`` raises as the direct call does instead of being served ``solve(8)``'s entry."""
    alpha = package.module("alpha")
    x64 = np.array([1.0, 2.0, 4.0])
    x32 = x64.astype(np.float32)
    assert alpha.divide.key(x64) != alpha.divide.key(x32)
    assert alpha.divide(x64) == alpha.divide.__wrapped__(x64)
    assert alpha.divide(x32) == alpha.divide.__wrapped__(x32)
    assert alpha.divide(x64) != alpha.divide(x32), "the precondition: the function tells the dtypes apart"
    alpha.solve(8)
    with pytest.raises(TypeError, match="cannot be interpreted as an integer"):
        alpha.solve(8.0)


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_process.py::test_m3_3_the_signature_binds_one_call_to_one_key')
def test_m3_5_a_declared_canonical_form_is_both_the_key_and_what_the_child_receives(package):
    """M3.5: ``total`` declares ``x`` canonical as a float64 array: a list and an array of equal values share
    one key and one generation, and the generator receives the canonical form (its answer carries +0.5 only for
    an ndarray), so the shared entry is the answer to what was keyed, never to the caller's spelling."""
    alpha = package.module("alpha")
    assert alpha.total.key([1, 2, 3]) == alpha.total.key(np.array([1.0, 2.0, 3.0]))
    assert alpha.total([1, 2, 3]) == 6.5
    assert alpha.total(np.array([1.0, 2.0, 3.0])) == 6.5
    assert package.generations("total") == 1


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_process.py::test_m3_3_the_signature_binds_one_call_to_one_key')
def test_m3_6_an_argument_with_no_content_is_refused_by_name(package):
    """M3.6: an argument with no content identity cannot be keyed: ``ContentlessError`` naming the parameter,
    before any generation."""
    alpha = package.module("alpha")
    with pytest.raises(api.name("ContentlessError") if hasattr(api.module(), "ContentlessError") else TypeError, match="'x'"):
        alpha.generate(object())
    assert package.generations() == 0


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_process.py::test_m3_3_the_signature_binds_one_call_to_one_key')
def test_m3_7_the_function_and_the_platform_are_in_the_key(package, monkeypatch):
    """M3.7: two functions given the same arguments are two keys; a different platform tag is a different key
    (an entry restored from another machine is never looked up here)."""
    alpha, beta = package.module("alpha"), package.module("beta")
    assert alpha.generate.key(2.0) != beta.inner.key(2.0)
    before = alpha.generate.key(2.0)
    monkeypatch.setattr(api.module(), "platform_tag", lambda: "another-machine")
    assert alpha.generate.key(2.0) != before


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_process.py::test_m3_1_a_miss_generates_once_and_a_hit_never')
def test_m3_8_the_child_is_born_clean_and_imports_from_the_parents_path(package, tmp_path):
    """M3.8: the package root is on this process's ``sys.path`` only (not ``PYTHONPATH``, not the working
    directory), and the child still imports it: it inherits the path, so a worktree parent never spawns a
    main-tree child (lessons L99). The entry's manifest holds the generator module's skeleton and
    ``orpheus/numerics/content.py``'s: tracing began before the first first-party import."""
    assert str(package.root) not in os.environ.get("PYTHONPATH", "") and os.getcwd() != str(package.root)
    alpha = package.module("alpha")
    alpha.generate(1.5)
    (entry,) = _entry_files(tmp_path)
    modules = {rel for rel, _ in api.manifest_rows(api.read_entry(entry), "ModulePin")}
    assert package.relpath("alpha") in modules, modules
    assert "orpheus/numerics/content.py" in modules, sorted(modules)[:10]


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_process.py::test_m3_1_a_miss_generates_once_and_a_hit_never')
def test_m3_9_arguments_cross_through_their_constructor(package, tmp_path):
    """M3.9, finding F1: a frozen dataclass argument is rebuilt in the child by its constructor, so its
    ``__post_init__`` is in the entry's manifest, and an edit to it makes the entry stale. (Default pickling
    restores the fields without running ``__post_init__``: ``[M]`` census F1, 2 of 2 workloads.)"""
    alpha = package.module("alpha")
    alpha.read_spec(alpha.Spec(3))
    (entry,) = _entry_files(tmp_path)
    functions = {q for _, q, *_ in api.manifest_rows(api.read_entry(entry), "DefPin")}
    assert "Spec.__post_init__" in functions, sorted(functions)
    package.edit("alpha", 'raise ValueError("Spec: n is a count")', 'raise ValueError("Spec: n counts")')
    assert api.verdict_kind(alpha.read_spec.lookup(alpha.Spec(3))) == "Stale"


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_process.py::test_m3_1_a_miss_generates_once_and_a_hit_never')
def test_m3_10_a_patch_in_the_parent_never_reaches_an_entry_and_bypass_reads_and_writes_nothing(package, tmp_path, monkeypatch):
    """M3.10: (a) a parent-side patch of a helper does not reach the child: the entry and the answer are
    honest; (b) under ``bypass()`` the patched function runs HERE and returns the decoy, a warm honest entry is
    NOT served (nothing is read), and no entry is written; (c) the control: outside the bypass the warm entry is
    a ``Hit`` and serves the honest value."""
    alpha = package.module("alpha")
    honest = alpha.generate.__wrapped__(1.5)
    monkeypatch.setattr(alpha, "helper", lambda x: -1000.0)
    decoy = alpha.generate.__wrapped__(1.5)
    assert decoy != honest
    assert alpha.generate(1.5) == honest  # (a)
    entries = _entry_files(tmp_path)
    with api.bypass():
        assert alpha.generate(1.5) == decoy  # (b): the warm entry was not read
        assert alpha.generate(2.5) == alpha.generate.__wrapped__(2.5)
    assert _entry_files(tmp_path) == entries, "the bypass wrote an entry"
    assert api.verdict_kind(alpha.generate.lookup(2.5)) == "Absent"
    assert alpha.generate(1.5) == honest  # (c)


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_process.py::test_m3_1_a_miss_generates_once_and_a_hit_never')
def test_m3_11_a_raising_generator_writes_nothing_and_raises_its_own_error(package, tmp_path):
    """M3.11: the child's exception reaches the caller with its type and message; no entry is written, so the
    next call generates (and raises) again: a refusal is never cached."""
    alpha = package.module("alpha")
    for attempt in (1, 2):
        with pytest.raises(ValueError, match="boom at 2.0"):
            alpha.boom(2.0)
        assert package.generations("boom") == attempt
    assert _entry_files(tmp_path) == []


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_manifest.py::test_m1_3c_signature_scopes_are_pinned_by_the_skeleton_not_by_a_body', 'tests/gates/numerics/test_traced_memo_process.py::test_m3_1_a_miss_generates_once_and_a_hit_never')
def test_m3_12a_a_parent_entry_pins_its_child_by_reference(package, tmp_path):
    """M3.12 (a): ``outer`` calls the memo ``inner``. The outer entry's manifest records the inner entry (its
    key and payload digest) and does NOT list ``inner``'s def, which ran only in the grandchild: the parent
    pins the child by reference, not by re-tracing it."""
    beta = package.module("beta")
    beta.outer(2.0)
    assert package.generations("outer") == 1 and package.generations("inner") == 1
    outer = api.read_entry(next(f for f in _entry_files(tmp_path) if api.entry_function(f).endswith(":outer")))
    assert [c[0].endswith(":inner") for c in api.manifest_rows(outer, "ChildPin")] == [True]
    assert "inner" not in {q for _, q, *_ in api.manifest_rows(outer, "DefPin")}


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_process.py::test_m3_12a_a_parent_entry_pins_its_child_by_reference')
@pytest.mark.parametrize("warm_child", [False, True], ids=["cold-child", "warm-child"])
def test_m3_12b_a_mutation_only_a_child_traced_makes_the_parent_stale(package, warm_child):
    """M3.12 (b), recursive validation: an edit to ``inner``'s body, which only the CHILD entry traced, makes
    the parent's lookup ``Stale`` with a reason naming the child. Both when the child was generated inside the
    parent's generation (cold) and when the parent's generation found it warm (a hit inside a generating
    process must still record the reference)."""
    beta = package.module("beta")
    if warm_child:
        beta.inner(2.0)
    beta.outer(2.0)
    assert api.verdict_kind(beta.outer.lookup(2.0)) == "Hit"
    package.edit("beta", "return x * FACTOR", "return x * FACTOR + 0.0")
    verdict = beta.outer.lookup(2.0)
    assert api.verdict_kind(verdict) == "Stale" and any("child" in r for r in verdict.reasons), verdict


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_process.py::test_m3_12b_a_mutation_only_a_child_traced_makes_the_parent_stale')
def test_m3_12c_a_child_regenerated_with_another_answer_makes_the_parent_stale(package):
    """M3.12 (c): the child entry is regenerated after an edit that changes its answer; the parent's recorded
    payload digest no longer matches, and the reason says so (the child is a valid ``Hit``, so only the digest
    comparison can see this)."""
    beta = package.module("beta")
    beta.outer(2.0)
    package.edit("beta", "return x * FACTOR", "return x * FACTOR * 2.0")
    beta.inner(2.0)
    verdict = beta.outer.lookup(2.0)
    assert api.verdict_kind(verdict) == "Stale" and any("payload" in r for r in verdict.reasons), verdict


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_process.py::test_m3_12a_a_parent_entry_pins_its_child_by_reference')
def test_m3_12d_two_calls_of_one_child_in_one_generation_are_one_entry(package, tmp_path):
    """M3.12 (d): ``outer_twice`` calls ``inner(x)`` twice: one inner generation, one recorded reference."""
    beta = package.module("beta")
    beta.outer_twice(2.0)
    assert package.generations("inner") == 1
    outer = api.read_entry(next(f for f in _entry_files(tmp_path) if api.entry_function(f).endswith(":outer_twice")))
    assert len(api.manifest_rows(outer, "ChildPin")) == 1


#: (row, the corruption)
CORRUPTIONS = [
    ("byte-flipped", lambda f: _flip(f, f.read_bytes().index(b'"payload_digest"') + 3)),
    ("truncated", lambda f: f.write_bytes(f.read_bytes()[:30])),
    ("entry-json-truncated", lambda f: _edit_json(f, None)),
    ("payload-value-edited", lambda f: _edit_json(f, lambda e: e["payload"]["fields"]["k"].update(f=float.hex(9.0)))),
]


def _flip(path, offset):
    data = bytearray(path.read_bytes())
    data[offset] ^= 0xFF
    path.write_bytes(bytes(data))


def _edit_json(path, change):
    """Rewrite the entry's JSON member: ``change`` edits it in place; ``None`` truncates it."""
    entry = api.read_entry(path)
    if change is None:
        import numpy as np

        with np.load(path, allow_pickle=False) as npz:
            arrays = {name: npz[name] for name in npz.files if name != "__entry__"}
        with open(path, "wb") as handle:
            np.savez(handle, **arrays, __entry__=np.frombuffer(json.dumps(entry).encode()[:50], np.uint8))
        return
    change(entry)
    api.write_entry(path, entry)


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_process.py::test_m3_1_a_miss_generates_once_and_a_hit_never')
@pytest.mark.parametrize("row,corrupt", CORRUPTIONS, ids=[c[0] for c in CORRUPTIONS])
def test_m3_13_a_corrupted_entry_is_refused_and_regenerated_never_served(package, tmp_path, row, corrupt):
    """M3.13: a flipped byte, a truncated file, a truncated entry JSON, and a payload value edited in place
    (the JSON parses; only the digest can see it) are each ``Corrupt``; the next call generates again and
    returns the honest answer. (A deleted entry is ``Absent``: one file per entry has no half to lose.)"""
    alpha = package.module("alpha")
    honest = alpha.solve(5)
    (entry,) = _entry_files(tmp_path)
    corrupt(entry)
    assert api.verdict_kind(alpha.solve.lookup(5)) == "Corrupt", row
    again = alpha.solve(5)
    assert package.generations("solve") == 2
    assert again.k == honest.k and np.array_equal(again.field, honest.field)


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_process.py::test_m3_1_a_miss_generates_once_and_a_hit_never')
def test_m3_13b_an_edited_manifest_is_stale(package, tmp_path):
    """M3.13 (b): a manifest edited on disk (a function digest replaced) is ``Stale``, never served."""
    alpha = package.module("alpha")
    alpha.solve(5)
    (entry,) = _entry_files(tmp_path)
    _edit_json(entry, lambda e: next(row for row in e["manifest"] if row[0] == "DefPin").__setitem__(4, "0" * 64))
    assert api.verdict_kind(alpha.solve.lookup(5)) == "Stale"


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_validation.py::test_m2_7_arrays_keep_dtype_shape_and_bits_and_load_fresh_and_read_only', 'tests/gates/numerics/test_traced_memo_process.py::test_m3_1_a_miss_generates_once_and_a_hit_never')
def test_m3_14_a_served_payload_is_fresh_read_only_and_bit_identical(package):
    """M3.14: the miss's return and every hit's return are the same bits as the direct call; each served
    array is read-only and shares no memory with another load (a consumer that writes into a hit raises)."""
    alpha = package.module("alpha")
    miss = alpha.solve(6, 1.5, offset=0.25)
    hit = alpha.solve(6, 1.5, offset=0.25)
    direct = alpha.solve.__wrapped__(6, 1.5, offset=0.25)
    assert miss.k.hex() == hit.k.hex() == direct.k.hex()
    assert miss.field.tobytes() == hit.field.tobytes() == direct.field.tobytes()
    assert type(hit) is alpha.Result and hit.n == 6 and hit.ok is True
    assert not hit.field.flags.writeable and not np.shares_memory(hit.field, miss.field)
    with pytest.raises(ValueError, match="read-only"):
        hit.field[0] = 0.0


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_process.py::test_m3_1_a_miss_generates_once_and_a_hit_never')
def test_m3_15_two_processes_racing_on_one_miss_leave_one_valid_entry(package, tmp_path):
    """M3.15: two parent processes miss the same key at once; both return the same answer and the entry left
    behind is a ``Hit`` (the write is atomic: a reader never sees a half-written entry)."""
    script = tmp_path / "race.py"
    script.write_text(textwrap.dedent(f'''
        import sys
        sys.path[:0] = {[str(package.root), os.getcwd()]!r}
        from {api.MODULE} import cache_root
        from {package.name}.alpha import solve
        with cache_root({str(tmp_path / "cache")!r}):
            print("K", solve(7).k.hex())
    '''))
    procs = [subprocess.Popen([sys.executable, "-O", str(script)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
             for _ in range(2)]
    outs = [p.communicate(timeout=120) for p in procs]
    answers = {line for out, _ in outs for line in out.split() if line != "K"}
    assert len(answers) == 1, outs
    alpha = package.module("alpha")
    assert api.verdict_kind(alpha.solve.lookup(7)) == "Hit"


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_process.py::test_m3_1_a_miss_generates_once_and_a_hit_never')
def test_m3_16_the_generating_process_sees_no_orpheus_environment(package, monkeypatch):
    """M3.16, the runtime backstop of M5.1: the child's environment carries no ``ORPHEUS_*`` variable, so a
    generator that reads one answers its default under the memo whatever the caller set, and a withdrawn
    generator (``ORPHEUS_RUN_WITHDRAWN``) refuses in the child rather than writing an entry a later caller
    without the switch would be served. The control: under ``bypass()`` the same call sees the variable."""
    alpha = package.module("alpha")
    monkeypatch.setenv("ORPHEUS_P3_PROBE", "set-by-the-caller")
    monkeypatch.setenv("P3_PROBE_OTHER", "kept")
    assert alpha.environment_probe("ORPHEUS_P3_PROBE") == "<unset>"
    assert alpha.environment_probe("P3_PROBE_OTHER") == "kept"
    with api.bypass():
        assert alpha.environment_probe("ORPHEUS_P3_PROBE") == "set-by-the-caller"
