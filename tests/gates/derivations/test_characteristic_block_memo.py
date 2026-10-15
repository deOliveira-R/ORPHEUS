r"""C1 of #405 P3 step 2c: one group's transport block as a traced memo (``group_transport``).

Spec ``scratch/characteristic_architecture/p3/c1/gate_spec.md``, rows B1-B6. ``group_transport(basis, walls, sigma_t,
resolution, support)`` equals ``LineRule.of(basis, walls, sigma_t, resolution.line_points).transport(support)``;
it is keyed on the CONTENT of its arguments (``PanelBasis``, ``ConcentricPartition`` and ``RigidMotion`` are
``ContentIdentity`` values since C1) and generated in a fresh interpreter whose manifest pins every def it ran.

Every fixture block costs under 0.2 s in process (``[M]`` 2026-10-14, ``ta/probes/p1_cost.log``, ``p2_cyl.log``); a
miss costs about 3 s (the generating interpreter, its imports and the recorder; ``ta/probes/p4_miss_cost.log``), a hit
about 0.013 s. Each row reads through its own ``cache_root(tmp_path)``, so no row reads an entry another row or an
earlier session wrote.
"""
from __future__ import annotations

import ast
import dataclasses
import inspect
import os
import subprocess
import sys
import textwrap
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from orpheus.derivations.continuous.characteristic import LineRule, PanelBasis, TransportResolution, Wall, Walls
from orpheus.derivations.continuous.characteristic import assembly as assembly_module
from orpheus.derivations.continuous.characteristic import system as system_module
from orpheus.derivations.continuous.characteristic.assembly import group_transport
from orpheus.geometry.chart import Chart
from orpheus.geometry.chord import ConcentricPartition
from orpheus.geometry.coord import CoordSystem
from tests.gates.numerics import _traced_memo_api as api

pytestmark = pytest.mark.foundation

REPO = Path(__file__).resolve().parents[3]
_HERE = "tests/gates/derivations/test_characteristic_block_memo.py::"
_MEMO = "tests/gates/numerics/test_traced_memo_process.py::"

_COORD = {"sphere": CoordSystem.SPHERICAL, "cylinder": CoordSystem.CYLINDRICAL, "slab": CoordSystem.CARTESIAN}
_MR3, _SLB3 = (0.0, 0.5, 1.5, 2.0), (0.0, 0.4, 1.5, 2.3)
_MIRROR, _WHITE = (1.0, 0.0), (0.0, 1.0)
#: Dyadic Sigma_t (exact in float32, so a float32 spelling is the same value: B2) and its void-middle twin (B3).
_SIG, _SIG_VOID = (0.625, 1.25, 0.4375), (0.625, 0.0, 0.4375)
#: Sigma_t NOT exact in float32 (the assembly rows' g0): B1's fixtures, so a rounding canonical form moves the block.
_SIG_G0, _SIG_G0_VOID = (0.6, 1.3, 0.45), (0.6, 0.0, 0.45)
#: The fixtures' resolution: four points per piece.
_R4 = TransportResolution(4)


def _require(cond: object, msg: str) -> None:
    if not cond:
        raise AssertionError(msg)


def _basis(chart: str, breakpoints, degree: int = 3, layers: int = 2) -> PanelBasis:
    return PanelBasis.of(ConcentricPartition(Chart(_COORD[chart]), tuple(breakpoints)), degree, layers, 0.4)


def _walls(chart: str, breakpoints, laws) -> Walls:
    n = len(breakpoints) - 1
    ends = (n,) if chart != "slab" and breakpoints[0] == 0.0 else (0, n)
    return Walls(tuple(Wall(k, law[0], law[1], k) for k, law in zip(ends, laws, strict=True)), n, Chart(_COORD[chart]))


def _emitting(sigma_t) -> np.ndarray:
    """The support a caller resolves when no source is posed: the emitting regions, ``sigma_t > 0``."""
    return np.asarray(sigma_t) > 0.0


#: (id, basis, walls, sigma_t, points, support): the three line shapes (impact, impact x polar, cosine), a void region.
_BLOCKS = {
    "sphere3-white": lambda: (_basis("sphere", _MR3), _walls("sphere", _MR3, (_WHITE,)), _SIG_G0, _R4, _emitting(_SIG_G0)),
    "sphere3-void-white": lambda: (_basis("sphere", _MR3), _walls("sphere", _MR3, (_WHITE,)), _SIG_G0_VOID, _R4,
                                   _emitting(_SIG_G0_VOID)),
    "slab3-mirror-white": lambda: (_basis("slab", _SLB3), _walls("slab", _SLB3, (_MIRROR, _WHITE)), _SIG_G0, _R4,
                                   _emitting(_SIG_G0)),
    "cylinder1-white-p1": lambda: (_basis("cylinder", (0.0, 1.0), 1, 0), _walls("cylinder", (0.0, 1.0), (_WHITE,)), (0.7,), TransportResolution(2),
                                   _emitting((0.7,))),
}


def _arrays(value: Any, path: str = "") -> dict[str, Any]:
    """Every array and scalar leaf of a block, by path (``GroupTransport`` -> ``WallCoupling`` -> ``DiffuseWalls``)."""
    if isinstance(value, np.ndarray):
        return {path: value}
    if dataclasses.is_dataclass(value):
        found: dict[str, Any] = {}
        for f in dataclasses.fields(value):
            found |= _arrays(getattr(value, f.name), f"{path}.{f.name}")
        return found
    return {path: value}


# ── B1: a served block is the in-process block ──────────────────────────────


@pytest.mark.parametrize("name", list(_BLOCKS))
@pytest.mark.rests_on(_MEMO + "test_m3_14_a_served_payload_is_fresh_read_only_and_bit_identical")
def test_b1_a_served_block_is_the_in_process_block_bit_for_bit(name: str, tmp_path, monkeypatch) -> None:
    """[B1] On a cold root the first call generates (one interpreter started), the second is a ``Hit`` (none started),
    and the served block equals ``group_transport`` under ``bypass()`` leaf for leaf: every array's dtype, shape and
    bytes, every scalar; the served arrays are read-only; ``.block`` agrees bitwise.

    First red (arm ``canonical-float32``, ``gate_spec.md`` §6): a canonical form for ``sigma_t`` that is not the body's
    own first conversion (float32): the key and the generation see the rounded column, the bypass call the given one,
    and the line part differs. Activation: the leaf count is printed and at least 9 (line, support, four coupling
    arrays, three diffuse-wall arrays)."""
    args = _BLOCKS[name]()
    spawns = api.SpawnCounter(monkeypatch)
    with api.cache_root(tmp_path):
        group_transport(*args)
        cold = spawns.count
        _require(api.verdict_kind(group_transport.lookup(*args)) == "Hit", f"{name}: the entry just written does not validate")
        served = group_transport(*args)
    _require((cold, spawns.count) == (1, 1), f"{name}: interpreters started cold {cold}, after the hit {spawns.count}")
    with api.bypass():
        fresh = group_transport(*args)
    leaves_served, leaves_fresh = _arrays(served), _arrays(fresh)
    _require(set(leaves_served) == set(leaves_fresh), f"{name}: leaves {sorted(leaves_served)} vs {sorted(leaves_fresh)}")
    _require(len(leaves_served) >= 9, f"{name}: activation: {len(leaves_served)} leaves")
    for path, a in leaves_served.items():
        b = leaves_fresh[path]
        if isinstance(a, np.ndarray):
            _require(a.dtype == b.dtype and a.shape == b.shape, f"{name}{path}: {a.dtype}{a.shape} vs {b.dtype}{b.shape}")
            _require(a.tobytes() == b.tobytes(), f"{name}{path}: max|a-b|/max|b| = "
                     f"{np.max(np.abs(a - b)) / max(np.max(np.abs(b)), 1e-300) if a.size else 0.0}")
            _require(not a.flags.writeable, f"{name}{path}: a served array is writeable")
        else:
            _require(type(a) is type(b) and a == b, f"{name}{path}: {a!r} vs {b!r}")
    _require(served.block.tobytes() == fresh.block.tobytes(), f"{name}: the block")
    print(f"B1 {name}: {len(leaves_served)} leaves compared")


# ── B2: one block, one key; every argument is in it ─────────────────────────

#: The reference call's support: every region of ``_SIG`` emits.
_EMIT = _emitting(_SIG)
#: Spellings of ONE call: each must key as the reference call does (the canonical forms, and content identity).
_SPELLINGS = {
    "sigma-as-ndarray": lambda b, w: (b, w, np.asarray(_SIG), _R4, _EMIT),
    "sigma-as-list": lambda b, w: (b, w, list(_SIG), _R4, _EMIT),
    "sigma-as-float32": lambda b, w: (b, w, np.asarray(_SIG, dtype=np.float32), _R4, _EMIT),
    # a numpy integer is a legitimate count the constructor normalises (``operator.index``); the field is typed ``int``
    "resolution-of-a-numpy-int": lambda b, w: (b, w, _SIG, TransportResolution(np.int64(4)), _EMIT),  # pyright: ignore[reportArgumentType]
    "support-as-int-list": lambda b, w: (b, w, _SIG, _R4, [1, 1, 1]),
    "basis-built-again-from-int-breakpoints": lambda b, w: (_basis("sphere", (0, 0.5, 1.5, 2)), w, _SIG, _R4, _EMIT),
    # a numpy integer degree, which the basis normalises (``operator.index``); the field is typed ``int``
    "basis-of-a-numpy-int-degree": lambda b, w: (_basis("sphere", _MR3, degree=np.int64(3)), w, _SIG, _R4, _EMIT),  # pyright: ignore[reportArgumentType]
    "walls-built-again": lambda b, w: (b, _walls("sphere", _MR3, (_WHITE,)), _SIG, _R4, _EMIT),
    # a wall's amplitudes as integers (``Wall`` parses them as reals) and its indices as numpy integers
    "walls-of-integer-amplitudes": lambda b, w: (
        b, Walls((Wall(np.int64(3), 0, 1, np.int64(3)),), 3, Chart(_COORD["sphere"])),  # pyright: ignore[reportArgumentType]
        _SIG, _R4, _EMIT),
}
#: One leg per PARAMETER of ``group_transport``: each moves that argument alone, and the key must move.
_SEPARATIONS = {
    "basis": lambda b, w: dict(basis=_basis("sphere", _MR3, degree=2)),
    "walls": lambda b, w: dict(walls=_walls("sphere", _MR3, (_MIRROR,))),
    "sigma_t": lambda b, w: dict(sigma_t=(0.625, 1.25, float(np.nextafter(0.4375, 1.0)))),
    "resolution": lambda b, w: dict(resolution=TransportResolution(5)),
    "support": lambda b, w: dict(support=np.array([True, False, True])),
}


@pytest.mark.rests_on("tests/gates/derivations/test_characteristic_reference.py::test_c1_equal_fields_are_one_value",
                      "tests/gates/geometry/test_content_identity_geometry.py::test_s5_2_equal_content_is_one_value")
def test_b2_one_call_is_one_key_and_every_argument_is_in_it() -> None:
    """[B2] Seven spellings of one call (``sigma_t`` as a tuple, an array, a list, float32 of exactly representable values;
    the resolution of an ``np.int64``; ``support`` as a list of integers; the basis and the walls built again, the breakpoints as
    integers) give the reference
    call's key. The separation legs are quantified over the SIGNATURE (``inspect.signature(group_transport)``,
    ``instrument-doctrine`` X2 STRUCTURAL-DENOMINATOR): one leg per parameter moves that argument alone and the key
    moves. The line rule's chunk and budget re-order the sum over lines (``test_the_block_does_not_depend_on_the_piece_
    budget``: up to 12.9 eps, not bitwise); they are not parameters, and the entry pins their defaults with
    ``LineRule.of``'s source (the elegance review, 2026-10-15).

    The void-region leg (a ``-0.0`` in sigma_t) and the integer-amplitude wall are the witnesses of the folds added after
    qa's part-two review (arms ``sigma-fold-removed``, ``wall-parse-removed``). First reds: the canonical forms dropped (``sigma-as-ndarray``, ``-list``, ``-float32`` and ``points-as-numpy-int``
    key apart: the exact key tags a tuple, a list and an array apart and an int from an ``np.int64``); the basis
    compared by identity (``ContentlessError`` at the key before C1); a parameter mapped to a constant by its canonical
    form (the ``support`` leg)."""
    b, w = _basis("sphere", _MR3), _walls("sphere", _MR3, (_WHITE,))
    reference = group_transport.key(b, w, _SIG, _R4, _EMIT)
    apart = {name: group_transport.key(*spell(b, w)) for name, spell in _SPELLINGS.items()}
    _require(all(k == reference for k in apart.values()), f"spellings keyed apart: {sorted(n for n, k in apart.items() if k != reference)}")
    void, signed_void = (0.625, 0.0, 0.4375), (0.625, -0.0, 0.4375)     # a void region spelled with a signed zero
    _require(group_transport.key(b, w, void, _R4, _emitting(void)) == group_transport.key(b, w, signed_void, _R4, _emitting(signed_void)),
             "a -0.0 total cross section keys apart from 0.0 (the canonical form must fold it)")
    parameters = set(inspect.signature(group_transport).parameters)
    _require(set(_SEPARATIONS) == parameters, f"separation legs {sorted(_SEPARATIONS)} vs parameters {sorted(parameters)}")
    base = dict(basis=b, walls=w, sigma_t=_SIG, resolution=_R4, support=_EMIT)
    same = [p for p, moved in _SEPARATIONS.items() if group_transport.key(**(base | moved(b, w))) == reference]
    _require(not same, f"arguments that do not move the key: {same}")


# ── B3: the caller resolves the support; the emitting regions are never every region ────


@pytest.mark.rests_on(_HERE + "test_b1_a_served_block_is_the_in_process_block_bit_for_bit")
def test_b3_the_support_is_required_and_the_emitting_regions_never_key_as_every_region() -> None:
    """[B3] ``support`` has no default, in ``group_transport`` and in ``LineRule.transport`` (the user's ruling of
    2026-10-14, ``gate_spec.md`` §3 B3; the elegance review of 2026-10-15 for the transport): the caller resolves it,
    so "no mask given" is neither a value the memo can key nor a second, silent definition of where emission lives
    (``sigma_t > 0`` is not the system's ``emission_support() | source_regions``). On a sphere with a void middle
    region the block of the emitting regions has 32 columns and the block of every region 52: the value law, read in
    process (``bypass()``; B1 is the served-equals-in-process row). The keys: the two masks key apart, and a mask keys
    as its own integer-list spelling.

    First reds: a default restored on either signature (the signature legs); ``_region_mask`` mapping every mask to
    the all-true mask (the value and key legs)."""
    for name, function in (("group_transport", group_transport), ("LineRule.transport", LineRule.transport)):
        default = inspect.signature(function).parameters["support"].default
        _require(default is inspect.Parameter.empty, f"{name}: support defaults to {default!r}")
    b, w = _basis("sphere", _MR3), _walls("sphere", _MR3, (_WHITE,))
    emitting, every = _emitting(_SIG_VOID), np.ones(3, dtype=bool)
    with api.bypass():
        resolved = group_transport(b, w, _SIG_VOID, _R4, emitting)
        everywhere = group_transport(b, w, _SIG_VOID, _R4, every)
    _require(resolved.support.shape == (32,) and everywhere.support.shape == (52,),
             f"value: supports {resolved.support.shape} and {everywhere.support.shape}")
    k_resolved, k_every = (group_transport.key(b, w, _SIG_VOID, _R4, s) for s in (emitting, every))
    _require(k_resolved != k_every, "the emitting regions and every region share a key")
    _require(k_resolved == group_transport.key(b, w, list(_SIG_VOID), _R4, emitting.astype(int)), "a resolved mask keys as itself")
    # two masks of one count, one a permutation of the other: a canonical form that sorts, or counts, the mask keys them
    # as one (qa, 2026-10-15: only B1's void fixture saw that collision)
    permuted = group_transport.key(b, w, _SIG_VOID, _R4, np.array([False, True, True]))
    _require(k_resolved != permuted, "two permutations of one support mask share a key")


# ── B4: the block entry pins the kernel, not its consumers ──────────────────

#: Defs the block runs (``[M]`` on the prototype, ``ta/probes/p5_pins.log``), one per kernel module the block reads.
_KERNEL = {"group_transport", "LineRule.of", "LineRule.transport", "Lines.chunks", "TraversalRule.on_chord",
           "AttenuationRule.of", "ConcentricPartition.chord", "PanelBasis.values", "DiffuseWalls.of",
           "StackedSources.of", "LinePeriod.inflow", "_leggauss"}
#: Defs of the block's consumers, none of which the block may run.
_CONSUMER_PREFIXES = ("GalerkinSystem.", "PointRule.", "CharacteristicDerivation.", "DensePencil.", "RegionCrossSections.",
                      "EmissionSpace.", "angular_flux")


@pytest.mark.rests_on(_HERE + "test_b1_a_served_block_is_the_in_process_block_bit_for_bit",
                      "tests/gates/numerics/test_traced_memo_manifest.py::test_m1_1_and_m1_2_the_digests_move_exactly_with_what_can_change_an_answer")
def test_b4_a_block_entry_pins_the_kernel_and_none_of_its_consumers(tmp_path) -> None:
    """[B4] A sphere block's entry pins the kernel defs it ran (``_KERNEL``, one per kernel module) and no def of a
    consumer (the system, the point rule, the door, the dense pencil, the cross sections), and no child entry: the
    manifest is the decoupling the user asked for (step 2b: "decoupling the impact"). The cylinder's ``polar_rule`` is
    absent on a sphere (the per-chart negative leg, a def the cylinder's block does pin).

    First red (arm ``door-keyed``, the census's C4 fallback): a block computed through the system
    (``GalerkinSystem(...).groups[g]``), which pins ``GalerkinSystem.__post_init__`` and ``RegionCrossSections.of``."""
    with api.cache_root(tmp_path):
        group_transport(*_BLOCKS["sphere3-white"]())
    (entry_file,) = api.entry_files(tmp_path)
    entry = api.read_entry(entry_file)
    functions = {q for _, q, *_ in api.manifest_rows(entry, "DefPin")}
    _require(_KERNEL <= functions, f"kernel defs not pinned: {sorted(_KERNEL - functions)}")
    consumers = sorted(q for q in functions if q.startswith(_CONSUMER_PREFIXES))
    _require(not consumers, f"the block ran consumer defs: {consumers}")
    _require("polar_rule" not in functions, "a sphere block pins the cylinder's polar rule")
    _require(api.manifest_rows(entry, "ChildPin") == [], "a block reads no other memo entry")
    print(f"B4: {len(functions)} defs pinned")


# ── B5: an edit misses exactly where the block ran (on a copy of the tree) ──

from tests.gates.numerics.test_traced_memo_clients import _copy_tree, _neutral_edit  # noqa: E402

_LOOKUP_DRIVER = textwrap.dedent("""
    import os, sys
    copy, cache, action = sys.argv[1:4]
    sys.path[:0] = [copy, sys.argv[4]]
    import orpheus
    assert orpheus.__file__.startswith(copy), orpheus.__file__
    if os.environ.get("C1_ARM") == "validate-off":          # the Stale rows' witness: a memo serving unvalidated entries
        import orpheus.numerics.traced_memo as memo
        memo.validate = lambda manifest, store=None: ()
    from tests.gates.derivations.test_characteristic_block_memo import _BLOCKS
    from tests.gates.numerics.test_traced_memo_clients import _tiny_characteristic_sphere
    from orpheus.derivations.continuous.characteristic.assembly import group_transport
    from orpheus.numerics.observable import PointValue
    from orpheus.numerics.traced_memo import cache_root
    with cache_root(cache):
        args = _BLOCKS["sphere3-white"]()
        reading, observable = _tiny_characteristic_sphere(), PointValue(position=0.4, group=1)
        if action == "generate":
            group_transport(*args)
            reading.evaluate(observable)
        print("VERDICT", type(group_transport.lookup(*args)).__name__, type(type(reading).evaluate.lookup(reading, observable)).__name__)
""")

#: (row, file under orpheus/, anchor, method, the BLOCK's verdict, the READING's verdict) after a neutral body edit. The
#: reading is the tiny sphere's point value, whose entry pins the solve and the blocks as children and runs the point
#: rule itself: every edit below is pinned by it (``Stale``), the activation leg proving each edit is a real pinned
#: edit, so a block ``Hit`` beside it is decoupling and not a no-op. The last row is a RECORD of the declared limit (step
#: 2b): the package ``__init__`` imports every submodule, so a block pins all 13 module skeletons; it flips to ``Hit``
#: if ``__init__`` becomes lazy, and then the flip is the point.
_EDITS = [
    ("kernel-traversal", "derivations/continuous/characteristic/transport.py", "class TraversalRule:", "on_chord", "Stale", "Stale"),
    ("kernel-attenuation", "derivations/continuous/characteristic/attenuation_rule.py", "class AttenuationRule:", "of", "Stale", "Stale"),
    ("kernel-chord", "geometry/chord.py", "class ConcentricPartition(", "chord", "Stale", "Stale"),
    ("kernel-basis", "derivations/continuous/characteristic/basis.py", "class PanelBasis(", "values", "Stale", "Stale"),
    ("kernel-assembly", "derivations/continuous/characteristic/assembly.py", "class LineRule:", "transport", "Stale", "Stale"),
    ("kernel-gauss", "derivations/common/quadrature.py", "def _leggauss(", None, "Stale", "Stale"),
    ("the-memo-itself", "derivations/continuous/characteristic/assembly.py", "def group_transport(", None, "Stale", "Stale"),
    # ``chunk`` and ``budget`` are no parameters of the memo: they order the sum over lines (content), and their only
    # witness is ``LineRule.of``'s signature, pinned through its def and its module's skeleton. Not a neutral edit:
    # a different default is a different answer, so the block must miss.
    ("knob-chunk-default", "derivations/continuous/characteristic/assembly.py", "        chunk: int = 512,\n",
     "replace:        chunk: int = 256,\n", "Stale", "Stale"),
    ("consumer-system", "derivations/continuous/characteristic/system.py", "class GalerkinSystem:", "pencil", "Hit", "Stale"),
    ("consumer-reading", "derivations/continuous/characteristic/reading.py", "class PointRule:", "row", "Hit", "Stale"),
    ("consumer-door", "derivations/continuous/characteristic/reference.py", "class CharacteristicDerivation(", "solve", "Hit", "Stale"),
    ("consumer-pencil", "derivations/common/dense_pencil.py", "    def fundamental(", None, "Hit", "Stale"),
    ("consumer-cross-sections", "derivations/continuous/characteristic/cross_sections.py", "class RegionCrossSections", "of", "Hit", "Stale"),
    ("record-system-skeleton", "derivations/continuous/characteristic/system.py",
     "@dataclass(frozen=True, eq=False)\nclass GalerkinSystem:", "constant-before", "Stale", "Stale"),
]


def _drive(copy: Path, cache: Path, action: str) -> tuple[str, str]:
    run = subprocess.run(
        [sys.executable, "-O", "-c", _LOOKUP_DRIVER, str(copy), str(cache), action, str(REPO)],
        capture_output=True, text=True, cwd=str(copy.parent), env={k: v for k, v in os.environ.items() if k != "PYTHONPATH"},
        timeout=600,
    )
    lines = [line for line in run.stdout.splitlines() if line.startswith("VERDICT ")]
    if not lines:
        raise AssertionError(f"the driver failed:\n{run.stderr[-3000:]}")
    _, block, reading = lines[-1].split(" ")
    return block, reading


@pytest.fixture(scope="module")
def _generated_copy(tmp_path_factory):
    root = tmp_path_factory.mktemp("block_copy")
    copy = _copy_tree(root)
    cache = root / "cache"
    return copy, cache, _drive(copy, cache, "generate")


@pytest.mark.parametrize("row,file,anchor,method,block,reading", _EDITS, ids=[e[0] for e in _EDITS])
@pytest.mark.rests_on("tests/gates/numerics/test_traced_memo_process.py::test_m3_2_an_edit_to_what_ran_regenerates_and_nothing_else_does",
                      _HERE + "test_b4_a_block_entry_pins_the_kernel_and_none_of_its_consumers")
def test_b5_an_edit_to_the_copied_tree_misses_exactly_where_the_block_ran(_generated_copy, row, file, anchor, method,
                                                                          block, reading) -> None:
    """[B5] On ONE copy of ``orpheus/`` whose cache holds the sphere block and the tiny sphere's point reading (the
    activation leg: both ``Hit`` before any edit), each row applies a behaviour-neutral edit, looks both up in a fresh
    interpreter, and restores the file. Kernel bodies make the block ``Stale``; consumer bodies leave it a ``Hit`` while
    the reading, which runs them, goes ``Stale`` (each edit is a pinned edit); a new constant in the system module's
    skeleton makes the block ``Stale`` (RECORD of the declared limit); a changed ``chunk`` default of ``LineRule.of``
    makes it ``Stale`` (the knob's only witness since ``chunk`` and ``budget`` left the memo's parameters).

    First reds: the nine block-Stale rows, a memo serving without validation (arm ``validate-off``,
    the five Hit rows green); the reading column, the same arm. The Hit rows' block column has no in-tree arm: B4's
    consumer leg is the static form of the same claim, and its arm ``door-keyed`` reds it."""
    copy, cache, generated = _generated_copy
    _require(generated == ("Hit", "Hit"), f"activation: before any edit the block and the reading read {generated}")
    path = copy / "orpheus" / file
    original = path.read_bytes()
    try:
        if method is not None and method.startswith("replace:"):
            text = path.read_text()
            _require(text.count(anchor) == 1, f"{anchor!r} occurs {text.count(anchor)} times in {path}")
            path.write_text(text.replace(anchor, method.removeprefix("replace:")))
        elif method == "constant-before":
            text = path.read_text()
            _require(text.count(anchor) == 1, f"{anchor!r} occurs {text.count(anchor)} times in {path}")
            path.write_text(text.replace(anchor, "_NEUTRAL_CONSTANT = None\n\n\n" + anchor))
        else:
            _neutral_edit(path, anchor, method)
        _require(path.read_bytes() != original, f"{row}: the edit did not change the file")
        verdicts = _drive(copy, cache, "lookup")
    finally:
        path.write_bytes(original)
    _require(verdicts == (block, reading), f"{row}: block {verdicts[0]} (expected {block}), reading {verdicts[1]} (expected {reading})")


# ── B6: the system reads every group through the memo ────────────────────────


def _fresh_system():
    """A two-group, three-region white sphere, built anew on every call. Not ``test_characteristic_system._system``: that
    builder is ``lru_cache``d, so two of its calls are one object whose ``groups`` (a ``cached_property``) keep whichever
    route computed them first, which made the served-versus-in-process comparison below one object compared with itself
    (``[M]`` 2026-10-14, arm ``consumer-writes`` green in the battery, red alone: ``gate_spec.md`` §6).

    The middle region's mixture transfers nothing into group 1 (``_SELF_ONLY``), so the groups' supports DIFFER
    ([T, T, T] and [T, F, T]): with every region emitting in every group, feeding group 0's support to every group left
    B6 green (qa, 2026-10-15)."""
    from orpheus.derivations.continuous.characteristic import GalerkinSystem, RegionCrossSections, TransportResolution
    from tests.gates.derivations.test_characteristic_system import _ABS, _SELF_ONLY, _UP2N  # two-group mixtures

    return GalerkinSystem(_basis("sphere", _MR3), _walls("sphere", _MR3, (_WHITE,)),
                          RegionCrossSections.of([_UP2N, _SELF_ONLY, _ABS]), TransportResolution(4))


def _production_calls(attribute: str) -> list[str]:
    """Every call ``<x>.<attribute>(`` in ``orpheus/``, by an AST walk (``file:line``)."""
    found = []
    for path in sorted((REPO / "orpheus").rglob("*.py")):
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and node.func.attr == attribute \
                    and isinstance(node.func.value, ast.Name) and node.func.value.id == "LineRule":
                found.append(f"{path.relative_to(REPO)}:{node.lineno}")
    return found


@pytest.mark.rests_on(_HERE + "test_b1_a_served_block_is_the_in_process_block_bit_for_bit")
def test_b6_the_system_reads_every_group_through_the_block_memo(monkeypatch) -> None:
    """[B6, ROUTE] ``GalerkinSystem.groups`` calls ``group_transport`` (as the system module binds it) once per group, with
    the system's basis and walls, the group's column of total cross sections, the resolution's line points and the
    group's support column; ``LineRule.of`` is called at exactly one site in ``orpheus/``, inside ``group_transport``
    (an AST census, its positive control the site itself).

    First red: today's tree (``groups`` calls ``LineRule.of`` directly: the spy counts 0, and the census finds the site in
    ``system.py``)."""
    calls: list[tuple] = []
    original = system_module.group_transport

    def spy(*args, **kwargs):
        calls.append((args, kwargs))
        return original(*args, **kwargs)

    monkeypatch.setattr(system_module, "group_transport", spy)
    system = _fresh_system()
    with api.bypass():
        groups = system.groups
    _require(len(calls) == system.n_groups == len(groups) >= 2, f"{len(calls)} calls for {system.n_groups} groups")
    support = system.cross_sections.emission_support() | system.source_regions
    _require(not np.array_equal(support[:, 0], support[:, 1]), f"activation: every group has one support {support.tolist()}")
    for g, (args, kwargs) in enumerate(calls):
        bound = inspect.signature(original).bind(*args, **kwargs)
        a = bound.arguments
        _require(a["basis"] is system.basis and a["walls"] is system.walls, f"group {g}: not the system's basis and walls")
        _require(np.array_equal(a["sigma_t"], system.cross_sections.total[:, g]), f"group {g}: sigma_t")
        _require(a["resolution"] is system.resolution, f"group {g}: not the system's resolution")
        _require(np.array_equal(a["support"], support[:, g]), f"group {g}: support")
    sites = _production_calls("of")
    _require(len(sites) == 1, f"LineRule.of( call sites in orpheus/: {sites}")
    body, start = inspect.getsourcelines(assembly_module.group_transport.__wrapped__)
    path, line = sites[0].rsplit(":", 1)
    _require(path == "orpheus/derivations/continuous/characteristic/assembly.py" and start <= int(line) < start + len(body),
             f"the one LineRule.of( site {sites[0]} is not in group_transport's body (lines {start}-{start + len(body) - 1})")
    unresolved = _transport_calls_without_support()
    _require(not unresolved, f".transport( calls in orpheus/ that pass no support: {unresolved}")


def _transport_calls_without_support() -> list[str]:
    """Every ``<x>.transport(...)`` call in ``orpheus/`` that passes neither a positional argument nor ``support=``
    (``file:line``): ``LineRule.transport`` has no default, so each caller resolves its mask. Positive control: the
    call in ``group_transport``'s body is found (it passes one)."""
    found, seen = [], 0
    for path in sorted((REPO / "orpheus").rglob("*.py")):
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and node.func.attr == "transport":
                seen += 1
                if not node.args and not any(k.arg == "support" for k in node.keywords):
                    found.append(f"{path.relative_to(REPO)}:{node.lineno}")
    _require(seen >= 1, "positive control: no .transport( call found in orpheus/")
    return found


@pytest.mark.rests_on(_HERE + "test_b6_the_system_reads_every_group_through_the_block_memo")
def test_b6_a_system_on_served_blocks_is_the_in_process_system(tmp_path) -> None:
    """[B6, the consumers] A two-group system whose blocks are served (a warm root) has the transport matrix and the k
    of the in-process system (``bypass()``), bitwise (``float.hex``): no consumer writes into a block (a served array is
    read-only, an in-process one is not). First red (arm ``consumer-writes``): ``GalerkinSystem.transport`` scaling a
    block's line part in place (``ValueError: assignment destination is read-only`` on the served side only)."""
    with api.cache_root(tmp_path):
        _fresh_system().groups
        warm = _fresh_system()
        _require(not any(g.line.flags.writeable for g in warm.groups), "activation: the warm system's blocks are not served")
        served_k, served_t = warm.pencil.fundamental().k, warm.transport
    with api.bypass():
        fresh = _fresh_system()
        fresh_k, fresh_t = fresh.pencil.fundamental().k, fresh.transport
    _require(served_t.tobytes() == fresh_t.tobytes(), "the transport matrices differ")
    _require(float(served_k).hex() == float(fresh_k).hex(), f"k {served_k!r} vs {fresh_k!r}")
