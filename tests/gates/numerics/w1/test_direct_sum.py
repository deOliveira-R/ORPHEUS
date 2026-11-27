r"""G05a, G05b, G05e, G05f, G05h, G05i, G05k, G05n, G05o and G05s: the direct sum of spaces (N05, spec 6.6 and 7.6).

One direct sum, ``DirectSumSpace`` (``numerics/direct_sum.py``), replaces ``FullFieldSpace`` and ``CoupledSpace``. The rows that build only ``numerics`` objects live here; the rows
that build a harness problem live in ``tests/gates/cross_method/w1/test_direct_sum_on_problems.py``. Every row names the input that reddens it (instrument doctrine X1) in its
own docstring; the mutations are written for the reviewers to run on a copy and are not run by this file.

Each name that exists only after the unit is read at call time (``_direct_sum_probe``), so on the base tip a row fails in its body: "import-only" below means the body's first line
raises ``ModuleNotFoundError`` for ``orpheus.numerics.direct_sum`` and no other red is available at the base tip. Rows with a second, behavioural red say so.

Host exactness (CI runs ubuntu with numpy 2.5.5, this host macOS arm64 with Accelerate and numpy 2.4.4): every row here is exact on any host. The metric rows compare a
diagonal metric applied by elementwise products in a fixed slot order, a copy, a concatenation or an index; the two rows that compare a reduction (the adjoint law and
``inner_product`` of the sum against the sum of the summands') state a relative tolerance of 1e-13 taken from the reduction depth (at most 24 terms per summand, so a few ulp), and no
row reads the sign a decomposition gives a zero eigenvalue.
"""

from __future__ import annotations

import ast

import numpy as np
import pytest

from tests.gates.numerics.w1 import _direct_sum_probe as dp
from tests.gates.numerics.w1 import _tensor_probe as tp

pytestmark = [
    pytest.mark.foundation,
    pytest.mark.rests_on(
        "tests/gates/numerics/w1/test_tensor_legality.py",
        "tests/gates/numerics/w1/test_one_tensor_type.py",
        "tests/gates/numerics/w1/test_coproduct_layout.py",
    ),
]

RETIRED_CLASSES = frozenset({"FullFieldSpace", "CoupledSpace"})
RETIRED_CONSTRUCTORS = frozenset({"FullFieldSpace", "CoupledSpace", "from_blocks", "from_systems"})
RETIRED_ATTRIBUTES = frozenset({"interior_space", "trace_space", "system_slices", "zeros_factory", "dual_zeros_factory"})


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def _a():
    return dp.space_of(dp.A_ARGS)


def _b():
    return dp.space_of(dp.B_ARGS)


def _c():
    return dp.space_of(dp.C_ARGS)


# ──────────────────────────────── G05a ────────────────────────────────


def census_one_sum(sources: list[tuple[str, str]]) -> tuple[list[str], list[str], list[str]]:
    """``(files defining a retired class, files calling a retired constructor, files defining DirectSumSpace)`` over ``(path, source)`` pairs."""
    defining = [p for p, s in sources if RETIRED_CLASSES & set(dp.class_defs(s))]
    calling = [p for p, s in sources if RETIRED_CONSTRUCTORS & set(dp.called_names(s))]
    direct = [p for p, s in sources if "DirectSumSpace" in dp.class_defs(s)]
    return defining, calling, direct


def test_one_direct_sum() -> None:
    """G05a. ``ast`` over ``orpheus/`` (derivations excluded): no ``ClassDef`` named ``FullFieldSpace`` or ``CoupledSpace``, none of the four retired constructors is called (0 of the 12
    production construction sites), and ``DirectSumSpace`` is a ``ClassDef`` (the positive control of the predicate). First red at the base tip: both classes are defined (``full_field_space.py``,
    ``coupled_system.py``) and ``DirectSumSpace`` is not, so the positive control is red as well. Mutation: append ``class FullFieldSpace(FunctionSpace): pass`` to a copy of a module source; the
    planted-source row below reads it. Runtime half (the production shapes) is ``test_direct_sum_on_problems.py::test_the_production_sums_have_the_declared_shape``."""
    planted = "class FullFieldSpace(FunctionSpace):\n    pass\n\ndef f():\n    return CoupledSpace.from_systems(())\n"
    defining, calling, _ = census_one_sum([("planted.py", planted)])
    _require(defining == ["planted.py"] and calling == ["planted.py"], f"the census predicate does not read a planted retired class and constructor: {defining}, {calling}")
    _, _, direct = census_one_sum([("planted.py", "class DirectSumSpace:\n    pass\n")])
    _require(direct == ["planted.py"], "the positive control does not read a planted DirectSumSpace class")
    sources = dp.orpheus_sources()
    _require(len(sources) > 100, f"the census read {len(sources)} modules of orpheus/")
    defining, calling, direct = census_one_sum(sources)
    _require(direct == ["orpheus/numerics/direct_sum.py"], f"DirectSumSpace is defined in {direct}, expected only orpheus/numerics/direct_sum.py")
    _require(not defining, f"a retired class is still defined: {defining}")
    _require(not calling, f"a retired constructor or class is still called in {len(calling)} modules: {calling[:6]}")


# ──────────────────────────────── G05b ────────────────────────────────


def test_labels_and_order_are_identity() -> None:
    """G05b. Three arms on tensor summands ``A != B`` of equal shape (``A`` and ``C``): (i) ``of(x=A, y=B) != of(u=A, v=B)`` (the labels differ); (ii) ``of(x=A, y=B) != of(y=B, x=A)``
    (the labels are the same and the order differs, so the offsets differ); (iii) equal sums built twice from separately built summands are equal and hash equal, and a nested sum is not
    equal to its flattening. First red: import-only. Mutations: (i) drop the labels from ``identity_key`` (measured on the prototype: ``of(x=A, y=B) != of(u=A, v=B)`` reads False, while
    ``of(x=A, y=B) != of(x=B, y=A)`` stays True, so the arm-(i) comparison is the one that reddens); (ii) sort the summands by label in ``identity_key``; (iii) flatten nested sums in ``of``."""
    S = dp.direct_sum()
    a, c = _a(), _c()
    _require(a != c and a.shape == c.shape, "the two summands the row compares must be different spaces of one shape")
    _require(S.of(x=a, y=c) != S.of(u=a, v=c), "(i) two sums over the same summands with different labels are equal")
    _require(S.of(x=a, y=c) != S.of(y=c, x=a), "(ii) two sums over the same labelled summands in two orders are equal")
    first, second = S.of(x=a, y=c), S.of(x=_a(), y=_c())
    _require(first == second and hash(first) == hash(second), "(iii) equal sums built from separately built summands are not equal or do not hash equal")
    b = _b()
    nested, flat = S.of(a=S.of(b=a, c=b), d=c), S.of(b=a, c=b, d=c)
    _require(nested != flat, "(iii) a nested sum is equal to its flattening")


# ──────────────────────────────── G05e ────────────────────────────────

ACCUMULATORS = frozenset({"cursor", "start", "offset"})
OFFSET_HOMES = frozenset({"orpheus/numerics/coproduct_layout.py", "orpheus/numerics/metric.py"})


def offset_loops(sources: list[tuple[str, str]]) -> dict[str, list[int]]:
    return {p: lines for p, s in sources if (lines := dp.accumulator_loops(s, ACCUMULATORS))}


def offset_slices(source: str) -> list[int]:
    """The line of every subscript slice one of whose bounds reads an ``.offset`` attribute (a hand cut of a block's region: ``buf[slot.offset:slot.offset + size]``)."""
    hits = []
    for n in ast.walk(ast.parse(source)):
        if isinstance(n, ast.Slice):
            bounds = [b for b in (n.lower, n.upper) if b is not None]
            if any(isinstance(m, ast.Attribute) and m.attr == "offset" for b in bounds for m in ast.walk(b)):
                hits.append(n.lineno)
    return hits


def offset_reads(source: str) -> list[int]:
    """The line of every read of an ``.offset`` attribute (ruling V2: no consumer reads a block's offset, it reads ``region``, ``slice_of`` or ``slice_view``)."""
    return dp.attribute_reads(source, frozenset({"offset"}))


def block_start_reads(source: str) -> list[int]:
    """The line of every ``.start`` read of a block's region: an attribute ``start`` whose value is a ``slice_of(...)`` call or a ``region`` attribute (ruling X5: the offset under another spelling)."""
    hits = []
    for n in ast.walk(ast.parse(source)):
        if isinstance(n, ast.Attribute) and n.attr == "start":
            v = n.value
            if (isinstance(v, ast.Call) and isinstance(v.func, ast.Attribute) and v.func.attr == "slice_of") or (isinstance(v, ast.Attribute) and v.attr == "region"):
                hits.append(n.lineno)
    return hits


def test_offsets_once() -> None:
    """G05e. ``ast`` over ``orpheus/`` (derivations excluded): the accumulator pattern ``NAME += <expr>``, ``NAME`` in {cursor, start, offset}, inside a ``for`` appears in no file but the layout
    (``numerics/coproduct_layout.py``) and ``numerics/metric.py`` (``FactoredMetric``'s slot walk, exempt by name: a tensor slot offset). The predicate reads this spelling only: a ``np.cumsum`` or
    ``itertools.accumulate`` in ``direct_sum.py`` passes it, and G05d's slice laws carry the rest. Positive controls: a planted ``for`` loop is read, a planted ``np.cumsum`` is not, and the real
    ``metric.py`` slot walk is read (a known member of the shape). First red at the base tip: 6 loops, 3 of them summand loops in ``coupled_system.py`` (2) and ``space.py`` (1), so the row is
    behaviourally red there. Mutation: add ``cursor += size`` in a loop over summand sizes to a copy of any other module (the planted row reads it)."""
    planted = "def f(xs):\n    cursor = 0\n    for x in xs:\n        cursor += x\n    return cursor\n"
    _require(offset_loops([("p.py", planted)]) == {"p.py": [4]}, "the predicate does not read a planted accumulator loop")
    _require(offset_loops([("p.py", "import numpy as np\ndef f(xs):\n    return np.cumsum(xs)\n")]) == {}, "the predicate reads a np.cumsum as a loop")
    found = offset_loops(dp.orpheus_sources())
    _require("orpheus/numerics/metric.py" in found, "the predicate does not read the known slot walk of orpheus/numerics/metric.py")
    stray = {p: v for p, v in found.items() if p not in OFFSET_HOMES}
    _require(not stray, f"{sum(len(v) for v in stray.values())} offset accumulator loops outside the layout and the slot walk: {stray}")
    # The second half of law L13: no consumer slices a block by its offset; it reads ``slot.region``, ``layout.slice_of`` or ``slot.slice_view``.
    planted_slice = "def f(buf, slot):\n    buf[slot.offset : slot.offset + slot.flat_size] = 0\n"
    _require(offset_slices(planted_slice) == [2], "the offset-slice predicate does not read a planted hand cut")
    _require(offset_slices("def f(buf, slot):\n    buf[slot.region] = 0\n") == [], "the offset-slice predicate reads slot.region as a hand cut")
    cuts = {p: v for p, s in dp.orpheus_sources() if (v := offset_slices(s))}
    _require(not cuts, f"{sum(len(v) for v in cuts.values())} slices bounded by a block's .offset: {cuts}")
    # Ruling V2 widens it: no read of ``.offset`` at all outside the layout's own module (an index computed from it, ``start = slot.offset + k * stride``, is the same hand cut with no slice).
    _require(offset_reads("def f(slot):\n    return slot.offset + 3\n") == [2], "the offset-read predicate does not read a planted index computation")
    _require(offset_reads("def f(slot):\n    return slot.region.start\n") == [], "the offset-read predicate reads region.start as an offset read")
    # Ruling X5 widens it again: the same datum spelled ``space.slice_of(K).start`` or ``slot.region.start`` is the offset read, whatever shifts by it.
    _require(block_start_reads("def f(space):\n    return space.slice_of(K).start\n") == [2], "the start-read predicate does not read a planted slice_of(...).start")
    _require(block_start_reads("def f(slot):\n    return slot.region.start + 1\n") == [2], "the start-read predicate does not read a planted region.start")
    _require(block_start_reads("def f(buf):\n    return buf.start\n") == [], "the start-read predicate reads a bare .start")
    layout_home = "orpheus/numerics/coproduct_layout.py"
    _require(bool(offset_reads(next(s for p, s in dp.orpheus_sources() if p == layout_home))), "the offset-read predicate does not read the layout's own module, a known member")
    starts = {p: v for p, s in dp.orpheus_sources() if (v := block_start_reads(s))}
    _require(not starts, f"{sum(len(v) for v in starts.values())} .start reads of a block's region or slice_of in orpheus/: {starts}")
    reads = {p: v for p, s in dp.orpheus_sources() if p != layout_home and (v := offset_reads(s))}
    _require(not reads, f"{sum(len(v) for v in reads.values())} reads of a block's .offset outside the layout's module: {reads}")


# ──────────────────────────────── G05f ────────────────────────────────


def test_dual_summandwise() -> None:
    """G05f, the space half. ``D.dual().labels == D.labels``; ``D.dual().summand(l) == D.summand(l).dual()``; ``D.dual().dual() == D`` and ``is D`` (memoised). The carrier half (the dual's pair is the
    swapped pair: ``D.dual().carrier.zeros is D.carrier.dual_zeros``) reads the production sums and is in ``test_direct_sum_on_problems.py::test_carrier_pair_is_swapped_by_the_dual``. MB4
    (``tests/gates/cross_method/test_member_seam_builders.py``) stays a strict xfail (W1-P6). First red: import-only. Mutations: return the dual's summands in reversed order (the label arm reddens);
    do not memoise the dual (the ``is`` arm reddens)."""
    S = dp.direct_sum()
    d = S.of(x=_a(), y=_b())
    dual = d.dual()
    _require(dual.labels == d.labels, f"the dual's labels {dual.labels} differ from {d.labels}")
    for label in d.labels:
        _require(dual.summand(label) == d.summand(label).dual(), f"the dual's summand {label!r} is not the dual of the summand")
    _require(dual.dual() == d and dual.dual() is d, "the dual of the dual is not the sum itself (memoised)")


# ──────────────────────────────── G05h ────────────────────────────────


def _kernel_bearing() -> object:
    """A tensor space with a nonzero kernel: the second factor of ``signed_pair_space`` takes the quotient by ``Coord(4, (3,))``."""
    return tp.signed_pair_space()


def test_biproduct_laws() -> None:
    """G05h. On a tensor-only sum with a non-trivial diagonal form per summand and seeded vectors: (a) ``pi_l iota_m = delta_lm`` (exact: a copy), ``sum_l iota_l pi_l = id`` (exact) and the adjoint law
    ``<pi_l x, y>_{S_l} = <x, iota_l y>_D`` under the block-diagonal metric, within 1e-13 relative (a reduction of at most 24 terms); the injection and projection are written through ``pack`` and
    ``unpack`` (a mapping label to the part of the summand's shape); (b) nesting: ``pack`` of ``of(a=of(b=X, c=Y), d=Z)`` equals ``pack`` of ``of(b=X, c=Y, d=Z)`` byte for byte while the two spaces
    are unequal; (c) ``inner_product`` is the sum of the summands' (through the typed element of a ``ToyCarrier``), ``rank`` is the sum of the ranks, ``kernel_basis`` is the block diagonal and ``negative``
    the sum of the summands'. First red: import-only. Mutations: (a) couple two summands in ``apply_metric`` (add each part's mean to the next; the adjoint law and the metric reddens); (b) flatten
    nested sums in ``of`` (the unequal arm reddens); (c) take the maximum over summands in ``rank`` (``rank(D)`` is then the larger summand's rank, which differs from the sum as both are non-zero)."""
    S = dp.direct_sum()
    a, b = _a(), _b()
    d = S.of(a=a, b=b)
    rng = np.random.default_rng(5)
    xa, xb = rng.standard_normal(a.shape), rng.standard_normal(b.shape)
    zero = {"a": np.zeros(a.shape), "b": np.zeros(b.shape)}

    def iota(label: str, part: np.ndarray) -> np.ndarray:
        return d.pack({**zero, label: part})

    def pi(label: str, flat: np.ndarray) -> np.ndarray:
        return np.asarray(d.unpack(flat)[label])

    _require(np.array_equal(pi("a", iota("a", xa)), xa) and np.array_equal(pi("b", iota("b", xb)), xb), "(a) pi_l iota_l is not the identity on a summand")
    _require(not pi("b", iota("a", xa)).any() and not pi("a", iota("b", xb)).any(), "(a) pi_m iota_l is not zero across summands")
    v = rng.standard_normal(d.n_coordinates)
    _require(np.array_equal(iota("a", pi("a", v)) + iota("b", pi("b", v)), v), "(a) the sum of iota_l pi_l is not the identity")
    # adjoint: <pi_l x, y>_{S_l} = <x, iota_l y>_D, both inner products through the metric of production
    for label, s, y in (("a", a, xa), ("b", b, xb)):
        lhs = float(np.dot(pi(label, v).reshape(-1), np.asarray(s.apply_metric(y)).reshape(-1)))
        rhs = float(np.dot(v, np.asarray(d.apply_metric(iota(label, y))).reshape(-1)))
        _require(abs(lhs - rhs) <= 1e-13 * max(abs(lhs), abs(rhs), 1.0), f"(a) the adjoint law fails on summand {label!r}: {lhs!r} against {rhs!r}")
    # (b) nesting
    c = _c()
    nested, flat = S.of(a=S.of(b=a, c=b), d=c), S.of(b=a, c=b, d=c)
    xc = rng.standard_normal(c.shape)
    inner = np.concatenate([xa.reshape(-1), xb.reshape(-1)])
    _require(nested != flat, "(b) a nested sum equals its flattening")
    _require(nested.pack({"a": inner, "d": xc}).tobytes() == flat.pack({"b": xa, "c": xb, "d": xc}).tobytes(), "(b) pack of a nested sum differs from pack of the flattened sum")
    # (c) block sums
    pair = (a, b)
    x, y = dp.toy_element(d, 11), dp.toy_element(d, 12)
    expected = sum(s.inner_product(bx, by) for s, bx, by in zip(pair, x.blocks, y.blocks))
    got = d.inner_product(x, y)
    _require(abs(got - expected) <= 1e-13 * max(abs(expected), 1.0), f"(c) inner_product of the sum {got!r} is not the sum of the summands' {expected!r}")
    kb = _kernel_bearing()
    dk = S.of(k=kb, a=a)
    _require(dk.rank == kb.rank + a.rank, f"(c) rank of the sum {dk.rank} is not {kb.rank} + {a.rank}")
    basis = np.asarray(dk.kernel_basis())
    kernel_k, kernel_a = np.asarray(kb.kernel_basis()), np.asarray(a.kernel_basis())
    hand = np.zeros((kb.n_coordinates + a.n_coordinates, kernel_k.shape[1] + kernel_a.shape[1]))
    hand[: kb.n_coordinates, : kernel_k.shape[1]] = kernel_k
    hand[kb.n_coordinates :, kernel_k.shape[1] :] = kernel_a
    _require(basis.shape == hand.shape and kernel_k.shape[1] > 0 and np.array_equal(basis, hand), f"(c) kernel_basis of the sum is not the block diagonal of the summands' (shapes {basis.shape} and {hand.shape})")
    signed = S.of(s=dp.signed_space(), a=a)
    _require(signed.negative == dp.signed_space().negative + a.negative == 1, f"(c) negative of the sum is {signed.negative}, expected the sum of the summands' (1)")


# ──────────────────────────── V1: one door to the blocks ────────────────────────────


def test_a_carrier_without_labels_is_refused_by_every_typed_verb() -> None:
    """Ruling V1 (the build's verify stage, F3). A typed element reaches the sum's verbs through ``_blocks_of``, which reads ``x.labels`` directly: an object with ``blocks`` and ``with_blocks`` and NO ``labels`` is not
    a carrier of a sum, and each typed verb (``apply_metric``, ``apply_inverse_metric``, ``inner_product`` on either operand, ``with_coordinates``, ``member``, ``coordinates_of``) refuses it by name; the same
    blocks labelled in the sum's order are accepted (the positive control), and labelled in the swapped order are refused by label (the old swapped-block defect of the verify probe, which a ``getattr``
    default let through). Mutation: restore ``labels = getattr(x, "labels", None)`` with the ``is not None`` guard in ``_blocks_of`` (the label-less rows and the hasattr arm redden)."""
    S = dp.direct_sum()
    from dataclasses import dataclass
    from typing import Any

    from orpheus.numerics.direct_sum import SumCarrier

    @dataclass(frozen=True)
    class Bare:
        blocks: tuple[Any, ...]

        def with_blocks(self, blocks: tuple[Any, ...]) -> "Bare":
            return Bare(tuple(blocks))

    a, b = _a(), _b()
    d = S.of_pairs((("a", a), ("b", b)), carrier=SumCarrier(zeros=lambda: Bare(tuple(dp.member_of(s, np.zeros(s.shape)) for _, s in d.summands))))
    good = dp.toy_element(d, 1)
    bare = Bare(good.blocks)
    flat = np.ones(d.n_coordinates)
    verbs = {
        "apply_metric": lambda x: d.apply_metric(x),
        "apply_inverse_metric": lambda x: d.apply_inverse_metric(x),
        "inner_product (first)": lambda x: d.inner_product(x, good),
        "inner_product (second)": lambda x: d.inner_product(good, x),
        "with_coordinates": lambda x: d.with_coordinates(x, flat),
    }
    refusals = {**verbs, "coordinates_of": lambda x: d.coordinates_of(x)}  # no positive control: the toy names no space, which coordinates_of reads after the door
    for name, call in verbs.items():
        call(good)  # the positive control: the labelled carrier is accepted by the same verb
    for name, call in refusals.items():
        with pytest.raises(AttributeError, match="labels"):
            call(bare)
    with pytest.raises(AttributeError, match="labels"):
        d.member(flat)
    swapped = dp.ToyCarrier((good.blocks[1], good.blocks[0]), ("b", "a"))
    for name, call in refusals.items():
        with pytest.raises(ValueError, match="admitted by label"):
            call(swapped)


# ──────────────────────────────── G05i ────────────────────────────────


def hand_slices(source: str) -> list[int]:
    """The line of every read of ``n_coordinates`` or ``shape`` on a summand taken by name (``x.interior_space``, ``x.trace_space`` or ``x.summand(..)``): the size a hand cut of a flat chart at a summand reads."""
    hits = []
    for n in ast.walk(ast.parse(source)):
        if isinstance(n, ast.Attribute) and n.attr in {"n_coordinates", "shape"}:
            v = n.value
            if (isinstance(v, ast.Attribute) and v.attr in {"interior_space", "trace_space"}) or (isinstance(v, ast.Call) and isinstance(v.func, ast.Attribute) and v.func.attr == "summand"):
                hits.append(n.lineno)
    return hits


def size_cuts(source: str) -> list[int]:
    """The line of every subscript slice bounded by a ``.values.size`` read, directly or through a name assigned from one in the same function: a flat chart cut at a carrier block's size."""
    hits = []
    for fn in (n for n in ast.walk(ast.parse(source)) if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))):
        def is_size(e: ast.AST) -> bool:
            return isinstance(e, ast.Attribute) and e.attr == "size" and isinstance(e.value, ast.Attribute) and e.value.attr == "values"
        sizes = {t.id for a in ast.walk(fn) if isinstance(a, ast.Assign) and is_size(a.value) for t in a.targets if isinstance(t, ast.Name)}
        for n in ast.walk(fn):
            if isinstance(n, ast.Slice):
                for b in (x for x in (n.lower, n.upper) if x is not None):
                    if any(is_size(m) or (isinstance(m, ast.Name) and m.id in sizes) for m in ast.walk(b)):
                        hits.append(n.lineno)
    return sorted(set(hits))


def test_retired_members_have_successors() -> None:
    """G05i. ``ast`` over ``orpheus/`` (derivations excluded): no attribute read named ``interior_space``, ``trace_space``, ``system_slices``, ``zeros_factory`` or ``dual_zeros_factory`` (the five names with one
    defining class each; they cannot collide with the carrier pair's field names ``zeros`` and ``dual_zeros``); ``act_on_direct_sum_coordinates`` is not defined; no hand cut of a flat chart by a
    summand's ``n_coordinates`` (the three of ``sn/operators/loss_kernel_gauge.py``); every read of ``.dual_zeros`` or ``.zeros`` of a sum goes through ``carrier`` (by receiver text: an attribute
    ``dual_zeros`` read on any receiver but ``carrier``, ``self`` or a name ``carrier``); a runtime arm on a sum built in the gate: ``hasattr`` is False for ``interior_space``, ``trace_space``, ``systems``,
    ``n_systems``, ``system_slices``, ``zeros`` and ``dual_zeros`` (the successors are ``summand``, ``summands``, ``labels``, ``slice_of`` and ``carrier``), and the sum has no ``axes`` or one that raises
    only ``AttributeError`` (the dual record's ``getattr(space, "axes", None)``). Positive controls: the planted-source reads, and the known reads of the tip (22, 17, 0, 7 and 4 nodes of the five names
    in 7, 7, 0, 2 and 1 files of ``orpheus/``) are the base tip's red. First red at the base tip: those reads. Mutation: re-add ``interior_space`` as a property of the sum (the runtime arm reddens)."""
    planted = "x = s.interior_space\ny = s.trace_space.n_coordinates\nz = s.summand('b').n_coordinates\nw = t.dual_zeros\n"
    _require(len(dp.attribute_reads(planted, RETIRED_ATTRIBUTES)) == 2, "the attribute predicate does not read the planted retired names")
    _require(len(hand_slices(planted)) == 2, "the hand-slice predicate does not read the planted cuts")
    sources = dp.orpheus_sources()
    reads = {p: lines for p, s in sources if (lines := dp.attribute_reads(s, RETIRED_ATTRIBUTES))}
    _require(not reads, f"{sum(len(v) for v in reads.values())} reads of retired member names in {len(reads)} files: {sorted(reads)[:6]}")
    defined = [p for p, s in sources if any(isinstance(n, ast.FunctionDef) and n.name == "act_on_direct_sum_coordinates" for n in ast.walk(ast.parse(s)))]
    _require(not defined, f"act_on_direct_sum_coordinates is still defined in {defined}")
    cuts = {p: v for p, s in sources if p != "orpheus/numerics/coproduct_layout.py" and (v := hand_slices(s))}
    _require(not cuts, f"hand cuts of a flat chart at a summand's size: {cuts}")
    planted_cut = "def f(psi, flat):\n    n = psi.interior.values.size\n    return flat[:n], flat[n:]\n"
    _require(size_cuts(planted_cut) == [3], "the size-cut predicate does not read a planted cut at a block's .values.size")
    _require(size_cuts("def f(psi, flat):\n    n = psi.interior.values.shape\n    return flat[:3]\n") == [], "the size-cut predicate reads a constant slice as a cut")
    size_hits = {p: v for p, s in sources if (v := size_cuts(s))}
    _require(not size_hits, f"{sum(len(v) for v in size_hits.values())} flat-chart cuts bounded by a block's .values.size: {size_hits}")
    stray = []
    for p, s in sources:
        for n in ast.walk(ast.parse(s)):
            if isinstance(n, ast.Attribute) and n.attr == "dual_zeros":
                r = n.value
                ok = (isinstance(r, ast.Attribute) and r.attr == "carrier") or (isinstance(r, ast.Name) and r.id in {"self", "carrier"})
                if not ok:
                    stray.append((p, n.lineno))
    _require(not stray, f"reads of .dual_zeros whose receiver is not a carrier: {stray[:6]}")
    d = dp.direct_sum().of(x=_a(), y=_b())
    for name in ("interior_space", "trace_space", "systems", "n_systems", "system_slices", "zeros", "dual_zeros"):
        _require(not hasattr(d, name), f"a direct sum still has the retired member {name!r}")
    for name in ("summand", "summands", "labels", "slice_of", "carrier"):
        _require(hasattr(d, name), f"a direct sum has no member {name!r}")
    try:
        axes = d.axes
    except AttributeError:
        axes = None
    _require(axes is None, "a direct sum defines axes: the dual record's getattr(space, 'axes', None) would read it as an axis-built space")


# ──────────────────────────────── G05k ────────────────────────────────


def test_a_signed_summand_refuses_norm() -> None:
    """G05k. ``negative`` of a sum is the sum of its summands'; a sum with a signed summand (diagonal form ``[1, -1, 2]``) refuses ``norm`` with the message of N03b ("needs a positive form", naming #706),
    whichever summand is signed and at any nesting depth; the all-positive control returns ``sqrt(sum_l <x_l, x_l>)``, the base body over the sum's ``inner_product`` (ON05-22); ``inner_product`` of the signed
    sum returns its signed value (the hand sum of ``sum(g_l x_l^2)``). First red: import-only (today both retired classes override ``negative`` and refuse: 2 of 2 orders). Mutation: delete the override of
    ``negative`` on the sum (the base ``negative`` reads ``self.metric``, which a sum refuses with ``AttributeError``, so ``norm`` raises the wrong error: ``pytest.raises(ValueError)`` reddens)."""
    S = dp.direct_sum()
    signed, a, p = dp.signed_space(), _a(), dp.positive_space()
    for name, d, pair in (
        ("signed first", S.of(s=signed, a=a), (signed, a)),
        ("signed second", S.of(a=a, s=signed), (a, signed)),
        ("nested", S.of(n=S.of(s=signed, a=a), p=p), ((signed, a), p)),
    ):
        _require(d.negative == 1, f"{name}: negative of the sum is {d.negative}, expected 1")
        x = dp.toy_element(d, 3)
        with pytest.raises(ValueError, match="needs a positive form"):
            d.norm(x)
    d = S.of(s=signed, a=a)
    x = dp.toy_element(d, 4)
    hand = float(np.sum(np.array([1.0, -1.0, 2.0]) * x.blocks[0].values.reshape(-1) ** 2) + np.sum(dp.hand_diagonal(dp.A_ARGS) * x.blocks[1].values.reshape(-1) ** 2))
    got = d.inner_product(x, x)
    _require(abs(got - hand) <= 1e-13 * max(abs(hand), 1.0), f"inner_product of the signed sum is {got!r}, the hand sum is {hand!r}")
    positive = S.of(p=p, a=a)
    y = dp.toy_element(positive, 5)
    square = float(np.sum(np.array([1.0, 2.0, 3.0]) * y.blocks[0].values.reshape(-1) ** 2) + np.sum(dp.hand_diagonal(dp.A_ARGS) * y.blocks[1].values.reshape(-1) ** 2))
    _require(abs(positive.norm(y) - np.sqrt(square)) <= 1e-13 * np.sqrt(square), f"norm of the all-positive sum is {positive.norm(y)!r}, the hand value is {np.sqrt(square)!r}")


# ──────────────────────────────── G05n ────────────────────────────────

DOOR_FILES = (
    "orpheus/transport/operators/isotropic_transfer.py",
    "orpheus/transport/operators/lift.py",
    "orpheus/transport/operators/multiplication_operator.py",
    "orpheus/transport/operators/angular_lift.py",
)


def is_direct_sum_door(source: str) -> list[int]:
    """The line of every ``isinstance(x, DirectSumSpace)`` (or ``FullFieldSpace``) test of a source."""
    hits = []
    for n in ast.walk(ast.parse(source)):
        if isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and n.func.id == "isinstance" and len(n.args) == 2:
            names = {m.id for m in ast.walk(n.args[1]) if isinstance(m, ast.Name)} | {m.attr for m in ast.walk(n.args[1]) if isinstance(m, ast.Attribute)}
            if names & {"DirectSumSpace", "FullFieldSpace"}:
                hits.append(n.lineno)
    return hits


def test_bulk_and_trace() -> None:
    """G05n. ``bulk_and_trace(space, owner="x")`` returns the two summands of a full-field sum (labels ``bulk`` and ``trace``); any other space raises ``TypeError`` with the fragment "composite
    DirectSumSpace": a ray instance (``interior``, ``boundary``), a two-summand sum labelled ``full_field`` and ``radial_characteristic`` (the 6 two-summand coupled spaces of the harness are the same case)
    and a tensor space, because acceptance is by the labels and never by arity; the four doors call it and test no ``isinstance`` on a sum. First red: import-only; the 3 pinned lines read the old fragment
    "composite FullFieldSpace" today (``test_bulk_lift.py:179``, ``:216``, ``test_n2n_operator.py:407``). Mutations: return the old fragment; accept any two-summand sum (the arity reading: the
    ``full_field`` and ``radial_characteristic`` arm reddens)."""
    from orpheus.numerics import direct_sum as m

    S = dp.direct_sum()
    a, b = _a(), _b()
    got = m.bulk_and_trace(S.of_pairs(((m.BULK, a), (m.TRACE, b))), owner="x")
    _require(got[0] == a and got[1] == b, "bulk_and_trace does not return the bulk and the trace summands")
    _require((m.BULK, m.TRACE) == ("bulk", "trace"), f"the full-field labels are {(m.BULK, m.TRACE)}")
    for name, space in (("ray instance", S.of(interior=a, boundary=b)), ("coupled labels", S.of(full_field=a, radial_characteristic=b)), ("tensor", a)):
        with pytest.raises(TypeError, match="composite DirectSumSpace"):
            m.bulk_and_trace(space, owner="x")
        del name
    for rel in DOOR_FILES:
        source = (dp.ROOT / rel).read_text()
        _require("bulk_and_trace" in dp.called_names(source), f"{rel} does not call bulk_and_trace")
        _require(not is_direct_sum_door(source), f"{rel} still tests isinstance on a sum at lines {is_direct_sum_door(source)}")
    _require(is_direct_sum_door("if isinstance(x, DirectSumSpace):\n    pass\n") == [1], "the door predicate does not read a planted isinstance door")


# ──────────────────────────────── G05o ────────────────────────────────


def test_inherited_verbs_are_the_sums_own() -> None:
    """G05o. On a sum whose bulk is a modal head (``ff.basis_space(ff.angular_frame(Quadrature.folded_product(2, 4), 2))``) beside a 3-point tensor summand: ``unit`` is the concatenation of the summands'
    units in label order, ``integral`` is the concatenation of the summands' integrals (it raises on ``FullFieldSpace`` today), ``has_coordinate_cone`` is False; on a sum of two nodal tensor summands it is
    True; ``flat`` equals ``apply_metric`` on the flat chart; each leaf-arm axis verb (``axis``, ``without_axis``, ``with_axis``, ``positions_spanned_by``, ``_axis_index``, ``_axis_block``) raises the sum's own
    ``TypeError`` naming ``summand(label)``. The arms with a trace leaf (None, and False before None) are ``test_direct_sum_on_problems.py::test_cone_conjunction_with_a_trace_leaf``. Red on today's tree on
    the production class, by behaviour: ``FullFieldSpace.unit`` is ``np.ones(shape)``, which differs from the concatenation by 11.566 on this input (positive control in the row: the modal head's own
    unit is not all ones), ``integral`` raises ``AttributeError`` on 21 of 21 full fields, and ``has_coordinate_cone`` reads None on the composite where its head reads False; the gate is written on
    ``DirectSumSpace``, where it is import-only. Mutations: inherit the base's ``np.ones(shape)`` for ``unit`` (reddens by 11.566); evaluate None before False in the conjunction (the cone row of the companion file)."""
    S = dp.direct_sum()
    head, tail = dp.modal_head_and_tail()
    head_unit, tail_unit = np.asarray(head.unit).reshape(-1), np.asarray(tail.unit).reshape(-1)
    _require(not np.array_equal(head_unit, np.ones_like(head_unit)), "positive control: the modal head's own unit is all ones, so the row cannot see the base's np.ones")
    d = S.of(bulk=head, trace=tail)
    concat = np.concatenate([head_unit, tail_unit])
    _require(np.array_equal(np.asarray(d.unit).reshape(-1), concat), f"unit of the sum is not the concatenation of the summands' units (max difference {np.max(np.abs(np.asarray(d.unit).reshape(-1) - concat)):.6g})")
    integral = np.concatenate([np.asarray(head.integral.values).reshape(-1), np.asarray(tail.integral.values).reshape(-1)])
    _require(np.array_equal(np.asarray(d.integral.values).reshape(-1), integral), "integral of the sum is not the concatenation of the summands' integrals")
    _require(d.has_coordinate_cone is False, f"has_coordinate_cone of a modal head beside a nodal tensor is {d.has_coordinate_cone!r}, expected False")
    nodal = S.of(a=_a(), b=_b())
    _require(nodal.has_coordinate_cone is True, f"has_coordinate_cone of two nodal tensor summands is {nodal.has_coordinate_cone!r}, expected True")
    density = np.random.default_rng(9).standard_normal(d.n_coordinates)
    _require(np.array_equal(np.asarray(d.flat(density).values).reshape(-1), np.asarray(d.apply_metric(density)).reshape(-1)), "flat of the sum is not apply_metric on the flat chart")
    factor = tp.triple_factors()[0]
    for verb, args in (
        ("axis", ("x",)),
        ("without_axis", ("x",)),
        ("with_axis", (factor,)),
        ("positions_spanned_by", (_a(),)),
        ("_axis_index", ("x",)),
        ("_axis_block", (0,)),
    ):
        with pytest.raises(TypeError, match=r"summand\(label\)"):
            getattr(d, verb)(*args)


# ──────────────────────────────── G05s ────────────────────────────────


def named_toy(source: str) -> list[int]:
    """The line of every class definition, name and import alias spelled ``_ToySpace`` of a source."""
    tree = ast.parse(source)
    lines = [n.lineno for n in ast.walk(tree) if isinstance(n, ast.ClassDef) and n.name == "_ToySpace"]
    lines += [n.lineno for n in ast.walk(tree) if isinstance(n, ast.Name) and n.id == "_ToySpace"]
    lines += [n.lineno for n in ast.walk(tree) if isinstance(n, ast.alias) and n.name == "_ToySpace"]
    return sorted(set(lines))


def isinstance_calls_in(function: ast.FunctionDef) -> list[int]:
    return [n.lineno for n in ast.walk(function) if isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and n.func.id == "isinstance"]


def test_summands_are_function_spaces() -> None:
    """G05s. ``ast`` over ``tests/``: 0 names ``_ToySpace`` (the retirement witness of the double, ON05-19; positive control: the planted class is read, and the same predicate reads the definition at
    ``tests/gates/numerics/test_coupled_operator.py:169`` on the base tip); ``ast`` over ``DirectSumSpace.of_pairs``: no ``isinstance`` call (there is no runtime refusal, by the user's standing ruling of
    2026-10-30); a ``TensorSpace`` summand is accepted and the sum's summands are the ones passed. First red at the base tip: the double is defined and used in ``test_coupled_operator.py``
    (behavioural). Mutations: re-add the class to a copy of the test file; add an ``isinstance`` refusal to ``of_pairs``."""
    _require(named_toy("class _ToySpace:\n    pass\nx = _ToySpace()\n") == [1, 3], "the predicate does not read a planted _ToySpace")
    found = {p: v for p, s in dp.tree_sources("tests") if (v := named_toy(s))}
    _require(not found, f"_ToySpace is still named in {len(found)} test files: {sorted(found)}")
    source = (dp.ROOT / "orpheus/numerics/direct_sum.py").read_text()
    tree = ast.parse(source)
    cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == "DirectSumSpace")
    of_pairs = next(n for n in cls.body if isinstance(n, ast.FunctionDef) and n.name == "of_pairs")
    _require(not isinstance_calls_in(of_pairs), f"of_pairs refuses by isinstance at lines {isinstance_calls_in(of_pairs)}")
    a, b = _a(), _b()
    d = dp.direct_sum().of_pairs((("a", a), ("b", b)))
    _require(d.summand("a") == a and d.summand("b") == b, "of_pairs does not keep the summands passed")
