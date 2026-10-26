r"""The content encoder and the ``ContentIdentity`` mixin (#405 P1 step 5).

Spec: ``.claude/plans/reference_p1_spec.md`` §1.5. This module carries the
gates on the ENCODER, witnessed on values of this module or built from
literals:

* S5.1 — digests and hashes are the same in every process (``PYTHONHASHSEED``
  1 and 2), with a control that the harness sees salting;
* S5.4 — the canonical forms: one double per real scalar (``True == 1 ==
  1.0``), ``-0.0`` is ``+0.0``, NaN and an integer beyond 2**53 refused with
  the path, containers keep the distinctions ``==`` keeps, a sparse matrix is
  the matrix, anything else is ``ContentlessError``;
* S5.6 — the ROUTE gate: one encoder. Rebinding it moves the digest, the hash
  and the derived name of every way a synthetic content value of this module
  reaches a key;
* S5.7 — the RECORD fingerprint of four synthetic values;
* S5.8 — the schema tag (ruling of 2026-10-01): a field added, renamed or
  reordered, a version bump (on a new class and on a live one), a class
  renamed, each changes the digest; a class MOVED to another module does not
  (layer DAG ruling 10, 2026-10-22: the relocation gates are
  ``test_content_identity_relocation.py``).

The censuses of the same contract over every production type (S5.9, the
conformance of every ``ContentIdentity`` class and the population against the
union of the rosters; S5.6 and S5.7 over the production types; the qualname
uniqueness gate) range over packages above ``values`` and live at the root,
``tests/gates/test_content_identity_census.py`` (layer DAG ruling 28,
2026-10-26). The per-type legs of S5.2, S5.3 and S5.4 live in the owning
trees (``tests/gates/{data,geometry,mesh}/test_content_identity_*.py``).
S5.1's subprocess builds production values of five layers from a source
string, which the test-layer gate does not read (its declared blind spot).

First red, measured on ``main`` ``1dc31163``: collection fails with
``ModuleNotFoundError: No module named 'orpheus.values.content'``.
"""

from __future__ import annotations

import dataclasses
import enum
import os
import pickle
import subprocess
import sys
import textwrap
from collections import Counter
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from scipy.sparse import csr_matrix

import orpheus.values.content as content_module
from orpheus.values.content import ContentIdentity, ContentlessError, FrozenMapping, content_digest, encode, name_digest
from tests.gates._content_identity_helpers import require

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/values/test_content_identity.py"
_ROOT = Path(__file__).resolve().parents[3]


# ═════════════════════════════════════════════════════════════════════════════
# S5.1 — the same digests and hashes in every process
# ═════════════════════════════════════════════════════════════════════════════

_SEED_SCRIPT = textwrap.dedent(
    """
    import numpy as np
    import orpheus
    from orpheus.data.materials import Materials
    from orpheus.data.macro_xs.synthetic import get_mixture
    from orpheus.geometry import CoordSystem, StructuredGeometry
    from orpheus.mesh import Mesh2D
    from orpheus.numerics.axis import EnergyAxis
    from orpheus.values.content import content_digest
    from orpheus.numerics.space import FunctionSpace

    from orpheus.derivations.continuous.characteristic import PanelBasis
    from orpheus.geometry.chart import Chart
    from orpheus.geometry.chord import ConcentricPartition

    print("FILE", orpheus.__file__)
    values = {
        "panel basis": PanelBasis.of(ConcentricPartition(Chart(CoordSystem.SPHERICAL), (0.0, 0.5, 2.0)), 3, 2, 0.4),
        "mixture A 4g": get_mixture("A", "4g"),
        "materials": Materials({0: get_mixture("A", "2g"), 3: get_mixture("B", "2g")}),
        "hollow sphere": StructuredGeometry(
            coord=CoordSystem.SPHERICAL, breakpoints=(0.5, 1.0, 2.0), mat_ids=(0, 1),
            face_roles=(1, 0)),
        "mesh2d": Mesh2D(np.array([0.0, 1.0, 2.0]), np.array([0.0, 1.0]), np.array([[0], [1]]),
            face_roles={"xmin": 0, "xmax": 1, "ymin": 0, "ymax": 1}),
    }
    for name, value in values.items():
        print("DIGEST", name, content_digest(value).hex())
        try:
            print("HASH", name, hash(value))
        except TypeError as err:
            print("HASH", name, "UNHASHABLE", err)
    print("NAME", FunctionSpace.of_axes(EnergyAxis.synthetic(2)).name)
    print("CONTROL", hash("a salted str hash"))
    """
)


def _run_seeded(seed: str) -> list[str]:
    env = {**os.environ, "PYTHONHASHSEED": seed, "PYTHONPATH": str(_ROOT)}
    out = subprocess.run(
        [sys.executable, "-O", "-c", _SEED_SCRIPT], cwd=_ROOT, env=env,
        capture_output=True, text=True, timeout=300,
    )
    require(out.returncode == 0, f"seed {seed}: the subprocess failed:\n{out.stderr[-3000:]}")
    return out.stdout.splitlines()


def test_s5_1_digests_and_hashes_are_seed_stable() -> None:
    """Two processes, ``PYTHONHASHSEED`` 1 and 2, print the same digests and
    the same hashes for five values (one per layer), and the same derived
    space name. The CONTROL line (a ``str`` hash) must DIFFER, which proves
    the harness sees salting (X1). ``[M]`` at 1dc31163 ``hash(mixture A 2g)``
    was 2365733758199073423 under seed 1 and -2362898479209405550 under 2."""
    one, two = _run_seeded("1"), _run_seeded("2")
    files = [line.split(" ", 1)[1] for line in one if line.startswith("FILE ")]
    require(files and files[0].startswith(str(_ROOT)), f"the subprocess imported {files} (L22)")
    pick = lambda lines, kind: [l for l in lines if l.startswith(kind + " ")]  # noqa: E731
    for kind, n in (("DIGEST", 5), ("HASH", 5), ("NAME", 1)):
        require(len(pick(one, kind)) == n, f"{kind}: {len(pick(one, kind))} lines, expected {n}")
        require(pick(one, kind) == pick(two, kind), f"{kind} lines differ across seeds")
    unhashable = [l for l in pick(one, "HASH") if " UNHASHABLE " in l]
    require(not unhashable, f"values that do not hash: {unhashable}")
    require(pick(one, "CONTROL") != pick(two, "CONTROL"), "control: str hashes did not differ")


# ═════════════════════════════════════════════════════════════════════════════
# S5.4 — the canonical forms of the encoder
# ═════════════════════════════════════════════════════════════════════════════


class _Colour(enum.Enum):
    RED = "red"


class _Shade(enum.Enum):
    RED = "red"


@dataclass
class _Mutable:
    x: float = 1.0


class TestS54EncoderCanonicalForms:
    """Each row pairs a canonical form with the mutation that breaks it (spec
    §1.5 battery table, arms E1-E12)."""

    @pytest.mark.parametrize(
        "values",
        [
            pytest.param((True, 1, 1.0, np.int64(1), np.float32(1.0), np.float64(1.0)), id="one"),
            pytest.param((0, 0.0, -0.0, np.float64(-0.0), False), id="zero"),
            pytest.param((2**53, float(2**53)), id="2**53"),
        ],
    )
    def test_a_real_scalar_is_one_double(self, values: tuple) -> None:
        encodings = {encode(v) for v in values}
        require(len(encodings) == 1, f"{values} encode as {len(encodings)} values")

    def test_an_array_is_its_values_on_its_shape(self) -> None:
        require(encode(np.array([0.0, -0.0])) == encode(np.array([0.0, 0.0])), "signed zero in an array")
        require(encode(np.array([1, 2])) == encode(np.array([1.0, 2.0])), "an int array is its values")
        require(encode(np.zeros(4)) != encode(np.zeros((2, 2))), "the shape is content")
        require(encode(np.arange(3.0)) == encode(np.arange(3.0)[::1].copy(order="F")), "layout is storage")

    @pytest.mark.parametrize(
        "value,fragment",
        [
            pytest.param(float("nan"), r"^float.*NaN", id="a_NaN_scalar"),
            pytest.param({"a": [1.0, float("nan")]}, r"\['a'\]\[1\].*NaN", id="a_NaN_in_a_container"),
            pytest.param(np.array([1.0, np.nan]), r"ndarray.*NaN", id="a_NaN_array"),
            pytest.param(2**53 + 1, r"2\*\*53", id="an_int_beyond_2**53"),
            pytest.param(np.array([2**53 + 1]), r"2\*\*53", id="an_int_array_beyond_2**53"),
            pytest.param(csr_matrix(np.array([[np.nan, 0.0], [0.0, 1.0]])), r"NaN", id="a_NaN_in_a_CSR"),
        ],
    )
    def test_nan_and_wide_ints_are_refused_with_the_path(self, value: Any, fragment: str) -> None:
        with pytest.raises(ValueError, match=fragment):
            encode(value)

    @pytest.mark.parametrize(
        "a,b",
        [
            pytest.param((1,), [1], id="tuple_vs_list"),
            pytest.param({"1": 0}, {1: 0}, id="str_key_vs_int_key"),
            pytest.param("1", 1, id="str_vs_int"),
            pytest.param(b"a", "a", id="bytes_vs_str"),
            pytest.param(None, (), id="None_vs_empty_tuple"),
            pytest.param("", (), id="empty_str_vs_empty_tuple"),
            # Without the length prefix both read ``S a S S b``: the tags alone
            # cannot delimit a payload that contains a tag byte.
            pytest.param(("a", "Sb"), ("aS", "b"), id="a_split_moved_(length_prefix)"),
            pytest.param(((1,),), (1,), id="nesting"),
            pytest.param(_Colour.RED, _Shade.RED, id="an_Enum's_class"),
            pytest.param(_Colour.RED, "RED", id="an_Enum_vs_its_name"),
            pytest.param({1: 2.0, 2: 1.0}, {1: 1.0, 2: 2.0}, id="values_swapped_between_keys"),
        ],
    )
    def test_distinct_values_encode_distinctly(self, a: Any, b: Any) -> None:
        require(encode(a) != encode(b), f"{a!r} and {b!r} encode alike")

    @pytest.mark.parametrize(
        "a,b",
        [
            pytest.param({"x": 1.0, "y": 2.0}, {"y": 2.0, "x": 1.0}, id="mapping_insertion_order"),
            # 1 and 9 share a slot of an 8-slot table, so the insertion order
            # decides the iteration order ([M] [1, 9] vs [9, 1]); asserted below.
            pytest.param(frozenset([1, 9]), frozenset([9, 1]), id="frozenset"),
            pytest.param({1: "a"}, {np.int64(1): "a"}, id="numpy_integer_key"),
        ],
    )
    def test_equal_values_encode_alike(self, a: Any, b: Any) -> None:
        if isinstance(a, (dict, frozenset)) and len(a) > 1:
            require(list(a) != list(b), "activation: the two iterate in different orders")
        require(a == b, "activation: the pair is equal under ==")
        require(encode(a) == encode(b), f"{a!r} and {b!r} encode differently")

    def test_a_sparse_matrix_is_the_matrix(self) -> None:
        dense = np.array([[1.0, 0.0], [2.0, 3.0]])
        base = csr_matrix(dense)
        stored_zero = csr_matrix((np.array([1.0, 0.0, 2.0, 3.0]), np.array([0, 1, 0, 1]), np.array([0, 2, 4])), shape=(2, 2))
        # Row 1 stores column 0 twice (1 + 1): a CSR built from the raw triple
        # keeps the duplicate (``coo.tocsr()`` would sum it on the way in).
        duplicates = csr_matrix((np.array([1.0, 1.0, 1.0, 3.0]), np.array([0, 0, 0, 1]), np.array([0, 1, 4])), shape=(2, 2))
        int64 = csr_matrix((base.data.copy(), base.indices.astype(np.int64), base.indptr.astype(np.int64)), shape=(2, 2))
        unsorted = csr_matrix((np.array([1.0, 3.0, 2.0]), np.array([0, 1, 0]), np.array([0, 1, 3])), shape=(2, 2))
        require(stored_zero.nnz == 4, "activation: the zero is stored")
        require(not unsorted.has_sorted_indices, "activation: the indices are unsorted")
        require(not duplicates.has_canonical_format, "activation: the duplicates are stored")
        for name, other in (("a stored zero", stored_zero), ("duplicates", duplicates), ("int64 indices", int64), ("unsorted indices", unsorted)):
            require(np.array_equal(other.toarray(), dense), f"activation: {name} is the same matrix")
            require(encode(other) == encode(base), f"{name}: the encoding reads storage")
        require(encode(csr_matrix(np.array([[1.0, 0.0], [2.0, 4.0]]))) != encode(base), "a value is content")
        require(encode(csr_matrix(np.zeros((2, 3)))) != encode(csr_matrix(np.zeros((3, 2)))), "the shape is content")

    @pytest.mark.parametrize(
        "value,fragment",
        [
            pytest.param(lambda: 0, r"function", id="a_function"),
            pytest.param(object(), r"object", id="a_plain_object"),
            pytest.param(np.array([object()]), r"dtype object", id="an_object_array"),
            pytest.param(_Mutable(), r"mutable dataclass", id="a_mutable_dataclass"),
            pytest.param({"k": [lambda: 0]}, r"\['k'\]\[0\].*function", id="the_path_to_it"),
        ],
    )
    def test_anything_else_is_contentless(self, value: Any, fragment: str) -> None:
        require(issubclass(ContentlessError, TypeError), "ContentlessError is a TypeError")
        with pytest.raises(ContentlessError, match=fragment):
            encode(value)

    def test_a_mapping_view_is_not_an_immutable_part(self) -> None:
        """A ``MappingProxyType`` is a read-only VIEW of a dict its caller may
        still hold, so a frozen owner cannot carry one: writing through the
        dict after keying would move the content under a cached digest. First
        red (``[M]`` 2026-10-02, the step-7 elegance review): the frozen
        encoder admitted every mappingproxy, because ``FrozenMapping`` handed
        it one over its own storage; that part is now a private items view
        that encodes byte-identically (the RECORD pins did not move)."""
        from types import MappingProxyType

        from orpheus.values.content import FrozenMapping

        @dataclasses.dataclass(frozen=True, eq=False)
        class _HoldsAView(ContentIdentity):
            items: Any

        with pytest.raises(ContentlessError, match="mappingproxy mapping"):
            content_digest(_HoldsAView(MappingProxyType({"k": 1.0})))
        require(encode(FrozenMapping({"k": 1.0})) == encode(FrozenMapping({"k": 1})), "activation: a frozen mapping encodes")

    def test_the_digest_is_blake2b_256_of_the_encoding(self) -> None:
        import hashlib

        value = {"a": (1.0, [2, 3]), "b": np.arange(3.0)}
        require(content_digest(value) == hashlib.blake2b(encode(value), digest_size=32).digest(), "the digest")


# ═════════════════════════════════════════════════════════════════════════════
# S5.8 — the schema tag
# ═════════════════════════════════════════════════════════════════════════════

_MODULE_A = "tests.gates.numerics._schema_a"


def _schema_class(qualname: str, fields: tuple[str, ...], *, version: int = 1, mixin: bool = True) -> type:
    """A frozen dataclass ``qualname`` of the synthetic module ``_MODULE_A``, with float ``fields``."""
    namespace: dict[str, Any] = {
        # Real types, not strings: a string annotation sends ``dataclass`` to
        # ``sys.modules[_MODULE_A]``, which a synthetic module name is not in.
        "__annotations__": {name: float for name in fields},
        "__module__": _MODULE_A, "__qualname__": qualname,
    }
    if mixin:
        namespace["__content_version__"] = version
    cls = type(qualname, (ContentIdentity,) if mixin else (), namespace)
    return dataclass(frozen=True, eq=not mixin)(cls)


def _digest(cls: type, **values: float) -> bytes:
    return content_digest(cls(**values))


@dataclass(frozen=True, eq=False)
class _Versioned(ContentIdentity):
    """A content class of this module whose ``__content_version__`` a row raises in place."""

    x: float


class TestS58Schema:
    """The key covers the class's schema, so an entry written under an older
    schema misses (ruling of 2026-10-01; no ``__setstate__`` guard). Every
    leg holds the shared field VALUES equal, so only the tag can move."""

    def test_the_tag_is_by_name_not_by_class_object(self) -> None:
        """Control: two class objects with one schema digest alike (a class
        re-created by a reload, or in another process, is the same schema)."""
        a, b = _schema_class("Point", ("x",)), _schema_class("Point", ("x",))
        require(a is not b, "activation: two class objects")
        require(_digest(a, x=1.0) == _digest(b, x=1.0), "one schema, two digests")

    @pytest.mark.parametrize(
        "other,values",
        [
            pytest.param(lambda: _schema_class("Point", ("x", "y")), {"x": 1.0, "y": 0.0}, id="a_field_added"),
            pytest.param(lambda: _schema_class("Point", ("z",)), {"z": 1.0}, id="a_field_renamed"),
            pytest.param(lambda: _schema_class("Point", ("x",), version=2), {"x": 1.0}, id="the_version_bumped"),
            pytest.param(lambda: _schema_class("Spot", ("x",)), {"x": 1.0}, id="the_class_renamed"),
        ],
    )
    def test_a_schema_change_moves_the_digest(self, other: Callable[[], type], values: dict) -> None:
        require(_digest(other(), **values) != _digest(_schema_class("Point", ("x",)), x=1.0), "the digest did not move")

    def test_the_field_order_is_in_the_tag(self) -> None:
        """Equal values in both fields: only the ORDER of the names differs."""
        xy, yx = _schema_class("Pair", ("x", "y")), _schema_class("Pair", ("y", "x"))
        require(_digest(xy, x=1.0, y=1.0) != _digest(yx, y=1.0, x=1.0), "the field order is not in the tag")

    def test_a_plain_frozen_dataclass_is_tagged_too(self) -> None:
        """A law's part (a reemission closure, a source) is content through the
        frozen-dataclass rule, and carries the tag the same way."""
        a = _schema_class("Part", ("x",), mixin=False)
        b = _schema_class("Other", ("x",), mixin=False)
        require(encode(a(x=1.0)) != encode(b(x=1.0)), "two part classes encode alike")

    def test_a_live_version_bump_moves_the_digest(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """On a class that already exists and has already been keyed (not a class
        built with the new version, as the parametrized rows above are), through a
        FRESH instance: the version is read when a digest is computed, never
        cached with the class (a cached digest or a cached tag would mask the
        bump: the activation half of a route gate). Re-posed 2026-10-26 (layer DAG
        ruling 28) from a production ``Mixture`` onto the synthetic ``_Versioned``."""
        before = content_digest(_Versioned(1.0))
        require(content_digest(_Versioned(1.0)) == before, "activation: two fresh builds digest alike")
        monkeypatch.setattr(_Versioned, "__content_version__", 2, raising=False)
        require(content_digest(_Versioned(1.0)) != before, "a version bump did not move the digest")



# ═════════════════════════════════════════════════════════════════════════════
# S5.6 — the ROUTE gate: one encoder
# ═════════════════════════════════════════════════════════════════════════════


# ═════════════════════════════════════════════════════════════════════════════
# S5.7 — the RECORD fingerprint
# ═════════════════════════════════════════════════════════════════════════════


# ═════════════════════════════════════════════════════════════════════════════
# The synthetic content values of S5.6 and S5.7
# ═════════════════════════════════════════════════════════════════════════════
#
# Re-posed 2026-10-26 (layer DAG ruling 28): S5.6 and S5.7 are claims about values' encoder, so they are witnessed
# here on content dataclasses of this module; the same two claims over every production type (the fifteen
# rosters and the derived space names, and the digests of a mixture, a geometry and a material set) are a census
# of the contract over its consumers and live at the root, ``tests/gates/test_content_identity_census.py``.


class _Tint(enum.Enum):
    WARM = "warm"
    COOL = "cool"


@dataclass(frozen=True)
class _Leaf:
    """A part keyed through the frozen-dataclass rule (no mixin)."""

    weight: float
    tint: _Tint


@dataclass(frozen=True, eq=False)
class _Sample(ContentIdentity):
    """A content value whose parts reach each branch of the encoder a production type reaches: a real
    array, a sparse matrix, a frozen-dataclass part holding an ``Enum``, a ``FrozenMapping``, a string and a
    tuple of reals."""

    weights: np.ndarray
    transfer: csr_matrix
    leaf: _Leaf
    table: FrozenMapping
    label: str
    edges: tuple[float, ...]


def _read_only(array: np.ndarray) -> np.ndarray:
    array = np.array(array, dtype=float)
    array.setflags(write=False)
    return array


def _read_only_csr(dense: np.ndarray) -> csr_matrix:
    matrix = csr_matrix(dense)
    for part in (matrix.data, matrix.indices, matrix.indptr):
        part.setflags(write=False)
    return matrix


def _sample() -> _Sample:
    """A fresh ``_Sample`` of literals (IEEE-exact: no libm or LAPACK output enters, so its bytes are
    platform-free, AP38)."""
    return _Sample(
        weights=_read_only(np.array([0.5, 1.25, -0.0, 3.0])),
        transfer=_read_only_csr(np.array([[1.0, 0.0], [0.25, 2.0]])),
        leaf=_Leaf(weight=0.75, tint=_Tint.COOL),
        table=FrozenMapping({"fast": 1.5, "thermal": 2}),
        label="u235",
        edges=(0.0, 0.5, 2.0),
    )


def _synthetic_route_subjects() -> dict[str, Callable[[], object]]:
    """Every way a value of this module reaches a key: the mixin's digest and hash, a frozen-dataclass part,
    a ``FrozenMapping`` (a mixin type overriding ``content_parts``), an ``Enum`` member, a plain container and
    a derived name (``name_digest``, the spelling a space name carries). Each reading builds its subject fresh."""
    return {
        "digest _Sample": lambda: content_digest(_sample()),
        "hash _Sample": lambda: hash(_sample()),
        "digest _Leaf": lambda: content_digest(_Leaf(weight=0.75, tint=_Tint.WARM)),
        "digest FrozenMapping": lambda: content_digest(FrozenMapping({"a": 1.0, "b": (2.0, 3.0)})),
        "hash FrozenMapping": lambda: hash(FrozenMapping({"a": 1.0, "b": (2.0, 3.0)})),
        "digest Enum member": lambda: content_digest(_Tint.WARM),
        "digest container": lambda: content_digest({"w": np.arange(3.0), "t": (1.0, "x")}),
        "name of _Sample": lambda: name_digest(_sample()),
    }


# ═════════════════════════════════════════════════════════════════════════════
# S5.6 — the ROUTE gate: one encoder
# ═════════════════════════════════════════════════════════════════════════════


@pytest.mark.rests_on(f"{_HERE}::TestS54EncoderCanonicalForms")
def test_s5_6_one_encoder_is_the_only_route(monkeypatch: pytest.MonkeyPatch) -> None:
    """Rebind the encoder's one recursive entry, ``content._encode``, in every
    module that binds it, to a decoy that prefixes its output; every subject
    (the mixin's digest and hash, a part, a mapping, an enum member, a container,
    a derived name) must MOVE, and must have called the decoy (activation: a
    cached digest or a name computed before the swap would read unmoved for the
    wrong reason, so every subject is built fresh under the decoy). A route
    that hashes its own bytes stays unmoved: that is a second encoder (X4). The
    same law over every production type is the root census's
    ``test_s5_6_every_roster_type_routes_through_the_one_encoder``."""
    subjects = _synthetic_route_subjects()
    honest: dict[str, object] = {}
    unreadable: dict[str, str] = {}
    for name, read in subjects.items():
        try:
            honest[name] = read()
            if read() != honest[name]:
                unreadable[name] = "unstable: two honest builds read differently"
        except Exception as err:  # reported, never swallowed: the row reds on it below
            unreadable[name] = f"{type(err).__name__}: {err}"
    original = content_module._encode
    calls: Counter[str] = Counter()
    current = [""]

    def decoy(value: Any, path: str, *rest: Any) -> bytes:
        calls[current[0]] += 1
        return b"\x00decoy" + original(value, path, *rest)

    rebound = 0
    for module in list(sys.modules.values()):
        namespace = getattr(module, "__dict__", None)
        if not isinstance(namespace, dict):
            continue
        for key, bound in list(namespace.items()):
            if bound is original:
                monkeypatch.setattr(module, key, decoy)
                rebound += 1
    require(rebound >= 1, "the decoy was installed nowhere")
    print(f"S5.6: decoy rebound in {rebound} module binding(s); {len(subjects)} subjects")
    unmoved, idle = [], []
    for name, read in subjects.items():
        if name in unreadable:
            continue
        current[0] = name
        moved = read()
        if calls[name] == 0:
            idle.append(name)
        if moved == honest[name]:
            unmoved.append(name)
    require(
        not (unreadable or idle or unmoved),
        f"subjects with no honest reading: {unreadable}\n"
        f"subjects that never reached the encoder: {idle}\n"
        f"subjects unmoved by the decoy (a second encoder): {unmoved}",
    )


# ═════════════════════════════════════════════════════════════════════════════
# S5.7 — the RECORD fingerprint
# ═════════════════════════════════════════════════════════════════════════════

#: Pinned 2026-10-26 (layer DAG ruling 28), from literals and IEEE-exact arithmetic, so the pin is platform-free
#: (AP38). The type names inside these keys carry this module's name (a type outside ``orpheus`` is keyed as
#: ``module:qualname``), so moving this file moves the pins. The production fingerprint, the digests of mixture A
#: 2g, a hollow sphere and a material set, is the root census's
#: ``test_s5_7_the_production_digest_bytes_are_recorded``.
_PINNED: dict[str, str] = {
    "sample": "5cbb880052e098acf0fa271a04814863e161b81b6ae2a3f90966e3c64cd6102c",
    "leaf": "3d3b315a9a8f3c12813fdb79ea2df8c95dcdbf9d650eff82e47313ce90e7343b",
    "frozen mapping": "9a946144a1267c78e72c977c76701ac26b48350619622c7350e8219e1360c422",
    "enum member": "a7badf09c1c8ed91b5b5233017d32ae8a4eb8cbc43998131463dd4c18f00ece7",
}


def _fingerprint_values() -> dict[str, object]:
    return {
        "sample": _sample(),
        "leaf": _Leaf(weight=0.75, tint=_Tint.WARM),
        "frozen mapping": FrozenMapping({"a": 1.0, "b": (2.0, 3.0)}),
        "enum member": _Tint.COOL,
    }


@pytest.mark.rests_on(_HERE + "::test_s5_6_one_encoder_is_the_only_route")
def test_s5_7_the_digest_bytes_are_recorded() -> None:
    """RECORD (designed to red on any encoder edit): the producer fingerprint
    owed to every frozen-byte consumer of the digests (``vv-principles``,
    bit-identity). Re-pin only with the reason, in the commit message."""
    got = {name: content_digest(value).hex() for name, value in _fingerprint_values().items()}
    for name, pin in _PINNED.items():
        require(
            got[name] == pin,
            f"{name}: the digest bytes moved ({got[name]} != {pin}): every cache key of P3 "
            f"is invalidated; re-pin only with the reason",
        )
