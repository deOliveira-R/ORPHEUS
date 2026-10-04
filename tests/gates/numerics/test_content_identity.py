r"""The content encoder and the ``ContentIdentity`` mixin (#405 P1 step 5).

Spec: ``.claude/plans/reference_p1_spec.md`` §1.5. This module carries the
gates on the ENCODER and on the type POPULATION:

* S5.1 — digests and hashes are the same in every process (``PYTHONHASHSEED``
  1 and 2), with a control that the harness sees salting;
* S5.4 — the canonical forms: one double per real scalar (``True == 1 ==
  1.0``), ``-0.0`` is ``+0.0``, NaN and an integer beyond 2**53 refused with
  the path, containers keep the distinctions ``==`` keeps, a sparse matrix is
  the matrix, anything else is ``ContentlessError``;
* S5.6 — the ROUTE gate: one encoder. Rebinding it moves the digest of every
  type in the three rosters and every derived space name (the axes product,
  the three trace spaces, the full-field direct sum);
* S5.7 — the RECORD fingerprint (pinned after 5b);
* S5.8 — the schema tag (ruling of 2026-10-01): a field added, renamed or
  reordered, a version bump, a class moved or renamed, each changes the digest;
* S5.9 — conformance: every ``ContentIdentity`` class takes ``__eq__`` and
  ``__hash__`` from the mixin (or a declared override), is ``eq=False`` if a
  dataclass, and has a roster entry (the population S5.3 quantifies over).

The per-type legs of S5.2, S5.3 and S5.4 live in the owning trees
(``tests/gates/{data,geometry,mesh}/test_content_identity_*.py``).

First red, measured on ``main`` ``1dc31163``: collection fails with
``ModuleNotFoundError: No module named 'orpheus.numerics.content'``.
"""

from __future__ import annotations

import dataclasses
import enum
import inspect
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

import orpheus.numerics.content as content_module
from orpheus.numerics.content import ContentIdentity, ContentlessError, content_digest, encode
from tests.gates._content_identity_helpers import Entry, require
from tests.gates.data.test_content_identity_data import ROSTER as DATA_ROSTER
from tests.gates.geometry.test_content_identity_geometry import ROSTER as GEOMETRY_ROSTER
from tests.gates.mesh.test_content_identity_mesh import ROSTER as MESH_ROSTER
from tests.gates.numerics.test_content_identity_axis import ROSTER as AXIS_ROSTER
from tests.gates.numerics.test_content_identity_mesh_free import ROSTER as MESH_FREE_ROSTER
from tests.gates.numerics.test_content_identity_question import ROSTER as QUESTION_ROSTER
from tests.gates.numerics.test_content_identity_observable import ROSTER as OBSERVABLE_ROSTER
from tests.gates.numerics.test_enclosure import ROSTER as ENCLOSURE_ROSTER
from tests.gates.reference.test_readings import ROSTER as READING_ROSTER
from tests.gates.reference.test_published import ROSTER as PUBLISHED_ROSTER
from tests.gates.reference.test_reference_certificate import ROSTER as CERTIFICATE_ROSTER
from tests.gates.specification.test_content_identity_specification import ROSTER as SPECIFICATION_ROSTER

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/numerics/test_content_identity.py"
_ROOT = Path(__file__).resolve().parents[3]
ROSTER: tuple[Entry, ...] = DATA_ROSTER + GEOMETRY_ROSTER + MESH_ROSTER + AXIS_ROSTER + MESH_FREE_ROSTER + QUESTION_ROSTER + SPECIFICATION_ROSTER + ENCLOSURE_ROSTER + READING_ROSTER + OBSERVABLE_ROSTER + PUBLISHED_ROSTER + CERTIFICATE_ROSTER


# ═════════════════════════════════════════════════════════════════════════════
# S5.1 — the same digests and hashes in every process
# ═════════════════════════════════════════════════════════════════════════════

_SEED_SCRIPT = textwrap.dedent(
    """
    import numpy as np
    import orpheus
    from orpheus.data.materials import Materials
    from orpheus.derivations.common.xs_library import get_mixture
    from orpheus.geometry import BC, CoordSystem, StructuredGeometry
    from orpheus.geometry.boundary import AlbedoBoundary, IsotropicReturn
    from orpheus.mesh import Mesh2D
    from orpheus.numerics.axis import EnergyAxis
    from orpheus.numerics.content import content_digest
    from orpheus.numerics.space import FunctionSpace

    print("FILE", orpheus.__file__)
    values = {
        "mixture A 4g": get_mixture("A", "4g"),
        "materials": Materials({0: get_mixture("A", "2g"), 3: get_mixture("B", "2g")}),
        "bc albedo": BC("albedo", {"albedo": 0.3}),
        "hollow sphere": StructuredGeometry(
            coord=CoordSystem.SPHERICAL, breakpoints=(0.5, 1.0, 2.0), mat_ids=(0, 1),
            boundaries=(AlbedoBoundary(0.3, IsotropicReturn("x", -1)), BC.vacuum)),
        "mesh2d": Mesh2D(np.array([0.0, 1.0, 2.0]), np.array([0.0, 1.0]), np.array([[0], [1]]),
            face_laws={"xmin": BC.reflective, "xmax": BC.vacuum,
                       "ymin": BC.reflective, "ymax": BC.vacuum}),
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

        from orpheus.numerics.content import FrozenMapping

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

_MODULE_A, _MODULE_B = "tests.gates.numerics._schema_a", "tests.gates.numerics._schema_b"


def _schema_class(qualname: str, fields: tuple[str, ...], *, version: int = 1,
                  module: str = _MODULE_A, mixin: bool = True) -> type:
    """A frozen dataclass named ``module.qualname`` with float ``fields``."""
    namespace: dict[str, Any] = {
        # Real types, not strings: a string annotation sends ``dataclass`` to
        # ``sys.modules[module]``, which a synthetic module name is not in.
        "__annotations__": {name: float for name in fields},
        "__module__": module, "__qualname__": qualname,
    }
    if mixin:
        namespace["__content_version__"] = version
    cls = type(qualname, (ContentIdentity,) if mixin else (), namespace)
    return dataclass(frozen=True, eq=not mixin)(cls)


def _digest(cls: type, **values: float) -> bytes:
    return content_digest(cls(**values))


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
            pytest.param(lambda: _schema_class("Point", ("x",), module=_MODULE_B), {"x": 1.0}, id="the_class_moved"),
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

    def test_a_production_version_bump_moves_the_digest(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """On a production type, through a FRESH instance (a cached digest
        would mask the bump: the activation half of a route gate)."""
        from orpheus.data.macro_xs.mixture import Mixture
        from orpheus.derivations.common.xs_library import get_mixture

        before = content_digest(get_mixture("A", "2g"))
        monkeypatch.setattr(Mixture, "__content_version__", 2, raising=False)
        require(content_digest(get_mixture("A", "2g")) != before, "a version bump did not move the digest")


# ═════════════════════════════════════════════════════════════════════════════
# S5.9 — conformance: the mixin's dunders, and the population
# ═════════════════════════════════════════════════════════════════════════════

#: The classes the plan rules may override ``__eq__`` (the string arm, kept by
#: the user's ruling); they still take ``__hash__`` from the mixin.
_DECLARED_EQ_OVERRIDES = frozenset({"VacuumInflow", "ReflectiveBoundary"})

#: The packages whose classes S5.9 walks.
_PACKAGES = ("orpheus.data", "orpheus.geometry", "orpheus.mesh", "orpheus.numerics", "orpheus.specification", "orpheus.reference", "orpheus.transport")


def _content_classes() -> list[type]:
    import importlib
    import pkgutil

    for name in _PACKAGES:
        package = importlib.import_module(name)
        for info in pkgutil.walk_packages(package.__path__, prefix=f"{name}."):
            importlib.import_module(info.name)
    found: list[type] = []
    stack = list(ContentIdentity.__subclasses__())
    while stack:
        cls = stack.pop()
        if cls not in found:
            found.append(cls)
            stack.extend(cls.__subclasses__())
    return [c for c in found if c.__module__.startswith("orpheus.")]


def _owner(cls: type, name: str) -> type:
    return next(k for k in cls.__mro__ if name in vars(k))


@pytest.mark.rests_on(f"{_HERE}::TestS58Schema")
def test_s5_9_every_content_class_takes_the_mixins_dunders() -> None:
    """First red (the carve's): ``Mixture``, ``Mesh1D`` and ``CellEdges`` define
    ``__eq__``/``__hash__`` in their own class dicts today; ``FaceLaws``
    inherits ``Mapping.__eq__``; a frozen dataclass declared with ``eq=True``
    generates both and shadows the mixin."""
    classes = _content_classes()
    require(len(classes) >= 17, f"activation: {len(classes)} ContentIdentity classes found")
    problems: list[str] = []
    for cls in classes:
        # The FUNCTION, not its owner: a class that overrides ``__eq__`` must
        # re-bind ``__hash__ = ContentIdentity.__hash__`` (defining ``__eq__``
        # nulls an inherited hash), which is the mixin's own hash.
        if cls.__hash__ is not ContentIdentity.__hash__:
            problems.append(f"{cls.__qualname__}.__hash__ comes from {_owner(cls, '__hash__').__qualname__}")
        eq_owner = _owner(cls, "__eq__")
        if eq_owner is not ContentIdentity and not (
            cls.__qualname__ in _DECLARED_EQ_OVERRIDES and eq_owner is cls
        ):
            problems.append(f"{cls.__qualname__}.__eq__ comes from {eq_owner.__qualname__}")
        params = getattr(cls, "__dataclass_params__", None)
        if params is not None and (params.eq or not params.frozen):
            problems.append(f"{cls.__qualname__} is a dataclass with eq={params.eq}, frozen={params.frozen}")
    require(not problems, "\n".join(problems))


@pytest.mark.rests_on(f"{_HERE}::test_s5_9_every_content_class_takes_the_mixins_dunders")
def test_s5_9_the_population_is_the_roster() -> None:
    """Every concrete ``ContentIdentity`` class has a roster entry, so S5.3
    perturbs its every part; every roster entry marked content-identity is one.
    A new content type reds here until it joins a roster."""
    found = {c for c in _content_classes() if not inspect.isabstract(c)}
    rostered = {e.cls for e in ROSTER if e.content_identity}
    require(found - rostered == set(), f"content types with no roster entry: {sorted(c.__qualname__ for c in found - rostered)}")
    require(rostered - found == set(), f"roster entries without the mixin: {sorted(c.__qualname__ for c in rostered - found)}")
    print(f"S5.9 population: {len(found)} concrete ContentIdentity classes")


# ═════════════════════════════════════════════════════════════════════════════
# S5.6 — the ROUTE gate: one encoder
# ═════════════════════════════════════════════════════════════════════════════


def _space_names() -> dict[str, Callable[[], str]]:
    """Every derived space name the plan moves onto the encoder (5a)."""
    from orpheus.numerics.axis import Axis, BasisKind, EnergyAxis
    from orpheus.numerics.face_layout import FaceLayout
    from orpheus.numerics.quadrature import Quadrature
    from orpheus.numerics.space import FunctionSpace
    from orpheus.numerics.spaces.angular_trace_space import AngularTraceSpace
    from orpheus.numerics.spaces.full_field_space import FullFieldSpace
    from orpheus.numerics.spaces.radial_characteristic_space import (
        RadialCharacteristicBoundarySpace,
        RadialCharacteristicInteriorSpace,
    )
    from orpheus.numerics.spaces.scalar_trace_space import ScalarTraceSpace

    def angular() -> AngularTraceSpace:
        quad = Quadrature.gauss_legendre(4)
        layout = FaceLayout.from_named_shapes([("xmin", (quad.N, 1)), ("xmax", (quad.N, 1))])
        return AngularTraceSpace.from_quadrature_and_layout(quad, layout)

    def scalar() -> ScalarTraceSpace:
        return ScalarTraceSpace.for_faces([("xmin", ()), ("xmax", ())], 2, {"xmin": 1.0, "xmax": 1.0})

    levels, cell_volumes = (0, 2, 5), np.array([1.5, 2.5, 3.5, 4.5])
    return {
        "of_axes": lambda: FunctionSpace.of_axes(
            EnergyAxis.synthetic(2), Axis("harmonic", (3,), kind=BasisKind.MODAL)).name,
        "angular trace": lambda: angular().name,
        "scalar trace": lambda: scalar().name,
        "radial interior": lambda: RadialCharacteristicInteriorSpace.for_levels(
            levels, ng=2, nx=4, cell_volumes=cell_volumes).name,
        "radial boundary": lambda: RadialCharacteristicBoundarySpace.for_levels(
            levels, ng=2, nx=4, cell_volumes=cell_volumes).name,
        "full field over content-named blocks": lambda: FullFieldSpace.from_blocks(
            FunctionSpace.of_axes(EnergyAxis.synthetic(2)), scalar()).name,
        "full field over raw blocks": lambda: FullFieldSpace.from_blocks(
            FunctionSpace(name="bulk", shape=(4,)), FunctionSpace(name="trace", shape=(2,))).name,
    }


def _route_subjects() -> dict[str, Callable[[], object]]:
    subjects: dict[str, Callable[[], object]] = {}
    for entry in ROSTER:
        subjects[f"digest {entry.id}"] = lambda e=entry: content_digest(e.base())
        if entry.content_identity:
            subjects[f"hash {entry.id}"] = lambda e=entry: hash(e.base())
    subjects.update({f"name {k}": v for k, v in _space_names().items()})
    return subjects


@pytest.mark.rests_on(f"{_HERE}::test_s5_9_the_population_is_the_roster")
def test_s5_6_one_encoder_is_the_only_route(monkeypatch: pytest.MonkeyPatch) -> None:
    """Rebind the encoder's one recursive entry, ``content._encode``, in every
    module that binds it, to a decoy that prefixes its output; every subject
    (each roster type's digest and hash, each derived space name) must MOVE,
    and must have called the decoy (activation: a cached digest or a name
    computed before the swap would read unmoved for the wrong reason, so every
    subject is built fresh under the decoy). A type or a space still hashing
    its own bytes stays unmoved: that is a second encoder (X4)."""
    subjects = _route_subjects()
    honest: dict[str, object] = {}
    unreadable: dict[str, str] = {}
    for name, read in subjects.items():
        try:
            honest[name] = read()
            # Stability control: a reading that differs between two honest
            # builds (an identity hash) would read "moved" for the wrong reason.
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

#: Pinned 2026-10-02 at #405 P1 step 5 (macOS, Python 3.14). The three values are
#: built from literals and IEEE-exact arithmetic (mixture A is ``xs_library``
#: literals; no libm or LAPACK output enters), so the pin is platform-free (AP38).
_PIN_AFTER_5B = "PIN-AFTER-5B"
_PINNED: dict[str, str] = {
    "mixture A 2g": "016a3c66cc2169eaa3e6dd8a090909de7cabb95bcbc04bee5e03bd8f9779f48f",
    "hollow sphere, three intervals": "4c4b7298d06faa05c80dcdce74612a113a24694d798e3a6401f050936fb6ee6a",
    "materials of two mixtures": "343e62741213b34e1b615fb0f5f110adbc83dd5bf85936594c3124e95281d946",
}


def _fingerprint_values() -> dict[str, object]:
    from orpheus.data.materials import Materials
    from orpheus.derivations.common.xs_library import get_mixture
    from orpheus.geometry import BC, CoordSystem, StructuredGeometry
    from orpheus.geometry.boundary import AlbedoBoundary, SpecularReturn

    return {
        "mixture A 2g": get_mixture("A", "2g"),
        "hollow sphere, three intervals": StructuredGeometry(
            coord=CoordSystem.SPHERICAL, breakpoints=(0.25, 0.5, 1.0, 2.0), mat_ids=(0, 1, 0),
            boundaries=(AlbedoBoundary(0.5, SpecularReturn("x")), BC.vacuum)),
        "materials of two mixtures": Materials({0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")}),
    }


@pytest.mark.rests_on(_HERE + "::test_s5_6_one_encoder_is_the_only_route")
def test_s5_7_the_digest_bytes_are_recorded() -> None:
    """RECORD (designed to red on any encoder edit): the producer fingerprint
    owed to every frozen-byte consumer of the digests (``vv-principles``,
    bit-identity). Re-pin only with the reason, in the commit message."""
    got = {name: content_digest(value).hex() for name, value in _fingerprint_values().items()}
    if any(pin == _PIN_AFTER_5B for pin in _PINNED.values()):
        pytest.fail(f"pin after 5b: replace the placeholders of _PINNED with {got}")
    for name, pin in _PINNED.items():
        require(
            got[name] == pin,
            f"{name}: the digest bytes moved ({got[name]} != {pin}): every cache key of P3 "
            f"is invalidated; re-pin only with the reason",
        )
