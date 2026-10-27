r"""Censuses of the content-identity contract over every ``orpheus`` package (layer DAG ruling 28, 2026-10-26).

``orpheus.values.content`` defines what a value's key is; these rows check that every class of the tree that
the encoder names keeps the contract, so they range over thirteen packages above ``values`` and live at the root
of ``tests/gates/``, beside ``test_layer_imports.py`` (the root is not judged by the test-layer gate):

* UNIQUENESS. A key names an ``orpheus`` type by its qualname, so no two nameable ``orpheus`` classes may share
  one. Moved from ``tests/gates/values/test_content_identity_relocation.py``, which keeps the relocation and
  foreign-type rules and the synthetic two-module loader these controls reuse. The population is every class
  defined in a module ``pkgutil.walk_packages`` finds under ``orpheus`` (all imported; ``tests/`` is not walked)
  that some encoder names: a ``ContentIdentity`` subclass, an ``Enum`` subclass, a frozen dataclass, or a
  subclass of ``str``, ``bytes``, ``int``, ``float``, ``np.ndarray``, ``tuple``, ``list``, ``dict``, ``set``,
  ``frozenset`` or a scipy sparse base. ``[M]`` 2026-10-22: 447 such classes of 697, 0 sharing a qualname. The
  positive controls: the predicate over ALL ``orpheus`` classes finds the 3 shared qualnames (``DerivationError``,
  ``THParams``, ``THResult``), and it reports a synthetic pair.
"""

from __future__ import annotations

import enum
import importlib
import inspect
import pkgutil
import sys
from collections import Counter
from collections.abc import Callable, Iterable
from pathlib import Path
from typing import Any

import numpy as np
import pytest
import scipy.sparse as sp

import orpheus
import orpheus.data
import orpheus.derivations.continuous.characteristic
import orpheus.diffusion
import orpheus.geometry
import orpheus.mesh
import orpheus.numerics
import orpheus.reference
import orpheus.specification
import orpheus.transport
import orpheus.values
import orpheus.values.content as content_module
from orpheus.diffusion import ZeroFluxCondition  # the one content type of orpheus.diffusion; it has no fields
from orpheus.numerics.axis import Declared, Modal
from orpheus.values.content import ContentIdentity, content_digest
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
from tests.gates.transport.test_content_identity_transport import ROSTER as TRANSPORT_ROSTER
from tests.gates.derivations.test_characteristic_reference import ROSTER as CHARACTERISTIC_ROSTER
from tests.gates.values.test_content_identity_relocation import (
    _ORPHEUS_A,
    _ORPHEUS_B,
    _load,
    _no_synthetic_class_outlives_its_test,  # noqa: F401  (autouse: the synthetic pair's classes die with its test)
    _subclasses,
)

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/test_content_identity_census.py"
_VALUES = "tests/gates/values/test_content_identity.py"

DIFFUSION_ROSTER: tuple[Entry, ...] = (
    Entry(
        cls=ZeroFluxCondition, base=ZeroFluxCondition, parts=(), perturb={},
        pairs=(("two builds", ZeroFluxCondition, ZeroFluxCondition),),
    ),
)

ROSTER: tuple[Entry, ...] = DATA_ROSTER + GEOMETRY_ROSTER + MESH_ROSTER + AXIS_ROSTER + MESH_FREE_ROSTER + QUESTION_ROSTER + SPECIFICATION_ROSTER + TRANSPORT_ROSTER + ENCLOSURE_ROSTER + READING_ROSTER + OBSERVABLE_ROSTER + PUBLISHED_ROSTER + CERTIFICATE_ROSTER + CHARACTERISTIC_ROSTER + DIFFUSION_ROSTER


# ═════════════════════════════════════════════════════════════════════════════
# UNIQUENESS — no two nameable orpheus classes share a qualname
# ═════════════════════════════════════════════════════════════════════════════


def _resolved(cls: type) -> type:
    """The class its module holds under its qualname (a ``slots=True`` dataclass's discarded original resolves to
    the class that replaced it); a class its module does not hold (``<locals>``) is itself."""
    holder: Any = sys.modules.get(cls.__module__)
    for part in cls.__qualname__.split("."):
        holder = getattr(holder, part, None)
    return holder if isinstance(holder, type) else cls


def _orpheus_classes() -> set[type]:
    """Every class defined in a module ``pkgutil.walk_packages`` finds under ``orpheus``, after importing each
    (a module this file registers under ``orpheus.`` is not walked, so it is not counted)."""
    walked = {orpheus.__name__}
    for info in pkgutil.walk_packages(orpheus.__path__, prefix="orpheus."):
        importlib.import_module(info.name)
        walked.add(info.name)
    defined = {cls for cls in _subclasses(object) if vars(cls).get("__module__") in walked}
    return {_resolved(cls) for cls in defined}


_NAMED_BASES: tuple[type, ...] = (str, bytes, int, float, np.ndarray, tuple, list, dict, set, frozenset,
                                  sp.spmatrix, sp.sparray)


def _kind(cls: type) -> str | None:
    """The kind under which some encoder names ``cls``, or ``None`` when no encoder does."""
    if issubclass(cls, ContentIdentity):
        return "ContentIdentity"
    if issubclass(cls, enum.Enum):
        return "Enum"
    params = getattr(cls, "__dataclass_params__", None)
    if params is not None and params.frozen:
        return "frozen dataclass"
    if issubclass(cls, _NAMED_BASES):
        return "builtin, array or sparse subclass"
    return None


def shared_qualnames(classes: Iterable[type]) -> dict[str, list[str]]:
    """Each qualname held by two or more classes, with the modules defining them."""
    by_name: dict[str, list[str]] = {}
    for cls in classes:
        by_name.setdefault(cls.__qualname__, []).append(cls.__module__)
    return {name: sorted(modules) for name, modules in by_name.items() if len(modules) > 1}


def test_no_two_nameable_orpheus_classes_share_a_qualname() -> None:
    """A key names an orpheus type by its qualname, so a shared one would let two types' values share a key."""
    named = [cls for cls in _orpheus_classes() if _kind(cls) is not None]
    counts = Counter(_kind(cls) for cls in named)
    require(counts["ContentIdentity"] >= 60 and counts["Enum"] >= 15 and counts["frozen dataclass"] >= 300
            and counts["builtin, array or sparse subclass"] >= 1, f"activation: the population is {dict(counts)}")
    shared = shared_qualnames(named)
    require(not shared, f"nameable orpheus classes sharing a qualname (rename one; the key names it): {shared}")


@pytest.mark.rests_on(f"{_HERE}::test_no_two_nameable_orpheus_classes_share_a_qualname")
class TestTheUniquenessGateSeesACollision:
    """Positive controls: the predicate that reads 0 above reads a collision where one exists."""

    def test_over_every_orpheus_class(self) -> None:
        """Not every class is nameable; among all of them three qualnames are shared."""
        shared = shared_qualnames(_orpheus_classes())
        known = {"DerivationError", "THParams", "THResult"}
        require(known <= shared.keys(), f"the known shared qualnames {sorted(known)} not all found: {shared}")

    def test_over_a_synthetic_pair(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        a = _load(tmp_path, monkeypatch, _ORPHEUS_A)
        b = _load(tmp_path, monkeypatch, _ORPHEUS_B)
        pair = [a.Whole, b.Whole, a.Part, a.Shade, a.Stencil]
        require(all(_kind(c) is not None for c in pair), "activation: the synthetic classes are nameable")
        require(shared_qualnames(pair) == {"Whole": [_ORPHEUS_A, _ORPHEUS_B]},
                f"the synthetic pair is not reported: {shared_qualnames(pair)}")


# ═════════════════════════════════════════════════════════════════════════════
# S5.9 — conformance: the mixin's dunders, and the population
# ═════════════════════════════════════════════════════════════════════════════

# No class overrides ``__eq__``: the two that did (the string arm of the retired ``VacuumInflow`` and ``ReflectiveBoundary``) retired in WF2.4,
# so the walk below admits the mixin's ``__eq__`` alone.

#: The packages whose classes S5.9 walks. Of ``orpheus.derivations``, the package whose derivation is a reference
#: solution (its receiver is a traced-memo key); the rest of ``orpheus.derivations`` is NOT walked. The
#: trajectory-resolvent package was the second until step (e2) of the characteristic-reference campaign deleted it.
_PACKAGES = (orpheus.values, orpheus.data, orpheus.geometry, orpheus.mesh, orpheus.numerics, orpheus.specification, orpheus.reference,
             orpheus.transport, orpheus.derivations.continuous.characteristic, orpheus.diffusion)


def _content_classes() -> list[type]:
    import importlib
    import pkgutil

    for package in _PACKAGES:
        for info in pkgutil.walk_packages(package.__path__, prefix=f"{package.__name__}."):
            importlib.import_module(info.name)
    found: list[type] = []
    stack = list(ContentIdentity.__subclasses__())
    while stack:
        cls = stack.pop()
        if cls not in found:
            found.append(cls)
            stack.extend(cls.__subclasses__())
    # the population is STATED: the classes defined under the walked packages, never whatever else the run has imported
    return [c for c in found if c.__module__.startswith(tuple(f"{package.__name__}." for package in _PACKAGES))]


def _owner(cls: type, name: str) -> type:
    return next(k for k in cls.__mro__ if name in vars(k))


@pytest.mark.rests_on(f"{_VALUES}::TestS58Schema")
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
        if eq_owner is not ContentIdentity:
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
# S5.6 over the tree — every roster type and every derived space name routes through the one encoder
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
        layout = FaceLayout.from_named_shapes([("xmin", (quad.n_ordinates, 1)), ("xmax", (quad.n_ordinates, 1))])
        return AngularTraceSpace.from_quadrature_and_layout(quad, layout)

    def scalar() -> ScalarTraceSpace:
        return ScalarTraceSpace.for_faces([("xmin", ()), ("xmax", ())], 2, {"xmin": 1.0, "xmax": 1.0})

    levels, cell_volumes = (0, 2, 5), np.array([1.5, 2.5, 3.5, 4.5])
    return {
        "of_axes": lambda: FunctionSpace.of_axes(
            EnergyAxis.synthetic(2), Axis("harmonic", (3,), origin=Declared(Modal(np.eye(3)[0])))).name,
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


@pytest.mark.rests_on(f"{_HERE}::test_s5_9_the_population_is_the_roster", f"{_VALUES}::test_s5_6_one_encoder_is_the_only_route")
def test_s5_6_every_roster_type_routes_through_the_one_encoder(monkeypatch: pytest.MonkeyPatch) -> None:
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
# S5.7 over the tree — the production RECORD fingerprint
# ═════════════════════════════════════════════════════════════════════════════

#: Pinned 2026-10-02 at #405 P1 step 5 (macOS, Python 3.14). The three values are
#: built from literals and IEEE-exact arithmetic (mixture A is ``macro_xs.synthetic``
#: literals; no libm or LAPACK output enters), so the pin is platform-free (AP38).
_PIN_AFTER_5B = "PIN-AFTER-5B"
# Re-pinned 2026-10-22 (layer DAG ruling 10): the schema tag and the Enum encoding dropped the module of an orpheus
# type, so every pin moved; the old _schema_tag and the old Enum encoding swapped back in-process together reproduce
# each old pin byte for byte (the tag alone reproduces 11 of the 14; the commit body lists old -> new).
_PINNED: dict[str, str] = {
    "mixture A 2g": "905ce741c85384eefb289905ea4dbd8591499002889a16b1e481e73992d845e6",
    "hollow sphere, three intervals": "42c18f89853f810317774beeb060bc4ddc7186bfcaa8c843c062e3a39fcdb2d9",
    "materials of two mixtures": "3ae54bab38bf8af73ec713740aa69963d2b86aa73ae8b489b6e63b0eae6f34e6",
}


def _fingerprint_values() -> dict[str, object]:
    from orpheus.data.materials import Materials
    from orpheus.data.macro_xs.synthetic import get_mixture
    from orpheus.geometry import CoordSystem, StructuredGeometry

    return {
        "mixture A 2g": get_mixture("A", "2g"),
        "hollow sphere, three intervals": StructuredGeometry(
            coord=CoordSystem.SPHERICAL, breakpoints=(0.25, 0.5, 1.0, 2.0), mat_ids=(0, 1, 0),
            face_roles=(1, 0)),
        "materials of two mixtures": Materials({0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")}),
    }


@pytest.mark.rests_on(_HERE + "::test_s5_6_every_roster_type_routes_through_the_one_encoder", f"{_VALUES}::test_s5_7_the_digest_bytes_are_recorded")
def test_s5_7_the_production_digest_bytes_are_recorded() -> None:
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
