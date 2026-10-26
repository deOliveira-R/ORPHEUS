r"""A type's name inside a key: an ``orpheus`` type by its qualname, any other by ``module:qualname`` (layer DAG ruling 10).

A value's identity is its type and its parts; the module a class is defined in is neither. So
:func:`orpheus.values.content._type_name`, the one spelling of a type's name inside a key, names a type of
``orpheus`` by its ``__qualname__`` and any other type by ``module:qualname``; it is used at all five sites a key
names a type (the schema tag, an ``Enum`` member, and the exact encoder's container, subclass and sparse tags).
Three gates; the third, a census over every ``orpheus`` package, lives at the root since 2026-10-26
(``tests/gates/test_content_identity_census.py``, layer DAG ruling 28):

* RELOCATION. One source defines a frozen dataclass, a ``ContentIdentity`` subclass, an ``Enum``, a
  ``NamedTuple``, a ``str`` and an ``int`` subclass and a scipy sparse subclass; it is imported as two modules
  under ``orpheus.`` from two files, and every value encodes alike across them (``encode``, ``content_digest`` and
  ``encode_exact``), so a carve that moves a file between packages re-keys nothing. Negative controls, on the
  relocated classes: a field value, a part name, a ``__content_version__`` change and a different enum member
  each move the digest, so the invariance is not a digest blind to the class.
* FOREIGN TYPES KEEP THEIR MODULE. The same source imported as two modules OUTSIDE ``orpheus`` keys apart, one row
  per site, because no gate can enumerate foreign types: without the module a foreign ``Enum`` or frozen
  dataclass sharing a qualname with an ``orpheus`` one would share its key.
* UNIQUENESS. A qualname is an identity only if no two nameable ``orpheus`` classes share one. The population is
  every class defined in a module that ``pkgutil.walk_packages`` finds under ``orpheus`` (all of them are
  imported; ``tests/`` is NOT walked, since its types are foreign to the rule) that some encoder names: a
  ``ContentIdentity`` subclass, an ``Enum`` subclass, a frozen dataclass, or a subclass of ``str``, ``bytes``,
  ``int``, ``float``, ``np.ndarray``, ``tuple``, ``list``, ``dict``, ``set``, ``frozenset`` (a ``NamedTuple``
  included) or of a scipy sparse base (334 frozen dataclasses, 62 ``ContentIdentity``, 18 ``Enum``, 33 builtin
  subclasses, mostly ``NamedTuple``). ``[M]`` 2026-10-22: 447 such classes of 697 ``orpheus`` classes (391
  modules), 0 sharing a qualname. The positive controls: the same predicate over ALL ``orpheus`` classes finds
  the 3 shared qualnames (``DerivationError``, ``THParams``, ``THResult``), and it reports a synthetic pair.

A ``@dataclass(slots=True)`` returns a NEW class and leaves the one it was given alive among the base's
subclasses until it is collected, so a raw subclass walk sees two objects under one module and qualname; the
population resolves each class through its module (``module.qualname``) and counts the object found there. A
class whose qualname holds ``<locals>`` cannot be resolved and is counted as itself: two calls of one factory
are then two classes with one qualname, which is exactly the collision the gate refuses.

First red: the tree before ruling 10 (``e1bad58d``, the module in every tag) reds the relocation legs; a mutation
of each of the five sites reds a row here (the module put back at the schema tag, the ``Enum`` branch, the
container, subclass and sparse tags; the module taken away at the schema tag, the ``Enum`` branch, the subclass
and the sparse tag). Taking the module away at the container tag alone reds nothing: a container whose type is
not its base also carries the subclass tag, which names the type again.
"""

from __future__ import annotations

import gc
import importlib
import importlib.util
import sys
import textwrap
import warnings
from collections.abc import Iterator
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pytest

from orpheus.values.content import content_digest, encode, encode_exact
from tests.gates._content_identity_helpers import require

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/values/test_content_identity_relocation.py"

#: Two module names under ``orpheus`` (the rule names their types by qualname) and two outside it.
_ORPHEUS_A, _ORPHEUS_B = "orpheus._wp0_relocation_a", "orpheus._wp0_relocation_b"
_FOREIGN_A, _FOREIGN_B = "wp0_foreign_a", "wp0_foreign_b"


# ═════════════════════════════════════════════════════════════════════════════
# The source, imported as two modules
# ═════════════════════════════════════════════════════════════════════════════

_SOURCE = '''
from __future__ import annotations

import enum
from dataclasses import dataclass
from typing import NamedTuple

import scipy.sparse as sp

from orpheus.values.content import ContentIdentity


class Shade(enum.Enum):
    DARK = 1
    LIGHT = 2


@dataclass(frozen=True)
class Part:
    {field}: float


@dataclass(frozen=True, eq=False)
class Whole(ContentIdentity):
    __content_version__ = {version}
    {field}: float
    part: Part
    shade: Shade


class Pair(NamedTuple):
    left: int
    right: int


class Label(str):
    pass


class Count(int):
    pass


class Stencil(sp.csr_matrix):
    pass
'''


@pytest.fixture(autouse=True)
def _no_synthetic_class_outlives_its_test() -> Iterator[None]:
    """The synthetic ``orpheus._wp0_*`` classes are collected once their modules leave ``sys.modules``, so no
    population walk of a later test (the root census ``tests/gates/test_content_identity_census.py``: S5.9's
    roster and the qualname uniqueness gate) meets them as ``orpheus`` classes. Autouse fixtures are set up before ``monkeypatch``, so this teardown runs after it restores
    ``sys.modules``."""
    yield
    gc.collect()
    left = sorted(f"{c.__module__}.{c.__qualname__}" for c in _subclasses(object)
                  if str(vars(c).get("__module__", "")).startswith("orpheus._wp0_"))
    require(not left, f"synthetic orpheus classes outlived their test: {left}")


def _load(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, name: str, source: str = _SOURCE, *,
          field: str = "x", version: int = 1) -> ModuleType:
    """``source`` imported from its own file as module ``name`` (registered, as an import registers it)."""
    path = tmp_path / f"{name.rsplit('.', 1)[-1]}.py"
    path.write_text(textwrap.dedent(source.format(field=field, version=version)))
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None  # type narrowing
    module = importlib.util.module_from_spec(spec)
    monkeypatch.setitem(sys.modules, name, module)
    spec.loader.exec_module(module)
    return module


def _values(m: ModuleType, x: float = 1.5, *, field: str = "x") -> dict[str, Any]:
    """One value per site a key names a type at (the site in the key)."""
    part = m.Part(**{field: x})
    return {
        "Enum member": m.Shade.DARK,
        "schema tag (frozen dataclass)": part,
        "schema tag (ContentIdentity)": m.Whole(**{field: x}, part=part, shade=m.Shade.LIGHT),
        "container tag (NamedTuple)": m.Pair(1, 2),
        "subclass tag (str)": m.Label("u235"),
        "subclass tag (int)": m.Count(92),
        "sparse tag": m.Stencil(np.array([[1.0, 0.0], [0.0, 2.0]])),
    }


_ENCODINGS = {"encode": encode, "content_digest": content_digest, "encode_exact": encode_exact}


# ═════════════════════════════════════════════════════════════════════════════
# RELOCATION — an orpheus class moved between modules keeps every key
# ═════════════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("encoding", sorted(_ENCODINGS))
def test_a_relocated_class_keeps_its_digest(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, encoding: str) -> None:
    """One source in two files, two ``orpheus`` modules: every value encodes alike in both."""
    a = _load(tmp_path, monkeypatch, _ORPHEUS_A)
    b = _load(tmp_path, monkeypatch, _ORPHEUS_B)
    require(a.Whole is not b.Whole and a.Whole.__module__ != b.Whole.__module__, "activation: two modules, two classes")
    moved = []
    for (kind, va), vb in zip(_values(a).items(), _values(b).values()):
        require(type(va).__qualname__ == type(vb).__qualname__, f"activation: {kind} shares one qualname")
        if _ENCODINGS[encoding](va) != _ENCODINGS[encoding](vb):
            moved.append(kind)
    require(not moved, f"moving the class moved the {encoding} of: {moved}")


@pytest.mark.rests_on(f"{_HERE}::test_a_relocated_class_keeps_its_digest")
class TestTheRelocatedDigestStillReadsTheClass:
    """Negative controls: the invariance is not a digest blind to the class. Each leg holds the module fixed
    against the other leg's module (the two are relocated), so only the named change separates them."""

    _OBJECTS = ("schema tag (frozen dataclass)", "schema tag (ContentIdentity)")

    def test_a_field_value_moves_it(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        a = _load(tmp_path, monkeypatch, _ORPHEUS_A)
        b = _load(tmp_path, monkeypatch, _ORPHEUS_B)
        va, vb = _values(a, 1.5), _values(b, 2.5)
        for kind in self._OBJECTS:
            require(content_digest(va[kind]) != content_digest(vb[kind]), f"{kind}: a field value did not move the digest")

    def test_a_part_name_moves_it(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        a = _load(tmp_path, monkeypatch, _ORPHEUS_A)
        b = _load(tmp_path, monkeypatch, _ORPHEUS_B, field="y")
        va, vb = _values(a), _values(b, field="y")
        for kind in self._OBJECTS:
            require(content_digest(va[kind]) != content_digest(vb[kind]), f"{kind}: a part name did not move the digest")

    def test_a_content_version_bump_moves_it(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        a = _load(tmp_path, monkeypatch, _ORPHEUS_A)
        b = _load(tmp_path, monkeypatch, _ORPHEUS_B, version=2)
        kind = "schema tag (ContentIdentity)"
        require(content_digest(_values(a)[kind]) != content_digest(_values(b)[kind]),
                "a __content_version__ bump did not move the digest")

    def test_a_renamed_enum_member_moves_it(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        a = _load(tmp_path, monkeypatch, _ORPHEUS_A)
        b = _load(tmp_path, monkeypatch, _ORPHEUS_B)
        require(content_digest(a.Shade.DARK) != content_digest(b.Shade.LIGHT), "two enum members digest alike")


# ═════════════════════════════════════════════════════════════════════════════
# FOREIGN TYPES — the same qualname in two modules outside orpheus keys apart
# ═════════════════════════════════════════════════════════════════════════════

#: Per site, the encoder whose key names the type there (the by-equality encoder names no container, subclass
#: or sparse type: a NamedTuple equals its tuple under ``==``).
_FOREIGN_ROWS = {
    "Enum member": content_digest,
    "schema tag (frozen dataclass)": content_digest,
    "schema tag (ContentIdentity)": content_digest,
    "container tag (NamedTuple)": encode_exact,
    "subclass tag (str)": encode_exact,
    "subclass tag (int)": encode_exact,
    "sparse tag": encode_exact,
}


@pytest.mark.parametrize("kind", list(_FOREIGN_ROWS))
def test_a_foreign_type_keeps_its_module_in_the_key(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, kind: str) -> None:
    """Two same-qualname foreign types in two modules key apart: no gate can enumerate foreign types."""
    a = _load(tmp_path, monkeypatch, _FOREIGN_A)
    b = _load(tmp_path, monkeypatch, _FOREIGN_B)
    va, vb = _values(a)[kind], _values(b)[kind]
    require(type(va).__qualname__ == type(vb).__qualname__ and type(va).__module__ != type(vb).__module__,
            "activation: one qualname, two foreign modules")
    key = _FOREIGN_ROWS[kind]
    require(key(va) != key(vb), f"{kind}: two foreign modules' same-qualname types share a {key.__name__} key")


# ═════════════════════════════════════════════════════════════════════════════
# The subclass walk (the fixture above; the uniqueness census at the root imports it)
# ═════════════════════════════════════════════════════════════════════════════


def _subclasses(root: type) -> set[type]:
    """Every class below ``root``. Hashing ``typing``'s deprecated ``_UnionGenericAlias`` warns; the walk only
    visits it, so the warning is silenced here."""
    found: set[type] = set()
    stack = [root]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        while stack:
            for cls in type.__subclasses__(stack.pop()):
                if cls not in found:
                    found.add(cls)
                    stack.append(cls)
    return found
