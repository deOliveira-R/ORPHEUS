r"""Content identity is invariant under relocation, and the qualname that keys it is unique (layer DAG ruling 10).

A value's identity is its type and its parts; the module a class is defined in is neither. The schema tag of
:mod:`orpheus.numerics.content` therefore names a class by its ``__qualname__``, its ``__content_version__`` and the
ordered names of its parts, and an ``Enum`` member by its class's qualname and its own name, so moving a class
between modules (the layer-DAG carves move whole files) keeps every digest that encodes it. Two gates:

* RELOCATION. One source defines a frozen dataclass, a ``ContentIdentity`` subclass, an ``Enum`` and an object
  holding all three; it is imported as two modules from two files, and every value digests alike across them
  (``encode``, ``content_digest`` and ``encode_exact``). Negative controls, each on the relocated classes: a field
  value, a part name and a ``__content_version__`` change each move the digest, so the invariance is not a
  digest that ignores the class.
* UNIQUENESS. A qualname is an identity only if no two content-keyed classes share one. The population is every
  class defined in a module of the ``orpheus`` package (``pkgutil.walk_packages`` over all of it; ``tests/`` is
  NOT walked, since a test class never reaches a persistent key) that is content-keyed: a ``ContentIdentity``
  subclass, an ``Enum`` subclass, or a frozen dataclass. ``[M]`` 2026-10-22: 414 such classes (62
  ``ContentIdentity``, 18 ``Enum``, the rest frozen dataclasses), 0 shared qualnames. The positive control runs
  the same predicate over ALL ``orpheus`` classes, where 3 qualnames are shared (``DerivationError``,
  ``THParams``, ``THResult``), and over a synthetic pair built in this file.

A ``@dataclass(slots=True)`` returns a NEW class and leaves the one it was given alive among the base's
subclasses until it is collected, so a raw subclass walk sees two objects under one module and qualname; the
population resolves each class through its module (``module.qualname``) and counts the object found there. A
class whose qualname holds ``<locals>`` cannot be resolved and is counted as itself: two calls of one factory
are then two classes with one qualname, which is exactly the collision the gate refuses.

First red (the tree before ruling 10, ``refactor/layer-dag-round1`` ``e1bad58d``): every relocation leg reds,
the schema tag spelling ``module.qualname``.
"""

from __future__ import annotations

import dataclasses
import enum
import importlib
import importlib.util
import pkgutil
import sys
import textwrap
from collections import Counter
from collections.abc import Iterable
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

import orpheus
from orpheus.numerics.content import ContentIdentity, content_digest, encode, encode_exact
from tests.gates._content_identity_helpers import require

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/numerics/test_content_identity_relocation.py"


# ═════════════════════════════════════════════════════════════════════════════
# RELOCATION — one source, two modules
# ═════════════════════════════════════════════════════════════════════════════

_SOURCE = '''
from __future__ import annotations

import enum
from dataclasses import dataclass

from orpheus.numerics.content import ContentIdentity


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
'''


def _load(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, name: str, *, field: str = "x", version: int = 1) -> ModuleType:
    """The source imported from its own file as module ``name`` (registered, as an import registers it)."""
    path = tmp_path / f"{name}.py"
    path.write_text(textwrap.dedent(_SOURCE.format(field=field, version=version)))
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None  # type narrowing
    module = importlib.util.module_from_spec(spec)
    monkeypatch.setitem(sys.modules, name, module)
    spec.loader.exec_module(module)
    return module


def _values(m: ModuleType, x: float = 1.5, *, field: str = "x") -> dict[str, Any]:
    part = m.Part(**{field: x})
    return {
        "enum member": m.Shade.DARK,
        "frozen dataclass": part,
        "ContentIdentity": m.Whole(**{field: x}, part=part, shade=m.Shade.LIGHT),
    }


_ENCODINGS = {"encode": encode, "content_digest": content_digest, "encode_exact": encode_exact}


@pytest.mark.parametrize("encoding", sorted(_ENCODINGS))
def test_a_relocated_class_keeps_its_digest(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, encoding: str) -> None:
    """One source in two files, two modules: every value encodes alike in both."""
    a = _load(tmp_path, monkeypatch, "wp0_relocation_a")
    b = _load(tmp_path, monkeypatch, "wp0_relocation_b")
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

    def test_a_field_value_moves_it(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        a = _load(tmp_path, monkeypatch, "wp0_relocation_a")
        b = _load(tmp_path, monkeypatch, "wp0_relocation_b")
        for (kind, va), vb in zip(_values(a, 1.5).items(), _values(b, 2.5).values()):
            if kind != "enum member":
                require(content_digest(va) != content_digest(vb), f"{kind}: a field value did not move the digest")

    def test_a_part_name_moves_it(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        a = _load(tmp_path, monkeypatch, "wp0_relocation_a")
        b = _load(tmp_path, monkeypatch, "wp0_relocation_b", field="y")
        for (kind, va), vb in zip(_values(a).items(), _values(b, field="y").values()):
            if kind != "enum member":
                require(content_digest(va) != content_digest(vb), f"{kind}: a part name did not move the digest")

    def test_a_content_version_bump_moves_it(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        a = _load(tmp_path, monkeypatch, "wp0_relocation_a")
        b = _load(tmp_path, monkeypatch, "wp0_relocation_b", version=2)
        require(content_digest(_values(a)["ContentIdentity"]) != content_digest(_values(b)["ContentIdentity"]),
                "a __content_version__ bump did not move the digest")

    def test_a_renamed_enum_member_moves_it(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        a = _load(tmp_path, monkeypatch, "wp0_relocation_a")
        b = _load(tmp_path, monkeypatch, "wp0_relocation_b")
        require(content_digest(a.Shade.DARK) != content_digest(b.Shade.LIGHT), "two enum members digest alike")


# ═════════════════════════════════════════════════════════════════════════════
# UNIQUENESS — no two content-keyed classes share a qualname
# ═════════════════════════════════════════════════════════════════════════════


def _subclasses(root: type) -> set[type]:
    found: set[type] = set()
    stack = [root]
    while stack:
        for cls in type.__subclasses__(stack.pop()):
            if cls not in found:
                found.add(cls)
                stack.append(cls)
    return found


def _resolved(cls: type) -> type:
    """The class its module holds under its qualname (a ``slots=True`` dataclass's discarded original resolves to
    the class that replaced it); a class its module does not hold (``<locals>``) is itself."""
    holder: Any = sys.modules.get(cls.__module__)
    for part in cls.__qualname__.split("."):
        holder = getattr(holder, part, None)
    return holder if isinstance(holder, type) else cls


def _orpheus_classes() -> set[type]:
    """Every class defined in a module of the ``orpheus`` package, after importing every such module."""
    for info in pkgutil.walk_packages(orpheus.__path__, prefix="orpheus."):
        importlib.import_module(info.name)
    defined = {cls for cls in _subclasses(object)
               if isinstance(vars(cls).get("__module__"), str) and cls.__module__.split(".")[0] == "orpheus"}
    return {_resolved(cls) for cls in defined}


def is_content_keyed(cls: type) -> bool:
    """A class the content encoder names by its qualname: a ``ContentIdentity``, an ``Enum`` or a frozen dataclass."""
    params = getattr(cls, "__dataclass_params__", None)
    return issubclass(cls, (ContentIdentity, enum.Enum)) or (params is not None and params.frozen)


def shared_qualnames(classes: Iterable[type]) -> dict[str, list[str]]:
    """Each qualname held by two or more classes, with the modules defining them."""
    by_name: dict[str, list[str]] = {}
    for cls in classes:
        by_name.setdefault(cls.__qualname__, []).append(cls.__module__)
    return {name: sorted(modules) for name, modules in by_name.items() if len(modules) > 1}


def test_no_two_content_keyed_classes_share_a_qualname() -> None:
    """The schema tag keys on the qualname, so a shared one would let two types' values share a digest."""
    keyed = [cls for cls in _orpheus_classes() if is_content_keyed(cls)]
    counts = Counter("ContentIdentity" if issubclass(c, ContentIdentity) else "Enum" if issubclass(c, enum.Enum)
                     else "frozen dataclass" for c in keyed)
    require(counts["ContentIdentity"] >= 60 and counts["Enum"] >= 15 and counts["frozen dataclass"] >= 300,
            f"activation: the population is {dict(counts)}")
    shared = shared_qualnames(keyed)
    require(not shared, f"content-keyed classes sharing a qualname (rename one; the digest keys on it): {shared}")


@pytest.mark.rests_on(f"{_HERE}::test_no_two_content_keyed_classes_share_a_qualname")
class TestTheUniquenessGateSeesACollision:
    """Positive controls: the predicate that reads 0 above reads a collision where one exists."""

    def test_over_every_orpheus_class(self) -> None:
        """Not every class is content-keyed; among all of them three qualnames are shared."""
        shared = shared_qualnames(_orpheus_classes())
        known = {"DerivationError", "THParams", "THResult"}
        require(known <= shared.keys(), f"the known shared qualnames {sorted(known)} not all found: {shared}")

    def test_over_a_synthetic_pair(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        a = _load(tmp_path, monkeypatch, "wp0_relocation_a")
        b = _load(tmp_path, monkeypatch, "wp0_relocation_b")
        pair = [a.Whole, b.Whole, a.Part, a.Shade]
        require(all(is_content_keyed(c) for c in pair), "activation: the synthetic classes are content-keyed")
        require(shared_qualnames(pair) == {"Whole": ["wp0_relocation_a", "wp0_relocation_b"]},
                f"the synthetic pair is not reported: {shared_qualnames(pair)}")


# ═════════════════════════════════════════════════════════════════════════════
# THE ASYMMETRY — the exact encoder's foreign type tags keep the module
# ═════════════════════════════════════════════════════════════════════════════

_FOREIGN_SOURCE = '''
class Label(str):
    pass


class Count(int):
    pass
'''


@pytest.mark.parametrize("kind", ["str subclass", "int subclass"])
def test_a_foreign_like_subclass_keeps_its_module_in_the_exact_key(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, kind: str
) -> None:
    """A ``str`` or ``int`` subclass is a type no uniqueness gate enumerates (it may be foreign), so its tag in
    ``encode_exact`` keeps ``module.qualname``: two same-qualname subclasses in two modules key differently."""
    modules = []
    for name in ("wp0_foreign_a", "wp0_foreign_b"):
        path = tmp_path / f"{name}.py"
        path.write_text(_FOREIGN_SOURCE)
        spec = importlib.util.spec_from_file_location(name, path)
        assert spec is not None and spec.loader is not None  # type narrowing
        module = importlib.util.module_from_spec(spec)
        monkeypatch.setitem(sys.modules, name, module)
        spec.loader.exec_module(module)
        modules.append(module)
    a, b = (m.Label("u235") if kind == "str subclass" else m.Count(92) for m in modules)
    require(type(a).__qualname__ == type(b).__qualname__ and type(a).__module__ != type(b).__module__,
            "activation: one qualname, two modules")
    require(a == b, "activation: the two values are equal as builtins")
    require(encode_exact(a) != encode_exact(b),
            f"the exact key of a {kind} dropped its module: two modules' types share a key")
