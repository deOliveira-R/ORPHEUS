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
import pkgutil
import sys
from collections import Counter
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import numpy as np
import pytest
import scipy.sparse as sp

import orpheus
from orpheus.values.content import ContentIdentity
from tests.gates._content_identity_helpers import require
from tests.gates.values.test_content_identity_relocation import (
    _ORPHEUS_A,
    _ORPHEUS_B,
    _load,
    _no_synthetic_class_outlives_its_test,  # noqa: F401  (autouse: the synthetic pair's classes die with its test)
    _subclasses,
)

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/test_content_identity_census.py"


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
