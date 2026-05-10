"""Foundation tests for :mod:`orpheus.numerics.registry`.

Tests the RegistryMixin self-registration mechanism shipped in Phase
B / Issue 9.6 — ``__init_subclass__`` registration under the ``key=``
class-creation kwarg, ``create()`` factory, duplicate-key detection,
and registry-root isolation across independent hierarchies.

The pre-existing :mod:`tests.numerics.test_registry` file covers a
different concept (the quadrature-rule selection registry); this file
covers the new :class:`RegistryMixin` primitive itself.
"""

from __future__ import annotations

from typing import ClassVar

import pytest

from orpheus.numerics.registry import RegistryMixin


# Define each test's hierarchy in its own module-level scope so the
# ``__init_subclass__`` side effects are visible to the test.


class _RootA(RegistryMixin):
    registry: ClassVar[dict[str, type]] = {}

    @classmethod
    def _registry_base(cls) -> type:
        return _RootA


class _ConcreteA1(_RootA, key="a1"):
    pass


class _ConcreteA2(_RootA, key="a2"):
    pass


class _RootB(RegistryMixin):
    registry: ClassVar[dict[str, type]] = {}

    @classmethod
    def _registry_base(cls) -> type:
        return _RootB


class _ConcreteB1(_RootB, key="b1"):
    pass


# ─────────────────────────────────────────────────────────────────────
# Self-registration on class creation
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_subclass_registers_under_key():
    assert _RootA.registry["a1"] is _ConcreteA1
    assert _RootA.registry["a2"] is _ConcreteA2


@pytest.mark.foundation
def test_class_attribute_key_set():
    """``cls.key`` is set by ``__init_subclass__`` from the kwarg."""
    assert _ConcreteA1.key == "a1"
    assert _ConcreteA2.key == "a2"


@pytest.mark.foundation
def test_abstract_intermediate_not_registered():
    """A subclass without ``key=`` is not registered (allows abstract layers)."""

    class _Intermediate(_RootA):
        pass

    assert _Intermediate not in _RootA.registry.values()
    # And key remains None on the intermediate.
    assert _Intermediate.key is None


# ─────────────────────────────────────────────────────────────────────
# Factory
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_create_returns_concrete_instance():
    instance = _RootA.create("a1")
    assert isinstance(instance, _ConcreteA1)


@pytest.mark.foundation
def test_create_unknown_key_raises():
    with pytest.raises(KeyError, match="unknown"):
        _RootA.create("does_not_exist")


@pytest.mark.foundation
def test_create_unknown_key_lists_available():
    """Error message lists available keys to aid debugging."""
    with pytest.raises(KeyError) as exc:
        _RootA.create("typo")
    msg = str(exc.value)
    assert "a1" in msg
    assert "a2" in msg


@pytest.mark.foundation
def test_create_forwards_kwargs():
    """Kwargs passed to ``create`` reach the concrete constructor."""

    class _Root(RegistryMixin):
        registry: ClassVar[dict[str, type]] = {}

        @classmethod
        def _registry_base(cls) -> type:
            return _Root

    class _Concrete(_Root, key="cfg"):
        def __init__(self, x: int, y: int = 0) -> None:
            self.x = x
            self.y = y

    instance = _Root.create("cfg", x=1, y=2)
    assert instance.x == 1
    assert instance.y == 2


# ─────────────────────────────────────────────────────────────────────
# Duplicate-key detection
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_duplicate_key_raises():
    """Registering twice under the same key is a hard error."""

    class _Root(RegistryMixin):
        registry: ClassVar[dict[str, type]] = {}

        @classmethod
        def _registry_base(cls) -> type:
            return _Root

    class _C1(_Root, key="dup"):
        pass

    with pytest.raises(KeyError, match="duplicate"):

        class _C2(_Root, key="dup"):  # noqa: F841
            pass


# ─────────────────────────────────────────────────────────────────────
# Registry-root isolation
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_registries_are_isolated():
    """Distinct registry roots have independent dicts."""
    assert "a1" in _RootA.registry
    assert "a1" not in _RootB.registry
    assert "b1" in _RootB.registry
    assert "b1" not in _RootA.registry


@pytest.mark.foundation
def test_registries_disjoint_across_codebase():
    """BoundaryOperator and CellUpdateBase registries don't share keys
    just because they both use RegistryMixin."""
    from orpheus.geometry.boundary import BoundaryOperator
    from orpheus.sn.spatial.cell_update import CellUpdateBase

    # The two registries are distinct dicts.
    assert BoundaryOperator.registry is not CellUpdateBase.registry
    # And they hold their own subclasses, not each other's.
    for cls in BoundaryOperator.registry.values():
        assert cls not in CellUpdateBase.registry.values()
