r"""Self-registering subclass mixin.

Many extensible families in ORPHEUS — boundary conditions,
quadrature rules, cell-update strategies, solvers — share a common
shape: an abstract base type plus a name-keyed registry of concrete
subtypes that consumers select by string ID. Historically these
registries were maintained by hand: each concrete subtype's module
ran an explicit ``REGISTRY["my_key"] = MyConcrete`` statement.

Hand-maintained registries drift. A new subtype can be added without
the registry insert; a registry key can be reused by typo. The
mechanical fix is :class:`RegistryMixin`: subclasses register
themselves at class-creation time via the ``key=`` class-creation
keyword. Adding a subclass without giving it a registry slot is now
syntactically incomplete (the ``key=`` is right next to the base
class in the ``class`` statement), and duplicate keys raise
:class:`KeyError` at import time rather than silently shadowing.

This pattern follows Grand Report v3 §4. It is the ONE registry
primitive used across the operator-algebra families; every
extensible family in this codebase should mount on it rather than
maintaining a parallel dict.

Usage
-----

A registry root inherits :class:`RegistryMixin` and overrides
:meth:`_registry_base` to return itself; concrete subclasses pass
``key=...`` in their class statement::

    class BoundaryOperator(LinearOperatorMixin, RegistryMixin, ABC):
        registry: ClassVar[dict[str, type["BoundaryOperator"]]] = {}

        @classmethod
        def _registry_base(cls) -> type:
            return BoundaryOperator

        # ... rest of the abstract base ...

    class VacuumBoundaryOperator(BoundaryOperator, key="vacuum"):
        ...
    class SpecularBoundaryOperator(BoundaryOperator, key="reflective"):
        ...

After import, ``BoundaryOperator.registry`` contains
``{"vacuum": VacuumBoundaryOperator, "reflective":
SpecularBoundaryOperator, ...}`` and
``BoundaryOperator.create("vacuum")`` returns a fresh
``VacuumBoundaryOperator``.

Abstract intermediate classes (those that should not be reachable
via :meth:`create`) omit the ``key=`` kwarg — :meth:`__init_subclass__`
skips registration when ``key`` is ``None``.
"""

from __future__ import annotations

from typing import Any, ClassVar, Optional

__all__ = ["RegistryMixin"]


class RegistryMixin:
    """Self-registering subclass mixin.

    Each registry root re-declares its own :pydata:`registry` ClassVar
    (so different roots have independent dicts) and overrides
    :meth:`_registry_base` to return itself.

    Class attributes
    ----------------
    registry : dict[str, type]
        The registry — populated by ``__init_subclass__`` at class
        creation time. Each registry root MUST re-declare this as an
        empty dict; otherwise subtypes from sibling roots would leak
        into each other's registries via the inherited attribute.
    key : str | None
        The key under which the concrete subclass is registered.
        ``None`` on the base / abstract intermediate classes. Set
        automatically by :meth:`__init_subclass__` when ``key=`` is
        passed to the class statement.
    """

    registry: ClassVar[dict[str, type]] = {}
    key: ClassVar[Optional[str]] = None

    def __init_subclass__(
        cls,
        *,
        key: Optional[str] = None,
        **kwargs: Any,
    ) -> None:
        super().__init_subclass__(**kwargs)

        # Tolerate the ``@dataclass(slots=True)`` replacement pattern.
        # When ``@dataclass(slots=True)`` replaces a class, it creates
        # a new class object inheriting from the original — and the
        # new class re-runs ``__init_subclass__`` (with ``key=None``,
        # since the dataclass internal ``type(...)`` call doesn't
        # forward the original ``key=`` kwarg). Detect this case by
        # checking whether ``cls`` carries an inherited ``key`` and an
        # existing registry entry whose qualname matches ``cls``: the
        # replacement is a strict descendant of the original and the
        # registry entry should point at the replacement.
        if key is None:
            inherited_key = getattr(cls, "key", None)
            if inherited_key is None:
                # Genuine abstract-intermediate class: do not register.
                return
            base = cls._registry_base()
            existing = base.registry.get(inherited_key)
            if (
                existing is not None
                and existing is not cls
                and existing.__qualname__ == cls.__qualname__
                and existing.__module__ == cls.__module__
            ):
                base.registry[inherited_key] = cls
            return

        cls.key = key
        base = cls._registry_base()
        if key in base.registry:
            existing = base.registry[key]
            raise KeyError(
                f"duplicate registry key {key!r} on {base.__name__}: "
                f"already registered {existing.__name__}, attempted to "
                f"register {cls.__name__}"
            )
        base.registry[key] = cls

    @classmethod
    def _registry_base(cls) -> type:
        """Return the registry-root class for this hierarchy.

        Default returns ``cls`` itself; each registry-root must
        override this classmethod to return the root, so that
        subclasses many levels deep find the correct registry to
        register into.
        """
        return cls

    @classmethod
    def create(cls, key: str, **kwargs: Any) -> "RegistryMixin":
        """Construct a registered subclass instance by key.

        Parameters
        ----------
        key : str
            Registry key — one of the keys present in
            ``cls.registry``.
        **kwargs
            Forwarded to the concrete subclass constructor.

        Returns
        -------
        RegistryMixin
            A new instance of the registered concrete class.

        Raises
        ------
        KeyError
            If ``key`` is not registered. The error message lists
            the available keys to aid debugging.
        """
        base = cls._registry_base()
        if key not in base.registry:
            available = sorted(base.registry.keys())
            raise KeyError(
                f"unknown {base.__name__} key {key!r}; "
                f"available keys: {available}"
            )
        concrete = base.registry[key]
        return concrete(**kwargs)
