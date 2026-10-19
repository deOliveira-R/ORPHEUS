r"""Values built once per host instance and then returned as the same object.

A downstream constructor (an operator or a space built FROM an upstream
object, rulings 22 and 34 of the numerics-untangling plan) must often hand
out one object per host instance: the faces of a frame, the lift over one
axis of a space. Identity tests rest on it (``R.codomain is E.domain``; a
Galerkin frame's test space ``is`` its trial space). The value is stored on
the host, in its instance ``__dict__`` (which also bypasses a frozen
dataclass's ``__setattr__``), under ONE private slot named for the module
that owns the values (``"_frame_faces"``, ``"_axis_marginal"``), so no cached
value sits under a public attribute name, where it could shadow or stand in
for a retired attribute. A cache is not a concept: the host does not know
what is cached on it, and the module that builds the value owns the slot.
"""

from __future__ import annotations

from collections.abc import Callable, Hashable
from typing import Any, TypeVar, cast

__all__ = ["cached_on"]

_V = TypeVar("_V")


def cached_on(host: object, slot: str, key: Hashable, build: Callable[[], _V]) -> _V:
    r"""The value ``key`` in ``host``'s private ``slot``, built once by ``build``.

    ``slot`` names the owning module's cache (``"_frame_faces"``,
    ``"_axis_marginal"``); it is a dict keyed by ``key``.
    """
    cache = cast("dict[Hashable, Any]", vars(host).setdefault(slot, {}))
    if key not in cache:
        cache[key] = build()
    return cache[key]
