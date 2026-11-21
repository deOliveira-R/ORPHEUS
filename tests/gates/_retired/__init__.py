"""The row types of the shared retirement engine (specification ``ontology_w0_spec.md`` section 1.5).

Five ``NamedTuple``s: four row types (:class:`Retired`, :class:`Moved`, :class:`RetiredName`, :class:`RetiredModule`), the members of
a rows module's ``ROWS``, and the entry type :class:`Kept` of its ``KEPT``. The walk (``_walk.walk``) dispatches on the type.

A ``NamedTuple`` cannot state an invariant at construction (``__post_init__`` is never called on one, measured on Python 3.14), so
the invalid rows (an empty field, a ``file`` or ``home`` that does not end in ``.py``, an empty ``keeper``, a row listed twice, a
``RetiredName`` that is the last component of a ``Kept`` qualname) are refused by the walk, as a ``GateError``.

Rows live in the modules ``rows_*.py`` of this package, each defining ``ROWS`` and ``KEPT`` and nothing else the engine reads; the
engine reads them from this directory and never from the root it scans. The package is excluded from the scan (its rows name the
retired symbols by design), so every done-when ``git grep`` over a retired name excludes ``tests/gates/_retired``.
"""

from __future__ import annotations

from typing import NamedTuple


class Retired(NamedTuple):
    """The definition ``qualname`` is gone from ``file`` and nothing imports, reads or sets it by its definition site."""

    file: str
    qualname: str


class Moved(NamedTuple):
    """As :class:`Retired`, and ``home_name`` (a top-level name) is defined at the repository path ``home``."""

    file: str
    qualname: str
    home: str
    home_name: str


class RetiredName(NamedTuple):
    """A bare identifier that appears nowhere in the five trees; used only where the name has no surviving homonym."""

    name: str


class RetiredModule(NamedTuple):
    """A dotted module path that no longer exists: no import resolves to it, and neither its ``.py`` file nor its package does."""

    module: str


class Kept(NamedTuple):
    """The entry type of ``KEPT``: the definition ``(file, qualname)`` must STILL exist; ``keeper`` is the wave or consumer that keeps it."""

    file: str
    qualname: str
    keeper: str


Row = Retired | Moved | RetiredName | RetiredModule
