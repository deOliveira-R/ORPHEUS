r"""Content identity: one encoder for every value a persistent key covers.

A reference cache keys every entry on the CONTENT of the question it answers
(the materials, the geometry with its boundary laws, the discretisation), so
two values that are the same physics must give the same key in every process,
on every platform, and two values that differ in anything a consumer reads
must give different keys. Python's ``hash`` serves neither: it is salted per
process, and a value's ``__eq__`` is whatever its author wrote. This module is
the one definition of "same content":

* :func:`encode` turns a value into bytes, recursively, every chunk
  TYPE-TAGGED and LENGTH-PREFIXED, so no two different values share a byte
  stream (the encoding is injective on the values it admits);
* :func:`content_digest` is the blake2b-256 of those bytes;
* :class:`ContentIdentity` derives ``==`` and ``hash`` from the digest, so a
  type's equality and its cache key cannot drift (X4: one definition);
* :class:`FrozenMapping` is the one frozen, picklable mapping value the
  content types hold (a declaration's materials, a tag's parameters, a mesh's
  face laws).

The canonical forms, each with its reason
=========================================

**A real scalar is its value.** ``bool``, ``int``, ``float`` and the numpy
real scalars encode as one canonical IEEE-754 double, because Python's ``==``
already identifies them (``True == 1 == 1.0``; the user's ruling of
2026-10-02: the digest follows ``==``, so ``BC("albedo", {"albedo": 1})`` and
``{"albedo": 1.0}`` are one value). ``-0.0`` becomes ``+0.0`` (``-0.0 ==
0.0``). NaN is refused, with the path to it: NaN is not equal to itself, so a
value holding one has no content equality to encode. An integer beyond
:math:`2^{53}` is refused, since it has no exact double and ``==`` would then
separate it from the double it rounds to.

**A numeric array is its values on its shape.** Any boolean, integer or real
array encodes as its shape and its entries as canonical doubles, following
``np.array_equal`` (an integer material map equals its float twin); the
dtype is storage, not content. **A sparse matrix is the matrix**, not its
storage: zeros eliminated, duplicates summed, indices sorted, so an explicitly
stored zero equals its absence. A complex array or matrix is refused (no
content type holds one, and dropping the imaginary part would not be
injective).

**Containers keep the distinctions ``==`` keeps.** A ``tuple`` and a
``list`` carry different tags (``(1,) != [1]``); a ``Mapping`` is ordered by
its encoded keys, so insertion order is never content; a ``str`` key and an
``int`` key stay distinct (``"1" != 1``). An ``Enum`` member is its
module-qualified class and its name.

**An object is its schema and its parts.** A :class:`ContentIdentity` value,
and any frozen dataclass whose own equality is by value (``eq=True``),
encodes as its SCHEMA TAG followed by its parts, the dataclass fields with
``compare=True`` unless the class says otherwise. The tag is the class's
``__module__.__qualname__``, its ``__content_version__`` (a ``ClassVar[int]``,
default 1) and the ordered names of its parts. So adding, removing or
renaming a field, bumping the version (for a change in what a field MEANS),
and moving or renaming the class all change every digest of that class. That
is the user's ruling of 2026-10-01 (#405): a cache entry keyed under an older
schema misses, it never loads; there is no per-class ``__setstate__`` guard,
because the key covers the schema once for every class and every future
carve. A field that is not content (provenance, a derived quantity) is
declared ``field(compare=False)``, the one spelling of "not content".

**A part of an object is immutable.** A list, a dict, a set, a writeable
array or a sparse matrix whose arrays are writeable is refused as a part (it
could change after the object was keyed, under a cached digest), with
:class:`ContentlessError`, and so is a ``MappingProxyType``, a read-only view
over a ``dict`` its caller may still hold and write; a ``tuple``, a
``frozenset``, a :class:`FrozenMapping` and a read-only array are admitted.
The array check reads the array's own flag, so a read-only VIEW of a
writeable array passes it (#561). A value
handed to :func:`content_digest` directly is not a part and is not held to
this.

**Anything else has no content** and is refused with
:class:`ContentlessError`: a function, a plain object or a dataclass compared
by identity (``eq=False``), a mutable dataclass, an object array. A
:class:`ContentIdentity` value holding one has no digest; it is equal only to
itself, and it is unhashable.

**One limitation, by design.** The schema tag names a class by its module
and qualified name, so two classes defined under one name (in one function,
on two calls) share a tag. A class defined inside a function is therefore
not a persistent type, and nothing in the package defines one.

This module imports only :mod:`orpheus.numerics.scalars` (the one definition
of a real number, itself a leaf), so the data, geometry, mesh and numerics
layers all use it with no new layer edge.
"""

from __future__ import annotations

import dataclasses
import enum
import hashlib
import struct
import weakref
from collections.abc import Iterable, Iterator, Mapping
from typing import Any, ClassVar, Generic, TypeVar

import numpy as np

from orpheus.numerics.scalars import canonical_reals, exact_double

__all__ = [
    "ContentIdentity",
    "ContentlessError",
    "FrozenMapping",
    "constructor_arguments",
    "content_digest",
    "encode",
    "encode_exact",
    "name_digest",
]

#: Content digests are blake2b-256: wide enough that a collision between two
#: specifications is not a failure mode a cache needs to consider.
_DIGEST_SIZE = 32

#: The width, in bytes, of the digest a derived space NAME carries. A name is
#: a label and a diagnostic, not a persistent key (space identity is
#: structural since CS4c step 6), so eight bytes suffice.
_NAME_DIGEST_SIZE = 8

class ContentlessError(TypeError):
    r"""A value with no content identity reached the encoder.

    The message names the value's type and its PATH from the root value
    (``Materials.mixtures[3].SigT``), so the refusal says where the
    contentless part sits, not only that one exists.
    """


def _chunk(tag: bytes, payload: bytes) -> bytes:
    """One self-delimiting chunk: a one-byte tag, an 8-byte length, the bytes."""
    return tag + len(payload).to_bytes(8, "little") + payload


class _ByEquality:
    """Leaves encoded by their value under ``==``: one canonical double per real (``True == 1 == 1.0``,
    ``-0.0 == 0.0``), NaN and integers beyond :math:`2^{53}` refused; an array is its values on its shape, its
    dtype storage. The content identity of the module docstring."""

    def real(self, value: Any, path: str) -> bytes:
        return _chunk(b"R", struct.pack("<d", exact_double(value, path)))

    def entries(self, value: np.ndarray, path: str) -> bytes:
        return np.ascontiguousarray(canonical_reals(value, path), dtype="<f8").tobytes()


class _Exactly:
    """Leaves encoded exactly, with their types: ``True``, ``1``, ``1.0`` and ``np.float32(1)`` are four
    leaves, ``-0.0`` is not ``0.0``, a NaN is its bits; an array is its dtype and its bytes. A memo key's
    encoding (:func:`encode_exact`): a function can tell apart what ``==`` identifies, so a key that followed
    ``==`` could serve one call another's answer (#405 P3, the user's ruling of 2026-10-04, Q2)."""

    def real(self, value: Any, path: str) -> bytes:
        if isinstance(value, np.generic):
            return _chunk(b"n", value.dtype.str.encode() + b"|" + value.tobytes())
        if isinstance(value, bool):
            return _chunk(b"b", b"\x01" if value else b"\x00")
        if isinstance(value, int):
            return _chunk(b"i", str(value).encode())
        return _chunk(b"f", struct.pack("<d", value))

    def entries(self, value: np.ndarray, path: str) -> bytes:
        return value.dtype.str.encode() + b"|" + np.ascontiguousarray(value).tobytes()


_BY_EQUALITY = _ByEquality()
_EXACTLY = _Exactly()
_Leaves = _ByEquality | _Exactly


def _immutable(value: np.ndarray, path: str, frozen: bool) -> None:
    if frozen and value.flags.writeable:
        raise ContentlessError(
            f"{path}: a writeable array is not an immutable part (it could "
            f"change after its owner was keyed); store it read-only."
        )


def _array(value: np.ndarray, path: str, frozen: bool, leaves: _Leaves = _BY_EQUALITY) -> bytes:
    """Shape plus the entries as ``leaves`` encodes them; refuses a non-real array."""
    if value.dtype.kind not in "biuf":
        raise ContentlessError(
            f"{path}: an array of dtype {value.dtype} has no content "
            f"encoding (only boolean, integer and real arrays do)."
        )
    _immutable(value, path, frozen)
    shape = b"".join(int(n).to_bytes(8, "little") for n in value.shape)
    return _chunk(b"A", _chunk(b"s", shape) + _chunk(b"d", leaves.entries(value, path)))


def _sparse(value: Any, path: str, frozen: bool, leaves: _Leaves = _BY_EQUALITY) -> bytes:
    """The matrix, not its storage: CSR, zeros eliminated, indices sorted."""
    if frozen:
        for name in ("data", "indices", "indptr"):
            arr = getattr(value, name, None)
            if isinstance(arr, np.ndarray):
                _immutable(arr, f"{path}.{name}", frozen)
    m = value.tocsr(copy=True)
    # The data's own checks (dtype, the 2**53 bound, NaN) before any cast.
    _array(np.asarray(m.data), f"{path}.data", frozen=False)
    m.data = canonical_reals(np.asarray(m.data), f"{path}.data")
    m.sum_duplicates()
    m.eliminate_zeros()
    m.sort_indices()
    parts = (
        _array(np.asarray(m.shape), f"{path}.shape", False),
        _array(np.asarray(m.indptr, dtype=np.int64), f"{path}.indptr", False),
        _array(np.asarray(m.indices, dtype=np.int64), f"{path}.indices", False),
        _array(m.data, f"{path}.data", False, leaves),
    )
    return _chunk(b"C", b"".join(parts))


def _schema_tag(cls: type, names: tuple[str, ...]) -> bytes:
    """``module.qualname``, the content version and the ordered part names."""
    version = getattr(cls, "__content_version__", 1)
    text = f"{cls.__module__}.{cls.__qualname__}|v{version}|" + ",".join(names)
    return _chunk(b"h", text.encode())


def _object(value: Any, parts: tuple[tuple[str, Any], ...], path: str, leaves: _Leaves = _BY_EQUALITY) -> bytes:
    names = tuple(name for name, _ in parts)
    body = _schema_tag(type(value), names) + b"".join(
        _encode(part, f"{path}.{name}", True, leaves) for name, part in parts
    )
    return _chunk(b"O", body)


def _dataclass_parts(value: Any) -> tuple[tuple[str, Any], ...]:
    return tuple(
        (f.name, getattr(value, f.name))
        for f in dataclasses.fields(value)
        if f.compare
    )


def _mutable(kind: str, path: str) -> ContentlessError:
    return ContentlessError(
        f"{path}: a {kind} is not an immutable part (it could change after "
        f"its owner was keyed); store a tuple, a frozenset or a FrozenMapping."
    )


def _encode(value: Any, path: str, frozen: bool, leaves: _Leaves = _BY_EQUALITY) -> bytes:
    # Order: the content objects first (a FrozenMapping is a Mapping, and a
    # ContentIdentity dataclass is a dataclass), then the Enum before the
    # scalars (an IntEnum member is an int).
    if isinstance(value, ContentIdentity):
        return _object(value, value.content_parts(), path, leaves)
    if value is None:
        return _chunk(b"N", b"")
    if isinstance(value, enum.Enum):
        cls = type(value)
        return _chunk(
            b"E", f"{cls.__module__}.{cls.__qualname__}.{value.name}".encode()
        )
    if isinstance(value, (bool, int, float, np.bool_, np.integer, np.floating)):
        return leaves.real(value, path)
    if isinstance(value, str):
        return _chunk(b"S", value.encode())
    if isinstance(value, bytes):
        return _chunk(b"B", value)
    if isinstance(value, np.ndarray):
        return _array(value, path, frozen, leaves)
    if _is_sparse(value):
        return _sparse(value, path, frozen, leaves)
    if isinstance(value, tuple):
        return _chunk(
            b"T",
            b"".join(_encode(v, f"{path}[{i}]", frozen, leaves) for i, v in enumerate(value)),
        )
    if isinstance(value, list):
        if frozen:
            raise _mutable("list", path)
        return _chunk(
            b"L",
            b"".join(_encode(v, f"{path}[{i}]", frozen, leaves) for i, v in enumerate(value)),
        )
    if isinstance(value, Mapping):
        if frozen and not isinstance(value, _OwnedItems):
            raise _mutable(f"{type(value).__qualname__} mapping", path)
        items = sorted(
            (_encode(k, f"{path}<key>", frozen, leaves), _encode(v, f"{path}[{k!r}]", frozen, leaves))
            for k, v in value.items()
        )
        return _chunk(b"M", b"".join(k + v for k, v in items))
    if isinstance(value, (frozenset, set)):
        if frozen and isinstance(value, set):
            raise _mutable("set", path)
        elements = sorted(_encode(v, f"{path}{{}}", frozen, leaves) for v in value)
        return _chunk(b"F", b"".join(elements))
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        params = type(value).__dataclass_params__  # type: ignore[attr-defined]
        if not params.frozen:
            raise ContentlessError(
                f"{path}: {type(value).__qualname__} is a mutable dataclass, "
                f"whose content can change after it is keyed."
            )
        if not params.eq:
            raise ContentlessError(
                f"{path}: {type(value).__qualname__} is a dataclass compared by "
                f"identity (eq=False), so two of them with equal fields are "
                f"unequal values; its content is not its identity."
            )
        return _object(value, _dataclass_parts(value), path, leaves)
    raise ContentlessError(
        f"{path}: a {type(value).__module__}.{type(value).__qualname__} has no "
        f"content identity (it is neither a value type the encoder knows, a "
        f"frozen dataclass compared by value, nor a ContentIdentity). A "
        f"function or a plain object is compared by identity, which no "
        f"persistent key can carry."
    )


def _is_sparse(value: Any) -> bool:
    # scipy is imported lazily: the module must stay a leaf for the layers
    # that never build a sparse matrix.
    try:
        import scipy.sparse as sp
    except ImportError:  # pragma: no cover - scipy is a hard dependency
        return False
    return bool(sp.issparse(value))


def encode(value: Any) -> bytes:
    r"""The injective content encoding of ``value`` (module docstring).

    Raises
    ------
    ContentlessError
        When ``value`` holds a part with no content identity, or a mutable
        part inside an object.
    ValueError
        When it holds NaN, or an integer beyond :math:`2^{53}`.
    """
    return _encode(value, type(value).__qualname__, False)


def encode_exact(value: Any) -> bytes:
    r"""The EXACT encoding of ``value``: the walk of :func:`encode`, each leaf with its type and its bits
    (``1``, ``1.0`` and ``True`` are three values; ``-0.0`` is not ``0.0``; an array is its dtype and its bytes;
    a NaN is admitted as its bits). An object is still its schema and its parts. It is finer than content
    equality by design: the key of a call a function computes from its arguments (the traced memo's).

    Raises
    ------
    ContentlessError
        When ``value`` holds a part with no content identity, or a mutable part inside an object.
    """
    return _encode(value, type(value).__qualname__, False, _EXACTLY)


def content_digest(value: Any) -> bytes:
    """The blake2b-256 digest of :func:`encode` of ``value``."""
    if isinstance(value, ContentIdentity):
        return value.content_digest
    return hashlib.blake2b(encode(value), digest_size=_DIGEST_SIZE).digest()


def name_digest(value: Any) -> str:
    """The short hex digest a derived space NAME carries (eight bytes)."""
    return content_digest(value)[:_NAME_DIGEST_SIZE].hex()


#: Digests cached by object id, dropped when the object dies. Not stored on
#: the instance: a digest in ``__dict__`` would be pickled with it and read
#: back under a later schema: the stale key the schema-tag ruling of
#: 2026-10-01 (#405) forbids.
_DIGESTS: dict[int, bytes] = {}


class ContentIdentity:
    r"""Equality and hash derived from the content digest.

    A subclass is a frozen value: a frozen dataclass declared ``eq=False``
    (so the generated dunders cannot shadow these), or a class overriding
    :meth:`content_parts`. Two values are equal iff they are of one type and
    their digests agree; the hash is the digest's first eight bytes, the same
    in every process.

    A value whose content holds something contentless (a
    :class:`ContentlessError` from the encoder) is equal only to itself and
    unhashable: identity is the honest equality of a value with no content,
    and an identity hash would let it into a persistent key.
    """

    __slots__ = ()
    __content_version__: ClassVar[int] = 1

    def content_parts(self) -> tuple[tuple[str, Any], ...]:
        """The named parts the content is, in schema order.

        Default: the dataclass fields with ``compare=True``.
        """
        if not dataclasses.is_dataclass(self):
            raise TypeError(
                f"{type(self).__qualname__} is not a dataclass, so it must "
                f"override content_parts() to say what its content is."
            )
        return _dataclass_parts(self)

    @property
    def content_digest(self) -> bytes:
        """The blake2b-256 digest of this value's content (module docstring)."""
        key = id(self)
        cached = _DIGESTS.get(key)
        if cached is not None:
            return cached
        # Through ``_encode``, the one route, like every other value (a
        # field-less type would otherwise reach only the schema tag).
        digest = hashlib.blake2b(
            _encode(self, type(self).__qualname__, False),
            digest_size=_DIGEST_SIZE,
        ).digest()
        try:
            weakref.finalize(self, _DIGESTS.pop, key, None)
        except TypeError:  # a slotted class without __weakref__: no cache
            return digest
        _DIGESTS[key] = digest
        return digest

    def __reduce__(self):
        """Pickle through the constructor, so every law re-runs on load.

        A default unpickle restores the attributes without ``__post_init__``,
        so a value's arrays come back WRITEABLE and its admission checks never
        run; through the constructor they are read-only copies again and the
        laws hold. A pickle written under an older schema (a removed field)
        then fails to load with a ``TypeError`` naming the field, rather than
        loading as a value it is not.
        """
        if not dataclasses.is_dataclass(self):
            return super().__reduce__()
        return (_rebuild, (type(self), constructor_arguments(self)))

    def __eq__(self, other: object) -> bool:
        if self is other:
            return True
        if type(other) is not type(self):
            return NotImplemented
        try:
            return self.content_digest == other.content_digest  # type: ignore[attr-defined]
        except ContentlessError:
            return NotImplemented

    def __hash__(self) -> int:
        try:
            digest = self.content_digest
        except ContentlessError as err:
            raise TypeError(
                f"unhashable {type(self).__qualname__}: it has no content "
                f"identity ({err})"
            ) from err
        return int.from_bytes(digest[:8], "little", signed=True)


def constructor_arguments(value: Any) -> dict[str, Any]:
    """A dataclass value's constructor form, its ``init`` fields by name: what rebuilds it through its
    constructor (:func:`_rebuild`), so its laws re-run. The one spelling for every route that carries a value
    by its construction (a pickle, the traced memo's process boundary and its payload)."""
    return {f.name: getattr(value, f.name) for f in dataclasses.fields(value) if f.init}


def _rebuild(cls: type, init: dict[str, Any]) -> Any:
    """The unpickling constructor of a dataclass value (``ContentIdentity.__reduce__``)."""
    return cls(**init)


_K = TypeVar("_K")
_V = TypeVar("_V")


class FrozenMapping(ContentIdentity, Mapping[_K, _V], Generic[_K, _V]):
    r"""A frozen, picklable mapping VALUE with content identity.

    Keeps the declared order for iteration (order is behaviour: a caller may
    document it) while its content, and so its equality and hash, is
    order-free. Equal only to another mapping of its own type with equal
    items, never to a plain ``dict``. The one spelling of a frozen mapping
    the content types hold; a subclass that adds behaviour (a mesh's face
    laws) keeps the storage.
    """

    __slots__ = ("_items", "_index", "__weakref__")

    _items: tuple[tuple[_K, _V], ...]
    _index: dict[_K, _V]

    def __init__(self, items: "Mapping[_K, _V] | Iterable[tuple[_K, _V]]" = ()) -> None:
        pairs = tuple(items.items() if isinstance(items, Mapping) else items)
        index = dict(pairs)
        if len(index) != len(pairs):
            raise ValueError(f"{type(self).__name__}: a key is given twice in {pairs!r}")
        object.__setattr__(self, "_items", pairs)
        object.__setattr__(self, "_index", index)

    def __getitem__(self, key: _K) -> _V:
        return self._index[key]

    def __iter__(self) -> Iterator[_K]:
        return (key for key, _ in self._items)

    def __len__(self) -> int:
        return len(self._items)

    def __setattr__(self, name: str, value: object) -> None:
        raise AttributeError(f"{type(self).__name__} is immutable")

    def __reduce__(self):
        return (type(self), (self._items,))

    def __repr__(self) -> str:
        body = ", ".join(f"{k!r}: {v!r}" for k, v in self._items)
        return f"{type(self).__name__}({{{body}}})"

    def content_parts(self) -> tuple[tuple[str, Any], ...]:
        return (("items", _OwnedItems(self._index)),)


class _OwnedItems(Mapping[_K, _V]):
    """The items of a :class:`FrozenMapping`, as its content part.

    A read-only view over storage only the frozen mapping holds, so the
    frozen encoder admits it. A ``MappingProxyType`` is refused instead: it is
    a view over a ``dict`` its caller may still hold and write after the
    owner was keyed, so its content is not fixed.
    """

    __slots__ = ("_index",)

    def __init__(self, index: dict[_K, _V]) -> None:
        self._index = index

    def __getitem__(self, key: _K) -> _V:
        return self._index[key]

    def __iter__(self) -> Iterator[_K]:
        return iter(self._index)

    def __len__(self) -> int:
        return len(self._index)
