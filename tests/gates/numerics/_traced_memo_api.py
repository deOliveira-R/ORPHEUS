r"""One spelling of the traced memo's API for the M1–M5 gates (#405 P3, ``.claude/plans/reference_p3_spec.md``).

Every production name the gates touch is resolved HERE, at call time, off the module imported at the top
(lessons ``L100``: a battery that rebinds a module attribute must reach every row; layer DAG ruling 27: the module is
imported by statement, not by string): if the build lands a different module path or name, this file is the one
edit. The memoised clients (``CLIENTS``) live in ``test_traced_memo_clients.py``, their one consumer. The declared API (spec §0.2):

* ``traced_memo(function)`` / ``traced_memo(canonical={param: f})``: the decorator, returning a ``TracedMemo``
  with ``__call__``, ``lookup(*args, **kwargs) -> Verdict`` (no generation), ``key(*args, **kwargs) -> str``,
  ``__wrapped__`` and ``function_id``;
* ``Verdict = Hit | Absent | Stale | Corrupt`` (``Stale.reasons``, ``Corrupt.reason``);
* ``cache_root(path)`` and ``bypass()``: context managers;
* ``trace_call(f, *args, **kwargs) -> (value, Manifest)``: one traced call in THIS process;
* ``Manifest`` with ``functions`` ``(relpath, qualname, digest)``, ``modules`` ``(relpath, digest)``,
  ``distributions`` ``(name, version)``, ``python``, ``children`` ``(function id, key, payload digest)``;
* ``validate(manifest) -> tuple[str, ...]`` (empty: valid); ``function_digest``, ``skeleton_digest``;
* ``classify(filename) -> (kind, detail)``, read off ``origin(filename)``; ``Unpinnable``; ``Unencodable``;
* ``encode_payload(value, allowed) -> (tree, arrays)`` and ``decode_payload(tree, arrays, allowed)``;
* ``platform_tag()``.
"""
from __future__ import annotations

from typing import Any

from orpheus.numerics import traced_memo as _traced_memo

#: The module's name, for the sources the gates write (a synthetic package's or a subprocess's import line).
MODULE = _traced_memo.__name__


def module() -> Any:
    return _traced_memo


def name(attribute: str) -> Any:
    return getattr(module(), attribute)


def traced_memo(*args: Any, **kwargs: Any) -> Any:
    return name("traced_memo")(*args, **kwargs)


def trace_call(function: Any, *args: Any, **kwargs: Any) -> Any:
    return name("trace_call")(function, *args, **kwargs)


def validate(manifest: Any) -> tuple[str, ...]:
    return tuple(name("validate")(manifest))


def function_digest(source: str | bytes, qualname: str) -> Any:
    return name("function_digest")(source, qualname)


def skeleton_digest(source: str | bytes) -> Any:
    return name("skeleton_digest")(source)


#: The kind of each origin, as the gates spell it.
_ORIGIN_KIND = {"Dropped": "dropped", "Interpreter": "python", "Distribution": "distribution", "Source": "source"}


def classify(filename: str) -> tuple[str, str]:
    """``(kind, detail)`` of a traced file's origin (``kind`` one of dropped, python, distribution, source)."""
    found = name("origin")(filename, name("_Installed")())
    return _ORIGIN_KIND[type(found).__name__], found[0]


def encode_payload(value: Any, allowed: frozenset[str] = frozenset()) -> Any:
    """The payload of ``value``; ``allowed`` names the dataclasses it may hold (a function's return annotation)."""
    return name("encode_payload")(value, set(allowed))


def interpreter_identity(manifest: Any) -> str:
    """The interpreter a manifest pins (its one ``InterpreterPin``)."""
    (pin,) = manifest.of(name("InterpreterPin"))
    return pin.identity


def decode_payload(tree: Any, arrays: Any, allowed: set[str]) -> Any:
    return name("decode_payload")(tree, arrays, allowed)


def cache_root(path: Any) -> Any:
    return name("cache_root")(path)


def bypass() -> Any:
    return name("bypass")()


def platform_tag() -> str:
    return name("platform_tag")()


# ── an entry on disk: one ``.npz`` per entry, its JSON in the member ``__entry__`` ──────────────────────


def entry_files(root: Any) -> list[Any]:
    """Every entry under a store root (``<root>/<function id>/<key>.npz``), sorted."""
    from pathlib import Path

    return sorted(Path(root).rglob("*.npz"))


def entry_function(path: Any) -> str:
    """The function id an entry file belongs to."""
    return path.parent.name


def read_entry(path: Any) -> dict:
    """An entry's JSON (the key, the manifest rows, the payload tree and its digest)."""
    import json

    import numpy as np

    with np.load(path, allow_pickle=False) as npz:
        return json.loads(npz["__entry__"].tobytes())


def write_entry(path: Any, entry: dict) -> None:
    """Rewrite an entry's JSON in place, keeping its arrays (a corruption the payload digest must see)."""
    import json

    import numpy as np

    with np.load(path, allow_pickle=False) as npz:
        arrays = {name: npz[name] for name in npz.files if name != "__entry__"}
    with open(path, "wb") as handle:
        np.savez(handle, **arrays, __entry__=np.frombuffer(json.dumps(entry).encode(), np.uint8))


def manifest_rows(entry: dict, kind: str) -> list[list[Any]]:
    """The fields of every pin of ``kind`` (``"DefPin"``, ``"ModulePin"``, ``"ChildPin"``, ...) in an entry."""
    return [fields for tag, *fields in entry["manifest"] if tag == kind]


def verdict_kind(verdict: Any) -> str:
    """``"Hit" | "Absent" | "Stale" | "Corrupt"``: the verdict's member, by its class name."""
    kind = type(verdict).__name__
    if kind not in ("Hit", "Absent", "Stale", "Corrupt"):
        raise AssertionError(f"not a verdict: {verdict!r}")
    return kind



class GenerationCounter:
    """Counts the generations THIS process asks for, forked from a template or spawned fresh (``#602`` cut E: a
    miss no longer always starts an interpreter, so an interpreter count is not a generation count): every
    generation passes through ``TracedMemo._generate``, which this wraps. Observed on the route (``vv-principles``
    #26), not reported by the memo; a hit asks for none, which each user of this counter asserts."""

    def __init__(self, monkeypatch: Any) -> None:
        self.count = 0
        memo_class = name("TracedMemo")
        original = memo_class._generate

        def counting(memo: Any, job: Any) -> None:
            self.count += 1
            return original(memo, job)

        monkeypatch.setattr(memo_class, "_generate", counting)
