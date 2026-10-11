r"""One spelling of the traced memo's API for the M1–M5 gates (#405 P3, ``.claude/plans/reference_p3_spec.md``).

Every production name the gates touch is resolved HERE, at call time (lessons ``L100``: a battery that
rebinds a module attribute must reach every row): if the build lands a different module path or name, this
file is the one edit. The declared API (spec §0.2):

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

import importlib
from typing import Any

MODULE = "orpheus.numerics.traced_memo"


def module() -> Any:
    return importlib.import_module(MODULE)


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
    found = name("origin")(filename)
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


# ── the P3 clients (step 4): where each memo is bound ────────────────────────────

#: (module, dotted attribute) of each memoised client: the characteristic reference's reading (its ``evaluate``, keyed on
#: the derivation's content and the observable) and its solve (its ``solve``, keyed on the derivation, a child entry of
#: every reading since #592). The trajectory resolvent's three clients (its two multi-region solvers and its reading)
#: were deleted with the family in step (e2) of the characteristic-reference campaign. The
#: exact infinite medium is NOT a client: its reading costs less than an interpreter start (the user's ruling of
#: 2026-10-04, ``reference_cache.md`` "P3 specification ruled", Q5).
CLIENTS = {
    "characteristic_reading": ("orpheus.derivations.continuous.characteristic.reference", "CharacteristicDerivation.evaluate"),
    "characteristic_solve": ("orpheus.derivations.continuous.characteristic.reference", "CharacteristicDerivation.solve"),
}


def spelled_function_id(name_: str) -> str:
    """The function id a client WOULD have (``module:qualname``, ``TracedMemo.function_id``'s spelling), read off
    :data:`CLIENTS` without importing the attribute: a gate can then name a client the tree does not yet bind and
    red on its claim, not on an ``AttributeError``. ``test_m4_1`` checks it equals each bound client's own id."""
    module_name, dotted = CLIENTS[name_]
    return f"{module_name}:{dotted}"


def client(name_: str) -> Any:
    module_name, dotted = CLIENTS[name_]
    value: Any = importlib.import_module(module_name)
    for part in dotted.split("."):
        value = getattr(value, part)
    return value


def function_id(name_: str) -> str:
    return client(name_).function_id


class SpawnCounter:
    """Counts the interpreters started from THIS process (a ``subprocess.Popen`` of ``sys.executable``): the
    route a generation takes, observed on the route (``vv-principles`` #26), not reported by the memo."""

    def __init__(self, monkeypatch: Any) -> None:
        import subprocess
        import sys

        self.count = 0
        original = subprocess.Popen.__init__

        def counting(popen: Any, args: Any, *a: Any, **k: Any) -> None:
            if isinstance(args, (list, tuple)) and args and str(args[0]) == sys.executable:
                self.count += 1
            original(popen, args, *a, **k)

        monkeypatch.setattr(subprocess.Popen, "__init__", counting)
