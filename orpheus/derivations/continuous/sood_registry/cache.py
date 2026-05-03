r"""Persistent solver-output cache for Sood-family benchmark cases.

The :mod:`sood_registry` package owns the **truth values** for Sood
benchmark cases (k_eff, k_inf, critical dimensions, flux ratios). Those
values are hard-coded data — no caching infrastructure needed; the
registry IS the cache for ground truth.

This module addresses the **complementary** problem: when a test asks a
slow solver to *recompute* its answer on a Sood case (e.g.
``solve_fn_sphere_bare_critical(c=1.30, n_modes=15)`` or a Variant-α
cross-check at fine quadrature), we want to memoize the result across
Python sessions so subsequent test runs are fast.

Design
------

* **Pure addition.** Existing tests do not import or call the cache.
  Cache adoption is opt-in: a test that wants memoization wraps a
  solver call with the :func:`sood_cache` decorator (or with the
  context-managed :class:`SoodResultCache` for finer control).
* **Hash-stable keys.** The cache key is a tuple
  ``(solver_qualified_name, params_hash)`` where ``params_hash`` is a
  SHA-256 over a canonicalized JSON representation of the keyword
  arguments. The same call across Python sessions produces the same
  key bit-for-bit.
* **Versioning.** Solver output may change if the solver implementation
  changes. Cache entries carry a ``version`` string; on read, a
  mismatch is treated as a miss and the entry is rewritten. The
  default version is auto-detected from the current git HEAD SHA;
  callers may pin an explicit string instead.
* **Storage.** Pickle files under ``.cache/sood_registry/`` at the
  project root. Each entry lives in its own file
  ``<solver-slug>_<hash>.pkl`` so a corrupt entry does not poison the
  whole cache. Pickle is chosen over JSON because solver result types
  are dataclasses carrying numpy arrays — JSON would need a
  serialization shim per type. Inspecting a cache entry from a
  debugger is one ``pickle.load`` call.
* **API ergonomics.** The decorator form is the recommended path; the
  class form is for tests that need to invalidate or inspect entries
  programmatically.

The cache is local-machine state. The cache directory is gitignored;
nothing should ever be committed.

Caveats
-------

* **Single-process safety only.** No file locking. Concurrent test
  workers writing to the same key may produce a "last writer wins"
  outcome, which is harmless for this use (deterministic solvers
  produce identical bytes) but is not transactional.
* **Cache poisoning.** If a solver has a bug and the cache stores a
  wrong answer, the wrong answer survives across runs until either
  the version bumps or the user calls :func:`clear_cache`. Tests that
  actively verify a solver against an external reference are
  themselves the protection — a poisoned cache fails the test.

Usage — decorator
-----------------

.. code-block:: python

    from orpheus.derivations.continuous.sood_registry.cache import sood_cache

    @sood_cache()
    def fn_sphere_at_n(c: float, n_modes: int):
        return solve_fn_sphere_bare_critical(c=c, n_modes=n_modes)

    res = fn_sphere_at_n(c=1.30, n_modes=15)  # first call: computes
    res = fn_sphere_at_n(c=1.30, n_modes=15)  # second call: cache hit

Usage — class form
------------------

.. code-block:: python

    from orpheus.derivations.continuous.sood_registry.cache import (
        SoodResultCache,
    )

    cache = SoodResultCache()
    res = cache.get_or_compute(
        solver_name="solve_fn_sphere_bare_critical",
        params={"c": 1.30, "n_modes": 15},
        compute=lambda: solve_fn_sphere_bare_critical(c=1.30, n_modes=15),
    )
"""
from __future__ import annotations

import functools
import hashlib
import json
import os
import pickle
import re
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Iterator, Mapping

import numpy as np


# ═══════════════════════════════════════════════════════════════════
# Cache directory + version detection
# ═══════════════════════════════════════════════════════════════════


def _project_root() -> Path:
    """Return the project root (the directory containing pyproject.toml).

    Walks upward from this file's location. Does not depend on cwd —
    sub-agents and pytest both reset cwd between calls, so any path
    derived from cwd would be unreliable.
    """
    here = Path(__file__).resolve()
    for parent in (here, *here.parents):
        if (parent / "pyproject.toml").exists():
            return parent
    # Fallback: this file lives at
    # ``orpheus/derivations/continuous/sood_registry/cache.py`` —
    # five parents up is the project root.
    return here.parents[4]


_DEFAULT_CACHE_DIR = _project_root() / ".cache" / "sood_registry"


def _detect_git_sha() -> str:
    """Return the current git HEAD SHA (short form), or ``"nogit"``.

    Used as the default cache version. A version mismatch on read
    invalidates the entry, so a code change that bumps the SHA
    naturally retires stale cache entries.
    """
    try:
        out = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=str(_project_root()),
            stderr=subprocess.DEVNULL,
            text=True,
        ).strip()
        return out or "nogit"
    except (subprocess.SubprocessError, FileNotFoundError, OSError):
        return "nogit"


# ═══════════════════════════════════════════════════════════════════
# Hash + canonicalization
# ═══════════════════════════════════════════════════════════════════


def _canonicalize(obj: Any) -> Any:
    """Recursively convert ``obj`` to a JSON-serializable canonical form.

    Handles numpy scalars/arrays, tuples, sets, frozensets, mappings,
    and arbitrary nested structures. Non-serializable values raise
    :class:`TypeError` so callers see the failure at hash time, not at
    cache-read time.
    """
    if isinstance(obj, (str, bool, type(None))):
        return obj
    if isinstance(obj, (int, float)):
        return obj
    if isinstance(obj, np.ndarray):
        return {"__nd__": True, "shape": list(obj.shape),
                "dtype": str(obj.dtype), "data": obj.tolist()}
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, Mapping):
        # Sort keys for deterministic ordering.
        return {str(k): _canonicalize(obj[k]) for k in sorted(obj, key=str)}
    if isinstance(obj, (list, tuple)):
        return [_canonicalize(x) for x in obj]
    if isinstance(obj, (set, frozenset)):
        return {"__set__": sorted([_canonicalize(x) for x in obj], key=repr)}
    raise TypeError(
        f"sood_cache: cannot canonicalize object of type {type(obj).__name__}"
    )


def _hash_params(params: Mapping[str, Any]) -> str:
    """Deterministic 16-char SHA-256 prefix of canonicalized params.

    The hash collapses a nested dict of solver kwargs to a fixed-width
    string suitable for a filename. 16 chars of SHA-256 give 2^64
    keyspace — collisions are not a practical concern for our scale
    (a few thousand cache entries at most).
    """
    canon = _canonicalize(dict(params))
    raw = json.dumps(canon, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()[:16]


_SLUG_RE = re.compile(r"[^A-Za-z0-9_]+")


def _solver_slug(solver_name: str) -> str:
    """Convert a qualified solver name to a filesystem-safe slug."""
    slug = _SLUG_RE.sub("_", solver_name).strip("_")
    return slug or "solver"


# ═══════════════════════════════════════════════════════════════════
# Cache entry on-disk format
# ═══════════════════════════════════════════════════════════════════


@dataclass(frozen=True)
class CacheEntry:
    """Pickled metadata + payload of a single cache entry.

    Stored verbatim to disk. The wrapper dataclass ensures we always
    know which version produced a given payload, regardless of how the
    payload type evolves.
    """

    solver_name: str
    params: dict[str, Any]
    version: str
    timestamp: float
    payload: Any


# ═══════════════════════════════════════════════════════════════════
# Class-form cache
# ═══════════════════════════════════════════════════════════════════


class SoodResultCache:
    """File-backed solver-output cache for Sood benchmark cases.

    Parameters
    ----------
    version :
        Cache-invalidation tag. Default ``"auto"`` — auto-detect via
        ``git rev-parse --short HEAD``; callers may pin an explicit
        string for reproducible offline runs. Falls back to
        ``"nogit"`` outside a git checkout.
    cache_dir :
        Directory holding pickled entries. Defaults to
        ``<project-root>/.cache/sood_registry/``. Created lazily on
        first write.

    Notes
    -----
    The cache is single-process. It is correct under the SECT
    (single-execution-context, terminating) test runner that ORPHEUS
    uses; concurrent writers to the same key will produce a
    "last-writer-wins" outcome (harmless for deterministic solvers).
    """

    def __init__(
        self,
        version: str = "auto",
        cache_dir: Path | str | None = None,
    ) -> None:
        self.version: str = _detect_git_sha() if version == "auto" else version
        self.cache_dir: Path = (
            Path(cache_dir) if cache_dir is not None else _DEFAULT_CACHE_DIR
        )

    # -----------------------------------------------------------------
    # Path helpers
    # -----------------------------------------------------------------

    def _entry_path(self, solver_name: str, params_hash: str) -> Path:
        return self.cache_dir / f"{_solver_slug(solver_name)}_{params_hash}.pkl"

    # -----------------------------------------------------------------
    # Read / write
    # -----------------------------------------------------------------

    def get(
        self,
        solver_name: str,
        params: Mapping[str, Any],
    ) -> Any | None:
        """Return the cached payload, or ``None`` on miss / version mismatch.

        A corrupt entry (unpicklable bytes, malformed wrapper) is
        treated as a miss and the file is removed so the next call
        recomputes cleanly. Filesystem read errors propagate.
        """
        path = self._entry_path(solver_name, _hash_params(params))
        if not path.exists():
            return None
        try:
            with path.open("rb") as fh:
                entry = pickle.load(fh)
        except (pickle.UnpicklingError, EOFError, AttributeError, TypeError):
            # Corrupt entry — evict and return miss.
            try:
                path.unlink()
            except OSError:
                pass
            return None
        if not isinstance(entry, CacheEntry):
            try:
                path.unlink()
            except OSError:
                pass
            return None
        if entry.version != self.version:
            return None
        return entry.payload

    def put(
        self,
        solver_name: str,
        params: Mapping[str, Any],
        payload: Any,
    ) -> Path:
        """Write a payload to the cache. Returns the entry path."""
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        path = self._entry_path(solver_name, _hash_params(params))
        entry = CacheEntry(
            solver_name=solver_name,
            params=dict(params),
            version=self.version,
            timestamp=time.time(),
            payload=payload,
        )
        # Write to a sibling tempfile then rename — avoids leaving a
        # half-written .pkl if the process is killed mid-write.
        tmp = path.with_suffix(path.suffix + ".tmp")
        with tmp.open("wb") as fh:
            pickle.dump(entry, fh, protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(tmp, path)
        return path

    def get_or_compute(
        self,
        solver_name: str,
        params: Mapping[str, Any],
        compute: Callable[[], Any],
    ) -> Any:
        """Return the cached payload, computing+storing it on miss."""
        cached = self.get(solver_name, params)
        if cached is not None:
            return cached
        payload = compute()
        self.put(solver_name, params, payload)
        return payload

    # -----------------------------------------------------------------
    # Maintenance / introspection
    # -----------------------------------------------------------------

    def clear(self) -> int:
        """Delete every cache entry. Returns the number of files removed.

        Idempotent: clearing an empty cache is a no-op. Used by tests
        to enforce a clean slate.
        """
        if not self.cache_dir.exists():
            return 0
        n = 0
        for entry in self.cache_dir.iterdir():
            if entry.is_file() and entry.suffix in (".pkl", ".tmp"):
                try:
                    entry.unlink()
                    n += 1
                except OSError:
                    pass
        return n

    def info(self) -> list[dict[str, Any]]:
        """List every cache entry as a dict: solver_name, params,
        version, timestamp, path. Useful for debugging.

        Skips corrupt entries silently (they will be evicted on next
        :meth:`get`).
        """
        if not self.cache_dir.exists():
            return []
        out: list[dict[str, Any]] = []
        for path in sorted(self.cache_dir.glob("*.pkl")):
            try:
                with path.open("rb") as fh:
                    entry = pickle.load(fh)
            except Exception:
                continue
            if not isinstance(entry, CacheEntry):
                continue
            out.append({
                "solver_name": entry.solver_name,
                "params": entry.params,
                "version": entry.version,
                "timestamp": entry.timestamp,
                "path": str(path),
            })
        return out

    def __iter__(self) -> Iterator[dict[str, Any]]:
        return iter(self.info())


# ═══════════════════════════════════════════════════════════════════
# Decorator-form cache
# ═══════════════════════════════════════════════════════════════════


def sood_cache(
    *,
    solver_name: str | None = None,
    version: str = "auto",
    cache_dir: Path | str | None = None,
    key_args: list[str] | None = None,
) -> Callable[[Callable[..., Any]], Callable[..., Any]]:
    """Persist a solver function's outputs across Python sessions.

    Parameters
    ----------
    solver_name :
        Override the cache key's ``solver_name`` component. Defaults
        to the wrapped function's qualified name (``module.qualname``).
        Pin this explicitly when the function lives in a test file
        that may be renamed without changing the underlying solver.
    version :
        Cache-invalidation tag — ``"auto"`` for git SHA, or any
        explicit string. See :class:`SoodResultCache`.
    cache_dir :
        Override the default cache directory.
    key_args :
        Whitelist of kwarg names that contribute to the cache key.
        If ``None`` (default), every kwarg is hashed. Use this to
        ignore non-load-bearing kwargs (e.g. a ``verbose`` flag) that
        should not invalidate the cache.

    Returns
    -------
    decorator
        Wraps a function whose return value is picklable. The wrapped
        function is callable with the same signature as the original.

    Notes
    -----
    The decorator hashes only **keyword arguments**. Positional
    arguments are forwarded to the wrapped function but do NOT
    contribute to the cache key. This is the intentional ergonomic
    rule: tests should call cached solvers with kwargs only, so the
    cache hash is unambiguous regardless of how the caller phrases
    the call.

    Examples
    --------
    >>> @sood_cache()
    ... def slow_solver(*, c, n_modes):
    ...     return expensive_computation(c, n_modes)
    >>> slow_solver(c=1.30, n_modes=10)  # computes
    >>> slow_solver(c=1.30, n_modes=10)  # cache hit
    """

    def decorator(fn: Callable[..., Any]) -> Callable[..., Any]:
        cache = SoodResultCache(version=version, cache_dir=cache_dir)
        sname = solver_name or f"{fn.__module__}.{fn.__qualname__}"

        @functools.wraps(fn)
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            if args:
                raise TypeError(
                    f"sood_cache-wrapped function {sname!r} must be called "
                    f"with keyword arguments only (got {len(args)} positional)."
                )
            params = (
                {k: kwargs[k] for k in key_args if k in kwargs}
                if key_args is not None
                else dict(kwargs)
            )
            return cache.get_or_compute(
                solver_name=sname,
                params=params,
                compute=lambda: fn(**kwargs),
            )

        # Expose the underlying cache for test inspection / clearing.
        wrapper.cache = cache  # type: ignore[attr-defined]
        wrapper.solver_name = sname  # type: ignore[attr-defined]
        return wrapper

    return decorator


# ═══════════════════════════════════════════════════════════════════
# Module-level conveniences
# ═══════════════════════════════════════════════════════════════════


def clear_cache(cache_dir: Path | str | None = None) -> int:
    """Clear the default cache (or one at ``cache_dir`` if given).

    Returns the number of entries removed. Convenience wrapper over
    :meth:`SoodResultCache.clear` for tests / CLI use.
    """
    return SoodResultCache(cache_dir=cache_dir).clear()


def cache_info(
    cache_dir: Path | str | None = None,
) -> list[dict[str, Any]]:
    """List entries in the default cache (or ``cache_dir`` if given)."""
    return SoodResultCache(cache_dir=cache_dir).info()


__all__ = [
    "CacheEntry",
    "SoodResultCache",
    "sood_cache",
    "clear_cache",
    "cache_info",
]
