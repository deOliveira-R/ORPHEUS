r"""Foundation tests for the Sood-registry solver-output cache.

The cache is a pure-addition convenience for tests that wrap a slow
solver call with :func:`sood_cache` (decorator) or
:class:`SoodResultCache` (class form). The tests below pin the
load-bearing invariants:

* Round-trip: write then read returns the stored payload.
* Miss: an unseen ``(solver_name, params)`` returns ``None``.
* Version invalidation: writing at version X and reading at version
  Y returns a miss.
* Hash stability: equal kwargs (modulo order, modulo numpy/Python
  scalar type) produce the same cache key.
* ``clear()`` empties the cache.
* ``info()`` lists every stored entry with metadata.
* Decorator integration: a wrapped function recomputes once, then
  hits the cache on every subsequent call.
* Production-solver integration: cache a real F_N sphere solver call
  and verify byte-equal payload on second invocation (L1 smoke).

All tests use a per-test temporary cache directory; the production
``.cache/sood_registry/`` is never touched. Cache tests are
``@pytest.mark.foundation`` because they verify a software invariant
of the cache primitive, not a physics claim.
"""
from __future__ import annotations

import pickle
import tempfile
from pathlib import Path

import numpy as np
import pytest

from orpheus.derivations.continuous.sood_registry import (
    PUA_1_0_IN,
    SoodResultCache,
    cache_info,
    clear_cache,
    sood_cache,
)
from orpheus.derivations.continuous.sood_registry.cache import (
    CacheEntry,
    _hash_params,
)


# ═══════════════════════════════════════════════════════════════════
# Fixtures
# ═══════════════════════════════════════════════════════════════════


@pytest.fixture
def tmp_cache_dir(tmp_path: Path) -> Path:
    """Per-test scratch cache dir. pytest cleans up automatically."""
    d = tmp_path / "sood_cache"
    d.mkdir()
    return d


# ═══════════════════════════════════════════════════════════════════
# Round-trip / miss / clear
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_cache_write_then_read_returns_stored_payload(tmp_cache_dir: Path):
    """Basic round trip: ``put`` then ``get`` returns the same payload."""
    cache = SoodResultCache(version="v1", cache_dir=tmp_cache_dir)
    payload = {"r_c": 2.4248249802, "coefficients": np.array([1.0, 0.5, 0.25])}
    cache.put("solve_fn_sphere", {"c": 1.30, "n_modes": 10}, payload)
    got = cache.get("solve_fn_sphere", {"c": 1.30, "n_modes": 10})
    assert got is not None
    assert got["r_c"] == payload["r_c"]
    np.testing.assert_array_equal(got["coefficients"], payload["coefficients"])


@pytest.mark.foundation
def test_cache_miss_returns_none(tmp_cache_dir: Path):
    """An unseen (solver, params) returns ``None`` cleanly."""
    cache = SoodResultCache(version="v1", cache_dir=tmp_cache_dir)
    assert cache.get("never_called", {"x": 1}) is None
    cache.put("solve_fn_sphere", {"c": 1.30, "n_modes": 10}, {"r": 1.0})
    # Different params = miss.
    assert cache.get("solve_fn_sphere", {"c": 1.30, "n_modes": 11}) is None
    # Different solver name = miss.
    assert cache.get("solve_other", {"c": 1.30, "n_modes": 10}) is None


@pytest.mark.foundation
def test_cache_clear_empties_cache(tmp_cache_dir: Path):
    """``clear()`` removes every entry and is idempotent."""
    cache = SoodResultCache(version="v1", cache_dir=tmp_cache_dir)
    cache.put("a", {"x": 1}, "alpha")
    cache.put("b", {"x": 2}, "beta")
    assert len(cache.info()) == 2
    n = cache.clear()
    assert n == 2
    assert cache.info() == []
    # Idempotent: clearing an empty cache is a no-op.
    assert cache.clear() == 0


@pytest.mark.foundation
def test_cache_info_lists_entries_with_metadata(tmp_cache_dir: Path):
    """``info()`` reports solver_name, params, version, timestamp, path."""
    cache = SoodResultCache(version="v1", cache_dir=tmp_cache_dir)
    cache.put("solver_a", {"c": 1.30, "n": 10}, {"r": 2.42})
    info = cache.info()
    assert len(info) == 1
    entry = info[0]
    assert entry["solver_name"] == "solver_a"
    assert entry["params"] == {"c": 1.30, "n": 10}
    assert entry["version"] == "v1"
    assert isinstance(entry["timestamp"], float)
    assert Path(entry["path"]).exists()


# ═══════════════════════════════════════════════════════════════════
# Version invalidation
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_cache_version_mismatch_is_a_miss(tmp_cache_dir: Path):
    """Writing at version X and reading at version Y returns ``None``.

    Version invalidation is the load-bearing protection against stale
    entries surviving a solver-implementation change. The test models
    a code-change-induced SHA bump.
    """
    write_cache = SoodResultCache(version="sha-aaaaa", cache_dir=tmp_cache_dir)
    write_cache.put("solve_fn_sphere", {"c": 1.30, "n_modes": 10}, "old-result")
    # Same dir, different version.
    read_cache = SoodResultCache(version="sha-bbbbb", cache_dir=tmp_cache_dir)
    assert read_cache.get("solve_fn_sphere", {"c": 1.30, "n_modes": 10}) is None
    # The original-version cache still hits (the file is unchanged).
    assert (
        write_cache.get("solve_fn_sphere", {"c": 1.30, "n_modes": 10})
        == "old-result"
    )


# ═══════════════════════════════════════════════════════════════════
# Hash stability
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_hash_stable_across_kwarg_order():
    """``{a:1,b:2}`` and ``{b:2,a:1}`` produce the same cache key."""
    h1 = _hash_params({"c": 1.30, "n_modes": 10})
    h2 = _hash_params({"n_modes": 10, "c": 1.30})
    assert h1 == h2


@pytest.mark.foundation
def test_hash_stable_for_equivalent_numpy_and_python_scalars():
    """``np.float64(1.30)`` and ``1.30`` produce the same key.

    Solver tests pass either flavor depending on how the c-value was
    computed; the cache must not split the entry on this distinction.
    """
    h_python = _hash_params({"c": 1.30, "n": 10})
    h_numpy = _hash_params({"c": np.float64(1.30), "n": np.int64(10)})
    assert h_python == h_numpy


@pytest.mark.foundation
def test_hash_distinguishes_different_params():
    """Different params produce different keys."""
    h1 = _hash_params({"c": 1.30, "n_modes": 10})
    h2 = _hash_params({"c": 1.30, "n_modes": 11})  # one int different
    h3 = _hash_params({"c": 1.40, "n_modes": 10})  # one float different
    assert h1 != h2 != h3 != h1


@pytest.mark.foundation
def test_hash_handles_numpy_arrays():
    """Numpy arrays are canonicalized to a stable representation."""
    h1 = _hash_params({"xs": np.array([1.0, 2.0, 3.0])})
    h2 = _hash_params({"xs": np.array([1.0, 2.0, 3.0])})
    assert h1 == h2
    h3 = _hash_params({"xs": np.array([1.0, 2.0, 3.1])})
    assert h1 != h3


# ═══════════════════════════════════════════════════════════════════
# Decorator form
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_decorator_first_call_computes_second_hits(tmp_cache_dir: Path):
    """A decorated function computes once, then hits cache forever."""
    counter = {"n": 0}

    @sood_cache(version="v1", cache_dir=tmp_cache_dir)
    def expensive(*, c: float, n_modes: int) -> dict:
        counter["n"] += 1
        return {"r": c * n_modes, "n": n_modes}

    r1 = expensive(c=1.30, n_modes=10)
    r2 = expensive(c=1.30, n_modes=10)
    r3 = expensive(c=1.30, n_modes=10)
    assert r1 == r2 == r3
    assert counter["n"] == 1, "Cache should have prevented recomputation"

    # Different params → recompute.
    r4 = expensive(c=1.30, n_modes=15)
    assert counter["n"] == 2
    assert r4 == {"r": 1.30 * 15, "n": 15}


@pytest.mark.foundation
def test_decorator_rejects_positional_arguments(tmp_cache_dir: Path):
    """Positional args are rejected — kwargs only for unambiguous hashing."""

    @sood_cache(version="v1", cache_dir=tmp_cache_dir)
    def fn(*, x: int) -> int:
        return x

    fn(x=1)  # OK
    with pytest.raises(TypeError, match="keyword arguments only"):
        fn(1)


@pytest.mark.foundation
def test_decorator_key_args_filter(tmp_cache_dir: Path):
    """``key_args`` whitelist excludes non-load-bearing kwargs from the hash."""
    counter = {"n": 0}

    @sood_cache(
        version="v1",
        cache_dir=tmp_cache_dir,
        key_args=["c", "n_modes"],
    )
    def fn(*, c: float, n_modes: int, verbose: bool = False) -> float:
        counter["n"] += 1
        return c * n_modes

    fn(c=1.30, n_modes=10, verbose=False)
    fn(c=1.30, n_modes=10, verbose=True)  # cache hit despite verbose change
    assert counter["n"] == 1


# ═══════════════════════════════════════════════════════════════════
# Corruption recovery
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_corrupt_cache_entry_is_evicted_and_returns_miss(tmp_cache_dir: Path):
    """A corrupt .pkl is treated as a miss and evicted on next read.

    Protects test runs from a partially-written or unparseable cache
    file (e.g. truncated by a prior crashed process).
    """
    cache = SoodResultCache(version="v1", cache_dir=tmp_cache_dir)
    # Manufacture a corrupt entry at the path the cache would write.
    cache.put("solver_a", {"x": 1}, "valid_payload")
    info = cache.info()
    assert len(info) == 1
    path = Path(info[0]["path"])
    path.write_bytes(b"not-a-pickle")
    # Read returns None (miss); the file is removed.
    assert cache.get("solver_a", {"x": 1}) is None
    assert not path.exists()


# ═══════════════════════════════════════════════════════════════════
# Module-level conveniences
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_clear_cache_module_function(tmp_cache_dir: Path):
    """Module-level ``clear_cache`` clears the cache at the given dir."""
    cache = SoodResultCache(version="v1", cache_dir=tmp_cache_dir)
    cache.put("a", {"x": 1}, "alpha")
    cache.put("b", {"x": 2}, "beta")
    n = clear_cache(cache_dir=tmp_cache_dir)
    assert n == 2
    assert cache_info(cache_dir=tmp_cache_dir) == []


# ═══════════════════════════════════════════════════════════════════
# L1 smoke: cache a real solver call
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_cache_with_real_kinf_homogeneous_solver(tmp_cache_dir: Path):
    """Wrap a real solver call (transfer-matrix k_inf on Sood PUa-1-0-IN)
    and verify the cached payload byte-equals the recomputed payload.

    The k_inf solver is fast (microseconds) — the test exists to
    confirm the cache works end-to-end with an actual production
    solver, not for a wall-clock speedup. The same wiring scales to
    slow solvers (sphere F_N at N≥15, Variant α at fine quadrature)
    where the speedup would be 1-2 orders of magnitude.
    """
    from orpheus.derivations.common.eigenvalue import kinf_homogeneous

    counter = {"n": 0}

    @sood_cache(version="v1", cache_dir=tmp_cache_dir)
    def solve(*, case_id: str) -> float:
        counter["n"] += 1
        mix = PUA_1_0_IN.materials[0]
        st = np.asarray(mix.SigT, dtype=float)
        ss = mix.SigS[0].toarray().astype(float)
        nsf = np.asarray(mix.SigP, dtype=float)
        ch = np.asarray(mix.chi, dtype=float)
        return float(kinf_homogeneous(st, ss, nsf, ch))

    k1 = solve(case_id="PUa-1-0-IN")
    k2 = solve(case_id="PUa-1-0-IN")
    assert counter["n"] == 1
    # Byte-equal float — pickle round-trip preserves IEEE-754.
    assert k1 == k2
    # Sanity-check Sood truth (not the gate; a smoke check).
    assert abs(k1 - 2.612903) < 1e-5
