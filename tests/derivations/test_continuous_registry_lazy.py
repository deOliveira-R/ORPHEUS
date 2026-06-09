"""Regression gates for the lazy continuous-reference registry (Issue #212).

Before the lazy-by-name registry, :func:`continuous_get` built the **entire**
continuous-reference registry on first access — including 13 Peierls hollow
cyl/sph/slab references, each an adaptive-``mpmath`` eigenvalue solve
(O(minutes) in aggregate). Fetching ANY reference (e.g. an unrelated SN MMS
case) therefore appeared to *hang*. The fix records Peierls references as
``name -> thunk`` via the opt-in ``continuous_case_builders()`` contract and
materialises only the requested reference (memoised).

These gates pin:

* **the fast path** — an SN fetch registers the Peierls names lazily but does
  NOT build them (the #212 win);
* **cheap enumeration** — ``continuous_all_names()`` lists every Peierls name
  without building;
* **the name contract** — the cheap builder keys match the names the inner
  builders stamp onto the BUILT references (drift guard, slow);
* **lazy build correctness** — fetching a Peierls name builds exactly that ref.

The two ``slow`` gates materialise the Peierls references and so pay the
O(minutes) cost; the rest are sub-second.
"""
from __future__ import annotations

import pytest

from orpheus.derivations import reference_values as rv
from orpheus.derivations.continuous.peierls_nystrom import cases as peierls_cases


@pytest.fixture
def fresh_continuous_registry():
    """Snapshot/restore the module-global registry so a test can force a clean
    (re)build without leaking state into the rest of the session."""
    saved = (rv._CONTINUOUS, rv._CONTINUOUS_BUILDERS)
    rv._CONTINUOUS = None
    rv._CONTINUOUS_BUILDERS = None
    try:
        yield
    finally:
        rv._CONTINUOUS, rv._CONTINUOUS_BUILDERS = saved


def _peierls_names() -> set[str]:
    """The shipped Class-A Peierls names, read cheaply from the builder keys."""
    return set(peierls_cases.continuous_case_builders())


@pytest.mark.foundation
def test_sn_fetch_does_not_build_peierls(fresh_continuous_registry):
    """The #212 regression: fetching an SN reference registers the Peierls
    names lazily but does NOT pay their O(minutes) build cost."""
    ref = rv.continuous_get("sn_slab_1eg_2rg_S8")
    assert ref.name == "sn_slab_1eg_2rg_S8"

    peierls = _peierls_names()
    assert peierls, "expected the Peierls module to register lazy builders"

    # Peierls names are registered as lazy builders ...
    assert peierls <= set(rv._CONTINUOUS_BUILDERS or {})
    # ... but NONE were materialised by the SN fetch (the win).
    built = set(rv._CONTINUOUS or {})
    assert not (built & peierls), (
        "fetching an SN reference must not build any Peierls reference (#212)"
    )


@pytest.mark.foundation
def test_continuous_all_names_enumerates_peierls_cheaply(fresh_continuous_registry):
    """All Peierls names enumerate without building (cheap audit path)."""
    names = set(rv.continuous_all_names())
    peierls = _peierls_names()
    assert peierls <= names
    # Enumeration must not have materialised them.
    assert not (set(rv._CONTINUOUS or {}) & peierls)


@pytest.mark.foundation
def test_builder_keyset_is_the_shipped_class_a_inventory():
    """The lazy builder exposes exactly the shipped Class-A inventory:
    slab 2G-2rg + hollow cyl/sph at r0/R ∈ {0.10, 0.20, 0.30}, 1G & 2G."""
    expected = {
        "peierls_slab_2eg_2rg",
        *(
            f"peierls_{shape}_hollow_{ng}eg_1rg_r0_{tag}"
            for shape in ("cyl1D", "sph1D")
            for ng in (1, 2)
            for tag in ("10", "20", "30")
        ),
    }
    assert set(peierls_cases.continuous_case_builders()) == expected


@pytest.mark.slow
@pytest.mark.foundation
def test_builder_keys_match_built_names():
    """Drift guard: the cheap builder keys equal the names the inner builders
    stamp onto the BUILT references. This is the contract that lets
    :func:`continuous_get` trust the lazy keys without ever building. Slow —
    materialises all Peierls references (the O(minutes) cost)."""
    keys = set(peierls_cases.continuous_case_builders())
    built = {c.name for c in peierls_cases.continuous_cases()}
    assert keys == built, (
        f"builder keys drifted from built names; "
        f"missing={sorted(built - keys)}, extra={sorted(keys - built)}"
    )


@pytest.mark.slow
@pytest.mark.foundation
def test_lazy_peierls_fetch_builds_requested_ref(fresh_continuous_registry):
    """Fetching a Peierls name materialises exactly that reference (the
    returned object's name matches the request), and memoises it."""
    ref = rv.continuous_get("peierls_slab_2eg_2rg")
    assert ref.name == "peierls_slab_2eg_2rg"
    assert rv._CONTINUOUS["peierls_slab_2eg_2rg"] is ref
