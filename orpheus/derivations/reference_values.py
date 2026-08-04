"""Unified registry of all verification reference solutions.

Holds two parallel registries during the Phase-0 → Phase-6 migration:

1. ``_CASES`` — the legacy :class:`VerificationCase` registry keyed
   by name, populated from existing ``derive_*`` functions. Every
   currently-green test consumes this.

2. ``_CONTINUOUS`` — the new :class:`ContinuousReferenceSolution`
   registry. Populated by derivations that have been retrofitted to
   the Phase-0 contract. Keys are identical to ``_CASES`` when a
   derivation has been upgraded, so a test can migrate one call site
   at a time.

Retrieval functions:

- :func:`get` — legacy, returns a :class:`VerificationCase`.
- :func:`continuous_get` — returns a :class:`ContinuousReferenceSolution`
  if registered; raises :class:`KeyError` otherwise. Tests that need
  a continuous reference call this; tests that only need a scalar
  ``k_inf`` can stay on :func:`get`.
- :func:`continuous_all_names`, :func:`continuous_all` — enumerate
  the upgraded derivations.

The two registries are additive — upgrading a derivation registers
the new :class:`ContinuousReferenceSolution` **without** removing
the legacy case (the new one can call ``.as_verification_case()`` to
produce a backward-compatible bridge), so the migration is incremental.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from typing import TypeVar

from .common.continuous_reference import ContinuousReferenceSolution
from .common.verification_case import VerificationCase

# Legacy registry populated lazily on first access
_CASES: dict[str, VerificationCase] | None = None

# Phase-0 continuous-reference registry. ``_CONTINUOUS`` holds eagerly-built
# references (cheap producers) plus any lazily-materialised ones (memoised on
# first access). ``_CONTINUOUS_BUILDERS`` holds ``name -> thunk`` for producers
# that opt into the lazy ``continuous_case_builders()`` contract, so an
# O(minutes) reference (e.g. the Peierls adaptive-mpmath solves) is built only
# when that exact name is requested — Issue #212 (fetching ANY reference used to
# pay the full Peierls build cost, which looked like a solver hang).
_CONTINUOUS: dict[str, ContinuousReferenceSolution] | None = None
_CONTINUOUS_BUILDERS: dict[str, Callable[[], ContinuousReferenceSolution]] | None = None


def _build_registry() -> dict[str, VerificationCase]:
    """Import all derivation modules and collect analytical/semi-analytical cases.

    The legacy T3 ``solver_cases()`` side-channel (Richardson-extrapolated
    heterogeneous references computed BY the solvers under test) is fully
    retired: ``sn.py`` deleted its arm in Phase 2.1a (superseded by the
    MMS continuous references), ``moc.py`` never regrew one, and
    ``diffusion.py`` retired the last one at #290 P6 together with the
    MATLAB-port island solver that computed it. Every registry entry is
    now analytical or semi-analytical — no reference is derived from a
    production solver.
    """
    cases: dict[str, VerificationCase] = {}

    from .continuous.analytical import homogeneous
    from .continuous.cases import diffusion, sn, moc, mc
    from .continuous.flat_source_cp import slab as cp_slab
    from .continuous.flat_source_cp import cylinder as cp_cylinder
    from .continuous.flat_source_cp import sphere as cp_sphere
    for module in [homogeneous, sn, cp_slab, cp_cylinder, cp_sphere, moc, mc, diffusion]:
        for case in module.all_cases():
            cases[case.name] = case

    return cases


def _ensure_loaded() -> dict[str, VerificationCase]:
    global _CASES
    if _CASES is None:
        _CASES = _build_registry()
    return _CASES


def get(name: str) -> VerificationCase:
    """Get a verification case by name.

    Raises KeyError if not found.
    """
    return _ensure_loaded()[name]


def all_names() -> list[str]:
    """List all available verification case names."""
    return sorted(_ensure_loaded().keys())


def all_cases() -> list[VerificationCase]:
    """Return all verification cases."""
    return list(_ensure_loaded().values())


def by_geometry(geometry: str) -> list[VerificationCase]:
    """Filter cases by geometry type."""
    return [c for c in _ensure_loaded().values() if c.geometry == geometry]


def by_groups(n_groups: int) -> list[VerificationCase]:
    """Filter cases by number of energy groups."""
    return [c for c in _ensure_loaded().values() if c.n_groups == n_groups]


def by_method(method: str) -> list[VerificationCase]:
    """Filter cases by solver method (homo, cp, sn, moc, mc, dif)."""
    return [c for c in _ensure_loaded().values() if c.method == method]


# ═══════════════════════════════════════════════════════════════════════
# Phase-0 continuous-reference registry
# ═══════════════════════════════════════════════════════════════════════

def _build_continuous_registry() -> tuple[
    dict[str, ContinuousReferenceSolution],
    dict[str, Callable[[], ContinuousReferenceSolution]],
]:
    """Import retrofitted derivation modules and collect their continuous references.

    As each module in :mod:`orpheus.derivations` is upgraded to the
    Phase-0 contract (a ``continuous_cases() -> list[ContinuousReferenceSolution]``
    entry point), it is **auto-discovered** — no manual registration
    required. This replaces the previous manual ``_continuous_modules``
    list that required incremental updates and was a known source of
    "module added but not registered" bugs (e.g., the cylinder / sphere
    hollow F.4 registration this required two follow-up commits).

    The walker examines every ``orpheus.derivations.*`` module, checks
    for a callable ``continuous_cases`` attribute, and ingests its
    results. Modules with no ``continuous_cases`` are silently skipped
    (the contract is opt-in via the presence of the function).

    What the walk costs, and why the three guards below exist
    --------------------------------------------------------

    Auto-discovery buys exactly one thing — it **cannot go stale by
    omission**, which is the bug the manual list kept producing. It pays
    for that with a hazard class the rest of this codebase deliberately
    does NOT carry: the consumer edge is created by a *property* (living
    under this package root), never by a name, so it is invisible to a
    call graph, a text grep, and a constructor search alike. Measured
    2026-08-03: this walk **imports 134 modules to reach 6** that carry a
    hook — 96 % of the imports exist only for their side effects, which
    is how a retired module's ``DeprecationWarning`` surfaced inside
    unrelated SN eigenvalue tests.

    :mod:`orpheus.transport.method` faced the same trade for boundary
    laws and ruled the other way — per-method ``BOUNDARY_OPERATOR_REGISTRY``
    literals, "no central registration step", on the grounds that
    dissolving the indirection *deletes* the hazard class rather than
    gating it. That ruling is right for a set whose membership is a
    semantic decision by an owner. It does not transfer here: this set is
    open, grows module-by-module as derivations are retrofitted, and has
    no single owner to decide admission — which is what keeps the walk.

    So the walk stays and its three silent failure modes are closed
    instead. Each was silent because the walk's *aggregate output* was
    never a declared artifact; the gate in
    ``tests/derivations/test_continuous_registry_lazy.py`` makes it one.
    """
    import importlib
    import pkgutil

    from . import __name__ as pkg_name
    from . import __path__ as pkg_path

    refs: dict[str, ContinuousReferenceSolution] = {}
    builders: dict[str, Callable[[], ContinuousReferenceSolution]] = {}

    # Walk every ``orpheus.derivations.*`` module recursively (including
    # the ``common``, ``discrete`` and ``continuous`` sub-packages added
    # by the three-path reorganisation). Skip private modules (leading
    # underscore) since the contract is "opt-in public".
    for module_info in pkgutil.walk_packages(pkg_path, prefix=f"{pkg_name}."):
        # Skip any path component starting with an underscore — keeps
        # private utility modules and ``__pycache__`` out of the walk.
        if any(part.startswith("_") for part in module_info.name.split(".")):
            continue
        try:
            module = importlib.import_module(module_info.name)
        except ImportError as exc:
            # GUARD 1 — do not swallow a break in our OWN tree.
            #
            # This was a bare ``continue``, justified as tolerating "a docs
            # build with solvers off the path". That tolerance is correct
            # for a MISSING THIRD-PARTY dependency and wrong for a module
            # of ours that fails to import: the reference it publishes
            # simply vanishes from the registry, every test that requests
            # it fails with a confusing "unknown reference", and nothing
            # names the real cause. Discriminate on WHICH module was not
            # found — ``ImportError.name`` is the missing one, not the one
            # being imported.
            missing = getattr(exc, "name", None) or ""
            if missing == pkg_name or missing.startswith(f"{pkg_name}."):
                raise ImportError(
                    f"{module_info.name} is part of {pkg_name} but failed to "
                    f"import ({exc}). A broken module in our own tree is a "
                    f"bug, not an optional dependency: it would silently "
                    f"remove every reference that module publishes from the "
                    f"registry. Fix the import, or move the module under a "
                    f"leading-underscore path to opt it out of discovery."
                ) from exc
            continue
        # Prefer the lazy ``continuous_case_builders()`` contract (Issue #212):
        # record ``name -> thunk`` WITHOUT building, so an expensive producer's
        # eigenvalue solves run only when that exact reference is requested. A
        # module exposing builders is NOT also eagerly built. Modules that only
        # expose ``continuous_cases()`` (cheap producers) are built immediately.
        builders_fn = getattr(module, "continuous_case_builders", None)
        if callable(builders_fn):
            _claim(builders, refs, builders_fn(), module_info.name)
            continue
        cases_fn = getattr(module, "continuous_cases", None)
        if not callable(cases_fn):
            continue
        _claim(refs, builders, {ref.name: ref for ref in cases_fn()},
               module_info.name)

    return refs, builders


#: Value type of whichever registry ``_claim`` is inserting into — an eager
#: :class:`ContinuousReferenceSolution` or a lazy thunk producing one.
_Claimed = TypeVar("_Claimed")


def _claim(
    target: dict[str, _Claimed],
    other: Mapping[str, object],
    incoming: Mapping[str, _Claimed],
    module_name: str,
) -> None:
    """Insert ``incoming`` into ``target``, refusing to overwrite a claimed name.

    GUARD 2 — a reference name has exactly one producer.

    Both halves of the walk used to be plain assignment (``refs[name] = ref``
    / ``builders.update(...)``), which is last-wins. Two modules publishing
    the same name would therefore resolve by ``walk_packages`` traversal
    order, silently, with no consumer anywhere to explain which one won —
    and the losing producer's math would simply never be exercised. That is
    precisely the shape a deprecation shim would have taken had it
    re-exported the discovery hooks alongside its successor.

    ``other`` is the sibling registry: a name must not be claimed as an
    eager case in one module and a lazy builder in another, since
    :func:`get` resolves the two in a fixed order and the loser would be
    dead weight.
    """
    for name, value in incoming.items():
        if name in target or name in other:
            raise ValueError(
                f"{module_name} publishes continuous reference {name!r}, "
                f"which is already claimed. A reference name has exactly one "
                f"producer — resolving this by walk order would silently pick "
                f"one implementation and never exercise the other. Rename one "
                f"of them, or delete the duplicate producer."
            )
        target[name] = value


def _ensure_continuous_loaded() -> dict[str, ContinuousReferenceSolution]:
    global _CONTINUOUS, _CONTINUOUS_BUILDERS
    if _CONTINUOUS is None:
        _CONTINUOUS, _CONTINUOUS_BUILDERS = _build_continuous_registry()
    return _CONTINUOUS


def _materialise_all_continuous() -> dict[str, ContinuousReferenceSolution]:
    """Force every lazy builder into ``_CONTINUOUS`` (Issue #212).

    Audit-only path: pays the full build cost (e.g. the Peierls mpmath
    solves). Used by :func:`continuous_all` / :func:`continuous_by_operator_form`
    which need the materialised objects, not just the names.
    """
    refs = _ensure_continuous_loaded()
    for name, build in list((_CONTINUOUS_BUILDERS or {}).items()):
        if name not in refs:
            refs[name] = build()
    return refs


def continuous_register(ref: ContinuousReferenceSolution) -> None:
    """Register a :class:`ContinuousReferenceSolution` into the registry.

    The preferred registration path is to add the producing module
    to the ``_continuous_modules`` list inside
    :func:`_build_continuous_registry`. This explicit entry point
    exists for tests and one-off derivations that need to inject
    a reference at import time without touching the registry source.
    """
    refs = _ensure_continuous_loaded()
    refs[ref.name] = ref


def continuous_get(name: str) -> ContinuousReferenceSolution:
    """Retrieve a continuous reference solution by name.

    Raises :class:`KeyError` if ``name`` has not been upgraded yet.
    Tests that need a continuous reference should call this directly;
    tests that only need a scalar ``k_inf`` should stay on
    :func:`get` until Phase 2 of the migration.
    """
    refs = _ensure_continuous_loaded()
    if name in refs:
        return refs[name]
    # Materialise a lazily-registered reference on first request, then memoise
    # so subsequent fetches are O(1) (Issue #212).
    if _CONTINUOUS_BUILDERS and name in _CONTINUOUS_BUILDERS:
        refs[name] = _CONTINUOUS_BUILDERS[name]()
        return refs[name]
    return refs[name]  # not registered → standard KeyError with the name


def continuous_all_names() -> list[str]:
    """List every registered continuous reference solution name.

    Cheap by construction: enumerates eager references and lazy-builder
    names WITHOUT triggering any O(minutes) build (Issue #212).
    """
    _ensure_continuous_loaded()
    names = set(_CONTINUOUS or {}) | set(_CONTINUOUS_BUILDERS or {})
    return sorted(names)


def continuous_all() -> list[ContinuousReferenceSolution]:
    """Return every registered continuous reference solution.

    Forces all lazy builders (Issue #212) — the audit path, which pays the
    full build cost. Prefer :func:`continuous_all_names` when only the names
    are needed.
    """
    return list(_materialise_all_continuous().values())


def continuous_by_operator_form(form: str) -> list[ContinuousReferenceSolution]:
    """Filter continuous references by :attr:`ContinuousReferenceSolution.operator_form`.

    Used by the verification audit tool to group references by the
    equation form they commit to, and by tests that want to pull
    every reference valid for their solver's operator.
    """
    return [
        r for r in _materialise_all_continuous().values()
        if r.operator_form == form
    ]
