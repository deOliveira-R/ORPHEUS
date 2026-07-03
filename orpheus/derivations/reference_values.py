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

from collections.abc import Callable

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
    # underscore) since the contract is "opt-in public", and skip any
    # module whose import fails (e.g. on a docs build with solvers off
    # the path).
    for module_info in pkgutil.walk_packages(pkg_path, prefix=f"{pkg_name}."):
        # Skip any path component starting with an underscore — keeps
        # private utility modules and ``__pycache__`` out of the walk.
        if any(part.startswith("_") for part in module_info.name.split(".")):
            continue
        try:
            module = importlib.import_module(module_info.name)
        except ImportError:
            continue
        # Prefer the lazy ``continuous_case_builders()`` contract (Issue #212):
        # record ``name -> thunk`` WITHOUT building, so an expensive producer's
        # eigenvalue solves run only when that exact reference is requested. A
        # module exposing builders is NOT also eagerly built. Modules that only
        # expose ``continuous_cases()`` (cheap producers) are built immediately.
        builders_fn = getattr(module, "continuous_case_builders", None)
        if callable(builders_fn):
            builders.update(builders_fn())
            continue
        cases_fn = getattr(module, "continuous_cases", None)
        if not callable(cases_fn):
            continue
        for ref in cases_fn():
            refs[ref.name] = ref

    return refs, builders


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
