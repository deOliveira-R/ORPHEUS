"""Unified registry of all analytical verification cases.

Imports all derivation modules and collects their VerificationCase
objects into a single dict for easy lookup by tests and documentation.
"""

from __future__ import annotations

from ._types import VerificationCase

# Registry populated lazily on first access
_CASES: dict[str, VerificationCase] | None = None


def _build_registry() -> dict[str, VerificationCase]:
    """Import all derivation modules and collect cases."""
    cases: dict[str, VerificationCase] = {}

    # Phase 1: analytical derivations (own k_inf)
    from . import homogeneous, cp_slab, cp_cylinder, diffusion, sn, moc, mc
    for module in [homogeneous, sn, cp_slab, cp_cylinder, moc, mc, diffusion]:
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
    cases = _ensure_loaded()
    return cases[name]


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
