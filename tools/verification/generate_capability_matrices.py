"""Meta-generator for capability matrices via auto-discovery.

Discovers ``cases.py:capability_rows()`` across every package in
``orpheus.derivations.continuous`` and emits one
``docs/theory/_<package_name>_capability_matrix.inc.rst`` per
discovered package. Replaces per-method generators
(``generate_peierls_nystrom_matrix.py``,
``generate_fn_method_matrix.py``) with a single tool that scales
linearly with new methods at zero per-method generator cost.

Schema contract
---------------

Each row dict returned by ``capability_rows()`` must carry the
**required** keys:

- ``name`` — registry / human-readable case name (rendered as the
  first column).
- ``geometry`` — ``"slab" | "cylinder-1d" | "sphere-1d" | "infinite"``.
- ``n_groups`` — energy-group count (``int``).
- ``n_regions`` — spatial-region count (``int``; ``0`` for infinite
  medium).
- ``bc`` — boundary-condition descriptor (RST-ready string).
- ``status`` — production-readiness label (``"shipped" | "stub" | ...``).

Optional keys are auto-detected: a column appears in the rendered
table iff at least one row carries the key. The recognised optional
keys are:

- ``r0_over_R`` — :math:`r_0/R` for hollow curvilinear (rendered as
  ``"—"`` when ``None`` or absent).
- ``closure`` — closure scheme label (matrix-Galerkin specific).
- ``accuracy`` — accuracy-class string.
- ``scattering_order`` — Legendre order (``int``).
- ``multiplying`` — ``True`` for ``c > 1`` cases (``bool``).
- ``orbit_space_class`` — orbit-space M/G class label. peierls_nystrom
  uses ``"A"`` (two-endpoint) / ``"B"`` (one-endpoint); other methods
  may use the values returned by
  :attr:`~orpheus.derivations.continuous.peierls_nystrom.geometry.CurvilinearGeometry.orbit_space_class`
  (``"two-endpoint"`` / ``"one-endpoint"``).
- ``extra`` — free-form ``dict`` of additional metadata; ignored by
  the renderer (carried for downstream consumers).

Discovery rules
---------------

For each ``pkg`` in ``pkgutil.iter_modules(orpheus.derivations.continuous)``:

1. Try ``importlib.import_module(f"{pkg.name}.cases")``.
2. If the import fails (``ImportError``), skip silently.
3. If the module imports but has no ``capability_rows`` attribute,
   skip silently.
4. Otherwise call ``capability_rows()`` and render to
   ``docs/theory/_{pkg.name}_capability_matrix.inc.rst``.

At the start of each run, all existing
``docs/theory/_*_capability_matrix.inc.rst`` files are wiped — stale
matrices for now-removed packages are removed automatically.

Modes
-----

``--check`` : compare what *would* be written against what is
on disk; exit 1 if any matrix differs (used by Sphinx CI to ensure
matrices are committed). Otherwise exit 0.

Without ``--check``, matrices are written to disk in place.

Usage::

    python -m tools.verification.generate_capability_matrices [--check]
"""

from __future__ import annotations

import importlib
import pkgutil
import sys
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]
DOCS_THEORY_DIR = REPO_ROOT / "docs" / "theory"
MATRIX_GLOB = "_*_capability_matrix.inc.rst"


# Required keys every row must carry. Missing keys are a hard error
# (raised before any output is written).
_REQUIRED_KEYS: tuple[str, ...] = (
    "name",
    "geometry",
    "n_groups",
    "n_regions",
    "bc",
    "status",
)

# Optional keys auto-detected per row presence; the column appears in
# the rendered table iff at least one row carries the key.
_OPTIONAL_KEYS: tuple[str, ...] = (
    "r0_over_R",
    "closure",
    "accuracy",
    "scattering_order",
    "multiplying",
    "orbit_space_class",
)

# Display order (left-to-right) for the rendered columns. Required
# keys come first, then optional keys in this priority order.
_COLUMN_ORDER: tuple[str, ...] = (
    "name",
    "geometry",
    "n_groups",
    "n_regions",
    "r0_over_R",
    "orbit_space_class",
    "bc",
    "closure",
    "scattering_order",
    "multiplying",
    "status",
    "accuracy",
)

# Column header labels (mapping field name → RST header text).
_HEADERS: dict[str, str] = {
    "name": "Reference name",
    "geometry": "Geometry",
    "n_groups": "n_g",
    "n_regions": "n_reg",
    "r0_over_R": r":math:`r_0 / R`",
    "orbit_space_class": "Orbit-space M/G",
    "bc": "BC",
    "closure": "Closure",
    "scattering_order": "Scattering order",
    "multiplying": "Multiplying",
    "status": "Status",
    "accuracy": "Accuracy class",
}

# Per-column ``:widths:`` (rough proportional weights).
_WIDTHS: dict[str, int] = {
    "name": 30,
    "geometry": 10,
    "n_groups": 6,
    "n_regions": 6,
    "r0_over_R": 8,
    "orbit_space_class": 8,
    "bc": 14,
    "closure": 14,
    "scattering_order": 8,
    "multiplying": 8,
    "status": 24,
    "accuracy": 18,
}

_GEOM_LABELS: dict[str, str] = {
    "slab": "slab",
    "cylinder-1d": "cylinder-1d",
    "sphere-1d": "sphere-1d",
    "infinite": "infinite",
}


def _validate_rows(pkg_name: str, rows: list[dict[str, Any]]) -> None:
    """Raise ``ValueError`` if any row violates the schema contract."""
    if not isinstance(rows, list):
        raise ValueError(
            f"{pkg_name}.capability_rows() must return list[dict]; "
            f"got {type(rows).__name__}"
        )
    for i, row in enumerate(rows):
        if not isinstance(row, dict):
            raise ValueError(
                f"{pkg_name}.capability_rows()[{i}] must be a dict; "
                f"got {type(row).__name__}"
            )
        missing = [k for k in _REQUIRED_KEYS if k not in row]
        if missing:
            raise ValueError(
                f"{pkg_name}.capability_rows()[{i}] missing required "
                f"keys {missing!r}; row={row!r}"
            )


def _format_cell(field: str, value: Any) -> str:
    """Render a single cell value, with field-aware formatting."""
    if field == "geometry":
        return _GEOM_LABELS.get(str(value), str(value))
    if field == "name":
        # Names are pre-formatted: peierls uses ``backtick`` quoting,
        # fn_method uses plain prose. Pass through verbatim.
        return str(value)
    if field == "r0_over_R":
        if value is None:
            return "—"
        return f"{float(value):.1f}"
    if field == "multiplying":
        return "yes" if bool(value) else "no"
    if field is None or value is None:
        return "—"
    return str(value)


def _list_table(
    title: str,
    columns: list[str],
    rows: list[dict[str, Any]],
) -> str:
    """Render a ``.. list-table::`` block from row dicts."""
    headers = [_HEADERS[c] for c in columns]
    widths = [_WIDTHS[c] for c in columns]
    lines: list[str] = []
    lines.append(f".. list-table:: {title}")
    lines.append("   :header-rows: 1")
    lines.append(f"   :widths: {' '.join(str(w) for w in widths)}")
    lines.append("")
    lines.append(f"   * - {headers[0]}")
    for h in headers[1:]:
        lines.append(f"     - {h}")
    for row in rows:
        cells = [_format_cell(c, row.get(c)) for c in columns]
        lines.append(f"   * - {cells[0]}")
        for cell in cells[1:]:
            lines.append(f"     - {cell}")
    lines.append("")
    return "\n".join(lines)


def _select_columns(rows: list[dict[str, Any]]) -> list[str]:
    """Pick which columns appear in the rendered table.

    Required columns always appear. Optional columns appear iff at
    least one row carries the key. Order follows ``_COLUMN_ORDER``.
    """
    present_optional: set[str] = set()
    for row in rows:
        for key in _OPTIONAL_KEYS:
            if key in row and row[key] is not None:
                present_optional.add(key)
    selected = [
        c for c in _COLUMN_ORDER
        if c in _REQUIRED_KEYS or c in present_optional
    ]
    return selected


def _render(pkg_name: str, rows: list[dict[str, Any]]) -> str:
    """Render the full include file for one package."""
    columns = _select_columns(rows)
    title = f"Shipped {pkg_name} continuous references"

    lines: list[str] = []
    lines.append(
        "..\n"
        "   AUTO-GENERATED by "
        "``tools/verification/generate_capability_matrices.py``.\n"
        "   Do not edit by hand — changes will be overwritten on the\n"
        "   next Sphinx build. The single source of truth is\n"
        f"   ``orpheus.derivations.continuous.{pkg_name}.cases.capability_rows()``.\n"
    )
    lines.append("")
    lines.append(_list_table(title=title, columns=columns, rows=rows))
    return "\n".join(lines) + "\n"


def _discover_packages() -> list[tuple[str, list[dict[str, Any]]]]:
    """Walk ``orpheus.derivations.continuous`` for ``cases.capability_rows``.

    Returns a sorted list of ``(pkg_name, rows)`` tuples for every
    package that exposes the metadata function. Packages without
    ``cases.py`` or without ``capability_rows`` are skipped silently.

    Sorted by ``pkg_name`` for deterministic output ordering.
    """
    import orpheus.derivations.continuous as continuous_pkg

    discovered: list[tuple[str, list[dict[str, Any]]]] = []
    for pkg in pkgutil.iter_modules(continuous_pkg.__path__):
        if not pkg.ispkg:
            continue
        full_name = f"orpheus.derivations.continuous.{pkg.name}.cases"
        try:
            cases_mod = importlib.import_module(full_name)
        except ImportError:
            continue
        fn = getattr(cases_mod, "capability_rows", None)
        if fn is None:
            continue
        rows = fn()
        _validate_rows(pkg.name, rows)
        discovered.append((pkg.name, rows))

    discovered.sort(key=lambda pair: pair[0])
    return discovered


def _wipe_existing_matrices() -> None:
    """Remove every ``_*_capability_matrix.inc.rst`` under docs/theory/."""
    for path in sorted(DOCS_THEORY_DIR.glob(MATRIX_GLOB)):
        path.unlink()


def main(argv: list[str] | None = None) -> int:
    argv = argv if argv is not None else sys.argv[1:]
    check_mode = "--check" in argv

    discovered = _discover_packages()
    if not discovered:
        # No packages — nothing to do, but not an error per se.
        print("no capability_rows() found in any package — nothing to write")
        return 0

    # Build the (path, content) pairs first; only touch disk after all
    # rendering succeeds, so a half-written state cannot leak.
    rendered: list[tuple[Path, str]] = []
    for pkg_name, rows in discovered:
        out_path = (
            DOCS_THEORY_DIR
            / f"_{pkg_name}_capability_matrix.inc.rst"
        )
        rst = _render(pkg_name, rows)
        rendered.append((out_path, rst))

    if check_mode:
        # Compare every (path, content) against on-disk; report and
        # exit 1 if any differ. Also flag stale on-disk matrices that
        # do not correspond to any discovered package.
        rendered_paths = {p for p, _ in rendered}
        existing_paths = set(DOCS_THEORY_DIR.glob(MATRIX_GLOB))
        stale = sorted(existing_paths - rendered_paths)
        differences: list[str] = []
        for stale_path in stale:
            differences.append(
                f"stale (no matching package): "
                f"{stale_path.relative_to(REPO_ROOT)}"
            )
        for out_path, expected in rendered:
            if not out_path.exists():
                differences.append(
                    f"missing: {out_path.relative_to(REPO_ROOT)}"
                )
                continue
            actual = out_path.read_text(encoding="utf-8")
            if actual != expected:
                differences.append(
                    f"out-of-date: {out_path.relative_to(REPO_ROOT)}"
                )
        if differences:
            print("capability matrices out of sync with cases.py:")
            for d in differences:
                print(f"  - {d}")
            print(
                "run `python -m tools.verification.generate_capability_matrices` "
                "to regenerate"
            )
            return 1
        print(
            f"all {len(rendered)} capability matrices are up to date "
            f"({', '.join(name for name, _ in discovered)})"
        )
        return 0

    # Write mode: wipe existing then write rendered.
    _wipe_existing_matrices()
    for out_path, rst in rendered:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(rst, encoding="utf-8")
        print(
            f"wrote {out_path.relative_to(REPO_ROOT)} "
            f"({len(_rows_for_path(rendered, out_path))} rows)"
        )
    return 0


def _rows_for_path(
    rendered: list[tuple[Path, str]],
    out_path: Path,
) -> list[Any]:
    """Tiny helper for the print line — counts rows by re-reading
    the rendered RST. Cheaper than threading the count back."""
    # Approximate row count: count occurrences of "   * - " minus the
    # header row. Good enough for the stdout summary.
    for path, rst in rendered:
        if path == out_path:
            count = rst.count("\n   * - ")
            return list(range(max(0, count - 1)))
    return []


if __name__ == "__main__":
    sys.exit(main())
