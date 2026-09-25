"""Generate ``docs/theory/verification/matrix.rst`` from the pytest test registry.

Runs ``python -m tests._harness.audit --json`` under the hood to
populate :data:`tests._harness.registry.TEST_REGISTRY`, then emits a
Sphinx RST page with:

- Overall V&V level distribution (L0/L1/L2/L3/foundation/unmarked counts)
- Per-module level × count grid
- Equation coverage table (label → number of declared tests)
- Orphan equations (`.. math:: :label:` blocks with zero declared tests,
  excluding ``.. vv-status: <label> documented`` labels)
- Documented-only labels (excluded from the orphan gate)
- Phantom verifies-targets (tests naming a ``:label:`` that no longer
  exists anywhere under ``docs/`` — the inverse orphan gate, issue #224)
- Claims held by a withdrawal (labels whose every carrier is a test
  marked ``@pytest.mark.withdrawn``; neither covered nor orphan)
- a pointer to the L0 error catalogue (its own page and generator since #308)
- Unmarked tests listing

The ``foundation`` bucket is orthogonal to the L0..L3 physics ladder
— foundation tests verify software invariants (data-structure
contracts, numerical primitives, factory outputs) rather than physics
equations. They appear in their own column in the module grid and
never contribute to the equation-coverage or orphan-equation tables.
See the ``vv-foundation-tests`` section of the harness architecture
page.

The page is built every time Sphinx rebuilds (the generator is
invoked as a ``pre-build`` step, see ``docs/conf.py`` and the
harness architecture page). It closes ORPHEUS issue #79
("Sphinx L0 verification page generated from test docstrings").

Usage::

    python -m tools.verification.generate_matrix [OUT_RST]

``OUT_RST`` defaults to ``docs/theory/verification/matrix.rst``.
"""

from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

from tests._harness.audit import (
    AUDIT_SNAPSHOT,
    VV_AUDIT_SKIP_MARKER,
    AuditFailed,
    audit_payload,
    write_audit_snapshot,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = REPO_ROOT / "docs" / "theory" / "verification" / "matrix.rst"


def _run_audit() -> dict:
    """Run the harness audit and return its JSON payload.

    On audit failure (e.g. vv-status sentinel violations, collection
    error) the audit's stderr is re-emitted on THIS process's stderr
    before exiting non-zero, so the Sphinx ``builder-inited`` hook's
    ``CalledProcessError.stderr`` carries the actual diagnostic (the
    hook logs it as a build warning — fatal under ``-W``).
    """
    try:
        return audit_payload()
    except AuditFailed as exc:
        sys.stderr.write(exc.stderr)
        raise SystemExit(exc.returncode) from exc


def _render(payload: dict) -> str:
    total = payload["total"]
    by_level = payload["by_level"]
    by_source = payload["by_source"]
    grid = payload["grid"]
    coverage = payload["equation_coverage"]
    orphans = payload["orphan_equations"]
    # ``documented_equations`` was added by the :vv-status: directive
    # work (Phase B.0 of issue #87). Older audit payloads may not
    # include it, so fall back to an empty list for robustness.
    documented = payload.get("documented_equations", [])
    untagged = payload["untagged"]

    lines: list[str] = []

    # Header. The page lives INSIDE the audit's theory tree, so it
    # must opt out of the label/sentinel scan: its label mentions are
    # prose about the census, not declarations. The marker string is
    # imported from the audit (the parser side) so emit and parse
    # cannot drift; the audit reports the exclusion.
    lines.append(f"{VV_AUDIT_SKIP_MARKER}\n\n")
    lines.append("Verification Matrix\n")
    lines.append("===================\n\n")
    lines.append(
        ".. note::\n\n"
        "   Auto-generated from ``tests._harness.registry.TEST_REGISTRY``\n"
        "   by ``tools/verification/generate_matrix.py``. Do not edit by\n"
        "   hand — changes will be overwritten on the next rebuild.\n\n"
    )

    lines.append(f"Total tests collected: **{total}**\n\n")

    # V&V level distribution. ``foundation`` is orthogonal to the
    # L0..L3 ladder and reported alongside it for visibility.
    lines.append("V&V level distribution\n")
    lines.append("----------------------\n\n")
    level_rows = []
    for lvl in ("L0", "L1", "L2", "L3", "foundation", "unmarked"):
        count = by_level.get(lvl, 0)
        pct = f"{100 * count / total:.1f}%" if total else "0.0%"
        level_rows.append([lvl, str(count), pct])
    lines.append(".. csv-table::\n")
    lines.append("   :header: Level, Count, Share\n")
    lines.append("   :widths: 15, 10, 10\n\n")
    for row in level_rows:
        lines.append(f"   {row[0]}, {row[1]}, {row[2]}\n")
    lines.append("\n")

    # Tagging source distribution
    lines.append("Tagging source\n")
    lines.append("--------------\n\n")
    lines.append(
        "How each test acquired its V&V level "
        "(see ``tests/conftest.py`` for the precedence chain).\n\n"
    )
    src_rows = []
    for src in (
        "explicit",
        "class-name",
        "func-name",
        "case",
        "unmarked",
    ):
        count = by_source.get(src, 0)
        src_rows.append([src, str(count)])
    lines.append(".. csv-table::\n")
    lines.append("   :header: Source, Count\n")
    lines.append("   :widths: 20, 10\n\n")
    for row in src_rows:
        lines.append(f"   {row[0]}, {row[1]}\n")
    lines.append("\n")

    # Module × level grid. ``FD`` is the foundation column.
    lines.append("Module × level grid\n")
    lines.append("-------------------\n\n")
    mod_rows = []
    for module in sorted(grid):
        row = grid[module]
        mod_rows.append(
            [
                module,
                str(row.get("L0", 0)),
                str(row.get("L1", 0)),
                str(row.get("L2", 0)),
                str(row.get("L3", 0)),
                str(row.get("foundation", 0)),
                str(row.get("unmarked", 0)),
            ]
        )
    lines.append(".. csv-table::\n")
    lines.append("   :header: Module, L0, L1, L2, L3, FD, ??\n")
    lines.append("   :widths: 40, 6, 6, 6, 6, 6, 6\n\n")
    for row in mod_rows:
        lines.append(f"   {', '.join(row)}\n")
    lines.append("\n")

    # Equation coverage
    lines.append("Equation coverage\n")
    lines.append("-----------------\n\n")
    lines.append(
        "Every Sphinx ``.. math:: :label:`` block declared under "
        "``docs/theory/**/*.rst`` (recursive) and the number of RUNNING "
        "tests carrying ``@pytest.mark.verifies(\"label\")`` that "
        "reference it. A test marked ``@pytest.mark.withdrawn`` is not "
        "counted here; see \"Claims held by a withdrawal\".\n\n"
    )
    if coverage:
        lines.append(".. csv-table::\n")
        lines.append("   :header: Equation label, Tests\n")
        lines.append("   :widths: 50, 10\n\n")
        for eq in sorted(coverage, key=lambda e: (-len(coverage[e]), e)):
            lines.append(f"   ``{eq}``, {len(coverage[eq])}\n")
    else:
        lines.append("*(no equations declared)*\n")
    lines.append("\n")

    # Orphan equations
    lines.append("Orphan equations\n")
    lines.append("----------------\n\n")
    lines.append(
        f"Equations with zero tests carrying "
        f"``@pytest.mark.verifies(\"label\")``, excluding labels "
        f"explicitly marked ``.. vv-status: <label> documented`` and "
        f"labels held by a withdrawal. "
        f"**{len(orphans)}** of the testable equations found on "
        f"theory pages are orphan.\n\n"
    )
    if orphans:
        for eq in sorted(orphans):
            lines.append(f"- ``{eq}``\n")
    else:
        lines.append("*(none — every testable theory equation has at "
                     "least one verifying test)*\n")
    lines.append("\n")

    # Documented-only equations (excluded from the orphan gate)
    lines.append("Documented-only equations\n")
    lines.append("-------------------------\n\n")
    lines.append(
        f"Theory labels marked ``.. vv-status: <label> documented`` in "
        f"their RST source. These are excluded from the orphan-equation "
        f"gate because they are either definitional (no single "
        f"implementing function — e.g. ``boltzmann``), describe a "
        f"module whose Python port does not yet exist (e.g. the "
        f"thermal-hydraulics / fuel-behaviour / reactor-kinetics "
        f"equations), or have a deliberately deferred test paired with "
        f"a tracking issue. **{len(documented)}** labels carry the "
        f"sentinel. See :ref:`vv-status-documented` for the full "
        f"taxonomy.\n\n"
    )
    if documented:
        for eq in sorted(documented):
            lines.append(f"- ``{eq}``\n")
    else:
        lines.append("*(none)*\n")
    lines.append("\n")

    # Scan-exempt files (the ``.. vv-audit: skip-file`` marker)
    skipped = payload.get("skipped_theory_files", [])
    lines.append("Scan-exempt files\n")
    lines.append("-----------------\n\n")
    lines.append(
        f"Files under the theory tree excluded from the label/sentinel "
        f"scan by an explicit ``.. vv-audit: skip-file`` marker — "
        f"syntax-teaching and generated pages whose label mentions are "
        f"not declarations (see the harness architecture page). "
        f"**{len(skipped)}** file(s).\n\n"
    )
    if skipped:
        for f in sorted(skipped):
            lines.append(f"- ``{f}``\n")
    else:
        lines.append("*(none)*\n")
    lines.append("\n")

    # Phantom verifies-targets (the inverse orphan gate, issue #224)
    phantoms = payload.get("phantom_verifies", {})
    lines.append("Phantom verifies targets\n")
    lines.append("------------------------\n\n")
    lines.append(
        f"Labels declared by ``@pytest.mark.verifies(\"label\")`` with "
        f"NO matching ``:label:`` anywhere under ``docs/`` — the "
        f"inverse of the orphan gate (issue #224): a theory-page label "
        f"rename or removal that is not migrated into its tests "
        f"silently drops those tests from the coverage table above. "
        f"**{len(phantoms)}** phantom label(s).\n\n"
    )
    if phantoms:
        lines.append(".. csv-table::\n")
        lines.append("   :header: Phantom label, Tests\n")
        lines.append("   :widths: 50, 10\n\n")
        for eq in sorted(phantoms):
            lines.append(f"   ``{eq}``, {len(phantoms[eq])}\n")
    else:
        lines.append(
            "*(none — every verifies-target resolves to a live "
            "``:label:``)*\n"
        )
    lines.append("\n")

    # Claims held by a withdrawal: labels whose EVERY carrier is a test
    # marked ``@pytest.mark.withdrawn`` (skipped by default), so they are
    # neither covered above nor orphan below.
    held = payload.get("held_by_withdrawal", {})
    lines.append("Claims held by a withdrawal\n")
    lines.append("---------------------------\n\n")
    lines.append(
        f"Labels whose every carrier is a test marked "
        f"``@pytest.mark.withdrawn(reason, issue=N)``: the test consumes a "
        f"withdrawn reference generator and is skipped unless "
        f"``ORPHEUS_RUN_WITHDRAWN`` names its issue, so the label is "
        f"neither covered (the table above counts running carriers only) "
        f"nor orphan (the list below excludes it). Closing the issue "
        f"returns the carriers. **{len(held)}** label(s).\n\n"
    )
    if held:
        lines.append(".. csv-table::\n")
        lines.append("   :header: Equation label, Withdrawn carriers, Issue\n")
        lines.append("   :widths: 50, 10, 10\n\n")
        for eq in sorted(held):
            issues = " ".join(f"#{n}" for n in held[eq]["issues"])
            lines.append(f"   ``{eq}``, {held[eq]['carriers']}, {issues}\n")
    else:
        lines.append("*(none)*\n")
    lines.append("\n")

    # The L0 error catalogue used to be tabulated here. It is a record
    # of DEFECTS that happened in this codebase; this page is a
    # statement about the TEST REGISTRY. Two concepts, so the catalogue
    # got its own generator (2026-08-17, #308) and its own page — the
    # execution is still consolidated in ``docs/conf.py``'s
    # ``_GENERATORS``, only the concepts are separated.
    lines.append("L0 error-catalog coverage\n")
    lines.append("-------------------------\n\n")
    lines.append(
        "Lives in :doc:`error_catalog` — one ``.. error-entry::`` per\n"
        "defect, each a graph node that ``@pytest.mark.catches`` resolves\n"
        "onto. ``nexus errors`` lists them with their catcher counts,\n"
        "uncaught first; the same table is generated into the\n"
        "``vv-principles`` skill index by\n"
        "``tools/verification/generate_error_index.py``.\n\n"
    )

    # Untagged tests summary
    lines.append("Unmarked tests\n")
    lines.append("--------------\n\n")
    untagged_count = len(untagged)
    if untagged_count:
        lines.append(
            f"**{untagged_count} tests** have no V&V level marker.\n"
            "This is a gap — every test in the tree should carry either\n"
            "a physics-ladder marker (``l0``..``l3``) or the orthogonal\n"
            "``foundation`` marker (``@pytest.mark.foundation``) for\n"
            "tests that verify software invariants rather than physics\n"
            "equations. See :ref:`vv-foundation-tests` for the\n"
            "taxonomy.\n\n"
        )
        by_file: Counter[str] = Counter()
        for nodeid in untagged:
            by_file[nodeid.split("::", 1)[0]] += 1
        lines.append(".. csv-table::\n")
        lines.append("   :header: File, Unmarked tests\n")
        lines.append("   :widths: 60, 10\n\n")
        for f, c in by_file.most_common():
            lines.append(f"   ``{f}``, {c}\n")
    else:
        lines.append(
            "*(none — every test carries an L0/L1/L2/L3 or "
            "foundation marker)*\n"
        )
    lines.append("\n")

    return "".join(lines)


def main(argv: list[str] | None = None) -> int:
    argv = argv if argv is not None else sys.argv[1:]
    out_path = Path(argv[0]) if argv else DEFAULT_OUT
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Persist the payload once per build: the error-catalogue index
    # (``generate_error_index``, at ``build-finished``) reads its dormant
    # column from it rather than collecting the suite a second time. The
    # previous build's copy is deleted FIRST, so a failed audit leaves no
    # stale snapshot behind, and the new one carries its tree stamp.
    AUDIT_SNAPSHOT.unlink(missing_ok=True)
    payload = _run_audit()
    write_audit_snapshot(payload)
    rst = _render(payload)
    out_path.write_text(rst, encoding="utf-8")
    print(f"wrote {out_path.relative_to(REPO_ROOT)} "
          f"({payload['total']} tests, {len(payload['equation_coverage'])} "
          f"equations covered)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
