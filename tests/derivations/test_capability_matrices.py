"""Foundation tests for the capability-matrix meta-generator.

Pins the software invariants of
:mod:`tools.verification.generate_capability_matrices`:

1. ``--check`` mode exits 0 when matrices are in sync with the
   ``cases.py:capability_rows()`` registries (regression gate
   for hand-edited matrices vs registry drift).
2. Generation is deterministic — running the meta-generator twice
   produces byte-identical output.
3. At least the two existing packages (``peierls_nystrom`` and
   ``fn_method``) are auto-discovered. Catches accidental breakage
   of the discovery loop.
4. **Every published row is identifiable** — no two rows of one matrix
   share a name, and none share a full identity (#345).
5. **A registry-keyed package's rows resolve, ALL of them** — and a
   package is registry-keyed all-or-nothing, so a *partial* join reds.

What #345 actually changed, and why item 5 is not what it first looks like
=========================================================================

``peierls_nystrom.rst`` used to carry a ``.. note::`` conceding that
nothing cross-checked ``capability_rows()`` against the continuous-reference
registry, "so a row whose metadata drifts from the registered
``ContinuousReferenceSolution`` will not be caught here". Measuring that gap
found the drift was already real: the matrix computed its :math:`r_0` name
tag as ``round(r0 * 100)`` while the registry used
``round(r0/R_out * 100)``, agreeing only because every shipped outer radius
is ``1.0``.

So the fix was **not** to add the cross-check. A name join between two
hand-written enumerations detects drift *after* someone ships it; hoisting
the enumeration to one grid
(:data:`~orpheus.derivations.continuous.peierls_nystrom.cases.SHIPPED_CLASS_A`)
and the name to one rule
(:func:`~orpheus.derivations.continuous.peierls_nystrom.naming.reference_name`)
means the two **cannot** disagree. Prevention over detection —
`coding-elegance` Patterns 2 and 4.

⚠ That makes the row→registry name comparison partly tautological, and per
`coding-standards`' rewire-demotion clause this file says so rather than
keeping an authoritative name for a check that can no longer fail. What item
5 still genuinely tests is the **discovery and registration path**: the grid
is walked by ``continuous_case_builders()``, which the registry walker in
``reference_values`` must find and register. Rename that contract, break the
``pkgutil`` traversal, or ship a package whose rows are *half* registry-keyed,
and this reds.

**Still ungated, deliberately**: per-*field* agreement (does the row's
``n_groups`` match the built reference's?). That needs the reference OBJECT,
and each Peierls reference is an O(minutes) adaptive-``mpmath`` eigenvalue
solve — Issue #212's documented "hang". The note's own proposed mechanism,
``continuous_all()`` filtered by ``operator_form``, forces all of them; it is
not affordable at any tier. The grid removes the drift this would have caught
for the identity fields; the prose fields (accuracy class, status) have no
registry counterpart to disagree with.

These tests are tagged ``@pytest.mark.foundation`` because they
verify a software invariant of the documentation infrastructure
(the meta-generator), not a numerical solver claim — so they do not
carry a ``vv_level`` or ``verifies(...)`` mark.
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import pytest

from orpheus.derivations.reference_values import continuous_all_names
from tools.verification.generate_capability_matrices import (
    MATRIX_GLOB,
    MATRIX_OUTPUT_DIR,
    _discover_packages,
)

REPO_ROOT = Path(__file__).resolve().parents[2]

#: Packages whose ``capability_rows()`` ``name`` IS the continuous-reference
#: registry key, so every row must resolve.
REGISTRY_KEYED = frozenset({"peierls_nystrom"})

#: Packages whose rows are deliberately NOT registry keys, with the reason.
#: `[M]` 2026-08-09: ``fn_method`` registers **nothing** into the continuous
#: registry — its 9 rows enumerate method *families* (one row covers the 12
#: LA-13511 1G ``k_inf`` cases, another the 6 2G cases, …), and the cases
#: themselves live in the method-agnostic
#: :mod:`orpheus.derivations.continuous.sood_registry` keyed by Sood case
#: name. Its row names are therefore human prose by design. That is a
#: different granularity, not drift — and it is declared here rather than
#: left as an untested silence, so a NEW non-joining package forces the same
#: decision explicitly.
NOT_REGISTRY_KEYED = frozenset({"fn_method"})


def _bare(name: object) -> str:
    """Strip the RST literal markup a rendered row name carries."""
    return re.sub(r"``([^`]+)``", r"\1", str(name)).strip()


def _run_generator(*extra_args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [
            sys.executable,
            "-m",
            "tools.verification.generate_capability_matrices",
            *extra_args,
        ],
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )


@pytest.mark.foundation
def test_check_mode_exits_zero_when_in_sync():
    """``--check`` must report exit 0 when matrices match cases.py.

    First force a fresh write so the on-disk matrices match the
    registry, then run ``--check``. Any drift between the two would
    surface as a non-zero exit — that is exactly what this test
    forbids in the green state.
    """
    write_result = _run_generator()
    assert write_result.returncode == 0, (
        f"meta-generator write failed:\n"
        f"stdout={write_result.stdout!r}\nstderr={write_result.stderr!r}"
    )
    check_result = _run_generator("--check")
    assert check_result.returncode == 0, (
        f"--check failed after a clean write — meta-generator is "
        f"non-deterministic or there is on-disk corruption:\n"
        f"stdout={check_result.stdout!r}\nstderr={check_result.stderr!r}"
    )


@pytest.mark.foundation
def test_generation_is_deterministic():
    """Running the meta-generator twice must produce identical output.

    Reads each ``_<pkg>_capability_matrix.inc.rst`` after the first
    write, runs the generator again, and asserts byte equality.
    Catches non-deterministic dict iteration / set ordering bugs.
    """
    first_write = _run_generator()
    assert first_write.returncode == 0, first_write.stderr

    first_contents = {
        path.name: path.read_text(encoding="utf-8")
        for path in sorted(MATRIX_OUTPUT_DIR.glob(MATRIX_GLOB))
    }
    assert first_contents, "no capability matrices found after first write"

    second_write = _run_generator()
    assert second_write.returncode == 0, second_write.stderr

    second_contents = {
        path.name: path.read_text(encoding="utf-8")
        for path in sorted(MATRIX_OUTPUT_DIR.glob(MATRIX_GLOB))
    }

    assert first_contents.keys() == second_contents.keys(), (
        f"matrix file set changed between runs:\n"
        f"  first:  {sorted(first_contents)}\n"
        f"  second: {sorted(second_contents)}"
    )
    for name in first_contents:
        assert first_contents[name] == second_contents[name], (
            f"matrix {name} differs between two consecutive runs — "
            f"meta-generator is non-deterministic"
        )


@pytest.mark.foundation
def test_at_least_peierls_nystrom_and_fn_method_discovered():
    """The two existing packages must be auto-discovered.

    Asserts that running the meta-generator emits matrices for at
    least ``peierls_nystrom`` and ``fn_method``. This is the
    minimum-discovery floor — any future regression of the
    discovery loop (broken ``pkgutil.iter_modules`` traversal,
    accidental rename of ``cases.capability_rows``, etc.) trips
    this test.
    """
    write_result = _run_generator()
    assert write_result.returncode == 0, write_result.stderr

    found = {
        path.name for path in MATRIX_OUTPUT_DIR.glob(MATRIX_GLOB)
    }
    expected = {
        "_peierls_nystrom_capability_matrix.inc.rst",
        "_fn_method_capability_matrix.inc.rst",
    }
    missing = expected - found
    assert not missing, (
        f"meta-generator failed to discover expected packages: "
        f"missing={sorted(missing)}; found={sorted(found)}"
    )


# ---------------------------------------------------------------------
# #345 — the published rows and the registry
# ---------------------------------------------------------------------


@pytest.mark.foundation
def test_every_discovered_package_declares_whether_its_rows_are_registry_keys():
    """A new package must not silently join the "not checked" bucket.

    The two buckets are a *decision* about what a package's ``name``
    column means. Leaving a third package undeclared is how the original
    #345 gap arose — nothing said whether the join was expected to hold,
    so nobody noticed it did not exist.
    """
    discovered = {pkg for pkg, _rows in _discover_packages()}
    declared = REGISTRY_KEYED | NOT_REGISTRY_KEYED
    assert discovered == declared, (
        f"undeclared package(s) {sorted(discovered - declared)}; "
        f"stale declaration(s) {sorted(declared - discovered)}.\n"
        f"Add each new package to REGISTRY_KEYED (its row ``name`` is a "
        f"continuous-registry key, so every row must resolve) or to "
        f"NOT_REGISTRY_KEYED (its rows are a different granularity — say "
        f"which, and why, in the constant's comment)."
    )


@pytest.mark.foundation
@pytest.mark.parametrize("pkg_name", sorted(REGISTRY_KEYED))
def test_every_registry_keyed_row_names_a_registered_solution(pkg_name):
    """Every published row of a registry-keyed package resolves.

    The doc's promised cross-check, in its affordable form. It is by
    NAME, not by field: materialising a Peierls reference is an
    O(minutes) mpmath solve (#212), so ``continuous_all_names()`` — cheap
    by construction — is the only registry view a default-tier gate can
    consult.

    ⚠ Since #345 the row name and the registry key are both produced by
    :func:`~orpheus.derivations.continuous.peierls_nystrom.naming.reference_name`
    from one grid, so this cannot fail on a *spelling* divergence any
    more — that class is now unspellable. It still fails if the package
    stops being registered at all (a broken ``pkgutil`` walk, a renamed
    ``continuous_case_builders`` contract, a builder that raises during
    registration), which is a live failure mode with no other catcher.
    """
    rows = dict(_discover_packages())[pkg_name]
    registered = set(continuous_all_names())
    unresolved = sorted(
        _bare(row["name"]) for row in rows if _bare(row["name"]) not in registered
    )
    assert not unresolved, (
        f"{pkg_name}.capability_rows() publishes {len(unresolved)} row(s) "
        f"naming no registered continuous reference: {unresolved}\n"
        f"The matrix is rendered into the theory corpus, so a row naming a "
        f"reference a reader cannot fetch is a published falsehood."
    )


@pytest.mark.foundation
@pytest.mark.parametrize("pkg_name", sorted(NOT_REGISTRY_KEYED))
def test_a_non_registry_keyed_package_joins_ZERO_rows(pkg_name):
    """All-or-nothing: a PARTIAL join is the drift signal.

    A package declared non-registry-keyed whose rows start resolving is
    mid-migration — someone registered part of the family. That state is
    exactly where a reader's "these names are registry keys" inference
    silently becomes half-true, so it reds here and forces the package
    into :data:`REGISTRY_KEYED` (finishing the job) or the new rows back
    out.
    """
    rows = dict(_discover_packages())[pkg_name]
    registered = set(continuous_all_names())
    resolved = sorted(
        _bare(row["name"]) for row in rows if _bare(row["name"]) in registered
    )
    assert not resolved, (
        f"{pkg_name} is declared NOT_REGISTRY_KEYED but {len(resolved)} of "
        f"its {len(rows)} rows now resolve against the continuous registry: "
        f"{resolved}\nFinish the migration and move the package to "
        f"REGISTRY_KEYED, or drop the newly-registry-named rows."
    )


@pytest.mark.foundation
@pytest.mark.parametrize("pkg_name", sorted(REGISTRY_KEYED | NOT_REGISTRY_KEYED))
def test_published_row_names_are_unique_within_a_matrix(pkg_name):
    """The ``name`` column is what a reader scans; it must identify a row.

    `[M]` 2026-08-09: ``fn_method`` shipped two rows both labelled
    ``k_inf — mG general (Sood Eq 76)``, distinguishable only by the
    ``n_g`` column — off-pattern against every sibling row, all of which
    name their order (``1G isotropic``, ``2G no-upscatter``). Renamed to
    ``3G`` / ``6G`` at #345.
    """
    rows = dict(_discover_packages())[pkg_name]
    seen: dict[str, int] = {}
    for row in rows:
        seen[_bare(row["name"])] = seen.get(_bare(row["name"]), 0) + 1
    duplicates = sorted(name for name, n in seen.items() if n > 1)
    assert not duplicates, (
        f"{pkg_name}.capability_rows() publishes duplicate row name(s): "
        f"{duplicates}. Two rows under one label are indistinguishable to "
        f"a reader scanning the matrix; add the discriminator (group count, "
        f"geometry, r_0/R) to the name itself."
    )


@pytest.mark.foundation
@pytest.mark.parametrize("pkg_name", sorted(REGISTRY_KEYED | NOT_REGISTRY_KEYED))
def test_no_two_published_rows_are_the_same_row(pkg_name):
    """Full-identity uniqueness — the copy-paste catcher.

    Weaker than name-uniqueness above and kept separately on purpose:
    this one cannot false-red on a legitimate same-family pair that
    differs in some column, so it stays valid even if a future matrix
    genuinely wants a repeated label.
    """
    rows = dict(_discover_packages())[pkg_name]
    identities = [tuple(sorted((k, str(v)) for k, v in row.items())) for row in rows]
    assert len(set(identities)) == len(identities), (
        f"{pkg_name}.capability_rows() contains rows that are identical in "
        f"every field — a copy-paste duplicate, not a distinct capability."
    )
