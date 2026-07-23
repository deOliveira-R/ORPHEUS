"""Unit tests for the V&V audit CLI scanner.

Covers the ``:label:`` / ``.. vv-status: <label> documented`` parsing
that feeds the orphan-equation gate. See ``tests/_harness/audit.py``
for the full behaviour contract; this file fixes the invariants:

* every ``:label:`` in a theory RST becomes a labelled equation,
* a ``.. vv-status: <label> documented`` comment excludes its
  label from the orphan set even when no test verifies it,
* the sentinel schema is FAIL-LOUD (2026-07 single-status ruling,
  task #10): ``documented`` is the only status — an unknown status
  word, a sentinel whose label is absent from the SAME file (dead or
  misplaced), or a malformed line each produce a violation entry that
  hard-fails the audit, never a silent drop. ``tested``/``verified``
  are derived from ``@pytest.mark.verifies``, not hand-asserted.

The tests write fixture RST files into a tmp_path so they are
hermetic and independent of the real ``docs/theory`` tree.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from tests._harness.audit import _scan_theory_equations

pytestmark = pytest.mark.l0


def _write_rst(path: Path, name: str, body: str) -> None:
    (path / f"{name}.rst").write_text(body, encoding="utf-8")


def test_scanner_collects_labels(tmp_path: Path) -> None:
    _write_rst(tmp_path, "sample", """
.. math::
   :label: foo

   \\int f(x)\\,dx = 1

.. math::
   :label: bar

   a + b = c
""")
    labels, documented, violations = _scan_theory_equations(tmp_path)
    assert labels == {"foo", "bar"}
    assert documented == set()
    assert violations == []


def test_scanner_marks_documented(tmp_path: Path) -> None:
    """A same-file ``.. vv-status: <label> documented`` excludes the label.

    Positive leg of the sentinel contract: a correct sentinel MUST
    parse cleanly (zero violations) and populate the documented set.
    """
    _write_rst(tmp_path, "sample", """
.. math::
   :label: boltzmann

   \\partial_t \\psi + \\Omega \\cdot \\nabla \\psi = S

.. vv-status: boltzmann documented

.. math::
   :label: testable

   x = 1
""")
    labels, documented, violations = _scan_theory_equations(tmp_path)
    assert labels == {"boltzmann", "testable"}
    assert documented == {"boltzmann"}
    assert violations == []


def test_scanner_flags_dead_sentinel(tmp_path: Path) -> None:
    """A sentinel naming a label that exists NOWHERE is a violation.

    The pre-2026-07 behaviour silently dropped it (fail-closed for the
    orphan gate but invisible to the author — three dead sentinels
    accumulated in the tree that way). The single-status ruling makes
    it a hard audit error instead.
    """
    _write_rst(tmp_path, "typo", """
.. math::
   :label: real_label

   x = 1

.. vv-status: rael_label documented
""")
    labels, documented, violations = _scan_theory_equations(tmp_path)
    assert labels == {"real_label"}
    assert documented == set()
    assert len(violations) == 1
    assert "rael_label" in violations[0]
    assert "dead sentinel" in violations[0]


def test_scanner_flags_cross_file_sentinel(tmp_path: Path) -> None:
    """A sentinel in a different file from its ``:label:`` is a violation.

    The same-file rule (a sentinel is point-of-use metadata) was
    documented but unenforced; the violation message distinguishes the
    misplaced case from the dead case so the fix is obvious.
    """
    _write_rst(tmp_path, "owner", """
.. math::
   :label: homed_label

   x = 1
""")
    _write_rst(tmp_path, "stray", """
.. vv-status: homed_label documented
""")
    labels, documented, violations = _scan_theory_equations(tmp_path)
    assert labels == {"homed_label"}
    assert documented == set()
    assert len(violations) == 1
    assert "another file" in violations[0]


def test_scanner_flags_unknown_status(tmp_path: Path) -> None:
    """Any status other than ``documented`` is a violation, not a no-op.

    ``tested`` / ``verified`` are DERIVED from ``@pytest.mark.verifies``
    — a hand-written coverage claim would be a second source of truth
    that can lie. 24 such inert sentinels existed before the ruling;
    the parser now rejects the whole class.
    """
    _write_rst(tmp_path, "misc", """
.. math::
   :label: alpha

   x = 1

.. vv-status: alpha verified
""")
    labels, documented, violations = _scan_theory_equations(tmp_path)
    assert labels == {"alpha"}
    assert documented == set()
    assert len(violations) == 1
    assert "'verified'" in violations[0]
    assert "verifies" in violations[0]


def test_scanner_flags_malformed_sentinel(tmp_path: Path) -> None:
    """A vv-status line that is not ``<label> <status>`` is a violation."""
    _write_rst(tmp_path, "broken", """
.. math::
   :label: beta

   x = 1

.. vv-status: beta documented (see issue #42)
""")
    labels, documented, violations = _scan_theory_equations(tmp_path)
    assert labels == {"beta"}
    assert documented == set()
    assert len(violations) == 1
    assert "malformed" in violations[0]


def test_audit_main_exits_2_on_sentinel_violation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """``main()`` returns 2 on a violation, BEFORE collection runs.

    The whole fail-loud chain (audit exit 2 → generator SystemExit →
    conf.py hook warning → ``-W`` fatal) rests on this exit code; a
    refactor that moved the check after collection or softened the
    return would pass every scanner test and silently defeat the
    Sphinx gate, so the propagation itself is pinned here.
    """
    _write_rst(tmp_path, "bad", """
.. math::
   :label: gamma

   x = 1

.. vv-status: gamma tested
""")
    from tests._harness import audit as audit_mod

    collection_calls: list[bool] = []
    monkeypatch.setattr(
        audit_mod, "_run_collection", lambda: collection_calls.append(True) or 0
    )
    rc = audit_mod.main(["--theory-dir", str(tmp_path), "--json"])
    assert rc == 2
    assert collection_calls == [], (
        "sentinel violations must short-circuit BEFORE pytest collection"
    )
    err = capsys.readouterr().err
    assert "sentinel violations" in err
    assert "gamma" in err


def test_scanner_returns_empty_for_missing_dir(tmp_path: Path) -> None:
    labels, documented, violations = _scan_theory_equations(tmp_path / "nope")
    assert labels == set()
    assert documented == set()
    assert violations == []


def test_scanner_handles_subdirectories(tmp_path: Path) -> None:
    """``rglob`` picks up labels in nested RST files."""
    sub = tmp_path / "nested"
    sub.mkdir()
    _write_rst(tmp_path, "top", """
.. math::
   :label: top_label

   1 = 1
""")
    _write_rst(sub, "deep", """
.. math::
   :label: deep_label

   2 = 2

.. vv-status: deep_label documented
""")
    labels, documented, violations = _scan_theory_equations(tmp_path)
    assert labels == {"top_label", "deep_label"}
    assert documented == {"deep_label"}
    assert violations == []


# ── Foundation marker — software-invariant classification ───────────

def test_foundation_marker_is_registered_in_pyproject() -> None:
    """``pytest.mark.foundation`` must be declared in pyproject.toml
    so pytest doesn't emit ``PytestUnknownMarkWarning`` for it.

    Reads the raw pyproject.toml text to keep the test hermetic —
    no TOML parser needed and the failure mode is a clear string
    miss rather than a nested lookup error.
    """
    pyproject = Path("pyproject.toml").read_text(encoding="utf-8")
    assert '"foundation:' in pyproject, (
        "pyproject.toml markers list must contain a 'foundation:' entry "
        "so pytest registers the marker; otherwise strict-markers mode "
        "would reject it and ORPHEUS cannot enforce the foundation bucket."
    )


def test_marker_to_level_mapping_is_orthogonal() -> None:
    """The conftest marker-to-level helper must classify ``l0``..``l3``
    as uppercase physics levels and ``foundation`` as lowercase.

    Pins the asymmetry that makes ``foundation`` sort below every L<N>
    in the conflicting-marker tiebreak, so a test accidentally carrying
    both ``l1`` and ``foundation`` markers is resolved to ``L1`` (the
    stronger physics claim), not ``foundation``.
    """
    from tests.conftest import _marker_to_level

    assert _marker_to_level("l0") == "L0"
    assert _marker_to_level("l1") == "L1"
    assert _marker_to_level("l2") == "L2"
    assert _marker_to_level("l3") == "L3"
    assert _marker_to_level("foundation") == "foundation"
    # Tiebreak invariant: sorted order must put ``foundation`` below
    # every L<N>, so sorted(...)[ -1] picks the physics level.
    assert sorted(["foundation", "l1"])[-1] == "l1"
    assert sorted(["foundation", "l0"])[-1] == "l0"


def test_vvlevel_literal_includes_foundation() -> None:
    """The ``VVLevel`` type alias must advertise ``foundation`` as a
    legal value so downstream tooling (Sphinx matrix generator, Nexus
    ingest, audit JSON output) accepts foundation tests without
    falling back to ``unmarked``."""
    from typing import get_args

    from tests._harness.registry import VVLevel

    legal = set(get_args(VVLevel))
    assert legal == {"L0", "L1", "L2", "L3", "foundation"}, (
        f"VVLevel Literal mismatch: got {legal}"
    )
