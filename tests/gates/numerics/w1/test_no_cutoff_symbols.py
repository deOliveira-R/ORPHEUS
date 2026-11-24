r"""G02c: no cutoff on a direction cosine survives in ``orpheus/`` (N02, spec 7.3).

Two ``ast`` arms over source files, with no import of the modules read (the arm reads ``sn/`` and ``mesh/`` sources by path):

(i) none of the five cutoff names is defined, imported or read anywhere in ``orpheus/`` (``derivations/`` excluded);
(ii) in the 8 production files that held a cutoff, the ``Compare`` nodes that carry a tolerance operand equal the allow-list (one
comparison, ``reading <= _GAUGE_AUDIBLE_FLOOR`` of ``loss_kernel_gauge.py``), keyed on the file and the source segment, never a line.

What it cannot see: a cutoff under a new name in a file outside the eight, a cutoff held in a parameter or a local, a tolerance passed
to ``np.isclose`` or ``np.allclose``, a name not ending ``EPS`` or ``TOL``, a literal above ``1e-6``, ``moc/geometry.py:143`` and
``:153`` (they stay), and an exact comparison of a cosine (G02g and issue #752 are about those).

Exact on any host: parse trees only. First red today: (i) 5 of 5 names are defined; (ii) 13 comparisons carry a tolerance operand in
the 8 files, 12 of them not on the allow-list. Each arm is run on a mutant source (a copy of ``mesh/axis.py`` with a cutoff appended) in
a row that is green on any tree.
"""

from __future__ import annotations

import pytest

from tests.gates import _half_range_helpers as H

pytestmark = pytest.mark.foundation


def _package_sources() -> list[tuple[str, str]]:
    root = H.orpheus_root()
    return [(str(p.relative_to(root)), p.read_text()) for p in sorted(root.rglob("*.py")) if "derivations" not in p.relative_to(root).parts]


def test_no_cutoff_symbols() -> None:
    sources = _package_sources()
    assert len(sources) > 200, f"the census reads {len(sources)} files: the filter dropped the tree"
    assert any(rel == "mesh/axis.py" for rel, _ in sources), "mesh/axis.py is not in the population"
    hits = [(rel, name, line) for rel, src in sources for name, line in H.cutoff_name_hits(src)]
    assert not hits, f"{len(hits)} occurrences of the five cutoff names in {len(sources)} files (file, name, line): {hits}"


def test_no_comparison_with_a_tolerance_survives_in_the_files_that_held_a_cutoff() -> None:
    root = H.orpheus_root()
    found = {
        (rel, seg) for rel in H.CUTOFF_FILES for _, seg in H.tolerance_compares((root / rel).read_text())
    }
    assert found == H.TOLERANCE_ALLOW_LIST, (
        f"{len(found - H.TOLERANCE_ALLOW_LIST)} comparisons carry a tolerance beyond the allow-list "
        f"(file, segment): {sorted(found - H.TOLERANCE_ALLOW_LIST)}; allow-list entries not found: "
        f"{sorted(H.TOLERANCE_ALLOW_LIST - found)}"
    )


def test_the_name_arm_reads_an_appended_cutoff() -> None:
    r"""Mutation (i): ``_OUTFLOW_COSINE_TOL = 1e-15`` appended to a copy of ``mesh/axis.py`` reads as a definition; the read form and the import form too."""
    src = (H.orpheus_root() / "mesh/axis.py").read_text()
    base = H.cutoff_name_hits(src)
    assert H.cutoff_name_hits(src + "\n_OUTFLOW_COSINE_TOL = 1e-15\n") != base
    assert H.cutoff_name_hits("y = c > TANGENTIAL_EPS\n") == [("TANGENTIAL_EPS", 1)]
    assert H.cutoff_name_hits("from a import _MU_DIRECTION_EPS\n") == [("_MU_DIRECTION_EPS", 1)]
    assert H.cutoff_name_hits("y = m._OCTANT_SIGN_EPS\n") == [("_OCTANT_SIGN_EPS", 1)]
    assert H.cutoff_name_hits("y = c > 0.0\n") == []


def test_the_comparison_arm_reads_a_cutoff_under_a_new_name() -> None:
    r"""Mutation (ii): ``_NEW_TOL = 1e-15`` and ``y = c > _NEW_TOL`` appended to a copy of ``mesh/axis.py`` add one comparison to the arm."""
    src = (H.orpheus_root() / "mesh/axis.py").read_text()
    before = H.tolerance_compares(src)
    after = H.tolerance_compares(src + "\n_NEW_TOL = 1e-15\ny = c > _NEW_TOL\n")
    assert len(after) == len(before) + 1, (len(before), len(after))
    assert H.tolerance_compares("y = c > 1e-15\n") != []
    assert H.tolerance_compares("y = c > np.finfo(float).eps\n") != [], "the inline .eps form is read"
    assert H.tolerance_compares("E = 4 * np.finfo(np.float64).eps\ny = c > E\n") != []
    assert H.tolerance_compares("y = c > 0.0\n") == []
