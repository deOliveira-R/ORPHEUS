r"""#711 (W8a): within a process one ``GeneratingMeasure`` name is one content.

``GeneratingMeasure`` compares by ``name`` alone (``support`` and ``exact_recurrence`` are outside equality, the
latter ``compare=False``); nothing checks the implementation against the name. The refusal is an intern table
``_NAME_CONTENT`` from the name to ``_MeasureContent(support, alpha, beta)`` (the bytes of the first 8 float64
recurrence coefficients) and ONE ``raise ValueError`` in ``GeneratingMeasure.__post_init__`` tagged
``ELEGANCE-DEBT[guard] #711 (W8a)``. Specification: ``.claude/plans/ontology_w0_spec.md`` 5.1 to 5.3.

Law: name implies content. The converse FAILS and is stated: ``LEGENDRE.on(-1.0, 1.0)`` has the support and the
bit-equal ``recurrence(8)`` of ``LEGENDRE`` and is not equal to it. Limit: two measures that agree on 8
coefficients and differ later are not told apart.

Activation: ``dataclasses.replace(CHEBYSHEV_T, name=LEGENDRE.name)`` differs from LEGENDRE at ``beta[0]``, the
coefficient the probe reads first, and its ``gauss(3)`` weights differ (the first red also reads equal ``==`` and
equal hashes today). The three classical Jacobi pairs are the admitted legs that stay green; they are built through
``_jacobi.__wrapped__`` because ``jacobi`` is ``lru_cache``d (the second ``jacobi(0.0, 0.0)`` of a process returns
the first object without running ``__post_init__``).

Mutations: the refusal replaced by ``assert False, (`` in a copy (leg 1 constructs under ``-O``); the content
record made to carry ``id(self.exact_recurrence)`` in a copy (the classical-pair admitted legs are refused). The
content is built inline, so an in-process patch could only replace the whole method, a twin of it.

Every assertion is a function call or an ``if ... pytest.fail`` (``python -O`` safe).
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from orpheus.numerics import generating_measure as gm
from orpheus.numerics.generating_measure import CHEBYSHEV_T, CHEBYSHEV_U, LEGENDRE
from tests.gates import _refusal_gates as rg

pytestmark = pytest.mark.foundation

FRAGMENT = "already names a different measure"
_FILE = "orpheus/numerics/generating_measure.py"

_CHILD_IMPORTS = (
    "import dataclasses\n"
    "from orpheus.numerics import generating_measure as gm\n"
    "from orpheus.numerics.generating_measure import CHEBYSHEV_T, LEGENDRE\n"
)
_CHILD_INPUT = "dataclasses.replace(CHEBYSHEV_T, name=LEGENDRE.name)"
_CHILD_CLASSICAL = (
    "for a, b in ((0.0, 0.0), (-0.5, -0.5), (0.5, 0.5)):\n"
    "    gm._jacobi.__wrapped__(a, b)\n"
)

#: ``(a, b)`` of the three classical pairs, and the module-level constant each equals.
_CLASSICAL = [(0.0, 0.0, LEGENDRE), (-0.5, -0.5, CHEBYSHEV_T), (0.5, 0.5, CHEBYSHEV_U)]


@pytest.fixture
def isolated_name_table(monkeypatch):
    """A copy of the process-global table, so no record leaks between tests. (Absent before the refusal lands: nothing to isolate.)"""
    table = getattr(gm, "_NAME_CONTENT", None)
    if table is not None:
        monkeypatch.setattr(gm, "_NAME_CONTENT", dict(table))


# ----------------------------------------------------------------- admitted legs


@pytest.mark.parametrize(("a", "b", "constant"), _CLASSICAL, ids=["jacobi(0,0)", "jacobi(-1/2,-1/2)", "jacobi(1/2,1/2)"])
def test_the_three_classical_pairs_agree_with_their_constants_and_are_admitted(a, b, constant, isolated_name_table):
    """``jacobi(0, 0)``, ``jacobi(-1/2, -1/2)`` and ``jacobi(1/2, 1/2)`` agree with ``LEGENDRE``, ``CHEBYSHEV_T`` and
    ``CHEBYSHEV_U`` on the name, the support and the bytes of ``recurrence(8)``: admitted by the law."""
    built = gm._jacobi.__wrapped__(a, b)
    if built.name != constant.name or built.support != constant.support:
        pytest.fail(f"{built.name!r} on {built.support!r} is not {constant.name!r} on {constant.support!r}")
    for mine, theirs in zip(built.recurrence(8), constant.recurrence(8)):
        if np.asarray(mine, dtype=float).tobytes() != np.asarray(theirs, dtype=float).tobytes():
            pytest.fail("the classical pair and its constant differ in the first 8 recurrence coefficients")


def test_one_name_two_contents_control_the_two_measures_really_differ(isolated_name_table):
    """CONTROL: the refused input is a different measure, not an equal one in disguise."""
    _, beta_t = CHEBYSHEV_T.recurrence(8)
    _, beta_l = LEGENDRE.recurrence(8)
    if np.array_equal(beta_t, beta_l):
        pytest.fail("CONTROL INVALID: CHEBYSHEV_T and LEGENDRE share their recurrence")


def test_content_does_not_imply_name_the_remap_is_unequal_to_its_source(isolated_name_table):
    """The converse FAILS, stated and gated: equality is strictly finer than content (W8a is where identity becomes content)."""
    remapped = LEGENDRE.on(-1.0, 1.0)
    if remapped == LEGENDRE:
        pytest.fail("the remap onto its own interval became equal to its source: the converse now holds, retire this leg")
    for mine, theirs in zip(remapped.recurrence(8), LEGENDRE.recurrence(8)):
        if np.asarray(mine, dtype=float).tobytes() != np.asarray(theirs, dtype=float).tobytes():
            pytest.fail("CONTROL INVALID: the remap's recurrence is not bit-equal to its source's")


# ----------------------------------------------------------------- refused input, in process


def test_one_name_for_two_measures_is_refused(isolated_name_table):
    """First red: ``replace(CHEBYSHEV_T, name=LEGENDRE.name) == LEGENDRE`` is True today (equal hashes, different
    ``gauss(3)`` weights)."""
    with pytest.raises(ValueError, match=FRAGMENT) as info:
        dataclasses.replace(CHEBYSHEV_T, name=LEGENDRE.name)
    if LEGENDRE.name not in str(info.value):
        pytest.fail(f"the refusal does not name the name {LEGENDRE.name!r}: {info.value}")


def test_the_refused_input_is_a_different_measure_today_by_its_gauss_weights():
    """ACTIVATION: ``gauss(3)`` of the two differ, so the identity the name claims is false; this reads the
    implementation, not the refusal, and is green before and after."""
    if np.array_equal(CHEBYSHEV_T.gauss(3).weights, LEGENDRE.gauss(3).weights):
        pytest.fail("CONTROL INVALID: the two measures give one 3-point rule")


# ----------------------------------------------------------------- cold child, -O


def test_the_unmutated_child_builds_the_classical_pairs():
    """CONTROL of the child legs."""
    child = rg.run_child(_CHILD_IMPORTS, _CHILD_CLASSICAL)
    if not rg.imported_from(child, rg.import_root()):
        pytest.fail(f"the child did not import orpheus from {rg.import_root()}:\n{child.stdout}\n{child.stderr}")
    if not rg.constructed(child, FRAGMENT):
        pytest.fail(f"a classical pair was refused in the child:\n{child.stderr}")


def test_the_refusal_survives_python_O():
    child = rg.run_child(_CHILD_IMPORTS, _CHILD_INPUT)
    if not rg.imported_from(child, rg.import_root()):
        pytest.fail(f"the child did not import orpheus from {rg.import_root()}:\n{child.stdout}\n{child.stderr}")
    if not rg.refused(child, FRAGMENT):
        pytest.fail(f"the renamed measure was not refused under -O (exit {child.returncode}):\n{child.stderr}")


# ----------------------------------------------------------------- mutations


def test_mutation_the_refusal_replaced_by_an_assert_constructs_under_python_O(tmp_path):
    """Leg (1): the ``raise`` after the #711 tag becomes ``assert False, (``; the renamed measure constructs."""
    honest = rg.run_child(_CHILD_IMPORTS, _CHILD_INPUT)
    if not rg.refused(honest, FRAGMENT):
        pytest.fail(f"the unmutated tree does not refuse the input, so the mutation proves nothing:\n{honest.stderr}")
    root = rg.copy_orpheus_py(rg.import_root(), tmp_path)
    rg.remove_refusal(root, _FILE, 711)
    mutant = rg.run_child(_CHILD_IMPORTS, _CHILD_INPUT, root=root)
    if not rg.imported_from(mutant, root):
        pytest.fail(f"the mutant did not import from the copy (never red evidence):\n{mutant.stdout}\n{mutant.stderr}")
    if not rg.constructed(mutant, FRAGMENT):
        pytest.fail(f"the mutant still refused the input (exit {mutant.returncode}):\n{mutant.stderr}")


def test_mutation_a_content_that_carries_the_callable_refuses_the_classical_pairs(tmp_path):
    """``_MeasureContent(self.support, `` replaced by ``_MeasureContent((self.support, id(self.exact_recurrence)), ``:
    the content now distinguishes two builds of one measure, so the admitted ``jacobi(0, 0)`` against ``LEGENDRE``
    is refused with the fragment. The honest tree admits the three pairs (read first)."""
    honest = rg.run_child(_CHILD_IMPORTS, _CHILD_CLASSICAL)
    if not rg.constructed(honest, FRAGMENT):
        pytest.fail(f"the unmutated tree refuses a classical pair, so the mutation proves nothing:\n{honest.stderr}")
    root = rg.copy_orpheus_py(rg.import_root(), tmp_path)
    rg.replace_once(
        root, _FILE, "_MeasureContent(self.support, ", "_MeasureContent((self.support, id(self.exact_recurrence)), ",
    )
    mutant = rg.run_child(_CHILD_IMPORTS, _CHILD_CLASSICAL, root=root)
    if not rg.imported_from(mutant, root):
        pytest.fail(f"the mutant did not import from the copy (never red evidence):\n{mutant.stdout}\n{mutant.stderr}")
    if not rg.refused(mutant, FRAGMENT):
        pytest.fail(f"the mutant still admits the classical pairs (exit {mutant.returncode}):\n{mutant.stderr}")
