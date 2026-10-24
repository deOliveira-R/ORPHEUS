r"""The realizer table types its two verbs on two spaces (wave 2, WF2.3 slice S3; the hand-off (a) of slice S2b, plan section 6.3).

Spec: ``.claude/plans/wave2_boundary_carve.md`` section 6.3, "Handed over by S2b (2026-10-23)", item (a): ``RealizerTable`` has ONE
space parameter for two verbs. ``realize(declared, site, space)`` takes the METHOD space of the face (``SNMethodSpace`` for SN) and
``sample_inflow(terms, site, space)`` takes the trace's INFLOW HALF-SPACE (SN passes ``trace.inflow_space(face)``), the space
``Gamma_-`` on which a source ``q`` lives. One type parameter cannot type both, so a typed call of the source verb needs
``# type: ignore[arg-type]`` (two of them stand at ``tests/gates/transport/test_boundary_source_from_specs.py``, lines 106 and 108),
and ``SNSource.sample`` re-checks the space's ``role`` at run time (ERR-047). The done-when is that the source verb is typed on its
own space and the two ignores are gone, with ``SN_REALIZERS`` and every S2 gate unchanged.

What is pinned, and the shape it leaves free. The split's SHAPE (a third type parameter, a protocol, a narrower annotation) is the
implementer's: this file pins the two observable consequences that hold for every shape that works and fails for the present one.

* the two ``type: ignore[arg-type]`` comments of ``test_boundary_source_from_specs.py`` are gone (a source-level scan of the file's
  comment tokens, validated by a planted control; the type checker is the instrument that proves the call typed, and it is the
  orchestrator's at integration, because the checker is vacuous in worktrees (#640) except through the documented workaround);
* the annotation of the ``space`` parameter of ``RealizerTable.sample_inflow`` is not the annotation of the ``space`` parameter of
  ``RealizerTable.realize`` (with one type parameter for both they are the same name; after the split they cannot be).

If the implementer finds the split is not the right shape (for example, diffusion samples on its method space, so that only SN
needs a distinct source space), this gate is deleted with that reason in ``NEEDS``, not repaired by weakening.

Red at the base tip (``0757a8d9``): both legs (the file holds two ignores; the two annotations are both ``SpaceT``).

Mutations for the reviewers: the two ignores restored in a copy of the test file; ``sample_inflow``'s ``space`` annotated with the
same type variable as ``realize``'s; a ``type: ignore`` moved onto the line above (a comment token on another line is still found:
the scan reads every comment of the file, not the two lines).
"""

from __future__ import annotations

import inspect
import io
import tokenize
from pathlib import Path

import pytest

from tests.gates.transport._boundary_conditions_helpers import bc

pytestmark = pytest.mark.foundation

_THIS = "tests/gates/transport/test_realizer_table_source_space.py"
_SOURCE_TEST = Path(__file__).with_name("test_boundary_source_from_specs.py")


def arg_type_ignores(source: str) -> list[int]:
    """The line numbers of every comment token carrying ``type: ignore[arg-type]`` (a comment scan, so a string holding the text is not found)."""
    tokens = tokenize.generate_tokens(io.StringIO(source).readline)
    return [t.start[0] for t in tokens if t.type == tokenize.COMMENT and "type: ignore[arg-type]" in t.string]


def test_the_scanner_sees_what_it_must() -> None:
    r"""Positive controls: a planted comment is found on its line, a string literal with the same text is not, and a different
    ignore code is not. Green today; reddens on a scanner that reads nothing or matches strings."""
    planted = 'x = f(a)  # type: ignore[arg-type]  # reason\ny = "type: ignore[arg-type]"\nz = g(b)  # type: ignore[misc]\n'
    assert arg_type_ignores(planted) == [1]


@pytest.mark.rests_on(f"{_THIS}::test_the_scanner_sees_what_it_must")
def test_the_source_verb_test_needs_no_arg_type_ignore() -> None:
    r"""``test_boundary_source_from_specs.py`` types its two calls of ``SN_REALIZERS.sample_inflow`` on the trace half-space without
    ``# type: ignore[arg-type]``. The population is the file's comment tokens (read: the file exists and holds the two calls).

    Red input at the base tip: lines 106 and 108 carry the ignore. Mutation: the ignores restored.
    """
    source = _SOURCE_TEST.read_text()
    assert source.count("sample_inflow(") >= 2, "the file no longer calls the source verb: the gate reads the wrong file"
    assert arg_type_ignores(source) == [], f"{_SOURCE_TEST.name} still carries type: ignore[arg-type] at lines {arg_type_ignores(source)}"


@pytest.mark.rests_on(f"{_THIS}::test_the_scanner_sees_what_it_must")
def test_the_two_verbs_are_not_typed_on_one_space() -> None:
    r"""``RealizerTable.realize`` and ``RealizerTable.sample_inflow`` annotate their ``space`` parameter differently: the method space and
    the source space are two types. (With one type parameter both read ``SpaceT``.)

    Red input at the base tip: both annotations are ``SpaceT``. Mutation: the same type variable on both.
    """
    table = bc.RealizerTable
    realize = inspect.signature(table.realize).parameters["space"].annotation
    sample = inspect.signature(table.sample_inflow).parameters["space"].annotation
    assert realize and sample and realize != sample, f"realize and sample_inflow both annotate space as {realize!r}"
