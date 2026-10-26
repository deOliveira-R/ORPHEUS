"""``realize_boundary`` completes no boundary declaration by default (the transport layer's half of S3.10).

Moved from ``tests/gates/mesh/test_mesh1d.py`` by claim in the layer-DAG re-filing (Phase 1b,
2026-10-26): the row reads only ``orpheus.transport.method``.
"""

from __future__ import annotations

import pytest

pytestmark = pytest.mark.foundation


class TestTheDefaultsAreGone:
    """S3.10, step 3c (the user's ruling of 2026-09-29: ``None`` retires as a
    boundary declaration everywhere), re-posed at wave 2 S2b. No entry takes a
    singular ``boundary_condition`` (the declared conditions travel as the keyword
    ``boundary_conditions``, keyed by tag), ``_apply_default_bcs`` is gone, and
    ``realize_boundary`` completes no declaration by default. Claim kind: THEOREM
    (the retirement's witness, ``retirement-audit`` D.18).

    The resolver rows of the earlier spelling (a ``None`` law slot built around the
    constructor, and an AST scan for ``BC``) are re-posed below on the new
    vocabulary; the keyed refusal of a ``None`` role is gated at construction
    (``tests/gates/mesh/test_mesh2d_face_roles.py::TestTheOneElementParser``).

    First red: a default ``WallCondition`` constructed in ``realize_boundary``
    (the structural row). The signature row and the name row (no SN entry takes
    ``boundary_condition``; ``_apply_default_bcs`` is gone) are SN claims and moved
    to ``tests/gates/sn/solve/test_no_entry_takes_a_boundary_condition.py``; this
    structural row is a claim about the transport layer's ``realize_boundary`` and
    moved here from ``tests/gates/mesh/test_mesh1d.py`` (the layer-DAG re-filing,
    Phase 1b, 2026-10-26).
    """

    # ``test_the_resolver_completes_no_undeclared_face`` is RETIRED with its subject (wave 2 S2b). It built an axis around its
    # constructor holding ``None`` in a law slot and required ``resolve_boundary_conditions`` to realise no law on that face. The
    # resolver no longer reads the axes: ``realize_boundary`` reads ``method.binding``, which ``bind`` parses from the carrier's roles
    # at construction (a ``None`` or a non-role is refused there, gated by ``TestTheOneElementParser``), so the state the row built
    # cannot reach the realization body and a stand-in that fed it ``None`` would be testing a path no caller has. The structural
    # leg below (the body completes nothing by default) is kept, re-spelled on the new default-fill vocabulary.

    def test_the_realization_body_constructs_no_declaration(self):
        """``realize_boundary`` spells no ``WallCondition`` and no deck: a default fill (``... or WallCondition.vacuum()``, the spelling
        every earlier default took, in the vocabulary of the day) would construct one inside the body.

        The scan is run first on a planted body that does spell the default, so a zero below is a measurement and not a blind pattern."""
        import ast
        import inspect
        import textwrap

        from orpheus.transport.method import realize_boundary

        names = {"WallCondition", "MirrorDeck", "TranslationDeck"}

        def spelled(source: str) -> list[str]:
            return [
                ast.unparse(node) for node in ast.walk(ast.parse(textwrap.dedent(source)))
                if isinstance(node, ast.Name) and node.id in names
            ]

        assert spelled("def f(b):\n    return b.get('xmin') or WallCondition.vacuum()") == ["WallCondition"]  # the positive control
        assert spelled(inspect.getsource(realize_boundary)) == []
