r"""The method-mesh's twin binding duties live ONCE, in the shared method layer (wave 2, slice S4, hand-off (a) of slice S3).

Spec: ``.claude/plans/wave2_boundary_carve.md`` section 6.4, "Handed over by S3 (2026-10-23)" (a), and 11.4 ruling 5.

What is gated
=============

A method-mesh (a :class:`~orpheus.transport.mesh.material_mesh.MaterialMesh` that owns a ``binding``: the declared exterior
conditions joined to its faces) owes the identity and the re-spelling machinery two duties, and until this commit both were
spelled twice, identically, in ``orpheus/sn/problem.py`` and ``orpheus/diffusion/augmented_mesh.py``:

* **the key**: ``_contractibility_key`` extends the material mesh's key by ``conditions_key(self.binding.declared)``, so two
  problems with the same mesh and different conditions are different phase spaces;
* **the re-spelling**: the methods that rebuild the problem (``with_cross_sections``, ``with_scattering_order``) re-pass
  ``boundary_conditions=self.binding.declared``, because a rebuilt carrier whose tags are unbound is refused.

CP's method-mesh would have been the third copy (and MoC's and MC's the fourth and fifth), which is a stop signal under the
cardinal rule on architecture. The lift puts both duties in ONE base of the method layer, ``orpheus.transport.method.MethodMesh``;
``SNProblem`` and ``DiffusionMesh`` derive from it and spell neither.

The ladder
==========

::

    test_the_scanner_sees_what_it_must                     foundation: the AST census reads planted members
        ^
    test_the_duties_are_spelled_once_in_the_method_layer   structural: 0 sites in sn/diffusion, 1 key site in the base
        ^
    test_the_key_carries_the_conditions[sn|diffusion]      behaviour: different conditions, different phase space
    test_a_respelled_problem_keeps_its_conditions[...]     behaviour: with_cross_sections re-passes the declared map

Red at the base tip (``d4da4328``)
==================================

* ``test_the_duties_are_spelled_once_in_the_method_layer``: the census reads 2 ``conditions_key(`` calls (one per method module) and 5 reads of
  ``self.binding.declared`` (SN 3, diffusion 2) in the two method modules, and no definition of ``MethodMesh`` in ``method.py``. This is the gate's own red input, in the tree today.
* The four behavioural rows are green at the base (SN and diffusion carry the duty today): they are the "keep every SN and
  diffusion gate green" leg, and the instruments for the mutation below.

Mutation for the reviewers (M-key-drops-conditions): in the shared base, ``_contractibility_key`` returns the material mesh's
key WITHOUT ``conditions_key(...)``. Red: both ``test_the_key_carries_the_conditions`` rows (``vac != white`` becomes False),
beside the diffusion row ``tests/gates/diffusion/test_augmented_mesh.py::TestConditionsAreThePhaseSpace`` that already pinned it
for diffusion alone. Second mutation (M-respell-drops-conditions): the shared re-spelling keywords return ``{}``: both
``test_a_respelled_problem_keeps_its_conditions`` rows raise ``UnboundTagError``.

Moved from tests/gates/transport/test_method_mesh_binding.py by layer DAG Phase 4 (2026-11-01): these rows build an
SNProblem and a DiffusionMesh, so they are cross-method; the key row rests on the transport file's scanner row.
"""

from __future__ import annotations

from typing import Any

import pytest

_THIS = "tests/gates/cross_method/test_method_mesh_binding.py"
#: The scanner row the key row rests on stayed in transport.
_SCANNER = "tests/gates/transport/test_method_mesh_binding.py"


def _sn(condition: Any) -> Any:
    from orpheus.data.macro_xs.synthetic import get_mixture
    from orpheus.geometry import StructuredGeometry
    from orpheus.geometry.deck import MirrorDeck
    from orpheus.mesh import CellsByCount, Mesher
    from orpheus.numerics.quadrature import Quadrature
    from orpheus.sn.problem import SNProblem

    mesh = Mesher(StructuredGeometry.slab((0.0, 2.0), (0,), left=MirrorDeck(), right=0)).partition(CellsByCount.uniform_width(4)).mesh
    return SNProblem(mesh, Quadrature.gauss_legendre(4), {0: get_mixture("B", "2g")}, boundary_conditions={0: condition})


def _diffusion(condition: Any) -> Any:
    from orpheus.data.macro_xs.synthetic import get_mixture
    from orpheus.diffusion.augmented_mesh import DiffusionMesh
    from orpheus.geometry import StructuredGeometry
    from orpheus.geometry.deck import MirrorDeck
    from orpheus.mesh import CellsByCount, Mesher

    mesh = Mesher(StructuredGeometry.slab((0.0, 2.0), (0,), left=MirrorDeck(), right=0)).partition(CellsByCount.uniform_width(4)).mesh
    return DiffusionMesh(mesh, {0: get_mixture("B", "2g")}, boundary_conditions={0: condition})


_BUILDERS = {"sn": _sn, "diffusion": _diffusion}


@pytest.mark.l1
@pytest.mark.rests_on(f"{_SCANNER}::test_the_scanner_sees_what_it_must")
@pytest.mark.parametrize("method", list(_BUILDERS))
def test_the_key_carries_the_conditions(method: str) -> None:
    r"""Two problems over the same mesh with different conditions are different phase spaces, and two with the same
    condition are the same one (the control: the key is not merely unequal for everything).

    Red input: none at the base tip (both methods carry the duty today); the instrument of the mutation M-key-drops-conditions in
    the shared base, which turns ``vac != white`` into False for BOTH rows.
    """
    from orpheus.transport.boundary_conditions import WallCondition

    build = _BUILDERS[method]
    vac, white, white_again = build(WallCondition.vacuum()), build(WallCondition.white(0.5)), build(WallCondition.white(0.5))
    assert vac != white and not vac.same_phase_space(white), f"{method}: vacuum and white(0.5) are one phase space"
    assert white == white_again and white.same_phase_space(white_again)


@pytest.mark.l1
@pytest.mark.rests_on(f"{_THIS}::test_the_key_carries_the_conditions")
@pytest.mark.parametrize("method", list(_BUILDERS))
def test_a_respelled_problem_keeps_its_conditions(method: str) -> None:
    r"""``with_cross_sections`` rebuilds the problem over the same phase space: the declared conditions are carried (without them
    the rebuilt carrier's tag is unbound, ``UnboundTagError``), and the rebuilt problem has the original's key (same sigma_t
    datum).

    Red input: none at the base tip. Mutation M-respell-drops-conditions (the shared re-spelling keywords return ``{}``): both
    rows raise ``UnboundTagError``.
    """
    from orpheus.transport.boundary_conditions import WallCondition

    problem = _BUILDERS[method](WallCondition.white(0.5))
    again = problem.with_cross_sections(problem.sigma_t_cell.copy())
    assert again.binding.declared == problem.binding.declared
    assert again._contractibility_key == problem._contractibility_key
