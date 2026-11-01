r"""Carrier-admission gates: ``mesh is None`` has ONE meaning (S7 G7.1–G7.3, re-posed at the CS4c coda).

History. CS4b S7 (2026-08-24) un-welded a sentinel that carried TWO
meanings — the mesh-less infinite-medium 1-cell carrier
(``MaterialMesh.from_materials``) and the d≥3 axis-native carrier — with
typed refusals (G7.1 the SN promotion refusal as a real ``raise`` under the
canonical ``-O``; G7.2 one ``areas`` message per arm; G7.3 the two states
discriminated by ``ndim``). The CS4c coda (C1 ``5caad3d6``, C2 this commit)
retired the fabricated carrier and its factory — the infinite-medium
problem poses on ``HomogeneousProblem.space`` and builds no carrier — so
the arms that served it (the SN promotion refusal, the ``areas`` "no faces
at all" arm, the diffusion bounded-geometry refusal) had NO reachable
input left and retired with it (memo H-1: `[M]` no producer of
``mesh is None and ndim == 1`` exists — ``SNProblem.from_axes`` synthesizes
a legacy adapter at d≤2 and leaves ``mesh = None`` only at d≥3;
``MaterialMesh.__init__`` always carries a mesh).

What is gated now:

G7.1 — a real carrier promotes (the positive leg survives; its negative
leg's subject is gone).
G7.2 — the two surviving ``areas`` arms each name their own case, pinned
by their shortest distinctive fragments; collapsing them reddens a row.
G7.3 — the SINGLETON law that replaced the discrimination law: every
d≤2 constructor carries a mesh, and ``mesh is None`` ⟹ ``ndim ≥ 3``.
This is the THEOREM the retirement rests on (``coding-standards``: what
kept the duplicate arms distinguishable is now the sole guarantor, and
it owes a witness) — a positive control (the d=3 carrier constructs with
``mesh is None``) and the d=1 / d=2 axis-native constructions that
synthesize their adapter.
R1 — the retirement is UNSPELLABLE: no class in the carrier hierarchy
has a ``from_materials``. (The energy-arm rule that once shared the name,
``EnergyAxis.from_materials``, is now ``orpheus.data.materials.energy_axis_of``,
gated at ``test_kernels.py``; with the homonym gone there is no name left to
confuse, so this file no longer pins it.)

Moved from tests/gates/transport/test_material_mesh_admission.py by layer DAG Phase 4 (2026-11-01): this file holds
R1's carrier-hierarchy row only, which names MaterialMesh, SNProblem and DiffusionMesh, so it is cross-method; the
other rows are in transport and sn/mesh/test_problem_admission.py.
"""

from __future__ import annotations

import pytest

from orpheus.diffusion import DiffusionMesh
from orpheus.sn.problem import SNProblem
from orpheus.transport.mesh.material_mesh import MaterialMesh

pytestmark = [pytest.mark.foundation]


def _require(condition: bool, message: str) -> None:
    if not condition:
        pytest.fail(message)


class TestR1TheFactoryIsUnspellable:
    def test_no_carrier_class_has_a_from_materials(self):
        """The retirement is structural: the name is absent from every class in
        the carrier hierarchy (a subclass could not have inherited it either).
        Re-adding the factory anywhere in the hierarchy reddens this."""
        for cls in (MaterialMesh, SNProblem, DiffusionMesh):
            _require(
                not hasattr(cls, "from_materials"),
                f"{cls.__name__}.from_materials exists — the fabricated carrier's factory is back",
            )

