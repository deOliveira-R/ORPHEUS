r"""SN ``compute_production_rate`` = IRR(νΣf) + (n,2n): the (n,2n) term genuinely adds.

Moved from ``tests/gates/transport/test_integrated_reaction_rate.py`` (leg C) by layer DAG Phase 4: its subject is
``SNSolver.compute_production_rate``, an SN claim; the transport file keeps the integrated reaction rate's own laws.
Every ``macro_xs.synthetic`` mixture has ``Sig2 = 0``, so a library fixture would null the (n,2n) term vacuously; this
row uses ``mixture_with_n2n`` (``tests/gates/transport/_materials.py``) to activate it.

``foundation`` — L0 value correctness.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry.deck import MirrorDeck
from orpheus.mesh import Mesh2D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.problem import SNProblem
from orpheus.sn.solver import SNSolver
from orpheus.transport.boundary_conditions import WallCondition
from orpheus.transport.fields.cross_section_field import CrossSectionField
from orpheus.transport.reaction_rate_functional import IntegratedReactionRate
from tests.gates._exterior import conditions_for
from tests.gates.transport._functional_helpers import require
from tests.gates.transport._materials import mixture_with_n2n

pytestmark = pytest.mark.foundation

_CONDITIONS = {0: WallCondition.vacuum()}


class TestN2NActivationInProductionRate:
    """SN ``compute_production_rate`` = IRR(νΣf) + (n,2n) — the (n,2n) genuinely adds."""

    def test_production_includes_n2n_when_sig2_nonzero(self):
        mat = mixture_with_n2n("2g")  # Sig2 ≠ 0 (macro_xs.synthetic mixtures are all Sig2=0)
        mesh = Mesh2D(
            edges_x=np.linspace(0.0, 0.8, 5), edges_y=np.linspace(0.0, 0.6, 4),
            mat_map=np.zeros((4, 3), dtype=int),
            face_roles={"xmin": MirrorDeck(), "xmax": MirrorDeck(), "ymin": MirrorDeck(), "ymax": MirrorDeck()},
        )
        solver = SNSolver(SNProblem(mesh, Quadrature.lebedev(order=17), {0: mat}, boundary_conditions=conditions_for(mesh, _CONDITIONS)))
        ng = solver.ng
        nx, ny = solver.problem.spatial_shape
        flux = np.random.default_rng(3).uniform(0.1, 1.0, size=(ng, nx, ny))

        fission_only = IntegratedReactionRate(
            CrossSectionField.fission_production_of(solver.problem.mat_xs)
        ).evaluate(flux)
        total = solver.compute_production_rate(flux)
        require(
            total > fission_only * (1.0 + 1e-9),
            f"(n,2n) must ADD to production when Sig2≠0: total={total!r}, "
            f"fission-only IRR={fission_only!r}. Guards the Mode-10 vacuity — every "
            f"macro_xs.synthetic mixture has Sig2=0, so a library fixture would null the term.",
        )
