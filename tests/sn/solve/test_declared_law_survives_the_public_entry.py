r"""A boundary law DECLARED on the geometry reaches a public fixed-source solve.

The declaration channel, landed alongside campaign phase **P4** of
`.claude/plans/archive/affine_boundary_source_channel.md`. Two layers had to be crossed,
and the second is the one that makes the fix worth having:

1. **A ``BC`` tag cannot express a law that carries a function.**
   ``BC`` is ``(kind: str, params: dict[str, float])``
   (``orpheus/geometry/mesh.py``), so a
   :class:`~orpheus.geometry.boundary.PrescribedInflow` whose source is a
   manufactured solution restricted to a face has no tag spelling and never
   will. The channel admits an already-typed
   :class:`~orpheus.geometry.boundary.BoundaryTraceLaw` wherever a tag is
   admitted, and :func:`~orpheus.transport.method._law_from_tag` passes it
   through unchanged.

2. ⭐ **Every public solver entry point rebuilds the method mesh, so anything
   installed on a constructed one is DISCARDED.**
   :func:`~orpheus.sn.solver.solve_sn_fixed_source` takes a raw geometry and
   calls ``_as_sn_mesh(...)``, which constructs a fresh
   :class:`~orpheus.sn.mesh.augmented_mesh.SNMesh`. Before the channel, the only
   way to install a non-tag-expressible law was to mutate a constructed mesh's
   already-resolved ``bc`` dict — and that mesh was not the one the solver used.
   ⟹ the law had to ride the **geometry**, which is what the public API already
   takes. Declaring there makes the rebuild a no-op for the declaration instead
   of a discard, and no ``sn_mesh=`` parameter is needed anywhere.

⚠ **What this module does NOT claim.** It pins that the declaration ARRIVES and
changes the answer — not that the delivered ``q`` is correct in detail. The
value claim is
``test_declared_inflow_reaches_the_rhs.py::test_the_declared_boundary_law_holds_on_the_answer``
(``γ₋ψ|_f == q_f``, ERR-075), which reads the boundary law off the converged
answer. Kept separate because they fail for different reasons: that one for a
mis-delivered source, this one for a declaration that never arrives.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.geometry.boundary import (
    ConstantInflowSource,
    PrescribedInflow,
    VacuumInflow,
)
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solver import solve_sn_fixed_source
from tests.sn._test_helpers import placeholder_materials

pytestmark = pytest.mark.l1

_VALUE = 2.5
_N_ORD = 8
_NG = 2
_NX = 8


def _solve(xmin_declaration, inner: str = "source_iteration") -> np.ndarray:
    """Drive the FULLY PUBLIC entry point — no private helper anywhere."""
    geom = StructuredGeometry(
        geometry="SLB",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(xmin_declaration, BC.vacuum),
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=_NX),))
    quad = Quadrature.gauss_legendre(n_ordinates=_N_ORD)
    source = np.ones((quad.N, _NG, len(mesh.centers)))
    solution = solve_sn_fixed_source(
        placeholder_materials(ng=_NG), mesh, quad, source,
        inner_solver=inner,
    )
    return np.asarray(solution.scalar_flux.values)


class TestTheDeclarationSurvivesTheRebuild:
    r"""⭐ The keystone: a geometry-declared law changes a public solve."""

    @pytest.mark.parametrize("inner", ["source_iteration", "krylov"])
    def test_a_declared_inflow_changes_the_answer(self, inner: str) -> None:
        r"""`[M]` ``|φ_declared − φ_vacuum|_inf = 1.877290`` on this fixture.

        The row that fails if the declaration is dropped anywhere along
        geometry → mesh → ``_as_sn_mesh`` → ``resolve_boundary_conditions`` →
        ``_law_from_tag`` → realizer → ``from_mesh_laws`` → RHS. Before the
        channel it was dropped at the first step (``StructuredGeometry`` refused
        the law outright) and would have been dropped again at ``_as_sn_mesh``.

        Both inner solvers, because they consume the composite RHS through
        different bodies — and because Krylov was unreachable with a declared
        prescribed inflow at all until P3 retired the affine operator.
        """
        declared = _solve(
            PrescribedInflow(source=ConstantInflowSource(value=_VALUE)), inner,
        )
        vacuum = _solve(BC.vacuum, inner)
        moved = float(np.max(np.abs(declared - vacuum)))
        assert moved > 1e-6, (
            f"[{inner}] a DECLARED prescribed inflow left the answer identical "
            f"to vacuum (|Δφ| = {moved:.3e}) — the declaration was discarded "
            f"somewhere between the geometry and the solve, which is exactly "
            f"the defect this channel closes."
        )

    def test_the_tag_and_the_law_agree_where_BOTH_can_express_the_same_thing(
        self,
    ) -> None:
        r"""⭐ The two declaration arms are one channel, not two.

        ``BC.vacuum`` (a tag) and :class:`VacuumInflow` (the law it parses to)
        name the same physics, so they must produce the same float program.
        BIT-identical, not merely close: the tag arm's only extra work is the
        registry lookup that *builds* this very law, so any difference means the
        arms diverged after the parse.

        This is the row that would catch a law arm which quietly bypassed some
        step the tag arm performs — the failure mode a "does it change the
        answer?" test cannot see.
        """
        np.testing.assert_array_equal(
            _solve(BC.vacuum), _solve(VacuumInflow()),
            err_msg=(
                "declaring the TAG and declaring the LAW it parses to gave "
                "different answers — the two arms of the declaration channel "
                "have diverged"
            ),
        )


class TestTheChannelDoesNotWeakenTheTagPath:
    """The tag remains first-class; the law arm is an addition, not a detour."""

    @pytest.mark.parametrize("tag", [BC.vacuum, BC.reflective])
    def test_a_tag_still_resolves(self, tag) -> None:
        """Both shipped tags still drive a solve to a finite answer.

        Cheap, but it is the row that reddens if the pass-through arm were
        placed so that it swallowed tags too (e.g. an ``isinstance`` test broad
        enough to match ``BC``).
        """
        phi = _solve(tag)
        assert np.all(np.isfinite(phi))
        assert float(np.max(np.abs(phi))) > 0.0

    def test_an_unsupported_tag_still_refuses(self) -> None:
        r"""⭐ The law arm must not become a hole in the admission table.

        ``BOUNDARY_OPERATOR_REGISTRY`` gates which TAGS a method admits, and
        ``BC("white")`` is not among SN's (#189). The pass-through arm returns
        early for law OBJECTS only — if it were reached for tags, an unsupported
        tag would sail through the registry check and fail later, or worse,
        silently.
        """
        with pytest.raises(ValueError, match="does not support boundary"):
            _solve(BC("white"))
