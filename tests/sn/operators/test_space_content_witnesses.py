r"""Manufactured witnesses for the S3 identity arms that had NO catcher.

CS4b verification plan §1.2/§11.1: 8 of the 22 identity-guard arms were
measured 0-red over a ≥3936-row denominator — a rewrite of an unwitnessed
guard is indistinguishable from its deletion (plan-authoring §6c), so the
S3 re-key OWES each arm a witness in the same step. This file carries the
arms whose operators had no natural refusal row anywhere (O4, O5, O8, O11,
O13, F5); O12's row rides the windowing suite's own fixture and O17's the
streaming-collision file, next to their siblings.

Every row is the ADMISSION pattern (vv #11, both legs where cheap): the
content-differing operand refuses with the ``space-content`` family token
(or the guard's own named message), and the twin-carrier leg — the
correctly-blind half of the F2 doctrine — is covered by the re-derived
sibling rows in each operator's home file.

Discriminators come from §4's table: stretched volumes for bulk content,
a graded radial mesh for ray content.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.operators.boundary import SNBoundaryOperator
from orpheus.sn.operators.radial_characteristic import (
    RadialCharacteristicReconstruction,
    RadialCharacteristicSeeding,
)
from orpheus.sn.coupled_system import build_within_group_system
from orpheus.sn.solver import (
    SNSolver,
    evaluate_residual,
    solve_sn_adjoint_fixed_source,
)
from orpheus.sn.solution import Solution
from orpheus.sn.operators.streaming import StreamingOperator
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.scalar_flux import ScalarFlux
from orpheus.transport.full_field import FullField
from orpheus.transport.source_sinks import (
    AngularBoundarySourceSink,
    AngularSourceSink,
)
from orpheus.transport.timed_full_field import TimedFullField
from tests.sn._test_helpers import placeholder_materials

pytestmark = pytest.mark.foundation


def _slab(*, width: float = 1.0, nx: int = 4, ng: int = 2) -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.0, width, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    return SNMesh(
        mesh, Quadrature.gauss_legendre(4), placeholder_materials(ng=ng)
    )


def _sphere(*, power: float = 1.0, nx: int = 5, ng: int = 2) -> SNMesh:
    """Seed-carrying sphere; ``power != 1`` grades the radii (the ray-content
    discriminator — same shape, different Δr, different ray metric)."""
    radii = 4.0 * (np.arange(nx + 1) / nx) ** power
    mesh = Mesh1D(
        edges=radii,
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_right=BC("vacuum"),
    )
    return SNMesh(
        mesh, Quadrature.gauss_legendre(4), placeholder_materials(ng=ng)
    )


def _composite(sn: SNMesh) -> FullField:
    return FullField(
        interior=AngularFlux.zeros_on(sn),
        boundary=AngularBoundaryFlux.zeros_on(sn),
    )


class TestO4DiffusionBoundaryOperator:
    def test_apply_refuses_content_differing_interior(self):
        from orpheus.diffusion import DiffusionBoundaryOperator, DiffusionMesh
        from orpheus.derivations.common.xs_library import get_mixture
        from orpheus.transport.fields.scalar_boundary_flux import (
            ScalarBoundaryFlux,
        )

        def _dmesh(width: float) -> DiffusionMesh:
            m = Mesh1D(
                np.linspace(0.0, width, 5), np.zeros(4, dtype=int),
            )
            return DiffusionMesh(m, {0: get_mixture("A", "2g")})

        op = DiffusionBoundaryOperator(_dmesh(10.0))
        stretched = _dmesh(20.0)
        foreign = FullField(
            interior=ScalarFlux.zeros_on(stretched),
            boundary=ScalarBoundaryFlux.zeros_on(stretched),
        )
        with pytest.raises(ValueError, match="space-content"):
            op.apply(foreign)


class TestO5AdjointEntryDetector:
    def test_composite_detector_on_foreign_content_refuses(self):
        sn_foreign = _slab(width=2.0)
        q_star = FullField(
            interior=AngularSourceSink.zeros_on(sn_foreign),
            boundary=AngularBoundarySourceSink.zeros_on(sn_foreign),
        )
        base = _slab()
        with pytest.raises(ValueError, match="space-content"):
            solve_sn_adjoint_fixed_source(
                base.materials,
                Mesh1D(
                    edges=np.linspace(0.0, 1.0, 5),
                    mat_ids=np.zeros(4, dtype=int),
                    coord=CoordSystem.CARTESIAN,
                    bc_left=BC("vacuum"),
                    bc_right=BC("vacuum"),
                ),
                Quadrature.gauss_legendre(4),
                q_star,
            )


class TestO8SolutionRayMember:
    def test_ray_member_with_foreign_ray_content_refuses(self):
        sn = _sphere()
        graded = _sphere(power=1.5)
        from orpheus.transport.radial_characteristic_field import (
            RadialCharacteristicField,
        )

        member = RadialCharacteristicField.source_zeros_on(graded)
        # A flux-role member is what a Solution carries; the source-zeros
        # factory gives the cheapest structurally-complete member — the
        # guard reads only its block SPACES, which is the claim.
        psi = TimedFullField.zeros(
            interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn,
        )
        phi = ScalarFlux.zeros_on(sn)
        with pytest.raises(ValueError, match="ray spaces"):
            Solution(
                angular_flux=psi, scalar_flux=phi, mesh=sn,
                radial_characteristic=member,  # type: ignore[arg-type]
            )


class TestO11ReconstructionTranspose:
    def test_apply_transpose_refuses_foreign_cotangent(self):
        sn = _sphere()
        graded = _sphere(power=1.5)
        cot = _composite(graded)
        with pytest.raises(ValueError, match="space-content"):
            RadialCharacteristicReconstruction(sn).apply_transpose(cot)

    def test_forward_arm_sibling_still_witnessed(self):
        """O10's forward twin (the classic forward/transpose asymmetry the
        plan flagged): the seeding guard's graded-refusal rows live in
        test_psi_half_coupling; this sibling pins the SAME guard through
        the seeding entry so the pair is visible side by side."""
        sn = _sphere()
        graded = _sphere(power=1.5)
        with pytest.raises(ValueError, match="space-content"):
            RadialCharacteristicSeeding(sn).apply_transpose(
                _composite(graded)
            )


class TestO12WindowingAnalysis:
    def test_apply_refuses_content_differing_composite(self):
        from orpheus.sn.operators.windowing import BulkAnalysisOperator
        from orpheus.transport.operators.scattering import ScatteringOperator
        from tests.sn._test_helpers import material_xs_from_raw

        sn = _slab()
        p0 = np.array([[0.2, 0.0], [0.05, 0.18]])
        p1 = np.array([[0.02, 0.0], [0.01, 0.015]])
        mat = material_xs_from_raw(
            sig_s={0: [p0, p1]},
            cells_by_mat={0: (np.arange(4), np.zeros(4, dtype=int))},
            ng=2, nx=4,
        )
        S = ScatteringOperator(
            mat_xs=mat, quadrature=sn.quad, scattering_order=1,
            space=sn.full_field_space,  # S4-amendment: .frame demands the pose
        )
        op = BulkAnalysisOperator(S.frame, sn)
        stretched = _slab(width=2.0)
        with pytest.raises(ValueError, match="space-content"):
            op.apply(_composite(stretched))


class TestO13BoundaryOperator:
    def test_apply_refuses_content_differing_composite(self):
        sn = _slab()
        stretched = _slab(width=2.0)
        with pytest.raises(ValueError, match="space-content"):
            SNBoundaryOperator(sn).apply(_composite(stretched))


class TestF5ResidualPoseGuard:
    def test_bare_composite_on_a_carrying_system_refuses_by_arity(self):
        """R1's Mode-12(b) guard, re-keyed at CS4b S4 (the F5 ruling): the
        pose travels with the call as the POSED system, and a bare
        System-A state on a 2×2 system — which would silently drop System
        B's defect — is refused by the system's ARITY, never by a mesh
        read through the field."""
        sn = _sphere()
        assert sn.radial_characteristic_field_space is not None
        solver = SNSolver(sn, inner_solver="source_iteration")
        system = build_within_group_system(
            sn, solver.mat_xs, scattering_op=solver.scattering_op,
        )
        assert system.loss.n_cols == 2  # the carrying pose IS the arity
        psi = _composite(sn)
        q_ext = FullField(
            interior=AngularSourceSink.zeros_on(sn),
            boundary=AngularBoundarySourceSink.zeros_on(sn),
        )
        with pytest.raises(ValueError, match="starting-direction levels"):
            evaluate_residual(system, psi, q_ext)

    def test_a_bare_operator_cannot_make_the_full_system_claim(self):
        """The retired signature refuses loudly: a caller-supplied
        operator carries no pose, so the full-system claim is unspellable
        with one (arm-level equations use the module-private
        _typed_balance, which claims nothing)."""
        sn = _slab()
        psi = _composite(sn)
        with pytest.raises(TypeError, match="POSED system"):
            evaluate_residual(
                StreamingOperator(sn), psi, psi,  # type: ignore[arg-type]
            )
