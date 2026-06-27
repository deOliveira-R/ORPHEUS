r"""Wave O (Issue #208) step O.4a.2 Commit 1 — the whole-trace boundary
operator ``B`` (:class:`~orpheus.sn.operators.boundary.SNBoundaryOperator`).

``B`` is the ``A_ss`` block of the SN algebra ``(L_full + C − S − F − B)``: a
BOUNDARY-block leaf on the :class:`TimedFullField` carrier that applies each
true boundary face's realized law to that face's trace slot, with zero bulk
action. These foundation tests pin the assembly BEFORE anything consumes ``B``
(the ``−B`` wiring + the bare-``L_full`` flip is O.4a.2 Commit 2):

* the role / domain / capabilities contract;
* zero bulk action;
* **per-face wiring** — ``B`` applies the RIGHT face's law to the RIGHT slot,
  emitting on the **inflow row only** (``B`` is the ``A_ss`` block
  ``V_outflow → V_inflow``; the discriminating case uses asymmetric BCs so a
  face↔face swap is caught);
* **block-diagonal over faces** — a single-face perturbation stays on that face;
* the ``apply_transpose`` capability intersection (advertised iff every face law
  honours it — white would drop it; see the stub negative).
"""
from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.numerics.operator import (
    CAP_APPLY,
    CAP_APPLY_TRANSPOSE,
    BlockRole,
    BoundaryOperator,
    BulkOperator,
    FullOperator,
)
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.operators.boundary import SNBoundaryOperator
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.boundary_flux import BoundaryFlux
from orpheus.transport.source_sinks import BoundarySourceSink
from orpheus.transport.timed_full_field import TimedFullField
from tests.sn._test_helpers import placeholder_materials

pytestmark = [pytest.mark.foundation]


def _sn(geometry: str, bcs: tuple, nx: int = 4, ng: int = 1) -> SNMesh:
    geom = StructuredGeometry(
        geometry=geometry,
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=bcs,
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=nx),))
    # Cylinder's angular redistribution needs a level-structured quadrature;
    # slab / sphere accept the 1-D Gauss–Legendre set.
    quad = (
        Quadrature.product(n_mu=2, n_phi=4)
        if geometry == "CYL"
        else Quadrature.gauss_legendre(n_ordinates=4)
    )
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _random_state(sn: SNMesh, seed: int = 7) -> TimedFullField:
    rng = np.random.default_rng(seed)
    z = TimedFullField.zeros(bulk=AngularFlux, boundary=BoundaryFlux, mesh=sn)
    return replace(
        z,
        bulk=replace(z.bulk, values=rng.uniform(0.5, 2.0, size=z.bulk.values.shape)),
        boundary=replace(
            z.boundary, values=rng.uniform(0.5, 2.0, size=z.boundary.values.shape),
        ),
    )


# Geometry × BC cases reachable through SNMesh (1-D faces support only
# reflective / vacuum). Slab uses ASYMMETRIC BCs so the per-face wiring test
# discriminates a face↔face swap.
_CASES = {
    "slab_vacuum_reflective": ("SLB", (BC.vacuum, BC.reflective)),
    "slab_reflective_reflective": ("SLB", (BC.reflective, BC.reflective)),
    "sphere_reflective": ("SPH", (BC.reflective,)),
    "cyl_reflective": ("CYL", (BC.reflective,)),
}


class TestContract:
    def test_block_role_is_boundary_and_exclusive(self) -> None:
        B = SNBoundaryOperator(_sn("SLB", (BC.vacuum, BC.reflective)))
        assert B.block_role is BlockRole.BOUNDARY
        assert isinstance(B, BoundaryOperator)
        assert not isinstance(B, BulkOperator)
        assert not isinstance(B, FullOperator)

    def test_domain_and_codomain_are_the_full_field_space(self) -> None:
        # Wave O / O.2b R5: ``B`` is an endomorphism on the composite
        # carrier (bulk ⊕ trace) — ``B.apply`` consumes / emits a full
        # ``TimedFullField`` (zero bulk + reflected trace), so it advertises
        # ``sn.full_field_space``, the SAME space L/C/S/F report (so the
        # OperatorSum guard accepts ``L + C - S - F - B``). The trace metric
        # ``B.H`` reads lives on the composite's trace block.
        sn = _sn("SLB", (BC.vacuum, BC.reflective))
        B = SNBoundaryOperator(sn)
        assert B.domain is sn.full_field_space
        assert B.codomain is sn.full_field_space
        # the composite trace block IS the mesh trace space (block identity)
        assert B.domain.trace_space is sn.trace

    def test_apply_advertised(self) -> None:
        B = SNBoundaryOperator(_sn("SLB", (BC.vacuum, BC.reflective)))
        assert CAP_APPLY in B.capabilities


class TestApply:
    @pytest.mark.parametrize("case_id", list(_CASES))
    def test_apply_is_zero_bulk(self, case_id: str) -> None:
        """``B`` is ``A_ss`` only — no bulk action."""
        sn = _sn(*_CASES[case_id])
        out = SNBoundaryOperator(sn).apply(_random_state(sn))
        assert not out.bulk.values.any()
        # B.5.2: B.apply emits B·ψ.outflow — the operator output is Aψ (a
        # source/sink), NOT a residual.  Its boundary is the source/sink role
        # leaf (mirrors the bulk's AngularSourceSink); the residual only arises
        # from from_balance, never straight off the operator output.
        assert isinstance(out.boundary, BoundarySourceSink)

    @pytest.mark.parametrize("case_id", list(_CASES))
    def test_apply_per_face_equals_legacy_bc_apply_on_inflow_row(
        self, case_id: str,
    ) -> None:
        """Per-face wiring: ``B`` applies the face's OWN law to the face's slot,
        emitting on the **inflow row only**.

        ``B`` is the ``A_ss`` block ``V_outflow → V_inflow`` — its codomain is
        the inflow trace, so ``B.apply`` must agree with the realized per-face
        ``bc.apply`` on the inflow ordinate slots and be ZERO on the outflow
        slots (the realized law is a full-face operator whose specular
        permutation also maps the input inflow onto the output outflow slots;
        that ``R·ψ.inflow`` term is spurious for the block and would corrupt the
        outflow-definition residual of ``−B`` — see
        :meth:`SNBoundaryOperator._apply_faces`).

        Asymmetric-BC slab makes a face↔face swap observable (vacuum xmin vs
        reflective xmax produce different inflow outputs) — vv L11 catches
        Failure mode #5 (index/face swap)."""
        sn = _sn(*_CASES[case_id])
        psi = _random_state(sn)
        out = SNBoundaryOperator(sn).apply(psi)
        for face in sn.trace.layout.faces:
            bc = sn.bc[face]
            inflow = sn.trace.inflow_indices_for_face(face)
            outflow = sn.trace.outflow_indices_for_face(face)
            got = out.boundary.face_view(face)
            expected_full = bc.apply(psi.boundary.face_view(face))
            # Inflow row: agrees with the realized law (the consistency action).
            np.testing.assert_array_equal(got[inflow], expected_full[inflow])
            # Outflow row: zero — ``B`` contributes nothing to the outflow
            # definition residual (the row carries only ``I·ψ.outflow − streamed``).
            assert not got[outflow].any(), (
                f"{case_id} face {face!r}: B emitted non-zero on the outflow "
                f"row — it is not a clean A_ss (V_outflow → V_inflow) block."
            )

    def test_block_diagonal_no_face_mixing(self) -> None:
        """A perturbation on ONE face's input slot affects ONLY that face's
        output (``B`` is block-diagonal over faces — it never mixes faces)."""
        sn = _sn("SLB", (BC.reflective, BC.reflective))
        B = SNBoundaryOperator(sn)
        psi = _random_state(sn, seed=1)
        other = _random_state(sn, seed=2)
        # psi3 = psi with ONLY the xmin slot replaced (xmax slot identical to psi).
        b3 = replace(psi.boundary, values=psi.boundary.values.copy())
        psi3 = replace(psi, boundary=b3)
        psi3.boundary.face_view("xmin")[:] = other.boundary.face_view("xmin")

        out = B.apply(psi)
        out3 = B.apply(psi3)
        # xmax output unchanged — it depends only on the (identical) xmax input.
        np.testing.assert_array_equal(
            out.boundary.face_view("xmax"),
            out3.boundary.face_view("xmax"),
        )
        # Sanity: the xmin perturbation actually changed the xmin output (else
        # the block-diagonal claim would be vacuous).
        assert not np.array_equal(
            out.boundary.face_view("xmin"),
            out3.boundary.face_view("xmin"),
        )


class TestApplyTransposeCapability:
    @pytest.mark.parametrize("case_id", list(_CASES))
    def test_capabilities_include_apply_transpose_when_all_faces_support(
        self, case_id: str,
    ) -> None:
        """Reflective / vacuum faces all advertise apply_transpose, so ``B`` does."""
        sn = _sn(*_CASES[case_id])
        B = SNBoundaryOperator(sn)
        assert CAP_APPLY_TRANSPOSE in B.capabilities
        # The Euclidean transpose ``Bᵀ: V_inflow → V_outflow`` emits on the
        # OUTFLOW row only (the transpose of the inflow-row projection), and
        # agrees there with the per-face ``bc.apply_transpose``.
        psi = _random_state(sn)
        out = B.apply_transpose(psi)
        for face in sn.trace.layout.faces:
            bc = sn.bc[face]
            inflow = sn.trace.inflow_indices_for_face(face)
            outflow = sn.trace.outflow_indices_for_face(face)
            got = out.boundary.face_view(face)
            expected_full = bc.apply_transpose(psi.boundary.face_view(face))
            np.testing.assert_array_equal(got[outflow], expected_full[outflow])
            assert not got[inflow].any(), (
                f"{case_id} face {face!r}: Bᵀ emitted non-zero on the inflow "
                f"row — the transpose of an A_ss (V_outflow → V_inflow) block "
                f"must land on the outflow row."
            )

    def test_capabilities_drop_apply_transpose_when_a_face_lacks_it(self) -> None:
        """The capability is an INTERSECTION — if any face law cannot transpose
        (e.g. the white BC, self-adjoint only under the |Ω·n|·w metric), ``B``
        must NOT advertise apply_transpose (vv L11 negative; prevents a silent
        wrong/raising adjoint in a Krylov adjoint solve)."""

        class _NoTransposeLaw:
            capabilities = frozenset({CAP_APPLY})

            def apply(self, x):  # noqa: D401 - stub
                return x

        sn = _sn("SLB", (BC.vacuum, BC.reflective))

        class _BWithStubFace(SNBoundaryOperator):
            @property
            def _face_laws(self):
                laws = dict(super()._face_laws)
                laws[next(iter(laws))] = _NoTransposeLaw()
                return laws

        B = _BWithStubFace(sn)
        assert CAP_APPLY in B.capabilities
        assert CAP_APPLY_TRANSPOSE not in B.capabilities


class TestFaceRestrictedReflect:
    """Phase 3 sub-step 2: ``reflect_into_inflow(faces=...)`` restricts the
    trace reflection to a face subset — the octant-group Gauss-Seidel schedule
    reflects ONLY the just-swept group's reflective outgoing faces between
    octant-group sweeps (the ``(L+C−B_lower)⁻¹`` forward substitution).

    ``B`` is block-diagonal over faces, so the subset MUST be the EXACT
    restriction: the per-face inflow output is identical whether reflected
    alone or as part of the whole trace, and the per-face reflects PARTITION
    the whole-trace reflect (no coupling dropped, no term double-counted).
    """

    def test_subset_reflects_only_selected_faces(self) -> None:
        """``faces=['xmax']`` emits reflected inflow on xmax and leaves the
        unselected xmin face untouched (zero)."""
        sn = _sn("SLB", (BC.reflective, BC.reflective))
        B = SNBoundaryOperator(sn)
        boundary = _random_state(sn).boundary
        only_xmax = B.reflect_into_inflow(boundary, faces=["xmax"])
        full = B.reflect_into_inflow(boundary)  # faces=None -> whole trace
        xmax_inflow = sn.trace.inflow_indices_for_face("xmax")
        # Selected face: bit-identical to the whole-trace reflect (the exact
        # restriction — face-block-diagonal, so xmax's output is independent
        # of whether xmin was also reflected).
        np.testing.assert_array_equal(
            only_xmax.face_view("xmax"), full.face_view("xmax"),
        )
        # Unselected face: untouched -> stays zero.
        assert not only_xmax.face_view("xmin").any(), (
            "reflect_into_inflow(faces=['xmax']) emitted on the unselected "
            "xmin face — the restriction is not clean."
        )
        # Sanity: the selected face actually carries non-zero reflected inflow
        # (else the restriction claim would be vacuous).
        assert only_xmax.face_view("xmax")[xmax_inflow].any()

    def test_subset_partitions_the_whole_trace_reflect(self) -> None:
        """Single-face reflects sum EXACTLY to the whole-trace reflect — ``B``
        never couples faces, so the per-face restrictions are a clean partition
        (vv L11: catches a face↔face coupling leak that the
        reflect-only-selected test alone would miss)."""
        sn = _sn("SLB", (BC.reflective, BC.reflective))
        B = SNBoundaryOperator(sn)
        boundary = _random_state(sn).boundary
        full = B.reflect_into_inflow(boundary)
        xmin_only = B.reflect_into_inflow(boundary, faces=["xmin"])
        xmax_only = B.reflect_into_inflow(boundary, faces=["xmax"])
        np.testing.assert_array_equal(
            full.values, xmin_only.values + xmax_only.values,
        )

    def test_faces_none_equals_explicit_all_faces(self) -> None:
        """The default (``faces=None``) is the whole-trace reflect — identical
        to passing every boundary face explicitly."""
        sn = _sn("SLB", (BC.reflective, BC.reflective))
        B = SNBoundaryOperator(sn)
        boundary = _random_state(sn).boundary
        all_faces = list(sn.trace.layout.faces)
        np.testing.assert_array_equal(
            B.reflect_into_inflow(boundary).values,
            B.reflect_into_inflow(boundary, faces=all_faces).values,
        )

    def test_unknown_face_raises(self) -> None:
        """A face not on the mesh boundary is a caller error — raise, do not
        silently skip (illegal states unrepresentable)."""
        sn = _sn("SLB", (BC.reflective, BC.reflective))
        B = SNBoundaryOperator(sn)
        boundary = _random_state(sn).boundary
        with pytest.raises(ValueError, match="boundary faces"):
            B.reflect_into_inflow(boundary, faces=["bogus_face"])
