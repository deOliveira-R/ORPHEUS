r"""Foundation tests for the :class:`HarmonicFrame` factory and its minted faces.

The F-1 contract (``frame_square_recarve.md``):

* the frame is the SHARED ``(basis, measure)`` operator factory —
  :meth:`from_galerkin` upgrades a generic ``frame_faces.angular_frame(quadrature, L)``
  sharing its basis + measure (table bit-identical), and identity is the
  table's identity (same-pairing frames compare equal);
* the mint verbs (:meth:`flux_analysis_on` / :meth:`source_reconstruction_on`)
  emit BOUND, carrier-typed operator faces — members of the numerics
  :class:`~orpheus.numerics.projection.AnalysisOperator` /
  :class:`~orpheus.numerics.projection.ReconstructionOperator` role family —
  taking at mint time the one input the frame cannot know (the angular field
  space) and deriving the moment codomain once (moment = f(angular, L),
  never the reverse), with the F-0 Parseval-dressed ``basis_space(frame)`` as
  the codomain's SH factor (single source);
* each face's ``apply`` is **bit-identical** to the generic ``np.ndarray``
  face it wraps — the typed seam adds carriers and binding, never a number;
* a wrong carrier raises ``TypeError``; a right-role carrier on a
  content-different space is REFUSED (the bound-operator admission, keeping
  the space-content vocabulary);
* a minted face's ``.H`` is the PHYSICAL Hilbert adjoint on the F-0 metrics —
  the S6-precursor pin ``M.H = R/W`` runs through the TRANSPORT face.

Foundation mark: these verify software invariants (the typed wrapper equals
the raw-face path; which object is bound where); the ``M.H = R/W`` closure is
the F-0 physics relation re-pinned through the transport seam.

Layer DAG Phase 4 (2026-11-01): this file holds the one row of tests/gates/transport/frames/test_harmonic_frame.py
that compares a face's moment codomain with the SN problem's own moment-field factory; the frame's laws stay in
transport.
"""
from __future__ import annotations


import pytest

from orpheus.geometry import StructuredGeometry
from orpheus.mesh import CellsByCount, Mesher
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.problem import SNProblem
from orpheus.transport.fields.harmonic_moment_flux import HarmonicMomentFlux
from orpheus.transport.frames import HarmonicFrame

from tests.gates.transport._materials import placeholder_materials
from orpheus.numerics import frame_faces
from orpheus.transport.boundary_conditions import WallCondition
from tests.gates._exterior import conditions_for

#: Exterior conditions of the tags this file's meshes carry (wave 2: a face's role is a deck or an int tag; the condition a tag
#: stands for is declared at the solve, keyword ``boundary_conditions``).
_CONDITIONS = {0: WallCondition.vacuum()}


pytestmark = pytest.mark.foundation

_L = 2


def _slab_mesh(nx: int = 4, ng: int = 2) -> SNProblem:
    mesh = Mesher(StructuredGeometry.from_homogeneous(1.0, 0)).partition(
        CellsByCount.uniform_width(nx),
    ).mesh
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    return SNProblem(mesh, quad, placeholder_materials(ng=ng), boundary_conditions=conditions_for(mesh, _CONDITIONS))


def _frame(m: SNProblem, L: int = _L) -> HarmonicFrame:
    return HarmonicFrame.from_galerkin(frame_faces.angular_frame(m.quad, L))


# ── from_galerkin: a behaviour-identical upgrade; the frame is SHARED ──


# ── the mint: bound faces, derived codomain, refusals ──────────────────


class TestMint:

    def test_moment_codomain_content_equals_the_carrier_mint(self) -> None:
        """moment = f(angular, L), derived ONCE at mint time — and it
        content-equals the carrier mint's own space (the admission seam is
        metric-blind, so the F-0 dressing does not disturb it).

        ⚠ Since #429 tracker 2.5 (2026-09-02) BOTH sides derive from ONE
        source — the quadrature's frame basis — so the equality below can
        no longer catch a drift between a hand-minted carrier space and
        the face's codomain (coding-standards: single-sourcing a duplicate
        demotes every gate that compared its copies). What it still tests
        is the DISCOVERY path: the carrier's mint and the face's mint reach
        the SAME frame object through the HUB, asserted by identity. The
        external pin the demotion needs — a hand-written literal shape —
        is the next test."""
        m = _slab_mesh()
        face = _frame(m).flux_analysis_on(m.angular_bulk_space)
        assert (
            face.codomain
            == HarmonicMomentFlux.zeros_for_problem_and_L(m, _L).space
        )
        assert face.frame is HarmonicFrame.for_space(m.angular_bulk_space, _L)


# ── the flux-analysis face (M ⊗ I) ─────────────────────────────────────


# ── the source-reconstruction face (R ⊗ I) ─────────────────────────────


# ── shareability: one factory, many faces, one table ───────────────────


# ── CS4c §14.4 — the HUB: interning + the blessed for_space chain ──────


# ── moved here by the layer DAG Phase 2 refresh (specs_2r §5) ──


