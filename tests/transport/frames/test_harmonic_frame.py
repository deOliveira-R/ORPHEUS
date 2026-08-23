r"""Foundation tests for the :class:`HarmonicFrame` factory and its minted faces.

The F-1 contract (``frame_square_recarve.md``):

* the frame is the SHARED ``(basis, measure)`` operator factory —
  :meth:`from_galerkin` upgrades a generic ``quadrature.angular_frame(L)``
  sharing its basis + measure (table bit-identical), and identity is the
  table's identity (same-pairing frames compare equal);
* the mint verbs (:meth:`flux_analysis_on` / :meth:`source_reconstruction_on`)
  emit BOUND, carrier-typed operator faces — members of the numerics
  :class:`~orpheus.numerics.projection.AnalysisOperator` /
  :class:`~orpheus.numerics.projection.ReconstructionOperator` role family —
  taking at mint time the one input the frame cannot know (the angular field
  space) and deriving the moment codomain once (moment = f(angular, L),
  never the reverse), with the F-0 Parseval-dressed ``frame.basis_space`` as
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
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D, Mesh2D
from orpheus.numerics.basis.indicator_basis import IndicatorBasis
from orpheus.numerics.frame import GalerkinFrame
from orpheus.numerics.measure import DiscreteMeasure
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.spaces.spherical_harmonic_space import SphericalHarmonicSpace
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.harmonic_moment_flux import HarmonicMomentFlux
from orpheus.transport.frames import (
    HarmonicAnalysisOperator,
    HarmonicFrame,
    HarmonicReconstructionOperator,
)
from orpheus.transport.source_sinks import (
    AngularSourceSink,
    HarmonicMomentSourceSink,
)

from tests.sn._test_helpers import placeholder_materials

pytestmark = pytest.mark.foundation

_L = 2


def _slab_mesh(nx: int = 4, ng: int = 2) -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _2d_mesh(nx: int = 3, ny: int = 3, ng: int = 1) -> SNMesh:
    mesh = Mesh2D(
        edges_x=np.linspace(0, 1, nx + 1),
        edges_y=np.linspace(0, 1, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
    )
    quad = Quadrature.level_symmetric(sn_order=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _frame(m: SNMesh, L: int = _L) -> HarmonicFrame:
    return HarmonicFrame.from_galerkin(m.quad.angular_frame(L))


def _angular_values(m: SNMesh, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.standard_normal((m.quad.N, m.ng, *m.spatial_shape))


def _moment_values(m: SNMesh, seed: int, L: int = _L) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.standard_normal((L + 1, 2 * L + 1, m.ng, *m.spatial_shape))


# ── from_galerkin: a behaviour-identical upgrade; the frame is SHARED ──


class TestFromGalerkin:
    def test_is_a_galerkin_frame_sharing_basis_measure(self) -> None:
        m = _slab_mesh()
        gf = m.quad.angular_frame(_L)
        frame = HarmonicFrame.from_galerkin(gf)
        assert isinstance(frame, type(gf))  # Liskov: HarmonicFrame IS-A GalerkinFrame
        assert frame.basis is gf.basis
        assert frame.measure is gf.measure
        # the trial table (Yₗᵐ tabulation) is bit-identical
        np.testing.assert_array_equal(frame.table, gf.table)

    def test_inherited_faces_unchanged(self) -> None:
        """The generic ndarray faces are inherited untouched (0-ULP-safe)."""
        m = _slab_mesh()
        gf = m.quad.angular_frame(_L)
        frame = HarmonicFrame.from_galerkin(gf)
        vals = _angular_values(m, 1)
        np.testing.assert_array_equal(
            frame.analysis.apply(vals), gf.analysis.apply(vals),
        )

    def test_non_sh_trial_basis_is_rejected_at_both_doors(self) -> None:
        """The SH narrowing is from_galerkin's ONLY remaining job (F-1) —
        and the direct constructor carries the same guard, so a harmonic
        frame over an indicator trial is unspellable, loudly, early."""
        indicator = IndicatorBasis((np.array([0.0, 1.0, 2.0]),))
        measure = DiscreteMeasure(
            nodes=np.array([0.5, 1.5]), weights=np.ones(2), support="spatial_R1",
        )
        with pytest.raises(TypeError, match="spherical-harmonic trial"):
            HarmonicFrame.from_galerkin(GalerkinFrame(indicator, measure))
        with pytest.raises(TypeError, match="spherical-harmonic trial"):
            HarmonicFrame(indicator, measure)  # type: ignore[arg-type]

    def test_same_pairing_frames_compare_equal(self) -> None:
        """The shareability half A1 destroyed: identity is the TABLE's
        identity again — two frames over the same (basis, measure) are the
        same projection, and no field-space binding distinguishes them."""
        m = _slab_mesh()
        gf = m.quad.angular_frame(_L)
        assert HarmonicFrame.from_galerkin(gf) == HarmonicFrame.from_galerkin(gf)


# ── the mint: bound faces, derived codomain, refusals ──────────────────


class TestMint:
    def test_flux_analysis_face_is_a_bound_role_member(self) -> None:
        m = _slab_mesh()
        frame = _frame(m)
        face = frame.flux_analysis_on(m.angular_bulk_space)
        from orpheus.numerics.projection import AnalysisOperator

        assert isinstance(face, HarmonicAnalysisOperator)
        assert isinstance(face, AnalysisOperator)
        assert face.domain is m.angular_bulk_space   # bound instance, not a copy
        assert face.frame is frame                   # the shared factory

    def test_moment_codomain_content_equals_the_carrier_mint(self) -> None:
        """moment = f(angular, L), derived ONCE at mint time — and it
        content-equals the carrier mint's own space (the admission seam is
        metric-blind, so the F-0 dressing does not disturb it)."""
        m = _slab_mesh()
        face = _frame(m).flux_analysis_on(m.angular_bulk_space)
        assert (
            face.codomain
            == HarmonicMomentFlux.zeros_for_mesh_and_L(m, _L).space
        )

    def test_codomain_sh_factor_is_the_frames_parseval_space(self) -> None:
        """The single-source pin: the minted codomain's SH factor IS the
        frame's own F-0-dressed ``basis_space`` object — the Parseval metric
        rides into the product with no re-derivation. (LS4 measures a
        DIAGONAL discrete Gram, so the factor is genuinely dressed here;
        on the slab's DENSE verdict the same read returns the undressed
        continuum space — honest either way, one source both ways.)"""
        m = _2d_mesh()
        frame = _frame(m)
        face = frame.flux_analysis_on(m.angular_bulk_space)
        from orpheus.numerics.space import TensorProductSpace

        assert isinstance(face.codomain, TensorProductSpace)
        sh = face.codomain.find_factor(SphericalHarmonicSpace)
        assert sh is frame.basis_space

    def test_mint_refuses_a_shape_only_angular_space(self) -> None:
        from orpheus.numerics.space import FunctionSpace

        m = _slab_mesh()
        shape_only = FunctionSpace("bare", (m.quad.N, m.ng, *m.spatial_shape))
        with pytest.raises(TypeError, match="axis-built"):
            _frame(m).flux_analysis_on(shape_only)

    def test_widened_angular_space_mints_widened_faces(self) -> None:
        """A spatial-moment-widened angular domain (the LD iterate width)
        derives a moment codomain carrying the densified moment factor —
        read off the SPACE at mint time, not off any operand — and the two
        minted faces close the square on the SAME derived moment space."""
        from orpheus.transport.spatial import LinearDiscontinuous

        mesh = Mesh1D(
            edges=np.linspace(0.0, 1.0, 5),
            mat_ids=np.zeros(4, dtype=int),
            coord=CoordSystem.CARTESIAN,
            bc_left=BC("vacuum"),
            bc_right=BC("vacuum"),
        )
        m = SNMesh(
            mesh,
            Quadrature.gauss_legendre(n_ordinates=4),
            placeholder_materials(ng=2),
            scheme=LinearDiscontinuous(),
        )
        widened = AngularFlux.zeros_on(m, spatial_moments=2)
        frame = _frame(m)
        analysis = frame.flux_analysis_on(widened.space)
        assert analysis.domain is widened.space
        moments = analysis.apply(widened)
        assert isinstance(moments, HarmonicMomentFlux)
        assert moments.space is analysis.codomain
        assert moments.spatial_moments_per_axis == 2
        # the square's other minted edge derives the SAME moment space
        recon = frame.source_reconstruction_on(widened.space)
        assert recon.domain == analysis.codomain
        assert recon.codomain is widened.space


# ── the flux-analysis face (M ⊗ I) ─────────────────────────────────────


class TestFluxAnalysisFace:
    def test_apply_bit_identical_and_rides_the_bound_codomain(self) -> None:
        m = _slab_mesh()
        frame = _frame(m)
        face = frame.flux_analysis_on(m.angular_bulk_space)
        vals = _angular_values(m, 2)
        psi = AngularFlux.from_mesh(vals, m)
        moments = face.apply(psi)
        assert isinstance(moments, HarmonicMomentFlux)
        assert moments.space is face.codomain    # bound codomain, no re-mint
        assert moments.L == _L
        np.testing.assert_array_equal(
            moments.values, frame.analysis.apply(vals),
        )

    def test_wrong_carrier_raises(self) -> None:
        m = _slab_mesh()
        face = _frame(m).flux_analysis_on(m.angular_bulk_space)
        moment = HarmonicMomentFlux.from_mesh_and_L(_moment_values(m, 4), m, _L)
        with pytest.raises(TypeError, match="unsupported carrier"):
            face.apply(moment)  # type: ignore[arg-type]

    def test_content_different_operand_is_refused(self) -> None:
        """The bound-operator admission (the space-content invariant): a
        right-role carrier riding a content-DIFFERENT space is refused,
        naming both spaces."""
        m = _slab_mesh(nx=4)
        other = _slab_mesh(nx=5)  # different spatial axis → different space
        face = _frame(m).flux_analysis_on(m.angular_bulk_space)
        psi = AngularFlux.from_mesh(_angular_values(other, 11), other)
        with pytest.raises(TypeError, match="bound to angular domain"):
            face.apply(psi)

    def test_transpose_round_trips_the_carriers_and_accepts_values(self) -> None:
        """``apply_transpose`` is the representation transpose: codomain
        carrier → domain carrier, bit-identical to the numerics face; the
        raw-``ndarray`` seam (what ``_AdjointOperator`` drives) returns raw
        values."""
        m = _slab_mesh()
        frame = _frame(m)
        face = frame.flux_analysis_on(m.angular_bulk_space)
        moment = HarmonicMomentFlux(
            values=_moment_values(m, 5), space=face.codomain, L=_L,
            spatial_moments=1,
        )
        back = face.apply_transpose(moment)
        assert isinstance(back, AngularFlux)
        assert back.space is face.domain
        np.testing.assert_array_equal(
            back.values, frame.analysis.apply_transpose(moment.values),
        )
        np.testing.assert_array_equal(
            face.apply_transpose(moment.values),
            frame.analysis.apply_transpose(moment.values),
        )

    def test_s6_precursor_the_faces_H_is_the_physical_adjoint(self) -> None:
        r"""``face.H = R/W`` THROUGH the transport face — the S6-precursor pin.

        The minted face's ``.H`` routes the generic ``_AdjointOperator``
        sandwich through the BOUND product spaces: the F-0 Parseval SH factor
        on the moment side, the per-axis (quadrature ⊗ cell) measures on the
        angular side. The cell metrics cancel interior (the composite law)
        and the SH closure collapses to the one scalar ``W`` — so the
        transport face's adjoint equals ``frame.reconstruction/W`` exactly
        as the numerics face's does. Runs on the LS4 2-D mesh (DIAGONAL
        verdict — the dressing is live; the slab's DENSE arm has no such
        claim and is pinned separately in ``test_frame.py``).
        """
        m = _2d_mesh()
        frame = _frame(m)
        face = frame.flux_analysis_on(m.angular_bulk_space)
        W = float(frame.measure.weights.sum())
        y = _moment_values(m, 6)
        np.testing.assert_allclose(
            face.H.apply(y),
            frame.reconstruction.apply(y) / W,
            rtol=1e-12, atol=1e-13,
        )


# ── the source-reconstruction face (R ⊗ I) ─────────────────────────────


class TestSourceReconstructionFace:
    def test_apply_bit_identical_and_rides_the_bound_codomain(self) -> None:
        m = _slab_mesh()
        frame = _frame(m)
        face = frame.source_reconstruction_on(m.angular_bulk_space)
        vals = _moment_values(m, 6)
        q = HarmonicMomentSourceSink(
            values=vals, space=face.domain, L=_L, spatial_moments=1,
        )
        out = face.apply(q)
        assert isinstance(out, AngularSourceSink)
        assert out.space is m.angular_bulk_space
        np.testing.assert_array_equal(
            out.values, frame.reconstruction.apply(vals),
        )

    def test_wrong_carrier_raises(self) -> None:
        m = _slab_mesh()
        face = _frame(m).source_reconstruction_on(m.angular_bulk_space)
        psi = AngularFlux.from_mesh(_angular_values(m, 7), m)
        with pytest.raises(TypeError, match="unsupported carrier"):
            face.apply(psi)  # type: ignore[arg-type]

    def test_content_different_operand_is_refused(self) -> None:
        m = _slab_mesh(nx=4)
        other = _slab_mesh(nx=5)
        face = _frame(m).source_reconstruction_on(m.angular_bulk_space)
        q = HarmonicMomentSourceSink.from_mesh_and_L(
            _moment_values(other, 12), other, _L,
        )
        with pytest.raises(TypeError, match="bound to moment domain"):
            face.apply(q)

    def test_is_adjointable_and_H_is_W_times_M(self) -> None:
        r"""The face ADDS the transpose over the apply-only role, so ``R.H``
        is free — and on the F-0 metrics it is the physical ``W·M``
        (the mirror of the analysis face's S6-precursor pin)."""
        m = _2d_mesh()
        frame = _frame(m)
        face = frame.source_reconstruction_on(m.angular_bulk_space)
        assert face.is_adjointable
        W = float(frame.measure.weights.sum())
        v = _angular_values(m, 8)
        np.testing.assert_allclose(
            face.H.apply(v),
            W * frame.analysis.apply(v),
            rtol=1e-12, atol=1e-13,
        )

    def test_2d_typed_round(self) -> None:
        """2-D smoke on the minted pair: analyse a flux, synthesise a source
        — the two faces the production windowed path composes."""
        m = _2d_mesh()
        frame = _frame(m)
        analysis = frame.flux_analysis_on(m.angular_bulk_space)
        recon = frame.source_reconstruction_on(m.angular_bulk_space)
        psi = AngularFlux.from_mesh(_angular_values(m, 10), m)
        moments = analysis.apply(psi)
        assert moments.values.shape == (_L + 1, 2 * _L + 1, m.ng, *m.spatial_shape)
        q = HarmonicMomentSourceSink(
            values=moments.values, space=recon.domain, L=_L, spatial_moments=1,
        )
        assert isinstance(recon.apply(q), AngularSourceSink)


# ── shareability: one factory, many faces, one table ───────────────────


class TestShareability:
    def test_two_mints_share_the_frame_and_its_table(self) -> None:
        """The 'derived, not independent' guarantee at the array level: every
        face minted from one frame reads the SAME cached table (the F-1
        repair of A1's unshareable factory)."""
        m = _slab_mesh()
        frame = _frame(m)
        a = frame.flux_analysis_on(m.angular_bulk_space)
        r = frame.source_reconstruction_on(m.angular_bulk_space)
        assert a.frame is r.frame is frame
        assert a.frame.table is r.frame.table
