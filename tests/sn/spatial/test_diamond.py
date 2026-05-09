"""Tests for the :class:`DiamondDifference` cell-update strategy.

Round 2 of Wave C of the SN reshape campaign.  These tests pin the
**bit-identical** contract between :class:`DiamondDifference` and
the existing inlined sweep math at :mod:`orpheus.sn.sweep` — they
assert ``np.array_equal`` (not ``np.allclose``) against per-cell
scalar formulas that mirror the sweep's operation order verbatim.

The tests are tagged ``@pytest.mark.foundation`` because:

* Diamond Difference's L1 transport accuracy is verified
  transitively via the existing sweep MMS suite at
  ``tests/sn/l1_analytical/`` and the bit-identical regression
  snapshot tests at ``tests/sn/regression/`` — Round 2 does not
  modify the sweep, so those gates remain green by construction.
* These tests verify the **software-contract** invariants of the
  strategy itself: that its per-branch algebra reproduces the
  sweep's per-cell algebra to the bit, that traits are class-level
  constants, and that the cylindrical-degenerate branch correctly
  signals "no spatial face flow" via
  ``CellResult.outgoing_spatial_flux=None``.

Test classes
============

* :class:`TestTraits` — class-level ``is_linear`` /
  ``is_positivity_preserving`` attributes have the expected
  values.
* :class:`TestBitIdenticalSlab` — DD's slab branch reproduces
  :func:`orpheus.sn.sweep._sweep_1d_cumprod` /
  :func:`orpheus.sn.sweep._solve_recurrence` lines 117-123 +
  208-222 bit-for-bit on synthetic per-cell inputs.
* :class:`TestBitIdenticalCurvilinear` — DD's curvilinear branch
  reproduces :func:`orpheus.sn.sweep._sweep_1d_spherical` lines
  350-361 bit-for-bit.
* :class:`TestCylindricalDegenerate` — DD's degenerate branch
  reproduces :func:`orpheus.sn.sweep._sweep_1d_cylindrical` lines
  533-546 bit-for-bit and signals via
  ``outgoing_spatial_flux=None``.
* :class:`TestPositivityFailure` — DD can produce negative
  outgoing fluxes from positive inputs in thin / large-source
  cells (Lewis & Miller §5.3 counter-example) — gates the
  ``is_positivity_preserving = False`` trait against the actual
  numerical behaviour.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import (
    BC,
    CoordSystem,
    Mesh1D,
    cylindrical_streaming,
    slab_streaming,
    spherical_streaming,
)
from orpheus.geometry.reduced_operator import StreamingTerms
from orpheus.sn.quadrature import GaussLegendre1D, ProductQuadrature
from orpheus.sn.spatial import DiamondDifference, UpstreamState
from orpheus.sn.spatial.cell_update import CellVisit


# ═══════════════════════════════════════════════════════════════════════
# Mesh fixtures
# ═══════════════════════════════════════════════════════════════════════

def _slab_mesh(nx: int = 5, length: float = 1.0) -> Mesh1D:
    """Slab mesh with vacuum BCs and uniform spacing."""
    return Mesh1D(
        edges=np.linspace(0.0, length, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )


def _spherical_mesh(nx: int = 5, radius: float = 1.0) -> Mesh1D:
    """Spherical mesh with reflective inner / vacuum outer BCs."""
    return Mesh1D(
        edges=np.linspace(0.0, radius, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )


def _cylindrical_mesh(nx: int = 5, radius: float = 1.0) -> Mesh1D:
    """Cylindrical mesh with reflective inner / vacuum outer BCs."""
    return Mesh1D(
        edges=np.linspace(0.0, radius, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )


# ═══════════════════════════════════════════════════════════════════════
# TestTraits
# ═══════════════════════════════════════════════════════════════════════

class TestTraits:
    """Class-level traits of :class:`DiamondDifference`.

    Traits are declared as ``ClassVar`` so they live on the class,
    not the instance — verified here by class-level access.  Lewis &
    Miller §5.3 motivates the ``is_positivity_preserving = False``
    value: DD's WDD spatial closure can produce negative outgoing
    flux from positive inputs in thin cells with large source.
    """

    @pytest.mark.foundation
    def test_is_linear_true(self):
        assert DiamondDifference.is_linear is True

    @pytest.mark.foundation
    def test_is_positivity_preserving_false(self):
        assert DiamondDifference.is_positivity_preserving is False

    @pytest.mark.foundation
    def test_traits_accessible_on_instance(self):
        s = DiamondDifference()
        assert s.is_linear is True
        assert s.is_positivity_preserving is False


# ═══════════════════════════════════════════════════════════════════════
# TestBitIdenticalSlab
# ═══════════════════════════════════════════════════════════════════════

class TestBitIdenticalSlab:
    """Slab DD reproduces :func:`orpheus.sn.sweep._solve_recurrence`
    bit-for-bit.

    The sweep solves the slab DD recurrence over a whole row of
    cells via cumulative products; here we mirror the **per-cell
    scalar form** of that recurrence and assert ``np.array_equal``
    against :class:`DiamondDifference`'s slab branch.

    The reference per-cell scalar form below mirrors sweep.py:117-123
    (where ``denom``, ``stream_coeff``, ``source_coeff``, ``bQ`` are
    set up vectorised) followed by sweep.py:208-222 (where
    ``_solve_recurrence`` applies them per cell).  At ``cell_idx=0``
    the recurrence reduces to ``psi_in[0] = psi0``, ``psi_out[0] =
    a[0] * psi_in[0] + s[0]``, ``psi_avg[0] = 0.5 * (psi_in[0] +
    psi_out[0])``.  Operation order is preserved exactly.
    """

    @pytest.mark.foundation
    @pytest.mark.verifies("dd-slab-scalar")
    def test_slab_first_cell_bit_identical(self):
        """First-cell DD update reproduces the sweep's scalar form bit-for-bit.

        Mirrors ``_solve_recurrence`` at ``i = 0`` (the entry-cell of
        a forward sweep), which is the cleanest discriminator: at
        ``i = 0`` the cumulative-product chain has not yet
        accumulated, so the per-cell scalar form is directly visible
        without unrolling.
        """
        mesh = _slab_mesh(nx=5, length=1.0)
        quad = GaussLegendre1D.create(4)
        op = slab_streaming(mesh, quad)

        # Pull the per-cell streaming terms for cell 0, the most
        # positive ordinate (forward-sweep entry).
        ng = 2
        cell_idx = 0
        # Most-positive ordinate is the last index for GL quadrature.
        direction_idx = quad.N - 1
        st = op.streaming_terms(cell_idx, direction_idx)

        # Synthetic ng-shaped inputs.
        total_xs = np.array([1.5, 0.7])
        # ``source`` is per-contract Q · V · weight_norm — for slab
        # V == chord_length, so we pre-multiply by chord_length and
        # weight_norm here as the sweep would.
        weight_norm = 1.0 / quad.weights.sum()
        Q = np.array([2.5, 0.4])
        source = Q * st.chord_length * weight_norm
        psi_in = np.array([0.3, 0.1])
        upstream = UpstreamState(
            spatial_upstream=psi_in,
            angular_upstream=None,
        )

        # Reference scalar form — mirrors sweep.py:119-123 + 221-222.
        chord = st.chord_length
        abs_mu = st.abs_mu
        ref_denom = 2.0 * abs_mu + chord * total_xs
        ref_a = (2.0 * abs_mu - chord * total_xs) / ref_denom
        ref_s = 2.0 * source / ref_denom
        ref_psi_out = ref_a * psi_in + ref_s
        ref_psi_avg = 0.5 * (psi_in + ref_psi_out)

        # Slab visit: face_area_downstream is None — slab DD does not
        # read face areas (the recurrence is built around chord /
        # |μ| only).
        visit = CellVisit(cell_idx=cell_idx, streaming_terms=st)
        strat = DiamondDifference()
        result = strat.update(visit, total_xs, source, upstream)

        # Bit-identical: np.array_equal, not np.allclose.
        assert np.array_equal(result.cell_average_flux, ref_psi_avg)
        assert np.array_equal(result.outgoing_spatial_flux, ref_psi_out)
        # Slab has no angular redistribution.
        assert result.outgoing_angular_state is None

    @pytest.mark.foundation
    @pytest.mark.verifies("dd-slab-scalar")
    def test_slab_interior_cell_bit_identical(self):
        """Interior-cell DD update reproduces the sweep's scalar form
        bit-for-bit.

        Mirrors ``_solve_recurrence`` at ``i = 2`` of a 5-cell mesh —
        a non-degenerate interior position with a non-trivial
        ``psi_in`` that is the would-be cumulative result of upstream
        cells.  The per-cell algebra is the same as the entry-cell
        case; this test guards against accidental special-casing of
        ``cell_idx == 0``.
        """
        mesh = _slab_mesh(nx=5, length=1.0)
        quad = GaussLegendre1D.create(4)
        op = slab_streaming(mesh, quad)

        ng = 3
        cell_idx = 2
        direction_idx = quad.N - 2  # second-most-positive ordinate
        st = op.streaming_terms(cell_idx, direction_idx)

        total_xs = np.array([2.0, 0.5, 1.1])
        weight_norm = 1.0 / quad.weights.sum()
        Q = np.array([0.9, 1.7, 0.2])
        source = Q * st.chord_length * weight_norm
        psi_in = np.array([0.42, 0.31, 0.55])  # synthetic upstream
        upstream = UpstreamState(
            spatial_upstream=psi_in,
            angular_upstream=None,
        )

        chord = st.chord_length
        abs_mu = st.abs_mu
        ref_denom = 2.0 * abs_mu + chord * total_xs
        ref_a = (2.0 * abs_mu - chord * total_xs) / ref_denom
        ref_s = 2.0 * source / ref_denom
        ref_psi_out = ref_a * psi_in + ref_s
        ref_psi_avg = 0.5 * (psi_in + ref_psi_out)

        visit = CellVisit(cell_idx=cell_idx, streaming_terms=st)
        strat = DiamondDifference()
        result = strat.update(visit, total_xs, source, upstream)

        assert np.array_equal(result.cell_average_flux, ref_psi_avg)
        assert np.array_equal(result.outgoing_spatial_flux, ref_psi_out)
        assert result.outgoing_angular_state is None

    @pytest.mark.foundation
    @pytest.mark.verifies("dd-slab-scalar")
    def test_slab_negative_ordinate_bit_identical(self):
        """Backward-sweep ordinate (negative μ) — abs_mu carries the
        magnitude.

        The sweep's slab cumprod path takes ``abs(quad.mu_x)``; the
        per-cell strategy reads ``streaming_terms.abs_mu`` directly.
        Both paths use the same magnitude, so bit-equality holds.
        """
        mesh = _slab_mesh(nx=5, length=1.0)
        quad = GaussLegendre1D.create(4)
        op = slab_streaming(mesh, quad)

        ng = 2
        cell_idx = 4  # last cell, where backward sweep enters
        direction_idx = 0  # most-negative ordinate
        st = op.streaming_terms(cell_idx, direction_idx)

        # Verify we got a negative-μ ordinate so this test isn't
        # just a re-run of the positive case.
        assert st.mu is not None and st.mu < 0
        assert st.abs_mu == abs(st.mu)

        total_xs = np.array([1.2, 0.8])
        weight_norm = 1.0 / quad.weights.sum()
        Q = np.array([1.5, 0.6])
        source = Q * st.chord_length * weight_norm
        psi_in = np.array([0.0, 0.0])  # vacuum BC at outer face
        upstream = UpstreamState(
            spatial_upstream=psi_in,
            angular_upstream=None,
        )

        chord = st.chord_length
        abs_mu = st.abs_mu
        ref_denom = 2.0 * abs_mu + chord * total_xs
        ref_a = (2.0 * abs_mu - chord * total_xs) / ref_denom
        ref_s = 2.0 * source / ref_denom
        ref_psi_out = ref_a * psi_in + ref_s
        ref_psi_avg = 0.5 * (psi_in + ref_psi_out)

        visit = CellVisit(cell_idx=cell_idx, streaming_terms=st)
        strat = DiamondDifference()
        result = strat.update(visit, total_xs, source, upstream)

        assert np.array_equal(result.cell_average_flux, ref_psi_avg)
        assert np.array_equal(result.outgoing_spatial_flux, ref_psi_out)


# ═══════════════════════════════════════════════════════════════════════
# TestBitIdenticalCurvilinear
# ═══════════════════════════════════════════════════════════════════════

class TestBitIdenticalCurvilinear:
    """Curvilinear DD reproduces :func:`orpheus.sn.sweep._sweep_1d_spherical`
    bit-for-bit.

    The reference scalar form below mirrors sweep.py:328-329
    (closure constants ``c_out``, ``c_in``) and sweep.py:350-361
    (denominator, numerator, WDD + M-M closures) verbatim.  The
    structurally identical cylindrical inward / outward branches at
    sweep.py:511-531 / :548-575 use the same algebra; spherical
    coverage suffices to gate the curvilinear branch.
    """

    @pytest.mark.foundation
    @pytest.mark.verifies("dd-curvilinear-scalar", "dd-mm-closure-constants")
    def test_spherical_outward_bit_identical(self):
        """Outward-sweep cell (positive μ, away from r=0).

        Mirrors sweep.py:368-394 for a single cell with synthetic
        ``psi_spatial_in``, ``psi_angle_in``.
        """
        mesh = _spherical_mesh(nx=5, radius=1.0)
        quad = GaussLegendre1D.create(8)
        op = spherical_streaming(mesh, quad)

        ng = 2
        cell_idx = 2  # interior
        direction_idx = quad.N - 2  # positive μ, not extremal
        st = op.streaming_terms(cell_idx, direction_idx)
        # Confirm we have the curvilinear bundle populated.
        assert st.alpha_in is not None
        assert st.abs_mu >= 1e-15  # non-degenerate

        total_xs = np.array([1.3, 0.9])
        weight_norm = 1.0 / quad.weights.sum()
        Q = np.array([2.1, 0.4])
        source = Q * st.volume * weight_norm
        psi_spat_in = np.array([0.21, 0.11])
        psi_angle_in = np.array([0.05, 0.03])
        upstream = UpstreamState(
            spatial_upstream=psi_spat_in,
            angular_upstream=psi_angle_in,
        )

        # Reference scalar form — mirrors sweep.py:328-329 + 350-361.
        # Outward (μ > 0): downstream face is the OUTER face.
        abs_mu = st.abs_mu
        A_inner = st.face_area_inner
        A_outer = st.face_area_outer
        A_downstream = A_outer  # outward sweep
        dA_w = st.delta_A_over_w
        alpha_in = st.alpha_in
        alpha_out = st.alpha_out
        tau = st.tau_mm
        V = st.volume
        ref_c_out = alpha_out / tau
        ref_c_in = (1.0 - tau) / tau * alpha_out + alpha_in
        ref_denom = (
            2.0 * abs_mu * A_downstream + dA_w * ref_c_out + total_xs * V
        )
        ref_numer = (
            source
            + abs_mu * (A_inner + A_outer) * psi_spat_in
            + dA_w * ref_c_in * psi_angle_in
        )
        ref_psi_avg = ref_numer / ref_denom
        ref_psi_spat_out = 2.0 * ref_psi_avg - psi_spat_in
        ref_psi_angle_out = (ref_psi_avg - (1.0 - tau) * psi_angle_in) / tau

        # Outward visit: face_area_downstream = outer face.
        visit = CellVisit(
            cell_idx=cell_idx,
            streaming_terms=st,
            face_area_downstream=A_outer,
        )
        strat = DiamondDifference()
        result = strat.update(visit, total_xs, source, upstream)

        # Bit-identical: np.array_equal, not np.allclose.
        assert np.array_equal(result.cell_average_flux, ref_psi_avg)
        assert np.array_equal(result.outgoing_spatial_flux, ref_psi_spat_out)
        assert np.array_equal(result.outgoing_angular_state, ref_psi_angle_out)

    @pytest.mark.foundation
    @pytest.mark.verifies("dd-curvilinear-scalar", "dd-mm-closure-constants")
    def test_spherical_inward_bit_identical(self):
        """Inward-sweep cell (negative μ, marching outer → inner).

        Mirrors sweep.py:336-366 — the inward branch differs only in
        which face is "downstream" (the inner face for inward, the
        outer face for outward).  The cell-update algebra is
        identical; the sweep orchestrator (now via
        :meth:`SNMesh.iter_cell_visits`) resolves the downstream face
        before issuing the visit, so the strategy sees no
        sign-of-:math:`\\mu` branching.
        """
        mesh = _spherical_mesh(nx=5, radius=1.0)
        quad = GaussLegendre1D.create(8)
        op = spherical_streaming(mesh, quad)

        ng = 2
        cell_idx = 3
        direction_idx = 1  # negative μ, second-most-negative
        st = op.streaming_terms(cell_idx, direction_idx)
        assert st.abs_mu >= 1e-15
        # Confirm the streaming terms are direction-independent
        # (geometric labels): inner == A[i], outer == A[i+1].
        assert st.face_area_inner == float(op.face_areas[cell_idx])
        assert st.face_area_outer == float(op.face_areas[cell_idx + 1])
        # Signed mu is the direction discriminator.
        assert st.mu < 0

        total_xs = np.array([0.8, 1.4])
        weight_norm = 1.0 / quad.weights.sum()
        Q = np.array([1.1, 0.7])
        source = Q * st.volume * weight_norm
        psi_spat_in = np.array([0.05, 0.02])
        psi_angle_in = np.array([0.18, 0.09])
        upstream = UpstreamState(
            spatial_upstream=psi_spat_in,
            angular_upstream=psi_angle_in,
        )

        abs_mu = st.abs_mu
        # Inward (μ < 0): downstream face is the INNER face.
        A_downstream = st.face_area_inner
        ref_c_out = st.alpha_out / st.tau_mm
        ref_c_in = (1.0 - st.tau_mm) / st.tau_mm * st.alpha_out + st.alpha_in
        ref_denom = (
            2.0 * abs_mu * A_downstream
            + st.delta_A_over_w * ref_c_out
            + total_xs * st.volume
        )
        ref_numer = (
            source
            + abs_mu * (st.face_area_inner + st.face_area_outer) * psi_spat_in
            + st.delta_A_over_w * ref_c_in * psi_angle_in
        )
        ref_psi_avg = ref_numer / ref_denom
        ref_psi_spat_out = 2.0 * ref_psi_avg - psi_spat_in
        ref_psi_angle_out = (
            (ref_psi_avg - (1.0 - st.tau_mm) * psi_angle_in) / st.tau_mm
        )

        # Inward visit: face_area_downstream = inner face.
        visit = CellVisit(
            cell_idx=cell_idx,
            streaming_terms=st,
            face_area_downstream=st.face_area_inner,
        )
        strat = DiamondDifference()
        result = strat.update(visit, total_xs, source, upstream)

        assert np.array_equal(result.cell_average_flux, ref_psi_avg)
        assert np.array_equal(result.outgoing_spatial_flux, ref_psi_spat_out)
        assert np.array_equal(result.outgoing_angular_state, ref_psi_angle_out)


# ═══════════════════════════════════════════════════════════════════════
# TestCylindricalDegenerate
# ═══════════════════════════════════════════════════════════════════════

class TestCylindricalDegenerate:
    """Cylindrical pure-azimuthal degenerate (``abs_mu < 1e-15``).

    For cylindrical 1-D radial sweeps with a product quadrature, an
    ordinate with axial direction cosine :math:`|\\mu_z| \\to 1` has
    radial direction cosine :math:`|\\eta| \\to 0`.  The cell has
    no radial face flow — the streaming term
    :math:`\\mu_x \\partial_r` vanishes — and the cell-update
    algebra collapses to the redistribution-only form at
    sweep.py:533-546.

    The strategy returns ``CellResult(outgoing_spatial_flux=None,
    ...)`` to signal "no face-flux write" to the sweep driver.
    The Morel--Montry angular closure remains active.
    """

    @pytest.mark.foundation
    @pytest.mark.verifies("dd-curvilinear-scalar", "dd-mm-closure-constants")
    def test_degenerate_cell_synthetic(self):
        """Synthetic StreamingTerms with abs_mu < 1e-15.

        Constructed directly via the
        :class:`~orpheus.geometry.reduced_operator.StreamingTerms`
        constructor with a tiny ``abs_mu`` — gates the branch
        without depending on the choice of quadrature.  The other
        fields (``alpha_in``, ``alpha_out``, ``tau_mm``,
        ``delta_A_over_w``, ``volume``) come from a real cylindrical
        streaming-terms instance so the algebra exercises a non-
        synthetic geometry.
        """
        mesh = _cylindrical_mesh(nx=4, radius=1.0)
        quad = ProductQuadrature.create(n_mu=4, n_phi=4)
        op = cylindrical_streaming(mesh, quad)

        # Pull a real cylindrical streaming-terms packet, then
        # construct a synthetic one with abs_mu below the threshold.
        # Cylindrical streaming_terms requires mu_level_idx.
        st_real = op.streaming_terms(
            cell_idx=1, direction_idx=0, mu_level_idx=0,
        )
        # Construct degenerate variant with abs_mu = 1e-16 (well
        # below the 1e-15 threshold).
        st = StreamingTerms(
            chord_length=st_real.chord_length,
            mu=0.0,
            face_area_inner=st_real.face_area_inner,
            face_area_outer=st_real.face_area_outer,
            delta_A_over_w=st_real.delta_A_over_w,
            alpha_in=st_real.alpha_in,
            alpha_out=st_real.alpha_out,
            tau_mm=st_real.tau_mm,
            volume=st_real.volume,
            abs_mu=1e-16,
        )

        ng = 2
        total_xs = np.array([0.6, 1.0])
        # ``source`` for this synthetic cell stands in for
        # ``QV[i]`` in sweep.py:539.
        source = np.array([0.05, 0.02])
        psi_spat_in = np.array([0.10, 0.04])  # ignored by branch
        psi_angle_in = np.array([0.07, 0.03])
        upstream = UpstreamState(
            spatial_upstream=psi_spat_in,
            angular_upstream=psi_angle_in,
        )

        # Reference scalar form — mirrors sweep.py:533-543.
        ref_c_out = st.alpha_out / st.tau_mm
        ref_c_in = (1.0 - st.tau_mm) / st.tau_mm * st.alpha_out + st.alpha_in
        ref_denom = st.delta_A_over_w * ref_c_out + total_xs * st.volume
        ref_numer = source + st.delta_A_over_w * ref_c_in * psi_angle_in
        ref_psi_avg = ref_numer / ref_denom
        ref_psi_angle_out = (
            (ref_psi_avg - (1.0 - st.tau_mm) * psi_angle_in) / st.tau_mm
        )

        # Degenerate visit: no spatial face flow ⇒
        # face_area_downstream is None.
        visit = CellVisit(
            cell_idx=1,
            streaming_terms=st,
            face_area_downstream=None,
        )
        strat = DiamondDifference()
        result = strat.update(visit, total_xs, source, upstream)

        # cell_average_flux bit-identical to the sweep's degenerate scalar.
        assert np.array_equal(result.cell_average_flux, ref_psi_avg)
        # outgoing_spatial_flux signals "no face flow" via None.
        assert result.outgoing_spatial_flux is None
        # outgoing_angular_state still produced via the M-M closure.
        assert result.outgoing_angular_state is not None
        assert np.array_equal(
            result.outgoing_angular_state, ref_psi_angle_out
        )

    @pytest.mark.foundation
    def test_degenerate_does_not_consume_psi_spatial_in(self):
        """Verify the degenerate branch ignores
        ``upstream_state.spatial_upstream``.

        Per sweep.py:533-546, the degenerate cell has no radial face
        flow — the ``|μ| (A_in + A_out) ψ^s_in`` term drops out.  We
        gate this by re-running the cell update with two distinct
        ``psi_spatial_in`` values; the cell-average flux must be
        identical.
        """
        mesh = _cylindrical_mesh(nx=4, radius=1.0)
        quad = ProductQuadrature.create(n_mu=4, n_phi=4)
        op = cylindrical_streaming(mesh, quad)
        st_real = op.streaming_terms(
            cell_idx=1, direction_idx=0, mu_level_idx=0,
        )
        st = StreamingTerms(
            chord_length=st_real.chord_length,
            mu=0.0,
            face_area_inner=st_real.face_area_inner,
            face_area_outer=st_real.face_area_outer,
            delta_A_over_w=st_real.delta_A_over_w,
            alpha_in=st_real.alpha_in,
            alpha_out=st_real.alpha_out,
            tau_mm=st_real.tau_mm,
            volume=st_real.volume,
            abs_mu=1e-16,
        )

        total_xs = np.array([0.6, 1.0])
        source = np.array([0.05, 0.02])
        psi_angle_in = np.array([0.07, 0.03])

        upstream_a = UpstreamState(
            spatial_upstream=np.array([0.0, 0.0]),
            angular_upstream=psi_angle_in,
        )
        upstream_b = UpstreamState(
            spatial_upstream=np.array([99.0, -42.0]),  # wildly different
            angular_upstream=psi_angle_in,
        )

        visit = CellVisit(
            cell_idx=1,
            streaming_terms=st,
            face_area_downstream=None,
        )
        strat = DiamondDifference()
        result_a = strat.update(visit, total_xs, source, upstream_a)
        result_b = strat.update(visit, total_xs, source, upstream_b)

        # cell_average_flux insensitive to spatial_upstream
        # (no radial face flow on this cell).
        assert np.array_equal(
            result_a.cell_average_flux, result_b.cell_average_flux,
        )
        # Both signal None for outgoing_spatial_flux.
        assert result_a.outgoing_spatial_flux is None
        assert result_b.outgoing_spatial_flux is None


# ═══════════════════════════════════════════════════════════════════════
# TestPositivityFailure
# ═══════════════════════════════════════════════════════════════════════

class TestPositivityFailure:
    """:class:`DiamondDifference` is **not** positivity-preserving.

    The DD spatial closure :math:`\\psi_{\\rm out} = 2\\overline{\\psi}
    - \\psi_{\\rm in}` can produce negative outgoing flux from
    positive inputs in thin / large-source cells (Lewis & Miller
    §5.3).  This test is the canonical numerical witness that
    :attr:`DiamondDifference.is_positivity_preserving` correctly
    advertises ``False``.

    Cell-thickness criterion (Lewis & Miller §5.3): negative outgoing
    flux occurs when the optical thickness :math:`\\Delta x \\Sigma_t`
    is large compared to :math:`2|\\mu|` while the source term
    dominates the upstream-flux contribution.  Wave C-extension's
    :class:`Step` and :class:`ExponentialCharacteristic` strategies
    will be positivity-preserving by construction.
    """

    @pytest.mark.foundation
    def test_thin_cell_large_source_can_produce_negative_outgoing(self):
        """Thin slab cell with large source produces negative ``psi_out``.

        Construct a slab cell with ``Δx · Σ_t >> 2|μ|`` and a large
        source.  In this regime ``a = (2|μ| - Δx Σ_t) / denom < 0``,
        and if ``psi_in`` is large compared to the source term the
        ``a · psi_in`` contribution drives ``psi_out`` negative even
        when both inputs are positive.
        """
        # Slab cell with optically thick mesh (Δx · Σ_t ≫ 2|μ|)
        mesh = _slab_mesh(nx=2, length=2.0)  # Δx = 1.0
        quad = GaussLegendre1D.create(2)  # |μ| ≈ 0.577
        op = slab_streaming(mesh, quad)
        st = op.streaming_terms(cell_idx=0, direction_idx=quad.N - 1)

        # Σ_t * Δx = 10 vs 2|μ| ≈ 1.15 — strongly thick cell.
        total_xs = np.array([10.0])
        weight_norm = 1.0 / quad.weights.sum()
        # Small source — let the upstream-flux contribution dominate.
        Q = np.array([0.01])
        source = Q * st.chord_length * weight_norm
        # Positive but bounded upstream flux.
        psi_in = np.array([1.0])
        upstream = UpstreamState(
            spatial_upstream=psi_in,
            angular_upstream=None,
        )

        visit = CellVisit(cell_idx=0, streaming_terms=st)
        strat = DiamondDifference()
        result = strat.update(visit, total_xs, source, upstream)

        # The outgoing flux should be negative — DD's positivity
        # failure mode.  This is a pre-condition check that the
        # constructed inputs hit the failure regime; if it does not,
        # the test is vacuous and needs a stiffer cell.
        assert result.outgoing_spatial_flux is not None
        assert np.any(result.outgoing_spatial_flux < 0.0), (
            "Thin-cell DD positivity-failure test did not exhibit "
            "the failure mode.  Increase optical thickness "
            "(total_xs * chord_length) or decrease source / |mu|."
        )
