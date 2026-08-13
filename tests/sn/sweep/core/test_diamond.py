"""Tests for the :class:`DiamondDifference` cell-update strategy.

Round 2 of Wave C of the SN reshape campaign.  These tests pin the
**bit-identical** contract between :class:`DiamondDifference` and
the existing inlined sweep math at ``orpheus.sn.loss_representation`` (the dissolved ``sweep.py``) — they
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
  ``_sweep_1d_cumprod`` (the dissolved ``sweep.py``) /
  ``_solve_recurrence`` (the dissolved ``sweep.py``) lines 117-123 +
  208-222 bit-for-bit on synthetic per-cell inputs.
* :class:`TestBitIdenticalCurvilinear` — DD's curvilinear branch
  reproduces ``_sweep_1d_spherical`` (the dissolved ``sweep.py``) lines
  350-361 bit-for-bit.
* :class:`TestCylindricalDegenerate` — DD's degenerate branch
  reproduces ``_sweep_1d_cylindrical`` (the dissolved ``sweep.py``) lines
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
from orpheus.numerics.quadrature import Quadrature
from orpheus.transport.spatial import DiamondDifference, UpstreamState
from orpheus.transport.spatial.scheme import CellVisit

# Issue #236 Phase 2 B2 Fix 3 / Step C — the ONE shared independent
# surrogate for the M-M ``(c_in, c_out)`` constants and the ``(τ, α)``
# triple (was a private byte-identical copy here; unified with
# ``test_cell_balance_for_streaming.py`` and the production-stamp catcher).
# Step C: the geometry-side τ producer is retired; α comes from the
# operator's surviving dome.  τ comes from production
# (``pole_angular_closure.morel_montry_tau_per_level``) — this said "the
# structurally-independent ``angular_differencing.morel_montry_weights``"
# until 2026-08-12, which was wrong twice over: that wrapper was retired,
# and it had already stopped being independent (it delegated).  The τ leg
# of this surrogate is TAUTOLOGICAL; see ``_c_surrogate``'s header.
from tests.sn.sweep.core._c_surrogate import (
    c_from_constants,
    mm_constants_for_ordinate,
)


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
    """Slab DD reproduces ``_solve_recurrence`` (the dissolved ``sweep.py``)
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
        quad = Quadrature.gauss_legendre(4)
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

        # Slab visit: face_area_downstream = 1.0 (Issue #196 Phase G
        # Step 2.5 — neutral curvature; slab and curvilinear share the
        # unified cell-balance helper).
        visit = CellVisit(
            cell_idx=cell_idx, streaming_terms=st,
            face_area_downstream=1.0,
        )
        strat = DiamondDifference()
        result = strat.update(visit, total_xs, source, upstream)

        # Algebraically-identical to the cumprod recurrence
        # ``a·ψ_in + 2q/denom; ½(ψ_in + ψ_out)`` — different IEEE-754
        # operation order than the unified ``(q + 2|μ|·ψ_in)/denom;
        # 2·ψ_avg − ψ_in``.  Re-baselined to ``np.allclose(rtol=1e-13)``
        # per the migration-endpoint clause documented in the
        # pre-Step-2.5 diamond.py module docstring and the plan
        # ``.claude/plans/issue_196_phase_g_replan.md`` §"Step 2.5".
        np.testing.assert_allclose(
            result.cell_average_flux, ref_psi_avg, rtol=1e-13,
        )
        np.testing.assert_allclose(
            result.outgoing_spatial_flux, ref_psi_out, rtol=1e-13,
        )
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
        quad = Quadrature.gauss_legendre(4)
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

        visit = CellVisit(
            cell_idx=cell_idx, streaming_terms=st,
            face_area_downstream=1.0,
        )
        strat = DiamondDifference()
        result = strat.update(visit, total_xs, source, upstream)

        # Issue #196 Step 2.5: re-baseline to np.allclose(rtol=1e-13)
        # per migration-endpoint clause; see test 1 docstring.
        np.testing.assert_allclose(
            result.cell_average_flux, ref_psi_avg, rtol=1e-13,
        )
        np.testing.assert_allclose(
            result.outgoing_spatial_flux, ref_psi_out, rtol=1e-13,
        )
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
        quad = Quadrature.gauss_legendre(4)
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

        visit = CellVisit(
            cell_idx=cell_idx, streaming_terms=st,
            face_area_downstream=1.0,
        )
        strat = DiamondDifference()
        result = strat.update(visit, total_xs, source, upstream)

        # Issue #196 Step 2.5: re-baseline to np.allclose(rtol=1e-13).
        np.testing.assert_allclose(
            result.cell_average_flux, ref_psi_avg, rtol=1e-13,
        )
        np.testing.assert_allclose(
            result.outgoing_spatial_flux, ref_psi_out, rtol=1e-13,
        )


# ═══════════════════════════════════════════════════════════════════════
# TestBitIdenticalCurvilinear
# ═══════════════════════════════════════════════════════════════════════

class TestBitIdenticalCurvilinear:
    """Curvilinear DD reproduces ``_sweep_1d_spherical`` (the dissolved ``sweep.py``)
    bit-for-bit.

    The reference scalar form below mirrors sweep.py:328-329
    (closure constants ``c_out``, ``c_in``) and sweep.py:350-361
    (denominator, numerator, WDD + M-M closures) verbatim.  The
    structurally identical cylindrical inward / outward branches at
    sweep.py:511-531 / :548-575 use the same algebra; spherical
    coverage suffices to gate the curvilinear branch.
    """

    @pytest.mark.sentinel
    @pytest.mark.foundation
    @pytest.mark.verifies("dd-curvilinear-scalar", "dd-mm-closure-constants")
    def test_spherical_outward_bit_identical(self):
        """Outward-sweep cell (positive μ, away from r=0).

        Mirrors sweep.py:368-394 for a single cell with synthetic
        ``psi_spatial_in``, ``psi_angle_in``.
        """
        mesh = _spherical_mesh(nx=5, radius=1.0)
        quad = Quadrature.gauss_legendre(8)
        op = spherical_streaming(mesh, quad)

        ng = 2
        cell_idx = 2  # interior
        direction_idx = quad.N - 2  # positive μ, not extremal
        st = op.streaming_terms(cell_idx, direction_idx)
        # Confirm we have the curvilinear bundle populated.
        assert st.delta_A_over_w is not None
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
        # Outward (μ > 0): downstream face is the OUTER face.  Issue #236
        # Step C: τ / α come from the independent surrogate (geometry-τ
        # retired); c_* from the hand-transcribed formula.
        tau, alpha_in, alpha_out = mm_constants_for_ordinate(
            op, cell_idx, direction_idx,
        )
        abs_mu = st.abs_mu
        A_inner = st.face_area_inner
        A_outer = st.face_area_outer
        A_downstream = A_outer  # outward sweep
        dA_w = st.delta_A_over_w
        V = st.volume
        ref_c_in, ref_c_out = c_from_constants(tau, alpha_in, alpha_out)
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

        # Outward visit: face_area_downstream = outer face.  Issue #236
        # Phase 2 B3: DD.update reads the M-M c_in / c_out / τ off the visit;
        # stamp them with the same closure-equivalent values the reference uses.
        visit = CellVisit(
            cell_idx=cell_idx,
            streaming_terms=st,
            face_area_downstream=A_outer,
            c_in=ref_c_in, c_out=ref_c_out, tau=tau,
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
        :meth:`SNMesh.dag_walk`) resolves the downstream face
        before issuing the visit, so the strategy sees no
        sign-of-:math:`\\mu` branching.

        Why this row asserts nULP and not ``array_equal``
        -------------------------------------------------
        The reference below sums the numerator left-to-right,
        ``(source + spatial) + angular``; production groups the two
        inflow terms, ``source + (spatial + angular)``
        (``cell_balance.py`` → ``diamond.py``).  Identical in exact
        arithmetic, and the ONLY difference between the two sides.

        ⭐ The bit-identity this row used to assert was a **numerical
        coincidence, never a guarantee**: with the pre-``579d5eaf``
        Gauss-Legendre nodes the two associations happened to round
        the same way.  ``579d5eaf`` moved the GL nodes/weights by
        3/17 ULP (a verified improvement — the new GL8 integrates the
        exact rational moments 5.8× better than ``numpy.leggauss``,
        winning 7 of 7 non-trivial even moments), and the coincidence
        stopped holding.  Nothing about the cell-update algebra moved.

        ``[M]`` 2026-08-12, this fixture: the association hypothesis is
        PROVEN, not inferred — ``array_equal(production, right-assoc)``
        is ``True`` while ``array_equal(production, left-assoc)`` is
        ``False``.  Drift, and its propagation from the single 1-ULP
        numerator difference:

        =========================  ====  ==========
        quantity                   ULP   rel
        =========================  ====  ==========
        ``cell_average_flux``      1     2.059e-16
        ``outgoing_spatial_flux``  2     2.528e-16
        ``outgoing_angular_state`` 4     6.807e-16
        =========================  ====  ==========

        The amplification is dimensionally explainable:
        ``2·avg − ψ_in`` carries a mild cancellation factor
        (≈1.23), and ``(avg − (1−τ)ψ)/τ`` divides by τ ≈ ½.
        ``nulp=8`` is 2× the worst measured value — headroom for
        platform FP variation, and the same order as the
        ``< 8 × ulp`` contract ``tests/numerics/test_rules_1d.py``
        already pins Gauss-Legendre against ``numpy`` with.

        ⚠ This relaxation costs the gate **no catch power**: it exists
        to catch a wrong cell-update algebra (sign flip, wrong
        downstream face, wrong M-M closure constant), every one of
        which is an O(1) error, not an O(ULP) one.  The outward
        sibling still asserts strict ``array_equal`` — it is
        unaffected only by ordinate-index luck, so if it ever reddens
        at this scale, read this note before touching production.
        """
        mesh = _spherical_mesh(nx=5, radius=1.0)
        quad = Quadrature.gauss_legendre(8)
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

        # Issue #236 Step C: τ / α from the independent surrogate.
        tau, alpha_in, alpha_out = mm_constants_for_ordinate(
            op, cell_idx, direction_idx,
        )
        abs_mu = st.abs_mu
        # Inward (μ < 0): downstream face is the INNER face.
        A_downstream = st.face_area_inner
        ref_c_in, ref_c_out = c_from_constants(tau, alpha_in, alpha_out)
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
            (ref_psi_avg - (1.0 - tau) * psi_angle_in) / tau
        )

        # Inward visit: face_area_downstream = inner face.  Issue #236
        # Phase 2 B3: stamp the M-M c_in / c_out / τ DD.update reads off the
        # visit with the same closure-equivalent values the reference uses.
        visit = CellVisit(
            cell_idx=cell_idx,
            streaming_terms=st,
            face_area_downstream=st.face_area_inner,
            c_in=ref_c_in, c_out=ref_c_out, tau=tau,
        )
        strat = DiamondDifference()
        result = strat.update(visit, total_xs, source, upstream)

        # An inward curvilinear visit MUST produce all three outputs; a None
        # here is a contract breach, not a tolerance question. Narrowed
        # explicitly because ``assert_array_almost_equal_nulp`` types its
        # arguments more strictly than ``np.array_equal`` did, and
        # ``CellResult.outgoing_*`` are ``ndarray | None``.
        spat_out = result.outgoing_spatial_flux
        ang_out = result.outgoing_angular_state
        if spat_out is None or ang_out is None:
            raise AssertionError(
                "inward curvilinear visit returned no outgoing state "
                f"(spatial={spat_out!r}, angular={ang_out!r})"
            )

        # nULP, not array_equal — see "Why this row asserts nULP" above.
        # np.testing.* raises unconditionally, so this survives -O.
        np.testing.assert_array_almost_equal_nulp(
            result.cell_average_flux, ref_psi_avg, nulp=8,
        )
        np.testing.assert_array_almost_equal_nulp(
            spat_out, ref_psi_spat_out, nulp=8,
        )
        np.testing.assert_array_almost_equal_nulp(
            ang_out, ref_psi_angle_out, nulp=8,
        )


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
    @pytest.mark.verifies(
        "dd-curvilinear-scalar", "dd-mm-closure-constants",
        "dd-cylindrical-degenerate",
    )
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
        quad = Quadrature.folded_product(n_mu=4, n_phi=4)
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
            mu_start=st_real.mu_start,
            volume=st_real.volume,
            abs_mu=1e-16,
        )
        # Issue #236 Step C: τ / α (clamped cylinder τ) for the same
        # ordinate the real packet came from, via the independent surrogate.
        tau, alpha_in, alpha_out = mm_constants_for_ordinate(
            op, cell_idx=1, direction_idx=0, mu_level_idx=0,
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
        ref_c_in, ref_c_out = c_from_constants(tau, alpha_in, alpha_out)
        ref_denom = st.delta_A_over_w * ref_c_out + total_xs * st.volume
        ref_numer = source + st.delta_A_over_w * ref_c_in * psi_angle_in
        ref_psi_avg = ref_numer / ref_denom
        ref_psi_angle_out = (
            (ref_psi_avg - (1.0 - tau) * psi_angle_in) / tau
        )

        # Degenerate visit: no spatial face flow ⇒
        # face_area_downstream = 0.0 (geometric truth; Issue #196
        # Step 2.5 replaced the ``None`` sentinel with float 0.0).
        # Issue #236 Phase 2 B3: stamp the M-M c_in / c_out / τ DD.update
        # reads off the visit (clamped cylinder τ from the surrogate).
        visit = CellVisit(
            cell_idx=1,
            streaming_terms=st,
            face_area_downstream=0.0,
            c_in=ref_c_in, c_out=ref_c_out, tau=tau,
        )
        strat = DiamondDifference()
        result = strat.update(visit, total_xs, source, upstream)

        # Issue #196 Step 2.5: the unified cell-balance helper retains
        # the ``|μ|·A_total·ψ^s_in`` term naturally — it vanishes as
        # ``|μ| → 0`` rather than via an explicit-drop branch.  For
        # synthetic ``abs_mu = 1e-16``, the residual drift is
        # ``1e-16 * O(A_total) * O(ψ^s_in) ≈ 1e-17``, well within
        # FP-non-associativity per ``vv-principles`` §"Bit-identity vs
        # principled-equivalence".
        np.testing.assert_allclose(
            result.cell_average_flux, ref_psi_avg, rtol=1e-13, atol=1e-15,
        )
        # outgoing_spatial_flux signals "no face flow" via None.
        assert result.outgoing_spatial_flux is None
        # outgoing_angular_state still produced via the M-M closure.
        assert result.outgoing_angular_state is not None
        np.testing.assert_allclose(
            result.outgoing_angular_state, ref_psi_angle_out,
            rtol=1e-13, atol=1e-15,
        )

    @pytest.mark.foundation
    def test_degenerate_does_not_consume_psi_spatial_in(self):
        """Verify the degenerate branch's spatial-upstream sensitivity
        is FP-noise level (Issue #196 Step 2.5).

        Pre-Step-2.5 the degenerate helper explicitly dropped the
        ``|μ|·A_total·ψ^s_in`` term so the cell-average flux was
        BIT-identical for any ``psi_spatial_in``.  Step 2.5's unified
        helper keeps the term but it vanishes as ``|μ|→0`` — for
        ``abs_mu = 1e-16`` the spatial-upstream sensitivity is bounded
        by ``1e-16 · A_total · |ψ^s_in|`` ≈ ``1e-14`` for the wildly-
        scaled probe.  Principled-equivalence per ``vv-principles``.
        """
        mesh = _cylindrical_mesh(nx=4, radius=1.0)
        quad = Quadrature.folded_product(n_mu=4, n_phi=4)
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
            mu_start=st_real.mu_start,
            volume=st_real.volume,
            abs_mu=1e-16,
        )
        # Issue #236 Step C: τ / α from the independent surrogate.
        tau, alpha_in, alpha_out = mm_constants_for_ordinate(
            op, cell_idx=1, direction_idx=0, mu_level_idx=0,
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

        # Issue #236 Phase 2 B3: stamp the M-M c_in / c_out / τ DD.update
        # reads off the visit (the production denom carries dA_w·c_out, which
        # sets the |μ|→0 spatial-upstream sensitivity floor this test bounds).
        c_in_v, c_out_v = c_from_constants(tau, alpha_in, alpha_out)
        visit = CellVisit(
            cell_idx=1,
            streaming_terms=st,
            face_area_downstream=0.0,
            c_in=c_in_v, c_out=c_out_v, tau=tau,
        )
        strat = DiamondDifference()
        result_a = strat.update(visit, total_xs, source, upstream_a)
        result_b = strat.update(visit, total_xs, source, upstream_b)

        # cell_average_flux insensitive to spatial_upstream to
        # FP-noise level (no radial face flow on this cell; the
        # |μ|·A_total·ψ^s_in term vanishes as |μ|→0).
        np.testing.assert_allclose(
            result_a.cell_average_flux, result_b.cell_average_flux,
            rtol=0, atol=1e-13,
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
        quad = Quadrature.gauss_legendre(2)  # |μ| ≈ 0.577
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

        visit = CellVisit(
            cell_idx=0, streaming_terms=st, face_area_downstream=1.0,
        )
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


# ═══════════════════════════════════════════════════════════════════════
# Issue 9.6 — DiscretizationSchemeBase registry membership
# ═══════════════════════════════════════════════════════════════════════


class TestDiscretizationSchemeBaseRegistry:
    """``DiamondDifference`` self-registers under ``key="diamond_difference"``."""

    @pytest.mark.foundation
    def test_diamond_difference_registered(self) -> None:
        from orpheus.transport.spatial.scheme import DiscretizationSchemeBase
        from orpheus.transport.spatial.diamond import DiamondDifference

        assert "diamond_difference" in DiscretizationSchemeBase.registry
        assert (
            DiscretizationSchemeBase.registry["diamond_difference"] is DiamondDifference
        )

    @pytest.mark.foundation
    def test_diamond_difference_factory_returns_concrete(self) -> None:
        from orpheus.transport.spatial.scheme import DiscretizationSchemeBase
        from orpheus.transport.spatial.diamond import DiamondDifference

        instance = DiscretizationSchemeBase.create("diamond_difference")
        assert isinstance(instance, DiamondDifference)


# ═══════════════════════════════════════════════════════════════════════
# Issue #196 Phase G Step 1 replan — TestResidual
# ═══════════════════════════════════════════════════════════════════════


def _slab_visit_inputs(
    nx: int = 5,
    cell_idx: int = 2,
    direction_idx: int | None = None,
    n_groups: int = 2,
    *,
    total_xs: np.ndarray | None = None,
    Q: np.ndarray | None = None,
    psi_in: np.ndarray | None = None,
) -> tuple[
    CellVisit, np.ndarray, np.ndarray, UpstreamState,
]:
    """Slab visit + inputs for residual round-trip / linearity tests.

    Returns ``(visit, total_xs, source, upstream)``.  ``source`` is
    already weight-normalised (``Q · chord · weight_norm``) per the
    :class:`DiscretizationScheme` contract — the helper does the bookkeeping
    here so tests stay focused on the residual contract.
    """
    mesh = _slab_mesh(nx=nx, length=1.0)
    quad = Quadrature.gauss_legendre(4)
    op = slab_streaming(mesh, quad)
    if direction_idx is None:
        direction_idx = quad.N - 1
    st = op.streaming_terms(cell_idx, direction_idx)
    if total_xs is None:
        total_xs = np.linspace(0.6, 1.5, n_groups)
    if Q is None:
        Q = np.linspace(0.2, 2.4, n_groups)
    if psi_in is None:
        psi_in = np.linspace(0.05, 0.35, n_groups)
    weight_norm = 1.0 / quad.weights.sum()
    source = Q * st.chord_length * weight_norm
    upstream = UpstreamState(spatial_upstream=psi_in, angular_upstream=None)
    # Slab visit: face_area_downstream = 1.0 (Issue #196 Step 2.5
    # neutral curvature) so the unified DD body's spatial-closure
    # branch runs (slab DOES have a downstream face — the cell-edge).
    # Slab α = 0 / τ = 1 → c_in = c_out = 0.0 (#236 Phase 2 B2).
    tau, alpha_in, alpha_out = mm_constants_for_ordinate(
        op, cell_idx, direction_idx,
    )
    c_in, c_out = c_from_constants(tau, alpha_in, alpha_out)
    visit = CellVisit(
        cell_idx=cell_idx, streaming_terms=st, face_area_downstream=1.0,
        c_in=c_in, c_out=c_out,
    )
    return visit, total_xs, source, upstream


def _sphere_visit_inputs(
    nx: int = 5,
    cell_idx: int = 2,
    direction_idx: int | None = None,
    n_groups: int = 2,
    *,
    outward: bool = True,
    total_xs: np.ndarray | None = None,
    Q: np.ndarray | None = None,
    psi_spat_in: np.ndarray | None = None,
    psi_angle_in: np.ndarray | None = None,
) -> tuple[
    CellVisit, np.ndarray, np.ndarray, UpstreamState,
]:
    """Sphere visit + inputs for residual contract tests."""
    mesh = _spherical_mesh(nx=nx, radius=1.0)
    quad = Quadrature.gauss_legendre(8)
    op = spherical_streaming(mesh, quad)
    if direction_idx is None:
        direction_idx = quad.N - 2 if outward else 1
    st = op.streaming_terms(cell_idx, direction_idx)
    if total_xs is None:
        total_xs = np.linspace(0.7, 1.4, n_groups)
    if Q is None:
        Q = np.linspace(0.3, 2.1, n_groups)
    if psi_spat_in is None:
        psi_spat_in = np.linspace(0.05, 0.27, n_groups)
    if psi_angle_in is None:
        psi_angle_in = np.linspace(0.02, 0.11, n_groups)
    weight_norm = 1.0 / quad.weights.sum()
    source = Q * st.volume * weight_norm
    upstream = UpstreamState(
        spatial_upstream=psi_spat_in,
        angular_upstream=psi_angle_in,
    )
    A_down = st.face_area_outer if outward else st.face_area_inner
    tau, alpha_in, alpha_out = mm_constants_for_ordinate(
        op, cell_idx, direction_idx,
    )
    c_in, c_out = c_from_constants(tau, alpha_in, alpha_out)
    visit = CellVisit(
        cell_idx=cell_idx, streaming_terms=st, face_area_downstream=A_down,
        c_in=c_in, c_out=c_out,
    )
    return visit, total_xs, source, upstream


def _cylinder_visit_inputs(
    nx: int = 4,
    cell_idx: int = 1,
    direction_idx: int = 0,
    mu_level_idx: int = 0,
    n_groups: int = 2,
    *,
    total_xs: np.ndarray | None = None,
    Q: np.ndarray | None = None,
    psi_spat_in: np.ndarray | None = None,
    psi_angle_in: np.ndarray | None = None,
) -> tuple[
    CellVisit, np.ndarray, np.ndarray, UpstreamState,
]:
    """Cylinder (non-degenerate) visit + inputs for residual contract tests."""
    mesh = _cylindrical_mesh(nx=nx, radius=1.0)
    quad = Quadrature.folded_product(n_mu=4, n_phi=4)
    op = cylindrical_streaming(mesh, quad)
    st = op.streaming_terms(
        cell_idx=cell_idx,
        direction_idx=direction_idx,
        mu_level_idx=mu_level_idx,
    )
    if total_xs is None:
        total_xs = np.linspace(0.55, 1.25, n_groups)
    if Q is None:
        Q = np.linspace(0.2, 1.8, n_groups)
    if psi_spat_in is None:
        psi_spat_in = np.linspace(0.04, 0.18, n_groups)
    if psi_angle_in is None:
        psi_angle_in = np.linspace(0.03, 0.13, n_groups)
    weight_norm = 1.0 / quad.weights.sum()
    source = Q * st.volume * weight_norm
    upstream = UpstreamState(
        spatial_upstream=psi_spat_in,
        angular_upstream=psi_angle_in,
    )
    # Determine sweep direction from signed μ.
    A_down = (
        st.face_area_outer if (st.mu is not None and st.mu >= 0.0)
        else st.face_area_inner
    )
    tau, alpha_in, alpha_out = mm_constants_for_ordinate(
        op, cell_idx, direction_idx, mu_level_idx=mu_level_idx,
    )
    c_in, c_out = c_from_constants(tau, alpha_in, alpha_out)
    visit = CellVisit(
        cell_idx=cell_idx, streaming_terms=st, face_area_downstream=A_down,
        c_in=c_in, c_out=c_out,
    )
    return visit, total_xs, source, upstream


def _cylinder_degenerate_visit_inputs(
    nx: int = 4,
    cell_idx: int = 1,
    n_groups: int = 2,
    *,
    total_xs: np.ndarray | None = None,
    source: np.ndarray | None = None,
    psi_angle_in: np.ndarray | None = None,
) -> tuple[
    CellVisit, np.ndarray, np.ndarray, UpstreamState,
]:
    """Cylindrical pure-azimuthal degenerate visit + inputs.

    Constructed by hand from a real streaming_terms instance with
    ``abs_mu`` overridden to ``1e-16`` — the canonical pattern from
    :class:`TestCylindricalDegenerate` above.
    """
    mesh = _cylindrical_mesh(nx=nx, radius=1.0)
    quad = Quadrature.folded_product(n_mu=4, n_phi=4)
    op = cylindrical_streaming(mesh, quad)
    st_real = op.streaming_terms(
        cell_idx=cell_idx, direction_idx=0, mu_level_idx=0,
    )
    st = StreamingTerms(
        chord_length=st_real.chord_length,
        mu=0.0,
        face_area_inner=st_real.face_area_inner,
        face_area_outer=st_real.face_area_outer,
        delta_A_over_w=st_real.delta_A_over_w,
        mu_start=st_real.mu_start,
        volume=st_real.volume,
        abs_mu=1e-16,
    )
    # Issue #236 Step C: τ / α from the independent surrogate.
    tau, alpha_in, alpha_out = mm_constants_for_ordinate(
        op, cell_idx=cell_idx, direction_idx=0, mu_level_idx=0,
    )
    if total_xs is None:
        total_xs = np.linspace(0.5, 1.2, n_groups)
    if source is None:
        source = np.linspace(0.02, 0.07, n_groups)
    if psi_angle_in is None:
        psi_angle_in = np.linspace(0.03, 0.13, n_groups)
    upstream = UpstreamState(
        spatial_upstream=np.zeros(n_groups),
        angular_upstream=psi_angle_in,
    )
    c_in, c_out = c_from_constants(tau, alpha_in, alpha_out)
    visit = CellVisit(
        cell_idx=cell_idx, streaming_terms=st, face_area_downstream=0.0,
        c_in=c_in, c_out=c_out,
    )
    return visit, total_xs, source, upstream


# Geometry-keyed visit factory for parametrized tests.
_GEOMETRY_FACTORIES = {
    "slab": _slab_visit_inputs,
    "sphere_outward": lambda **kw: _sphere_visit_inputs(outward=True, **kw),
    "sphere_inward": lambda **kw: _sphere_visit_inputs(outward=False, **kw),
    "cylinder": _cylinder_visit_inputs,
    "cylinder_degenerate": _cylinder_degenerate_visit_inputs,
}


class TestResidual:
    r""":class:`DiamondDifference` apply-direction contract.

    Issue #196 Phase G Step 1 replan — pins the contract on
    :meth:`DiamondDifference.residual` (the strategy-layer apply
    direction).  These tests are the salvageable invariants from the
    reverted ``test_sncell_operator.py``, ported to test the strategy
    directly per the corrected architecture (Pattern 6 — defer
    abstraction; the per-cell strategy is not a LinearOperator).

    Coverage matrix
    ---------------

    Five geometry classes × varied multi-group bias × varied source
    activation × heterogeneous total_xs:

    * **slab** — Cartesian DD, closed-form residual.
    * **sphere_outward** — curvilinear non-degenerate, outward sweep.
    * **sphere_inward** — curvilinear non-degenerate, inward sweep
      (downstream face = inner).
    * **cylinder** — curvilinear non-degenerate (1-D radial product
      quadrature, μ-level 0).
    * **cylinder_degenerate** — ``abs_mu < 1e-15`` branch (axial
      cosine :math:`|\mu_z| \to 1`).

    Mode 7 (`vv-principles` MMS simplification bias): each test
    declares which terms are activated AND which are nulled, so
    ansatz-driven cancellation is visible by inspection.
    """

    GEOMETRIES = list(_GEOMETRY_FACTORIES.keys())

    # ── 1. Round-trip: residual at solved cell_avg is zero ─────────

    @pytest.mark.foundation
    @pytest.mark.parametrize("geometry", GEOMETRIES)
    def test_residual_zero_at_solved_cell_avg(self, geometry: str) -> None:
        r"""Apply-vs-solve round-trip: ``residual(update(q).cell_avg, q) ≡ 0``.

        Strategy-layer Pattern 2 contract — :meth:`update` and
        :meth:`residual` describe the same per-cell linear system;
        evaluating the residual at the solved cell-average flux must
        return zero to FP rounding.  Tolerance ``atol=1e-13`` reflects
        the principled-equivalence band of a single division ULP.

        Term activation: every active solver term in each branch is
        exercised — slab has streaming + collision; sphere /
        cylinder have streaming + collision + Bailey
        :math:`\Delta A / w` + M-M (via ``c_in``, ``c_out``);
        cylindrical-degenerate has collision + Bailey + M-M only
        (no radial streaming).  No term is nulled by ansatz.
        """
        visit, total_xs, source, upstream = _GEOMETRY_FACTORIES[geometry]()
        strat = DiamondDifference()

        result = strat.update(visit, total_xs, source, upstream)
        residual = strat.residual(
            result.cell_average_flux,
            visit, total_xs, source, upstream,
        )

        # Round-trip identity to FP rounding — one division ULP band.
        np.testing.assert_allclose(residual, 0.0, atol=1e-13)

    @pytest.mark.sentinel
    @pytest.mark.foundation
    @pytest.mark.parametrize("geometry", GEOMETRIES)
    @pytest.mark.parametrize("n_groups", [1, 2, 4])
    def test_residual_zero_multi_group_heterogeneous(
        self, geometry: str, n_groups: int,
    ) -> None:
        """Round-trip holds across 1G / 2G / 4G with heterogeneous XS.

        Per `vv-principles` hygiene rule H1, 1-group is a degenerate
        regime; H2 forbids accepting homogeneous-only verification.
        This test exercises both axes: parametrised over n_groups
        and constructs *heterogeneous* per-group total_xs (varying
        magnitudes 0.6 → 1.5) so the multi-group bias is non-trivial.
        """
        total_xs = np.linspace(0.6, 1.5, n_groups)
        Q = np.linspace(0.2, 2.4, n_groups)
        if geometry == "slab":
            visit, total_xs, source, upstream = _slab_visit_inputs(
                n_groups=n_groups, total_xs=total_xs, Q=Q,
                psi_in=np.linspace(0.05, 0.35, n_groups),
            )
        elif geometry == "sphere_outward":
            visit, total_xs, source, upstream = _sphere_visit_inputs(
                n_groups=n_groups, outward=True,
                total_xs=total_xs, Q=Q,
                psi_spat_in=np.linspace(0.05, 0.27, n_groups),
                psi_angle_in=np.linspace(0.02, 0.11, n_groups),
            )
        elif geometry == "sphere_inward":
            visit, total_xs, source, upstream = _sphere_visit_inputs(
                n_groups=n_groups, outward=False,
                total_xs=total_xs, Q=Q,
                psi_spat_in=np.linspace(0.05, 0.27, n_groups),
                psi_angle_in=np.linspace(0.02, 0.11, n_groups),
            )
        elif geometry == "cylinder":
            visit, total_xs, source, upstream = _cylinder_visit_inputs(
                n_groups=n_groups,
                total_xs=total_xs, Q=Q,
                psi_spat_in=np.linspace(0.04, 0.18, n_groups),
                psi_angle_in=np.linspace(0.03, 0.13, n_groups),
            )
        elif geometry == "cylinder_degenerate":
            visit, total_xs, source, upstream = (
                _cylinder_degenerate_visit_inputs(
                    n_groups=n_groups, total_xs=total_xs,
                    source=np.linspace(0.02, 0.07, n_groups),
                    psi_angle_in=np.linspace(0.03, 0.13, n_groups),
                )
            )
        else:  # pragma: no cover — exhaustive parametrize
            raise ValueError(geometry)

        strat = DiamondDifference()
        result = strat.update(visit, total_xs, source, upstream)
        residual = strat.residual(
            result.cell_average_flux,
            visit, total_xs, source, upstream,
        )

        np.testing.assert_allclose(residual, 0.0, atol=1e-13)

    # ── 2. Linearity in cell_avg ───────────────────────────────────

    @pytest.mark.sentinel
    @pytest.mark.foundation
    @pytest.mark.parametrize("geometry", GEOMETRIES)
    def test_residual_linear_in_cell_avg(self, geometry: str) -> None:
        r"""Residual is linear in ``cell_avg`` (Diamond Difference closure
        is linear; M-M and Bailey terms enter only through fixed
        upstream-state contributions, not through ``cell_avg``).

        Sentinel (Phase S2) for ``DiamondDifference.residual`` — the
        apply-direction. The cheap sentinels for the sweep nodes did NOT
        call ``residual()`` at all (0 % of its mutants); this linearity
        identity kills the interior-arithmetic mutants (e.g. the lone
        S0 ``-``→``//`` FloorDiv survivor breaks linearity).

        Mathematical identity (DD ``is_linear = True``):

        .. math::

            r(\lambda a + (1-\lambda) b)
              \;=\; \lambda\,r(a) + (1-\lambda)\,r(b).

        Failure here signals an accidental non-linearity in the
        residual implementation (e.g. wrong dependence on
        ``upstream_state``, mistaken ``cell_avg`` usage in a
        denominator).
        """
        visit, total_xs, source, upstream = _GEOMETRY_FACTORIES[geometry]()
        strat = DiamondDifference()

        rng = np.random.default_rng(seed=42)
        n_groups = source.shape[0]
        probe_a = rng.normal(loc=1.0, scale=0.5, size=n_groups)
        probe_b = rng.normal(loc=2.0, scale=0.5, size=n_groups)

        lam = 0.37
        probe_mix = lam * probe_a + (1.0 - lam) * probe_b

        r_a = strat.residual(probe_a, visit, total_xs, source, upstream)
        r_b = strat.residual(probe_b, visit, total_xs, source, upstream)
        r_mix = strat.residual(
            probe_mix, visit, total_xs, source, upstream,
        )

        np.testing.assert_allclose(
            r_mix, lam * r_a + (1.0 - lam) * r_b, rtol=1e-12, atol=1e-13,
        )

    # ── 3. Affine in source ────────────────────────────────────────

    @pytest.mark.foundation
    @pytest.mark.parametrize("geometry", GEOMETRIES)
    def test_residual_affine_in_source(self, geometry: str) -> None:
        r"""Residual is affine in ``source`` — ``∂r/∂q = -1``.

        Mathematical identity:

        .. math::

            r(\bar\psi;\, q + \delta q) \;=\;
                r(\bar\psi;\, q) \;-\; \delta q.

        The cell residual carries the source on the RHS of the cell
        balance.  Shifting the source shifts the residual by minus
        the shift — a property that holds in all three branches
        because :func:`cell_balance_terms` and the slab closed form
        both treat the source as an affine term added on the outside.
        """
        visit, total_xs, source, upstream = _GEOMETRY_FACTORIES[geometry]()
        strat = DiamondDifference()

        rng = np.random.default_rng(seed=7)
        n_groups = source.shape[0]
        probe = rng.normal(loc=1.0, scale=0.5, size=n_groups)
        ds = rng.normal(loc=0.0, scale=0.3, size=n_groups)

        r0 = strat.residual(probe, visit, total_xs, source, upstream)
        r1 = strat.residual(probe, visit, total_xs, source + ds, upstream)

        np.testing.assert_allclose(r1, r0 - ds, rtol=1e-12, atol=1e-13)

    # ── 4. Bit-identity vs analytic zero at zero source ────────────

    @pytest.mark.foundation
    @pytest.mark.parametrize("geometry", GEOMETRIES)
    def test_residual_bit_identity_at_zero_source(
        self, geometry: str,
    ) -> None:
        r"""At ``source = 0`` and ``cell_avg = update(q=0).cell_avg``,
        the residual is zero to a tight FP band.

        This is the strongest round-trip variant — null the volumetric
        source so the only non-zero RHS contribution is
        ``numer_upstream`` (curvilinear) or ``2|μ| psi_in`` (slab).
        :meth:`update` solves the homogeneous-source balance; the
        round-trip ``residual(update(q=0).cell_avg, q=0)`` returns
        zero to a single division ULP.

        Term activation: streaming + collision + redistribution
        (curvilinear) — all active except ``source`` itself.  Tests
        that the residual algebra correctly handles the
        source-nulled limit (the historical homogeneous-only fallacy
        the brief warns about — see ``vv-principles`` hygiene rule
        H2).
        """
        visit, total_xs, _, upstream = _GEOMETRY_FACTORIES[geometry]()
        n_groups = total_xs.shape[0]
        source = np.zeros(n_groups)
        strat = DiamondDifference()

        result = strat.update(visit, total_xs, source, upstream)
        residual = strat.residual(
            result.cell_average_flux,
            visit, total_xs, source, upstream,
        )

        # Single-division ULP band.
        np.testing.assert_allclose(residual, 0.0, atol=1e-13)

    # ── 5. Slab residual matches closed-form analytic ──────────────

    @pytest.mark.foundation
    @pytest.mark.verifies("dd-slab-scalar")
    def test_slab_residual_closed_form(self) -> None:
        r"""Slab residual is exactly
        :math:`2|\mu|(\bar\psi - \psi_{\rm in})
        + \mathrm{chord}\,\Sigma_t\,\bar\psi - q`.

        Hand-calc cross-check against the closed-form expression.
        Issue #196 Step 2.5: re-baselined to ``np.allclose(rtol=1e-13)``
        because the unified residual builds via
        ``denom·cell_avg − (source + numer_upstream)`` (operation
        order in the unified helper) instead of the explicit
        ``2|μ|(cell_avg − ψ_in) + chord·Σ_t·cell_avg − source``
        closed form below.  The two forms are algebraically identical
        but ULP-different at IEEE-754.
        """
        visit, total_xs, source, upstream = _slab_visit_inputs()
        strat = DiamondDifference()

        # Probe at a non-converged cell_avg so the residual is non-trivial.
        rng = np.random.default_rng(seed=2026)
        n_groups = source.shape[0]
        cell_avg = rng.normal(loc=1.0, scale=0.5, size=n_groups)

        # Closed-form reference.
        st = visit.streaming_terms
        psi_in = upstream.spatial_upstream
        ref = (
            2.0 * st.abs_mu * (cell_avg - psi_in)
            + st.chord_length * total_xs * cell_avg
            - source
        )

        computed = strat.residual(
            cell_avg, visit, total_xs, source, upstream,
        )

        # Issue #196 Step 2.5: re-baselined per migration-endpoint
        # clause; see test docstring.
        np.testing.assert_allclose(computed, ref, rtol=1e-13)

    # ── 6. Curvilinear residual matches CellBalanceTerms by composition

    @pytest.mark.foundation
    @pytest.mark.verifies("dd-curvilinear-scalar")
    def test_curvilinear_residual_matches_cell_balance(self) -> None:
        r"""Sphere curvilinear residual equals
        ``denom · cell_avg - (source + numer_upstream)`` from the
        shared :func:`cell_balance_terms` helper — bit-identical.

        Verifies the Pattern 2 contract directly: both
        :meth:`update` and :meth:`residual` consume the same
        :func:`cell_balance_terms` helper (Issue #196 Step 2.5
        unified body — no geometry dispatch), so the residual at
        any probe point is the closed-form rearrangement of the
        balance equation the helper produces.
        """
        from orpheus.transport.spatial.cell_balance import cell_balance_terms

        visit, total_xs, source, upstream = _sphere_visit_inputs()
        strat = DiamondDifference()

        rng = np.random.default_rng(seed=11)
        n_groups = source.shape[0]
        cell_avg = rng.normal(loc=1.0, scale=0.5, size=n_groups)

        terms = cell_balance_terms(
            visit.streaming_terms,
            visit.face_area_downstream,
            total_xs,
            upstream,
            c_in=visit.c_in,
            c_out=visit.c_out,
        )
        ref = terms.denom * cell_avg - (source + terms.numer_upstream)

        computed = strat.residual(
            cell_avg, visit, total_xs, source, upstream,
        )

        # Bit-identical via shared-helper composition.
        assert np.array_equal(computed, ref)

    # ── 7. Cylindrical degenerate residual matches degenerate helper

    @pytest.mark.foundation
    @pytest.mark.verifies("dd-curvilinear-scalar")
    def test_cylindrical_degenerate_residual_matches_cell_balance(
        self,
    ) -> None:
        r"""Cylindrical degenerate residual equals
        ``denom · cell_avg - (source + numer_upstream)`` from
        the unified :func:`cell_balance_terms` (Issue #196 Step 2.5:
        the degenerate helper was retired in favour of the unified
        helper, which handles ``2|μ|·A_down = 0`` via
        ``visit.face_area_downstream = 0.0``).
        """
        from orpheus.transport.spatial.cell_balance import cell_balance_terms

        visit, total_xs, source, upstream = (
            _cylinder_degenerate_visit_inputs()
        )
        strat = DiamondDifference()

        rng = np.random.default_rng(seed=29)
        n_groups = source.shape[0]
        cell_avg = rng.normal(loc=1.0, scale=0.5, size=n_groups)

        terms = cell_balance_terms(
            visit.streaming_terms,
            visit.face_area_downstream,
            total_xs,
            upstream,
            c_in=visit.c_in,
            c_out=visit.c_out,
        )
        ref = terms.denom * cell_avg - (source + terms.numer_upstream)

        computed = strat.residual(
            cell_avg, visit, total_xs, source, upstream,
        )

        assert np.array_equal(computed, ref)

    # ── 8. Negative-μ ordinate (sphere inward) ─────────────────────

    @pytest.mark.foundation
    def test_residual_inward_sweep_uses_inner_face(self) -> None:
        r"""Inward-sphere residual differs from outward-sphere on the
        same streaming_terms (only ``face_area_downstream`` changes).

        Negative-μ guard: the residual must consume
        ``visit.face_area_downstream`` (= inner face on inward sweeps)
        and NOT the outer face.  Two calls on the same underlying
        streaming_terms but with different ``face_area_downstream``
        must produce different ``denom`` and thus a different residual.
        """
        # Use a cell with non-trivial geometric asymmetry between
        # inner and outer face areas.
        mesh = _spherical_mesh(nx=8, radius=1.0)
        quad = Quadrature.gauss_legendre(8)
        op = spherical_streaming(mesh, quad)
        st = op.streaming_terms(cell_idx=3, direction_idx=1)  # μ<0
        tau, alpha_in, alpha_out = mm_constants_for_ordinate(
            op, cell_idx=3, direction_idx=1,
        )
        total_xs = np.array([1.0, 0.5])
        source = np.array([0.1, 0.05])
        upstream = UpstreamState(
            spatial_upstream=np.array([0.2, 0.1]),
            angular_upstream=np.array([0.05, 0.02]),
        )
        cell_avg = np.array([0.7, 0.3])

        c_in, c_out = c_from_constants(tau, alpha_in, alpha_out)
        visit_inward = CellVisit(
            cell_idx=3, streaming_terms=st,
            face_area_downstream=st.face_area_inner,
            c_in=c_in, c_out=c_out,
        )
        visit_outward = CellVisit(
            cell_idx=3, streaming_terms=st,
            face_area_downstream=st.face_area_outer,
            c_in=c_in, c_out=c_out,
        )

        strat = DiamondDifference()
        r_in = strat.residual(
            cell_avg, visit_inward, total_xs, source, upstream,
        )
        r_out = strat.residual(
            cell_avg, visit_outward, total_xs, source, upstream,
        )

        # Different downstream face areas ⇒ different residuals.
        assert not np.allclose(r_in, r_out)

    # ── 9. Cylindrical-degenerate residual independent of psi_spat_in

    @pytest.mark.foundation
    def test_degenerate_residual_independent_of_spatial_upstream(self) -> None:
        """Degenerate residual's spatial-upstream sensitivity is FP-noise.

        Issue #196 Step 2.5: the unified cell-balance helper retains
        the ``|μ|·A_total·ψ^s_in`` term naturally (no explicit-drop
        branch); for cyl-degenerate ``abs_mu ≈ 1e-16`` it vanishes
        physically as ``|μ|→0``.  Two probes with wildly different
        ``spatial_upstream`` produce the same residual to FP-noise
        bound ``|μ|·A_total·|ψ^s_in| ≈ 1e-14`` for the test's wild
        probe.
        """
        visit, total_xs, source, upstream_a = (
            _cylinder_degenerate_visit_inputs()
        )
        psi_angle_in = upstream_a.angular_upstream
        upstream_b = UpstreamState(
            spatial_upstream=np.array([99.0, -42.0]),
            angular_upstream=psi_angle_in,
        )
        cell_avg = np.array([0.7, 0.3])
        strat = DiamondDifference()

        r_a = strat.residual(cell_avg, visit, total_xs, source, upstream_a)
        r_b = strat.residual(cell_avg, visit, total_xs, source, upstream_b)

        np.testing.assert_allclose(r_a, r_b, rtol=0, atol=1e-13)
