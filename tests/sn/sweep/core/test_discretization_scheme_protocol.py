"""Software-contract tests for the DiscretizationScheme protocol.

These tests exercise the **strategy contract** shipped in Round 1 of
Wave C of the SN reshape campaign — the
:class:`~orpheus.transport.spatial.scheme.DiscretizationScheme` ``Protocol`` plus
the :class:`~orpheus.transport.spatial.scheme.UpstreamState` and
:class:`~orpheus.transport.spatial.scheme.CellResult` dataclasses.  The
L1 transport math is verified transitively via the existing sweep
MMS suite; these tests are software-contract claims (runtime
checkability, immutability, slab vs curvilinear discrimination) and
are tagged ``@pytest.mark.foundation``.

No concrete strategies (e.g. ``DiamondDifference``) ship in Round 1
— they appear in Round 2 (Issue #158).  The synthetic strategies in
this file are minimal stand-ins to drive the protocol-conformance
tests.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

import numpy as np
import pytest

from orpheus.geometry import (
    BC,
    CoordSystem,
    Mesh1D,
    slab_streaming,
    spherical_streaming,
)
from orpheus.numerics.quadrature import Quadrature
from orpheus.transport.spatial.scheme import (
    CellResult,
    DiscretizationScheme,
    CellVisit,
    UpstreamState,
)
from tests.sn.sweep.core._c_surrogate import mm_constants_for_ordinate


# ═══════════════════════════════════════════════════════════════════════
# Synthetic strategies — minimal stand-ins for protocol-conformance tests
# ═══════════════════════════════════════════════════════════════════════

@dataclass(frozen=True, slots=True)
class IdentityDiscretizationScheme:
    """Synthetic strategy: trivial closure ``ψ_avg = source / Σ_t``.

    Returns ``outgoing_spatial_flux=None`` and
    ``outgoing_angular_state=None``.  Has no physical meaning — exists
    purely to drive isinstance() and trait-attribute checks against the
    Protocol.

    Traits are declared as ``ClassVar`` so they are genuine
    class-level constants (not turned into per-instance slots by
    ``@dataclass(slots=True)``).
    """

    is_linear: ClassVar[bool] = True
    is_positivity_preserving: ClassVar[bool] = True
    is_affine_scannable: ClassVar[bool] = False
    # #240 D5-0: the Protocol gained a fourth class-level trait
    # (cross-axis separability); a conforming scheme must declare it.
    transverse_coupling_is_facewise: ClassVar[bool] = False
    # #240 D5b-S2: the Protocol gained a fifth class-level trait (the per-axis
    # spatial-moment basis size); a slopeless cell-average scheme declares 1.
    spatial_basis_per_axis: ClassVar[int] = 1
    # #236 Phase 1a: the Protocol gained a sixth class-level trait (the spatial
    # diffusion-limit condition); a conforming scheme must declare it.  This
    # synthetic stand-in makes no physics claim, so the conservative ``False``.
    diffusion_limit_consistent: ClassVar[bool] = False
    # #236 Phase 1b: the Protocol gained a seventh trait (curvilinear
    # capability); this slab-only stand-in (returns angular=None) declares False.
    supports_curvilinear: ClassVar[bool] = False

    def update(
        self,
        visit: CellVisit,
        total_xs: np.ndarray,
        source: np.ndarray,
        upstream_state: UpstreamState,
    ) -> CellResult:
        del visit, upstream_state  # unused by IdentityDiscretizationScheme
        return CellResult(
            cell_average_flux=source / total_xs,
            outgoing_spatial_flux=None,
            outgoing_angular_state=None,
        )

    def residual(
        self,
        cell_avg: np.ndarray,
        visit: CellVisit,
        total_xs: np.ndarray,
        source: np.ndarray,
        upstream_state: UpstreamState,
    ) -> np.ndarray:
        """Apply-direction companion to :meth:`update`.

        ``IdentityDiscretizationScheme`` solves ``Σ_t · ψ = q``, so the residual
        is ``Σ_t · cell_avg − source``.  At ``cell_avg = source / Σ_t``
        the residual is zero by construction.
        """
        del visit, upstream_state  # unused by IdentityDiscretizationScheme
        return total_xs * cell_avg - source

@dataclass(frozen=True, slots=True)
class BadDiscretizationScheme:
    """Synthetic non-strategy: missing the ``update`` method.

    Used to verify that the runtime-checkable Protocol correctly
    rejects a class lacking the contract method.
    """

    is_linear: ClassVar[bool] = True
    is_positivity_preserving: ClassVar[bool] = True


@dataclass(frozen=True, slots=True)
class FakeCurvilinearStrategy:
    """Synthetic strategy that asserts a curvilinear-shaped input arrives.

    Verifies the shape contract:
    ``visit.streaming_terms.delta_A_over_w is not None``,
    ``upstream_state.angular_upstream.shape == (ng,)``.  Returns
    ``CellResult`` with all three fields populated (shape ``(ng,)``).
    """

    is_linear: ClassVar[bool] = True
    is_positivity_preserving: ClassVar[bool] = False
    is_affine_scannable: ClassVar[bool] = False
    # #240 D5-0: the Protocol gained a fourth class-level trait
    # (cross-axis separability); a conforming scheme must declare it.
    transverse_coupling_is_facewise: ClassVar[bool] = False
    # #240 D5b-S2: the Protocol gained a fifth class-level trait (the per-axis
    # spatial-moment basis size); a slopeless cell-average scheme declares 1.
    spatial_basis_per_axis: ClassVar[int] = 1
    # #236 Phase 1a: the Protocol gained a sixth class-level trait (the spatial
    # diffusion-limit condition); a conforming scheme must declare it.
    diffusion_limit_consistent: ClassVar[bool] = False
    # #236 Phase 1b: the Protocol gained a seventh trait (curvilinear
    # capability); this curvilinear stand-in handles angular_upstream → True.
    supports_curvilinear: ClassVar[bool] = True

    def update(
        self,
        visit: CellVisit,
        total_xs: np.ndarray,
        source: np.ndarray,
        upstream_state: UpstreamState,
    ) -> CellResult:
        st = visit.streaming_terms
        # Curvilinear shape check.  Issue #236 Step C: the M-M τ / α packing
        # on StreamingTerms was retired; the angular weight τ is now read off
        # the CellVisit (closure-owned).  The curvilinear-shape discriminator
        # is the surviving geometry redistribution factor.
        assert st.delta_A_over_w is not None, (
            "FakeCurvilinearStrategy expects curvilinear streaming terms "
            "(delta_A_over_w must be populated)."
        )
        assert st.face_area_inner is not None
        assert st.face_area_outer is not None
        assert visit.tau is not None
        assert st.volume is not None
        assert st.abs_mu is not None
        # Curvilinear non-degenerate visits carry a positive
        # resolved downstream face area (Issue #196 Step 2.5:
        # float-always; "> 0.0" replaces "is not None").
        assert visit.face_area_downstream > 0.0

        ng = total_xs.shape[0]
        assert source.shape == (ng,)
        assert upstream_state.spatial_upstream.shape == (ng,)
        assert upstream_state.angular_upstream is not None, (
            "Curvilinear cell update needs an upstream angular state."
        )
        assert upstream_state.angular_upstream.shape == (ng,)

        # Stand-in math — not physically meaningful, just a shape-correct
        # CellResult that exercises every output channel.
        avg = source / total_xs
        out_spatial = 2.0 * avg - upstream_state.spatial_upstream
        tau = visit.tau
        out_angular = (
            avg - (1.0 - tau) * upstream_state.angular_upstream
        ) / tau
        return CellResult(
            cell_average_flux=avg,
            outgoing_spatial_flux=out_spatial,
            outgoing_angular_state=out_angular,
        )

    def residual(
        self,
        cell_avg: np.ndarray,
        visit: CellVisit,
        total_xs: np.ndarray,
        source: np.ndarray,
        upstream_state: UpstreamState,
    ) -> np.ndarray:
        """Apply-direction companion to :meth:`update`.

        Matches the stand-in math: ``Σ_t · cell_avg − source``.
        Not physically meaningful — exists purely to satisfy the
        Protocol contract for protocol-conformance tests.
        """
        del visit, upstream_state  # unused by FakeCurvilinearStrategy
        return total_xs * cell_avg - source

# ═══════════════════════════════════════════════════════════════════════
# Mesh fixtures (re-used from the geometry test suite pattern)
# ═══════════════════════════════════════════════════════════════════════

def _slab_mesh() -> Mesh1D:
    return Mesh1D(
        edges=np.linspace(0.0, 1.0, 6),
        mat_ids=np.zeros(5, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )


def _spherical_mesh() -> Mesh1D:
    return Mesh1D(
        edges=np.linspace(0.0, 1.0, 6),
        mat_ids=np.zeros(5, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )


# ═══════════════════════════════════════════════════════════════════════
# Tests
# ═══════════════════════════════════════════════════════════════════════

class TestProtocolConformance:
    """``isinstance`` runtime-checkable semantics on the DiscretizationScheme Protocol."""

    @pytest.mark.foundation
    def test_identity_strategy_is_recognized(self):
        strat = IdentityDiscretizationScheme()
        assert isinstance(strat, DiscretizationScheme)

    @pytest.mark.foundation
    def test_bad_strategy_is_rejected(self):
        bad = BadDiscretizationScheme()
        assert not isinstance(bad, DiscretizationScheme)

    @pytest.mark.foundation
    def test_curvilinear_strategy_is_recognized(self):
        strat = FakeCurvilinearStrategy()
        assert isinstance(strat, DiscretizationScheme)


class TestTraitAttributes:
    """Class-level traits ``is_linear`` / ``is_positivity_preserving``."""

    @pytest.mark.foundation
    def test_identity_traits(self):
        assert IdentityDiscretizationScheme.is_linear is True
        assert IdentityDiscretizationScheme.is_positivity_preserving is True
        # Also accessible on the instance
        s = IdentityDiscretizationScheme()
        assert s.is_linear is True
        assert s.is_positivity_preserving is True

    @pytest.mark.foundation
    def test_curvilinear_strategy_traits(self):
        assert FakeCurvilinearStrategy.is_linear is True
        assert FakeCurvilinearStrategy.is_positivity_preserving is False


class TestCapabilityTraitsAreGenuineBools:
    """The four capability traits MUST be genuine ``bool``, not merely truthy.

    ``DiscretizationScheme`` is ``@runtime_checkable``; on Python 3.12+
    ``isinstance`` validates member *presence*, NOT type — so a scheme
    declaring e.g. ``transverse_coupling_is_facewise = "yes"`` would pass
    ``isinstance(x, DiscretizationScheme)`` and then read *truthy* in
    ``ScanMarch.supports``, silently mis-routing (a narrower instance of the
    very #240 D5-0 silent-DD bug the trait was minted to close).  The
    presence-only check has no teeth against this; THIS test is the teeth:
    every REGISTERED production scheme declares all four capability traits as
    genuine ``bool``.  (#240 D5-0 elegance review — the gatekeeper's
    symmetric-on-the-Protocol resolution: keep the traits symmetric, close the
    footgun with a registry-wide type assertion rather than dropping the
    member from the Protocol.)
    """

    _CAPABILITY_TRAITS = (
        "is_linear",
        "is_positivity_preserving",
        "is_affine_scannable",
        "transverse_coupling_is_facewise",
        # #236 Phase 1a: the spatial diffusion-limit condition is a bool
        # capability trait too — subject to the same #240 D5-0 truthy-non-bool
        # footgun, so the genuine-bool teeth must cover it.
        "diffusion_limit_consistent",
        # #236 Phase 1b: curvilinear capability is read in a boolean context
        # by the selection layer (`_curvilinear_capability`) — same footgun.
        "supports_curvilinear",
    )

    @pytest.mark.foundation
    def test_registered_schemes_declare_genuine_bool_traits(self):
        # Force-import the concrete schemes so they populate the registry.
        import orpheus.transport.spatial.diamond  # noqa: F401
        import orpheus.transport.spatial.linear_discontinuous  # noqa: F401
        from orpheus.transport.spatial.scheme import DiscretizationSchemeBase

        registry = DiscretizationSchemeBase.registry
        if not registry:
            pytest.fail(
                "DiscretizationSchemeBase.registry is empty — the concrete "
                "schemes did not register (import side effect missing); the "
                "bool-trait teeth would be vacuous."
            )
        for key, scheme_cls in registry.items():
            for trait in self._CAPABILITY_TRAITS:
                value = getattr(scheme_cls, trait)
                # ``isinstance(value, bool)`` is exactly "is a genuine bool":
                # ``isinstance(1, bool)`` is False (int is not bool), and
                # np.bool_ is not a Python-bool subclass — both rejected, as
                # intended (capability traits are Python-level declarations).
                if not isinstance(value, bool):
                    pytest.fail(
                        f"{scheme_cls.__name__} (key={key!r}).{trait} = "
                        f"{value!r} is {type(value).__name__}, not bool. "
                        f"Capability traits MUST be genuine bool — a truthy "
                        f"non-bool passes the presence-only @runtime_checkable "
                        f"isinstance and then silently mis-routes in "
                        f"supports() (#240 D5-0)."
                    )


class TestDataclassImmutability:
    """``UpstreamState`` and ``CellResult`` are frozen dataclasses."""

    @pytest.mark.foundation
    def test_upstream_state_holds_reference(self):
        spatial = np.array([1.0, 2.0])
        angular = np.array([3.0, 4.0])
        st = UpstreamState(
            spatial_upstream=spatial, angular_upstream=angular,
        )
        assert st.spatial_upstream is spatial
        assert st.angular_upstream is angular

    @pytest.mark.foundation
    def test_upstream_state_is_frozen(self):
        st = UpstreamState(
            spatial_upstream=np.array([1.0]),
            angular_upstream=None,
        )
        with pytest.raises(AttributeError):
            st.spatial_upstream = np.array([99.0])  # type: ignore[misc]

    @pytest.mark.foundation
    def test_upstream_state_slab_default(self):
        """Slab UpstreamState has angular_upstream defaulting to None."""
        st = UpstreamState(spatial_upstream=np.array([1.0]))
        assert st.angular_upstream is None

    @pytest.mark.foundation
    def test_cell_result_is_frozen(self):
        r = CellResult(cell_average_flux=np.array([1.0]))
        with pytest.raises(AttributeError):
            r.cell_average_flux = np.array([99.0])  # type: ignore[misc]

    @pytest.mark.foundation
    def test_cell_result_default_outputs_none(self):
        r = CellResult(cell_average_flux=np.array([1.0]))
        assert r.outgoing_spatial_flux is None
        assert r.outgoing_angular_state is None


class TestSlabVsCurvilinearDiscrimination:
    """Geometry-as-data — slab carries neutral curvature (Issue #196 Step 2.5).

    Pre-Step-2.5 the protocol's strategies branched on ``delta_A_over_w
    is None`` to discriminate slab from curvilinear.  Step 2.5 retires
    that branch: slab populates neutral values (``delta_A_over_w = 0``,
    ``face_area_inner = face_area_outer = 1``) so the unified
    cell-balance helper consumes one packet for all geometries.  This
    test pins the neutral-value contract.

    Issue #236 Step C: the M-M ``alpha_in`` / ``alpha_out`` / ``tau_mm``
    packing on ``StreamingTerms`` was retired (the angular weight is now
    closure-owned, stamped on ``CellVisit``); the surviving neutral
    geometry fields are pinned here.
    """

    @pytest.mark.foundation
    def test_slab_streaming_terms_neutral_curvature(self):
        mesh = _slab_mesh()
        quad = Quadrature.gauss_legendre(8)
        op = slab_streaming(mesh, quad)
        st = op.streaming_terms(cell_idx=0, direction_idx=0)
        assert st.delta_A_over_w == 0.0
        assert st.face_area_inner == 1.0
        assert st.face_area_outer == 1.0

    @pytest.mark.foundation
    def test_spherical_streaming_terms_have_curvature(self):
        mesh = _spherical_mesh()
        quad = Quadrature.gauss_legendre(8)
        op = spherical_streaming(mesh, quad)
        st = op.streaming_terms(cell_idx=0, direction_idx=0)
        assert st.delta_A_over_w is not None
        assert st.face_area_inner is not None
        assert st.face_area_outer is not None


class TestCurvilinearStrategyDriven:
    """End-to-end shape check: real spherical streaming_terms feed a strategy."""

    @pytest.mark.foundation
    def test_spherical_streaming_terms_drive_strategy(self):
        mesh = _spherical_mesh()
        quad = Quadrature.gauss_legendre(8)
        op = spherical_streaming(mesh, quad)
        n = 4  # μ > 0 (outward sweep) for GL with n_half=4
        st = op.streaming_terms(cell_idx=2, direction_idx=n)
        # Wrap into a CellVisit packet with sweep-direction-resolved
        # downstream face area (outward → outer face).  Issue #236 Step C:
        # the angular weight τ is stamped on the visit (closure-owned); the
        # fake strategy reads visit.tau, so stamp the independent surrogate τ.
        tau, _, _ = mm_constants_for_ordinate(op, 2, n)
        visit = CellVisit(
            cell_idx=2,
            streaming_terms=st,
            face_area_downstream=st.face_area_outer,
            tau=tau,
        )
        ng = 3
        total_xs = np.full(ng, 1.5)
        source = np.full(ng, 0.7)
        upstream = UpstreamState(
            spatial_upstream=np.full(ng, 0.4),
            angular_upstream=np.full(ng, 0.3),
        )
        strat = FakeCurvilinearStrategy()
        result = strat.update(visit, total_xs, source, upstream)

        assert isinstance(result, CellResult)
        assert result.cell_average_flux.shape == (ng,)
        assert result.outgoing_spatial_flux is not None
        assert result.outgoing_spatial_flux.shape == (ng,)
        assert result.outgoing_angular_state is not None
        assert result.outgoing_angular_state.shape == (ng,)


class TestCellVisitPacket:
    """Foundation-level checks on the new :class:`CellVisit` dataclass."""

    @pytest.mark.foundation
    def test_cell_visit_is_frozen(self):
        mesh = _slab_mesh()
        quad = Quadrature.gauss_legendre(8)
        op = slab_streaming(mesh, quad)
        st = op.streaming_terms(cell_idx=0, direction_idx=0)
        v = CellVisit(cell_idx=0, streaming_terms=st)
        with pytest.raises(AttributeError):
            v.cell_idx = 99  # type: ignore[misc]

    @pytest.mark.foundation
    def test_cell_visit_default_downstream_zero(self):
        """Default face_area_downstream is 0.0 (Issue #196 Step 2.5).

        Pre-Step-2.5 the default was ``None``; Step 2.5 flipped to
        float-always with default ``0.0`` (signalling "no downstream
        face" via geometric truth, not by sentinel value).  The
        cell-update strategies test ``> 0.0`` rather than ``is not
        None`` to decide whether the spatial closure produces an
        output.
        """
        mesh = _slab_mesh()
        quad = Quadrature.gauss_legendre(8)
        op = slab_streaming(mesh, quad)
        st = op.streaming_terms(cell_idx=0, direction_idx=0)
        v = CellVisit(cell_idx=0, streaming_terms=st)
        assert v.face_area_downstream == 0.0
        assert isinstance(v.face_area_downstream, float)

    @pytest.mark.foundation
    def test_cell_visit_curvilinear_downstream(self):
        """Outward sphere visit: downstream = outer face by convention."""
        mesh = _spherical_mesh()
        quad = Quadrature.gauss_legendre(8)
        op = spherical_streaming(mesh, quad)
        n = 4  # μ > 0 outward
        st = op.streaming_terms(cell_idx=2, direction_idx=n)
        v = CellVisit(
            cell_idx=2,
            streaming_terms=st,
            face_area_downstream=st.face_area_outer,
        )
        assert v.face_area_downstream == st.face_area_outer
        assert v.face_area_downstream != st.face_area_inner
