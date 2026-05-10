"""Foundation tests for the BoundaryFaceFlux Protocol + concrete strategies.

These tests pin the **software contract** of the Phase A boundary-face-flux
strategies introduced for Issue #168 (the curvilinear FD operator's
outer-face truncation defects).  They are foundation-tagged because the
claims are about the strategy API (Protocol conformance, registry
self-registration, hand-calc closure algebra) rather than transport-equation
identities — those are verified transitively via the operator-level tests
at :file:`tests/sn/test_snstreamingoperator.py` and the curvilinear MMS
suite.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.sn.spatial.boundary_face_flux import (
    BoundaryFaceFlux,
    BoundaryFaceFluxBase,
    CellCenter,
    DDExtrapolation,
)

pytestmark = pytest.mark.foundation


# ═══════════════════════════════════════════════════════════════════════
# Protocol conformance
# ═══════════════════════════════════════════════════════════════════════


class TestProtocolConformance:
    """The :class:`BoundaryFaceFlux` Protocol shape is honoured by both
    concrete strategies via structural typing — no inheritance needed."""

    def test_dd_extrapolation_satisfies_protocol(self) -> None:
        """``DDExtrapolation`` is a structural ``BoundaryFaceFlux``."""
        assert isinstance(DDExtrapolation(), BoundaryFaceFlux)

    def test_cell_center_satisfies_protocol(self) -> None:
        """``CellCenter`` is a structural ``BoundaryFaceFlux``."""
        assert isinstance(CellCenter(), BoundaryFaceFlux)

    def test_is_linear_class_attr_advertised(self) -> None:
        """Both Phase-A strategies advertise ``is_linear = True``."""
        assert DDExtrapolation.is_linear is True
        assert CellCenter.is_linear is True

    def test_dd_extrapolation_repr(self) -> None:
        """``__repr__`` is the bare class name (no parameter clutter)."""
        assert repr(DDExtrapolation()) == "DDExtrapolation()"

    def test_cell_center_repr(self) -> None:
        assert repr(CellCenter()) == "CellCenter()"


# ═══════════════════════════════════════════════════════════════════════
# DD extrapolation hand-calc — outer face
# ═══════════════════════════════════════════════════════════════════════


class TestDDExtrapolationOuterFace:
    """Hand-calc tests for :class:`DDExtrapolation` at the outer face.

    Formula: ``psi_face_out = 1.5 * psi_cells[:, nx-1] - 0.5 * psi_cells[:, nx-2]``.
    """

    def test_outer_face_single_group(self) -> None:
        """Single-group (ng=1) outer-face DD extrapolation."""
        psi = np.array([[1.0, 2.0, 3.0, 5.0]])  # (ng=1, nx=4)
        # Outer at cell_idx = 3: 1.5*5 - 0.5*3 = 6.0
        result = DDExtrapolation()(psi, 3)
        np.testing.assert_array_equal(result, np.array([6.0]))

    def test_outer_face_multi_group(self) -> None:
        """Multi-group (ng=3) outer-face DD extrapolation."""
        psi = np.array([
            [1.0, 2.0, 3.0],
            [10.0, 20.0, 30.0],
            [0.5, 0.7, 0.9],
        ])  # (ng=3, nx=3)
        # Outer at cell_idx = 2:
        # group 0: 1.5*3 - 0.5*2 = 3.5
        # group 1: 1.5*30 - 0.5*20 = 35.0
        # group 2: 1.5*0.9 - 0.5*0.7 = 1.0
        result = DDExtrapolation()(psi, 2)
        np.testing.assert_array_equal(result, np.array([3.5, 35.0, 1.0]))

    def test_outer_face_constant_solution_exact(self) -> None:
        """For a constant solution, DD extrapolation reproduces the
        constant to round-off (``1.5·c - 0.5·c == c`` algebraically;
        in FP the multiply-add chain costs at most 1 ULP).  This is
        the Phase-A invariant: constants pass through every boundary
        closure unchanged up to FP non-associativity."""
        c = 7.42
        psi = np.full((2, 5), c)
        result = DDExtrapolation()(psi, 4)
        np.testing.assert_allclose(result, np.full(2, c), rtol=1e-15)

    def test_outer_face_linear_solution_exact(self) -> None:
        """For a linear solution :math:`\\psi(r) = a + b\\,r` evaluated at
        cell centres ``r_i = (i + 0.5)·h``, DD extrapolation evaluates
        the boundary face flux at ``r = nx·h`` exactly.  This is the
        load-bearing claim for O(h²) truncation."""
        nx = 6
        h = 0.5
        a, b = 1.3, 2.7
        # Cell-centres: r_i = (i + 0.5) * h
        r_centres = (np.arange(nx) + 0.5) * h
        psi_cells = (a + b * r_centres).reshape(1, nx)
        # Outer face is at r = nx * h.
        r_face = nx * h
        expected = np.array([a + b * r_face])
        result = DDExtrapolation()(psi_cells, nx - 1)
        np.testing.assert_allclose(result, expected, rtol=1e-14)


# ═══════════════════════════════════════════════════════════════════════
# DD extrapolation pole — Phase A preservation
# ═══════════════════════════════════════════════════════════════════════


class TestDDExtrapolationPole:
    """At the pole (cell_idx=0), Phase A preserves the cell-centre value.
    Defect 3 (sphere-pole stencil) is deferred to Phase B."""

    def test_pole_returns_cell_centre(self) -> None:
        """At ``cell_idx = 0`` DD-extrapolation returns ``psi_cells[:, 0]``
        unchanged (the Phase-A pole closure)."""
        psi = np.array([[1.0, 2.0, 3.0, 5.0]])
        result = DDExtrapolation()(psi, 0)
        np.testing.assert_array_equal(result, np.array([1.0]))

    def test_pole_multi_group(self) -> None:
        psi = np.array([
            [4.0, 5.0, 6.0],
            [40.0, 50.0, 60.0],
        ])
        result = DDExtrapolation()(psi, 0)
        np.testing.assert_array_equal(result, np.array([4.0, 40.0]))


# ═══════════════════════════════════════════════════════════════════════
# CellCenter (legacy reproducer) hand-calc
# ═══════════════════════════════════════════════════════════════════════


class TestCellCenterReproducer:
    """:class:`CellCenter` returns the cell-centre value at any
    boundary cell — reproduces the pre-Phase-A Defect-1 behaviour for
    ablation studies."""

    def test_outer_face_returns_cell_centre(self) -> None:
        psi = np.array([[1.0, 2.0, 3.0, 5.0]])
        result = CellCenter()(psi, 3)
        np.testing.assert_array_equal(result, np.array([5.0]))

    def test_pole_returns_cell_centre(self) -> None:
        psi = np.array([[1.0, 2.0, 3.0, 5.0]])
        result = CellCenter()(psi, 0)
        np.testing.assert_array_equal(result, np.array([1.0]))

    def test_cell_center_does_NOT_match_dd_at_outer_face(self) -> None:
        """Sanity check: on a non-constant input, ``CellCenter`` and
        ``DDExtrapolation`` give different answers at the outer face.
        This is the whole point of Phase A — Defect-1 only manifests
        on non-constant solutions."""
        psi = np.array([[1.0, 2.0, 3.0, 5.0]])
        cc_value = CellCenter()(psi, 3)
        dd_value = DDExtrapolation()(psi, 3)
        assert not np.array_equal(cc_value, dd_value)
        # cc gives 5.0, dd gives 6.0 — the difference IS Defect 1.
        np.testing.assert_array_equal(cc_value, [5.0])
        np.testing.assert_array_equal(dd_value, [6.0])


# ═══════════════════════════════════════════════════════════════════════
# Registry / self-registration
# ═══════════════════════════════════════════════════════════════════════


class TestRegistry:
    """RegistryMixin self-registration for the BoundaryFaceFlux family.

    Mirrors the CellUpdate registry contract — both DD-extrapolation
    and cell-centre strategies self-register at import time and are
    name-keyed addressable via ``BoundaryFaceFluxBase.create(...)``.
    """

    def test_dd_extrapolation_registered(self) -> None:
        """DD-extrapolation is registered under ``"dd_extrapolation"``."""
        assert "dd_extrapolation" in BoundaryFaceFluxBase.registry
        assert BoundaryFaceFluxBase.registry["dd_extrapolation"] is DDExtrapolation

    def test_cell_center_registered(self) -> None:
        assert "cell_center" in BoundaryFaceFluxBase.registry
        assert BoundaryFaceFluxBase.registry["cell_center"] is CellCenter

    def test_create_dd_extrapolation_returns_instance(self) -> None:
        """``BoundaryFaceFluxBase.create("dd_extrapolation")`` returns a
        fresh ``DDExtrapolation`` instance."""
        instance = BoundaryFaceFluxBase.create("dd_extrapolation")
        assert isinstance(instance, DDExtrapolation)

    def test_create_cell_center_returns_instance(self) -> None:
        instance = BoundaryFaceFluxBase.create("cell_center")
        assert isinstance(instance, CellCenter)

    def test_create_unknown_key_raises(self) -> None:
        """Unknown registry keys raise ``KeyError`` listing available keys."""
        with pytest.raises(KeyError, match="unknown BoundaryFaceFluxBase key"):
            BoundaryFaceFluxBase.create("unknown_strategy")


# ═══════════════════════════════════════════════════════════════════════
# Frozen + slotted dataclass invariants
# ═══════════════════════════════════════════════════════════════════════


class TestImmutability:
    """Both strategies are ``@dataclass(frozen=True, slots=True)``."""

    def test_dd_extrapolation_is_frozen(self) -> None:
        """Cannot mutate attributes on a ``DDExtrapolation`` instance."""
        strategy = DDExtrapolation()
        with pytest.raises((AttributeError, TypeError)):
            strategy.is_linear = False  # type: ignore[misc]

    def test_cell_center_is_frozen(self) -> None:
        strategy = CellCenter()
        with pytest.raises((AttributeError, TypeError)):
            strategy.is_linear = False  # type: ignore[misc]
