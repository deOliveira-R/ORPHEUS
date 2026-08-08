r"""The Quadrature-tier fold (Q5.6): the ``quotient`` lift and ``folded_product``.

Q5.1 landed the measure fold (orbit-stabilizer weights, mass preserved,
symmetry consumed) and Q5.3 the structure descent (charts bit-copied,
order by restriction) — both gated at their own tiers
(``tests/numerics``).  This file gates the THIN LIFT to the
:class:`Quadrature` tier (the lift must delegate and add nothing) and
the :meth:`Quadrature.folded_product` factory — the σ_y-QUOTIENT of the
staggered product, the cylindrical fold's named object, whose offset is
DERIVED by requiring the quotient to be free (T25: staggered, even
``n_phi``).

The factory's output is in the CARRYING class on every level
(``march_start_structure_per_level``: no edge-node start, no η-tie) —
the exact admissibility contract ``SNMesh(CYLINDRICAL)`` demands at the
wiring step, so this gate is the factory half of that handshake.

The pole-map gate discharges the owed check from the plan's G/P
reconciliation note: the r = 0 pole continuation needs the σ_x ordinate
permutation, and folding by σ_y must not destroy it (geometrically:
σ_x maps the ξ > 0 arc onto itself) — gated, not assumed.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry.transformation import RigidMotion
from orpheus.numerics.quadrature import (
    STAGGERED,
    Quadrature,
    gauss_legendre_on_mu,
    periodic_trapezoid,
    spherical_product,
)
from orpheus.numerics.symmetry import SubgroupOfO3
from orpheus.geometry import CoordSystem
from orpheus.sn.sweep.pole_angular_closure import (
    march_start_structure_per_level,
)

pytestmark = pytest.mark.foundation

_MIRROR_Y = SubgroupOfO3.Mirror("y")


class TestQuotientLift:
    """``Quadrature.quotient`` delegates both halves and adds nothing."""

    def test_the_lift_delegates_both_halves(self):
        """The lift == the hand-spelled seam fold, array for array."""
        measure, structure = spherical_product(
            gauss_legendre_on_mu(4), periodic_trapezoid(8, shift=STAGGERED)
        )
        folded_measure = measure.quotient(_MIRROR_Y)
        folded_structure = structure.quotient(
            parent=measure, onto=folded_measure
        )

        lifted = Quadrature(
            measure=measure, level_structure=structure
        ).quotient(_MIRROR_Y)

        np.testing.assert_array_equal(
            lifted.measure.nodes, folded_measure.nodes
        )
        np.testing.assert_array_equal(
            lifted.measure.weights, folded_measure.weights
        )
        assert lifted.level_structure is not None
        for lifted_level, expected_level in zip(
            lifted.level_indices,
            folded_structure.level_indices,
            strict=True,
        ):
            np.testing.assert_array_equal(lifted_level, expected_level)

    def test_a_rule_without_structure_folds_its_measure_alone(self):
        """Lebedev (no level structure) folds; ``level_structure`` stays None.

        O_h ⊇ σ_y, so the fold is defined; Lebedev HAS nodes on the
        mirror (axis points), so this also exercises the orbit-
        stabilizer mixed-weight arm — total mass must still be 4π.
        """
        full = Quadrature.lebedev(17)
        folded = full.quotient(_MIRROR_Y)
        assert folded.level_structure is None
        assert folded.N < full.N
        np.testing.assert_allclose(
            folded.weights.sum(), 4.0 * np.pi, rtol=1e-14
        )

    def test_refolding_by_the_consumed_mirror_is_refused(self):
        """The fold consumes ITS OWN symmetry: the representatives all sit
        on the ξ > 0 side, so σ_y no longer carries the node set onto
        itself and the certificate cannot be built.  (A DIFFERENT
        surviving symmetry may still fold — σ_z's certificate builds
        from the nodes; on the product it is then the STRUCTURE that
        refuses, because σ_z merges the signed-μ_z levels — the
        level-merging arm gated at Q5.3.)"""
        folded = Quadrature.folded_product(4, 8)
        with pytest.raises(ValueError, match="quotient is defined only"):
            folded.quotient(_MIRROR_Y)

    def test_the_slab_rule_refuses_through_the_primitive(self):
        """A 1-D measure has no 3-D realization — the measure-level refusal
        passes through the lift unchanged."""
        with pytest.raises(ValueError, match="quotient is defined only"):
            Quadrature.gauss_legendre(4).quotient(SubgroupOfO3.Mirror("x"))


class TestFoldedProduct:
    """The factory: the derived-offset product, folded free."""

    @pytest.mark.parametrize("n_phi", [8, 16, 32])
    def test_every_level_is_a_carrying_arc(self, n_phi: int):
        """The admissibility contract the cylindrical mesh will demand.

        No level starts on an edge node (the fold is free — no node on
        Σ) and no level carries an η-tie (one representative per mirror
        orbit) — so every level consumes an independent ψ½ seed.
        """
        quad = Quadrature.folded_product(4, n_phi)
        starts = march_start_structure_per_level(
            quad, CoordSystem.CYLINDRICAL
        )
        for p, start in enumerate(starts):
            assert not start.on_edge_node, (
                f"level {p}: edge-node start — the fold was not free"
            )
            assert not start.degenerate, (
                f"level {p}: η-tie — a double cover survived the fold"
            )
            assert start.consumes_independent_seed

    def test_count_and_mass(self):
        """One representative per orbit; the orbit weight is the whole 2w."""
        quad = Quadrature.folded_product(4, 8)
        assert quad.N == 4 * 8 // 2
        assert len(quad.level_indices) == 4
        assert float(quad.weights.sum()) == 4.0 * np.pi  # bit-exact

    def test_odd_n_phi_is_refused(self):
        """Odd staggered puts one node per level ON the mirror (φ = π ∈ Σ)."""
        with pytest.raises(ValueError, match="even n_phi"):
            Quadrature.folded_product(4, 5)

    def test_the_pole_map_survives_the_fold(self):
        """σ_x still permutes the folded ordinates — the r = 0 pole
        continuation stays derivable (owed by the reconciliation note;
        gated, not assumed).  The permutation is an involution."""
        quad = Quadrature.folded_product(4, 8)
        pi = quad.ordinate_permutation(
            RigidMotion.reflection(normal=[1.0, 0.0, 0.0])
        )
        assert pi is not None, (
            "the σ_y fold destroyed σ_x-closure — the cylindrical pole "
            "map cannot be derived on this rule"
        )
        np.testing.assert_array_equal(
            pi.indices[pi.indices], np.arange(quad.N)
        )
