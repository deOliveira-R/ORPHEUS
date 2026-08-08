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


class TestFoldedHarmonics:
    """The folded rule's harmonic machinery binds the σ-EVEN sub-basis (Q5.6).

    On the quotient, σ_y-odd harmonics are not in the function space and
    their raw ``YᵀW`` moments are GARBAGE, not zero (`[M]` pre-fix the
    ξ-carrying l = 1 slot read +6.49 for a FLAT flux — the scattering
    kernel has no Gram division anywhere to absorb it, so the garbage
    reconstructs straight into the P1 source).  The restricted basis
    keeps the rectangular (L+1, 2L+1) layout and structurally zeroes the
    odd columns, so odd moments come out EXACT 0.0 and every shape
    contract survives; the even block stays exactly orthogonal on the
    quotient (even × even is even, which the fold integrates exactly).
    """

    def test_flat_moments_are_the_isotropic_moment_alone(self):
        """Σ w Y ψ_flat = [4π at (0,0), EXACT 0.0 elsewhere] — the +6.49
        garbage channel is structurally closed."""
        fold = Quadrature.folded_product(4, 8)
        Y = fold.spherical_harmonics(2)
        M = np.einsum("n,nlm->lm", fold.weights, np.ones(fold.N)[:, None, None] * Y)
        np.testing.assert_allclose(M[0, 0], 4.0 * np.pi, rtol=1e-14)
        M_rest = M.copy()
        M_rest[0, 0] = 0.0
        # the masked (odd) slots are EXACT zeros; the surviving even
        # slots vanish to quadrature exactness on the quotient.
        assert float(np.abs(M_rest).max()) < 1e-12
        assert float(M[1, 2]) == 0.0  # the measured garbage carrier, bit-zero

    def test_even_slots_match_the_unrestricted_table(self):
        """The mask touches ONLY the σ-odd columns — the even sector is
        the parent basis bit-for-bit."""
        from orpheus.numerics.basis.spherical_harmonic_basis import (
            MirrorEvenSphericalHarmonicBasis,
            SphericalHarmonicBasis,
        )

        fold = Quadrature.folded_product(4, 8)
        restricted = MirrorEvenSphericalHarmonicBasis(L=3, mirror_axis=1)
        full = SphericalHarmonicBasis(L=3)
        Yr = restricted.evaluate(fold.measure.nodes)
        Yf = full.evaluate(fold.measure.nodes)
        mask = restricted.even_slot_mask.astype(bool)
        np.testing.assert_array_equal(Yr[:, mask], Yf[:, mask])
        np.testing.assert_array_equal(
            Yr[:, ~mask], np.zeros_like(Yr[:, ~mask])
        )
        # the derived odd-count per l is exactly l (the σ_y parity split);
        # the l=1 masked slot is the measured ξ carrier [1, 2].
        odd = ~mask
        for l in range(4):
            assert int(odd[l].sum()) == l, f"l={l}: odd-slot count != l"
        assert odd[1, 2] and not odd[1, 0] and not odd[1, 1]

    def test_even_block_gram_is_the_continuum_metric_on_the_quotient(self):
        """The discriminating gate the route-equivalent P1 hand-ref cannot
        be: the restricted basis's discrete Gram equals the CONTINUUM
        metric diag(4π/(2l+1)) on the quotient (nonzero-diagonal block),
        which is false for the unrestricted basis (the even-odd cross
        block carries the +6.49)."""
        from orpheus.numerics.basis.spherical_harmonic_basis import (
            MirrorEvenSphericalHarmonicBasis,
        )

        fold = Quadrature.folded_product(4, 8)
        basis = MirrorEvenSphericalHarmonicBasis(L=2, mirror_axis=1)
        g = basis.mass_matrix(fold.measure).reshape(15, 15)
        mask = basis.even_slot_mask.reshape(-1).astype(bool)
        g_even = g[np.ix_(mask, mask)]
        ell = np.repeat(np.arange(3), 5)[mask]
        expected = np.diag(4.0 * np.pi / (2.0 * ell + 1.0))
        live = np.diag(g_even) > 1e-10          # drop |m|>l padding rows
        np.testing.assert_allclose(
            g_even[np.ix_(live, live)],
            expected[np.ix_(live, live)],
            rtol=0.0, atol=1e-12,
        )

    def test_the_folded_frame_analysis_is_isotropic_on_a_flat_flux(self):
        """The PRODUCTION seam: ``ScatteringOperator.frame`` wraps exactly
        ``quadrature.angular_frame(L)``, whose raw analysis face
        (``YᵀW``, no Gram division) is where the +6.49 garbage entered.
        On the folded rule that frame now binds the restricted basis:
        a flat flux analyses to the isotropic moment ALONE, the
        ξ-carrier slot bit-zero."""
        fold = Quadrature.folded_product(4, 8)
        frame = fold.angular_frame(1)
        moments = frame.analysis.apply(np.ones(fold.N))
        assert moments.shape == (2, 3)
        np.testing.assert_allclose(moments[0, 0], 4.0 * np.pi, rtol=1e-14)
        assert float(moments[1, 2]) == 0.0      # the garbage carrier
        rest = moments.copy()
        rest[0, 0] = 0.0
        assert float(np.abs(rest).max()) < 1e-12
