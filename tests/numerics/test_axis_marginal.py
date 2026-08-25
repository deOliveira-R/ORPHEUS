r"""S6 gates for the axis collapse pair (CS4b S6.0 + S6.0b).

G6.1–G6.6 of the CS4b verification plan §11 S6, realized on
:meth:`FunctionSpace.retraction` / :meth:`FunctionSpace.section` — the
split epi/mono pair :math:`R = \pi_*` (fiber integration) and :math:`E`
(its measure-normalized section), ruled onto :class:`FunctionSpace` by
the user 2026-08-24; names ratified canonical (retraction / section, Mac
Lane CWM §I.5) the same day. S6.0b re-carved the REALIZATION: the pair
is the single-region indicator frame's induced output, minted at one
site (:func:`orpheus.numerics.frame._collapse_pair`), memoized on the
space, with the frame discarded (the stage-2 generator discipline) —
``TestFrameInduction`` gates that induction (tightness, gram-derivation,
the clause-2 energy refusal, the one-mint memo).

The convention discipline under test (the anti-ERR-051 design): the two
arrows differ by exactly the axis's total weight,
:math:`R^\dagger = \Sigma w \cdot E`, and carry different NAMES and
TYPES so no call site can silently swap them. Each law row names the
mutation that reddens it; the metric-loaded row (G6.3) carries its
NEGATIVE leg per vv #19 (a positive reading alone cannot discriminate
loaded from blind).

Fixtures: the synthetic 3-axis product (weights chosen non-uniform so no
cancellation flatters a law) + the REAL SN carrier for the two
equivalence gates (G6.5/G6.6 pin ``np.array_equal`` against the shipped
kernels — the implementation deliberately spells the same float ops).

vv Mode-8 discipline: assertions are function calls (``np.testing.*`` /
``pytest.raises`` / ``pytest.fail``) — canonical invocation is
``python -O``.
"""

from __future__ import annotations

import numpy as np
import numpy.testing as npt
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.axis import Axis, BasisKind
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.space import FunctionSpace
from orpheus.sn.mesh.augmented_mesh import SNMesh
from tests.sn._test_helpers import placeholder_materials

pytestmark = pytest.mark.foundation


# ── fixtures ──────────────────────────────────────────────────────────

_W_ANG = np.array([0.3, 0.7, 0.5, 0.5])          # non-uniform, Σw = 2.0
_W_SPA = np.array([0.2, 0.3, 0.4, 0.7, 1.4])     # a genuine V vector


def _product() -> FunctionSpace:
    """angular(w) ⊗ energy(counting) ⊗ spatial(V) — the synthetic carrier."""
    return FunctionSpace.of_axes(
        Axis("angular", (4,), weights=_W_ANG, kind=BasisKind.NODAL),
        Axis("energy", (2,), kind=BasisKind.NODAL),
        Axis("spatial", (5,), weights=_W_SPA, kind=BasisKind.NODAL),
    )


def _sn() -> SNMesh:
    mesh = Mesh1D(
        edges=np.array([0.0, 0.2, 0.5, 0.9, 1.6, 3.0]),
        mat_ids=np.zeros(5, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"), bc_right=BC("vacuum"),
    )
    return SNMesh(mesh, Quadrature.gauss_legendre(4), placeholder_materials(ng=2))


def _rand(shape, seed):
    return np.random.default_rng(seed).standard_normal(shape)


# ── G6.1 / G6.2 — the section laws (LAW, reference-free) ─────────────


class TestSectionLaws:
    def test_g61_retraction_of_section_is_the_identity(self):
        """R ∘ E = id on the marginal — nulp-tier, NOT bit-exact in
        general: ⚠ [M] 2026-08-24 seed sweep on THIS fixture,
        array_equal fails 844/2000 seeds (worst rel 1.5e-16 — float
        re-association of Σ w_n/Σw); the originally-pinned seed
        happened to land in the exact set, so the earlier "BIT-EXACT"
        row was seed-fragile (S7 docs-audit finding, verified). On the
        shipped SN carrier (GL4): 200/200 draws exact. Reddened by
        dropping ÷Σw from E (an O(Σw) error, far above nulp)."""
        V = _product()
        R, E = V.retraction("angular"), V.section("angular")
        phi = _rand((2, 5), 1)
        npt.assert_array_almost_equal_nulp(R.apply(E.apply(phi)), phi, nulp=1)

    def test_g62_projection_is_idempotent(self):
        """P = E ∘ R is idempotent — nulp-tier (⚠ [M] 2026-08-24 sweep:
        array_equal fails 57/200 seeds on this fixture; same
        re-association as G6.1, same S7 audit finding). Same mutation
        as G6.1 (a mis-scaled section makes P(P) = Σw²·P)."""
        V = _product()
        R, E = V.retraction("angular"), V.section("angular")
        x = _rand(V.shape, 2)
        p = E.apply(R.apply(x))
        npt.assert_array_almost_equal_nulp(E.apply(R.apply(p)), p, nulp=1)

    def test_the_marginal_space_keeps_the_remaining_measures(self):
        """R.codomain / E.domain are the SAME content — the remaining
        axes verbatim, measures intact (the marginal's metric stays
        physical)."""
        V = _product()
        R, E = V.retraction("angular"), V.section("angular")
        marginal = FunctionSpace.of_axes(*V.axes[1:])
        if R.codomain != marginal or E.domain != marginal:
            pytest.fail("marginal space must be the remaining axes verbatim")


# ── G6.3 — the metric-loaded adjoint pairing (LAW + negative leg) ────


class TestAdjointPairing:
    def test_g63_pairing_holds_on_the_physical_metrics(self):
        """⟨Rψ, φ⟩_marginal = ⟨ψ, R.Hφ⟩_product — nulp-tier ([M]
        5.6e-16 on this fixture). Reddened by dropping V or w from
        either space (the negative leg below IS that mutation, kept
        permanently per vv #19)."""
        V = _product()
        R = V.retraction("angular")
        x, phi = _rand(V.shape, 3), _rand((2, 5), 4)
        lhs = R.codomain.inner_product(R.apply(x), phi)
        rhs = V.inner_product(x, np.asarray(R.H.apply(phi)))
        npt.assert_allclose(lhs, rhs, rtol=1e-13)

    def test_g63_negative_leg_the_pairing_is_metric_loaded(self):
        """vv #19: the positive reading cannot discriminate loaded from
        blind — cite the reading under the DELIBERATELY-WRONG structure.
        Strip the spatial measure from the domain only: the pairing must
        break at O(1)."""
        stripped = FunctionSpace.of_axes(
            Axis("angular", (4,), weights=_W_ANG, kind=BasisKind.NODAL),
            Axis("energy", (2,), kind=BasisKind.NODAL),
            Axis("spatial", (5,), kind=BasisKind.NODAL),  # V dropped
        )
        V = _product()
        R = V.retraction("angular")
        x, phi = _rand(V.shape, 3), _rand((2, 5), 4)
        lhs = R.codomain.inner_product(R.apply(x), phi)          # keeps V
        rhs = stripped.inner_product(x, np.asarray(R.H.apply(phi)))
        rel = abs(lhs - rhs) / max(abs(lhs), 1e-30)
        if rel < 1e-2:
            pytest.fail(
                f"stripping V moved the pairing by only {rel:.2e} — the "
                f"gate is NOT metric-loaded on this fixture"
            )


# ── G6.4 — the two-arrow discrimination (LAW) ────────────────────────


class TestTwoArrows:
    def test_g64_hilbert_adjoint_is_sigma_w_times_the_section(self):
        """R.H == Σw · E — nulp-tier (the metric sandwich (w·V·φ)/(w·V)
        against Σw·(φ/Σw); [M] 1.1e-16 here, exact on the §9 production
        fixture). THE anti-ERR-051 row: swapping R.H for E at a call
        site is a Σw-sized error this gate names."""
        V = _product()
        R, E = V.retraction("angular"), V.section("angular")
        phi = _rand((2, 5), 5)
        npt.assert_array_almost_equal_nulp(
            np.asarray(R.H.apply(phi)), _W_ANG.sum() * E.apply(phi), nulp=4,
        )

    def test_the_two_arrows_are_different_types(self):
        """The convention lives in the TYPE system: R.H is not an
        AxisSectionOperator and E is not an adjoint — a swapped call
        site cannot type-narrow its way through."""
        from orpheus.numerics.operator import (
            AxisSectionOperator,
            AxisRetractionOperator,
        )

        V = _product()
        R, E = V.retraction("angular"), V.section("angular")
        if not isinstance(R, AxisRetractionOperator):
            pytest.fail("retraction() must mint the retraction type")
        if not isinstance(E, AxisSectionOperator):
            pytest.fail("section() must mint the section type")
        if isinstance(R.H, AxisSectionOperator):
            pytest.fail("R.H must NOT be the section type (Σw apart)")


# ── G6.5 / G6.6 — equivalence with the shipped kernels (LAW, real SN) ─


class TestShippedKernelEquivalence:
    def test_g65_retraction_is_the_angular_reduction_bit_identical(self):
        """R over the angular axis == the canonical reduction PROGRAM —
        the hand-spelled ``einsum("n,ng...->g...", w, ψ)`` written HERE,
        independent of production. np.array_equal. ⚠ Re-scoped at S6.2:
        the pre-S6.2 target (``_integrate_angular_values``) now ROUTES
        through this very retraction (the single-sourcing carve made a
        production comparison tautological), so the hand spelling is the
        surviving independent pin — plus the external value anchors in
        ``test_angular_flux.py::TestIntegrateAngular``. Reddened by
        changing the contraction order or dropping w."""
        sn = _sn()
        psi_values = _rand((sn.quad.N, sn.ng, *sn.spatial_shape), 6)
        R = sn.angular_bulk_space.retraction("angular")
        npt.assert_array_equal(
            R.apply(psi_values),
            np.einsum("n,ng...->g...", sn.quad.weights, psi_values),
        )

    def test_g66_section_is_the_from_isotropic_kernel_bit_identical(self):
        """E over the angular axis == the hand-spelled iso-projection
        kernel (÷Σw then broadcast), np.array_equal ON THIS GL4 FIXTURE
        — written HERE, independent of production. Pre-S6.2 this was
        the licence for re-keying ``from_isotropic`` through E; S6.2
        consumed it (the factory now routes through this very section),
        so the inline spelling is the surviving independent pin. ⚠ The
        identity is NOT universal: at GL8 the induced divisor is 1 ULP
        off ``weights.sum()`` (the frame-induction gate's GL8 row) —
        principled-over-bit-identical, ruled. Reddened by dropping ÷Σw
        (G6.1's mutation)."""
        sn = _sn()
        Q = _rand((sn.ng, *sn.spatial_shape), 7)
        E = sn.angular_bulk_space.section("angular")
        sum_w = float(sn.quad.weights.sum())
        expected = np.broadcast_to(
            (Q / sum_w)[None], (sn.quad.N, sn.ng, *sn.spatial_shape),
        ).copy()
        npt.assert_array_equal(E.apply(Q), expected)

    def test_the_marginal_of_the_angular_bulk_is_the_scalar_bulk(self):
        """R.codomain content-equals the carrier's scalar bulk_space —
        the two mints share the energy/spatial axes verbatim (S1), so
        the marginal IS the scalar carrier, not a lookalike."""
        sn = _sn()
        R = sn.angular_bulk_space.retraction("angular")
        if R.codomain != sn.bulk_space:
            pytest.fail(
                f"the angular marginal must BE the scalar bulk space; "
                f"got {R.codomain!r} vs {sn.bulk_space!r}"
            )


# ── axis-generic legs (the verbs are not angular-only) ───────────────


class TestAxisGeneric:
    def test_untyped_axis_is_admitted_whatever_its_label(self):
        """The clause gate reads the axis TYPE, never the label string
        (stringly dispatch rejected): a GENERIC axis labeled "energy" —
        a synthetic test factor — is admitted with clause-3 semantics.
        Counting measure (weights None IS all-ones): R is the plain
        sum, E the uniform 1/n disaggregation, R∘E=id bit-exact. The
        typed EnergyAxis refusal is TestFrameInduction's clause-2 gate."""
        V = _product()
        R, E = V.retraction("energy"), V.section("energy")
        x = _rand(V.shape, 8)
        npt.assert_array_equal(R.apply(x), x.sum(axis=1))
        phi = _rand((4, 5), 9)
        npt.assert_array_equal(R.apply(E.apply(phi)), phi)

    def test_multidim_axis_contracts_all_its_dims(self):
        """A 2-D spatial axis (one factor, shape (nx, ny), 2-D V) is
        contracted over BOTH dims with its own weights — the volume
        integral, checked against the hand einsum."""
        Vw = np.outer(np.array([0.2, 0.5, 0.9]), np.array([1.0, 2.0]))
        V = FunctionSpace.of_axes(
            Axis("energy", (2,), kind=BasisKind.NODAL),
            Axis("spatial", (3, 2), weights=Vw, kind=BasisKind.NODAL),
        )
        R = V.retraction("spatial")
        x = _rand(V.shape, 10)
        npt.assert_allclose(
            R.apply(x), np.einsum("xy,gxy->g", Vw, x), rtol=1e-15,
        )
        E = V.section("spatial")
        phi = _rand((2,), 11)
        # The section law is ULP-tier here, NOT bit-exact: the flattened
        # 2-D measure sums 6 terms of φ/W whose reassociation wobbles 1
        # ULP ([M] 2.2e-16 on this fixture) — the bit-exact G6.1 leg is
        # the 1-D angular case's, where Σ w_n/W re-associates cleanly.
        npt.assert_array_almost_equal_nulp(R.apply(E.apply(phi)), phi, nulp=4)


# ── admission (typed refusals, both legs where a positive exists) ────


class TestAdmission:
    def test_modal_axis_is_refused_with_the_average_slot_pointer(self):
        V = FunctionSpace.of_axes(
            Axis("angular", (4,), weights=_W_ANG, kind=BasisKind.NODAL),
            Axis("moment", (2,), weights=np.array([1.0, 3.0]),
                 kind=BasisKind.MODAL),
        )
        with pytest.raises(TypeError, match="MODAL.*average slot"):
            V.retraction("moment")

    def test_unknown_label_is_refused_naming_the_axes(self):
        with pytest.raises(ValueError, match="names 0 axes"):
            _product().retraction("azimuthal")

    def test_non_axis_built_space_is_refused(self):
        bare = FunctionSpace(name="legacy", shape=(4, 2))
        with pytest.raises(TypeError, match="not axis-built"):
            bare.retraction("angular")

    def test_single_axis_space_is_refused(self):
        lone = FunctionSpace.of_axes(
            Axis("angular", (4,), weights=_W_ANG, kind=BasisKind.NODAL),
        )
        with pytest.raises(ValueError, match="bare scalar"):
            lone.retraction("angular")

    def test_zero_total_weight_refuses_the_section_only(self):
        """A signed measure summing to zero has no section — the mint's
        rank-one Gram is singular (no canonical dual), so the section
        arm is unminted and the verb refuses. The RETRACTION over the
        same axis is legal (a contraction needs no division); the
        asymmetry is the point."""
        signed = FunctionSpace.of_axes(
            Axis("angular", (2,), weights=np.array([1.0, -1.0]),
                 kind=BasisKind.NODAL),
            Axis("spatial", (3,), kind=BasisKind.NODAL),
        )
        signed.retraction("angular")  # legal
        with pytest.raises(ValueError, match="zero total weight"):
            signed.section("angular")

    def test_wrong_shape_inputs_are_refused_both_directions(self):
        V = _product()
        R, E = V.retraction("angular"), V.section("angular")
        with pytest.raises(ValueError, match="full space"):
            R.apply(np.zeros((2, 5)))
        with pytest.raises(ValueError, match="marginal space"):
            E.apply(np.zeros(V.shape))


# ── S6.0b — the frame induction (the stage-2 generator's gates) ──────


class TestFrameInduction:
    """The pair is the single-region indicator frame's induced output.

    The mint (:func:`orpheus.numerics.frame._collapse_pair`) builds the
    literal ``GalerkinFrame(IndicatorBasis, axis measure)``, reads the
    induced data off it, and discards it. These gates re-build that frame
    INDEPENDENTLY and pin the minted operators against its face contents
    — the generator discipline's consistency (tightness) gate — plus the
    gram-derivation of the section divisor, the clause-2 energy refusal,
    and the one-mint memoization.
    """

    @staticmethod
    def _literal_frame():
        """The rank-one frame over the angular axis, spelled independently
        of the mint (same construction law, re-derived here so the gate
        compares two spellings, not one object with itself)."""
        from orpheus.numerics.basis.indicator_basis import IndicatorBasis
        from orpheus.numerics.frame import GalerkinFrame
        from orpheus.numerics.measure import DiscreteMeasure

        n = _W_ANG.shape[0]
        return GalerkinFrame(
            basis=IndicatorBasis(
                edges_per_axis=(np.array([-0.5, n - 0.5]),),
            ),
            measure=DiscreteMeasure(
                nodes=np.arange(n, dtype=float),
                weights=_W_ANG,
                support="index(angular)",
            ),
        )

    def test_tightness_the_pair_is_the_single_region_frames_faces(self):
        """The minted kernels ≡ the literal frame's face contents, all
        BIT-EXACT ([M] 2026-08-24 on this fixture): R vs the analysis
        content (coefficient slot 0), R.T vs analyze_transpose (the
        w-scatter), E vs reconstruction ∘ G⁻¹ (the canonical-dual
        normalization). The operator kernels are hand einsums and the
        frame path runs the basis's table einsums — different programs,
        so agreement is a real claim. Reddened by any drift between the
        operator kernels and the frame's contraction law."""
        frame = self._literal_frame()
        V = _product()
        R, E = V.retraction("angular"), V.section("angular")
        x, phi = _rand(V.shape, 12), _rand((2, 5), 13)
        analysis = frame.basis.analyze(x, frame.table, frame.measure.weights)
        npt.assert_array_equal(R.apply(x), analysis[0])
        npt.assert_array_equal(
            R.apply_transpose(phi),
            frame.basis.analyze_transpose(
                phi[None], frame.table, frame.measure.weights,
            ),
        )
        gram = frame.discrete_gram
        npt.assert_array_equal(
            E.apply(phi),
            frame.basis.reconstruct((phi / gram[0, 0])[None], frame.table),
        )

    def test_the_section_divisor_is_the_frames_discrete_gram(self):
        """The section's divisor IS the rank-one Parseval metric — the
        1×1 ``discrete_gram`` entry of the literal frame (F-0's
        inverse-discrete-Gram theorem at K=1), pinned EXACTLY against
        the frame. ⚠ Its agreement with the OLD ``weights.sum()``
        spelling is fixture-dependent: [M] 2026-08-24 (corrected by the
        S7 docs audit — the first probe's 8 fixtures skipped GL8) exact
        at GL{2,4,5,6,12,16,32,64} and 1 ULP off at GL8; the GL8 row
        below records that divergence at its honest tier. The second
        equality here pins the exact coincidence on THIS fixture
        (n=4)."""
        frame = self._literal_frame()
        E = _product().section("angular")
        if frame.discrete_gram.shape != (1, 1):
            pytest.fail(
                f"the single-region frame's Gram must be 1×1; got "
                f"{frame.discrete_gram.shape}"
            )
        npt.assert_array_equal(E.total_weight, frame.discrete_gram[0, 0])
        npt.assert_array_equal(E.total_weight, _W_ANG.sum())

    def test_gl8_divisor_is_ulp_equivalent_to_the_sum_spelling(self):
        """The principled-over-bit-identical record, falsifiable: at GL8
        the induced divisor (the frame's gram einsum) and the old
        ``weights.sum()`` spelling differ — [M] 2026-08-24, gram
        1.9999999999999998 vs sum 2.0 (1 ULP), and the section then
        differs from the pre-S6.0b iso kernel by 2.07e-16 max rel.
        This row pins the BOUND (≤ 1 ulp on the divisor; ≤ 4 nulp on
        the kernel), not the inequality — a future numpy that closes
        the gap tightens silently, which is the correct direction."""
        from orpheus.geometry import BC, CoordSystem, Mesh1D
        from orpheus.sn.mesh.augmented_mesh import SNMesh

        mesh = Mesh1D(
            edges=np.array([0.0, 0.2, 0.5, 0.9, 1.6, 3.0]),
            mat_ids=np.zeros(5, dtype=int),
            coord=CoordSystem.CARTESIAN,
            bc_left=BC("vacuum"), bc_right=BC("vacuum"),
        )
        sn = SNMesh(
            mesh, Quadrature.gauss_legendre(8), placeholder_materials(ng=2),
        )
        E = sn.angular_bulk_space.section("angular")
        sum_w = float(sn.quad.weights.sum())
        npt.assert_array_almost_equal_nulp(E.total_weight, sum_w, nulp=1)
        Q = _rand((sn.ng, *sn.spatial_shape), 40)
        old_kernel = np.broadcast_to(
            (Q / sum_w)[None], (sn.quad.N, sn.ng, *sn.spatial_shape),
        ).copy()
        npt.assert_array_almost_equal_nulp(E.apply(Q), old_kernel, nulp=4)

    def test_typed_energy_axis_is_refused_with_the_condensation_pointer(self):
        """Collapse doctrine clause 2 (partition-integration of an L¹
        class): the ENERGY axis persists at its one-cell member —
        ⟨σ̄,φ⟩ consumes the partition — so a drop-form marginal is
        refused, pointing at the machinery that owns the energy
        collapse (EnergyGrid.overlap_to, the PG condensation frames).
        The refusal reads the TYPE; TestAxisGeneric's untyped-label row
        is the other half of the clause gate."""
        from orpheus.numerics.axis import EnergyAxis

        V = FunctionSpace.of_axes(
            EnergyAxis.synthetic(2),
            Axis("spatial", (5,), weights=_W_SPA, kind=BasisKind.NODAL),
        )
        with pytest.raises(TypeError, match="condensation"):
            V.retraction("energy")
        with pytest.raises(TypeError, match="condensation"):
            V.section("energy")

    def test_both_verbs_share_one_memoized_mint(self):
        """§5.3 of the induction plan: one mint per space per axis —
        repeated verb calls return the SAME operators (identity), and
        the two arrows share ONE marginal-space instance (both
        inductions minted together at one site)."""
        V = _product()
        R1, R2 = V.retraction("angular"), V.retraction("angular")
        E1, E2 = V.section("angular"), V.section("angular")
        if R1 is not R2 or E1 is not E2:
            pytest.fail("the collapse pair must be memoized per axis label")
        if R1.codomain is not E1.domain:
            pytest.fail("the pair must share ONE marginal space (one mint)")
