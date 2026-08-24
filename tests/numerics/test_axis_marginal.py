r"""S6 gates for the axis retraction / embedding pair (CS4b S6.0).

G6.1–G6.6 of the CS4b verification plan §11 S6, realized on
:meth:`FunctionSpace.retraction` / :meth:`FunctionSpace.embedding` —
the space-level marginal pair :math:`R` (measure contraction) and
:math:`E` (its section), ruled onto :class:`FunctionSpace` by the user
2026-08-24.

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
    def test_g61_retraction_of_embedding_is_the_identity(self):
        """R ∘ E = id on the marginal — [M] BIT-EXACT (§9 2026-08-24 and
        re-measured on this fixture). Reddened by dropping ÷Σw from E
        (the section becomes the plain broadcast, R∘E = Σw·id)."""
        V = _product()
        R, E = V.retraction("angular"), V.embedding("angular")
        phi = _rand((2, 5), 1)
        npt.assert_array_equal(R.apply(E.apply(phi)), phi)

    def test_g62_projection_is_idempotent(self):
        """P = E ∘ R is idempotent — [M] BIT-EXACT. Same mutation as
        G6.1 (a mis-scaled section makes P(P) = Σw²·P)."""
        V = _product()
        R, E = V.retraction("angular"), V.embedding("angular")
        x = _rand(V.shape, 2)
        p = E.apply(R.apply(x))
        npt.assert_array_equal(E.apply(R.apply(p)), p)

    def test_the_marginal_space_keeps_the_remaining_measures(self):
        """R.codomain / E.domain are the SAME content — the remaining
        axes verbatim, measures intact (the marginal's metric stays
        physical)."""
        V = _product()
        R, E = V.retraction("angular"), V.embedding("angular")
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
        R, E = V.retraction("angular"), V.embedding("angular")
        phi = _rand((2, 5), 5)
        npt.assert_array_almost_equal_nulp(
            np.asarray(R.H.apply(phi)), _W_ANG.sum() * E.apply(phi), nulp=4,
        )

    def test_the_two_arrows_are_different_types(self):
        """The convention lives in the TYPE system: R.H is not an
        AxisEmbeddingOperator and E is not an adjoint — a swapped call
        site cannot type-narrow its way through."""
        from orpheus.numerics.operator import (
            AxisEmbeddingOperator,
            AxisRetractionOperator,
        )

        V = _product()
        R, E = V.retraction("angular"), V.embedding("angular")
        if not isinstance(R, AxisRetractionOperator):
            pytest.fail("retraction() must mint the retraction type")
        if not isinstance(E, AxisEmbeddingOperator):
            pytest.fail("embedding() must mint the embedding type")
        if isinstance(R.H, AxisEmbeddingOperator):
            pytest.fail("R.H must NOT be the section type (Σw apart)")


# ── G6.5 / G6.6 — equivalence with the shipped kernels (LAW, real SN) ─


class TestShippedKernelEquivalence:
    def test_g65_retraction_is_the_angular_reduction_bit_identical(self):
        """R over the angular axis == _integrate_angular_values' einsum,
        np.array_equal (the implementation spells the same contraction;
        the plan pinned the array_equal leg — vv bit-identity criterion
        3). Reddened by changing the contraction order."""
        from orpheus.transport.fields.angular_flux import AngularFlux

        sn = _sn()
        psi = AngularFlux(
            values=_rand((sn.quad.N, sn.ng, *sn.spatial_shape), 6),
            space=sn.angular_bulk_space,
        )
        R = sn.angular_bulk_space.retraction("angular")
        npt.assert_array_equal(
            R.apply(psi.values), psi._integrate_angular_values(),
        )

    def test_g66_embedding_is_the_from_isotropic_kernel_bit_identical(self):
        """E over the angular axis == the from_isotropic kernel
        (÷Σw then broadcast), np.array_equal — the S6 prototype's
        licence: re-spelling the factory through E is a pure
        re-spelling. Reddened by dropping ÷Σw (G6.1's mutation)."""
        sn = _sn()
        Q = _rand((sn.ng, *sn.spatial_shape), 7)
        E = sn.angular_bulk_space.embedding("angular")
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
    def test_energy_marginal_degrades_to_the_group_sum(self):
        """The energy axis carries the counting measure (weights None IS
        all-ones) — R is the plain group collapse, E the uniform 1/ng
        disaggregation, and R∘E=id still holds bit-exactly."""
        V = _product()
        R, E = V.retraction("energy"), V.embedding("energy")
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
        E = V.embedding("spatial")
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

    def test_zero_total_weight_refuses_the_embedding_only(self):
        """A signed measure summing to zero has no section (E divides by
        Σw) — but the RETRACTION over it is legal (a contraction needs
        no division). The asymmetry is the point."""
        signed = FunctionSpace.of_axes(
            Axis("angular", (2,), weights=np.array([1.0, -1.0]),
                 kind=BasisKind.NODAL),
            Axis("spatial", (3,), kind=BasisKind.NODAL),
        )
        signed.retraction("angular")  # legal
        with pytest.raises(ValueError, match="zero total weight"):
            signed.embedding("angular")

    def test_wrong_shape_inputs_are_refused_both_directions(self):
        V = _product()
        R, E = V.retraction("angular"), V.embedding("angular")
        with pytest.raises(ValueError, match="full space"):
            R.apply(np.zeros((2, 5)))
        with pytest.raises(ValueError, match="marginal space"):
            E.apply(np.zeros(V.shape))
