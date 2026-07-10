r"""§16.A — the starting-direction CARRIER intrinsic gates (2.5d d1, #282 route (a)).

Algebraic carrier laws for the ψ½ block of the augmented composite
:math:`V = V_{\rm bulk} \oplus V_{\rm trace} \oplus V_{\rm sd}` — the
gate spec is ``a3_solve_transpose_verification.md`` §16.A, with A1
re-ruled by **R12a** (presence iff the level's first-ordinate raw M-M
weight ``τ_raw ∈ (0, 1)`` exclusive — supersedes both the R12 letter's
μ-node-membership spelling AND §16.A A1's original cylinder-HAS leg):

* **A1** — R12a presence pinned BOTH ways: the sphere-GL composite HAS
  the block (1 level); the cylinder (BOTH production rules — the
  level-symmetric one is the R12a discriminator: μ_start ∉ its μ-nodes,
  yet τ_raw = 1.0 bit-exact makes the seed dead) and Cartesian (1-D and
  2-D) MUST NOT. The chosen absent-spelling is pinned positively:
  mesh predicate → ``None`` space → ``zeros`` omits the block → the
  leaf factory RAISES.
* **A2** — to_flat/from_flat round-trip INCLUDING the seed tail +
  the flat-length pin + the drop-slice mutation teeth (a serializer
  that silently truncates the block reddens BOTH legs).
* **A3** — zeros / algebra closure threads the block through every
  composite operation (additive source algebra; the flux torsor triple
  ``ψ⊖ψ → Δ`` / ``ψ⊕Δ → ψ`` / ``ψ+ψ → ⊥``; mixed presence raises —
  the anti-silent-drop law; ``advance`` presence law; ``copy``).
* **A4** — the STATE-METRIC gates (INVERTED from the former
  ghost-metric scope): the seed weights ARE ``G_sd = V_cell`` (SPD),
  ``apply_metric`` SCALES the block by ``V_cell`` and
  ``apply_inverse_metric`` DIVIDES by it (no masking — empty null
  space) while leaving bulk/trace bit-identical to the unseeded twin,
  and ``inner_product`` receives a NONZERO ``Σ V_cell·x·y`` seed term.

.. note:: **Mode-12 is CLOSED for the seed rows.** The ghost
   ``G_sd = 0`` put the seed rows inside the G-reciprocity functional's
   invariance group (identically blind to a seed-row transpose error —
   the false-green this A4 class used to pin); the ``G_sd = V_cell``
   fix moves them OUT of it. The Mode-12-CLOSURE gate — a seed-row flip
   now REDS G-reciprocity while the UNMUTATED nonzero-seed reciprocity
   holds — lives at
   :func:`tests.sn.sweep.curvilinear.test_282_direct_seed_fixed_point.test_mode12_g_reciprocity_catches_a_seed_row_flip`
   (production path) and the promoted ``diag_gsd_02`` dense sweep
   (canonical). See the derivation-of-record
   ``radial_characteristic_metric_gauge_derivation.md`` and
   :mod:`orpheus.numerics.spaces.radial_characteristic_space`.

vv Mode-8 discipline: ``np.testing.assert_*`` / ``pytest.raises`` only
(no bare ``assert`` on numerical claims — the canonical invocation is
``python -O``).
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.spaces.radial_characteristic_space import RadialCharacteristicSpace
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.fields.radial_characteristic_flux import RadialCharacteristicFlux
from orpheus.transport.displacements import RadialCharacteristicDisplacement
from orpheus.transport.full_field import FullField
from orpheus.transport.source_sinks import (
    AngularBoundarySourceSink,
    AngularSourceSink,
    RadialCharacteristicSourceSink,
)
from orpheus.transport.timed_full_field import TimedFullField
from tests.sn._test_helpers import cart2d_2g_nonsquare, placeholder_materials

pytestmark = pytest.mark.foundation


# ═══════════════════════════════════════════════════════════════════════
# Builders — one per R12a trichotomy class + the Cartesian controls
# ═══════════════════════════════════════════════════════════════════════


def _mesh_1d(coord: CoordSystem, quad, *, nx: int = 4, ng: int = 2) -> SNMesh:
    kwargs = {"bc_right": BC("reflective")}
    if coord is CoordSystem.CARTESIAN:
        kwargs["bc_left"] = BC("reflective")
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=coord,
        **kwargs,
    )
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _sphere(ng: int = 2, nx: int = 4) -> SNMesh:
    """Sphere-GL: τ_raw,0 ≈ 0.42 ∈ (0,1) — the CARRYING instance (1 level)."""
    return _mesh_1d(
        CoordSystem.SPHERICAL, Quadrature.gauss_legendre(4), nx=nx, ng=ng,
    )


def _cyl_level_symmetric() -> SNMesh:
    """Cylinder-LS: the R12a DISCRIMINATOR — μ_start ∉ the level μ-nodes
    (the R12 letter would carry), but duplicate-η midpoint edges collapse
    onto η₀ so τ_raw,0 = 1.0 bit-exact: the seed's thread weight (1−τ₀)
    vanishes — dead state, NO block (measured 0.0-bit solve insensitivity)."""
    return _mesh_1d(CoordSystem.CYLINDRICAL, Quadrature.level_symmetric(4))


def _cyl_product() -> SNMesh:
    """Cylinder-product: η₀ = η_{1/2} = −sinθ bit-exact (#229) → τ_raw,0 = 0
    — the seed would be a rank-duplicate of ψ₀, NO block."""
    return _mesh_1d(
        CoordSystem.CYLINDRICAL, Quadrature.product(n_mu=2, n_phi=4),
    )


def _slab() -> SNMesh:
    return _mesh_1d(CoordSystem.CARTESIAN, Quadrature.gauss_legendre(4))


def _seeded_flux_composite(sn: SNMesh, seed: int = 20260704) -> TimedFullField:
    """A random seed-carrying flux composite on ``sn`` (must be sphere)."""
    z = TimedFullField.zeros(
        interior=AngularFlux,
        boundary=AngularBoundaryFlux,
        mesh=sn,
        radial_characteristic=RadialCharacteristicFlux,
    )
    rng = np.random.default_rng(seed)
    assert z.radial_characteristic is not None  # builder precondition, not a gate
    return z._recombine(
        interior=replace(z.interior, values=rng.normal(size=z.interior.values.shape)),
        boundary=replace(
            z.boundary, values=rng.normal(size=z.boundary.values.shape),
        ),
        radial_characteristic=replace(
            z.radial_characteristic,
            values=rng.normal(size=z.radial_characteristic.values.shape),
        ),
    )


def _seeded_source_composite(sn: SNMesh, seed: int) -> FullField:
    """A random seed-carrying SOURCE composite (additive role algebra)."""
    z = FullField.zeros(
        interior=AngularSourceSink,
        boundary=AngularBoundarySourceSink,
        mesh=sn,
        radial_characteristic=RadialCharacteristicSourceSink,
    )
    rng = np.random.default_rng(seed)
    assert z.radial_characteristic is not None
    return z._recombine(
        interior=replace(z.interior, values=rng.normal(size=z.interior.values.shape)),
        boundary=replace(
            z.boundary, values=rng.normal(size=z.boundary.values.shape),
        ),
        radial_characteristic=replace(
            z.radial_characteristic,
            values=rng.normal(size=z.radial_characteristic.values.shape),
        ),
    )


# ═══════════════════════════════════════════════════════════════════════
# A1 — R12a presence, pinned BOTH ways
# ═══════════════════════════════════════════════════════════════════════


class TestA1PresenceR12a:
    def test_sphere_carries_one_seed_level(self):
        sn = _sphere()
        np.testing.assert_array_equal(sn.radial_characteristic_levels, (0,))
        space = sn.radial_characteristic_space
        if space is None:
            pytest.fail("sphere-GL mesh must carry a starting-direction space")
        # 1 level × 2 signs × (ng·nx cells + ng corner)
        expected = 1 * 2 * (sn.ng * 4 + sn.ng)
        np.testing.assert_array_equal(space.shape, (expected,))
        np.testing.assert_array_equal(space.levels, (0,))

    @pytest.mark.parametrize(
        "builder",
        [_cyl_level_symmetric, _cyl_product, _slab, cart2d_2g_nonsquare],
        ids=["cyl_level_symmetric", "cyl_product", "slab", "cart2d"],
    )
    def test_non_carrying_meshes_have_no_block(self, builder):
        """Absent-on-non-carrying, pinned positively (a future leak REDS here).

        The cylinder-LS row is the R12a discriminator: the R12 letter
        (μ_start ∉ μ-nodes) fires on it, so THIS row is what pins the
        refined predicate — a re-spelling back to node membership grows
        dead seed blocks on every level-symmetric cylinder and reds
        this gate.
        """
        sn = builder()
        np.testing.assert_array_equal(sn.radial_characteristic_levels, ())
        if sn.radial_characteristic_space is not None:
            pytest.fail(
                f"{builder.__name__}: non-carrying mesh must have "
                f"radial_characteristic_space None (R12a)"
            )
        z = TimedFullField.zeros(
            interior=AngularFlux,
            boundary=AngularBoundaryFlux,
            mesh=sn,
            radial_characteristic=RadialCharacteristicFlux,  # passed UNIFORMLY
        )
        if z.radial_characteristic is not None:
            pytest.fail(
                "composite factory must omit the ψ½ block on a "
                "non-carrying mesh even when the leaf class is passed"
            )
        with pytest.raises(ValueError, match="no starting-direction"):
            RadialCharacteristicFlux.zeros_on(sn)

    def test_sphere_composite_zeros_carry_typed_block(self):
        sn = _sphere()
        z = TimedFullField.zeros(
            interior=AngularFlux,
            boundary=AngularBoundaryFlux,
            mesh=sn,
            radial_characteristic=RadialCharacteristicFlux,
        )
        sd = z.radial_characteristic
        if not isinstance(sd, RadialCharacteristicFlux):
            pytest.fail(f"expected RadialCharacteristicFlux, got {type(sd)}")
        np.testing.assert_array_equal(sd.values, 0.0)
        np.testing.assert_array_equal(sd.cells(0, -1).shape, (sn.ng, 4))
        np.testing.assert_array_equal(sd.cells(0, +1).shape, (sn.ng, 4))
        np.testing.assert_array_equal(sd.corner(0, -1).shape, (sn.ng,))
        np.testing.assert_array_equal(sd.corner(0, +1).shape, (sn.ng,))

    def test_full_field_space_third_block_keyed_by_presence(self):
        """The composite SPACE grows by exactly the seed size on the sphere
        and is unchanged on non-carrying meshes."""
        sn = _sphere()
        seed_space = sn.radial_characteristic_space
        assert seed_space is not None  # narrowing (A1 sphere gate above pins it)
        n_interior = sn.quad.N * sn.ng * 4
        n_trace = sn.angular_trace.shape[0]
        np.testing.assert_array_equal(
            sn.full_field_space.shape,
            (n_interior + n_trace + seed_space.shape[0],),
        )
        slab = _slab()
        np.testing.assert_array_equal(
            slab.full_field_space.shape,
            (slab.quad.N * slab.ng * 4 + slab.angular_trace.shape[0],),
        )

    def test_views_alias_one_flat_backing(self):
        """cells/corner are slice VIEWS — writing through them mutates the
        flat buffer (the AngularBoundaryFlux flat-storage discipline)."""
        sn = _sphere(ng=1)
        space = sn.radial_characteristic_space
        assert space is not None
        buf = np.zeros(space.shape)
        space.cells_view(buf, 0, -1)[:] = 1.0
        space.corner_view(buf, 0, +1)[:] = 2.0
        np.testing.assert_allclose(np.sum(buf == 1.0), space.ng * space.nx)
        np.testing.assert_allclose(np.sum(buf == 2.0), space.ng)
        # Layout order (documented): (−1 cells, −1 corner, +1 cells, +1 corner)
        per_sign = space.ng * space.nx + space.ng
        np.testing.assert_array_equal(buf[: space.ng * space.nx], 1.0)
        np.testing.assert_array_equal(buf[per_sign + space.ng * space.nx :], 2.0)

    def test_unknown_level_and_sign_rejected(self):
        sn = _sphere()
        space = sn.radial_characteristic_space
        assert space is not None
        buf = np.zeros(space.shape)
        with pytest.raises(KeyError, match="no starting-direction block"):
            space.cells_view(buf, 1, -1)
        with pytest.raises(ValueError, match="sign must be"):
            space.corner_view(buf, 0, 0)

    def test_empty_levels_space_unrepresentable(self):
        """Absence is spelled None, never a zero-DOF phantom space.

        ``cell_volumes`` is now a REQUIRED kwarg (the state metric G_sd =
        V_cell), but irrelevant here — the empty-``levels`` guard fires
        first; a valid dummy keeps the call well-formed so the ValueError
        (not a TypeError from the missing kwarg) is what we pin.
        """
        with pytest.raises(ValueError, match="levels is empty"):
            RadialCharacteristicSpace.for_levels(
                (), ng=2, nx=4, cell_volumes=np.ones(4),
            )


# ═══════════════════════════════════════════════════════════════════════
# A2 — flat round-trip including the seed tail (+ the mutation teeth)
# ═══════════════════════════════════════════════════════════════════════


class TestA2FlatRoundTrip:
    def test_round_trip_identity_and_length(self):
        sn = _sphere()
        x = _seeded_flux_composite(sn)
        assert x.radial_characteristic is not None
        flat = x.to_flat()
        expected_len = (
            x.interior.values.size
            + x.boundary.values.size
            + x.radial_characteristic.values.size
        )
        np.testing.assert_array_equal(flat.size, expected_len)
        y = FullField.from_flat(flat, x)
        np.testing.assert_array_equal(y.interior.values, x.interior.values)
        np.testing.assert_array_equal(y.boundary.values, x.boundary.values)
        assert y.radial_characteristic is not None
        np.testing.assert_array_equal(
            y.radial_characteristic.values, x.radial_characteristic.values,
        )
        if type(y) is not TimedFullField:
            pytest.fail("template type must be preserved through from_flat")

    def test_seed_tail_layout_is_the_trailing_block(self):
        """The flat layout is [bulk, trace, seed] — pin the tail slice."""
        sn = _sphere()
        x = _seeded_flux_composite(sn)
        assert x.radial_characteristic is not None
        flat = x.to_flat()
        n_head = x.interior.values.size + x.boundary.values.size
        np.testing.assert_array_equal(
            flat[n_head:], x.radial_characteristic.values,
        )

    def test_mutation_teeth_dropped_seed_slice_reds_both_legs(self, monkeypatch):
        """MUTATION (in-process, per the mutation-hygiene rule): a to_flat
        that silently truncates the ψ½ block breaks (a) the length pin and
        (b) from_flat's size guard — proving the block is genuinely
        serialized, not decorative."""
        sn = _sphere()
        x = _seeded_flux_composite(sn)
        assert x.radial_characteristic is not None

        def _two_block_to_flat(self):  # the pre-2.5d serializer
            return np.concatenate([
                self.interior.values.ravel(), self.boundary.values,
            ])

        monkeypatch.setattr(FullField, "to_flat", _two_block_to_flat)
        mutated = x.to_flat()
        expected_len = (
            x.interior.values.size
            + x.boundary.values.size
            + x.radial_characteristic.values.size
        )
        if mutated.size == expected_len:
            pytest.fail("mutation did not bite — the teeth probe is broken")
        with pytest.raises(ValueError, match="does not match template"):
            FullField.from_flat(mutated, x)


# ═══════════════════════════════════════════════════════════════════════
# A3 — zeros / algebra closure threads the block
# ═══════════════════════════════════════════════════════════════════════


class TestA3AlgebraClosure:
    def test_source_composite_additive_algebra(self):
        """SourceSink composites are a genuine vector space — every
        operation threads the seed block exactly leaf-wise."""
        sn = _sphere()
        a = _seeded_source_composite(sn, 1)
        b = _seeded_source_composite(sn, 2)
        assert a.radial_characteristic is not None
        assert b.radial_characteristic is not None
        av, bv = a.radial_characteristic.values, b.radial_characteristic.values
        s = a + b
        assert s.radial_characteristic is not None
        np.testing.assert_array_equal(s.radial_characteristic.values, av + bv)
        d = a - b
        assert d.radial_characteristic is not None
        np.testing.assert_array_equal(d.radial_characteristic.values, av - bv)
        n = -a
        assert n.radial_characteristic is not None
        np.testing.assert_array_equal(n.radial_characteristic.values, -av)
        m = 2.5 * a
        assert m.radial_characteristic is not None
        np.testing.assert_array_equal(m.radial_characteristic.values, 2.5 * av)
        q = a / 4.0
        assert q.radial_characteristic is not None
        np.testing.assert_array_equal(q.radial_characteristic.values, av / 4.0)
        # role closure: the summed block keeps the SourceSink role
        if type(s.radial_characteristic) is not RadialCharacteristicSourceSink:
            pytest.fail("source algebra must stay in the SourceSink role")

    def test_flux_torsor_triple_threads_the_seed(self):
        sn = _sphere()
        psi1 = _seeded_flux_composite(sn, 11)
        psi2 = _seeded_flux_composite(sn, 12)
        assert psi1.radial_characteristic is not None
        assert psi2.radial_characteristic is not None
        # ψ ⊖ ψ → displacement composite (seed block mints its sibling)
        delta = psi2 - psi1
        if type(delta.radial_characteristic) is not RadialCharacteristicDisplacement:
            pytest.fail(
                f"flux ⊖ flux must mint RadialCharacteristicDisplacement; got "
                f"{type(delta.radial_characteristic).__name__}"
            )
        np.testing.assert_array_equal(
            delta.radial_characteristic.values,
            psi2.radial_characteristic.values - psi1.radial_characteristic.values,
        )
        # ψ ⊕ Δ → flux. The THREADING pin is bitwise (the block computes
        # exactly ψ₁ ⊕ Δ leaf-wise); the torsor LAW ψ₁ ⊕ (ψ₂ ⊖ ψ₁) = ψ₂
        # holds to FP round-off only (a + (b − a) ≠ b bitwise).
        recovered = psi1 + delta
        if type(recovered.radial_characteristic) is not RadialCharacteristicFlux:
            pytest.fail("torsor action must return the flux role")
        np.testing.assert_array_equal(
            recovered.radial_characteristic.values,
            psi1.radial_characteristic.values + delta.radial_characteristic.values,
        )
        np.testing.assert_allclose(
            recovered.radial_characteristic.values,
            psi2.radial_characteristic.values,
            rtol=0.0, atol=1e-15,
        )
        # ψ + ψ → ⊥ (the affine gate — raised by the leaf algebra)
        with pytest.raises(TypeError, match="cannot add two"):
            psi1 + psi2

    @pytest.mark.parametrize("op", ["add", "sub"], ids=["add", "sub"])
    def test_mixed_presence_raises(self, op):
        """The anti-silent-drop law: a seeded ⊕ unseeded pair is an ERROR,
        never a silent drop of the ψ½ block (§16.A A3's feared bug)."""
        sn = _sphere()
        a = _seeded_source_composite(sn, 21)
        unseeded = FullField(interior=a.interior, boundary=a.boundary)
        with pytest.raises(ValueError, match="MIXED starting-direction"):
            (a + unseeded) if op == "add" else (a - unseeded)
        with pytest.raises(ValueError, match="MIXED starting-direction"):
            (unseeded + a) if op == "add" else (unseeded - a)

    def test_advance_presence_law_and_snapshot(self):
        sn = _sphere()
        psi = _seeded_flux_composite(sn, 31)
        assert psi.radial_characteristic is not None
        new = _seeded_flux_composite(sn, 32)
        assert new.radial_characteristic is not None
        adv = psi.advance(new.interior, new.boundary, new.radial_characteristic)
        assert adv.radial_characteristic is not None
        np.testing.assert_array_equal(
            adv.radial_characteristic.values, new.radial_characteristic.values,
        )
        lag1 = adv.at_lag(1)
        assert lag1.radial_characteristic is not None
        np.testing.assert_array_equal(
            lag1.radial_characteristic.values, psi.radial_characteristic.values,
        )
        with pytest.raises(TypeError, match="presence must match"):
            psi.advance(new.interior, new.boundary)  # dropping the block

    def test_copy_owns_the_seed_block(self):
        sn = _sphere()
        x = _seeded_flux_composite(sn, 41)
        assert x.radial_characteristic is not None
        c = x.copy()
        assert c.radial_characteristic is not None
        np.testing.assert_array_equal(
            c.radial_characteristic.values, x.radial_characteristic.values,
        )
        if np.shares_memory(c.radial_characteristic.values, x.radial_characteristic.values):
            pytest.fail("copy must own its seed ndarray (no aliasing)")

    def test_cross_mesh_seed_arithmetic_rejected(self):
        """Mesh-bound discipline on the new leaf (the BulkField guard law)."""
        a = RadialCharacteristicFlux.zeros_on(_sphere())
        b = RadialCharacteristicFlux.zeros_on(_sphere())
        with pytest.raises(ValueError, match="distinct SNMesh"):
            a - b


# ═══════════════════════════════════════════════════════════════════════
# A4 — the ghost metric: honest-scope gates (d1-landable part)
# ═══════════════════════════════════════════════════════════════════════


class TestA4SeedStateMetricVcell:
    r"""A4 (INVERTED) — the ψ½ STATE metric is ``G_sd = V_cell`` (SPD), NOT
    the ghost zero.  The starting-direction ray is a first-class radial STATE
    field (its self-block ``A_ss`` is a banded radial transport operator, not
    a grazing angular face — diag_gsd_05), so its Hilbert metric is the radial
    cell volume, mirroring the bulk ``G_bulk = V_cell·w_n``.  The ghost
    ``G_sd = 0`` was the single FORBIDDEN value: it puts the seed in null(G),
    severing the seed→bulk ``A_bs`` coupling from ``A.H`` — a WRONG adjoint
    for any nonzero seed (production defect 1.3e-2; derivation-of-record +
    diag_gsd_02/03).  These four gates were the green-but-blind ghost gates;
    inverted to the V_cell reality.

    Object-level (vv Mode-12): assert on the weight ARRAY and the scaled
    block values directly — NEVER via a downstream scalar (a scalar
    functional can lie in the metric's invariance group).
    """

    def test_seed_weights_are_v_cell_spd(self):
        r"""The state metric IS ``V_cell``: every ``(level, sign)`` leg's
        cells == ``V_cell`` (group-broadcast), corner == ``V[-1]`` (the gauge
        slot); strictly positive (SPD — the ghost 0 is the forbidden value)."""
        sn = _sphere()
        space = sn.radial_characteristic_space
        assert space is not None
        w = space.inner_product_weights
        if w is None:
            pytest.fail("the state metric must be EXPLICIT V_cell weights, "
                        "not None (None selects the Euclidean inner product)")
        V = np.asarray(sn.volumes, dtype=float).ravel()
        for p in space.levels:
            for sign in (-1, +1):
                np.testing.assert_array_equal(
                    space.cells_view(w, p, sign),
                    np.broadcast_to(V, (space.ng, space.nx)),
                )
                np.testing.assert_array_equal(space.corner_view(w, p, sign), V[-1])
        if not np.all(w > 0.0):
            pytest.fail(f"state metric must be SPD (strictly positive); got "
                        f"min {float(w.min()):.3e} — the forbidden ghost value")

    def test_apply_metric_scales_seed_by_v_cell_and_leaves_bulk_trace(self):
        r"""``G⊙x`` SCALES the ψ½ block by ``V_cell`` (not zero); bulk/trace
        stay BIT-IDENTICAL to the unseeded twin's (the metric is
        block-diagonal — the seed block changes nothing outside itself)."""
        sn = _sphere()
        ffs = sn.full_field_space
        x = _seeded_flux_composite(sn, 51)
        twin = FullField(interior=x.interior, boundary=x.boundary)
        gx = ffs.apply_metric(x)
        g_twin = ffs.apply_metric(twin)
        assert gx.radial_characteristic is not None
        w = sn.radial_characteristic_space.inner_product_weights
        np.testing.assert_array_equal(
            gx.radial_characteristic.values, w * x.radial_characteristic.values,
        )
        if not np.any(gx.radial_characteristic.values != 0.0):
            pytest.fail("G⊙x zeroed the seed block — the ghost metric lives")
        np.testing.assert_array_equal(gx.interior.values, g_twin.interior.values)
        np.testing.assert_array_equal(gx.boundary.values, g_twin.boundary.values)

    def test_apply_inverse_metric_divides_seed_by_v_cell(self):
        r"""``G⁺⊙x`` DIVIDES the ψ½ block by ``V_cell`` (no masking — the SPD
        metric has empty null space); the round-trip ``G⁺(G⊙x) = x`` holds on
        the seed (inverted from the ghost's masked-zero whole-block)."""
        sn = _sphere()
        ffs = sn.full_field_space
        x = _seeded_flux_composite(sn, 52)
        w = sn.radial_characteristic_space.inner_product_weights
        ginv_x = ffs.apply_inverse_metric(x)
        assert ginv_x.radial_characteristic is not None
        np.testing.assert_array_equal(
            ginv_x.radial_characteristic.values, x.radial_characteristic.values / w,
        )
        if not np.all(np.isfinite(ginv_x.to_flat())):
            pytest.fail("inverse metric produced non-finite values")
        round_trip = ffs.apply_inverse_metric(ffs.apply_metric(x))
        assert round_trip.radial_characteristic is not None
        np.testing.assert_allclose(
            round_trip.radial_characteristic.values,
            x.radial_characteristic.values, rtol=0.0, atol=1e-14,
        )

    def test_inner_product_seed_contribution_is_v_cell_weighted(self):
        r"""``⟨x,y⟩_G`` now CHANGES with the seed values (inverted from the
        ghost's seed-blindness): the seed term = ``Σ V_cell·x_seed·y_seed``,
        NONZERO for a nonzero seed — the Mode-12 CLOSURE mechanism (the seed
        rows carry metric weight, so a seed-transpose error is visible)."""
        sn = _sphere()
        ffs = sn.full_field_space
        x = _seeded_flux_composite(sn, 53)
        y = _seeded_flux_composite(sn, 54)
        x_twin = FullField(interior=x.interior, boundary=x.boundary)
        y_twin = FullField(interior=y.interior, boundary=y.boundary)
        w = np.asarray(
            sn.radial_characteristic_space.inner_product_weights, dtype=float,
        )
        assert x.radial_characteristic is not None and y.radial_characteristic is not None
        seed_term = float(np.sum(
            w * x.radial_characteristic.values * y.radial_characteristic.values
        ))
        full = ffs.inner_product(x, y)
        bulk_trace = ffs.inner_product(x_twin, y_twin)
        # the composite routes EXACTLY the V_cell-weighted seed term into
        # the total (inverted from "seed contributes 0"):
        np.testing.assert_allclose(
            full - bulk_trace, seed_term, rtol=1e-11, atol=1e-12,
        )
        if abs(seed_term) < 1e-9:
            pytest.fail(f"seed inner-product term ~0 ({seed_term:.2e}) — the "
                        f"ghost metric lives (expected Σ V_cell·x·y ≠ 0)")

    def test_inner_product_mixed_presence_raises(self):
        sn = _sphere()
        ffs = sn.full_field_space
        x = _seeded_flux_composite(sn, 55)
        y_twin = FullField(interior=x.interior, boundary=x.boundary)
        with pytest.raises(ValueError, match="mixed starting-direction"):
            ffs.inner_product(x, y_twin)

    def test_seeded_field_on_seedless_space_rejected(self):
        """The illegal pairing: a ψ½-carrying composite through a composite
        space with no seed block space (the transitional discipline's one
        forbidden quadrant)."""
        sn = _sphere()
        slab = _slab()
        x = _seeded_flux_composite(sn, 56)
        with pytest.raises(RuntimeError, match="no.*radial_characteristic_space"):
            slab.full_field_space.apply_metric(x)
