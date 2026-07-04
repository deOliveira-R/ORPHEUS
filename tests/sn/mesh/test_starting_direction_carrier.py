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
* **A4** — the ghost-metric HONEST-SCOPE gates: all seed weights are
  exactly zero, ``apply_metric`` / ``apply_inverse_metric`` zero the
  block while leaving bulk/trace bit-identical to the unseeded twin,
  and ``inner_product`` receives EXACTLY zero seed contribution.

.. note:: **A4's positive Mode-12 pin is deferred to d3** (gate spec
   §16.A A4 final leg): the "sign-flip the seed-row transpose → G3
   reciprocity STAYS GREEN while the C(i) round-trip REDS" asymmetry
   proof needs a transpose that READS the block — that lands with the
   d3 walk triple. Until then the blindness is documented here and on
   :mod:`orpheus.numerics.spaces.starting_direction_space`: the zero
   G-weight puts the seed rows inside the reciprocity functional's
   invariance group — G3 can NEVER be credited for seed-row
   correctness (the catchers are §16.B B1, §16.C C(i), and 2.5b's
   Euclidean :math:`M^{\mathsf T}` oracle).

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
from orpheus.numerics.spaces.starting_direction_space import StartingDirectionSpace
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.fields.starting_direction_flux import StartingDirectionFlux
from orpheus.transport.displacements import StartingDirectionDisplacement
from orpheus.transport.full_field import FullField
from orpheus.transport.source_sinks import (
    AngularBoundarySourceSink,
    AngularSourceSink,
    StartingDirectionSourceSink,
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
        bulk=AngularFlux,
        boundary=AngularBoundaryFlux,
        mesh=sn,
        starting_direction=StartingDirectionFlux,
    )
    rng = np.random.default_rng(seed)
    assert z.starting_direction is not None  # builder precondition, not a gate
    return z._recombine(
        bulk=replace(z.bulk, values=rng.normal(size=z.bulk.values.shape)),
        boundary=replace(
            z.boundary, values=rng.normal(size=z.boundary.values.shape),
        ),
        starting_direction=replace(
            z.starting_direction,
            values=rng.normal(size=z.starting_direction.values.shape),
        ),
    )


def _seeded_source_composite(sn: SNMesh, seed: int) -> FullField:
    """A random seed-carrying SOURCE composite (additive role algebra)."""
    z = FullField.zeros(
        bulk=AngularSourceSink,
        boundary=AngularBoundarySourceSink,
        mesh=sn,
        starting_direction=StartingDirectionSourceSink,
    )
    rng = np.random.default_rng(seed)
    assert z.starting_direction is not None
    return z._recombine(
        bulk=replace(z.bulk, values=rng.normal(size=z.bulk.values.shape)),
        boundary=replace(
            z.boundary, values=rng.normal(size=z.boundary.values.shape),
        ),
        starting_direction=replace(
            z.starting_direction,
            values=rng.normal(size=z.starting_direction.values.shape),
        ),
    )


# ═══════════════════════════════════════════════════════════════════════
# A1 — R12a presence, pinned BOTH ways
# ═══════════════════════════════════════════════════════════════════════


class TestA1PresenceR12a:
    def test_sphere_carries_one_seed_level(self):
        sn = _sphere()
        np.testing.assert_array_equal(sn.starting_direction_levels, (0,))
        space = sn.starting_direction_space
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
        np.testing.assert_array_equal(sn.starting_direction_levels, ())
        if sn.starting_direction_space is not None:
            pytest.fail(
                f"{builder.__name__}: non-carrying mesh must have "
                f"starting_direction_space None (R12a)"
            )
        z = TimedFullField.zeros(
            bulk=AngularFlux,
            boundary=AngularBoundaryFlux,
            mesh=sn,
            starting_direction=StartingDirectionFlux,  # passed UNIFORMLY
        )
        if z.starting_direction is not None:
            pytest.fail(
                "composite factory must omit the ψ½ block on a "
                "non-carrying mesh even when the leaf class is passed"
            )
        with pytest.raises(ValueError, match="no starting-direction"):
            StartingDirectionFlux.zeros_on(sn)

    def test_sphere_composite_zeros_carry_typed_block(self):
        sn = _sphere()
        z = TimedFullField.zeros(
            bulk=AngularFlux,
            boundary=AngularBoundaryFlux,
            mesh=sn,
            starting_direction=StartingDirectionFlux,
        )
        sd = z.starting_direction
        if not isinstance(sd, StartingDirectionFlux):
            pytest.fail(f"expected StartingDirectionFlux, got {type(sd)}")
        np.testing.assert_array_equal(sd.values, 0.0)
        np.testing.assert_array_equal(sd.cells(0, -1).shape, (sn.ng, 4))
        np.testing.assert_array_equal(sd.cells(0, +1).shape, (sn.ng, 4))
        np.testing.assert_array_equal(sd.corner(0, -1).shape, (sn.ng,))
        np.testing.assert_array_equal(sd.corner(0, +1).shape, (sn.ng,))

    def test_full_field_space_third_block_keyed_by_presence(self):
        """The composite SPACE grows by exactly the seed size on the sphere
        and is unchanged on non-carrying meshes."""
        sn = _sphere()
        seed_space = sn.starting_direction_space
        assert seed_space is not None  # narrowing (A1 sphere gate above pins it)
        n_bulk = sn.quad.N * sn.ng * 4
        n_trace = sn.angular_trace.shape[0]
        np.testing.assert_array_equal(
            sn.full_field_space.shape,
            (n_bulk + n_trace + seed_space.shape[0],),
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
        space = sn.starting_direction_space
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
        space = sn.starting_direction_space
        assert space is not None
        buf = np.zeros(space.shape)
        with pytest.raises(KeyError, match="no starting-direction block"):
            space.cells_view(buf, 1, -1)
        with pytest.raises(ValueError, match="sign must be"):
            space.corner_view(buf, 0, 0)

    def test_empty_levels_space_unrepresentable(self):
        """Absence is spelled None, never a zero-DOF phantom space."""
        with pytest.raises(ValueError, match="levels is empty"):
            StartingDirectionSpace.for_levels((), ng=2, nx=4)


# ═══════════════════════════════════════════════════════════════════════
# A2 — flat round-trip including the seed tail (+ the mutation teeth)
# ═══════════════════════════════════════════════════════════════════════


class TestA2FlatRoundTrip:
    def test_round_trip_identity_and_length(self):
        sn = _sphere()
        x = _seeded_flux_composite(sn)
        assert x.starting_direction is not None
        flat = x.to_flat()
        expected_len = (
            x.bulk.values.size
            + x.boundary.values.size
            + x.starting_direction.values.size
        )
        np.testing.assert_array_equal(flat.size, expected_len)
        y = FullField.from_flat(flat, x)
        np.testing.assert_array_equal(y.bulk.values, x.bulk.values)
        np.testing.assert_array_equal(y.boundary.values, x.boundary.values)
        assert y.starting_direction is not None
        np.testing.assert_array_equal(
            y.starting_direction.values, x.starting_direction.values,
        )
        if type(y) is not TimedFullField:
            pytest.fail("template type must be preserved through from_flat")

    def test_seed_tail_layout_is_the_trailing_block(self):
        """The flat layout is [bulk, trace, seed] — pin the tail slice."""
        sn = _sphere()
        x = _seeded_flux_composite(sn)
        assert x.starting_direction is not None
        flat = x.to_flat()
        n_head = x.bulk.values.size + x.boundary.values.size
        np.testing.assert_array_equal(
            flat[n_head:], x.starting_direction.values,
        )

    def test_mutation_teeth_dropped_seed_slice_reds_both_legs(self, monkeypatch):
        """MUTATION (in-process, per the mutation-hygiene rule): a to_flat
        that silently truncates the ψ½ block breaks (a) the length pin and
        (b) from_flat's size guard — proving the block is genuinely
        serialized, not decorative."""
        sn = _sphere()
        x = _seeded_flux_composite(sn)
        assert x.starting_direction is not None

        def _two_block_to_flat(self):  # the pre-2.5d serializer
            return np.concatenate([
                self.bulk.values.ravel(), self.boundary.values,
            ])

        monkeypatch.setattr(FullField, "to_flat", _two_block_to_flat)
        mutated = x.to_flat()
        expected_len = (
            x.bulk.values.size
            + x.boundary.values.size
            + x.starting_direction.values.size
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
        assert a.starting_direction is not None
        assert b.starting_direction is not None
        av, bv = a.starting_direction.values, b.starting_direction.values
        s = a + b
        assert s.starting_direction is not None
        np.testing.assert_array_equal(s.starting_direction.values, av + bv)
        d = a - b
        assert d.starting_direction is not None
        np.testing.assert_array_equal(d.starting_direction.values, av - bv)
        n = -a
        assert n.starting_direction is not None
        np.testing.assert_array_equal(n.starting_direction.values, -av)
        m = 2.5 * a
        assert m.starting_direction is not None
        np.testing.assert_array_equal(m.starting_direction.values, 2.5 * av)
        q = a / 4.0
        assert q.starting_direction is not None
        np.testing.assert_array_equal(q.starting_direction.values, av / 4.0)
        # role closure: the summed block keeps the SourceSink role
        if type(s.starting_direction) is not StartingDirectionSourceSink:
            pytest.fail("source algebra must stay in the SourceSink role")

    def test_flux_torsor_triple_threads_the_seed(self):
        sn = _sphere()
        psi1 = _seeded_flux_composite(sn, 11)
        psi2 = _seeded_flux_composite(sn, 12)
        assert psi1.starting_direction is not None
        assert psi2.starting_direction is not None
        # ψ ⊖ ψ → displacement composite (seed block mints its sibling)
        delta = psi2 - psi1
        if type(delta.starting_direction) is not StartingDirectionDisplacement:
            pytest.fail(
                f"flux ⊖ flux must mint StartingDirectionDisplacement; got "
                f"{type(delta.starting_direction).__name__}"
            )
        np.testing.assert_array_equal(
            delta.starting_direction.values,
            psi2.starting_direction.values - psi1.starting_direction.values,
        )
        # ψ ⊕ Δ → flux. The THREADING pin is bitwise (the block computes
        # exactly ψ₁ ⊕ Δ leaf-wise); the torsor LAW ψ₁ ⊕ (ψ₂ ⊖ ψ₁) = ψ₂
        # holds to FP round-off only (a + (b − a) ≠ b bitwise).
        recovered = psi1 + delta
        if type(recovered.starting_direction) is not StartingDirectionFlux:
            pytest.fail("torsor action must return the flux role")
        np.testing.assert_array_equal(
            recovered.starting_direction.values,
            psi1.starting_direction.values + delta.starting_direction.values,
        )
        np.testing.assert_allclose(
            recovered.starting_direction.values,
            psi2.starting_direction.values,
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
        unseeded = FullField(bulk=a.bulk, boundary=a.boundary)
        with pytest.raises(ValueError, match="MIXED starting-direction"):
            (a + unseeded) if op == "add" else (a - unseeded)
        with pytest.raises(ValueError, match="MIXED starting-direction"):
            (unseeded + a) if op == "add" else (unseeded - a)

    def test_advance_presence_law_and_snapshot(self):
        sn = _sphere()
        psi = _seeded_flux_composite(sn, 31)
        assert psi.starting_direction is not None
        new = _seeded_flux_composite(sn, 32)
        assert new.starting_direction is not None
        adv = psi.advance(new.bulk, new.boundary, new.starting_direction)
        assert adv.starting_direction is not None
        np.testing.assert_array_equal(
            adv.starting_direction.values, new.starting_direction.values,
        )
        lag1 = adv.at_lag(1)
        assert lag1.starting_direction is not None
        np.testing.assert_array_equal(
            lag1.starting_direction.values, psi.starting_direction.values,
        )
        with pytest.raises(TypeError, match="presence must match"):
            psi.advance(new.bulk, new.boundary)  # dropping the block

    def test_copy_owns_the_seed_block(self):
        sn = _sphere()
        x = _seeded_flux_composite(sn, 41)
        assert x.starting_direction is not None
        c = x.copy()
        assert c.starting_direction is not None
        np.testing.assert_array_equal(
            c.starting_direction.values, x.starting_direction.values,
        )
        if np.shares_memory(c.starting_direction.values, x.starting_direction.values):
            pytest.fail("copy must own its seed ndarray (no aliasing)")

    def test_cross_mesh_seed_arithmetic_rejected(self):
        """Mesh-bound discipline on the new leaf (the BulkField guard law)."""
        a = StartingDirectionFlux.zeros_on(_sphere())
        b = StartingDirectionFlux.zeros_on(_sphere())
        with pytest.raises(ValueError, match="distinct SNMesh"):
            a - b


# ═══════════════════════════════════════════════════════════════════════
# A4 — the ghost metric: honest-scope gates (d1-landable part)
# ═══════════════════════════════════════════════════════════════════════


class TestA4GhostMetricHonestScope:
    def test_all_seed_weights_exactly_zero(self):
        sn = _sphere()
        space = sn.starting_direction_space
        assert space is not None
        w = space.inner_product_weights
        if w is None:
            pytest.fail("the ghost metric must be EXPLICIT zeros, not None "
                        "(None selects the Euclidean inner product)")
        np.testing.assert_array_equal(w, 0.0)

    def test_apply_metric_zeroes_seed_and_matches_unseeded_twin(self):
        """G⊙x zeroes the ψ½ block; bulk/trace are BIT-IDENTICAL to the
        unseeded twin's — the seed changes nothing outside its block."""
        sn = _sphere()
        ffs = sn.full_field_space
        x = _seeded_flux_composite(sn, 51)
        twin = FullField(bulk=x.bulk, boundary=x.boundary)
        gx = ffs.apply_metric(x)
        g_twin = ffs.apply_metric(twin)
        assert gx.starting_direction is not None
        np.testing.assert_array_equal(gx.starting_direction.values, 0.0)
        np.testing.assert_array_equal(gx.bulk.values, g_twin.bulk.values)
        np.testing.assert_array_equal(gx.boundary.values, g_twin.boundary.values)

    def test_apply_inverse_metric_masked_zero_on_seed(self):
        """The Moore–Penrose masked inverse: the whole seed block is metric
        null space — zeroed, never divided."""
        sn = _sphere()
        ffs = sn.full_field_space
        x = _seeded_flux_composite(sn, 52)
        gx = ffs.apply_inverse_metric(x)
        assert gx.starting_direction is not None
        np.testing.assert_array_equal(gx.starting_direction.values, 0.0)
        if not np.all(np.isfinite(gx.to_flat())):
            pytest.fail("masked pseudo-inverse must never divide by the "
                        "zero seed weights")

    def test_inner_product_seed_contribution_exactly_zero(self):
        """⟨x, y⟩_G is UNCHANGED by arbitrary seed values — the Mode-12
        blindness mechanism, asserted as an exact identity."""
        sn = _sphere()
        ffs = sn.full_field_space
        x = _seeded_flux_composite(sn, 53)
        y = _seeded_flux_composite(sn, 54)
        x_twin = FullField(bulk=x.bulk, boundary=x.boundary)
        y_twin = FullField(bulk=y.bulk, boundary=y.boundary)
        np.testing.assert_array_equal(
            ffs.inner_product(x, y), ffs.inner_product(x_twin, y_twin),
        )

    def test_inner_product_mixed_presence_raises(self):
        sn = _sphere()
        ffs = sn.full_field_space
        x = _seeded_flux_composite(sn, 55)
        y_twin = FullField(bulk=x.bulk, boundary=x.boundary)
        with pytest.raises(ValueError, match="mixed starting-direction"):
            ffs.inner_product(x, y_twin)

    def test_seeded_field_on_seedless_space_rejected(self):
        """The illegal pairing: a ψ½-carrying composite through a composite
        space with no seed block space (the transitional discipline's one
        forbidden quadrant)."""
        sn = _sphere()
        slab = _slab()
        x = _seeded_flux_composite(sn, 56)
        with pytest.raises(RuntimeError, match="no.*starting_direction_space"):
            slab.full_field_space.apply_metric(x)
