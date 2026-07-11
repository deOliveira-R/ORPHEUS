r"""§16.A — the starting-direction PRESENCE + space intrinsic gates (post-B.2d).

The 2.5d-era §16.A carrier-law suite pinned the ψ½ block's life ON the
3-block augmented composite (presence-keyed ``zeros``, the mixed-presence
law, the flat seed tail, the seed metric arms of ``FullFieldSpace``). The
B.2d ray eviction dissolved that carrier: ``FullField`` is a pure 2-block
composite and System B is its own
:class:`~orpheus.transport.radial_characteristic_composite.RadialCharacteristicComposite`,
so the LAW tests retired with the law and the surviving content here is the
**presence FACTS** (the R12a trichotomy, re-homed as System-B EXISTENCE) and
the **unified-space intrinsics** (layout views, the ``G_sd = V_cell`` state
metric) that the walk currency still rides:

* **A1** — R12a presence pinned BOTH ways, as System-B existence: the
  sphere-GL mesh carries the space (1 level) and CONSTRUCTS System B; the
  cylinder (BOTH production rules — the level-symmetric one is the R12a
  discriminator: μ_start ∉ its μ-nodes, yet τ_raw = 1.0 bit-exact makes the
  seed dead) and Cartesian (1-D and 2-D) MUST NOT — the space is ``None``,
  the leaf factory raises, System B is unconstructable, and the 2-block
  composite factories REJECT the retired ``radial_characteristic=`` kwarg
  outright (the G-d2.2 mesh-side tooth: an eviction that left a seed arm
  behind reds here).
* **A4** — the unified space's STATE metric ``G_sd = V_cell`` (SPD) — the
  space-level fact, unchanged by the eviction. Its COMPOSITE realization
  (the member-wise metric dispatch) lives on System B's own composite space
  (``sn.radial_characteristic_composite_space``), gated below against the
  same ``V_cell``/corner-gauge numbers.

Retirement ledger (B.2d d2 — where each retired class's claim went):

* **A2 (flat round-trip incl. the seed tail)** — the 3-block layout is
  unrepresentable; the 2-block round trip is the base composite law
  (``tests/transport/test_composite.py``) and the COUPLED round trip +
  coverage is G-c2.4 (``test_psi_half_coupling::TestCoupledBuilder``).
* **A3 (algebra closure threading the seed; the mixed-presence law; the
  ``advance`` presence law)** — the seed-threading arms are unspellable;
  System B's own 2-block algebra closure lives in
  ``tests/transport/test_radial_characteristic_composite.py``; the
  mixed-presence and advance presence LAWS retired with the block (their
  raise sites are gone — presence is System-B existence, G-c2.2).
* **A4's FullFieldSpace seed arms** (apply_metric/apply_inverse_metric/
  inner_product seed terms + the mixed-presence raise + the illegal
  pairing) — the space is 2-block; System B's metric rides its OWN
  composite-space instance (gated here) and the coupled G-reciprocity is
  G-c2.5.

vv Mode-8 discipline: ``np.testing.assert_*`` / ``pytest.raises`` only
(no bare ``assert`` on numerical claims — the canonical invocation is
``python -O``).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.spaces.radial_characteristic_space import RadialCharacteristicSpace
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.fields.radial_characteristic_flux import RadialCharacteristicFlux
from orpheus.transport.full_field import FullField
from orpheus.transport.radial_characteristic_composite import (
    RadialCharacteristicComposite,
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


# ═══════════════════════════════════════════════════════════════════════
# A1 — R12a presence, pinned BOTH ways (as System-B existence, B.2d)
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

    def test_sphere_constructs_system_b(self):
        """The presence FACT, re-homed (B.2d): a carrying mesh's System B
        EXISTS — the member composite constructs presence-gated, its
        composite member space is cached and 2-block-consistent."""
        sn = _sphere()
        ray = RadialCharacteristicComposite.from_mesh(sn)
        np.testing.assert_array_equal(ray.to_flat(), 0.0)
        cspace = sn.radial_characteristic_composite_space
        if cspace is None:
            pytest.fail("carrying mesh must cache the System-B member space")
        np.testing.assert_array_equal(
            ray.to_flat().size, int(np.prod(cspace.shape)),
        )

    @pytest.mark.parametrize(
        "builder",
        [_cyl_level_symmetric, _cyl_product, _slab, cart2d_2g_nonsquare],
        ids=["cyl_level_symmetric", "cyl_product", "slab", "cart2d"],
    )
    def test_non_carrying_meshes_have_no_system_b(self, builder):
        """Absent-on-non-carrying, pinned positively (a future leak REDS here).

        The cylinder-LS row is the R12a discriminator: the R12 letter
        (μ_start ∉ μ-nodes) fires on it, so THIS row is what pins the
        refined predicate — a re-spelling back to node membership grows
        dead seed blocks on every level-symmetric cylinder and reds
        this gate. B.2d: absence is SYSTEM-B NON-EXISTENCE (space None,
        leaf factory raises, the composite unconstructable).
        """
        sn = builder()
        np.testing.assert_array_equal(sn.radial_characteristic_levels, ())
        if sn.radial_characteristic_space is not None:
            pytest.fail(
                f"{builder.__name__}: non-carrying mesh must have "
                f"radial_characteristic_space None (R12a)"
            )
        if sn.radial_characteristic_composite_space is not None:
            pytest.fail(
                f"{builder.__name__}: non-carrying mesh must have NO "
                f"System-B member space"
            )
        with pytest.raises(ValueError, match="no starting-direction"):
            RadialCharacteristicFlux.zeros_on(sn)
        with pytest.raises(
            ValueError, match="no radial_characteristic_interior_space",
        ):
            RadialCharacteristicComposite.from_mesh(sn)

    @pytest.mark.parametrize(
        "builder", [_sphere, _slab], ids=["sphere", "slab"],
    )
    def test_composite_factories_reject_the_retired_seed_kwarg(self, builder):
        """G-d2.2 (mesh side): the 2-block factories REJECT a
        ``radial_characteristic=`` kwarg outright — on EVERY mesh, carrying
        included (the eviction left no seed arm; a revived arm reds here)."""
        sn = builder()
        with pytest.raises(TypeError):
            FullField.zeros(
                interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn,
                radial_characteristic=RadialCharacteristicFlux,
            )
        with pytest.raises(TypeError):
            TimedFullField.zeros(
                interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn,
                radial_characteristic=RadialCharacteristicFlux,
            )

    def test_full_field_space_is_two_block_everywhere(self):
        """G-d2.2 (space side): System A's composite space is the honest
        bulk ⊕ trace sum on EVERY mesh — the sphere's seed DOFs live on
        System B's OWN member space, never as a third block here."""
        for sn in (_sphere(), _slab()):
            n_interior = sn.quad.N * sn.ng * 4
            n_trace = sn.angular_trace.shape[0]
            np.testing.assert_array_equal(
                sn.full_field_space.shape, (n_interior + n_trace,),
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
# A4 — the ψ½ STATE metric G_sd = V_cell (SPD): the space-level fact +
# its System-B composite realization
# ═══════════════════════════════════════════════════════════════════════


class TestA4SeedStateMetricVcell:
    r"""A4 — the ψ½ STATE metric is ``G_sd = V_cell`` (SPD), NOT the ghost
    zero.  The starting-direction ray is a first-class radial STATE field
    (its self-block ``A_BB`` is a banded radial transport operator, not a
    grazing angular face — diag_gsd_05), so its Hilbert metric is the radial
    cell volume, mirroring the bulk ``G_bulk = V_cell·w_n``.  The ghost
    ``G_sd = 0`` was the single FORBIDDEN value: it puts the seed in null(G),
    severing the seed→bulk coupling from ``A.H`` — a WRONG adjoint for any
    nonzero seed (production defect 1.3e-2; derivation-of-record +
    diag_gsd_02/03).

    B.2d: the metric's COMPOSITE realization rides System B's OWN member
    space (the member-wise dispatch of the family-blind composite-space
    class); the coupled G-reciprocity it feeds is G-c2.5.

    Object-level (vv Mode-12): assert on the weight ARRAY directly — NEVER
    via a downstream scalar (a scalar functional can lie in the metric's
    invariance group).
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

    def test_composite_member_metric_matches_the_unified_gauge(self):
        r"""System B's OWN composite-space inner product carries the SAME
        ``Σ V_cell·x·y`` numbers as the unified gauge — the member-wise
        realization of the state metric (the B.2d home of the retired
        ``FullFieldSpace`` seed arms; the split↔unified re-label is exact,
        so the two spellings must agree to reassociation)."""
        sn = _sphere()
        cspace = sn.radial_characteristic_composite_space
        uspace = sn.radial_characteristic_space
        assert cspace is not None and uspace is not None
        rng = np.random.default_rng(53)
        x_u = RadialCharacteristicFlux.zeros_on(sn)
        y_u = RadialCharacteristicFlux.zeros_on(sn)
        x_u.values[...] = rng.standard_normal(x_u.values.shape)
        y_u.values[...] = rng.standard_normal(y_u.values.shape)
        w = np.asarray(uspace.inner_product_weights, dtype=float)
        unified_term = float(np.sum(w * x_u.values * y_u.values))
        composite_term = cspace.inner_product(
            RadialCharacteristicComposite.from_unified(x_u),
            RadialCharacteristicComposite.from_unified(y_u),
        )
        np.testing.assert_allclose(
            composite_term, unified_term, rtol=1e-12, atol=1e-13,
        )
        if abs(unified_term) < 1e-9:
            pytest.fail(f"seed inner-product term ~0 ({unified_term:.2e}) — "
                        f"the metric gate is vacuous (expected Σ V_cell·x·y ≠ 0)")
