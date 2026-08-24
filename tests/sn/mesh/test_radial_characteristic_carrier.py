r"""§16.A — the starting-direction PRESENCE + space intrinsic gates (post-B.2d).

The 2.5d-era §16.A carrier-law suite pinned the ψ½ block's life ON the
3-block augmented composite (presence-keyed ``zeros``, the mixed-presence
law, the flat seed tail, the seed metric arms of ``FullFieldSpace``). The
B.2d ray eviction dissolved that carrier: ``FullField`` is a pure 2-block
composite and System B is its own
:class:`~orpheus.transport.radial_characteristic_field.RadialCharacteristicField`,
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
  (``sn.radial_characteristic_field_space``), gated below against the
  same ``V_cell``/corner-gauge numbers.

Retirement ledger (B.2d d2 — where each retired class's claim went):

* **A2 (flat round-trip incl. the seed tail)** — the 3-block layout is
  unrepresentable; the 2-block round trip is the base composite law
  (``tests/transport/test_composite.py``) and the COUPLED round trip +
  coverage is G-c2.4 (``test_psi_half_coupling::TestCoupledBuilder``).
* **A3 (algebra closure threading the seed; the mixed-presence law; the
  ``advance`` presence law)** — the seed-threading arms are unspellable;
  System B's own 2-block algebra closure lives in
  ``tests/transport/test_radial_characteristic_field.py``; the
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
from orpheus.numerics.spaces.radial_characteristic_space import (
    RadialCharacteristicInteriorSpace,
)
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.full_field import FullField
from orpheus.transport.radial_characteristic_field import (
    RadialCharacteristicField,
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


def _cyl_folded() -> SNMesh:
    """Cylinder on the admitted folded family — CARRYING on every level.

    Until Q5.6.3 this file held two NON-carrying cylinder fixtures
    (`_cyl_level_symmetric`: the dead-seed τ_raw,0 = 1 discriminator;
    `_cyl_product`: the edge-node τ_raw,0 = 0 rank-duplicate).  The
    admission flip made both UNCONSTRUCTIBLE — their refusals, with the
    per-family reasons, are gated one tier up in
    ``test_cylindrical_quadrature_admission.py``; presence here covers
    only meshes that exist."""
    return _mesh_1d(
        CoordSystem.CYLINDRICAL, Quadrature.folded_product(n_mu=2, n_phi=4),
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
        space = sn.radial_characteristic_field_space
        if space is None:
            pytest.fail("sphere-GL mesh must carry a starting-direction space")
        # 1 level × 2 signs × (ng·nx cells + ng corner) — the composite space's
        # flat size (the successor of the retired unified space total, 4e).
        expected = 1 * 2 * (sn.ng * 4 + sn.ng)
        np.testing.assert_array_equal(space.shape, (expected,))
        np.testing.assert_array_equal(
            sn.radial_characteristic_interior_space.levels, (0,),
        )

    def test_sphere_constructs_system_b(self):
        """The presence FACT, re-homed (B.2d): a carrying mesh's System B
        EXISTS — the member composite constructs presence-gated, its
        composite member space is cached and 2-block-consistent."""
        sn = _sphere()
        ray = RadialCharacteristicField.flux_zeros(sn.radial_characteristic_field_space)
        np.testing.assert_array_equal(ray.to_flat(), 0.0)
        cspace = sn.radial_characteristic_field_space
        if cspace is None:
            pytest.fail("carrying mesh must cache the System-B member space")
        np.testing.assert_array_equal(
            ray.to_flat().size, int(np.prod(cspace.shape)),
        )

    def test_folded_cylinder_carries_every_level_and_constructs_system_b(self):
        """Q5.6.3 — the admitted cylinder is a MULTI-LEVEL carrier: every
        level position carries, the composite space exists, and the
        System-B member constructs (the cylinder twin of the two sphere
        presence gates above)."""
        sn = _cyl_folded()
        n_levels = len(sn.quad.level_indices)
        np.testing.assert_array_equal(
            sn.radial_characteristic_levels, tuple(range(n_levels)),
        )
        space = sn.radial_characteristic_field_space
        if space is None:
            pytest.fail("folded cylinder must carry a starting-direction space")
        ray = RadialCharacteristicField.flux_zeros(sn.radial_characteristic_field_space)
        np.testing.assert_array_equal(ray.to_flat(), 0.0)
        np.testing.assert_array_equal(
            ray.to_flat().size, int(np.prod(space.shape)),
        )

    @pytest.mark.parametrize(
        "builder",
        [_slab, cart2d_2g_nonsquare],
        ids=["slab", "cart2d"],
    )
    def test_non_carrying_meshes_have_no_system_b(self, builder):
        """Absent-on-non-carrying, pinned positively (a future leak REDS here).

        Since Q5.6.3 the Cartesian geometries are the ONLY admitted
        non-carrying meshes (the non-carrying cylinder rows this test
        held became unconstructible at the admission flip — their
        refusal reasons are pinned in
        ``test_cylindrical_quadrature_admission.py``). B.2d: absence is
        SYSTEM-B NON-EXISTENCE (space None, the member-space parse
        refuses — S5's composite seam — the composite unconstructable).
        """
        sn = builder()
        np.testing.assert_array_equal(sn.radial_characteristic_levels, ())
        if sn.radial_characteristic_field_space is not None:
            pytest.fail(
                f"{builder.__name__}: non-carrying mesh must have NO "
                f"System-B member space"
            )
        if sn.radial_characteristic_interior_space is not None:
            pytest.fail(
                f"{builder.__name__}: non-carrying mesh must have NO "
                f"interior split space (R12a)"
            )
        with pytest.raises(ValueError, match="System B is absent"):
            RadialCharacteristicField.flux_zeros(sn.radial_characteristic_field_space)

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
                interior=AngularFlux, boundary=AngularBoundaryFlux, space=sn.full_field_space,
                radial_characteristic=None,
            )
        with pytest.raises(TypeError):
            TimedFullField.zeros(
                interior=AngularFlux, boundary=AngularBoundaryFlux, space=sn.full_field_space,
                radial_characteristic=None,
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
        # 4e: the split spaces carry their own flat backings (cells in the
        # interior space, corners in the boundary space); slot_view is a slice
        # VIEW into each, so writes mutate the flat buffer.
        interior = sn.radial_characteristic_interior_space
        boundary = sn.radial_characteristic_boundary_space
        assert interior is not None and boundary is not None
        ibuf = np.zeros(interior.shape)
        bbuf = np.zeros(boundary.shape)
        interior.slot_view(ibuf, 0, -1)[:] = 1.0
        boundary.slot_view(bbuf, 0, +1)[:] = 2.0
        np.testing.assert_allclose(np.sum(ibuf == 1.0), interior.ng * interior.nx)
        np.testing.assert_allclose(np.sum(bbuf == 2.0), boundary.ng)
        # Interior layout order: (−1 cells, +1 cells) ⇒ the −1 cells are the
        # first slot; boundary: (−1 corner, +1 corner) ⇒ the +1 corner is last.
        np.testing.assert_array_equal(ibuf[: interior.ng * interior.nx], 1.0)
        np.testing.assert_array_equal(bbuf[boundary.ng :], 2.0)

    def test_unknown_level_and_sign_rejected(self):
        sn = _sphere()
        interior = sn.radial_characteristic_interior_space
        boundary = sn.radial_characteristic_boundary_space
        assert interior is not None and boundary is not None
        ibuf = np.zeros(interior.shape)
        bbuf = np.zeros(boundary.shape)
        with pytest.raises(KeyError, match="no starting-direction block"):
            interior.slot_view(ibuf, 1, -1)
        with pytest.raises(ValueError, match="sign must be"):
            boundary.slot_view(bbuf, 0, 0)

    def test_empty_levels_space_unrepresentable(self):
        """Absence is spelled None, never a zero-DOF phantom space.

        ``cell_volumes`` is now a REQUIRED kwarg (the state metric G_sd =
        V_cell), but irrelevant here — the empty-``levels`` guard fires
        first; a valid dummy keeps the call well-formed so the ValueError
        (not a TypeError from the missing kwarg) is what we pin.
        """
        with pytest.raises(ValueError, match="levels is empty"):
            RadialCharacteristicInteriorSpace.for_levels(
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
        interior = sn.radial_characteristic_interior_space
        boundary = sn.radial_characteristic_boundary_space
        assert interior is not None and boundary is not None
        iw = interior.inner_product_weights
        bw = boundary.inner_product_weights
        if iw is None or bw is None:
            pytest.fail("the state metric must be EXPLICIT V_cell weights, "
                        "not None (None selects the Euclidean inner product)")
        V = np.asarray(sn.volumes, dtype=float).ravel()
        for p in interior.levels:
            for sign in (-1, +1):
                np.testing.assert_array_equal(
                    interior.slot_view(iw, p, sign),
                    np.broadcast_to(V, (interior.ng, interior.nx)),
                )
                np.testing.assert_array_equal(boundary.slot_view(bw, p, sign), V[-1])
        if not (np.all(iw > 0.0) and np.all(bw > 0.0)):
            pytest.fail("state metric must be SPD (strictly positive) — the "
                        "forbidden ghost value")

    def test_composite_member_metric_is_the_v_cell_gauge(self):
        r"""System B's OWN composite-space inner product carries the
        ``Σ V_cell·x·y`` numbers — the member-wise realization of the state
        metric — checked against a HAND-BUILT ``V_cell`` gauge from
        ``sn.volumes`` (4e RE-POSE: the retired unified gauge is gone, so the
        reference is built directly from the cell volumes, structurally
        independent of the space's stored weights)."""
        sn = _sphere()
        cspace = sn.radial_characteristic_field_space
        interior = sn.radial_characteristic_interior_space
        boundary = sn.radial_characteristic_boundary_space
        assert cspace is not None and interior is not None and boundary is not None
        rng = np.random.default_rng(53)
        ns = int(cspace.shape[0])
        x = RadialCharacteristicField.from_flat(
            rng.standard_normal(ns), RadialCharacteristicField.flux_zeros(sn.radial_characteristic_field_space),
        )
        y = RadialCharacteristicField.from_flat(
            rng.standard_normal(ns), RadialCharacteristicField.flux_zeros(sn.radial_characteristic_field_space),
        )
        # Hand-built V_cell gauge (composite to_flat order = interior ⊕ boundary),
        # from sn.volumes directly — the structurally-independent reference.
        V = np.asarray(sn.volumes, dtype=float).ravel()
        iw = np.zeros(interior.shape[0])
        bw = np.zeros(boundary.shape[0])
        for p in interior.levels:
            for sign in (-1, +1):
                interior.slot_view(iw, p, sign)[:] = V[None, :]
                boundary.slot_view(bw, p, sign)[:] = V[-1]
        w = np.concatenate([iw, bw])
        hand_term = float(np.sum(w * x.to_flat() * y.to_flat()))
        composite_term = cspace.inner_product(x, y)
        np.testing.assert_allclose(
            composite_term, hand_term, rtol=1e-12, atol=1e-13,
        )
        if abs(hand_term) < 1e-9:
            pytest.fail(f"seed inner-product term ~0 ({hand_term:.2e}) — "
                        f"the metric gate is vacuous (expected Σ V_cell·x·y ≠ 0)")
