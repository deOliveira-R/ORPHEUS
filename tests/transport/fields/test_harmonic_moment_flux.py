r"""Member contracts of :class:`HarmonicMomentFlux` (#399 / CS4b S6.1 — G6.7).

The derived-view members — :meth:`truncate`, :meth:`isotropic_part`,
:meth:`anisotropic_part`, :meth:`l_block` — plus :meth:`scalar_flux`'s
self-derive contract, gated at BOTH widths (``spatial_moments`` 1 and 2).

#399's finding (census 2026-08-21) was that the members were
spatial-moment-blind. Re-measured at HEAD 2026-08-24 before designing
(the plan-authoring shelf-life check): the parts were ALREADY
tail-correct — the S4 space-passing rework repaired them, ungated — and
``truncate`` carried a loud widened ``NotImplementedError`` defer, not
the census's latent broadcast crash. S6.1 lifted the defer with a
space-derived rebuild (swap the spherical-harmonic head factor, keep
every remaining factor verbatim) and this module is the first gate
coverage the members have at ANY width — `[M]` 2026-08-24, zero callers
of all four members in orpheus/ or tests/ before it.

Mutation record (G6.7's row, `[M]` in-process 2026-08-24): restoring the
un-tailed ``new_shape = (L+1, 2L+1, *spatial-only)`` in ``truncate``
reddens the widened truncation row; dropping the tail factor from the
rebuilt space (the pre-S6.1 ``factors[1]``-only rebuild) reddens the
space-content row; zeroing the partition's ℓ=0 copy reddens the
partition law.

vv Mode-8 discipline: assertions are function calls (``np.testing.*`` /
``pytest.raises`` / ``pytest.fail``) — canonical invocation is
``python -O``.
"""

from __future__ import annotations

import numpy as np
import numpy.testing as npt
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.transport.fields.harmonic_moment_flux import HarmonicMomentFlux
from tests.sn._test_helpers import placeholder_materials

pytestmark = [pytest.mark.foundation]

_L, _NG, _NX = 2, 2, 5


def _sn() -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, _NX + 1),
        mat_ids=np.zeros(_NX, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"), bc_right=BC("vacuum"),
    )
    return SNMesh(mesh, Quadrature.gauss_legendre(4), placeholder_materials(ng=_NG))


def _field(sn: SNMesh, spatial_moments: int, seed: int) -> HarmonicMomentFlux:
    tail = (spatial_moments,) if spatial_moments > 1 else ()
    values = np.random.default_rng(seed).standard_normal(
        (_L + 1, 2 * _L + 1, _NG, _NX, *tail)
    )
    return HarmonicMomentFlux.from_mesh_and_L(
        values, sn, _L, spatial_moments=spatial_moments,
    )


# ── truncate — the space-derived rebuild (both widths) ───────────────


class TestTruncate:
    @pytest.mark.parametrize("sm", [1, 2])
    def test_truncate_matches_the_factory_mint_at_both_widths(self, sm):
        """The truncated field's space content-equals the factory's OWN
        mint at (mesh, L_new, spatial_moments) — the single-source
        done-when of the space-derived rebuild. The sm=2 row is #399's
        FLIPPED witness (pre-S6.1 it raised the widened defer; pre-S4 it
        was the census's broadcast ValueError). Reddened by dropping the
        tail factor from the rebuild (the pre-S6.1 factors[1]-only
        spelling) or by re-introducing the un-tailed new_shape."""
        sn = _sn()
        f = _field(sn, sm, seed=20 + sm)
        L_new = 1
        g = f.truncate(L_new)
        tail = (sm,) if sm > 1 else ()
        expected_shape = (L_new + 1, 2 * L_new + 1, _NG, _NX, *tail)
        if g.values.shape != expected_shape:
            pytest.fail(
                f"truncate(sm={sm}) shape {g.values.shape} != "
                f"{expected_shape}"
            )
        factory_space = HarmonicMomentFlux.from_mesh_and_L(
            np.zeros(expected_shape), sn, L_new, spatial_moments=sm,
        ).space
        if g.space != factory_space:
            pytest.fail(
                f"truncated space must content-equal the factory mint at "
                f"(mesh, L_new={L_new}, sm={sm}); got {g.space!r}"
            )
        if g.L != L_new or g.spatial_moments != sm:
            pytest.fail(
                f"truncate must thread L and spatial_moments: got "
                f"L={g.L}, spatial_moments={g.spatial_moments}"
            )
        npt.assert_array_equal(
            g.values[: L_new + 1, : 2 * L_new + 1],
            f.values[: L_new + 1, : 2 * L_new + 1],
        )

    def test_truncate_widened_reads_the_width_off_its_space(self):
        """The truncated widened field's SPACE carries the moment factor
        (the single source): spatial_moments_per_axis reads 2 off it."""
        g = _field(_sn(), 2, seed=23).truncate(0)
        if g.spatial_moments_per_axis != 2:
            pytest.fail(
                f"the truncated space must carry the moment factor; "
                f"per-axis width read {g.spatial_moments_per_axis}"
            )

    @pytest.mark.parametrize("sm", [1, 2])
    def test_truncate_at_full_order_is_the_identity(self, sm):
        """truncate(self.L) reproduces the field bit-exactly, same-space
        content — the no-op leg of the truncation contract."""
        f = _field(_sn(), sm, seed=25 + sm)
        g = f.truncate(f.L)
        npt.assert_array_equal(g.values, f.values)
        if g.space != f.space:
            pytest.fail("full-order truncation must keep the space content")

    def test_truncate_range_refusals(self):
        f = _field(_sn(), 1, seed=27)
        with pytest.raises(ValueError, match="L_new=3 >"):
            f.truncate(_L + 1)
        with pytest.raises(ValueError, match="< 0"):
            f.truncate(-1)


# ── the ℓ-partition (both widths) ────────────────────────────────────


class TestPartition:
    @pytest.mark.parametrize("sm", [1, 2])
    def test_parts_partition_bit_exactly(self, sm):
        """iso + aniso == self, values bit-exact (the partition law the
        docstring names), at both widths — the sm=2 legs are #399's
        other two flipped witnesses (already repaired by S4, first
        GATED here). Reddened by zeroing the partition's ℓ=0 copy."""
        f = _field(_sn(), sm, seed=30 + sm)
        iso, aniso = f.isotropic_part(), f.anisotropic_part()
        npt.assert_array_equal(iso.values + aniso.values, f.values)
        if np.any(iso.values[1:]) or np.any(iso.values[0, 1:]):
            pytest.fail("isotropic_part must zero every ℓ ≥ 1 / m ≠ 0 slot")
        if np.any(aniso.values[0, 0]):
            pytest.fail("anisotropic_part must zero the (ℓ=0, m=0) slot")

    @pytest.mark.parametrize("sm", [1, 2])
    def test_parts_share_the_space_instance(self, sm):
        """The parts are replace-style derivations: SAME space instance
        (not a re-mint), same L, same width — the space-derived
        contract's cheapest observable."""
        f = _field(_sn(), sm, seed=33 + sm)
        for part in (f.isotropic_part(), f.anisotropic_part()):
            if part.space is not f.space:
                pytest.fail("a part must carry the parent's space instance")
            if part.L != f.L or part.spatial_moments != sm:
                pytest.fail("a part must thread L and spatial_moments")


# ── l_block (view) + scalar_flux (self-derive contract) ──────────────


class TestViews:
    def test_l_block_shapes_carry_the_tail(self):
        """l_block(l) is a raw view of the legitimate m-entries — shape
        (2l+1, ng, nx, tail) on a widened field, tail riding."""
        f = _field(_sn(), 2, seed=36)
        for l in range(_L + 1):
            blk = f.l_block(l)
            if blk.shape != (2 * l + 1, _NG, _NX, 2):
                pytest.fail(f"l_block({l}) shape {blk.shape} wrong")
        with pytest.raises(ValueError, match="out of range"):
            f.l_block(_L + 1)

    def test_scalar_flux_width1_self_derives_the_bulk_factor(self):
        """Width-1 scalar_flux derives its target from the product's
        cell-group factor — the carrier's shared bulk mint."""
        f = _field(_sn(), 1, seed=37)
        s = f.scalar_flux()
        npt.assert_array_equal(s.values, f.values[0, 0])
        if s.space != f.space.factors[1]:  # type: ignore[union-attr]
            pytest.fail("scalar target must be the product's bulk factor")

    def test_scalar_flux_widened_self_derive_refuses(self):
        """The widened self-derive is REFUSED by contract (S4): the
        widened target carries the scheme's mass-bearing moment axis,
        which the densified SpatialMomentSpace factor does not hold —
        the caller holding the pose passes space=. NOT a #399 defect."""
        f = _field(_sn(), 2, seed=38)
        with pytest.raises(TypeError, match="cannot self-derive"):
            f.scalar_flux()
