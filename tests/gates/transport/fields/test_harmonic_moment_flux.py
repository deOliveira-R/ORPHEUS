r"""Member contracts of :class:`HarmonicMomentFlux` (#399 / CS4b S6.1 — G6.7).

The :meth:`scalar_flux` self-derive contract, gated at BOTH widths (``spatial_moments`` 1 and 2). (The derived-view
members ``truncate``, ``isotropic_part``, ``anisotropic_part`` and ``l_block``, which this module first gated, were
retired with their zero production callers, 2026-11-20.)

vv Mode-8 discipline: assertions are function calls (``np.testing.*`` /
``pytest.raises`` / ``pytest.fail``) — canonical invocation is
``python -O``.

Layer DAG Phase 4 (2026-11-01): posed on a material mesh with the frame's moment_space_on.
"""

from __future__ import annotations

import numpy as np
import numpy.testing as npt
import pytest

from orpheus.geometry import StructuredGeometry
from orpheus.mesh import CellsByCount, Mesher
from orpheus.numerics.quadrature import Quadrature

from tests.gates.transport._phase_space_cases import Posed, moment_space
from orpheus.transport.mesh.material_mesh import MaterialMesh
from orpheus.transport.fields.harmonic_moment_flux import HarmonicMomentFlux
from orpheus.transport.spatial import LinearDiscontinuous
from tests.gates.transport._materials import placeholder_materials
from orpheus.numerics.coefficient_layout import DegreeLayout
from orpheus.numerics.spaces.moment_head import head_layout
from orpheus.numerics import frame_faces


pytestmark = [pytest.mark.foundation]

_L, _NG, _NX = 2, 2, 5

#: ⛔ RE-KEYED 2026-09-02 (#429). The whole module rode ONE 1-D fixture and
#: built its values as ``(L+1, 2L+1, …)`` — the RECTANGULAR spherical-harmonic
#: head — on a rule that now binds the FLAT Legendre head ``(L+1,)``. The
#: claims here (the scalar self-derive) are layout-GENERIC, so
#: the repair is to read the layout off the frame and to run every row on
#: BOTH families: a gate keyed on one certifies the other by accident.
_FAMILIES: tuple[str, ...] = ("flat", "rectangular")


def _sn(family: str = "flat") -> Posed:
    """The 1-D mesh, with the quadrature that induces the requested angular head.

    ``"flat"``    — ``gauss_legendre(4)``: its measure lives on
    :math:`S^2/O(2)_x`, so its frame binds the Legendre basis and the head is
    ``(L+1,)``.
    ``"rectangular"`` — ``level_symmetric(4)``: a full-sphere rule, so the head
    is the harmonics' ``(L+1, 2L+1)``.
    """
    mesh = Mesher(StructuredGeometry.from_homogeneous(1.0, 0)).partition(
        CellsByCount.uniform_width(_NX),
    ).mesh
    quadrature = (
        Quadrature.gauss_legendre(4) if family == "flat"
        else Quadrature.level_symmetric(4)
    )
    # Linear-Discontinuous: the widened (spatial_moments = 2) rows need a
    # scheme that MINTS a moment axis — since CS4c step 6 item 6.2c-iii the
    # moment product's tail is the scheme's own mass-weighted axis (the
    # angular side's rule since CS4b S4), so a widened request on a DD
    # carrier is refused rather than given a Euclidean tail. Width-1 rows
    # are scheme-blind (no tail).
    # the moment space is the frame's ``moment_space_on`` of the angular bulk, widened by the scheme's ``trial_space``
    # at width 2 (``_phase_space_cases.moment_space``)
    return Posed(MaterialMesh(mesh, placeholder_materials(ng=_NG)), quadrature, LinearDiscontinuous())


def _head(sn: Posed, L: int = _L) -> DegreeLayout:
    """The layout of the angular head this mesh's frame induces — read off the frame's Parseval-dressed head (item 6.2c-ii) through the one door (S7)."""
    return head_layout(frame_faces.basis_space(frame_faces.angular_frame(sn.quadrature, L)))


def _field(sn: Posed, spatial_moments: int, seed: int) -> HarmonicMomentFlux:
    tail = (spatial_moments,) if spatial_moments > 1 else ()
    values = np.random.default_rng(seed).standard_normal(
        (*_head(sn).shape, _NG, _NX, *tail)
    )
    return HarmonicMomentFlux(values=values, space=moment_space(sn.mesh, sn.quadrature, _L, sn.scheme if spatial_moments > 1 else None), L=_L, spatial_moments=spatial_moments)


# ── scalar_flux (self-derive contract) ──────────────


class TestViews:
    @pytest.mark.parametrize("family", _FAMILIES)
    def test_scalar_flux_width1_self_derives_the_bulk_factor(self, family):
        """Width-1 scalar_flux derives its target from the product's
        cell-group factor — the carrier's shared bulk mint."""
        sn = _sn(family)
        f = _field(sn, 1, seed=37)
        s = f.scalar_flux()
        npt.assert_array_equal(s.values, f.values[_head(sn).isotropic_slot])
        if s.space != f.space.sub(*f.space.labels[1:]):  # type: ignore[union-attr]
            pytest.fail("scalar target must be the product's bulk factor")

    @pytest.mark.parametrize("family", _FAMILIES)
    def test_scalar_flux_widened_self_derives_the_widened_bulk(self, family):
        """The widened self-derive is HONEST since CS4c step 6 item
        6.2c-iii: the moment product's cell-group factor is the carrier's
        widened bulk carrying the scheme's mass-weighted moment axis (the
        S4 refusal's reason — a Euclidean ``SpatialMomentSpace`` tail
        without that axis — is gone), so the target is the product's own
        factor, equal to the widened angular space's scalar marginal a
        caller used to have to pass."""
        sn = _sn(family)
        f = _field(sn, 2, seed=38)
        s = f.scalar_flux()
        npt.assert_array_equal(s.values, f.values[_head(sn).isotropic_slot])
        if s.space != f.space.sub(*f.space.labels[1:]):  # type: ignore[union-attr]
            pytest.fail("scalar target must be the product's widened bulk factor")
        from orpheus.numerics.basis import SPATIAL_MOMENT_AXIS_LABEL
        assert s.space.axes is not None
        tail = [ax for ax in s.space.axes if ax.label == SPATIAL_MOMENT_AXIS_LABEL]
        if len(tail) != 1 or tail[0] != sn.scheme.moment_axis(sn.mesh.axes):
            pytest.fail("the widened scalar target must carry the scheme's own moment axis")
        # and it equals the explicit-target path
        explicit = f.scalar_flux(space=s.space)
        if explicit.space != s.space:
            pytest.fail("explicit and self-derived targets must be one space")
