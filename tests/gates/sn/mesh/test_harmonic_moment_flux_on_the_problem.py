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

Layer DAG Phase 4 (2026-11-01): this file holds the one row of
tests/gates/transport/fields/test_harmonic_moment_flux.py whose subject is the SN problem's moment-field factory
(truncate's space equals the factory mint on the problem); the source file keeps the field's laws, posed on a
material mesh and the frame's moment_space_on.
"""

from __future__ import annotations

import numpy as np
import numpy.testing as npt
import pytest

from orpheus.geometry import StructuredGeometry
from orpheus.mesh import CellsByCount, Mesher
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.problem import SNProblem
from orpheus.transport.fields.harmonic_moment_flux import HarmonicMomentFlux
from orpheus.transport.spatial import LinearDiscontinuous
from tests.gates.transport._materials import placeholder_materials
from orpheus.numerics.coefficient_layout import DegreeLayout
from orpheus.numerics.spaces.moment_head import head_layout
from orpheus.numerics import frame_faces
from orpheus.transport.boundary_conditions import WallCondition
from tests.gates._exterior import conditions_for

#: Exterior conditions of the tags this file's meshes carry (wave 2: a face's role is a deck or an int tag; the condition a tag
#: stands for is declared at the solve, keyword ``boundary_conditions``).
_CONDITIONS = {0: WallCondition.vacuum()}


pytestmark = [pytest.mark.foundation]

_L, _NG, _NX = 2, 2, 5

#: ⛔ RE-KEYED 2026-09-02 (#429). The whole module rode ONE 1-D fixture and
#: built its values as ``(L+1, 2L+1, …)`` — the RECTANGULAR spherical-harmonic
#: head — on a rule that now binds the FLAT Legendre head ``(L+1,)``. The
#: claims here (truncate, the ℓ-partition, the views) are layout-GENERIC, so
#: the repair is to read the layout off the frame and to run every row on
#: BOTH families: a gate keyed on one certifies the other by accident.
_FAMILIES: tuple[str, ...] = ("flat", "rectangular")


def _sn(family: str = "flat") -> SNProblem:
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
    return SNProblem(mesh, quadrature, placeholder_materials(ng=_NG), scheme=LinearDiscontinuous(), boundary_conditions=conditions_for(mesh, _CONDITIONS))


def _head(sn: SNProblem, L: int = _L) -> DegreeLayout:
    """The layout of the angular head this mesh's frame induces — read off the frame's Parseval-dressed head (item 6.2c-ii) through the one door (S7)."""
    return head_layout(frame_faces.basis_space(frame_faces.angular_frame(sn.quad, L)))


def _field(sn: SNProblem, spatial_moments: int, seed: int) -> HarmonicMomentFlux:
    tail = (spatial_moments,) if spatial_moments > 1 else ()
    values = np.random.default_rng(seed).standard_normal(
        (*_head(sn).shape, _NG, _NX, *tail)
    )
    return HarmonicMomentFlux.from_problem_and_L(
        values, sn, _L, spatial_moments=spatial_moments,
    )


# ── truncate — the space-derived rebuild (both widths) ───────────────


class TestTruncate:
    @pytest.mark.parametrize("family", _FAMILIES)
    @pytest.mark.parametrize("sm", [1, 2])
    def test_truncate_matches_the_factory_mint_at_both_widths(self, sm, family):
        """The truncated field's space content-equals the carrier's OWN
        mint at (mesh, L_new, spatial_moments) — ``SNProblem.moment_space``
        since CS4c step 6 item 6.2b, the cached object the factory reads
        — the single-source done-when of the space-derived rebuild. The sm=2 row is #399's
        FLIPPED witness (pre-S6.1 it raised the widened defer; pre-S4 it
        was the census's broadcast ValueError). Reddened by dropping the
        tail factor from the rebuild (the pre-S6.1 factors[1]-only
        spelling) or by re-introducing the un-tailed new_shape."""
        sn = _sn(family)
        f = _field(sn, sm, seed=20 + sm)
        L_new = 1
        g = f.truncate(L_new)
        tail = (sm,) if sm > 1 else ()
        expected_shape = (*_head(sn, L_new).shape, _NG, _NX, *tail)
        if g.values.shape != expected_shape:
            pytest.fail(
                f"truncate(sm={sm}) shape {g.values.shape} != "
                f"{expected_shape}"
            )
        factory_space = HarmonicMomentFlux.from_problem_and_L(
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
        head_new = _head(sn, L_new)
        for l in range(L_new + 1):
            npt.assert_array_equal(
                g.values[head_new.degree_block(l)],
                f.values[_head(sn).degree_block(l)],
            )


# ── the ℓ-partition (both widths) ────────────────────────────────────


# ── l_block (view) + scalar_flux (self-derive contract) ──────────────


