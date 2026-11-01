r"""The infinite-medium problem's quotient point and a genuine one-cell mesh are different spaces.

Moved from ``tests/gates/numerics/test_space_of_axes.py`` (row B9) by layer DAG Phase 4 (2026-11-01): it compares the
homogeneous problem's pose (``orpheus.homogeneous.solver._pose_space``) with a transport ``MaterialMesh``'s bulk, so
its home is the homogeneous gates (which may import transport); the numerics file keeps ``of_axes``'s own laws.
"""
from __future__ import annotations

import numpy as np
import pytest

pytestmark = pytest.mark.foundation


def _require(condition: bool, message: str) -> None:
    if not condition:
        pytest.fail(message)


def test_quotient_point_and_a_genuine_one_cell_mesh_are_DIFFERENT_spaces() -> None:
    r"""B9 ⭐ — the collapse doctrine's retrodiction (A.5 row 4), mechanized.

    The infinite-medium problem's pose (``_pose_space``: Energy ⊗ the
    counting point — the normalized "per unit volume" density convention;
    until the CS4c coda a fabricated one-cell carrier's ``bulk_space``
    minted the same space) is the QUOTIENT side. A genuine one-cell slab of
    width 2 keeps ``V = 2`` BY THE DATA. Both spaces have shape
    ``(ng, 1)`` — shape carries nothing; the derived NAME (hence space
    identity) is the only discriminator.

    ⚠ Per F2 (measured): this distinction is provably INVISIBLE to
    ``.H`` — a scalar metric commutes with every operator, so ``V = 2``
    reads bit-identical to ``V = 1`` on every adjoint — which is exactly
    why this identity gate exists and why no adjoint-flavoured gate can
    replace it (mutation M17's MUST-STAY-GREEN column is the proof).

    # CS1.5 re-point: ``bulk_space`` moves to Medium; this gate re-points
    # with it.
    """
    from orpheus.data.macro_xs.synthetic import get_mixture
    from orpheus.geometry import StructuredGeometry
    from orpheus.geometry.deck import MirrorDeck
    from orpheus.homogeneous.solver import _pose_space
    from orpheus.mesh import CellsByCount, Mesher
    from orpheus.transport.mesh.material_mesh import MaterialMesh

    mix = get_mixture("A", "2g")
    quotient = _pose_space(mix)
    one_cell_mesh = MaterialMesh(
        Mesher(StructuredGeometry.from_homogeneous(2.0, MirrorDeck())).partition(
            CellsByCount.uniform_width(1),
        ).mesh,
        {0: mix},
    )
    _require(
        bool(np.array_equal(one_cell_mesh.volumes, [2.0])),
        "precondition lost: the one-cell slab no longer has V = 2",
    )
    one_cell = one_cell_mesh.bulk_space
    _require(quotient.shape == one_cell.shape == (2, 1), "precondition: same shape")
    _require(
        quotient != one_cell,
        "the quotient point and a genuine one-cell mesh collapsed to ONE "
        "space — the measure fell out of the identity",
    )
