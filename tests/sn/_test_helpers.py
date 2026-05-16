"""Shared test helpers for the SN test suite.

Issue #197 PR-TYPED-0 introduced ``materials`` as a REQUIRED parameter
on :class:`SNMesh`.  Many geometry-only tests don't consume cross-
section values — they exercise sweep DAGs, BC realisation,
quadrature, cache structure, etc.  ``placeholder_materials`` provides
a minimal-but-valid :class:`Mixture` dict that those tests can hand
to :class:`SNMesh` so the construction succeeds without inviting any
real cross-section semantics into the test.

Tests that DO need realistic cross sections continue to use
``orpheus.derivations.common.xs_library.get_mixture`` etc. — this
helper is for the geometry-only call sites.
"""
from __future__ import annotations

import numpy as np
from scipy.sparse import csr_matrix


def placeholder_materials(
    ng: int = 1, mat_ids: tuple[int, ...] = (0,),
) -> dict:
    """Build a placeholder ``{mat_id: Mixture}`` dict for SNMesh tests.

    Parameters
    ----------
    ng
        Number of energy groups.  All Mixtures will report this value
        via :attr:`Mixture.ng`.
    mat_ids
        Material ids to include in the dict.  Default ``(0,)`` covers
        the common case where the mesh's ``mat_ids`` / ``mat_map`` is
        zeros.

    Returns
    -------
    dict[int, Mixture]
        Each entry has ``SigT = ones(ng)`` and all other cross sections
        zero.  Suitable for SNMesh tests that don't compute physical
        quantities from the materials.
    """
    from orpheus.data.macro_xs.mixture import Mixture
    z = np.zeros(ng)
    z_mat = csr_matrix(np.zeros((ng, ng)))
    return {
        int(mid): Mixture(
            SigC=z.copy(),
            SigL=z.copy(),
            SigF=z.copy(),
            SigP=z.copy(),
            SigT=np.ones(ng),
            SigS=[z_mat],
            Sig2=z_mat,
            chi=z.copy(),
        )
        for mid in mat_ids
    }
