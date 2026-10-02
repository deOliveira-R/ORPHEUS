r"""Fixtures of the step-8 gates (#405 P1 step 8, spec §1.8): mixtures that carry
exactly the channels a row needs, geometries and functions.

Not a test module. Each mixture is built with ``Mixture(...)`` DIRECTLY (never
``make_mixture``, which nulls ``Sig2`` and ``SigL``: test-architect lessons §3),
and each states the emission channels it carries, so a row that reads one
channel's carriage has a fixture on which that channel alone is non-zero.
"""

from __future__ import annotations

import numpy as np
from scipy.sparse import csr_matrix

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.geometry import BC, CoordSystem, StructuredGeometry


def _mixture(*, fission: bool, scattering: bool, n2n: bool, ng: int = 2) -> Mixture:
    """A mixture whose emission cells are non-zero exactly where asked."""
    if ng == 2:
        s0 = np.array([[0.30, 0.05], [0.02, 0.80]])  # [g_from, g_to], with upscatter
        n2 = np.array([[0.004, 0.003], [0.0, 0.0]])
        sig_c = np.array([0.01, 0.06])
        sig_f = np.array([0.008, 0.10])
        chi = np.array([0.97, 0.03])
    else:
        s0 = np.full((ng, ng), 0.1 / ng)
        n2 = np.full((ng, ng), 0.001)
        sig_c = np.full(ng, 0.05)
        sig_f = np.full(ng, 0.04)
        chi = np.full(ng, 1.0 / ng)
    s0 = s0 if scattering else np.zeros((ng, ng))
    n2 = n2 if n2n else np.zeros((ng, ng))
    sig_f = sig_f if fission else np.zeros(ng)
    sig_l = np.zeros(ng)
    sig_t = sig_c + sig_l + sig_f + s0.sum(1) + n2.sum(1)
    return Mixture(
        SigC=sig_c, SigL=sig_l, SigF=sig_f, SigP=2.4 * sig_f, SigT=sig_t,
        SigS=(csr_matrix(s0), csr_matrix(0.2 * s0)), Sig2=(csr_matrix(n2),),
        chi=chi if fission else np.zeros(ng),
    )


def fuel(ng: int = 2) -> Mixture:
    """Carries all three emission channels; 2 groups with upscatter and (n,2n)."""
    return _mixture(fission=True, scattering=True, n2n=True, ng=ng)


def moderator() -> Mixture:
    """Carries the scattering emission only."""
    return _mixture(fission=False, scattering=True, n2n=False)


def n2n_only() -> Mixture:
    """Carries the (n,2n) emission only."""
    return _mixture(fission=False, scattering=False, n2n=True)


def fission_only() -> Mixture:
    """Carries the fission emission only."""
    return _mixture(fission=True, scattering=False, n2n=False)


def slab2(mat_ids: tuple[int, int] = (0, 1)) -> StructuredGeometry:
    """A 2-interval slab, reflective | vacuum."""
    return StructuredGeometry.slab((0.0, 1.0, 3.0), mat_ids, left=BC.reflective, right=BC.vacuum)


def body(coord: CoordSystem, mat_ids: tuple[int, ...] = (0, 1)) -> StructuredGeometry:
    """A 2-interval solid body in ``coord`` (a slab gets reflective | vacuum)."""
    breakpoints = (0.0, 1.0, 3.0)
    if coord is CoordSystem.CARTESIAN:
        return StructuredGeometry.slab(breakpoints, mat_ids, left=BC.reflective, right=BC.vacuum)
    if coord is CoordSystem.CYLINDRICAL:
        return StructuredGeometry.cylinder(breakpoints, mat_ids, outer=BC.vacuum)
    return StructuredGeometry.sphere(breakpoints, mat_ids, outer=BC.vacuum)
