r"""The A|B|A cross-check problem, defined once (#405 P2 step 7b.2, spec §1.7b.2).

Fuel A | moderator B | fuel A at outer radii 0.5, 1.5, 2.0 cm, 2 groups,
reflective at r = R, on a sphere and on a cylinder: the problem of the
trajectory-resolvent cross-check rows. It is posed as a
:class:`~orpheus.specification.specification.GeometrySpecification` from
ISOTROPIC (P0) mixtures: the trajectory resolvent solves isotropic scattering
only, and the SN rows run at ``scattering_order=0``, while the xs_library's
``get_mixture("B", "2g")`` carries a P1 moment (mean cosine 0.6). A
specification built from those would pose a problem neither side solves.

This module is the home step 7b.2.3 migrates the retiring
``_certified_agreement.py``'s ABA definitions into; until then it holds only
the specification and its materials.
"""
from __future__ import annotations

import functools

import numpy as np

from orpheus.data.cells import CellCoefficient, Channel
from orpheus.data.macro_xs.mixture import Mixture
from orpheus.data.materials import Materials
from orpheus.derivations.common.xs_library import get_xs, make_mixture
from orpheus.geometry import BC, CoordSystem, StructuredGeometry
from orpheus.numerics.question import Eigen
from orpheus.specification.specification import GeometrySpecification

#: The outer radius of each region, cm.
ABA_RADII = (0.5, 1.5, 2.0)
#: The material of each region: fuel A (id 0), moderator B (id 1), fuel A.
ABA_MATERIAL_IDS = (0, 1, 0)


def isotropic_mixture(key: str) -> Mixture:
    """The xs_library's 2-group mixture ``key`` with its P0 scattering only (no ``sig_s1``)."""
    xs = get_xs(key, "2g")
    return make_mixture(sig_t=xs["sig_t"], sig_c=xs["sig_c"], sig_f=xs["sig_f"], nu=xs["nu"], chi=xs["chi"], sig_s=xs["sig_s"])


@functools.cache
def aba_materials() -> Materials:
    """Fuel A as material 0, moderator B as material 1, both isotropic."""
    return Materials({0: isotropic_mixture("A"), 1: isotropic_mixture("B")})


def aba_geometry(coord: CoordSystem) -> StructuredGeometry:
    """The A|B|A body on ``coord`` (spherical or cylindrical), reflective at r = R."""
    return StructuredGeometry.from_thicknesses(
        coord=coord, thicknesses=tuple(np.diff((0.0, *ABA_RADII))), mat_ids=ABA_MATERIAL_IDS,
        boundaries=(BC.reflective,),
    )


@functools.cache
def aba_specification(coord: CoordSystem) -> GeometrySpecification:
    """The k question on the A|B|A body: the specification the trajectory-resolvent reference answers."""
    return GeometrySpecification(aba_materials(), aba_geometry(coord), Eigen(CellCoefficient.every(Channel.FISSION_EMISSION)))
