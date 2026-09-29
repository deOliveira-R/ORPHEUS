"""The one-material reference generators read a geometry through one refusal.

``Spectrum``, ``MomentSpace``, ``BasisSpace`` and ``Billiard`` solve on one
material filling a slab, or a solid cylinder or sphere. Before P1 step 2
each read that body off the geometry on its own: a geometry of several
materials was read as its first material, and (once step 2 admitted hollow bodies) a hollow
sphere would have been read as a solid one whose radius is the shell
thickness. :func:`~orpheus.derivations.common.homogeneous_body.homogeneous_body`
is now the one reading, and it refuses both (the user's ruling of
2026-09-29, ``.claude/plans/reference_cache.md``, "P1 step 2 opened").
"""
from __future__ import annotations

import pytest

from orpheus.derivations.common.homogeneous_body import (
    HomogeneousBody,
    homogeneous_body,
)
from orpheus.derivations.common.xs_library import get_mixture
from orpheus.derivations.continuous.fn_method.moment_space import MomentSpace
from orpheus.derivations.continuous.galerkin_spectral.basis_space import BasisSpace
from orpheus.derivations.continuous.singular_eigenfunction.spectrum import Spectrum
from orpheus.derivations.continuous.trajectory_resolvent.billiard import Billiard
from orpheus.geometry import BC, CoordSystem, StructuredGeometry

pytestmark = pytest.mark.foundation


def _hollow_sphere() -> StructuredGeometry:
    return StructuredGeometry(
        coord=CoordSystem.SPHERICAL, breakpoints=(0.5, 2.0), mat_ids=(0,),
        boundaries=(BC.reflective, BC.vacuum),
    )


def _two_material_slab() -> StructuredGeometry:
    return StructuredGeometry(
        coord=CoordSystem.CARTESIAN, breakpoints=(0.0, 1.0, 2.0), mat_ids=(0, 1),
        boundaries=(BC.vacuum, BC.vacuum),
    )


class TestHomogeneousBody:
    @pytest.mark.parametrize(
        "coord, breakpoints, boundaries",
        [
            (CoordSystem.CARTESIAN, (-1.0, 2.0), (BC.vacuum, BC.vacuum)),
            (CoordSystem.CYLINDRICAL, (0.0, 2.0), (BC.vacuum,)),
            (CoordSystem.SPHERICAL, (0.0, 2.0), (BC.vacuum,)),
        ],
        ids=["slab-any-origin", "solid-cylinder", "solid-sphere"],
    )
    def test_a_body_is_read(self, coord, breakpoints, boundaries):
        g = StructuredGeometry(
            coord=coord, breakpoints=breakpoints, mat_ids=(7,), boundaries=boundaries,
        )
        assert homogeneous_body(g, owner="test") == HomogeneousBody(
            coord=coord, extent_cm=breakpoints[1] - breakpoints[0], mat_id=7,
        )

    def test_several_materials_are_refused(self):
        with pytest.raises(ValueError, match="holds the materials \\[0, 1\\]"):
            homogeneous_body(_two_material_slab(), owner="test")

    def test_several_intervals_of_one_material_are_one_body(self):
        """The two-sided leg: interior breakpoints between intervals of one
        material are not interfaces, so the body is read, over its full width."""
        g = StructuredGeometry(
            coord=CoordSystem.SPHERICAL, breakpoints=(0.0, 1.0, 2.5), mat_ids=(3, 3),
            boundaries=(BC.vacuum,),
        )
        assert homogeneous_body(g, owner="test") == HomogeneousBody(
            coord=CoordSystem.SPHERICAL, extent_cm=2.5, mat_id=3,
        )

    def test_a_hollow_body_is_refused(self):
        with pytest.raises(ValueError, match="solid spherical body centred at r = 0"):
            homogeneous_body(_hollow_sphere(), owner="test")


#: Each generator, built on a geometry; its refusal names it.
_GENERATORS = {
    "MomentSpace": lambda g: MomentSpace(geometry=g, materials={0: get_mixture("A", "1g")}),
    "BasisSpace": lambda g: BasisSpace(geometry=g, materials={0: get_mixture("A", "1g")}),
    "Spectrum": lambda g: Spectrum(geometry=g, materials={0: get_mixture("A", "1g")}),
    "Billiard": lambda g: Billiard(geometry=g, materials={0: get_mixture("A", "1g")}, alpha=0.0),
}


@pytest.mark.parametrize("owner", _GENERATORS)
@pytest.mark.parametrize(
    "geometry, fragment",
    [(_hollow_sphere, "is hollow"), (_two_material_slab, "holds the materials")],
    ids=["hollow-sphere", "two-material-slab"],
)
def test_every_generator_refuses_what_it_cannot_read(owner, geometry, fragment):
    """A ROUTE gate: every generator reaches the one refusal, named."""
    with pytest.raises(ValueError, match=fragment) as caught:
        _GENERATORS[owner](geometry())
    assert str(caught.value).startswith(owner)


def test_billiard_reads_the_body_material():
    """Billiard takes its cross sections from the body's material id, not
    from key 0 (qa, 2026-09-29: with mat_ids (3,) and materials {0: A, 3: C}
    it returned A's k_inf)."""
    def sphere(mat_id: int) -> StructuredGeometry:
        return StructuredGeometry(
            coord=CoordSystem.SPHERICAL, breakpoints=(0.0, 2.0), mat_ids=(mat_id,),
            boundaries=(BC.reflective,),
        )

    decoy, body = get_mixture("A", "1g"), get_mixture("B", "1g")
    read = Billiard(geometry=sphere(3), materials={0: decoy, 3: body}, alpha=0.0)
    reference = Billiard(geometry=sphere(0), materials={0: body}, alpha=0.0)
    assert read.xs_payload == reference.xs_payload
    with pytest.raises(ValueError, match="the body's material id 3"):
        Billiard(geometry=sphere(3), materials={0: body}, alpha=0.0)
