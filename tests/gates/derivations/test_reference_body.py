"""The reference generators read a geometry as one of four body shapes, and each
serves exactly the shapes its solvers solve (P1 step 2b).

:func:`~orpheus.derivations.common.reference_body.reference_body` reads a
:class:`~orpheus.geometry.StructuredGeometry` as a homogeneous body, a hollow
body of one material, a symmetric reflected slab, or a layered body; it is
total and knows no solver. Each generator matches on the shape and refuses
the ones it does not serve through the one door,
:func:`~orpheus.derivations.common.reference_body.refuse_unserved`. Before
step 2b one shared refusal blocked every multi-material body, including the
ones a solver under the generator already solves: the trajectory resolvent's
(``Billiard``, deleted in step (e2) of the characteristic-reference campaign)
multi-region sphere and cylinder, its hollow sphere and annulus (#190, #421), and
MomentSpace's reflected slab (Neshat and Maiorino 1980, Sood's problem 4,
which only a test adapter reached). The user's rulings of 2026-09-29 are in
``.claude/plans/reference_cache.md``, "P1 step 2b".
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.reference_body import (
    HollowBody,
    HomogeneousBody,
    LayeredBody,
    ReflectedSlab,
    reference_body,
)
from orpheus.derivations.common.xs_library import get_mixture, make_mixture
from orpheus.derivations.continuous.fn_method.moment_space import MomentSpace
from orpheus.derivations.continuous.galerkin_spectral.basis_space import BasisSpace
from orpheus.derivations.continuous.singular_eigenfunction.spectrum import Spectrum
from orpheus.geometry import CoordSystem, StructuredGeometry
from tests.gates.derivations import _posed_walls as walls
from tests.gates.derivations._posed_walls import MIRROR, VACUUM, WHITE, partial

pytestmark = pytest.mark.foundation

_SLAB, _CYL, _SPH = CoordSystem.CARTESIAN, CoordSystem.CYLINDRICAL, CoordSystem.SPHERICAL


def _body(coord, breakpoints, mat_ids, declared=None) -> walls.Posed:
    """The geometry and its exterior; ``None`` states vacuum on every boundary face."""
    return walls.posed(coord, breakpoints, mat_ids, declared)


def _geometry(coord, breakpoints, mat_ids, declared=None) -> StructuredGeometry:
    return _body(coord, breakpoints, mat_ids, declared).geometry


def _thick(thicknesses, mat_ids, declared=None) -> walls.Posed:
    return walls.from_thicknesses(_SLAB, thicknesses, mat_ids, declared)


def _unit_sigma_t_mixture(c: float, sigma_t: float = 1.0):
    """One group, total cross section ``sigma_t`` (1 by default: the unit
    convention of Sood's two-media slabs), secondaries ratio ``c``: a
    scatterer, multiplying by fission when c > 0.9. Every cross section
    scales with ``sigma_t``, so ``c`` does not depend on it."""
    sig_s = min(c, 0.9)
    nu_sig_f = c - sig_s
    nu = 2.5 if nu_sig_f else 0.0
    sig_f = nu_sig_f / nu if nu_sig_f else 0.0
    return make_mixture(
        sig_t=np.array([1.0]) * sigma_t,
        sig_c=np.array([1.0 - sig_s - sig_f]) * sigma_t,
        sig_f=np.array([sig_f]) * sigma_t,
        nu=np.array([nu]),
        chi=np.array([1.0 if nu_sig_f else 0.0]),
        sig_s=np.array([[sig_s]]) * sigma_t,
    )


# ─────────────────────────────────────────────────────────────────────
# The classification: total, one shape per geometry, knows no solver
# ─────────────────────────────────────────────────────────────────────


class TestTheClassification:
    @pytest.mark.parametrize("coord", [_SLAB, _CYL, _SPH], ids=lambda c: c.name.lower())
    def test_one_material_solid_is_homogeneous(self, coord):
        body = reference_body(_geometry(coord, (0.0, 2.0), (4,)))
        assert body == HomogeneousBody(coord=coord, extent_cm=2.0, mat_id=4)

    def test_intervals_of_one_material_merge(self):
        """An interior breakpoint between intervals of one material is not an
        interface: the body is read over its full width."""
        body = reference_body(_geometry(_SPH, (0.0, 1.0, 2.5), (3, 3)))
        assert body == HomogeneousBody(coord=_SPH, extent_cm=2.5, mat_id=3)

    @pytest.mark.parametrize("coord", [_CYL, _SPH], ids=lambda c: c.name.lower())
    def test_one_material_hollow_is_hollow(self, coord):
        g = _geometry(coord, (0.5, 2.0), (1,), declared=(MIRROR, VACUUM))
        assert reference_body(g) == HollowBody(
            coord=coord, inner_radius_cm=0.5, outer_radius_cm=2.0, mat_id=1,
        )

    def test_a_symmetric_three_run_slab_is_a_reflected_slab(self):
        """Built from thicknesses, as a registry states Sood's slabs: the
        right reflector's width carries the breakpoints' rounding, and the
        classification still reads it as symmetric."""
        g = _thick((0.5, 0.8603, 0.5), (1, 0, 1)).geometry
        body = reference_body(g)
        assert isinstance(body, ReflectedSlab)
        assert (body.core_mat_id, body.reflector_mat_id) == (0, 1)
        assert body.reflector_width_cm == 0.5
        assert body.coord is _SLAB

    @pytest.mark.parametrize(
        "thicknesses, mat_ids",
        [
            ((0.5, 1.0, 0.6), (1, 0, 1)),
            ((0.5, 1.0, 0.5), (1, 0, 2)),
            ((1.0, 0.5), (0, 1)),
            ((0.1, 0.5, 0.1, 0.2), (2, 0, 2, 3)),
        ],
        ids=["asymmetric-widths", "different-reflectors", "one-sided", "four-run"],
    )
    def test_any_other_layered_slab_is_layered(self, thicknesses, mat_ids):
        """The discriminating rows: Sood's problem 3 (one-sided) and 30 (the
        Fe/U/Fe/Na slab) are layered slabs, not reflected slabs."""
        g = _thick(thicknesses, mat_ids).geometry
        assert isinstance(reference_body(g), LayeredBody)

    def test_a_layered_solid_sphere(self):
        body = reference_body(_geometry(_SPH, (0.0, 1.0, 3.0), (0, 1)))
        assert isinstance(body, LayeredBody)
        assert body == LayeredBody(
            coord=_SPH, breakpoints=(0.0, 1.0, 3.0), mat_ids=(0, 1), is_hollow=False,
        )
        assert not body.is_hollow

    def test_a_layered_hollow_cylinder(self):
        g = _geometry(_CYL, (0.5, 1.0, 3.0), (0, 1), declared=(MIRROR, VACUUM))
        body = reference_body(g)
        assert isinstance(body, LayeredBody) and body.is_hollow


# ─────────────────────────────────────────────────────────────────────
# Each owner serves exactly its shapes, and refuses the rest through the door
# ─────────────────────────────────────────────────────────────────────

_ONE_GROUP = {0: get_mixture("A", "1g"), 1: get_mixture("B", "1g")}


_SHAPES = {
    "hollow-sphere": lambda: _body(_SPH, (0.5, 2.0), (0,), (MIRROR, VACUUM)),
    "layered-sphere": lambda: _body(_SPH, (0.0, 1.0, 2.0), (0, 1)),
    "layered-slab": lambda: _thick((1.0, 0.5), (0, 1)),
    "reflected-slab": lambda: _thick((0.5, 1.0, 0.5), (1, 0, 1)),
    "hollow-layered-sphere": lambda: _body(_SPH, (0.5, 1.0, 2.0), (0, 1), (MIRROR, VACUUM)),
}

#: (owner, the shapes it refuses): the shapes it serves are tested by route below.
_REFUSALS = [
    ("Spectrum", lambda g: Spectrum(**g._asdict(), materials=_ONE_GROUP),
     ["hollow-sphere", "layered-sphere", "layered-slab", "reflected-slab", "hollow-layered-sphere"]),
    ("BasisSpace", lambda g: BasisSpace(**g._asdict(), materials=_ONE_GROUP),
     ["hollow-sphere", "layered-sphere", "layered-slab", "reflected-slab", "hollow-layered-sphere"]),
    ("MomentSpace", lambda g: MomentSpace(**g._asdict(), materials=_ONE_GROUP),
     ["hollow-sphere", "layered-sphere", "layered-slab", "hollow-layered-sphere"]),
]


@pytest.mark.parametrize(
    "owner, build, shape",
    [pytest.param(o, b, s, id=f"{o}-{s}") for o, b, shapes in _REFUSALS for s in shapes],
)
def test_an_unserved_shape_is_refused_through_the_door(owner, build, shape):
    with pytest.raises(NotImplementedError, match="#536") as caught:
        build(_SHAPES[shape]())
    assert str(caught.value).startswith(f"{owner} does not solve a ")


class TestMomentSpaceReflectedSlabPreconditions:
    def test_unequal_total_cross_sections_are_refused(self):
        g = _SHAPES["reflected-slab"]()
        materials = {0: _unit_sigma_t_mixture(1.5), 1: get_mixture("B", "1g")}
        with pytest.raises(NotImplementedError, match="share one total cross section"):
            MomentSpace(**g._asdict(), materials=materials)

    def test_several_groups_are_refused(self):
        g = _SHAPES["reflected-slab"]()
        materials = {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")}
        with pytest.raises(NotImplementedError, match="one-group reflected slab"):
            MomentSpace(**g._asdict(), materials=materials)

    def test_its_flux_reconstruction_is_refused(self):
        g = _SHAPES["reflected-slab"]()
        materials = {0: _unit_sigma_t_mixture(1.5), 1: _unit_sigma_t_mixture(0.9)}
        with pytest.raises(NotImplementedError, match="flux reconstruction"):
            MomentSpace(**g._asdict(), materials=materials).reconstruct_flux()


# ─────────────────────────────────────────────────────────────────────
# Each route is bit-identical to the solver it reaches
# ─────────────────────────────────────────────────────────────────────


def test_moment_space_routes_a_reflected_slab():
    """Sood's problem 4 configuration (c_core 1.50, c_reflector 0.90,
    reflector 0.5 mfp): the facade and the bare F_N function agree bit for
    bit. The published value is gated in
    ``tests/gates/cross_method/test_eigenvalue.py::test_fn_reflected_slab_matches_truth``,
    whose adapter now goes through this facade."""
    from orpheus.derivations.continuous.fn_method.slab.reflected import (
        solve_fn_slab_reflected_critical,
    )

    g = _thick((0.5, 0.8603, 0.5), (1, 0, 1))
    materials = {0: _unit_sigma_t_mixture(1.5), 1: _unit_sigma_t_mixture(0.9)}
    solution = MomentSpace(**g._asdict(), materials=materials, fn_order=7).solve_critical()
    bare = solve_fn_slab_reflected_critical(
        c_core=float(materials[0].scattering_ratio[0]),
        c_reflector=float(materials[1].scattering_ratio[0]),
        reflector_half_thickness=0.5,
        n_modes=7,
    )
    assert solution.parameter_value == bare.tau_critical_mfp
    assert solution.parameter_kind == "core_half_thickness_mfp"


# ─────────────────────────────────────────────────────────────────────
# The laws are part of what an owner serves (the user's ruling, 2026-09-29)
# ─────────────────────────────────────────────────────────────────────

from orpheus.derivations.common.reference_body import specular_albedos  # noqa: E402


class TestSpecularAlbedo:
    @pytest.mark.parametrize(
        "wall, albedo",
        [
            (VACUUM, 0.0),
            (MIRROR, 1.0),
            (partial(0.3), 0.3),
            (partial(0.7), 0.7),
        ],
        ids=["vacuum", "mirror", "partial-0.3", "partial-0.7"],
    )
    @pytest.mark.parametrize("coord", [_SLAB, _SPH], ids=lambda c: c.name.lower())
    def test_a_specular_wall_reads_as_its_albedo(self, coord, wall, albedo):
        """A mirror is a deck on a slab and the tag bound to ``(1, 0)`` on a sphere; both read 1."""
        body = _body(coord, (0.0, 2.0), (0,), (wall,) if coord is _SPH else (wall, wall))
        got = specular_albedos(body.geometry, body.exterior, owner="test")
        assert got == (albedo,) * len(body.geometry.face_roles)

    @pytest.mark.parametrize("wall", [WHITE, walls.white(0.7), walls.WRAP], ids=["white", "partial-white", "wrap"])
    def test_any_other_wall_is_refused(self, wall):
        body = _body(_SLAB, (0.0, 2.0), (0,), (wall, wall))
        with pytest.raises(NotImplementedError, match="#536"):
            specular_albedos(body.geometry, body.exterior, owner="test")


def _reflected_slab(faces=(VACUUM, VACUUM), thicknesses=(0.5, 1.0, 0.5)) -> walls.Posed:
    return _thick(thicknesses, (1, 0, 1), faces)


class TestTheLawsAreServed:
    def test_moment_space_refuses_a_reflected_slab_with_mirror_faces(self):
        """The F1 witness: before the laws joined the served pattern, a
        reflected slab with reflective faces returned the vacuum answer, bit
        for bit (the elegance review, 2026-09-29)."""
        materials = {0: _unit_sigma_t_mixture(1.5), 1: _unit_sigma_t_mixture(0.9)}
        with pytest.raises(NotImplementedError, match="reflecting boundaries"):
            MomentSpace(**_reflected_slab((MIRROR, MIRROR))._asdict(), materials=materials)

    @pytest.mark.parametrize(
        "owner, build",
        [("MomentSpace", lambda g: MomentSpace(**g._asdict(), materials=_ONE_GROUP)),
         ("BasisSpace", lambda g: BasisSpace(**g._asdict(), materials=_ONE_GROUP))],
    )
    def test_the_bare_solvers_refuse_a_reflecting_boundary(self, owner, build):
        g = _body(_SPH, (0.0, 2.0), (0,), (partial(0.5),))
        with pytest.raises(NotImplementedError, match="reflecting boundaries") as caught:
            build(g)
        assert str(caught.value).startswith(owner)

    def test_spectrum_refuses_a_slab_whose_faces_differ(self):
        """Atalay's slab puts one R on both faces."""
        g = _body(_SLAB, (0.0, 2.0), (0,), (VACUUM, partial(0.5)))
        with pytest.raises(NotImplementedError, match="unequal face albedos"):
            Spectrum(**g._asdict(), materials=_ONE_GROUP)

    def test_spectrum_refuses_a_reflected_cylinder(self):
        g = _body(_CYL, (0.0, 2.0), (0,), (partial(0.5),))
        with pytest.raises(NotImplementedError, match="reflected cylinder"):
            Spectrum(**g._asdict(), materials=_ONE_GROUP)

    @pytest.mark.parametrize(
        "build",
        [lambda g: Spectrum(**g._asdict(), materials=_ONE_GROUP)],
        ids=["Spectrum"],
    )
    def test_a_white_wall_is_refused(self, build):
        with pytest.raises(NotImplementedError, match="another angular shape"):
            build(_body(_SPH, (0.0, 2.0), (0,), (WHITE,)))


class TestReflectedSlabRoute:
    def test_the_symmetry_tolerance_meets_a_real_rounding(self):
        """F9: from the thicknesses (0.3, 1.1, 0.3) the right reflector's
        width differs from the left's by a fraction of an ulp (the fold's
        rounding); the slab is still symmetric. Control: a 5-ulp asymmetry
        is a layered slab (the tolerance is 4 ulp)."""
        import math

        g = _reflected_slab(thicknesses=(0.3, 1.1, 0.3)).geometry
        bp = g.breakpoints
        assert (bp[1] - bp[0]) != (bp[3] - bp[2])
        assert isinstance(reference_body(g), ReflectedSlab)
        r_3 = bp[3] + 5 * math.ulp(bp[3])
        skewed = _geometry(_SLAB, (bp[0], bp[1], bp[2], r_3), (1, 0, 1))
        assert isinstance(reference_body(skewed), LayeredBody)

    def test_the_reflector_width_is_converted_to_mean_free_paths(self):
        """qa's 2026-09-29 row: every reflected fixture had Sigma_t = 1, so
        the cm -> mfp conversion was never exercised. With Sigma_t = 2 and a
        0.25 cm reflector the body is the Sigma_t = 1, 0.5 cm body in mean
        free paths: the same tau, bit for bit."""
        unit = MomentSpace(
            **_reflected_slab(thicknesses=(0.5, 1.0, 0.5))._asdict(),
            materials={0: _unit_sigma_t_mixture(1.5), 1: _unit_sigma_t_mixture(0.9)},
        ).solve_critical()
        doubled = MomentSpace(
            **_reflected_slab(thicknesses=(0.25, 0.5, 0.25))._asdict(),
            materials={0: _unit_sigma_t_mixture(1.5, 2.0), 1: _unit_sigma_t_mixture(0.9, 2.0)},
        ).solve_critical()
        assert doubled.metadata["reflector_half_thickness_mfp"] == 0.5
        assert doubled.parameter_value == unit.parameter_value


