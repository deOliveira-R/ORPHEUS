"""The reference generators read a geometry as one of four body shapes, and each
serves exactly the shapes its solvers solve (P1 step 2b).

:func:`~orpheus.derivations.common.reference_body.reference_body` reads a
:class:`~orpheus.geometry.StructuredGeometry` as a homogeneous body, a hollow
body of one material, a symmetric reflected slab, or a layered body; it is
total and knows no solver. Each generator matches on the shape and refuses
the ones it does not serve through the one door,
:func:`~orpheus.derivations.common.reference_body.refuse_unserved`. Before
step 2b one shared refusal blocked every multi-material body, including the
ones a solver under the generator already solves: Billiard's multi-region
sphere and cylinder, its hollow sphere and annulus (#190, #421), and
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
from orpheus.derivations.common.xs_library import get_mixture, get_xs, make_mixture
from orpheus.derivations.continuous.fn_method.moment_space import MomentSpace
from orpheus.derivations.continuous.galerkin_spectral.basis_space import BasisSpace
from orpheus.derivations.continuous.singular_eigenfunction.spectrum import Spectrum
from orpheus.derivations.continuous.trajectory_resolvent.billiard import Billiard
from orpheus.geometry import BC, CoordSystem, StructuredGeometry

pytestmark = pytest.mark.foundation

_SLAB, _CYL, _SPH = CoordSystem.CARTESIAN, CoordSystem.CYLINDRICAL, CoordSystem.SPHERICAL


def _geometry(coord, breakpoints, mat_ids, boundaries=None) -> StructuredGeometry:
    if boundaries is None:
        solid = coord is _SLAB or breakpoints[0] == 0.0
        boundaries = (BC.vacuum,) if (coord is not _SLAB and solid) else (BC.vacuum, BC.vacuum)
    return StructuredGeometry(
        coord=coord, breakpoints=breakpoints, mat_ids=mat_ids, boundaries=boundaries,
    )


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
        g = _geometry(coord, (0.5, 2.0), (1,), boundaries=(BC.reflective, BC.vacuum))
        assert reference_body(g) == HollowBody(
            coord=coord, inner_radius_cm=0.5, outer_radius_cm=2.0, mat_id=1,
        )

    def test_a_symmetric_three_run_slab_is_a_reflected_slab(self):
        """Built from thicknesses, as a registry states Sood's slabs: the
        right reflector's width carries the breakpoints' rounding, and the
        classification still reads it as symmetric."""
        g = StructuredGeometry.from_thicknesses(
            coord=_SLAB, thicknesses=(0.5, 0.8603, 0.5), mat_ids=(1, 0, 1),
            boundaries=(BC.vacuum, BC.vacuum),
        )
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
        g = StructuredGeometry.from_thicknesses(
            coord=_SLAB, thicknesses=thicknesses, mat_ids=mat_ids,
            boundaries=(BC.vacuum, BC.vacuum),
        )
        assert isinstance(reference_body(g), LayeredBody)

    def test_a_layered_solid_sphere(self):
        body = reference_body(_geometry(_SPH, (0.0, 1.0, 3.0), (0, 1)))
        assert isinstance(body, LayeredBody)
        assert body == LayeredBody(
            coord=_SPH, breakpoints=(0.0, 1.0, 3.0), mat_ids=(0, 1), is_hollow=False,
        )
        assert not body.is_hollow

    def test_a_layered_hollow_cylinder(self):
        g = _geometry(_CYL, (0.5, 1.0, 3.0), (0, 1), boundaries=(BC.reflective, BC.vacuum))
        body = reference_body(g)
        assert isinstance(body, LayeredBody) and body.is_hollow


# ─────────────────────────────────────────────────────────────────────
# Each owner serves exactly its shapes, and refuses the rest through the door
# ─────────────────────────────────────────────────────────────────────

_ONE_GROUP = {0: get_mixture("A", "1g"), 1: get_mixture("B", "1g")}


def _isotropic(key: str, groups: str = "1g"):
    """The library mixture with its P0 scattering only: ``Billiard`` refuses a higher moment (#405 P2 step 7b.2.3)."""
    xs = get_xs(key, groups)
    return make_mixture(sig_t=xs["sig_t"], sig_c=xs["sig_c"], sig_f=xs["sig_f"], nu=xs["nu"], chi=xs["chi"], sig_s=xs["sig_s"])


#: The 1-group materials for ``Billiard``, which solves isotropic scattering only.
_ONE_GROUP_ISOTROPIC = {0: _isotropic("A"), 1: _isotropic("B")}

_SHAPES = {
    "hollow-sphere": lambda: _geometry(_SPH, (0.5, 2.0), (0,), (BC.reflective, BC.vacuum)),
    "layered-sphere": lambda: _geometry(_SPH, (0.0, 1.0, 2.0), (0, 1)),
    "layered-slab": lambda: StructuredGeometry.from_thicknesses(
        coord=_SLAB, thicknesses=(1.0, 0.5), mat_ids=(0, 1), boundaries=(BC.vacuum, BC.vacuum),
    ),
    "reflected-slab": lambda: StructuredGeometry.from_thicknesses(
        coord=_SLAB, thicknesses=(0.5, 1.0, 0.5), mat_ids=(1, 0, 1), boundaries=(BC.vacuum, BC.vacuum),
    ),
    "hollow-layered-sphere": lambda: _geometry(
        _SPH, (0.5, 1.0, 2.0), (0, 1), (BC.reflective, BC.vacuum),
    ),
}

#: (owner, the shapes it refuses): the shapes it serves are tested by route below.
_REFUSALS = [
    ("Spectrum", lambda g: Spectrum(geometry=g, materials=_ONE_GROUP),
     ["hollow-sphere", "layered-sphere", "layered-slab", "reflected-slab", "hollow-layered-sphere"]),
    ("BasisSpace", lambda g: BasisSpace(geometry=g, materials=_ONE_GROUP),
     ["hollow-sphere", "layered-sphere", "layered-slab", "reflected-slab", "hollow-layered-sphere"]),
    ("MomentSpace", lambda g: MomentSpace(geometry=g, materials=_ONE_GROUP),
     ["hollow-sphere", "layered-sphere", "layered-slab", "hollow-layered-sphere"]),
    ("Billiard", lambda g: Billiard(geometry=g, materials=_ONE_GROUP_ISOTROPIC),
     ["layered-slab", "reflected-slab", "hollow-layered-sphere"]),
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
            MomentSpace(geometry=g, materials=materials)

    def test_several_groups_are_refused(self):
        g = _SHAPES["reflected-slab"]()
        materials = {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")}
        with pytest.raises(NotImplementedError, match="one-group reflected slab"):
            MomentSpace(geometry=g, materials=materials)

    def test_its_flux_reconstruction_is_refused(self):
        g = _SHAPES["reflected-slab"]()
        materials = {0: _unit_sigma_t_mixture(1.5), 1: _unit_sigma_t_mixture(0.9)}
        with pytest.raises(NotImplementedError, match="flux reconstruction"):
            MomentSpace(geometry=g, materials=materials).reconstruct_flux()


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

    g = StructuredGeometry.from_thicknesses(
        coord=_SLAB, thicknesses=(0.5, 0.8603, 0.5), mat_ids=(1, 0, 1),
        boundaries=(BC.vacuum, BC.vacuum),
    )
    materials = {0: _unit_sigma_t_mixture(1.5), 1: _unit_sigma_t_mixture(0.9)}
    solution = MomentSpace(geometry=g, materials=materials, fn_order=7).solve_critical()
    bare = solve_fn_slab_reflected_critical(
        c_core=float(materials[0].scattering_ratio[0]),
        c_reflector=float(materials[1].scattering_ratio[0]),
        reflector_half_thickness=0.5,
        n_modes=7,
    )
    assert solution.parameter_value == bare.tau_critical_mfp
    assert solution.parameter_kind == "core_half_thickness_mfp"


_SMALL_SPHERE = {"n_r": 8, "n_mu": 8, "n_traj_quad": 16}
_SMALL_CYLINDER = {"n_r": 6, "n_mu_axial": 4, "n_phi_az": 8, "n_traj_quad": 16}


#: A two-group fuel | moderator layering written out as arrays ([from, to]
#: scattering, the solvers' and ``make_mixture``'s one orientation). The
#: mixtures are built from these arrays, and the bare solver is handed the
#: SAME arrays directly, so the gate does not share production's stacking.
_FUEL = dict(sig_t=[0.5, 1.2], sig_s=[[0.30, 0.15], [0.0, 0.90]],
             sig_f=[0.01, 0.20], nu=[2.5, 2.5], chi=[1.0, 0.0])
_MODERATOR = dict(sig_t=[0.6, 1.5], sig_s=[[0.40, 0.18], [0.0, 1.45]],
                  sig_f=[0.0, 0.0], nu=[0.0, 0.0], chi=[0.0, 0.0])


def _mixture_of(xs):
    sig_t, sig_s, sig_f = (np.array(xs[k], dtype=float) for k in ("sig_t", "sig_s", "sig_f"))
    return make_mixture(
        sig_t=sig_t, sig_c=sig_t - sig_f - sig_s.sum(axis=1), sig_f=sig_f,
        nu=np.array(xs["nu"], dtype=float), chi=np.array(xs["chi"], dtype=float),
        sig_s=sig_s,
    )


@pytest.mark.parametrize(
    "coord, quadrature, solver_name",
    [
        (_SPH, _SMALL_SPHERE, "sphere_mr"),
        (_CYL, _SMALL_CYLINDER, "cylinder_mr"),
    ],
    ids=["layered-sphere", "layered-cylinder"],
)
def test_billiard_routes_a_layered_solid_body(coord, quadrature, solver_name):
    """#190's done-when: the facade's k is the bare multi-region solver's,
    bit for bit, on a two-group fuel / moderator body."""
    from orpheus.derivations.continuous.trajectory_resolvent import (
        greens_function as gf_sphere,
        greens_function_cylinder as gf_cyl,
    )

    materials = {0: _mixture_of(_FUEL), 1: _mixture_of(_MODERATOR)}
    g = _geometry(coord, (0.0, 1.0, 2.0), (0, 1), (BC.reflective,))
    b = Billiard(geometry=g, materials=materials, quadrature=dict(quadrature))
    assert b.geometry_kind == solver_name
    k_facade = b.solve_critical().eigenvalue
    bare = gf_sphere.solve_greens_function_sphere_mr if coord is _SPH else gf_cyl.solve_greens_function_cylinder_mr
    k_bare = bare(
        radii=np.array([1.0, 2.0]),
        sigma_t=np.array([_FUEL["sig_t"], _MODERATOR["sig_t"]]),
        sigma_s=np.array([_FUEL["sig_s"], _MODERATOR["sig_s"]]),
        nu_sigma_f=np.array(_FUEL["sig_f"]) * np.array(_FUEL["nu"])[None, :] * np.array([[1.0], [0.0]]),
        chi=np.array([_FUEL["chi"], _MODERATOR["chi"]]),
        alpha=1.0, **quadrature,
    ).k_eff
    assert k_facade == k_bare


@pytest.mark.parametrize(
    "coord, kind", [(_SPH, "hollow_sphere"), (_CYL, "annulus")], ids=["hollow-sphere", "annulus"],
)
def test_billiard_routes_a_hollow_body(coord, kind):
    """#421: the hollow arms are reachable from the constructor."""
    g = _geometry(coord, (0.5, 2.0), (0,), (BC.reflective, BC.vacuum))
    b = Billiard(geometry=g, materials={0: _isotropic("A")})
    assert b.geometry_kind == kind
    assert b.geometry_payload == {"R_in": 0.5, "R_out": 2.0}


def test_billiard_reads_the_body_material():
    """Billiard takes its cross sections from the body's material id, not
    from key 0 (qa, 2026-09-29: with mat_ids (3,) and materials {0: A, 3: C}
    it returned A's k_inf)."""
    def sphere(mat_id: int) -> StructuredGeometry:
        return _geometry(_SPH, (0.0, 2.0), (mat_id,), (BC.reflective,))

    decoy, body = _isotropic("A"), _isotropic("B")
    read = Billiard(geometry=sphere(3), materials={0: decoy, 3: body})
    reference = Billiard(geometry=sphere(0), materials={0: body})
    assert read.xs_payload == reference.xs_payload
    with pytest.raises(ValueError, match="the body's material id 3"):
        Billiard(geometry=sphere(3), materials={0: body})


# ─────────────────────────────────────────────────────────────────────
# The laws are part of what an owner serves (the user's ruling, 2026-09-29)
# ─────────────────────────────────────────────────────────────────────

from orpheus.derivations.common.reference_body import specular_albedo  # noqa: E402
from orpheus.geometry.boundary import (  # noqa: E402
    AlbedoBoundary,
    ReflectiveBoundary,
    SpecularReturn,
    VacuumInflow,
)


class TestSpecularAlbedo:
    @pytest.mark.parametrize(
        "law, albedo",
        [
            (BC.vacuum, 0.0),
            (VacuumInflow(), 0.0),
            (BC.reflective, 1.0),
            (ReflectiveBoundary(axis="x"), 1.0),
            (BC("partial", {"albedo": 0.3}), 0.3),
            (AlbedoBoundary(0.7, SpecularReturn(axis="x")), 0.7),
        ],
        ids=["vacuum-tag", "vacuum-law", "mirror-tag", "mirror-law", "partial-tag", "specular-albedo-law"],
    )
    def test_a_specular_law_reads_as_its_albedo(self, law, albedo):
        assert specular_albedo(law, owner="test") == albedo

    @pytest.mark.parametrize(
        "law", [BC("white"), AlbedoBoundary(0.7)], ids=["white", "unstated-reemission"],
    )
    def test_any_other_law_is_refused(self, law):
        with pytest.raises(NotImplementedError, match="another angular shape"):
            specular_albedo(law, owner="test")


def _reflected_slab(faces=(BC.vacuum, BC.vacuum), thicknesses=(0.5, 1.0, 0.5)):
    return StructuredGeometry.from_thicknesses(
        coord=_SLAB, thicknesses=thicknesses, mat_ids=(1, 0, 1), boundaries=faces,
    )


class TestTheLawsAreServed:
    def test_moment_space_refuses_a_reflected_slab_with_mirror_faces(self):
        """The F1 witness: before the laws joined the served pattern, a
        reflected slab with reflective faces returned the vacuum answer, bit
        for bit (the elegance review, 2026-09-29)."""
        materials = {0: _unit_sigma_t_mixture(1.5), 1: _unit_sigma_t_mixture(0.9)}
        with pytest.raises(NotImplementedError, match="reflecting boundaries"):
            MomentSpace(geometry=_reflected_slab((BC.reflective, BC.reflective)), materials=materials)

    @pytest.mark.parametrize(
        "owner, build",
        [("MomentSpace", lambda g: MomentSpace(geometry=g, materials=_ONE_GROUP)),
         ("BasisSpace", lambda g: BasisSpace(geometry=g, materials=_ONE_GROUP))],
    )
    def test_the_bare_solvers_refuse_a_reflecting_boundary(self, owner, build):
        g = _geometry(_SPH, (0.0, 2.0), (0,), (BC("partial", {"albedo": 0.5}),))
        with pytest.raises(NotImplementedError, match="reflecting boundaries") as caught:
            build(g)
        assert str(caught.value).startswith(owner)

    def test_spectrum_refuses_a_slab_whose_faces_differ(self):
        """Atalay's slab puts one R on both faces."""
        g = _geometry(_SLAB, (0.0, 2.0), (0,), (BC.vacuum, BC("partial", {"albedo": 0.5})))
        with pytest.raises(NotImplementedError, match="unequal face albedos"):
            Spectrum(geometry=g, materials=_ONE_GROUP)

    def test_spectrum_refuses_a_reflected_cylinder(self):
        g = _geometry(_CYL, (0.0, 2.0), (0,), (BC("partial", {"albedo": 0.5}),))
        with pytest.raises(NotImplementedError, match="reflected cylinder"):
            Spectrum(geometry=g, materials=_ONE_GROUP)

    def test_billiard_reads_a_slab_with_different_faces_as_two_surface(self):
        g = _geometry(_SLAB, (0.0, 2.0), (0,), (BC("partial", {"albedo": 0.7}), BC.reflective))
        b = Billiard(geometry=g, materials=_ONE_GROUP_ISOTROPIC)
        assert b.geometry_kind == "slab_asymmetric"
        assert b.alpha_payload == {"alpha_left": 0.7, "alpha_right": 1.0}
        assert b.closure_rank == 2

    def test_billiard_reads_a_hollow_body_s_two_laws(self):
        """qa's 2026-09-29 row: a declared (reflective, vacuum) hollow sphere
        is solved with those albedos, never with a separate parameter."""
        g = _geometry(_SPH, (0.5, 2.0), (0,), (BC.reflective, BC.vacuum))
        assert Billiard(geometry=g, materials=_ONE_GROUP_ISOTROPIC).alpha_payload == {
            "alpha_in": 1.0, "alpha_out": 0.0,
        }

    @pytest.mark.parametrize(
        "build",
        [lambda g: Spectrum(geometry=g, materials=_ONE_GROUP),
         lambda g: Billiard(geometry=g, materials=_ONE_GROUP_ISOTROPIC)],
        ids=["Spectrum", "Billiard"],
    )
    def test_a_white_law_is_refused(self, build):
        with pytest.raises(NotImplementedError, match="another angular shape"):
            build(_geometry(_SPH, (0.0, 2.0), (0,), (BC("white"),)))


class TestReflectedSlabRoute:
    def test_the_symmetry_tolerance_meets_a_real_rounding(self):
        """F9: from the thicknesses (0.3, 1.1, 0.3) the right reflector's
        width differs from the left's by a fraction of an ulp (the fold's
        rounding); the slab is still symmetric. Control: a 5-ulp asymmetry
        is a layered slab (the tolerance is 4 ulp)."""
        import math

        g = _reflected_slab(thicknesses=(0.3, 1.1, 0.3))
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
            geometry=_reflected_slab(thicknesses=(0.5, 1.0, 0.5)),
            materials={0: _unit_sigma_t_mixture(1.5), 1: _unit_sigma_t_mixture(0.9)},
        ).solve_critical()
        doubled = MomentSpace(
            geometry=_reflected_slab(thicknesses=(0.25, 0.5, 0.25)),
            materials={0: _unit_sigma_t_mixture(1.5, 2.0), 1: _unit_sigma_t_mixture(0.9, 2.0)},
        ).solve_critical()
        assert doubled.metadata["reflector_half_thickness_mfp"] == 0.5
        assert doubled.parameter_value == unit.parameter_value


@pytest.mark.catches("ERR-091")
def test_billiard_sphere_mr_fixed_source_reports_every_group():
    """ERR-091: the multi-region sphere's fixed-source arm read the group
    count as ``reshape(-1, 1).shape[1]`` (always 1) and returned group 0 as
    the scalar flux. On a two-group source the count is 2 and the scalar
    flux is the total, the sum over groups."""
    materials = {0: _mixture_of(_FUEL), 1: _mixture_of(_MODERATOR)}
    g = _geometry(_SPH, (0.0, 1.0, 2.0), (0, 1), (BC.vacuum,))
    b = Billiard(geometry=g, materials=materials, quadrature=dict(_SMALL_SPHERE))
    flux = b.solve_fixed_source(external_source=np.array([[1.0, 0.5], [0.0, 0.0]]))
    phi_g = flux.metadata["phi_g"]
    assert flux.metadata["n_groups"] == 2
    assert phi_g.shape[0] == 2
    np.testing.assert_array_equal(flux.scalar_flux, phi_g.sum(axis=0))
