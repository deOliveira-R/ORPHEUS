"""DRY-RUN of the proposed B3.4c gates, verbatim from the plan.

Run from the repo root with:
    .venv/bin/python -O -m pytest /Users/rodrigo/.claude/jobs/c30e4f25/tmp/test_b34c_proposed.py -q

Expected TODAY (step 4 NOT done): the C4/C5/C7 rows RED, the C1/C2/C3/C6/E.2
rows GREEN.
"""
from __future__ import annotations

import itertools
from dataclasses import replace

import numpy as np
import pytest

from orpheus.geometry.boundary import PeriodicBoundary
from orpheus.geometry.boundary._errors import BoundaryError
from orpheus.geometry.boundary._factors import (
    IdentityMap, SpatialWrap, SpecularMirror,
)
from orpheus.geometry.mesh import Mesh1D
from orpheus.numerics.face_layout import (
    AXIS_NAMES, FACE_NAMES, face_name, face_normal, face_opposite,
)
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.boundary.realizer import SNBoundaryRealizer
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.operators.boundary import SNBoundaryOperator
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.timed_full_field import TimedFullField
from tests.sn._test_helpers import face_method_space, placeholder_materials

pytestmark = [pytest.mark.foundation]

_AXES = range(len(AXIS_NAMES))
_SIGNS = (-1, +1)


# ── C1 ────────────────────────────────────────────────────────────────
@pytest.mark.parametrize("axis,sign", list(itertools.product(_AXES, _SIGNS)))
def test_render_then_parse_is_the_identity(axis, sign):
    assert face_normal(face_name(axis, sign)) == (axis, sign)


@pytest.mark.parametrize("face", FACE_NAMES)
def test_parse_then_render_is_the_identity(face):
    assert face_name(*face_normal(face)) == face


def test_the_inventory_is_the_whole_product():
    assert set(FACE_NAMES) == {f"{a}{s}" for a in AXIS_NAMES for s in ("min", "max")}
    assert len(FACE_NAMES) == len(set(FACE_NAMES)) == 2 * len(AXIS_NAMES)


@pytest.mark.parametrize("axis,sign", list(itertools.product(_AXES, _SIGNS)))
def test_the_bijection_reproduces_the_retired_transcriptions(axis, sign):
    legacy = f"{AXIS_NAMES[axis]}{'max' if sign == +1 else 'min'}"
    assert face_name(axis, sign) == legacy
    assert face_normal(legacy) == (axis, sign)


@pytest.mark.parametrize("bad", ["wmin", "xmid", "x", "", "minx", "xminmax"])
def test_face_normal_refuses_a_non_canonical_name(bad):
    with pytest.raises(ValueError):
        face_normal(bad)


# ── C2 ────────────────────────────────────────────────────────────────
@pytest.mark.parametrize("face", FACE_NAMES)
def test_opposite_is_an_involution(face):
    assert face_opposite(face_opposite(face)) == face


@pytest.mark.parametrize("face", FACE_NAMES)
def test_opposite_preserves_the_axis_and_flips_the_normal(face):
    axis, sign = face_normal(face)
    assert face_normal(face_opposite(face)) == (axis, -sign)


@pytest.mark.parametrize("face", FACE_NAMES)
def test_opposite_has_no_fixed_point(face):
    assert face_opposite(face) != face


# ── C2b ───────────────────────────────────────────────────────────────
@pytest.mark.parametrize("axis,ep,sign",
                        [(a, ep, s) for a in _AXES for ep, s in (("min", -1), ("max", +1))])
def test_facelabel_render_agrees_with_the_bijection(axis, ep, sign):
    from orpheus.transport.mesh.axis import FaceLabel
    assert FaceLabel(axis, ep).face_name == face_name(axis, sign)


# ── C3 ────────────────────────────────────────────────────────────────
@pytest.mark.parametrize("gmap", [IdentityMap(), SpecularMirror(axis="x"),
                                  SpecularMirror(axis="y")])
@pytest.mark.parametrize("face", FACE_NAMES)
def test_a_non_quotient_map_consumes_the_face_it_is_installed_on(gmap, face):
    assert gmap.domain_face(face) == face


@pytest.mark.parametrize("axis", AXIS_NAMES)
def test_the_wrap_consumes_the_OPPOSITE_face(axis):
    for sign in (-1, +1):
        face = face_name(AXIS_NAMES.index(axis), sign)
        assert SpatialWrap(axis=axis).domain_face(face) == face_opposite(face)
        assert SpatialWrap(axis=axis).domain_face(face) != face


@pytest.mark.parametrize("wrap_axis,face", [
    ("y", "xmin"), ("y", "xmax"), ("x", "ymin"), ("x", "ymax"),
    ("z", "xmin"), ("x", "zmax"),
])
def test_a_wrap_on_a_foreign_axis_REFUSES(wrap_axis, face):
    with pytest.raises(BoundaryError) as exc:
        SpatialWrap(axis=wrap_axis).domain_face(face)
    assert exc.value.law == "periodic"
    assert wrap_axis in str(exc.value) and face in str(exc.value)


def test_an_unparseable_face_REFUSES_as_a_BoundaryError():
    with pytest.raises(BoundaryError):
        SpatialWrap(axis="x").domain_face("bogus")


def test_every_geometry_map_answers_domain_face():
    for gmap in (IdentityMap(), SpecularMirror(), SpatialWrap()):
        assert isinstance(gmap.domain_face("xmax"), str)


# ── C6 ────────────────────────────────────────────────────────────────
def test_the_wrap_declares_an_honest_transpose():
    assert SpatialWrap(axis="x").is_adjointable


def test_every_deck_transformation_is_adjointable():
    for gmap in (IdentityMap(), SpecularMirror(), SpatialWrap()):
        assert gmap.is_adjointable


# ── E.2 ───────────────────────────────────────────────────────────────
def test_the_identification_compares_SETS_not_sizes():
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    space = face_method_space(quad, face="xmin", faces=("xmin", "xmax"))
    poisoned = replace(space, inflow_indices=np.asarray(space.outflow_indices))
    assert np.asarray(poisoned.inflow_indices).size == np.asarray(space.inflow_indices).size
    with pytest.raises(BoundaryError) as exc:
        SNBoundaryRealizer().realize(PeriodicBoundary(axis="x"), poisoned)
    assert exc.value.law == "periodic"


def test_the_honest_space_realizes():
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    space = face_method_space(quad, face="xmin", faces=("xmin", "xmax"))
    SNBoundaryRealizer().realize(PeriodicBoundary(axis="x"), space)


# ── C4 / C5 / C7 fixtures ─────────────────────────────────────────────
def _periodic_slab(n_ordinates: int = 8) -> SNMesh:
    quad = Quadrature.gauss_legendre(n_ordinates=n_ordinates)
    mesh = Mesh1D(edges=np.linspace(0.0, 2.0, 5), mat_ids=np.zeros(4, dtype=int))
    sn = SNMesh(mesh, quad, placeholder_materials())
    for face in ("xmin", "xmax"):
        sn.bc[face] = sn.realize_boundary_law(PeriodicBoundary(axis="x"), face)
    return sn


def _independently_seeded_trace(sn: SNMesh, seed: int):
    z = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn)
    rng = np.random.default_rng(seed)
    return replace(z, boundary=replace(
        z.boundary, values=rng.uniform(0.5, 2.0, size=z.boundary.values.shape)))


class TestPeriodicIsACrossFaceCoupling:
    def test_the_two_faces_carry_different_data(self):
        sn = _periodic_slab()
        psi = _independently_seeded_trace(sn, seed=20260801)
        tr = sn.angular_trace
        assert not np.array_equal(
            psi.boundary.face_view("xmin")[tr.outflow_indices_for_face("xmin")],
            psi.boundary.face_view("xmax")[tr.outflow_indices_for_face("xmax")],
        )

    @pytest.mark.parametrize("face,partner", [("xmin", "xmax"), ("xmax", "xmin")])
    def test_the_inflow_is_the_PARTNERS_outflow(self, face, partner):
        sn = _periodic_slab()
        psi = _independently_seeded_trace(sn, seed=20260801)
        tr = sn.angular_trace
        got = SNBoundaryOperator(sn).apply(psi).boundary.face_view(face)[
            tr.inflow_indices_for_face(face)]
        want = psi.boundary.face_view(partner)[tr.outflow_indices_for_face(partner)]
        np.testing.assert_array_equal(got, want)

    @pytest.mark.parametrize("face", ["xmin", "xmax"])
    def test_it_is_NOT_the_faces_own_outflow(self, face):
        sn = _periodic_slab()
        psi = _independently_seeded_trace(sn, seed=20260801)
        tr = sn.angular_trace
        got = SNBoundaryOperator(sn).apply(psi).boundary.face_view(face)[
            tr.inflow_indices_for_face(face)]
        own = psi.boundary.face_view(face)[tr.outflow_indices_for_face(face)]
        assert not np.array_equal(got, own)

    @pytest.mark.parametrize("face", ["xmin", "xmax"])
    def test_nothing_lands_off_gamma_minus(self, face):
        sn = _periodic_slab()
        psi = _independently_seeded_trace(sn, seed=20260801)
        tr = sn.angular_trace
        out = SNBoundaryOperator(sn).apply(psi).boundary.face_view(face)
        off = np.setdiff1d(np.arange(out.shape[0]), tr.inflow_indices_for_face(face))
        np.testing.assert_array_equal(out[off], np.zeros_like(out[off]))


@pytest.mark.parametrize("n_ordinates", [4, 8])
def test_euclidean_reciprocity_on_a_periodic_slab(n_ordinates):
    sn = _periodic_slab(n_ordinates=n_ordinates)
    B = SNBoundaryOperator(sn)
    x = _independently_seeded_trace(sn, seed=1)
    y = _independently_seeded_trace(sn, seed=2)
    lhs = float(np.sum(B.apply(x).boundary.values * y.boundary.values))
    rhs = float(np.sum(x.boundary.values * B.apply_transpose(y).boundary.values))
    n_terms = x.boundary.values.size
    np.testing.assert_allclose(
        lhs, rhs, rtol=0.0,
        atol=n_terms * float(np.finfo(np.float64).eps) * max(abs(lhs), abs(rhs)))


def test_the_transpose_scatters_over_the_PARTNERS_outflow():
    sn = _periodic_slab()
    tr = sn.angular_trace
    z = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn)
    y = replace(z, boundary=replace(z.boundary, values=z.boundary.values.copy()))
    y.boundary.face_view("xmin")[tr.inflow_indices_for_face("xmin")] = 1.0
    out = SNBoundaryOperator(sn).apply_transpose(y).boundary
    target = out.face_view("xmax")[tr.outflow_indices_for_face("xmax")]
    np.testing.assert_array_equal(target, np.ones_like(target))
    np.testing.assert_array_equal(out.face_view("xmin"),
                                  np.zeros_like(out.face_view("xmin")))


def test_B_is_block_structured_not_block_diagonal():
    """ZERO BASE — the differencing form is NOT bit-exact ([M] 1 ULP)."""
    sn = _periodic_slab()
    tr = sn.angular_trace
    z = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn)
    bumped = replace(z, boundary=replace(z.boundary, values=z.boundary.values.copy()))
    bumped.boundary.face_view("xmax")[tr.outflow_indices_for_face("xmax")] = 1.0
    out = SNBoundaryOperator(sn).apply(bumped).boundary
    sel = out.face_view("xmin")[tr.inflow_indices_for_face("xmin")]
    np.testing.assert_array_equal(sel, np.ones_like(sel))
    off = np.setdiff1d(np.arange(out.face_view("xmin").shape[0]),
                       tr.inflow_indices_for_face("xmin"))
    np.testing.assert_array_equal(out.face_view("xmin")[off],
                                  np.zeros_like(out.face_view("xmin")[off]))
    # DIRECTION: xmax's Γ₋ reads xmin's Γ₊, which is zero here.
    np.testing.assert_array_equal(out.face_view("xmax"),
                                  np.zeros_like(out.face_view("xmax")))


# ── C8 (re-posed: the refusal's home is the COMPOSITE) ────────────────
def test_a_wrap_whose_partner_face_is_absent_is_REFUSED_at_the_COMPOSITE():
    """[M] realization is SILENT on a one-faced mesh; `_face_domains`'s
    permutation certification is what catches it."""
    from orpheus.geometry.mesh import CoordSystem
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    sph = Mesh1D(edges=np.linspace(0.1, 1.0, 6), mat_ids=np.zeros(5, dtype=int),
                 coord=CoordSystem.SPHERICAL)
    sn = SNMesh(sph, quad, placeholder_materials())
    assert tuple(sn.angular_trace.layout.faces) == ("xmax",)
    sn.bc["xmax"] = sn.realize_boundary_law(PeriodicBoundary(axis="x"), "xmax")
    with pytest.raises(ValueError):
        SNBoundaryOperator(sn)._face_domains


# ── B.5 re-pose ───────────────────────────────────────────────────────
def test_the_composition_supplies_the_partner_half_trace():
    from pathlib import Path

    from tests.geometry import _generate_bc_equivalence_snapshots as _gen
    case = _gen.case_by_id("periodic_lebedev17")
    snap = np.load(Path(_gen.__file__).parent / "snapshots"
                   / "bc_equivalence_periodic_lebedev17.npz")
    space = face_method_space(case.build_quadrature(), face=case.face, faces=case.faces)
    law = PeriodicBoundary(axis="x")
    assert law.geometry_map.domain_face(case.face) == str(snap["partner_face"])
    source_face = law.geometry_map.domain_face(case.face)
    probe = {case.face: snap["psi_out"], source_face: snap["psi_out_partner"]}[source_face]
    actual = SNBoundaryRealizer().realize(law, space).apply(probe)
    np.testing.assert_array_equal(actual, snap["psi_in"])
