r"""The carry of wave 2 for SN: every number the carve must not move (law L-K; plan sections 6.0 and 8; S2b gate 3).

``tests/gates/transport/carry_head.npz`` was captured FROM THE HEAD TREE (``main`` ``5d04c296``, a detached worktree,
before any edit) by ``tests/gates/transport/_capture_carry_head.py``: for the SN bodies (a slab, a cylinder, a sphere and a
hollow cylinder, two regions, two groups) the eigenvalue ``k`` and the scalar flux of ``solve_sn``, the scalar flux of
``solve_sn_fixed_source`` under a bulk source, the slab's fixed-source flux under a uniform inflow, the dense matrix of
the realized operator of both faces of an S16 slab; and, for each of the seven ``bc_equivalence`` snapshots, the image the
HEAD realizer gives of the frozen probe. The same problems, spelled in the NEW vocabulary (a tag bound to a
``WallCondition``, or a deck: ``tests/gates/sn/architecture/_sn_consumers_helpers.py``), must give the same answers
BIT-EQUAL (``array_equal``; a non-bit difference is a finding to investigate, never a tolerance to widen). The declared
exceptions are exactly those of section 8: R-11 (the Lambertian divisor, bit-identical on every measured pair: the mass
census of ``test_realized_boundary.py``) and nothing else.

The fixture also carries the other methods' ``k`` (diffusion on three charts, CP, MoC, MC with the fixed seed): slices S3
and S4 gate them. This file checks only that they were captured (``test_the_fixture_is_the_head_capture``).

What each row can and cannot see
================================

* The solves read every layer between the declaration and the answer: the conditions package, the binding, the realizer
  table, the typed operators, the sweep schedule, the corner swap, the leakage term (``compute_keff``), DSA. A move of the
  mirror's ordinate permutation, of the Lambertian divisor, of a return fraction, or of the face a deck reads moves a
  ``k`` here. Blind: a defect that is zero on these 2-group two-region bodies (a one-group problem would be blind to the
  whole spatial operator, ``vv-principles`` anti-pattern 3; a homogeneous one to redistribution; hence two regions with a
  fuel and a moderator).
* The face matrices see the realized operator alone (no sweep): they are the finest witness of the typed operators and
  localise a red to the realizer.
* The snapshots are independent derivations from the law's equation (``_generate_bc_equivalence_snapshots.py``); the
  fixture's recorded sha256 proves they were NOT re-captured by this carve (R-10).

First red at the base tip: ``ImportError`` from ``import orpheus.sn`` for every row that executes SN; the two
fixture rows (``test_the_fixture_is_the_head_capture``, ``test_the_snapshots_were_not_rebaselined``) read files only and
are green today (they pin the instrument, not the carve).

Mutations for the reviewers: an orientation-blind Lambertian (reddens the white rows and the ``white_*`` snapshots); a
mirror pairing taken from the wrong face (reddens every specular row); a divisor by the OUTFLOW mass on a rule with
``M_+ != M_-`` (none is shipped, so it reddens nothing here: the asymmetric row of ``test_realized_boundary.py`` is its
witness); ``SNSpecularChain`` reporting ``alpha = 1.0`` (reddens the partial-specular ``k`` rows); a leakage term that
sums every face (the ERR-064/094 family).
"""

from __future__ import annotations

import hashlib
from functools import cache
from typing import Any

import numpy as np
import pytest

from tests.gates.sn.architecture import _sn_consumers_helpers as S
from tests.gates.transport import _carry_cases as carry

_THIS = "tests/gates/sn/architecture/test_wave2_carry.py"
_CASES = carry.declared_cases()
_IDS = [f"{body}-{carry.tag(d)}" for body, d in _CASES]

CARRY_FIXTURE = carry.CARRY_FIXTURE
SNAPSHOT_DIR = carry.SNAPSHOT_DIR


@cache
def _fixture() -> Any:
    assert CARRY_FIXTURE.exists(), f"{CARRY_FIXTURE} is missing: the carry fixture is the instrument of this file"
    return np.load(CARRY_FIXTURE)


# ═════════════════════════════════════════════════════════════════════════════
# The instrument
# ═════════════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_the_fixture_is_the_head_capture() -> None:
    r"""The fixture was captured at ``5d04c296`` and holds its populations: 17 eigenvalues and 17 + 17 flux arrays of the
    SN bodies, the slab inflow flux, 14 face matrices, 7 snapshot images, 12 diffusion, 6 CP, 1 MoC and 1 MC eigenvalues.
    The positive control is physical: on every body the eigenvalue is ordered vacuum < specular 0.3 < specular 0.7 <
    mirror, so a fixture of a flat or inverted column (a broken capture) fails here.

    Red input: none today (the capture is complete); it reddens on a fixture captured from the wrong tree or truncated.
    """
    fixture = _fixture()
    assert str(fixture["meta/commit"]) == carry.HEAD_COMMIT
    names = list(fixture.files)

    def count(prefix: str) -> int:
        return sum(1 for n in names if n.startswith(prefix))

    assert (count("keff/"), count("phi_eigen/"), count("phi_fixed/")) == (17, 17, 17 + 1)
    assert count("face_matrix/") == 2 * len(carry.DECLARATIONS["slab"]) and count("snapshot_image/") == 7 and count("snapshot_sha256/") == 7
    assert (count("keff_diffusion/"), count("keff_cp/"), count("keff_moc/"), count("keff_mc/")) == (12, 6, 1, 1)
    for body in ("slab", "cylinder", "sphere"):
        ordered = [float(fixture[carry.key(body, d, "keff")]) for d in (("vacuum",), ("specular", 0.3), ("specular", 0.7), ("mirror",))]
        assert ordered == sorted(ordered) and len(set(ordered)) == 4, f"{body}: k is not ordered vacuum < 0.3 < 0.7 < mirror: {ordered}"


@pytest.mark.foundation
def test_the_snapshots_were_not_rebaselined() -> None:
    r"""The seven ``bc_equivalence_*.npz`` references are the carry fixture's independent half (R-10): their sha256 equals
    the one recorded at HEAD. A re-baseline here would be a carry failure dressed as a fix.

    Red input: none today. It reddens on a regenerated snapshot (the generator is not to be re-run to make a red go away).
    """
    fixture = _fixture()
    files = sorted(SNAPSHOT_DIR.glob("bc_equivalence_*.npz"))
    assert len(files) == 7, f"{len(files)} snapshot files under {SNAPSHOT_DIR}"
    for path in files:
        case_id = path.stem.removeprefix("bc_equivalence_")
        assert hashlib.sha256(path.read_bytes()).hexdigest() == str(fixture[f"snapshot_sha256/{case_id}"]), f"{path.name} differs from the HEAD snapshot"


# ═════════════════════════════════════════════════════════════════════════════
# The solves
# ═════════════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.rests_on(f"{_THIS}::test_the_fixture_is_the_head_capture")
@pytest.mark.parametrize("body,declaration", _CASES, ids=_IDS)
def test_the_carry_of_the_solves(body: str, declaration: tuple) -> None:
    r"""The k-eigenvalue and the scalar flux of ``solve_sn``, and the scalar flux of ``solve_sn_fixed_source`` under the
    bulk source of the fixture, equal the HEAD answers bit-for-bit, for each of the 17 (body, declaration) pairs: vacuum,
    mirror (a deck on the slab, ``specular(1.0)`` on a curved body), partial specular 0.3 and 0.7, isotropic 0.5 and 1.0
    (Cartesian only), and the hollow cylinder's reflective cavity under a vacuum and a 0.7 specular outer wall.

    Red input: at the base tip ``boundary_conditions`` is not a keyword of any entry. A defect that moves any layer between
    the declaration and the answer reddens it; the blind spots are in the module docstring.
    """
    from orpheus.sn.solver import solve_sn, solve_sn_fixed_source

    fixture = _fixture()
    mesh, quad, conditions = S.build(body, declaration)
    mats = carry.materials()
    eigen = solve_sn(mats, mesh, quad, boundary_conditions=conditions, **carry.EIGEN)
    assert eigen.outcome.keff == float(fixture[carry.key(body, declaration, "keff")]), f"{body}/{carry.tag(declaration)}: k moved"
    np.testing.assert_array_equal(np.asarray(eigen.scalar_flux.values), fixture[carry.key(body, declaration, "phi_eigen")])

    problem = S.problem(mesh, quad, conditions)
    source = carry.bulk_source(quad.n_ordinates, problem.ng, tuple(problem.spatial_shape))
    fixed = solve_sn_fixed_source(mats, mesh, quad, source, boundary_conditions=conditions, **carry.FIXED)
    np.testing.assert_array_equal(np.asarray(fixed.scalar_flux.values), fixture[carry.key(body, declaration, "phi_fixed")])


@pytest.mark.l1
@pytest.mark.rests_on(f"{_THIS}::test_the_fixture_is_the_head_capture")
def test_the_carry_of_the_slab_inflow() -> None:
    r"""The slab fixed-source solve under a UNIFORM inflow of ``1.0`` on the left face (tag 0 bound to
    ``WallCondition.uniform(1.0)``, the right face a vacuum tag) equals the HEAD flux bit-for-bit: the boundary-source
    channel end to end (the binding's ``sources``, ``from_binding``, the sweep's inflow seed).

    Red input: at the base tip the keyword does not exist. Mutation: a ``from_binding`` that drops a face's source.
    """
    from orpheus.sn.solver import solve_sn_fixed_source

    B = S.bc()
    fixture = _fixture()
    body = carry.BODIES["slab"]
    mesh = S.mesh_of("slab", left=0, right=1)
    quad = S.quadrature_of("slab")
    conditions = {0: B.WallCondition.uniform(carry.INFLOW_VALUE), 1: B.WallCondition.vacuum()}
    problem = S.problem(mesh, quad, conditions)
    zero = np.zeros((quad.n_ordinates, problem.ng, *problem.spatial_shape))
    fixed = solve_sn_fixed_source(carry.materials(), mesh, quad, zero, boundary_conditions=conditions, **carry.FIXED)
    flux = np.asarray(fixed.scalar_flux.values)
    assert body.name == "slab" and float(np.max(flux)) > 0.0, "the control: the inflow must populate the slab"
    np.testing.assert_array_equal(flux, fixture[carry.inflow_key("phi_fixed")])


@pytest.mark.l1
@pytest.mark.rests_on(f"{_THIS}::test_the_fixture_is_the_head_capture")
@pytest.mark.parametrize("face", ("xmin", "xmax"))
@pytest.mark.parametrize("declaration", carry.DECLARATIONS["slab"], ids=[carry.tag(d) for d in carry.DECLARATIONS["slab"]])
def test_the_carry_of_the_face_operators(declaration: tuple, face: str) -> None:
    r"""The dense matrix of the realized operator of each face of an S16 slab equals the HEAD matrix bit-for-bit, for every
    slab declaration: the finest witness of the typed operators, with no sweep between the realizer and the number. A red
    here and green solves is impossible; a red solve with green matrices localises the defect AFTER the realizer.

    Red input: at the base tip ``boundary_conditions`` is not a keyword. Mutation: a mirror permutation read from the
    opposite face; the Lambertian divisor changed.
    """
    from tests.gates.sn.architecture import _realized_boundary_helpers as H

    fixture = _fixture()
    roles, conditions = S.roles_and_conditions("slab", declaration)
    mesh = S.mesh_of("slab", **roles)
    quad = S.quadrature_of("slab", carry.S16_ORDINATES)
    problem = S.problem(mesh, quad, conditions)
    n = int(len(problem.angular_trace.outflow_indices_for_face(face)))
    np.testing.assert_array_equal(H.dense(problem.bc[face], n), fixture[carry.face_key(declaration, face)])


# ═════════════════════════════════════════════════════════════════════════════
# The seven equivalence snapshots (R-10): the bit-equal carry of the realizer's images
# ═════════════════════════════════════════════════════════════════════════════

_FOUR = ("xmin", "xmax", "ymin", "ymax")
_TWO = ("xmin", "xmax")


def _snapshot_cases() -> dict[str, tuple]:
    """``case_id -> (quadrature builder, face, faces, declared builder)``, typed from the snapshot registry's descriptions."""
    from orpheus.geometry.deck import TranslationDeck

    def q() -> Any:
        from orpheus.numerics.quadrature import Quadrature

        return Quadrature

    def w() -> Any:
        return S.bc().WallCondition

    return {
        "vacuum_lebedev17": (lambda: q().lebedev(17), "xmin", _TWO, lambda: w().vacuum()),
        "specular_x_lebedev17": (lambda: q().lebedev(17), "xmax", _FOUR, lambda: w().specular(1.0)),
        "specular_y_partial_07_LS6": (lambda: q().level_symmetric(sn_order=6), "ymax", _FOUR, lambda: w().specular(0.7)),
        "white_xmax_LS4": (lambda: q().level_symmetric(sn_order=4), "xmax", _FOUR, lambda: w().white(1.0)),
        "white_xmin_partial_03_GL": (lambda: q().gauss_legendre(n_ordinates=8), "xmin", _TWO, lambda: w().white(0.3)),
        "mixed_30spec_70white_LS4": (lambda: q().level_symmetric(sn_order=4), "xmax", _FOUR, lambda: w().from_amplitudes(0.3, 0.7)),
        "periodic_lebedev17": (lambda: q().lebedev(17), "xmin", _TWO, lambda: TranslationDeck()),
    }


@pytest.mark.l1
@pytest.mark.rests_on(f"{_THIS}::test_the_snapshots_were_not_rebaselined")
@pytest.mark.parametrize("case_id", list(_snapshot_cases()))
def test_the_equivalence_snapshots_carry(case_id: str) -> None:
    r"""Each of the seven snapshots, built on the NEW types through the real SN table (the mixed one is the ``RealizedSum``
    witness: ``from_amplitudes(0.3, 0.7)``, the periodic one a ``TranslationDeck`` read on the PARTNER face's outflow),
    gives the image the HEAD realizer gave of the frozen probe BIT-EQUAL, and agrees with the independent frozen reference
    ``psi_in`` at the old derived tolerance: ``array_equal`` for the vacuum, specular and periodic cases; ``|Gamma_+| eps``
    for the white cases and ``(|Gamma_+| + 2) eps`` for the mixed one (a reduction of positive terms, condition number 1).

    The specular cases are also checked as a mirror DECK on the Cartesian face (rule R2: a reflective declaration on a
    Cartesian face is a deck): the same image, bit-equal.

    Red input at the base tip: the table does not exist. Mutation: the Lambertian taken over the wrong hemisphere (white
    rows); the periodic image read on the installation face's own outflow (the probe routes by ``domain_face``).
    """
    from orpheus.geometry.deck import MirrorDeck
    from tests.gates.sn.architecture import _realized_boundary_helpers as H

    build_quad, face, faces, build_declared = _snapshot_cases()[case_id]
    quad = build_quad()
    snapshot = np.load(SNAPSHOT_DIR / f"bc_equivalence_{case_id}.npz")
    trace = H.trace_of(quad, faces)
    space = H.method_space(quad, face, trace)
    where = H.site(face, dimension=len(faces) // 2)
    declared = build_declared()
    op = H.sn_table().realize(declared, where, space)
    probe = snapshot["psi_out_partner"] if op.domain_face != face else snapshot["psi_out"]
    image = np.asarray(op.apply(probe))
    np.testing.assert_array_equal(image, _fixture()[f"snapshot_image/{case_id}"], err_msg=f"{case_id}: the image moved against HEAD")

    eps = float(np.finfo(np.float64).eps)
    reference = snapshot["psi_in"]
    n_plus = int(snapshot["outflow_indices"].size)
    if case_id.startswith(("white",)):
        np.testing.assert_allclose(image, reference, rtol=n_plus * eps, atol=0.0)
    elif case_id.startswith("mixed"):
        np.testing.assert_allclose(image, reference, rtol=(n_plus + 2) * eps, atol=0.0)
    else:
        np.testing.assert_array_equal(image, reference)

    if case_id.startswith("specular_x"):
        deck_image = np.asarray(H.sn_table().realize(MirrorDeck(), where, space).apply(probe))
        np.testing.assert_array_equal(deck_image, image)
