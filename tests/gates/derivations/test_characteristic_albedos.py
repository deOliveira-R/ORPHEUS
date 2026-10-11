"""How the characteristic reference's eigenvalue answers to the walls' albedos: the method of images and the ordering theorem.

P1 step (e1b) of the characteristic-reference campaign
(``.claude/plans/characteristic_reference_architecture.md``, "Step (e): the audit
and four rulings"): the successors, on the characteristic reference, of the old
trajectory-resolvent family's albedo rows, which step (e2) deleted with the
family. Verification spec ``scratch/characteristic_architecture/p1_verification_spec.md``,
rows B4c (the method of images) and D8 (the ordering theorem). The row-by-row
map from the old rows is ``scratch/characteristic_architecture/p1_step_e/ta_e1b/README.md``.

The ladder, bottom up (``rests_on`` on each row):

1. the closed bodies read k_inf (D5, ``test_characteristic_system.py``) and the
   door reads the system's k (``test_characteristic_reference.py``);
2. the method of images [B4c]: an albedo-1 wall is the symmetry plane of the
   doubled body, an identity between two specifications of the same reference
   that share no wall law;
3. the ordering theorem [D8]: k is strictly increasing in each wall's albedo,
   from the vacuum body to k_inf, and tends to k_inf as the body thickens.

Every row reads through the door (:class:`~orpheus.derivations.continuous.characteristic.CharacteristicDerivation`)
in this process (:func:`~orpheus.numerics.traced_memo.bypass`), at rungs of the
joint ladder (:func:`tests.gates.derivations._characteristic_ladders.rung`).
Bands are measured (``[M]`` 2026-10-10, ``.venv/bin/python -O``, probes in
``scratch/characteristic_architecture/p1_step_e/ta_e1b/probes/``); each row's first
red is an arm of ``scratch/characteristic_architecture/p1_step_e/ta_e1b/battery/``.
"""
from __future__ import annotations

from functools import lru_cache
from typing import Any

import numpy as np
import pytest

from orpheus.data.cells import CellCoefficient, Channel
from orpheus.data.materials import Materials
from orpheus.derivations.continuous.characteristic import CharacteristicDerivation, Resolution
from orpheus.geometry import BC, StructuredGeometry
from orpheus.geometry.boundary import ReflectiveBoundary
from orpheus.numerics.observable import Eigenvalue, PointValue
from orpheus.numerics.question import Eigen
from orpheus.numerics.traced_memo import bypass
from orpheus.reference.reading import Uncertified
from orpheus.specification.specification import GeometrySpecification
from tests.gates.derivations._characteristic_ladders import rung
from tests.gates.derivations.test_characteristic_system import _PU2, _mixture

pytestmark = pytest.mark.filterwarnings("error::RuntimeWarning")

_HERE = "tests/gates/derivations/test_characteristic_albedos.py::"
_SYSTEM = "tests/gates/derivations/test_characteristic_system.py::"
_REFERENCE = "tests/gates/derivations/test_characteristic_reference.py::"
_D5 = _SYSTEM + "test_a_closed_body_reads_k_inf_with_a_flat_flux_in_the_infinite_mediums_group_ratio"
_DOOR_K = _REFERENCE + "test_the_doors_k_is_the_fundamental_of_the_system_written_by_hand"

_K = CellCoefficient.every(Channel.FISSION_EMISSION)


@lru_cache(maxsize=None)
def _derivation(specification: GeometrySpecification, resolution: Resolution) -> CharacteristicDerivation:
    return CharacteristicDerivation(specification, resolution)


def _read(specification: GeometrySpecification, observable: Any, resolution: Resolution) -> float:
    with bypass():
        reading = _derivation(specification, resolution).evaluate(observable)
    assert type(reading) is Uncertified
    return float(reading.value)


def _eigen(geometry: StructuredGeometry, mixture) -> GeometrySpecification:
    return GeometrySpecification(Materials({m: mixture for m in set(geometry.mat_ids)}), geometry, Eigen(_K))


# ── B4c, the method of images ────────────────────────────────────────────

#: The old rows' one-group medium (Sigma_t 0.5, Sigma_s 0.38, nu Sigma_f 0.025: k about 0.0585 on [0, 2]), and
#: PU-2-0's two-group mixture (downscatter, chi in both groups).
_IMAGE_MIXTURES = {"1g": _mixture([0.5], [[0.38]], None, [0.025], [1.0]), "PU2": _PU2}
_HALF = 1.0
#: Rung 4 of the joint ladder. ``[M]`` 2026-10-10 (``probes/probe_moi.py``): the relative k difference between the
#: two specifications is 3.5e-8 (1G) and 3.6e-8 (PU2) at p = 3, 4.6e-10 and 2.3e-10 at p = 4, 3.1e-10 and 1.8e-10 at
#: p = 5, 1.6e-14 and 1.2e-14 at p = 6; the point ratio's distance from 2 is 2.1e-5, 3.9e-6, 3.9e-6, 1.8e-7 (1G) and
#: 4.6e-5, 2.8e-6, 2.8e-6, 1.3e-7 (PU2). p = 4 costs 6 s (1G) and 21 s (PU2) for both bodies.
_IMAGE_RESOLUTION = rung(4)
#: 20x the larger measured k difference at p = 4; 25x the larger measured ratio distance at p = 4.
_IMAGE_K_BAND = 1e-8
_IMAGE_RATIO_BAND = 1e-4


def _images(name: str) -> tuple[GeometrySpecification, GeometrySpecification]:
    """The slab [0, L] mirror | vacuum, and the doubled slab [0, 2L] vacuum | vacuum, of one mixture."""
    mixture = _IMAGE_MIXTURES[name]
    half = StructuredGeometry.slab((0.0, _HALF), (0,), left=ReflectiveBoundary(axis="x"), right=BC.vacuum)
    doubled = StructuredGeometry.slab((0.0, 2.0 * _HALF), (0,), left=BC.vacuum, right=BC.vacuum)
    return _eigen(half, mixture), _eigen(doubled, mixture)


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-slab-asym-method-of-images")
@pytest.mark.catches("ERR-034")
@pytest.mark.parametrize("name", sorted(_IMAGE_MIXTURES))
@pytest.mark.rests_on(_D5, _DOOR_K)
def test_a_mirror_is_the_symmetry_plane_of_the_doubled_vacuum_slab(name: str) -> None:
    """[B4c] k of the slab [0, 1] under a mirror at 0 and vacuum at 1 equals k of the vacuum slab [0, 2], to 1e-8.

    The fundamental mode of the symmetric vacuum slab is even about its
    centre, so its centre plane carries exactly the mirror's condition: the
    two problems have one eigenvalue (``peierls-greens-slab-asym-method-of-images``).
    The two specifications share no wall law and no breakpoint but 0, so the
    identity holds between two posings of the reference, not inside one; the
    band is the two posings' discretisation difference at rung 4 (module
    constants). It succeeds the old family's
    ``test_method_of_images_reflective_vacuum_equals_double_vacuum``.

    First red, ERR-034 re-dropped (``[M]`` 2026-10-10): the position along
    a slab line advanced by the arc length instead of the arc length times
    the cosine, the attenuation left honest. Three arms. ``err034``, the
    defect as catalogued (x - s), carries positions out of the body, and
    every row reds by a singular loss form (a red for a structural reason,
    ``vv`` anti-#18). ``err034-inclass`` (x - mu |mu| s) reds 2 of 4 by value
    (the PU2 rows: k 4.68 against 13.96) and 2 of 4 by ``NoFundamentalMode``
    (the 1G rows: a negative dominant eigenvalue). The value witness is
    qa's ``err034-scaled`` (the position advanced by 0.999 mu s): 4 of 4 red
    by value (PU2 k 0.66648 against 0.66662, 1G 0.058650 against 0.058669).
    Unlike the old family's collocation, the Galerkin assembly transports
    each basis function, so a flat emission does not hide the defect: the
    closed slabs of D5 red under the in-class arm (6 of 6 slab rows).
    """
    half, doubled = _images(name)
    k_half = _read(half, Eigenvalue(), _IMAGE_RESOLUTION)
    k_doubled = _read(doubled, Eigenvalue(), _IMAGE_RESOLUTION)
    assert abs(k_half - k_doubled) / k_doubled < _IMAGE_K_BAND, (k_half, k_doubled)


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-slab-asym-method-of-images")
@pytest.mark.catches("ERR-034")
@pytest.mark.parametrize("name", sorted(_IMAGE_MIXTURES))
@pytest.mark.rests_on(_HERE + "test_a_mirror_is_the_symmetry_plane_of_the_doubled_vacuum_slab")
def test_the_mirrored_slabs_flux_is_twice_the_doubled_slabs_right_half(name: str) -> None:
    """[B4c, the flux] phi_half(x) = 2 phi_doubled(1 + x) at x in {0, 0.3, 0.7, 1}, every group, to 1e-4 relative.

    Both fluxes are scaled to a total production of 1 over their body, and the
    doubled body produces twice over its two symmetric halves, so the factor
    is exactly 2: the row needs no normalisation of its own (the old row
    normalised both by their maxima, which is blind to a uniform scale). The
    point readings transport the converged emission, so the two sides share
    no node. It succeeds the old family's
    ``test_method_of_images_flux_shape_matches_doubled_slab_right_half``.

    First red (``[M]`` 2026-10-10): ERR-034 (arm ``err034``).
    """
    half, doubled = _images(name)
    groups = _IMAGE_MIXTURES[name].SigT.size
    for x in (0.0, 0.3 * _HALF, 0.7 * _HALF, _HALF):
        for g in range(groups):
            ratio = _read(half, PointValue(x, g), _IMAGE_RESOLUTION) / _read(doubled, PointValue(_HALF + x, g), _IMAGE_RESOLUTION)
            assert abs(ratio / 2.0 - 1.0) < _IMAGE_RATIO_BAND, (x, g, ratio)


# ── D8, the ordering theorem ─────────────────────────────────────────────

#: PU-2-0's two-group mixture (downscatter, chi in both groups), k_inf 2.684: every D8 row is two-group (``vv``,
#: 1-group degeneracy), and every body is homogeneous so that its closed form is k_inf of the mixture.
_ORDER_MIXTURE = _PU2
_ALBEDOS = (0.0, 0.3, 0.7, 1.0)
#: The values are read at rung 3 (the door default); the error estimate of each is its step from rung 2, an
#: over-estimate on a geometrically converging ladder. ``[M]`` 2026-10-10 (``probes/probe_d8.log``,
#: ``probe_d8b.log``): steps at most 4.3e-6 relative (sphere, slab, hollow sphere), gaps between neighbouring
#: albedos at least 1.3e-2 relative.
_ORDER_VALUE, _ORDER_STEP = rung(3), rung(2)
#: The margin over the step that a strict inequality must clear (spec D8: 10 x the D9 step).
_ORDER_MARGIN = 10.0


def _law(albedo: float):
    from orpheus.geometry.boundary import AlbedoBoundary, SpecularReturn, VacuumInflow

    return VacuumInflow() if albedo == 0.0 else AlbedoBoundary(albedo, SpecularReturn(axis="x"))


def _body(kind: str, inner: float, outer: float) -> StructuredGeometry:
    """The homogeneous bodies of D8: solid sphere and cylinder of radius 1; slab [0, 1]; hollow sphere and annulus
    (0.4, 1.4). ``inner`` is the slab's left wall; it is not read on a solid body."""
    match kind:
        case "sphere":
            return StructuredGeometry.sphere((0.0, 1.0), (0,), outer=_law(outer))
        case "cylinder":
            return StructuredGeometry.cylinder((0.0, 1.0), (0,), outer=_law(outer))
        case "slab":
            return StructuredGeometry.slab((0.0, 1.0), (0,), left=_law(inner), right=_law(outer))
        case "hollow":
            return StructuredGeometry.sphere((0.4, 1.4), (0,), inner=_law(inner), outer=_law(outer))
        case "annulus":
            return StructuredGeometry.cylinder((0.4, 1.4), (0,), inner=_law(inner), outer=_law(outer))
    raise AssertionError(kind)


def _k_inf(mixture) -> float:
    """rho((diag Sigma_t - Sigma_s^T)^-1 chi nu Sigma_f), written here from the mixture's arrays."""
    loss = np.diag(np.asarray(mixture.SigT, dtype=float)) - mixture.SigS[0].toarray().T
    production = np.outer(mixture.chi, mixture.SigP)
    return float(np.max(np.abs(np.linalg.eigvals(np.linalg.solve(loss, production)))))


def _four_group():
    """Fuel A at four groups, P0 (``_aba_reference.isotropic_mixture``): the old family's 4G rows' mixture."""
    from tests.gates.sn.verification.analytical._aba_reference import isotropic_mixture

    return isotropic_mixture("A", "4g")


def _mixture_of(kind: str):
    """``kind`` is a body, or a body and ``:A4`` for fuel A at four groups; PU2 otherwise."""
    return _four_group() if kind.endswith(":A4") else _ORDER_MIXTURE


def _k(kind: str, inner: float, outer: float, resolution: Resolution) -> float:
    return _read(_eigen(_body(kind.split(":")[0], inner, outer), _mixture_of(kind)), Eigenvalue(), resolution)


#: (id, body, the albedo path: "inner" varies the inner (left) wall with the outer at 0, "outer" the reverse,
#: "both" sets both walls to the albedo). A path whose last point closes every wall ends at k_inf.
_PATHS = [
    pytest.param("sphere", "outer", id="sphere-outer"),
    pytest.param("sphere:A4", "outer", id="sphere-outer-4g"),
    pytest.param("slab", "inner", id="slab-left"),
    pytest.param("slab", "outer", id="slab-right"),
    pytest.param("slab", "both", id="slab-both"),
    pytest.param("hollow", "inner", id="hollow-inner"),
    pytest.param("hollow", "outer", id="hollow-outer"),
    pytest.param("hollow", "both", id="hollow-both"),
    pytest.param("cylinder", "outer", id="cylinder-outer", marks=pytest.mark.slow),
    pytest.param("annulus", "inner", id="annulus-inner", marks=pytest.mark.slow),
    pytest.param("annulus", "outer", id="annulus-outer", marks=pytest.mark.slow),
]

#: The cylinder and the annulus cost about 60 s and 190 s a solve at rung 3 (``[M]`` 2026-10-10,
#: ``probes/probe_d8b.log``; 4 s and 11 s at rung 2), so their values are read at rung 2 and each value's error
#: estimate is its step to rung 3 sampled at the path's two ends (the vacuum body, shared by the annulus's two paths,
#: and the varied wall at albedo 1); elsewhere every value carries its own step. The annulus's both-walls path is not
#: posed: its closing end is the annulus of D5.
_ORDER_RUNGS = {"cylinder": (rung(2), rung(3)), "annulus": (rung(2), rung(3))}


def _walls_on(path: str, albedo: float) -> tuple[float, float]:
    return {"inner": (albedo, 0.0), "outer": (0.0, albedo), "both": (albedo, albedo)}[path]


@pytest.mark.l1
@pytest.mark.parametrize(("kind", "path"), _PATHS)
@pytest.mark.rests_on(_D5, _DOOR_K)
def test_k_increases_strictly_with_each_walls_albedo_from_the_vacuum_body_toward_k_inf(kind: str, path: str) -> None:
    """[D8] Along each albedo path 0 < 0.3 < 0.7 < 1, k increases strictly, by more than 10x its rung step, and
    stays below k_inf unless the path closes every wall, where it reads k_inf to 1e-10 or to 10x its step if larger.

    The theorem: the transport operator is positive and its boundary
    resolvent is monotone in each albedo (a returned neutron can only add
    to the flux), so by Perron-Frobenius the dominant eigenvalue is strictly
    increasing in each albedo, from the vacuum body to the closed one, whose
    k is the medium's k_inf (D5). Each value is read at rung 3 and its error
    estimated by its step from rung 2 (module constants); the strict
    inequalities hold with ``_ORDER_MARGIN`` over the largest step.
    It succeeds the old family's per-geometry Q rows (the vacuum body leaks,
    the mirror|vacuum and vacuum|mirror shells, the symmetric intermediate
    albedo, k_vacuum < k_closed, the 2G and 4G vacuum bodies below k_inf),
    which asserted bands or single inequalities with no margin. Every path is
    PU2 except ``sphere-outer-4g``, fuel A at four groups (the old 4G rows).

    Declared blind: the albedo pairing (which wall's albedo multiplies which
    traversal) keeps monotonicity in each albedo; B1, T11 and the
    unfolded-path rows catch it. First reds (``[M]`` 2026-10-10, battery
    ``scratch/characteristic_architecture/p1_step_e/ta_e1b/battery``): an albedo
    entering the closure as 1 - alpha (arm ``albedo-complement``).
    """
    k_inf = _k_inf(_mixture_of(kind))
    value_rung, step_rung = _ORDER_RUNGS.get(kind, (_ORDER_VALUE, _ORDER_STEP))
    sampled = _ALBEDOS if kind not in _ORDER_RUNGS else (_ALBEDOS[0], _ALBEDOS[-1])
    values, steps = [], []
    for albedo in _ALBEDOS:
        inner, outer = _walls_on(path, albedo)
        value = _k(kind, inner, outer, value_rung)
        values.append(value)
        if albedo in sampled:
            steps.append(abs(value - _k(kind, inner, outer, step_rung)) / value)
    step = max(steps)
    for low, high, a in zip(values, values[1:], _ALBEDOS[1:], strict=False):
        assert (high - low) / high > _ORDER_MARGIN * step, (kind, path, a, values, step)
    closes = kind.split(":")[0] in ("sphere", "cylinder") or path == "both"
    if closes:
        # the closed body reads k_inf to its own error estimate (``[M]``: the cylinder at rung 2 misses by 9.7e-8)
        assert abs(values[-1] - k_inf) / k_inf < max(1e-10, _ORDER_MARGIN * step), (values[-1], k_inf, step)
    else:
        assert (k_inf - values[-1]) / k_inf > _ORDER_MARGIN * step, (values[-1], k_inf, step)


@pytest.mark.l1
@pytest.mark.parametrize("albedo", [0.3, 0.7, 1.0])
@pytest.mark.rests_on(_HERE + "test_k_increases_strictly_with_each_walls_albedo_from_the_vacuum_body_toward_k_inf")
def test_a_homogeneous_slabs_k_is_blind_to_which_wall_carries_the_albedo(albedo: float) -> None:
    """[D8, the reflection] A homogeneous slab and its mirror image are one problem: k(alpha | 0) = k(0 | alpha) to 1e-12.

    The reflection x -> L - x maps the body onto itself and swaps its walls.
    ``[M]`` 2026-10-10 (``probes/probe_d8.log``): 8.5e-16 at most at rungs 3
    and 4. It is the declared blindness of the ordering row to the albedo
    pairing, made a row: on a homogeneous slab no k-valued row can see a swap
    of the two walls. First red: a per-wall convention drift, the left wall
    reading its partial mirror's albedo squared (arm ``left-wall-albedo-squared``).
    """
    left = _k("slab", albedo, 0.0, _ORDER_VALUE)
    right = _k("slab", 0.0, albedo, _ORDER_VALUE)
    assert abs(left - right) / right < 1e-12, (left, right)


#: The thickening ladder of the vacuum bodies, in cm (PU2: Sigma_t 0.1 to 0.3 per cm in its two groups, so 0.5 cm
#: is thin and 32 cm thick). ``[M]`` 2026-10-10 (``probes/probe_d8b.log``): sphere k 0.106, 0.212, 0.416, 0.794,
#: 1.386 at 0.5 to 8 cm against k_inf 2.684, each at most 0.1 s at rung 3.
_THICKENING = (0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0)


@pytest.mark.l1
@pytest.mark.parametrize("kind", ["sphere", "slab"])
@pytest.mark.rests_on(_HERE + "test_k_increases_strictly_with_each_walls_albedo_from_the_vacuum_body_toward_k_inf")
def test_a_vacuum_bodys_k_rises_toward_k_inf_as_it_thickens(kind: str) -> None:
    """[D8, the thickening] On the vacuum sphere and slab of PU2, k increases strictly with the size and its distance
    from k_inf decreases strictly, each step clearing 10x the rung step; every k stays below k_inf.

    The leakage fraction of a vacuum body falls monotonically as it grows,
    so k rises toward the infinite medium's k_inf. A finite ladder cannot
    prove the limit; the row asserts the monotone approach over sizes from
    thin to 32 cm (the old rows: ``test_a1_vacuum_thick_sphere_approaches_k_inf``,
    ``test_a2_thick_sphere_approaches_k_inf``, which banded one thick sphere
    at 0.95 k_inf). First red (``[M]`` 2026-10-10): a vacuum wall read as a
    mirror (arm ``vacuum-as-mirror``), every k then k_inf. Declared blind,
    measured green: a vacuum wall returning half (``vacuum-half-mirror``),
    the scattering transposed, Sigma_t read 1e-9 large; each keeps the
    approach monotone. A monotone law is a weak instrument; the value rows
    (D5, D2', PS-1982) carry the size.
    """
    k_inf = _k_inf(_ORDER_MIXTURE)
    values, steps = [], []
    for size in _THICKENING:
        geometry = (StructuredGeometry.sphere((0.0, size), (0,), outer=_law(0.0)) if kind == "sphere"
                    else StructuredGeometry.slab((0.0, size), (0,), left=_law(0.0), right=_law(0.0)))
        spec = _eigen(geometry, _ORDER_MIXTURE)
        value = _read(spec, Eigenvalue(), _ORDER_VALUE)
        values.append(value)
        steps.append(abs(value - _read(spec, Eigenvalue(), _ORDER_STEP)) / value)
    step = max(steps)
    gaps = [(k_inf - v) / k_inf for v in values]
    assert min(gaps) > _ORDER_MARGIN * step, (gaps, step)
    for low, high in zip(gaps, gaps[1:], strict=False):
        assert low - high > _ORDER_MARGIN * step, (gaps, step)


# ── D10, the fundamental mode in the positive cone ───────────────────────

#: The vacuum bodies whose mode D10 reads: a sphere of radius 2 and a slab [0, 2] of PU2, two groups.
_MODE_BODIES = {
    "sphere": StructuredGeometry.sphere((0.0, 2.0), (0,), outer=BC.vacuum),
    "slab": StructuredGeometry.slab((0.0, 2.0), (0,), left=BC.vacuum, right=BC.vacuum),
    "cylinder": StructuredGeometry.cylinder((0.0, 2.0), (0,), outer=BC.vacuum),
}
_MODE_POINTS = (0.0, 0.25, 0.5, 0.75, 1.0)   # fractions of the half-width, centre first


@pytest.mark.l1
@pytest.mark.parametrize("kind", ["sphere", "slab", pytest.param("cylinder", marks=pytest.mark.slow)])
@pytest.mark.rests_on(_HERE + "test_k_increases_strictly_with_each_walls_albedo_from_the_vacuum_body_toward_k_inf")
def test_a_vacuum_bodys_mode_is_positive_and_falls_from_its_centre_to_its_wall(kind: str) -> None:
    """[D10] The fundamental mode of a homogeneous vacuum body, read at points: positive in every group, strictly
    decreasing from the centre to the wall, and on the slab even about its centre to 1e-10 relative.

    The theorem: the dominant mode of a positive operator is in the positive
    cone (Perron-Frobenius), and on a homogeneous body symmetric about its
    centre it is symmetric and peaked there (the leakage is largest at the
    wall). It succeeds the old family's shape rows
    (``test_a1_vacuum_spatial_mode_nontrivial``, ``test_mr_issue132_spatial_mode_physical``,
    ``test_a2_phi_shape_qualitative_agreement``'s monotone leg), which banded a
    spread and a centre-to-surface ratio; the PS-1982 shape itself is in
    ``test_characteristic_independent_references.py``. First red
    (``[M]`` 2026-10-10): the eigen answer's emission returned with its sign
    flipped (arm ``mode-sign-flipped``).
    """
    geometry = _MODE_BODIES[kind]
    spec = _eigen(geometry, _ORDER_MIXTURE)
    lo, hi = geometry.breakpoints[0], geometry.breakpoints[-1]
    centre, half = (0.0, hi) if kind != "slab" else (0.5 * (lo + hi), 0.5 * (hi - lo))
    for g in range(_ORDER_MIXTURE.SigT.size):
        profile = [_read(spec, PointValue(centre + f * half, g), _ORDER_VALUE) for f in _MODE_POINTS]
        assert min(profile) > 0.0, (g, profile)
        assert np.all(np.diff(profile) < 0.0), (g, profile)
        if kind == "slab":
            mirrored = [_read(spec, PointValue(centre - f * half, g), _ORDER_VALUE) for f in _MODE_POINTS]
            np.testing.assert_allclose(mirrored, profile, rtol=1e-10, atol=0.0)


# ── the hollow body as its cavity closes ─────────────────────────────────

#: The old rows' one-group medium on a body of radius 5 cm, vacuum at every wall.
_CAVITY_MIXTURE = _IMAGE_MIXTURES["1g"]
_CAVITY_FRACTIONS = (1e-1, 1e-2, 1e-3)
#: How much the distance to the solid body must shrink per decade of the cavity's radius: a vacuum cavity absorbs
#: what enters it, a fraction of the lines that scales as its cross section, (R_in / R)^2 on the sphere and R_in / R
#: on the cylinder. ``[M]`` 2026-10-10 (``probes/probe_rin.py``): sphere 1.9e-2, 2.2e-4, 1.8e-6 (factors 88 and
#: 122 per decade); cylinder at rung 2 9.4e-2, 1.1e-2, 1.2e-3 (factors 8.4 and 9.6). The bound is half the order's
#: factor.
_CAVITY_FACTOR = {"sphere": 50.0, "cylinder": 5.0}
_CAVITY_RUNG = {"sphere": rung(3), "cylinder": rung(2)}


@pytest.mark.l1
@pytest.mark.parametrize("kind", ["sphere", pytest.param("cylinder", marks=pytest.mark.slow)])
@pytest.mark.rests_on(_D5, _DOOR_K)
def test_a_vacuum_shells_k_tends_to_the_solid_bodys_as_its_cavity_closes(kind: str) -> None:
    """[X, the edge R_in -> 0] The hollow vacuum body's k tends to the solid vacuum body's as R_in / R falls from
    1e-1 to 1e-3, the distance shrinking per decade by at least half the factor its order predicts (100 on the
    sphere, 10 on the cylinder) and by at most twice it; and k(hollow) < k(solid) at every R_in.

    The edge of the hollow body equal to a foundation, the solid body: the
    impact parameter partitions the lines at b = R_in (the partition itself
    is verified by ``test_characteristic_closure.py::test_the_period_matches_the_hand_counted_table``,
    which carries the two ``*-impact-parameter-partition`` labels; this row
    asserts a limit, not the partition), and as R_in -> 0
    the lines through the cavity, rank 2, vanish as its cross section. It
    succeeds the old family's ``test_R_in_to_zero_limit_matches_solid_sphere_vacuum``
    and ``test_R_in_to_zero_limit_matches_solid_cylinder_vacuum``, which
    asserted 1e-7 and 1e-4 at R_in = 1e-3 R; the new reference reads 1.8e-6 on
    the sphere there, the second-order distance a vacuum cavity of that size
    absorbs, so the old row's 1e-7 was below the physics and is not inherited
    (a finding in the README). First red (``[M]`` 2026-10-10): the inner
    vacuum wall read as a mirror (arm ``inner-vacuum-as-mirror``), which makes
    the cavity transparent: at R_in = 1e-3 R the hollow sphere's k then
    exceeds the solid one's by 4.2e-7 relative (the material removed absorbed
    more than it produced), and the k(hollow) < k(solid) leg reds.
    """
    make = StructuredGeometry.sphere if kind == "sphere" else StructuredGeometry.cylinder
    resolution = _CAVITY_RUNG[kind]
    solid = _read(_eigen(make((0.0, 5.0), (0,), outer=BC.vacuum), _CAVITY_MIXTURE), Eigenvalue(), resolution)
    distances = []
    for fraction in _CAVITY_FRACTIONS:
        hollow = make((5.0 * fraction, 5.0), (0,), inner=BC.vacuum, outer=BC.vacuum)
        k = _read(_eigen(hollow, _CAVITY_MIXTURE), Eigenvalue(), resolution)
        assert k < solid, (fraction, k, solid)
        distances.append((solid - k) / solid)
    for wide, narrow in zip(distances, distances[1:], strict=False):
        assert narrow * _CAVITY_FACTOR[kind] < wide < narrow * 4.0 * _CAVITY_FACTOR[kind], (kind, distances)


#: A closed three-region one-group body with a 4x total cross-section contrast at each interface (the old row's
#: 10x contrast fixture's spirit: Sigma_t 2.0 | 0.5 | 1.5, fissile outer regions), radii (1, 2.5, 4).
_CONTRAST_MATERIALS = {0: _mixture([2.0], [[1.5]], None, [0.6], [1.0]), 1: _mixture([0.5], [[0.45]]),
                       2: _mixture([1.5], [[1.0]], None, [0.6], [1.0])}
_CONTRAST_BREAKPOINTS = (0.0, 1.0, 2.5, 4.0)


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-cylinder-mr-interface-continuity")
@pytest.mark.parametrize("kind", ["sphere", pytest.param("cylinder", marks=pytest.mark.slow)])
@pytest.mark.rests_on(_HERE + "test_a_vacuum_bodys_mode_is_positive_and_falls_from_its_centre_to_its_wall")
def test_the_modes_scalar_flux_is_continuous_across_a_material_interface(kind: str) -> None:
    """[D10, interface continuity] The fundamental mode's scalar flux read at r_k (1 - eps) and r_k (1 + eps) on
    each interface of a closed contrasted body: the jump falls by more than 30x as eps falls from 1e-3 to 1e-5
    and again to 1e-7, and is below 1e-5 relative at 1e-7.

    The angular flux has a kink at an interface (the source jumps with
    Sigma_t) but the angle-integrated flux is continuous: a theorem of the
    integral equation (``peierls-greens-cylinder-mr-interface-continuity``).
    The point reading transports the emission to each point, so the two
    one-sided readings share nothing but the emission. ``[M]`` 2026-10-10
    (``probes/probe_d10b.py``, sphere, rung 3): jumps 6.2e-4, 1.1e-5, 1.6e-7
    (r = 1) and 1.2e-3, 1.5e-5, 1.9e-7 (r = 2.5). It succeeds the old family's
    ``test_mr_interface_continuity_3region`` (a 5e-2 gross-indexing band on
    spline limits). First red (``[M]`` 2026-10-10): the point value read from
    the Galerkin flux coefficients instead of transported (arm
    ``point-from-galerkin``; the projection is discontinuous at panel ends).
    Declared blind: any error in the transport itself (a wrong Sigma_t, a
    wrong closure) keeps the transported reading continuous; the value rows
    carry those.
    """
    make = StructuredGeometry.sphere if kind == "sphere" else StructuredGeometry.cylinder
    spec = GeometrySpecification(Materials(_CONTRAST_MATERIALS), make(_CONTRAST_BREAKPOINTS, (0, 1, 2), outer=BC.reflective),
                                 Eigen(_K))
    resolution = rung(3) if kind == "sphere" else rung(2)
    for r_k in _CONTRAST_BREAKPOINTS[1:-1]:
        jumps = []
        for eps in (1e-3, 1e-5, 1e-7):
            inside = _read(spec, PointValue(r_k * (1.0 - eps), 0), resolution)
            outside = _read(spec, PointValue(r_k * (1.0 + eps), 0), resolution)
            jumps.append(abs(inside - outside) / abs(outside))
        assert jumps[1] * 30.0 < jumps[0] and jumps[2] * 30.0 < jumps[1] and jumps[2] < 1e-5, (r_k, jumps)


@pytest.mark.l1
@pytest.mark.rests_on(_D5, _HERE + "test_a_vacuum_bodys_mode_is_positive_and_falls_from_its_centre_to_its_wall")
def test_a_closed_one_group_heterogeneous_sphere_balances_and_stays_below_the_fuels_k_inf() -> None:
    """[D6/D5, the closed heterogeneous body] Fuel A | moderator B at one group in a mirrored sphere (0, 0.5, 1):
    the fundamental is certified (the door answers: real, simple, single-signed, ``DensePencil``'s refusals); its
    k equals the body's production over its absorption, both read as flux integrals, to 1e-12; and
    0 < k_vacuum < k < k_inf(A).

    In one group a closed body's k is the ratio of its production to its
    absorption, so it lies below the largest region's nu Sigma_f / Sigma_a
    (here fuel A's 1.5; the moderator produces nothing). It succeeds the old
    family's issue-132 rows (``test_mr_issue132_no_catastrophe_closed_sphere``:
    0.5 < k < 0.95; ``test_mr_issue132_vacuum_below_closed``). ``[M]``
    2026-10-10: k 0.73679519, the balance 2.0e-15 relative. First red
    (``[M]`` 2026-10-10): Sigma_t read 1e-9 large (arm ``total-scaled``: the
    balance reads the true Sigma_a and moves by 1e-8 relative).

    Declared blind: the balance holds for ANY flux the pencil returns, so an
    error that only redistributes neutrons in space or between the regions
    (a wrong transport kernel, a wrong closure, a transposed region pairing
    that conserves) keeps production / absorption = k and leaves this row
    green; the bound k < k_inf(A) is equally loose. The shape and value rows
    (D10, PS-1982, D2') carry the redistribution.
    """
    from tests.gates.sn.verification.analytical._aba_reference import isotropic_mixture

    from orpheus.numerics.mesh_free_function import RegionwiseConstant
    from orpheus.numerics.observable import FluxIntegral

    fuel, moderator = isotropic_mixture("A", "1g"), isotropic_mixture("B", "1g")
    absorption = [float(m.SigT[0] - m.SigS[0].toarray().sum()) for m in (fuel, moderator)]

    def body(outer: BC) -> GeometrySpecification:
        return GeometrySpecification(Materials({0: fuel, 1: moderator}),
                                     StructuredGeometry.sphere((0.0, 0.5, 1.0), (0, 1), outer=outer), Eigen(_K))

    closed = body(BC.reflective)
    k = _read(closed, Eigenvalue(), rung(3))
    production = _read(closed, FluxIntegral(RegionwiseConstant(np.array([[fuel.SigP[0]], [moderator.SigP[0]]]))), rung(3))
    absorbed = _read(closed, FluxIntegral(RegionwiseConstant(np.array([[absorption[0]], [absorption[1]]]))), rung(3))
    assert abs(production / absorbed / k - 1.0) < 1e-12, (k, production / absorbed)
    k_vacuum = _read(body(BC.vacuum), Eigenvalue(), rung(3))
    assert 0.0 < k_vacuum < k < fuel.SigP[0] / absorption[0], (k_vacuum, k)



#: The fast tier's cylinder and annulus (qa's ruling of 2026-10-10: a retirement must not move coverage from the
#: default tier to slow-only). Rung 2 of the joint ladder, about 4 s a cylinder solve and 10 s an annulus one
#: (``[M]`` ``probes/probe_d8b.log``); the closed bodies read k_inf there to 9.7e-8 (``[M]`` the slow D8 run).
_FAST_RUNG = rung(2)
_FAST_CLOSED_BAND = 1e-6


@pytest.mark.l1
@pytest.mark.parametrize(("kind", "path"), [("cylinder", "outer"), ("annulus", "inner")], ids=["cylinder-outer", "annulus-inner"])
@pytest.mark.rests_on(_D5, _DOOR_K)
def test_a_cylinders_k_increases_from_its_vacuum_body_to_k_inf_in_the_fast_tier(kind: str, path: str) -> None:
    """[D8, the fast tier] The cylinder (outer wall) and the annulus (inner wall, outer vacuum) of PU2 at rung 2:
    k(0) < k(0.5) < k(1), each gap above 1e-3 relative, the cylinder closed at albedo 1 reading k_inf to 1e-6, the
    annulus at (1, 0) below k_inf by more than 1e-3 and the annulus closed at (1, 1) reading k_inf to 1e-6.

    The default-tier successor of the old family's fast cylinder and annulus
    rows (the vacuum cylinder leaks, the annulus mirror|vacuum and
    vacuum|mirror shells, their closed bodies); the full ordering with each
    value's own step is the slow ``test_k_increases_strictly_with_each_walls_albedo_from_the_vacuum_body_toward_k_inf``.
    ``[M]`` 2026-10-10: the gaps are 0.1 relative or more; the closed
    cylinder's distance to k_inf at rung 2 is 9.7e-8. First reds by value
    (``[M]`` 2026-10-10): the in-plane speed read 1e-4 large in the chord (qa's
    ``obliquity-dropped`` at QSPEED 1.0001, battery arm ``cylinder-speed-scaled``):
    the closed cylinder reads 2.683010 against k_inf 2.683767, the cylinder row
    red, the annulus row green (its legs are inequalities); the albedo read as
    its complement (``albedo-complement``): both rows.
    """
    k_inf = _k_inf(_ORDER_MIXTURE)
    values = [_k(kind, *_walls_on(path, a), _FAST_RUNG) for a in (0.0, 0.5, 1.0)]
    for low, high in zip(values, values[1:], strict=False):
        assert (high - low) / high > 1e-3, (kind, values)
    if kind == "cylinder":
        assert abs(values[-1] - k_inf) / k_inf < _FAST_CLOSED_BAND, (values[-1], k_inf)
    else:
        assert (k_inf - values[-1]) / k_inf > 1e-3, (values[-1], k_inf)
        closed = _k(kind, 1.0, 1.0, _FAST_RUNG)
        assert abs(closed - k_inf) / k_inf < _FAST_CLOSED_BAND, (closed, k_inf)
