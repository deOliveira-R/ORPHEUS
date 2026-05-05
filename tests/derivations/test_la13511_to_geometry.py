r"""Foundation tests for :meth:`La13511Case.to_geometry`.

The Phase B adapter on
:mod:`orpheus.derivations.continuous.sood_registry.la13511` lifts the
published critical dimension off ``GeometrySpec`` (a registry-truth
artefact) onto :class:`La13511Truth` and provides
:meth:`La13511Case.to_geometry` to build a new-API
:class:`~orpheus.geometry.structured_geometry.StructuredGeometry`
on demand.

These tests cover:

* Slab (1G): full-width cm + vacuum-vacuum BCs, geometry tag
  ``"SLB"`` with two endpoints.
* Sphere (1G): radius-cm region + single vacuum BC, tag ``"SPH"``.
* Cylinder (1G): radius-cm region + single vacuum BC, tag ``"CYL"``.
* Infinite-medium: raises :class:`ValueError` pointing the caller at
  ``MomentSpace.solve_kinf`` / ``solve_homogeneous_infinite``.
* Round-trip consistency: the new-API ``domain_extent_cm`` matches
  the legacy ``geometry_spec.domain_extent_cm`` to machine precision
  modulo author transcription rounding.
* Provenance: the structured ``provenance`` field is populated on
  every case in the registry.
"""
from __future__ import annotations

import math

import pytest

from orpheus.derivations.continuous.sood_registry import (
    LA13511_CASES,
    PUA_1_0_IN,
    PUA_1_0_SL,
    PUB_1_0_CY_STUB,
    PU_2_0_IN,
    Provenance,
    UA_1_0_CY_STUB,
    UA_1_0_SL_STUB,
    UA_1_0_SP_STUB,
)
from orpheus.derivations.continuous.sood_registry.atalay1997 import (
    ATALAY_ALL_CASES,
)
from orpheus.derivations.continuous.sood_registry.la13511 import _ALL_CASES
from orpheus.geometry.mesh import BC
from orpheus.geometry.structured_geometry import (
    Region as StructuredRegion,
    StructuredGeometry,
)


# ═══════════════════════════════════════════════════════════════════
# Slab adapter — full-width cm, vacuum-vacuum BCs
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_slab_to_geometry_full_width_with_vacuum_bcs() -> None:
    r"""1G slab case: ``to_geometry()`` returns SLB with full-width cm.

    For ``Ua-1-0-SL`` (U-235 (a), 1G isotropic, ``c=1.30``):

    * ``truth.critical_dimension_mfp = 0.93772556`` (half-thickness in mfp).
    * ``Σ_t = 0.32640``.
    * Expected full slab width: :math:`2 \cdot 0.93772556 / 0.32640
      = 5.74586815...` cm.

    The new-API geometry must be ``"SLB"`` with one region of that
    full-width cm thickness and vacuum-vacuum BCs.
    """
    case = UA_1_0_SL_STUB
    geom = case.to_geometry()

    sigma_t = float(case.materials[0].SigT[0])
    expected_full_width_cm = 2.0 * case.truth.critical_dimension_mfp / sigma_t

    assert isinstance(geom, StructuredGeometry)
    assert geom.geometry == "SLB"
    assert len(geom.regions) == 1
    region = geom.regions[0]
    assert isinstance(region, StructuredRegion)
    assert region.mat_id == case.geometry_spec.mat_id
    assert math.isclose(
        region.outer_thickness_cm, expected_full_width_cm, rel_tol=1e-12,
    )
    assert math.isclose(
        geom.domain_extent_cm, expected_full_width_cm, rel_tol=1e-12,
    )

    # SLB carries 2 BCs (left, right); both vacuum.
    assert geom.n_endpoints == 2
    assert geom.bcs == (BC.vacuum, BC.vacuum)


# ═══════════════════════════════════════════════════════════════════
# Sphere adapter — radius cm, single vacuum BC
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_sphere_to_geometry_radius_with_vacuum_outer_bc() -> None:
    r"""1G sphere case: ``to_geometry()`` returns SPH with radius cm.

    For ``Ua-1-0-SP`` (U-235 (a), 1G isotropic, ``c=1.30``):

    * ``truth.critical_dimension_mfp = 2.4248249802`` (radius in mfp).
    * ``Σ_t = 0.32640``.
    * Expected radius: :math:`2.4248249802 / 0.32640 \approx 7.4290`
      cm.
    """
    case = UA_1_0_SP_STUB
    geom = case.to_geometry()

    sigma_t = float(case.materials[0].SigT[0])
    expected_R_cm = case.truth.critical_dimension_mfp / sigma_t

    assert isinstance(geom, StructuredGeometry)
    assert geom.geometry == "SPH"
    assert len(geom.regions) == 1
    region = geom.regions[0]
    assert isinstance(region, StructuredRegion)
    assert region.mat_id == case.geometry_spec.mat_id
    assert math.isclose(
        region.outer_thickness_cm, expected_R_cm, rel_tol=1e-12,
    )
    assert math.isclose(geom.domain_extent_cm, expected_R_cm, rel_tol=1e-12)

    # SPH carries 1 BC (outer); centreline reflective is implicit.
    assert geom.n_endpoints == 1
    assert geom.bcs == (BC.vacuum,)


# ═══════════════════════════════════════════════════════════════════
# Cylinder adapter — radius cm, single vacuum BC
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_cylinder_to_geometry_radius_with_vacuum_outer_bc() -> None:
    r"""1G cylinder case: ``to_geometry()`` returns CYL with radius cm.

    For ``Ua-1-0-CY`` (U-235 (a), 1G isotropic, ``c=1.30``):

    * ``truth.critical_dimension_mfp = 1.72500292`` (radius in mfp).
    * ``Σ_t = 0.32640``.
    * Expected radius: :math:`1.72500292 / 0.32640 \approx 5.2849` cm.
    """
    case = UA_1_0_CY_STUB
    geom = case.to_geometry()

    sigma_t = float(case.materials[0].SigT[0])
    expected_R_cm = case.truth.critical_dimension_mfp / sigma_t

    assert isinstance(geom, StructuredGeometry)
    assert geom.geometry == "CYL"
    assert len(geom.regions) == 1
    region = geom.regions[0]
    assert isinstance(region, StructuredRegion)
    assert region.mat_id == case.geometry_spec.mat_id
    assert math.isclose(
        region.outer_thickness_cm, expected_R_cm, rel_tol=1e-12,
    )
    assert math.isclose(geom.domain_extent_cm, expected_R_cm, rel_tol=1e-12)

    # CYL carries 1 BC (outer); centreline reflective is implicit.
    assert geom.n_endpoints == 1
    assert geom.bcs == (BC.vacuum,)


# ═══════════════════════════════════════════════════════════════════
# Infinite-medium — must raise with a directive error message
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_infinite_to_geometry_raises_with_kinf_pointer() -> None:
    """Infinite-medium cases raise ValueError pointing at the kinf APIs.

    ``PUa-1-0-IN`` is a bare infinite medium — there is no geometry
    to build. The error message must point the caller at
    ``MomentSpace.solve_kinf`` / ``solve_homogeneous_infinite``
    so a fresh consumer can recover.
    """
    case = PUA_1_0_IN
    with pytest.raises(ValueError) as excinfo:
        case.to_geometry()

    msg = str(excinfo.value)
    assert "infinite" in msg.lower()
    # Point at one of the kinf entry points by name.
    assert (
        "solve_kinf" in msg or "solve_homogeneous_infinite" in msg
    ), f"error message must mention a k_inf entry point; got: {msg}"


@pytest.mark.foundation
def test_infinite_2g_to_geometry_raises() -> None:
    """2G infinite case (``PU-2-0-IN``) also raises."""
    with pytest.raises(ValueError):
        PU_2_0_IN.to_geometry()


# ═══════════════════════════════════════════════════════════════════
# Round-trip — new-API extent matches legacy GeometrySpec extent
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.parametrize(
    "case_id",
    [
        c.case_id
        for c in _ALL_CASES
        if c.geometry_spec.geometry != "infinite"
    ],
)
def test_to_geometry_extent_matches_geometry_spec(case_id: str) -> None:
    r"""``geom.domain_extent_cm`` matches ``case.geometry_spec.domain_extent_cm``.

    Both quantities are derived from the same published mfp value and
    the same :math:`\Sigma_t`. Author transcription rounding (the
    catalogue stores published-mfp and published-cm separately, both
    transcribed verbatim from LA-13511) puts an upper bound on the
    achievable agreement at ~3e-7 relative — well below the 1e-5
    benchmark tolerance.
    """
    case = LA13511_CASES[case_id]
    geom = case.to_geometry()
    expected_cm = case.geometry_spec.domain_extent_cm
    actual_cm = geom.domain_extent_cm
    assert math.isclose(actual_cm, expected_cm, rel_tol=1e-5), (
        f"{case_id}: legacy extent={expected_cm}, new extent={actual_cm}"
    )


# ═══════════════════════════════════════════════════════════════════
# Provenance — every case carries structured citation metadata
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.parametrize("case_id", list(LA13511_CASES.keys()))
def test_la13511_case_has_provenance(case_id: str) -> None:
    """Every LA-13511 registry case has a populated Provenance."""
    case = LA13511_CASES[case_id]
    assert case.provenance is not None, (
        f"{case_id}: provenance must be populated post-Phase-B"
    )
    assert isinstance(case.provenance, Provenance)
    assert case.provenance.paper_id == "LA-13511", (
        f"{case_id}: paper_id should be 'LA-13511', "
        f"got {case.provenance.paper_id!r}"
    )
    # Mirror invariants vs the legacy flat fields (Phase F drops
    # the flat fields; until then they must agree).
    assert case.provenance.paper_table == case.sood_table
    assert case.provenance.primary_reference == case.primary_reference
    assert case.provenance.notes == case.notes


@pytest.mark.foundation
@pytest.mark.parametrize(
    "case_id",
    [c.case_id for c in ATALAY_ALL_CASES],
)
def test_atalay_case_has_provenance(case_id: str) -> None:
    """Every Atalay 1997 case has a Provenance with paper_id="Atalay-1997"."""
    case_by_id = {c.case_id: c for c in ATALAY_ALL_CASES}
    case = case_by_id[case_id]
    assert case.provenance is not None
    assert isinstance(case.provenance, Provenance)
    assert case.provenance.paper_id == "Atalay-1997"
    assert case.provenance.primary_reference == case.primary_reference
    assert case.provenance.notes == case.notes


# ═══════════════════════════════════════════════════════════════════
# Truth.critical_dimension_mfp — populated for every finite case
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.parametrize("case_id", list(LA13511_CASES.keys()))
def test_truth_critical_dimension_mfp_mirrors_geometry_spec(
    case_id: str,
) -> None:
    """``truth.critical_dimension_mfp`` matches the legacy
    ``geometry_spec.critical_dimension_mfp`` exactly.

    Both fields are populated during Phase B (the lift). They will
    diverge only if a post-Phase-F migration discovers a
    transcription error — until then they must agree by construction.
    """
    case = LA13511_CASES[case_id]
    assert case.truth.critical_dimension_mfp == (
        case.geometry_spec.critical_dimension_mfp
    ), (
        f"{case_id}: truth.critical_dimension_mfp"
        f"={case.truth.critical_dimension_mfp} disagrees with "
        f"geometry_spec.critical_dimension_mfp"
        f"={case.geometry_spec.critical_dimension_mfp}"
    )
