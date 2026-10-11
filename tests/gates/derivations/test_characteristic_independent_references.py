"""The characteristic reference against structurally independent references: PS-1982, Sood's and WM-72's cylinder, Nystrom.

P1 step (e1b) of the characteristic-reference campaign
(``.claude/plans/characteristic_reference_architecture.md``, "Step (e): the audit
and four rulings"): the old trajectory-resolvent family's cross-reference rows,
re-pointed at the characteristic reference before step (e2) deleted the family.
The row-by-row map is ``scratch/characteristic_architecture/p1_step_e/ta_e1b/README.md``.

Each reference here shares no transport primitive with the characteristic
reference above the trusted-library line (``vv-principles``, structural
independence):

- **PS-1982** (:func:`~orpheus.derivations.continuous.peierls_nystrom.ps1982_reference.solve_ps1982_vacuum_sphere`):
  the vacuum sphere's scalar-flux integral equation in its E_1 kernel, power
  iteration with QUADPACK on the log singularity; no lines, no closure;
- **Sood 2003 Table 10** (Ua-1-0-CY, a printed critical radius) and **WM-72**
  (:func:`~orpheus.derivations.continuous.singular_eigenfunction.solve_singular_eigenfunction_cylinder_bare_critical`,
  Case singular eigenfunctions and a Wiener-Hopf factorisation): no
  Bickley-Naylor function anywhere;
- the **Peierls Nystrom** slab and sphere: withdrawn under #506 with their
  generators; their re-pointed rows are ``test_characteristic_nystrom_withdrawn.py``
  (a file of withdrawn rows only, so the withdrawal's M5 run over the files that
  hold them stays cheap).

Every characteristic reading runs in this process (:func:`~orpheus.numerics.traced_memo.bypass`).
Bands are measured (``[M]`` 2026-10-10, probes in
``scratch/characteristic_architecture/p1_step_e/ta_e1b/probes/``); first reds are arms of
``scratch/characteristic_architecture/p1_step_e/ta_e1b/battery/``.
"""
from __future__ import annotations

from decimal import Decimal
from functools import lru_cache

import numpy as np
import pytest

from orpheus.data.cells import CellCoefficient, Channel
from orpheus.data.materials import Materials
from orpheus.derivations.continuous.characteristic import CharacteristicDerivation, Resolution
from orpheus.derivations.continuous.sood_registry import sood2003 as sood
from orpheus.geometry import BC, StructuredGeometry
from orpheus.numerics.observable import Eigenvalue, PointValue
from orpheus.numerics.question import Eigen
from orpheus.numerics.traced_memo import bypass
from orpheus.specification.specification import GeometrySpecification
from tests.gates.derivations._characteristic_ladders import rung
from tests.gates.derivations._ladder_rules import geometric_error
from tests.gates.derivations.test_characteristic_system import _mixture

pytestmark = pytest.mark.filterwarnings("error::RuntimeWarning")

_HERE = "tests/gates/derivations/test_characteristic_independent_references.py::"
_SYSTEM = "tests/gates/derivations/test_characteristic_system.py::"
_ALBEDOS = "tests/gates/derivations/test_characteristic_albedos.py::"
_D5 = _SYSTEM + "test_a_closed_body_reads_k_inf_with_a_flat_flux_in_the_infinite_mediums_group_ratio"
_D2 = _SYSTEM + "test_k_is_one_at_soods_two_group_critical_sizes"
_D8 = _ALBEDOS + "test_k_increases_strictly_with_each_walls_albedo_from_the_vacuum_body_toward_k_inf"

_K = CellCoefficient.every(Channel.FISSION_EMISSION)


@lru_cache(maxsize=None)
def _derivation(specification: GeometrySpecification, resolution: Resolution) -> CharacteristicDerivation:
    return CharacteristicDerivation(specification, resolution)


def _read(specification: GeometrySpecification, observable, resolution: Resolution) -> float:
    with bypass():
        return float(_derivation(specification, resolution).evaluate(observable).value)


def _one_group(sigma_t: float, sigma_s: float, nu_sigma_f: float):
    return _mixture([sigma_t], [[sigma_s]], None, [nu_sigma_f], [1.0])


def _vacuum_sphere(radius: float, mixture) -> GeometrySpecification:
    return GeometrySpecification(Materials({0: mixture}), StructuredGeometry.sphere((0.0, radius), (0,), outer=BC.vacuum),
                                 Eigen(_K))


# ── PS-1982, the vacuum sphere ───────────────────────────────────────────

#: The old rows' four vacuum spheres (R cm, Sigma_t, Sigma_s, nu Sigma_f): optical radii 2.5 and 5, scattering
#: ratios 0.4 and 0.6.
_PS_CASES = [
    pytest.param((5.0, 0.5, 0.20, 0.025), id="tauR2.5-c0.4"),
    pytest.param((10.0, 0.5, 0.20, 0.025), id="tauR5-c0.4", marks=pytest.mark.slow),
    pytest.param((5.0, 0.5, 0.30, 0.025), id="tauR2.5-c0.6", marks=pytest.mark.slow),
    pytest.param((10.0, 0.5, 0.30, 0.025), id="tauR5-c0.6", marks=pytest.mark.slow),
]
#: PS-1982 at 40 Gauss-Legendre nodes, converged to 1e-11 in its power iteration (59 to 110 iterations, 12 to 22 s).
_PS_NODES, _PS_TOL = 40, 1e-11
#: The characteristic reference at rung 4. ``[M]`` 2026-10-10 (``probes/probe_ps1982.log``): k misses PS-1982 at
#: 40 nodes by 2.8e-7 to 4.2e-7 relative at rung 4 (6.2e-7 to 8.8e-7 at rung 3); PS-1982's own step from 30 to 40
#: nodes is 5.4e-7 to 7.0e-7. The band is 10x the larger of the two sides' steps.
_PS_RESOLUTION = rung(4)
_PS_K_BAND = 1e-5
#: The shape: phi at PS-1982's nodes over phi at its first node, both sides. ``[M]``: 3.0e-6 (tauR 2.5) and 5.2e-6
#: (tauR 5) at rung 4 against PS-1982 at 40 nodes; 3.4e-5 and 5.5e-5 at rung 3. The band is 20x the larger.
_PS_SHAPE_BAND = 1e-4


@lru_cache(maxsize=None)
def _ps1982(case: tuple[float, float, float, float]):
    from orpheus.derivations.continuous.peierls_nystrom.ps1982_reference import solve_ps1982_vacuum_sphere

    radius, sigma_t, sigma_s, nu_sigma_f = case
    result = solve_ps1982_vacuum_sphere(R=radius, sigma_t=sigma_t, sigma_s=sigma_s, nu_sigma_f=nu_sigma_f,
                                        n_quad=_PS_NODES, max_iter=400, tol=_PS_TOL)
    assert result.converged, f"PS-1982 did not converge on {case} in {result.iterations} iterations"
    return result


@pytest.mark.l1
@pytest.mark.parametrize("case", _PS_CASES)
@pytest.mark.rests_on(_D5, _D8)
def test_the_vacuum_spheres_k_is_ps1982s(case) -> None:
    """[V] The one-group vacuum sphere's k is PS-1982's (Pomraning and Siewert 1982, Eq. (21)), to 1e-5 relative.

    PS-1982 solves the scalar-flux integral equation in its E_1 kernel by
    a Nystrom product quadrature: no lines, no boundary closure, no panel
    basis. It succeeds the old family's ``test_a2_variant_alpha_agrees_with_ps1982``
    (1e-4); the band is 10x the larger of the two sides' measured steps
    (module constants). First reds (``[M]`` 2026-10-10): the total cross
    section read 1e-5 large (arm ``total-scaled-1e-5``); the vacuum wall read
    as returning half of what reaches it (arm ``vacuum-half-mirror``).
    """
    radius, sigma_t, sigma_s, nu_sigma_f = case
    spec = _vacuum_sphere(radius, _one_group(sigma_t, sigma_s, nu_sigma_f))
    k = _read(spec, Eigenvalue(), _PS_RESOLUTION)
    k_ps = _ps1982(case).k_eff
    assert abs(k - k_ps) / k_ps < _PS_K_BAND, (k, k_ps)


@pytest.mark.l1
@pytest.mark.parametrize("case", _PS_CASES[:2])
@pytest.mark.rests_on(_HERE + "test_the_vacuum_spheres_k_is_ps1982s")
def test_the_vacuum_spheres_mode_is_ps1982s_shape(case) -> None:
    """[D10, against PS-1982] phi(r) / phi(r_0) at PS-1982's 40 nodes equals PS-1982's own, to 1e-4.

    The characteristic side transports its converged emission to each of
    PS-1982's nodes (in cm, the node's optical radius over Sigma_t), so the
    two shapes are compared point for point with no interpolation. It
    succeeds the old family's ``test_a2_phi_shape_qualitative_agreement``,
    which compared the centre-to-surface ratio to 20 %. First red
    (``[M]`` 2026-10-10): the vacuum wall returning half (arm
    ``vacuum-half-mirror``), which flattens the shape.
    """
    radius, sigma_t, sigma_s, nu_sigma_f = case
    spec = _vacuum_sphere(radius, _one_group(sigma_t, sigma_s, nu_sigma_f))
    ps = _ps1982(case)
    radii = np.asarray(ps.r_nodes_optical) / sigma_t
    phi = np.array([_read(spec, PointValue(float(r), 0), _PS_RESOLUTION) for r in radii])
    np.testing.assert_allclose(phi / phi[0], np.asarray(ps.phi) / ps.phi[0], rtol=0.0, atol=_PS_SHAPE_BAND)


@pytest.mark.l1
@pytest.mark.rests_on(_D8)
def test_a_thick_vacuum_spheres_k_is_just_below_k_inf() -> None:
    """[Q, the thick sphere] At an optical radius of 25 (R = 50 cm, c = 0.4), 0.95 k_inf < k < k_inf.

    The old row (``test_a2_thick_sphere_approaches_k_inf``) also compared PS-1982
    at 2e-3. That comparison is not re-pointed: PS-1982 does not converge
    there in 200 power iterations at 30 or 40 nodes (``[M]`` 2026-10-10,
    ``probes/probe_ps1982.log``: ``it=200``, and the characteristic reference
    at rungs 3 and 4 agree to 1.6e-7 while missing PS-1982 by 9.2e-4), so the
    old band measured PS-1982's iteration, not the transport. The monotone
    approach to k_inf is D8's thickening row. First red: a vacuum wall read as a
    mirror (arm ``vacuum-as-mirror``, k = k_inf).
    """
    sigma_t, sigma_s, nu_sigma_f = 0.5, 0.2, 0.025
    k_inf = nu_sigma_f / (sigma_t - sigma_s)
    k = _read(_vacuum_sphere(50.0, _one_group(sigma_t, sigma_s, nu_sigma_f)), Eigenvalue(), rung(3))
    assert 0.95 * k_inf < k < k_inf * (1.0 - 1e-6), (k, k_inf)


# ── Ua-1-0-CY: Sood's and WM-72's bare critical cylinder ─────────────────


def _last_digit(value: float) -> float:
    return 10.0 ** int(Decimal(repr(value)).as_tuple().exponent)


#: The fast leg's band: the reference's measured error at rung 2 on Ua-1-0-CY is -8.4e-6 (``[M]`` 2026-10-10,
#: ``probes/probe_ua_cy.log``: k - 1 = -8.4e-6, -2.7e-7 at rungs 2, 3, the slow row's estimate converging below 1e-6),
#: plus the truth's resolution (|dk/dmfp| 0.45 x 0.5e-8, 2.3e-9): 3e-5 is 3.6x the rung-2 error.
_UA_CY_FAST_BAND = 3e-5


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-cylinder-mr-wm72-vacuum")
@pytest.mark.rests_on(_D2)
def test_k_is_one_at_the_one_group_cylinders_published_critical_radius_in_the_fast_tier() -> None:
    """[V, the fast tier] Ua-1-0-CY (Sood 2003 Table 10, c = 1.30) at its printed critical radius reads k = 1 to 3e-5
    at rung 2 (about 3 s), the default-tier leg of the slow row below (qa's ruling of 2026-10-10: the cylinder keeps a
    default-tier eigenvalue row). The band is the measured rung-2 error with a 3.6x margin, a RECORD of the rung's
    accuracy against an independent truth, not the slow row's estimated bound. First red by value (``[M]``
    2026-10-10): the in-plane speed read 1e-4 large (qa's ``obliquity-dropped`` at QSPEED 1.0001, battery arm
    ``cylinder-speed-scaled``), which moves k by about 1e-4."""
    case = sood.UA_1_0_CY_STUB
    mixture = case.materials[0]
    printed = case.truth.critical_dimension_mfp
    assert printed is not None
    geometry = StructuredGeometry.cylinder((0.0, printed / float(mixture.SigT[0])), (0,), outer=BC.vacuum)
    k = _read(GeometrySpecification(Materials({0: mixture}), geometry, Eigen(_K)), Eigenvalue(), rung(2))
    assert abs(k - 1.0) < _UA_CY_FAST_BAND, k - 1.0


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("peierls-greens-cylinder-mr-wm72-vacuum")
@pytest.mark.rests_on(_D2)
def test_k_is_one_at_the_one_group_cylinders_published_critical_radius() -> None:
    """[V; ``peierls-greens-cylinder-mr-wm72-vacuum``] Sood 2003 Table 10, Ua-1-0-CY (c = 1.30): k = 1 at the
    printed critical radius 1.72500292 mfp, within the truth's resolution plus the reference's error estimate;
    and WM-72's singular-eigenfunction critical radius is Sood's to 2e-6.

    Two independent truths of one quantity: Sood's printed digits, and
    WM-72's F_N solve (no Bickley-Naylor function anywhere). The band is
    |dk/dmfp| x half a unit of the printed last digit (the slope measured
    here) plus the reference's error at rung 3, its step to rung 4 over one
    minus the steps' ratio (``_ladder_rules.geometric_error``): a printed
    truth decides only beyond its resolution (lessons, the published-truth
    rule). It succeeds the old family's ``test_mr_single_region_vacuum_matches_wm72``
    and ``test_a2_variant_alpha_agrees_with_sood2003_cylinder`` (1e-5 each).
    ``[M]`` 2026-10-10 (``probes/probe_ua_cy.log``): k - 1 = -8.4e-6, -2.7e-7
    at rungs 2, 3; the slope 0.45 per mfp. First red: the cylinder's
    obliquity dropped from the in-plane lengths (arm ``cylinder-no-obliquity``).
    """
    from orpheus.derivations.continuous.singular_eigenfunction import (
        solve_singular_eigenfunction_cylinder_bare_critical,
    )

    case = sood.UA_1_0_CY_STUB
    mixture = case.materials[0]
    sigma_t = float(mixture.SigT[0])
    printed = case.truth.critical_dimension_mfp
    wm = solve_singular_eigenfunction_cylinder_bare_critical(c=1.30, sigma_t=sigma_t, n_grid=24)
    assert printed is not None and wm.converged and wm.r_c_cm is not None
    assert abs(wm.r_c_cm * sigma_t / printed - 1.0) < 2e-6, (wm.r_c_cm * sigma_t, printed)

    def k_at(radius_mfp: float, p: int) -> float:
        geometry = StructuredGeometry.cylinder((0.0, radius_mfp / sigma_t), (0,), outer=BC.vacuum)
        return _read(GeometrySpecification(Materials({0: mixture}), geometry, Eigen(_K)), Eigenvalue(), rung(p))

    k2, k3, k4 = (k_at(printed, p) for p in (2, 3, 4))
    slope = (k_at(printed * (1.0 + 1e-4), 2) - k2) / (printed * 1e-4)
    band = abs(slope) * 0.5 * _last_digit(printed) + geometric_error(k4 - k3, k3 - k2)
    assert abs(k3 - 1.0) < band, (k3 - 1.0, band)
