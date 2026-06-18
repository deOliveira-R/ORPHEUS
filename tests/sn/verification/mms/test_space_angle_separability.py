r"""Characterization gate — (spatial :math:`\otimes` angular) discretization-error
separability of ORPHEUS SN.  Issue #236 Phase 3 (sub-task ST5).

This is the LAST piece of the #236 (spatial :math:`\otimes` angular product)
campaign (Phase 1 pairing-validity + Phase 2 :math:`\tau`-carve already shipped).
It pins, in permanent test form, the campaign's headline claim — backed by the
literature (LMM-1987 spatial / BMC-2010 angular diffusion-limit split) AND by the
empirical probing of ``diag_sep_*.py`` — that the SN discretization is a tensor
product of two independently-selectable axes (a SPATIAL closure, mesh refinement
``h`` = ``n_cells``; an ANGULAR factor, quadrature order ``N`` = ``n_ordinates`` /
``n_phi``), whose error STRUCTURE differs by geometry:

* **Cartesian (slab / slab-iso):** the :math:`(1-\mu^2)/r` angular-redistribution
  term is ABSENT, so space :math:`\perp` angle.  The error SEPARATES additively:
  :math:`E(h,N) \approx E_{\rm space}(h) + E_{\rm angle}(N)`.  Signature: the
  spatial convergence RATE is the SAME at every quadrature order (N-independent),
  and the **mixed second difference**
  :math:`M = E[h_1,N_1] - E[h_1,N_2] - E[h_2,N_1] + E[h_2,N_2]` vanishes relative
  to the individual deltas.

* **Curvilinear (sphere / cylinder):** the Morel–Montry angular thread couples
  space and angle through the shared cell-balance denominator (see the coupling
  memo + ``docs/theory/discrete_ordinates.rst`` :eq:`dd-curvilinear-scalar` /
  :eq:`mm-weights`).  The error does NOT separate additively — it GATES:
  :math:`E(h,N) \approx \max(E_{\rm space}(h), E_{\rm angle}(N))`.  You cannot
  harvest fine-``h`` accuracy without also refining ``N``.  Signature: at a coarse
  quadrature the spatial refinement SATURATES (the h-ratio collapses toward 1.0)
  while at a fine quadrature the same spatial ladder keeps O(h²); + a LARGE mixed
  second difference.

What this gate is and is NOT
============================

This is a **characterization gate**, modelled on
``test_curvilinear_pole_cell_characterization.py`` (#233).  It pins what is TRUE
about the error STRUCTURE (the regime DISCRIMINATION) WITHOUT calcifying the
limitation: the curvilinear gating is the standing behaviour today; if a future
2-D angular closure (#229) or a higher-order spatial scheme (#158 / #6) lifts the
gating, the discrimination assertions are designed to stay meaningful (the
Cartesian leg pins separability as a positive CLAIM; the curvilinear leg pins the
gating as a positive SIGNATURE — neither is an xfail-pending-fix).

It is **NOT** a convergence-order verification of the solver value (that lives in
``test_curvilinear_aniso_convergence.py`` + ``test_mms.py``).  It is NOT an
eigenvalue claim (vv-principles: MMS does not reach the eigenvalue layer; the
measurand here is a discretization-error-surface STRUCTURE, a math claim).  The
V&V level is **L1** (an MMS-convergence-structure / math claim, the same level as
the pole-cell characterization gate it mirrors).

The L27 scalar-vs-per-ordinate decision (documented)
====================================================

The signature is built from the SCALAR (weight-summed) flux
``res.scalar_flux.values[0, :]``.  vv-principles L27 / Mode-8 H3 warns that a
weight-summed metric is BLIND BY CONSTRUCTION to a wrong angular closure (the M-M
:math:`\alpha`-dome telescopes under the weight sum:
:math:`\sum_n w_n \psi_n` holds even when per-ordinate balance is wrong).  Because
the curvilinear GATING is itself an angular-closure phenomenon, this gate adds a
PER-ORDINATE leg (``test_curvilinear_gating_per_ordinate_not_blind``) that
reproduces the sphere gating signature from the max-over-ordinates per-ordinate
L2, so the gate cannot be blind to a future angular-closure regression that the
scalar reconstruction would telescope away.

⭐ CONVENTION TRAP (measured 2026-06-18): ``case.psi_exact(r, mu_n)`` returns
:math:`A(r) + B(r)\mu_n` **WITHOUT** the :math:`1/W` factor (``W = \sum_n w_n``),
by its own documented design — "Returned without the 1/W factor for caller
convenience".  The solver stores the per-ordinate flux WITH the :math:`1/W` factor
(``\sum_n w_n \psi_n^{\rm solver} = \phi`` exactly).  The per-ordinate metric MUST
divide ``psi_exact`` by ``W`` before comparing, else a 2× (= W) normalization
mismatch swamps the metric (a naive per-ordinate L2 reads ~8.1 with FLAT h-ratios
— a measured trap this docstring records so it is never re-introduced).

Measured evidence (this branch, ``feature/sn-spatial-angular-product`` @
``cdb8f07``, ``.venv/bin/python``, 2026-06-18; the four ``diag_sep_*.py`` probes
reproduced the documented cross-terms bit-for-bit post the Phase-2 :math:`\tau`
carve — slab-iso ``|M|/max`` 0.000–0.005, cylinder 0.019/0.083, sphere 0.373
@ N8→16):

Reduced grid ``nx \in {20,40,80}``, total wall-clock ~2.0 s (NO ``@slow`` needed,
but marked ``@slow`` anyway — the curvilinear solves dominate and the gate is a
characterization gate, not a fast canary):

* SPHERE  scalar L2:  N=8  [1.4665e-2, 5.4016e-3, 4.6851e-3]  h-ratios [2.71, 1.15]
                       N=32 [1.5006e-2, 3.7136e-3, 9.2864e-4]  h-ratios [4.04, 4.00]
                       cross-term |M|/max (nx20→80 × N8→32) = 0.411
          per-ord L2:  N=8  [2.2897e-2, 1.4310e-2, 1.2383e-2]  h-ratios [1.60, 1.16]
                       N=32 [3.3372e-2, 9.8571e-3, 3.3207e-3]  h-ratios [3.39, 2.97]
* CYLINDER scalar L2:  n_phi=8  [1.9504e-2, 1.9116e-2, 1.9038e-2]  h-ratios [1.02, 1.00]
                       n_phi=16 [8.0504e-3, 7.4726e-3, 7.3719e-3]
                       azimuthal floor drop n_phi 8→16 @ nx=80 = 2.58×
                       cross-term |M|/max = 0.019 (small only because E≈E_angle swamps)
* SLAB-ISO scalar L2:  N=4  [5.4249e-4, 1.3540e-4, 3.3837e-5]  h-ratios [4.01, 4.00]
                       N=16 [5.3992e-4, 1.3477e-4, 3.3679e-5]  h-ratios [4.01, 4.00]
                       cross-term |M|/max = 0.0047
* SLAB-P1 scalar L2:   N=4  [6.7449e-3, 6.7837e-3, 6.7973e-3]  (flat — angular floor)
                       N=16 [6.7325e-3, 6.7715e-3, 6.7851e-3]
                       cross-term |M|/max = 0.0038

References
----------
* ``.claude/agent-memory/numerics-investigator/sn_space_angle_discretization_coupling.md``
  — the measured coupling (Cartesian separable, curvilinear gates).
* ``.claude/agent-memory/literature-researcher/space_angle_discretization_separability.md``
  — LMM-1987 (spatial) × BMC-2010 (angular) diffusion-limit split; WHY Cartesian
  separates and curvilinear gates.
* Issue #229 (azimuthal floor), #233 (pole-cell O(h)), #236 (the campaign),
  ERR-026 (the curvilinear M-M closure the curvilinear legs activate).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.mms.sn import (
    build_1d_slab_mms_case,
    build_cylindrical_anisotropic_mms_case,
    build_p1_aniso_mms_case,
    build_spherical_anisotropic_mms_case,
)
from orpheus.sn import solve_sn_fixed_source

from tests.sn._test_helpers import scalar_flux_l2_ladder, volume_weighted_l2

# Every leg verifies the space⊗angle decomposition law (eq. sn-space-angle-
# separability) from a different angle; the cross-term legs ALSO verify the
# mixed-second-difference discriminator (eq. sn-space-angle-cross-term).
pytestmark = [
    pytest.mark.l1,
    pytest.mark.verifies("sn-space-angle-separability"),
]

# Reduced grid: nx ∈ {20,40,80} resolves both the gated (collapse-to-1) and the
# O(h²) (ratio≈4) regimes cleanly while keeping the whole module ~2 s wall-clock.
_N_CELLS = (20, 40, 80)


# ──────────────────────────────────────────────────────────────────────────
# Shared machinery (mirrors test_curvilinear_aniso_convergence.py idioms).
# The error norm is the single-source ``volume_weighted_l2`` primitive
# (tests/sn/_test_helpers.py) — NOT a 7th private copy of the mms/ ``_l2_1d``.
# ──────────────────────────────────────────────────────────────────────────

def _sphere_per_ordinate_max_l2_ladder(case, n_cells) -> np.ndarray:
    r"""Max-over-ordinates per-ordinate volume-weighted L2 error ladder (the L27
    angular-aware metric).

    ⭐ ``case.psi_exact(r, mu_n)`` returns :math:`A + B\mu_n` WITHOUT the
    :math:`1/W` factor (its documented contract); the solver stores the
    per-ordinate flux WITH :math:`1/W`.  Divide the reference by ``W`` here, else
    a 2× normalization mismatch swamps the metric (see the module docstring trap).
    """
    mu = case.quadrature.mu_x
    sum_w = float(case.quadrature.weights.sum())
    errors = []
    for nc in n_cells:
        mesh = case.build_mesh(nc)
        Q = case.external_source(mesh)
        result = solve_sn_fixed_source(
            case.materials, mesh, case.quadrature, Q,
            max_inner=500, inner_tol=1e-13,
        )
        psi = result.angular_flux.bulk.values  # (N, ng, nx)
        centers = mesh.centers
        V = mesh.volumes
        worst = 0.0
        for n in range(psi.shape[0]):
            psi_ref_n = case.psi_exact(centers, float(mu[n])) / sum_w
            worst = max(worst, volume_weighted_l2(psi[n, 0, :], psi_ref_n, V))
        errors.append(worst)
    return np.asarray(errors)


def _h_ratios(errors: np.ndarray) -> np.ndarray:
    """E[h]/E[h/2] per consecutive pair — ~4 for O(h²), ~1 for a saturated floor."""
    return errors[:-1] / errors[1:]


def _mixed_second_difference(L2: np.ndarray) -> tuple[float, float, float, float]:
    r"""The cross-term :math:`M = E[h_1,N_1]-E[h_1,N_2]-E[h_2,N_1]+E[h_2,N_2]` over
    the COARSEST-h × FINEST-h corners (rows 0 and -1) and the two columns of a
    2-column error table.  Returns ``(M, dEh, dEN, |M|/max(dEh,dEN))``.

    ``dEh`` = spatial delta at the coarse-N column; ``dEN`` = angular delta at the
    coarse-h row.  ``|M|/max`` ≪ 1 ⟺ additively separable (no interaction).
    """
    assert L2.shape[1] == 2, (
        f"_mixed_second_difference needs a 2-column (two-quadrature) error "
        f"table; got shape {L2.shape} — a 3+-column caller would silently drop "
        f"a column."
    )
    i1, i2 = 0, L2.shape[0] - 1
    M = L2[i1, 0] - L2[i1, 1] - L2[i2, 0] + L2[i2, 1]
    dEh = abs(L2[i1, 0] - L2[i2, 0])
    dEN = abs(L2[i1, 0] - L2[i1, 1])
    rel = abs(M) / max(dEh, dEN, 1e-300)
    return float(M), float(dEh), float(dEN), float(rel)


# A clean discrimination band.  The Cartesian cross-term is measured ≤ 0.005;
# the curvilinear gating cross-term is ≥ 0.41 (sphere) — three orders of
# headroom.  The bounds below sit comfortably between, so the gate DISCRIMINATES
# separable from gated rather than pinning a brittle exact number.
_CARTESIAN_SEPARABLE_MAX = 0.05   # |M|/max ABOVE which Cartesian is no longer separable
_CURVILINEAR_GATING_MIN = 0.15    # |M|/max BELOW which the sphere is no longer gating


# ══════════════════════════════════════════════════════════════════════════
# Cartesian legs — space ⊥ angle SEPARATES additively (the robust easy leg)
# ══════════════════════════════════════════════════════════════════════════

@pytest.mark.slow
@pytest.mark.verifies("sn-space-angle-cross-term")
def test_cartesian_slab_iso_space_angle_separable():
    r"""[Cartesian, separable] The slab ISOTROPIC MMS spatial O(h²) RATE is
    N-INDEPENDENT and the mixed second difference vanishes.

    This is the cleanest separability signature: the slab-iso case has a real
    O(h²) window, so the spatial h-ratio is ≈4 at EVERY quadrature order — the
    angular axis does not touch the spatial convergence.  The
    :math:`(1-\mu^2)/r` term is absent in Cartesian geometry, so
    :math:`E(h,N) \approx E_{\rm space}(h) + E_{\rm angle}(N)` and
    :math:`\partial^2 E/\partial h\,\partial N \approx 0`.

    Measured 2026-06-18 (nx[20,40,80]):
      N=4  h-ratios [4.01, 4.00];  N=16 h-ratios [4.01, 4.00];  |M|/max = 0.0047.
    A violation (N-dependent spatial rate, or |M|/max ≳ 0.05) would mean the
    Cartesian path acquired a space–angle coupling it must not have.
    """
    L2 = np.column_stack([
        scalar_flux_l2_ladder(build_1d_slab_mms_case(n_ordinates=4), _N_CELLS),
        scalar_flux_l2_ladder(build_1d_slab_mms_case(n_ordinates=16), _N_CELLS),
    ])
    r_n4 = _h_ratios(L2[:, 0])
    r_n16 = _h_ratios(L2[:, 1])
    M, dEh, dEN, rel = _mixed_second_difference(L2)
    print(f"slab-iso N4 ratios={r_n4} N16 ratios={r_n16} |M|/max={rel:.4f}")

    # Spatial rate is O(h²) and the SAME at both quadrature orders (N-independent).
    assert np.all(r_n4 > 3.5), (
        f"slab-iso N=4 spatial rate not O(h²): h-ratios={r_n4} (want >3.5). "
        f"A coupled Cartesian path would degrade the rate at coarse N."
    )
    assert np.all(r_n16 > 3.5), (
        f"slab-iso N=16 spatial rate not O(h²): h-ratios={r_n16} (want >3.5)."
    )
    assert np.allclose(r_n4, r_n16, rtol=0.05), (
        f"slab-iso spatial rate is N-DEPENDENT (N4={r_n4} vs N16={r_n16}) — "
        f"Cartesian space and angle must be independent; a rate that moves with "
        f"N means a space–angle coupling crept into the Cartesian path."
    )
    assert rel < _CARTESIAN_SEPARABLE_MAX, (
        f"slab-iso cross-term |M|/max={rel:.4f} ≥ {_CARTESIAN_SEPARABLE_MAX} — "
        f"Cartesian space×angle is no longer additively separable "
        f"(M={M:+.3e}, dEh={dEh:.3e}, dEN={dEN:.3e})."
    )


@pytest.mark.slow
@pytest.mark.verifies("sn-space-angle-cross-term")
def test_cartesian_slab_p1_aniso_floor_n_independent():
    r"""[Cartesian, separable — angular axis ACTIVE] The slab P1-ANISOTROPIC MMS
    spatial-error floor is N-independent and the mixed second difference vanishes.

    Companion to the slab-iso leg with the ANGULAR axis genuinely active (the
    P1 :math:`\mu`-anisotropic ansatz).  This P1-aniso MMS has NO h-convergence
    window — the error sits at a flat ~6.8e-3 angular/MMS floor at every ``h`` —
    so the separability signature is that this floor is the SAME at every
    quadrature order (the angular anisotropy does not induce an h-dependence) and
    the cross-term stays ≈0.  This confirms separability survives an active
    angular term, not just the isotropic degenerate case.

    Measured 2026-06-18 (nx[20,40,80]):
      N=4  [6.745e-3, 6.784e-3, 6.797e-3];  N=16 [6.733e-3, 6.772e-3, 6.785e-3];
      |M|/max = 0.0038.  Floor differs N4 vs N16 by < 0.3 %.
    """
    L2 = np.column_stack([
        scalar_flux_l2_ladder(build_p1_aniso_mms_case(n_ordinates=4), _N_CELLS),
        scalar_flux_l2_ladder(build_p1_aniso_mms_case(n_ordinates=16), _N_CELLS),
    ])
    M, dEh, dEN, rel = _mixed_second_difference(L2)
    floor_n4 = L2[-1, 0]
    floor_n16 = L2[-1, 1]
    print(f"slab-P1 floor N4={floor_n4:.4e} N16={floor_n16:.4e} |M|/max={rel:.4f}")

    assert abs(floor_n4 - floor_n16) < 0.05 * floor_n4, (
        f"slab-P1 spatial floor is N-DEPENDENT (N4={floor_n4:.4e} vs "
        f"N16={floor_n16:.4e}, >5 % apart) — the active P1 angular axis must not "
        f"induce an h-dependence in Cartesian geometry."
    )
    assert rel < _CARTESIAN_SEPARABLE_MAX, (
        f"slab-P1 cross-term |M|/max={rel:.4f} ≥ {_CARTESIAN_SEPARABLE_MAX} — "
        f"the active angular axis induced a space–angle coupling in Cartesian "
        f"(M={M:+.3e}, dEh={dEh:.3e}, dEN={dEN:.3e})."
    )


# ══════════════════════════════════════════════════════════════════════════
# Curvilinear legs — space ⊗ angle GATES (the hard discriminating leg)
# ══════════════════════════════════════════════════════════════════════════

@pytest.mark.slow
@pytest.mark.catches("ERR-026")
def test_sphere_spatial_rate_is_quadrature_gated():
    r"""[Curvilinear, gating — THE discriminator] On the spherical anisotropic
    MMS the spatial convergence RATE is GATED by the quadrature order: at a
    COARSE quadrature (N=8) the spatial h-ratio collapses toward 1 (refinement
    SATURATES at the angular floor); at a FINE quadrature (N=32) the SAME spatial
    ladder recovers O(h²).

    This is the robust POSITIVE signature of gating (option (i) of the ST5
    brief): not "the cross-term is large" (a fragile negative claim) but "the
    spatial rate is N-gated" (a direct, mechanism-anchored measurement).  The M-M
    angular thread feeds the spatial cell-balance denominator, so the angular
    interpolation floor CAPS what the spatial scheme can resolve —
    :math:`E(h,N) \approx \max(E_{\rm space}(h), E_{\rm angle}(N))`.  You cannot
    harvest fine-``h`` accuracy at coarse ``N``; that is the campaign's headline.

    Measured 2026-06-18 (nx[20,40,80]):
      N=8  scalar [1.4665e-2, 5.4016e-3, 4.6851e-3]  h-ratios [2.71, 1.15]
      N=32 scalar [1.5006e-2, 3.7136e-3, 9.2864e-4]  h-ratios [4.04, 4.00]
    The FINEST coarse-N ratio (1.15) is the saturation; the fine-N ratios (≈4)
    are O(h²).  ``@catches("ERR-026")``: the curvilinear M-M closure is active.

    NOTE (vv anti-pattern #5/#17): the gating bounds are designed to stay
    meaningful if #229 lifts the floor — a future 2-D angular closure would
    raise the coarse-N saturated ratio toward 4, at which point THIS test would
    redden and must be RE-TUNED to the new (better) regime, not deleted.  That
    re-tune is the intended signal that the gating was lifted, not a regression.
    """
    L2_n8 = scalar_flux_l2_ladder(build_spherical_anisotropic_mms_case(n_ordinates=8), _N_CELLS)
    L2_n32 = scalar_flux_l2_ladder(build_spherical_anisotropic_mms_case(n_ordinates=32), _N_CELLS)
    r_n8 = _h_ratios(L2_n8)
    r_n32 = _h_ratios(L2_n32)
    print(f"sphere gating: N8 err={L2_n8} ratios={r_n8}  N32 err={L2_n32} ratios={r_n32}")

    # Coarse-N spatial refinement SATURATES: the finest h-ratio collapses toward 1.
    assert r_n8[-1] < 1.6, (
        f"sphere N=8 spatial refinement did NOT saturate: finest h-ratio "
        f"{r_n8[-1]:.2f} ≥ 1.6 (measured 1.15).  If this rose toward 4, the "
        f"#229 angular floor was LIFTED — the gating is GONE; re-tune this gate "
        f"to the new regime (it is not a regression)."
    )
    # Fine-N recovers O(h²): the SAME spatial ladder converges when N is refined.
    assert np.all(r_n32 > 3.0), (
        f"sphere N=32 spatial rate not O(h²): h-ratios={r_n32} (want >3.0).  "
        f"At a fine quadrature the angular floor drops below the spatial error, "
        f"so the spatial ladder MUST recover its O(h²) rate; if it does not, the "
        f"spatial closure itself regressed (not a gating phenomenon)."
    )
    # DISCRIMINATION: the rate genuinely depends on N (the defining gating fact).
    assert r_n32[-1] > 2.0 * r_n8[-1], (
        f"sphere spatial rate is NOT quadrature-gated: finest h-ratio N32="
        f"{r_n32[-1]:.2f} is not ≳2× N8={r_n8[-1]:.2f} — the space/angle error "
        f"is separable here, which contradicts the curvilinear coupling claim."
    )


@pytest.mark.slow
@pytest.mark.verifies("sn-space-angle-cross-term")
def test_sphere_cross_term_large_discriminates_from_cartesian():
    r"""[Curvilinear, gating] The spherical anisotropic MMS mixed second
    difference :math:`|M|/\max(\Delta E_h, \Delta E_N)` is LARGE — comfortably
    above the Cartesian separable ceiling — so the cross-term DISCRIMINATES
    gated curvilinear from separable Cartesian.

    This corroborates :func:`test_sphere_spatial_rate_is_quadrature_gated` with
    the cross-term estimator directly.  The bound is a LOWER bound (``> 0.15``)
    well above the measured Cartesian upper bound (``≤ 0.005``) and well below
    the measured sphere value (``0.41``), so it is a DISCRIMINATION threshold,
    not a brittle "this is exactly 0.41" pin.

    Measured 2026-06-18 (nx20→80 corners × N8→32): |M|/max = 0.411
    (M=-6.13e-3-scale; dEh=9.98e-3, dEN=3.41e-4).

    NO ``@catches("ERR-026")``: the M-M closure is active here, but a re-emerged
    ERR-026 (a divergent wrong fixed point) does not have a MUTATION-PROVEN
    direct effect on this cross-term assertion — the proven ERR-026 catcher is
    the O(h²)-recovery assertion in
    :func:`test_sphere_spatial_rate_is_quadrature_gated` (a divergent fine-N
    ladder makes ``r_n32 > 3.0`` go RED; verified 2026-06-18).  Attaching the
    marker here would be a phantom coverage claim (vv-principles: a ``catches``
    marker is a coverage claim, not a topic tag).
    """
    L2 = np.column_stack([
        scalar_flux_l2_ladder(build_spherical_anisotropic_mms_case(n_ordinates=8), _N_CELLS),
        scalar_flux_l2_ladder(build_spherical_anisotropic_mms_case(n_ordinates=32), _N_CELLS),
    ])
    M, dEh, dEN, rel = _mixed_second_difference(L2)
    print(f"sphere cross-term |M|/max={rel:.3f} (M={M:+.3e} dEh={dEh:.3e} dEN={dEN:.3e})")

    assert rel > _CURVILINEAR_GATING_MIN, (
        f"sphere cross-term |M|/max={rel:.3f} ≤ {_CURVILINEAR_GATING_MIN} — the "
        f"space×angle error is too SEPARABLE for a gating regime (Cartesian sits "
        f"at ≤0.005, gating sphere at ≈0.41).  Either the gating was lifted "
        f"(#229 — re-tune) or the curvilinear coupling regressed."
    )


@pytest.mark.slow
def test_cylinder_spatial_saturates_at_azimuthal_floor():
    r"""[Curvilinear, gating — cylinder] On the cylindrical anisotropic MMS the
    spatial refinement SATURATES at the azimuthal angular floor: at fixed
    ``n_phi`` the spatial h-ratio is ≈1 (no benefit from refining ``h``), and the
    floor DROPS when the AZIMUTHAL quadrature ``n_phi`` doubles.

    The cylinder is the extreme of the gating regime — it has NO pre-floor O(h²)
    window at any practical quadrature (the (η,φ) angular variation a 1-D η-march
    cannot thread exactly, #229), so ``E ≈ E_angle(n_phi)`` and the spatial axis
    is fully gated.  The positive signatures: (a) spatial saturation (h-ratio≈1
    at fixed n_phi), (b) the floor scales with the AZIMUTHAL order.  Together
    they show space is slaved to the azimuthal angular floor.

    Measured 2026-06-18 (n_mu=4, nx[20,40,80]):
      n_phi=8  [1.9504e-2, 1.9116e-2, 1.9038e-2]  h-ratios [1.02, 1.00]
      n_phi=16 [8.0504e-3, 7.4726e-3, 7.3719e-3]
      azimuthal floor drop n_phi 8→16 @ nx=80 = 2.58×.

    NO ``@catches("ERR-026")``: the cylinder M-M closure is active, but the
    MUTATION-PROVEN ERR-026 catcher is the sphere O(h²)-recovery assertion (see
    :func:`test_sphere_spatial_rate_is_quadrature_gated`); the cylinder's
    band-level ERR-026 catcher already lives in
    ``test_curvilinear_aniso_convergence.py``
    (``test_sn_cylindrical_aniso_mms_converges_second_order``).  Attaching the
    marker here without a proven direct reddening would inflate the per-ERR
    coverage count with a non-catcher (vv-principles).
    """
    case8 = build_cylindrical_anisotropic_mms_case(n_mu=4, n_phi=8)
    case16 = build_cylindrical_anisotropic_mms_case(n_mu=4, n_phi=16)
    L2_p8 = scalar_flux_l2_ladder(case8, _N_CELLS)
    L2_p16 = scalar_flux_l2_ladder(case16, _N_CELLS)
    r_p8 = _h_ratios(L2_p8)
    floor_drop = L2_p8[-1] / L2_p16[-1]
    print(f"cyl: n_phi8 err={L2_p8} ratios={r_p8}  n_phi16 err={L2_p16}  "
          f"floor_drop={floor_drop:.2f}x")

    # Spatial refinement is GATED: at fixed azimuthal order it does not converge.
    assert np.all(r_p8 < 1.5), (
        f"cyl n_phi=8 spatial refinement did NOT saturate: h-ratios={r_p8} "
        f"(want all <1.5, measured [1.02, 1.00]).  A ratio rising toward 4 would "
        f"mean the azimuthal gating was lifted (#229 — re-tune)."
    )
    # The floor IS the azimuthal angular floor: it drops when n_phi doubles.
    assert floor_drop > 2.0, (
        f"cyl azimuthal floor did NOT scale with n_phi: drop {floor_drop:.2f}× "
        f"(< 2.0) when n_phi 8→16 — the saturated floor is not the #229 "
        f"azimuthal-thread floor; investigate a fixed closure artefact."
    )


@pytest.mark.slow
@pytest.mark.catches("ERR-026")
def test_curvilinear_gating_per_ordinate_not_blind():
    r"""[Curvilinear, gating — L27 angular-aware] The sphere gating signature
    reproduces from the PER-ORDINATE (max-over-ordinates) L2, so the gate is NOT
    blind to a wrong angular closure that the scalar weight-sum would telescope
    away.

    vv-principles L27 / H3: the scalar (weight-summed) metric is blind by
    construction to a per-ordinate closure error (the M-M :math:`\alpha`-dome
    telescopes under :math:`\sum_n w_n \psi_n`).  Because the curvilinear gating
    is ITSELF an angular-closure phenomenon, this leg measures the gating from
    the per-ordinate flux — the metric the scalar legs are blind to — and
    confirms the SAME N-gated spatial rate: coarse N saturates, fine N recovers.

    ⭐ Uses the :math:`1/W`-corrected reference (``psi_exact / sum_w``); see the
    module docstring convention trap.

    Measured 2026-06-18 (nx[20,40,80], per-ord max-over-ordinates L2):
      N=8  [2.2897e-2, 1.4310e-2, 1.2383e-2]  h-ratios [1.60, 1.16]
      N=32 [3.3372e-2, 9.8571e-3, 3.3207e-3]  h-ratios [3.39, 2.97]
    The per-ordinate finest coarse-N ratio (1.16) saturates; fine-N (≈3) converges
    — the gating is visible in the angular-aware metric, not only the scalar.
    ``@catches("ERR-026")``: the curvilinear M-M closure is active.
    """
    case8 = build_spherical_anisotropic_mms_case(n_ordinates=8)
    case32 = build_spherical_anisotropic_mms_case(n_ordinates=32)
    po_n8 = _sphere_per_ordinate_max_l2_ladder(case8, _N_CELLS)
    po_n32 = _sphere_per_ordinate_max_l2_ladder(case32, _N_CELLS)
    r_n8 = _h_ratios(po_n8)
    r_n32 = _h_ratios(po_n32)
    print(f"sphere per-ord gating: N8 err={po_n8} ratios={r_n8}  "
          f"N32 err={po_n32} ratios={r_n32}")

    # Sanity: the per-ordinate metric is correctly normalized (O(1e-2), not O(8)).
    assert po_n8[0] < 1.0, (
        f"per-ordinate L2 {po_n8[0]:.3e} is O(1) not O(1e-2) — the 1/W "
        f"normalization of psi_exact was NOT applied (the convention trap)."
    )
    # Coarse-N per-ordinate refinement saturates.
    assert r_n8[-1] < 1.6, (
        f"sphere per-ordinate N=8 did NOT saturate: finest h-ratio {r_n8[-1]:.2f} "
        f"≥ 1.6 (measured 1.16) — gating not visible in the angular-aware metric."
    )
    # Fine-N per-ordinate converges (toward O(h²)); the gating is N-dependent.
    assert r_n32[-1] > 2.0, (
        f"sphere per-ordinate N=32 finest h-ratio {r_n32[-1]:.2f} ≤ 2.0 "
        f"(measured 2.97) — the per-ordinate spatial rate must recover when the "
        f"quadrature is refined, the angular-aware proof of gating."
    )
    assert r_n32[-1] > 1.5 * r_n8[-1], (
        f"sphere per-ordinate spatial rate not quadrature-gated: N32="
        f"{r_n32[-1]:.2f} not ≳1.5× N8={r_n8[-1]:.2f} — angular-aware metric "
        f"does not see the gating."
    )
