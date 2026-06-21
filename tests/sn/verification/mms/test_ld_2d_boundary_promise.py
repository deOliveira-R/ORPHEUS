r"""#257 S9 — the LD 2-D Cartesian BOUNDARY coherent-promise gate + sub-floor verdict pins.

Does the 2-D Cartesian LD boundary transverse-SLOPE moment improve the
converged-flux error magnitude or spatial-convergence ORDER above the O(h^2)
floor, in ANY physically-meaningful regime — INCLUDING the optically-thin /
streaming-dominated limit?

VERDICT (see the test bodies for the verbatim convergence tables): NO. Including
the boundary transverse-slope moment (``mom``) NEVER beats the average-only drop
(``flat``) above the O(h^2) floor — not at Sigma_t*L = 0.1 (streaming), 0.5,
1.0, 2.0, nor with grazing-heavy S8 quadrature, nor with the boundary slope
amplified 20x relative to the average.  In every regime flat == mom == flip at
O(h^2) globally; the slope is a sub-floor perturbation that, if anything, makes
the converged flux SLIGHTLY WORSE (it adds its own O(h^2) discretization noise to
a boundary representation the average moment already resolves to O(h^2)).

MECHANISM (the decisive scaling analysis, test_first_cell_row_already_second_order):
the FIRST-CELL-ROW error (the cells that directly consume the boundary inflow) is
ALREADY O(h^2) with the average-only inflow.  There is NO O(h) deficiency in the
converged flux for the slope to fix.  The slope DOES lift the *inflow
representation* O(h)->O(h^2), but that representation refinement is a second-order
correction to an ALREADY-second-order face balance, so it cannot move the
converged flux above the bulk O(h^2) floor.  This holds in the streaming limit
too: even when the inflow propagates ballistically across the domain, the LD cell
balance integrates the inflow against the cell's own linear basis, and the
cell-AVERAGE moment is what that integral needs to O(h^2).

THE COHERENT PROMISE (#257 S9 / #263).  "LD gives 2nd order EVERYWHERE including
the boundary" is TRUE and ALREADY DELIVERED — by the AVERAGE moment alone (the
prescribed inflow is exact at the face cells; the bulk LD closure carries it
inward at O(h^2)).  The transverse boundary-slope moment is a representation
refinement (O(h)->O(h^2) on the inflow trace) that is genuinely consumed (a flip
moves the flux ~4.1e-3 near-bdy >> tol, ``test_mms_ld_2d`` consumption gate) and
genuinely threaded (machine precision, ``test_mms_ld_2d`` threading gate) but its
converged-flux contribution is sub-floor.  S9 LOCKS this promise with a gate; it
does NOT fix it.  This is the canonical Mode-10 resolution where the
companion-isolating value gate is UNAVAILABLE (the boundary transverse-slope has
NO regime where it is the dominant forcing in the converged flux — the boundary
is codimension-1, measure-zero in the refinement limit).  The boundary moment is
a PROPERTY, not a type (#263 defers the first-class ``SpatialMomentField`` to the
collocation trigger).

If a verdict pin ever FAILS (flat and mom diverge in order, or mom beats flat
above the floor by more than the documented sub-floor band), the LD boundary
closure has changed and the S9 verdict must be revisited.

Promoted from ``derivations/diagnostics/diag_s9_ld_boundary_slope_optical_sweep.py``
(numerics-investigator, 2026-06-20).  Mode-8: the order-band checks are
``np.testing.*`` / ``pytest.fail`` (function calls that fire under ``python -O``),
NEVER bare ``assert`` (which ``-O`` strips to a no-op).  ``_solve_with_boundary_slope``
is imported as a SIBLING from ``test_mms_ld_2d`` (the public-API toggle the #251
Leg B gates already drive); the verdict harness is the same module's
``_face_transverse_buffers`` + ``leggauss`` (L11, never a production projector).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.mms.sn import (
    SN2DCartesianLDStressMMSCase,
    build_2d_cartesian_ld_stress_mms_case,
)
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn import loss_representation as LR

from .test_mms_ld_2d import _solve_with_boundary_slope


# ── Custom optical-depth-tunable streaming MMS case ──────────────────────────


def _streaming_case(sigt, c_ratio, *, length=1.0, order=4, bc_scale=1.0):
    """An ``SN2DCartesianLDStressMMSCase`` with UNIFORM, scalable cross sections.

    ``sigt`` sets the optical scale (``Sigma_t * length`` is the optical
    thickness), ``c_ratio`` the within-group scattering ratio (``c = Sigma_s /
    Sigma_t``; ``c_ratio = 0`` is a pure absorber / streaming medium where the
    prescribed inflow propagates with minimal redistribution — the regime where
    an inflow-representation error is LEAST boundary-confined).  ``bc_scale``
    amplifies the mu-dependent B,C drivers, inflating the boundary
    transverse-slope-to-average ratio (the user's "non-trivial slope at the
    boundary" hypothesis, stressed to the limit).  The MMS source machinery
    (the PDE residual) stays exact for any Sigma_t / Sigma_s by construction.
    """
    def sigma_t_fn(x, y, g):
        x = np.asarray(x, dtype=float)
        return np.full_like(x, sigt)

    def sigma_s_fn(x, y, g_from, g_to):
        x = np.asarray(x, dtype=float)
        if g_from == g_to:
            return np.full_like(x, c_ratio * sigt)
        return np.zeros_like(x)

    case = SN2DCartesianLDStressMMSCase(
        name="s9_streaming",
        length_x=length,
        length_y=length,
        c_spectrum=np.asarray((1.0, 0.4), dtype=float),
        sigma_t_fn=sigma_t_fn,
        sigma_s_fn=sigma_s_fn,
        quadrature=Quadrature.level_symmetric(order),
    )
    if bc_scale == 1.0:
        return case
    return _amplify_boundary_slope(case, bc_scale)


def _amplify_boundary_slope(base_case, bc_scale):
    """Return a case whose mu-dependent (B,C) drivers — the carriers of the
    boundary transverse slope — are scaled by ``bc_scale``.  The MMS stays
    exact because BOTH the manufactured source and the inflow trace read the
    SAME ``_drivers`` (scaling them consistently re-derives a valid source)."""
    class _AmpCase(type(base_case)):
        def _drivers(self, x, y, g):
            (A, dAx, dAy, B, dBx, dBy, C, dCx, dCy) = super()._drivers(x, y, g)
            s = bc_scale
            return (A, dAx, dAy, s * B, s * dBx, s * dBy, s * C, s * dCx, s * dCy)

    # Rebuild as the subclass (frozen dataclass — copy fields).  NOTE: this
    # field list must track ``SN2DCartesianLDStressMMSCase``'s dataclass fields;
    # a new field added there is silently dropped here until mirrored.
    return _AmpCase(
        name=base_case.name,
        length_x=base_case.length_x,
        length_y=base_case.length_y,
        c_spectrum=base_case.c_spectrum,
        sigma_t_fn=base_case.sigma_t_fn,
        sigma_s_fn=base_case.sigma_s_fn,
        quadrature=base_case.quadrature,
    )


def _errors(case, nc, sign):
    """Volume-weighted GLOBAL L2 and boundary-RING L2 of the group-0 scalar-flux
    error vs the MMS exact (``phi_exact == A_0``).  ``sign``: None = average-only
    (slope dropped); +1 = moment-resolved with the projected transverse slope;
    -1 = flipped slope."""
    phi, mesh = _solve_with_boundary_slope(case, nc, slope_sign=sign)
    ref = case.phi_exact(mesh.centers_x, mesh.centers_y, 0)
    num = phi[0]
    err = num - ref
    hx = mesh.edges_x[1] - mesh.edges_x[0]
    hy = mesh.edges_y[1] - mesh.edges_y[0]
    vol = hx * hy
    glob = float(np.sqrt((err ** 2 * vol).sum()))
    ring = np.zeros_like(err, dtype=bool)
    ring[0, :] = ring[-1, :] = ring[:, 0] = ring[:, -1] = True
    edge = float(np.sqrt((err[ring] ** 2 * vol).sum()))
    return glob, edge


def _order(seq):
    return [float(np.log2(seq[i] / seq[i + 1])) for i in range(len(seq) - 1)]


# A generous "second-order band": the measured order must lie in [1.7, 2.4] for
# every region/treatment (the boundary ring is super-convergent ~2.4; the global
# is ~2.0).  Documented sub-floor: mom may DIFFER from flat in magnitude by up to
# this fraction without changing the order — that is the sub-floor signature.
_SUBFLOOR_FRACTION = 0.30  # mom within 30% of flat magnitude (always observed <20%)


# ── The COHERENT-PROMISE gate — first-cell-row is ALREADY O(h^2) ─────────────


@pytest.mark.l1
@pytest.mark.verifies("ld-cartesian-2d")
def test_first_cell_row_already_second_order():
    """#257 S9 COHERENT-PROMISE gate — the FIRST-CELL-ROW (i=0, the cells that
    directly consume the xmin inflow) is ALREADY O(h^2) with the average-only
    inflow: LD is 2nd-order AT the boundary, the no-asterisk guarantee.

    This is the decisive scaling argument behind the coherent promise.  The slope
    lifts the *inflow representation* O(h)->O(h^2) (a real effect), but the
    converged flux at the boundary cell is ALREADY second-order from the average
    moment, so the slope is a sub-floor refinement — there is NO O(h) deficiency
    in the converged flux for it to repair.  Both flat and mom first-row orders
    sit at ~2.0.  ``-O``-safe (``np.testing`` / ``pytest.fail``)."""
    case = _streaming_case(0.1, 0.0, length=1.0)
    ncs = [16, 32, 64, 128]
    rows = {"flat": [], "mom": []}
    for nc in ncs:
        for sign, name in [(None, "flat"), (+1.0, "mom")]:
            phi, mesh = _solve_with_boundary_slope(case, nc, slope_sign=sign)
            ref = case.phi_exact(mesh.centers_x, mesh.centers_y, 0)
            err = phi[0] - ref
            hy = mesh.edges_y[1] - mesh.edges_y[0]
            # 1-D L2 along the first column i=0 (the boundary-cell-row norm).
            rows[name].append(float(np.sqrt((err[0, :] ** 2 * hy).sum())))

    print("\n[first-cell-row (i=0) L2 — flat already O(h^2)?]")
    for name in ("flat", "mom"):
        v = rows[name]
        print(f"  {name}: {[f'{x:.3e}' for x in v]} order={[f'{x:.3f}' for x in _order(v)]}")

    flat_orders = _order(rows["flat"])
    for o in flat_orders:
        if not (o >= 1.85):
            pytest.fail(
                f"first-cell-row flat order {o:.3f} < 1.85 — the average inflow "
                f"would then be O(h)-DEFICIENT at the boundary cell, and the "
                f"slope COULD repair it (the coherent promise would carry an "
                f"asterisk and the verdict would flip); it is in fact O(h^2)"
            )
    mom_orders = _order(rows["mom"])
    for o in mom_orders:
        if not (o >= 1.85):
            pytest.fail(
                f"first-cell-row mom order {o:.3f} < 1.85 (regression)"
            )


# ── The SUB-FLOOR VERDICT PINS — guard the "no value gate" conclusion ────────


@pytest.mark.slow
@pytest.mark.l1
@pytest.mark.parametrize(
    "sigt, c_ratio, label",
    [
        (0.1, 0.0, "thin-streaming"),       # Sigma_t*L = 0.1, pure absorber
        (0.5, 0.0, "moderate-streaming"),   # Sigma_t*L = 0.5, pure absorber
        (1.0, 0.0, "unit-streaming"),       # Sigma_t*L = 1.0, pure absorber
        (2.0, 0.0, "thick-streaming"),      # Sigma_t*L = 2.0, pure absorber
        (0.1, 0.5, "thin-scatter"),         # Sigma_t*L = 0.1, c = 0.5
        (1.0, 0.5, "unit-scatter"),         # Sigma_t*L = 1.0, c = 0.5
    ],
)
def test_optical_sweep_slope_never_beats_floor(sigt, c_ratio, label):
    """#257 S9 verdict pin — across the optical-depth axis the boundary
    transverse-slope moment NEVER beats the average-only drop above the O(h^2)
    floor.  flat == mom == flip at O(h^2) globally, in EVERY regime including the
    streaming limit (where the inflow propagates ballistically and an
    inflow-representation error is least boundary-confined).

    Guards the "no value gate / sub-floor" S9 conclusion: if the LD boundary
    closure ever changes so the slope matters above the floor, this reddens and
    the verdict is revisited (#257 S9 / #263).  ``-O``-safe."""
    case = _streaming_case(sigt, c_ratio, length=1.0)
    ncs = [16, 32, 64]
    series = {}
    for sign, name in [(None, "flat"), (+1.0, "mom"), (-1.0, "flip")]:
        gl, ed = zip(*[_errors(case, nc, sign) for nc in ncs])
        series[name] = {"glob": list(gl), "edge": list(ed)}

    print(f"\n[Sigma_t*L={sigt:.1f}  c={c_ratio}  {label}]")
    for name in ("flat", "mom", "flip"):
        g, e = series[name]["glob"], series[name]["edge"]
        print(
            f"  {name:4s} glob={[f'{v:.4e}' for v in g]} ord={[f'{v:.3f}' for v in _order(g)]}"
            f" | edge={[f'{v:.4e}' for v in e]} ord={[f'{v:.3f}' for v in _order(e)]}"
        )

    # Every treatment converges in the second-order band globally.
    for name in ("flat", "mom", "flip"):
        for o in _order(series[name]["glob"]):
            if not (1.7 <= o <= 2.4):
                pytest.fail(
                    f"{label}/{name} global order {o:.3f} left the 2nd-order band "
                    f"[1.7, 2.4] — the LD discretization regressed"
                )
    # mom does NOT beat flat above the floor: the magnitudes track within the
    # sub-floor band (the slope is a sub-floor perturbation, not a repair).
    g_flat = np.array(series["flat"]["glob"])
    g_mom = np.array(series["mom"]["glob"])
    rel = np.abs(g_mom - g_flat) / g_flat
    if not (rel.max() <= _SUBFLOOR_FRACTION):
        pytest.fail(
            f"{label}: |mom - flat|/flat = {rel.max():.3f} exceeded the documented "
            f"sub-floor band {_SUBFLOOR_FRACTION} — the slope would then be "
            f"above-floor and a value gate WOULD be achievable (revisit the verdict)"
        )


@pytest.mark.slow
@pytest.mark.l1
@pytest.mark.parametrize("bc_scale", [1.0, 5.0, 20.0])
def test_amplified_boundary_slope_still_subfloor(bc_scale):
    """#257 S9 verdict pin — the user's strongest hypothesis (a boundary slope
    large relative to the average, in the streaming limit) does NOT unlock a value
    gate.  Amplifying the mu-dependent (slope-carrying) drivers up to 20x leaves
    flat == mom at O(h^2); the slope makes the converged flux MONOTONICALLY
    WORSE, never better.  The sub-floor wall is FUNDAMENTAL to a boundary-trace
    moment, not regime-specific (#257 S9 / #263).  ``-O``-safe."""
    case = _streaming_case(0.1, 0.0, length=1.0, bc_scale=bc_scale)
    ncs = [16, 32, 64]
    gl_flat = [_errors(case, nc, None)[0] for nc in ncs]
    gl_mom = [_errors(case, nc, +1.0)[0] for nc in ncs]

    print(f"\n[bc_scale={bc_scale}  Sigma_t*L=0.1 streaming]")
    print(f"  flat glob={[f'{v:.4e}' for v in gl_flat]} ord={[f'{v:.3f}' for v in _order(gl_flat)]}")
    print(f"  mom  glob={[f'{v:.4e}' for v in gl_mom]} ord={[f'{v:.3f}' for v in _order(gl_mom)]}")
    improves = [m < f for m, f in zip(gl_mom, gl_flat)]
    print(f"  mom improves on flat? {improves} (expected all False — sub-floor)")

    for o in _order(gl_flat):
        if not (1.7 <= o <= 2.4):
            pytest.fail(f"flat order {o:.3f} regressed at bc_scale={bc_scale}")
    for o in _order(gl_mom):
        if not (1.7 <= o <= 2.4):
            pytest.fail(f"mom order {o:.3f} regressed at bc_scale={bc_scale}")
    # The decisive claim: mom never beats flat (no value gate possible).
    if any(improves):
        pytest.fail(
            f"bc_scale={bc_scale}: mom beat flat in the converged flux — the value "
            f"gate WOULD be achievable; the S9 verdict must be revisited"
        )


# ── The Mode-11 sentinel — the toggle reaches the production consumer ────────


@pytest.mark.foundation
def test_slope_toggle_reaches_inflow_to_moments():
    """#257 S9 Mode-11 guard: the flat/mom/flip toggle ACTUALLY changes the
    boundary moment slot-1 that the production LD closure consumes (proves the
    sub-floor null result is real, not a vacuous toggle).  Sentinel-instruments
    the production ``_LossRepresentation._inflow_to_moments`` and confirms:
    flat/zero-ctrl see slot-1 == 0, mom/flip see slot-1 != 0, and the converged
    flux DIFFERS (the slope is genuinely consumed).  ``-O``-safe."""
    case = build_2d_cartesian_ld_stress_mms_case()
    orig = LR._LossRepresentation._inflow_to_moments
    seen = {"slot1_max": 0.0, "n": 0}

    def spy(self, inflow):
        out = orig(self, inflow)
        for face in out:
            if face.ndim >= 1 and face.shape[-1] == 2:
                seen["slot1_max"] = max(seen["slot1_max"], float(np.abs(face[..., 1]).max()))
                seen["n"] += 1
        return out

    results = {}
    try:
        LR._LossRepresentation._inflow_to_moments = spy
        for sign, name in [(None, "flat"), (+1.0, "mom"), (-1.0, "flip"), (0.0, "zero")]:
            seen["slot1_max"] = 0.0
            seen["n"] = 0
            phi, _ = _solve_with_boundary_slope(case, 8, slope_sign=sign)
            results[name] = (seen["slot1_max"], seen["n"], float(phi.sum()))
    finally:
        LR._LossRepresentation._inflow_to_moments = orig

    print("\n[Mode-11 sentinel: slot-1 reaching production _inflow_to_moments]")
    for name, (s1, n, ps) in results.items():
        print(f"  {name:5s} slot1_max={s1:.4e} n_calls={n} phi_sum={ps:.8e}")

    # The production consumer is reached on every solve.
    for name in results:
        if not (results[name][1] > 0):
            pytest.fail(f"{name}: _inflow_to_moments never saw a moment slot")
    # flat and zero-ctrl drop the slope (slot-1 == 0); mom/flip carry it.
    np.testing.assert_array_equal(results["flat"][0], 0.0)
    np.testing.assert_array_equal(results["zero"][0], 0.0)
    if not (results["mom"][0] > 1e-4):
        pytest.fail("mom: slope slot-1 did not reach the consumer")
    if not (results["flip"][0] > 1e-4):
        pytest.fail("flip: slope slot-1 did not reach the consumer")
    # The slope is genuinely consumed — the converged flux differs.
    if not (abs(results["mom"][2] - results["flat"][2]) > 1e-3):
        pytest.fail(
            "mom converged flux == flat — the slope was carried but NOT consumed "
            "(the toggle would be vacuous and the null result meaningless)"
        )
    # zero-ctrl is byte-identical to flat (the Leg-B asymmetry / no-op control).
    np.testing.assert_array_equal(results["zero"][2], results["flat"][2])
