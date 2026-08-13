r"""``D`` — the over-determination residual of the per-level angular march.

A curvilinear :math:`\mu`-level's angular march has **two** endpoints, not
one.  The redistribution coefficient vanishes at each
(:math:`\alpha_{1/2} = \alpha_{M+1/2} = 0`), so the cell balance decouples
there into a plain radial DD ODE, and production solves **both** with one
engine — ``carlson_inward_sweep_from_source``, called twice from
``RadialCharacteristicOperator.solve`` (the inward leg becomes the M-M
recurrence's seed, the outward leg is stored as ``cells(p, +1)``).

The recurrence, marched from the seed across all :math:`M` ordinates, then
predicts the far end a *second* time, as ``_MMHalfGrid.trailing_face``.  Only
one of the two is imposed.  ``D`` is what the unimposed one is violated by:

.. math::

   D_p \;:=\; \psi_{M+1/2}\big|_p \;-\; \psi^{\text{marched}}_p(+1)

Two structurally different discretizations of ONE point of phase space,
computed on every curvilinear solve since the ψ½ route landed and — until
these gates — compared by nothing.

⛔ WHAT ``D`` IS NOT
===================

**It is not an error estimator, and it may not vote on** :math:`\tau`.  This
is a measurement, not a caution.  `[M]` 2026-08-12, against the **analytic**
anisotropic cylindrical MMS (``build_cylindrical_anisotropic_mms_case`` —
an analytic reference, not another ORPHEUS solver, so not the L49
reference-limited trap; nx = 80, ``inner_tol = 1e-13``, all 12 solves
converged), the Pearson correlation of :math:`\log D` against :math:`\log`
of the true MMS error across four :math:`\tau` variants:

=============  ==============  ==========================
n_φ            ranks agreeing  Pearson r on log(metric)
=============  ==============  ==========================
8              2 / 4           ``+0.7515``
16             0 / 4           ``+0.2608``
32             0 / 4           ``+0.0630``
=============  ==============  ==========================

The correlation **degrades monotonically to zero** as angle refines.
Structurally it must: :math:`D = e_1 - e_2` is a DIFFERENCE of two truncation
errors, hence small exactly when both are large and equal.

``D`` ranks the shipped Q5.6.4 partition first by 2.6–45× over garbage
:math:`\tau` — and that ranking is **not** evidence for the partition,
because the instrument is uncorrelated with accuracy.  The campaign still
has NO reference-free instrument that can rank :math:`\tau`
(``scratch/q65_endpoint_defect_findings.md`` §F7/§F10, and the memo's
closing section).  Any future ``D``-based :math:`\tau` argument must cite
that first.

⛔ **``D`` must also not be used to CORRECT the seed.**  The march's linear
part is exactly :math:`(-1)^M I` — it follows from
:math:`\prod_m (1-\tau_m)/\tau_m = 1`, which is gated on both arms in
``test_psi_half_positivity.py`` — and BOTH endpoint values come from
physics, so imposing both is an over-determination (a constraint on the
interior solution), not an equation for a free parameter.

WHAT THESE GATES PIN
====================

Only what ``D`` legitimately is: a cheap, reference-free, pointwise
**consistency** residual, whose behaviour under refinement is a property of
the scheme.

FIXTURE DISCIPLINE — the blindness is structural
================================================

Every row here uses a **heterogeneous, two-group, VACUUM** curvilinear
problem (`vv-principles` #3/#4).  Not decoration:

* a **flat-in-angle** fixture is *provably* blind — the recurrence's
  flat-flux fixed point makes both endpoints coincide exactly, which
  :func:`test_the_defect_is_exactly_zero_when_the_endpoints_agree` pins as
  the positive control (`vv-principles` #11) *and* as the reason no flat
  gate anywhere in the tree could ever have caught this;
* a **1-group / homogeneous** fixture nulls the redistribution terms the
  march is made of.

⚠ A **reflective** outer face is NOT blind — `[M]` §F6 measures ``max|D|``
at ``9.2e-2 … 1.6e-1`` under reflection, comparable to vacuum.  What
reflection removes is the *divergence* of ``L∞(D)`` (issue #360), not the
defect.

MUTATION LEDGER — `[M]` 2026-08-13, and what each gate CANNOT see
=================================================================

In-process monkeypatch only (never a git-level revert — the tree carries
uncommitted state).  ``M3`` is the **positive control** (`vv-principles`
#17): a battery with no control cannot tell "nothing to catch" from "the
harness is dead", and it reports the reassuring one.

===========================================  ==  ==  ==
mutation                                     M1  M2  M3
===========================================  ==  ==  ==
``trailing_face -> faces[:, -2, :]``         ✓
``D vs cells(p, -1)``, i.e. the SEED             ✓
``trailing_face -> zeros`` (CONTROL)                 ✓
===========================================  ==  ==  ==

==============================================  ===  ===  ===
gate                                            M1   M2   M3
==============================================  ===  ===  ===
``…trailing_face_is_the_slice…``                RED   .   RED
``…defect_is_the_trailing_face_minus…``          .   RED   .
``…defect_is_exactly_zero…``                     .    .   RED
``…cylinder…falls_second_order…``               RED  RED  RED
``…sphere…turnover_moves_out…``                 RED  RED  RED
``…sees_the_ANGULAR_PARTITION…``                RED  RED  RED
**total reddened**                              **4** **4** **5**
==============================================  ===  ===  ===

⚠ **The two survivals are structural, not gaps — and saying so is the point**
(`vv-principles` #19/#20: a row that cannot see a property must not be
counted as covering it).

* ``…defect_is_the_trailing_face_minus…`` survives M1 and M3 because **both
  of its sides read** ``trailing_face``.  It is a tautology with respect to
  *which slice* the property returns, and can only ever test the
  *subtraction*.  The slice itself is covered by
  ``…trailing_face_is_the_slice…``, which reddens for exactly those two.
  Neither gate is redundant; neither is sufficient alone.
* ``…defect_is_exactly_zero…`` survives M1 and M2 because a flat flux
  collapses **every** face to the same constant — so no face-index error and
  no seed-vs-far-endpoint confusion is observable through it.  That is the
  same blindness it exists to document, now measured rather than argued.

NORM CONVENTION
===============

``_rms`` below is the plain RMS over every ``(level, group, cell)`` entry.
⚠ The findings memo's "L2" figures are :math:`\sqrt{n_g}` LARGER — it
normalised per ``(level, cell)`` while summing over groups, and never said
so (`plan-authoring` §4: a number without its convention is not reusable).
`[M]` the two agree to machine precision after that factor; ``max|D|``, the
signed mean, and every RATIO are convention-free and match the memo's
digits exactly.
"""
from __future__ import annotations

from typing import Any, NamedTuple

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn import solve_sn_fixed_source
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.sweep.pole_angular_closure import MorelMontryAngularSweep

pytestmark = pytest.mark.foundation


# ═══════════════════════════════════════════════════════════════════════
# Fixtures — heterogeneous, 2-group, vacuum, both curvilinear arms
# ═══════════════════════════════════════════════════════════════════════


def _two_region_2g(coord: CoordSystem, quad: Quadrature, *, nx: int) -> SNMesh:
    """Fuel-like A inside r < 1, moderator-like B outside, vacuum at R = 2."""
    edges = np.linspace(0.0, 2.0, nx + 1)
    r_mid = 0.5 * (edges[:-1] + edges[1:])
    mesh = Mesh1D(
        edges=edges,
        mat_ids=np.where(r_mid < 1.0, 0, 1).astype(int),
        coord=coord,
        bc_right=BC("vacuum"),
    )
    materials = {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")}
    return SNMesh(mesh, quad, materials)


def _cylinder(n_phi: int, *, n_mu: int = 4, nx: int = 12) -> SNMesh:
    return _two_region_2g(
        CoordSystem.CYLINDRICAL,
        Quadrature.folded_product(n_mu=n_mu, n_phi=n_phi),
        nx=nx,
    )


def _sphere(n_ordinates: int, *, nx: int) -> SNMesh:
    return _two_region_2g(
        CoordSystem.SPHERICAL,
        Quadrature.gauss_legendre(n_ordinates),
        nx=nx,
    )


class _Marched(NamedTuple):
    """Everything the endpoint defect needs from one converged solve."""

    closure: MorelMontryAngularSweep
    psi: np.ndarray
    interior: Any                       # RadialCharacteristicInteriorField
    state: tuple


def _march(sn: SNMesh) -> _Marched:
    """Solve, then drive the PRODUCTION recurrence — the one solve site.

    Fails loudly on a starved solve (#340) and on a seedless mesh: ``D``'s
    magnitude has exactly the signature a starved inner would fake, so a
    silent degradation here would be read as a scheme property.
    """
    mesh = sn.mesh
    if mesh is None:
        pytest.fail("the fixture built no Mesh1D — nothing to solve")
    solution = solve_sn_fixed_source(
        sn.materials, mesh, sn.quad,
        external_source=np.ones((sn.quad.N, 2, sn.nx)),
        inner_tol=1e-12,
    )
    if not solution.converged():
        pytest.fail(
            "the fixed-source solve did not fully converge, so any D "
            "measured from it is an iteration artefact, not a property of "
            "the angular march (#340)"
        )
    rc = solution.radial_characteristic
    if rc is None:
        pytest.fail(
            "no radial_characteristic member: this mesh marched no psi-half "
            "state, so there is no SECOND endpoint and D is undefined rather "
            "than zero — the fixture is not exercising the two-ended march"
        )
    closure = sn.pole_angular_closure
    # The defect is a property of a MARCHING closure. The Cartesian identity
    # closure has no angular march and therefore no endpoints at all, so the
    # method lives on M-M alone (never a base-class stub) and consumers narrow.
    assert isinstance(closure, MorelMontryAngularSweep)
    psi = np.asarray(solution.angular_flux.interior.values, dtype=float)
    state = closure.precompute_psi_state(psi, radial_characteristic=rc.interior)
    return _Marched(closure, psi, rc.interior, state)


def _solved_defect(sn: SNMesh) -> dict[int, np.ndarray]:
    """``{level: D_p}`` from a CONVERGED fixed-source solve."""
    m = _march(sn)
    return m.closure.angular_endpoint_defect_per_level(m.state, m.interior)


def _rms(defect: dict[int, np.ndarray]) -> float:
    """RMS over every ``(level, group, cell)`` entry — see the module norm note."""
    return float(
        np.sqrt(np.mean(np.concatenate([d.ravel() for d in defect.values()]) ** 2))
    )


# ═══════════════════════════════════════════════════════════════════════
# The CONTRACT — what the trailing face is, and what D subtracts
# ═══════════════════════════════════════════════════════════════════════


def test_the_trailing_face_is_the_slice_upstream_per_ordinate_drops():
    r"""⭐ ``trailing_face`` is exactly the one face no consumer reads.

    The off-by-one contract :class:`_MMHalfGrid` exists to enforce, stated
    from the other side: ``upstream_per_ordinate`` is ``faces[:, :-1, :]``
    and ``trailing_face`` is ``faces[:, -1, :]``, so together they partition
    the grid with no overlap and nothing left over.

    Solve-free, and it is the gate that would redden if a future carve
    re-indexed either accessor.
    """
    sn = _cylinder(8)
    _closure, _psi, _interior, state = _march(sn)

    for grid in state:
        m_total = grid.n_ordinates
        np.testing.assert_array_equal(grid.trailing_face, grid.faces[:, m_total, :])
        np.testing.assert_array_equal(grid.trailing_face, grid.faces[:, -1, :])
        # Partition: upstream slice + trailing face reconstruct the grid.
        assert grid.upstream_per_ordinate.shape[1] == m_total
        np.testing.assert_array_equal(
            np.concatenate(
                [grid.upstream_per_ordinate, grid.trailing_face[:, None, :]], axis=1,
            ),
            grid.faces,
        )


@pytest.mark.verifies("sn-angular-endpoint-defect-eq")
def test_the_defect_is_the_trailing_face_minus_the_marched_endpoint():
    r"""``D_p = trailing_face - cells(p, +1)``, over the CARRYING levels only.

    Pins the definition against a hand-written difference, and pins the
    key set: a non-carrying level seeds from ``edge_extrapolated_seed``,
    which is an *interpolation* and not a solved endpoint, so no comparison
    is defined for it and it must not appear.

    ⚠ **This row is a TAUTOLOGY with respect to which slice**
    ``trailing_face`` **returns** — both of its sides read that property, so
    re-indexing it moves them together.  `[M]` it stays green under both the
    off-by-one and the zeroing mutation (module ledger, M1/M3).  It tests
    the SUBTRACTION and the key set, nothing more; the slice is covered by
    :func:`test_the_trailing_face_is_the_slice_upstream_per_ordinate_drops`.
    """
    sn = _cylinder(8)
    closure, _psi, interior, state = _march(sn)
    defect = closure.angular_endpoint_defect_per_level(state, interior)

    carrying = sorted(closure._carrying_levels)
    assert sorted(defect) == carrying
    # A folded cylinder carries on EVERY level since Q5.6.3 — if this ever
    # reads fewer, the fixture stopped exercising the two-ended march.
    assert len(carrying) == sn.quad.n_levels

    for p in carrying:
        np.testing.assert_array_equal(
            defect[p],
            state[p].trailing_face - interior.cells(p, +1),
        )
        assert defect[p].shape == (2, sn.nx)


def test_the_defect_is_exactly_zero_when_the_endpoints_agree():
    r"""⭐ POSITIVE CONTROL — and the proof that flat fixtures are blind.

    `vv-principles` #11: a residual gate tested only against configurations
    that produce a non-zero reading validates the arithmetic, not the claim.
    Here the zero is reachable exactly, with no tolerance argument.

    On a flux that is CONSTANT in angle, seeded with that same constant, the
    M-M recurrence has :math:`\psi` as an exact fixed point:

    .. math::

       \psi_{m+1/2} = \frac{\psi_m - (1-\tau_m)\psi_{m-1/2}}{\tau_m}
                    = \frac{c - (1-\tau_m)c}{\tau_m} = c

    for **every** :math:`\tau_m`.  So the trailing face is :math:`c`, and if
    the marched outward endpoint is also :math:`c` the defect vanishes.

    ⚠ **Algebraically exact, not bit-exact.**  Floating point makes each of
    the :math:`M` steps cost about a ULP (one subtraction, one division),
    and the map is neutrally stable (:math:`|{\rm gain}| = 1`), so the
    round-off neither grows nor cancels.  The bound below is therefore
    derived — 4 ULP of the constant per step — and not fitted.  `[M]`
    2026-08-13 the worst entry is ``5.551e-17``, i.e. **0.67 ULP** of
    ``0.375`` and ~24× inside the bound, while a REAL endpoint disagreement
    on this fixture is :math:`O(10^{-1})` — fourteen orders larger, so the
    bound has no chance of hiding one.

    ⟹ This is precisely why no flat-flux gate in the tree could ever have
    caught the endpoint disagreement, and why every measuring row in this
    module is heterogeneous.

    ⚠ **What this row cannot see, measured rather than argued.** A flat flux
    collapses every face to the same constant, so it is blind to
    :math:`\tau` (the fixed point holds for all of them), blind to a
    face-INDEX error, and blind to reading the near endpoint in place of the
    far one — `[M]` it stays green under both M1 and M2 (module ledger).
    Do not read it as coverage of the partition, of the slice, or of the
    operand choice.  It covers exactly one thing: that the defect can reach
    zero at all, which is what makes a non-zero reading elsewhere meaningful.
    """
    sn = _cylinder(8)
    closure, _psi, interior, _state = _march(sn)

    constant = 0.375
    psi_flat = np.full(np.asarray(_psi).shape, constant, dtype=float)
    for p in sorted(closure._carrying_levels):
        interior.cells(p, -1)[...] = constant   # the seed
        interior.cells(p, +1)[...] = constant   # the marched far endpoint

    state = closure.precompute_psi_state(psi_flat, radial_characteristic=interior)
    defect = closure.angular_endpoint_defect_per_level(state, interior)

    assert defect, "no carrying level — the control proves nothing"
    eps = float(np.finfo(float).eps)
    for p, d_p in defect.items():
        # 4 ULP of the constant per recurrence step — derived, not fitted.
        bound = 4.0 * state[p].n_ordinates * eps * constant
        worst = float(np.abs(d_p).max())
        assert worst <= bound, (
            f"level {p}: flat flux with a consistent seed left D = {worst:.4e}, "
            f"above the {bound:.4e} round-off bound. The recurrence's flat-flux "
            f"fixed point holds for EVERY tau, so this is not a partition "
            f"question — either the seed stopped reaching the grid or the "
            f"marched endpoint is no longer read from the same state."
        )


# ═══════════════════════════════════════════════════════════════════════
# The CHARACTERIZATION — how D behaves under refinement
# ═══════════════════════════════════════════════════════════════════════


def test_cylinder_endpoint_defect_falls_second_order_in_angle():
    r"""⭐ ``D`` converges at ~2nd order in :math:`n_\varphi` on the cylinder.

    `[M]` 2026-08-13, ``folded_product(n_mu=4, n_phi=·)``, 2-region 2-group,
    ``R_out = 2``, **vacuum**, ``nx = 12``, ``external_source = 1``,
    ``inner_tol = 1e-12``, every solve ``converged()`` — RMS over all
    ``(level, group, cell)``:

    =====  ============  =======  =================
    n_φ    RMS(D)        ratio    observed order
    =====  ============  =======  =================
    8      ``5.8276e-2``  —       —
    16     ``1.3109e-2``  4.45×   ``2.15``
    32     ``3.4039e-3``  3.85×   ``1.95``
    =====  ============  =======  =================

    ⚠ **The ORDER is not a property of the shipped partition.**  A two-point
    slope once read as one — from {8, 16} alone the shipped :math:`\tau`
    gave 2.15 against a diamond :math:`\tau \equiv \tfrac12`'s 1.53, which
    looks like a recovered order.  With :math:`n_\varphi = 32` the diamond
    reaches 2.15 and a *shuffled* :math:`\tau` reaches 2.19: `[M]` **all
    four variants are asymptotically ~2nd order**.  What the shipped
    partition buys is a smaller error CONSTANT and reaching the rate at the
    coarsest rung tested — recorded in
    ``scratch/q65_endpoint_defect_findings.md`` §F7, refuted in place there.

    The bound below is deliberately loose on the low side (this is a
    characterization, not an accuracy claim) and tight enough that a
    first-order or non-convergent regression reddens it.
    """
    rms = {n_phi: _rms(_solved_defect(_cylinder(n_phi))) for n_phi in (8, 16, 32)}

    orders = [
        float(np.log2(rms[coarse] / rms[fine]))
        for coarse, fine in ((8, 16), (16, 32))
    ]
    assert all(o > 1.7 for o in orders), (
        f"D stopped converging at ~2nd order in angle: RMS {rms}, orders {orders}"
    )
    assert all(o < 2.6 for o in orders), (
        f"observed order above the plausible band — suspect the fixture, not "
        f"a scheme improvement: RMS {rms}, orders {orders}"
    )
    # The ladder itself, loosely pinned so a 2x drift in the CONSTANT is audible.
    for n_phi, expected in ((8, 5.8276e-2), (16, 1.3109e-2), (32, 3.4039e-3)):
        assert 0.5 * expected < rms[n_phi] < 2.0 * expected, (
            f"n_phi={n_phi}: RMS(D) = {rms[n_phi]:.4e} against the recorded "
            f"{expected:.4e} — the error CONSTANT moved even if the order held"
        )


def test_sphere_endpoint_defect_turnover_moves_out_with_the_spatial_grid():
    r"""⭐⭐ The sphere's ``D`` DOES converge in angle — the floor is SPATIAL.

    This gate exists because the opposite was believed, measured, and
    reported — twice, independently — before being refuted.

    `[M]` at ``nx = 12`` the sphere's ``D`` falls from N=4 to N=8 and then
    **rises** (``7.4973e-2 → 3.0749e-2 → 3.4389e-2 → 4.4163e-2`` at
    N = 4/8/16/32).  Read alone that says the angular march does not
    converge in angular order.  It is an artefact: at ``nx = 12`` the
    SPATIAL error is already the floor by N=8.

    `[M]` at ``nx = 48`` the same ladder falls monotonically
    (``7.5456e-2 → 2.1956e-2 → 1.0125e-2 → 8.2441e-3``, ratios
    ``3.44 / 2.17 / 1.23``) — the turnover has moved OUT past N=32.

    ⟹ **The turnover point moving out with** ``nx`` **is the textbook
    signature of a spatially-set floor**, and that is what this gate pins —
    the mechanism, not a ladder.  A regression that genuinely broke angular
    convergence would break the nx=48 monotonicity; a regression that only
    coarsened the spatial arm would move the turnover back IN.

    (``scratch/q65_endpoint_defect_findings.md`` §F8 records the refuted
    reading in place, and §F8b separately kills the pole-consistency
    hypothesis that was offered to explain it.)
    """
    coarse = [_rms(_solved_defect(_sphere(n, nx=12))) for n in (4, 8, 16, 32)]
    fine = [_rms(_solved_defect(_sphere(n, nx=48))) for n in (4, 8, 16, 32)]

    # nx = 48: angular convergence is real and monotone across the ladder.
    assert all(b < a for a, b in zip(fine, fine[1:])), (
        f"the sphere's endpoint defect stopped falling with angular order at "
        f"nx=48, where no spatial floor should yet bite: {fine}"
    )
    # nx = 12: the spatial floor bites, and the defect turns back up.
    assert coarse[-1] > coarse[1], (
        f"nx=12 no longer shows the spatial floor turning D back up; the "
        f"refuted-reading record in F8 no longer reproduces: {coarse}"
    )
    # The floor is SPATIAL: refining space lowers the fine end, not the coarse.
    assert fine[-1] < 0.5 * coarse[-1], (
        f"refining nx did not lower D at the finest angle, so the N=32 floor "
        f"is not spatial after all: nx12={coarse[-1]:.4e} nx48={fine[-1]:.4e}"
    )


def test_the_endpoint_defect_sees_the_ANGULAR_PARTITION_not_just_its_product():
    r"""⭐ ``D`` is loaded on :math:`\tau` — the per-step partition, not the product.

    The tree's neutral-stability gates pin
    :math:`\prod_m (1-\tau_m)/\tau_m = 1`
    (``test_psi_half_positivity.py``), which is a statement about the
    PRODUCT.  A permutation of the shipped :math:`\tau` multiset preserves
    that product **exactly**, so those gates are GREEN for it — they are
    structurally incapable of seeing the per-step partition.

    ``D`` sees it.  This row is the teeth: an in-class mutation
    (`vv-principles` #18 — reversing a permutation keeps the recurrence
    linear, keeps the product identity, and changes only the thing under
    test) must move ``D`` substantially.

    ⛔ **This is a SENSITIVITY claim and nothing more.**  It does NOT say the
    shipped order is better: `[M]` ``D`` is uncorrelated with accuracy
    (module docstring), so the direction of the change carries no verdict.
    The gate asserts only that the instrument is not :math:`\tau`-blind, and
    it deliberately does not assert a sign.
    """
    sn = _cylinder(16)
    closure, psi, interior, _state = _march(sn)

    shipped_tau = closure._tau_per_level
    reversed_tau = tuple(t[::-1].copy() for t in shipped_tau)

    # The product identity is preserved EXACTLY by the permutation — so the
    # existing neutral-stability gates cannot distinguish these two.
    for a, b in zip(shipped_tau, reversed_tau):
        np.testing.assert_allclose(
            np.prod((1.0 - a) / a), np.prod((1.0 - b) / b), rtol=0.0, atol=1e-12,
        )

    state = closure.precompute_psi_state(psi, radial_characteristic=interior)
    shipped = _rms(closure.angular_endpoint_defect_per_level(state, interior))

    object.__setattr__(closure, "_tau_per_level", reversed_tau)
    try:
        state_mut = closure.precompute_psi_state(psi, radial_characteristic=interior)
        mutated = _rms(
            closure.angular_endpoint_defect_per_level(state_mut, interior)
        )
    finally:
        object.__setattr__(closure, "_tau_per_level", shipped_tau)

    assert mutated > 3.0 * shipped, (
        f"D barely moved under a tau permutation that the product-identity "
        f"gates cannot see ({shipped:.4e} -> {mutated:.4e}); the instrument "
        f"has gone tau-blind and no longer adds anything those gates lack"
    )
