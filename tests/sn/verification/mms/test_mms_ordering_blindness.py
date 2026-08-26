r"""The curvilinear MMS suite CANNOT adjudicate the per-level ordinate
ordering — this file is that Mode-7 declaration made executable.

Promoted 2026-08-01 from
``derivations/diagnostics/diag_326_mms_ordering_blindness.py``
(numerics-investigator, GitHub issue **#326**).

WHY IT LIVES HERE
-----------------
``test_mms_curvilinear.py`` and ``test_curvilinear_aniso_convergence.py`` are
the gates a fresh session will reach for when asked "does the curvilinear MMS
see X?".  For X = *the per-level ordinate ordering* the answer is **no, and
exactly no** — but nothing in those files says so, and their Mode-7 honest-scope
note exists only as prose.  This module is the executable form of it, so the
claim cannot decay silently.

THE SITUATION
-------------
``rules_product.py`` orders each mu-level by ``np.argsort(mu_x, kind="stable")``.
``eta = mu_x = sin(theta) cos(phi)`` is 2-to-1 over ``phi`` in ``[0, 2*pi)``:
the azimuthal mirror pair ``(phi, 2*pi - phi)`` shares it.  The level was never
totally ordered, and the tie-break moves a heterogeneous cylindrical solve by
0.6 %-7.2 %.

UPDATED 2026-08-02 — the tie-break used to be decided BY ROUNDING NOISE.  The
azimuths are now roots of unity (#325's generator half, landed with the
exactness carve), so the ties are exact and the tie-break is a named rule:
``kind="stable"``, eta ascending with ties broken by increasing phi.  ``[M]``
What that substitution actually moved, isotropic cylindrical MMS at
``nx=40, n_mu=4, n_phi=8``::

    MMS L2 error   5.389276744519e-04 -> 5.389276744525e-04   (1.2e-12)
    scalar flux    max rel diff 6.7e-16
    per-ordinate   max rel diff 3.0e-05

i.e. machine-level in the scalar flux, solver-tolerance in the ladder, and
real but small per-ordinate — exactly the sector this file proves the MMS is
blind to, and in which ``test_azimuthal_mirror_symmetry.py`` measures the
defect magnitude to be ordering-INVARIANT.  So the substitution did not make
the solver more or less correct; it replaced a noise-decided label with a
named one.

WHAT THIS FILE PROVES
---------------------
1. ``alpha`` and ``tau`` are BIT-IDENTICAL across tie-breaks.  ``alpha`` is a
   partial sum of ``w_m * eta_m`` and ``w`` is constant within a product-rule
   level, so permuting equal ``eta`` cannot move any partial sum; ``tau`` reads
   only the sorted ``eta`` sequence.  The tie-break moves ONLY which ordinate
   sits at which position — i.e. which member of each mirror pair receives
   ``tau = 1`` and which the clamped ``tau = 1/2``.

2. **BOTH production cylindrical MMS ansatzes are functions of ``eta_n`` and
   ``xi_n**2`` ONLY** (``mms/sn.py``: ``psi_n = (A + B*eta_n)/W``;
   ``Q_n = eta A' + eta**2 B' + xi**2 B/r + (SigT-SigS) A + SigT eta B``).  The
   mirror pair shares ``eta`` — that IS the sort-key degeneracy — and shares
   ``xi**2``.  So it carries identical ``psi_ref``, identical ``Q``, identical
   ``w``: swapping it is a PURE RELABELING of two ordinates carrying identical
   data, and ``phi = sum_n w_n psi_n`` is invariant.

   **Therefore the curvilinear MMS is structurally blind to the entire
   ``xi``-ODD sector**, of which the ordering degeneracy is one instance.  That
   is ``vv-principles`` Mode 12 in exact form — the measured functional's
   invariance group CONTAINS the mutation class.  No tolerance, no mesh
   refinement, no quadrature order can expose it through these ansatzes.

3. A ``xi``-ODD companion ansatz ``psi_n = (A(r) + B(r) xi_n)/W`` leaves the
   symmetric sector and DOES see the tie-break (20.6 % apart at nx=10) — but it
   does NOT adjudicate: both orderings converge to the SAME angular floor from
   opposite sides (0.28 % apart at nx=80) at zero spatial order.  There is no
   "right" ordering to find; the closure itself is what is under-determined.

WHAT THIS FILE IS BLIND TO
--------------------------
It verifies no transport equation and carries no ``verifies(...)`` edge.  Every
row is a claim about the TEST DESIGN (an algebraic reduction invariant, a
property of an MMS case object, a blindness measurement), which is what
``foundation`` means.  The criterion that DOES bite the underlying defect is
the geometric ``xi``-mirror symmetry —
``tests/sn/sweep/curvilinear/test_azimuthal_mirror_symmetry.py``.

MEASURED 2026-08-01 (branch ``refactor/operator-strategy-layers``)
------------------------------------------------------------------
Isotropic cylindrical MMS, nx [20,40,80,160], the production gate's ladder::

    lexsort  [2.160863081428e-03 5.389276744521e-04 1.346506256551e-04 3.365754646309e-05]
    stable   [2.160863081428e-03 5.389276744522e-04 1.346506256550e-04 3.365754646321e-05]
    max rel diff 3.4e-12, orders identical [2.0034 2.0009 2.0002]

Anisotropic cylindrical MMS, nx [10,20,40,80]::

    both orderings  [0.022096038325 0.019504300884 0.019116118385 0.019037806484]
    max rel diff 8.6e-15
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
import pytest

from orpheus.derivations.continuous.mms.sn import (
    _make_1g_mixture,
    build_cylindrical_mms_case,
    build_cylindrical_anisotropic_mms_case,
)
from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.geometry.boundary import SelfPairedDeck
from orpheus.numerics.quadrature import (
    STAGGERED,
    Quadrature,
    gauss_legendre_on_mu,
    periodic_trapezoid,
    spherical_product,
)
from orpheus.sn import solve_sn_fixed_source


def _xi_mirror_pairing(quad: Quadrature) -> np.ndarray:
    """The σ_y (ξ-mirror) ordinate pairing, from production's own source."""
    pi = quad.ordinate_permutation(SelfPairedDeck.mirror(axis="y").motion)
    assert pi is not None, "every fixture here is σ_y-closed (k = 0 plane)"
    return pi.indices

from tests.sn._test_helpers import product_level_ordering, volume_weighted_l2

# Every claim here is about test design / quadrature algebra, not about a
# transport equation — so the whole module is `foundation` and NOTHING carries
# a `verifies(...)` edge.
pytestmark = pytest.mark.foundation


# ── 0. the swap is NON-VACUOUS (control legs) ───────────────────────────

def test_the_two_tie_breaks_really_order_the_level_differently():
    """Control: without this the invariance results below prove nothing."""
    got = {}
    for tie_break in ("lexsort", "stable"):
        with product_level_ordering(tie_break):
            got[tie_break] = [
                list(map(int, a))
                for a in Quadrature.product(n_mu=4, n_phi=8).level_indices
            ]
    print(f"level0 lexsort={got['lexsort'][0]}  stable={got['stable'][0]}")
    assert got["lexsort"] != got["stable"], (
        f"the two tie-breaks produced the SAME ordering {got} — the "
        f"blindness measurements below would be vacuous"
    )


def test_trig_nodes_hide_the_tie_break_behind_rounding_noise():
    """#325 <-> #326 coupling — HALF OF THIS RESOLVED 2026-08-02.

    With ``np.cos(np.linspace(0, 2pi, n_phi, endpoint=False))`` the mirror
    pair's ``eta`` differ by ~1 ULP, so EVERY sort rule agrees and the
    degeneracy is invisible.  That half is unchanged and is asserted below
    against the ``exact_nodes=False`` arm, which is now a HISTORICAL control
    rather than a description of production.

    Production moved: ``rules_product`` builds its azimuths with
    ``periodic_trapezoid`` (roots of unity, #325's generator half), so the
    ties are now EXACT and the tie-break IS the decision.  The companion
    ``test_production_now_has_exact_ties_and_a_NAMED_tie_break`` is the
    other half of this pair.

    On "#326 blocks #325": the blocking condition was that exact nodes hand
    the ordering decision to a tie-break rule, so one has to be CHOSEN
    rather than inherited from noise.  It was — ``kind="stable"``, i.e. eta
    ascending with ties broken by increasing phi — so the condition is
    discharged, not violated.  Choosing is not the same as *ranking*: this
    file's own headline result stands, that no ordering is more correct
    than another, and Q5.3 may yet re-pose it as the azimuthal march.
    """
    got = {}
    for tie_break in ("lexsort", "stable"):
        with product_level_ordering(tie_break, exact_nodes=False):
            got[tie_break] = [
                list(map(int, a))
                for a in Quadrature.product(n_mu=4, n_phi=8).level_indices
            ]
    assert got["lexsort"] == got["stable"], (
        f"trig-evaluated nodes unexpectedly left a REACHABLE tie: {got}"
    )


def test_production_now_has_exact_ties_and_a_NAMED_tie_break():
    r"""The inversion of the row above, and the reason it is a pair.

    ``eta = sin(theta) cos(phi)`` and ``cos(phi_m) == cos(phi_{n-m})``
    bit-exactly for roots of unity, so a production level carries only
    ``n_phi//2 + 1`` distinct eta among ``n_phi`` ordinates and the
    interior gaps are EXACTLY zero — where they used to be ~1 ULP of
    rounding noise standing in for a tie.

    Keeping both rows is the point: the historical claim is what makes
    the present one legible, and a lone "gaps are zero" assertion would
    read as an unremarkable property rather than as a defect removed.
    """
    production = Quadrature.product(n_mu=4, n_phi=8)
    level = np.asarray(production.level_indices[0])
    eta = production.eta[level]
    gaps = np.abs(np.diff(eta))
    tie_gaps = gaps[gaps < 1e-12]
    print(f"production eta gaps at the mirror pairs: {tie_gaps}")
    assert tie_gaps.size > 0, "no near-ties at all — fixture stopped working"
    assert np.all(tie_gaps == 0.0), (
        f"expected EXACT ties on roots-of-unity azimuths, got {tie_gaps}"
    )
    assert len(np.unique(eta)) == 8 // 2 + 1

    # The tie-break is now reachable, which is the whole consequence.
    got = {}
    for tie_break in ("lexsort", "stable"):
        with product_level_ordering(tie_break, exact_nodes=True):
            got[tie_break] = [
                list(map(int, a))
                for a in Quadrature.product(n_mu=4, n_phi=8).level_indices
            ]
    assert got["lexsort"] != got["stable"], (
        f"exact nodes left the tie-break UNREACHABLE: {got}"
    )


# ── 1. the tie-break cannot move alpha or tau ───────────────────────────

def test_alpha_and_tau_are_bit_identical_across_tie_breaks():
    r"""The geometry coefficients do not see the tie-break at all.

    ⛔ **Narrowed at Q5.6.4 (2026-08-11).** The τ half used to compare two
    τ ARRAYS on this full-circle fixture.  The cell partition is now the
    ω-midpoint, and a full-circle level has no angular cells, so there is
    no τ array to compare — the producer refuses.  The claim survives in
    the form that is still meaningful and is arguably sharper: **the
    REFUSAL is itself tie-break-invariant**.  If a tie-break could turn a
    refusal into an acceptance (or vice versa), the ordering would be
    reaching the closure, which is exactly what this file exists to deny.

    α is unaffected and still compared bit-for-bit — it is computed
    inline here from ``(eta, w)``, needs no cell partition, and is the
    coefficient the tie-break could most plausibly move.
    """
    from orpheus.sn.angular.closure import morel_montry_tau_per_level

    def coefficients(tie_break):
        with product_level_ordering(tie_break):
            quad = Quadrature.product(n_mu=4, n_phi=8)
            alphas = []
            for level in quad.level_indices:
                eta, w = quad.mu_x[level], quad.weights[level]
                alpha = np.zeros(len(level) + 1)
                for m in range(len(level)):
                    alpha[m + 1] = alpha[m] - w[m] * eta[m]
                alphas.append(alpha)
            try:
                morel_montry_tau_per_level(quad, CoordSystem.CYLINDRICAL)
                tau_verdict = "accepted"
            except ValueError as exc:
                tau_verdict = f"refused: {'monotone arc' in str(exc)}"
        return np.array(alphas), tau_verdict

    alpha_lex, tau_lex = coefficients("lexsort")
    alpha_stable, tau_stable = coefficients("stable")
    np.testing.assert_array_equal(alpha_lex, alpha_stable)   # bit-identical
    assert tau_lex == tau_stable == "refused: True", (
        f"the cell partition's verdict is tie-break DEPENDENT — lexsort "
        f"gave {tau_lex!r}, stable gave {tau_stable!r}. The per-level "
        f"ordering must not be able to decide whether a level has angular "
        f"cells."
    )


def test_the_full_circle_double_cover_is_REFUSED_by_the_cell_partition():
    r"""The full-circle level has no angular cells, and saying so is the fix.

    ⛔ **Retired-with-tombstone at Q5.6.4 (2026-08-11):** this row used to
    be ``test_the_degenerate_eta_pair_splits_tau_into_one_and_one_half``,
    which asserted that a mirror pair's shared ``eta`` collapses its
    angular cell (``eta_edge[m+1] = (eta[m]+eta[m+1])/2 = eta[m]``,
    giving ``tau_raw`` the ``[0,1,0,1,…]`` double-cover fingerprint) and
    that the ``[1/2, 1]`` absorber then laundered it into
    ``{1, 1/2}``.  **Both halves of that thesis are gone:** the absorber
    RETIRED, and the cell partition is no longer taken in ``eta`` at all
    (:func:`~orpheus.sn.angular.closure.angular_cell_edges_per_level`
    takes the midpoint in :math:`\omega`), so a shared ``eta`` no longer
    collapses anything — the two mirror partners have DIFFERENT
    :math:`\omega`.

    The successor claim is strictly stronger, and it is what this file is
    about: the double cover is now **refused at the partition**, with a
    message naming the cause, rather than being silently laundered into a
    finite wrong answer 400 lines downstream.  A full-circle level
    carries :math:`\omega` of both signs, so "the midpoint in
    :math:`\omega`" is undefined for it — that is not a tolerance
    question, and no epsilon appears below.

    ⚠ Note ``product(2, 8)`` is ALSO refused earlier, at cylindrical
    ``SNMesh`` admission (``assert_carrying_quadrature``).  This row
    exercises the partition producer DIRECTLY, so it still pins the
    inner guard for any caller that reaches it without a mesh.
    """
    from orpheus.sn.angular.closure import (
        angular_cell_edges_per_level,
        morel_montry_tau_per_level,
    )
    quad = Quadrature.product(n_mu=2, n_phi=8)
    omega = np.arctan2(quad.mu_y[quad.level_indices[0]],
                       quad.mu_x[quad.level_indices[0]])
    assert np.any(omega < 0.0) and np.any(omega > 0.0), (
        f"fixture is not a double cover any more: omega={omega} has one "
        f"sign, so it cannot exercise the refusal this row pins"
    )
    for producer in (angular_cell_edges_per_level, morel_montry_tau_per_level):
        with pytest.raises(ValueError, match="not a monotone arc in omega"):
            producer(quad, CoordSystem.CYLINDRICAL)


# ── 2. the ansatzes live inside the symmetric sector ────────────────────

def _fold_parent(n_mu: int = 4, n_phi: int = 8) -> Quadrature:
    """The σ_y-closed PARENT of the production default — the same
    spelling :meth:`Quadrature.folded_product` quotients (GL(n_mu) ×
    staggered trapezoid(n_phi)), built full so the mirror pairing
    exists.  Staggered ⟹ Σ = ∅ ⟹ the pairing is a PERFECT matching
    (no fixed points) — strictly stronger than the pre-6.3
    NODE_ALIGNED fixture, whose Σ nodes were self-paired."""
    measure, structure = spherical_product(
        gauss_legendre_on_mu(n_mu),
        periodic_trapezoid(n_phi, shift=STAGGERED),
    )
    return Quadrature(measure=measure, level_structure=structure)


@pytest.mark.parametrize(
    "builder",
    [build_cylindrical_mms_case, build_cylindrical_anisotropic_mms_case],
    ids=["isotropic", "anisotropic"],
)
def test_production_mms_source_is_identical_on_the_mirror_pair(builder):
    """``Q_n`` and ``psi_ref`` depend on ``eta`` and ``xi**2`` only.

    This is the Mode-7 declaration made executable: the ansatz does not merely
    make the ordering error small, it makes the two ordinates the ordering
    permutes CARRY THE SAME DATA.  Everything downstream follows from here.

    Since the 6.3 flip the builders default to the σ_y-FOLDED rule, on
    which no mirror partner exists to compare against (the quotient
    keeps one representative per orbit) — so the parity claim is
    evaluated on the fold's PARENT rule, through the same production
    ``external_source`` path, on the same frozen case
    (``dataclasses.replace`` swaps only the quadrature, re-running
    construction).  A second leg pins the RESTRICTION: the folded
    default evaluates bit-identically to the parent at the kept nodes,
    so migrating the builders lost nothing.
    """
    case = builder()                                  # folded default
    parent_case = replace(case, quadrature=_fold_parent())
    quad = parent_case.quadrature
    partner = _xi_mirror_pairing(quad)
    np.testing.assert_array_equal(
        np.sort(partner), np.arange(len(partner)),
    )  # a permutation…
    assert not np.any(partner == np.arange(len(partner))), (
        "staggered parent must have NO σ_y-fixed ordinate (Σ = ∅)"
    )  # …with no fixed points
    mesh = case.build_mesh(20)
    Q_parent = parent_case.external_source(mesh)      # (2N, ng, nx)
    np.testing.assert_allclose(Q_parent, Q_parent[partner], rtol=0, atol=1e-15)
    np.testing.assert_array_equal(quad.weights, quad.weights[partner])

    # Restriction leg: each folded node IS a parent node (bit-exact
    # representative selection, Q5.1/Q5.3), and the folded case's
    # source rows equal the parent case's rows there — the σ_y-even
    # closed form restricted, nothing recomputed differently.
    folded = case.quadrature
    parent_nodes = np.column_stack([quad.mu_x, quad.mu_y, quad.mu_z])
    folded_nodes = np.column_stack(
        [folded.mu_x, folded.mu_y, folded.mu_z]
    )
    fold_to_parent = np.array([
        np.flatnonzero((parent_nodes == node).all(axis=1))[0]
        for node in folded_nodes
    ])
    Q_folded = case.external_source(mesh)             # (N, ng, nx)
    np.testing.assert_array_equal(Q_folded, Q_parent[fold_to_parent])


# ── 3. the MMS ladders are invariant (the headline negative result) ─────

def _ladder(case, n_cells):
    errors = []
    for nc in n_cells:
        mesh = case.build_mesh(nc)
        result = solve_sn_fixed_source(
            case.materials, mesh, case.quadrature,
            case.external_source(mesh), max_inner=800, inner_tol=1e-13,
        )
        errors.append(volume_weighted_l2(
            np.asarray(result.scalar_flux.values)[0, :],
            case.phi_exact(mesh.centers), mesh.volumes,
        ))
    return np.asarray(errors)


@pytest.mark.slow
@pytest.mark.parametrize(
    "builder,n_cells",
    [
        (build_cylindrical_mms_case, [20, 40, 80]),
        (build_cylindrical_anisotropic_mms_case, [10, 20, 40]),
    ],
    ids=["isotropic", "anisotropic"],
)
def test_curvilinear_mms_ladder_is_blind_to_the_tie_break(builder, n_cells):
    """The production MMS gate cannot adjudicate issue #326's ordering."""
    got = {}
    for tie_break in ("lexsort", "stable"):
        with product_level_ordering(tie_break):
            got[tie_break] = _ladder(builder(), n_cells)
    rel = float(np.max(np.abs(got["lexsort"] - got["stable"])
                       / np.abs(got["lexsort"])))
    print(f"lexsort={got['lexsort']}\nstable ={got['stable']}\nrel={rel:.3e}")
    assert rel < 1e-9, (
        f"the curvilinear MMS unexpectedly SEES the tie-break (rel {rel:.3e}) "
        f"— the Mode-7 blindness declaration this module exists to pin is no "
        f"longer true; re-adjudicate #326 against it"
    )


# ── 4. the xi-ODD companion: sees it, but still cannot adjudicate ───────

@dataclass
class _CylXiOddMMSCase:
    r"""``psi_n(r) = (A(r) + B(r) xi_n)/W`` — the same ``A``, ``B`` shapes as
    the production anisotropic case, but ODD in the azimuthal cosine, so the
    mirror pair carries DIFFERENT data.

    Substituting into ``eta d_r psi - (xi/r) d_omega psi + SigT psi
    = Qext + SigS phi/W`` with ``d_omega xi = eta`` and
    ``phi = sum_n w_n psi_n = A`` (since ``sum_n w_n xi_n = 0``)::

        Qext_n = [eta A' + eta xi B' - eta xi B/r
                  + (SigT - SigS) A + SigT xi B] / W

    ``A(R) = B(R) = 0`` and ``A(0) = B(0) = 0``, so psi vanishes at both r = 0
    (reflective) and r = R (vacuum) exactly as the production case does.
    """

    sigma_t: float
    sigma_s: float
    radius: float
    materials: dict
    mat_id: int
    quadrature: Quadrature

    def A(self, r):
        return np.sin(np.pi * np.asarray(r) / self.radius)

    def Ap(self, r):
        R = self.radius
        return (np.pi / R) * np.cos(np.pi * np.asarray(r) / R)

    def B(self, r):
        rr = np.asarray(r) / self.radius
        return rr * (1.0 - rr) * np.cos(np.pi * rr)

    def Bp(self, r):
        R = self.radius
        rr = np.asarray(r) / R
        return ((1.0 - 2.0 * rr) * np.cos(np.pi * rr) / R
                - (np.pi / R) * rr * (1.0 - rr) * np.sin(np.pi * rr))

    def phi_exact(self, r):
        return self.A(r)

    def build_mesh(self, n_cells: int) -> Mesh1D:
        return Mesh1D(
            edges=np.linspace(0.0, self.radius, n_cells + 1),
            mat_ids=np.full(n_cells, self.mat_id, dtype=int),
            coord=CoordSystem.CYLINDRICAL,
            bc_left=BC("reflective"), bc_right=BC("vacuum"),
        )

    def external_source(self, mesh: Mesh1D) -> np.ndarray:
        r = mesh.centers
        A_, Ap_, B_, Bp_ = self.A(r), self.Ap(r), self.B(r), self.Bp(r)
        eta, xi = self.quadrature.eta, self.quadrature.xi
        W = float(self.quadrature.weights.sum())
        Q = (
            eta[:, None] * Ap_[None, :]
            + (eta * xi)[:, None] * Bp_[None, :]
            - (eta * xi)[:, None] * (B_ / r)[None, :]
            + (self.sigma_t - self.sigma_s) * A_[None, :]
            + self.sigma_t * xi[:, None] * B_[None, :]
        ) / W
        return Q[:, None, :]


def _build_xi_odd_case(n_mu: int = 4, n_phi: int = 8) -> _CylXiOddMMSCase:
    return _CylXiOddMMSCase(
        sigma_t=1.0, sigma_s=0.5, radius=5.0,
        materials={1: _make_1g_mixture(1.0, 0.5)}, mat_id=1,
        quadrature=Quadrature.product(n_mu=n_mu, n_phi=n_phi),
    )


def test_xi_odd_companion_source_DOES_differ_on_the_mirror_pair():
    """Control: the companion ansatz really leaves the symmetric sector."""
    case = _build_xi_odd_case()
    partner = _xi_mirror_pairing(case.quadrature)
    Q = case.external_source(case.build_mesh(20))
    rel = float(np.max(np.abs(Q - Q[partner])) / np.max(np.abs(Q)))
    print(f"xi-odd companion: mirror-pair source difference {rel:.3e}")
    assert rel > 1e-2, (
        f"the xi-odd companion source is (nearly) mirror-symmetric "
        f"({rel:.3e}) — it would be as blind as the production ansatz"
    )


def test_the_xi_odd_companion_fixture_is_INADMISSIBLE_since_the_63_flip():
    r"""⛔ RETIRED-WITH-TOMBSTONE — the adjudicator outlived its question.

    This row was
    ``test_xi_odd_companion_sees_the_tie_break_but_does_not_adjudicate_it``
    (``@pytest.mark.slow``), and it RAN A SOLVE on ``product(4, 8)`` — a
    FULL-CIRCLE rule.  Its record, kept because the measurement was real
    and is the reason #326's ordering question closed the way it did:

        Measured 2026-08-01 (``product(4, 8)``, nx [10,20,40,80]):
          lexsort [0.0727, 0.0790, 0.0806, 0.0810]  orders [-0.12, -0.03, -0.007]
          stable  [0.0877, 0.0827, 0.0815, 0.0812]  orders [+0.08, +0.02, +0.005]
        20.6 % apart at nx=10, 0.28 % apart at nx=80, spatial order ~0 for
        both: they converge to the SAME angular floor from opposite sides.
        Neither is "correct"; the floor is the defect.

    **Two independent reasons it is retired, and it should have ridden the
    6.3 flip commit** (task #22 scheduled exactly this class — "the 3
    ``_XFAIL_326`` rows + mechanism tests RIDE THE 6.3 FLIP COMMIT as
    retire-with-tombstone"); it was missed because it is ``slow``-marked
    and the `-m "not slow"` ledger never selected it:

    1. `[M]` its fixture is **inadmissible**: a cylindrical ``SNMesh``
       refuses a full-circle rule at construction
       (``assert_carrying_quadrature``, the 6.3 flip), so the solve raises
       before any coefficient is computed.  ⚠ This has been RED in the
       slow tier since the flip, uncounted by the #51 ledger — found
       2026-08-11 while re-running the τ slice without the marker filter.
    2. Its question is **closed**: T22 adjudicated the per-level ordering
       (the fold makes the η-tie unspellable on arcs), so there is nothing
       left for an adjudicator to decide.

    The successor claim is the refusal itself — which is the honest
    statement, and it is fast, so the ``slow`` marker goes too.
    """
    with pytest.raises(ValueError, match="CARRYING"):
        _ladder(_build_xi_odd_case(), [10])
