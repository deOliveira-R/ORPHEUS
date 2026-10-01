r"""CHARACTERISATION: the sign of the half-angle flux :math:`\hat\psi`
under the Morel–Montry angular recurrence.

The M-M angular recurrence :eq:`pole-mm-recurrence`

.. math::
   \hat\psi_{m+1/2} \;=\;
   \bigl(\psi_m - (1-\tau_m)\,\hat\psi_{m-1/2}\bigr)\,/\,\tau_m

is a **first-order linear recursion with amplification factor**
:math:`-(1-\tau_m)/\tau_m`.  Since the shipped cylinder chart puts
:math:`\tau` in :math:`[\tfrac14, \tfrac34]`, individual steps can
amplify by up to 3, with alternating sign.  Nothing in the literature
read for this seam — Reed & Lathrop 1970, Bailey–Morel–Chang 2010,
Lathrop 2000, Hébert 2009 — states a positivity condition for the
ANGULAR recurrence (the positivity literature is about the SPATIAL
closure).  **These rows therefore CHARACTERISE, they do not certify.**
No row here carries ``verifies(...)``: there is no equation whose truth
they establish.

⛔ **What this module measured, and it corrects the framing the seam
review carried in** (`[M]` 2026-08-11, all numbers below on the fixture
built by :func:`_heterogeneous_2g_cylinder`):

* on the **production value path** — a converged flux with the marched
  :math:`\psi_{1/2}` STATE as the seed — :math:`\hat\psi` is **strictly
  POSITIVE**, ``min = +0.1337 / +0.1286 / +0.1287`` at
  :math:`n_\varphi = 6/8/16`, i.e. within 12 % of ``min psi`` itself.
  The recurrence is well behaved when the seed is *consistent* with the
  flux it marches, because the fixed point of the recursion at flat
  :math:`\psi` is :math:`\hat\psi = \psi` and the converged state sits
  near it.
* with a **zero seed** — the legitimate ψ-independent COEFFICIENT use
  (the transpose walk's ``denom``-only state, where these values are
  never read as fluxes) — the SAME converged flux gives
  ``min = -12.09 / -16.35 / -25.89``, and the excursion GROWS with
  :math:`n_\varphi`.
* with a random non-negative iterate and a zero seed, ``-48 / -77 /
  -133``.

⟹ the observed :math:`\hat\psi` sign is a property of the **seed's
consistency**, not of the scheme; a reported "``min psi_hat ≈ -77``"
is a statement about an inconsistent seed, not about the production
value path.  Both regimes are pinned below so a change to either is
audible, and the *mechanism* behind both — the worst partial
amplification :math:`A(M) = \max_m \prod_{k \le m}(1-\tau_k)/\tau_k` —
is pinned solve-free against a structurally independent reference.

Fixture discipline: the solve-based rows use a **heterogeneous,
two-group** cylinder (`vv-principles` #3/#4, `AGENT.md` §0.6) at
:math:`n_\varphi \in \{6, 8, 16\}` — :math:`M = 3, 4, 8`, **both
parities**, and never :math:`n_\varphi = 4`, whose :math:`M = 2` is the
partition-blind degenerate case (see
``tests/gates/sn/sweep/test_angular_cell_partition.py::test_the_M2_fold_is_BLIND_to_the_partition_choice``).

Mutation-proof — measured 2026-08-11
====================================

Same 13-mutation battery and 298-row scope as
``tests/gates/sn/sweep/test_angular_cell_partition.py`` (see its docstring for
the mutation definitions).  Rows reddened per gate:

======================================================  == == == == === == == == == == == === ===
gate                                                    MC M1 M2 M3 M3b M4 M5 M6 M7 M8 M9 M10 M11
======================================================  == == == == === == == == == == == === ===
``cylinder_recurrence_amplification…closed_form``        7  .  .  .   .  7  8  8  8  8  8   .   .
``sphere_recurrence_is_neutrally_stable…``               4  5  5  4   .  .  .  .  .  2  5   .   .
``the_production_value_path_keeps_psi_half_POSITIVE``    1  .  .  .   .  .  3  3  .  .  .   3   .
``a_zero_seed_drives_psi_half_NEGATIVE…``                1  .  .  .   .  .  3  3  3  3  3   3   3
======================================================  == == == == === == == == == == == === ===

M10 divides the recurrence by :math:`1-\tau` instead of :math:`\tau`;
M11 drops the :math:`(1-\tau_m)` thread weight.

⛔ **Two measured blindnesses, stated so no audit counts them as
coverage** (`vv-principles` #20):

* ``the_production_value_path_keeps_psi_half_POSITIVE`` is **green under
  M4, M7, M8 and M9** — four partition/chart mutations.  A *consistent*
  seed keeps the recurrence positive across a wide family of charts, so
  this row is a regression detector for the seam's END state, NOT a
  chart gate.  The chart gates are the amplification row here and the
  value rows in ``test_tau_producer_equivalence.py`` /
  ``test_angular_cell_partition.py``.
* ``sphere_recurrence_is_neutrally_stable…`` catches M8
  (:math:`\tau \to 1-\tau`) at only **2 of 5** orders, and never through
  its product leg — :math:`\tau \to 1-\tau` INVERTS every ratio, leaving
  :math:`\prod (1-\tau)/\tau = 1` exactly.  Only the amplification-band
  leg sees it, and only where the inverted product's peak leaves
  :math:`(1, 10)`.  The dedicated orientation catcher is
  ``test_angular_cell_partition.py::test_the_march_orientation_sign_is_pinned``
  (12 of 12 rows red under M8).
"""

from __future__ import annotations

from fractions import Fraction

from orpheus.sn.angular.closure import angular_cell_edges_per_level
from tests._harness.float_bounds import Tracked, U, gamma

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, CoordSystem, StructuredGeometry
from orpheus.mesh import CellsByCount, Mesher
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn import solve_sn_fixed_source
from orpheus.sn.problem import SNProblem
from orpheus.sn.angular.closure import morel_montry_tau_per_level

pytestmark = pytest.mark.foundation

_EPS = float(np.finfo(float).eps)

# M = n_phi/2 per level: 3 (odd), 4 (even), 8 (even).  n_phi = 4 is
# excluded on purpose — see the module docstring.
_SOLVE_N_PHI = (6, 8, 16)


def _heterogeneous_2g_cylinder(n_phi: int, *, n_mu: int = 4, nx: int = 12):
    r"""A 2-region, 2-group, vacuum-outer cylinder on the folded rule.

    Fuel-like mixture A for :math:`r < 1`, moderator-like B outside, so
    the flux carries a real spatial gradient AND a real group shape —
    a homogeneous or 1-group fixture would null exactly the terms the
    angular recurrence redistributes (`vv-principles` #3/#4).
    """
    mesh = Mesher(StructuredGeometry.cylinder(
        (0.0, 1.0, 2.0), (0, 1), outer=BC("vacuum"),
    )).partition(CellsByCount.uniform_width(nx // 2)).mesh
    quad = Quadrature.folded_product(n_mu=n_mu, n_phi=n_phi)
    materials = {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")}
    return SNProblem(mesh, quad, materials)


def _converged_flux(sn):
    """``(psi, radial_characteristic)`` from a converged fixed-source solve.

    The ``radial_characteristic`` member is ``Optional`` on the Solution
    (it is ``None`` on a seedless mesh — the slab).  A folded cylinder
    CARRIES on every level since Q5.6.3, so ``None`` here would mean the
    marched-seed regime this module contrasts against does not exist and
    every row below would be measuring the zero-seed path twice.  Fail
    loudly rather than silently degrade.
    """
    solution = solve_sn_fixed_source(
        sn.materials, sn.mesh, sn.quad,
        external_source=np.ones((sn.quad.N, 2, sn.nx)),
        inner_tol=1e-12,
    )
    psi = np.asarray(solution.angular_flux.interior.values, dtype=float)
    radial_characteristic = solution.radial_characteristic
    if radial_characteristic is None:
        pytest.fail(
            "the folded-cylinder solve returned no radial_characteristic "
            "member, so there is no marched psi_1/2 seed to contrast the "
            "zero seed against — every row in this module would silently "
            "collapse onto the coefficient path"
        )
    return psi, radial_characteristic.interior


def _min_psi_half(closure, psi, radial_characteristic) -> float:
    """Drive the PRODUCTION recurrence and return the global ``min``."""
    state = closure.precompute_psi_state(
        psi, radial_characteristic=radial_characteristic
    )
    return min(float(grid.faces.min()) for grid in state)


def _amplification(tau_level: np.ndarray) -> float:
    r""":math:`A = \max_m \prod_{k\le m}\,(1-\tau_k)/\tau_k` — the worst
    partial amplification the recurrence applies to a seed error."""
    return float(np.max(np.cumprod((1.0 - tau_level) / tau_level)))


# ═══════════════════════════════════════════════════════════════════════
# The MECHANISM — solve-free, structurally independent, both arms
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("n_phi", (4, 6, 8, 10, 16, 18, 32, 64))
def test_cylinder_recurrence_amplification_matches_the_closed_form(n_phi):
    r"""⭐ The recurrence's worst partial amplification :math:`A(M)`,
    against the analytic τ closed form.

    **Pins** the quantity that explains every :math:`\hat\psi` excursion
    in this module — and it does so **without a solve**, so it is the
    cheap regression detector for the seam.  A seed error :math:`\delta`
    reaches half-face :math:`m` multiplied by
    :math:`\prod_{k\le m}-(1-\tau_k)/\tau_k`; its worst magnitude is
    :math:`A(M)`.  `[M]` 2026-08-11:

    ====  ===============  ====================
    M     A (production)   A (closed-form ref)
    ====  ===============  ====================
    2     2.414213562      2.414213562
    3     2.732050808      2.732050808
    4     3.359160854      3.359160854
    5     3.656875757      3.656875757
    8     4.728870031      4.728870031
    9     4.976120036      4.976120036
    16    6.679689249      6.679689249
    32    9.443672213      9.443672213
    ====  ===============  ====================

    The reference is built from
    :math:`\tau_m = \tfrac12 + \tfrac12\cot\omega_m\tan(\Delta\omega/4)`
    — the hand-derived arc closed form, sharing no code path with the
    producer (`vv-principles` L11).  Relative agreement `[M]`
    :math:`\le 1.2\times 10^{-13}` at M = 32; asserted at ``1e-11``.

    ⚠ The value ``9.44`` at :math:`M = 32` is the same number the Q5.6.4
    partition adjudication quotes as the ω-midpoint chart's
    "recurrence error-amplification" (`vv-principles` #24b) — this row
    is the committed home of that figure, which until now lived only in
    a skill and a plan.

    **Second leg — the END-TO-END identity.** :math:`\prod_{k}
    (1-\tau_k)/\tau_k = 1` **exactly**, because the reversal identity
    :math:`\tau_m + \tau_{M-1-m} = 1` makes
    :math:`(1-\tau_k) = \tau_{M-1-k}`, so the numerators are the
    denominators re-ordered.  The recurrence is therefore *neutrally
    stable across a whole level* while amplifying by :math:`A(M)` in the
    middle — which is why the marched (consistent) seed stays positive
    and a zero seed does not.  `[M]` residual :math:`\le 8.6e{-14}`;
    asserted at ``1e-12``.

    **Cannot catch**: (a) an error in the *seed*, the *α*-dome, or
    :math:`\Delta A/w` — this row sees only τ; (b) a τ error that
    preserves the multiset :math:`\{(1-\tau_k)/\tau_k\}` and its worst
    prefix — in particular :math:`\tau \to 1-\tau` **inverts every
    ratio**, so the second leg (total product) is invariant under it and
    only the first leg reds.  The dedicated orientation catcher is
    ``tests/gates/sn/sweep/test_angular_cell_partition.py::test_the_march_orientation_sign_is_pinned``.
    """
    quad = Quadrature.folded_product(n_mu=4, n_phi=n_phi)
    tau_levels = morel_montry_tau_per_level(quad, CoordSystem.CYLINDRICAL)

    for p, level_idx in enumerate(quad.level_indices):
        M = len(level_idx)
        omega = np.arctan2(quad.mu_y[level_idx], quad.mu_x[level_idx])
        tau_reference = (
            0.5 + 0.5 / np.tan(omega) * np.tan(np.pi / M / 4.0)
        )

        np.testing.assert_allclose(
            _amplification(tau_levels[p]),
            _amplification(tau_reference),
            rtol=1e-11, atol=0.0,
            err_msg=(
                f"n_phi={n_phi} level {p} (M={M}): the recurrence's worst "
                f"partial amplification A = max_m prod (1-tau)/tau departs "
                f"from the analytic-arc reference. Every psi_hat excursion "
                f"characterised in this module scales with A."
            ),
        )
        ratios = (1.0 - tau_levels[p]) / tau_levels[p]
        np.testing.assert_allclose(
            float(np.prod(ratios)), 1.0, rtol=0.0, atol=1e-12,
            err_msg=(
                f"n_phi={n_phi} level {p} (M={M}): the end-to-end product "
                f"of (1-tau)/tau is not 1 — the recurrence is no longer "
                f"neutrally stable across a level, which means the "
                f"reversal identity tau_m + tau_{{M-1-m}} = 1 has broken"
            ),
        )


@pytest.mark.parametrize("N", (4, 8, 16, 32, 64))
def test_sphere_recurrence_is_neutrally_stable_across_the_level(N):
    r"""The sphere arm of the end-to-end identity, and its :math:`A(N)` — every float claim carries a DERIVED bound.

    ``[M]`` 2026-08-11 :math:`A` = 1.639121 (N=4) / 2.015931 (8) / 2.510677
    (16) / 3.145681 (32) / 3.952582 (64). (Unchanged.)

    **The identity, exactly.** For a rule symmetric to the bit
    (:math:`\mu_i = -\mu_{N-1-i}`, :math:`w_i = w_{N-1-i}`) with exact
    weight sum :math:`S = \sum w_i` (the float weights read as rationals),
    the cumulative-weight edges :math:`e_n = -1 + \sum_{i<n} w_i` satisfy
    :math:`e_{N-1-m} = S - 2 - e_{m+1}`, so the exact Morel–Montry weights
    :math:`\tau_m = (\mu_m - e_m)/w_m` obey

    .. math::

       \tau_m + \tau_{N-1-m} = 1 - \frac{S-2}{w_m},
       \qquad
       \prod_k \frac{1-\tau_k}{\tau_k} = \prod_k\Bigl(1 + \frac{S-2}{w_k\,\tau_{N-1-k}}\Bigr),

    the reversal identity with its mass-defect term. Asserted with ``==`` in
    rational arithmetic: a theorem, so no tolerance exists to pick. At
    :math:`S = 2` it is the old statement, product exactly 1.

    **The float values, within derived radii.** Production marches the
    edges left to right, :math:`\hat e_n = \mathrm{fl}(\hat e_{n-1} +
    w_{n-1})`, so :math:`|\hat e_n - e_n| \le R_n = R_{n-1} +
    u|\hat e_n|/(1-u)` (running error analysis, Higham §3.3, from the
    computed edges), and :math:`\hat\tau_m =
    \mathrm{fl}(\mathrm{fl}(\mu_m - \hat e_m)/\mathrm{fl}(\hat e_{m+1} -
    \hat e_m))` carries the radius :class:`Tracked` propagates from
    :math:`R_m, R_{m+1}`. The ratios :math:`(1-\hat\tau)/\hat\tau` are
    tracked the same way, and ``np.prod`` adds :math:`\gamma_{N-1}`
    relative in any order (Lemma 3.1).

    Until 2026-10 this row asserted :math:`|\hat P - 1| \le 10^{-12}`: the
    float product's own rounding reaches 1.8e-12 at N = 64 on the shipped
    rule and 7.4e-12 on another eigensolver's, so the tolerance sat inside
    the instrument's noise band (red in 7 of 10 ±3-ULP perturbations). The
    derived radius at N = 64 is printed in the failure message.

    **Cannot catch** (unchanged in kind): a τ error that preserves the
    multiset of ratios — but :math:`\tau \to 1-\tau`, the inversion the
    product is blind to, is now caught by the per-ordinate value leg.
    """
    quad = Quadrature.gauss_legendre(N)
    (tau,) = morel_montry_tau_per_level(quad, CoordSystem.SPHERICAL)
    (edges,) = angular_cell_edges_per_level(quad, CoordSystem.SPHERICAL)
    mu = [Fraction(float(m)) for m in quad.mu_x]
    w = [Fraction(float(x)) for x in quad.weights]

    # The rule's symmetry, to the bit — the identity's hypothesis.
    if not (all(mu[i] == -mu[N - 1 - i] for i in range(N)) and all(w[i] == w[N - 1 - i] for i in range(N))):
        pytest.fail(f"GL-{N} is not exactly symmetric: the reversal identity has no hypothesis")

    # Exact τ and the exact identity.
    S = sum(w)
    e = [Fraction(-1)]
    for wi in w:
        e.append(e[-1] + wi)
    tau_x = [(mu[m] - e[m]) / w[m] for m in range(N)]
    for m in range(N):
        if tau_x[m] + tau_x[N - 1 - m] != 1 - (S - 2) / w[m]:
            pytest.fail(f"sphere N={N}: the exact reversal identity fails at m={m}")
    ratio_x = [(1 - t) / t for t in tau_x]
    P_x = Fraction(1)
    for r in ratio_x:
        P_x *= r
    P_identity = Fraction(1)
    for k in range(N):
        P_identity *= 1 + (S - 2) / (w[k] * tau_x[N - 1 - k])
    if P_x != P_identity:
        pytest.fail(f"sphere N={N}: prod (1-tau)/tau != the mass-defect identity, exactly")

    # Production τ̂ against the exact τ, within the running-error radius.
    R = [Fraction(0)]
    for n in range(1, N + 1):
        R.append(R[-1] + U * abs(Fraction(float(edges[n]))) / (1 - U))
    E = [Tracked(e[n], R[n]) for n in range(N + 1)]
    ratio_t = []
    for m in range(N):
        t = (Tracked.of(quad.mu_x[m]) - E[m]) / (E[m + 1] - E[m])
        if abs(Fraction(float(tau[m])) - t.v) > t.r:
            pytest.fail(
                f"sphere N={N}: tau[{m}] = {float(tau[m])!r} is not the Morel–Montry "
                f"weight {float(t.v)!r} within its derived radius {float(t.r):.3e}"
            )
        ratio_t.append((Tracked.const(1) - t) / t)

    # The float product, within the derived bound of the exact one.
    P_hat = Fraction(float(np.prod((1.0 - tau) / tau)))
    upper = Fraction(1)
    for rt in ratio_t:
        upper *= abs(rt.v) + rt.r
    bound = (upper - abs(P_x)) + gamma(N - 1) * upper
    if abs(P_hat - P_x) > bound:
        pytest.fail(
            f"sphere N={N}: prod (1-tau)/tau = {float(P_hat)!r} is "
            f"{float(abs(P_hat - P_x)):.3e} from the exact {float(P_x)!r}, outside the "
            f"derived bound {float(bound):.3e} — the M-M recurrence is no longer "
            f"neutrally stable across the level"
        )
    if bound > Fraction(1, 10**8):
        pytest.fail(f"sphere N={N}: derived bound {float(bound):.3e} is not a rounding-scale statement")

    amplification = _amplification(tau)
    if not 1.0 < amplification < 10.0:
        pytest.fail(
            f"sphere N={N}: the recurrence's worst partial amplification "
            f"A = {amplification!r} left the characterised band (1, 10)."
        )


# ═══════════════════════════════════════════════════════════════════════
# The two SEED regimes on a real heterogeneous 2G cylinder
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("n_phi", _SOLVE_N_PHI)
def test_the_production_value_path_keeps_psi_half_POSITIVE(n_phi: int):
    r"""⭐ With the MARCHED seed, :math:`\hat\psi > 0` on a converged
    heterogeneous 2G cylinder.

    This is the path every value-carrying matvec and sweep takes: the
    seed is the composite's first-class :math:`\psi_{1/2}` state
    (#282 route (a)), marched from the true :math:`q_{1/2}` source, so
    it is *consistent* with the flux the recurrence then walks.

    `[M]` 2026-08-11, ``min psi_hat`` and its ratio to ``min psi``:

    =======  ====  =============  ==================
    n_φ      M     min psi_hat    ÷ min psi
    =======  ====  =============  ==================
    6        3     ``+0.133705``  0.8841
    8        4     ``+0.128600``  0.9348
    16       8     ``+0.128651``  0.9837
    =======  ====  =============  ==================

    Asserted as ``min psi_hat > 0`` AND
    ``min psi_hat >= 0.5 * min psi`` — the ratio floor is the
    regression detector (≥ 1.77× margin at the tightest row), the strict
    positivity is the characterised fact.  Both are ONE-SIDED: an
    improvement (a larger ratio) keeps them green.

    ⚠ **CHARACTERISATION, not a contract.** No source read for this seam
    prescribes positivity for the angular recurrence, and it is NOT a
    theorem — the identical fixture with an inconsistent seed goes to
    ``-26`` (next row).  What this row pins is that the *production*
    combination of chart, seed and converged flux lands positive, so a
    change that makes it negative is audible instead of silent.

    **Cannot catch**: (a) a τ chart that is wrong but still yields a
    positive :math:`\hat\psi` — the chart's value gates live in
    ``test_tau_producer_equivalence.py``; (b) anything about the
    magnitude of :math:`\hat\psi` beyond the ``0.5 × min psi`` floor;
    (c) anything on the spherical arm.
    """
    sn = _heterogeneous_2g_cylinder(n_phi)
    psi, radial_characteristic = _converged_flux(sn)
    min_psi = float(psi.min())
    if min_psi <= 0.0:
        pytest.fail(
            f"n_phi={n_phi}: the FIXTURE is degenerate — the converged "
            f"angular flux itself has min {min_psi!r} <= 0, so a "
            f"positivity claim about psi_hat says nothing"
        )

    min_psi_half = _min_psi_half(
        sn.angular_closure, psi, radial_characteristic
    )
    if min_psi_half <= 0.0:
        pytest.fail(
            f"n_phi={n_phi}: the production value path (marched psi_1/2 "
            f"seed, converged heterogeneous 2G flux) now produces a "
            f"NEGATIVE half-angle flux: min psi_hat = {min_psi_half!r} "
            f"(min psi = {min_psi!r}).\n"
            f"`[M]` 2026-08-11 this was +0.1337/+0.1286/+0.1287 at "
            f"n_phi = 6/8/16. Something changed in the tau chart, the "
            f"cell partition, or the seed march. This is a "
            f"CHARACTERISATION, not a theorem — but the change should be "
            f"deliberate and recorded."
        )
    if min_psi_half < 0.5 * min_psi:
        pytest.fail(
            f"n_phi={n_phi}: min psi_hat = {min_psi_half!r} has fallen "
            f"below half of min psi = {min_psi!r} (ratio "
            f"{min_psi_half / min_psi:.4f}; `[M]` 2026-08-11 the ratio "
            f"was 0.8841/0.9348/0.9837 at n_phi = 6/8/16). The angular "
            f"recurrence is amplifying where it used not to."
        )


@pytest.mark.parametrize("n_phi", _SOLVE_N_PHI)
def test_a_zero_seed_drives_psi_half_NEGATIVE_by_the_amplification(n_phi):
    r"""The other regime: with the ψ-independent ZERO seed the SAME
    converged flux gives :math:`\hat\psi \ll 0`, bounded by
    :math:`A(M)\,\max\psi`.

    ``radial_characteristic=None`` is a legitimate production call — the
    ψ-independent COEFFICIENT state the transpose walk builds, where
    these half-faces are never read as fluxes.  It is used HERE as the
    controlled *inconsistent-seed* experiment: same chart, same flux,
    only the seed changed.

    `[M]` 2026-08-11:

    =======  ====  ==============  ========  ==========================
    n_φ      M     min psi_hat     A(M)      \|min\| ÷ (A·max psi)
    =======  ====  ==============  ========  ==========================
    6        3     ``-12.089129``  2.732051  0.620
    8        4     ``-16.351438``  3.359161  0.687
    16       8     ``-25.890124``  4.728870  0.770
    =======  ====  ==============  ========  ==========================

    Two assertions, and together they say *the sign is real and the
    magnitude is the amplification, not something worse*:

    * ``min psi_hat < 0`` — the characterised fact;
    * ``|min psi_hat| <= A(M) * max psi`` — the derived envelope
      (a seed error of at most ``max psi`` amplified by at most
      :math:`A(M)`), which `[M]` holds with ≥ 1.3× margin.

    ⚠ **If a future change makes this row RED by turning the recurrence
    positivity-preserving, that is a WELCOME event** — re-pose the row
    (and the module's framing) rather than relaxing it.  The row is
    written as an assert-the-defect gate on purpose: a non-strict xfail
    would flip to xpass silently and an imperative ``pytest.xfail``
    could never flip at all (`vv-principles` Mode 8, classes 4 and 9).

    **Cannot catch**: anything about the production VALUE path (the row
    above owns that); the envelope leg cannot catch a sign error inside
    :math:`A` because it takes an absolute value.
    """
    sn = _heterogeneous_2g_cylinder(n_phi)
    psi, _radial_characteristic = _converged_flux(sn)
    closure = sn.angular_closure

    min_psi_half = _min_psi_half(closure, psi, None)
    max_psi = float(np.abs(psi).max())
    amplification = max(
        _amplification(tau) for tau in closure._tau_per_level
    )

    if min_psi_half >= 0.0:
        pytest.fail(
            f"n_phi={n_phi}: the ZERO-seed (coefficient-path) half-angle "
            f"flux is no longer negative — min psi_hat = "
            f"{min_psi_half!r}, `[M]` 2026-08-11 it was "
            f"-12.09/-16.35/-25.89 at n_phi = 6/8/16.\n"
            f"This module's framing rests on the recurrence NOT being "
            f"positivity-preserving under an inconsistent seed. If the "
            f"scheme has gained that property, that is good news: "
            f"re-pose this row and the module docstring, do not relax "
            f"the assertion."
        )
    envelope = amplification * max_psi
    if abs(min_psi_half) > envelope:
        pytest.fail(
            f"n_phi={n_phi}: |min psi_hat| = {abs(min_psi_half)!r} has "
            f"escaped the derived amplification envelope A*max|psi| = "
            f"{amplification!r} * {max_psi!r} = {envelope!r}. The "
            f"excursion is no longer explained by the recurrence's own "
            f"worst partial amplification, so a second mechanism is at "
            f"work (a wrong seed shape, a wrong level partition, or a "
            f"tau outside its chart)."
        )
