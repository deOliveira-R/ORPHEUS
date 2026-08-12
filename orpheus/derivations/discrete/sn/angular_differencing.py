r"""Analysis and diagnostics for SN **angular differencing schemes**.

An angular differencing scheme is a triple:

* the **quadrature** — nodes and weights on the direction sphere;
* the **angular cell partition** — the cell each ordinate owns, produced
  in production by
  :func:`~orpheus.sn.sweep.pole_angular_closure.angular_cell_edges_per_level`;
* the **closure weight** :math:`\tau` — how the ordinate flux relates to
  the two half-angle fluxes bounding its cell.

This module is where that triple gets *interrogated*: the predicates a
scheme must satisfy, the diagnostics that quantify how badly it fails, and
— importantly — a written record of which diagnostics are **blind** on
which rules. It is analysis, not production: it CONSUMES the production
partition rather than re-deriving it (Cardinal Rule 2), so a "reference"
here can never silently drift into a second definition of the cell
boundary. The one deliberate exception is flagged inline.

.. warning::

   Replaces ``derivations/discrete/sn/contamination.py`` (retired
   2026-08-11, Q5.6.4), whose cylindrical arm had become **present-tense
   wrong**: it built the retired η-midpoint edges, so its τ disagreed with
   production by up to 6.8e-2. It had one accurate and load-bearing
   comment, preserved as :ref:`the incompatibility note <why-not-weights>`.

-------------------------------------------------------------------------
NOMENCLATURE — three different τ, two different β
-------------------------------------------------------------------------

Both letters are overloaded in this codebase AND in the literature, and in
both cases the collision has already cost real time. Spell out which one
you mean, every time.

``tau`` / :math:`\tau`
    1. **The angular closure weight** — THIS module's τ. A dimensionless
       barycentric coordinate in :math:`[0, 1]`, one per ordinate,
       relating the ordinate flux to its two half-angle fluxes. Produced
       by
       :func:`~orpheus.sn.sweep.pole_angular_closure.morel_montry_tau_per_level`.
    2. **Optical depth** :math:`\Sigma_t s` — the path integral of the
       total cross section. Unrelated; lives in ``peierls_nystrom``, MoC,
       and the ``transport`` spatial schemes.
    3. **A critical half-thickness in mean free paths** — the ``fn_method``
       slab/reflector unknown.

    ⚠ And ``tau_inv`` in ``sn/sweep/cache.py`` is :math:`1/\tau` of sense
    **(1)**, sitting a few files from code using sense (2). When touching
    either, read the surrounding units before assuming.

``beta`` / :math:`\beta`
    1. **The BMC contamination coefficient** (Bailey–Morel–Chang 2010
       Eq. 41 sphere / Eq. 75 cylinder) — ONE SCALAR per level, the
       :math:`J^{(2)}` contamination of the asymptotic diffusion limit.
       Zero iff :math:`\tau` is the Morel–Montry weight.
       ⟹ :func:`contamination_beta`.
    2. **The Lathrop α-defect** (Lathrop 2000 Eq. 25) — A SEQUENCE,
       :math:`\alpha_{m+1/2} = 1 - \eta^2_{m+1/2} + \beta_{m+1/2}`, i.e.
       α's pointwise departure from the exact tangential target. Zero iff
       :math:`\delta \equiv 0`, i.e. :math:`\tau \equiv \tfrac12`.
       ⟹ :func:`alpha_defect_beta`.

    ⛔ **These two are near-opposites and the confusion is not
    hypothetical.** (1) vanishes for the Morel–Montry τ; (2) vanishes for
    the *diamond* τ. Reading a claim about one as a claim about the other
    produced a full design cycle's worth of wrong conclusions on
    2026-08-11 before the two definitions were put side by side. The
    conditions also live at *different orders*: (1) is a leading-order
    diffusion-limit statement, (2) a first-order truncation statement.

    Lathrop's :math:`\delta` relates to τ affinely: :math:`\tau =
    (1 + \delta)/2`, so :math:`\delta = 0 \iff \tau = \tfrac12`.

-------------------------------------------------------------------------
THE PREDICATE LADDER — what a scheme must satisfy, rule-agnostically
-------------------------------------------------------------------------

Extracted from Lathrop 2000 (Eqs. 23, 25–26, 29, 53–54) and BMC 2010
(Eqs. 12, 41–43, 52, 74–75). Stated so they can be evaluated for ANY 1-D
angular march, because *the literature's particular edge recursion is one
solution to this system, not the system itself*:

======  ===============================================================  =================
 name   condition                                                        order
======  ===============================================================  =================
 P0     :math:`\sum w = |{\rm range}|`, :math:`\sum w\eta = 0`           α closes
 P1     :math:`c \equiv \sum_m w_m \eta_m^2 = \tfrac23`                  LEADING
 P2     :math:`\eta_m = \tau_m e_{m+1/2} + (1-\tau_m) e_{m-1/2}`         FIRST
 P3     :math:`\eta_m \in [e_{m-1/2}, e_{m+1/2}] \iff \tau \in [0,1]`    well-posedness
 P4     :math:`\alpha_{m+1/2} = \alpha_{m-1/2} - k\,w_m \eta_m`          conservation
======  ===============================================================  =================

**P1 is the diffusion condition** (the :math:`P_1` pair is recovered iff
:math:`3c/2 = 1`) and it is a property of **the quadrature alone** — no
edges, no τ. **P2 DETERMINES τ** once a partition is chosen: τ is the
barycentric coordinate, equivalently the unique closure weight exact for a
flux *affine in the radial cosine*. **P3 is the predicate the retired
``[½,1]`` absorber was groping for** — and note it is the honest one:
`[M]` at :math:`S_8` Gauss–Legendre four of eight Morel–Montry τ lie below
:math:`\tfrac12`, so ``[½,1]`` was never the admissible range in either
arm.

.. _why-not-weights:

**Why the cylinder cannot use BMC's cumulative-weight edges.** Inherited
verbatim from the retired module, because it is correct and load-bearing:
*"the weight-sum approach is wrong for cylindrical because weights are
uniform in* :math:`\varphi` *-space, not* :math:`\eta` *-space."* Sharpened
with the measurement: an arc cell's η-measure is
:math:`2\sin\theta \sin\omega_m \sin(\Delta\omega/2) \propto \sin\omega_m`,
**not** constant (spread 0.30→1.53 across a level), while a trapezoid
weight is. Accumulating weights in η therefore violates P3, and it gets
WORSE with refinement: `[M]` ordinates outside their own cell go
0/4 → 4/8 → 12/16 → 28/32 at :math:`n_\varphi = 8/16/32/64`, and the solve
diverges (NaN) from :math:`n_\varphi \ge 16`. ⟹ BMC Eq. (52) is not a law;
it is the statement that *in their* quadrature the weight equals the
cell's η-measure. Ours does not, so we satisfy the same **predicate** with
a different partition (the ω-midpoint).

-------------------------------------------------------------------------
⛔ WHICH DIAGNOSTICS ARE BLIND, AND WHERE
-------------------------------------------------------------------------

Recording this is the whole reason this module has a docstring this long.
A diagnostic that reads "fine" on a rule it cannot see is worse than no
diagnostic, because an audit trusts it.

* **:func:`contamination_beta` is IDENTICALLY ZERO on a σ_y-folded arc,
  for ANY antisymmetric edge set — including random garbage.** The fold
  makes the nodes antisymmetric and the α dome symmetric, so
  :math:`{\rm term}_{M-1-m} = -{\rm term}_m` and the sum cancels pairwise.
  `[M]` measured on ``folded_product(4, 16)`` level 0: chord
  ``+6.94e-18``, edges scaled 0.5× ``+3.47e-18``, edges CUBED
  ``+1.73e-18``, random antisymmetrised ``−3.47e-18``; only breaking
  antisymmetry moves it (``−3.53e-03``). ⟹ **never gate a cylinder
  partition on β.** Gate it on the unfolded parent, or gate the
  antisymmetry and say that is what you tested.
* **:func:`nu_closure_residual` is the diagnostic that DOES discriminate**
  on a folded arc, and it needs no solve: the march implied by τ must land
  on the level's far endpoint. `[M]` ``1.000000`` for any derived τ,
  ``1.016389`` for the retired clamp, ``1.164784`` for :math:`\tau \equiv
  \tfrac12` at :math:`n_\varphi = 8` — i.e. those two correspond to no
  partition of the level at all.
* **A fixture with M = 2 ordinates per level cannot distinguish the
  ω-midpoint partition from the retired chord one.** `[M]` they are
  **bit-identical** at M = 2 (the interior chord midpoint
  :math:`(\eta_0+\eta_1)/2 = 0` IS the arc edge at
  :math:`\omega = \pi/2`), and diverge only from M = 3
  (``3.17e-02``, ``4.46e-02``, ``6.82e-02`` at M = 3/4/6). Any gate whose
  cylinder fixture is ``folded_product(·, 4)`` is silent on the partition
  choice — a Mode-12 blindness, not a tolerance question.
* **The MMS L2 norm is the wrong instrument for a closure re-pose.** It
  measures truncation order, which is what :math:`\tau = \tfrac12`
  optimises; it says nothing about the diffusion limit, which is what τ
  exists to fix. Principled ≠ more accurate.
"""

from __future__ import annotations

import numpy as np

from orpheus.geometry import CoordSystem
from orpheus.geometry.reduced_operator import alpha_dome as _production_alpha_dome
from orpheus.sn.sweep.pole_angular_closure import (
    angular_cell_edges_per_level,
    morel_montry_tau_per_level,
)

__all__ = [
    "alpha_defect_beta",
    "alpha_dome",
    "contamination_beta",
    "diffusion_limit_c",
    "morel_montry_beta",
    "morel_montry_weights",
    "nu_closure_residual",
]


def _coord_of(geometry: str) -> CoordSystem:
    """Map this module's legacy ``"spherical"``/``"cylindrical"`` strings."""
    try:
        return {
            "spherical": CoordSystem.SPHERICAL,
            "cylindrical": CoordSystem.CYLINDRICAL,
        }[geometry]
    except KeyError:
        raise ValueError(
            f"Unknown geometry: {geometry!r}; expected 'spherical' or "
            f"'cylindrical'."
        ) from None


def _levels(quad, coord: CoordSystem) -> tuple[np.ndarray, ...]:
    """The radial cosine per level, matching the partition's indexing."""
    if coord is CoordSystem.SPHERICAL:
        return (np.asarray(quad.mu_x),)
    return tuple(np.asarray(quad.mu_x)[idx] for idx in quad.level_indices)


def _weights(quad, coord: CoordSystem) -> tuple[np.ndarray, ...]:
    if coord is CoordSystem.SPHERICAL:
        return (np.asarray(quad.weights),)
    return tuple(np.asarray(quad.weights)[idx] for idx in quad.level_indices)


def alpha_dome(mu: np.ndarray, w: np.ndarray) -> np.ndarray:
    r"""P4's recursion: :math:`\alpha_{m+1/2} = \alpha_{m-1/2} - w_m\mu_m`.

    Returns ``(M+1,)`` with :math:`\alpha_0 = 0` and, when P0 holds,
    :math:`\alpha_M \approx 0` — the dome closes.

    ⭐ This is the analysis-facing name for the PRODUCTION recursion
    (:func:`orpheus.geometry.reduced_operator.alpha_dome`), which both
    curvilinear streaming factories run.  It was an independent third
    spelling of the same arithmetic until 2026-08-12; nothing compared
    the copies, so collapsing them demoted no gate.

    ⚠ It deliberately does NOT carry the production ADMISSION contract
    (``_assert_alpha_dome_closes``): this module's P0/P4 predicate ladder
    exists precisely to characterise measures whose dome does **not**
    close, and a guard welded into the recursion would make that
    analysis unspellable.  Ask :func:`alpha_defect_beta` for the defect.
    """
    return _production_alpha_dome(mu, w)


def diffusion_limit_c(quad, geometry: str = "spherical") -> np.ndarray:
    r"""**P1** — :math:`c = \sum_m w_m \eta_m^2`, one value per level.

    Lathrop 2000 Eq. 29 / Eqs. 53–54: the :math:`P_1` pair is recovered
    iff :math:`3c/2 = 1`, i.e. :math:`c = \tfrac23`.

    ⚠ **That constant belongs to the SPHERE's normalisation** — where the
    march covers :math:`\mu \in [-1, 1]` and :math:`\sum w = 2`, so
    :math:`c = \int_{-1}^{1}\mu^2\,d\mu = \tfrac23` exactly. `[M]`
    Gauss–Legendre delivers ``0.666666666667`` at N = 8 and 16.
    In general the target is :math:`|{\rm range}| / 3`.

    ⛔ Do NOT compare a CYLINDER level against ``2/3``: a level's radial
    cosine spans :math:`[-\sin\theta, \sin\theta]` and its weights sum to
    that level's azimuthal share, not to the full range — `[M]` the raw
    per-level moment reads ``0.282433`` on ``folded_product(4, ·)``,
    which is not a failure of P1 but a different normalisation. This
    function returns the RAW moment; scaling it is the caller's job, and
    the caller must say which range it is scaling by.

    A property of the QUADRATURE ONLY — no partition, no τ.  That is why
    a τ re-pose cannot break the diffusion limit on a fixed node set, and
    why "τ = ½ has the wrong diffusion limit" is a statement about the
    NODES that ``δ = 0`` implies (cell midpoints), not about τ.
    """
    coord = _coord_of(geometry)
    return np.array([
        float(np.sum(w * mu ** 2))
        for mu, w in zip(_levels(quad, coord), _weights(quad, coord))
    ])


def morel_montry_weights(quad, geometry: str = "spherical"):
    r"""**P2** — the closure weight τ, per level.

    ⚠ **This DELEGATES to production** (Cardinal Rule 2). The retired
    ``contamination.morel_montry_weights`` re-derived the edges itself and
    thereby became a second, divergent definition of the angular cell —
    `[M]` disagreeing with production by up to 6.8e-2 once the partition
    moved to ω. A "reference" that can drift is not a reference.

    ⟹ **Consequence for gates, stated so nobody has to rediscover it:**
    this is NO LONGER an independent reference for τ. A gate needing
    independence must use one of

    * the **analytic closed form** on a folded arc (structural
      independence, no shared code path)::

          tau_m = 1/2 + 1/2 * cot(omega_m) * tan(d_omega / 4)

    * or, on the sphere, a hand-written cumulative-weight expression.

    It remains useful as the τ *input* to the diagnostics below, which is
    what the surviving consumers want.

    Returns ``(N,)`` for the sphere, a list of ``(M_p,)`` for the cylinder
    (list, not tuple, to match the retired module's contract).
    """
    coord = _coord_of(geometry)
    tau = morel_montry_tau_per_level(quad, coord)
    return tau[0] if coord is CoordSystem.SPHERICAL else list(tau)


def contamination_beta(quad, geometry: str = "spherical"):
    r"""**β₁ — the BMC contamination coefficient.** One scalar per level.

    BMC 2010 Eq. 41 (sphere) / Eq. 75 (cylinder)::

        beta = sum_m eta_m [ alpha_{m+1/2} e_{m+1/2}
                             - alpha_{m-1/2} e_{m-1/2} ]

    (the spherical form carries a ½.)  :math:`\beta = 0` is
    diffusion-limit consistency — no Morel–Montry flux dip.

    ⛔ **BLIND on a σ_y-folded arc** — identically zero for ANY
    antisymmetric edge set, random values included. See the module
    docstring. Use :func:`nu_closure_residual` there.

    Returns a float for the sphere, ``(n_levels,)`` for the cylinder —
    the retired module's contract, preserved.
    """
    coord = _coord_of(geometry)
    edges = angular_cell_edges_per_level(quad, coord)
    out = []
    for mu, w, e in zip(
        _levels(quad, coord), _weights(quad, coord), edges
    ):
        alpha = alpha_dome(mu, w)
        scale = 0.5 if coord is CoordSystem.SPHERICAL else 1.0
        out.append(float(sum(
            scale * mu[m] * (alpha[m + 1] * e[m + 1] - alpha[m] * e[m])
            for m in range(mu.size)
        )))
    return out[0] if coord is CoordSystem.SPHERICAL else np.array(out)


def morel_montry_beta(
    mu: np.ndarray, w: np.ndarray, tau: np.ndarray,
) -> float:
    r"""**β₃ — M&M Eq. (6a) on the cell edges IMPLIED BY τ.** One scalar.

    .. math::

        \tilde\mu_{m+1/2} = \frac{\mu_m - (1-\tau_m)\,\tilde\mu_{m-1/2}}
                                 {\tau_m},
        \qquad \tilde\mu_{1/2} = -1

        \beta = 3 \sum_m \mu_m\bigl(\alpha_{m+1/2}\tilde\mu_{m+1/2}
                                  - \alpha_{m-1/2}\tilde\mu_{m-1/2}\bigr)

    with the standard recursion :math:`\alpha_{m+1/2} = \alpha_{m-1/2}
    - \mu_m W_m`, :math:`\alpha_{1/2} = 0`.  It is the :math:`2\beta/r`
    corruption of the S\ :sub:`N` diffusion coefficient
    :math:`D = 1/(3(\sigma_t + 2\beta/r))` (M&M Eq. 7a) — the flux dip.

    ``mu`` ascending; weights are renormalised here to :math:`\sum W = 1`
    (M&M's convention, in which the diffusion condition reads
    :math:`\sum \mu^2 W = 1/3`).

    ⭐⭐ **THE DIFFERENCE FROM ITS SIBLING** :func:`contamination_beta`
    **is the EDGE SET, and it is the whole reason one is τ-blind and this
    one is not.** Feeding the STANDARD weight-partition edges
    :math:`\mu_{m+1/2} = \mu_{m-1/2} + 2W_m` into Eq. (6a) **is** M&M's
    own proof that β = 0 (their Eqs. 17–19) — so a β built that way is
    blind to τ *by construction*, not by accident, and that is exactly
    why it sits in the campaign's instrument graveyard. The τ-IMPLIED
    edges above are what make this one read τ.

    `[M]` 2026-08-12: exactly zero (``< 1e-13``) for the shipped τ at
    Gauss S2/4/8/16/32, and non-zero for τ ≡ ½ — falling ~4 orders per
    doubling of N. It **catches** the ``τ → 1−τ`` reflection that the
    membership guard, the fold box and the reversal gates are all blind
    to (a Mode-12 stabiliser hole).

    ⛔ **Sphere only.** `[M]` on a σ_y-folded cylinder arc it is
    identically zero for BOTH the shipped τ and τ ≡ ½ at every
    :math:`n_\varphi` — refuted there, for the same antisymmetry reason
    that kills :func:`contamination_beta`.

    Gated by ``tests/sn/sweep/curvilinear/test_angular_beta_identity.py``
    (positive + negative legs); evidence
    ``scratch/q68_flux_dip_discriminator.md``.
    """
    mu = np.asarray(mu, dtype=float)
    tau = np.asarray(tau, dtype=float)
    normalized = np.asarray(w, dtype=float)
    normalized = normalized / normalized.sum()
    alpha = alpha_dome(mu, normalized)
    implied = np.zeros(mu.size + 1)
    implied[0] = -1.0
    for m in range(mu.size):
        implied[m + 1] = (mu[m] - (1.0 - tau[m]) * implied[m]) / tau[m]
    return float(
        3.0 * np.sum(mu * (alpha[1:] * implied[1:]
                           - alpha[:-1] * implied[:-1]))
    )


def alpha_defect_beta(quad, geometry: str = "spherical"):
    r"""**β₂ — the Lathrop α-defect.** A SEQUENCE per level, not a scalar.

    Lathrop 2000 Eq. 25: :math:`\alpha_{m+1/2} = 1 - \eta^2_{m+1/2}
    + \beta_{m+1/2}`, so this returns

    .. math:: \beta_{m+1/2} = \alpha_{m+1/2} - (1 - \eta_{m+1/2}^2)

    — α's pointwise departure from its exact tangential target at the
    cell edges.  Zero iff :math:`\delta \equiv 0`, i.e.
    :math:`\tau \equiv \tfrac12` (Eq. 26).

    ⛔ **This is NOT :func:`contamination_beta`.** That one vanishes for
    the Morel–Montry τ; this one vanishes for the *diamond* τ. They are
    near-opposites — see the module docstring's nomenclature section.

    Returns one ``(M+1,)`` array per level.
    """
    coord = _coord_of(geometry)
    return tuple(
        alpha_dome(mu, w) - (1.0 - np.asarray(e) ** 2)
        for mu, w, e in zip(
            _levels(quad, coord),
            _weights(quad, coord),
            angular_cell_edges_per_level(quad, coord),
        )
    )


def nu_closure_residual(quad, geometry: str = "spherical") -> np.ndarray:
    r"""⭐ The discriminator that survives the fold. One value per level.

    BMC Eq. 43 read as a RECURSION: given τ, the cell edges it implies are

    .. math::
        \nu_{1/2} = \eta_{\rm start}, \qquad
        \nu_{m+1/2} = \frac{\eta_m - (1-\tau_m)\,\nu_{m-1/2}}{\tau_m}

    and a τ that came from an actual partition must march this back onto
    the level's FAR endpoint.  Returns
    :math:`\nu_{M+1/2} / \eta_{\rm end}`, so **1.0 is exact closure**.

    Solve-free, needs no manufactured solution, and unlike
    :func:`contamination_beta` it is NOT annihilated by the fold's
    antisymmetry.  `[M]` on ``folded_product(4, n_phi)``:

    ============================  ========  ========  ========
    τ                             n_φ = 8   16        64
    ============================  ========  ========  ========
    derived (ω or chord)          1.000000  1.000000  1.000000
    the retired ``[½,1]`` clamp   1.016389  1.001930  1.000030
    :math:`\tau \equiv \tfrac12`  1.164784  1.039182  1.002412
    ============================  ========  ========  ========

    ⟹ the clamp and the diamond weight correspond to **no partition of
    the level**; that is the principled condemnation of the absorber, and
    it references no MMS.
    """
    coord = _coord_of(geometry)
    residuals = []
    for mu, tau, e in zip(
        _levels(quad, coord),
        (morel_montry_tau_per_level(quad, coord)),
        angular_cell_edges_per_level(quad, coord),
    ):
        nu = float(e[0])
        for m in range(mu.size):
            nu = (mu[m] - (1.0 - tau[m]) * nu) / tau[m]
        residuals.append(nu / float(e[-1]))
    return np.array(residuals)
