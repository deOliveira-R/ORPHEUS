r"""Every production rule's ADVERTISED degree of exactness, measured.

Issue **#327**, acceptance criterion 1. A
:class:`~orpheus.numerics.exactness.ExactnessClaim` is a promise about a
function space: *every polynomial up to this degree integrates exactly*. Nothing
in the tree measured whether the shipped rules keep it — the existing coverage
gates the claim TYPE's algebra
(:mod:`tests.numerics.test_exactness`) and the operators' floating-point
realization (:mod:`tests.numerics.test_symmetry_exactness`), neither of which
can see a rule whose promise is simply false.

`[M]` it was false. ``level_symmetric`` advertises :math:`N-1` and delivers
**3 at every order** — an over-claim of 12 at :math:`S_{16}`.

The reference is the CLOSED FORM, not another rule
==================================================

On :math:`S^2`, in gamma functions:

.. math::

    \int_{S^2} x^a y^b z^c \,\mathrm{d}\Omega =
    \begin{cases}
      2\,\dfrac{\Gamma(\frac{a+1}{2})\Gamma(\frac{b+1}{2})\Gamma(\frac{c+1}{2})}
               {\Gamma(\frac{a+b+c+3}{2})} & a,b,c \text{ all even} \\[2mm]
      0 & \text{otherwise}
    \end{cases}

and on :math:`[-1,1]`, :math:`\int \mu^k \mathrm{d}\mu = \frac{2}{k+1}` for even
:math:`k`, else 0. Both are structurally independent of every quadrature in the
tree (``vv-principles`` L11) — which matters, because the defect being caught is
*a rule that agrees with the other rules about what it can do*.

⭐ Two claims, not one — and they fail in opposite directions
============================================================

A bare ``measured == advertised`` conflates two consequences that a consumer
experiences completely differently, so they are separate rows:

==========================  ==================================================
``measured >= advertised``  **SAFETY.** The promise is kept. Breaking this is
                            a live wrong answer for anyone who trusted the tag
                            to size a calculation.
``measured == advertised``  **ACCURACY.** The promise is tight. Breaking this
                            upward is safe but wasteful — a selector picks a
                            bigger rule than it needs, or refuses a job it
                            could do.
==========================  ==================================================

`[M]` ``level_symmetric`` fails BOTH, in both directions: it over-claims at
:math:`N \ge 6` (unsafe) and **under-claims at** :math:`N = 2` (advertised 1,
measured 3). That is the tell that the tag is not a mis-calibrated measurement
but a **formula describing a different rule** — ``N-1`` is the Carlson–Lathrop
moment-matched construction's degree, and this rule assigns one equal weight to
every ordinate (`[M]` exactly ONE distinct weight at every order). The two
coincide only at :math:`N = 4`.

The controls are in the SAME parametrization, deliberately
==========================================================

``gauss_legendre``, ``lebedev`` and ``product`` are swept by the same body that
sweeps ``level_symmetric``. ``vv`` #17: a probe that reports a failure must be
shown capable of reporting a pass on the same code path, or "0 caught" and "the
instrument is dead" are indistinguishable. `[M]` all three controls measure
EXACTLY their advertised degree, which is what licenses reading
``level_symmetric``'s 3 as a property of the rule rather than of this file.

Why the failing rows are ``xfail(strict=True)`` and not deleted or relaxed
=========================================================================

The defect is LIVE and #327 is open. A strict xfail pins it: the row cannot
silently start passing (that is an ``XPASS`` failure, which forces the tag to be
corrected in the same change that fixes the rule), and it cannot be forgotten,
because the xfail list IS the todo list. Relaxing the assertion to ``>= 3``
would encode the defect as the contract.
"""

from __future__ import annotations

from math import gamma

import numpy as np
import pytest

from orpheus.numerics.quadrature import Quadrature

pytestmark = [pytest.mark.foundation]

#: Absolute floor for "exact". Well above FP noise on a degree-31 monomial sum
#: and far below any real degree deficit — `[M]` the level_symmetric failures
#: are O(1) relative, and #327 measured them robust across five orders of
#: magnitude of this tolerance (1e-12 … 1e-2 all give degree 3).
_ATOL = 1e-10


def _sphere_monomial(a: int, b: int, c: int) -> float:
    r"""``∫_{S²} xᵃyᵇzᶜ dΩ`` in closed form."""
    if a % 2 or b % 2 or c % 2:
        return 0.0
    return 2.0 * (
        gamma((a + 1) / 2) * gamma((b + 1) / 2) * gamma((c + 1) / 2)
        / gamma((a + b + c + 3) / 2)
    )


def _interval_monomial(k: int) -> float:
    r"""``∫₋₁¹ μᵏ dμ`` in closed form."""
    return 2.0 / (k + 1) if k % 2 == 0 else 0.0


def _measured_degree(quad: "Quadrature", *, ceiling: int) -> int:
    """The largest ``d`` with EVERY monomial of total degree ``≤ d`` exact.

    Swept from 0 upward and stopped at the first inexact degree — the claim is
    "every polynomial up to d", so a single failure at degree d caps the answer
    at d-1 regardless of what happens higher up. (A rule can be exact at some
    isolated higher degree by symmetry; that is not a degree of exactness.)
    """
    nodes = np.asarray(quad.nodes)
    weights = np.asarray(quad.weights)
    one_dimensional = nodes.ndim == 1

    best = -1
    for degree in range(ceiling + 1):
        for a in range(degree + 1):
            for b in range(degree + 1 - a):
                c = degree - a - b
                if one_dimensional:
                    if b or c:
                        continue
                    got = float(np.sum(weights * nodes**a))
                    want = _interval_monomial(a)
                else:
                    got = float(np.sum(
                        weights * nodes[:, 0] ** a
                        * nodes[:, 1] ** b * nodes[:, 2] ** c
                    ))
                    want = _sphere_monomial(a, b, c)
                if abs(got - want) > _ATOL * max(abs(want), 1.0):
                    return best
        best = degree
    return best


def _xfail(reason: str):
    """A STRICT xfail mark.

    ⛔ Deliberately the MARK and not the imperative ``pytest.xfail(...)``. The
    imperative form aborts the test unconditionally, so it can never report
    ``XPASS`` — a rule that got fixed would keep reporting ``XFAIL`` and the
    stale tag would sail through. That is the "gate that cannot redden" class
    this whole module exists to catch, and the first draft of this file had it.
    """
    return pytest.mark.xfail(strict=True, reason=reason)


_OVER_CLAIM_REASON = (
    "#327 OPEN: {rid} over-claims by {delta} — it is an equal-weight O_h orbit "
    "rule (ONE distinct weight at every order), not the Carlson-Lathrop "
    "moment-matched construction its advertised N-1 describes."
)

#: ``(id, build, ceiling)`` — swept to ~2N so a rule cannot pass by the sweep
#: stopping below its own claim.
_RULES = [
    ("gauss_legendre(4)", lambda: Quadrature.gauss_legendre(4), 12),
    ("gauss_legendre(8)", lambda: Quadrature.gauss_legendre(8), 20),
    ("gauss_legendre(16)", lambda: Quadrature.gauss_legendre(16), 36),
    ("lebedev(17)", lambda: Quadrature.lebedev(17), 22),
    ("product(4,8)", lambda: Quadrature.product(n_mu=4, n_phi=8), 12),
    ("product(6,12)", lambda: Quadrature.product(n_mu=6, n_phi=12), 16),
    ("level_symmetric(2)", lambda: Quadrature.level_symmetric(2), 8),
    ("level_symmetric(4)", lambda: Quadrature.level_symmetric(4), 10),
    ("level_symmetric(6)", lambda: Quadrature.level_symmetric(6), 12),
    ("level_symmetric(8)", lambda: Quadrature.level_symmetric(8), 14),
    ("level_symmetric(12)", lambda: Quadrature.level_symmetric(12), 18),
    ("level_symmetric(16)", lambda: Quadrature.level_symmetric(16), 22),
]


def _advertised(quad: "Quadrature") -> int:
    claim = quad.measure.exactness
    assert claim is not None, "a production rule must carry an exactness claim"
    return int(claim.degree)


#: `[M]` 2026-08-06 — the rules whose PROMISE is broken today (#327). Strict, so
#: fixing the rule without correcting the tag is an ``XPASS`` failure.
_OVER_CLAIMS = {
    "level_symmetric(6)": 2,
    "level_symmetric(8)": 4,
    "level_symmetric(12)": 8,
    "level_symmetric(16)": 12,
}
#: …and the one that under-claims. Safe, but the same root cause: the tag is a
#: formula for a rule this is not.
_UNDER_CLAIMS = {"level_symmetric(2)": -2}


_SAFETY_PARAMS = [
    pytest.param(
        rid, build, ceiling, id=rid,
        marks=(
            _xfail(_OVER_CLAIM_REASON.format(rid=rid, delta=_OVER_CLAIMS[rid]))
            if rid in _OVER_CLAIMS else ()
        ),
    )
    for rid, build, ceiling in _RULES
]

_TIGHTNESS_PARAMS = [
    pytest.param(
        rid, build, ceiling, id=rid,
        marks=(
            _xfail(
                f"#327 OPEN: {rid} advertised - measured = "
                f"{ {**_OVER_CLAIMS, **_UNDER_CLAIMS}[rid] }"
            )
            if rid in _OVER_CLAIMS or rid in _UNDER_CLAIMS else ()
        ),
    )
    for rid, build, ceiling in _RULES
]


@pytest.mark.parametrize("rule_id,build,ceiling", _SAFETY_PARAMS)
def test_the_promise_is_KEPT(rule_id: str, build, ceiling: int) -> None:
    r"""⭐ SAFETY: ``measured ≥ advertised``. Breaking this is a wrong answer.

    A consumer that reads ``degree_of_exactness`` to size a calculation — the
    selector does exactly this — gets a silently under-resolved result when the
    promise is broken, and refining the order does not converge it away because
    the deficit is order-independent.
    """
    quad = build()
    advertised, measured = _advertised(quad), _measured_degree(quad, ceiling=ceiling)
    assert measured >= advertised, (
        f"{rule_id} PROMISES degree {advertised} and delivers {measured} — an "
        f"over-claim of {advertised - measured}. Every consumer that sized a "
        f"calculation from the tag is under-resolved, and refining the order "
        f"will not fix it if the deficit is order-independent."
    )


@pytest.mark.parametrize("rule_id,build,ceiling", _TIGHTNESS_PARAMS)
def test_the_promise_is_TIGHT(rule_id: str, build, ceiling: int) -> None:
    r"""ACCURACY: ``measured == advertised``. Breaking this upward is waste.

    Separate from the safety row because the consequence is different in kind:
    an under-claiming rule is correct, but a selector reading it picks a larger
    rule than needed — or declines a job it could have done.
    """
    quad = build()
    advertised, measured = _advertised(quad), _measured_degree(quad, ceiling=ceiling)
    assert measured == advertised, (
        f"{rule_id} advertises {advertised}, measures {measured}"
    )


def test_the_probe_itself_is_sound() -> None:
    r"""⭐ The positive control (``vv`` #17), as ONE row with the numbers in it.

    The per-rule rows above would all read "pass" for a probe that could not
    detect anything — so this states, in one place and with values, that the
    probe DOES resolve a degree and that the three control families land exactly
    on their claims. `[M]` 2026-08-06:

    ==================  ==========  ========
    rule                advertised  measured
    ==================  ==========  ========
    gauss_legendre(8)   15          15
    lebedev(17)         17          17
    product(6,12)       11          11
    level_symmetric(8)  7           **3**
    ==================  ==========  ========

    Without the first three, ``level_symmetric``'s 3 is a statement about this
    file. With them, it is a statement about the rule.
    """
    for rule_id, build, ceiling, expected in (
        ("gauss_legendre(8)", lambda: Quadrature.gauss_legendre(8), 20, 15),
        ("lebedev(17)", lambda: Quadrature.lebedev(17), 22, 17),
        ("product(6,12)", lambda: Quadrature.product(n_mu=6, n_phi=12), 16, 11),
        ("level_symmetric(8)", lambda: Quadrature.level_symmetric(8), 14, 3),
    ):
        got = _measured_degree(build(), ceiling=ceiling)
        assert got == expected, (
            f"the probe measured {rule_id} at degree {got}, expected "
            f"{expected}. If a CONTROL moved, the probe is broken and every "
            f"verdict in this module is void; if only level_symmetric moved, "
            f"the rule changed and #327 wants re-reading."
        )


def test_level_symmetric_has_exactly_ONE_distinct_weight() -> None:
    r"""The mechanism behind the degree, pinned separately from the symptom.

    ``rules_sphere.py`` assigns ``w_octant = 4π/(8·n_octant)`` to every ordinate.
    Classical Carlson–Lathrop solves a moment-matching system for **per-level**
    weights; the degree-3 this rule does reach is free — it comes from the
    :math:`O_h` orbit symmetry, which any equal-weight orbit set with
    :math:`\sum w = 4\pi` achieves.

    Pinning the mechanism rather than only the degree is what makes the eventual
    fix legible: the day this row reddens, the weights became per-level and the
    degree rows above should XPASS in the same commit.
    """
    for n in (4, 8, 16):
        weights = np.asarray(Quadrature.level_symmetric(n).weights)
        distinct = np.unique(np.round(weights, 15))
        assert distinct.size == 1, (
            f"level_symmetric({n}) now has {distinct.size} distinct weights — "
            f"if these are moment-matched per-level weights, #327's fix has "
            f"landed and the exactness rows above must be re-adjudicated."
        )
        assert float(np.sum(weights)) == pytest.approx(4.0 * np.pi), (
            f"level_symmetric({n}) weight sum is not 4π"
        )
