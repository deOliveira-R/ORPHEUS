r"""Every production rule's ADVERTISED degree of exactness, measured.

Issue **#327**, acceptance criterion 1. A
:class:`~orpheus.numerics.exactness.ExactnessClaim` is a promise about a
function space: *every polynomial up to this degree integrates exactly*. Nothing
in the tree measured whether the shipped rules keep it — the existing coverage
gates the claim TYPE's algebra
(:mod:`tests.numerics.test_exactness`) and the operators' floating-point
realization (:mod:`tests.numerics.test_symmetry_exactness`), neither of which
can see a rule whose promise is simply false.

`[M]` it was false. ``level_symmetric`` advertised :math:`N-1` and delivered
**3 at every order** — an over-claim of 12 at :math:`S_{16}`. Fixed by #327 on
2026-08-06 (per-:math:`O_h`-orbit moment-matched weights); **every production
rule now measures exactly its advertised degree**, and the xfail set below is
empty because this module's own strict marks refused to let the fix land
without the tag correction.

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

`[M]` **before #327** ``level_symmetric`` failed BOTH, in opposite directions:
it over-claimed at :math:`N \ge 6` (unsafe) and **under-claimed at**
:math:`N = 2` (advertised 1, measured 3). Failing in both directions is what
identified the tag as a **formula describing a different rule** rather than a
mis-calibrated measurement of this one — ``N-1`` is the moment-matched
construction's degree, and the rule then assigned one equal weight to every
ordinate. They coincided only at :math:`N = 4`. A bare
``measured == advertised`` would have reddened identically and hidden the
direction, which is why the split is worth its extra row.

The controls are in the SAME parametrization, deliberately
==========================================================

``gauss_legendre``, ``lebedev`` and ``product`` are swept by the same body that
sweeps ``level_symmetric``. ``vv`` #17: a probe that reports a failure must be
shown capable of reporting a pass on the same code path, or "0 caught" and "the
instrument is dead" are indistinguishable. `[M]` all three controls measure
EXACTLY their advertised degree, which is what licenses reading
``level_symmetric``'s 3 as a property of the rule rather than of this file.

The xfail mechanism, and what it actually did
=============================================

The failing rows shipped as ``xfail(strict=True)`` while #327 was open. When the
fix landed, all nine became ``XPASS(strict)`` **failures** — the suite would not
go green until the advertised degree was corrected in the same change as the
weights. That is the mechanism working end to end, not a hypothetical: the tag
and the rule were physically unable to drift apart across the fix.

⛔ The first draft used the IMPERATIVE ``pytest.xfail(...)``, which aborts
unconditionally and can never report ``XPASS``. Under it, the fix would have
landed with all nine still reporting ``XFAIL`` and the stale ``N-1`` intact —
the exact "gate that cannot redden" class this module exists to catch, inside
the module itself.
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
    ("level_symmetric(10)", lambda: Quadrature.level_symmetric(10), 16),
    ("level_symmetric(12)", lambda: Quadrature.level_symmetric(12), 18),
]


def _advertised(quad: "Quadrature") -> int:
    claim = quad.measure.exactness
    assert claim is not None, "a production rule must carry an exactness claim"
    return int(claim.degree)


#: ⭐ **EMPTY since #327 landed (2026-08-06) — and the emptiness was FORCED.**
#:
#: These held the four over-claims and the one under-claim as strict xfails.
#: When the moment-matched weights landed, all nine rows turned into
#: ``XPASS(strict)`` failures and the suite refused to go green until the tag
#: was corrected alongside the rule. That is the whole reason the marks were
#: strict rather than imperative ``pytest.xfail()``, which would have kept
#: reporting ``XFAIL`` and let a stale tag through.
#:
#: Kept as empty dicts rather than deleted: they are the seam a future
#: over-claim is recorded in, and an empty seam states "nothing is
#: over-claiming today" where an absent one states nothing at all.
_OVER_CLAIMS: "dict[str, int]" = {}
_UNDER_CLAIMS: "dict[str, int]" = {}


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

    ==================  ==========  ====================  ===============
    rule                advertised  measured (post-#327)  before #327
    ==================  ==========  ====================  ===============
    gauss_legendre(8)   15          15                    15
    lebedev(17)         17          17                    17
    product(6,12)       11          11                    11
    level_symmetric(8)  7           **7**                 **3**
    ==================  ==========  ====================  ===============

    Without the first three, ``level_symmetric``'s number is a statement about
    this file. With them, it is a statement about the rule — and the last
    column is why the controls had to be here BEFORE the fix, not added with
    it: they are what licensed reading the pre-fix 3 as real.
    """
    for rule_id, build, ceiling, expected in (
        ("gauss_legendre(8)", lambda: Quadrature.gauss_legendre(8), 20, 15),
        ("lebedev(17)", lambda: Quadrature.lebedev(17), 22, 17),
        ("product(6,12)", lambda: Quadrature.product(n_mu=6, n_phi=12), 16, 11),
        ("level_symmetric(8)", lambda: Quadrature.level_symmetric(8), 14, 7),
    ):
        got = _measured_degree(build(), ceiling=ceiling)
        assert got == expected, (
            f"the probe measured {rule_id} at degree {got}, expected "
            f"{expected}. If a CONTROL moved, the probe is broken and every "
            f"verdict in this module is void; if only level_symmetric moved, "
            f"the rule changed and #327 wants re-reading."
        )


#: `[M]` 2026-08-06 — distinct weights == number of :math:`O_h` orbits.
_ORBIT_COUNT = {2: 1, 4: 1, 6: 2, 8: 3, 10: 4, 12: 5}


@pytest.mark.parametrize("n,orbits", sorted(_ORBIT_COUNT.items()))
def test_the_weights_are_ONE_PER_ORBIT_and_positive(n: int, orbits: int) -> None:
    r"""⭐ RE-POSED at #327 — this row asserted the OPPOSITE until 2026-08-06.

    It read ``test_level_symmetric_has_exactly_ONE_distinct_weight`` and pinned
    the equal-weight construction as the mechanism behind the degree-3, with the
    note: *"the day this row reddens, the weights became per-level and the
    degree rows above should XPASS in the same commit."* That is exactly what
    happened, in this commit — so the row is re-posed to the new mechanism
    rather than deleted, and the old claim is recorded here because it is how
    the fix was recognised.

    The claim now: :math:`O_h`-invariance forces the weight to be constant on
    each orbit, so the number of DISTINCT weights is the number of orbits — no
    fewer (that would be the old equal-weight collapse) and no more (that would
    be a weight varying within an orbit, i.e. a rule that is not
    :math:`O_h`-invariant at all, whatever it advertises).

    ⛔ **Positivity is asserted, not assumed.** ``φ = Σ wₙ ψₙ`` must be
    non-negative for a non-negative angular flux, and the boundary response
    kernels assert it. `[M]` the moment-matched solve goes NEGATIVE from
    :math:`S_{14}` (−0.027) — which is why the family refuses above
    :math:`S_{12}` rather than trading positivity for degree.
    """
    weights = np.asarray(Quadrature.level_symmetric(n).weights)
    distinct = np.unique(np.round(weights, 15))

    assert distinct.size == orbits, (
        f"level_symmetric({n}) has {distinct.size} distinct weights, expected "
        f"{orbits} (one per O_h orbit). Fewer means the equal-weight collapse "
        f"is back and the degree has dropped to 3; more means a weight varies "
        f"WITHIN an orbit, so the rule is not O_h-invariant."
    )
    assert float(np.min(weights)) > 0.0, (
        f"level_symmetric({n}) has a non-positive weight {float(np.min(weights))} "
        f"— a positive angular flux can then integrate to a negative scalar "
        f"flux, and the boundary response kernels assert positivity."
    )
    assert float(np.sum(weights)) == pytest.approx(4.0 * np.pi), (
        f"level_symmetric({n}) weight sum is not 4π"
    )


def test_the_family_REFUSES_where_it_has_no_positive_solution() -> None:
    r"""⭐ Above :math:`S_{12}` the construction does not exist, and says so.

    `[M]` the moment-matched solve yields a negative weight from :math:`S_{14}`
    up (−0.027 at 14, −0.018 at 16, −0.142 at 20). Rather than ship a second
    silent construction under one name, or a rule whose scalar flux can go
    negative, the family is **defined exactly where it is valid**.

    ⭐ The frontier is COMPUTED, never hardcoded — the builder reads the sign of
    its own solution — so it tracks the node set instead of going stale beside
    it. This row therefore pins the BEHAVIOUR (refusal, with a message naming
    positivity and an alternative), not the number 12.
    """
    with pytest.raises(ValueError, match="no POSITIVE solution"):
        Quadrature.level_symmetric(14)
    with pytest.raises(ValueError, match="lebedev|product"):
        Quadrature.level_symmetric(16)
