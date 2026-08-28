r"""The α-dome recursion and its admission contract.

:func:`~orpheus.sn.angular.redistribution.alpha_dome` is the
Lathrop--Carlson angular redistribution recursion
:math:`\alpha_{m+1/2} = \alpha_{m-1/2} - w_m\mu_m`, and
:func:`~orpheus.sn.angular.redistribution._assert_alpha_dome_closes` is
its admission contract :math:`\alpha_{M+1/2} = 0` — a *consequence* of
the measure's antisymmetry, not an axiom of the one-sided recursion.

The contract is tested as a real ``raise``, not a bare ``assert``: the
canonical runner is ``python -O``, which strips ``assert`` at compile
time, so a contract spelled that way does not run in the suite that
matters (`coding-standards`).  Both legs are here — a closing measure is
ACCEPTED, a non-closing one is REFUSED — per ``vv-principles`` #11.

⚠ **This file has moved twice, following its system under test.**  It was
``tests/geometry/test_reduced_operator.py``.  P4.4 (2026-08-28) took the
factory-binding and packet-plumbing tests out to
``tests/sn/mesh/test_reduced_operator.py`` with the chart connection; P4.2, the
same day, brought these 7 α cases here with
:mod:`orpheus.sn.angular.redistribution` itself.

These tests are tagged ``@pytest.mark.foundation``.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import CoordSystem
from orpheus.sn.angular.redistribution import (
    _ALPHA_CLOSURE_ATOL,
    _assert_alpha_dome_closes,
    alpha_dome,
)
from orpheus.numerics.quadrature import Quadrature



# ═══════════════════════════════════════════════════════════════════════
# The α-dome: ONE recursion, and a contract that survives ``-O``
# ═══════════════════════════════════════════════════════════════════════


class TestAlphaDomeClosureContract:
    r"""The admission contract :math:`\alpha_{M+1/2} = 0`, on both arms.

    ⛔ **What these rows exist to prevent, measured 2026-08-12.** The
    recursion :math:`\alpha_{m+1/2} = \alpha_{m-1/2} - w_m\mu_m` is
    strictly ONE-SIDED — :math:`\alpha_{1/2} = 0` is an axiom, the far-end
    closure is *not*; it is a consequence of the measure's antisymmetry
    and therefore a real contract on every admitted quadrature. Before
    this batch it was enforced:

    * on the SPHERE by a **bare** ``assert abs(alpha[N]) < 1e-12``;
    * on the CYLINDER **not at all**.

    And ``.claude/rules/vv-testing.md`` makes ``python -O -m pytest``
    canonical, which strips every ``assert``. `[M]` a measure closing at
    ``alpha[N] = +0.2`` was REFUSED under plain ``python`` and **ACCEPTED**
    under ``python -O`` — so the one check that existed did not run in the
    suite that matters.

    ⭐ **These rows are therefore only meaningful because they use
    ``pytest.raises``, not ``assert``.** A future contributor who
    "simplifies" the guard back to an ``assert`` reddens them under the
    canonical runner; that is the whole point of the class.
    """

    #: Every curvilinear rule the tree ships through these factories.
    _SPHERE_N = (4, 8, 16, 32, 64, 128)
    _CYL_RULES = ((4, 8), (4, 16), (4, 32), (8, 16), (6, 12), (2, 6))

    # ---- positive leg: the shipped rules are admitted --------------------
    # (`vv-principles` #11 — a negative-only guard test validates the
    # RAISING behaviour but never the claim that the contract is right.)

    @pytest.mark.foundation
    @pytest.mark.parametrize("N", _SPHERE_N)
    def test_every_shipped_gauss_legendre_dome_closes(self, N):
        """`[M]` residual 5.6e-17 (N=4) … 8.2e-17 (N=64) — ≈1 ULP of the
        dome peak, and NOT drifting with N."""
        quad = Quadrature.gauss_legendre(N)
        alpha = alpha_dome(quad.mu_x, quad.weights)
        _assert_alpha_dome_closes(alpha, coord=CoordSystem.SPHERICAL)
        assert abs(alpha[-1]) < 1e-15, (
            f"GL{N} dome residual {alpha[-1]!r} is far above the ~1e-16 "
            f"floor these rules have always delivered — the guard's 1e-12 "
            f"admission band would still pass it, so this row is the "
            f"tighter statement about the PRODUCERS (`vv-principles` #16: "
            f"the type permits 1e-12; Gauss-Legendre achieves 1e-16)."
        )

    @pytest.mark.foundation
    @pytest.mark.parametrize("n_mu,n_phi", _CYL_RULES)
    def test_every_shipped_folded_product_dome_closes_on_every_level(
        self, n_mu, n_phi,
    ):
        """`[M]` worst residual over the shipped rules is 2.78e-16 at
        ``folded_product(4, 32)``; ``(2, 6)`` closes at exactly 0.0."""
        quad = Quadrature.folded_product(n_mu=n_mu, n_phi=n_phi)
        for p, idx in enumerate(quad.level_indices):
            alpha = alpha_dome(quad.mu_x[idx], quad.weights[idx])
            _assert_alpha_dome_closes(
                alpha, coord=CoordSystem.CYLINDRICAL, level=p,
            )
            assert abs(alpha[-1]) < 1e-15, (
                f"folded_product({n_mu}, {n_phi}) level {p} residual "
                f"{alpha[-1]!r} is far above the ~1e-16 floor"
            )

    # ---- negative leg: an inadmissible dome is REFUSED, under -O ---------

    @pytest.mark.foundation
    @pytest.mark.parametrize("coord,level", [
        (CoordSystem.SPHERICAL, None),
        (CoordSystem.CYLINDRICAL, 3),
    ])
    def test_a_dome_that_does_not_close_is_refused(self, coord, level):
        """The row that the retired ``assert`` could not carry: this must
        raise under ``python -O`` as well as plain ``python``."""
        bad = np.zeros(5)
        bad[-1] = 0.2
        with pytest.raises(ValueError, match="does not close"):
            _assert_alpha_dome_closes(bad, coord=coord, level=level)

    @pytest.mark.foundation
    def test_the_cylinder_message_names_the_offending_level(self):
        """A folded rule closes on every level or none, but a
        level-symmetric rule can lose antisymmetry on ONE level — a
        whole-measure check would only say "somewhere"."""
        bad = np.zeros(5)
        bad[-1] = 0.2
        with pytest.raises(ValueError, match="on level 2"):
            _assert_alpha_dome_closes(
                bad, coord=CoordSystem.CYLINDRICAL, level=2,
            )

    @pytest.mark.foundation
    def test_the_guard_admits_exactly_up_to_its_own_stated_tolerance(self):
        """`vv-principles` #16 — the gate quotes the guard's OWN threshold
        rather than a tighter one the guard does not promise."""
        ok, bad = np.zeros(3), np.zeros(3)
        ok[-1] = 0.5 * _ALPHA_CLOSURE_ATOL
        bad[-1] = 2.0 * _ALPHA_CLOSURE_ATOL
        _assert_alpha_dome_closes(ok, coord=CoordSystem.SPHERICAL)  # no raise
        with pytest.raises(ValueError, match="does not close"):
            _assert_alpha_dome_closes(bad, coord=CoordSystem.SPHERICAL)

    # ---- the single-source claim ----------------------------------------

    @pytest.mark.foundation
    def test_the_derivations_module_does_not_SPELL_the_recursion(self):
        r"""Cardinal Rule 2 — the recursion had THREE spellings until
        2026-08-12 (sphere factory, cylinder factory, and the derivations
        analysis module), which is why the contract could be enforced on
        one arm only.

        ⭐ This gate was PROMOTED on 2026-08-28 (P4.2) and its claim class
        changed, so its description changed with it (``coding-standards``'
        mirror clause).  It used to assert that the derivations name
        *delegated* to production — a comparison no input could fail, since
        the body was ``return _production_alpha_dome(mu, w)``.  When
        :math:`\alpha` moved to ``orpheus.sn.angular.redistribution``, that
        L0 module could no longer import it at all (``derivations`` is L0,
        ``sn`` is L3), so the wrapper retired and the graders now ACCEPT
        :math:`\alpha` as a keyword.

        ⟹ what is pinned now is that the twin is **unspellable**, not merely
        delegating: the module must define no recursion of its own, and its
        three graders must be *handed* the dome.  A fourth copy
        re-introduced under any name would have to come back through one of
        these two doors.
        """
        import inspect

        from orpheus.derivations.discrete.sn import angular_differencing as ad

        assert not hasattr(ad, "alpha_dome"), (
            "the L0 analysis module has re-acquired an ``alpha_dome`` — that "
            "is the Pattern-2 twin the 2026-08-12 collapse retired, and it "
            "cannot import the production one (derivations is L0, sn is L3), "
            "so a local spelling is the only way it could have come back"
        )
        for fn in (ad.contamination_beta, ad.alpha_defect_beta,
                   ad.morel_montry_beta):
            param = inspect.signature(fn).parameters.get("alpha")
            assert param is not None, (
                f"{fn.__name__} no longer accepts ``alpha`` — it must be "
                "handed the dome it grades, not compute one"
            )
            assert param.kind is inspect.Parameter.KEYWORD_ONLY, (
                f"{fn.__name__}'s ``alpha`` must stay keyword-only, like "
                "``tau`` and ``edges``: a positional slot invites a caller "
                "to pass the wrong per-level sequence silently"
            )
            assert param.default is inspect.Parameter.empty, (
                f"{fn.__name__}'s ``alpha`` must stay REQUIRED — a default "
                "that fetched would reinstate the forbidden edge, and one "
                "that guessed would grade the wrong dome"
            )

    def test_the_derivations_recursion_stays_UNGUARDED(self):
        r"""The contract lives at the production ADMISSION point, not
        inside the recursion — deliberately.

        ``derivations.discrete.sn.angular_differencing``'s P0/P4 predicate
        ladder exists precisely to characterise measures whose dome does
        NOT close. A guard welded into the shared recursion would make
        that analysis unspellable, so the split into
        :func:`alpha_dome` (pure) + :func:`_assert_alpha_dome_closes`
        (admission) is load-bearing, not stylistic.

        ⚠ Re-pointed at the PRODUCTION recursion on 2026-08-28 (P4.2): the
        L0 spelling this used to call has retired, and the L0 ladder now
        runs on a dome its caller builds with *this* function.  So the
        claim is unchanged and its subject is now the only recursion there
        is — still two independently-evaluated functions, so welding the
        guard into :func:`alpha_dome` reddens the first assertion.
        """
        # A deliberately non-antisymmetric measure: the recursion must
        # still RETURN it (open dome and all), not raise.
        mu = np.array([-0.8, -0.3, 0.3, 0.8])
        w = np.array([0.50, 0.25, 0.25, 0.25])
        alpha = alpha_dome(mu, w)
        assert alpha[-1] == pytest.approx(0.2), (
            "the analysis recursion must report the OPEN dome, not refuse it"
        )
        # …and the production admission point must refuse that same dome.
        with pytest.raises(ValueError, match="does not close"):
            _assert_alpha_dome_closes(alpha, coord=CoordSystem.SPHERICAL)
