r"""The α-dome recursion and its admission contract.

:func:`~orpheus.geometry.reduced_operator.alpha_dome` is the
Lathrop--Carlson angular redistribution recursion
:math:`\alpha_{m+1/2} = \alpha_{m-1/2} - w_m\mu_m`, and
:func:`~orpheus.geometry.reduced_operator._assert_alpha_dome_closes` is
its admission contract :math:`\alpha_{M+1/2} = 0` — a *consequence* of
the measure's antisymmetry, not an axiom of the one-sided recursion.

The contract is tested as a real ``raise``, not a bare ``assert``: the
canonical runner is ``python -O``, which strips ``assert`` at compile
time, so a contract spelled that way does not run in the suite that
matters (`coding-standards`).  Both legs are here — a closing measure is
ACCEPTED, a non-closing one is REFUSED — per ``vv-principles`` #11.

⚠ **This file is the residue of a split.**  Until 2026-08-28 it also
carried the factory-binding and packet-plumbing tests; those moved with
their system under test to ``tests/sn/mesh/test_reduced_operator.py``
(un-weld arc P4.4).  These 7 cases stay with
:mod:`orpheus.geometry.reduced_operator`, and travel again to
``tests/sn/angular/`` at P4.2 when :math:`\alpha` itself does.

These tests are tagged ``@pytest.mark.foundation``.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import CoordSystem
from orpheus.geometry.reduced_operator import (
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
    def test_the_derivations_alpha_dome_IS_the_production_one(self):
        r"""Cardinal Rule 2 — the recursion had THREE spellings until
        2026-08-12 (sphere factory, cylinder factory, and the derivations
        analysis module), which is why the contract could be enforced on
        one arm only.

        ⚠ This row is **not** a two-implementation agreement gate and must
        never be described as one: the derivations name now delegates, so
        no input can make the two disagree (``coding-standards``' "single-
        sourcing a duplicate demotes every gate that compared its copies"
        — the demotion is CORRECT here and the fix stays). What it pins is
        that the delegation is still in place, i.e. that a fourth copy has
        not been re-introduced under the old name.
        """
        from orpheus.derivations.discrete.sn.angular_differencing import (
            alpha_dome as derivations_alpha_dome,
        )

        quad = Quadrature.gauss_legendre(8)
        np.testing.assert_array_equal(
            derivations_alpha_dome(quad.mu_x, quad.weights),
            alpha_dome(quad.mu_x, quad.weights),
        )

    @pytest.mark.foundation
    def test_the_derivations_recursion_stays_UNGUARDED(self):
        r"""The contract lives at the production ADMISSION point, not
        inside the recursion — deliberately.

        ``derivations.discrete.sn.angular_differencing``'s P0/P4 predicate
        ladder exists precisely to characterise measures whose dome does
        NOT close. A guard welded into the shared recursion would make
        that analysis unspellable, so the split into
        :func:`alpha_dome` (pure) + :func:`_assert_alpha_dome_closes`
        (admission) is load-bearing, not stylistic.
        """
        from orpheus.derivations.discrete.sn.angular_differencing import (
            alpha_dome as derivations_alpha_dome,
        )

        # A deliberately non-antisymmetric measure: the recursion must
        # still RETURN it (open dome and all), not raise.
        mu = np.array([-0.8, -0.3, 0.3, 0.8])
        w = np.array([0.50, 0.25, 0.25, 0.25])
        alpha = derivations_alpha_dome(mu, w)
        assert alpha[-1] == pytest.approx(0.2), (
            "the analysis recursion must report the OPEN dome, not refuse it"
        )
        # …and the production admission point must refuse that same dome.
        with pytest.raises(ValueError, match="does not close"):
            _assert_alpha_dome_closes(alpha, coord=CoordSystem.SPHERICAL)
