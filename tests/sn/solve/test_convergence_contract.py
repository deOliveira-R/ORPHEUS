"""The convergence CONTRACT: a best-effort exit is distinguishable and audible.

Issues #340 (a truncated exit was indistinguishable from a certified one at
the public entries) and #342 (``solve_sn`` hardcoded ``converged=True``).

The defect these gates exist to prevent is not a wrong number — it is a
**silently arbitrary** one.  A solve that runs out of iterations returns a
mid-descent iterate; if nothing distinguishes it from the converged answer,
a caller asserts physics against whatever the loop happened to be holding
when the budget expired.  That is how
``test_d3_pure_absorber_per_ordinate_psi_exact`` came to assert a
closed-form identity against the 999th of 1631 needed sweeps and pass, for
months, by luck (``scratch/d3_absorber_diagnosis.md``).

The contract has three parts, and each gets a gate below:

1. **The loop records why it stopped** — ``power_iteration`` returns a
   :class:`~orpheus.numerics.eigenvalue.PowerIterationOutcome` carrying
   ``converged``, instead of a triple that discarded it.  With the fact at
   the source, ``converged=True`` is not something a consumer can assert.
2. **The flag is honest at every entry** — forward and adjoint, eigenvalue
   and fixed-source, all read that fact (or the one shared residual
   predicate) rather than transcribing their own.
3. **A truncated exit is audible** — it emits
   :class:`~orpheus.numerics.convergence.ConvergenceWarning`, escalatable
   to a hard failure with
   :data:`~orpheus.numerics.convergence.ESCALATION_FLAG`.

⭐ Every gate here is written to be REDDENABLE, which for a boolean-flag
contract needs care: asserting ``converged is True`` on a solve that
converges is satisfied by the very hardcoded ``True`` this file exists to
forbid.  So each honesty gate is a PAIR — a converging configuration and a
deliberately-starved one — and it is the *starved* leg that has teeth.  A
regression to ``converged=True`` reds the starved leg of every pair.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture, make_mixture
from orpheus.geometry import BC
from orpheus.numerics.convergence import ESCALATION_FLAG, ConvergenceWarning
from orpheus.numerics.eigenvalue import PowerIterationOutcome
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solver import solve_sn, solve_sn_fixed_source
from orpheus.transport.mesh.axis import AxisMesh

_REFL = BC("reflective")


def _absorber_2g():
    """Pure absorber, 2 groups — no scattering iteration, so the reflective
    boundary coupling is the ONLY loop and its budget is the only knob."""
    return make_mixture(
        sig_t=np.array([0.8, 1.6]), sig_c=np.array([0.8, 1.6]),
        sig_f=np.array([0.0, 0.0]), nu=np.array([0.0, 0.0]),
        chi=np.zeros(2), sig_s=np.zeros((2, 2)),
    )


def _d3_reflective_axes():
    """The SI-hard configuration: all-reflective 3-D, [M] 1631 sweeps needed.

    Zero leakage plus zero scattering leaves only the boundary coupling, and
    only absorption damps it — the slow mode is the DD face sawtooth.  This
    is what makes the box a reliable *starvation* fixture: a modest
    ``max_inner`` cannot converge it, by measured construction rather than
    by a magic number.
    """
    return tuple(
        AxisMesh(edges=np.linspace(0.0, e, n + 1), bc_low=_REFL, bc_high=_REFL)
        for e, n in zip((1.0, 2.0, 3.0), (3, 4, 5))
    )


def _uniform_source(quad, cells, Q_g=(1.0, 0.5)):
    W = float(np.sum(quad.weights))
    return np.broadcast_to(
        (np.asarray(Q_g) / W)[None, :, *([None] * len(cells))],
        (quad.N, 2, *cells),
    ).copy()


def _eigenvalue(*, max_outer: int, keff_tol: float):
    """A 2-D all-reflective k-eigenvalue solve on a FISSILE mixture.

    Mixture "A" deliberately: "B" and "C" have ``SigF = 0``, so a k-solve on
    them divides by a zero production rate and yields ``nan`` — which the
    starved leg would then "pass" for entirely the wrong reason.  A control
    that is green because the physics is degenerate is not a control.
    """
    return solve_sn(
        {0: get_mixture("A", "2g")},
        tuple(
            AxisMesh(edges=np.linspace(0.0, e, n + 1),
                     bc_low=_REFL, bc_high=_REFL)
            for e, n in zip((2.0, 1.5), (4, 3))
        ),
        Quadrature.level_symmetric(sn_order=4),
        max_outer=max_outer, keff_tol=keff_tol, inner_tol=1e-10,
    )


def _fixed_source(max_inner: int, **kw):
    quad = Quadrature.level_symmetric(sn_order=4)
    return solve_sn_fixed_source(
        {0: _absorber_2g()}, _d3_reflective_axes(), quad,
        external_source=_uniform_source(quad, (3, 4, 5)),
        boundary_condition="reflective",
        inner_tol=1e-13, max_inner=max_inner, **kw,
    )


# ─── 1. The primitive records the fact ───────────────────────────────────


@pytest.mark.foundation
class TestPowerIterationCarriesItsOutcome:
    """``power_iteration`` returns WHY it stopped, not just what it found."""

    def test_outcome_is_not_tuple_unpackable(self) -> None:
        """The old ``keff, hist, flux = power_iteration(...)`` idiom is DEAD.

        Not pedantry: that idiom is exactly how the flag got discarded at
        five call sites.  If the outcome destructured like the triple it
        replaced, every one of them would still compile while ignoring
        ``converged``, and the retirement would be cosmetic.
        """
        outcome = PowerIterationOutcome(
            keff=1.0, keff_history=[1.0], flux_distribution=np.ones(3),
            converged=True,
        )
        with pytest.raises(TypeError):
            _a, _b, _c = outcome           # type: ignore[misc]

    def test_a_starved_budget_reports_not_converged(self) -> None:
        """max_outer=1 cannot satisfy a criterion needing >=3 outers.

        The teeth of the whole file: ``SNSolver.converged`` returns False for
        ``iteration <= 2`` by construction, so this outcome is structurally
        unconverged and any hardcoded ``True`` anywhere on the path reds
        here.
        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConvergenceWarning)
            sol = _eigenvalue(max_outer=1, keff_tol=1e-10)
        assert sol.history is not None
        assert sol.history.converged is False

    def test_a_generous_budget_reports_converged(self) -> None:
        """The positive control — the flag is not simply pinned False."""
        sol = _eigenvalue(max_outer=200, keff_tol=1e-8)
        assert sol.history is not None
        assert sol.history.converged is True


# ─── 2. The flag is honest at the fixed-source entry ─────────────────────


@pytest.mark.foundation
class TestFixedSourceFlagIsHonest:
    """The pair that would have caught the d3 defect the day it appeared."""

    def test_starved_budget_reports_not_converged(self) -> None:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConvergenceWarning)
            sol = _fixed_source(max_inner=50)
        assert sol.history is not None
        assert sol.history.converged is False

    def test_sufficient_budget_reports_converged(self) -> None:
        """[M] 1631 sweeps needed; 4000 is ~2.5x headroom."""
        sol = _fixed_source(max_inner=4000)
        assert sol.history is not None
        assert sol.history.converged is True


# ─── 3. A truncated exit is AUDIBLE ──────────────────────────────────────


@pytest.mark.foundation
class TestTruncationIsAudible:
    """The warning exists, fires exactly when the flag is False, and says
    enough to act on."""

    def test_starved_solve_warns(self) -> None:
        with pytest.warns(ConvergenceWarning) as record:
            _fixed_source(max_inner=50)
        assert len(record) >= 1

    def test_converged_solve_is_SILENT(self) -> None:
        """The other half of the pair — a warning that always fires is noise,
        and noise gets filtered, which is how the next truncation goes
        unnoticed."""
        with warnings.catch_warnings():
            warnings.simplefilter("error", ConvergenceWarning)
            _fixed_source(max_inner=4000)   # must not raise

    def test_the_message_carries_budget_tolerance_and_distance(self) -> None:
        """A bare "did not converge" cannot be acted on.

        The reader needs to tell "one more sweep" from "diverging", so the
        message must name the budget that ran out, the tolerance that was
        missed, and how far the last iterate actually was.
        """
        with pytest.warns(ConvergenceWarning) as record:
            _fixed_source(max_inner=50)
        msg = str(record[0].message)
        assert "max_inner=50" in msg
        assert "tol=" in msg
        assert "last residual" in msg
        assert "BEST-EFFORT" in msg

    def test_it_is_escalatable_to_an_error(self) -> None:
        """Prove the CATEGORY escalates rather than merely being emitted.

        ⚠ This proves the mechanism, NOT the published recipe — it installs
        the filter through the in-process API, so it is green for any
        spelling whatsoever.  The string is a separate claim and has its own
        gate below; read the two together.
        """
        with warnings.catch_warnings():
            warnings.simplefilter("error", ConvergenceWarning)
            with pytest.raises(ConvergenceWarning):
                _fixed_source(max_inner=50)

    def test_the_published_escalation_flag_actually_parses(self) -> None:
        """The CI contract is a STRING, and a string can be unspellable.

        ⛔ 2026-08-09 (#340).  The recipe shipped as
        ``-W error::ConvergenceWarning`` at four sites — including inside
        the emitted warning message — and **it does not parse**: ``-W``
        resolves an undotted category against ``builtins``, so pytest dies
        at startup with ``AttributeError: module 'builtins' has no
        attribute 'ConvergenceWarning'`` and collects ZERO tests.  The CI
        contract was imaginary for as long as it was published, and the
        sibling gate above stayed green throughout because it never touched
        the string.

        So this gate consumes the STRING, through pytest's own parser.
        """
        from _pytest.config import UsageError, parse_warning_filter

        assert ESCALATION_FLAG.startswith("-W ")
        spec = ESCALATION_FLAG.removeprefix("-W ")

        # It parses, and resolves to the category we actually emit.
        assert parse_warning_filter(spec, escape=False)[2] is ConvergenceWarning

        # And the teeth: the short spelling that shipped must still be
        # rejected, so a future "simplification" back to it reds here
        # instead of silently disarming CI.
        with pytest.raises(UsageError):
            parse_warning_filter("error::ConvergenceWarning", escape=False)
