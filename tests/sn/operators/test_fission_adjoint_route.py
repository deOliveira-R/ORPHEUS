r"""**G5.7** — the one-level witness: ``IsotropicFission.apply_transpose`` is REACHED.

A claim about the **PATH**, not about a value (`vv` Mode-11 / anti-#26): a
correct implementation and one that computes the same numbers by routing
*around* the energy binding are indistinguishable in the return value, so the
instrument is a call counter, not an assertion on the answer.

**What the carve changed.** Until CS4c step 5 both fission-adjoint routes
reached the KERNEL and stepped over the operator that owns it:

* ``sn/solver.py`` posed the carrying-mesh ray emission as
  ``RadialCharacteristicEmission(F.kernel, …)`` — the
  :class:`~orpheus.numerics.operator.TensorProductOperator`, one level below
  the binding;
* ``FissionMomentOperator.apply_transpose`` read ``self.energy.kernel.apply_transpose``
  — the same step-over, in the moment factor.

`[M]` (the step-5 plan's census, over ``test_sn_adjoint_certification.py``):
``IsotropicScattering.apply_transpose`` **5309** calls,
``IsotropicN2N.apply_transpose`` **5309**, ``IsotropicFission.apply_transpose``
**0** — a verb that reads DEAD while its arithmetic is the hottest thing in the
adjoint. That is `vv` #29's fifth way (*a consumer fed the operator's KERNEL
instead of the operator*), and the honest repair is one-level delegation, not a
retirement. Both sites now pass the OPERATOR
(``fission.py`` ``FissionMomentOperator`` calls ``self.energy.apply_transpose``;
``sn/solver.py`` ``_adjoint_posing_parts`` builds
``RadialCharacteristicEmission(F.isotropic_energy, …)``).

**The gate, and why it is attributable.** A bare "count > 0" cannot say WHICH
route fired — and the two repairs are independent, so a single count would let
one of them silently regress. The spy therefore records the CALLER's frame, and
the rows assert a count per route:

============================ ============ ============
caller frame                 sphere       slab
============================ ============ ============
``FissionMomentOperator``    ``30``       ``9``
``RadialCharacteristicEmission`` ``30``   ``0``
============================ ============ ============

`[M]` 2026-09-04, ``get_mixture("A", "2g")``, 4 cells over ``[0, 5]``,
``gauss_legendre(4)``, ``keff_tol=1e-6 / flux_tol=1e-5 / inner_tol=1e-7``;
sphere ``k = 0.4845397047151718``, slab ``k = 1.8750000113975396``; **0.65 s**
and **0.41 s** respectively.

The SLAB is the negative control for the emission route (``solve_sn_adjoint``
takes the ``radial_characteristic_field_space is None`` early return, so no
``RadialCharacteristicEmission`` is ever built) — without it a green reading on
the sphere could not tell "the emission route fired" from "the moment factor
fired twice".

**The positive control** (`vv` #17): a counter says a call HAPPENED, not that
its result is consumed. Scaling the verb's output by 2 in-process moves the
adjoint eigenvalue by exactly 2× (`[M]` ``0.4845397047151718`` →
``0.9690791590419177``) — so the counted call is load-bearing, and a route that
counted but discarded would fail this row.

⚠ **Mutate in-process** (``monkeypatch``), never ``git checkout`` — the working
tree carries uncommitted edits (``process-discipline``).

`vv` Mode-8: ``pytest.fail`` / ``np.testing.*`` only (fire under ``python -O``).
"""

from __future__ import annotations

import sys
import warnings
from collections import Counter

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.sn.solver import solve_sn_adjoint
from orpheus.numerics.quadrature import Quadrature
from orpheus.transport.operators.isotropic_transfer import (
    IsotropicFission,
    IsotropicN2N,
    IsotropicScattering,
)

pytestmark = pytest.mark.foundation

_EMISSION = "RadialCharacteristicEmission.apply_transpose"
_MOMENT = "FissionMomentOperator.apply_transpose"

#: The run config's own convergence tolerance — the SINGLE SOURCE both the
#: solve and every tolerance below read (a k pinned tighter than the iteration
#: that produced it is pinning iteration noise, not a claim).
_KEFF_TOL = 1e-6
#: SAFETY × the run tolerance; the mutated solve carries its own increment
#: slack, so the doubled comparison gets 2 × the band.
_SAFETY = 10.0

#: `[M]` 2026-09-04 on the fixture below — the value leg's anchor for the
#: positive control (the shipped adjoint records are the campaign's real value
#: gate; this number exists only so the control can be shown to MOVE it).
_SPHERE_K = 0.4845397047151718


def _mesh(coord: CoordSystem) -> Mesh1D:
    """4 cells over ``[0, 5]`` — the cheapest carrying/non-carrying pair."""
    curvilinear = coord is not CoordSystem.CARTESIAN
    return Mesh1D(
        edges=np.linspace(0.0, 5.0, 5),
        mat_ids=np.zeros(4, dtype=int),
        coord=coord,
        **(
            {"bc_left": BC("reflective"), "bc_right": BC("vacuum")}
            if curvilinear
            else {}
        ),
    )


def _solve(coord: CoordSystem):
    return solve_sn_adjoint(
        {0: get_mixture("A", "2g")},
        _mesh(coord),
        Quadrature.gauss_legendre(n_ordinates=4),
        keff_tol=_KEFF_TOL,
        flux_tol=1e-5,
        inner_tol=1e-7,
        max_outer=50,
    )


def _install_counter(monkeypatch, cls, counter: Counter) -> None:
    r"""Wrap ``cls.apply_transpose`` through the descriptor protocol, keying the
    count on the CALLER's qualified frame name.

    Patched on the CLASS, so a ``cached_property``-held instance warmed by an
    earlier solve cannot mask it (`lessons L68b`): method lookup is dynamic, so
    the wrap bites for every instance, however it was minted.
    """
    original = cls.apply_transpose

    def spying(self, chi, /):
        counter[sys._getframe(1).f_code.co_qualname] += 1
        return original(self, chi)

    monkeypatch.setattr(cls, "apply_transpose", spying)


def test_the_carrying_adjoint_reaches_the_energy_binding(monkeypatch):
    r"""Both repaired routes enter ``IsotropicFission.apply_transpose``.

    The sphere adjoint runs the ray emission (``RadialCharacteristicEmission``,
    posed with ``F.isotropic_energy``) AND the moment factor
    (``FissionMomentOperator``, delegating to ``self.energy.apply_transpose``),
    so both frames must appear. **FIRST RED:** at ``f90f7914`` this counter
    reads 0 on both, on the very suite that runs 5309 calls of the sibling
    verb.
    """
    counter: Counter[str] = Counter()
    _install_counter(monkeypatch, IsotropicFission, counter)

    outcome = _solve(CoordSystem.SPHERICAL)
    if outcome.keff is None:
        pytest.fail("the sphere adjoint returned no eigenvalue")

    for frame in (_EMISSION, _MOMENT):
        if counter[frame] == 0:
            pytest.fail(
                f"IsotropicFission.apply_transpose was never entered from "
                f"{frame} — that route reads the KERNEL instead of the "
                f"OPERATOR again (`vv` #29's fifth way), and the energy "
                f"binding's transpose reads dead while its arithmetic runs. "
                f"Observed frames: {dict(counter)}"
            )


def test_the_sibling_verbs_are_untouched(monkeypatch):
    r"""ATTRIBUTION (`vv` #17): the loss-side siblings' counts are the row's
    denominator, not its subject.

    ``IsotropicScattering`` / ``IsotropicN2N`` already ran 5309 transposes each
    before the carve — so a gate that only asserted "some energy binding's
    transpose fired" would have been green all along. This row states the
    contrast explicitly: the siblings fire (they always did), and the row above
    is about the FISSION verb specifically.
    """
    counter: Counter[str] = Counter()
    _install_counter(monkeypatch, IsotropicScattering, counter)
    _install_counter(monkeypatch, IsotropicN2N, counter)

    _solve(CoordSystem.SPHERICAL)
    if sum(counter.values()) == 0:
        pytest.fail(
            "neither loss-side energy binding's transpose fired on the sphere "
            "adjoint — the fixture no longer exercises the adjoint gain, so "
            "the fission row above is not being compared against anything."
        )


def test_the_slab_is_the_emission_routes_negative_control(monkeypatch):
    r"""The NON-carrying mesh must NOT enter the emission route.

    ``solve_sn_adjoint`` branches on ``radial_characteristic_field_space is
    None`` (``sn/solver.py``), so a slab poses ``F`` directly and builds no
    ``RadialCharacteristicEmission`` at all. Without this row a green reading
    on the sphere cannot discriminate "the emission route fired" from "the
    moment factor fired twice" — the two frames would be indistinguishable in
    a bare total.

    `[M]` slab: ``FissionMomentOperator`` **9**, ``RadialCharacteristicEmission``
    **0**.
    """
    counter: Counter[str] = Counter()
    _install_counter(monkeypatch, IsotropicFission, counter)

    _solve(CoordSystem.CARTESIAN)

    if counter[_EMISSION] != 0:
        pytest.fail(
            f"the SLAB adjoint entered {_EMISSION} — the carrying branch is "
            f"being taken on a non-carrying mesh, so the sphere row's "
            f"attribution is void. Observed: {dict(counter)}"
        )
    if counter[_MOMENT] == 0:
        pytest.fail(
            f"the slab adjoint never entered {_MOMENT} — the moment factor's "
            f"one-level delegation regressed to reading the kernel. "
            f"Observed: {dict(counter)}"
        )


def test_the_counted_call_is_load_bearing(monkeypatch):
    r"""POSITIVE CONTROL (`vv` #17 / anti-#26): scaling the verb moves the answer.

    A counter proves a call HAPPENED; it cannot prove the result is CONSUMED.
    Doubling ``IsotropicFission.apply_transpose``'s output in-process must
    double the adjoint eigenvalue — the fission operator enters the adjoint
    pencil linearly, so ``2F†`` scales ``k`` exactly.

    `[M]` ``0.4845397047151718`` → ``0.9690791590419177``; the residual
    ``2.5e-07`` is the two solves' own increment slack (both stop at
    ``keff_tol``), which is why the band is derived from the RUN CONFIG
    (``SAFETY × keff_tol``, doubled for the doubled solve) and not invented —
    a k pinned tighter than the iteration that produced it pins noise, not a
    claim. The control's DISCRIMINATION is unaffected: a route whose result is
    discarded leaves ``k`` unchanged, i.e. off by ``0.48``, five orders above
    the band. Mutation is in-process (``monkeypatch``); the
    ``ConvergenceWarning`` a mutated pencil can raise is silenced HERE only,
    so the canonical ``-W`` wall stays intact elsewhere.
    """
    original = IsotropicFission.apply_transpose

    def doubled(self, chi, /):
        return 2.0 * np.asarray(original(self, chi))

    baseline = _solve(CoordSystem.SPHERICAL)
    if baseline.keff is None:
        pytest.fail("the sphere adjoint returned no eigenvalue")
    np.testing.assert_allclose(
        float(baseline.keff), _SPHERE_K, rtol=0.0, atol=_SAFETY * _KEFF_TOL,
        err_msg="the fixture's adjoint k moved — re-derive the control's "
                "expected ratio before trusting the mutation below.",
    )

    monkeypatch.setattr(IsotropicFission, "apply_transpose", doubled)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mutated = _solve(CoordSystem.SPHERICAL)
    if mutated.keff is None:
        pytest.fail("the mutated sphere adjoint returned no eigenvalue")
    np.testing.assert_allclose(
        float(mutated.keff), 2.0 * float(baseline.keff),
        rtol=0.0, atol=2.0 * _SAFETY * _KEFF_TOL,
        err_msg="doubling IsotropicFission.apply_transpose did NOT double the "
                "adjoint eigenvalue — the counted call's result is not the "
                "one the pencil consumes, so the route rows above count a "
                "call that does not matter.",
    )
