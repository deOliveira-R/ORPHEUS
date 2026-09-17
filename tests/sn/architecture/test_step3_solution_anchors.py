r"""Consumers-campaign step 3 — PRE-carve anchors for "the Solution carries
its posing".

Campaign: ``.claude/plans/consumers_step3_design.md`` (§3 the design, §5 the
rulings, §6.4 the landable units U1–U7).
Verification plan: ``scratch/_consumers/step3/test_architect_step3.md``.
Sibling step-2 anchors (the format precedent):
``tests/sn/architecture/test_step2_terminal_object_anchors.py`` and
``tests/sn/operators/test_step2_posed_fission_anchors.py``.

What step 3 does, in the domain's terms
=======================================

A Solution is the pair (Problem, **posing**) plus the Strategy that produced
it and the records.  Its KIND is the posing's TYPE — an eigen Solution carries
a ray representative, λ and the normalisation **gauge** that picked the
representative; a source Solution carries the coset representative and the
kernel gauge.  Never ``keff is not None``.

Today the kind is a PROPERTY of a nullable field
(:meth:`~orpheus.sn.solution.SolutionBase.is_eigenvalue` = ``keff is not
None``), the normalisation that picked the representative is recorded
NOWHERE, the multiplying-source entry measures its admissibility ``k`` and
drops it, and the two adjoint arms pose two different objects.  These rows
freeze that, exactly, so the carve is LOUD.

The two row kinds, and why each exists
======================================

**RECORD** (``TestRecord*``) — green at the time of writing, describing the
pre-carve tree exactly.  Each is *designed to RED at the carve* and is then
**DELETED, never repaired**: its job is to make the API change loud, because a
``strict`` xfail only flips on XPASS and is therefore SILENT when the carve
lands the API with the wrong semantics (``vv-principles`` Mode-8, fourth
class; the step-2 module's own pairing shape).

**RULED** (``TestRuled*``, ``xfail(strict=True)``) — the ruled post-carve
behaviour, RED until its unit lands; each is paired with the RECORD row that
states the pre-carve answer.  Every one is asserted over
``dataclasses.fields`` / TYPES / counting spies, **never over a name the
carve gets to choose** — step 2's O-2/O-4 lesson.

⚠ Mode-12 / config blindness, MEASURED, and it decides two fixtures
===================================================================

``[M]`` ``scratch/_consumers/step3/probes/p6_anchor_facts.py`` on the
2-group fuel|moderator slab: SN's ``compute_production_rate`` (fission **+**
the (n,2n) emission) and the fission-only ``IntegratedReactionRate`` read
**the same number** — ``0.9999999999999999`` both — because every shipped
library mixture ships ``Sig2 = 0``.  On the finalize module's
``_LIBRARY_N2N`` they differ by **1.218606e-01, rel 1.2186e-01**.  ⟹ *a
gauge row on a Σ₂-free mixture is structurally blind to half of the
functional it names*, so :class:`TestRecordTheFourProductionRateFunctionals`
carries BOTH populations and says which is the discriminating one.

``[M]`` same probe, the homogeneous family: ``IntegratedReactionRate.evaluate``
silently accepts the ``(ng,)`` flux ``HomogeneousResult`` stores and returns
``200.0`` where the ``(ng, 1)`` column reads ``100.0`` (2g; 411.27 vs 100.0 at
4g; **equal at 1g**).  ⟹ a 1-group gauge row cannot see it — the Cardinal
Rule's ≥2G demand, at the gauge tier.

Marks
=====

``foundation`` — architecture invariants of the Solution carrier and its
posing; no theory ``:label:``, hence no ``verifies(...)`` (the verifies ⊥
level doctrine).
"""

from __future__ import annotations

import dataclasses
import functools
import typing
import warnings

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry.mesh import BC, Mesh1D
from orpheus.geometry.coord import CoordSystem
from orpheus.numerics.coupled_system import CoupledField
from orpheus.numerics.eigenvalue import ProductionRateSolver, power_iteration
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solution import AdjointSolution, Solution, SolutionBase
from orpheus.sn.solver import (
    SNSolver,
    _adjoint_posing_parts,
    _as_sn_mesh,
    _build_fixed_source_rhs,
    _balance_projection,
    solve_sn,
    solve_sn_adjoint,
    solve_sn_fixed_source,
    solve_sn_multiplying_source,
)
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.full_field import FullField
from orpheus.transport.reaction_rate_functional import IntegratedReactionRate

pytestmark = pytest.mark.foundation


def _require(condition: object, message: str) -> None:
    """A ``-O``-firing assertion (NOT a bare ``assert``; ``vv`` Mode 8)."""
    if not condition:
        pytest.fail(message)


# ═══════════════════════════════════════════════════════════════════════
# Fixtures — fissile, ≥2G, heterogeneous (AGENT.md §0.6)
# ═══════════════════════════════════════════════════════════════════════
#
# ``lessons`` L26 / L7: a brief's "reuse the existing fixture" is a
# hypothesis.  Every ``tests/sn/architecture/_config`` mesh is built from
# ``_two_region_materials``, which carries NO fission — ``solve_sn`` on one
# raises.  Every row here that SOLVES builds its own fissile mesh from the
# shipped cross-section library, as the step-2 anchors do.

_QUAD_N = 8
#: ``[M]`` probes/p6: k = 0.907457573 at L = 4.0 — SUBCRITICAL, so the
#: multiplying-source entry is admitted, and 1/(1−k) = 10.8 makes the lagged
#: fission term the dominant one (the strong discriminator; the L = 2.0 slab's
#: k = 0.435 makes it a 1.8× perturbation and hides an ``A`` vs ``A − F``
#: difference).
_SUBCRITICAL_LENGTH = 4.0


def _slab(length: float = _SUBCRITICAL_LENGTH, bc_right: str = "vacuum"):
    """Fuel | moderator 2-G slab, GL-8, reflective | vacuum, 4 + 4 cells.

    The ``tests/sn/solve/test_subcritical_multiplying_source.py`` fixture,
    spelled once here (its own ``_slab``); the two modules pin the same k.
    """
    materials = {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")}
    mesh = Mesh1D(
        edges=np.linspace(0.0, length, 9),
        mat_ids=np.array([0] * 4 + [1] * 4, dtype=int),
        bc_left=BC("reflective"), bc_right=BC(bc_right),
    )
    return materials, mesh, Quadrature.gauss_legendre(_QUAD_N)


def _carrying_sphere():
    """A CARRYING (System-B bearing) hub — the arm whose adjoint iterates on
    the coupled carrier while the seedless one iterates on a bare field."""
    materials = {0: get_mixture("A", "2g")}
    mesh = Mesh1D(
        edges=np.linspace(0.0, 3.0, 7), mat_ids=np.zeros(6, dtype=int),
        coord=CoordSystem.SPHERICAL, bc_right=BC("vacuum"),
    )
    return materials, mesh, Quadrature.gauss_legendre(_QUAD_N)


def _uniform_source(quadrature: Quadrature, ng: int, nx: int) -> np.ndarray:
    return np.ones((quadrature.N, ng, nx))


# ═══════════════════════════════════════════════════════════════════════
# The two SCANNERS — shape-keyed, so the carve names its members freely
# ═══════════════════════════════════════════════════════════════════════
#
# Both ruled rows below assert over the answer's SHAPE rather than over a
# member NAME the carve gets to choose (step 2's O-2/O-4 lesson).  A filter
# is a claim, and a wrong one fails SILENTLY and FLATTERINGLY (it finds
# nothing, which reads as "the tree does not have this") — so each is
# hoisted here and validated against a PLANTED member by
# :class:`TestFilterTheScannersFindAPlantedMember`, which is a permanent
# THEOREM row (green before AND after the carve; never delete it).
#
# ⚠ Both walk ``dataclasses.fields`` / ``NamedTuple._fields``, i.e. DECLARED
# members — an attribute merely stapled onto an instance is invisible to
# them, and that is deliberate: the ruling is that the Solution *carries*
# these, not that something can be attached to it.

_SKIP_FIELDS = frozenset({"mesh", "angular_flux", "scalar_flux"})


@typing.runtime_checkable
class _SectionLike(typing.Protocol):
    """``ScaleGauge``'s ruled member set (§3.2), as a structural type — so the
    scanner's return type-checks without a suppression (``coding-elegance``
    anti-pattern #19)."""

    functional: typing.Callable[[typing.Any], float]
    target: float

    def apply(self, x: typing.Any) -> typing.Any: ...


def _declared_member_names(node: object) -> list[str]:
    if dataclasses.is_dataclass(node) and not isinstance(node, type):
        return [f.name for f in dataclasses.fields(node)]
    fields = getattr(node, "_fields", None)
    return list(fields) if isinstance(fields, tuple) else []


def _section_shaped_members(
    root: object, max_depth: int = 3,
) -> list[_SectionLike]:
    r"""Every declared member reachable from ``root`` that looks like a
    normalisation SECTION: it exposes ``functional``, a float ``target`` and a
    callable ``apply`` (``ScaleGauge``'s ruled member set, §3.2)."""
    found: list[_SectionLike] = []
    seen: set[int] = set()
    frontier: list[object] = [root]
    for _ in range(max_depth):
        nxt: list[object] = []
        for node in frontier:
            if node is None or id(node) in seen:
                continue
            seen.add(id(node))
            for name in _declared_member_names(node):
                member = getattr(node, name, None)
                if member is None:
                    continue
                if (
                    callable(getattr(member, "apply", None))
                    and getattr(member, "functional", None) is not None
                    and isinstance(getattr(member, "target", None), float)
                ):
                    found.append(member)
                nxt.append(member)
        frontier = nxt
    return found


def _declared_numeric_data(root: object, max_depth: int = 4) -> list[float]:
    """Every scalar reachable through DECLARED members, skipping the flux
    carriers and the hub (whose numbers are the Problem's, not the answer's)."""
    seen: set[int] = set()
    out: list[float] = []

    def walk(node: object, depth: int) -> None:
        if depth > max_depth or node is None or id(node) in seen:
            return
        seen.add(id(node))
        if isinstance(node, bool):
            return
        if isinstance(node, (float, int, np.floating, np.integer)):
            out.append(float(node))
            return
        if isinstance(node, (str, bytes, np.ndarray)):
            return
        if isinstance(node, (tuple, list)):
            for item in node[:64]:
                walk(item, depth + 1)
            return
        for name in _declared_member_names(node):
            if name in _SKIP_FIELDS:
                continue
            walk(getattr(node, name, None), depth + 1)

    walk(root, 0)
    return out


# ── the planted members the filter validation uses ────────────────────


@dataclasses.dataclass(frozen=True)
class _PlantedSection:
    functional: object
    target: float

    def apply(self, x: object) -> object:
        return x


@dataclasses.dataclass(frozen=True)
class _PlantedCertificate:
    bound: float
    by: str


@dataclasses.dataclass(frozen=True)
class _PlantedOutcome:
    state: object
    #: a name the carve need NOT choose — the scanner keys on shape
    picked_by: _PlantedSection


@dataclasses.dataclass(frozen=True)
class _PlantedSolution:
    mesh: object
    outcome: _PlantedOutcome
    certificate: _PlantedCertificate


class TestFilterTheScannersFindAPlantedMember:
    """THEOREM — green before AND after the carve.  Never delete.

    ``nexus-tools``/``vv`` #17: a filter is validated against a known member
    before any of its negatives is believed.  The two ruled rows below report
    *"nothing found"*, and *nothing found* is exactly what a broken filter
    reports — these rows are what separate the two readings.
    """

    @staticmethod
    def _planted() -> _PlantedSolution:
        section = _PlantedSection(functional=float, target=1.0)
        return _PlantedSolution(
            mesh=object(),
            outcome=_PlantedOutcome(state=np.ones(3), picked_by=section),
            certificate=_PlantedCertificate(bound=0.907457573, by="k-solve"),
        )

    def test_the_section_scanner_finds_a_planted_section(self) -> None:
        found = _section_shaped_members(self._planted())
        _require(
            len(found) == 1 and isinstance(found[0], _PlantedSection),
            f"the section scanner found {found!r} on a shape that plants "
            f"exactly one section — its negative on the real tree carries no "
            f"information until this passes.",
        )

    def test_the_section_scanner_finds_NOTHING_on_todays_carrier(self) -> None:
        """The negative control: a member-less object must read empty, or the
        scanner is matching on something incidental."""
        found = _section_shaped_members(_PlantedCertificate(1.0, "x"))
        _require(
            found == [],
            f"the section scanner matched {found!r} on a certificate-shaped "
            f"object with no section — it is over-matching.",
        )

    def test_the_numeric_scanner_finds_a_planted_scalar(self) -> None:
        values = _declared_numeric_data(self._planted())
        _require(
            any(abs(v - 0.907457573) < 1e-12 for v in values),
            f"the numeric scanner missed the planted k in {values!r}.",
        )

    def test_every_xfail_in_this_module_is_STRICT(self) -> None:
        """THEOREM — ``lessons`` L61g: a non-strict xfail reports ``x`` and is
        SILENT when its subject lands, so the marker set stops being the
        campaign's todo list.  ``pyproject.toml`` sets no ``xfail_strict``, so
        the default is non-strict and only this row catches a dropped flag."""
        import sys

        module = sys.modules[__name__]
        lax: list[str] = []
        for name in dir(module):
            obj = getattr(module, name)
            for mark in getattr(obj, "pytestmark", ()):
                if mark.name == "xfail" and mark.kwargs.get("strict") is not True:
                    lax.append(f"{name}: {mark.kwargs}")
        _require(
            lax == [],
            f"non-strict xfail marker(s) in this module: {lax}. A strict "
            f"marker XPASS-fails when its unit lands, which is the whole "
            f"reason the RULED rows are committed red.",
        )

    def test_the_numeric_scanner_skips_the_flux_carriers(self) -> None:
        """``_SKIP_FIELDS`` is part of the predicate — state it, and prove it
        bites, so a future widening cannot silently make the k-row vacuous."""
        planted = self._planted()
        _require(
            _SKIP_FIELDS == {"mesh", "angular_flux", "scalar_flux"},
            f"the scanner's skip set moved to {sorted(_SKIP_FIELDS)}.",
        )
        _require(
            len(_declared_numeric_data(planted)) < 8,
            "the numeric scanner returned a large set on a 3-scalar shape — "
            "it is walking into carriers it should skip.",
        )


@functools.lru_cache(maxsize=4)
def _hub_k() -> float:
    """``[M]`` the subcritical slab's k, through the hub's own eigen solve."""
    materials, mesh, quadrature = _slab()
    keff = solve_sn(materials, mesh, quadrature).keff
    if keff is None:  # pragma: no cover - the fixture is fissile by choice
        raise RuntimeError("the anchor fixture stopped being an eigenproblem")
    return float(keff)


@functools.lru_cache(maxsize=4)
def _multiplying_solution_truncated() -> Solution:
    """A deliberately TRUNCATED multiplying-source solve.

    Truncated because every ``None`` the certificate retires is reachable only
    off the converged path: ``_exit_balance_defect`` returns ``None`` when
    ``record.fully_converged`` (``solver.py:640``), so a converged fixture
    cannot see the balance number at all.
    """
    materials, mesh, quadrature = _slab()
    return solve_sn_multiplying_source(
        materials, mesh, quadrature,
        _uniform_source(quadrature, 2, 8),
        inner_tol=1e-12, max_inner=5,
    )


# ═══════════════════════════════════════════════════════════════════════
# RECORD — today's tree, stated exactly.  DELETE these at the carve.
# ═══════════════════════════════════════════════════════════════════════


class TestRecordTheCarriersShape:
    """The six fields, the two role leaves, and the three Optionals.

    ``[M]`` ``probes/p6_anchor_facts.py`` at ``d946ba9d``::

        angular_flux           'TimedFullField'
        scalar_flux            'ScalarFlux'
        mesh                   'SNMesh'
        keff                   float | None
        history                IterationHistory | None
        radial_characteristic  'RadialCharacteristicField | None'

    subclasses: ``['Solution', 'AdjointSolution']`` — ⚠ the campaign charter's
    *"every SolutionBase subclass across the six method families"* was a wrong
    premise; there are exactly TWO and both are SN (census §2.6).
    """

    _EXPECTED = (
        "angular_flux", "scalar_flux", "mesh", "keff", "history",
        "radial_characteristic",
    )

    def test_record_the_six_fields(self) -> None:
        """RECORD — the field NAMES and their order."""
        names = tuple(f.name for f in dataclasses.fields(SolutionBase))
        _require(
            names == self._EXPECTED,
            f"the carrier's field list moved: {names} != {self._EXPECTED}. "
            f"If step 3 landed, DELETE this row (its ruled successor is "
            f"TestRuledTheKindIsTheOutcomesType).",
        )

    def test_record_three_fields_are_optional_by_KIND(self) -> None:
        """RECORD — ``keff``/``history``/``radial_characteristic`` are
        ``| None``, and that Optionality is what encodes the kind today."""
        optional = tuple(
            f.name for f in dataclasses.fields(SolutionBase)
            if "None" in str(f.type)
        )
        _require(
            optional == ("keff", "history", "radial_characteristic"),
            f"the Optional-by-kind field set moved: {optional}.",
        )

    def test_record_exactly_two_role_leaves(self) -> None:
        """RECORD — the ROLE is a type (#276 A5); the KIND is not."""
        leaves = sorted(c.__name__ for c in SolutionBase.__subclasses__())
        _require(
            leaves == ["AdjointSolution", "Solution"],
            f"the role leaves moved: {leaves} — step 3 keeps TWO leaves "
            f"(the role axis), it does not add kind leaves (§3.2).",
        )
        _require(
            issubclass(Solution, SolutionBase)
            and issubclass(AdjointSolution, SolutionBase),
            "the two leaves must still be SolutionBase subclasses.",
        )

    def test_record_the_kind_is_read_off_the_nullable_keff(self) -> None:
        """RECORD — ``is_eigenvalue()`` IS ``keff is not None``, on a real
        pair of solves rather than by reading the source."""
        materials, mesh, quadrature = _slab()
        eigen = solve_sn(materials, mesh, quadrature)
        fixed = solve_sn_fixed_source(
            materials, mesh, quadrature, _uniform_source(quadrature, 2, 8),
        )
        _require(
            eigen.is_eigenvalue() and eigen.keff is not None,
            "the eigen entry stopped reporting a keff.",
        )
        _require(
            fixed.is_fixed_source() and fixed.keff is None,
            "the fixed-source entry stopped reporting keff=None.",
        )
        _require(
            type(eigen) is type(fixed) is Solution,
            f"today BOTH kinds are the same TYPE ({type(eigen).__name__} / "
            f"{type(fixed).__name__}) — that is exactly what step 3 ends.",
        )


class TestRecordTheMultiplyingEntryDropsItsK:
    r"""``solve_sn_multiplying_source`` measures its admissibility ``k`` and
    throws it away.

    ``[M]`` ``probes/p6_anchor_facts.py``: the returned ``Solution.keff`` is
    ``None`` and ``is_fixed_source()`` is ``True`` while the hub's own k-solve
    — which the entry RAN, at ``solver.py:3546``, to certify
    :math:`\rho(A^{-1}F) < 1` — reads ``0.9074575729668507``.  ⟹ the
    ``(M, q)`` cell's Solution is indistinguishable from a pure-transport one
    **by its own data**.
    """

    def test_record_the_admissibility_k_is_not_on_the_solution(self) -> None:
        """RECORD — the k is measurable only by RE-RUNNING the eigen solve."""
        solution = _multiplying_solution_truncated()
        _require(
            solution.keff is None,
            f"the multiplying entry now reports keff={solution.keff!r} — if "
            f"step 3 landed, DELETE this row.",
        )
        k = _hub_k()
        _require(
            0.0 < k < 1.0,
            f"non-vacuity: the fixture's k = {k!r} must be SUBCRITICAL, or "
            f"the entry refuses and this row measures a raise.",
        )
        _require(
            abs(k - 0.907457573) < 1e-6,
            f"the fixture's k moved: {k!r} vs the recorded 0.907457573.",
        )


class TestRecordTheMultiplyingEntryIsSilent:
    """The fifth entry bypasses the hoisted convergence warning.

    ``[M]`` ``probes/p6_anchor_facts.py`` on the SAME truncated configuration
    (``inner_tol=1e-12``, ``max_inner=5``):
    ``solve_sn_multiplying_source`` → ``converged=False``, **0 warnings**;
    ``solve_sn_fixed_source`` → ``converged=False``, **1 ConvergenceWarning**.
    The silence is gated as CORRECT today by the 7-site pin
    (``tests/numerics/test_family_convergence_contract.py:685``).
    """

    def test_record_the_multiplying_entry_warns_on_NOTHING(self) -> None:
        """RECORD — a truncated multiplying solve is silent; its sibling is
        not.  The paired reading is the point: a bare "0 warnings" could be a
        fixture that converged."""
        materials, mesh, quadrature = _slab()
        source = _uniform_source(quadrature, 2, 8)
        caught_by_entry: dict[str, list[str]] = {}
        converged: dict[str, bool] = {}
        for label, entry in (
            ("multiplying", solve_sn_multiplying_source),
            ("fixed_source", solve_sn_fixed_source),
        ):
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                solution = entry(
                    materials, mesh, quadrature, source,
                    inner_tol=1e-12, max_inner=5,
                )
            caught_by_entry[label] = sorted(
                type(w.message).__name__ for w in caught
            )
            converged[label] = solution.converged()
        _require(
            converged == {"multiplying": False, "fixed_source": False},
            f"non-vacuity: BOTH solves must be truncated, else 'no warning' "
            f"is the correct answer; got {converged}.",
        )
        _require(
            caught_by_entry["multiplying"] == [],
            f"the multiplying entry now warns ({caught_by_entry['multiplying']}) "
            f"— if step 3's U2 landed the hoisted warnings, DELETE this row.",
        )
        _require(
            caught_by_entry["fixed_source"] == ["ConvergenceWarning"],
            f"the SIBLING entry must warn on the same truncation, or this row "
            f"is measuring a fixture that cannot warn; got "
            f"{caught_by_entry['fixed_source']}.",
        )


class TestRecordTheMultiplyingEntryIsSilentOnAGAUGESINGULARHub:
    r"""The sharper half of the same silence — and the §6c witness U2's new
    ``warn_if_gauge_freedom`` row will need.

    A gauge-singular SUBCRITICAL FISSILE hub is CONSTRUCTIBLE (``[M]``
    ``probes/p7_gauge_singular_subcritical.py``: the ledger's 2-D
    all-reflective ``(3, 4)`` box with a dilute fissile mixture reads
    ``gauge_freedom(hub).present = True`` at every dilution and
    ``k = 0.003 … 0.15``), so the multiplying entry both RUNS and GAUGES there:

    =========================  =================  ==========
    entry                       gauge_correction   warnings
    =========================  =================  ==========
    ``…_multiplying_source``    6.080482952671906e-02   **none**
    ``…_fixed_source``          6.080482952672409e-02   ``GaugeFreedomWarning``
    =========================  =================  ==========

    ⟹ the repair FIRES on the fifth entry and nobody is told.  The paired
    reading is what makes the row non-vacuous: a bare "0 warnings" could be a
    hub with no kernel.
    """

    _CELLS = (3, 4)

    @staticmethod
    def _dilute_fissile(ng: int = 2):
        """Fissile enough to be admitted, dilute enough to stay SUBCRITICAL —
        an all-reflective box of any library mixture is supercritical
        (``k_inf = 1.875`` for ``A``) and the entry would REFUSE."""
        from orpheus.derivations.common.xs_library import make_mixture

        sig_t = np.linspace(0.8, 1.6, ng)
        sig_f = 0.05 * np.ones(ng)
        return make_mixture(
            sig_t=sig_t, sig_c=sig_t - sig_f, sig_f=sig_f,
            nu=2.4 * np.ones(ng),
            chi=np.array([1.0] + [0.0] * (ng - 1)),
            sig_s=np.zeros((ng, ng)),
        )

    def _hub_and_source(self):
        from orpheus.geometry import Mesh2D

        quadrature = Quadrature.level_symmetric(sn_order=4)
        reflective = BC("reflective")
        mesh = Mesh2D(
            edges_x=np.linspace(0.0, 1.0, self._CELLS[0] + 1),
            edges_y=np.linspace(0.0, 2.0, self._CELLS[1] + 1),
            mat_map=np.zeros(self._CELLS, dtype=int),
            bc_xmin=reflective, bc_xmax=reflective,
            bc_ymin=reflective, bc_ymax=reflective,
        )
        materials = {0: self._dilute_fissile()}
        source = np.full(
            (quadrature.weights.size, 2) + tuple(self._CELLS),
            1.0 / float(quadrature.weights.sum()),
        )
        return materials, mesh, quadrature, source

    def test_record_the_gauge_fires_and_the_entry_says_nothing(self) -> None:
        """RECORD — the gauge repair runs on the fifth entry; the warning the
        other four emit does not."""
        from orpheus.sn.operators.loss_kernel_gauge import gauge_freedom

        materials, mesh, quadrature, source = self._hub_and_source()
        hub = _as_sn_mesh(mesh, quadrature, materials, None)
        _require(
            gauge_freedom(hub).present,
            "non-vacuity: this hub is NOT gauge-singular, so 'no warning' is "
            "the correct answer and the row proves nothing (an ODD first axis "
            "and ≥2 reflective axis pairs are what excite the kernel).",
        )
        # ``lessons`` L44k: a ``**kwargs`` splat from an untyped dict is the
        # anti-#4 stringly-typed shape pyright cannot read — the two entries
        # are called EXPLICITLY, which also makes their differing signatures
        # (the fifth entry takes no ``inner_solver``) visible.
        def _run(is_multiplying: bool) -> tuple[list[str], float | None]:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                if is_multiplying:
                    solution: Solution = solve_sn_multiplying_source(
                        materials, mesh, quadrature, source,
                        boundary_condition=None,
                        inner_schedule="gauss_seidel",
                        inner_tol=1e-13, max_inner=400_000,
                    )
                else:
                    solution = solve_sn_fixed_source(
                        materials, mesh, quadrature, source,
                        boundary_condition=None,
                        inner_solver="source_iteration",
                        inner_schedule="gauss_seidel",
                        inner_tol=1e-13, max_inner=400_000,
                    )
            history = solution.history
            return (
                sorted(type(w.message).__name__ for w in caught),
                None if history is None else history.gauge_correction,
            )

        mult_warnings, mult_gauge = _run(True)
        fixed_warnings, fixed_gauge = _run(False)
        _require(
            mult_gauge is not None and mult_gauge > 1e-3,
            f"non-vacuity: the multiplying entry's gauge did NOT fire "
            f"({mult_gauge!r}) — the silence would then be correct.",
        )
        assert mult_gauge is not None  # narrowing for pyright
        _require(
            "GaugeFreedomWarning" in fixed_warnings,
            f"the SIBLING entry did not warn on this hub ({fixed_warnings}) — "
            f"the fixture cannot see the asymmetry this row is about.",
        )
        _require(
            mult_warnings == [],
            f"the multiplying entry now warns ({mult_warnings}) — if U2 "
            f"hoisted the warnings onto it, DELETE this row and keep the "
            f"fixture: it is the §6c witness for the new row.",
        )
        _require(
            fixed_gauge is not None
            and abs(mult_gauge - fixed_gauge) < 1e-9 * abs(fixed_gauge),
            f"the two entries gauge DIFFERENT amounts ({mult_gauge!r} vs "
            f"{fixed_gauge!r}) — they no longer share the exit path, so the "
            f"ledger's exemption prose for the fifth entry is stale.",
        )


class TestRecordTheFourProductionRateFunctionals:
    r"""The eigen gauge ships as several functionals under ONE word, and the
    Σ₂-free fixture cannot tell two of them apart.

    ``[M]`` ``probes/p6_anchor_facts.py``:

    ========================  =========================  =====================
    library                    SN ``compute_production``  fission-only ``R``
    ========================  =========================  =====================
    ``{A, B}`` 2g (Σ₂ = 0)     0.9999999999999999          0.9999999999999999
    finalize ``_LIBRARY_N2N``  1.0                         0.8781393737492378
    ========================  =========================  =====================

    ⟹ the SN forward gauge's (n,2n) half is INVISIBLE on every shipped
    library mixture; only the manufactured Σ₂ stack activates it (rel
    **1.2186e-01**).  Any step-3 gauge row written on the plain 2-G slab is a
    provable non-catcher for "which production rate?" — ``vv`` Mode 12 at the
    fixture.
    """

    def test_record_the_two_functionals_agree_on_a_SIGMA2_FREE_mixture(
        self,
    ) -> None:
        """RECORD — the blindness itself, stated as a measurement."""
        materials, mesh, quadrature = _slab()
        sn_mesh = _as_sn_mesh(mesh, quadrature, materials)
        solver = SNSolver(sn_mesh)
        phi = np.asarray(
            solve_sn(materials, mesh, quadrature).scalar_flux.values,
            dtype=float,
        )
        total = solver.compute_production_rate(phi)
        fission_only = float(
            IntegratedReactionRate(
                sn_mesh.mat_xs.fission_production_field
            ).evaluate(phi)
        )
        _require(
            total == fission_only,
            f"the Σ₂-free slab no longer makes the two functionals agree "
            f"({total!r} vs {fission_only!r}) — the library grew an (n,2n) "
            f"channel and this blindness record is stale.",
        )

    def test_record_a_SIGMA2_CARRYING_mixture_separates_them(self) -> None:
        """RECORD — the discriminating population, from the finalize module's
        own manufactured stack (ONE source, not a second copy)."""
        from tests.sn.solve.test_eigenvalue_finalize_reconstruction import (
            _LIBRARY_N2N,
            _slab as _finalize_slab,
        )

        mesh = _finalize_slab(BC.vacuum)
        quadrature = Quadrature.gauss_legendre(n_ordinates=_QUAD_N)
        sn_mesh = _as_sn_mesh(mesh, quadrature, _LIBRARY_N2N)
        solver = SNSolver(sn_mesh)
        phi = np.asarray(
            solve_sn(_LIBRARY_N2N, mesh, quadrature).scalar_flux.values,
            dtype=float,
        )
        total = solver.compute_production_rate(phi)
        fission_only = float(
            IntegratedReactionRate(
                sn_mesh.mat_xs.fission_production_field
            ).evaluate(phi)
        )
        relative = abs(total - fission_only) / abs(total)
        _require(
            relative > 1e-2,
            f"the (n,2n)-carrying fixture no longer separates the two "
            f"functionals (rel {relative:.4e} ≤ 1e-2) — the gauge rows that "
            f"cite it have lost their teeth (recorded: 1.2186e-01).",
        )


class TestRecordNothingRecordsWhichGaugeApplied:
    r"""Two normalisation conventions, one hub, no record.

    ``eigenvalue.py:437`` renormalises to unit production rate **only** when
    the solver ``isinstance``-conforms to
    :class:`~orpheus.numerics.eigenvalue.ProductionRateSolver`; a
    non-conforming solver keeps the legacy un-normalised trajectory.  ``[M]``
    ``probes/p6``-family measurement on the subcritical slab: the two
    conventions agree on k to **1.03e-09** and return fluxes whose sums differ
    by a factor **0.408356**.  Nothing on either answer says which ran.

    ⚠ **The brief's stronger form is NOT constructible through** ``solve_sn``:
    that entry always builds an :class:`~orpheus.sn.solver.SNSolver`, which IS
    a ``ProductionRateSolver``, so no pair of *Solutions* can differ by gauge
    convention.  The constructible pair is at the ``power_iteration`` tier, and
    the gauge-value distinguishability law belongs with ``ScaleGauge`` in
    ``tests/numerics/test_gauge.py`` (U1).  This row records the SN half.
    """

    class _NoProductionRate:
        """Delegate everything EXCEPT ``compute_production_rate`` — the legacy
        (non-``ProductionRateSolver``) convention the driver falls back to."""

        def __init__(self, inner: object) -> None:
            object.__setattr__(self, "_inner_solver", inner)

        def __getattr__(self, name: str) -> object:
            if name == "compute_production_rate":
                raise AttributeError(name)
            return getattr(self._inner_solver, name)

    def test_record_two_conventions_one_hub_no_record(self) -> None:
        """RECORD — the conventions are distinguishable by the FLUX and by
        nothing the answer carries."""
        materials, mesh, quadrature = _slab()
        sn_mesh = _as_sn_mesh(mesh, quadrature, materials)
        kwargs = dict(
            inner_solver="source_iteration", keff_tol=1e-10, flux_tol=1e-9,
        )
        gauged = SNSolver(sn_mesh, **kwargs)  # type: ignore[arg-type]
        legacy = self._NoProductionRate(SNSolver(sn_mesh, **kwargs))  # type: ignore[arg-type]
        _require(
            isinstance(gauged, ProductionRateSolver),
            "the SN solver stopped conforming to ProductionRateSolver — the "
            "gauged leg of this row is measuring the fallback twice.",
        )
        _require(
            not isinstance(legacy, ProductionRateSolver),
            "the surrogate still conforms — the two legs are the SAME "
            "convention and this row proves nothing (vv #17's AIM trap).",
        )
        a = power_iteration(gauged, max_iter=500, budget_name="max_outer")
        b = power_iteration(legacy, max_iter=500, budget_name="max_outer")  # type: ignore[arg-type]
        phi_a = np.asarray(a.flux_distribution, dtype=float)
        phi_b = np.asarray(b.flux_distribution, dtype=float)
        _require(
            abs(float(a.keff) - float(b.keff)) < 1e-7,
            f"non-vacuity: the two conventions must agree on k (the gauge is "
            f"a scale, k is scale-invariant); got {a.keff!r} vs {b.keff!r}.",
        )
        ratio = float(phi_b.sum()) / float(phi_a.sum())
        _require(
            abs(ratio - 1.0) > 0.1,
            f"the two conventions now return the SAME scale (ratio "
            f"{ratio:.6f}) — the record's discriminator is gone (recorded: "
            f"0.408356).",
        )
        # And the whole point: neither ANSWER says which functional ran.
        recorded = {
            name for outcome in (a, b)
            for name in dir(outcome) if "gauge" in name or "production" in name
        }
        _require(
            recorded == set(),
            f"a gauge/production member appeared on the power-iteration "
            f"outcome ({sorted(recorded)}) — if step 3 landed, DELETE this row.",
        )


class TestRecordTheAdjointArmsAreAsymmetric:
    r"""The seedless adjoint iterates on a BARE carrier, the carrying one on
    the coupled carrier — #467's subject.

    ``[M]`` ``probes/p2_adjoint_1x1_lift.py`` and ``probes/p3_typed_residual.py``:
    ``_adjoint_posing_parts`` returns a ``FullField`` template on the slab and
    a ``CoupledField`` template on the carrying sphere, while the hub's own
    ``system.space`` is a ONE-system ``CoupledSpace`` in both cases.
    """

    def test_record_the_two_arms_return_different_carrier_types(self) -> None:
        """RECORD — the arm asymmetry, read off the shipped helper."""
        slab_materials, slab_mesh, quadrature = _slab()
        sphere_materials, sphere_mesh, _ = _carrying_sphere()
        seedless = _as_sn_mesh(slab_mesh, quadrature, slab_materials)
        carrying = _as_sn_mesh(sphere_mesh, quadrature, sphere_materials)
        _require(
            seedless.radial_characteristic_field_space is None,
            "the slab hub started carrying System B — the arm split moved.",
        )
        _require(
            carrying.radial_characteristic_field_space is not None,
            "the sphere hub stopped carrying System B — this row's carrying "
            "arm is measuring a seedless mesh.",
        )
        _, _, _, seedless_template = _adjoint_posing_parts(seedless)
        _, _, _, carrying_template = _adjoint_posing_parts(carrying)
        _require(
            isinstance(seedless_template, FullField)
            and not isinstance(seedless_template, CoupledField),
            f"the seedless adjoint template is "
            f"{type(seedless_template).__name__}, not a bare FullField — if "
            f"U5 landed the 1×1 lift, DELETE this row.",
        )
        _require(
            isinstance(carrying_template, CoupledField),
            f"the carrying adjoint template is "
            f"{type(carrying_template).__name__}, not a CoupledField.",
        )
        _require(
            len(seedless.system.space._require_systems()) == 1,
            "the seedless hub's own system space is a ONE-system "
            "CoupledSpace — that is what makes the 1×1 lift spellable.",
        )


class TestRecordTheEigenScalarFluxIsNotTheAngularIntegral:
    r"""The eigen exit packages the POWER ITERATION's φ, not :math:`\int\psi
    d\Omega` of the ψ it returns.

    ``[M]`` ``probes/p1_scalar_flux_vs_integral.py`` over the 16 finalize
    cases: ``array_equal`` **0 of 16**, worst relative gap **7.4094e-11** (the
    ``cart2d_L0`` arm), all well inside the pins' band
    (``SAFETY(10) × flux_tol = 1e-8``) — so U4's derivation is a principled
    re-read whose pins hold, and whose ``DriftWarning`` tripwire will fire on
    every one of them.

    This row runs TWO of the sixteen (a slab and a CARRYING sphere — the
    marginal-axes and ray guards differ) so the record is cheap and the
    campaign's own probe carries the full table.
    """

    @pytest.mark.parametrize(
        "fixture", ["slab", "carrying_sphere"],
    )
    def test_record_the_two_flux_members_disagree(self, fixture: str) -> None:
        """RECORD — the two members are NOT the same object, and the gap is
        the convergence residual, not a bug."""
        materials, mesh, quadrature = (
            _slab() if fixture == "slab" else _carrying_sphere()
        )
        solution = solve_sn(
            materials, mesh, quadrature,
            keff_tol=1e-10, flux_tol=1e-9, inner_tol=1e-11,
        )
        stored = np.asarray(solution.scalar_flux.values, dtype=float)
        interior = solution.angular_flux.interior
        _require(
            isinstance(interior, AngularFlux),
            f"the eigen carrier's interior is {type(interior).__name__}, not "
            f"an AngularFlux — ∫ψ dΩ is not spellable on it and U4's "
            f"derivation would need a different reduction.",
        )
        assert isinstance(interior, AngularFlux)  # narrowing for pyright
        derived = np.asarray(
            interior.integrate_angular().values, dtype=float,
        )
        _require(
            stored.shape == derived.shape,
            f"the two members' shapes diverged: {stored.shape} vs "
            f"{derived.shape} — the derivation is not a re-read.",
        )
        _require(
            not np.array_equal(stored, derived),
            "the stored φ became bit-identical to ∫ψ dΩ — if U4 landed the "
            "derivation, DELETE this row.",
        )
        relative = float(
            np.max(np.abs(stored - derived) / np.maximum(np.abs(stored), 1e-300))
        )
        _require(
            relative < 1e-8,
            f"the gap {relative:.4e} exceeds the pins' own band "
            f"(SAFETY×flux_tol = 1e-8) — U4's 'principled re-read' claim does "
            f"not hold on this fixture.",
        )


class TestRecordTheMultiplyingEntrysBalanceReadsTheLOSS:
    r"""⛔ The certificate's number MOVES on the fifth entry — the plan's
    *"NO reported number changes"* holds for four of five.

    Today the fixed-source arms hand ``_exit_balance_defect`` the hub's LOSS
    (``system.loss`` / its bare arm, ``solver.py:3792``), so the multiplying
    entry — whose own equation is :math:`(A - F)\psi = q` — reports the
    imbalance of a DIFFERENT equation.  Step 3 records
    ``hub.source_posing(q)`` on that entry, whose operator is
    ``pencil.at(1.0) = A − F``.

    ``[M]`` ``probes/p6_anchor_facts.py`` on the truncated subcritical slab:

    * today, against :math:`A`   — **0.8294593510371534** (bit-identical to
      the number ``IterationHistory.balance_defect`` actually carries);
    * the design, against :math:`A - F` — **0.8758249879057027** (+5.59 %).

    The design's number is the HONEST one; the point of the row is that the
    change is real, predicted, and must be gated rather than asserted away.
    """

    def test_record_the_reported_defect_is_the_loss_imbalance(self) -> None:
        """RECORD — reproduce the reported number from ``A`` alone, and show
        ``A − F`` gives a different one on the same iterate."""
        solution = _multiplying_solution_truncated()
        history = solution.history
        _require(history is not None, "the entry stopped carrying a history.")
        assert history is not None  # narrowing for the type checker
        reported = history.balance_defect
        _require(
            reported is not None,
            "non-vacuity: the fixture must be TRUNCATED, else the defect is "
            "None by design and this row compares two Nones.",
        )
        materials, mesh, quadrature = _slab()
        sn_mesh = _as_sn_mesh(mesh, quadrature, materials)
        source = _build_fixed_source_rhs(
            _uniform_source(quadrature, 2, 8), sn_mesh,
        )
        state = CoupledField(systems=(solution.angular_flux,))
        source_rate = _balance_projection(source, sn_mesh=sn_mesh)
        denominator = float(np.linalg.norm(np.asarray(source_rate)))
        _require(denominator > 0.0, "the source rate integrates to zero.")

        def defect(operator) -> float:
            applied = operator.apply(state).systems[0]
            per_group = (
                _balance_projection(applied, sn_mesh=sn_mesh) - source_rate
            )
            return float(np.linalg.norm(np.asarray(per_group))) / denominator

        against_loss = defect(sn_mesh.pencil.at(0.0))
        against_posing = defect(sn_mesh.pencil.at(1.0))
        _require(
            against_loss == reported,
            f"the reported defect {reported!r} is no longer the LOSS "
            f"imbalance {against_loss!r} — the instrument this row reproduces "
            f"has moved (recorded: 0.8294593510371534).",
        )
        _require(
            against_posing != against_loss,
            f"non-vacuity: the fixture's fission term vanished, so A and A−F "
            f"agree ({against_loss!r}) and this row cannot see the change. "
            f"Use a SUBCRITICAL fissile hub with k near 1.",
        )
        shift = abs(against_posing - against_loss) / against_loss
        _require(
            shift > 1e-3,
            f"the A vs A−F shift collapsed to {shift:.4e}; recorded 5.59e-02 "
            f"at k = 0.907. A near-critical hub is what makes this visible.",
        )


class TestRecordThePureTransportPosingIsTheLOSS:
    r"""F11's identity read, and the Pattern-2 reconciliation, both measured.

    ``[M]`` ``probes/p3_typed_residual.py``: ``pencil.at(0.0)`` **is**
    ``pencil.lhs`` **is** ``system.loss`` (object identity — ``pencil.py:69``),
    so *"was fission suppressed?"* is an ``is`` test and needs no datum; and
    ``SourcePosing(pencil.at(0), q).residual(ψ)`` is **bit-exactly**
    ``−evaluate_residual(system, ψ, q)`` on BOTH arms (``array_equal`` True,
    ``max|Δ| = 0.000000e+00``), so the sign flip U2c lands is exact.

    ⚠ The same claim is FALSE for ``hub.source_posing(q)`` (= ``at(1.0)``):
    ``[M]`` the two differ by :math:`|F\psi| = 2.25\mathrm{e}{-1}` on the
    carrying sphere — see the sibling record above.
    """

    @pytest.mark.parametrize("fixture", ["slab", "carrying_sphere"])
    def test_record_at_zero_is_the_loss_by_IDENTITY(self, fixture: str) -> None:
        """RECORD — identity, not equality; the F11 read is free."""
        materials, mesh, quadrature = (
            _slab() if fixture == "slab" else _carrying_sphere()
        )
        sn_mesh = _as_sn_mesh(mesh, quadrature, materials)
        _require(
            sn_mesh.pencil.at(0.0) is sn_mesh.pencil.lhs,
            "pencil.at(0.0) stopped returning lhs itself — F11's identity "
            "read ('was fission suppressed?') needs a datum again.",
        )
        _require(
            sn_mesh.pencil.lhs is sn_mesh.system.loss,
            "the pencil's lhs is no longer the record's loss object.",
        )

    @pytest.mark.parametrize("fixture", ["slab", "carrying_sphere"])
    def test_record_the_residual_sign_twin_is_EXACT(self, fixture: str) -> None:
        """RECORD — ``q − Aψ`` and ``Aψ − q`` agree to the bit under the
        PURE-TRANSPORT posing (the one U2 flips the sign of)."""
        from orpheus.numerics.posing import SourcePosing
        from orpheus.sn.solver import evaluate_residual

        materials, mesh, quadrature = (
            _slab() if fixture == "slab" else _carrying_sphere()
        )
        sn_mesh = _as_sn_mesh(mesh, quadrature, materials)
        ng, nx = 2, len(mesh.mat_ids)
        source = _build_fixed_source_rhs(
            _uniform_source(quadrature, ng, nx), sn_mesh,
        )
        lifted_source = (
            source if isinstance(source, CoupledField)
            else CoupledField(systems=(source,))
        )
        template = sn_mesh.system.space.zeros()
        state = CoupledField.from_flat(
            np.ones(template.to_flat().size), template,
        )
        posing = SourcePosing(sn_mesh.pencil.at(0.0), lifted_source)
        mine = np.asarray(posing.residual(state).to_flat(), dtype=float)
        # ``evaluate_residual``'s arity guard: the SEEDLESS system refuses a
        # coupled wrapper (`[M]` "this system is 1×1 (seedless) — pass the
        # bare FullField pair").  That guard is a §6b member of U4.
        carrying = sn_mesh.radial_characteristic_field_space is not None
        probe: FullField | CoupledField
        if carrying:
            probe = state
        else:
            member = state.systems[0]
            _require(
                isinstance(member, FullField),
                f"the seedless system's member is {type(member).__name__}, "
                f"not a FullField — evaluate_residual's 1×1 arity guard "
                f"cannot be fed.",
            )
            assert isinstance(member, FullField)  # narrowing for pyright
            probe = member
        theirs_field = evaluate_residual(sn_mesh.system, probe, source)
        theirs = np.asarray(theirs_field.to_flat(), dtype=float)
        _require(
            np.array_equal(-mine, theirs),
            f"the sign twin is no longer exact: max|(-r) - er| = "
            f"{np.max(np.abs(-mine - theirs)):.6e} (recorded 0.000000e+00).",
        )


# ═══════════════════════════════════════════════════════════════════════
# XFAIL(strict) — the ruled post-carve behaviour.  RED today, by design.
#
# Every body below is structured so EXACTLY ONE statement can fail and it is
# the documented one (``vv`` Mode-8 fourth class): preconditions are gathered
# into ``evidence`` and reported inside the message, never asserted.
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.xfail(
    strict=True,
    reason="U2 (F4): the kind becomes the OUTCOME's type; `keff` retires "
           "from the base and a source Solution has no `.keff` at all.",
)
class TestRuledTheKindIsTheOutcomesType:
    """§3.1 F4 / §3.2 — asserted over ``dataclasses.fields``, never over the
    name the carve picks for the outcome member."""

    def test_ruled_no_field_is_optional_by_kind(self) -> None:
        """The carrier's field set carries no ``| None`` — every Optional it
        has today encodes the kind, and the kind becomes a TYPE."""
        optional = sorted(
            f.name for f in dataclasses.fields(SolutionBase)
            if "None" in str(f.type)
        )
        _require(
            optional == [],
            f"SolutionBase still carries Optional-by-kind fields {optional}; "
            f"the kind is still read off a nullable datum (§3.1 F4).",
        )

    def test_ruled_a_source_solution_has_no_keff_attribute(self) -> None:
        """Not ``keff is None`` — **no attribute at all**, so the old
        discrimination is UNSPELLABLE rather than merely discouraged."""
        materials, mesh, quadrature = _slab()
        evidence: list[str] = []
        try:
            fixed = solve_sn_fixed_source(
                materials, mesh, quadrature,
                _uniform_source(quadrature, 2, 8),
            )
            evidence.append(f"type={type(fixed).__name__}")
            present = hasattr(fixed, "keff")
        except Exception as exc:  # noqa: BLE001 — evidence, not the verdict
            evidence.append(f"the fixed-source solve raised {exc!r}")
            present = True
        _require(
            not present,
            f"a fixed-source Solution still answers `.keff` "
            f"({'; '.join(evidence)}) — the kind is still a property.",
        )


@pytest.mark.xfail(
    strict=True,
    reason="U2 (F5): the eigen Solution records the SECTION that picked its "
           "representative — a functional and the target it was scaled to.",
)
class TestRuledTheEigenSolutionRecordsItsGauge:
    r"""§3.1 F5 — the gauge is an OBJECT (``ScaleGauge(functional, target)``,
    ``apply``/``displacement``), and the functional is *the object that ran*,
    not a label: ``[M]`` four functionals ship under the one word "production
    rate" and they differ by **1.2186e-01** on a Σ₂-carrying mixture.

    Asserted structurally — a member reachable from the Solution that exposes
    ``functional``, ``target`` and ``apply``, and whose section law
    ``functional(state) ≈ target`` holds on the RETURNED state.  The member's
    NAME is the carve's to choose; the scanner
    (:func:`_section_shaped_members`) is validated against a planted member by
    :class:`TestFilterTheScannersFindAPlantedMember`.
    """

    def test_ruled_the_section_is_recoverable_and_its_law_holds(self) -> None:
        """One statement: a section exists whose law closes on the state."""
        materials, mesh, quadrature = _slab()
        evidence: list[str] = []
        holds = False
        try:
            solution = solve_sn(materials, mesh, quadrature)
            sections = _section_shaped_members(solution)
            evidence.append(f"{len(sections)} section-shaped member(s)")
            for section in sections:
                state = getattr(
                    getattr(solution, "outcome", None), "state", None,
                )
                probe = state if state is not None else solution.scalar_flux
                value = float(section.functional(probe))
                evidence.append(
                    f"functional(state)={value!r} target={section.target!r}"
                )
                if np.isclose(value, section.target, rtol=1e-9, atol=0.0):
                    holds = True
        except Exception as exc:  # noqa: BLE001 — evidence, not the verdict
            evidence.append(f"raised {exc!r}")
        _require(
            holds,
            f"no recorded section reproduces its own normalisation on the "
            f"returned state ({'; '.join(evidence) or 'nothing found'}) — the "
            f"gauge that picked the representative is still unrecorded.",
        )


@pytest.mark.xfail(
    strict=True,
    reason="U2 (§3.3): the multiplying entry records the admissibility k it "
           "measured, as Certified(k, keff_tol).",
)
class TestRuledTheMultiplyingSolutionRecordsItsAdmissibility:
    """§3.1 / §3.3 — asserted by SEARCHING the answer for the number
    (:func:`_declared_numeric_data`, validated against a planted scalar), so
    the carve is free to name the certificate member."""

    def test_ruled_the_k_is_recoverable_from_the_answer(self) -> None:
        """One statement: some non-array datum on the Solution reproduces the
        hub's k."""
        k = _hub_k()
        solution = _multiplying_solution_truncated()
        seen = _declared_numeric_data(solution)
        recovered = [v for v in seen if abs(v - k) <= 1e-6 * abs(k)]
        _require(
            recovered,
            f"the admissibility k = {k!r} the entry MEASURED is nowhere on "
            f"its answer; the scanned numeric data was {sorted(set(seen))[:20]}"
            f" (… {len(seen)} values).",
        )


@pytest.mark.xfail(
    strict=True,
    reason="U5 (#467): both adjoint arms pose ONE object — the hub's own "
           "eigen_posing, daggered.",
)
class TestRuledBothAdjointArmsPoseTheHubsDaggeredQuestion:
    r"""§3.1 / §3.3 — a COUNTING SPY on
    :meth:`~orpheus.numerics.posing.EigenPosing.H`, because the claim is about
    the ROUTE (which object the entry poses), and no value functional states
    it: ``[M]`` ``probes/p2`` the three candidate routes agree on ``k_adj`` to
    **1.1e-15 / 1.6e-15** relative, so a value gate is a provable non-catcher.

    ⚠ ``k^† = k`` makes the k-equality rows Mode-12 blind to an UNDAGGERED
    posing (battery arm A12) — this row's instrument is the spy, not the
    number.
    """

    @pytest.mark.parametrize("fixture", ["slab", "carrying_sphere"])
    def test_ruled_the_entry_calls_EigenPosing_H(
        self, fixture: str, monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """One statement: the spy counted at least one dagger per solve."""
        import orpheus.numerics.posing as posing_module

        calls: list[str] = []
        original = posing_module.EigenPosing.H

        def spy(self):  # type: ignore[no-untyped-def]
            calls.append(type(self.pencil.lhs).__name__)
            return original(self)

        monkeypatch.setattr(posing_module.EigenPosing, "H", spy)
        _require(
            posing_module.EigenPosing.H is spy,
            "the spy did not bind — every count below is a false zero "
            "(lessons L46e).",
        )

        materials, mesh, quadrature = (
            _slab() if fixture == "slab" else _carrying_sphere()
        )
        evidence: list[str] = []
        try:
            adjoint = solve_sn_adjoint(
                materials, mesh, quadrature, keff_tol=1e-9, flux_tol=1e-8,
            )
            evidence.append(f"k_adj={adjoint.keff!r}")
        except Exception as exc:  # noqa: BLE001 — evidence, not the verdict
            evidence.append(f"the adjoint solve raised {exc!r}")
        _require(
            calls,
            f"solve_sn_adjoint posed its question without daggering the hub's "
            f"eigen_posing ({'; '.join(evidence)}) — the two arms still build "
            f"their own EigenPosing (#467).",
        )


@pytest.mark.xfail(
    strict=True,
    reason="U4 (F3): the state is stored WHOLE; the flux members and the ray "
           "member become derived, so the presence biconditional is deleted.",
)
class TestRuledTheStateIsStoredWhole:
    """§3.1 F3 — asserted over the FIELD set plus the descriptor type; the
    arity derivation itself is U4's own gate."""

    def test_ruled_the_flux_members_are_derived_not_stored(self) -> None:
        """``angular_flux``/``scalar_flux``/``radial_characteristic`` keep
        their reader NAMES (0 reader edits) and stop being fields."""
        names = frozenset(f.name for f in dataclasses.fields(SolutionBase))
        stored = sorted(
            names & {"angular_flux", "scalar_flux", "radial_characteristic"}
        )
        _require(
            stored == [],
            f"{stored} are still FIELDS of the carrier — the state is not "
            f"stored whole (§3.1 F3), so the presence biconditional and the "
            f"marginal-axes guard still have inputs to refuse.",
        )

    def test_ruled_the_ray_member_is_a_derived_reader(self) -> None:
        """The reader survives as a property so the 15 ``.radial_characteristic``
        call sites need no edit."""
        member = typing.cast(
            object, getattr(SolutionBase, "radial_characteristic", None),
        )
        _require(
            isinstance(member, property),
            f"SolutionBase.radial_characteristic is {type(member).__name__}, "
            f"not a property — the ray member is not yet derived from the "
            f"state's ARITY.",
        )
