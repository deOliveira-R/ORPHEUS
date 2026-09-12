r"""Consumers campaign step 1 — the PRE-CARVE anchors for the SN Problem's identity.

Landed on the UNMODIFIED tree, before the first production edit of the
identity step (plan ``.claude/plans/cs4c_binding_design.md`` §27, rulings
**R-cc3 / R-cc4 / R-cc8 / R-cc9**; GitHub **#459**), so every row here is a
measurement of the tree the carve starts from.

**The two ruled predicates (R-cc8), and why they are two.**

* ``same_phase_space`` — CONTRACTIBILITY. *May two solutions' FIELDS be
  paired?* Geometry × quadrature × layout × scheme type, by CONTENT. Serves
  ``SolutionBase.compare``, ``Solution.homogenize(adjoint=…)`` and
  ``Solution.condense(adjoint=…)``. The closure class and the truncation
  order are EXCLUDED — and that exclusion is measured, not assumed:
  ``[M]`` 2026-09-12 the returned ``angular_flux`` is the ORDINATE carrier
  at every order on both arms (1-D ``TimedFullField(4, 2, 8)`` at L = 0 and
  L = 1; 2-D windowed ``TimedFullField(24, 2, 4, 4)`` at both), so no field
  layout moves with ``L``.
* ``__eq__`` / ``__hash__`` — full Problem IDENTITY. *Is this the same
  problem?* Every generating datum by content, closure class and truncation
  order INCLUDED. Serves save-state and ``Solution`` provenance. A P0
  forward and a P3 adjoint share a phase space and are DIFFERENT problems.

**Three claim kinds live here with different fates — read the class
docstring before touching a row.**

* ``TestTodays…`` classes — **RECORD** of states the carve DELETES or
  FLIPS. Green today, designed to RED at the carve, deleted (not repaired)
  in the carve's own commit. They exist so the API change cannot land
  quietly: an ``xfail`` row is SILENT when the new code lands wrong, and a
  RECORD row is not (``vv`` Mode 8, the misattributed-strict-xfail class).
* ``xfail(strict=True)`` classes — the RULED post-carve gates, a
  self-retiring todo list. Their XPASS is a failure, so the carve's commit
  must delete the marker.
* :class:`TestTheClosureExclusionSurvives` — **MUST STAY GREEN**: the one
  place R-cc8 deliberately keeps today's behaviour, so a carve that
  "strengthens" ``same_phase_space`` by folding the closure in reds here.

⚠ **Two activation facts every order row asserts in-test.** ``[M]`` all
twelve ``xs_library`` mixtures ship ``len(SigS) == 2``, so the clamp
``min(L, len(SigS) − 1)`` maps **every** request ≥ 1 to **1**: a row spelled
"P2 vs P3" asserts ``1 != 1`` and is a false green. The only non-vacuous
abstract-library pair is ``(0, ≥ 1)``. And ``[M]`` the P1 block is non-zero
there — ``solve_sn`` reads ``k = 1.2122522010124397`` at L = 0 against
``1.2180192347287149`` at L ≥ 1 on this module's two-region slab, a 476 pcm
separation, so the channel is live rather than merely present.
"""

from __future__ import annotations

import inspect

import numpy as np
import pytest

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry.coord import CoordSystem
from orpheus.geometry.mesh import BC, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.angular.closure import IdentityAngularClosure
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.solution import Solution, SolutionBase
from orpheus.transport.mesh.axis import AxisMesh
from orpheus.transport.spatial.linear_discontinuous import LinearDiscontinuous

pytestmark = [pytest.mark.foundation, pytest.mark.catches("ERR-084")]

_RULING_ID = "#459 / R-cc8 — same_phase_space (contractibility) and __eq__ (identity)"
_RULING_ORDER = "#459 / R-cc9 — the clamped truncation order is a hub constructor datum"


def _require(cond: object, msg: str) -> None:
    """``-O``-safe assertion (``coding-standards``: the canonical runner strips ``assert``)."""
    if not cond:
        raise AssertionError(msg)


# ── fixtures: every call builds FRESH constituent objects ─────────────


def _mats(region_b: str = "B") -> dict[int, Mixture]:
    """FRESH ``Mixture`` objects — equal data, distinct identity."""
    return {0: get_mixture("A", "2g"), 1: get_mixture(region_b, "2g")}


#: SHARED constituents. Rows that isolate the SPATIAL leg must hold the
#: quadrature and the materials fixed by ``is`` identity — otherwise today's
#: per-mixture ``is`` leg refuses the pair for a reason the row is not about,
#: and the RECORD reads ``False`` for the wrong cause. ⛔ I shipped this file
#: once with ``_mats()`` called per hub and two rows failed exactly that way.
_SHARED_QUAD = Quadrature.level_symmetric(sn_order=4)
_SHARED_MATS = _mats()


def _axes(d: int, cells=(3, 4, 5), extents=(1.0, 2.0, 3.0), bc=(None, None)):
    return tuple(
        AxisMesh(edges=np.linspace(0.0, e, n + 1), bc_low=bc[0], bc_high=bc[1])
        for e, n in list(zip(extents, cells))[:d]
    )


def _hub(d: int = 3, *, cells=(3, 4, 5), extents=(1.0, 2.0, 3.0), bc=(None, None),
         quad=None, mats=None, scheme=None, mat_map=None) -> SNMesh:
    """A hub over the SHARED quadrature and materials — the spatial leg isolated."""
    return SNMesh.from_axes(
        _axes(d, cells, extents, bc),
        _SHARED_QUAD if quad is None else quad,
        _SHARED_MATS if mats is None else mats,
        mat_map=mat_map,
        scheme=scheme,
    )


def _hub_independent(d: int = 3, **kw) -> SNMesh:
    """A hub sharing NOTHING by identity — the honest CONTENT positive control.

    A content key that silently fell back to identity would pass a positive
    control built over shared constituents and fail here, which is the whole
    reason the two helpers are separate (``vv`` #19 — the reading that
    cannot change).
    """
    return _hub(d, quad=Quadrature.level_symmetric(sn_order=4), mats=_mats(), **kw)


def _slab_two_region() -> tuple[dict[int, Mixture], Mesh1D, Quadrature]:
    """The two-region 2-group slab the order rows solve on."""
    mesh = Mesh1D(edges=np.linspace(0.0, 2.0, 9),
                  mat_ids=np.array([0] * 4 + [1] * 4, dtype=int))
    return _mats(), mesh, Quadrature.gauss_legendre(4)


#: The closure pair's shared geometry and quadrature — shared by ``is`` so the
#: pair isolates the CLOSURE and nothing else.
_SPHERE_MESH = Mesh1D(edges=np.linspace(0.0, 1.0, 5), mat_ids=np.zeros(4, dtype=int),
                      coord=CoordSystem.SPHERICAL, bc_right=BC("vacuum"))
_SPHERE_QUAD = Quadrature.gauss_legendre(4)


def _sphere(closure=None) -> SNMesh:
    """A SPHERICAL hub — the only chart where the closure override CONSTRUCTS.

    ``[M]`` 2026-09-12: a slab REFUSES ``MorelMontryAngularSweep``
    (``ValueError: march_start_structure_per_level supports SPHERICAL or
    CYLINDRICAL …``), so the closure witness must be built the other way
    round — a sphere carrying the Cartesian ``IdentityAngularClosure``.
    """
    return SNMesh(_SPHERE_MESH, _SPHERE_QUAD, _SHARED_MATS, angular_closure=closure)


class TestTheClosureExclusionSurvives:
    r"""R-cc8 keeps the contractibility ruling FOR ``same_phase_space``.

    The predicate's own docstring (``augmented_mesh.py:571-577``) and
    ``docs/theory/methods/sn/index.rst:849-856`` say *"do not strengthen
    that predicate by adding the closure"* — and R-cc8 agrees, for the
    CONTRACTIBILITY predicate. This class is the gate that makes a
    well-meaning "fold the closure in" red.

    ``[M]`` 2026-09-12 the witness CONSTRUCTS: a spherical hub whose
    default closure is ``MorelMontryAngularSweep`` accepts
    ``angular_closure=IdentityAngularClosure``, and today's predicate reads
    ``True`` across the pair. ⚠ The reverse pair is UNCONSTRUCTIBLE — a
    slab refuses ``MorelMontryAngularSweep`` outright — so this chart is
    not a preference, it is the only one available.
    """

    def test_two_closures_over_one_geometry_share_the_phase_space(self) -> None:
        default, overridden = _sphere(), _sphere(closure=IdentityAngularClosure)
        _require(
            type(default.angular_closure) is not type(overridden.angular_closure),
            "activation: the two hubs must carry DIFFERENT closure classes",
        )
        _require(
            default.same_phase_space(overridden),
            "R-cc8: two closures over one geometry remain contractible",
        )


class TestCrossClassComparisonIsAlreadySafe:
    """MUST STAY GREEN — a property the carve must PRESERVE, not create.

    ⛔ I shipped this as an ``xfail`` row and the harness refuted it:
    ``[M]`` 2026-09-12 it ``XPASS(strict)``-ed on the unmodified tree.
    ``SNMesh`` inherits ``object.__eq__``, so a foreign comparison is
    ``False`` and never raises.

    ⟹ not a gap the carve closes but a REGRESSION PIN on the hand-written
    ``__eq__`` that replaces it. It is load-bearing because ``Solution`` is a
    frozen dataclass with a GENERATED ``__eq__`` over a field tuple that
    INCLUDES ``mesh`` (``solution.py:346/667/1082``, no hand-written
    ``__eq__``): a hub whose ``__eq__`` raises on a foreign operand would
    surface at a ``Solution == Solution`` call site that never mentions
    ``SNMesh``.
    """

    def test_comparison_across_classes_never_raises(self) -> None:
        a = _hub(3)
        _require(not (a == object()), "a hub is not an arbitrary object")
        _require(not (a == None), "a hub is not None")  # noqa: E711
        _require(a != object(), "the negated form must agree")


# ════════════════════════════════════════════════════════════════════
# The RULED post-carve gates — strict xfail, a self-retiring todo list
# ════════════════════════════════════════════════════════════════════


class TestSamePhaseSpaceIsContractibilityByContent:
    """R-cc8, predicate 1 — ``same_phase_space`` over CONTENT."""

    @pytest.mark.parametrize("d", [1, 2, 3])
    def test_same_data_two_constructions_pair(self, d: int) -> None:
        """POSITIVE CONTROL at every rank, over INDEPENDENT constituents."""
        a, b = _hub_independent(d), _hub_independent(d)
        _require(a is not b, "activation: two distinct hubs")
        _require(a.same_phase_space(b), "equal generating data ⟹ one phase space")

    @pytest.mark.parametrize(
        "label,other",
        [
            ("edges", dict(extents=(1.0, 2.0, 9.0))),
            ("cell_count", dict(cells=(3, 4, 6))),
            ("bc_tag", dict(bc=(BC("vacuum"), BC("reflective")))),
            ("quadrature", dict(quad=Quadrature.level_symmetric(sn_order=6))),
            ("materials", dict(mats={0: get_mixture("A", "2g"), 1: get_mixture("C", "2g")})),
            ("scheme_type", dict(scheme=LinearDiscontinuous())),
        ],
    )
    def test_one_moved_datum_refuses(self, label: str, other: dict) -> None:
        """The per-datum NEGATIVE legs — one flip per contractibility datum."""
        a = _hub(3)
        b = _hub(3, **other)
        _require(
            not a.same_phase_space(b),
            f"moving {label} must make the two phase spaces different",
        )

    def test_the_material_MAP_is_a_contractibility_datum(self) -> None:
        """``mat_map`` — the datum with no home on an ``AxisMesh``.

        It rides the ``from_axes`` keyword, so it is invisible to any key
        built by walking ``self.axes`` alone.
        """
        shape = (3, 4, 5)
        flat = np.zeros(shape, dtype=int)
        other = flat.copy()
        other[0, 0, 0] = 1
        a = _hub(3, mat_map=flat)
        b = _hub(3, mat_map=other)
        _require(not a.same_phase_space(b), "a different material assignment is a different space")

    def test_the_truncation_order_is_NOT_a_contractibility_datum(self) -> None:
        r"""R-cc8's insensitivity leg — the one the fields' shapes license.

        ``[M]`` the returned ``angular_flux`` is the ORDINATE carrier at
        every order (1-D ``(4, 2, 8)`` at L = 0 and L = 1; 2-D windowed
        ``(24, 2, 4, 4)`` at both), so a P0 forward and a P3 adjoint pair
        and ``condense(adjoint=…)`` must proceed.
        """
        mats, mesh, quad = _slab_two_region()
        a = SNMesh(mesh, quad, mats, scattering_order=0) 
        b = SNMesh(mesh, quad, mats, scattering_order=3) 
        _require(a.scattering_order != b.scattering_order,
                 "activation: the two hubs must retain DIFFERENT orders")
        _require(a.same_phase_space(b), "R-cc8: the order does not move the field layout")


class TestProblemIdentityIsEveryGeneratingDatum:
    """R-cc8, predicate 2 — ``__eq__`` / ``__hash__`` over the FULL data."""

    @pytest.mark.parametrize("d", [1, 2, 3])
    def test_same_data_compares_equal_and_hashes_equal(self, d: int) -> None:
        """POSITIVE CONTROL over INDEPENDENT constituents — nothing shared by ``is``."""
        a, b = _hub_independent(d), _hub_independent(d)
        _require(a is not b, "activation: two distinct hubs")
        _require(a == b, "equal generating data ⟹ the same problem")
        _require(hash(a) == hash(b), "equal problems must hash equal")
        _require(len({a, b}) == 1, "a set of two equal problems holds one member")

    def test_the_closure_class_IS_an_identity_datum(self) -> None:
        r"""The fork's discriminating row — ``same_phase_space`` True, ``==`` False.

        This pair is the ONLY place the two ruled predicates disagree on a
        constructible input, so it is what makes "two predicates" a
        measurable design rather than a naming choice.
        """
        default, overridden = _sphere(), _sphere(closure=IdentityAngularClosure)
        _require(
            type(default.angular_closure) is not type(overridden.angular_closure),
            "activation: the two hubs must carry DIFFERENT closure classes",
        )
        _require(default.same_phase_space(overridden), "they still contract (R-cc8)")
        _require(not (default == overridden), "R-cc8: the closure is generating data")

    def test_the_truncation_order_IS_an_identity_datum(self) -> None:
        r"""R-cc4's headline: a P0 forward and a P3 adjoint are DIFFERENT problems.

        ⚠ The pair is ``(0, 3)`` and NOT ``(2, 3)`` for a measured reason:
        ``[M]`` every ``xs_library`` mixture has ``len(SigS) == 2``, so the
        clamp maps both 2 and 3 to **1** and a ``(2, 3)`` row would assert
        ``1 != 1`` — a false green. The clamped values are asserted below so
        the row cannot decay into one.
        """
        mats, mesh, quad = _slab_two_region()
        a = SNMesh(mesh, quad, mats, scattering_order=0) 
        b = SNMesh(mesh, quad, mats, scattering_order=3) 
        _require(a.scattering_order == 0, "activation: the P0 hub retains 0") 
        _require(b.scattering_order == 1, "activation: P3 CLAMPS to 1 on this library")
        _require(not (a == b), "R-cc4: two retained orders are two problems")

    def test_identity_is_strictly_finer_than_contractibility(self) -> None:
        """The law relating the two ruled predicates (``vv`` #15).

        ``a == b`` ⟹ ``a.same_phase_space(b)``, and the containment is
        STRICT — the closure pair witnesses a contractible pair that is not
        one problem. A carve that got the containment backwards passes both
        per-datum tables and fails here.

        ⚠ The non-vacuity guard is load-bearing and it is what makes this a
        gate: ``[M]`` 2026-09-12 the row PASSES on the unmodified tree
        without it, because ``SNMesh.__eq__`` is ``object.__eq__`` and no
        pair is ever ``==``, so the implication holds over an empty set. I
        shipped it that way once and the harness reported ``XPASS(strict)``.
        """
        equal_pairs = [
            (_hub_independent(3), _hub_independent(3)),
            (_sphere(), _sphere()),
        ]
        strict_pairs = [(_sphere(), _sphere(closure=IdentityAngularClosure))]
        _require(
            all(a == b for a, b in equal_pairs),
            "NON-VACUITY: the implication needs at least one EQUAL pair to range over",
        )
        for a, b in equal_pairs + strict_pairs:
            if a == b:
                _require(
                    a.same_phase_space(b),
                    "identity must imply contractibility on every pair",
                )
        for a, b in strict_pairs:
            _require(
                a.same_phase_space(b) and not (a == b),
                "STRICTNESS: a contractible pair that is not one problem must exist",
            )



class TestTheHubOwnsTheClampedOrder:
    """R-cc9 — one clamp, at construction, read by every entry."""

    @pytest.mark.parametrize("requested,retained", [(0, 0), (1, 1), (3, 1), (5, 1)])
    def test_the_clamp_runs_once_at_construction(self, requested: int, retained: int) -> None:
        r"""``min(L, len(SigS) − 1)``, evaluated on the hub.

        ``[M]`` the abstract library ships ``len(SigS) == 2`` on all twelve
        mixtures, so the retained order saturates at 1 — the row's
        expectation table IS that measurement, and the activation leg
        asserts the premise so a library change reds here rather than
        silently making the table meaningless.
        """
        mats, mesh, quad = _slab_two_region()
        _require(
            min(len(m.SigS) for m in mats.values()) - 1 == 1,
            "activation: this expectation table assumes a P1 library",
        )
        hub = SNMesh(mesh, quad, mats, scattering_order=requested)
        _require(hub.scattering_order == retained,
                 f"request {requested} must retain {retained}")

    def test_the_adjoint_entry_never_exceeds_the_hubs_order(self) -> None:
        r"""The ROUTE claim (Mode 11), and the row that measures R-cc9's whole point.

        Today ``solve_sn_adjoint(scattering_order=3)`` reaches
        ``TransferKernel.at_order(3)`` on a P1 library (``[M]`` above). After
        the carve the adjoint posing reads the HUB's clamped value, so the
        largest order the kernel ever sees is the hub's.

        ⭐ **What this gate measures, stated because a battery arm measured
        it.** The spy observes the order the CALLER PASSES, not the order the
        kernel serves. ``[M]`` 2026-09-12, battery arm ``A5`` — a clamp
        pushed DOWN into ``at_order`` itself (bite confirmed: ``at_order(3)``
        honest → order 3, mutant → order 1) leaves this row and its RECORD
        sibling **0 reds**. That is the CORRECT refusal, not a blindness: a
        clamp inside the kernel would be a FOURTH spelling of "the retained
        order", which is the defect R-cc9 exists to remove. The gate is
        keyed on the caller because R-cc9's claim is about the caller.
        """
        import orpheus.transport.kernels as kernels
        from orpheus.sn.solver import solve_sn_adjoint

        mats, mesh, quad = _slab_two_region()
        seen: set[int] = set()
        original = kernels.TransferKernel.at_order

        def _spy(self, order):  # noqa: ANN001, ANN202
            seen.add(int(order))
            return original(self, order)

        kernels.TransferKernel.at_order = _spy
        try:
            solve_sn_adjoint(mats, mesh, quad, scattering_order=3,
                             keff_tol=1e-9, inner_tol=1e-10, max_outer=400)
        finally:
            kernels.TransferKernel.at_order = original
        _require(seen, "activation: the spy must observe at least one at_order call")
        _require(
            max(seen) <= 1,
            f"the adjoint must never exceed the hub's clamped order; saw {sorted(seen)}",
        )


class TestTheInternedGeometryCacheUnderContentIdentity:
    """Open ruling **O-3** — what happens to ``_GEOM_CACHE_INTERN``.

    ``[M]`` 2026-09-12, two content-equal, distinct, simultaneously-live
    hubs and six alternating ``geometry_cache_for`` calls:

    ================================  ======  ============
    regime                            builds  intern size
    ================================  ======  ============
    identity keys (today)             2       2
    simulated content ``__eq__``      **6**   **1**
    ================================  ======  ============

    Mechanism: the intern validates ``entry[0] is angular_closure`` and each
    hub binds its OWN closure instance (``cls(...)``), so two content-equal
    hubs ping-pong, each overwriting the other's Stratum-1 table. The
    committed cache gates cannot see it: ``test_cache.py:303`` counts
    ``CollisionCache._build_count`` (Stratum **2**) and ``:361-369`` asserts
    identity on ONE mesh.

    The fix is a WIN rather than a patch — ``[M]`` the table's eight fields
    are bare arrays with no mesh and no closure reference, and it is σ-free,
    so two content-equal hubs may legitimately SHARE one entry once the
    validation is re-keyed from the closure INSTANCE to its CLASS; the count
    then goes 2 → 1. The RECORD row below is the pre-carve reading; the
    xfail row states the ruled outcome.
    """

    @staticmethod
    def _count_builds(a: SNMesh, b: SNMesh, rounds: int = 3) -> int:
        import orpheus.sn.loss_representation as loss_representation

        loss_representation._GEOM_CACHE_INTERN.clear()
        builds: list[int] = []
        cache_cls = loss_representation.StreamingCoefficientCache
        original = cache_cls.from_mesh_and_quad

        def _spy(mesh: SNMesh):  # noqa: ANN202
            builds.append(id(mesh))
            return original(mesh)

        # ``vv`` #29 / L77e: rebind the BOUND classmethod on the class object, so
        # the wrapper is what every caller resolves. A monkeypatch-only battery
        # is crash-safe by construction (L63h) — nothing on disk to restore.
        setattr(cache_cls, "from_mesh_and_quad", _spy)
        try:
            for _ in range(rounds):
                loss_representation.geometry_cache_for(a, a.angular_closure)
                loss_representation.geometry_cache_for(b, b.angular_closure)
        finally:
            setattr(cache_cls, "from_mesh_and_quad", original)
        _require(builds, "POSITIVE CONTROL: the spy must observe at least one build")
        return len(builds)

    def test_content_equal_hubs_share_ONE_table(self) -> None:
        """O-3 RULED (2026-09-12): a content key + closure-CLASS validation ⟹ two
        content-equal live hubs share ONE Stratum-1 table. (The pre-carve RECORD
        read 2 — one table per hub under identity keys; the rejected instance
        validation would have read 6, the ping-pong.)"""
        mats, mesh, quad = _slab_two_region()
        a, b = SNMesh(mesh, quad, mats), SNMesh(mesh, quad, mats)
        _require(a is not b and a == b, "activation: two distinct, content-equal, live hubs")
        _require(
            self._count_builds(a, b) == 1,
            "content-equal hubs must share one Stratum-1 table (2 was the identity-key reading; 6 the ping-pong)",
        )

    def test_content_equal_hubs_do_not_pingpong(self) -> None:
        """Whichever way O-3 is ruled, the ping-pong must not ship.

        ``2`` = the intern stays identity-keyed; ``1`` = it shares by
        content with a closure-CLASS validation. ``6`` is the failure this
        row exists to make impossible.
        """
        mats, mesh, quad = _slab_two_region()
        a, b = SNMesh(mesh, quad, mats), SNMesh(mesh, quad, mats)
        _require(a == b, "activation: the two hubs must be content-equal after the carve")
        builds = self._count_builds(a, b)
        _require(
            builds <= 2,
            f"content-equal hubs rebuilt the Stratum-1 table {builds} times in 6 calls",
        )
