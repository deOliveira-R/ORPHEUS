r"""B2 regression floor — production ASKS the factors, and gets the old answer.

Campaign phase B1 populated
:class:`~orpheus.geometry.boundary.BoundaryGeometryMap` /
:class:`~orpheus.geometry.boundary.BoundaryResponseKernel` on all seven laws
(pinned in :mod:`tests.geometry.test_boundary_factors`). Phase **B2** made them
reachable — B2.0 stopped discarding the law at realization — and repointed the
five production sites that had been answering structural questions by comparing
``law.kind`` strings against literals.

This file is the floor under that repoint. It carries two claims that must be
kept apart:

**Equivalence** (:class:`TestTagEquivalence`) — for every registered law, each
repointed predicate returns exactly what the retired tag expression returned.
The retired expressions are reproduced verbatim below, so this is a genuine
old-vs-new comparison rather than a restatement of the new code. It is what
makes B2 a *repoint* rather than a redefinition, and it reddens if a law is
ever added whose factors disagree with its own tag.

**Consumption** (:class:`TestSpecMutationsPropagate`) — mutating a law's spec
object must MOVE each dependent site. Equivalence alone is Mode-11 blind: a
site that ignores the factors entirely and hard-codes the old answer would pass
every row above. Only a red mutation proves the site reads the spec.

Tagged ``@pytest.mark.foundation``: structural contracts about dispatch, not
discretization claims, so no ``verifies(...)`` label (the verifies⊥level
doctrine).
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from orpheus.diffusion.boundary_realizer import DiffusionBoundaryRealizer
from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.geometry.boundary import (
    AlbedoBoundary,
    BoundaryError,
    BoundaryTraceLaw,
    IdentityMap,
    PrescribedInflow,
    ReflectiveBoundary,
    ScalarResponse,
    VacuumInflow,
    WhiteBoundary,
    ZeroFluxBoundary,
)
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.acceleration.dsa import DSALowOrderSystem
from orpheus.sn.loss_representation.sweep_schedule import _reflective_faces
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.operators.boundary import (
    RadialCharacteristicBoundaryOperator,
    _has_ruled_corner_action,
)
from tests.sn._test_helpers import placeholder_materials


pytestmark = pytest.mark.foundation


#: The retired tag frozensets, reproduced VERBATIM from the constants B2
#: deleted — ``_RULED_CORNER_KINDS`` (``sn/operators/boundary.py``) and
#: ``_SUPPORTED_BC`` (``sn/acceleration/dsa.py``). Both were
#: ``frozenset({"vacuum", "reflective"})``. Keeping the literal here is what
#: makes the equivalence rows below an old-vs-new comparison; importing a
#: shared constant would make them tautologies.
_RETIRED_RULED_CORNER_KINDS = frozenset({"vacuum", "reflective"})
_RETIRED_SUPPORTED_BC = frozenset({"vacuum", "reflective"})


def _every_registered_law() -> list[BoundaryTraceLaw]:
    """One instance of each registered law, with a NON-DEFAULT albedo.

    ``albedo=0.42`` matters: at the default 1.0 a response-vs-geometry mix-up
    would be invisible, since ``ScalarResponse(1.0)`` and "is a mirror" both
    read as truthy. Off-default, the two are distinguishable.
    """
    laws = []
    for _key, cls in sorted(BoundaryTraceLaw.registry.items()):
        fields = {f.name for f in dataclasses.fields(cls)}
        laws.append(cls(albedo=0.42) if "albedo" in fields else cls())
    return laws


def _retired_diffusion_albedo(law: BoundaryTraceLaw) -> object:
    """``DiffusionBoundaryRealizer._partial_current_albedo``, as it was.

    The five-arm ``isinstance`` ladder B2 collapsed into
    ``law.response_kernel.amplitude``, transcribed from the pre-B2 source.
    """
    if isinstance(law, VacuumInflow):
        return 0.0
    if isinstance(law, ZeroFluxBoundary):
        return -1.0
    if isinstance(law, (ReflectiveBoundary, WhiteBoundary)):
        return float(law.albedo)
    if isinstance(law, AlbedoBoundary):
        return float(law.albedo)
    return "REFUSE"


def _live_diffusion_albedo(law: BoundaryTraceLaw) -> object:
    try:
        return DiffusionBoundaryRealizer._partial_current_albedo(law)
    except BoundaryError:
        return "REFUSE"


#: ``(site, retired tag expression, live expression, calls_production)``.
#:
#: **``calls_production`` is load-bearing, not decoration.** Two of the four
#: live sides invoke the production reader directly, so a bug inside it reddens
#: the row. The other two — the sweep schedule and the DSA admission — need a
#: whole ``SNMesh`` to reach production, which cannot be built for the five
#: laws SN does not admit; their live side is the predicate EXPRESSION, so the
#: row proves the expression is equivalent but NOT that production evaluates
#: it. What ties those two to production is
#: :class:`TestSpecMutationsPropagate` (measured: breaking
#: ``_reflective_faces`` reddens the sweep mutation test and leaves this table
#: green). Stating the split is the point — an unlabelled restatement row is a
#: Mode-11 gate wearing an equivalence gate's name.
#:
#: The one deliberate divergence — the solver's leakage list — is NOT here; it
#: has its own test below, which states the divergence rather than hiding it.
_EQUIVALENCE_SITES = [
    (
        "sweep_schedule._reflective_faces",
        lambda law: type(law).key == "reflective",
        lambda law: law.geometry_map.permutes_ordinates,
        False,
    ),
    (
        "RadialCharacteristicBoundaryOperator ruled-corner admission",
        lambda law: type(law).key in _RETIRED_RULED_CORNER_KINDS,
        _has_ruled_corner_action,
        True,
    ),
    (
        "DSALowOrderSystem low-order-row admission",
        lambda law: type(law).key in _RETIRED_SUPPORTED_BC,
        lambda law: not isinstance(law, PrescribedInflow)
        and (
            law.response_kernel.is_zero
            or law.geometry_map.permutes_ordinates
        ),
        False,
    ),
    (
        "DiffusionBoundaryRealizer._partial_current_albedo",
        _retired_diffusion_albedo,
        _live_diffusion_albedo,
        True,
    ),
]


class TestTagEquivalence:
    """Old tag expression == new structural expression, law by law."""

    @pytest.mark.parametrize(
        "site,retired,live,calls_production",
        _EQUIVALENCE_SITES,
        ids=[s.split(".")[0].split(" ")[0] for s, _, _, _ in _EQUIVALENCE_SITES],
    )
    def test_structural_predicate_matches_the_retired_tag_expression(
        self, site, retired, live, calls_production,
    ) -> None:
        divergences = [
            (type(law).__name__, retired(law), live(law))
            for law in _every_registered_law()
            if retired(law) != live(law)
        ]
        if divergences:
            rows = "\n".join(
                f"    {name}: retired={old!r}  live={new!r}"
                for name, old, new in divergences
            )
            reach = (
                "this row CALLS the production reader, so the bug may be in it"
                if calls_production else
                "this row checks the predicate EXPRESSION only (production "
                "needs an SNMesh) — if the expression is right, look at "
                "TestSpecMutationsPropagate for the production wiring"
            )
            pytest.fail(
                f"{site} answers differently than the tag expression it "
                f"replaced:\n{rows}\n"
                f"({reach}.)\n"
                f"B2 is a REPOINT — the structural form must return what the "
                f"string comparison returned. A divergence here is either a "
                f"repoint bug or a deliberate semantic change that belongs in "
                f"its own test, stated, the way the solver-leakage one is."
            )

    def test_the_one_deliberate_divergence_is_stated_not_hidden(self) -> None:
        r"""``solver.compute_keff``'s leakage list DOES diverge — on one law.

        The retired test was ``op.kind == "vacuum"``; the live one is
        ``op.law.response_kernel.is_zero``. They agree on six of seven laws and
        differ on :class:`PrescribedInflow`, which returns nothing to the domain
        (:math:`R = 0`) and therefore leaks its whole outflow — so the live
        predicate includes it and the tag test missed it.

        The new answer is the correct one, and it is unreachable through the
        production path (SN's admission table is ``{reflective, vacuum}``, so no
        BC tag installs a prescribed law on an SN face). This test exists so the
        divergence stays a recorded decision instead of drifting into folklore.
        """
        retired, live = [], []
        for law in _every_registered_law():
            if type(law).key == "vacuum":
                retired.append(type(law).__name__)
            if law.response_kernel.is_zero:
                live.append(type(law).__name__)
        assert retired == ["VacuumInflow"]
        assert live == ["PrescribedInflow", "VacuumInflow"]
        # And the reachability claim the decision rests on.
        assert set(SNMesh.BOUNDARY_OPERATOR_REGISTRY) == {"reflective", "vacuum"}


def _slab(left: str, right: str) -> SNMesh:
    return SNMesh(
        Mesh1D(
            edges=np.linspace(0.0, 1.0, 5),
            mat_ids=np.zeros(4, dtype=int),
            coord=CoordSystem.CARTESIAN,
            bc_left=BC(left),
            bc_right=BC(right),
        ),
        Quadrature.gauss_legendre(4),
        placeholder_materials(),
    )


def _sphere(outer: str) -> SNMesh:
    return SNMesh(
        Mesh1D(
            edges=np.linspace(0.0, 1.0, 5),
            mat_ids=np.zeros(4, dtype=int),
            coord=CoordSystem.SPHERICAL,
            bc_left=BC("reflective"),
            bc_right=BC(outer),
        ),
        Quadrature.gauss_legendre(8),
        placeholder_materials(),
    )


class TestSpecMutationsPropagate:
    """Mutate a law's spec — every dependent site must MOVE.

    Equivalence is necessary but Mode-11 blind on its own: a site that never
    reads the factors and hard-codes the old answer passes every row above.
    These are the reddening mutations that prove consumption.
    """

    def test_sweep_schedule_moves_when_the_mirror_stops_permuting(
        self, monkeypatch,
    ) -> None:
        sn = _slab("reflective", "vacuum")
        assert _reflective_faces(sn) == frozenset({"xmin"})
        monkeypatch.setattr(
            ReflectiveBoundary, "geometry_map",
            property(lambda self: IdentityMap()),
        )
        mutated = _reflective_faces(_slab("reflective", "vacuum"))
        if mutated:
            pytest.fail(
                f"the octant Gauss-Seidel schedule still sees {sorted(mutated)} "
                f"as reflective after the law's geometry map stopped permuting "
                f"ordinates — it is not reading `geometry_map`."
            )

    def test_corner_adjointability_moves_when_the_mirror_stops_permuting(
        self, monkeypatch,
    ) -> None:
        assert RadialCharacteristicBoundaryOperator(
            _sphere("reflective")).is_adjointable
        monkeypatch.setattr(
            ReflectiveBoundary, "geometry_map",
            property(lambda self: IdentityMap()),
        )
        if RadialCharacteristicBoundaryOperator(
                _sphere("reflective")).is_adjointable:
            pytest.fail(
                "B_b still advertises a transpose after its outer law stopped "
                "permuting ordinates — the corner swap it would transpose is "
                "no longer expressible, so this is a promise it cannot keep."
            )

    def test_corner_adjointability_moves_when_vacuum_stops_returning_nothing(
        self, monkeypatch,
    ) -> None:
        assert RadialCharacteristicBoundaryOperator(
            _sphere("vacuum")).is_adjointable
        monkeypatch.setattr(
            VacuumInflow, "response_kernel",
            property(lambda self: ScalarResponse(0.5)),
        )
        if RadialCharacteristicBoundaryOperator(
                _sphere("vacuum")).is_adjointable:
            pytest.fail(
                "B_b still advertises a transpose after its outer law grew a "
                "non-zero response — a half-returning corner has no ruled "
                "action, so it is not reading `response_kernel`."
            )

    @pytest.mark.parametrize(
        "target,attr,value,label",
        [
            (VacuumInflow, "response_kernel", ScalarResponse(0.5),
             "vacuum's response stops being zero"),
            (ReflectiveBoundary, "geometry_map", IdentityMap(),
             "reflective's geometry stops permuting"),
        ],
        ids=["vacuum-response", "reflective-geometry"],
    )
    def test_dsa_admission_moves_when_a_proven_row_loses_its_property(
        self, monkeypatch, target, attr, value, label,
    ) -> None:
        DSALowOrderSystem.from_sn_mesh(_slab("vacuum", "reflective"))
        monkeypatch.setattr(target, attr, property(lambda self: value))
        with pytest.raises(NotImplementedError, match="no proven low-order row"):
            DSALowOrderSystem.from_sn_mesh(_slab("vacuum", "reflective"))

    def test_dsa_boundary_ROW_selection_moves_too_not_just_admission(
        self, monkeypatch,
    ) -> None:
        r"""The row selector reads the spec, not only the admission gate.

        Distinct from the admission mutations above: this one keeps the law
        ADMITTED (a mirror still permutes) while flipping the property the row
        choice turns on, so only ``_build``'s
        ``if law.response_kernel.is_zero`` can react. Reflective at
        :math:`R = 0` must take the Marshak (38) row instead of the (39) row,
        which changes the assembled edge operator.
        """
        base = DSALowOrderSystem.from_sn_mesh(
            _slab("vacuum", "reflective")).a_low.copy()
        monkeypatch.setattr(
            ReflectiveBoundary, "response_kernel",
            property(lambda self: ScalarResponse(0.0)),
        )
        mutated = DSALowOrderSystem.from_sn_mesh(
            _slab("vacuum", "reflective")).a_low
        delta = float(np.abs(base - mutated).max())
        if delta == 0.0:
            pytest.fail(
                "the low-order boundary row is unchanged after the reflective "
                "law's response went to zero — `_build` is not selecting the "
                "Marshak vs (39) row from `response_kernel`, so a future law "
                "would silently get the wrong closure."
            )

    @pytest.mark.parametrize(
        "law_cls,attr,value",
        [
            (VacuumInflow, "response_kernel", ScalarResponse(0.75)),
            (AlbedoBoundary, "response_kernel", ScalarResponse(0.11)),
        ],
        ids=["vacuum", "albedo"],
    )
    def test_diffusion_albedo_follows_the_response_factor(
        self, monkeypatch, law_cls, attr, value,
    ) -> None:
        monkeypatch.setattr(law_cls, attr, property(lambda self: value))
        got = DiffusionBoundaryRealizer._partial_current_albedo(law_cls())
        assert got == value.amplitude, (
            f"diffusion's 𝒜 for {law_cls.__name__} is {got!r}, not the law's "
            f"response {value.amplitude!r} — the ladder was collapsed onto "
            f"`response_kernel`, so 𝒜 must track it."
        )
