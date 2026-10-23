r"""The declaration -> snapshot bridge: ``AngularBoundarySourceSink.from_binding`` (file name kept: it was ``from_specs``).

Phase **P1** of `.claude/plans/archive/affine_boundary_source_channel.md`, re-posed at wave 2 S2b. The affine boundary law is
:math:`\gamma_-\psi = L\,\gamma_+\psi + q`, and there are two routes to :math:`q`, related as *declaration -> snapshot* rather
than as duplicates:

* :meth:`~orpheus.transport.fields.AngularBoundarySourceSink.prescribed_inflow` — the **known-per-face-array** route;
* :meth:`~orpheus.transport.fields.AngularBoundarySourceSink.from_binding` — the **declared-condition** route: it reads the
  sources of a ``Binding`` (the tagged faces whose bound ``WallCondition`` has a nonzero source), samples each on that face's inflow
  half-trace through the method's realizer table, and packs the result.

The earlier bridge ``from_specs`` took a per-face mapping of ``InflowSourceSpec`` recipes; the recipe is now the condition's
``InflowTerms`` (``WallCondition.uniform(v)`` for a constant, ``WallCondition.framed(FaceFrameInflow(sampler))`` for a profile),
the mapping is the binding (so an unknown face is unspellable), and the sample goes through ``SN_REALIZERS.sample_inflow``.

⭐ **The load-bearing claim is the last class in this module:** the bridge evaluates each recipe at the face's inflow half-trace and
writes the result to the inflow rows **in order**. The oracle is the sampler called DIRECTLY on the half-trace space.

⚠ **What these gates deliberately do NOT claim.** They pin the *materialisation*, not the *wiring*: that a declared inflow reaches a
solve is gated in ``tests/gates/sn/solve/test_declared_inflow_reaches_the_rhs.py``.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import StructuredGeometry
from orpheus.mesh import CellsByCount, Mesher
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.boundary.errors import BoundarySourceNotOnIncomingTraceError
from orpheus.sn.boundary.realizer import SN_REALIZERS
from orpheus.sn.problem import SNProblem
from orpheus.transport.boundary_conditions import FaceFrameInflow, FaceSite, WallCondition
from orpheus.transport.fields import AngularBoundarySourceSink
from tests.gates.sn._test_helpers import placeholder_materials

pytestmark = pytest.mark.l1


def _slab(xmin: WallCondition | None = None, xmax: WallCondition | None = None, n_ord: int = 8, ng: int = 2, nx: int = 4) -> SNProblem:
    """A two-face slab, tag 0 on ``xmin`` and tag 1 on ``xmax`` (vacuum unless a condition is given). ``ng = 2`` so the trailing-axis
    broadcast is exercised."""
    mesh = Mesher(StructuredGeometry.slab((0.0, 2.0), (0,), left=0, right=1)).partition(CellsByCount.uniform_volume(nx)).mesh
    return SNProblem(
        mesh, Quadrature.gauss_legendre(n_ordinates=n_ord), placeholder_materials(ng=ng),
        boundary_conditions={0: xmin or WallCondition.vacuum(), 1: xmax or WallCondition.vacuum()},
    )


@pytest.mark.catches("ERR-047")
class TestTheRecipeLandsOnTheInflowSlotsOnly:
    r""":math:`q \in \Gamma_-` — by construction, not by an erasure.

    ⭐ **The ERR-047 marker MIGRATED here at P3** (2026-08-05) from the law-tier rows of
    ``tests/gates/geometry/test_bc_universal_invariants.py``, and at wave 2 S2b the refusal that the SN side keeps is
    ``SNSource.sample``'s: a source sampled on a space that is NOT the inflow half-trace of its face is refused
    (``BoundarySourceNotOnIncomingTraceError``). The catalogued hazard is unchanged in substance: an inflow source written into
    slots the sweep discards leaves the total inflow SHORT of intent. What changed is which code could commit it.
    """

    def test_the_inflow_rows_carry_the_recipes_value(self) -> None:
        """A constant recipe fills every inflow ordinate of the named face."""
        sn = _slab(xmin=WallCondition.uniform(2.5))
        q = sn.declared_boundary_source()
        inflow = sn.angular_trace.inflow_indices_for_face("xmin")
        np.testing.assert_array_equal(np.asarray(q.face_view("xmin"))[inflow], 2.5)

    def test_every_other_row_of_that_face_is_zero(self) -> None:
        r"""⭐ The outflow AND tangential rows stay zero.

        This is the claim that :math:`q` lives on :math:`\Gamma_-` — the ERR-047 mask having dissolved into the type. It is a
        separate leg from the row above because a recipe that filled the WHOLE slot would satisfy that one and fail this one.
        """
        sn = _slab(xmin=WallCondition.uniform(2.5))
        inflow = sn.angular_trace.inflow_indices_for_face("xmin")
        view = np.asarray(sn.declared_boundary_source().face_view("xmin"))
        others = np.ones(view.shape[0], dtype=bool)
        others[inflow] = False
        assert others.any(), (
            "this fixture has no non-inflow ordinates, so the claim is vacuous — pick a quadrature/face where Γ₊ or the "
            "tangential band is non-empty"
        )
        np.testing.assert_array_equal(view[others], 0.0)

    def test_faces_absent_from_the_declarations_stay_vacuum(self) -> None:
        """An undeclared source is untouched — the declaration is not broadcast."""
        sn = _slab(xmin=WallCondition.uniform(2.5))
        np.testing.assert_array_equal(np.asarray(sn.declared_boundary_source().face_view("xmax")), 0.0)

    @pytest.mark.parametrize("wall", [WallCondition.vacuum(), WallCondition.uniform(0.0)], ids=["vacuum", "zero-source"])
    def test_a_sourceless_boundary_equals_the_zero_trace(self, wall) -> None:
        """Both spellings of "no inflow" reproduce the canonical zero field."""
        sn = _slab(xmin=wall)
        assert sn.declared_boundary_source().linf == AngularBoundarySourceSink.zeros(sn.angular_trace).linf == 0.0

    def test_a_source_sampled_off_the_inflow_tier_is_refused(self) -> None:
        r"""⭐ ERR-047 at the realizer: the sampler is handed a space that is not :math:`\Gamma_-` of the face and the source realizer
        refuses it, naming the tier. Run on the face's OUTFLOW half-trace (a real space of the same size, so no shape check
        could tell) with the inflow half-trace as the positive control."""
        sn = _slab()
        site = FaceSite("xmin", sn.mesh.coord, 1)
        terms = WallCondition.uniform(1.0).source
        inflow, outflow = sn.angular_trace.inflow_space("xmin"), sn.angular_trace.outflow_space("xmin")
        np.testing.assert_array_equal(SN_REALIZERS.sample_inflow(terms, site, inflow), 1.0)
        with pytest.raises(BoundarySourceNotOnIncomingTraceError):
            SN_REALIZERS.sample_inflow(terms, site, outflow)


class TestTheContractIsEnforced:
    """The one invariant a sampler promises: the array has the requested shape."""

    def test_a_sampler_returning_the_wrong_shape_is_refused(self) -> None:
        r"""Caught at materialisation because it cannot be caught later — an over-long array would raise deep in the packing with a
        message about slots, and an array of the wrong TRAILING shape could broadcast silently into the wrong group."""

        class _Liar:
            def sample(self, space):
                shape = tuple(space.shape)
                return np.zeros((shape[0] + 1,) + shape[1:])

        sn = _slab(xmin=WallCondition.framed(FaceFrameInflow(_Liar())))
        with pytest.raises(ValueError, match="not exactly space.shape"):
            sn.declared_boundary_source()

    def test_the_positive_leg_of_the_refusal(self) -> None:
        """An honest sampler on a real face constructs — so the row above pins the DISAGREEMENT, not "everything is rejected"."""

        class _Honest:
            def sample(self, space):
                return np.ones(space.shape)

        q = _slab(xmin=WallCondition.framed(FaceFrameInflow(_Honest()))).declared_boundary_source()
        assert isinstance(q, AngularBoundarySourceSink)


class TestTheTwoRoutesAgreeByConstruction:
    r"""⭐ The keystone: the bridge delivers exactly what the SAMPLER produces.

    ``from_binding`` asks each recipe for ``space.shape`` on the face's inflow half-trace and writes the answer to that face's inflow
    rows. The oracle asks the same sampler for the same space, derived by hand from the trace. Independence is intact: one side is the
    bridge, the other is the user's own sampler object; neither calls the other.
    """

    # (the value -1.0 of the earlier recipe row is gone: a negative inflow is outside the positive cone and ``WallCondition.uniform``
    # refuses it at construction, ``BoundaryResponseNotPositiveError``, gated in ``test_wall_condition.py``; 7.25 takes its place)
    @pytest.mark.parametrize("value", [2.5, 7.25, 0.0])
    def test_the_bridge_reproduces_the_recipes_own_values(self, value: float) -> None:
        """Bit-identical to the recipe's own output, on the inflow rows."""
        sn = _slab(xmin=WallCondition.uniform(value))
        trace = sn.angular_trace
        inflow = trace.inflow_indices_for_face("xmin")
        bridged = np.asarray(sn.declared_boundary_source().face_view("xmin"))[inflow]
        direct = value * np.ones(trace.inflow_space("xmin").shape)
        np.testing.assert_array_equal(bridged, direct)

    def test_a_nonconstant_recipe_also_agrees(self) -> None:
        r"""⭐ The negative leg — a constant recipe cannot see an ORDERING bug.

        Every assertion above uses a recipe whose output is uniform, so a bridge that shuffled the inflow rows, or wrote them
        transposed, would pass all of them. This sampler's value varies along the ordinate axis AND the group axis, so the comparison
        pins the row ORDER as well as the values.
        """

        class _Ramp:
            """``q[n, g] = n + 10·g`` — distinct in BOTH axes, deliberately."""

            def sample(self, space):
                shape = tuple(space.shape)
                n, ng = shape[0], shape[1]
                return (np.arange(n, dtype=float)[:, None] + 10.0 * np.arange(ng, dtype=float)[None, :]).reshape(shape)

        sn = _slab(xmin=WallCondition.framed(FaceFrameInflow(_Ramp())))
        trace = sn.angular_trace
        inflow = trace.inflow_indices_for_face("xmin")
        bridged = np.asarray(sn.declared_boundary_source().face_view("xmin"))[inflow]
        direct = _Ramp().sample(trace.inflow_space("xmin"))

        np.testing.assert_array_equal(bridged, direct)
        # Activation guard: the fixture must actually vary, or the ordering claim above is vacuous.
        assert bridged.min() != bridged.max(), (
            "the ramp collapsed to a constant — this row's whole purpose is a recipe that DISCRIMINATES row order, so a constant "
            "makes it blind"
        )


# ``TestTheTraceDoorsRefuseANonSNCarrier`` (Layer-DAG P3) is RETIRED with its subject (wave 2 S2b): it pinned the typed doors
# ``from_problem_laws`` / ``prescribed_inflow`` opened on a ``MaterialMesh`` carrier (``TransportMethod.bc`` and
# ``mesh.angular_trace``), refusing a bare material mesh and a diffusion mesh by name ("owns no angular trace"). ``from_binding``
# takes its parts explicitly (the binding, the trace, the realizer table, the face-site hook) and no carrier, so the doors have no
# subject: a caller that has no trace has nothing to pass for ``trace``. The diffusion half of the control (``DiffusionMesh``
# owns no angular trace) is slice S3's. ``prescribed_inflow`` keeps its own door, gated by its callers.
