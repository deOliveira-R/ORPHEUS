r"""The recipe → snapshot bridge: ``AngularBoundarySourceSink.from_specs``.

Phase **P1** of `.claude/plans/affine_boundary_source_channel.md`. The affine
boundary law is :math:`\gamma_-\psi = L\,\gamma_+\psi + q`, and there are two
routes to :math:`q`, related as *recipe → snapshot* rather than as duplicates:

* :meth:`~orpheus.transport.source_sinks.AngularBoundarySourceSink.prescribed_inflow`
  — the **known-per-face-array** route, already shipped;
* :meth:`~orpheus.transport.source_sinks.AngularBoundarySourceSink.from_specs`
  — the **lazy-recipe** route built here, materialising per-face
  :class:`~orpheus.geometry.boundary._source.InflowSourceSpec` generators onto
  the trace.

⭐ **The load-bearing claim is the last test in this module:** ``from_specs``
evaluates each spec at ``(|Γ₋|,) + trailing`` and writes the result to the
inflow rows **in order**. The oracle is the spec called DIRECTLY on a
hand-derived shape — the independence lives in deriving that shape and placing
those rows, which is the whole of what ``from_specs`` does.

**P3 note (2026-08-05).** Until P3 the oracle was
``IncomingSourceOperator.apply``, whose entire body was
``source.evaluate((n_inflow,) + trailing)``. Retiring that operator did NOT cost
this module an independent reference: the operator was a one-line adapter, never
the source of the independence, so the migration is a straight inlining of the
call it made. Both sides are still produced independently — one by the bridge,
one by asking the user's own spec object.

⚠ **What these gates deliberately do NOT claim.** They pin the *materialisation*,
not the *wiring*: nothing here asserts that a declared
:class:`~orpheus.geometry.boundary.PrescribedInflow` reaches a solve. That is
P2′/P3's claim, gated in ``tests/sn/solve/``; a gate here implying otherwise
would be the "honest metadata that gates nothing" shape this campaign keeps
finding.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.geometry.boundary import ConstantInflowSource, NoSource
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.transport.source_sinks import AngularBoundarySourceSink
from tests.sn._test_helpers import placeholder_materials

pytestmark = pytest.mark.l1


def _slab(n_ord: int = 8, ng: int = 2, nx: int = 4) -> SNMesh:
    """A two-face slab. ``ng = 2`` so the trailing-axis broadcast is exercised."""
    geom = StructuredGeometry(
        geometry="SLB",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(BC.vacuum, BC.vacuum),
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=nx),))
    quad = Quadrature.gauss_legendre(n_ordinates=n_ord)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


@pytest.mark.catches("ERR-047")
class TestTheRecipeLandsOnTheInflowSlotsOnly:
    r""":math:`q \in \Gamma_-` — by construction, not by an erasure.

    ⭐ **The ERR-047 marker MIGRATED here at P3** (2026-08-05), from
    ``tests/geometry/test_bc_universal_invariants.py::
    TestSourceLivesOnIncomingTraceInvariant``, whose operator-side rows asked
    the REALIZED law to deliver ``q`` and then checked which rows it landed on.
    Post-P3 the realized law delivers nothing — it is the zero morphism — so
    the postcondition had to travel with the delivery, and this is the channel
    that now delivers. (The law-tier rows there keep their own copy of the
    marker: they exercise ``assert_source_is_placeable`` itself,
    which is what the catalog entry names.)

    The catalogued hazard is unchanged in substance: an inflow source written
    into slots the sweep discards leaves the total inflow SHORT of intent. What
    changed is which code could commit it.
    """

    def test_the_inflow_rows_carry_the_recipes_value(self) -> None:
        """A constant recipe fills every inflow ordinate of the named face."""
        sn = _slab()
        q = AngularBoundarySourceSink.from_specs(
            sn, {"xmin": ConstantInflowSource(value=2.5)}
        )
        inflow = sn.angular_trace.inflow_indices_for_face("xmin")
        view = np.asarray(q.face_view("xmin"))
        np.testing.assert_array_equal(view[inflow], 2.5)

    def test_every_other_row_of_that_face_is_zero(self) -> None:
        r"""⭐ The outflow AND tangential rows stay zero.

        This is the claim that :math:`q` lives on :math:`\Gamma_-` — the
        ERR-047 mask having dissolved into the type. It is a separate leg from
        the row above because a spec that filled the WHOLE slot would satisfy
        that one and fail this one.
        """
        sn = _slab()
        q = AngularBoundarySourceSink.from_specs(
            sn, {"xmin": ConstantInflowSource(value=2.5)}
        )
        inflow = sn.angular_trace.inflow_indices_for_face("xmin")
        view = np.asarray(q.face_view("xmin"))
        others = np.ones(view.shape[0], dtype=bool)
        others[inflow] = False
        assert others.any(), (
            "this fixture has no non-inflow ordinates, so the claim is vacuous "
            "— pick a quadrature/face where Γ₊ or the tangential band is "
            "non-empty"
        )
        np.testing.assert_array_equal(view[others], 0.0)

    def test_faces_absent_from_the_mapping_stay_vacuum(self) -> None:
        """An unnamed face is untouched — the mapping is not broadcast."""
        sn = _slab()
        q = AngularBoundarySourceSink.from_specs(
            sn, {"xmin": ConstantInflowSource(value=2.5)}
        )
        np.testing.assert_array_equal(np.asarray(q.face_view("xmax")), 0.0)

    @pytest.mark.parametrize(
        "specs", [{}, {"xmin": NoSource()}], ids=["empty", "NoSource"]
    )
    def test_a_sourceless_boundary_equals_zeros_on(self, specs) -> None:
        """Both spellings of "no inflow" reproduce the canonical zero field."""
        sn = _slab()
        built = AngularBoundarySourceSink.from_specs(sn, specs)
        assert built.linf == AngularBoundarySourceSink.zeros_on(sn).linf == 0.0


class TestTheContractIsEnforced:
    """The one invariant ``InflowSourceSpec`` promises, and the face names."""

    def test_a_spec_returning_the_wrong_shape_is_refused(self) -> None:
        r"""The Protocol's ONE invariant: the array has the requested shape.

        Caught at materialisation because it cannot be caught later — an
        over-long array would raise deep in the packing with a message about
        slots, and an array of the wrong TRAILING shape could broadcast
        silently into the wrong group.
        """
        sn = _slab()

        class _Liar:
            def evaluate(self, space):
                shape = tuple(space.shape)
                return np.zeros((shape[0] + 1,) + shape[1:])

        with pytest.raises(ValueError, match="expected"):
            AngularBoundarySourceSink.from_specs(sn, {"xmin": _Liar()})

    def test_an_unknown_face_is_refused(self) -> None:
        """A typo'd face name is a silent no-op unless it raises."""
        sn = _slab()
        with pytest.raises(ValueError, match="not a face of the layout"):
            AngularBoundarySourceSink.from_specs(
                sn, {"xmid": ConstantInflowSource(value=1.0)}
            )

    def test_the_positive_leg_of_both_refusals(self) -> None:
        """An honest spec on a real face constructs — so the rows above pin
        the DISAGREEMENT, not "from_specs rejects everything"."""
        sn = _slab()
        q = AngularBoundarySourceSink.from_specs(
            sn, {"xmin": ConstantInflowSource(value=1.0)}
        )
        assert isinstance(q, AngularBoundarySourceSink)


class TestTheTwoRoutesAgreeByConstruction:
    r"""⭐ The keystone: the bridge delivers exactly what the SPEC produces.

    ``from_specs`` asks each spec for ``(|Γ₋|,) + slot.shape[1:]`` and writes
    the answer to that face's inflow rows. The oracle asks the same spec for
    ``(|Γ₋|, ng)`` derived by hand from the trace. Independence is intact: one
    side is the bridge, the other is the user's own spec object; neither calls
    the other.

    Until **P3** the oracle read ``IncomingSourceOperator(spec,
    n_inflow=…).apply(probe)``, whose body was precisely
    ``spec.evaluate((n_inflow,) + probe.shape[1:])`` — so retiring the operator
    inlines the call rather than removing a reference. See the module docstring.
    """

    @pytest.mark.parametrize("value", [2.5, -1.0, 0.0])
    def test_the_bridge_reproduces_the_specs_own_values(
        self, value: float
    ) -> None:
        """Bit-identical to the spec's own output, on the inflow rows."""
        sn = _slab()
        trace = sn.angular_trace
        inflow = trace.inflow_indices_for_face("xmin")
        spec = ConstantInflowSource(value=value)

        bridged = np.asarray(
            AngularBoundarySourceSink.from_specs(sn, {"xmin": spec})
            .face_view("xmin")
        )[inflow]

        direct = spec.evaluate(trace.inflow_space("xmin"))
        np.testing.assert_array_equal(bridged, direct)

    def test_a_nonconstant_recipe_also_agrees(self) -> None:
        r"""⭐ The negative leg — a constant recipe cannot see an ORDERING bug.

        Every assertion above uses a spec whose output is uniform, so a bridge
        that shuffled the inflow rows, or wrote them transposed, would pass
        all of them. This recipe's value varies along the ordinate axis, so
        the comparison pins the row ORDER as well as the values.
        """
        sn = _slab()
        trace = sn.angular_trace
        inflow = trace.inflow_indices_for_face("xmin")

        class _Ramp:
            """``q[n, g] = n + 10·g`` — distinct in BOTH axes, deliberately."""

            def evaluate(self, space):
                shape = tuple(space.shape)
                n, ng = shape[0], shape[1]
                return (
                    np.arange(n, dtype=float)[:, None]
                    + 10.0 * np.arange(ng, dtype=float)[None, :]
                ).reshape(shape)

        bridged = np.asarray(
            AngularBoundarySourceSink.from_specs(sn, {"xmin": _Ramp()})
            .face_view("xmin")
        )[inflow]
        direct = _Ramp().evaluate(trace.inflow_space("xmin"))

        np.testing.assert_array_equal(bridged, direct)
        # Activation guard: the fixture must actually vary, or the ordering
        # claim above is vacuous.
        assert bridged.min() != bridged.max(), (
            "the ramp collapsed to a constant — this row's whole purpose is a "
            "spec that DISCRIMINATES row order, so a constant makes it blind"
        )
