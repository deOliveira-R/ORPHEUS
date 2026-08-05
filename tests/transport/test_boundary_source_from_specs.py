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

⭐ **The load-bearing claim is the last test in this module:**
``from_specs`` evaluates each spec at ``(|Γ₋|,) + trailing``, which is *exactly*
the shape :meth:`~orpheus.sn.boundary.angular.IncomingSourceOperator.apply` asks
for. So the two routes to :math:`q` agree **by construction** rather than by a
transcription that could drift. That property is what the later phases rely on
when they retire the inline path — without it, "the bridge computes the same
thing" would be a claim needing its own perpetual regression gate.

⚠ **What these gates deliberately do NOT claim.** They pin the *materialisation*,
not the *wiring*: nothing here asserts that a declared
:class:`~orpheus.geometry.boundary.PrescribedInflow` reaches a solve. It does not
(P2′ — the mesh-BC bridge), and a gate implying otherwise would be the
"honest metadata that gates nothing" shape this campaign keeps finding.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.geometry.boundary import ConstantInflowSource, NoSource
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.boundary.angular import IncomingSourceOperator
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


class TestTheRecipeLandsOnTheInflowSlotsOnly:
    r""":math:`q \in \Gamma_-` — by construction, not by an erasure."""

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
            def evaluate(self, shape):
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
    r"""⭐ The keystone: the bridge and the inline path evaluate the SAME shape.

    :meth:`IncomingSourceOperator.apply` asks the spec for
    ``(|Γ₋|,) + psi_out.shape[1:]``; ``from_specs`` asks for
    ``(|Γ₋|,) + slot.shape[1:]``. On a face those trailing axes are the same
    axes, so the two calls are the same call — and the equality below is a
    consequence rather than a coincidence that needs watching.
    """

    @pytest.mark.parametrize("value", [2.5, -1.0, 0.0])
    def test_the_bridge_reproduces_the_inline_operators_values(
        self, value: float
    ) -> None:
        """Bit-identical, on the inflow rows the inline path produces."""
        sn = _slab()
        trace = sn.angular_trace
        inflow = trace.inflow_indices_for_face("xmin")
        outflow = trace.outflow_indices_for_face("xmin")
        spec = ConstantInflowSource(value=value)

        bridged = np.asarray(
            AngularBoundarySourceSink.from_specs(sn, {"xmin": spec})
            .face_view("xmin")
        )[inflow]

        inline = IncomingSourceOperator(
            spec, n_inflow=int(np.size(inflow))
        ).apply(
            # apply IGNORES its input except for the trailing axes, which is
            # exactly why a zero probe of the right shape is sufficient.
            np.zeros((int(np.size(outflow)), sn.ng))
        )
        np.testing.assert_array_equal(bridged, inline)

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
        outflow = trace.outflow_indices_for_face("xmin")

        class _Ramp:
            """``q[n, g] = n + 10·g`` — distinct in BOTH axes, deliberately."""

            def evaluate(self, shape):
                n, ng = shape[0], shape[1]
                return (
                    np.arange(n, dtype=float)[:, None]
                    + 10.0 * np.arange(ng, dtype=float)[None, :]
                ).reshape(shape)

        bridged = np.asarray(
            AngularBoundarySourceSink.from_specs(sn, {"xmin": _Ramp()})
            .face_view("xmin")
        )[inflow]
        inline = IncomingSourceOperator(
            _Ramp(), n_inflow=int(np.size(inflow))
        ).apply(np.zeros((int(np.size(outflow)), sn.ng)))

        np.testing.assert_array_equal(bridged, inline)
        # Activation guard: the fixture must actually vary, or the ordering
        # claim above is vacuous.
        assert bridged.min() != bridged.max(), (
            "the ramp collapsed to a constant — this row's whole purpose is a "
            "spec that DISCRIMINATES row order, so a constant makes it blind"
        )
