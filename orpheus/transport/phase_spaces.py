r"""The transport phase spaces, built from their ingredients: a material mesh and the method's generators.

Every space a transport field lives on is a product or a sum of factors that are already typed below the method
layer: the quadrature's ordinate axis (numerics), the material mesh's scalar bulk :math:`V_E \otimes V_r` (this
package), the scheme's within-cell moment axis (:mod:`orpheus.transport.spatial`) and the traced faces of the mesh
(:mod:`orpheus.mesh`). This module spells each space once, as a raw function of those ingredients:

* :func:`angular_bulk_space` — :math:`V_\Omega \otimes V_E \otimes V_r`;
* :func:`angular_face_layout` and :func:`angular_trace_space` — the angular boundary trace, one slot per traced face
  shaped ``(N, ng, *face_shape, *moment_tail)``;
* :func:`scalar_trace_space` — the scalar ``(J⁺, J⁻)`` trace, the angular trace with the angular measure integrated out.

The angular TRIAL space is the scheme's verb,
:meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.trial_space`, applied to the angular bulk; the
full-field space is :meth:`FullFieldSpace.from_blocks
<orpheus.numerics.spaces.full_field_space.FullFieldSpace.from_blocks>` of an interior and a trace; the harmonic moment
space is :meth:`HarmonicFrame.moment_space_on
<orpheus.transport.frames.harmonic_frame.HarmonicFrame.moment_space_on>`. None of them needs a function here.

A method's problem owns these spaces as cached properties (:class:`~orpheus.sn.problem.SNProblem`,
:class:`~orpheus.diffusion.augmented_mesh.DiffusionMesh`), and each property is one call to the function below, so the
rule has one spelling and the hub keeps one instance per problem (layer DAG ruling 24). A consumer that holds a material
mesh and a quadrature, and no problem, builds the same space here: equal content, a fresh instance.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

from orpheus.mesh import trace_face_shapes
from orpheus.numerics.axis import Axis
from orpheus.numerics.face_layout import FaceLayout
from orpheus.numerics.space import FunctionSpace, factor_tail
from orpheus.numerics.spaces.angular_trace_space import AngularTraceSpace
from orpheus.numerics.spaces.scalar_trace_space import ScalarTraceSpace

if TYPE_CHECKING:
    from orpheus.numerics.quadrature.directional import Quadrature
    from orpheus.transport.mesh.material_mesh import MaterialMesh
    from orpheus.transport.spatial.scheme import DiscretizationSchemeBase

__all__ = ["angular_bulk_space", "angular_face_layout", "angular_trace_space", "scalar_trace_space"]


def angular_bulk_space(mesh: MaterialMesh, quadrature: Quadrature) -> FunctionSpace:
    r"""The per-ordinate bulk space :math:`V_\Omega \otimes V_E \otimes V_r`, axes ``(angular, energy, spatial)``.

    The ordinate axis is :meth:`Axis.of_quadrature <orpheus.numerics.axis.Axis.of_quadrature>` (the weights
    :math:`w_n`, NODAL, carrying the quadrature as its origin); the energy and spatial axes are the material mesh's
    :attr:`~orpheus.transport.mesh.material_mesh.MaterialMesh.bulk_space` axes, reused, so the scalar bulk is this
    space minus axis 0 and the two cannot disagree on the shared factors.

    No pairing of chart and quadrature is checked here. Today the only refusal of a quadrature that does not serve its
    chart is in :func:`~orpheus.sn.mesh.reduced_operator.cylindrical_streaming`, which refuses a quadrature without
    level structure. Whether suitability is a property of this space or of the streaming operator is not yet derived
    (#686).
    """
    scalar = mesh.bulk_space
    if scalar.axes is None:
        raise ValueError(f"angular_bulk_space: the material mesh's bulk space is not axis-built: {scalar!r}.")
    return FunctionSpace.of_axes(Axis.of_quadrature(quadrature), *scalar.axes)


def angular_face_layout(
    mesh: MaterialMesh, quadrature: Quadrature, scheme: DiscretizationSchemeBase,
) -> FaceLayout[str]:
    r"""The angular boundary trace's face layout: one slot per traced face, shaped ``(N, ng, *face_shape, *tail)``.

    The faces and their spatial shapes are :func:`~orpheus.mesh.axis.trace_face_shapes`, in the canonical
    :func:`~orpheus.mesh.axis.face_labels` order. The tail is the scheme's per-face transverse moment count
    :math:`p^{d-1}` through :func:`~orpheus.numerics.space.factor_tail` ("append iff > 1"), so a slopeless closure and
    every 1-D face (a point) carry no tail (#251).
    """
    n_ordinates = quadrature.n_ordinates
    moment_tail = factor_tail(scheme.face_moment_count(mesh.ndim))
    return FaceLayout.from_named_shapes([
        (label.face_name, (n_ordinates, mesh.ng, *shape, *moment_tail))
        for label, shape in trace_face_shapes(mesh.axes)
    ])


def angular_trace_space(
    mesh: MaterialMesh, quadrature: Quadrature, scheme: DiscretizationSchemeBase,
) -> AngularTraceSpace:
    r"""The angular boundary trace over :func:`angular_face_layout`, carrying the signed :math:`\Omega\cdot\hat n_f`.

    Its metric is the partial-current surface measure :math:`|\Omega\cdot\hat n_f|\,w_n`
    (:class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`).
    """
    return AngularTraceSpace.from_quadrature_and_layout(quadrature, angular_face_layout(mesh, quadrature, scheme))


def scalar_trace_space(mesh: MaterialMesh) -> ScalarTraceSpace:
    r"""The scalar boundary trace: per traced face a ``(J⁺, J⁻)`` pair per group, under the face-area metric.

    The faces are :func:`~orpheus.mesh.axis.trace_face_shapes`; each face's surface measure is the material mesh's
    :attr:`~orpheus.transport.mesh.material_mesh.MaterialMesh.areas` at its endpoint (slab 1, cylinder
    :math:`2\pi R`, sphere :math:`4\pi R^2`).

    A mesh of two or more axes is refused with :class:`ValueError`: its face measure is a per-face-cell array (the
    product of the transverse edge widths), and :class:`ScalarTraceSpace` takes one scalar per face.

    SCOPE-BOUNDARY[guard] machinery: one ``face_measure(mesh, label)`` for every dimension, with a per-face-cell measure
    on the scalar trace (#685).
    ruling: the elegance review of layer DAG Phase 4 (2026-11-01, ``.claude/plans/layer_dag.md``), filed as #685.
    revisit: when #685 lands, the 2-D scalar trace is built and this refusal retires.
    """
    if mesh.ndim != 1:
        raise ValueError(
            f"scalar_trace_space: a {mesh.ndim}-D mesh has a per-face-cell surface measure, which the scalar trace "
            "does not carry yet (#685)."
        )
    faces = trace_face_shapes(mesh.axes)
    face_areas = {
        label.face_name: float(mesh.areas[0] if label.endpoint == "min" else mesh.areas[-1])
        for label, _ in faces
    }
    return ScalarTraceSpace.for_faces([(label.face_name, shape) for label, shape in faces], mesh.ng, face_areas)
