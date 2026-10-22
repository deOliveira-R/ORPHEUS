r"""L2 typed transport fields — every class inherits from
:class:`~orpheus.numerics.field.Field`.

The Depth-B unification (plan §3.3): every typed transport field is
``(values, space) + algebra`` (the L1 ABC) plus domain-specific
fields (``mesh``, ``boundary``, …). The dunder algebra is inherited
verbatim; concrete subclasses add construction factories and
diagnostics tailored to their physical role.

This package holds both halves of every role pair: the flux leaves and the
source/sink leaves. Importing it imports every leaf that declares a pair, so
the flux ↔ source/sink map of
:class:`~orpheus.transport.fields._bases.RolePair` is complete whenever the
package is. The pairs themselves are the registry's, not a list here:
:class:`~orpheus.transport.fields._bases.RolePair` holds them and
``tests/gates/transport/fields/test_role_partners.py`` enumerates them.

Migration sequence (Depth B plan §6, with the as-built outcome):

* D-D — :class:`ScalarFlux` (simplest case).
* D-E — :class:`HarmonicMomentFlux` (cleanest gap).
* D-F — :class:`AngularSourceSink` / :class:`AngularBoundarySourceSink`
  (the source/sink role leaves; they live in this package, below).
* D-G — :class:`AngularBoundaryFlux` (the boundary trace cochain
  :math:`C^1_\partial`, on a
  :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`).
* D-H — :class:`AngularFlux` (the phase-space directional flux, on the
  ``TimedFullField`` composite).

Why ``SourceSink`` (not ``Source``)
===================================

A **source** (production: external :math:`q_{\rm ext}`, in-scattering
:math:`S\psi`, fission :math:`F\psi`) and a **sink** (loss: an operator
output such as the streaming+collision rate :math:`L\psi`) are *the same
physical quantity with opposite sign* — a signed rate-density term in the
balance :math:`L\psi + C\psi - S\psi - F\psi = q`. They share one
dimensional signature (:data:`~orpheus.numerics.units.ANGULAR_RATE_UNITS`
/ :data:`~orpheus.numerics.units.SCALAR_RATE_UNITS` in the bulk; flux
units on the trace) and one storage shape, and they add freely. So the
role is one leaf, and the name says so: ``SourceSink``. (The original
``Source`` name read as production-only and hid that operator ``.apply``
outputs — sinks — land in the *same* type; the rename, B.5, removes that
confusion at the signature level.)

These are Field-shaped objects with DIFFERENT physical units than
fluxes. Per Issue #207's architectural principle and Depth B plan §3.6,
cross-class arithmetic between a source/sink and a flux is forbidden by
the Field ABC's Layer 1 class-identity gate — even when shape and units
might align. Cross-class same-units operations (e.g. combining an
isotropic source/sink with a per-ordinate one) go through the canonical
subspace-containment injection or a NAMED composition factory, not an
ad-hoc cross-class dunder.

Interior face access — NO per-face type
========================================

An early sketch of D-G paired :class:`AngularBoundaryFlux` with a per-face
``BoundaryFaceFlux`` object. That object was **rejected** (Wave O
#205, cross-domain-attacker frame memo
``field_role_typing_faceflux_frames.md``): a cochain is an assignment
of values to faces, and the load-bearing vectorisation constraint (the
``(N_oct, ng, n_diag)`` wavefront batch is the unit of operation)
forbids a per-face Python object — it would reintroduce the per-cell
fold that caused a 10–20× regression. Per-face access is therefore a
**zero-copy view** (``face_view``), exactly as :class:`AngularBoundaryFlux`
does. The *interior* cell-face cochain ``C^1_int`` was, through S6.4,
the typed ``WavefrontFlux`` (ephemeral, flux-only; with
:class:`AngularBoundaryFlux` ``C^1_∂`` it biproduct-decomposed the full face
cochain ``C^1 = C^1_int ⊕ C^1_∂``). RETIRED at S6.4(f): the concept
lives on in its two native realizations — the rolling front
(``orpheus.sn.loss_representation.sweep_graph._MovingFrontier``, the values AT the moving
wavefront, per-level ι_*-seeded / ι*-shed) and the per-octant full
cochain history (``FullFieldWavefront._octant_face_cochain``); the
whole-trace boundary exchange is the shared ``_OctantWalk`` frame's
job. See the theory page ``operator_algebra.rst``
§ "The interior face-flux cochain" for the succession history.
"""

from __future__ import annotations

from orpheus.transport.fields.cross_section_field import CrossSectionField
from orpheus.transport.fields.harmonic_moment_flux import HarmonicMomentFlux
from orpheus.transport.fields.scalar_boundary_flux import ScalarBoundaryFlux
from orpheus.transport.fields.scalar_flux import ScalarFlux
from orpheus.transport.fields.radial_characteristic_interior_flux import (
    RadialCharacteristicInteriorFlux,
)
from orpheus.transport.fields.radial_characteristic_boundary_flux import (
    RadialCharacteristicBoundaryFlux,
)
from orpheus.transport.fields.scalar_source_sink import ScalarSourceSink
from orpheus.transport.fields.angular_source_sink import AngularSourceSink
from orpheus.transport.fields.angular_boundary_source_sink import AngularBoundarySourceSink
from orpheus.transport.fields.scalar_boundary_source_sink import (
    ScalarBoundarySourceSink,
)
from orpheus.transport.fields.harmonic_moment_source_sink import (
    HarmonicMomentSourceSink,
)
from orpheus.transport.fields.radial_characteristic_interior_source_sink import (
    RadialCharacteristicInteriorSourceSink,
)
from orpheus.transport.fields.radial_characteristic_boundary_source_sink import (
    RadialCharacteristicBoundarySourceSink,
)

__all__ = [
    "CrossSectionField",
    "HarmonicMomentFlux",
    "ScalarBoundaryFlux",
    "ScalarFlux",
    "RadialCharacteristicInteriorFlux",
    "RadialCharacteristicBoundaryFlux",
    "ScalarSourceSink",
    "AngularSourceSink",
    "AngularBoundarySourceSink",
    "ScalarBoundarySourceSink",
    "HarmonicMomentSourceSink",
    "RadialCharacteristicInteriorSourceSink",
    "RadialCharacteristicBoundarySourceSink",
]
