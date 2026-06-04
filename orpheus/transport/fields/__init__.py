r"""L2 typed transport fields — every class inherits from
:class:`~orpheus.numerics.field.Field`.

The Depth-B unification (plan §3.3): every typed transport field is
``(values, space) + algebra`` (the L1 ABC) plus domain-specific
fields (``mesh``, ``boundary``, …). The dunder algebra is inherited
verbatim; concrete subclasses add construction factories and
diagnostics tailored to their physical role.

Migration sequence (Depth B plan §6, with the as-built outcome):

* D-D — :class:`ScalarFlux` (simplest case).
* D-E — :class:`HarmonicMomentField` (cleanest gap).
* D-F — :class:`AngularSourceSink` / :class:`BoundarySourceSink`
  (the source/sink role leaves; see ``orpheus.transport.source_sinks``).
* D-G — :class:`BoundaryFlux` (the boundary trace cochain
  :math:`C^1_\partial`, on a
  :class:`~orpheus.numerics.spaces.trace_space.TraceSpace`).
* D-H — :class:`AngularFlux` (the phase-space directional flux, on the
  ``TimedFullField`` composite).

Interior face access — NO per-face type
========================================

An early sketch of D-G paired :class:`BoundaryFlux` with a per-face
``BoundaryFaceFlux`` object. That object was **rejected** (Wave O
#205, cross-domain-attacker frame memo
``field_role_typing_faceflux_frames.md``): a cochain is an assignment
of values to faces, and the load-bearing vectorisation constraint (the
``(N_oct, ng, n_diag)`` wavefront batch is the unit of operation)
forbids a per-face Python object — it would reintroduce the per-cell
fold that caused a 10–20× regression. Per-face access is therefore a
**zero-copy view** (``face_view`` / ``face(axis)``), exactly as
:class:`BoundaryFlux` already does. The *interior* cell-face fluxes
are the named cochain
:class:`~orpheus.transport.fields.wavefront_flux.WavefrontFlux`
(``C^1_int``, ephemeral, flux-only) with ``face`` / ``edge_view``
zero-copy views; together with :class:`BoundaryFlux` (``C^1_∂``) they
biproduct-decompose the full face cochain ``C^1 = C^1_int ⊕ C^1_∂``.
See the theory page ``operator_algebra.rst`` § "The interior face-flux
cochain — WavefrontFlux".
"""

from __future__ import annotations

from orpheus.transport.fields.harmonic_moment_field import HarmonicMomentField
from orpheus.transport.fields.scalar_flux import ScalarFlux

__all__ = ["HarmonicMomentField", "ScalarFlux"]
