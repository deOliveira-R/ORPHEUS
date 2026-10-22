---
name: transport-carrier-surface-typing
description: How transport fields read method-carrier surfaces (angular_bulk_space, moment_space, angular_trace, bc) without naming SNProblem; the P3 census of 2026-10-22
metadata:
  type: project
---

`[M]` 2026-10-22 at f672a2d5 (scripts `scratch/layer_dag/p3/`): the gate's
TYPE_CHECKING tolerance admitted exactly 7 imports, all transport -> sn/diffusion,
0 from numerics (control: a synthetic transport TC import is listed). 3 of 7 were
dead (0 AST uses: `_bases.py` DiffusionMesh, `angular_flux.py` and
`scalar_source_sink.py` SNProblem).

Durable shape: transport already had THREE spellings of "the carrier owns this
family's space": a runtime_checkable Protocol + isinstance (`_CarriesMomentSpace`),
a getattr-None door (`_face_space_of` for angular/scalar traces), and an `SNProblem`
annotation (`AngularField._space_for_problem`). `TransportMethod`
(transport/method.py) is the transport Protocol carrying `bc`; it deliberately
excludes trace spaces.

On Python 3.14 a runtime_checkable Protocol isinstance uses getattr_static: it does
NOT evaluate a cached_property (`angular_bulk_space` stays uncached), and costs
~0.25 us. On a bare `object.__new__` instance an instance-attribute member (`bc`)
reads absent, so probe Protocol conformance on a REAL instance.

`RadialCharacteristicField.source_from_angular` reads quad levels and the angular
closure: SN machinery living in transport (the P3 design moves it to
`sn/operators/radial_characteristic.py`). Related: [[geometry-numerics-package-cycle]].
