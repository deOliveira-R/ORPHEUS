---
name: question-and-source-vocabulary
description: Where ORPHEUS's question (posing, spectral map, adjoint arity) and source (intensional spec vs extensional sink) concepts live; identity/digest precedents (#405 P1, 2026-09-25)
metadata:
  type: project
---

Tree at 99ac3e66 (2026-09-25). Re-verify before repeating a count.

- Posing TYPES are L1 (`numerics/posing.py`, `pencil.py`); every production MINT is L3 (sn, homogeneous); L2 mints none. A brief saying "posings at L2" is wrong.
- Adjoint is an ARITY, not a flag: `EigenPosing.H()` nullary, `SourcePosing.H(detector)` unary. A stored forward/adjoint bool makes the source adjoint spellable without its detector.
- k vs α differ only by `SpectralMap` on one pencil; c-eigen needs a different pencil split (no production poser). `SpectralMap` holds lambdas: no content identity.
- FixedSource carries a σ point: `at(1)` physical (hub `source_posing`), `at(0)` pure transport (`solve_sn_fixed_source`, ruling F11 = the entry's choice).
- Intensional source precedent: `geometry/boundary/_source.py` `InflowSourceSpec.evaluate(space)`; projection bridge `AngularBoundarySourceSink.from_specs`. A bulk `Source` is its twin (X4). The volumetric path is extensional only (caller projects; `from_isotropic` = /Σw).
- MMS sources are phase-space Q(x, Ω, g) on the method's quadrature at cell centres; SymPy only derives/self-checks (srepr 0).
- The question is stored today as tags ON THE ANSWER: `CriticalSolution.(eigenvalue_kind, parameter_kind)`.
- Stable digest precedent: `Axis._structural_bytes` + blake2b. `Mixture._identity_key` is stable across seeds; only `hash` is salted. `Materials` is eq=False by design.
- Census traps: registry functions are re-exported from `orpheus/derivations/__init__`, so a resolver that skips re-exports undercounts `continuous_get` (8 vs 22). `_CASES` has 6 test-local homonyms. The regex `rossi` matches "crossing".

Related: [[reference-producer-landscape]].
