r"""Method-agnostic Sood/Forster/Parsons benchmark case registry.

The :mod:`sood_registry` package is the **single source of truth** for
benchmark case configurations from the Sood-family literature
(LA-13511, LA-UR-03-1987 cylinder benchmarks, future KLL Tables, etc.).

Design intent
-------------

Before this package, each method (F_N, Variant α, etc.) carried its
own per-case data alongside its solvers. That coupling is the wrong
factoring: the same Sood case (XS, geometry, published reference
value) must be consumable by **every** method that wants to benchmark
itself against it — semi-analytical reference solvers (F_N, PS-1982),
production discrete solvers (CP, SN, MOC), Monte Carlo, etc.

Each :class:`La13511Case` carries:

* **Production-protocol XS + geometry**: ``materials: dict[int, Mixture]``
  + ``mesh_template: MeshTemplate``. The same objects production
  solvers (:func:`orpheus.cp.solver.solve_cp`,
  :func:`orpheus.sn.solver.solve_sn`) consume directly.
* **Tabulated truth values**: :math:`k_{\rm eff}` / :math:`k_\infty`,
  flux ratios, critical dimensions, etc. — whatever the published
  reference tabulates.
* **Provenance**: Sood table number + primary reference + notes.

Tests live in ``tests/derivations/`` and import case + solver(s),
producing the value to compare. See e.g.
``tests/derivations/test_fn_la13511_kinf.py`` for the F_N consumer
or ``tests/derivations/test_sood_registry_compatibility.py`` for the
production-protocol smoke gates.

Module layout
-------------

* :mod:`.la13511` — Sood/Forster/Parsons LA-13511 (1999) cases
  (5 first-slice cases; more added as Phase B lands).
* :mod:`.builders` — case → ``(materials, mesh, params)`` helpers
  for production-solver consumers.
* :mod:`.extractors` — :class:`Mixture` → numpy-array extractors for
  semi-analytical-reference consumers (F_N, etc.).

References
----------

* Sood, Forster & Parsons (1999), *Analytical Benchmark Test Set for
  Criticality Code Verification*, LA-13511, Los Alamos National
  Laboratory.
"""
from __future__ import annotations

from .la13511 import (
    ALL_FIRST_SLICE,
    LA13511_CASES,
    La13511Case,
    La13511Truth,
    MeshTemplate,
    PU_2_0_IN,
    PUA_1_0_IN,
    UA_1_0_CY_STUB,
    UA_1_0_SL_STUB,
    UA_1_0_SP_STUB,
)
from .builders import build_cp_params, build_materials, build_mesh
from .cache import SoodResultCache, cache_info, clear_cache, sood_cache
from .extractors import mixture_to_fn_arrays

__all__ = [
    # Core schema
    "La13511Case",
    "La13511Truth",
    "MeshTemplate",
    # Cases
    "PUA_1_0_IN",
    "PU_2_0_IN",
    "UA_1_0_SL_STUB",
    "UA_1_0_CY_STUB",
    "UA_1_0_SP_STUB",
    "ALL_FIRST_SLICE",
    "LA13511_CASES",
    # Builders / extractors
    "build_materials",
    "build_mesh",
    "build_cp_params",
    "mixture_to_fn_arrays",
    # Cache (Phase B4)
    "SoodResultCache",
    "sood_cache",
    "clear_cache",
    "cache_info",
]
