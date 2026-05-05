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
  + ``geometry_spec: GeometrySpec``. The same objects production
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
    # Core schema
    La13511Case,
    La13511Truth,
    Provenance,
    GeometrySpec,
    # Phase A first slice (5)
    ALL_FIRST_SLICE,
    LA13511_CASES,
    PU_2_0_IN,
    PUA_1_0_IN,
    UA_1_0_CY_STUB,
    UA_1_0_SL_STUB,
    UA_1_0_SP_STUB,
    # Phase B3 — 1G k_inf (8 + 3 anisotropic)
    PUB_1_0_IN,
    UA_1_0_IN,
    UB_1_0_IN,
    UC_1_0_IN,
    UD_1_0_IN,
    UD2O_1_0_IN,
    UE_1_0_IN,
    PU_1_1_IN,
    UD2OA_1_1_IN,
    UD2OB_1_1_IN,
    UD2OC_1_1_IN,
    # Phase B3 — 2G k_inf (7)
    U_2_0_IN,
    UAL_2_0_IN,
    URRA_2_0_IN,
    URRB_2_0_IN,
    URRC_2_0_IN,
    URRD_2_0_IN,
    UD2O_2_0_IN,
    # Phase B3 — 3G/6G k_inf (2)
    URR_3_0_IN,
    URR_6_0_IN,
    # Phase B3 — 1G bare-critical (5)
    PUA_1_0_SL,
    PUB_1_0_SL,
    UD2O_1_0_SL,
    PUB_1_0_SP,
    UD2O_1_0_SP,
    # Wave 2-C — 1G P_1 anisotropic bare-critical (5)
    PUA_1_1_SL,
    PUB_1_1_SL,
    UD2OA_1_1_SP,
    UD2OB_1_1_SP,
    UD2OC_1_1_SP,
    WIDE_SLICE_BARE_CRITICAL_1G_P1,
    # Phase B3 — STUBS (cylinder + 2G bare-critical)
    PUB_1_0_CY_STUB,
    UD2O_1_0_CY_STUB,
    PU_2_0_SL_STUB,
    PU_2_0_SP_STUB,
    U_2_0_SL_STUB,
    U_2_0_SP_STUB,
    UAL_2_0_SL_STUB,
    UAL_2_0_SP_STUB,
    URRA_2_0_SL_STUB,
    URRA_2_0_SP_STUB,
    UD2O_2_0_SL_STUB,
    UD2O_2_0_SP_STUB,
    # Slice tuples
    WIDE_SLICE_KINF,
    WIDE_SLICE_BARE_CRITICAL_1G,
    WIDE_SLICE_STUBS,
)
from .atalay1997 import (
    ATALAY_ALL_CASES,
    ATALAY_SLAB_CASES,
    ATALAY_SLAB_C130_R000_F0,
    ATALAY_SLAB_C130_R000_F010,
    ATALAY_SLAB_C130_R025_F0,
    ATALAY_SLAB_C130_R050_F0,
    ATALAY_SLAB_C130_R050_F010,
    ATALAY_SLAB_C130_R075_F0,
    ATALAY_SPHERE_C130_R000_F0,
    ATALAY_SPHERE_CASES,
)
from .builders import build_cp_params, build_materials, build_mesh
from .cache import SoodResultCache, cache_info, clear_cache, sood_cache
from .extractors import mixture_to_fn_arrays

__all__ = [
    # Core schema
    "La13511Case",
    "La13511Truth",
    "Provenance",
    "GeometrySpec",
    # Phase A first slice
    "PUA_1_0_IN",
    "PU_2_0_IN",
    "UA_1_0_SL_STUB",
    "UA_1_0_CY_STUB",
    "UA_1_0_SP_STUB",
    "ALL_FIRST_SLICE",
    # Phase B3 — k_inf cases
    "PUB_1_0_IN",
    "UA_1_0_IN",
    "UB_1_0_IN",
    "UC_1_0_IN",
    "UD_1_0_IN",
    "UD2O_1_0_IN",
    "UE_1_0_IN",
    "PU_1_1_IN",
    "UD2OA_1_1_IN",
    "UD2OB_1_1_IN",
    "UD2OC_1_1_IN",
    "U_2_0_IN",
    "UAL_2_0_IN",
    "URRA_2_0_IN",
    "URRB_2_0_IN",
    "URRC_2_0_IN",
    "URRD_2_0_IN",
    "UD2O_2_0_IN",
    "URR_3_0_IN",
    "URR_6_0_IN",
    # Phase B3 — bare-critical 1G
    "PUA_1_0_SL",
    "PUB_1_0_SL",
    "UD2O_1_0_SL",
    "PUB_1_0_SP",
    "UD2O_1_0_SP",
    # Wave 2-C — P_1 anisotropic bare-critical
    "PUA_1_1_SL",
    "PUB_1_1_SL",
    "UD2OA_1_1_SP",
    "UD2OB_1_1_SP",
    "UD2OC_1_1_SP",
    "WIDE_SLICE_BARE_CRITICAL_1G_P1",
    # Phase B3 — stubs
    "PUB_1_0_CY_STUB",
    "UD2O_1_0_CY_STUB",
    "PU_2_0_SL_STUB",
    "PU_2_0_SP_STUB",
    "U_2_0_SL_STUB",
    "U_2_0_SP_STUB",
    "UAL_2_0_SL_STUB",
    "UAL_2_0_SP_STUB",
    "URRA_2_0_SL_STUB",
    "URRA_2_0_SP_STUB",
    "UD2O_2_0_SL_STUB",
    "UD2O_2_0_SP_STUB",
    # Slice tuples
    "WIDE_SLICE_KINF",
    "WIDE_SLICE_BARE_CRITICAL_1G",
    "WIDE_SLICE_STUBS",
    # Top-level registry
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
    # Atalay 1997 reflected slab + sphere catalogue (Wave 2-B)
    "ATALAY_SLAB_C130_R000_F0",
    "ATALAY_SLAB_C130_R025_F0",
    "ATALAY_SLAB_C130_R050_F0",
    "ATALAY_SLAB_C130_R075_F0",
    "ATALAY_SLAB_C130_R000_F010",
    "ATALAY_SLAB_C130_R050_F010",
    "ATALAY_SPHERE_C130_R000_F0",
    "ATALAY_SLAB_CASES",
    "ATALAY_SPHERE_CASES",
    "ATALAY_ALL_CASES",
]
