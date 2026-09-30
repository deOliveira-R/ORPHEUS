r"""Method-agnostic Sood/Forster/Parsons benchmark case registry.

The :mod:`sood_registry` package is the **single source of truth** for
benchmark case configurations from the Sood-family literature
(Sood, Forster & Parsons 2003, Atalay 1997, future KLL Tables, etc.).

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
  + ``geometry_kind: str`` (``"slab"`` / ``"sphere"`` /
  ``"cylinder"`` / ``"infinite"``). Use
  :meth:`La13511Case.to_geometry` to materialise a
  :class:`StructuredGeometry`; production solvers consume
  ``materials`` + a :class:`Mesh1D` built by a
  :class:`~orpheus.mesh.mesher.Mesher`.
* **Tabulated truth values**: :math:`k_{\rm eff}` / :math:`k_\infty`,
  flux ratios, critical dimensions, etc. — whatever the published
  reference tabulates.
* **Citations**: the problem's source on the case (``problem``) and the
  published values' sources on the truth (``truth.sources``), each a
  :class:`~orpheus.data.citation.Citation` whose key is in ``docs/refs.bib``.

Tests live in ``tests/gates/derivations/`` and import case + solver(s),
producing the value to compare. See e.g.
``tests/gates/derivations/test_fn_sood2003_kinf.py`` for the F_N consumer
or ``tests/gates/derivations/test_sood_registry_compatibility.py`` for the
production-protocol smoke gates.

Module layout
-------------

* :mod:`.sood2003` — the 47 Sood, Forster & Parsons cases, cited to
  the 2003 edition.
* :mod:`.atalay1997` — the 6 Atalay 1997 reflected-slab cases.
* :mod:`.case` — the case schema both are built on.
* :mod:`.builders` — case → ``(materials, mesh, params)`` helpers
  for production-solver consumers.

References
----------

* Sood, A., Forster, R.A. & Parsons, D.K. (2003), "Analytical
  benchmark test set for criticality code verification", *Progress in
  Nuclear Energy* 42(1), 55-106: the edition the cases cite.
* The same authors' 1999 report LA-13511 (Los Alamos National
  Laboratory) is the earlier edition, with other table, equation and
  reference numbers.
"""
from __future__ import annotations

from .case import La13511Case, La13511Truth
from .sood2003 import (
    # Phase A first slice (5)
    ALL_FIRST_SLICE,
    SOOD2003_CASES,
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
    ATALAY1997_CASES,
    ATALAY_SLAB_C130_R000_F0,
    ATALAY_SLAB_C130_R000_F010,
    ATALAY_SLAB_C130_R025_F0,
    ATALAY_SLAB_C130_R050_F0,
    ATALAY_SLAB_C130_R050_F010,
    ATALAY_SLAB_C130_R075_F0,
)
from .builders import build_cp_params, build_materials, build_mesh
from .cache import SoodResultCache, cache_info, clear_cache, sood_cache

__all__ = [
    # Core schema
    "La13511Case",
    "La13511Truth",
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
    "SOOD2003_CASES",
    # Builders / extractors
    "build_materials",
    "build_mesh",
    "build_cp_params",
    # Cache (Phase B4)
    "SoodResultCache",
    "sood_cache",
    "clear_cache",
    "cache_info",
    # Atalay 1997 reflected-slab catalogue (Wave 2-B)
    "ATALAY_SLAB_C130_R000_F0",
    "ATALAY_SLAB_C130_R025_F0",
    "ATALAY_SLAB_C130_R050_F0",
    "ATALAY_SLAB_C130_R075_F0",
    "ATALAY_SLAB_C130_R000_F010",
    "ATALAY_SLAB_C130_R050_F010",
    "ATALAY1997_CASES",
]
