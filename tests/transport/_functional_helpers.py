r"""Shared fixtures + the SUT import guard for the #257 S5 ``Functional``
category verification suite.

S5 was PRE-IMPLEMENTATION at the time these specs were written — the
system-under-test (``orpheus.numerics.functional.Functional`` and the
concrete production-rate functional, now
:class:`~orpheus.transport.reaction_rate_functional.ReactionRateFunctional`)
did not yet exist on the tree. Every S5 test file routes its SUT import
through :func:`require_functional` / :func:`require_production_rate_functional`
so the files **collect** (green test session, no ``ImportError`` at
collection time) and **skip** with a clear reason if the production code
is absent. Now that it has landed, the skips resolve to real assertions
automatically — no edit to the test files is required.

The SUT-coupling is deliberately minimal (the test-architect brief
leaves the exact class layout to the method-implementer): the specs
probe the SUT's *structural surface* (does it satisfy the ``Functional``
Protocol? does it NOT satisfy ``LinearOperator``? does ``evaluate``
return the hand-derived value?) rather than its concrete construction
internals. The ONE construction assumption — the positional constructor
``ReactionRateFunctional(<CrossSectionField>)`` — is isolated
in :func:`build_production_rate_functional` so a different chosen layout
is a single-line change here, not a sweep across four files.

``foundation`` — these are software invariants on the algebra's TYPE
SURFACE (Protocol membership, category distinctness) plus an
L0 value-correctness check against a hand-derived reference. No
theory-page ``:label:`` is claimed (so no ``verifies(...)``).

vv Mode-8: the canonical ORPHEUS invocation is ``python -O``, which
strips bare ``assert`` to a NO-OP. Every structural gate in the S5
suite routes through :func:`require` (a function call, fires under
``-O``) or ``np.testing.*`` — NEVER a bare ``assert``.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D, Mesh2D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.geometry import SNMesh
from orpheus.transport.fields.cross_section_field import CrossSectionField
from tests.sn._test_helpers import placeholder_materials


# ───────────────────────────────────────────────────────────────────────
# vv Mode-8: -O-firing assertion (NOT a bare ``assert``).
# ───────────────────────────────────────────────────────────────────────


def require(condition: bool, message: str) -> None:
    """Fail with ``message`` if ``condition`` is false. Fires under ``-O``."""
    if not condition:
        pytest.fail(message)


# ───────────────────────────────────────────────────────────────────────
# Not-yet-written-SUT import guards (pytest.importorskip).
# ───────────────────────────────────────────────────────────────────────


def require_functional():
    """The ``Functional`` Protocol (S5) — skip the file if not yet landed."""
    mod = pytest.importorskip(
        "orpheus.numerics.functional",
        reason="#257 S5 PRE-IMPL: orpheus.numerics.functional.Functional "
        "not yet written by the method-implementer.",
    )
    if not hasattr(mod, "Functional"):
        pytest.skip(
            "#257 S5 PRE-IMPL: module exists but `Functional` symbol "
            "not yet defined."
        )
    return mod.Functional


def require_production_rate_functional():
    """The concrete production-rate functional — now ``ReactionRateFunctional``.

    The original ``ProductionRateFunctional`` (S5) was RETIRED:
    :class:`~orpheus.transport.reaction_rate_functional.ReactionRateFunctional`
    generalises it (production = ``ReactionRateFunctional(νΣf)``, absorption =
    ``ReactionRateFunctional(Σa)``). The surviving foil consumers
    (``test_functional_category`` / ``test_integral_kernel_category``) only
    need a concrete ``Functional`` carrying ``νΣf`` to foil against; the
    reaction-rate functional IS that object, so this probe now resolves the
    successor. The ``importorskip`` shape is kept (it always resolves on a
    landed tree) so the structural surface of these helpers is unchanged.
    """
    mod = pytest.importorskip(
        "orpheus.transport.reaction_rate_functional",
        reason="orpheus.transport.reaction_rate_functional.ReactionRateFunctional "
        "not importable.",
    )
    cls = getattr(mod, "ReactionRateFunctional", None)
    if cls is None:
        pytest.skip(
            "module exists but `ReactionRateFunctional` symbol not defined."
        )
    return cls


# ───────────────────────────────────────────────────────────────────────
# Meshes.
# ───────────────────────────────────────────────────────────────────────


def slab_mesh(nx: int = 4, ng: int = 2) -> SNMesh:
    """A 1-D slab, ``ng`` groups, ``nx`` cells, GL-4 quadrature."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def cartesian_2d_mesh(nx: int = 5, ny: int = 3, ng: int = 2) -> SNMesh:
    """The axis-discriminating regime: nx ≠ ny, ng = 2.

    A wrong contraction axis (group-vs-spatial variable-swap, mode #2)
    mis-scales or raises ONLY when the carrier has ≥2 complement axes of
    DIFFERENT sizes; a square mesh would silently agree.
    """
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 2.0, nx + 1),
        edges_y=np.linspace(0.0, 1.0, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_xmin=BC("vacuum"), bc_xmax=BC("vacuum"),
        bc_ymin=BC("vacuum"), bc_ymax=BC("vacuum"),
    )
    quad = Quadrature.level_symmetric(sn_order=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


# ───────────────────────────────────────────────────────────────────────
# Cross-section + flux fields. ASYMMETRIC, HETEROGENEOUS, ≥2G — so a
# group/variable swap (νΣf ↔ φ) or a wrong-axis contraction cannot null.
# ───────────────────────────────────────────────────────────────────────


def asymmetric_nu_sigma_f(ng: int, spatial_shape: tuple[int, ...]) -> np.ndarray:
    """νΣf ``(ng, *spatial)`` — distinct per group AND per cell, no symmetry.

    The leading group axis carries a strong per-group offset (group 0 is
    ~10× group 1) so the contraction ``Σ_g νΣf_g · φ_g`` is dominated by
    a different group than a νΣf↔φ swap would produce. The per-cell
    variation is a smooth ramp, distinct per group, so the result is a
    non-flat field in EVERY cell — a flat or symmetric νΣf would null the
    redistribution-style discrimination (H2).
    """
    base = np.arange(1, int(np.prod(spatial_shape)) + 1, dtype=float)
    base = base.reshape(spatial_shape)
    rows = []
    for g in range(ng):
        # group 0: ~[2.0 .. ], group 1: ~[0.2 .. ], etc. — strongly asymmetric.
        scale = 2.0 / (10.0**g)
        rows.append(scale * (1.0 + 0.13 * base + 0.07 * g))
    return np.stack(rows, axis=0)


def asymmetric_phi(ng: int, spatial_shape: tuple[int, ...]) -> np.ndarray:
    """φ ``(ng, *spatial)`` — distinct per group AND per cell, ≠ νΣf shape ramp.

    Built with a DIFFERENT per-cell variation pattern than
    :func:`asymmetric_nu_sigma_f` so that swapping the two arrays
    (νΣf↔φ, mode #2) changes the contraction value — a symmetric pair
    would agree under the swap and the test would be blind.
    """
    base = np.arange(1, int(np.prod(spatial_shape)) + 1, dtype=float)
    base = base.reshape(spatial_shape)
    rows = []
    for g in range(ng):
        # Inverted group trend (group 0 small, group 1 large) + a
        # cosine ripple → distinct from the νΣf ramp.
        scale = 0.5 * (g + 1)
        rows.append(scale * (3.0 - 0.05 * base + 0.4 * np.cos(base + g)))
    return np.stack(rows, axis=0)


def cross_section_field(values: np.ndarray, sn_mesh: SNMesh) -> CrossSectionField:
    """Wrap raw ``(ng, *spatial)`` into a typed CrossSectionField."""
    return CrossSectionField.from_mesh(values, sn_mesh)


def build_production_rate_functional(nu_sigma_f_field: CrossSectionField):
    """Construct the production-rate functional foil.

    THE single construction assumption for the surviving foil consumers —
    the positional constructor
    ``ReactionRateFunctional(<CrossSectionField>)`` (the retired
    ``ProductionRateFunctional(nu_sigma_f=...)`` is superseded). If the
    construction surface changes, change this ONE function, not the spec
    files.
    """
    cls = require_production_rate_functional()
    return cls(nu_sigma_f_field)


# ───────────────────────────────────────────────────────────────────────
# The STRUCTURALLY-INDEPENDENT hand-derived reference (B-leg correctness).
# ───────────────────────────────────────────────────────────────────────


def hand_derived_production_density(
    nu_sigma_f: np.ndarray, phi: np.ndarray,
) -> np.ndarray:
    r"""Per-cell fission emission density by an EXPLICIT Python double-loop.

    .. math::

        p(\vec r) \;=\; \sum_{g'} \nu\Sigma_{f,g'}(\vec r)\,\phi_{g'}(\vec r)

    Computed by an explicit nested Python ``for`` loop with a running
    scalar accumulator per cell — NOT ``numpy.einsum``, NOT
    ``.sum(axis=...)``, NOT a numpy broadcast-multiply. This is the
    structurally-independent reference of the B-leg: it shares NO
    reduction primitive with the production
    ``(self.right * x).sum(axis=self.axis, keepdims=True)`` extraction
    inside ``RankOneOperator.apply`` and NONE of the ORPHEUS operator
    algebra. A wrong-axis contraction, a group/variable swap, or a
    measure-fold in the production code disagrees with this loop.

    Returns the density ``(1, *spatial)`` (group axis collapsed,
    ``keepdims``-style leading-1) matching the production
    ``inner = (right * x).sum(axis=0, keepdims=True)`` shape. The S5
    spec asserts the SUT's output equals this — the SUT decides whether
    it returns ``(1, *spatial)`` or ``(*spatial,)``; the value
    comparison in the spec file squeezes both to a common shape.
    """
    ng = nu_sigma_f.shape[0]
    spatial_shape = nu_sigma_f.shape[1:]
    out = np.zeros(spatial_shape, dtype=float)
    # Explicit per-cell, per-group scalar accumulation. The flat index
    # walk over the spatial product keeps this a genuinely different
    # structural path (no numpy reduction) for ANY spatial rank.
    for flat in range(int(np.prod(spatial_shape))):
        idx = np.unravel_index(flat, spatial_shape)
        acc = 0.0
        for g in range(ng):
            acc += float(nu_sigma_f[(g, *idx)]) * float(phi[(g, *idx)])
        out[idx] = acc
    return out.reshape((1, *spatial_shape))


def squeeze_density(values: np.ndarray) -> np.ndarray:
    """Collapse a leading length-1 group axis if present, else return as-is.

    The SUT may return the density as ``(1, *spatial)`` (keepdims, like
    the production ``inner``) or ``(*spatial,)`` (squeezed). The value
    comparison is agnostic to that choice — both squeeze to ``(*spatial,)``.
    """
    arr = np.asarray(values)
    if arr.ndim >= 1 and arr.shape[0] == 1:
        return arr.reshape(arr.shape[1:])
    return arr
