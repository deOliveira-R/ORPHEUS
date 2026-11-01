r"""#246 — the named moment-axis predicates: positive + negative (anti-pattern #11).

SKELETON (test-architect, PRE-IMPL for #246). The predicates are NOT written.
These tests pin the CONTRACT and will ERROR (AttributeError) until
``DiscretizationSchemeBase.is_multi_moment`` and
``BulkField.has_spatial_moment_axis`` land (see the spec
``.claude/agent-memory/test-architect/issue_246_moment_axis_predicate_spec.md``).

Two predicates, two pairs:

* ``DiscretizationSchemeBase.is_multi_moment`` (= ``spatial_basis_per_axis > 1``)
  — the SCHEME-level UNCONDITIONAL truth used at the inner-walk sites. LD → True,
  DD → False, regardless of any field provenance.
* ``BulkField.spatial_moments_per_axis`` (the pre-existing moment-width count) —
  the FIELD-level provenance check. (#246 TRIMMED the speculative
  ``has_spatial_moment_axis`` boolean — zero production consumers, Pattern 6; the
  count already answers "does this field carry the moment axis".)
  PROVENANCE-DEPENDENT: ``> 1`` only on output that came through a moment-aware
  producer (the factor is set by the OUTPUT wrap, default 1 everywhere, NOT
  auto-read from the mesh — ``transport/fields/_bases.py:183-194``). A hand-built
  bare field on an LD mesh is ``1`` BY DESIGN — the load-bearing
  construct-general pin (P4'). The rows P3 / P4 read an SN solve's output and live in
  ``tests/gates/sn/solve/test_solutions_carry_the_moment_axis.py`` since layer DAG Phase 4.

vv-principles anti-pattern #11: a predicate needs a POSITIVE (real multi-moment →
True, MUST NOT mis-report) AND a NEGATIVE (real single-moment → False). Negative-
only validates the raising/return behaviour but NOT the invariant claim.

``foundation`` level (software invariant — a typed query over a data structure,
no equation ``:label:``; NO ``verifies``). ``-O``-safe: assert via
``pytest.fail`` on the bool, never bare ``assert``.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import StructuredGeometry
from orpheus.mesh import CellsByCount, Mesher
from orpheus.numerics.quadrature import Quadrature
from orpheus.transport.spatial import DiamondDifference, LinearDiscontinuous
from orpheus.transport.fields.angular_flux import AngularFlux


pytestmark = pytest.mark.foundation


def _require(predicate: bool, msg: str) -> None:
    r"""-O-safe boolean assert (a bare ``assert`` would be stripped under -O)."""
    if not predicate:
        pytest.fail(msg)


# ═══════════════════════════════════════════════════════════════════════
# P1/P2 — DiscretizationSchemeBase.is_multi_moment (scheme-level, UNCONDITIONAL)
# ═══════════════════════════════════════════════════════════════════════
#
# MUTATION CHECK (run once at impl, document the result in the PR): redefine
# is_multi_moment to a constant True → P2 must redden; constant False → P1 must
# redden. If a mutation does NOT redden, the body is not `spatial_basis_per_axis
# > 1` as claimed (Mode-10 activated-but-unconstrained).


def test_is_multi_moment_true_for_linear_discontinuous() -> None:
    r"""P1 (positive): LD is multi-moment, and the predicate equals the
    underlying ``spatial_basis_per_axis > 1`` (not a hardcoded True)."""
    ld = LinearDiscontinuous()
    _require(
        ld.is_multi_moment is True,
        "LinearDiscontinuous().is_multi_moment should be True",
    )
    _require(
        ld.is_multi_moment == (ld.spatial_basis_per_axis > 1),
        "is_multi_moment must equal spatial_basis_per_axis > 1 (single source)",
    )


def test_is_multi_moment_false_for_diamond_difference() -> None:
    r"""P2 (negative): DD is the slopeless cell-average closure → False."""
    dd = DiamondDifference()
    _require(
        dd.is_multi_moment is False,
        "DiamondDifference().is_multi_moment should be False",
    )
    _require(
        dd.is_multi_moment == (dd.spatial_basis_per_axis > 1),
        "is_multi_moment must equal spatial_basis_per_axis > 1 (single source)",
    )


@pytest.mark.skip(
    reason="Step scheme is a docstring stub today (no live class — see "
    "d5_nd_polymorphism_verification). UN-SKIP when Step lands; it must be "
    "is_multi_moment False (spatial_basis_per_axis == 1)."
)
def test_is_multi_moment_false_for_step() -> None:
    r"""P2' (negative): the Step scheme (slopeless) → False. Gated until Step
    is instantiable."""
    raise NotImplementedError  # pragma: no cover — un-skip when Step exists


# ═══════════════════════════════════════════════════════════════════════
# P3/P4 — BulkField.spatial_moments_per_axis (field-level, PROVENANCE-dependent)
# ═══════════════════════════════════════════════════════════════════════


def test_spatial_moment_axis_absent_on_hand_built_ld_field() -> None:
    r"""P4' (the CONSTRUCT-GENERAL pin — load-bearing): a hand-built bare field
    on an LD mesh LACKS the factor (``spatial_moments_per_axis == 1``), DESPITE
    the mesh being LD.

    This is INTENTIONAL: the spatial-moment axis is set by the OUTPUT
    wrap (``spatial_moments=per_axis``), default 1 everywhere, NOT auto-read from
    ``mesh.scheme`` (``transport/fields/_bases.py:183-194``). WITHOUT this pin a
    future change that auto-reads the scheme into the field default would pass
    silently and break LD byte-identity. This test is WHY the inner walk uses the
    SCHEME-level ``is_multi_moment``, not a field-level query.

    Build a bare width-1 field on the width-1 angular bulk (``angular_bulk_space(mesh, quadrature)``, NOT the LD
    closure's widened ``trial_space``) and assert the field does NOT carry the spatial-moment axis even though the
    closure IS multi-moment: widening is the CALLER's selection (layer DAG Phase 4: posed on a material mesh and the
    raw constructors; the SN solve's output rows P3 / P4 moved to
    ``tests/gates/sn/solve/test_solutions_carry_the_moment_axis.py``).
    """
    from orpheus.transport.mesh.material_mesh import MaterialMesh
    from orpheus.transport.phase_spaces import angular_bulk_space
    from tests.gates.transport._materials import placeholder_materials

    nx = 5
    mesh = Mesher(StructuredGeometry.from_homogeneous(1.0, 0)).partition(
        CellsByCount.uniform_width(nx),
    ).mesh
    quad = Quadrature.gauss_legendre(4)
    scheme = LinearDiscontinuous()
    material_mesh = MaterialMesh(mesh, placeholder_materials(ng=2))
    # Premise: the closure IS multi-moment — the field's False is NOT a trivial DD result; it documents the provenance
    # discipline.
    _require(scheme.is_multi_moment is True, "premise broken: the closure must be LD (multi-moment)")
    # Bare field — NO spatial_moments=, the construct-general default 1.
    values = np.zeros((quad.n_ordinates, material_mesh.ng, nx))
    bulk = AngularFlux(values=values, space=angular_bulk_space(material_mesh, quad))
    _require(
        bulk.spatial_moments_per_axis == 1,
        "a hand-built bare field on an LD mesh must NOT carry the moment axis "
        "(the factor is producer-set, not mesh-derived; construct-general pin)",
    )
