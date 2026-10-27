"""The test-layer gate: a test file imports only what its layer's package may import (layer DAG ruling 22).

The rule. A file under ``tests/gates/<dir>/`` (at any depth below it) whose
``<dir>`` names an ``orpheus`` package with a row in the production contract,
``FORBIDDEN_EDGES`` of ``tests/gates/test_layer_imports.py``, may import only
what that package may import. The contract is READ from that module, never
restated here, so a change to the layers moves both gates at once. An edge is
an import statement at any depth: module top level, a function or class body
(a lazy import), and an ``if TYPE_CHECKING:`` or ``if typing.TYPE_CHECKING:``
block. Two kinds of target are judged:

- ``orpheus.<pkg>``, judged as the package ``<pkg>``;
- ``tests.gates.<other_dir>`` (a helper or fixture module of another gate
  directory), judged as the layer of ``<other_dir>``: a ``numerics`` gate
  importing ``tests.gates.sn._fixtures`` reaches SN machinery as surely as
  one importing ``orpheus.sn``.

The scanner is the production linter's own (``_imports_of_source``), with the
importing package set to ``tests.gates`` so a relative import resolves under
the gate tree.

Exempt, because their directory has no row in ``FORBIDDEN_EDGES``:
``tests/gates/cross_method/`` (the cross-method comparisons, which import
every method by design), ``tests/gates/tools/``, ``tests/gates/_mutation/``
and the files at the root of ``tests/gates/`` (this one included). A
``tests.gates`` target that is not a layered directory (a root helper such as
``tests.gates._exterior``, or ``cross_method``) is not judged.

The contract's former gap, closed 2026-10-26 (``layer_dag.md`` ruling 29):
``FORBIDDEN_EDGES`` had no ``derivations`` (L0) entry in the rows of
``numerics``, ``data``, ``mesh`` and ``geometry``, and no ``plotting`` (L4)
entry below the values row, so a test there importing a derivation read as
allowed (``[M]`` 2026-10-25, qa: 22 files, numerics 15, data 6, mesh 1). Most
of those edges were the synthetic cross-section catalogue, which moved to
``orpheus.data.macro_xs.synthetic``; the rest are listed below.

An import by string is refused, so the scanner sees every edge it judges
(layer DAG ruling 27, ``test_no_layered_gate_file_imports_a_project_module_by_string``): in a layered directory, a call to
``import_module`` or ``__import__`` whose first argument is a constant
string, or an f-string whose constant head names ``orpheus`` or
``tests.gates.``, or a module-level name bound once to such a string
(``_MODULE = "orpheus.numerics.enclosure"``), is red; the spelling is an
import statement of the MODULE, its attribute read at call time. What stays
dynamic, and is not judged: a walk over discovered modules (any other
variable or an expression argument, or an f-string that begins with a
substitution, ``f"{package.__name__}.{name}"``);
a synthetic package written at test time (the same shapes); an absence check,
spelled ``importlib.util.find_spec(name) is None``, since an import statement
cannot assert that a module is absent; the standard library; and the root
files.

Not seen, the declared blind spots: a walk's targets (discovered, not
written), and so also a known name reached through a parameter, a loop
variable, a container or an expression (``for name in ("orpheus.a", ...):
import_module(name)``): ``[M]`` 2026-10-26, 20 of the 87 calls at
``3efcd39d`` had that shape and 0 after ruling 27, so the rule, not this
gate, keeps it at 0; an import spelled inside a string literal (a subprocess source, or
a synthetic package's module text); and the TRANSITIVE cone of a root helper
(a ``numerics`` gate importing ``tests.gates._exterior``, which imports
``orpheus.mesh``). The census of 2026-10-25
(``scratch/layer_dag/refiling/census/test_imports.json``) counts the second
separately and lists the third under ``forbidden_via_root_helper``.

Today's violations are listed EXACTLY in ``OPEN_TEST_LAYER_VIOLATIONS``,
keyed by (file, imported package), the way ``_OPEN_INTRA_PACKAGE_CYCLES``
lists the open module cycles: a new violation is red, and a listed row that no
longer violates (or whose file is gone) is red too, so the list only shrinks.
A test move of the re-filing campaign (ruling 23) is measured as rows deleted
from it.

The second gate here is the ``rests_on`` resolution: every target of a
``rests_on`` marker names a collected test (``vv-principles``, "A capability's
tests form a ladder"). No other check exists (the census: 0 readers of the
marker in ``tests/_harness``, the conftests and ``tools/``).
"""
from __future__ import annotations

import ast
import json
import os
import pathlib
import subprocess
import sys
from collections.abc import Iterable, Mapping

import pytest

from tests._harness.rests_on_manifest import MANIFEST_ENV
from tests.gates.test_layer_imports import FORBIDDEN_EDGES, ORPHEUS_ROOT, _imports_of_source

REPO_ROOT = ORPHEUS_ROOT.parent
GATES_ROOT = REPO_ROOT / "tests" / "gates"

#: The package a gate file's path is relative to, for relative-import resolution.
_GATES_PACKAGE: tuple[str, ...] = ("tests", "gates")


def _layer_of(rel: pathlib.PurePosixPath) -> str | None:
    """The layer a gate file is judged as: its first directory under ``tests/gates/``, when it has a contract row."""
    if len(rel.parts) < 2:
        return None
    directory = rel.parts[0]
    return directory if directory in FORBIDDEN_EDGES else None


def _target_layer(module_name: str) -> tuple[str, str] | None:
    """``(imported package, layer)`` of an import target, ``None`` when the target is not judged.

    ``orpheus.sn.solver`` is ``("orpheus.sn", "sn")``; ``tests.gates.sn._fixtures``
    is ``("tests.gates.sn", "sn")``; ``tests.gates._exterior`` (a root helper) and
    ``tests.gates.cross_method`` (an exempt directory) are not judged.
    """
    parts = module_name.split(".")
    if parts[0] == "orpheus" and len(parts) >= 2:
        return f"orpheus.{parts[1]}", parts[1]
    if parts[:2] == ["tests", "gates"] and len(parts) >= 3 and parts[2] in FORBIDDEN_EDGES:
        return f"tests.gates.{parts[2]}", parts[2]
    return None


def _test_layer_violations(rel: pathlib.PurePosixPath, src: str) -> set[tuple[str, str]]:
    """The (file, imported package) violations of one gate file, ``rel`` relative to ``tests/gates/``."""
    layer = _layer_of(rel)
    if layer is None:
        return set()
    forbidden = FORBIDDEN_EDGES[layer]
    found: set[tuple[str, str]] = set()
    for module_name, _is_tc in _imports_of_source(rel, src, _GATES_PACKAGE):
        target = _target_layer(module_name)
        if target is not None and target[1] in forbidden:
            found.add((rel.as_posix(), target[0]))
    return found


def _gate_files() -> list[pathlib.Path]:
    return sorted(p for p in GATES_ROOT.rglob("*.py") if "__pycache__" not in p.parts)


def _list_diff(
    found: set[tuple[str, str]], listed: frozenset[tuple[str, str]], existing: set[str],
) -> tuple[list[str], list[str]]:
    """(new, stale) rows: found and not listed; listed and not found (the file named, or gone)."""
    new = [f"{file} imports {package}" for file, package in sorted(found - listed)]
    stale = [
        f"{file} imports {package}" + ("" if file in existing else " (no such file)")
        for file, package in sorted(listed - found)
    ]
    return new, stale


# Today's violations, exact: minted by this gate's scanner on 2026-10-25, after the four Phase 1 moves
# of the re-filing campaign (geometry -> mesh, geometry -> numerics, numerics -> cross_method,
# numerics -> transport), and re-minted on 2026-10-26 when the contract closed its derivations and plotting
# gap (ruling 29): the synthetic catalogue's numerics importers now read numerics -> data. The list only shrinks.
OPEN_TEST_LAYER_VIOLATIONS: frozenset[tuple[str, str]] = frozenset(
    {
        ("cp/test_cp_vacuum_not_above_white.py", "orpheus.sn"),
        ("data/test_cells.py", "orpheus.mesh"),
        ("data/test_cells.py", "orpheus.specification"),
        ("data/test_cells.py", "orpheus.transport"),
        ("data/test_cells.py", "tests.gates.specification"),
        ("data/test_content_identity_data.py", "orpheus.mesh"),
        ("data/test_content_identity_data.py", "orpheus.transport"),
        ("data/test_n2n_yield_convention.py", "orpheus.transport"),
        ("derivations/_characteristic_ladders.py", "tests.gates.sn"),
        ("derivations/test_characteristic_albedos.py", "tests.gates.sn"),
        ("derivations/test_characteristic_convergence.py", "tests.gates.sn"),
        ("derivations/test_characteristic_reference.py", "tests.gates.sn"),
        ("derivations/test_characteristic_system.py", "tests.gates.sn"),
        ("derivations/test_dsa_production_tie.py", "orpheus.sn"),
        ("derivations/test_dsa_production_tie.py", "orpheus.transport"),
        ("derivations/test_reference_pairing_loop.py", "orpheus.transport"),
        ("derivations/test_sn_mms_ld_2d_stress_symbolic.py", "orpheus.sn"),
        ("derivations/test_sn_mms_ld_2d_stress_symbolic.py", "orpheus.transport"),
        ("derivations/test_sn_mms_ld_2d_stress_symbolic.py", "tests.gates.sn"),
        ("derivations/test_sood_registry_compatibility.py", "orpheus.cp"),
        ("derivations/test_sood_registry_compatibility.py", "orpheus.sn"),
        ("diffusion/test_augmented_mesh.py", "orpheus.sn"),
        ("diffusion/test_solver.py", "orpheus.homogeneous"),
        ("geometry/test_angular_chart.py", "orpheus.numerics"),
        ("geometry/test_face_roles.py", "orpheus.transport"),
        ("geometry/test_geometry_extent.py", "tests.gates.specification"),
        ("geometry/test_kernel_corroboration.py", "orpheus.mc"),
        ("geometry/test_kernel_corroboration.py", "orpheus.moc"),
        ("geometry/test_manifold.py", "orpheus.numerics"),
        ("geometry/test_symmetry.py", "orpheus.numerics"),
        ("geometry/test_transformation.py", "orpheus.numerics"),
        ("mc/test_cross_verification.py", "orpheus.cp"),
        ("moc/test_verification.py", "orpheus.cp"),
        ("numerics/_axis_marginal_lift_support.py", "orpheus.data"),
        ("numerics/_axis_marginal_lift_support.py", "orpheus.mesh"),
        ("numerics/_axis_marginal_lift_support.py", "orpheus.sn"),
        ("numerics/_axis_marginal_lift_support.py", "orpheus.transport"),
        ("numerics/test_adjointability_is_honest.py", "orpheus.diffusion"),
        ("numerics/test_adjointability_is_honest.py", "orpheus.sn"),
        ("numerics/test_adjointability_is_honest.py", "orpheus.transport"),
        ("numerics/test_adjointability_is_honest.py", "tests.gates.diffusion"),
        ("numerics/test_adjointability_is_honest.py", "tests.gates.sn"),
        ("numerics/test_angular_symmetry.py", "orpheus.data"),
        ("numerics/test_angular_symmetry.py", "orpheus.mesh"),
        ("numerics/test_angular_symmetry.py", "orpheus.sn"),
        ("numerics/test_angular_symmetry.py", "orpheus.transport"),
        ("numerics/test_angular_trace_space.py", "orpheus.mesh"),
        ("numerics/test_assembly_functor.py", "orpheus.derivations"),
        ("numerics/test_assembly_functor.py", "orpheus.diffusion"),
        ("numerics/test_assembly_functor.py", "orpheus.homogeneous"),
        ("numerics/test_assembly_functor.py", "orpheus.transport"),
        ("numerics/test_assembly_functor.py", "tests.gates.diffusion"),
        ("numerics/test_assembly_functor.py", "tests.gates.homogeneous"),
        ("numerics/test_assembly_functor.py", "tests.gates.transport"),
        ("numerics/test_axis.py", "orpheus.data"),
        ("numerics/test_axis_identity.py", "orpheus.sn"),
        ("numerics/test_axis_identity.py", "orpheus.transport"),
        ("numerics/test_axis_identity.py", "tests.gates.sn"),
        ("numerics/test_axis_marginal.py", "orpheus.mesh"),
        ("numerics/test_axis_marginal.py", "orpheus.sn"),
        ("numerics/test_axis_marginal.py", "orpheus.transport"),
        ("numerics/test_axis_marginal.py", "tests.gates.sn"),
        ("numerics/test_axis_metric_is_a_derived_object.py", "tests.gates.sn"),
        ("numerics/test_basis_domain.py", "orpheus.data"),
        ("numerics/test_basis_domain.py", "orpheus.mesh"),
        ("numerics/test_basis_domain.py", "orpheus.sn"),
        ("numerics/test_dual_metric.py", "orpheus.transport"),
        ("numerics/test_gauge.py", "orpheus.data"),
        ("numerics/test_gauge.py", "orpheus.homogeneous"),
        ("numerics/test_gauge.py", "orpheus.mesh"),
        ("numerics/test_gauge.py", "orpheus.sn"),
        ("numerics/test_gauge.py", "orpheus.transport"),
        ("numerics/test_indicator_basis.py", "orpheus.mesh"),
        ("numerics/test_interval_legendre_basis.py", "orpheus.transport"),
        ("numerics/test_inverse_metric_operator.py", "orpheus.data"),
        ("numerics/test_inverse_metric_operator.py", "orpheus.mesh"),
        ("numerics/test_inverse_metric_operator.py", "orpheus.sn"),
        ("numerics/test_inverse_universal.py", "orpheus.sn"),
        ("numerics/test_inverse_universal.py", "tests.gates.sn"),
        ("numerics/test_iteration.py", "orpheus.data"),
        ("numerics/test_iteration.py", "orpheus.derivations"),
        ("numerics/test_iteration.py", "orpheus.mesh"),
        ("numerics/test_iteration.py", "orpheus.sn"),
        ("numerics/test_iteration.py", "orpheus.transport"),
        ("numerics/test_iteration.py", "tests.gates.sn"),
        ("numerics/test_linear_combination.py", "orpheus.sn"),
        ("numerics/test_linear_combination.py", "orpheus.transport"),
        ("numerics/test_linear_combination.py", "tests.gates.sn"),
        ("numerics/test_matrix_inverse_operator.py", "orpheus.data"),
        ("numerics/test_matrix_inverse_operator.py", "orpheus.derivations"),
        ("numerics/test_matrix_inverse_operator.py", "orpheus.transport"),
        ("numerics/test_matrix_inverse_operator.py", "tests.gates.transport"),
        ("numerics/test_moment_head_axis_built_premise.py", "orpheus.transport"),
        ("numerics/test_outcome.py", "orpheus.data"),
        ("numerics/test_outcome.py", "orpheus.homogeneous"),
        ("numerics/test_outcome.py", "orpheus.mesh"),
        ("numerics/test_outcome.py", "orpheus.sn"),
        ("numerics/test_outcome.py", "orpheus.transport"),
        ("numerics/test_pencil.py", "orpheus.data"),
        ("numerics/test_pencil.py", "orpheus.homogeneous"),
        ("numerics/test_pencil.py", "orpheus.sn"),
        ("numerics/test_pencil.py", "tests.gates.sn"),
        ("numerics/test_posing.py", "orpheus.data"),
        ("numerics/test_posing.py", "orpheus.homogeneous"),
        ("numerics/test_registry_mixin.py", "orpheus.sn"),
        ("numerics/test_registry_mixin.py", "orpheus.transport"),
        ("numerics/test_retraction_adjoint_is_the_pullback.py", "orpheus.data"),
        ("numerics/test_retraction_adjoint_is_the_pullback.py", "orpheus.mesh"),
        ("numerics/test_retraction_adjoint_is_the_pullback.py", "orpheus.sn"),
        ("numerics/test_retraction_adjoint_is_the_pullback.py", "orpheus.transport"),
        ("numerics/test_si_diagnostic_trajectory.py", "orpheus.data"),
        ("numerics/test_si_diagnostic_trajectory.py", "orpheus.mesh"),
        ("numerics/test_si_diagnostic_trajectory.py", "orpheus.sn"),
        ("numerics/test_si_diagnostic_trajectory.py", "orpheus.transport"),
        ("numerics/test_space_of_axes.py", "orpheus.data"),
        ("numerics/test_space_of_axes.py", "orpheus.diffusion"),
        ("numerics/test_space_of_axes.py", "orpheus.homogeneous"),
        ("numerics/test_space_of_axes.py", "orpheus.mesh"),
        ("numerics/test_space_of_axes.py", "orpheus.sn"),
        ("numerics/test_space_of_axes.py", "orpheus.transport"),
        ("numerics/test_space_of_axes.py", "tests.gates.sn"),
        ("numerics/test_space_of_axes.py", "tests.gates.transport"),
        ("numerics/test_symbolic_steps.py", "orpheus.mesh"),
        ("numerics/test_symbolic_steps.py", "tests.gates.sn"),
        ("numerics/test_tensor_product_metric_is_factored.py", "orpheus.transport"),
        ("numerics/test_tensor_product_metric_is_factored.py", "tests.gates.sn"),
        ("numerics/test_traced_memo_clients.py", "orpheus.data"),
        ("numerics/test_traced_memo_clients.py", "orpheus.derivations"),
        ("numerics/test_traced_memo_clients.py", "orpheus.specification"),
        ("numerics/test_traced_memo_findings.py", "orpheus.data"),
        ("numerics/test_traced_memo_ownership.py", "orpheus.derivations"),
        ("numerics/test_traced_memo_source_index.py", "orpheus.derivations"),
        ("numerics/test_traced_memo_source_index.py", "orpheus.specification"),
        ("sn/mesh/test_hub_owns_the_moment_space.py", "orpheus.diffusion"),
        ("sn/operators/test_isotropic_fission.py", "orpheus.diffusion"),
        ("sn/operators/test_space_content_witnesses.py", "orpheus.diffusion"),
        ("sn/solve/test_subcritical_multiplying_source.py", "orpheus.homogeneous"),
        ("sn/verification/analytical/test_cp_standoff_curvilinear.py", "orpheus.cp"),
        ("specification/test_exterior_on_the_specification.py", "orpheus.transport"),
        ("transport/_functional_helpers.py", "orpheus.sn"),
        ("transport/_functional_helpers.py", "tests.gates.sn"),
        ("transport/_wave2_s4_helpers.py", "tests.gates.diffusion"),
        ("transport/fields/test_angular_boundary_flux.py", "orpheus.sn"),
        ("transport/fields/test_angular_boundary_flux.py", "tests.gates.sn"),
        ("transport/fields/test_angular_boundary_source_sink_residual.py", "orpheus.sn"),
        ("transport/fields/test_angular_boundary_source_sink_residual.py", "tests.gates.sn"),
        ("transport/fields/test_angular_flux.py", "orpheus.sn"),
        ("transport/fields/test_angular_flux.py", "tests.gates.sn"),
        ("transport/fields/test_coefficient_fields.py", "orpheus.sn"),
        ("transport/fields/test_coefficient_fields.py", "tests.gates.sn"),
        ("transport/fields/test_facefield_hierarchy.py", "tests.gates.sn"),
        ("transport/fields/test_field_space_identity.py", "orpheus.sn"),
        ("transport/fields/test_field_space_identity.py", "tests.gates.sn"),
        ("transport/fields/test_flux_vector_algebra.py", "orpheus.sn"),
        ("transport/fields/test_flux_vector_algebra.py", "tests.gates.sn"),
        ("transport/fields/test_harmonic_moment_flux.py", "orpheus.sn"),
        ("transport/fields/test_harmonic_moment_flux.py", "tests.gates.sn"),
        ("transport/fields/test_meshless_construction.py", "orpheus.sn"),
        ("transport/fields/test_meshless_construction.py", "tests.gates.sn"),
        ("transport/fields/test_scalar_boundary_flux.py", "orpheus.diffusion"),
        ("transport/fields/test_spatial_moment_field_space.py", "orpheus.sn"),
        ("transport/fields/test_spatial_moment_field_space.py", "tests.gates.sn"),
        ("transport/frames/test_frame_mints_refuse_a_moment_space.py", "tests.gates.sn"),
        ("transport/frames/test_harmonic_frame.py", "orpheus.sn"),
        ("transport/frames/test_harmonic_frame.py", "tests.gates.sn"),
        ("transport/frames/test_moment_metric_fork_premise.py", "tests.gates.sn"),
        ("transport/frames/test_moment_space_is_read_off_the_frame.py", "orpheus.sn"),
        ("transport/frames/test_moment_space_is_read_off_the_frame.py", "tests.gates.sn"),
        ("transport/residuals/test_typed_residuals.py", "orpheus.sn"),
        ("transport/residuals/test_typed_residuals.py", "tests.gates.sn"),
        ("transport/spatial/test_cell_basis_head.py", "orpheus.sn"),
        ("transport/spatial/test_cell_basis_head.py", "tests.gates.sn"),
        ("transport/spatial/test_face_transmission_xverif.py", "orpheus.sn"),
        ("transport/spatial/test_ld_slope_frame.py", "orpheus.sn"),
        ("transport/spatial/test_ld_ubld_primitive.py", "orpheus.sn"),
        ("transport/spatial/test_ld_ubld_symbolic.py", "orpheus.sn"),
        ("transport/spatial/test_linear_discontinuous.py", "orpheus.sn"),
        ("transport/spatial/test_moment_axis_predicates.py", "orpheus.sn"),
        ("transport/spatial/test_moment_axis_predicates.py", "tests.gates.sn"),
        ("transport/spatial/test_spatial_moment_tail_is_the_schemes_axis.py", "orpheus.sn"),
        ("transport/spatial/test_spatial_moment_tail_is_the_schemes_axis.py", "tests.gates.sn"),
        ("transport/test_admission_matrix.py", "tests.gates.diffusion"),
        ("transport/test_admission_matrix.py", "tests.gates.sn"),
        ("transport/test_angular_lift.py", "orpheus.sn"),
        ("transport/test_boundary_source_from_specs.py", "orpheus.sn"),
        ("transport/test_boundary_source_from_specs.py", "tests.gates.sn"),
        ("transport/test_bulk_lift.py", "orpheus.diffusion"),
        ("transport/test_composite.py", "orpheus.diffusion"),
        ("transport/test_field_units.py", "orpheus.sn"),
        ("transport/test_field_units.py", "tests.gates.sn"),
        ("transport/test_full_field.py", "orpheus.sn"),
        ("transport/test_full_field.py", "tests.gates.sn"),
        ("transport/test_integral_kernel_category.py", "orpheus.sn"),
        ("transport/test_integrated_reaction_rate.py", "orpheus.sn"),
        ("transport/test_integrated_reaction_rate.py", "tests.gates.cp"),
        ("transport/test_integrated_reaction_rate.py", "tests.gates.sn"),
        ("transport/test_kernels.py", "orpheus.diffusion"),
        ("transport/test_kernels.py", "orpheus.sn"),
        ("transport/test_kernels.py", "tests.gates.sn"),
        ("transport/test_material_field.py", "tests.gates.sn"),
        ("transport/test_material_mesh.py", "orpheus.sn"),
        ("transport/test_material_mesh_admission.py", "orpheus.diffusion"),
        ("transport/test_material_mesh_admission.py", "orpheus.sn"),
        ("transport/test_material_mesh_admission.py", "tests.gates.sn"),
        ("transport/test_method.py", "orpheus.diffusion"),
        ("transport/test_method.py", "orpheus.sn"),
        ("transport/test_method_mesh_binding.py", "orpheus.diffusion"),
        ("transport/test_method_mesh_binding.py", "orpheus.sn"),
        ("transport/test_multiplication_operator.py", "orpheus.sn"),
        ("transport/test_multiplication_operator.py", "tests.gates.sn"),
        ("transport/test_radial_characteristic_field.py", "orpheus.sn"),
        ("transport/test_radial_characteristic_field.py", "tests.gates.sn"),
        ("transport/test_tier2_equivalence_s_family.py", "orpheus.sn"),
        ("transport/test_tier2_equivalence_s_family.py", "tests.gates.sn"),
        ("transport/test_timed_full_field.py", "orpheus.sn"),
        ("transport/test_timed_full_field.py", "tests.gates.sn"),
    }
)


# ---------------------------------------------------------------------------
# The positive controls: the predicate on synthetic sources (X1)
# ---------------------------------------------------------------------------

_FORBIDDEN_SOURCES = [
    pytest.param("from orpheus.sn.solver import solve_sn\n", "orpheus.sn", id="top"),
    pytest.param("import orpheus.sn.solver\n", "orpheus.sn", id="import"),
    pytest.param("from orpheus import sn\n", "orpheus.sn", id="from-orpheus-import-package"),
    pytest.param("def f():\n    from orpheus.sn.solver import solve_sn\n", "orpheus.sn", id="lazy"),
    pytest.param("class C:\n    from orpheus.sn import solver\n", "orpheus.sn", id="class-body"),
    pytest.param(
        "from typing import TYPE_CHECKING\nif TYPE_CHECKING:\n    from orpheus.sn.problem import SNProblem\n",
        "orpheus.sn", id="type-checking",
    ),
    pytest.param(
        "import typing\nif typing.TYPE_CHECKING:\n    from orpheus.sn.problem import SNProblem\n",
        "orpheus.sn", id="typing-type-checking",
    ),
    pytest.param("from tests.gates.sn._fixtures import fixture\n", "tests.gates.sn", id="helper"),
    pytest.param("from tests.gates import sn\n", "tests.gates.sn", id="from-gates-import-dir"),
    pytest.param("from ..sn import _fixtures\n", "tests.gates.sn", id="relative-helper"),
]


@pytest.mark.foundation
@pytest.mark.parametrize(("source", "package"), _FORBIDDEN_SOURCES)
def test_the_predicate_finds_each_import_shape(source: str, package: str) -> None:
    """A ``numerics`` gate importing SN, at every import kind and through a helper of the ``sn`` directory,
    is a violation; the same source in an ``sn`` gate is clean (the negative leg), and so is the same source
    in the exempt ``cross_method`` directory."""
    assert "sn" in FORBIDDEN_EDGES["numerics"], sorted(FORBIDDEN_EDGES["numerics"])
    rel = pathlib.PurePosixPath("numerics/test_probe.py")
    assert _test_layer_violations(rel, source) == {(rel.as_posix(), package)}
    assert _test_layer_violations(pathlib.PurePosixPath("sn/test_probe.py"), source) == set()
    assert _test_layer_violations(pathlib.PurePosixPath("cross_method/test_probe.py"), source) == set()
    assert _test_layer_violations(pathlib.PurePosixPath("test_probe.py"), source) == set()


@pytest.mark.foundation
def test_a_nested_file_is_judged_as_its_first_directory() -> None:
    """``tests/gates/sn/regression/x.py`` is an ``sn`` gate; ``tests/gates/numerics/sub/x.py`` a ``numerics`` one."""
    source = "from orpheus.sn.solver import solve_sn\n"
    assert _test_layer_violations(pathlib.PurePosixPath("sn/regression/test_probe.py"), source) == set()
    deep = pathlib.PurePosixPath("numerics/sub/test_probe.py")
    assert _test_layer_violations(deep, source) == {(deep.as_posix(), "orpheus.sn")}
    # a relative import climbs from the nested package: `...sn` from numerics/sub/ is tests.gates.sn
    assert _test_layer_violations(deep, "from ...sn import _fixtures\n") == {(deep.as_posix(), "tests.gates.sn")}


@pytest.mark.foundation
def test_the_diff_reports_a_new_and_a_stale_row() -> None:
    found = {("numerics/test_a.py", "orpheus.sn")}
    listed = frozenset({("numerics/test_b.py", "orpheus.sn"), ("numerics/test_gone.py", "orpheus.sn")})
    new, stale = _list_diff(found, listed, {"numerics/test_a.py", "numerics/test_b.py"})
    assert new == ["numerics/test_a.py imports orpheus.sn"]
    assert stale == ["numerics/test_b.py imports orpheus.sn", "numerics/test_gone.py imports orpheus.sn (no such file)"]
    assert _list_diff(found, frozenset(found), {"numerics/test_a.py"}) == ([], [])


# ---------------------------------------------------------------------------
# The tree
# ---------------------------------------------------------------------------


@pytest.mark.foundation
def test_every_gate_file_imports_only_its_layer_but_the_listed_rows() -> None:
    files = _gate_files()
    rels = [pathlib.PurePosixPath(p.relative_to(GATES_ROOT).as_posix()) for p in files]
    judged = [rel for rel in rels if _layer_of(rel) is not None]
    # The denominators are stated, so an empty walk cannot read as a clean tree.
    assert len(files) > 500, len(files)
    assert len(judged) > 400, len(judged)
    assert any(rel.parts[0] == "cross_method" for rel in rels), "the walk reached no exempt directory"
    assert any(len(rel.parts) > 2 for rel in judged), "the walk reached no nested gate file"
    found: set[tuple[str, str]] = set()
    for path, rel in zip(files, rels):
        found |= _test_layer_violations(rel, path.read_text(encoding="utf-8"))
    new, stale = _list_diff(found, OPEN_TEST_LAYER_VIOLATIONS, {rel.as_posix() for rel in rels})
    assert not new, (
        f"{len(new)} gate file(s) import above their layer (move the test to the layer its imports name, "
        "or re-base it onto fixtures of its own layer): " + "; ".join(new)
    )
    assert not stale, (
        f"{len(stale)} listed row(s) no longer violate; retire them from OPEN_TEST_LAYER_VIOLATIONS: "
        + "; ".join(stale)
    )


# ---------------------------------------------------------------------------
# No import by string (layer DAG ruling 27)
# ---------------------------------------------------------------------------

_RULING_27 = (
    "an import whose target is known when the test is written is an import STATEMENT, so this gate sees it "
    "(layer_dag.md ruling 27): write `from orpheus.numerics import enclosure` and read `enclosure.Enclosure` at "
    "call time (a monkeypatch of the module attribute still reaches every row); an absence check is "
    "`importlib.util.find_spec(name) is None`"
)


def _is_dynamic_import(func: ast.expr) -> bool:
    """``import_module`` (on any receiver, or bare) or the builtin ``__import__``."""
    if isinstance(func, ast.Attribute):
        return func.attr == "import_module"
    return isinstance(func, ast.Name) and func.id in {"import_module", "__import__"}


def _leading_literal(argument: ast.expr) -> str | None:
    """The module name a call spells at writing time: a constant string, or the constant head of an f-string;
    ``None`` for a variable, an expression, or an f-string that begins with a substitution (a walk's name)."""
    if isinstance(argument, ast.Constant) and isinstance(argument.value, str):
        return argument.value
    if isinstance(argument, ast.JoinedStr) and argument.values:
        head = argument.values[0]
        if isinstance(head, ast.Constant) and isinstance(head.value, str):
            return head.value
    return None


def _module_constants(tree: ast.Module) -> dict[str, ast.expr]:
    """The module-level names bound once to a constant string (``_MODULE = "orpheus.numerics.enclosure"``): a name
    handed to ``import_module`` is a string import when it is one of these, and a walk's when it is not."""
    bound: dict[str, list[ast.expr]] = {}
    for statement in tree.body:
        if isinstance(statement, ast.Assign):
            for target in statement.targets:
                if isinstance(target, ast.Name):
                    bound.setdefault(target.id, []).append(statement.value)
        elif isinstance(statement, ast.AnnAssign) and isinstance(statement.target, ast.Name) and statement.value is not None:
            bound.setdefault(statement.target.id, []).append(statement.value)
    return {name: values[0] for name, values in bound.items() if len(values) == 1 and _leading_literal(values[0]) is not None}


def _string_imports(rel: pathlib.PurePosixPath, src: str) -> list[str]:
    """``file:line: argument`` of every string import of a project module in a layered gate file."""
    if _layer_of(rel) is None:
        return []
    tree = ast.parse(src)
    constants = _module_constants(tree)
    found: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Call) and node.args and _is_dynamic_import(node.func):
            argument = node.args[0]
            if isinstance(argument, ast.Name) and argument.id in constants:
                argument = constants[argument.id]
            literal = _leading_literal(argument)
            if literal is not None and (literal.split(".")[0] == "orpheus" or literal.startswith("tests.gates.")):
                found.append(f"{rel.as_posix()}:{node.lineno}: {ast.unparse(node.args[0])}")
    return found


_STRING_IMPORTS = [
    pytest.param('import importlib\nimportlib.import_module("orpheus.cp.solver")\n', id="literal"),
    pytest.param('import importlib as _il\ndef f():\n    return _il.import_module("orpheus.cp.solver")\n', id="aliased-lazy"),
    pytest.param('import importlib\ndef f(n):\n    return importlib.import_module(f"orpheus.sn.{n}")\n', id="f-string"),
    pytest.param('m = __import__("orpheus.transport.fields._bases", fromlist=["FieldRole"])\n', id="dunder-import"),
    pytest.param('x = __import__("orpheus").__file__\n', id="the-package"),
    pytest.param('import importlib\nimportlib.import_module(f"tests.gates.mc.{name}")\n', id="gate-helper"),
    pytest.param('import importlib\n_MODULE = "orpheus.numerics.enclosure"\ndef f():\n    return importlib.import_module(_MODULE)\n',
                 id="module-constant"),
    pytest.param('import importlib\nMODULE: str = "orpheus.numerics.traced_memo"\nimportlib.import_module(MODULE)\n', id="annotated-constant"),
]

_NOT_STRING_IMPORTS = [
    pytest.param('import importlib\nimportlib.util.find_spec("orpheus.numerics.invariance")\n', id="absence-find-spec"),
    pytest.param('import importlib\nfor name in names:\n    importlib.import_module(name)\n', id="walk-variable"),
    pytest.param('import importlib\nimportlib.import_module(m.name)\n', id="walk-expression"),
    pytest.param('import importlib\nimportlib.import_module(f"{self.name}.{module}")\n', id="synthetic-package"),
    pytest.param('x = __import__("json")\n', id="stdlib"),
    pytest.param('from orpheus.numerics import enclosure\nenclosure.Enclosure\n', id="the-statement"),
    pytest.param('import importlib\n_MODULE = "orpheus.a"\n_MODULE = "orpheus.b"\nimportlib.import_module(_MODULE)\n', id="rebound-name"),
    pytest.param('import importlib\ndef f(_MODULE):\n    return importlib.import_module(_MODULE)\n', id="a-parameter-not-the-constant"),
]


@pytest.mark.foundation
@pytest.mark.parametrize("source", _STRING_IMPORTS)
def test_a_string_import_of_a_project_module_is_found(source: str) -> None:
    """Each planted shape is one finding in a layered file (``numerics``) and none at the gate root."""
    assert len(_string_imports(pathlib.PurePosixPath("numerics/test_probe.py"), source)) == 1, source
    assert _string_imports(pathlib.PurePosixPath("test_probe.py"), source) == []


@pytest.mark.foundation
@pytest.mark.parametrize("source", _NOT_STRING_IMPORTS)
def test_what_stays_dynamic_is_not_a_string_import(source: str) -> None:
    """An absence check, a walk, a synthetic package, the standard library and the statement itself read 0."""
    assert _string_imports(pathlib.PurePosixPath("numerics/test_probe.py"), source) == [], source


#: One kept dynamic import of each kind ruling 27 keeps, in the tree: (file, kind). Each must carry a dynamic
#: call (read by the census below) and read 0 under the gate.
_KEPT = (
    ("numerics/test_axis_identity.py", "walk: pkgutil.walk_packages over orpheus"),
    ("transport/fields/test_role_partners.py", "walk: the transport packages' iter_modules"),
    ("numerics/_traced_memo_synthetic.py", "synthetic package: written under tmp_path at test time"),
    ("numerics/test_invariance_composition.py", "absence: find_spec of the deleted module"),
    ("mesh/test_module_layout.py", "absence: find_spec of the old homes"),
)


def _dynamic_calls(src: str) -> list[ast.Call]:
    names = {"import_module", "__import__", "find_spec"}
    return [
        node for node in ast.walk(ast.parse(src))
        if isinstance(node, ast.Call) and node.args and (
            (isinstance(node.func, ast.Attribute) and node.func.attr in names)
            or (isinstance(node.func, ast.Name) and node.func.id in names))
    ]


@pytest.mark.foundation
@pytest.mark.parametrize(("file", "kind"), _KEPT, ids=[k.split(":")[0] + "-" + f.split("/")[-1] for f, k in _KEPT])
def test_a_kept_dynamic_import_stays_green(file: str, kind: str) -> None:
    src = (GATES_ROOT / file).read_text(encoding="utf-8")
    assert _dynamic_calls(src), f"{file} no longer holds the dynamic call it is named for ({kind})"
    assert _string_imports(pathlib.PurePosixPath(file), src) == []


@pytest.mark.foundation
def test_no_layered_gate_file_imports_a_project_module_by_string() -> None:
    files = _gate_files()
    rels = [pathlib.PurePosixPath(p.relative_to(GATES_ROOT).as_posix()) for p in files]
    judged = [rel for rel in rels if _layer_of(rel) is not None]
    assert len(judged) > 400, len(judged)
    found = [hit for path, rel in zip(files, rels) for hit in _string_imports(rel, path.read_text(encoding="utf-8"))]
    assert not found, f"{len(found)} string import(s) of a project module; {_RULING_27}: " + "; ".join(found)


# ---------------------------------------------------------------------------
# rests_on resolution
# ---------------------------------------------------------------------------

# Targets that resolve to no collected test today, exact (empty since the
# 2026-10-25 repairs: two re-pointed renames, and the two targets of the retired
# geometry/test_reemission_closure.py re-pointed to the successors its own
# retirement note named). The machinery stays: a new dangling target is red.
OPEN_UNRESOLVED_RESTS_ON: frozenset[str] = frozenset()


def _resolves(target: str, collected: set[str], prefixes: set[str]) -> bool:
    """A target with a ``[...]`` case id names one case and must be collected exactly; a target without one
    names a test, a class or a parametrised function, and resolves when it is a collected id or a prefix
    (``target::`` or ``target[``) of one."""
    if "[" in target:
        return target in collected
    return target in collected or target in prefixes


def _prefixes(collected: Iterable[str]) -> set[str]:
    """Every ``::``-prefix of every collected id, case ids stripped."""
    out: set[str] = set()
    for nodeid in collected:
        parts = nodeid.split("[", 1)[0].split("::")
        out.update("::".join(parts[:i]) for i in range(1, len(parts) + 1))
    return out


def _unresolved(manifest: list[Mapping[str, object]]) -> tuple[set[str], set[str]]:
    """(the cited targets, the cited targets that resolve to nothing)."""
    collected = {str(row["nodeid"]) for row in manifest}
    prefixes = _prefixes(collected)
    cited = {str(t) for row in manifest for t in row["rests_on"]}  # type: ignore[union-attr]
    return cited, {t for t in cited if not _resolves(t, collected, prefixes)}


@pytest.mark.foundation
def test_the_resolver_resolves_a_known_id_and_refuses_a_fabricated_one() -> None:
    manifest = [
        {"nodeid": "tests/gates/a/test_x.py::TestC::test_f[case-1]", "rests_on": []},
        {"nodeid": "tests/gates/a/test_x.py::test_g", "rests_on": [
            "tests/gates/a/test_x.py::TestC",  # a class
            "tests/gates/a/test_x.py::TestC::test_f",  # a parametrised function
            "tests/gates/a/test_x.py::TestC::test_f[case-1]",  # one case, exact
            "tests/gates/a/test_x.py::TestC::test_f[case-2]",  # a stale case id
            "tests/gates/nowhere.py::nope",  # fabricated
            "tests/gates/a/test_x.py::test_no_such_test",  # fabricated, real file
        ]},
    ]
    cited, unresolved = _unresolved(manifest)
    assert len(cited) == 6
    assert unresolved == {
        "tests/gates/a/test_x.py::TestC::test_f[case-2]",
        "tests/gates/nowhere.py::nope",
        "tests/gates/a/test_x.py::test_no_such_test",
    }


@pytest.mark.foundation
def test_every_rests_on_target_resolves_but_the_listed_rows(tmp_path: pathlib.Path) -> None:
    """Collect the whole suite in a fresh interpreter (the in-process session holds only what THIS run
    selected) and resolve every ``rests_on`` target against the collected ids."""
    manifest_path = tmp_path / "rests_on_manifest.json"
    result = subprocess.run(
        [sys.executable, "-O", "-m", "pytest", "--collect-only", "-q", "-p", "no:cacheprovider",
         "-p", "tests._harness.rests_on_manifest", "tests"],
        capture_output=True, text=True, cwd=REPO_ROOT,
        env={**os.environ, MANIFEST_ENV: str(manifest_path)},
    )
    assert result.returncode == 0, result.stdout[-4000:] + result.stderr[-4000:]
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    cited, unresolved = _unresolved(manifest)
    # The denominators: the collection and the cited targets (2026-10-25: 19 859 items, 526 targets).
    assert len(manifest) > 15000, len(manifest)
    assert len(cited) > 400, len(cited)
    new = sorted(unresolved - OPEN_UNRESOLVED_RESTS_ON)
    stale = sorted(OPEN_UNRESOLVED_RESTS_ON - unresolved)
    assert not new, f"{len(new)} rests_on target(s) name no collected test: " + "; ".join(new)
    assert not stale, (
        f"{len(stale)} listed target(s) now resolve or are no longer cited; retire them from "
        "OPEN_UNRESOLVED_RESTS_ON: " + "; ".join(stale)
    )
