r"""#707 Step A1, the member seam on the production builders' composites: each builder wires a zero element, so the
composite space mints its members.

The numerics laws (leaf, dual, tensor product, the kernel basis) are
``tests/gates/numerics/test_member_seam_laws.py``. This file poses the composites that only a builder above the numerics
layer can wire (``DirectSumSpace.carrier.zeros`` in ``sn/problem.py``, ``transport/full_field.py`` and
``diffusion/augmented_mesh.py``; ``DirectSumSpace``'s ``zeros`` in ``sn/coupled_system.py``), on the smallest production
problems: a 3-cell 2-group slab (vacuum | mirror, LS4), the 2-cell reflective sphere (GL2, which carries System B,
the ψ½ ray), and a 4-cell 2-group diffusion slab (vacuum | mirror). It lives in ``cross_method/`` because it imports SN
and diffusion together.

Rows (THEOREM; ``foundation``):

* MB1: ``coordinates_of(member(c)) == c`` bit-exactly; the member is of the carrier class the builder names; its
  block spaces ARE the composite's blocks (``is``); the member's derived space equals V by content (``==``, see below).
* MB2: the other leg. ``x`` is built by the carrier classes' own constructors (``AngularFlux(values=…, space=…)``
  and so on, never ``member``); ``coordinates_of(x)`` is the hand concatenation of the block values (interior, then
  trace; system order for a coupled space), and ``member(coordinates_of(x))`` is ``x`` again in coordinates and in
  every leaf class.
* MB3: every production composite here is full rank, so its kernel basis has zero columns (the block-wise law with
  live columns is the numerics file's KC row).
* MB4 (strict xfail, A3): the dual coupled space's members live on the dual. Today the builder's cotangent zero is
  ``*SourceSink`` blocks on the PRIMAL leaf spaces, so ``member(c).space`` equals the primal coupled space and not
  ``V.dual()``.

Why ``==`` and not ``is`` for a composite: ``Composite.space`` and ``CoupledField.space`` are derived per instance from
the blocks (the F4 ruling, ``transport/full_field.py``: "Not ``is``: the wrapper is minted per composite").
``FunctionSpace.member``'s docstring states the same: the same object on a leaf and its dual, an equal space on a
composite. ``[M]`` 6 of 6 composites here read ``is`` False and ``==`` True.
"""

from __future__ import annotations

from typing import Any, Callable

import numpy as np
import pytest

from orpheus.data.macro_xs.synthetic import get_mixture
from orpheus.diffusion import DiffusionMesh
from orpheus.geometry import StructuredGeometry
from orpheus.geometry.deck import MirrorDeck
from orpheus.mesh import CellsByCount, Mesher
from orpheus.numerics.coupled_system import CoupledField
from orpheus.numerics.direct_sum import DirectSumSpace
from orpheus.sn.coupled_system import build_within_group_system
from orpheus.sn.solver import SNSolver
from orpheus.transport.boundary_conditions import WallCondition
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.scalar_boundary_flux import ScalarBoundaryFlux
from orpheus.transport.fields.scalar_flux import ScalarFlux
from orpheus.transport.full_field import FullField
from orpheus.transport.radial_characteristic_field import (
    RadialCharacteristicBoundaryFlux,
    RadialCharacteristicField,
    RadialCharacteristicInteriorFlux,
)
from tests.gates._exterior import declare_sides
from tests.gates.sn._singular_loss_box import REFLECTIVE, VACUUM, build, scatterer
from tests.gates.sn._test_helpers import make_tiny_spherical_sn_mesh

pytestmark = pytest.mark.foundation

_NM1 = "tests/gates/numerics/test_true_ends_member_seam.py::test_nm1_a_space_mints_its_members"
_M1 = "tests/gates/numerics/test_member_seam_laws.py::test_m1_coordinates_of_member_is_the_identity"
_M2 = "tests/gates/numerics/test_member_seam_laws.py::test_m2_member_of_coordinates_is_the_identity"
_KC = "tests/gates/numerics/test_member_seam_laws.py::test_kc_a_direct_sum_kernel_basis_is_block_wise"


def _slab():
    problem, system, _ = build((3,), ((VACUUM, REFLECTIVE),), scatterer(2))
    return problem, system


def _sphere():
    problem = make_tiny_spherical_sn_mesh()
    system = build_within_group_system(problem, SNSolver(problem, inner_solver="source_iteration").problem.mat_xs)
    return problem, system


def _diffusion() -> DiffusionMesh:
    (left, right), conditions = declare_sides(WallCondition.vacuum(), MirrorDeck())
    geometry = StructuredGeometry.slab((0.0, 10.0), (0,), left=left, right=right)
    mesh = Mesher(geometry).partition(CellsByCount.uniform_width(4)).mesh
    return DiffusionMesh(mesh, {0: get_mixture("A", "2g")}, boundary_conditions=conditions)


#: The leaf classes of each composite carrier, named here (not read off the builder's zero) so MB2's ``x`` is
#: independent of the factory ``member`` reads.
_SN_FULL = (FullField, AngularFlux, AngularBoundaryFlux)
_RC = (RadialCharacteristicField, RadialCharacteristicInteriorFlux, RadialCharacteristicBoundaryFlux)
_DIFF_FULL = (FullField, ScalarFlux, ScalarBoundaryFlux)


def _spaces() -> dict[str, tuple[Callable[[], Any], tuple]]:
    """name -> (builder of the composite space, its carrier shape): a 3-tuple of classes for a ``DirectSumSpace``, a
    tuple of such 3-tuples (system order) for a ``DirectSumSpace``."""
    return {
        "SN slab full_field_space": (lambda: _slab()[0].full_field_space, _SN_FULL),
        "SN sphere full_field_space": (lambda: _sphere()[0].full_field_space, _SN_FULL),
        "SN sphere radial_characteristic_field_space": (lambda: _sphere()[0].radial_characteristic_field_space, _RC),
        "SN slab coupled space (System A)": (lambda: _slab()[1].space, (_SN_FULL,)),
        "SN sphere coupled space (Systems A, B)": (lambda: _sphere()[1].space, (_SN_FULL, _RC)),
        "diffusion full_field_space": (lambda: _diffusion().full_field_space, _DIFF_FULL),
    }


_CASES = list(_spaces())


def _is_coupled(V: Any) -> bool:
    """Whether ``V`` is a coupled sum: a sum whose labels are neither the full field's nor the ray instance's (a hand rule on the production labels)."""
    return V.labels not in (('bulk', 'trace'), ('interior', 'boundary'))


def _check_composite(case: str, V: Any, x: Any, classes: tuple) -> None:
    """``x`` is a member of the DirectSumSpace ``V``: carrier and leaf classes, blocks ``is`` V's, space ``==`` V."""
    outer, interior_cls, boundary_cls = classes
    if (type(x), type(x.interior), type(x.boundary)) != (outer, interior_cls, boundary_cls):
        pytest.fail(
            f"[{case}] member classes {(type(x).__name__, type(x.interior).__name__, type(x.boundary).__name__)}, "
            f"the builder names {(outer.__name__, interior_cls.__name__, boundary_cls.__name__)}"
        )
    if x.interior.space is not V.summand(V.labels[0]) or x.boundary.space is not V.summand(V.labels[1]):
        pytest.fail(f"[{case}] the member's block spaces are not the composite's blocks")
    if not x.space == V:
        pytest.fail(f"[{case}] the member's derived space differs from V by content")


def _check_member(case: str, V: Any, x: Any, classes: tuple) -> None:
    if _is_coupled(V):
        if type(x) is not CoupledField or len(x.systems) != len(classes):
            pytest.fail(f"[{case}] member is a {type(x).__name__} of arity {len(getattr(x, 'systems', ()))}")
        for i, (member_space, member, member_classes) in enumerate(zip((s for _, s in V.summands), x.systems, classes)):
            _check_composite(f"{case} system {i}", member_space, member, member_classes)
        if not x.space == V:
            pytest.fail(f"[{case}] the coupled member's derived space differs from V by content")
    else:
        _check_composite(case, V, x, classes)


def _build_composite(V: Any, classes: tuple, rng: np.random.Generator) -> tuple[Any, np.ndarray]:
    """``x`` from its classes' constructors on V's blocks, and its hand coordinates (interior then trace)."""
    outer, interior_cls, boundary_cls = classes
    bulk = rng.standard_normal(V.summand(V.labels[0]).shape)
    trace = rng.standard_normal(V.summand(V.labels[1]).shape)
    x = outer(
        interior=interior_cls(values=bulk, space=V.summand(V.labels[0])),
        boundary=boundary_cls(values=trace, space=V.summand(V.labels[1])),
    )
    return x, np.concatenate([bulk.ravel(), trace.ravel()])


def _build(V: Any, classes: tuple, rng: np.random.Generator) -> tuple[Any, np.ndarray]:
    if _is_coupled(V):
        built = [_build_composite(s, c, rng) for s, c in zip((s for _, s in V.summands), classes)]
        return CoupledField(systems=tuple(b[0] for b in built), labels=V.labels), np.concatenate([b[1] for b in built])
    return _build_composite(V, classes, rng)


@pytest.mark.rests_on(_NM1, _M1)
@pytest.mark.parametrize("case", _CASES)
def test_mb1_coordinates_of_member_is_the_identity(case: str) -> None:
    """MB1. ``coordinates_of(member(c)) == c`` bit-exactly, the member of the builder's carrier classes on V's blocks.
    First reds: ``DirectSumSpace.member`` splitting ``c`` at the wrong index; ``DirectSumSpace.member`` reversing the
    systems; an unwired ``carrier.zeros`` (the loud refusal)."""
    builder, classes = _spaces()[case]
    V = builder()
    c = np.random.default_rng(3).standard_normal(V.n_coordinates)
    x = V.member(c)
    _check_member(case, V, x, classes)
    np.testing.assert_array_equal(np.asarray(V.coordinates_of(x)).ravel(), c, err_msg=f"[{case}] coordinates_of(member(c)) != c")


@pytest.mark.rests_on(_NM1, _M2)
@pytest.mark.parametrize("case", _CASES)
def test_mb2_member_of_coordinates_is_the_identity(case: str) -> None:
    """MB2. ``x`` built by its carrier classes' constructors: ``coordinates_of(x)`` is the hand concatenation of its
    block values, and ``member(coordinates_of(x))`` is ``x`` (coordinates bit-exact, every class).
    First red: ``coordinates_of`` placing the trace block before the interior."""
    builder, classes = _spaces()[case]
    V = builder()
    x, hand = _build(V, classes, np.random.default_rng(5))
    coords = np.asarray(V.coordinates_of(x)).ravel()
    np.testing.assert_array_equal(coords, hand, err_msg=f"[{case}] coordinates_of(x) is not interior then trace")
    y = V.member(coords)
    _check_member(case, V, y, classes)
    np.testing.assert_array_equal(np.asarray(V.coordinates_of(y)).ravel(), hand, err_msg=f"[{case}] member(coordinates_of(x)) != x")


@pytest.mark.rests_on(_NM1)
@pytest.mark.parametrize("case", _CASES)
def test_mb5_the_dual_refuses_a_primal_composite(case: str) -> None:
    """MB5, the identity leg. A composite member of the ``DirectSumSpace`` V is not a member of ``V.dual()``, an
    equal-shape space with the inverse metric: ``coordinates_of`` compares the composite's derived space with its own,
    as the base does for a ``Member`` (``space.py``, the identity leg), and refuses with the base's message.
    The coupled sums are in the population (``CoupledField`` carries blocks and labels, and the sum's identity leg serves it too).
    First red: the composite override running only the shape leg (the A1 elegance review, finding 3: accepted on 4 of
    4 composites before the leg was restored; with the identity leg removed from ``DirectSumSpace.coordinates_of``, 6 of 6 cases accept)."""
    builder, classes = _spaces()[case]
    V = builder()
    x, _ = _build(V, classes, np.random.default_rng(7))
    with pytest.raises(ValueError, match="is not a member of"):
        V.dual().coordinates_of(x)


@pytest.mark.rests_on(_KC)
@pytest.mark.parametrize("case", _CASES)
def test_mb3_a_full_rank_production_composite_has_no_kernel_columns(case: str) -> None:
    """MB3. The production composites here are full rank (``rank == n``), so ``kernel_basis`` is ``(n, 0)``.
    First red: a block-wise override that returns a block's identity instead of its kernel."""
    V = _spaces()[case][0]()
    K = np.asarray(V.kernel_basis())
    if V.rank != V.n_coordinates or K.shape != (V.n_coordinates, 0):
        pytest.fail(f"[{case}] rank {V.rank} of {V.n_coordinates}, kernel basis {K.shape}")


@pytest.mark.xfail(
    strict=True,
    reason="#707 A3 (the integrated data side): the cotangent zero is *SourceSink on the PRIMAL leaves, so a dual "
    "coupled member's space equals the primal coupled space, not V.dual()",
)
@pytest.mark.rests_on(_NM1)
def test_mb4_a_dual_coupled_member_lives_on_the_dual() -> None:
    """MB4 (strict xfail until A3). ``V.dual().member(c).space == V.dual()`` on the slab's coupled space. ``[M]``
    2026-11-12: the member's blocks are ``AngularSourceSink`` / ``AngularBoundarySourceSink`` on the primal
    ``angular_trial_space`` / ``angular_trace``; ``== V`` True, ``== V.dual()`` False."""
    V = _slab()[1].space
    D = V.dual()
    y = D.member(np.arange(D.n_coordinates, dtype=float))
    if not y.space == D:
        pytest.fail(f"the dual coupled member's space is not V.dual() (equals the primal: {y.space == V})")
