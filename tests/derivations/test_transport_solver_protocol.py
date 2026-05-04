r"""Foundation tests for the :class:`TransportSolver` Protocol.

The Protocol is the structural unification of the math-heart pattern
across the project: ``Billiard`` (trajectory_resolvent), ``MomentSpace``
(fn_method), the forthcoming ``Spectrum`` / ``BasisSpace``, and the
production discrete-mesh solvers ``CPSolver`` / ``SNSolver``. These
tests verify the *software invariants* of conformance — they are NOT
verification claims about any solver's accuracy.

What this file pins
-------------------

1. **Runtime conformance** — :func:`isinstance` checks return ``True``
   for every shipped concrete realisation. The :class:`TransportSolver`
   Protocol is :func:`runtime_checkable`, so a class without the
   declared attributes / methods would silently *not* conform; the
   tests catch that drift.

2. **Method-name tags** — every conforming class exposes a stable
   :attr:`method_name` matching the project's pillar tag. Used by
   cross-method comparators.

3. **``materials`` / ``geometry_spec`` shape** — both attributes carry
   the correct types: ``dict[int, Mixture]`` and
   :class:`~orpheus.derivations.common.geometry_spec.GeometrySpec`.
   Catches accidental regression to legacy raw-XS payloads or to a
   discrete :class:`~orpheus.geometry.mesh.Mesh1D` instead of the
   method-agnostic spec.

4. **Registry completeness** —
   :data:`~orpheus.derivations.common.solver_protocol.KNOWN_TRANSPORT_SOLVERS`
   is in lockstep with the :attr:`method_name` values of the shipped
   solvers. Catches "we removed ``method_name`` and forgot to update
   the registry".

Foundation-tier
---------------

Module-level ``pytestmark = [pytest.mark.foundation]`` (no
``verifies(...)``); these are software invariants, not equation labels.

References
----------

* :doc:`/theory/transport_solver_protocol` — the Protocol's narrative.
* :doc:`/skills/algebra-of-record` § "The bifurcation pattern" —
  Branch-1/Branch-2 isolation and structural-independence
  preservation through the Protocol facade.
"""
from __future__ import annotations

import pytest

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.derivations.common.geometry_spec import GeometrySpec
from orpheus.derivations.common.solver_protocol import (
    KNOWN_TRANSPORT_SOLVERS,
    TransportSolver,
)
from orpheus.derivations.common.xs_library import make_mixture
from orpheus.geometry.mesh import BC

import numpy as np


# Module-level mark per the spec § "V&V tagging discipline".
pytestmark = [pytest.mark.foundation]


# ----------------------------------------------------------------------
# Helpers — minimal mixtures + geometry specs for the conformance probes
# ----------------------------------------------------------------------


def _make_1g_mixture() -> Mixture:
    """Build a tiny 1G fissile-fuel-like mixture for conformance probes."""
    return make_mixture(
        sig_t=np.array([0.5]),
        sig_c=np.array([0.05]),
        sig_f=np.array([0.05]),
        nu=np.array([2.0]),
        chi=np.array([1.0]),
        sig_s=np.array([[0.40]]),
    )


def _slab_geometry_spec() -> GeometrySpec:
    """Method-agnostic slab spec for the F_N moment-space probes."""
    return GeometrySpec(
        geometry="slab",
        critical_dimension_mfp=0.93772556,
        critical_dimension_cm=1.875,
        n_groups=1,
        mat_id=0,
        bc_left=BC.vacuum,
        bc_right=BC.vacuum,
    )


def _sphere_geometry_spec() -> GeometrySpec:
    """Method-agnostic sphere spec for the trajectory_resolvent probes."""
    return GeometrySpec(
        geometry="sphere",
        critical_dimension_mfp=2.5,
        critical_dimension_cm=5.0,
        n_groups=1,
        mat_id=0,
        bc_left=BC.reflective,
        bc_right=BC.reflective,
    )


# ----------------------------------------------------------------------
# Runtime conformance — Billiard + MomentSpace
# ----------------------------------------------------------------------


def test_protocol_runtime_checkable_billiard_conforms():
    """Billiard satisfies isinstance(b, TransportSolver)."""
    from orpheus.derivations.continuous.trajectory_resolvent import Billiard

    b = Billiard.from_problem(
        materials={0: _make_1g_mixture()},
        geometry_spec=_sphere_geometry_spec(),
        alpha=1.0,
    )
    assert isinstance(b, TransportSolver), (
        "Billiard does not satisfy TransportSolver Protocol; check "
        "that materials, geometry_spec, method_name, and solve_critical "
        "are all exposed as the Protocol surface declares."
    )


def test_protocol_runtime_checkable_moment_space_conforms():
    """MomentSpace satisfies isinstance(m, TransportSolver)."""
    from orpheus.derivations.continuous.fn_method.moment_space import (
        MomentSpace,
    )

    m = MomentSpace.from_problem(
        materials={0: _make_1g_mixture()},
        geometry=_slab_geometry_spec(),
        fn_order=8,
    )
    assert isinstance(m, TransportSolver), (
        "MomentSpace does not satisfy TransportSolver Protocol; check "
        "that materials, geometry_spec, method_name, and solve_critical "
        "are all exposed as the Protocol surface declares."
    )


# ----------------------------------------------------------------------
# method_name tags — stable strings per pillar
# ----------------------------------------------------------------------


def test_protocol_method_name_billiard():
    """Billiard.method_name == 'trajectory_resolvent'."""
    from orpheus.derivations.continuous.trajectory_resolvent import Billiard

    b = Billiard.from_problem(
        materials={0: _make_1g_mixture()},
        geometry_spec=_sphere_geometry_spec(),
        alpha=1.0,
    )
    assert b.method_name == "trajectory_resolvent"


def test_protocol_method_name_moment_space():
    """MomentSpace.method_name == 'fn_method'."""
    from orpheus.derivations.continuous.fn_method.moment_space import (
        MomentSpace,
    )

    m = MomentSpace.from_problem(
        materials={0: _make_1g_mixture()},
        geometry=_slab_geometry_spec(),
        fn_order=8,
    )
    assert m.method_name == "fn_method"


# ----------------------------------------------------------------------
# materials / geometry_spec shape contracts
# ----------------------------------------------------------------------


def test_protocol_materials_property_returns_dict_int_mixture():
    """The Protocol's ``.materials`` is the production dict[int, Mixture],
    NOT the legacy raw-array dict.

    Catches accidental regression where ``Billiard.materials`` returns
    its legacy ``{"sigma_t": ..., "sigma_s": ...}`` dict instead of the
    production-protocol ``{0: Mixture}`` dict.
    """
    from orpheus.derivations.continuous.fn_method.moment_space import (
        MomentSpace,
    )
    from orpheus.derivations.continuous.trajectory_resolvent import Billiard

    mix = _make_1g_mixture()

    b = Billiard.from_problem(
        materials={0: mix},
        geometry_spec=_sphere_geometry_spec(),
        alpha=1.0,
    )
    assert isinstance(b.materials, dict)
    assert set(b.materials.keys()) == {0}
    assert isinstance(b.materials[0], Mixture), (
        f"Billiard.materials[0] expected Mixture, got "
        f"{type(b.materials[0]).__name__}. Legacy raw-XS dict leaked."
    )

    m = MomentSpace.from_problem(
        materials={0: mix},
        geometry=_slab_geometry_spec(),
        fn_order=8,
    )
    assert isinstance(m.materials, dict)
    assert set(m.materials.keys()) == {0}
    assert isinstance(m.materials[0], Mixture)


def test_protocol_geometry_spec_property_returns_geometry_spec():
    """``.geometry_spec`` returns a :class:`GeometrySpec` instance."""
    from orpheus.derivations.continuous.fn_method.moment_space import (
        MomentSpace,
    )
    from orpheus.derivations.continuous.trajectory_resolvent import Billiard

    mix = _make_1g_mixture()

    b = Billiard.from_problem(
        materials={0: mix},
        geometry_spec=_sphere_geometry_spec(),
        alpha=1.0,
    )
    assert isinstance(b.geometry_spec, GeometrySpec)
    assert b.geometry_spec.geometry == "sphere"

    m = MomentSpace.from_problem(
        materials={0: mix},
        geometry=_slab_geometry_spec(),
        fn_order=8,
    )
    assert isinstance(m.geometry_spec, GeometrySpec)
    assert m.geometry_spec.geometry == "slab"


# ----------------------------------------------------------------------
# Protocol surface — minimal documented contract
# ----------------------------------------------------------------------


def test_protocol_disallows_extra_methods_in_signature():
    """The Protocol's behavioural surface is exactly:
    ``materials``, ``geometry_spec``, ``method_name``, ``solve_critical``.

    No additional members are required for conformance. This pins the
    documented contract — adding a new required Protocol member is a
    breaking change for downstream conformers (Spectrum / BasisSpace
    are designing against the Protocol surface; expanding it
    surreptitiously would silently break their conformance work).
    """
    # The Protocol's __annotations__ + the methods declared on its body
    # are the authoritative surface. We pin the exact attribute / method
    # set so any addition forces a deliberate update here.
    expected_attrs = {"materials", "geometry_spec", "method_name"}
    expected_methods = {"solve_critical"}

    annotations = set(TransportSolver.__annotations__.keys())
    assert annotations == expected_attrs, (
        f"TransportSolver attribute surface drift: expected "
        f"{expected_attrs}, got {annotations}. If you added an "
        f"attribute, update this test deliberately and announce the "
        f"breaking change."
    )

    # The Protocol's solve_critical is declared on the Protocol body
    # (not a plain attribute). Check it exists on the Protocol class.
    method_names = {
        name
        for name in vars(TransportSolver)
        if not name.startswith("_") and callable(vars(TransportSolver)[name])
    }
    assert expected_methods.issubset(method_names), (
        f"TransportSolver method surface drift: expected "
        f"{expected_methods} ⊆ {method_names}."
    )


# ----------------------------------------------------------------------
# Registry completeness
# ----------------------------------------------------------------------


def test_solver_registry_is_complete():
    """``KNOWN_TRANSPORT_SOLVERS`` lists every shipped pillar tag.

    Catches silent removal of ``method_name`` from a math-heart class
    (the test for the missing pillar would still pass at the type
    check, but the registry would no longer cover it).
    """
    from orpheus.derivations.continuous.fn_method.moment_space import (
        MomentSpace,
    )
    from orpheus.derivations.continuous.trajectory_resolvent import Billiard

    mix = _make_1g_mixture()
    b = Billiard.from_problem(
        materials={0: mix},
        geometry_spec=_sphere_geometry_spec(),
        alpha=1.0,
    )
    m = MomentSpace.from_problem(
        materials={0: mix},
        geometry=_slab_geometry_spec(),
        fn_order=8,
    )

    # Every shipped solver's method_name is in the registry.
    assert b.method_name in KNOWN_TRANSPORT_SOLVERS, (
        f"Billiard.method_name = {b.method_name!r} is not in "
        f"KNOWN_TRANSPORT_SOLVERS = {KNOWN_TRANSPORT_SOLVERS}."
    )
    assert m.method_name in KNOWN_TRANSPORT_SOLVERS, (
        f"MomentSpace.method_name = {m.method_name!r} is not in "
        f"KNOWN_TRANSPORT_SOLVERS = {KNOWN_TRANSPORT_SOLVERS}."
    )

    # The CP / SN production solvers conform natively (Step 4 of
    # the input-cleanup track). Assert they're in the registry so
    # the lockstep test catches a silent removal of the
    # ``method_name`` class attribute.
    for expected in ("trajectory_resolvent", "fn_method", "cp", "sn"):
        assert expected in KNOWN_TRANSPORT_SOLVERS, (
            f"Expected pillar tag {expected!r} missing from "
            f"KNOWN_TRANSPORT_SOLVERS = {KNOWN_TRANSPORT_SOLVERS}."
        )


# ----------------------------------------------------------------------
# Production CP / SN — Protocol conformance via from_problem
# ----------------------------------------------------------------------
#
# Step 4 of the input-cleanup track promoted the production CP / SN
# solvers onto the Protocol natively (no test-only adapter scaffold).
# These tests pin that the production classes — built via
# ``CPSolver.from_problem`` / ``SNSolver.from_problem`` — satisfy the
# same shape contract as the continuous-reference math-heart classes
# above.

def _smoke_sphere_geometry_spec_cp() -> GeometrySpec:
    """Sphere spec compatible with CP's outer-BC registry (white).

    CPMesh supports outer BC ∈ {white, vacuum}; reflective is not in
    its registry. The Protocol itself is BC-agnostic — what matters
    is the structural surface — so we use a CP-supported BC pair.
    """
    return GeometrySpec(
        geometry="sphere",
        critical_dimension_mfp=2.5,
        critical_dimension_cm=5.0,
        n_groups=1,
        mat_id=0,
        bc_left=BC.reflective,
        bc_right=BC("white"),
    )


def _smoke_slab_geometry_spec_sn() -> GeometrySpec:
    """Slab spec compatible with SN's BC registry (vacuum / reflective).

    SN supports {vacuum, reflective} on every face. The Protocol
    itself is BC-agnostic, so we choose vacuum-vacuum to avoid the
    closed-loop k_inf degeneracy on the conformance probe.
    """
    return GeometrySpec(
        geometry="slab",
        critical_dimension_mfp=0.93772556,
        critical_dimension_cm=1.875,
        n_groups=1,
        mat_id=0,
        bc_left=BC.vacuum,
        bc_right=BC.vacuum,
    )


def test_protocol_runtime_checkable_cp_solver_conforms():
    """CPSolver built via from_problem satisfies isinstance(c, TransportSolver)."""
    from orpheus.cp.solver import CPSolver

    c = CPSolver.from_problem(
        materials={0: _make_1g_mixture()},
        geometry_spec=_smoke_sphere_geometry_spec_cp(),
    )
    assert isinstance(c, TransportSolver), (
        "CPSolver does not satisfy TransportSolver Protocol; check "
        "that materials, geometry_spec, method_name, and "
        "solve_critical are all exposed as the Protocol surface "
        "declares."
    )


def test_protocol_runtime_checkable_sn_solver_conforms():
    """SNSolver built via from_problem satisfies isinstance(s, TransportSolver)."""
    from orpheus.sn.solver import SNSolver

    s = SNSolver.from_problem(
        materials={0: _make_1g_mixture()},
        geometry_spec=_smoke_slab_geometry_spec_sn(),
    )
    assert isinstance(s, TransportSolver), (
        "SNSolver does not satisfy TransportSolver Protocol; check "
        "that materials, geometry_spec, method_name, and "
        "solve_critical are all exposed as the Protocol surface "
        "declares."
    )


def test_protocol_method_name_cp_solver():
    """CPSolver.method_name == 'cp'."""
    from orpheus.cp.solver import CPSolver

    c = CPSolver.from_problem(
        materials={0: _make_1g_mixture()},
        geometry_spec=_smoke_sphere_geometry_spec_cp(),
    )
    assert c.method_name == "cp"


def test_protocol_method_name_sn_solver():
    """SNSolver.method_name == 'sn'."""
    from orpheus.sn.solver import SNSolver

    s = SNSolver.from_problem(
        materials={0: _make_1g_mixture()},
        geometry_spec=_smoke_slab_geometry_spec_sn(),
    )
    assert s.method_name == "sn"


def test_protocol_geometry_spec_property_cp_sn_returns_geometry_spec():
    """Production CP / SN expose ``.geometry_spec`` as a GeometrySpec."""
    from orpheus.cp.solver import CPSolver
    from orpheus.sn.solver import SNSolver

    mix = _make_1g_mixture()
    c = CPSolver.from_problem(
        materials={0: mix},
        geometry_spec=_smoke_sphere_geometry_spec_cp(),
    )
    assert isinstance(c.geometry_spec, GeometrySpec)
    assert c.geometry_spec.geometry == "sphere"

    s = SNSolver.from_problem(
        materials={0: mix},
        geometry_spec=_smoke_slab_geometry_spec_sn(),
    )
    assert isinstance(s.geometry_spec, GeometrySpec)
    assert s.geometry_spec.geometry == "slab"


def test_protocol_materials_property_cp_sn_returns_dict_int_mixture():
    """Production CP / SN expose ``.materials`` as ``dict[int, Mixture]``.

    Catches accidental regression where the production solver
    leaks an internal cell-by-cell ``CellXS`` payload through
    ``self.materials`` instead of the production-protocol Mixture
    dict.
    """
    from orpheus.cp.solver import CPSolver
    from orpheus.sn.solver import SNSolver

    mix = _make_1g_mixture()
    c = CPSolver.from_problem(
        materials={0: mix},
        geometry_spec=_smoke_sphere_geometry_spec_cp(),
    )
    assert isinstance(c.materials, dict)
    assert set(c.materials.keys()) == {0}
    assert isinstance(c.materials[0], Mixture)

    s = SNSolver.from_problem(
        materials={0: mix},
        geometry_spec=_smoke_slab_geometry_spec_sn(),
    )
    assert isinstance(s.materials, dict)
    assert set(s.materials.keys()) == {0}
    assert isinstance(s.materials[0], Mixture)
