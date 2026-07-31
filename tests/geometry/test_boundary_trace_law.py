r"""Foundation tests for the :class:`BoundaryTraceLaw` ABC.

This file pins the Issue-#186 (B3 + β2) descriptor-model contract:

* :class:`BoundaryTraceLaw` is a **pure descriptor** — no
  ``apply`` method, no :class:`LinearOperator` inheritance.
* Concrete subclasses MAY add their own ``apply`` for unit-test
  convenience, but the canonical realisation path is via
  :class:`SNBoundaryRealizer.realize`.
* The five ``assert_*`` invariants default to no-ops.
* :meth:`source` returns a :class:`NoSource` sentinel by default.
* :meth:`realize` raises :class:`NotImplementedError` directing
  callers to a method-specific realizer.
* The :class:`BoundaryTraceLaw` registry self-populates via
  ``__init_subclass__(key=...)``.
* :attr:`geometry_map` and :attr:`response_kernel` default to
  ``None``.
* Descriptor-tree algebra (``+``, ``-``, ``*``, ``/``, ``-``)
  returns :class:`LawSum` / :class:`LawScaled` nodes (Issue #186).
  Composition assertions are pinned in
  :mod:`tests.geometry.test_law_composition`.

Wave 7 merged the legacy ``BoundaryOperator`` ABC into
:class:`BoundaryTraceLaw`; Wave O step O.4a.1-β retired the deprecated
alias entirely. The pre-Wave-7 registry-disjointness tests were
replaced with a UNIFIED-registry assertion — all rank-1 BCs
(``vacuum`` / ``reflective`` / ``white`` / ``periodic`` / ``albedo`` /
``prescribed_inflow``) now live in :pyattr:`BoundaryTraceLaw.registry`.

Tagged ``@pytest.mark.foundation`` per :mod:`tests._harness`.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from orpheus.geometry.boundary import (
    InflowSourceSpec,
    BoundaryTraceLaw,
    ConstantInflowSource,
    NoSource,
)


# ─────────────────────────────────────────────────────────────────────
# Stub concrete law used across tests
# ─────────────────────────────────────────────────────────────────────


class _StubLaw(BoundaryTraceLaw, key="_stub_for_test"):
    """Minimal concrete law.

    Used to exercise the registry, the default property surface
    (``geometry_map``, ``response_kernel``, ``source``), and the
    descriptor-tree algebra. The descriptor itself is not callable —
    realisation happens through a method-specific realizer.
    """


# ─────────────────────────────────────────────────────────────────────
# Descriptor surface (Issue #186 / B3 + β2)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_descriptor_has_no_apply_method() -> None:
    """The ABC ships no ``apply`` and no :class:`LinearOperator`.

    Issue #186 / B3 + β2: descriptors are not callable operators.
    Realisation is the sole bridge — route through
    ``SNBoundaryRealizer().realize(law, method_space)`` or
    ``realize_recursively`` for rank-N composition.
    """
    assert not hasattr(BoundaryTraceLaw, "apply")
    # _StubLaw inherits — no apply method present.
    assert not hasattr(_StubLaw(), "apply")


@pytest.mark.foundation
def test_minimal_concrete_subclass_constructs() -> None:
    """A trivial subclass under the descriptor model constructs."""
    law = _StubLaw()
    assert isinstance(law, BoundaryTraceLaw)


# ─────────────────────────────────────────────────────────────────────
# Default-property contract
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_geometry_map_default_is_none() -> None:
    assert _StubLaw().geometry_map is None


@pytest.mark.foundation
def test_response_kernel_default_is_none() -> None:
    assert _StubLaw().response_kernel is None


@pytest.mark.foundation
def test_source_default_is_no_source() -> None:
    law = _StubLaw()
    src = law.source
    assert isinstance(src, NoSource)
    out = src.evaluate((8, 1))
    assert out.shape == (8, 1)
    assert np.all(out == 0.0)


@pytest.mark.foundation
def test_descriptor_has_no_capabilities_attribute() -> None:
    """Capabilities are an operator-tree concept; descriptors don't
    advertise them. Realised operators (from
    :class:`SNBoundaryRealizer`) carry the capability set."""
    assert not hasattr(_StubLaw, "capabilities")
    assert not hasattr(BoundaryTraceLaw, "capabilities")


@pytest.mark.foundation
def test_descriptor_has_no_domain_or_codomain() -> None:
    """``domain`` / ``codomain`` are :class:`LinearOperator`
    attributes; descriptors drop the mixin in B3 + β2."""
    law = _StubLaw()
    assert not hasattr(law, "domain")
    assert not hasattr(law, "codomain")


# ─────────────────────────────────────────────────────────────────────
# assert_* invariants default to no-ops
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_assert_invariants_default_behaviour() -> None:
    """Four ``assert_*`` methods on the ABC default to no-ops; the
    fifth carries the real ERR-047 body (#52).

    Concrete BCs override the relevant invariants with logic that
    raises the typed errors from
    :mod:`orpheus.geometry.boundary._errors`. The first four defaults
    must not raise on ANY input;
    ``assert_source_lives_on_incoming_trace`` probes ``self.source``
    against the quadrature's ordinate count — for the inherited
    :class:`NoSource` (:math:`q \\equiv 0`) it certifies trivially,
    with or without an inflow set. Its no-op era ended at #52 (the
    d3 audit's "a marker today would be blind" finding)."""
    law = _StubLaw()
    # ``None`` as the quadrature stand-in -- the four no-op defaults
    # don't access it.
    assert law.assert_inflow_outflow_classification(None) is None  # type: ignore[arg-type]
    assert law.assert_outgoing_leakage_unconstrained(None) is None  # type: ignore[arg-type]
    assert law.assert_geometry_map_measure_preserving(None) is None  # type: ignore[arg-type]
    assert law.assert_response_positive_if_declared() is None

    class _QuadStub:
        N = 6

    # The real body reads quadrature.N + self.source only; NoSource
    # certifies masklessly (both index states).
    assert law.assert_source_lives_on_incoming_trace(_QuadStub()) is None  # type: ignore[arg-type]
    assert (
        law.assert_source_lives_on_incoming_trace(
            _QuadStub(), np.arange(3),  # type: ignore[arg-type]
        )
        is None
    )


# ─────────────────────────────────────────────────────────────────────
# realize hook
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_realize_default_raises_not_implemented() -> None:
    """The base ``realize`` default raises :class:`NotImplementedError`,
    directing callers to a method-specific realizer
    (e.g. :class:`SNBoundaryRealizer.realize`) or
    :func:`realize_recursively` for descriptor trees."""
    law = _StubLaw()
    with pytest.raises(NotImplementedError, match="method-specific realizer"):
        law.realize(method_space=None)


# ─────────────────────────────────────────────────────────────────────
# Registry contracts
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_stub_law_self_registers_under_its_key() -> None:
    """The ``key="_stub_for_test"`` kwarg on the class statement
    populates :attr:`BoundaryTraceLaw.registry`."""
    assert BoundaryTraceLaw.registry.get("_stub_for_test") is _StubLaw
    # And we can construct via the registry's ``create`` classmethod.
    inst = BoundaryTraceLaw.create("_stub_for_test")
    assert isinstance(inst, _StubLaw)


@pytest.mark.foundation
def test_unified_registry_holds_all_concretes() -> None:
    """Post-Wave-7, the ONE registry on :class:`BoundaryTraceLaw`
    holds every concrete BC.

    Before Wave 7 the legacy ``BoundaryOperator`` ABC and
    :class:`BoundaryTraceLaw` held two disjoint registry dicts. Wave
    7 merged the ABCs into :class:`BoundaryTraceLaw` (Wave O O.4a.1-β
    retired the deprecated ``BoundaryOperator`` alias), and the
    registry merged with them — one dict containing the rank-1
    concretes plus the test stub plus any Wave-7 additions
    (``prescribed_inflow``). Wave 11 removed ``MixedBoundaryOperator``
    and its ``"mixed"`` registry key — rank-N compositions are now
    Wave-0 ``OperatorSum``-algebra over realised leaves, not a
    registered concrete BC.
    """
    # Every Wave-4 rank-1 concrete BC lives in the unified registry.
    keys = set(BoundaryTraceLaw.registry.keys())
    expected_concretes = {
        "vacuum",
        "reflective",
        "white",
        "periodic",
        "albedo",
    }
    assert expected_concretes <= keys, (
        f"missing concretes in unified registry: "
        f"{expected_concretes - keys}"
    )
    # Wave 11 — ``"mixed"`` MUST no longer be a registered key.
    assert "mixed" not in keys, (
        "Wave 11 removed MixedBoundaryOperator; the 'mixed' key MUST "
        "no longer be in BoundaryTraceLaw.registry."
    )
    # The test stub registered in this module is also present.
    assert "_stub_for_test" in keys


# ─────────────────────────────────────────────────────────────────────
# Source helpers
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_no_source_evaluate_returns_zeros_of_requested_shape() -> None:
    src = NoSource()
    out = src.evaluate((12, 2))
    assert out.shape == (12, 2)
    assert out.dtype == np.float64
    assert np.all(out == 0.0)


@pytest.mark.foundation
def test_constant_inflow_source_returns_value_everywhere() -> None:
    src = ConstantInflowSource(value=2.5)
    out = src.evaluate((6, 3))
    assert out.shape == (6, 3)
    assert out.dtype == np.float64
    assert np.all(out == 2.5)


@pytest.mark.foundation
def test_boundary_source_protocol_is_runtime_checkable() -> None:
    """A bare instance of :class:`NoSource` or
    :class:`ConstantInflowSource` satisfies the Protocol via
    duck-typing (the ``evaluate`` method is present)."""
    assert isinstance(NoSource(), InflowSourceSpec)
    assert isinstance(ConstantInflowSource(1.0), InflowSourceSpec)
