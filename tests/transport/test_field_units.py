r"""Foundation tests for the B.4 ``UNITS`` class constant + units diagnostics.

Pins the View-G units-on-field commitment (step B.4):

* every role leaf declares a :class:`pint.Unit` ``UNITS`` class constant
  drawn from the four signatures in :mod:`orpheus.numerics.units`;
* the four signatures are pairwise distinct under **exact** unit
  comparison — but collapse to **two** SI ``dimensionality`` classes
  because ``sr`` is dimensionless (the documented reason we compare units
  exactly, not by dimensionality: so a missing-``sr`` / missing-angular-
  integration bug stays visible);
* ``UNITS`` is metadata, NOT the arithmetic gate — same-role pairs
  (e.g. ``AngularSourceSink`` / ``AngularResidual``, and the all-flux boundary
  trio) share a unit signature yet are distinct classes, so it is **class
  identity** that separates them;
* ``UNITS`` is a ``ClassVar`` (absent from the dataclass fields), so it is
  never per-instance state and never touched on the arithmetic hot path
  (the ``python -O`` zero-cost structural guarantee);
* the diagnostics surface ``UNITS`` in :meth:`Field.__repr__` and in the
  cross-class :class:`TypeError`.

``foundation`` — software invariants, not physics claims.

References
----------

* ``.claude/plans/field_role_typing_view_g.md`` — Phase B step B.4.
* ``orpheus/numerics/units.py`` — the registry + the four signatures +
  the eV-free / exact-``sr`` convention rationale.
"""
from __future__ import annotations

import numpy as np
import pint
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.field import Field
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.units import (
    ANGULAR_FLUX_UNITS,
    ANGULAR_RATE_UNITS,
    SCALAR_FLUX_UNITS,
    SCALAR_RATE_UNITS,
    UREG,
)
from orpheus.sn.geometry import SNMesh
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.boundary_flux import BoundaryFlux
from orpheus.transport.fields.harmonic_moment_field import HarmonicMomentField
from orpheus.transport.fields.scalar_flux import ScalarFlux
from orpheus.transport.residuals import (
    AngularResidual,
    BoundaryResidual,
    ScalarResidual,
)
from orpheus.transport.source_sinks import AngularSourceSink, BoundarySourceSink, ScalarSourceSink

from tests.sn._test_helpers import placeholder_materials

pytestmark = pytest.mark.foundation

# ── The role grid → its expected dimensional-class signature ─────────
#: (leaf class, expected UNITS constant). The single source of truth the
#: rest of the module checks against.
LEAF_UNITS = [
    (AngularFlux, ANGULAR_FLUX_UNITS),
    (BoundaryFlux, ANGULAR_FLUX_UNITS),
    (BoundarySourceSink, ANGULAR_FLUX_UNITS),
    (BoundaryResidual, ANGULAR_FLUX_UNITS),
    (ScalarFlux, SCALAR_FLUX_UNITS),
    (HarmonicMomentField, SCALAR_FLUX_UNITS),
    (AngularSourceSink, ANGULAR_RATE_UNITS),
    (AngularResidual, ANGULAR_RATE_UNITS),
    (ScalarSourceSink, SCALAR_RATE_UNITS),
    (ScalarResidual, SCALAR_RATE_UNITS),
]

THE_FOUR = [
    ANGULAR_FLUX_UNITS, SCALAR_FLUX_UNITS, ANGULAR_RATE_UNITS, SCALAR_RATE_UNITS,
]


def _slab_mesh(nx: int = 4, ng: int = 2) -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


# ════════════════════════════════════════════════════════════════════
# Every leaf declares the right UNITS, from the one shared registry.
# ════════════════════════════════════════════════════════════════════


class TestUnitsDeclaration:
    @pytest.mark.parametrize(
        "leaf, expected", LEAF_UNITS, ids=lambda x: getattr(x, "__name__", "")
    )
    def test_leaf_units_match_expected(self, leaf, expected) -> None:
        assert leaf.UNITS == expected

    @pytest.mark.parametrize(
        "leaf, _expected", LEAF_UNITS, ids=lambda x: getattr(x, "__name__", "")
    )
    def test_leaf_units_is_pint_unit_from_shared_registry(
        self, leaf, _expected
    ) -> None:
        assert isinstance(leaf.UNITS, pint.Unit)
        # Single-registry rule: every UNITS comes from the one UREG (pint
        # forbids cross-registry mixing).
        assert leaf.UNITS._REGISTRY is UREG

    def test_field_declares_units_contract_but_sets_no_value(self) -> None:
        """``Field`` declares the bare ``UNITS: ClassVar`` contract (the
        units-side companion of ``_SPACE_NAME``) so every concrete leaf's
        obligation to set it is typed and visible — but the abstract base
        assigns no value. Diagnostics fall back to ``"<no units>"``; a
        consumer that reads ``.UNITS`` directly (e.g. #208's unit-gain
        check) gets an ``AttributeError`` on a forgetful leaf rather than a
        silent ``None``."""
        assert "UNITS" in Field.__annotations__          # contract declared
        assert getattr(Field, "UNITS", None) is None      # ... but unset on base


# ════════════════════════════════════════════════════════════════════
# The four signatures: distinct exactly, but only two SI dim classes.
# ════════════════════════════════════════════════════════════════════


class TestFourSignatures:
    def test_four_signatures_pairwise_distinct_exact(self) -> None:
        import itertools
        assert all(a != b for a, b in itertools.combinations(THE_FOUR, 2))

    def test_collapse_to_two_si_dimensionality_classes(self) -> None:
        """``sr`` is dimensionless (SI), so the four unit signatures share
        only TWO ``.dimensionality`` classes (areal vs volumetric). This is
        WHY we compare units exactly, not by dimensionality — exact
        comparison keeps the missing-``sr`` bug class visible."""
        dim_classes = {str(u.dimensionality) for u in THE_FOUR}
        assert len(dim_classes) == 2
        # The angular/scalar pair within each areal/volumetric group is
        # dimensionally identical (sr dimensionless) yet exactly distinct.
        assert ANGULAR_FLUX_UNITS.dimensionality == SCALAR_FLUX_UNITS.dimensionality
        assert ANGULAR_FLUX_UNITS != SCALAR_FLUX_UNITS
        assert ANGULAR_RATE_UNITS.dimensionality == SCALAR_RATE_UNITS.dimensionality
        assert ANGULAR_RATE_UNITS != SCALAR_RATE_UNITS


# ════════════════════════════════════════════════════════════════════
# UNITS is metadata, not the gate: same-units / different-class pairs.
# ════════════════════════════════════════════════════════════════════


class TestSameUnitsDifferentRole:
    @pytest.mark.parametrize(
        "a, b",
        [
            (AngularSourceSink, AngularResidual),   # rate density, Ω
            (ScalarSourceSink, ScalarResidual),     # rate density, scalar
            (AngularFlux, BoundaryFlux),        # flux density, Ω
            (BoundaryFlux, BoundarySourceSink),     # boundary all-flux
            (BoundarySourceSink, BoundaryResidual),
            (ScalarFlux, HarmonicMomentField),  # areal scalar (ℓ=0 IS φ)
        ],
    )
    def test_same_units_but_distinct_class(self, a, b) -> None:
        assert a.UNITS == b.UNITS       # identical dimensional signature
        assert a is not b               # ... but distinct classes
        # => the arithmetic gate must be class identity, not units.


# ════════════════════════════════════════════════════════════════════
# `python -O` zero-cost: UNITS is a ClassVar, not per-instance state.
# ════════════════════════════════════════════════════════════════════


class TestUnitsIsClassVarNotField:
    @pytest.mark.parametrize(
        "leaf, _expected", LEAF_UNITS, ids=lambda x: getattr(x, "__name__", "")
    )
    def test_units_is_classvar_not_a_real_field(self, leaf, _expected) -> None:
        """``UNITS`` is a genuine ``ClassVar`` — absent from the *public*
        dataclass fields and from ``__init__``, with internal
        ``_field_type == _FIELD_CLASSVAR``. So it is never per-instance
        state and never touched in the values-arithmetic dunders (the
        ``-O`` zero-cost guarantee). (Note: ClassVar pseudo-fields DO appear
        in the raw ``__dataclass_fields__`` mapping — the public
        :func:`dataclasses.fields` is what filters them out.)"""
        import dataclasses
        from dataclasses import _FIELD_CLASSVAR

        assert "UNITS" not in {f.name for f in dataclasses.fields(leaf)}
        assert leaf.__dataclass_fields__["UNITS"]._field_type is _FIELD_CLASSVAR

    def test_units_is_class_level_shared(self) -> None:
        """The constant lives on the class and is shared by every instance
        (one object, constructed once at import)."""
        m = _slab_mesh()
        a = AngularSourceSink.from_mesh(np.zeros((m.quad.N, m.ng, m.nx, m.ny)), m)
        b = AngularSourceSink.from_mesh(np.ones((m.quad.N, m.ng, m.nx, m.ny)), m)
        assert a.UNITS is b.UNITS is AngularSourceSink.UNITS


# ════════════════════════════════════════════════════════════════════
# Diagnostics: repr + the cross-class error surface UNITS.
# ════════════════════════════════════════════════════════════════════


class TestUnitsDiagnostics:
    def test_repr_is_concise_and_shows_units(self) -> None:
        m = _slab_mesh()
        src = AngularSourceSink.from_mesh(np.zeros((m.quad.N, m.ng, m.nx, m.ny)), m)
        r = repr(src)
        assert r.startswith("AngularSourceSink(")
        assert "shape=(4, 2, 4, 1)" in r
        assert "units=1/cm³/s/sr" in r
        # No raw ndarray dump (the dataclass auto-repr wart is gone).
        assert "array(" not in r

    def test_cross_class_error_shows_both_units_and_guidance(self) -> None:
        m = _slab_mesh()
        shp = (m.quad.N, m.ng, m.nx, m.ny)
        src = AngularSourceSink.from_mesh(np.ones(shp), m)
        res = AngularResidual.from_mesh(np.ones(shp), m)
        with pytest.raises(TypeError) as ei:
            _ = res - src
        msg = str(ei.value)
        # both operands' (identical) units are surfaced...
        assert msg.count("1/cm³/s/sr") == 2
        # ... with the "same units ≠ meaning" guidance + the named route,
        # and the legacy "same-class" token existing matchers rely on.
        assert "same-class" in msg
        assert "not meaning" in msg
        assert "from_balance" in msg

    def test_error_label_for_non_field_partner(self) -> None:
        """A non-Field operand has no UNITS — the label degrades gracefully."""
        m = _slab_mesh()
        src = AngularSourceSink.from_mesh(np.ones((m.quad.N, m.ng, m.nx, m.ny)), m)
        with pytest.raises(TypeError, match="<no units>"):
            _ = src + 42  # type: ignore[operator]
