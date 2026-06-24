r"""Shared unit registry + the field dimensional-class constants (View-G, B.4).

This module is the **single source of truth for units** in ORPHEUS. Under
the "View-G" decision (issues #205 / #207) a
:class:`~orpheus.numerics.space.FunctionSpace` is pure geometry and carries
NO units; units are a property of the *quantity* — they live on the field
role-leaf as the class constant
:attr:`~orpheus.numerics.field.Field.UNITS` (this module) and, in #208, on
the operator as its unit-gain. The two composition chains check
independently: geometric (space identity, at compose/apply) and dimensional
(units, at operator construction).

Why a units constant at all (it is NOT the arithmetic gate)
===========================================================

The arithmetic gate is **class identity**, not units — see
:meth:`~orpheus.numerics.field.Field._check_partner`. This is structural,
not incidental: the ten role leaves fall into only **four distinct unit
signatures** (below), so units are far coarser than classes
(e.g. ``AngularFlux``, ``BoundaryFlux``, ``BoundarySourceSink`` and
``BoundaryResidual`` are all :data:`ANGULAR_FLUX_UNITS`). "Same units" grants
permission to add in linear algebra; it does NOT grant meaning. ``UNITS`` is
therefore **machine-readable metadata**, not a gate. Its consumers are:

* **diagnostics (now):** :meth:`Field._check_partner` and
  :meth:`Field.__repr__` surface ``UNITS`` so a cross-class error reads
  ``"AngularResidual [1/cm³/s/sr] vs AngularSourceSink [1/cm³/s/sr]: same units,
  different role — use AngularResidual.from_balance"``;
* **named composition (B.5):** each residual leaf's ``from_balance``
  factory sanity-checks that its two same-class operands share its own
  unit signature (``sr``-exact) before forming the residual;
* **operator unit-gain (#208):** the construction-time dimensional check
  composes operator gains against these field units.

The four unit signatures (a 2×2: areal/volumetric × angular/scalar)
===================================================================

============================  =================  =====================================
Constant                      Unit               Leaves
============================  =================  =====================================
:data:`ANGULAR_FLUX_UNITS`    ``1/(cm²·s·sr)``   AngularFlux, BoundaryFlux,
                                                 BoundarySourceSink, BoundaryResidual
:data:`SCALAR_FLUX_UNITS`     ``1/(cm²·s)``      ScalarFlux, HarmonicMomentFlux
:data:`ANGULAR_RATE_UNITS`    ``1/(cm³·s·sr)``   AngularSourceSink, AngularResidual
:data:`SCALAR_RATE_UNITS`     ``1/(cm³·s)``      ScalarSourceSink, ScalarResidual
============================  =================  =====================================

The *areal vs volumetric* axis (``cm⁻²`` vs ``cm⁻³``) is *flux* vs
*rate-density* — physical, not a discretization choice. The *angular vs
scalar* axis (``sr`` present vs absent) is per-solid-angle density vs
angle-integrated. Note the boundary leaves are **all-flux** (``1/(cm²·s·sr)``):
on the trace, the prescribed inflow ``q`` and the affine-BC residual are both
flux-typed (added to / differenced from ``γψ``), so — unlike the bulk —
source and residual do NOT pick up the volumetric ``cm⁻³``. This is the
sharpest "units ≠ meaning" case and the reason the gate must be class.

The coefficient signature (#257 — the multiplier-algebra coefficient)
=====================================================================

One further signature carries the cross-section COEFFICIENT — the
:class:`~orpheus.transport.fields._coefficient_role.CoefficientRole` leaf,
distinct from the ten state leaves above. A coefficient is not a flux or a
rate density; it is the *symbol* of a multiplication operator
(``C = M[σ_t]``, the §5.7 promotion). The #208 operator unit-gain is exactly
the product of a state signature and a coefficient signature:

============================  ===============  =====================================
Constant                      Unit             Leaves
============================  ===============  =====================================
:data:`CROSS_SECTION_UNITS`   ``1/cm``         CrossSectionField (σ_t, σ_a, νΣ_f)
============================  ===============  =====================================

``ANGULAR_FLUX_UNITS`` (``cm⁻²·s⁻¹·sr⁻¹``) × ``CROSS_SECTION_UNITS`` (``cm⁻¹``)
= ``ANGULAR_RATE_UNITS`` (``cm⁻³·s⁻¹·sr⁻¹``): multiplying a flux by a cross
section yields the matching rate density — which is why the collision operator
maps the flux signature to the rate signature. (The fission emission spectrum
χ is a dimensionless coefficient too, but its probability-simplex invariant
lives on the source ``Mixture.chi``, not on a field — see #257 — so it adds no
unit signature here.)

Two conventions baked in here (both deliberate)
===============================================

**1. Energy is eV-free.** A stored flux estimate is always integrated over an
energy *bin* — a multigroup group (deterministic) OR a Monte-Carlo tally bin
(continuous-energy). So it is group-integrated by construction:
``φ_g = ∫_{Eg} φ(E) dE`` has units ``1/(cm²·s)``, never ``1/(cm²·s·eV)``.
"Continuous energy" lives in the **cross-section data and collision kernel**
(σ(E) sampled at continuous E), NOT in the flux field — which is exactly what
lets a CE-MC tally be compared against an MG-deterministic solve at the same
dimensional signature (a V&V win). ``eV`` is *defined* in the registry (pint
built-in) so a future per-eV *spectral-density* field type — a derived
``φ(E)`` with ``1/(...·eV)`` — can exist; the current ten leaves carry
exponent-zero on it. Neutron count is dimensionless.

**2. Steradian is kept EXPLICIT and compared EXACTLY.** pint (correctly, per
SI) treats ``sr`` as *dimensionless* — so ``1/(cm²·s·sr)`` and ``1/(cm²·s)``
have identical ``.dimensionality`` and the four signatures collapse to two
SI dimensionality classes. That would make a missing-angular-integration /
missing-``/4π`` bug (the ERR-039 normalization class) dimensionally
invisible. We therefore keep ``sr`` in the unit string and compare units by
**exact equality** (``unit_a == unit_b``), NOT by ``.dimensionality`` — exact
comparison is ``sr``-sensitive (``1/(cm²·s·sr) != 1/(cm²·s)``) and so catches
the normalization-bug class. (Promoting ``sr`` to a base *dimension* via
``UnitRegistry.define`` was evaluated and rejected — it produced an
inconsistent registry in pint 0.25.3.) Consequence for #208: the operator
unit-gain check must likewise compare ``sr``-sensitively (exact unit), not by
SI dimensionality.

Single-registry rule
====================

pint forbids mixing quantities/units across registries, so **every** ``UNITS``
constant MUST come from this one :data:`UREG`. Construction happens once at
import; ``UNITS`` is never touched in the dunder hot path, so ``python -O``
per-op cost is zero (units appear only at construction, in ``-O``-excluded
asserts, and on the error/repr diagnostic paths).

References
----------

* ``.claude/plans/field_role_typing_view_g.md`` — Phase B step B.4.
* ``orpheus/numerics/field.py`` — the View-G enforcement story; ``UNITS`` is
  the field-side half (the operator-side unit-gain is #208).
* Pre-A.1 ``FunctionSpace.units`` (commit ``1e0bb98``, removed ``20dc7d3``) —
  the prior pint-based representation this re-homes onto the field leaf.
"""

from __future__ import annotations

import pint
from pint import Unit

__all__ = [
    "UREG",
    "Unit",
    "ANGULAR_FLUX_UNITS",
    "SCALAR_FLUX_UNITS",
    "ANGULAR_RATE_UNITS",
    "SCALAR_RATE_UNITS",
    "CROSS_SECTION_UNITS",
]

#: The ONE shared unit registry for ORPHEUS. Every ``UNITS`` constant is built
#: from this instance — pint raises on cross-registry mixing.
UREG = pint.UnitRegistry()

# The four field dimensional-class signatures (built once, at import).
# Operator precedence: ``**`` binds before ``/``, so each reads left-to-right
# as the displayed unit.

#: ``1/(cm²·s·sr)`` — areal per-solid-angle flux density. Angular flux and
#: every boundary-trace leaf (flux, source, residual are all flux-typed on
#: the trace).
ANGULAR_FLUX_UNITS: Unit = UREG.cm**-2 / UREG.s / UREG.sr

#: ``1/(cm²·s)`` — areal angle-integrated flux. Scalar flux and the harmonic
#: moments (a moment is angle-integrated; the ``sr`` of ``dΩ`` cancels, and
#: the ``ℓ=0`` moment IS the scalar flux — verified against the
#: ``SphericalHarmonicSpace`` no-prefactor convention).
SCALAR_FLUX_UNITS: Unit = UREG.cm**-2 / UREG.s

#: ``1/(cm³·s·sr)`` — volumetric per-solid-angle rate density. Per-ordinate
#: sources and residuals (operator outputs / external drives).
ANGULAR_RATE_UNITS: Unit = UREG.cm**-3 / UREG.s / UREG.sr

#: ``1/(cm³·s)`` — volumetric angle-integrated rate density. Scalar sources
#: and residuals.
SCALAR_RATE_UNITS: Unit = UREG.cm**-3 / UREG.s

# The coefficient signature (#257 — the multiplier-algebra coefficient).
# NOT a state quantity (flux / rate density) but the cross-section COEFFICIENT
# that promotes to a multiplication operator (C = M[σ_t], the §5.7 promotion).

#: ``1/cm`` — macroscopic cross section :math:`\Sigma` (σ_t, σ_a, νΣ_f, the
#: σ_s diagonal). The coefficient of a collision / production
#: :class:`~orpheus.numerics.operator.LinearOperator` promoted from a
#: :class:`~orpheus.transport.fields.cross_section_field.CrossSectionField`.
#: Multiplied into a flux it yields the matching rate density — the #208
#: operator unit-gain ``ANGULAR_FLUX_UNITS × CROSS_SECTION_UNITS =
#: ANGULAR_RATE_UNITS``.
CROSS_SECTION_UNITS: Unit = UREG.cm**-1
