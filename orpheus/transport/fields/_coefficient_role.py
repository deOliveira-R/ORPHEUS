r"""The :class:`CoefficientRole` mixin — the coefficient algebra role marker.

A cross section is not a *state* — it is a **coefficient**, the symbol of a
zeroth-order operator. The grand report (§5.5–5.7) names cross-section fields
``CoefficientField``\ s and the operators they become *promotions*
(``C = M[σ_t]``, the multiplier-algebra embedding ``M: L^∞ → B(L²)``).
:class:`CoefficientRole` is the role marker on the field side of that
promotion, the **exact complement** of
:class:`~orpheus.transport.fields._flux_role.FluxRole`:

==========================  ====================================  ===============================
Role                        Algebra                               Mints
==========================  ====================================  ===============================
``FluxRole`` (state)        affine torsor — NO origin             ``flux ⊖ flux → displacement``
                            (``flux + flux`` forbidden)
``CoefficientRole`` (coef)  plain vector space — HAS an origin    nothing new (keeps ``Field``)
                            (``σ + σ`` legitimate)
==========================  ====================================  ===============================

The flux leaves form an affine space (no natural origin — the iterate has no
"zero flux" baseline that arithmetic respects), so :class:`FluxRole`
*overrides* the additive algebra (removes ``flux + flux``, retypes
``flux − flux`` to a displacement). The coefficients are the opposite: they
form a genuine **vector space with an origin** (the zero cross section
``Σ = 0`` IS a coefficient — it promotes to ``M_0 = ZeroOperator``), closed
under ``+`` (homogenisation ``Σ_mix = Σ_m N_m Σ_m`` is a number-density-weighted
sum), unary ``−``, and scalar ``·`` — i.e. exactly the plain
:class:`~orpheus.numerics.field.Field` vector-space dunders. So
:class:`CoefficientRole` overrides **nothing**: the *absence* of an affine gate
IS its content relative to :class:`FluxRole`. ``σ − σ′`` stays the same class
(a coefficient difference is itself a coefficient — possibly signed; physical
nonnegativity is the *cone*, a tested property of physical inputs, not a type
invariant on every intermediate — keeping the multiplier-algebra domain a full
vector space).

Why a marker mixin and not a bare ``Field`` subclass
----------------------------------------------------

Two reasons it earns its name even while dunder-empty:

* **Taxonomy.** It marks the coefficient leaves —
  :class:`~orpheus.transport.fields.cross_section_field.CrossSectionField`
  today, and the inverse-velocity ``1/v`` time-mass coefficient (units
  ``s/cm``, the ``TimeMassOperator``'s multiplier) imminent next — as one
  family. ``isinstance(x, CoefficientRole)`` discriminates a coefficient from
  a flux state without a unit check, mirroring how ``FluxRole`` tags the state
  leaves.
* **The future multiplier product.** When ``MultiplicationOperator`` lands
  (#257 S3), the law ``M_f @ M_g = M_{f·g}`` needs the pointwise *field·field*
  product ``f·g`` — the commutative-algebra multiplication that the flux never
  has. That product is dimensionally graded (``σ·σ′`` has units ``cm⁻²``, not
  ``cm⁻¹``), so its closure is a promotion-layer (S3) concern, not S1; this
  mixin is its designated home, deferred until it does work
  (``coding-elegance`` Pattern 6).

Mixed in BEFORE the storage base
(``CrossSectionField(CoefficientRole, ScalarField)``) so it occupies the same
MRO slot as :class:`FluxRole` does for the state leaves — even though, in S1,
it contributes no override there.

References
----------

* ``.claude/agent-memory/cross-domain-attacker/coefficient_field_promotion_frames.md``
  — Frame 2 (the coefficient cone algebra is NOT the flux torsor).
* ``.claude/plans/issue_257_coefficient_field_promotion.md`` — S1.
* ``coding-elegance`` Pattern 1 (the role reads like the domain) + Pattern 6
  (the field·field multiplier product is deferred to S3, where it does work).
"""
from __future__ import annotations


class CoefficientRole:
    r"""Role marker for the cross-section coefficient leaves.

    The complement of :class:`~orpheus.transport.fields._flux_role.FluxRole`:
    where ``FluxRole`` removes ``flux + flux`` and retypes ``flux − flux``,
    ``CoefficientRole`` keeps the plain :class:`~orpheus.numerics.field.Field`
    vector-space algebra (``σ + σ`` legitimate, ``Σ = 0`` is the origin) and
    adds **no** affine gate. It overrides nothing in S1 — the *absence* of the
    torsor gate is its content. It is the designated home for the future
    pointwise field·field product (the multiplier-algebra multiplication, #257
    S3), deferred until that product is built.
    """
