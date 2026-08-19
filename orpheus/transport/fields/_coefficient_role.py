r"""The :class:`CoefficientRole` mixin — the coefficient algebra role marker.

A cross section is not a *state* — it is a **coefficient**, the symbol of a
zeroth-order operator. The grand report (§5.5–5.7) names cross-section fields
``CoefficientField``\ s and the operators they become *promotions*
(``C = M[σ_t]``, the multiplier-algebra embedding ``M: L^∞ → B(L²)``).
:class:`CoefficientRole` is the role marker on the field side of that
promotion. (Until campaign 1 CS3, 2026-08-19, it was the **exact
complement** of the retired ``FluxRole`` affine mixin:)

==========================  ====================================  ===============================
Role                        Algebra                               Mints
==========================  ====================================  ===============================
``FluxRole`` (state)        ⛔ RETIRED (campaign 1 CS3,           nothing — ``±`` return the
                            2026-08-19): flux lives in V,         same class (fiber-guarded)
                            same vector algebra as below
``CoefficientRole`` (coef)  plain vector space — HAS an origin    nothing new (keeps ``Field``)
                            (``σ + σ`` legitimate)
==========================  ====================================  ===============================

(The overturned doctrine held that flux leaves form an affine space with no
natural origin, ``FluxRole`` overriding their additive algebra. The CS3
ruling inverted it: flux lives in the positive cone K ⊂ V, and THIS
module's posture — cone as tested property, algebra as plain vector
space — became every family's.) The coefficients always did
form a genuine **vector space with an origin** (the zero cross section
``Σ = 0`` IS a coefficient — it promotes to ``M_0 = ZeroOperator``), closed
under ``+`` (homogenisation ``Σ_mix = Σ_m N_m Σ_m`` is a number-density-weighted
sum), unary ``−``, and scalar ``·`` — i.e. exactly the plain
:class:`~orpheus.numerics.field.Field` vector-space dunders. So
:class:`CoefficientRole` overrides **nothing** — historically the *absence*
of an affine gate was its content relative to ``FluxRole``; today it is
simply the family's shared shape. ``σ − σ′`` stays the same class
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
  a flux state without a unit check (the state leaves lost their role mixin
  at CS3; class identity itself now tags them).
* **The future multiplier product.** When ``MultiplicationOperator`` lands
  (#257 S3), the law ``M_f @ M_g = M_{f·g}`` needs the pointwise *field·field*
  product ``f·g`` — the commutative-algebra multiplication that the flux never
  has. That product is dimensionally graded (``σ·σ′`` has units ``cm⁻²``, not
  ``cm⁻¹``), so its closure is a promotion-layer (S3) concern, not S1; this
  mixin is its designated home, deferred until it does work
  (``coding-elegance`` Pattern 6).

Mixed in BEFORE the storage base
(``CrossSectionField(CoefficientRole, ScalarField)``) — the MRO slot the
retired ``FluxRole`` used to occupy on the state leaves — even though, in
S1, it contributes no override there.

References
----------

* ``.claude/agent-memory/cross-domain-attacker/coefficient_field_promotion_frames.md``
  — Frame 2 (the coefficient cone algebra vs the then-flux-torsor; the
  torsor was overturned at campaign 1 CS3 — this doctrine generalised).
* ``.claude/plans/issue_257_coefficient_field_promotion.md`` — S1.
* ``coding-elegance`` Pattern 1 (the role reads like the domain) + Pattern 6
  (the field·field multiplier product is deferred to S3, where it does work).
"""
from __future__ import annotations


class CoefficientRole:
    r"""Role marker for the cross-section coefficient leaves.

    Keeps the plain :class:`~orpheus.numerics.field.Field` vector-space
    algebra (``σ + σ`` legitimate, ``Σ = 0`` is the origin) and adds **no**
    gate. (Historically the complement of the retired ``FluxRole`` affine
    mixin — the *absence* of the torsor gate was its content, and at
    campaign 1 CS3 that shape became every family's.) It overrides nothing
    in S1. It is the designated home for the future
    pointwise field·field product (the multiplier-algebra multiplication, #257
    S3), deferred until that product is built.
    """
