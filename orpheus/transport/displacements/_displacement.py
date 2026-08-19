r"""The :class:`Displacement` marker + convergence-diagnostics mixin.

A **flux displacement** :math:`\Delta\psi = \psi^{(i)} \ominus \psi^{(i-1)}` is
the iterate INCREMENT — an element of the difference vector space :math:`V` of
the affine flux space :math:`A` (see
:class:`~orpheus.transport.fields._flux_role.FluxRole`). It shares its flux
sibling's storage family, :class:`~orpheus.numerics.space.FunctionSpace`, and
metric (it IS the tangent space of that same function space), but is a DISTINCT
physical kind: a *step*, not a *state*.

Why a distinct type (the affine-axiom argument)
===============================================

Flux states form an affine space with no natural origin: ``state + state`` is
undefined (it lands off the affine subspace; treating :math:`\psi=0` as an
origin is a category error). The legitimate operations are ``state ⊖ state →
displacement`` (the difference, in :math:`V`) and ``state ⊕ displacement →
state`` (the torsor action :math:`A \times V \to A`, e.g. the relaxation
update). Typing the increment as the flux STATE type both admits the illegal
``state + state`` AND strands the contraction data (:math:`\rho`, the
a-posteriori bound) with no home — the displacement is the ONLY object that
knows "previous" / "step". This mixin gives it the methods
:class:`~orpheus.numerics.field.Field` (a state) structurally cannot carry.

The vector-space algebra (``+`` between displacements, scalar ``·``, unary
``-``) is INHERITED from :class:`~orpheus.numerics.field.Field` — :math:`V` is a
genuine vector space (unlike :math:`A`). This mixin adds ONLY the diagnostics
and serves as the family marker (``isinstance(x, Displacement)``).

References
----------

* ``.claude/agent-memory/cross-domain-attacker/issue_208_delta_psi_affine_frames.md``
  — the affine / torsor / Banach-fixed-point / Krylov-dual triad.
* ``.claude/agent-memory/numerics-investigator/issue_208_flux_displacement_residual_typing_debug_value.md``
  — the method catalogue + the :math:`c\to 1` false-convergence finding.
* Adams, M.L. & Larsen, E.W. (2002). *Fast iterative methods for
  discrete-ordinates particle transport calculations.* Prog. Nucl. Energy 40(1)
  — the SI contraction factor :math:`\rho \approx \max \Sigma_s/\Sigma_t`.
"""
from __future__ import annotations

from typing import ClassVar


class Displacement:
    r"""Marker for difference-space (tangent) fields.

    Mixed in BEFORE the storage base on each displacement leaf
    (``AngularDisplacement(Displacement, AngularField)``, …). Carries no
    storage and no dunder algebra — those come from the storage base /
    :class:`~orpheus.numerics.field.Field`.

    ⛔ RETIRED SURFACE (campaign 1 CS3 step 1, 2026-08-19): the convergence
    diagnostics this marker used to carry (``contraction_ratio`` /
    ``true_error_estimate`` / ``where_largest``) RELOCATED — ρ and the
    :math:`c\to 1` geometric-tail estimate now derive from
    :attr:`~orpheus.numerics.convergence.IterationRecord.increment_norms`
    on the iteration record (space-norm convention, user ruling), and the
    per-entry magnitude map is
    :meth:`~orpheus.numerics.field.Field.where_largest`, a property of any
    field. What remains here is the Rep-keyed flux↔displacement pairing
    (:meth:`sibling_of`), which retires with the type family at CS3 step 3.
    """

    #: Registry ``Rep → Displacement`` populated by :meth:`__init_subclass__` as
    #: each displacement leaf is defined. The ``displacements`` package
    #: ``__init__`` imports every leaf eagerly, so the registry is complete once
    #: that package is imported (W-H: the flux↔displacement pairing is DERIVED
    #: from the shared Rep, not hand-set on each flux leaf).
    _BY_REP: ClassVar[dict[type, type["Displacement"]]] = {}

    def __init_subclass__(cls, **kwargs) -> None:
        super().__init_subclass__(**kwargs)
        Displacement._BY_REP[_carrier_rep(cls)] = cls

    @classmethod
    def sibling_of(cls, carrier_cls: type) -> type["Displacement"]:
        r"""The displacement role-class sharing ``carrier_cls``'s Rep.

        A flux and its displacement are two role-classes over the SAME
        Field-family Rep (``ScalarFlux`` / ``ScalarDisplacement`` both over
        ``ScalarField``); the lookup is keyed by that shared Rep, so the pairing
        is STRUCTURAL, not nominal — it handles the ``HarmonicMomentFlux`` ↔
        ``MomentDisplacement`` name asymmetry a name-mangling derive could not.
        """
        return Displacement._BY_REP[_carrier_rep(carrier_cls)]

def _carrier_rep(cls: type) -> type:
    r"""The Field-family Rep of a flux / displacement role-leaf.

    Both ``XFlux(FluxRole, XField)`` and ``XDisplacement(Displacement, XField)``
    carry exactly one :class:`~orpheus.numerics.field.Field`-family base — the
    storage representation ``XField`` they share. The role mixin ``FluxRole`` and
    the ``Displacement`` marker are not ``Field`` subclasses, so they are
    excluded; that shared Rep is the key the flux↔displacement pairing is derived
    from (W-H).
    """
    from orpheus.numerics.field import Field

    reps = [
        base
        for base in cls.__bases__
        if isinstance(base, type)
        and issubclass(base, Field)
        and not issubclass(base, Displacement)
    ]
    if len(reps) != 1:
        raise TypeError(
            f"{cls.__name__}: a flux/displacement role-leaf must derive from "
            f"exactly one Field-family Rep; found {[r.__name__ for r in reps]}."
        )
    return reps[0]
