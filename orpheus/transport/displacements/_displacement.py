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

import numpy as np


class Displacement:
    r"""Marker + convergence-diagnostics for difference-space (tangent) fields.

    Mixed in BEFORE the storage base on each displacement leaf
    (``AngularDisplacement(Displacement, AngularField)``, …). Carries no
    storage and no dunder algebra — those come from the storage base /
    :class:`~orpheus.numerics.field.Field`. It provides the three diagnostics a
    flux state cannot:

    * :meth:`contraction_ratio` — the Banach factor :math:`\rho`.
    * :meth:`true_error_estimate` — :math:`\lVert\Delta\psi\rVert/(1-\rho)`, the
      fix for the :math:`c\to 1` false-convergence.
    * :meth:`where_largest` — the per-cell / per-ordinate convergence map.

    All three use the space-induced metric norm :attr:`Field.l2` (consistent
    with the flux because the displacement shares the flux's space). NOTE these
    are DIAGNOSTICS — the SI stopping test stays the flat
    ``np.linalg.norm(to_flat())`` for bit-identity; switching it to ``.l2``
    would re-interpret ``conv_tol``.
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

    def contraction_ratio(self, previous: "Displacement") -> float:
        r"""The Banach contraction factor :math:`\rho \approx
        \lVert\Delta\psi^{(i)}\rVert / \lVert\Delta\psi^{(i-1)}\rVert`.

        :math:`\rho > 1` diverges (wrong fixed point / unstable scheme);
        :math:`\rho \approx 1` stalled (the :math:`c\to 1` slow mode, curvilinear
        / reflective slow modes); :math:`\rho < 1` healthy. Turns the
        :math:`\rho`-blind :math:`\lVert\Delta\psi\rVert` stopping test honest.

        Parameters
        ----------
        previous : Displacement
            The previous iterate's displacement :math:`\Delta\psi^{(i-1)}`.

        Raises
        ------
        ZeroDivisionError
            If the previous displacement has zero norm (the iteration was
            already at the fixed point — the ratio is undefined).
        """
        prev_norm = previous.l2  # type: ignore[attr-defined]
        if prev_norm == 0.0:
            raise ZeroDivisionError(
                "contraction_ratio: the previous displacement has zero norm "
                "(already at the fixed point) — the ratio is undefined."
            )
        return self.l2 / prev_norm  # type: ignore[attr-defined]

    def true_error_estimate(self, contraction_ratio: float) -> float:
        r"""The geometric-tail true error :math:`\lVert\Delta\psi\rVert/(1-\rho)`.

        The core :math:`c\to 1` false-convergence fix: the bare increment
        :math:`\lVert\Delta\psi\rVert` understates the distance to the fixed
        point by :math:`1/(1-\rho)`. At :math:`\rho\approx 0.99` (e.g.
        :math:`c=0.99`) the true error is ~100× the increment, so a
        :math:`\lVert\Delta\psi\rVert < \text{tol}` "convergence" can be ~100·tol
        from the solution. This estimate surfaces that stall.

        Parameters
        ----------
        contraction_ratio : float
            The Banach factor :math:`\rho \in [0, 1)`, typically from
            :meth:`contraction_ratio`.

        Raises
        ------
        ValueError
            If ``contraction_ratio`` is not in :math:`[0, 1)` (a non-contracting
            iteration has no finite geometric-tail error estimate).
        """
        if not (0.0 <= contraction_ratio < 1.0):
            raise ValueError(
                f"true_error_estimate: contraction_ratio must be in [0, 1) for "
                f"a finite geometric-tail estimate; got {contraction_ratio!r} "
                f"(>= 1 means the iteration is not contracting)."
            )
        return self.l2 / (1.0 - contraction_ratio)  # type: ignore[attr-defined]

    def where_largest(self, k: int = 1) -> list[tuple[int, ...]]:
        r"""The ``k`` index tuples with the largest :math:`|\Delta\psi|`.

        The per-cell / per-group / per-ordinate convergence map: WHICH entries
        are not converging (pole-cell resonance, material-interface slow modes,
        a lagging group). Indices are into the leaf's ``values`` layout (e.g.
        ``(n, g, ix, iy)`` for an angular displacement), largest first.

        Parameters
        ----------
        k : int
            How many of the largest-magnitude entries to return (default 1).
        """
        flat = np.abs(np.asarray(self.values)).ravel()  # type: ignore[attr-defined]
        k = max(1, min(int(k), flat.size))
        top = np.argpartition(flat, -k)[-k:]
        top = top[np.argsort(flat[top])[::-1]]  # largest first
        shape = np.asarray(self.values).shape  # type: ignore[attr-defined]
        return [tuple(int(i) for i in np.unravel_index(j, shape)) for j in top]


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
