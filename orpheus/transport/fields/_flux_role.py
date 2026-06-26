r"""The :class:`FluxRole` mixin — affine-point algebra for flux-STATE leaves.

Flux states (:class:`AngularFlux`, :class:`ScalarFlux`,
:class:`HarmonicMomentFlux`, :class:`BoundaryFlux`) form an **affine space**
:math:`A` over a distinct **difference vector space** :math:`V` (the
displacements — :mod:`orpheus.transport.displacements`). An affine space has no
natural origin, so the additive algebra of a flux leaf is NOT the plain
vector-space algebra a generic :class:`~orpheus.numerics.field.Field` inherits.
This mixin overrides it to the torsor algebra:

.. math::

    \psi_2 \ominus \psi_1 &\;\to\; \Delta\psi \in V
        \quad&&\text{(the iterate increment; the ONLY mint of a displacement)}\\
    \psi \oplus \Delta\psi &\;\to\; \psi' \in A
        \quad&&\text{(the torsor action } A\times V\to A\text{; the update step)}\\
    \psi_1 \oplus \psi_2 &\;\to\; \bot
        \quad&&\text{(no origin — } \textstyle\sum\lambda_i = 2 \text{ lands off } A)\\
    \textstyle\sum_i \lambda_i \psi_i,\ \sum_i\lambda_i=1 &\;\to\; \psi \in A
        \quad&&\text{(affine combination — the only legal multi-flux blend)}

The relaxation :math:`\omega\,\psi_{\rm new} + (1-\omega)\,\psi_{\rm old}` IS an
affine combination (:math:`\sum\lambda=1`), equivalently the torsor form
:math:`\psi_{\rm old} \oplus \omega\,(\psi_{\rm new}\ominus\psi_{\rm old})`.

**Scalar scaling is UNTOUCHED.** ``__mul__`` / ``__rmul__`` / ``__truediv__`` /
``__neg__`` take a SCALAR and are inherited from :class:`Field` — eigenvalue
normalisation (``ψ / k``), :math:`\omega\cdot\Delta\psi`, etc. all survive. The
gate forbids ONLY the binary ``flux + flux``; ``__sub__`` is RETYPED (→
displacement) and the cross-class / cross-mesh guards are inherited unchanged.

The sibling displacement class is DERIVED from the shared Rep (W-H): a flux
``XFlux(FluxRole, XField)`` and its displacement ``XDisplacement(Displacement,
XField)`` are two roles over the SAME ``XField`` Rep, so ``_DISPLACEMENT_CLS`` is
looked up by that Rep (:meth:`Displacement.sibling_of`), not hand-set on each
leaf. The mint copies the flux's dataclass init-fields (``values`` / ``space`` /
``mesh``, plus ``L`` for the moment family) into that class, swapping ``values``
— uniform across families and handling the moment leaf's ``L`` that a
``cls.from_mesh(...)`` engine could not. The displacement SHARES the flux's
:class:`~orpheus.numerics.space.FunctionSpace` (the tangent vector lives in the
flux's own function space); the CLASS identity is the role gate.

References
----------

* ``.claude/agent-memory/cross-domain-attacker/issue_208_delta_psi_affine_frames.md``
  — the affine / torsor / Banach / Krylov triad (the design source).
* ``.claude/plans/issue_208_residual_typing_closeout.md`` §3a — the locked
  mechanism (field-copy mint, share-the-flux-space, TimedFullField transparency).
* ``coding-elegance`` Pattern 4 (illegal states unrepresentable) — ``flux+flux``
  is now unspellable; Pattern 1 (the torsor algebra reads like the math).
"""
from __future__ import annotations

from dataclasses import fields as _dataclass_fields, replace
from typing import Sequence


class FluxRole:
    r"""Affine-point additive algebra for the flux-state role leaves.

    Mixed in BEFORE the storage base (``AngularFlux(FluxRole, AngularField)``)
    so its ``__add__`` / ``__sub__`` take MRO precedence over
    :class:`~orpheus.numerics.field.Field`'s vector-space versions. Provides
    nothing else — storage, ``_check_partner``, scalar scaling, and diagnostics
    come from the storage base / :class:`Field`.
    """

    @property
    def _DISPLACEMENT_CLS(self) -> type:
        r"""The sibling displacement class, DERIVED from the shared Rep (W-H).

        A flux and its displacement are two role-classes over the SAME
        Field-family Rep (``ScalarFlux`` / ``ScalarDisplacement`` both over
        ``ScalarField``), so the pairing is structural — looked up via
        :meth:`Displacement.sibling_of` keyed by this flux's Rep, NOT hand-set on
        each leaf. Handles the ``HarmonicMomentFlux`` ↔ ``MomentDisplacement``
        name asymmetry a name-mangling derive could not.
        """
        from orpheus.transport.displacements._displacement import Displacement

        return Displacement.sibling_of(type(self))

    # ── The affine algebra ──────────────────────────────────────────────

    def __sub__(self, other):
        r"""``flux ⊖ flux → displacement`` (the iterate increment).

        Validates via the flux's own :meth:`_check_partner` (same class, space,
        mesh, and ``L`` for the moment leaf), then mints the sibling
        displacement carrying the signed difference ``self − other``.
        """
        self._check_partner(other)  # type: ignore[attr-defined]
        return self._mint_displacement(self.values - other.values)  # type: ignore[attr-defined]

    def __add__(self, other):
        r"""``flux ⊕ displacement → flux`` (torsor); ``flux ⊕ flux → TypeError``.

        The torsor action recovers a flux from a flux + a displacement (the
        update step). Adding two fluxes is forbidden — an affine space has no
        origin (:math:`\sum\lambda=2` lands off :math:`A`); the error names the
        legal alternatives.
        """
        if isinstance(other, self._DISPLACEMENT_CLS):
            self._check_torsor_partner(other)
            return replace(self, values=self.values + other.values)  # type: ignore[attr-defined]
        if type(self) is type(other):
            raise TypeError(
                f"cannot add two {type(self).__name__} states: flux states form "
                f"an affine space with no origin, so '+' between two fluxes is "
                f"undefined (Σλ = 2 lands off the affine subspace). Use "
                f"{type(self).__name__}.affine_combination([(λᵢ, ψᵢ)]) with "
                f"Σλᵢ = 1 for a relaxation / convex blend, ψ + Δψ "
                f"(flux + displacement) for the update step, or add the "
                f"rate-density L-actions (sources / residuals)."
            )
        # Cross-class / cross-role: defer to the standard gate, which raises the
        # cross-class TypeError (same units ≠ same meaning).
        self._check_partner(other)  # type: ignore[attr-defined]
        raise TypeError(  # pragma: no cover — _check_partner always raises here
            f"{type(self).__name__} + {type(other).__name__} is undefined."
        )

    @classmethod
    def affine_combination(cls, pairs: Sequence[tuple[float, "FluxRole"]]):
        r"""``Σλᵢ ψᵢ`` with ``Σλᵢ = 1`` → flux (the only legal multi-flux blend).

        The partition-of-unity constraint :math:`\sum\lambda_i = 1` IS the affine
        constraint (it keeps the result on the affine subspace). The relaxation
        :math:`\omega\psi_{\rm new} + (1-\omega)\psi_{\rm old}` is the canonical
        instance. All operand fluxes must be the same class / space / mesh.

        Parameters
        ----------
        pairs : sequence of (float, flux)
            The ``(λᵢ, ψᵢ)`` pairs. ``Σλᵢ`` must equal 1 (within 1e-12).

        Raises
        ------
        ValueError
            If ``pairs`` is empty or :math:`\sum\lambda_i \neq 1`.
        """
        if len(pairs) == 0:
            raise ValueError(
                "affine_combination requires at least one (λ, ψ) pair."
            )
        coeffs = [float(c) for c, _ in pairs]
        if abs(sum(coeffs) - 1.0) > 1e-12:
            raise ValueError(
                f"affine_combination requires Σλᵢ = 1 (the partition-of-unity "
                f"affine constraint — keeps the blend on the affine subspace); "
                f"got Σλ = {sum(coeffs)!r}. A relaxation ω·ψ_new + (1−ω)·ψ_old "
                f"satisfies it; an arbitrary weighted sum does not."
            )
        fluxes = [f for _, f in pairs]
        first = fluxes[0]
        for f in fluxes[1:]:
            first._check_partner(f)  # type: ignore[attr-defined]
        # The affine blend Σλᵢ·ψᵢ (the partition of unity Σλ=1 is enforced
        # above). Seed the reduction with the first scaled flux's values — NOT
        # the builtin ``sum``'s int ``0`` start — so the accumulation is over
        # float arrays throughout and the intermediate is the named blend.
        blended_values = sum(
            (c * f.values for c, f in zip(coeffs[1:], fluxes[1:])),  # type: ignore[attr-defined]
            start=coeffs[0] * first.values,  # type: ignore[attr-defined]
        )
        return replace(first, values=blended_values)  # type: ignore[arg-type]

    # ── Mint / torsor helpers ───────────────────────────────────────────

    def _mint_displacement(self, new_values):
        r"""Build the sibling displacement carrying ``new_values``.

        Copies this flux's dataclass init-fields (``values`` / ``space`` /
        ``mesh``, plus ``L`` for the moment leaf) into ``_DISPLACEMENT_CLS``,
        swapping ``values``. The displacement shares this flux's ``space``
        object, so it lives in the flux's own function space (same metric).
        """
        init = {
            f.name: getattr(self, f.name)
            for f in _dataclass_fields(self) if f.init  # type: ignore[arg-type]
        }
        init["values"] = new_values
        return self._DISPLACEMENT_CLS(**init)

    def _check_torsor_partner(self, displacement) -> None:
        r"""Reject an ill-formed torsor action ``flux ⊕ displacement``.

        The partner is already known to be the sibling displacement class. It
        shares the flux's space; verify mesh identity (the field is mesh-bound)
        and shape compatibility — a displacement from a different mesh cannot
        act on this flux.
        """
        self_mesh = getattr(self, "mesh", None)
        other_mesh = getattr(displacement, "mesh", None)
        if self_mesh is not other_mesh:
            raise ValueError(
                f"{type(self).__name__} ⊕ {type(displacement).__name__} across "
                f"distinct SNMesh instances is forbidden — the field is "
                f"mesh-bound."
            )
        if self.values.shape != displacement.values.shape:  # type: ignore[attr-defined]
            raise ValueError(
                f"{type(self).__name__} ⊕ {type(displacement).__name__} shape "
                f"mismatch: {self.values.shape} vs {displacement.values.shape}."  # type: ignore[attr-defined]
            )
