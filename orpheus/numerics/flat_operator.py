r"""The flattened realization of a typed-carrier operator.

:class:`FlattenedOperator` conjugates a
:class:`~orpheus.numerics.operator.LinearOperator` whose carrier is a
**flat-vector-like composite** (anything exposing the duck-typed
``to_flat()`` / ``from_flat(flat, template)`` pair — the
:data:`~orpheus.numerics.iteration.FlatVectorLike` protocol; the
canonical instance is the transport
``FullField``) into an operator on **bare 1-D ndarrays**:

.. math::

    A_{\rm flat} \;=\; \operatorname{flat} \circ A \circ
        \operatorname{unflat} .

Why it exists (#290 P4 — the diffusion resolvent)
=================================================

:meth:`~orpheus.numerics.operator.LinearOperator.as_matrix` feeds its
apply-to-basis loop with **ndarray** basis columns, so a typed-carrier
operator cannot materialize directly (the honest-scope note on
``as_matrix``: SN composites "stay matrix-free" — they are also far
above any sane size gate). The diffusion loss
:math:`A = L + C - S - B` is the counterexample the note anticipated:
typed on the scalar composite AND small enough for the dense direct
resolvent. This adapter is the one bridge — the flat-layout knowledge
stays where it lives (the carrier's own ``to_flat`` / ``from_flat``,
Pattern 2), and the resolvent is spelled

.. code-block:: python

    A_flat = FlattenedOperator(A, template)     # template: a zero composite
    A_inv  = MatrixInverseOperator(A_flat)      # eager LU — the #290 resolvent
    x_flat = A_inv.apply(q.to_flat())
    x      = type(template).from_flat(x_flat, template)

``domain`` / ``codomain`` need no re-derivation: a composite space
(:class:`~orpheus.numerics.spaces.full_field_space.FullFieldSpace`)
already carries the flat direct-sum dimension as its ``shape``
identity, so the inner operator's spaces are reported verbatim and
``as_matrix`` derives ``basis_shape = (n_interior + n_trace,)`` for free.

The TEMPLATE is input-side metadata only: basis vectors are unflattened
through it (shapes, spaces, meshes, concrete composite class), while the
output side just re-flattens whatever composite ``inner.apply`` emits —
role changes (flux in, source out) pass through untouched.

Structural axes: apply-only. The flattened operator is a
*serialization* of ``inner``, not new structure — inversion goes
through an explicit
:class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`
over it (values, not structure), and the metric-correct adjoint lives
on the TYPED operator where the block metrics are visible.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Optional

import numpy as np

from orpheus.numerics.operator import LinearOperator

if TYPE_CHECKING:
    from orpheus.numerics.space import FunctionSpace


__all__ = ["FlattenedOperator"]


class FlattenedOperator(LinearOperator):
    r"""``flat ∘ A ∘ unflat`` — a typed-carrier operator on bare 1-D vectors.

    (Deliberately unparameterized, the ndarray-carrier family
    convention — the ``IsotropicScattering`` / ``DiagonalOperator``
    precedent: the carrier TypeVars are bounded by the typed-vector
    protocol a bare ``ndarray`` does not satisfy.)

    Parameters
    ----------
    inner : LinearOperator
        The typed-carrier operator (e.g. the diffusion loss
        ``L + C − S − B`` on the scalar composite).
    template : FlatVectorLike
        A composite instance supplying the unflatten metadata (shapes,
        spaces, meshes, concrete class) — canonically a zero composite
        on the operator's mesh. Must satisfy
        ``from_flat(template.to_flat(), template) ≡ template``
        structurally (the round-trip invariant of the flat protocol).

    Notes
    -----
    ``apply(x)`` accepts any array raveling to the template's flat size
    and returns the 1-D flattened action. The input basis element is
    reconstructed as the SAME concrete class as ``template`` (via the
    carrier's ``from_flat``); the output is whatever composite
    ``inner.apply`` emits, re-flattened via ITS ``to_flat``.
    """

    def __init__(self, inner: "LinearOperator", template: Any) -> None:
        self.inner = inner
        self.template = template
        #: The flat dimension — frozen from the template at construction
        #: (one ``to_flat`` allocation; also the input-size guard).
        self.n_flat: int = int(np.asarray(template.to_flat()).size)

    @property
    def domain(self) -> "Optional[FunctionSpace]":
        # The inner operator's space identity IS flat-compatible: a
        # composite FullFieldSpace's ``shape`` is the direct-sum flat
        # dimension, so as_matrix derives basis_shape from it directly.
        return getattr(self.inner, "domain", None)

    @property
    def codomain(self) -> "Optional[FunctionSpace]":
        return getattr(self.inner, "codomain", None)

    def apply(self, x: Any, /) -> np.ndarray:
        r"""Return ``inner.apply(unflat(x)).to_flat()``."""
        flat = np.asarray(x, dtype=float).ravel()
        if flat.size != self.n_flat:
            raise ValueError(
                f"FlattenedOperator.apply: input size {flat.size} does "
                f"not match the template's flat size {self.n_flat}."
            )
        carrier = type(self.template).from_flat(flat, self.template)
        return np.asarray(self.inner.apply(carrier).to_flat())

    # ── Assembly passthrough (stencil-assembly 2b) ────────────────────
    #
    # The flat serialization is TRANSPARENT to structural assembly: a
    # typed operator's sparse emission is ALREADY in the flat layout
    # (the emitters scatter into the carrier's own ``to_flat`` DOF
    # numbering — the 2-P0 never-re-derive rule), so the flattened
    # view's assembly IS the inner operator's, dimension-checked against
    # the template. Through this passthrough the ``as_matrix``
    # delegation (R2) reaches the diffusion resolvent: once every loss
    # leaf emits, ``MatrixInverseOperator(FlattenedOperator(A, t))``
    # LU-factors the assembled matrix instead of probing n columns.

    @property
    def is_assemblable(self) -> bool:
        return self.inner.is_assemblable

    def assemble(self) -> "Any":
        r"""Return the inner operator's assembly (the flat layouts coincide).

        Raises
        ------
        MissingAssembly
            Propagated from composer recursion when a leg of ``inner``
            cannot emit (guarded here via the checked bridge).
        ValueError
            The inner emission's dimensions contradict the template's
            flat size — a layout drift between an emitter and the
            carrier's ``to_flat`` (must be impossible by construction;
            checked because silence here would corrupt every consumer).
        """
        from orpheus.numerics.operator import MissingAssembly, assemblable

        if not assemblable(self.inner):
            raise MissingAssembly(
                f"FlattenedOperator.assemble: the wrapped "
                f"{type(self.inner).__name__} is not assemblable."
            )
        assembled = self.inner.assemble()
        rows, cols = assembled.shape
        if rows != self.n_flat or cols != self.n_flat:
            raise ValueError(
                f"FlattenedOperator.assemble: the inner emission is "
                f"{rows}×{cols} but the template's flat size is "
                f"{self.n_flat} — an emitter drifted from the carrier's "
                f"to_flat layout."
            )
        return assembled

    def __repr__(self) -> str:
        return (
            f"FlattenedOperator({self.inner!r}, n_flat={self.n_flat})"
        )
