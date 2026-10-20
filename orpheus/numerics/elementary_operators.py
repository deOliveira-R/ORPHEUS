r"""The elementary operators: concrete leaves built on the operator algebra.

:mod:`orpheus.numerics.operator` holds the ALGEBRA of linear operators:
:class:`~orpheus.numerics.operator.LinearOperator`, its sums, products,
tensor products, the inverse, the adjoint
:math:`\sharp\circ\mathrm{dual}\circ\flat`, the dual and the Riesz operators.
This module holds the concrete operators and constructors that are
INSTANCES built on that algebra and which the algebra itself does not need:

- the structural index maps :class:`PermutationOperator` and
  :class:`TraceRestrictionOperator`, with their shared extent check
  :func:`checked_space_extent`;
- the pointwise multiplier :class:`DiagonalOperator`;
- the rank-one dyad :class:`RankOneOperator` and its constructor :func:`outer`.

In the concept order of :mod:`orpheus.numerics` this module sits one step
downstream of ``operator``: it imports the algebra, and the algebra never
imports it (``tests/gates/numerics/test_concept_order.py``). The definitions
moved here verbatim from ``operator.py``; behaviour is unchanged. The sum of
tensor products that lived here retired at package S2: a sum of tensor
products is an :class:`~orpheus.numerics.operator.OperatorSum` whose words
are :class:`~orpheus.numerics.operator.TensorProductOperator` (#260 decides
whether separability earns its own type).

The zero maps :class:`~orpheus.numerics.operator.ZeroOperator` and
:class:`~orpheus.numerics.operator.ZeroMorphism` moved here too, and went
back the same day (the user's ruling, 2026-10-19): the zero of
:math:`\mathrm{Hom}(V, W)` is the additive identity of the algebra, as
:class:`~orpheus.numerics.operator.IdentityOperator` is its multiplicative
one, so it belongs to the algebra.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Optional, cast

import numpy as np

from orpheus.numerics.operator import (
    Domain,
    InverseOperator,
    LinearOperator,
    MissingAdjoint,
    NotInvertible,
    PointwiseOperator,
    TensorProductOperator,
    adjointable,
)
from orpheus.numerics.vector import Vector

if TYPE_CHECKING:
    from orpheus.numerics.functional import Functional
    from orpheus.numerics.space import FunctionSpace

__all__ = [
    "PermutationOperator",
    "TraceRestrictionOperator",
    "checked_space_extent",
    "DiagonalOperator",
    "RankOneOperator",
    "outer",
]


class PermutationOperator(LinearOperator):
    r"""Index permutation along a configurable axis: :math:`(P x)_i = x_{\pi(i)}`.

    For a permutation :math:`\pi : \{0, \ldots, N-1\} \to \{0, \ldots, N-1\}`
    represented as an integer array ``perm`` of length :math:`N`, the
    apply action gathers entries along ``axis`` according to ``perm``:

    .. math::

        (P\,x)_{i_0 \ldots i_{a-1}\,j\,i_{a+1} \ldots} \;=\;
        x_{i_0 \ldots i_{a-1}\,\pi(j)\,i_{a+1} \ldots}

    The transpose is the inverse permutation :math:`\pi^{-1}`, computed
    once at construction via ``np.argsort(perm)``.

    ⭐ **Ask the algebra whether it is an involution; there is no flag.**
    Until G6.3 step 5 this class carried an ``is_involution`` attribute,
    set at construction from ``perm[perm] == np.arange(N)``. It is retired,
    and the retirement is the *point* rather than a tidy-up: an involution
    is a claim about :math:`P \circ P`, so it needs that composition to
    exist, and the composition is exactly what the operator algebra already
    knows how to form and to refuse. Spell it::

        (P @ P).apply(x)        # same-space: compare against x
        P @ P                   # cross-space: IncompatibleOperatorComposition

    The second line is the whole argument. A permutation bound
    :math:`\Gamma_+ \to \Gamma_-` — every specular kernel the SN realizer
    builds — has no square, and a stored ``bool`` had to answer *something*
    anyway. `[M]` for ONE physical law (a mirror about ``x`` on ``xmin``)
    the raw index test answered ``True`` on ``gauss_legendre(4/8)``,
    ``product(4,4)`` and ``level_symmetric(6)`` and ``False`` on
    ``lebedev(17)`` — an answer tracking the quadrature's local index
    ordering rather than the mirror, and unfalsifiable at every row because
    the flag's documented purpose (self-adjointness in the unweighted inner
    product) is undefined between two different spaces. Routing the
    question through ``@`` replaces a value that could be wrong with a
    composition that cannot be formed. A full-space rule like the mirror
    pairing
    :meth:`~orpheus.numerics.quadrature.Quadrature.ordinate_permutation`
    derives for a reflection IS an involution on
    :math:`\{0, \ldots, N-1\}`, and asks and answers as one; periodic
    shifts and rotational reorderings are not.

    A permutation is always invertible (:math:`P^{-1} = P^T`), and its
    inverse is ALGEBRA-CLOSED: :meth:`inverse` returns the inverse
    permutation as a first-class :class:`PermutationOperator` (#226
    taxonomy step 1) whose ``apply`` is the same
    ``np.take(·, inverse_perm)`` gather as :meth:`apply_transpose`.

    Parameters
    ----------
    perm
        1-D integer array of length :math:`N` whose entries are a
        permutation of :math:`\{0, \ldots, N-1\}`. Validated at
        construction; rejecting duplicates and out-of-range entries
        with :class:`ValueError`.
    axis
        Tensor axis along which the permutation acts. The operator
        broadcasts on every other axis.
    domain, codomain : FunctionSpace, optional
        The spaces this permutation maps between. For the SN deck
        transformation — a specular mirror realized as a **length-1 chain**
        (campaign step G6.3, issue **#330**) — these are
        :math:`\Gamma_+(f)` and :math:`\Gamma_-(f)`, and binding them is
        what makes ``@`` refuse a mis-composed chain and ``.H`` the
        metric-aware Hilbert adjoint rather than the bare transpose.
        Both extents are checked against :attr:`n` along :attr:`axis` at
        construction, because a mis-bound space is SILENT at apply-time.

    Attributes
    ----------
    perm : np.ndarray
        Forward permutation :math:`\pi`, as 1-D ``intp`` array.
    inverse_perm : np.ndarray
        Inverse permutation :math:`\pi^{-1}`, precomputed via
        :func:`numpy.argsort`.
    axis : int
        The tagged tensor axis.
    n : int
        Length of the permuted axis.
    """

    def __init__(
        self, perm: np.ndarray, axis: int = 0,
        *,
        domain: Optional[FunctionSpace] = None,
        codomain: Optional[FunctionSpace] = None,
    ) -> None:
        perm = np.asarray(perm, dtype=np.intp)
        if perm.ndim != 1:
            raise ValueError(
                f"PermutationOperator perm must be 1-D; got shape {perm.shape}"
            )
        n = perm.size
        # Validate: perm is a true permutation of {0, ..., n-1}.
        if n == 0 or not (
            perm.min() == 0
            and perm.max() == n - 1
            and np.unique(perm).size == n
        ):
            raise ValueError(
                "PermutationOperator perm must be a permutation of "
                f"{{0, ..., {n - 1}}}; got {perm!r}."
            )
        self.perm = perm
        self.inverse_perm = np.argsort(perm)
        self.axis = int(axis)
        self.n = n
        self._domain = checked_space_extent(
            domain, n, axis=self.axis,
            owner="PermutationOperator", role="domain",
        )
        self._codomain = checked_space_extent(
            codomain, n, axis=self.axis,
            owner="PermutationOperator", role="codomain",
        )

    @property
    def domain(self) -> Optional[FunctionSpace]:
        return self._domain

    @property
    def codomain(self) -> Optional[FunctionSpace]:
        return self._codomain

    def apply(self, x: np.ndarray) -> np.ndarray:
        return np.take(x, self.perm, axis=self.axis)

    def apply_transpose(self, x: np.ndarray) -> np.ndarray:
        return np.take(x, self.inverse_perm, axis=self.axis)

    # NO ``solve``: the inverse is ALGEBRA-CLOSED —
    # :meth:`inverse` returns the inverse permutation as a first-class
    # forward whose ``apply`` is the SAME ``np.take(·, inverse_perm)``
    # gather (P^{-1} = P^T, bit-identical) — so solving is
    # ``.inverse().apply(b)``; ``apply_transpose`` keeps the gather as
    # the Euclidean-transpose verb.

    @property
    def is_invertible(self) -> bool:
        return True  # P^{-1} = P^T — a permutation is always invertible

    def inverse(self) -> "PermutationOperator":
        r"""Return :math:`P^{-1}` as a first-class :class:`PermutationOperator`.

        Built on the precomputed :attr:`inverse_perm` — its :meth:`apply`
        is the SAME integer gather as this operator's :meth:`solve` /
        :meth:`apply_transpose` (bit-identical: no arithmetic at all), and
        ``argsort`` of a permutation is exactly involutive in integer math,
        so :math:`(P^{-1})^{-1}` reproduces :attr:`perm` EXACTLY (§13 I2).

        ⭐ **The binding is INVERTED, not carried and not dropped.**
        :math:`P : V \to W` has :math:`P^{-1} : W \to V`, so the spaces
        swap. Dropping them instead — which is what this returned before
        the ends could be bound (G6.3 step 5) — would be the quieter bug of
        the two: the inverse would compose with anything, and the
        deck transformation's return leg would lose exactly the typing the
        forward leg had just been given.
        """
        return PermutationOperator(
            self.inverse_perm, axis=self.axis,
            domain=self._codomain, codomain=self._domain,
        )

    @property
    def is_adjointable(self) -> bool:
        return True


def checked_space_extent(
    space: Optional[FunctionSpace],
    expected: int,
    *,
    axis: int = 0,
    owner: str,
    role: str,
) -> Optional[FunctionSpace]:
    r"""Refuse a bound space whose extent contradicts a stored length.

    ⭐ **Why this exists as a shared primitive.** An operator that carries BOTH
    a length (``n_total``, ``n``, ``n_inflow``, ``len(cos_w)``) and a bound
    space is describing the SAME fact twice, and the two can disagree. That
    disagreement
    is **invisible at apply-time** — the arrays still broadcast, so the operator
    computes a plausible wrong answer — which makes construction the only place
    it can be caught.

    The redundancy is transitional: the tree-wide mandate (**#330**) retires
    the lengths in favour of MANDATORY spaces — G6.5 (2026-08-07) measured why
    the retirement cannot land sooner (binding is optional this era, and the
    space-less arms must still name their extents). Until then this keeps the
    pair honest, and it is one routine rather than four because four operators
    now need it
    (:class:`TraceRestrictionOperator`, :class:`PermutationOperator`,
    :class:`~orpheus.sn.boundary.angular.PartialCurrentOperator`,
    :class:`~orpheus.sn.boundary.angular.IsotropicEmissionOperator`) — well
    past the threshold at which a repeated check earns extraction.

    Parameters
    ----------
    space : FunctionSpace or None
        The bound space, or ``None`` while binding remains optional.
    expected : int
        The length the operator stores for that end.
    axis : int
        Which of the space's axes the length describes.
    owner, role : str
        Names for the error message — the class and which end (``"domain"`` /
        ``"codomain"``).

    Returns
    -------
    FunctionSpace or None
        ``space`` unchanged when consistent, so callers can assign the result.
    """
    if space is None:
        return None
    if not 0 <= axis < len(space.shape):
        raise ValueError(
            f"{owner} {role}={space!r} has {len(space.shape)} axes, so "
            f"axis={axis} is out of range."
        )
    actual = space.shape[axis]
    if actual != expected:
        raise ValueError(
            f"{owner} {role}={space!r} has extent {actual} along axis {axis}, "
            f"but this operator's {role} is {expected} rows. The space and the "
            f"length describe the SAME fact and disagree; a mis-bound space is "
            f"SILENT at apply-time because the arrays still broadcast."
        )
    return space


class TraceRestrictionOperator(LinearOperator):
    r"""Restriction onto an index subset along an axis:
    :math:`(\gamma_S x)_i = x_{S(i)}`.

    Given a **sorted, unique** index set
    :math:`S \subset \{0, \ldots, N-1\}` of size :math:`m < N`, the apply
    action gathers those entries along ``axis``, producing an array whose
    length along that axis is :math:`m` rather than :math:`N`:

    .. math::

        \gamma_S : \mathbb{R}^N \to \mathbb{R}^m, \qquad
        (\gamma_S\,x)_j \;=\; x_{S(j)}.

    The transpose is the **scatter** — zeros everywhere, the input written
    back into the selected rows:

    .. math::

        \iota_S = \gamma_S^{\mathsf T} : \mathbb{R}^m \to \mathbb{R}^N,
        \qquad
        (\iota_S\,y)_i \;=\;
        \begin{cases} y_j & i = S(j) \\ 0 & i \notin S. \end{cases}

    Relation to :class:`PermutationOperator` — a **sibling, not a
    subclass**. The mechanism is the same ``np.take`` gather with a
    non-square index array, but the algebra is different in kind: a
    permutation is a bijection (invertible, with an algebra-closed
    :meth:`~PermutationOperator.inverse`), whereas a
    restriction is **rank-deficient by construction** and has a scatter
    transpose rather than an inverse. Inheriting would have promised
    guarantees this type cannot honour.

    Why it exists
    -------------

    These are the **trace operators** :math:`\gamma_\pm` of the affine
    boundary form :math:`\gamma_-\psi = R\,G\,\gamma_+\psi + q` — the very
    maps the theory page names, which the codebase had spelled three
    different ways and typed as none of them. Every one of those spellings
    is a composition of this operator and its transpose:

    ================================================  ==================
    spelling found in the boundary subsystem          is
    ================================================  ==================
    a slice-write ``out[sel] = full[sel]``             :math:`\iota_S \circ \gamma_S`
    a dense diagonal multiply by an inflow mask        :math:`\iota_S \circ \gamma_S`
    a sparse tensor zeroing the inflow rows            :math:`I - \iota_S \circ \gamma_S`
    ================================================  ==================

    (That third spelling was a class of its own, ``IncomingOrdinateMaskTensor``,
    until campaign phase **B3.3** retired it: once B3.2 narrowed the SN law's
    domain to :math:`\Gamma_+`, the rows it preserved left the operator's
    domain and it lost its last construction site. Note it is
    :math:`I - \iota_S \circ \gamma_S`, whose range is
    :math:`\Gamma_+ \oplus \Gamma_{\text{tan}}` — the face partition is
    THREE-way, so "not inflow" was never "outflow".)

    and the observation that ``P_in ∘ P_out = 0`` stops being a curiosity:
    it is :math:`\gamma_- \circ \iota_+ = 0`, true because two disjoint
    index sets have nothing to hand each other.

    Guards, and what each one prevents
    ----------------------------------

    **Sorted** is canonical form, and it is load-bearing one tier up: the
    local↔global remap (``to_local``) lives on the half-trace SPACE since
    G6.5 — the embedding data is the space's — and its ``searchsorted``
    haystack carries, on the canonical trace-builder path, the same VALUES
    as this gather's index array (two objects; the agreement is CHECKED
    elementwise at binding — ERR-077 — not assumed). One ascending
    spelling per subspace means the gather, the space's row order, and the
    metric restriction are three views of one selection that cannot drift.
    (The 1-D-prefix-vs-2-D lesson that used to justify the remap living
    here travels with it — see
    :meth:`~orpheus.numerics.spaces.angular_trace_space.AngularFaceTraceSpace.to_local`.)

    **Unique** because a repeated row is not a restriction — the transpose
    would silently drop all but one contribution, and the pair would stop
    satisfying :math:`\gamma_S \circ \iota_S = I`.

    **The shape guard on apply** is the one a hand-rolled reduced
    permutation does not have: fed a full-length input it returns a
    same-shaped array of *wrong values*, with no raise. Measured. Here a
    domain mismatch is a loud ``ValueError`` naming both lengths.

    Structurally NON-invertible — the restriction discards rows — so it
    declares no ``inverse()``; misuse is a static error rather than a
    silent wrong answer. (A restriction covering *every* row is the
    identity in disguise; use :class:`~orpheus.numerics.operator.IdentityOperator` for that.)

    Parameters
    ----------
    indices : np.ndarray
        Sorted, unique row indices in ``[0, n_total)`` — the subspace this
        operator restricts onto.
    n_total : int
        Length of the FULL axis, i.e. the domain's extent along ``axis``.
    axis : int
        Axis to gather along. Trailing axes broadcast, so this composes as
        a :class:`~orpheus.numerics.operator.TensorProductOperator` factor.
    domain, codomain : FunctionSpace, optional
        The spaces this restriction maps between — for a boundary trace,
        :math:`\gamma_\pm : \Gamma(f) \to \Gamma_\pm(f)`. Binding them is what
        makes ``.H`` the **Hilbert** adjoint :math:`G_V^{-1}\gamma^{\mathsf T}G_W`
        rather than the bare Euclidean transpose, and what lets ``@`` refuse a
        mis-composed chain; unbound, both metric applications are silently
        skipped (campaign step G6.3, issue **#330**).

        ⭐ **Checked against** ``n_total`` **and** ``indices`` **at
        construction**, so the space and the length cannot disagree. That
        matters because they are redundant *today*: ``n_total`` /
        :attr:`n_restricted` are lengths standing where a space belongs, and
        the guard is what keeps the duplication honest until the tree-wide
        mandate (#330) retires them in favour of MANDATORY spaces. (G6.5,
        2026-08-07, measured why the retirement cannot land sooner: binding
        is optional this era, and the trace-less arm must still name its
        extent. What G6.5 did move is the local↔global remap — ``to_local``
        now lives on the half-trace space, where the index data does.) A
        binding that contradicts its own lengths is the one failure this
        class cannot detect any other way — the shapes still broadcast, so
        the wrong answer would be silent.
    """

    def __init__(
        self, indices: np.ndarray, n_total: int, axis: int = 0,
        *,
        domain: Optional[FunctionSpace] = None,
        codomain: Optional[FunctionSpace] = None,
    ) -> None:
        idx = np.asarray(indices, dtype=np.intp)
        if idx.ndim != 1:
            raise ValueError(
                f"TraceRestrictionOperator indices must be 1-D; got shape "
                f"{idx.shape}."
            )
        n_total = int(n_total)
        if n_total <= 0:
            raise ValueError(
                f"TraceRestrictionOperator n_total must be positive; got "
                f"{n_total}."
            )
        if idx.size and (idx.min() < 0 or idx.max() >= n_total):
            raise ValueError(
                f"TraceRestrictionOperator indices must lie in "
                f"[0, {n_total}); got min={idx.min()}, max={idx.max()}."
            )
        if np.unique(idx).size != idx.size:
            raise ValueError(
                "TraceRestrictionOperator indices must be unique — a repeated "
                "row is not a restriction: the scatter transpose would drop "
                f"all but one contribution. Got {idx!r}."
            )
        if idx.size > 1 and not np.all(idx[1:] > idx[:-1]):
            raise ValueError(
                "TraceRestrictionOperator indices must be SORTED ascending — "
                "one canonical spelling per subspace, matching the half-trace "
                "space's own row order (whose `to_local` searchsorts the same "
                f"array). Got {idx!r}; pass `np.sort(indices)`."
            )
        self.indices = idx
        self.n_total = n_total
        self.axis = int(axis)
        self._domain = self._checked_space(domain, n_total, "domain")
        self._codomain = self._checked_space(codomain, idx.size, "codomain")

    def _checked_space(
        self, space: Optional[FunctionSpace], expected: int, role: str,
    ) -> Optional[FunctionSpace]:
        """Refuse a space whose extent along :attr:`axis` contradicts the length.

        For the CODOMAIN, extent alone is not enough (**ERR-077**): the gather
        emits rows in :attr:`indices` order, and a space that declares its own
        ambient rows (``ordinate_indices``) must declare the SAME ones — or
        the space's ``to_local`` answers a row order this gather never
        produces. Pre-G6.5 the operator-side ``to_local`` was implicitly
        closed under the gather's own array; relocating the remap to the
        space made the agreement a separate fact, so it is CHECKED where the
        pair is bound rather than assumed. Duck-typed deliberately: numerics
        cannot import the spaces module (circular), and any space declaring
        ambient rows owes the agreement regardless of its class.
        """
        space = checked_space_extent(
            space, expected,
            axis=self.axis, owner="TraceRestrictionOperator", role=role,
        )
        declared = getattr(space, "ordinate_indices", None)
        if role == "codomain" and declared is not None:
            declared = np.asarray(declared, dtype=np.intp)
            if not np.array_equal(declared, self.indices):
                raise ValueError(
                    f"TraceRestrictionOperator: the bound codomain "
                    f"{getattr(space, 'name', space)!r} declares "
                    f"ordinate_indices {declared.tolist()} but this gather "
                    f"reads rows {self.indices.tolist()} — same extent, "
                    f"different rows. The space's to_local would answer a row "
                    f"order the gather does not emit (ERR-077)."
                )
        return space

    @property
    def domain(self) -> Optional[FunctionSpace]:
        r""":math:`\Gamma(f)` — the full space restricted FROM, when bound."""
        return self._domain

    @property
    def codomain(self) -> Optional[FunctionSpace]:
        r""":math:`\Gamma_\pm(f)` — the subspace restricted ONTO, when bound."""
        return self._codomain

    @property
    def n_restricted(self) -> int:
        """Extent of the codomain along :attr:`axis`."""
        return int(self.indices.size)

    def apply(self, x: np.ndarray) -> np.ndarray:
        x = np.asarray(x)
        if x.shape[self.axis] != self.n_total:
            raise ValueError(
                f"TraceRestrictionOperator.apply: input has "
                f"{x.shape[self.axis]} rows along axis {self.axis}, but this "
                f"restriction's DOMAIN is the full space of {self.n_total}. "
                f"Passing the restricted space back in ("
                f"{self.n_restricted} rows) is the likely mistake — the "
                f"restriction consumes the full trace and emits the subspace, "
                f"not the reverse."
            )
        return np.take(x, self.indices, axis=self.axis)

    def apply_transpose(self, x: np.ndarray) -> np.ndarray:
        r"""The scatter :math:`\iota_S` — zeros, then write the selected rows."""
        x = np.asarray(x)
        if x.shape[self.axis] != self.n_restricted:
            raise ValueError(
                f"TraceRestrictionOperator.apply_transpose: input has "
                f"{x.shape[self.axis]} rows along axis {self.axis}, but the "
                f"scatter's DOMAIN is the restricted space of "
                f"{self.n_restricted}. Passing the full space in "
                f"({self.n_total} rows) is the likely mistake."
            )
        full_shape = list(x.shape)
        full_shape[self.axis] = self.n_total
        out = np.zeros(full_shape, dtype=x.dtype)
        sel: list = [slice(None)] * x.ndim
        sel[self.axis] = self.indices
        out[tuple(sel)] = x
        return out

    # NO ``to_local``: the local↔global remap moved to the half-trace SPACE
    # at G6.5 (`AngularFaceTraceSpace.to_local`) — the embedding data is the
    # space's, and on the canonical builder path the space's row order and
    # this operator's gather are the same array by construction.

    # NO ``inverse()``: a restriction is rank-deficient by construction.
    # The transpose is the scatter, NOT an inverse — ``ι ∘ γ`` is the
    # projector onto the subspace, not the identity on the full space.

    @property
    def is_invertible(self) -> bool:
        return False

    @property
    def is_adjointable(self) -> bool:
        return True


class DiagonalOperator(PointwiseOperator):
    r"""Diagonal (pointwise) multiplication by a coefficient field.

    The operator multiplies a carrier tensor :math:`x` by a coefficient
    array :math:`c` that occupies a **sub-product** of the carrier's
    axes and is **constant** over the complementary ``broadcast_axes``:

    .. math::

        (D x)_{\mathbf{i}} \;=\;
        c_{\,\mathbf{i}\setminus\mathrm{bcast}} \; x_{\mathbf{i}}

    i.e. ``D.apply(x) == np.expand_dims(c, broadcast_axes) * x``. The
    coefficient's rank equals ``x.ndim - len(broadcast_axes)`` and its
    axes map, in order, onto the carrier axes NOT in ``broadcast_axes``.

    This is the canonical "diagonal in some basis" / pointwise-multiply
    operator. Two regimes it must express:

    * **1-D special case** — a 1-D coefficient on ONE carrier axis,
      broadcast over all others. This is the Grand Report v3 §9
      :math:`W` (``AngularWeightMatrix``) and the "multiply-by-weights
      along one axis" primitive (MoC track-weight diagonal, CP
      region-volume weighting, MC importance weighting). Spell it
      ``DiagonalOperator(w, axis=k)``; the action is rank-agnostic
      (broadcasts over however many other axes the carrier has).
    * **The multigroup-collision case** — a coefficient of shape
      ``(ng, *spatial)`` broadcast over the LEADING ordinate axis of a
      ``(N, ng, *spatial)`` angular flux, so that
      ``D.apply(psi) == sigma[None] * psi``. This is the broadcast
      engine the transport-layer ``MultiplicationOperator(σ_t)``
      delegates to. Spell it
      ``DiagonalOperator(sigma, broadcast_axes=(0,))``.

    Self-adjoint by construction (real-valued coefficient), so
    :meth:`apply_transpose` is the same code path as :meth:`apply`.
    Invertible iff every coefficient entry is non-zero — the
    VALUE-dependent arm of the split: :meth:`inverse` is declared and
    :meth:`solve` divides, both refusing eagerly
    (:class:`~orpheus.numerics.operator.NotInvertible`) on a zero entry.

    Parameters
    ----------
    coefficient : np.ndarray
        The coefficient field :math:`c`. Its rank determines how many
        carrier axes it occupies (``x.ndim - len(broadcast_axes)``).
    broadcast_axes : tuple of int, optional
        The carrier axes over which the coefficient is constant — the
        positions :func:`numpy.expand_dims` inserts singleton dims.
        When omitted, the **1-D special case** applies: ``coefficient``
        MUST be 1-D and ``axis`` selects the single carrier axis it
        occupies (rank-agnostic broadcast over the rest).
    axis : int, default 0
        Used ONLY in the 1-D special case (``broadcast_axes is None``):
        the single carrier axis the 1-D coefficient occupies. Ignored
        when ``broadcast_axes`` is given.

    Notes
    -----

    Construction does NOT materialise a dense diagonal matrix; the
    action is a single broadcast-multiply
    (``self._broadcast(x.ndim) * x``) so memory cost is
    :math:`O(\mathrm{size}(c))` regardless of the carrier's shape.

    The two construction modes are unified through one broadcast helper
    (:meth:`_broadcast`): the 1-D ``axis`` mode is rank-agnostic (the
    carrier's rank is read at apply-time, so the same operator acts on
    a 1-D, 2-D, or N-D carrier), whereas the explicit ``broadcast_axes``
    mode pins both the coefficient rank and the complement, which a
    multi-axis coefficient requires.

    The ``weights``/``axis`` attributes remain available in the 1-D
    case (``coefficient`` is exposed for both modes) as the back-compat
    alias for ``from_measure`` ergonomics and the existing 1-D call
    sites. Composition (``&`` / :class:`~orpheus.numerics.operator.TensorProductOperator`,
    ``@`` / :class:`~orpheus.numerics.operator.OperatorProduct`) does NOT read them — every composer
    routes purely through ``apply`` / ``solve``.

    Use :meth:`from_measure` when a 1-D coefficient lives on a
    :class:`~orpheus.numerics.measure.DiscreteMeasure` — common for the
    angular axis of an SN field, where the operator is built from
    ``quad.weights``.
    """

    def __init__(
        self,
        coefficient: np.ndarray,
        broadcast_axes: tuple[int, ...] | None = None,
        *,
        axis: int = 0,
    ) -> None:
        coeff = np.asarray(coefficient, dtype=float)

        if broadcast_axes is None:
            # 1-D special case: a single-axis coefficient broadcasting,
            # rank-agnostically, over every other carrier axis. The
            # carrier rank is unknown at construction, so the broadcast
            # placement is deferred to apply-time (see _broadcast).
            if coeff.ndim != 1:
                raise ValueError(
                    f"DiagonalOperator without broadcast_axes is the 1-D "
                    f"special case; coefficient must be 1-D, got shape "
                    f"{coeff.shape}. For an N-D coefficient pass "
                    f"broadcast_axes=(...)."
                )
            self.broadcast_axes: tuple[int, ...] | None = None
        else:
            # General case: an N-D coefficient pinned onto an explicit
            # complement of carrier axes.
            bcast = tuple(int(a) for a in broadcast_axes)
            if len(set(bcast)) != len(bcast):
                raise ValueError(
                    f"DiagonalOperator broadcast_axes must be distinct, "
                    f"got {bcast}."
                )
            self.broadcast_axes = bcast

        # ``axis`` is consulted ONLY in the 1-D special case; storing it
        # as a plain int (default 0) keeps the attribute well-typed and
        # harmless in broadcast mode.
        self.axis = int(axis)
        self.coefficient = coeff

    @classmethod
    def from_measure(
        cls, measure, axis: int = 0,
    ) -> "DiagonalOperator":
        """Construct from the weights of a :class:`DiscreteMeasure`.

        Convenience constructor for the canonical 1-D case where the
        diagonal IS the discrete measure's weights — e.g.
        ``DiagonalOperator.from_measure(quad.measure, axis=0)`` is the
        Grand Report v3 §9 ``AngularWeightMatrix``.
        """
        return cls(measure.weights, axis=axis)

    @property
    def weights(self) -> np.ndarray:
        """The 1-D coefficient vector (the historical ``weights`` name).

        Available ONLY in the 1-D special case (``broadcast_axes is
        None``); it is the back-compat alias for ``from_measure`` and
        the existing 1-D call sites. Reading it on a multi-axis-
        coefficient instance is an illegal state and raises (Pattern 4)
        rather than returning an N-D array under a 1-D name.
        """
        if self.broadcast_axes is not None:
            raise AttributeError(
                "DiagonalOperator.weights is the 1-D special case's "
                "coefficient vector; this operator has an N-D coefficient "
                f"of shape {self.coefficient.shape} on broadcast_axes="
                f"{self.broadcast_axes}. Use .coefficient instead."
            )
        return self.coefficient

    def _broadcast(self, ndim: int) -> np.ndarray:
        """Reshape the coefficient to broadcast over an ``ndim`` carrier.

        Single source of truth for both construction modes:

        * 1-D ``axis`` mode — return a view of shape
          ``(1, ..., 1, N, 1, ..., 1)`` with ``N`` at ``self.axis``
          (rank-agnostic: built fresh for the carrier's actual ``ndim``).
        * explicit ``broadcast_axes`` mode — return
          ``np.expand_dims(coefficient, broadcast_axes)``, inserting a
          singleton at each broadcast axis so the coefficient occupies
          the complementary axes in order.
        """
        if self.broadcast_axes is None:
            shape = [1] * ndim
            shape[self.axis] = -1
            return self.coefficient.reshape(shape)
        return np.expand_dims(self.coefficient, self.broadcast_axes)

    def _check_shape(self, x: np.ndarray) -> None:
        """Validate the carrier's rank/axis sizes against the coefficient."""
        if self.broadcast_axes is None:
            if x.shape[self.axis] != self.coefficient.shape[0]:
                raise ValueError(
                    f"DiagonalOperator(axis={self.axis}) expects axis size "
                    f"{self.coefficient.shape[0]}; got {x.shape[self.axis]} "
                    f"in input of shape {x.shape}."
                )
            return
        expected_rank = x.ndim - len(self.broadcast_axes)
        if self.coefficient.ndim != expected_rank:
            raise ValueError(
                f"DiagonalOperator(broadcast_axes={self.broadcast_axes}) "
                f"expects a rank-{expected_rank} coefficient for a "
                f"{x.ndim}-D carrier; got rank-{self.coefficient.ndim} "
                f"coefficient of shape {self.coefficient.shape}."
            )

    def apply(self, x: np.ndarray) -> np.ndarray:
        x_arr = np.asarray(x)
        self._check_shape(x_arr)
        return self._broadcast(x_arr.ndim) * x_arr

    def apply_transpose(self, x: np.ndarray) -> np.ndarray:
        # Real-valued diagonal is self-adjoint.
        return self.apply(x)

    def solve(self, b_vec: np.ndarray) -> np.ndarray:
        # The value-dependent guard (Pattern 4 heritage: never a silent
        # IEEE NaN on a σ=0 division — the legacy bare-σ path had no gate).
        if not self.is_invertible:
            raise NotInvertible(
                "DiagonalOperator.solve requires non-zero coefficient "
                "entries; this operator has at least one zero entry."
            )
        b_arr = np.asarray(b_vec)
        self._check_shape(b_arr)
        return b_arr / self._broadcast(b_arr.ndim)

    def solve_transpose(self, b_vec: np.ndarray) -> np.ndarray:
        # A real diagonal is its own transpose, so (D^T)^{-1} = D^{-1}.
        return self.solve(b_vec)

    @property
    def is_invertible(self) -> bool:
        # Invertible iff every coefficient entry is non-zero (D^{-1} = 1/c).
        return bool(np.all(self.coefficient != 0.0))

    @property
    def is_transpose_invertible(self) -> bool:
        # D^T = D, so the transposed solve exists exactly when the solve does.
        return self.is_invertible

    def inverse(self) -> "InverseOperator":
        r"""Return :math:`D^{-1}` as an :class:`~orpheus.numerics.operator.InverseOperator` over this leaf.

        Delegation, NOT a reciprocal-coefficient twin: the returned
        object's ``apply`` IS :meth:`solve` (the division ``b / c``),
        bit-identical — whereas ``DiagonalOperator(1/c)`` would multiply
        by a rounded reciprocal and drift by an ulp. The generic name is
        the honest one (round-trip alone earns exactly
        "InverseOperator"; a diagonal division carries no distinguishing
        invariant beyond it).
        """
        if not self.is_invertible:
            raise NotInvertible(
                "DiagonalOperator.inverse requires non-zero coefficient "
                "entries; this operator has at least one zero entry."
            )
        return InverseOperator(self)



class RankOneOperator(LinearOperator):
    r"""The rank-1 dyad :math:`|v\rangle\langle w|` — a reconstruction column ⊗ a functional row.

    A :class:`RankOneOperator` is the outer product of a **reconstruction**
    vector :math:`v` (the column, an output direction) and a
    :class:`~orpheus.numerics.functional.Functional` :math:`\langle w|` (the
    row, the contraction). Its action on a carrier :math:`x` is

    .. math::

        (\,|v\rangle\langle w|\,)\,x \;=\; v \,\cdot\, \langle w, x\rangle ,

    i.e. ``reconstruction * functional.evaluate(x)``: the functional contracts
    :math:`x` to the inner product :math:`\langle w, x\rangle` (with
    ``keepdims`` on the contracted axis), and the reconstruction broadcasts back
    over that length-1 axis. The functional OWNS the contraction (its weight and
    axis); the operator only broadcasts — so the matvec routes THROUGH
    ``functional.evaluate`` and there is no parallel inline reduction to drift
    from it.

    Build instances with :func:`outer` (the readable verb,
    ``outer(reconstruction, functional)``). A genuine ``M × K`` rank-1 operator
    (``v ∈ ℝ^M``, ``w ∈ ℝ^K``, ``M ≠ K``) is legal — there is no same-shape
    constraint between the column and the row (the old ``left.shape ==
    right.shape`` check was an artifact of the square-only legacy form).

    Native to the multigroup fission emission
    :math:`F = |\chi\rangle\langle\nu\Sigma_f| =
    \texttt{outer}(\chi,\ \mathrm{ReactionRateFunctional}(\nu\Sigma_f))`
    (Grand Report v3 §15): the production-rate co-vector
    :math:`\langle\nu\Sigma_f, \phi\rangle = \sum_{g'}\nu\Sigma_{f,g'}\phi_{g'}`
    is the
    :class:`~orpheus.transport.reaction_rate_functional.ReactionRateFunctional`
    row, and the emission spectrum :math:`\chi` is the reconstruction column.
    Fission is the **rank-1 (single-mode) degenerate** of the multi-mode
    scattering kernel ``R∘Λ∘M`` (a :class:`~orpheus.numerics.frame.FrameBase`
    manages the analogous *stack* of dyads); see
    :mod:`orpheus.transport.operators.integral_kernel_operator`.

    Relation to :class:`~orpheus.numerics.operator.TensorProductOperator`
    -------------------------------------------

    A :class:`RankOneOperator` satisfies the TP-factor contract (it acts on the
    functional's contracted axis and broadcasts on the others), so it composes
    as a TP factor when the algebra wants the type-visible separable form:

    .. code-block:: python

        F_kernel = outer(chi, InnerProductFunctional(nu_sigma_f, axis=0)) & IdentityOperator()

    The :class:`~orpheus.numerics.operator.IdentityOperator` factor advertises the spatial-axis broadcast;
    the TP fold reduces to :meth:`RankOneOperator.apply` bit-identically
    (``IdentityOperator.apply`` returns ``x``).

    Adjointable exactly when the row is an
    :class:`~orpheus.numerics.functional.InnerProductFunctional` (the usual
    case, including its
    :class:`~orpheus.transport.reaction_rate_functional.ReactionRateFunctional`
    specialisation) — the VALUE-dependent arm on the adjoint axis. Rank-1
    operators are **structurally non-invertible** (no ``inverse()`` declared —
    the kernel is the orthogonal complement of the row along the contracted axis),
    but they DO have a **transpose**: :meth:`apply_transpose` is the dual dyad
    :math:`|w\rangle\langle v|` — swap the column :math:`v` with the row's
    weight :math:`w`, contracting :math:`\langle v,\cdot\rangle` over the same
    axis. This is the Euclidean transpose :math:`A^{T}`; the metric-correct
    Hilbert adjoint :math:`A^\dagger = G^{-1}A^{T}G` is the
    :attr:`~LinearOperator.H` wrapper's job. The fission adjoint
    :math:`F^\dagger\psi^* = \nu\Sigma_f\,(\chi\cdot\psi^*)` is exactly this
    dyad-swap (#276). A nonlinear / opaque functional has no dual
    column, so such a dyad advertises ``apply`` only. See
    :ref:`operator-adjoint`.

    Parameters
    ----------
    reconstruction : Vector | numpy.ndarray
        The column :math:`v` — the output direction the inner product is
        broadcast against. Aligns with the carrier on the complement of the
        functional's contracted axis; its size on that axis is the output
        dimension ``M``.
    functional : Functional
        The row co-vector :math:`\langle w|` — contracts the carrier to
        :math:`\langle w, x\rangle` over its own axis (typically the leading
        group axis for the multigroup reaction rate). Usually an
        :class:`~orpheus.numerics.functional.InnerProductFunctional` (generic)
        or a
        :class:`~orpheus.transport.reaction_rate_functional.ReactionRateFunctional`
        (the domain-typed reaction rate).

    Notes
    -----
    Bit-identity with the legacy ``(right * x).sum(axis, keepdims) * left``
    formulation is preserved because
    :meth:`~orpheus.numerics.functional.InnerProductFunctional.evaluate` performs
    that exact ``(w * x).sum(axis, keepdims=True)`` reduction — the same numpy
    primitive, the same axis, the same order — and the reconstruction broadcast
    is elementwise. IEEE-754 pairwise-reduction order is preserved.
    """

    def __init__(
        self,
        reconstruction: "Vector | np.ndarray",
        functional: "Functional",
    ) -> None:
        # The dyad |v⟩⟨w|: ``reconstruction`` is the column (output) vector v;
        # ``functional`` is the row co-vector ⟨w| that OWNS the contraction (its
        # weight and axis). NO same-shape guard — a genuine M×K rank-1 operator
        # (M ≠ K) is legal; the functional validates its own contraction against
        # the carrier at apply time (the old left.shape == right.shape check was
        # an artifact of the square-only legacy form, not a real constraint).
        self.reconstruction = reconstruction
        self.functional = functional

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        r"""``None`` — a DOCUMENTED Optional (the S4-amendment's fourth
        answer): the dyad's two vectors are bare arrays that carry no
        spaces today, so there is nothing honest to derive. The owning
        flip is CS4c's kernel binding (``FissionKernel``'s dyad bound to
        its Energy space), where the row/column gain their spaces and
        this override narrows.
        """
        return None

    @property
    def codomain(self) -> Optional["FunctionSpace"]:
        r"""``None`` — see :attr:`domain` (the CS4c-owned debt)."""
        return None

    @property
    def is_adjointable(self) -> bool:
        # The VALUE-dependent arm on the ADJOINT axis: the dual dyad
        # |w⟩⟨v| exists iff the row ⟨w| is a genuine co-vector whose
        # weight IS the dual column — an InnerProductFunctional (the
        # ReactionRateFunctional specialisation included). A nonlinear /
        # opaque functional has no dual column. Mirrors — and gates —
        # the apply_transpose realization below. Whether ``.H`` is POSED is a
        # separate question (:func:`~orpheus.numerics.operator.posed`): this dyad's ends are ``None``
        # (the CS4c debt), so ``.H`` refuses while ``.dual()`` runs.
        from orpheus.numerics.functional import InnerProductFunctional

        return isinstance(self.functional, InnerProductFunctional)

    def apply(self, x: np.ndarray) -> np.ndarray:
        # |v⟩⟨w| x = v · ⟨w, x⟩. The functional IS the contraction — it returns
        # the inner product ⟨w, x⟩ with ``keepdims`` on the contracted axis, and
        # the reconstruction broadcasts back over that length-1 axis. Routing the
        # matvec THROUGH ``functional.evaluate`` is what makes the row-factor a
        # first-class object (no parallel inline reduction to drift from it).
        recon = np.asarray(getattr(self.reconstruction, "values", self.reconstruction))
        return recon * self.functional.evaluate(x)

    def apply_transpose(self, x: np.ndarray) -> np.ndarray:
        # The dual dyad: (|v⟩⟨w|)ᵀ = |w⟩⟨v|. The transpose swaps the column and
        # the row — the new column is the old row's weight w (the dual column),
        # the new row is ⟨v| (the old reconstruction as a co-vector on the SAME
        # contracted axis). So Aᵀx = w · ⟨v, x⟩, the Euclidean transpose (the .H
        # wrapper adds the metric). Reuses the InnerProductFunctional contraction
        # primitive — single source of truth with the forward `apply` row-factor.
        from orpheus.numerics.functional import InnerProductFunctional

        # The isinstance IS the narrowing (the same fact is_adjointable
        # advertises): the body reads the IPF-typed row's .weight/.axis.
        if not isinstance(self.functional, InnerProductFunctional):
            raise MissingAdjoint(
                "RankOneOperator.apply_transpose requires the row functional to "
                "be an InnerProductFunctional (a co-vector with a dual column); "
                f"got {type(self.functional).__name__} — a nonlinear functional "
                "has no dual column."
            )
        column = np.asarray(
            getattr(self.functional.weight, "values", self.functional.weight)
        )
        dual_row = InnerProductFunctional(
            np.asarray(getattr(self.reconstruction, "values", self.reconstruction)),
            axis=self.functional.axis,
        )
        return column * dual_row.evaluate(x)


def outer(
    reconstruction: "Vector | np.ndarray",
    functional: "Functional",
) -> RankOneOperator:
    r"""Build the rank-1 dyad :math:`|v\rangle\langle w|` from a column and a co-vector.

    The universal constructor for a rank-1 :class:`~orpheus.numerics.operator.LinearOperator`: a
    :class:`~orpheus.numerics.vector.Vector` (or ``ndarray``) ``reconstruction``
    :math:`v` (the column, the output direction) tensored with a
    :class:`~orpheus.numerics.functional.Functional` ``functional`` :math:`\langle w|`
    (the row, the contraction). The action is

    .. math::

        (\,|v\rangle\langle w|\,)\,x \;=\; v \,\cdot\, \langle w, x\rangle ,

    i.e. ``reconstruction * functional.evaluate(x)``. Every separable rank-1
    operator in the algebra is one of these; the multi-mode generalisation is a
    sum of dyads managed by a :class:`~orpheus.numerics.frame.FrameBase` (the
    spectral / block-term decomposition). The canonical transport instance is
    the fission emission kernel
    :math:`F = \texttt{outer}(\chi,\ \mathrm{ReactionRateFunctional}(\nu\Sigma_f))`
    (see :class:`~orpheus.transport.operators.fission.FissionOperator`).
    """
    return RankOneOperator(reconstruction, functional)
