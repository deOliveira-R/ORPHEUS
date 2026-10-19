r"""The marginal over one axis of a product space: retraction, section, pullback.

A space built from axes, :math:`V = V_k \otimes V_{\rm rest}`, carries the
product metric :math:`G = G_k \otimes G_{\rm rest}` (an axis-built space
admits no other: its metric is derived from its axes, and a positioned form
sits only on a measure-less axis's own block). Marginalising over the axis
:math:`k` is integration along it, which is the pairing with the constant
function :math:`1` on that factor.

**The unit and its covector** (ruling 36 of the numerics-untangling plan).
The axis owns its unit :math:`u`, the coordinates of the constant function
:math:`1` in the axis's own representation
(:attr:`~orpheus.numerics.axis.Axis.unit`; minted with the axis, ruling 37):
all ones on a nodal axis, :math:`e_0` on a continuum harmonic or Legendre
head, the moments :math:`\Phi^{\mathsf T} w` on a head a frame has dressed,
:math:`e_0` on the LD scheme's spatial-moment axis. The covector of
integration is :math:`m = G_k\,u`, the axis's metric block applied to the
unit (the weights, a positioned form, or the identity of the counting
measure). Then :math:`\langle m, c\rangle = \langle u, c\rangle_{G_k}` is
the integral of the function :math:`c` represents, against whatever measure
the axis carries, and :math:`\langle m, u\rangle` is that measure's MASS
(:math:`\Sigma w`, :math:`4\pi`, the reference cell's volume 1).

**The three arrows**, minted together, once per (space, axis label):

* the retraction :math:`R = m^{\mathsf T} \otimes \mathrm{id}_{\rm rest}`
  (:class:`AxisRetractionOperator`), fibre integration, the pushforward
  :math:`\pi_*`;
* the pullback :math:`\pi^* = u \otimes \mathrm{id}_{\rm rest}`
  (:class:`AxisPullbackOperator`), the Hilbert adjoint :math:`R^\dagger`
  exactly, because under the product metric
  :math:`\langle R\psi, \phi\rangle_{\rm rest}
  = \langle \psi, u \otimes \phi\rangle_G`: no metric sandwich, and
  :math:`G` is never inverted;
* the section :math:`E = (u / \langle m, u\rangle) \otimes
  \mathrm{id}_{\rm rest}` (:class:`AxisSectionOperator`), the Moore–Penrose
  right inverse :math:`R^\dagger (R R^\dagger)^{-1}`, so :math:`R \circ E =
  \mathrm{id}`. It does not exist when the mass is zero (a signed measure
  summing to zero).

The nodal case is the pair the tree shipped before the unit existed:
:math:`u = 1`, :math:`m = w`, :math:`\langle m, u\rangle = \Sigma w`, and
the three arrows are bit-identical to it (gated). The divisor is the
correctly rounded mass, ``math.fsum`` (ruling 38).

**What is refused, and why** (each at the one mint, :func:`_lift`):

* an axes-less space, or a label naming not exactly one axis (the space's
  by-label resolver, ``FunctionSpace._axis_index``);
* a typed :class:`~orpheus.numerics.axis.EnergyAxis`: the energy axis
  persists at its one-cell member, and its collapse is condensation;
* a single-axis space: its marginal would be a bare scalar, not a space;
* the section over an axis of zero mass.

Before ruling 36 a MODAL axis was refused ("the modal average is the
coefficient at the average slot") and so was a positioned form on the
collapsed axis ("two measures on one block"). Both refusals are retired: the
unit makes the modal marginal an integral, and a form on the axis's block IS
:math:`G_k`, which :math:`m = G_k u` uses.

The arrows are downstream constructors of the space (``.of(space, label)``,
ruling 34's form); the space itself carries no lift verb.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, NamedTuple

import numpy as np

from orpheus.numerics.axis import EnergyAxis
from orpheus.numerics.instance_cache import cached_on
from orpheus.numerics.operator import LinearOperator

if TYPE_CHECKING:
    from orpheus.numerics.space import FunctionSpace

__all__ = ["AxisPullbackOperator", "AxisRetractionOperator", "AxisSectionOperator"]


class _AxisMarginalBase(LinearOperator):
    r"""What the three arrows of one axis marginal retain.

    The bound product spaces, the ndarray dims the axis occupies, the axis's
    shape, and one flat vector over the axis (the covector :math:`m` on the
    retraction, the unit :math:`u` on the pullback and the section). Nothing
    else: no axis, no frame, no generator.

    An axis may span several ndarray dims (the 2-D spatial axis is one
    factor of shape ``(nx, ny)``); the arrows contract and scatter over all
    of them. Every arrow is born bound: an operator without its two spaces
    does not exist here.
    """

    def __init__(
        self,
        *,
        full_space: "FunctionSpace",
        marginal_space: "FunctionSpace",
        axis_shape: tuple[int, ...],
        dims: tuple[int, ...],
    ) -> None:
        self._full_space = full_space
        self._marginal_space = marginal_space
        self._axis_shape = axis_shape
        self._dims = dims

    # ── the two kernels (one spelling per direction) ────────────────────

    def _contract(self, x: np.ndarray, vector: np.ndarray) -> np.ndarray:
        r"""``Σ_axis vector · x`` — the pairing with a flat vector over the axis.

        Spelled to be BIT-IDENTICAL with the shipped angular reduction
        (``einsum("n,ng...->g...", w, values)``) on its case, a leading
        one-dim axis: ``moveaxis`` to the front is then the identity view
        (gated ``np.array_equal``). A multi-dim axis is fused to one flat
        dim first.
        """
        x = np.asarray(x)
        if len(self._dims) == 1:
            moved = np.moveaxis(x, self._dims[0], 0)
        else:
            moved = np.moveaxis(x, self._dims, range(len(self._dims)))
            moved = moved.reshape(-1, *moved.shape[len(self._dims):])
        return np.einsum("n,n...->...", vector, moved)

    def _scatter(self, x: np.ndarray, vector: np.ndarray) -> np.ndarray:
        r"""``vector ⊗ x`` — ``x`` placed across the axis dims, slot ``n``
        scaled by ``vector[n]`` (a fresh C-contiguous array)."""
        out = np.multiply.outer(vector.reshape(self._axis_shape), np.asarray(x))
        return np.ascontiguousarray(np.moveaxis(out, range(len(self._dims)), self._dims))

    def _require(self, x: np.ndarray, space: "FunctionSpace", which: str, verb: str) -> np.ndarray:
        x = np.asarray(x)
        if x.shape != space.shape:
            raise ValueError(
                f"{type(self).__name__}.{verb}: input shape {x.shape} does "
                f"not match the {which} space {space.shape}."
            )
        return x

    @property
    def is_invertible(self) -> bool:
        return False

    @property
    def is_adjointable(self) -> bool:
        return True


class AxisRetractionOperator(_AxisMarginalBase):
    r"""The retraction :math:`R = \pi_* = m^{\mathsf T} \otimes \mathrm{id}`:
    fibre integration over one named axis,
    :math:`(R\psi)(\cdot) = \sum_n m_n\, \psi(n, \cdot)`.

    **Canonical names.** :math:`R \circ E = \mathrm{id}` (up to the
    re-association of :math:`\Sigma_n m_n u_n / \langle m, u\rangle`:
    `[M]` at most :math:`N/3` ULP over 120 shipped rules, :math:`N` the
    axis's slot count) makes :math:`R` a *retraction* (split epimorphism)
    and :math:`E` its *section* (Mac Lane, CWM §I.5). Its Hilbert adjoint is
    the pullback :math:`\pi^*` (:class:`AxisPullbackOperator`), returned by
    :meth:`adjoint` as an object minted with the retraction, so
    ``R.H.H is R``. The generic metric sandwich :math:`\sharp\,
    R^{\mathsf T} \flat` is the same operator in exact arithmetic on a
    nonsingular metric; it re-associates and so rounds a few ULP away, and
    where the metric is singular (a zero weight, a frame head of deficient
    rank) the adjoint is not unique and the two differ while both satisfy
    reciprocity. The closed form is the law, and the sandwich is its
    independent witness.

    On the angular axis this is the scalar flux
    :math:`\phi = \sum_n w_n \psi_n` (`[M]` bit-identical with the shipped
    einsum), on a spatial axis the volume integral, and on a moment head the
    integral of the represented function: :math:`\langle m, \varphi\rangle =
    \int \psi\,d\mu` on a frame-dressed head (:math:`m = G^{+} B e_0`),
    :math:`4\pi c_{00}` on a continuum harmonic head.

    The two lifts differ by the mass: :math:`R^\dagger = \langle m,
    u\rangle\, E` up to rounding. Naming both arrows is the anti-ERR-051
    move: a single verb would have had to choose a convention, and a
    re-pointed call site would silently change a source by the mass.

    Structurally rank-deficient (the marginal discards the axis); the
    Euclidean transpose is the scatter :math:`m \otimes \phi`.
    """

    def __init__(
        self,
        *,
        full_space: "FunctionSpace",
        marginal_space: "FunctionSpace",
        axis_shape: tuple[int, ...],
        dims: tuple[int, ...],
        covector: np.ndarray,
        unit: np.ndarray,
    ) -> None:
        super().__init__(
            full_space=full_space, marginal_space=marginal_space,
            axis_shape=axis_shape, dims=dims,
        )
        self._covector = covector
        self._unit = unit
        # The adjoint is minted with the retraction from the same induced
        # data; the pair refer to each other, so ``R.H.H is R``.
        self._pullback = AxisPullbackOperator(self)

    @classmethod
    def of(cls, space: "FunctionSpace", axis_label: str) -> "AxisRetractionOperator":
        r"""The retraction over ``axis_label`` of ``space``, minted once per
        (space, label) together with its section and pullback (:func:`_lift`)."""
        return _lift(space, axis_label).retraction

    @property
    def covector(self) -> np.ndarray:
        r"""The covector of integration :math:`m = G_k u`, flat over the axis."""
        return self._covector

    def adjoint(self) -> "AxisPullbackOperator":
        r"""The Hilbert adjoint :math:`R^\dagger = \pi^*`, the pullback."""
        return self._pullback

    @property
    def domain(self) -> "FunctionSpace":
        r"""The full space, the product carrying the integrated axis."""
        return self._full_space

    @property
    def codomain(self) -> "FunctionSpace":
        r"""The marginal space: the remaining axes, measures and forms intact."""
        return self._marginal_space

    def apply(self, x: np.ndarray) -> np.ndarray:
        return self._contract(self._require(x, self._full_space, "full", "apply"), self._covector)

    def apply_transpose(self, x: np.ndarray) -> np.ndarray:
        r"""The scatter :math:`m \otimes \phi`, the Euclidean transpose."""
        return self._scatter(
            self._require(x, self._marginal_space, "marginal", "apply_transpose"), self._covector,
        )


class AxisSectionOperator(_AxisMarginalBase):
    r"""The section of the axis retraction,
    :math:`E = (u / \langle m, u\rangle) \otimes \mathrm{id}`:
    :math:`(E\phi)(n, \cdot) = u_n\, \phi(\cdot) / \langle m, u\rangle`.

    **Canonical name.** :math:`E` is DEFINED by :math:`R \circ E =
    \mathrm{id}`, the right inverse of the retraction (split monomorphism;
    Mac Lane CWM §I.5); it is the Moore–Penrose right inverse
    :math:`R^\dagger (R R^\dagger)^{-1}`. "Embedding" was rejected as
    non-canonical (ratified 2026-08-24): the pullback is an embedding too,
    so that name cannot discriminate the two arrows. The composite
    :math:`P = E \circ R` is the conditional expectation onto the functions
    constant along the axis.

    **The divisor is the mass** :math:`\langle m, u\rangle` of the axis's
    measure, summed correctly rounded (``math.fsum``, ruling 38): a named
    quantity with a platform-free value, where the retired spelling (the
    1×1 Gram entry of a throw-away indicator frame) was the least accurate
    of the candidates (`[M]` up to 39 ULP off the correctly rounded sum, on
    ``lebedev29``).

    On the angular axis this is the isotropic-source projection
    :math:`Q / \Sigma w` broadcast across the ordinates (``u = 1``; divide,
    then scatter, the shipped kernel's float operations). It is NOT the
    adjoint of :class:`AxisRetractionOperator`: that is the pullback, which
    does not divide.

    Does not exist on an axis of zero mass: :meth:`of` refuses there, while
    the retraction over the same axis stays legal.
    """

    def __init__(
        self,
        *,
        full_space: "FunctionSpace",
        marginal_space: "FunctionSpace",
        axis_shape: tuple[int, ...],
        dims: tuple[int, ...],
        unit: np.ndarray,
        total_weight: float,
    ) -> None:
        super().__init__(
            full_space=full_space, marginal_space=marginal_space,
            axis_shape=axis_shape, dims=dims,
        )
        self._unit = unit
        self._total_weight = total_weight

    @classmethod
    def of(cls, space: "FunctionSpace", axis_label: str) -> "AxisSectionOperator":
        r"""The section over ``axis_label`` of ``space``, minted once per
        (space, label) with the retraction (:func:`_lift`).

        Raises
        ------
        ValueError
            If the axis's measure has zero mass (a signed measure summing to
            zero): the section divides by it, so none exists.
        """
        section = _lift(space, axis_label).section
        if section is None:
            raise ValueError(
                f"AxisSectionOperator.of: axis {axis_label!r} of {space!r} has "
                f"zero total weight (a signed measure whose mass <m, u> is 0) "
                f"— the section divides by the mass, so none exists. The "
                f"retraction over this axis is still legal."
            )
        return section

    @property
    def total_weight(self) -> float:
        r"""The mass :math:`\langle m, u\rangle` of the axis's measure (``math.fsum``)."""
        return self._total_weight

    @property
    def domain(self) -> "FunctionSpace":
        r"""The marginal space the section lifts FROM."""
        return self._marginal_space

    @property
    def codomain(self) -> "FunctionSpace":
        r"""The full space the section lifts INTO."""
        return self._full_space

    def apply(self, x: np.ndarray) -> np.ndarray:
        x = self._require(x, self._marginal_space, "marginal", "apply")
        # Divide THEN scatter: on a nodal axis (u = 1) these are the shipped
        # isotropic-source kernel's float operations, so the equivalence
        # gate pins np.array_equal.
        return self._scatter(x / self._total_weight, self._unit)

    def apply_transpose(self, x: np.ndarray) -> np.ndarray:
        r"""The pairing with :math:`u` over the mass, the Euclidean transpose."""
        x = self._require(x, self._full_space, "full", "apply_transpose")
        return self._contract(x, self._unit) / self._total_weight


class AxisPullbackOperator(_AxisMarginalBase):
    r"""The pullback :math:`\pi^* = u \otimes \mathrm{id}` along the
    projection that forgets one axis: :math:`(\pi^*\phi)(n, \cdot) = u_n\,
    \phi(\cdot)`, the constant function along the axis times :math:`\phi`.
    The Hilbert adjoint of :class:`AxisRetractionOperator`.

    **Why it is the adjoint.** For :math:`R\psi = \sum_n m_n \psi_n` with
    :math:`m = G_k u`, on a full space whose metric is the product
    :math:`G = G_k \otimes G_{\rm rest}`,

    .. math::

       \langle R\psi, \phi\rangle_{\rm rest}
       = \sum_{n,n'} (G_k)_{n n'} u_{n'} \langle \psi_n, \phi\rangle_{\rm rest}
       = \langle \psi, u \otimes \phi\rangle_G,

    so :math:`R^\dagger = \pi^*`, with no metric in it. Every axis-built
    space carries the product metric, so the closed form is the adjoint of
    every retraction :func:`_lift` can mint, and ``.H`` carries no guard.

    **The detector lift.** On the angular axis (:math:`u = 1`) this is how a
    detector response :math:`\Sigma_d(\vec r, g)` enters phase space: the
    functional :math:`\psi \mapsto \langle \Sigma_d, R\psi\rangle` has the
    Riesz representative :math:`R^\dagger \Sigma_d = \pi^*\Sigma_d`, with no
    division by the mass. A source rate enters through the section
    :class:`AxisSectionOperator` instead. The two lifts differ by
    :math:`R \circ R^\dagger = \langle m, u\rangle`: :math:`4\pi` for a rule
    on the whole sphere, 2 for a 1-D rule.

    Minted only by the retraction it is the adjoint of; the two refer to
    each other, so ``R.H.H is R``.
    """

    def __init__(self, retraction: AxisRetractionOperator) -> None:
        if hasattr(retraction, "_pullback"):
            raise TypeError(
                "AxisPullbackOperator is minted by its retraction, once; read "
                "it as retraction.H (a second pullback would break R.H.H is R)"
            )
        super().__init__(
            full_space=retraction.domain,
            marginal_space=retraction.codomain,
            axis_shape=retraction._axis_shape,
            dims=retraction._dims,
        )
        self._retraction = retraction
        self._unit = retraction._unit

    @property
    def domain(self) -> "FunctionSpace":
        r"""The marginal space the pullback lifts FROM."""
        return self._marginal_space

    @property
    def codomain(self) -> "FunctionSpace":
        r"""The full space the pullback lifts INTO."""
        return self._full_space

    def apply(self, x: np.ndarray) -> np.ndarray:
        return self._scatter(self._require(x, self._marginal_space, "marginal", "apply"), self._unit)

    def apply_transpose(self, x: np.ndarray) -> np.ndarray:
        r"""The pairing with :math:`u`, the Euclidean transpose."""
        return self._contract(self._require(x, self._full_space, "full", "apply_transpose"), self._unit)

    def adjoint(self) -> AxisRetractionOperator:
        r"""The involution as an object identity: ``R.H.H is R``."""
        return self._retraction


class _AxisMarginal(NamedTuple):
    """One axis's marginal arrows (``section`` is ``None`` iff the mass is 0)."""

    retraction: AxisRetractionOperator
    section: AxisSectionOperator | None


def _lift(space: "FunctionSpace", axis_label: str) -> _AxisMarginal:
    r"""The marginal over ``axis_label``, minted once per (space, label).

    The one mint of the three arrows (the stage-2 discipline: what the
    operators retain is induced data copied out of the axis and the space,
    and nothing holds the generator). Memoised in the space's private slot
    ``"_axis_marginal"``, so the retraction, the section and the pullback
    share ONE marginal-space instance (``R.codomain is E.domain``) and a
    carrier that caches its spaces gets warm operators.
    """
    return cached_on(space, "_axis_marginal", axis_label, lambda: _mint(space, axis_label))


def _mint(space: "FunctionSpace", axis_label: str) -> _AxisMarginal:
    # Axis-built, exactly one label hit: the space's by-label resolver, the
    # one home of that refusal vocabulary.
    k = space._axis_index(axis_label)
    axes = space.axes
    assert axes is not None  # narrowing only: _axis_index refused the None case
    axis = axes[k]
    if isinstance(axis, EnergyAxis):
        raise TypeError(
            f"axis marginal: axis {axis_label!r} is a typed EnergyAxis, which "
            f"PERSISTS at its one-cell member (collapse doctrine clause 2 — "
            f"the one-group limit keeps its edges and spectrum, because "
            f"⟨σ̄,φ⟩ consumes the partition). The energy collapse is "
            f"condensation: use EnergyGrid.overlap_to / the Petrov-Galerkin "
            f"condensation frame, not a drop-form marginal."
        )
    if len(axes) == 1:
        raise ValueError(
            f"axis marginal: {space!r} has only the {axis_label!r} axis — its "
            f"marginal would be a bare scalar, which is not a FunctionSpace. "
            f"Contract with the space's inner product instead."
        )

    # The ndarray dims this axis occupies: axes map to dims by cumulative
    # rank (an axis's shape may span several dims).
    start = sum(len(ax.shape) for ax in axes[:k])
    dims = tuple(range(start, start + len(axis.shape)))
    # The marginal keeps the remaining axes' measures and forms, so the full
    # metric is exactly G_k ⊗ (the marginal's metric).
    marginal_space = space.without_axis(axis_label)

    # m = G_k u: the axis's metric block is a positioned form (on a
    # measure-less axis), else the axis's weights, else the counting measure.
    unit = axis.unit
    assert unit is not None  # narrowing only: the axis fills its unit at construction
    form = space._overlay_forms()[k]
    if form is not None:
        covector = np.asarray(form.apply(unit))
    elif axis.weights is not None:
        covector = axis.weights * unit
    else:
        covector = unit
    flat_unit = unit.ravel()
    flat_covector = covector.ravel()
    mass = math.fsum(flat_covector * flat_unit)

    retraction = AxisRetractionOperator(
        full_space=space, marginal_space=marginal_space,
        axis_shape=axis.shape, dims=dims,
        covector=flat_covector, unit=flat_unit,
    )
    section = None if mass == 0.0 else AxisSectionOperator(
        full_space=space, marginal_space=marginal_space,
        axis_shape=axis.shape, dims=dims,
        unit=flat_unit, total_weight=mass,
    )
    return _AxisMarginal(retraction, section)
