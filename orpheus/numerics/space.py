r"""Function spaces for the operator-algebra framework.

A *function space* is the domain or codomain of a linear operator. In
the matrix-free transport algebra, every operator :math:`A : V \to W`
acts on discrete flux distributions that live in a definite
:class:`FunctionSpace` — angular-flux space, scalar-flux space, a
boundary-trace space, and so on. Tagging operators with their
:attr:`domain` and :attr:`range` lets the composition machinery in
:mod:`orpheus.numerics.operator` reject ill-formed compositions at
construction time (raising :class:`IncompatibleOperatorComposition`),
preventing the harmful-stub anti-pattern where a downstream Krylov
consumer hits a shape mismatch mid-iteration.

The *Hilbert-adjoint* refinement (:meth:`LinearOperator.adjoint`,
``A.H``) further requires the domain and range to carry an inner
product. :class:`FunctionSpace` stores the inner-product weights as
metadata: the L² inner product
:math:`\langle x, y \rangle_w = \sum_i w_i \, x_i \, y_i` reduces to
the Euclidean :math:`\sum_i x_i \, y_i` when the weights are absent,
which is the natural default for spaces that have no canonical
quadrature (cell-flux, region-flux). For angular-flux and
boundary-trace spaces the quadrature weights ARE the canonical
inner-product weights, so the adjoint identity
:math:`\langle A x, y \rangle_W = \langle x, A^* y \rangle_V`
becomes a non-trivial consistency check (see test
``tests/numerics/test_operator.py::test_hilbert_adjoint_weighted_identity``).

Future direction (Grand Report v3 §5.3 + §6.1)
==============================================

The Grand Report v3 anticipates a richer Space ontology with these
specialisations layered on top of :class:`FunctionSpace`:

* **MeshFunctionSpace** — functions on a structured mesh; carries a
  reference to the :class:`~orpheus.geometry.mesh.Mesh1D` /
  :class:`~orpheus.geometry.mesh.Mesh2D` instance.
* **TraceSpace** — functions on a boundary face; the domain/range
  for :class:`~orpheus.geometry.boundary.BoundaryOperator`. Carries
  a directional tag (``"in"`` vs ``"out"``).
* **RegionSpace** — region-piecewise constant fields (one value per
  homogenised region); used by region-collapsed CP / homogenisation.
* **EnergyGroupSpace** — multi-group flux space; tensored with a
  spatial space to form the full state.
* **DiscreteAngularSpace** — quadrature-tagged angular space carrying
  invariance-group metadata and the underlying
  :class:`~orpheus.numerics.measure.DiscreteMeasure`.

And these compositional dunders:

* ``S * T`` — tensor product; produces a :class:`FunctionSpace` whose
  shape is the concatenation of factor shapes.
* ``S + T`` — direct sum; concatenated dimension on a shared abstract
  type tag.
* ``S.dual()`` — dual space; for inner-product-bearing spaces this is
  isomorphic to ``S`` itself but carries a covariance tag for
  bra-ket-style composition checks.

These are NOT shipped in 9.6 — the file is structured so they can be
slotted in additively without disturbing the :class:`FunctionSpace`
base.

References
----------

* Trefethen, L.N. & Bau, D. (1997). *Numerical Linear Algebra*. SIAM.
  §1 — vector spaces, inner products, the Hilbert adjoint vs.
  representation transpose distinction.
* Reed, M. & Simon, B. (1980). *Methods of Modern Mathematical
  Physics I: Functional Analysis*, §III.6 (Hilbert spaces, Riesz
  representation, the adjoint operator).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal, Optional

import numpy as np
import pint
from numpy.typing import NDArray

__all__ = [
    "DualSpace",
    "FunctionSpace",
    "TensorProductSpace",
    "angular_flux_space",
    "scalar_flux_space",
    "boundary_trace_space",
]


@dataclass(frozen=True)
class FunctionSpace:
    r"""A finite-dimensional vector space of discrete fields.

    Parameters
    ----------
    name : str
        Human-readable identifier. Used by :meth:`__repr__` and by
        :class:`IncompatibleOperatorComposition` error messages. Two
        spaces with the same ``name`` and ``shape`` compare equal even
        if they were constructed via different factory functions —
        ``name`` is the **identity** of the space, not a description
        of its contents.
    shape : tuple[int, ...]
        Tensor shape of the elements. A ``shape=(n_cells,
        n_ordinates, n_groups)`` space contains 3-D arrays; a
        ``shape=(n_cells,)`` space contains 1-D arrays.
    inner_product_weights : NDArray | None, default None
        Diagonal weights for the L² inner product
        :math:`\langle x, y \rangle = \sum_i w_i \, x_i \, y_i`. The
        weights array MUST be broadcast-compatible with ``shape``;
        most commonly a 1-D array along one axis (e.g. quadrature
        weights along the ordinate axis of an angular-flux space).
        ``None`` selects the Euclidean inner product
        :math:`\sum_i x_i \, y_i`.
    units : pint.Unit | None, default None
        Dimensional label of every :class:`~orpheus.numerics.field.Field`
        that lives on this space. Two spaces with the same
        ``(name, shape)`` but different ``units.dimensionality`` are
        structurally different mathematical objects — :meth:`__eq__`
        returns ``False`` in that case so :class:`Field` arithmetic
        across them raises. ``None`` is a permissive sentinel that
        opts out of dimensional checking (treated as compatible with
        any other units in :meth:`__eq__`).

    Notes
    -----
    The class is **frozen**. Two function spaces are considered
    identical iff their ``(name, shape)`` tuple matches AND their
    ``units.dimensionality`` matches (when both sides set units).
    Inner-product weights are metadata that affect the inner product
    but not the identity of the space — this matches the abstract-
    vector-space framing where two copies of :math:`\mathbb{R}^n` are
    "the same" space regardless of which inner product is installed.

    Hashing is on ``(name, shape)`` only; units do not participate.
    Two spaces with the same ``(name, shape)`` but incompatible
    dimensionality will collide in dict keys but resolve via the
    ``__eq__`` dimensionality check — Python's hashable contract
    (equal-implies-same-hash) is preserved.
    """

    name: str
    shape: tuple[int, ...]
    inner_product_weights: Optional[NDArray] = field(default=None, repr=False)
    units: Optional[pint.Unit] = field(default=None)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, FunctionSpace):
            return NotImplemented
        if self.name != other.name or self.shape != other.shape:
            return False
        if self.units is None or other.units is None:
            return True
        return self.units.dimensionality == other.units.dimensionality

    def __hash__(self) -> int:
        return hash((self.name, self.shape))

    def __repr__(self) -> str:
        if self.units is None:
            return f"FunctionSpace({self.name!r}, shape={self.shape})"
        return f"FunctionSpace({self.name!r}, shape={self.shape}, units={self.units!s})"

    # ------------------------------------------------------------------
    # Inner product / norm
    # ------------------------------------------------------------------

    def inner_product(self, x: NDArray, y: NDArray) -> float:
        r"""Return :math:`\langle x, y \rangle`.

        With diagonal weights ``w`` the inner product is the weighted
        sum :math:`\sum_i w_i \, x_i \, y_i`. Without weights it
        reduces to the Euclidean :math:`\sum_i x_i \, y_i`. The
        weights array is broadcast against ``x * y`` so a 1-D weight
        vector along (say) the ordinate axis acts on the full
        ``(n_cells, n_ordinates, n_groups)`` tensor without manual
        reshaping.

        Parameters
        ----------
        x, y : NDArray
            Arrays of shape :attr:`shape`.

        Returns
        -------
        float
            Scalar inner product.
        """
        if self.inner_product_weights is None:
            return float(np.sum(x * y))
        return float(np.sum(self.inner_product_weights * x * y))

    def norm(self, x: NDArray) -> float:
        r"""Return the induced :math:`L^2` norm
        :math:`\sqrt{\langle x, x \rangle}`."""
        return float(np.sqrt(self.inner_product(x, x)))

    # ------------------------------------------------------------------
    # Space algebra (Depth B step D-B)
    # ------------------------------------------------------------------

    def __mul__(self, other: "FunctionSpace") -> "TensorProductSpace":
        r"""Return the tensor product :math:`V \otimes W` of this space
        with ``other``.

        Implements ``V * W`` per grand-report §6.1: the resulting
        :class:`TensorProductSpace` carries the concatenated shape
        ``self.shape + other.shape``, the outer-product inner-product
        weights, and the multiplied units. Associative on its inputs:
        ``(A * B) * C`` and ``A * (B * C)`` both produce a flat
        3-factor :class:`TensorProductSpace`.

        Loadbearing for the Wave T tensor-network rewires per the
        grand report §15.1 (streaming as
        :math:`L = \sum_{\text{axis}} D_{\text{axis}} \otimes \Omega_{\text{axis}} \otimes I_g`)
        and §16A.10 (boundary as
        :math:`B = G_{\text{patch}} \otimes K_\omega \otimes K_g`).
        See ``.claude/plans/depth_b_field_on_function_space.md`` §6
        step D-B for the design.
        """
        if not isinstance(other, FunctionSpace):
            return NotImplemented
        self_factors = (
            self.factors if isinstance(self, TensorProductSpace) else (self,)
        )
        other_factors = (
            other.factors if isinstance(other, TensorProductSpace) else (other,)
        )
        return TensorProductSpace.from_factors(self_factors + other_factors)

    def dual(self) -> "DualSpace":
        r"""Return the dual space :math:`V^*`.

        Under L²-Riesz identification (the standard ORPHEUS setting
        where every :class:`FunctionSpace` carries an inner product),
        :math:`V^*` is isomorphic to :math:`V` itself with a covariance
        tag for bra-ket-style composition tracking. The dual carries
        the same shape, weights, and units as the primal; its
        ``primal`` attribute holds a reference back.

        Used by the Hilbert-adjoint machinery
        (:meth:`~orpheus.numerics.operator.LinearOperator.adjoint`,
        ``A.H``) to track which spaces are codomain-sourced vs
        domain-sourced through operator composition.
        """
        return DualSpace.of(self)


# ---------------------------------------------------------------------------
# Tensor-product and dual space (Depth B step D-B)
# ---------------------------------------------------------------------------


def _tensor_product_inner_weights(
    factors: tuple["FunctionSpace", ...],
) -> Optional[NDArray]:
    r"""Compute the outer-product inner-product weights of a tensor product.

    For factors with weights :math:`w_1, w_2, \ldots, w_k`, the tensor-
    product weights tensor has shape ``factors[0].shape + factors[1].shape
    + ...`` with entries
    :math:`W[i_1, i_2, \ldots, i_k] = w_1[i_1] \cdot w_2[i_2] \cdots w_k[i_k]`.
    Factor weights ``None`` (Euclidean) contribute identity (ones broadcast
    to the factor shape). If ALL factors are Euclidean, the result is
    ``None`` (preserving the Euclidean default — no allocation).
    """
    if all(f.inner_product_weights is None for f in factors):
        return None
    result: Optional[NDArray] = None
    for f in factors:
        w = (
            f.inner_product_weights
            if f.inner_product_weights is not None
            else np.ones(f.shape)
        )
        w = np.broadcast_to(w, f.shape)
        result = w if result is None else np.multiply.outer(result, w)
    return result


def _tensor_product_units(
    factors: tuple["FunctionSpace", ...],
) -> "pint.Unit | None":
    r"""Multiply the units of the factors via pint's unit algebra.

    :math:`V_1 \otimes V_2 \otimes \cdots` lives in units that are the
    product of the factor units. If any factor has ``units=None``, the
    product is ``None`` (cannot infer; documented gap).
    """
    if any(f.units is None for f in factors):
        return None
    result = factors[0].units
    for f in factors[1:]:
        result = result * f.units
    return result


@dataclass(frozen=True)
class TensorProductSpace(FunctionSpace):
    r"""A function space that decomposes as
    :math:`V = V_1 \otimes V_2 \otimes \cdots \otimes V_k`.

    The tensor-product structure makes algebraic identities of operators
    on this space (adjoint distributivity, composition distributivity,
    representation polymorphism) checkable at the type level. See grand-
    report §5.3, §15, §32.4 for the L1 motivation; see
    ``.claude/plans/wave_t_tensor_network.md`` for the production
    consumers being wired in Wave T.

    Construction
    ------------
    Two equivalent paths:

    * **Operator-algebra dispatch** — ``A * B`` where ``A`` and ``B``
      are :class:`FunctionSpace` instances returns a
      :class:`TensorProductSpace`. The dunder is associative on its
      inputs (``(A * B) * C`` flattens to a 3-factor product, never
      nests).
    * **Explicit factory** — :meth:`from_factors` with a tuple of
      factor spaces.

    Notes
    -----
    The class is **frozen**. Identity is the inherited
    ``(name, shape, units.dimensionality)`` tuple, where ``name`` and
    ``shape`` are derived from the factors. Two
    :class:`TensorProductSpace` instances with the same factor sequence
    compare equal even if reached via different composition paths.

    The ``factors`` field is metadata that supports introspection (e.g.,
    operator-algebra factor matching for ``(A & B) ∘ (C & D)`` rewriting)
    and is not part of the identity.

    Parameters
    ----------
    factors : tuple[FunctionSpace, ...]
        The factor spaces, in order. Should have ``len >= 2`` for a
        meaningful tensor product; trivial 1-factor TensorProductSpaces
        are permitted by the dataclass but produced only via the
        :meth:`__mul__` flattening edge cases.
    """

    factors: tuple["FunctionSpace", ...] = field(default=(), compare=False, repr=False)

    @classmethod
    def from_factors(
        cls, factors: tuple["FunctionSpace", ...],
    ) -> "TensorProductSpace":
        r"""Construct a :class:`TensorProductSpace` from a tuple of factor
        spaces.

        Derives:
        * ``name`` from ``" ⊗ ".join(f.name for f in factors)``
        * ``shape`` from concatenated factor shapes
        * ``inner_product_weights`` from the outer product of factor
          weights (``None`` if all factors are Euclidean)
        * ``units`` from the product of factor units (``None`` if any
          factor is unitless)
        """
        if len(factors) < 2:
            raise ValueError(
                f"TensorProductSpace.from_factors requires at least 2 "
                f"factors; got {len(factors)}"
            )
        name = " ⊗ ".join(f.name for f in factors)
        shape: tuple[int, ...] = ()
        for f in factors:
            shape = shape + f.shape
        return cls(
            name=name,
            shape=shape,
            inner_product_weights=_tensor_product_inner_weights(factors),
            units=_tensor_product_units(factors),
            factors=factors,
        )

    def __repr__(self) -> str:
        if self.units is None:
            return f"TensorProductSpace({self.name!r}, shape={self.shape})"
        return (
            f"TensorProductSpace({self.name!r}, shape={self.shape}, "
            f"units={self.units!s})"
        )


@dataclass(frozen=True)
class DualSpace(FunctionSpace):
    r"""The dual :math:`V^*` of a :class:`FunctionSpace`.

    Under L²-Riesz identification (the standard ORPHEUS setting where
    every :class:`FunctionSpace` carries an inner product),
    :math:`V^*` is isomorphic to :math:`V` itself but carries a
    covariance tag that the operator-algebra layer reads through to
    track which spaces participate as bras vs. kets in composition
    chains. The dual carries the same ``shape``, ``inner_product_weights``,
    and ``units`` as the primal; the ``primal`` field is the
    introspection link.

    Used by the Hilbert-adjoint machinery
    (:meth:`~orpheus.numerics.operator.LinearOperator.adjoint`,
    ``A.H``) — taking ``A.H`` swaps domain ↔ codomain AND flips both
    to their duals, so the adjoint's domain is the original's codomain
    dual.

    Notes
    -----
    ``V.dual().dual() == V`` is enforced by :meth:`of` recognising a
    :class:`DualSpace` argument and returning its primal (idempotency).

    Parameters
    ----------
    primal : FunctionSpace
        The primal space :math:`V` of which this is the dual.
    """

    primal: Optional["FunctionSpace"] = field(default=None, compare=False, repr=False)

    @classmethod
    def of(cls, primal: "FunctionSpace") -> "FunctionSpace":
        r"""Construct the dual of ``primal``.

        Idempotent: ``of(of(V)) == V`` (returns the primal of a passed
        :class:`DualSpace`, never wraps twice).
        """
        if isinstance(primal, DualSpace):
            return primal.primal
        return cls(
            name=f"{primal.name}*",
            shape=primal.shape,
            inner_product_weights=primal.inner_product_weights,
            units=primal.units,
            primal=primal,
        )

    def __repr__(self) -> str:
        return f"DualSpace({self.name!r}, shape={self.shape})"


# ---------------------------------------------------------------------------
# Pre-populated common-space factories
# ---------------------------------------------------------------------------
#
# These are *factory functions*, not module-level singletons, because
# the shape of an angular-flux space depends on the mesh dimensions —
# every solver instance carries its own (n_cells, n_ordinates, n_groups)
# triple. Factories give the caller named, well-typed construction
# sites without forcing premature commitment to a global instance.


def angular_flux_space(
    n_cells: int,
    n_ordinates: int,
    n_groups: int,
    *,
    quadrature_weights: NDArray | None = None,
) -> FunctionSpace:
    r"""Construct the angular-flux space for an SN solve.

    Shape is ``(n_cells, n_ordinates, n_groups)``. When
    ``quadrature_weights`` is provided, it is broadcast along the
    ordinate axis as the inner-product metadata so the canonical
    angular inner product
    :math:`\langle \psi, \varphi \rangle_\Omega = \sum_n w_n \,
    \psi_n \, \varphi_n` (summed over cells / groups too) becomes
    :meth:`FunctionSpace.inner_product`.

    Parameters
    ----------
    n_cells, n_ordinates, n_groups : int
        Tensor dimensions.
    quadrature_weights : NDArray, optional
        Shape ``(n_ordinates,)`` quadrature weights along the
        ordinate axis. Reshaped to ``(1, n_ordinates, 1)`` so it
        broadcasts against the full angular-flux tensor.
    """
    weights: NDArray | None = None
    if quadrature_weights is not None:
        w = np.asarray(quadrature_weights)
        if w.shape != (n_ordinates,):
            raise ValueError(
                f"quadrature_weights must have shape ({n_ordinates},), "
                f"got {w.shape}"
            )
        weights = w.reshape(1, n_ordinates, 1)
    return FunctionSpace(
        name="angular_flux",
        shape=(n_cells, n_ordinates, n_groups),
        inner_product_weights=weights,
    )


def scalar_flux_space(n_cells: int, n_groups: int) -> FunctionSpace:
    r"""Construct the scalar-flux space.

    Shape is ``(n_cells, n_groups)``. No inner-product weights — the
    canonical inner product on scalar-flux space is the Euclidean
    sum (or, equivalently, the volume-weighted L² inner product when
    the cell volumes are absorbed into the operator that produced
    :math:`\phi`; see :meth:`Mesh1D.volume_measure` for the
    volume-weighted variant).
    """
    return FunctionSpace(
        name="scalar_flux",
        shape=(n_cells, n_groups),
    )


def boundary_trace_space(
    direction: Literal["in", "out"],
    n_ordinates: int,
    n_groups: int,
    *,
    quadrature_weights: NDArray | None = None,
) -> FunctionSpace:
    r"""Construct the directional boundary-trace space.

    Shape is ``(n_ordinates, n_groups)``. The ``direction`` tag
    distinguishes the *incoming* and *outgoing* faces of a boundary
    operator: a :class:`BoundaryOperator` ``B`` has
    ``B.domain = boundary_trace_space("out", ...)`` and
    ``B.codomain = boundary_trace_space("in", ...)`` (it consumes
    outgoing flux, produces incoming flux).

    Parameters
    ----------
    direction : ``"in"`` | ``"out"``
        Whether this trace describes flux flowing **into** the
        domain (incoming) or **out of** the domain (outgoing).
    n_ordinates, n_groups : int
        Tensor dimensions.
    quadrature_weights : NDArray, optional
        Shape ``(n_ordinates,)`` weights along the ordinate axis.
    """
    weights: NDArray | None = None
    if quadrature_weights is not None:
        w = np.asarray(quadrature_weights)
        if w.shape != (n_ordinates,):
            raise ValueError(
                f"quadrature_weights must have shape ({n_ordinates},), "
                f"got {w.shape}"
            )
        weights = w.reshape(n_ordinates, 1)
    return FunctionSpace(
        name=f"boundary_trace_{direction}",
        shape=(n_ordinates, n_groups),
        inner_product_weights=weights,
    )
