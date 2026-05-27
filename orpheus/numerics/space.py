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
from numpy.typing import NDArray

__all__ = [
    "FunctionSpace",
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

    Notes
    -----
    The class is **frozen**. Two function spaces are considered
    identical iff their ``(name, shape)`` tuple matches; weights are
    metadata that affect the inner product but not the identity of
    the space. This convention matches the abstract-vector-space
    framing where two copies of :math:`\mathbb{R}^n` are "the same"
    space regardless of which inner product is installed.
    """

    name: str
    shape: tuple[int, ...]
    inner_product_weights: Optional[NDArray] = field(default=None, repr=False)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, FunctionSpace):
            return NotImplemented
        return self.name == other.name and self.shape == other.shape

    def __hash__(self) -> int:
        return hash((self.name, self.shape))

    def __repr__(self) -> str:
        return f"FunctionSpace({self.name!r}, shape={self.shape})"

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
