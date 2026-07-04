r"""Function space for within-cell spatial-moment (tensor-Legendre DG) coefficients.

This module mints the typed home of the **spatial** moment basis — the
within-cell tensor-product Legendre Discontinuous-Galerkin (DG) basis the
Linear-Discontinuous (LD) closure carries: how the angular flux
:math:`\psi` varies in space *inside one mesh cell*. It is the spatial
sibling of :class:`~orpheus.numerics.spaces.spherical_harmonic_space.SphericalHarmonicSpace`
(the **angular** moment basis), and the naming collision between the two
"moment" notions is precisely the trap this type exists to dispel:

================  ==================================  ========================
Moment kind       basis                               varies over
================  ==================================  ========================
**angular**       real spherical harmonics            direction :math:`\Omega`
                  :math:`\{Y_\ell^m\}`
                  (:class:`SphericalHarmonicSpace`)
**spatial**       tensor-Legendre :math:`\{1, P_1\}`  space :math:`x` (within a
                  per axis                            cell)
                  (:class:`SpatialMomentSpace`)
================  ==================================  ========================

The two axes are **orthogonal**: a flux representation can carry a
spatial-moment axis on TOP of any angular representation (full angular
:math:`\psi`, or angular moments :math:`\phi_\ell^m`). The spatial-moment
axis is the within-cell DG slope structure that the LD scheme produces and
that the diffusion-limit-consistent scattering source
:math:`S_{\rm full} = \Sigma_s \otimes I_{\rm spatial}` requires to travel
between source-iteration sweeps (#240 D5b-S3).

The Kronecker / tensor-product basis
=====================================

The per-cell DG basis is the tensor (Kronecker) product of a 1-D basis
applied per spatial axis. The 1-D basis size is :attr:`per_axis`:

* ``per_axis == 1`` — the cell-average basis :math:`\{1\}`. This is what
  Diamond-Difference / Step carry (a slopeless, cell-average closure).
* ``per_axis == 2`` — the linear basis :math:`\{1, P_1\}` (constant +
  first Legendre). This is what Linear-Discontinuous carries.

For ``ndim`` spatial axes the total within-cell moment count is
:attr:`per_axis` :sup:`ndim` (DD: ``1``; LD-1D: ``2``; LD-2D: ``4``;
LD-3D: ``8``). The ordering is the Kronecker product of the per-axis
bases, **x-outer / y-inner** (matching
:func:`orpheus.transport.spatial._ubld.assemble_ubld`), so the slot-0 entry is
always the all-:math:`P_0` cell average:

* d = 1 (``per_axis = 2``): :math:`[\bar\psi,\ \hat\psi_x]`
* d = 2 (``per_axis = 2``): :math:`[\bar\psi,\ \hat\psi_y,\ \hat\psi_x,\ \hat\psi_{xy}]`
* d = 3 (``per_axis = 2``): the 8-vector
  :math:`[\bar\psi,\ \hat\psi_z,\ \hat\psi_y,\ \hat\psi_{yz},\ \hat\psi_x,\ \hat\psi_{xz},\ \hat\psi_{xy},\ \hat\psi_{xyz}]`

The slot-0 (cell-average) convention is single-sourced from
:data:`orpheus.numerics.moment_layout.AVERAGE_MOMENT`, the same constant every
moment consumer reduces on — this space does not re-spell the literal.

Construct-general / select-narrow (#240 D5b-S3-A0)
=================================================

This type is the CAPABILITY half of the spatial-moment iterate. It is
composed (optionally) into the bulk-field spaces via the field-space
factories, gated so that ``per_axis == 1`` appends NO factor — the
DD/Step field space stays byte-identical to its pre-S3 shape. No
production field selects the axis yet; the iterate / cell-emit / source
seams that DO select it are the next sub-step (S3-A). The "append iff
> 1" policy is single-sourced from
:func:`orpheus.numerics.moment_layout.face_moment_tail` (the cell analogue,
:func:`~orpheus.numerics.spaces.spatial_moment_space.spatial_moment_tail`,
delegates to it so the policy lives in exactly one place).

References
----------

* :mod:`orpheus.numerics.moment_layout` — the physics-free moment-layout
  policy (``AVERAGE_MOMENT`` slot-0 + the "append iff > 1" tail) this space
  surfaces; single-sourced there (#245).
* :mod:`orpheus.transport.spatial._ubld` — the UBLD cell assembler; the Kronecker
  moment ordering this space's slot layout mirrors.
* :class:`orpheus.transport.spatial.scheme.DiscretizationSchemeBase` — carries
  ``spatial_basis_per_axis`` (DD/Step = 1, LD = 2), the per-axis basis
  size this space's :attr:`per_axis` is derived from.
* :class:`~orpheus.numerics.spaces.spherical_harmonic_space.SphericalHarmonicSpace`
  — the angular-moment sibling and the construction mold (``from_*``
  factory, size in ``shape``, ``find_factor``-queryable).
* ``.claude/plans/issue_240_d5b_s3_crosswalk.md`` — the S3 design record
  (angular-vs-spatial moment distinction, the FP resolution).
"""

from __future__ import annotations

from dataclasses import dataclass

from orpheus.numerics.space import FunctionSpace

# The slot-0 cell-average index + the "append iff > 1" tail policy are the
# physics-free moment-layout primitives, now homed in the leaf numerics
# module ``orpheus.numerics.moment_layout`` (#245).  A top-level import is
# safe — ``moment_layout`` is leaf (stdlib only), so it cannot re-introduce
# the numerics → SN cycle the old deferred ``_ubld`` import worked around.
from orpheus.numerics.moment_layout import (
    AVERAGE_MOMENT,
    cell_moment_count,
    face_moment_tail,
)


__all__ = ["SpatialMomentSpace", "spatial_moment_tail"]


# ─────────────────────────────────────────────────────────────────────
# The "append iff > 1" policy — the CELL analogue of face_moment_tail
# ─────────────────────────────────────────────────────────────────────


def spatial_moment_tail(n_cell_moments: int) -> tuple[int, ...]:
    r"""Trailing spatial-moment-axis shape suffix for a bulk-field buffer.

    The CELL analogue of :func:`orpheus.numerics.moment_layout.face_moment_tail`
    (which sizes the per-FACE transverse cochain). A multi-moment closure
    (LD, ``n_cell_moments == per_axis**ndim > 1``) carries a trailing
    spatial-moment axis on the bulk field; a cell-average closure (DD/Step,
    ``== 1``) leaves the field rank untouched (``()`` — NO length-1 axis
    appended) so its buffers / spaces stay byte-identical (#240 D5b — the
    backward-compat invariant).

    Delegates to :func:`~orpheus.numerics.moment_layout.face_moment_tail` so
    the "append iff > 1" decision lives in EXACTLY ONE place — the cell-moment
    tail and the face-cochain tail must never disagree on the policy
    (``coding-elegance`` Pattern 7: normalise the convention at one site).
    """
    return face_moment_tail(n_cell_moments)


# ─────────────────────────────────────────────────────────────────────
# Space class
# ─────────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class SpatialMomentSpace(FunctionSpace):
    r"""Function space of within-cell tensor-Legendre DG moment coefficients.

    A peer of
    :class:`~orpheus.numerics.spaces.spherical_harmonic_space.SphericalHarmonicSpace`
    in the Space hierarchy: a tensor-product moment carrier whose size is
    encoded in :attr:`shape` and whose metadata
    (:attr:`per_axis`, :attr:`ndim`) describe the tensor-Legendre basis.
    Composes via ``*`` into a :class:`~orpheus.numerics.space.TensorProductSpace`
    and is ``find_factor``-queryable (see
    :meth:`~orpheus.numerics.space.TensorProductSpace.find_factor`).

    Parameters
    ----------
    name : str
        Inherited from :class:`FunctionSpace`. Convention:
        ``"spatial_moment_space"``.
    shape : tuple[int, ...]
        Inherited from :class:`FunctionSpace`. MUST equal
        ``(per_axis ** ndim,)`` — a single trailing moment axis of length
        :math:`(\text{per\_axis})^{\text{ndim}}`. ``__post_init__`` checks.
    per_axis : int, default 1
        The 1-D basis size (per spatial axis): ``1`` for the cell-average
        closures (DD / Step, basis :math:`\{1\}`); ``2`` for
        Linear-Discontinuous (basis :math:`\{1, P_1\}`). Mirrors
        :attr:`DiscretizationSchemeBase.spatial_basis_per_axis`.
    ndim : int, default 0
        Number of spatial dimensions over which the tensor product is
        taken. Must satisfy ``shape == (per_axis ** ndim,)``.

    Notes
    -----
    Frozen dataclass with the same subclassing pattern as
    :class:`SphericalHarmonicSpace`. The defaults ``per_axis = 1``,
    ``ndim = 0`` are required by dataclass-inheritance rules (they follow
    :attr:`FunctionSpace.inner_product_weights` which has a default), and
    they describe the degenerate cell-average / no-spatial-axes space of
    size :math:`1^0 = 1`.

    Equality and hashing are by ``(name, shape)`` alone (inherited from
    :class:`FunctionSpace`). ``shape == (per_axis ** ndim,)`` encodes the
    SIZE but NOT the ``(per_axis, ndim)`` factorisation — a
    ``(per_axis=4, ndim=1)`` space and a ``(per_axis=2, ndim=2)`` space
    both have ``shape == (4,)`` and so compare equal. The two never
    coexist on one mesh (``per_axis`` is the scheme's, ``ndim`` the
    mesh's), so the size-identity convention is the right one — it matches
    the abstract-vector-space framing (identity = type tag + dimension)
    that :class:`SphericalHarmonicSpace` and
    :class:`~orpheus.numerics.space.TensorProductSpace` already follow.

    The **inner product is Euclidean** (``inner_product_weights = None``):
    a within-cell DG moment carries no canonical diagonal metric of its
    own here (the cell-volume / mass-matrix weighting lives on the UBLD
    cell operator, not on the field space — issue #207's "units do not
    live on the space" framing).
    """

    per_axis: int = 1
    ndim: int = 0

    def __post_init__(self) -> None:
        if self.per_axis < 1:
            raise ValueError(
                f"SpatialMomentSpace: per_axis={self.per_axis} must be >= 1 "
                f"(1 = cell-average DD/Step; 2 = linear-discontinuous)."
            )
        if self.ndim < 0:
            raise ValueError(
                f"SpatialMomentSpace: ndim={self.ndim} must be >= 0."
            )
        expected = (cell_moment_count(self.per_axis, self.ndim),)
        if self.shape != expected:
            raise ValueError(
                f"SpatialMomentSpace: shape={self.shape} inconsistent with "
                f"per_axis={self.per_axis}, ndim={self.ndim}; expected "
                f"shape={expected}."
            )

    # ── Equality / hashing inherited from FunctionSpace ───────────────
    #
    # The @dataclass(frozen=True) decorator on this subclass would
    # otherwise auto-generate an __eq__/__hash__ comparing ALL fields
    # (including the ndarray inner_product_weights, whose == raises the
    # ambiguous-truth-value error).  Delegating restores the (name,
    # shape) identity convention shared by every FunctionSpace subclass.

    def __eq__(self, other: object) -> bool:
        return FunctionSpace.__eq__(self, other)

    def __hash__(self) -> int:
        return FunctionSpace.__hash__(self)

    # ── Constructor ──────────────────────────────────────────────────

    @classmethod
    def from_per_axis(cls, per_axis: int, ndim: int) -> "SpatialMomentSpace":
        r"""Construct the spatial-moment space for a per-axis basis size and dimension.

        The :class:`SpatialMomentSpace` analogue of
        :meth:`SphericalHarmonicSpace.from_L` — derives the
        :math:`(\text{per\_axis})^{\text{ndim}}` total moment count and
        encodes it in ``shape``.

        Parameters
        ----------
        per_axis : int
            The 1-D basis size per spatial axis (DD/Step = 1, LD = 2).
            Typically sourced from
            :attr:`DiscretizationSchemeBase.spatial_basis_per_axis`.
        ndim : int
            Number of spatial dimensions (the mesh's ``ndim``).

        Returns
        -------
        SpatialMomentSpace
            With ``name="spatial_moment_space"``,
            ``shape=(per_axis ** ndim,)``.
        """
        return cls(
            name="spatial_moment_space",
            shape=(cell_moment_count(per_axis, ndim),),
            per_axis=per_axis,
            ndim=ndim,
        )

    # ── Derived metadata ─────────────────────────────────────────────

    @property
    def n_moments(self) -> int:
        r"""The total within-cell moment count :math:`(\text{per\_axis})^{\text{ndim}}`.

        Equal to ``shape[0]``. DD/Step → ``1``; LD-1D → ``2``; LD-2D →
        ``4``; LD-3D → ``8``.
        """
        return self.shape[0]

    @property
    def average_moment_index(self) -> int:
        r"""The slot of the cell-average (all-:math:`P_0`) moment — always 0.

        Single-sourced from :data:`orpheus.numerics.moment_layout.AVERAGE_MOMENT`
        (the Kronecker layout puts the all-:math:`P_0` moment first); this
        property surfaces the convention on the space without re-spelling
        the literal.
        """
        return AVERAGE_MOMENT
