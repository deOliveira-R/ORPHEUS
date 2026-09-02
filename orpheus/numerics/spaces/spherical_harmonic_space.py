r"""Function space for spherical-harmonic moment coefficients.

This module ships the typed home of the ERR-039 Gram matrix
:math:`g_C = \mathrm{diag}(4\pi/(2\ell+1))`. Pre-frame the SH Gram lived as
a free literal (the ``two_l_plus_one`` array, which is
:math:`4\pi \cdot g_C^{-1}` wearing a disguise) on the since-retired
harmonic-reconstruction operator, and as a prose warning on the moment
projection's representation transpose. It now lives here: the space carries
the metric as :attr:`inner_product_weights`, and the metric formula itself is
sourced from :class:`~orpheus.numerics.basis.SphericalHarmonicBasis` so the
:math:`(2\ell+1)` literal exists in exactly one place.

What this enables
=================

The discrete spherical-harmonic :class:`~orpheus.numerics.frame.GalerkinFrame`'s
codomain is this space **re-dressed by the frame** (F-0,
``frame_square_recarve.md``): :meth:`SphericalHarmonicSpace.from_L` carries the
CONTINUUM Gram :math:`g_C = 4\pi/(2\ell+1)`, and :attr:`FrameBase.basis_space
<orpheus.numerics.frame.FrameBase.basis_space>` REPLACES the metric with the
PARSEVAL metric — the inverse of the frame's discrete Gram,
:math:`(2\ell+1)/4\pi` on a degree-exact sphere rule — because the carried
moments are COVARIANT (:math:`\varphi = Gc`) and only :math:`G^{-1}` makes
analysis an isometry onto its image (`[M]` ``scratch/probe_f1_parseval.py``,
2026-08-24: continuum-side Parseval ratio 118.7 vs 1.000). The generic
``AdjointOperator`` machinery then computes ``frame.analysis.H`` as the
physical :math:`S_0 \circ G^{-1} = R/W` with no bespoke code.

ERR-039 in one sentence: the addition-theorem reconstruction
:math:`R = (2\ell+1) S_0` and the analysis face's Hilbert adjoint
:math:`\Pi^* = R/W` differ by exactly the total weight :math:`W = 4\pi`. With
the metric carried on the space and the basis providing the convention, they
fall out as two derived expressions from a common ground.

References
----------

* Grand Report v3 §5.3 — Space hierarchy.
* Grand Report v3 §6.3 — ``.T`` (representation transpose) vs ``.H``
  (Hilbert adjoint).
* Grand Report v3 §19 — Harmonic projection.
* :mod:`orpheus.numerics.basis.spherical_harmonic_basis` — the basis,
  the SH convention, and the source of
  :attr:`SphericalHarmonicSpace.metric_per_ell`.
* ERR-039 entry: ``docs/theory/verification/error_catalog.rst``.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Optional

import numpy as np
from numpy.typing import NDArray

from orpheus.numerics.basis.spherical_harmonic_basis import SphericalHarmonicBasis
from orpheus.numerics.space import FunctionSpace


__all__ = ["SphericalHarmonicSpace"]


# ─────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────


def _padded_metric_tensor(L: int, metric_per_ell: NDArray) -> NDArray:
    r"""Broadcast the per-:math:`\ell` metric to the padded ``(L+1, 2L+1)`` storage layout.

    The :class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`
    storage convention is to allocate ``(L+1, 2L+1)`` slots with the
    addition-theorem-shifted :math:`m`-index ``[l + m]`` — entries
    outside :math:`|m| \le \ell` are zero by construction. The
    inner-product weights MUST share that layout so the broadcast
    :math:`\sum_{\ell m} w_{\ell m} \, x_{\ell m} \, y_{\ell m}` reduces
    correctly to :math:`\sum_\ell (4\pi/(2\ell+1)) \sum_{m=-\ell}^{\ell}
    x_\ell^m y_\ell^m` (the padded slots contribute zero either way).

    Parameters
    ----------
    L : int
        Maximum harmonic degree.
    metric_per_ell : NDArray, shape ``(L+1,)``
        The :math:`4\pi/(2\ell+1)` diagonal per :math:`\ell`. Sourced
        from :attr:`SphericalHarmonicBasis.metric_per_ell` so the formula
        lives in exactly one place.

    Returns
    -------
    NDArray, shape ``(L+1, 2L+1)``
        Padded metric: row :math:`\ell` holds ``metric_per_ell[ell]`` in
        the :math:`2\ell+1` valid slots ``[0, 2*ell]``, zero elsewhere.
    """
    cols = np.arange(2 * L + 1)
    rows = np.arange(L + 1)
    valid_mask = cols[None, :] <= 2 * rows[:, None]   # (L+1, 2L+1) bool
    return np.where(valid_mask, metric_per_ell[:, None], 0.0)


# ─────────────────────────────────────────────────────────────────────
# Space class
# ─────────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class SphericalHarmonicSpace(FunctionSpace):
    r"""Function space of real-spherical-harmonic moment coefficients up to degree :math:`L`.

    Parameters
    ----------
    name : str
        Inherited from :class:`FunctionSpace`. Convention:
        ``"spherical_harmonic_space"``.
    shape : tuple[int, ...]
        Inherited from :class:`FunctionSpace`. MUST equal
        ``(L + 1, 2 * L + 1)``; ``__post_init__`` checks.
    inner_product_weights : NDArray, optional
        Inherited from :class:`FunctionSpace`. The padded ``(L+1, 2L+1)``
        metric tensor, zero in the :math:`|m| > \ell` padding. WHICH metric
        depends on the instance's provenance: :meth:`from_L` installs the
        CONTINUUM Gram (row :math:`\ell` holds :math:`4\pi/(2\ell+1)`),
        while a frame's dressed ``basis_space`` — built from this same class
        via :func:`dataclasses.replace` — carries the PARSEVAL inverse
        (:math:`(2\ell+1)/4\pi` on a degree-exact rule; F-0,
        :attr:`FrameBase.basis_space
        <orpheus.numerics.frame.FrameBase.basis_space>`).
    L : int, default 0
        Maximum harmonic degree retained. Must satisfy
        ``shape == (L + 1, 2 * L + 1)``.

    Notes
    -----
    Frozen dataclass with the same subclassing pattern as
    :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`. The default
    ``L = 0`` is required by dataclass-inheritance rules: it follows
    :attr:`FunctionSpace.inner_product_weights` which has a default,
    so every subsequent field must also have one.

    Equality and hashing are by ``(name, shape)`` alone (inherited from
    :class:`FunctionSpace`); ``shape == (L+1, 2L+1)`` already encodes
    :math:`L`, so equal-:math:`L` spaces compare equal even when their
    ``inner_product_weights`` arrays are distinct ``ndarray`` objects.
    """

    L: int = 0

    def __post_init__(self) -> None:
        # The base guards first (one-metric-source exclusivity + metric
        # admission — P7): an override that skips them silently opts the
        # whole subclass out of the construction contract.
        super().__post_init__()
        expected = (self.L + 1, 2 * self.L + 1)
        if self.shape != expected:
            raise ValueError(
                f"SphericalHarmonicSpace: shape={self.shape} inconsistent with "
                f"L={self.L}; expected shape={expected}."
            )

    # ── Equality / hashing inherited from FunctionSpace ───────────────
    #
    # FunctionSpace defines explicit __eq__ and __hash__ on (name, shape).
    # The @dataclass(frozen=True) decorator on this subclass would
    # otherwise auto-generate its own __eq__ / __hash__ that compare ALL
    # fields (including the ndarray inner_product_weights) — and ndarray
    # equality returns an array, raising at use time.  Explicitly
    # delegating restores the (name, shape) identity convention; ``L``
    # is already encoded in ``shape`` so excluding it from equality is
    # consistent with the abstract-vector-space framing where the
    # identity of a space is its type tag + dimension.

    def __eq__(self, other: object) -> bool:
        return FunctionSpace.__eq__(self, other)

    def __hash__(self) -> int:
        return FunctionSpace.__hash__(self)

    # ── Constructor ──────────────────────────────────────────────────

    @classmethod
    def from_L(cls, L: int) -> "SphericalHarmonicSpace":
        r"""Construct the canonical SH space for truncation degree :math:`L`.

        Builds the metric tensor from :class:`SphericalHarmonicBasis` so
        the :math:`(2\ell+1)` / :math:`4\pi/(2\ell+1)` formulas live in
        exactly one place (the basis).

        Parameters
        ----------
        L : int
            Maximum harmonic degree. Must be non-negative.

        Returns
        -------
        SphericalHarmonicSpace
            With ``name="spherical_harmonic_space"``,
            ``shape=(L+1, 2L+1)``, and ``inner_product_weights`` carrying
            the padded :math:`4\pi/(2\ell+1)` CONTINUUM Gram. (A frame's
            ``basis_space`` is this object re-dressed with the discrete
            Parseval inverse — :attr:`FrameBase.basis_space
            <orpheus.numerics.frame.FrameBase.basis_space>`, F-0.)
        """
        basis = SphericalHarmonicBasis(L=L)
        weights = _padded_metric_tensor(L, basis.metric_per_ell)
        return cls(
            name="spherical_harmonic_space",
            shape=(L + 1, 2 * L + 1),
            inner_product_weights=weights,
            L=L,
        )

    # ── Delegated properties (single source of truth in the basis) ───

    @property
    def basis(self) -> SphericalHarmonicBasis:
        r"""The associated :class:`SphericalHarmonicBasis`.

        Uniquely determined by :attr:`L`; constructed on access.
        """
        return SphericalHarmonicBasis(L=self.L)

    @property
    def metric_per_ell(self) -> NDArray:
        r"""The :math:`4\pi/(2\ell+1)` diagonal per :math:`\ell`, shape ``(L+1,)``.

        Delegated to :attr:`SphericalHarmonicBasis.metric_per_ell` so the
        SH convention's formula has a single canonical home.
        """
        return self.basis.metric_per_ell

    def truncated(self, L_new: int) -> "SphericalHarmonicSpace":
        r"""This family's space at the lower order ``L_new`` — the head a moment field truncates TO.

        The moment carrier truncates by asking its angular HEAD factor for
        the same family one order down, never by re-minting a family from
        an integer (#429 tracker 2.5): a Legendre head truncates to a
        Legendre head, a spherical-harmonic head to this. Returns the
        CONTINUUM-metric member (:meth:`from_L`) under THIS head's own
        name — the head keeps its identity and only its order moves — so
        a frame-dressed head drops its Parseval dressing on truncation,
        exactly as it always did: the metric is the frame's to install,
        never the field's to carry across orders. (The production name
        is the class's canonical one, so this is bit-identical to the
        ``from_L(L_new)`` it replaced; what it refuses to do is hand a
        renamed head back the default name — the tell of an integer mint.)
        """
        if not 0 <= L_new <= self.L:
            raise ValueError(
                f"SphericalHarmonicSpace.truncated: L_new={L_new} must lie "
                f"in [0, {self.L}]."
            )
        return replace(type(self).from_L(L_new), name=self.name)

    @property
    def addition_theorem_factor(self) -> NDArray:
        r"""The :math:`(2\ell+1)` factor per :math:`\ell`, shape ``(L+1,)``.

        Delegated to :attr:`SphericalHarmonicBasis.addition_theorem_factor` —
        the canonical source of the addition-theorem factor that the discrete
        spherical-harmonic :class:`~orpheus.numerics.frame.GalerkinFrame`'s
        reconstruction face reads.
        """
        return self.basis.addition_theorem_factor
