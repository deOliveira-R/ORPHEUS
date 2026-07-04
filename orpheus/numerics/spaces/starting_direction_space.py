r"""The starting-direction space — per-level ψ½ carrier with the ghost metric.

The function space of the **starting-direction block** of the augmented
transport composite :math:`V = V_{\rm bulk} \oplus V_{\rm trace} \oplus
V_{\rm sd}` (#282 route (a), campaign #280 phase 2.5d). Its elements are
the per-μ-level half-angle fluxes :math:`\psi_{1/2}` — the angular flux
evaluated at a level's starting-direction edge :math:`\mu_{\rm start}`
(Hébert §3.9.4: the :math:`\mu = \pm 1` rays of a sphere are straight
characteristics) — promoted from a lagged solver-internal estimate to
**typed state** so the sweep is a genuine triangular solve (the #282
back edge dissolves).

Which levels carry a block (ruling R12a)
========================================

Presence is keyed PER LEVEL by the structural predicate

    **the level's first-ordinate raw Morel–Montry weight
    τ_raw ∈ (0, 1) exclusive**

— i.e. the M-M half-angle recurrence *genuinely consumes independent
starting-direction state*. The trichotomy, bit-exact on the production
quadratures (τ_raw from
:func:`~orpheus.sn.spatial.pole_angular_closure.morel_montry_tau_raw_per_level`):

* ``τ_raw = 0`` — the starting direction coincides with the level's
  first node (cylinder *product* rules: :math:`\eta_0 = \eta_{1/2} =
  -\sin\theta` bit-exactly, the #229 clamp fact). The seed would be a
  rank-duplicate of :math:`\psi_0` — NO block.
* ``τ_raw = 1`` — the first node sits ON the level's second edge
  (cylinder *level-symmetric* rules: duplicate-η nodes collapse the
  midpoint edge onto :math:`\eta_0`), so the seed's only consumption
  path — the recurrence weight :math:`(1-\tau_0)` — vanishes. The
  seed would be dead state — NO block. (This is why the measured
  cylinder-LS seed sensitivity is 0.0-bit.)
* ``τ_raw ∈ (0,1)`` — the recurrence consumes the seed with a genuine
  weight (sphere Gauss–Legendre: :math:`\tau_{\rm raw,0} \approx
  0.39\text{–}0.42`). The level CARRIES the block.

R12a refines the R12 letter ("μ_start ∉ the level's μ-nodes"), whose
claimed equivalence to ``τ_raw ≠ 0`` is empirically false on
level-symmetric cylinder rules (μ_start ∉ nodes there, yet τ_raw = 1 —
dead). The predicate is evaluated by
:attr:`~orpheus.sn.mesh.augmented_mesh.SNMesh.starting_direction_levels`;
this space is deliberately quadrature-blind (pure layout + metric).

Layout — one flat backing buffer, shaped views (R13)
====================================================

Per ruling R13 both directions are STATE (the outward corner's defect
row must be linear in state) and the (R, μ = ∓1) corner slots land with
the carrier. The flat layout, mirroring the
:class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux`
flat-backing precedent::

    for level in levels (ascending):          # the DAG order
        for sign in (-1, +1):                 # seed⁻ ≺ seed⁺ (pole continuation)
            cells  (ng, nx)   — ψ½ at every radial cell
            corner (ng,)      — ψ½ at r = R (inflow: given data / outflow: defect)

Access goes through :meth:`cells_view` / :meth:`corner_view` — slice
views into the flat buffer, no copies (the exact
:class:`~orpheus.numerics.face_layout.FaceLayout` discipline, with the
``(level, sign)`` key typed instead of stringly).

The ghost metric — ALL weights zero (structural, not fabricated)
================================================================

Every inner-product weight of this space is **exactly zero**: the
angular through-flux factor :math:`(1-\mu^2)` vanishes at
:math:`\mu = \pm 1` — the SAME structural fact as the α-dome endpoints
:math:`\alpha_{1/2} = \alpha_{N+1/2} = 0`. The starting-direction rays
carry zero quadrature measure; a nonzero weight would be a fabricated
volume. Consequences (all inherited from
:class:`~orpheus.numerics.space.FunctionSpace` with zero weights — no
overrides):

* :meth:`~orpheus.numerics.space.FunctionSpace.apply_metric` zeroes the
  block (:math:`G\odot x = 0`);
* :meth:`~orpheus.numerics.space.FunctionSpace.apply_inverse_metric` is
  the masked Moore–Penrose pseudo-inverse — zero on the whole block;
* :meth:`~orpheus.numerics.space.FunctionSpace.inner_product`
  contributes exactly ``0.0``.

**Honest-scope note (vv-principles Mode 12, gate spec §16.A A4):** the
zero weight puts the seed rows inside the G-adjoint reciprocity
functional's invariance group — the G3 reciprocity gate is IDENTICALLY
blind to any error in the seed-row transpose, at every tolerance, in
every regime. The seed rows are metric-invisible YET ACTIVE; they are
constrained by the §16.B direct-solver closed-form pin, the §16.C
solve∘apply residual over the FULL augmented field, and 2.5b's
Euclidean :math:`M^{\mathsf T}` oracle — NEVER by G3.

References
==========

* Hébert, A. (2009). *Applied Reactor Physics*. §3.9.4 (the
  starting-direction equation, Eqs. 3.432–3.435).
* GH #282 (the spherical seed lag), #280 (the walk unification), #229
  (the cylinder τ clamp fact).
* ``.claude/plans/stencil_assembly_dsa_roadmap.md`` — rulings R12/R12a,
  R13; ``a3_solve_transpose_verification.md`` §16.A (the carrier gates).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

from orpheus.numerics.space import FunctionSpace

__all__ = ["StartingDirectionSpace", "fold_moments_to_starting_direction"]

#: The two starting-direction legs, in flat-layout (and DAG) order:
#: the inward leg first (seed⁻ ≺ seed⁺ — the outward leg is
#: pole-continued from the inward one, ψ½⁺(0) = ψ½⁻(0)).
_SIGNS: tuple[int, int] = (-1, +1)


@dataclass(frozen=True)
class StartingDirectionSpace(FunctionSpace):
    r"""Flat per-level ψ½ space with typed ``(level, sign)`` views and zero metric.

    Parameters
    ----------
    name, shape, inner_product_weights
        Inherited from :class:`~orpheus.numerics.space.FunctionSpace`.
        ``name`` is ``"starting_direction"``; ``shape`` is the flat
        total ``(n_levels · 2 · (ng·nx + ng),)``; the weights are
        **all-zero** (the ghost metric — see the module docstring).
        Build via :meth:`for_levels`, never the bare constructor.
    levels : tuple[int, ...]
        The seed-carrying μ-level indices (ascending), as selected by
        the R12a predicate on the owning mesh. Metadata
        (``compare=False`` — space identity stays ``(name, shape)``).
    ng, nx : int
        Group and radial-cell counts — the per-``(level, sign)`` block
        shapes ``(ng, nx)`` (cells) and ``(ng,)`` (corner). Metadata.
    """

    levels: tuple[int, ...] = field(kw_only=True, compare=False, repr=False)
    ng: int = field(kw_only=True, compare=False, repr=False)
    nx: int = field(kw_only=True, compare=False, repr=False)

    # ── Equality / hashing inherited from FunctionSpace ───────────────
    #
    # Same rationale as AngularTraceSpace / TensorProductSpace / FullFieldSpace:
    # restore the (name, shape) identity convention over the dataclass
    # auto-generation (the metadata fields are already ``compare=False``).

    def __eq__(self, other: object) -> bool:
        return FunctionSpace.__eq__(self, other)

    def __hash__(self) -> int:
        return FunctionSpace.__hash__(self)

    def __repr__(self) -> str:
        return (
            f"StartingDirectionSpace(levels={self.levels}, ng={self.ng}, "
            f"nx={self.nx}, shape={self.shape})"
        )

    # ── Construction ──────────────────────────────────────────────────

    @classmethod
    def for_levels(
        cls, levels: tuple[int, ...], *, ng: int, nx: int,
    ) -> "StartingDirectionSpace":
        r"""Build the space for the given seed-carrying levels.

        The ONLY construction path (the bare constructor cannot derive
        the flat shape or the zero weights). An EMPTY ``levels`` is
        rejected: absence of the block is spelled ``None`` at the
        carrier/mesh layer, never a zero-DOF phantom space.

        Parameters
        ----------
        levels : tuple[int, ...]
            Seed-carrying μ-level indices, strictly ascending (the
            canonical layout order; also the DAG order).
        ng, nx : int
            Group / radial-cell counts of the owning mesh.
        """
        levels = tuple(int(p) for p in levels)
        if not levels:
            raise ValueError(
                "StartingDirectionSpace.for_levels: levels is empty — a mesh "
                "with no seed-carrying levels has NO starting-direction "
                "space (spelled None), not an empty one."
            )
        if any(b <= a for a, b in zip(levels, levels[1:])):
            raise ValueError(
                f"StartingDirectionSpace.for_levels: levels must be strictly "
                f"ascending (the canonical layout order); got {levels!r}."
            )
        if ng <= 0 or nx <= 0:
            raise ValueError(
                f"StartingDirectionSpace.for_levels: ng and nx must be "
                f"positive; got ng={ng}, nx={nx}."
            )
        total = len(levels) * 2 * (ng * nx + ng)
        return cls(
            name="starting_direction",
            shape=(total,),
            # The ghost metric: every seed DOF carries zero measure —
            # (1-μ²)|_{μ=±1} = 0, the α_{1/2} = 0 fact (module docstring).
            inner_product_weights=np.zeros(total),
            levels=levels,
            ng=int(ng),
            nx=int(nx),
        )

    # ── Layout arithmetic (single source of the flat offsets) ─────────

    @property
    def n_levels(self) -> int:
        r"""Number of seed-carrying levels."""
        return len(self.levels)

    @property
    def _per_sign(self) -> int:
        r"""Flat size of one ``(level, sign)`` leg: cells ⊕ corner."""
        return self.ng * self.nx + self.ng

    def _leg_offset(self, level: int, sign: int) -> int:
        r"""Flat offset of the ``(level, sign)`` leg (cells first, corner after).

        Raises
        ------
        KeyError
            If ``level`` is not a seed-carrying level of this space.
        ValueError
            If ``sign`` is not ``-1`` or ``+1``.
        """
        if sign not in _SIGNS:
            raise ValueError(
                f"StartingDirectionSpace: sign must be -1 (inward) or +1 "
                f"(outward); got {sign!r}."
            )
        try:
            pos = self.levels.index(level)
        except ValueError:
            raise KeyError(
                f"StartingDirectionSpace: level {level!r} carries no "
                f"starting-direction block; seed-carrying levels are "
                f"{self.levels!r} (R12a predicate)."
            ) from None
        sign_pos = _SIGNS.index(sign)
        return (2 * pos + sign_pos) * self._per_sign

    def cells_slice(self, level: int, sign: int) -> slice:
        r"""Flat slice of the ``(level, sign)`` cells block (``ng·nx`` entries)."""
        base = self._leg_offset(level, sign)
        return slice(base, base + self.ng * self.nx)

    def corner_slice(self, level: int, sign: int) -> slice:
        r"""Flat slice of the ``(level, sign)`` corner slot (``ng`` entries)."""
        base = self._leg_offset(level, sign) + self.ng * self.nx
        return slice(base, base + self.ng)

    # ── Shaped views (no copies) ──────────────────────────────────────

    def cells_view(self, buffer: NDArray, level: int, sign: int) -> NDArray:
        r"""The ``(ng, nx)`` ψ½ cells view of ``buffer`` for ``(level, sign)``.

        A reshaped slice view — shares memory with ``buffer``.
        """
        return buffer[self.cells_slice(level, sign)].reshape(self.ng, self.nx)

    def corner_view(self, buffer: NDArray, level: int, sign: int) -> NDArray:
        r"""The ``(ng,)`` r = R corner view of ``buffer`` for ``(level, sign)``.

        Inflow corner (``sign = -1``): the given-data / identity row.
        Outflow corner (``sign = +1``): the defect row (ruling R13).
        Shares memory with ``buffer``.
        """
        return buffer[self.corner_slice(level, sign)]


def fold_moments_to_starting_direction(
    moments: np.ndarray,
    sign: int,
) -> np.ndarray:
    r"""Fold angular source moments to a starting direction: :math:`\bar Q(\mu=\pm 1)`.

    The Hébert (3.432) Legendre fold evaluated at the closed rays
    :math:`\mu = \pm 1` (ruling R14 — the FULL :math:`(-1)^\ell` fold
    from day one):

    .. math::

        \bar Q(\mu = \pm 1)
        \;=\; \sum_{\ell} \frac{2\ell + 1}{2}\, Q_\ell \, P_\ell(\pm 1)
        \;=\; \sum_{\ell} \frac{2\ell + 1}{2}\, Q_\ell \, (\pm 1)^\ell ,

    the exact 1-D addition-theorem weight :math:`(2\ell+1)/2` with
    :math:`P_\ell(\pm 1) = (\pm 1)^\ell`. The single source of the q½
    source construction for the starting-direction carrier blocks
    (#282 route (a)): the seed-arm emitters (the S composite arm, the
    F arm, the q_ext factories) all fold through HERE — the
    :math:`P_1(-1)` sign is spelled ONCE (vv Mode 1/6; the §16.B B2b
    2-term pin is live on this function).

    At :math:`\ell = 0` this reduces to :math:`\bar Q = \tfrac12 Q_0`
    — the same convention :func:`carlson_inward_sweep_from_source`'s
    ``Q_bar`` parameter documents (§16.B B2a). The helper accepts any
    number of moments; the CURRENT production emitters feed
    :math:`\ell = 0` only (the curvilinear operator's reach — the
    anisotropic seed fold activates with the >linear-in-μ companion
    gate), so :math:`\ell \ge 1` is manufactured-before-needed per
    §0.6's isotropic-snapshot-blindness discipline.

    Parameters
    ----------
    moments : np.ndarray, shape ``(n_moments, ...)``
        The source moments :math:`Q_\ell`, ℓ-leading (``moments[l]`` is
        :math:`Q_\ell`; any trailing shape — typically ``(ng, nx)``).
    sign : int
        ``-1`` (the inward starting direction) or ``+1`` (the
        pole-continued outward leg, R13).

    Returns
    -------
    np.ndarray
        :math:`\bar Q(\mu = \mathrm{sign})`, the trailing shape of
        ``moments``.
    """
    if sign not in (-1, +1):
        raise ValueError(
            f"fold_moments_to_starting_direction: sign must be -1 or +1; "
            f"got {sign!r}."
        )
    moments = np.asarray(moments)
    if moments.ndim < 1 or moments.shape[0] < 1:
        raise ValueError(
            f"fold_moments_to_starting_direction: moments must carry a "
            f"leading ℓ axis with at least the ℓ = 0 moment; got shape "
            f"{moments.shape!r}."
        )
    n_moments = moments.shape[0]
    ell = np.arange(n_moments)
    # (2ℓ+1)/2 · sign^ℓ, broadcast down the trailing axes.
    coeff = ((2.0 * ell + 1.0) / 2.0) * np.float64(sign) ** ell
    return np.tensordot(coeff, moments, axes=(0, 0))
