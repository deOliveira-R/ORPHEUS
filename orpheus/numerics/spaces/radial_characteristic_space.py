r"""The starting-direction space — per-level ψ½ carrier with the V_cell state metric.

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
:attr:`~orpheus.sn.mesh.augmented_mesh.SNMesh.radial_characteristic_levels`;
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
views into the flat buffer, no copies. The offsets are sourced from a
real :class:`~orpheus.numerics.face_layout.FaceLayout` carried on the
space (:attr:`~RadialCharacteristicSpace.layout`), keyed by the structured
``(level, sign, part)`` tuple — the SAME flat-buffer discipline the
spatial trace uses with ``str`` face-name keys (``FaceLayout[str]``), the
key merely typed instead of stringly.

The state metric — radial cell volume :math:`V_{\rm cell}` (SPD)
================================================================

The inner-product weight of this space is the **state metric**
:math:`G_{\rm sd} = V_{\rm cell}` — the radial cell-volume measure,
strictly positive and therefore symmetric positive-definite (SPD),
group-broadcast across the cells and gauge-extended to the corners. It
mirrors the bulk metric :math:`G_{\rm bulk} = V_{\rm cell}\,w_n`: the
SAME spatial measure, restricted to the single :math:`\mu = \pm 1` ray
(the angular factor :math:`w` is dropped — see the gauge argument).

**ψ½ is a first-class radial STATE field, not a face trace.** Its
operator self-block :math:`A_{\rm ss}` is a *banded radial transport
operator* :math:`\mu\,\partial_r + \sigma_t` (Hébert Eqs. 3.434–3.435;
measured off-diagonal norm :math:`\approx 71` on a 6-cell sphere — a
genuine interior radial dynamics, unlike the spatial trace's pure
restriction map). **A state's Hilbert metric is set by its OPERATOR
ROLE, not by its integration weight** — so ψ½'s metric is the bulk's
spatial measure, made nonzero by the radial-field volume.

Three pole-vanishing quantities, historically conflated — keep them apart:

* **M1 — the moment / output weight** :math:`= 0`. The *open*
  Gauss–Legendre rule has no node at the pole, so ψ½ carries zero weight
  in the scalar-flux moment :math:`\phi = \sum_n w_n\psi_n`. This
  correctly *excludes* ψ½ from the flux; it lives in the moment reducer,
  not in this space.
* **M2 — the angular through-flux coefficient**
  :math:`(1-\mu^2)\big|_{\mu=\pm 1} = 0` (the α-dome endpoints
  :math:`\alpha_{1/2} = \alpha_{N+1/2} = 0`). This is an **operator
  coefficient** *inside* :math:`A` — the strength of the
  angular-derivative redistribution term — governing the
  straight-characteristic / triangular structure. Correctly zero, and it
  lives in the operator, not in the state inner product.
* **M3 — the state metric** :math:`G_{\rm sd} = V_{\rm cell} \neq 0`.
  *This space's* inner product; it governs the G-adjoint reciprocity
  :math:`\langle A\psi,\chi\rangle_G = \langle\psi, A^\dagger\chi\rangle_G`.

The retired **"ghost metric" bug** installed **M2** (an operator
coefficient) as **M3** (the state metric): it read the angular
through-flux weight :math:`(1-\mu^2)|_{\rm pole} = 0` as the Hilbert
metric and set :math:`G_{\rm sd} \equiv 0`.

**Why 0 is the one forbidden value, and why** :math:`V_{\rm cell}`.
Because :meth:`~orpheus.numerics.operator.SupportsAdjoint.apply_transpose`
is the *exact* Euclidean transpose (:math:`T = A^{\mathsf T}`, measured
:math:`\lVert T - A^{\mathsf T}\rVert = 3.6\times10^{-16}`), the
determining relation :math:`A^{\mathsf T}G = G A^\dagger` behind
:math:`A^\dagger = G^{-1}A^{\mathsf T}G` holds for **every** SPD
:math:`G_{\rm sd}` — the reciprocity is gauge-free among SPD choices, so
the metric is fixed only up to an SPD gauge. The single **forbidden**
value is :math:`0`: it puts the seed rows in :math:`\ker G`, which
severs the seed → bulk coupling :math:`A_{\rm bs}` from
:math:`A^\dagger = G^{-1}A^{\mathsf T}G` (the adjoint's seed block goes
identically zero). The shipped :math:`G_{\rm sd} = 0` was therefore a
**wrong adjoint the instant the seed carries data** — a production
reciprocity defect of :math:`1.3\times10^{-2}`, green only because the
gate fed a present-but-zero seed. We gauge-fix to
:math:`G_{\rm sd} = V_{\rm cell}` (dropping :math:`w` — a single
:math:`\mu=\pm 1` ray has no canonical quadrature weight) so the
adjoint's seed block is the **physical backward radial march**, and all
bulk/trace observables (:math:`\phi^\dagger`, adjoint reaction rates)
are **bitwise gauge-invariant** (the block-upper-triangular
:math:`A^\dagger` seats the seed at the top, so only
:math:`\phi^\dagger_{\rm seed}` moves with the gauge). Derivation of
record:
``.claude/agent-memory/numerics-investigator/radial_characteristic_metric_gauge_derivation.md``.

Consequences (inherited from
:class:`~orpheus.numerics.space.FunctionSpace` with the SPD weight
vector — no overrides):
:meth:`~orpheus.numerics.space.FunctionSpace.apply_metric` **scales** the
block by :math:`V_{\rm cell}`;
:meth:`~orpheus.numerics.space.FunctionSpace.apply_inverse_metric`
**divides** by it (no masking — the null space is empty); and
:meth:`~orpheus.numerics.space.FunctionSpace.inner_product` contributes
:math:`\sum V_{\rm cell}\,x\,y`.

**vv-principles Mode 12 — CLOSED (ERR-067).** Under the ghost
:math:`G_{\rm sd} = 0` the seed rows lay inside the G-adjoint
reciprocity functional's invariance group, so the G3 reciprocity gate
was *identically blind* to any seed-row transpose error (a false green,
at every tolerance, in every regime). With
:math:`G_{\rm sd} = V_{\rm cell}` the seed rows carry metric weight and
move **out** of that invariance group: a seed-row (:math:`A_{\rm ss}`)
sign flip now REDs G-reciprocity — the Mode-12-closure gate
:func:`tests.sn.sweep.curvilinear.test_282_direct_seed_fixed_point.test_mode12_g_reciprocity_catches_a_seed_row_flip`
— while the unmutated nonzero-seed reciprocity holds
:math:`< 10^{-12}`. The direct-solver closed-form pin (§16.B), the
solve∘apply residual over the full augmented field (§16.C), and 2.5b's
Euclidean :math:`M^{\mathsf T}` oracle pin the seed *coefficients* —
orthogonal to, and now complemented by, this metric-level reciprocity
catch.

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

from orpheus.numerics.face_layout import FaceLayout
from orpheus.numerics.space import FunctionSpace

__all__ = ["RadialCharacteristicSpace", "fold_moments_to_radial_characteristic"]

#: The two starting-direction legs, in flat-layout (and DAG) order:
#: the inward leg first (seed⁻ ≺ seed⁺ — the outward leg is
#: pole-continued from the inward one, ψ½⁺(0) = ψ½⁻(0)).
_SIGNS: tuple[int, int] = (-1, +1)


@dataclass(frozen=True)
class RadialCharacteristicSpace(FunctionSpace):
    r"""Flat per-level ψ½ space with typed ``(level, sign)`` views and the V_cell metric.

    Parameters
    ----------
    name, shape, inner_product_weights
        Inherited from :class:`~orpheus.numerics.space.FunctionSpace`.
        ``name`` is ``"radial_characteristic"``; ``shape`` is the flat
        total ``(n_levels · 2 · (ng·nx + ng),)``; the weights are the
        SPD **state metric** :math:`G_{\rm sd} = V_{\rm cell}` (see the
        module docstring's "The state metric" section).
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
    #: The flat-buffer :class:`~orpheus.numerics.face_layout.FaceLayout`
    #: keyed by the structured ``(level, sign, part)`` tuple — the single
    #: source of the per-``(level, sign)`` cells/corner offsets (the same
    #: discipline the spatial trace carries via ``str`` keys). Metadata
    #: (``compare=False`` — space identity stays ``(name, shape)``); built
    #: unconditionally by :meth:`for_levels`.
    layout: FaceLayout[tuple[int, int, str]] = field(
        kw_only=True, compare=False, repr=False,
    )

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
            f"RadialCharacteristicSpace(levels={self.levels}, ng={self.ng}, "
            f"nx={self.nx}, shape={self.shape})"
        )

    # ── Construction ──────────────────────────────────────────────────

    @classmethod
    def for_levels(
        cls,
        levels: tuple[int, ...],
        *,
        ng: int,
        nx: int,
        cell_volumes: NDArray,
    ) -> "RadialCharacteristicSpace":
        r"""Build the space for the given seed-carrying levels.

        The ONLY construction path (the bare constructor cannot derive
        the flat shape or the state metric). An EMPTY ``levels`` is
        rejected: absence of the block is spelled ``None`` at the
        carrier/mesh layer, never a zero-DOF phantom space.

        Parameters
        ----------
        levels : tuple[int, ...]
            Seed-carrying μ-level indices, strictly ascending (the
            canonical layout order; also the DAG order).
        ng, nx : int
            Group / radial-cell counts of the owning mesh.
        cell_volumes : NDArray, shape ``(nx,)``
            The radial cell volumes :math:`V_{\rm cell}` of the owning
            mesh — the state metric :math:`G_{\rm sd} = V_{\rm cell}` (see
            the module docstring "The state metric" section). The SAME
            array the bulk metric :math:`G_{\rm bulk} = V_{\rm cell}\,w_n`
            reads (``SNMesh.full_field_space``).
        """
        levels = tuple(int(p) for p in levels)
        if not levels:
            raise ValueError(
                "RadialCharacteristicSpace.for_levels: levels is empty — a mesh "
                "with no seed-carrying levels has NO starting-direction "
                "space (spelled None), not an empty one."
            )
        if any(b <= a for a, b in zip(levels, levels[1:])):
            raise ValueError(
                f"RadialCharacteristicSpace.for_levels: levels must be strictly "
                f"ascending (the canonical layout order); got {levels!r}."
            )
        if ng <= 0 or nx <= 0:
            raise ValueError(
                f"RadialCharacteristicSpace.for_levels: ng and nx must be "
                f"positive; got ng={ng}, nx={nx}."
            )
        cell_volumes = np.asarray(cell_volumes, dtype=float)
        if cell_volumes.shape != (nx,):
            raise ValueError(
                f"RadialCharacteristicSpace.for_levels: cell_volumes must have "
                f"shape ({nx},) to match nx; got {cell_volumes.shape}."
            )
        if np.any(cell_volumes <= 0.0):
            raise ValueError(
                "RadialCharacteristicSpace.for_levels: cell_volumes must be "
                "strictly positive — the state metric G_sd = V_cell must be "
                "SPD (a zero weight is the forbidden ghost metric that severs "
                "the seed from the G-adjoint; see the module docstring)."
            )
        # The STATE metric G_sd = V_cell — the radial cell-volume measure.
        # The pole ray ψ½ is a first-class radial STATE field (its self-block
        # A_ss is a banded radial transport operator μ∂_r + σ_t — Hébert
        # 3.434–3.435), so its Hilbert metric is the bulk's spatial measure
        # restricted to the ray: G_sd = V_cell, mirroring the bulk
        # G_bulk = V_cell·w_n. The angular factor w is a FREE GAUGE (a single
        # μ=±1 ray carries no quadrature weight); gauge-fixed to V_cell so
        # A.H's seed block is the physical backward radial march, and so
        # bulk/trace observables are bitwise gauge-invariant. This is NOT the
        # angular through-flux (1−μ²)|_pole = 0 — that is an OPERATOR
        # coefficient (it correctly excludes ψ½ from the moment φ = Σ w_n ψ_n),
        # not the state inner product. Setting G_sd = 0 severs the seed→bulk
        # A_bs coupling from A.H, making it a WRONG adjoint for any nonzero
        # seed (derivation: [[starting-direction-metric-gauge-derivation]]).
        # ONE walk over the (level, sign, part) legs emits BOTH the layout
        # shapes AND the state-metric weights, so the flat-buffer offsets and
        # the G_sd = V_cell weights cannot diverge — a single source of the leg
        # order (retiring the hand-rolled _leg_offset). Cells then corner within
        # each leg; legs in (level ascending, sign ∈ (-1, +1)) order. Per leg:
        # V_cell per group over the cells, and V(r = R) = cell_volumes[-1] on
        # the corner (the free-gauge weight).
        corner_gauge = float(cell_volumes[-1])
        named_shapes: list[tuple[tuple[int, int, str], tuple[int, ...]]] = []
        metric_pieces: list[NDArray] = []
        for level in levels:
            for sign in _SIGNS:
                named_shapes.append(((level, sign, "cells"), (ng, nx)))
                metric_pieces.append(np.tile(cell_volumes, ng))
                named_shapes.append(((level, sign, "corner"), (ng,)))
                metric_pieces.append(np.full(ng, corner_gauge))
        layout = FaceLayout.from_named_shapes(named_shapes)
        metric = np.concatenate(metric_pieces)
        return cls(
            name="radial_characteristic",
            shape=(layout.total_size,),
            inner_product_weights=metric,
            levels=levels,
            ng=int(ng),
            nx=int(nx),
            layout=layout,
        )

    # ── Layout access (offsets sourced from the FaceLayout) ───────────

    @property
    def n_levels(self) -> int:
        r"""Number of seed-carrying levels."""
        return len(self.levels)

    def _slot_key(
        self, level: int, sign: int, part: str,
    ) -> tuple[int, int, str]:
        r"""Validate ``(level, sign)`` and return the ``(level, sign, part)`` layout key.

        The flat offsets themselves live on :attr:`layout` (the
        single-source :class:`~orpheus.numerics.face_layout.FaceLayout`);
        this reproduces the retired ``_leg_offset``'s error contract so the
        view accessors reject the same misuse.

        Raises
        ------
        ValueError
            If ``sign`` is not ``-1`` or ``+1``.
        KeyError
            If ``level`` is not a seed-carrying level of this space.
        """
        if sign not in _SIGNS:
            raise ValueError(
                f"RadialCharacteristicSpace: sign must be -1 (inward) or +1 "
                f"(outward); got {sign!r}."
            )
        if level not in self.levels:
            raise KeyError(
                f"RadialCharacteristicSpace: level {level!r} carries no "
                f"starting-direction block; seed-carrying levels are "
                f"{self.levels!r} (R12a predicate)."
            )
        return (level, sign, part)

    # ── Shaped views (no copies — slice views into the flat buffer) ───

    def cells_view(self, buffer: NDArray, level: int, sign: int) -> NDArray:
        r"""The ``(ng, nx)`` ψ½ cells view of ``buffer`` for ``(level, sign)``.

        A reshaped slice view through the :attr:`layout` slot — shares
        memory with ``buffer``.
        """
        slot = self.layout.faces[self._slot_key(level, sign, "cells")]
        return slot.slice_view(buffer)

    def corner_view(self, buffer: NDArray, level: int, sign: int) -> NDArray:
        r"""The ``(ng,)`` r = R corner view of ``buffer`` for ``(level, sign)``.

        Inflow corner (``sign = -1``): the given-data / identity row.
        Outflow corner (``sign = +1``): the defect row (ruling R13).
        Shares memory with ``buffer``.
        """
        slot = self.layout.faces[self._slot_key(level, sign, "corner")]
        return slot.slice_view(buffer)


def fold_moments_to_radial_characteristic(
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
            f"fold_moments_to_radial_characteristic: sign must be -1 or +1; "
            f"got {sign!r}."
        )
    moments = np.asarray(moments)
    if moments.ndim < 1 or moments.shape[0] < 1:
        raise ValueError(
            f"fold_moments_to_radial_characteristic: moments must carry a "
            f"leading ℓ axis with at least the ℓ = 0 moment; got shape "
            f"{moments.shape!r}."
        )
    n_moments = moments.shape[0]
    ell = np.arange(n_moments)
    # (2ℓ+1)/2 · sign^ℓ, broadcast down the trailing axes.
    coeff = ((2.0 * ell + 1.0) / 2.0) * np.float64(sign) ** ell
    return np.tensordot(coeff, moments, axes=(0, 0))
