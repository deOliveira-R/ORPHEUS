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

Presence is keyed PER LEVEL by two structural facts about the level's
march-start edge
(:class:`~orpheus.sn.angular.closure.MarchStart`, Q5.4/T26):
a level carries a block **iff the M-M half-angle recurrence genuinely
consumes independent starting-direction state**, i.e. NEITHER

* ``on_edge_node`` — the start edge IS an ordinate (cylinder *product*
  NODE_ALIGNED rules: the η-minimum node lies on :math:`\Sigma`,
  :math:`\eta_0 = \eta_{1/2} = -\sin\theta` bit-exactly, the #229
  fact). The seed would be a rank-duplicate of :math:`\psi_0` — NO
  block; NOR
* ``degenerate`` — the η-minimum is shared (cylinder *level-symmetric*
  rules: duplicate-η hemisphere partners collapse the midpoint edge
  onto :math:`\eta_0`), so the seed's only consumption path — the
  recurrence weight :math:`(1-\tau_0)` — vanishes. The seed would be
  dead state — NO block. (This is why the measured cylinder-LS seed
  sensitivity is 0.0-bit.)

Sphere Gauss–Legendre is the carrying instance (a genuine off-node
start), as is every level of a σ_y-folded product rule (the arc, T22b).
The former encoding — the raw first-ordinate M-M weight ``τ_raw,0 ∈
(0, 1)`` exclusive — is a bit-exact gated *consequence* of the facts
(``on_edge_node ⟹ 0``, ``degenerate ⟹ 1``, neither ⟹ strict interior),
no longer the predicate itself.

R12a refines the R12 letter ("μ_start ∉ the level's μ-nodes"), whose
claimed equivalence to ``τ_raw ≠ 0`` is empirically false on
level-symmetric cylinder rules (μ_start ∉ nodes there, yet the seed is
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
space layouts, keyed by the structured
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
  (the cylinder azimuthal floor — its structural content is the level's
  edge-inclusion, a property of the CIRCLE, not of the
  :math:`[\tfrac12, 1]` absorber retired at Q5.6.4).
* ``.claude/plans/archive/stencil_assembly_dsa_roadmap.md`` — rulings R12/R12a,
  R13; ``a3_solve_transpose_verification.md`` §16.A (the carrier gates).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import ClassVar, NamedTuple, Self

import numpy as np
from numpy.typing import NDArray

from orpheus.numerics.face_layout import FaceLayout
from orpheus.numerics.space import FunctionSpace

__all__ = [
    "RadialCharacteristicInteriorSpace",
    "RadialCharacteristicBoundarySpace",
    "fold_moments_to_radial_characteristic",
    "fold_moments_to_radial_characteristic_transpose",
]

#: The two starting-direction legs, in flat-layout (and DAG) order:
#: the inward leg first (seed⁻ ≺ seed⁺ — the outward leg is
#: pole-continued from the inward one, ψ½⁺(0) = ψ½⁻(0)).
_SIGNS: tuple[int, int] = (-1, +1)


def _validate_ray_space_inputs(
    cls_name: str,
    levels: tuple[int, ...],
    ng: int,
    nx: int,
    cell_volumes: NDArray,
) -> tuple[tuple[int, ...], NDArray]:
    r"""The shared ``for_levels`` precondition check for every ψ½ ray space.

    The two split
    :class:`RadialCharacteristicInteriorSpace` / :class:`RadialCharacteristicBoundarySpace`
    share ONE construction contract: non-empty strictly-ascending seed levels,
    positive ``ng`` / ``nx``, and a strictly-positive ``cell_volumes`` of shape
    ``(nx,)`` (the SPD ``G_sd = V_cell`` state metric — a zero weight is the
    forbidden ghost metric that severs the seed from the G-adjoint, ERR-067). The
    check is single-sourced here (``cls_name`` selects the message prefix) so the
    three spaces cannot validate differently. Returns the cleaned ``(levels,
    cell_volumes)`` (int-tuple levels; float ``cell_volumes``).
    """
    levels = tuple(int(p) for p in levels)
    if not levels:
        raise ValueError(
            f"{cls_name}.for_levels: levels is empty — a mesh "
            f"with no seed-carrying levels has NO starting-direction "
            f"space (spelled None), not an empty one."
        )
    if any(b <= a for a, b in zip(levels, levels[1:])):
        raise ValueError(
            f"{cls_name}.for_levels: levels must be strictly "
            f"ascending (the canonical layout order); got {levels!r}."
        )
    if ng <= 0 or nx <= 0:
        raise ValueError(
            f"{cls_name}.for_levels: ng and nx must be "
            f"positive; got ng={ng}, nx={nx}."
        )
    cell_volumes = np.asarray(cell_volumes, dtype=float)
    if cell_volumes.shape != (nx,):
        raise ValueError(
            f"{cls_name}.for_levels: cell_volumes must have "
            f"shape ({nx},) to match nx; got {cell_volumes.shape}."
        )
    if np.any(cell_volumes <= 0.0):
        raise ValueError(
            f"{cls_name}.for_levels: cell_volumes must be "
            f"strictly positive — the state metric G_sd = V_cell must be "
            f"SPD (a zero weight is the forbidden ghost metric that severs "
            f"the seed from the G-adjoint; see the module docstring)."
        )
    return levels, cell_volumes


class _RayLeg(NamedTuple):
    r"""One shared leg of the ψ½ layout walk (see :func:`_radial_characteristic_legs`)."""

    #: ``"cells"`` (interior) or ``"corner"`` (boundary).
    part: str
    #: ``(level, sign)`` — the split interior / boundary spaces' layout key.
    split_key: tuple[int, int]
    #: ``(level, sign, part)`` — the unified space's 3-tuple layout key.
    #: The per-leg block shape (``(ng, nx)`` cells / ``(ng,)`` corner).
    shape: tuple[int, ...]
    #: The ``G_sd = V_cell`` state-metric weights for this leg.
    metric: NDArray


def _radial_characteristic_legs(
    levels: tuple[int, ...],
    ng: int,
    nx: int,
    cell_volumes: NDArray,
) -> list[_RayLeg]:
    r"""The single-source ``(level, sign)`` leg walk shared by every ψ½ ray space.

    Emits, in the canonical flat-layout (= DAG) order — ``level`` ascending,
    ``sign ∈ (-1, +1)``, ``cells`` then ``corner`` — one :data:`_RayLeg` per leg.
    Each split space filters on its ``part`` (keying
    by ``split_key``). Because all three spaces read their layout order AND their
    ``G_sd = V_cell`` state metric from THIS one walk, a split space's
    layout/metric cannot drift from the unified (Pattern 2 — the twin-path is
    unspellable; the split-fidelity gate is then a safety net, not the only
    guard). Per leg: ``V_cell`` per group over the ``(ng, nx)`` cells, and
    ``V(r = R) = cell_volumes[-1]`` over the ``(ng,)`` corner (the free-gauge
    weight — see the module docstring's "state metric" section for why the corner
    gauge is the outer cell volume and why ``w`` is dropped).
    """
    corner_gauge = float(cell_volumes[-1])
    legs: list[_RayLeg] = []
    for level in levels:
        for sign in _SIGNS:
            legs.append(
                _RayLeg(
                    part="cells",
                    split_key=(level, sign),
                    shape=(ng, nx),
                    metric=np.tile(cell_volumes, ng),
                )
            )
            legs.append(
                _RayLeg(
                    part="corner",
                    split_key=(level, sign),
                    shape=(ng,),
                    metric=np.full(ng, corner_gauge),
                )
            )
    return legs


# The historical UNIFIED ``RadialCharacteristicSpace`` (one flat buffer
# keyed by the 3-tuple ``(level, sign, part)``) RETIRED at Phase C 4e —
# the fused walk marches System B's split composite natively, so the split
# interior / boundary spaces below are the ONLY ψ½ layouts. The shared
# ``(level, sign)`` leg walk + metric (``_radial_characteristic_legs``)
# and validation (``_validate_ray_space_inputs``) single-source both.


# ═══════════════════════════════════════════════════════════════════════
# The SPLIT ψ½ spaces — System B's interior ⊕ boundary (Phase B)
# ═══════════════════════════════════════════════════════════════════════


@dataclass(frozen=True)
class _RadialCharacteristicSubSpace(FunctionSpace):
    r"""Shared base of the two SPLIT ψ½ spaces (interior cells / boundary corner).

    Phase B (the coupled-block campaign) poses the ψ½ ray as **System B** — its
    OWN ``interior ⊕ boundary`` composite — by splitting the unified
    the retired unified space (keyed by the 3-tuple ``(level, sign,
    part)``) along its ``part`` axis into two per-``(level, sign)`` spaces:

    * :class:`RadialCharacteristicInteriorSpace` — the ``cells`` legs
      (``(ng, nx)`` ψ½ at every radial cell), A_BB's marched interior state,
      metric ``G_sd = V_cell``;
    * :class:`RadialCharacteristicBoundarySpace` — the ``corner`` legs
      (``(ng,)`` ψ½ at r = R), B_b's BC locus (inflow = data / outflow = defect),
      metric ``G = V(r = R)``.

    Both are ``FaceLayout[tuple[int, int]]`` spaces keyed by ``(level, sign)``,
    differing ONLY in their per-leg block and their name; each concrete subclass
    fixes those via the ``_PART`` / ``_SPACE_NAME`` ClassVars. The ``(level,
    sign)`` leg order AND the metric are single-sourced from
    :func:`_radial_characteristic_legs` (shared with the unified space), so a
    split space cannot drift from the unified — the split-fidelity gate is a
    safety net, not the only guard. Space identity stays the ``(name, shape)``
    tuple (the metadata fields are ``compare=False``).
    """

    #: The ``part`` this space projects the shared leg walk onto — ``"cells"``
    #: (interior) or ``"corner"`` (boundary). Set by each concrete subclass; the
    #: base is abstract (a bare ``_RadialCharacteristicSubSpace`` has no ``_PART``).
    _PART: ClassVar[str]
    #: The ``(name, shape)`` identity name. Set by each concrete subclass.
    _SPACE_NAME: ClassVar[str]

    levels: tuple[int, ...] = field(kw_only=True, compare=False, repr=False)
    ng: int = field(kw_only=True, compare=False, repr=False)
    nx: int = field(kw_only=True, compare=False, repr=False)
    #: The flat-buffer :class:`~orpheus.numerics.face_layout.FaceLayout` keyed by
    #: ``(level, sign)`` (the ``part`` is fixed by ``_PART``). Metadata
    #: (``compare=False`` — identity stays ``(name, shape)``); built by
    #: :meth:`for_levels`.
    layout: FaceLayout[tuple[int, int]] = field(
        kw_only=True, compare=False, repr=False,
    )

    # ── Equality / hashing inherited from FunctionSpace ───────────────
    def __eq__(self, other: object) -> bool:
        return FunctionSpace.__eq__(self, other)

    def __hash__(self) -> int:
        return FunctionSpace.__hash__(self)

    def __repr__(self) -> str:
        return (
            f"{type(self).__name__}(levels={self.levels}, ng={self.ng}, "
            f"nx={self.nx}, shape={self.shape})"
        )

    # ── Construction (a projection of the unified for_levels) ─────────
    @classmethod
    def for_levels(
        cls,
        levels: tuple[int, ...],
        *,
        ng: int,
        nx: int,
        cell_volumes: NDArray,
    ) -> Self:
        r"""Build the split space for the given seed-carrying levels.

        The SAME contract + inputs as
        the retired unified ``for_levels`` (the split spaces project the same
        leg walk): validated by the shared
        :func:`_validate_ray_space_inputs`, then the shared
        :func:`_radial_characteristic_legs` walk FILTERED to this space's
        ``_PART``. ``nx`` / ``cell_volumes`` are required by BOTH spaces even
        though the boundary corner is ``(ng,)`` — the corner's ``V(r = R)`` gauge
        is ``cell_volumes[-1]``, and the shared walk / validation are uniform
        across the three spaces.
        """
        levels, cell_volumes = _validate_ray_space_inputs(
            cls.__name__, levels, ng, nx, cell_volumes,
        )
        named_shapes: list[tuple[tuple[int, int], tuple[int, ...]]] = []
        metric_pieces: list[NDArray] = []
        for leg in _radial_characteristic_legs(levels, ng, nx, cell_volumes):
            if leg.part != cls._PART:
                continue
            named_shapes.append((leg.split_key, leg.shape))
            metric_pieces.append(leg.metric)
        layout = FaceLayout.from_named_shapes(named_shapes)
        # CS4b S3 (F2): the name carries a CONTENT digest — the level set,
        # the split layout, and the ray metric — so the inherited
        # ``(name, shape)`` equality IS content equality (mirrors the
        # trace-space mints). Two carriers with the same ray structure mint
        # EQUAL spaces; a different level set, grid, or metric refuses.
        import hashlib

        metric = np.concatenate(metric_pieces)
        payload = b"".join((
            repr([
                (str(k), int(s.offset), int(s.flat_size))
                for k, s in layout.faces.items()
            ]).encode(),
            repr(tuple(int(lv) for lv in levels)).encode(),
            metric.tobytes(),
        ))
        digest = hashlib.blake2b(payload, digest_size=8).hexdigest()
        return cls(
            name=f"{cls._SPACE_NAME}#{digest}",
            shape=(layout.total_size,),
            inner_product_weights=metric,
            levels=levels,
            ng=int(ng),
            nx=int(nx),
            layout=layout,
        )

    # ── Layout access ─────────────────────────────────────────────────
    @property
    def n_levels(self) -> int:
        r"""Number of seed-carrying levels."""
        return len(self.levels)

    def _slot_key(self, level: int, sign: int) -> tuple[int, int]:
        r"""Validate ``(level, sign)`` and return the layout key.

        Reproduces the retired unified ``_slot_key``'s error contract,
        minus the ``part`` (fixed per space by ``_PART``).
        """
        if sign not in _SIGNS:
            raise ValueError(
                f"{type(self).__name__}: sign must be -1 (inward) or +1 "
                f"(outward); got {sign!r}."
            )
        if level not in self.levels:
            raise KeyError(
                f"{type(self).__name__}: level {level!r} carries no "
                f"starting-direction block; seed-carrying levels are "
                f"{self.levels!r} (R12a predicate)."
            )
        return (level, sign)

    def slot_view(self, buffer: NDArray, level: int, sign: int) -> NDArray:
        r"""The ``(level, sign)`` block view of ``buffer`` — memory-shared.

        A reshaped slice view through the :attr:`layout` slot (interior:
        ``(ng, nx)`` cells; boundary: ``(ng,)`` corner). Shares memory with
        ``buffer``. The field leaves expose it under their own name
        (``cells`` / ``corner``).
        """
        slot = self.layout.faces[self._slot_key(level, sign)]
        return slot.slice_view(buffer)


@dataclass(frozen=True)
class RadialCharacteristicInteriorSpace(_RadialCharacteristicSubSpace):
    r"""The ψ½ interior (``cells``) space — System B's interior block.

    The ``(ng, nx)`` half-angle flux at every radial cell, per seed-carrying
    ``(level, sign)`` leg, under the SPD ``G_sd = V_cell`` state metric.
    :class:`~orpheus.sn.operators.radial_characteristic.RadialCharacteristicOperator`
    (A_BB) marches this block (μ∂_r + σ_t); A_AB reads its inward leg; A_BA
    writes it. See :class:`_RadialCharacteristicSubSpace`.
    """

    _PART: ClassVar[str] = "cells"
    _SPACE_NAME: ClassVar[str] = "radial_characteristic_interior"


@dataclass(frozen=True)
class RadialCharacteristicBoundarySpace(_RadialCharacteristicSubSpace):
    r"""The ψ½ boundary (``corner``) space — System B's boundary block.

    The ``(ng,)`` half-angle flux at r = R, per seed-carrying ``(level, sign)``
    leg, under the ``G = V(r = R)`` corner gauge. Inflow corner (``sign = -1``) =
    given BC data; outflow corner (``sign = +1``) = the defect row (ruling R13).
    :class:`~orpheus.sn.operators.boundary.RadialCharacteristicBoundaryOperator`
    (B_b) acts on this block. See :class:`_RadialCharacteristicSubSpace`.
    """

    _PART: ClassVar[str] = "corner"
    _SPACE_NAME: ClassVar[str] = "radial_characteristic_boundary"


def _radial_characteristic_reconstruction_weights(
    n_moments: int,
    sign: int,
) -> np.ndarray:
    r"""The 1-D Legendre reconstruction weights at the closed ray :math:`\mu=\pm 1`.

    .. math::

        w_\ell \;=\; \frac{2\ell + 1}{2}\,(\pm 1)^\ell ,
        \qquad \ell = 0 \ldots n_{\rm moments} - 1,

    the exact addition-theorem weight :math:`(2\ell+1)/2` with
    :math:`P_\ell(\pm 1) = (\pm 1)^\ell` (Hébert 3.432). The **single source**
    of the fold coefficient: both the forward reconstruction
    :func:`fold_moments_to_radial_characteristic` (which contracts these
    weights against the moments) and its Euclidean transpose
    :func:`fold_moments_to_radial_characteristic_transpose` (which expands a
    ray cotangent onto them) read this — so the :math:`P_1(-1)` sign is
    spelled ONCE (vv Mode 1/6; the §16.B B2b 2-term pin lives on the forward).
    """
    if sign not in (-1, +1):
        raise ValueError(
            f"_radial_characteristic_reconstruction_weights: sign must be "
            f"-1 or +1; got {sign!r}."
        )
    if n_moments < 1:
        raise ValueError(
            f"_radial_characteristic_reconstruction_weights: n_moments must "
            f"be ≥ 1 (at least the ℓ = 0 moment); got {n_moments!r}."
        )
    ell = np.arange(n_moments)
    # (2ℓ+1)/2 · sign^ℓ.
    return ((2.0 * ell + 1.0) / 2.0) * np.float64(sign) ** ell


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
    moments = np.asarray(moments)
    if moments.ndim < 1 or moments.shape[0] < 1:
        raise ValueError(
            f"fold_moments_to_radial_characteristic: moments must carry a "
            f"leading ℓ axis with at least the ℓ = 0 moment; got shape "
            f"{moments.shape!r}."
        )
    # (2ℓ+1)/2 · sign^ℓ from the single-source weights, contracted down the
    # ℓ axis (the sign is validated inside the weights helper).
    weights = _radial_characteristic_reconstruction_weights(
        moments.shape[0], sign,
    )
    return np.tensordot(weights, moments, axes=(0, 0))


def fold_moments_to_radial_characteristic_transpose(
    rays_bar: np.ndarray,
    sign: int,
    n_moments: int,
) -> np.ndarray:
    r"""Euclidean transpose of :func:`fold_moments_to_radial_characteristic`.

    The forward reconstruction contracts the ray value out of the moments,
    :math:`\bar Q(\mu=\mathrm{sign}) = \sum_\ell w_\ell\,Q_\ell`; its adjoint
    with respect to the moments therefore EXPANDS a cotangent on that ray
    value back onto the moments,

    .. math::

        \bar Q_\ell \;=\; w_\ell(\mathrm{sign})\,\overline{\bar Q},
        \qquad \ell = 0 \ldots n_{\rm moments} - 1,

    with the SAME reconstruction weights
    :math:`w_\ell = \tfrac{2\ell+1}{2}(\pm 1)^\ell`
    (:func:`_radial_characteristic_reconstruction_weights` — single source).
    This is the transpose that single-sources the scattering seed adjoint
    (:math:`\partial S / \partial \psi_{1/2}` cotangent → bulk moment): it
    retires the hand-rolled :math:`\ell = 0` factor :math:`\tfrac12` that arm
    previously carried, keeping the sign/weight in ONE place (vv Mode 1).

    Parameters
    ----------
    rays_bar : np.ndarray, shape ``(...)``
        The cotangent on :math:`\bar Q(\mu = \mathrm{sign})` — typically
        ``(ng, nx)`` (a per-level, per-sign ray-cells cotangent).
    sign : int
        ``-1`` (inward starting direction) or ``+1`` (outward, pole-continued).
    n_moments : int
        The number of moments to expand onto — the reconstruction operator's
        domain dimension (≥ 1). For the isotropic production reach this is 1.

    Returns
    -------
    np.ndarray, shape ``(n_moments, *rays_bar.shape)``
        :math:`\bar Q_\ell`, the moment-space cotangent (ℓ-leading).
    """
    rays_bar = np.asarray(rays_bar)
    weights = _radial_characteristic_reconstruction_weights(n_moments, sign)
    # outer product: (n_moments,) ⊗ (…) → (n_moments, …), the transpose of
    # the forward tensordot's ℓ-contraction.
    return weights.reshape((n_moments,) + (1,) * rays_bar.ndim) * rays_bar[None]
