r"""Half-angle face flux seed strategies for the M-M angular recurrence.

Issue #168 Phase D — closes ERR-026 on sphere SN MMS by replacing the
hardcoded ``psi_half_left = 0`` seed in
:meth:`~orpheus.sn.spatial.pole_angular_closure.MorelMontryAngularSweep._weighted_angular_recurrence_single_level`
with the canonical Hébert §3.9.4 Carlson coupled-pole inward μ = −1
sweep.

The bug
=======

The Phase B :class:`~orpheus.sn.spatial.pole_angular_closure.MorelMontryAngularSweep`
recurrence runs Hébert Eqs. 3.437 / 3.439 forward across ordinates
:math:`m = 1, \ldots, N`:

.. math::

   \phi_{m+1/2,i,g} \;=\;
   \frac{\phi_{m,i,g} \;-\; (1 - \tau_m)\,\phi_{m-1/2,i,g}}{\tau_m}

The recurrence needs a seed :math:`\phi_{1/2,i,g}` at the auxiliary
inward starting direction μ = −1.  Phase B hardcoded
:math:`\phi_{1/2,i,g} = 0` (the "Carlson zero-weight starting
direction" convention with α\ :sub:`1/2` = 0 by Hébert Eq. 3.423),
reasoning that the term ``α_{1/2}·ψ_{1/2,i,g}`` in the redistribution
formula vanishes regardless of the seed value because
α\ :sub:`1/2` = 0 by GL antisymmetry.

That reasoning is WRONG.  The seed ``ψ_{1/2,i,g}`` enters the
recurrence's denominator through the propagation chain — every
subsequent half-angle face flux ``ψ_{m+1/2,i,g}`` depends on
``ψ_{m-1/2,i,g}`` recursively, and the chain inherits the seed
through the M-M weighting ``(1 - τ_m)``.  Setting the seed to zero
when the canonical Hébert §3.9.4 form says ``ψ_{1/2,i,g} = φ̄_{1/2,i}``
(the cell-centred output of the inward μ = −1 sweep) is a **wrong
term initialization** (Mode 3 per :doc:`/.claude/skills/numerical-bug-signatures/SKILL`).

This wrong seed survives Phase B's flat-ψ-identity test
(``tests/sn/spatial/test_pole_angular_closure.py``) because that
test only compares the three strategies against each other on flat
ψ — it does NOT compare against the closed-form fixed-point
identity ``L·ψ = Σ_t·ψ``.  Cylindrical Gate 1.1 MMS passes despite
the wrong seed because each μ-level's α-dome telescopes back to
zero at the level's edges, absorbing the wrong seed via cancellation.
Spherical Gate 1.1 MMS FAILS because the sphere cascade is more
rigid (no per-level telescoping) — the wrong seed propagates to a
wrong fixed point.  See ERR-026 in
:file:`.claude/skills/vv-principles/error_catalog.md`.

The fix
=======

Hébert §3.9.4 Eqs. (3.432)–(3.435) give the **inward μ = −1
"Carlson coupled-pole" sweep**.  At μ = −1 the angular redistribution
coefficient :math:`1 - \mu^2 = 0`, so the streaming–collision balance
decouples from the α-cascade and reduces to a plain DD inward
recurrence in radius:

.. math::
   :label: hebert-3-434

   \bar\phi_i \;=\; \frac{\Delta r_i \cdot \bar Q_i
                            + 2 \cdot \bar\phi_{i+1/2}}
                          {\Delta r_i \cdot \Sigma_i + 2}
   \quad\text{Hébert Eq. (3.434)}

.. math::
   :label: hebert-3-435

   \bar\phi_{i-1/2} \;=\; 2 \cdot \bar\phi_i - \bar\phi_{i+1/2}
   \quad\text{Hébert Eq. (3.435)}

The cell-averaged source at μ = −1 is the Legendre-moment-folded
source built from the input ψ's angular moments evaluated at μ = −1
via :math:`P_\ell(-1) = (-1)^\ell`:

.. math::
   :label: hebert-3-432-source

   \bar Q_i \;=\; \sum_\ell \frac{2\ell+1}{2}\,Q_\ell(r_i)\,(-1)^\ell

The boundary value :math:`\bar\phi_{nx+1/2}` is the angular flux at
the outer face at μ = −1.  For vacuum BC this is 0; for reflective
BC it is the outgoing-face value at μ = +1 (mirror reflection).

The sweep proceeds **INWARD** from ``i = nx-1`` down to ``i = 0``,
returning the cell-centred profile :math:`\bar\phi_i \equiv
\phi_{1/2,i}` per cell.  The resulting array seeds the M-M angular
recurrence's ``psi_half_left`` initialisation — replacing the
hardcoded zero with the structurally-correct Hébert value.

Two concrete strategies
=======================

* :class:`ZeroSeed` (key ``"zero"``) — reproduces Phase B's hardcoded
  ``psi_half_left = 0`` behaviour.  Used as the regression-safety
  ablation against the Carlson canonical form, and to A/B test the
  Phase D fix.
* :class:`CarlsonInwardSweep` (key ``"carlson_inward_sweep"``) —
  the canonical Hébert §3.9.4 (3.432)-(3.435) inward sweep.  Phase D
  default.

Architectural choice — Option α (composition, not sibling Protocol)
====================================================================

The seed is **M-M-specific**: the
:class:`~orpheus.sn.spatial.pole_angular_closure.LegacyTauSymmetricInterpolation`
and :class:`~orpheus.sn.spatial.pole_angular_closure.BaileyFlatFluxRedist`
closures do not have a ``psi_half_left`` variable to seed — their
half-angle face flux evaluation collapses to cell-centre values
unconditionally.  Composing the seed strategy as a field on
:class:`MorelMontryAngularSweep` keeps the abstraction local to where
the seed is consumed.  Phase D's Step 3a evaluation against an
alternative "sibling Protocol on ``SNMesh``" architecture (per the
original plan) showed that the sibling Protocol would force every
geometry consumer to handle an "irrelevant for non-M-M" Protocol,
violating SRP.

Linear-operator preservation
============================

Both strategies are LINEAR in the input ``psi_cells`` (verified by
:attr:`is_linear` ClassVar = True).  This is the load-bearing
property: the apply matvec must be a linear operator, otherwise
the operator-algebra capabilities of
:class:`~orpheus.sn.operator.SNStreamingOperator` (apply,
apply_transpose, dense matrix probing) break.  The
:class:`CarlsonInwardSweep` is linear because:

* The moment-folded source :math:`Q_\ell` is linear in ``psi_cells``
  (Legendre projection is linear).
* The BC trace value :math:`\bar\phi_{nx+1/2}` is linear in ``psi_cells``
  (BC realisation is linear).
* The inward recurrence is an affine function of (source, BC) with
  constant coefficients (depends on ``Σ_t``, ``Δr`` only).

References
==========

* Hébert, A. (2009).  *Applied Reactor Physics*.  §3.9.4 pp. 141-144,
  Eqs. (3.432)-(3.435).  Local copy:
  ``scratch/literature/Hebert(2009)Chapter3.pdf``.
* Pomraning, G.C. (1989).  *The Transport Equation in General
  Geometry*.  NSE 101:330-340 p. 339.  Cross-reference for the
  structural singularity at r = 0 in curvilinear geometry.
* Literature memo: ``.claude/agent-memory/literature-researcher/phase_d_carlson_coupled_pole.md``.
* Diagnosis memo: ``.claude/agent-memory/numerics-investigator/phase_d_gate_1_1_sphere_mms_diagnosis.md``.
* Diagnostic script: ``tests/sn/diagnostics/gate_1_1_sphere_mms_failure.py``.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import ClassVar, Protocol, runtime_checkable

import numpy as np

from orpheus.numerics.registry import RegistryMixin


# ═══════════════════════════════════════════════════════════════════════
# CarlsonSweepContext — bundle of inputs the Carlson seed strategy needs
# ═══════════════════════════════════════════════════════════════════════


@dataclass(frozen=True)
class CarlsonSweepContext:
    r"""Bundle of external inputs the Carlson inward sweep needs.

    The :class:`CarlsonInwardSweep` strategy operates on the
    Hébert §3.9.4 inward μ = −1 sweep equation, which requires
    information NOT present in the standard
    :class:`~orpheus.sn.spatial.pole_angular_closure.PoleAngularClosure`
    Protocol's call signature (cross-sections, BC at the outer face,
    radial widths).  This context object bundles those inputs so the
    seed strategy receives them via a single kwarg.

    The :class:`ZeroSeed` strategy IGNORES the context entirely
    (Protocol-shape compatibility only); :class:`CarlsonInwardSweep`
    READS every field.

    Parameters
    ----------
    sigma_t : np.ndarray
        Cell-centred total cross-section, shape ``(ng, nx)``.
        Per-group, per-cell.
    dr : np.ndarray
        Radial cell widths, shape ``(nx,)``.  For uniform meshes this
        is constant; for non-uniform meshes the per-cell value enters
        the recurrence's denominator.
    mu_quad : np.ndarray
        Quadrature direction cosines, shape ``(N,)``.  Used to build
        Legendre moments of the input ψ.  GL1D for sphere; for
        cylindrical the level's μ-axis cosines.
    weights : np.ndarray
        Quadrature weights, shape ``(N,)``.  Used in moment integration.
    bc_outer_value : np.ndarray
        Outer-face angular flux value at μ = −1, shape ``(ng,)``.
        For vacuum BC: zeros.  For reflective BC: the outgoing-face
        flux at μ = +1 (mirrored).  For white / albedo: see literature
        memo §8 design decisions; default to 0 with documented
        treatment.

    Notes
    -----
    The context is constructed once per matvec at the call site
    (:func:`~orpheus.sn.operator.transport_operator_matvec_spherical`
    and ``_cylindrical``) and passed through to the M-M closure.
    Carrying the context as a dataclass keeps the call-signature
    expansion local: only one new kwarg on the pole-angular-closure
    Protocol's ``__call__`` instead of four.
    """

    sigma_t: np.ndarray
    dr: np.ndarray
    mu_quad: np.ndarray
    weights: np.ndarray
    bc_outer_value: np.ndarray


# ═══════════════════════════════════════════════════════════════════════
# PsiHalfAngleSeed Protocol
# ═══════════════════════════════════════════════════════════════════════


@runtime_checkable
class PsiHalfAngleSeed(Protocol):
    r"""Strategy contract for the M-M recurrence's half-angle face flux seed.

    Produces the initial ``ψ_{1/2,i,g}`` value for
    :meth:`~orpheus.sn.spatial.pole_angular_closure.MorelMontryAngularSweep._weighted_angular_recurrence_single_level`.

    The recurrence consumes this seed as ``psi_half_left`` at
    ``m = 0`` and propagates it forward across ordinates.

    Parameters
    ----------
    psi_level : np.ndarray
        Cell-centred angular flux at the level's ordinates, shape
        ``(ng, M, nx)``.  Linear in the matvec's INPUT ψ vector;
        used by :class:`CarlsonInwardSweep` to build the moment-folded
        source at μ = −1.
    context : CarlsonSweepContext
        Bundle of external inputs (Σ_t, Δr, μ_quad, weights,
        bc_outer_value).  See :class:`CarlsonSweepContext` for shapes.

    Returns
    -------
    np.ndarray
        Cell-centred half-angle face flux seed, shape ``(ng, nx)``.
        Returned per the M-M recurrence's expectation: one value per
        radial cell, per group.

    Attributes
    ----------
    is_linear : bool
        Whether the seed is linear in ``psi_level``.  Both ZeroSeed
        and CarlsonInwardSweep are linear; required by the operator
        algebra (the M-M closure is linear, and a non-linear seed
        would propagate non-linearity through the recurrence).
    """

    is_linear: bool

    def __call__(
        self,
        psi_level: np.ndarray,
        context: CarlsonSweepContext,
    ) -> np.ndarray:
        ...


# ═══════════════════════════════════════════════════════════════════════
# PsiHalfAngleSeedBase — concrete ABC with self-registration
# ═══════════════════════════════════════════════════════════════════════


class PsiHalfAngleSeedBase(RegistryMixin, ABC):
    r"""Concrete abstract base for self-registering half-angle seed strategies.

    Subclasses inherit this ABC and pass ``key="..."`` in the class
    statement to self-register; the registry is consulted via
    :meth:`PsiHalfAngleSeedBase.create("carlson_inward_sweep")` (or
    any other registered key).

    Mirrors
    :class:`~orpheus.sn.spatial.pole_angular_closure.PoleAngularClosureBase`'s
    Protocol + ABC + RegistryMixin design.

    Subclasses MUST declare:

    * ``is_linear: ClassVar[bool]`` — whether the seed is linear in
      ``psi_level``.
    * :meth:`__call__` — the seed-computation algorithm.
    """

    registry: ClassVar[dict[str, type["PsiHalfAngleSeedBase"]]] = {}

    is_linear: ClassVar[bool]

    @classmethod
    def _registry_base(cls) -> type:
        return PsiHalfAngleSeedBase

    @abstractmethod
    def __call__(
        self,
        psi_level: np.ndarray,
        context: CarlsonSweepContext,
    ) -> np.ndarray:
        ...


# ═══════════════════════════════════════════════════════════════════════
# ZeroSeed — regression-safety ablation
# ═══════════════════════════════════════════════════════════════════════


@dataclass(frozen=True, slots=True)
class ZeroSeed(PsiHalfAngleSeedBase, key="zero"):
    r"""Hardcoded zero seed — reproduces Phase B's pre-fix behaviour.

    Returns ``np.zeros((ng, nx))`` regardless of input.  Exists as
    the regression-safety ablation against
    :class:`CarlsonInwardSweep` (the Phase D canonical fix), and to
    A/B test the Phase D fix's empirical evidence on Gate 1.1.

    **Not a production strategy.**  This reproduces the wrong
    behaviour ERR-026 diagnoses (per
    :doc:`/.claude/skills/numerical-bug-signatures/SKILL` Signature 1
    and Mode 3 failure mode).  Phase D's default-flip activates
    :class:`CarlsonInwardSweep`; :class:`ZeroSeed` is reachable only
    by explicit user opt-in.

    Notes
    -----
    Frozen + slotted: instances are immutable and lightweight.
    """

    is_linear: ClassVar[bool] = True
    """Trivially linear: constant-zero function of input is linear."""

    def __call__(
        self,
        psi_level: np.ndarray,
        context: CarlsonSweepContext,
    ) -> np.ndarray:
        """Return ``np.zeros((ng, nx))`` — the Phase B hardcoded seed."""
        del context  # unused — Protocol-shape compatibility only
        ng, _M, nx = psi_level.shape
        return np.zeros((ng, nx), dtype=psi_level.dtype)

    def __repr__(self) -> str:
        return "ZeroSeed()"


# ═══════════════════════════════════════════════════════════════════════
# CarlsonInwardSweep — Hébert §3.9.4 (3.432)-(3.435) canonical seed
# ═══════════════════════════════════════════════════════════════════════


def carlson_inward_sweep_from_source(
    Q_bar: np.ndarray,
    sigma_t: np.ndarray,
    dr: np.ndarray,
    bc_outer_value: np.ndarray,
) -> np.ndarray:
    r"""Run the Hébert §3.9.4 (3.434)-(3.435) inward μ = −1 DD sweep.

    Source-driven entry point for the Carlson coupled-pole sweep —
    factored from :meth:`CarlsonInwardSweep.__call__` so the SI/sweep
    path (Phase F) can invoke the same recurrence math without going
    through the ψ-input Legendre-moment fold.  At the converged
    eigenmode ``Σ_t·φ_0 = Q̄``, so the matvec path (which builds
    ``Q̄`` from ``Σ_t·Σ_n w_n ψ_n``) and the sweep path (which uses
    ``Q̄`` directly from the within-group source) produce the same
    seed on the fixed point.

    Parameters
    ----------
    Q_bar : np.ndarray, shape ``(ng, nx)``
        Cell-averaged source at :math:`\mu = -1`.  For an L = 0
        isotropic operator: ``Q_bar = (1/2) · Q̄_iso`` where ``Q̄_iso``
        is the angle-averaged within-group source.  The factor ``1/2``
        comes from the Hébert (3.432) Legendre fold at ``ℓ = 0``
        (``(2·0 + 1)/2 = 1/2``) with ``P_0(-1) = 1``.
    sigma_t : np.ndarray, shape ``(ng, nx)``
        Cell-centred total cross-section, per-group, per-cell.
    dr : np.ndarray, shape ``(nx,)``
        Radial cell widths.
    bc_outer_value : np.ndarray, shape ``(ng,)``
        Outer-face angular flux at :math:`\mu = -1`.  Vacuum BC: 0.
        Reflective BC: the outgoing-face flux at :math:`\mu = +1`
        (mirrored), evaluated via the BC operator on the current
        outflow estimate.

    Returns
    -------
    np.ndarray, shape ``(ng, nx)``
        Cell-centred half-angle face flux seed
        :math:`\bar\phi_i \equiv \phi_{1/2,i}`.

    Notes
    -----
    The recurrence is sequential in cells (each step depends on its
    outer neighbour's face value) and vectorised across groups via
    NumPy broadcasting.
    """
    ng, nx = Q_bar.shape
    phi_aux = np.zeros((ng, nx), dtype=Q_bar.dtype)
    phi_face = bc_outer_value.copy()  # (ng,) — outer face at i = nx
    for k in range(nx - 1, -1, -1):
        denom = dr[k] * sigma_t[:, k] + 2.0       # (ng,)
        phi_cell = (dr[k] * Q_bar[:, k] + 2.0 * phi_face) / denom
        phi_aux[:, k] = phi_cell
        # Hébert (3.435): step inward to the next face.
        phi_face = 2.0 * phi_cell - phi_face
    return phi_aux


@dataclass(frozen=True, slots=True)
class CarlsonInwardSweep(
    PsiHalfAngleSeedBase, key="carlson_inward_sweep",
):
    r"""Hébert §3.9.4 Eqs. (3.432)-(3.435) inward μ = −1 sweep.

    Canonical Phase D fix for the M-M angular recurrence's half-angle
    seed.  Solves the streaming–collision balance at μ = −1
    (where the angular redistribution coefficient :math:`1 - \mu^2`
    vanishes and the equation decouples from the α-cascade) and
    returns the cell-centred output as the M-M recurrence seed.

    Algorithm
    ---------

    1. Build Legendre moments of the input ψ at L = 0 (isotropic
       only — the operator's collision term is isotropic in the
       current implementation):

       .. math::

          \phi_0[g, i] \;=\; \sum_n w_n \cdot \psi[g, n, i]

    2. Build the cell-averaged source at μ = −1:

       .. math::

          \bar Q[g, i] \;=\; \frac{1}{2} \cdot \Sigma_t[g, i]
                              \cdot \phi_0[g, i]

       (only the ℓ = 0 term survives for an isotropic operator;
       :math:`P_0(-1) = 1`.)

    3. Run the Hébert (3.434)-(3.435) recurrence INWARD from
       :math:`i = n_x - 1` to :math:`i = 0`:

       .. math::

          \bar\phi_i &= \frac{\Delta r_i \cdot \bar Q_i
                                + 2 \cdot \bar\phi_{i+1/2}}
                              {\Delta r_i \cdot \Sigma_i + 2} \\
          \bar\phi_{i-1/2} &= 2 \cdot \bar\phi_i - \bar\phi_{i+1/2}

    The result :math:`\{\bar\phi_i\}_{i=0..n_x-1}` is returned as
    the M-M recurrence's seed.

    Why this is linear in ``psi_level``
    -----------------------------------

    * :math:`\phi_0` is a linear projection of ``psi_level``.
    * :math:`\bar Q` is :math:`\Sigma_t / 2 \cdot \phi_0` — linear
      in ``psi_level`` (Σ_t is constant input).
    * The Hébert (3.434)-(3.435) recurrence is an affine function of
      :math:`(\bar Q, \bar\phi_{i+1/2})` with constant coefficients
      depending only on (Σ_t, Δr).  Composition of linear maps is
      linear.
    * ``bc_outer_value`` is linear in ``psi_level`` (the matvec
      builds it via cell-centred extraction at the outer face).

    Multi-group: per-group sweep is independent; vectorised across
    groups via broadcasting.  The per-cell inward sweep is genuinely
    sequential (each cell depends on its outer neighbour's face
    value), but the per-group axis is vectorised.

    Notes
    -----
    For multi-region σ_t step discontinuities (different Σ_t per
    cell), the recurrence's denominator ``dr · Σ_t + 2`` is computed
    per cell so the discontinuity is handled trivially (no special
    treatment required as long as the radial mesh respects material
    boundaries).

    The current implementation evaluates only the L = 0 (isotropic)
    Legendre moment.  Anisotropic operators (P1+ scattering) would
    require additional moments :math:`\phi_\ell` and the full
    moment-folded sum :math:`\bar Q = \Sigma_\ell (2\ell+1)/2 \cdot
    \Sigma_t \cdot \phi_\ell \cdot (-1)^\ell`.  Since the
    :class:`~orpheus.sn.operator.SNStreamingOperator` apply matvec
    currently carries an isotropic collision term, the L = 0
    truncation is consistent with the operator's structure.

    .. warning::

       **L = 0 isotropic-only assumption is load-bearing.** This
       strategy assumes the apply matvec's L operator carries ONLY
       the isotropic collision term :math:`\Sigma_t \psi` — scattering
       (including anisotropic P1+ moments) is composed externally
       via the ``SCATTERING`` operator, NOT included in L.  If a
       future refactor moves scattering INTO L (e.g., to enable
       a "monolithic" SN apply that includes within-group scattering),
       this Carlson seed becomes WRONG: the source at :math:`\mu = -1`
       (Hébert Eq. 3.432) needs the full Legendre-moment sum
       :math:`\bar Q_i = \sum_\ell (2\ell+1)/2 \cdot Q_\ell(r) \cdot (-1)^\ell`,
       not just :math:`\Sigma_t \phi_0`.  This is a Mode-6
       convention-drift risk per :ref:`vv-principles
       <vv-principles>`.  A foundation test pinning the
       isotropic-only assumption (e.g., asserting the apply matvec
       does NOT couple to ``self_scattering``) would catch a future
       drift; in its absence, this WARNING block is the only
       safeguard.

    See Also
    --------
    :func:`tests/sn/diagnostics/gate_1_1_sphere_mms_failure.py
    <tests.sn.diagnostics.gate_1_1_sphere_mms_failure>` —
    the diagnostic script that empirically validates this
    intervention (intervention ``[B]``) collapses the Gate 1.1 sphere
    MMS residual to machine precision.
    """

    is_linear: ClassVar[bool] = True
    """Linear in ``psi_level`` — see class docstring §"Why this is linear"."""

    def __call__(
        self,
        psi_level: np.ndarray,
        context: CarlsonSweepContext,
    ) -> np.ndarray:
        r"""Run the Hébert §3.9.4 inward μ = −1 sweep.

        Parameters
        ----------
        psi_level : np.ndarray, shape ``(ng, M, nx)``
            Cell-centred angular flux at the level's ordinates.
        context : CarlsonSweepContext
            Bundle of external inputs (σ_t, Δr, μ_quad, weights,
            bc_outer_value).

        Returns
        -------
        np.ndarray, shape ``(ng, nx)``
            Cell-centred half-angle face flux seed φ̄_{1/2,i}.
        """
        sigma_t = context.sigma_t          # (ng, nx)
        dr = context.dr                    # (nx,)
        weights = context.weights          # (M,)
        bc_outer = context.bc_outer_value  # (ng,)

        # ── Step 1: Legendre ℓ = 0 moment (scalar flux) ─────────────
        # φ_0[g, i] = Σ_n w_n · ψ[g, n, i]
        # einsum keeps the operation vectorised across (g, i); the
        # M axis is contracted.
        phi_0 = np.einsum("m,gmx->gx", weights, psi_level)  # (ng, nx)

        # ── Step 2: cell-averaged source at μ = −1 ──────────────────
        # Q̄[g, i] = (1/Σw_level) · Σ_t[g, i] · φ_0[g, i] · P_0(−1)
        #
        # P_0(−1) = 1; the normalisation factor is the *per-level
        # quadrature weight sum* ``weights.sum() = Σw_level`` so that
        # the isotropic Carlson seed reproduces ``ψ_flat = φ_0 / Σw``
        # under Pomraning structural-singularity isotropy.
        #
        # Phase G Step 2 (Issue #196): the previous form ``0.5 · Σ_t
        # · φ_0`` hardcoded the sphere quadrature convention `Σw = 2`
        # (Hébert §3.9.4's `(2ℓ+1)/2 = 1/2` for ℓ = 0 on the
        # sphere-1D angular measure `dμ` on [−1, 1]).  For 3D
        # quadratures (ProductQuadrature, LevelSymmetric, Lebedev)
        # the per-level weight sum is `2π` and the factor `0.5` is
        # wrong by `2π`, producing a 580% L0 error on the cylinder
        # homogeneous streaming-equilibrium gauntlet (ERR-048
        # manifestation #3).
        #
        # ``1.0 / weights.sum()`` reduces bit-identically to ``0.5``
        # for sphere (single level over `(−1, 1)`, GL-N weights sum
        # to 2) and to ``1/(2π)`` for cylinder ProductQuadrature
        # azimuthal levels — the geometry-general normalisation that
        # restores Pomraning isotropy on flat ψ for any quadrature
        # family.  Per ``coding-elegance`` Pattern 7 (normalise at
        # the definition site, not at every consumer) the literal
        # is replaced with the typed source ``context.weights``.
        Q_bar = sigma_t * phi_0 / weights.sum()  # (ng, nx)

        # ── Step 3: Hébert (3.434)-(3.435) inward DD recurrence ─────
        # Delegated to :func:`carlson_inward_sweep_from_source` so the
        # SI/sweep path (Phase F) can reuse the same recurrence.
        return carlson_inward_sweep_from_source(
            Q_bar=Q_bar,
            sigma_t=sigma_t,
            dr=dr,
            bc_outer_value=bc_outer,
        )

    def __repr__(self) -> str:
        return "CarlsonInwardSweep()"


__all__ = [
    "CarlsonInwardSweep",
    "CarlsonSweepContext",
    "PsiHalfAngleSeed",
    "PsiHalfAngleSeedBase",
    "ZeroSeed",
    "carlson_inward_sweep_from_source",
]
