r"""Shared cross-method solution types.

This module defines the **first cross-method shared vocabulary** for
continuous reference solvers. Two very different mathematical attacks
on the one-speed neutron transport equation —

* the F_N method (Galerkin half-range moment projection of the
  angular flux, realised in
  :mod:`orpheus.derivations.continuous.fn_method`), and
* trajectory_resolvent (the Variant α Green's function carried along
  bouncing characteristics, realised in
  :mod:`orpheus.derivations.continuous.trajectory_resolvent`)

— produce solutions to the same boundary-value problem and therefore
should produce **structurally comparable result containers**. The
math-heart classes that own each method (``MomentSpace`` for F_N,
``Billiard`` for trajectory_resolvent) deliberately return the same
result types so a downstream consumer can hold a ``CriticalSolution``
without knowing which pillar produced it.

The two result types correspond to the two canonical questions one
asks of a transport reference solver:

* :class:`CriticalSolution` — *"What is the critical configuration?"*
  Carries the eigenvalue (``k_eff`` or ``k_inf``) plus the
  configuration parameter the eigenvalue is associated with (the
  critical half-thickness for slab F_N; the critical radius for
  sphere F_N; the eigenvalue at a given ``L`` for trajectory_resolvent
  power iteration). The ``parameter_kind`` field disambiguates which.
* :class:`FluxSolution` — *"Given a configuration, what is the flux
  shape?"* Carries the scalar flux (and optionally the angular flux),
  spatial / angular grids, and the eigenvalue at which the flux
  shape was extracted.

Why this lives in ``derivations/common/``
=========================================

Per the algebra-of-record skill (``algebra-of-record/SKILL.md``), the
``common`` package is the **shared layer above the trusted-library
line** — it carries primitives that both Branch-1 (SymPy /
SymPy+mpmath) and Branch-2 (numpy / scipy) implementations of all
methods share. Solution types belong here because they describe
*what is being computed* in a method-agnostic way, without bringing
any method's machinery into scope. (The retired ``GeometrySpec``
carrier lived here for the same reason on the *input* side; Phase F.3
deleted it in favour of the geometry-layer
:class:`~orpheus.geometry.structured_geometry.StructuredGeometry`,
which every reference solver now consumes directly.)

Why this is **not** a Protocol
==============================

Per the project's "unify after two instances" memory
(``feedback_unify_after_two_instances.md``), the unifying Protocol
across math-heart classes can be designed only AFTER ≥2 working
instances exist. ``MomentSpace`` (F_N) and ``Billiard``
(trajectory_resolvent) ARE the first two instances. Their solution
types — defined here — are the first concrete unification: a
*structural* contract on what every math-heart returns, before any
behavioural Protocol is posited.

The deliberate choice is to UNIFY the result types eagerly (they're
small, stable dataclasses; the cost of premature unification is
trivial) but DEFER the Protocol over the math-heart classes
themselves until both are working and the patterns of variation are
empirically observed.

Field philosophy
================

These dataclasses are **frozen** and **slot-friendly** — no behaviour,
no methods beyond ``__repr__`` and equality. They are containers,
not actors. Method-specific diagnostics (number of F_N modes, number
of trajectory quadrature nodes, determinant residuals, Atkinson
panel counts) live on the math-heart class instance and can be
queried via the ``metadata`` dict if reproducibility from the
result alone is required.

The ``metadata`` field is open-ended on purpose: every solver carries
its own diagnostic vocabulary, and forcing a closed schema would
trade structural simplicity for premature standardisation. The
canonical names expected by the cross-method protocol
(:mod:`tests.cross_method.adapters`) are documented per-adapter; new
solvers SHOULD follow precedent but aren't structurally constrained
to.

References
----------

* :doc:`/skills/algebra-of-record` — the bifurcation discipline that
  this module's "Branch-1 + Branch-2 share the same return type"
  pattern slots into. Both branches return the same shape because
  they should be substitutable at the call site (``MomentSpace`` can
  be backed by either; ``Billiard`` can be backed by either).
* :doc:`/skills/vv-principles` § "The three pillars of verification"
  — these dataclasses are the *interface* across pillars. A
  closed-form analytical reference, a semi-analytical solver, and an
  MMS-derived reference all populate the same ``CriticalSolution`` /
  ``FluxSolution``, which is what makes them substitutable in
  cross-check tests.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal, Optional

import numpy as np


# Type alias for the parameter-kind tag on CriticalSolution. The
# enumeration is closed (these are the only configuration parameters
# any 1-D reference solver in this codebase reports), but stays a
# Literal rather than an Enum to keep the dataclasses pickleable
# without extra import boilerplate downstream.
ParameterKind = Literal[
    "half_thickness_mfp",   # bare slab: half-thickness in mean free paths
    "radius_mfp",           # sphere/cylinder: radius in mean free paths
    "half_thickness_cm",    # bare slab: half-thickness in cm
    "radius_cm",            # sphere/cylinder: radius in cm
    "domain_extent_cm",     # generic: full domain extent in cm (slab L, sphere R)
    "k_inf_only",           # infinite medium: no spatial parameter
    "fixed_geometry",       # power iteration at a given fixed config
]


@dataclass(frozen=True)
class CriticalSolution:
    r"""Solution of a critical-configuration problem.

    A reference solver in *critical-configuration mode* reports both:

    * an **eigenvalue** — typically ``k_eff = 1`` for a critical
      bare-multiplying configuration (the critical-radius / critical-
      thickness root-find), or ``k_inf`` for an infinite medium, or
      ``k_eff`` evaluated at a given fixed geometry (Variant α
      power iteration), and
    * a **configuration parameter** — the half-thickness, radius,
      domain extent, or other geometric scalar that the eigenvalue
      is associated with.

    The two together completely identify the (shape, configuration)
    pair the solver converged to. Detail beyond this — internal
    expansion coefficients, iteration count, residual norms — lives
    in :attr:`metadata`.

    Parameters
    ----------
    eigenvalue : float
        The reported eigenvalue. Naming convention is method-specific
        and recorded in :attr:`eigenvalue_kind` (``"k_eff"``,
        ``"k_inf"``, ``"sigma_a"`` for fixed-source absorption-rate
        eigenvalue, etc.).
    eigenvalue_kind : str
        One of ``"k_eff"``, ``"k_inf"``, or any other documented tag.
        Free-form to allow novel conventions; comparators across
        methods should check this field before comparing eigenvalues.
    parameter_value : float
        The configuration parameter (typically the critical
        half-thickness / radius). Units encoded in
        :attr:`parameter_kind`.
    parameter_kind : ParameterKind
        Tag for what :attr:`parameter_value` represents and what
        units it carries. See the :data:`ParameterKind` literal.
    converged : bool
        Whether the iteration / root-find converged to its tolerance.
    metadata : dict[str, Any]
        Method-specific diagnostics. Free-form. Conventions used in
        the cross-method test protocol:

        * ``"n_modes"`` (F_N method) — F_N order :math:`N`.
        * ``"determinant_residual"`` (F_N) — :math:`\det M` at the
          converged configuration.
        * ``"iterations"`` (trajectory_resolvent) — power-iteration
          step count.
        * ``"n_traj_quad"`` (trajectory_resolvent) — trajectory
          quadrature density.

    Notes
    -----
    The choice to expose both ``eigenvalue`` and ``parameter_value``
    as plain floats (rather than separate sub-dataclasses or a
    method-specific Result type) is deliberate: the cross-method
    protocol consumes both as scalars. Users who want the rich
    method-specific result type (e.g., ``SlabFNResult`` with its
    F_N expansion coefficients) should call the underlying solver
    directly; ``MomentSpace.solve_critical()`` and
    ``Billiard.solve_critical()`` populate this dataclass for
    cross-method comparability and stash the rich result in
    ``metadata["raw_result"]`` for callers that want both.
    """
    eigenvalue: float
    eigenvalue_kind: str
    parameter_value: float
    parameter_kind: ParameterKind
    converged: bool
    """Whether the producing solve met its tolerance — **required, no default**.

    ⛔ This defaulted to ``True`` until 2026-08-09.  That is the same defect
    as #342 with the assertion moved into the type: a producer that never
    thought about convergence claimed it anyway.  A field whose SAFE value
    is the optimistic one lies by omission.

    `[M]` removing the default cost **zero** churn — all 33 construction
    sites already passed it explicitly, so the default protected nothing
    and only stood ready to hide the next producer that forgot.  Twin of
    :attr:`orpheus.sn.solution.IterationHistory.converged`.
    """
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class FluxSolution:
    r"""Solution of a flux-reconstruction problem.

    A reference solver in *flux-reconstruction mode* reports the flux
    shape on a known configuration. Carries:

    * **scalar flux** :math:`\phi(r)` on a 1-D spatial grid (always),
    * **angular flux** :math:`\psi(r, \mu)` on a (spatial, polar)
      grid (optional — ``None`` if the method only computes the
      scalar flux),
    * the **eigenvalue** at which the flux was extracted (since flux
      shape and eigenvalue are tied in a critical-configuration
      problem).

    Parameters
    ----------
    spatial_nodes : np.ndarray
        ``(n_x,)`` array of spatial coordinates. Units are method-
        specific (mfp or cm) and recorded in :attr:`spatial_units`.
    scalar_flux : np.ndarray
        ``(n_x,)`` scalar flux :math:`\phi`. Normalised so that the
        peak value is 1, unless :attr:`metadata` says otherwise
        (the cross-method gates always work with ratios; absolute
        normalisation is a method-specific convention).
    angular_flux : np.ndarray | None
        ``(n_x, n_mu)`` angular flux :math:`\psi(r, \mu)` on a
        product grid, or ``None`` if the method does not return
        angular flux. Same normalisation convention as
        :attr:`scalar_flux`.
    angular_nodes : np.ndarray | None
        ``(n_mu,)`` polar quadrature nodes :math:`\mu_k` if
        :attr:`angular_flux` is set; ``None`` otherwise.
    eigenvalue : float
        The eigenvalue at which the flux shape was extracted. For
        bare-critical configurations this is ``k_eff = 1`` to the
        solver's tolerance; for fixed-source / sub-critical
        configurations this can be any positive real.
    eigenvalue_kind : str
        ``"k_eff"`` / ``"k_inf"`` / ``"none"`` — the same convention
        as :class:`CriticalSolution`. ``"none"`` for fixed-source
        solutions where no eigenvalue is computed.
    spatial_units : Literal["mfp", "cm"]
        Units of :attr:`spatial_nodes`.
    metadata : dict[str, Any]
        Method-specific diagnostics — same convention as
        :class:`CriticalSolution`.

    Notes
    -----
    The grid layout (cell-centred vs cell-face, uniform vs
    Chebyshev-distributed, half-range vs full-range :math:`\mu`) is
    method-specific and not normalised here — users should always
    interpolate to a target grid before comparing flux shapes across
    methods. The cross-method gates compare flux ratios at a small
    set of well-defined points (e.g., :math:`\phi(0)/\phi(R)`) rather
    than full-grid agreement, sidestepping the normalisation question.
    """
    spatial_nodes: np.ndarray
    scalar_flux: np.ndarray
    angular_flux: Optional[np.ndarray]
    angular_nodes: Optional[np.ndarray]
    eigenvalue: float
    eigenvalue_kind: str
    spatial_units: Literal["mfp", "cm"]
    metadata: dict[str, Any] = field(default_factory=dict)


__all__ = [
    "ParameterKind",
    "CriticalSolution",
    "FluxSolution",
]
