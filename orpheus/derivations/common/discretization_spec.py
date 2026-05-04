r""":class:`DiscretizationSpec` — discrete-mesh-adapter parameter dataclass.

This module defines a small, immutable specification of *discretization*
parameters for the production CP / SN / MOC solvers when wrapped behind
the :class:`~orpheus.derivations.common.solver_protocol.TransportSolver`
Protocol. It is the discrete-mesh counterpart to the method-specific
kwargs that ``MomentSpace`` (``fn_order``) and ``Billiard``
(``quadrature``) expose.

Where it sits
-------------

The three-axis carve from the Protocol theory page is:

1. **What** — :class:`~orpheus.derivations.common.geometry_spec.GeometrySpec`
   + ``materials`` (method-agnostic).
2. **How** — method-specific parameters: ``fn_order`` for F_N,
   :class:`DiscretizationSpec` for production CP/SN, etc.
3. **What is asked** —
   :class:`~orpheus.derivations.common.solution_types.CriticalSolution`
   / :class:`~orpheus.derivations.common.solution_types.FluxSolution`
   (method-agnostic).

This dataclass lives at the (2) layer for production discrete solvers.
A continuous-reference math-heart class does NOT consume it — its
"how" is its own dataclass field (``fn_order`` / ``quadrature``).
The shared layer is the (1) and (3) shape, which is enforced by the
Protocol.

What it does NOT cover
----------------------

* Solver tolerances / iteration counts (``max_outer``, ``keff_tol``,
  ``flux_tol``) — those are CP/SN solver concerns, not
  *discretization* concerns. Pass them through the production
  solver's ``CPParams`` / SN kwargs instead.
* Boundary conditions — those live on
  :class:`GeometrySpec.bc_left` / :attr:`GeometrySpec.bc_right`. The
  discretization-spec layer is BC-agnostic.
* Quadrature *rule* selection (Gauss-Legendre vs Lebedev vs
  level-symmetric for SN). That is owned by the SN adapter's choice
  of :class:`AngularQuadrature` factory at build time. We expose a
  scalar ``n_angular`` here that downstream maps to the natural
  family for the geometry — slab GL of order ``n_angular``,
  spherical Lebedev of order :math:`\ge` ``n_angular``.

Pure dataclass, no behaviour
----------------------------

The class is frozen and carries no methods beyond the dataclass-
generated ``__repr__`` / ``__eq__`` / ``__hash__``. Construction
validation is minimal: positive-int checks on each field. The
adapter consumes the spec at ``from_problem`` time and translates
each scalar into the production solver's native parameter shape.

References
----------

* :doc:`/theory/transport_solver_protocol` § "Discretisation
  parameters for discrete adapters" — the rich-narrative theory page
  that frames this dataclass within the three-axis carve.
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class DiscretizationSpec:
    r"""Discretization parameters for discrete-mesh transport adapters.

    The fields enumerate the smallest set of scalars that pin a
    production CP / SN solve from the discrete-adapter scaffold:

    * :attr:`n_cells` — the number of equal-volume sub-cells in the
      1-D production mesh built via
      :meth:`GeometrySpec.build(n_cells=...)`.
    * :attr:`n_angular` — the angular quadrature order. Maps to:

      - Slab SN: Gauss-Legendre of order ``n_angular``.
      - Spherical / cylindrical SN: Lebedev rule of order
        :math:`\ge` ``n_angular`` (rounded up to a Lebedev-supported
        order).
      - CP: ignored at the discretization-spec layer (CP integrates
        analytically over angle via ``E_3`` / ``Ki_3``); the SN
        adapter is the one that consumes ``n_angular``.

    * :attr:`n_chord_quad` — the chord-integration quadrature order
      for CP curvilinear geometries (``CPParams.n_quad_y``). Pinned
      to a default that matches :class:`CPParams.n_quad_y = 64`. For
      slab CP the chord integral is closed-form; the field is
      ignored for slabs.

    Defaults match the production-solver-side defaults (``n_cells =
    20``, ``n_chord_quad = 32``) so a bare ``DiscretizationSpec()``
    constructs to a usable smoke-test configuration.

    Examples
    --------

    Construct a discretisation spec for a small smoke test::

        from orpheus.derivations.common.discretization_spec import (
            DiscretizationSpec,
        )

        spec = DiscretizationSpec(
            n_cells=20,
            n_angular=8,
            n_chord_quad=32,
        )

    Pass to a production discrete solver via its
    ``from_problem`` Protocol entry point::

        from orpheus.cp.solver import CPSolver

        cp = CPSolver.from_problem(
            materials={0: pu_mixture},
            geometry_spec=GeometrySpec(geometry="sphere", ...),
            discretization=spec,
        )
        sol = cp.solve_critical()

    Why a dedicated dataclass and not a dict
    -----------------------------------------

    A dict would have hidden the shape contract and made downstream
    typing fuzzy. The dataclass:

    * Documents the field set in one place (this docstring + IDE
      hover).
    * Catches typos at construction time (``n_celles=20`` raises
      :class:`TypeError`).
    * Hashable + frozen, so adapters can cache solves keyed by the
      spec without surprising mutation issues.

    Attributes
    ----------
    n_cells : int, default 20
        Number of equal-volume sub-cells passed to
        :meth:`GeometrySpec.build`.
    n_angular : int, default 16
        Angular quadrature order. SN-only.
    n_chord_quad : int, default 32
        CP curvilinear chord integration order. CP-only.
    """

    n_cells: int = 20
    n_angular: int = 16
    n_chord_quad: int = 32

    def __post_init__(self) -> None:
        """Validate that all fields are positive integers."""
        for field_name, value in (
            ("n_cells", self.n_cells),
            ("n_angular", self.n_angular),
            ("n_chord_quad", self.n_chord_quad),
        ):
            if not isinstance(value, int):
                raise TypeError(
                    f"DiscretizationSpec.{field_name} must be int, "
                    f"got {type(value).__name__}={value!r}"
                )
            if value <= 0:
                raise ValueError(
                    f"DiscretizationSpec.{field_name} must be positive, "
                    f"got {value}"
                )


__all__ = ["DiscretizationSpec"]
