r"""The GAUGES — the section that picked the representative a solve returns.

A converged solve does not return *the* answer; it returns one member of the
answer's SOLUTION SET, and the set carries a symmetry group Γ that acts freely
and transitively on it (a Γ-torsor):

* the eigen kind — :math:`A\psi = \mu M\psi` — has a RAY in the cone plus a
  scalar; the ray is a torsor under :math:`(\mathbb{R}_+, \times)` (and, on the
  dense engines, its sign is a :math:`\mathbb{Z}/2` leg on top);
* the source kind — :math:`A\psi = q` — has a COSET :math:`\psi_0 + \ker A`, a
  torsor under :math:`(\ker A, +)`, non-trivial on an all-reflective Cartesian
  box closed by diamond differencing (the loss operator is exactly singular
  there, #344).

Returning a representative means choosing a SECTION of the quotient
:math:`S \to S/\Gamma`; that choice is data — it says which member you hold
and how far it was moved — and until step 3 of the consumers campaign
(2026-09-17) no result type in the tree recorded it.  ``[M]`` the eigen
normalization shipped as FOUR functionals under the one name "production
rate" (SN's ``compute_production_rate`` INCLUDING the (n,2n) emission,
``KEigenvalue``'s fission-only member, diffusion's :math:`\langle\nu\Sigma_f,
\phi\rangle`, the homogeneous solver's rescale to 100), applied
conditionally (``isinstance(solver, ProductionRateSolver)`` inside
``power_iteration``), and recorded nowhere: a consumer holding a Solution
could not determine the scale of the flux it held.

Two realizations of ONE concept:

* :class:`ScaleGauge` — a degree-1 homogeneous functional ``n`` and a target
  ``t``: the section :math:`\psi \mapsto \psi \cdot t / n(\psi)`.  The
  ``functional`` is stored as the OBJECT that ran (a bound method, a
  co-vector's ``evaluate``), never as a label — a label would record a
  falsehood the day two solvers spell "production rate" differently, which
  ``[M]`` they do.
* :class:`KernelGauge` — a PROTOCOL: the :math:`G`-orthogonal projector
  :math:`\Pi` onto :math:`\ker A` with the section :math:`\psi \mapsto \psi -
  \Pi\psi` (the minimum-:math:`G`-norm member, where the exact solution sits
  — a theorem, not a convention: every kernel mode is mirror-odd).  Satisfied
  structurally by the SN tier's ``LossKernelGauge`` — the kernel is a TRACE
  object there, so ``gauge`` acts on the boundary-trace values and the bulk is
  untouched by construction; this tier never imports the SN class (the layer
  contract forbids ``numerics → sn``).

The defining laws each realization ships a test of (``tests/gates/numerics/``):
the section lands on the target (``n(apply(ψ)) ≈ t``), it is IDEMPOTENT
(``apply(apply(ψ)) ≈ apply(ψ)`` — ``allclose``, not bit-equal: a rescale
multiplies by :math:`1/(1 \pm \varepsilon)`), it is Γ-INVARIANT
(``apply(γ·ψ) ≈ apply(ψ)`` for every group element — for the kernel gauge
the element is drawn from :math:`\Pi`'s own range on the TRACE space), and,
for the kernel gauge, it is RESIDUAL-NEUTRAL (:math:`A(\psi - \Pi\psi) =
A\psi` because :math:`\Pi\psi \in \ker A`).  The negative leg that motivated
the type: two solves of ONE problem under two functionals return fluxes that
differ by a positive scalar, and their recorded gauges must be
DISTINGUISHABLE values.

See ``.claude/plans/consumers_step3_design.md`` §3.1–§3.2 and the frame pass
``scratch/_consumers/step3/attacker_solution_articulation.md`` §4 (the
torsor / projective frames).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Protocol, runtime_checkable

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True)
class ScaleGauge:
    r"""The section of the :math:`\mathbb{R}_+`-quotient a degree-1 functional fixes.

    Parameters
    ----------
    functional : callable
        A degree-1 homogeneous functional :math:`n` of the state —
        :math:`n(c\psi) = c\,n(\psi)` — stored as the OBJECT that produced the
        returned scale (``SNSolver.compute_production_rate``,
        ``KEigenvalue.compute_production_rate``, a reaction-rate co-vector's
        ``evaluate``).  Never a label.
    target : float
        The value :math:`t` the section lands on: ``1.0`` for the
        ``power_iteration`` family, ``100.0`` for the homogeneous infinite
        medium.

    Notes
    -----
    A negative functional value is a legal input: the displacement is then
    negative and :meth:`apply` flips the sign — the :math:`\mathbb{Z}/2` leg
    the dense engines fix by hand (``_sign_normalised``) is one case of this
    section, not a second convention.  A ZERO functional value has no section
    (the state is orthogonal to the functional — nothing to scale against) and
    is refused with a ``ValueError`` naming the functional.  ``[M]``
    ``power_iteration``'s own rescale is guarded by ``if p > 0.0``
    (``eigenvalue.py:437-447``), so it silently skips BOTH the zero reading
    (where no section exists) and a negative one (where this section exists
    and flips the sign) — two cases the guard cannot tell apart.
    """

    functional: Callable[[Any], float]
    target: float

    def displacement(self, state: Any) -> float:
        r"""The group element :math:`t / n(\psi)` that carries ``state`` onto the section."""
        n = float(self.functional(state))
        if n == 0.0:
            raise ValueError(
                f"ScaleGauge: the functional {self.functional!r} reads 0 on the "
                f"state — no section of the scale quotient exists for it "
                f"(the state is orthogonal to the gauge functional)."
            )
        return self.target / n

    def apply(self, state: Any) -> Any:
        r"""The canonical representative :math:`\psi \cdot t / n(\psi)` — idempotent to ``allclose``."""
        return state * self.displacement(state)


@runtime_checkable
class KernelGauge(Protocol):
    r"""The section of the :math:`(\ker A, +)`-quotient: :math:`\psi \mapsto \psi - \Pi\psi`.

    A structural contract — the SN tier's ``LossKernelGauge`` satisfies it
    without inheriting (the layer contract forbids ``numerics → sn``, and a
    projector is what the concept IS, whichever tier builds it).

    ``gauge(trace)`` returns the canonical member of the trace's coset (the
    minimum-:math:`G`-norm one); ``dimension`` is :math:`\dim \ker A` on the
    component the projector spans — ``0`` when the configuration has no
    freedom, in which case ``gauge`` is the identity and NO caller needs a
    ``None`` branch (the zero-block projector is the honest value, not an
    absence).
    """

    def gauge(self, trace: NDArray[np.floating]) -> NDArray[np.floating]: ...

    @property
    def dimension(self) -> int: ...
