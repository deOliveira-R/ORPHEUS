r"""The OUTCOMES — what a solve ANSWERED, fused with the question it answered
and the gauge that picked the representative; and the exit CERTIFICATE.

Step 3 of the consumers campaign (``.claude/plans/consumers_step3_design.md``,
RULED 2026-09-14): a Solution is the pair (Problem, posing) plus the Strategy
that produced it and the records, and **its KIND is the posing's TYPE** —
never ``keff is not None``, the answer-derived discrimination the SN
``Solution`` carried until this step (``[M]`` ``sn/solution.py:543-549``, a
kind read one tier too late: the multiplying-source Solution was
indistinguishable from a pure-transport one by its data, and the eigen gauge
was recorded nowhere).

Shape (B) at the Solution tier — TWO kind-typed outcomes, ZERO Optionals, ONE
fused member so that every illegal pairing is UNCONSTRUCTIBLE rather than
merely refused (``coding-elegance`` Pattern 4):

* :class:`EigenOutcome` ``(posing: EigenPosing, state, lam, trajectory,
  gauge: ScaleGauge)`` — the eigen kind: a RAY representative plus the
  physical eigenvalue λ (the spectral map's reading of the pencil's μ) plus
  the scale section that fixed the representative.  Its adjoint question is
  NULLARY (``posing.H()``, :math:`k^\dagger = k`).
* :class:`SourceOutcome` ``(posing: SourcePosing, state, gauge: KernelGauge)``
  — the source kind: a COSET representative plus the kernel section.  It has
  no λ, no ``keff``, no nullary adjoint (a source problem's dual needs a
  DETECTOR — ``posing.H(detector)`` stays on the posing).

The pairings ``(EigenPosing, KernelGauge)``, ``(SourcePosing, λ)`` and "a
state beside a question of the other kind" cannot be spelled: the constructor
of each outcome PARSES its posing's type once, at the boundary, and the
kind-specific verbs (``keff``, ``rayleigh``, ``adjoint_posing``) exist only on
the kind that has them — a ``SourceOutcome.keff`` is a type error, not
``None``.

**The certificate.** What was MEASURED about the returned state, and — when
nothing was — WHY NOT, as a typed sum instead of a ``None`` with five
documented meanings (``[M]`` ``IterationHistory.balance_defect`` until this
step): :class:`Measured` (a number), :class:`Certified` (a bound the exit
certificate asserted, so no number was needed), :class:`NotApplicable` (the
question does not arise — a zero source, no kernel freedom, a pure-transport
posing has no admissibility to check), :class:`NotYet` (open work, by issue —
the carrying eigen exit's balance #354, the daggered eigen exit's #353, the
LD residual mint #310).  A :class:`ExitCertificate` carries four members —
``balance`` (the per-group balance defect ratio, a DIAGNOSTIC never a gate:
#340 N5), ``gauge`` (the kernel-gauge displacement :math:`\lVert\Pi\psi\rVert
/ \lVert\psi\rVert`), ``rayleigh_gap`` (:math:`|\lambda - \lambda_{\rm
Rayleigh}(\psi)|` — the reference-class agreement between the method-tier
estimator and the posing's own quotient, RECORDED, never asserted at
construction: a bare ``assert`` is inert under ``-O`` and a ``raise`` on a
magnitude is the refuted N5), and ``admissibility`` (the multiplying source
problem's certificate: the hub's own :math:`k_{\rm eff} < 1`, with the
tolerance it was measured at).  There is deliberately no ``value_or_none``
accessor — it would re-import the leak the sum retires.

The chain, in one line: ``Problem.pencil → posing → (Strategy) → outcome +
certificate + record → Solution``.  The outcome types live at the numerics
tier because they know nothing of SN: a ``HomogeneousResult`` carries an
:class:`EigenOutcome` too (the 0-D pencil under the k map, the ``=100``
production-rate gauge — the one shipped deliberately-named target).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Generic

from orpheus.numerics.gauge import KernelGauge, ScaleGauge
from orpheus.numerics.posing import EigenPosing, SourcePosing
from orpheus.numerics.vector import Carrier

# The outcomes are generic in their CARRIER through the tree's deliberately
# UNBOUNDED ``Carrier`` (``numerics/vector.py``, the ``power_iteration``
# precedent): numpy's stubs cannot prove ``ndarray ⊨ Vector``, and the 0-D
# family's state IS an ndarray column.  The posing's own carrier parameter is
# erased here for the same reason; the posing checks its ends at runtime.


# ── the evidence sum ──────────────────────────────────────────────────────


@dataclass(frozen=True)
class Measured:
    """A number was measured on the returned state."""

    value: float


@dataclass(frozen=True)
class Certified:
    """No number was needed: the exit certificate ASSERTED the bound (raising on a defect beyond it)."""

    bound: float
    by: str


@dataclass(frozen=True)
class NotApplicable:
    """The question does not arise for this posing / configuration — say why."""

    reason: str


@dataclass(frozen=True)
class NotYet:
    """Open work — the issue that will make it measurable, and what is missing."""

    issue: int
    reason: str


Evidence = Measured | Certified | NotApplicable | NotYet


@dataclass(frozen=True)
class ExitCertificate:
    r"""What the exit measured about the RETURNED state, member by member (see the module docstring)."""

    balance: Evidence
    gauge: Evidence
    rayleigh_gap: Evidence
    admissibility: Evidence


# ── the two kinds ─────────────────────────────────────────────────────────


@dataclass(frozen=True)
class EigenOutcome(Generic[Carrier]):
    r"""The eigen kind's answer: a ray representative, λ, the λ-trajectory, the scale section.

    Parameters
    ----------
    posing : EigenPosing
        The question — the pencil :math:`(A, M)` with its spectral map.
    state : Carrier
        The returned representative of the ray — the driver's iterate AFTER the
        gauge (production rate on the section's target).
    lam : float
        The PHYSICAL eigenvalue the driver produced (k under ``K_MAP``, α under
        ``ALPHA_MAP``) — the method-tier estimator's value, bit-for-bit what
        the driver converged.  ``[M]`` the posing's own Rayleigh quotient
        agrees with it to the convergence residual (a REFERENCE-class pair);
        the gap is the certificate's ``rayleigh_gap``, not this field's law.
    trajectory : tuple of float
        λ at every outer step, the last one ``lam`` (co-indexing enforced at
        construction; ``(lam,)`` for a direct solve).
    gauge : ScaleGauge
        The section that fixed the representative — the functional AS THE
        OBJECT THAT RAN and its target.
    """

    posing: EigenPosing[Any]
    state: Carrier
    lam: float
    trajectory: tuple[float, ...]
    gauge: ScaleGauge

    def __post_init__(self) -> None:
        # The ONE kind parse, at the boundary (parse, don't validate): a source
        # question cannot be paired with a λ and a scale section.  pyright sees
        # the annotation; this is the untyped path's refusal.
        if not isinstance(self.posing, EigenPosing):
            raise TypeError(
                f"EigenOutcome: the posing must be an EigenPosing (the eigen "
                f"question); got {type(self.posing).__name__} — a source "
                f"question has no eigenvalue and no scale section; build a "
                f"SourceOutcome."
            )
        if self.trajectory and float(self.trajectory[-1]) != float(self.lam):
            raise ValueError(
                f"EigenOutcome: the trajectory's last entry ({self.trajectory[-1]!r}) "
                f"must be lam ({self.lam!r}) — the trajectory is λ at every outer "
                f"step and lam is the last of them."
            )

    # ── the verbs, all delegating to the question ────────────────────────

    def residual(self) -> Any:
        r"""The eigen-residual :math:`A\psi - \mu(\lambda) M\psi` of the returned state."""
        return self.posing.residual(self.state, self.lam)

    def balance(self, w: Any = 1.0) -> float:
        r""":math:`\beta(\psi, \lambda) = \langle w, \mathcal{A}(\mu(\lambda))\psi\rangle` — zero at the solution."""
        return self.posing.balance(self.state, self.lam, w)

    def rayleigh(self, w: Any = 1.0) -> float:
        r"""The posing's own eigenvalue estimate on the returned state (``w = 1``: the balance member; ``w = ψ†``: the stationary one)."""
        return self.posing.rayleigh(self.state, w)

    @property
    def keff(self) -> float:
        r"""λ read as the multiplication factor — legal only under the k map."""
        if self.posing.spectral_map.name != "k":
            raise ValueError(
                f"EigenOutcome.keff: the posing's spectral map is "
                f"{self.posing.spectral_map.name!r}, not 'k' — read `lam` "
                f"(the physical eigenvalue under that map) instead."
            )
        return self.lam

    def dominance_ratio(self) -> float | None:
        r"""The empirical dominance ratio :math:`|\lambda_n - \lambda_{n-1}| / |\lambda_{n-1}|` — ``None`` for a single-point trajectory (undefined, not unmeasured)."""
        if len(self.trajectory) < 2:
            return None
        prev = self.trajectory[-2]
        if prev == 0.0:
            return None
        return abs(self.trajectory[-1] - prev) / abs(prev)

    @property
    def adjoint_posing(self) -> EigenPosing[Any]:
        r"""NULLARY — the daggered question :math:`(A^\dagger, M^\dagger)` shares the spectrum; no datum is needed."""
        return self.posing.H()


@dataclass(frozen=True)
class SourceOutcome(Generic[Carrier]):
    r"""The source kind's answer: a coset representative and the kernel section that picked it.

    Parameters
    ----------
    posing : SourcePosing
        The affine question :math:`A\psi = q` — over the pure transport operator
        (``pencil.at(0)``, the driver's own point in Λ on a fissile hub), or
        over the physical multiplying member (``hub.source_posing(q)``,
        ``pencil.at(1)``).
    state : Carrier
        The returned representative of the coset — the driver's iterate AFTER
        the kernel gauge.
    gauge : KernelGauge
        The projector :math:`\Pi` onto :math:`\ker A` whose section returned
        the minimum-norm member — the zero-block projector when the
        configuration has no freedom (a value, not an absence).

    No ``keff``, no ``lam``, no nullary adjoint: the adjoint source problem needs
    a DETECTOR (``posing.H(detector)``).
    """

    posing: SourcePosing[Any]
    state: Carrier
    gauge: KernelGauge

    def __post_init__(self) -> None:
        if not isinstance(self.posing, SourcePosing):
            raise TypeError(
                f"SourceOutcome: the posing must be a SourcePosing (the affine "
                f"question); got {type(self.posing).__name__} — an eigen "
                f"question has no source and its answer carries λ; build an "
                f"EigenOutcome."
            )

    def residual(self) -> Any:
        r""":math:`A\psi - q` on the returned state (the loss-sign convention every SN residual uses)."""
        return self.posing.residual(self.state)

    def balance(self, w: Any = 1.0) -> float:
        r""":math:`\beta(\psi) = \langle w, A\psi - q\rangle` — the exit balance defect, signed, zero at the solution."""
        return self.posing.balance(self.state, w)
