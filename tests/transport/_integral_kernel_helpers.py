r"""Shared fixtures + not-yet-written-SUT import guards for the #257 S6
``IntegralKernelOperator`` category verification suite.

S6 is PRE-IMPLEMENTATION when these specs are written. The
system-under-test does NOT yet exist on the tree:

* the Protocol ``orpheus.transport.operators.integral_kernel_operator.IntegralKernelOperator``;
* the new ``FissionOperator.production_rate`` property (the S5 functional);
* the new ``ScatteringOperator.kernel`` property (the typed ``R∘Λ∘M``).

Every S6 spec routes its NOT-YET-LANDED SUT access through a probe guard
so the files **collect** (green session, no ``ImportError``) and **skip**
with a clear reason until the method-implementer lands the code. Once
landed, the skips become real assertions automatically — no edit to the
spec files is required. This mirrors the S5 ``_functional_helpers.py``
pattern verbatim.

What S6 adds (the system-under-test — additive, bit-identical)
--------------------------------------------------------------

S6 introduces the §5.6 *Kernel* suffix of the operator-algebra category
taxonomy (Operator / Kernel / Functional / Projection / Reconstruction).
An :class:`IntegralKernelOperator` is a :class:`LinearOperator` whose
action is **nonlocal** — it integrates a field against a measure on ≥1
axis — distinct from a :class:`MultiplicationOperator` (local / diagonal,
S3b) and from a :class:`Functional` (field → scalar, S5). The sole
required member is a ``kernel`` property returning a :class:`LinearOperator`.

* ``FissionOperator`` gains a ``production_rate`` property → the S5
  :class:`ProductionRateFunctional` over ``νΣf`` — and so satisfies the
  Protocol (it already exposes ``kernel`` since Wave T). The matvec arms
  are UNCHANGED.
* ``ScatteringOperator`` gains a ``kernel`` property → the typed
  ``OperatorProduct(R, OperatorProduct(Λ, M))`` reproducing the existing
  ``_aniso_source_from_moment_values(M·ψ)`` chain. The 5 dispatch arms
  are UNCHANGED.

vv claim layer (1.5 gate)
-------------------------

Every row across the S6 suite is one of:

* a CATEGORY / structural claim (Protocol membership) — the reference is
  the Protocol definitions themselves (positive + negative + discriminator
  make the gate non-self-referential, L11);
* a BIT-IDENTICAL equivalence/de-risk claim (the new ``kernel`` /
  ``production_rate`` property reproduces the EXISTING fused realization)
  — 0 ULP, ``assert_array_equal``; clearly demarcated as equivalence, NOT
  correctness (the structurally-independent reference for the scattering
  *physics* is the existing aniso MMS gate, NOT this equivalence).

Zero eigenvalue claims, zero MMS, zero convergence-order — so no
structurally-independent NUMERICAL reference is asserted here.

vv Mode-8: the canonical ORPHEUS invocation is ``python -O``, which
strips bare ``assert`` to a NO-OP. Every structural gate routes through
:func:`require` (a function call, fires under ``-O``) or ``np.testing.*``
— NEVER a bare ``assert``.

vv Mode-11: the cross-check gates MUST execute the NEW ``kernel`` /
``production_rate`` property, not a sibling path. Each cross-check reads
the property OFF the live operator (``op.production_rate`` /
``op.kernel``) and compares against an independently-built reference —
so mutating the property reddens the gate.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest

if TYPE_CHECKING:
    from orpheus.numerics.functional import Functional
    from orpheus.numerics.operator import LinearOperator


# ───────────────────────────────────────────────────────────────────────
# vv Mode-8: -O-firing assertion (NOT a bare ``assert``).
# ───────────────────────────────────────────────────────────────────────


def require(condition: bool, message: str) -> None:
    """Fail with ``message`` if ``condition`` is false. Fires under ``-O``."""
    if not condition:
        pytest.fail(message)


# ───────────────────────────────────────────────────────────────────────
# Not-yet-written-SUT import guards (probe + skip).
# ───────────────────────────────────────────────────────────────────────


def require_integral_kernel_operator():
    """The ``IntegralKernelOperator`` Protocol (S6) — skip if not yet landed.

    The brief sketches ``orpheus.transport.operators.integral_kernel_operator`` (the
    numerics/transport-floor home alongside the S5 ``Functional``). The
    exact module path is the method-implementer's latitude; probe the
    likely homes and skip with a clear reason otherwise.
    """
    for path in (
        "orpheus.transport.operators.integral_kernel_operator",
        "orpheus.numerics.integral_kernel_operator",
        "orpheus.transport.kernel_operator",
    ):
        try:
            mod = __import__(path, fromlist=["IntegralKernelOperator"])
        except ImportError:
            continue
        cls = getattr(mod, "IntegralKernelOperator", None)
        if cls is not None:
            return cls
    pytest.skip(
        "#257 S6 PRE-IMPL: IntegralKernelOperator Protocol not yet written "
        "(probed transport.integral_kernel_operator / "
        "numerics.integral_kernel_operator / transport.kernel_operator)."
    )


def require_production_rate_property(fission_op) -> "Functional":
    """The NEW ``FissionOperator.production_rate`` property (S6) — skip if absent.

    THE single construction assumption for the fission ``production_rate``
    surface — a property named ``production_rate`` returning a
    :class:`~orpheus.numerics.functional.Functional`. If the
    method-implementer chooses a different name, change this ONE function.
    The probe is duck-typed (``hasattr``) so the gate skips (PRE-IMPL),
    rather than raising, until the property lands.
    """
    if not hasattr(fission_op, "production_rate"):
        pytest.skip(
            "#257 S6 PRE-IMPL: FissionOperator.production_rate property not "
            "yet written by the method-implementer."
        )
    return fission_op.production_rate


def require_scattering_kernel_property(scattering_op) -> "LinearOperator":
    """The NEW ``ScatteringOperator.kernel`` property (S6) — skip if absent.

    THE single construction assumption for the scattering ``kernel``
    surface — a property named ``kernel`` returning a
    :class:`~orpheus.numerics.operator.LinearOperator` (the typed
    ``R∘Λ∘M`` composition). If the method-implementer chooses a different
    name, change this ONE function.
    """
    if not hasattr(scattering_op, "kernel"):
        pytest.skip(
            "#257 S6 PRE-IMPL: ScatteringOperator.kernel property not yet "
            "written by the method-implementer."
        )
    return scattering_op.kernel


# ───────────────────────────────────────────────────────────────────────
# The STRUCTURALLY-INDEPENDENT reference for the FISSION kernel cross-check
# (Part B). χ · (νΣf · φ contracted over groups), by an EXPLICIT Python
# double-loop — shares NO numpy reduction primitive with the production
# ``RankOneOperator.apply`` / ``ProductionRateFunctional.evaluate`` path.
# ───────────────────────────────────────────────────────────────────────


def hand_derived_fission_emission(
    chi: np.ndarray, nu_sigma_f: np.ndarray, phi: np.ndarray,
) -> np.ndarray:
    r"""Per-group fission emission :math:`(F\phi)_g` by an EXPLICIT Python loop.

    .. math::

        (F\,\phi)_g(\vec r) \;=\; \chi_g(\vec r)\,
            \sum_{g'} \nu\Sigma_{f,g'}(\vec r)\,\phi_{g'}(\vec r)

    Computed by an explicit nested Python ``for`` loop with a running
    scalar accumulator per cell — NOT ``numpy.einsum``, NOT
    ``.sum(axis=...)``, NOT a numpy broadcast-multiply. This is the
    structurally-independent reference for the fission cross-check: it
    shares NO reduction primitive with the production
    ``inner = (right * x).sum(axis=0, keepdims=True)`` extraction inside
    :meth:`RankOneOperator.apply` / :meth:`ProductionRateFunctional.evaluate`
    and NONE of the ORPHEUS operator algebra. A νΣf↔φ swap, a wrong-axis
    contraction, or a dropped χ broadcast in the production code disagrees
    with this loop.

    Parameters
    ----------
    chi, nu_sigma_f, phi : np.ndarray
        ``(ng, *spatial)`` per-group per-cell arrays. ``chi`` is the
        emission spectrum, ``nu_sigma_f`` the production cross section,
        ``phi`` the flux.

    Returns
    -------
    np.ndarray
        Per-group emission ``(ng, *spatial)`` (the production density
        re-broadcast across groups by χ).
    """
    ng = nu_sigma_f.shape[0]
    spatial_shape = nu_sigma_f.shape[1:]
    out = np.zeros((ng, *spatial_shape), dtype=float)
    for flat in range(int(np.prod(spatial_shape))):
        idx = np.unravel_index(flat, spatial_shape)
        # production density at this cell (group contraction), by hand.
        acc = 0.0
        for gp in range(ng):
            acc += float(nu_sigma_f[(gp, *idx)]) * float(phi[(gp, *idx)])
        # re-broadcast across groups weighted by the emission spectrum.
        for g in range(ng):
            out[(g, *idx)] = float(chi[(g, *idx)]) * acc
    return out
