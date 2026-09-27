r"""The verification floor for an SN-against-reference comparison, and the heterogeneous cylinder reference it is applied to.

A comparison of a production solve with a reference certifies the solve only
when the reference's own error bound on the compared observable is at most a
tenth of the tolerance (the floor of ``.claude/plans/reference_cache.md``,
"The architecture as it stands", point 6). :func:`certify_agreement` returns
that comparison as a value, an :class:`AgreementCertificate`, whose
:meth:`~AgreementCertificate.require` raises for the floor first and the
agreement second, so a row whose reference cannot certify its bound is red
for that reason and not for the reading. It is the test-side stand-in for the
comparison verb that phase P2 of that plan gives the ``ReferenceSolution``,
and retires when that verb lands.

This module also owns the ONE definition of the heterogeneous closed-cylinder
reference problem the SN rows compare against (fuel A | moderator B | fuel A
at outer radii 0.5, 1.5, 2.0 cm, 2 groups, reflective), its cached solve, its
(absent) bound, the strict-xfail mark every row that needs the bound carries,
and the RECORD values those rows' companions pin. The ladders every number
comes from, and the command that re-measures them, are
:mod:`tests.gates.derivations._trajectory_resolvent_ladders`.
"""

from __future__ import annotations

import functools
from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_xs
from orpheus.derivations.continuous.trajectory_resolvent.greens_function_cylinder import (
    CylinderGreensMRResult,
    solve_greens_function_cylinder_mr,
)

# ── the verification floor ───────────────────────────────────────────────


@dataclass(frozen=True)
class AgreementCertificate:
    """One comparison of a production reading with a reference: the values, the bound, the tolerance, the floor verdict."""

    observable: str
    reading: float
    tolerance: float
    reference_bound: float | None

    @property
    def floor_holds(self) -> bool:
        """The reference carries a bound, and it is at most a tenth of the tolerance."""
        return self.reference_bound is not None and self.reference_bound <= self.tolerance / 10

    @property
    def agrees(self) -> bool:
        """The floor holds and the reading is below the tolerance."""
        return self.floor_holds and self.reading < self.tolerance

    def require(self) -> None:
        """Raise ``AssertionError`` unless the floor holds, then unless the reading agrees."""
        if self.reference_bound is None:
            raise AssertionError(
                f"{self.observable}: the reference carries no certified error bound, so this "
                f"comparison cannot certify the SN solve (today's reading {self.reading:.3e}, "
                f"tolerance {self.tolerance:.0e})"
            )
        if not self.floor_holds:
            raise AssertionError(
                f"{self.observable}: the reference's error bound {self.reference_bound:.1e} exceeds "
                f"a tenth of the tolerance {self.tolerance:.0e}, so this comparison cannot certify "
                f"the SN solve (today's reading {self.reading:.3e})"
            )
        if not self.agrees:
            raise AssertionError(
                f"{self.observable}: SN against the reference reads {self.reading:.3e}, above the "
                f"tolerance {self.tolerance:.0e} (reference bound {self.reference_bound:.1e})"
            )


def certify_agreement(
    observable: str, reading: float, tolerance: float, reference_bound: float | None,
) -> AgreementCertificate:
    """The comparison as a value; the row calls :meth:`AgreementCertificate.require`."""
    return AgreementCertificate(observable, float(reading), float(tolerance), reference_bound)


def assert_record(readings: Mapping[str, float], recorded: Mapping[str, float], band: float, *, relative: frozenset[str]) -> None:
    """RECORD: each recorded reading still reads within ``band`` (relative for the names in ``relative``, absolute otherwise).

    Not verification: a RECORD row pins what the code prints today so that a
    change on either side of a comparison that cannot yet certify reddens it.
    It raises (never a bare ``assert``: this module is a helper, not a
    collected test module, so ``python -O`` would strip one).
    """
    for name, value in recorded.items():
        scale = abs(value) if name in relative else 1.0
        if not abs(readings[name] - value) <= band * scale:
            raise AssertionError(
            f"record {name!r} moved: {readings[name]!r} against the recorded {value!r}. If the "
            f"reference was repaired (#516), re-derive its bound, re-measure the records "
            f"(python -O -m tests.gates.derivations._trajectory_resolvent_ladders records) and "
            f"lift the xfails; otherwise an SN change moved the solve"
            )


# ── the heterogeneous closed cylinder ─────────────────────────────────────

ABA_RADII = np.array([0.5, 1.5, 2.0])
_ABA_KEYS = ("A", "B", "A")


def aba_xs_2g() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """``sigma_t (3,2)``, ``sigma_s (3,2,2)`` (``[region, from, to]``), ``nu_sigma_f (3,2)``, ``chi (3,2)``."""
    parts = [get_xs(k, "2g") for k in _ABA_KEYS]
    return (
        np.stack([p["sig_t"] for p in parts]),
        np.stack([p["sig_s"] for p in parts]),
        np.stack([p["nu"] * p["sig_f"] for p in parts]),
        np.stack([p["chi"] for p in parts]),
    )


@functools.cache
def cylinder_3reg_reference() -> CylinderGreensMRResult:
    """The cylinder reference at (n_r, n_mu_axial, n_phi_az) = (24, 16, 32), solved once per session."""
    sigma_t, sigma_s, nu_sigma_f, chi = aba_xs_2g()
    return solve_greens_function_cylinder_mr(
        radii=ABA_RADII, sigma_t=sigma_t, sigma_s=sigma_s,
        nu_sigma_f=nu_sigma_f, chi=chi, alpha=1.0,
        n_r=24, n_mu_axial=16, n_phi_az=32, n_traj_quad=64,
        max_iter=2000, tol=1e-9, initial_k=1.23,
    )


#: NO certified bound (#516): the reference's ladder
#: (:data:`~tests.gates.derivations._trajectory_resolvent_ladders.CYLINDER_3REG_REFERENCE_K`)
#: is not monotone in the azimuthal order and couples to the radial one, so
#: no finite ladder bounds it. The rows whose bound needs it therefore fail
#: the floor, and they XPASS only when this ``None`` is replaced by a bound
#: derived from a ladder that converges (after #516's repair), not when the
#: repair alone lands.
CYLINDER_3REG_REFERENCE_BOUND: dict[str, float | None] = {"k": None, "shape": None}

#: The mark every row that needs the cylinder bound carries: strict (an XPASS
#: fails) and expecting only the floor's ``AssertionError`` (a crash before it
#: is not absorbed; ``vv-principles`` mode 8(4)).
awaits_cylinder_bound = pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason=(
        "#516: the cylinder trajectory-resolvent reference carries no certified error "
        "bound (its azimuthal Gauss-Legendre rule meets the tangency kinks of the "
        "interior interfaces); the floor assertion is the expected failure"
    ),
)

#: RECORD values of the cylinder comparisons, ``[M]`` 2026-09-26 with
#: ``python -O -m tests.gates.derivations._trajectory_resolvent_ladders records``.
#: ``k_ref`` is :func:`cylinder_3reg_reference`; ``phase_c_*`` are the live
#: folded-16x32 SN solve and its gaps (``test_phase_c_crosscheck``);
#: ``unified_k`` is the folded-4x8 Krylov-on-(L + C) solve
#: (``test_unified_matvec_cylinder``); ``standoff_sweep_k_nx40`` the folded-4x8
#: source-iteration solve on 40 uniform cells (``test_l1_standoff_slab_cylinder``).
CYLINDER_3REG_RECORD: dict[str, float] = {
    "k_ref": 1.231036749830859,
    "phase_c_k_sn": 1.2317720792844793,
    "phase_c_k_gap": 0.000596968762311497,
    "phase_c_shape_fast": 0.005314982320667484,
    "phase_c_shape_thermal": 0.0008352414999330727,
    "unified_k": 1.2310184196907399,
    "standoff_sweep_k_nx40": 1.23101841974857,
}
#: The eigenvalues are compared relatively, the gaps absolutely.
CYLINDER_3REG_RECORD_RELATIVE = frozenset({"k_ref", "phase_c_k_sn", "unified_k", "standoff_sweep_k_nx40"})
#: Ten times the looser iterative floor: the 4x8 SN solves stop on a relative
#: k increment of 1e-7 with a dominance ratio near 0.95, so each k is within
#: about 2e-6 of its fixed point (the reference at tol 1e-9 and the 16x32 SN
#: solve at keff_tol 1e-12 are far tighter).
CYLINDER_3REG_RECORD_BAND = 2e-5
