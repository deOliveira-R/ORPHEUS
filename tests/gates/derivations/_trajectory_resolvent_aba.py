r"""The OLD reference's spelling on the SN rows' problems, for the corroboration rows until step (e) deletes the family.

Until P1 step (d) of ``.claude/plans/characteristic_reference_architecture.md``
the SN rows read the trajectory resolvent through
``tests/gates/sn/verification/analytical/_aba_reference.py`` (the A|B|A
bodies) and ``test_partial_reflector_resolvent.py`` (the ERR-094 bodies). Step
(d) re-pointed both onto the characteristic reference, so the old spelling
lives here, on the OLD side of
``tests/gates/derivations/test_characteristic_reference_corroboration.py``:
this module imports the old family and nothing that reaches the new
reference, which that file's independence legs assert (it is one of its
``_OLD_HELPERS``). The A|B|A specification is therefore an ARGUMENT, built by
the caller from ``_aba_reference`` (the problem is the shared input; ``spec
§9``, X4), never imported here.
"""
from __future__ import annotations

from collections.abc import Mapping

import numpy as np

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.derivations.continuous.trajectory_resolvent.reference import trajectory_resolvent_reference
from orpheus.geometry import CoordSystem
from orpheus.reference.solution import ReferenceSolution
from orpheus.specification.specification import GeometrySpecification

#: The old reference's resolution on each A|B|A body (the fixtures of the 2026-09-26 ladders,
#: :mod:`tests.gates.derivations._trajectory_resolvent_ladders`).
ABA_REFERENCE_QUADRATURE = {
    CoordSystem.SPHERICAL: {"n_r": 36, "n_mu": 96, "n_traj_quad": 64},
    CoordSystem.CYLINDRICAL: {"n_r": 24, "n_mu_axial": 16, "n_phi_az": 32, "n_traj_quad": 64},
}
_INITIAL_K = {CoordSystem.SPHERICAL: 1.38, CoordSystem.CYLINDRICAL: 1.23}
#: The old power iteration's settings on the A|B|A bodies.
ABA_SOLVE_TOL = 1e-9
_ABA_MAX_ITER = 2000


def aba_reference_at(specification: GeometrySpecification, quadrature: Mapping[str, int], coord: CoordSystem) -> ReferenceSolution:
    """The trajectory-resolvent reference on an A|B|A ``specification`` at ``quadrature``: lazy, uncertified (#566;
    the cylinder also #516). The one spelling of its power iteration's settings, for the rows and the ladders alike."""
    return trajectory_resolvent_reference(
        specification, quadrature, max_iter=_ABA_MAX_ITER, tol=ABA_SOLVE_TOL, initial_k=_INITIAL_K[coord],
    )


def partial_reflector_cross_sections() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Mixture A as the old ERR-094 solvers read it: the library's two-group mixture's P0 transfer ``[g_from, g_to]``,
    with its total, production and spectrum (the SN rows solve it at scattering order 0)."""
    fuel = get_mixture("A", "2g")
    return (
        np.asarray(fuel.SigT, dtype=float), fuel.SigS[0].toarray(),
        np.asarray(fuel.SigP, dtype=float), np.asarray(fuel.chi, dtype=float),
    )
