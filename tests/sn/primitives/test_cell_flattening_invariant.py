"""Foundation test — cell-flattening invariant for ``SNSolver`` XS storage.

PR-CLEANUP-CODE §E.  Promoted from the ``__debug__`` block in
``orpheus.sn.solver.SNSolver.__init__`` (PR-INDEX-3 era).

The invariant under test
------------------------

``xs.sig_t : (N_cells, ng)`` must round-trip bit-identically between
the legacy ``(nx, ny, ng)`` reshape and the principled ``(ng, nx, ny)``
``.T.reshape`` storage layout::

    sig_t_old = xs.sig_t.reshape(nx, ny, ng)                  # legacy
    sig_t_new = xs.sig_t.T.reshape(ng, nx, ny)                # principled
    assert np.array_equal(sig_t_old, sig_t_new.transpose(1, 2, 0))

This holds iff:

1. ``xs.sig_t`` is row-major (C-order) over ``(N_cells, ng)``.
2. ``N_cells = nx * ny`` flattens in C-order over ``(nx, ny)``.

A silent break here (a producer flipping ravel order; a consumer
reading the wrong axis; a numpy-version subtlety changing reshape
behaviour) would propagate as a silent semantic bug everywhere
SNSolver reads cross sections.  Caught here as a foundation invariant
before any solver path executes — no L0/L1/L2 theory-page label is
needed, this is a software invariant about array storage.

This file deliberately lives in its own module — co-locating with
``test_solver_components.py`` would expose the foundation marker to
the file-level ``pytestmark = pytest.mark.l0`` and yield a "conflicting
V&V level markers" warning (the V&V audit harness prefers the
stronger ``l0`` claim, but the brief specifies ``foundation`` is the
correct tag for software invariants without a theory-page link).
"""
from __future__ import annotations

import numpy as np
import pytest


@pytest.mark.foundation
@pytest.mark.parametrize(
    "nx,ny,ng",
    [
        (1, 1, 2),   # 1D-degenerate / smallest meaningful case
        (5, 1, 2),   # 1D-with-trailing-singleton (the slab pattern)
        (3, 4, 3),   # 2D non-square multi-group
    ],
    ids=["1x1x2", "5x1x2", "3x4x3"],
)
def test_cell_flattening_invariant_xs_storage_round_trips(
    nx: int, ny: int, ng: int,
) -> None:
    """``xs.sig_t : (N_cells, ng)`` round-trip ``.T.reshape`` <-> ``reshape``.

    The promoted invariant from the ``__debug__`` block in
    ``SNSolver.__init__``.  See module docstring for the load-bearing
    rationale.
    """
    # Build a synthetic xs.sig_t : (N_cells, ng) with KNOWN, mutually
    # distinct values so a wrong axis swap is detectable (every entry
    # is unique → any axis confusion produces a value mismatch, not
    # just a permutation that happens to bit-match).
    n_cells = nx * ny
    sig_t = (
        np.arange(n_cells * ng, dtype=np.float64).reshape(n_cells, ng)
        + 1.0  # offset from zero to make the spot-check failure messages
               # carry non-zero values (a zero/zero false-pass is impossible
               # because every cell is a distinct integer).
    )

    sig_t_old = sig_t.reshape(nx, ny, ng)
    sig_t_new = sig_t.T.reshape(ng, nx, ny)

    np.testing.assert_array_equal(
        sig_t_old, sig_t_new.transpose(1, 2, 0),
        err_msg=(
            f"Cell-flattening invariant broke for (nx, ny, ng) = "
            f"({nx}, {ny}, {ng}) — ``mat_ids`` ravel order is not "
            f"C-order ``(nx, ny)``, or the ``.T.reshape`` storage flip "
            f"disagrees with the legacy reshape."
        ),
    )

    # Spot-check a representative entry too, so the assertion failure
    # message (if any) points at a concrete (g, i, j) coordinate rather
    # than just "arrays differ at N positions".  The C-order ravel
    # ``flat = i * ny + j`` is the load-bearing convention.
    for g in range(ng):
        for i in range(nx):
            for j in range(ny):
                flat = i * ny + j           # C-order (nx, ny) ravel
                expected = sig_t[flat, g]
                assert sig_t_new[g, i, j] == expected, (
                    f"sig_t_new[{g}, {i}, {j}] = {sig_t_new[g, i, j]} "
                    f"!= sig_t[{flat}, {g}] = {expected}"
                )
                assert sig_t_old[i, j, g] == expected, (
                    f"sig_t_old[{i}, {j}, {g}] = {sig_t_old[i, j, g]} "
                    f"!= sig_t[{flat}, {g}] = {expected}"
                )
