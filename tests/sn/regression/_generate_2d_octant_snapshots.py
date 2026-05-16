"""Generate the 2-D octant-sweep equivalence snapshots from LEGACY code.

Run::

    python -m tests.sn.regression._generate_2d_octant_snapshots
    python -m tests.sn.regression._generate_2d_octant_snapshots --case 03_l7_trap_mixedBC_2g_het_LS4
    python -m tests.sn.regression._generate_2d_octant_snapshots --list

Each snapshot writes ``snapshots/2d_octant_equivalence_<case_id>.npz``
containing (Issue #196 PR-INDEX-5 — principled ``g`` after ``N``):

* ``angular_flux`` — ``(N, ng, nx, ny)`` float64
* ``scalar_flux`` — ``(ng, nx, ny)`` float64
* ``psi_x_post`` — ``(N, ng, nx+1, ny)`` float64
* ``psi_y_post`` — ``(N, ng, nx, ny+1)`` float64
* ``case_id`` — np.array(case.case_id)
* ``case_description`` — np.array(case.description)
* ``failure_mode`` — np.array(case.failure_mode)
* ``generator_commit`` — short SHA

**This script MUST be run on the LEGACY (pre-Wave-2) commit**, before
the ``_sweep_2d_wavefront`` refactor lands.  The companion test at
:file:`tests/sn/test_2d_octant_sweep_equivalence.py` then re-runs each
case against whatever implementation is current and asserts
bit-for-bit (within ``nulp=64``) agreement.

Parity with :mod:`tests.sn.regression._generate_snapshots`:

* Same ``snapshots/`` directory.
* Same ``--case`` / ``--list`` CLI.
* Same ``generator_commit`` metadata.

Drift protocol — when this generator's output legitimately changes
(e.g. an upstream Mesh2D refactor changes the canonical edges/center
spacing): (1) audit why the new output is correct, (2) re-run the
generator, (3) commit both new snapshots AND the audit narrative in
the same commit.
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import numpy as np

from tests.sn.test_2d_octant_sweep_equivalence import (
    CASES,
    SNAPSHOT_DIR,
    OctantEquivalenceCase,
    _snapshot_path,
)
from orpheus.sn.sweep import _sweep_2d_wavefront


def _git_short_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short=12", "HEAD"],
            stderr=subprocess.DEVNULL,
        ).decode().strip()
    except Exception:
        return "unknown"


def generate_one(
    case: OctantEquivalenceCase, *, sha: str | None = None,
) -> Path:
    """Run the case under the LEGACY sweep and write the .npz snapshot."""
    inputs = case.builder()

    # The legacy ``_sweep_2d_wavefront`` mutates the persistent
    # buffers in-place; we capture the post-sweep buffers so the test
    # can assert agreement on the stateful output as well as on the
    # angular/scalar flux.  Issue #197 PR-TYPED-2: the typed
    # :class:`BoundaryFlux` exposes the two persistent buffers as
    # ``xmin_xmax_buf`` and ``ymin_ymax_buf``.
    angular_flux = scalar_flux = None
    for _ in range(case.n_sweeps):
        angular_flux, scalar_flux = _sweep_2d_wavefront(
            inputs.Q, inputs.sig_t, inputs.sn_mesh, inputs.boundary_flux,
            Q_aniso=inputs.aniso_source,
        )
    psi_x_post = inputs.boundary_flux.xmin_xmax_buf
    psi_y_post = inputs.boundary_flux.ymin_ymax_buf

    SNAPSHOT_DIR.mkdir(parents=True, exist_ok=True)
    out = _snapshot_path(case.case_id)

    payload = dict(
        angular_flux=np.asarray(angular_flux, dtype=np.float64),
        scalar_flux=np.asarray(scalar_flux, dtype=np.float64),
        psi_x_post=np.asarray(psi_x_post, dtype=np.float64),
        psi_y_post=np.asarray(psi_y_post, dtype=np.float64),
        case_id=np.array(case.case_id),
        case_description=np.array(case.description),
        failure_mode=np.array(case.failure_mode),
        generator_commit=np.array(sha or _git_short_sha()),
    )
    np.savez_compressed(out, **payload)
    return out


def generate_all(case_ids: list[str] | None = None) -> list[Path]:
    sha = _git_short_sha()
    written = []
    for case in CASES:
        if case_ids and case.case_id not in case_ids:
            continue
        path = generate_one(case, sha=sha)
        written.append(path)
        print(f"wrote  {path.relative_to(Path.cwd())}")
    return written


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Generate 2-D octant-sweep equivalence snapshots. "
            "Run on the LEGACY commit BEFORE the Wave-2 refactor."
        ),
    )
    parser.add_argument(
        "--case", action="append", default=None,
        help="Restrict to a specific case_id (repeatable).",
    )
    parser.add_argument(
        "--list", action="store_true",
        help="List available cases and exit.",
    )
    args = parser.parse_args()
    if args.list:
        for case in CASES:
            print(f"  {case.case_id:50s}  {case.description}")
        return
    written = generate_all(case_ids=args.case)
    print(f"generated {len(written)} snapshot(s) in {SNAPSHOT_DIR}")


if __name__ == "__main__":
    main()
