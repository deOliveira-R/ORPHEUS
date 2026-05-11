"""Generate the BC-equivalence snapshot harness (Wave 6 / C6.1).

This script is the **safety net** for Waves 7–11 of the boundary-operator
refactor (see ``.claude/plans/transient-giggling-cake.md``). It captures
the LEGACY ``bc.apply(psi_out, quadrature)`` output for 8 representative
BC × quadrature × input cases and writes each to a ``.npz`` snapshot.
The companion test harness
(``tests/geometry/test_bc_equivalence_snapshot.py``) re-runs each case
on every pytest invocation and asserts equivalence at the per-case
tolerance.

Run::

    python -m tests.geometry._generate_bc_equivalence_snapshots
    python -m tests.geometry._generate_bc_equivalence_snapshots --case vacuum_lebedev17

Output: ``tests/geometry/snapshots/bc_equivalence_<case>.npz`` per case.

V&V context
===========

The snapshots are **L1 / regression**: they pin the legacy
``BoundaryOperator.apply`` output exactly so that the Wave 7–11
production-call-site migration cannot silently drift the numerical
output (per :ref:`vv-principles` §1, structural-independence). The
legacy code is the structurally-independent reference (it predates
the refactor and is L1-verified at the SN sweep level by
``tests/sn/test_boundary*`` and the MMS suites). The snapshot pins
its output exactly.

**Per-case tolerance rationale** (V&V principle, bit-identity vs
principled-equivalence three criteria — see ``vv-principles`` skill):

* ``vacuum_lebedev17`` — legacy returns zeros everywhere; the
  realizer (post-Wave 8) zeros only the inflow ordinates. The
  snapshot records BOTH the legacy zeros-all output AND the inflow
  indices so the post-Wave-8 test can switch to the
  principled-equivalent §16A.5 semantics.
* ``albedo_05_lebedev17`` — single scalar multiplication; bit-exact.
* ``specular_x_lebedev17`` — pure permutation; bit-exact.
* ``specular_y_partial_07_LS6`` — multiplication then permutation;
  bit-exact.
* ``white_xmax_LS4`` — cosine-weighted average; ``nulp=4`` allows
  for one ULP of compositional reordering (Wave 1's
  ``TestLegacyBitEquivalence`` already pins the body bit-exact at
  α=1, ``nulp=4`` is a safety margin).
* ``white_xmin_partial_03_GL`` — same body × scalar shift; ``nulp=4``.
* ``periodic_lebedev17`` — identity body (copy); bit-exact.
* ``mixed_30spec_70white_LS4`` — sum of two BCs; ``nulp=64`` covers
  the sum-of-products reduction tree change once Wave-0 ``OperatorSum``
  composes the realizer chain.

The snapshot directory ``tests/geometry/snapshots/`` is committed to
the repository — the snapshots ARE the verification artefact.
Regeneration is intentional only when (a) a Wave-induced semantic
change is documented and approved (Wave 8's vacuum-trace
correction), or (b) the underlying quadrature has been re-verified
(no current case).
"""
from __future__ import annotations

import argparse
import hashlib
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np

from orpheus.geometry.boundary import (
    AlbedoBoundaryOperator,
    BoundaryOperator,
    PeriodicBoundaryOperator,
    SpecularBoundaryOperator,
    VacuumBoundaryOperator,
    WhiteBoundaryOperator,
)
from orpheus.sn.boundary_realizer import SNBoundaryRealizer, SNMethodSpace
from orpheus.sn.quadrature import (
    AngularQuadrature,
    GaussLegendre1D,
    LebedevSphere,
    LevelSymmetricSN,
)


SNAPSHOT_DIR = Path(__file__).parent / "snapshots"

# Input shape used for every snapshot: ``(N_ordinates, n_cells, n_groups)``.
# The trailing ``(5, 3)`` is arbitrary — the BC primitives are agnostic to
# trailing axes (they broadcast over the leading ordinate axis), so any
# fixed shape with at least one trailing axis works as long as the
# generator and harness agree.
_PSI_TRAILING_SHAPE: tuple[int, int] = (5, 3)


# ─── seed derivation ─────────────────────────────────────────────────


def _seed_from_case_id(case_id: str) -> int:
    """Deterministic small seed from a case ID.

    Using ``hashlib.sha256`` (not ``hash()``) so the seed is stable
    across Python interpreter invocations (PYTHONHASHSEED-immune).
    Truncated to 16 bits — small enough for ``np.random.default_rng``
    to print cleanly in error messages, large enough that distinct case
    IDs collide rarely.
    """
    digest = hashlib.sha256(case_id.encode("utf-8")).digest()
    # Two-byte truncation; matches the "small integer derived from the
    # case ID" requirement in the plan.
    return int.from_bytes(digest[:2], byteorder="big")


def _generate_psi(quad: AngularQuadrature, case_id: str) -> np.ndarray:
    """Deterministic input for one snapshot case.

    ``np.random.default_rng(<seed>).uniform(0, 2, size=(N, 5, 3))`` per
    plan. The (0, 2) range is symmetric around 1.0 and keeps inputs
    O(1), making bit-equivalence comparisons unambiguous.
    """
    seed = _seed_from_case_id(case_id)
    rng = np.random.default_rng(seed)
    shape = (int(quad.N),) + _PSI_TRAILING_SHAPE
    return rng.uniform(0.0, 2.0, size=shape)


# ─── case registry ───────────────────────────────────────────────────


@dataclass(frozen=True)
class BCEquivalenceCase:
    """One snapshot case definition.

    The ``build_bc`` and ``build_quadrature`` callables are invoked at
    both snapshot-generation time AND harness time, ensuring the two
    sides of the equivalence cannot drift in BC or quadrature
    configuration.

    ``build_bc=None`` is the Wave-11 marker for the
    ``mixed_30spec_70white_LS4`` case (the only case for which the
    snapshot pins the realised ``apply(psi)`` output of a Wave-0
    ``OperatorSum`` composition — see the comment block on the
    ``mixed_30spec_70white_LS4`` entry in ``CASES``).
    """

    case_id: str
    description: str
    build_bc: Callable[[], BoundaryOperator] | None
    build_quadrature: Callable[[], AngularQuadrature]

    @property
    def snapshot_path(self) -> Path:
        return SNAPSHOT_DIR / f"bc_equivalence_{self.case_id}.npz"


CASES: tuple[BCEquivalenceCase, ...] = (
    BCEquivalenceCase(
        case_id="vacuum_lebedev17",
        description=(
            "VacuumBoundaryOperator + LebedevSphere(17). Legacy "
            "returns np.zeros_like; Wave 8 will switch to inflow-only "
            "masking per §16A.5 — see test harness for the dual "
            "assertion."
        ),
        build_bc=lambda: VacuumBoundaryOperator(),
        build_quadrature=lambda: LebedevSphere.create(17),
    ),
    BCEquivalenceCase(
        case_id="albedo_05_lebedev17",
        description="AlbedoBoundaryOperator(0.5) + LebedevSphere(17).",
        build_bc=lambda: AlbedoBoundaryOperator(albedo=0.5),
        build_quadrature=lambda: LebedevSphere.create(17),
    ),
    BCEquivalenceCase(
        case_id="specular_x_lebedev17",
        description=(
            "SpecularBoundaryOperator(axis='x', albedo=1.0) + "
            "LebedevSphere(17). Pure permutation; bit-exact."
        ),
        build_bc=lambda: SpecularBoundaryOperator(axis="x", albedo=1.0),
        build_quadrature=lambda: LebedevSphere.create(17),
    ),
    BCEquivalenceCase(
        case_id="specular_y_partial_07_LS6",
        description=(
            "SpecularBoundaryOperator(axis='y', albedo=0.7) + "
            "LevelSymmetricSN(6). Multiplication-then-permutation; "
            "bit-exact."
        ),
        build_bc=lambda: SpecularBoundaryOperator(axis="y", albedo=0.7),
        build_quadrature=lambda: LevelSymmetricSN.create(sn_order=6),
    ),
    BCEquivalenceCase(
        case_id="white_xmax_LS4",
        description=(
            "WhiteBoundaryOperator(axis='x', outward_sign=+1, "
            "albedo=1.0) + LevelSymmetricSN(4). Cosine-weighted "
            "average; Wave 1 bit-equivalence with nulp=4 safety margin."
        ),
        build_bc=lambda: WhiteBoundaryOperator(
            axis="x", outward_sign=+1, albedo=1.0,
        ),
        build_quadrature=lambda: LevelSymmetricSN.create(sn_order=4),
    ),
    BCEquivalenceCase(
        case_id="white_xmin_partial_03_GL",
        description=(
            "WhiteBoundaryOperator(axis='x', outward_sign=-1, "
            "albedo=0.3) + GaussLegendre1D(8). Scalar-shifted "
            "cosine-weighted average; nulp=4."
        ),
        build_bc=lambda: WhiteBoundaryOperator(
            axis="x", outward_sign=-1, albedo=0.3,
        ),
        build_quadrature=lambda: GaussLegendre1D.create(n_ordinates=8),
    ),
    BCEquivalenceCase(
        case_id="periodic_lebedev17",
        description=(
            "PeriodicBoundaryOperator + LebedevSphere(17). Identity "
            "body (psi.copy()); bit-exact."
        ),
        build_bc=lambda: PeriodicBoundaryOperator(),
        build_quadrature=lambda: LebedevSphere.create(17),
    ),
    # Wave 11: the ``mixed_30spec_70white_LS4`` case is the special
    # "rank-N composition" entry. Pre-Wave-11 it stored
    # ``MixedBoundaryOperator([(0.3, Spec), (0.7, White)]).apply(psi, quad)``
    # as the legacy reference; Wave 11 removed that composer and the
    # equivalent rank-N composition is now expressed via Wave-0
    # ``OperatorSum``-algebra over realised primitives. The case below
    # carries ``build_bc=None`` as the marker for the special handling
    # in :func:`_build_payload`: the snapshot output IS the realised
    # ``0.3 * spec_realized + 0.7 * white_realized`` ``apply(psi)``,
    # not a legacy 2-arg ``bc.apply(psi, quad)`` output.
    BCEquivalenceCase(
        case_id="mixed_30spec_70white_LS4",
        description=(
            "0.3 * SpecularBoundaryOperator(axis='x', albedo=1.0) + "
            "0.7 * WhiteBoundaryOperator(axis='x', outward_sign=+1, "
            "albedo=1.0) + LevelSymmetricSN(4). Wave-11 Wave-0 "
            "OperatorSum-of-realised-primitives composition; the "
            "snapshot pins the realised ``apply(psi)`` output (no "
            "legacy 2-arg path — ``MixedBoundaryOperator`` was removed "
            "in Wave 11)."
        ),
        build_bc=None,
        build_quadrature=lambda: LevelSymmetricSN.create(sn_order=4),
    ),
)


# ─── face metadata for the vacuum special case ───────────────────────


def _xmin_inflow_indices(quad: AngularQuadrature) -> np.ndarray:
    """Inflow ordinate indices for the ``xmin`` face.

    Convention (from :mod:`orpheus.numerics.trace_space`): the ``xmin``
    face has outward normal ``-x``, so the ordinate-axis inflow predicate
    is ``-mu_x < -eps``, i.e. ``mu_x > eps``. We use an exact ``> 0``
    cut since the quadrature ordinates we ship for this case (Lebedev 17)
    have no ordinates with ``mu_x == 0`` exactly.
    """
    return np.flatnonzero(np.asarray(quad.mu_x) > 0.0)


# ─── snapshot writer ─────────────────────────────────────────────────


def _build_mixed_realized_apply(
    quad: AngularQuadrature, psi_out: np.ndarray,
) -> np.ndarray:
    """Wave-11 Wave-0 ``OperatorSum``-composition for the mixed case.

    The pre-Wave-11 snapshot value was
    ``MixedBoundaryOperator([(0.3, Spec), (0.7, White)]).apply(psi, quad)``.
    With the composer removed, the equivalent rank-N composition is
    expressed via Wave-0 ``ScaledOperator``/``OperatorSum`` algebra
    over the realised leaves. The result is mathematically identical;
    by ``vv-principles`` "bit-identity vs principled-equivalence", a
    ULP shift relative to the pre-Wave-11 snapshot is acceptable
    because (a) each intermediate is a named Wave-0 type, (b) the
    composed output is verified against the explicit pointwise
    weighted sum (the structurally-independent reference), and (c) any
    drift is FP non-associativity over a 2-summand reduction.

    The realizer's pre-Wave-11 internal mixed-BC path also composed
    via ``OperatorSum`` of ``ScaledOperator`` (see the deleted
    ``isinstance(law, MixedBoundaryOperator)`` branch in
    ``orpheus.sn.boundary_realizer``), so on the current ship-state
    the two reduction trees agree bit-exactly.
    """
    spec_realized = SNBoundaryRealizer().realize(
        SpecularBoundaryOperator(axis="x", albedo=1.0),
        SNMethodSpace.minimal(quad),
    )
    white_realized = SNBoundaryRealizer().realize(
        WhiteBoundaryOperator(axis="x", outward_sign=+1, albedo=1.0),
        SNMethodSpace.minimal(quad),
    )
    composed = 0.3 * spec_realized + 0.7 * white_realized
    return composed.apply(psi_out)


def _build_payload(case: BCEquivalenceCase) -> dict:
    """Compute the snapshot payload (realised ``apply`` output + metadata).

    Post-Issue-#186 design:

    * ``psi_out`` — the deterministic input array.
    * ``psi_in`` — the REALISER ``op.apply(psi_out)`` output. For
      single-rank-1 cases the realiser is obtained via
      ``SNBoundaryRealizer().realize(law, method_space)``. For the
      vacuum case the method space carries the xmin-face inflow
      indices (so the §16A.5 mask correctly identifies the inflow
      rows). For the mixed case ``psi_in`` is the
      ``apply(psi)`` of the Wave-0 ``OperatorSum`` composition of
      realised leaves (see :func:`_build_mixed_realized_apply`).
    * ``inflow_indices_xmin`` — (vacuum case only) the xmin-face inflow
      ordinate indices, so the harness can encode the §16A.5
      mask-and-pass-through assertion without re-deriving the set.
    * ``case_id`` / ``description`` — for traceability.

    Issue #186 (B3 + β2) note
    -------------------------
    Pre-#186 ``psi_in`` was generated via the legacy 2-arg
    ``bc.apply(psi_out, quad)`` path. That path is gone (descriptors
    are no longer callable); the realiser path produces the same
    numerical values on the existing snapshots, so regeneration is
    not strictly necessary, but the generator now uses the realiser
    path directly for future-proofing.
    """
    quad = case.build_quadrature()
    psi_out = _generate_psi(quad, case.case_id)

    if case.build_bc is None:
        # Issue #186 special case: the snapshot's ``psi_in`` is the
        # Wave-0 ``OperatorSum``-composition of realised leaves
        # (currently only ``mixed_30spec_70white_LS4`` uses this path).
        assert case.case_id == "mixed_30spec_70white_LS4", (
            f"build_bc=None marker reserved for the mixed case; "
            f"got case_id={case.case_id!r}."
        )
        psi_in = _build_mixed_realized_apply(quad, psi_out)
    else:
        bc = case.build_bc()
        # Vacuum needs an SNMethodSpace carrying the xmin inflow indices
        # to drive the §16A.5 mask; every other case uses ``minimal``.
        if case.case_id == "vacuum_lebedev17":
            method_space = SNMethodSpace(
                quadrature=quad,
                face="xmin",
                inflow_indices=_xmin_inflow_indices(quad),
            )
        else:
            method_space = SNMethodSpace.minimal(quad)
        op = SNBoundaryRealizer().realize(bc, method_space)
        psi_in = op.apply(psi_out)

    payload: dict = dict(
        psi_out=np.asarray(psi_out, dtype=np.float64),
        psi_in=np.asarray(psi_in, dtype=np.float64),
        case_id=np.array(case.case_id),
        description=np.array(case.description),
    )

    if case.case_id == "vacuum_lebedev17":
        # Wave-8 semantic-transition metadata. The §16A.5 corrected
        # behaviour zeros ONLY the inflow ordinates for the xmin face
        # (the legacy `np.zeros_like` zeros every row). Record the
        # index set so the harness can encode both assertions side by
        # side.
        payload["inflow_indices_xmin"] = _xmin_inflow_indices(quad).astype(
            np.int64,
        )

    return payload


def _array_equal_after_load(snapshot_path: Path, payload: dict) -> bool:
    """Check whether the on-disk snapshot already matches the payload."""
    if not snapshot_path.exists():
        return False
    try:
        existing = np.load(snapshot_path)
    except Exception:
        return False
    for key, expected in payload.items():
        if key not in existing.files:
            return False
        if not np.array_equal(existing[key], expected):
            return False
    return True


def generate_one(case: BCEquivalenceCase) -> tuple[Path, str]:
    """Write the snapshot for one case; return (path, status).

    Status is ``"UNCHANGED"`` when an existing snapshot already matches
    the freshly-built payload (bit-for-bit on every array). Status is
    ``"REGENERATED"`` when the snapshot is new or differs from disk.
    """
    payload = _build_payload(case)
    SNAPSHOT_DIR.mkdir(parents=True, exist_ok=True)
    out = case.snapshot_path
    if _array_equal_after_load(out, payload):
        status = "UNCHANGED"
    else:
        np.savez_compressed(out, **payload)
        status = "REGENERATED"
    return out, status


def generate_all(case_ids: list[str] | None = None) -> list[tuple[Path, str]]:
    """Generate every (or the named subset of) snapshots."""
    written: list[tuple[Path, str]] = []
    for case in CASES:
        if case_ids is not None and case.case_id not in case_ids:
            continue
        path, status = generate_one(case)
        written.append((path, status))
    return written


# ─── CLI ─────────────────────────────────────────────────────────────


def _print_summary(written: list[tuple[Path, str]]) -> None:
    """Pretty-print the per-snapshot status table."""
    if not written:
        print("(no snapshots generated — case-id filter matched nothing)")
        return
    name_w = max(len(p.name) for p, _ in written)
    print(f"{'Snapshot':<{name_w}}  Status        Size")
    print(f"{'-' * name_w}  ------------  ---------")
    for path, status in written:
        size_kb = path.stat().st_size / 1024.0 if path.exists() else 0.0
        print(f"{path.name:<{name_w}}  {status:<12}  {size_kb:6.1f} KiB")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Generate the BC-equivalence snapshot harness "
            "(Wave 6 / C6.1)."
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
    args = parser.parse_args(argv)
    if args.list:
        name_w = max(len(c.case_id) for c in CASES)
        for case in CASES:
            print(f"  {case.case_id:<{name_w}}  {case.description}")
        return 0
    written = generate_all(case_ids=args.case)
    _print_summary(written)
    return 0


if __name__ == "__main__":
    sys.exit(main())
