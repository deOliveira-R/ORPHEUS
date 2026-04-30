"""Diagnostic: F.4 structural-floor baseline across the 6-point reference grid.

Created by numerics-investigator on 2026-04-22 (Direction N, Issue #123).

Purpose
-------
L17-L19 (research log, 2026-04-22 late session) established empirically that
F.4's own RICH reference is quadrature-limited at sigma_t * R >= 10:
  * anchor (tau=10, rho=0.3): RICH err = 0.0033%, ULTRA err = 0.0010% (3x lower).
  * (tau=20, rho=0.5):       RICH err = +0.0054%, RICH+panels err = -0.0039%.
The second point is a SIGN FLIP of F.4's own signed error under one-panel
refinement -- F.4's true structural floor sits below RICH's quadrature noise.

This scanner pins an L19-compliant reference quadrature per point on the
canonical 6-point grid (sigma_t*R in {5, 10, 20}; rho = r_0/R in {0.3, 0.5}).
For each point it runs F.4 at RICH, RICH+panels, RICH+pp, and ULTRA (or as
far as a 2-min/point budget allows), records (quadrature, signed err, wall),
and prints a table that is exactly the input for the L19 protocol test's
reference baseline.

A point is declared RESOLVED when there exists a quadrature Q_i in the
scanned ladder such that:
  (a) |signed_err(Q_i)| < 0.005%,
  (b) |signed_err(Q_i)| < |signed_err(Q_{i-1})|      (monotone DOWN), and
  (c) sign(signed_err(Q_i)) == sign(signed_err(Q_{i-1}))  (no late flip).
The earliest such Q_i is the per-point REFERENCE quadrature.

If no quadrature in the ladder satisfies (a-c) within 120 s, the point is
flagged "unresolved at 2-min budget". The signed-err trajectory is still
printed so the extrapolation cost can be estimated.

Process isolation
-----------------
Each (point, quad) invocation runs in its OWN python subprocess with
subprocess.run(..., timeout=120). A single runaway must not wedge the whole
scan. This is the lesson from Direction C's first kill (3 hr hang, rescued
by per-point subprocess timeout). F.4 is much lighter than that workload,
but the pattern is the same.

Usage
-----
    python scratch/derivations/diagnostics/diag_f4_structural_floor_baseline.py

Output: a plain-text Markdown-friendly table at the end.

NOT a pytest test -- this is a scanner. The distilled baseline is copied
into tests/cp/test_peierls_rank_n_protocol.py as the F.4 per-point reference.
"""
from __future__ import annotations

import json
import subprocess
import sys
import textwrap
import time
from pathlib import Path

# -------------------------------------------------------------------------
# Grid and quadrature ladder
# -------------------------------------------------------------------------

GRID = [
    (5.0, 0.3),
    (5.0, 0.5),
    (10.0, 0.3),
    (10.0, 0.5),
    (20.0, 0.3),
    (20.0, 0.5),
]

# Ladder from Direction C (diag_pca_sectors_rich_adaptive.py).
# RICH is the historical reference; RICH+panels revealed L17's sign flip
# at (20, 0.5). RICH+pp and ULTRA further isolate the structural floor.
QUAD_LADDER = [
    ("RICH",         dict(n_panels=4, p_order=8,  n_ang=64)),
    ("RICH+panels",  dict(n_panels=5, p_order=8,  n_ang=64)),
    ("RICH+pp",      dict(n_panels=5, p_order=10, n_ang=64)),
    ("ULTRA",        dict(n_panels=5, p_order=10, n_ang=96)),
]

K_INF = 1.5  # nu_sig_f / (sig_t - sig_s) = 1.0 / (1.0 - 1/3) = 1.5
WALL_BUDGET_S = 120  # per-point/per-quad subprocess timeout

RESOLUTION_THRESHOLD = 5.0e-5  # 0.005%


# -------------------------------------------------------------------------
# Subprocess worker -- runs a single (tau, rho, quad) triple, prints JSON.
# -------------------------------------------------------------------------

_WORKER_CODE = textwrap.dedent(
    r"""
    import json, sys, time
    sys.path.insert(0, '/workspaces/ORPHEUS/scratch/derivations/diagnostics')
    from diag_cin_aware_split_basis_keff import run_scalar_f4

    tau = float(sys.argv[1])
    rho = float(sys.argv[2])
    n_panels = int(sys.argv[3])
    p_order  = int(sys.argv[4])
    n_ang    = int(sys.argv[5])
    K_INF    = 1.5

    R = tau  # sigma_t = 1 so R = sigma_t * R
    r_0 = rho * R

    t0 = time.time()
    k = run_scalar_f4(r_0, R, 1.0, 1.0/3.0, 1.0,
                      n_panels=n_panels, p_order=p_order, n_ang=n_ang)
    dt = time.time() - t0
    signed = (k - K_INF) / K_INF  # SIGNED fractional residual
    print(json.dumps({"k": k, "signed": signed, "wall": dt}))
    """
).strip()


def run_point_quad(tau: float, rho: float, quad: dict, wall_budget: float = WALL_BUDGET_S):
    """Run F.4 at (tau, rho) with the given quadrature in a subprocess.

    Returns
    -------
    dict with keys:
      status : "ok" | "timeout" | "error"
      signed : signed fractional residual (k - K_INF) / K_INF
               (None on timeout/error)
      wall   : wall-time seconds (None on timeout/error)
      msg    : short message
    """
    args = [
        sys.executable, "-c", _WORKER_CODE,
        str(tau), str(rho),
        str(quad["n_panels"]), str(quad["p_order"]), str(quad["n_ang"]),
    ]
    t0 = time.time()
    try:
        proc = subprocess.run(
            args, capture_output=True, text=True, timeout=wall_budget, check=False,
        )
    except subprocess.TimeoutExpired:
        return {"status": "timeout", "signed": None, "wall": wall_budget,
                "msg": f"timed out at {wall_budget:.0f}s"}
    total_wall = time.time() - t0
    if proc.returncode != 0:
        return {"status": "error", "signed": None, "wall": total_wall,
                "msg": (proc.stderr or proc.stdout or "")[-200:]}
    try:
        out = json.loads(proc.stdout.strip().splitlines()[-1])
    except (json.JSONDecodeError, IndexError) as exc:
        return {"status": "error", "signed": None, "wall": total_wall,
                "msg": f"parse: {exc}; stdout={proc.stdout!r}"}
    return {"status": "ok", "signed": out["signed"], "wall": out["wall"],
            "msg": ""}


# -------------------------------------------------------------------------
# Reference-quadrature selection per point
# -------------------------------------------------------------------------


def select_reference_quadrature(results: list[dict]):
    """Pick the earliest quadrature satisfying the L19 resolution criteria.

    results: list of dicts in ladder order with keys
             {'label', 'signed', 'status', 'wall'}

    Returns the dict of the chosen reference, or None if unresolved.
    """
    # Walk forward; need two consecutive OK points to evaluate monotonicity
    # and sign stability.
    ok_points = [(i, r) for i, r in enumerate(results) if r["status"] == "ok"]
    if len(ok_points) < 2:
        return None
    for idx in range(1, len(ok_points)):
        i_prev, prev = ok_points[idx - 1]
        i_curr, curr = ok_points[idx]
        # abs smaller than threshold
        if abs(curr["signed"]) >= RESOLUTION_THRESHOLD:
            continue
        # monotone DOWN in magnitude
        if abs(curr["signed"]) >= abs(prev["signed"]):
            continue
        # sign-stable between the two consecutive quads
        if (prev["signed"] > 0) != (curr["signed"] > 0) and prev["signed"] != 0:
            continue
        return curr
    return None


# -------------------------------------------------------------------------
# Scanner
# -------------------------------------------------------------------------


def scan():
    print("=" * 78)
    print("F.4 structural-floor baseline -- Direction N / Issue #123")
    print("=" * 78)
    print(f"Grid: {GRID}")
    print(f"Ladder: {[label for label, _ in QUAD_LADDER]}")
    print(f"Per-run wall budget: {WALL_BUDGET_S}s")
    print(f"Resolution threshold: |signed err| < {RESOLUTION_THRESHOLD*100:.3f}%"
          f" AND monotone DOWN AND sign-stable")
    print()

    grand_t0 = time.time()
    all_results = {}

    for (tau, rho) in GRID:
        print(f"--- Point (sigma_t*R = {tau}, rho = {rho}) ---")
        point_results = []
        for label, quad in QUAD_LADDER:
            res = run_point_quad(tau, rho, quad)
            res["label"] = label
            res["quad"] = quad
            point_results.append(res)
            if res["status"] == "ok":
                print(f"  {label:<12} {quad}: "
                      f"signed = {res['signed']*100:+.6f}%  "
                      f"wall = {res['wall']:.1f}s")
            else:
                print(f"  {label:<12} {quad}: [{res['status']}] {res['msg']}")
        ref = select_reference_quadrature(point_results)
        all_results[(tau, rho)] = {"ladder": point_results, "reference": ref}
        if ref:
            print(f"  REFERENCE: {ref['label']}  signed = {ref['signed']*100:+.6f}%")
        else:
            print(f"  UNRESOLVED at 2-min budget.")
        print()

    grand_wall = time.time() - grand_t0
    print(f"Total scan wall: {grand_wall:.1f}s")
    print()

    # -------------------------------------------------------------------------
    # Final table (plain-text)
    # -------------------------------------------------------------------------
    print("=" * 78)
    print("SUMMARY TABLE")
    print("=" * 78)
    header = (
        f"| {'sigma_t*R':<9} | {'rho':<4} | {'RICH signed':<13} | "
        f"{'R+pan signed':<13} | {'R+pp signed':<13} | {'ULTRA signed':<13} | "
        f"{'reference':<12} |"
    )
    print(header)
    print("|" + "-" * (len(header) - 2) + "|")
    for (tau, rho), rec in all_results.items():
        ladder = {r["label"]: r for r in rec["ladder"]}
        cells = [f"{tau:<9}", f"{rho:<4}"]
        for lbl in ["RICH", "RICH+panels", "RICH+pp", "ULTRA"]:
            r = ladder[lbl]
            if r["status"] == "ok":
                cells.append(f"{r['signed']*100:+.6f}%".ljust(13))
            elif r["status"] == "timeout":
                cells.append("timeout".ljust(13))
            else:
                cells.append("error".ljust(13))
        ref = rec["reference"]
        cells.append((ref["label"] if ref else "UNRESOLVED").ljust(12))
        print("| " + " | ".join(cells) + " |")

    # -------------------------------------------------------------------------
    # Machine-readable dump
    # -------------------------------------------------------------------------
    dump = []
    for (tau, rho), rec in all_results.items():
        ladder = [
            {"label": r["label"], "quad": r["quad"],
             "status": r["status"], "signed": r["signed"], "wall": r["wall"]}
            for r in rec["ladder"]
        ]
        ref = rec["reference"]
        dump.append({
            "tau": tau, "rho": rho,
            "ladder": ladder,
            "reference": (
                {"label": ref["label"], "quad": ref["quad"],
                 "signed": ref["signed"], "wall": ref["wall"]}
                if ref else None
            ),
        })
    out_path = Path("/tmp/diag_f4_structural_floor_baseline.json")
    out_path.write_text(json.dumps(dump, indent=2))
    print(f"\nMachine-readable dump written to {out_path}")

    return all_results


if __name__ == "__main__":
    scan()
