"""Diagnostic: randomized QMC (Owen-scrambled Sobol') angular quadrature for F.4.

Created by numerics-investigator on 2026-04-22 (Frame 5, cross-domain attack).

Purpose
-------
L17/L19/L20 (research log, 2026-04-22 late) established that F.4's own
signed error under the product-Gauss panel refinement (RICH 4 panels vs
RICH+panels 5 panels) flips sign at the anchor (sigma_t*R=10, rho=0.3)
and pathology point (sigma_t*R=20, rho=0.5). The integrand exp(-tau*d)
has bounded Hardy-Krause variation, so Owen-scrambled Sobol' should
give O(N^-1 (log N)^d) convergence with a tau-BOUNDED constant,
whereas product-Gauss bias at tau scales as tau^{2p} in the constant.

Option chosen: **A** (angular-only QMC).

Why A over B:
  - Option A: replace only the angular GL (n_ang knob, driving
    compute_P_esc_outer, compute_P_esc_inner, compute_G_bc_outer,
    compute_G_bc_inner) with Owen-scrambled Sobol' on [0, 1] mapped
    to the native angular range (sphere: [0, pi]). K_vol (which uses
    gl_nodes_weights directly, NOT gl_float) is left untouched, as is
    the radial composite-GL from composite_gl_r.
  - Option B (full 3D Sobol' over (r, r', mu)): requires unstacking
    the entire F.4 assembly — ~200 LOC of infrastructure. L17's
    root pathology is radial-panel-driven, so Option A may not fully
    kill L17 but IS the correct isolation of the angular contribution.
    If Option A passes at anchor AND pathology, Option B is
    unnecessary.

Implementation: monkey-patch peierls_geometry.gl_float to return Sobol'
nodes+weights on demand. Radial composite-GL uses gl_nodes_weights
(a different function) and is NOT affected.

Load-bearing tests (spec from frame attack memo):
  - Anchor: sigma_t*R=10, rho=0.3. RICH signed = +0.00329%.
  - Pathology: sigma_t*R=20, rho=0.5. RICH signed = +0.00543%,
    RICH+panels signed = -0.00394% (sign flip).

Pass criterion:
  - Anchor: 95% CI width < 0.003% AND sign stable across 32 scrambles.
  - Pathology: 95% CI of the bootstrap mean does NOT cross zero.

Usage
-----
    python scratch/derivations/diagnostics/diag_f4_qmc_quadrature.py

Subprocess isolation: each (point, quadrature-or-scramble) triple runs
in its own python subprocess with subprocess.run(..., timeout=600).
Lesson from L20 (Direction N): unbounded Brent fallbacks can wedge the
driver for hours if a single point hangs.

NOT a pytest test — this is a scanner. The conclusions feed back into
the L19 protocol proposal in the final memo.
"""
from __future__ import annotations

import json
import subprocess
import sys
import textwrap
import time
from pathlib import Path

# -------------------------------------------------------------------------
# Grid and QMC parameters
# -------------------------------------------------------------------------

ANCHOR = (10.0, 0.3)
PATHOLOGY = (20.0, 0.5)

FULL_GRID = [
    (5.0, 0.3),
    (5.0, 0.5),
    (10.0, 0.3),
    (10.0, 0.5),
    (20.0, 0.3),
    (20.0, 0.5),
]

QUAD_REFS = [
    ("RICH",         dict(n_panels=4, p_order=8,  n_ang=64)),
    ("RICH+panels",  dict(n_panels=5, p_order=8,  n_ang=64)),
]

# QMC sample count for the angular integral.
# Each scramble uses N_QMC Sobol' points in the angular variable.
# Radial composite-GL is (n_panels=4, p_order=8), fixed.
N_QMC = 4096
N_SCRAMBLES = 32
QMC_N_PANELS = 4  # match RICH radial
QMC_P_ORDER = 8   # match RICH radial

K_INF = 1.5  # nu_sig_f / (sig_t - sig_s) = 1.0 / (1.0 - 1/3)
WALL_BUDGET_S = 600

# -------------------------------------------------------------------------
# QMC worker (runs F.4 with gl_float monkey-patched to scrambled Sobol')
# -------------------------------------------------------------------------

_QMC_WORKER_CODE = textwrap.dedent(
    r"""
    # Batched QMC worker: build K_vol ONCE, then run n_scrambles closure
    # solves. Amortizes K_vol (dominant ~45s cost) across all scrambles.
    import json, sys, time
    import numpy as np
    from scipy.stats import qmc

    sys.path.insert(0, '/workspaces/ORPHEUS/scratch/derivations/diagnostics')
    from orpheus.derivations import peierls_geometry as pg
    from orpheus.derivations.peierls_geometry import (
        CurvilinearGeometry, composite_gl_r, build_volume_kernel,
        build_closure_operator,
    )
    from diag_cin_aware_split_basis_keff import solve_k_eff

    tau         = float(sys.argv[1])
    rho         = float(sys.argv[2])
    n_qmc       = int(sys.argv[3])
    seed_lo     = int(sys.argv[4])
    seed_hi     = int(sys.argv[5])  # exclusive
    n_pan       = int(sys.argv[6])
    p_ord       = int(sys.argv[7])

    K_INF       = 1.5

    if (n_qmc & (n_qmc - 1)) != 0:
        raise ValueError(f'n_qmc={n_qmc} must be a power of 2 for Sobol balance')
    m_log2 = int(round(np.log2(n_qmc)))

    # --- Build geometry + K_vol ONCE (unaffected by gl_float monkey-patch) ---
    R = tau
    r_0 = rho * R
    geom = CurvilinearGeometry(kind='sphere-1d', inner_radius=r_0)
    radii = np.array([R])
    sig_t = np.array([1.0])
    r_nodes, r_wts, panels = composite_gl_r(
        radii, n_pan, p_ord, dps=15, inner_radius=r_0,
    )
    t_kvol = time.time()
    K_vol = build_volume_kernel(
        geom, r_nodes, panels, radii, sig_t,
        n_angular=64, n_rho=64, dps=15,
    )
    kvol_wall = time.time() - t_kvol

    # --- Monkey-patch gl_float; it's rebound per scramble ---
    _orig_gl_float = pg.gl_float
    results = []
    total_closure_wall = 0.0

    try:
        for seed in range(seed_lo, seed_hi):
            sobol = qmc.Sobol(d=1, scramble=True, seed=seed)
            sobol_u = sobol.random_base2(m=m_log2).ravel()

            def _gl_float_qmc(n, a, b, dps=30, _u=sobol_u, _N=n_qmc):
                nodes = a + (b - a) * _u
                wts = np.full(_N, (b - a) / _N)
                return nodes, wts

            pg.gl_float = _gl_float_qmc

            t0 = time.time()
            op = build_closure_operator(
                geom, r_nodes, r_wts, radii, sig_t,
                reflection='white', n_angular=n_qmc, n_surf_quad=n_qmc, dps=15,
            )
            K = K_vol + op.as_matrix()
            k = solve_k_eff(K, 1.0, 1.0/3.0, 1.0)
            dt = time.time() - t0
            total_closure_wall += dt
            signed = (k - K_INF) / K_INF
            results.append({'seed': seed, 'k': k, 'signed': signed, 'wall': dt})
    finally:
        pg.gl_float = _orig_gl_float

    print(json.dumps({
        'results': results,
        'kvol_wall': kvol_wall,
        'total_closure_wall': total_closure_wall,
        'mean_closure_wall': total_closure_wall / max(1, len(results)),
    }))
    """
).strip()

# -------------------------------------------------------------------------
# Product-Gauss reference worker (reproduce RICH, RICH+panels).
# -------------------------------------------------------------------------

_PG_WORKER_CODE = textwrap.dedent(
    r"""
    import json, sys, time
    sys.path.insert(0, '/workspaces/ORPHEUS/scratch/derivations/diagnostics')
    from diag_cin_aware_split_basis_keff import run_scalar_f4

    tau      = float(sys.argv[1])
    rho      = float(sys.argv[2])
    n_panels = int(sys.argv[3])
    p_order  = int(sys.argv[4])
    n_ang    = int(sys.argv[5])
    K_INF    = 1.5

    R = tau
    r_0 = rho * R
    t0 = time.time()
    k = run_scalar_f4(r_0, R, 1.0, 1.0/3.0, 1.0,
                      n_panels=n_panels, p_order=p_order, n_ang=n_ang)
    dt = time.time() - t0
    signed = (k - K_INF) / K_INF
    print(json.dumps({'k': k, 'signed': signed, 'wall': dt}))
    """
).strip()


def run_pg(tau, rho, quad, timeout=WALL_BUDGET_S):
    args = [sys.executable, "-c", _PG_WORKER_CODE,
            str(tau), str(rho),
            str(quad["n_panels"]), str(quad["p_order"]), str(quad["n_ang"])]
    t0 = time.time()
    try:
        proc = subprocess.run(args, capture_output=True, text=True,
                              timeout=timeout, check=False)
    except subprocess.TimeoutExpired:
        return {"status": "timeout", "signed": None, "wall": timeout}
    if proc.returncode != 0:
        return {"status": "error", "signed": None, "wall": time.time()-t0,
                "msg": (proc.stderr or proc.stdout or "")[-300:]}
    try:
        out = json.loads(proc.stdout.strip().splitlines()[-1])
    except Exception as e:
        return {"status": "error", "signed": None, "wall": time.time()-t0,
                "msg": f"parse: {e}; stdout={proc.stdout!r}"}
    return {"status": "ok", **out}


def run_qmc_batch(tau, rho, n_qmc, seed_lo, seed_hi, timeout=WALL_BUDGET_S):
    """Run a batch of scrambles in one subprocess (amortizes K_vol build)."""
    args = [sys.executable, "-c", _QMC_WORKER_CODE,
            str(tau), str(rho), str(n_qmc),
            str(seed_lo), str(seed_hi),
            str(QMC_N_PANELS), str(QMC_P_ORDER)]
    t0 = time.time()
    try:
        proc = subprocess.run(args, capture_output=True, text=True,
                              timeout=timeout, check=False)
    except subprocess.TimeoutExpired:
        return {"status": "timeout", "results": [], "wall": timeout}
    if proc.returncode != 0:
        return {"status": "error", "results": [], "wall": time.time()-t0,
                "msg": (proc.stderr or proc.stdout or "")[-500:]}
    try:
        out = json.loads(proc.stdout.strip().splitlines()[-1])
    except Exception as e:
        return {"status": "error", "results": [], "wall": time.time()-t0,
                "msg": f"parse: {e}; stdout={proc.stdout!r}"}
    return {"status": "ok", **out, "subprocess_wall": time.time()-t0}


# -------------------------------------------------------------------------
# Bootstrap 95% CI for the mean signed error
# -------------------------------------------------------------------------


def bootstrap_mean_ci(samples, n_boot=10_000, alpha=0.05, rng_seed=0):
    import numpy as np
    samples = np.asarray(samples, dtype=float)
    n = len(samples)
    rng = np.random.default_rng(rng_seed)
    boot_means = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, n, size=n)
        boot_means[b] = samples[idx].mean()
    lo = float(np.percentile(boot_means, 100 * alpha/2))
    hi = float(np.percentile(boot_means, 100 * (1 - alpha/2)))
    return lo, hi, float(samples.mean()), float(samples.std(ddof=1))


# -------------------------------------------------------------------------
# Scan drivers
# -------------------------------------------------------------------------


def run_point(tau, rho, *, n_scrambles=N_SCRAMBLES, n_qmc=N_QMC,
              quad_refs=QUAD_REFS, verbose=True):
    out = {"tau": tau, "rho": rho, "pg": {}, "qmc": {}}
    if verbose:
        print(f"\n=== Point (sigma_t*R = {tau}, rho = {rho}) ===")

    # Product-Gauss references
    for label, quad in quad_refs:
        t0 = time.time()
        res = run_pg(tau, rho, quad)
        if verbose:
            if res["status"] == "ok":
                print(f"  PG {label:<12}: signed = {res['signed']*100:+.6f}%  "
                      f"wall = {res['wall']:.1f}s")
            else:
                print(f"  PG {label:<12}: [{res['status']}] "
                      f"{res.get('msg', '')[:80]}")
        out["pg"][label] = {"quad": quad, **res}

    # QMC scrambles — batched in one subprocess to amortize K_vol.
    if verbose:
        print(f"  QMC N={n_qmc}, {n_scrambles} scrambles (batched)...")
    qmc_t0 = time.time()
    batch = run_qmc_batch(tau, rho, n_qmc, 0, n_scrambles, timeout=WALL_BUDGET_S)
    qmc_total_wall = time.time() - qmc_t0
    qmc_signed = []
    qmc_walls = []
    if batch["status"] == "ok":
        for r in batch["results"]:
            qmc_signed.append(r["signed"])
            qmc_walls.append(r["wall"])
        if verbose:
            print(f"    K_vol build: {batch.get('kvol_wall', 0):.1f}s  "
                  f"closure mean: {batch.get('mean_closure_wall', 0):.2f}s  "
                  f"total: {qmc_total_wall:.1f}s")
    else:
        if verbose:
            print(f"    batch {batch['status']}: "
                  f"{batch.get('msg', '')[:300]}")

    if len(qmc_signed) >= 2:
        import numpy as np
        samples = np.array(qmc_signed)
        lo, hi, mean, std = bootstrap_mean_ci(samples)
        ci_width = (hi - lo)
        sign_stable = all(s > 0 for s in samples) or all(s < 0 for s in samples)
        sign_of_mean = "+" if mean > 0 else "-"
        ci_crosses_zero = (lo < 0.0 < hi)
        out["qmc"] = {
            "n_scrambles_ok": len(qmc_signed),
            "samples": qmc_signed,
            "mean": mean, "std": std,
            "ci_lo": lo, "ci_hi": hi, "ci_width": ci_width,
            "sign_stable": sign_stable,
            "ci_crosses_zero": ci_crosses_zero,
            "total_wall": qmc_total_wall,
            "mean_wall_per_scramble": (sum(qmc_walls) / len(qmc_walls)),
        }
        if verbose:
            print(f"  QMC bootstrap: mean = {mean*100:+.6f}% "
                  f"std = {std*100:.6f}%")
            print(f"  95% CI = [{lo*100:+.6f}%, {hi*100:+.6f}%]  "
                  f"width = {ci_width*100:.6f}%")
            print(f"  sign stable across scrambles: {sign_stable}")
            print(f"  CI crosses zero: {ci_crosses_zero}")
            print(f"  Total QMC wall: {qmc_total_wall:.1f}s "
                  f"(mean per scramble: {out['qmc']['mean_wall_per_scramble']:.1f}s)")
    else:
        out["qmc"] = {"n_scrambles_ok": len(qmc_signed),
                      "total_wall": qmc_total_wall, "samples": qmc_signed}
        if verbose:
            print(f"  QMC: only {len(qmc_signed)} scrambles succeeded.")
    return out


def main():
    print("=" * 78)
    print("Frame 5 — randomized Sobol' angular quadrature for F.4")
    print("=" * 78)
    print(f"N_QMC = {N_QMC}  (power of 2: {N_QMC & (N_QMC-1) == 0})")
    print(f"N_SCRAMBLES = {N_SCRAMBLES}")
    print(f"Radial composite-GL: n_panels={QMC_N_PANELS}, p_order={QMC_P_ORDER}")
    print(f"K_vol angular GL: n_angular=64 (fixed; uses gl_nodes_weights "
          f"directly, NOT monkey-patched)")
    print()

    grand_t0 = time.time()
    all_results = {}

    # Determine points to scan from env / CLI.
    # Default: anchor + pathology (fast).
    import os
    mode = os.environ.get("QMC_SCAN_MODE", "anchor").lower()
    if mode == "full":
        points = FULL_GRID
        print(f"Scan mode: FULL 6-point grid.")
    else:
        points = [ANCHOR, PATHOLOGY]
        print(f"Scan mode: anchor + pathology only (set QMC_SCAN_MODE=full for 6-point).")

    for pt in points:
        all_results[pt] = run_point(*pt)

    grand_wall = time.time() - grand_t0
    print(f"\nTotal scan wall: {grand_wall:.1f}s")

    # Dump JSON
    # Coerce keys to strings for JSON.
    def _coerce(v):
        if isinstance(v, float):
            return v
        return v
    dump = {}
    for (tau, rho), rec in all_results.items():
        dump[f"{tau}_{rho}"] = rec
    out_path = Path("/tmp/diag_f4_qmc_quadrature.json")
    out_path.write_text(json.dumps(dump, indent=2, default=float))
    print(f"Machine-readable dump: {out_path}")

    # Final summary table
    print("\n" + "=" * 78)
    print("SUMMARY — anchor & pathology")
    print("=" * 78)
    print(f"{'point':<18} {'PG RICH':<18} {'PG RICH+panels':<22} "
          f"{'QMC mean':<18} {'QMC 95% CI':<32} {'verdict':<10}")
    for (tau, rho), rec in all_results.items():
        pg_rich = rec["pg"].get("RICH", {}).get("signed")
        pg_rp = rec["pg"].get("RICH+panels", {}).get("signed")
        q = rec.get("qmc", {})
        mean = q.get("mean")
        lo, hi = q.get("ci_lo"), q.get("ci_hi")
        crosses = q.get("ci_crosses_zero")
        if mean is None:
            verdict = "NO QMC"
        elif crosses:
            verdict = "CI0"
        elif q.get("sign_stable", False):
            verdict = "PASS"
        else:
            verdict = "LOOSE"
        print(f"({tau},{rho})".ljust(18)
              + (f"{pg_rich*100:+.6f}%" if pg_rich is not None else "err").ljust(18)
              + (f"{pg_rp*100:+.6f}%" if pg_rp is not None else "err").ljust(22)
              + (f"{mean*100:+.6f}%" if mean is not None else "err").ljust(18)
              + (f"[{lo*100:+.5f}%, {hi*100:+.5f}%]"
                 if lo is not None else "err").ljust(32)
              + verdict)

    return all_results


if __name__ == "__main__":
    main()
