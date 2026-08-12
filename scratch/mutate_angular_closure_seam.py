"""In-process mutation plugin for the SN curvilinear angular-closure seam.

Usage:
    MUT=<name> PYTHONPATH=<this dir> .venv/bin/python -O -m pytest \
        -p mutate_seam <targets>

Discipline (test-architect lessons L44i / L46e / vv #17):
  * mutants are built by TRANSFORMING ``inspect.getsource``, never by
    hand-copying a production body (a hand copy is a twin path);
  * a ``str.replace`` whose target is ABSENT raises, so the instrument
    asserts its own installation instead of printing a banner;
  * every ``sys.modules`` entry whose attribute ``is`` the original is
    rebound, and the rebind COUNT is asserted > 0.
"""
from __future__ import annotations

import inspect
import os
import sys
import textwrap

import numpy as np

import orpheus.sn.sweep.pole_angular_closure as P
from orpheus.geometry import CoordSystem

_TARGET = os.environ.get("MUT", "")


def _recompile(fn, pairs):
    """Return a mutated copy of ``fn`` built from its own source."""
    src = textwrap.dedent(inspect.getsource(fn))
    for old, new in pairs:
        if old not in src:
            raise RuntimeError(
                f"mutation target absent in {fn.__name__}: {old!r} "
                f"— the instrument is stale, not the code"
            )
        src = src.replace(old, new, 1)
    ns = dict(P.__dict__)
    exec(compile(src, f"<mutant:{fn.__name__}>", "exec"), ns)
    return ns[fn.__name__]


def _rebind(name, new_obj):
    original = getattr(P, name)
    n = 0
    for mod in list(sys.modules.values()):
        if mod is None:
            continue
        try:
            if getattr(mod, name, None) is original:
                setattr(mod, name, new_obj)
                n += 1
        except Exception:  # pragma: no cover - defensive
            continue
    if n == 0:
        raise RuntimeError(f"rebound 0 bindings of {name!r} — instrument dead")
    return n


# ── mutant definitions ────────────────────────────────────────────────
# Each entry: (symbol name, list of (old, new) source replacements)
# or a callable(name) -> object for a wholesale replacement.

_EDGES = "angular_cell_edges_per_level"
_TAU = "morel_montry_tau_per_level"
_REC = "_psi_half_grid_single_level"

_SOURCE_MUTANTS = {
    # ── POSITIVE CONTROL ────────────────────────────────────────────
    "MC_control_uniform_partition": (
        _EDGES,
        [("return (mu_edge,)", "return (np.linspace(-1.0, 1.0, N + 1),)"),
         ("edges.append(sin_theta * np.cos(edge_omega))",
          "edges.append(np.linspace(-sin_theta, sin_theta, M + 1))")],
    ),
    # ── sphere partition ────────────────────────────────────────────
    "M1_sphere_seed_zero": (_EDGES, [("mu_edge[0] = -1.0", "mu_edge[0] = 0.0")]),
    "M2_sphere_half_weight": (
        _EDGES, [("mu_edge[n + 1] = mu_edge[n] + w[n]",
                  "mu_edge[n + 1] = mu_edge[n] + 0.5 * w[n]")],
    ),
    "M3_sphere_no_accumulation": (
        _EDGES, [("mu_edge[n + 1] = mu_edge[n] + w[n]",
                  "mu_edge[n + 1] = -1.0 + (n + 1) * (2.0 / N)")],
    ),
    "M3b_sphere_palindromic_weights": (
        _EDGES, [("mu_edge[n + 1] = mu_edge[n] + w[n]",
                  "mu_edge[n + 1] = mu_edge[n] + w[N - 1 - n]")],
    ),
    # ── cylinder partition ──────────────────────────────────────────
    "M4_cyl_chord_partition": (
        _EDGES,
        [("edges.append(sin_theta * np.cos(edge_omega))",
          "edge_eta = np.empty(M + 1)\n"
          "            edge_eta[0] = -sin_theta\n"
          "            edge_eta[M] = sin_theta\n"
          "            edge_eta[1:M] = 0.5 * (eta[:-1] + eta[1:])\n"
          "            edges.append(edge_eta)")],
    ),
    "M5_cyl_drop_sin_theta": (
        _EDGES, [("edges.append(sin_theta * np.cos(edge_omega))",
                  "edges.append(np.cos(edge_omega))")],
    ),
    "M6_cyl_endpoint_swap": (
        _EDGES, [("edge_omega[0] = np.pi\n", "edge_omega[0] = 0.0\n")],
    ),
    "M7_cyl_left_node_edge": (
        _EDGES, [("edge_omega[1:M] = 0.5 * (omega[:-1] + omega[1:])",
                  "edge_omega[1:M] = ("
                  "0.75 * omega[:-1] + 0.25 * omega[1:])")],
    ),
    # ── tau chart ───────────────────────────────────────────────────
    "M8_tau_one_minus_tau": (
        _TAU, [("(mu[m] - mu_edge[m]) / dmu",
                "(mu_edge[m + 1] - mu[m]) / dmu")],
    ),
    "M9_tau_plain_diamond": (
        _TAU, [("(mu[m] - mu_edge[m]) / dmu if abs(dmu) > 1e-15 else 0.5",
                "0.5")],
    ),
    # ── psi-half recurrence ─────────────────────────────────────────
    "M10_recurrence_wrong_divisor": (
        _REC, [(") / tau_m", ") / (1.0 - tau_m)")],
    ),
    "M11_recurrence_drop_thread_weight": (
        _REC, [("(1.0 - tau_m) * psi_half[:, m, :]", "0.0 * psi_half[:, m, :]")],
    ),
}


def pytest_configure(config):
    if not _TARGET:
        return
    if _TARGET not in _SOURCE_MUTANTS:
        raise RuntimeError(f"unknown mutation {_TARGET!r}")
    name, pairs = _SOURCE_MUTANTS[_TARGET]
    mutant = _recompile(getattr(P, name), pairs)
    n = _rebind(name, mutant)
    print(f"\n[MUTATE] {_TARGET}: rebound {name!r} in {n} module(s)")


def pytest_report_header(config):
    return f"MUTATION = {_TARGET or '(none)'}"


# ── smoke-test entry point: prove each mutant's OUTPUT is the defect ──
if __name__ == "__main__":
    from orpheus.numerics.quadrature import Quadrature

    sph = Quadrature.gauss_legendre(4)
    cyl = Quadrature.folded_product(n_mu=4, n_phi=8)
    base_e_s = P.angular_cell_edges_per_level(sph, CoordSystem.SPHERICAL)[0]
    base_e_c = P.angular_cell_edges_per_level(cyl, CoordSystem.CYLINDRICAL)[0]
    base_t_s = P.morel_montry_tau_per_level(sph, CoordSystem.SPHERICAL)[0]
    base_t_c = P.morel_montry_tau_per_level(cyl, CoordSystem.CYLINDRICAL)[0]
    print("BASE sphere edges :", np.round(base_e_s, 6))
    print("BASE cyl    edges :", np.round(base_e_c, 6))
    print("BASE sphere tau   :", np.round(base_t_s, 6))
    print("BASE cyl    tau   :", np.round(base_t_c, 6))
    print()
    for label, (name, pairs) in _SOURCE_MUTANTS.items():
        mutant = _recompile(getattr(P, name), pairs)
        if name == _EDGES:
            try:
                s = mutant(sph, CoordSystem.SPHERICAL)[0]
                ds = float(np.max(np.abs(s - base_e_s)))
            except Exception as exc:
                ds = f"{type(exc).__name__}"
            try:
                c = mutant(cyl, CoordSystem.CYLINDRICAL)[0]
                dc = float(np.max(np.abs(c - base_e_c)))
            except Exception as exc:
                dc = f"{type(exc).__name__}"
            print(f"{label:36s} d(sphere edges)={ds}  d(cyl edges)={dc}")
        elif name == _TAU:
            try:
                s = mutant(sph, CoordSystem.SPHERICAL)[0]
                ds = float(np.max(np.abs(s - base_t_s)))
            except Exception as exc:
                ds = f"{type(exc).__name__}"
            try:
                c = mutant(cyl, CoordSystem.CYLINDRICAL)[0]
                dc = float(np.max(np.abs(c - base_t_c)))
            except Exception as exc:
                dc = f"{type(exc).__name__}"
            print(f"{label:36s} d(sphere tau)={ds}  d(cyl tau)={dc}")
        else:
            psi = np.arange(1.0, 1.0 + 2 * 4 * 3).reshape(2, 4, 3)
            base = P._psi_half_grid_single_level(psi, base_t_c)
            got = mutant(psi, base_t_c)
            print(f"{label:36s} d(psi_half)="
                  f"{float(np.max(np.abs(got - base))):.4g}")
