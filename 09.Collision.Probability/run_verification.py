#!/usr/bin/env python3
"""Formal verification of transport solvers using synthetic benchmarks.

Runs 1-group, 2-group, and 4-group problems with known analytical
solutions against the homogeneous, CP, SN, diffusion, MOC, and MC solvers.
"""

import sys
_root = str(__import__('pathlib').Path(__file__).resolve().parent.parent)
sys.path.insert(0, _root)
sys.path.insert(0, str(__import__('pathlib').Path(_root) / '01.Homogeneous.Reactors'))
sys.path.insert(0, str(__import__('pathlib').Path(_root) / '02.Discrete.Ordinates'))
sys.path.insert(0, str(__import__('pathlib').Path(_root) / '03.Method.Of.Characteristics'))
sys.path.insert(0, str(__import__('pathlib').Path(_root) / '04.Monte.Carlo'))
sys.path.insert(0, str(__import__('pathlib').Path(_root) / '05.Diffusion.1D'))

import numpy as np
from data.macro_xs.benchmarks import (
    benchmark_1g_homogeneous,
    benchmark_2g_homogeneous,
    benchmark_4g_homogeneous,
    benchmark_1g_slab,
    benchmark_2g_slab,
    benchmark_1g_cylinder,
    benchmark_2g_cylinder,
    _kinf_from_cp,
    _make_mixture,
)
from homogeneous import solve_homogeneous_infinite
from collision_probability_slab import SlabGeometry, solve_slab_cp, _compute_slab_cp_group
from collision_probability import CPGeometry, solve_collision_probability
from discrete_ordinates import (
    DOParams, PinCellGeometry, Quadrature, solve_discrete_ordinates,
)


def run_homogeneous_benchmarks():
    """Verify the infinite-medium eigenvalue solver."""
    print("=" * 65)
    print("HOMOGENEOUS INFINITE MEDIUM BENCHMARKS")
    print("=" * 65)
    print(f"{'Benchmark':<22s}  {'Groups':>6s}  {'Analytical':>10s}  "
          f"{'Solver':>10s}  {'Error':>10s}  {'OK':>3s}")
    print("-" * 65)

    for label, bench_fn in [
        ("1G homogeneous", benchmark_1g_homogeneous),
        ("2G homogeneous", benchmark_2g_homogeneous),
        ("4G homogeneous", benchmark_4g_homogeneous),
    ]:
        mix, k_analytical = bench_fn()
        result = solve_homogeneous_infinite(mix)
        err = abs(result.k_inf - k_analytical)
        ok = "✓" if err < 1e-4 else "✗"
        print(f"{label:<22s}  {mix.ng:>6d}  {k_analytical:10.6f}  "
              f"{result.k_inf:10.6f}  {err:10.2e}  {ok:>3s}")

    print()


def run_slab_benchmarks():
    """Verify the slab CP solver against analytical CP eigenvalues."""
    print("=" * 65)
    print("HETEROGENEOUS SLAB CP BENCHMARKS")
    print("=" * 65)
    print(f"{'Benchmark':<22s}  {'Groups':>6s}  {'Analytical':>10s}  "
          f"{'Solver':>10s}  {'Error':>10s}  {'OK':>3s}")
    print("-" * 65)

    for label, bench_fn in [
        ("1G 2-region slab", benchmark_1g_slab),
        ("2G 2-region slab", benchmark_2g_slab),
    ]:
        materials, geom_params, k_analytical = bench_fn()
        ng = materials[2].ng

        geom = SlabGeometry.default_pwr(**geom_params)
        result = solve_slab_cp(materials, geom,
                               keff_tol=1e-7, flux_tol=1e-6)

        err = abs(result.keff - k_analytical)
        ok = "✓" if err < 1e-3 else "✗"
        print(f"{label:<22s}  {ng:>6d}  {k_analytical:10.6f}  "
              f"{result.keff:10.6f}  {err:10.2e}  {ok:>3s}")

    print()


def run_cylinder_benchmarks():
    """Verify the cylindrical (Wigner-Seitz) CP solver."""
    print("=" * 65)
    print("HETEROGENEOUS CYLINDRICAL CP BENCHMARKS")
    print("=" * 65)
    print(f"{'Benchmark':<22s}  {'Groups':>6s}  {'Analytical':>10s}  "
          f"{'Solver':>10s}  {'Error':>10s}  {'OK':>3s}")
    print("-" * 65)

    for label, bench_fn in [
        ("1G 2-region cyl", benchmark_1g_cylinder),
        ("2G 2-region cyl", benchmark_2g_cylinder),
    ]:
        materials, geom_params, k_analytical = bench_fn()
        ng = materials[2].ng

        geom = CPGeometry.default_pwr(**geom_params)
        result = solve_collision_probability(
            materials, geom,
        )

        err = abs(result.keff - k_analytical)
        ok = "✓" if err < 1e-3 else "✗"
        print(f"{label:<22s}  {ng:>6d}  {k_analytical:10.6f}  "
              f"{result.keff:10.6f}  {err:10.2e}  {ok:>3s}")

    print()


def _slab_geom_for_do(n_fuel, n_mod, delta, ny=2):
    """Build an SN PinCellGeometry for a 2-region slab benchmark.

    Returns the geometry AND the effective material thicknesses (accounting
    for boundary half-volumes) so that the analytical CP reference can be
    computed for the exact same dimensions.
    """
    nx = n_fuel + n_mod
    vol = np.full((nx, ny), delta**2)
    vol[0, :] /= 2
    vol[-1, :] /= 2
    vol[:, 0] /= 2
    vol[:, -1] /= 2

    mat = np.zeros((nx, ny), dtype=int)
    mat[:n_fuel, :] = 2   # fuel
    mat[n_fuel:, :] = 0   # moderator

    geom = PinCellGeometry(nx=nx, ny=ny, delta=delta, mat_map=mat, volume=vol)

    # Effective thicknesses (sum of x-widths per material)
    x_widths = np.full(nx, delta)
    x_widths[0] /= 2
    x_widths[-1] /= 2
    t_fuel_eff = x_widths[:n_fuel].sum()
    t_mod_eff = x_widths[n_fuel:].sum()

    return geom, t_fuel_eff, t_mod_eff


def _analytical_slab_kinf(fuel_mix, mod_mix, t_fuel, t_mod):
    """Compute the analytical CP eigenvalue for a 2-region slab."""
    ng = fuel_mix.ng
    sig_t_all = np.array([
        fuel_mix.SigT,
        mod_mix.SigT,
    ])  # (2, ng)

    # Build P_inf for each group using the slab E_3 formula
    from scipy.special import expn

    def e3(x):
        return float(expn(3, max(x, 0.0)))

    t_arr = np.array([t_fuel, t_mod])
    P_inf_g = np.zeros((2, 2, ng))

    for g in range(ng):
        sig_t_g = sig_t_all[:, g]
        tau = sig_t_g * t_arr
        bnd_pos = np.array([0.0, tau[0], tau[0] + tau[1]])

        rcp = np.zeros((2, 2))
        for i in range(2):
            sti, tau_i = sig_t_g[i], tau[i]
            rcp[i, i] += 0.5 * sti * (2 * t_arr[i] - (2.0 / sti) * (0.5 - e3(tau_i)))
            for j in range(2):
                tau_j = tau[j]
                if j > i:
                    gap_d = bnd_pos[j] - bnd_pos[i + 1]
                elif j < i:
                    gap_d = bnd_pos[i] - bnd_pos[j + 1]
                else:
                    gap_d = None
                if gap_d is not None:
                    gap_d = max(gap_d, 0.0)
                    dd = e3(gap_d) - e3(gap_d + tau_i) - e3(gap_d + tau_j) \
                         + e3(gap_d + tau_i + tau_j)
                else:
                    dd = 0.0
                gap_c = bnd_pos[i] + bnd_pos[j]
                dc = e3(gap_c) - e3(gap_c + tau_i) - e3(gap_c + tau_j) \
                     + e3(gap_c + tau_i + tau_j)
                rcp[i, j] += 0.5 * (dd + dc)

        P_cell = np.zeros((2, 2))
        for i in range(2):
            P_cell[i, :] = rcp[i, :] / (sig_t_g[i] * t_arr[i])
        P_out = np.maximum(1.0 - P_cell.sum(axis=1), 0.0)
        P_in = sig_t_g * t_arr * P_out
        P_inout = max(1.0 - P_in.sum(), 0.0)
        P_inf_g[:, :, g] = P_cell + np.outer(P_out, P_in) / (1.0 - P_inout)

    # Build scattering / fission data for the eigenvalue solver
    sig_s_fuel = fuel_mix.SigS[0].toarray()
    sig_s_mod = mod_mix.SigS[0].toarray()
    nu_sigf_fuel = fuel_mix.SigP
    nu_sigf_mod = mod_mix.SigP

    return _kinf_from_cp(
        P_inf_g=P_inf_g, sig_t_all=sig_t_all, V_arr=t_arr,
        sig_s_mats=[sig_s_fuel, sig_s_mod],
        nu_sig_f_mats=[nu_sigf_fuel, nu_sigf_mod],
        chi_mats=[fuel_mix.chi, mod_mix.chi],
    )


def run_do_slab_benchmarks():
    """Verify the SN solver on slab benchmarks.

    Builds SN meshes at several resolutions and compares against the
    analytical CP eigenvalue for the EXACT same effective geometry.
    Shows mesh convergence.
    """
    print("=" * 65)
    print("DISCRETE ORDINATES SLAB BENCHMARKS")
    print("=" * 65)

    # --- 1-group benchmark ---
    fuel_1g = _make_mixture(
        sig_t=np.array([1.0]), sig_c=np.array([0.2]),
        sig_f=np.array([0.3]), nu=np.array([2.5]),
        chi=np.array([1.0]), sig_s=np.array([[0.5]]),
    )
    mod_1g = _make_mixture(
        sig_t=np.array([2.0]), sig_c=np.array([0.1]),
        sig_f=np.array([0.0]), nu=np.array([0.0]),
        chi=np.array([1.0]), sig_s=np.array([[1.9]]),
    )

    # --- 2-group benchmark ---
    fuel_2g = _make_mixture(
        sig_t=np.array([0.50, 1.00]),
        sig_c=np.array([0.01, 0.02]),
        sig_f=np.array([0.01, 0.08]),
        nu=np.array([2.50, 2.50]),
        chi=np.array([1.00, 0.00]),
        sig_s=np.array([[0.38, 0.10], [0.00, 0.90]]),
    )
    mod_2g = _make_mixture(
        sig_t=np.array([0.60, 2.00]),
        sig_c=np.array([0.02, 0.05]),
        sig_f=np.array([0.00, 0.00]),
        nu=np.array([0.00, 0.00]),
        chi=np.array([1.00, 0.00]),
        sig_s=np.array([[0.40, 0.18], [0.00, 1.95]]),
    )

    for label, fuel, mod in [
        ("1G slab SN", fuel_1g, mod_1g),
        ("2G slab SN", fuel_2g, mod_2g),
    ]:
        ng = fuel.ng
        materials = {2: fuel, 0: mod}

        print(f"\n  {label} (mesh convergence):")
        print(f"  {'delta':>8s}  {'nx':>4s}  {'t_fuel':>7s}  {'t_mod':>7s}  "
              f"{'k_analytical':>12s}  {'k_SN':>12s}  {'Error':>10s}")
        print("  " + "-" * 70)

        for delta in [0.1, 0.05, 0.02, 0.01]:
            n_fuel = max(2, round(0.5 / delta))
            n_mod = max(2, round(0.5 / delta))

            geom, t_fuel, t_mod = _slab_geom_for_do(n_fuel, n_mod, delta)

            k_ref = _analytical_slab_kinf(fuel, mod, t_fuel, t_mod)

            result = solve_discrete_ordinates(
                materials, geom,
                params=DOParams(max_outer=300, bicgstab_tol=1e-6),
            )

            err = abs(result.keff - k_ref)
            print(f"  {delta:8.3f}  {geom.nx:4d}  {t_fuel:7.4f}  {t_mod:7.4f}  "
                  f"{k_ref:12.6f}  {result.keff:12.6f}  {err:10.2e}")

    print()


def run_sn_1d_benchmarks():
    """Verify the 1D SN solver with Gauss-Legendre quadrature.

    Tests:
    1. Homogeneous benchmarks (exact analytical reference)
    2. Heterogeneous slab convergence (spatial and angular)
       with observed order-of-accuracy
    """
    from sn_1d import (
        GaussLegendreQuadrature, Slab1DGeometry, solve_sn_1d,
    )

    print("=" * 78)
    print("1D SN (GAUSS-LEGENDRE) BENCHMARKS")
    print("=" * 78)

    # --- 1. Homogeneous benchmarks (exact reference) ---
    print("\n  Homogeneous infinite medium (SN with reflective BCs):")
    print(f"  {'Benchmark':<22s}  {'Groups':>6s}  {'Analytical':>10s}  "
          f"{'SN 1D':>10s}  {'Error':>10s}  {'OK':>3s}")
    print("  " + "-" * 68)

    for label, bench_fn in [
        ("1G homogeneous", benchmark_1g_homogeneous),
        ("2G homogeneous", benchmark_2g_homogeneous),
        ("4G homogeneous", benchmark_4g_homogeneous),
    ]:
        mix, k_analytical = bench_fn()
        ng = mix.ng
        materials = {i: mix for i in range(3)}
        geom = Slab1DGeometry.homogeneous(20, 2.0, mat_id=0)
        quad = GaussLegendreQuadrature.gauss_legendre(8)
        result = solve_sn_1d(materials, geom, quad)
        err = abs(result.keff - k_analytical)
        ok = "✓" if err < 1e-4 else "✗"
        print(f"  {label:<22s}  {ng:>6d}  {k_analytical:10.6f}  "
              f"{result.keff:10.6f}  {err:10.2e}  {ok:>3s}")

    # --- 2. Heterogeneous slab convergence ---
    # Cross-section data (same as benchmark_1g_slab / benchmark_2g_slab)
    fuel_1g = _make_mixture(
        sig_t=np.array([1.0]), sig_c=np.array([0.2]),
        sig_f=np.array([0.3]), nu=np.array([2.5]),
        chi=np.array([1.0]), sig_s=np.array([[0.5]]),
    )
    mod_1g = _make_mixture(
        sig_t=np.array([2.0]), sig_c=np.array([0.1]),
        sig_f=np.array([0.0]), nu=np.array([0.0]),
        chi=np.array([1.0]), sig_s=np.array([[1.9]]),
    )

    fuel_2g = _make_mixture(
        sig_t=np.array([0.50, 1.00]),
        sig_c=np.array([0.01, 0.02]),
        sig_f=np.array([0.01, 0.08]),
        nu=np.array([2.50, 2.50]),
        chi=np.array([1.00, 0.00]),
        sig_s=np.array([[0.38, 0.10], [0.00, 0.90]]),
    )
    mod_2g = _make_mixture(
        sig_t=np.array([0.60, 2.00]),
        sig_c=np.array([0.02, 0.05]),
        sig_f=np.array([0.00, 0.00]),
        nu=np.array([0.00, 0.00]),
        chi=np.array([1.00, 0.00]),
        sig_s=np.array([[0.40, 0.18], [0.00, 1.95]]),
    )

    t_fuel, t_mod = 0.5, 0.5
    sn_params = dict(max_outer=300, max_inner=500, inner_tol=1e-10)

    # --- Spatial convergence study (1G only, S16 fixed) ---
    label, fuel, mod = "1G slab", fuel_1g, mod_1g
    materials = {2: fuel, 0: mod}

    print(f"\n  {label} — spatial convergence (S16, diamond-difference O(h²)):")
    print(f"  {'dx':>8s}  {'N_cells':>7s}  {'keff':>12s}  {'Error':>10s}  {'Order':>6s}")
    print("  " + "-" * 50)

    keffs_sp = []
    dxs = []
    for n_per in [5, 10, 20, 40]:
        geom = Slab1DGeometry.from_benchmark(
            n_fuel=n_per, n_mod=n_per,
            t_fuel=t_fuel, t_mod=t_mod,
        )
        quad = GaussLegendreQuadrature.gauss_legendre(16)
        result = solve_sn_1d(materials, geom, quad, **sn_params)
        keffs_sp.append(result.keff)
        dxs.append(t_fuel / n_per)

    # Richardson extrapolation (O(h²), ratio 2) using two finest meshes
    k_ref_sp = keffs_sp[-1] + (keffs_sp[-1] - keffs_sp[-2]) / 3.0
    for i, (dx, k) in enumerate(zip(dxs, keffs_sp)):
        err = abs(k - k_ref_sp)
        if i > 0 and abs(keffs_sp[i - 1] - k_ref_sp) > 0 and err > 0:
            p = np.log(abs(keffs_sp[i - 1] - k_ref_sp) / err) / np.log(dxs[i - 1] / dx)
            order_str = f"{p:6.2f}"
        else:
            order_str = "   ---"
        print(f"  {dx:8.4f}  {2 * round(t_fuel / dx):>7d}  "
              f"{k:12.8f}  {err:10.2e}  {order_str}")
    print(f"  Richardson ref = {k_ref_sp:.8f}")

    # --- Angular convergence study (1G only, 40 cells fixed) ---
    print(f"\n  {label} — angular convergence (40 cells/region, GL spectral):")
    print(f"  {'N_ord':>6s}  {'keff':>12s}  {'Error':>10s}  {'Order':>6s}")
    print("  " + "-" * 42)

    keffs_ang = []
    n_ords = [4, 8, 16, 32]
    for N_ord in n_ords:
        geom = Slab1DGeometry.from_benchmark(
            n_fuel=40, n_mod=40,
            t_fuel=t_fuel, t_mod=t_mod,
        )
        quad = GaussLegendreQuadrature.gauss_legendre(N_ord)
        result = solve_sn_1d(materials, geom, quad, **sn_params)
        keffs_ang.append(result.keff)

    k_ref_ang = keffs_ang[-1]
    for i, (N, k) in enumerate(zip(n_ords, keffs_ang)):
        err = abs(k - k_ref_ang)
        if i > 0 and abs(keffs_ang[i - 1] - k_ref_ang) > 0 and err > 0:
            p = np.log(abs(keffs_ang[i - 1] - k_ref_ang) / err) \
                / np.log(N / n_ords[i - 1])
            order_str = f"{p:6.2f}"
        else:
            order_str = "   ---"
        print(f"  {N:>6d}  {k:12.8f}  {err:10.2e}  {order_str}")

    # CP white-BC reference for comparison
    k_cp = _analytical_slab_kinf(fuel, mod, t_fuel, t_mod)
    print(f"\n  CP (white BC) ref = {k_cp:.8f}  "
          f"(SN-CP gap = {abs(k_ref_ang - k_cp):.2e}, due to white-BC approx.)")

    # --- 2G slab: single point at moderate resolution ---
    label2, fuel2, mod2 = "2G slab", fuel_2g, mod_2g
    materials2 = {2: fuel2, 0: mod2}
    geom2 = Slab1DGeometry.from_benchmark(
        n_fuel=20, n_mod=20, t_fuel=t_fuel, t_mod=t_mod,
    )
    quad2 = GaussLegendreQuadrature.gauss_legendre(16)
    result2 = solve_sn_1d(materials2, geom2, quad2, **sn_params)
    k_cp2 = _analytical_slab_kinf(fuel2, mod2, t_fuel, t_mod)
    print(f"\n  2G slab SN-1D (S16, 40 cells): keff = {result2.keff:.8f}  "
          f"CP ref = {k_cp2:.8f}  gap = {abs(result2.keff - k_cp2):.2e}")

    print()


def run_moc_homogeneous_benchmarks():
    """Verify the MOC solver against homogeneous benchmarks."""
    from method_of_characteristics import MoCGeometry, solve_moc

    print("=" * 65)
    print("METHOD OF CHARACTERISTICS HOMOGENEOUS BENCHMARKS")
    print("=" * 65)
    print(f"  {'Benchmark':<22s}  {'Groups':>6s}  {'Analytical':>10s}  "
          f"{'MOC':>10s}  {'Error':>10s}  {'OK':>3s}")
    print("  " + "-" * 62)

    for label, bench_fn in [
        ("1G homogeneous", benchmark_1g_homogeneous),
        ("2G homogeneous", benchmark_2g_homogeneous),
    ]:
        mix, k_analytical = bench_fn()
        materials = {0: mix, 1: mix, 2: mix}
        geom = MoCGeometry.default_pwr()
        result = solve_moc(materials, geom, max_outer=200)
        err = abs(result.keff - k_analytical)
        ok = "✓" if err < 1e-2 else "✗"
        print(f"  {label:<22s}  {mix.ng:>6d}  {k_analytical:10.6f}  "
              f"{result.keff:10.6f}  {err:10.2e}  {ok:>3s}")

    print()


def run_mc_homogeneous_benchmarks():
    """Verify the Monte Carlo solver against homogeneous benchmarks."""
    from monte_carlo import MCParams, solve_monte_carlo

    print("=" * 65)
    print("MONTE CARLO HOMOGENEOUS BENCHMARKS")
    print("=" * 65)
    print(f"  {'Benchmark':<16s}  {'Groups':>6s}  {'Analytical':>10s}  "
          f"{'MC':>10s}  {'sigma':>8s}  {'|z|':>6s}  {'OK':>3s}")
    print("  " + "-" * 65)

    for label, bench_fn in [
        ("1G homogeneous", benchmark_1g_homogeneous),
        ("2G homogeneous", benchmark_2g_homogeneous),
    ]:
        mix, k_analytical = bench_fn()
        materials = {0: mix, 1: mix, 2: mix}
        # Fewer histories for 2G (slower per neutron due to scattering)
        n_active = 200 if mix.ng > 1 else 500
        params = MCParams(n_neutrons=200, n_inactive=50,
                          n_active=n_active, seed=42)
        result = solve_monte_carlo(materials, params)
        z_score = abs(result.keff - k_analytical) / max(result.sigma, 1e-10)
        ok = "✓" if z_score < 5.0 else "✗"
        print(f"  {label:<16s}  {mix.ng:>6d}  {k_analytical:10.6f}  "
              f"{result.keff:10.6f}  {result.sigma:8.5f}  {z_score:6.2f}  {ok:>3s}")

    print()


def run_diffusion_benchmarks():
    """Verify the 1D diffusion solver against analytical buckling eigenvalue."""
    from diffusion_1d import CoreGeometry, TwoGroupXS, solve_diffusion_1d

    print("=" * 78)
    print("1D DIFFUSION BARE-SLAB BENCHMARKS")
    print("=" * 78)

    # Use fuel XS from the default diffusion problem
    fuel_xs = TwoGroupXS(
        transport=np.array([0.2181, 0.7850]),
        absorption=np.array([0.0096, 0.0959]),
        fission=np.array([0.0024, 0.0489]),
        production=np.array([0.0061, 0.1211]),
        chi=np.array([1.0, 0.0]),
        scattering=np.array([0.0160, 0.0]),
    )

    fuel_height = 50.0  # cm — short slab so FD error is visible

    # Analytical keff: 2-group diffusion with buckling B² = (π/H)²
    # The FD scheme sets zero flux at the physical boundary.
    D_coeff = 1.0 / (3.0 * fuel_xs.transport)
    B2 = (np.pi / fuel_height) ** 2
    A_an = np.diag(D_coeff * B2 + fuel_xs.absorption + fuel_xs.scattering) \
        - np.array([[0.0, 0.0], [fuel_xs.scattering[0], 0.0]])
    F_an = np.outer(fuel_xs.chi, fuel_xs.production)
    M_an = np.linalg.solve(A_an, F_an)
    k_analytical = float(np.max(np.real(np.linalg.eigvals(M_an))))

    # Spatial convergence study
    print(f"\n  Bare fuel slab (H={fuel_height} cm, vacuum BCs, 2-group)")
    print(f"  Analytical keff = {k_analytical:.6f}")
    print(f"  Spatial convergence (O(h²) expected):")
    print(f"  {'dz':>8s}  {'N_cells':>7s}  "
          f"{'k_diffusion':>12s}  {'Error':>10s}  {'Order':>6s}")
    print("  " + "-" * 52)

    keffs = []
    dzs = [5.0, 2.5, 1.25, 0.625]
    for dz in dzs:
        geom = CoreGeometry(
            bot_refl_height=0.0, fuel_height=fuel_height,
            top_refl_height=0.0, dz=dz,
        )
        result = solve_diffusion_1d(
            geom=geom, reflector_xs=fuel_xs, fuel_xs=fuel_xs,
        )
        keffs.append(result.keff)

    for i, (dz, k) in enumerate(zip(dzs, keffs)):
        err = abs(k - k_analytical)
        if i > 0 and abs(keffs[i - 1] - k_analytical) > 0 and err > 0:
            p = np.log(abs(keffs[i - 1] - k_analytical) / err) / np.log(dzs[i - 1] / dz)
            order_str = f"{p:6.2f}"
        else:
            order_str = "   ---"
        nc = int(fuel_height / dz)
        print(f"  {dz:8.2f}  {nc:>7d}  {k:12.6f}  {err:10.2e}  {order_str}")

    print()


def main():
    print()
    run_homogeneous_benchmarks()
    run_slab_benchmarks()
    run_cylinder_benchmarks()
    run_do_slab_benchmarks()
    run_sn_1d_benchmarks()
    run_diffusion_benchmarks()
    run_moc_homogeneous_benchmarks()
    run_mc_homogeneous_benchmarks()

    print("=" * 65)
    print("VERIFICATION COMPLETE")
    print("=" * 65)


if __name__ == "__main__":
    main()
