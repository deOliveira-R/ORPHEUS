#!/usr/bin/env python3
"""1-D two-group diffusion for a PWR subassembly — the modern operator path.

The CORE1D.m reference problem (reflector | fuel | reflector, dz = 5 cm)
posed through the #290 operator algebra: materials are ``Mixture``
objects, the core layout is the mesh's own ``mat_ids``, and the
boundary condition is the honestly-named ``zero_flux`` Dirichlet law
(what the MATLAB port called "vacuum" — in ORPHEUS, *vacuum* MEANS zero
incoming current, the Marshak condition; the MATLAB reference used the
φ_Γ = 0 idealization).

Reference MATLAB result: keff = 1.022173.

The legacy cross sections encode bit-identically (campaign ruling 4)
through :func:`~orpheus.derivations.common.xs_library.
mixture_from_diffusion_tables` — ``sig_t := transport`` with no P1
moment, so ``Σ_tr = Σ_t`` exactly and ``D = 1/(3·transport)``; see its
docstring for the full encoding. The physical-P1 re-baseline is the
campaign's close-out follow-up.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from orpheus.derivations.common.xs_library import mixture_from_diffusion_tables
from orpheus.diffusion import solve_diffusion_1d
from orpheus.geometry import BC
from orpheus.geometry.mesh import Mesh1D

OUTPUT = Path("results")

# ── CORE1D.m two-group tables (legacy transport-XS encoding) ─────────────

def _core1d_materials():
    """Reflector (id 0) + fuel (id 1) from MATLAB CORE1D.m."""
    reflector = mixture_from_diffusion_tables(dict(
        transport=np.array([0.3416, 0.9431]),
        absorption=np.array([0.0029, 0.0933]),
        fission=np.array([0.0, 0.0]),
        production=np.array([0.0, 0.0]),
        chi=np.array([0.0, 0.0]),
        scattering=np.array([2.4673e-04, 0.0]),
    ))
    fuel = mixture_from_diffusion_tables(dict(
        transport=np.array([0.2181, 0.7850]),
        absorption=np.array([0.0096, 0.0959]),
        fission=np.array([0.0024, 0.0489]),
        production=np.array([0.0061, 0.1211]),
        chi=np.array([1.0, 0.0]),
        scattering=np.array([0.0160, 0.0]),
    ))
    return {0: reflector, 1: fuel}


def main():
    print("=" * 70)
    print("1D TWO-GROUP DIFFUSION — PWR SUBASSEMBLY (operator algebra)")
    print("=" * 70)

    # 50 cm reflector | 300 cm fuel | 50 cm reflector, dz = 5 cm.
    dz = 5.0
    n_refl, n_fuel = 10, 60
    edges = np.arange(0.0, (2 * n_refl + n_fuel) * dz + dz / 2, dz)
    mat_ids = np.array([0] * n_refl + [1] * n_fuel + [0] * n_refl)
    mesh = Mesh1D(
        edges=edges, mat_ids=mat_ids,
        bc_left=BC("zero_flux"), bc_right=BC("zero_flux"),
    )

    result = solve_diffusion_1d(_core1d_materials(), mesh)

    print(f"\n  keff = {result.keff:.6f}  (MATLAB reference: 1.022173)")
    match = "YES" if abs(result.keff - 1.022173) < 1e-4 else "NO"
    print(f"  Match: {match}   ({len(result.keff_history)} outer iterations)")

    # The eigenmode is production-normalized (∫νΣf·φ dV = 1); rescaling
    # to an absolute flux at a target power is one multiplication.
    z_cells = 0.5 * (edges[:-1] + edges[1:])
    z_faces = edges
    flux = result.flux.bulk.values          # (2, n_cells)
    current = result.current                # (2, n_faces), axis-signed

    OUTPUT.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots()
    ax.plot(z_cells, flux[0, :], "-or", label="Fast", markersize=3)
    ax.plot(z_cells, flux[1, :], "-ob", label="Thermal", markersize=3)
    ax.set_xlabel("z (cm)")
    ax.set_ylabel(r"Neutron flux (arb. units, $\int\nu\Sigma_f\phi\,dV=1$)")
    ax.legend()
    ax.grid(True)
    fig.tight_layout()
    fig.savefig(OUTPUT / "DIF_01_flux.pdf")
    plt.close(fig)

    fig, ax = plt.subplots()
    ax.plot(z_faces, current[0, :], "-or", label="Fast", markersize=3)
    ax.plot(z_faces, current[1, :], "-ob", label="Thermal", markersize=3)
    ax.set_xlabel("z (cm)")
    ax.set_ylabel("Neutron net current (arb. units)")
    ax.legend()
    ax.grid(True)
    fig.tight_layout()
    fig.savefig(OUTPUT / "DIF_02_current.pdf")
    plt.close(fig)

    print(f"\n  Plots saved to {OUTPUT.resolve()}/")


if __name__ == "__main__":
    main()
