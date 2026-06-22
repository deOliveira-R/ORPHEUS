"""Shared plotting utilities for ORPHEUS reactor physics solvers.

Common functions used across multiple solver modules. Module-specific
plots (e.g. DO spatial flux) remain in their respective directories.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Patch, Rectangle

if TYPE_CHECKING:
    from orpheus.homogeneous.solver import HomogeneousResult
    from orpheus.moc.solver import MoCResult
    from orpheus.moc.geometry import MOCMesh
    from orpheus.mc.solver import MCResult


# ═══════════════════════════════════════════════════════════════════════
# Generic helpers
# ═══════════════════════════════════════════════════════════════════════

def plot_2d_field(
    nx: int, ny: int, delta: float,
    field: np.ndarray,
    title: str,
    filepath: Path,
) -> None:
    """2D colored-rectangle plot of a scalar field on a Cartesian mesh."""
    aspect_ratio = max(ny / nx, 0.3)
    fig, ax = plt.subplots(figsize=(8, max(3, 8 * aspect_ratio)))
    ax.set_aspect("equal")
    ax.set_xlim(0, delta * (nx - 1))
    ax.set_ylim(-delta * (ny - 1), 0)
    ax.set_xlabel("Distance from the cell centre (cm)")
    ax.set_ylabel("Distance from the cell centre (cm)")
    ax.set_title(title)

    fmin, fmax = field.min(), field.max()

    y = 0.0
    for iy in range(ny):
        height = delta - delta / 2 * (iy == 0 or iy == ny - 1)
        y -= height
        x = 0.0
        for ix in range(nx):
            width = delta - delta / 2 * (ix == 0 or ix == nx - 1)
            c = (field[ix, iy] - fmin) / max(fmax - fmin, 1e-30)
            ax.add_patch(Rectangle(
                (x, y), width, height,
                facecolor=(c, 0, 1 - c), edgecolor=(0.5, 0.5, 0.5), linewidth=0.2,
            ))
            x += width

    fig.tight_layout()
    fig.savefig(filepath)
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════
# Homogeneous reactor
# ═══════════════════════════════════════════════════════════════════════

def plot_spectrum(
    result: HomogeneousResult,
    title: str = "",
    output_dir: Path | str = ".",
    prefix: str = "spectrum",
) -> None:
    """Generate flux-per-energy and flux-per-lethargy plots.

    Skipped (with a warning) when the result has no physical energy
    grid (``result.eg_mid is None``). Per-energy plotting is meaningless
    for synthetic / Sood-style XS — consumers should plot only
    ``result.flux`` directly in that case.
    """
    if result.eg_mid is None:
        import warnings
        warnings.warn(
            "plot_spectrum: result has no energy grid (synthetic XS); "
            "skipping per-energy and per-lethargy plots.",
            stacklevel=2,
        )
        return
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.loglog(result.eg_mid, result.flux_per_energy)
    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel(r"Neutron flux per unit energy (cm$^{-2}$s$^{-1}$eV$^{-1}$)")
    ax.set_title(f"{title} — Flux per unit energy" if title else "Flux per unit energy")
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_dir / f"{prefix}_flux_energy.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.semilogx(result.eg_mid, result.flux_per_lethargy)
    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel(r"Neutron flux per unit lethargy (cm$^{-2}$s$^{-1}$)")
    ax.set_title(f"{title} — Flux per unit lethargy" if title else "Flux per unit lethargy")
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_dir / f"{prefix}_flux_lethargy.pdf")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════
# Method of Characteristics
# ═══════════════════════════════════════════════════════════════════════

def plot_moc_rays(
    geom: MOCMesh,
    output_dir: Path | str = ".",
    filename: str = "MOC_01_tracks.pdf",
) -> None:
    """Plot the characteristic ray tracks across the pin cell.

    Each track is the straight chord from its entry point to its exit
    point on the equal-area square boundary (``pitch = R_outer * sqrt(pi)``,
    the inverse Wigner-Seitz transform), coloured by azimuthal angle.  The
    concentric material rings (``geom.radii``) are overlaid for context.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    pitch = geom.pitch
    cx, cy = pitch / 2.0, pitch / 2.0

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_aspect("equal")
    ax.set_xlim(0.0, pitch)
    ax.set_ylim(0.0, pitch)
    ax.set_xlabel("x (cm)")
    ax.set_ylabel("y (cm)")
    ax.set_title("Characteristic ray tracks")

    cmap = plt.get_cmap("hsv")
    n_azi = max(1, geom.quad.n_azi)
    for track in geom.tracks:
        (x0, y0), (x1, y1) = track.entry_point, track.exit_point
        ax.plot(
            [x0, x1], [y0, y1],
            color=cmap(track.azi_index / n_azi),
            linewidth=0.4, alpha=0.6,
        )

    # Material-ring + cell-boundary overlay (outline only).
    for r in geom.radii:
        ax.add_patch(Circle((cx, cy), float(r), fill=False,
                            edgecolor="black", linewidth=0.6))
    ax.add_patch(Rectangle((0.0, 0.0), pitch, pitch, fill=False,
                          edgecolor="black"))

    fig.tight_layout()
    fig.savefig(output_dir / filename)
    plt.close(fig)


def plot_moc_mesh(
    geom: MOCMesh,
    output_dir: Path | str = ".",
    filename: str = "MOC_02_mesh.pdf",
) -> None:
    """Plot the pin-cell material map.

    The cylindrical Wigner-Seitz pin cell is a set of concentric material
    rings (``geom.radii``) inside the equal-area square outer boundary; the
    outermost region (moderator) fills the square corners.  Regions are
    painted largest-radius first (painter's algorithm) and coloured by
    material id.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    pitch = geom.pitch
    cx, cy = pitch / 2.0, pitch / 2.0
    mat_ids = geom.region_mat_ids

    unique_mats = sorted({int(m) for m in mat_ids})
    cmap = plt.get_cmap("tab10")
    color_of = {m: cmap(i % 10) for i, m in enumerate(unique_mats)}

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_aspect("equal")
    ax.set_xlim(0.0, pitch)
    ax.set_ylim(0.0, pitch)
    ax.set_xlabel("x (cm)")
    ax.set_ylabel("y (cm)")
    ax.set_title("Pin-cell material map")

    # Outermost region fills the equal-area square; inner rings paint on
    # top, largest radius first, leaving each annulus its own colour.
    ax.add_patch(Rectangle((0.0, 0.0), pitch, pitch,
                          facecolor=color_of[int(mat_ids[-1])],
                          edgecolor="black"))
    for k in range(geom.n_regions - 2, -1, -1):
        ax.add_patch(Circle((cx, cy), float(geom.radii[k]),
                            facecolor=color_of[int(mat_ids[k])],
                            edgecolor="black", linewidth=0.5))

    handles = [
        Patch(facecolor=color_of[m], edgecolor="black", label=f"material {m}")
        for m in unique_mats
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=8)

    fig.tight_layout()
    fig.savefig(output_dir / filename)
    plt.close(fig)


def plot_moc_convergence(
    result: MoCResult,
    output_dir: Path | str = ".",
) -> None:
    """Plot keff convergence for MoC."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots()
    ax.plot(range(1, len(result.keff_history) + 1), result.keff_history, "-or", markersize=3)
    ax.set_xlabel("Iteration number")
    ax.set_ylabel("k-effective")
    ax.grid(True)
    fig.tight_layout()
    fig.savefig(output_dir / "MOC_03_keff.pdf")
    plt.close(fig)


def plot_moc_spectra(
    result: MoCResult,
    output_dir: Path | str = ".",
) -> None:
    """Plot neutron spectra per unit lethargy in fuel, cladding, coolant.

    Skipped (with a warning) when the result has no physical energy
    grid (``result.eg is None``).  Per-lethargy plotting requires a
    real grid; synthetic Sood-style mixtures don't carry one.
    """
    if result.eg is None:
        import warnings
        warnings.warn(
            "plot_moc_spectra: result has no energy grid (synthetic XS); "
            "skipping per-lethargy plot.",
            stacklevel=2,
        )
        return
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    eg = result.eg
    eg_mid = 0.5 * (eg[:-1] + eg[1:])
    du = np.abs(np.log(eg[1:] / eg[:-1]))

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.semilogx(eg_mid, result.flux_fuel / du, "-r", label="Fuel")
    ax.semilogx(eg_mid, result.flux_clad / du, "-g", label="Cladding")
    ax.semilogx(eg_mid, result.flux_cool / du, "-b", label="Coolant")
    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel("Neutron flux per unit lethargy (a.u.)")
    ax.legend(loc="upper left")
    ax.grid(True)
    fig.tight_layout()
    fig.savefig(output_dir / "MOC_04_flux_lethargy.pdf")
    plt.close(fig)


def plot_moc_spatial_flux(
    result: MoCResult,
    output_dir: Path | str = ".",
) -> None:
    """Plot thermal/resonance/fast flux along cell centerline and 2D."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    geom = result.geometry
    n = geom.n_cells
    x = np.arange(n) * geom.delta

    FI_T = result.scalar_flux[:, 0, :50].sum(axis=1)
    FI_R = result.scalar_flux[:, 0, 50:287].sum(axis=1)
    FI_F = result.scalar_flux[:, 0, 287:].sum(axis=1)

    fig, ax = plt.subplots()
    ax.plot(x, FI_F, "-or", label="Fast", markersize=3)
    ax.plot(x, FI_R, "-og", label="Resonance", markersize=3)
    ax.plot(x, FI_T, "-ob", label="Thermal", markersize=3)
    ax.set_xlabel("Distance from the cell centre (cm)")
    ax.set_ylabel("Neutron flux (a.u.)")
    ax.legend(loc="upper left")
    ax.grid(True)
    fig.tight_layout()
    fig.savefig(output_dir / "MOC_05_flux_cell.pdf")
    plt.close(fig)

    fun_T = result.scalar_flux[:, :, :50].sum(axis=2)
    fun_R = result.scalar_flux[:, :, 50:355].sum(axis=2)
    fun_F = result.scalar_flux[:, :, 355:].sum(axis=2)

    plot_2d_field(n, n, geom.delta, fun_T,
                  "Thermal flux distribution", output_dir / "MOC_06_flux_thermal.pdf")
    plot_2d_field(n, n, geom.delta, fun_R,
                  "Resonance flux distribution", output_dir / "MOC_07_flux_resonance.pdf")
    plot_2d_field(n, n, geom.delta, fun_F,
                  "Fast flux distribution", output_dir / "MOC_08_flux_fast.pdf")


# ═══════════════════════════════════════════════════════════════════════
# Monte Carlo
# ═══════════════════════════════════════════════════════════════════════

def plot_mc_keff(
    result: MCResult,
    output_dir: Path | str = ".",
) -> None:
    """Plot keff convergence with uncertainty bands."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    iters = np.arange(1, len(result.keff_history) + 1)
    fig, ax = plt.subplots()
    ax.plot(iters, result.keff_history, "-r", label=r"$k_{eff}$")
    ax.plot(iters, result.keff_history + result.sigma_history, "--b",
            label=r"$k_{eff} \pm \sigma$")
    ax.plot(iters, result.keff_history - result.sigma_history, "--b")
    ax.set_xlabel("Active cycle number")
    ax.set_ylabel("k-effective")
    ax.legend()
    ax.grid(True)
    fig.tight_layout()
    fig.savefig(output_dir / "MC_01_keff.pdf")
    plt.close(fig)


def plot_mc_spectrum(
    result: MCResult,
    output_dir: Path | str = ".",
) -> None:
    """Plot cell-averaged neutron flux per unit lethargy.

    Skipped (with a warning) when the result has no physical energy
    grid (``result.eg_mid is None``).  Per-lethargy plotting is
    undefined for synthetic / Sood-style XS.
    """
    if result.eg_mid is None or result.flux_per_lethargy is None:
        import warnings
        warnings.warn(
            "plot_mc_spectrum: result has no energy grid (synthetic XS); "
            "skipping per-lethargy plot.",
            stacklevel=2,
        )
        return
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.semilogx(result.eg_mid, result.flux_per_lethargy)
    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel("Neutron flux per unit lethargy (a.u.)")
    ax.grid(True)
    fig.tight_layout()
    fig.savefig(output_dir / "MC_02_flux_lethargy.pdf")
    plt.close(fig)
