# ORPHEUS — Open Reactor Physics Educational University System

Provides practical tasks covering the main steps of nuclear reactor analysis:
cross section processing, neutron transport, diffusion, fuel behaviour,
thermal hydraulics, and reactor kinetics.

## Getting Started

### Prerequisites

- Python 3.12+ (3.14 recommended)
- git-lfs (for nuclear data files)

### Installation

```bash
git clone git@github.com:deOliveira-R/reactor_physics_lessons.git
cd reactor_physics_lessons
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### Convert nuclear data

The repository ships GENDF (.GXS) cross section files via git-lfs.
Convert them to HDF5 before running any calculations:

```bash
cd data/micro_xs
python convert_gxs_to_hdf5.py
```

This reads the 12 `.GXS` files and produces one `.h5` file per isotope
containing all available temperatures.

## Architecture

All deterministic eigenvalue solvers satisfy the `EigenvalueSolver`
protocol defined in `numerics/eigenvalue.py` and share a generic
`power_iteration()` function.  Power iteration converges to the
**dominant eigenvalue** (k_eff) and its eigenvector — the
**fundamental mode** of the flux distribution.  This is the unique
non-negative solution (Perron-Frobenius); higher harmonics change sign
in space and decay with each iteration.  Absolute flux normalization
(by power, integral flux, etc.) is a separate post-processing step.

```
numerics/
    eigenvalue.py        EigenvalueSolver protocol + power_iteration()

data/
    micro_xs/            421-group microscopic XS (GENDF/HDF5), Isotope dataclass
    macro_xs/
        mixture.py       Mixture dataclass (SigT, SigS, absorption_xs, ...)
        cell_xs.py       CellXS dataclass + assemble_cell_xs() — shared by all spatial solvers
        benchmarks.py    Synthetic benchmarks with analytical eigenvalues
        recipes.py       Material recipes (UO2, Zircaloy, borated water)
        sigma_zeros.py   Sigma-zero self-shielding iteration
        interpolation.py XS interpolation at converged background XS
    materials/           MATPRO correlations + water/steam properties (pyXSteam)

01-08.*                  Lecture modules (one per physics domain)
09.Collision.Probability Collision probability solvers + formal verification suite
```

### Terminology

The codebase uses 3D mesh terminology consistently, even in reduced dimensions (1D and 2D are spatially degenerate 3D cases):

- **cell / volume** — where material properties and volume-averaged fields live
- **face** — interface between cells, where diffusion coefficients, currents, and gradients are defined

## Modules

Each numbered folder is a self-contained lecture module.  Run from the
repository root:

```bash
python 01.Homogeneous.Reactors/run_homogeneous.py
```

| Module | Description | Solver class | Entry script |
|--------|-------------|-------------|-------------|
| **00.Demo** | Central Limit Theorem and spherical harmonics demos | — | `central_limit_theorem.py`, `spherical_harmonics.py` |
| **01.Homogeneous.Reactors** | Infinite medium eigenvalue problem (421 groups) | `HomogeneousSolver` | `run_homogeneous.py` |
| **02.Discrete.Ordinates** | 2D SN transport with Lebedev quadrature (110 ordinates) | `DO2DSolver` | `run_discrete_ordinates.py` |
| | 1D SN transport with Gauss-Legendre quadrature | `SN1DSolver` | — (used via verification) |
| **03.Method.Of.Characteristics** | 2D MoC transport with 8 ray directions | `MoCSolver` | `run_moc.py` |
| **04.Monte.Carlo** | Monte Carlo with Woodcock delta tracking | — (stochastic) | `run_monte_carlo.py` |
| **05.Diffusion.1D** | Two-group 1D axial diffusion for a PWR subassembly | `DiffusionSolver` | `run_diffusion_1d.py` |
| **06.Fuel.Behaviour** | 1D radial thermo-mechanical fuel rod analysis (6 years) | — | `run_fuel_behaviour.py` |
| **07.Thermal.Hydraulics** | Coupled TH + fuel mechanics under LOCA conditions (600 s) | — | `run_thermal_hydraulics.py` |
| **08.Reactor.Kinetics.0D** | Point kinetics + TH + fuel mechanics under RIA conditions | — | `run_reactor_kinetics.py` |
| **09.Collision.Probability** | CP method for slab and concentric cylindrical geometries | `CPSolver` | `run_cp_slab.py`, `run_cp_concentric.py` |

### Collision Probability

Module 09 provides two CP solvers sharing a single `CPSolver` class:

- **Slab** (`solve_cp_slab`) — 1D half-cell with E_3 exponential-integral kernel
- **Concentric** (`solve_cp_concentric`) — Wigner-Seitz cylindrical cell with Ki_3/Ki_4 Bickley-Naylor kernel

The CP matrices are geometry-specific; the eigenvalue iteration is
geometry-independent.

## Nuclear Data

Cross sections are in the IAEA 421-group GENDF format, downloaded from:
https://www-nds.iaea.org/ads/adsgendf.html

Isotopes included: H-1, B-10, B-11, O-16, Na-23, U-235, U-238,
Zr-90, Zr-91, Zr-92, Zr-94, Zr-96.

## SQA

### MATLAB comparison

Modules 01-07 are benchmarked against the original MATLAB results:

| Module | Quantity | Python | MATLAB | Match |
|--------|----------|--------|--------|-------|
| 01 Homogeneous (aqueous) | k_inf | 1.03596 | 1.03596 | exact |
| 01 Homogeneous (PWR) | k_inf | 1.01357 | 1.01357 | exact |
| 02 Discrete Ordinates | keff | 1.04190 | 1.04188 | 2e-5 |
| 03 Method of Characteristics | keff | 1.04923 | 1.04923 | exact |
| 04 Monte Carlo | keff | 1.038 +/- 0.002 | 1.035 +/- 0.002 | stochastic |
| 05 1D Diffusion | keff | 1.022171 | 1.022173 | 2e-6 |
| 06 Fuel Behaviour (t=1d) | Fuel center T | 1017.80 C | 1017.86 C | 0.006 C |

### Formal verification

A formal verification suite (`09.Collision.Probability/run_verification.py`)
tests all deterministic solvers against synthetic benchmarks with known
analytical eigenvalues and reports observed order of accuracy:

```bash
python 09.Collision.Probability/run_verification.py
```

| Solver | Benchmark | Reference | Result |
|--------|-----------|-----------|--------|
| Homogeneous | 1G/2G/4G infinite medium | Analytical eigenvalue | Machine precision |
| CP Slab | 1G/2G 2-region slab | E_3 CP matrix eigenvalue | < 2e-7 |
| CP Concentric | 1G/2G Wigner-Seitz | Ki_4 CP matrix eigenvalue | < 2e-6 |
| SN 1D (Gauss-Legendre) | Homogeneous slab | Analytical | Machine precision |
| SN 1D | Heterogeneous slab | Richardson extrapolation | O(h^2) spatial, spectral angular |
| Diffusion 1D | 2-group bare slab | Buckling eigenvalue | O(h^2) spatial |
| MOC | 1G/2G homogeneous pin cell | Analytical | Machine precision |
| Monte Carlo | 1G/2G homogeneous | Analytical | |z| < 5 sigma |

## Attribution

Based on the MATLAB educational system by Konstantin Mikityuk (PSI).
Python port and augmentation by Rodrigo de Oliveira.
