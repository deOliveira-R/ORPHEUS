# ORPHEUS — Open Reactor Physics Educational University System

A rigorous Python implementation of reactor physics solvers with a
comprehensive verification suite, backed by an educational example collection
covering the main steps of nuclear reactor analysis: cross section processing,
neutron transport, diffusion, fuel behaviour, thermal hydraulics, and reactor
kinetics.

## Getting Started

### Prerequisites

- Python 3.11+ (3.14 recommended)
- [`uv`](https://docs.astral.sh/uv/) (recommended) — fast, disk-efficient venvs:
  `brew install uv` (or `pipx install uv`, or `curl -LsSf https://astral.sh/uv/install.sh | sh`).
  Plain `python3 -m venv` works too.
- git-lfs (for nuclear data files)
- [GitHub CLI](https://cli.github.com/) (`gh`) — for issue tracking

### Installation

One command sets up a checkout (uses `uv` if present, else stdlib `venv`):

```bash
git clone git@github.com:deOliveira-R/ORPHEUS.git
cd ORPHEUS
./scripts/setup.sh            # .venv + editable .[dev] install + Sphinx/Nexus graph
                              #   add --no-docs to skip the (slow) docs build
source .venv/bin/activate     # optional — for bare python / pytest / nexus
```

…or by hand:

```bash
uv venv && uv pip install -e ".[dev]"                        # with uv (recommended)
python3 -m venv .venv && .venv/bin/pip install -e ".[dev]"   # stdlib fallback
```

**Working in a git worktree?** Run `./scripts/setup.sh` *inside* the worktree —
a venv per checkout. `orpheus` is an editable install, so the venv must be created
from the checkout for `import orpheus` to resolve to *that* checkout's code; a
shared or symlinked `.venv` would import the main clone instead. `uv` makes the
per-worktree venv cheap (deps are hardlinked, not recopied).

### Editor / type checking (optional)

The repo ships a Pyright config (`[tool.pyright]` in `pyproject.toml`, pointed at
the project `.venv`). For in-editor diagnostics or the Pyright LSP:

```bash
npm install -g pyright       # provides pyright + pyright-langserver
```

Run the install/setup above first — Pyright resolves third-party imports
(numpy/scipy/…) from `.venv`.

### Convert nuclear data

The repository ships GENDF (.GXS) cross section files via git-lfs.
Convert them to HDF5 before running any calculations:

```bash
cd orpheus/data/micro_xs
.venv/bin/python convert_gxs_to_hdf5.py
```

### Run tests

```bash
.venv/bin/python -m pytest                  # non-slow tests
.venv/bin/python -m pytest -m slow          # slow tests (Richardson, MC high-stats)
.venv/bin/python -m pytest -v               # all ~5000 tests
```

### Build documentation

```bash
.venv/bin/python -m orpheus.derivations.generate_rst
.venv/bin/python -m sphinx -b html docs docs/_build/html
open docs/_build/html/index.html
```

## Project Structure

The project is organized into two main areas:

### `orpheus/` — The Python Package

Importable, testable, pip-installable solver library.

```
orpheus/
    __init__.py
    plotting.py              Shared plotting utilities

    # ── Solvers ──────────────────────────────────────────────
    homogeneous/             Infinite medium eigenvalue problem (421 groups)
    sn/                      Discrete ordinates (SN) transport
        solver.py            Unified SN solver: 1D slab/curvilinear + 2D Cartesian; source iteration + Krylov
        sweep.py             Diamond-difference transport sweep
        sweep_graph.py       Upwind sweep dependency graph (octant wavefronts)
        operator.py          Typed transport operator algebra (L + C − S − F)
        geometry.py          Augmented SN mesh (SNMesh)
        axis.py              Dimension-agnostic axis primitives (Axis1D, FaceLabel)
    moc/                     Method of characteristics (2D)
    mc/                      Monte Carlo with Woodcock delta tracking
    cp/                      Collision probability (slab + cylindrical)
    diffusion/               Two-group 1D axial diffusion
    fuel/                    1D radial thermo-mechanical fuel rod analysis
    thermal_hydraulics/      Coupled TH + fuel mechanics (ODE + DAE)
    kinetics/                Point kinetics + TH + fuel mechanics

    # ── Infrastructure ───────────────────────────────────────
    numerics/
        eigenvalue.py        EigenvalueSolver protocol + power_iteration()
        quadrature/          Angular quadratures (GL, Lebedev, level-symmetric, product)
        field.py             Typed Field algebra (flux / source / residual roles)
        iteration.py         Source iteration + Krylov acceleration
        face_layout.py       Boundary face-layout descriptor
    data/
        micro_xs/            421-group microscopic XS (GENDF/HDF5), Isotope dataclass
        macro_xs/            Mixture, CellXS, recipes, self-shielding
        materials/           MATPRO correlations + water/steam properties
    geometry/                Mesh1D/Mesh2D, coordinate systems, factory functions
    derivations/             SymPy analytical derivations (verification truth)
```

All deterministic eigenvalue solvers satisfy the `EigenvalueSolver`
protocol defined in `orpheus.numerics.eigenvalue` and share a generic
`power_iteration()` function.

### `examples/` — Educational Entry Points

Demo scripts that teach reactor physics concepts. Each subdirectory is a
self-contained lesson:

```bash
cd examples/homogeneous
python demo_homogeneous.py
```

| Example | Description | Solver |
|---------|-------------|--------|
| `demo/` | Central Limit Theorem, spherical harmonics | — |
| `homogeneous/` | Infinite medium eigenvalue | `orpheus.homogeneous` |
| `discrete_ordinates/` | SN transport (1D slab, 2D Cartesian) | `orpheus.sn` |
| `method_of_characteristics/` | 2D MoC with ray tracing | `orpheus.moc` |
| `monte_carlo/` | Delta tracking, pluggable geometry | `orpheus.mc` |
| `diffusion/` | 1D axial diffusion for PWR subassembly | `orpheus.diffusion` |
| `fuel_behaviour/` | Fuel rod thermo-mechanics | `orpheus.fuel` |
| `thermal_hydraulics/` | LOCA transient analysis | `orpheus.thermal_hydraulics` |
| `reactor_kinetics/` | RIA transient (kinetics + TH + fuel) | `orpheus.kinetics` |
| `collision_probability/` | CP method (slab + cylindrical) | `orpheus.cp` |

### Other directories

```
tests/               pytest verification suite (~500 tests)
tools/research/      Literature search utilities (arXiv, OSTI, Scopus, ...)
matlab_archive/      Original MATLAB code by K. Mikityuk (PSI)
docs/                Sphinx documentation (theory + API)
```

## Verification

The verification suite uses SymPy-derived analytical references as the
single source of truth. Each solver has its own derivation from its own
equations — no cross-verification.

```bash
pytest tests/ -v    # run all verification tests
```

| Method | Geometry | Groups × Regions | Reference type |
|--------|----------|-------------------|---------------|
| Homogeneous | — | 1/2/4 × 1 | Analytical (matrix eigenvalue) |
| SN 1D | Slab | 1/2/4 × 1,2,4 | Analytical + MMS O(h²) |
| CP Slab | Slab | 1/2/4 × 1,2,4 | Semi-analytical (E₃ eigenvalue) |
| CP Cylinder | Cyl1D | 1/2/4 × 1,2,4 | Semi-analytical (Ki₄ eigenvalue) |
| MOC | Cyl1D | 1/2/4 × 1,2,4 | Analytical + Richardson |
| MC | Cyl1D | 1/2/4 × 1,2,4 | Analytical + CP reference |
| Diffusion | Slab | 2 × 1,2 | Analytical (buckling) + Richardson |

Unit tests verify structural properties: CP conservation/reciprocity,
SN particle balance/flux symmetry, diffusion vacuum BCs.

## Nuclear Data

Cross sections are in the IAEA 421-group GENDF format, downloaded from:
https://www-nds.iaea.org/ads/adsgendf.html

Isotopes included: H-1, B-10, B-11, O-16, Na-23, U-235, U-238,
Zr-90, Zr-91, Zr-92, Zr-94, Zr-96.

## Attribution

Based on the MATLAB educational system by Konstantin Mikityuk (PSI).
Python port and augmentation by Rodrigo de Oliveira.
