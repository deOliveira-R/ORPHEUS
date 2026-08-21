r"""CS4a K1 — the interaction kernels (gates G1.1–G1.9).

The kernels (:mod:`orpheus.transport.kernels`) are representation-free
per-material physics data — the datum the CS4a bindings realize as bound
operators. These gates are software invariants of that posing
(``foundation``, the ``test_operator_spaces.py`` placement argument);
none carries ``verifies(...)``.

**Fixture discipline** (verification plan §3, ``lessons`` L1): every
mixture is built DIRECTLY through the campaign's one builder
(:func:`tests.sn.architecture._config.anisotropic_mixture`) — never
``make_mixture``, which hardcodes ``SigL = 0``, defaults ``Sig2`` to
zero, and offers no P1 channel, making the ℓ≥1 and (n,2n) mutation arms
vacuously green. Fixtures are function-scoped builders (fresh per call)
so the G1.5 carrier-mutation leg cannot pollute a shared object.
"""

from __future__ import annotations

import ast
import dataclasses
import inspect
from pathlib import Path

import numpy as np
import pytest

import orpheus.transport.kernels as kernels_module
from orpheus.data.energy_grid import EnergyGrid
from orpheus.data.macro_xs.mixture import Mixture
from orpheus.derivations.common.xs_library import get_mixture
from orpheus.numerics.axis import EnergyAxis
from orpheus.transport.kernels import FissionKernel, N2NKernel, ScatteringKernel
from orpheus.transport.mesh.material_mesh import MaterialMesh
from orpheus.transport.operators.isotropic_scattering import (
    IsotropicN2N,
    IsotropicScattering,
)
from tests.sn.architecture._config import anisotropic_mixture

pytestmark = pytest.mark.foundation


# ═════════════════════════════════════════════════════════════════════════
# Fixtures — direct builds, one home (lessons L1)
# ═════════════════════════════════════════════════════════════════════════

#: The shipped Sood-style library — every (region, ng) pair, all ``eg=None``.
_SHIPPED_PAIRS = [
    (region, ng_key)
    for region in ("A", "B", "C", "D")
    for ng_key in ("1g", "2g", "4g")
]


def _asymmetric_fissile_2g() -> Mixture:
    """The CATCHER fixture: asymmetric ``SigS``/``Sig2``, P1 stack, χ ∦ νΣf.

    Asymmetry makes a group-axis transpose observable (Mode 6); the P1
    moment makes an ℓ-truncation observable (M1.3's requirement); the
    non-parallel (χ, νΣf) pair makes a dyad factor swap observable
    (G1.7's requirement).
    """
    return anisotropic_mixture(
        [1.1, 2.3],
        [[0.38, 0.10], [0.05, 0.90]],
        [[0.02, 0.01], [0.00, 0.04]],
        sig_f=[0.02, 0.31], chi=[0.95, 0.05],
        sig_l=[0.004, 0.011],
        sig_2=[[0.0, 0.03], [0.01, 0.0]],
    )


def _second_fissile_2g() -> Mixture:
    """A second, different material — the per-``mid`` iteration witness."""
    return anisotropic_mixture(
        [0.9, 2.6],
        [[0.22, 0.03], [0.12, 1.10]],
        [[0.05, 0.02], [0.01, 0.03]],
        sig_f=[0.05, 0.12], chi=[0.80, 0.20],
        sig_l=[0.002, 0.007],
        sig_2=[[0.0, 0.02], [0.005, 0.0]],
    )


def _symmetric_2g() -> Mixture:
    """The declared NON-CATCHER control: symmetric ``SigS``/``Sig2``.

    A symmetric matrix annihilates every transpose/orientation mutation
    (``vv`` Mode 6) — this row exists to DOCUMENT that blindness beside
    the asymmetric catcher, never to be counted as coverage of it.
    """
    return anisotropic_mixture(
        [1.0, 1.5],
        [[0.30, 0.08], [0.08, 0.45]],
        [[0.02, 0.01], [0.01, 0.04]],
        sig_l=[0.003, 0.009],
        sig_2=[[0.010, 0.006], [0.006, 0.020]],
    )


def _two_material_carrier():
    """A fresh two-material carrier + its mixtures, keyed by ``mid``.

    ``from_materials`` retains the extra entry (the single cell uses id
    0), and ``MaterialXSField``'s dense caches cover EVERY materials-dict
    entry — so one degenerate carrier exercises the per-``mid`` accessor
    surface for both ids.
    """
    mixtures = {0: _asymmetric_fissile_2g(), 1: _second_fissile_2g()}
    return MaterialMesh.from_materials(mixtures).material_xs_field(), mixtures


# ═════════════════════════════════════════════════════════════════════════
# G1.1 — ng is the Mixture's, for every shipped pair
# ═════════════════════════════════════════════════════════════════════════

@pytest.mark.parametrize(("region", "ng_key"), _SHIPPED_PAIRS)
def test_kernel_ng_matches_the_mixture(region, ng_key):
    """**G1.1** — all three kernels report the source mixture's ``ng``.

    Reddened by reading a wrong length (``len(chi)`` where a matrix edge
    is meant, an off-by-one in the stack) anywhere in a kernel's shape
    plumbing.
    """
    mixture = get_mixture(region, ng_key)
    assert ScatteringKernel.from_mixture(mixture).ng == mixture.ng
    assert N2NKernel.from_mixture(mixture).ng == mixture.ng
    assert FissionKernel.from_mixture(mixture).ng == mixture.ng


# ═════════════════════════════════════════════════════════════════════════
# G1.2 — truncation is exact, the identity at L, and a REFUSAL beyond
# ═════════════════════════════════════════════════════════════════════════

def test_truncated_returns_exactly_the_requested_stack():
    """**G1.2** — ``truncated(order)`` is exact (positive + negative, vv#11).

    The negative leg is the load-bearing one: moments beyond the stored
    order do not exist, and silent zero-padding would misreport the
    material's anisotropy content as a measured zero (a fabricated
    datum — the campaign's O1 tell).
    """
    kernel = ScatteringKernel.from_mixture(_asymmetric_fissile_2g())
    assert kernel.order == 1  # the fixture ships a P1 stack

    p0_only = kernel.truncated(0)
    assert p0_only.order == 0
    assert len(p0_only.moments) == 1
    np.testing.assert_array_equal(p0_only.p0, kernel.p0)

    identity = kernel.truncated(kernel.order)
    assert identity.order == kernel.order
    for ours, theirs in zip(identity.moments, kernel.moments, strict=True):
        np.testing.assert_array_equal(ours, theirs)

    with pytest.raises(ValueError, match="not invented"):
        kernel.truncated(kernel.order + 1)
    with pytest.raises(ValueError, match="not invented"):
        kernel.truncated(-1)


# ═════════════════════════════════════════════════════════════════════════
# G1.3 — kernel ≡ carrier cache, 0 ULP (bit-identity, never view-identity)
# ═════════════════════════════════════════════════════════════════════════

@pytest.mark.parametrize(
    "build_mixture",
    [_asymmetric_fissile_2g, _symmetric_2g],
    ids=["asymmetric-catcher", "symmetric-declared-non-catcher"],
)
def test_kernel_equals_carrier_cache_bit_identical(build_mixture):
    """**G1.3** — every kernel datum equals the carrier's cache at 0 ULP.

    The claim's licence (verification plan §2(h).5): the two sides are
    **independently assembled from the Mixture's sparse data** — the
    kernel densifies ``SigS``/``Sig2``/reads ``chi``/``SigP`` in
    ``from_mixture``; the carrier densifies through its own
    ``_build_dense_caches`` path — so their agreement is a real claim
    whatever the later absorb-vs-delegate fate of the caches, and
    ``array_equal`` (never ``allclose``) is right because both sides are
    ``todense()`` of one sparse source at reduction depth 0.

    The ``symmetric-declared-non-catcher`` row is exactly that (Mode 6):
    a transpose mutation is invisible on it, and it is shipped to
    document the blindness beside the catcher, not as coverage — which
    is why its carrier holds ONLY the symmetric material ([M] first
    battery run: pairing it with the asymmetric second material made the
    "blind" row catch, falsifying its own declaration). The asymmetric
    row keeps the two-material carrier for the per-``mid`` iteration.
    """
    if build_mixture is _symmetric_2g:
        mixtures = {0: _symmetric_2g()}
    else:
        mixtures = {0: build_mixture(), 1: _second_fissile_2g()}
    mat_xs = MaterialMesh.from_materials(mixtures).material_xs_field()

    for mid, mixture in mixtures.items():
        scattering = ScatteringKernel.from_mixture(mixture)
        cached = mat_xs.sig_s_legendre(mid)
        assert len(scattering.moments) == len(cached)
        for l, cache_matrix in enumerate(cached):
            np.testing.assert_array_equal(scattering.moments[l], cache_matrix)

        np.testing.assert_array_equal(
            N2NKernel.from_mixture(mixture).matrix, mat_xs.n2n_matrix(mid)
        )

        fission = FissionKernel.from_mixture(mixture)
        np.testing.assert_array_equal(fission.chi, mat_xs.chi_per_material(mid))
        np.testing.assert_array_equal(
            fission.nu_sig_f, mat_xs.fission_production_per_material(mid)
        )


# ═════════════════════════════════════════════════════════════════════════
# G1.4 — the ℓ=0 slice IS what the iso pair consumes
# ═════════════════════════════════════════════════════════════════════════

def test_p0_and_emission_are_what_the_iso_pair_consumes():
    """**G1.4** — ``p0``/``emission_matrix`` against the storage-side oracle.

    ``dense_per_material`` is the deliberately independent
    transpose-copy view (its own docstring's stated job), so these rows
    pin the kernel's slice AND orientation against an oracle that shares
    no realization with it: ``p0`` is the ``[g_from, g_to]`` transpose
    of the iso operator matrix, and ``emission_matrix()`` carries the
    multiplicity the raw ``matrix`` deliberately does not (M1.5's
    two-tier separation from G1.3).

    The CS4a-constructible half of the done-when's "slice-consistency
    crosscheck" (F8): the ANGULAR ℓ=0-block agreement is CS4c's, when S
    itself re-points at the kernel.
    """
    mat_xs, mixtures = _two_material_carrier()
    iso_scatter = IsotropicScattering(mat_xs).dense_per_material()
    iso_n2n = IsotropicN2N(mat_xs).dense_per_material()

    for mid, mixture in mixtures.items():
        np.testing.assert_array_equal(
            ScatteringKernel.from_mixture(mixture).p0, iso_scatter[mid].T
        )
        np.testing.assert_array_equal(
            N2NKernel.from_mixture(mixture).emission_matrix(), iso_n2n[mid]
        )


# ═════════════════════════════════════════════════════════════════════════
# G1.5 — the kernel does NOT alias the carrier cache (the F4 hazard, closed)
# ═════════════════════════════════════════════════════════════════════════

def test_kernel_does_not_alias_the_carrier_cache():
    """**G1.5** — non-aliasing, read-only, and carrier-mutation isolation.

    The hazard this pins (CS4a fact F4): the shipped
    ``sig_s_legendre`` returns the production cache object itself,
    writable — a consumer mutation reaches the loss matrix. Four legs,
    each a separate arm of the guard (vv#17 granularity —
    M1.6 reddens the identity+isolation legs, M1.7 the flags leg ALONE):

    1. the kernel array is a different object than the cache;
    2. every kernel buffer is write-protected;
    3. writing through the kernel RAISES;
    4. mutating the CARRIER's cache leaves the kernel datum unchanged.
    """
    mat_xs, mixtures = _two_material_carrier()
    kernel = ScatteringKernel.from_mixture(mixtures[0])
    cache = mat_xs.sig_s_legendre(0)

    # Identity legs FIRST (an aliasing mutation fails HERE; a
    # copies-but-unfrozen mutation fails on the flags loop below — the
    # vv#17 per-arm distinction, made legible by the leg order).
    assert kernel.moments[0] is not cache[0]

    # The MIXTURE-side identity legs — the aliasing surface
    # ``from_mixture`` actually has: its matrix data goes through
    # ``todense()`` (fresh by construction), but ``chi``/``SigP`` are
    # the mixture's own dense arrays, handed in by reference. The
    # carrier legs pin the F4 direction; these pin the provenance
    # direction.
    mixture = mixtures[0]
    fission = FissionKernel.from_mixture(mixture)
    assert fission.nu_sig_f is not mixture.SigP
    assert fission.chi is not mixture.chi

    for array in (
        *kernel.moments,
        N2NKernel.from_mixture(mixtures[0]).matrix,
        fission.chi,
        fission.nu_sig_f,
    ):
        assert array.flags.writeable is False

    with pytest.raises(ValueError):
        kernel.moments[0][0, 0] = 999.0

    before = kernel.p0[0, 0]
    cache[0][0, 0] += 999.0  # the measured F4 reach — must NOT propagate
    assert kernel.p0[0, 0] == before

    sig_p_before = fission.nu_sig_f[0]
    mixture.SigP[0] += 999.0  # mutate the SOURCE — must not reach the kernel
    assert fission.nu_sig_f[0] == sig_p_before


# ═════════════════════════════════════════════════════════════════════════
# G1.6 — the hoisted energy-arm rule (one home: EnergyAxis.from_materials)
# ═════════════════════════════════════════════════════════════════════════

_EDGES_2G = np.array([2.0e7, 1.0e5, 1.0e-3])  # descending, fast-first


@pytest.mark.parametrize(("region", "ng_key"), _SHIPPED_PAIRS)
def test_energy_arm_all_absent_is_synthetic(region, ng_key):
    """**G1.6** (all-absent, 12 witnesses) — ``eg=None`` everywhere ⟹ synthetic."""
    mixture = get_mixture(region, ng_key)
    assert mixture.eg is None  # the row's own precondition, asserted
    axis = EnergyAxis.from_materials([mixture])
    assert axis == EnergyAxis.synthetic(mixture.ng)
    assert axis.edges is None


def test_energy_arm_content_equal_edges_is_from_grid():
    """**G1.6** (the ONE from_grid witness) — content-equal edges ⟹ the grid axis.

    Equality is content (edges BYTES), never object identity: the two
    materials carry separately-constructed edge arrays.
    """
    first = dataclasses.replace(_asymmetric_fissile_2g(), eg=_EDGES_2G)
    second = dataclasses.replace(_second_fissile_2g(), eg=_EDGES_2G.copy())
    axis = EnergyAxis.from_materials([first, second])
    assert axis == EnergyAxis.from_grid(EnergyGrid(_EDGES_2G))
    assert axis.edges is not None
    assert axis != EnergyAxis.synthetic(2)


def test_energy_arm_differing_or_mixed_edges_are_synthetic():
    """**G1.6** (differing / mixed / empty) — anything short of unanimity is synthetic.

    Reddened by flipping the arm (M1.8) or by weakening unanimity to
    majority/first-wins.
    """
    with_edges = dataclasses.replace(_asymmetric_fissile_2g(), eg=_EDGES_2G)
    other_edges = dataclasses.replace(
        _second_fissile_2g(), eg=np.array([1.0e7, 5.0e4, 1.0e-3])
    )
    absent = _symmetric_2g()

    differing = EnergyAxis.from_materials([with_edges, other_edges])
    assert differing == EnergyAxis.synthetic(2)

    mixed = EnergyAxis.from_materials([with_edges, absent])
    assert mixed == EnergyAxis.synthetic(2)

    with pytest.raises(ValueError, match="at least one material"):
        EnergyAxis.from_materials([])


def test_bulk_space_energy_arm_reads_only_reachable_materials():
    """**G1.6** (the call-site denominator) — ``bulk_space`` passes REACHABLE
    materials, not the whole dict.

    The leak principle's witness: a retained SPECTATOR entry with
    ``eg=None`` must not flip the axis identity of the problem the
    single cell (material 0) actually poses. Reddened by widening the
    ``bulk_space`` call site's denominator from reachable to all
    materials.
    """
    with_edges = dataclasses.replace(_asymmetric_fissile_2g(), eg=_EDGES_2G)
    spectator = _symmetric_2g()
    assert spectator.eg is None
    carrier = MaterialMesh.from_materials({0: with_edges, 1: spectator})
    axes = carrier.bulk_space.axes
    assert axes is not None
    assert axes[0] == EnergyAxis.from_grid(EnergyGrid(_EDGES_2G))


# ═════════════════════════════════════════════════════════════════════════
# G1.7 — the fission dyad: direction pinned, transpose = the factor swap
# ═════════════════════════════════════════════════════════════════════════

def test_fission_dyad_direction_and_transpose_theorem():
    """**G1.7** — ``dyad()`` is |χ⟩⟨νΣf| and its transpose is the factor swap.

    Both rows read the RAW mixture factors on the right-hand side, so a
    kernel that swaps the factors in the forward direction reds here on
    the χ ∦ νΣf fixture (and only there — parallel factors annihilate
    the swap, which is why the fixture asserts its own non-parallelism).

    Scope (verification plan §2(h).4): this is a THEOREM about the dyad,
    gated at the kernel tier with no operator in the room. The
    production fission transpose (``FissionOperator``'s fresh
    ``TensorProductOperator``) is UNCHANGED at CS4a — it rebinds at
    CS4c.
    """
    mixture = _asymmetric_fissile_2g()
    chi = np.asarray(mixture.chi, dtype=float)
    nu_sig_f = np.asarray(mixture.SigP, dtype=float)
    cosine = float(chi @ nu_sig_f / (np.linalg.norm(chi) * np.linalg.norm(nu_sig_f)))
    assert cosine < 0.999  # χ ∦ νΣf — the swap is observable

    kernel = FissionKernel.from_mixture(mixture)
    np.testing.assert_array_equal(kernel.dyad(), np.outer(chi, nu_sig_f))
    np.testing.assert_array_equal(kernel.dyad().T, np.outer(nu_sig_f, chi))


# ═════════════════════════════════════════════════════════════════════════
# G1.8 — frozen, and replace re-validates
# ═════════════════════════════════════════════════════════════════════════

def test_kernels_are_frozen_and_replace_revalidates():
    """**G1.8** — immutability + the route-through-replace invariant.

    The design's ruling, stated per the gate table: ``ng`` (and
    ``order``) are DERIVED properties, not fields — so
    ``replace(kernel, ng=...)`` raises ``TypeError`` by construction
    (there is no field to set), and every legal ``replace`` re-runs
    ``__post_init__``, re-establishing shape coherence, the read-only
    buffers, and the χ law.
    """
    scattering = ScatteringKernel.from_mixture(_asymmetric_fissile_2g())
    n2n = N2NKernel.from_mixture(_asymmetric_fissile_2g())
    fission = FissionKernel.from_mixture(_asymmetric_fissile_2g())

    for kernel, field_name in (
        (scattering, "moments"), (n2n, "matrix"), (fission, "chi"),
    ):
        with pytest.raises(dataclasses.FrozenInstanceError):
            setattr(kernel, field_name, None)

    with pytest.raises(TypeError):
        dataclasses.replace(scattering, ng=3)  # type: ignore[call-arg]

    with pytest.raises(ValueError, match="same square"):
        dataclasses.replace(
            scattering,
            moments=(scattering.moments[0], np.zeros((3, 3))),
        )

    # The χ law re-fires on replace: a producing kernel refuses a
    # non-simplex spectrum (the one law, enforce_emission_spectrum).
    with pytest.raises(ValueError):
        dataclasses.replace(fission, chi=np.array([0.4, 0.3]))


# ═════════════════════════════════════════════════════════════════════════
# G1.9 — the C8 import fence, as an import-list assertion
# ═════════════════════════════════════════════════════════════════════════

def test_module_imports_nothing_from_scattering_or_frames():
    """**G1.9** — ``kernels.py`` never imports scattering/frame machinery.

    The fence's direction is the doctrine: ``ScatteringOperator``
    re-points at :class:`ScatteringKernel` (CS4c), never the reverse.
    The walk covers EVERY import statement in the module — late
    function-body imports included — so the fence cannot be tunneled
    under.
    """
    source_path = inspect.getsourcefile(kernels_module)
    assert source_path is not None
    source = Path(source_path).read_text()
    imported: list[str] = []
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Import):
            imported.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            imported.append(node.module or "")
    offenders = [
        name for name in imported
        if "scattering" in name or "frame" in name
    ]
    assert offenders == [], (
        f"orpheus/transport/kernels.py imports {offenders} — the C8 fence "
        f"forbids the kernel module from reaching scattering/frame "
        f"machinery (the dependency points the other way)."
    )
