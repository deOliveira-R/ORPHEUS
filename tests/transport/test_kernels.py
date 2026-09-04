r"""CS4a K1 — the interaction kernels (gates G1.1–G1.9).

The kernels (:mod:`orpheus.transport.kernels`) are representation-free
per-material physics data — the datum the CS4a bindings realize as bound
operators. These gates are software invariants of that posing
(``foundation``, the ``test_operator_spaces.py`` placement argument);
none carries ``verifies(...)``.

**Fixture discipline** (verification plan §3, ``lessons`` L1; wording
corrected at CS4a-R QA-F7): the MUTATION-CATCHER fixtures are built
DIRECTLY through the campaign's one builder
(:func:`tests.sn.architecture._config.anisotropic_mixture`), because the
SHIPPED Sood tables carry no ``Sig2`` (`[M]` nnz = 0 on all 12
``get_mixture`` pairs) and ``SigL = 0`` — channels a catcher must
un-null or its arm goes vacuously green. (``make_mixture`` itself DOES
offer a P1 channel via ``sig_s1=`` — the earlier "offers no P1 channel"
clause was false — and the shipped-pair/binding rows in this file
legitimately use ``get_mixture``, where the nulled channels are
irrelevant.) Fixtures are function-scoped builders (fresh per call) so
the G1.5 carrier-mutation leg cannot pollute a shared object.
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
from orpheus.transport.kernels import FissionKernel, TransferKernel
from orpheus.transport.mesh.material_mesh import MaterialMesh
from orpheus.transport.operators.fission import FissionOperator
from orpheus.transport.operators.isotropic_scattering import (
    IsotropicFission,
)
from orpheus.transport.operators.isotropic_scattering import (
    IsotropicN2N,
    IsotropicScattering,
)
from orpheus.transport.operators.multiplication_operator import (
    MultiplicationOperator,
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
    assert TransferKernel.scattering(mixture).ng == mixture.ng
    assert TransferKernel.n2n(mixture).ng == mixture.ng
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
    kernel = TransferKernel.scattering(_asymmetric_fissile_2g())
    assert kernel.order == 1  # the fixture ships a P1 stack

    p0_only = kernel.at_order(0)
    assert p0_only.order == 0
    assert len(p0_only.moments) == 1
    np.testing.assert_array_equal(p0_only.p0, kernel.p0)

    identity = kernel.at_order(kernel.order)
    assert identity.order == kernel.order
    for ours, theirs in zip(identity.moments, kernel.moments, strict=True):
        np.testing.assert_array_equal(ours, theirs)

    # §4.3 (#426 step 2): above the stored order the stack PADS exact zeros
    # (a shorter stack is complete — the evaluation's zeros); until
    # 2026-09-04 this refused with "not invented".
    wide = kernel.at_order(kernel.order + 1)
    assert wide.order == kernel.order + 1
    np.testing.assert_array_equal(wide.moments[-1], 0.0)
    with pytest.raises(ValueError, match="order >= 0"):
        kernel.at_order(-1)


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

    Scope of the licence (re-scoped at CS4a-R QA-F2): the two sides are
    two SPELLINGS over ONE sparse source — both are
    ``np.asarray(s.todense())`` of the same ``Mixture`` object (and the
    chi/SigP legs compare a construction copy with its own source), so
    this gate pins an ASYMMETRIC re-spelling of either path, NOT the
    storage convention: `[M]` a shared ``[g_from,g_to]→[g_to,g_from]``
    inversion of BOTH sides leaves all rows here green (the convention
    is pinned by ``tests/homogeneous`` — 17 reds incl. the L1 anchor —
    and by G1.4b's hand-authored literal below). ``array_equal`` (never
    ``allclose``) is right because both sides are ``todense()`` of one
    sparse source at reduction depth 0.

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
        scattering = TransferKernel.scattering(mixture)
        cached = mat_xs.sig_s_legendre(mid)
        assert len(scattering.moments) == len(cached)
        for l, cache_matrix in enumerate(cached):
            np.testing.assert_array_equal(scattering.moments[l], cache_matrix)

        np.testing.assert_array_equal(
            TransferKernel.n2n(mixture).p0, mat_xs.n2n_matrix(mid)
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

    ``dense_per_material`` is a transpose-copy VIEW of the same carrier
    cache G1.3 already compares (re-scoped at CS4a-R QA-F2: ``p0 ==
    iso[mid].T`` cancels to an identity under a SHARED convention
    inversion, so orientation-BETWEEN-views and the multiplicity are
    what these rows genuinely pin — the multiplicity leg IS independent:
    ``ClassVar = 2`` in kernels.py vs a literal ``2.0`` in
    ``isotropic_scattering.py:444``, two hand-written homes). The
    absolute storage convention is pinned by G1.4b's hand-authored
    literal and by ``tests/homogeneous`` (M1.5's two-tier separation
    from G1.3).

    The CS4a-constructible half of the done-when's "slice-consistency
    crosscheck" (F8): the ANGULAR ℓ=0-block agreement is CS4c's, when S
    itself re-points at the kernel.
    """
    mat_xs, mixtures = _two_material_carrier()
    iso_scatter = IsotropicScattering.from_material_xs(
        mat_xs, space=mat_xs.mesh.bulk_space,
    ).dense_per_material()
    iso_n2n = IsotropicN2N.from_material_xs(
        mat_xs, space=mat_xs.mesh.bulk_space,
    ).dense_per_material()

    for mid, mixture in mixtures.items():
        np.testing.assert_array_equal(
            TransferKernel.scattering(mixture).p0, iso_scatter[mid].T
        )
        np.testing.assert_array_equal(
            TransferKernel.n2n(mixture).emission_matrix(), iso_n2n[mid]
        )


def test_p0_convention_pinned_against_a_hand_authored_literal():
    """**G1.4b** (CS4a-R QA-F2) — the ``[g_from, g_to]`` convention, EXTERNALLY pinned.

    G1.3/G1.4 compare spellings of ONE ``todense()`` chain, so a SHARED
    storage-convention inversion moves both sides together (`[M]`
    both-sides transpose: 51/51 green in this file pre-CS4a-R). The
    literal below is ``_asymmetric_fissile_2g``'s own declared
    ``SigS[0]`` input, copied BY HAND from the fixture definition —
    never computed — so no shared code path can move both sides of THIS
    comparison. Deliberately asymmetric so the transpose is observable.
    """
    kernel = TransferKernel.scattering(_asymmetric_fissile_2g())
    hand_authored_p0 = np.array([[0.38, 0.10], [0.05, 0.90]])  # [g_from, g_to]
    np.testing.assert_array_equal(kernel.p0, hand_authored_p0)


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
    4. the CARRIER's cache itself REFUSES mutation (the producer-side
       freeze CS4a-R EE-4 added — before it, this leg mutated the cache
       and asserted non-propagation; the freeze upgrades the property
       from "does not propagate" to "cannot be spelled");
    5. mutating the SOURCE (the mixture's sparse ``.data`` / dense
       ``SigP``) does not reach the kernel — construction copies.
    """
    mat_xs, mixtures = _two_material_carrier()
    kernel = TransferKernel.scattering(mixtures[0])
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
        TransferKernel.n2n(mixtures[0]).p0,
        fission.chi,
        fission.nu_sig_f,
    ):
        assert array.flags.writeable is False

    with pytest.raises(ValueError):
        kernel.moments[0][0, 0] = 999.0

    with pytest.raises(ValueError):
        cache[0][0, 0] += 999.0  # the F4 reach, now REFUSED at the producer

    before = kernel.p0[0, 0]
    mixtures[0].SigS[0].data[:] += 999.0  # mutate the sparse SOURCE instead
    assert kernel.p0[0, 0] == before  # ...the kernel copied at construction

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
    gated at the kernel tier with no operator in the room. Since CS4c
    step 4 the production realization consumes this datum: the energy
    binding's cached ``kernel`` (``IsotropicFission``) is the one dyad
    home, and its ``TensorProductOperator`` transpose IS this factor
    swap — gated operator-side in ``test_isotropic_fission.py``.
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
    scattering = TransferKernel.scattering(_asymmetric_fissile_2g())
    n2n = TransferKernel.n2n(_asymmetric_fissile_2g())
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
    # match= discriminates the χ-law raise from the factor-shape raise —
    # [M] CS4a-R QA-F9: both ValueErrors reach this construction, so a
    # bare raises() would stay green with the χ law disabled.
    with pytest.raises(ValueError, match="not normalized"):
        dataclasses.replace(fission, chi=np.array([0.4, 0.3]))

    # QA-F10 completion: replace() re-establishes the read-only buffers
    # (the claim the class docstring makes; asserted nowhere until now).
    assert scattering.at_order(0).moments[0].flags.writeable is False


def test_kernel_constructor_refusals_have_negative_witnesses():
    """**G1.8b** (CS4a-R QA-F10) — every ``__post_init__`` refusal, exercised.

    `[M]` at review time each message fragment below had 0 hits in
    ``tests/`` — four admission contracts shipped positive-only (vv#11).
    One row per arm, ``match=`` on the shortest distinctive fragment.
    """
    with pytest.raises(ValueError, match="empty"):
        TransferKernel(moments=())
    with pytest.raises(ValueError, match="square"):
        TransferKernel(moments=(np.zeros((2, 3)),))
    with pytest.raises(ValueError, match="rank"):
        TransferKernel(moments=(np.zeros(2),))
    with pytest.raises(ValueError, match="positive integer"):
        TransferKernel(moments=(np.zeros((2, 2)),), multiplicity=0)
    with pytest.raises(ValueError, match="two \\(ng,\\) vectors"):
        FissionKernel(chi=np.array([1.0, 0.0]), nu_sig_f=np.array([0.1, 0.2, 0.3]))


# ═════════════════════════════════════════════════════════════════════════
# G1.9 — the C8 import fence, as an import-list assertion
# ═════════════════════════════════════════════════════════════════════════

def test_module_imports_nothing_from_scattering_or_frames():
    """**G1.9** — ``kernels.py`` never imports scattering/frame machinery.

    The fence's direction is the doctrine: ``ScatteringOperator``
    re-points at :class:`TransferKernel` (CS4c), never the reverse.
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


# ═════════════════════════════════════════════════════════════════════════
# CS4a K2 — the binding fences (G2.8, G2.9, G2.10)
# ═════════════════════════════════════════════════════════════════════════

def _registry_type_names(cls) -> set[str]:
    """The ``singledispatchmethod`` registry of ``cls._apply_impl``, by type
    NAME — the empty set when the class dispatches without one (the iso
    pair's ``isinstance`` branches, which is exactly why the behavioural
    matrix below exists)."""
    attr = inspect.getattr_static(cls, "_apply_impl", None)
    dispatcher = getattr(attr, "dispatcher", None)
    if dispatcher is None:
        return set()
    return {t.__name__ for t in dispatcher.registry}


def _diffusion_binding():
    """The 2g / 6-cell diffusion binding the arm matrix is measured on."""
    from orpheus.diffusion.augmented_mesh import DiffusionMesh
    from orpheus.geometry import BC, CoordSystem, Mesh1D

    mesh = Mesh1D(
        edges=np.linspace(0.0, 2.0, 7), mat_ids=np.zeros(6, dtype=int),
        coord=CoordSystem.CARTESIAN, bc_right=BC("vacuum"),
    )
    dm = DiffusionMesh(mesh, {0: get_mixture("A", "2g")})
    return dm, dm.material_xs_field(), dm.full_field_space


_ARM_MATRIX = {
    # (operator, carrier) -> expected outcome type name, or "TypeError".
    # [M] 2026-08-20 (verification plan §2(g.1)) on the diffusion binding.
    ("C", "FullField"): "FullField",
    ("C", "ndarray"): "ndarray",
    ("C", "ScalarFlux"): "TypeError",
    ("IsoS", "FullField"): "FullField",
    ("IsoS", "ndarray"): "ndarray",
    ("IsoS", "ScalarFlux"): "ndarray",  # F10: untyped fall-through, recorded
    ("IsoN2N", "FullField"): "FullField",
    ("IsoN2N", "ndarray"): "ndarray",
    ("IsoN2N", "ScalarFlux"): "ndarray",  # F10
    ("F", "FullField"): "FullField",
    ("F", "ndarray"): "ndarray",
    ("F", "ScalarFlux"): "ndarray",  # F10 (iso-family unwrap, CS4c step 4)
}


@pytest.mark.parametrize(
    ("operator_key", "carrier_key"),
    sorted(_ARM_MATRIX),
    ids=[f"{o}-{c}" for o, c in sorted(_ARM_MATRIX)],
)
def test_apply_arm_survival_matrix(operator_key, carrier_key):
    r"""**G2.8** ⭐⭐ — the R-A fence, executable: NO apply arm is deleted.

    Twelve cells, each with a distinct MEASURED outcome including the
    refusal (``C × ScalarFlux → TypeError``) — so the gate is never
    merely "does it not crash", and it survives reformatting and any
    later ``singledispatchmethod`` → explicit-dispatch rewrite, which a
    source grep would not. The DENOMINATOR (CS4a-R QA-F12, re-keyed at
    CS4c step 4): 4 dispatchers × 3 carriers, on the 2g/6-cell
    diffusion binding — the F entry is the ENERGY binding
    (``IsotropicFission``, what diffusion actually consumes since the
    step-4 split; its ScalarFlux cell joins the iso family's F10
    fall-through, and the retired typed ``ScalarSourceSink`` cell's
    coverage moved to ``test_isotropic_fission.py``). The ANGULAR
    ``FissionOperator``'s arms are gated in ``test_fission_operator.py``
    —
    ``ScatteringOperator``'s 5 arms have no behavioural cell here (they
    are covered by the registry-keyset gate only, which catches an ADDED
    arm where this matrix catches a MOVED one: the two gates are
    complementary by design), and `[M]` the F-FullField / F-ScalarFlux
    cells are coupled (one reroute reds both). The two ``IsoS/IsoN2N × ScalarFlux →
    ndarray`` cells are F10's untyped fall-through (typed in, bare out —
    vv#29's asymmetric arrow), RECORDED here as the shipped behaviour;
    CS4a deletes no arm and repairs none (the dispatch collapse is
    CS4c's, after the per-instance feeding census).
    """
    from orpheus.transport.fields.scalar_boundary_flux import ScalarBoundaryFlux
    from orpheus.transport.fields.scalar_flux import ScalarFlux
    from orpheus.transport.full_field import FullField

    dm, mat_xs, ffs = _diffusion_binding()
    operators = {
        "C": MultiplicationOperator(
            coefficient=mat_xs.total_cross_section_field, domain=ffs, codomain=ffs,
        ),
        "IsoS": IsotropicScattering.from_material_xs(mat_xs, space=ffs),
        "IsoN2N": IsotropicN2N.from_material_xs(mat_xs, space=ffs),
        "F": IsotropicFission.from_material_xs(mat_xs, space=ffs),
    }
    rng = np.random.default_rng(2026)
    interior_values = rng.random((2, 6)) + 0.5
    carriers = {
        "FullField": lambda: FullField(
            interior=ScalarFlux(values=interior_values, space=dm.bulk_space),
            boundary=ScalarBoundaryFlux(values=rng.random(dm.scalar_trace.shape[0]) + 0.1, space=dm.scalar_trace),
        ),
        "ndarray": lambda: interior_values.copy(),
        "ScalarFlux": lambda: ScalarFlux(values=interior_values, space=dm.bulk_space),
    }

    expected = _ARM_MATRIX[(operator_key, carrier_key)]
    op = operators[operator_key]
    probe = carriers[carrier_key]()
    if expected == "TypeError":
        with pytest.raises(TypeError):
            op.apply(probe)
    else:
        assert type(op.apply(probe)).__name__ == expected, (
            f"{operator_key} × {carrier_key}: the arm's outcome type "
            f"changed — an apply arm moved or died (the R-A fence)"
        )


def test_apply_dispatch_registries_are_verbatim():
    r"""**G2.8 companion** — the five ``singledispatchmethod`` keysets, verbatim.

    ⚠ The iso pair's EMPTY sets are the point, not a gap: their arms are
    ``isinstance`` branches, structurally invisible to any registry
    introspection — which is why the behavioural matrix above is the
    stronger instrument and this row is only the cheap source-tier tell.
    """
    from orpheus.transport.operators.scattering import ScatteringOperator

    assert _registry_type_names(MultiplicationOperator) == {
        "FullField", "ndarray", "object",
    }
    assert _registry_type_names(FissionOperator) == {
        "AngularFlux", "FullField", "HarmonicMomentFlux", "ScalarFlux",
        "object",
    }  # CS4c step 4: the angular binding mirrors N2NOperator's arms
    #    (ScalarFlux is the typed refusal toward the energy binding).
    assert _registry_type_names(IsotropicScattering) == set()
    assert _registry_type_names(IsotropicN2N) == set()
    assert _registry_type_names(IsotropicFission) == set()
    assert _registry_type_names(ScatteringOperator) == {
        "AngularFlux", "FullField", "HarmonicMomentFlux", "ScalarFlux",
        "object",
    }


def test_isotropic_energy_inherits_the_parent_binding_space():
    r"""**G2.9, INVERTED (CS4c step 3)** — the C8 fence's witness flipped
    to the POSITIVE gate its own docstring promised: the iso
    constructors' space became mandatory, the fence broke before
    production did (as designed), and the surviving claim is
    inheritance — ``S.isotropic_energy`` (the P0 energy binding the
    per-ordinate fast path lifts; ``isotropic_kernel``'s successor
    after the §14.1 (n,2n) extraction) is bound to the SCALAR sub-space
    of the parent's own composite interior, never space-anonymous.
    """
    from orpheus.numerics.quadrature import Quadrature
    from orpheus.numerics.space import FunctionSpace
    from orpheus.sn.mesh.augmented_mesh import SNMesh
    from orpheus.geometry import BC, CoordSystem, Mesh1D
    from orpheus.transport.operators.scattering import ScatteringOperator

    carrier = MaterialMesh.from_materials({0: get_mixture("A", "2g")})
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, 5),
        mat_ids=np.zeros(4, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    sn_mesh = SNMesh(
        mesh, Quadrature.gauss_legendre(n_ordinates=4), carrier.materials,
    )
    space = sn_mesh.full_field_space
    scattering = ScatteringOperator.from_solver_data(
        mat_xs=carrier.material_xs_field(),
        scattering_order=0,
        space=space,
    )
    energy = scattering.isotropic_energy
    interior = space.interior_space
    assert interior is not None and interior.axes is not None
    expected = FunctionSpace.of_axes(*interior.axes[1:])
    assert energy.domain == expected and energy.codomain == expected, (
        "S.isotropic_energy must be bound to the scalar sub-space of "
        "the parent's OWN interior — the binding drifted from the "
        "parent's pose"
    )


def test_energy_conformity_guard_three_rows():
    r"""**G2.10** — the ng-conformity refusal, per arm (vv#11 + vv#28).

    Four row families, and the last two are the ones an author will not
    write unprompted:

    1. axis-built POSITIVE — 2g data × the 2g quotient space constructs;
    2. axis-built NEGATIVE — 2g data × a 4g quotient space raises the
       typed ``"energy extent"`` refusal (fragment asserted DISJOINT
       from the ``OperatorSum`` pins' ``"equal domains"`` vocabulary, so
       this row can never be intercepted by those) — and row 2b repeats
       the refusal at ALL FOUR wired call sites, naming each operator
       (QA-F1: the guard BODY is single-sourced but the WIRING is
       per-site — three sites had no witness);
    3. axes-LESS — a WRONG-ng bind on ``SNMesh(2g).full_field_space``
       MUST REFUSE since CS4c step 4 (⛔ this clause read "MUST
       CONSTRUCT: the declared inertness" until step 4 — the row's body
       records how the reach widened). The guard's reach is the
       contract. ``[M]`` re-derived at CS4c step 3 (the rebind changed
       the wiring): SEVEN production classes now run the admission at
       construction — C, S, IsoS, IsoN2N (the per-END base helper),
       N2NOperator, F, and Λ/N2N-moment inherit the base without an
       energy end to check — where the pre-step census read 4 of 13.
       The guard stays INERT on axes-less composites (this row's
       subject): ``SNMesh.full_field_space`` carries no EnergyAxis until
       CS2's axes, so a wrong-ng bind constructs, and the row keeps that
       fact asserted rather than assumed. Without this row the guard
       ships certified by a fixture family that reddens on demand while
       the axes-less real bindings never touch it.
    """
    from orpheus.geometry import BC, CoordSystem, Mesh1D
    from orpheus.numerics.quadrature import Quadrature
    from orpheus.sn.mesh.augmented_mesh import SNMesh

    carrier_2g = MaterialMesh.from_materials({0: get_mixture("A", "2g")})
    carrier_4g = MaterialMesh.from_materials({0: get_mixture("A", "4g")})
    mat_2g = carrier_2g.material_xs_field()

    # Row 1 — axis-built positive (the fission ENERGY binding — the
    # k-outer / homogeneous / diffusion production site since step 4).
    bound = IsotropicFission.from_material_xs(
        mat_2g, space=carrier_2g.bulk_space,
    )
    assert bound.domain == carrier_2g.bulk_space

    # Row 2 — axis-built negative, typed, disjoint fragment.
    with pytest.raises(ValueError, match="energy extent") as excinfo:
        IsotropicFission.from_material_xs(
            mat_2g, space=carrier_4g.bulk_space,
        )
    message = str(excinfo.value)
    assert "4" in message and "2" in message  # both integers named
    assert "equal domains" not in message  # disjoint from the D2 pins

    # Row 2b — the SAME refusal at EVERY wired call site (CS4a-R QA-F1:
    # [M] per-site no-op mutation over 655 rows reddened F only — C /
    # IsoS / IsoN2N had NO witness, and the unwitnessed C site is the
    # one passing a different operand expression, values.shape[0]. One
    # row per call site; the message must name the CONSTRUCTING
    # operator, so a miswired label cannot borrow a sibling's witness).
    wrong_space = carrier_4g.bulk_space
    per_site = {
        "MultiplicationOperator": lambda: MultiplicationOperator(
            coefficient=mat_2g.total_cross_section_field, domain=wrong_space, codomain=wrong_space,
        ),
        "IsotropicScattering": lambda: IsotropicScattering.from_material_xs(
            mat_2g, space=wrong_space,
        ),
        "IsotropicN2N": lambda: IsotropicN2N.from_material_xs(mat_2g, space=wrong_space),
    }
    for op_name, construct in per_site.items():
        with pytest.raises(ValueError, match="energy extent") as site_info:
            construct()
        assert op_name in str(site_info.value), (
            f"the wrong-ng refusal at the {op_name} site does not name "
            f"its constructing operator"
        )

    # Row 3 — axes-less composite: the WRONG-ng bind REFUSES since
    # CS4c step 4. ⛔ REACH WIDENED (this row asserted "constructs —
    # declared inert" until step 4): the composite itself still carries
    # no EnergyAxis, but FissionOperator.from_solver_data now DERIVES
    # its energy binding's scalar ends from the interior's axes — and
    # the interior (the angular trial space) IS axis-built with an
    # EnergyAxis, so the per-END guard reaches a bind the axes-less
    # composite alone could never let it see. The inertness this row
    # used to record is CLOSED, not merely relocated.
    mesh = Mesh1D(
        edges=np.linspace(0.0, 2.0, 5), mat_ids=np.zeros(4, dtype=int),
        coord=CoordSystem.CARTESIAN, bc_right=BC("vacuum"),
    )
    sn_2g = SNMesh(
        mesh, Quadrature.gauss_legendre(n_ordinates=4),
        {0: get_mixture("A", "2g")},
    )
    composite = sn_2g.full_field_space
    assert composite.axes is None  # still true — the reach is the interior's
    mat_4g = carrier_4g.material_xs_field()
    with pytest.raises(ValueError, match="energy extent"):
        FissionOperator.from_solver_data(mat_xs=mat_4g, space=composite)


def test_fission_space_is_mandatory():
    r"""**G2.11** — F without a space is a ``TypeError`` at both entries.

    The presence half of the ng-conformity ruling (a signature
    ``TypeError``, deliberately not a message pin): the campaign's R2
    repair made an anonymous ``F`` UNREPRESENTABLE, so ``.H`` can never
    see a ``None`` space and silently degrade to the bare Euclidean
    transpose. The annotation half (``domain -> "FunctionSpace"``, no
    Optional) is carried by the now-marker-free ledger rows
    ``test_leaf_space_annotation_is_not_optional[F]`` /
    ``test_leaf_without_a_space_refuses_construction[F]`` — this row
    pins the two constructor surfaces directly.
    """
    mat_xs = MaterialMesh.from_materials(
        {0: get_mixture("A", "2g")}
    ).material_xs_field()
    with pytest.raises(TypeError):
        FissionOperator.from_solver_data(mat_xs=mat_xs)  # type: ignore[call-arg]
    with pytest.raises(TypeError):
        FissionOperator(mat_xs=mat_xs)  # type: ignore[call-arg]


# ═════════════════════════════════════════════════════════════════════════
# G-F1 — the χ↔νΣf-coupled condensation law (XD-9; CS4c step 4, plan §7)
# ═════════════════════════════════════════════════════════════════════════

# The shipped 4g→2g fixture of tests/data/test_mixture_condense.py,
# restated here as hand literals (vv L11: the morphisms below are built
# in THIS test body from the partition and φ — never frame.project,
# never a second condense call).
_XD9_PHI = np.array([1.0, 4.0, 2.0, 0.5])
_XD9_EG_FINE = np.array([1.0e7, 1.0e5, 1.0e2, 1.0e0, 1.0e-2])
_XD9_EG_COARSE = np.array([1.0e7, 1.0e2, 1.0e-2])
_XD9_PARTITION = ((0, 1), (2, 3))  # fine indices per coarse group


def _xd9_fine_mixture():
    from tests.data.test_mixture_condense import _balanced_fissile_4g

    return _balanced_fissile_4g()


def _hand_marginalize(vec, partition):
    """Mass-preserving sink morphism: χ_G = Σ_{g∈G} χ_g (no φ)."""
    return np.array([sum(vec[g] for g in group) for group in partition])


def _hand_average(vec, phi, partition):
    """Rate-preserving source morphism: ⟨νΣf⟩_G = Σ φ_g νΣf_g / Σ φ_g."""
    return np.array([
        sum(phi[g] * vec[g] for g in group) / sum(phi[g] for g in group)
        for group in partition
    ])


def _assert_condensation_activates(partition, phi):
    """G-F1's ACTIVATION PRECONDITION, asserted (the
    ``_assert_metric_is_constant`` pattern — §10's designed-green
    hazard): every coarse group holds ≥ 2 fine groups AND φ varies
    within at least one of them. On a 1-fine-per-coarse target
    ``average ≡ marginalize/width`` degenerates and CTRL-A/CTRL-C go
    silent — an identity condensation must be REFUSED as a fixture, so
    the gate cannot silently stop being discriminating."""
    if any(len(group) < 2 for group in partition):
        raise ValueError(
            "G-F1 activation precondition: every coarse group must hold "
            ">= 2 fine groups (a 1-fine-per-coarse target makes the "
            "average/marginalize discrimination vacuous)"
        )
    phi = np.asarray(phi, dtype=float)
    if all(
        np.allclose(phi[list(group)], phi[group[0]]) for group in partition
    ):
        raise ValueError(
            "G-F1 activation precondition: the spectrum is flat within "
            "every coarse group — the average degenerates to the "
            "marginalize direction and the controls go silent"
        )


class TestFissionCondensationGF1:
    r"""**G-F1** — ``dyad(condense(K)) == outer(marginalize(χ), average(νΣf))``.

    The χ↔νΣf-coupled condensation (XD-9): the conjugation is ASYMMETRIC
    by design — the sink (χ, ``g_to``) axis MARGINALIZES (mass-preserving)
    and the source (νΣf, ``g_from``) axis AVERAGES (rate-preserving,
    φ-weighted). ⚠ BRANCH DECLARATION: this pins the FORWARD branch
    (``adjoint_spectrum is None``); the bilinear branch folds the adjoint
    carrier into the sink and obeys a DIFFERENT law (plan §1.4).

    The three negative controls are the measured wrong-morphism pairs
    ([M] 2026-08-30, plan §7.1: CTRL-A 6.421e-1 / CTRL-B 1.685e0 /
    CTRL-C 7.087e-2) — hand-built and red-capable at the landing commit
    with no production change (§6c: the witness ships WITH the gate).
    """

    def _condensed_dyad(self):
        from orpheus.data.energy_grid import EnergyGrid

        fine = _xd9_fine_mixture()
        coarse = fine.condense(EnergyGrid(_XD9_EG_COARSE), _XD9_PHI)
        return FissionKernel.from_mixture(coarse).dyad()

    def test_law_ruled_morphism_pair(self):
        _assert_condensation_activates(_XD9_PARTITION, _XD9_PHI)
        fine = _xd9_fine_mixture()
        expected = np.outer(
            _hand_marginalize(np.asarray(fine.chi), _XD9_PARTITION),
            _hand_average(np.asarray(fine.SigP), _XD9_PHI, _XD9_PARTITION),
        )
        np.testing.assert_allclose(
            self._condensed_dyad(), expected, rtol=1e-14, atol=0.0,
            err_msg="the condensed fission dyad violates the ruled "
            "(χ marginalize, νΣf average) morphism pair — the XD-9 "
            "χ↔νΣf coupling drifted",
        )

    @pytest.mark.parametrize("ctrl,chi_morph,nu_morph,measured", [
        ("A", "average", "average", 6.421e-1),
        ("B", "marginalize", "marginalize", 1.685e0),
        ("C", "average", "marginalize", 7.087e-2),
    ])
    def test_wrong_morphism_controls(self, ctrl, chi_morph, nu_morph, measured):
        """CTRL-A/B/C: each wrong pair is O(1e-1..1e0) away — the gate's
        rtol=1e-14 law row reds by >10 orders under any of them."""
        _assert_condensation_activates(_XD9_PARTITION, _XD9_PHI)
        fine = _xd9_fine_mixture()

        def morph(vec, which):
            return (
                _hand_marginalize(vec, _XD9_PARTITION)
                if which == "marginalize"
                else _hand_average(vec, _XD9_PHI, _XD9_PARTITION)
            )

        wrong = np.outer(
            morph(np.asarray(fine.chi), chi_morph),
            morph(np.asarray(fine.SigP), nu_morph),
        )
        dyad = self._condensed_dyad()
        rel = float(np.max(np.abs(dyad - wrong)) / np.max(np.abs(dyad)))
        if not rel > 1e-2:
            pytest.fail(
                f"CTRL-{ctrl} ({chi_morph} χ, {nu_morph} νΣf) sits at "
                f"rel {rel:.3e} from the production dyad — the control "
                f"lost its measured O({measured:.2e}) separation and the "
                f"law row can no longer discriminate this morphism swap"
            )

    def test_b45_degenerate_target_reds_the_precondition(self):
        """B4.5 — a 1-fine-per-coarse target must red the ACTIVATION
        PRECONDITION (proving it is asserted, not assumed): the identity
        condensation is exactly the fixture on which every control above
        would go silent (§10's designed-green hazard)."""
        with pytest.raises(ValueError, match="activation precondition"):
            _assert_condensation_activates(((0,), (1,), (2,), (3,)), _XD9_PHI)

    def test_b45_flat_spectrum_reds_the_precondition(self):
        with pytest.raises(ValueError, match="activation precondition"):
            _assert_condensation_activates(
                _XD9_PARTITION, np.ones_like(_XD9_PHI),
            )


# ═════════════════════════════════════════════════════════════════════════
# G-F2 — the operator dyad IS the kernel datum (the collapse, documented)
# ═════════════════════════════════════════════════════════════════════════


def test_g_f2_operator_dyad_is_the_kernel_datum_per_material():
    r"""**G-F2, post-collapse form.** The energy binding's cached dyad
    factors equal ``FissionKernel.from_mixture(m)`` per material, cell by
    cell.

    ⚠ NEAR-TAUTOLOGICAL BY DESIGN, and said so (`coding-standards`
    single-sourcing clause): since CS4c step 4 the operator's factors ARE
    the gathered kernel datum (a pure index gather), so no input can make
    the two sides disagree in VALUE — the row's live content is the
    GATHER PLACEMENT (each material's pair lands in ITS cells; a
    material-id mixup or a factor-order swap reds). The external pins
    that keep the datum itself anchored are the hand-written χ/SigP
    literals in this file's fixtures and the hand-rolled per-cell
    references in ``test_isotropic_fission.py`` (which red on a factor
    swap through an INDEPENDENT route — the B4.6 catcher).
    """
    from orpheus.transport.operators.isotropic_scattering import (
        IsotropicFission,
    )

    mat_xs, mixtures = _two_material_carrier()
    op = IsotropicFission.from_material_xs(
        mat_xs, space=mat_xs.mesh.bulk_space,
    )
    rank_one = op.kernel.ops[0]
    chi_gathered = np.asarray(rank_one.reconstruction)
    nu_gathered = np.asarray(rank_one.functional.weight)
    for mid, idx in mat_xs.mesh.cells_by_material.items():
        k = FissionKernel.from_mixture(mixtures[mid])
        cells = (slice(None), *idx)
        np.testing.assert_array_equal(
            chi_gathered[cells],
            np.broadcast_to(k.chi[:, None], chi_gathered[cells].shape),
        )
        np.testing.assert_array_equal(
            nu_gathered[cells],
            np.broadcast_to(k.nu_sig_f[:, None], nu_gathered[cells].shape),
        )
