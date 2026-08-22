r"""A declared boundary law leaves the OPERATOR linear — and leaves it alone.

Phase **P5** of `.claude/plans/archive/affine_boundary_source_channel.md`. The affine
boundary law is

.. math::

    \gamma_-\psi = L\,\gamma_+\psi + q

and the whole campaign is the claim that the realizer tier carries **only**
:math:`L`, while :math:`q` travels the composite-source channel. Prescribed
inflow is the :math:`L = 0` case, so its realized operator is the *zero
morphism* — bound to :math:`\Gamma_+(f) \to \Gamma_-(f)`, stamped
``BlockRole.BOUNDARY``, and contributing nothing to ``B``.

⛔ **P3 fixed a LIVE bug here, not a hypothetical.** Before ``8d552395`` a
declared ``PrescribedInflow`` realized to an *affine* operator carrying its own
source, and :meth:`SNBoundaryOperator._face_laws` collects every face's law with
**no** ``block_role`` filter — so it reached ``B``. Measured then: ``|B(0)| =
2.5``; the inflow was delivered TWICE through SI (``γ₋ψ = 5.0``, ratio exactly
``2.000000``); and on Krylov an affine ``A(x) = A_lin(x) − c`` breaks GMRES's
Arnoldi relation ``A V_k = V_{k+1} H_k``, so the solve RAISED
``ConvergenceCertificateError``. These rows are the operator-tier statement of
that fix (ERR-075).

⭐ What "linear" MEANS on this carrier
======================================

Since campaign 1 CS3 (2026-08-19) the carrier IS a vector space: flux lives in
V (the affine/torsor field algebra was overturned by user ruling), so the
textbook laws are directly spellable and this module states them directly:

============================  ================================================
law                           states
============================  ================================================
``B(0) = 0``                  no constant term. The zero flux is V's origin.
``B(c·ψ) = c·B(ψ)``           homogeneity (scalar scaling — how ``ψ/k`` is
                              legal).
``B(ψ₁ + ψ₂) =``              additivity, DIRECT — and applied to a
``B(ψ₁) + B(ψ₂)``             difference, ``B(ψ₁ − ψ₂) = B(ψ₁) − B(ψ₂)``,
                              which is the **#331 closure**: the operator acts
                              on any element of V, increments included.
============================  ================================================

**History (kept per coding-standards — it records why the old spelling
existed).** Until 2026-08-19 the domain was an *affine space* over a distinct
displacement space (the #208 torsor carve): ``psi + psi`` raised ``TypeError``
("no origin"), so additivity could only be spelled as base-point independence —
``B(ψ₁ ⊕ σ) − B(ψ₂ ⊕ σ) == B(ψ₁) − B(ψ₂)`` with σ a displacement — and #331
recorded that ``B`` refused a displacement argument outright. Writing the
direct row was the first draft's mistake then; it is the correct row now.

⚠ **Attributability** (``vv`` anti-pattern #18): all three rows now red for an
affine ``B(x) = Lx + q`` — including additivity, where ``q`` survives once on
the left and twice on the right. The retired base-point spelling could NOT
catch the affine regression (``B(x₁)−B(x₂) = L(x₁−x₂)`` whatever ``q``), which
its docstring said honestly; the direct spelling is strictly sharper, so the
CS3 re-spell UPGRADED this battery's affine coverage by two rows.

`[M]` **The battery** (``vv`` #17 + #18). Three mutations of the realizer's
prescribed arm, over this module + the reciprocity module + the P3 trace gate
(**69 rows**, baseline all green):

==============================  =======  ===========================================
mutation                        reds     what it establishes
==============================  =======  ===========================================
``control`` — ``B.apply``       27       the positive control. It reddens every
scaled by ``|ψ_trace|_inf``              homogeneity row AND both base-point rows,
(NONLINEAR)                              so a later "0 caught" verdict means the
                                         harness died, not that the gates are weak.
``affine`` — the P3 regression  19       ⭐ ``B(0)``, ``A(0)``, all 8 homogeneity
(``+ q`` inside the operator)            rows, §0's two rows, §3, the 3 NEW
                                         reciprocity cases, and the P3 trace gate.
                                         The 9 pre-existing reflective reciprocity
                                         cases stay GREEN — the new cases are
                                         attributable to the declaration.
``identity`` — ``L := I`` for   6        every linearity row stays **green**, as
prescribed (perfectly LINEAR)            :func:`test_B_vanishes_at_zero` predicts.
                                         Caught by §0's two rows, by §3, and by the
                                         P3 trace gate.
==============================  =======  ===========================================

The ``affine`` column's *misses* were the informative part: neither
base-point-independence row reddened, exactly as their docstrings said. (That
battery ran against the pre-CS3 spelling; the direct additivity rows that
replaced the base-point rows at the cone carve DO red on the affine mutation —
re-measured at the carve, see their docstrings.)

⚠ The battery harness itself has now failed **three** times in this campaign,
every time reading as "your tests are fine" — twice by mis-parsing pytest's
output, once by patching a module object pytest never collected. The third was a
reviewer's independent run of THIS battery: an unquoted ``$VAR`` of three paths
became one bogus argument, pytest collected nothing and exited 0 in 0.01 s.
Re-run it with ``--color=no -rf`` and **read the collected count** before
believing any verdict it gives.

⛔ **Why this fixture and not the §4.6 MMS one.** The P4 MMS slab declares
``PrescribedInflow`` on **both** faces, so after P3 ``B`` is *identically* the
zero morphism there — ``|B(x)|_inf = 0.0`` for random ``x``. Every row above
posed on it would hold because both sides are structurally zero, and **no input
could make it red**: the ``vv`` Mode-8 tautological-companion class. This module
therefore pairs the declared prescribed face with a **reflective** partner,
whose live specular law is what makes ``B`` a non-trivial operator at all
(:func:`test_the_fixture_activates_B_on_the_reflective_face_alone` measures both
halves of that precondition, and every row below rests on it).

The leaf-tier statement of §3's theorem is already landed
(``test_operator_block_role.py::test_prescribed_inflow_realizes_the_same_object_vacuum_does``
pins type + both spaces off the realizer). What it cannot say is what the
ASSEMBLED ``B`` and the full matvec DO — and the assembled tier is precisely
where the P3 defect lived, because the missing stamp was never read.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.geometry.boundary import ConstantInflowSource, PrescribedInflow
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.coupled_system import build_within_group_system
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.operators.boundary import SNBoundaryOperator
from orpheus.sn.solver import SNSolver
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.source_sinks import AngularBoundarySourceSink
from orpheus.transport.timed_full_field import TimedFullField

pytestmark = [pytest.mark.foundation]

_VALUE = 2.5
_N_ORD = 8

#: The scalars homogeneity is asserted at. ``-2.0`` is a power of two (exact in
#: IEEE, so it isolates a *structural* break from rounding); ``3.7`` is generic;
#: ``1e3`` / ``1e-3`` span six decades so a defect scaling with the INPUT rather
#: than with ``c`` cannot hide inside one amplitude.
_SCALARS = (3.7, -2.0, 1e3, 1e-3)

#: Relative slack for the FULL matvec's rows, measured. ``B`` alone is exact
#: (``0.0`` at every scalar), but ``(L+C) − S`` accumulates the sweep's
#: floating-point non-associativity: `[M]` the worst relative defect over
#: ``_SCALARS`` and over the base-point-independence row is ``3.7e-16`` — 1.6
#: machine epsilons, stable across six decades of ``c``. ``1e-14`` leaves 27×
#: headroom on that while keeping **twelve orders** of margin against the
#: regression it must see: `[M]` under the ``affine`` mutation this row's own
#: message reads ``relative 4.274e-02`` at ``c = 3.7`` (and 5.9e-2 … 5.9e+1
#: across the other scalars).
#:
#: ⚠ This comment said ``≈ 1.6e-1`` until 2026-08-06. That number is not the
#: regression's reading at all — ``1.603e-01`` is the fixture's
#: boundary/interior magnitude ratio (``6.844 / 42.685``), i.e. a value read off
#: the wrong line of a probe. No verdict changes (12.6 orders of margin either
#: way), which is exactly why it survived: a plausible constant in the one
#: module whose thesis is that plausible constants must be measured.
_MATVEC_RTOL = 1e-14

#: ULP budget for ``B``'s base-point-independence row. `[M]` 4 nulp for one
#: displacement, 8 for a second, on the boundary block; the interior is exact.
#: The defect is NOT an operator defect: the row re-bases both arguments, and
#: ``fl(a+σ) − fl(b+σ) ≠ a − b`` for the same reason any FP subtraction is
#: inexact. 32 leaves 4× headroom over the worst reading while staying far
#: below anything a real linearity break could produce.
_INCREMENT_NULP = 32


# ═══════════════════════════════════════════════════════════════════════
# The fixture — declared on the GEOMETRY, which is the user path
# ═══════════════════════════════════════════════════════════════════════


#: The declaration under test, spelled once.
_PRESCRIBED = PrescribedInflow(source=ConstantInflowSource(value=_VALUE))


def _slab(xmin=_PRESCRIBED, xmax=None) -> SNMesh:
    r"""Het 2G slab, GL-8. Default: ``PrescribedInflow`` at ``xmin``, reflective
    at ``xmax`` — the fixture every row below uses.

    ONE construction site for every mesh in this module, so a row can vary the
    declaration and nothing else. The declaration rides the **geometry**, not a
    poked-at mesh: that is the channel landed at ``985497b5`` and the only route
    that survives ``_as_sn_mesh``'s rebuild inside a public solver.

    Heterogeneous 2G with ``c ≈ 0.90–0.96`` so ``S`` is genuinely live in the
    full-matvec rows (``placeholder_materials`` has ``SigS ≡ 0``, which would
    make ``− S`` the zero operator and quietly drop a third of the algebra
    under test).
    """
    from orpheus.derivations.common.xs_library import get_mixture

    geom = StructuredGeometry(
        geometry="SLB",
        regions=(Region(mat_id=0, outer_thickness_cm=1.0),
                 Region(mat_id=1, outer_thickness_cm=2.0)),
        bcs=(xmin, BC.reflective if xmax is None else xmax),
    )
    mesh = Mesh1D.from_geometry(
        geom, region_meshes=(RegionMesh(n_cells=6), RegionMesh(n_cells=6)),
    )
    return SNMesh(
        mesh, Quadrature.gauss_legendre(n_ordinates=_N_ORD),
        {0: get_mixture("A", "2g"), 1: get_mixture("D", "2g")},
    )


def _zero_flux(sn: SNMesh) -> TimedFullField:
    return TimedFullField.zeros(
        interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn,
    )


def _random_flux(sn: SNMesh, seed: int) -> TimedFullField:
    """A random state on BOTH blocks — a trace-only or bulk-only probe would
    leave half of ``B``'s domain untouched."""
    rng = np.random.default_rng(seed)
    shape = _zero_flux(sn)
    return TimedFullField(
        interior=AngularFlux.from_mesh(
            rng.uniform(0.5, 2.0, size=shape.interior.values.shape), sn,
        ),
        boundary=AngularBoundaryFlux(
            values=rng.uniform(0.5, 2.0, size=shape.boundary.values.shape),
            space=sn.angular_trace,
        ),
        _history=(), history_depth=2,
    )


def _blocks(field):
    """The ``(interior, boundary)`` value arrays of a composite or coupled state."""
    inner = field.systems[0] if hasattr(field, "systems") else field
    return (
        np.asarray(inner.interior.values), np.asarray(inner.boundary.values),
    )


def _linf(field) -> float:
    return max(float(np.max(np.abs(a))) for a in _blocks(field))


def _assert_vanishes(field, what: str) -> None:
    """Every block is exactly zero.

    Asserted against the literal ``0.0`` rather than against a zero *state*:
    ``B`` maps a flux to a RATE DENSITY, so a state-vs-state comparison would
    silently compare an ``AngularSourceSink`` with an ``AngularFlux`` through
    their raw arrays and read green on a role confusion.
    """
    for name, values in zip(("interior", "boundary"), _blocks(field)):
        np.testing.assert_array_equal(
            values, 0.0, err_msg=f"{what} — the {name} block is non-zero",
        )


def _assert_identical(got, want, what: str) -> None:
    """Bit-identity on both blocks, block-named in the message.

    ``assert_array_equal`` and not ``allclose``: ``B``'s action is a gather /
    permutation / scatter with no arithmetic that could round, and `[M]` every
    defect measured here is exactly ``0.0``. A tolerance would be slack this
    claim has not earned — and slack is what lets an affine term back in.
    """
    for name, (g, w) in zip(("interior", "boundary"),
                            zip(_blocks(got), _blocks(want))):
        np.testing.assert_array_equal(
            g, w, err_msg=f"{what} — the {name} block is not bit-identical",
        )


def _assert_within_nulp(got, want, what: str) -> None:
    """``|got − want| ≤ _INCREMENT_NULP`` ULP elementwise, reading reported.

    The nulp is computed here rather than delegated to
    ``np.testing.assert_array_almost_equal_nulp``, which takes no ``err_msg``:
    a caller wanting a message has to append one after the call, which silently
    builds a tuple and discards it. (That trap cost a debugging cycle earlier in
    this campaign; spelling the comparison out makes it unspellable here.)
    """
    for name, (g, w) in zip(("interior", "boundary"),
                            zip(_blocks(got), _blocks(want))):
        spacing = np.spacing(np.maximum(np.abs(g), np.abs(w)))
        reading = float(np.max(
            np.abs(g - w) / np.where(spacing == 0.0, np.spacing(1.0), spacing)
        ))
        assert reading <= _INCREMENT_NULP, (
            f"{what} — the {name} block differs by {reading:.1f} ULP against a "
            f"budget of {_INCREMENT_NULP} (absolute "
            f"{float(np.max(np.abs(g - w))):.3e} at scale "
            f"{float(np.max(np.abs(w))):.6f}). A few ULP is the re-basing's own "
            f"rounding; a large reading means the increment genuinely depends "
            f"on the base point, i.e. the operator is not affine."
        )


def _assert_close_relative(got, want, scale: float, what: str) -> None:
    """``|got − want|_inf ≤ _MATVEC_RTOL · scale``, with the reading reported.

    Normalised against a norm rather than asserted elementwise in ULP: the
    matvec's output is a DIFFERENCE of streaming and collision terms, so its
    small entries carry the cancellation of large ones and an elementwise ULP
    budget there measures the cancellation, not the linearity. `[M]` on this
    fixture the elementwise nulp reading peaks at **213** for homogeneity and
    **256** for the base-point row, while the norm-relative readings are
    ``3.6e-16`` and ``3.7e-16`` — the second pair is the honest statement of the
    same fact.
    """
    defect = max(
        float(np.max(np.abs(g - w)))
        for g, w in zip(_blocks(got), _blocks(want))
    )
    assert defect <= _MATVEC_RTOL * scale, (
        f"{what}: |defect|_inf = {defect:.6e} against a budget of "
        f"{_MATVEC_RTOL:.0e}·{scale:.6e} = {_MATVEC_RTOL * scale:.6e} "
        f"(relative {defect / scale:.3e}). An AFFINE B reads ~1e-1 relative "
        f"here, so a reading in that range means the source is back in the "
        f"operator; a reading just over budget is a rounding question."
    )


def _full_matvec(sn: SNMesh):
    """``A = (L + C) − S − B`` — the object the drivers actually consume.

    Built through :func:`build_within_group_system`, the single construction
    site, so this module cannot drift from the posed system by re-spelling the
    sum locally. On a seedless slab the loss grid is the 1×1 ``[[A_AA]]``, i.e.
    exactly ``(L+C) − S − B_a``.
    """
    solver = SNSolver(sn, inner_solver="source_iteration", scattering_order=0)
    system = build_within_group_system(sn, solver.mat_xs, scattering_order=0)
    return system.loss, system.space


def _coupled(space, field: TimedFullField):
    """Lift a composite state into the 1-system coupled vector the grid takes."""
    return replace(space.zeros(), systems=(field,))


# ═══════════════════════════════════════════════════════════════════════
# §0 — the precondition every row below rests on
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.catches("ERR-075")
def test_an_all_prescribed_mesh_makes_B_the_zero_morphism() -> None:
    r"""⭐ This module's own PREMISE, committed instead of merely asserted in prose.

    The docstring rejects the §4.6 MMS fixture because ``B`` is identically zero
    there, which would make every row below a Mode-8 tautology. That reasoning
    is load-bearing — it is why this module builds a different mesh — and until
    now it lived only as a sentence. `[M]` on the all-prescribed het 2G slab:
    both blocks exactly ``0.0``, ``array_equal(boundary, 0) is True``, and
    ``B.is_adjointable is True``.

    ⚠ The failure this closes is a *joint* one, which is why neither module
    catches it alone: if a regression re-armed a prescribed face, the P4 MMS
    would silently regain its double delivery **and** this module's stated
    rationale would become false — while this module, which deliberately avoids
    that fixture, stayed green throughout. The premise has to be checked on the
    fixture it is a premise ABOUT.

    The adjointability leg rides along because it is free here and states the
    P3 widening in its purest form: a mesh whose every face is prescribed used
    to have no transpose at all.
    """
    sn = _slab(xmin=_PRESCRIBED, xmax=_PRESCRIBED)
    B = SNBoundaryOperator(sn)
    _assert_vanishes(
        B.apply(_random_flux(sn, 7)),
        "an all-prescribed mesh's B is NOT the zero morphism, so the §4.6 MMS "
        "fixture has a live boundary operator again — the P4 module's rows are "
        "no longer measuring what they claim, and THIS module's reason for "
        "using a different fixture has become false",
    )
    assert B.is_adjointable, (
        "an all-prescribed B is not adjointable — the P3 capability widening "
        "has regressed at its purest point (every face is the zero morphism, "
        "which transposes to itself under every metric)"
    )


def test_the_fixture_activates_B_on_the_reflective_face_alone() -> None:
    r"""⭐ The activation guard, decomposed PER FACE so it is attributable.

    A bare ``|B(x)|_inf > 0`` would establish that ``B`` is not the zero
    operator, and nothing more — in particular it would be satisfied by the
    reflective face while the prescribed face did something arbitrary. The two
    halves must be measured separately:

    * ``xmin`` (prescribed) — the WHOLE face slot is exactly zero. This is the
      P3 claim at the assembled tier: the realized law is the zero morphism, so
      the face contributes nothing to ``B``, not merely something small.
    * ``xmax`` (reflective) — the inflow rows are non-zero, which is what makes
      every row below non-vacuous.

    `[M]` on this fixture at seed 7: ``|B(x)|`` on ``xmin`` = ``0.0`` exactly,
    on ``xmax``'s inflow rows = ``1.8239310798528774``.
    """
    sn = _slab()
    out = SNBoundaryOperator(sn).apply(_random_flux(sn, 7))

    np.testing.assert_array_equal(
        np.asarray(out.boundary.face_view("xmin")), 0.0,
        err_msg=(
            "the PRESCRIBED face contributes to B. Its realized law must be "
            "the zero morphism — a non-zero reading is the affine operator "
            "back in the boundary block (ERR-075)."
        ),
    )
    reflective_rows = sn.angular_trace.inflow_indices_for_face("xmax")
    live = float(np.max(np.abs(
        np.asarray(out.boundary.face_view("xmax"))[reflective_rows]
    )))
    assert live > 1e-3, (
        f"the REFLECTIVE face contributes |B(x)| = {live:.3e} ≈ 0, so B is "
        f"the zero operator on this fixture and every row below is a "
        f"tautology (vv Mode 8). Fix the fixture, not the rows."
    )


# ═══════════════════════════════════════════════════════════════════════
# §1 — B is LINEAR (the leaf the campaign changed)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.catches("ERR-075")
def test_B_vanishes_at_zero() -> None:
    r"""⭐ ``B(0) = 0`` — the one row that reads the affine term directly.

    An affine operator is ``x ↦ Lx + q``; at ``x = 0`` it returns exactly the
    constant it should not be carrying. `[M]` pre-P3 on a prescribed face:
    ``|B(0)|_inf = 2.5``, i.e. the declared source itself. Post-P3: ``0.0``.

    ⚠ Necessary but NOT sufficient, and the gap matters: setting
    ``L := IdentityOperator`` for prescribed inflow is perfectly linear, keeps
    ``B(0) = 0``, and this row stays green. The row that reddens for THAT
    mutation is
    ``tests/sn/solve/test_declared_inflow_reaches_the_rhs.py::test_the_declared_boundary_law_holds_on_the_answer``
    (``γ₋ψ|_f == q_f`` on the converged answer). The two are complementary and
    neither replaces the other.
    """
    sn = _slab()
    _assert_vanishes(
        SNBoundaryOperator(sn).apply(_zero_flux(sn)),
        "B(0) ≠ 0, so B is AFFINE: the boundary source is inside the operator "
        "instead of in q_∂ (ERR-075)",
    )


@pytest.mark.parametrize("c", _SCALARS)
def test_B_is_homogeneous(c: float) -> None:
    r"""``B(c·ψ) = c·B(ψ)``, bit-identically. `[M]` defect ``0.0`` at every ``c``.

    Reds for an affine ``B``: ``L(cψ) + q`` differs from ``c(Lψ + q)`` by
    ``(1−c)·q``, which is non-zero at every ``c ≠ 1`` in ``_SCALARS``.
    """
    sn = _slab()
    B, psi = SNBoundaryOperator(sn), _random_flux(sn, 7)
    _assert_identical(
        B.apply(c * psi), c * B.apply(psi), f"B is not homogeneous at c = {c}",
    )


def test_B_is_additive_and_acts_on_differences() -> None:
    r"""``B(ψ₁ + ψ₂) = B(ψ₁) + B(ψ₂)`` and ``B(ψ₁ − ψ₂) = B(ψ₁) − B(ψ₂)`` —
    the direct additivity law, spellable since the CS3 cone carve.

    The second identity is the **#331 closure**: a difference of iterates is
    an ordinary element of V and ``B`` acts on it (until 2026-08-19 ``B``
    refused a displacement argument, and this module could only spell the
    base-point-independence surrogate — see the module docstring's history).

    Unlike the retired surrogate, this row DOES catch the affine regression:
    for ``B(x) = Lx + q`` the left side carries ``q`` once and the right side
    twice.

    Activation: the fixture's reflective face makes ``B`` non-trivial
    (:func:`test_the_fixture_activates_B_on_the_reflective_face_alone`), and
    the in-test guard below re-asserts it — a linearity row on a zero
    morphism holds structurally and can never red (lessons L40c).
    """
    sn = _slab()
    B = SNBoundaryOperator(sn)
    psi_1, psi_2 = (_random_flux(sn, s) for s in (7, 11))
    B_1, B_2 = B.apply(psi_1), B.apply(psi_2)
    if not _linf(B_1) > 0.0:
        raise AssertionError("B is the zero morphism on this fixture — "
                             "the additivity rows below are vacuous")
    _assert_within_nulp(
        B.apply(psi_1 + psi_2), B_1 + B_2,
        "B(ψ₁+ψ₂) ≠ B(ψ₁)+B(ψ₂) — the direct additivity law failed",
    )
    _assert_within_nulp(
        B.apply(psi_1 - psi_2), B_1 - B_2,
        "B(ψ₁−ψ₂) ≠ B(ψ₁)−B(ψ₂) — the #331 closure failed: B does not act "
        "linearly on an iterate difference",
    )


# ═══════════════════════════════════════════════════════════════════════
# §2 — and the FULL matvec stays linear with the declaration installed
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.catches("ERR-075")
def test_the_full_matvec_vanishes_at_zero() -> None:
    r"""``A(0) = 0`` for ``A = (L+C) − S − B``. `[M]` exactly ``0.0``.

    ⭐ **Why the composite row is not redundant with** :func:`test_B_vanishes_at_zero`.
    The drivers never apply ``B``; they apply the posed loss grid, and an
    ``OperatorSum`` is where a leaf's affineness becomes the *solver's* problem
    — that is the mechanism by which the P3 defect took down Krylov
    (``A V_k = V_{k+1} H_k`` is an identity about a LINEAR ``A``; scipy's
    tracked residual is meaningless without it, and
    ``_certify_within_group_exit`` caught the divergence at
    ``‖Aψ − q‖/‖q‖ = 1.718``). This row is the algebraic precondition for
    GMRES being applicable at all on a mesh with a declared inflow.
    """
    sn = _slab()
    A, space = _full_matvec(sn)
    _assert_vanishes(
        A.apply(_coupled(space, _zero_flux(sn))),
        "A(0) ≠ 0 — the within-group loss operator is AFFINE, which breaks "
        "GMRES's Arnoldi relation and makes its residual meaningless",
    )


@pytest.mark.parametrize("c", _SCALARS)
def test_the_full_matvec_is_homogeneous(c: float) -> None:
    r"""``A(c·ψ) = c·A(ψ)`` to ``_MATVEC_RTOL``. `[M]` worst relative ``3.6e-16``."""
    sn = _slab()
    A, space = _full_matvec(sn)
    psi = _coupled(space, _random_flux(sn, 7))
    A_psi = A.apply(psi)
    _assert_close_relative(
        A.apply(c * psi), c * A_psi, abs(c) * _linf(A_psi),
        f"A is not homogeneous at c = {c}",
    )


def test_the_full_matvec_is_additive_and_acts_on_differences() -> None:
    r"""``A(x₁ + x₂) = A(x₁) + A(x₂)`` and ``A(x₁ − x₂) = A(x₁) − A(x₂)`` for
    the full within-group matvec ``A = (L+C) − S − B`` — direct, typed.

    The composite sibling of :func:`test_B_is_additive_and_acts_on_differences`.
    What it adds over the ``B``-only row is that ``(L+C)``, ``S`` and ``B``
    COMPOSE into the property — the difference leg drives ``S.apply`` and
    ``B.apply`` with an iterate difference through the assembled sum, which is
    exactly the spelling #331 reported as refused (``S.apply: unsupported
    input type AngularDisplacement``; both gains now act on any element of V).
    """
    sn = _slab()
    A, space = _full_matvec(sn)
    psi_1, psi_2 = (_random_flux(sn, s) for s in (7, 11))
    A_1 = A.apply(_coupled(space, psi_1))
    A_2 = A.apply(_coupled(space, psi_2))
    if not _linf(A_1) > 0.0:
        raise AssertionError("A is the zero morphism on this fixture — "
                             "the additivity rows below are vacuous")
    scale = max(_linf(A_1), _linf(A_2))
    _assert_close_relative(
        A.apply(_coupled(space, psi_1 + psi_2)), A_1 + A_2, scale,
        "A(x₁+x₂) ≠ A(x₁)+A(x₂) — the assembled matvec lost additivity",
    )
    _assert_close_relative(
        A.apply(_coupled(space, psi_1 - psi_2)), A_1 - A_2, scale,
        "A(x₁−x₂) ≠ A(x₁)−A(x₂) — the #331 closure failed at the assembled "
        "tier: a gain refused or mishandled an iterate difference",
    )


# ═══════════════════════════════════════════════════════════════════════
# §3 — ⭐⭐ the campaign theorem, at the assembled tier
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.catches("ERR-075")
def test_declaring_prescribed_moves_q_and_leaves_the_operator_untouched() -> None:
    r"""⭐⭐ Vacuum and prescribed inflow differ **only** in :math:`q`.

    The sharpest statement the campaign has, and the one row here that is fully
    attributable to the prescribed declaration. Two meshes identical except for
    the ``xmin`` declaration — ``PrescribedInflow(2.5)`` versus ``BC.vacuum`` —
    must satisfy BOTH halves:

    * **the operators are bit-identical**: ``B`` and the whole loss grid
      ``(L+C) − S − B`` give the same floats on the same input, because both
      laws realize to the SAME bound zero morphism. `[M]` ``|ΔB|_inf = 0.0``
      and ``|ΔA|_inf = 0.0``, ``array_equal = True``, at ``|A(ψ)|_inf = 42.685``.
    * **the sources differ**: ``q_∂`` reads ``2.5`` against vacuum's ``0.0``.

    Neither half is optional. Without the second, a channel that dropped the
    declaration entirely would satisfy the first perfectly — which is exactly
    the pre-P2′ state (a silently-inert declaration), so a one-sided row here
    would certify the campaign's *opening* defect as correct.

    ⛔ ``assert_array_equal``, never a tolerance. The two meshes do not run two
    computations that happen to agree; after P3 they run ONE float program with
    two different ``q``. Any non-zero difference means a second, law-dependent
    path has reappeared in the operator, and a ``2.9e-14`` relative divergence
    would sail through any ``rtol`` scaled off a solver tolerance.

    ⚠ The reflective partner face is load-bearing. On an all-prescribed mesh
    both ``B``'s would be identically zero and this row would hold for a ``B``
    that ignored its faces altogether.

    ⭐ `[M]` This row catches both ``affine`` AND ``identity`` in the module
    docstring's three-mutation battery — including ``L := IdentityOperator``,
    which every linearity row above is blind to by construction. It is the only
    operator-tier row here that sees a prescribed law whose linear factor is
    wrong but *linear*. (It is green under ``control``, which mutates
    ``SNBoundaryOperator.apply`` identically on both meshes and so cancels in a
    mesh-vs-mesh comparison — correct behaviour, not a gap.)
    """
    prescribed, vacuum = _slab(), _slab(xmin=BC.vacuum)

    # Half 1 — the OPERATORS agree, at both tiers.
    psi_p, psi_v = _random_flux(prescribed, 7), _random_flux(vacuum, 7)
    _assert_identical(
        SNBoundaryOperator(prescribed).apply(psi_p),
        SNBoundaryOperator(vacuum).apply(psi_v),
        "B differs between a declared prescribed face and a vacuum face — the "
        "LINEAR factor is reading the source, which is the P3 defect",
    )
    A_p, space_p = _full_matvec(prescribed)
    A_v, space_v = _full_matvec(vacuum)
    _assert_identical(
        A_p.apply(_coupled(space_p, psi_p)), A_v.apply(_coupled(space_v, psi_v)),
        "the within-group loss differs between prescribed and vacuum — the "
        "declaration reached the OPERATOR, not just q",
    )

    # Half 2 — …and the SOURCES do not. Without this the row is vacuous.
    assert AngularBoundarySourceSink.from_mesh_laws(prescribed).linf == _VALUE, (
        "the declared inflow did not reach q_∂, so half 1 above is comparing "
        "two vacuum problems and proves nothing"
    )
    assert AngularBoundarySourceSink.from_mesh_laws(vacuum).linf == 0.0
