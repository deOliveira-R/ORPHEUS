r"""BC reference-image harness — the realized law vs an independent expression.

For each of seven BC × quadrature × face cases, this file realizes the SN
boundary operator and asserts it reproduces a **frozen, independently-derived
reference image** at the per-case tolerance. The reference lives in
``tests/geometry/snapshots/bc_equivalence_<case>.npz`` and is produced by
:mod:`tests.geometry._generate_bc_equivalence_snapshots`.

.. _bc-snapshot-harness-inversion:

What changed on 2026-08-01, and why the claim got stronger
==========================================================

Until this rewrite the artefacts were **recordings of production output**: the
generator called ``SNBoundaryRealizer().realize(...)`` and froze
``op.apply(psi_out)``. The assertion was therefore ``production == a recording
of production`` — a regression LOCK worth exactly what the recorded code was
right, and (mechanically) a thing that breaks whenever production's signature
changes. It broke twice: at **B3.2**, when the SN law narrowed to
:math:`\Gamma_+ \to \Gamma_-`, and again at **B3.4a**.

The interim repair had been to keep the frozen full-face artefacts and
*restrict* them, on the argument that a pre-narrowing artefact is an
independent statement. That argument is only **procedural** independence
(``vv-principles`` §structural independence, L11): it certifies *"the new path
agrees with the old path"*, and this campaign's premise is that the old path
read the wrong half-trace.

The user's ruling was to migrate the schema and **re-anchor every case against
an expression derived from the mathematics** — the affine form
:eq:`affine-bc-form`, the mirror isometry
:math:`\Omega \mapsto \Omega - 2(\Omega\cdot\hat n)\hat n`, and the Lambertian
partial-current balance :math:`J^- = \alpha J^+`. All seven artefacts were
regenerated. The gate now says *"the operator does what the law's equation
says"*, not *"the operator still does what it did"* — a promotion in the
evidence hierarchy, from drift lock to correctness claim.

What the migration measured (2026-08-01, and worth keeping)
-----------------------------------------------------------

The derived reference was compared against the retired pre-B3.2 recordings,
restricted to :math:`\Gamma_-`:

* ``vacuum``, ``specular_x``, ``specular_y`` — **bit-identical**. The
  independent expression reproduces the frozen artefact exactly.
* ``white_xmax``, ``white_xmin``, ``mixed`` — **1 to 2 ULP** (reduction order:
  the reference contracts with :func:`numpy.tensordot`, production runs a
  broadcast-multiply-then-``sum``).
* ``periodic`` — **98 % relative**, by design; see that case below.

So the re-anchoring did not quietly move six of the seven claims. It changed
where they come from.

Why the artefact is still frozen, now that the reference is computable
----------------------------------------------------------------------

Because a reference the harness recomputed would let the generator's
expression and production drift **together**. The file on disk is the barrier:
this module reads ``psi_in`` off disk and **never imports the generator's
``_*_image`` helpers**, so changing an expression moves the gate only through a
regeneration — a reviewable diff of committed artefacts.

The registry (``CASES``) is imported, and only the registry: it is the single
source of truth for which law, quadrature and face each case is about. Before
this migration every case's law and quadrature were spelled twice — once in the
generator, once here — with nothing to keep them in step.

V&V tags
--------

* ``@pytest.mark.foundation`` — a software invariant about operator
  construction, with no theory-page ``:label:`` to point at.
* ``@pytest.mark.regression`` — the reference is a committed frozen artefact.
* No ``@pytest.mark.verifies`` — stacking it on a ``foundation`` test is silent
  level conflation.

**B0.3 RELABEL (2026-07-30), still standing.** This file was once marked
``l1`` on the claim that the realizer path is a structurally-independent
re-implementation. That was false while the baseline was self-generated. It is
now *less* false — the reference IS structurally independent — but the level
marker stays ``foundation``, because what these rows verify is operator
construction, not an equation on a theory page. The honest home for the
correctness content is the expression's derivation, which lives in the
generator's docstring.

Regenerating
------------

::

    python -m tests.geometry._generate_bc_equivalence_snapshots

Regeneration is intentional only when the reference EXPRESSION or the probe
changes, and must carry an inline justification. It is no longer coupled to
production: the generator imports nothing from
:mod:`orpheus.sn.boundary.realizer`.
"""
from __future__ import annotations

import ast
from pathlib import Path

import numpy as np
import pytest

from orpheus.geometry.boundary import (
    AlbedoBoundary,
    SpecularReturn,
)
from orpheus.numerics.operator import LinearOperator
from orpheus.sn.boundary.realizer import SNBoundaryRealizer
from tests.geometry._generate_bc_equivalence_snapshots import (
    BCEquivalenceCase,
    CASES,
    case_by_id,
)
from tests.sn._test_helpers import face_method_space


pytestmark = [pytest.mark.foundation, pytest.mark.regression]


_SNAPSHOT_DIR = Path(__file__).parent / "snapshots"
_EPS = float(np.finfo(np.float64).eps)
_GENERATOR_HINT = (
    "Regenerate with `python -m tests.geometry._generate_bc_equivalence_snapshots`. "
    "Since the 2026-08-01 re-anchoring that generator computes each reference "
    "from the law's EQUATION and imports nothing from the SN realizer, so it "
    "cannot be broken by a production signature change (it was, twice: B3.2 "
    "and B3.4a). Regeneration is legitimate only when the reference expression "
    "or the probe changes — never to make a red go away, which would re-anchor "
    "the gate on the very code it gates."
)


def _load_snapshot(case_id: str) -> np.lib.npyio.NpzFile:
    """Load the frozen reference for ``case_id``; a missing file is a FAILURE.

    B0.3 REPAIR — this helper used to ``pytest.skip`` on a missing snapshot.
    A skip is invisible in a summary line, so that was a way for the widest net
    in the subsystem to be silently disabled by a deleted or renamed baseline
    (``vv-principles`` Mode-8, the SKIP-SWALLOWED class). A missing reference is
    a broken gate, so it reds.
    """
    path = _SNAPSHOT_DIR / f"bc_equivalence_{case_id}.npz"
    if not path.exists():
        pytest.fail(
            f"Reference {path.name} is MISSING — this gate is the widest "
            f"mutation net in the boundary subsystem and a missing baseline "
            f"silently disables it. {_GENERATOR_HINT}",
        )
    return np.load(path)


def _verified_space(case: BCEquivalenceCase, snapshot: np.lib.npyio.NpzFile):
    r"""The live method space, cross-checked against the artefact's Layer-1 metadata.

    Every artefact stores its face and BOTH half-traces, so the classification
    :math:`\Omega\cdot\hat n \lessgtr \mp\epsilon` is itself pinned by a frozen
    independent derivation (the generator writes the index sets from the face
    convention, never from the trace space). This is the row that would catch
    the B3.4a ``> 0.0`` classification twin returning: it mis-admits tangential
    ordinates into :math:`\Gamma_+`, which moves the index SET.

    It matters that this runs before every value comparison: a value gate on
    the wrong half-trace is a Mode-12 blind spot waiting to happen, because
    ``|Γ₊| == |Γ₋|`` on every quadrature in the tree — so shapes agree even
    when the sets do not.

    **Honest scope, MEASURED 2026-08-01.** This row is BLIND to the ``> 0.0``
    twin *on the case set as it stands*: mutating ``TANGENTIAL_EPS`` to ``0.0``
    leaves all seven cases green, because every quadrature here
    (``lebedev(17)``, ``level_symmetric(4)`` / ``(6)``, ``gauss_legendre(8)``)
    carries its tangential cosines as EXACT ``0.0``, which a strict ``> 0.0``
    cut excludes just as the band does. The quadrature that discriminates is
    ``product(2, 4)``, whose tangential cosines are round-off (±5e-17), and it
    is covered by
    :class:`tests.geometry.test_reemission_closure.TestDiffuseAgainstAnIndependentExpression`
    — deliberately, since duplicating it here would add a case and no
    coverage. What this row DOES catch is any change that moves an index set
    at all, which a flipped outward normal does (measured: 7 reds).
    """
    space = face_method_space(
        case.build_quadrature(), face=case.face, faces=case.faces,
    )
    assert str(snapshot["face"]) == case.face
    np.testing.assert_array_equal(
        np.sort(np.asarray(space.inflow_indices)), snapshot["inflow_indices"],
        err_msg=(
            f"{case.case_id}: the live Γ₋ at face {case.face!r} differs from "
            f"the frozen one. The artefact's index sets come from the "
            f"documented Ω·n̂ convention; a divergence is either an ordinate "
            f"re-classification (the > 0.0 tangential twin) or a changed "
            f"quadrature — never something to re-baseline away."
        ),
    )
    np.testing.assert_array_equal(
        np.sort(np.asarray(space.outflow_indices)), snapshot["outflow_indices"],
        err_msg=f"{case.case_id}: the live Γ₊ differs from the frozen one.",
    )
    return space


def _realize(case: BCEquivalenceCase, space) -> LinearOperator:
    """Build the case's operator through the registry's own composition."""
    realizer = SNBoundaryRealizer()
    return case.compose(lambda law: realizer.realize(law, space))


def _reduction_order_rtol(snapshot: np.lib.npyio.NpzFile, extra: int = 0) -> float:
    r"""The reduction-order bound for a :math:`|\Gamma_+|`-term positive sum.

    ``(|Γ₊| + extra) · ε``. Recursive summation of :math:`n` terms carries a
    relative error :math:`\le (n-1)\,u\,\kappa` with :math:`u = \varepsilon/2`;
    every summand here is positive (the probe is ``U(0,2)``, the quadrature
    weights and :math:`|\Omega\cdot\hat n|` are positive on :math:`\Gamma_+`),
    so :math:`\kappa = 1` and the bound is clean. ``|Γ₊| · ε`` is that bound
    with roughly a factor two of room, which absorbs the trailing division;
    ``extra`` buys one rounding per additional operation downstream.

    This is DERIVED, not a measurement rounded up — the point of replacing the
    ``nulp=4`` / ``nulp=64`` constants this file used to carry, which encoded no
    claim at all.
    """
    return float(int(snapshot["outflow_indices"].size) + extra) * _EPS


# ─────────────────────────────────────────────────────────────────────
# Layer-0 — the inversion itself, enforced rather than documented
# ─────────────────────────────────────────────────────────────────────


def _imported_names(path: Path) -> set[str]:
    """Every module and symbol name an ``import`` statement in ``path`` binds."""
    names: set[str] = set()
    for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
        if isinstance(node, ast.Import):
            names.update(a.name for a in node.names)
        elif isinstance(node, ast.ImportFrom):
            module = node.module or ""
            names.add(module)
            names.update(f"{module}.{a.name}" for a in node.names)
    return names


class TestTheInversionIsStructural:
    r"""The generator computes a reference; it does not record production.

    This is the whole design decision of the 2026-08-01 re-anchoring, and a
    docstring cannot hold it: re-adding one import would silently turn the
    artefacts back into recordings, and every row in this file would keep
    passing while claiming something it no longer does. So it is asserted.

    Both directions are load-bearing. If the generator imported the realizer,
    the artefact would be a recording. If this harness imported the
    generator's ``_*_image`` expressions, the frozen file would stop being the
    barrier that keeps the reference and production from drifting *together* —
    the gate would recompute both sides from one source and go green on a
    shared mistake.
    """

    _GENERATOR = Path(__file__).parent / "_generate_bc_equivalence_snapshots.py"

    def test_the_generator_imports_no_realization_code(self) -> None:
        offenders = sorted(
            name for name in _imported_names(self._GENERATOR)
            if "realizer" in name
            or "method_space" in name
            or name.startswith("orpheus.sn")
            or name.endswith("SNBoundaryRealizer")
        )
        assert not offenders, (
            f"the reference generator imports {offenders} — it is a RECORDER "
            f"again, and its artefacts certify only that production still does "
            f"what it did. That coupling is what broke this harness at B3.2 "
            f"and B3.4a, and what the re-anchoring removed."
        )

    def test_this_harness_imports_only_the_registry(self) -> None:
        generator_module = "tests.geometry._generate_bc_equivalence_snapshots"
        pulled = {
            name.rsplit(".", 1)[-1]
            for name in _imported_names(Path(__file__))
            if name.startswith(generator_module + ".")
        }
        allowed = {"BCEquivalenceCase", "CASES", "case_by_id"}
        assert pulled <= allowed, (
            f"this harness pulls {sorted(pulled - allowed)} from the "
            f"generator. Only the case REGISTRY may cross: importing an "
            f"expression helper would make the reference computable here, and "
            f"the frozen artefact — the one thing standing between a wrong "
            f"expression and a green gate — would stop mattering."
        )

    def test_the_registry_covers_every_committed_artefact(self) -> None:
        """No orphaned ``.npz``, no case without one.

        An artefact whose case was deleted is dead weight that looks like
        coverage; a case whose artefact is missing reds elsewhere, and this
        row names the asymmetry directly. (``albedo_05_lebedev17`` was retired
        at B3.4b together with its baseline, which is what this pins.)
        """
        on_disk = {
            p.stem.removeprefix("bc_equivalence_")
            for p in _SNAPSHOT_DIR.glob("bc_equivalence_*.npz")
        }
        registered = {c.case_id for c in CASES}
        assert on_disk == registered, (
            f"orphaned artefacts {sorted(on_disk - registered)}; "
            f"cases without an artefact {sorted(registered - on_disk)}."
        )


# ─────────────────────────────────────────────────────────────────────
# Layer-1 — the schema and the classification, every case
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("case", CASES, ids=[c.case_id for c in CASES])
class TestEveryArtefactIsSelfDescribing:
    """The frozen metadata is complete, and it agrees with the live geometry.

    Uniform across all seven cases (the pre-migration schema carried the index
    set for the vacuum case only, so six artefacts could not say which
    half-traces they were about). These rows are cheap and they run first: a
    value comparison against the wrong half-trace is worse than no comparison.
    """

    def test_the_schema_is_complete(self, case: BCEquivalenceCase) -> None:
        snapshot = _load_snapshot(case.case_id)
        required = {
            "psi_out", "psi_in", "face", "inflow_indices", "outflow_indices",
            "case_id", "description",
        }
        missing = required - set(snapshot.files)
        assert not missing, (
            f"{case.case_id}: artefact is missing {sorted(missing)}. "
            f"{_GENERATOR_HINT}"
        )
        assert str(snapshot["case_id"]) == case.case_id

    def test_the_narrowed_shapes_are_the_two_half_traces(
        self, case: BCEquivalenceCase,
    ) -> None:
        r"""``psi_out`` lives on :math:`\Gamma_+`, ``psi_in`` on :math:`\Gamma_-`.

        The B3.2 contract, read off the artefact. It cannot stand alone —
        :math:`|\Gamma_+| = |\Gamma_-|` on every quadrature in the tree, so a
        law that stayed an endomorphism would satisfy it (``vv-principles``
        Mode 12) — which is why every case below pairs it with a value.

        The trailing ``(5, 3)`` deliberately RESTATES the generator's
        ``_PSI_TRAILING_SHAPE`` rather than importing it: a silent change to
        the probe shape should red here and force a look, and the allowed-
        import gate above keeps the registry the only thing that crosses. It
        is not square, so a group ↔ spatial axis transpose cannot hide behind
        it.
        """
        snapshot = _load_snapshot(case.case_id)
        n_in = int(snapshot["inflow_indices"].size)
        n_out = int(snapshot["outflow_indices"].size)
        assert snapshot["psi_out"].shape == (n_out, 5, 3)
        assert snapshot["psi_in"].shape == (n_in, 5, 3)

    def test_the_frozen_half_traces_match_the_live_trace_space(
        self, case: BCEquivalenceCase,
    ) -> None:
        """The classification gate; see :func:`_verified_space`."""
        _verified_space(case, _load_snapshot(case.case_id))

    def test_the_probe_activates_the_law(self, case: BCEquivalenceCase) -> None:
        """NON-VACUITY of the input, for every case including vacuum.

        Vacuum's reference is zero, so its value row cannot assert a non-zero
        image; what makes that row a statement is that the operator was handed
        something non-zero and returned nothing. Asserting it here keeps the
        property in one place for all seven cases, and keeps it from being
        confused with "the answer is non-zero".
        """
        snapshot = _load_snapshot(case.case_id)
        assert np.count_nonzero(snapshot["psi_out"]), (
            f"{case.case_id}: the probe is all zeros, so every comparison "
            f"below is 0 == 0."
        )


# ─────────────────────────────────────────────────────────────────────
# Case 1 — Vacuum × Lebedev(17) at xmin
# ─────────────────────────────────────────────────────────────────────


class TestVacuumLebedev17Snapshot:
    r"""Vacuum: :math:`R = q = 0`, so :math:`\gamma_-\psi = 0`.

    Since B3.2 the realized law is typed :math:`\Gamma_+ \to \Gamma_-`, so
    vacuum is the honest ZERO MAP and there is no "outflow pass-through" half
    left to check — those rows are outside the operator's domain. (Pre-B3.2
    this case carried a zero-mask check AND a pass-through check, because the
    §16A.5 realization was a full-face projector that preserved the outflow
    rows.)

    ``assert_array_equal``: the zero map performs no arithmetic.
    """

    case_id = "vacuum_lebedev17"

    def test_realizes_the_zero_map_onto_gamma_minus(self) -> None:
        case = case_by_id(self.case_id)
        snapshot = _load_snapshot(self.case_id)
        space = _verified_space(case, snapshot)
        actual = _realize(case, space).apply(snapshot["psi_out"])
        assert actual.shape == snapshot["psi_in"].shape, (
            f"vacuum emitted {actual.shape}; the narrowed codomain is Γ₋, "
            f"which has {snapshot['psi_in'].shape[0]} rows. A zero map that "
            f"echoes its input's shape is the endomorphism this narrowing "
            f"removed — and it looks correct here only because |Γ₊| == |Γ₋|."
        )
        np.testing.assert_array_equal(
            actual, snapshot["psi_in"],
            err_msg="the narrowed vacuum law is the ZERO map Γ₊ → Γ₋",
        )


# ─────────────────────────────────────────────────────────────────────
# Case 2 — Specular(x, α=1) × Lebedev(17) at xmax
# ─────────────────────────────────────────────────────────────────────


class TestSpecularXLebedev17Snapshot:
    r"""Specular at :math:`\alpha = 1`: the mirror gather.

    Reference: :math:`\psi^-(\Omega) = \psi^+(\Omega')` with
    :math:`\Omega' = \Omega - 2(\Omega\cdot\hat n)\hat n`, the pairing derived
    from the isometry and matched onto the ordinate table by exact float
    equality — never read from production's pairing derivation
    (``quadrature.ordinate_permutation``; until §7d.3, the
    ``reflection_index`` table).

    ``assert_array_equal``: a gather introduces no re-association, so the
    predicted drift is EXACTLY zero and any tolerance would admit the failure
    mode the row exists to catch. MEASURED bit-identical, and additionally
    bit-identical to the retired pre-B3.2 recording restricted to
    :math:`\Gamma_-`.
    """

    case_id = "specular_x_lebedev17"

    def test_matches_the_frozen_mirror_gather(self) -> None:
        case = case_by_id(self.case_id)
        snapshot = _load_snapshot(self.case_id)
        space = _verified_space(case, snapshot)
        actual = _realize(case, space).apply(snapshot["psi_out"])
        np.testing.assert_array_equal(actual, snapshot["psi_in"])
        assert np.count_nonzero(snapshot["psi_in"]), "vacuous comparison"

    def test_albedo_with_specular_closure_folds_alpha_onto_the_same_image(
        self,
    ) -> None:
        r"""**B3.4b** — the ``≡`` theorem, pinned against the MATHEMATICS.

        ``AlbedoBoundary(α, SpecularReturn("x")) ≡ ReflectiveBoundary("x", α)``
        as realized operators, so the albedo route's image is the α = 1 mirror
        gather scaled by α.

        Before the re-anchoring this row's argument was provenance: the
        reference was an artefact frozen pre-B3.2 by a *different* law. The
        migration replaced that with something stronger — the reference is now
        the mirror isometry itself, so this checks the albedo route against the
        equation rather than against a sibling implementation's recorded
        output. (:class:`tests.geometry.test_reemission_closure.TestEquivalenceTheorems`
        checks the two routes against EACH OTHER and is structurally blind to a
        bug in the body they share; this row and the anchors next door are what
        catch that.)

        ``lebedev(17)`` is one of the two quadratures where the specular
        pairing and the retired array-position pairing genuinely DISAGREE, so
        the row discriminates.
        """
        case = case_by_id(self.case_id)
        snapshot = _load_snapshot(self.case_id)
        space = _verified_space(case, snapshot)
        op = SNBoundaryRealizer().realize(
            AlbedoBoundary(0.5, SpecularReturn(axis="x")), space,
        )
        actual = op.apply(snapshot["psi_out"])
        expected = 0.5 * snapshot["psi_in"]
        assert np.count_nonzero(expected), "vacuous comparison"
        np.testing.assert_array_equal(actual, expected)


# ─────────────────────────────────────────────────────────────────────
# Case 3 — Specular(y, α=0.7) × LevelSymmetric(6) at ymax
# ─────────────────────────────────────────────────────────────────────


class TestSpecularYPartial07LS6Snapshot:
    r"""Specular about **y**, at :math:`\alpha = 0.7`.

    The α-fold row, and the only case on a non-``x`` axis and a ``y`` face — so
    an ``axis`` or ``outward_sign`` hard-coded to the x-convention reds here.

    ``assert_array_equal``: scale-then-gather is one multiplication of the same
    two floats followed by a permutation; reduction depth 0.
    """

    case_id = "specular_y_partial_07_LS6"

    def test_matches_the_frozen_scaled_mirror_gather(self) -> None:
        case = case_by_id(self.case_id)
        snapshot = _load_snapshot(self.case_id)
        space = _verified_space(case, snapshot)
        actual = _realize(case, space).apply(snapshot["psi_out"])
        np.testing.assert_array_equal(actual, snapshot["psi_in"])
        assert np.count_nonzero(snapshot["psi_in"]), "vacuous comparison"


# ─────────────────────────────────────────────────────────────────────
# Case 4 — White(xmax, α=1) × LevelSymmetric(4)
# ─────────────────────────────────────────────────────────────────────


class TestWhiteXmaxLS4Snapshot:
    r"""White (Lambertian) at :math:`\alpha = 1`.

    Reference: isotropic re-emission carrying the incident partial current,
    :math:`J^- = J^+` in the boundary measure
    :math:`w_n|\Omega_n\cdot\hat n|`. The generator derives it from that
    balance; production spells the same number as *the cosine-weighted average
    of the outflow* (it normalises over :math:`\Gamma_+` where the balance
    normalises over :math:`\Gamma_-`). The two agree because the mirror is an
    exact measure-preserving bijection between the half-traces, which the
    generator asserts rather than assumes.

    ``rtol = |Γ₊| · ε`` — see :func:`_reduction_order_rtol`. MEASURED max
    relative error 1.9e-16 against a bound of 2.7e-15.
    """

    case_id = "white_xmax_LS4"

    def test_matches_the_frozen_lambertian(self) -> None:
        case = case_by_id(self.case_id)
        snapshot = _load_snapshot(self.case_id)
        space = _verified_space(case, snapshot)
        actual = _realize(case, space).apply(snapshot["psi_out"])
        np.testing.assert_allclose(
            actual, snapshot["psi_in"],
            rtol=_reduction_order_rtol(snapshot), atol=0.0,
            err_msg=(
                "the realized Lambertian disagrees with the frozen "
                "partial-current-balance reference beyond the |Γ₊|·ε "
                "reduction-order bound"
            ),
        )
        assert np.count_nonzero(snapshot["psi_in"]), "vacuous comparison"

    def test_the_frozen_image_is_isotropic_over_gamma_minus(self) -> None:
        """The defining Lambertian law, as an intrinsic property of the artefact.

        A re-emission that is not constant over the incoming hemisphere is not
        Lambertian however well it matches an average. Asserted on the frozen
        reference because that is what the value row above compares against: if
        the artefact ever stopped being isotropic, the row would still pass
        while pinning something that is not a white boundary.
        """
        snapshot = _load_snapshot(self.case_id)
        psi_in = snapshot["psi_in"]
        np.testing.assert_array_equal(
            psi_in, np.broadcast_to(psi_in[0][None, ...], psi_in.shape),
        )


# ─────────────────────────────────────────────────────────────────────
# Case 5 — White(xmin, α=0.3) × GaussLegendre1D(8)
# ─────────────────────────────────────────────────────────────────────


class TestWhiteXminPartial03GLSnapshot:
    r"""White at :math:`\alpha = 0.3` on a 1-D Gauss-Legendre quadrature.

    Same law as Case 4 with the amplitude folded in, on a quadrature whose
    :math:`\sum w` is 2 rather than :math:`4\pi` and on the ``xmin`` face
    (outward normal :math:`-\hat x`). It is the canary against a hard-coded
    normalisation constant — ``numerical-bug-signatures`` Signature 4, the
    ERR-004 family — which a spherical-only case set cannot see.

    ``rtol = |Γ₊| · ε`` with :math:`|\Gamma_+| = 4`, i.e. 8.9e-16. MEASURED
    1.2e-16.
    """

    case_id = "white_xmin_partial_03_GL"

    def test_matches_the_frozen_scaled_lambertian(self) -> None:
        case = case_by_id(self.case_id)
        snapshot = _load_snapshot(self.case_id)
        space = _verified_space(case, snapshot)
        actual = _realize(case, space).apply(snapshot["psi_out"])
        np.testing.assert_allclose(
            actual, snapshot["psi_in"],
            rtol=_reduction_order_rtol(snapshot), atol=0.0,
        )
        assert np.count_nonzero(snapshot["psi_in"]), "vacuous comparison"


# ─────────────────────────────────────────────────────────────────────
# Case 6 — Mixed (0.3 Specular + 0.7 White) × LevelSymmetric(4) at xmax
# ─────────────────────────────────────────────────────────────────────


class TestMixed30Spec70WhiteLS4Snapshot:
    r"""The rank-N Marshak mix, realized through the Wave-0 operator algebra.

    .. math::

        R = 0.3\,R_{\text{spec}} + 0.7\,R_{\text{diff}}

    Each leaf is realized independently against the SAME face-ful method space
    and the composition is an
    :class:`~orpheus.numerics.operator.OperatorSum` of
    :class:`~orpheus.numerics.operator.ScaledOperator` wrappers — stated once,
    in the case registry, so the generator and this harness cannot describe
    different mixtures.

    The reference is the pointwise convex combination of the two independent
    images. Before the re-anchoring the generator composed
    ``0.3 * spec_realized + 0.7 * white_realized`` out of REALIZED leaves,
    which made this the most self-referential row in the file: the snapshot and
    the test were the same expression.

    Both leaves are narrowed (B3.2 reflective, B3.4a white), which is what
    makes the sum well-typed. It must NOT be re-posed over albedo or periodic:
    those are still full-face, and because :math:`|\Gamma_+| = |\Gamma_-|` such
    a sum does not RAISE — it runs silently wrong (``vv-principles`` Mode 12),
    so the gate would be green and worthless.

    ``rtol = (|Γ₊| + 2) · ε``: the white term's bound, plus one rounding for
    its 0.7 scaling and one for the sum. No cancellation (both weights and both
    images are positive) so the condition number stays 1. MEASURED 3.2e-16
    against a bound of 3.1e-15.
    """

    case_id = "mixed_30spec_70white_LS4"

    def test_matches_the_frozen_convex_combination(self) -> None:
        case = case_by_id(self.case_id)
        snapshot = _load_snapshot(self.case_id)
        space = _verified_space(case, snapshot)
        actual = _realize(case, space).apply(snapshot["psi_out"])
        np.testing.assert_allclose(
            actual, snapshot["psi_in"],
            rtol=_reduction_order_rtol(snapshot, extra=2), atol=0.0,
        )
        assert np.count_nonzero(snapshot["psi_in"]), "vacuous comparison"

    def test_the_mix_is_neither_leaf(self) -> None:
        r"""ACTIVATION: the frozen image is not one of its own summands.

        A mixing weight dropped to ``(1.0, 0.0)`` or ``(0.0, 1.0)`` would make
        the row above a restatement of Case 2 or Case 4 — green, and blind to
        the composition it exists to test. The specular leaf varies over
        :math:`\Gamma_-` and the Lambertian one does not, so "the image is
        neither constant nor a pure gather" separates the mix from both.
        """
        psi_in = _load_snapshot(self.case_id)["psi_in"]
        assert not np.allclose(psi_in, psi_in[0][None, ...]), (
            "the mixed image is CONSTANT over Γ₋, i.e. purely Lambertian — "
            "the specular weight has been lost."
        )
        spread = psi_in.max(axis=0) - psi_in.min(axis=0)
        assert np.all(spread > 0.0), (
            "the mixed image varies over Γ₋ in some entries but not others"
        )


# ─────────────────────────────────────────────────────────────────────
# Case 7 — Periodic × Lebedev(17), xmin fed by xmax
# ─────────────────────────────────────────────────────────────────────


# This case carried ``xfail(strict=True)`` from the re-anchoring until
# **B3.4c** landed. The documented defect: periodic is a TWO-FACE coupling and
# production realized it as a per-face angular identity, so the composition
# ``ι₋ ∘ law ∘ γ₊`` handed it this face's own Γ₊ = {μ_x < 0} and it fed
# outgoing-left flux back in as incoming-right flux — MEASURED 98 % relative
# error, an O(1) wrong answer rather than a tolerance question.
#
# B3.4c built the partner-face channel and the marker is DELETED. Note what
# replaced it, because a mechanical "swap the `apply` line to the partner
# probe" would have produced a VERBATIM DUPLICATE of the live body-pin below:
# B3.4c did not change the operator's body, it changed which half-trace the
# composition supplies. So the re-posed gate asks the LAW which face to read
# (``geometry_map.domain_face``) and feeds the probe that answer names. A
# regression of ``domain_face`` back to the installation face now selects the
# wrong probe and reds it, which is precisely the defect the xfail described.


class TestPeriodicLebedev17Snapshot:
    r"""Periodic: the partner face's outflow IS this face's inflow.

    .. math::

        \psi^-\!\left(\mathbf{r}_f, \Omega\right)
        \;=\; \psi^+\!\left(T\mathbf{r}_f, \Omega\right)

    where :math:`T` is the translation that identifies the two faces. The
    artefact therefore carries TWO probes — this face's :math:`\Gamma_+` and
    the partner's — drawn from independent seeds. That is not decoration: with
    one shared draw, a per-face endomorphism would look correct on the rows
    that coincide, and the whole content of a periodic law is that it reads the
    *other* face.

    Before the re-anchoring this case pinned ``PeriodicWrapOperator``'s
    identity body against a recording of itself — the most vacuous row in the
    file, and green throughout. Re-anchoring it to the partner-face identity
    converts an artefact that pinned a placeholder into a gate that pins the
    requirement.
    """

    case_id = "periodic_lebedev17"

    def test_the_two_faces_are_a_periodic_pair(self) -> None:
        r"""LIVE, and it is the premise everything else here rests on.

        :math:`\Gamma_-(\text{xmin})` and :math:`\Gamma_+(\text{xmax})` must be
        the SAME ordinate set, because the translation acts on space and not on
        direction. If they ever differ, the two faces are not a periodic pair
        and the xfail below would be red for the wrong reason.
        """
        case = case_by_id(self.case_id)
        snapshot = _load_snapshot(self.case_id)
        _verified_space(case, snapshot)
        np.testing.assert_array_equal(
            snapshot["inflow_indices"], snapshot["partner_outflow_indices"],
            err_msg=(
                "Γ₋(xmin) != Γ₊(xmax): a periodic pair must agree on the "
                "ordinate set."
            ),
        )
        assert str(snapshot["partner_face"]) == "xmax"

    def test_the_two_probes_discriminate(self) -> None:
        """LIVE. The fixture ACTIVATES the defect, asserted rather than assumed.

        The xfail below is a claim that production reads the wrong face. That
        claim is only observable if the two faces carry different data on the
        rows in question — a shared draw would null it exactly (the
        config-blindness discipline: the convenient fixture nulls the term).
        """
        snapshot = _load_snapshot(self.case_id)
        assert not np.array_equal(
            snapshot["psi_out"], snapshot["psi_out_partner"],
        ), (
            "the two faces carry identical data, so 'production reads the "
            "wrong face' is unobservable and the xfail below is vacuous."
        )

    def test_the_identity_body_is_correct_on_the_partner_half_trace(self) -> None:
        """LIVE, and it is the flip-proof for the xfail below.

        The realized operator is an angular identity. Fed the PARTNER face's
        :math:`\\Gamma_+` — the law's true domain — it already reproduces the
        reference exactly. So the defect is entirely in WHICH half-trace the
        composition supplies, which is what B3.4c changes; the body needs no
        repair.

        It also pins the body: an operator that scaled, permuted or averaged
        would fail here, so this is not merely ``x == x``.
        """
        case = case_by_id(self.case_id)
        snapshot = _load_snapshot(self.case_id)
        space = _verified_space(case, snapshot)
        actual = _realize(case, space).apply(snapshot["psi_out_partner"])
        np.testing.assert_array_equal(actual, snapshot["psi_in"])

    def test_delivers_the_partner_faces_outflow(self) -> None:
        r"""The REQUIREMENT, live since **B3.4c**:
        :math:`\gamma_-\psi|_{\rm xmin} = \gamma_+\psi|_{\rm xmax}`.

        The probe is chosen by the LAW, not by this test: it asks
        ``geometry_map.domain_face`` which face's :math:`\Gamma_+` the law
        consumes and feeds the half-trace that answer names. That is the whole
        thing B3.4c built, and asking it here is what keeps this from being a
        second copy of the body-pin above — that one hard-codes the partner
        probe to pin the *body*, this one lets the law route and so pins the
        *routing*. If ``domain_face`` ever regresses to the installation face,
        this selects ``psi_out`` and reds at 98 % relative error, which is the
        exact defect its retired ``xfail(strict=True)`` documented.
        """
        case = case_by_id(self.case_id)
        snapshot = _load_snapshot(self.case_id)
        space = face_method_space(
            case.build_quadrature(), face=case.face, faces=case.faces,
        )
        # The law comes from the REGISTRY's own composition, captured on the
        # way through, never re-spelled here: a second spelling could name a
        # different axis than the case realizes and the two would agree by
        # coincidence (the Pattern-2 twin the case dataclass exists to kill).
        composed_laws: list = []

        def _capture(law):
            composed_laws.append(law)
            return SNBoundaryRealizer().realize(law, space)

        operator = case.compose(_capture)
        (law,) = composed_laws
        probes = {
            case.face: snapshot["psi_out"],
            str(snapshot["partner_face"]): snapshot["psi_out_partner"],
        }
        domain_face = law.geometry_map.domain_face(case.face)
        assert domain_face in probes, (
            f"the law names {domain_face!r} as its domain face, which is "
            f"neither this face nor its recorded partner — the artefact "
            f"carries no probe for it."
        )
        actual = operator.apply(probes[domain_face])
        np.testing.assert_array_equal(actual, snapshot["psi_in"])
