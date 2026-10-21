r"""S4 — the Green owns its realisation (numerics rows of the S4 gate spec).

Spec: ``scratch/numerics_untangling/packages/S4_gate_spec.md`` (rows G4–G6;
the SN route rows G1–G3 are ``tests/gates/sn/solve/test_strategy_mints_green_sn.py``).

* G4 — ruling 67, one seed law for a scaled inverse :math:`(\alpha A)^{-1}`:
  the forward face, the transpose face and the driver slot each hand the word
  :math:`\alpha x_0`, keyed by ONE capability predicate
  (:func:`~orpheus.numerics.iteration.accepts_seed`).
* G5 — the promise face over a Krylov realisation (ruling 68).
* G6 — ruling 66: the Green's adjoint stays the base ``♯ ∘ dual ∘ ♭``, which is
  what keeps the honesty registry's ``_sn_green`` an independent witness.
"""
from __future__ import annotations

import ast
import inspect
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from orpheus.numerics import iteration as iteration_mod
from orpheus.numerics.green_operator import GreenOperator
from orpheus.numerics.iteration import KrylovAcceleration, SourceIteration, seeded_inverse
from orpheus.numerics.operator import IdentityOperator, InverseWrapMixin, LinearOperator
from orpheus.numerics.elementary_operators import DiagonalOperator

pytestmark = [pytest.mark.foundation]

_ORPHEUS = Path(iteration_mod.__file__).resolve().parents[1]
_GREEN_FILE = (_ORPHEUS / "numerics" / "green_operator.py").resolve()


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


from tests.gates.numerics.test_linear_combination import Mat, V  # noqa: E402  (the S2 doubles)


class _RecordingInverse(InverseWrapMixin, LinearOperator):
    """A wrap-delegate inverse that RECORDS every seed it is handed (forward
    and transpose). Exact, so a seed never moves its value: the record is the
    only observable of the seed path ([M] the 26-configuration baseline is
    value-blind to the driver-slot seed, spec finding F3)."""

    def __init__(self, fwd) -> None:
        super().__init__(fwd)
        self.seeds: list[Any] = []
        self.transpose_seeds: list[Any] = []

    def apply(self, b, /, *, initial_guess=None):
        self.seeds.append(None if initial_guess is None else np.array(initial_guess))
        return np.linalg.solve(self.inner.M, b)

    @property
    def is_adjointable(self) -> bool:
        return True

    def apply_transpose(self, b, /, *, initial_guess=None):
        self.transpose_seeds.append(None if initial_guess is None else np.array(initial_guess))
        return np.linalg.solve(self.inner.M.T, b)


class _WrapForward(Mat):
    """An invertible bound diagonal whose ``inverse()`` is the recorder."""

    def __init__(self) -> None:
        super().__init__(np.diag([2.0, 3.0, 4.0, 5.0]), V, V, invertible=True)
        self.rec = _RecordingInverse(self)

    def inverse(self):
        return self.rec


_X = np.array([1.0, -2.0, 0.5, 4.0])
_X0 = np.array([3.0, 1.0, -1.0, 2.0])
_ALPHA = 2.0


# ═══════════════════════════════════════════════════════════════════════
# G4 — ruling 67: one seed law for a scaled inverse (numerics)
# ═══════════════════════════════════════════════════════════════════════


def test_g4a_forward_face_seeds_the_word_with_alpha_x0():
    """PRE-CARVE-RUNNABLE. ``(αA)⁻¹.apply(b, initial_guess=x₀)`` hands A⁻¹ the
    seed α·x₀ (the seed lives in the output space of (αA)⁻¹ = α⁻¹A⁻¹).
    Re-states c30 on the recording double so G4a/G4b/G4c read one fixture.
    First red: ``seed / c`` → ``seed`` in ``OperatorSum._apply_scaled``."""
    fwd = _WrapForward()
    inv = fwd.rec
    (_ALPHA * fwd).inverse().apply(_X, initial_guess=_X0)
    _require(len(inv.seeds) == 1 and inv.seeds[0] is not None,
             "the forward face dropped the seed of a scaled inverse")
    np.testing.assert_array_equal(inv.seeds[0], _ALPHA * _X0)


def test_g4b_transpose_face_seeds_the_word_with_alpha_x0():
    """PRE-CARVE-RUNNABLE, NEW (0 tests today, census §3). The transpose face of
    the one-term combination threads α·x₀ into the word's ``apply_transpose``.
    First red: ``seed / c`` → ``seed`` at the one-term ``apply_transpose``
    (operator.py, the qa-F3 fix), or the seed dropped there."""
    fwd = _WrapForward()
    inv = fwd.rec
    (_ALPHA * fwd).inverse().apply_transpose(_X, initial_guess=_X0)
    _require(len(inv.transpose_seeds) == 1 and inv.transpose_seeds[0] is not None,
             "the transpose face dropped the seed of a scaled inverse")
    np.testing.assert_array_equal(inv.transpose_seeds[0], _ALPHA * _X0)


def test_g4c_driver_slot_threads_the_seed_of_a_scaled_inverse():
    """RED ON ae6bdbeb (its first red: ``_SeededExactApply`` drops the seed, the
    recorder sees None). The driver slot: ``seeded_inverse(αA)`` placed in a
    SourceIteration hands A⁻¹ the iterate's seed scaled by α, the same law as
    the direct face (ruling 67). One step from x₀, so the seed is x₀."""
    fwd = _WrapForward()
    inv = fwd.rec
    slot = seeded_inverse(_ALPHA * fwd)
    zero_gain = 0.0 * fwd
    SourceIteration(slot, zero_gain, max_iter=1, tol=0.0).solve(_X, initial_guess=_X0)
    _require(inv.seeds and inv.seeds[0] is not None,
             "the driver slot dropped the seed of (αA)⁻¹ (ruling 67: one seed law)")
    np.testing.assert_array_equal(inv.seeds[0], _ALPHA * _X0)


def test_g4d_a_word_that_refuses_a_seed_drops_it_in_the_driver_slot():
    """PRE-CARVE-RUNNABLE. (2I)⁻¹ in the driver slot must not raise and returns
    the exact value (the seed is dropped, Q1); (2D)⁻¹ is a one-term combination
    over a wrap member and is threaded (its pair row is G4c). Control for the
    capability predicate: an implementation that passes every one-term
    combination as-is reds HERE with ``TypeError`` from IdentityOperator."""
    slot = seeded_inverse(_ALPHA * IdentityOperator())
    np.testing.assert_array_equal(slot.apply(_X, initial_guess=_X0), _X / _ALPHA)


def test_g4e_a_combination_over_a_seedless_word_refuses_the_seed_directly():
    """PRE-CARVE-RUNNABLE, CHARACTERISATION (Q1: "as is"; the orchestrator may
    change it in the carve and will say so). ``(2I)⁻¹.apply(b, initial_guess=x)``
    raises TypeError today: the one-term combination forwards the keyword to a
    word whose apply has no seed slot."""
    with pytest.raises(TypeError, match="initial_guess"):
        (_ALPHA * IdentityOperator()).inverse().apply(_X, initial_guess=_X0)


def test_g4f_the_seed_capability_is_one_named_predicate():
    """POST-CARVE. Q1: the key is ONE named predicate read by ``seeded_inverse``,
    never an inline isinstance chain. AST of ``seeded_inverse``: no
    ``isinstance`` call in its body; it calls ``accepts_seed`` (name a
    hypothesis); and the predicate itself answers the four members
    (wrap → True; (2D)⁻¹ → True; I → False; (2I)⁻¹ → False)."""
    tree = ast.parse(inspect.getsource(iteration_mod.seeded_inverse))
    calls = [n.func for n in ast.walk(tree) if isinstance(n, ast.Call)]
    names = {getattr(f, "id", getattr(f, "attr", None)) for f in calls}
    _require("isinstance" not in names, "seeded_inverse branches by isinstance inline (Q1)")
    _require("accepts_seed" in names, f"seeded_inverse does not read the named predicate; calls {names}")
    pred = iteration_mod.accepts_seed
    d = DiagonalOperator(np.array([1.0, 2.0, 4.0, 8.0]))
    _require(pred(d.inverse()) is True, "a wrap member must accept a seed")
    _require(pred((_ALPHA * d).inverse()) is True, "(2D)⁻¹ must accept a seed")
    _require(pred(IdentityOperator().inverse()) is False, "I⁻¹ has no seed slot")
    _require(pred((_ALPHA * IdentityOperator()).inverse()) is False, "(2I)⁻¹'s word has no seed slot")


# ═══════════════════════════════════════════════════════════════════════
# G5 — the promise face over a Krylov realisation (numerics, NEW)
# ═══════════════════════════════════════════════════════════════════════

_RNG = np.random.default_rng(2026)
_n = 12
_Mmat = np.diag(np.linspace(2.0, 5.0, _n))
_Nmat = 0.45 * np.abs(_RNG.standard_normal((_n, _n))) / _n  # ρ(M⁻¹N) < 1
_Q = _RNG.standard_normal(_n)


def _split_ops():
    from tests.gates.numerics.test_iteration import MatrixOperator

    return MatrixOperator(_Mmat, can_solve=True), MatrixOperator(_Nmat)


def _green(realisation_kind: str, *, tol: float = 1e-10, max_iter: int = 200):
    from orpheus.numerics.green_operator import Krylov

    M, N = _split_ops()
    if realisation_kind == "krylov":
        real = Krylov(M.inverse(), restart=_n)
    else:  # the default realisation: Richardson on M's own seeded inverse
        return GreenOperator(M, N, tol=tol, max_iter=max_iter)
    return GreenOperator(M, N, realisation=real, tol=tol, max_iter=max_iter)


def test_g5a_a_krylov_realised_green_keeps_the_promise():
    """POST-CARVE, REFERENCE. ``G.apply(q)`` over a Krylov realisation returns
    ψ with TRUE relative residual ‖(M−N)ψ − q‖/‖q‖ < tol, and ψ agrees with
    ``np.linalg.solve(M−N, q)`` within κ·tol (κ the dense 2-norm condition
    number, measured at landing). First red: the promise check skipped for
    the Krylov arm (apply returns the driver's raw iterate)."""
    g = _green("krylov")
    psi = g.apply(_Q)
    A = _Mmat - _Nmat
    res = np.linalg.norm(A @ psi - _Q) / np.linalg.norm(_Q)
    _require(res < g.tol, f"Krylov-realised Green broke its promise: residual {res:.2e}")
    ref = np.linalg.solve(A, _Q)
    kappa = np.linalg.cond(A)
    err = np.linalg.norm(psi - ref) / np.linalg.norm(ref)
    _require(err <= 10.0 * kappa * g.tol, f"Krylov Green value off by {err:.2e} (κ={kappa:.1f})")


def test_g5b_a_krylov_realised_green_refuses_an_unmet_promise():
    """POST-CARVE. ``tol=0`` is unsatisfiable: ``apply`` raises
    ``ConvergenceFailure`` naming the TRUE relative residual (fragment pinned
    by test_green_operator.py:271, B.11). First red: the Krylov arm returns
    its best effort silently."""
    from orpheus.numerics.green_operator import ConvergenceFailure

    with pytest.raises(ConvergenceFailure, match="TRUE relative residual"):
        _green("krylov", tol=0.0, max_iter=20).apply(_Q)


def test_g5c_the_promise_is_written_over_the_relaxation():
    """POST-CARVE, ROUTE. ``apply`` drives the realisation through ``self.driver``
    (ruling 68: the relaxation is constructed in one place): a spy on
    ``GreenOperator.driver`` counts ≥ 1 call per ``apply`` for BOTH
    realisations, and the drivers built are KrylovAcceleration / SourceIteration
    respectively. First red: apply constructing its driver directly."""
    calls: list[str] = []
    real = GreenOperator.driver

    def spy(self, **kw):
        d = real(self, **kw)
        calls.append(type(d).__name__)
        return d

    mp = pytest.MonkeyPatch()
    mp.setattr(GreenOperator, "driver", spy)
    try:
        _green("krylov").apply(_Q)
        _green("richardson").apply(_Q)
    finally:
        mp.undo()
    _require("KrylovAcceleration" in calls and "SourceIteration" in calls,
             f"apply did not build its drivers through G.driver: {calls}")


def test_g5d_two_realisations_one_green():
    """POST-CARVE, THEOREM (the module docstring's "a Richardson-realized Green
    and a Krylov-realized Green are the SAME Green operator"). The two
    answers agree within 2·κ·tol. First red: a realisation iterating a
    different operator (e.g. the Krylov matvec reading M only)."""
    a, b = _green("krylov").apply(_Q), _green("richardson").apply(_Q)
    kappa = np.linalg.cond(_Mmat - _Nmat)
    _require(np.linalg.norm(a - b) / np.linalg.norm(b) <= 2 * 10.0 * kappa * 1e-10,
             "the two realisations of one Green disagree beyond the promise")


# ═══════════════════════════════════════════════════════════════════════
# G6 — ruling 66: the Green's adjoint stays the base ♯∘dual∘♭
# ═══════════════════════════════════════════════════════════════════════


def test_g6a_green_does_not_override_its_adjoint():
    """PRE- AND POST-CARVE. AST of green_operator.py: the GreenOperator class
    defines no ``H``/``adjoint`` (it may define ``dual``/``apply_transpose``).
    First red: re-introducing ``def adjoint(self): return GreenOperator(M.H, N.H)``,
    the override removed on 2026-10-18 that made the honesty witness a
    tautology (``_sn_green``, test_adjointability_is_honest.py)."""
    tree = ast.parse(_GREEN_FILE.read_text())
    cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == "GreenOperator")
    defined = {n.name for n in cls.body if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))}
    _require("dual" in defined, "control: the filter must find GreenOperator.dual")
    _require(not ({"H", "adjoint"} & defined), f"GreenOperator overrides its adjoint: {defined}")


def test_g6b_no_production_green_is_built_on_daggered_parts_outside_the_strategy():
    """POST-CARVE. AST over orpheus/: every ``GreenOperator(`` call whose first
    argument is an ``.H`` attribute lives in ``orpheus/sn/solver.py`` (the
    adjoint entries, design §2: the Green of the ADJOINT SPLIT minted by the
    Strategy) — none in numerics/, where it would be a second route to G.H.
    Population printed; control: the dual's own ``GreenOperator(`` in
    green_operator.py is found by the filter (its args are ``.dual()`` calls)."""
    hits, all_calls = [], []
    for f in sorted(_ORPHEUS.rglob("*.py")):
        for n in ast.walk(ast.parse(f.read_text())):
            if isinstance(n, ast.Call) and getattr(n.func, "id", getattr(n.func, "attr", None)) == "GreenOperator":
                all_calls.append(f)
                a0 = n.args[0] if n.args else None
                if isinstance(a0, ast.Attribute) and a0.attr == "H":
                    hits.append(f.relative_to(_ORPHEUS).as_posix())
    _require(_GREEN_FILE in [p.resolve() for p in all_calls], "control: dual's GreenOperator( not found")
    _require(all(h == "sn/solver.py" for h in hits), f"daggered-part Greens outside the Strategy: {hits}")




# ═══════════════════════════════════════════════════════════════════════
# G7 — KEigenvalue's inner is its Green's relaxation face
# ═══════════════════════════════════════════════════════════════════════


def test_g7_keigenvalue_builds_its_inner_through_the_greens_driver(monkeypatch):
    """``KEigenvalue(posing, green, max_inner=7, inner_tol=1e-9)`` builds its
    inner solve through ``green.driver`` with exactly that budget, tolerance and
    the ``max_inner`` budget name, and holds the Green it was handed (S4,
    ruling 68). Reds when the budget is not forwarded, or the inner is built
    by hand beside the Green."""
    from orpheus.numerics.iteration import KEigenvalue
    from orpheus.numerics.pencil import OperatorPencil
    from orpheus.numerics.posing import K_MAP, EigenPosing
    from tests.gates.numerics.test_iteration import MatrixOperator

    A = MatrixOperator(np.diag([3.0, 4.0]), can_solve=True)
    S = MatrixOperator(np.array([[0.5, 0.1], [0.2, 0.4]]))
    F = MatrixOperator(np.array([[1.0, 0.5], [0.3, 0.2]]))
    green = GreenOperator(A, S)
    seen: list[dict] = []
    real = GreenOperator.driver

    def spy(self, **kw):
        seen.append({"green": self, **kw})
        return real(self, **kw)

    monkeypatch.setattr(GreenOperator, "driver", spy)
    ke = KEigenvalue(EigenPosing(OperatorPencil(A - S, F), K_MAP), green, max_inner=7, inner_tol=1e-9)
    _require(ke.green is green, "KEigenvalue does not hold the Green it was handed")
    _require(len(seen) == 1, f"green.driver called {len(seen)}× at construction, expected 1")
    call = seen[0]
    _require(call["green"] is green and call["max_iter"] == 7 and call["tol"] == 1e-9
             and call["budget_name"] == "max_inner", f"the inner's budget was not forwarded: {call}")


# ═══════════════════════════════════════════════════════════════════════
# G9 — ruling 70: Richardson derives its preconditioner from M
# ═══════════════════════════════════════════════════════════════════════


def test_g9a_richardson_accepts_no_preconditioner():
    """The X1 mismatch ruling 70 closes cannot be spelled: ``Richardson`` takes
    no preconditioner, positionally or by keyword, so a source iteration over
    ANOTHER split's inverse is not a value (`[M]` the S4 elegance review: a
    free preconditioner of another split's diagonal gave a record "converged"
    15.5 % off). Reds on any re-introduced preconditioner parameter."""
    from orpheus.numerics.green_operator import Richardson

    params = inspect.signature(Richardson).parameters
    _require("preconditioner" not in params, f"Richardson takes a preconditioner again: {list(params)}")
    _require(all(p.kind is inspect.Parameter.KEYWORD_ONLY for p in params.values()),
             f"Richardson takes a positional argument: {list(params)}")
    M, _ = _split_ops()
    with pytest.raises(TypeError):
        Richardson(M.inverse())  # type: ignore[call-arg]  # the refusal under test
    with pytest.raises(TypeError):
        Richardson(preconditioner=M.inverse())  # type: ignore[call-arg]  # the refusal under test


def test_g9b_the_native_richardson_iterates_with_m_s_own_inverse():
    """A Richardson Green's preconditioner IS its ``implicit_inverse`` (the
    native representation: ``precondition`` is the identity map), and the
    driver it builds holds exactly that object and the realisation's corrector.
    Reds when the preconditioner is re-derived (a second inverse object) or
    the corrector is dropped at the driver."""
    from orpheus.numerics.green_operator import Richardson

    M, N = _split_ops()
    from tests.gates.numerics.test_iteration import MatrixOperator

    corrector = MatrixOperator(np.zeros((_n, _n)))
    g = GreenOperator(M, N, realisation=Richardson(corrector=corrector))
    _require(g.preconditioner is g.implicit_inverse, "the native preconditioner is not M's own inverse")
    drv = g.driver(max_iter=3, tol=1e-12)
    _require(drv.A_inv is g.preconditioner, "the driver iterates with another inverse")
    _require(drv.corrector is corrector, "the driver dropped the realisation's corrector")


def test_g9c_a_krylov_green_s_preconditioner_is_the_strategy_s():
    """GMRES's left preconditioner is free (it moves the trajectory, never the
    fixed point), so a Krylov realisation keeps the one handed in: the Green's
    ``preconditioner`` and its driver's are ``left_preconditioner``. Reds when a
    Krylov Green silently swaps in ``implicit_inverse``."""
    from orpheus.numerics.green_operator import Krylov

    M, N = _split_ops()
    left = 1.0 * M.inverse()
    g = GreenOperator(M, N, realisation=Krylov(left, restart=_n))
    _require(g.preconditioner is left, "the Green replaced the Strategy's left preconditioner")
    _require(g.driver(max_iter=3, tol=1e-12).preconditioner is left, "the GMRES driver does not use it")
