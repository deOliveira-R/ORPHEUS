r"""B3.4a mutation harness — an in-process pytest plugin proving gate teeth.

Sibling of ``scratch/b3_2_mutations.py``, same contract. Usage (SERIAL,
canonical ``-O``)::

    ORPHEUS_B34A=W1 .venv/bin/python -O -m pytest \
        tests/sn/operators/test_angular_average_operator.py \
        tests/sn/operators/test_sn_boundary_realizer.py \
        tests/sn/operators/test_capability_survival.py \
        tests/geometry/test_bc_universal_invariants.py \
        -p no:randomly -p scratch.b3_4a_mutations -q

``ORPHEUS_B34A`` unset ⇒ the CONTROL leg: nothing is patched, everything must
be GREEN. Any other value installs the named mutation by monkeypatching
production IN PROCESS — no file on disk is touched, so this can never destroy
uncommitted work (``process-discipline``: never ``git checkout`` a file you
hold edits in).

Each mutation is the plausible transcription of a real B3.4a hazard, not an
arbitrary perturbation:

======  =====================================================================
id      the wrong code it simulates
======  =====================================================================
W1      ``from_quadrature`` classifies Γ₊ with the RETIRED private
        ``(outward_sign * mu) > 0.0`` test instead of ``TANGENTIAL_EPS`` —
        the twin the narrowing dissolved. Silent on every quadrature whose
        tangential cosines are exactly 0.0; wrong on ``product(2, 4)``
W2      the normalization drops ``|Ω·n̂|``: ``norm = Σ w`` instead of
        ``Σ w |Ω·n̂|`` (the ERR-004 family — a quadrature-dependent constant
        standing in for a computed one)
W3      the codomain is sized from Γ₊ instead of Γ₋. VALUE-INVISIBLE by
        construction: ``|Γ₊| == |Γ₋|`` on every reachable face (vv Mode 12).
        Shipped as the DESIGNED-GREEN control — if any gate reds on it, the
        gate is reading something other than what it claims
W4      ``apply`` broadcasts over the INPUT's leading axis instead of the
        codomain — the codomain echoes the domain
W5      ``apply`` drops its domain shape guard, so a caller that has not
        migrated can still hand it the full face slot
W6      the ``cos_w`` positivity guard reverts to the pre-B3.4a ``>= 0``,
        re-admitting a masked full-face array
W7      the normalization is face-sign dependent (``×0.9`` on lower faces) —
        a mutation that breaks RECIPROCITY while leaving the ``+1``-face
        conservation legs exactly green, i.e. the discriminator proving the
        reciprocity gate carries teeth the conservation gate does not
Q1      ``IncomingSourceOperator.apply`` sizes ``q`` from the INPUT rather
        than from ``|Γ₋|`` — the ERR-047 full-face shape returning
Q2      the ``n_inflow`` non-negativity guard is dropped
R1      the white realizer skips ``_checked_angular_average`` — no
        orientation cross-check (the pre-B3.4a state)
R2      the white realizer skips ``_outflow_restriction`` — a faceless method
        space realizes white again
R3      the prescribed-inflow realizer skips ``_outflow_restriction``
======  =====================================================================

A mutation is only evidence if it BIT (``vv`` Mode-8 METHOD WARNING). Every
run prints a fingerprint — a value hash plus the constructed half-trace sizes
on the discriminating ``product(2, 4)`` fixture — so a no-op patch is visible
in the run's first line rather than being read as "the gate is blind". ``W3``
is the deliberate exception: its fingerprint MUST show the size change while
the value hash stays put.
"""
from __future__ import annotations

import hashlib
import os

import numpy as np

_MUT = os.environ.get("ORPHEUS_B34A", "").strip()


# ─────────────────────────────────────────────────────────────────────
# AngularAverageOperator
# ─────────────────────────────────────────────────────────────────────


def _mutated_from_quadrature(mut: str):
    """``from_quadrature`` carrying mutation ``mut`` (W1 / W2 / W3 / W7).

    Transcribed from ``orpheus/sn/boundary/angular.py`` at B3.4a so each
    mutation is a one-line difference; the mutation points are marked ``# MUT``.
    """
    from orpheus.numerics.face_layout import AXIS_NAMES
    from orpheus.numerics.spaces.angular_trace_space import (
        TANGENTIAL_EPS,
        build_omega_dot_n,
    )

    def from_quadrature(cls, quadrature, axis, outward_sign):
        if axis not in AXIS_NAMES:
            raise ValueError(
                f"Unknown axis: {axis!r}; expected one of {AXIS_NAMES}"
            )
        if outward_sign not in (+1, -1):
            raise ValueError(
                f"outward_sign must be +1 or -1; got {outward_sign}"
            )
        face = f"{axis}{'max' if outward_sign == +1 else 'min'}"
        omega_dot_n = build_omega_dot_n(quadrature, (face,))[0]
        if mut == "W1":
            outflow = np.flatnonzero(omega_dot_n > 0.0)      # MUT retired test
        else:
            outflow = np.flatnonzero(omega_dot_n > +TANGENTIAL_EPS)
        inflow = np.flatnonzero(omega_dot_n < -TANGENTIAL_EPS)
        if outflow.size == 0:
            raise ValueError(
                f"AngularAverageOperator.from_quadrature: no outgoing "
                f"ordinates on axis={axis!r}, outward_sign={outward_sign}. "
                f"Quadrature is degenerate for this face."
            )
        if inflow.size == 0:
            raise ValueError(
                f"AngularAverageOperator.from_quadrature: no incoming "
                f"ordinates on axis={axis!r}, outward_sign={outward_sign} "
                f"— the re-emitted flux has no Γ₋ to land on. Quadrature "
                f"is degenerate for this face."
            )
        w = np.asarray(quadrature.weights, dtype=float)
        cos_w = w[outflow] * omega_dot_n[outflow]
        if mut == "W2":
            norm = float(w[outflow].sum())                   # MUT drops |Ω·n̂|
        elif mut == "W7":
            norm = float(cos_w.sum()) * (1.0 if outward_sign == +1 else 0.9)  # MUT
        else:
            norm = float(cos_w.sum())
        n_inflow = outflow.size if mut == "W3" else inflow.size  # MUT (W3)
        return cls(cos_w=cos_w, norm=norm, n_inflow=int(n_inflow))

    return classmethod(from_quadrature)


def _mutated_apply(mut: str):
    """``AngularAverageOperator.apply`` carrying W4 / W5."""

    def apply(self, psi):
        psi = np.asarray(psi)
        if mut != "W5" and psi.shape[0] != self.n_outflow:  # MUT (W5 drops it)
            raise ValueError(
                f"AngularAverageOperator.apply: psi.shape[0] = "
                f"{psi.shape[0]}, expected |Γ₊| = {self.n_outflow}."
            )
        cos_w = self._cos_w
        if mut == "W5" and psi.shape[0] != cos_w.size:
            # A dropped guard does not make the contraction legal; the
            # pre-B3.4a code would have used a full-face cos_w. Emulate the
            # un-guarded read by truncating/padding to the input's length.
            cos_w = np.resize(cos_w, psi.shape[0])
        psi_avg = (
            cos_w.reshape((-1,) + (1,) * (psi.ndim - 1)) * psi
        ).sum(axis=0) / self._norm
        n_rows = psi.shape[0] if mut == "W4" else self.n_inflow  # MUT (W4)
        return np.broadcast_to(
            psi_avg[None, ...], (n_rows,) + psi.shape[1:],
        ).copy()

    return apply


def _mutated_init():
    """W6 — the pre-B3.4a ``>= 0`` positivity guard."""

    def __init__(self, cos_w, norm, n_inflow):
        cos_w = np.asarray(cos_w, dtype=float)
        if cos_w.ndim != 1:
            raise ValueError("cos_w must be 1-D")
        if (cos_w < 0).any():                                # MUT was (<= 0)
            raise ValueError("cos_w must be non-negative")
        norm = float(norm)
        if norm <= 0.0:
            raise ValueError("AngularAverageOperator norm must be positive")
        n_inflow = int(n_inflow)
        if n_inflow <= 0:
            raise ValueError(
                "AngularAverageOperator n_inflow must be positive"
            )
        self._cos_w = cos_w.copy()
        self._norm = norm
        self.n_outflow = int(cos_w.size)
        self.n_inflow = n_inflow

    return __init__


# ─────────────────────────────────────────────────────────────────────
# IncomingSourceOperator
# ─────────────────────────────────────────────────────────────────────


def _mutated_source_apply():
    """Q1 — the delivered q is sized from the INPUT, not from Γ₋."""

    def apply(self, psi_out):
        psi_out = np.asarray(psi_out)
        return self.source.evaluate(tuple(int(s) for s in psi_out.shape))  # MUT

    return apply


def _mutated_source_init():
    """Q2 — the codomain-size guard is dropped."""

    def __init__(self, source, *, n_inflow):
        self.source = source
        self.n_inflow = int(n_inflow)                        # MUT no guard

    return __init__


# ─────────────────────────────────────────────────────────────────────
# SNBoundaryRealizer
# ─────────────────────────────────────────────────────────────────────


def _mutate_realizer(mut: str) -> None:
    from orpheus.sn.boundary import realizer as _r

    if mut == "R1":
        # MUT: no orientation cross-check — build straight from the law.
        # Signature updated at B3.4b: the helper now takes the declared
        # orientation as primitives (so albedo+IsotropicReturn can reach the
        # same body) instead of a WhiteBoundary.
        _r._checked_angular_average = (
            lambda quadrature, method_space, gamma_out, *,
            axis, outward_sign, law_key: (
                _r.AngularAverageOperator.from_quadrature(
                    quadrature, axis, outward_sign,
                )
            )
        )
    elif mut in {"R2", "R3"}:
        # MUT: the domain guard is skipped for the named law, so a faceless
        # method space realizes it again (the pre-narrowing state).
        law_key = "white" if mut == "R2" else "prescribed_inflow"
        original = _r._outflow_restriction

        def patched(method_space, key, _orig=original, _skip=law_key):
            if key == _skip and method_space.outflow_indices is None:
                from orpheus.numerics.operator import TraceRestrictionOperator

                return TraceRestrictionOperator(
                    np.arange(method_space.quadrature.N, dtype=np.intp),
                    n_total=method_space.quadrature.N, axis=0,
                )
            return _orig(method_space, key)

        _r._outflow_restriction = patched


# ─────────────────────────────────────────────────────────────────────
# Plugin hook + bite check
# ─────────────────────────────────────────────────────────────────────


def pytest_configure(config) -> None:  # noqa: D401 - pytest hook
    if not _MUT:
        print("\n[b3_4a_mutations] CONTROL leg — no mutation installed.")
        return
    from orpheus.sn.boundary.angular import (
        AngularAverageOperator,
        IncomingSourceOperator,
    )

    if _MUT in {"W1", "W2", "W3", "W7"}:
        AngularAverageOperator.from_quadrature = _mutated_from_quadrature(_MUT)
    elif _MUT in {"W4", "W5"}:
        AngularAverageOperator.apply = _mutated_apply(_MUT)
    elif _MUT == "W6":
        AngularAverageOperator.__init__ = _mutated_init()
    elif _MUT == "Q1":
        IncomingSourceOperator.apply = _mutated_source_apply()
    elif _MUT == "Q2":
        IncomingSourceOperator.__init__ = _mutated_source_init()
    elif _MUT in {"R1", "R2", "R3"}:
        _mutate_realizer(_MUT)
    else:
        raise ValueError(f"unknown ORPHEUS_B34A variant {_MUT!r}")
    print(f"\n[b3_4a_mutations] variant={_MUT} installed; bite check: "
          f"{_bite_fingerprint()}")


def _bite_fingerprint() -> str:
    """Positive control — prove the patch CHANGED production.

    Reports BOTH a value hash and the constructed half-trace sizes on the
    discriminating ``product(2, 4)`` xmax fixture. W3 is expected to move the
    SIZES only; W1 moves the sizes AND the values; W2/W7 move the values only.
    A variant that moves neither did not install.
    """
    try:
        from orpheus.numerics.quadrature import Quadrature
        from orpheus.sn.boundary.angular import (
            AngularAverageOperator,
            IncomingSourceOperator,
        )
        from orpheus.geometry.boundary import ConstantInflowSource

        parts = []
        quad = Quadrature.product(n_mu=2, n_phi=4)
        try:
            op = AngularAverageOperator.from_quadrature(quad, "x", +1)
            probe = np.arange(1, op.n_outflow + 1, dtype=float)
            out = op.apply(probe)
            parts.append(
                f"product(2,4) dom/cod={op.n_outflow}/{op.n_inflow} "
                f"sha={hashlib.sha256(out.tobytes()).hexdigest()[:16]}"
            )
        except Exception as exc:
            parts.append(f"<AAO: {type(exc).__name__}: {exc}>")
        leb = Quadrature.lebedev(17)
        try:
            op2 = AngularAverageOperator.from_quadrature(leb, "x", +1)
            rng = np.random.default_rng(5)
            out2 = op2.apply(rng.uniform(0.5, 1.5, size=op2.n_outflow))
            parts.append(
                f"lebedev17 dom/cod={op2.n_outflow}/{op2.n_inflow} "
                f"sha={hashlib.sha256(out2.tobytes()).hexdigest()[:16]}"
            )
        except Exception as exc:
            parts.append(f"<AAO-leb: {type(exc).__name__}: {exc}>")
        try:
            q = IncomingSourceOperator(ConstantInflowSource(value=2.0), n_inflow=3)
            parts.append(f"q rows={q.apply(np.ones((7, 2))).shape}")
        except Exception as exc:
            parts.append(f"<ISO: {type(exc).__name__}: {exc}>")
        return " | ".join(parts)
    except Exception as exc:  # pragma: no cover - diagnostic only
        return f"<bite check unavailable: {type(exc).__name__}: {exc}>"
