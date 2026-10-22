r"""Generator of ``s5_invariance_baseline.json`` — FROZEN at ``c27830f5``; never re-run after S5.

S5 of the numerics untangling (``scratch/numerics_untangling/packages/S5_design.md``,
gate 3; the spec is ``scratch/numerics_untangling/packages/S5_gate_spec.md``) deletes
``orpheus.numerics.invariance`` and re-spells every invariance verb as a composition
``relabelling_to(pushforward(support.action(g)))``.  This script records, on the tree
BEFORE that change, what the old kernel answered, so that the gates in
``tests/gates/numerics/test_invariance_composition.py`` can pin the new path against it
bit for bit.

⛔ **Do not re-run it after the carve.**  It imports the free functions S5 retires, so
it will not import; and a generator re-pointed at the new verbs would record the new
code against itself (``retirement-audit`` D.14, demotion).  The fixture is a RECORD
(what the code printed on 2026-10-19); the independent anchor beside it is the
Hungarian-matcher row in the gate module, which builds every permutation from raw
3-D rotation matrices and ``scipy.optimize.linear_sum_assignment``.

Population (each member rebuilt from a label at gate time, so a battery reaches it):

* every ``(registry spec, parameter set)`` that ``test_registry_rule_list._grid``
  sweeps (degree inversion at targets 0..21, plus that file's ``_EXTRA_GRID``),
  built through the registry's factory (the ``Quadrature`` classmethods);
* ``Quadrature.folded_product(n_mu, n_phi)`` for ``n_mu`` in {2, 4, 6, 8} and even
  ``n_phi`` in {4, 6, 8, 16, 32} (the tested range, by a census of ``tests/``);
* the polar orbits ``gauss_legendre_on_polar_orbit(n, axis)`` (S²/O2_a) for ``n`` in
  {1, 2, 3, 4, 8, 16} and every axis;
* three non-σ_y folds: ``lebedev(7)``, ``level_symmetric(4)`` and ``product(4, 8)``
  quotiented by each of ``Mirror(x|y|z)`` the rule admits;
* the chart μ-rule ``gauss_legendre_on_mu(n)`` (support ``Interval[-1,1]``), n in
  {2, 4, 8}: its verdicts are the ZERO-PADDING's, recorded so the spec can name the
  rows S5 turns into refusals — no gate pins them as answers.

Per rule: a fingerprint of the nodes and weights (so a moved rule is told apart from a
moved verdict), ``is_invariant_under`` for every candidate (finite and continuous), the
certificate's permutations (``indices`` of each element, in the realization's element
order) or ``None`` for every FINITE candidate, ``singular_set_under`` where a
certificate exists, ``symmetry_groups`` under both methods, ``ordinate_permutation``
for the 3-D x/y/z mirrors (``geometry.boundary._factors._mirror_motion``, the motion
every specular face uses) plus the pole deck ``SelfPairedDeck.mirror("x").motion`` and
one NON-normalising motion (a quarter turn about z), and on every product rule the
three ``_derived_product_group`` factor verdicts with their permutations.

Run: ``PYTHONPATH=<worktree> .venv/bin/python <this file>`` from outside the
repository (the script prints ``orpheus.__file__`` first).
"""

from __future__ import annotations

import hashlib
import json
import pathlib
import sys

import numpy as np

import orpheus
from orpheus.geometry.boundary import SelfPairedDeck
from orpheus.geometry.boundary._factors import _mirror_motion
from orpheus.geometry.transformation import RigidMotion
from orpheus.numerics import invariance as inv
from orpheus.numerics.quadrature import Quadrature, quadrature_registry
from orpheus.numerics.quadrature.rules_1d import (
    gauss_legendre_on_mu,
    gauss_legendre_on_polar_orbit,
)
from orpheus.numerics.quadrature.rules_circle import (
    NODE_ALIGNED,
    STAGGERED,
    periodic_trapezoid,
)
from orpheus.geometry.symmetry import SubgroupOfO3

COMMIT = "c27830f5"
OUT = pathlib.Path(__file__).with_name("s5_invariance_baseline.json")

TARGET_CAP = 21
EXTRA_GRID = {
    "LevelSymmetricSN": [{"sn_order": n} for n in (12, 16)],
    "GaussLegendre1D": [{"n_ordinates": n} for n in (1, 3, 5, 16, 24, 33, 64)],
}
FOLDED = [(m, p) for m in (2, 4, 6, 8) for p in (4, 6, 8, 16, 32)]
POLAR = [(n, a) for n in (1, 2, 3, 4, 8, 16) for a in ("x", "y", "z")]
MIRROR_FOLDS = [("lebedev", {"order": 7}), ("level_symmetric", {"sn_order": 4}),
                ("product", {"n_mu": 4, "n_phi": 8})]
CHART_MU = (2, 4, 8)

#: Finite groups asked of every rule beside its own ``candidate_groups``.
NAMED_FINITE = (
    [SubgroupOfO3.Trivial, SubgroupOfO3.OctahedralOh, SubgroupOfO3.IcosahedralIh]
    + [SubgroupOfO3.Mirror(a) for a in "xyz"]
    + [SubgroupOfO3.Cn(n) for n in range(2, 9)]
    + [SubgroupOfO3.Dnh(n) for n in range(1, 9)]
)
CONTINUOUS = (
    [SubgroupOfO3.Dinfh, SubgroupOfO3.SO3, SubgroupOfO3.O3]
    + [SubgroupOfO3.SO2(a) for a in "xyz"]
    + [SubgroupOfO3.O2(a) for a in "xyz"]
)


def grid(spec) -> list[dict]:
    seen: list[dict] = []
    for target in range(TARGET_CAP + 1):
        params = spec.degree_of_exactness_for(target)
        if params is not None and params not in seen:
            seen.append(params)
    for extra in EXTRA_GRID.get(spec.name, []):
        if extra not in seen:
            seen.append(extra)
    return seen


def fingerprint(measure) -> str:
    h = hashlib.sha256()
    h.update(np.ascontiguousarray(measure.nodes, dtype=float).tobytes())
    h.update(np.ascontiguousarray(measure.weights, dtype=float).tobytes())
    return h.hexdigest()[:16]


def perm_list(pi) -> list[int] | None:
    return None if pi is None else [int(i) for i in pi.indices]


def quarter_turn_z() -> RigidMotion:
    return RigidMotion.rotation(
        plane=(np.array([1.0, 0, 0]), np.array([0, 1.0, 0])), angle=np.pi / 2,
    )


def record_measure(measure) -> dict:
    groups: dict[str, SubgroupOfO3] = {}
    for g in list(inv.candidate_groups(measure)) + NAMED_FINITE + CONTINUOUS:
        groups.setdefault(g.name, g)
    invariant = {}
    certificates = {}
    singular = {}
    for name in sorted(groups):
        g = groups[name]
        invariant[name] = bool(measure.is_invariant_under(g))
        if g.realization.is_finite:
            cert = measure.certificate_under(g)
            certificates[name] = (
                None if cert is None
                else [[int(i) for i in p.indices] for p in cert.permutations]
            )
            if cert is not None:
                singular[name] = [int(i) for i in measure.singular_set_under(g)]
    return {
        "support": measure.support.name,
        "n_points": int(measure.n_points),
        "fingerprint": fingerprint(measure),
        "candidates": sorted(g.name for g in inv.candidate_groups(measure)),
        "is_invariant_under": invariant,
        "certificate_permutations": certificates,
        "singular_set_under": singular,
        "symmetry_groups_walk": sorted(g.name for g in measure.symmetry_groups(method="walk")),
        "symmetry_groups_bruteforce": sorted(
            g.name for g in measure.symmetry_groups(method="bruteforce")
        ),
        "permutation_under": {
            k: perm_list(measure.permutation_under(m.linear_part))
            for k, m in sorted(motions().items())
        },
    }


def motions() -> dict[str, RigidMotion]:
    out = {f"mirror_{a}": _mirror_motion(a, 3) for a in "xyz"}
    out["pole_deck_x"] = SelfPairedDeck.mirror(axis="x").motion
    out["quarter_turn_z"] = quarter_turn_z()
    return out


def record_quadrature(q: Quadrature) -> dict:
    row = record_measure(q.measure)
    row["ordinate_permutation"] = {
        k: perm_list(q.ordinate_permutation(m)) for k, m in sorted(motions().items())
    }
    return row


def product_factor_verdicts(polar, azimuthal, atol: float) -> dict:
    """The three factor checks ``rules_product._derived_product_group`` makes, verbatim."""
    n_phi = azimuthal.n_points
    e1, e2 = np.array([1.0, 0.0]), np.array([0.0, 1.0])
    step = RigidMotion.rotation(plane=(e1, e2), angle=2.0 * np.pi / n_phi)
    vertical_mirror = RigidMotion.reflection(normal=e2)
    polar_mirror = RigidMotion.reflection(normal=np.array([1.0]))
    rot = step.preserves(azimuthal.nodes, azimuthal.weights, atol=atol, weight_atol=atol)
    mir = vertical_mirror.preserves(
        azimuthal.nodes, azimuthal.weights, atol=atol, weight_atol=atol
    )
    pol = polar_mirror.preserves(
        polar.nodes.reshape(-1, 1), polar.weights, atol=atol, weight_atol=atol
    )
    return {
        "atol": atol,
        "rotation": perm_list(rot),
        "vertical_mirror": perm_list(mir),
        "polar_mirror": perm_list(pol),
    }


def main() -> None:
    print("orpheus from", orpheus.__file__, file=sys.stderr)
    rules: dict[str, dict] = {}
    specs = {s.name: s for s in quadrature_registry}
    for name in sorted(specs):
        for params in grid(specs[name]):
            label = f"registry|{name}|" + ",".join(f"{k}={v}" for k, v in sorted(params.items()))
            q = specs[name].factory(**params)
            rules[label] = record_quadrature(q)
    for n_mu, n_phi in FOLDED:
        rules[f"folded_product|n_mu={n_mu},n_phi={n_phi}"] = record_quadrature(
            Quadrature.folded_product(n_mu, n_phi)
        )
    for n, axis in POLAR:
        m = gauss_legendre_on_polar_orbit(n, axis)
        rules[f"polar_orbit|n={n},axis={axis}"] = record_measure(m)
    for builder, params in MIRROR_FOLDS:
        q = getattr(Quadrature, builder)(**params)
        for a in "xyz":
            g = SubgroupOfO3.Mirror(a)
            if q.measure.certificate_under(g) is None:
                continue
            label = f"mirror_fold|{builder}|" + ",".join(
                f"{k}={v}" for k, v in sorted(params.items())
            ) + f"|by=sigma_{a}"
            # The MEASURE's fold: a mu-mirror fold of a level-structured rule is
            # refused by the level structure (an orbit crosses a level), not by
            # the measure, so the measure-level quotient is the shipped object here.
            rules[label] = record_measure(q.measure.quotient(g))
    for n in CHART_MU:
        rules[f"chart_mu|n={n}"] = record_measure(gauss_legendre_on_mu(n))

    # The product rules' factor checks — every product rule in the population.
    factors: dict[str, dict] = {}
    atol = _product_atol()
    product_params = grid(specs["ProductQuadrature"])
    for p in product_params:
        n_mu, n_phi = p["n_mu"], p["n_phi"]
        factors[f"product|n_mu={n_mu},n_phi={n_phi}"] = product_factor_verdicts(
            gauss_legendre_on_mu(n_mu), periodic_trapezoid(n_phi, shift=NODE_ALIGNED), atol,
        )
    for n_mu, n_phi in FOLDED:
        factors[f"folded_product_parent|n_mu={n_mu},n_phi={n_phi}"] = product_factor_verdicts(
            gauss_legendre_on_mu(n_mu), periodic_trapezoid(n_phi, shift=STAGGERED), atol,
        )

    payload = {
        "commit": COMMIT,
        "generator": "tests/gates/numerics/data/generate_s5_invariance_baseline.py",
        "weight_atol": inv.WEIGHT_ATOL,
        "node_window_factor": inv._NODE_WINDOW_FACTOR,
        "n_rules": len(rules),
        "n_product_factor_rows": len(factors),
        "rules": rules,
        "product_factors": factors,
    }
    OUT.write_text(json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n")
    n_cert = sum(len(r["certificate_permutations"]) for r in rules.values())
    n_inv = sum(len(r["is_invariant_under"]) for r in rules.values())
    n_true = sum(sum(r["is_invariant_under"].values()) for r in rules.values())
    n_cert_some = sum(
        sum(v is not None for v in r["certificate_permutations"].values()) for r in rules.values()
    )
    print(
        f"rules={len(rules)} verdict_rows={n_inv} (True {n_true}) "
        f"finite_rows={n_cert} (certificate {n_cert_some}) "
        f"product_factor_rows={len(factors)} -> {OUT}",
        file=sys.stderr,
    )


def _product_atol() -> float:
    """The ``atol`` ``spherical_product`` hands ``_derived_product_group`` by default."""
    import inspect

    from orpheus.numerics.quadrature.rules_product import spherical_product

    return float(inspect.signature(spherical_product).parameters["atol"].default)


if __name__ == "__main__":
    main()
