"""Shared fixtures for the singular-loss gates (space-membership U0, the Fredholm verdict). Values `[M]` 2026-10-31 on the
membership worktree, by the test-architect (``scratch/membership_covector/u0/probes/_fx.py`` is the probe twin)."""
from __future__ import annotations

import numpy as np

from orpheus.data.macro_xs.synthetic import make_mixture
from orpheus.geometry import StructuredGeometry
from orpheus.geometry.deck import MirrorDeck
from orpheus.mesh import CellsByCount, Mesher


def mix(sig_t, sig_c, sig_s, sig_f=None, nu=None, chi=None):
    ng = len(sig_t)
    z = np.zeros(ng)
    return make_mixture(
        sig_t=np.asarray(sig_t, float), sig_c=np.asarray(sig_c, float),
        sig_f=z if sig_f is None else np.asarray(sig_f, float), nu=z if nu is None else np.asarray(nu, float),
        chi=z if chi is None else np.asarray(chi, float), sig_s=np.asarray(sig_s, float),
    )


PURE_1G = mix([0.5], [0.0], [[0.5]])
BOTH_2G = mix([0.5, 0.5], [0.0, 0.0], [[0.4, 0.1], [0.0, 0.5]])
G2_ONLY = mix([0.5, 0.5], [0.1, 0.0], [[0.3, 0.1], [0.0, 0.5]])
G2_UPSCATTER = mix([0.5, 0.5], [0.1, 0.0], [[0.3, 0.1], [0.05, 0.45]])
NEAR_1E14 = mix([0.5], [1e-14], [[0.5 - 1e-14]])
NEAR_1E8 = mix([0.5], [1e-8], [[0.5 - 1e-8]])
NEAR_2G_1E10 = mix([0.5, 0.5], [0.1, 1e-10], [[0.3, 0.1], [0.0, 0.5 - 1e-10]])   # homogeneous rcond: 2g (1g is 1x1, rcond = 1)
FISSILE_NO_CAPTURE = mix([0.5], [0.0], [[0.45]], sig_f=[0.05], nu=[2.5], chi=[1.0])
ABSORBER = mix([0.5], [0.1], [[0.4]])
#: fission in g1 emits into g1, g1 downscatters into the absorption-free closed class {g2}: <w, chi> = 0.5 (w = (0.5, 1)).
FIS_FEEDS = mix([0.5, 0.5], [0.05, 0.0], [[0.3, 0.1], [0.0, 0.5]], sig_f=[0.05, 0.0], nu=[2.5, 0.0], chi=[1.0, 0.0])
#: g2 is an absorption-free class nobody feeds: <w, chi> = 0 (w = (0, 1)); k_inf = nu Sf1 / (St1 - Ss11) = 0.125/0.1.
FIS_ISOLATED = mix([0.5, 0.5], [0.05, 0.0], [[0.4, 0.0], [0.0, 0.5]], sig_f=[0.05, 0.0], nu=[2.5, 0.0], chi=[1.0, 0.0])
FIS_ISOLATED_KINF = 0.125 / 0.1
FIS_REGULAR = mix([0.5, 0.5], [0.05, 0.02], [[0.3, 0.1], [0.0, 0.48]], sig_f=[0.05, 0.0], nu=[2.5, 0.0], chi=[1.0, 0.0])


def slab(left=None, right=None, bps=(0.0, 10.0), mats=(0,), n=(8,)):
    left = MirrorDeck() if left is None else left
    right = MirrorDeck() if right is None else right
    g = StructuredGeometry.slab(bps, mats, left=left, right=right)
    return Mesher(g).partition(tuple(CellsByCount.uniform_width(k) for k in n)).mesh


#: The gauge row's mesh: two zones of one material, 3 + 5 cells, so the bulk metric is NON-uniform (1.0 and 1.4).
#: On the uniform 8-cell mesh an antisymmetric balanced source makes the Euclidean and G projections coincide [M].
TWO_ZONE = dict(bps=(0.0, 3.0, 10.0), mats=(0, 0), n=(3, 5))


def dense_of(apply, n):
    A = np.empty((n, n))
    e = np.zeros(n)
    for j in range(n):
        e[j] = 1.0
        A[:, j] = apply(e)
        e[j] = 0.0
    return A


def null(A, rel=1e-10):
    U, s, Vt = np.linalg.svd(A)
    k = int(np.sum(s < rel * s[0]))
    n = A.shape[0]
    return U[:, n - k:], Vt[n - k:].T


def bordered_member(A, q, right, constraint_rows):
    """The member of {x : A x = q} with constraint_rows @ x = 0, by one dense bordered solve (the oracle; it never asks
    the SUT's projector). For the minimum-G-norm member, constraint_rows = (G * right).T."""
    k = right.shape[1]
    B = np.block([[A, right], [constraint_rows, np.zeros((k, k))]])
    return np.linalg.solve(B, np.concatenate([q, np.zeros(k)]))[: A.shape[0]]
