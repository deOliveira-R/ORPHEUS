"""Dense oracles for the near-kernel witness rows (independent of the witness: the loss assembled through the production matvec, then numpy)."""
from __future__ import annotations

import numpy as np


def residual(span: np.ndarray, vecs: np.ndarray) -> float:
    if vecs.shape[1] == 0:
        return 0.0
    if span.shape[1] == 0:
        return 1.0
    Q, sv, _ = np.linalg.svd(span, full_matrices=False)
    Q = Q[:, sv > 1e-10 * sv[0]]
    return float(np.linalg.norm(vecs - Q @ (Q.T @ vecs)) / np.linalg.norm(vecs))


def dense_kernels(A: np.ndarray, g: np.ndarray):
    """(s, right kernel R, left kernel in .H coordinates G^+ L, k modulo ker G)."""
    U, s, Vt = np.linalg.svd(A)
    n = A.shape[0]
    raw = int(np.sum(s < 1e-10 * s[0]))
    R, L = Vt[n - raw:].T, U[:, n - raw:]
    pos = g > 0.0
    k = int(np.sum(np.linalg.svd(R[pos], compute_uv=False) > 1e-8)) if raw else 0
    gp = np.divide(1.0, g, out=np.zeros_like(g), where=pos)
    return s, R, gp[:, None] * L, k


def derived_containment_tol(stop: float, s: np.ndarray, k: int) -> float:
    """Davis-Kahan / Wedin: sin(theta) <= r sigma_max / sigma_{k+1}; asserted with a factor 10."""
    n = s.size
    return 10.0 * stop * s[0] / s[n - k - 1]


def metric_pseudo_inverse_rcond(A: np.ndarray, g: np.ndarray, R: np.ndarray, LH: np.ndarray) -> float:
    """1 / (||A||_1 ||S||_1), S = (I - Pi_G) A^+ (I - Pi^H_G): the G-metric pseudo-inverse (the ruled quantity)."""
    n = A.shape[0]
    U, s, Vt = np.linalg.svd(A)
    k = R.shape[1]
    pinv = Vt[: n - k].T @ np.diag(1.0 / s[: n - k]) @ U[:, : n - k].T
    G = np.diag(g)
    Pg = np.eye(n) - (R @ np.linalg.solve(R.T @ G @ R, R.T @ G) if k else 0.0)
    Ph = np.eye(n) - (LH @ np.linalg.solve(LH.T @ G @ LH, LH.T @ G) if k else 0.0)
    S = Pg @ pinv @ Ph
    return 1.0 / (np.abs(A).sum(axis=0).max() * np.abs(S).sum(axis=0).max())
