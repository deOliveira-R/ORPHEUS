"""The white-BC slab flux rebuilt from the Peierls equation and the partial-
current balance, every volume integral by mpmath.quad (no antiderivative)."""
import time, mpmath
from orpheus.derivations.continuous.peierls_nystrom.reference import slab_uniform_source_white_bc_analytical as f
import importlib, sys
sys.path.insert(0, "scratch/reference_architecture/p0probe")
from err032_arm import phi_wrong
def phi_quad(x, L, s, dps=30):
    with mpmath.workdps(dps):
        x, L, s = map(mpmath.mpf, (x, L, s)); tL = s*L
        jplus_vol = mpmath.quad(lambda xp: mpmath.expint(2, s*(L-xp))/2, [0, L])   # outgoing at L, uniform S=1
        jminus = jplus_vol / (1 - 2*mpmath.expint(3, tL))                            # Mark closure J- = J+
        pts = [0, x, L] if 0 < x < L else [0, L]
        vol = mpmath.quad(lambda xp: mpmath.expint(1, s*abs(x-xp))/2, pts)
        return vol + 2*jminus*(mpmath.expint(2, s*x) + mpmath.expint(2, s*(L-x)))
t0=time.perf_counter(); worst=0; worstw=1e9
for L, s in [(0.1,1.0),(1.0,1.0),(5.0,2.0),(100.0,0.5)]:
    for fr in (0.0,0.2,0.5,0.8,1.0):
        q = phi_quad(fr*L, L, s)
        r = abs(q - f(fr*L, L, s, dps=30))/abs(q); worst=max(worst, r)
        rw = abs(q - phi_wrong(fr*L, L, s, dps=30))/abs(q); worstw=min(worstw, rw)
print(f"max rel honest {float(worst):.2e}  min rel mutant {float(worstw):.2e}  {time.perf_counter()-t0:.2f}s")
