"""-p plugin: ERR-032 re-introduced — the white-BC slab closed form replaced
by the catalogue's phi_wrong, rebound in every sys.modules binding."""
import sys, mpmath, pytest
TARGET = "orpheus.derivations.continuous.peierls_nystrom.reference"
def phi_wrong(x, L, sig_t, *, dps=50):
    with mpmath.workdps(dps):
        s = mpmath.mpf(sig_t); tL = s*L
        beta = (1 - mpmath.expint(3, tL)) / (1 - 2*mpmath.expint(3, tL))
        return (1/(2*s)) * (2 + (2*beta - 1)*(mpmath.expint(2, s*x) + mpmath.expint(2, s*(L - x))))
N = [0]
def pytest_configure(config):
    import importlib
    mod = importlib.import_module(TARGET); orig = mod.slab_uniform_source_white_bc_analytical
    honest = float(orig(0.3, 1.0, 1.0)); mutant = float(phi_wrong(0.3, 1.0, 1.0))
    if honest == mutant: raise pytest.UsageError("Uninstallable: mutant equals honest")
    config._err032 = (orig, honest, mutant)
def _sweep(orig):
    for m in list(sys.modules.values()):
        d = getattr(m, "__dict__", None)
        if isinstance(d, dict):
            for k, v in list(d.items()):
                if v is orig: setattr(m, k, phi_wrong); N[0] += 1
def pytest_collection_finish(session):
    _sweep(session.config._err032[0])
    o, h, m = session.config._err032
    print(f"\n[err032_arm] rebinds={N[0]} honest={h!r} mutant={m!r}", file=sys.stderr)
    if N[0] < 2: raise pytest.UsageError(f"Uninstallable: {N[0]} rebinds")
