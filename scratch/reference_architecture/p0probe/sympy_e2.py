import time, sympy as sp
t0=time.perf_counter()
u,tau=sp.symbols("u tau", positive=True)
I = sp.integrate(sp.expint(2,u),(u,0,tau))
print("integrate ->", I, f"{time.perf_counter()-t0:.2f}s")
d = sp.simplify(I - (sp.Rational(1,2) - sp.expint(3,tau)))
print("diff vs 1/2-E3:", d)
dw = sp.simplify(I - (1 - sp.expint(3,tau)))
print("diff vs 1-E3 (the ERR-032 identity):", dw)
# derivative route independent of integrate(): d/dtau(1/2 - E3) == E2 and value at 0
print("deriv check:", sp.simplify(sp.diff(sp.Rational(1,2)-sp.expint(3,tau),tau) - sp.expint(2,tau)),
      "at0:", sp.limit(sp.Rational(1,2)-sp.expint(3,tau), tau, 0), "wrong at0:", sp.limit(1-sp.expint(3,tau),tau,0))
print(f"{time.perf_counter()-t0:.2f}s")
