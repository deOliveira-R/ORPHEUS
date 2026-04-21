---
name: CP/MOC moment-integral recursions — canonical references
description: Where to cite J_k(z)=∫u^k E_1(u)du (slab) and T_k(n,z)=∫u^k Ki_n(u)du (cylinder). Carlvik 1965/1967 do NOT have them; Hébert 2009/2020 ch.3 and Stamm'ler-Abbate ch.4-6 are the textbook homes.
type: reference
---

**Recursion 1 — slab Peierls polynomial-source moments**:
J_k(z) = ∫₀^z u^k E_1(u) du = z^(k+1)/(k+1) · E_1(z) + γ(k+1, z)/(k+1).
This is integration by parts on E_1; it is **not** a transport-physics
result, it is just the antiderivative of u^k E_1(u). The cleanest
canonical citations are:
- **Abramowitz & Stegun (1964) §5.1.32 and §5.1.36** — recursion for
  ∫E_n and the relation E_n(z) = z^(n-1) Γ(1-n, z)/Γ(n) which
  collapses to the γ(k+1,z) form above when integrated by parts.
- **Lewis & Miller (1984/1993) "Computational Methods of Neutron
  Transport", App. C** — restates the slab kernel moments in the
  exact form used for higher-order Peierls / MOC sources.
- **Hébert (2009/2020) "Applied Reactor Physics" §3.2-3.3** — derives
  the polynomial-source slab CP from this recursion.

**Recursion 2 — cylindrical Bickley-moment integrals**:
T_k(1, z) = ∫₀^z u^k Ki_1(u) du = -Σ_{m=0}^k k!/(k-m)! · z^(k-m) · Ki_(2+m)(z) + k! · Ki_(k+2)(0).
Derivation is mechanical: Ki_n' = -Ki_(n-1), ∫Ki_n du = -Ki_(n+1),
plus k integrations by parts. Closed-form values at zero use
Ki_n(0) = (√π/2) Γ(n/2)/Γ((n+1)/2) — see Bickley & Nayler (1935)
Phil. Mag. 20:343 (the original definition paper, also recapped in
Carlvik 1965 eq. 3.1.2). The recursion itself is **not in Carlvik
1965 or 1967** (both stop at flat-source Ki₃ kernels). The textbook
homes are:
- **Stamm'ler & Abbate (1983) "Methods of Steady-State Reactor
  Physics in Nuclear Design", Ch. IV** (pp. 105-141) — the standard
  book-length treatment of cylindrical CP integration. Tightened
  from the earlier "Ch. 4-6" range: Ch.V is P_L multi-D, Ch.VI is
  the SN method — neither is CP. **Stamm'ler Ch.IV uses only
  flat-source (rank-0) CP**; the Bickley-moment recursion for
  higher-order polynomial sources is NOT in Stamm'ler either.
  Primary home for higher-order cylindrical moments remains Hébert
  §3.4-3.5 and Abramowitz-Stegun §5.1 (slab).
- **Hébert (2009/2020) "Applied Reactor Physics" §3.4-3.5** — the
  modern restatement; the 3rd ed. is DOI 10.1515/9782553017445 (not
  open-access, but standard library holding).

**What is NOT a primary source**:
- Carlvik 1965 (Geneva A/CONF.28/P/681) — flat-source CP only.
- Carlvik 1967 (AE-279) — flat-source CP + DIT; no Ki_n moment table.
- Sanchez & McCormick 1982 NSE 80:481 — review article, gives the
  *form* of polynomial-flux multifunction expansion (Sec III) but
  does NOT tabulate the closed-form moment integrals; it cites
  primary sources instead.

**How to apply**: when implementing high-order Peierls/CP/MOC moment
solvers, cite Abramowitz-Stegun (slab) + Hébert §3.4 or
Stamm'ler-Abbate Ch. 5 (cylinder). Do NOT cite Carlvik 1965 or
Sanchez-McCormick 1982 for the closed-form recursions — they only
cover the flat-source kernel and the formulation framework
respectively.
