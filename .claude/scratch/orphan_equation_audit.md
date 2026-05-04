# Orphan Equation Inventory

Generated from `docs/verification/matrix.rst` § Orphan equations.
Each entry shows: theory page, line, section, and the equation block (directive + label + body until blank line).

## `c-in-remapping`

- File: `docs/theory/peierls_nystrom.rst`
- Line: 3111
- Section: 'The structural obstruction — the :math:`c_{\\rm in}` remapping'

```rst
.. math::
   :label: c-in-remapping

   c_{\rm in}(\mu_{\rm emit}) \;=\;
      \sqrt{1 - \left(\frac{R}{r_0}\right)^{\!2}\!
         \bigl(1 - \mu_{\rm emit}^2\bigr)},
```

## `singular-eigenfunction-eq40`

- File: `docs/theory/singular_eigenfunction.rst`
- Line: 488
- Section: 'X-function (Atalay Eq 40)'

```rst
.. math::
   :label: singular-eigenfunction-eq40

   X(\mu) = \exp\!\Bigg\{ -\frac{c}{2} \int_0^1 d\nu\, g_1(c,\nu)\,
       \Big[d^2(\nu^2)\Big(1 + \frac{c\nu^2}{1-\nu^2}\Big)
            + 3 f_1 (1-c)^2 \nu^2 d(-\nu^2)\Big]\,\ln(\nu - \mu) \Bigg\}
```

## `singular-eigenfunction-eq42`

- File: `docs/theory/singular_eigenfunction.rst`
- Line: 508
- Section: 'Extrapolated endpoint :math:`z_0` (Atalay Eq 42)'

```rst
.. math::
   :label: singular-eigenfunction-eq42

   z_0 = -\frac{\nu_0}{2} \ln\!\frac{d(-\nu_0\bar\nu)}{d(\nu_0\bar\nu)}
        + \frac{c\,\nu_0}{4} \int_0^1 d\mu\, g_1(c,\mu)\,
                                \Big[d^2(\mu^2)\Big(1 + \frac{c\mu^2}{1-\mu^2}\Big)
                                     + 3 f_1 (1-c)^2 \mu^2 d(-\mu^2)\Big]\,
                                \ln\!\frac{\nu_0 + \mu}{\nu_0 - \mu}
```

## `singular-eigenfunction-eq46`

- File: `docs/theory/singular_eigenfunction.rst`
- Line: 209
- Section: 'Slab criticality (Eq 46)'

```rst
.. math::
   :label: singular-eigenfunction-eq46

   \pm\frac{\pi}{2} - \arctan\!\frac{R \sin[(d-z_0)/|\nu_0|] + \sin[(d+z_0)/|\nu_0|]}
                                       {R \cos[(d-z_0)/|\nu_0|] - \cos[(d+z_0)/|\nu_0|]}
   = \arctan\!\frac{\big(K_0 \bar\nu - 3 f_1 (1-c) \bar\nu
                       [K_1 \bar\nu d(\nu_0^2) - K_0 \nu_0^2 d(\bar\nu^2)]\big) |\nu_0|}
                  {(1+K_2) d(\nu_0 \bar\nu) d(-\nu_0 \bar\nu)
                   + K_1 \bar\nu d(\nu_0^2) - K_0 \nu_0^2 d(\bar\nu^2)}
```

## `singular-eigenfunction-eq5`

- File: `docs/theory/singular_eigenfunction.rst`
- Line: 534
- Section: 'Validity bound (Atalay Eq 5)'

```rst
.. math::
   :label: singular-eigenfunction-eq5

   c \le 1 + \frac{1}{3 f_1}
```

## `singular-eigenfunction-eq54`

- File: `docs/theory/singular_eigenfunction.rst`
- Line: 244
- Section: 'Sphere criticality (Eq 54) via parity flip'

```rst
.. math::
   :label: singular-eigenfunction-eq54

   \arctan\!\frac{\sin[(d+z_0)/|\nu_0|] - R \sin[(d-z_0)/|\nu_0|]}
                {\cos[(d+z_0)/|\nu_0|] + R \cos[(d-z_0)/|\nu_0|]}
   = \arctan\!\frac{\big(L_0 \bar\nu - 3 f_1 (1-c) \bar\nu
                       [L_1 \bar\nu d(\nu_0^2) - L_0 \nu_0^2 d(\bar\nu^2)]\big) |\nu_0|}
                  {(1+L_2) d(\nu_0 \bar\nu) d(-\nu_0 \bar\nu)
                   + L_1 \bar\nu d(\nu_0^2) - L_0 \nu_0^2 d(\bar\nu^2)}
```

## `e1-decomposition`

- File: `docs/theory/collision_probability.rst`
- Line: 2310
- Section: 'Nystrom method details'

```rst
.. math::
   :label: e1-decomposition

   E_1(z) = \bigl[-\ln z - \gamma\bigr] + R(z),
   \qquad R(z) \equiv E_1(z) + \ln z + \gamma,\quad R(0) = 0.
```

## `mode-conservation-target`

- File: `docs/theory/peierls_nystrom.rst`
- Line: 3131
- Section: 'The structural obstruction — the :math:`c_{\\rm in}` remapping'

```rst
.. math::
   :label: mode-conservation-target

   W_{\rm oo}[n, n] \;+\; W_{\rm io}[n, n] \;=\; \delta_{n, 0}
   \quad\text{(naive per-mode conservation)}.
```

## `peierls-3d`

- File: `docs/theory/peierls.rst`
- Line: 172
- Section: 'The transport equation in integral form'

```rst
.. math::
   :label: peierls-3d

   \phi(\mathbf r) = \int Q(\mathbf r')\,
       \frac{e^{-\Sigt{}\,|\mathbf r-\mathbf r'|}}
            {4\pi |\mathbf r-\mathbf r'|^2}\,\mathrm d^3\mathbf r',
```

## `peierls-M-rank-1`

- File: `docs/theory/peierls_nystrom.rst`
- Line: 3290
- Section: 'Why F.4 works at rank-1 but does not generalise (Issue #122 close-out)'

```rst
.. math::
   :label: peierls-M-rank-1

   M^{(1)} \;=\; \tfrac{\sqrt{2}}{2}\;\approx\;0.7071,
   \qquad
   (B^{\mu})^{(1)} \;=\; \tfrac{1}{2}.
```

## `peierls-M-rank-2`

- File: `docs/theory/peierls_nystrom.rst`
- Line: 3297
- Section: 'Why F.4 works at rank-1 but does not generalise (Issue #122 close-out)'

```rst
.. math::
   :label: peierls-M-rank-2

   M^{(2)} \;=\;
   \begin{pmatrix}
     \tfrac{\sqrt{2}}{2} & \tfrac{\sqrt{6}}{6} \\[4pt]
     0                   & \tfrac{\sqrt{3}}{3}
   \end{pmatrix},
   \qquad
   (B^{\mu})^{(2)} \;=\;
   \begin{pmatrix}
     \tfrac{1}{2} & \tfrac{\sqrt{3}}{6} \\[4pt]
     \tfrac{\sqrt{3}}{6} & \tfrac{1}{2}
   \end{pmatrix}.
```

## `peierls-WM-WL-asymmetric`

- File: `docs/theory/peierls_nystrom.rst`
- Line: 3324
- Section: 'Why F.4 works at rank-1 but does not generalise (Issue #122 close-out)'

```rst
.. math::
   :label: peierls-WM-WL-asymmetric

   W_M \;=\; B^{\mu}\,W_L,
   \qquad
   B^{\mu} \;=\; M^{\!\top} M,
```

## `peierls-bc-general`

- File: `docs/theory/peierls.rst`
- Line: 296
- Section: 'Boundary conditions parametrised by α + β (Sanchez 1986 Eq. A3.a)'

```rst
.. math::
   :label: peierls-bc-general

   \psi(\mathbf r_b,\Omega) = K(\Omega) +
       \alpha\,\psi(\mathbf r_b,\Omega_R) +
       \beta\,\chi(\Omega)\!\int_{\Omega'\cdot n>0}
            \psi(\mathbf r_b,\Omega')\,(\Omega'\cdot n)\,\mathrm d\Omega',
       \qquad \Omega\cdot n \le 0,
```

## `peierls-bc-operator`

- File: `docs/theory/peierls_nystrom.rst`
- Line: 8342
- Section: 'Section 24 — From continuous integral equation to finite tensor'

```rst
.. math::
   :label: peierls-bc-operator

   S_{\rm bc}(r) \;=\; (T_{\rm bc}\,q)(r)
     \;=\; \int_V K_{\rm bc}(r, r')\,q(r')\,\mathrm d V'.
```

## `peierls-boltzmann`

- File: `docs/theory/peierls.rst`
- Line: 139
- Section: 'The transport equation in integral form'

```rst
.. math::
   :label: peierls-boltzmann

   \Omega\cdot\nabla\psi(\mathbf r,\Omega) + \Sigt{}\,
       \psi(\mathbf r,\Omega) = \frac{1}{4\pi}\,Q(\mathbf r),
```

## `peierls-change-of-basis`

- File: `docs/theory/peierls_nystrom.rst`
- Line: 3279
- Section: 'Why F.4 works at rank-1 but does not generalise (Issue #122 close-out)'

```rst
.. math::
   :label: peierls-change-of-basis

   M_{nm} \;=\; \langle \psi_n^M, \phi_m^L \rangle_M,
   \qquad
   B^{\mu}_{mn} \;=\; \langle \phi_m^L, \phi_n^L \rangle_M
                \;=\; (M^{\!\top} M)_{mn}.
```

## `peierls-class-b-Jn-canonical`

- File: `docs/theory/peierls_nystrom.rst`
- Line: 3742
- Section: 'Root cause — Probe G normalisation mismatch'

```rst
.. math::
   :label: peierls-class-b-Jn-canonical

   J^{+}_n(r_i) \;=\; \frac{1}{A_d}\,
                      \int_\Omega \tilde P_n(\mu_{\rm exit})\,
                      \rho_{\max}^{\,2}(r_i,\Omega)\,e^{-\tau}\,
                      \mathrm d\Omega,
```

## `peierls-class-b-Pss-homogeneous`

- File: `docs/theory/peierls_nystrom.rst`
- Line: 3940
- Section: 'Surface-to-surface probability'

```rst
.. math::
   :label: peierls-class-b-Pss-homogeneous

   P_{ss}(\Sigma_t, R) \;=\;
       2\!\int_0^{\pi/2}\cos\theta'\,\sin\theta'\,
                       e^{-2\Sigma_t R\cos\theta'}\,d\theta'
       \;=\;
       \frac{1 - (1 + 2\tau_R)\,e^{-2\tau_R}}{2\,\tau_R^{\,2}}
```

## `peierls-class-b-hebert-closure`

- File: `docs/theory/peierls_nystrom.rst`
- Line: 3917
- Section: 'Class B sphere — Hébert (1−P_ss)⁻¹ resolution (Issue #132 partial fix)'

```rst
.. math::
   :label: peierls-class-b-hebert-closure

   \mathbb{P}_{\rm white} \;=\;
       \mathbb{P}_{\rm vac}
       \;+\;
       \frac{\beta^{+}}{1 - \beta^{+}\,P_{ss}}\;
       \mathbf{P}_{iS}\,\mathbf{P}_{Sj}^{\top}
```

## `peierls-cyl-Pss-derivation`

- File: `docs/theory/peierls_nystrom.rst`
- Line: 4134
- Section: 'Follow-up directions'

```rst
.. math::
   :label: peierls-cyl-Pss-derivation

   P_{ss}^{\rm cyl}(\Sigma_t, R) = \frac{4}{\pi}\!\int_0^{\pi/2}
       \cos\alpha\;\mathrm{Ki}_3\!\bigl(2\Sigma_t R\cos\alpha\bigr)\,d\alpha
```

## `peierls-cyl-foundations`

- File: `docs/theory/peierls.rst`
- Line: 217
- Section: 'Cylinder reduction (Bickley-Naylor :math:`\\mathrm{Ki}_n` machinery)'

```rst
.. math::
   :label: peierls-cyl-foundations

   \Sigt{}\,\phi(r) = \int_0^R Q(r')\,\frac{r'}{|r-r'|}\,
       \mathrm{Ki}_1(\Sigt{}\,|r-r'|)\,\mathrm d r'\,\cdots
       \;+\; \text{angular details},
```

## `peierls-factored-kernel`

- File: `docs/theory/peierls_nystrom.rst`
- Line: 8471
- Section: 'Section 25 — The factored form :math:`K_{\\rm bc} = G\\,R\\,P`'

```rst
.. math::
   :label: peierls-factored-kernel

   \boxed{\;
     \mathbf K_{\rm bc} \;=\; G \cdot R \cdot P
   \;}
   \qquad
   \Longleftrightarrow
   \qquad
   (\mathbf K_{\rm bc})^{i}{}_{j}
     \;=\; G^{i}{}_{n}\,R^{n}{}_{m}\,P^{m}{}_{j},
```

## `peierls-greens-A1-split`

- File: `docs/theory/trajectory_resolvent.rst`
- Line: 436
- Section: "The angle-resolved Green's function with BC absorbed"

```rst
.. math::
   :label: peierls-greens-A1-split

   t(r',\Omega' \to r,\Omega) = \bar t(r',\Omega' \to r,\Omega) +
       t_h(r',\Omega' \to r,\Omega).
```

## `peierls-greens-A5-specular`

- File: `docs/theory/trajectory_resolvent.rst`
- Line: 454
- Section: "The angle-resolved Green's function with BC absorbed"

```rst
.. math::
   :label: peierls-greens-A5-specular

   t_h(r' \to r,\mu) = \frac{1}{2\pi A}\,e^{-\tau_+ - \tau_-}\,
       T(\mu_+) \cdot \frac{\delta(\mu_- - \mu_+)}{\mu_+}
```

## `peierls-greens-L0`

- File: `docs/theory/trajectory_resolvent.rst`
- Line: 510
- Section: 'Trajectory geometry'

```rst
.. math::
   :label: peierls-greens-L0

   L_0(r_i, \mu_q) = \sqrt{R^2 - r_i^2 (1 - \mu_q^2)} - r_i \mu_q,
```

## `peierls-greens-Lp`

- File: `docs/theory/trajectory_resolvent.rst`
- Line: 528
- Section: 'Trajectory geometry'

```rst
.. math::
   :label: peierls-greens-Lp

   L_p(\mu_{\rm surf}) = 2 R \mu_{\rm surf}.
```

## `peierls-greens-T-alpha`

- File: `docs/theory/trajectory_resolvent.rst`
- Line: 1055
- Section: 'Bounce-sum closure with α'

```rst
.. math::
   :label: peierls-greens-T-alpha

   \psi_{\rm surf}(\mu_{\rm surf})
   = \frac{\alpha\,B(\mu_{\rm surf})}
          {1 - \alpha\,e^{-\Sigt{}\,L_p(\mu_{\rm surf})}}.
```

## `peierls-greens-T-mu-surf`

- File: `docs/theory/trajectory_resolvent.rst`
- Line: 606
- Section: 'Geometric bounce sum'

```rst
.. math::
   :label: peierls-greens-T-mu-surf

   \psi_{\rm surf}(\mu_{\rm surf}) = T(\mu_{\rm surf})\,B(\mu_{\rm surf}),
   \qquad
   T(\mu_{\rm surf}) = \frac{1}{1 - e^{-\Sigt{}\,L_p(\mu_{\rm surf})}}.
```

## `peierls-greens-T00-integrand`

- File: `docs/theory/trajectory_resolvent.rst`
- Line: 786
- Section: 'V_α2 — rank-1 equivalence to Hébert white BC'

```rst
.. math::
   :label: peierls-greens-T00-integrand

   T_{00}^{\rm sphere} = 2 \int_0^1 \mu\,\tilde P_0(\mu)^2\,
       e^{-2\Sigt{} R \mu}\,\mathrm d\mu
   = 2 \int_0^1 \mu\,e^{-2\Sigt{} R \mu}\,\mathrm d\mu,
```

## `peierls-greens-V-alpha-1`

- File: `docs/theory/trajectory_resolvent.rst`
- Line: 656
- Section: 'V_α1 — closed-sphere k_inf identity'

```rst
.. math::
   :label: peierls-greens-V-alpha-1

   (K \cdot 1)(r, \mu) = \omega_0,
   \qquad \omega_0 = \frac{\Sigs{}}{\Sigt{}},
```

## `peierls-greens-V-alpha-2`

- File: `docs/theory/trajectory_resolvent.rst`
- Line: 803
- Section: 'V_α2 — rank-1 equivalence to Hébert white BC'

```rst
.. math::
   :label: peierls-greens-V-alpha-2

   T_{00}^{\rm sphere} = P_{ss}^{\rm sphere} =
       \frac{1 - (1 + 2\tau_R)\,e^{-2\tau_R}}{2\,\tau_R^{2}},
   \qquad \tau_R = \Sigt{}\,R.
```

## `peierls-greens-V-alpha-3`

- File: `docs/theory/trajectory_resolvent.rst`
- Line: 835
- Section: 'V_α3 — vacuum reduction at α = 0'

```rst
.. math::
   :label: peierls-greens-V-alpha-3

   g_h(\rho' \to \rho)\bigr|_{\alpha = 0} = 0,
```

## `peierls-greens-bounce-period-integral`

- File: `docs/theory/trajectory_resolvent.rst`
- Line: 569
- Section: 'Bounce-period integral'

```rst
.. math::
   :label: peierls-greens-bounce-period-integral

   B(\mu_{\rm surf}) = \int_0^{L_p(\mu_{\rm surf})}
       q\bigl(|r_{\rm chord}(s)|\bigr)\,
       e^{-\Sigt{}\,s}\,\mathrm d s,
```

## `peierls-greens-bounce-sum-alpha`

- File: `docs/theory/trajectory_resolvent.rst`
- Line: 1047
- Section: 'Bounce-sum closure with α'

```rst
.. math::
   :label: peierls-greens-bounce-sum-alpha

   \psi_{\rm surf}(\mu_{\rm surf})
   = B(\mu_{\rm surf}) + \alpha\,e^{-\Sigt{}\,L_p}\,\psi_{\rm surf},
```

## `peierls-greens-defining-bvp`

- File: `docs/theory/trajectory_resolvent.rst`
- Line: 415
- Section: "The angle-resolved Green's function with BC absorbed"

```rst
.. math::
   :label: peierls-greens-defining-bvp

   \begin{aligned}
     (\Omega \cdot \nabla + \Sigt{})\,t &= \delta(r - r')\,
        \delta_2(\Omega \cdot \Omega'),
        \quad a_- \le \rho < a, \\
     t(r',\Omega' \to r,\Omega) &= \alpha\,t(r',\Omega' \to r,\Omega_R),
        \quad |\rho| = a,\;\Omega \cdot n \le 0,
   \end{aligned}
```

## `peierls-greens-fixed-source-iteration`

- File: `docs/theory/trajectory_resolvent.rst`
- Line: 1548
- Section: 'Fixed-source iteration'

```rst
.. math::
   :label: peierls-greens-fixed-source-iteration

   \psi_g^{(n+1)}(r,\mu) \;=\; K_g\!\left[\,
       \tfrac{1}{4\pi}\!\left(\,\sum_{g'}\Sigs{g'\to g,\,k(r)}\,
           \phi_{g'}^{(n)}(r) + Q_{{\rm ext},k(r),g}\right)\,
       \right] (r,\mu),
```

## `peierls-greens-function-architecture`

- File: `docs/theory/trajectory_resolvent.rst`
- Line: 628
- Section: 'Total angular flux'

```rst
.. math::
   :label: peierls-greens-function-architecture

   \boxed{\;
   \psi(r_i, \mu_q) = F(r_i, \mu_q) +
       e^{-\Sigt{}\,L_0(r_i,\mu_q)}\,
       T(\mu_{\rm surf})\,B(\mu_{\rm surf})
   \;}
```

## `peierls-greens-garcia-convention`

- File: `docs/theory/trajectory_resolvent.rst`
- Line: 1528
- Section: 'Convention conversion to Garcia table 5'

```rst
.. math::
   :label: peierls-greens-garcia-convention

   \phi_{\rm mine}(r) \;=\; \frac{2\pi}{4\pi}\,\phi_{\rm Garcia}(r)
                       \;=\; \tfrac{1}{2}\,\phi_{\rm Garcia}(r).
```

## `peierls-greens-hollow-sph-outer-only-resolvent`

- File: `docs/theory/trajectory_resolvent.rst`
- Line: 3467
- Section: 'Outer-only resolvent (rank-1)'

```rst
.. math::
   :label: peierls-greens-hollow-sph-outer-only-resolvent

   \psi_{\rm surf}(b) = \frac{\alpha_{\rm out}\,B(b; q)}
                              {1 - \alpha_{\rm out}\,
                               e^{-\Sigma_t\,L_{\rm period}(b)}},
```

## `peierls-greens-k-inf`

- File: `docs/theory/trajectory_resolvent.rst`
- Line: 670
- Section: 'V_α1 — closed-sphere k_inf identity'

```rst
.. math::
   :label: peierls-greens-k-inf

   k_{\rm eff} = \kinf = \frac{\nSigf{}}{\Siga{}},
   \qquad \Siga{} = \Sigt{} - \Sigs{}.
```

## `peierls-greens-mg-source`

- File: `docs/theory/trajectory_resolvent.rst`
- Line: 1204
- Section: 'Multi-group operator action'

```rst
.. math::
   :label: peierls-greens-mg-source

   q_g(r) \;=\; \sum_{g'=1}^{G} \Sigs{g'\to g}\,\phi_{g'}(r) \;+\;
       \frac{\chi_g}{k_{\rm eff}}\,\sum_{g'=1}^{G}\nSigf{g'}\,\phi_{g'}(r),
```

## `peierls-greens-mr-piecewise-tau`

- File: `docs/theory/trajectory_resolvent.rst`
- Line: 1357
- Section: 'Piecewise optical depth'

```rst
.. math::
   :label: peierls-greens-mr-piecewise-tau

   \tau(s) \;=\; \sum_{(s_a, s_b, k)\,\subset\,[0,s]}
       \Sigt{,k}\,(\min(s, s_b) - s_a),
```

## `peierls-greens-mr-trajectory-segments`

- File: `docs/theory/trajectory_resolvent.rst`
- Line: 1324
- Section: 'Trajectory and bounce-period segment decomposition'

```rst
.. math::
   :label: peierls-greens-mr-trajectory-segments

   r(s)^2 \;=\; r_i^2 - 2 r_i \mu_q s + s^2,
   \qquad s \in [0, L_{\rm back}],
```

## `peierls-greens-mu-surf`

- File: `docs/theory/trajectory_resolvent.rst`
- Line: 517
- Section: 'Trajectory geometry'

```rst
.. math::
   :label: peierls-greens-mu-surf

   \mu_{\rm surf}(r_i, \mu_q) = \frac{1}{R}\sqrt{R^2 - r_i^2 (1 - \mu_q^2)}.
```

## `peierls-greens-sanchez-A6`

- File: `docs/theory/trajectory_resolvent.rst`
- Line: 362
- Section: 'Motivation: the Phase 5 retreat made angle-resolution structurally necessary'

```rst
.. math::
   :label: peierls-greens-sanchez-A6

   g_h(\rho'\to\rho) = 2\alpha \int_{\mu_0}^{1} T(\mu_-)\,\mu_*^{-1}\,
       \cosh(\rho\mu)\,\cosh(\rho'\mu_*)\,e^{-2 a \mu_-}\,\mathrm d\mu
```

## `peierls-greens-slab-asym-monodromy`

- File: `docs/theory/trajectory_resolvent.rst`
- Line: 2947
- Section: 'Rank-2 monodromy and resolvent'

```rst
.. math::
   :label: peierls-greens-slab-asym-monodromy

   S(\alpha_L, \alpha_R, \tau) = \begin{pmatrix}
       0                              & \alpha_L\,e^{-\tau} \\
       \alpha_R\,e^{-\tau}            & 0
   \end{pmatrix}, \qquad
   \tau = \Sigma_t\,L /|\mu|
```

## `peierls-greens-slab-bounce-period`

- File: `docs/theory/trajectory_resolvent.rst`
- Line: 2622
- Section: 'Slab phase-space and 2-bounce-per-period structure'

```rst
.. math::
   :label: peierls-greens-slab-bounce-period

   L_{\rm period}^{\rm slab}(\mu) = \frac{2 L}{|\mu|}.
```

## `peierls-greens-surface-fixed-point`

- File: `docs/theory/trajectory_resolvent.rst`
- Line: 595
- Section: 'Geometric bounce sum'

```rst
.. math::
   :label: peierls-greens-surface-fixed-point

   \psi_{\rm surf} = B(\mu_{\rm surf}) +
       e^{-\Sigt{}\,L_p}\,\psi_{\rm surf}.
```

## `peierls-greens-trajectory-integral`

- File: `docs/theory/trajectory_resolvent.rst`
- Line: 544
- Section: 'First-leg trajectory integral'

```rst
.. math::
   :label: peierls-greens-trajectory-integral

   F(r_i, \mu_q) = \int_0^{L_0(r_i,\mu_q)}
       q\bigl(|r_i - s\,\Omega_\mu|\bigr)\,
       e^{-\Sigt{}\,s}\,\mathrm d s.
```

## `peierls-greens-unification-resolvent`

- File: `docs/theory/trajectory_resolvent.rst`
- Line: 2299
- Section: 'The cross-domain frame'

```rst
.. math::
   :label: peierls-greens-unification-resolvent

   T \;=\; (I - S)^{-1}, \qquad
   S = \alpha \cdot R_{\rm chord},
```

## `peierls-half-range-inner-products`

- File: `docs/theory/peierls_nystrom.rst`
- Line: 3266
- Section: 'Why F.4 works at rank-1 but does not generalise (Issue #122 close-out)'

```rst
.. math::
   :label: peierls-half-range-inner-products

   \langle f, g \rangle_L \;=\; \int_0^1 f(\mu)\,g(\mu)\,\mathrm d\mu,
   \qquad
   \langle f, g \rangle_M \;=\; \int_0^1 f(\mu)\,g(\mu)\,\mu\,\mathrm d\mu.
```

## `peierls-integral-form`

- File: `docs/theory/peierls.rst`
- Line: 162
- Section: 'The transport equation in integral form'

```rst
.. math::
   :label: peierls-integral-form

   \psi(\mathbf r,\Omega) = \int_0^\infty
       \frac{Q(\mathbf r - s\Omega)}{4\pi}\,
       e^{-\Sigt{}\,s}\,\mathrm d s.
```

## `peierls-mg-operator`

- File: `docs/theory/peierls_nystrom.rst`
- Line: 1088
- Section: 'Subsection — The multi-group operator form'

```rst
.. math::
   :label: peierls-mg-operator

   \Sigma_{t,g_{\rm out}}(r_i)\,\varphi_{g_{\rm out}}(r_i)
     \;=\;\sum_{j}\!K^{(g_{\rm out})}_{ij}\!\!
         \sum_{g_{\rm in}}\!\Biggl[\,
           \Sigma_{s,\,g_{\rm in} \to g_{\rm out}}(r_j)\,
           \varphi_{g_{\rm in}}(r_j)
           + \frac{1}{k}\,\chi_{g_{\rm out}}(r_i)\,
             \nu\Sigma_{f,g_{\rm in}}(r_j)\,
             \varphi_{g_{\rm in}}(r_j)
         \,\Biggr],
```

## `peierls-operator-factorisation`

- File: `docs/theory/peierls_nystrom.rst`
- Line: 8446
- Section: 'Section 25 — The factored form :math:`K_{\\rm bc} = G\\,R\\,P`'

```rst
.. math::
   :label: peierls-operator-factorisation

   T_{\rm bc} \;=\; G_\infty \;\circ\; R_\infty \;\circ\; P_\infty.
```

## `peierls-operator-form`

- File: `docs/theory/peierls_nystrom.rst`
- Line: 8314
- Section: 'Section 24 — From continuous integral equation to finite tensor'

```rst
.. math::
   :label: peierls-operator-form

   \Sigma_t\,\varphi \;=\; T_{\rm vol}\,q \;+\; S_{\rm bc},
```

## `peierls-rank-n-P-esc-moment`

- File: `docs/theory/peierls_nystrom.rst`
- Line: 2961
- Section: 'Surface-to-observer Jacobian :math:`(\\rho_{\\max}/R)^2`'

```rst
.. math::
   :label: peierls-rank-n-P-esc-moment

   P_{\rm esc}^{(n)}(r_i) \;=\; C_d\!\int_0^\pi W_\Omega(\Omega)\,
                    \Bigl(\tfrac{\rho_{\max}(r_i,\Omega)}{R}\Bigr)^{\!2}\,
                    \tilde P_n\!\bigl(\mu_{\rm exit}\bigr)\,
                    K_{\rm esc}(\tau)\,\mathrm d\Omega.
```

## `peierls-rank-n-jacobian-derivation`

- File: `docs/theory/peierls_nystrom.rst`
- Line: 2947
- Section: 'Surface-to-observer Jacobian :math:`(\\rho_{\\max}/R)^2`'

```rst
.. math::
   :label: peierls-rank-n-jacobian-derivation

   J^{+}_n(r_i) \;=\; \frac{1}{A_d}\,
                      \int_\Omega \tilde P_n(\mu_{\rm exit})\,
                      \rho_{\max}^{\,2}(r_i,\Omega)\,e^{-\tau}\,
                      \mathrm d\Omega.
```

## `peierls-slab-foundations`

- File: `docs/theory/peierls.rst`
- Line: 196
- Section: 'Slab reduction (:math:`E_n` machinery)'

```rst
.. math::
   :label: peierls-slab-foundations

   \Sigt{}\,\phi(x) = \frac{1}{2}\int_0^L
       Q(x')\,E_1(\Sigt{}\,|x-x'|)\,\mathrm d x',
```

## `peierls-slab-polar`

- File: `docs/theory/peierls_nystrom.rst`
- Line: 516
- Section: 'Subsection — The slab polar-form equation (recap)'

```rst
.. math::
   :label: peierls-slab-polar

   \Sigma_t(x_i)\,\varphi(x_i)
     \;=\;\frac{1}{2}
     \int_{-1}^{1}\!\mathrm d\mu
     \int_{0}^{\rho_{\max}(x_i,\mu)}
       e^{-\int_0^\rho \Sigma_t(x_i + s\mu)\,\mathrm ds}\,
       q\bigl(x_i + \rho\mu\bigr)\,\mathrm d\rho,
```

## `peierls-specular-M-tridiagonal`

- File: `docs/theory/peierls_nystrom.rst`
- Line: 4385
- Section: 'The R_specular operator — closed form in the half-range Legendre basis'

```rst
.. math::
   :label: peierls-specular-M-tridiagonal

   M_{nn} \;=\; \frac{1}{2(2n+1)},
   \qquad
   M_{n,n+1} \;=\; M_{n+1,n} \;=\;
   \frac{n+1}{2(2n+1)(2n+3)}.
```

## `peierls-specular-R-formula`

- File: `docs/theory/peierls_nystrom.rst`
- Line: 4374
- Section: 'The R_specular operator — closed form in the half-range Legendre basis'

```rst
.. math::
   :label: peierls-specular-R-formula

   R_{\rm spec} \;=\; \tfrac{1}{2}\,M^{-1},
   \qquad
   M_{nm} \;=\; \int_0^1 \mu\,\tilde P_n(\mu)\,\tilde P_m(\mu)\,
   \mathrm d\mu.
```

## `peierls-specular-T-matrix`

- File: `docs/theory/peierls_nystrom.rst`
- Line: 4965
- Section: 'Per-geometry derivation of :math:`T`'

```rst
.. math::
   :label: peierls-specular-T-matrix

   T_{mn}^{\rm sph} \;=\; 2\!\int_0^1\!\mu\,\tilde P_m(\mu)\,\tilde P_n(\mu)\,
                  e^{-\tau(\mu)}\,\mathrm d\mu,
```

## `peierls-specular-multibounce-formula`

- File: `docs/theory/peierls_nystrom.rst`
- Line: 4939
- Section: 'Multi-bounce-corrected specular (sphere / cylinder / slab — Phase 4)'

```rst
.. math::
   :label: peierls-specular-multibounce-formula

   K_{\rm bc}^{\rm spec,mb} \;=\; G \cdot R \cdot
                                  (I - T\,R)^{-1} \cdot P,
```

## `peierls-sph-ps1982-foundations`

- File: `docs/theory/peierls.rst`
- Line: 246
- Section: 'Sphere reduction (radial-integration / cosh-extension)'

```rst
  .. math::
     :label: peierls-sph-ps1982-foundations

     r\phi(r) = \int_0^R x\,Q(x)\,
         \bigl[\,E_1(\Sigt{}\,|r-x|) - E_1(\Sigt{}\,(r+x))\,\bigr]\,
         \mathrm d x,
```

## `peierls-sphere-G-bc`

- File: `docs/theory/collision_probability.rst`
- Line: 3794
- Section: "Surface-to-volume Green's function :math:`G_{\\rm bc}`"

```rst
.. math::
   :label: peierls-sphere-G-bc

   G_{\rm bc}^{\rm sph}(r_i)
     \;=\; 2\!\int_{0}^{\pi}\!\sin\theta\,
       e^{-\tau_{\rm surf}(r_i, \theta)}\,\mathrm d\theta.
```

## `peierls-sphere-equation`

- File: `docs/theory/collision_probability.rst`
- Line: 3446
- Section: 'Observer-centred polar form and Jacobian cancellation'

```rst
.. math::
   :label: peierls-sphere-equation

   \Sigma_t(r_i)\,\varphi(r_i)
     \;=\; \frac{\Sigma_t(r_i)}{2}\!
       \int_{0}^{\pi}\!\sin\theta\,\mathrm d\theta\!
       \int_{0}^{\rho_{\max}(r_i, \theta)}\!\!
         e^{-\tau(r_i, \rho, \theta)}\,
         q\!\bigl(r'(r_i, \rho, \theta)\bigr)\,\mathrm d\rho
     \;+\; S_{\rm bc}(r_i).
```

## `peierls-sphere-nystrom`

- File: `docs/theory/collision_probability.rst`
- Line: 3616
- Section: 'Nyström assembly in polar coordinates'

```rst
.. math::
   :label: peierls-sphere-nystrom

   \Sigma_t(r_i)\,\varphi_i
     \;=\; \sum_{j=1}^{N} K_{ij}\,q_j + S_{\rm bc}(r_i),
```

## `peierls-sphere-ray-optical-depth`

- File: `docs/theory/collision_probability.rst`
- Line: 3648
- Section: 'Ray optical-depth walker'

```rst
.. math::
   :label: peierls-sphere-ray-optical-depth

   \tau(r_i, \rho, \theta)
     \;=\; \int_{0}^{\rho}\!
       \Sigma_t\!\bigl(r'(r_i, s, \theta)\bigr)\,\mathrm ds,
```

## `peierls-svd`

- File: `docs/theory/peierls_nystrom.rst`
- Line: 8368
- Section: 'Section 24 — From continuous integral equation to finite tensor'

```rst
  .. math::
     :label: peierls-svd

     K(r, r') \;=\; \sum_{k=1}^{\infty} \sigma_k\,u_k(r)\,v_k(r'),
     \qquad \sigma_k \to 0,
```

## `peierls-tensor-G-definition`

- File: `docs/theory/peierls_nystrom.rst`
- Line: 8596
- Section: 'Section 27 — The escape and response tensors: where geometry lives'

```rst
.. math::
   :label: peierls-tensor-G-definition

   G^{i}{}_{n} \;\propto\;
     \frac{\Sigma_t(r_i)\,\mathcal G^{(n)}(r_i)}{A_d^{\rm divisor}},
```

## `peierls-tensor-P-definition`

- File: `docs/theory/peierls_nystrom.rst`
- Line: 8590
- Section: 'Section 27 — The escape and response tensors: where geometry lives'

```rst
.. math::
   :label: peierls-tensor-P-definition

   P^{n}{}_{j} \;\propto\;
     r_j^{d-1}\,w_j\,\mathcal P^{(n)}(r_j),
```

## `peierls-vacuum-bc-cylinder`

- File: `docs/theory/peierls_nystrom.rst`
- Line: 2474
- Section: 'Cylinder (infinite axial extent)'

```rst
.. math::
   :label: peierls-vacuum-bc-cylinder

   \varphi_{\rm cyl}(r) \;=\; \frac{1}{\pi\,\Sigma_t}\!\int_0^\pi\!
     \Bigl[\,1 - \mathrm{Ki}_2\!\bigl(\Sigma_t\,L_{2D}(r,\theta')\bigr)\,\Bigr]
     \,\mathrm d\theta'.
```

## `peierls-vacuum-bc-flux`

- File: `docs/theory/peierls_nystrom.rst`
- Line: 2417
- Section: 'Vacuum-BC analytical flux references (2026-04-20 milestone)'

```rst
.. math::
   :label: peierls-vacuum-bc-flux

   \varphi_d(r) \;=\; \int_{\mathcal V}\!K_d(r, r')\,q(r')\,\mathrm dV'
```

## `peierls-vacuum-bc-row-sum-gate`

- File: `docs/theory/peierls_nystrom.rst`
- Line: 2431
- Section: 'Vacuum-BC analytical flux references (2026-04-20 milestone)'

```rst
.. math::
   :label: peierls-vacuum-bc-row-sum-gate

   \sum_{j=1}^{N} K_{ij}\cdot 1 \;\stackrel{!}{=}\; \Sigma_t(r_i)\,\varphi_d(r_i).
```

## `peierls-vacuum-bc-slab`

- File: `docs/theory/peierls_nystrom.rst`
- Line: 2444
- Section: 'Slab'

```rst
.. math::
   :label: peierls-vacuum-bc-slab

   \varphi_{\rm slab}(x) \;=\; \frac{1}{2\,\Sigma_t}
     \Bigl[\,2 - E_2(\Sigma_t\,x) - E_2(\Sigma_t\,(L - x))\,\Bigr].
```

## `peierls-vacuum-bc-sphere`

- File: `docs/theory/peierls_nystrom.rst`
- Line: 2514
- Section: 'Sphere'

```rst
.. math::
   :label: peierls-vacuum-bc-sphere

   \varphi_{\rm sph}(r) \;=\; \frac{1}{2\,\Sigma_t}\!\left[\,
       2 - \int_{-1}^{1}\!\exp\!\Bigl(
           -\Sigma_t\bigl[-r\mu + \sqrt{R^{2} - r^{2} + r^{2}\mu^{2}}\bigr]
         \Bigr)\,\mathrm d\mu\,\right].
```

## `peierls-white-bc-slab`

- File: `docs/theory/peierls_nystrom.rst`
- Line: 2568
- Section: 'White-BC analytical flux — slab (Wigner-Seitz identity)'

```rst
.. math::
   :label: peierls-white-bc-slab

   \varphi_{\rm white}(x) \;\equiv\; \frac{1}{\Sigma_t}
   \qquad (\text{for all } x, \text{ any } L).
```

## `vacuum-bc`

- File: `docs/theory/discrete_ordinates.rst`
- Line: 1189
- Section: 'Supported Types'

```rst
.. math::
   :label: vacuum-bc

   \psi_n^{\rm in} = 0
```
