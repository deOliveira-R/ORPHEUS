# Q5.6.4 second attempt — raw probe outputs (2026-08-11)

Probes live in `$CLAUDE_JOB_DIR/tmp/q64_probe{A,B,C,D,E}_*.py` (EPHEMERAL).
Regenerate from the descriptions in the memo if the job dir is gone.


```text
==============================================================================
A1/A2 — is alpha a DISPLACED exact dome (an edge) or an INFLATED one?
==============================================================================

  exact continuous dome:  alpha_exact(omega) = -w_gl * xi(omega)
    [from the conservative angular term -(1/r) d(xi*psi)/domega]
  production (recursion): alpha_{m+1/2} = alpha_{m-1/2} - w_m eta_m

  n_phi=  8  M= 4  d_omega=0.785398  kappa=1.026172153  (kappa-1 = 2.617e-02)
      max|alpha_prod - kappa*alpha_exact_on_arc| = 7.259e-01   <-- A1
      max|alpha_prod -       alpha_exact_on_arc| = 7.166e-01
      A2: edges with NO real recursion-defined position: 1 of 5   (min RHS = -1.371e-02)

  n_phi= 16  M= 8  d_omega=0.392699  kappa=1.006454543  (kappa-1 = 6.455e-03)
      max|alpha_prod - kappa*alpha_exact_on_arc| = 7.119e-01   <-- A1
      max|alpha_prod -       alpha_exact_on_arc| = 7.096e-01
      A2: edges with NO real recursion-defined position: 1 of 9   (min RHS = -3.347e-03)

  n_phi= 32  M=16  d_omega=0.196350  kappa=1.001608189  (kappa-1 = 1.608e-03)
      max|alpha_prod - kappa*alpha_exact_on_arc| = 7.085e-01   <-- A1
      max|alpha_prod -       alpha_exact_on_arc| = 7.079e-01
      A2: edges with NO real recursion-defined position: 1 of 17   (min RHS = -8.319e-04)

  n_phi= 64  M=32  d_omega=0.098175  kappa=1.000401708  (kappa-1 = 4.017e-04)
      max|alpha_prod - kappa*alpha_exact_on_arc| = 7.076e-01   <-- A1
      max|alpha_prod -       alpha_exact_on_arc| = 7.075e-01
      A2: edges with NO real recursion-defined position: 1 of 33   (min RHS = -2.077e-04)

  READ: A1 == 0 to machine precision means alpha IS the exact dome
  scaled by kappa.  A scale factor is an AMPLITUDE, not a position:
  there is no edge at which the exact dome takes the value alpha has,
  wherever kappa*xi > sin_theta -- which A2 counts.

==============================================================================
A3 — nu-closure, run in ETA (the memo's table) and in OMEGA
==============================================================================

  nu_{1/2} = start ;  nu_{m+1/2} = (x_m - (1-tau_m) nu_{m-1/2}) / tau_m
  A partition-consistent tau lands nu exactly on the far endpoint.

  n_phi convention                   nu_close/end (ETA)  nu_close/end (OMEGA)
  ---------------------------------------------------------------------------
      8 chord (retired prod.)                  1.000000              0.159076
      8 chord + [1/2,1] absorber               1.016389              0.079538
      8 arc / omega-mid (LANDED)               1.000000              0.148847
      8 tau == 1/2 (diamond)                   1.164784              0.000000

     16 chord (retired prod.)                  1.000000              0.081139
     16 chord + [1/2,1] absorber               1.001930              0.040570
     16 arc / omega-mid (LANDED)               1.000000              0.074765
     16 tau == 1/2 (diamond)                   1.039182             -0.000000

     32 chord (retired prod.)                  1.000000              0.040772
     32 chord + [1/2,1] absorber               1.000238              0.020386
     32 arc / omega-mid (LANDED)               1.000000              0.037427
     32 tau == 1/2 (diamond)                   1.009677              0.000000

     64 chord (retired prod.)                  1.000000              0.020411
     64 chord + [1/2,1] absorber               1.000030              0.010206
     64 arc / omega-mid (LANDED)               1.000000              0.018719
     64 tau == 1/2 (diamond)                   1.002412              0.000000

  READ: the ETA column is the memo's table (1.000000 == exact).
  The OMEGA column is the same test in the MARCH variable; 0.000000
  == exact.  A convention that closes in one and not the other is
  affine in that variable and not the other -- which is the whole
  content of 'which partition?'.
```

```text
==============================================================================
B1 — alpha's SHAPE: is alpha_prod proportional to xi at the arc edges?
==============================================================================

  n_phi=  8  M= 4  kappa=1.026172153
      alpha/xi(e_arc) over interior edges: rel spread = 1.555e-16   (constant <=> alpha ~ xi)
      fitted constant c = 0.713917910716 = kappa * 0.695709690275
      => alpha is the arc-edge dome INFLATED by kappa (2.617% too tall)
  n_phi= 16  M= 8  kappa=1.006454543
      alpha/xi(e_arc) over interior edges: rel spread = 1.110e-15   (constant <=> alpha ~ xi)
      fitted constant c = 0.700200178247 = kappa * 0.695709690275
      => alpha is the arc-edge dome INFLATED by kappa (0.645% too tall)
  n_phi= 32  M=16  kappa=1.001608189
      alpha/xi(e_arc) over interior edges: rel spread = 1.912e-15   (constant <=> alpha ~ xi)
      fitted constant c = 0.696828523004 = kappa * 0.695709690275
      => alpha is the arc-edge dome INFLATED by kappa (0.161% too tall)
  n_phi= 64  M=32  kappa=1.000401708
      alpha/xi(e_arc) over interior edges: rel spread = 2.313e-14   (constant <=> alpha ~ xi)
      fitted constant c = 0.695989162531 = kappa * 0.695709690275
      => alpha is the arc-edge dome INFLATED by kappa (0.040% too tall)

  READ: a constant ratio means alpha's SHAPE is exactly xi at the arc
  edges and its AMPLITUDE is kappa too large.  An amplitude error has
  no edge to sit at: probe A's A2 counted 1 edge per level (the
  mid-level one, where xi is maximal) at which kappa*xi > sin_theta,
  so NO real direction carries the value alpha has.  Candidate (3) is
  ILL-POSED, at every quadrature order.

==============================================================================
B2 — ERROR AMPLIFICATION of the angular recurrence
==============================================================================

  per-step factor |(1 - tau_m)/tau_m| ;  tau >= 1/2  <=>  factor <= 1

  n_phi  convention                       min tau  max factor  worst running prod  #tau<1/2
  -----------------------------------------------------------------------------------------
      8  (1) chord                       0.219545      3.5549        5.027339e+00         2
      8  (1c) chord + absorber           0.500000      1.0000        1.000000e+00         0
      8  (2) arc / omega-mid  [LANDED]   0.259892      2.8478        3.359161e+00         2
      8  (C) tau == 1/2                  0.500000      1.0000        1.000000e+00         0

     16  (1) chord                       0.204689      3.8855        1.015317e+01         4
     16  (1c) chord + absorber           0.500000      1.0000        1.000000e+00         0
     16  (2) arc / omega-mid  [LANDED]   0.252425      2.9616        4.728870e+00         4
     16  (C) tau == 1/2                  0.500000      1.0000        1.000000e+00         0

     32  (1) chord                       0.201161      3.9712        2.035547e+01         8
     32  (1c) chord + absorber           0.500000      1.0000        1.000000e+00         0
     32  (2) arc / omega-mid  [LANDED]   0.250603      2.9904        6.679689e+00         8
     32  (C) tau == 1/2                  0.500000      1.0000        1.000000e+00         0

     64  (1) chord                       0.200289      3.9928        4.073548e+01        16
     64  (1c) chord + absorber           0.500000      1.0000        1.000000e+00         0
     64  (2) arc / omega-mid  [LANDED]   0.250151      2.9976        9.443672e+00        16
     64  (C) tau == 1/2                  0.500000      1.0000        1.000000e+00         0

  READ: a worst running product >> 1 means the recurrence AMPLIFIES an
  upstream error before it reaches the outward directions.  This is
  what the retired [1/2, 1] absorber enforced (factor <= 1 by
  construction) -- a STABILITY guard on the recurrence, which is a
  different claim from 'it compensates a wrong partition'.

==============================================================================
B3 — the P1 closure defect at the arc-edge faces
==============================================================================

  psi affine in the direction cosines (diffusion limit).  Feed eta and
  xi through the closure; compare the produced face value to the TRUE
  value at the arc edge.  Weighted defect uses xi(e_arc), the
  coefficient alpha multiplies in the balance.

  n_phi  convention                      max|d eta|   max|d xi|  xi-weighted
  --------------------------------------------------------------------------
      8  (1) chord                        2.736e-02   6.277e-01    3.191e-01
      8  (1c) chord + absorber            7.150e-02   1.332e-01    3.011e-02
      8  (2) arc / omega-mid  [LANDED]    1.110e-16   3.891e-01    1.869e-01
      8  (C) tau == 1/2                   8.377e-02   4.189e-02    3.011e-02

     16  (1) chord                        9.025e-03   6.541e-01    3.325e-01
     16  (1c) chord + absorber            1.916e-02   6.540e-02    7.887e-03
     16  (2) arc / omega-mid  [LANDED]    1.670e-16   2.827e-01    1.391e-01
     16  (C) tau == 1/2                   1.992e-02   9.960e-03    7.779e-03

     32  (1) chord                        2.401e-03   6.617e-01    3.364e-01
     32  (1c) chord + absorber            4.872e-03   3.263e-02    2.211e-03
     32  (2) arc / omega-mid  [LANDED]    5.274e-16   2.004e-01    1.003e-01
     32  (C) tau == 1/2                   4.920e-03   2.460e-03    1.921e-03

     64  (1) chord                        6.094e-04   6.637e-01    3.374e-01
     64  (1c) chord + absorber            1.223e-03   1.631e-02    5.790e-04
     64  (2) arc / omega-mid  [LANDED]    6.939e-16   1.415e-01    7.157e-02
     64  (C) tau == 1/2                   1.226e-03   6.131e-04    4.791e-04

  READ: no single tau can be exact for BOTH eta and xi -- requiring it
  gives  sin p + sin q = sin(p + q)  for the two half-cell angles,
  whose only solutions are p = 0 or q = 0 (the node ON an edge).  So
  every convention trades one against the other; the xi-weighted
  column is the trade as the balance equation actually sees it.
```

```text
==============================================================================
tau-BLINDNESS AUDIT — quadrature and partition held FIXED
  folded_product(n_mu=4, n_phi=16), level 0, M = 8
==============================================================================

                                                    production (arc)                tau == 1/2        GARBAGE tau == 0.7   GARBAGE random in (0,1)
  ------------------------------------------------------------------------------------------------------------------------------------------------
  P1   c = sum w eta^2  (Lathrop 29/53-54)            0.282432589924            0.282432589924            0.282432589924            0.282432589924
  BMC  beta  (contamination, Eq. 41/75)                            0                         0                         0                         0
  Lath beta  (alpha defect, Eq. 25)                   0.741555747146            0.741555747146            0.741555747146            0.741555747146
  nu-closure ratio (BMC Eq. 43 forward)                            1             1.03918231642             1.01319215547            0.972815198734
  * amplification  max prod |(1-t)/t|                  4.72887003107                         1            0.428571428571             2.18309886184
  * P1 closure defect at the arc faces                0.139106457701          0.00777890941286           0.0211225798756            0.274253487045

  Instruments marked * are the two this session added.  Every
  UNMARKED row is CONSTANT across all four columns -- including the
  two garbage taus -- so it is tau-blind by construction: it reads
  the quadrature and the partition only.

  Consequence: 'the diffusion limit is why tau exists' is true of the
  literature's DERIVATION of tau, but the published beta / c
  functionals cannot GRADE a tau.  A suite of tau-blind instruments
  will certify any tau, which is what happened at 6.4.
```

```text
==============================================================================
E1 — SPHERE: the march variable IS mu, so nothing moves
==============================================================================

  N=  4   max|P2_in_march - production| = 0.000e+00   (tau range [0.399200, 0.600800])
  N=  8   max|P2_in_march - production| = 0.000e+00   (tau range [0.392282, 0.607718])
  N= 16   max|P2_in_march - production| = 0.000e+00   (tau range [0.390354, 0.609646])
  N= 32   max|P2_in_march - production| = 0.000e+00   (tau range [0.389840, 0.610160])
  N= 64   max|P2_in_march - production| = 0.000e+00   (tau range [0.389708, 0.610292])

  READ: the sphere arm is BMC Eq. 12 verbatim and it stays verbatim.
  The re-pose is a cylinder-only change of march variable.

==============================================================================
E2 — CYLINDER: P2 in omega == 1/2 as a CONSEQUENCE, not a constant
==============================================================================

  folded_product(2,  4)   max|P2_in_omega - 1/2| over all levels = 0.000e+00
  folded_product(2,  8)   max|P2_in_omega - 1/2| over all levels = 0.000e+00
  folded_product(2, 16)   max|P2_in_omega - 1/2| over all levels = 5.551e-16
  folded_product(2, 32)   max|P2_in_omega - 1/2| over all levels = 1.110e-15
  folded_product(2, 64)   max|P2_in_omega - 1/2| over all levels = 4.552e-15

  folded_product(4,  4)   max|P2_in_omega - 1/2| over all levels = 0.000e+00
  folded_product(4,  8)   max|P2_in_omega - 1/2| over all levels = 0.000e+00
  folded_product(4, 16)   max|P2_in_omega - 1/2| over all levels = 5.551e-16
  folded_product(4, 32)   max|P2_in_omega - 1/2| over all levels = 1.110e-15
  folded_product(4, 64)   max|P2_in_omega - 1/2| over all levels = 8.993e-15

  control — a NON-equispaced monotone arc (the body does real work):
    omega     = [2.9  2.2  1.9  1.   0.35]
    tau       = [0.408377, 0.700000, 0.250000, 0.580645, 0.481481]
    max|tau - 1/2| = 0.2500   (non-trivial => not a hardcoded constant)
    all in [0,1]? True   -- P3 still a theorem on a monotone arc

==============================================================================
E3 — the guards and the R12a facts under tau == 1/2
==============================================================================

  _assert_tau_within_unit_interval: 1/2 in [0,1] trivially.
  The tau_0 = 0 endpoint pathology of the CHORD partition (a node on
  Sigma divides by zero) cannot occur: tau is 1/2, never 0.

  folded_product(2,  4)  march-start facts (on_edge_node, degenerate, carrying) per level:
      [(False, False, True), (False, False, True)]
  folded_product(4,  8)  march-start facts (on_edge_node, degenerate, carrying) per level:
      [(False, False, True), (False, False, True), (False, False, True), (False, False, True)]
  folded_product(4, 16)  march-start facts (on_edge_node, degenerate, carrying) per level:
      [(False, False, True), (False, False, True), (False, False, True), (False, False, True)]
  folded_product(4, 32)  march-start facts (on_edge_node, degenerate, carrying) per level:
      [(False, False, True), (False, False, True), (False, False, True), (False, False, True)]

  READ: these are integer facts read off xi == 0 and eta_0 == eta_1
  bit-exactly -- they never consulted tau (Q5.4/T26).  A tau re-pose
  cannot move them, so R12a admission is untouched.
```
