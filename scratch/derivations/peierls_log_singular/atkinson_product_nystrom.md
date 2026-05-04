# Atkinson's product-Nyström method for log-singular Fredholm kernels

**Memo for**: numerics-investigator (will implement in
`orpheus/derivations/continuous/fn_method/<slab|sphere>/flux_reconstruction.py`)
and archivist (will document in Sphinx).

**Sources used (both local in `scratch/literature/`)**:

- Atkinson, K.E. (1997). *The Numerical Solution of Integral Equations
  of the Second Kind*. Cambridge Monographs on Applied and Computational
  Mathematics, Vol. 4. Cambridge University Press.
  Chapter 4 §4.2 "Product integration methods" (pp. 116–135).
- Atkinson, K.E. (1972). "The Numerical Solution of Fredholm Integral
  Equations of the Second Kind with Singular Kernels." *Numerische
  Mathematik* 19, 248–259. (Cited in the book as "the kernel was
  assumed capable of being split as ..." — origin of equation
  (4.2.80) below.)

The numbering used below follows the **1997 book** (`(4.2.x)`) when both
exist; the 1972 paper is cited as "[1972]" with its own equation labels.

---

## 0. Notation bridge: Atkinson ↔ ORPHEUS Bickley-Naylor

Atkinson formulates everything for the canonical 1D Fredholm form

$$
\lambda x(t) - \int_a^b K(t, s)\, x(s)\, ds = y(t), \quad a \le t \le b.
\tag{4.2.60}
$$

The two singular cases he treats explicitly are

$$
H(t,s) = \log|t-s|, \qquad H(t,s) = |t-s|^{\gamma-1}, \quad \gamma > 0.
$$

Our slab Peierls kernel is

$$
K_{\text{slab}}(x, y) = \frac{1}{2}\, \mathrm{Ki}_1\!\bigl(\Sigma\,|x-y|\bigr),
$$

where (Abramowitz & Stegun 5.1.30, with our sign convention)

$$
\mathrm{Ki}_1(\tau) = \int_0^{\pi/2} e^{-\tau/\cos\theta}\, d\theta
= \int_0^\infty \mathrm{sech}\,u \cdot e^{-\tau\,\cosh u}\, du.
$$

The small-argument expansion of $\mathrm{Ki}_1$ is

$$
\mathrm{Ki}_1(\tau) = \frac{\pi}{2} + \tau\bigl(\log\tau + \gamma_E - 1 - \log 2\bigr) + O(\tau^2),
\qquad \tau \to 0^+,
$$

where $\gamma_E$ is the Euler–Mascheroni constant. Hence

$$
\boxed{\;
\mathrm{Ki}_1(\Sigma|x-y|) =
\underbrace{\frac{\pi}{2}}_{\text{regular}}
\;+\; \Sigma|x-y|\,\bigl[\log(\Sigma|x-y|)\bigr]
\;+\; \text{regular}\;}
$$

This kernel has the form Atkinson's machinery is *most powerful for*:
the leading term is **bounded** (so $\mathrm{Ki}_1 \in L^\infty$),
but its derivative at $|x-y|=0$ blows up logarithmically, which is what
kills polynomial Gauss-Legendre. Concretely, the singular component
matches Atkinson's general decomposition (4.2.80) with

$$
K_{\text{slab}}(x,y) = \tfrac12 L_1(x,y) H_1(x,y) + \tfrac12 L_2(x,y) H_2(x,y),
$$

| Atkinson | Concrete form | Role |
|----------|---------------|------|
| $L_1(x,y)$ | $\equiv 1$ | smooth weight |
| $H_1(x,y)$ | $\pi/2$ + smooth Taylor remainder of $\mathrm{Ki}_1$ minus its log piece | bounded, integrable |
| $L_2(x,y)$ | $\Sigma\,|x-y|$ | smooth weight ($C^\infty$ except at $x=y$) |
| $H_2(x,y)$ | $\log\bigl(\Sigma|x-y|\bigr)$ | the bare logarithmic singularity |

The map to ORPHEUS variables: Atkinson's $t,s$ are positions $x,y$ on
the slab; the integration interval $[a,b]$ is $[-X,X]$ (full slab) or
$[0,X]$ (half-slab with reflective BC). Atkinson's eigenvalue parameter
$\lambda$ corresponds to our $1$ in the second-kind form
$\phi - \mathcal K\phi = q$ (so $\lambda = 1$); for the criticality
problem $1/k_{\text{eff}}$ multiplies the fission source.

For the **sphere case** the Peierls kernel after radial integration
behaves as $|r-r'|^{-1}\,$Ki$_1(\Sigma|r-r'|)$ in 1D. Despite the
$1/|r-r'|$ prefactor this **also** falls under Atkinson's framework
with $\gamma = 0$ in the $|t-s|^{\gamma-1}$ family — see §3 below.

---

## 1. Method specification: product trapezoidal & product Simpson rules

### 1.1 The model log-kernel equation

Atkinson introduces product integration on the **prototype**

$$
\lambda x(t) - \int_a^b L(t,s)\,\log|s-t|\,x(s)\, ds = y(t), \quad a \le t \le b.
\tag{4.2.61}
$$

The kernel splits as $K(t,s) = L(t,s)\,\log|s-t|$ with $L$ smooth
(several times continuously differentiable) and $x(t)$ assumed smooth
**at first** (graded mesh in §4.2.5 lifts this assumption).

### 1.2 Product trapezoidal rule (the entry-level method)

Place $n+1$ nodes $t_j = a + jh$, $j = 0, \ldots, n$, $h = (b-a)/n$.
**Interpolate $L(t,s)\,x(s)$ piecewise linearly in $s$** at the nodes:

$$
[L(t,s)\,x(s)]_n =
\frac{1}{h}\Bigl[
(t_j - s)\,L(t, t_{j-1})\,x(t_{j-1}) +
(s - t_{j-1})\,L(t, t_j)\,x(t_j)
\Bigr],
\quad t_{j-1} \le s \le t_j.
\tag{4.2.63}
$$

Then **integrate the singular part exactly** against the linear basis:

$$
\mathcal K_n x(t) = \int_a^b [L(t,s)\,x(s)]_n\,\log|s-t|\, ds
= \sum_{j=0}^n w_j(t)\, L(t, t_j)\, x(t_j).
\tag{4.2.64–4.2.65}
$$

The weights $w_j(t)$ are the **load-bearing piece** of the method.
They are *not* tabulated quadrature weights of any classical rule —
they are computed from the analytic integral of the linear basis
times the singular kernel $\log|s-t|$:

$$
\begin{aligned}
w_0(t) &= \frac{1}{h}\int_{t_0}^{t_1} (t_1 - s)\,\log|t-s|\, ds, \\
w_n(t) &= \frac{1}{h}\int_{t_{n-1}}^{t_n} (s - t_{n-1})\,\log|t-s|\, ds, \\
w_j(t) &= \frac{1}{h}\int_{t_{j-1}}^{t_j}(s - t_{j-1})\log|t-s|\, ds
       + \frac{1}{h}\int_{t_j}^{t_{j+1}}(t_{j+1} - s)\log|t-s|\, ds,
\end{aligned}
\tag{4.2.66, 4.2.67}
$$

for the interior nodes $j = 1, \ldots, n-1$.

The discretized integral equation is

$$
\lambda x_n(t) - \sum_{j=0}^n w_j(t)\, L(t, t_j)\, x_n(t_j) = y(t),
\quad a \le t \le b.
\tag{4.2.68}
$$

The Nyström-style implementation: solve at the collocation nodes
first

$$
\lambda x_n(t_i) - \sum_{j=0}^n w_j(t_i)\, L(t_i, t_j)\, x_n(t_j)
= y(t_i),
\quad i = 0, \ldots, n,
\tag{4.2.69}
$$

then reconstruct off-node values via the **Nyström interpolation
formula**

$$
x_n(t) = \frac{1}{\lambda}\Bigl[\,y(t) +
\sum_{j=0}^n w_j(t)\, L(t, t_j)\, x_n(t_j)\,\Bigr],
\quad a \le t \le b.
\tag{4.2.70}
$$

Atkinson's exact words (p. 117): *"With this method, we approximate
those parts of the integrand in (4.2.61) that can be well-approximated
by piecewise linear interpolation, and we integrate exactly the
remaining more singular parts of the integrand."*

### 1.3 Product Simpson rule

Replace piecewise linear with **piecewise quadratic** interpolation
on subintervals $[t_{2j-2}, t_{2j}]$ of length $2h$. The 1972 paper
([1972] §III, Eqs. (16)–(19)) gives the three cardinal-function
weights $\delta_0, \delta_1, \delta_2$:

$$
\delta_0(t) = \frac{1}{2h^2}(t-h)(t-2h),\quad
\delta_1(t) = -\frac{1}{h^2}\,t(t-2h),\quad
\delta_2(t) = \delta_0(2h-t),
$$

and the singular weights are

$$
\int_{t_{2j-2}}^{t_{2j}} K(t_i, s)\,\delta_l(s - t_{2j-2})\, ds,
\quad l = 0, 1, 2.
$$

Atkinson notes (p. 122) that all the analysis of §4.2.2 generalizes
to higher-degree piecewise polynomials, with rates upgraded as in
equation (4.2.79) below.

### 1.4 The 1972 modification: hybrid product-Simpson + standard Simpson

The 1972 paper introduces an important **practical modification** for
when the kernel is singular **only** at $s = t$ but $L$ itself is
expensive to evaluate (which is exactly our Peierls situation —
$\mathrm{Ki}_1$ is not cheap):

> *Define $d_j(t_i)$ to be the distance from $t_i$ to
> $[t_{2j-2}, t_{2j}]$. Choose $\delta > 0$. On subintervals far
> from the singularity ($d_j(t_i) > \delta$) use **standard** Simpson.
> On subintervals near the singularity ($d_j(t_i) \le \delta$) use
> **product** Simpson.* ([1972] Eqs. (16), (17))

This is *exactly* the right thing for our log-singular Peierls reconstruction:
within $\delta \approx 0.2$ mfp of the diagonal, use product weights;
elsewhere, use cheap classical Simpson against the smooth tail of
$\mathrm{Ki}_1$. [1972] Eq. (17) is the master discrete equation:

$$
\lambda x_n(t_i)
- \!\!\!\sum_{d_j(t_i)\le\delta}\!\!\! \int_{t_{2j-2}}^{t_{2j}}
   K(t_i, t)\, \Lambda x(t)\, dt
\;-\!\!\sum_{d_j(t_i) > \delta}\!\!\! \Omega_j x(t_i)
= y(t_i),
\quad i = 0, \ldots, n,
\tag{[1972] (17)}
$$

where $\Lambda x(t)$ is the piecewise quadratic interpolant of $x$
and $\Omega_j x$ is the standard Simpson rule.

---

## 2. Computation of the quadrature weights

The book's §4.2.1 gives a **closed-form, $O(n)$ table-based**
algorithm for $H(t,s) = \log|t-s|$.

Define the half-weights

$$
\alpha_j(t_i) = \frac{1}{h}\int_{t_{j-1}}^{t_j}(t_j - s)\,\log|t_i - s|\, ds,
\quad
\beta_j(t_i) = \frac{1}{h}\int_{t_{j-1}}^{t_j}(s - t_{j-1})\,\log|t_i - s|\, ds,
$$

so that $w_0 = \alpha_1$, $w_n = \beta_n$, and
$w_j = \beta_j + \alpha_{j+1}$ for interior nodes. (Top of book p. 119.)

### 2.1 The $\psi$-table (key implementation trick)

Atkinson tabulates the **dimensionless** integrals

$$
\psi_0(k) = \int_0^1 \log|k - \mu|\, d\mu,
\qquad
\psi_1(k) = \int_0^1 \mu\,\log|k - \mu|\, d\mu,
\quad k = 0, \pm 1, \pm 2, \ldots
\tag{4.2.74}
$$

Both have closed forms (integrate by parts; see A&S 4.1.55). Then

$$
\boxed{
\begin{aligned}
\alpha_j(t_i) &= \frac{h}{2}\log h + h\bigl[\psi_0(i-j+1) - \psi_1(i-j+1)\bigr], \\
\beta_j(t_i)  &= \frac{h}{2}\log h + h\,\psi_1(i-j+1).
\end{aligned}
}
\tag{book p. 119}
$$

Build a table of $\psi_0(k), \psi_1(k)$ for $k = -(N-1), \ldots, N$
**once** at $O(N)$ cost. Then every weight $w_j(t_i)$ for any pair
$(i, j)$ is just two table lookups + a multiply-add. Setup cost is
$O(N)$, **not** the naïve $O(N^2)$ for $N^2$ integrations.

> **Implementation gotcha** (Atkinson p. 119, last paragraph):
> *"The standard formulas for (4.2.74) often lead to loss of
> significance errors in their calculation when $k$ becomes large,
> and care must be taken to deal with this."*
> Use the asymptotic series for large $|k|$ rather than the
> finite-difference formula for the closed form. (For us this means:
> when implementing $\psi_0, \psi_1$, branch on $|k| \gtrsim 50$ and
> use the asymptotic expansion of $\log|k-\mu|$ around $\mu = 1/2$.)

### 2.2 For other $H(t,s) = |t-s|^{\gamma-1}$ kernels

Atkinson (book p. 119, last sentence): *"It should be clear that a
similar discussion can be given for any case in which $H(t,s)$ depends
only on $|t-s|$"*. For us this means: replace $\psi_0, \psi_1$ with
their analogs for the kernel of interest. For the Bickley-Naylor case
the analog is

$$
\Psi_0^{\mathrm{Ki}_1}(k) = \int_0^1 \mathrm{Ki}_1(\Sigma h |k - \mu|)\, d\mu,
\quad
\Psi_1^{\mathrm{Ki}_1}(k) = \int_0^1 \mu\,\mathrm{Ki}_1(\Sigma h |k - \mu|)\, d\mu.
$$

These have **no closed form** but are smooth in $k$ (only the integrand
is singular as a function of $\mu$ near $\mu = k$), so they can be
tabulated once via adaptive quadrature with relative error $10^{-12}$
or so. Cost: $O(N)$ in $N$, independent of the number of nodes used
to discretise the FN equation.

---

## 3. Convergence theory

### 3.1 The clean-solution rate (4.2.78)–(4.2.79)

**Theorem 4.2.1** (book p. 120). Let $L(t,s)$ be continuous and
$H(t,s)$ satisfy

$$
c_H \equiv \sup_{a \le t \le b} \int_a^b |H(t,s)|\, ds < \infty,
\qquad
\lim_{h \to 0^+} \omega_H(h) = 0,
\tag{4.2.75, 4.2.76}
$$

where $\omega_H(h) = \sup_{|t-\tau|\le h} \int_a^b |H(t,s) - H(\tau,s)|\,ds$
(modulus of continuity of $H$ in the integrated norm). Both
$\log|t-s|$ and $|t-s|^{\gamma-1}$ ($\gamma > 0$) satisfy
(4.2.75)–(4.2.76).

If the integral equation is uniquely solvable then for sufficiently
large $n$, the discretised equation (4.2.73) is uniquely solvable
with uniformly bounded inverses, and

$$
\|x - x_n\|_\infty \le c\,\|\mathcal K x - \mathcal K_n x\|_\infty,
\quad n \ge N.
\tag{4.2.77}
$$

This is the **abstract** error bound. The explicit rates come from
(4.2.78) and (4.2.79):

- **Product trapezoidal** (book Eq. 4.2.78, requires $L(t,\cdot) \in C^2[a,b]$
  and $x \in C^2[a,b]$):
$$
\|x - x_n\|_\infty \le \frac{c\,h^2}{8}
  \max_{a \le s \le b}\Bigl|\frac{\partial^2 L(t,s)\,x(s)}{\partial s^2}\Bigr|.
$$

- **Product polynomial degree $m$** (book Eq. 4.2.79, requires
  $L(t,\cdot) \in C^{m+1}$ and $x \in C^{m+1}$):
$$
\|x - x_n\|_\infty \le c\, h^{m+1}
  \max_{a \le s \le b}\Bigl|\frac{\partial^{m+1} L(t,s)\,x(s)}{\partial s^{m+1}}\Bigr|.
$$

For product Simpson ($m=2$), this gives **at least** $O(h^3)$ — but
see §3.2.

### 3.2 The de Hoog–Weiss superconvergence (4.2.83)

Atkinson highlights an **observed** product-Simpson rate of $O(h^4 \log h)$
in his Table 4.5, **better** than the theoretical $O(h^3)$ from (4.2.79).
The improvement comes from cancellation in the product integration of
errors, as for ordinary Simpson's rule (book p. 124).

de Hoog & Weiss [book ref. 164] prove (Eq. 4.2.83): for $x \in C^4[a,b]$,

$$
\boxed{
\|\mathcal K x - \mathcal K_n x\|_\infty \le
\begin{cases}
c\,h^4\,\log h, & H(t,s) = \log|t-s|, \\
c\,h^{3+\gamma}, & H(t,s) = |t-s|^{\gamma - 1}, \quad 0 < \gamma < 1.
\end{cases}
}
\tag{4.2.83}
$$

For our Bickley-Naylor slab kernel: $H_2 = \log|x-y|$, so the
**predicted product-Simpson rate is $O(h^4 \log h)$**.

The convergence ratio in Atkinson's Table 4.5 (p. 124) is **15.6 for $h \to h/2$**
— that is, error reduction factor $\approx 16 = 2^4$ for product
Simpson on a clean test problem. This is the empirical signature
to look for in our convergence tests.

### 3.3 The reality check: solution is **not** smooth at endpoints

The product-Simpson $O(h^{4} \log h)$ rate **assumes $x \in C^4[a,b]$**.
For an integral equation with a logarithmic kernel, this is *generically false*
(Atkinson p. 125, Eq. 4.2.84):

$$
\int_0^1 \log|t-s|\, ds = t \log t + (1-t)\log(1-t) - 1, \quad 0 \le t \le 1.
$$

This $\mathcal K x_0$ is **not** in $C^1[0,1]$ even when $x_0 \equiv 1$.
The solution $x$ of $(\lambda - \mathcal K)x = y$ inherits this
endpoint singularity. **Schneider's Theorem 4.2.2** (book p. 126)
makes this precise: with $H(t,s) = g_\gamma(t-s)$ where $g_\gamma(u) =
u^{\gamma-1}$ ($0 < \gamma < 1$) or $\log|u|$ ($\gamma = 1$), the
solution behaves like $(t-a)^\gamma$ near $a$ and $(b-t)^\gamma$
near $b$. For our log case ($\gamma = 1$):

$$
|x^{(i)}(t)| \le c_i\, (t-a)^{1 - \epsilon - i},
\quad a < t \le \tfrac12(a+b),
\tag{4.2.89}
$$

for any $\epsilon \in (0,1)$. This is **fatal** to high-order
convergence on a uniform mesh.

### 3.4 Graded-mesh repair (Theorem 4.2.3, p. 132)

Atkinson's solution: **graded mesh** clustering nodes near the
endpoints. For the interval $[0,1]$ with singularity at $0$ define
nodes (Rice 1969 [book ref. 461]):

$$
t_j = \Bigl(\frac{j}{n}\Bigr)^q, \quad j = 0, 1, \ldots, n,
\tag{4.2.90}
$$

with **grading exponent** $q \ge 1$. For symmetric singularities at
both endpoints (our slab case), use double-sided grading:

$$
t_j = a + \Bigl(\frac{2j}{n}\Bigr)^q\!\cdot\frac{b-a}{2},
\quad t_{n-j} = b + a - t_j,\quad j = 0,\ldots,\tfrac{n}{2},
\quad n \text{ even}.
\tag{book p. 132}
$$

**Theorem 4.2.3** (book p. 132). Under the assumptions of Schneider's
theorem with $k = m+1$, choose $q$ to satisfy

$$
\boxed{\;
q \ge \frac{m+1}{\gamma}\;\text{ for } 0 < \gamma < 1,\qquad
q > m+1 \;\text{ for } \gamma = 1\;
}
\tag{4.2.107, 4.2.108}
$$

Then

$$
\|x - x_n\|_\infty \le \frac{c}{n^{m+1}}.
\tag{4.2.109}
$$

**For our log-kernel case ($\gamma = 1$) with product Simpson ($m = 2$)**:
the grading parameter must satisfy $q > 3$; Atkinson's example uses
$q = 4$. The empirical rates in Table 4.6 (book p. 133) confirm:

| $n$ | $q=1$ | $q=2$ | $q=3$ | $q=4$ |
|----:|------:|------:|------:|------:|
|   8 | 7.98e-2 | 3.95e-2 | 2.20e-2 | 3.35e-2 |
|  16 | 5.46e-2 | 1.91e-2 | 7.34e-3 | 9.97e-3 |
|  32 | 3.79e-2 | 9.38e-3 | 2.51e-3 | 2.80e-3 |
|  64 | 2.65e-2 | 4.64e-3 | 8.64e-4 | 7.50e-4 |
| 128 | 1.86e-2 | 2.31e-3 | 3.00e-4 | 1.95e-4 |
| 256 | 1.31e-2 | 1.15e-3 | 1.05e-4 | 4.99e-5 |
| empirical order | 0.51 | 1.01 | 1.51 | 1.97 |

**Read this carefully**: the *uniform* mesh ($q=1$) gives only
**half-order** convergence on the test problem (4.2.110)
— the worst-case behavior our slab solver is **already exhibiting** at
1–3% accuracy. The $q=4$ graded mesh recovers the full
$O(n^{-2})$ rate.

> Note: Atkinson's test in Table 4.6 uses $\gamma = 1/2$ (kernel
> $|t-s|^{-1/2}$). For our log kernel ($\gamma = 1$) the required
> grading from (4.2.108) is even harsher ($q > m+1$, not $q \ge (m+1)/\gamma$).

---

## 4. Implementation algorithm (concrete recipe)

Given: an FN-method angular flux $\psi_n$ already converged to
$\sim 10^{-5}$. Reconstruct scalar flux $\phi(x) = \int K(x,y)\,q(y)\,dy +
\text{angular contributions}$, where the integral has the slab Peierls
form.

### Step 1: Mesh

Build a graded mesh on $[a, b]$ with grading exponent $q$ chosen
according to (4.2.108):

```python
# Symmetric graded mesh on [a, b] for log-singular endpoints,
# using product-Simpson (m = 2) so grading must be q > 3; pick q = 4.
n = ...  # number of subintervals (must be even for product-Simpson)
q = 4
half = n // 2
mid = 0.5 * (a + b)
mu = np.arange(half + 1) / half   # 0, 1/half, ..., 1
left_half  = a + (mu**q) * (mid - a)
right_half = b - (mu[::-1]**q) * (b - mid)
nodes = np.concatenate([left_half, right_half[1:]])  # length n+1
```

For *interior* singularities (which we have at every node $t_i$
in the kernel diagonal $|x-y|=0$), the grading is **endpoint-only**;
the diagonal singularity is handled by the product-rule analytic
weights, **not** by mesh refinement.

### Step 2: Tabulate $\psi$ functions

```python
# For our Bickley-Naylor case, replace psi_0/psi_1 with:
def Psi_0_KiN(k, h, Sigma):
    # int_0^1 Ki_1(Sigma * h * |k - mu|) dmu
    integrand = lambda mu: Ki1(Sigma * h * abs(k - mu))
    return scipy.integrate.quad(integrand, 0, 1,
                                points=[k] if 0 <= k <= 1 else [],
                                epsrel=1e-12)[0]

def Psi_1_KiN(k, h, Sigma):
    integrand = lambda mu: mu * Ki1(Sigma * h * abs(k - mu))
    return scipy.integrate.quad(integrand, 0, 1,
                                points=[k] if 0 <= k <= 1 else [],
                                epsrel=1e-12)[0]
```

The `points=[k]` argument tells `quad` to subdivide AT the singular
abscissa even when $k \in [0,1]$ — this is the workaround for QUADPACK's
inability to detect interior endpoint singularities automatically.

Build a 2D table $\Psi_0[i, j], \Psi_1[i, j]$ of size $(n+1) \times (n+1)$
keyed on the *index difference* $k = i - j + 1$ if mesh is uniform,
or on *both indices* on a graded mesh (since $h_j$ varies).

> **Graded-mesh subtlety**: on a graded mesh the subinterval lengths
> $h_j = t_j - t_{j-1}$ vary, so the dimensionless variable is
> different on each subinterval. Replace the single $\psi(k)$ table
> by a per-subinterval recipe: for the $(i,j)$ pair, change variables
> $s = t_{j-1} + h_j \mu$ and compute
>
> $$
> \alpha_j(t_i) = \frac{1}{h_j}\int_0^{h_j}(h_j - u)\,\mathrm{Ki}_1(\Sigma|t_i - t_{j-1} - u|)\,du
> $$
>
> directly via adaptive 1D quadrature with a forced subdivision at
> $u^* = t_i - t_{j-1}$ if it lies in $[0, h_j]$.

### Step 3: Assemble the linear system

For product **Simpson** (m=2), each pair of subintervals
$[t_{2j-2}, t_{2j}]$ contributes three weights:

```python
def assemble_product_simpson_log(nodes, Sigma, L_smooth, lam):
    """
    Assemble (lam I - K_n) for the slab Peierls equation
        lam phi(t) - int K_slab(t,s) phi(s) ds = q(t)
    with K_slab = (1/2) Ki_1(Sigma |t-s|).

    Decomposition: K_slab = L_1 H_1 + L_2 H_2 with
        L_1 = 1/2,   H_1 = Ki_1_smooth(...)
        L_2 = (Sigma/2)|t-s|,  H_2 = log(...)

    L_smooth(t,s) = (1/2) [Ki_1(Sigma|t-s|) - Sigma|t-s| log(Sigma|t-s|)]
                    is in C^infty (the Taylor remainder).
    """
    n = len(nodes) - 1
    A = lam * np.eye(n + 1)
    for i, t_i in enumerate(nodes):
        for j_pair in range(0, n, 2):  # 0, 2, 4, ..., n-2
            t_lo, t_md, t_hi = nodes[j_pair], nodes[j_pair+1], nodes[j_pair+2]
            # Three product-Simpson weights against H_2 = log|t_i - s|
            # and the Lagrange basis on [t_lo, t_hi]:
            w0 = product_simpson_weight_log(t_i, t_lo, t_md, t_hi, vertex='lo')
            w1 = product_simpson_weight_log(t_i, t_lo, t_md, t_hi, vertex='md')
            w2 = product_simpson_weight_log(t_i, t_lo, t_md, t_hi, vertex='hi')
            # Plus three standard-Simpson weights against L_smooth * H_1:
            ...
            A[i, j_pair  ] -= w0 * L_singular(t_i, t_lo) + ...
            A[i, j_pair+1] -= w1 * L_singular(t_i, t_md) + ...
            A[i, j_pair+2] -= w2 * L_singular(t_i, t_hi) + ...
    return A
```

The `product_simpson_weight_log` evaluates an integral of the form

$$
\int_{t_{2j-2}}^{t_{2j}} \log|t_i - s|\,\delta_l(s - t_{2j-2})\, ds, \quad l \in \{0, 1, 2\},
$$

with $\delta_l$ the Lagrange cardinal of degree 2 — **closed forms exist**
(integrate $s^k \log|t-s|$ for $k = 0, 1, 2$), and they should be
hard-coded for numerical safety.

### Step 4: Solve and reconstruct off-node

Once $\phi_n(t_i)$ are known at the nodes, use Eq. (4.2.70) for any
off-node target:

$$
\phi_n(t) = \frac{1}{\lambda}\Bigl[\,y(t) +
\sum_j w_j(t)\, L_{\text{smooth}}(t, t_j)\,\phi_n(t_j)\,\Bigr].
$$

This is the **Nyström interpolation formula** — accurate everywhere,
not just at nodes. Atkinson notes (book p. 154): *"The Nyström
interpolation formula is quite accurate, as it generally maintains
at all points of the interval the accuracy obtained at the node
points."*

---

## 5. Where this lives in Atkinson's reference network

Atkinson's "Discussion of the literature" (book pp. 154–156) credits
the foundational papers:

- **Young (1954) [book ref. 585]**: the original product-integration
  idea (`Proc. Roy. Soc. London A 224, 561–573`).
- **de Hoog & Weiss (1973) [book ref. 164]**: improved error bounds
  for log/algebraic singular kernels — the source of (4.2.83).
- **Schneider (1979) [book ref. 493, 494]**: solution-regularity theorem
  4.2.2 and the proof that graded meshes restore high-order rates.
- **Rice (1969) [book ref. 461]**: introduced the $(j/n)^q$ graded mesh
  (4.2.90).
- **Graham (1981) [book ref. 227]**: extensions to Hölder spaces.

For our specific Bickley-Naylor kernel context:

- **Atkinson (1972)** (the local PDF) is the most directly relevant —
  it is the worked-out hybrid product/standard Simpson scheme for
  *exactly* the situation where evaluating $K$ off-diagonal is
  expensive. It is the best procedural template for our
  implementation.

- The book's §4.2 is the **clean general theory**; the 1972 paper is
  the **practical algorithm** with adaptive subdivision and "near-/far-
  field" splitting at threshold $\delta$.

Atkinson does **not** discuss Bickley-Naylor functions explicitly
anywhere in the book or paper. The notation bridge in §0 of this memo
is the connector — readers should treat it as the single load-bearing
piece of authorship in this memo. Verify it independently.

---

## 6. Variants relevant to our case

### 6.1 The general split (4.2.80)

**This is the equation the numerics-investigator should pin to the wall.**

When $K(t,s)$ does not factor cleanly as $L(t,s) H(t,s)$ with one
smooth $L$ and one singular $H$, write the kernel as a finite sum:

$$
K(t,s) = \sum_{j=1}^r L_j(t,s)\, H_j(t,s),
\tag{4.2.80}
$$

with each $L_j$ smooth and each $H_j$ from the catalog of
"weight-computable" singular functions ($\log|t-s|$, $|t-s|^{\gamma-1}$).
Atkinson (book pp. 122–123, Eq. 4.2.82): the equation
$x(t) - \int_0^\pi x(s)\log|\cos t - \cos s|\,ds = 1$ is decomposed as

$$
\log|\cos t - \cos s| = \log\Bigl|\frac{2\sin\tfrac12(t-s)\,\sin\tfrac12(t+s)}{(t-s)(t+s)(2\pi-t-s)}\Bigr|
+ \log|t-s| + \log(t+s) + \log(2\pi-t-s),
$$

i.e. one $C^\infty$ piece plus three $\log|\cdot|$ pieces. Each
$\log|\cdot|$ piece gets product-integration weights; the $C^\infty$
piece gets ordinary Simpson. **This is the pattern we will follow** for
the slab Peierls kernel: extract the $\log$ piece analytically via the
$\mathrm{Ki}_1$ Taylor expansion and run product-Simpson on it; the
remainder $\mathrm{Ki}_1(\tau) - (\pi/2) - \tau(\log\tau + \cdots)$
is $C^\infty$ and gets standard Simpson.

### 6.2 Sphere kernel ($1/|r-r'|$ prefactor)

The sphere Peierls kernel (after the angular reduction)
$\sim |r-r'|^{-1}\,\mathrm{Ki}_1(\Sigma|r-r'|)$ has an apparent
$|r-r'|^{-1}$ singularity, but the $\mathrm{Ki}_1(\tau) \sim \tau\log\tau$
small-argument behavior makes the **product** $|r-r'|^{-1} \mathrm{Ki}_1(\Sigma|r-r'|)$
go as $\log|r-r'|$ as $|r-r'|\to 0$ — i.e. exactly the same log
singularity as the slab. The decomposition (4.2.80) then writes

$$
\frac{\mathrm{Ki}_1(\Sigma|r-r'|)}{|r-r'|}
= L(r,r') \log|r-r'| + L'(r,r'),
$$

with both $L$ and $L'$ smooth. Product-trapezoidal/Simpson applies
identically, with the **same** $\psi$-table machinery.

### 6.3 Graded mesh on the slab is two-sided

Our slab problem has log-induced solution singularities at **both**
endpoints (the vacuum boundaries). Use the symmetric grading from
Atkinson p. 132 (an "even $n$" mesh with reflection about the
midpoint).

### 6.4 The eigenvalue problem (relevant to FN k_eff iteration)

Atkinson does not develop the eigenvalue case in detail in §4.2 (his
main eigenvalue treatment is in Anselone 1971 ref. [16]). However, the
Nyström-Kantorovich theory of compact-operator approximations applies
unchanged: if $\mathcal K_n \to \mathcal K$ in the *collectively
compact* sense (book §4.1.2, Theorem 4.1.2), then the eigenvalues
of $\mathcal K_n$ converge to those of $\mathcal K$ at the same rate
as the resolvent equation. Practically, the same product-Simpson
matrix $\mathcal K_n$ used in the resolvent solve is what one feeds
to a power-iteration / Arnoldi for the criticality eigenvalue.

---

## 7. Procedural summary (what numerics-investigator should do)

1. **Validate the diagnosis**: confirm the current 1–3 % stall is in
   `flux_reconstruction.py` and not in the F_N moment computation.
   The marker is: *moments converge, scalar flux saturates*. Test with
   a manufactured solution where $\phi(x)$ is known analytically;
   plot $\|\phi - \phi_n\|_\infty$ vs $n$ and confirm half-order
   (slope $\approx -1/2$ on log-log) — this is the Atkinson-Schneider
   $C^{(0,1)}$ endpoint-singularity signature.

2. **Implement product Simpson on a uniform mesh first.** Use the
   $\psi$-table approach of §2.1, adapted to $\mathrm{Ki}_1$. Re-test
   the manufactured solution. Expect: rate improves from $O(n^{-1/2})$
   to $O(n^{-?})$ where $?$ is **still bounded above by the solution
   regularity** (so still half-order on the slab, but with much
   smaller constant — Atkinson's product trapezoidal Table 4.6 column
   $q=1$ shows the same pattern: ratio 1.05 not 1.97).

3. **Add graded mesh with $q = 4$**. This is the actual fix for
   the rate. Re-test. Expect: rate $\to O(n^{-3})$ (since product Simpson
   is $m=2$, but Atkinson's de Hoog–Weiss superconvergence kicks in
   for log kernels giving $O(h^4 \log h)$ in the *operator* norm; the
   *solution* norm is bounded by the operator norm via (4.2.77),
   modulo the regularity constant).

4. **Cross-validate against an MoC reference**. The semi-analytical
   F_N moments give a structurally-independent validation. Confirm
   slab and sphere both converge to the same $\phi$ at three different
   problem sizes and three different $\Sigma$ values.

5. **Document with full equation references** (this memo is a starting
   point — the archivist will expand it into Sphinx with the
   notation-bridge derivation made fully rigorous).

### Numerical-stability gotchas to watch

- **$\psi_0, \psi_1$ table**: use the asymptotic series for
  $|k| \gtrsim 50$ to avoid catastrophic cancellation in
  $\log|k-\mu|$ for $\mu \in [0,1]$ when $k$ is large.
- **Product-Simpson weight closed forms**: integrate $s^k \log|t-s|$
  ($k = 0, 1, 2$) symbolically, factor the log, simplify;
  do **not** quadrature-evaluate at runtime.
- **Off-diagonal far-field**: cut over to standard Simpson once
  $d_j(t_i) > \delta$ where $\delta \sim 0.2$ mfp (Atkinson 1972
  empirical), to avoid evaluating $\mathrm{Ki}_1$ on $\sim n^2/2$
  points per matrix.
- **Graded-mesh corner**: at the interior breakpoint of two grading
  domains the basis is $C^0$ but **not** $C^1$. The product-Simpson
  weights are computed correctly subinterval-by-subinterval, so this
  is a non-issue for the matrix assembly, but it does mean
  Nyström interpolation (4.2.70) is only $C^0$ across the breakpoint.
  This is consistent with the underlying solution regularity and
  not a bug.

---

## Cross-references for the archivist

When this is documented in Sphinx:

- **Pillar of theory**: integral form of the transport equation for
  homogeneous slab → reduces to Fredholm second-kind with kernel
  $\frac12 \mathrm{Ki}_1(\Sigma|x-y|)$ (Lewis-Miller, Bell-Glasstone,
  Stamm'ler-Abbate Ch. III).
- **Pillar of method**: Nyström approximation (book §4.1) +
  product integration for singular kernels (book §4.2) + graded mesh
  (book §4.2.5) = the full reconstruction prescription.
- **Convergence-rate citation**: de Hoog & Weiss 1973 +
  Schneider 1979 (graded-mesh) + Atkinson 1997 book Eq. 4.2.107–4.2.109
  for the explicit rate guarantee.
- **Implementation citation**: Atkinson 1972 *Numer. Math.* 19 for the
  hybrid product/standard Simpson template with near-/far-field
  switching.

---

## Memo authored by literature-researcher

- Sources: both Atkinson references read in full from
  `/workspaces/ORPHEUS/scratch/literature/`.
- Coverage: book §4.2 (pp. 116–135), 1972 paper (all 12 pp.),
  index, and "Discussion of the literature" pp. 154–156.
- Bickley-Naylor mapping (§0 of this memo): the connector is *not in
  Atkinson*. It is constructed here from A&S 5.1.30 + the small-
  argument expansion of $\mathrm{Ki}_1$. **Verify before claiming
  it as Atkinson-cited**.
- One open question for which Atkinson is silent: whether the
  $\mathrm{Ki}_1$-specific $\Psi$ tables admit a closed-form recurrence
  the way $\psi_0, \psi_1$ do. Conjecture: yes, via integration by parts
  using the recurrence $\mathrm{Ki}_n'(\tau) = -\mathrm{Ki}_{n-1}(\tau)$.
  Recommend deferring this until empirical performance forces the
  question.
