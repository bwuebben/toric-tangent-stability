# The prism reduction and zero-sum stability

*2026-07-14, updated 2026-08-12. Derivation, proofs, and exact validation behind the "prism
structure" subsection of §7. Validation harness: `src/prism_check.py` (98/98 cases, dims 3–8,
all exact). The August updates prove every zero-sum case, record the unbalanced family, and
remove the false nonzero-sum converse, refuted already by the stable four-fold
$X(\alpha,\beta)$.*

## 1. The prism lemma

Coordinates: $N_f = \mathbb{Z}\alpha \oplus \mathbb{Z}\beta$, $\alpha=(1,0)$, $\beta=(0,1)$,
$\ell_t(x) = \langle x, t\rangle$ on $M_f$. The dP₃ hexagon is
$Q = \{x : \ell_\rho(x) \ge -1 \ \forall \text{ roots } \rho\} = \{|x_1|,|x_2|,|x_1+x_2|\le 1\}$,
area 3, six lattice-length-1 edges $E_\rho \subset \{\ell_\rho = -1\}$.

For $X(t_1..t_k)$ (dimension $n = 2+k$) the anticanonical polytope is, directly from the ray
list of Definition 7.1,
$$P = \{(x,y) : x \in Q,\ -1-\ell_{t_i}(x) \le y_i \le 1\},$$
a **generalized prism**: over each $x \in Q$ the fibre is a box with side lengths
$w_i(x) = 2 + \ell_{t_i}(x) \in [1,3]$. With $F(x) = \prod_i w_i(x)$:

- $(-K)^n = n! \int_Q F$, hence $\mu(T_X) = (n-1)!\int_Q F$;
- $\deg D_{f_i+t_i} = \deg D_{-f_i} = (n-1)!\int_Q F/w_i$ (the two base facets of factor $i$
  are identified by the unimodular shear $y_i \mapsto -y_i - \ell_{t_i}(x)$ — this also
  re-proves the degree symmetry that the fan automorphism $f_i \mapsto -f_i - t_i$ gives);
- $\deg D_\rho = (n-1)!\int_{E_\rho} F$ (lattice measure on the edge).

Hand-checked against the paper's canonical example $X(\alpha,\beta,-(\alpha+\beta))$:
$\int_Q F = 43/2$ ⟹ $(-K)^5 = 2580$ ✓; $b = 24\int F/(2+x_1) = 278$ ✓ (the "278×6" of §7 are
the six *base* degrees); edge integrals give 148 on the twist orbit $\{\alpha,\beta,-(\alpha{+}\beta)\}$
and 156 on the opposite orbit ✓ (the "148×3, 156×3").

## 2. The divergence identity and Φ

For a reflexive polygon (facets on $\langle x,\rho\rangle = -1$, primitive normals) and any
$C^1$ function $g$, applying the divergence theorem to the field $x\,g$:
$$\int_{\partial Q} g \;=\; 2\int_Q g + \int_Q \langle x, \nabla g\rangle$$
(lattice boundary measure; the Euclidean factor $|\rho|$ cancels against
$\langle x,\nu\rangle = 1/|\rho|$). With $g = F$:
$\langle x, \nabla F\rangle = F\sum_i \ell_{t_i}/(2+\ell_{t_i})$, so defining
$$\Phi(t_\bullet) := \int_Q F \sum_i \frac{\ell_{t_i}}{2+\ell_{t_i}}\,,$$
the fibre subsheaf (relative tangent of $\pi: X \to (\mathbb{P}^1)^k$, $V = N_f$) satisfies
$$\mu(\mathcal{F}_{N_f}) - \mu(T_X) = \tfrac{(n-1)!}{2}\,\Phi.$$
Validated exactly on all 98 cases. (Canonical dim-5: $\Phi = -5$ ⟹ slope gap $-60$: $456$ vs $516$ ✓.)

## 3. Proved statements (now Proposition 7.4 in the paper)

**(i) Some $t_i = 0$ ⟹ not stable.** $V = \mathbb{Q}f_i$ contains exactly the rays $\pm f_i$,
each of degree $(n-1)!\int F/2$, so $\mu(\mathcal{F}_V) = (n-1)!\int F = \mu(T_X)$ exactly.
This proves that an untwisted factor obstructs stability (equality subsheaf = relative tangent
of the $i$-th $\mathbb{P}^1$ direction — the $F_1$-type mechanism of §5).

**(ii) All $t_i$ equal to one nonzero root ⟹ not stable; unstable for $k \ge 2$.**
$\Phi = k\int_Q (2+\ell_t)^{k-1}\ell_t$; folding $Q$ along $x \mapsto -x$:
$\Phi = k\int_{Q \cap \{\ell_t > 0\}} \ell_t[(2+\ell_t)^{k-1} - (2-\ell_t)^{k-1}] \ge 0$,
strict for $k\ge2$. So the fibre tangent destabilizes. (For $k=1$: $\Phi = 0$, equality.)
This gives a simple slice of unbalanced unstable examples.

**(iii) Zero-sum nonzero roots ⟹ the fibre subsheaf is strictly subcritical.**
$u \mapsto u/(2+u)$ is strictly concave on $u > -2$; Jensen pointwise gives
$\sum_i \ell_{t_i}/(2+\ell_{t_i}) \le k\,\bar\ell/(2+\bar\ell)$ with
$\bar\ell = \tfrac1k\ell_{\sum t_i} = 0$, equality only where all $\ell_{t_i}(x)$ coincide.
A zero-sum multiset of nonzero roots contains two distinct roots, so the equality locus lies in
a proper linear subspace (measure zero) and $\Phi < 0$ strictly. **The subsheaf that makes the untwisted product
merely polystable — and a frequent destabilizer in the census (§8.1) — moves strictly below
critical.**

## 4. The full reduction (validated, mostly open)

Every ray-spanned subspace $V$ has $W = V \cap N_{f,\mathbb{C}} \in \{0, \text{root line}, N_f\}$
(single lifts $u_i^\pm$ contribute no fibre directions; a pair $\{u_i^+,u_i^-\}$ contributes
exactly $\mathbb{Q}t_i$). Hence the max subsheaf slope is attained on the finite list:
$W \in \{0, \mathbb{Q}\rho, N_f\}$ plus, per factor $i$, one of {skip, single lift ($+b_i$, dim+1),
pair ($+2b_i$, dim+1; needs $t_i \in W$ or $t_i = 0$)}. `prism_check.py` confirms this list
reproduces the harness max slope **exactly in 98/98 cases** — the "finitely many combinatorial
types uniform in $k$" claim is now a verified theorem-shaped statement.

Via the divergence identity, the $W = N_f$ family reduces to: for every sub-multiset
$T \subseteq \{t_1..t_k\}$ (the twists *not* absorbed as pairs),
$$\int_Q F \sum_{t \in T} \frac{\ell_t}{2+\ell_t} < 0.$$
$T$ = full multiset is §3(iii). $|T| = 1$ explains the observed near-critical hyperplanes:
canonical dim-5 has max slope 506 from $V = N_f \oplus \langle f_1, f_2\rangle$, and indeed
$516 - 506 = 10 = -\tfrac{4!}{4}\int F\,\ell_{t_3}/(2+\ell_{t_3}) = -6\cdot(-5/3)$ ✓.

Equality mechanism for $\sum t_i \ne 0$, seen concretely: for $(\alpha,\beta,-(\alpha+\beta),\alpha)$,
$\int_Q F_3\,\ell_\alpha = 0$ by the $\rho$-orbit symmetry of the zero-sum triple $F_3$
(the three orbit integrals are equal and sum to 0), giving slope exactly $\mu$ for
$V = N_f \oplus \langle f_1,f_2,f_3\rangle$ — matching the harness verdict "strictly ss".

**At this stage of the derivation, the open steps were:** (a) the $W = $ root-line inequalities
$\int_{E_t}F + \int_{E_{-t}}F + \ldots < \int_Q F$ (edge-vs-interior comparisons; the constant-
field divergence identities pick up slant-edge terms that don't cancel — needs a better test
field or a direct estimate); (b) the sub-multiset inequalities for general $T$ (Jensen gives
$\le |T|\,f(\ell_{\sigma_T}/|T|)$ with $\sigma_T \ne 0$, sign not pointwise-determined);
The later exact criterion below resolves (a)--(b). A general nonzero-sum obstruction is false:
$X(\alpha,\beta)$ is stable.

## 5. KE moment criterion and its solution

The barycenter of $P$: base components $\propto -\tfrac12\int_Q F\,\ell_{t_i}$, fibre components
$\propto \int_Q x F$. Since each $\ell_{t_i}$ is a combination of $\ell_\alpha, \ell_\beta$:
$$X(t_\bullet) \text{ is K\"ahler--Einstein} \iff \int_Q x\,F(x)\,dx = 0.$$
Validated on all 98 cases.  For a zero-sum multiset, use the normal form
$$F=H^d\prod_{j=1}^3(4-x_j^2)^{q_j},\qquad H=\prod_{j=1}^3(2+x_j).$$
If $d=0$, the multiset is negation-invariant and the density is centrally even.  If $d\ge1$,
the sector-ordering argument in the paper proves the sharp implication
$$q_i=\max q_\nu>\min q_\nu=q_j\quad\Longrightarrow\quad
  \int_Qx_iF<\int_Qx_jF.$$
Thus the moment can vanish only when $q_1=q_2=q_3$, which is precisely invariance under the
order-three rotation.  Consequently, among zero-sum twists, $X(t_\bullet)$ is
Kähler--Einstein exactly when the multiset is invariant under negation or rotation.  Claim 7
in `src/inequality_check.py` verifies the permutation identity and this ordering on an exact grid.

## 6. RESOLVED: the reduction criterion and zero-sum stability

The open items of §4 are now resolved; verification: `src/inequality_check.py`
(claims 1–7, all exact, ALL VERIFIED).

**Exact criterion (Prop 7.5 in the paper).** For nonzero-root twists with $X$ Fano, define
$\varphi_i = \int_Q F\ell_{t_i}/(2+\ell_{t_i})$ and line shares
$S_\rho = \int_{E_\rho}F + \int_{E_{-\rho}}F + \sum_{t_i=\pm\rho}\max(0,-\varphi_i)$. Then
**$T_X$ stable ⟺ (a) all $\varphi_i < 0$ and (b) all $S_\rho < \int_Q F$.** Key discoveries:

1. **The (A)-family collapses termwise**: $\sum_{i\in T}\varphi_i < 0$ for all sub-multisets
   $T$ ⟺ each $\varphi_i < 0$ — because each $\varphi_i$ is itself a full member of the
   family (take $A$ = complement of $\{i\}$).
2. **The share identity**: given (a), $S_\alpha + S_\beta + S_{\alpha+\beta} = 2\int_Q F$
   exactly (six edges = $\int_{\partial Q}F = 2\int F + \Phi$ by divergence; corrections
   $= -\Phi$). Hence **(b) ⟺ the three shares satisfy the strict triangle inequality** —
   the strict triangle inequality. The paper's centered-box argument proves it for every
   zero-sum multiset; the dim-7 asymmetric example has shares
   207/301, 395/602, 395/602 of $\int F$.

**Proofs achieved (Thm 7.6):**
- (a) for **negation-symmetric** multisets: pair $t$ with $-t$:
  $\varphi = -\int_Q F\ell_t^2/(4-\ell_t^2) < 0$. One line.
- (a) for **pure rotation-triples** $F = H^m$: cyclic-averaging gives
  $3\varphi = \int(3e-2p)H^{m-1}$ ($p = \sum\ell_{\gamma_j}^2$, $e = \prod\ell_{\gamma_j}$);
  split $H^{m-1} = (A+B)^{m-1}$, $A = 8-p$, $B = e$; pointwise
  $3B\Sigma_- \le 3|e|\Sigma_+ \le \frac32 p\Sigma_+$ (via $|e| \le p/2$, $|\Sigma_-|\le\Sigma_+$)
  ⟹ $3\varphi \le -\frac12\int p\,\Sigma_+ < 0$.
- (b) for **triples**: rotation permutes lines, $F$ invariant ⟹ all shares $= \frac23\int F$.
  Margin exactly $\int F/3$ (matches the observed constant-1/3 margins).
- (b) for **one-line pairs** $F = q^r$, $q = 4-\ell_\alpha^2$: slice over $\ell_\alpha = s$;
  $\int_Q F = 4K - 2L$, $S_\beta = S_{\alpha+\beta} = 2K$, $S_\alpha = 2\int F - 4K$ (share
  identity!); triangle ⟺ $0 < L < K$, trivial. Margin on the twisted line
  $= 2L_r = (4^{r+1}-3^{r+1})/(r+1) = o(\int F)$ — **asymptotically critical**, explaining
  the near-critical slopes in the census. (The K-recursion $(2r+1)K_r = 3^r + 8rK_{r-1}$
  route gives the same margin; the share identity makes it unnecessary.)

**Theorem (now in paper): $X((\alpha,-\alpha)^{\times r})$ (dim $2r+2$) and
$X((\alpha,\beta,-(\alpha+\beta))^{\times m})$ (dim $3m+2$) are smooth toric Fano, stable, KE —
all $r, m \ge 1$.** Smooth-Fano via integral prism vertices + unimodular vertex normal cones;
KE via symmetry. Stable KE examples of Picard rank $n+2$ in every even dim and every
dim ≡ 2 (mod 3). Independently corroborated: reduced criterion evaluated exactly for
$r \le 6$ (dim 14), $m \le 4$ (dim 14) — all stable, all KE (Claim 4).

**August 2026 addition: the unbalanced family.** For every $s\ge0$,
$Y_s=X((\alpha,-\alpha)^{\times s},\alpha,\beta)$ is stable and non-KE. Slicing in
$x=\ell_\alpha$ gives all singleton integrals with a negative sign and all three line shares
strictly below $\int_QF$; meanwhile $\int_QxF>0$. Claim 6 checks every displayed identity
through $s=20$, and an explicit $GL(4,\mathbb Z)$ map identifies $Y_0$ with census entry
`F.4D.0061`. Thus zero twist sum is sufficient for stability but not necessary.

**August 2026 completion of the zero-sum case.** Every zero-sum multiset has density
$H^d\prod_j(4-x_j^2)^{q_j}$. A positive layer decomposition of the pair factors reduces both
the singleton inequalities and the three line-share inequalities to centered box sections.
The singleton signs follow from a radial comparison and a sharp convolution ratio; a derivative
version of that comparison proves every line-share inequality. Hence every nonempty zero-sum
multiset of nonzero roots is stable.

**August 2026 completion of the Kähler--Einstein case.**  Central pairing isolates the odd part
of $H^d$, and a decomposition of the hexagon into twelve permutation sectors proves that the
moment belonging to a maximal pair exponent is strictly smaller than that belonging to a minimal
one.  Section 5 then gives the exact negation/rotation criterion.  The zero-sum root-twist theory
has no remaining conjectural clause; what remains open is the classification of stability and
vanishing moments for arbitrary, not necessarily zero-sum, root multisets.

## 7. Validation harness

`src/prism_check.py`: 98 cases = 8 canonical (dims 4–8, from `output/root_twist.json`) + all
k=1,2 tuples over roots∪{0} + 34 random k=3,4 tuples; checks per case: 6 hexagon degrees,
2k base degrees, total degree, reduced-list max slope == harness max slope, divergence
identity, KE-moment criterion. All exact (Fraction). 98/98.
