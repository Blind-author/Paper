# Mathematical Guarantees for the Koopman-Enhanced Objective (K-TreeROM)

**Complete proofs of the performance guarantees** referenced in Section *Performance Guarantees* of the paper. The full document is reproduced below,. Supplementary experimental results and figures are in [`RESULTS.md`](RESULTS.md).

**Contents**

- [1. Purpose of This Section](#1-purpose-of-this-section)
- [2. Definitions](#2-definitions)
- [Lemma 1 — Orthogonal Procrustes alignment](#lemma-1-orthogonal-procrustes-alignment)
- [Lemma 2 — Equivariance of reduced Koopman estimators](#lemma-2-equivariance-of-reduced-koopman-estimators)
- [Theorem 1 — The aligned Frobenius discrepancy is coordinate-invariant](#theorem-1-the-aligned-frobenius-discrepancy-is-coordinate-invariant)
- [Proposition 1 — The π/4 alignment certificate is the exact closeness boundary](#proposition-1-the-π4-alignment-certificate-is-the-exact-closeness-boundary)
- [Assumption 1 — Idealized local rollout setting](#assumption-1-idealized-local-rollout-setting)
- [Theorem 2 — Finite-horizon rollout-error upper bound](#theorem-2-finite-horizon-rollout-error-upper-bound)
- [Theorem 3 — K-TreeROM gives a bound-controlled objective for any ω > 0](#theorem-3-k-treerom-gives-a-bound-controlled-objective-for-any-ω--0)
- [Corollary 1 — Comparison with the geometry-only objective](#corollary-1-comparison-with-the-geometry-only-objective)
- [Proposition 2 — Geometry-only TreeROM can be non-identifiable](#proposition-2-geometry-only-treerom-can-be-non-identifiable)
- [Proposition 3 — The Koopman term restores a positive separation signal](#proposition-3-the-koopman-term-restores-a-positive-separation-signal)
- [3. Summary](#3-summary)
- [References](#references)

---

## 1. Purpose of This Section

This document is the complete mathematical guarantee section for the Koopman-enhanced Grassmannian regression tree objective. It incorporates the following clarifications:

1. the Procrustes result is stated as a standard alignment lemma rather than a novel contribution;
2. the coordinate-invariance result explicitly covers equivariant reduced-operator estimators, including least-squares and ridge/DMD estimators;
3. the reference frame is taken to be the node basis $\hat\Phi_m$, which matches the local tree construction;
4. the finite-horizon rollout bound includes the coordinate-mismatch term $\lVert \hat\Phi_m^T \tilde\Phi_i - I \rVert$, which is absorbed into the Grassmannian term, producing a factor of two in the geometric constant;
5. the upper-bound theorem is stated for arbitrary $\omega > 0$, with a separate special case showing when the objective is exactly proportional to the bound proxy;
6. the finite-horizon Frobenius bound is explicitly limited to the aligned-Frobenius branch $d^\star = d_F$. The spectral branch is still useful for frame-free regime separation, but it does not by itself control $\lVert \tilde A^r_i - \tilde A^r_m \rVert_F$.

## 2. Definitions

Let $\Phi_i \in \mathbb{R}^{n\times r}$ denote a local POD basis with orthonormal columns,

$$\Phi_i^T \Phi_i = I_r.$$

For a candidate node or region $R_m$, let $\hat\Phi_m \in \mathbb{R}^{n\times r}$ denote the locally aggregated POD basis, also with orthonormal columns,

$$\hat\Phi_m^T \hat\Phi_m = I_r.$$

The basis $\Phi_i$ represents a point on the Grassmann manifold $G(r,n)$, and is unique only up to a right orthogonal transformation:

$$\Phi_i \sim \Phi_i S, \qquad S \in O(r).$$

Consequently, a reduced Koopman operator is coordinate-dependent. If

$$A^r_i = \Phi_i^T A_i \Phi_i,$$

then under $\Phi_i \mapsto \Phi_i S$, one obtains

$$A^r_i \mapsto S^T A^r_i S.$$

Thus, Frobenius-norm comparisons between reduced Koopman operators are meaningful only after the operators are expressed in a common reduced coordinate system.

For the node $R_m$, set the reference basis to be $\hat\Phi_m$. Define

$$P_{i,m} = \Phi_i^T \hat\Phi_m,$$

and let

$$P_{i,m} = U_{i,m} \Sigma_{i,m} V_{i,m}^T$$

be its singular value decomposition. The Procrustes alignment from $\Phi_i$ to $\hat\Phi_m$ is

$$Q^\star_{i,m} = U_{i,m} V_{i,m}^T.$$

The aligned local basis and aligned reduced Koopman operator are

$$\tilde\Phi_i = \Phi_i Q^\star_{i,m}, \qquad \tilde A^r_i = (Q^\star_{i,m})^T A^r_i Q^\star_{i,m}.$$

The node Koopman operator, fitted directly in the node basis $\hat\Phi_m$, is denoted $\hat A^r_m$. Since it is already expressed in the node coordinates, write

$$\tilde A^r_m = \hat A^r_m.$$

The **aligned Frobenius discrepancy** is

$$d_F(i,m) = \lVert \tilde A^r_i - \tilde A^r_m \rVert_F.$$

The **spectral discrepancy** is

$$d_\sigma(A,B) = \min_{\pi \in S_r} \max_{1\le k\le r} \left| \mu_k(A) - \mu_{\pi(k)}(B) \right|,$$

where $\lbrace \mu_k(A) \rbrace_{k=1}^r$ are the eigenvalues of $A$. The **hybrid discrepancy** is

$$d^\star(i,m) = \begin{cases} \lVert \tilde A^r_i - \tilde A^r_m \rVert_F, & \sigma_{\min}(\Phi_i^T \hat\Phi_m) \ge \cos(\pi/4), \\ d_\sigma(A^r_i, \hat A^r_m), & \sigma_{\min}(\Phi_i^T \hat\Phi_m) < \cos(\pi/4). \end{cases}$$

For a candidate split $s = (j,l)$, define the **Grassmannian dispersion**

$$G(s) = \sum_{i \in R_L(s)} \delta(\Phi_i, \hat\Phi_L) + \sum_{i \in R_R(s)} \delta(\Phi_i, \hat\Phi_R),$$

and the **Koopman dispersion**

$$K(s) = \sum_{i \in R_L(s)} d^\star(i,L) + \sum_{i \in R_R(s)} d^\star(i,R).$$

The geometry-only TreeROM objective minimizes $G(s)$, whereas the Koopman-enhanced objective is

$$J_\omega(s) = \frac{G(s)}{\bar\delta_0} + \omega \frac{K(s)}{\bar d_0}, \qquad \omega > 0.$$

---

## Lemma 1 (Orthogonal Procrustes alignment)

**Statement.** Let $\Phi_i, \hat\Phi_m \in \mathbb{R}^{n\times r}$ have orthonormal columns. Then

$$Q^\star_{i,m} = \underset{Q \in O(r)}{\arg\min} \lVert \Phi_i Q - \hat\Phi_m \rVert_F$$

is given by

$$Q^\star_{i,m} = U_{i,m} V_{i,m}^T, \qquad \text{where} \qquad \Phi_i^T \hat\Phi_m = U_{i,m} \Sigma_{i,m} V_{i,m}^T.$$

If $\Phi_i^T\hat\Phi_m$ is nonsingular and the usual Procrustes nondegeneracy conditions hold, the minimizer is unique.

**Proof.** Expand the objective:

$$\lVert \Phi_i Q - \hat\Phi_m \rVert_F^2 = \mathrm{tr}\left( (\Phi_i Q - \hat\Phi_m)^T (\Phi_i Q - \hat\Phi_m) \right).$$

Since both bases are orthonormal,

$$\lVert \Phi_i Q - \hat\Phi_m \rVert_F^2 = 2r - 2\thinspace\mathrm{tr}(Q^T \Phi_i^T \hat\Phi_m).$$

Thus the problem is equivalent to maximizing

$$\mathrm{tr}(Q^T P_{i,m}), \qquad P_{i,m} = \Phi_i^T \hat\Phi_m.$$

Using the SVD $P_{i,m} = U_{i,m}\Sigma_{i,m}V_{i,m}^T$,

$$\mathrm{tr}(Q^T P_{i,m}) = \mathrm{tr}(V_{i,m}^T Q^T U_{i,m} \Sigma_{i,m}).$$

Let

$$M = V_{i,m}^T Q^T U_{i,m}.$$

Since $Q, U_{i,m}, V_{i,m} \in O(r)$, we have $M \in O(r)$. Hence

$$\mathrm{tr}(Q^T P_{i,m}) = \sum_{k=1}^r M_{kk} \sigma_k \le \sum_{k=1}^r \sigma_k,$$

because $|M_{kk}| \le 1$. Equality is achieved when $M = I_r$, which gives

$$V_{i,m}^T Q^T U_{i,m} = I_r.$$

Therefore,

$$Q = U_{i,m} V_{i,m}^T. \qquad \blacksquare$$

---

## Lemma 2 (Equivariance of reduced Koopman estimators)

**Statement.** Let $Z = \Phi_i^T X$ and $Z^+ = \Phi_i^T X^+$ denote reduced time-shifted data in the local coordinates. Consider the ridge-regularized reduced operator estimator

$$A^r(Z, Z^+) = Z^+ Z^T (Z Z^T + \gamma I_r)^{-1}, \qquad \gamma \ge 0,$$

with the Moore–Penrose interpretation when $\gamma = 0$. If the local basis is rotated by $S \in O(r)$, so that

$$Z' = S^T Z, \qquad Z^{+\prime} = S^T Z^+,$$

then

$$A^r(Z', Z^{+\prime}) = S^T A^r(Z, Z^+) S.$$

The same equivariance holds for the projected full-order operator $A^r_i = \Phi_i^T A_i \Phi_i$.

**Proof.** For $\gamma > 0$,

$$A^r(Z', Z^{+\prime}) = S^T Z^+ Z^T S \left( S^T Z Z^T S + \gamma I_r \right)^{-1}.$$

Since $S$ is orthogonal,

$$S^T Z Z^T S + \gamma I_r = S^T (Z Z^T + \gamma I_r) S,$$

and therefore

$$\left( S^T (Z Z^T + \gamma I_r) S \right)^{-1} = S^T (Z Z^T + \gamma I_r)^{-1} S.$$

Thus,

$$A^r(Z', Z^{+\prime}) = S^T Z^+ Z^T S S^T (Z Z^T + \gamma I_r)^{-1} S = S^T A^r(Z, Z^+) S.$$

For $\gamma = 0$, the same statement follows from orthogonal invariance of the Moore–Penrose pseudoinverse. For a projected full-order operator,

$$(\Phi_i S)^T A_i (\Phi_i S) = S^T (\Phi_i^T A_i \Phi_i) S. \qquad \blacksquare$$

---

## Theorem 1 (The aligned Frobenius discrepancy is coordinate-invariant)

**Statement.** Assume the reduced Koopman estimator is equivariant under right-orthogonal basis rotations, as in Lemma 2. If the local POD basis is replaced by another basis for the same subspace,

$$\Phi'_i = \Phi_i S_i, \qquad S_i \in O(r),$$

then the aligned Koopman operator is unchanged:

$$\tilde A^{r\prime}_i = \tilde A^r_i.$$

Consequently, $\lVert \tilde A^r_i - \tilde A^r_m \rVert_F$ is independent of arbitrary POD basis orientations.

**Proof.** Under the basis rotation $\Phi'_i = \Phi_i S_i$, equivariance gives

$$A^{r\prime}_i = S_i^T A^r_i S_i.$$

Moreover,

$$(\Phi'_i)^T \hat\Phi_m = S_i^T \Phi_i^T \hat\Phi_m.$$

If

$$\Phi_i^T \hat\Phi_m = U_{i,m}\Sigma_{i,m}V_{i,m}^T,$$

then

$$(\Phi'_i)^T \hat\Phi_m = (S_i^T U_{i,m}) \Sigma_{i,m} V_{i,m}^T.$$

Hence the new Procrustes factor is

$$Q'_{i,m} = S_i^T U_{i,m} V_{i,m}^T = S_i^T Q^\star_{i,m}.$$

Therefore,

$$\tilde A^{r\prime}_i = (Q'_{i,m})^T A^{r\prime}_i Q'_{i,m} = (S_i^T Q^\star_{i,m})^T (S_i^T A^r_i S_i)(S_i^T Q^\star_{i,m}) = (Q^\star_{i,m})^T S_i S_i^T A^r_i S_i S_i^T Q^\star_{i,m} = (Q^\star_{i,m})^T A^r_i Q^\star_{i,m} = \tilde A^r_i.$$

Thus the aligned operator depends only on the subspace and the node reference frame, not on the arbitrary numerical POD basis returned by the SVD. $\blacksquare$

---

## Proposition 1 (The π/4 alignment certificate is the exact closeness boundary)

**Statement.** Let $\theta_{\max}$ be the largest principal angle between $\mathrm{span}(\Phi_i)$ and $\mathrm{span}(\hat\Phi_m)$. Then

$$\sigma_{\min}(\Phi_i^T \hat\Phi_m) = \cos\theta_{\max}.$$

Moreover,

$$\theta_{\max} \le \frac{\pi}{4}$$

if and only if every aligned principal direction is at least as close to the node subspace as to its orthogonal complement. Equivalently,

$$\sigma_{\min}(\Phi_i^T \hat\Phi_m) \ge \cos\left(\frac{\pi}{4}\right) = \frac{1}{\sqrt{2}}.$$

**Proof.** Let $u \in \mathrm{span}(\Phi_i)$ be a principal direction with principal angle $\theta$. Then

$$u = \cos\theta \cdot v + \sin\theta \cdot w,$$

where

$$v \in \mathrm{span}(\hat\Phi_m), \qquad w \in \mathrm{span}(\hat\Phi_m)^\perp, \qquad \lVert v\rVert = \lVert w\rVert = 1.$$

The distance from $u$ to the node subspace is

$$\mathrm{dist}\left(u, \mathrm{span}(\hat\Phi_m)\right) = \sin\theta,$$

whereas the distance from $u$ to the orthogonal complement is

$$\mathrm{dist}\left(u, \mathrm{span}(\hat\Phi_m)^\perp\right) = \cos\theta.$$

Thus $u$ is at least as close to the node subspace as to its orthogonal complement if and only if

$$\sin\theta \le \cos\theta,$$

which is equivalent to $\theta \le \pi/4$. Applying this to every principal direction gives $\theta_{\max} \le \pi/4$.

Finally, the singular values of $\Phi_i^T \hat\Phi_m$ are

$$\sigma_k(\Phi_i^T \hat\Phi_m) = \cos\theta_k.$$

Therefore $\sigma_{\min}(\Phi_i^T \hat\Phi_m) = \cos\theta_{\max}$, and the result follows. $\blacksquare$

---

## Assumption 1 (Idealized local rollout setting)

For the finite-horizon bound below, assume the following idealized setting.

1. The sample trajectory is exactly represented in the aligned local basis:

$$x_k^{(i)} = \tilde\Phi_i z_k, \qquad \tilde\Phi_i = \Phi_i Q^\star_{i,m}.$$

2. The local and node reduced rollouts use the same initial reduced coordinate:

$$z_0 = \hat z_0, \qquad \lVert z_0 \rVert \le M.$$

3. The aligned local and node reduced operators are power-bounded by

$$\lVert \tilde A^r_i \rVert_2 \le \rho, \qquad \lVert \tilde A^r_m \rVert_2 \le \rho.$$

4. The alignment certificate holds for the sample relative to the node basis:

$$\sigma_{\min}(\Phi_i^T \hat\Phi_m) \ge \cos(\pi/4).$$

Therefore, the discrepancy branch is $d^\star = d_F$.

---

## Theorem 2 (Finite-horizon rollout-error upper bound)

**Statement.** Under Assumption 1, let

$$z_{k+1} = \tilde A^r_i z_k, \qquad \hat z_{k+1} = \tilde A^r_m \hat z_k.$$

Then, for every $k \ge 1$,

$$\left\lVert x_k^{(i)} - \hat x_k^{(m)} \right\rVert \le 2 M \rho^k \delta(\Phi_i, \hat\Phi_m) + M k \rho^{k-1} \lVert \tilde A^r_i - \tilde A^r_m \rVert_F,$$

where $\hat x_k^{(m)} = \hat\Phi_m \hat z_k$. Consequently, over a finite horizon $T$,

$$E_T^{(i,m)} \le C_\Phi \delta(\Phi_i, \hat\Phi_m) + C_A \lVert \tilde A^r_i - \tilde A^r_m \rVert_F,$$

with

$$C_\Phi = 2M \sum_{k=0}^{T} \rho^k, \qquad C_A = M \sum_{k=1}^{T} k \rho^{k-1}.$$

**Proof.** The aligned local trajectory and the node trajectory are

$$z_k = (\tilde A^r_i)^k z_0, \qquad \hat z_k = (\tilde A^r_m)^k z_0.$$

Let

$$\hat\Pi_m = \hat\Phi_m \hat\Phi_m^T, \qquad C_{i,m} = \hat\Phi_m^T \tilde\Phi_i.$$

Since $x_k^{(i)} = \tilde\Phi_i z_k$ and $\hat x_k^{(m)} = \hat\Phi_m \hat z_k$, decompose the state error as

$$x_k^{(i)} - \hat x_k^{(m)} = \tilde\Phi_i z_k - \hat\Phi_m \hat z_k = (I - \hat\Pi_m)\tilde\Phi_i z_k + \hat\Phi_m\left( \hat\Phi_m^T \tilde\Phi_i z_k - \hat z_k \right) = (I - \hat\Pi_m)\tilde\Phi_i z_k + \hat\Phi_m \left( (C_{i,m} - I) z_k + (z_k - \hat z_k) \right).$$

Taking norms and using $\lVert \hat\Phi_m \rVert_2 = 1$,

$$\left\lVert x_k^{(i)} - \hat x_k^{(m)} \right\rVert \le \lVert (I - \hat\Pi_m)\tilde\Phi_i z_k \rVert + \lVert (C_{i,m} - I) z_k \rVert + \lVert z_k - \hat z_k \rVert.$$

**First**, the orthogonal projection term is controlled by the largest principal angle:

$$\lVert (I - \hat\Pi_m)\tilde\Phi_i z_k \rVert \le \sin\theta_{\max}(\Phi_i, \hat\Phi_m) \lVert z_k \rVert.$$

Since $\lVert z_k \rVert \le M\rho^k$ and

$$\sin\theta_{\max}(\Phi_i, \hat\Phi_m) \le \left( \sum_{\ell=1}^r \theta_\ell^2 \right)^{1/2} = \delta(\Phi_i, \hat\Phi_m),$$

we have

$$\lVert (I - \hat\Pi_m)\tilde\Phi_i z_k \rVert \le M \rho^k \delta(\Phi_i, \hat\Phi_m).$$

**Second**, the coordinate mismatch term is not zero in general. The Procrustes alignment implies

$$C_{i,m} = \hat\Phi_m^T \Phi_i Q^\star_{i,m} = V_{i,m} \Sigma_{i,m} V_{i,m}^T.$$

Therefore,

$$\lVert C_{i,m} - I \rVert_2 = \max_\ell |\cos\theta_\ell - 1| = 1 - \cos\theta_{\max}.$$

For $0 \le \theta \le \pi/2$, one has $1 - \cos\theta \le \sin\theta$. Hence,

$$\lVert C_{i,m} - I \rVert_2 \le \sin\theta_{\max} \le \delta(\Phi_i, \hat\Phi_m),$$

and thus

$$\lVert (C_{i,m} - I) z_k \rVert \le M \rho^k \delta(\Phi_i, \hat\Phi_m).$$

This is the additional term that must be included; it is absorbed into the Grassmannian contribution, doubling the geometric constant.

**Third**, the dynamical term satisfies

$$\lVert z_k - \hat z_k \rVert = \left\lVert \left[ (\tilde A^r_i)^k - (\tilde A^r_m)^k \right] z_0 \right\rVert.$$

Using the telescoping identity

$$A^k - B^k = \sum_{q=0}^{k-1} A^{k-1-q}(A - B)B^{q},$$

we obtain

$$\left\lVert (\tilde A^r_i)^k - (\tilde A^r_m)^k \right\rVert_2 \le \sum_{q=0}^{k-1} \lVert \tilde A^r_i \rVert_2^{k-1-q} \lVert \tilde A^r_i - \tilde A^r_m \rVert_2 \lVert \tilde A^r_m \rVert_2^{q} \le k \rho^{k-1} \lVert \tilde A^r_i - \tilde A^r_m \rVert_2 \le k \rho^{k-1} \lVert \tilde A^r_i - \tilde A^r_m \rVert_F.$$

Therefore,

$$\lVert z_k - \hat z_k \rVert \le M k \rho^{k-1} \lVert \tilde A^r_i - \tilde A^r_m \rVert_F.$$

Combining the three terms gives

$$\left\lVert x_k^{(i)} - \hat x_k^{(m)} \right\rVert \le 2M\rho^k\delta(\Phi_i,\hat\Phi_m) + Mk\rho^{k-1}\lVert \tilde A^r_i - \tilde A^r_m \rVert_F.$$

Summing over $k = 0, \dots, T$ gives the finite-horizon result with

$$C_\Phi = 2M\sum_{k=0}^{T}\rho^k, \qquad C_A = M\sum_{k=1}^{T} k\rho^{k-1}. \qquad \blacksquare$$

> **Remark 1 (Scope of the Frobenius rollout bound).** Theorem 2 uses $\lVert \tilde A^r_i - \tilde A^r_m \rVert_F$. Therefore, it applies directly on the aligned-Frobenius branch $d^\star = d_F$, i.e., when the alignment certificate is satisfied. The spectral branch $d_\sigma$ is frame-free and useful for detecting dynamically distinct regimes, but $d_\sigma$ alone does not control the Frobenius-norm difference of nonnormal operators. Thus, the finite-horizon Frobenius rollout bound should be stated only for node/sample pairs passing the alignment certificate.

---

## Theorem 3 (K-TreeROM gives a bound-controlled objective for any ω > 0)

**Statement.** Assume all sample/node pairs involved in a candidate split satisfy the aligned-Frobenius certificate, so that $d^\star = d_F$. Then the finite-horizon split error satisfies

$$E_T(s) \le C_\Phi G(s) + C_A K(s).$$

For any $\omega > 0$, define

$$C_\omega = \max\left\lbrace C_\Phi \bar\delta_0, \frac{C_A \bar d_0}{\omega} \right\rbrace.$$

Then

$$E_T(s) \le C_\omega J_\omega(s).$$

Thus, every positive value of $\omega$ gives a valid scaled upper-bound proxy. Moreover, if

$$\omega^\star = \frac{C_A \bar d_0}{C_\Phi \bar\delta_0},$$

then

$$J_{\omega^\star}(s) = \frac{1}{C_\Phi \bar\delta_0}\left[ C_\Phi G(s) + C_A K(s) \right],$$

so minimizing $J_{\omega^\star}$ is exactly equivalent to minimizing the finite-horizon bound proxy

$$B(s) = C_\Phi G(s) + C_A K(s).$$

**Proof.** From Theorem 2 summed over all samples assigned by split $s$,

$$E_T(s) \le C_\Phi G(s) + C_A K(s).$$

Rewrite

$$G(s) = \bar\delta_0 \frac{G(s)}{\bar\delta_0}, \qquad K(s) = \frac{\bar d_0}{\omega} \cdot \omega\frac{K(s)}{\bar d_0}.$$

Therefore,

$$C_\Phi G(s) + C_A K(s) = C_\Phi \bar\delta_0 \frac{G(s)}{\bar\delta_0} + \frac{C_A \bar d_0}{\omega} \cdot \omega \frac{K(s)}{\bar d_0} \le \max\left\lbrace C_\Phi\bar\delta_0, \frac{C_A \bar d_0}{\omega} \right\rbrace \left( \frac{G(s)}{\bar\delta_0} + \omega\frac{K(s)}{\bar d_0} \right) = C_\omega J_\omega(s).$$

Thus $E_T(s) \le C_\omega J_\omega(s)$. For the special choice

$$\omega^\star = \frac{C_A \bar d_0}{C_\Phi \bar\delta_0},$$

we get

$$J_{\omega^\star}(s) = \frac{G(s)}{\bar\delta_0} + \frac{C_A \bar d_0}{C_\Phi \bar\delta_0} \cdot \frac{K(s)}{\bar d_0} = \frac{1}{C_\Phi\bar\delta_0}\left[ C_\Phi G(s) + C_A K(s) \right].$$

Since the scaling constant is positive, minimizing $J_{\omega^\star}$ is equivalent to minimizing $B(s) = C_\Phi G(s) + C_A K(s)$. $\blacksquare$

---

## Corollary 1 (Comparison with the geometry-only objective)

**Statement.** The geometry-only objective controls only $G(s)$. The Koopman-enhanced objective controls both $G(s)$ and $K(s)$. Therefore, whenever operator mismatch contributes to rollout error and $K(s)$ varies across candidate splits, the Koopman-enhanced objective is a strictly more informative error proxy than the geometry-only objective.

> **Remark 2 (How to phrase the result without overclaiming).** Theorem 3 should not be read as saying that the experimental value of $\omega$ is necessarily the exact optimal theoretical value $\omega^\star$. Rather, it shows that for any $\omega > 0$, the normalized hybrid objective upper-bounds the same two physical contributors to rollout error: spatial subspace mismatch and reduced-operator mismatch. The special value $\omega^\star$ only gives the exact proportionality between the objective and the idealized bound proxy.

---

## Proposition 2 (Geometry-only TreeROM can be non-identifiable)

**Statement.** Assume an idealized two-regime system with identical subspaces but different dynamics:

$$\Phi(\lambda) \equiv \Phi,$$

and

$$A(\lambda) = \begin{cases} A_1, & \lambda \le \lambda_s, \\ A_2, & \lambda > \lambda_s. \end{cases}$$

Then the geometry-only objective is flat:

$$G(s) = 0 \quad \text{for all candidate splits } s.$$

Thus, the true dynamical regime boundary is not identifiable from Grassmannian dispersion alone.

**Proof.** Since $\Phi_i = \Phi$ for every sample, the locally aggregated basis in any region is also $\hat\Phi_m = \Phi$, in the idealized noise-free setting. Hence

$$\delta(\Phi_i, \hat\Phi_m) = 0$$

for every sample and every region. Therefore $G(s) = 0$ for all candidate splits $s$. The geometry-only objective has no separating signal for the true boundary. $\blacksquare$

---

## Proposition 3 (The Koopman term restores a positive separation signal)

**Statement.** In the same two-regime setting, assume the operator spectra are separated by

$$g = d_\sigma(A_1, A_2) > 0.$$

For any mixed region with representative $\hat A$, the triangle inequality gives

$$\max\left\lbrace d_\sigma(A_1, \hat A), d_\sigma(A_2, \hat A) \right\rbrace \ge \frac{g}{2}.$$

Thus, any split that mixes spectrally distinct regimes pays a positive Koopman penalty, whereas a pure split does not in the idealized noise-free case.

**Proof.** By the triangle inequality for the bottleneck spectral matching distance,

$$d_\sigma(A_1, A_2) \le d_\sigma(A_1, \hat A) + d_\sigma(\hat A, A_2).$$

Since $d_\sigma(A_1, A_2) = g$,

$$g \le d_\sigma(A_1, \hat A) + d_\sigma(A_2, \hat A).$$

Therefore, at least one of the two terms must be at least $g/2$:

$$\max\left\lbrace d_\sigma(A_1, \hat A), d_\sigma(A_2, \hat A) \right\rbrace \ge \frac{g}{2}.$$

This gives a positive dynamical separation signal on mixed regions. $\blacksquare$

---

## 3. Summary

The theoretical results above justify the Koopman-enhanced splitting criterion as follows. The Procrustes alignment places locally computed reduced operators into the node coordinate frame, and the resulting aligned Frobenius discrepancy is invariant to arbitrary right-orthogonal rotations of the POD basis. The certificate $\sigma_{\min}(\Phi_i^T\hat\Phi_m) \ge \cos(\pi/4)$ is the exact condition under which every aligned principal direction remains at least as close to the node subspace as to its orthogonal complement. Under this certificate and bounded reduced dynamics, the finite-horizon rollout error is bounded by a spatial term proportional to $\delta(\Phi_i, \hat\Phi_m)$ and a dynamical term proportional to $\lVert \tilde A^r_i - \tilde A^r_m \rVert_F$. Hence the proposed objective controls both subspace projection mismatch and accumulated reduced-dynamics mismatch, while the original Grassmannian objective controls only the spatial component. When subspaces are identical or nearly identical but the reduced dynamics differ, the geometry-only objective can be non-identifiable, whereas the Koopman term provides a positive separation signal proportional to the spectral gap.

## References

1. P. H. Schönemann, "A generalized solution of the orthogonal Procrustes problem," *Psychometrika*, vol. 31, pp. 1–10, 1966.
2. R. Bhatia, *Matrix Analysis*. Springer, 1997.
3. C. Davis and W. M. Kahan, "The rotation of eigenvectors by a perturbation. III," *SIAM Journal on Numerical Analysis*, vol. 7, no. 1, pp. 1–46, 1970.
