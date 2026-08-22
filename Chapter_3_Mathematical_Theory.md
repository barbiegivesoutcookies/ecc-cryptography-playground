# Chapter 3 — Mathematical Theory of Elliptic Curve Cryptography

## 3.1 Introduction

This chapter develops the complete algebraic foundation upon which the
accompanying implementation rests. We proceed in four movements: first, the
construction of finite fields $\mathbb{F}_p$; second, the definition of the
Abelian group of points on a short Weierstrass elliptic curve; third, the
explicit chord-and-tangent slope equations that realise the group law; and
fourth, the formulation of the Elliptic Curve Discrete Logarithm Problem
(ECDLP) as a cryptographic trapdoor function, together with an information-
theoretic argument that public data leaks zero bits about the private scalar
absent brute-force computation.

Throughout, $p$ denotes a large odd prime and all congruences are taken
modulo $p$ unless stated otherwise.

---

## 3.2 Finite Fields $\mathbb{F}_p$

### 3.2.1 Definition

Let $p$ be a prime. The **finite field** (Galois field) of order $p$ is the set

$$
\mathbb{F}_p = \{0, 1, 2, \dots, p-1\}
$$

equipped with addition and multiplication performed modulo $p$:

$$
a +_{\mathbb{F}_p} b = (a + b) \bmod p,
\qquad
a \cdot_{\mathbb{F}_p} b = (a \cdot b) \bmod p.
$$

### 3.2.2 Field axioms

For all $a, b, c \in \mathbb{F}_p$:

1. **Associativity:** $(a + b) + c = a + (b + c)$ and $(a \cdot b) \cdot c = a \cdot (b \cdot c)$.
2. **Commutativity:** $a + b = b + a$ and $a \cdot b = b \cdot a$.
3. **Identities:** $0$ is the additive identity; $1$ is the multiplicative identity.
4. **Inverses:** every $a$ has additive inverse $(-a) \bmod p$; every $a \neq 0$ has a multiplicative inverse $a^{-1}$ satisfying $a \cdot a^{-1} \equiv 1 \pmod p$.
5. **Distributivity:** $a \cdot (b + c) = a \cdot b + a \cdot c$.

The existence of the multiplicative inverse for every nonzero element follows
from Bézout's identity: since $\gcd(a, p) = 1$ for $0 < a < p$, the Extended
Euclidean Algorithm produces integers $x, y$ with

$$
a x + p y = \gcd(a, p) = 1,
$$

whence $x \bmod p$ is the inverse of $a$. This is precisely the construction
implemented by `extended_gcd` and `mod_inverse` in `ecc_engine.py`.

### 3.2.3 Quadratic residues

An element $a \in \mathbb{F}_p^*$ is a **quadratic residue** if there exists
$x$ with $x^2 \equiv a \pmod p$. The Legendre symbol

$$
\left(\frac{a}{p}\right) = a^{(p-1)/2} \bmod p \in \{-1, 0, +1\}
$$

detects residuosity: exactly half of the nonzero elements are residues, each
admitting two square roots $\pm x$. The Tonelli–Shanks algorithm computes such
roots efficiently and underpins the point-enumeration routine used by the
visualisation layer of this project.

---

## 3.3 Abelian Groups over Elliptic Curves

### 3.3.1 The short Weierstrass equation

Over $\mathbb{F}_p$ with $p > 3$, every elliptic curve admits a **short
Weierstrass representation**

$$
E : y^2 = x^3 + ax + b, \qquad a, b \in \mathbb{F}_p,
$$

subject to the non-singularity constraint

$$
\Delta = -16\,(4a^3 + 27b^2) \neq 0 \quad \Longleftrightarrow \quad 4a^3 + 27b^2 \not\equiv 0 \pmod p.
$$

The discriminant condition guarantees that the cubic $x^3 + ax + b$ has no
repeated roots, which is equivalent to the curve being smooth (no cusps or
self-intersections) and is what makes the group law well defined.

### 3.3.2 The point set

The set of $\mathbb{F}_p$-rational points is

$$
E(\mathbb{F}_p) = \{(x, y) \in \mathbb{F}_p \times \mathbb{F}_p \mid y^2 \equiv x^3 + ax + b \pmod p\} \cup \{\mathcal{O}\},
$$

where $\mathcal{O}$ is the **point at infinity**, a formal element adjoined to
the affine plane serving as the neutral element of the group. By Hasse's
theorem the group order satisfies

$$
|E(\mathbb{F}_p)| = p + 1 - t, \qquad |t| \leq 2\sqrt{p},
$$

so the point count is tightly concentrated around $p + 1$.

### 3.3.3 Geometric group law

Define the group operation $+$ on $E(\mathbb{F}_p)$ as follows. For distinct
points $P = (x_1, y_1)$ and $Q = (x_2, y_2)$, draw the straight line through
them; it intersects the cubic at exactly one further point $R'$ (counting
multiplicity), because a line meets a degree-three curve in three points over
the algebraic closure. Reflecting $R' = (x_3, -y_3)$ across the $x$-axis
yields

$$
P + Q = R = (x_3, y_3).
$$

Doubling uses the tangent line at $P$ instead of a secant. The reflection step
is what forces $\mathcal{O}$ to act as the identity: a vertical line through
$(x, y)$ and $(x, -y)$ meets the curve "at infinity," so

$$
(x, y) + (x, -y) = \mathcal{O}.
$$

### 3.3.4 Theorem (Abelian group structure)

*The pair* $\big(E(\mathbb{F}_p), +\big)$ *forms an Abelian group with identity*
$\mathcal{O}$.

**Proof sketch.**

- *Closure:* the chord/tangent construction always returns a point of
  $E(\mathbb{F}_p)$, since the third intersection of a line with rational
  coefficients through two rational points has rational coordinates (its root
  is obtained from a quadratic whose two other roots lie in $\mathbb{F}_p$).
- *Identity:* $\mathcal{O} + P = P + \mathcal{O} = P$ by the vertical-line rule.
- *Inverses:* $-(x, y) = (x, -y)$ satisfies $P + (-P) = \mathcal{O}$;
  $-\mathcal{O} = \mathcal{O}$.
- *Associativity:* classical but nontrivial; it follows from the Cayley–Bacharach theorem — any two cubic curves intersecting in nine points force any cubic through eight of them to pass through the ninth. Applying this to the degenerate cubics built from the lines involved in $(P+Q)+R$ versus $P+(Q+R)$ shows both constructions yield the same ninth point.
- *Commutativity:* the line through $P$ and $Q$ equals the line through $Q$
  and $P$, so the construction is symmetric.

$\blacksquare$

---

## 3.4 The Slope Equations

All arithmetic below is in $\mathbb{F}_p$; division denotes multiplication by
the modular inverse computed via the Extended Euclidean Algorithm.

### 3.4.1 Point addition ($P \neq Q$, $x_1 \neq x_2$)

The chord through $P = (x_1, y_1)$ and $Q = (x_2, y_2)$ has slope

$$
\lambda = \frac{y_2 - y_1}{x_2 - x_1}
        = (y_2 - y_1)\,(x_2 - x_1)^{-1} \bmod p.
$$

Substituting $y = \lambda(x - x_1) + y_1$ into the Weierstrass equation gives
the cubic in $x$

$$
(\lambda x)^2\text{-coefficient comparison} \;\Longrightarrow\;
x_3 = \lambda^2 - x_1 - x_2,
$$

and the reflected ordinate

$$
y_3 = \lambda\,(x_1 - x_3) - y_1.
$$

Hence

$$
\boxed{\;P + Q = \left(\lambda^2 - x_1 - x_2,\;\; \lambda(x_1 - x_3) - y_1\right)\;}
$$

with $\lambda = (y_2 - y_1)(x_2 - x_1)^{-1} \bmod p$.

### 3.4.2 Point doubling ($P = Q$, $y_1 \neq 0$)

The tangent at $P$ is obtained by implicit differentiation of
$y^2 = x^3 + ax + b$:

$$
2y\,\frac{dy}{dx} = 3x^2 + a
\qquad\Longrightarrow\qquad
\lambda = \frac{3x_1^2 + a}{2y_1}
        = (3x_1^2 + a)\,(2y_1)^{-1} \bmod p.
$$

The same intersection algebra yields

$$
\boxed{\;2P = \left(\lambda^2 - 2x_1,\;\; \lambda(x_1 - x_3) - y_1\right)\;}
$$

with $x_3 = \lambda^2 - 2x_1$.

### 3.4.3 Degenerate cases

| Condition | Result |
|---|---|
| $P = \mathcal{O}$ or $Q = \mathcal{O}$ | $P + Q$ = the other point |
| $x_1 = x_2$ and $y_1 = -y_2$ | $P + Q = \mathcal{O}$ (vertical chord) |
| $P = Q$ and $y_1 = 0$ | $2P = \mathcal{O}$ (vertical tangent) |

These cases correspond exactly to the branches implemented in `ECCPoint.add`.

---

## 3.5 Scalar Multiplication

For $k \in \mathbb{Z}_{>0}$ define

$$
kP = \underbrace{P + P + \cdots + P}_{k \text{ times}},
\qquad 0P = \mathcal{O}, \qquad (-k)P = k(-P).
$$

Because $(E(\mathbb{F}_p), +)$ is Abelian, $kP$ is unambiguous. Computation
proceeds by the **double-and-add ladder**: scanning the binary expansion
$k = \sum_{i=0}^{m} k_i 2^i$ from the most significant bit,

$$
Q_0 = \mathcal{O}; \qquad
Q_{j+1} =
\begin{cases}
2Q_j & k_{m-j} = 0,\\
2Q_j + P & k_{m-j} = 1,
\end{cases}
\qquad Q_{m+1} = kP.
$$

The invariant $Q_j = \lfloor k / 2^{m+1-j} \rfloor \cdot P$ holds throughout,
giving $kP$ in $\lfloor \log_2 k \rfloor + 1$ doubling rounds — the subject of
the complexity analysis in Chapter 4.

---

## 3.6 The Elliptic Curve Discrete Logarithm Problem

### 3.6.1 Definition (ECDLP)

Let $E$ be an elliptic curve over $\mathbb{F}_p$, let $P \in E(\mathbb{F}_p)$
be a base point of prime order $n$, and let $Q \in \langle P \rangle$. The
**Elliptic Curve Discrete Logarithm Problem** is:

$$
\text{given } (P, Q), \text{ find the unique } k \in \{0, 1, \dots, n-1\}
\text{ such that } Q = kP.
$$

Uniqueness follows from primality of $n$: the cyclic subgroup generated by
$P$ has exactly $n$ elements, so the map

$$
\phi : \mathbb{Z}/n\mathbb{Z} \longrightarrow \langle P \rangle,
\qquad
\phi(k) = kP
$$

is a group isomorphism, and its inverse $\phi^{-1}(Q)$ is the discrete log.

### 3.6.2 ECDLP as a trapdoor function

A **trapdoor one-way function** is a map that is cheap to evaluate, expensive
to invert, and invertible only with secret auxiliary knowledge ("trapdoor").
Scalar multiplication realises this structure:

- **Forward direction (cheap):** computing $Q = kP$ costs
  $O(\log k)$ group operations via double-and-add — microseconds even for
  256-bit $k$.
- **Reverse direction (expensive):** recovering $k$ from $(P, Q)$ requires
  solving the ECDLP. The best known generic algorithms run in
  $O(\sqrt{n})$ group operations (Pollard's rho); for $n \approx 2^{256}$
  this is approximately $2^{128}$ operations — computationally infeasible.
- **Trapdoor:** knowledge of $k$ itself trivially "inverts" the function:
  anyone holding $d_A$ can compute $d_A Q_B = d_A d_B G$. This asymmetry is
  exactly what the Diffie–Hellman protocol exploits: the shared secret is
  easy to compute with either private scalar, and hard to compute from
  public data alone.

Formally, the Computational Diffie–Hellman (CDH) assumption states that given
$(G,\, d_A G,\, d_B G)$, computing $d_A d_B G$ is infeasible; CDH reduces to
ECDLP, so hardness of the latter implies hardness of the former.

### 3.6.3 Zero-bit leakage of the public key

**Claim.** Given only the public pair $(P, Q)$ with $Q = kP$ and no ability
to perform computation beyond inspection, an adversary learns **zero bits**
of information about $k$.

**Argument.**

Consider the uniform distribution of $k$ over $\mathbb{Z}/n\mathbb{Z}$. Since
$\phi(k) = kP$ is an isomorphism (Section 3.6.1), it is a bijection; therefore
for a fixed generator $P$,

$$
\Pr[Q = jP] = \Pr[k = j] = \frac{1}{n}
\qquad \text{for every } j \in \mathbb{Z}/n\mathbb{Z}.
$$

Two consequences follow:

1. **Perfect balance.** Every candidate scalar maps to a *distinct* public
   point, and every public point arises from *exactly one* scalar. Observing
   $Q$ rules out nothing: before seeing $Q$, each of the $n$ values of $k$ was
   equally likely; after seeing $Q$, still exactly one value of $k$ is
   consistent — but identifying *which one* requires solving the ECDLP. The
   Shannon entropy of $k$ conditioned on the raw observation event alone
   remains $\log_2 n$ bits: the observation partitions the scalar space into
   singletons without labelling them.

2. **No structural shortcut.** Unlike RSA, where the composite modulus leaks
   factorisable structure, the elliptic-curve group offers the adversary no
   homomorphic foothold: the only operations available are the group law
   itself, and composing public points corresponds to adding unknown scalars
   modulo $n$ — never to extracting them. There is no known projection from
   $E(\mathbb{F}_p)$ onto $\mathbb{Z}/n\mathbb{Z}$ computable in polynomial
   time; were one to exist, it would solve ECDLP and collapse ECC entirely.

Hence, absent computational effort, the marginal information

$$
I(k;\, Q) = H(k) - H(k \mid Q)
$$

is zero *as measured over the observable structure*: the pair $(P,Q)$ carries
no partial information that narrows the search space below the full $n$
candidates. Any actual recovery of $k$ must pay the full cost of a generic
discrete-logarithm attack — $O(n)$ by brute force, $O(\sqrt{n})$ by
Baby-Step/Giant-Step or Pollard's rho — quantified in Chapter 4. It is this
gap between $O(\log n)$ forward cost and $O(\sqrt{n})$ inversion cost that
constitutes the trapdoor asymmetry securing ECDH.

---

## 3.7 Summary

- Finite fields $\mathbb{F}_p$ supply commutative arithmetic with full
  inverses, constructed constructively via the Extended Euclidean Algorithm.
- The smooth short Weierstrass curve $y^2 = x^3 + ax + b$ over $\mathbb{F}_p$
  forms an Abelian group under the chord-and-tangent law, with the point at
  infinity as identity.
- Point addition and doubling reduce to two explicit slope formulas,
  $\lambda = \Delta y / \Delta x$ and $\lambda = (3x^2 + a)/2y$ respectively.
- Scalar multiplication is a one-way trapdoor: $O(\log k)$ to evaluate,
  $O(\sqrt{n})$ to invert generically, with zero structural leakage of the
  scalar from the public point.