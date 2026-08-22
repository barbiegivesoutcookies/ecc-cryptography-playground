# Chapter 4 — Performance Analysis of Elliptic Curve Scalar Multiplication

## 4.1 Introduction

This chapter provides a formal computer-science analysis of the computational
costs that govern elliptic curve cryptography. The central object of study is
the **double-and-add ladder** used to evaluate scalar multiplication $kP$,
whose running time we prove to be $\Theta(\log k)$ group operations. We then
contrast this efficient forward direction against the asymptotic costs of the
cryptanalytic attacks that attempt to invert it — exhaustive brute force,
Baby-Step/Giant-Step (BSGS), and Pollard's rho — and close with an evaluation
table mapping key-size dimensions to memory consumption scales.

Cost model: throughout, one **group operation** denotes either a point
addition or a point doubling on $E(\mathbb{F}_p)$, each comprising a constant
number of modular multiplications and one modular inversion (or, in projective
implementations, a fixed larger batch of multiplications). Bit operations on
$\ell$-bit integers cost $O(\ell)$ to $O(\ell^2)$ depending on multiplication
strategy; since $\ell$ is fixed by the security parameter, we count group
operations as the atomic unit.

---

## 4.2 The Double-and-Add Ladder

### 4.2.1 Algorithm

Given a point $P$ and a positive integer $k$ with binary expansion

$$
k = \sum_{i=0}^{m} k_i \, 2^i, \qquad k_m = 1, \qquad m = \lfloor \log_2 k \rfloor,
$$

the most-significant-bit-first ladder computes:

```
Q ← O                      # accumulator
for j = m down to 0:
    Q ← 2Q                 # doubling step
    if k_j = 1:
        Q ← Q + P          # addition step
return Q
```

### 4.2.2 Correctness invariant

**Lemma 4.1.** *At the start of iteration* $j$ *(processing bit* $k_j$*)*, the
accumulator satisfies

$$
Q = \left\lfloor \frac{k}{2^{j+1}} \right\rfloor P.
$$

**Proof.** By induction on the iteration index. Initially ($j = m$),
$Q = \mathcal{O} = 0 \cdot P$ and $\lfloor k / 2^{m+1} \rfloor = 0$ since
$k < 2^{m+1}$. Inductive step: entering iteration $j$, the loop doubles the
accumulator and conditionally adds $P$:

$$
Q' = 2\left\lfloor \frac{k}{2^{j+1}} \right\rfloor P + k_j P
   = \left(2\left\lfloor \frac{k}{2^{j+1}} \right\rfloor + k_j\right) P.
$$

By the identity $2\lfloor x/2 \rfloor + (x \bmod 2) = x$ applied to
$x = \lfloor k / 2^{j} \rfloor$,

$$
2\left\lfloor \frac{k}{2^{j+1}} \right\rfloor + k_j = \left\lfloor \frac{k}{2^{j}} \right\rfloor,
$$

which re-establishes the invariant for iteration $j-1$. After the final
iteration ($j = 0$), $Q = \lfloor k / 1 \rfloor P = kP$. $\blacksquare$

### 4.2.3 Time complexity theorem

**Theorem 4.2.** *Double-and-add computes* $kP$ *in exactly* $m$ doublings
*and between* $0$ *and* $m$ *additions, i.e.,*

$$
T(k) \in \Theta(\log k) \text{ group operations}.
$$

**Proof.**

*Upper bound.* The loop body executes once per bit of $k$, i.e., $m + 1 =
\lfloor \log_2 k \rfloor + 1$ times. Each execution performs at most two
group operations (one doubling, one conditional addition). Hence

$$
T(k) \leq 2(\lfloor \log_2 k \rfloor + 1) \in O(\log k).
$$

*Lower bound.* Every iteration performs at least one doubling, so

$$
T(k) \geq \lfloor \log_2 k \rfloor + 1 \in \Omega(\log k).
$$

Combining both bounds yields $T(k) \in \Theta(\log k)$. More precisely, if
$k$ has Hamming weight $w(k)$ (the number of set bits), then

$$
T(k) = (\lfloor \log_2 k \rfloor + 1) + w(k),
$$

and for uniformly random $k$ the expected weight is $\tfrac{1}{2}\log_2 k$,
giving an average-case cost of approximately $1.5 \log_2 k$ group
operations. $\blacksquare$

### 4.2.4 Contrast with naive repeated addition

The naive approach $P + P + \cdots + P$ performs exactly $k - 1$ additions:

$$
T_{\text{naive}}(k) = k - 1 \in \Theta(k).
$$

Since $\Theta(\log k) \subsetneq \Theta(k)$, the ladder is exponentially
faster. For secp256k1 with $k \approx 2^{256}$, the naive method would need
on the order of $10^{77}$ group operations, whereas the ladder needs roughly
$384$ — a speedup factor of about $2^{256} / (1.5 \cdot 256) \approx
2^{248}$. This exponential gap is precisely what makes the forward direction
of the trapdoor function cheap.

---

## 4.3 Cryptanalytic Attack Complexities

Let $n$ denote the prime order of the base point $G$; the ECDLP instance is
$(G,\, Q = kG)$ with $k \in [0, n)$.

### 4.3.1 Brute force — $O(n)$

Exhaustively compute $G, 2G, 3G, \dots$ until $Q$ appears. Each candidate
requires one group operation (incremental addition), so the worst case costs

$$
T_{\text{BF}}(n) = n - 1 \in O(n)
$$

group operations and $O(1)$ memory. For $n \approx 2^{256}$ this is
$\approx 1.16 \times 10^{77}$ operations — utterly infeasible: at $10^{18}$
operations per second (exascale), the search would take roughly
$3.7 \times 10^{51}$ years, versus $\approx 1.4 \times 10^{10}$ years for the
age of the universe.

### 4.3.2 Baby-Step Giant-Step — $O(\sqrt{n})$ time, $O(\sqrt{n})$ memory

Shanks' BSGS exploits the birthday-style decomposition $k = im + j$ with
$m = \lceil \sqrt{n} \rceil$, $0 \leq i, j < m$:

1. **Baby steps:** tabulate $jG$ for all $j \in [0, m)$ — $m$ group
   operations and $O(m)$ stored points.
2. **Giant steps:** compute $Q - imG$ for $i = 0, 1, \dots$ until a match is
   found in the table — at most $m$ further group operations.

By the pigeonhole principle some pair $(i, j)$ must satisfy $im + j = k$
since $k < n \leq m^2$, so the algorithm terminates within

$$
T_{\text{BSGS}}(n) \leq 2\lceil \sqrt{n} \rceil \in O(\sqrt{n})
$$

group operations, but requires storing $\Theta(\sqrt{n})$ group elements.
For $n \approx 2^{256}$: $\approx 2^{128}$ operations *and*
$\approx 2^{128} \times 32$ bytes $\approx 10^{39}$ bytes of memory —
exceeding the estimated total data capacity of civilization by dozens of
orders of magnitude. BSGS is therefore memory-bound long before it is
time-bound.

### 4.3.3 Pollard's rho — $O(\sqrt{n})$ time, $O(1)$ memory

Pollard's rho replaces BSGS's table with pseudorandom iteration inside the
finite group. Define a partition of $E(\mathbb{F}_p)$ into three classes and
an iterating function

$$
f(R) =
\begin{cases}
G + R & R \in S_0,\\
2R    & R \in S_1,\\
Q + R & R \in S_2,
\end{cases}
$$

tracking affine coefficients $(a_i, b_i)$ with $R_i = a_i G + b_i Q$
throughout. Iteration from $R_0$ eventually enters a cycle ("rho" shape);
by Floyd's cycle detection (tortoise-and-hare) or Brent's variant, a
collision $R_i = R_j$ is found with matching coefficient pairs:

$$
a_i G + b_i Q = a_j G + b_j Q
\;\Longrightarrow\;
(a_i - a_j)\,G = (b_j - b_i)\,k\,G
\;\Longrightarrow\;
k = (a_i - a_j)(b_j - b_i)^{-1} \bmod n.
$$

Because the iteration behaves as a random walk on a set of size $n$, the
birthday paradox guarantees a collision after an expected

$$
\mathbb{E}[T_{\rho}(n)] = \sqrt{\frac{\pi n}{2}} \approx 1.2533\,\sqrt{n}
\in O(\sqrt{n})
$$

iterations, using only $O(1)$ memory. Parallelisation across $r$ processors
yields a linear speedup to $O(\sqrt{n}/r)$ (van Oorschot–Wiener), which is
why Pollard's rho — not brute force or BSGS — defines the true generic
security level of ECC: an $m$-bit curve offers approximately $m/2$ bits of
generic security.

### 4.3.4 Comparative summary

| Method | Time (group ops) | Memory | Generic? | Notes |
|---|---|---|---|---|
| Double-and-add (forward) | $\Theta(\log k)$ | $O(1)$ | — | Legitimate evaluation path |
| Naive addition (forward) | $\Theta(k)$ | $O(1)$ | — | Pedagogical baseline only |
| Brute force | $O(n)$ | $O(1)$ | Yes | Dominated by rho; never optimal |
| Baby-Step Giant-Step | $O(\sqrt{n})$ | $O(\sqrt{n})$ | Yes | Deterministic; memory-infeasible at scale |
| Pollard's rho | $O(\sqrt{n})$ expected | $O(1)$ | Yes | Best known generic attack |

The asymmetry is stark: the legitimate user pays $O(\log n)$; the attacker
pays $O(\sqrt{n})$ — a factor of $\sqrt{n}/\log n$. For $n \approx 2^{256}$,
that ratio is approximately $2^{128}/256 = 2^{120}$.

---

## 4.4 Key Size Dimensions and Memory Consumption

Security level is governed by the bit length $b$ of the group order $n$: the
attacker's work is $\approx 2^{b/2}$ while storage per point is proportional
to the coordinate size $\approx b$ bits. The table below maps representative
key dimensions to attack-work and memory scales.

| Curve class | Order size $b$ (bits) | Forward cost $\approx 1.5b$ ops | Attack work $\approx 2^{b/2}$ ops | BSGS memory (32-byte points) | Security level |
|---|---|---|---|---|---|
| Toy visualisation (this project) | ~7 ($p = 97$) | ~10 ops | ~$2^{3.5} \approx 11$ ops | ~352 bytes | None (demonstration) |
| Legacy toy | 64 | ~96 ops | $2^{32} \approx 4.3 \times 10^{9}$ ops | ~137 GB | Broken trivially |
| Deprecated | 112 | ~168 ops | $2^{56} \approx 7.2 \times 10^{16}$ ops | ~$2.9 \times 10^{18}$ GB | Insufficient |
| Standard | 192 | ~288 ops | $2^{96} \approx 7.9 \times 10^{28}$ ops | ~$3.2 \times 10^{30}$ GB | 96-bit |
| **secp256k1 / P-256** | **256** | **~384 ops** | **$2^{128} \approx 3.4 \times 10^{38}$ ops** | **~$5.4 \times 10^{39}$ GB** | **128-bit** |
| High security | 384 | ~576 ops | $2^{192} \approx 6.3 \times 10^{57}$ ops | ~$10^{60}$ GB | 192-bit |
| Extreme | 512 | ~768 ops | $2^{256} \approx 1.2 \times 10^{77}$ ops | ~$10^{79}$ GB | 256-bit |

Reading the table: each additional bit of key length squares the attacker's
work while adding only a constant amount to the legitimate user's cost —
the defining economic property of exponential-security cryptosystems. Note
also that the BSGS memory column grows as $2^{b/2} \cdot b/8$ bytes, which is
why even moderately sized legacy parameters fall to storage exhaustion before
compute exhaustion.

---

## 4.5 Implementation-Level Considerations

The asymptotic analysis assumes constant-time group operations. Two practical
caveats apply to the reference implementation in `ecc_engine.py`:

1. **Affine inversion cost.** Each `add` call invokes `mod_inverse`
   (Extended Euclidean). Its cost depends on operand bit length, which is
   bounded by $p$'s bit length — constant per curve — so the
   $\Theta(\log k)$ operation count stands. Production libraries avoid
   inversions entirely via projective (Jacobian) coordinates.
2. **Side channels.** The textbook double-and-add branches on secret bits
   ($k_j$), enabling timing/power side-channel attacks. Constant-time ladders
   (Montgomery ladder, always-double-and-add) perform identical operation
   sequences regardless of $k$; this hardening is out of scope for the
   educational codebase but mandatory in production.

Neither caveat alters the asymptotic conclusions of Sections 4.2–4.3.

---

## 4.6 Summary

- Double-and-add evaluates $kP$ in exactly $\lfloor \log_2 k \rfloor + 1$
  doublings plus $w(k)$ additions: $\Theta(\log k)$ total, proven via the
  invariant $Q_j = \lfloor k/2^{j+1} \rfloor P$.
- Generic ECDLP attacks require $O(n)$ (brute force) or $O(\sqrt{n})$
  (BSGS, Pollard's rho) operations; Pollard's rho dominates with $O(1)$
  memory and defines the generic security level $b/2$ bits for an $b$-bit
  order.
- The defender/attacker cost ratio $\sqrt{n}/\log n$ reaches roughly
  $2^{120}$ at the 256-bit secp256k1 scale, underpinning ECC's practical
  security with minimal key sizes.