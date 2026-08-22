"""
ecc_engine.py
=============

A library-free implementation of Elliptic Curve Cryptography (ECC) over prime
finite fields F_p.

This module provides:

1. Arbitrary-precision modular arithmetic primitives built directly on top of
   Python's native ``int`` type (which is itself arbitrary precision).
2. The Extended Euclidean Algorithm and modular inversion.
3. An ``ECCPoint`` class modelling points on a short Weierstrass curve:

       y^2 = x^3 + a*x + b  (mod p)

   including the neutral element, the Point at Infinity.
4. Algebraic group operations: point addition, point doubling, negation,
   scalar multiplication via the Double-and-Add ladder in O(log k) steps.
5. A complete Elliptic Curve Diffie-Hellman (ECDH) key exchange protocol.

No third-party cryptography libraries are used anywhere in this module.

Author: Senior Cryptographer / Software Architect
License: MIT
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import secrets
from typing import ClassVar, Dict, Final, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Modular arithmetic primitives (arbitrary-precision integer arithmetic)
# ---------------------------------------------------------------------------


def mod_add(a: int, b: int, p: int) -> int:
    """Compute ``(a + b) mod p`` for arbitrary-precision integers.

    Args:
        a: First addend. May be any integer (negative values allowed).
        b: Second addend. May be any integer (negative values allowed).
        p: The modulus. Must be a positive integer.

    Returns:
        The residue of ``a + b`` in the canonical range ``[0, p)``.

    Raises:
        ValueError: If ``p`` is not a positive integer.
    """
    if p <= 0:
        raise ValueError("Modulus p must be a positive integer.")
    return (a + b) % p


def mod_sub(a: int, b: int, p: int) -> int:
    """Compute ``(a - b) mod p`` for arbitrary-precision integers.

    Args:
        a: Minuend. May be any integer.
        b: Subtrahend. May be any integer.
        p: The modulus. Must be a positive integer.

    Returns:
        The residue of ``a - b`` in the canonical range ``[0, p)``.

    Raises:
        ValueError: If ``p`` is not a positive integer.
    """
    if p <= 0:
        raise ValueError("Modulus p must be a positive integer.")
    return (a - b) % p


def mod_mul(a: int, b: int, p: int) -> int:
    """Compute ``(a * b) mod p`` for arbitrary-precision integers.

    Args:
        a: First factor. May be any integer.
        b: Second factor. May be any integer.
        p: The modulus. Must be a positive integer.

    Returns:
        The residue of ``a * b`` in the canonical range ``[0, p)``.

    Raises:
        ValueError: If ``p`` is not a positive integer.
    """
    if p <= 0:
        raise ValueError("Modulus p must be a positive integer.")
    return (a * b) % p


def mod_pow(base: int, exponent: int, p: int) -> int:
    """Compute ``base**exponent mod p`` using fast modular exponentiation.

    Python's built-in three-argument ``pow`` implements square-and-multiply
    internally on arbitrary-precision integers; it is used here as the
    canonical primitive (no third-party libraries involved).

    Args:
        base: The base value. May be negative; results are normalised.
        exponent: A non-negative integer exponent.
        p: The modulus. Must be a positive integer.

    Returns:
        The residue of ``base**exponent`` in the canonical range ``[0, p)``.

    Raises:
        ValueError: If ``p`` is not positive or the exponent is negative.
    """
    if p <= 0:
        raise ValueError("Modulus p must be a positive integer.")
    if exponent < 0:
        raise ValueError("Exponent must be non-negative.")
    return pow(base % p, exponent, p)


def extended_gcd(a: int, b: int) -> Tuple[int, int, int]:
    """Extended Euclidean Algorithm.

    Iteratively computes integers ``(g, x, y)`` satisfying Bezout's identity:

        a*x + b*y = g = gcd(a, b)

    The iterative formulation avoids Python recursion limits even for very
    large operands (e.g., 256-bit cryptographic primes).

    Args:
        a: First integer operand.
        b: Second integer operand.

    Returns:
        A tuple ``(g, x, y)`` where ``g = gcd(a, b)`` and
        ``a*x + b*y == g``.
    """
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1

    while r != 0:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_s, s = s, old_s - quotient * s
        old_t, t = t, old_t - quotient * t

    # Normalise so that gcd is non-negative.
    if old_r < 0:
        old_r = -old_r
        old_s = -old_s
        old_t = -old_t

    return old_r, old_s, old_t


def mod_inverse(n: int, p: int) -> int:
    """Compute the multiplicative inverse of ``n`` modulo ``p``.

    Uses the Extended Euclidean Algorithm to find ``x`` such that:

        n * x === 1 (mod p)

    This requires ``gcd(n, p) == 1``; when ``p`` is prime this holds for
    every ``n`` not congruent to zero.

    Args:
        n: The integer to invert. May be negative.
        p: The modulus. Must be a positive integer.

    Returns:
        The inverse of ``n`` modulo ``p`` in the range ``[0, p)``.

    Raises:
        ValueError: If ``p <= 0`` or if ``n`` is not invertible mod ``p``.
    """
    if p <= 0:
        raise ValueError("Modulus p must be a positive integer.")

    g, x, _ = extended_gcd(n % p, p)
    if g != 1:
        raise ValueError(
            f"{n} has no multiplicative inverse modulo {p} "
            f"(gcd = {g})."
        )
    return x % p


def legendre_symbol(a: int, p: int) -> int:
    """Compute the Legendre symbol ``(a / p)`` for an odd prime ``p``.

    Returns:
        +1 if ``a`` is a quadratic residue mod ``p``,
        -1 if ``a`` is a quadratic non-residue,
         0 if ``a === 0 (mod p)``.
    """
    if p <= 2:
        raise ValueError("Legendre symbol requires an odd prime modulus.")
    residue: int = mod_pow(a, (p - 1) // 2, p)
    if residue == 0:
        return 0
    if residue == 1:
        return 1
    return -1


def tonelli_shanks_sqrt(a: int, p: int) -> Optional[int]:
    """Compute a square root of ``a`` modulo an odd prime ``p``.

    Implements the full Tonelli-Shanks algorithm, which handles every odd
    prime including those where ``p === 3 (mod 4)`` (fast path) and
    ``p === 1 (mod 4)`` (general path).

    Args:
        a: The value whose square root is sought.
        p: An odd prime modulus.

    Returns:
        One square root of ``a`` mod ``p``, or ``None`` if ``a`` is a
        quadratic non-residue.
    """
    a %= p
    if a == 0:
        return 0
    if p == 2:
        return a
    if legendre_symbol(a, p) != 1:
        return None

    # Fast path: p === 3 (mod 4) => root = a^((p+1)/4).
    if p % 4 == 3:
        return mod_pow(a, (p + 1) // 4, p)

    # General Tonelli-Shanks: write p - 1 = q * 2^s with q odd.
    q: int = p - 1
    s: int = 0
    while q % 2 == 0:
        q //= 2
        s += 1

    # Find a quadratic non-residue z by trial over small candidates.
    z: int = 2
    while legendre_symbol(z, p) != -1:
        z += 1

    m: int = s
    c: int = mod_pow(z, q, p)
    t: int = mod_pow(a, q, p)
    r: int = mod_pow(a, (q + 1) // 2, p)

    while t != 1:
        # Find the least i in [1, m) with t^(2^i) === 1.
        i: int = 0
        t2i: int = t
        while t2i != 1:
            t2i = mod_mul(t2i, t2i, p)
            i += 1
            if i >= m:
                raise ArithmeticError(
                    "Tonelli-Shanks failed: input p may not be prime."
                )
        b: int = c
        for _ in range(m - i - 1):
            b = mod_mul(b, b, p)
        m = i
        c = mod_mul(b, b, p)
        t = mod_mul(t, c, p)
        r = mod_mul(r, b, p)

    return r


def is_on_curve_x(x: int, curve: "EllipticCurve") -> bool:
    """Check whether ``x`` admits at least one valid point on the curve."""
    rhs: int = curve.evaluate_rhs(x)
    return legendre_symbol(rhs, curve.p) != -1


# ---------------------------------------------------------------------------
# Curve definition
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class EllipticCurve:
    """A short Weierstrass elliptic curve over a prime field F_p.

    The curve equation is::

        y^2 = x^3 + a*x + b   (mod p)

    subject to the non-singularity constraint ``4*a^3 + 27*b !== 0 (mod p)``.

    Attributes:
        a: The linear coefficient of the Weierstrass equation.
        b: The constant term of the Weierstrass equation.
        p: The prime field characteristic.
        gx: The x-coordinate of the standard generator (base point), if any.
        gy: The y-coordinate of the standard generator (base point), if any.
        n: The order of the base point (group order), if known.
        name: A human-readable identifier for the curve.
    """

    a: int
    b: int
    p: int
    gx: Optional[int] = None
    gy: Optional[int] = None
    n: Optional[int] = None
    name: str = "custom"

    def __post_init__(self) -> None:
        """Validate the curve parameters upon construction."""
        if self.p <= 2:
            raise ValueError("Field characteristic p must be an odd prime.")
        discriminant: int = (
            4 * pow(self.a, 3, self.p) + 27 * pow(self.b, 2, self.p)
        ) % self.p
        if discriminant == 0:
            raise ValueError(
                "Curve parameters are singular: 4a^3 + 27b^2 === 0 (mod p)."
            )

    @property
    def generator(self) -> "ECCPoint":
        """Return the standard base point as an :class:`ECCPoint`."""
        if self.gx is None or self.gy is None:
            raise ValueError(f"Curve '{self.name}' defines no base point.")
        return ECCPoint(self.gx, self.gy, self)

    def evaluate_rhs(self, x: int) -> int:
        """Evaluate the right-hand side ``x^3 + a*x + b (mod p)``."""
        x_mod: int = x % self.p
        return (pow(x_mod, 3, self.p)
                + self.a * x_mod + self.b) % self.p

    def contains_point(self, x: int, y: int) -> bool:
        """Return True iff ``(x, y)`` satisfies the curve equation."""
        lhs: int = pow(y % self.p, 2, self.p)
        return lhs == self.evaluate_rhs(x)


# ---------------------------------------------------------------------------
# Standardised curves
# ---------------------------------------------------------------------------

SECP256K1_P: Final[int] = (
    0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
)
SECP256K1_A: Final[int] = 0
SECP256K1_B: Final[int] = 7
SECP256K1_GX: Final[int] = (
    0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
)
SECP256K1_GY: Final[int] = (
    0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8
)
SECP256K1_N: Final[int] = (
    0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
)

CURVE_SECP256K1: Final[EllipticCurve] = EllipticCurve(
    a=SECP256K1_A,
    b=SECP256K1_B,
    p=SECP256K1_P,
    gx=SECP256K1_GX,
    gy=SECP256K1_GY,
    n=SECP256K1_N,
    name="secp256k1",
)

# A tiny toy curve used for visualisation: y^2 = x^3 + 2x + 3 (mod 97).
TOY_CURVE_P: Final[int] = 97
TOY_CURVE_A: Final[int] = 2
TOY_CURVE_B: Final[int] = 3

CURVE_TOY_97: Final[EllipticCurve] = EllipticCurve(
    a=TOY_CURVE_A,
    b=TOY_CURVE_B,
    p=TOY_CURVE_P,
    gx=None,
    gy=None,
    n=None,
    name="toy-p97",
)


# ---------------------------------------------------------------------------
# Point representation and group law
# ---------------------------------------------------------------------------


class ECCPoint:
    """An affine point on a short Weierstrass elliptic curve over F_p.

    The Point at Infinity (the group identity element) is represented by the
    singleton flag ``is_infinity=True`` with sentinel coordinates ``(0, 0)``.

    Attributes:
        x: The affine x-coordinate (meaningless for infinity).
        y: The affine y-coordinate (meaningless for infinity).
        curve: The parent :class:`EllipticCurve`.
        is_infinity: True iff this instance is the Point at Infinity.
    """

    __slots__ = ("x", "y", "curve", "is_infinity")

    #: Sentinel coordinates used by the Point at Infinity.
    INFINITY_COORDS: ClassVar[Tuple[int, int]] = (0, 0)

    def __init__(
        self,
        x: Optional[int],
        y: Optional[int],
        curve: EllipticCurve,
        is_infinity: bool = False,
    ) -> None:
        """Construct a curve point.

        Args:
            x: Affine x-coordinate, or ``None`` for infinity.
            y: Affine y-coordinate, or ``None`` for infinity.
            curve: The parent curve.
            is_infinity: Explicitly mark this point as the identity element.

        Raises:
            ValueError: If finite coordinates do not satisfy the curve
                equation.
        """
        self.curve: EllipticCurve = curve
        if is_infinity or x is None or y is None:
            self.is_infinity: bool = True
            self.x: int = ECCPoint.INFINITY_COORDS[0]
            self.y: int = ECCPoint.INFINITY_COORDS[1]
            return

        x_val: int = x % curve.p
        y_val: int = y % curve.p
        if not curve.contains_point(x_val, y_val):
            raise ValueError(
                f"Point ({x_val}, {y_val}) does not lie on curve "
                f"'{curve.name}': y^2 = x^3 + {curve.a}x + {curve.b} "
                f"(mod {curve.p})."
            )
        self.is_infinity = False
        self.x = x_val
        self.y = y_val

    # -- Factory helpers ----------------------------------------------------

    @staticmethod
    def infinity(curve: EllipticCurve) -> "ECCPoint":
        """Return the Point at Infinity (identity) on ``curve``."""
        return ECCPoint(None, None, curve, is_infinity=True)

    @classmethod
    def from_y_parity(cls, x: int, y_root: int, curve: EllipticCurve,
                      parity: str) -> "ECCPoint":
        """Construct a point from an x-coordinate and a chosen y parity.

        Args:
            x: The x-coordinate.
            y_root: Any square root of the curve RHS at ``x``.
            curve: The parent curve.
            parity: Either ``"even"`` or ``"odd"`` selecting which of the two
                mirror points to construct.

        Returns:
            The selected :class:`ECCPoint`.

        Raises:
            ValueError: If ``parity`` is not recognised or no root exists.
        """
        if y_root is None:
            raise ValueError("No square root exists for the given x.")
        if parity == "even":
            y: int = y_root if y_root % 2 == 0 else (-y_root) % curve.p
        elif parity == "odd":
            y = y_root if y_root % 2 == 1 else (-y_root) % curve.p
        else:
            raise ValueError("parity must be 'even' or 'odd'.")
        return cls(x, y, curve)

    # -- Group operations ---------------------------------------------------

    def negate(self) -> "ECCPoint":
        """Return ``-P``: the additive inverse of this point.

        For a finite point ``P = (x, y)``, the inverse is ``(x, -y mod p)``;
        the inverse of the Point at Infinity is itself.
        """
        if self.is_infinity:
            return ECCPoint.infinity(self.curve)
        return ECCPoint(self.x, (-self.y) % self.curve.p, self.curve)

    def __neg__(self) -> "ECCPoint":
        """Unary minus operator alias for :meth:`negate`."""
        return self.negate()

    def add(self, other: "ECCPoint") -> "ECCPoint":
        """Add two points using the chord-and-tangent group law.

        Rules implemented (all arithmetic mod ``p``):

        * ``O + P = P + O = P``           (identity rule)
        * ``P + (-P) = O``                (inverse rule)
        * ``P + Q, P.x != Q.x``:          chord slope
              lambda = (yQ - yP) / (xQ - xP)
              xR = lambda^2 - xP - xQ ; yR = lambda*(xP - xR) - yP
        * ``2P, P.y != 0``:               tangent slope
              lambda = (3*xP^2 + a) / (2*yP)
              xR = lambda^2 - 2*xP     ; yR = lambda*(xP - xR) - yP
        * ``P + P`` with ``P.y == 0``:    result is ``O``

        Args:
            other: The right-hand addend.

        Returns:
            The sum as a new :class:`ECCPoint`.

        Raises:
            TypeError: If ``other`` lies on a different curve.
        """
        if other.curve is not self.curve:
            raise TypeError("Cannot add points defined on different curves.")

        p: int = self.curve.p
        a: int = self.curve.a

        # Identity rules.
        if self.is_infinity:
            return other
        if other.is_infinity:
            return self

        # Inverse rule: vertical line intersects only at infinity.
        if self.x == other.x and (self.y + other.y) % p == 0:
            return ECCPoint.infinity(self.curve)

        if self.x == other.x and self.y == other.y:
            # Point doubling (tangent case).
            if self.y == 0:
                return ECCPoint.infinity(self.curve)
            numerator: int = (3 * pow(self.x, 2, p) + a) % p
            denominator: int = (2 * self.y) % p
            lam: int = mod_mul(numerator, mod_inverse(denominator, p), p)
        else:
            # General point addition (chord case).
            numerator = (other.y - self.y) % p
            denominator = (other.x - self.x) % p
            lam = mod_mul(numerator, mod_inverse(denominator, p), p)

        x_r: int = (pow(lam, 2, p) - self.x - other.x) % p
        y_r: int = (lam * (self.x - x_r) - self.y) % p
        return ECCPoint(x_r, y_r, self.curve)

    def double(self) -> "ECCPoint":
        """Return ``2P`` via the tangent construction."""
        return self.add(self)

    def __add__(self, other: object) -> "ECCPoint":
        """Operator overload for ``+``."""
        if not isinstance(other, ECCPoint):
            return NotImplemented
        return self.add(other)

    def __sub__(self, other: object) -> "ECCPoint":
        """Operator overload for ``-``."""
        if not isinstance(other, ECCPoint):
            return NotImplemented
        return self.add(other.negate())

    def __eq__(self, other: object) -> bool:
        """Structural equality: same curve and same coordinates/infinity."""
        if not isinstance(other, ECCPoint):
            return NotImplemented
        if self.curve is not other.curve:
            return False
        if self.is_infinity or other.is_infinity:
            return self.is_infinity and other.is_infinity
        return self.x == other.x and self.y == other.y

    def __hash__(self) -> int:
        """Hash consistent with :meth:`__eq__`."""
        if self.is_infinity:
            return hash(("infinity", id(self.curve)))
        return hash((self.x, self.y, id(self.curve)))

    def scalar_multiply(self, k: int) -> "ECCPoint":
        """Compute ``k*P`` using the Double-and-Add ladder in O(log k).

        The binary expansion of ``k`` is scanned from the most significant
        bit downward; each step doubles the accumulator and conditionally
        adds the base point. This yields exactly ``floor(log2(k)) + 1``
        doubling rounds and at most that many additions.

        Negative scalars are handled via ``k*P = (-k)*(-P)``.

        Args:
            k: The integer scalar multiplier.

        Returns:
            The resulting :class:`ECCPoint` (possibly the Point at Infinity
            when ``k === 0`` or the point order divides ``k``).
        """
        if k < 0:
            return self.negate().scalar_multiply(-k)

        result: ECCPoint = ECCPoint.infinity(self.curve)
        addend: ECCPoint = self
        remaining: int = k

        while remaining > 0:
            if remaining & 1:
                result = result.add(addend)
            addend = addend.double()
            remaining >>= 1

        return result

    def __mul__(self, k: object) -> "ECCPoint":
        """Operator overload for ``P * k``."""
        if not isinstance(k, int):
            return NotImplemented
        return self.scalar_multiply(k)

    def __rmul__(self, k: object) -> "ECCPoint":
        """Operator overload for ``k * P``."""
        if not isinstance(k, int):
            return NotImplemented
        return self.scalar_multiply(k)

    # -- Introspection -------------------------------------------------------

    def __repr__(self) -> str:
        """Developer-facing representation."""
        if self.is_infinity:
            return f"ECCPoint(infinity, curve='{self.curve.name}')"
        return (
            f"ECCPoint(x={hex(self.x)}, y={hex(self.y)}, "
            f"curve='{self.curve.name}')"
        )

    def __str__(self) -> str:
        """Human-readable representation."""
        if self.is_infinity:
            return "Point at Infinity (O)"
        return f"({self.x}, {self.y})"

    def to_affine_tuple(self) -> Optional[Tuple[int, int]]:
        """Return ``(x, y)`` or ``None`` for the Point at Infinity."""
        if self.is_infinity:
            return None
        return (self.x, self.y)


# ---------------------------------------------------------------------------
# ECDH key exchange
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ECDHPublicShare:
    """One party's public contribution to an ECDH exchange.

    Attributes:
        public_point: The public key ``Q = d*G``.
        party_name: Label such as ``"Alice"`` or ``"Bob"``.
    """

    public_point: ECCPoint
    party_name: str


@dataclass(frozen=True)
class ECDHResult:
    """The outcome of a completed ECDH handshake.

    Attributes:
        shared_point: The identical raw shared point ``d_A * Q_B = d_B * Q_A``.
        shared_secret_hex: A hex-encoded symmetric key derived from the
            shared point via SHA-256 (used purely for demonstration).
        alice_public: Alice's public point.
        bob_public: Bob's public point.
        alice_private_bits: Bit length of Alice's private scalar.
        bob_private_bits: Bit length of Bob's private scalar.
    """

    shared_point: ECCPoint
    shared_secret_hex: str
    alice_public: ECCPoint
    bob_public: ECCPoint
    alice_private_bits: int
    bob_private_bits: int


class ECDHExchange:
    """Elliptic Curve Diffie-Hellman key exchange orchestrator.

    Protocol outline (with base point ``G`` of prime order ``n``):

    1. Alice selects secret ``d_A`` uniformly from ``[1, n-1]`` and publishes
       ``Q_A = d_A * G``.
    2. Bob selects secret ``d_B`` uniformly from ``[1, n-1]`` and publishes
       ``Q_B = d_B * G``.
    3. Alice computes ``S = d_A * Q_B``; Bob computes ``S = d_B * Q_A``.
    4. By associativity both equal ``d_A * d_B * G`` — the shared secret.

    Example:
        >>> exchange = ECDHExchange(CURVE_SECP256K1)
        >>> result = exchange.run_exchange()
        >>> len(result.shared_secret_hex)
        64
    """

    def __init__(
        self,
        curve: EllipticCurve,
        generator: Optional[ECCPoint] = None,
        rng_seed_private_a: Optional[int] = None,
        rng_seed_private_b: Optional[int] = None,
    ) -> None:
        """Initialise the exchange on a given curve.

        Args:
            curve: The elliptic curve domain parameters.
            generator: Override the base point; defaults to the curve's
                standard generator.
            rng_seed_private_a: Deterministic private scalar for Alice
                (testing only). ``None`` draws cryptographically secure
                randomness.
            rng_seed_private_b: Deterministic private scalar for Bob
                (testing only). ``None`` draws cryptographically secure
                randomness.

        Raises:
            ValueError: If the curve lacks a usable generator.
        """
        self.curve: EllipticCurve = curve
        if generator is None:
            generator = curve.generator
        if generator.is_infinity:
            raise ValueError("Generator cannot be the Point at Infinity.")
        self.generator: ECCPoint = generator

        order: Optional[int] = curve.n
        if order is not None and order > 2:
            upper_bound: int = order - 1
        else:
            # Toy curves without a declared order: bound by field size.
            upper_bound = curve.p - 1

        self._upper_bound: int = upper_bound
        self.alice_private: Optional[int] = rng_seed_private_a
        self.bob_private: Optional[int] = rng_seed_private_b
        self.alice_public: Optional[ECCPoint] = None
        self.bob_public: Optional[ECCPoint] = None
        self.shared_secret_point: Optional[ECCPoint] = None

    @property
    def order_bound(self) -> int:
        """The inclusive upper bound used for private scalar sampling."""
        return self._upper_bound

    def generate_private_scalar(self, forced: Optional[int] = None) -> int:
        """Sample a private scalar in ``[1, bound]``.

        Uses ``secrets.randbelow`` (CSPRNG backed by OS entropy) unless a
        deterministic value is supplied for reproducible testing.
        """
        if forced is not None:
            if not 1 <= forced <= self._upper_bound:
                raise ValueError(
                    f"Forced scalar must lie in [1, {self._upper_bound}]."
                )
            return forced
        return secrets.randbelow(self._upper_bound) + 1

    def compute_public_key(self, private_scalar: int) -> ECCPoint:
        """Derive the public key ``Q = d*G`` from a private scalar."""
        return self.generator.scalar_multiply(private_scalar)

    def derive_shared_secret(
        self,
        own_private: int,
        peer_public: ECDHPublicShare,
    ) -> ECCPoint:
        """Compute the raw shared point ``own_private * peer_public``.

        Args:
            own_private: This party's private scalar.
            peer_public: The counterparty's published public share.

        Returns:
            The shared :class:`ECCPoint`. Callers should reject an all-
            infinity result (a negligible-probability invalid-peer event).
        """
        if peer_public.public_point.curve is not self.curve:
            raise TypeError("Peer public key lies on a different curve.")
        return peer_public.public_point.scalar_multiply(own_private)

    @staticmethod
    def derive_symmetric_key(shared_point: ECCPoint) -> str:
        """Map the raw shared point to a 256-bit hex symmetric key.

        The canonical x-coordinate (or the literal string ``"infinity"``)
        is hashed with SHA-256. In production systems a KDF such as
        HKDF would be preferred; SHA-256 suffices for demonstration.
        """
        if shared_point.is_infinity:
            payload: bytes = b"infinity"
        else:
            byte_length: int = (shared_point.curve.p.bit_length() + 7) // 8
            payload = shared_point.x.to_bytes(byte_length, byteorder="big")
        digest: bytes = hashlib.sha256(payload).digest()
        return digest.hex()

    def run_exchange(
        self,
        alice_scalar: Optional[int] = None,
        bob_scalar: Optional[int] = None,
    ) -> ECDHResult:
        """Execute the full four-step ECDH handshake.

        Args:
            alice_scalar: Optional deterministic scalar for Alice.
            bob_scalar: Optional deterministic scalar for Bob.

        Returns:
            An :class:`ECDHResult` containing both public keys, the identical
            shared point computed independently by each party, and the
            derived hex symmetric key.

        Raises:
            RuntimeError: If the two parties derive divergent shared points
                (indicative of an internal invariant violation).
        """
        d_a: int = self.generate_private_scalar(alice_scalar)
        d_b: int = self.generate_private_scalar(bob_scalar)

        # Step 1 & 2: publish public points.
        q_a: ECCPoint = self.compute_public_key(d_a)
        q_b: ECCPoint = self.compute_public_key(d_b)

        share_a: ECDHPublicShare = ECDHPublicShare(q_a, "Alice")
        share_b: ECDHPublicShare = ECDHPublicShare(q_b, "Bob")

        # Step 3: each party mixes its private scalar with the peer's
        # public point.
        s_alice: ECCPoint = self.derive_shared_secret(d_a, share_b)
        s_bob: ECCPoint = self.derive_shared_secret(d_b, share_a)

        if s_alice != s_bob:
            raise RuntimeError(
                "ECDH invariant violated: parties derived different "
                "shared points."
            )

        self.alice_private = d_a
        self.bob_private = d_b
        self.alice_public = q_a
        self.bob_public = q_b
        self.shared_secret_point = s_alice

        return ECDHResult(
            shared_point=s_alice,
            shared_secret_hex=ECDHExchange.derive_symmetric_key(s_alice),
            alice_public=q_a,
            bob_public=q_b,
            alice_private_bits=d_a.bit_length(),
            bob_private_bits=d_b.bit_length(),
        )


# ---------------------------------------------------------------------------
# Discrete-logarithm utilities (for pedagogy / toy fields only)
# ---------------------------------------------------------------------------


def brute_force_dlog(target: ECCPoint, generator: ECCPoint,
                     max_k: int) -> Optional[int]:
    """Recover ``k`` with ``k*G == target`` by exhaustive search up to max_k.

    WARNING: Exponential-time routine intended exclusively for tiny toy
    fields to illustrate why ECDLP is hard at cryptographic scale.
    """
    if target.is_infinity:
        return 0
    accumulator: ECCPoint = ECCPoint.infinity(generator.curve)
    for candidate in range(1, max_k + 1):
        accumulator = accumulator.add(generator)
        if accumulator == target:
            return candidate
    return None


def enumerate_curve_points(curve: EllipticCurve) -> List[Tuple[int, int]]:
    """Enumerate every affine point on ``curve`` by scanning all x in F_p.

    For each x, the right-hand side is evaluated and a square root is
    extracted with Tonelli-Shanks; both mirror roots yield distinct points.

    Complexity: O(p) field evaluations — suitable only for visualisation
    fields (p on the order of hundreds).
    """
    points: List[Tuple[int, int]] = []
    for x in range(curve.p):
        rhs: int = curve.evaluate_rhs(x)
        root: Optional[int] = tonelli_shanks_sqrt(rhs, curve.p)
        if root is None:
            continue
        points.append((x, root))
        mirrored: int = (-root) % curve.p
        if mirrored != root:
            points.append((x, mirrored))
    points.sort()
    return points


# ---------------------------------------------------------------------------
# Module self-test hook
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    demo_curve: EllipticCurve = CURVE_TOY_97
    demo_gen: ECCPoint = ECCPoint(3, 6, demo_curve)
    print(f"Curve: {demo_curve.name}: "
          f"y^2 = x^3 + {demo_curve.a}x + {demo_curve.b} (mod {demo_curve.p})")
    print(f"Base point G = {demo_gen}")
    print(f"2G  = {demo_gen.double()}")
    print(f"5G  = {demo_gen.scalar_multiply(5)}")
    print(f"10G = {demo_gen.scalar_multiply(10)}")

    secp_demo: ECDHExchange = ECDHExchange(CURVE_SECP256K1)
    secp_result: ECDHResult = secp_demo.run_exchange(alice_scalar=7, bob_scalar=13)
    print("")
    print("secp256k1 ECDH demonstration:")
    print(f"  Alice public Q_A = {secp_result.alice_public!r}")
    print(f"  Bob   public Q_B = {secp_result.bob_public!r}")
    print(f"  Shared secret (SHA-256 derived): "
          f"{secp_result.shared_secret_hex[:32]}...")