"""
test_ecc.py
===========

Native ``unittest`` test suite for the library-free ECC engine.

Coverage
--------
1. Modular arithmetic primitives: extended GCD / Bezout identity, modular
   inversion, and round-trip properties.
2. Algebraic group axioms on elliptic curve points:
   - Identity:  P + O = O + P = P
   - Inverse:   P + (-P) = O  and  -(-P) = P
   - Closure:   every sum of two curve points remains a valid curve point
   - Commutativity and associativity (supporting axioms)
3. Scalar multiplication correctness against naive repeated addition.
4. ECDH key exchange: 100 randomized iterations asserting that Alice and Bob
   independently derive identical shared secrets, plus determinism checks
   with forced scalars.

Run with:
    python -m unittest test_ecc -v
"""

from __future__ import annotations

import random
import unittest

from ecc_engine import (
    CURVE_SECP256K1,
    CURVE_TOY_97,
    ECCPoint,
    ECDHExchange,
    ECDHResult,
    EllipticCurve,
    extended_gcd,
    mod_inverse,
    tonelli_shanks_sqrt,
)


class ModularArithmeticTests(unittest.TestCase):
    """Verify the arbitrary-precision modular arithmetic primitives."""

    def test_extended_gcd_bezout_identity(self) -> None:
        """a*x + b*y must equal gcd(a, b) for assorted operand pairs."""
        cases: list[tuple[int, int]] = [
            (240, 46),
            (46, 240),
            (17, 5),
            (0, 7),
            (7, 0),
            (270, 192),
            (2**64 + 3, 2**32 + 9),
            (-24, 60),
        ]
        for a, b in cases:
            with self.subTest(a=a, b=b):
                g, x, y = extended_gcd(a, b)
                self.assertEqual(a * x + b * y, g)
                self.assertGreaterEqual(g, 0)
                self.assertEqual(g, extended_gcd(abs(a), abs(b))[0])

    def test_mod_inverse_round_trip(self) -> None:
        """n * mod_inverse(n, p) mod p must equal 1."""
        prime: int = CURVE_TOY_97.p
        for n in range(1, prime):
            with self.subTest(n=n):
                inv: int = mod_inverse(n, prime)
                self.assertEqual((n * inv) % prime, 1)

    def test_mod_inverse_large_prime(self) -> None:
        """Inversion must succeed on the 256-bit secp256k1 prime."""
        p: int = CURVE_SECP256K1.p
        n: int = 0x123456789ABCDEF0123456789ABCDEF0123456789ABCDEF0
        inv: int = mod_inverse(n, p)
        self.assertEqual((n % p) * inv % p, 1)

    def test_mod_inverse_non_invertible_raises(self) -> None:
        """Zero has no inverse; requesting one must raise ValueError."""
        with self.assertRaises(ValueError):
            mod_inverse(0, 97)

    def test_tonelli_shanks_square_roots(self) -> None:
        """Every returned root must square back to the original value."""
        prime: int = 97  # 97 === 1 (mod 4): exercises general Tonelli-Shanks.
        for x in range(prime):
            square: int = pow(x, 2, prime)
            root = tonelli_shanks_sqrt(square, prime)
            with self.subTest(x=x):
                self.assertIsNotNone(root)
                assert root is not None
                self.assertEqual(pow(root, 2, prime), square)


class GroupAxiomTests(unittest.TestCase):
    """Verify the algebraic group axioms of the elliptic curve point set."""

    def setUp(self) -> None:
        """Enumerate all affine points of the toy curve for reuse."""
        self.curve: EllipticCurve = CURVE_TOY_97
        from ecc_engine import enumerate_curve_points

        coords: list[tuple[int, int]] = enumerate_curve_points(self.curve)
        self.points: list[ECCPoint] = [
            ECCPoint(x, y, self.curve) for x, y in coords
        ]
        self.infinity: ECCPoint = ECCPoint.infinity(self.curve)

    def test_identity_property(self) -> None:
        """P + O = O + P = P for every point P, including O itself."""
        for point in self.points + [self.infinity]:
            with self.subTest(point=str(point)):
                self.assertEqual(point.add(self.infinity), point)
                self.assertEqual(self.infinity.add(point), point)

    def test_inverse_property(self) -> None:
        """P + (-P) = O and -(-P) = P for every point P."""
        for point in self.points:
            with self.subTest(point=str(point)):
                negated: ECCPoint = point.negate()
                self.assertTrue(
                    negated.is_infinity
                    or self.curve.contains_point(negated.x, negated.y)
                )
                self.assertEqual(point.add(negated), self.infinity)
                self.assertEqual(negated.negate(), point)
        # The inverse of infinity is infinity.
        self.assertEqual(self.infinity.negate(), self.infinity)

    def test_closure_property(self) -> None:
        """The sum of any two curve points is again a valid curve point."""
        enumerated: set[tuple[int, int]] = {
            (p.x, p.y) for p in self.points
        }
        for left in self.points:
            for right in self.points:
                total: ECCPoint = left.add(right)
                if total.is_infinity:
                    continue
                key: tuple[int, int] = (total.x, total.y)
                self.assertIn(
                    key,
                    enumerated,
                    f"Sum {left} + {right} = {total} escaped the point set.",
                )

    def test_commutativity(self) -> None:
        """P + Q = Q + P for sampled point pairs."""
        rng: random.Random = random.Random(20240601)
        for _ in range(200):
            left: ECCPoint = rng.choice(self.points)
            right: ECCPoint = rng.choice(self.points)
            with self.subTest(left=str(left), right=str(right)):
                self.assertEqual(left.add(right), right.add(left))

    def test_associativity(self) -> None:
        """(P + Q) + R = P + (Q + R) for sampled point triples."""
        rng: random.Random = random.Random(20240602)
        for _ in range(100):
            first: ECCPoint = rng.choice(self.points)
            second: ECCPoint = rng.choice(self.points)
            third: ECCPoint = rng.choice(self.points)
            with self.subTest(
                first=str(first), second=str(second), third=str(third)
            ):
                self.assertEqual(
                    first.add(second).add(third),
                    first.add(second.add(third)),
                )


class ScalarMultiplicationTests(unittest.TestCase):
    """Verify double-and-add scalar multiplication against naive addition."""

    def setUp(self) -> None:
        """Fix the toy-curve generator used across these tests."""
        self.curve: EllipticCurve = CURVE_TOY_97
        self.generator: ECCPoint = ECCPoint(3, 6, self.curve)

    def test_scalar_multiply_matches_repeated_addition(self) -> None:
        """k*G via double-and-add equals k sequential additions of G."""
        accumulator: ECCPoint = ECCPoint.infinity(self.curve)
        for k in range(0, 40):
            with self.subTest(k=k):
                self.assertEqual(
                    self.generator.scalar_multiply(k), accumulator
                )
            accumulator = accumulator.add(self.generator)

    def test_scalar_multiply_zero_is_infinity(self) -> None:
        """0*G must be the Point at Infinity."""
        result: ECCPoint = self.generator.scalar_multiply(0)
        self.assertTrue(result.is_infinity)

    def test_scalar_multiply_negative_scalar(self) -> None:
        """(-k)*G must equal -(k*G)."""
        positive: ECCPoint = self.generator.scalar_multiply(13)
        negative: ECCPoint = self.generator.scalar_multiply(-13)
        self.assertEqual(negative, positive.negate())

    def test_secp256k1_generator_order(self) -> None:
        """n*G must be the Point at Infinity on secp256k1."""
        n: int = CURVE_SECP256K1.n
        if n is not None:
            result: ECCPoint = CURVE_SECP256K1.generator.scalar_multiply(n)
            self.assertTrue(result.is_infinity)


class ECDHExchangeTests(unittest.TestCase):
    """Randomized end-to-end verification of the ECDH protocol."""

    ITERATIONS: int = 100

    def test_randomized_ecdh_shared_secret_agreement(self) -> None:
        """Across 100 randomized runs Alice and Bob must agree exactly.

        Each iteration draws fresh private scalars from the CSPRNG-backed
        generator, executes the full handshake, and asserts that both
        parties derive the same shared point and symmetric key.
        """
        exchange: ECDHExchange = ECDHExchange(CURVE_SECP256K1)
        for iteration in range(ECDHExchangeTests.ITERATIONS):
            with self.subTest(iteration=iteration):
                result: ECDHResult = exchange.run_exchange()

                # Independent recomputation of each party's view.
                alice_view: ECCPoint = result.bob_public.scalar_multiply(
                    exchange.alice_private  # type: ignore[arg-type]
                )
                bob_view: ECCPoint = result.alice_public.scalar_multiply(
                    exchange.bob_private  # type: ignore[arg-type]
                )

                self.assertFalse(result.shared_point.is_infinity)
                self.assertEqual(alice_view, result.shared_point)
                self.assertEqual(bob_view, result.shared_point)
                self.assertEqual(alice_view, bob_view)
                self.assertEqual(len(result.shared_secret_hex), 64)

    def test_deterministic_scalars_produce_stable_secret(self) -> None:
        """Forced scalars must yield reproducible shared secrets."""
        first: ECDHResult = ECDHExchange(CURVE_SECP256K1).run_exchange(
            alice_scalar=42, bob_scalar=1337
        )
        second: ECDHResult = ECDHExchange(CURVE_SECP256K1).run_exchange(
            alice_scalar=42, bob_scalar=1337
        )
        self.assertEqual(first.shared_point, second.shared_point)
        self.assertEqual(first.shared_secret_hex, second.shared_secret_hex)

    def test_swapped_scalars_change_secret(self) -> None:
        """Different scalar pairs must (overwhelmingly) differ."""
        base: ECDHResult = ECDHExchange(CURVE_SECP256K1).run_exchange(
            alice_scalar=42, bob_scalar=1337
        )
        swapped: ECDHResult = ECDHExchange(CURVE_SECP256K1).run_exchange(
            alice_scalar=43, bob_scalar=1337
        )
        self.assertNotEqual(base.shared_point, swapped.shared_point)

    def test_forced_scalar_bounds_enforced(self) -> None:
        """Out-of-range forced scalars must be rejected."""
        exchange: ECDHExchange = ECDHExchange(CURVE_SECP256K1)
        with self.assertRaises(ValueError):
            exchange.generate_private_scalar(forced=0)
        with self.assertRaises(ValueError):
            exchange.generate_private_scalar(
                forced=CURVE_SECP256K1.n  # type: ignore[arg-type]
            )

    def test_toy_curve_handshake_agreement(self) -> None:
        """ECDH must also agree on the small visualisation field."""
        curve: EllipticCurve = CURVE_TOY_97
        generator: ECCPoint = ECCPoint(3, 6, curve)
        exchange: ECDHExchange = ECDHExchange(curve, generator=generator)
        for iteration in range(25):
            with self.subTest(iteration=iteration):
                result: ECDHResult = exchange.run_exchange()
                alice_view: ECCPoint = result.bob_public.scalar_multiply(
                    exchange.alice_private  # type: ignore[arg-type]
                )
                self.assertEqual(alice_view, result.shared_point)


if __name__ == "__main__":
    unittest.main(verbosity=2)