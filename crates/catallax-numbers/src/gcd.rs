//! The greatest common divisor of two integers.
//!
//! A rational is kept in lowest terms, so every sum and product of two
//! rationals ends in one or two of these. The work meter prices a rational
//! operation by the limbs of its operands, multiplied or added: the limb
//! operations of the schoolbook product or sum (`decider.rs`, "Logical
//! arithmetic costs"). The reduction has to cost about that much too, or time
//! and counted work part ways on small rationals.
//!
//! `num-integer`'s algorithm for a `BigUint` is Stein's binary one. It makes
//! one pass over both operands for about every two bits, so 32 n passes for n
//! limbs, each a shift and a subtraction, against the n passes of a schoolbook
//! product. It does that even when one operand is 1.
//!
//! Three cases here:
//!
//! 1. An operand that is 0 or 1 needs no arithmetic.
//! 2. Two operands of at most 128 bits use the binary algorithm on machine
//!    words (Stein 1967; Knuth, The Art of Computer Programming, volume 2,
//!    section 4.5.2, Algorithm B).
//! 3. Longer operands use Lehmer's algorithm (Lehmer 1938; Knuth, the same
//!    section, Algorithm L). The leading bits of both operands, held in
//!    machine words, give the quotients of several steps of Euclid's algorithm
//!    at once, and one pass over the operands applies them.
//!
//! The divisor of two integers is unique, so every caller gets the value it
//! got before.

use num_bigint::{BigInt, BigUint};
use num_traits::{One, ToPrimitive, Zero};

/// Non-negative greatest common divisor; `gcd(0, 0)` is 0.
pub(crate) fn gcd(a: &BigInt, b: &BigInt) -> BigInt {
    BigInt::from(of_magnitudes(a.magnitude(), b.magnitude()))
}

fn of_magnitudes(a: &BigUint, b: &BigUint) -> BigUint {
    if a.is_zero() {
        return b.clone();
    }
    if b.is_zero() {
        return a.clone();
    }
    if a.is_one() || b.is_one() {
        return BigUint::one();
    }
    if let (Some(a), Some(b)) = (a.to_u128(), b.to_u128()) {
        return BigUint::from(binary(a, b));
    }
    lehmer(a, b)
}

// Algorithm B on machine words. `of(0, b)` is `b`.
macro_rules! algorithm_b {
    ($(#[$doc:meta])* $name:ident, $word:ty) => {
        $(#[$doc])*
        pub(crate) fn $name(mut a: $word, mut b: $word) -> $word {
            if a == 0 || b == 0 {
                return a | b;
            }
            let twos = (a | b).trailing_zeros();
            a >>= a.trailing_zeros();
            loop {
                b >>= b.trailing_zeros();
                if a > b {
                    std::mem::swap(&mut a, &mut b);
                }
                // Only the shared twos divide 1 and another number.
                if a == 1 {
                    return 1 << twos;
                }
                b -= a;
                if b == 0 {
                    return a << twos;
                }
            }
        }
    };
}
algorithm_b!(
    /// The divisor of two double words.
    binary,
    u128
);
algorithm_b!(
    /// The divisor of two single words.
    word,
    u64
);

// The bits Algorithm L reads at the head of the longer operand. Knuth's
// bounds in step L2 keep the four cofactors of the simulated steps at or
// below 2^LEAD in magnitude, so 62 leaves each one inside a 64-bit word and
// every intermediate inside a signed 128-bit one.
const LEAD: u64 = 62;

/// Algorithm L for two non-zero operands, at least one longer than 128 bits.
fn lehmer(x: &BigUint, y: &BigUint) -> BigUint {
    // `u >= v` here and after every round: each round replaces the two by a
    // later pair of consecutive remainders of Euclid's algorithm on them.
    let (mut u, mut v) = if x >= y {
        (x.clone(), y.clone())
    } else {
        (y.clone(), x.clone())
    };
    while v.bits() > 128 {
        // Step L1: the leading parts, cut at the same bit of both.
        let shift = u.bits() - LEAD;
        let (mut head_u, mut head_v) = (lead(&u, shift), lead(&v, shift));
        let (mut a, mut b, mut c, mut d) = (1i128, 0i128, 0i128, 1i128);
        // Steps L2 and L3: a quotient is certain while both bounds of the
        // true ratio give the same one.
        while head_v + c != 0 && head_v + d != 0 {
            let q = (head_u + a) / (head_v + c);
            if q != (head_u + b) / (head_v + d) {
                break;
            }
            (a, c) = (c, a - q * c);
            (b, d) = (d, b - q * d);
            (head_u, head_v) = (head_v, head_u - q * head_v);
        }
        // Step L4.
        if b == 0 {
            // The leading parts decided no step. One full division.
            let r = &u % &v;
            u = std::mem::replace(&mut v, r);
        } else {
            let next = combine(a, &u, b, &v);
            v = combine(c, &u, d, &v);
            u = next;
        }
    }
    if v.is_zero() {
        return u;
    }
    let r = &u % &v;
    let words = |x: &BigUint| x.to_u128().expect("at most 128 bits");
    BigUint::from(binary(words(&v), words(&r)))
}

/// Bits `shift` to `shift + LEAD` of `x`, which has none above them.
fn lead(x: &BigUint, shift: u64) -> i128 {
    let skip = usize::try_from(shift / 64).expect("a limb index");
    let mut limbs = x.iter_u64_digits().skip(skip);
    let low = u128::from(limbs.next().unwrap_or(0));
    let high = u128::from(limbs.next().unwrap_or(0));
    i128::try_from(((high << 64) | low) >> (shift % 64)).expect("no more than LEAD bits")
}

/// `a u + b v`, a remainder of Euclid's algorithm and so not negative. The
/// two cofactors never have the same sign.
fn combine(a: i128, u: &BigUint, b: i128, v: &BigUint) -> BigUint {
    let word = |k: i128| u64::try_from(k.unsigned_abs()).expect("a cofactor fits a word");
    if b <= 0 {
        u * word(a) - v * word(b)
    } else {
        v * word(b) - u * word(a)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_integer::Integer as _;

    // A fixed sequence (Marsaglia's xorshift64), so a failure repeats.
    struct Sequence(u64);
    impl Sequence {
        fn word(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
        fn below(&mut self, n: u64) -> u64 {
            self.word() % n
        }
        fn integer(&mut self, bits: u64) -> BigUint {
            let mut value = BigUint::zero();
            for _ in 0..bits.div_ceil(64) {
                value = (value << 64u32) | BigUint::from(self.word());
            }
            value >> (bits.div_ceil(64) * 64 - bits)
        }
    }

    fn agrees(a: &BigUint, b: &BigUint) {
        let expected = a.gcd(b);
        assert_eq!(of_magnitudes(a, b), expected, "gcd({a}, {b})");
        assert_eq!(of_magnitudes(b, a), expected, "gcd({b}, {a})");
    }

    #[test]
    fn the_divisor_is_stein_s_on_operands_of_every_length() {
        let mut sequence = Sequence(0x9E37_79B9_7F4A_7C15);
        for _ in 0..4000 {
            let bits = [
                1 + sequence.below(700),
                1 + sequence.below(700),
                sequence.below(300),
            ];
            let shared = sequence.integer(bits[2]) + 1u32;
            let a = sequence.integer(bits[0]) * &shared;
            let b = sequence.integer(bits[1]) * &shared;
            agrees(&a, &b);
        }
    }

    #[test]
    fn the_divisor_of_operands_that_meet_a_case_boundary() {
        let one = BigUint::one();
        let long = (&one << 300u32) + 12_345u32;
        let edge = &one << 128u32;
        let cases = [
            (BigUint::zero(), BigUint::zero()),
            (BigUint::zero(), long.clone()),
            (one.clone(), long.clone()),
            (long.clone(), long.clone()),
            (long.clone(), &long * 7u32),
            (long.clone(), &long + 1u32),
            (&edge - 1u32, edge.clone()),
            (edge.clone(), &edge + 1u32),
            (&edge * 3u32, &edge * 9u32),
            (&one << 500u32, &one << 200u32),
            (&long * &edge, BigUint::from(6u32)),
            // The same leading bits: no step is certain from them.
            ((&one << 400u32) + 5u32, (&one << 400u32) + 3u32),
        ];
        for (a, b) in &cases {
            agrees(a, b);
        }
    }

    #[test]
    fn the_divisor_of_neighbours_of_the_slowest_sequence() {
        // Consecutive Fibonacci numbers make every quotient 1, the most steps
        // Euclid's algorithm takes for operands of their length (Lame 1844).
        let (mut a, mut b) = (BigUint::one(), BigUint::one());
        for _ in 0..1200 {
            let next = &a + &b;
            a = std::mem::replace(&mut b, next);
            agrees(&a, &b);
            agrees(&(&a * 6u32), &(&b * 10u32));
        }
    }

    #[test]
    fn the_divisor_of_machine_words_is_stein_s() {
        let mut sequence = Sequence(0xD1B5_4A32_D192_ED03);
        for _ in 0..4000 {
            let bits = sequence.below(20);
            let shared = 1 + sequence.below(1 << bits);
            let mut multiple = || {
                let cut = sequence.below(64);
                (sequence.word() >> cut) / shared * shared
            };
            let (a, b) = (multiple(), multiple());
            assert_eq!(word(a, b), a.gcd(&b), "gcd({a}, {b})");
            let (wide_a, wide_b) = (u128::from(a) << 40, u128::from(b) * 3);
            assert_eq!(binary(wide_a, wide_b), wide_a.gcd(&wide_b));
        }
        for (a, b) in [
            (0, 0),
            (0, 9),
            (9, 0),
            (1, 1 << 40),
            (1 << 40, 1 << 12),
            (3, 3),
        ] {
            assert_eq!(word(a, b), a.gcd(&b), "gcd({a}, {b})");
        }
    }

    #[test]
    fn the_sign_of_an_operand_does_not_reach_the_divisor() {
        let a = BigInt::from(-84);
        let b = BigInt::from(36);
        assert_eq!(gcd(&a, &b), BigInt::from(12));
        assert_eq!(gcd(&b, &a), BigInt::from(12));
        assert_eq!(gcd(&-&b, &a), BigInt::from(12));
        assert_eq!(gcd(&a, &BigInt::zero()), BigInt::from(84));
    }
}
