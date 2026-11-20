//! A second seeded family, kept apart from `fresh`.
//!
//! `fresh` is frozen: its generator and its random stream are the ones the
//! known-answer baseline was recorded with, so a seed replays the same edges.
//! Its surface family draws every control height of the patch at random, and
//! some of those controls sit inside a knot cluster as narrow as 2^-20. That
//! makes slopes of millions of millimetres per unit parameter, which is a
//! stress case for rounding and not a surface any builder produces.
//!
//! This family has the same hidden bump on a smooth patch. It draws from its
//! own stream, so adding or changing it cannot move a `fresh` edge.

use super::exact::*;
use super::*;

/// SplitMix64, as in `fresh` (Steele, Lea and Flood, OOPSLA 2014; Vigna's
/// public-domain `splitmix64.c`). A separate copy, because the frozen
/// module's methods are private to it.
struct Stream(u64);

impl Stream {
    /// The complement of the run's seed starts a stream apart from `fresh`.
    fn new(seed: u64) -> Self {
        Self(!seed)
    }
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    /// A whole number in `lo..=hi`.
    fn between(&mut self, lo: i64, hi: i64) -> i64 {
        lo + (self.next() % (hi - lo + 1) as u64) as i64
    }
    /// A dyadic fraction strictly inside `(0, 1)` with `bits` binary places.
    fn share(&mut self, bits: i64) -> Rational {
        qi(self.between(1, (1 << bits) - 1)).mul_pow2(-bits)
    }
}

/// `least..=most` distinct sorted dyadic knots strictly inside `(0, 1)`.
fn random_knots(random: &mut Stream, least: i64, most: i64) -> Vec<Rational> {
    let count = random.between(least, most) as usize;
    let mut knots = std::collections::BTreeSet::new();
    while knots.len() < count {
        knots.insert(random.share(10));
    }
    knots.into_iter().collect()
}

/// Five distinct knots of random spacing in a window 2^-3 to 2^-20 wide that
/// lies strictly inside one span of the interior `knots`: the support of one
/// cubic basis function.
fn cluster(random: &mut Stream, knots: &[Rational]) -> Vec<Rational> {
    loop {
        let width = pow2(-random.between(3, 20));
        let start = random.share(30);
        let end = &start + &width;
        let clear = end < Rational::one() && knots.iter().all(|knot| knot < &start || knot > &end);
        if !clear {
            continue;
        }
        let mut shares = std::collections::BTreeSet::new();
        while shares.len() < 3 {
            shares.insert(random.share(7));
        }
        let mut cluster = vec![start.clone()];
        cluster.extend(shares.iter().map(|share| &start + &(&width * share)));
        cluster.push(end);
        return cluster;
    }
}

/// A smooth random spline patch and its own iso-curve as the 3D copy. One
/// control of the patch is then lifted, so the bump lives on the surface: its
/// size is the lift times the v basis value at the iso-parameter times the u
/// basis function along the edge.
fn surface_bump(random: &mut Stream) -> (Curve, CaseSurface, FaceCopy, Truth) {
    let degree_v = random.between(1, 3) as usize;
    let mut interior = random_knots(random, 1, 6);
    let cluster = cluster(random, &interior);
    interior.extend(cluster.iter().cloned());
    interior.sort();
    let knots_u = clamped(3, &interior, &Rational::zero(), &Rational::one());
    let knots_v = clamped(
        degree_v,
        &random_knots(random, 0, 3),
        &Rational::zero(),
        &Rational::one(),
    );
    let (us, vs) = (greville(&knots_u, 3), greville(&knots_v, degree_v));
    let reach = [10, 1_000, 1_000_000][random.between(0, 2) as usize];
    let origin: [Rational; 3] = std::array::from_fn(|_| qi(random.between(-reach, reach)));
    // The heights follow a random quadratic in the Greville abscissae, so the
    // patch is smooth and only the hidden bump is narrow.
    let bend: [Rational; 6] = std::array::from_fn(|_| ratio(random.between(-64, 64), 8));
    let mut net: Vec<Vec<[Rational; 3]>> = us
        .iter()
        .map(|u| {
            vs.iter()
                .map(|v| {
                    let terms = [Rational::one(), u.clone(), v.clone(), u * v, u * u, v * v];
                    let height = bend
                        .iter()
                        .zip(&terms)
                        .fold(Rational::zero(), |sum, (c, term)| sum + c * term);
                    [
                        &origin[0] + &(qi(100) * u.clone()),
                        &origin[1] + &(qi(20) * v.clone()),
                        &origin[2] + &height,
                    ]
                })
                .collect()
        })
        .collect();
    let v = loop {
        let v = random.share(12);
        if !knots_v.contains(&v) {
            break v;
        }
    };
    let shares: Vec<_> = (0..vs.len())
        .map(|j| basis(&knots_v, degree_v, j, &v))
        .collect();
    let row = net
        .iter()
        .map(|column| {
            std::array::from_fn(|axis| {
                column
                    .iter()
                    .zip(&shares)
                    .fold(Rational::zero(), |sum, (p, share)| sum + share * &p[axis])
            })
        })
        .collect();
    let i = knots_u
        .iter()
        .position(|knot| knot == &cluster[0])
        .expect("cluster");
    let j = (0..vs.len())
        .max_by(|a, b| shares[*a].cmp(&shares[*b]))
        .expect("columns");
    // A bump size between 1e-7 and 1e-2 mm.
    let lift =
        qi(random.between(10, 99)) * qi(10).pow(-(random.between(4, 8) as i32)).expect("power");
    net[i][j][2] = &net[i][j][2] + &lift;
    let mut along = vec![Rational::zero(); us.len()];
    along[i] = &lift * &shares[j];
    (
        Curve::BSpline(spline(3, knots_u.clone(), row, None)),
        CaseSurface::Carrier(carrier(
            [3, degree_v],
            [knots_u.clone(), knots_v],
            net.into_iter().flatten().collect(),
        )),
        FaceCopy::Native(straight(
            [Rational::zero(), v.clone()],
            [Rational::one(), v],
        )),
        spline_peak(&knots_u, 3, &along),
    )
}

/// `count` smooth-surface edges for a seed.
pub(super) fn cases(seed: u64, count: usize) -> Vec<EdgeCase> {
    let mut random = Stream::new(seed);
    (0..count)
        .map(|round| {
            let (curve, surface, face_copy, truth) = surface_bump(&mut random);
            EdgeCase {
                name: format!("smooth-surface-bump-{round}"),
                source: "seeded: hidden bump on one control of a smooth surface",
                geometry: EdgeGeometry {
                    curve,
                    parameters: [Rational::zero(), Rational::one()],
                    surface,
                    face_copy,
                },
                truth,
                domain_exit: None,
                zero: None,
                zero_ceiling: Rational::zero(),
            }
        })
        .collect()
}
