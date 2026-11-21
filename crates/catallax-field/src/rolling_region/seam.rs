//! A ball rolling over an earlier mitre of two equal rounds (4.2a step 5).
//!
//! An earlier fillet rounds the two edges where a plane `T` meets planes `A`
//! and `B` at one radius `R`, and the rounds meet in a mitre. A later straight
//! run between `A` and `B`, with a ball of radius `r < R` of the rounds'
//! convexity, reaches both rounds at once at its end `p`. RECIPE 4.5.1: "The
//! ball rolls between the two rounds to the mitre's apex, with its centre on
//! the ellipse where their offsets meet, and sweeps one exact rational band."
//!
//! Frame. `n_a` and `n_b` are the unit normals of `A` and `B` toward the ball
//! and `d` the run's unit direction past `p`, pairwise orthogonal. The rounds'
//! axes run along `n_b` and `n_a` and meet at `q = p + R (n_a + n_b)`. With
//! `x - q = a n_a + b n_b + h d`, a point lies `sqrt(a^2 + h^2)` from the
//! first axis and `sqrt(b^2 + h^2)` from the second. The ball lies inside both
//! cylinders, so its centre lies `k = R - r` from both axes: `a^2 = b^2`. The
//! branch through the run's last centre `p + r (n_a + n_b)` is `a = b`, the
//! ellipse
//!
//! ```text
//! c(t) = q - k cos t (n_a + n_b) + k sin t d,      t in [0, pi/2],
//! ```
//!
//! from that centre to the apex `q + k d`, whose ball touches `T` where both
//! rounds' rails meet. It is step 3's crossing of two equal cylinders on
//! meeting axes (their offsets meet in the planes bisecting the axes), cut to
//! one quarter of one ellipse. The contact rays `m_1 = -cos t n_a + sin t d`
//! and `m_2 = -cos t n_b + sin t d` are radial from each axis and square to
//! `c'(t) = k (sin t (n_a + n_b) + cos t d)`, so the face is the arc between
//! them in the normal plane. With `tau = tan(t/2) = (u + 1) / 2`, every
//! station of a rational frame is rational.
//!
//! The fan. A point `c(t) + s_1 m_1 + s_2 m_2` with `s_1, s_2 >= 0` has
//! `a = -(k + s_1) cos t`, `b = -(k + s_2) cos t` and
//! `h = (k + s_1 + s_2) sin t`. In the open octant `a, b < 0 < h`, eliminating
//! `s_1` and `s_2` leaves `g(t) = -(a + b) / cos t - h / sin t = k`, and `g`
//! rises from minus to plus infinity on `(0, pi/2)`, so exactly one station
//! holds the point's section. By that rise, `s_1 >= 0` (`cos t <= -a / k`) is
//! `a^2 (b^2 + h^2) >= k^2 b^2`, and `s_2 >= 0` is
//! `b^2 (a^2 + h^2) >= k^2 a^2`. These two polynomial walls, `a <= 0`,
//! `b <= 0` and the ball side of both supports (`a^2 + h^2 <= R^2`,
//! `b^2 + h^2 <= R^2`) bound the region, all exact. `h >= 0` is the junction
//! with the run's straight piece, whose section is the one at `t = 0`. Off its
//! centres the fan lies outside both offsets, `a^2 + h^2 > k^2` and
//! `b^2 + h^2 > k^2`: two bands, which keep out the segment `a = b = 0`,
//! `h < k`, where both polynomial walls, `a` and `b` vanish.
//!
//! No closed form gives the station (a quartic in `tau`), so membership near
//! the face and distance take the whole-interval panel cover of step 1, as
//! the saddle does. The field derives the path from the run's two planes, the
//! rounds' radius and the ball's; nothing here reads a B-rep fillet module or
//! a fillet definition.

use super::elementary::{balls, local, point_balls, position, positive, reach, vector};
use super::path::{Cover, Frame, Path};
use super::{BITS, interval, q};
use crate::geometry::{enclose, sign};
use catallax_numbers::{Ball, Decider, Rational, Real, Refusal, Sign};
use catallax_shapes::{Isometry, Point3, Vector3};

/// The ball path over a mitre of two equal rounds: a quarter of the ellipse
/// where their offsets meet.
#[derive(Clone, Debug)]
pub struct Seam {
    // Where the rounds' axes meet.
    meet: Point3,
    // The unit normals of the run's planes toward the ball; the first
    // round's axis runs along the second normal and the reverse.
    normals: [Vector3; 2],
    // The run's unit direction past its end, toward the plane the rounds
    // roll on.
    out: Vector3,
    // The rounds' radius.
    round: Real,
    radius: Real,
    // The box `a, b in [-R, 0]`, `h in [0, R]` in world coordinates.
    bounds: [[Rational; 2]; 3],
}

/// Interval data of a [`Seam`].
#[derive(Clone, Debug)]
pub struct SeamCover {
    meet: [Ball; 3],
    normals: [[Ball; 3]; 2],
    // `n_a + n_b`.
    sum: [Ball; 3],
    out: [Ball; 3],
    // `k = R - r`.
    offset: Ball,
    radius: Ball,
}

/// `cos t` and `sin t` at `tan(t/2) = (u + 1) / 2`, exactly.
fn half_angle(u: &Rational) -> (Rational, Rational) {
    let tau = &(u + &q(1)) / &q(2);
    let squared = &tau * &tau;
    let weight = &q(1) + &squared;
    (&(&q(1) - &squared) / &weight, &(&q(2) * &tau) / &weight)
}

impl SeamCover {
    /// Enclosures of `cos t` and `sin t` over a panel of `u`. On `[-1, 1]`
    /// the angle runs over `[0, pi/2]`, where the cosine falls and the sine
    /// rises, so the panel's ends bound both.
    fn trig(u: &Ball) -> (Ball, Ball) {
        let clamp = |x: Rational| x.max(q(-1)).min(q(1));
        let (cos_lo, sin_lo) = half_angle(&clamp(u.lo()));
        let (cos_hi, sin_hi) = half_angle(&clamp(u.hi()));
        (interval(&cos_hi, &cos_lo), interval(&sin_lo, &sin_hi))
    }
    fn centre_at(&self, cos: &Ball, sin: &Ball) -> [Ball; 3] {
        let across = &self.offset * cos;
        let up = &self.offset * sin;
        std::array::from_fn(|i| {
            &(&self.meet[i] - &(&across * &self.sum[i])) + &(&up * &self.out[i])
        })
    }
}

impl Cover for SeamCover {
    fn radius(&self) -> &Ball {
        &self.radius
    }
    fn centre(&self, u: &Ball) -> Result<[Ball; 3], Refusal> {
        let (cos, sin) = Self::trig(u);
        Ok(self.centre_at(&cos, &sin))
    }
    fn frame(&self, u: &Ball) -> Result<Frame, Refusal> {
        let (cos, sin) = Self::trig(u);
        let centre = self.centre_at(&cos, &sin);
        let tangent = std::array::from_fn(|i| &(&sin * &self.sum[i]) + &(&cos * &self.out[i]));
        let feet = [0, 1].map(|j| {
            std::array::from_fn(|i| {
                let ray = &(&sin * &self.out[i]) - &(&cos * &self.normals[j][i]);
                &centre[i] + &(&self.radius * &ray)
            })
        });
        Ok(Frame {
            centre,
            tangent,
            feet,
            side: None,
        })
    }
}

impl Seam {
    /// The path over the mitre whose rounds' axes meet at `meet`, from the
    /// run's plane normals toward the ball, its direction past its end, the
    /// rounds' radius and the ball's. The three directions must be pairwise
    /// orthogonal units and the ball narrower than the rounds.
    pub fn new(
        meet: Point3,
        normals: [Vector3; 2],
        out: Vector3,
        round: Real,
        radius: Real,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        dec.charge_work(1)?;
        if !positive(&radius, dec)? || !positive(&(&round - &radius), dec)? {
            return Err(reach(
                "the ball must be narrower than the rounds it rolls between",
            ));
        }
        let directions = [&normals[0], &normals[1], &out];
        for (i, v) in directions.iter().enumerate() {
            if sign(&(&v.norm_squared() - &Real::one()), dec)? != Sign::Zero {
                return Err(reach("a mitre's frame needs unit directions"));
            }
            for w in &directions[i + 1..] {
                if sign(&v.dot(w), dec)? != Sign::Zero {
                    return Err(reach("a mitre's frame needs orthogonal directions"));
                }
            }
        }
        let minus = -&round;
        let mut corners = Vec::with_capacity(8);
        for corner in 0..8 {
            let mut v = meet.clone();
            for (bit, (direction, length)) in
                [(&normals[0], &minus), (&normals[1], &minus), (&out, &round)]
                    .into_iter()
                    .enumerate()
            {
                if (corner >> bit) & 1 == 1 {
                    v = v.plus(&direction.scaled(length));
                }
            }
            corners.push(point_balls(&v, dec)?);
        }
        let mut bounds: [[Rational; 2]; 3] = std::array::from_fn(|i| {
            let x = &corners[0][i];
            [x.lo(), x.hi()]
        });
        for corner in &corners {
            for (axis, x) in corner.iter().enumerate() {
                bounds[axis][0] = bounds[axis][0].clone().min(x.lo());
                bounds[axis][1] = bounds[axis][1].clone().max(x.hi());
            }
        }
        Ok(Self {
            meet,
            normals,
            out,
            round,
            radius,
            bounds,
        })
    }

    /// Ball radius in mm.
    pub fn radius(&self) -> &Real {
        &self.radius
    }

    fn offset(&self) -> Real {
        &self.round - &self.radius
    }

    /// The point's coordinates `a`, `b`, `h` in the mitre's frame.
    fn coordinates(&self, local: &[Real; 3]) -> [Real; 3] {
        let v = vector(local).minus(&position(&self.meet));
        [
            v.dot(&self.normals[0]),
            v.dot(&self.normals[1]),
            v.dot(&self.out),
        ]
    }

    /// Applies an exact rigid motion, reflections included; the box moves by
    /// its corners.
    pub fn transform(&mut self, pose: &Isometry) {
        self.meet = pose.apply_point(&self.meet);
        for n in &mut self.normals {
            *n = pose.apply_vector(n);
        }
        self.out = pose.apply_vector(&self.out);
        self.bounds = super::moved_box(pose, &self.bounds);
    }

    /// An axis-aligned box holding the region.
    pub fn hull(&self, _dec: &mut Decider) -> Result<[Ball; 3], Refusal> {
        Ok(self.bounds.each_ref().map(|[lo, hi]| interval(lo, hi)))
    }

    /// Boxes that together hold the region: its one box.
    pub fn boxes(&self, dec: &mut Decider) -> Result<Vec<[Ball; 3]>, Refusal> {
        Ok(vec![self.hull(dec)?])
    }

    /// The path in binary64, for display values and starting guesses.
    fn rough(&self) -> Rough {
        let f = |v: &Vector3| [v.x.to_f64(), v.y.to_f64(), v.z.to_f64()];
        Rough {
            meet: [
                self.meet.x.to_f64(),
                self.meet.y.to_f64(),
                self.meet.z.to_f64(),
            ],
            normals: [f(&self.normals[0]), f(&self.normals[1])],
            out: f(&self.out),
            offset: self.offset().to_f64(),
        }
    }

    /// Tentative binary64 value for display: negative inside the region. The
    /// square root is the caller's. The nearest centre comes from 512 samples
    /// of the path, so near the face it is only a display value.
    pub(crate) fn fast(&self, p: [f64; 3], sqrt: fn(f64) -> f64) -> f64 {
        let rough = self.rough();
        let r = self.radius.to_f64();
        let big = self.round.to_f64();
        let k = rough.offset;
        let mut nearest = f64::INFINITY;
        const SAMPLES: usize = 512;
        for i in 0..SAMPLES {
            let u = -1.0 + 2.0 * (i as f64 + 0.5) / SAMPLES as f64;
            let c = rough.centre(u);
            nearest = nearest.min(sqrt(dot3(sub3(p, c), sub3(p, c))));
        }
        let v = sub3(p, rough.meet);
        let [a, b, h] = [
            dot3(v, rough.normals[0]),
            dot3(v, rough.normals[1]),
            dot3(v, rough.out),
        ];
        let (first, second) = (a * a + h * h, b * b + h * h);
        let s1 = big - sqrt(first);
        let s2 = big - sqrt(second);
        let fan = a <= 0.0
            && b <= 0.0
            && h >= 0.0
            && b * b * first >= k * k * a * a
            && a * a * second >= k * k * b * b;
        if fan && s1 >= 0.0 && s2 >= 0.0 && nearest > r {
            (-s1).max(-s2).max(r - nearest)
        } else {
            (nearest - r).abs().max(f64::MIN_POSITIVE)
        }
    }
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// A [`Seam`] in binary64.
struct Rough {
    meet: [f64; 3],
    normals: [[f64; 3]; 2],
    out: [f64; 3],
    offset: f64,
}
impl Rough {
    /// The centre at `u`, mapped as the cover maps it.
    fn centre(&self, u: f64) -> [f64; 3] {
        let tau = (u + 1.0) / 2.0;
        let weight = 1.0 + tau * tau;
        let (cos, sin) = ((1.0 - tau * tau) / weight, 2.0 * tau / weight);
        std::array::from_fn(|i| {
            self.meet[i] - self.offset * cos * (self.normals[0][i] + self.normals[1][i])
                + self.offset * sin * self.out[i]
        })
    }

    /// The squared distance's slope in `u` at `u`, up to a positive factor:
    /// `(c - p) . c'`, with `c'` the tangent `sin t (n_a + n_b) + cos t d`
    /// over `k`, and `t` rising with `u`.
    fn slope(&self, u: f64, p: [f64; 3]) -> f64 {
        let tau = (u + 1.0) / 2.0;
        let weight = 1.0 + tau * tau;
        let (cos, sin) = ((1.0 - tau * tau) / weight, 2.0 * tau / weight);
        let tangent: [f64; 3] = std::array::from_fn(|i| {
            sin * (self.normals[0][i] + self.normals[1][i]) + cos * self.out[i]
        });
        dot3(sub3(self.centre(u), p), tangent)
    }

    /// The parameter of the centre nearest `p`: the best of 1024 samples,
    /// polished by a golden-section search on its neighbourhood, as for the
    /// saddle. A golden-section search settles only near the square root of
    /// the unit roundoff, because the squared distance is flat at its least,
    /// so where the slope changes sign across the final bracket, bisection
    /// on the slope's sign finishes it.
    fn nearest(&self, p: [f64; 3]) -> f64 {
        const SAMPLES: usize = 1024;
        let gap = |u: f64| {
            let c = self.centre(u);
            dot3(sub3(p, c), sub3(p, c))
        };
        let step = 2.0 / SAMPLES as f64;
        let mut best = (f64::INFINITY, 0.0);
        for i in 0..SAMPLES {
            let u = -1.0 + step * (i as f64 + 0.5);
            let g = gap(u);
            if g < best.0 {
                best = (g, u);
            }
        }
        let ratio = (5f64.sqrt() - 1.0) / 2.0;
        let (mut lo, mut hi) = ((best.1 - step).max(-1.0), (best.1 + step).min(1.0));
        for _ in 0..80 {
            let a = hi - ratio * (hi - lo);
            let b = lo + ratio * (hi - lo);
            if gap(a) < gap(b) {
                hi = b;
            } else {
                lo = a;
            }
        }
        let golden = (lo + hi) / 2.0;
        // Near an end where the distance only grows inward, the end itself
        // is the least.
        if golden - step <= -1.0 && self.slope(-1.0, p) >= 0.0 {
            return -1.0;
        }
        if golden + step >= 1.0 && self.slope(1.0, p) <= 0.0 {
            return 1.0;
        }
        let (mut lo, mut hi) = ((lo - step).max(-1.0), (hi + step).min(1.0));
        if !(self.slope(lo, p) < 0.0 && self.slope(hi, p) > 0.0) {
            // The slope's sign is lost in roundoff: the golden-section
            // station stands.
            return golden;
        }
        for _ in 0..64 {
            let mid = (lo + hi) / 2.0;
            if self.slope(mid, p) < 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        (lo + hi) / 2.0
    }
}

impl Path for Seam {
    type Cover = SeamCover;
    fn local(&self, point: &Point3) -> [Real; 3] {
        local(point)
    }
    fn walls(&self, local: &[Real; 3], _: &mut Decider) -> Result<Vec<Real>, Refusal> {
        let [a, b, h] = self.coordinates(local);
        let (aa, bb, hh) = (a.square(), b.square(), h.square());
        let first = &aa + &hh;
        let second = &bb + &hh;
        let rr = self.round.square();
        let kk = self.offset().square();
        Ok(vec![
            &rr - &first,
            &rr - &second,
            &(&aa * &second) - &(&kk * &bb),
            &(&bb * &first) - &(&kk * &aa),
            -&a,
            -&b,
        ])
    }
    fn walls_are_faces(&self) -> bool {
        true
    }
    fn bands(&self, local: &[Real; 3], _: &mut Decider) -> Result<Vec<Real>, Refusal> {
        let [a, b, h] = self.coordinates(local);
        let hh = h.square();
        let kk = self.offset().square();
        Ok(vec![&(&a.square() + &hh) - &kk, &(&b.square() + &hh) - &kk])
    }
    fn junctions(&self, local: &[Real; 3], _: &mut Decider) -> Result<Vec<Real>, Refusal> {
        let [_, _, h] = self.coordinates(local);
        Ok(vec![h])
    }
    fn cover(&self, dec: &mut Decider) -> Result<SeamCover, Refusal> {
        let normals = [balls(&self.normals[0], dec)?, balls(&self.normals[1], dec)?];
        let sum = std::array::from_fn(|i| &normals[0][i] + &normals[1][i]);
        Ok(SeamCover {
            meet: point_balls(&self.meet, dec)?,
            normals,
            sum,
            out: balls(&self.out, dec)?,
            offset: enclose(&self.offset(), BITS, dec)?,
            radius: enclose(&self.radius, BITS, dec)?,
        })
    }
    fn guess(&self, local: &[Real; 3]) -> Option<Rational> {
        let p = local.each_ref().map(|x| x.to_f64());
        Rational::from_f64_exact(self.rough().nearest(p))
    }
}
