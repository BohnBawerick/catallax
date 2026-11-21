//! The filtered predicate of a decision with a tolerant participant.
//!
//! Checkpoint 2 decisions, section 4.5: "Any decision with a tolerant
//! participant is 'distance against the sum of tolerances' through a
//! filtered predicate: f64, then interval, then exact (the CGAL and Shewchuk
//! pattern; Catallax already uses it in mesh and field)."
//!
//! A decision here is the sign of a polynomial of rational operands: a
//! squared distance less a squared reach, where the reach is a size plus the
//! tolerances of the participants. The sign is read in three stages, and a
//! stage is asked only where the one before it could not tell:
//!
//! 1. Binary64, with a bound of its own rounding error kept beside every
//!    value (the running error of Higham, "Accuracy and Stability of
//!    Numerical Algorithms", second edition, section 3.3). The sign is taken
//!    where the value is larger than its bound. This is the first stage of
//!    Shewchuk, "Adaptive Precision Floating-Point Arithmetic and Fast Robust
//!    Geometric Predicates" (Discrete and Computational Geometry 18, 1997),
//!    with the bound computed with the value instead of ahead of it.
//! 2. Outward-rounded binary64 intervals ([`Interval`]), the dynamic filter
//!    of Bronnimann, Burnikel and Pion, "Interval arithmetic yields efficient
//!    dynamic filters for computational geometry" (Discrete Applied
//!    Mathematics 109, 2001), which CGAL's filtered kernel uses.
//! 3. Rational arithmetic, which gives the sign of the polynomial itself.
//!
//! So the sign returned is always the exact one, zero included. The stages
//! change what it costs and never what it is. OpenCASCADE decides the same
//! kind of question, a distance against a sum of tolerances, in binary64
//! alone (`IntTools_Context::ComputeVE` and the vertex and edge interferences
//! of `BOPAlgo_PaveFiller`); the stages after the first are what the
//! decisions add to it.
//!
//! Work: one unit is charged before each of the first two stages, and the
//! third reserves every rational operation by the limb rule before it runs,
//! as the decider does.
//!
//! Only operands that are rational have a third stage. A curve or a surface
//! that is not a line is not decided here: the decisions rule out the exact
//! stage for it ("the exact fallback is impractical on splines"), and its
//! caller bounds it by enclosures and decides the points it is left with.

use catallax_numbers::{Decider, Interval, Rational, Refusal, Sign, StrictSign};

/// How many signs each stage settled.
#[derive(Clone, Copy, Default, PartialEq, Eq, Debug)]
pub struct Stages {
    /// Settled in binary64 with its error bound.
    pub float: u64,
    /// Settled in outward-rounded intervals.
    pub interval: u64,
    /// Settled in rational arithmetic.
    pub exact: u64,
}

impl Stages {
    /// Every sign taken.
    pub fn total(&self) -> u64 {
        self.float + self.interval + self.exact
    }
}

/// The number type of one stage.
trait Arithmetic {
    type Number: Clone;
    fn lift(&mut self, value: &Rational) -> Self::Number;
    fn sum(&mut self, a: &Self::Number, b: &Self::Number) -> Result<Self::Number, Refusal>;
    fn difference(&mut self, a: &Self::Number, b: &Self::Number) -> Result<Self::Number, Refusal>;
    fn product(&mut self, a: &Self::Number, b: &Self::Number) -> Result<Self::Number, Refusal>;
}

/// A binary64 value and a bound of its distance from the true value.
#[derive(Clone, Copy)]
struct Rounded {
    value: f64,
    error: f64,
}

/// The first stage.
///
/// One operation rounds its result by at most half a unit in the last place,
/// `2^-53` of the result. Each operation adds `2^-52` of its result to the
/// bound, twice that: the other half covers the rounding of the bound's own
/// sums and products, which are of the same relative size (Higham, section
/// 3.3). The smallest positive binary64 is added for a result that rounds
/// into the subnormal range. A rational that is read is taken as rounded
/// the same way.
struct Floats;

impl Floats {
    fn rounded(value: f64, carried: f64) -> Rounded {
        Rounded {
            value,
            error: carried + value.abs() * f64::EPSILON + f64::MIN_POSITIVE,
        }
    }
}

impl Arithmetic for Floats {
    type Number = Rounded;
    fn lift(&mut self, value: &Rational) -> Rounded {
        Self::rounded(value.to_f64_nearest(), 0.0)
    }
    fn sum(&mut self, a: &Rounded, b: &Rounded) -> Result<Rounded, Refusal> {
        Ok(Self::rounded(a.value + b.value, a.error + b.error))
    }
    fn difference(&mut self, a: &Rounded, b: &Rounded) -> Result<Rounded, Refusal> {
        Ok(Self::rounded(a.value - b.value, a.error + b.error))
    }
    fn product(&mut self, a: &Rounded, b: &Rounded) -> Result<Rounded, Refusal> {
        Ok(Self::rounded(
            a.value * b.value,
            a.value.abs() * b.error + b.value.abs() * a.error + a.error * b.error,
        ))
    }
}

/// The second stage.
struct Intervals;

impl Arithmetic for Intervals {
    type Number = Interval;
    fn lift(&mut self, value: &Rational) -> Interval {
        Interval::from_rational(value)
    }
    fn sum(&mut self, a: &Interval, b: &Interval) -> Result<Interval, Refusal> {
        Ok(*a + *b)
    }
    fn difference(&mut self, a: &Interval, b: &Interval) -> Result<Interval, Refusal> {
        Ok(*a - *b)
    }
    fn product(&mut self, a: &Interval, b: &Interval) -> Result<Interval, Refusal> {
        Ok(*a * *b)
    }
}

/// The third stage. Each operation is reserved by the limb rule before it
/// runs, and counts itself when it does.
struct Exact<'a> {
    dec: &'a mut Decider,
}

impl Arithmetic for Exact<'_> {
    type Number = Rational;
    fn lift(&mut self, value: &Rational) -> Rational {
        value.clone()
    }
    fn sum(&mut self, a: &Rational, b: &Rational) -> Result<Rational, Refusal> {
        self.dec.reserve_exact_arithmetic(a, b, false)?;
        Ok(a + b)
    }
    fn difference(&mut self, a: &Rational, b: &Rational) -> Result<Rational, Refusal> {
        self.dec.reserve_exact_arithmetic(a, b, false)?;
        Ok(a - b)
    }
    fn product(&mut self, a: &Rational, b: &Rational) -> Result<Rational, Refusal> {
        self.dec.reserve_exact_arithmetic(a, b, true)?;
        Ok(a * b)
    }
}

type Triple<N> = [N; 3];

fn lifted<A: Arithmetic>(a: &mut A, v: &[Rational; 3]) -> Triple<A::Number> {
    [a.lift(&v[0]), a.lift(&v[1]), a.lift(&v[2])]
}

fn minus<A: Arithmetic>(
    a: &mut A,
    x: &Triple<A::Number>,
    y: &Triple<A::Number>,
) -> Result<Triple<A::Number>, Refusal> {
    Ok([
        a.difference(&x[0], &y[0])?,
        a.difference(&x[1], &y[1])?,
        a.difference(&x[2], &y[2])?,
    ])
}

fn dot<A: Arithmetic>(
    a: &mut A,
    x: &Triple<A::Number>,
    y: &Triple<A::Number>,
) -> Result<A::Number, Refusal> {
    let first = a.product(&x[0], &y[0])?;
    let second = a.product(&x[1], &y[1])?;
    let third = a.product(&x[2], &y[2])?;
    let sum = a.sum(&first, &second)?;
    a.sum(&sum, &third)
}

fn cross<A: Arithmetic>(
    a: &mut A,
    x: &Triple<A::Number>,
    y: &Triple<A::Number>,
) -> Result<Triple<A::Number>, Refusal> {
    let mut out = Vec::with_capacity(3);
    for (i, j) in [(1, 2), (2, 0), (0, 1)] {
        let left = a.product(&x[i], &y[j])?;
        let right = a.product(&x[j], &y[i])?;
        out.push(a.difference(&left, &right)?);
    }
    let [p, q, r]: [A::Number; 3] = out.try_into().ok().expect("three components");
    Ok([p, q, r])
}

/// A polynomial whose sign is a decision.
enum Polynomial<'a> {
    /// `|(p - o) x d|^2 - reach^2 (d . d)`: the squared distance of `p`
    /// from the line less the squared reach, times `d . d`.
    PointLine {
        line: &'a Line,
        point: &'a [Rational; 3],
        reach: &'a Rational,
    },
    /// `w . w` with `w = (b - a) x d`: zero where the segment runs along
    /// the line's direction.
    Across {
        line: &'a Line,
        a: &'a [Rational; 3],
        b: &'a [Rational; 3],
    },
    /// `u . w + from_b (w . w)` with `u = (a - o) x d`: the side of the
    /// segment's end `a`, or of its end `b`, on which the foot of the line
    /// lies.
    Foot {
        line: &'a Line,
        a: &'a [Rational; 3],
        b: &'a [Rational; 3],
        from_b: bool,
    },
    /// `(u . u)(w . w) - (u . w)^2 - reach^2 (d . d)(w . w)`: the squared
    /// distance between the two lines less the squared reach, times
    /// `(d . d)(w . w)`.
    LineLine {
        line: &'a Line,
        a: &'a [Rational; 3],
        b: &'a [Rational; 3],
        reach: &'a Rational,
    },
    /// `p.x q.y - p.y q.x`: the turn from `p` to `q` about the origin.
    Turn {
        p: &'a [Rational; 2],
        q: &'a [Rational; 2],
    },
    /// `|p - q|^2 - reach^2`.
    Gap {
        p: &'a [Rational; 3],
        q: &'a [Rational; 3],
        reach: &'a Rational,
    },
    /// `(p - a) . (b - a)`, less `|b - a|^2` from `b`: the side of the
    /// segment's end on which the foot of `p` lies.
    FootOfPoint {
        p: &'a [Rational; 3],
        a: &'a [Rational; 3],
        b: &'a [Rational; 3],
        from_b: bool,
    },
    /// `|u|^2 |w|^2 - (u . w)^2 - reach^2 |w|^2` with `u = p - a` and
    /// `w = b - a`: the squared distance of `p` from the line through the
    /// segment less the squared reach, times `|w|^2`.
    PointFromLine {
        p: &'a [Rational; 3],
        a: &'a [Rational; 3],
        b: &'a [Rational; 3],
        reach: &'a Rational,
    },
    /// A quantity of the two segments `a b` and `c d`, with `d1 = b - a`,
    /// `d2 = d - c` and `r = a - c` (Ericson, "Real-Time Collision
    /// Detection", 2005, section 5.1.9).
    Segments {
        ends: [&'a [Rational; 3]; 4],
        what: OfSegments<'a>,
    },
}

/// What is asked of two segments.
#[derive(Clone, Copy)]
enum OfSegments<'a> {
    /// `|d1|^2 |d2|^2 - (d1 . d2)^2`: zero where they are parallel.
    Skew,
    /// The numerator of the common perpendicular's foot on the first or the
    /// second segment, less the denominator [`Self::Skew`] from its far end.
    Foot { second: bool, from_end: bool },
    /// `(r . n)^2 - reach^2 |n|^2` with `n = d1 x d2`: the squared distance
    /// between the two lines less the squared reach, times `|n|^2`.
    Apart(&'a Rational),
}

impl Polynomial<'_> {
    fn value<A: Arithmetic>(&self, a: &mut A) -> Result<A::Number, Refusal> {
        // `u = (a - o) x d` and `w = (b - a) x d`.
        fn frame<A: Arithmetic>(
            a: &mut A,
            line: &Line,
            from: &[Rational; 3],
            to: &[Rational; 3],
        ) -> Result<[Triple<A::Number>; 3], Refusal> {
            let d = lifted(a, &line.direction);
            let origin = lifted(a, &line.origin);
            let from = lifted(a, from);
            let to = lifted(a, to);
            let offset = minus(a, &from, &origin)?;
            let along = minus(a, &to, &from)?;
            Ok([cross(a, &offset, &d)?, cross(a, &along, &d)?, d])
        }
        match self {
            Self::PointLine { line, point, reach } => {
                let d = lifted(a, &line.direction);
                let origin = lifted(a, &line.origin);
                let point = lifted(a, point);
                let offset = minus(a, &point, &origin)?;
                let w = cross(a, &offset, &d)?;
                let distance = dot(a, &w, &w)?;
                let reach = a.lift(reach);
                let reach = a.product(&reach, &reach)?;
                let length = dot(a, &d, &d)?;
                let reach = a.product(&reach, &length)?;
                a.difference(&distance, &reach)
            }
            Self::Across {
                line,
                a: from,
                b: to,
            } => {
                let [_, w, _] = frame(a, line, from, to)?;
                dot(a, &w, &w)
            }
            Self::Foot {
                line,
                a: from,
                b: to,
                from_b,
            } => {
                let [u, w, _] = frame(a, line, from, to)?;
                let foot = dot(a, &u, &w)?;
                if *from_b {
                    let whole = dot(a, &w, &w)?;
                    a.sum(&foot, &whole)
                } else {
                    Ok(foot)
                }
            }
            Self::LineLine {
                line,
                a: from,
                b: to,
                reach,
            } => {
                let [u, w, d] = frame(a, line, from, to)?;
                let uu = dot(a, &u, &u)?;
                let ww = dot(a, &w, &w)?;
                let uw = dot(a, &u, &w)?;
                let apart = a.product(&uu, &ww)?;
                let foot = a.product(&uw, &uw)?;
                let apart = a.difference(&apart, &foot)?;
                let reach = a.lift(reach);
                let reach = a.product(&reach, &reach)?;
                let length = dot(a, &d, &d)?;
                let reach = a.product(&reach, &length)?;
                let reach = a.product(&reach, &ww)?;
                a.difference(&apart, &reach)
            }
            Self::Turn { p, q } => {
                let [px, py] = [a.lift(&p[0]), a.lift(&p[1])];
                let [qx, qy] = [a.lift(&q[0]), a.lift(&q[1])];
                let left = a.product(&px, &qy)?;
                let right = a.product(&py, &qx)?;
                a.difference(&left, &right)
            }
            Self::Gap { p, q, reach } => {
                let (p, q) = (lifted(a, p), lifted(a, q));
                let between = minus(a, &p, &q)?;
                let distance = dot(a, &between, &between)?;
                let reach = a.lift(reach);
                let reach = a.product(&reach, &reach)?;
                a.difference(&distance, &reach)
            }
            Self::FootOfPoint {
                p,
                a: from,
                b: to,
                from_b,
            } => {
                let (p, from, to) = (lifted(a, p), lifted(a, from), lifted(a, to));
                let u = minus(a, &p, &from)?;
                let w = minus(a, &to, &from)?;
                let foot = dot(a, &u, &w)?;
                if *from_b {
                    let whole = dot(a, &w, &w)?;
                    a.difference(&foot, &whole)
                } else {
                    Ok(foot)
                }
            }
            Self::PointFromLine {
                p,
                a: from,
                b: to,
                reach,
            } => {
                let (p, from, to) = (lifted(a, p), lifted(a, from), lifted(a, to));
                let u = minus(a, &p, &from)?;
                let w = minus(a, &to, &from)?;
                let uu = dot(a, &u, &u)?;
                let ww = dot(a, &w, &w)?;
                let uw = dot(a, &u, &w)?;
                let apart = a.product(&uu, &ww)?;
                let foot = a.product(&uw, &uw)?;
                let apart = a.difference(&apart, &foot)?;
                let reach = a.lift(reach);
                let reach = a.product(&reach, &reach)?;
                let reach = a.product(&reach, &ww)?;
                a.difference(&apart, &reach)
            }
            Self::Segments { ends, what } => {
                let [from, to, other, end] = ends.map(|point| lifted(a, point));
                let d1 = minus(a, &to, &from)?;
                let d2 = minus(a, &end, &other)?;
                let r = minus(a, &from, &other)?;
                let skew = |a: &mut A| -> Result<A::Number, Refusal> {
                    let first = dot(a, &d1, &d1)?;
                    let second = dot(a, &d2, &d2)?;
                    let both = dot(a, &d1, &d2)?;
                    let lengths = a.product(&first, &second)?;
                    let square = a.product(&both, &both)?;
                    a.difference(&lengths, &square)
                };
                match what {
                    OfSegments::Skew => skew(a),
                    OfSegments::Foot { second, from_end } => {
                        let both = dot(a, &d1, &d2)?;
                        let on_first = dot(a, &d1, &r)?;
                        let on_second = dot(a, &d2, &r)?;
                        let (left, right) = if *second {
                            // (d1 . d1)(d2 . r) - (d1 . d2)(d1 . r)
                            let first = dot(a, &d1, &d1)?;
                            (a.product(&first, &on_second)?, a.product(&both, &on_first)?)
                        } else {
                            // (d1 . d2)(d2 . r) - (d1 . r)(d2 . d2)
                            let second = dot(a, &d2, &d2)?;
                            (
                                a.product(&both, &on_second)?,
                                a.product(&on_first, &second)?,
                            )
                        };
                        let foot = a.difference(&left, &right)?;
                        if *from_end {
                            let whole = skew(a)?;
                            a.difference(&foot, &whole)
                        } else {
                            Ok(foot)
                        }
                    }
                    OfSegments::Apart(reach) => {
                        let n = cross(a, &d1, &d2)?;
                        let across = dot(a, &r, &n)?;
                        let apart = a.product(&across, &across)?;
                        let reach = a.lift(reach);
                        let reach = a.product(&reach, &reach)?;
                        let length = dot(a, &n, &n)?;
                        let reach = a.product(&reach, &length)?;
                        a.difference(&apart, &reach)
                    }
                }
            }
        }
    }

    /// The sign, by the first stage that can tell it.
    fn sign(&self, stages: &mut Stages, dec: &mut Decider) -> Result<Sign, Refusal> {
        dec.charge_work(1)?;
        let Rounded { value, error } = self.value(&mut Floats)?;
        if value.is_finite() && error.is_finite() && value.abs() > error {
            stages.float += 1;
            return Ok(if value > 0.0 {
                Sign::Positive
            } else {
                Sign::Negative
            });
        }
        dec.charge_work(1)?;
        if let Some(sign) = self.value(&mut Intervals)?.sign() {
            stages.interval += 1;
            return Ok(match sign {
                StrictSign::Positive => Sign::Positive,
                StrictSign::Negative => Sign::Negative,
            });
        }
        let exact = self.value(&mut Exact { dec })?;
        stages.exact += 1;
        Ok(exact.sign())
    }
}

/// A line through `origin` along `direction`, of any length but zero.
#[derive(Clone, Debug)]
pub struct Line {
    /// A point of the line.
    pub origin: [Rational; 3],
    /// Its direction.
    pub direction: [Rational; 3],
}

impl Line {
    /// The sign of the distance from `point` to this line less `reach`:
    /// positive where the point is farther than the reach, zero where it is
    /// at the reach, negative inside it. `reach` must not be negative.
    pub fn point_beyond(
        &self,
        point: &[Rational; 3],
        reach: &Rational,
        stages: &mut Stages,
        dec: &mut Decider,
    ) -> Result<Sign, Refusal> {
        Polynomial::PointLine {
            line: self,
            point,
            reach,
        }
        .sign(stages, dec)
    }

    /// The sign of the least distance from the segment `a b` to this line
    /// less `reach`: positive where every point of the segment is farther
    /// than the reach.
    ///
    /// The nearest point of the segment is an end, or the foot of the common
    /// perpendicular where that lies between the ends. Which it is, is two
    /// signs of the same kind, and each goes through the stages.
    pub fn segment_beyond(
        &self,
        a: &[Rational; 3],
        b: &[Rational; 3],
        reach: &Rational,
        stages: &mut Stages,
        dec: &mut Decider,
    ) -> Result<Sign, Refusal> {
        let at = |point, stages: &mut Stages, dec: &mut Decider| {
            Polynomial::PointLine {
                line: self,
                point,
                reach,
            }
            .sign(stages, dec)
        };
        if (Polynomial::Across { line: self, a, b }).sign(stages, dec)? == Sign::Zero {
            // Along the line's direction: every point is as far as `a`.
            return at(a, stages, dec);
        }
        let foot = |from_b, stages: &mut Stages, dec: &mut Decider| {
            Polynomial::Foot {
                line: self,
                a,
                b,
                from_b,
            }
            .sign(stages, dec)
        };
        if foot(false, stages, dec)? != Sign::Negative {
            return at(a, stages, dec);
        }
        if foot(true, stages, dec)? != Sign::Positive {
            return at(b, stages, dec);
        }
        Polynomial::LineLine {
            line: self,
            a,
            b,
            reach,
        }
        .sign(stages, dec)
    }
}

/// The sign of the distance between two points less `reach`.
pub fn points_beyond(
    p: &[Rational; 3],
    q: &[Rational; 3],
    reach: &Rational,
    stages: &mut Stages,
    dec: &mut Decider,
) -> Result<Sign, Refusal> {
    Polynomial::Gap { p, q, reach }.sign(stages, dec)
}

/// The sign of the least distance from `p` to the segment `a b` less
/// `reach`. The nearest point of the segment is an end, or the foot of `p`
/// where that lies between the ends, and which it is, is two signs.
pub fn point_segment_beyond(
    p: &[Rational; 3],
    a: &[Rational; 3],
    b: &[Rational; 3],
    reach: &Rational,
    stages: &mut Stages,
    dec: &mut Decider,
) -> Result<Sign, Refusal> {
    let foot = |from_b, stages: &mut Stages, dec: &mut Decider| {
        Polynomial::FootOfPoint { p, a, b, from_b }.sign(stages, dec)
    };
    if foot(false, stages, dec)? != Sign::Positive {
        return points_beyond(p, a, reach, stages, dec);
    }
    if foot(true, stages, dec)? != Sign::Negative {
        return points_beyond(p, b, reach, stages, dec);
    }
    Polynomial::PointFromLine { p, a, b, reach }.sign(stages, dec)
}

/// The sign of the least distance between the segments `a b` and `c d` less
/// `reach`.
///
/// The least distance is between the feet of the common perpendicular where
/// both lie inside their segments, and from an end of one segment to the
/// other segment otherwise (Ericson, "Real-Time Collision Detection", 2005,
/// section 5.1.9). The four ends are always read, and the perpendicular
/// where its feet are inside, and the least of the signs is returned.
pub fn segments_beyond(
    [a, b]: [&[Rational; 3]; 2],
    [c, d]: [&[Rational; 3]; 2],
    reach: &Rational,
    stages: &mut Stages,
    dec: &mut Decider,
) -> Result<Sign, Refusal> {
    let order = |sign: Sign| match sign {
        Sign::Negative => 0,
        Sign::Zero => 1,
        Sign::Positive => 2,
    };
    let mut least = Sign::Positive;
    for (point, from, to) in [(a, c, d), (b, c, d), (c, a, b), (d, a, b)] {
        let sign = point_segment_beyond(point, from, to, reach, stages, dec)?;
        if order(sign) < order(least) {
            least = sign;
        }
    }
    let of = |what, stages: &mut Stages, dec: &mut Decider| {
        Polynomial::Segments {
            ends: [a, b, c, d],
            what,
        }
        .sign(stages, dec)
    };
    let inside = |second, stages: &mut Stages, dec: &mut Decider| -> Result<bool, Refusal> {
        Ok(of(
            OfSegments::Foot {
                second,
                from_end: false,
            },
            stages,
            dec,
        )? == Sign::Positive
            && of(
                OfSegments::Foot {
                    second,
                    from_end: true,
                },
                stages,
                dec,
            )? == Sign::Negative)
    };
    if least != Sign::Negative
        && of(OfSegments::Skew, stages, dec)? == Sign::Positive
        && inside(false, stages, dec)?
        && inside(true, stages, dec)?
    {
        let sign = of(OfSegments::Apart(reach), stages, dec)?;
        if order(sign) < order(least) {
            least = sign;
        }
    }
    Ok(least)
}

/// The sign of the turn from `p` to `q` about the origin, `p.x q.y - p.y
/// q.x`: positive where `q` is to the left of the ray through `p`.
pub fn turn(
    p: &[Rational; 2],
    q: &[Rational; 2],
    stages: &mut Stages,
    dec: &mut Decider,
) -> Result<Sign, Refusal> {
    Polynomial::Turn { p, q }.sign(stages, dec)
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::{current_operation_work, with_operation_budget};

    fn q(text: &str) -> Rational {
        Rational::parse(text).unwrap()
    }
    fn p(point: [&str; 3]) -> [Rational; 3] {
        point.map(q)
    }
    fn z_axis() -> Line {
        Line {
            origin: p(["1", "2", "0"]),
            direction: p(["0", "0", "3"]),
        }
    }

    /// The signs the stages give, straight from rational arithmetic.
    fn exact_point(line: &Line, point: &[Rational; 3], reach: &Rational) -> Sign {
        let d = &line.direction;
        let o: [Rational; 3] = std::array::from_fn(|i| &point[i] - &line.origin[i]);
        let w = [
            &o[1] * &d[2] - &o[2] * &d[1],
            &o[2] * &d[0] - &o[0] * &d[2],
            &o[0] * &d[1] - &o[1] * &d[0],
        ];
        let distance = &w[0] * &w[0] + &w[1] * &w[1] + &w[2] * &w[2];
        let length = &d[0] * &d[0] + &d[1] * &d[1] + &d[2] * &d[2];
        (distance - reach * reach * length).sign()
    }

    #[test]
    fn a_point_far_from_the_reach_is_settled_in_binary64() {
        let mut stages = Stages::default();
        let mut dec = Decider::default();
        let line = z_axis();
        let far = p(["11", "2", "7"]);
        assert_eq!(
            line.point_beyond(&far, &q("3"), &mut stages, &mut dec)
                .unwrap(),
            Sign::Positive
        );
        let inside = p(["2", "2", "-40"]);
        assert_eq!(
            line.point_beyond(&inside, &q("3"), &mut stages, &mut dec)
                .unwrap(),
            Sign::Negative
        );
        assert_eq!(
            stages,
            Stages {
                float: 2,
                interval: 0,
                exact: 0
            }
        );
    }

    #[test]
    fn a_point_at_the_reach_is_settled_exactly_and_reads_zero() {
        let mut stages = Stages::default();
        let mut dec = Decider::default();
        let line = z_axis();
        // Three and four along x and y from the line: five away.
        let on = p(["4", "6", "1/3"]);
        assert_eq!(
            line.point_beyond(&on, &q("5"), &mut stages, &mut dec)
                .unwrap(),
            Sign::Zero
        );
        assert_eq!(stages.exact, 1);
        // A reach longer by one part in 2^80 is still told apart.
        let hair = Rational::one().mul_pow2(-80);
        let longer = &q("5") + &hair;
        let shorter = &q("5") - &hair;
        assert_eq!(
            line.point_beyond(&on, &longer, &mut stages, &mut dec)
                .unwrap(),
            Sign::Negative
        );
        assert_eq!(
            line.point_beyond(&on, &shorter, &mut stages, &mut dec)
                .unwrap(),
            Sign::Positive
        );
        assert_eq!(
            stages,
            Stages {
                float: 0,
                interval: 0,
                exact: 3
            }
        );
    }

    #[test]
    fn the_interval_settles_what_the_binary64_bound_cannot() {
        // Operands that are binary64 values are read by the interval with no
        // width, and by the first stage with its bound for a rounded value.
        let mut stages = Stages::default();
        let mut dec = Decider::default();
        let one = Rational::one();
        let above = &one + &Rational::one().mul_pow2(-50);
        let sign = turn(
            &[one.clone(), one.clone()],
            &[above, one.clone()],
            &mut stages,
            &mut dec,
        )
        .unwrap();
        assert_eq!(sign, Sign::Negative);
        assert_eq!(
            stages,
            Stages {
                float: 0,
                interval: 1,
                exact: 0
            }
        );
    }

    #[test]
    fn every_stage_gives_the_sign_rational_arithmetic_gives() {
        // Points on a grid about the line, at reaches that put them outside,
        // inside, at the reach and a hair from it.
        let mut stages = Stages::default();
        let mut dec = Decider::default();
        let line = Line {
            origin: p(["1/3", "-2/7", "5"]),
            direction: p(["2", "1", "-2"]),
        };
        let hair = Rational::one().mul_pow2(-70);
        let mut checked = 0;
        for i in -3..=3 {
            for j in -3..=3 {
                let point = [
                    Rational::from_ratio(i, 3).unwrap(),
                    Rational::from_ratio(j * 5, 7).unwrap(),
                    Rational::from_ratio(i + j, 2).unwrap(),
                ];
                for reach in [q("1/2"), q("1"), q("7/3"), q("1") + hair.clone()] {
                    let want = exact_point(&line, &point, &reach);
                    let got = line
                        .point_beyond(&point, &reach, &mut stages, &mut dec)
                        .unwrap();
                    assert_eq!(got, want, "{point:?} at {reach}");
                    checked += 1;
                }
            }
        }
        assert_eq!(stages.total(), checked);
        assert!(stages.float > 0, "{stages:?}");
    }

    #[test]
    fn a_segment_is_read_at_its_nearest_point() {
        let mut stages = Stages::default();
        let mut dec = Decider::default();
        let line = z_axis();
        let reach = q("1");
        let mut beyond = |a: [&str; 3], b: [&str; 3]| {
            line.segment_beyond(&p(a), &p(b), &reach, &mut stages, &mut dec)
                .unwrap()
        };
        // Both ends far, the middle through the line.
        assert_eq!(beyond(["-9", "2", "1"], ["11", "2", "5"]), Sign::Negative);
        // Both ends far, the middle two away.
        assert_eq!(beyond(["-9", "4", "1"], ["11", "4", "5"]), Sign::Positive);
        // The middle exactly at the reach.
        assert_eq!(beyond(["-9", "3", "1"], ["11", "3", "5"]), Sign::Zero);
        // The foot beyond an end: the end decides.
        assert_eq!(beyond(["3", "2", "0"], ["9", "2", "0"]), Sign::Positive);
        assert_eq!(beyond(["9", "2", "0"], ["3/2", "2", "0"]), Sign::Negative);
        assert_eq!(beyond(["2", "2", "0"], ["9", "2", "0"]), Sign::Zero);
        // Along the line's direction.
        assert_eq!(beyond(["3", "2", "0"], ["3", "2", "8"]), Sign::Positive);
        assert_eq!(beyond(["1", "5/2", "0"], ["1", "5/2", "8"]), Sign::Negative);
    }

    #[test]
    fn a_point_against_a_point_and_against_a_segment() {
        let mut stages = Stages::default();
        let mut dec = Decider::default();
        let (a, b) = (p(["0", "0", "0"]), p(["4", "0", "0"]));
        // Three and four apart: five.
        let far = p(["3", "4", "0"]);
        let gap = |reach: &str, stages: &mut Stages, dec: &mut Decider| {
            points_beyond(&a, &far, &q(reach), stages, dec).unwrap()
        };
        assert_eq!(gap("4", &mut stages, &mut dec), Sign::Positive);
        assert_eq!(gap("5", &mut stages, &mut dec), Sign::Zero);
        assert_eq!(gap("6", &mut stages, &mut dec), Sign::Negative);
        let mut beyond = |point: [&str; 3], reach: &str| {
            point_segment_beyond(&p(point), &a, &b, &q(reach), &mut stages, &mut dec).unwrap()
        };
        // Above the middle: the foot is inside, and the height decides.
        assert_eq!(beyond(["2", "3", "0"], "5/2"), Sign::Positive);
        assert_eq!(beyond(["2", "3", "0"], "3"), Sign::Zero);
        assert_eq!(beyond(["2", "3", "0"], "7/2"), Sign::Negative);
        // Past either end: the end decides, five away.
        assert_eq!(beyond(["-3", "4", "0"], "5"), Sign::Zero);
        assert_eq!(beyond(["7", "0", "4"], "5"), Sign::Zero);
        assert_eq!(beyond(["7", "0", "4"], "4"), Sign::Positive);
    }

    #[test]
    fn two_segments_are_read_at_their_nearest_points() {
        let mut stages = Stages::default();
        let mut dec = Decider::default();
        let mut beyond = |first: [[&str; 3]; 2], second: [[&str; 3]; 2], reach: &str| {
            let [a, b] = first.map(p);
            let [c, d] = second.map(p);
            let forward =
                segments_beyond([&a, &b], [&c, &d], &q(reach), &mut stages, &mut dec).unwrap();
            let backward =
                segments_beyond([&d, &c], [&b, &a], &q(reach), &mut stages, &mut dec).unwrap();
            assert_eq!(forward, backward);
            forward
        };
        let along_x = [["-2", "0", "0"], ["2", "0", "0"]];
        // Across it, three above: the common perpendicular, between interiors.
        let across = [["0", "-2", "3"], ["0", "2", "3"]];
        assert_eq!(beyond(along_x, across, "2"), Sign::Positive);
        assert_eq!(beyond(along_x, across, "3"), Sign::Zero);
        assert_eq!(beyond(along_x, across, "4"), Sign::Negative);
        // The same line moved off the end: an end against an end, five away.
        let off = [["5", "4", "0"], ["5", "9", "0"]];
        assert_eq!(beyond(along_x, off, "5"), Sign::Zero);
        assert_eq!(beyond(along_x, off, "49/10"), Sign::Positive);
        // Parallel, two apart, overlapping in part.
        let beside = [["1", "2", "0"], ["9", "2", "0"]];
        assert_eq!(beyond(along_x, beside, "2"), Sign::Zero);
        assert_eq!(beyond(along_x, beside, "21/10"), Sign::Negative);
        // Crossing lines whose segments do not reach each other: the lines
        // meet at the origin, and the second segment starts one away.
        let short = [["0", "1", "0"], ["0", "5", "0"]];
        assert_eq!(beyond(along_x, short, "1"), Sign::Zero);
        assert_eq!(beyond(along_x, short, "1/2"), Sign::Positive);
        // Two segments that cross.
        let crossing = [["0", "-1", "0"], ["0", "1", "0"]];
        assert_eq!(beyond(along_x, crossing, "0"), Sign::Zero);
        assert_eq!(beyond(along_x, crossing, "1/1000"), Sign::Negative);
    }

    #[test]
    fn the_exact_stage_is_charged_by_the_size_of_its_operands() {
        let line = z_axis();
        let reach = q("5");
        let charged = |point: [Rational; 3]| {
            let (stages, work) = with_operation_budget(1_000_000, || {
                let mut stages = Stages::default();
                let before = current_operation_work().used;
                line.point_beyond(&point, &reach, &mut stages, &mut Decider::default())
                    .unwrap();
                (stages, current_operation_work().used - before)
            })
            .0;
            (stages, work)
        };
        let (far, cheap) = charged(p(["11", "2", "7"]));
        assert_eq!(far.float, 1);
        // One unit, and the unit of its deadline check.
        assert_eq!(cheap, 2);
        let (on, short) = charged(p(["4", "6", "1/3"]));
        assert_eq!(on.exact, 1);
        assert!(short > 4, "{short}");
        // The same point with a height of 400 bits: the cross products with
        // the line's direction multiply it, and cost more by the limb rule.
        let tall = Rational::one().mul_pow2(400) + Rational::one();
        let (long, dear) = charged([q("4"), q("6"), tall]);
        assert_eq!(long.exact, 1);
        assert!(dear > short, "{dear} against {short}");
    }

    #[test]
    fn a_stage_the_allowance_cannot_pay_for_is_not_made() {
        let line = z_axis();
        let ((result, stages), work) = with_operation_budget(3, || {
            let mut stages = Stages::default();
            let result = line.point_beyond(
                &p(["4", "6", "1/3"]),
                &q("5"),
                &mut stages,
                &mut Decider::default(),
            );
            (result, stages)
        });
        assert_eq!(result.unwrap_err().code, "numbers.operation-budget");
        assert_eq!(stages.total(), 0);
        assert!(work.exhausted);
    }
}
