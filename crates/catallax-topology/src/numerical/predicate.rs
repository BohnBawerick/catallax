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
