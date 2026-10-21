//! Piecewise degree-seven Taylor polynomials, converted to Bezier spans.
//! Each span carries a uniform Taylor remainder and a convex-hull bound for
//! decimal rounding, enclosure width, and endpoint replacement.

use catallax_intersect::{CurveTrace, TracePiece};
use catallax_numbers::{Ball, Decider, Rational, Refusal};

use crate::StepOptions;
use crate::number::{self, Point};

const DEGREE: usize = 7;

pub(crate) struct Fit {
    pub(crate) points: Vec<Point>,
    pub(crate) degree: usize,
    pub(crate) knots: Vec<Rational>,
    pub(crate) multiplicities: Vec<usize>,
    pub(crate) error: Rational,
    pub(crate) bounds: [[Rational; 2]; 3],
}

struct Span {
    points: Vec<Point>,
    ideals: Vec<[[Rational; 2]; 3]>,
    remainder: Rational,
    bounds: [[Rational; 2]; 3],
}

impl Span {
    fn error(&self) -> Rational {
        let controls = self
            .points
            .iter()
            .zip(&self.ideals)
            .map(|(p, b)| {
                (0..3).fold(Rational::zero(), |s, i| {
                    s + (&p.xyz[i].value - &b[i][0])
                        .abs()
                        .max((&p.xyz[i].value - &b[i][1]).abs())
                })
            })
            .max()
            .unwrap();
        &self.remainder + &controls
    }
}

pub(crate) fn fit(
    trace: &CurveTrace,
    start: &Point,
    end: &Point,
    reverse: bool,
    options: &StepOptions,
    dec: &Decider,
) -> Result<Fit, Refusal> {
    dec.check_deadline()?;
    let policy = dec.policy();
    assert!(!trace.pieces().is_empty(), "empty curve trace");
    // TracePiece's internal arithmetic is excluded, as allowed by piece 6's
    // convention. Each span constructs five input balls and one hull.
    let pieces = u64::try_from(trace.pieces().len()).unwrap();
    let initial_operations = pieces
        .checked_mul(6)
        .ok_or_else(|| number::budget("too many initial spline spans"))?;
    let p = number::work_precision(
        policy,
        options
            .significant_digits
            .max(1)
            .saturating_mul(4)
            .saturating_add(32),
        initial_operations,
    )
    .ok_or_else(|| number::budget("no work available for an initial spline fit"))?;
    let cost = u64::from(p.div_ceil(64)).pow(2) * 6;
    // Reserve every initial span before spending anything on refinement.
    let mut remaining = policy.max_work - cost * pieces;
    let mut spans = Vec::new();
    for piece in trace.pieces() {
        let (lo, hi) = piece.range();
        let whole = Ball::from_rational(&lo, p).hull(&Ball::from_rational(&hi, p));
        if whole.lo() < lo || whole.hi() > hi {
            return Err(number::budget(
                "trace range cannot be enclosed at working precision",
            ));
        }
        let initial = span(piece, &lo, &hi, p, options.significant_digits, &whole, dec)?;
        // Its work was reserved above.
        let mut pending = vec![(lo, hi, initial)];
        while let Some((lo, hi, current)) = pending.pop() {
            dec.check_deadline()?;
            let error = current.error();
            // Once decimal/enclosure error dominates, subdivision cannot help.
            let useful = &current.remainder * &Rational::from_i64(4)
                > error.clone().max(options.curve_tolerance.clone());
            if error > options.curve_tolerance && useful && remaining >= cost.saturating_mul(2) {
                remaining -= cost * 2;
                let mid = (&lo + &hi).mul_pow2(-1);
                let left = span(piece, &lo, &mid, p, options.significant_digits, &whole, dec)?;
                let right = span(piece, &mid, &hi, p, options.significant_digits, &whole, dec)?;
                if left.error().max(right.error()) < error {
                    pending.push((mid.clone(), hi, right));
                    pending.push((lo, mid, left));
                } else {
                    spans.push(current);
                }
            } else {
                spans.push(current);
            }
        }
    }
    let (first, last) = if reverse { (end, start) } else { (start, end) };
    spans.first_mut().unwrap().points[0] = first.clone();
    spans.last_mut().unwrap().points[DEGREE] = last.clone();
    let mut points: Vec<Point> = Vec::with_capacity(spans.len() * DEGREE + 1);
    let mut error = Rational::zero();
    let mut bounds = spans[0].bounds.clone();
    for mut span in spans {
        if let Some(previous) = points.last() {
            span.points[0] = previous.clone();
        }
        error = error.max(span.error());
        for (b, other) in bounds.iter_mut().zip(&span.bounds) {
            b[0] = b[0].clone().min(other[0].clone());
            b[1] = b[1].clone().max(other[1].clone());
        }
        let skip = usize::from(!points.is_empty());
        points.extend(span.points.into_iter().skip(skip));
    }
    if reverse {
        points.reverse();
    }
    let count = (points.len() - 1) / DEGREE;
    let mut knots = Vec::with_capacity(count + 1);
    let mut knot = Rational::zero();
    knots.push(knot.clone());
    for _ in 0..count {
        // A repeated Bezier knot isolates each span's parameterization. Choose
        // the next exactly representable decimal so low digit counts cannot
        // collapse distinct knots and change the fitted curve.
        let next = (&knot + &Rational::one())
            .to_decimal_bounds(options.significant_digits.max(1))
            .1;
        knot = number::formatted_value(&next);
        knots.push(knot.clone());
    }
    let multiplicities = (0..=count)
        .map(|i| {
            if i == 0 || i == count {
                DEGREE + 1
            } else {
                DEGREE
            }
        })
        .collect();
    Ok(Fit {
        points,
        degree: DEGREE,
        knots,
        multiplicities,
        error,
        bounds,
    })
}

fn span(
    piece: &TracePiece,
    lo: &Rational,
    hi: &Rational,
    precision: u32,
    digits: u32,
    whole: &Ball,
    dec: &Decider,
) -> Result<Span, Refusal> {
    let mid = (lo + hi).mul_pow2(-1);
    let width = hi - lo;
    let lower = whole.lo();
    let upper = whole.hi();
    let inside = |ball: &Ball| ball.lo() >= lower && ball.hi() <= upper;
    let raw_center = Ball::from_rational(&mid, precision);
    let center = if inside(&raw_center) {
        &raw_center
    } else {
        whole
    };
    // Outward rounding at low precision may leave a piece's range. Fall
    // back to the certified whole range rather than evaluate beyond it.
    let candidate = (lo != &lower || hi != &upper)
        .then(|| Ball::from_rational(lo, precision).hull(&Ball::from_rational(hi, precision)));
    let domain = candidate
        .as_ref()
        .filter(|ball| inside(ball))
        .unwrap_or(whole);
    let coefficients = piece.taylor(center, DEGREE, dec)?;
    let remainder_coefficient = piece.taylor(domain, DEGREE + 1, dec)?.pop().unwrap();
    let remainder = remainder_coefficient
        .iter()
        .fold(Rational::zero(), |s, b| s + number::magnitude(b))
        * width.mul_pow2(-1).pow((DEGREE + 1) as i32).unwrap();
    let mut ideals = Vec::with_capacity(DEGREE + 1);
    let mut points = Vec::with_capacity(DEGREE + 1);
    for j in 0..=DEGREE {
        let mut control: [[Rational; 2]; 3] =
            std::array::from_fn(|_| [Rational::zero(), Rational::zero()]);
        for (k, coefficient) in coefficients.iter().enumerate() {
            let factor = (0..=k.min(j)).fold(Rational::zero(), |s, l| {
                let term =
                    Rational::from_ratio(binomial(k, l) * binomial(j, l), binomial(DEGREE, l))
                        .unwrap()
                        * Rational::from_ratio(-1, 2)
                            .unwrap()
                            .pow((k - l) as i32)
                            .unwrap();
                s + term
            }) * width.pow(k as i32).unwrap();
            for i in 0..3 {
                let a = &coefficient[i].lo() * &factor;
                let b = &coefficient[i].hi() * &factor;
                control[i][0] = &control[i][0] + &a.clone().min(b.clone());
                control[i][1] = &control[i][1] + &a.max(b);
            }
        }
        points.push(Point {
            xyz: control
                .each_ref()
                .map(|b| number::rounded(&b[0], &b[1], digits)),
        });
        ideals.push(control);
    }
    for (index, t) in [(0, lo), (DEGREE, hi)] {
        let raw = Ball::from_rational(t, precision);
        let parameter = if inside(&raw) { &raw } else { whole };
        points[index] = Point {
            xyz: piece
                .evaluate(parameter, dec)?
                .each_ref()
                .map(|b| number::ball_number(b, digits)),
        };
    }
    let bounds = piece
        .evaluate(domain, dec)?
        .each_ref()
        .map(|b| [b.lo(), b.hi()]);
    Ok(Span {
        points,
        ideals,
        remainder,
        bounds,
    })
}

fn binomial(n: usize, k: usize) -> i64 {
    (0..k).fold(1, |v, i| v * (n - i) as i64 / (i + 1) as i64)
}
