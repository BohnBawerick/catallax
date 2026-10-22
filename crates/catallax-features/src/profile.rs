//! Exact Jordan-loop checks and odd-nesting regions in sketch coordinates.

use catallax_numbers::{Ball, Decider, Outcome, Rational, Real, Refusal, Sign, StrictSign};
use catallax_recipe::{ElementId, Loop, Point2, Segment, Turn};
use catallax_shapes::{Pose, RPoint3, RVector3, Rotation};

type P = [Rational; 2];
type X = [Real; 2];

/// A validated sketch, with maximal lines and circular arcs and odd-nesting material.
#[derive(Clone, Debug)]
pub struct Profile {
    pub(crate) pose: Pose,
    pub(crate) loops: Vec<Boundary>,
}

#[derive(Clone, Debug)]
pub(crate) struct Boundary {
    pub segments: Vec<Arc>,
    pub depth: usize,
    pub parent: Option<usize>,
}

#[derive(Clone, Debug)]
pub(crate) struct Arc {
    pub ids: Vec<ElementId>,
    pub start: P,
    pub end: P,
    pub circle: Option<(P, Rational, Turn)>,
    pub full: bool,
}

pub(crate) fn sign(value: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    match dec.sign(value) {
        Outcome::Proven(StrictSign::Positive) => Ok(Sign::Positive),
        Outcome::Proven(StrictSign::Negative) => Ok(Sign::Negative),
        Outcome::ProvenCoincident(_) => Ok(Sign::Zero),
        Outcome::NotExact(_) => Err(Refusal::new(
            "features.precision",
            "an exact geometric decision exhausted its budget",
        )),
        Outcome::Refused(e) => Err(e),
    }
}

fn bad(code: &'static str, detail: &str) -> Refusal {
    Refusal::new(code, detail)
}

pub(crate) fn point(p: &Point2) -> P {
    [p.0[0].value().clone(), p.0[1].value().clone()]
}

fn sub(a: &P, b: &P) -> P {
    [&a[0] - &b[0], &a[1] - &b[1]]
}
fn cross(a: &P, b: &P) -> Rational {
    &a[0] * &b[1] - &a[1] * &b[0]
}
fn dot(a: &P, b: &P) -> Rational {
    &a[0] * &b[0] + &a[1] * &b[1]
}
fn real(p: &P) -> X {
    p.clone().map(Real::from)
}
fn xsub(a: &X, b: &X) -> X {
    [&a[0] - &b[0], &a[1] - &b[1]]
}
fn xcross(a: &X, b: &X) -> Real {
    &a[0] * &b[1] - &a[1] * &b[0]
}
fn xdot(a: &X, b: &X) -> Real {
    &a[0] * &b[0] + &a[1] * &b[1]
}
fn same(a: &X, b: &P, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(
        sign(&(&a[0] - &Real::from(b[0].clone())), dec)? == Sign::Zero
            && sign(&(&a[1] - &Real::from(b[1].clone())), dec)? == Sign::Zero,
    )
}

impl Arc {
    fn reversed(&mut self) {
        std::mem::swap(&mut self.start, &mut self.end);
        if let Some((_, _, turn)) = &mut self.circle {
            *turn = match turn {
                Turn::Ccw => Turn::Cw,
                Turn::Cw => Turn::Ccw,
            };
        }
    }

    pub(crate) fn contains(&self, p: &X, dec: &mut Decider) -> Result<bool, Refusal> {
        if let Some((c, _, turn)) = &self.circle {
            if self.full {
                return Ok(true);
            }
            let (start, end) = if *turn == Turn::Ccw {
                (&self.start, &self.end)
            } else {
                (&self.end, &self.start)
            };
            let a = sub(start, c);
            let b = sub(end, c);
            let v = xsub(p, &real(c));
            let av = sign(&xcross(&real(&a), &v), dec)?;
            let vb = sign(&xcross(&v, &real(&b)), dec)?;
            Ok(match cross(&a, &b).sign() {
                Sign::Positive => av != Sign::Negative && vb != Sign::Negative,
                Sign::Negative => av != Sign::Negative || vb != Sign::Negative,
                Sign::Zero => av != Sign::Negative,
            })
        } else {
            let d = sub(&self.end, &self.start);
            let q = xsub(p, &real(&self.start));
            if sign(&xcross(&real(&d), &q), dec)? != Sign::Zero {
                return Ok(false);
            }
            let t = xdot(&q, &real(&d));
            Ok(sign(&t, dec)? != Sign::Negative
                && sign(&(t - Real::from(dot(&d, &d))), dec)? != Sign::Positive)
        }
    }

    pub(crate) fn midpoint(&self) -> X {
        if let Some((c, radius_sq, turn)) = &self.circle {
            let a = sub(&self.start, c);
            if self.full {
                return real(&self.start);
            }
            let b = sub(&self.end, c);
            let mut v = [&a[0] + &b[0], &a[1] + &b[1]];
            if v.iter().all(Rational::is_zero) {
                v = [-&a[1], a[0].clone()];
                if *turn == Turn::Cw {
                    v = v.map(|x| -x);
                }
            } else {
                let major = match cross(&a, &b).sign() {
                    Sign::Positive => *turn == Turn::Cw,
                    Sign::Negative => *turn == Turn::Ccw,
                    Sign::Zero => false,
                };
                if major {
                    v = v.map(|x| -x);
                }
            }
            let scale = (Real::from(radius_sq.clone()) / Real::from(dot(&v, &v))).sqrt();
            std::array::from_fn(|i| Real::from(c[i].clone()) + &Real::from(v[i].clone()) * &scale)
        } else {
            std::array::from_fn(|i| {
                Real::from((&self.start[i] + &self.end[i]) / Rational::from_i64(2))
            })
        }
    }
}

fn line_circle(a: &P, d: &P, c: &P, radius_sq: &Rational) -> Vec<X> {
    let v = sub(a, c);
    let aa = dot(d, d);
    let b = dot(&v, d);
    let discriminant = &b * &b - &aa * &(dot(&v, &v) - radius_sq.clone());
    if discriminant.sign() == Sign::Negative {
        return vec![];
    }
    let root = Real::from(discriminant.clone()).sqrt();
    let mut points = vec![];
    for flip in if discriminant.is_zero() {
        &[1][..]
    } else {
        &[-1, 1][..]
    } {
        let t = (Real::from(-&b) + &Real::from_i64(*flip) * &root) / Real::from(aa.clone());
        points.push(std::array::from_fn(|i| {
            Real::from(a[i].clone()) + &Real::from(d[i].clone()) * &t
        }));
    }
    points
}

fn intersections(a: &Arc, b: &Arc) -> Vec<X> {
    match (&a.circle, &b.circle) {
        (None, None) => {
            let d = sub(&a.end, &a.start);
            let e = sub(&b.end, &b.start);
            let delta = sub(&b.start, &a.start);
            let det = cross(&d, &e);
            if det.is_zero() {
                if cross(&delta, &d).is_zero() {
                    vec![
                        real(&a.start),
                        real(&a.end),
                        real(&b.start),
                        real(&b.end),
                        a.midpoint(),
                        b.midpoint(),
                    ]
                } else {
                    vec![]
                }
            } else {
                let t = cross(&delta, &e) / det;
                vec![std::array::from_fn(|i| {
                    Real::from(&a.start[i] + &(&d[i] * &t))
                })]
            }
        }
        (None, Some((c, r, _))) => line_circle(&a.start, &sub(&a.end, &a.start), c, r),
        (Some((c, r, _)), None) => line_circle(&b.start, &sub(&b.end, &b.start), c, r),
        (Some((c, r, _)), Some((q, s, _))) => {
            let d = sub(q, c);
            let n = dot(&d, &d);
            if n.is_zero() {
                if r == s {
                    vec![
                        real(&a.start),
                        real(&a.end),
                        real(&b.start),
                        real(&b.end),
                        a.midpoint(),
                        b.midpoint(),
                    ]
                } else {
                    vec![]
                }
            } else {
                let t = ((r - s) + n.clone()) / (Rational::from_i64(2) * n);
                let base = std::array::from_fn(|i| &c[i] + &(&d[i] * &t));
                let v = [-&d[1], d[0].clone()];
                line_circle(&base, &v, c, r)
            }
        }
    }
}

fn enclosed(segments: &[Arc], p: &P, dec: &mut Decider) -> Result<bool, Refusal> {
    // A rational ray avoids every vertex and every circle tangency exactly.
    let mut slope = 0_i64;
    let direction = loop {
        let d = [Rational::one(), Rational::from_i64(slope)];
        let good = segments.iter().all(|s| {
            if !s.full
                && [&s.start, &s.end]
                    .iter()
                    .any(|v| cross(&sub(v, p), &d).is_zero())
            {
                return false;
            }
            if let Some((c, r, _)) = &s.circle {
                let cp = sub(c, p);
                let det = cross(&cp, &d);
                if &det * &det == r * &dot(&d, &d) {
                    return false;
                }
            }
            true
        });
        if good {
            break d;
        }
        slope += 1;
    };
    let mut crossings = 0;
    for s in segments {
        let points = if let Some((c, r, _)) = &s.circle {
            line_circle(p, &direction, c, r)
        } else {
            let e = sub(&s.end, &s.start);
            let det = cross(&direction, &e);
            if det.is_zero() {
                continue;
            }
            let t = cross(&sub(&s.start, p), &e) / det;
            vec![std::array::from_fn(|i| {
                Real::from(&p[i] + &(&direction[i] * &t))
            })]
        };
        for q in points {
            if sign(&xdot(&xsub(&q, &real(p)), &real(&direction)), dec)? == Sign::Positive
                && s.contains(&q, dec)?
            {
                crossings += 1;
            }
        }
    }
    Ok(crossings % 2 == 1)
}

fn ccw(segments: &[Arc], dec: &Decider) -> Result<bool, Refusal> {
    // A simple Jordan loop has nonzero signed area. Refinement decides its sign
    // without rounding an arc angle or sampling its boundary into chords.
    let mut precision = 64;
    loop {
        dec.check_deadline()?;
        let mut area = Ball::from_rational(&Rational::zero(), precision);
        for s in segments {
            if let Some((c, r, turn)) = &s.circle {
                let a = sub(&s.start, c);
                let b = sub(&s.end, c);
                let determinant = cross(&a, &b);
                let scalar = dot(&a, &b);
                let pi = Ball::pi(precision);
                let period = &pi * &Ball::from_rational(&Rational::from_i64(2), precision);
                let angle = if s.full {
                    period.clone()
                } else if determinant.is_zero() {
                    pi
                } else {
                    let y = Ball::from_rational(&determinant, precision);
                    let x = Ball::from_rational(&scalar, precision);
                    let principal = Ball::atan2(&y, &x).expect("nonzero circle vectors");
                    match (turn, determinant.sign()) {
                        (Turn::Ccw, Sign::Negative) => principal + period,
                        (Turn::Cw, Sign::Positive) => principal - period,
                        _ => principal,
                    }
                };
                let angle = if (s.full || determinant.is_zero()) && *turn == Turn::Cw {
                    -angle
                } else {
                    angle
                };
                area = area
                    + Ball::from_rational(&cross(c, &sub(&s.end, &s.start)), precision)
                    + Ball::from_rational(r, precision) * angle;
            } else {
                area = area + Ball::from_rational(&cross(&s.start, &s.end), precision);
            }
        }
        if let Some(s) = area.sign() {
            return Ok(s == StrictSign::Positive);
        }
        if precision >= dec.policy().max_precision_bits {
            return Err(bad(
                "features.precision",
                "signed profile area could not be separated from zero within the precision budget",
            ));
        }
        precision = precision
            .saturating_mul(2)
            .min(dec.policy().max_precision_bits);
    }
}

pub(crate) fn placement(plane: &catallax_recipe::Pose) -> Result<Pose, Refusal> {
    let [w, x, y, z] = plane.rotation.each_ref().map(|x| x.value().clone());
    Ok(Pose::new(
        Rotation::from_quaternion(w, x, y, z)?,
        RVector3::new(
            plane.translation.0[0].value().clone(),
            plane.translation.0[1].value().clone(),
            plane.translation.0[2].value().clone(),
        ),
    ))
}

impl Profile {
    /// Validates a solved sketch. Crossing, tangency and nesting decisions are exact.
    pub fn new(
        plane: &catallax_recipe::Pose,
        loops: &[Loop],
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let pose = placement(plane)?;
        let mut boundaries = Vec::with_capacity(loops.len());
        for lp in loops {
            if lp.0.is_empty() {
                continue;
            }
            let mut segments = Vec::with_capacity(lp.0.len());
            for (i, s) in lp.0.iter().enumerate() {
                let previous = &lp.0[(i + lp.0.len() - 1) % lp.0.len()];
                let end_of = |s: &Segment| match s {
                    Segment::Line { to, .. } | Segment::Arc { to, .. } => Some(point(to)),
                    _ => None,
                };
                let arc = match s {
                    Segment::Circle { id, centre, radius } => {
                        if lp.0.len() != 1 {
                            return Err(bad(
                                "sketch.circle-loop",
                                "a full circle must be alone in its loop",
                            ));
                        }
                        if radius.value().sign() != Sign::Positive {
                            return Err(bad(
                                "sketch.non-positive-radius",
                                "circle radius must be positive",
                            ));
                        }
                        let c = point(centre);
                        let start = [&c[0] + radius.value(), c[1].clone()];
                        Arc {
                            ids: vec![id.clone()],
                            start: start.clone(),
                            end: start,
                            circle: Some((c, radius.value() * radius.value(), Turn::Ccw)),
                            full: true,
                        }
                    }
                    Segment::Line { id, to } | Segment::Arc { id, to, .. } => {
                        let start = end_of(previous).ok_or_else(|| {
                            bad("sketch.circle-loop", "a circle cannot join another segment")
                        })?;
                        let end = point(to);
                        if start == end {
                            return Err(bad(
                                "sketch.zero-segment",
                                "segment endpoints must differ",
                            ));
                        }
                        let circle = if let Segment::Arc { centre, turn, .. } = s {
                            let c = point(centre);
                            let a = sub(&start, &c);
                            let b = sub(&end, &c);
                            let r = dot(&a, &a);
                            if r != dot(&b, &b) {
                                return Err(bad(
                                    "sketch.arc-radius",
                                    "arc endpoints have different distances from the centre",
                                ));
                            }
                            Some((c, r, *turn))
                        } else {
                            None
                        };
                        Arc {
                            ids: vec![id.clone()],
                            start,
                            end,
                            circle,
                            full: false,
                        }
                    }
                };
                segments.push(arc);
            }
            boundaries.push(Boundary {
                segments,
                depth: 0,
                parent: None,
            });
        }
        for i in 0..boundaries.len() {
            for j in i..boundaries.len() {
                let a = &boundaries[i].segments;
                let b = &boundaries[j].segments;
                for (k, first) in a.iter().enumerate() {
                    for (l, second) in b.iter().enumerate() {
                        if i == j && l <= k {
                            continue;
                        }
                        let adjacent = i == j && ((k + 1) % a.len() == l || (l + 1) % a.len() == k);
                        for p in intersections(first, second) {
                            if !first.contains(&p, dec)? || !second.contains(&p, dec)? {
                                continue;
                            }
                            let endpoint = adjacent
                                && (same(&p, &first.start, dec)? || same(&p, &first.end, dec)?)
                                && (same(&p, &second.start, dec)? || same(&p, &second.end, dec)?);
                            if !endpoint {
                                return Err(bad(
                                    if i == j {
                                        "sketch.self-intersection"
                                    } else {
                                        "sketch.loop-intersection"
                                    },
                                    "sketch loops must neither cross nor touch",
                                ));
                            }
                        }
                    }
                }
            }
        }
        let mut enclosing = vec![Vec::new(); boundaries.len()];
        for i in 0..boundaries.len() {
            for j in 0..boundaries.len() {
                if i != j
                    && enclosed(
                        &boundaries[j].segments,
                        &boundaries[i].segments[0].start,
                        dec,
                    )?
                {
                    enclosing[i].push(j);
                }
            }
        }
        for i in 0..boundaries.len() {
            boundaries[i].depth = enclosing[i].len();
            boundaries[i].parent = enclosing[i]
                .iter()
                .copied()
                .max_by_key(|j| enclosing[*j].len());
            if ccw(&boundaries[i].segments, dec)? != (boundaries[i].depth % 2 == 0) {
                boundaries[i].segments.reverse();
                for s in &mut boundaries[i].segments {
                    s.reversed();
                }
            }
            merge(&mut boundaries[i].segments);
        }
        Ok(Self {
            pose,
            loops: boundaries,
        })
    }

    /// Whether the profile has no material.
    pub fn is_empty(&self) -> bool {
        self.loops.is_empty()
    }

    pub(crate) fn world(&self, p: &P, offset: &RVector3) -> RPoint3 {
        self.pose
            .apply_point(&RPoint3::new(p[0].clone(), p[1].clone(), Rational::zero()))
            .plus(offset)
    }
}

fn merge(segments: &mut Vec<Arc>) {
    loop {
        let mut found = None;
        for i in 0..segments.len() {
            let j = (i + 1) % segments.len();
            if i == j {
                break;
            }
            let a = &segments[i];
            let b = &segments[j];
            let compatible = match (&a.circle, &b.circle) {
                (Some(x), Some(y)) => x == y,
                (None, None) => {
                    let u = sub(&a.end, &a.start);
                    let v = sub(&b.end, &b.start);
                    cross(&u, &v).is_zero() && dot(&u, &v).sign() == Sign::Positive
                }
                _ => false,
            };
            if compatible {
                found = Some((i, j));
                break;
            }
        }
        let Some((i, j)) = found else {
            break;
        };
        let b = segments[j].clone();
        segments[i].end = b.end;
        segments[i].ids.extend(b.ids);
        segments[i].full = segments[i].start == segments[i].end;
        segments.remove(j);
    }
}
