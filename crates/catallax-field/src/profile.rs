//! Field trim predicates over validated recipe loops. Validation stays in features.

use crate::Location;
use catallax_numbers::{
    Decider, Integer, Outcome, Rational, Real, Refusal, Sign, StrictSign, UPoly,
};
use catallax_recipe::{ElementId, Loop, Segment, Turn};
use catallax_shapes::{Pose, RPoint3, RVector3, Rotation};

pub(crate) type P = [Rational; 2];
pub(crate) type X = [Real; 2];

#[derive(Clone, Debug)]
pub(crate) struct Edge {
    pub ids: Vec<ElementId>,
    pub start: P,
    pub end: P,
    pub circle: Option<(P, Rational, Turn)>,
    pub full: bool,
}

#[derive(Clone, Debug)]
pub(crate) struct Profile {
    pub pose: Pose,
    pub loops: Vec<Vec<Edge>>,
}

pub(crate) fn sign(x: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    match dec.sign(x) {
        Outcome::Proven(StrictSign::Positive) => Ok(Sign::Positive),
        Outcome::Proven(StrictSign::Negative) => Ok(Sign::Negative),
        Outcome::ProvenCoincident(_) => Ok(Sign::Zero),
        Outcome::NotExact(_) => Err(Refusal::new(
            "features.precision",
            "unresolved exact profile predicate",
        )),
        Outcome::Refused(e) => Err(e),
    }
}
pub(crate) fn sqrt(n: &Rational) -> Option<Rational> {
    fn root(n: &Integer) -> Option<Integer> {
        if n.sign() == Sign::Negative {
            return None;
        }
        if n.is_zero() {
            return Some(Integer::zero());
        }
        let mut x = Integer::one().mul_pow2(n.bit_length().div_ceil(2));
        loop {
            let y = (&x + &n.div_rem_floor(&x)?.0)
                .div_rem_floor(&Integer::from_i64(2))?
                .0;
            if y >= x {
                return (&x * &x == *n).then_some(x);
            }
            x = y;
        }
    }
    Rational::new(root(n.numerator())?, root(n.denominator())?)
}

/// Construct the nonnegative rational radicand's exact root without a floating
/// square-root filter. Root isolation and its enclosures use shared arithmetic.
pub(crate) fn positive_sqrt(n: &Rational, dec: &mut Decider) -> Result<Real, Refusal> {
    if let Some(root) = sqrt(n) {
        return Ok(Real::from(root));
    }
    match dec.real_roots(&UPoly::new(vec![-n, Rational::zero(), Rational::one()])) {
        Outcome::Proven(mut roots) => Ok(Real::from_algebraic(
            roots.pop().expect("positive radicand").value,
        )),
        Outcome::Refused(e) => Err(e),
        _ => Err(Refusal::new(
            "features.precision",
            "unresolved exact radial root",
        )),
    }
}

// Sign of a*sqrt(d)+b for d>0. Opposite signs can be compared by squaring
// without introducing another radical or a platform square root.
fn root_sign(a: &Real, b: &Real, d: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    let sa = sign(a, dec)?;
    let sb = sign(b, dec)?;
    if sa == Sign::Zero {
        return Ok(sb);
    }
    if sb == Sign::Zero || sa == sb {
        return Ok(sa);
    }
    Ok(match sign(&(a.square() * d.clone() - b.square()), dec)? {
        Sign::Positive => sa,
        Sign::Negative => sb,
        Sign::Zero => Sign::Zero,
    })
}
pub(crate) fn sub(a: &P, b: &P) -> P {
    [&a[0] - &b[0], &a[1] - &b[1]]
}
pub(crate) fn cross(a: &P, b: &P) -> Rational {
    &a[0] * &b[1] - &a[1] * &b[0]
}
pub(crate) fn dot(a: &P, b: &P) -> Rational {
    &a[0] * &b[0] + &a[1] * &b[1]
}
pub(crate) fn real(p: &P) -> X {
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

impl Edge {
    pub fn contains(&self, p: &X, dec: &mut Decider) -> Result<bool, Refusal> {
        if let Some((c, r, turn)) = &self.circle {
            let v = xsub(p, &real(c));
            if sign(&(xdot(&v, &v) - Real::from(r.clone())), dec)? != Sign::Zero {
                return Ok(false);
            }
            if self.full {
                return Ok(true);
            }
            let (a, b) = if *turn == Turn::Ccw {
                (&self.start, &self.end)
            } else {
                (&self.end, &self.start)
            };
            let a = sub(a, c);
            let b = sub(b, c);
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
}

impl Profile {
    pub fn new(plane: &catallax_recipe::Pose, loops: &[Loop]) -> Result<Self, Refusal> {
        // One validator, including tangencies, odd nesting and impossible-input codes.
        catallax_features::Profile::new(plane, loops, &mut Decider::default())?;
        let [w, x, y, z] = plane.rotation.each_ref().map(|n| n.value().clone());
        let [tx, ty, tz] = plane.translation.0.each_ref().map(|n| n.value().clone());
        let pose = Pose::new(
            Rotation::from_quaternion(w, x, y, z)?,
            RVector3::new(tx, ty, tz),
        );
        let mut boundaries = Vec::new();
        for lp in loops.iter().filter(|lp| !lp.0.is_empty()) {
            let mut edges = Vec::new();
            for (i, segment) in lp.0.iter().enumerate() {
                let p = |p: &catallax_recipe::Point2| p.0.each_ref().map(|n| n.value().clone());
                let edge = match segment {
                    Segment::Circle { id, centre, radius } => {
                        let c = p(centre);
                        let start = [&c[0] + radius.value(), c[1].clone()];
                        Edge {
                            ids: vec![id.clone()],
                            start: start.clone(),
                            end: start,
                            circle: Some((c, radius.value() * radius.value(), Turn::Ccw)),
                            full: true,
                        }
                    }
                    Segment::Line { id, to } | Segment::Arc { id, to, .. } => {
                        let previous = &lp.0[(i + lp.0.len() - 1) % lp.0.len()];
                        let start = match previous {
                            Segment::Line { to, .. } | Segment::Arc { to, .. } => p(to),
                            _ => unreachable!("validated loop"),
                        };
                        let end = p(to);
                        let circle = if let Segment::Arc { centre, turn, .. } = segment {
                            let c = p(centre);
                            let d = sub(&start, &c);
                            Some((c, dot(&d, &d), *turn))
                        } else {
                            None
                        };
                        Edge {
                            ids: vec![id.clone()],
                            start,
                            end,
                            circle,
                            full: false,
                        }
                    }
                };
                edges.push(edge);
            }
            // The field retains all aliases of a maximal straight or circular run.
            loop {
                let pair = (0..edges.len()).find_map(|i| {
                    let j = (i + 1) % edges.len();
                    if i == j {
                        return None;
                    }
                    let (a, b) = (&edges[i], &edges[j]);
                    let compatible = match (&a.circle, &b.circle) {
                        (Some(a), Some(b)) => a == b,
                        (None, None) => {
                            let u = sub(&a.end, &a.start);
                            let v = sub(&b.end, &b.start);
                            cross(&u, &v).is_zero() && dot(&u, &v).sign() == Sign::Positive
                        }
                        _ => false,
                    };
                    compatible.then_some((i, j))
                });
                let Some((i, j)) = pair else {
                    break;
                };
                let b = edges[j].clone();
                edges[i].end = b.end;
                edges[i].ids.extend(b.ids);
                edges[i].full = edges[i].start == edges[i].end;
                edges.remove(j);
            }
            boundaries.push(edges);
        }
        Ok(Self {
            pose,
            loops: boundaries,
        })
    }
    pub fn world(&self, p: &P) -> RPoint3 {
        self.pose
            .apply_point(&RPoint3::new(p[0].clone(), p[1].clone(), Rational::zero()))
    }
}

pub(crate) fn location(loops: &[Vec<Edge>], p: &X, dec: &mut Decider) -> Result<Location, Refusal> {
    let mut inside = false;
    for edges in loops {
        match loop_location(edges, p, dec)? {
            Location::Boundary => return Ok(Location::Boundary),
            Location::Inside => inside = !inside,
            _ => (),
        }
    }
    Ok(if inside {
        Location::Inside
    } else {
        Location::Outside
    })
}

pub(crate) fn loop_location(edges: &[Edge], p: &X, dec: &mut Decider) -> Result<Location, Refusal> {
    for e in edges {
        if e.contains(p, dec)? {
            return Ok(Location::Boundary);
        }
    }
    // Choose a ray missing every vertex and circle tangency. The finite forbidden
    // slopes cannot exhaust the deterministic nonnegative integer enumeration.
    let mut slope = 0_i64;
    let direction = loop {
        let d = [Rational::one(), Rational::from_i64(slope)];
        let mut good = true;
        for e in edges {
            if !e.full {
                for v in [&e.start, &e.end] {
                    if sign(&xcross(&xsub(&real(v), p), &real(&d)), dec)? == Sign::Zero {
                        good = false;
                    }
                }
            }
            if let Some((c, r, _)) = &e.circle {
                let det = xcross(&xsub(&real(c), p), &real(&d));
                if sign(&(&det * &det - Real::from(r * &dot(&d, &d))), dec)? == Sign::Zero {
                    good = false;
                }
            }
        }
        if good {
            break d;
        }
        slope += 1;
    };
    let d = real(&direction);
    let mut crossings = 0;
    for e in edges {
        if let Some((c, r, _)) = &e.circle {
            let v = xsub(p, &real(c));
            let aa = Real::from(dot(&direction, &direction));
            let b = xdot(&v, &d);
            let discriminant = &b * &b - &aa * &(xdot(&v, &v) - Real::from(r.clone()));
            if sign(&discriminant, dec)? != Sign::Positive {
                continue;
            }
            for flip in [-1, 1] {
                let flip = Real::from_i64(flip);
                if root_sign(&flip, &(-&b), &discriminant, dec)? != Sign::Positive {
                    continue;
                }
                if e.full {
                    crossings += 1;
                    continue;
                }
                let base: X = std::array::from_fn(|i| &v[i] - &((&d[i] * &b) / aa.clone()));
                let coefficient: X = std::array::from_fn(|i| (&d[i] * &flip) / aa.clone());
                let (_, _, turn) = e.circle.as_ref().expect("circle branch");
                let (start, end) = if *turn == Turn::Ccw {
                    (&e.start, &e.end)
                } else {
                    (&e.end, &e.start)
                };
                let start = sub(start, c);
                let end = sub(end, c);
                let av = root_sign(
                    &xcross(&real(&start), &coefficient),
                    &xcross(&real(&start), &base),
                    &discriminant,
                    dec,
                )?;
                let vb = root_sign(
                    &xcross(&coefficient, &real(&end)),
                    &xcross(&base, &real(&end)),
                    &discriminant,
                    dec,
                )?;
                let contains = match cross(&start, &end).sign() {
                    Sign::Positive => av != Sign::Negative && vb != Sign::Negative,
                    Sign::Negative => av != Sign::Negative || vb != Sign::Negative,
                    Sign::Zero => av != Sign::Negative,
                };
                if contains {
                    crossings += 1;
                }
            }
        } else {
            let edge = sub(&e.end, &e.start);
            let det = cross(&direction, &edge);
            if det.is_zero() {
                continue;
            }
            let t = xcross(&xsub(&real(&e.start), p), &real(&edge)) / Real::from(det);
            if sign(&t, dec)? != Sign::Positive {
                continue;
            }
            let q = std::array::from_fn(|i| &p[i] + &(&d[i] * &t));
            if e.contains(&q, dec)? {
                crossings += 1;
            }
        }
    }
    Ok(if crossings % 2 == 1 {
        Location::Inside
    } else {
        Location::Outside
    })
}
