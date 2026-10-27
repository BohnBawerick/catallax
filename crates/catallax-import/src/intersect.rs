//! Exact self-intersection test for the planar faces of a recognised solid.
//!
//! Two faces intersect when their interiors share a positive-length segment that
//! is interior to at least one of them, or, if coplanar and facing the same way,
//! an area. Coplanar faces facing opposite ways touch and are allowed. Contact only
//! along boundaries, at a single point, or at a vertex is allowed. Everything is
//! rational arithmetic on the file's coordinates. Faces with a curved surface or
//! a circular edge are not compared.

use std::collections::BTreeSet;

use catallax_numbers::{Rational, Sign};

use crate::exact::{cross, dot, is_zero, sub};
use crate::{Curve, Face, ImportRefusal, Solid, StepFile, Surface, Vec3};

#[derive(PartialEq, Eq, Clone, Copy, Debug)]
enum Loc {
    Out,
    Boundary,
    Interior,
}

/// A planar polygon face with holes: directed boundary segments in travel order,
/// interior on the left of the face normal.
struct Poly<'a> {
    name: &'a str,
    /// Outward normal, not normalised.
    normal: Vec3,
    point: Vec3,
    segments: Vec<(Vec3, Vec3)>,
}

fn poly(face: &Face) -> Option<Poly<'_>> {
    let Surface::Plane(p) = &face.surface else {
        return None;
    };
    let mut segments = Vec::new();
    for bound in &face.loops {
        for e in &bound.edges {
            if !matches!(e.curve, Curve::Line(..)) {
                return None;
            }
            let travel = e.forward == bound.orientation;
            let (a, b) = (e.start.point.clone(), e.end.point.clone());
            segments.push(if travel { (a, b) } else { (b, a) });
        }
    }
    let normal = if face.same_sense {
        p.axis.clone()
    } else {
        p.axis.clone().map(|c| -c)
    };
    Some(Poly {
        name: &face.name,
        normal,
        point: p.origin.clone(),
        segments,
    })
}

fn on_segment(m: &Vec3, a: &Vec3, b: &Vec3) -> bool {
    is_zero(&cross(&sub(m, a), &sub(b, a))) && dot(&sub(m, a), &sub(m, b)).sign() != Sign::Positive
}

fn locate(p: &Poly, m: &Vec3) -> Loc {
    if p.segments.iter().any(|(a, b)| on_segment(m, a, b)) {
        return Loc::Boundary;
    }
    // Drop the axis along which the normal is largest, then count crossings.
    let drop = (0..3)
        .max_by_key(|&i| p.normal[i].abs())
        .expect("three axes");
    let (u, v) = ((drop + 1) % 3, (drop + 2) % 3);
    let mut inside = false;
    for (a, b) in &p.segments {
        if (a[v] > m[v]) != (b[v] > m[v]) {
            let slope = &(&m[v] - &a[v]) * &(&b[u] - &a[u]);
            let x = &a[u] + &(&slope / &(&b[v] - &a[v]));
            if m[u] < x {
                inside = !inside;
            }
        }
    }
    if inside { Loc::Interior } else { Loc::Out }
}

fn at(a: &Vec3, d: &Vec3, t: &Rational) -> Vec3 {
    std::array::from_fn(|i| &a[i] + &(&d[i] * t))
}

fn half() -> Rational {
    Rational::from_ratio(1, 2).expect("nonzero")
}

/// Whether the two faces' interiors intersect.
fn faces_intersect(a: &Poly, b: &Poly) -> bool {
    let n = cross(&a.normal, &b.normal);
    if is_zero(&n) {
        let coplanar = dot(&a.normal, &sub(&b.point, &a.point)).is_zero();
        // Faces on opposite sides of one plane touch; only faces looking the same
        // way overlap material.
        let same_way = dot(&a.normal, &b.normal).sign() == Sign::Positive;
        return coplanar && same_way && (coplanar_overlap(a, b) || coplanar_overlap(b, a));
    }
    let nn = dot(&n, &n);
    let (da, db) = (dot(&a.normal, &a.point), dot(&b.normal, &b.point));
    let (c1, c2) = (cross(&b.normal, &n), cross(&n, &a.normal));
    let p0: Vec3 = std::array::from_fn(|i| &(&(&da * &c1[i]) + &(&db * &c2[i])) / &nn);
    let param = |pt: &Vec3| &dot(&sub(pt, &p0), &n) / &nn;
    let mut ts: BTreeSet<Rational> = BTreeSet::new();
    for (poly, other) in [(a, b), (b, a)] {
        let d = dot(&other.normal, &other.point);
        for (s, e) in &poly.segments {
            let (ss, se) = (&dot(&other.normal, s) - &d, &dot(&other.normal, e) - &d);
            if ss.is_zero() && se.is_zero() {
                ts.insert(param(s));
                ts.insert(param(e));
            } else if ss.sign() != se.sign() || ss.is_zero() || se.is_zero() {
                let k = &ss / &(&ss - &se);
                ts.insert(param(&at(s, &sub(e, s), &k)));
            }
        }
    }
    let ts: Vec<Rational> = ts.into_iter().collect();
    ts.windows(2).any(|w| {
        let m = at(&p0, &n, &(&(&w[0] + &w[1]) * &half()));
        let (la, lb) = (locate(a, &m), locate(b, &m));
        la != Loc::Out && lb != Loc::Out && (la == Loc::Interior || lb == Loc::Interior)
    })
}

/// Pieces of `a`'s boundary, split where they meet `b`'s, that lie inside `b`'s
/// interior, or along a shared boundary with both interiors on the same side.
fn coplanar_overlap(a: &Poly, b: &Poly) -> bool {
    for (s, e) in &a.segments {
        let r = sub(e, s);
        let mut us: BTreeSet<Rational> = BTreeSet::new();
        us.insert(Rational::zero());
        us.insert(Rational::one());
        let rr = dot(&r, &r);
        for (c, d) in &b.segments {
            let q = sub(d, c);
            let rq = cross(&r, &q);
            if is_zero(&rq) {
                if is_zero(&cross(&sub(c, s), &r)) {
                    for pt in [c, d] {
                        let u = &dot(&sub(pt, s), &r) / &rr;
                        if u > Rational::zero() && u < Rational::one() {
                            us.insert(u);
                        }
                    }
                }
            } else {
                let den = dot(&rq, &rq);
                let cs = sub(c, s);
                let u = &dot(&cross(&cs, &q), &rq) / &den;
                let v = &dot(&cross(&cs, &r), &rq) / &den;
                if u > Rational::zero()
                    && u < Rational::one()
                    && v >= Rational::zero()
                    && v <= Rational::one()
                {
                    us.insert(u);
                }
            }
        }
        let us: Vec<Rational> = us.into_iter().collect();
        for w in us.windows(2) {
            let m = at(s, &r, &(&(&w[0] + &w[1]) * &half()));
            match locate(b, &m) {
                Loc::Interior => return true,
                Loc::Boundary => {
                    let Some((c, d)) = b.segments.iter().find(|(c, d)| on_segment(&m, c, d)) else {
                        continue;
                    };
                    let side_a = cross(&a.normal, &r);
                    let side_b = cross(&b.normal, &sub(d, c));
                    // Interiors on the same side of a shared piece overlap.
                    if is_zero(&cross(&r, &sub(d, c)))
                        && dot(&side_a, &side_b).sign() == Sign::Positive
                    {
                        return true;
                    }
                }
                Loc::Out => {}
            }
        }
    }
    false
}

/// The first pair of planar faces of the solid whose interiors intersect.
fn solid_pair(solid: &Solid) -> Option<(&str, &str)> {
    let polys: Vec<Poly> = solid.faces.iter().filter_map(poly).collect();
    for i in 0..polys.len() {
        for j in i + 1..polys.len() {
            if faces_intersect(&polys[i], &polys[j]) {
                return Some((polys[i].name, polys[j].name));
            }
        }
    }
    None
}

/// Refuses a file when two planar faces of one solid intersect. Pairs involving a
/// curved surface or circular edge are not compared. Separate solids are not compared.
pub fn check_self_intersection(file: &StepFile) -> Result<(), ImportRefusal> {
    for solid in &file.solids {
        if let Some((a, b)) = solid_pair(solid) {
            return Err(ImportRefusal::SelfIntersecting {
                solid: solid.name.clone(),
                faces: [a.to_string(), b.to_string()],
            });
        }
    }
    Ok(())
}

/// The same test for one solid.
pub(crate) fn check_solid(solid: &Solid) -> Result<(), ImportRefusal> {
    match solid_pair(solid) {
        Some((a, b)) => Err(ImportRefusal::SelfIntersecting {
            solid: solid.name.clone(),
            faces: [a.to_string(), b.to_string()],
        }),
        None => Ok(()),
    }
}
