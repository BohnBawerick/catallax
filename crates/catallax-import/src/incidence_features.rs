//! Boundary feature bounds for complete simple planar regions, including holes.

use catallax_numbers::{Rational, Sign};

use crate::exact::{cross, dot, is_zero, sub};
use crate::{Curve, Face, Surface, Vec3};

fn available() -> Result<(), String> {
    if catallax_numbers::current_operation_work().exhausted {
        Err("the feature proof exhausted the operation allowance".into())
    } else {
        Ok(())
    }
}

fn point_segment2(p: &Vec3, a: &Vec3, b: &Vec3) -> Rational {
    let d = sub(b, a);
    let t = (dot(&sub(p, a), &d) / dot(&d, &d))
        .max(Rational::zero())
        .min(Rational::one());
    let v: Vec3 = std::array::from_fn(|i| &p[i] - &(&a[i] + &(&t * &d[i])));
    dot(&v, &v)
}

fn segment_distance2(a: &Vec3, b: &Vec3, c: &Vec3, d: &Vec3) -> Rational {
    let mut answer = point_segment2(a, c, d)
        .min(point_segment2(b, c, d))
        .min(point_segment2(c, a, b))
        .min(point_segment2(d, a, b));
    let u = sub(b, a);
    let v = sub(d, c);
    let w = sub(a, c);
    let aa = dot(&u, &u);
    let bb = dot(&u, &v);
    let cc = dot(&v, &v);
    let dd = dot(&u, &w);
    let ee = dot(&v, &w);
    let det = &aa * &cc - &bb * &bb;
    if det > Rational::zero() {
        let s = (&bb * &ee - &cc * &dd) / det.clone();
        let t = (&aa * &ee - &bb * &dd) / det;
        if s >= Rational::zero()
            && s <= Rational::one()
            && t >= Rational::zero()
            && t <= Rational::one()
        {
            let r: Vec3 = std::array::from_fn(|i| &(&w[i] + &(&s * &u[i])) - &(&t * &v[i]));
            answer = answer.min(dot(&r, &r));
        }
    }
    answer
}

fn inside(point: &Vec3, ring: &[Vec3], drop: usize) -> bool {
    let (u, v) = ((drop + 1) % 3, (drop + 2) % 3);
    let mut inside = false;
    for i in 0..ring.len() {
        let a = &ring[i];
        let b = &ring[(i + 1) % ring.len()];
        if (a[v] > point[v]) != (b[v] > point[v]) {
            let x = &a[u] + &(&(&point[v] - &a[v]) * &(&b[u] - &a[u]) / (&b[v] - &a[v]));
            if point[u] < x {
                inside = !inside;
            }
        }
    }
    inside
}

struct Segment {
    a: Vec3,
    b: Vec3,
    ring: usize,
    index: usize,
    length2: Rational,
}

/// Proves simple, consistently oriented boundaries with one outer region and
/// disjoint holes. The lower bound includes every edge and separated boundary
/// feature. Distances use finite segments, so extensions cannot introduce a
/// fictitious narrow feature outside the repaired face.
pub(crate) fn feature_size2(face: &Face) -> Result<Rational, String> {
    let Surface::Plane(plane) = &face.surface else {
        return Err("the boundary feature proof needs a plane".into());
    };
    if is_zero(&plane.axis) {
        return Err("the boundary feature proof needs a nonzero normal".into());
    }
    let mut rings = Vec::new();
    for bound in &face.loops {
        if bound.vertex.is_some() {
            return Err("a planar feature proof cannot cover a vertex loop".into());
        }
        for chain in crate::build::chains(&face.name, bound.orientation, &bound.edges)
            .map_err(|e| format!("boundary feature proof: {e:?}"))?
        {
            if chain.len() < 3 {
                return Err("a planar feature proof needs closed nondegenerate polygons".into());
            }
            let mut points = Vec::new();
            for (edge, forward) in chain {
                available()?;
                let Curve::Line(o, d) = &edge.curve else {
                    return Err("the feature proof needs certified straight repaired edges".into());
                };
                let (a, b) = if forward {
                    (&edge.start.point, &edge.end.point)
                } else {
                    (&edge.end.point, &edge.start.point)
                };
                if a == b
                    || is_zero(d)
                    || !is_zero(&cross(&sub(a, o), d))
                    || !is_zero(&cross(&sub(b, o), d))
                    || !dot(&sub(a, &plane.origin), &plane.axis).is_zero()
                    || !dot(&sub(b, &plane.origin), &plane.axis).is_zero()
                {
                    return Err("a boundary feature lacks exact finite planar incidence".into());
                }
                points.push(a.clone());
            }
            rings.push(points);
        }
    }
    if rings.is_empty() {
        return Err("a planar feature proof needs a finite boundary".into());
    }
    let mut segments = Vec::new();
    for (ring, points) in rings.iter().enumerate() {
        for index in 0..points.len() {
            let a = points[index].clone();
            let b = points[(index + 1) % points.len()].clone();
            let d = sub(&b, &a);
            segments.push(Segment {
                a,
                b,
                ring,
                index,
                length2: dot(&d, &d),
            });
        }
    }
    let mut feature = segments
        .iter()
        .map(|s| s.length2.clone())
        .min()
        .ok_or("empty boundary feature set")?;
    for i in 0..segments.len() {
        for j in i + 1..segments.len() {
            available()?;
            let a = &segments[i];
            let b = &segments[j];
            let n = rings[a.ring].len();
            let adjacent =
                a.ring == b.ring && ((a.index + 1) % n == b.index || (b.index + 1) % n == a.index);
            if adjacent {
                let shared = if a.a == b.a || a.a == b.b {
                    &a.a
                } else if a.b == b.a || a.b == b.b {
                    &a.b
                } else {
                    return Err("adjacent repaired boundary edges do not meet".into());
                };
                let u = sub(if &a.a == shared { &a.b } else { &a.a }, shared);
                let v = sub(if &b.a == shared { &b.b } else { &b.a }, shared);
                if is_zero(&cross(&u, &v)) && dot(&u, &v) > Rational::zero() {
                    return Err("adjacent repaired boundary edges overlap".into());
                }
            } else {
                let distance = segment_distance2(&a.a, &a.b, &b.a, &b.b);
                if distance.is_zero() {
                    return Err("repaired boundary loops cross or touch".into());
                }
                feature = feature.min(distance);
            }
        }
    }
    let drop = (0..3).max_by_key(|&i| plane.axis[i].abs()).unwrap();
    let mut outer = 0;
    for (i, ring) in rings.iter().enumerate() {
        available()?;
        let depth = rings
            .iter()
            .enumerate()
            .filter(|(j, other)| *j != i && inside(&ring[0], other, drop))
            .count();
        if depth == 0 {
            outer += 1;
        } else if depth != 1 {
            return Err("the planar face has nested islands rather than disjoint holes".into());
        }
        let area = ring
            .iter()
            .enumerate()
            .fold(Rational::zero(), |sum, (k, a)| {
                sum + dot(
                    &cross(
                        &sub(a, &plane.origin),
                        &sub(&ring[(k + 1) % ring.len()], &plane.origin),
                    ),
                    &plane.axis,
                )
            });
        let outward = if face.same_sense { area } else { -area };
        if outward.sign()
            != if depth == 0 {
                Sign::Positive
            } else {
                Sign::Negative
            }
        {
            return Err("a repaired boundary has unproved face orientation".into());
        }
        for point in ring {
            for edge in &segments {
                available()?;
                if point == &edge.a || point == &edge.b {
                    continue;
                }
                feature = feature.min(point_segment2(point, &edge.a, &edge.b));
            }
        }
    }
    if outer != 1 || feature <= Rational::zero() {
        return Err(
            "the planar face lacks a single covered outer region and positive feature bound".into(),
        );
    }
    Ok(feature)
}
