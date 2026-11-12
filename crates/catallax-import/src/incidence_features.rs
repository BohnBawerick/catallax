//! Boundary feature bounds for complete simple planar regions, including holes.

use catallax_numbers::{Decider, Rational, Sign};

use crate::exact::is_zero;
use crate::incidence_arithmetic::{
    arithmetic, checkpoint, compare, cross, determinant, dot, minimum, same, sub,
};
use crate::{Curve, Face, Surface, Vec3};

fn point_segment2(p: &Vec3, a: &Vec3, b: &Vec3, dec: &mut Decider) -> Result<Rational, String> {
    let d = sub(b, a, dec)?;
    let delta = sub(p, a, dec)?;
    let numerator = dot(&delta, &d, dec)?;
    let denominator = dot(&d, &d, dec)?;
    let t = arithmetic(&numerator, &denominator, '/', dec)?;
    let t = if t.sign() == Sign::Negative {
        Rational::zero()
    } else {
        minimum(t, Rational::one(), dec)?
    };
    let mut v = delta;
    for i in 0..3 {
        let scaled = arithmetic(&t, &d[i], '*', dec)?;
        v[i] = arithmetic(&v[i], &scaled, '-', dec)?;
    }
    dot(&v, &v, dec)
}

fn segment_distance2(
    a: &Vec3,
    b: &Vec3,
    c: &Vec3,
    d: &Vec3,
    dec: &mut Decider,
) -> Result<Rational, String> {
    let mut answer = point_segment2(a, c, d, dec)?;
    for (p, x, y) in [(b, c, d), (c, a, b), (d, a, b)] {
        let distance = point_segment2(p, x, y, dec)?;
        answer = minimum(answer, distance, dec)?;
    }
    let u = sub(b, a, dec)?;
    let v = sub(d, c, dec)?;
    let w = sub(a, c, dec)?;
    let aa = dot(&u, &u, dec)?;
    let bb = dot(&u, &v, dec)?;
    let cc = dot(&v, &v, dec)?;
    let dd = dot(&u, &w, dec)?;
    let ee = dot(&v, &w, dec)?;
    let det = determinant(&aa, &cc, &bb, &bb, dec)?;
    if det.sign() == Sign::Positive {
        let sn = determinant(&bb, &ee, &cc, &dd, dec)?;
        let tn = determinant(&aa, &ee, &bb, &dd, dec)?;
        let s = arithmetic(&sn, &det, '/', dec)?;
        let t = arithmetic(&tn, &det, '/', dec)?;
        if s.sign() != Sign::Negative
            && compare(&s, &Rational::one(), dec)?.is_le()
            && t.sign() != Sign::Negative
            && compare(&t, &Rational::one(), dec)?.is_le()
        {
            let mut r = w;
            for i in 0..3 {
                let su = arithmetic(&s, &u[i], '*', dec)?;
                let tv = arithmetic(&t, &v[i], '*', dec)?;
                r[i] = arithmetic(&r[i], &su, '+', dec)?;
                r[i] = arithmetic(&r[i], &tv, '-', dec)?;
            }
            let distance = dot(&r, &r, dec)?;
            answer = minimum(answer, distance, dec)?;
        }
    }
    Ok(answer)
}

fn inside(point: &Vec3, ring: &[Vec3], drop: usize, dec: &mut Decider) -> Result<bool, String> {
    let (u, v) = ((drop + 1) % 3, (drop + 2) % 3);
    let mut inside = false;
    for i in 0..ring.len() {
        checkpoint(dec)?;
        let a = &ring[i];
        let b = &ring[(i + 1) % ring.len()];
        if compare(&a[v], &point[v], dec)?.is_gt() != compare(&b[v], &point[v], dec)?.is_gt() {
            let height = arithmetic(&point[v], &a[v], '-', dec)?;
            let width = arithmetic(&b[u], &a[u], '-', dec)?;
            let rise = arithmetic(&b[v], &a[v], '-', dec)?;
            let offset = arithmetic(&height, &width, '*', dec)?;
            let offset = arithmetic(&offset, &rise, '/', dec)?;
            let x = arithmetic(&a[u], &offset, '+', dec)?;
            if compare(&point[u], &x, dec)?.is_lt() {
                inside = !inside;
            }
        }
    }
    Ok(inside)
}

struct Segment {
    a: Vec3,
    b: Vec3,
    ring: usize,
    index: usize,
}

/// Proves simple, consistently oriented boundaries with one outer region and
/// disjoint holes. Distances use finite segments, so extensions cannot introduce
/// a fictitious narrow feature outside the repaired face.
pub(crate) fn feature_size2(face: &Face) -> Result<Rational, String> {
    let mut dec = Decider::default();
    checkpoint(&mut dec)?;
    let Surface::Plane(plane) = &face.surface else {
        return Err("the boundary feature proof needs a plane".into());
    };
    if is_zero(&plane.axis) {
        return Err("the boundary feature proof needs a nonzero normal".into());
    }
    let mut rings = Vec::new();
    for bound in &face.loops {
        checkpoint(&mut dec)?;
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
                checkpoint(&mut dec)?;
                let Curve::Line(o, d) = &edge.curve else {
                    return Err("the feature proof needs certified straight repaired edges".into());
                };
                let (a, b) = if forward {
                    (&edge.start.point, &edge.end.point)
                } else {
                    (&edge.end.point, &edge.start.point)
                };
                let ao = sub(a, o, &mut dec)?;
                let bo = sub(b, o, &mut dec)?;
                let ap = sub(a, &plane.origin, &mut dec)?;
                let bp = sub(b, &plane.origin, &mut dec)?;
                if same(a, b, &mut dec)?
                    || is_zero(d)
                    || !is_zero(&cross(&ao, d, &mut dec)?)
                    || !is_zero(&cross(&bo, d, &mut dec)?)
                    || !dot(&ap, &plane.axis, &mut dec)?.is_zero()
                    || !dot(&bp, &plane.axis, &mut dec)?.is_zero()
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
    let mut feature: Option<Rational> = None;
    for (ring, points) in rings.iter().enumerate() {
        for index in 0..points.len() {
            checkpoint(&mut dec)?;
            let a = points[index].clone();
            let b = points[(index + 1) % points.len()].clone();
            let d = sub(&b, &a, &mut dec)?;
            let length2 = dot(&d, &d, &mut dec)?;
            feature = Some(if let Some(previous) = feature {
                minimum(previous, length2.clone(), &mut dec)?
            } else {
                length2.clone()
            });
            segments.push(Segment { a, b, ring, index });
        }
    }
    let mut feature = feature.ok_or("empty boundary feature set")?;
    for i in 0..segments.len() {
        for j in i + 1..segments.len() {
            checkpoint(&mut dec)?;
            let a = &segments[i];
            let b = &segments[j];
            let n = rings[a.ring].len();
            let adjacent =
                a.ring == b.ring && ((a.index + 1) % n == b.index || (b.index + 1) % n == a.index);
            if adjacent {
                let shared = if same(&a.a, &b.a, &mut dec)? || same(&a.a, &b.b, &mut dec)? {
                    &a.a
                } else if same(&a.b, &b.a, &mut dec)? || same(&a.b, &b.b, &mut dec)? {
                    &a.b
                } else {
                    return Err("adjacent repaired boundary edges do not meet".into());
                };
                let u = sub(
                    if same(&a.a, shared, &mut dec)? {
                        &a.b
                    } else {
                        &a.a
                    },
                    shared,
                    &mut dec,
                )?;
                let v = sub(
                    if same(&b.a, shared, &mut dec)? {
                        &b.b
                    } else {
                        &b.a
                    },
                    shared,
                    &mut dec,
                )?;
                if is_zero(&cross(&u, &v, &mut dec)?)
                    && dot(&u, &v, &mut dec)?.sign() == Sign::Positive
                {
                    return Err("adjacent repaired boundary edges overlap".into());
                }
            } else {
                let distance = segment_distance2(&a.a, &a.b, &b.a, &b.b, &mut dec)?;
                if distance.is_zero() {
                    return Err("repaired boundary loops cross or touch".into());
                }
                feature = minimum(feature, distance, &mut dec)?;
            }
        }
    }
    let mut drop = 0;
    for i in 1..3 {
        if compare(&plane.axis[i].abs(), &plane.axis[drop].abs(), &mut dec)?.is_ge() {
            drop = i;
        }
    }
    let mut outer = 0;
    for (i, ring) in rings.iter().enumerate() {
        checkpoint(&mut dec)?;
        let mut depth = 0;
        for (j, other) in rings.iter().enumerate() {
            if j != i && inside(&ring[0], other, drop, &mut dec)? {
                depth += 1;
            }
        }
        if depth == 0 {
            outer += 1;
        } else if depth != 1 {
            return Err("the planar face has nested islands rather than disjoint holes".into());
        }
        let mut area = Rational::zero();
        for (k, a) in ring.iter().enumerate() {
            checkpoint(&mut dec)?;
            let u = sub(a, &plane.origin, &mut dec)?;
            let v = sub(&ring[(k + 1) % ring.len()], &plane.origin, &mut dec)?;
            let cross = cross(&u, &v, &mut dec)?;
            let product = dot(&cross, &plane.axis, &mut dec)?;
            area = arithmetic(&area, &product, '+', &mut dec)?;
        }
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
                checkpoint(&mut dec)?;
                if same(point, &edge.a, &mut dec)? || same(point, &edge.b, &mut dec)? {
                    continue;
                }
                let distance = point_segment2(point, &edge.a, &edge.b, &mut dec)?;
                feature = minimum(feature, distance, &mut dec)?;
            }
        }
    }
    if outer != 1 || feature.sign() != Sign::Positive {
        return Err(
            "the planar face lacks a single covered outer region and positive feature bound".into(),
        );
    }
    Ok(feature)
}
#[cfg(test)]
#[path = "../tests/common/mod.rs"]
mod common;

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::with_operation_budget;

    #[test]
    fn exact_boundary_work_is_charged_and_stops_at_the_operation_allowance() {
        let text = common::box_step(common::BoxFixture::new(["2", "3", "4"]));
        let file = crate::read_step(text.as_bytes()).unwrap();
        let face = &file.solids[0].faces[0];
        let (result, work) = with_operation_budget(20_000_000, || feature_size2(face));
        assert!(result.is_ok(), "{result:?}");
        assert!(
            work.used > 0,
            "exact boundary predicates must charge their work"
        );
        assert!(!work.exhausted);
        let (refused, limited) = with_operation_budget(work.used / 2, || feature_size2(face));
        assert!(
            refused.is_err(),
            "the proof must stop before returning a bound"
        );
        assert!(refused.unwrap_err().contains("numbers.operation-budget"));
        assert!(limited.exhausted);
        assert!(limited.used <= work.used / 2);
    }

    #[test]
    fn exact_boundary_work_accounts_for_coordinate_bit_size() {
        let make_face = |fixture| {
            let text = common::box_step(fixture);
            crate::read_step(text.as_bytes())
                .unwrap()
                .solids
                .remove(0)
                .faces
                .remove(0)
        };
        let small = make_face(common::BoxFixture::new(["2", "3", "4"]));
        let mut fixture = common::BoxFixture::new(["2", "3", "4"]);
        let offset = Rational::parse(&format!("0.{}1", "0".repeat(300))).unwrap();
        fixture.origin = [offset.clone(), offset.clone(), offset];
        let large = make_face(fixture);
        let (a, small_work) = with_operation_budget(20_000_000, || feature_size2(&small));
        let (b, large_work) = with_operation_budget(20_000_000, || feature_size2(&large));
        assert_eq!(
            a.unwrap(),
            b.unwrap(),
            "translation preserves physical feature size"
        );
        assert!(!large_work.exhausted);
        assert!(
            large_work.used > small_work.used,
            "larger exact operands must cost more"
        );
        let (refused, limited) = with_operation_budget(small_work.used, || feature_size2(&large));
        assert!(refused.is_err());
        assert!(limited.exhausted);
    }
}
