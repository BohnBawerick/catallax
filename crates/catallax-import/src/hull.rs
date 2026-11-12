//! Spline face pairs kept apart by their control hulls.
//!
//! A rational B-spline surface whose weights are all positive lies in the
//! convex hull of its controls, and so does every face on it. So two spline
//! faces with no common vertex cannot meet when a plane strictly separates
//! their controls.
//!
//! Two faces that share one edge, and no vertex off it, meet only along that
//! edge when a plane through the edge's two vertices has every control of one
//! face strictly on one side and every control of the other strictly on the
//! other side, except controls on the closed segment between the vertices.
//! Each hull then meets the plane only inside that segment, so the faces share
//! no point off the segment. The shared edge joins the two vertices inside
//! that segment, so it is the segment, and the faces meet only along it.
//!
//! The plane is searched exactly among finitely many candidates. Controls with
//! no common point are tried against the separating-axis candidates: normals
//! of planes through three controls of one face, and cross products of a
//! direction between two controls of each face. For a shared edge, the
//! candidates are normals of planes through the edge and one control, summed
//! in pairs, and each control's direction across the edge. A search that finds
//! no plane proves nothing, and the pair is left to the surface intersection.

use catallax_numbers::{Decider, Rational, Sign};
use std::cmp::Ordering;

use crate::curved::NotVerified;
use crate::{Face, Surface, Vec3};

/// The most controls a face may have for the search. Larger nets are left
/// to the surface intersection, so the search stays small.
const MAX_CONTROLS: usize = 9;

/// Whether two faces of one solid are proved to meet nowhere, or only along
/// their one shared edge. `Ok(false)` proves nothing.
pub(crate) fn apart(a: &Face, b: &Face, dec: &mut Decider) -> Result<bool, NotVerified> {
    let (Some(first), Some(second)) = (controls(a), controls(b)) else {
        return Ok(false);
    };
    let mut shared: Vec<&Vec3> = Vec::new();
    for p in vertices(a) {
        for q in vertices(b) {
            if same(p, q, dec)? && !contains(&shared, p, dec)? {
                shared.push(p);
            }
        }
    }
    match shared[..] {
        [] => separated(&first, &second, dec),
        [p, q] if shares_edge(a, b, p, q, dec)? => hinged(&first, &second, p, q, dec),
        _ => Ok(false),
    }
}

/// A spline face's controls, when its weights are all positive and the net
/// is small enough to search.
fn controls(face: &Face) -> Option<Vec<Vec3>> {
    let Surface::BSpline(spline) = &face.surface else {
        return None;
    };
    let geometry = &spline.geometry;
    if spline.self_intersect == Some(true)
        || geometry.controls().len() > MAX_CONTROLS
        || geometry
            .weights()
            .iter()
            .any(|w| w.sign() != Sign::Positive)
    {
        return None;
    }
    Some(geometry.controls().to_vec())
}

/// Every vertex point on a face's loops.
fn vertices(face: &Face) -> impl Iterator<Item = &Vec3> {
    face.loops.iter().flat_map(|lp| {
        lp.edges
            .iter()
            .flat_map(|e| [&e.start.point, &e.end.point])
            .chain(lp.vertex.iter().map(|v| &v.point))
    })
}

fn contains(points: &[&Vec3], p: &Vec3, dec: &mut Decider) -> Result<bool, NotVerified> {
    for q in points {
        if same(q, p, dec)? {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Whether one edge of both faces runs between `p` and `q`.
fn shares_edge(
    a: &Face,
    b: &Face,
    p: &Vec3,
    q: &Vec3,
    dec: &mut Decider,
) -> Result<bool, NotVerified> {
    for ea in a.loops.iter().flat_map(|lp| &lp.edges) {
        let joins = (same(&ea.start.point, p, dec)? && same(&ea.end.point, q, dec)?)
            || (same(&ea.start.point, q, dec)? && same(&ea.end.point, p, dec)?);
        if joins
            && b.loops
                .iter()
                .flat_map(|lp| &lp.edges)
                .any(|eb| eb.edge == ea.edge)
        {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Whether a plane strictly separates two sets of controls.
fn separated(a: &[Vec3], b: &[Vec3], dec: &mut Decider) -> Result<bool, NotVerified> {
    let mut normals = Vec::new();
    for points in [a, b] {
        for i in 0..points.len() {
            for j in i + 1..points.len() {
                for k in j + 1..points.len() {
                    let u = sub(&points[j], &points[i], dec)?;
                    let v = sub(&points[k], &points[i], dec)?;
                    normals.push(cross(&u, &v, dec)?);
                }
            }
        }
    }
    let directions = |points: &[Vec3], dec: &mut Decider| -> Result<Vec<Vec3>, NotVerified> {
        let mut out = Vec::new();
        for i in 0..points.len() {
            for j in i + 1..points.len() {
                out.push(sub(&points[j], &points[i], dec)?);
            }
        }
        Ok(out)
    };
    let (along_a, along_b) = (directions(a, dec)?, directions(b, dec)?);
    for u in &along_a {
        for v in &along_b {
            normals.push(cross(u, v, dec)?);
        }
    }
    for normal in &normals {
        dec.check_deadline()?;
        if normal.iter().all(Rational::is_zero) {
            continue;
        }
        let (a_low, a_high) = extent(normal, a, dec)?;
        let (b_low, b_high) = extent(normal, b, dec)?;
        if compare(&a_high, &b_low, dec)?.is_lt() || compare(&b_high, &a_low, dec)?.is_lt() {
            return Ok(true);
        }
    }
    Ok(false)
}

/// The least and greatest projection of the points on a normal.
fn extent(
    normal: &Vec3,
    points: &[Vec3],
    dec: &mut Decider,
) -> Result<(Rational, Rational), NotVerified> {
    let mut values = Vec::with_capacity(points.len());
    for p in points {
        values.push(dot(normal, p, dec)?);
    }
    let mut low = values[0].clone();
    let mut high = values[0].clone();
    for value in &values[1..] {
        if compare(value, &low, dec)?.is_lt() {
            low = value.clone();
        }
        if compare(value, &high, dec)?.is_gt() {
            high = value.clone();
        }
    }
    Ok((low, high))
}

/// Whether a plane through `p` and `q` has the controls of `a` strictly on
/// one side and those of `b` strictly on the other, except controls on the
/// closed segment from `p` to `q`.
fn hinged(
    a: &[Vec3],
    b: &[Vec3],
    p: &Vec3,
    q: &Vec3,
    dec: &mut Decider,
) -> Result<bool, NotVerified> {
    let d = sub(q, p, dec)?;
    let length = dot(&d, &d, dec)?;
    // Each control off the segment: its face, its offset from `p` and the
    // normal of the plane through the edge's line and the control.
    let mut off: Vec<(bool, Vec3, Vec3)> = Vec::new();
    for (first, points) in [(true, a), (false, b)] {
        for x in points {
            let v = sub(x, p, dec)?;
            let normal = cross(&d, &v, dec)?;
            if normal.iter().all(Rational::is_zero) {
                // On the edge's line: it must lie on the closed segment.
                let t = dot(&v, &d, dec)?;
                if t.sign() == Sign::Negative || compare(&t, &length, dec)?.is_gt() {
                    return Ok(false);
                }
                continue;
            }
            off.push((first, v, normal));
        }
    }
    let mut candidates = Vec::new();
    for (_, v, _) in &off {
        // `v` without its part along the edge, scaled by the edge's length.
        let along = dot(v, &d, dec)?;
        let mut across = [Rational::zero(), Rational::zero(), Rational::zero()];
        for i in 0..3 {
            let scaled = mul(&v[i], &length, dec)?;
            let part = mul(&d[i], &along, dec)?;
            across[i] = minus(&scaled, &part, dec)?;
        }
        candidates.push(across);
    }
    for i in 0..off.len() {
        for j in i + 1..off.len() {
            let (u, w) = (&off[i].2, &off[j].2);
            let mut sum = [Rational::zero(), Rational::zero(), Rational::zero()];
            let mut difference = sum.clone();
            for k in 0..3 {
                sum[k] = add(&u[k], &w[k], dec)?;
                difference[k] = minus(&u[k], &w[k], dec)?;
            }
            candidates.push(sum);
            candidates.push(difference);
        }
    }
    'candidates: for normal in &candidates {
        dec.check_deadline()?;
        // The side of the first face, once a control has shown it.
        let mut side: Option<Sign> = None;
        for (first, v, _) in &off {
            let sign = dot(normal, v, dec)?.sign();
            if sign == Sign::Zero {
                continue 'candidates;
            }
            let wanted = if *first { sign } else { flip(sign) };
            match side {
                None => side = Some(wanted),
                Some(s) if s == wanted => {}
                Some(_) => continue 'candidates,
            }
        }
        return Ok(true);
    }
    Ok(false)
}

fn flip(sign: Sign) -> Sign {
    match sign {
        Sign::Positive => Sign::Negative,
        Sign::Negative => Sign::Positive,
        Sign::Zero => Sign::Zero,
    }
}

// Each rational operation reserves its cost first.
fn mul(a: &Rational, b: &Rational, dec: &mut Decider) -> Result<Rational, NotVerified> {
    dec.reserve_exact_arithmetic(a, b, true)?;
    Ok(a * b)
}

fn add(a: &Rational, b: &Rational, dec: &mut Decider) -> Result<Rational, NotVerified> {
    dec.reserve_exact_arithmetic(a, b, false)?;
    Ok(a + b)
}

fn minus(a: &Rational, b: &Rational, dec: &mut Decider) -> Result<Rational, NotVerified> {
    dec.reserve_exact_arithmetic(a, b, false)?;
    Ok(a - b)
}

fn compare(a: &Rational, b: &Rational, dec: &mut Decider) -> Result<Ordering, NotVerified> {
    // Fraction ordering can cross-multiply twice.
    dec.reserve_exact_arithmetic(a, b, true)?;
    dec.reserve_exact_arithmetic(a, b, true)?;
    Ok(a.cmp(b))
}

fn same(a: &Vec3, b: &Vec3, dec: &mut Decider) -> Result<bool, NotVerified> {
    for i in 0..3 {
        dec.reserve_exact_arithmetic(&a[i], &b[i], false)?;
        if a[i] != b[i] {
            return Ok(false);
        }
    }
    Ok(true)
}

fn sub(a: &Vec3, b: &Vec3, dec: &mut Decider) -> Result<Vec3, NotVerified> {
    Ok([
        minus(&a[0], &b[0], dec)?,
        minus(&a[1], &b[1], dec)?,
        minus(&a[2], &b[2], dec)?,
    ])
}

fn dot(a: &Vec3, b: &Vec3, dec: &mut Decider) -> Result<Rational, NotVerified> {
    let mut result = Rational::zero();
    for i in 0..3 {
        let product = mul(&a[i], &b[i], dec)?;
        result = add(&result, &product, dec)?;
    }
    Ok(result)
}

fn cross(a: &Vec3, b: &Vec3, dec: &mut Decider) -> Result<Vec3, NotVerified> {
    let term = |i: usize, j: usize, dec: &mut Decider| -> Result<Rational, NotVerified> {
        let left = mul(&a[i], &b[j], dec)?;
        let right = mul(&a[j], &b[i], dec)?;
        minus(&left, &right, dec)
    };
    Ok([term(1, 2, dec)?, term(2, 0, dec)?, term(0, 1, dec)?])
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(x: i64, y: i64, z: i64) -> Vec3 {
        [x, y, z].map(Rational::from)
    }

    /// The M27 loft's first two side faces: bilinear nets sharing the edge
    /// from `(1, -1, 0)` to `(1, 0, 2)`.
    fn loft_sides() -> (Vec<Vec3>, Vec<Vec3>) {
        (
            vec![p(-1, -1, 0), p(0, -1, 2), p(1, -1, 0), p(1, 0, 2)],
            vec![p(1, -1, 0), p(1, 0, 2), p(1, 1, 0), p(0, 2, 2)],
        )
    }

    #[test]
    fn adjacent_loft_sides_meet_only_along_their_edge() {
        let (a, b) = loft_sides();
        let mut dec = Decider::default();
        assert_eq!(
            hinged(&a, &b, &p(1, -1, 0), &p(1, 0, 2), &mut dec),
            Ok(true)
        );
    }

    #[test]
    fn a_control_across_the_edge_proves_nothing() {
        let (a, mut b) = loft_sides();
        // The second face folds back over the first.
        b[3] = p(-1, -2, 2);
        let mut dec = Decider::default();
        assert_eq!(
            hinged(&a, &b, &p(1, -1, 0), &p(1, 0, 2), &mut dec),
            Ok(false)
        );
    }

    #[test]
    fn a_control_on_the_line_past_the_edge_proves_nothing() {
        let (mut a, b) = loft_sides();
        a.push(p(1, 1, 4));
        let mut dec = Decider::default();
        assert_eq!(
            hinged(&a, &b, &p(1, -1, 0), &p(1, 0, 2), &mut dec),
            Ok(false)
        );
    }

    #[test]
    fn opposite_loft_sides_are_separated() {
        let first = vec![p(1, -1, 0), p(1, 0, 2), p(1, 1, 0), p(0, 2, 2)];
        let third = vec![p(-1, 1, 0), p(-2, 0, 2), p(-1, -1, 0), p(0, -1, 2)];
        let mut dec = Decider::default();
        assert_eq!(separated(&first, &third, &mut dec), Ok(true));
    }

    #[test]
    fn crossing_nets_are_not_separated() {
        let first = vec![p(0, 0, 0), p(2, 0, 0), p(0, 2, 0), p(2, 2, 0)];
        let second = vec![p(1, 1, -1), p(1, 1, 1), p(3, 3, -1), p(3, 3, 1)];
        let mut dec = Decider::default();
        assert_eq!(separated(&first, &second, &mut dec), Ok(false));
    }

    #[test]
    fn a_small_allowance_stops_the_search() {
        let (a, b) = loft_sides();
        let (result, work) = catallax_numbers::with_operation_budget(10, || {
            hinged(&a, &b, &p(1, -1, 0), &p(1, 0, 2), &mut Decider::default())
        });
        assert!(work.exhausted);
        let refusal = result.unwrap_err();
        assert_eq!(
            (refusal.code, refusal.detail.as_str()),
            (
                "numbers.operation-budget",
                "exact geometric proof arithmetic stopped"
            )
        );
    }
}
