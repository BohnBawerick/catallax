//! Complete oriented boundary moments of a rectangle intersected with disks.
use super::{Circle, Counter, q};
use crate::geometry;
use catallax_numbers::{Ball, Decider, Rational, Real, Refusal, Sign};
use std::cmp::Ordering;

type Point = [Real; 2];
const POWERS: [[usize; 2]; 6] = [[0, 0], [1, 0], [0, 1], [2, 0], [0, 2], [1, 1]];
fn r(value: &Rational) -> Real {
    Real::from(value.clone())
}
fn fraction(n: i64, d: i64) -> Rational {
    q(n).checked_div(&q(d)).unwrap()
}
fn sign(value: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    geometry::sign(value, dec)
}
fn cmp(a: &Real, b: &Real, dec: &mut Decider) -> Result<Ordering, Refusal> {
    Ok(match sign(&(a - b), dec)? {
        Sign::Negative => Ordering::Less,
        Sign::Zero => Ordering::Equal,
        Sign::Positive => Ordering::Greater,
    })
}
fn offset(point: &Point, circle: &Circle) -> Point {
    [
        &point[0] - &r(&circle.centre[0]),
        &point[1] - &r(&circle.centre[1]),
    ]
}
fn cross(a: &Point, b: &Point) -> Real {
    &a[0] * &b[1] - &a[1] * &b[0]
}
fn dot(a: &Point, b: &Point) -> Real {
    &a[0] * &b[0] + &a[1] * &b[1]
}
fn half(v: &Point, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(match sign(&v[1], dec)? {
        Sign::Positive => false,
        Sign::Negative => true,
        Sign::Zero => sign(&v[0], dec)? == Sign::Negative,
    })
}
fn angle_cmp(a: &Point, b: &Point, c: &Circle, dec: &mut Decider) -> Result<Ordering, Refusal> {
    let (a, b) = (offset(a, c), offset(b, c));
    let (ha, hb) = (half(&a, dec)?, half(&b, dec)?);
    if ha != hb {
        return Ok(ha.cmp(&hb));
    }
    Ok(match sign(&cross(&a, &b), dec)? {
        Sign::Positive => Ordering::Less,
        Sign::Negative => Ordering::Greater,
        Sign::Zero => Ordering::Equal,
    })
}
fn sorted_points(
    mut points: Vec<Point>,
    circle: &Circle,
    dec: &mut Decider,
) -> Result<Vec<Point>, Refusal> {
    for i in 1..points.len() {
        let mut j = i;
        while j > 0 && angle_cmp(&points[j - 1], &points[j], circle, dec)? == Ordering::Greater {
            points.swap(j - 1, j);
            j -= 1;
        }
    }
    let mut unique: Vec<Point> = Vec::new();
    for point in points {
        if unique.last().is_none()
            || angle_cmp(unique.last().unwrap(), &point, circle, dec)? != Ordering::Equal
        {
            unique.push(point);
        }
    }
    Ok(unique)
}
fn roots(circle: &Circle, axis: usize, coordinate: &Rational) -> Vec<Point> {
    let d = coordinate - &circle.centre[axis];
    let square = &circle.radius * &circle.radius - &d * &d;
    if square < q(0) {
        return Vec::new();
    }
    let root = r(&square).sqrt();
    [1, -1]
        .into_iter()
        .map(|s| {
            let mut p = circle.centre.each_ref().map(r);
            p[axis] = r(coordinate);
            p[1 - axis] = &p[1 - axis] + &(&Real::from_i64(s) * &root);
            p
        })
        .collect()
}
fn meetings(a: &Circle, b: &Circle) -> Vec<Point> {
    let d = [&b.centre[0] - &a.centre[0], &b.centre[1] - &a.centre[1]];
    let square = &d[0] * &d[0] + &d[1] * &d[1];
    if square.is_zero()
        || square > (&a.radius + &b.radius).pow(2).unwrap()
        || square < (&a.radius - &b.radius).pow(2).unwrap()
    {
        return Vec::new();
    }
    let t = (&(&a.radius * &a.radius) - &(&b.radius * &b.radius) + square.clone())
        .checked_div(&square.mul_pow2(1))
        .unwrap();
    let k = (&a.radius * &a.radius).checked_div(&square).unwrap() - &t * &t;
    let root = r(&k).sqrt();
    let base = [&a.centre[0] + &(&d[0] * &t), &a.centre[1] + &(&d[1] * &t)];
    let normal = [-&d[1], d[0].clone()];
    [1, -1]
        .into_iter()
        .map(|s| {
            std::array::from_fn(|i| r(&base[i]) + (&r(&normal[i]) * &root) * Real::from_i64(s))
        })
        .collect()
}
fn distance(point: &Point, circle: &Circle) -> Real {
    let d = offset(point, circle);
    dot(&d, &d) - r(&(&circle.radius * &circle.radius))
}
fn disks_inside(
    point: &Point,
    circles: &[Circle],
    except: Option<usize>,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let mut inside = true;
    for (i, c) in circles.iter().enumerate() {
        if except == Some(i) {
            continue;
        }
        match sign(&distance(point, c), dec)? {
            Sign::Positive => inside = false,
            Sign::Negative => {}
            Sign::Zero => {
                return Err(Refusal::new(
                    "field.mass-prismatic-cell",
                    format!("open boundary-cell point lies on another carrier: {point:?}"),
                ));
            }
        }
    }
    Ok(inside)
}
fn arc_inside(
    point: &Point,
    circles: &[Circle],
    except: usize,
    rectangle: &[[Rational; 2]; 2],
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let mut inside = disks_inside(point, circles, Some(except), dec)?;
    for i in 0..2 {
        let lo = cmp(&point[i], &r(&rectangle[i][0]), dec)?;
        let hi = cmp(&point[i], &r(&rectangle[i][1]), dec)?;
        if lo == Ordering::Equal || hi == Ordering::Equal {
            return Err(Refusal::new(
                "field.mass-prismatic-cell",
                format!("open circular cell meets a rectangle carrier: {point:?}"),
            ));
        }
        inside &= lo == Ordering::Greater && hi == Ordering::Less;
    }
    Ok(inside)
}
fn midpoint(
    a: &Point,
    b: &Point,
    c: &Circle,
    full: bool,
    dec: &mut Decider,
) -> Result<Point, Refusal> {
    let (a, b) = (offset(a, c), offset(b, c));
    let v = if full {
        [-&a[0], -&a[1]]
    } else {
        let s = sign(&cross(&a, &b), dec)?;
        if s == Sign::Zero {
            [-&a[1], a[0].clone()]
        } else {
            let v = [&a[0] + &b[0], &a[1] + &b[1]];
            if s == Sign::Negative {
                [-&v[0], -&v[1]]
            } else {
                v
            }
        }
    };
    let norm = dot(&v, &v).sqrt();
    Ok(std::array::from_fn(|i| {
        r(&c.centre[i]) + r(&c.radius) * (&v[i] / &norm)
    }))
}
fn delta(
    a: &Point,
    b: &Point,
    c: &Circle,
    full: bool,
    bits: u32,
    dec: &mut Decider,
) -> Result<Ball, Refusal> {
    let pi = Ball::pi_checked(bits, dec)?;
    if full {
        return Ok(&pi * &Ball::from_rational(&q(2), bits));
    }
    let (a, b) = (offset(a, c), offset(b, c));
    let (y, x) = (cross(&a, &b), dot(&a, &b));
    let (sy, sx) = (sign(&y, dec)?, sign(&x, dec)?);
    if sy == Sign::Zero {
        return if sx == Sign::Negative {
            Ok(pi)
        } else {
            Err(Refusal::new(
                "field.mass-prismatic-cell",
                "duplicate endpoints of an open circular cell",
            ))
        };
    }
    if sx == Sign::Zero {
        return Ok(&pi
            * &Ball::from_rational(&fraction(if sy == Sign::Positive { 1 } else { 3 }, 2), bits));
    }
    let x = geometry::enclose(&x, bits, dec)?;
    let y = geometry::enclose(&y, bits, dec)?;
    let quotient = y.checked_div(&x).ok_or_else(|| {
        Refusal::new(
            "field.mass-prismatic-refine",
            format!(
                "arc-angle dot enclosure [{},{}] mm2 at {bits} bits needs refinement",
                x.lo().to_decimal_bounds(10).0,
                x.hi().to_decimal_bounds(10).1
            ),
        )
    })?;
    let angle = quotient.atan_checked(dec)?;
    Ok(if sx == Sign::Negative {
        &angle + &pi
    } else if sy == Sign::Negative {
        &angle + &(&pi * &Ball::from_rational(&q(2), bits))
    } else {
        angle
    })
}

fn primitive(c: &Real, s: &Real, cp: usize, sp: usize) -> (Rational, Real) {
    let cs = c * s;
    let fourth = &cs * &(&c.pow_u32(2) - &s.pow_u32(2));
    match (cp, sp) {
        (1, 0) => (q(0), s.clone()),
        (2, 0) => (fraction(1, 2), cs / Real::from_i64(2)),
        (3, 0) => (q(0), s - &(&s.pow_u32(3) / &Real::from_i64(3))),
        (4, 0) => (
            fraction(3, 8),
            cs / Real::from_i64(2) + fourth / Real::from_i64(8),
        ),
        (k, 1) => (
            q(0),
            -&c.pow_u32((k + 1) as u32) / Real::from_i64((k + 1) as i64),
        ),
        (1, 2) => (q(0), s.pow_u32(3) / Real::from_i64(3)),
        (2, 2) => (fraction(1, 8), -fourth / Real::from_i64(8)),
        _ => unreachable!("moments have degree at most two"),
    }
}
fn choose(n: usize, k: usize) -> i64 {
    const ROWS: [&[i64]; 4] = [&[1], &[1, 1], &[1, 2, 1], &[1, 3, 3, 1]];
    ROWS[n][k]
}
fn arc_moments(
    a: &Point,
    b: &Point,
    c: &Circle,
    full: bool,
    bits: u32,
    dec: &mut Decider,
) -> Result<[Ball; 6], Refusal> {
    let delta = delta(a, b, c, full, bits, dec)?;
    let radius = r(&c.radius);
    let a = offset(a, c).map(|v| v / radius.clone());
    let b = offset(b, c).map(|v| v / radius.clone());
    let mut trig = Vec::new();
    for pair in [
        (1, 0),
        (2, 0),
        (3, 0),
        (4, 0),
        (1, 1),
        (2, 1),
        (3, 1),
        (1, 2),
        (2, 2),
    ] {
        dec.charge_work(1)?;
        let (factor, fa) = primitive(&a[0], &a[1], pair.0, pair.1);
        let (_, fb) = primitive(&b[0], &b[1], pair.0, pair.1);
        let value = &geometry::enclose(&(fb - fa), bits, dec)?
            + &(&delta * &Ball::from_rational(&factor, bits));
        trig.push((pair, value));
    }
    Ok(std::array::from_fn(|index| {
        let [p, qexp] = POWERS[index];
        let mut sum = Ball::from_rational(&q(0), bits);
        for i in 0..=p + 1 {
            for j in 0..=qexp {
                let factor = (&q(choose(p + 1, i) * choose(qexp, j))
                    * &(&c.centre[0].pow((p + 1 - i) as i32).unwrap()
                        * &c.centre[1].pow((qexp - j) as i32).unwrap())
                    * c.radius.pow((i + j + 1) as i32).unwrap())
                .checked_div(&q((p + 1) as i64))
                .unwrap();
                if factor.is_zero() {
                    continue;
                }
                let value = &trig.iter().find(|(pair, _)| *pair == (i + 1, j)).unwrap().1;
                sum = &sum + &(value * &Ball::from_rational(&factor, bits));
            }
        }
        sum
    }))
}

pub(super) fn moments(
    circles: &[Circle],
    rectangle: &[[Rational; 2]; 2],
    bits: u32,
    counter: &mut Counter,
    dec: &mut Decider,
) -> Result<[Ball; 6], Refusal> {
    let mut points: Vec<Vec<Point>> = vec![Vec::new(); circles.len()];
    for (i, c) in circles.iter().enumerate() {
        for axis in 0..2 {
            for bound in &rectangle[axis] {
                points[i].extend(roots(c, axis, bound));
            }
        }
        for j in 0..i {
            counter.step(dec)?;
            let meeting = meetings(c, &circles[j]);
            points[i].extend(meeting.iter().cloned());
            points[j].extend(meeting);
        }
    }
    let mut result: [Ball; 6] = std::array::from_fn(|_| Ball::from_rational(&q(0), bits));
    for (i, c) in circles.iter().enumerate() {
        let points = sorted_points(std::mem::take(&mut points[i]), c, dec)?;
        if points.is_empty() {
            return Err(Refusal::new(
                "field.mass-prismatic-cell",
                "no vertices for a bounded disk intersection",
            ));
        }
        for k in 0..points.len() {
            counter.step(dec)?;
            let (a, b) = (&points[k], &points[(k + 1) % points.len()]);
            let full = points.len() == 1;
            if !arc_inside(&midpoint(a, b, c, full, dec)?, circles, i, rectangle, dec)? {
                continue;
            }
            let moments = arc_moments(a, b, c, full, bits, dec)?;
            for j in 0..6 {
                result[j] = &result[j] + &moments[j];
            }
        }
    }
    // Horizontal rectangle edges have dy=0 in Green's integral.
    for (side, x) in rectangle[0].iter().enumerate() {
        let mut nodes = vec![r(&rectangle[1][0]), r(&rectangle[1][1])];
        for c in circles {
            for point in roots(c, 0, x) {
                if cmp(&point[1], &r(&rectangle[1][0]), dec)? != Ordering::Less
                    && cmp(&point[1], &r(&rectangle[1][1]), dec)? != Ordering::Greater
                {
                    nodes.push(point[1].clone());
                }
            }
        }
        for i in 1..nodes.len() {
            let mut j = i;
            while j > 0 && cmp(&nodes[j - 1], &nodes[j], dec)? == Ordering::Greater {
                nodes.swap(j - 1, j);
                j -= 1;
            }
        }
        let mut unique = Vec::new();
        for node in nodes {
            if unique.last().is_none()
                || cmp(unique.last().unwrap(), &node, dec)? != Ordering::Equal
            {
                unique.push(node);
            }
        }
        for pair in unique.windows(2) {
            counter.step(dec)?;
            let point = [r(x), (&pair[0] + &pair[1]) / Real::from_i64(2)];
            if !disks_inside(&point, circles, None, dec)? {
                continue;
            }
            for index in 0..6 {
                let [p, qexp] = POWERS[index];
                let value = (r(x).pow_u32((p + 1) as u32)
                    * (pair[1].pow_u32((qexp + 1) as u32) - pair[0].pow_u32((qexp + 1) as u32)))
                    / Real::from_i64(((p + 1) * (qexp + 1)) as i64);
                let value = if side == 0 { -value } else { value };
                result[index] = &result[index] + &geometry::enclose(&value, bits, dec)?;
            }
        }
    }
    Ok(result)
}
