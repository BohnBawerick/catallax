use super::{equal, section_data, sign};
use crate::BezierSectionPair;
use catallax_numbers::{Decider, Rational, Real, Refusal, Sign};
use catallax_shapes::{BSplineCurve, CurveId, GeomStore, Point3, Vector3};

fn unsupported(detail: &'static str) -> Refusal {
    Refusal::new("intersect.bezier-curve-pair", detail)
}

#[derive(Clone)]
struct Span {
    controls: Vec<Point3>,
    domain: [Rational; 2],
}
impl Span {
    fn whole(curve: &BSplineCurve) -> Self {
        Self {
            controls: (0..curve.control_count())
                .map(|i| curve.control(i).clone())
                .collect(),
            domain: [Rational::zero(), Rational::one()],
        }
    }
    fn endpoints(&self) -> [&Point3; 2] {
        [&self.controls[0], self.controls.last().unwrap()]
    }
    fn split(&self) -> [Self; 2] {
        let mut row = self.controls.clone();
        let mut left = vec![row[0].clone()];
        let mut right = vec![row.last().unwrap().clone()];
        let half = Real::from(Rational::parse("1/2").unwrap());
        for count in (1..row.len()).rev() {
            for i in 0..count {
                row[i] = row[i].plus(&row[i + 1].minus(&row[i]).scaled(&half));
            }
            left.push(row[0].clone());
            right.push(row[count - 1].clone());
        }
        right.reverse();
        let mid = (&self.domain[0] + &self.domain[1]) / Rational::from_i64(2);
        [
            Self {
                controls: left,
                domain: [self.domain[0].clone(), mid.clone()],
            },
            Self {
                controls: right,
                domain: [mid, self.domain[1].clone()],
            },
        ]
    }
}

fn projection(point: &Point3, axis: &Vector3) -> Real {
    Vector3::new(point.x.clone(), point.y.clone(), point.z.clone()).dot(axis)
}
fn bounds(span: &Span, axis: &Vector3, dec: &mut Decider) -> Result<(Real, Real), Refusal> {
    let mut lo = projection(&span.controls[0], axis);
    let mut hi = lo.clone();
    for point in &span.controls[1..] {
        let value = projection(point, axis);
        if sign(&(&value - &lo), dec)? == Sign::Negative {
            lo = value.clone();
        }
        if sign(&(&value - &hi), dec)? == Sign::Positive {
            hi = value;
        }
    }
    Ok((lo, hi))
}
fn allowed(point: &Point3, shared: &[Point3], dec: &mut Decider) -> Result<bool, Refusal> {
    for other in shared {
        if equal(point, other, dec)? {
            return Ok(true);
        }
    }
    Ok(false)
}
fn flat(span: &Span, axis: &Vector3, cut: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
    for point in &span.controls {
        if sign(&(&projection(point, axis) - cut), dec)? != Sign::Zero {
            return Ok(false);
        }
    }
    Ok(true)
}

// Positive Bernstein weights keep a nonconstant projection off its supporting
// plane in the open span. Only endpoint contacts remain on that plane.
fn contact(
    spans: [&Span; 2],
    curves: [&BSplineCurve; 2],
    axis: &Vector3,
    cut: &Real,
    shared: &[Point3],
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let is_flat = [
        flat(spans[0], axis, cut, dec)?,
        flat(spans[1], axis, cut, dec)?,
    ];
    if is_flat == [true, true] {
        return Ok(false);
    }
    if is_flat[0] || is_flat[1] {
        let f = usize::from(is_flat[1]);
        let curved = 1 - f;
        for point in spans[curved].endpoints() {
            if sign(&(&projection(point, axis) - cut), dec)? != Sign::Zero
                || allowed(point, shared, dec)?
            {
                continue;
            }
            if let Some(parameter) = curves[f].unique_parameter_of(point, dec)? {
                let value = parameter.value();
                if sign(&(value - &Real::from(spans[f].domain[0].clone())), dec)? != Sign::Negative
                    && sign(&(&Real::from(spans[f].domain[1].clone()) - value), dec)?
                        != Sign::Negative
                {
                    return Err(unsupported(
                        "finite Bezier sections meet away from shared original endpoints",
                    ));
                }
            }
        }
        return Ok(true);
    }
    for point in spans[0].endpoints() {
        if sign(&(&projection(point, axis) - cut), dec)? != Sign::Zero {
            continue;
        }
        for other in spans[1].endpoints() {
            if equal(point, other, dec)? && !allowed(point, shared, dec)? {
                return Err(unsupported(
                    "finite Bezier sections meet away from shared original endpoints",
                ));
            }
        }
    }
    Ok(true)
}

fn separated(
    spans: [&Span; 2],
    curves: [&BSplineCurve; 2],
    normals: [&Vector3; 2],
    shared: &[Point3],
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let mut axes = vec![
        Vector3::new(Real::one(), Real::zero(), Real::zero()),
        Vector3::new(Real::zero(), Real::one(), Real::zero()),
        Vector3::new(Real::zero(), Real::zero(), Real::one()),
        normals[0].clone(),
        normals[1].clone(),
    ];
    let edges: [Vec<Vector3>; 2] = spans.map(|span| {
        let mut result = Vec::new();
        for i in 0..span.controls.len() {
            for j in i + 1..span.controls.len() {
                result.push(span.controls[j].minus(&span.controls[i]));
            }
        }
        result
    });
    for edge in edges.iter().flatten() {
        for normal in normals {
            axes.push(normal.cross(edge));
        }
    }
    for a in &edges[0] {
        for b in &edges[1] {
            axes.push(a.cross(b));
        }
    }
    for axis in axes {
        if sign(&axis.norm_squared(), dec)? == Sign::Zero {
            continue;
        }
        let (a_lo, a_hi) = bounds(spans[0], &axis, dec)?;
        let (b_lo, b_hi) = bounds(spans[1], &axis, dec)?;
        let ab = sign(&(&a_hi - &b_lo), dec)?;
        let ba = sign(&(&b_hi - &a_lo), dec)?;
        if ab == Sign::Negative || ba == Sign::Negative {
            return Ok(true);
        }
        if ab == Sign::Zero && contact(spans, curves, &axis, &a_hi, shared, dec)?
            || ba == Sign::Zero && contact(spans, curves, &axis, &a_lo, shared, dec)?
        {
            return Ok(true);
        }
    }
    Ok(false)
}

pub(crate) fn meetings(
    geom: &GeomStore,
    a: CurveId,
    b: CurveId,
    dec: &mut Decider,
) -> Result<Option<BezierSectionPair>, Refusal> {
    let Some((first, normal_a)) = section_data(geom, a, dec)? else {
        return Ok(None);
    };
    let Some((second, normal_b)) = section_data(geom, b, dec)? else {
        return Ok(None);
    };
    prove(&first, &second, &normal_a, &normal_b, dec).map(Some)
}

pub(super) fn prove(
    first: &BSplineCurve,
    second: &BSplineCurve,
    normal_a: &Vector3,
    normal_b: &Vector3,
    dec: &mut Decider,
) -> Result<BezierSectionPair, Refusal> {
    let before = dec.stats().coincident_separation;
    let a = Span::whole(first);
    let b = Span::whole(second);
    if a.controls.len() == b.controls.len() {
        for reverse in [false, true] {
            let mut same = true;
            for i in 0..a.controls.len() {
                let j = if reverse { b.controls.len() - 1 - i } else { i };
                if !equal(&a.controls[i], &b.controls[j], dec)? {
                    same = false;
                    break;
                }
            }
            if same {
                return Ok(BezierSectionPair::Coincident);
            }
        }
    }
    let mut shared = Vec::new();
    for point in a.endpoints() {
        for other in b.endpoints() {
            if equal(point, other, dec)? && !allowed(point, &shared, dec)? {
                shared.push(point.clone());
            }
        }
    }
    let mut pending = vec![(a, b, 0)];
    let mut cells = 0;
    while let Some((a, b, depth)) = pending.pop() {
        dec.check_deadline()?;
        if separated(
            [&a, &b],
            [first, second],
            [normal_a, normal_b],
            &shared,
            dec,
        )? {
            continue;
        }
        cells += 1;
        if cells > 4096 || depth >= 48 {
            return Err(unsupported(
                "control hull subdivision did not prove complete finite Bezier meetings",
            ));
        }
        if depth % 2 == 0 {
            let [left, right] = a.split();
            pending.push((right, b.clone(), depth + 1));
            pending.push((left, b, depth + 1));
        } else {
            let [left, right] = b.split();
            pending.push((a.clone(), right, depth + 1));
            pending.push((a, left, depth + 1));
        }
    }
    if dec.stats().coincident_separation != before {
        return Err(unsupported(
            "finite Bezier meetings need exact incidence without a separation-bound zero",
        ));
    }
    Ok(BezierSectionPair::Points(shared))
}
