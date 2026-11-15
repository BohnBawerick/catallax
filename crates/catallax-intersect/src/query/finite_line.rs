//! Exact equality of complete finite straight spline point sets.

use crate::spline_arithmetic::{arithmetic, difference, sign};
use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_shapes::{BSplineCurve, Point3};

const CODE: &str = "intersect.finite-line-equivalence";

fn coordinates(point: &Point3) -> [&Real; 3] {
    [&point.x, &point.y, &point.z]
}

fn equal(a: &Point3, b: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
    for (a, b) in coordinates(a).into_iter().zip(coordinates(b)) {
        if sign(&difference(a, b, dec)?, dec, CODE)? != Sign::Zero {
            return Ok(false);
        }
    }
    Ok(true)
}

struct Segment {
    ends: [Point3; 2],
    direction: [Real; 3],
    axis: usize,
    sense: Sign,
}

impl Segment {
    fn collinear(&self, point: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        let at = coordinates(point);
        let origin = coordinates(&self.ends[0]);
        let delta = difference(at[self.axis], origin[self.axis], dec)?;
        for other in (0..3).filter(|axis| *axis != self.axis) {
            let offset = difference(at[other], origin[other], dec)?;
            let a = arithmetic(&self.direction[self.axis], &offset, true, dec)?;
            let b = arithmetic(&self.direction[other], &delta, true, dec)?;
            if sign(&difference(&a, &b, dec)?, dec, CODE)? != Sign::Zero {
                return Ok(false);
            }
        }
        Ok(true)
    }

    fn contains(&self, point: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        let value = coordinates(point)[self.axis];
        let lo = coordinates(&self.ends[0])[self.axis];
        let hi = coordinates(&self.ends[1])[self.axis];
        for delta in [difference(value, lo, dec)?, difference(hi, value, dec)?] {
            let sign = sign(&delta, dec, CODE)?;
            if sign != Sign::Zero && sign != self.sense {
                return Ok(false);
            }
        }
        Ok(true)
    }
}

fn segment(curve: &BSplineCurve, dec: &mut Decider) -> Result<Option<Segment>, Refusal> {
    let ends = [
        curve.control(0).clone(),
        curve.control(curve.control_count() - 1).clone(),
    ];
    let a = coordinates(&ends[0]);
    let b = coordinates(&ends[1]);
    let direction = [
        difference(b[0], a[0], dec)?,
        difference(b[1], a[1], dec)?,
        difference(b[2], a[2], dec)?,
    ];
    let mut selected = None;
    for (axis, value) in direction.iter().enumerate() {
        let sense = sign(value, dec, CODE)?;
        if sense != Sign::Zero {
            selected = Some((axis, sense));
            break;
        }
    }
    let Some((axis, sense)) = selected else {
        return Ok(None);
    };
    let segment = Segment {
        ends,
        direction,
        axis,
        sense,
    };
    for control in curve.controls() {
        if !segment.collinear(control, dec)? || !segment.contains(control, dec)? {
            return Ok(None);
        }
    }
    // Admission guarantees positive weights and C0 continuity. The convex
    // hull stays in the segment; the endpoint values cover it continuously.
    Ok(Some(segment))
}

pub(super) fn same(
    a: &BSplineCurve,
    b: &BSplineCurve,
    dec: &mut Decider,
) -> Result<Option<bool>, Refusal> {
    let first = segment(a, dec)?;
    let second = segment(b, dec)?;
    let first_is_segment = first.is_some();
    match (first, second) {
        (Some(a), Some(b)) => Ok(Some(
            (equal(&a.ends[0], &b.ends[0], dec)? && equal(&a.ends[1], &b.ends[1], dec)?)
                || (equal(&a.ends[0], &b.ends[1], dec)? && equal(&a.ends[1], &b.ends[0], dec)?),
        )),
        (Some(line), None) | (None, Some(line)) => {
            let other = if first_is_segment { b } else { a };
            // The admitted B-spline basis is linearly independent. A control
            // off the line proves the rational curve is not wholly on it.
            for control in other.controls() {
                if !line.collinear(control, dec)? {
                    return Ok(Some(false));
                }
            }
            for end in [other.control(0), other.control(other.control_count() - 1)] {
                if !line.contains(end, dec)? {
                    return Ok(Some(false));
                }
            }
            Ok(None)
        }
        (None, None) => Ok(None),
    }
}
