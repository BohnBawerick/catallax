//! Complete transverse line crossings of a finite polynomial extrusion.

use super::*;
use catallax_numbers::RealPoly;

fn unresolved(detail: &str) -> Refusal {
    Refusal::new("intersect.parametric-surface", detail)
}

fn in_unit(value: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(exact_sign(value, dec)? != Sign::Negative
        && exact_sign(&(&Real::one() - value), dec)? != Sign::Negative)
}

impl BezierExtrusion {
    /// Every transverse meeting of a whole line with the finite patch, in line order.
    /// A tangent meeting or a contained generator returns a named refusal.
    pub fn line_hits(
        &self,
        origin: &Point3,
        direction: &Vector3,
        dec: &mut Decider,
    ) -> Result<Vec<Point3>, Refusal> {
        math::validate_point(origin, dec)?;
        math::validate_vector(direction, dec)?;
        math::positive(&direction.norm_squared(), "shapes.zero-direction", dec)?;
        let axis = self.direction();
        let axis_sq = axis.norm_squared();
        let shadow = direction.cross(axis);
        let curve = self.lower_edge();
        if exact_sign(&shadow.norm_squared(), dec)? == Sign::Zero {
            let height = &origin.minus(&self.controls()[0]).dot(axis) / &axis_sq;
            let projected = origin.plus(&axis.scaled(&(-height)));
            return if curve.unique_parameter_of(&projected, dec)?.is_some() {
                Err(unresolved(
                    "the line contains a generator of the finite Bezier patch",
                ))
            } else {
                Ok(Vec::new())
            };
        }
        let along = direction.minus(&axis.scaled(&(&direction.dot(axis) / &axis_sq)));
        let along_sq = along.norm_squared();
        let coordinates = curve.query_span_polynomials(curve.degree(), dec)?;
        let targets = [&origin.x, &origin.y, &origin.z];
        let normal = [&shadow.x, &shadow.y, &shadow.z];
        let mut equation = RealPoly::constant(1, Real::zero());
        for i in 0..3 {
            equation =
                equation + (&coordinates[i] - &coordinates[3].scale(targets[i])).scale(normal[i]);
        }
        if equation.is_zero(dec)? {
            return Err(unresolved(
                "the line is coplanar with an entire Bezier generator family",
            ));
        }
        let mut hits: Vec<(Real, Point3)> = Vec::new();
        for parameter in equation.real_roots(dec)? {
            if !in_unit(&parameter, dec)? {
                continue;
            }
            let base = curve.point_at(&parameter, dec)?;
            let station = &base.minus(origin).dot(&along) / &along_sq;
            let point = origin.plus(&direction.scaled(&station));
            let height = &point.minus(&base).dot(axis) / &axis_sq;
            if !in_unit(&height, dec)? {
                continue;
            }
            let tangent = curve.first_at(&parameter, dec)?;
            if exact_sign(&tangent.dot(&shadow), dec)? == Sign::Zero {
                return Err(unresolved("the line is tangent to the finite Bezier patch"));
            }
            let mut position = hits.len();
            let mut duplicate = false;
            for (i, (previous, _)) in hits.iter().enumerate() {
                match exact_sign(&(&station - previous), dec)? {
                    Sign::Negative => {
                        position = i;
                        break;
                    }
                    Sign::Zero => {
                        duplicate = true;
                        break;
                    }
                    Sign::Positive => {}
                }
            }
            if !duplicate {
                hits.insert(position, (station, point));
            }
        }
        Ok(hits.into_iter().map(|(_, point)| point).collect())
    }
}
