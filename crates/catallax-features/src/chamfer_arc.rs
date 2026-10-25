//! Finite circular domains, with no angular rounding or seam-dependent order.

use catallax_numbers::{Decider, Rational, Real, Refusal, Sign};
use catallax_shapes::{Circle, Curve, Plane, Point3, Vector3};

use super::{ChamferDefinition, bad, root, unit};
use crate::profile::sign;

/// A nonempty proper counterclockwise arc in its circle's natural direction.
#[derive(Clone, Debug)]
pub struct BoundedArc {
    /// Exact circular support.
    pub circle: Circle,
    /// Natural start and end, independent of a topological edge's sense.
    pub ends: [Point3; 2],
}

impl BoundedArc {
    /// Checks incidence and rejects a collapsed or whole-circle domain.
    pub fn new(circle: Circle, ends: [Point3; 2], dec: &mut Decider) -> Result<Self, Refusal> {
        if sign(&circle.normal.norm_squared(), dec)? == Sign::Zero
            || sign(&circle.ref_dir.norm_squared(), dec)? == Sign::Zero
            || sign(&circle.ref_dir.dot(&circle.normal), dec)? != Sign::Zero
            || sign(&circle.radius_sq, dec)? != Sign::Positive
            || sign(&ends[0].minus(&ends[1]).norm_squared(), dec)? == Sign::Zero
        {
            return Err(bad(
                "a bounded circular run needs distinct regular endpoints",
            ));
        }
        for p in &ends {
            let v = p.minus(&circle.center);
            if sign(&v.dot(&circle.normal), dec)? != Sign::Zero
                || sign(&(&v.norm_squared() - &circle.radius_sq), dec)? != Sign::Zero
            {
                return Err(bad("circular run endpoint misses its support"));
            }
        }
        Ok(Self { circle, ends })
    }

    /// Angular membership, including endpoints; the point may have another radius.
    pub fn contains_direction(&self, p: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        let [a, b] = self.ends.each_ref().map(|p| p.minus(&self.circle.center));
        let v = p.minus(&self.circle.center);
        if sign(&v.norm_squared(), dec)? == Sign::Zero {
            return Ok(false);
        }
        let cross = |a: &Vector3, b: &Vector3| a.cross(b).dot(&self.circle.normal);
        let av = sign(&cross(&a, &v), dec)?;
        let vb = sign(&cross(&v, &b), dec)?;
        Ok(match sign(&cross(&a, &b), dec)? {
            Sign::Positive => av != Sign::Negative && vb != Sign::Negative,
            Sign::Negative => av != Sign::Negative || vb != Sign::Negative,
            Sign::Zero => av != Sign::Negative,
        })
    }

    /// Exact interior direction, valid for minor, major and seam-crossing arcs.
    pub fn middle(&self, dec: &mut Decider) -> Result<Point3, Refusal> {
        let [a, b] = self.ends.each_ref().map(|p| p.minus(&self.circle.center));
        let mut v = a.plus(&b);
        if sign(&v.norm_squared(), dec)? == Sign::Zero {
            v = unit(&self.circle.normal).cross(&a);
        } else if sign(&a.cross(&b).dot(&self.circle.normal), dec)? == Sign::Negative {
            v = v.scaled(&Real::from_i64(-1));
        }
        Ok(self
            .circle
            .center
            .plus(&v.scaled(&root(&(&self.circle.radius_sq / &v.norm_squared())))))
    }

    /// A rational point strictly inside a rationally parameterised circle arc.
    /// Approximation proposes a parameter only; exact membership accepts it.
    pub fn rational_sample(&self, dec: &mut Decider) -> Result<Point3, Refusal> {
        let middle = self.middle(dec)?.minus(&self.circle.center);
        let x = unit(&self.circle.ref_dir);
        let y = unit(&self.circle.normal).cross(&x);
        let radius = root(&self.circle.radius_sq);
        let cosine = &middle.dot(&x) / &radius;
        let sine = &middle.dot(&y) / &radius;
        // Choose the stereographic chart whose denominator is bounded from zero.
        let flip = if sign(&cosine, dec)? == Sign::Negative {
            -1
        } else {
            1
        };
        let k = Real::from_i64(flip);
        let parameter = &sine / &(&Real::one() + &(&k * &cosine));
        for bits in [32, 64, 128, 256, 512] {
            let ball = dec.enclose_budgeted_at(&parameter, bits)?;
            let t = Real::from((ball.lo() + ball.hi()) / Rational::from_i64(2));
            let denominator = Real::one() + t.square();
            let c = &k * &(&(Real::one() - t.square()) / &denominator);
            let s = &(&Real::from_i64(2) * &t) / &denominator;
            let p = self
                .circle
                .center
                .plus(&x.scaled(&(&radius * &c)).plus(&y.scaled(&(&radius * &s))));
            if p.as_rational().is_some()
                && self.contains_direction(&p, dec)?
                && sign(&p.minus(&self.ends[0]).norm_squared(), dec)? != Sign::Zero
                && sign(&p.minus(&self.ends[1]).norm_squared(), dec)? != Sign::Zero
            {
                return Ok(p);
            }
        }
        Err(bad(
            "circular material probes need a rational interior point",
        ))
    }

    /// Continues the selected transverse circle/plane root onto a coaxial
    /// contact circle. The derivative sign fixes the branch exactly.
    pub fn contact_at(
        &self,
        circle: &Circle,
        end: usize,
        plane: &Plane,
        dec: &mut Decider,
    ) -> Result<Point3, Refusal> {
        if sign(&plane.normal.norm_squared(), dec)? == Sign::Zero
            || sign(&circle.radius_sq, dec)? != Sign::Positive
            || sign(
                &circle.normal.cross(&self.circle.normal).norm_squared(),
                dec,
            )? != Sign::Zero
            || sign(
                &circle
                    .center
                    .minus(&self.circle.center)
                    .cross(&self.circle.normal)
                    .norm_squared(),
                dec,
            )? != Sign::Zero
            || sign(&plane.normal.dot(&self.circle.normal), dec)? != Sign::Zero
            || sign(&self.ends[end].minus(&plane.origin).dot(&plane.normal), dec)? != Sign::Zero
        {
            return Err(bad(
                "finite circular ends need transverse planes parallel to the cylinder axis",
            ));
        }
        let u = unit(&plane.normal);
        let v = unit(&self.circle.normal).cross(&u);
        let side = sign(&self.ends[end].minus(&self.circle.center).dot(&v), dec)?;
        if side == Sign::Zero {
            return Err(bad("a tangent circular end needs continuation"));
        }
        let offset = plane.origin.minus(&circle.center).dot(&u);
        let discriminant = &circle.radius_sq - &offset.square();
        if sign(&discriminant, dec)? != Sign::Positive {
            return Err(Refusal::new(
                "chamfer.unsupported-overflow",
                "an end plane consumes or touches a contact circle",
            ));
        }
        let along =
            root(&discriminant) * Real::from_i64(if side == Sign::Positive { 1 } else { -1 });
        Ok(circle
            .center
            .plus(&u.scaled(&offset))
            .plus(&v.scaled(&along)))
    }

    /// Proves that the two end cuts remain transverse and do not meet anywhere
    /// in the swept annular band. Returns oriented end halfspaces and their
    /// intersection/union rule, certified on the complete original circle.
    pub fn trim(
        &self,
        definition: &ChamferDefinition,
        planes: &[Plane; 2],
        dec: &mut Decider,
    ) -> Result<([Plane; 2], bool), Refusal> {
        let mut lo = self.circle.radius_sq.clone();
        let mut hi = lo.clone();
        for contact in &definition.contacts {
            let Curve::Circle(c) = contact else {
                return Err(bad("circular end needs circular contacts"));
            };
            for (end, plane) in planes.iter().enumerate() {
                self.contact_at(c, end, plane, dec)?;
            }
            if sign(&(&c.radius_sq - &lo), dec)? == Sign::Negative {
                lo = c.radius_sq.clone();
            }
            if sign(&(&c.radius_sq - &hi), dec)? == Sign::Positive {
                hi = c.radius_sq.clone();
            }
        }
        for plane in planes {
            let offset = plane
                .origin
                .minus(&self.circle.center)
                .dot(&unit(&plane.normal));
            if sign(&(&lo - &offset.square()), dec)? != Sign::Positive {
                return Err(Refusal::new(
                    "chamfer.unsupported-overflow",
                    "an end cut loses transversality inside the affected band",
                ));
            }
        }
        let [a, b] = planes.each_ref().map(|p| &p.normal);
        let aa = a.norm_squared();
        let bb = b.norm_squared();
        let ab = a.dot(b);
        let determinant = &(&aa * &bb) - &ab.square();
        if sign(&determinant, dec)? != Sign::Zero {
            let da = planes[0].origin.minus(&self.circle.center).dot(a);
            let db = planes[1].origin.minus(&self.circle.center).dot(b);
            let crossing = a.scaled(&(&(&da * &bb) - &(&db * &ab)));
            // Solve the two linear end equations in the circle plane.
            let crossing = crossing.scaled(&(&Real::one() / &determinant)).plus(
                &b.scaled(&(&(&db * &aa) - &(&da * &ab)))
                    .scaled(&(&Real::one() / &determinant)),
            );
            let radius_sq = crossing.norm_squared();
            if sign(&(&radius_sq - &lo), dec)? != Sign::Negative
                && sign(&(&radius_sq - &hi), dec)? != Sign::Positive
            {
                return Err(Refusal::new(
                    "chamfer.unsupported-overflow",
                    "the two end cuts cross inside the affected band",
                ));
            }
        }
        let middle = self.middle(dec)?;
        let mut gates = planes.clone();
        for plane in &mut gates {
            match sign(&middle.minus(&plane.origin).dot(&plane.normal), dec)? {
                Sign::Positive => plane.normal = plane.normal.scaled(&Real::from_i64(-1)),
                Sign::Negative => (),
                Sign::Zero => {
                    return Err(bad("end halfspaces do not select a regular arc interior"));
                }
            }
        }
        let mut points = self.ends.to_vec();
        for (end, plane) in planes.iter().enumerate() {
            let p = self.contact_at(&self.circle, end, plane, dec)?;
            let u = unit(&plane.normal);
            let foot = self
                .circle
                .center
                .plus(&u.scaled(&plane.origin.minus(&self.circle.center).dot(&u)));
            points.push(foot.plus(&foot.minus(&p)));
        }
        let points = ordered_points(&self.circle, points, dec)?;
        for major in [false, true] {
            let mut matches = true;
            for i in 0..points.len() {
                let arc = Self::new(
                    self.circle.clone(),
                    [points[i].clone(), points[(i + 1) % points.len()].clone()],
                    dec,
                )?;
                let p = arc.middle(dec)?;
                let mut inside = [false; 2];
                for (i, plane) in gates.iter().enumerate() {
                    inside[i] =
                        sign(&p.minus(&plane.origin).dot(&plane.normal), dec)? != Sign::Positive;
                }
                let clipped = if major {
                    inside[0] || inside[1]
                } else {
                    inside[0] && inside[1]
                };
                matches &= clipped == self.contains_direction(&p, dec)?;
            }
            if matches {
                return Ok((gates, major));
            }
        }
        Err(bad(
            "the end halfspaces do not bound one circular input component",
        ))
    }
}

/// Sorts and deduplicates exact points around a circle, with its reference as seam.
pub fn ordered_points(
    circle: &Circle,
    points: Vec<Point3>,
    dec: &mut Decider,
) -> Result<Vec<Point3>, Refusal> {
    let x = &circle.ref_dir;
    let y = circle.normal.cross(x);
    let half = |v: &Vector3, dec: &mut Decider| -> Result<bool, Refusal> {
        let sy = sign(&v.dot(&y), dec)?;
        Ok(sy == Sign::Negative || (sy == Sign::Zero && sign(&v.dot(x), dec)? == Sign::Negative))
    };
    let mut ordered: Vec<Point3> = Vec::new();
    for point in points {
        let v = point.minus(&circle.center);
        let mut index = ordered.len();
        let mut duplicate = false;
        for (i, old) in ordered.iter().enumerate() {
            if sign(&point.minus(old).norm_squared(), dec)? == Sign::Zero {
                duplicate = true;
                break;
            }
            let w = old.minus(&circle.center);
            let (a, b) = (half(&v, dec)?, half(&w, dec)?);
            if (!a && b)
                || (a == b && sign(&v.cross(&w).dot(&circle.normal), dec)? == Sign::Positive)
            {
                index = i;
                break;
            }
        }
        if !duplicate {
            ordered.insert(index, point);
        }
    }
    Ok(ordered)
}
