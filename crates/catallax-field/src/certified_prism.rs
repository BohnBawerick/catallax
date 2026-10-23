//! Exact-target queries tied to an independently checked spline-body certificate.

use catallax_features::{CertifiedBuilt, PrismTarget};
use catallax_numbers::{Ball, Decider, Interval, Rational, Real, RealPoly, Refusal, Sign};
use catallax_shapes::{EllipticCylinder, RPoint3, Surface, Vector3};
use catallax_topology::FaceId;

use crate::{Direction, DistanceBound, FaceQuery, FaceState, Location, Query, Value};
use crate::geometry::{enclose, sign};
use crate::metric::{PRECISION, Scalar};

/// Nominal field of a certified offset prism. Construction checks the actual
/// body and its complete certificate, and borrows that immutable result.
/// Queries concern the original exact target, not the fitted control net;
/// the body's recorded uniform bound relates the two geometries separately.
/// This consumer does not give authority to a diagnostic mark or raw spline.
#[derive(Debug)]
pub struct CertifiedPrismField<'a> {
    built: &'a CertifiedBuilt,
    frame: Option<Frame<'a>>,
}

#[derive(Debug)]
struct Frame<'a> {
    original: &'a EllipticCylinder,
    axis: Vector3,
    major: Vector3,
    minor: Vector3,
    raw_minor: Vector3,
    reference_squared: Real,
}

struct Cross {
    signed: Real,
    sign: Sign,
    normal: Option<Vector3>,
}

struct Sample {
    query: Query,
    cross: Sign,
    axial: Sign,
    cap: usize,
    side_normal: Option<Vector3>,
}

impl<'a> CertifiedPrismField<'a> {
    /// Checks the geometry-bound record before retaining the immutable body.
    pub fn new(built: &'a CertifiedBuilt, dec: &mut Decider) -> Result<Self, Refusal> {
        built.verify(dec)?;
        let target = built.result_certificate().target();
        let frame = if target.is_empty() { None } else {
            if target.definition().side() == catallax_shapes::OffsetSide::Outward
                && sign(target.definition().distance(), dec)? != Sign::Zero
            {
                return Err(Refusal::new("field.certified-prism-branch", "this consumer requires an inward regular offset or its exact zero-distance source"));
            }
            let Surface::EllipticCylinder(original) = target.definition().original() else {
                unreachable!("validated elliptic-offset target");
            };
            let axis = original.axis.scaled(&(&Real::one() / &original.axis.norm_squared().sqrt()));
            let reference_squared = original.major_dir.norm_squared();
            let inverse_length = &Real::one() / &reference_squared.sqrt();
            let raw_minor = axis.cross(&original.major_dir);
            Some(Frame {
                original, axis, major: original.major_dir.scaled(&inverse_length),
                minor: raw_minor.scaled(&inverse_length), raw_minor, reference_squared,
            })
        };
        Ok(Self { built, frame })
    }

    /// Authoritative nominal classification and true Euclidean distance at an
    /// exact rational point, under the caller's unchanged numeric policy.
    pub fn proven_checked(&self, point: &RPoint3, dec: &mut Decider) -> Result<Query, Refusal> {
        Ok(self.sample(point, dec)?.query)
    }

    /// Membership and outward direction of one corresponding nominal face.
    /// Handles are the actual certificate's lower cap, upper cap and side IDs.
    /// A fitted mesh point need not lie on the exact target; no tolerance silently
    /// changes that point into an exact nominal boundary point.
    pub fn face_query(&self, face: FaceId, point: &RPoint3, dec: &mut Decider) -> Result<FaceQuery, Refusal> {
        let certificate = self.built.result_certificate();
        let Some(faces) = certificate.faces() else { return Ok(not_on_face()); };
        let Some(index) = faces.iter().position(|candidate| *candidate == face) else { return Ok(not_on_face()); };
        let sample = self.sample(point, dec)?;
        let normal = if index == 2 && sample.cross == Sign::Zero && sample.axial != Sign::Positive {
            sample.side_normal
        } else if index < 2 && sample.axial == Sign::Zero && sample.cross != Sign::Positive && sample.cap == index {
            let axis = &self.frame.as_ref().unwrap().axis;
            Some(if index == 0 { axis.scaled(&Real::from_i64(-1)) } else { axis.clone() })
        } else { None };
        match normal {
            Some(normal) => Ok(FaceQuery { state: FaceState::OnFace, direction: crate::unit_direction(normal, dec)? }),
            None => Ok(not_on_face()),
        }
    }

    fn sample(&self, point: &RPoint3, dec: &mut Decider) -> Result<Sample, Refusal> {
        dec.check_deadline()?;
        let Some(frame) = &self.frame else {
            let infinity = Interval::point(f64::INFINITY).unwrap();
            return Ok(Sample {
                query: Query { location: Location::Outside, value: Value::TrueDistance(infinity),
                    distance: DistanceBound { millimetres: infinity }, direction: Direction::Undefined },
                cross: Sign::Positive, axial: Sign::Positive, cap: 0, side_normal: None,
            });
        };
        let target = self.built.result_certificate().target();
        let displacement = point.to_real().minus(&frame.original.origin);
        let z = displacement.dot(&frame.axis);
        let axial = target.axial();
        let midpoint = (axial[0] + axial[1]) / Real::from_i64(2);
        let cap = usize::from(sign(&(&z - &midpoint), dec)? != Sign::Negative);
        let axial_distance = if cap == 0 { axial[0] - &z } else { &z - axial[1] };
        let axial_sign = sign(&axial_distance, dec)?;
        let cross = frame.cross(target, &displacement, axial_sign != Sign::Positive, dec)?;
        let location = if cross.sign == Sign::Positive || axial_sign == Sign::Positive {
            Location::Outside
        } else if cross.sign == Sign::Zero || axial_sign == Sign::Zero {
            Location::Boundary
        } else { Location::Inside };
        let a = enclose(&cross.signed, PRECISION, dec)?;
        let b = enclose(&axial_distance, PRECISION, dec)?;
        let zero = Ball::from_rational(&Rational::zero(), PRECISION);
        let outside2 = a.max(&zero).square() + b.max(&zero).square();
        let outside = outside2.sqrt().ok_or_else(|| Refusal::new("field.distance-domain", "nonnegative product-distance square failed enclosure"))?;
        let signed = outside + a.max(&b).min(&zero);
        let value = if location == Location::Boundary { Interval::point(0.0).unwrap() } else { signed.to_interval() };
        let direction = if location != Location::Boundary { Direction::Undefined }
        else if cross.sign == Sign::Zero && axial_sign == Sign::Zero { Direction::Ambiguous }
        else if cross.sign == Sign::Zero { crate::unit_direction(cross.normal.as_ref().unwrap().clone(), dec)? }
        else {
            crate::unit_direction(if cap == 0 { frame.axis.scaled(&Real::from_i64(-1)) } else { frame.axis.clone() }, dec)?
        };
        Ok(Sample {
            query: Query { location, value: Value::TrueDistance(value),
                distance: DistanceBound { millimetres: crate::unsigned(value) }, direction },
            cross: cross.sign, axial: axial_sign, cap, side_normal: cross.normal,
        })
    }
}

fn not_on_face() -> FaceQuery {
    FaceQuery { state: FaceState::NotOnFace, direction: Direction::Undefined }
}

impl Frame<'_> {
    fn normal(&self, u: &Real, v: &Real) -> Vector3 {
        let a = self.major.scaled(u);
        let b = self.minor.scaled(v);
        Vector3::new(&a.x + &b.x, &a.y + &b.y, &a.z + &b.z)
    }

    fn cross(&self, target: &PrismTarget, delta: &Vector3, need_normal: bool, dec: &mut Decider) -> Result<Cross, Refusal> {
        let a2 = &self.original.major_sq;
        let b2 = &self.original.minor_sq;
        // Square projections in the original rational directions, rather than
        // eliminate dependent normalising radicals from posed coordinates.
        let raw_x = delta.dot(&self.original.major_dir);
        let raw_y = delta.dot(&self.raw_minor);
        let x2 = &(&raw_x * &raw_x) / &self.reference_squared;
        let y2 = &(&raw_y * &raw_y) / &self.reference_squared;
        let d = target.definition().distance();
        let gap = &(&(&x2 / a2) + &(&y2 / b2)) - &Real::one();
        let original_sign = sign(&gap, dec)?;
        if original_sign == Sign::Zero {
            let class = sign(d, dec)?;
            let normal = if need_normal && class == Sign::Zero {
                let u = delta.dot(&self.major);
                let v = delta.dot(&self.minor);
                Some(self.normal(&(&u / a2), &(&v / b2)))
            } else { None };
            return Ok(Cross { signed: d.clone(), sign: class, normal });
        }
        let mut normal = None;
        let mut normal_coordinates = None;
        let squared = if sign(&raw_x, dec)? == Sign::Zero {
            let distance = y2.sqrt() - b2.sqrt();
            let signed = &distance + d;
            let class = sign(&signed, dec)?;
            if need_normal && class == Sign::Zero {
                normal = Some(self.minor.scaled(&Real::from_i64(if sign(&raw_y, dec)? == Sign::Negative { -1 } else { 1 })));
            }
            return Ok(Cross { signed, sign: class, normal });
        } else if sign(&raw_y, dec)? == Sign::Zero {
            let difference = a2 - b2;
            let threshold = &(&difference * &difference) / a2;
            if sign(&(&x2 - &threshold), dec)? != Sign::Positive {
                b2 * &(&Real::one() - &(&x2 / &difference))
            } else {
                let signed = &(&x2.sqrt() - &a2.sqrt()) + d;
                let class = sign(&signed, dec)?;
                if need_normal && class == Sign::Zero {
                    normal = Some(self.major.scaled(&Real::from_i64(if sign(&raw_x, dec)? == Sign::Negative { -1 } else { 1 })));
                }
                return Ok(Cross { signed, sign: class, normal });
            }
        } else {
            let variable = RealPoly::var(1, 0);
            let da = &variable + &RealPoly::constant(1, a2.clone());
            let db = &variable + &RealPoly::constant(1, b2.clone());
            let da2 = da.pow_u32(2);
            let db2 = db.pow_u32(2);
            let polynomial = &(&da2 * &db2) - &(&db2.scale(&(a2 * &x2)) + &da2.scale(&(b2 * &y2)));
            // For nonzero minor projection, the nearest normal's multiplier is
            // the unique root above -b^2: its normal equation decreases strictly
            // from +infinity to zero there. Other stationary normals are not the
            // nearest point. Axis singularities are the two branches above.
            let mut nearest = None;
            for candidate in polynomial.real_roots(dec)? {
                if sign(&(&candidate + b2), dec)? == Sign::Positive { nearest = Some(candidate); break; }
            }
            let lambda = nearest.ok_or_else(|| Refusal::new("field.ellipse-distance-root", "nearest normal multiplier was not certified"))?;
            let du = a2 + &lambda;
            let dv = b2 + &lambda;
            let u = &lambda / &du;
            let v = &lambda / &dv;
            let squared = &(&(&u * &u) * &x2) + &(&(&v * &v) * &y2);
            if need_normal && original_sign == Sign::Negative { normal_coordinates = Some([du, dv]); }
            squared
        };
        let d2 = d * d;
        let class = if original_sign == Sign::Positive { Sign::Positive }
        else { match sign(&(&squared - &d2), dec)? {
            Sign::Positive => Sign::Negative, Sign::Zero => Sign::Zero, Sign::Negative => Sign::Positive,
        }};
        if class == Sign::Zero && let Some([du, dv]) = normal_coordinates {
            let u = delta.dot(&self.major);
            let v = delta.dot(&self.minor);
            normal = Some(self.normal(&(&u / &du), &(&v / &dv)));
        }
        // A point on an inward regular offset is within the original normal
        // reach. The original nearest point is unique there, so its normal is
        // also the offset's outward normal. A major-axis medial point cannot
        // be on this branch because its distance is at least the reach.
        let signed = if class == Sign::Zero { Real::zero() }
        else if original_sign == Sign::Negative { d - &squared.sqrt() }
        else { &squared.sqrt() + d };
        Ok(Cross { signed, sign: class, normal })
    }
}
