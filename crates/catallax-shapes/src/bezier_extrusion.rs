//! Exact regular normal extrusions of monotone polynomial Bezier spans.

use crate::{
    BSplineCurve, BSplineSurface, Isometry, Point3, Pose, RPoint3, RVector3, SplineBoundary,
    Vector3, math,
};
use catallax_numbers::{Decider, Outcome, Rational, Real, Refusal, Sign, StrictSign, ZeroProof};
use std::sync::Arc;

mod inverse;
mod line;

#[derive(Debug)]
struct Data {
    controls: Vec<Point3>,
    direction: Vector3,
    axis_sq: Real,
    chord: Vector3,
    chord_sq: Real,
    projection: Vec<Rational>,
    carrier: BSplineSurface,
}

/// A finite quadratic or cubic Bezier span extruded normally to its plane.
/// Monotone projection proves injectivity; endpoint tangents prove regularity
/// where the interior projection derivative may tend to zero.
#[derive(Clone, Debug)]
pub struct BezierExtrusion {
    data: Arc<Data>,
}

impl BezierExtrusion {
    /// Checks coplanarity, normal extrusion, regularity and embedding exactly.
    /// A span whose chord projection cannot prove embedding is refused.
    pub fn new(
        controls: Vec<RPoint3>,
        direction: RVector3,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        dec.check_deadline()?;
        if !(3..=4).contains(&controls.len()) {
            return Err(Refusal::new(
                "shapes.bezier-extrusion-controls",
                "a quadratic or cubic span needs three or four controls",
            ));
        }
        if direction.is_zero() {
            return Err(Refusal::new(
                "shapes.bezier-extrusion-direction",
                "the extrusion vector must be nonzero",
            ));
        }
        let first = &controls[0];
        for point in &controls[1..] {
            if !point.minus(first).dot(&direction).is_zero() {
                return Err(Refusal::new(
                    "shapes.bezier-extrusion-plane",
                    "the controls must lie in a plane normal to the extrusion vector",
                ));
            }
        }
        let last = controls.last().expect("checked control count");
        let chord = last.minus(first);
        if chord.is_zero()
            || controls[1].minus(first).is_zero()
            || last.minus(&controls[controls.len() - 2]).is_zero()
        {
            return Err(Refusal::new(
                "shapes.bezier-extrusion-singular",
                "span endpoints and their first tangents must be distinct and nonzero",
            ));
        }
        // Nonnegative Bernstein coefficients and one positive coefficient make
        // the chord projection derivative positive throughout the open span.
        let slopes: Vec<_> = controls
            .windows(2)
            .map(|pair| pair[1].minus(&pair[0]).dot(&chord))
            .collect();
        if slopes.iter().any(|s| s.sign() == Sign::Negative) || slopes.iter().all(Rational::is_zero)
        {
            return Err(Refusal::new(
                "shapes.bezier-extrusion-embedding",
                "the span has no proved monotone chord projection",
            ));
        }
        let count = controls.len();
        let knots = std::iter::repeat_n(Rational::zero(), count)
            .chain(std::iter::repeat_n(Rational::one(), count))
            .collect();
        let net = controls
            .iter()
            .flat_map(|point| [point.to_real(), point.plus(&direction).to_real()])
            .collect();
        let carrier = BSplineSurface::new(
            [count - 1, 1],
            [
                knots,
                vec![
                    Rational::zero(),
                    Rational::zero(),
                    Rational::one(),
                    Rational::one(),
                ],
            ],
            [count, 2],
            net,
            vec![Real::one(); count * 2],
            dec,
        )?;
        Ok(Self {
            data: Arc::new(Data {
                controls: controls.iter().map(RPoint3::to_real).collect(),
                direction: direction.to_real(),
                axis_sq: direction.norm_squared().into(),
                chord: chord.to_real(),
                chord_sq: chord.norm_squared().into(),
                projection: controls
                    .iter()
                    .map(|point| point.minus(first).dot(&chord) / chord.norm_squared())
                    .collect(),
                carrier,
            }),
        })
    }

    /// The exact polynomial ruled carrier on the unit parameter square.
    pub fn carrier(&self) -> &BSplineSurface {
        &self.data.carrier
    }

    /// The retained lower-boundary controls in parameter order.
    pub fn controls(&self) -> &[Point3] {
        &self.data.controls
    }

    /// The displacement from the lower boundary to the upper boundary.
    pub fn direction(&self) -> &Vector3 {
        &self.data.direction
    }

    /// A boundary view sharing the lower control row with the ruled carrier.
    pub fn lower_edge(&self) -> BSplineCurve {
        BSplineCurve::boundary(self.carrier(), SplineBoundary::VMin)
    }

    /// A boundary view sharing the upper control row with the ruled carrier.
    pub fn upper_edge(&self) -> BSplineCurve {
        BSplineCurve::boundary(self.carrier(), SplineBoundary::VMax)
    }

    fn parameters_at(
        &self,
        point: &Point3,
        dec: &mut Decider,
    ) -> Result<Option<(Real, Real)>, Refusal> {
        math::validate_point(point, dec)?;
        let direction = self.direction();
        let height = &point.minus(&self.controls()[0]).dot(direction) / &self.data.axis_sq;
        if exact_sign(&height, dec)? == Sign::Negative
            || exact_sign(&(&Real::one() - &height), dec)? == Sign::Negative
        {
            return Ok(None);
        }
        let projected = point.plus(&direction.scaled(&(-&height)));
        let curve = self.lower_edge();
        Ok(curve
            .unique_parameter_of(&projected, dec)?
            .map(|parameter| (parameter.value().clone(), height)))
    }

    /// The natural chart normal at a point on the finite regular patch.
    /// Points off the patch return a domain refusal.
    pub fn normal_at(&self, point: &Point3, dec: &mut Decider) -> Result<Vector3, Refusal> {
        let Some((parameter, _)) = self.parameters_at(point, dec)? else {
            return Err(Refusal::new(
                "shapes.bezier-extrusion-domain",
                "the normal point is not on the finite Bezier extrusion",
            ));
        };
        Ok(self
            .lower_edge()
            .first_at(&parameter, dec)?
            .cross(self.direction()))
    }

    /// Natural oriented side inside the finite chord-projection chart.
    /// Off-patch points on a singular endpoint chart keep a domain refusal.
    pub fn side_at(&self, point: &Point3, dec: &mut Decider) -> Result<Sign, Refusal> {
        use catallax_numbers::RealPoly;
        math::validate_point(point, dec)?;
        let offset = point.minus(&self.controls()[0]);
        let height = &offset.dot(self.direction()) / &self.data.axis_sq;
        let along = offset.dot(&self.data.chord);
        if !in_unit(&height, dec)?
            || exact_sign(&along, dec)? == Sign::Negative
            || exact_sign(&(&self.data.chord_sq - &along), dec)? == Sign::Negative
        {
            return Err(query_domain(
                "the side point lies outside the finite Bezier chart",
            ));
        }
        let curve = self.lower_edge();
        let coordinates = curve.query_span_polynomials(curve.degree(), dec)?;
        let targets = [&point.x, &point.y, &point.z];
        let chord = [&self.data.chord.x, &self.data.chord.y, &self.data.chord.z];
        let mut equation = RealPoly::constant(1, Real::zero());
        for i in 0..3 {
            equation =
                equation + (&coordinates[i] - &coordinates[3].scale(targets[i])).scale(chord[i]);
        }
        let mut found: Option<Real> = None;
        for parameter in equation.real_roots(dec)? {
            if !in_unit(&parameter, dec)? {
                continue;
            }
            if let Some(previous) = &found {
                if exact_sign(&(&parameter - previous), dec)? != Sign::Zero {
                    return Err(query_domain("the finite graph parameter is ambiguous"));
                }
            } else {
                found = Some(parameter);
            }
        }
        let parameter =
            found.ok_or_else(|| query_domain("no finite graph parameter was proved"))?;
        let base = curve.point_at(&parameter, dec)?;
        let across = self.direction().cross(&self.data.chord);
        let residual = point.minus(&base).dot(&across);
        let sign = exact_sign(&residual, dec)?;
        if sign == Sign::Zero {
            return Ok(Sign::Zero);
        }
        let tangent = curve.first_at(&parameter, dec)?;
        if exact_sign(&tangent.dot(&self.data.chord), dec)? != Sign::Positive {
            return Err(query_domain(
                "the endpoint graph chart has no off-surface side proof",
            ));
        }
        Ok(match sign {
            Sign::Positive => Sign::Negative,
            Sign::Negative => Sign::Positive,
            Sign::Zero => unreachable!("incidence handled above"),
        })
    }

    /// Exact rigid transport of the carrier and its retained embedding proof.
    pub fn transformed(&self, pose: &Pose) -> Self {
        self.mapped(&Isometry::from(pose))
    }

    pub(crate) fn mapped(&self, map: &Isometry) -> Self {
        Self {
            data: Arc::new(Data {
                controls: self
                    .controls()
                    .iter()
                    .map(|point| map.apply_point(point))
                    .collect(),
                direction: map.apply_vector(self.direction()),
                axis_sq: self.data.axis_sq.clone(),
                chord: map.apply_vector(&self.data.chord),
                chord_sq: self.data.chord_sq.clone(),
                projection: self.data.projection.clone(),
                carrier: self.carrier().mapped(map),
            }),
        }
    }

    /// Positive retained-definition identity for geometric query caches.
    pub fn shares_definition_nodes(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.data, &other.data)
    }
}

fn exact_sign(value: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    match dec.sign(value) {
        Outcome::Proven(StrictSign::Positive) => Ok(Sign::Positive),
        Outcome::Proven(StrictSign::Negative) => Ok(Sign::Negative),
        Outcome::ProvenCoincident(proof) if proof.proof != ZeroProof::SeparationBound => {
            Ok(Sign::Zero)
        }
        Outcome::ProvenCoincident(_) => Err(Refusal::new(
            "numbers.geometric-gap-unproven",
            "a finite Bezier query needs an exact geometric decision",
        )),
        Outcome::NotExact(error) => Err(Refusal::new(error.code, error.detail)),
        Outcome::Refused(error) => Err(error),
    }
}

fn query_domain(detail: &str) -> Refusal {
    Refusal::new("shapes.bezier-extrusion-domain", detail)
}

fn in_unit(value: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(exact_sign(value, dec)? != Sign::Negative
        && exact_sign(&(&Real::one() - value), dec)? != Sign::Negative)
}
