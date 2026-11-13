//! Definition-native queries for a finite oblique cylindrical draft cell.

use super::*;
use catallax_features::draft::oblique::ObliqueCell;

#[derive(Clone, Debug)]
pub(super) struct FastOblique {
    origin: [f64; 3],
    axis: [f64; 3],
    across: [f64; 3],
    side: [f64; 3],
    heights: [f64; 2],
    amplitude: f64,
    tangent: f64,
    lean: [f64; 2],
}
impl FastOblique {
    pub(super) fn new(cell: &ObliqueCell) -> Self {
        let s = &cell.support;
        Self {
            origin: point(s.origin()),
            axis: vector(s.axis()),
            across: vector(s.across()),
            side: vector(s.side()),
            heights: s.heights().each_ref().map(Real::to_f64),
            amplitude: s.amplitude(&Real::zero()).to_f64(),
            tangent: s.tangent().to_f64(),
            lean: s.lean().each_ref().map(Real::to_f64),
        }
    }
    pub(super) fn wgsl(&self) -> String {
        use crate::unproved_wgsl::{number, vector};
        format!(
            "let d=p-{}; let z=dot(d,{}); let x=dot(d,{}); let y=dot(d,{}); \
             let rho=length(vec2<f32>(x,y)); let a={}-{}*z; var radial=-a; \
             if(rho!=0.0) {{radial=rho-a-({}*x+{}*y)/rho;}} \
             return max(max(radial,{}-z),z-{});",
            vector(&self.origin),
            vector(&self.axis),
            vector(&self.across),
            vector(&self.side),
            number(self.amplitude),
            number(self.tangent),
            number(self.lean[0]),
            number(self.lean[1]),
            number(self.heights[0]),
            number(self.heights[1]),
        )
    }
    pub(super) fn residual(&self, p: [f64; 3]) -> f64 {
        let d = sub(p, self.origin);
        let z = dot(d, self.axis);
        let x = dot(d, self.across);
        let y = dot(d, self.side);
        let rho = x.hypot(y);
        let a = self.amplitude - self.tangent * z;
        let radial = if rho == 0.0 {
            -a
        } else {
            rho - a - (self.lean[0] * x + self.lean[1] * y) / rho
        };
        radial.max(self.heights[0] - z).max(z - self.heights[1])
    }
}

pub(super) fn location(
    cell: &ObliqueCell,
    p: &Point3,
    dec: &mut Decider,
) -> Result<Location, Refusal> {
    Ok(match cell.locate(p, dec)? {
        Sign::Negative => Location::Inside,
        Sign::Zero => Location::Boundary,
        Sign::Positive => Location::Outside,
    })
}

pub(super) fn distance(
    cell: &ObliqueCell,
    p: &Point3,
    l: Location,
    dec: &mut Decider,
) -> Result<Ball, Refusal> {
    // Every chosen point on the cap ring is an actual boundary point. The
    // nearest boundary distance is nonnegative and at most its distance to
    // that point. This is a geometric enclosure, not a radial residual or a
    // claim that this convenient point is the nearest point.
    if l == Location::Boundary {
        return Ok(Ball::from_rational(&Rational::zero(), PRECISION));
    }
    let boundary = cell
        .support
        .point([&Real::one(), &Real::zero()], &cell.support.heights()[0]);
    let upper = geometry::enclose(&p.minus(&boundary).norm_squared().sqrt(), PRECISION, dec)?.hi();
    let ball = Ball::from_rational(&Rational::zero(), PRECISION)
        .hull(&Ball::from_rational(&upper, PRECISION));
    Ok(if l == Location::Inside { -ball } else { ball })
}

pub(super) fn on_face(
    cell: &ObliqueCell,
    origin: &Origin,
    p: &Point3,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    if cell.locate(p, dec)? != Sign::Zero {
        return Ok(false);
    }
    if cell.origins.contains(origin) {
        return Ok(cell.support.side_at(p, dec)? == Sign::Zero);
    }
    for (plane, names) in &cell.caps {
        if names.contains(origin) {
            return Ok(
                geometry::sign(&p.minus(&plane.origin).dot(&plane.normal), dec)? == Sign::Zero,
            );
        }
    }
    Ok(false)
}

pub(super) fn normal(
    cell: &ObliqueCell,
    origin: &Origin,
    p: &Point3,
    dec: &mut Decider,
) -> Result<Option<Vector3>, Refusal> {
    if !on_face(cell, origin, p, dec)? {
        return Ok(None);
    }
    let z = cell.support.coordinates(p)[2].clone();
    if cell.origins.contains(origin) {
        if geometry::sign(&(&z - &cell.support.heights()[0]), dec)? == Sign::Zero
            || geometry::sign(&(&z - &cell.support.heights()[1]), dec)? == Sign::Zero
        {
            return Ok(None);
        }
        return cell.support.normal_at(p, dec).map(Some);
    }
    if cell.support.side_at(p, dec)? == Sign::Zero {
        return Ok(None);
    }
    Ok(cell
        .caps
        .iter()
        .find(|(_, names)| names.contains(origin))
        .map(|(p, _)| p.normal.clone()))
}
