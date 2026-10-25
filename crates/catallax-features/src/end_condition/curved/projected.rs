//! Whole-rectangle quadratic branch certificates, never sample acceptance.

use super::*;
use catallax_numbers::{Ball, RealPoly};

/// The exact affine chart of the profile plane and raw authored ray direction.
#[derive(Clone, Debug)]
pub struct Frame {
    /// World-space chart origin.
    pub origin: Point3,
    /// Chart u direction.
    pub u: Vector3,
    /// Chart v direction.
    pub v: Vector3,
    /// Unnormalised authored sweep direction.
    pub direction: Vector3,
}

impl Frame {
    /// An exact point in the profile plane.
    pub fn point(&self, q: &[Real; 2]) -> Point3 {
        self.origin
            .plus(&self.u.scaled(&q[0]))
            .plus(&self.v.scaled(&q[1]))
    }
}

pub(super) struct Certificate {
    pub branch: Branch,
    pub midpoint: Plane,
    pub extent: Rational,
}

pub(super) fn coefficients(surface: &Surface, frame: &Frame) -> Result<[RealPoly; 3], Refusal> {
    // Variables are u,v,t. Substitution retains coefficients and root selectors.
    let vars: [RealPoly; 3] = std::array::from_fn(|i| RealPoly::var(3, i));
    let o = [&frame.origin.x, &frame.origin.y, &frame.origin.z];
    let u = [&frame.u.x, &frame.u.y, &frame.u.z];
    let v = [&frame.v.x, &frame.v.y, &frame.v.z];
    let w = [&frame.direction.x, &frame.direction.y, &frame.direction.z];
    let coordinates: [RealPoly; 3] = std::array::from_fn(|i| {
        RealPoly::constant(3, o[i].clone())
            + vars[0].scale(u[i])
            + vars[1].scale(v[i])
            + vars[2].scale(w[i])
    });
    let mut polynomial = RealPoly::zero(3);
    let implicit = surface
        .implicit()
        .ok_or_else(|| coverage("limiting carrier has no implicit polynomial"))?;
    for (m, coefficient) in implicit.terms() {
        let mut term = RealPoly::constant(3, coefficient);
        for (i, coordinate) in coordinates.iter().enumerate() {
            term = term * coordinate.pow_u32(m.exponents()[i]);
        }
        polynomial = polynomial + term;
    }
    let mut coefficients: [Vec<_>; 3] = std::array::from_fn(|_| Vec::new());
    for (m, coefficient) in polynomial.terms() {
        let degree = m.exponents()[2] as usize;
        if degree > 2 {
            return Err(coverage(
                "nonquadratic ray equation requires a certified first-hit branch, not a truncated quadratic",
            ));
        }
        coefficients[degree].push((m.exponents()[..2].to_vec(), coefficient));
    }
    Ok(coefficients.map(|terms| RealPoly::from_terms(2, terms)))
}

pub(super) fn certify(
    scene: &Scene,
    frame: &Frame,
    bounds: &[[Rational; 2]; 2],
    dec: &mut Decider,
) -> Result<Certificate, Refusal> {
    dec.check_deadline()?;
    if sign(&frame.direction.norm_squared(), dec)? == Sign::Zero {
        return Err(Refusal::new(
            "extrude.zero-direction",
            "extrusion direction is zero",
        ));
    }
    let n = frame.u.cross(&frame.v);
    let projection = n.dot(&frame.direction);
    if sign(&projection, dec)? == Sign::Zero {
        return Err(Refusal::new(
            "extrude.zero-thickness",
            "ray lies in a degenerate profile plane",
        ));
    }
    for [lo, hi] in bounds {
        if lo >= hi {
            return Err(coverage(
                "profile enclosure has no two-dimensional interior",
            ));
        }
    }
    let [c, b, a] = coefficients(&scene.surface, frame)?;
    if a.total_degree(dec)? != Some(0) {
        return Err(coverage("ray equation is not a constant-leading quadratic"));
    }
    let av = a.eval_real(&[Real::zero(), Real::zero()]);
    if sign(&av, dec)? != Sign::Positive {
        return Err(coverage(
            "ray does not have a positive quadratic leading coefficient",
        ));
    }
    let cells = bounds.each_ref().map(|[lo, hi]| interval(lo, hi));
    let aa = a.eval_ball(&cells, dec)?;
    let bb = b.eval_ball(&cells, dec)?;
    let cc = c.eval_ball(&cells, dec)?;
    let discriminant = &b * &b - (&a * &c).scale(&Real::from_i64(4));
    let dd = discriminant.eval_ball(&cells, dec)?;
    if dd.lo().sign() != Sign::Positive {
        return Err(coverage(
            "whole-profile positive discriminant is not certified; silhouette subdivision is unfinished",
        ));
    }
    let branch = if cc.hi().sign() == Sign::Negative {
        Branch::Exit
    } else if cc.lo().sign() == Sign::Positive && bb.hi().sign() == Sign::Negative {
        Branch::Entry
    } else {
        return Err(coverage(
            "the profile crosses the target boundary or has no uniformly positive branch",
        ));
    };
    let roots = dd
        .sqrt()
        .ok_or_else(|| coverage("positive discriminant enclosure has no square root"))?;
    let denominator = &aa * &Ball::from_rational(&Rational::from_i64(2), 128);
    let upper = (-&bb + roots)
        .checked_div(&denominator)
        .ok_or_else(|| coverage("quadratic denominator is not separated"))?
        .hi();
    if upper.sign() != Sign::Positive {
        return Err(coverage("positive root bound is not certified"));
    }
    let extent = upper + Rational::one();
    // Every intermediate point from the start through both positive hits must
    // lie strictly inside every target trim. This excludes an earlier cap hit.
    let time = interval(&Rational::zero(), &extent);
    let world: [Ball; 3] = std::array::from_fn::<_, 3, _>(|i| {
        let o = [&frame.origin.x, &frame.origin.y, &frame.origin.z][i];
        let u = [&frame.u.x, &frame.u.y, &frame.u.z][i];
        let v = [&frame.v.x, &frame.v.y, &frame.v.z][i];
        let w = [&frame.direction.x, &frame.direction.y, &frame.direction.z][i];
        // Frames compiled from recipe poses are rational; algebraic caller
        // frames are enclosed below through the same budgeted evaluator.
        Ok::<_, Refusal>(
            enclose(o, dec)?
                + &enclose(u, dec)? * &cells[0]
                + &enclose(v, dec)? * &cells[1]
                + &enclose(w, dec)? * &time,
        )
    })
    .into_iter()
    .collect::<Result<Vec<_>, _>>()?
    .try_into()
    .expect("three coordinates");
    for trim in &scene.trims {
        let value = Surface::Plane(trim.clone())
            .implicit()
            .ok_or_else(|| coverage("stock trim has no implicit polynomial"))?
            .eval_ball(&world, dec)?;
        if value.hi().sign() != Sign::Negative {
            return Err(coverage(
                "the ray enclosure reaches a target trim; projected trim arrangement is unfinished",
            ));
        }
    }
    let at_origin = scene
        .surface
        .gradient(&frame.origin)
        .ok_or_else(|| coverage("limiting carrier has no analytic gradient"))?
        .dot(&frame.direction);
    let normal = scene
        .surface
        .gradient(&frame.origin.plus(&frame.direction))
        .ok_or_else(|| coverage("limiting carrier has no analytic gradient"))?
        .minus(
            &scene
                .surface
                .gradient(&frame.origin)
                .ok_or_else(|| coverage("limiting carrier has no analytic gradient"))?,
        );
    // grad(F).w = b(q)+2*a*t, and its world gradient is 2 Q w.
    let norm = normal.norm_squared();
    if sign(&norm, dec)? != Sign::Positive {
        return Err(coverage("quadratic midpoint plane is degenerate"));
    }
    let point = frame.origin.plus(&normal.scaled(&(-at_origin / norm)));
    let reference = perpendicular(&normal, dec)?;
    let midpoint = Plane::new(point, normal, reference, dec)?;
    Ok(Certificate {
        branch,
        midpoint,
        extent,
    })
}

pub(super) fn enclose(x: &Real, dec: &mut Decider) -> Result<Ball, Refusal> {
    if let Some(q) = x.rational_ref() {
        Ok(Ball::from_rational(q, 128))
    } else {
        dec.enclose_budgeted_at(x, 128)
    }
}

pub(super) fn perpendicular(n: &Vector3, dec: &mut Decider) -> Result<Vector3, Refusal> {
    for i in 0..3 {
        let v = n.cross(&axis(i));
        if sign(&v.norm_squared(), dec)? == Sign::Positive {
            return Ok(v);
        }
    }
    Err(coverage("normal is zero"))
}

fn interval(lo: &Rational, hi: &Rational) -> Ball {
    Ball::from_rational(lo, 128).hull(&Ball::from_rational(hi, 128))
}

#[cfg(test)]
#[path = "projected_tests.rs"]
mod tests;
