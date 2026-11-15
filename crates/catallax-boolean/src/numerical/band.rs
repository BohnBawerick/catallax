//! Retained circular ring ownership and bounded finite cone-height rounding.

use super::*;
use catallax_topology::Loop;

pub(super) fn read(
    model: &Model,
    face: FaceId,
    cone: &Cone,
    tolerance: &Rational,
    dec: &mut Decider,
) -> Result<Option<(Vec<NumericalConeRing>, Vec<NumericalCheck>)>, Refusal> {
    if model.face(face).loops.len() != 2 {
        return Ok(None);
    }
    let axis = unit(&cone.axis, dec)?;
    let mut rings = Vec::new();
    let mut samples = Vec::new();
    for &boundary in &model.face(face).loops {
        dec.charge_work(1)?;
        let Loop::Coedges(uses) = model.loop_(boundary) else {
            return Ok(None);
        };
        let [use_id] = uses.as_slice() else {
            return Ok(None);
        };
        let coedge = model.coedge(*use_id);
        let edge = model.edge(coedge.edge);
        let Curve::Circle(circle) = model.geom().curve(edge.curve) else {
            return Ok(None);
        };
        if edge.start.is_some() || edge.end.is_some() {
            return Ok(None);
        }
        let delta = circle.center.minus(&cone.apex);
        let height = delta.dot(&axis);
        let bound = interval(&height, dec)?;
        let rounded = bound.mid();
        let radius = circle.radius_sq.sqrt();
        let normal = unit(&circle.normal, dec)?;
        let centre_error = delta.minus(&axis.scaled(&height)).norm_squared().sqrt();
        let tilt = normal.cross(&axis).norm_squared().sqrt();
        let radial_error = &radius - &(&cone.slope * &height);
        let rounding = upper(&(&bound - &Ball::from_rational(&rounded, BITS)));
        let error = interval(&centre_error, dec)?.hi()
            + interval(&(&(&Real::from(2) * &radius) * &tilt), dec)?.hi()
            + upper(&interval(&radial_error, dec)?)
            + (q(1) + upper(&interval(&cone.slope, dec)?)) * rounding;
        if error > tolerance / &q(16) {
            return Err(bad(
                "original ring differs from the finite cone-band boundary beyond its tolerance share",
            ));
        }
        samples.push(NumericalSample {
            parameters: vec![rounded.clone()],
            distance_upper_mm: error,
        });
        rings.push(NumericalConeRing {
            original_loop: boundary,
            coedge: *use_id,
            edge: coedge.edge,
            sense: coedge.sense,
            edge_sense: edge.sense,
            circle: circle.clone(),
            height_mm: rounded,
            original_height_mm: [bound.lo(), bound.hi()],
        });
    }
    rings.sort_by(|a, b| a.height_mm.cmp(&b.height_mm));
    if rings[0].height_mm.sign() != Sign::Positive || rings[0].height_mm >= rings[1].height_mm {
        return Err(bad(
            "finite cone ring heights must be distinct and strictly positive",
        ));
    }
    Ok(Some((
        rings,
        vec![NumericalCheck::Samples {
            name: "finite.ring-cone-band-deviation-bound".into(),
            scope: format!("{face:?}"),
            samples,
        }],
    )))
}

pub(super) fn angle(
    ring: &NumericalConeRing,
    cone: &Cone,
    u: &Rational,
    dec: &mut Decider,
) -> Result<[Rational; 2], Refusal> {
    let axis = unit(&cone.axis, dec)?;
    let x = unit(
        &cone.ref_dir.minus(&axis.scaled(&cone.ref_dir.dot(&axis))),
        dec,
    )?;
    let normal = unit(&ring.circle.normal, dec)?;
    let rx = unit(
        &ring
            .circle
            .ref_dir
            .minus(&normal.scaled(&ring.circle.ref_dir.dot(&normal))),
        dec,
    )?;
    let ry = normal.cross(&rx);
    let cosine = interval(&x.dot(&rx), dec)?;
    let sine = interval(&x.dot(&ry), dec)?;
    let pi = Ball::pi_checked(BITS, dec)?;
    let phase = if cosine.hi().sign() == Sign::Negative {
        sine.checked_div(&cosine)
            .ok_or_else(|| bad("ring phase denominator vanished"))?
            .atan_checked(dec)?
            + pi.clone()
    } else {
        Ball::atan2(&sine, &cosine).ok_or_else(|| bad("ring phase is unresolved"))?
    };
    let orientation = interval(&axis.dot(&normal), dec)?;
    let direction = if orientation.lo().sign() == Sign::Positive {
        q(1)
    } else if orientation.hi().sign() == Sign::Negative {
        q(-1)
    } else {
        return Err(bad("ring parameter orientation is unresolved"));
    };
    let quarter = (0..4)
        .find(|i| u >= &q(*i) && u <= &q(*i + 1))
        .ok_or_else(|| bad("ring source parameter leaves the retained full circle"))?;
    let local = u - &q(quarter);
    let h = (&q(7) * &local) / (q(10) - &q(3) * &local);
    let angle = Ball::from_rational(&h, BITS).atan_checked(dec)? * Ball::from_rational(&q(2), BITS)
        + pi.clone()
            * Ball::from_rational(
                &Rational::from_ratio(quarter, 2).expect("fixed denominator"),
                BITS,
            );
    let result = phase + angle * Ball::from_rational(&direction, BITS);
    Ok([result.lo(), result.hi()])
}
