//! Whole-band integrals of the exact band a ball sweeps over a mitre.
//!
//! A band point is `P = c + r e`, with the centre `c(t) = O + alpha cos t x +
//! k sin t y` for `alpha = k / sa`, and `e = cos h m1 + sin h e2` running from
//! `m1` to `m2` over the angle `f` between them. With `sc^2 + sa^2 = 1`,
//!
//! ```text
//! C = m1 . m2 = 1 - 2 sc^2 cos^2 t,   Q = sqrt(1 - sc^2 cos^2 t) >= sa,
//! S = sin f = 2 sc cos t Q,           e2 = (z + sc cos t m1) / Q,
//! ```
//!
//! so nothing divides by `S`, which vanishes at the apex, and the section's
//! angle is `f = pi - 2 phi` for `cos phi = sc cos t`, `sin phi = Q`. Both rays meet the
//! unit tangent's derivative in `T' . m = -k / |c'|`, because
//! `alpha sa = k`, which gives the surface element
//! `dA = r (|c'| + r k (cos h + sin h sc cos t / Q) / |c'|) dh dt`. Its
//! integrands are polynomials in `cos h` and `sin h` with exact primitives
//! over the section, which leaves one integral in `tau = tan(t/2)` over
//! `[0, 1]`, with `dt = 2 dtau / (1 + tau^2)`, for the Taylor engine. The
//! chart's own normal points along `e`, away from the ball.

use std::collections::BTreeMap;

use catallax_numbers::{Decider, Real, Refusal};
use catallax_shapes::{RPoint3, RollingMitre};
use catallax_topology::{FaceId, Model, Sense};

use crate::integrate::{Builder, Expr, Program};
use crate::pipe::{Primitives, Section};
use crate::planar::{real_components as components, real_direction as direction};
use crate::surface::Property;

type Vector<'a> = [Expr<'a>; 3];

fn unsupported(detail: &'static str) -> Refusal {
    Refusal::new("measure.unsupported-surface", detail)
}

/// A complete band: one loop of its end section and its two rails.
pub(crate) fn prepare(
    model: &Model,
    face: FaceId,
    band: &RollingMitre,
    anchor: &RPoint3,
    property: Property,
    dec: &mut Decider,
) -> Result<crate::pipe::Band, Refusal> {
    catallax_topology::RollingMitreBand::of_face(model, face, dec)?
        .ok_or_else(|| unsupported("a trimmed rolling mitre band has no measure"))?;
    Ok(crate::pipe::Band::parameter(forms(
        band,
        anchor,
        property,
        model.face(face).sense,
    )))
}

fn forms(band: &RollingMitre, anchor: &RPoint3, property: Property, sense: Sense) -> Program {
    let k = Builder::default();
    let exact = |value: &Real| k.exact(value);
    let [out, up, side] = band
        .frame()
        .each_ref()
        .map(|v| direction(v).map(|x| k.exact(x)));
    let offset = band.round() - band.radius();
    let half = Real::from(catallax_numbers::Rational::from_ratio(1, 2).expect("nonzero"));
    let cosine = band.axes()[0].dot(&band.axes()[1]);
    let along = exact(&(&(&Real::one() - &cosine) * &half).sqrt());
    let across = exact(&(&(&Real::one() + &cosine) * &half).sqrt());
    let (radius, offset) = (exact(band.radius()), exact(&offset));
    let alpha = offset / along;
    let origin = components(band.origin()).map(|x| k.exact(x));
    let anchor = components(&anchor.to_real()).map(|x| k.exact(x));
    let tau = k.input(0);
    let wide = k.k(1) + tau.square();
    let cos = (k.k(1) - tau.square()) / wide;
    let sin = k.k(2) * tau / wide;
    let turn = k.k(2) / wide;
    let zero = k.k(0);
    let centre: Vector = std::array::from_fn(|i| {
        origin[i] - anchor[i] + alpha * cos * out[i] + offset * sin * up[i]
    });
    let speed = (alpha.square() * sin.square() + offset.square() * cos.square()).sqrt();
    let first: Vector =
        std::array::from_fn(|i| along * cos * out[i] + sin * up[i] - across * cos * side[i]);
    let c = k.k(1) - k.k(2) * across.square() * cos.square();
    let q = (k.k(1) - across.square() * cos.square()).sqrt();
    let s = k.k(2) * across * cos * q;
    // `f = pi - 2 phi` for `cos phi = sc cos t` and `sin phi = Q >= sa`, so
    // the angle's sine stays clear of zero where `S` vanishes at the apex.
    let quarter = k.k(0).angle_to(k.k(1));
    let angle = k.k(2) * (quarter - (across * cos).angle_to(q));
    let second: Vector = std::array::from_fn(|i| (side[i] + across * cos * first[i]) / q);
    let position: [Section; 3] =
        std::array::from_fn(|i| Section::linear(centre[i], radius * first[i], radius * second[i]));
    let ray: [Section; 3] = std::array::from_fn(|i| Section::linear(zero, first[i], second[i]));
    // The surface element per `dtau`.
    let reach = radius.square() * offset / speed;
    let density = Section::linear(
        radius * speed * turn,
        reach * turn,
        reach * across * cos / q * turn,
    );
    // `P . e = r + c . e`, which keeps the flux linear in the section.
    let [towards_first, towards_second] =
        [first, second].map(|ray| (0..3).fold(zero, |sum, i| sum + centre[i] * ray[i]));
    let flux = Section::linear(radius, towards_first, towards_second);
    let carried = flux.times(&density);
    let side_sign = k.k(if sense == Sense::Same { 1 } else { -1 });
    let mut primitives = Primitives {
        builder: &k,
        angle,
        cosine: c,
        sine: s,
        known: BTreeMap::new(),
    };
    let outputs = match property {
        Property::Area => vec![primitives.integral(&density)],
        Property::Volume => vec![side_sign * primitives.integral(&carried) / k.k(3)],
        Property::Centroid => {
            let volume = side_sign * primitives.integral(&carried) / k.k(3);
            let mut outputs = vec![volume];
            for i in 0..3 {
                let form = position[i]
                    .times(&position[i])
                    .times(&ray[i].times(&density));
                outputs.push(side_sign * primitives.integral(&form) / k.k(2));
            }
            outputs
        }
        Property::FaceCentroid => {
            let mut outputs = vec![primitives.integral(&density)];
            for at in &position {
                outputs.push(primitives.integral(&at.times(&density)));
            }
            outputs
        }
        Property::SecondMoment(i, j) => {
            let form = position[i].times(&position[j]).times(&carried);
            vec![side_sign * primitives.integral(&form) / k.k(5)]
        }
    };
    k.finish(&outputs)
}
