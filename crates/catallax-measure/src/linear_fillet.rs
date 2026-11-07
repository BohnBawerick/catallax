//! Whole-face integrals of exact affine-radius contact sections.
//!
//! `S(t,v) = E0 + v L + r(v) h(t)`, with rational quadratic `h`.
//! Thus `St x Sv = r(v) h'(t) x (L + dr h(t))`. Surface-area density is
//! affine in `v`; all volume, centroid and second-moment fluxes have degree
//! at most four in `v`. Three-point Gauss integration is exact in that
//! coordinate, including a collapsed zero-radius boundary. The remaining
//! rational/square-root integral uses the validated whole-panel Taylor engine.

use catallax_numbers::{Decider, Refusal};
use catallax_shapes::{LinearFillet, RPoint3};
use catallax_topology::{FaceId, LinearFilletBand, Model, Sense};

use crate::integrate::{Builder, Expr, Program};
use crate::planar::{real_components as components, real_direction as direction};
use crate::surface::Property;

type Vector<'a> = [Expr<'a>; 3];

pub(crate) fn prepare(
    model: &Model,
    face: FaceId,
    law: &LinearFillet,
    anchor: &RPoint3,
    property: Property,
    dec: &mut Decider,
) -> Result<crate::pipe::Band, Refusal> {
    LinearFilletBand::of_face(model, face, dec)?.ok_or_else(|| {
        Refusal::new(
            "measure.unsupported-surface",
            "a variable fillet integral requires its complete checked band",
        )
    })?;
    Ok(crate::pipe::Band::parameter(forms(
        law,
        anchor,
        property,
        model.face(face).sense,
    )))
}

fn forms(law: &LinearFillet, anchor: &RPoint3, property: Property, sense: Sense) -> Program {
    let k = Builder::default();
    let zero = k.k(0);
    let one = k.k(1);
    let two = k.k(2);
    let t = k.input(0);
    let other = one - t;
    let weight = k.exact(&law.carrier().weights()[2]);
    let denominator = other.square() + two * weight * t * other + t.square();
    let derivative = two * (-other + weight * (one - two * t) + t);
    let first = direction(&law.shift().plus(&law.rays()[0])).map(|x| k.exact(x));
    let last = direction(&law.shift().plus(&law.rays()[1])).map(|x| k.exact(x));
    let numerator: Vector =
        std::array::from_fn(|i| first[i] * other.square() + last[i] * t.square());
    let h: Vector = numerator.map(|n| n / denominator);
    let dh: Vector = std::array::from_fn(|i| {
        (two * (-first[i] * other + last[i] * t) * denominator - numerator[i] * derivative)
            / denominator.square()
    });
    let along = direction(law.along()).map(|x| k.exact(x));
    let base = components(&law.edge()[0]).map(|x| k.exact(x));
    let anchor = components(&anchor.to_real()).map(|x| k.exact(x));
    let radius = law.radii().each_ref().map(|r| k.exact(r));
    let dr = radius[1] - radius[0];
    let station: Vector = std::array::from_fn(|i| along[i] + dr * h[i]);
    let cross: Vector = [
        dh[1] * station[2] - dh[2] * station[1],
        dh[2] * station[0] - dh[0] * station[2],
        dh[0] * station[1] - dh[1] * station[0],
    ];
    let area = if matches!(property, Property::Area | Property::FaceCentroid) {
        cross.iter().fold(zero, |sum, x| sum + x.square()).sqrt()
    } else {
        zero
    };
    let orientation = k.k(if sense == Sense::Same { 1 } else { -1 });
    let root = (k.k(3) / k.k(5)).sqrt();
    let nodes = [(one - root) / two, one / two, (one + root) / two];
    let weights = [k.k(5) / k.k(18), k.k(4) / k.k(9), k.k(5) / k.k(18)];
    let count = if matches!(property, Property::Centroid | Property::FaceCentroid) {
        4
    } else {
        1
    };
    let mut outputs = [zero; 4];
    for (v, weight) in nodes.into_iter().zip(weights) {
        let r = radius[0] + v * dr;
        let p: Vector = std::array::from_fn(|i| base[i] - anchor[i] + v * along[i] + r * h[i]);
        let flux = p.iter().zip(cross).fold(zero, |sum, (p, n)| sum + *p * n) * r * orientation;
        let values = match property {
            Property::Area => [r * area, zero, zero, zero],
            Property::Volume => [flux / k.k(3), zero, zero, zero],
            Property::Centroid => std::array::from_fn(|i| {
                if i == 0 {
                    flux / k.k(3)
                } else {
                    p[i - 1].square() * cross[i - 1] * r * orientation / two
                }
            }),
            Property::FaceCentroid => std::array::from_fn(|i| {
                if i == 0 {
                    r * area
                } else {
                    p[i - 1] * r * area
                }
            }),
            Property::SecondMoment(i, j) => [p[i] * p[j] * flux / k.k(5), zero, zero, zero],
        };
        for (sum, value) in outputs.iter_mut().zip(values) {
            *sum = *sum + weight * value;
        }
    }
    k.finish(&outputs[..count])
}
