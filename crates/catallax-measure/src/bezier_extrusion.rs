//! Finite ruled Bezier bands, with exact polynomial integration along the axis.

use crate::integrate::{Builder, Expr, Program};
use crate::planar::{real_components as components, real_direction as direction};
use crate::surface::Property;
use catallax_numbers::{Decider, Refusal};
use catallax_shapes::RPoint3;
use catallax_topology::{BezierExtrusionBand, FaceId, Model, Sense};

type Vector<'a> = [Expr<'a>; 3];

fn evaluate<'a>(mut controls: Vec<Vector<'a>>, u: Expr<'a>, one: Expr<'a>) -> Vector<'a> {
    for count in (1..controls.len()).rev() {
        for i in 0..count {
            controls[i] =
                std::array::from_fn(|j| (one - u) * controls[i][j] + u * controls[i + 1][j]);
        }
    }
    controls[0]
}

pub(crate) fn prepare(
    model: &Model,
    face: FaceId,
    anchor: &RPoint3,
    property: Property,
    dec: &mut Decider,
) -> Result<crate::pipe::Band, Refusal> {
    let band = BezierExtrusionBand::of_face(model, face, dec)?.ok_or_else(|| {
        Refusal::new(
            "measure.unsupported-surface",
            "a Bezier integral needs its checked finite band",
        )
    })?;
    Ok(crate::pipe::Band::parameter(forms(
        &band,
        anchor,
        property,
        model.face(face).sense,
    )))
}

fn forms(
    band: &BezierExtrusionBand<'_>,
    anchor: &RPoint3,
    property: Property,
    sense: Sense,
) -> Program {
    let k = Builder::default();
    let zero = k.k(0);
    let one = k.k(1);
    let two = k.k(2);
    let u = k.input(0);
    let controls: Vec<Vector> = band
        .support()
        .controls()
        .iter()
        .map(|p| components(p).map(|x| k.exact(x)))
        .collect();
    let degree = k.k((controls.len() - 1) as i64);
    let tangent = evaluate(
        controls
            .windows(2)
            .map(|p| std::array::from_fn(|i| degree * (p[1][i] - p[0][i])))
            .collect(),
        u,
        one,
    );
    let curve = evaluate(controls, u, one);
    let [lo, hi] = band.axial().map(|v| k.exact(v));
    let axis = direction(band.support().direction()).map(|x| k.exact(x));
    let displacement: Vector = axis.map(|x| (hi - lo) * x);
    let anchor = components(&anchor.to_real()).map(|x| k.exact(x));
    let base: Vector = std::array::from_fn(|i| curve[i] - anchor[i] + lo * axis[i]);
    let normal: Vector = [
        tangent[1] * displacement[2] - tangent[2] * displacement[1],
        tangent[2] * displacement[0] - tangent[0] * displacement[2],
        tangent[0] * displacement[1] - tangent[1] * displacement[0],
    ];
    let area = if matches!(property, Property::Area | Property::FaceCentroid) {
        normal.iter().fold(zero, |sum, x| sum + x.square()).sqrt()
    } else {
        zero
    };
    let orientation = k.k(if sense == Sense::Same { 1 } else { -1 });
    // Every flux and moment has axial degree at most three. Three-point
    // Gauss integration is exact there; scalar-parameter integration proves
    // the remaining polynomial or square-root integral along the curve.
    let root = (k.k(3) / k.k(5)).sqrt();
    let nodes = [(one - root) / two, one / two, (one + root) / two];
    let weights = [k.k(5) / k.k(18), k.k(4) / k.k(9), k.k(5) / k.k(18)];
    let mut outputs = [zero; 4];
    for (v, weight) in nodes.into_iter().zip(weights) {
        let p: Vector = std::array::from_fn(|i| base[i] + v * displacement[i]);
        let flux = p.iter().zip(normal).fold(zero, |sum, (p, n)| sum + *p * n) * orientation;
        let values = match property {
            Property::Volume => [flux / k.k(3), zero, zero, zero],
            Property::Area => [area, zero, zero, zero],
            Property::Centroid => std::array::from_fn(|i| {
                if i == 0 {
                    flux / k.k(3)
                } else {
                    p[i - 1].square() * normal[i - 1] * orientation / two
                }
            }),
            Property::FaceCentroid => {
                std::array::from_fn(|i| if i == 0 { area } else { p[i - 1] * area })
            }
            Property::SecondMoment(i, j) => [p[i] * p[j] * flux / k.k(5), zero, zero, zero],
        };
        for (sum, value) in outputs.iter_mut().zip(values) {
            *sum = *sum + weight * value;
        }
    }
    let count = if matches!(property, Property::Centroid | Property::FaceCentroid) {
        4
    } else {
        1
    };
    k.finish(&outputs[..count])
}
