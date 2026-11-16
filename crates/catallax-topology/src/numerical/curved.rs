//! Sampled parameter walks for native spline edges on planar faces.

use super::*;

pub(super) fn projected(face: &NumericalFace, edges: &[NumericalEdge]) -> bool {
    matches!(face.surface, Surface::Plane(_))
        && face.loops.iter().flatten().any(|trim| {
            edges
                .get(trim.edge)
                .is_some_and(|edge| matches!(edge.curve, Curve::BSpline(_)))
        })
}

pub(super) fn parameters(
    surface: &Surface,
    point: &[Ball; 3],
    dec: &Decider,
) -> Result<[Rational; 2], Refusal> {
    let at = surface
        .parameters_of(point, dec)?
        .ok_or_else(|| bad("trim", "planar curve projection is unresolved"))?;
    Ok(at.each_ref().map(Ball::mid))
}

pub(super) fn check_loops(
    face: &NumericalFace,
    edges: &[NumericalEdge],
    dec: &mut Decider,
) -> Result<(), Refusal> {
    if !projected(face, edges) {
        return super::check_loops(face, dec);
    }
    let mut sampled = face.clone();
    sampled.loops.clear();
    for boundary in &face.loops {
        let mut polygon = Vec::new();
        for trim in boundary {
            let edge = edges
                .get(trim.edge)
                .ok_or_else(|| bad("trim", "trim names an unknown edge"))?;
            let [a, b] = if trim.sense == Sense::Same {
                edge.parameters.each_ref()
            } else {
                [&edge.parameters[1], &edge.parameters[0]]
            };
            let mut previous = trim.uv[0].clone();
            for i in 1..=SUBDIVISIONS {
                dec.charge_work(1)?;
                let uv = if i == SUBDIVISIONS {
                    trim.uv[1].clone()
                } else {
                    let share = Rational::from_ratio(i, SUBDIVISIONS).unwrap();
                    let t = a + &(&share * &(b - a));
                    let point = crate::numerical::curve_at(&edge.curve, &t, dec)?
                        .ok_or_else(|| bad("unsupported", "edge has no parameter evaluator"))?;
                    parameters(&face.surface, &point, dec)?
                };
                polygon.push(NumericalTrim {
                    edge: trim.edge,
                    sense: trim.sense,
                    uv: [previous, uv.clone()],
                });
                previous = uv;
            }
        }
        sampled.loops.push(polygon);
    }
    super::check_loops(&sampled, dec)
}
