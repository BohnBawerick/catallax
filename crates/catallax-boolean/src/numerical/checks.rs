//! Actual generated edge endpoints and per-use native pcurve joins.

use super::*;
use std::collections::BTreeMap;

fn on_surface(
    surface: &Surface,
    uv: &[Rational; 2],
    dec: &mut Decider,
) -> Result<[Ball; 3], Refusal> {
    let Surface::BSpline(SplineSupport::RationalPatch(patch)) = surface else {
        return Err(bad(
            "finite cone partition lost its retained rational patch",
        ));
    };
    patch
        .carrier()
        .evaluate_box([&uv[0], &uv[0]], [&uv[1], &uv[1]], BITS, dec)
}

pub(super) fn open(
    definition: &NumericalDefinition,
    native: &[NumericalNativeTrim],
    dec: &mut Decider,
) -> Result<Vec<NumericalCheck>, Refusal> {
    let bindings: BTreeMap<_, _> = native
        .iter()
        .map(|t| ((t.face.as_str(), t.boundary, t.use_index), t))
        .collect();
    if bindings.len() != native.len() {
        return Err(bad("duplicate generated native binding address"));
    }
    let mut vertices = Vec::new();
    let mut coupled = Vec::new();
    let mut declared = Vec::new();
    let mut joins = Vec::new();
    for face in &definition.faces {
        for (boundary, row) in face.loops.iter().enumerate() {
            if row.is_empty() {
                return Err(bad("generated face contains an empty loop"));
            }
            let mut actual = Vec::new();
            for (use_index, trim) in row.iter().enumerate() {
                let edge = &definition.edges[trim.edge];
                let mut parameters = edge.parameters.each_ref();
                let mut endpoints = edge.vertices;
                if trim.sense == Sense::Reversed {
                    parameters.swap(0, 1);
                    endpoints.swap(0, 1)
                }
                let binding = bindings.get(&(face.key.as_str(), boundary, use_index));
                if binding.is_some_and(|binding| binding.edge != edge.key) {
                    return Err(bad("generated native binding names another edge"));
                }
                let mut physical = Vec::new();
                for end in 0..2 {
                    dec.charge_work(1)?;
                    let t = parameters[end];
                    let xyz = if let Curve::BSpline(curve) = &edge.curve {
                        curve.evaluate_box([t, t], BITS, dec)?
                    } else {
                        edge.curve
                            .evaluate(&Ball::from_rational(t, BITS), dec)?
                            .ok_or_else(|| bad("generated edge has no finite evaluator"))?
                    };
                    let uv = if let Some(binding) = binding {
                        let point = events::rational_point(
                            binding.curve.point_at(&Real::from(t.clone()), dec)?,
                        )?;
                        [point[0].clone(), point[1].clone()]
                    } else {
                        trim.uv[end].clone()
                    };
                    let on = on_surface(&face.surface, &uv, dec)?;
                    let logical = on_surface(&face.surface, &trim.uv[end], dec)?;
                    vertices.push(NumericalSample {
                        parameters: vec![t.clone()],
                        distance_upper_mm: distance(
                            &xyz,
                            &events::box_point(&definition.vertices[endpoints[end]].position),
                        )?,
                    });
                    coupled.push(NumericalSample {
                        parameters: vec![t.clone(), uv[0].clone(), uv[1].clone()],
                        distance_upper_mm: distance(&xyz, &on)?,
                    });
                    declared.push(NumericalSample {
                        parameters: vec![t.clone(), uv[0].clone(), uv[1].clone()],
                        distance_upper_mm: distance(&on, &logical)?,
                    });
                    physical.push(on);
                }
                actual.push(physical);
            }
            for i in 0..actual.len() {
                joins.push(NumericalSample {
                    parameters: vec![q(boundary as i64), q(i as i64)],
                    distance_upper_mm: distance(&actual[i][1], &actual[(i + 1) % actual.len()][0])?,
                });
            }
        }
    }
    if definition.faces.is_empty() {
        return Ok(vec![NumericalCheck::Passed {
            name: "finite.no-retained-cone-regions".into(),
            scope: definition.name.clone(),
        }]);
    }
    Ok([
        ("finite.generated-edge-vertex-endpoints", vertices),
        ("finite.generated-native-endpoint-coupling", coupled),
        (
            "finite.generated-declared-endpoint-correspondence",
            declared,
        ),
        ("finite.actual-native-pcurve-joins", joins),
    ]
    .into_iter()
    .map(|(name, samples)| NumericalCheck::Samples {
        name: name.into(),
        scope: definition.name.clone(),
        samples,
    })
    .collect())
}
