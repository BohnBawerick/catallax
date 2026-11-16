//! Planar cap wiring from retained stock rings and tool sections.

use std::collections::{BTreeMap, BTreeSet};

use catallax_numbers::{Decider, Rational, Real, Refusal};
use catallax_shapes::{BSplineCurve, Curve, Plane, Point3, Surface};
use catallax_topology::{
    Sense,
    numerical::{NumericalDefinition, NumericalFace, NumericalNativeTrim, NumericalTrim},
};

fn bad(detail: &str) -> Refusal {
    Refusal::new("hole.cap-boundary", detail)
}

/// Appends one generating cap with an outer boundary and optional bore boundary.
/// Each boundary supplies unordered, already oriented shared edge uses. All uses
/// must form one closed walk. Branches and extra disconnected walks refuse.
/// The original plane, 3D carriers and finite parameters are retained. Projected
/// native pcurves require subsequent numerical body admission under its preset
/// tolerance; this operation grants no incidence or STEP authority.
pub fn append_planar_cap(
    definition: &mut NumericalDefinition,
    native: &mut Vec<NumericalNativeTrim>,
    key: &str,
    plane: &Plane,
    sense: Sense,
    boundaries: &[Vec<(usize, Sense)>],
    dec: &mut Decider,
) -> Result<usize, Refusal> {
    dec.charge_work(1)?;
    if key.trim().is_empty() || boundaries.is_empty() || boundaries.len() > 2 {
        return Err(bad(
            "a cap needs a key, an outer boundary and at most one bore",
        ));
    }
    for face in &definition.faces {
        dec.charge_work(1)?;
        if face.key == key {
            return Err(bad("cap face key is already in use"));
        }
    }
    for binding in native.iter() {
        dec.charge_work(1)?;
        if binding.face == key {
            return Err(bad("cap face key already has native bindings"));
        }
    }
    Plane::new(
        plane.origin.clone(),
        plane.normal.clone(),
        plane.ref_dir.clone(),
        dec,
    )?;
    let normal = plane
        .normal
        .scaled(&(&Real::one() / &plane.normal.norm_squared().sqrt()));
    let x = plane
        .ref_dir
        .scaled(&(&Real::one() / &plane.ref_dir.norm_squared().sqrt()));
    let y = normal.cross(&x);
    let project = |point: &Point3| {
        let relative = point.minus(&plane.origin);
        Point3::new(relative.dot(&x), relative.dot(&y), Real::zero())
    };
    let mut loops = Vec::new();
    let mut bindings = Vec::new();
    let mut used = BTreeSet::new();
    let mut parameters = BTreeMap::<usize, [Rational; 2]>::new();
    for (boundary, uses) in boundaries.iter().enumerate() {
        let ordered = walk(definition, uses, dec)?;
        let mut trims = Vec::new();
        for (use_index, (index, sense)) in ordered.into_iter().enumerate() {
            dec.charge_work(1)?;
            if !used.insert(index) {
                return Err(bad("a cap boundary repeats an edge"));
            }
            let edge = &definition.edges[index];
            let curve = match &edge.curve {
                Curve::BSpline(curve) => {
                    let controls = curve
                        .controls()
                        .map(|point| {
                            dec.charge_work(1)?;
                            Ok(project(point))
                        })
                        .collect::<Result<Vec<_>, Refusal>>()?;
                    let weights = curve
                        .weights()
                        .map(|weight| {
                            dec.charge_work(1)?;
                            Ok(weight.clone())
                        })
                        .collect::<Result<Vec<_>, Refusal>>()?;
                    dec.charge_work(curve.knots().len() as u64)?;
                    BSplineCurve::new(
                        curve.degree(),
                        curve.knots().to_vec(),
                        controls,
                        weights,
                        dec,
                    )?
                }
                Curve::Line(line) => {
                    let mut domain = edge.parameters.clone();
                    domain.sort();
                    BSplineCurve::new(
                        1,
                        vec![
                            domain[0].clone(),
                            domain[0].clone(),
                            domain[1].clone(),
                            domain[1].clone(),
                        ],
                        domain
                            .map(|t| project(&line.origin.plus(&line.dir.scaled(&Real::from(t)))))
                            .to_vec(),
                        vec![Real::one(); 2],
                        dec,
                    )?
                }
                _ => return Err(bad("cap edges require retained splines or finite lines")),
            };
            let mut vertices = edge.vertices;
            if sense == Sense::Reversed {
                vertices.reverse();
            }
            let mut uv = Vec::new();
            for vertex in vertices {
                if let Some(at) = parameters.get(&vertex) {
                    uv.push(at.clone());
                    continue;
                }
                let [x, y, z] = definition.vertices[vertex].position.clone().map(Real::from);
                let point = project(&Point3::new(x, y, z));
                let at = [
                    dec.enclose_budgeted_at(&point.x, 128)?.mid(),
                    dec.enclose_budgeted_at(&point.y, 128)?.mid(),
                ];
                uv.push(at.clone());
                parameters.insert(vertex, at);
            }
            trims.push(NumericalTrim {
                edge: index,
                sense,
                uv: uv.try_into().expect("two endpoints"),
            });
            bindings.push(NumericalNativeTrim {
                face: key.into(),
                boundary,
                use_index,
                edge: edge.key.clone(),
                curve,
            });
        }
        loops.push(trims);
    }
    let index = definition.faces.len();
    definition.faces.push(NumericalFace {
        key: key.into(),
        surface: Surface::Plane(plane.clone()),
        sense,
        loops,
    });
    native.extend(bindings);
    Ok(index)
}

fn walk(
    definition: &NumericalDefinition,
    uses: &[(usize, Sense)],
    dec: &mut Decider,
) -> Result<Vec<(usize, Sense)>, Refusal> {
    let mut outgoing = BTreeMap::new();
    let mut incoming = BTreeSet::new();
    for (index, sense) in uses {
        dec.charge_work(1)?;
        let edge = definition
            .edges
            .get(*index)
            .ok_or_else(|| bad("cap use names an unknown shared edge"))?;
        if edge
            .vertices
            .iter()
            .any(|v| *v >= definition.vertices.len())
        {
            return Err(bad("cap edge names an unknown vertex"));
        }
        let mut vertices = edge.vertices;
        if *sense == Sense::Reversed {
            vertices.reverse();
        }
        if outgoing
            .insert(vertices[0], (*index, *sense, vertices[1]))
            .is_some()
            || !incoming.insert(vertices[1])
        {
            return Err(bad("cap boundary has a branch or repeated vertex"));
        }
    }
    let Some((&start, _)) = outgoing.first_key_value() else {
        return Err(bad("cap boundary is empty"));
    };
    let mut at = start;
    let mut ordered = Vec::new();
    while let Some((edge, sense, end)) = outgoing.remove(&at) {
        dec.charge_work(1)?;
        ordered.push((edge, sense));
        at = end;
    }
    if at != start || !outgoing.is_empty() {
        return Err(bad("cap boundary is open or contains another closed walk"));
    }
    Ok(ordered)
}
