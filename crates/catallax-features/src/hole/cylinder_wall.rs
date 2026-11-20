//! Cylinder quarter walls assembled from owned shared boundary curves.

use std::collections::{BTreeMap, BTreeSet};

use catallax_boolean::numerical::NumericalCylinderSource;
use catallax_numbers::{Ball, Decider, Rational, Real, Refusal, Sign};
use catallax_shapes::{BSplineCurve, Curve, Point3};
use catallax_topology::{
    Sense,
    numerical::{
        NumericalDefinition, NumericalEdge, NumericalFace, NumericalNativeTrim, NumericalTrim,
    },
};

/// One retained shared edge, oriented with material on the face's interior side.
#[derive(Clone, Debug)]
pub struct CylinderWallBoundary {
    /// Index in the caller's common edge table.
    pub edge: usize,
    /// Original canonical cylinder quarter.
    pub quarter: usize,
    /// Traversal relative to the unchanged 3D edge parameter.
    pub sense: Sense,
    /// Coupled native UV definition at that same parameter.
    pub curve: BSplineCurve,
}

fn bad(detail: impl Into<String>) -> Refusal {
    Refusal::new("hole.cylinder-boundary", detail)
}
fn q(value: i64) -> Rational {
    Rational::from_i64(value)
}
fn at(curve: &BSplineCurve, t: &Rational, dec: &mut Decider) -> Result<[Rational; 2], Refusal> {
    let point = curve.evaluate(&Ball::from_rational(t, 128), dec)?;
    Ok([point[0].mid(), point[1].mid()])
}
fn ends(definition: &NumericalDefinition, trim: &NumericalTrim) -> [usize; 2] {
    let mut vertices = definition.edges[trim.edge].vertices;
    if trim.sense == Sense::Reversed {
        vertices.reverse();
    }
    vertices
}
fn point(
    source: &NumericalCylinderSource,
    uv: &[Rational; 2],
    dec: &Decider,
) -> Result<[Rational; 3], Refusal> {
    Ok(source
        .source()
        .evaluate(
            &Ball::from_rational(&uv[0], 128),
            &Ball::from_rational(&uv[1], 128),
            dec,
        )?
        .map(|v| v.mid()))
}

/// Wires the supplied boundaries and fills only unmatched native quarter seams.
/// Shared 3D curves, finite parameters and native maps are retained. All uses
/// must be consumed by closed directed walks. The generating result still needs
/// complete numerical body admission, including composition and material checks.
pub fn append_cylinder_wall(
    definition: &mut NumericalDefinition,
    native: &mut Vec<NumericalNativeTrim>,
    key: &str,
    source: &NumericalCylinderSource,
    sense: Sense,
    boundaries: &[CylinderWallBoundary],
    dec: &mut Decider,
) -> Result<Vec<usize>, Refusal> {
    dec.charge_work(
        (definition.vertices.len()
            + definition.edges.len()
            + definition.faces.len()
            + native.len()
            + boundaries.len()) as u64,
    )?;
    if key.trim().is_empty() || boundaries.iter().any(|boundary| boundary.quarter >= 4) {
        return Err(bad(
            "cylinder boundary needs a key and canonical quarter IDs",
        ));
    }
    let mut candidate = definition.clone();
    let mut added = Vec::new();
    let mut faces = Vec::new();
    let tolerance = source.evidence().tolerance_mm();
    let mut seams = BTreeMap::new();
    for quarter in 0..4 {
        dec.charge_work(boundaries.len() as u64)?;
        let mut uses = Vec::new();
        let mut curves = BTreeMap::new();
        let mut vertices = BTreeMap::<usize, [Rational; 2]>::new();
        for boundary in boundaries.iter().filter(|b| b.quarter == quarter) {
            dec.charge_work(1)?;
            let edge = candidate
                .edges
                .get(boundary.edge)
                .ok_or_else(|| bad("unknown cylinder boundary edge"))?;
            if edge
                .vertices
                .iter()
                .any(|vertex| *vertex >= candidate.vertices.len())
            {
                return Err(bad("cylinder edge names an unknown shared vertex"));
            }
            let [lo, hi] = boundary.curve.domain();
            if edge
                .parameters
                .iter()
                .any(|parameter| parameter < lo || parameter > hi)
            {
                return Err(bad(
                    "native cylinder map does not cover the shared edge interval",
                ));
            }
            let mut uv = edge
                .parameters
                .each_ref()
                .map(|t| at(&boundary.curve, t, dec))
                .into_iter()
                .collect::<Result<Vec<_>, _>>()?;
            for (vertex, at) in edge.vertices.into_iter().zip(&mut uv) {
                for side in [q(quarter as i64), q(quarter as i64 + 1)] {
                    if (&at[0] - &side).abs() <= *tolerance {
                        let snapped = [side, at[1].clone()];
                        let position = point(source, &snapped, dec)?;
                        if position
                            .iter()
                            .zip(&candidate.vertices[vertex].position)
                            .all(|(a, b)| (a - b).abs() <= tolerance.clone() / q(2))
                        {
                            *at = snapped;
                        }
                    }
                }
                if let Some(previous) = vertices.get(&vertex) {
                    let actual = point(source, at, dec)?;
                    let declared = point(source, previous, dec)?;
                    if actual
                        .iter()
                        .zip(&declared)
                        .any(|(a, b)| (a - b).abs() > tolerance.clone() / q(2))
                    {
                        return Err(bad("shared cylinder vertex has incompatible native maps"));
                    }
                    *at = previous.clone();
                } else {
                    vertices.insert(vertex, at.clone());
                }
            }
            if boundary.sense == Sense::Reversed {
                uv.reverse();
            }
            uses.push(NumericalTrim {
                tag: catallax_topology::numerical::CopyTag::Independent,
                edge: boundary.edge,
                sense: boundary.sense,
                uv: uv.try_into().expect("two endpoints"),
            });
            if curves
                .insert(boundary.edge, boundary.curve.clone())
                .is_some()
            {
                return Err(bad("a cylinder quarter repeats a boundary edge"));
            }
        }
        if uses.is_empty() {
            continue;
        }
        for right in [false, true] {
            let u = q(quarter as i64 + i64::from(right));
            dec.charge_work(vertices.len() as u64)?;
            let mut nodes: Vec<_> = vertices
                .iter()
                .filter(|(_, uv)| uv[0] == u)
                .map(|(id, uv)| (*id, uv.clone()))
                .collect();
            dec.charge_work((nodes.len() * nodes.len()) as u64)?;
            nodes.sort_by(|a, b| a.1[1].cmp(&b.1[1]));
            if right == (sense == Sense::Reversed) {
                nodes.reverse();
            }
            for pair in nodes.windows(2) {
                dec.charge_work(1)?;
                let a = pair[0].0;
                let b = pair[1].0;
                dec.charge_work((uses.len() as u64).saturating_mul(2))?;
                let outgoing = uses.iter().filter(|t| ends(&candidate, t)[0] == a).count();
                let incoming = uses.iter().filter(|t| ends(&candidate, t)[1] == b).count();
                if outgoing != 0 || incoming != 0 {
                    continue;
                }
                let mut ordered = [pair[0].clone(), pair[1].clone()];
                ordered.sort_by(|a, b| a.1[1].cmp(&b.1[1]));
                let seam = (
                    if u == q(4) { q(0) } else { u.clone() },
                    ordered[0].0,
                    ordered[1].0,
                );
                let edge = if let Some(index) = seams.get(&seam) {
                    *index
                } else {
                    let pole = seam.0.to_f64_nearest() as usize * 2;
                    let controls = source.source().controls()[pole * 2..pole * 2 + 2].to_vec();
                    let weights = source.source().weights()[pole * 2..pole * 2 + 2].to_vec();
                    let index = candidate.edges.len();
                    candidate.edges.push(NumericalEdge {
                        key: format!("{key}/seam-{}-{}-{}", seam.0, seam.1, seam.2),
                        curve: Curve::BSpline(BSplineCurve::new(
                            1,
                            vec![q(0), q(0), q(1), q(1)],
                            controls,
                            weights,
                            dec,
                        )?),
                        vertices: [ordered[0].0, ordered[1].0],
                        parameters: [ordered[0].1[1].clone(), ordered[1].1[1].clone()],
                    });
                    seams.insert(seam, index);
                    index
                };
                let native_curve = BSplineCurve::new(
                    1,
                    vec![q(0), q(0), q(1), q(1)],
                    vec![
                        Point3::new(u.clone().into(), Real::zero(), Real::zero()),
                        Point3::new(u.clone().into(), Real::one(), Real::zero()),
                    ],
                    vec![Real::one(); 2],
                    dec,
                )?;
                curves.insert(edge, native_curve);
                uses.push(NumericalTrim {
                    tag: catallax_topology::numerical::CopyTag::Independent,
                    edge,
                    sense: if a == ordered[0].0 {
                        Sense::Same
                    } else {
                        Sense::Reversed
                    },
                    uv: [pair[0].1.clone(), pair[1].1.clone()],
                });
            }
        }
        let loops = walk(&candidate, uses, dec)?;
        let mut outers = Vec::new();
        let mut holes = Vec::new();
        for boundary in loops {
            let polygon = polygon(&candidate, &boundary, &curves, dec)?;
            dec.charge_work(polygon.len() as u64)?;
            let area = polygon
                .iter()
                .zip(polygon.iter().cycle().skip(1))
                .take(polygon.len())
                .fold(q(0), |sum, (a, b)| sum + &a[0] * &b[1] - &a[1] * &b[0]);
            if area.sign() == Sign::Zero {
                return Err(bad("a cylinder loop has zero sampled area"));
            }
            if (area.sign() == Sign::Positive) == (sense == Sense::Same) {
                outers.push((boundary, polygon, Vec::new()));
            } else {
                holes.push((boundary, polygon));
            }
        }
        for (boundary, polygon) in holes {
            let mut owners = Vec::new();
            for (index, (_, outer, _)) in outers.iter().enumerate() {
                dec.charge_work(1)?;
                if inside(&polygon[0], outer, dec)? {
                    owners.push(index);
                }
            }
            let [owner] = owners.as_slice() else {
                return Err(bad("a cylinder hole has no unique enclosing region"));
            };
            outers[*owner].2.push(boundary);
        }
        for (region, (outer, _, holes)) in outers.into_iter().enumerate() {
            let face_key = format!("{key}/quarter-{quarter}/region-{region}");
            let loops: Vec<_> = std::iter::once(outer).chain(holes).collect();
            for (boundary, row) in loops.iter().enumerate() {
                for (use_index, trim) in row.iter().enumerate() {
                    added.push(NumericalNativeTrim {
                        tag: catallax_topology::numerical::CopyTag::Independent,
                        face: face_key.clone(),
                        boundary,
                        use_index,
                        edge: candidate.edges[trim.edge].key.clone(),
                        curve: curves[&trim.edge].clone(),
                    });
                }
            }
            faces.push(candidate.faces.len());
            candidate.faces.push(NumericalFace {
                key: face_key,
                surface: source.quarters()[quarter].surface().clone(),
                sense,
                loops,
            });
        }
    }
    for face in &candidate.faces[definition.faces.len()..] {
        dec.charge_work(definition.faces.len() as u64)?;
        if definition.faces.iter().any(|old| old.key == face.key) {
            return Err(bad("cylinder face key is already in use"));
        }
    }
    *definition = candidate;
    native.extend(added);
    Ok(faces)
}

fn walk(
    definition: &NumericalDefinition,
    uses: Vec<NumericalTrim>,
    dec: &mut Decider,
) -> Result<Vec<Vec<NumericalTrim>>, Refusal> {
    let mut outgoing = BTreeMap::new();
    let mut incoming = BTreeSet::new();
    for (i, trim) in uses.iter().enumerate() {
        dec.charge_work(1)?;
        let [a, b] = ends(definition, trim);
        if outgoing.insert(a, i).is_some() || !incoming.insert(b) {
            return Err(bad(format!("cylinder boundary branches at {a}/{b}")));
        }
    }
    let mut loops = Vec::new();
    while let Some((&start, _)) = outgoing.first_key_value() {
        let mut current = start;
        let mut row = Vec::new();
        while let Some(index) = outgoing.remove(&current) {
            dec.charge_work(1)?;
            row.push(uses[index].clone());
            current = ends(definition, &uses[index])[1];
        }
        if current != start {
            return Err(bad(format!(
                "cylinder boundary remains open at vertex {current}"
            )));
        }
        loops.push(row);
    }
    Ok(loops)
}
fn polygon(
    definition: &NumericalDefinition,
    boundary: &[NumericalTrim],
    curves: &BTreeMap<usize, BSplineCurve>,
    dec: &mut Decider,
) -> Result<Vec<[Rational; 2]>, Refusal> {
    let mut result = Vec::new();
    for trim in boundary {
        let mut domain = definition.edges[trim.edge].parameters.clone();
        if trim.sense == Sense::Reversed {
            domain.reverse()
        }
        for i in 0..16 {
            result.push(at(
                &curves[&trim.edge],
                &(domain[0].clone()
                    + (&domain[1] - &domain[0]) * Rational::from_ratio(i, 16).unwrap()),
                dec,
            )?);
        }
    }
    Ok(result)
}
fn inside(
    point: &[Rational; 2],
    polygon: &[[Rational; 2]],
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let mut parity = false;
    for (a, b) in polygon
        .iter()
        .zip(polygon.iter().cycle().skip(1))
        .take(polygon.len())
    {
        dec.charge_work(1)?;
        if (a[1] > point[1]) != (b[1] > point[1])
            && point[0] < a[0].clone() + (&point[1] - &a[1]) * (&b[0] - &a[0]) / (&b[1] - &a[1])
        {
            parity = !parity
        }
    }
    Ok(parity)
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::with_operation_budget;

    #[test]
    fn polygon_ownership_scan_obeys_the_work_limit() {
        let polygon = vec![[q(0), q(0)]; 10_000];
        let (result, work) = with_operation_budget(16, || {
            inside(&[q(1), q(1)], &polygon, &mut Decider::default())
        });
        assert!(work.exhausted, "unmetered polygon scan: {work:?}");
        assert_eq!(result.unwrap_err().code, "numbers.operation-budget");
    }
}
