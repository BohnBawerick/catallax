//! Checked normal-strip partitions of finite polynomial extrusion faces.

use catallax_numbers::{Decider, Outcome, Real, Refusal, Sign, StrictSign, ZeroProof};
use catallax_shapes::{Curve, SplineSupport, Surface};
use catallax_topology::{BezierExtrusionBand, EdgeId, Face, FaceId, Loop, LoopId, Model, Sense};
use std::collections::{BTreeMap, BTreeSet};

fn coverage(detail: &str) -> Refusal {
    Refusal::new("boolean.bezier-arrangement-coverage", detail)
}
fn sign(value: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    match dec.sign(value) {
        Outcome::Proven(StrictSign::Positive) => Ok(Sign::Positive),
        Outcome::Proven(StrictSign::Negative) => Ok(Sign::Negative),
        Outcome::ProvenCoincident(p) if p.proof != ZeroProof::SeparationBound => Ok(Sign::Zero),
        Outcome::ProvenCoincident(_) => Err(coverage("a strip boundary needs exact incidence")),
        Outcome::NotExact(error) => Err(Refusal::new(error.code, error.detail)),
        Outcome::Refused(error) => Err(error),
    }
}
struct Cell {
    interval: [Real; 2],
    boundary: LoopId,
}

pub(super) fn pieces(
    model: &mut Model,
    face: FaceId,
    loops: &[LoopId],
    boundary: &BTreeSet<EdgeId>,
    sections: &[(EdgeId, Sense)],
    vertices: &[LoopId],
    dec: &mut Decider,
) -> Result<Option<Vec<FaceId>>, Refusal> {
    let original = model.face(face).clone();
    let Surface::BSpline(SplineSupport::BezierExtrusion(support)) =
        model.geom().surface(original.surface)
    else {
        return Ok(None);
    };
    let support = support.clone();
    if !vertices.is_empty() {
        return Err(Refusal::new(
            "boolean.bezier-vertex-trim",
            "normal Bezier strip arrangements do not accept isolated vertex trims",
        ));
    }
    let checked = BezierExtrusionBand::of_face(model, face, dec)?
        .ok_or_else(|| coverage("the original spline face has no checked band"))?;
    let original_interval = checked.axial().map(Clone::clone);
    let mut interior = BTreeMap::new();
    for &(edge, _) in sections {
        if boundary.contains(&edge) {
            continue;
        }
        let Curve::BSpline(curve) = model.geom().curve(model.edge(edge).curve) else {
            return Err(Refusal::new(
                "boolean.bezier-strip-sections",
                "Bezier strip imprinting requires complete normal polynomial sections",
            ));
        };
        let (height, _) = support
            .section_parameter(curve, dec)?
            .ok_or_else(|| coverage("a section is not an exact complete rail of this carrier"))?;
        if sign(&(&height - &original_interval[0]), dec)? != Sign::Positive
            || sign(&(&height - &original_interval[1]), dec)? != Sign::Negative
        {
            return Err(coverage(
                "a non-boundary section does not lie strictly within the original band",
            ));
        }
        interior.insert(edge, height);
    }
    let mut cells: Vec<Cell> = Vec::new();
    for &boundary_loop in loops {
        let temporary = model.add_face(Face {
            surface: original.surface,
            sense: original.sense,
            loops: vec![boundary_loop],
        });
        let band = BezierExtrusionBand::of_face(model, temporary, dec)?
            .ok_or_else(|| coverage("an arrangement walk has no checked strip domain"))?;
        let interval = band.axial().map(Clone::clone);
        let inside = sign(&(&interval[0] - &original_interval[0]), dec)? != Sign::Negative
            && sign(&(&interval[1] - &original_interval[1]), dec)? != Sign::Positive;
        if !inside {
            if sign(&(&interval[1] - &original_interval[0]), dec)? != Sign::Positive
                || sign(&(&interval[0] - &original_interval[1]), dec)? != Sign::Negative
            {
                continue;
            }
            return Err(coverage(
                "an arrangement walk crosses an unsplit original strip boundary",
            ));
        }
        let mut index = 0;
        while index < cells.len()
            && sign(&(&interval[0] - &cells[index].interval[0]), dec)? != Sign::Negative
        {
            index += 1;
        }
        cells.insert(
            index,
            Cell {
                interval,
                boundary: boundary_loop,
            },
        );
    }
    let (Some(first), Some(last)) = (cells.first(), cells.last()) else {
        return Err(coverage("the arrangement does not cover the original band"));
    };
    if sign(&(&first.interval[0] - &original_interval[0]), dec)? != Sign::Zero
        || sign(&(&last.interval[1] - &original_interval[1]), dec)? != Sign::Zero
    {
        return Err(coverage("the strip partition omits an original boundary"));
    }
    for pair in cells.windows(2) {
        if sign(&(&pair[0].interval[1] - &pair[1].interval[0]), dec)? != Sign::Zero {
            return Err(coverage("the strip partition has a gap or overlap"));
        }
    }
    let mut uses: BTreeMap<EdgeId, [usize; 2]> = BTreeMap::new();
    for cell in &cells {
        let Loop::Coedges(coedges) = model.loop_(cell.boundary) else {
            return Err(coverage("a strip boundary has no coedges"));
        };
        for coedge in coedges {
            let coedge = model.coedge(*coedge);
            uses.entry(coedge.edge).or_default()[usize::from(coedge.sense == Sense::Reversed)] += 1;
        }
    }
    for edge in boundary {
        if uses
            .get(edge)
            .is_none_or(|counts| counts[0] + counts[1] != 1)
        {
            return Err(coverage("an original boundary edge is missing or repeated"));
        }
    }
    for (edge, height) in &interior {
        if uses.get(edge) != Some(&[1, 1]) {
            return Err(coverage(
                "an interior section is not sewn once to each adjacent strip",
            ));
        }
        let mut found = false;
        for pair in cells.windows(2) {
            if sign(&(&pair[0].interval[1] - height), dec)? == Sign::Zero {
                found = true;
                break;
            }
        }
        if !found {
            return Err(coverage(
                "an interior section is absent from the strip partition",
            ));
        }
    }
    if uses
        .keys()
        .any(|edge| !boundary.contains(edge) && !interior.contains_key(edge))
    {
        return Err(coverage(
            "the strip partition introduces an unproved boundary",
        ));
    }
    let origins = model.face_origins(face).to_vec();
    Ok(Some(
        cells
            .into_iter()
            .map(|cell| {
                model.add_named_face(
                    Face {
                        surface: original.surface,
                        sense: original.sense,
                        loops: vec![cell.boundary],
                    },
                    origins.clone(),
                )
            })
            .collect(),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::Rational;
    use catallax_shapes::{
        BSplineCurve, BezierExtrusion, Line, PointDef, RPoint3, RVector3, SurfaceId,
    };
    use catallax_topology::{Coedge, Edge, Vertex};

    fn surface(model: &mut Model) -> SurfaceId {
        let controls = [(0, 0), (0, 4), (4, 4), (4, 0)]
            .into_iter()
            .map(|(x, y)| RPoint3::new(x.into(), y.into(), 0.into()))
            .collect();
        let support = BezierExtrusion::new(
            controls,
            RVector3::new(0.into(), 0.into(), 2.into()),
            &mut Decider::default(),
        )
        .unwrap();
        model
            .geom_mut()
            .add_surface(Surface::BSpline(support.into()))
    }
    fn band(
        model: &mut Model,
        surface: SurfaceId,
        lo: &str,
        hi: &str,
    ) -> (FaceId, LoopId, BTreeSet<EdgeId>) {
        let Surface::BSpline(SplineSupport::BezierExtrusion(support)) =
            model.geom().surface(surface)
        else {
            panic!("polynomial carrier expected")
        };
        let support = support.clone();
        let lo = Real::from(Rational::parse(lo).unwrap());
        let hi = Real::from(Rational::parse(hi).unwrap());
        let positions = [
            support.controls()[0].plus(&support.direction().scaled(&lo)),
            support.controls()[3].plus(&support.direction().scaled(&lo)),
            support.controls()[3].plus(&support.direction().scaled(&hi)),
            support.controls()[0].plus(&support.direction().scaled(&hi)),
        ];
        let vertices = positions.clone().map(|point| {
            let point = model.geom_mut().add_point(PointDef::Explicit(point));
            model.add_vertex(Vertex { point })
        });
        let row = |height: &Real| {
            BSplineCurve::new(
                3,
                support.lower_edge().knots().to_vec(),
                support
                    .controls()
                    .iter()
                    .map(|p| p.plus(&support.direction().scaled(height)))
                    .collect(),
                vec![Real::one(); 4],
                &mut Decider::default(),
            )
            .unwrap()
        };
        let curves = [
            Curve::BSpline(row(&lo)),
            Curve::Line(
                Line::new(
                    positions[1].clone(),
                    positions[2].minus(&positions[1]),
                    &mut Decider::default(),
                )
                .unwrap(),
            ),
            Curve::BSpline(row(&hi)),
            Curve::Line(
                Line::new(
                    positions[0].clone(),
                    positions[3].minus(&positions[0]),
                    &mut Decider::default(),
                )
                .unwrap(),
            ),
        ]
        .map(|curve| model.geom_mut().add_curve(curve));
        let edges: [EdgeId; 4] = std::array::from_fn(|index| {
            let (from, to) = [(0, 1), (1, 2), (3, 2), (0, 3)][index];
            model.add_edge(Edge {
                curve: curves[index],
                sense: Sense::Same,
                start: Some(vertices[from]),
                end: Some(vertices[to]),
            })
        });
        let coedges = edges
            .into_iter()
            .enumerate()
            .map(|(index, edge)| {
                model.add_coedge(Coedge {
                    edge,
                    sense: if index < 2 {
                        Sense::Same
                    } else {
                        Sense::Reversed
                    },
                })
            })
            .collect();
        let boundary_loop = model.add_loop(Loop::Coedges(coedges));
        let face = model.add_face(Face {
            surface,
            sense: Sense::Same,
            loops: vec![boundary_loop],
        });
        (face, boundary_loop, edges.into_iter().collect())
    }
    #[test]
    fn complete_original_domains_are_retained_but_gaps_and_overlaps_refuse() {
        let mut model = Model::new();
        let surface = surface(&mut model);
        let (face, original, boundary) = band(&mut model, surface, "0", "1");
        let result = pieces(
            &mut model,
            face,
            &[original],
            &boundary,
            &[],
            &[],
            &mut Decider::default(),
        )
        .unwrap()
        .unwrap();
        assert_eq!(result.len(), 1);
        let checked = BezierExtrusionBand::of_face(&model, result[0], &mut Decider::default())
            .unwrap()
            .unwrap();
        assert!(checked.axial()[0].rational_ref().unwrap().is_zero());
        assert_eq!(checked.axial()[1].rational_ref().unwrap(), &Rational::one());
        for ((alo, ahi), (blo, bhi)) in [(("0", "1/4"), ("3/4", "1")), (("0", "3/4"), ("1/2", "1"))]
        {
            let (_, a, _) = band(&mut model, surface, alo, ahi);
            let (_, b, _) = band(&mut model, surface, blo, bhi);
            assert_eq!(
                pieces(
                    &mut model,
                    face,
                    &[a, b],
                    &boundary,
                    &[],
                    &[],
                    &mut Decider::default()
                )
                .unwrap_err()
                .code,
                "boolean.bezier-arrangement-coverage"
            );
        }
    }
    #[test]
    fn an_omitted_interior_section_cannot_pass_a_complete_outer_loop() {
        let mut model = Model::new();
        let surface = surface(&mut model);
        let (face, original, boundary) = band(&mut model, surface, "0", "1");
        let (_, middle, _) = band(&mut model, surface, "1/2", "3/4");
        let Loop::Coedges(coedges) = model.loop_(middle) else {
            panic!("coedges expected")
        };
        let section = model.coedge(coedges[0]).edge;
        assert_eq!(
            pieces(
                &mut model,
                face,
                &[original],
                &boundary,
                &[(section, Sense::Same)],
                &[],
                &mut Decider::default()
            )
            .unwrap_err()
            .code,
            "boolean.bezier-arrangement-coverage"
        );
    }
}
