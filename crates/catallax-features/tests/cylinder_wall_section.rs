//! The actual E07 wall section shares a retained cone/tool ring vertex.

use catallax_boolean::numerical::clip_cone_cylinder_face;
use catallax_features::{evaluate, hole::HoleTool};
use catallax_intersect::numerical_native_cylinder_intersection;
use catallax_numbers::{Ball, Decider, Rational, Real, with_operation_budget};
use catallax_recipe::{HoleDepth, HolePlacement, HoleTip, Recipe};
use catallax_shapes::{BSplineCurve, Curve, Surface};
use catallax_topology::{
    Loop,
    numerical::{NumericalEdge, NumericalOptions, NumericalVertex},
};
use serde_json::json;
use std::collections::BTreeSet;

fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}

fn source() -> Recipe {
    let side = |segment: &str| json!({"step":"ring","role":{"side":{"segment":segment}}});
    Recipe::parse(&json!({"format":"catallax-recipe/3","name":"E07 collar skew prefix","steps":[
        {"id":"section","op":"sketch","plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},"loops":[[
            {"kind":"line","id":"bottom","to":["12","0"]},
            {"kind":"line","id":"right","to":["12","12"]},
            {"kind":"line","id":"top","to":["6","12"]},
            {"kind":"line","id":"left","to":["6","0"]}]]},
        {"id":"ring","op":"revolve","profile":"section","axis":{"point":["0","0","0"],"direction":["0","1","0"]},"angle":"360"},
        {"id":"collar","op":"chamfer","body":"ring","distance":"1","edges":[
            {"between":[side("right"),side("bottom")]},{"between":[side("right"),side("top")]}]}
    ],"result":"collar"}).to_string()).unwrap()
}

fn at(curve: &BSplineCurve, parameter: &Rational, dec: &mut Decider) -> [Rational; 3] {
    curve
        .evaluate(&Ball::from_rational(parameter, 128), dec)
        .unwrap()
        .map(|v| v.mid())
}

fn distance(a: &[Rational; 3], b: &[Rational; 3]) -> Rational {
    a.iter()
        .zip(b)
        .fold(Rational::zero(), |sum, (a, b)| sum + (a - b).abs())
}

#[test]
fn skew_wall_section_joins_the_original_cone_ring_after_the_diameter_edit() {
    for diameter in ["2", "3/2"] {
        let ((), work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let built = evaluate(&source(), &mut dec).unwrap();
            let faces = built.model.faces_of_body(built.result);
            let wall = *faces.iter().find(|id| matches!(built.model.geom().surface(built.model.face(**id).surface),Surface::Cylinder(c) if c.radius_sq.as_rational()==Some(q("144")))).unwrap();
            let upper = *faces.iter().find(|id| matches!(built.model.geom().surface(built.model.face(**id).surface),Surface::Cone(c) if c.apex.y.to_f64()>12.)).unwrap();
            let Surface::Cylinder(stock) =
                built.model.geom().surface(built.model.face(wall).surface)
            else {
                unreachable!()
            };
            let axis = stock
                .axis
                .scaled(&(&Real::one() / &stock.axis.norm_squared().sqrt()));
            let mut heights = Vec::new();
            for boundary in &built.model.face(wall).loops {
                let Loop::Coedges(uses) = built.model.loop_(*boundary) else {
                    panic!("retained circular wall bounds")
                };
                for use_id in uses {
                    let edge = built.model.edge(built.model.coedge(*use_id).edge);
                    if let Curve::Circle(circle) = built.model.geom().curve(edge.curve) {
                        heights.push(
                            dec.enclose_at(&circle.center.minus(&stock.origin).dot(&axis), 128)
                                .unwrap()
                                .mid(),
                        );
                    }
                }
            }
            heights.sort();
            heights.dedup();
            let heights: [Rational; 2] = heights.try_into().unwrap();
            let placement: HolePlacement = serde_json::from_value(
                json!({"id":"radial-0","point":["21/2","12","0"],"direction":["1","-1","1"]}),
            )
            .unwrap();
            let span = HoleTool::new(
                &placement,
                &q(diameter),
                &HoleDepth::Through,
                &HoleTip::Flat,
                None,
                &q("26"),
                &mut dec,
            )
            .unwrap()
            .cylinder_span()
            .unwrap();
            assert_eq!(span.interval_mm, [q("-99/2"), q("99/2")]);
            let mut cone = clip_cone_cylinder_face(
                &built.model,
                upper,
                "upper",
                &span.cylinder,
                span.interval_mm.clone(),
                NumericalOptions {
                    tolerance_mm: q("1e-7"),
                    vertex_tolerance_mm: None,
                },
                &mut dec,
            )
            .unwrap()
            .unwrap();
            let wall = numerical_native_cylinder_intersection(
                [stock, &span.cylinder],
                [heights, span.interval_mm.clone()],
                &mut dec,
            )
            .unwrap();
            assert!(!wall.pieces().is_empty());
            assert_eq!(wall.evidence().tolerance_mm(), &q("1e-7"));
            assert!(
                Surface::Cylinder(wall.sources()[1].original().clone())
                    .same_definition(&Surface::Cylinder(span.cylinder.clone()))
            );
            let ring = cone
                .rings
                .iter()
                .position(|ring| ring.circle.radius_sq.as_rational() == Some(q("144")))
                .unwrap();
            let ring_vertices: BTreeSet<_> = cone
                .ring_links
                .iter()
                .filter(|link| link.ring == ring)
                .flat_map(|link| cone.definition.edges[link.edge].vertices)
                .collect();
            let mut joined = None;
            'events: for event in &cone.events {
                let Some(vertex) = event.vertex.filter(|v| ring_vertices.contains(v)) else {
                    continue;
                };
                let position = &cone.definition.vertices[vertex].position;
                for piece in wall.pieces() {
                    for (end, parameter) in piece.curve().domain().into_iter().enumerate() {
                        let endpoint = at(piece.curve(), parameter, &mut dec);
                        if distance(position, &endpoint) > q("1e-7") {
                            continue;
                        }
                        let uv = piece
                            .uv()
                            .each_ref()
                            .map(|curve| at(curve, parameter, &mut dec));
                        let same_tool_map = event.hits.iter().any(|hit| {
                            (uv[1][0].clone() - hit.uv[1][0].clone()).abs()
                                + (uv[1][1].clone() - hit.uv[1][1].clone()).abs()
                                <= q("1e-7")
                        });
                        if !same_tool_map {
                            continue;
                        }
                        for (source, native) in wall.sources().iter().zip(&uv) {
                            let composed = source
                                .source()
                                .evaluate(
                                    &Ball::from_rational(&native[0], 128),
                                    &Ball::from_rational(&native[1], 128),
                                    &dec,
                                )
                                .unwrap()
                                .map(|b| b.mid());
                            assert!(distance(position, &composed) <= q("1e-7"));
                        }
                        joined = Some((vertex, piece.clone(), end));
                        break 'events;
                    }
                }
            }
            let (shared, piece, end) =
                joined.expect("one actual cone/tool-to-wall/tool ring join with both native maps");
            assert!(
                cone.definition
                    .edges
                    .iter()
                    .any(|edge| edge.vertices.contains(&shared))
            );
            let other_parameter = piece.curve().domain()[1 - end];
            let other = cone.definition.vertices.len();
            cone.definition.vertices.push(NumericalVertex {
                key: "wall/other-end".into(),
                position: at(piece.curve(), other_parameter, &mut dec),
            });
            let mut vertices = [other, other];
            vertices[end] = shared;
            cone.definition.edges.push(NumericalEdge {
                key: "wall/coupled-section".into(),
                curve: Curve::BSpline(piece.curve().clone()),
                vertices,
                parameters: piece.curve().domain().map(Clone::clone),
            });
            assert_eq!(cone.definition.edges.last().unwrap().vertices[end], shared);
            println!(
                "CYLINDER_SECTION diameter={diameter} pieces={} work={:?} shared_vertex={shared}",
                wall.pieces().len(),
                wall.work_stages()
            );
        });
        assert!(!work.exhausted, "diameter={diameter}: {work:?}");
        println!("CYLINDER_JOIN diameter={diameter} total_work={}", work.used);
    }
}
