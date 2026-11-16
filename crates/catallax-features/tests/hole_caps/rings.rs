//! The existing skew component's sliver-side original-ring anchor.

use super::*;
use catallax_boolean::numerical::NumericalCylinderSource;
use catallax_features::hole::cylinder_ring_curve;
use catallax_intersect::numerical_native_cylinder_intersection;
use catallax_numbers::{NumericalEvidence, transfer};
use catallax_shapes::{BSplineCurve, Curve, Surface};
use catallax_topology::Loop;
use std::collections::BTreeSet;

fn at(curve: &BSplineCurve, t: &Rational, dec: &mut Decider) -> [Rational; 3] {
    let point = curve.point_at(&Real::from(t.clone()), dec).unwrap();
    [&point.x, &point.y, &point.z].map(|v| dec.enclose_budgeted_at(v, 128).unwrap().mid())
}

fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.into_iter()
        .zip(b)
        .map(|(a, b)| (a - b).powi(2))
        .sum::<f64>()
        .sqrt()
}

fn on(source: &BSplineSurface, uv: &[Rational; 3], dec: &Decider) -> [f64; 3] {
    floats(
        source
            .evaluate(
                &Ball::from_rational(&uv[0], 128),
                &Ball::from_rational(&uv[1], 128),
                dec,
            )
            .unwrap(),
    )
}

#[test]
fn shared_ring_anchor_keeps_original_native_maps_after_diameter_edit() {
    let (records, work) = with_operation_budget(20_000_000, || {
        let mut dec = Decider::default();
        let mut records = Vec::new();
        let mut previous_keys = None;
        for diameter in [q("2"), q("3/2")] {
            let component = super::skew::component(&diameter, &mut dec);
            assert_eq!(component.span.interval_mm, [q("-99/2"), q("99/2")]);
            assert_eq!(component.definition.faces.last().unwrap().key, "top");
            assert_eq!(
                component
                    .bindings
                    .iter()
                    .filter(|b| b.face == "top")
                    .count(),
                12
            );
            let built = &component.built;
            let faces = built.model.faces_of_body(built.result);
            let stock_face = *faces
                .iter()
                .find(|id| {
                    matches!(built.model.geom().surface(built.model.face(**id).surface),
                Surface::Cylinder(c) if c.radius_sq.as_rational()==Some(q("144")))
                })
                .unwrap();
            let Surface::Cylinder(stock) = built
                .model
                .geom()
                .surface(built.model.face(stock_face).surface)
            else {
                unreachable!()
            };
            let axis = stock
                .axis
                .scaled(&(&Real::one() / &stock.axis.norm_squared().sqrt()));
            let mut rings = BTreeSet::new();
            let mut heights = Vec::new();
            for boundary in &built.model.face(stock_face).loops {
                let Loop::Coedges(uses) = built.model.loop_(*boundary) else {
                    panic!("original circular stock boundary")
                };
                for use_id in uses {
                    let edge = built.model.coedge(*use_id).edge;
                    rings.insert(edge);
                    let Curve::Circle(circle) =
                        built.model.geom().curve(built.model.edge(edge).curve)
                    else {
                        panic!("original circular ring")
                    };
                    heights.push(
                        dec.enclose_budgeted_at(
                            &circle.center.minus(&stock.origin).dot(&axis),
                            128,
                        )
                        .unwrap()
                        .mid(),
                    );
                }
            }
            heights.sort();
            heights.dedup();
            let heights: [Rational; 2] = heights.try_into().unwrap();
            assert_eq!(heights, [q("1"), q("11")]);
            let options = NumericalOptions {
                tolerance_mm: q("1e-7"),
            };
            let stock =
                NumericalCylinderSource::new(stock, heights, options.clone(), &mut dec).unwrap();
            let part = &component.cone;
            let link = part
                .ring_links
                .iter()
                .find(|link| {
                    rings.contains(&part.rings[link.ring].edge)
                        && link.source_u[0] == q("0")
                        && link.source_u[1] < q("1/10")
                })
                .unwrap();
            let edge = &part.definition.edges[link.edge];
            let Curve::BSpline(original) = &edge.curve else {
                panic!("retained original ring")
            };
            let mapping = cylinder_ring_curve(
                &edge.key,
                original,
                edge.parameters.clone(),
                &stock,
                1,
                &options,
                &mut dec,
            )
            .unwrap()
            .unwrap();
            assert_eq!(mapping.quarter(), 3);
            assert_eq!(mapping.curve().domain(), [&q("0"), &q("1")]);
            assert_eq!(mapping.evidence().tolerance_mm(), &q("1e-7"));
            for check in mapping.evidence().checks() {
                if let NumericalCheck::Samples { name, samples, .. } = check {
                    assert_eq!(
                        samples.len(),
                        if name.ends_with("composition") { 17 } else { 2 }
                    );
                    assert!(
                        samples
                            .iter()
                            .all(|sample| sample.distance_upper_mm <= q("1e-7"))
                    );
                }
            }
            let vertex = edge.vertices[1];
            assert_eq!(
                component.definition.vertices[vertex].position,
                part.definition.vertices[vertex].position
            );
            assert_eq!(
                component.definition.vertices[vertex].key,
                part.definition.vertices[vertex].key
            );
            let event = part
                .events
                .iter()
                .find(|event| event.vertex == Some(vertex))
                .unwrap();
            let hit = event
                .hits
                .iter()
                .find(|hit| hit.uv[0][1] == q("1"))
                .unwrap();
            let section = part
                .cylinder_trims
                .iter()
                .find(|trim| {
                    trim.piece == hit.piece
                        && part.definition.edges[trim.edge].vertices.contains(&vertex)
                })
                .unwrap();
            let section_edge = &part.definition.edges[section.edge];
            let end = section_edge
                .vertices
                .iter()
                .position(|v| *v == vertex)
                .unwrap();
            let parameter = &section_edge.parameters[end];
            assert_eq!(parameter, &hit.parameter);
            let shared = &component.definition.vertices[vertex];
            assert_eq!(shared.key, part.definition.vertices[vertex].key);
            assert_eq!(shared.position, part.definition.vertices[vertex].position);
            assert_eq!(
                component.definition.edges[link.edge].vertices[1],
                component.definition.edges[section.edge].vertices[end]
            );
            let point = shared.position.each_ref().map(Rational::to_f64_nearest);
            let stock_uv = at(mapping.curve(), &edge.parameters[1], &mut dec);
            let tool_uv = at(&section.curve, parameter, &mut dec);
            let cone_uv = at(
                &part.intersection.pieces()[hit.piece].uv()[0],
                parameter,
                &mut dec,
            );
            let Surface::BSpline(stock_quarter) = stock.quarters()[mapping.quarter()].surface()
            else {
                unreachable!()
            };
            let stock_point = on(stock_quarter.carrier(), &stock_uv, &dec);
            let tool_point = on(part.cylinder_source.source(), &tool_uv, &dec);
            let cone_point = on(&part.intersection.windows().sources()[0], &cone_uv, &dec);
            for actual in [stock_point, tool_point, cone_point] {
                assert!(distance(actual, point) <= 1e-7);
            }
            assert_eq!(section.quarter, 1);
            let incoming = numerical_native_cylinder_intersection(
                [stock.original(), &component.span.cylinder],
                [
                    stock.axis_parameters_mm().clone(),
                    component.span.interval_mm.clone(),
                ],
                &mut dec,
            )
            .unwrap();
            assert!(
                incoming.sources()[0]
                    .source()
                    .same_definition(stock.source())
            );
            assert!(
                incoming.sources()[1]
                    .source()
                    .same_definition(part.cylinder_source.source())
            );
            let mut joins = Vec::new();
            for (index, piece) in incoming.pieces().iter().enumerate() {
                for (end, parameter) in piece.curve().domain().into_iter().enumerate() {
                    let endpoint = at(piece.curve(), parameter, &mut dec);
                    if distance(endpoint.each_ref().map(Rational::to_f64_nearest), point) > 1e-7 {
                        continue;
                    }
                    let native = piece
                        .uv()
                        .each_ref()
                        .map(|curve| at(curve, parameter, &mut dec));
                    assert!(
                        distance(
                            native[0].each_ref().map(Rational::to_f64_nearest),
                            stock_uv.each_ref().map(Rational::to_f64_nearest)
                        ) <= 1e-7
                    );
                    assert!(
                        distance(
                            native[1].each_ref().map(Rational::to_f64_nearest),
                            tool_uv.each_ref().map(Rational::to_f64_nearest)
                        ) <= 1e-7
                    );
                    for (source, uv) in incoming.sources().iter().zip(&native) {
                        assert!(distance(on(source.source(), uv, &dec), point) <= 1e-7);
                    }
                    joins.push((index, end));
                }
            }
            let [(piece_index, end)] = joins.as_slice() else {
                panic!("one incoming wall section at the owned original-ring anchor: {joins:?}")
            };
            let piece = &incoming.pieces()[*piece_index];
            let mut joined = component.definition.clone();
            let other = joined.vertices.len();
            joined
                .vertices
                .push(catallax_topology::numerical::NumericalVertex {
                    key: "stock-wall/other-end".into(),
                    position: at(piece.curve(), piece.curve().domain()[1 - end], &mut dec),
                });
            let mut vertices = [other, other];
            vertices[*end] = vertex;
            joined
                .edges
                .push(catallax_topology::numerical::NumericalEdge {
                    key: "stock-wall/coupled-section".into(),
                    curve: Curve::BSpline(piece.curve().clone()),
                    vertices,
                    parameters: piece.curve().domain().map(Clone::clone),
                });
            assert_eq!(joined.edges.last().unwrap().vertices[*end], vertex);
            assert_eq!(
                joined.vertices[vertex].key,
                part.definition.vertices[vertex].key
            );
            let naive_uv = [&q("4") - &edge.parameters[1], q("1"), q("0")];
            let naive_error = distance(on(stock_quarter.carrier(), &naive_uv, &dec), stock_point);
            assert!(naive_error > 1e-3);
            assert!(
                cylinder_ring_curve(
                    "wrong-stock-ring",
                    original,
                    edge.parameters.clone(),
                    &stock,
                    0,
                    &options,
                    &mut dec
                )
                .unwrap()
                .is_none()
            );
            let saved = transfer::encode(
                &(mapping.curve().clone(), mapping.evidence().clone()),
                &mut dec,
            )
            .unwrap();
            let (reopened, evidence): (BSplineCurve, NumericalEvidence) =
                transfer::decode(&saved, &mut dec).unwrap();
            assert!(reopened.same_definition(mapping.curve()));
            assert_eq!(evidence, *mapping.evidence());
            assert!(
                distance(
                    on(
                        stock_quarter.carrier(),
                        &at(&reopened, &edge.parameters[1], &mut dec),
                        &dec
                    ),
                    point
                ) <= 1e-7
            );
            let keys = (
                part.definition.vertices[vertex].key.clone(),
                edge.key.clone(),
                section_edge.key.clone(),
            );
            if let Some(previous) = previous_keys {
                assert_eq!(keys, previous);
            }
            previous_keys = Some(keys);
            records.push(json!({"diameter":diameter.to_string(),"vertex_key":part.definition.vertices[vertex].key,
                "original_ring":format!("{:?}",part.rings[link.ring].edge),"point_mm":point,
                "stock_quarter":mapping.quarter(),"stock_uv":stock_uv.each_ref().map(ToString::to_string),
                "tool_quarter":section.quarter,"tool_uv":tool_uv.each_ref().map(ToString::to_string),
                "cone_uv":cone_uv.each_ref().map(ToString::to_string),"naive_error_mm":naive_error,
                "incoming_piece":piece_index,"incoming_end":end,
                "incoming_pieces":incoming.pieces().iter().map(|piece| json!({
                    "domains":piece.source_domains().clone().map(|domain|domain.map(|axis|axis.map(|v|v.to_string()))),
                    "ends":piece.curve().domain().map(|t|at(piece.curve(),t,&mut dec).map(|v|v.to_f64_nearest())),
                    "uv":piece.uv().each_ref().map(|curve|piece.curve().domain().map(|t|at(curve,t,&mut dec).map(|v|v.to_f64_nearest())))
                })).collect::<Vec<_>>(),
                "status":"incoming cylinder/cylinder section shares the original cone/tool/ring vertex and both native maps"}));
        }
        records
    });
    assert!(!work.exhausted);
    if let Some(path) = std::env::var_os("CATALLAX_RING_JOIN_PROOF") {
        let path = std::path::PathBuf::from(path);
        std::fs::create_dir_all(&path).unwrap();
        std::fs::write(
            path.join("anchor.json"),
            serde_json::to_string_pretty(&json!({"work":work.used,"records":records})).unwrap(),
        )
        .unwrap();
    }
    println!("one shared ring anchor work={}", work.used);
}
