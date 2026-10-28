use super::*;
use catallax_numbers::{Rational, Real, StrictSign, with_operation_budget};
use catallax_recipe::{Recipe, Role};
use catallax_shapes::{Circle, Ellipse, Line, Vector3};
use serde_json::{Value, json};

const FIXTURE: &str = include_str!("../../tests/fixtures/trimmed-curves.json");

fn recipe(scale: i64, clockwise: bool) -> Recipe {
    let mut value: Value = serde_json::from_str(FIXTURE).unwrap();
    let y = if clockwise { -4 } else { 4 } * scale;
    value["steps"][0]["loops"][0][0]["to"] = json!([(3 * scale).to_string(), y.to_string()]);
    value["steps"][0]["loops"][0][1]["to"] = json!([(3 * scale).to_string(), (-y).to_string()]);
    value["steps"][0]["loops"][0][0]["turn"] = json!(if clockwise { "cw" } else { "ccw" });
    Recipe::parse(&value.to_string()).unwrap()
}

fn build(recipe: &Recipe) -> ApplicationBody {
    let result = crate::build(recipe, &Default::default(), &mut Decider::default()).unwrap();
    assert!(result.built.model.check(result.built.result).is_empty());
    assert!(result.display.marks.is_empty());
    result
}

fn top_arc(body: &ApplicationBody) -> EdgeCurve {
    body.names
        .edges
        .keys()
        .find_map(|edge| {
            let query = body.edge_curve(*edge, &mut Decider::default()).unwrap();
            (query.origins.iter().any(|o| o.role == Role::End)
                && matches!(query.span, CurveSpan::CircularArc { .. }))
            .then_some(query)
        })
        .unwrap()
}

fn point(p: &Point3, expected: [i64; 3]) {
    for (value, expected) in [&p.x, &p.y, &p.z].into_iter().zip(expected) {
        assert!(matches!(
            Decider::default().sign(&(value - &Real::from_i64(expected))),
            Outcome::ProvenCoincident(_)
        ));
    }
}

fn circle_radius(query: &EdgeCurve, expected: i64) {
    let CurveForm::Circle {
        center,
        normal,
        radius_sq,
        ref_dir,
    } = &query.form
    else {
        panic!("circle")
    };
    assert!(matches!(
        Decider::default().sign(&(radius_sq - &Real::from_i64(expected * expected))),
        Outcome::ProvenCoincident(_)
    ));
    let carrier = Curve::Circle(Circle {
        center: center.clone(),
        normal: normal.clone(),
        radius_sq: radius_sq.clone(),
        ref_dir: ref_dir.clone(),
    });
    assert!(query.definition.same_definition(&carrier));
}

#[test]
fn public_application_queries_keep_exact_finite_definitions_without_display() {
    let input = recipe(1, false);
    let mut body = build(&input);
    let saved = top_arc(&body);
    // Exact picks must survive even when a consumer discards all display data.
    body.display.value.positions.clear();
    body.display.value.triangles.clear();
    body.display.value.edges.clear();
    for &edge in body.names.edges.keys() {
        let query = body.edge_curve(edge, &mut Decider::default()).unwrap();
        let stored = body.built.model.edge(edge);
        assert_eq!(query.curve, stored.curve);
        assert_eq!(query.sense, stored.sense);
        assert_eq!(query.name, body.names.edges[&edge]);
        assert_eq!(
            query.constructions,
            body.built.model.edge_constructions(edge)
        );
        assert!(
            query
                .definition
                .same_definition(body.built.model.geom().curve(stored.curve))
        );
        let CurveTrim::Finite { start, end } = &query.trim else {
            panic!("finite extrusion")
        };
        assert_eq!(Some(start.vertex), stored.start);
        assert_eq!(Some(end.vertex), stored.end);
        for endpoint in [start, end] {
            let retained = body.built.model.geom().position(endpoint.point);
            assert!(endpoint.position.x.same_node(&retained.x));
            assert!(endpoint.position.y.same_node(&retained.y));
            assert!(endpoint.position.z.same_node(&retained.z));
        }
        assert!(matches!(
            (&query.form, query.span),
            (CurveForm::Line { .. }, CurveSpan::Line)
                | (CurveForm::Circle { .. }, CurveSpan::CircularArc { .. })
        ));
    }
    let query = body
        .edge_curve(saved.edge, &mut Decider::default())
        .unwrap();
    circle_radius(&query, 5);
    assert_eq!(body.built.recipe(), &input);
    let CurveTrim::Finite { start, end } = query.trim else {
        panic!()
    };
    point(&start.position, [3, -4, 2]);
    point(&end.position, [3, 4, 2]);
    let CurveSpan::CircularArc {
        interior,
        crosses_seam,
    } = query.span
    else {
        panic!()
    };
    assert!(
        !crosses_seam,
        "an endpoint on the seam is not an interior crossing"
    );
    assert!(matches!(
        Decider::default().sign(&(interior.x - Real::from_i64(3))),
        Outcome::Proven(StrictSign::Positive)
    ));
}

#[test]
fn reversed_and_complementary_arcs_preserve_the_selected_branch() {
    let forward = top_arc(&build(&recipe(1, false)));
    let mut reversed_body = build(&recipe(1, false));
    reverse_edges(&mut reversed_body);
    let reverse = top_arc(&reversed_body);
    assert_eq!(forward.sense, Sense::Same);
    assert_eq!(reverse.sense, Sense::Reversed);
    let CurveTrim::Finite { start, end } = &reverse.trim else {
        panic!()
    };
    point(&start.position, [3, 4, 2]);
    point(&end.position, [3, -4, 2]);
    assert!(matches!(
        reverse.span,
        CurveSpan::CircularArc {
            crosses_seam: false,
            ..
        }
    ));
    check_step_edge(&reversed_body, &reverse, 5);
    for &edge in reversed_body.names.edges.keys() {
        assert_eq!(
            reversed_body
                .edge_curve(edge, &mut Decider::default())
                .unwrap()
                .sense,
            Sense::Reversed
        );
    }

    let mut input: Value = serde_json::from_str(FIXTURE).unwrap();
    input["steps"][0]["loops"][0][0]["turn"] = json!("cw");
    let major = top_arc(&build(&Recipe::parse(&input.to_string()).unwrap()));
    let CurveSpan::CircularArc {
        interior,
        crosses_seam,
    } = major.span
    else {
        panic!()
    };
    assert!(!crosses_seam);
    // The complementary major arc lies on the opposite side of the chord.
    assert!(matches!(
        Decider::default().sign(&(interior.x - Real::from_i64(3))),
        Outcome::Proven(StrictSign::Negative)
    ));
}

#[test]
fn edit_resolves_stable_origin_then_matches_retained_definition_and_step() {
    let old = build(&recipe(1, false));
    let selected = top_arc(&old);
    let name = selected.name.clone().unwrap();
    for (scale, clockwise) in [(1, false), (2, false), (2, true)] {
        let input = recipe(scale, clockwise);
        let body = build(&input);
        let ids = body.built.model.edges_with_name(body.built.result, &name);
        assert_eq!(ids.len(), 1);
        let query = body.edge_curve(ids[0], &mut Decider::default()).unwrap();
        assert_eq!(query.name.as_ref(), Some(&name));
        assert_eq!(query.origins, selected.origins);
        circle_radius(&query, scale * 5);
        assert_eq!(body.built.recipe(), &input);
        check_step_edge(&body, &query, scale * 5);
    }
    // A retained query continues to describe its own revision after edits.
    circle_radius(&selected, 5);
    assert_eq!(old.built.recipe(), &recipe(1, false));
}

fn check_step_edge(body: &ApplicationBody, query: &EdgeCurve, radius: i64) {
    let options = catallax_step::StepOptions::default();
    let mut bytes = Vec::new();
    let report = catallax_step::write_step(
        &body.built.model,
        &[body.built.result],
        &options,
        &mut Decider::default(),
        &mut bytes,
    )
    .unwrap();
    assert_eq!(report.value.approximated_curves, 0);
    // Exact 3D circles still have certified auxiliary-pcurve fitting and rounding bounds.
    assert!(report.value.max_deviation <= options.curve_tolerance);
    let text = String::from_utf8(bytes).unwrap();
    let entities: std::collections::BTreeMap<_, _> = text
        .lines()
        .filter_map(|line| {
            let (id, value) = line.strip_prefix('#')?.split_once('=')?;
            Some((id.parse::<usize>().unwrap(), value))
        })
        .collect();
    let refs = |s: &str| -> Vec<usize> {
        s.split('#')
            .skip(1)
            .map(|s| {
                s.split(|c: char| !c.is_ascii_digit())
                    .next()
                    .unwrap()
                    .parse()
                    .unwrap()
            })
            .collect()
    };
    let vertex_point = |id| {
        let id = refs(entities[&id])[0];
        let s = entities[&id]
            .strip_prefix("CARTESIAN_POINT('',(")
            .unwrap()
            .split(')')
            .next()
            .unwrap();
        s.split(',')
            .map(|v| Rational::parse(v.trim_end_matches('.')).unwrap())
            .collect::<Vec<_>>()
    };
    let CurveTrim::Finite { start, end } = &query.trim else {
        panic!()
    };
    let coords = |p: &Point3| {
        [&p.x, &p.y, &p.z]
            .map(|v| v.as_rational().unwrap())
            .to_vec()
    };
    let start = coords(&start.position);
    let end = coords(&end.position);
    assert!(
        entities
            .values()
            .filter(|v| v.starts_with("EDGE_CURVE("))
            .any(|entity| {
                let ids = refs(entity);
                let curve = entities[&ids[2]];
                let curve =
                    if curve.starts_with("SURFACE_CURVE(") || curve.starts_with("SEAM_CURVE(") {
                        entities[&refs(curve)[0]]
                    } else {
                        curve
                    };
                vertex_point(ids[0]) == start
                    && vertex_point(ids[1]) == end
                    && curve.starts_with("CIRCLE(")
                    && curve.ends_with(&format!(",{radius}.);"))
                    && entity.ends_with(if query.sense == Sense::Same {
                        ",.T.);"
                    } else {
                        ",.F.);"
                    })
            }),
        "selected exact endpoints, radius and sense must match one exported circle edge"
    );
}

#[test]
fn full_circles_are_explicit_and_incremental_picks_use_the_checkpoint() {
    let input = Recipe::parse(r#"{"format":"catallax-recipe/3","name":"ring","steps":[{"id":"stock","op":"cylinder","radius":"5","height":"2","pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}}],"result":"stock"}"#).unwrap();
    let body = build(&input);
    for &edge in body.names.edges.keys() {
        let query = body.edge_curve(edge, &mut Decider::default()).unwrap();
        assert!(matches!(query.span, CurveSpan::FullCircle));
        assert!(matches!(query.trim, CurveTrim::Closed));
        circle_radius(&query, 5);
    }
    let checkpoint = crate::build_incremental(
        &recipe(1, false),
        None,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let edited = crate::build_incremental(
        &recipe(2, true),
        Some(checkpoint.checkpoint.clone()),
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap();
    for &edge in edited.names.edges.keys() {
        let query = edited.edge_curve(edge, &mut Decider::default()).unwrap();
        if matches!(query.span, CurveSpan::CircularArc { .. }) {
            circle_radius(&query, 10);
        }
    }
}

#[test]
fn seam_boolean_keeps_algebraic_endpoints_and_complete_origin_sets() {
    let body = build(&Recipe::parse(include_str!("../../../../parts/P0/seam.json")).unwrap());
    let mut arcs = 0;
    let mut seam_crossings = 0;
    for &edge in body.names.edges.keys() {
        let query = body.edge_curve(edge, &mut Decider::default()).unwrap();
        if let CurveSpan::CircularArc {
            interior,
            crosses_seam,
        } = query.span
        {
            arcs += 1;
            seam_crossings += usize::from(crosses_seam);
            let CurveForm::Circle {
                center,
                normal,
                radius_sq,
                ..
            } = query.form
            else {
                panic!()
            };
            assert!(matches!(
                Decider::default().sign(&(interior.minus(&center).norm_squared() - radius_sq)),
                Outcome::ProvenCoincident(_)
            ));
            assert!(matches!(
                Decider::default().sign(&interior.minus(&center).dot(&normal)),
                Outcome::ProvenCoincident(_)
            ));
            let CurveTrim::Finite { start, end } = query.trim else {
                panic!()
            };
            assert!(start.position.as_rational().is_none());
            assert!(end.position.as_rational().is_none());
            assert!(!query.origins.is_empty());
            for step in ["a", "b"] {
                assert!(
                    query
                        .origins
                        .iter()
                        .any(|origin| origin.step.as_str() == step),
                    "the merged cap retains both cylinders' origins"
                );
            }
        }
    }
    assert_eq!(arcs, 4);
    assert!(seam_crossings > 0);
}

// Reversing every edge and its coedge uses is the same valid boundary with
// opposite edge storage directions. This controls orientation independently
// of profile canonicalisation, which can reverse clockwise authored loops.
fn reverse_edges(body: &mut ApplicationBody) {
    use catallax_topology::{Body, Coedge, Face, Lump, Shell};
    use std::collections::BTreeMap;
    let model = &mut body.built.model;
    let reverse = |sense| {
        if sense == Sense::Same {
            Sense::Reversed
        } else {
            Sense::Same
        }
    };
    let mut edges = BTreeMap::new();
    for old in model.edges_of_body(body.built.result) {
        let mut definition = model.edge(old).clone();
        definition.sense = reverse(definition.sense);
        std::mem::swap(&mut definition.start, &mut definition.end);
        let new = model.add_edge(definition);
        edges.insert(old, new);
    }
    let mut lumps = Vec::new();
    for lump in model.body(body.built.result).lumps.clone() {
        let mut shells = Vec::new();
        for shell in model.lump(lump).shells.clone() {
            let mut faces = Vec::new();
            for face in model.shell(shell).faces.clone() {
                let definition = model.face(face).clone();
                let mut loops = Vec::new();
                for lp in definition.loops {
                    let Loop::Coedges(coedges) = model.loop_(lp).clone() else {
                        panic!()
                    };
                    let coedges = coedges
                        .into_iter()
                        .map(|coedge| {
                            let coedge = model.coedge(coedge).clone();
                            model.add_coedge(Coedge {
                                edge: edges[&coedge.edge],
                                sense: reverse(coedge.sense),
                            })
                        })
                        .collect();
                    loops.push(model.add_loop(Loop::Coedges(coedges)));
                }
                let origins = model.face_origins(face).to_vec();
                faces.push(model.add_named_face(
                    Face {
                        loops,
                        ..definition
                    },
                    origins,
                ));
            }
            shells.push(model.add_shell(Shell { faces }));
        }
        lumps.push(model.add_lump(Lump { shells }));
    }
    body.built.result = model.add_body(Body {
        lumps,
        ..model.body(body.built.result).clone()
    });
    assert!(model.check(body.built.result).is_empty());
    body.names = crate::result_names(model, body.built.result);
}

#[test]
fn unsupported_carriers_and_invalid_or_exhausted_queries_are_explicit() {
    let body = build(&recipe(1, false));
    let selected = top_arc(&body);
    let (result, work) = with_operation_budget(0, || {
        body.edge_curve(selected.edge, &mut Decider::default())
    });
    assert_eq!(result.unwrap_err().code, "numbers.operation-budget");
    assert!(work.exhausted);
    assert_eq!(top_arc(&body).name, selected.name);
    let mut model = body.built.model.clone();
    let edge = model.edge(selected.edge).clone();
    let outsider = model.add_edge(edge.clone());
    assert_eq!(
        edge_curve(&model, body.built.result, outsider, &mut Decider::default())
            .unwrap_err()
            .code,
        "app.curve-edge"
    );

    let p = Point3::new(Real::zero(), Real::zero(), Real::zero());
    let x = Vector3::new(Real::one(), Real::zero(), Real::zero());
    let z = Vector3::new(Real::zero(), Real::zero(), Real::one());
    let (model, result, edge) = isolated_edge(
        Curve::Line(Line {
            origin: p.clone(),
            dir: x.clone(),
        }),
        None,
    );
    assert_eq!(
        edge_curve(&model, result, edge, &mut Decider::default())
            .unwrap_err()
            .code,
        "app.curve-trim"
    );
    let (model, result, edge) = isolated_edge(
        Curve::Circle(Circle {
            center: p,
            normal: z,
            ref_dir: x,
            radius_sq: Real::from_i64(25),
        }),
        Some(Point3::new(Real::from_i64(5), Real::zero(), Real::zero())),
    );
    let query = edge_curve(&model, result, edge, &mut Decider::default()).unwrap();
    assert!(matches!(query.span, CurveSpan::FullCircle));
    assert!(matches!(query.trim, CurveTrim::Finite { .. }));
}

#[test]
fn valid_conic_and_spline_bodies_keep_the_unsupported_boundary_explicit() {
    let mut input: Value = serde_json::from_str(include_str!(
        "../../../catallax-features/tests/fixtures/cone-distance/oblique-inward.json"
    ))
    .unwrap();
    input["steps"].as_array_mut().unwrap().pop();
    input["result"] = json!("clipped");
    let body = build(&Recipe::parse(&input.to_string()).unwrap());
    let mut ellipses = 0;
    for &edge in body.names.edges.keys() {
        let query = body.edge_curve(edge, &mut Decider::default()).unwrap();
        if matches!(query.form, CurveForm::Ellipse { .. }) {
            ellipses += 1;
            assert!(matches!(
                query.span,
                CurveSpan::FullEllipse | CurveSpan::EllipticArc { .. }
            ));
        }
    }
    assert!(ellipses > 0);

    let wrapper: Value = serde_json::from_str(include_str!(
        "../../../../.claude/skills/verify-catallax/references/examples/certified-recipe.json"
    ))
    .unwrap();
    let recipe = Recipe::parse(&wrapper["recipe"].to_string()).unwrap();
    let certified = crate::build_certified(&recipe, None, &mut Decider::default()).unwrap();
    let built = certified.built().built();
    assert!(built.model.check(built.result).is_empty());
    let edges = built.model.edges_of_body(built.result);
    assert_eq!(edges.len(), 2);
    for edge in edges {
        let query = edge_curve(&built.model, built.result, edge, &mut Decider::default()).unwrap();
        assert!(matches!(query.form, CurveForm::BSpline(_)));
        assert!(matches!(query.span, CurveSpan::Unsupported));
        assert!(
            query
                .definition
                .same_definition(built.model.geom().curve(query.curve))
        );
    }
}

fn isolated_edge(curve: Curve, anchor: Option<Point3>) -> (Model, BodyId, EdgeId) {
    isolated_span(curve, anchor.clone().zip(anchor), Sense::Same)
}

fn isolated_span(
    curve: Curve,
    ends: Option<(Point3, Point3)>,
    sense: Sense,
) -> (Model, BodyId, EdgeId) {
    let mut model = Model::new();
    let curve = model.geom_mut().add_curve(curve);
    attach_span(model, curve, ends, sense)
}

/// Adds a one-edge body over an existing curve of `model`.
fn attach_span(
    mut model: Model,
    curve: CurveId,
    ends: Option<(Point3, Point3)>,
    sense: Sense,
) -> (Model, BodyId, EdgeId) {
    use catallax_shapes::{Plane, Surface};
    use catallax_topology::{Body, BodyKind, Coedge, Edge, Face, Lump, Shell, Vertex};
    let mut vertex = |position| {
        let point = model.geom_mut().add_point(PointDef::Explicit(position));
        model.add_vertex(Vertex { point })
    };
    let (start, end) = match ends {
        Some((a, b)) => (Some(vertex(a)), Some(vertex(b))),
        None => (None, None),
    };
    let edge = model.add_edge(Edge {
        curve,
        sense,
        start,
        end,
    });
    let coedge = model.add_coedge(Coedge {
        edge,
        sense: Sense::Same,
    });
    let lp = model.add_loop(Loop::Coedges(vec![coedge]));
    let surface = model.geom_mut().add_surface(Surface::Plane(Plane {
        origin: Point3::new(Real::zero(), Real::zero(), Real::zero()),
        normal: Vector3::new(Real::zero(), Real::zero(), Real::one()),
        ref_dir: Vector3::new(Real::one(), Real::zero(), Real::zero()),
    }));
    let face = model.add_face(Face {
        surface,
        sense: Sense::Same,
        loops: vec![lp],
    });
    let shell = model.add_shell(Shell { faces: vec![face] });
    let lump = model.add_lump(Lump {
        shells: vec![shell],
    });
    let body = model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks: Default::default(),
    });
    (model, body, edge)
}

fn ellipse() -> Curve {
    Curve::Ellipse(Ellipse {
        center: Point3::new(Real::zero(), Real::zero(), Real::zero()),
        normal: Vector3::new(Real::zero(), Real::zero(), Real::one()),
        major_dir: Vector3::new(Real::one(), Real::zero(), Real::zero()),
        major_sq: Real::from_i64(25),
        minor_sq: Real::from_i64(9),
    })
}

fn frac(n: i64, d: i64) -> Real {
    Real::from(Rational::from_ratio(n, d).unwrap())
}

fn xy(x: Real, y: Real) -> Point3 {
    Point3::new(x, y, Real::zero())
}

/// Exact ellipse x^2/25 + y^2/9 = 1 points: the seam, vertices and (+-3, +-12/5).
fn ep(name: &str) -> Point3 {
    let i = Real::from_i64;
    match name {
        "seam" => xy(i(5), i(0)),
        "top" => xy(i(0), i(3)),
        "left" => xy(i(-5), i(0)),
        "bottom" => xy(i(0), i(-3)),
        "ur" => xy(i(3), frac(12, 5)),
        "ul" => xy(i(-3), frac(12, 5)),
        "ll" => xy(i(-3), frac(-12, 5)),
        _ => panic!("{name}"),
    }
}

/// Queries an isolated ellipse edge and checks the exact interior point lies on
/// the ellipse and strictly between the traversal endpoints on the right branch.
fn ellipse_arc(start: &str, end: &str, sense: Sense) -> (EdgeCurve, Point3) {
    let (model, result, edge) = isolated_span(ellipse(), Some((ep(start), ep(end))), sense);
    let query = edge_curve(&model, result, edge, &mut Decider::default()).unwrap();
    assert_eq!(query.sense, sense);
    assert!(query.definition.same_definition(&ellipse()));
    assert!(matches!(query.form, CurveForm::Ellipse { .. }));
    let CurveSpan::EllipticArc { interior, .. } = &query.span else {
        panic!("elliptic arc")
    };
    let on = &interior.x * &interior.x * frac(1, 25) + &interior.y * &interior.y * frac(1, 9)
        - Real::one();
    assert!(matches!(
        Decider::default().sign(&on),
        Outcome::ProvenCoincident(_)
    ));
    // Independent control: ordered counter-clockwise from the carrier's view.
    let (from, to) = match sense {
        Sense::Same => (ep(start), ep(end)),
        Sense::Reversed => (ep(end), ep(start)),
    };
    let ccw = |a: &Point3, b: &Point3, c: &Point3| {
        Decider::default().sign(&b.minus(a).cross(&c.minus(a)).z)
    };
    assert!(matches!(
        ccw(&from, interior, &to),
        Outcome::Proven(StrictSign::Positive)
    ));
    let interior = interior.clone();
    (query, interior)
}

fn crosses(query: &EdgeCurve) -> bool {
    let CurveSpan::EllipticArc { crosses_seam, .. } = query.span else {
        panic!()
    };
    crosses_seam
}

#[test]
fn exact_ellipse_arcs_keep_branch_sense_and_strict_seam_flag() {
    // Minor arc over the top, no seam.
    let (minor, interior) = ellipse_arc("ur", "ul", Sense::Same);
    assert!(!crosses(&minor));
    assert!(matches!(
        Decider::default().sign(&(interior.y - frac(12, 5))),
        Outcome::Proven(StrictSign::Positive)
    ));
    // Complementary major arc passes through the seam and the lower half.
    let (major, interior) = ellipse_arc("ul", "ur", Sense::Same);
    assert!(crosses(&major));
    assert!(matches!(
        Decider::default().sign(&(interior.y - frac(12, 5))),
        Outcome::Proven(StrictSign::Negative)
    ));
    // Seam endpoints are excluded from the strict crossing flag.
    assert!(!crosses(&ellipse_arc("seam", "top", Sense::Same).0));
    assert!(!crosses(&ellipse_arc("bottom", "seam", Sense::Same).0));
    // A seam crossing arc between two quadrants.
    assert!(crosses(&ellipse_arc("bottom", "top", Sense::Same).0));
    // Reversed storage walks the carrier from end to start.
    let (reversed_minor, interior) = ellipse_arc("ul", "ur", Sense::Reversed);
    assert!(!crosses(&reversed_minor));
    assert!(matches!(
        Decider::default().sign(&(interior.y - frac(12, 5))),
        Outcome::Proven(StrictSign::Positive)
    ));
    let CurveTrim::Finite { start, end } = &reversed_minor.trim else {
        panic!()
    };
    for (endpoint, name) in [(start, "ul"), (end, "ur")] {
        assert!(matches!(
            Decider::default().sign(&(&endpoint.position.x - &ep(name).x)),
            Outcome::ProvenCoincident(_)
        ));
    }
    assert!(crosses(&ellipse_arc("ur", "ul", Sense::Reversed).0));
}

#[test]
fn full_ellipse_is_explicit_for_ring_and_coincident_endpoints() {
    let (model, result, edge) = isolated_span(ellipse(), None, Sense::Same);
    let query = edge_curve(&model, result, edge, &mut Decider::default()).unwrap();
    assert!(matches!(query.span, CurveSpan::FullEllipse));
    assert!(matches!(query.trim, CurveTrim::Closed));
    let (model, result, edge) = isolated_span(ellipse(), Some((ep("ul"), ep("ul"))), Sense::Same);
    let query = edge_curve(&model, result, edge, &mut Decider::default()).unwrap();
    assert!(matches!(query.span, CurveSpan::FullEllipse));
    assert!(matches!(query.trim, CurveTrim::Finite { .. }));
}

fn cone_recipe(base_radius: &str) -> Recipe {
    let mut input: Value = serde_json::from_str(include_str!(
        "../../../catallax-features/tests/fixtures/cone-distance/oblique-inward.json"
    ))
    .unwrap();
    input["steps"].as_array_mut().unwrap().pop();
    input["result"] = json!("clipped");
    input["steps"][0]["base_radius"] = json!(base_radius);
    Recipe::parse(&input.to_string()).unwrap()
}

#[test]
fn built_ellipse_edges_keep_origins_names_and_follow_an_upstream_edit() {
    let body = build(&cone_recipe("4"));
    let mut found = Vec::new();
    for &edge in body.names.edges.keys() {
        let query = body.edge_curve(edge, &mut Decider::default()).unwrap();
        if let CurveForm::Ellipse { .. } = query.form {
            if let CurveSpan::EllipticArc { interior, .. } = &query.span {
                let CurveForm::Ellipse {
                    center,
                    normal,
                    major_sq,
                    ..
                } = &query.form
                else {
                    unreachable!()
                };
                assert!(matches!(
                    Decider::default().sign(&interior.minus(center).dot(normal)),
                    Outcome::ProvenCoincident(_)
                ));
                assert!(matches!(
                    Decider::default()
                        .sign(&(interior.minus(center).norm_squared() - major_sq.clone())),
                    Outcome::Proven(StrictSign::Negative) | Outcome::ProvenCoincident(_)
                ));
            } else {
                assert!(matches!(query.span, CurveSpan::FullEllipse));
            }
            assert!(!query.origins.is_empty());
            assert!(!query.constructions.is_empty());
            found.push(query);
        }
    }
    assert!(!found.is_empty());
    // The same stable name resolves again after the cone's base radius changes.
    let name = found[0].name.clone().unwrap();
    let edited = build(&cone_recipe("5"));
    let ids = edited
        .built
        .model
        .edges_with_name(edited.built.result, &name);
    assert_eq!(ids.len(), 1);
    let again = edited.edge_curve(ids[0], &mut Decider::default()).unwrap();
    assert!(matches!(
        again.span,
        CurveSpan::EllipticArc { .. } | CurveSpan::FullEllipse
    ));
    assert_eq!(again.origins, found[0].origins);
    assert!(!again.definition.same_definition(&found[0].definition));
}

fn open_conic_recipe(name: &str, edit: impl FnOnce(&mut Value)) -> Recipe {
    let text = match name {
        "parabola" => include_str!("../../tests/fixtures/open-conics/parabola.json"),
        "hyperbola" => include_str!("../../tests/fixtures/open-conics/hyperbola.json"),
        _ => panic!("{name}"),
    };
    let mut value: Value = serde_json::from_str(text).unwrap();
    edit(&mut value);
    Recipe::parse(&value.to_string()).unwrap()
}

fn is_zero(value: Real) -> bool {
    matches!(
        Decider::default().sign(&value),
        Outcome::ProvenCoincident(_)
    )
}

fn strictly(value: Real, sign: StrictSign) -> bool {
    matches!(Decider::default().sign(&value), Outcome::Proven(s) if s == sign)
}

/// The carrier point at transverse coordinate `y`, from the conic equation alone.
fn carrier_point(form: &CurveForm, y: &Real) -> Point3 {
    let unit = |v: &Vector3| v.scaled(&(Real::one() / v.norm_squared().sqrt()));
    let (origin, normal, axis, x) = match form {
        CurveForm::Parabola {
            vertex,
            normal,
            axis_dir,
            focal,
        } => (
            vertex,
            normal,
            axis_dir,
            y * y / (Real::from_i64(4) * focal.clone()),
        ),
        CurveForm::Hyperbola {
            center,
            normal,
            axis_dir,
            transverse_sq,
            conjugate_sq,
            branch,
        } => {
            let x = (transverse_sq.clone() * (Real::one() + y * y / conjugate_sq.clone())).sqrt();
            let x = if *branch == StrictSign::Positive {
                x
            } else {
                -x
            };
            (center, normal, axis_dir, x)
        }
        _ => panic!("open conic"),
    };
    let across = normal.cross(axis);
    origin
        .plus(&unit(axis).scaled(&x))
        .plus(&unit(&across).scaled(y))
}

fn same_position(a: &Point3, b: &Point3) -> bool {
    is_zero(a.minus(b).norm_squared())
}

/// Queries every conic edge and checks it against the fixture's independent
/// equations: the stock cone `x^2 + y^2 = (4 - 3z/4)^2`, the cutting surface
/// (`plane`) and the finite height 0..4 of the stock.
fn check_open_arcs(
    body: &ApplicationBody,
    parabola: bool,
    plane: &dyn Fn(&Point3) -> Real,
) -> Vec<EdgeCurve> {
    let mut found = Vec::new();
    for &edge in body.names.edges.keys() {
        let query = body.edge_curve(edge, &mut Decider::default()).unwrap();
        let expected = match (&query.form, parabola) {
            (CurveForm::Parabola { .. }, true) | (CurveForm::Hyperbola { .. }, false) => true,
            (CurveForm::Parabola { .. } | CurveForm::Hyperbola { .. }, _) => {
                panic!("unexpected conic kind")
            }
            _ => false,
        };
        if !expected {
            continue;
        }
        let CurveSpan::OpenConicArc {
            interior,
            first,
            last,
        } = &query.span
        else {
            panic!("open conic arc")
        };
        assert!(
            query
                .definition
                .same_definition(body.built.model.geom().curve(query.curve))
        );
        let CurveTrim::Finite { start, end } = &query.trim else {
            panic!("finite conic edge")
        };
        let (from, to) = match query.sense {
            Sense::Same => (start, end),
            Sense::Reversed => (end, start),
        };
        // The interval reproduces the carrier-order endpoints and the interior.
        assert!(same_position(
            &carrier_point(&query.form, first),
            &from.position
        ));
        assert!(same_position(
            &carrier_point(&query.form, last),
            &to.position
        ));
        let origin = match &query.form {
            CurveForm::Parabola { vertex, .. } => vertex,
            CurveForm::Hyperbola { center, .. } => center,
            _ => unreachable!(),
        };
        let (CurveForm::Parabola {
            normal, axis_dir, ..
        }
        | CurveForm::Hyperbola {
            normal, axis_dir, ..
        }) = &query.form
        else {
            unreachable!()
        };
        let across = normal.cross(axis_dir);
        let across = across.scaled(&(Real::one() / across.norm_squared().sqrt()));
        let coordinate = interior.minus(origin).dot(&across);
        assert!(same_position(
            &carrier_point(&query.form, &coordinate),
            interior
        ));
        // Monotone interval, and the interior is strictly inside it.
        assert!(!is_zero(last - first));
        assert!(strictly(
            (&coordinate - first) * (&coordinate - last),
            StrictSign::Negative
        ));
        for point in [&from.position, &to.position, interior] {
            let r = Real::from_i64(4) - point.z.clone() * frac(3, 4);
            assert!(is_zero(&point.x * &point.x + &point.y * &point.y - &r * &r));
            assert!(is_zero(plane(point)));
        }
        for point in [&from.position, &to.position] {
            assert!(!strictly(point.z.clone(), StrictSign::Negative));
            assert!(!strictly(
                point.z.clone() - Real::from_i64(4),
                StrictSign::Positive
            ));
        }
        assert!(strictly(interior.z.clone(), StrictSign::Positive));
        assert!(strictly(
            interior.z.clone() - Real::from_i64(4),
            StrictSign::Negative
        ));
        // Finite-domain exclusion: carrier points just outside the interval
        // satisfy the conic equation but fall outside the stock's height.
        let (low, high) = if strictly(last - first, StrictSign::Positive) {
            (first, last)
        } else {
            (last, first)
        };
        let width = (high - low) * frac(1, 10);
        for outside in [low - &width, high + &width] {
            let z = carrier_point(&query.form, &outside).z;
            assert!(
                strictly(z.clone(), StrictSign::Negative)
                    || strictly(z - Real::from_i64(4), StrictSign::Positive)
            );
        }
        assert!(!query.origins.is_empty());
        found.push(query);
    }
    found
}

fn parabola_plane(d: i64) -> impl Fn(&Point3) -> Real {
    move |p| Real::from_i64(4) * p.x.clone() + Real::from_i64(3) * p.z.clone() - Real::from_i64(d)
}

fn hyperbola_plane(half: i64) -> impl Fn(&Point3) -> Real {
    hyperbola_plane_sq(Real::from_i64(half * half))
}

fn hyperbola_plane_sq(square: Real) -> impl Fn(&Point3) -> Real {
    move |p| &p.x * &p.x - square.clone()
}

#[test]
fn parabola_arcs_keep_exact_interval_branch_and_traversal() {
    let body = build(&open_conic_recipe("parabola", |_| {}));
    let arcs = check_open_arcs(&body, true, &parabola_plane(12));
    assert_eq!(arcs.len(), 2);
    assert!(arcs.iter().all(|q| q.sense == Sense::Same));
    assert!(arcs.iter().all(|q| !q.constructions.is_empty()));
    let mut reversed = build(&open_conic_recipe("parabola", |_| {}));
    reverse_edges(&mut reversed);
    let flipped = check_open_arcs(&reversed, true, &parabola_plane(12));
    assert_eq!(flipped.len(), 2);
    assert!(flipped.iter().all(|q| q.sense == Sense::Reversed));
    // Reversed storage walks the same carrier interval, trim swapped.
    for flip in &flipped {
        let twin = arcs
            .iter()
            .find(|q| {
                q.origins == flip.origins && q.name == flip.name && {
                    let (
                        CurveSpan::OpenConicArc { first: a, .. },
                        CurveSpan::OpenConicArc { first: b, .. },
                    ) = (&q.span, &flip.span)
                    else {
                        unreachable!()
                    };
                    is_zero(a - b)
                }
            })
            .expect("same arc after reversal");
        let (CurveTrim::Finite { start: s0, end: e0 }, CurveTrim::Finite { start: s1, end: e1 }) =
            (&twin.trim, &flip.trim)
        else {
            unreachable!()
        };
        assert!(same_position(&s0.position, &e1.position));
        assert!(same_position(&e0.position, &s1.position));
    }
}

#[test]
fn hyperbola_arcs_on_both_sides_keep_exact_interval_and_origins() {
    let body = build(&open_conic_recipe("hyperbola", |_| {}));
    let arcs = check_open_arcs(&body, false, &hyperbola_plane(2));
    assert_eq!(arcs.len(), 2);
    // The two cutting planes give the two mirrored carriers.
    let sides: Vec<_> = arcs
        .iter()
        .map(|q| {
            let CurveForm::Hyperbola { center, .. } = &q.form else {
                unreachable!()
            };
            Decider::default().sign(&center.x)
        })
        .collect();
    assert!(sides.contains(&Outcome::Proven(StrictSign::Positive)));
    assert!(sides.contains(&Outcome::Proven(StrictSign::Negative)));
    assert_ne!(arcs[0].origins, arcs[1].origins);
    let mut reversed = build(&open_conic_recipe("hyperbola", |_| {}));
    reverse_edges(&mut reversed);
    let flipped = check_open_arcs(&reversed, false, &hyperbola_plane(2));
    assert!(flipped.iter().all(|q| q.sense == Sense::Reversed));
    assert_eq!(flipped.len(), 2);
}

#[test]
fn edited_trims_resolve_by_origin_and_change_the_interval() {
    for (name, parabola) in [("parabola", true), ("hyperbola", false)] {
        let old = build(&open_conic_recipe(name, |_| {}));
        let before = check_open_arcs(
            &old,
            parabola,
            &if parabola {
                Box::new(parabola_plane(12)) as Box<dyn Fn(&Point3) -> Real>
            } else {
                Box::new(hyperbola_plane(2))
            },
        );
        let edited = build(&open_conic_recipe(name, |value| {
            if parabola {
                // The plane 4x + 3z = 10 is still parallel to a generator.
                value["steps"][1]["pose"]["translation"] = json!(["5/2", "-10", "0"]);
            } else {
                value["steps"][1]["size"] = json!(["5", "20", "20"]);
                value["steps"][1]["pose"]["translation"] = json!(["-5/2", "-10", "0"]);
            }
        }));
        let after = check_open_arcs(
            &edited,
            parabola,
            &if parabola {
                Box::new(parabola_plane(10)) as Box<dyn Fn(&Point3) -> Real>
            } else {
                Box::new(hyperbola_plane_sq(frac(25, 4)))
            },
        );
        assert_eq!(after.len(), before.len());
        for old_arc in &before {
            let name = old_arc.name.clone().unwrap();
            let again: Vec<_> = after
                .iter()
                .filter(|q| q.name.as_ref() == Some(&name))
                .collect();
            let same_name = before
                .iter()
                .filter(|q| q.name.as_ref() == Some(&name))
                .count();
            assert_eq!(again.len(), same_name);
            assert!(again.iter().any(|q| q.origins == old_arc.origins));
            // The edit moves the cutting surface, so no matching arc keeps the
            // old exact interval.
            let CurveSpan::OpenConicArc {
                first: old_first,
                last: old_last,
                ..
            } = &old_arc.span
            else {
                unreachable!()
            };
            assert!(again.iter().all(|q| {
                let CurveSpan::OpenConicArc { first, last, .. } = &q.span else {
                    unreachable!()
                };
                !(is_zero(first - old_first) && is_zero(last - old_last))
            }));
        }
    }
}

#[test]
fn open_conic_refusals_stay_typed() {
    let body = build(&open_conic_recipe("parabola", |_| {}));
    let arc = check_open_arcs(&body, true, &parabola_plane(12)).remove(0);
    let (CurveSpan::OpenConicArc { first, last, .. }, CurveTrim::Finite { .. }) =
        (&arc.span, &arc.trim)
    else {
        unreachable!()
    };
    let point = |y: &Real| carrier_point(&arc.form, y);
    let (a, b) = (point(first), point(last));
    let model = body.built.model.clone();
    let code = |ends: Option<(Point3, Point3)>| {
        let (model, result, edge) = attach_span(model.clone(), arc.curve, ends, Sense::Same);
        edge_curve(&model, result, edge, &mut Decider::default())
    };
    // The valid pair reproduces the original interval.
    let query = code(Some((a.clone(), b.clone()))).unwrap();
    let CurveSpan::OpenConicArc {
        first: f2,
        last: l2,
        ..
    } = &query.span
    else {
        panic!()
    };
    assert!(is_zero(f2 - first) && is_zero(l2 - last));
    // Reversed carrier order, coincident ends and a ring all refuse.
    assert_eq!(
        code(Some((b.clone(), a.clone()))).unwrap_err().code,
        "app.curve-trim"
    );
    assert_eq!(
        code(Some((a.clone(), a.clone()))).unwrap_err().code,
        "app.curve-trim"
    );
    assert_eq!(code(None).unwrap_err().code, "app.curve-trim");
    // An exhausted allowance refuses instead of weakening the answer.
    let (result, work) =
        with_operation_budget(0, || body.edge_curve(arc.edge, &mut Decider::default()));
    assert_eq!(result.unwrap_err().code, "numbers.operation-budget");
    assert!(work.exhausted);
}
