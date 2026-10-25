//! Public convex-profile caps: independent set oracles, names and consumers.

use catallax_features::end_condition::curved::{AuthoredCap, Branch, Location};
use catallax_field::{CurvedField, Plan, Workspace};
use catallax_numbers::{Ball, Decider, Outcome, Rational, Real};
use catallax_recipe::{Operation, Origin, Recipe, Role, StepId};
use catallax_shapes::RPoint3;
use serde_json::{Value, json};

const CASES: [&str; 9] = [
    include_str!("fixtures/polygon-limits/sphere-triangle.json"),
    include_str!("fixtures/polygon-limits/sphere-pentagon.json"),
    include_str!("fixtures/polygon-limits/sphere-rotated.json"),
    include_str!("fixtures/polygon-limits/cylinder-triangle.json"),
    include_str!("fixtures/polygon-limits/cylinder-pentagon.json"),
    include_str!("fixtures/polygon-limits/cylinder-rotated.json"),
    include_str!("fixtures/polygon-limits/cone-triangle.json"),
    include_str!("fixtures/polygon-limits/cone-pentagon.json"),
    include_str!("fixtures/polygon-limits/cone-rotated.json"),
];

fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn id(s: &str) -> StepId {
    StepId::new(s).unwrap()
}
fn recipe(v: &Value) -> Recipe {
    Recipe::parse(&v.to_string()).unwrap()
}
fn compile(recipe: &Recipe, dec: &mut Decider) -> AuthoredCap {
    let (direction, limit) = match &recipe.steps[2].op {
        Operation::ExtrudeNext { direction, .. } => (direction.as_ref(), None),
        Operation::ExtrudeTo {
            direction, limit, ..
        } => (direction.as_ref(), Some(limit)),
        _ => unreachable!(),
    };
    AuthoredCap::compile(recipe, &id("target"), &id("profile"), direction, limit, dec).unwrap()
}
fn point(cap: &AuthoredCap, uv: [Rational; 2], t: Rational) -> RPoint3 {
    let frame = cap.definition().frame();
    frame
        .point(&uv.map(Real::from))
        .plus(&frame.direction.scaled(&t.into()))
        .as_rational()
        .unwrap()
}
fn variant(text: &str, to: bool, exit: bool) -> Value {
    let mut v: Value = serde_json::from_str(text).unwrap();
    let sphere = v["steps"][1]["op"] == "sphere";
    let axis = if sphere { 2 } else { 0 };
    if exit {
        v["steps"][0]["plane"]["translation"][axis] = json!("5");
    }
    v["steps"][2]["op"] = json!(if to { "extrude_to" } else { "extrude_next" });
    if to {
        v["steps"][2]["limit"] =
            json!({"face":{"step":"target","role":{"side":{"segment":"wall"}}}});
        let mut direction = ["0"; 3];
        direction[axis] = "2";
        v["steps"][2]["direction"] = json!(direction);
    }
    v
}
fn vertices(v: &Value) -> Vec<[Rational; 2]> {
    v["steps"][0]["loops"][0]
        .as_array()
        .unwrap()
        .iter()
        .map(|s| std::array::from_fn(|i| q(s["to"][i].as_str().unwrap())))
        .collect()
}
fn equal(a: &Real, b: Rational, dec: &mut Decider) {
    assert!(matches!(
        dec.sign(&(a - &Real::from(b))),
        Outcome::ProvenCoincident(_)
    ));
}

#[test]
fn triangles_pentagons_and_rotated_polygons_keep_entry_exit_sides_in_both_modes_and_kernels() {
    for (case, text) in CASES.iter().enumerate() {
        for to in [false, true] {
            for exit in [false, true] {
                let v = variant(text, to, exit);
                let recipe = recipe(&v);
                let mut dec = Decider::default();
                let cap = compile(&recipe, &mut dec);
                assert_eq!(
                    cap.definition().branch(),
                    if exit { Branch::Exit } else { Branch::Entry }
                );
                let root = match (case / 3, exit) {
                    (0 | 1, false) => 3,
                    (0 | 1, true) => 2,
                    (2, false) => 2,
                    (2, true) => 3,
                    _ => unreachable!(),
                };
                let scale = if to { 2 } else { 1 };
                equal(
                    &cap.definition()
                        .height(&[Real::zero(), Real::zero()], &mut dec)
                        .unwrap(),
                    Rational::from_ratio(root, scale).unwrap(),
                    &mut dec,
                );
                let radial = Real::from(if case / 3 == 2 {
                    q("143/16")
                } else {
                    q("63/16")
                })
                .sqrt();
                let expected = if exit {
                    radial
                } else {
                    Real::from_i64(5) - radial
                } / Real::from_i64(scale);
                assert!(matches!(
                    dec.sign(
                        &(cap
                            .definition()
                            .height(&[q("1/4").into(), Real::zero()], &mut Decider::default())
                            .unwrap()
                            - expected)
                    ),
                    Outcome::ProvenCoincident(_)
                ));
                let mut built = catallax_features::evaluate(&recipe, &mut dec).unwrap();
                assert!(
                    built.model.check(built.result).is_empty(),
                    "{}",
                    recipe.name
                );
                assert!(built.model.body(built.result).marks.is_empty());
                let plan = Plan::compile(&recipe).unwrap();
                let query = CurvedField::new(cap.clone(), id("cap"));
                let mut workspace = Workspace::default();
                // All three domains contain zero, but exclude this point of
                // their bounding rectangle. This detects a rectangle takeover.
                for (uv, t, expected) in [
                    ([q("0"), q("0")], q("1/2"), catallax_field::Location::Inside),
                    (
                        [q("2/5"), q("2/5")],
                        q("1/2"),
                        catallax_field::Location::Outside,
                    ),
                    ([q("0"), q("0")], q("4"), catallax_field::Location::Outside),
                    ([q("0"), q("0")], q("-1"), catallax_field::Location::Outside),
                ] {
                    let p = point(&cap, uv.clone(), t.clone());
                    assert_eq!(
                        plan.proven(&p, &mut workspace).location,
                        expected,
                        "{} at {p:?}",
                        recipe.name
                    );
                    let definition = cap
                        .definition()
                        .location(&uv.map(Real::from), &t.into(), &mut dec)
                        .unwrap();
                    let want = if expected == catallax_field::Location::Inside {
                        Location::Inside
                    } else {
                        Location::Outside
                    };
                    assert_eq!(definition, want);
                    assert_eq!(query.location(&p, &mut dec).unwrap(), want);
                    let want = if expected == catallax_field::Location::Inside {
                        catallax_boolean::Containment::Inside
                    } else {
                        catallax_boolean::Containment::Outside
                    };
                    assert!(
                        matches!(catallax_boolean::classify_point(&mut built.model, built.result, &p, &mut dec), Outcome::Proven(found) if found == want)
                    );
                }
                let vertices = vertices(&v);
                let names: Vec<_> = built
                    .model
                    .faces_of_body(built.result)
                    .into_iter()
                    .flat_map(|f| built.model.face_origins(f))
                    .cloned()
                    .collect();
                assert!(names.iter().all(|o| o.step == id("cap")));
                for role in [Role::Start, Role::End] {
                    assert!(names.iter().any(|o| o.role == role));
                }
                for (i, b) in vertices.iter().enumerate() {
                    let a = &vertices[(i + vertices.len() - 1) % vertices.len()];
                    let uv = std::array::from_fn(|j| (&a[j] + &b[j]) / q("2"));
                    let p = point(&cap, uv, q("1/2"));
                    let origin = Origin {
                        step: id("cap"),
                        role: Role::Side {
                            segment: cap.side_ids()[i].clone(),
                            path: None,
                        },
                    };
                    assert_eq!(
                        cap.side_ids()[i].as_str(),
                        v["steps"][0]["loops"][0][i]["id"].as_str().unwrap()
                    );
                    assert!(names.contains(&origin), "missing side {origin:?}");
                    assert_eq!(query.origins(&p, &mut dec).unwrap(), vec![origin.clone()]);
                    assert_eq!(
                        plan.proven(&p, &mut workspace).location,
                        catallax_field::Location::Boundary
                    );
                    assert!(matches!(
                        catallax_boolean::classify_point(
                            &mut built.model,
                            built.result,
                            &p,
                            &mut dec
                        ),
                        Outcome::ProvenCoincident(_)
                    ));
                    let handle = plan
                        .faces()
                        .iter()
                        .position(|f| f.origin == origin)
                        .and_then(|i| plan.result_face(i))
                        .unwrap();
                    assert_eq!(
                        plan.face_query(handle, &p, &mut workspace).unwrap().state,
                        catallax_field::FaceState::OnFace
                    );
                }
            }
        }
    }
}

#[test]
fn independent_rational_grid_oracle_checks_the_polygon_and_first_crossing() {
    for (case, text) in CASES.iter().enumerate() {
        let v = variant(text, false, false);
        let recipe = recipe(&v);
        let cap = compile(&recipe, &mut Decider::default());
        let plan = Plan::compile(&recipe).unwrap();
        let mut workspace = Workspace::default();
        for i in -3..=3 {
            for j in -3..=3 {
                let (u, v) = (
                    Rational::from_ratio(i, 5).unwrap(),
                    Rational::from_ratio(j, 5).unwrap(),
                );
                // These inequalities are authored independently from the cap's
                // side construction and use no support equations or B-rep.
                let margins = match case % 3 {
                    0 => vec![&v + &q("1/2"), q("1/2") - (&v + &(&u.abs() * &q("2")))],
                    1 => vec![
                        &v + &q("1/2"),
                        q("1/2") - (&u.abs() - &(&v / &q("2"))),
                        q("1/2") - (&u.abs() + &v),
                    ],
                    2 => vec![q("1/2") - (&u.abs() + &v.abs())],
                    _ => unreachable!(),
                };
                let in_profile = margins.iter().all(|m| *m >= q("0"));
                let on_profile = margins.iter().any(|m| *m == q("0"));
                for t in [q("1"), q("5/2"), q("7/2"), q("8")] {
                    let radial_sq = match case / 3 {
                        0 => &u * &u + &v * &v,
                        1 | 2 => &u * &u,
                        _ => unreachable!(),
                    };
                    let radius = if case / 3 == 2 {
                        q("3") - &v / &q("5")
                    } else {
                        q("2")
                    };
                    let delta = &t - &q("5");
                    // Entry keeps only the exterior before the near wall;
                    // exterior after the far wall must also be excluded.
                    let inside =
                        in_profile && t < q("5") && radial_sq + &delta * &delta > &radius * &radius;
                    let p = point(&cap, [u.clone(), v.clone()], t);
                    let actual = plan.proven(&p, &mut workspace).location;
                    assert_eq!(
                        actual,
                        if inside && on_profile {
                            catallax_field::Location::Boundary
                        } else if inside {
                            catallax_field::Location::Inside
                        } else {
                            catallax_field::Location::Outside
                        },
                        "case {case}, {p:?}"
                    );
                }
            }
        }
    }
}

#[test]
fn enclosure_failure_remains_coverage_even_when_every_polygon_ray_has_a_hit() {
    let mut v = variant(CASES[2], false, false);
    v["steps"][0]["loops"][0] = json!([
        {"kind":"line","id":"a","to":["0","-3/2"]},
        {"kind":"line","id":"b","to":["3/2","0"]},
        {"kind":"line","id":"c","to":["0","3/2"]},
        {"kind":"line","id":"d","to":["-3/2","0"]}
    ]);
    // The convex diamond lies inside radius 3/2 < 2, so it has a positive
    // sphere entry everywhere. Its enclosing square crosses the silhouette.
    assert!(q("3/2") * q("3/2") < q("4"));
    let recipe = recipe(&v);
    let refusal = AuthoredCap::compile(
        &recipe,
        &id("target"),
        &id("profile"),
        None,
        None,
        &mut Decider::default(),
    )
    .unwrap_err();
    assert_eq!(refusal.code, "features.unsupported-curved-limit");
    assert!(refusal.detail.contains("discriminant"));
    assert_eq!(
        Plan::compile(&recipe).unwrap_err().code,
        "features.unsupported-curved-limit"
    );
}

#[test]
fn clockwise_transformed_and_reversed_scaled_rays_keep_exact_physical_sides() {
    for (start, w, expected_branch, expected_root) in [
        ("0", "3", Branch::Entry, "1"),
        ("5", "3", Branch::Exit, "2/3"),
        ("10", "-3", Branch::Entry, "1"),
        ("5", "-3", Branch::Exit, "2/3"),
    ] {
        let mut v = variant(CASES[1], false, false);
        // Reverse the traversal while keeping each id on its physical edge.
        let old = v["steps"][0]["loops"][0].as_array().unwrap().clone();
        let n = old.len();
        v["steps"][0]["loops"][0] = Value::Array(
            (0..n)
                .rev()
                .map(|i| json!({"kind":"line","id":old[(i+1)%n]["id"],"to":old[i]["to"]}))
                .collect(),
        );
        let origin_y = q("-11") - q(start);
        v["steps"][0]["plane"] =
            json!({"rotation":["1","1","0","0"],"translation":["7",origin_y.to_string(),"3"]});
        v["steps"][1]["pose"] = json!({"rotation":["1","1","0","0"],"translation":["7","-16","3"]});
        v["steps"][2]["direction"] = json!(["0", (-q(w)).to_string(), "0"]);
        let recipe = recipe(&v);
        let mut dec = Decider::default();
        let cap = compile(&recipe, &mut dec);
        assert_eq!(cap.definition().branch(), expected_branch);
        equal(
            &cap.definition()
                .height(&[Real::zero(), Real::zero()], &mut dec)
                .unwrap(),
            q(expected_root),
            &mut dec,
        );
        let plan = Plan::compile(&recipe).unwrap();
        let query = CurvedField::new(cap.clone(), id("cap"));
        let mut built = catallax_features::evaluate(&recipe, &mut dec).unwrap();
        let vertices = vertices(&v);
        for (i, b) in vertices.iter().enumerate() {
            let a = &vertices[(i + vertices.len() - 1) % vertices.len()];
            let uv = std::array::from_fn(|j| (&a[j] + &b[j]) / q("2"));
            let p = point(&cap, uv, q("1/4"));
            assert!(matches!(
                catallax_boolean::classify_point(&mut built.model, built.result, &p, &mut dec),
                Outcome::ProvenCoincident(_)
            ));
            assert_eq!(
                plan.proven(&p, &mut Workspace::default()).location,
                catallax_field::Location::Boundary
            );
            assert_eq!(
                query.origins(&p, &mut dec).unwrap(),
                vec![Origin {
                    step: id("cap"),
                    role: Role::Side {
                        segment: cap.side_ids()[i].clone(),
                        path: None
                    }
                }]
            );
        }
    }
}

#[test]
fn polygon_brep_display_and_step_are_complete_and_extrusion_owned() {
    for text in CASES {
        for to in [false, true] {
            for exit in [false, true] {
                let built = catallax_features::evaluate(
                    &recipe(&variant(text, to, exit)),
                    &mut Decider::default(),
                )
                .unwrap();
                let mesh = catallax_mesh::mesh_body(
                    &built.model,
                    built.result,
                    &Default::default(),
                    &mut Decider::default(),
                )
                .unwrap();
                assert!(mesh.marks.is_empty());
                assert!(!mesh.value.triangles.is_empty());
                let mut bytes = Vec::new();
                let report = catallax_step::write_step(
                    &built.model,
                    &[built.result],
                    &Default::default(),
                    &mut Decider::default(),
                    &mut bytes,
                )
                .unwrap();
                assert!(report.marks.is_empty());
                assert_eq!(report.value.solids, 1);
                assert!(bytes.starts_with(b"ISO-10303-21;"));
            }
        }
    }
}

#[test]
fn requested_cylinder_triangle_volume_encloses_an_independent_exact_integral() {
    let built =
        catallax_features::evaluate(&Recipe::parse(CASES[3]).unwrap(), &mut Decider::default())
            .unwrap();
    let measured = catallax_measure::volume(
        &built.model,
        built.result,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert!(measured.marks.is_empty());
    let b = |s| Ball::from_rational(&q(s), 192);
    // For this triangle the vertical width is 1-2|u| on [-1/2,1/2].
    // Integral sqrt(4-u²)*(1-2|u|) du = 11 sqrt(15)/4
    // + 4 asin(1/4) - 32/3. Its area is 1/2; entry height is
    // 5-sqrt(4-u²), independently of the construction cell decomposition.
    let root = b("15").sqrt().unwrap();
    let expected = b("5/2") - (b("11/4") * root + b("4") * b("1/4").asin().unwrap() - b("32/3"));
    assert!(measured.value.lo <= expected.lo() && measured.value.hi >= expected.hi());
    assert!(measured.value.width() <= q("1e-11"));
}

#[test]
fn every_stock_and_polygon_volume_keeps_the_exact_entry_exit_integral_identity() {
    for text in CASES {
        let v = variant(text, false, false);
        let points = vertices(&v);
        let mut twice_area = q("0");
        for (i, a) in points.iter().enumerate() {
            let b = &points[(i + 1) % points.len()];
            twice_area = twice_area + (&a[0] * &b[1] - &a[1] * &b[0]);
        }
        let expected = twice_area.abs() * q("5/2");
        let volumes: [catallax_measure::Enclosure; 2] = std::array::from_fn(|i| {
            let built = catallax_features::evaluate(
                &recipe(&variant(text, false, i == 1)),
                &mut Decider::default(),
            )
            .unwrap();
            catallax_measure::volume(
                &built.model,
                built.result,
                &Default::default(),
                &mut Decider::default(),
            )
            .unwrap()
            .value
        });
        // Every near-wall height 5-r(q) plus its corresponding exit r(q)
        // is exactly 5. This identity integrates the authored polygon itself.
        assert!(&volumes[0].lo + &volumes[1].lo <= expected);
        assert!(&volumes[0].hi + &volumes[1].hi >= expected);
        assert!(volumes.iter().all(|v| v.width() <= q("1e-11")));
    }
}

#[test]
fn collinear_authored_sides_keep_their_names_after_coplanar_face_merge() {
    let mut v = variant(CASES[0], false, false);
    v["steps"][0]["loops"][0] = json!([
        {"kind":"line","id":"bottom-a","to":["0","-1/2"]},
        {"kind":"line","id":"bottom-b","to":["1/2","-1/2"]},
        {"kind":"line","id":"right","to":["0","1/2"]},
        {"kind":"line","id":"left","to":["-1/2","-1/2"]}
    ]);
    let recipe = recipe(&v);
    let built = catallax_features::evaluate(&recipe, &mut Decider::default()).unwrap();
    let plan = Plan::compile(&recipe).unwrap();
    for segment in ["bottom-a", "bottom-b", "right", "left"] {
        let origin = Origin {
            step: id("cap"),
            role: Role::Side {
                segment: catallax_recipe::ElementId::new(segment).unwrap(),
                path: None,
            },
        };
        assert!(
            !built
                .model
                .faces_with_origin(built.result, &origin)
                .is_empty()
        );
        assert!(plan.faces().iter().any(|f| f.origin == origin));
    }
}

#[test]
fn radius_edits_resolve_every_original_polygon_face_and_physical_edge_name() {
    for text in CASES {
        let v: Value = serde_json::from_str(text).unwrap();
        let mut edited = v.clone();
        if edited["steps"][1]["op"] == "cone" {
            edited["steps"][1]["base_radius"] = json!("9/2");
            edited["steps"][1]["top_radius"] = json!("5/2");
        } else {
            edited["steps"][1]["radius"] = json!("5/2");
        }
        let before = catallax_features::evaluate(&recipe(&v), &mut Decider::default()).unwrap();
        let after = catallax_features::evaluate(&recipe(&edited), &mut Decider::default()).unwrap();
        for face in before.model.faces_of_body(before.result) {
            for origin in before.model.face_origins(face) {
                assert!(
                    !after
                        .model
                        .faces_with_origin(after.result, origin)
                        .is_empty(),
                    "{origin:?}"
                );
            }
        }
        for edge in before.model.edges_of_body(before.result) {
            let name = before
                .model
                .edge_name(before.result, edge)
                .expect("named physical edge");
            assert!(
                !after.model.edges_with_name(after.result, &name).is_empty(),
                "{name:?}"
            );
        }
        let plan = Plan::compile(&recipe(&edited)).unwrap();
        for origin in before
            .model
            .faces_of_body(before.result)
            .iter()
            .flat_map(|f| before.model.face_origins(*f))
        {
            assert!(plan.faces().iter().any(|f| &f.origin == origin));
        }
    }
}

#[test]
fn nonconvex_degenerate_trim_crossing_and_exhausted_requests_stay_fail_closed() {
    for points in [
        json!([
            ["-1/2", "-1/2"],
            ["1/2", "-1/2"],
            ["0", "0"],
            ["1/2", "1/2"],
            ["-1/2", "1/2"]
        ]),
        json!([
            ["-1/2", "-1/2"],
            ["1/2", "1/2"],
            ["1/2", "-1/2"],
            ["-1/2", "1/2"]
        ]),
        json!([["0", "0"], ["1/2", "0"], ["1", "0"]]),
        json!([["0", "0"], ["1/2", "0"], ["0", "0"], ["0", "1/2"]]),
    ] {
        let mut v: Value = serde_json::from_str(CASES[0]).unwrap();
        v["steps"][0]["loops"][0] = Value::Array(
            points
                .as_array()
                .unwrap()
                .iter()
                .enumerate()
                .map(|(i, p)| json!({"kind":"line","id":format!("side-{i}"),"to":p}))
                .collect(),
        );
        let recipe = recipe(&v);
        let refusal = AuthoredCap::compile(
            &recipe,
            &id("target"),
            &id("profile"),
            None,
            None,
            &mut Decider::default(),
        )
        .unwrap_err();
        assert_eq!(refusal.code, "features.unsupported-curved-limit");
    }
    let mut v: Value = serde_json::from_str(CASES[3]).unwrap();
    v["steps"][0]["plane"]["translation"][2] = json!("1/4");
    let recipe = recipe(&v);
    let refusal = AuthoredCap::compile(
        &recipe,
        &id("target"),
        &id("profile"),
        None,
        None,
        &mut Decider::default(),
    )
    .unwrap_err();
    assert_eq!(refusal.code, "features.unsupported-curved-limit");
    assert!(refusal.detail.contains("trim"));
    let recipe = Recipe::parse(CASES[0]).unwrap();
    let (result, work) = catallax_numbers::with_operation_budget(0, || {
        AuthoredCap::compile(
            &recipe,
            &id("target"),
            &id("profile"),
            None,
            None,
            &mut Decider::default(),
        )
    });
    assert!(work.exhausted);
    assert_eq!(result.unwrap_err().code, "numbers.operation-budget");
    let cap = compile(&recipe, &mut Decider::default());
    assert!(
        cap.definition()
            .construction(&[], &mut Decider::default())
            .is_err()
    );
    let query = CurvedField::new(cap, id("cap"));
    let (result, work) = catallax_numbers::with_operation_budget(0, || {
        query.location(&RPoint3::origin(), &mut Decider::default())
    });
    assert!(work.exhausted);
    assert_eq!(result.unwrap_err().code, "numbers.operation-budget");
}
