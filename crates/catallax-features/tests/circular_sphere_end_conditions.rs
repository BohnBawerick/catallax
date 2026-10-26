//! Exact circular first hits, checked against a separate rational set oracle.

use catallax_features::end_condition::curved::{AuthoredCap, Branch, Location};
use catallax_field::{CurvedField, Plan, Workspace};
use catallax_numbers::{Ball, Decider, Outcome, Rational, Real};
use catallax_recipe::{Operation, Recipe, Role, StepId};
use catallax_shapes::RPoint3;
use serde_json::{Value, json};

const SOURCE: &str = include_str!("fixtures/circular-sphere-limits/sphere-entry-next.json");

fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn id(s: &str) -> StepId {
    StepId::new(s).unwrap()
}
fn recipe(v: &Value) -> Recipe {
    Recipe::parse(&v.to_string()).unwrap()
}
fn variant(to: bool, exit: bool) -> Value {
    let mut v: Value = serde_json::from_str(SOURCE).unwrap();
    if exit {
        v["steps"][0]["plane"]["translation"][2] = json!("10");
    }
    if to {
        v["steps"][2]["op"] = json!("extrude_to");
        v["steps"][2]["limit"] =
            json!({"face":{"step":"target","role":{"side":{"segment":"wall"}}}});
    }
    v
}
fn compile(recipe: &Recipe, dec: &mut Decider) -> Result<AuthoredCap, catallax_numbers::Refusal> {
    let (direction, limit) = match &recipe.steps[2].op {
        Operation::ExtrudeNext { direction, .. } => (direction.as_ref(), None),
        Operation::ExtrudeTo {
            direction, limit, ..
        } => (direction.as_ref(), Some(limit)),
        _ => unreachable!(),
    };
    AuthoredCap::compile(recipe, &id("target"), &id("profile"), direction, limit, dec)
}
fn point(cap: &AuthoredCap, u: Rational, v: Rational, t: Rational) -> RPoint3 {
    let f = cap.definition().frame();
    f.point(&[u.into(), v.into()])
        .plus(&f.direction.scaled(&t.into()))
        .as_rational()
        .unwrap()
}
fn equal(a: &Real, b: &Real) {
    assert!(matches!(
        Decider::default().sign(&(a - b)),
        Outcome::ProvenCoincident(_)
    ));
}

#[test]
fn circle_first_hit_constructs_exact_disk_and_both_branches_in_both_kernels() {
    for to in [false, true] {
        for exit in [false, true] {
            let r = recipe(&variant(to, exit));
            let mut dec = Decider::default();
            let cap = compile(&r, &mut dec).unwrap();
            assert_eq!(
                cap.definition().branch(),
                if exit { Branch::Exit } else { Branch::Entry }
            );
            for (uv, h) in [
                (["0", "0"], "5"),
                (["3", "0"], if exit { "4" } else { "6" }),
                (["9/5", "12/5"], if exit { "4" } else { "6" }),
            ] {
                equal(
                    &cap.definition()
                        .height(&uv.map(|s| q(s).into()), &mut dec)
                        .unwrap(),
                    &q(h).into(),
                );
            }
            let mut built = catallax_features::evaluate(&r, &mut dec).unwrap();
            assert!(built.model.check(built.result).is_empty());
            assert!(built.model.body(built.result).marks.is_empty());
            assert_eq!(built.model.faces_of_body(built.result).len(), 3);
            let field = CurvedField::new(cap.clone(), id("cap"));
            let plan = Plan::compile(&r).unwrap();
            let mut workspace = Workspace::default();
            for (uv, t, expected) in [
                (["0", "0"], "1", Location::Inside),
                (["3", "0"], "1", Location::Boundary),
                (["0", "0"], "0", Location::Boundary),
                (["0", "0"], "5", Location::Boundary),
                (["3", "3"], "1", Location::Outside),
                (["0", "0"], "-1", Location::Outside),
                (["0", "0"], "16", Location::Outside),
            ] {
                let p = point(&cap, q(uv[0]), q(uv[1]), q(t));
                assert_eq!(field.location(&p, &mut dec).unwrap(), expected);
                assert_eq!(
                    cap.definition()
                        .location(&uv.map(|s| q(s).into()), &q(t).into(), &mut dec)
                        .unwrap(),
                    expected
                );
                let want = match expected {
                    Location::Inside => catallax_field::Location::Inside,
                    Location::Outside => catallax_field::Location::Outside,
                    Location::Boundary => catallax_field::Location::Boundary,
                };
                assert_eq!(plan.proven(&p, &mut workspace).location, want);
                match catallax_boolean::classify_point(&mut built.model, built.result, &p, &mut dec)
                {
                    Outcome::Proven(catallax_boolean::Containment::Inside) => {
                        assert_eq!(expected, Location::Inside)
                    }
                    Outcome::Proven(catallax_boolean::Containment::Outside) => {
                        assert_eq!(expected, Location::Outside)
                    }
                    Outcome::ProvenCoincident(_) => assert_eq!(expected, Location::Boundary),
                    other => panic!("uncertified containment: {other:?}"),
                }
            }
            for (uv, t, roles) in [
                (["0", "0"], "0", vec![Role::Start]),
                (["0", "0"], "5", vec![Role::End]),
                (
                    ["3", "0"],
                    "1",
                    vec![Role::Side {
                        segment: catallax_recipe::ElementId::new("rim").unwrap(),
                        path: None,
                    }],
                ),
            ] {
                let found = field
                    .origins(&point(&cap, q(uv[0]), q(uv[1]), q(t)), &mut dec)
                    .unwrap();
                assert!(found.iter().all(|o| o.step == id("cap")));
                assert_eq!(found.into_iter().map(|o| o.role).collect::<Vec<_>>(), roles);
            }
            let names: Vec<_> = built
                .model
                .faces_of_body(built.result)
                .into_iter()
                .flat_map(|f| built.model.face_origins(f))
                .collect();
            assert_eq!(names.len(), 3);
            assert!(names.iter().all(|o| o.step == id("cap")));
            assert!(
                names.iter().any(
                    |o| matches!(&o.role,Role::Side { segment,.. } if segment.as_str()=="rim")
                )
            );
            let rim = point(&cap, q("3"), q("0"), q(if exit { "4" } else { "6" }));
            let rim_origins = field.origins(&rim, &mut dec).unwrap();
            assert_eq!(rim_origins.len(), 2);
            assert!(rim_origins.iter().any(|o| o.role == Role::End));
            assert!(
                rim_origins
                    .iter()
                    .any(|o| matches!(&o.role, Role::Side {segment,..} if segment.as_str()=="rim"))
            );
            for role in [
                Role::Start,
                Role::End,
                Role::Side {
                    segment: catallax_recipe::ElementId::new("rim").unwrap(),
                    path: None,
                },
            ] {
                let p = if role == Role::Start {
                    point(&cap, q("0"), q("0"), q("0"))
                } else if role == Role::End {
                    point(&cap, q("0"), q("0"), q("5"))
                } else {
                    point(&cap, q("3"), q("0"), q("1"))
                };
                let handles: Vec<_> = plan
                    .faces()
                    .iter()
                    .enumerate()
                    .filter(|(_, f)| f.origin.step == id("cap") && f.origin.role == role)
                    .filter_map(|(i, _)| plan.result_face(i))
                    .collect();
                assert!(!handles.is_empty());
                for handle in handles {
                    assert_eq!(
                        plan.face_query(handle, &p, &mut workspace).unwrap().state,
                        catallax_field::FaceState::OnFace
                    );
                    assert_eq!(
                        plan.face_query(
                            handle,
                            &point(&cap, q("0"), q("0"), q("1")),
                            &mut workspace
                        )
                        .unwrap()
                        .state,
                        catallax_field::FaceState::NotOnFace
                    );
                }
            }
        }
    }
}

#[test]
fn rational_disk_oracle_excludes_square_corners_far_roots_and_start_extensions() {
    for to in [false, true] {
        for exit in [false, true] {
            let r = recipe(&variant(to, exit));
            let cap = compile(&r, &mut Decider::default()).unwrap();
            let field = CurvedField::new(cap.clone(), id("cap"));
            let plan = Plan::compile(&r).unwrap();
            let mut workspace = Workspace::default();
            let mut dec = Decider::default();
            for i in -7..=7 {
                for j in -7..=7 {
                    for k in [-1, 0, 1, 4, 5, 6, 10, 15, 16] {
                        let u = Rational::from_ratio(i, 2).unwrap();
                        let v = Rational::from_ratio(j, 2).unwrap();
                        let radial = &u * &u + &v * &v;
                        let z = k - if exit { 0 } else { 10 };
                        let residual = &radial + &q(&(z * z).to_string()) - q("25");
                        let admitted = radial <= q("9")
                            && k >= 0
                            && if exit {
                                residual <= q("0")
                            } else {
                                residual >= q("0") && z < 0
                            };
                        let want = if !admitted {
                            Location::Outside
                        } else if radial == q("9") || k == 0 || residual == q("0") {
                            Location::Boundary
                        } else {
                            Location::Inside
                        };
                        let p = point(&cap, u, v, Rational::from_i64(k));
                        assert_eq!(field.location(&p, &mut dec).unwrap(), want, "{p:?}");
                        let expected = match want {
                            Location::Inside => catallax_field::Location::Inside,
                            Location::Outside => catallax_field::Location::Outside,
                            Location::Boundary => catallax_field::Location::Boundary,
                        };
                        assert_eq!(plan.proven(&p, &mut workspace).location, expected, "{p:?}");
                    }
                }
            }
        }
    }
}

#[test]
fn circular_display_step_and_requested_exact_volume_follow_disk_integrals() {
    for to in [false, true] {
        for exit in [false, true] {
            for (radius, sphere, start) in [
                ("3", "5", if exit { "10" } else { "0" }),
                ("4", "5", if exit { "11" } else { "0" }),
                ("3", "6", if exit { "9" } else { "0" }),
            ] {
                let mut v = variant(to, exit);
                v["steps"][0]["loops"][0][0]["radius"] = json!(radius);
                v["steps"][0]["plane"]["translation"][2] = json!(start);
                v["steps"][1]["radius"] = json!(sphere);
                let built =
                    catallax_features::evaluate(&recipe(&v), &mut Decider::default()).unwrap();
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
                let step = catallax_step::write_step(
                    &built.model,
                    &[built.result],
                    &Default::default(),
                    &mut Decider::default(),
                    &mut bytes,
                )
                .unwrap();
                assert!(step.marks.is_empty());
                assert_eq!(step.value.solids, 1);
                assert!(bytes.starts_with(b"ISO-10303-21;"));
                let accuracy = catallax_measure::Accuracy {
                    relative: q("0"),
                    absolute: q("1e-11"),
                };
                let measured = catallax_measure::volume(
                    &built.model,
                    built.result,
                    &accuracy,
                    &mut Decider::default(),
                )
                .unwrap();
                assert!(measured.marks.is_empty());
                let b = |s| Ball::from_rational(&q(s), 192);
                let rad = b(radius);
                let r = b(sphere);
                let h = (&r * &r - &rad * &rad).sqrt().unwrap();
                // Integral over the disk of sqrt(R²-rho²): 2pi/3 (R³-h³).
                let integral = Ball::pi(192) * b("2/3") * ((&r * &r) * r - (&h * &h) * h);
                let offset = Ball::pi(192) * (&rad * &rad) * (b("10") - b(start));
                let expected = if exit {
                    offset + integral
                } else {
                    offset - integral
                };
                assert!(
                    measured.value.lo <= expected.lo() && measured.value.hi >= expected.hi(),
                    "{radius}, {sphere}, {start}: {:?} vs {expected:?}",
                    measured.value
                );
                assert!(measured.value.width() <= q("1e-11"));
            }
        }
    }
}

#[test]
fn silhouette_branch_change_offset_oblique_and_invalid_requests_keep_refusals() {
    for change in [
        "silhouette",
        "beyond",
        "branch",
        "axis-zero-hit",
        "rim-zero-hit",
        "offset",
        "oblique",
        "annulus",
        "cylinder",
        "zero",
        "parallel",
        "missing",
    ] {
        let mut v = variant(true, false);
        let mut expected = "features.unsupported-curved-limit";
        match change {
            "silhouette" => v["steps"][0]["loops"][0][0]["radius"] = json!("5"),
            "beyond" => v["steps"][0]["loops"][0][0]["radius"] = json!("6"),
            "branch" => v["steps"][0]["plane"]["translation"][2] = json!("11/2"),
            "axis-zero-hit" => v["steps"][0]["plane"]["translation"][2] = json!("5"),
            "rim-zero-hit" => v["steps"][0]["plane"]["translation"][2] = json!("14"),
            "offset" => v["steps"][0]["loops"][0][0]["centre"] = json!(["1/10", "0"]),
            "oblique" => v["steps"][2]["direction"] = json!(["1", "0", "2"]),
            "annulus" => v["steps"][0]["loops"]
                .as_array_mut()
                .unwrap()
                .push(json!([{"kind":"circle","id":"inner","centre":["0","0"],"radius":"1"}])),
            "cylinder" => {
                v["steps"][1]["op"] = json!("cylinder");
                v["steps"][1]["height"] = json!("3");
            }
            "zero" => {
                v["steps"][2]["direction"] = json!(["0", "0", "0"]);
                expected = "extrude.zero-direction";
            }
            "parallel" => {
                v["steps"][2]["direction"] = json!(["1", "0", "0"]);
                expected = "extrude.zero-thickness";
            }
            "missing" => {
                v["steps"][2]["limit"]["touching"] =
                    json!([{"step":"target","role":{"side":{"segment":"wall"}}}]);
                expected = "recipe.missing-reference";
            }
            _ => unreachable!(),
        }
        let r = recipe(&v);
        assert_eq!(
            compile(&r, &mut Decider::default()).unwrap_err().code,
            expected,
            "{change}"
        );
        let catallax_features::BuildError::Geometry { refusal, .. } =
            catallax_features::evaluate(&r, &mut Decider::default()).unwrap_err()
        else {
            panic!("expected geometry refusal")
        };
        assert_eq!(refusal.code, expected, "{change}");
        assert_eq!(Plan::compile(&r).unwrap_err().code, expected, "{change}");
    }
}

#[test]
fn circular_compilation_build_and_query_honour_the_original_zero_work_limit() {
    let r = recipe(&variant(false, false));
    let cap = compile(&r, &mut Decider::default()).unwrap();
    let field = CurvedField::new(cap.clone(), id("cap"));
    let (result, work) =
        catallax_numbers::with_operation_budget(0, || compile(&r, &mut Decider::default()));
    assert!(work.exhausted);
    assert_eq!(result.unwrap_err().code, "numbers.operation-budget");
    let (result, work) = catallax_numbers::with_operation_budget(0, || {
        cap.build(&mut Default::default(), &id("cap"), &mut Decider::default())
    });
    assert!(work.exhausted);
    assert_eq!(result.unwrap_err().code, "numbers.operation-budget");
    let (result, work) = catallax_numbers::with_operation_budget(0, || {
        field.location(
            &RPoint3::new(q("0"), q("0"), q("1")),
            &mut Decider::default(),
        )
    });
    assert!(work.exhausted);
    assert_eq!(result.unwrap_err().code, "numbers.operation-budget");
}

#[test]
fn raw_ray_scale_reversal_proper_world_frames_and_length_scales_preserve_first_hits() {
    for to in [false, true] {
        for exit in [false, true] {
            for reverse in [false, true] {
                for scale in ["1/100", "1", "7"] {
                    let s = q(scale);
                    for placed in [false, true] {
                        let mut v = variant(to, exit);
                        let start = if exit {
                            10
                        } else if reverse {
                            20
                        } else {
                            0
                        };
                        let signed = if reverse { -&s } else { s.clone() };
                        v["steps"][2]["direction"] = json!(["0", "0", signed.to_string()]);
                        v["steps"][0]["plane"]["translation"][2] = json!(start.to_string());
                        if placed {
                            // Proper rational rotation: z -> (0,-4/5,3/5).
                            // World translation is (7,11,-3).
                            for (index, key, z) in [(0, "plane", start), (1, "pose", 10)] {
                                v["steps"][index][key] = json!({"rotation":["2","1","0","0"],"translation":[
                                    "7",(q("11")-q("4/5")*Rational::from_i64(z)).to_string(),
                                    (q("-3")+q("3/5")*Rational::from_i64(z)).to_string()]});
                            }
                            v["steps"][2]["direction"] = json!([
                                "0",
                                (-&signed * q("4/5")).to_string(),
                                (&signed * &q("3/5")).to_string()
                            ]);
                        }
                        let r = recipe(&v);
                        let cap = compile(&r, &mut Decider::default()).unwrap();
                        equal(
                            &cap.definition()
                                .height(&[Real::zero(), Real::zero()], &mut Decider::default())
                                .unwrap(),
                            &Real::from(q("5") / s.clone()),
                        );
                        let plan = Plan::compile(&r).unwrap();
                        let mut built =
                            catallax_features::evaluate(&r, &mut Decider::default()).unwrap();
                        let p = point(&cap, q("0"), q("0"), q("1") / s.clone());
                        assert_eq!(
                            plan.proven(&p, &mut Workspace::default()).location,
                            catallax_field::Location::Inside
                        );
                        assert!(matches!(
                            catallax_boolean::classify_point(
                                &mut built.model,
                                built.result,
                                &p,
                                &mut Decider::default()
                            ),
                            Outcome::Proven(catallax_boolean::Containment::Inside)
                        ));
                    }
                }
            }
        }
    }
    for scale in ["1/10", "10"] {
        let s = q(scale);
        let mut v = variant(false, false);
        v["steps"][0]["loops"][0][0]["radius"] = json!((q("3") * s.clone()).to_string());
        v["steps"][1]["radius"] = json!((q("5") * s.clone()).to_string());
        v["steps"][1]["pose"]["translation"][2] = json!((q("10") * s.clone()).to_string());
        let cap = compile(&recipe(&v), &mut Decider::default()).unwrap();
        equal(
            &cap.definition()
                .height(&[Real::zero(), Real::zero()], &mut Decider::default())
                .unwrap(),
            &Real::from(q("5") * s),
        );
        catallax_features::evaluate(&recipe(&v), &mut Decider::default()).unwrap();
    }
}

#[test]
fn radius_and_start_height_edits_keep_circular_faces_and_physical_rim_names() {
    for to in [false, true] {
        for exit in [false, true] {
            let v = variant(to, exit);
            let before = catallax_features::evaluate(&recipe(&v), &mut Decider::default()).unwrap();
            for edit in ["radius", "sphere", "height", "chart-centre"] {
                let mut edited = v.clone();
                match edit {
                    "radius" => edited["steps"][0]["loops"][0][0]["radius"] = json!("4"),
                    "sphere" => edited["steps"][1]["radius"] = json!("6"),
                    "height" => {
                        edited["steps"][0]["plane"]["translation"][2] =
                            json!(if exit { "11" } else { "-1" })
                    }
                    "chart-centre" => {
                        edited["steps"][0]["loops"][0][0]["centre"] = json!(["1", "2"]);
                        edited["steps"][0]["plane"]["translation"][0] = json!("-1");
                        edited["steps"][0]["plane"]["translation"][1] = json!("-2");
                    }
                    _ => unreachable!(),
                }
                let r = recipe(&edited);
                let after = catallax_features::evaluate(&r, &mut Decider::default()).unwrap();
                let plan = Plan::compile(&r).unwrap();
                for face in before.model.faces_of_body(before.result) {
                    for origin in before.model.face_origins(face) {
                        assert_eq!(after.model.faces_with_origin(after.result, origin).len(), 1);
                        assert!(plan.faces().iter().any(|f| &f.origin == origin));
                    }
                }
                for edge in before.model.edges_of_body(before.result) {
                    let name = before.model.edge_name(before.result, edge).unwrap();
                    assert_eq!(after.model.edges_with_name(after.result, &name).len(), 1);
                }
            }
        }
    }
}
