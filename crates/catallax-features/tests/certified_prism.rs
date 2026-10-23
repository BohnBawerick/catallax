//! Complete recipe-target/B-rep correspondence, independent of diagnostic marks.
use catallax_features::{BodyCertificate, BuildError, CertifiedBuilt, evaluate, evaluate_certified};
use catallax_numbers::{Decider, EvidenceBook, EvidenceBound, EvidenceUnit, GeometricProof, Rational, Real, Sign, with_operation_budget};
use catallax_recipe::{Recipe, StepId};
use catallax_shapes::{CurveKind, Point3, SurfaceKind};
use catallax_topology::BodyKind;
use serde_json::{Value, json};
use std::collections::BTreeSet;

fn r(value: &str) -> Rational { Rational::parse(value).unwrap() }
fn recipe(distance: &str, height: &str) -> Recipe {
    Recipe::parse(&json!({
        "format":"catallax-recipe/3", "name":"Certified inward oblique-circle prism",
        "steps":[
            {"id":"profile","op":"sketch","plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},"loops":[[{"kind":"circle","id":"wall","centre":["0","0"],"radius":"2"}]]},
            {"id":"sweep","op":"extrude","profile":"profile","distance":"10","direction":["3","0","4"]},
            {"id":"window","op":"block","size":["8","8",height],"pose":{"rotation":["3","0","1","0"],"translation":["-2","-4","4"]}},
            {"id":"source","op":"boolean","kind":"intersect","target":"sweep","tool":"window"},
            {"id":"offset","op":"offset","body":"source","distance":distance}
        ], "result":"offset"
    }).to_string()).unwrap()
}
fn build(recipe: &Recipe) -> CertifiedBuilt {
    let (result, work) = with_operation_budget(100_000_000, || evaluate_certified(recipe, &mut Decider::default()));
    eprintln!("certified prism request: {work:?}");
    let built = result.unwrap();
    assert!(!work.exhausted);
    built
}
fn brackets_target_boundary(certificate: &BodyCertificate, expected: [&str; 3]) {
    let spline = certificate.spline().expect("approximated prism support");
    let [x, y, z] = expected.map(|value| Real::from_rational(r(value)));
    // Two micrometres exceeds the hard fitting cap, independently of knot
    // phase, interpolation nodes or the producer's refinement granularity.
    let margin = Real::from_rational(r("1/500"));
    let inside = Point3::new(x.clone(), &y - &margin, z.clone());
    let outside = Point3::new(x, &y + &margin, z);
    assert_eq!(spline.side_at(&inside, &mut Decider::default()).unwrap(), Sign::Negative);
    assert_eq!(spline.side_at(&outside, &mut Decider::default()).unwrap(), Sign::Positive);
}

#[test]
fn inward_offset_builds_one_valid_solid_with_complete_cap_exceptions_and_shared_rings() {
    let recipe = recipe("-1/4", "4");
    assert!(matches!(evaluate(&recipe, &mut Decider::default()), Err(BuildError::Geometry { refusal, .. }) if refusal.code == "features.unsupported-operation"));
    let certified = build(&recipe);
    let built = certified.built();
    let certificate = certified.result_certificate();
    let counts = built.model.counts(built.result);
    assert_eq!((counts.lumps, counts.shells, counts.faces, counts.edges, counts.vertices), (1, 1, 3, 2, 0));
    assert!(built.model.check(built.result).is_empty());
    assert_eq!(certificate.record().regions.len(), 3);
    assert!(certificate.record().regions.values().all(|region| region.exception.is_some()));
    let target = certificate.target();
    assert!(Decider::default().compare(target.axial()[0], &Real::from_rational(r("9/4"))).is_coincident());
    assert!(Decider::default().compare(target.axial()[1], &Real::from_rational(r("23/4"))).is_coincident());
    assert_eq!(target.distance_mm(), &r("1/4"));
    assert_eq!(counts.edges_by_kind.get(&CurveKind::BSpline), Some(&2));
    assert_eq!(counts.faces_by_kind.get(&SurfaceKind::BSpline), Some(&1));
    brackets_target_boundary(certificate, ["27/20", "7/4", "9/5"]);
    brackets_target_boundary(certificate, ["69/20", "7/4", "23/5"]);
}

#[test]
fn a_structurally_valid_zero_error_forgery_cannot_authorize_the_fitted_body() {
    let certified = build(&recipe("-1/4", "4"));
    let certificate = certified.result_certificate();
    let mut forged = certificate.record().clone();
    forged.uniform_bound = Some(EvidenceBound { upper: Rational::zero(), unit: EvidenceUnit::Millimetres });
    for region in forged.regions.values_mut() { region.exception = None; }
    let mut evidence = EvidenceBook::new(certified.evidence().tolerance().clone());
    for (key, original) in certified.evidence().iter() {
        let mut record = original.clone();
        if let Some(bound) = &mut record.bound { bound.upper = Rational::zero(); }
        if matches!(record.proof, GeometricProof::Approximation { .. }) {
            record.proof = GeometricProof::ExactDefinition(record.identity.target);
        }
        evidence.record(*key, record, &Decider::default()).unwrap();
    }
    let inventory: BTreeSet<_> = forged.regions.keys().copied().collect();
    evidence.check_structure(&forged, &inventory, &Decider::default()).unwrap();
    assert!(certificate.verify_record(&certified.built().model, &forged, &evidence, certified.decisions(), &mut Decider::default()).is_err());
    let mut missing_cap = certificate.record().clone();
    missing_cap.regions.pop_first();
    assert!(certificate.verify_record(&certified.built().model, &missing_cap, certified.evidence(), certified.decisions(), &mut Decider::default()).is_err());
}

#[test]
fn collapsed_caps_are_exact_empty_but_nonregular_material_remains_a_named_gap() {
    let empty = build(&recipe("-1/4", "1/2"));
    assert_eq!(empty.built().model.body(empty.built().result).kind, BodyKind::Empty);
    assert!(empty.result_certificate().target().is_empty());
    assert!(empty.result_certificate().record().regions.is_empty());
    assert!(empty.result_certificate().record().uniform_bound.is_none());
    let (nonregular, work) = with_operation_budget(100_000_000, || evaluate_certified(&recipe("-7/5", "4"), &mut Decider::default()));
    assert!(!work.exhausted);
    assert!(matches!(nonregular, Err(BuildError::Geometry { refusal, .. }) if refusal.code == "features.offset-nonregular-branch"));
}

#[test]
fn a_second_offset_refits_the_original_law_instead_of_offsetting_old_controls() {
    let first = recipe("-1/4", "4");
    let mut value: Value = serde_json::to_value(&first).unwrap();
    value["steps"].as_array_mut().unwrap().push(json!({"id":"second","op":"offset","body":"offset","distance":"-1/4"}));
    value["result"] = json!("second");
    let built = build(&Recipe::parse(&value.to_string()).unwrap());
    let target = built.result_certificate().target();
    assert_eq!(target.distance_mm(), &r("1/2"));
    assert!(Decider::default().compare(target.axial()[0], &Real::from_rational(r("5/2"))).is_coincident());
    assert!(Decider::default().compare(target.axial()[1], &Real::from_rational(r("11/2"))).is_coincident());
    brackets_target_boundary(built.result_certificate(), ["3/2", "3/2", "2"]);
    let source = built.certificate(&StepId::new("offset").unwrap()).unwrap();
    let previous_origins: Vec<_> = source.faces().unwrap().map(|face| built.built().model.face_origins(face).to_vec()).into();
    let current_origins: Vec<_> = built.result_certificate().faces().unwrap().map(|face| built.built().model.face_origins(face).to_vec()).into();
    assert_eq!(previous_origins, current_origins);
}

#[test]
fn axial_planar_trim_and_later_offset_keep_the_original_exact_profile() {
    let mut value: Value = serde_json::to_value(recipe("-1/4", "4")).unwrap();
    value["steps"].as_array_mut().unwrap().extend([
        json!({"id":"trim_window","op":"block","size":["8","8","1"],"pose":{"rotation":["3","0","1","0"],"translation":["-7/5","-4","24/5"]}}),
        json!({"id":"trimmed","op":"boolean","kind":"intersect","target":"offset","tool":"trim_window"}),
    ]);
    value["result"] = json!("trimmed");
    let trimmed = build(&Recipe::parse(&value.to_string()).unwrap());
    let certificate = trimmed.result_certificate();
    assert!(Decider::default().compare(certificate.target().axial()[0], &Real::from_i64(3)).is_coincident());
    assert!(Decider::default().compare(certificate.target().axial()[1], &Real::from_i64(4)).is_coincident());
    assert_eq!(certificate.target().distance_mm(), &r("1/4"));
    brackets_target_boundary(certificate, ["9/5", "7/4", "12/5"]);
    let faces = certificate.faces().unwrap();
    assert!(trimmed.built().model.face_origins(faces[0]).iter().all(|origin| origin.step.as_str() == "trim_window"));
    assert!(trimmed.built().model.face_origins(faces[1]).iter().all(|origin| origin.step.as_str() == "trim_window"));
    assert!(trimmed.built().model.face_origins(faces[2]).iter().all(|origin| origin.step.as_str() == "sweep"));
    value["steps"].as_array_mut().unwrap().push(json!({"id":"after_trim","op":"offset","body":"trimmed","distance":"-1/4"}));
    value["result"] = json!("after_trim");
    let later = build(&Recipe::parse(&value.to_string()).unwrap());
    let target = later.result_certificate().target();
    assert_eq!(target.distance_mm(), &r("1/2"));
    assert!(Decider::default().compare(target.axial()[0], &Real::from_rational(r("13/4"))).is_coincident());
    assert!(Decider::default().compare(target.axial()[1], &Real::from_rational(r("15/4"))).is_coincident());
    brackets_target_boundary(later.result_certificate(), ["39/20", "3/2", "13/5"]);
}

#[test]
fn scalar_edit_keeps_names_but_rejects_the_previous_revision_certificate() {
    let before_recipe = recipe("-1/4", "4");
    let before = build(&before_recipe);
    let mut value: Value = serde_json::to_value(&before_recipe).unwrap();
    value["steps"][0]["loops"][0][0]["radius"] = json!("21/10");
    let after = build(&Recipe::parse(&value.to_string()).unwrap());
    brackets_target_boundary(after.result_certificate(), ["27/20", "37/20", "9/5"]);
    let names = |built: &CertifiedBuilt| built.result_certificate().faces().unwrap()
        .map(|face| built.built().model.face_origins(face).to_vec());
    assert_eq!(names(&before), names(&after));
    assert_eq!(before.built().recipe(), &before_recipe);
    assert!(before.result_certificate().verify_record(
        &after.built().model, after.result_certificate().record(),
        after.evidence(), after.decisions(), &mut Decider::default(),
    ).is_err());
}

#[test]
fn exact_empty_boolean_identity_cannot_regain_material_under_another_offset() {
    let mut value: Value = serde_json::to_value(recipe("-1/4", "4")).unwrap();
    value["steps"].as_array_mut().unwrap().extend([
        json!({"id":"nothing","op":"boolean","kind":"subtract","target":"offset","tool":"offset"}),
        json!({"id":"still_nothing","op":"offset","body":"nothing","distance":"-1/10"}),
    ]);
    value["result"] = json!("still_nothing");
    let built = build(&Recipe::parse(&value.to_string()).unwrap());
    assert!(built.result_certificate().target().is_empty());
    assert_eq!(built.built().model.body(built.built().result).kind, BodyKind::Empty);
    assert!(built.result_certificate().faces().is_none());
    assert!(built.result_certificate().record().regions.is_empty());
}

#[test]
fn spline_display_is_watertight_shared_and_within_model_relative_accuracy() {
    let options = catallax_mesh::MeshOptions {
        chord_tolerance: r("1/1000"), max_vertices: 20_000,
    };
    let (result, work) = with_operation_budget(100_000_000, || {
        let mut dec = Decider::default();
        let built = evaluate_certified(&recipe("-1/4", "4"), &mut dec)?;
        let mesh = catallax_mesh::mesh_body(
            &built.built().model, built.built().result, &options, &mut dec,
        ).map_err(|refusal| BuildError::Geometry { step: StepId::new("offset").unwrap(), refusal })?;
        Ok::<_, BuildError>((built, mesh.value))
    });
    eprintln!("certified display request: {work:?}");
    let (built, mesh) = result.unwrap();
    assert!(!work.exhausted);
    let faces: BTreeSet<_> = built.result_certificate().faces().unwrap().into_iter().collect();
    let mut uses = std::collections::BTreeMap::<_, (usize, i32, BTreeSet<_>)>::new();
    for triangle in &mesh.triangles {
        assert!(faces.contains(&triangle.face));
        for index in 0..3 {
            let a = triangle.vertices[index];
            let b = triangle.vertices[(index + 1) % 3];
            let entry = uses.entry((a.min(b), a.max(b))).or_default();
            entry.0 += 1;
            entry.1 += if a < b { 1 } else { -1 };
            entry.2.insert(triangle.face);
        }
    }
    assert!(uses.values().all(|(count, winding, _)| *count == 2 && *winding == 0));
    for edge in &mesh.edges {
        assert_eq!(edge.vertices.first(), edge.vertices.last());
        assert!(edge.max_chord_error <= options.chord_tolerance);
        for pair in edge.vertices.windows(2) {
            assert_eq!(uses[&(pair[0].min(pair[1]), pair[0].max(pair[1]))].2.len(), 2);
        }
    }
    assert!(mesh.faces.iter().all(|face| face.max_chord_error <= options.chord_tolerance));
    let error = mesh.faces.iter().map(|face| &face.max_chord_error).max().unwrap()
        + &built.result_certificate().record().uniform_bound.as_ref().unwrap().upper;
    let expected = [("27/100", "453/100"), ("-7/4", "7/4"), ("99/100", "541/100")];
    for (axis, (lo, hi)) in expected.into_iter().enumerate() {
        let actual_lo = mesh.positions.iter().map(|p| Rational::from_f64_exact(p[axis]).unwrap()).min().unwrap();
        let actual_hi = mesh.positions.iter().map(|p| Rational::from_f64_exact(p[axis]).unwrap()).max().unwrap();
        assert!((&actual_lo - &r(lo)).abs() <= error);
        assert!((&actual_hi - &r(hi)).abs() <= error);
    }
}

#[test]
fn requested_properties_measure_the_original_offset_law_not_the_fitted_mesh() {
    let accuracy = catallax_measure::Accuracy { relative: Rational::zero(), absolute: r("1e-10") };
    let (result, work) = with_operation_budget(100_000_000, || {
        let mut dec = Decider::default();
        let built = evaluate_certified(&recipe("-1/4", "4"), &mut dec)?;
        let target = built.result_certificate().target();
        catallax_measure::elliptic_offset_prism_mass_properties(
            target.definition(), target.axial(), &accuracy, &mut dec,
        ).map_err(|refusal| BuildError::Geometry { step: StepId::new("offset").unwrap(), refusal })
    });
    eprintln!("certified original-law properties request: {work:?}");
    let properties = result.unwrap().value;
    assert!(!work.exhausted);
    let inertia = properties.inertia.as_ref().unwrap();
    // Independent theta-domain composite Simpson references (131072 panels,
    // compensated binary64 sums), not values from this Taylor integrator or
    // fitted control points. A 1e-9 reference bracket covers numerical rounding.
    let reference_error = r("1e-9");
    for (value, reference) in [
        (&properties.volume, "25.94647735203738"),
        (&properties.area, "49.03510638766064"),
        (&inertia[0][0], "41.043118539348235"),
        (&inertia[1][1], "38.29601731673853"),
        (&inertia[2][2], "36.93326716845365"),
        (&inertia[0][2], "-7.0454594929621415"),
    ] {
        assert!(value.width() <= accuracy.absolute, "{value:?}");
        assert!(value.lo <= &r(reference) + &reference_error && value.hi >= &r(reference) - &reference_error, "{value:?}");
    }
    for (value, exact) in properties.centroid.as_ref().unwrap().iter().zip(["12/5", "0", "16/5"]) {
        assert!(value.contains(&r(exact)));
        assert!(value.width() <= accuracy.absolute);
    }
    assert!(inertia[0][1].contains(&Rational::zero()));
    assert!(inertia[1][2].contains(&Rational::zero()));
    assert_eq!(inertia[2][0], inertia[0][2]);
}

#[test]
fn certified_nominal_field_proves_true_distances_caps_corners_and_selected_faces() {
    use catallax_field::{CertifiedPrismField, Direction, FaceState, Location, Value as FieldValue};
    let point = |values: [&str; 3]| {
        let [x, y, z] = values.map(r);
        catallax_shapes::RPoint3::new(x, y, z)
    };
    let (result, work) = with_operation_budget(100_000_000, || {
        let mut dec = Decider::default();
        let built = evaluate_certified(&recipe("-1/4", "4"), &mut dec)?;
        let field = CertifiedPrismField::new(&built, &mut dec)?;
        for (coordinates, location, distance) in [
            (["12/5", "0", "16/5"], Location::Inside, -1.35),
            (["12/5", "7/4", "16/5"], Location::Boundary, 0.0),
            (["27/20", "0", "9/5"], Location::Boundary, 0.0),
            (["27/20", "7/4", "9/5"], Location::Boundary, 0.0),
            // Outside both factors: true product distance, not max residual.
            (["6/5", "41/20", "8/5"], Location::Outside, 0.3905124837953327),
            // Independent theta-domain nearest-point minimisation reference.
            (["12/5", "1", "17/5"], Location::Inside, -0.726634244277628),
        ] {
            let query = field.proven_checked(&point(coordinates), &mut dec)?;
            assert_eq!(query.location, location);
            let FieldValue::TrueDistance(enclosure) = query.value else { panic!("nominal metric was replaced by a residual"); };
            assert!(enclosure.lo() <= distance + 1e-12 && enclosure.hi() >= distance - 1e-12);
            assert!(enclosure.hi() - enclosure.lo() < 1e-12);
        }
        let side_point = point(["12/5", "7/4", "16/5"]);
        let side = field.proven_checked(&side_point, &mut dec)?;
        let Direction::Unique(normal) = side.direction else { panic!("smooth side requires its nominal normal"); };
        assert!(normal[0].contains(0.0) && normal[1].contains(1.0) && normal[2].contains(0.0));
        let corner_point = point(["27/20", "7/4", "9/5"]);
        assert_eq!(field.proven_checked(&corner_point, &mut dec)?.direction, Direction::Ambiguous);
        let [lower, upper, side] = built.result_certificate().faces().unwrap();
        assert_eq!(field.face_query(lower, &corner_point, &mut dec)?.state, FaceState::OnFace);
        assert_eq!(field.face_query(side, &corner_point, &mut dec)?.state, FaceState::OnFace);
        assert_eq!(field.face_query(upper, &corner_point, &mut dec)?.state, FaceState::NotOnFace);
        assert_eq!(field.face_query(lower, &side_point, &mut dec)?.state, FaceState::NotOnFace);
        Ok::<_, Box<dyn std::error::Error>>(())
    });
    eprintln!("certified nominal field request: {work:?}");
    result.unwrap();
    assert!(!work.exhausted);
}
