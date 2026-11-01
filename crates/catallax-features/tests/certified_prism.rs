//! Complete recipe-target/B-rep correspondence, independent of diagnostic marks.
use catallax_boolean::{Containment, classify_point, face_contains};
use catallax_features::{BodyCertificate, BuildError, CertifiedBuilt, evaluate_certified};
use catallax_intersect::{Contact, CurveSurface, intersect_curve_surface};
use catallax_numbers::{
    BudgetPolicy, Decider, EvidenceBook, EvidenceBound, EvidenceUnit, FactKind, GeometricProof,
    Outcome, Rational, Real, Refusal, Sign, with_operation_budget,
};
use catallax_recipe::{Recipe, Role, StepId};
use catallax_shapes::{
    Curve, CurveId, CurveKind, GeomKey, Isometry, Line, Point3, Pose, RVector3, Rotation,
    SplineSupport, SplineSurface, Surface, SurfaceKind, Vector3,
};
use catallax_topology::{BodyId, BodyKind, FaceId, Model, SplinePrism, transformed_body};
use serde_json::{Value, json};
use std::collections::BTreeSet;

fn r(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}
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
    let (result, work) = with_operation_budget(100_000_000, || {
        evaluate_certified(recipe, &mut Decider::default())
    });
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
    assert_eq!(
        spline.side_at(&inside, &mut Decider::default()).unwrap(),
        Sign::Negative
    );
    assert_eq!(
        spline.side_at(&outside, &mut Decider::default()).unwrap(),
        Sign::Positive
    );
}

#[test]
fn inward_offset_builds_one_valid_solid_with_complete_cap_exceptions_and_shared_rings() {
    let recipe = recipe("-1/4", "4");
    let certified = build(&recipe);
    let built = certified.built();
    let certificate = certified.result_certificate();
    let counts = built.model.counts(built.result);
    assert_eq!(
        (
            counts.lumps,
            counts.shells,
            counts.faces,
            counts.edges,
            counts.vertices
        ),
        (1, 1, 3, 2, 0)
    );
    assert!(built.model.check(built.result).is_empty());
    assert_eq!(certificate.record().regions.len(), 3);
    assert!(
        certificate
            .record()
            .regions
            .values()
            .all(|region| region.exception.is_some())
    );
    let target = certificate.target();
    assert!(
        Decider::default()
            .compare(target.axial()[0], &Real::from_rational(r("9/4")))
            .is_coincident()
    );
    assert!(
        Decider::default()
            .compare(target.axial()[1], &Real::from_rational(r("23/4")))
            .is_coincident()
    );
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
    forged.uniform_bound = Some(EvidenceBound {
        upper: Rational::zero(),
        unit: EvidenceUnit::Millimetres,
    });
    for region in forged.regions.values_mut() {
        region.exception = None;
    }
    let mut evidence = EvidenceBook::new(certified.evidence().tolerance().clone());
    for (key, original) in certified.evidence().iter() {
        let mut record = original.clone();
        if let Some(bound) = &mut record.bound {
            bound.upper = Rational::zero();
        }
        if matches!(record.proof, GeometricProof::Approximation { .. }) {
            record.proof = GeometricProof::ExactDefinition(record.identity.target);
        }
        evidence.record(*key, record, &Decider::default()).unwrap();
    }
    let inventory: BTreeSet<_> = forged.regions.keys().copied().collect();
    evidence
        .check_structure(&forged, &inventory, &Decider::default())
        .unwrap();
    assert!(
        certificate
            .verify_record(
                &certified.built().model,
                &forged,
                &evidence,
                certified.decisions(),
                &mut Decider::default()
            )
            .is_err()
    );
    let mut missing_cap = certificate.record().clone();
    missing_cap.regions.pop_first();
    assert!(
        certificate
            .verify_record(
                &certified.built().model,
                &missing_cap,
                certified.evidence(),
                certified.decisions(),
                &mut Decider::default()
            )
            .is_err()
    );
}

#[test]
fn collapsed_caps_are_exact_empty_but_nonregular_material_remains_a_named_gap() {
    let empty = build(&recipe("-1/4", "1/2"));
    assert_eq!(
        empty.built().model.body(empty.built().result).kind,
        BodyKind::Empty
    );
    assert!(empty.result_certificate().target().is_empty());
    assert!(empty.result_certificate().record().regions.is_empty());
    assert!(empty.result_certificate().record().uniform_bound.is_none());
    let (nonregular, work) = with_operation_budget(100_000_000, || {
        evaluate_certified(&recipe("-7/5", "4"), &mut Decider::default())
    });
    assert!(!work.exhausted);
    assert!(
        matches!(nonregular, Err(BuildError::Geometry { refusal, .. }) if refusal.code == "features.offset-nonregular-branch")
    );
}

#[test]
fn a_second_offset_refits_the_original_law_instead_of_offsetting_old_controls() {
    let first = recipe("-1/4", "4");
    let mut value: Value = serde_json::to_value(&first).unwrap();
    value["steps"]
        .as_array_mut()
        .unwrap()
        .push(json!({"id":"second","op":"offset","body":"offset","distance":"-1/4"}));
    value["result"] = json!("second");
    let built = build(&Recipe::parse(&value.to_string()).unwrap());
    let target = built.result_certificate().target();
    assert_eq!(target.distance_mm(), &r("1/2"));
    assert!(
        Decider::default()
            .compare(target.axial()[0], &Real::from_rational(r("5/2")))
            .is_coincident()
    );
    assert!(
        Decider::default()
            .compare(target.axial()[1], &Real::from_rational(r("11/2")))
            .is_coincident()
    );
    brackets_target_boundary(built.result_certificate(), ["3/2", "3/2", "2"]);
    let source = built.certificate(&StepId::new("offset").unwrap()).unwrap();
    let previous_origins: Vec<_> = source
        .faces()
        .unwrap()
        .map(|face| built.built().model.face_origins(face).to_vec())
        .into();
    let current_origins: Vec<_> = built
        .result_certificate()
        .faces()
        .unwrap()
        .map(|face| built.built().model.face_origins(face).to_vec())
        .into();
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
    assert!(
        Decider::default()
            .compare(certificate.target().axial()[0], &Real::from_i64(3))
            .is_coincident()
    );
    assert!(
        Decider::default()
            .compare(certificate.target().axial()[1], &Real::from_i64(4))
            .is_coincident()
    );
    assert_eq!(certificate.target().distance_mm(), &r("1/4"));
    brackets_target_boundary(certificate, ["9/5", "7/4", "12/5"]);
    let faces = certificate.faces().unwrap();
    assert!(
        trimmed
            .built()
            .model
            .face_origins(faces[0])
            .iter()
            .all(|origin| origin.step.as_str() == "trim_window")
    );
    assert!(
        trimmed
            .built()
            .model
            .face_origins(faces[1])
            .iter()
            .all(|origin| origin.step.as_str() == "trim_window")
    );
    assert!(
        trimmed
            .built()
            .model
            .face_origins(faces[2])
            .iter()
            .all(|origin| origin.step.as_str() == "sweep")
    );
    value["steps"]
        .as_array_mut()
        .unwrap()
        .push(json!({"id":"after_trim","op":"offset","body":"trimmed","distance":"-1/4"}));
    value["result"] = json!("after_trim");
    let later = build(&Recipe::parse(&value.to_string()).unwrap());
    let target = later.result_certificate().target();
    assert_eq!(target.distance_mm(), &r("1/2"));
    assert!(
        Decider::default()
            .compare(target.axial()[0], &Real::from_rational(r("13/4")))
            .is_coincident()
    );
    assert!(
        Decider::default()
            .compare(target.axial()[1], &Real::from_rational(r("15/4")))
            .is_coincident()
    );
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
    let names = |built: &CertifiedBuilt| {
        built
            .result_certificate()
            .faces()
            .unwrap()
            .map(|face| built.built().model.face_origins(face).to_vec())
    };
    assert_eq!(names(&before), names(&after));
    assert_eq!(before.built().recipe(), &before_recipe);
    assert!(
        before
            .result_certificate()
            .verify_record(
                &after.built().model,
                after.result_certificate().record(),
                after.evidence(),
                after.decisions(),
                &mut Decider::default(),
            )
            .is_err()
    );
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
    assert_eq!(
        built.built().model.body(built.built().result).kind,
        BodyKind::Empty
    );
    assert!(built.result_certificate().faces().is_none());
    assert!(built.result_certificate().record().regions.is_empty());
}

#[test]
fn spline_display_is_watertight_shared_and_within_model_relative_accuracy() {
    let options = catallax_mesh::MeshOptions {
        chord_tolerance: r("1/1000"),
        max_vertices: 20_000,
        max_edge_length: Some(r("10")),
        max_aspect_ratio: Some(r("1")),
        ..Default::default()
    };
    let (result, work) = with_operation_budget(100_000_000, || {
        let mut dec = Decider::default();
        let built = evaluate_certified(&recipe("-1/4", "4"), &mut dec)?;
        let mesh = catallax_mesh::mesh_body(
            &built.built().model,
            built.built().result,
            &options,
            &mut dec,
        )
        .map_err(|refusal| BuildError::Geometry {
            step: StepId::new("offset").unwrap(),
            refusal,
        })?;
        Ok::<_, BuildError>((built, mesh.value))
    });
    eprintln!("certified display request: {work:?}");
    let (built, mesh) = result.unwrap();
    assert!(!work.exhausted);
    let faces: BTreeSet<_> = built
        .result_certificate()
        .faces()
        .unwrap()
        .into_iter()
        .collect();
    let unavailable_faces: BTreeSet<_> = mesh
        .unavailable_controls
        .iter()
        .map(|control| {
            assert_eq!(control.control, catallax_mesh::MeshControl::AspectRatio);
            control.face
        })
        .collect();
    assert_eq!(unavailable_faces, faces);
    assert!(mesh.faces.iter().all(|face| {
        face.max_aspect_ratio
            .as_ref()
            .is_some_and(|ratio| *ratio > Rational::one())
    }));
    let angular = catallax_mesh::MeshOptions {
        angular_deviation: Some(r("1/100")),
        ..options.clone()
    };
    let refusal = catallax_mesh::mesh_body(
        &built.built().model,
        built.built().result,
        &angular,
        &mut Decider::default(),
    )
    .unwrap_err();
    assert_eq!(refusal.code, "mesh.spline-angular-coverage");
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
    assert!(
        uses.values()
            .all(|(count, winding, _)| *count == 2 && *winding == 0)
    );
    for edge in &mesh.edges {
        assert_eq!(edge.vertices.first(), edge.vertices.last());
        assert!(edge.max_chord_error <= options.chord_tolerance);
        for pair in edge.vertices.windows(2) {
            assert_eq!(
                uses[&(pair[0].min(pair[1]), pair[0].max(pair[1]))].2.len(),
                2
            );
        }
    }
    assert!(
        mesh.faces
            .iter()
            .all(|face| face.max_chord_error <= options.chord_tolerance)
    );
    let error = mesh
        .faces
        .iter()
        .map(|face| &face.max_chord_error)
        .max()
        .unwrap()
        + &built
            .result_certificate()
            .record()
            .uniform_bound
            .as_ref()
            .unwrap()
            .upper;
    let expected = [
        ("27/100", "453/100"),
        ("-7/4", "7/4"),
        ("99/100", "541/100"),
    ];
    for (axis, (lo, hi)) in expected.into_iter().enumerate() {
        let actual_lo = mesh
            .positions
            .iter()
            .map(|p| Rational::from_f64_exact(p[axis]).unwrap())
            .min()
            .unwrap();
        let actual_hi = mesh
            .positions
            .iter()
            .map(|p| Rational::from_f64_exact(p[axis]).unwrap())
            .max()
            .unwrap();
        assert!((&actual_lo - &r(lo)).abs() <= error);
        assert!((&actual_hi - &r(hi)).abs() <= error);
    }
}

#[test]
fn requested_properties_measure_the_original_offset_law_not_the_fitted_mesh() {
    let accuracy = catallax_measure::Accuracy {
        relative: Rational::zero(),
        absolute: r("1e-10"),
    };
    let (result, work) = with_operation_budget(100_000_000, || {
        let mut dec = Decider::default();
        let built = evaluate_certified(&recipe("-1/4", "4"), &mut dec)?;
        let target = built.result_certificate().target();
        catallax_measure::elliptic_offset_prism_mass_properties(
            target.definition(),
            target.axial(),
            &accuracy,
            &mut dec,
        )
        .map_err(|refusal| BuildError::Geometry {
            step: StepId::new("offset").unwrap(),
            refusal,
        })
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
        assert!(
            value.lo <= &r(reference) + &reference_error
                && value.hi >= &r(reference) - &reference_error,
            "{value:?}"
        );
    }
    for (value, exact) in properties
        .centroid
        .as_ref()
        .unwrap()
        .iter()
        .zip(["12/5", "0", "16/5"])
    {
        assert!(value.contains(&r(exact)));
        assert!(value.width() <= accuracy.absolute);
    }
    assert!(inertia[0][1].contains(&Rational::zero()));
    assert!(inertia[1][2].contains(&Rational::zero()));
    assert_eq!(inertia[2][0], inertia[0][2]);
}

#[test]
fn certified_nominal_field_proves_true_distances_caps_corners_and_selected_faces() {
    use catallax_field::{
        CertifiedPrismField, Direction, FaceState, Location, Value as FieldValue,
    };
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
            (
                ["6/5", "41/20", "8/5"],
                Location::Outside,
                0.3905124837953327,
            ),
            // Independent theta-domain nearest-point minimisation reference.
            (["12/5", "1", "17/5"], Location::Inside, -0.726634244277628),
        ] {
            let query = field.proven_checked(&point(coordinates), &mut dec)?;
            assert_eq!(query.location, location);
            let FieldValue::TrueDistance(enclosure) = query.value else {
                panic!("nominal metric was replaced by a residual");
            };
            assert!(enclosure.lo() <= distance + 1e-12 && enclosure.hi() >= distance - 1e-12);
            assert!(enclosure.hi() - enclosure.lo() < 1e-12);
        }
        let side_point = point(["12/5", "7/4", "16/5"]);
        let side = field.proven_checked(&side_point, &mut dec)?;
        let Direction::Unique(normal) = side.direction else {
            panic!("smooth side requires its nominal normal");
        };
        assert!(normal[0].contains(0.0) && normal[1].contains(1.0) && normal[2].contains(0.0));
        let corner_point = point(["27/20", "7/4", "9/5"]);
        assert_eq!(
            field.proven_checked(&corner_point, &mut dec)?.direction,
            Direction::Ambiguous
        );
        let [lower, upper, side] = built.result_certificate().faces().unwrap();
        assert_eq!(
            field.face_query(lower, &corner_point, &mut dec)?.state,
            FaceState::OnFace
        );
        assert_eq!(
            field.face_query(side, &corner_point, &mut dec)?.state,
            FaceState::OnFace
        );
        assert_eq!(
            field.face_query(upper, &corner_point, &mut dec)?.state,
            FaceState::NotOnFace
        );
        assert_eq!(
            field.face_query(lower, &side_point, &mut dec)?.state,
            FaceState::NotOnFace
        );
        Ok::<_, Box<dyn std::error::Error>>(())
    });
    eprintln!("certified nominal field request: {work:?}");
    result.unwrap();
    assert!(!work.exhausted);
}

// Native queries. Every expectation below is the hand derivation in the batch Q
// design, fixed before any query answered: a=(3/5,0,4/5), e=(0,1,0) and
// f=(-4/5,0,3/5) frame the sweep, P=4a is its centre, the caps sit at 2+d and
// 6-d along a, and the cardinal radius along e is b=2-d.
const QUERY_WORK: u64 = 20_000_000;

type Local = [Rational; 3];

#[derive(Clone, Copy, PartialEq)]
enum Placement {
    Original,
    Moved,
    Half,
}

struct Part {
    label: &'static str,
    placement: Placement,
    distance: &'static str,
    recipe: Recipe,
}

fn parts() -> [Part; 4] {
    let edited = recipe("-1/2", "4");
    let moved = Recipe::parse(&json!({
        "format":"catallax-recipe/3", "name":"Q-R",
        "steps":[
            {"id":"profile","op":"sketch","plane":{"rotation":["1","2","3","4"],"translation":["7","11","-13"]},"loops":[[{"kind":"circle","id":"wall","centre":["0","0"],"radius":"2"}]]},
            {"id":"sweep","op":"extrude","profile":"profile","distance":"10","direction":["14/15","14/3","23/15"]},
            {"id":"window","op":"block","size":["8","8","4"],"pose":{"rotation":["0","2","10","14"],"translation":["161/15","41/3","-253/15"]}},
            {"id":"source","op":"boolean","kind":"intersect","target":"sweep","tool":"window"},
            {"id":"offset","op":"offset","body":"source","distance":"-1/4"}
        ], "result":"offset"
    }).to_string()).unwrap();
    let half = Recipe::parse(&json!({
        "format":"catallax-recipe/3", "name":"Q-S",
        "steps":[
            {"id":"profile","op":"sketch","plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},"loops":[[{"kind":"circle","id":"wall","centre":["0","0"],"radius":"1"}]]},
            {"id":"sweep","op":"extrude","profile":"profile","distance":"5","direction":["3","0","4"]},
            {"id":"window","op":"block","size":["4","4","2"],"pose":{"rotation":["3","0","1","0"],"translation":["-1","-2","2"]}},
            {"id":"source","op":"boolean","kind":"intersect","target":"sweep","tool":"window"},
            {"id":"offset","op":"offset","body":"source","distance":"-1/8"}
        ], "result":"offset"
    }).to_string()).unwrap();
    [
        Part {
            label: "Q-N",
            placement: Placement::Original,
            distance: "1/4",
            recipe: recipe("-1/4", "4"),
        },
        Part {
            label: "Q-E",
            placement: Placement::Original,
            distance: "1/2",
            recipe: edited,
        },
        Part {
            label: "Q-R",
            placement: Placement::Moved,
            distance: "1/4",
            recipe: moved,
        },
        Part {
            label: "Q-S",
            placement: Placement::Half,
            distance: "1/4",
            recipe: half,
        },
    ]
}

impl Part {
    fn turned(&self, local: &Local) -> Local {
        const ROTATION: [[&str; 3]; 3] = [
            ["-2/3", "2/15", "11/15"],
            ["2/3", "-1/3", "2/3"],
            ["1/3", "14/15", "2/15"],
        ];
        match self.placement {
            Placement::Original => local.clone(),
            Placement::Half => local.clone().map(|value| value.mul_pow2(-1)),
            Placement::Moved => ROTATION.map(|row| {
                (0..3).fold(Rational::zero(), |sum, i| {
                    sum + r(row[i]) * local[i].clone()
                })
            }),
        }
    }
    fn point(&self, local: &Local) -> Point3 {
        let shift = match self.placement {
            Placement::Moved => [r("7"), r("11"), r("-13")],
            _ => [r("0"), r("0"), r("0")],
        };
        let [x, y, z]: [Real; 3] = std::array::from_fn(|i| {
            Real::from_rational(self.turned(local)[i].clone() + shift[i].clone())
        });
        Point3::new(x, y, z)
    }
    fn vector(&self, local: &Local) -> Vector3 {
        let [x, y, z] = self.turned(local).map(Real::from_rational);
        Vector3::new(x, y, z)
    }
}

fn along_a(height: &str) -> Local {
    let height = r(height);
    [
        r("3/5") * height.clone(),
        Rational::zero(),
        r("4/5") * height,
    ]
}
fn centre() -> Local {
    along_a("4")
}
fn axis_a() -> Local {
    along_a("1")
}
fn unit_e() -> Local {
    [r("0"), r("1"), r("0")]
}
fn unit_f() -> Local {
    [r("-4/5"), r("0"), r("3/5")]
}
fn shifted(from: &Local, by: &Local, amount: &Rational) -> Local {
    std::array::from_fn(|i| from[i].clone() + by[i].clone() * amount.clone())
}
fn built_part(part: &Part) -> (CertifiedBuilt, Model) {
    let (built, work) = with_operation_budget(QUERY_WORK, || {
        evaluate_certified(&part.recipe, &mut Decider::default())
    });
    assert!(!work.exhausted, "{}", part.label);
    let built = built.unwrap();
    let model = built.built().model.clone();
    (built, model)
}
fn support(model: &Model, body: BodyId) -> (FaceId, SplineSurface) {
    let side = SplinePrism::check(model, body, &mut Decider::default())
        .unwrap()
        .faces()[2];
    match model.geom().surface(model.face(side).surface) {
        Surface::BSpline(SplineSupport::EllipticOffset(spline)) => (side, spline.clone()),
        other => panic!("checked side has another support: {other:?}"),
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum Place {
    In,
    Out,
    On,
}

fn native_points(part: &Part) -> [(&'static str, Local, Place); 7] {
    let d = r(part.distance);
    let b = r("2") - d.clone();
    let edited = part.label == "Q-E";
    [
        ("centre", centre(), Place::In),
        (
            "edit margin",
            shifted(&centre(), &unit_e(), &r("13/8")),
            if edited { Place::Out } else { Place::In },
        ),
        (
            "beyond the profile",
            shifted(&centre(), &unit_e(), &r("3")),
            Place::Out,
        ),
        ("below the lower cap", along_a("2"), Place::Out),
        ("above the upper cap", along_a("6"), Place::Out),
        (
            "cardinal side point",
            shifted(&centre(), &unit_e(), &b),
            Place::On,
        ),
        (
            "lower cap interior",
            shifted(&[r("0"), r("0"), r("0")], &axis_a(), &(r("2") + d)),
            Place::On,
        ),
    ]
}

fn assert_binding_survives(part: &Part, built: &CertifiedBuilt, queried: &Model) {
    let original = &built.built().model;
    let body = built.built().result;
    let label = part.label;
    assert_eq!(
        format!("{:?}", original.counts(body)),
        format!("{:?}", queried.counts(body)),
        "{label}"
    );
    assert!(queried.check(body).is_empty(), "{label}");
    let before = SplinePrism::check(original, body, &mut Decider::default()).unwrap();
    let after = SplinePrism::check(queried, body, &mut Decider::default()).unwrap();
    assert_eq!(before.faces(), after.faces(), "{label}");
    assert_eq!(before.rings(), after.rings(), "{label}");
    for face in before.faces() {
        assert_eq!(
            original.face_origins(face),
            queried.face_origins(face),
            "{label}"
        );
    }
    for edge in before.rings() {
        assert_eq!(
            original.edge_constructions(edge),
            queried.edge_constructions(edge),
            "{label}"
        );
    }
    assert!(
        after
            .spline()
            .carrier()
            .same_definition(before.spline().carrier()),
        "{label}"
    );
    built
        .result_certificate()
        .verify_record(
            queried,
            built.result_certificate().record(),
            built.evidence(),
            built.decisions(),
            &mut Decider::default(),
        )
        .unwrap();
}

#[test]
fn every_native_point_answers_on_the_built_edited_moved_and_scaled_prism() {
    let mut names = Vec::new();
    for part in &parts() {
        let (built, mut model) = built_part(part);
        let body = built.built().result;
        let label = part.label;
        let counts = model.counts(body);
        assert_eq!(
            (
                counts.lumps,
                counts.shells,
                counts.faces,
                counts.edges,
                counts.vertices
            ),
            (1, 1, 3, 2, 0),
            "{label}"
        );
        assert_eq!(counts.faces_by_kind.get(&SurfaceKind::Plane), Some(&2));
        assert_eq!(counts.faces_by_kind.get(&SurfaceKind::BSpline), Some(&1));
        assert_eq!(counts.edges_by_kind.get(&CurveKind::BSpline), Some(&2));
        let (_, spline) = support(&model, body);
        let scale = if part.placement == Placement::Half {
            r("1/2")
        } else {
            Rational::one()
        };
        assert!(
            Decider::default()
                .compare(
                    spline.definition().distance(),
                    &Real::from_rational(r(part.distance) * scale)
                )
                .is_coincident(),
            "{label}: a copied label on another support"
        );
        let prism = SplinePrism::check(&model, body, &mut Decider::default()).unwrap();
        let origins = prism.faces().map(|face| model.face_origins(face).to_vec());
        let is = |origins: &[catallax_recipe::Origin], step: &str, role: &Role| {
            origins
                .iter()
                .any(|origin| origin.step.as_str() == step && origin.role == *role)
        };
        assert!(is(&origins[0], "window", &Role::Start), "{label}");
        assert!(is(&origins[1], "window", &Role::End), "{label}");
        assert!(
            origins[2]
                .iter()
                .any(|origin| origin.step.as_str() == "sweep"
                    && matches!(&origin.role, Role::Side { .. })),
            "{label}"
        );
        let rings = prism
            .rings()
            .map(|edge| model.edge_constructions(edge).to_vec());
        names.push((origins, rings));
        for (name, local, expected) in native_points(part) {
            let point = part.point(&local);
            let (answer, work) = with_operation_budget(QUERY_WORK, || {
                classify_point(&mut model, body, &point, &mut Decider::default())
            });
            assert!(!work.exhausted, "{label} {name}");
            eprintln!("{label} {name}: {work:?}");
            match (expected, &answer) {
                (Place::In, Outcome::Proven(Containment::Inside))
                | (Place::Out, Outcome::Proven(Containment::Outside))
                | (Place::On, Outcome::ProvenCoincident(_)) => {}
                _ => panic!("{label} {name}: {answer:?}, expected {expected:?}"),
            }
        }
        assert_binding_survives(part, &built, &model);
    }
    // The offset edit changes the distance and the material, never the names.
    assert_eq!(names[0], names[1]);
    assert_eq!(names[0], names[3]);
}

#[test]
fn finite_faces_keep_their_retained_domains() {
    let [canonical, ..] = parts();
    let (built, mut model) = built_part(&canonical);
    let body = built.built().result;
    let [lower, _, side] = SplinePrism::check(&model, body, &mut Decider::default())
        .unwrap()
        .faces();
    assert!(
        model
            .face_origins(lower)
            .iter()
            .any(|origin| origin.step.as_str() == "window" && origin.role == Role::Start)
    );
    for (face, xyz, expected, reason) in [
        (lower, ["27/20", "0", "9/5"], true, "centre of the cap"),
        (lower, ["27/20", "3", "9/5"], false, "plane beyond the disk"),
        (side, ["12/5", "7/4", "16/5"], true, "cardinal side point"),
        (side, ["6/5", "7/4", "8/5"], false, "below the finite side"),
    ] {
        let [x, y, z] = xyz.map(|value| Real::from_rational(r(value)));
        let point = Point3::new(x, y, z);
        let (answer, work) = with_operation_budget(QUERY_WORK, || {
            face_contains(&mut model, face, &point, &mut Decider::default())
        });
        assert!(!work.exhausted, "{reason}");
        eprintln!("face {reason}: {work:?}");
        assert_eq!(answer.unwrap(), expected, "{reason}");
    }
    assert_binding_survives(&canonical, &built, &model);
}

// The carrier's section continues past the sweep, so a point can sit exactly on
// it beyond a cap. At heights 1 and 7 this one is on no face and in no material.
#[test]
fn a_point_on_the_side_section_beyond_a_cap_is_outside() {
    let [canonical, ..] = parts();
    let (built, mut model) = built_part(&canonical);
    let body = built.built().result;
    let (side, _) = support(&model, body);
    for height in ["1", "7"] {
        let point = canonical.point(&shifted(&along_a(height), &unit_e(), &r("7/4")));
        let answer = classify_point(&mut model, body, &point, &mut Decider::default());
        assert!(
            matches!(answer, Outcome::Proven(Containment::Outside)),
            "{answer:?}"
        );
        assert!(!face_contains(&mut model, side, &point, &mut Decider::default()).unwrap());
    }
}

fn meet(
    part: &Part,
    model: &mut Model,
    surface: catallax_shapes::SurfaceId,
    through: &Local,
    direction: &Local,
    dec: &mut Decider,
) -> Result<(CurveId, CurveSurface), Refusal> {
    let line = Line::new(part.point(through), part.vector(direction), dec)?;
    let curve = model.geom_mut().add_curve(Curve::Line(line));
    let answer = intersect_curve_surface(model.geom_mut(), curve, surface, dec)?;
    assert!(answer.marks.is_empty(), "{}", part.label);
    Ok((curve, answer.value))
}

fn line_answer(
    part: &Part,
    model: &mut Model,
    surface: catallax_shapes::SurfaceId,
    through: &Local,
    direction: &Local,
) -> Result<(CurveId, CurveSurface), Refusal> {
    let (answer, work) = with_operation_budget(QUERY_WORK, || {
        meet(
            part,
            model,
            surface,
            through,
            direction,
            &mut Decider::default(),
        )
    });
    assert!(!work.exhausted, "{}", part.label);
    eprintln!("{} line: {work:?}", part.label);
    answer
}

fn zero(value: &Real) -> bool {
    matches!(Decider::default().sign(value), Outcome::ProvenCoincident(_))
}

// Complete, ordered crossings of the line through `through`, as line parameters
// along `direction`. Each records its incidence on both the line and the support.
// A root leaf has no cheap zero test, so its residual is bounded here, while the
// rows with a known answer compare it exactly.
fn crossings(
    part: &Part,
    model: &mut Model,
    surface: catallax_shapes::SurfaceId,
    through: &Local,
    direction: &Local,
    label: &str,
) -> Vec<Real> {
    let (curve, answer) = line_answer(part, model, surface, through, direction)
        .unwrap_or_else(|refusal| panic!("{label}: {refusal:?}"));
    let CurveSurface::Meets(hits) = answer else {
        panic!("{label}: {answer:?}");
    };
    let origin = part.point(through);
    let along = part.vector(direction);
    let parameters: Vec<Real> = hits
        .iter()
        .map(|hit| {
            assert_eq!(hit.contact, Contact::Crossing, "{label}");
            let facts = model.geom().facts();
            let point = GeomKey::Point(hit.point);
            assert!(
                facts.holds(FactKind::Incident, point, GeomKey::Surface(surface)),
                "{label}"
            );
            assert!(
                facts.holds(FactKind::Incident, point, GeomKey::Curve(curve)),
                "{label}"
            );
            let relative = model.geom().position(hit.point).minus(&origin);
            let off_line = Decider::default()
                .enclose_at(&relative.cross(&along).norm_squared(), 128)
                .unwrap();
            assert!(off_line.hi() < Rational::one().mul_pow2(-80), "{label}");
            &relative.dot(&along) / &along.norm_squared()
        })
        .collect();
    for pair in parameters.windows(2) {
        assert!(
            matches!(
                Decider::default().sign(&(&pair[1] - &pair[0])),
                Outcome::Proven(catallax_numbers::StrictSign::Positive)
            ),
            "{label}: line order"
        );
    }
    parameters
}

fn is_exactly(value: &Real, expected: &Rational) -> bool {
    zero(&(value - &Real::from_rational(expected.clone())))
}

fn within(value: &Real, lower: &str, upper: &str) -> bool {
    let above = Decider::default().sign(&(value - &Real::from_rational(r(lower))));
    let below = Decider::default().sign(&(&Real::from_rational(r(upper)) - value));
    matches!(
        (above, below),
        (
            Outcome::Proven(catallax_numbers::StrictSign::Positive),
            Outcome::Proven(catallax_numbers::StrictSign::Positive)
        )
    )
}

// The hit's enclosure maps onto the actual carrier, not onto the offset target.
fn on_represented_carrier(spline: &SplineSurface, point: &Point3, label: &str) {
    let mut dec = Decider::default();
    let width = r("1/1000000000000");
    let at = [
        dec.enclose(&point.x, &width).unwrap(),
        dec.enclose(&point.y, &width).unwrap(),
        dec.enclose(&point.z, &width).unwrap(),
    ];
    let support = SplineSupport::EllipticOffset(spline.clone());
    let [u, v] = support.parameters_of(&at, &dec).unwrap().unwrap();
    let image = spline.carrier().evaluate(&u, &v, &dec).unwrap();
    assert!(
        (0..3).all(|i| image[i].lo() <= at[i].hi() && at[i].lo() <= image[i].hi()),
        "{label}"
    );
}

#[test]
fn central_lines_get_both_exact_cardinal_crossings_on_every_part() {
    for part in &parts() {
        let (built, mut model) = built_part(part);
        let body = built.built().result;
        let (side, spline) = support(&model, body);
        let surface = model.face(side).surface;
        let label = part.label;
        let b = r("2") - r(part.distance);
        let parameters = crossings(part, &mut model, surface, &centre(), &unit_e(), label);
        // A placed direction keeps its local length factor, so every part reads b.
        assert_eq!(parameters.len(), 2, "{label}");
        assert!(is_exactly(&parameters[0], &(-b.clone())), "{label}");
        assert!(is_exactly(&parameters[1], &b), "{label}");
        for (parameter, sign) in [(&parameters[0], r("-1")), (&parameters[1], r("1"))] {
            let expected = part.point(&shifted(&centre(), &unit_e(), &(b.clone() * sign)));
            let hit = part
                .point(&centre())
                .plus(&part.vector(&unit_e()).scaled(parameter));
            assert!(zero(&hit.minus(&expected).norm_squared()), "{label}");
            on_represented_carrier(&spline, &hit, label);
        }
        assert_binding_survives(part, &built, &model);
    }
}

#[test]
fn every_transverse_line_row_answers_completely_on_the_represented_carrier() {
    let [canonical, ..] = parts();
    let (built, mut model) = built_part(&canonical);
    let body = built.built().result;
    let (side, spline) = support(&model, body);
    let surface = model.face(side).surface;
    let b = r("7/4");
    // Off-axis secant: the exact target has no rational root, only brackets.
    let through = shifted(&centre(), &unit_e(), &r("1/2"));
    let secant = crossings(
        &canonical,
        &mut model,
        surface,
        &through,
        &unit_f(),
        "secant",
    );
    assert_eq!(secant.len(), 2);
    assert!(within(&secant[0], "-3/2", "-1"), "first secant bracket");
    assert!(within(&secant[1], "1", "3/2"), "second secant bracket");
    let line = |parameter: &Real| {
        canonical
            .point(&through)
            .plus(&canonical.vector(&unit_f()).scaled(parameter))
    };
    for (parameter, outward) in secant.iter().zip([r("-1"), r("1")]) {
        // A rounded rational would pass every enclosure below, but is not the root.
        assert!(parameter.rational_ref().is_none(), "secant");
        on_represented_carrier(&spline, &line(parameter), "secant");
        // Rational points a nanometre either way bracket the represented crossing.
        let at = Decider::default().enclose_at(parameter, 128).unwrap().mid();
        let step = outward * r("1/1000000");
        for (point, side) in [(&at - &step, Sign::Negative), (&at + &step, Sign::Positive)] {
            let point = line(&Real::from_rational(point));
            let answer = spline.side_at(&point, &mut Decider::default()).unwrap();
            assert_eq!(answer, side, "secant");
        }
    }
    // Oblique: lambda is exactly plus or minus b, inside the finite axial interval.
    let oblique = shifted(&unit_e(), &axis_a(), &r("1/4"));
    let parameters = crossings(
        &canonical,
        &mut model,
        surface,
        &centre(),
        &oblique,
        "oblique",
    );
    assert_eq!(parameters.len(), 2);
    assert!(is_exactly(&parameters[0], &(-b.clone())), "oblique low");
    assert!(is_exactly(&parameters[1], &b), "oblique high");
    // Below the lower cap the whole line misses the finite side.
    let (_, below) =
        line_answer(&canonical, &mut model, surface, &along_a("2"), &unit_e()).unwrap();
    assert!(matches!(below, CurveSurface::Disjoint), "{below:?}");
    // Reversed direction: the same positions in the opposite order.
    let reverse = [r("0"), r("-1"), r("0")];
    let parameters = crossings(
        &canonical,
        &mut model,
        surface,
        &centre(),
        &reverse,
        "reverse",
    );
    assert_eq!(parameters.len(), 2);
    assert!(is_exactly(&parameters[0], &(-b.clone())), "reverse first");
    assert!(is_exactly(&parameters[1], &b), "reverse second");
    // A shifted origin keeps every whole-line root, both behind it.
    let origin = shifted(&centre(), &unit_e(), &r("4"));
    let parameters = crossings(
        &canonical,
        &mut model,
        surface,
        &origin,
        &unit_e(),
        "shifted",
    );
    assert_eq!(parameters.len(), 2);
    assert!(
        is_exactly(&parameters[0], &(-r("4") - b.clone())),
        "shifted"
    );
    assert!(
        is_exactly(&parameters[1], &(-r("4") + b.clone())),
        "shifted"
    );
    assert_binding_survives(&canonical, &built, &model);
}

// A chord between two interpolation nodes of a convex section meets it just
// there. The first node is the periodic seam, at both ends of the parameter.
#[test]
fn a_crossing_on_the_periodic_seam_counts_once() {
    let [canonical, ..] = parts();
    let (built, model) = built_part(&canonical);
    let (_, spline) = support(&model, built.built().result);
    let curve = spline.lower_edge();
    let seam = curve.control(0).clone();
    let far = curve.control(curve.control_count() / 2).clone();
    let raised = spline.axis().scaled(&Real::from_rational(r("3/2")));
    let (origin, end) = (seam.plus(&raised), far.plus(&raised));
    let (hits, work) = with_operation_budget(QUERY_WORK, || {
        spline.line_hits(&origin, &end.minus(&origin), &mut Decider::default())
    });
    assert!(!work.exhausted);
    let hits = hits.unwrap();
    assert_eq!(hits.len(), 2);
    assert!(zero(&hits[0].minus(&origin).norm_squared()));
    assert!(zero(&hits[1].minus(&end).norm_squared()));
}

// With a=1 and e along the line, the central row's heights 5-7/4 and 5+7/4 land
// on either side of the upper cap at 23/4, so only the first is on the carrier.
#[test]
fn a_crossing_beyond_the_sweep_is_not_a_hit() {
    let [canonical, ..] = parts();
    let (built, mut model) = built_part(&canonical);
    let surface = model.face(support(&model, built.built().result).0).surface;
    let through = along_a("5");
    let direction = shifted(&unit_e(), &axis_a(), &r("1"));
    let parameters = crossings(
        &canonical, &mut model, surface, &through, &direction, "one end",
    );
    assert_eq!(parameters.len(), 1);
    assert!(is_exactly(&parameters[0], &r("-7/4")));
}

// A finished body carried rigidly keeps its source controls, so the query runs
// in the source frame and its answers come back through the same motion.
#[test]
fn a_transported_body_answers_in_its_new_place() {
    let [canonical, ..] = parts();
    let (built, mut model) = built_part(&canonical);
    let turn = Rotation::from_quaternion(r("1"), r("2"), r("3"), r("4")).unwrap();
    let shift = RVector3::new(r("7"), r("11"), r("-13"));
    let pose = Pose::new(turn, shift.into(), &mut Decider::default()).unwrap();
    let motion = Isometry::from(&pose);
    let body = transformed_body(&mut model, built.built().result, &motion, Clone::clone);
    let middle = pose.apply_point(&canonical.point(&centre()));
    let (answer, work) = with_operation_budget(QUERY_WORK, || {
        classify_point(&mut model, body, &middle, &mut Decider::default())
    });
    assert!(!work.exhausted);
    assert!(
        matches!(answer, Outcome::Proven(Containment::Inside)),
        "{answer:?}"
    );
    let (_, spline) = support(&model, body);
    let along = pose.apply_vector(&canonical.vector(&unit_e()));
    let hits = spline
        .line_hits(&middle, &along, &mut Decider::default())
        .unwrap();
    let expected = pose.apply_point(&canonical.point(&shifted(&centre(), &unit_e(), &r("7/4"))));
    assert_eq!(hits.len(), 2);
    assert!(zero(&hits[1].minus(&expected).norm_squared()));
}

#[test]
fn a_tangent_or_an_overlapping_line_stays_a_typed_refusal() {
    let [canonical, ..] = parts();
    let (built, mut model) = built_part(&canonical);
    let body = built.built().result;
    let surface = model.face(support(&model, body).0).surface;
    let through = shifted(&centre(), &unit_e(), &r("7/4"));
    for (name, direction) in [("tangent", unit_f()), ("generator", axis_a())] {
        let answer = line_answer(&canonical, &mut model, surface, &through, &direction);
        match answer {
            Err(refusal) => assert_eq!(refusal.code, "intersect.parametric-surface", "{name}"),
            Ok((_, value)) => panic!("{name}: {value:?}"),
        }
    }
    assert_binding_survives(&canonical, &built, &model);
}

#[test]
fn the_central_solver_stays_under_the_caller_precision_and_work() {
    let [canonical, ..] = parts();
    let (built, mut model) = built_part(&canonical);
    let body = built.built().result;
    let surface = model.face(support(&model, body).0).surface;
    let mut dec = Decider::new(BudgetPolicy {
        max_work: QUERY_WORK,
        max_precision_bits: 64,
    });
    let (answer, work) = with_operation_budget(QUERY_WORK, || {
        meet(
            &canonical,
            &mut model,
            surface,
            &centre(),
            &unit_e(),
            &mut dec,
        )
    });
    assert!(!work.exhausted);
    assert!(dec.stats().max_bits <= 64, "{:?}", dec.stats());
    assert!(dec.stats().refined > 0, "{:?}", dec.stats());
    match answer {
        Ok((_, CurveSurface::Meets(hits))) => assert_eq!(hits.len(), 2),
        Ok((_, other)) => panic!("a guessed answer: {other:?}"),
        Err(refusal) => assert!(
            refusal.code == "numbers.precision-limit" || refusal.code == "numbers.operation-budget",
            "{refusal:?}"
        ),
    }
    let (answer, work) = with_operation_budget(0, || {
        meet(
            &canonical,
            &mut model,
            surface,
            &centre(),
            &unit_e(),
            &mut Decider::default(),
        )
    });
    assert!(work.exhausted);
    assert!(
        matches!(&answer, Err(refusal) if refusal.code == "numbers.operation-budget"),
        "{answer:?}"
    );
}
