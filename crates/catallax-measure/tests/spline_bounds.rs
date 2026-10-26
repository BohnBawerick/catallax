//! Public application fixture regressions for represented finite spline bounds.

use catallax_measure::{Accuracy, BoundingBox, bounding_box, face_bounding_box};
use catallax_numbers::{Decider, MarkLog, Rational, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_shapes::{Curve, Isometry, Pose, RPoint3, RVector3, Rotation};
use catallax_topology::{
    Body, BodyId, BodyKind, Face, FaceId, Loop, Lump, Model, Shell, SplinePrism, transformed_body,
};

#[path = "spline_bounds/oracle.rs"]
mod oracle;

fn q(text: &str) -> Rational {
    Rational::parse(text).unwrap()
}

fn wrapper() -> serde_json::Value {
    serde_json::from_str(include_str!(
        "../../../.claude/skills/verify-catallax/references/examples/certified-recipe.json"
    ))
    .unwrap()
}

fn fixture() -> catallax_app::CertifiedApplicationBody {
    let wrapper = wrapper();
    let recipe = Recipe::parse(&wrapper["recipe"].to_string()).unwrap();
    catallax_app::build_certified(&recipe, None, &mut Decider::default()).unwrap()
}

fn accuracy() -> Accuracy {
    Accuracy {
        relative: q("0"),
        absolute: q("1e-10"),
    }
}

fn body_with_faces<const N: usize>(model: &mut Model, faces: [FaceId; N]) -> BodyId {
    let shell = model.add_shell(Shell {
        faces: faces.into(),
    });
    let lump = model.add_lump(Lump {
        shells: vec![shell],
    });
    model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks: MarkLog::new(),
    })
}

fn shared_cap_fixture() -> (Model, [BodyId; 2], [FaceId; 4]) {
    let application = fixture();
    let built = application.built().built();
    let mut model = built.model.clone();
    let [lower, upper, side] = SplinePrism::check(&model, built.result, &mut Decider::default())
        .unwrap()
        .faces();
    let mut cap = model.face(lower).clone();
    let Loop::Coedges(uses) = model.loop_(cap.loops[0]).clone() else {
        unreachable!()
    };
    let uses = uses
        .into_iter()
        .map(|id| {
            let coedge = model.coedge(id).clone();
            model.add_coedge(coedge)
        })
        .collect();
    cap.loops = vec![model.add_loop(Loop::Coedges(uses))];
    let origins = model.face_origins(lower).to_vec();
    let replacement = model.add_named_face(cap, origins);
    let other = body_with_faces(&mut model, [replacement, upper, side]);
    (
        model,
        [built.result, other],
        [lower, upper, side, replacement],
    )
}

fn assert_body_rim_oracle(model: &Model, body: BodyId, bounds: &BoundingBox, accuracy: &Accuracy) {
    let checked = SplinePrism::check(model, body, &mut Decider::default()).unwrap();
    let rims = checked.rings().map(|edge| {
        let Curve::BSpline(curve) = model.geom().curve(model.edge(edge).curve) else {
            unreachable!()
        };
        oracle::extrema(curve)
    });
    for axis in 0..3 {
        let min_lo = rims[0].min[axis]
            .lo
            .clone()
            .min(rims[1].min[axis].lo.clone());
        let min_hi = rims[0].min[axis]
            .hi
            .clone()
            .min(rims[1].min[axis].hi.clone());
        let max_lo = rims[0].max[axis]
            .lo
            .clone()
            .max(rims[1].max[axis].lo.clone());
        let max_hi = rims[0].max[axis]
            .hi
            .clone()
            .max(rims[1].max[axis].hi.clone());
        assert!(bounds.min[axis] <= min_lo && bounds.max[axis] >= max_hi);
        assert!(&min_hi - &bounds.min[axis] <= accuracy.absolute);
        assert!(&bounds.max[axis] - &max_lo <= accuracy.absolute);
    }
}

#[test]
fn body_bounds_accept_two_valid_prisms_sharing_a_cap_ring() {
    let (model, bodies, faces) = shared_cap_fixture();
    for body in bodies {
        assert!(model.check(body).is_empty());
        SplinePrism::check(&model, body, &mut Decider::default()).unwrap();
    }
    let answers =
        bodies.map(|body| bounding_box(&model, body, &accuracy(), &mut Decider::default()));
    eprintln!(
        "valid shared-cap body refusals: {:?}",
        answers.each_ref().map(|answer| answer.as_ref().err())
    );
    assert!(answers.iter().all(Result::is_ok));
    let [a, b] = answers.map(Result::unwrap);
    assert!(a.marks.is_empty() && b.marks.is_empty());
    assert_eq!(a.value, b.value);
    for body in bodies {
        assert_body_rim_oracle(&model, body, a.value.as_ref().unwrap(), &accuracy());
    }
    // A face alone does not select one of the model-wide prism covers.
    // Body queries above must not impose that global ambiguity on A or B.
    for face in faces {
        assert_eq!(
            face_bounding_box(&model, face, &accuracy(), &mut Decider::default())
                .unwrap_err()
                .code,
            "measure.unsupported-surface"
        );
    }
}

#[test]
fn body_scope_reads_shared_loop_uses_in_every_face_order() {
    let (mut model, [a, _], [lower, upper, side, _]) = shared_cap_fixture();
    let expected = bounding_box(&model, a, &accuracy(), &mut Decider::default())
        .unwrap()
        .value;
    let origins = model.face_origins(lower).to_vec();
    // Cross-body sharing permits reusing the loop and coedge themselves. Their
    // model-wide first owner remains A's lower cap, not this replacement face.
    let replacement = model.add_named_face(model.face(lower).clone(), origins);
    let Loop::Coedges(uses) = model.loop_(model.face(lower).loops[0]) else {
        unreachable!()
    };
    let dangling = model.coedge(uses[0]).clone();
    model.add_coedge(dangling); // Unfinished foreign topology must not affect A.
    for faces in [
        [replacement, upper, side],
        [replacement, side, upper],
        [upper, replacement, side],
        [upper, side, replacement],
        [side, replacement, upper],
        [side, upper, replacement],
    ] {
        let body = body_with_faces(&mut model, faces);
        assert!(model.check(body).is_empty());
        SplinePrism::check(&model, body, &mut Decider::default()).unwrap();
        let answer = bounding_box(&model, body, &accuracy(), &mut Decider::default()).unwrap();
        assert!(answer.marks.is_empty());
        assert_eq!(answer.value, expected);
    }
    assert_eq!(
        bounding_box(&model, a, &accuracy(), &mut Decider::default())
            .unwrap()
            .value,
        expected
    );
}

#[test]
fn body_scope_cannot_borrow_a_missing_cap_from_another_prism() {
    let (mut model, [a, _], [lower, upper, side, _]) = shared_cap_fixture();
    let expected = bounding_box(&model, a, &accuracy(), &mut Decider::default())
        .unwrap()
        .value;
    for faces in [[upper, side], [lower, side], [lower, upper]] {
        let incomplete = body_with_faces(&mut model, faces);
        assert!(SplinePrism::check(&model, incomplete, &mut Decider::default()).is_err());
        assert!(bounding_box(&model, incomplete, &accuracy(), &mut Decider::default()).is_err());
    }
    assert_eq!(
        bounding_box(&model, a, &accuracy(), &mut Decider::default())
            .unwrap()
            .value,
        expected
    );
}

#[test]
fn shared_body_queries_keep_accuracy_caches_and_model_immutable() {
    let (model, bodies, faces) = shared_cap_fixture();
    let snapshots = faces.map(|id| (model.face(id).clone(), model.face_origins(id).to_vec()));
    let counts = bodies.map(|body| model.counts(body));
    let mut tight = None;
    for tolerance in ["1e-5", "1e-12", "1e-5", "1e-12"] {
        let accuracy = Accuracy {
            relative: q("0"),
            absolute: q(tolerance),
        };
        for body in [bodies[1], bodies[0]] {
            let (answer, work) = with_operation_budget(20_000_000, || {
                bounding_box(&model, body, &accuracy, &mut Decider::default())
            });
            assert!(!work.exhausted);
            let answer = answer.unwrap();
            assert!(answer.marks.is_empty());
            let bounds = answer.value.unwrap();
            assert_body_rim_oracle(&model, body, &bounds, &accuracy);
            if tolerance == "1e-12" {
                if let Some(previous) = &tight {
                    assert_eq!(&bounds, previous);
                }
                tight = Some(bounds);
            }
        }
    }
    assert_eq!(bodies.map(|body| model.counts(body)), counts);
    assert_eq!(
        faces.map(|id| (model.face(id).clone(), model.face_origins(id).to_vec())),
        snapshots
    );
}

fn assert_prism(model: &Model, body: BodyId) -> [BoundingBox; 3] {
    let checked = SplinePrism::check(model, body, &mut Decider::default()).unwrap();
    let faces = checked.faces().map(|face| {
        let answer = face_bounding_box(model, face, &accuracy(), &mut Decider::default()).unwrap();
        assert!(answer.marks.is_empty());
        answer.value
    });
    for (index, edge) in checked.rings().into_iter().enumerate() {
        assert!(model.edge(edge).start.is_none() && model.edge(edge).end.is_none());
        let Curve::BSpline(curve) = model.geom().curve(model.edge(edge).curve) else {
            unreachable!()
        };
        oracle::check(&faces[index], &oracle::extrema(curve), &accuracy());
    }
    for axis in 0..3 {
        assert_eq!(
            faces[2].min[axis],
            faces[0].min[axis].clone().min(faces[1].min[axis].clone())
        );
        assert_eq!(
            faces[2].max[axis],
            faces[0].max[axis].clone().max(faces[1].max[axis].clone())
        );
    }
    let answer = bounding_box(model, body, &accuracy(), &mut Decider::default()).unwrap();
    assert!(answer.marks.is_empty());
    assert_eq!(answer.value.as_ref(), Some(&faces[2]));
    faces
}

#[test]
fn certified_application_body_and_all_finite_faces_answer_bounds() {
    let application = fixture();
    let built = application.built().built();
    let checked = SplinePrism::check(&built.model, built.result, &mut Decider::default()).unwrap();
    let accuracy = Accuracy {
        relative: Rational::zero(),
        absolute: Rational::parse("1e-8").unwrap(),
    };
    let (body, work) = with_operation_budget(20_000_000, || {
        bounding_box(
            &built.model,
            built.result,
            &accuracy,
            &mut Decider::default(),
        )
    });
    eprintln!("requested body bounds work: {work:?}");
    assert!(!work.exhausted);
    let faces = checked
        .faces()
        .map(|face| face_bounding_box(&built.model, face, &accuracy, &mut Decider::default()));
    eprintln!(
        "body refusal: {:?}; face refusals: {:?}",
        body.as_ref().err(),
        faces.each_ref().map(|face| face.as_ref().err())
    );
    assert!(body.is_ok() && faces.iter().all(Result::is_ok));
    assert!(body.unwrap().marks.is_empty());
    for answer in faces {
        assert!(answer.unwrap().marks.is_empty());
    }
}

#[test]
fn independent_cubic_oracle_covers_caps_and_the_represented_sweep() {
    let application = fixture();
    let built = application.built().built();
    let boxes = assert_prism(&built.model, built.result);
    // A cap's v trim excludes the extrema on the opposite end of its carrier.
    assert!(&boxes[2].max[0] - &boxes[0].max[0] > q("2"));
    assert!(&boxes[1].min[2] - &boxes[2].min[2] > q("2"));
    assert!(
        bounding_box(
            &built.model,
            built.result,
            &Accuracy::default(),
            &mut Decider::default()
        )
        .unwrap()
        .marks
        .is_empty()
    );
}

#[test]
fn exact_rotations_and_negative_coordinates_retain_interior_extrema() {
    let application = fixture();
    let built = application.built().built();
    let rotation = Rotation::from_quaternion(q("1"), q("2"), q("3"), q("4")).unwrap();
    let pose = Pose::from_rational(rotation, RVector3::new(q("-20"), q("-30"), q("-40")));
    let other = Pose::from_rational(
        Rotation::from_quaternion(q("2"), q("-3"), q("1"), q("5")).unwrap(),
        RVector3::new(q("-10"), q("-20"), q("-30")),
    );
    for map in [Isometry::from(&pose), Isometry::from(&other)] {
        let mut model = built.model.clone();
        let body = transformed_body(&mut model, built.result, &map, Clone::clone);
        assert!(model.check(body).is_empty());
        let boxes = assert_prism(&model, body);
        let checked = SplinePrism::check(&model, body, &mut Decider::default()).unwrap();
        let Curve::BSpline(curve) = model.geom().curve(model.edge(checked.rings()[0]).curve) else {
            unreachable!()
        };
        let exact = oracle::extrema(curve);
        assert!(
            (0..3).any(|axis| &exact.max[axis].lo - &exact.endpoint_max[axis].hi > q("1e-12")),
            "a global extremum must lie strictly inside a knot span"
        );
        let target = oracle::target_cap_max(checked.spline());
        assert!(
            (0..3).any(|i| target[i].hi < exact.max[i].lo || exact.max[i].hi < target[i].lo),
            "the represented and retained target extrema must be distinguished"
        );
        let error = &checked.spline().bound().upper;
        for (represented, target) in exact.max.iter().zip(target) {
            assert!(&represented.hi - &target.lo <= *error);
            assert!(&target.hi - &represented.lo <= *error);
        }
        assert!(boxes[2].max.iter().all(|value| value < &q("0")));
        for (source, mapped) in built
            .model
            .faces_of_body(built.result)
            .iter()
            .zip(model.faces_of_body(body))
        {
            assert_eq!(
                built.model.face_origins(*source),
                model.face_origins(mapped)
            );
        }
    }
}

#[test]
fn reflected_raw_copy_retains_the_complete_prism_checker_refusal() {
    let application = fixture();
    let built = application.built().built();
    let reflection = Isometry::reflection(
        &RPoint3::new(q("0"), q("0"), q("0")),
        &RVector3::new(q("1"), q("2"), q("3")),
    )
    .unwrap();
    let mut model = built.model.clone();
    let body = transformed_body(&mut model, built.result, &reflection, Clone::clone);
    assert!(model.check(body).is_empty());
    // This dependency rejects the reflected rim traversal before measurement.
    // Bounds must not manufacture a finite-domain/material certificate for it.
    assert_eq!(
        SplinePrism::check(&model, body, &mut Decider::default())
            .unwrap_err()
            .code,
        "numbers.geometric-gap-unproven"
    );
    assert!(bounding_box(&model, body, &accuracy(), &mut Decider::default()).is_err());
    for face in model.faces_of_body(body) {
        assert!(face_bounding_box(&model, face, &accuracy(), &mut Decider::default()).is_err());
    }
}

#[test]
fn scale_offset_edit_and_shorter_axial_trim_are_remeasured() {
    for (scale, distance, height) in [
        ("1/1000", "-1/4", "4"),
        ("1", "-3/10", "3"),
        ("1000", "-1/4", "4"),
    ] {
        let mut recipe = wrapper()["recipe"].clone();
        recipe["steps"][2]["size"][2] = height.into();
        recipe["steps"][4]["distance"] = distance.into();
        for path in [
            "/steps/0/loops/0/0/radius",
            "/steps/1/distance",
            "/steps/2/size/0",
            "/steps/2/size/1",
            "/steps/2/size/2",
            "/steps/2/pose/translation/0",
            "/steps/2/pose/translation/1",
            "/steps/2/pose/translation/2",
            "/steps/4/distance",
        ] {
            let value = recipe.pointer_mut(path).unwrap();
            *value = (&q(value.as_str().unwrap()) * &q(scale)).to_string().into();
        }
        let recipe = Recipe::parse(&recipe.to_string()).unwrap();
        let application =
            catallax_app::build_certified(&recipe, None, &mut Decider::default()).unwrap();
        let built = application.built().built();
        assert_prism(&built.model, built.result);
    }
}

#[test]
fn unproved_spline_face_trims_refuse_and_zero_work_preserves_the_body() {
    let application = fixture();
    let built = application.built().built();
    let mut model = built.model.clone();
    let checked = SplinePrism::check(&model, built.result, &mut Decider::default()).unwrap();
    let [cap, _, side] = checked.faces();
    let surface = model.face(side).surface;
    let bare = model.add_face(Face {
        surface,
        sense: model.face(side).sense,
        loops: vec![],
    });
    let mut partial = model.face(side).clone();
    partial.loops.pop();
    let partial = model.add_face(partial);
    // Reusing a cap loop does not bind a second face to the checked prism.
    let detached = model.add_face(model.face(cap).clone());
    for face in [bare, partial, detached] {
        assert!(face_bounding_box(&model, face, &accuracy(), &mut Decider::default()).is_err());
    }
    let counts = built.model.counts(built.result);
    let (answer, work) = with_operation_budget(0, || {
        bounding_box(
            &built.model,
            built.result,
            &accuracy(),
            &mut Decider::default(),
        )
    });
    assert_eq!(answer.unwrap_err().code, "numbers.operation-budget");
    assert!(work.exhausted);
    assert_eq!(built.model.counts(built.result), counts);
    assert!(
        bounding_box(
            &built.model,
            built.result,
            &accuracy(),
            &mut Decider::default()
        )
        .is_ok()
    );
}
