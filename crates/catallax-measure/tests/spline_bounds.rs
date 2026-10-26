//! Public application fixture regressions for represented finite spline bounds.

use catallax_measure::{Accuracy, BoundingBox, bounding_box, face_bounding_box};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_shapes::{Curve, Isometry, Pose, RPoint3, RVector3, Rotation};
use catallax_topology::{BodyId, Face, Model, SplinePrism, transformed_body};

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
