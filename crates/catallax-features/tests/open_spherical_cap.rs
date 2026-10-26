//! Finite kept-face distance controls for one open cylindrical/spherical cell.

use std::cmp::Ordering;

use catallax_boolean::{Containment, classify_point};
use catallax_features::{BuildError, evaluate};
use catallax_field::{Location, Plan, Workspace};
use catallax_numbers::{Ball, Decider, Rational};
use catallax_recipe::{Recipe, Role, StepId};
use catallax_shapes::{Curve, Pose, RPoint3, RVector3, Rotation, Surface};
use catallax_topology::{BodyKind, Sense};
use serde_json::{Value, json};

const FIXTURE: &str = include_str!("fixtures/open-spherical-cap/base-open.json");

fn r(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn value() -> Value {
    serde_json::from_str(FIXTURE).unwrap()
}
fn parse(v: &Value) -> Recipe {
    Recipe::parse(&v.to_string()).unwrap()
}
fn square(v: &Rational) -> Rational {
    v * v
}

// Independent distance to CLOSED FINITE source faces, without CircularSkin or
// moved supports. Outside a spherical projection's polar domain its closest
// point is on the join circle; a cylindrical projection clamps to [0,4].
fn finite_material(radial: &Rational, z: &Rational) -> Location {
    let q_sq = square(radial) + square(z);
    if radial > &r("3") || z < &r("0") || q_sq > r("25") {
        return Location::Outside;
    }
    let clamped = z.clone().max(r("0")).min(r("4"));
    let cylinder_distance_sq = square(&(radial - &r("3"))) + square(&(z - &clamped));
    let cylinder = cylinder_distance_sq.cmp(&r("1/4"));
    let sphere = if z >= &r("0") && r("25") * square(z) >= &r("16") * &q_sq {
        // For a point inside the radius-5 sphere, compare (5-|p|)^2
        // with 1/4 without extracting a square root.
        r("81/4").cmp(&q_sq)
    } else {
        (square(&(radial - &r("3"))) + square(&(z - &r("4")))).cmp(&r("1/4"))
    };
    match cylinder.min(sphere) {
        Ordering::Greater => Location::Outside,
        Ordering::Equal => Location::Boundary,
        Ordering::Less => {
            if radial == &r("3") || z == &r("0") || q_sq == r("25") {
                Location::Boundary
            } else {
                Location::Inside
            }
        }
    }
}

#[test]
fn planar_base_opening_builds_in_both_independent_kernels() {
    let recipe = Recipe::parse(FIXTURE).unwrap();
    let built = evaluate(&recipe, &mut Decider::default()).unwrap();
    assert!(built.model.check(built.result).is_empty());
    Plan::compile(&recipe).unwrap();
}

#[test]
fn finite_face_oracle_checks_opening_wall_join_pole_and_material_senses() {
    let recipe = Recipe::parse(FIXTURE).unwrap();
    let mut built = evaluate(&recipe, &mut Decider::default()).unwrap();
    let plan = Plan::compile(&recipe).unwrap();
    let mut workspace = Workspace::default();
    // 375 exact rational controls, including the open plane and both sides of
    // the outer join. Agreement here supplements the whole-domain lemma.
    for iz in -1..=23 {
        for ir in 0..=14 {
            let radial = Rational::from_ratio(ir, 4).unwrap();
            let z = Rational::from_ratio(iz, 4).unwrap();
            let point = RPoint3::new(radial.clone(), r("0"), z.clone());
            let expected = finite_material(&radial, &z);
            let actual = plan
                .proven_checked(&point, &mut workspace, &mut Decider::default())
                .unwrap();
            assert_eq!(actual.location, expected, "at {point:?}");
        }
    }
    for (radial, z, expected) in [
        ("0", "0", Location::Outside),
        ("2", "1/100", Location::Outside),
        ("11/4", "0", Location::Boundary),
        ("5/2", "1", Location::Boundary),
        ("249/100", "1", Location::Outside),
        ("251/100", "1", Location::Inside),
        ("3", "4", Location::Boundary),
        ("27/10", "18/5", Location::Inside),
        ("5/2", "19/5", Location::Inside),
        ("45/26", "54/13", Location::Boundary),
        ("0", "9/2", Location::Boundary),
        ("0", "449/100", Location::Outside),
        ("0", "451/100", Location::Inside),
        ("0", "5", Location::Boundary),
        ("11/4", "-1/100", Location::Outside),
    ] {
        let point = RPoint3::new(r(radial), r("0"), r(z));
        assert_eq!(finite_material(&point.x, &point.z), expected);
        assert_eq!(
            plan.proven_checked(&point, &mut workspace, &mut Decider::default())
                .unwrap()
                .location,
            expected,
            "field at {point:?}"
        );
        if expected != Location::Boundary {
            let actual = classify_point(
                &mut built.model,
                built.result,
                &point,
                &mut Decider::default(),
            );
            assert_eq!(
                actual.proven(),
                Some(&if expected == Location::Inside {
                    Containment::Inside
                } else {
                    Containment::Outside
                }),
                "B-rep at {point:?}: {actual:?}"
            );
        }
    }
    let source = built.bodies[&StepId::new("clipped").unwrap()];
    assert_eq!(built.model.faces_of_body(source).len(), 3);
    let source_circles: Vec<_> = built
        .model
        .edges_of_body(source)
        .into_iter()
        .map(|edge| {
            let Curve::Circle(circle) = built.model.geom().curve(built.model.edge(edge).curve)
            else {
                panic!("complete source section")
            };
            (
                circle.center.z.rational_ref().cloned(),
                circle.radius_sq.rational_ref().cloned(),
            )
        })
        .collect();
    assert!(source_circles.contains(&(Some(r("0")), Some(r("9")))));
    assert!(source_circles.contains(&(Some(r("4")), Some(r("9")))));
    let faces = built.model.faces_of_body(built.result);
    assert_eq!(faces.len(), 5);
    assert_eq!(built.model.edges_of_body(built.result).len(), 4);
    assert!(built.model.vertices_of_body(built.result).is_empty());
    let mut inner = 0;
    for face in faces {
        let definition = built.model.face(face);
        let origins = built.model.face_origins(face);
        assert_eq!(origins.len(), 1, "{origins:?}");
        if let Role::Inner { of } = &origins[0].role {
            inner += 1;
            assert_eq!(origins[0].step.as_str(), "result");
            assert!(matches!(of.step.as_str(), "stock" | "limit"));
            assert_eq!(definition.sense, Sense::Reversed);
            assert!(matches!(
                built.model.geom().surface(definition.surface),
                Surface::Cylinder(_) | Surface::Sphere(_)
            ));
        } else {
            assert_eq!(definition.sense, Sense::Same);
            assert!(matches!(origins[0].step.as_str(), "stock" | "limit"));
        }
    }
    assert_eq!(inner, 2);
}

fn ball(s: &str) -> Ball {
    Ball::from_rational(&r(s), 256)
}
fn contains(actual: &catallax_measure::Enclosure, expected: Ball) {
    assert!(
        actual.lo <= expected.lo() && actual.hi >= expected.hi(),
        "{actual:?} versus {expected:?}"
    );
    assert!(actual.width() <= r("1e-8"));
}

// Independent disk integrals, in the axial interval [0,R]. The join h is
// sqrt(R^2-a^2); integrate a^2 below h and R^2-z^2 above h.
fn disk_moments(a: Ball, radius: Ball) -> [Ball; 4] {
    let h = (&radius.square() - &a.square()).sqrt().unwrap();
    let r3 = &radius.square() * &radius;
    let h3 = &h.square() * &h;
    let r5 = &r3 * &radius.square();
    let h5 = &h3 * &h.square();
    let pi = Ball::pi(256);
    [
        &pi * &ball("2/3") * (&r3 - &h3),
        &pi * &ball("1/4") * (radius.square().square() - h.square().square()),
        &pi * &(&ball("2/15") * &r5 - (&ball("1/3") * &radius.square()) * h3 + &ball("1/5") * &h5),
        &pi * &ball("2/15") * (&r5 - &h5),
    ]
}

#[test]
fn default_properties_enclose_independent_finite_disk_integrals() {
    let built = evaluate(&Recipe::parse(FIXTURE).unwrap(), &mut Decider::default()).unwrap();
    let mass = catallax_measure::mass_properties(
        &built.model,
        built.result,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert!(mass.marks.is_empty(), "{:?}", mass.marks);
    let outer = disk_moments(ball("3"), ball("5"));
    let inner = disk_moments(ball("5/2"), ball("9/2"));
    let [volume, mz, xx, zz] = std::array::from_fn(|i| &outer[i] - &inner[i]);
    let cz = &mz * &volume.recip().unwrap();
    contains(&mass.value.volume, volume.clone());
    contains(
        &mass.value.area,
        Ball::pi(256) * (ball("309/4") - ball("4") * ball("14").sqrt().unwrap()),
    );
    let centroid = mass.value.centroid.unwrap();
    contains(&centroid[0], ball("0"));
    contains(&centroid[1], ball("0"));
    contains(&centroid[2], cz.clone());
    let inertia = mass.value.inertia.unwrap();
    for (i, row) in inertia.iter().enumerate() {
        for (j, actual) in row.iter().enumerate() {
            let expected = if i != j {
                ball("0")
            } else if i == 2 {
                &ball("2") * &xx
            } else {
                &xx + &zz - &volume * &cz.square()
            };
            contains(actual, expected);
        }
    }
}

#[test]
fn rational_poses_reflected_axis_and_scales_keep_the_same_finite_material() {
    for (scale, quaternion) in [
        ("1/10", [1, 2, 3, 4]),
        ("1", [0, 1, 0, 0]),
        ("10", [1, 0, 0, 0]),
    ] {
        let mut v = value();
        let scale = r(scale);
        v["steps"][0]["radius"] = (&r("3") * &scale).to_string().into();
        v["steps"][0]["height"] = (&r("6") * &scale).to_string().into();
        v["steps"][1]["radius"] = (&r("5") * &scale).to_string().into();
        v["steps"][3]["thickness"] = (&r("1/2") * &scale).to_string().into();
        for i in [0, 1] {
            v["steps"][i]["pose"] = json!({"rotation":quaternion.map(|n| n.to_string()), "translation":["7","-2","11"]});
        }
        let [w, x, y, z] = quaternion.map(Rational::from_i64);
        let pose = Pose::from_rational(
            Rotation::from_quaternion(w, x, y, z).unwrap(),
            RVector3::new(r("7"), r("-2"), r("11")),
        );
        let recipe = parse(&v);
        let mut built = evaluate(&recipe, &mut Decider::default()).unwrap();
        assert!(built.model.check(built.result).is_empty());
        let plan = Plan::compile(&recipe).unwrap();
        let mut workspace = Workspace::default();
        for (radial, z) in [
            ("0", "1"),
            ("11/4", "1"),
            ("0", "451/100"),
            ("0", "449/100"),
            ("5/2", "1"),
        ] {
            let expected = finite_material(&r(radial), &r(z));
            let local = RPoint3::new(&r(radial) * &scale, r("0"), &r(z) * &scale);
            let point = pose.apply_rational_point(&local).unwrap();
            assert_eq!(
                plan.proven_checked(&point, &mut workspace, &mut Decider::default())
                    .unwrap()
                    .location,
                expected
            );
            if expected != Location::Boundary {
                let actual = classify_point(
                    &mut built.model,
                    built.result,
                    &point,
                    &mut Decider::default(),
                );
                assert_eq!(
                    actual.proven(),
                    Some(&if expected == Location::Inside {
                        Containment::Inside
                    } else {
                        Containment::Outside
                    }),
                    "{actual:?}"
                );
            }
        }
    }
}

#[test]
fn consumed_cavity_and_all_open_faces_regularise_without_inventing_material() {
    for thickness in ["3", "7/2"] {
        let mut v = value();
        v["steps"][3]["thickness"] = thickness.into();
        let recipe = parse(&v);
        let built = evaluate(&recipe, &mut Decider::default()).unwrap();
        assert_eq!(built.result, built.bodies[&StepId::new("clipped").unwrap()]);
        assert_eq!(built.model.faces_of_body(built.result).len(), 3);
        let plan = Plan::compile(&recipe).unwrap();
        assert_eq!(
            plan.proven_checked(
                &RPoint3::new(r("0"), r("0"), r("1")),
                &mut Workspace::default(),
                &mut Decider::default()
            )
            .unwrap()
            .location,
            Location::Inside
        );
    }
    let mut v = value();
    v["steps"][3]["open"] = json!([{"face":{"step":"stock"}}, {"face":{"step":"limit"}}]);
    let recipe = parse(&v);
    let built = evaluate(&recipe, &mut Decider::default()).unwrap();
    assert_eq!(built.model.body(built.result).kind, BodyKind::Empty);
    assert!(built.model.faces_of_body(built.result).is_empty());
    assert_eq!(
        Plan::compile(&recipe)
            .unwrap()
            .proven_checked(
                &RPoint3::new(r("0"), r("0"), r("1")),
                &mut Workspace::default(),
                &mut Decider::default()
            )
            .unwrap()
            .location,
        Location::Outside
    );
}

#[test]
fn lower_base_and_curved_openings_keep_a_typed_finite_envelope_refusal() {
    let mut lower = value();
    lower["steps"][0]["pose"]["translation"][2] = "-4".into();
    lower["steps"][0]["height"] = "10".into();
    // At (0,0,-4), a full-sphere support is only distance 1 away;
    // the kept upper spherical cap is distance sqrt(73) away. A support
    // minimum would invent material for t=3/2. This case must stay refused.
    lower["steps"][3]["thickness"] = "3/2".into();
    let mut curved = value();
    curved["steps"][3]["open"] = json!([{"face":{"step":"limit"}}]);
    for v in [lower, curved] {
        let recipe = parse(&v);
        let BuildError::Geometry { refusal, .. } =
            evaluate(&recipe, &mut Decider::default()).unwrap_err()
        else {
            panic!("expected finite-domain coverage refusal")
        };
        assert_eq!(refusal.code, "distance.cylinder-domain");
        assert_eq!(
            Plan::compile(&recipe).unwrap_err().code,
            "distance.cylinder-domain"
        );
    }
}

#[test]
fn decreasing_axis_cap_uses_the_same_finite_projection_proof() {
    let mut v = value();
    v["steps"][0]["pose"]["translation"][2] = "-6".into();
    v["steps"][3]["open"] = json!([{"face":{"step":"stock","role":"end"}}]);
    let recipe = parse(&v);
    let mut built = evaluate(&recipe, &mut Decider::default()).unwrap();
    assert!(built.model.check(built.result).is_empty());
    let plan = Plan::compile(&recipe).unwrap();
    let mut workspace = Workspace::default();
    for (radial, z) in [
        ("0", "1"),
        ("11/4", "1"),
        ("0", "451/100"),
        ("0", "449/100"),
        ("5/2", "1"),
        ("11/4", "0"),
    ] {
        let expected = finite_material(&r(radial), &r(z));
        let point = RPoint3::new(r(radial), r("0"), -r(z));
        assert_eq!(
            plan.proven_checked(&point, &mut workspace, &mut Decider::default())
                .unwrap()
                .location,
            expected
        );
        if expected != Location::Boundary {
            let actual = classify_point(
                &mut built.model,
                built.result,
                &point,
                &mut Decider::default(),
            );
            assert_eq!(
                actual.proven(),
                Some(&if expected == Location::Inside {
                    Containment::Inside
                } else {
                    Containment::Outside
                }),
                "{actual:?}"
            );
        }
    }
}

#[test]
fn raised_base_handles_inner_band_loss_exact_tie_and_axial_cavity_consumption() {
    for (base, thickness, faces, point, expected) in [
        ("3", "7/4", 4, ["0", "0", "25/8"], Location::Outside),
        ("39/10", "1/2", 4, ["0", "0", "4"], Location::Outside),
        ("39/10", "3/2", 3, ["0", "0", "4"], Location::Inside),
    ] {
        let mut v = value();
        v["steps"][0]["pose"]["translation"][2] = base.into();
        v["steps"][3]["thickness"] = thickness.into();
        let recipe = parse(&v);
        let mut built = evaluate(&recipe, &mut Decider::default()).unwrap();
        assert!(built.model.check(built.result).is_empty());
        assert_eq!(built.model.faces_of_body(built.result).len(), faces);
        let plan = Plan::compile(&recipe).unwrap();
        let [x, y, z] = point.map(r);
        let point = RPoint3::new(x, y, z);
        assert_eq!(
            plan.proven_checked(&point, &mut Workspace::default(), &mut Decider::default())
                .unwrap()
                .location,
            expected
        );
        let actual = classify_point(
            &mut built.model,
            built.result,
            &point,
            &mut Decider::default(),
        );
        assert_eq!(
            actual.proven(),
            Some(&if expected == Location::Inside {
                Containment::Inside
            } else {
                Containment::Outside
            }),
            "{actual:?}"
        );
    }
}
