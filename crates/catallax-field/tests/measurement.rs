//! Known deviations, unsampled extrema, exact identities and refinement.

use catallax_field::measurement::{CentrePath, FaceDefinition, MeasurementOptions, measure_face};
use catallax_numbers::{Decider, Rational, Real};
use catallax_recipe::{Operation, Pose, Recipe, Step, StepId};
use catallax_shapes::{Circle, Line, Plane, Point3, Surface, SurfaceKind, Vector3};

fn r(x: i64) -> Rational {
    Rational::from_i64(x)
}
fn p(x: i64, y: i64, z: i64) -> Point3 {
    Point3::new(Real::from_i64(x), Real::from_i64(y), Real::from_i64(z))
}
fn v(x: i64, y: i64, z: i64) -> Vector3 {
    Vector3::new(Real::from_i64(x), Real::from_i64(y), Real::from_i64(z))
}
fn built(op: Operation) -> catallax_features::Built {
    let id = StepId::new("base").unwrap();
    catallax_features::evaluate(
        &Recipe::new("measurement", vec![Step { id: id.clone(), op }], id),
        &mut Decider::default(),
    )
    .unwrap()
}
fn circle() -> CentrePath {
    CentrePath::Circle(
        Circle::new(
            p(0, 0, 0),
            v(0, 0, 1),
            Real::from_i64(16),
            v(1, 0, 0),
            &mut Decider::default(),
        )
        .unwrap(),
    )
}
fn options(cells: usize) -> MeasurementOptions {
    MeasurementOptions {
        max_cells: cells,
        width: Rational::parse("0.00001").unwrap(),
        ..Default::default()
    }
}

#[test]
fn exact_and_perturbed_rolling_torus_have_two_sided_certificates() {
    // A plane at z=-1 and an outer cylinder at rho=3 have the circle
    // z=0,rho=4 as their outward unit-offset intersection.
    for (minor, deviation) in [("1", "0"), ("1.0002", "0.0002"), ("0.9998", "0.0002")] {
        let mut built = built(Operation::Torus {
            major: 4.into(),
            minor: catallax_recipe::Number::parse(minor).unwrap(),
            pose: Pose::identity(),
        });
        let face = built.model.faces_of_body(built.result)[0];
        let definition = FaceDefinition::Fillet {
            path: circle(),
            radius: Real::one(),
        };
        let result = measure_face(
            &mut built.model,
            face,
            Some(&definition),
            &options(1),
            &mut Decider::default(),
        )
        .unwrap();
        let truth = Rational::parse(deviation).unwrap();
        assert!(result.lower <= truth && truth <= result.upper);
        assert!(result.upper.clone() - result.lower.clone() < Rational::parse("1e-30").unwrap());
        if truth.is_zero() {
            assert!(result.upper.is_zero());
        }
    }
}

#[test]
fn cylinder_rolling_between_two_planes_and_exact_faces() {
    // x=-1 and y=-1 pushed outward by one meet on the z axis.
    let mut built = built(Operation::Cylinder {
        radius: 1.into(),
        height: 3.into(),
        pose: Pose::identity(),
    });
    let face = built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .find(|f| {
            built
                .model
                .geom()
                .surface(built.model.face(*f).surface)
                .kind()
                == SurfaceKind::Cylinder
        })
        .unwrap();
    let definition = FaceDefinition::Fillet {
        path: CentrePath::Line(Line::new(p(0, 0, 0), v(0, 0, 1), &mut Decider::default()).unwrap()),
        radius: Real::one(),
    };
    let result = measure_face(
        &mut built.model,
        face,
        Some(&definition),
        &options(1),
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!((result.lower, result.upper), (r(0), r(0)));
    let result = measure_face(
        &mut built.model,
        face,
        None,
        &options(1),
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!((result.lower, result.upper), (r(0), r(0)));
}

#[test]
fn whole_planar_face_bound_contains_boundary_extremum_and_tightens() {
    let mut built = built(Operation::Block {
        size: [4.into(), 2.into(), 1.into()],
        pose: Pose::identity(),
    });
    let face = built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .find(|f| {
            built
                .model
                .face_origins(*f)
                .iter()
                .any(|o| o.role == catallax_recipe::Role::Start)
        })
        .unwrap();
    let original = Surface::Plane(
        Plane::new(p(0, 0, 0), v(1, 0, 1), v(0, 1, 0), &mut Decider::default()).unwrap(),
    );
    let definition = FaceDefinition::Offset {
        original,
        distance: Real::zero(),
    };
    // On z=0, deviation is x/sqrt(2), with supremum sqrt(8) at x=4.
    // The initial interior witness alone cannot certify that maximum.
    let coarse = measure_face(
        &mut built.model,
        face,
        Some(&definition),
        &options(1),
        &mut Decider::default(),
    )
    .unwrap();
    let fine = measure_face(
        &mut built.model,
        face,
        Some(&definition),
        &options(64),
        &mut Decider::default(),
    )
    .unwrap();
    for bound in [&coarse, &fine] {
        assert!(&bound.lower * &bound.lower <= r(8));
        assert!(&bound.upper * &bound.upper >= r(8));
    }
    assert!(fine.lower > coarse.lower);
    assert!(fine.upper <= coarse.upper);
    assert!(&fine.upper - &fine.lower < &coarse.upper - &coarse.lower);
    assert!(!coarse.converged);
}

#[test]
fn translated_torus_bound_encloses_known_normal_displacement() {
    let mut pose = Pose::identity();
    pose.translation.0[2] = catallax_recipe::Number::parse("0.25").unwrap();
    let mut built = built(Operation::Torus {
        major: 4.into(),
        minor: 1.into(),
        pose,
    });
    let face = built.model.faces_of_body(built.result)[0];
    let definition = FaceDefinition::Fillet {
        path: circle(),
        radius: Real::one(),
    };
    let coarse = measure_face(
        &mut built.model,
        face,
        Some(&definition),
        &options(1),
        &mut Decider::default(),
    )
    .unwrap();
    let fine = measure_face(
        &mut built.model,
        face,
        Some(&definition),
        &options(64),
        &mut Decider::default(),
    )
    .unwrap();
    let truth = Rational::parse("0.25").unwrap();
    for bound in [&coarse, &fine] {
        assert!(bound.lower <= truth && truth <= bound.upper);
    }
    assert!(
        fine.upper < coarse.upper,
        "coarse={coarse:?}, fine={fine:?}"
    );
    assert!(fine.lower >= coarse.lower);
}
