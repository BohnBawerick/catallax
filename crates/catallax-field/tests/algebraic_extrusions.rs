//! Independent finite sweep membership and native display, without B-rep geometry.
use catallax_field::{Location, Plan, Workspace};
use catallax_numbers::{Decider, Rational};
use catallax_recipe::{Operation, Recipe};
use catallax_shapes::{Pose, RPoint3, RVector3, Rotation};

fn r(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn recipe(profile: &str, distance: i64) -> Recipe {
    let source = match profile {
        "rectangle" => include_str!(
            "../../catallax-features/tests/fixtures/algebraic-extrusions/rectangle.json"
        ),
        "circle" => {
            include_str!("../../catallax-features/tests/fixtures/algebraic-extrusions/circle.json")
        }
        "slot" => {
            include_str!("../../catallax-features/tests/fixtures/algebraic-extrusions/slot.json")
        }
        _ => unreachable!(),
    };
    let mut recipe = Recipe::parse(source).unwrap();
    let Operation::Extrude { distance: d, .. } = &mut recipe.steps[1].op else {
        unreachable!()
    };
    *d = distance.into();
    recipe
}
fn analytic(profile: &str, distance: i64, [x, y, z]: [i64; 3]) -> Location {
    if if distance > 0 {
        z <= 0 || z >= 4
    } else {
        z >= 0 || z <= -4
    } {
        return Location::Outside;
    }
    let q = x - z;
    let value = match profile {
        "rectangle" => {
            if !(0..=4).contains(&q) || !(0..=2).contains(&y) {
                return Location::Outside;
            }
            return if q == 0 || q == 4 || y == 0 || y == 2 {
                Location::Boundary
            } else {
                Location::Inside
            };
        }
        "circle" => q * q + y * y - 4,
        "slot" => {
            if q < 0 {
                q * q + y * y - 1
            } else if q > 4 {
                (q - 4) * (q - 4) + y * y - 1
            } else {
                y * y - 1
            }
        }
        _ => unreachable!(),
    };
    match value.cmp(&0) {
        std::cmp::Ordering::Less => Location::Inside,
        std::cmp::Ordering::Equal => Location::Boundary,
        std::cmp::Ordering::Greater => Location::Outside,
    }
}
#[test]
fn signed_profiles_and_rigid_copies_agree_with_exact_analytic_containment() {
    for profile in ["rectangle", "circle", "slot"] {
        for distance in [5, -5] {
            let plan =
                Plan::compile_checked(&recipe(profile, distance), &mut Decider::default()).unwrap();
            for rotation in [[1, 0, 0, 0], [1, 2, 3, 4]] {
                let [w, x, y, z] = rotation.map(r);
                let pose = Pose::from_rational(
                    Rotation::from_quaternion(w, x, y, z).unwrap(),
                    RVector3::new(r(7), r(-3), r(11)),
                );
                let copy = plan.transformed(&pose, &mut Decider::default()).unwrap();
                let mut workspace = Workspace::default();
                for x in [-2, 0, 2, 4, 6] {
                    for y in [-1, 0, 1, 2] {
                        for z in [-4, -1, 1, 4] {
                            let point = pose
                                .apply_rational_point(&RPoint3::new(r(x), r(y), r(z)))
                                .unwrap();
                            let result = copy
                                .proven_checked(&point, &mut workspace, &mut Decider::default())
                                .unwrap();
                            assert_eq!(
                                result.location,
                                analytic(profile, distance, [x, y, z]),
                                "{profile} {distance} {rotation:?} ({x},{y},{z})"
                            );
                        }
                    }
                }
            }
        }
    }
}
#[test]
fn signed_line_circle_and_arc_profiles_have_native_certified_display() {
    for profile in ["rectangle", "circle", "slot"] {
        for distance in [5, -5] {
            let ((), work) = catallax_numbers::with_operation_budget(100_000_000, || {
                let plan =
                    Plan::compile_checked(&recipe(profile, distance), &mut Decider::default())
                        .unwrap();
                let mesh = catallax_field::display::mesh_plan(
                    &plan,
                    &Default::default(),
                    &mut Decider::default(),
                )
                .unwrap();
                assert!(!mesh.triangles.is_empty());
                assert!(mesh.max_deviation <= Rational::parse("0.01").unwrap());
                for face in plan.faces() {
                    assert!(
                        mesh.triangles
                            .iter()
                            .any(|t| plan.faces()[t.face.index()].origin == face.origin)
                    );
                }
            });
            assert!(!work.exhausted, "{profile} {distance}");
        }
    }
}
#[test]
fn impossible_directions_keep_their_distinct_refusals() {
    for (direction, code) in [
        ([0, 0, 0], "extrude.zero-direction"),
        ([1, 1, 0], "extrude.zero-thickness"),
    ] {
        let mut source = recipe("rectangle", 5);
        let Operation::Extrude { direction: d, .. } = &mut source.steps[1].op else {
            unreachable!()
        };
        *d = Some(catallax_recipe::Vector3(direction.map(Into::into)));
        assert_eq!(
            Plan::compile_checked(&source, &mut Decider::default())
                .unwrap_err()
                .code,
            code
        );
    }
}

#[test]
fn reflected_sweep_facets_keep_outward_winding_after_a_rigid_copy() {
    let plan = Plan::compile_checked(&recipe("rectangle", -5), &mut Decider::default()).unwrap();
    let pose = Pose::from_rational(
        Rotation::from_quaternion(r(1), r(2), r(3), r(4)).unwrap(),
        RVector3::new(r(7), r(-3), r(11)),
    );
    let copy = plan.transformed(&pose, &mut Decider::default()).unwrap();
    let mesh =
        catallax_field::display::mesh_plan(&copy, &Default::default(), &mut Decider::default())
            .unwrap();
    for triangle in &mesh.triangles {
        let normal = match &copy.faces()[triangle.face.index()].origin.role {
            catallax_recipe::Role::Start => [0, 0, 1],
            catallax_recipe::Role::End => [0, 0, -1],
            catallax_recipe::Role::Side { segment, .. } => match segment.as_str() {
                "a" => [0, -1, 0],
                "b" => [1, 0, -1],
                "c" => [0, 1, 0],
                "d" => [-1, 0, 1],
                _ => unreachable!(),
            },
            _ => unreachable!(),
        };
        let [x, y, z] = normal.map(r);
        let expected = pose.apply_rational_vector(&RVector3::new(x, y, z)).unwrap();
        let points = triangle.vertices.map(|i| {
            let [x, y, z] = mesh.positions[i].map(|x| Rational::from_f64_exact(x).unwrap());
            RPoint3::new(x, y, z)
        });
        let actual = points[1]
            .minus(&points[0])
            .cross(&points[2].minus(&points[0]));
        assert!(actual.dot(&expected) > r(0));
    }
    assert!(mesh.max_deviation <= Rational::parse("0.01").unwrap());
}

#[test]
fn a_signed_distance_edit_preserves_independent_field_origins() {
    let before = Plan::compile_checked(&recipe("slot", 5), &mut Decider::default()).unwrap();
    let after = Plan::compile_checked(&recipe("slot", -7), &mut Decider::default()).unwrap();
    assert_eq!(
        before.faces().iter().map(|f| &f.origin).collect::<Vec<_>>(),
        after.faces().iter().map(|f| &f.origin).collect::<Vec<_>>()
    );
}
