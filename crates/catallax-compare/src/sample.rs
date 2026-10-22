//! Rational primitive parameter samples, independent of the B-rep mesh.
use crate::Revision;
use catallax_numbers::Rational;
use catallax_recipe::{Operation, Origin, Recipe, Role};
use catallax_shapes::{Pose, RPoint3, RVector3, Rotation};

fn r(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn fraction(i: u32, n: u32) -> Rational {
    Rational::from_ratio(i64::from(2 * i + 1), i64::from(2 * n)).unwrap()
}
fn circle(n: u32) -> Vec<[Rational; 2]> {
    let mut points = Vec::new();
    for quarter in 0..4 {
        for i in 0..n {
            let t = fraction(i, n);
            let t2 = &t * &t;
            let denominator = &r(1) + &t2;
            let x = &(r(1) - t2) / &denominator;
            let y = r(2) * t / denominator;
            points.push(match quarter {
                0 => [x, y],
                1 => [-y, x],
                2 => [-x, -y],
                _ => [y, -x],
            });
        }
    }
    points
}
fn placement(pose: &catallax_recipe::Pose) -> Pose {
    let [w, x, y, z] = pose.rotation.each_ref().map(|n| n.value().clone());
    let rotation = Rotation::from_quaternion(w, x, y, z).expect("compiled field placement");
    let [x, y, z] = pose.translation.0.each_ref().map(|n| n.value().clone());
    Pose::new(rotation, RVector3::new(x, y, z))
}

/// Calls `visit` for each independently generated supporting-face candidate.
/// Membership in the result is always established separately by the field query.
pub(crate) fn visit(
    recipe: &Recipe,
    n: u32,
    revision: &Revision,
    mut visit: impl FnMut(RPoint3, &Origin),
) {
    let angular = circle(n);
    for step in &recipe.steps {
        if revision.is_cancelled() {
            return;
        }
        let pose = match &step.op {
            Operation::Block { pose, .. }
            | Operation::Cylinder { pose, .. }
            | Operation::Sphere { pose, .. }
            | Operation::Cone { pose, .. }
            | Operation::Torus { pose, .. } => placement(pose),
            _ => continue,
        };
        let mut emit = |xyz: [Rational; 3], role: Role| {
            if revision.is_cancelled() {
                return false;
            }
            let [x, y, z] = xyz;
            visit(
                pose.apply_point(&RPoint3::new(x, y, z)),
                &Origin {
                    step: step.id.clone(),
                    role,
                },
            );
            !revision.is_cancelled()
        };
        let wall = || Role::Side {
            segment: catallax_recipe::ElementId::new("wall").unwrap(),
            path: None,
        };
        match &step.op {
            Operation::Block { size, .. } => {
                for axis in 0..3 {
                    for upper in [false, true] {
                        for i in 0..n {
                            for j in 0..n {
                                let mut xyz = std::array::from_fn(|_| r(0));
                                xyz[axis] = if upper {
                                    size[axis].value().clone()
                                } else {
                                    r(0)
                                };
                                xyz[(axis + 1) % 3] =
                                    size[(axis + 1) % 3].value() * &fraction(i, n);
                                xyz[(axis + 2) % 3] =
                                    size[(axis + 2) % 3].value() * &fraction(j, n);
                                let role = if axis == 2 {
                                    if upper { Role::End } else { Role::Start }
                                } else {
                                    let name = match (axis, upper) {
                                        (0, false) => "x_min",
                                        (0, true) => "x_max",
                                        (1, false) => "y_min",
                                        _ => "y_max",
                                    };
                                    Role::Side {
                                        segment: catallax_recipe::ElementId::new(name).unwrap(),
                                        path: None,
                                    }
                                };
                                if !emit(xyz, role) {
                                    return;
                                }
                            }
                        }
                    }
                }
            }
            Operation::Cylinder { radius, height, .. } => round(
                radius.value(),
                radius.value(),
                height.value(),
                n,
                &angular,
                &mut emit,
            ),
            Operation::Cone {
                base_radius,
                top_radius,
                height,
                ..
            } => round(
                base_radius.value(),
                top_radius.value(),
                height.value(),
                n,
                &angular,
                &mut emit,
            ),
            Operation::Sphere { radius, .. } => {
                for [c, s] in &angular {
                    if c <= &r(0) {
                        continue;
                    }
                    for [x, y] in &angular {
                        if !emit(
                            [
                                &(radius.value() * c) * x,
                                &(radius.value() * c) * y,
                                radius.value() * s,
                            ],
                            wall(),
                        ) {
                            return;
                        }
                    }
                }
            }
            Operation::Torus { major, minor, .. } => {
                for [c, s] in &angular {
                    for [x, y] in &angular {
                        let rho = major.value() + &(minor.value() * c);
                        if !emit([&rho * x, &rho * y, minor.value() * s], wall()) {
                            return;
                        }
                    }
                }
            }
            _ => (),
        }
    }
}
fn round(
    bottom: &Rational,
    top: &Rational,
    height: &Rational,
    n: u32,
    angular: &[[Rational; 2]],
    emit: &mut impl FnMut([Rational; 3], Role) -> bool,
) {
    for i in 0..n {
        let h = fraction(i, n);
        let radius = bottom + &(&(top - bottom) * &h);
        for [x, y] in angular {
            if !emit(
                [&radius * x, &radius * y, height * &h],
                Role::Side {
                    segment: catallax_recipe::ElementId::new("wall").unwrap(),
                    path: None,
                },
            ) {
                return;
            }
        }
    }
    for (radius, z, role) in [
        (bottom, r(0), Role::Start),
        (top, height.clone(), Role::End),
    ] {
        if radius.is_zero() {
            continue;
        }
        for i in 0..n {
            for [x, y] in angular {
                let rho = radius * &fraction(i, n);
                if !emit([&rho * x, &rho * y, z.clone()], role.clone()) {
                    return;
                }
            }
        }
    }
}
