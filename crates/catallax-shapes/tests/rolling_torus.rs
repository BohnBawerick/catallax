//! The exact torus a ball sweeps round a circle where two surfaces of
//! revolution about one axis meet. Each case checks the law against centres
//! and contacts worked by hand in the half-plane through the axis, and then
//! checks the published circles on the supports' own equations, so the two
//! checks share nothing with the law's profiles.

use catallax_numbers::{BudgetPolicy, Decider, Outcome, Rational, Real};
use catallax_shapes::{
    BallSide::{self, Inside, Outside},
    Circle, Cone, Cylinder, EllipticCylinder, Plane, Point3, Pose, RollingTorus, Rotation, Sphere,
    Surface, Torus, Vector3,
};

fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 100_000_000,
        ..Default::default()
    })
}
fn real(value: &str) -> Real {
    Real::from(Rational::parse(value).unwrap())
}
fn root(value: &str) -> Real {
    real(value).sqrt()
}
fn vector(v: [&str; 3]) -> Vector3 {
    let [x, y, z] = v.map(real);
    Vector3::new(x, y, z)
}
fn point(v: [&str; 3]) -> Point3 {
    let [x, y, z] = v.map(real);
    Point3::new(x, y, z)
}
fn up() -> Vector3 {
    vector(["0", "0", "1"])
}
fn plane(z: &str, normal: &str) -> Surface {
    Surface::Plane(Plane {
        origin: point(["0", "0", z]),
        normal: vector(["0", "0", normal]),
        ref_dir: vector(["1", "0", "0"]),
    })
}
fn cylinder(radius: &str) -> Surface {
    Surface::Cylinder(Cylinder {
        origin: point(["0", "0", "0"]),
        axis: up(),
        radius_sq: real(radius).square(),
        ref_dir: vector(["1", "0", "0"]),
    })
}
/// A cone with its apex at height `apex`, opening down when `down`.
fn cone(apex: &str, slope: &str, down: bool) -> Surface {
    Surface::Cone(Cone {
        apex: point(["0", "0", apex]),
        axis: vector(["0", "0", if down { "-1" } else { "1" }]),
        slope: real(slope),
        ref_dir: vector(["1", "0", "0"]),
    })
}
fn sphere(z: &str, radius: &str) -> Surface {
    Surface::Sphere(Sphere {
        center: point(["0", "0", z]),
        radius: real(radius),
        axis: up(),
        ref_dir: vector(["1", "0", "0"]),
    })
}
fn torus(major: &str, minor: &str) -> Surface {
    Surface::Torus(Torus {
        center: point(["0", "0", "0"]),
        axis: up(),
        major: real(major),
        minor: real(minor),
        ref_dir: vector(["1", "0", "0"]),
    })
}
fn ring(z: &str, radius_sq: &str) -> Circle {
    Circle {
        center: point(["0", "0", z]),
        normal: up(),
        radius_sq: real(radius_sq),
        ref_dir: vector(["1", "0", "0"]),
    }
}
fn same(a: &Real, b: &Real) -> bool {
    matches!(decider().sign(&(a - b)), Outcome::ProvenCoincident(_))
}

/// A hand-worked rolling ball: the ball centre and both contacts in the
/// half-plane through the axis, as distance from the axis and height.
struct Case {
    name: &'static str,
    supports: [Surface; 2],
    sides: [BallSide; 2],
    edge: Circle,
    radius: &'static str,
    centre: [Real; 2],
    contacts: [[Real; 2]; 2],
    ring: bool,
}

fn cases() -> Vec<Case> {
    let r10 = root("10");
    let r2 = root("2");
    let r7 = root("7");
    let r24 = root("24");
    let r5 = root("5");
    vec![
        // A frustum's top: 3 rho + z = 15 below the plane z = 6.
        Case {
            name: "frustum top",
            supports: [plane("6", "1"), cone("15", "1/3", true)],
            sides: [Inside, Inside],
            edge: ring("6", "9"),
            radius: "1",
            centre: [&(&real("10") - &r10) / &real("3"), real("5")],
            contacts: [
                [&(&real("10") - &r10) / &real("3"), real("6")],
                [
                    &(&(&real("10") - &r10) / &real("3")) + &(&real("3") / &r10),
                    &real("5") + &(&real("1") / &r10),
                ],
            ],
            ring: true,
        },
        // The top of a frustum so narrow that the centre circle is smaller
        // than the ball: the band stays on the spindle torus's outer sheet.
        Case {
            name: "narrow frustum top",
            supports: [plane("6", "1"), cone("8", "1/2", true)],
            sides: [Inside, Inside],
            edge: ring("6", "1"),
            radius: "1",
            centre: [&(&real("3") - &r5) / &real("2"), real("5")],
            contacts: [
                [&(&real("3") - &r5) / &real("2"), real("6")],
                [
                    &(&(&real("3") - &r5) / &real("2")) + &(&real("2") / &r5),
                    &real("5") + &(&real("1") / &r5),
                ],
            ],
            ring: false,
        },
        // A chamfer cone rho + z = 14 meeting the cylinder of radius 5.
        Case {
            name: "chamfer wall",
            supports: [cone("14", "1", true), cylinder("5")],
            sides: [Inside, Inside],
            edge: ring("9", "25"),
            radius: "1/2",
            centre: [real("9/2"), &real("19/2") - &(&r2 / &real("2"))],
            contacts: [
                [
                    &real("9/2") + &(&r2 / &real("4")),
                    &real("19/2") - &(&r2 / &real("4")),
                ],
                [real("5"), &real("19/2") - &(&r2 / &real("2"))],
            ],
            ring: true,
        },
        // A dome: the sphere of radius 5 above the plane z = 2.
        Case {
            name: "dome rim",
            supports: [plane("2", "-1"), sphere("0", "5")],
            sides: [Inside, Inside],
            edge: ring("2", "21"),
            radius: "1",
            centre: [r7.clone(), real("3")],
            contacts: [
                [r7.clone(), real("2")],
                [&(&r7 * &real("5")) / &real("4"), real("15/4")],
            ],
            ring: true,
        },
        // A ball seated where a sphere stands on a plane: a concave ring.
        Case {
            name: "sphere seat",
            supports: [plane("20", "1"), sphere("23", "5")],
            sides: [Outside, Outside],
            edge: ring("20", "16"),
            radius: "1",
            centre: [&real("4") * &r2, real("21")],
            contacts: [
                [&real("4") * &r2, real("20")],
                [&(&real("10") * &r2) / &real("3"), real("64/3")],
            ],
            ring: true,
        },
        // A cylinder of radius 6 capped by a sphere of radius 10.
        Case {
            name: "rod cap",
            supports: [cylinder("6"), sphere("0", "10")],
            sides: [Inside, Inside],
            edge: ring("8", "36"),
            radius: "1",
            centre: [real("5"), root("56")],
            contacts: [
                [real("6"), root("56")],
                [real("50/9"), &(&root("56") * &real("10")) / &real("9")],
            ],
            ring: true,
        },
        // A lens of two spheres of radius 5 whose centres are 6 apart.
        Case {
            name: "lens rim",
            supports: [sphere("3", "5"), sphere("-3", "5")],
            sides: [Inside, Inside],
            edge: ring("0", "16"),
            radius: "1",
            centre: [r7.clone(), real("0")],
            contacts: [
                [&(&r7 * &real("5")) / &real("4"), real("-3/4")],
                [&(&r7 * &real("5")) / &real("4"), real("3/4")],
            ],
            ring: true,
        },
        // A toroidal groove of tube radius 4 round the circle of radius 10:
        // its floor edge on the plane z = 0, and its wall edge at height 4.
        Case {
            name: "groove floor",
            supports: [plane("0", "-1"), torus("10", "4")],
            sides: [Inside, Outside],
            edge: ring("0", "36"),
            radius: "1",
            centre: [&real("10") - &r24, real("1")],
            contacts: [
                [&real("10") - &r24, real("0")],
                [
                    &real("10") - &(&(&r24 * &real("4")) / &real("5")),
                    real("4/5"),
                ],
            ],
            ring: true,
        },
        Case {
            name: "groove wall",
            supports: [torus("10", "4"), cylinder("10")],
            sides: [Outside, Inside],
            edge: ring("4", "100"),
            radius: "1",
            centre: [real("9"), r24.clone()],
            contacts: [
                [real("46/5"), &(&r24 * &real("4")) / &real("5")],
                [real("10"), r24.clone()],
            ],
            ring: true,
        },
        // The same plane meets the torus again outside it: the other crossing.
        Case {
            name: "groove outer floor",
            supports: [plane("0", "-1"), torus("10", "4")],
            sides: [Inside, Outside],
            edge: ring("0", "196"),
            radius: "1",
            centre: [&real("10") + &r24, real("1")],
            contacts: [
                [&real("10") + &r24, real("0")],
                [
                    &real("10") + &(&(&r24 * &real("4")) / &real("5")),
                    real("4/5"),
                ],
            ],
            ring: true,
        },
    ]
}

/// The point of a circle at a quarter turn from its seam direction.
fn quarter(circle: &Circle) -> Point3 {
    let across = circle.normal.cross(&circle.ref_dir);
    circle
        .center
        .plus(&across.scaled(&(&circle.radius_sq / &across.norm_squared()).sqrt()))
}
/// The point of a circle at its seam direction.
fn seam(circle: &Circle) -> Point3 {
    circle.center.plus(
        &circle
            .ref_dir
            .scaled(&(&circle.radius_sq / &circle.ref_dir.norm_squared()).sqrt()),
    )
}
/// A point of the half-plane through the axis at `x`, turned onto the x axis.
fn lifted(at: &[Real; 2]) -> Point3 {
    Point3::new(at[0].clone(), Real::zero(), at[1].clone())
}

fn check(case: &Case, built: &RollingTorus) {
    let name = case.name;
    let (center, axis, major, minor) = match (&built.surface, case.ring) {
        (Surface::Torus(t), true) => (&t.center, &t.axis, &t.major, &t.minor),
        (Surface::TorusSheet(t), false) => (&t.center, &t.axis, &t.major, &t.minor),
        (other, _) => panic!("{name}: unexpected carrier {:?}", other.kind()),
    };
    assert!(same(major, &case.centre[0]), "{name}: major");
    assert!(same(minor, &real(case.radius)), "{name}: minor");
    assert!(same(&center.z, &case.centre[1]), "{name}: centre height");
    assert!(same(&center.x, &Real::zero()) && same(&center.y, &Real::zero()));
    assert!(same(&axis.x, &Real::zero()) && same(&axis.y, &Real::zero()));
    for (index, contact) in built.contacts.iter().enumerate() {
        let want = lifted(&case.contacts[index]);
        let at = seam(contact);
        for (found, wanted) in [(&at.x, &want.x), (&at.y, &want.y), (&at.z, &want.z)] {
            assert!(
                same(found, wanted),
                "{name}: contact {index}: {} against {}",
                found.to_f64(),
                wanted.to_f64()
            );
        }
        // The contact lies on its support and on the ball about the centre
        // circle, and they share a normal there; the same holds a quarter
        // turn round, so the whole circle is a contact.
        let support = &case.supports[index];
        let equation = support.implicit().unwrap();
        for p in [seam(contact), quarter(contact)] {
            assert!(
                same(
                    &equation.eval_real(&[p.x.clone(), p.y.clone(), p.z.clone()]),
                    &Real::zero()
                ),
                "{name}: contact {index} off its support"
            );
            let radial = Vector3::new(p.x.clone(), p.y.clone(), Real::zero());
            let ball = center.plus(&radial.scaled(&(major / &radial.norm_squared().sqrt())));
            let reach = p.minus(&ball);
            assert!(
                same(&reach.norm_squared(), &real(case.radius).square()),
                "{name}: reach"
            );
            let normal = support.gradient(&p).unwrap();
            assert!(
                same(&reach.cross(&normal).norm_squared(), &Real::zero()),
                "{name}: contact {index} not tangent"
            );
        }
        assert!(
            same(&contact.normal.z, &case.edge.normal.z),
            "{name}: rail direction"
        );
    }
}

#[test]
fn coaxial_supports_roll_on_the_hand_worked_ball() {
    for case in cases() {
        let supports = [&case.supports[0], &case.supports[1]];
        let built = RollingTorus::new(
            supports,
            case.sides,
            &case.edge,
            &real(case.radius),
            &mut decider(),
        )
        .unwrap_or_else(|e| panic!("{}: {e:?}", case.name));
        check(&case, &built);
    }
}

// A moved copy of the frustum gives the moved torus: the law reads its axis
// from the edge circle, not from the coordinate axes.
#[test]
fn a_moved_pair_rolls_the_moved_ball() {
    let mut dec = decider();
    let rotation = Rotation::from_quaternion(
        Rational::from_i64(1),
        Rational::from_i64(2),
        Rational::from_i64(3),
        Rational::from_i64(4),
    )
    .unwrap();
    let pose = Pose::new(rotation, vector(["7/3", "-11/7", "13/5"]), &mut dec).unwrap();
    let case = &cases()[0];
    let supports = case.supports.each_ref().map(|s| s.transformed(&pose));
    let edge = Circle {
        center: pose.apply_point(&case.edge.center),
        normal: pose.apply_vector(&case.edge.normal),
        radius_sq: case.edge.radius_sq.clone(),
        ref_dir: pose.apply_vector(&case.edge.ref_dir),
    };
    let built = RollingTorus::new(
        [&supports[0], &supports[1]],
        case.sides,
        &edge,
        &real(case.radius),
        &mut dec,
    )
    .unwrap();
    let Surface::Torus(moved) = &built.surface else {
        panic!("expected a ring torus");
    };
    let centre = pose.apply_point(&lifted(&[Real::zero(), case.centre[1].clone()]));
    for (found, wanted) in [
        (&moved.center.x, &centre.x),
        (&moved.center.y, &centre.y),
        (&moved.center.z, &centre.z),
    ] {
        assert!(same(found, wanted));
    }
    assert!(same(&moved.major, &case.centre[0]));
    for (index, contact) in built.contacts.iter().enumerate() {
        let want = pose.apply_point(&lifted(&case.contacts[index]));
        let at = seam(contact);
        for (found, wanted) in [(&at.x, &want.x), (&at.y, &want.y), (&at.z, &want.z)] {
            assert!(same(found, wanted), "contact {index}");
        }
    }
}

// The torus's tube angle starts on its outer equator and turns toward its
// axis. A band ending on that equator keeps to the axis side of it whichever
// way the edge runs; a band across it keeps the edge circle's normal.
#[test]
fn a_band_on_the_outer_equator_turns_toward_the_axis() {
    let axis = |case: &Case, normal: &str| -> Real {
        let edge = Circle {
            normal: vector(["0", "0", normal]),
            ..case.edge.clone()
        };
        let built = RollingTorus::new(
            [&case.supports[0], &case.supports[1]],
            case.sides,
            &edge,
            &real(case.radius),
            &mut decider(),
        )
        .unwrap();
        assert!(
            built
                .contacts
                .iter()
                .all(|c| same(&c.normal.z, &real(normal)))
        );
        let Surface::Torus(torus) = &built.surface else {
            panic!("expected a ring torus");
        };
        torus.axis.z.clone()
    };
    let cases = cases();
    // The wall contact lies on the outer equator, the groove contact below.
    let wall = &cases[8];
    assert_eq!(wall.name, "groove wall");
    for normal in ["1", "-1"] {
        assert!(same(&axis(wall, normal), &real("-1")));
    }
    // The dome's band runs from below its centre to above it.
    let dome = &cases[3];
    assert_eq!(dome.name, "dome rim");
    for normal in ["1", "-1"] {
        assert!(same(&axis(dome, normal), &real(normal)));
    }
}

fn refusal(supports: [&Surface; 2], sides: [BallSide; 2], edge: &Circle, radius: &str) -> String {
    let refusal = RollingTorus::new(supports, sides, edge, &real(radius), &mut decider())
        .expect_err("expected a refusal");
    format!("{}: {}", refusal.code, refusal.detail)
}

#[test]
fn unsupported_pairs_refuse_with_their_reason() {
    let dome = [plane("2", "-1"), sphere("0", "5")];
    let dome = [&dome[0], &dome[1]];
    assert_eq!(
        refusal(dome, [Inside, Inside], &ring("2", "21"), "2"),
        "shapes.rolling-torus: the offsets of the supports stop crossing before this radius"
    );
    assert_eq!(
        refusal(dome, [Inside, Inside], &ring("2", "20"), "1"),
        "shapes.rolling-torus-edge: the edge circle does not lie on both supports"
    );
    assert_eq!(
        refusal(dome, [Inside, Inside], &ring("2", "21"), "0"),
        "shapes.rolling-torus: the radius is not positive"
    );
    let off_axis = Surface::Sphere(Sphere {
        center: point(["1", "0", "0"]),
        radius: real("5"),
        axis: up(),
        ref_dir: vector(["1", "0", "0"]),
    });
    assert_eq!(
        refusal(
            [dome[0], &off_axis],
            [Inside, Inside],
            &ring("2", "21"),
            "1"
        ),
        "shapes.rolling-torus: the supports do not turn about the edge circle's axis"
    );
    let oval = Surface::EllipticCylinder(EllipticCylinder {
        origin: point(["0", "0", "0"]),
        axis: up(),
        major_dir: vector(["1", "0", "0"]),
        major_sq: real("25"),
        minor_sq: real("16"),
    });
    assert_eq!(
        refusal([dome[0], &oval], [Inside, Inside], &ring("2", "21"), "1"),
        "shapes.rolling-torus: a support is not a surface of revolution"
    );
    // A sphere and a cylinder of its radius touch along the equator.
    let (wall, ball) = (cylinder("5"), sphere("0", "5"));
    assert_eq!(
        refusal([&wall, &ball], [Inside, Inside], &ring("0", "25"), "1"),
        "shapes.rolling-torus: the supports are tangent along the edge"
    );
    // A frustum too narrow for its ball: 3 rho + 2 z = 15 below z = 6.
    let narrow = [plane("6", "1"), cone("15/2", "2/3", true)];
    assert_eq!(
        refusal(
            [&narrow[0], &narrow[1]],
            [Inside, Inside],
            &ring("6", "1"),
            "2"
        ),
        "shapes.rolling-torus: the ball centre reaches the axis"
    );
    // A pointed cone standing on a plane: past the apex the ball's contact
    // with the cone would cross the axis.
    let spike = [plane("0", "1"), cone("4", "3/4", true)];
    assert_eq!(
        refusal(
            [&spike[0], &spike[1]],
            [Outside, Outside],
            &ring("0", "9"),
            "12"
        ),
        "shapes.rolling-torus: a contact reaches the axis"
    );
}

#[test]
fn failing_offsets_give_the_radius_limit() {
    let mut dec = decider();
    // The dome is 3 high, so no ball wider than it rests on both faces.
    let dome = [plane("2", "-1"), sphere("0", "5")];
    let limit = RollingTorus::radius_limit(
        [&dome[0], &dome[1]],
        [Inside, Inside],
        &ring("2", "21"),
        &mut dec,
    )
    .unwrap()
    .unwrap();
    assert!(same(&limit, &real("3/2")));
    // A lens 4 thick on its axis: the inner offsets keep crossing while the
    // spheres' radii less the ball's stay more than half their centres apart.
    let lens = [sphere("3", "5"), sphere("-3", "5")];
    let limit = RollingTorus::radius_limit(
        [&lens[0], &lens[1]],
        [Inside, Inside],
        &ring("0", "16"),
        &mut dec,
    )
    .unwrap()
    .unwrap();
    assert!(same(&limit, &real("2")));
    // A seated ball outside both faces always finds a centre.
    let seat = [plane("20", "1"), sphere("23", "5")];
    assert!(
        RollingTorus::radius_limit(
            [&seat[0], &seat[1]],
            [Outside, Outside],
            &ring("20", "16"),
            &mut dec
        )
        .unwrap()
        .is_none()
    );
    // A cone's offset rounds its apex, so its line proves no bound.
    let frustum = [plane("6", "1"), cone("15", "1/3", true)];
    assert!(
        RollingTorus::radius_limit(
            [&frustum[0], &frustum[1]],
            [Inside, Inside],
            &ring("6", "9"),
            &mut dec
        )
        .unwrap()
        .is_none()
    );
}

fn circle_at(rho: &Real, z: &Real) -> Circle {
    Circle {
        center: Point3::new(real("0"), real("0"), z.clone()),
        normal: up(),
        radius_sq: rho.square(),
        ref_dir: vector(["1", "0", "0"]),
    }
}
// The circle a hand-worked half-plane point turns, compared as a point set.
fn on(circle: &Circle, rho: &str, z: &str) -> bool {
    same(&circle.center.z, &real(z))
        && same(&circle.center.x, &real("0"))
        && same(&circle.center.y, &real("0"))
        && same(&circle.radius_sq, &real(rho).square())
}

/// The plane z = 1 cuts the tube of the torus R = 6, r = 2 in a segment of
/// height 1, whose inscribed circle has radius 1/2 about (6, 3/2). Both rings,
/// at rho = 6 + sqrt 3 and 6 - sqrt 3, reach that one centre with contacts at
/// (6, 1) and (6, 2), and each keeps the half tube on its own side.
#[test]
fn a_ball_at_the_inscribed_radius_keeps_the_half_on_its_edge_side() {
    let supports = [plane("1", "-1"), torus("6", "2")];
    let supports = [&supports[0], &supports[1]];
    let root3 = root("3");
    for (rho, middle) in [(&real("6") + &root3, "13/2"), (&real("6") - &root3, "11/2")] {
        let edge = circle_at(&rho, &real("1"));
        let band = RollingTorus::inscribed(
            supports,
            [Inside, Inside],
            &edge,
            &real("1/2"),
            &mut decider(),
        )
        .unwrap()
        .expect("the inscribed half tube");
        assert!(same(&band.torus.major, &real("6")));
        assert!(same(&band.torus.minor, &real("1/2")));
        assert!(same(&band.torus.center.z, &real("3/2")));
        assert!(on(&band.contacts[0], "6", "1"), "plane contact");
        assert!(on(&band.contacts[1], "6", "2"), "torus contact");
        assert!(on(&band.middle, middle, "3/2"), "band side");
        // Away from the limit the regular strip answers, and past it nothing.
        for radius in ["2/5", "3/5"] {
            let away = RollingTorus::inscribed(
                supports,
                [Inside, Inside],
                &edge,
                &real(radius),
                &mut decider(),
            );
            assert!(away.unwrap().is_none(), "{radius}");
        }
        assert!(
            RollingTorus::new(
                supports,
                [Inside, Inside],
                &edge,
                &real("1/2"),
                &mut decider()
            )
            .is_err()
        );
    }
}

/// Two circle profiles: the sphere |q - (0, 8)| = 10 and the tube
/// |q - (6, 0)| = 2 bound a lens whose offsets 10 - r and 2 - r touch at
/// r = 1, about (27/5, 4/5), with feet (6, 0) and (24/5, 8/5). The rings lie
/// where 3 rho - 4 z = 17 meets the tube, at z = (4 +- 9 sqrt 11) / 25.
#[test]
fn two_circle_profiles_touch_at_the_inscribed_radius() {
    let supports = [sphere("8", "10"), torus("6", "2")];
    let supports = [&supports[0], &supports[1]];
    let root11 = root("11");
    for (z, middle) in [
        (
            &(&real("4") + &(&real("9") * &root11)) / &real("25"),
            ["31/5", "7/5"],
        ),
        (
            &(&real("4") - &(&real("9") * &root11)) / &real("25"),
            ["23/5", "1/5"],
        ),
    ] {
        let rho = &(&real("17") + &(&real("4") * &z)) / &real("3");
        let edge = circle_at(&rho, &z);
        let band = RollingTorus::inscribed(
            supports,
            [Inside, Inside],
            &edge,
            &real("1"),
            &mut decider(),
        )
        .unwrap()
        .expect("the inscribed half tube");
        assert!(same(&band.torus.major, &real("27/5")));
        assert!(same(&band.torus.center.z, &real("4/5")));
        assert!(on(&band.contacts[0], "6", "0"), "sphere contact");
        assert!(on(&band.contacts[1], "24/5", "8/5"), "torus contact");
        assert!(on(&band.middle, middle[0], middle[1]), "band side");
    }
}
