//! The sphere sector that closes a junction of three contacts of one ball,
//! convex or concave, independent of the strips that reach it.

use catallax_features::fillet_corner::SphereSector;
use catallax_numbers::{BudgetPolicy, Decider, Rational, Real};
use catallax_shapes::{Curve, Point3, RPoint3, Surface, Vector3};

fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 1_000_000,
        ..Default::default()
    })
}
fn r(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn vector(v: [&str; 3]) -> Vector3 {
    let [x, y, z] = v.map(|c| Real::from(r(c)));
    Vector3::new(x, y, z)
}
fn point(p: [&str; 3]) -> Point3 {
    let [x, y, z] = p.map(r);
    RPoint3::new(x, y, z).to_real()
}
fn exact(p: &Point3) -> [Rational; 3] {
    let p = p.as_rational().expect("rational point");
    [p.x, p.y, p.z]
}
fn axes(sign: &str) -> [Vector3; 3] {
    [
        vector([sign, "0", "0"]),
        vector(["0", sign, "0"]),
        vector(["0", "0", sign]),
    ]
}

#[test]
fn a_convex_sector_spans_the_cone_of_its_contact_directions() {
    let mut dec = decider();
    let sector = SphereSector::new(
        point(["1", "2", "3"]),
        r("2").into(),
        axes("1"),
        true,
        &mut dec,
    )
    .unwrap();
    assert!(sector.convex);
    let Surface::Sphere(sphere) = &sector.surface else {
        panic!("a sector lies on a sphere");
    };
    assert_eq!(exact(&sphere.center), [r("1"), r("2"), r("3")]);
    assert_eq!(sphere.radius.as_rational(), Some(r("2")));
    let feet: Vec<_> = sector.feet.iter().map(exact).collect();
    assert_eq!(
        feet,
        [
            [r("3"), r("2"), r("3")],
            [r("1"), r("4"), r("3")],
            [r("1"), r("2"), r("5")],
        ]
    );
    // (1, 2, 2) / 3 has every coefficient positive; flipping its first
    // component leaves the cone.
    let inside = point(["5/3", "10/3", "13/3"]);
    let outside = point(["1/3", "10/3", "13/3"]);
    assert!(sector.in_sector(&inside, &mut dec).unwrap());
    assert!(!sector.in_sector(&outside, &mut dec).unwrap());
    // Each boundary is the great circle through two contacts, turning from
    // the first to the second, and its points are on the closed sector.
    let Curve::Circle(edge) = sector.boundary(0, 1) else {
        panic!("a sector boundary is a circle");
    };
    assert_eq!(exact(&edge.center), [r("1"), r("2"), r("3")]);
    assert_eq!(edge.radius_sq.as_rational(), Some(r("4")));
    assert_eq!(
        edge.normal.as_rational().map(|v| [v.x, v.y, v.z]),
        Some([r("0"), r("0"), r("1")])
    );
    assert_eq!(
        edge.ref_dir.as_rational().map(|v| [v.x, v.y, v.z]),
        Some([r("1"), r("0"), r("0")])
    );
    let on_boundary = point(["11/5", "18/5", "3"]);
    assert!(sector.in_sector(&on_boundary, &mut dec).unwrap());
}

#[test]
fn a_concave_sector_is_the_convex_sector_of_the_complement() {
    let mut dec = decider();
    let sector = SphereSector::new(
        point(["6", "6", "5"]),
        r("1").into(),
        axes("-1"),
        false,
        &mut dec,
    )
    .unwrap();
    assert!(!sector.convex);
    let feet: Vec<_> = sector.feet.iter().map(exact).collect();
    assert_eq!(
        feet,
        [
            [r("5"), r("6"), r("5")],
            [r("6"), r("5"), r("5")],
            [r("6"), r("6"), r("4")],
        ]
    );
    // Toward the pocket corner at (5, 5, 4), not away from it.
    let toward = point(["17/3", "16/3", "13/3"]);
    let away = point(["19/3", "20/3", "17/3"]);
    assert!(sector.in_sector(&toward, &mut dec).unwrap());
    assert!(!sector.in_sector(&away, &mut dec).unwrap());
    let Curve::Circle(edge) = sector.boundary(1, 2) else {
        panic!("a sector boundary is a circle");
    };
    // (0, -1, 0) x (0, 0, -1) = (1, 0, 0).
    assert_eq!(
        edge.normal.as_rational().map(|v| [v.x, v.y, v.z]),
        Some([r("1"), r("0"), r("0")])
    );
}

#[test]
fn a_skew_sector_tests_the_dual_basis_not_the_directions() {
    // The acute prism corner: centre (-3, 1, 1), both points one unit away.
    let mut dec = decider();
    let directions = [
        vector(["0", "-1", "0"]),
        vector(["3/5", "4/5", "0"]),
        vector(["0", "0", "-1"]),
    ];
    let sector = SphereSector::new(
        point(["-3", "1", "1"]),
        r("1").into(),
        directions,
        true,
        &mut dec,
    )
    .unwrap();
    // Its dot product with (0, -1, 0) is negative, yet it is in the cone.
    let inside = point(["-7/3", "1/3", "2/3"]);
    assert!(sector.in_sector(&inside, &mut dec).unwrap());
    let outside = point(["-11/3", "1/3", "2/3"]);
    assert!(!sector.in_sector(&outside, &mut dec).unwrap());
    let duals = sector.dual_basis();
    for (i, dual) in duals.iter().enumerate() {
        for (j, direction) in sector.directions.iter().enumerate() {
            let want = if i == j { r("1") } else { r("0") };
            assert_eq!(dual.dot(direction).as_rational(), Some(want));
        }
    }
}

#[test]
fn degenerate_sectors_are_refused_as_unsupported_corners() {
    let mut dec = decider();
    let flat = [
        vector(["1", "0", "0"]),
        vector(["0", "1", "0"]),
        vector(["3/5", "4/5", "0"]),
    ];
    let refusal =
        SphereSector::new(point(["0", "0", "0"]), r("1").into(), flat, true, &mut dec).unwrap_err();
    assert_eq!(refusal.code, "fillet.unsupported-corner");
    assert!(refusal.detail.contains("independent"), "{refusal:?}");
    let refusal = SphereSector::new(
        point(["0", "0", "0"]),
        r("0").into(),
        axes("1"),
        true,
        &mut dec,
    )
    .unwrap_err();
    assert_eq!(refusal.code, "fillet.unsupported-corner");
    assert!(refusal.detail.contains("positive radius"), "{refusal:?}");
}
