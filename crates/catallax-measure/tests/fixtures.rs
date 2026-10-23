//! These exercise finished pieces 1–3 and the independent oracle, not the stubs.
mod common;
use catallax_shapes::{Curve, Point3, PointDef, Surface, Vector3};
use common::oracle::{to_kernel, to_q};
use common::*;
use rug::Rational as Q;

fn xyz(p: &Point3) -> [Q; 3] {
    p.rational_components()
        .expect("rational measurement fixture")
        .map(to_q)
}
fn vector(p: &Vector3) -> [Q; 3] {
    p.rational_components()
        .expect("rational measurement fixture")
        .map(to_q)
}
fn sub(a: &[Q; 3], b: &[Q; 3]) -> [Q; 3] {
    std::array::from_fn(|i| a[i].clone() - &b[i])
}
fn dot(a: &[Q; 3], b: &[Q; 3]) -> Q {
    (0..3).map(|i| a[i].clone() * &b[i]).sum()
}
fn cross(a: &[Q; 3], b: &[Q; 3]) -> [Q; 3] {
    std::array::from_fn(|i| {
        a[(i + 1) % 3].clone() * &b[(i + 2) % 3] - a[(i + 2) % 3].clone() * &b[(i + 1) % 3]
    })
}
fn point_on_curve(p: &Point3, c: &Curve) {
    match c {
        Curve::Line(c) => {
            let d = sub(&xyz(p), &xyz(&c.origin));
            assert_eq!(cross(&d, &vector(&c.dir)), [Q::new(), Q::new(), Q::new()]);
        }
        Curve::Circle(c) => {
            let d = sub(&xyz(p), &xyz(&c.center));
            assert_eq!(dot(&d, &vector(&c.normal)), 0);
            assert_eq!(dot(&d, &d), to_q(c.radius_sq.rational_ref().unwrap()));
        }
        Curve::Ellipse(c) => {
            let d = sub(&xyz(p), &xyz(&c.center));
            let a = vector(&c.major_dir);
            let b = cross(&vector(&c.normal), &a);
            assert_eq!(dot(&d, &vector(&c.normal)), 0);
            let da = dot(&d, &a);
            let db = dot(&d, &b);
            assert_eq!(
                da.clone() * da / (dot(&a, &a) * to_q(c.major_sq.rational_ref().unwrap()))
                    + db.clone() * db / (dot(&b, &b) * to_q(c.minor_sq.rational_ref().unwrap())),
                1
            );
        }
        _ => panic!("analytic fixture expected"),
    }
}
fn point_on_surface(p: &Point3, s: &Surface) {
    let point = xyz(p);
    match s {
        Surface::Plane(s) => assert_eq!(dot(&sub(&point, &xyz(&s.origin)), &vector(&s.normal)), 0),
        Surface::Sphere(s) => {
            let d = sub(&point, &xyz(&s.center));
            let radius = to_q(s.radius.rational_ref().unwrap());
            assert_eq!(dot(&d, &d), radius.clone() * radius);
        }
        Surface::Cylinder(s) => {
            let d = sub(&point, &xyz(&s.origin));
            let axis = vector(&s.axis);
            let c = cross(&d, &axis);
            let radius = to_q(s.radius().rational_ref().unwrap());
            assert_eq!(dot(&c, &c), radius.clone() * radius * dot(&axis, &axis));
        }
        Surface::Cone(s) => {
            let d = sub(&point, &xyz(&s.apex));
            let axis = vector(&s.axis);
            let c = cross(&d, &axis);
            let da = dot(&d, &axis);
            let slope = to_q(s.slope.rational_ref().unwrap());
            assert!(da >= 0);
            assert_eq!(dot(&c, &c), slope.clone() * slope * da.clone() * da);
        }
        _ => panic!("loopless torus has no vertices"),
    }
}

#[test]
fn fixtures_have_valid_incidence_and_intended_topology_before_measurement_exists() {
    let counts = [
        (6, 12, 8),
        (6, 12, 8),
        (3, 2, 0),
        (1, 0, 0),
        (3, 2, 0),
        (2, 1, 1),
        (1, 0, 0),
        (0, 0, 0),
        (7, 14, 8),
        (12, 24, 16),
        (2, 1, 0),
        (2, 1, 0),
        (4, 6, 4),
        (4, 6, 4),
        (4, 6, 4),
        (4, 4, 2),
    ];
    for (case, expected) in CASES.into_iter().zip(counts) {
        for (scale, pose) in [
            (r(1), Pose::identity()),
            (ten(-30), poses().remove(0)),
            (ten(30), poses().remove(19)),
        ] {
            let f = build(case, &scale, &pose);
            let m = &f.model;
            let c = m.counts(f.body);
            assert_eq!((c.faces, c.edges, c.vertices), expected, "{case:?}");
            for edge in m.edges_of_body(f.body) {
                let e = m.edge(edge);
                for id in [e.start, e.end].into_iter().flatten() {
                    let PointDef::Explicit(point) = m.geom().point(m.vertex(id).point) else {
                        panic!("explicit fixture expected")
                    };
                    point_on_curve(point, m.geom().curve(e.curve));
                    for coedge in m.coedges_of_edge(edge) {
                        let (_, face) = m.owner_of_coedge(coedge).unwrap();
                        point_on_surface(point, m.geom().surface(m.face(face).surface));
                    }
                }
            }
            // Independently integrate each face's closed form before testing the
            // kernel sum. This catches mismatched cap/side/void oracle assignments.
            let sum = f
                .faces
                .iter()
                .fold(rat(0, 1), |s, (_, a)| s + a.clone())
                .bounds(512);
            let whole = f.expected.area.bounds(512);
            assert!(
                sum.lo <= whole.hi && whole.lo <= sum.hi,
                "{case:?}: face formulas do not sum to body area"
            );
        }
    }
}

#[test]
fn oracle_directed_rounding_handles_cancellation_division_and_negative_values() {
    // Squaring directed sqrt bounds must surround 5; subtraction and reciprocal
    // must reverse the correct endpoints. No catallax_numbers::Ball is an oracle.
    let b = sqrt(5).bounds(512);
    assert!(b.lo.clone() * &b.lo < 5 && b.hi.clone() * &b.hi > 5);
    for expr in [
        rat(5, 1) - sqrt(5) * sqrt(5),
        pi() - pi(),
        rat(1, 1) - pi() / pi(),
    ] {
        let b = expr.bounds(512);
        assert!(b.lo <= 0 && b.hi >= 0);
    }
    for expr in [rat(-8, 3) / pi(), rat(64, 1) - rat(4, 1) * pi()] {
        let coarse = expr.bounds(128);
        let fine = expr.bounds(512);
        assert!(coarse.lo <= fine.lo && fine.hi <= coarse.hi);
        let enclosure = Enclosure {
            lo: to_kernel(&coarse.lo),
            hi: to_kernel(&coarse.hi),
        };
        contains(&enclosure, &expr);
    }
}
