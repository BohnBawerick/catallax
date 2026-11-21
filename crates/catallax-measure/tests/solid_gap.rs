//! Minimum distance between two whole solids: face pairs far from the minimum
//! need no closed form, a line parallel to a circle's plane has one, turned
//! parallel faces answer from their edges, and pairs no closed form covers
//! settle by proven stationary pairs.
use catallax_measure::{Accuracy, solid_gap};
use catallax_numbers::{Decider, Rational, set_threads, with_operation_budget};
use catallax_shapes::{Isometry, Pose, RPoint3, RVector3, Rotation, Vector3};
use catallax_topology::{BodyId, Model, block, cone, cylinder, torus, transformed_body};
use std::num::NonZeroUsize;

fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn at(x: &str, y: &str, z: &str) -> Pose {
    Pose::from_rational(Rotation::identity(), RVector3::new(q(x), q(y), q(z)))
}
/// A turn by `degrees` about the axis through the origin, then a shift.
fn turned(
    model: &mut Model,
    body: BodyId,
    axis: [&str; 3],
    degrees: &str,
    to: [&str; 3],
) -> BodyId {
    let mut map = Isometry::around_axis(
        &RPoint3::new(q("0"), q("0"), q("0")),
        &RVector3::new(q(axis[0]), q(axis[1]), q(axis[2])),
        &q(degrees),
        &mut Decider::default(),
    )
    .unwrap();
    map.translation = map.translation.plus(&Vector3::new(
        q(to[0]).into(),
        q(to[1]).into(),
        q(to[2]).into(),
    ));
    transformed_body(model, body, &map, Clone::clone)
}
fn gap(
    model: &Model,
    a: BodyId,
    b: BodyId,
) -> catallax_numbers::Answer<catallax_measure::Enclosure> {
    solid_gap(model, a, b, &Accuracy::default(), &mut Decider::default()).unwrap()
}
fn within_accuracy(e: &catallax_measure::Enclosure) -> bool {
    e.width() <= &Accuracy::default().relative * &e.lo
}

#[test]
fn an_edge_parallel_to_a_rim_plane_meets_the_rim_in_closed_form() {
    // The block's lower front edge runs along x at y = 3, z = 2. The nearest
    // cylinder point is on its top rim at (0, 1, 1): the gap is sqrt(5).
    let mut model = Model::new();
    let pin = cylinder(&mut model, q("1"), q("1"), &at("0", "0", "0")).unwrap();
    let bar = block(&mut model, [q("10"), q("1"), q("1")], &at("-5", "3", "2")).unwrap();
    let answer = gap(&model, pin, bar);
    let e = answer.value;
    assert!(&e.lo * &e.lo <= q("5") && &e.hi * &e.hi >= q("5"), "{e:?}");
    assert!(within_accuracy(&e) && answer.marks.is_empty());
}

#[test]
fn faces_far_from_the_minimum_need_no_closed_form() {
    // Cone rims against block edges have no closed form, but every such pair
    // is farther than the base plane's 188 mm and is never visited.
    let mut model = Model::new();
    let base = block(&mut model, [q("10"), q("10"), q("10")], &at("0", "0", "0")).unwrap();
    let tip = cone(&mut model, q("2"), q("1"), q("4"), &at("200", "5", "3")).unwrap();
    let answer = gap(&model, base, tip);
    assert_eq!(answer.value.lo, q("188"));
    assert_eq!(answer.value.hi, q("188"));
    assert!(answer.marks.is_empty());
}

#[test]
fn turned_parallel_faces_answer_from_their_edges() {
    // A 4 mm cube turned 30 degrees about z: its top and bottom faces stay
    // parallel to the base block's, on planes with no rational frame. Its
    // nearest corner is 12 - sqrt(3) from the origin in x: the gap is 2 - sqrt(3).
    let mut model = Model::new();
    let base = block(&mut model, [q("10"), q("10"), q("10")], &at("0", "0", "0")).unwrap();
    let cube = block(&mut model, [q("4"), q("4"), q("4")], &at("-2", "-2", "-2")).unwrap();
    let cube = turned(&mut model, cube, ["0", "0", "1"], "30", ["13", "5", "5"]);
    let answer = gap(&model, base, cube);
    let e = answer.value;
    let (a, b) = (&q("2") - &e.lo, &q("2") - &e.hi);
    assert!(
        b > q("0") && &a * &a >= q("3") && &b * &b <= q("3"),
        "{e:?}"
    );
    assert!(within_accuracy(&e) && answer.marks.is_empty());
}

#[test]
fn rims_at_non_parallel_normals_settle_by_proven_stationary_pairs() {
    // The second cylinder lies along y. Its near rim faces the first
    // cylinder's side at sqrt(26) - 2, a line against a circle in a parallel
    // plane. The two cylinders' rims meet at non-parallel normals, which no
    // closed form covers; their stationary pairs are proven one by one, so the
    // stated gap is the closed form's value within the requested accuracy.
    let mut model = Model::new();
    let a = cylinder(&mut model, q("2"), q("4"), &at("0", "0", "-2")).unwrap();
    let b = cylinder(&mut model, q("2"), q("4"), &at("0", "0", "-2")).unwrap();
    let b = turned(&mut model, b, ["1", "0", "0"], "90", ["7", "3", "1"]);
    let answer = gap(&model, a, b);
    let e = answer.value;
    let (lo, hi) = (&e.lo + &q("2"), &e.hi + &q("2"));
    assert!(&lo * &lo <= q("26") && &hi * &hi >= q("26"), "{e:?}");
    assert!(
        within_accuracy(&e) && answer.marks.is_empty(),
        "{e:?} {:?}",
        answer.marks
    );
}

/// A face as a map of the unit square, in the solid's local frame.
type Face = Box<dyn Fn(f64, f64) -> [f64; 3]>;

fn round(r0: f64, r1: f64, z0: f64, z1: f64) -> Vec<Face> {
    let tau = std::f64::consts::TAU;
    let disc = move |r: f64, z: f64| -> Face {
        Box::new(move |u, v| [r * u * (tau * v).cos(), r * u * (tau * v).sin(), z])
    };
    vec![
        Box::new(move |u, v| {
            let r = r0 + (r1 - r0) * v;
            [r * (tau * u).cos(), r * (tau * u).sin(), z0 + (z1 - z0) * v]
        }),
        disc(r0, z0),
        disc(r1, z1),
    ]
}

fn ring(major: f64, minor: f64) -> Vec<Face> {
    let tau = std::f64::consts::TAU;
    vec![Box::new(move |u, v| {
        let r = major + minor * (tau * v).cos();
        [
            r * (tau * u).cos(),
            r * (tau * u).sin(),
            minor * (tau * v).sin(),
        ]
    })]
}

/// The faces moved as [`turned`] moves the solid: a turn by `degrees` about the
/// unit coordinate `axis` through the origin, then a shift.
fn placed(faces: Vec<Face>, axis: usize, degrees: f64, to: [f64; 3]) -> Vec<Face> {
    let (s, c) = degrees.to_radians().sin_cos();
    faces
        .into_iter()
        .map(|f| -> Face {
            Box::new(move |u, v| {
                let p = f(u, v);
                let (i, j) = ((axis + 1) % 3, (axis + 2) % 3);
                let mut out = p;
                out[i] = c * p[i] - s * p[j];
                out[j] = s * p[i] + c * p[j];
                [out[0] + to[0], out[1] + to[1], out[2] + to[2]]
            })
        })
        .collect()
}

fn apart(a: &Face, b: &Face, x: [f64; 4]) -> f64 {
    let (p, q) = (a(x[0], x[1]), b(x[2], x[3]));
    ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt()
}

/// An independent floating-point yardstick: the best pair of a 48 x 48 grid on
/// each face pair, refined by grids that shrink around the best pair found.
/// It is the distance of a real pair of boundary points, so the certified
/// lower end may not exceed it.
fn sampled_minimum(a: &[Face], b: &[Face]) -> f64 {
    const N: usize = 48;
    let grid = |f: &Face| -> Vec<([f64; 2], [f64; 3])> {
        let mut out = Vec::new();
        for i in 0..=N {
            for j in 0..=N {
                let (u, v) = (i as f64 / N as f64, j as f64 / N as f64);
                out.push(([u, v], f(u, v)));
            }
        }
        out
    };
    let mut best = f64::INFINITY;
    for fa in a {
        let ga = grid(fa);
        for fb in b {
            let gb = grid(fb);
            let mut start = ([0.0; 4], f64::INFINITY);
            for (x, p) in &ga {
                for (y, q) in &gb {
                    let d = (p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2);
                    if d < start.1 {
                        start = ([x[0], x[1], y[0], y[1]], d);
                    }
                }
            }
            let (mut x, mut h) = (start.0, 1.0 / N as f64);
            let mut d = apart(fa, fb, x);
            for _ in 0..48 {
                let centre = x;
                for k in 0..7usize.pow(4) {
                    let mut y = centre;
                    for (axis, value) in y.iter_mut().enumerate() {
                        let step = (k / 7usize.pow(axis as u32)) % 7;
                        *value = (*value + h * (step as f64 - 3.0) / 3.0).clamp(0.0, 1.0);
                    }
                    let e = apart(fa, fb, y);
                    if e < d {
                        (d, x) = (e, y);
                    }
                }
                h /= 2.0;
            }
            best = best.min(d);
        }
    }
    best
}

#[test]
fn cone_and_torus_pairs_settle_by_proven_stationary_pairs() {
    // Pairs with a cone, or two tori at non-parallel axes, have no closed form
    // at their minimum: a cylinder's side or rim against a cone, two cones, a
    // torus against a cone and two tori. Each answer is proven to the requested
    // accuracy, and agrees with an independent floating-point search.
    let cylinder_faces = || round(2.0, 2.0, -2.0, 2.0);
    let cone_faces = || round(2.0, 1.0, -2.0, 2.0);
    let torus_faces = || ring(1.5, 0.5);
    type Solid = (fn(&mut Model) -> BodyId, fn() -> Vec<Face>);
    let cyl: Solid = (
        |m| cylinder(m, q("2"), q("4"), &at("0", "0", "-2")).unwrap(),
        cylinder_faces,
    );
    let con: Solid = (
        |m| cone(m, q("2"), q("1"), q("4"), &at("0", "0", "-2")).unwrap(),
        cone_faces,
    );
    let tor: Solid = (
        |m| torus(m, q("3/2"), q("1/2"), &at("0", "0", "0")).unwrap(),
        torus_faces,
    );
    // The second solid turned 30 degrees about z, or tipped 90 about x, then
    // moved to (7, 3, 1).
    let cases: [(Solid, Solid, usize, &str); 4] = [
        (cyl, con, 2, "30"),
        (con, con, 0, "90"),
        (con, tor, 2, "30"),
        (tor, tor, 0, "90"),
    ];
    for ((make_a, faces_a), (make_b, faces_b), axis, degrees) in cases {
        let mut model = Model::new();
        let a = make_a(&mut model);
        let b = make_b(&mut model);
        let mut turn = ["0"; 3];
        turn[axis] = "1";
        let b = turned(&mut model, b, turn, degrees, ["7", "3", "1"]);
        let answer = gap(&model, a, b);
        let e = answer.value;
        let reference = sampled_minimum(
            &faces_a(),
            &placed(faces_b(), axis, degrees.parse().unwrap(), [7.0, 3.0, 1.0]),
        );
        let (lo, hi) = (e.lo.to_f64_nearest(), e.hi.to_f64_nearest());
        assert!(
            within_accuracy(&e) && answer.marks.is_empty(),
            "{degrees} about {axis}: {e:?} {:?}",
            answer.marks
        );
        assert!(
            lo <= reference + 1e-12 && (hi - reference).abs() <= 1e-9,
            "{degrees} about {axis}: [{lo}, {hi}] against {reference}"
        );
    }
}

#[test]
fn proven_stationary_pairs_answer_the_same_at_every_thread_count() {
    let mut model = Model::new();
    let a = cylinder(&mut model, q("2"), q("4"), &at("0", "0", "-2")).unwrap();
    let b = cone(&mut model, q("2"), q("1"), q("4"), &at("0", "0", "-2")).unwrap();
    let b = turned(&mut model, b, ["0", "0", "1"], "30", ["7", "3", "1"]);
    let records: Vec<_> = [1, 2, 4]
        .map(|n| {
            set_threads(NonZeroUsize::new(n));
            let (answer, work) = with_operation_budget(100_000_000, || gap(&model, a, b));
            (answer.value, answer.marks.len(), work.used)
        })
        .into();
    set_threads(None);
    assert!(records.windows(2).all(|w| w[0] == w[1]), "{records:?}");
}
