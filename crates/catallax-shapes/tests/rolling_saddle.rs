//! The fitted ball surface between a narrow cylinder and a wide one it passes
//! through at an angle. The checks below rebuild the ball's path with plain
//! bisection on the two distance conditions, so they share nothing with the law.

use catallax_numbers::{Ball, BudgetPolicy, Decider, Outcome, Rational, Real, Refusal};
use catallax_shapes::{
    BSplineSurface, BallSide, Cylinder, Point3, Pose, RollingSaddle, RollingSaddleDefinition,
    Rotation, SplineBoundary, Surface, Vector3,
};

fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 1_000_000_000,
        ..Default::default()
    })
}
fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}
fn real(value: &str) -> Real {
    Real::from(q(value))
}
fn vector(v: [&str; 3]) -> Vector3 {
    let [x, y, z] = v.map(real);
    Vector3::new(x, y, z)
}
fn point(v: [&str; 3]) -> Point3 {
    let [x, y, z] = v.map(real);
    Point3::new(x, y, z)
}
fn cylinder(origin: [&str; 3], axis: [&str; 3], radius: &str) -> Cylinder {
    let axis = vector(axis);
    let ref_dir = Vector3::new(axis.y.clone(), -&axis.x, Real::zero());
    let ref_dir = if axis.x.as_rational().is_some_and(|x| x.is_zero())
        && axis.y.as_rational().is_some_and(|y| y.is_zero())
    {
        vector(["1", "0", "0"])
    } else {
        ref_dir
    };
    Cylinder {
        origin: point(origin),
        axis,
        radius_sq: real(radius).square(),
        ref_dir,
    }
}

type V = [f64; 3];
fn sub(a: V, c: V) -> V {
    [a[0] - c[0], a[1] - c[1], a[2] - c[2]]
}
fn add(a: V, c: V, scale: f64) -> V {
    [
        a[0] + c[0] * scale,
        a[1] + c[1] * scale,
        a[2] + c[2] * scale,
    ]
}
fn dot(a: V, c: V) -> f64 {
    a[0] * c[0] + a[1] * c[1] + a[2] * c[2]
}
fn cross(a: V, c: V) -> V {
    [
        a[1] * c[2] - a[2] * c[1],
        a[2] * c[0] - a[0] * c[2],
        a[0] * c[1] - a[1] * c[0],
    ]
}
fn norm(a: V) -> f64 {
    dot(a, a).sqrt()
}
fn unit(a: V) -> V {
    let scale = 1.0 / norm(a);
    [a[0] * scale, a[1] * scale, a[2] * scale]
}
fn n(value: &Real) -> f64 {
    decider()
        .enclose_at(value, 128)
        .unwrap()
        .mid()
        .to_f64_nearest()
}
fn f(v: &Vector3) -> V {
    [&v.x, &v.y, &v.z].map(n)
}
fn g(v: &Point3) -> V {
    [&v.x, &v.y, &v.z].map(n)
}

// A cylinder as plain numbers.
#[derive(Clone, Copy)]
struct Tube {
    at: V,
    axis: V,
    radius: f64,
}
impl Tube {
    fn of(c: &Cylinder) -> Self {
        Self {
            at: g(&c.origin),
            axis: unit(f(&c.axis)),
            radius: n(&c.radius_sq).sqrt(),
        }
    }
    fn distance(&self, p: V) -> f64 {
        let d = sub(p, self.at);
        norm(add(d, self.axis, -dot(d, self.axis)))
    }
}

// The ball's path found from the two distance conditions alone: walk the narrow
// offset's generator at angle `t` until it reaches the wide offset's radius.
struct Path {
    wide: Tube,
    narrow: Tube,
    across: V,
    side: V,
    big: f64,
    small: f64,
    leaving: bool,
    radius: f64,
}
impl Path {
    fn new(wide: Tube, narrow: Tube, offsets: [f64; 2], leaving: bool, radius: f64) -> Self {
        let across = unit(cross(narrow.axis, wide.axis));
        let side = cross(narrow.axis, across);
        Self {
            wide,
            narrow,
            across,
            side,
            big: offsets[0],
            small: offsets[1],
            leaving,
            radius,
        }
    }
    fn centre(&self, t: f64) -> V {
        let ray = add(
            [
                self.across[0] * t.cos(),
                self.across[1] * t.cos(),
                self.across[2] * t.cos(),
            ],
            self.side,
            t.sin(),
        );
        let start = add(self.narrow.at, ray, self.small);
        let direction = if self.leaving { 1.0 } else { -1.0 };
        // The wide axis's closest station on this generator, then outward.
        let d = sub(start, self.wide.at);
        let flat = add(
            self.narrow.axis,
            self.wide.axis,
            -dot(self.narrow.axis, self.wide.axis),
        );
        let near = -dot(d, flat) / dot(flat, flat);
        let (mut low, mut high) = (0.0, 1000.0);
        for _ in 0..120 {
            let middle = 0.5 * (low + high);
            let at = add(start, self.narrow.axis, near + direction * middle);
            if self.wide.distance(at) < self.big {
                low = middle;
            } else {
                high = middle;
            }
        }
        add(start, self.narrow.axis, near + direction * low)
    }
    fn foot(&self, c: V) -> [V; 2] {
        let wide = {
            let d = sub(c, self.wide.at);
            add(
                sub(c, [0.0; 3]),
                unit(add(d, self.wide.axis, -dot(d, self.wide.axis))),
                self.wide.radius - self.big,
            )
        };
        let narrow = {
            let d = sub(c, self.narrow.at);
            add(
                c,
                unit(add(d, self.narrow.axis, -dot(d, self.narrow.axis))),
                self.narrow.radius - self.small,
            )
        };
        [wide, narrow]
    }
}

struct Case {
    name: &'static str,
    wide: Cylinder,
    narrow: Cylinder,
    radius: &'static str,
    sides: [BallSide; 2],
    leaving: bool,
}

fn cases() -> Vec<Case> {
    use BallSide::{Inside, Outside};
    let tilted = || cylinder(["1", "2", "3"], ["0", "3", "4"], "4");
    let wide = || cylinder(["0", "0", "0"], ["0", "0", "1"], "10");
    let square = || cylinder(["3/2", "0", "0"], ["0", "1", "0"], "4");
    vec![
        Case {
            name: "boss",
            wide: wide(),
            narrow: tilted(),
            radius: "1",
            sides: [Outside, Outside],
            leaving: true,
        },
        Case {
            name: "far boss",
            wide: wide(),
            narrow: tilted(),
            radius: "1",
            sides: [Outside, Outside],
            leaving: false,
        },
        Case {
            name: "bored rim",
            wide: wide(),
            narrow: tilted(),
            radius: "3/2",
            sides: [Inside, Outside],
            leaving: true,
        },
        Case {
            name: "inner rim",
            wide: wide(),
            narrow: tilted(),
            radius: "1",
            sides: [Outside, Inside],
            leaving: true,
        },
        Case {
            name: "cavity",
            wide: wide(),
            narrow: tilted(),
            radius: "1",
            sides: [Inside, Inside],
            leaving: false,
        },
        Case {
            name: "square tee off centre",
            wide: wide(),
            narrow: square(),
            radius: "1/2",
            sides: [Outside, Outside],
            leaving: true,
        },
    ]
}

fn saddle(case: &Case) -> RollingSaddle {
    let mut dec = decider();
    let definition = RollingSaddleDefinition::new(
        &case.wide,
        &case.narrow,
        &real(case.radius),
        case.sides,
        case.leaving,
        &mut dec,
    )
    .unwrap();
    RollingSaddle::fit(definition, &mut dec).unwrap()
}

fn zero(value: &Real, dec: &mut Decider) -> bool {
    matches!(dec.sign(value), Outcome::ProvenCoincident(_))
}

fn on_cylinder(point: &Point3, cylinder: &Cylinder, radius_sq: &Real, dec: &mut Decider) -> bool {
    let d = point.minus(&cylinder.origin);
    let flat = d.cross(&cylinder.axis).norm_squared();
    zero(&(&flat - &(radius_sq * &cylinder.axis.norm_squared())), dec)
}

fn offset_sq(case: &Case, which: usize) -> Real {
    let (cylinder, side) = if which == 0 {
        (&case.wide, case.sides[0])
    } else {
        (&case.narrow, case.sides[1])
    };
    let step = match side {
        BallSide::Outside => real(case.radius),
        BallSide::Inside => -&real(case.radius),
    };
    (&cylinder.radius_sq.sqrt() + &step).square()
}

fn chart_samples() -> Vec<Rational> {
    (0..=40)
        .map(|i| &q("-1") + &(q("1/10") * Rational::from_i64(i)))
        .collect()
}

#[test]
fn the_law_rests_the_ball_on_both_cylinders_exactly() {
    for case in cases() {
        let mut dec = decider();
        let definition = RollingSaddleDefinition::new(
            &case.wide,
            &case.narrow,
            &real(case.radius),
            case.sides,
            case.leaving,
            &mut dec,
        )
        .unwrap();
        let (big, small) = (offset_sq(&case, 0), offset_sq(&case, 1));
        for u in chart_samples().iter().step_by(5) {
            let [centre, wide, narrow] = definition.section(u).unwrap();
            assert!(
                on_cylinder(&centre, &case.wide, &big, &mut dec),
                "{}",
                case.name
            );
            assert!(
                on_cylinder(&centre, &case.narrow, &small, &mut dec),
                "{}",
                case.name
            );
            assert!(
                on_cylinder(&wide, &case.wide, &case.wide.radius_sq, &mut dec),
                "{}",
                case.name
            );
            assert!(
                on_cylinder(&narrow, &case.narrow, &case.narrow.radius_sq, &mut dec),
                "{}",
                case.name
            );
            let ball = real(case.radius).square();
            assert!(zero(
                &(&wide.minus(&centre).norm_squared() - &ball),
                &mut dec
            ));
            assert!(zero(
                &(&narrow.minus(&centre).norm_squared() - &ball),
                &mut dec
            ));
        }
    }
}

#[test]
fn the_rails_lie_exactly_on_their_cylinders_and_close() {
    for case in cases() {
        let fitted = saddle(&case);
        let mut dec = decider();
        let [wide, narrow] = fitted.rails();
        assert!(
            wide.is_boundary_of(fitted.carrier(), SplineBoundary::VMin),
            "{}",
            case.name
        );
        assert!(
            narrow.is_boundary_of(fitted.carrier(), SplineBoundary::VMax),
            "{}",
            case.name
        );
        for (rail, cylinder) in [(&wide, &case.wide), (&narrow, &case.narrow)] {
            assert_eq!(rail.is_closed(), Some(true), "{}", case.name);
            for u in chart_samples().iter().step_by(3) {
                let at = rail.point_at(&Real::from(u.clone()), &mut dec).unwrap();
                assert!(
                    on_cylinder(&at, cylinder, &cylinder.radius_sq, &mut dec),
                    "{} at {u:?}",
                    case.name
                );
            }
        }
        // Equal end columns are what a reader that closes the surface needs.
        let surface = fitted.carrier();
        let [columns, rows] = surface.control_count();
        for row in 0..rows {
            let first = row;
            let last = (columns - 1) * rows + row;
            assert!(zero(
                &(&surface.weights()[first] - &surface.weights()[last]),
                &mut dec
            ));
            let gap = surface.controls()[first].minus(&surface.controls()[last]);
            assert!(zero(&gap.norm_squared(), &mut dec), "{}", case.name);
        }
    }
}

fn sample(surface: &BSplineSurface, u: f64, v: f64, dec: &Decider) -> V {
    let ball = |x: f64| Ball::from_rational(&Rational::from_f64_exact(x).unwrap(), 128);
    surface
        .evaluate(&ball(u), &ball(v), dec)
        .unwrap()
        .map(|x| x.mid().to_f64_nearest())
}

// Each carrier point sits at the ball's radius from the path's centre at the
// same angle, on the short arc between the contacts, within the claimed bound.
#[test]
fn the_carrier_stays_within_its_bound_of_the_rolling_ball() {
    for case in cases() {
        let fitted = saddle(&case);
        let dec = decider();
        let bound = fitted.bound().upper.to_f64_nearest();
        assert!(bound <= 1e-4, "{}: {bound}", case.name);
        let radius = n(&real(case.radius));
        let (wide, narrow) = (Tube::of(&case.wide), Tube::of(&case.narrow));
        let offsets = [0, 1].map(|which| n(&offset_sq(&case, which)).sqrt());
        let path = Path::new(wide, narrow, offsets, case.leaving, radius);
        let mut worst: f64 = 0.0;
        for i in 0..240 {
            let u = -1.0 + 4.0 * (f64::from(i) + 0.37) / 240.0;
            let (chart, s) = if u <= 1.0 { (0.0, u) } else { (1.0, u - 2.0) };
            let t = chart * std::f64::consts::PI + 2.0 * s.atan();
            let centre = path.centre(t);
            let [first, second] = path.foot(centre);
            let (a, c) = (unit(sub(first, centre)), unit(sub(second, centre)));
            let normal = unit(cross(a, c));
            for j in 0..=8 {
                let v = f64::from(j) / 8.0;
                let p = sample(fitted.carrier(), u, v, &dec);
                let d = sub(p, centre);
                worst = worst.max((norm(d) - path.radius).abs());
                worst = worst.max(dot(d, normal).abs());
                // Inside the cone of the contact rays, by more than the bound.
                assert!(dot(cross(a, d), normal) > -2.0 * bound, "{}", case.name);
                assert!(dot(cross(d, c), normal) > -2.0 * bound, "{}", case.name);
            }
        }
        assert!(worst <= bound, "{}: {worst} over {bound}", case.name);
        eprintln!(
            "{}: {} spans, bound {bound:.3e}, worst seen {worst:.3e}",
            case.name,
            fitted.span_count()
        );
    }
}

#[test]
fn a_moved_copy_keeps_its_rails_on_the_moved_cylinders() {
    let case = &cases()[0];
    let fitted = saddle(case);
    let mut dec = decider();
    let rotation = Rotation::from_quaternion(q("1"), q("2"), q("3"), q("4")).unwrap();
    let pose = Pose::new(rotation, vector(["7/3", "-11/7", "13/5"]), &mut dec).unwrap();
    let moved = fitted.transformed(&pose);
    assert_eq!(moved.span_count(), fitted.span_count());
    let Surface::Cylinder(wide) = Surface::Cylinder(case.wide.clone()).transformed(&pose) else {
        unreachable!()
    };
    let at = &moved.rails()[0];
    let big = offset_sq(case, 0);
    for u in ["-1", "1/3", "2"] {
        let point = at.point_at(&real(u), &mut dec).unwrap();
        assert!(on_cylinder(&point, &wide, &wide.radius_sq, &mut dec), "{u}");
        let [centre, ..] = moved.definition().section(&q(u)).unwrap();
        assert!(on_cylinder(&centre, &wide, &big, &mut dec), "{u}");
    }
}

fn refused(wide: Cylinder, narrow: Cylinder, radius: &str, sides: [BallSide; 2]) -> Refusal {
    let mut dec = decider();
    RollingSaddleDefinition::new(&wide, &narrow, &real(radius), sides, true, &mut dec).unwrap_err()
}

#[test]
fn cases_outside_the_regular_ring_stay_typed_refusals() {
    use BallSide::{Inside, Outside};
    let wide = || cylinder(["0", "0", "0"], ["0", "0", "1"], "10");
    let cases = [
        // Parallel axes belong to the exact cylinder case.
        refused(
            wide(),
            cylinder(["1", "0", "0"], ["0", "0", "2"], "4"),
            "1",
            [Outside, Outside],
        ),
        // Equal radii meet at the singular crossing.
        refused(
            wide(),
            cylinder(["0", "0", "0"], ["0", "1", "0"], "10"),
            "1",
            [Outside, Outside],
        ),
        // The narrow offset reaches the wide offset's side, so one ring is gone.
        refused(
            wide(),
            cylinder(["7", "0", "0"], ["0", "1", "0"], "4"),
            "1",
            [Outside, Outside],
        ),
        // The ball does not fit in the bore.
        refused(
            wide(),
            cylinder(["0", "0", "0"], ["0", "1", "0"], "4"),
            "5",
            [Outside, Inside],
        ),
    ];
    for refusal in &cases {
        assert_eq!(refusal.code, "shapes.rolling-saddle", "{}", refusal.detail);
    }
}
