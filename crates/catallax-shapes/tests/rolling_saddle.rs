//! The fitted ball surface between a narrow cylinder and a wide one it passes
//! through at an angle. The checks below rebuild the ball's path with plain
//! bisection on the two distance conditions, so they share nothing with the law.

use std::sync::OnceLock;

use catallax_numbers::{Ball, BudgetPolicy, Decider, Outcome, Rational, Real, Refusal, Sign};
use catallax_shapes::{
    BSplineSurface, BallSide, Cylinder, Point3, Pose, RollingSaddle, RollingSaddleDefinition,
    Rotation, SplineBoundary, SplineSupport, Surface, Vector3,
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

fn fit(case: &Case) -> RollingSaddle {
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

// Each listed case's fit, shared so a run pays for it once.
fn fitted(index: usize) -> RollingSaddle {
    static FITS: OnceLock<Vec<OnceLock<RollingSaddle>>> = OnceLock::new();
    let all = FITS.get_or_init(|| cases().iter().map(|_| OnceLock::new()).collect());
    all[index].get_or_init(|| fit(&cases()[index])).clone()
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
    for (index, case) in cases().iter().enumerate() {
        let fitted = fitted(index);
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
    for (index, case) in cases().iter().enumerate() {
        let fitted = fitted(index);
        let dec = decider();
        let bound = fitted.bound().upper.to_f64_nearest();
        assert!(bound <= 1e-4, "{}: {bound}", case.name);
        let radius = n(&real(case.radius));
        let (wide, narrow) = (Tube::of(&case.wide), Tube::of(&case.narrow));
        let offsets = [0, 1].map(|which| n(&offset_sq(case, which)).sqrt());
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
    let fitted = fitted(0);
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

// A triangulated grid, for counting crossings with nothing from the line solver.
struct Soup {
    triangles: Vec<[V; 3]>,
}
impl Soup {
    fn of(fitted: &RollingSaddle) -> Self {
        let dec = decider();
        let grid = fitted.display_grid(&q("1/200"), 2_000_000, &dec).unwrap();
        assert!(grid.chord <= q("1/200"));
        let point = |column: usize, row: usize| -> V {
            let at = &grid.points[column * grid.rows + row];
            std::array::from_fn(|i| at[i].mid().to_f64_nearest())
        };
        let mut triangles = Vec::new();
        for column in 0..grid.columns - 1 {
            for row in 0..grid.rows - 1 {
                let [a, c, d, e] = [
                    point(column, row),
                    point(column + 1, row),
                    point(column, row + 1),
                    point(column + 1, row + 1),
                ];
                triangles.push([a, c, e]);
                triangles.push([a, e, d]);
            }
        }
        Self { triangles }
    }
    // Moller and Trumbore's test, in line parameter order.
    fn crossings(&self, origin: V, direction: V) -> Vec<f64> {
        let mut found = Vec::new();
        for [a, c, d] in &self.triangles {
            let (e1, e2) = (sub(*c, *a), sub(*d, *a));
            let p = cross(direction, e2);
            let det = dot(e1, p);
            if det.abs() < 1e-14 {
                continue;
            }
            let t = sub(origin, *a);
            let u = dot(t, p) / det;
            let q = cross(t, e1);
            let v = dot(direction, q) / det;
            if u >= 0.0 && v >= 0.0 && u + v <= 1.0 {
                found.push(dot(e2, q) / det);
            }
        }
        found.sort_by(f64::total_cmp);
        found
    }
}

fn rational_of(point: V) -> Point3 {
    let [x, y, z] = point.map(|value| {
        let snapped = (value * 1024.0).round() / 1024.0;
        Real::from(Rational::from_f64_exact(snapped).unwrap())
    });
    Point3::new(x, y, z)
}

fn check_lines(fitted: &RollingSaddle) {
    let soup = Soup::of(fitted);
    let dec = decider();
    let directions = [["3", "7", "-5"], ["-2", "1", "9"], ["11", "-4", "2"]];
    let mut seen = 0;
    for i in 0..9 {
        let u = -0.93 + 0.44 * f64::from(i);
        let v = 0.1 + 0.09 * f64::from(i);
        let near = sample(fitted.carrier(), u, v, &dec);
        let direction = vector(directions[i as usize % 3]);
        let origin = rational_of(near);
        let mut dec = decider();
        let hits = fitted.line_hits(&origin, &direction, &mut dec).unwrap();
        let start = g(&origin);
        let along = f(&direction);
        let found = soup.crossings(start, along);
        assert_eq!(hits.len(), found.len(), "at {u} {v}: {found:?}");
        seen += hits.len();
        for (hit, want) in hits.iter().zip(&found) {
            let parameter = dot(sub(g(hit), start), along) / dot(along, along);
            assert!((parameter - want).abs() < 0.01, "{parameter} {want}");
        }
    }
    assert!(seen >= 9, "the lines went through the carrier: {seen}");
}

#[test]
fn a_line_meets_the_carrier_where_a_triangulation_says() {
    for index in 0..3 {
        check_lines(&fitted(index));
    }
}

// Unit axes here have square roots, which the span equations carry as variables.
#[test]
fn a_line_meets_a_carrier_at_an_angle_with_no_rational_unit_axis() {
    use BallSide::Outside;
    let wide = cylinder(["0", "0", "0"], ["0", "0", "2"], "10");
    for axis in [["0", "1", "1"], ["1", "1", "1"]] {
        let case = Case {
            name: "leaning",
            wide: wide.clone(),
            narrow: cylinder(["1", "0", "0"], axis, "3"),
            radius: "1",
            sides: [Outside, Outside],
            leaving: true,
        };
        check_lines(&fit(&case));
    }
}

#[test]
fn a_line_through_a_rail_point_hits_it_exactly() {
    let fitted = fitted(0);
    let mut dec = decider();
    let [wide, narrow] = fitted.rails();
    for (rail, u) in [(&wide, "1/3"), (&narrow, "7/5")] {
        let at = rail.point_at(&real(u), &mut dec).unwrap();
        let direction = vector(["1", "2", "-2"]);
        let hits = fitted.line_hits(&at, &direction, &mut dec).unwrap();
        let on = hits
            .iter()
            .filter(|hit| zero(&hit.minus(&at).norm_squared(), &mut dec))
            .count();
        assert_eq!(on, 1, "{u}");
    }
}

#[test]
fn rails_answer_incidence_normals_and_sides() {
    let case = &cases()[0];
    let fitted = fitted(0);
    let mut dec = decider();
    let support = SplineSupport::from(fitted.clone());
    let [wide, narrow] = fitted.rails();
    let at = narrow.point_at(&real("1/2"), &mut dec).unwrap();
    assert!(fitted.on_boundary(&at, &mut dec).unwrap());
    let normal = fitted.normal_at(&at, &mut dec).unwrap();
    assert!(zero(&normal.dot(&case.narrow.axis), &mut dec));
    let [centre, ..] = fitted.definition().section(&q("1/2")).unwrap();
    assert!(!fitted.on_boundary(&centre, &mut dec).unwrap());
    let wide_at = wide.point_at(&real("2"), &mut dec).unwrap();
    assert_eq!(fitted.side_at(&wide_at, &mut dec).unwrap(), Sign::Zero);
    assert_eq!(fitted.side_at(&centre, &mut dec).unwrap(), Sign::Positive);
    assert!(
        support
            .misses_point(&point(["100", "0", "0"]), &mut dec)
            .unwrap()
    );
}

#[test]
fn the_display_grid_closes_round_the_ring() {
    let fitted = fitted(0);
    let dec = decider();
    let grid = fitted.display_grid(&q("1/100"), 2_000_000, &dec).unwrap();
    assert!(grid.chord <= q("1/100"));
    // Between grid lines the carrier stays within the chord of its bilinear cell.
    let chord = grid.chord.to_f64_nearest();
    let at = |column: usize, row: usize| -> V {
        let p = &grid.points[column * grid.rows + row];
        std::array::from_fn(|i| p[i].mid().to_f64_nearest())
    };
    for i in 0..200 {
        let (x, y) = (0.137 + 0.618 * f64::from(i), 0.291 + 0.414 * f64::from(i));
        let (column, row) = (
            (x % (grid.columns - 1) as f64),
            (y % (grid.rows - 1) as f64),
        );
        let (c, r) = (column.floor() as usize, row.floor() as usize);
        let (s, t) = (column - c as f64, row - r as f64);
        let u = -1.0 + 4.0 * column / (grid.columns - 1) as f64;
        let v = row / (grid.rows - 1) as f64;
        let on = sample(fitted.carrier(), u, v, &dec);
        let near = add(
            add([0.0; 3], at(c, r), (1.0 - s) * (1.0 - t)),
            at(c + 1, r),
            s * (1.0 - t),
        );
        let near = add(
            add(near, at(c, r + 1), (1.0 - s) * t),
            at(c + 1, r + 1),
            s * t,
        );
        assert!(norm(sub(on, near)) <= chord, "{u} {v}");
    }
    for row in 0..grid.rows {
        let first = &grid.points[row];
        let last = &grid.points[(grid.columns - 1) * grid.rows + row];
        for k in 0..3 {
            assert!((&first[k].mid() - &last[k].mid()).abs() < q("1/1000000000"));
        }
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
