//! The exact band of a ball rolling over a mitre between two equal rounds.

use catallax_numbers::{Ball, BudgetPolicy, Decider, Outcome, Rational, Real, StrictSign};
use catallax_shapes::{
    Isometry, Point3, Pose, RPoint3, RVector3, RollingMitre, Rotation, SplineSupport, Surface,
    Vector3,
};

fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 200_000_000,
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
fn zero(value: &Real, dec: &mut Decider) -> bool {
    matches!(dec.sign(value), Outcome::ProvenCoincident(_))
}

struct Case {
    origin: [&'static str; 3],
    axes: [[&'static str; 3]; 2],
    round: &'static str,
    radius: &'static str,
    convex: bool,
}

// The narrow-mitres corner, a corner whose rounds meet at another angle, and
// the concave pocket form of the first.
fn cases() -> [Case; 3] {
    [
        Case {
            origin: ["1", "1", "9"],
            axes: [["1", "0", "0"], ["0", "1", "0"]],
            round: "1",
            radius: "1/2",
            convex: true,
        },
        Case {
            origin: ["2", "-1", "3"],
            axes: [["0", "0", "2"], ["3/5", "0", "4/5"]],
            round: "2",
            radius: "1/3",
            convex: true,
        },
        Case {
            origin: ["1", "1", "1"],
            axes: [["1", "0", "0"], ["0", "1", "0"]],
            round: "1",
            radius: "1/2",
            convex: false,
        },
    ]
}

fn build(case: &Case, dec: &mut Decider) -> RollingMitre {
    RollingMitre::new(
        point(case.origin),
        case.axes.map(vector),
        real(case.round),
        real(case.radius),
        case.convex,
        dec,
    )
    .unwrap()
}

// Squared distance from a point to the line through `origin` along unit `axis`.
fn from_axis(p: &Point3, origin: &Point3, axis: &Vector3) -> Real {
    let gap = p.minus(origin);
    &gap.norm_squared() - &gap.dot(axis).square()
}

const STEPS: [&str; 5] = ["0", "1/4", "1/2", "3/4", "1"];

#[test]
fn every_chart_point_is_a_ball_on_both_rounds_offsets() {
    let mut dec = decider();
    for case in cases() {
        let band = build(&case, &mut dec);
        let round = real(case.round);
        let offset = &round - &real(case.radius);
        for u in STEPS.map(q) {
            let centre = band.centre(&u);
            for axis in band.axes() {
                let gap = &from_axis(&centre, band.origin(), axis) - &offset.square();
                assert!(zero(&gap, &mut dec), "centre off an offset at u = {u}");
            }
            for v in STEPS.map(q) {
                let p = band.point(&u, &v);
                let ball = &p.minus(&centre).norm_squared() - &real(case.radius).square();
                assert!(zero(&ball, &mut dec), "point off its ball at ({u}, {v})");
            }
            for (i, v) in ["0", "1"].map(q).iter().enumerate() {
                let p = band.point(&u, v);
                let on = &from_axis(&p, band.origin(), &band.axes()[i]) - &round.square();
                assert!(zero(&on, &mut dec), "row {v} off round {i} at u = {u}");
            }
        }
    }
}

#[test]
fn the_band_starts_on_its_end_section_and_closes_at_the_apex() {
    let mut dec = decider();
    for case in cases() {
        let band = build(&case, &mut dec);
        let section = band.section(&mut dec).unwrap();
        for v in STEPS.map(q) {
            let p = band.point(&Rational::zero(), &v);
            let gap = p.minus(&section.center);
            assert!(zero(&gap.dot(&section.normal), &mut dec));
            assert!(zero(&(&gap.norm_squared() - &section.radius_sq), &mut dec));
            let top = band.point(&Rational::one(), &v).minus(&band.apex());
            assert!(
                zero(&top.norm_squared(), &mut dec),
                "column u = 1 is not the apex"
            );
        }
    }
}

#[test]
fn each_contact_row_is_its_rail_ellipse() {
    let mut dec = decider();
    for case in cases() {
        let band = build(&case, &mut dec);
        for (i, v) in ["0", "1"].map(q).iter().enumerate() {
            let rail = band.rail(i, &mut dec).unwrap();
            let minor = rail.normal.cross(&rail.major_dir);
            for u in STEPS.map(q) {
                let gap = band.point(&u, v).minus(&rail.center);
                assert!(zero(&gap.dot(&rail.normal), &mut dec));
                let a = &gap.dot(&rail.major_dir).square()
                    / &(&rail.major_dir.norm_squared() * &rail.major_sq);
                let b = &gap.dot(&minor).square() / &(&minor.norm_squared() * &rail.minor_sq);
                assert!(
                    zero(&(&(&a + &b) - &Real::one()), &mut dec),
                    "row {v} off rail {i}"
                );
            }
        }
    }
}

#[test]
fn rails_and_section_run_from_the_strip_end_toward_the_apex() {
    let mut dec = decider();
    for case in cases() {
        let band = build(&case, &mut dec);
        turns(&band, &real(case.radius), &mut dec);
    }
}

// The rails run from their feet to the apex and the section from the first
// rail's foot through the band to the second's.
fn turns(band: &RollingMitre, radius: &Real, dec: &mut Decider) {
    for (i, v) in ["0", "1"].map(q).iter().enumerate() {
        let rail = band.rail(i, dec).unwrap();
        // The seam is the foot, and a quarter turn on is the apex.
        let seam = rail.center.plus(&rail.major_dir);
        let foot = band.point(&Rational::zero(), v);
        assert!(
            zero(&seam.minus(&foot).norm_squared(), dec),
            "rail {i} seam"
        );
        assert!(zero(
            &(&rail.major_dir.norm_squared() - &rail.major_sq),
            dec
        ));
        let turn = rail.normal.cross(&rail.major_dir);
        let quarter = rail
            .center
            .plus(&turn.scaled(&(&rail.minor_sq.sqrt() / &turn.norm_squared().sqrt())));
        assert!(
            zero(&quarter.minus(&band.apex()).norm_squared(), dec),
            "rail {i} apex"
        );
    }
    // The section turns from the first rail's foot toward the second
    // through the band, not around the ball's other side.
    let section = band.section(dec).unwrap();
    let first = band
        .point(&Rational::zero(), &Rational::zero())
        .minus(&section.center);
    let middle = band
        .point(&Rational::zero(), &q("1/2"))
        .minus(&section.center);
    let seam = section.ref_dir.scaled(radius);
    assert!(
        zero(&first.minus(&seam).norm_squared(), dec),
        "section seam"
    );
    let turning = section.normal.cross(&section.ref_dir);
    assert!(matches!(
        dec.sign(&middle.dot(&turning)),
        Outcome::Proven(StrictSign::Positive)
    ));
}

// The written chart curves: each conic at its angle `x` from its foot is the
// band at `u = tan(x / 2)` on a rail and at `v = tan(x / 2) / l` on the end
// section, `l = across / along`.
fn charts(band: &RollingMitre, dec: &mut Decider) {
    let ratio = band.across() / band.along();
    let turn = |tangent: &Real| {
        let wide = &Real::one() + &tangent.square();
        (
            &(&Real::one() - &tangent.square()) / &wide,
            &(&Real::from_i64(2) * tangent) / &wide,
        )
    };
    for (i, v) in ["0", "1"].map(q).iter().enumerate() {
        let rail = band.rail(i, dec).unwrap();
        let side = rail.normal.cross(&rail.major_dir);
        let side = side.scaled(&(&rail.minor_sq.sqrt() / &side.norm_squared().sqrt()));
        for u in STEPS.map(q) {
            let (cos, sin) = turn(&Real::from(u.clone()));
            let on = rail
                .center
                .plus(&rail.major_dir.scaled(&cos))
                .plus(&side.scaled(&sin));
            let gap = on.minus(&band.point(&u, v)).norm_squared();
            assert!(zero(&gap, dec), "rail {i} chart at u = {u}");
        }
    }
    let section = band.section(dec).unwrap();
    let radius = section.radius_sq.sqrt();
    let side = section.normal.cross(&section.ref_dir);
    let side = side.scaled(&(&radius / &side.norm_squared().sqrt()));
    let seam = section
        .ref_dir
        .scaled(&(&radius / &section.ref_dir.norm_squared().sqrt()));
    for v in STEPS.map(q) {
        let (cos, sin) = turn(&(&ratio * &Real::from(v.clone())));
        let on = section
            .center
            .plus(&seam.scaled(&cos))
            .plus(&side.scaled(&sin));
        let gap = on.minus(&band.point(&Rational::zero(), &v)).norm_squared();
        assert!(zero(&gap, dec), "section chart at v = {v}");
    }
}

// The sign of the chart normal `S_u x S_v` along the ray from the ball's
// centre at an interior chart point.
fn chart_side(band: &RollingMitre, dec: &Decider) -> StrictSign {
    let (u, v) = (q("1/2"), q("1/2"));
    let [du, dv] = band
        .carrier()
        .derivatives(
            &Ball::from_rational(&u, 128),
            &Ball::from_rational(&v, 128),
            dec,
        )
        .unwrap();
    let ray = band.point(&u, &v).minus(&band.centre(&u));
    let ray = [&ray.x, &ray.y, &ray.z].map(|c| dec.enclose_at(c, 128).unwrap());
    let normal = [
        &(&du[1] * &dv[2]) - &(&du[2] * &dv[1]),
        &(&du[2] * &dv[0]) - &(&du[0] * &dv[2]),
        &(&du[0] * &dv[1]) - &(&du[1] * &dv[0]),
    ];
    let along = &(&(&normal[0] * &ray[0]) + &(&normal[1] * &ray[1])) + &(&normal[2] * &ray[2]);
    along
        .sign()
        .expect("the chart normal is clear of the section plane")
}

fn reflected(band: &RollingMitre, map: &Isometry) -> RollingMitre {
    let Surface::BSpline(SplineSupport::RollingMitre(band)) =
        map.surface(&Surface::BSpline(SplineSupport::RollingMitre(band.clone())))
    else {
        unreachable!("a mapped band stays a band");
    };
    band
}

fn motions() -> [Isometry; 3] {
    let pose = Pose::from_rational(
        Rotation::from_quaternion(q("1"), q("2"), q("3"), q("4")).unwrap(),
        RVector3::new(q("1000"), q("-500"), q("10")),
    );
    let mirror =
        Isometry::reflection(&RPoint3::origin(), &RVector3::new(q("1"), q("2"), q("3"))).unwrap();
    let posed = Isometry::from(&pose);
    [posed.clone(), mirror.clone(), posed.then(&mirror)]
}

#[test]
fn moved_and_mirrored_bands_keep_their_turns_charts_and_side() {
    let mut dec = decider();
    for case in cases() {
        let band = build(&case, &mut dec);
        assert!(!band.reversing());
        charts(&band, &mut dec);
        assert_eq!(chart_side(&band, &dec), StrictSign::Positive);
        for map in motions() {
            let moved = reflected(&band, &map);
            assert_eq!(moved.reversing(), map.reversing());
            assert_eq!(moved.convex(), band.convex());
            for (u, v) in [("0", "0"), ("1/2", "1/3"), ("1", "1")] {
                let (u, v) = (q(u), q(v));
                let image = map.apply_point(&band.point(&u, &v));
                assert!(zero(
                    &image.minus(&moved.point(&u, &v)).norm_squared(),
                    &mut dec
                ));
            }
            turns(&moved, &real(case.radius), &mut dec);
            charts(&moved, &mut dec);
            // A reflection turns the chart's normal toward the ball.
            let side = if map.reversing() {
                StrictSign::Negative
            } else {
                StrictSign::Positive
            };
            assert_eq!(chart_side(&moved, &dec), side);
        }
    }
}

#[test]
fn the_carrier_is_the_law() {
    let mut dec = decider();
    for case in cases() {
        let band = build(&case, &mut dec);
        assert!(
            band.carrier()
                .weights()
                .iter()
                .all(|w| matches!(dec.sign(w), Outcome::Proven(StrictSign::Positive)))
        );
        for u in STEPS.map(q) {
            for v in STEPS.map(q) {
                let at = |r: &Rational| Ball::from_rational(r, 256);
                let enclosed = band.carrier().evaluate(&at(&u), &at(&v), &dec).unwrap();
                let exact = band.point(&u, &v);
                for (ball, coordinate) in enclosed.iter().zip([&exact.x, &exact.y, &exact.z]) {
                    let value = dec.enclose_at(coordinate, 256).unwrap();
                    let gap = ball - &value;
                    assert!(
                        gap.lo() <= Rational::zero() && gap.hi() >= Rational::zero(),
                        "carrier misses the law at ({u}, {v})"
                    );
                    assert!(gap.hi() - gap.lo() < q("1/1000000000000"));
                }
            }
        }
    }
}

#[test]
fn the_display_grid_meets_its_chord_bound() {
    let mut dec = decider();
    for case in cases() {
        let band = build(&case, &mut dec);
        let tolerance = q("1/1000");
        let grid = band.display_grid(&tolerance, 100_000, &dec).unwrap();
        assert!(grid.chord <= tolerance);
        assert_eq!(grid.points.len(), grid.columns * grid.rows);
    }
}

#[test]
fn a_ball_as_wide_as_the_rounds_and_parallel_axes_refuse() {
    let mut dec = decider();
    let wide = RollingMitre::new(
        point(["0", "0", "0"]),
        [vector(["1", "0", "0"]), vector(["0", "1", "0"])],
        real("1"),
        real("1"),
        true,
        &mut dec,
    );
    assert!(wide.is_err());
    let parallel = RollingMitre::new(
        point(["0", "0", "0"]),
        [vector(["1", "0", "0"]), vector(["2", "0", "0"])],
        real("1"),
        real("1/2"),
        true,
        &mut dec,
    );
    assert!(parallel.is_err());
}

#[test]
fn boundary_points_have_exact_normals_and_interior_points_refuse() {
    let mut dec = decider();
    for case in cases() {
        let band = build(&case, &mut dec);
        let radius = real(case.radius);
        for u in STEPS.map(q) {
            for v in ["0", "1"].map(q) {
                let p = band.point(&u, &v);
                assert!(band.on_boundary(&p, &mut dec).unwrap());
                let normal = band.normal_at(&p, &mut dec).unwrap();
                let expected = band.centre(&u).minus(&p).scaled(&(&Real::one() / &radius));
                assert!(zero(&normal.minus(&expected).norm_squared(), &mut dec));
                assert_eq!(
                    band.side_at(&p, &mut dec).unwrap(),
                    catallax_numbers::Sign::Zero
                );
            }
        }
        for v in STEPS.map(q) {
            let p = band.point(&Rational::zero(), &v);
            assert!(
                band.on_boundary(&p, &mut dec).unwrap(),
                "end section at v = {v}"
            );
        }
        let inner = band.point(&q("1/2"), &q("1/2"));
        assert!(!band.on_boundary(&inner, &mut dec).unwrap());
        assert!(band.normal_at(&inner, &mut dec).is_err());
    }
}

#[test]
fn side_reports_the_tube_near_the_band() {
    let mut dec = decider();
    for case in cases() {
        let band = build(&case, &mut dec);
        for (u, v) in [
            ("1/3", "1/2"),
            ("1/2", "1/4"),
            ("2/3", "3/4"),
            ("1/8", "1/2"),
        ] {
            let (u, v) = (q(u), q(v));
            let centre = band.centre(&u);
            let ray = band.point(&u, &v).minus(&centre);
            let inside = centre.plus(&ray.scaled(&real("9/10")));
            let outside = centre.plus(&ray.scaled(&real("11/10")));
            assert_eq!(
                band.side_at(&inside, &mut dec).unwrap(),
                catallax_numbers::Sign::Positive
            );
            assert_eq!(
                band.side_at(&outside, &mut dec).unwrap(),
                catallax_numbers::Sign::Negative
            );
        }
    }
}

#[test]
fn a_line_through_a_section_centre_meets_the_band_once() {
    let mut dec = decider();
    for case in cases() {
        let band = build(&case, &mut dec);
        for (u, v) in [("1/3", "1/2"), ("1/2", "1/4"), ("0", "1/2")] {
            let (u, v) = (q(u), q(v));
            let centre = band.centre(&u);
            let p = band.point(&u, &v);
            let hits = band
                .line_hits(&centre, &p.minus(&centre), &mut dec)
                .unwrap();
            assert_eq!(hits.len(), 1, "hits at ({u}, {v})");
            assert!(zero(&hits[0].minus(&p).norm_squared(), &mut dec));
        }
        let far = band.origin().plus(&vector(["100", "100", "100"]));
        assert!(
            band.line_hits(&far, &vector(["1", "-2", "3"]), &mut dec)
                .unwrap()
                .is_empty()
        );
    }
}

#[test]
fn the_extremes_bound_the_display_grid() {
    let mut dec = decider();
    for case in cases() {
        let band = build(&case, &mut dec);
        let extremes = band.axis_extrema(&mut dec).unwrap();
        let grid = band.display_grid(&q("1/100"), 100_000, &dec).unwrap();
        for axis in 0..3 {
            let values: Vec<Rational> = extremes
                .iter()
                .map(|p| {
                    let c = [&p.x, &p.y, &p.z][axis];
                    dec.enclose_at(c, 128).unwrap().hi()
                })
                .collect();
            let lows: Vec<Rational> = extremes
                .iter()
                .map(|p| dec.enclose_at([&p.x, &p.y, &p.z][axis], 128).unwrap().lo())
                .collect();
            let high = values.iter().max().unwrap().clone();
            let low = lows.iter().min().unwrap().clone();
            let slack = q("1/1000000000");
            for point in &grid.points {
                assert!(
                    point[axis].hi() <= &high + &slack,
                    "grid above the extremes"
                );
                assert!(point[axis].lo() >= &low - &slack, "grid below the extremes");
            }
        }
    }
}

#[test]
fn lines_through_the_collapsed_apex_section_meet_only_the_apex() {
    let mut dec = decider();
    let band = RollingMitre::new(
        point(["0", "0", "0"]),
        [vector(["1", "0", "0"]), vector(["0", "1", "0"])],
        real("1"),
        real("1/2"),
        true,
        &mut dec,
    )
    .unwrap();
    // Both contact rays are `y` at the apex section, whose circle about
    // `(0, 0, 1/2)` meets the band only at the apex `(0, 0, 1)`. The upright
    // line also meets that circle at the origin, and a line through its
    // centre tilted across the mitre plane crosses it twice, simply, at
    // `(0, 0, 1/2) -+ (-1, 1, 1) / (2 sqrt 3)`, with no other band point.
    let upright = (point(["0", "0", "0"]), vector(["0", "0", "1"]), 1);
    let across = (point(["0", "0", "1/2"]), vector(["-1", "1", "1"]), 0);
    let identity = Isometry::identity();
    for map in [identity].into_iter().chain(motions()) {
        let moved = reflected(&band, &map);
        for (line, (origin, direction, count)) in
            [upright.clone(), across.clone()].into_iter().enumerate()
        {
            let hits = moved
                .line_hits(
                    &map.apply_point(&origin),
                    &map.apply_vector(&direction),
                    &mut dec,
                )
                .unwrap();
            assert_eq!(hits.len(), count, "line {line}");
            for hit in &hits {
                assert!(zero(&hit.minus(&moved.apex()).norm_squared(), &mut dec));
            }
        }
    }
}
