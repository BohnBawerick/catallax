//! The geometric gap and the change of parameter against known answers
//! (decision Q9, and section 4.2 step 2 of the Checkpoint 2 decisions).

use super::a_plus::copies;
use super::exact::*;
use super::*;
use crate::numerical::gap::{self, GapMeasure, GeometricMeasure, Hints, Limits, refit};

fn read(geometry: &EdgeGeometry) -> (GapMeasure, GeometricMeasure) {
    let mut dec = Decider::default();
    let copies = copies(geometry);
    let parametric =
        gap::measure(&copies, &Hints::default(), Limits::default(), &mut dec).expect("measure");
    let geometric = gap::measure_geometric(&copies, &Hints::default(), Limits::default(), &mut dec)
        .expect("geometric")
        .expect("a copy with its own u/v");
    (parametric, geometric)
}

fn fixed_case(name: &str) -> EdgeCase {
    fixed::cases()
        .into_iter()
        .find(|case| case.name == name)
        .expect("case")
}

/// The shares of a cubic over `[0, 1]` that runs `s + e (27/4) s^2 (1 - s)`:
/// the identity pushed ahead by at most `e`, at `s = 2/3`, with both ends
/// kept. It is the stretch of known-answer edge 8.
fn stretched(e: &Rational) -> [Rational; 4] {
    [
        Rational::zero(),
        ratio(1, 3),
        ratio(2, 3) + ratio(9, 4) * e.clone(),
        Rational::one(),
    ]
}

fn unit_knots() -> Vec<Rational> {
    clamped(3, &[], &Rational::zero(), &Rational::one())
}

/// A circular arc of radius `radius` over one radian, on a cylinder of
/// radius `radius + offset` about the same axis. The face copy keeps the
/// arc's height and runs its angle `slide + s + e (27/4) s^2 (1 - s)`.
///
/// The lift is the arc of the wider circle, so with no slide the geometric
/// gap is `offset` exactly, both ways, whatever the stretch. The
/// same-parameter gap is the chord between the angles `t` on one circle and
/// `u(t)` on the other, largest where `u - t` is.
fn arc_on_wider_cylinder(offset: &str, stretch: &str, slide: &str) -> EdgeGeometry {
    let radius = qi(10);
    let controls = stretched(&q(stretch))
        .into_iter()
        .map(|share| [q(slide) + share, Rational::zero()])
        .collect();
    EdgeGeometry {
        curve: circle(
            [Rational::zero(), Rational::zero(), Rational::zero()],
            [0, 0, 1],
            &radius,
            [1, 0, 0],
        ),
        parameters: [Rational::zero(), Rational::one()],
        surface: CaseSurface::Kernel(cylinder(
            [Rational::zero(), Rational::zero(), Rational::zero()],
            [0, 0, 1],
            &(&radius + &q(offset)),
            [1, 0, 0],
        )),
        face_copy: FaceCopy::Native(pcurve(3, unit_knots(), controls, None)),
    }
}

/// The chord between a point of a circle of radius `r` and a point of a
/// concentric circle of radius `r + d`, an angle `a` apart.
fn chord(r: f64, d: f64, a: f64) -> f64 {
    (d * d + 4.0 * r * (r + d) * (a / 2.0).sin().powi(2)).sqrt()
}

/// Known-answer edge 8, the stretch of solidcore #570: the face copy traces
/// the same straight line at another speed. The same-parameter gap is
/// 5.7e-5 mm and the geometric gap is zero. One span gives the copy the
/// edge's parameter, and its gap is rounding.
#[test]
fn a_stretched_copy_has_no_geometric_gap_and_one_span_gives_it_the_edge_parameter() {
    let case = fixed_case("08-parameter-stretch");
    let (parametric, geometric) = read(&case.geometry);
    assert!((parametric.raw_mm - 5.7e-5).abs() < 1e-12, "{parametric:?}");
    assert!(geometric.reported_mm() < 1e-12, "{geometric:?}");
    let mut dec = Decider::default();
    let next = refit::reparametrise(&copies(&case.geometry), 1e-7, &mut dec)
        .expect("reparametrise")
        .expect("a change of parameter");
    println!(
        "edge 8: same-parameter {:.3e} mm, geometric {:.3e} mm, reparametrised {:.3e} mm stored in {} span",
        parametric.raw_mm,
        geometric.reported_mm(),
        next.stored_mm,
        next.spans
    );
    assert_eq!(next.spans, 1);
    assert!(next.stored_mm < 1e-12, "{}", next.stored_mm);
}

/// Known-answer edge 14, the rational circle of solidcraft #8502: a circular
/// arc written as a rational spline, against a face copy that is straight in
/// the cylinder's angle. Both trace the same arc, so the geometric gap is
/// zero, and the 0.0436 mm of the same-parameter gap is all speed. The copy
/// is given the spline's parameter in spans that double until it stores rule
/// T1's least floor.
#[test]
fn a_circle_against_its_rational_copy_is_a_difference_of_speed() {
    let case = fixed_case("14-circle-against-rational-copy");
    let (parametric, geometric) = read(&case.geometry);
    assert!(parametric.raw_mm > 0.0436 && parametric.raw_mm < 0.0437);
    assert!(geometric.reported_mm() < 1e-12, "{geometric:?}");
    let mut dec = Decider::default();
    let floor = tolerance::floor_minimum_mm().to_f64_nearest();
    let next = refit::reparametrise(&copies(&case.geometry), floor, &mut dec)
        .expect("reparametrise")
        .expect("a change of parameter");
    println!(
        "edge 14: same-parameter {:.6e} mm, geometric {:.3e} mm, reparametrised {:.3e} mm stored in {} spans",
        parametric.raw_mm,
        geometric.reported_mm(),
        next.stored_mm,
        next.spans
    );
    assert!(next.stored_mm <= floor, "{}", next.stored_mm);
}

/// A geometric gap with a known answer under a difference of speed. The arc
/// lies 1e-4 mm inside its cylinder, and its face copy runs up to 1e-4 rad
/// ahead. The geometric gap is the 1e-4 mm, both ways. The same-parameter
/// gap is the chord, ten times as much. The reparametrised copy stores rule
/// T4's tolerance of the 1e-4 mm and stays on the old copy's u/v line.
#[test]
fn a_stretched_arc_on_a_wider_cylinder_keeps_its_radial_gap() {
    let (offset, stretch) = (1e-4, 1e-4);
    let geometry = arc_on_wider_cylinder("1e-4", "1e-4", "0");
    let (parametric, geometric) = read(&geometry);
    let chord_mm = chord(10.0, offset, stretch);
    assert!(
        (parametric.raw_mm - chord_mm).abs() < 1e-9 * chord_mm,
        "{parametric:?}"
    );
    for way in [geometric.forward_mm, geometric.backward_mm] {
        assert!((way - offset).abs() < 1e-9 * offset, "{geometric:?}");
    }
    // The done-when line of the checker, for the geometric gap: at least the
    // true gap and at most 1.1 times it.
    assert!(geometric.reported_mm() >= offset && geometric.reported_mm() <= 1.1 * offset);
    assert!(parametric.reported_mm() > gap::MARGIN * geometric.reported_mm());

    let mut dec = Decider::default();
    let target = gap::MARGIN * gap::MARGIN * geometric.reported_mm();
    let next = refit::reparametrise(&copies(&geometry), target, &mut dec)
        .expect("reparametrise")
        .expect("a change of parameter");
    println!(
        "arc: same-parameter {:.6e} mm, geometric {:.6e} mm forward and {:.6e} mm backward, reparametrised {:.6e} mm stored in {} spans",
        parametric.raw_mm, geometric.forward_mm, geometric.backward_mm, next.stored_mm, next.spans
    );
    assert!(next.stored_mm <= target, "{}", next.stored_mm);
    assert!(next.stored_mm >= gap::MARGIN * offset, "{}", next.stored_mm);
    // The new copy runs along the old one: the same height, the same two
    // ends, and an angle that ascends.
    let mut at = |t: Rational| {
        let point = next
            .curve
            .point_at(&Real::from(t), &mut dec)
            .expect("point");
        [point.x, point.y].map(|x| {
            dec.enclose_at(&x, BITS)
                .expect("enclose")
                .mid()
                .to_f64_nearest()
        })
    };
    let mut before = -1.0;
    for k in 0..=64 {
        let [u, v] = at(ratio(k, 64));
        assert_eq!(v, 0.0);
        assert!(u > before, "{u} after {before}");
        before = u;
    }
    assert_eq!(at(Rational::zero()), [0.0, 0.0]);
    assert_eq!(at(Rational::one()), [1.0, 0.0]);
}

/// A copy that is slid along its curve is not a copy at another speed: its
/// two ends are off by the slide, the ends are paired with each other, and
/// the geometric gap is the same-parameter gap. Nothing is reparametrised.
#[test]
fn a_copy_slid_along_its_curve_shows_at_its_ends() {
    let slide = 1e-4;
    let geometry = arc_on_wider_cylinder("0", "0", "1e-4");
    let (parametric, geometric) = read(&geometry);
    let chord_mm = chord(10.0, 0.0, slide);
    assert!((parametric.raw_mm - chord_mm).abs() < 1e-9 * chord_mm);
    assert!((geometric.raw_mm() - chord_mm).abs() < 1e-9 * chord_mm);
    assert!(parametric.reported_mm() <= gap::MARGIN * geometric.reported_mm());
}

/// "Geometric gap <= parametric gap always" (`checker-menu.md`, section
/// 2.6), on every fixed edge and on fresh ones: each foot starts from the
/// same-parameter pair and takes no step that does not come nearer.
#[test]
fn the_geometric_gap_is_never_read_larger_than_the_same_parameter_gap() {
    let seed = fresh::seed();
    println!("geometric seed={seed}");
    let cases = fixed::cases()
        .into_iter()
        .chain(fresh::cases(seed, 2))
        .chain(smooth::cases(seed, 2));
    let mut text = String::new();
    for case in cases.filter(|case| case.domain_exit.is_none()) {
        let (parametric, geometric) = read(&case.geometry);
        text.push_str(&format!(
            "edge={} parametric_mm={:.9e} forward_mm={:.9e} backward_mm={:.9e} geometric_evaluations={} parametric_evaluations={}\n",
            case.name,
            parametric.raw_mm,
            geometric.forward_mm,
            geometric.backward_mm,
            geometric.evaluations,
            parametric.evaluations,
        ));
        assert!(
            geometric.raw_mm() <= parametric.raw_mm + parametric.error_mm + geometric.error_mm,
            "seed {seed}, {}: {geometric:?} over {parametric:?}",
            case.name
        );
    }
    publish("known-answers-geometric.txt", &text);
}

/// The geometric reading runs under the allowance it is given, refuses by
/// name where that runs out, and charges the operation what it spent: the
/// lift, every unit it counted, and one deadline poll per charge.
#[test]
fn the_geometric_reading_counts_its_work_and_refuses_at_the_cap() {
    let geometry = arc_on_wider_cylinder("1e-4", "1e-4", "0");
    let capped = |work: u64| {
        with_operation_budget(DEFAULT_WORK, || {
            gap::measure_geometric(
                &copies(&geometry),
                &Hints::default(),
                Limits { work },
                &mut Decider::default(),
            )
        })
    };
    // No allowance at all refuses before the search: the lift's work alone.
    let (none, lift) = capped(0);
    assert_eq!(none.expect_err("capped").code, "numerical.gap-work-cap");
    let (found, charged) = capped(gap::EDGE_WORK);
    let found = found.expect("default").expect("a copy with its own u/v");
    assert!(found.work >= found.evaluations);
    assert!(charged.used >= lift.used + found.work);
    assert!(charged.used <= lift.used + 2 * found.work);
    let (short, spent) = capped(found.work / 2);
    assert_eq!(short.expect_err("capped").code, "numerical.gap-work-cap");
    assert!(spent.used < charged.used);
    // The same bits and the same work on a second run.
    let (again, recharged) = capped(gap::EDGE_WORK);
    let again = again.expect("default").expect("read");
    assert_eq!(again.forward_mm.to_bits(), found.forward_mm.to_bits());
    assert_eq!(again.backward_mm.to_bits(), found.backward_mm.to_bits());
    assert_eq!(again.work, found.work);
    assert_eq!(recharged.used, charged.used);
}

/// A curve projected on a plane has no u/v of its own: there is nothing to
/// read geometrically and nothing to reparametrise.
#[test]
fn a_projected_copy_has_no_geometric_reading() {
    let curve = circle(
        [Rational::zero(), Rational::zero(), Rational::zero()],
        [0, 0, 1],
        &qi(10),
        [1, 0, 0],
    );
    let surface = plane(
        [Rational::zero(), Rational::zero(), Rational::zero()],
        [0, 0, 1],
        [1, 0, 0],
    );
    let parameters = [Rational::zero(), Rational::one()];
    let copies = gap::EdgeCopies {
        curve: &curve,
        parameters: parameters.each_ref(),
        surface: gap::GapSurface::Kernel(&surface),
        face_copy: gap::FaceCopyRef::Projected,
    };
    let mut dec = Decider::default();
    assert!(
        gap::measure_geometric(&copies, &Hints::default(), Limits::default(), &mut dec)
            .expect("read")
            .is_none()
    );
    assert!(
        refit::reparametrise(&copies, 1e-7, &mut dec)
            .expect("read")
            .is_none()
    );
}
