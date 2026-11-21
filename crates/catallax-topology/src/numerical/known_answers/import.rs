//! The geometric gap and the change of parameter against known answers
//! (decision Q9, and section 4.2 step 2 of the Checkpoint 2 decisions).

use super::a_plus::copies;
use super::exact::*;
use super::*;
use crate::numerical::edge_gap::Accounts;
use crate::numerical::gap::{self, GapMeasure, GeometricMeasure, Hints, Limits, refit};
use crate::numerical::import_gap::{self, ImportGapOutcome, UseReading};

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
    let next = refit::reparametrise(&copies(&case.geometry), 1e-7, Limits::default(), &mut dec)
        .expect("reparametrise")
        .copy
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
    let next = refit::reparametrise(&copies(&case.geometry), floor, Limits::default(), &mut dec)
        .expect("reparametrise")
        .copy
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
    let next = refit::reparametrise(&copies(&geometry), target, Limits::default(), &mut dec)
        .expect("reparametrise")
        .copy
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

fn floor_mm() -> f64 {
    tolerance::floor_minimum_mm().to_f64_nearest()
}

/// The gap rule for one edge use on an account of `work` units that holds no
/// reading yet, and the work charged to the operation. The rule then takes
/// the first reading itself, on the same account, as admission takes it in
/// an import.
fn rule(geometry: &EdgeGeometry, work: u64) -> (Result<Option<UseReading>, Refusal>, u64) {
    let (read, charged) = with_operation_budget(DEFAULT_WORK, || {
        let mut accounts = Accounts::of_one_use(Limits { work });
        let mut account = accounts.draw(0)?;
        import_gap::read_use(
            &copies(geometry),
            &mut account,
            floor_mm(),
            &mut Decider::default(),
        )
    });
    (read, charged.used)
}

/// The three uses whose copies run at different speeds: the stretch of edge
/// 8, the rational circle of edge 14, and the stretched arc on its wider
/// cylinder.
fn uses_at_different_speeds() -> [(&'static str, EdgeGeometry); 3] {
    [
        (
            "08-parameter-stretch",
            fixed_case("08-parameter-stretch").geometry,
        ),
        (
            "14-circle-against-rational-copy",
            fixed_case("14-circle-against-rational-copy").geometry,
        ),
        (
            "stretched-arc-on-a-wider-cylinder",
            arc_on_wider_cylinder("1e-4", "1e-4", "0"),
        ),
    ]
}

/// Section 4.2 step 8: "A hard per-edge work cap refuses by name." Three
/// things are done with a use: its same-parameter gap is read, its geometric
/// gap is read, and its parameter is changed. This is the small-allowance
/// control that the three draw on one account and that nothing stands in for
/// the refusal. The first reading is the one admission takes in an import;
/// here the account holds none, so the rule takes it on the same account.
///
/// Read alone, each of the three fits an allowance one unit short of their
/// sum, so an allowance that started again for any of them would let the use
/// through. Under that one allowance the rule refuses by name. It refuses by
/// name under every smaller allowance tried as well: it never answers with
/// the copy as imported, nor with a candidate it read before the allowance
/// ran out. Every smaller allowance is tried for a use of under 1,024 units,
/// and about 512 of them, evenly spread, for a larger one. Under exactly the
/// sum the rule passes, with the answer and the work of the full cap: the
/// last step is the checker on the last candidate, and a reading of the
/// checker spends all it sets aside.
#[test]
fn everything_the_rule_does_with_a_use_draws_on_one_allowance() {
    let mut text = String::new();
    for (name, geometry) in &uses_at_different_speeds() {
        let (full, charged) = rule(geometry, gap::EDGE_WORK);
        let full = full.expect("the full cap").expect("a use the rule reads");
        assert_eq!(full.outcome, ImportGapOutcome::Reparametrised, "{name}");
        let next = full.reparametrised.as_ref().expect("a new copy");

        // The three steps, each counted alone at the full cap. The account
        // holds the first reading once and the two steps after it.
        let copies = copies(geometry);
        let mut dec = Decider::default();
        let target = (gap::MARGIN * gap::MARGIN * full.geometric.reported_mm()).max(floor_mm());
        let change = |limits: Limits, dec: &mut Decider| {
            refit::reparametrise(&copies, target, limits, dec).expect("alone")
        };
        let steps = [
            full.parametric.work,
            full.geometric.work,
            change(Limits::default(), &mut dec).work,
        ];
        let total: u64 = steps.iter().sum();
        assert_eq!(full.work, total, "{name}");
        assert!(total < gap::EDGE_WORK);

        // Alone, each step fits an allowance one unit short of the total.
        let short = Limits { work: total - 1 };
        let parametric = gap::measure(&copies, &Hints::default(), short, &mut dec).expect("alone");
        assert_eq!(parametric.work, steps[0]);
        let geometric = gap::measure_geometric(&copies, &Hints::default(), short, &mut dec)
            .expect("alone")
            .expect("a copy with its own u/v");
        assert_eq!(geometric.work, steps[1]);
        let alone = change(short, &mut dec);
        assert!(alone.copy.is_some());
        assert_eq!(alone.work, steps[2]);

        // Together they do not, and no smaller allowance gives any answer
        // but the refusal: a spread of them, and the ones at which a step
        // would just fit if the steps before it had cost nothing.
        let spread = (0..total).step_by((total / 512).max(1) as usize);
        let edges = [steps[0], steps[0] + steps[1], steps[2], total - 1];
        for work in spread.chain(edges) {
            let (refused, spent) = rule(geometry, work);
            match refused {
                Err(refusal) => assert_eq!(
                    refusal.code, "numerical.gap-work-cap",
                    "{name}: {work} of {total} units"
                ),
                Ok(answer) => panic!(
                    "{name}: {work} of {total} units gave {:?}",
                    answer.map(|answer| answer.outcome)
                ),
            }
            // Work done before a refusal is charged, and it is less.
            assert!(spent < charged, "{name}: {spent} against {charged}");
        }

        let (same, recharged) = rule(geometry, total);
        let same = same
            .expect("the sum of the steps is enough")
            .expect("a use the rule reads");
        let again = same.reparametrised.as_ref().expect("a new copy");
        assert_eq!(same.outcome, full.outcome);
        assert_eq!(again.stored_mm.to_bits(), next.stored_mm.to_bits());
        assert_eq!(again.spans, next.spans);
        assert!(again.curve.same_definition(&next.curve));
        assert_eq!(same.work, full.work);
        assert_eq!(recharged, charged);
        text.push_str(&format!(
            "use={name} same_parameter_work={} geometric_work={} change_of_parameter_work={} least_allowance={total} spans={} operation_work={charged}\n",
            steps[0], steps[1], steps[2], next.spans
        ));
    }
    print!("{text}");
    publish("known-answers-import-allowance.txt", &text);
}

/// The change of parameter counts its projections and the checker's work on
/// every candidate, and nothing starts again when the spans double.
///
/// The candidates of edge 14 are stepped through one at a time: a target no
/// copy misses stops the fit at its first candidate, and each candidate's
/// stored tolerance, lowered by one unit in the last place, stops it at the
/// next. From one stop to the next the fit's count rises by more than the
/// checker's own work on the new candidate, which is that work and the
/// projections of the doubling. So the count at the last candidate holds the
/// projections and the checker's work on every candidate before it.
///
/// The checker reads each candidate alone under an allowance of exactly its
/// own work there. The largest of those allowances, which is all a cap that
/// started again per candidate would ask for, is not enough for the fit: it
/// refuses by name, and so it does one unit short of its count.
#[test]
fn a_change_of_parameter_counts_its_projections_and_every_candidate() {
    let geometry = fixed_case("14-circle-against-rational-copy").geometry;
    let copies = copies(&geometry);
    let floor = floor_mm();
    let fit = |target: f64, work: u64| {
        with_operation_budget(DEFAULT_WORK, || {
            refit::reparametrise(&copies, target, Limits { work }, &mut Decider::default())
        })
    };
    let (full, charged) = fit(floor, gap::EDGE_WORK);
    let full = full.expect("the full cap");
    let full_spans = full.copy.as_ref().expect("a new copy").spans;

    let mut target = f64::INFINITY;
    let mut counted = 0;
    let mut largest = 0;
    let mut candidates = 0;
    loop {
        let stopped = fit(target, gap::EDGE_WORK).0.expect("the full cap");
        let stop = stopped.copy.as_ref().expect("a candidate");
        let alone = gap::measure(
            &gap::EdgeCopies {
                face_copy: gap::FaceCopyRef::Native(&stop.curve),
                ..copies
            },
            &Hints::default(),
            Limits::default(),
            &mut Decider::default(),
        )
        .expect("the checker on one candidate");
        assert_eq!(alone.stored_mm().to_bits(), stop.stored_mm.to_bits());
        gap::measure(
            &gap::EdgeCopies {
                face_copy: gap::FaceCopyRef::Native(&stop.curve),
                ..copies
            },
            &Hints::default(),
            Limits { work: alone.work },
            &mut Decider::default(),
        )
        .expect("one candidate under its own work");
        // The reading the fit keeps of its candidate is the checker's.
        assert_eq!(stop.reading.work, alone.work);
        assert!(
            stopped.work > counted + alone.work,
            "{} spans: {} after {counted} with {} for the candidate",
            stop.spans,
            stopped.work,
            alone.work
        );
        counted = stopped.work;
        largest = largest.max(alone.work);
        candidates += 1;
        if stop.spans == full_spans {
            break;
        }
        target = stop.stored_mm.next_down();
    }
    assert!(candidates > 1, "edge 14 doubles its spans");
    assert_eq!(counted, full.work);
    // The operation is charged at least what the fit counted.
    assert!(charged.used > full.work);

    let (one_candidate, _) = fit(floor, largest);
    assert_eq!(
        one_candidate.expect_err("capped").code,
        "numerical.gap-work-cap"
    );
    let (short, spent) = fit(floor, full.work - 1);
    assert_eq!(short.expect_err("capped").code, "numerical.gap-work-cap");
    assert!(spent.used < charged.used);
    // Nothing at all refuses before the first projection.
    let (none, _) = fit(floor, 0);
    assert_eq!(none.expect_err("capped").code, "numerical.gap-work-cap");
    println!(
        "edge 14: {candidates} candidates, {} units counted, {largest} for the largest candidate alone, {} charged to the operation",
        full.work, charged.used
    );
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
        curve: gap::GapCurve::Kernel(&curve),
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
        refit::reparametrise(&copies, 1e-7, Limits::default(), &mut dec)
            .expect("read")
            .copy
            .is_none()
    );
}

/// The arc of [`arc_on_wider_cylinder`] and its cylinder as a reader of a
/// foreign file holds them: two binary64 evaluators in one frame, with the
/// shape stated beside them. A circle has no degree and no breaks; the
/// cylinder closes on itself in u.
fn sampled_arc(offset: f64) -> (gap::SampledCurve, gap::SampledSurface) {
    let radius = 10.0_f64;
    let rounding = 4.0 * f64::EPSILON * (radius + offset);
    (
        gap::SampledCurve {
            point: std::sync::Arc::new(move |t: f64| {
                ([radius * t.cos(), radius * t.sin(), 0.0], rounding)
            }),
            breaks: Vec::new(),
            degree: 0,
            work: 0,
        },
        gap::SampledSurface {
            point: std::sync::Arc::new(move |uv: [f64; 2]| {
                let r = radius + offset;
                ([r * uv[0].cos(), r * uv[0].sin(), uv[1]], rounding)
            }),
            lines: [Vec::new(), Vec::new()],
            periods: [Some(std::f64::consts::TAU), None],
            degree: 0,
            work: 0,
        },
    )
}

/// A carrier with no kernel form is read through its own evaluator, by the
/// same search. The known answer is the arc on its wider cylinder: the
/// same-parameter gap is the chord, the geometric gap the 1e-4 mm, and the
/// rule gives the copy the arc's parameter. The plan is the kernel form's,
/// because the stated shape is the same: no degree, no breaks, one seam.
#[test]
fn sampled_carriers_are_read_by_the_same_search_and_the_same_rule() {
    let (offset, stretch) = (1e-4, 1e-4);
    let geometry = arc_on_wider_cylinder("1e-4", "1e-4", "0");
    let kernel = copies(&geometry);
    let (curve, surface) = sampled_arc(offset);
    let sampled = gap::EdgeCopies {
        curve: gap::GapCurve::Sampled(&curve),
        surface: gap::GapSurface::Sampled(&surface),
        ..kernel
    };
    let mut dec = Decider::default();
    let by_kernel =
        gap::measure(&kernel, &Hints::default(), Limits::default(), &mut dec).expect("kernel");
    let by_evaluator =
        gap::measure(&sampled, &Hints::default(), Limits::default(), &mut dec).expect("sampled");
    let chord_mm = chord(10.0, offset, stretch);
    assert!(
        (by_evaluator.raw_mm - chord_mm).abs() < 1e-9 * chord_mm,
        "{by_evaluator:?}"
    );
    assert_eq!(by_evaluator.plan, by_kernel.plan);
    assert_eq!(by_evaluator.method, gap::Method::Sampled);

    let (read, charged) = with_operation_budget(DEFAULT_WORK, || {
        import_gap::read_import_use(&sampled, Limits::default(), &mut Decider::default())
    });
    let read = read.expect("the rule");
    let geometric = read.geometric.expect("a geometric reading");
    for way in [geometric.forward_mm, geometric.backward_mm] {
        assert!((way - offset).abs() < 1e-9 * offset, "{geometric:?}");
    }
    assert_eq!(read.outcome, Some(ImportGapOutcome::Reparametrised));
    let next = read.new_copy().expect("a new copy");
    assert!(next.stored_mm >= gap::MARGIN * offset, "{}", next.stored_mm);
    assert!(
        next.stored_mm <= gap::MARGIN * gap::MARGIN * geometric.reported_mm(),
        "{}",
        next.stored_mm
    );
    assert_eq!(read.stored_mm(), next.stored_mm);
    // Every unit the use counts is charged to the operation, with one
    // deadline poll per charge on top.
    assert!(charged.used >= read.work, "{charged:?} under {}", read.work);
    println!(
        "sampled arc: same-parameter {:.6e} mm, geometric {:.6e} mm, reparametrised {:.6e} mm stored in {} spans, {} units",
        read.parametric.raw_mm,
        geometric.reported_mm(),
        next.stored_mm,
        next.spans,
        read.work
    );

    // One allowance for the use: a unit short of what it took refuses by
    // name, and nothing read before that stands in for the answer.
    let short = import_gap::read_import_use(
        &sampled,
        Limits {
            work: read.work - 1,
        },
        &mut Decider::default(),
    );
    assert_eq!(short.unwrap_err().code, "numerical.gap-work-cap");
    let full = import_gap::read_import_use(
        &sampled,
        Limits { work: read.work },
        &mut Decider::default(),
    )
    .expect("the least allowance");
    assert_eq!(full.stored_mm().to_bits(), read.stored_mm().to_bits());
}

/// An evaluator works in its owner's frame and a kernel carrier is shifted
/// to a frame of the checker's: one of each is refused by name, not read.
#[test]
fn a_sampled_carrier_is_not_read_against_a_kernel_carrier() {
    let geometry = arc_on_wider_cylinder("1e-4", "1e-4", "0");
    let kernel = copies(&geometry);
    let (curve, surface) = sampled_arc(1e-4);
    for mixed in [
        gap::EdgeCopies {
            curve: gap::GapCurve::Sampled(&curve),
            ..kernel
        },
        gap::EdgeCopies {
            surface: gap::GapSurface::Sampled(&surface),
            ..kernel
        },
    ] {
        let refused = gap::measure(
            &mixed,
            &Hints::default(),
            Limits::default(),
            &mut Decider::default(),
        );
        assert_eq!(refused.unwrap_err().code, "numerical.gap-unsupported");
    }
}

/// A sampled surface states its own cut lines: a face copy that crosses one
/// is cut there, as it is at a knot line of a kernel spline (section 4.2
/// step 3). A stated break of the 3D curve is a cut too.
#[test]
fn a_sampled_carrier_is_cut_at_the_breaks_it_states() {
    let geometry = arc_on_wider_cylinder("1e-4", "1e-4", "0");
    let kernel = copies(&geometry);
    let (mut curve, mut surface) = sampled_arc(1e-4);
    let plain = gap::EdgeCopies {
        curve: gap::GapCurve::Sampled(&curve),
        surface: gap::GapSurface::Sampled(&surface),
        ..kernel
    };
    let mut dec = Decider::default();
    let before = gap::plan(&plain, &Hints::default(), Limits::default(), &mut dec).expect("plan");
    curve.breaks = vec![0.25, 7.0];
    surface.lines[0] = vec![0.5];
    let cut = gap::EdgeCopies {
        curve: gap::GapCurve::Sampled(&curve),
        surface: gap::GapSurface::Sampled(&surface),
        ..kernel
    };
    let after = gap::plan(&cut, &Hints::default(), Limits::default(), &mut dec).expect("plan");
    assert_eq!(after.len(), before.len() + 2, "{after:?}");
    assert!(
        after
            .iter()
            .any(|piece| piece.piece[1].to_bits() == 0.25_f64.to_bits())
    );
    // The copy's angle reaches 0.5 near, not at, the parameter 0.5: it runs
    // up to 1e-4 rad ahead.
    assert!(
        after
            .iter()
            .any(|piece| (piece.piece[1] - 0.5).abs() < 1e-3 && (piece.piece[1] - 0.5).abs() > 0.0),
        "{after:?}"
    );
    // Every piece has the count of a shape with no degree.
    assert!(after.iter().all(|piece| piece.points == before[0].points));
}

/// Section 4.2 step 8: "A hard per-edge work cap refuses by name." A sampled
/// carrier states the work of one evaluation, and the checker takes that
/// from the use's one allowance with its own unit, before the evaluation.
///
/// The same use is read with evaluators that state no work and with ones
/// that state three and five units. The answer is the same and the work is
/// more. Under an allowance that fits the checker's own work exactly, the
/// use with the stated work refuses by name, with the operation's allowance
/// hardly touched. A unit short of its own total it refuses too, and under
/// exactly its total it gives the answer of the full cap.
#[test]
fn the_work_a_sampled_carrier_states_is_inside_the_edges_one_allowance() {
    let geometry = arc_on_wider_cylinder("1e-4", "1e-4", "0");
    let kernel = copies(&geometry);
    let read = |stated: [u64; 2], allowance: u64| {
        let (mut curve, mut surface) = sampled_arc(1e-4);
        (curve.work, surface.work) = (stated[0], stated[1]);
        let sampled = gap::EdgeCopies {
            curve: gap::GapCurve::Sampled(&curve),
            surface: gap::GapSurface::Sampled(&surface),
            ..kernel
        };
        with_operation_budget(DEFAULT_WORK, || {
            import_gap::read_import_use(
                &sampled,
                Limits { work: allowance },
                &mut Decider::default(),
            )
        })
    };
    let (free, _) = read([0, 0], gap::EDGE_WORK);
    let free = free.expect("free evaluators");
    let (costly, charged) = read([3, 5], gap::EDGE_WORK);
    let costly = costly.expect("costly evaluators");
    assert_eq!(costly.outcome, Some(ImportGapOutcome::Reparametrised));
    assert_eq!(costly.stored_mm().to_bits(), free.stored_mm().to_bits());
    assert!(
        costly.work > 4 * free.work,
        "{} and {}",
        costly.work,
        free.work
    );
    // Every unit of the use is charged to the operation.
    assert!(charged.used >= costly.work, "{charged:?}");
    println!(
        "sampled arc: {} units with free evaluators, {} with 3 and 5 units an evaluation, {} charged",
        free.work, costly.work, charged.used
    );

    // The checker's own work fits this allowance, and the total does not.
    let (alone, _) = read([0, 0], free.work);
    assert_eq!(
        alone
            .expect("the checker's work alone")
            .stored_mm()
            .to_bits(),
        free.stored_mm().to_bits()
    );
    let (refused, request) = read([3, 5], free.work);
    assert_eq!(refused.unwrap_err().code, "numerical.gap-work-cap");
    assert!(
        !request.exhausted && request.used <= 2 * free.work,
        "{request:?}"
    );

    let (short, _) = read([3, 5], costly.work - 1);
    assert_eq!(short.unwrap_err().code, "numerical.gap-work-cap");
    let (full, _) = read([3, 5], costly.work);
    assert_eq!(
        full.expect("the least allowance").stored_mm().to_bits(),
        costly.stored_mm().to_bits()
    );
}
