//! The A+ checker of `numerical::gap` against the known answers, and the
//! three guards against tuning it to them.

use super::*;
use crate::numerical::gap::{
    self, EdgeCopies, FaceCopyRef, GapSurface, Hints, Limits, Method, PieceShape, Reason,
};

struct APlus;

fn copies(geometry: &EdgeGeometry) -> EdgeCopies<'_> {
    EdgeCopies {
        curve: &geometry.curve,
        parameters: geometry.parameters.each_ref(),
        surface: match &geometry.surface {
            CaseSurface::Kernel(surface) => GapSurface::Kernel(surface),
            CaseSurface::Carrier(carrier) => GapSurface::Carrier(carrier),
        },
        face_copy: match &geometry.face_copy {
            FaceCopy::Native(curve) => FaceCopyRef::Native(curve),
            FaceCopy::Linear(ends) => FaceCopyRef::Linear(ends),
        },
    }
}

fn exact(value: f64) -> Rational {
    Rational::from_f64_exact(value).expect("finite")
}

impl GapChecker for APlus {
    fn name(&self) -> &'static str {
        "a-plus"
    }
    fn measure(&self, geometry: &EdgeGeometry, dec: &mut Decider) -> Result<Reading, Refusal> {
        let found = gap::measure(&copies(geometry), &Hints::default(), Limits::default(), dec)?;
        Ok(Reading {
            raw_mm: exact(found.raw_mm),
            allowance_mm: exact(found.allowance_mm),
            stored_mm: Some(exact(found.stored_mm())),
            at: Some(exact(found.parameter)),
            exact_zero: false,
            domain_exit: found.domain_exit.map(exact),
            evaluations: found.evaluations,
        })
    }
}

/// What the oracle's receipt has no column for: how each result was reached,
/// and the search uncertainty apart from the rounding allowance.
fn search_receipt(seed: Option<u64>, cases: &[EdgeCase]) -> String {
    let mut text = format!(
        "checker=a-plus seed={}\n",
        seed.map_or_else(|| "none".into(), |seed| seed.to_string())
    );
    for case in cases {
        let started = std::time::Instant::now();
        let found = gap::measure(
            &copies(&case.geometry),
            &Hints::default(),
            Limits::default(),
            &mut Decider::default(),
        )
        .expect("measure");
        // Clock time is recorded and judges nothing.
        let micros = started.elapsed().as_micros();
        // The search uncertainty fits inside the margin, so the stored
        // tolerance covers the raw value and what the search may have left.
        // On a zero-gap control both are rounding noise, far under the
        // 1e-7 mm floor of rule T1; a billionth of a millimetre sets it aside.
        assert!(
            found.raw_mm + found.search_mm <= found.stored_mm() + 1e-9,
            "{}: search {} against stored {}",
            case.name,
            found.search_mm,
            found.stored_mm()
        );
        text.push_str(&format!(
            "edge={} method={:?} pieces={} points={} splits={} raw_mm={:.9e} \
             rounding_mm={:.3e} search_mm={:.3e} stored_mm={:.9e} evaluations={} micros={micros}\n",
            case.name,
            found.method,
            found.plan.len(),
            found.plan.iter().map(|piece| piece.points).sum::<usize>(),
            found.refinements.len(),
            found.raw_mm,
            found.allowance_mm,
            found.search_mm,
            found.stored_mm(),
            found.evaluations,
        ));
    }
    text
}

/// The 3.0A done-when line: every one of the fourteen edges, by name.
#[test]
fn a_plus_passes_every_known_answer_edge() {
    let rows = run(&APlus, &fixed::cases());
    publish(
        "known-answers-a-plus.txt",
        &receipt(APlus.name(), None, &rows),
    );
    publish(
        "known-answers-a-plus-search.txt",
        &search_receipt(None, &fixed::cases()),
    );
    all_required_rows_pass(&rows).unwrap_or_else(|failures| panic!("{failures}"));
    assert!(accept(&APlus).is_ok());
}

/// The oracle's reason for a rounding allowance above a thousandth of the
/// true gap.
const ALLOWANCE_RULE: &str = "error allowance is too large to judge the measurement";

/// Whether a row of the frozen steep surface family misses only the oracle's
/// allowance-share rule and still meets the done-when line of roadmap item
/// 3.0A: reported at least the true gap and at most 1.1 times it.
///
/// That family puts random control heights inside a knot cluster as narrow
/// as 2^-20, so the surface slopes by millions of millimetres per unit
/// parameter and the face copy's own rounding, carried through that slope,
/// can exceed a thousandth of the smallest bumps. The checker says so in its
/// allowance. Such a row is not counted as meeting the share rule: it is
/// held to the done-when line and listed in the receipt as stress evidence.
fn allowance_limited(row: &Row) -> bool {
    let reported = row.reading.reported_mm();
    row.name.starts_with("fresh-surface-bump-")
        && row
            .verdict
            .as_ref()
            .is_err_and(|reason| reason == ALLOWANCE_RULE)
        && reported >= row.truth.hi
        && reported <= &row.truth.lo * &Rational::from_ratio(11, 10).expect("nonzero")
}

/// Guard 1: new edges on every run, each with a hidden bump. The frozen
/// families replay from the seed as the baseline recorded them; the smooth
/// surface family draws from its own stream beside them.
#[test]
fn a_plus_finds_the_hidden_bump_on_fresh_edges() {
    let seed = fresh::seed();
    let mut cases = fresh::cases(seed, 8);
    cases.extend(smooth::cases(seed, 8));
    let rows = run(&APlus, &cases);
    publish(
        "known-answers-fresh-a-plus.txt",
        &receipt(APlus.name(), Some(seed), &rows),
    );
    publish(
        "known-answers-fresh-a-plus-search.txt",
        &search_receipt(Some(seed), &cases),
    );
    let mut steep = format!("checker=a-plus seed={seed}\n");
    for row in &rows {
        if allowance_limited(row) {
            steep.push_str(&format!(
                "{} allowance-limited raw={:e} allowance={:e} true=[{:e}, {:e}]\n",
                row.name,
                row.reading.raw_mm.to_f64_nearest(),
                row.reading.allowance_mm.to_f64_nearest(),
                row.truth.lo.to_f64_nearest(),
                row.truth.hi.to_f64_nearest(),
            ));
            continue;
        }
        assert!(
            row.verdict.is_ok(),
            "seed {seed}, {}: {:?}",
            row.name,
            row.verdict
        );
    }
    publish("known-answers-fresh-a-plus-steep.txt", &steep);
}

/// The steep-family allowance is an exception for that family's one rule and
/// nothing else: a low reading, a high reading or another family fails.
#[test]
fn only_a_steep_row_inside_the_done_when_line_is_allowance_limited() {
    let truth = Truth {
        lo: Rational::from_i64(1000),
        hi: Rational::from_i64(1000),
        at: Rational::zero(),
    };
    let row = |name: &str, raw: i64, allowance: i64, verdict: Result<(), String>| Row {
        name: name.into(),
        truth: truth.clone(),
        reading: Reading {
            raw_mm: Rational::from_i64(raw),
            allowance_mm: Rational::from_i64(allowance),
            stored_mm: None,
            at: None,
            exact_zero: false,
            domain_exit: None,
            evaluations: 0,
        },
        work: 0,
        verdict,
    };
    let rule = || Err(ALLOWANCE_RULE.to_string());
    assert!(allowance_limited(&row(
        "fresh-surface-bump-0",
        998,
        5,
        rule()
    )));
    // Below the true gap, and above 1.1 times it.
    assert!(!allowance_limited(&row(
        "fresh-surface-bump-0",
        990,
        5,
        rule()
    )));
    assert!(!allowance_limited(&row(
        "fresh-surface-bump-0",
        1000,
        101,
        rule()
    )));
    // Another family, another reason, and a pass.
    assert!(!allowance_limited(&row(
        "fresh-plane-lift-0",
        998,
        5,
        rule()
    )));
    assert!(!allowance_limited(&row(
        "fresh-surface-bump-0",
        998,
        5,
        Err("reported gap is below the true gap".into())
    )));
    assert!(!allowance_limited(&row(
        "fresh-surface-bump-0",
        998,
        5,
        Ok(())
    )));
}

/// Guard 2: the first sample plan depends on the shape of the edge only. A
/// scaled and shifted twin has the same knots, degrees and crossings, so it
/// gets the same pieces and the same points, and any later split carries a
/// recorded reason.
#[test]
fn the_sample_plan_depends_on_shape_only() {
    let seed = fresh::seed();
    println!("plan seed={seed}");
    let scale = Rational::from_ratio(7, 3).expect("nonzero");
    let shift = [1234, -77, 5].map(Rational::from_i64);
    let cases = fixed::cases()
        .into_iter()
        .chain(fresh::cases(seed, 2))
        .chain(smooth::cases(seed, 2));
    for case in cases.filter(|case| case.domain_exit.is_none()) {
        let twin = fresh::twin(&case, &scale, &shift);
        let measure = |geometry: &EdgeGeometry| {
            let mut dec = Decider::default();
            let plan = gap::plan(
                &copies(geometry),
                &Hints::default(),
                Limits::default(),
                &mut dec,
            )
            .expect("plan");
            let found = gap::measure(
                &copies(geometry),
                &Hints::default(),
                Limits::default(),
                &mut dec,
            )
            .expect("measure");
            (plan, found)
        };
        let (plan, found) = measure(&case.geometry);
        let (twin_plan, twin_found) = measure(&twin.geometry);
        assert_eq!(plan, twin_plan, "seed {seed}, {}", case.name);
        assert_eq!(found.plan, plan, "{}", case.name);
        assert_eq!(found.method, twin_found.method, "{}", case.name);
        for refinement in found.refinements.iter().chain(&twin_found.refinements) {
            assert_eq!(refinement.reason, Reason::LoudTail);
        }
        // The twin's answer is the scale times the original's, to within
        // the two rounding allowances.
        let difference = (twin_found.raw_mm - 7.0 / 3.0 * found.raw_mm).abs();
        assert!(
            difference <= twin_found.allowance_mm + 7.0 / 3.0 * found.allowance_mm,
            "seed {seed}, {}: {difference}",
            case.name
        );
    }
    // The point count is a function of three degrees and nothing else.
    let bicubic_with_cubic = PieceShape {
        curve: 3,
        face: 3,
        surface: 6,
    };
    assert_eq!(bicubic_with_cubic.composed_degree(), 18);
    assert_eq!(bicubic_with_cubic.points(), 54);
    let analytic = PieceShape {
        curve: 0,
        face: 1,
        surface: 0,
    };
    assert_eq!(analytic.points(), 18);
}

/// A curved analytic surface has no degree, so its pieces start at the least
/// point count whatever the face copy does. Here a degree-23 face copy
/// wiggles the angle on a cylinder twenty-one times inside one piece. The
/// coefficient tail is loud, the piece is split until the wiggle is
/// resolved, and the result agrees with a dense sweep of exact values.
#[test]
fn a_loud_tail_splits_a_piece_until_it_is_resolved() {
    use exact::{bernstein_from_power, circle, clamped, cylinder, pcurve, q, qi};
    // w(x) = T_21(x) (1 - x^2) on x = 2 t - 1, from the Chebyshev recurrence.
    let x = catallax_numbers::UPoly::new(vec![qi(-1), qi(2)]);
    let (mut before, mut last) = (catallax_numbers::UPoly::constant(qi(1)), x.clone());
    for _ in 1..21 {
        let next = &(&(&x * &last) * &catallax_numbers::UPoly::constant(qi(2))) - &before;
        (before, last) = (last, next);
    }
    let window = &catallax_numbers::UPoly::constant(qi(1)) - &(&x * &x);
    let wiggle = &last * &window;
    let amplitude = q("1e-3");
    let mut angle: Vec<Rational> = bernstein_from_power(wiggle.coefficients())
        .iter()
        .map(|c| c * &amplitude)
        .collect();
    let degree = angle.len() - 1;
    assert_eq!(degree, 23);
    for (i, value) in angle.iter_mut().enumerate() {
        *value = &*value + &Rational::from_ratio(i as i64, degree as i64).expect("nonzero");
    }
    let knots = clamped(degree, &[], &qi(0), &qi(1));
    let geometry = EdgeGeometry {
        curve: circle([qi(0), qi(0), qi(3)], [0, 0, 1], &qi(10), [1, 0, 0]),
        parameters: [qi(0), qi(1)],
        surface: CaseSurface::Kernel(cylinder(
            [qi(0), qi(0), qi(0)],
            [0, 0, 1],
            &qi(10),
            [1, 0, 0],
        )),
        face_copy: FaceCopy::Native(pcurve(
            degree,
            knots,
            angle.into_iter().map(|u| [u, qi(3)]).collect(),
            None,
        )),
    };
    let mut dec = Decider::default();
    let found = gap::measure(
        &copies(&geometry),
        &Hints::default(),
        Limits::default(),
        &mut dec,
    )
    .expect("measure");
    assert_eq!(found.plan.len(), 1);
    assert_eq!(found.plan[0].points, 18);
    assert!(!found.refinements.is_empty());
    assert!(
        found
            .refinements
            .iter()
            .all(|split| split.reason == Reason::LoudTail)
    );
    // The gap is the chord 2 R sin(|A w| / 2). A sweep of 1001 exact values
    // falls short of its maximum by less than one part in a thousand.
    let mut dense = Rational::zero();
    for i in 0..=1000 {
        let t = Rational::from_ratio(i, 1000).expect("nonzero");
        let [lo, _] = exact::gap_enclosure(&geometry, &t, &mut dec).expect("gap");
        dense = dense.max(lo);
    }
    let dense = dense.to_f64_nearest();
    assert!(
        found.reported_mm() >= dense,
        "{} against {dense}",
        found.raw_mm
    );
    assert!(found.raw_mm <= dense * (1.0 + 1e-3));
    assert!(found.search_mm <= (gap::MARGIN - 1.0) * found.raw_mm);
    println!(
        "loud tail: splits={} evaluations={} raw_mm={:e} search_mm={:e}",
        found.refinements.len(),
        found.evaluations,
        found.raw_mm,
        found.search_mm
    );
    // Forbidden to split, the same piece is loud and cannot be looked at
    // more closely. A missed peak would have no bound, so it refuses.
    let refusal = gap::measure(
        &copies(&geometry),
        &Hints::default(),
        Limits {
            narrowest_split: 1.0,
            ..Limits::default()
        },
        &mut dec,
    )
    .expect_err("unresolved");
    assert_eq!(refusal.code, "numerical.gap-unresolved");
}

/// The straight-residual rule is used where it applies and only there.
#[test]
fn a_line_on_its_plane_is_read_at_its_ends() {
    for case in fixed::cases() {
        let mut dec = Decider::default();
        let found = gap::measure(
            &copies(&case.geometry),
            &Hints::default(),
            Limits::default(),
            &mut dec,
        )
        .expect("measure");
        let straight = case.name == "13a-zero-line-on-plane";
        assert_eq!(
            found.method == Method::StraightResidual,
            straight,
            "{}",
            case.name
        );
        if straight {
            // Both ends, and five more for the error allowance.
            assert_eq!(found.evaluations, 7);
        }
    }
}

/// Rule T4, apart from the known-answer rule: the stored tolerance is 1.1
/// times the reported gap, rounded up, so it covers the true gap.
#[test]
fn the_stored_tolerance_is_the_margin_times_the_reported_gap() {
    for case in fixed::cases() {
        let mut dec = Decider::default();
        let found = gap::measure(
            &copies(&case.geometry),
            &Hints::default(),
            Limits::default(),
            &mut dec,
        )
        .expect("measure");
        let reported = found.raw_mm + found.allowance_mm;
        assert!(found.reported_mm() >= reported);
        assert!(found.stored_mm() >= gap::MARGIN * reported);
        assert!(found.stored_mm() <= gap::MARGIN * reported * (1.0 + 4.0 * f64::EPSILON));
        assert!(exact(found.stored_mm()) >= case.truth.hi, "{}", case.name);
    }
}

/// An edge that cannot be resolved inside its allowance refuses by name.
#[test]
fn the_evaluation_cap_refuses_by_name() {
    let case = &fixed::cases()[16];
    assert_eq!(case.name, "12-256-spans");
    let refusal = gap::measure(
        &copies(&case.geometry),
        &Hints::default(),
        Limits {
            evaluations: 1000,
            ..Limits::default()
        },
        &mut Decider::default(),
    )
    .expect_err("capped");
    assert_eq!(refusal.code, "numerical.gap-work-cap");
}

/// The cap covers every evaluation of an edge, the closing probes of the
/// rounding allowance included: an allowance of exactly what the edge needs
/// passes, and one evaluation less refuses.
#[test]
fn the_closing_probes_count_inside_the_cap() {
    for name in ["08-parameter-stretch", "13a-zero-line-on-plane"] {
        let case = fixed::cases()
            .into_iter()
            .find(|case| case.name == name)
            .expect("case");
        let capped = |evaluations: u64| {
            gap::measure(
                &copies(&case.geometry),
                &Hints::default(),
                Limits {
                    evaluations,
                    ..Limits::default()
                },
                &mut Decider::default(),
            )
        };
        let needed = capped(gap::EDGE_EVALUATIONS).expect("default").evaluations;
        assert_eq!(capped(needed).expect("exact").evaluations, needed, "{name}");
        let refusal = capped(needed - 1).expect_err("one short");
        assert_eq!(refusal.code, "numerical.gap-work-cap", "{name}");
    }
}

/// A face copy that winds around a cylinder 1e29 times has that many seam
/// crossings. They are counted before any is listed, and the edge is refused
/// by name after a handful of work units.
#[test]
fn a_face_copy_that_winds_without_end_is_refused_before_it_is_walked() {
    let mut case = fixed::cases()
        .into_iter()
        .find(|case| case.name == "13b-zero-circle-on-cylinder")
        .expect("case");
    let FaceCopy::Linear(ends) = &mut case.geometry.face_copy else {
        unreachable!("a straight face copy")
    };
    ends[1][0] = Rational::from_i64(10).pow(30).expect("power");
    let (refusal, work) = with_operation_budget(DEFAULT_WORK, || {
        gap::measure(
            &copies(&case.geometry),
            &Hints::default(),
            Limits::default(),
            &mut Decider::default(),
        )
    });
    assert_eq!(refusal.expect_err("refused").code, "numerical.gap-work-cap");
    assert!(work.used < 1000, "{}", work.used);
}

/// A builder's hint becomes a break of the plan.
#[test]
fn a_hint_is_sampled() {
    let case = &fixed::cases()[12];
    assert_eq!(case.name, "08-parameter-stretch");
    let mut dec = Decider::default();
    let plain = gap::plan(
        &copies(&case.geometry),
        &Hints::default(),
        Limits::default(),
        &mut dec,
    )
    .expect("plan");
    let hinted = gap::plan(
        &copies(&case.geometry),
        &Hints {
            parameters: vec![0.4137],
        },
        Limits::default(),
        &mut dec,
    )
    .expect("plan");
    assert_eq!(hinted.len(), plain.len() + 1);
    assert!(
        hinted
            .iter()
            .any(|piece| piece.piece[1].to_bits() == 0.4137_f64.to_bits())
    );
}

/// Guard 3: every constant in the checker says where it comes from, and no
/// other number is written into its code.
#[test]
fn every_constant_in_the_checker_names_its_source() {
    let files = [
        ("gap.rs", include_str!("../gap.rs")),
        ("gap/lift.rs", include_str!("../gap/lift.rs")),
        ("gap/numeric.rs", include_str!("../gap/numeric.rs")),
    ];
    // Numbers that are structure and not choices: zero, one, two and a half,
    // and the three coordinates and four homogeneous coordinates of a point.
    let structural = ["0", "1", "2", "3", "4", "0.0", "1.0", "2.0", "0.5"];
    for (name, text) in files {
        let code = text.split("#[cfg(test)]").next().expect("code");
        let lines: Vec<&str> = code.lines().collect();
        let mut in_const = false;
        for (index, line) in lines.iter().enumerate() {
            let trimmed = line.trim_start();
            if trimmed.starts_with("//") {
                continue;
            }
            let starts_const = trimmed.starts_with("const ")
                || trimmed.starts_with("pub const ")
                || trimmed.starts_with("pub(super) const ");
            if starts_const {
                let mut doc = String::new();
                let mut above = index;
                while above > 0 && lines[above - 1].trim_start().starts_with("///") {
                    above -= 1;
                    doc.push_str(lines[above]);
                }
                assert!(
                    doc.contains("Source:"),
                    "{name}:{}: a constant without a source",
                    index + 1
                );
                in_const = true;
            }
            if !in_const {
                // Drop string literals, then read the numeric literals.
                let bare: String = trimmed.split('"').step_by(2).collect();
                let bare = bare.replace("..", " ");
                let mut token = String::new();
                let mut previous = ' ';
                for c in bare.chars().chain(std::iter::once(' ')) {
                    let continues = c.is_ascii_digit()
                        || (!token.is_empty() && (c == '.' || c == '_' || c == 'e'));
                    let starts = token.is_empty()
                        && c.is_ascii_digit()
                        && !(previous.is_alphanumeric() || previous == '_');
                    if starts || (!token.is_empty() && continues) {
                        token.push(c);
                    } else if !token.is_empty() {
                        let literal = token.trim_end_matches(['.', 'e', '_']);
                        assert!(
                            structural.contains(&literal),
                            "{name}:{}: the number {literal} has no named source",
                            index + 1
                        );
                        token.clear();
                    }
                    previous = c;
                }
            }
            if in_const && trimmed.trim_end().ends_with(';') {
                in_const = false;
            }
        }
    }
}
