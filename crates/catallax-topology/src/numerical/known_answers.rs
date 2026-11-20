//! Known-answer edges for an edge-gap checker.
//!
//! Each case is one edge with two copies: a 3D curve, and a face copy that is
//! lifted through a surface at the same parameter. The largest distance
//! between the copies is known exactly by construction, so a checker can be
//! judged against it. Nothing here measures a gap for admission; this module
//! is the oracle and holds no checker. A checker under test is handed an
//! [`EdgeGeometry`] and nothing else: names, truths and expectations stay here.
//!
//! The fourteen fixed edges copy recorded failures of sampled checkers in
//! OpenCASCADE, SMLib and solidcraft. `docs/research/cp2-grill/checker-failures.md`
//! lists the sources. The seeded edges in [`fresh`] are built by different
//! code, so a checker tuned to the fixed edges still meets hidden bumps.

use super::*;
use catallax_numbers::{Real, with_operation_budget};
use catallax_shapes::{BSplineCurve, BSplineSurface};

mod a_plus;
mod exact;
mod fixed;
mod fresh;
mod smooth;

/// The surface the face copy is lifted through.
#[derive(Clone, Debug)]
pub(super) enum CaseSurface {
    /// An analytic kernel surface.
    Kernel(Surface),
    /// A raw spline carrier. A kernel spline support needs a chart proof that
    /// a degree-25 fixture cannot supply, so the carrier is used directly.
    Carrier(BSplineSurface),
}

/// The copy of the edge in the face's parameters.
#[derive(Clone, Debug)]
pub(super) enum FaceCopy {
    /// A retained u/v curve in the 3D curve's parameter.
    Native(BSplineCurve),
    /// Straight in u/v between the two end values, as an ordinary trim.
    Linear(Box<[[Rational; 2]; 2]>),
}

/// The exact largest same-parameter gap, as a rational enclosure.
#[derive(Clone, Debug)]
pub(super) struct Truth {
    /// Lower end of the enclosure in millimetres.
    pub lo: Rational,
    /// Upper end of the enclosure in millimetres.
    pub hi: Rational,
    /// Edge parameter where the gap is largest, within the enclosure's precision.
    pub at: Rational,
}

impl Truth {
    pub(super) fn exact(value: Rational, at: Rational) -> Self {
        Self {
            lo: value.clone(),
            hi: value,
            at,
        }
    }
    pub(super) fn is_zero(&self) -> bool {
        self.hi.is_zero()
    }
}

/// Why two copies of a zero-gap control agree.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Zero {
    /// The face copy follows from the 3D curve by a rule: a line or circle on
    /// its own analytic surface, or a boundary view of the carrier. A checker
    /// may report an exact zero here.
    Structural,
    /// The copies agree as functions but are stored independently. Only a
    /// measurement can show it, so an exact zero would be a false claim.
    Measured,
}

/// Everything a checker is given: the two copies and the surface. It holds
/// no name, no answer and no expectation, so an adapter cannot read the oracle.
#[derive(Clone, Debug)]
pub(super) struct EdgeGeometry {
    /// The 3D copy.
    pub curve: Curve,
    /// Edge parameters in traversal order.
    pub parameters: [Rational; 2],
    /// The surface the face copy is lifted through.
    pub surface: CaseSurface,
    /// The face copy.
    pub face_copy: FaceCopy,
}

/// One edge with a known largest gap. Only [`EdgeCase::geometry`] reaches a
/// checker; the rest stays in the harness.
#[derive(Clone, Debug)]
pub(super) struct EdgeCase {
    /// Stable name used in receipts.
    pub name: String,
    /// The recorded failure or construction this case stands for.
    pub source: &'static str,
    /// What a checker sees.
    pub geometry: EdgeGeometry,
    /// The exact largest gap.
    pub truth: Truth,
    /// The parameter where the face copy leaves the surface domain, if it does.
    /// The truth is then measured against the end span's polynomial
    /// continuation, which is what a receiving system evaluates.
    pub domain_exit: Option<Rational>,
    /// Set on zero-gap controls.
    pub zero: Option<Zero>,
    /// Largest acceptable reported gap on a zero-gap control.
    pub zero_ceiling: Rational,
}

/// What a checker reports for one edge.
#[derive(Clone, Debug)]
pub(super) struct Reading {
    /// The largest gap the checker found, before any allowance or margin.
    pub raw_mm: Rational,
    /// The checker's own numerical error allowance for that value.
    pub allowance_mm: Rational,
    /// The tolerance the checker would store, if it stores one.
    pub stored_mm: Option<Rational>,
    /// Edge parameter of the reported gap.
    pub at: Option<Rational>,
    /// Whether the checker claims the copies agree exactly.
    pub exact_zero: bool,
    /// Where the checker saw the face copy leave the surface domain.
    pub domain_exit: Option<Rational>,
    /// Geometry evaluations spent.
    pub evaluations: u64,
}

impl Reading {
    /// The reported gap: the raw estimate plus the error allowance, before
    /// the storage margin. The known-answer rule reads this one quantity on
    /// both sides. A raw binary64 estimate is an attained value and sits at or
    /// below the true maximum by its rounding error, so the raw value alone
    /// cannot be required to reach the truth; the checker's own bound on that
    /// error is part of what it reports, and of what the margin multiplies.
    pub(super) fn reported_mm(&self) -> Rational {
        &self.raw_mm + &self.allowance_mm
    }
}

/// A checker under test. It is handed the geometry only.
pub(super) trait GapChecker {
    /// Name used in receipts.
    fn name(&self) -> &'static str;
    /// Measures one edge.
    fn measure(&self, geometry: &EdgeGeometry, dec: &mut Decider) -> Result<Reading, Refusal>;
}

/// One receipt row: a reading, its charged work and the oracle's verdict.
#[derive(Clone, Debug)]
pub(super) struct Row {
    pub name: String,
    pub truth: Truth,
    pub reading: Reading,
    pub work: u64,
    pub verdict: Result<(), String>,
}

/// The smallest value of rule T1's floor for binary64 geometry (Checkpoint 2
/// decisions, section 4.3): `eps(R) = max(1e-7 mm, 1e-11 x R)`.
fn binary64_floor_mm() -> Rational {
    Rational::from_ratio(1, 10_000_000).expect("nonzero")
}

/// The done-when rule of roadmap item 3.0A: the reported gap is at least the
/// true gap and at most 1.1 times it. Both sides read [`Reading::reported_mm`].
/// A truth known as an enclosure is taken at its upper end on the lower side
/// and at its lower end on the upper side, so the enclosure's width never
/// helps a checker. Nothing else is judged.
///
/// A zero-gap control has a true gap of zero, and 1.1 times zero is zero. No
/// binary64 measurement of a curved copy returns an exact zero, so a control
/// passes when its reported gap is within rule T1's floor for binary64
/// geometry: it must report about zero and raise no false alarm.
pub(super) fn known_answer_holds(case: &EdgeCase, reading: &Reading) -> Result<(), String> {
    let truth = &case.truth;
    let reported = reading.reported_mm();
    if case.zero.is_some() {
        if reported > binary64_floor_mm() {
            return Err("false alarm on a zero-gap control".into());
        }
        return Ok(());
    }
    if reported < truth.hi {
        return Err("reported gap is below the true gap".into());
    }
    let margin = Rational::from_ratio(11, 10).expect("nonzero");
    if reported > &truth.lo * &margin {
        return Err("reported gap is above 1.1 times the true gap".into());
    }
    Ok(())
}

/// Measures every case under its own work allowance. No case is skipped: a
/// caller that wants a subset passes a subset.
pub(super) fn run(checker: &dyn GapChecker, cases: &[EdgeCase]) -> Vec<Row> {
    cases
        .iter()
        .map(|case| {
            let (reading, work) = with_operation_budget(DEFAULT_WORK, || {
                checker.measure(&case.geometry, &mut Decider::default())
            });
            assert!(!work.exhausted, "{} ran out of work", case.name);
            // A refusal is a result: the row fails with the refusal's code.
            let (reading, verdict) = match reading {
                Ok(reading) => {
                    let verdict = known_answer_holds(case, &reading);
                    (reading, verdict)
                }
                Err(refusal) => (
                    Reading {
                        raw_mm: Rational::zero(),
                        allowance_mm: Rational::zero(),
                        stored_mm: None,
                        at: None,
                        exact_zero: false,
                        domain_exit: None,
                        evaluations: 0,
                    },
                    Err(format!("refused {}", refusal.code)),
                ),
            };
            Row {
                name: case.name.clone(),
                truth: case.truth.clone(),
                verdict,
                reading,
                work: work.used,
            }
        })
        .collect()
}

/// The acceptance run for a checker that claims the 3.0A done-when line:
/// every fixed case, by exact name and count, with no filter. A row that
/// fails or a case that is missing makes this an error listing each one.
pub(super) fn accept(checker: &dyn GapChecker) -> Result<Vec<Row>, String> {
    let rows = run(checker, &fixed::cases());
    all_required_rows_pass(&rows).map(|()| rows)
}

/// Checks rows against [`fixed::REQUIRED`]: the same names in the same order,
/// and every verdict a pass.
pub(super) fn all_required_rows_pass(rows: &[Row]) -> Result<(), String> {
    let names: Vec<&str> = rows.iter().map(|row| row.name.as_str()).collect();
    if names != fixed::REQUIRED {
        return Err(format!(
            "the required cases are {:?}, the run has {names:?}",
            fixed::REQUIRED
        ));
    }
    let failures: Vec<String> = rows
        .iter()
        .filter_map(|row| {
            row.verdict
                .as_ref()
                .err()
                .map(|reason| format!("{}: {reason}", row.name))
        })
        .collect();
    if failures.is_empty() {
        Ok(())
    } else {
        Err(failures.join("; "))
    }
}

/// A receipt: the seed of any generated cases and, per edge, the true gap, the
/// raw estimate, the error allowance, their sum (the reported gap), the stored
/// tolerance and the work.
pub(super) fn receipt(checker: &str, seed: Option<u64>, rows: &[Row]) -> String {
    let number = |value: &Rational| format!("{:.9e}", value.to_f64_nearest());
    let mut text = format!(
        "checker={checker} seed={}\n",
        seed.map_or_else(|| "none".into(), |seed| seed.to_string())
    );
    for row in rows {
        let ratio = if row.truth.is_zero() {
            "-".into()
        } else {
            format!(
                "{:.6}",
                row.reading.raw_mm.to_f64_nearest() / row.truth.lo.to_f64_nearest()
            )
        };
        text.push_str(&format!(
            "edge={} true_mm={} raw_mm={} allowance_mm={} reported_mm={} stored_mm={} \
             raw_over_true={ratio} at={} exact_zero={} domain_exit={} evaluations={} work={} \
             verdict={}\n",
            row.name,
            number(&row.truth.lo),
            number(&row.reading.raw_mm),
            number(&row.reading.allowance_mm),
            number(&row.reading.reported_mm()),
            row.reading.stored_mm.as_ref().map_or("-".into(), number),
            // Recorded, not judged.
            row.reading.at.as_ref().map_or("-".into(), number),
            row.reading.exact_zero,
            row.reading.domain_exit.as_ref().map_or("-".into(), number),
            row.reading.evaluations,
            row.work,
            row.verdict
                .as_ref()
                .map_or_else(|reason| format!("FAIL({reason})"), |()| "pass".into()),
        ));
    }
    text
}

/// Prints a receipt and, when `CATALLAX_GAP_RECEIPT` names a directory,
/// writes it there as well.
pub(super) fn publish(file: &str, text: &str) {
    println!("{text}");
    if let Some(directory) = std::env::var_os("CATALLAX_GAP_RECEIPT") {
        let path = std::path::Path::new(&directory).join(file);
        std::fs::write(&path, text).unwrap_or_else(|error| panic!("{}: {error}", path.display()));
    }
}

/// The oracle's own check, independent of any checker: the stated truth is
/// attained at the stated parameter and no sampled point exceeds it. Values
/// are 128-bit enclosures of the exact definitions.
pub(super) fn truth_is_attained_and_bounds(case: &EdgeCase, samples: i64) {
    let mut dec = Decider::default();
    let ((), work) = with_operation_budget(10 * DEFAULT_WORK, || {
        let truth = &case.truth;
        let [a, b] = &case.geometry.parameters;
        // Past a domain exit the truth is defined on the end span's polynomial
        // continuation, which the kernel evaluators refuse. Those cases state
        // their continuation in closed form and are checked up to the exit.
        let inside = |t: &Rational| {
            case.domain_exit
                .as_ref()
                .is_none_or(|exit| (t - a).abs() <= (exit - a).abs())
        };
        // The stated parameter can miss an irrational maximum by the width of
        // its bisection bracket, which moves the value by far less than this.
        let slack = &truth.hi * &exact::pow2(-60) + exact::pow2(-110);
        if inside(&truth.at) {
            let [lo, hi] = exact::gap_enclosure(&case.geometry, &truth.at, &mut dec)
                .unwrap_or_else(|refusal| panic!("{}: {refusal:?}", case.name));
            assert!(
                hi + slack.clone() >= truth.lo && lo <= &truth.hi + &slack,
                "{}: truth is not attained at its parameter",
                case.name
            );
        }
        let step = (b - a) / Rational::from_i64(samples);
        let near = (b - a).abs();
        let mut parameters: Vec<Rational> = (0..=samples)
            .map(|i| a + &(&step * &Rational::from_i64(i)))
            .collect();
        for k in 4..40 {
            for side in [-1, 1] {
                let t = &truth.at + &(&near * &exact::pow2(-k) * Rational::from_i64(side));
                if &t >= a.min(b) && &t <= a.max(b) {
                    parameters.push(t);
                }
            }
        }
        for t in parameters.into_iter().filter(inside) {
            let [lo, _] = exact::gap_enclosure(&case.geometry, &t, &mut dec)
                .unwrap_or_else(|refusal| panic!("{}: {refusal:?}", case.name));
            assert!(
                lo <= &truth.hi + &slack,
                "{}: a sampled gap exceeds the stated truth",
                case.name
            );
        }
    });
    assert!(!work.exhausted, "{}: oracle ran out of work", case.name);
}

#[cfg(test)]
mod rule_tests {
    use super::*;

    fn reading(raw: &str, allowance: &str) -> Reading {
        Reading {
            raw_mm: exact::q(raw),
            allowance_mm: exact::q(allowance),
            stored_mm: None,
            at: None,
            exact_zero: false,
            domain_exit: None,
            evaluations: 0,
        }
    }

    fn with_truth(lo: &str, hi: &str) -> EdgeCase {
        let mut case = fixed::cases().remove(0);
        case.truth = Truth {
            lo: exact::q(lo),
            hi: exact::q(hi),
            at: Rational::zero(),
        };
        case
    }

    #[test]
    fn both_sides_of_the_rule_read_the_reported_gap() {
        let case = with_truth("1", "1");
        let holds = |raw, allowance| known_answer_holds(&case, &reading(raw, allowance));
        assert!(holds("1", "0").is_ok());
        assert!(holds("1.1", "0").is_ok());
        // A raw estimate short of the truth fails unless the checker's own
        // allowance closes the difference.
        assert!(holds("0.9999", "0").is_err());
        assert!(holds("0.9999", "0.00009").is_err());
        assert!(holds("0.9999", "0.0001").is_ok());
        // The allowance counts on the upper side as well.
        assert!(holds("1.1", "0.0001").is_err());
        assert!(holds("1.0995", "0.001").is_err());
        assert!(holds("1.11", "0").is_err());
        // The rule reads the sum and nothing else: a large allowance that
        // keeps the sum inside the band passes.
        assert!(holds("0.998", "0.002").is_ok());
        assert!(holds("1", "0.0011").is_ok());
    }

    #[test]
    fn an_enclosed_truth_never_widens_the_rule() {
        let case = with_truth("1", "1.01");
        let holds = |raw| known_answer_holds(&case, &reading(raw, "0"));
        // The lower side reads the enclosure's upper end.
        assert!(holds("1.005").is_err());
        assert!(holds("1.01").is_ok());
        // The upper side reads 1.1 times its lower end.
        assert!(holds("1.1").is_ok());
        assert!(holds("1.105").is_err());
    }

    #[test]
    fn a_zero_control_fails_on_a_false_alarm() {
        let cases = fixed::cases();
        for control in cases.iter().filter(|case| case.zero.is_some()) {
            let holds = |raw, allowance| known_answer_holds(control, &reading(raw, allowance));
            assert!(holds("0", "0").is_ok());
            assert!(holds("1e-17", "1e-15").is_ok());
            assert!(holds("0", "1e-7").is_ok());
            assert!(holds("1e-7", "1e-9").is_err());
            assert!(holds("1e-6", "0").is_err());
        }
    }

    #[test]
    fn acceptance_needs_every_required_case_by_name() {
        struct Blind;
        impl GapChecker for Blind {
            fn name(&self) -> &'static str {
                "blind"
            }
            fn measure(&self, _: &EdgeGeometry, _: &mut Decider) -> Result<Reading, Refusal> {
                Ok(reading("0", "0"))
            }
        }
        // A checker that reports nothing passes the zero-gap controls only.
        let rows = run(&Blind, &fixed::cases());
        let passing: Vec<_> = rows.iter().filter(|row| row.verdict.is_ok()).collect();
        assert_eq!(passing.len(), 4);
        assert!(accept(&Blind).is_err());
        // Passing rows alone are not acceptance: a missing case is an error.
        let mut all_pass = rows.clone();
        for row in &mut all_pass {
            row.verdict = Ok(());
        }
        assert!(all_required_rows_pass(&all_pass).is_ok());
        all_pass.remove(3);
        assert!(all_required_rows_pass(&all_pass).is_err());
        all_pass.retain(|row| row.name.starts_with("13"));
        assert!(all_required_rows_pass(&all_pass).is_err());
    }
}
