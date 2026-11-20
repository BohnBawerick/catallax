//! Known-answer edges for an edge-gap checker.
//!
//! Each case is one edge with two copies: a 3D curve, and a face copy that is
//! lifted through a surface at the same parameter. The largest distance
//! between the copies is known exactly by construction, so a checker can be
//! judged against it. Nothing here measures a gap for admission; this module
//! is the oracle and holds no checker.
//!
//! The fourteen fixed edges copy recorded failures of sampled checkers in
//! OpenCASCADE, SMLib and solidcraft. `docs/research/cp2-grill/checker-failures.md`
//! lists the sources. The seeded edges in [`fresh`] are built by different
//! code, so a checker tuned to the fixed edges still meets hidden bumps.

use super::*;
use catallax_numbers::{Real, with_operation_budget};
use catallax_shapes::{BSplineCurve, BSplineSurface};

mod exact;
mod fixed;
mod fresh;
mod legacy;

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

/// One edge with a known largest gap.
#[derive(Clone, Debug)]
pub(super) struct EdgeCase {
    /// Stable name used in receipts.
    pub name: String,
    /// The recorded failure or construction this case stands for.
    pub source: &'static str,
    /// The 3D copy.
    pub curve: Curve,
    /// Edge parameters in traversal order.
    pub parameters: [Rational; 2],
    /// The surface the face copy is lifted through.
    pub surface: CaseSurface,
    /// The face copy.
    pub face_copy: FaceCopy,
    /// The exact largest gap.
    pub truth: Truth,
    /// The parameter where the face copy leaves the surface domain, if it does.
    /// The truth is then measured against the end span's polynomial
    /// continuation, which is what a receiving system evaluates.
    pub domain_exit: Option<Rational>,
    /// Set on zero-gap controls.
    pub zero: Option<Zero>,
    /// Largest acceptable raw value plus allowance on a zero-gap control.
    pub zero_ceiling: Rational,
}

/// What a checker reports for one edge.
#[derive(Clone, Debug)]
pub(super) struct Reading {
    /// The largest gap the checker found, before any margin.
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

/// A checker under test. Implementations see only the geometry of the case.
pub(super) trait GapChecker {
    /// Name used in receipts.
    fn name(&self) -> &'static str;
    /// Whether this checker has a rule for the case's kind of face copy.
    fn applies(&self, case: &EdgeCase) -> bool;
    /// Measures one edge.
    fn measure(&self, case: &EdgeCase, dec: &mut Decider) -> Result<Reading, Refusal>;
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

/// Share of the true gap that a checker's error allowance may take on an edge
/// with a nonzero answer. The stored margin is ten percent; an allowance above
/// one hundredth of that margin would let the allowance, not the measurement,
/// cover a missed maximum.
fn allowance_share() -> Rational {
    Rational::from_ratio(1, 1000).expect("nonzero")
}

/// The done-when rule of roadmap item 3.0A, read on the raw estimate:
/// the reported gap is at least the true gap and at most 1.1 times it.
pub(super) fn known_answer_holds(case: &EdgeCase, reading: &Reading) -> Result<(), String> {
    let truth = &case.truth;
    if reading.raw_mm.sign() == Sign::Negative || reading.allowance_mm.sign() == Sign::Negative {
        return Err("negative gap or allowance".into());
    }
    match (&case.domain_exit, &reading.domain_exit) {
        (Some(_), None) => return Err("face copy leaves the surface domain unreported".into()),
        (None, Some(_)) => return Err("reported a domain exit that does not exist".into()),
        (Some(expected), Some(found)) => {
            let span = (&case.parameters[1] - &case.parameters[0]).abs();
            if (expected - found).abs() > &span * &exact::pow2(-20) {
                return Err("domain exit reported at the wrong parameter".into());
            }
        }
        (None, None) => (),
    }
    if let Some(zero) = case.zero {
        if reading.exact_zero && zero == Zero::Measured {
            return Err("measured coincidence reported as exact".into());
        }
        if reading.exact_zero && !reading.raw_mm.is_zero() {
            return Err("exact zero claimed with a nonzero gap".into());
        }
        if &reading.raw_mm + &reading.allowance_mm > case.zero_ceiling {
            return Err("false alarm on a zero-gap control".into());
        }
        return Ok(());
    }
    if reading.exact_zero {
        return Err("exact zero claimed on an edge with a gap".into());
    }
    if reading.allowance_mm > &truth.lo * &allowance_share() {
        return Err("error allowance is too large to judge the measurement".into());
    }
    if &reading.raw_mm + &reading.allowance_mm < truth.lo {
        return Err("reported gap is below the true gap".into());
    }
    let margin = Rational::from_ratio(11, 10).expect("nonzero");
    if reading.raw_mm > &truth.hi * &margin {
        return Err("reported gap is above 1.1 times the true gap".into());
    }
    Ok(())
}

/// Measures every applicable case under its own work allowance.
pub(super) fn run(checker: &dyn GapChecker, cases: &[EdgeCase]) -> Vec<Row> {
    cases
        .iter()
        .filter(|case| checker.applies(case))
        .map(|case| {
            let (reading, work) = with_operation_budget(DEFAULT_WORK, || {
                checker.measure(case, &mut Decider::default())
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

/// A receipt: the seed of any generated cases and, per edge, the true gap, the
/// raw estimate, the error allowance, the stored tolerance and the work.
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
            "edge={} true_mm={} raw_mm={} allowance_mm={} stored_mm={} found_over_true={ratio} \
             evaluations={} work={} verdict={}\n",
            row.name,
            number(&row.truth.lo),
            number(&row.reading.raw_mm),
            number(&row.reading.allowance_mm),
            row.reading.stored_mm.as_ref().map_or("-".into(), number),
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
        let [a, b] = &case.parameters;
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
            let [lo, hi] = exact::gap_enclosure(case, &truth.at, &mut dec)
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
            let [lo, _] = exact::gap_enclosure(case, &t, &mut dec)
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
