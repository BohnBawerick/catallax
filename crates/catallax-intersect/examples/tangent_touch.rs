//! R2's bounded plane experiment. Run explicitly with `cargo run --release
//! -p catallax-intersect --example tangent_touch`. It never builds a solid.

#[path = "tangent_touch/fixtures.rs"]
mod fixtures;
#[path = "tangent_touch/solver.rs"]
mod solver;

use catallax_numbers::{Decider, Refusal, with_operation_budget};
use fixtures::{CASES, COPIES, build};
use solver::{Kind, contacts};

const ALLOWANCE: u64 = 20_000_000;

fn request(case: &str, copy: &str) -> (String, u64, bool) {
    let (result, work) = with_operation_budget(ALLOWANCE, || -> Result<String, Refusal> {
        let mut dec = Decider::default();
        let (patch, plane) = build(case, copy, &mut dec)?;
        let found = contacts(&patch, &plane, &mut dec)?;
        let expected: &[Kind] = match case {
            "bowl" | "off-centre" | "steep-off-centre" | "hill" | "rational-bowl"
            | "touch-and-cut" => &[Kind::IsolatedTouch],
            "several" | "algebraic" => &[Kind::IsolatedTouch, Kind::IsolatedTouch],
            "saddle" => &[Kind::TwoBranchesCrossing],
            "flat" => &[Kind::DegenerateContact],
            "near-miss" | "tiny-crossing" | "transverse" => &[],
            _ => panic!("{case} must refuse"),
        };
        assert_eq!(
            found.iter().map(|c| c.kind).collect::<Vec<_>>(),
            expected,
            "{case}/{copy}"
        );
        let mut answer = format!("{:?}", expected);
        for contact in found {
            if let Err(refusal) = contact.kind.local_outcome() {
                answer.push_str(&format!(" refused:{}", refusal.code));
            }
            for parameter in contact.parameters {
                let bound = dec.enclose_budgeted_at(&parameter, 48)?;
                answer.push_str(&format!(" [{},{}]", bound.lo(), bound.hi()));
            }
        }
        Ok(answer)
    });
    (
        result.unwrap_or_else(|error| format!("refused:{}", error.code)),
        work.used,
        work.exhausted,
    )
}

fn round(workers: usize) -> Vec<(usize, String, u64, bool)> {
    let mut rows = std::thread::scope(|scope| {
        let jobs = (0..workers)
            .map(|worker| {
                scope.spawn(move || {
                    (worker..CASES.len() * COPIES.len())
                        .step_by(workers)
                        .map(|index| {
                            let (answer, work, exhausted) =
                                request(CASES[index / COPIES.len()], COPIES[index % COPIES.len()]);
                            (index, answer, work, exhausted)
                        })
                        .collect::<Vec<_>>()
                })
            })
            .collect::<Vec<_>>();
        jobs.into_iter()
            .flat_map(|job| job.join().unwrap())
            .collect::<Vec<_>>()
    });
    rows.sort_by_key(|r| r.0);
    rows
}

fn main() {
    let mut baseline = None;
    let mut copies_equal = true;
    for workers in [1, 2, 4] {
        let rows = round(workers);
        for (index, answer, work, exhausted) in &rows {
            println!(
                "workers={workers}\t{}\t{}\t{work}/{ALLOWANCE}\texhausted={exhausted}\t{answer}",
                CASES[index / COPIES.len()],
                COPIES[index % COPIES.len()]
            );
        }
        for copies in rows.as_chunks::<{ COPIES.len() }>().0 {
            let same = copies
                .iter()
                .all(|r| r.1 == copies[0].1 && r.3 == copies[0].3);
            copies_equal &= same;
            println!(
                "copies workers={workers}\t{}\tequal={same}",
                CASES[copies[0].0 / COPIES.len()]
            );
        }
        if let Some(baseline) = &baseline {
            assert_eq!(&rows, baseline);
        } else {
            baseline = Some(rows);
        }
    }
    println!("cylinder=refused:prototype.tangent-cylinder-unmeasured");
    println!("cone=refused:prototype.tangent-cone-unsupported");
    println!("sphere=refused:prototype.tangent-sphere-out-of-scope");
    println!("torus=refused:prototype.tangent-torus-unsupported");
    assert!(copies_equal, "transformed copies changed a contact answer");
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::Sign;

    #[test]
    fn exact_nodes_touches_and_degeneracies() {
        for (case, expected) in [
            ("bowl", "[IsolatedTouch]"),
            ("saddle", "[TwoBranchesCrossing]"),
            ("flat", "[DegenerateContact]"),
            ("algebraic", "[IsolatedTouch, IsolatedTouch]"),
            ("near-miss", "[]"),
            ("tiny-crossing", "[]"),
            ("flush", "refused:prototype.tangent-flush-span"),
            ("line", "refused:prototype.tangent-stationary-component"),
            ("boundary", "refused:prototype.tangent-boundary-contact"),
            ("crease", "refused:shapes.rational-patch-unproven"),
        ] {
            let (answer, _, exhausted) = request(case, "original");
            assert!(!exhausted, "{case}");
            assert!(answer.starts_with(expected), "{case}: {answer}");
        }
    }

    #[test]
    fn exhausted_allowance_never_returns_contacts() {
        let (result, work) = with_operation_budget(100, || {
            let mut dec = Decider::default();
            let (patch, plane) = build("bowl", "original", &mut dec)?;
            contacts(&patch, &plane, &mut dec)
        });
        assert!(work.exhausted);
        assert!(result.is_err());
    }

    #[test]
    fn isolated_points_really_lie_on_the_plane() {
        let mut dec = Decider::default();
        for case in [
            "bowl",
            "off-centre",
            "several",
            "algebraic",
            "rational-bowl",
        ] {
            let (patch, plane) = build(case, "combined", &mut dec).unwrap();
            for contact in contacts(&patch, &plane, &mut dec).unwrap() {
                let point = patch.point_at(&contact.parameters, &mut dec).unwrap();
                assert_eq!(
                    solver::sign(&point.minus(&plane.origin).dot(&plane.normal), &mut dec).unwrap(),
                    Sign::Zero
                );
            }
        }
    }
}
