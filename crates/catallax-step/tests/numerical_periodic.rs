//! One shared cylinder seam carries two explicit periodic pcurves.
#[path = "../../catallax-topology/tests/common/periodic.rs"]
mod fixture;
use catallax_numbers::{Ball, Decider, Rational};
use catallax_step::numerical::{self as output, NumericalStepOptions};
use catallax_topology::{Sense, numerical as n};

#[test]
fn a_complete_bore_emits_one_seam_curve_with_two_uses() {
    let body = n::build(
        fixture::bore(),
        n::NumericalOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let file = output::build(
        &body,
        &NumericalStepOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let text = std::str::from_utf8(file.bytes()).unwrap();
    assert_eq!(text.matches("=SEAM_CURVE(").count(), 1);
    assert_eq!(file.report().counts, [7, 15, 10]);
    file.check_written(file.bytes(), &mut Decider::default())
        .unwrap();
    if let Ok(path) = std::env::var("CATALLAX_ANALYTIC_EVIDENCE") {
        std::fs::write(
            std::path::Path::new(&path).join("periodic-bore.step"),
            file.bytes(),
        )
        .unwrap();
    }
}

#[test]
fn equivalent_closed_rim_directions_and_periods_emit_the_same_step() {
    let period = Ball::pi(128).mid() * Rational::from_i64(2);
    let export = |definition| {
        let body = n::build(
            definition,
            n::NumericalOptions::default(),
            &mut Decider::default(),
        )
        .unwrap();
        output::build(
            &body,
            &NumericalStepOptions::default(),
            &mut Decider::default(),
        )
        .unwrap()
    };
    // The negative direction is the real 00003390 failure. Period shifts also
    // exercise the natural range which a STEP reader reconstructs from vertices.
    for phased in [false, true] {
        let mut baseline = fixture::bore();
        let index = baseline
            .edges
            .iter()
            .position(|edge| edge.key == "rim-0")
            .unwrap();
        let start = if phased {
            let catallax_shapes::Curve::Circle(circle) = &mut baseline.edges[index].curve else {
                unreachable!()
            };
            circle.ref_dir = catallax_shapes::RVector3::new(0.into(), 1.into(), 0.into()).to_real();
            -(&period / &Rational::from_i64(4))
        } else {
            Rational::zero()
        };
        baseline.edges[index].parameters = [start.clone(), &start + &period];
        let expected = export(baseline.clone());
        if phased && let Ok(path) = std::env::var("CATALLAX_ANALYTIC_EVIDENCE") {
            std::fs::write(
                std::path::Path::new(&path).join("periodic-bore-phased.step"),
                expected.bytes(),
            )
            .unwrap();
        }
        for (direction, periods) in [(1, -1), (1, 1), (-1, -1), (-1, 0), (-1, 1)] {
            let mut definition = baseline.clone();
            let phase = &start + &(&period * &Rational::from_i64(periods));
            definition.edges[index].parameters = [
                phase.clone(),
                &phase + &(&period * &Rational::from_i64(direction)),
            ];
            if direction < 0 {
                for trim in definition
                    .faces
                    .iter_mut()
                    .flat_map(|face| face.loops.iter_mut().flatten())
                {
                    if trim.edge == index {
                        trim.sense = match trim.sense {
                            Sense::Same => Sense::Reversed,
                            Sense::Reversed => Sense::Same,
                        };
                    }
                }
            }
            let file = export(definition);
            assert!(
                file.bytes() == expected.bytes(),
                "phased={phased}, direction={direction}, periods={periods}"
            );
        }
    }
}
