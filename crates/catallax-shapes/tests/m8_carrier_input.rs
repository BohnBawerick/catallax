//! The retained P5 carrier posed to M10's existing checked input routes.
//! Neither route may round its algebraic controls or enlarge the allowance.

use catallax_numbers::{Decider, Rational, Real, Refusal, with_operation_budget};
use catallax_shapes::{
    Point3, ROLLING_ENVELOPE_TARGET_NM, RationalAtlas, RationalPatch, RollingEnvelope,
    RollingEnvelopeDefinition, Vector3,
};
use std::num::NonZeroUsize;

const ALLOWANCE: u64 = 20_000_000;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Input {
    Patch,
    Atlas,
}

#[derive(Debug, PartialEq, Eq)]
struct Record {
    input: Input,
    reached_input: bool,
    nonrational_coefficients: usize,
    answer: Result<(), (String, String)>,
    used: u64,
    exhausted: bool,
}

fn vector(x: i64, y: i64, z: i64) -> Vector3 {
    Vector3::new(x.into(), y.into(), z.into())
}

fn run(input: Input) -> Record {
    let mut reached_input = false;
    let mut nonrational_coefficients = 0;
    let (answer, work) = with_operation_budget(ALLOWANCE, || -> Result<(), Refusal> {
        let mut dec = Decider::default();
        let definition = RollingEnvelopeDefinition::new(
            Point3::new(0.into(), 5.into(), 0.into()),
            vector(0, 1, 0),
            vector(1, 0, 0),
            12.into(),
            4.into(),
            17.into(),
            2.into(),
            [Real::one(), Real::one()],
            &mut dec,
        )?;
        let envelope = RollingEnvelope::fit(definition, &mut dec)?;
        let target = Rational::from_ratio(ROLLING_ENVELOPE_TARGET_NM, 1_000_000).unwrap();
        assert!(envelope.bound().upper <= target);
        let carrier = envelope.carrier();
        assert_eq!(carrier.degrees(), [8, 2]);
        assert_eq!(carrier.control_count(), [257, 3]);
        nonrational_coefficients = carrier
            .controls()
            .iter()
            .flat_map(|point| [&point.x, &point.y, &point.z])
            .chain(carrier.weights())
            .filter(|value| value.rational_ref().is_none())
            .count();
        reached_input = true;
        match input {
            Input::Patch => RationalPatch::new(carrier.clone(), &mut dec).map(|_| ()),
            Input::Atlas => {
                RationalAtlas::partition(carrier.clone(), [false, false], &mut dec).map(|_| ())
            }
        }
    });
    Record {
        input,
        reached_input,
        nonrational_coefficients,
        answer: answer.map_err(|refusal| (refusal.code.to_owned(), refusal.detail)),
        used: work.used,
        exhausted: work.exhausted,
    }
}

#[test]
#[ignore = "M8 opt-in: retained fitted carrier at M10's input checks"]
fn retained_p5_carrier_posed_to_m10_input_checks() {
    let mut first = None;
    for threads in [Some(1), Some(2), Some(4), None] {
        catallax_numbers::set_threads(threads.and_then(NonZeroUsize::new));
        let records = [Input::Patch, Input::Atlas].map(run);
        println!("M8 threads={threads:?} {records:?}");
        for (record, used) in records.iter().zip([19_999_970, 19_999_995]) {
            assert!(record.reached_input);
            assert_eq!(record.nonrational_coefficients, 844);
            assert_eq!(record.used, used);
            assert!(record.exhausted);
            assert_eq!(
                record.answer.as_ref().unwrap_err().0,
                "numbers.operation-budget"
            );
        }
        match &first {
            Some(first) => assert_eq!(first, &records),
            None => first = Some(records),
        }
    }
    catallax_numbers::set_threads(None);
}
