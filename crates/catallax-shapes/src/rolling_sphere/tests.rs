use super::*;
use catallax_numbers::with_operation_budget;
use std::time::{Duration, Instant};

fn point(x: i64, y: i64, z: i64) -> Point3 {
    Point3::new(Real::from_i64(x), Real::from_i64(y), Real::from_i64(z))
}

fn definition(branch: i64, dec: &mut Decider) -> RollingSphereDefinition {
    let axis = Vector3::new(Real::zero(), Real::zero(), Real::one());
    let reference = Vector3::new(Real::one(), Real::zero(), Real::zero());
    RollingSphereDefinition::new(
        &Sphere {
            center: point(0, 0, 0),
            radius: Real::from_i64(5),
            axis: axis.clone(),
            ref_dir: reference.clone(),
        },
        &Cylinder {
            origin: point(1, 0, -10),
            axis,
            radius_sq: Real::from_i64(4),
            ref_dir: reference,
        },
        &Real::from(Rational::from_ratio(1, 2).unwrap()),
        [BallSide::Outside, BallSide::Outside],
        &point(3, 0, 4 * branch),
        dec,
    )
    .unwrap()
}

fn trial(branch: i64, reuse: bool) -> ((Vec<Node>, Rational), u64) {
    let (result, work) = with_operation_budget(1_000_000, || {
        let mut dec = Decider::default();
        let def = definition(branch, &mut dec);
        let k = Constants::new(&def, BITS, &dec).unwrap();
        let mut samples = FitSamples::default();
        let coarse = produce(&k, FIRST_SPANS * 2, BITS, &mut samples, &mut dec).unwrap();
        let _ = recompute_bound(&k, &coarse, BITS, &mut samples, None, &mut dec);
        if !reuse {
            samples = FitSamples::default();
        }
        let nodes = produce(&k, FIRST_SPANS * 4, BITS, &mut samples, &mut dec).unwrap();
        let bound = recompute_bound(&k, &nodes, BITS, &mut samples, None, &mut dec)
            .unwrap()
            .unwrap();
        assert!(samples.nodes.len() <= 2 * LAST_SPANS);
        assert!(samples.rows.len() <= 2 * (LAST_SPANS + 1));
        (nodes, bound)
    });
    assert!(!work.exhausted);
    (result, work.used)
}

#[test]
fn refinement_reuses_identical_nodes_and_bounds_on_both_rings() {
    for branch in [-1, 1] {
        let ((cold, cold_bound), cold_work) = trial(branch, false);
        let ((reused, reused_bound), reused_work) = trial(branch, true);
        assert_eq!(cold_bound, reused_bound);
        assert_eq!(cold.len(), reused.len());
        for (a, b) in cold.iter().zip(&reused) {
            assert_eq!(a.u, b.u);
            assert_eq!(a.values, b.values);
            assert_eq!(a.slopes, b.slopes);
            assert_eq!(a.arriving, b.arriving);
        }
        assert!(reused_work <= cold_work, "{reused_work} > {cold_work}");
    }
}

#[test]
fn early_rejection_preserves_the_full_bound_at_the_threshold() {
    let ((), work) = with_operation_budget(1_000_000, || {
        let mut dec = Decider::default();
        let def = definition(1, &mut dec);
        let k = Constants::new(&def, BITS, &dec).unwrap();
        let mut samples = FitSamples::default();
        let nodes = produce(&k, FIRST_SPANS * 4, BITS, &mut samples, &mut dec).unwrap();
        let upper = recompute_bound(&k, &nodes, BITS, &mut samples, None, &mut dec)
            .unwrap()
            .unwrap();
        assert!(upper > Rational::zero());
        let mut early = FitSamples::default();
        assert!(
            recompute_bound(
                &k,
                &nodes,
                BITS,
                &mut early,
                Some(&Rational::zero()),
                &mut dec
            )
            .unwrap()
            .is_none()
        );
        assert!(early.rows.len() < samples.rows.len());
        assert_eq!(
            recompute_bound(&k, &nodes, BITS, &mut early, Some(&upper), &mut dec).unwrap(),
            Some(upper)
        );
    });
    assert!(!work.exhausted);
}

#[test]
fn cached_samples_still_observe_deadlines() {
    let mut dec = Decider::default();
    let def = definition(1, &mut dec);
    let k = Constants::new(&def, BITS, &dec).unwrap();
    let mut samples = FitSamples::default();
    let nodes = produce(&k, FIRST_SPANS, BITS, &mut samples, &mut dec).unwrap();
    dec.set_deadline(Some(Instant::now() - Duration::from_secs(1)));
    assert_eq!(
        produce(&k, FIRST_SPANS, BITS, &mut samples, &mut dec)
            .unwrap_err()
            .code,
        "numbers.deadline"
    );
    assert_eq!(
        recompute_bound(&k, &nodes, BITS, &mut samples, None, &mut dec)
            .unwrap_err()
            .code,
        "numbers.deadline"
    );
}

#[test]
fn retained_homogeneous_controls_match_the_world_carrier_with_less_work() {
    use catallax_numbers::{Outcome, reuse_operation_enclosures};
    for branch in [-1, 1] {
        let mut dec = Decider::default();
        let def = definition(branch, &mut dec);
        let k = Constants::new(&def, BITS, &dec).unwrap();
        let nodes = produce(&k, FIRST_SPANS, BITS, &mut FitSamples::default(), &mut dec).unwrap();
        let (carrier, homogeneous) = build(&def, &nodes, &mut dec).unwrap();
        let (direct, direct_work) = with_operation_budget(1_000_000, || {
            reuse_operation_enclosures();
            Net::new(&homogeneous, &dec).unwrap()
        });
        let (round_trip, round_trip_work) = with_operation_budget(1_000_000, || {
            reuse_operation_enclosures();
            let frame = def.frame(BITS, &dec).unwrap();
            carrier
                .controls()
                .iter()
                .zip(carrier.weights())
                .map(|(point, weight)| {
                    let w = dec.enclose_at(weight, BITS).unwrap();
                    let [x, y, z] = frame.local(&math::point(point, BITS, &dec).unwrap());
                    [&x * &w, &y * &w, &z * &w, w]
                })
                .collect::<Vec<_>>()
        });
        assert!(!direct_work.exhausted && !round_trip_work.exhausted);
        assert!(direct_work.used < round_trip_work.used);
        for (a, b) in direct.controls.iter().zip(&round_trip) {
            for (a, b) in a.iter().zip(b) {
                assert!(a.lo() <= b.hi() && b.lo() <= a.hi());
                assert!(a.hi() - a.lo() <= b.hi() - b.lo());
            }
        }
        for index in [0, homogeneous.len() / 2, homogeneous.len() - 1] {
            let local = def.local(&carrier.controls()[index]);
            for (axis, value) in local.iter().enumerate() {
                let residual = &(value * &carrier.weights()[index]) - &homogeneous[index][axis];
                assert!(matches!(dec.sign(&residual), Outcome::ProvenCoincident(_)));
            }
        }
        dec.set_deadline(Some(Instant::now() - Duration::from_secs(1)));
        assert_eq!(
            Net::new(&homogeneous, &dec).unwrap_err().code,
            "numbers.deadline"
        );
    }
}
