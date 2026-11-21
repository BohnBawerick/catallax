//! Finite round-primitive clearance, and the rim piece past the mouth contact.
#[path = "../examples/m18/cases.rs"]
mod cases;
#[test]
fn finite_round_controls_preserve_full_parallel_answers_and_exhaustion() {
    for cap in [1_000_000, 100] {
        let mut previous = None;
        for workers in [Some(1), Some(2), Some(3), Some(4), Some(9), None] {
            let (records, spawned, charged) = cases::matrix(workers, cap);
            println!(
                "P5 cap={cap} requested={workers:?} spawned={spawned} charged={charged} work={:?} exhausted={:?}",
                records.iter().map(|r| r.work).collect::<Vec<_>>(),
                records.iter().map(|r| r.exhausted).collect::<Vec<_>>()
            );
            assert_eq!(
                spawned,
                workers
                    .unwrap_or_else(|| std::thread::available_parallelism().unwrap().get())
                    .min(8)
            );
            if cap == 1_000_000 {
                // Radius 1/2 in three placements keeps the elementary
                // torus; radii 1, 2 and 5/2 take the field's rim piece
                // (4.2a step 4); at 7/2 the drill's top reaches the fan.
                assert!(
                    records[..6]
                        .iter()
                        .all(|r| r.result.is_ok() && !r.exhausted)
                );
                assert_eq!(
                    records[6].result.as_ref().unwrap_err().code,
                    "fillet.unsupported-overflow"
                );
                assert_eq!(
                    records[7].result.as_ref().unwrap_err().code,
                    "fillet.non-positive-radius"
                );
            } else {
                assert!(records.iter().all(|r| r.exhausted && r.result.is_err()));
            }
            let semantic = (records, charged);
            if let Some(old) = &previous {
                assert_eq!(old, &semantic);
            } else {
                previous = Some(semantic);
            }
        }
    }
}
