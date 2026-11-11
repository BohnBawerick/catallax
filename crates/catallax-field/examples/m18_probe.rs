//! Public P5 field coverage and real parallel validation; no B-rep or oracle.
#[path = "m18/cases.rs"]
mod cases;
fn main() {
    for cap in [1_000_000, 100] {
        let mut previous = None;
        for workers in [Some(1), Some(2), Some(4), None] {
            let (records, spawned, charged) = cases::matrix(workers, cap);
            println!("requested={workers:?} spawned={spawned} job_cap={cap} charged={charged}");
            for r in &records {
                println!(
                    "radius={} placement={} outcome={} work={} exhausted={}",
                    r.radius,
                    r.placement,
                    r.result
                        .as_ref()
                        .map_or_else(|e| e.code, |_| "built-and-queried"),
                    r.work,
                    r.exhausted
                );
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
