#![allow(missing_docs)]
#[path = "../tests/common/mod.rs"]
mod common;
use common::boolean_cases::*;
use common::*;
use std::hint::black_box;
use std::time::Instant;

fn record<T>(
    case: &str,
    property: &str,
    digits: i32,
    run: impl FnOnce(&mut Decider) -> Answered<T>,
    widths: impl FnOnce(&T) -> Vec<Rational>,
) {
    let mut dec = Decider::default();
    let start = Instant::now();
    let result = black_box(run(&mut dec));
    let nanos = start.elapsed().as_nanos();
    let stats = dec.stats();
    match result {
        Ok(answer) => println!(
            "case={case} property={property} relative=1e-{digits} wall_ns={nanos} decider={stats:?} widths={:?} marks={:?}",
            widths(&answer.value),
            answer.marks
        ),
        Err(refusal) => println!(
            "case={case} property={property} relative=1e-{digits} wall_ns={nanos} decider={stats:?} refused={}",
            refusal.code
        ),
    }
}
fn mass_widths(m: &MassProperties) -> Vec<Rational> {
    let mut w = vec![width(&m.volume), width(&m.area)];
    if let Some(c) = &m.centroid {
        w.extend(c.iter().map(width));
    }
    if let Some(tensor) = &m.inertia {
        for (i, row) in tensor.iter().enumerate() {
            w.extend(row[i..].iter().map(width));
        }
    }
    w
}
fn main() {
    println!(
        "measurement_ball_work=unavailable measurement_max_precision=unavailable until an observation API exists; Decider counters exclude measurement's own Ball arithmetic"
    );
    if std::env::var_os("CATALLAX_MEASURE_RUN_BENCH").as_deref() != Some(std::ffi::OsStr::new("1"))
    {
        println!(
            "Disabled by default. Set CATALLAX_MEASURE_RUN_BENCH=1 to measure piece 6; add CATALLAX_MEASURE_BOOLEAN_BENCH=1 for piece 5 results."
        );
        return;
    }
    // No timing assertions, iteration deadline, or speed-based exit code.
    // Centroid and inertia time are included in mass_properties.
    for case in CASES {
        let f = build(case, &r(1), &Pose::identity());
        let name = format!("{case:?}");
        for digits in [6, 12, 18, 24, 30] {
            let a = accuracy(digits);
            record(
                &name,
                "volume",
                digits,
                |d| measure::volume(&f.model, f.body, &a, d),
                |e| vec![width(e)],
            );
            record(
                &name,
                "area",
                digits,
                |d| measure::area(&f.model, f.body, &a, d),
                |e| vec![width(e)],
            );
            record(
                &name,
                "mass_properties",
                digits,
                |d| measure::mass_properties(&f.model, f.body, &a, d),
                mass_widths,
            );
            for (i, (face, _)) in f.faces.iter().enumerate() {
                record(
                    &name,
                    &format!("face_area[{i}]"),
                    digits,
                    |d| measure::face_area(&f.model, *face, &a, d),
                    |e| vec![width(e)],
                );
            }
        }
    }
    if std::env::var_os("CATALLAX_MEASURE_BOOLEAN_BENCH").as_deref()
        == Some(std::ffi::OsStr::new("1"))
    {
        for case in BOOLEAN_CASES {
            // Boolean construction is outside the measurement timer.
            let f = build_boolean(case);
            let name = format!("boolean.{case:?}");
            for digits in [6, 12, 18, 24, 30] {
                let a = accuracy(digits);
                record(
                    &name,
                    "volume",
                    digits,
                    |d| measure::volume(&f.model, f.body, &a, d),
                    |e| vec![width(e)],
                );
                record(
                    &name,
                    "area",
                    digits,
                    |d| measure::area(&f.model, f.body, &a, d),
                    |e| vec![width(e)],
                );
                record(
                    &name,
                    "mass_properties",
                    digits,
                    |d| measure::mass_properties(&f.model, f.body, &a, d),
                    mass_widths,
                );
                for (i, face) in f.model.faces_of_body(f.body).into_iter().enumerate() {
                    record(
                        &name,
                        &format!("face_area[{i}]"),
                        digits,
                        |d| measure::face_area(&f.model, face, &a, d),
                        |e| vec![width(e)],
                    );
                }
            }
        }
    }
}
