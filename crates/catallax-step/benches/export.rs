//! Non-gating STEP export time, file size and fitting control-point measurements.
#[path = "../tests/common/mod.rs"]
mod common;
use common::*;
use std::time::Instant;

fn measure(label: &str, m: &Model, b: BodyId, tolerance: &str) {
    let o = StepOptions {
        curve_tolerance: r(tolerance),
        ..options(17)
    };
    let mut bytes = Vec::new();
    let mut dec = Decider::default();
    let start = Instant::now();
    let result = write_step(m, &[b], &o, &mut dec, &mut bytes);
    let elapsed = start.elapsed();
    match result {
        Ok(a) => {
            // Parsing and counting are outside the measured fitting/export time.
            let s = Step::parse(&bytes, &o.product_name);
            let cps: usize = s
                .ids("B_SPLINE_CURVE_WITH_KNOTS")
                .iter()
                .map(|id| list(&s.entities[id].args[2]).len())
                .sum();
            println!(
                "{label},{tolerance},{},{},{},{},{},{},{},answered",
                a.value.faces,
                bytes.len(),
                elapsed.as_nanos(),
                a.value.approximated_curves,
                cps,
                a.value.max_deviation,
                a.marks.len()
            );
        }
        Err(e) => println!(
            "{label},{tolerance},{},0,{},0,0,0,0,{}",
            m.faces_of_body(b).len(),
            elapsed.as_nanos(),
            e.code
        ),
    }
}
fn main() {
    println!(
        "case,tolerance,faces,bytes,export_ns,splines,control_points,max_deviation,marks,outcome"
    );
    for (i, p) in primitives().iter().enumerate() {
        let mut m = Model::new();
        let b = p.build(&mut m, &Pose::identity());
        measure(&format!("primitive-{i}"), &m, b, "1e-7");
    }
    for count in [1, 2, 4, 8, 16, 32, 64, 128] {
        let mut m = Model::new();
        let mut lumps = Vec::new();
        for i in 0..count {
            let b = block(
                &mut m,
                [n(2), n(3), n(5)],
                &translated(&(i * 4).to_string(), "0", "0"),
            )
            .unwrap();
            lumps.extend_from_slice(&m.body(b).lumps);
        }
        let b = body(&mut m, lumps);
        measure(&format!("blocks-{count}"), &m, b, "1e-7");
    }
    for case in 1..=40 {
        for op in [
            catallax_boolean::BoolOp::Unite,
            catallax_boolean::BoolOp::Subtract,
            catallax_boolean::BoolOp::Intersect,
        ] {
            let (mut m, a, b) = boolean_pair(case);
            let label = format!("boolean-{case}-{op:?}");
            match catallax_boolean::boolean(&mut m, a, b, op, &mut Decider::default()) {
                Ok(a) => measure(&label, &m, a.body, "1e-7"),
                Err(e) => println!("{label},1e-7,0,0,0,0,0,0,0,{}", e.refusal.code),
            }
        }
    }
    for case in 17..=26 {
        let (mut m, a, b) = branch_pair(case);
        match catallax_boolean::boolean(
            &mut m,
            a,
            b,
            catallax_boolean::BoolOp::Intersect,
            &mut Decider::default(),
        ) {
            Ok(a) => {
                for digits in 6..=12 {
                    measure(
                        &format!("intersection-{case}"),
                        &m,
                        a.body,
                        &format!("1e-{digits}"),
                    );
                }
            }
            Err(e) => println!("intersection-{case},0,0,0,0,0,0,0,0,{}", e.refusal.code),
        }
    }
}
