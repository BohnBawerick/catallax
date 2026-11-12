#![allow(missing_docs)]
//! Differential gate for the line and Hermite polar boxes.

mod common;

use catallax_trace::extent::Extent1;
use catallax_trace::hermite::Hermite;
use catallax_trace::line::Line;
use catallax_trace::polar_box::{PolarBox, PolarState};
use catallax_trace::vector::Vector2;
use common::{Input, Output, SUCCESS, assert_same_with_workers, run_gate};

fn polar(out: &mut Output, box_: PolarBox) {
    out.i(SUCCESS);
    out.i(box_.state as i64);
    if box_.state == PolarState::Uninitialized {
        return;
    }
    for axis in box_.basis {
        out.v3(axis);
    }
    out.v2(box_.min);
    out.v2(box_.max);
}

fn evaluate(op: &str, input: &mut Input, out: &mut Output) {
    match op {
        "line" => {
            let mode = input.i();
            let (p, v) = (input.v3(), input.v3());
            let mut curve = if mode == 0 {
                let (range, scale) = (input.interval(), input.d());
                let (dim, inside) = (input.i(), input.i());
                out.i(SUCCESS);
                Line::new(p, v, range, scale, dim == 2, inside != 0)
            } else {
                let (dim, has) = (input.i(), input.i());
                let range = input.interval();
                let Some(curve) =
                    out.status(Line::segment(p, v, dim == 2, (has != 0).then_some(range)))
                else {
                    return;
                };
                curve
            };
            input.d();
            input.d();
            input.i();
            input.d();
            let part = input.interval();
            let reparam = input.i() != 0;
            let next = input.interval();
            polar(out, curve.polar_box(part));
            if reparam {
                curve.edit_parameterization(next);
                out.i(SUCCESS);
                let now = curve.natural_interval();
                polar(
                    out,
                    curve.polar_box(Extent1::new(now.evaluate(0.2), now.evaluate(0.9))),
                );
            }
        }
        "hermite" => {
            let planar = input.i() != 0;
            let (p1, d1, p2, d2) = (input.v3(), input.v3(), input.v3(), input.v3());
            input.d();
            input.i();
            input.i();
            let flat = |v: catallax_trace::vector::Vector3| Vector2::new(v.x, v.y);
            let curve = if planar {
                Hermite::planar(flat(p1), flat(d1), flat(p2), flat(d2))
            } else {
                Hermite::new(p1, d1, p2, d2)
            };
            // SMLib ignores the interval, including a collapsed interval.
            for _ in 0..3 {
                polar(out, curve.polar_box());
            }
        }
        _ => panic!("unknown polar-box operation {op}"),
    }
}

#[test]
fn curve_polar_gate() {
    let fixture = include_str!("fixtures/curve-polar/cases.txt");
    let report = run_gate(fixture, |_| 16, evaluate);
    assert!(report.failures.is_empty(), "{:#?}", report.failures);
    assert_same_with_workers(fixture, |_| 16, evaluate);
}
