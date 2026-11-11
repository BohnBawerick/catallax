#![allow(missing_docs)]
//! Differential gate for `SmOffsetSurface` (`EvaluatePoint`, `EvaluateNormal`,
//! `Evaluate` and `EvaluateSimple`) against saved outputs of the pinned C++
//! SMLib.

mod common;

use catallax_trace::spline::{BSplineSurface, ControlPoint};
use catallax_trace::surface_offset::OffsetSurface;
use catallax_trace::vector::Vector2;
use common::{Input, Output, run_gate};

/// `du dv nu nv rational knots_u... knots_v... (x y z w)...`, as stored.
fn raw_surface(input: &mut Input) -> BSplineSurface {
    let degrees = [input.i() as usize, input.i() as usize];
    let counts = [input.i() as usize, input.i() as usize];
    let _rational = input.i();
    let knots = [0, 1].map(|axis| {
        (0..counts[axis] + degrees[axis] + 1)
            .map(|_| input.d())
            .collect()
    });
    let points = (0..counts[0] * counts[1])
        .map(|_| ControlPoint {
            coordinates: [input.d(), input.d(), input.d(), input.d()],
        })
        .collect();
    BSplineSurface::new(degrees, knots, counts, points).expect("valid test surface")
}

fn run(op: &str, input: &mut Input, out: &mut Output) {
    let distance = input.d();
    let surface = raw_surface(input);
    let uv = Vector2::new(input.d(), input.d());
    let offset = OffsetSurface::new(&surface, distance);
    match op {
        "offset-point" => {
            if let Some(point) = out.status(offset.evaluate_point(uv)) {
                out.v3(point);
            }
        }
        "offset-normal" => {
            let sides = [input.i() != 0, input.i() != 0];
            out.i(1000);
            out.v3(offset.evaluate_normal(uv, sides));
        }
        "offset-eval" => {
            let simple = input.i() != 0;
            let (hu, hv) = (input.i() as usize, input.i() as usize);
            let sides = [input.i() != 0, input.i() != 0];
            let upper = input.i() != 0;
            let non_zero = input.i() != 0;
            let cells = if simple {
                offset.evaluate_simple(uv, hu, hv, sides, upper)
            } else {
                offset.evaluate(uv, hu, hv, sides, upper, non_zero)
            };
            if let Some(cells) = out.status(cells) {
                for cell in cells {
                    out.v3(cell);
                }
            }
        }
        other => panic!("unknown operation {other}"),
    }
}

#[test]
fn surface_offset_matches_smlib() {
    let report = run_gate(
        include_str!("fixtures/surface-offset/cases.txt"),
        |_| 0,
        run,
    );
    report.assert_passed(&["tangent", "touching", "singular", "boundary"]);
}
