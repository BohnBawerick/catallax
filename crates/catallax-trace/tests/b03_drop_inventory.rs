#![allow(missing_docs)]
//! Replay only the two shipped refusal rows; no donor execution or new answers.
mod common;
use catallax_trace::{
    extent::Extent2,
    global_solver::SolverRefusal,
    spline::BSplineSurface,
    surface_drop_curve::{drop_and_trim_curve, drop_curve_onto},
    vector::Vector2,
};
fn domain(input: &mut common::Input) -> Extent2 {
    Extent2::new(
        Vector2::new(input.d(), input.d()),
        Vector2::new(input.d(), input.d()),
    )
}
#[test]
fn public_remaining_drop_rows_preserve_their_refusals() {
    let rows = include_str!("fixtures/drop-curve/cases.txt")
        .lines()
        .filter(|line| line.starts_with("drefuse "))
        .collect::<Vec<_>>();
    assert_eq!(rows.len(), 2);
    let fixture = rows.join("\n");
    let report = common::run_gate(
        &fixture,
        |_| 16,
        |op, input, output| {
            assert_eq!(op, "drefuse");
            let surface = if input.word() == "s" {
                input.surface()
            } else {
                let kind = input.i();
                let (origin, x, y) = (input.v3(), input.v3(), input.v3());
                let scale = input.v2();
                let bounds = domain(input);
                if kind == 1 {
                    BSplineSurface::canonical_plane(origin, x, y).unwrap()
                } else {
                    BSplineSurface::plane(origin, x, y, scale, bounds).unwrap()
                }
            };
            let domain = domain(input);
            let curve = input.curve();
            let interval = input.interval();
            let tolerance = input.d();
            let keep_all = input.i() != 0;
            let on_surface = input.i() != 0;
            if on_surface {
                let direct =
                    drop_curve_onto(&surface, domain, &curve, interval, tolerance, false).unwrap();
                println!(
                    "B03 OFF direct general drop keep_all=false: curves={} max_drop={} max_deviation={}",
                    direct.curves.len(),
                    direct.max_drop,
                    direct.max_deviation
                );
            }
            let answer = drop_and_trim_curve(
                &surface, domain, &curve, interval, tolerance, keep_all, on_surface,
            );
            println!(
                "B03 on_surface={on_surface} complete refusal={:?}",
                answer.as_ref().err()
            );
            assert!(matches!(answer, Err(SolverRefusal::NotPorted(_))));
            output.i(-1);
        },
    );
    report.assert_passed(&["regular"]);
}
