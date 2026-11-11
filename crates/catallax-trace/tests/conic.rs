#![allow(missing_docs)]
//! Differential gate for SMLib's cubic and quartic solvers
//! (`SmPolynomial`) and the planar conics of the Dupin indicatrix
//! (`SmConic`) against saved outputs of the pinned C++ SMLib.

mod common;

use catallax_trace::conic::{Conic, ConicKind};
use catallax_trace::polynomial::{solve_cubic, solve_quartic};
use catallax_trace::vector::Vector2;
use common::{Input, Output, run_gate};

fn conic(input: &mut Input) -> Conic {
    let kind = match input.i() {
        0 => ConicKind::Ellipse,
        1 => ConicKind::Hyperbola,
        2 => ConicKind::Parabola,
        _ => ConicKind::Line,
    };
    let (a, b) = (input.d(), input.d());
    let mut conic = Conic::new(kind, a, b);
    for _ in 0..input.i() {
        let (angle, x, y) = (input.d(), input.d(), input.d());
        conic.transform_2d(angle, Vector2::new(x, y));
    }
    conic
}

fn evaluate(op: &str, input: &mut Input, out: &mut Output) {
    match op {
        "cubic" => {
            let c = [input.d(), input.d(), input.d(), input.d()];
            let (count, roots) = solve_cubic(c, input.d());
            out.i(common::SUCCESS);
            out.i(count as i64);
            for &root in &roots[..count] {
                out.d(root);
            }
        }
        "quartic" => {
            let c = [input.d(), input.d(), input.d(), input.d(), input.d()];
            let (count, roots) = solve_quartic(c, input.d());
            out.i(common::SUCCESS);
            out.i(count as i64);
            for &root in &roots[..count] {
                out.d(root);
            }
        }
        "conic" | "conicx" => {
            let one = conic(input);
            let two = conic(input);
            let tolerance = input.d();
            for t in [-0.75, -0.25, 0.3, 0.9, 1.5, 2.6] {
                let p = one.evaluate_point(t);
                out.v2(p);
                let back = one.invert_point(p);
                out.flag(back.is_some());
                if let Some(back) = back {
                    out.d(back);
                }
            }
            for c in two.implicit_equation() {
                out.d(c);
            }
            let found = one.intersect_conic(&two, tolerance);
            out.i(common::SUCCESS);
            out.i(found.len() as i64);
            for (t, t2, distance) in found {
                out.d(t);
                out.d(t2);
                out.flag(distance <= tolerance);
            }
        }
        other => panic!("unknown operation {other}"),
    }
}

#[test]
fn conic_matches_smlib() {
    // `conicx` marks the two documented exceptions (README.md): one-unit
    // cube-root differences between `libm` and the C library, magnified by
    // the quartic. Every other case keeps the 16-unit gate.
    let budget = |op: &str| if op == "conicx" { 32 } else { 16 };
    let report = run_gate(include_str!("fixtures/conic/cases.txt"), budget, evaluate);
    report.assert_passed(&[
        "tangent",
        "touching",
        "singular",
        "near-singular",
        "boundary",
    ]);
}

/// Pins the cause of the `conicx` exceptions: on these resolvent-cubic
/// arguments from the two pairs, the port's cube root (`libm`) and the C
/// library's (as SMLib links it, recorded from the same machine) differ by
/// exactly one unit in the last place. If the port's cube root changes, this
/// test and the exceptions must be revisited.
#[test]
fn conic_exception_cube_roots() {
    use catallax_trace::scalar::cbrt;
    // (argument, port, C library)
    let pins: [(u64, u64, u64); 6] = [
        (0x3fae99fe6acb3d3a, 0x3fd905dd26a1019f, 0x3fd905dd26a1019e),
        (0xbf957bdaf4f195eb, 0xbfd1a6ce65be2cb8, 0xbfd1a6ce65be2cb7),
        (0x400fb51f9115d392, 0x3ff9521fc89cd3ed, 0x3ff9521fc89cd3ec),
        (0xbff642a058ea294b, 0xbff1dc99ff55c5a1, 0xbff1dc99ff55c5a2),
        (0x3fa7a15caff47e26, 0x3fd6f4f7389f9c44, 0x3fd6f4f7389f9c43),
        (0x3f379793922f07c0, 0x3fb23608c482303c, 0x3fb23608c482303b),
    ];
    for (argument, port, library) in pins {
        let root = cbrt(f64::from_bits(argument)).to_bits();
        assert_eq!(root, port, "cube root of {argument:#x}");
        assert_eq!(
            port.abs_diff(library),
            1,
            "one unit apart for {argument:#x}"
        );
    }
}
