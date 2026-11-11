#![allow(missing_docs)]
//! Differential gate for B-spline surface evaluation
//! (`SmBSplineSurface::Evaluate`, `EvaluatePoint`, `EvaluateSimple`) against
//! saved outputs of the pinned C++ SMLib.

mod common;

use catallax_trace::surface_eval::{evaluate, evaluate_normal, evaluate_point, evaluate_simple};
use catallax_trace::vector::Vector2;
use common::{Input, Output, run_gate};

fn run(op: &str, input: &mut Input, out: &mut Output) {
    let surface = input.surface();
    let uv = Vector2::new(input.d(), input.d());
    if op == "point" {
        if let Some(point) = out.status(evaluate_point(&surface, uv)) {
            out.v3(point);
        }
        return;
    }
    if op == "normal" {
        let sides = [input.i() != 0, input.i() != 0];
        out.status(Ok(()));
        out.v3(evaluate_normal(&surface, uv, sides));
        return;
    }
    let (hu, hv) = (input.i() as usize, input.i() as usize);
    let sides = [input.i() != 0, input.i() != 0];
    let (upper, nonzero) = (input.i() != 0, input.i() != 0);
    let result = match op {
        "eval" => evaluate(&surface, uv, hu, hv, sides, upper, nonzero, true),
        "simple" => evaluate_simple(&surface, uv, hu, hv, sides, upper),
        other => panic!("unknown operation {other}"),
    };
    if let Some(cells) = out.status(result) {
        for k in 0..=hu {
            for l in 0..=hv {
                if !upper || k + l <= hu {
                    out.v3(cells[k * (hv + 1) + l]);
                }
            }
        }
    }
}

#[test]
fn surface_evaluation_matches_smlib() {
    let report = run_gate(include_str!("fixtures/surface-eval/cases.txt"), |_| 4, run);
    report.assert_passed(&["touching", "singular", "boundary"]);
}

mod unsupported_forms {
    use catallax_trace::Fault;
    use catallax_trace::spline::{BSplineSurface, ControlPoint, NO_WEIGHT};
    use catallax_trace::surface_eval::{
        derivatives, eval_point, evaluate, evaluate_normal, evaluate_point, evaluate_simple,
    };
    use catallax_trace::vector::Vector2;

    /// A polynomial surface with `P[i][j] = (i, j, 0)`.
    fn surface(degrees: [usize; 2], knots: [Vec<f64>; 2]) -> BSplineSurface {
        let counts = [
            knots[0].len() - degrees[0] - 1,
            knots[1].len() - degrees[1] - 1,
        ];
        let points = (0..counts[0])
            .flat_map(|i| {
                (0..counts[1]).map(move |j| ControlPoint {
                    coordinates: [i as f64, j as f64, 0.0, NO_WEIGHT],
                })
            })
            .collect();
        BSplineSurface::new(degrees, knots, counts, points).expect("accepted by the constructor")
    }

    fn clamped(degree: usize) -> Vec<f64> {
        let mut knots = vec![0.0; degree + 1];
        knots.extend(vec![1.0; degree + 1]);
        knots
    }

    fn assert_refused(s: &BSplineSurface, uv: Vector2) {
        assert_eq!(evaluate_point(s, uv), Err(Fault::InvalidInput));
        assert_eq!(eval_point(s, uv, [false, false]), Err(Fault::InvalidInput));
        assert_eq!(
            derivatives(s, uv, [true, true], false, 1, 1),
            Err(Fault::InvalidInput)
        );
        assert_eq!(
            evaluate(s, uv, 2, 2, [true, true], false, true, true),
            Err(Fault::InvalidInput)
        );
        assert_eq!(
            evaluate_simple(s, uv, 1, 1, [true, true], false),
            Err(Fault::InvalidInput)
        );
        // The normal has no status; it must not panic.
        let _ = evaluate_normal(s, uv, [true, true]);
    }

    /// SMLib's NLib sizes its basis work arrays for degree 17
    /// (`NL_MAXDEG`); a higher degree the constructor accepts is refused
    /// rather than indexed past the work arrays.
    #[test]
    fn degree_above_nlib_limit_is_refused() {
        let s = surface([18, 1], [clamped(18), clamped(1)]);
        for uv in [(1.0, 0.5), (0.0, 0.5), (0.4, 0.5)] {
            assert_refused(&s, Vector2::new(uv.0, uv.1));
        }
    }

    /// Degree 17, NLib's limit, still evaluates: a polynomial net
    /// `P[i][j] = (i, j, 0)` reproduces `(17 u, v, 0)`.
    #[test]
    fn degree_at_nlib_limit_evaluates() {
        let s = surface([17, 1], [clamped(17), clamped(1)]);
        for u in [0.0, 0.3, 1.0] {
            let point = evaluate_point(&s, Vector2::new(u, 0.5)).expect("degree 17");
            assert!((point.x - 17.0 * u).abs() < 1e-12, "{point:?}");
            assert!((point.y - 0.5).abs() < 1e-12, "{point:?}");
            let cells = derivatives(&s, Vector2::new(u, 0.5), [true, true], false, 1, 1)
                .expect("degree 17 derivatives");
            assert!((cells[0].x - 17.0 * u).abs() < 1e-12, "{:?}", cells[0]);
        }
    }

    /// NLib's end shortcut assumes clamped knots, as SMLib's surfaces have.
    /// An unclamped knot vector the constructor accepts is refused, so the
    /// point and derivative paths cannot disagree on it.
    #[test]
    fn unclamped_knots_are_refused() {
        let s = surface([2, 1], [vec![-2.0, -1.0, 0.0, 1.0, 2.0, 3.0], clamped(1)]);
        for uv in [(0.0, 0.5), (1.0, 0.5), (0.5, 0.5)] {
            assert_refused(&s, Vector2::new(uv.0, uv.1));
        }
    }
}

mod unsupported_cache {
    use catallax_trace::Fault;
    use catallax_trace::spline::{BSplineSurface, ControlPoint, NO_WEIGHT};
    use catallax_trace::surface_cache::{CacheRefusal, SurfaceCache, SurfaceCacheSettings};

    /// The cache refuses the forms evaluation refuses, before decomposing.
    #[test]
    fn cache_refuses_unsupported_forms() {
        for (degrees, knots) in [
            (
                [18, 1],
                [
                    [vec![0.0; 19], vec![1.0; 19]].concat(),
                    vec![0.0, 0.0, 1.0, 1.0],
                ],
            ),
            (
                [2, 1],
                [
                    vec![-2.0, -1.0, 0.0, 1.0, 2.0, 3.0],
                    vec![0.0, 0.0, 1.0, 1.0],
                ],
            ),
        ] {
            let counts = [knots[0].len() - degrees[0] - 1, 2];
            let points = (0..counts[0] * 2)
                .map(|k| ControlPoint {
                    coordinates: [(k / 2) as f64, (k % 2) as f64, 0.0, NO_WEIGHT],
                })
                .collect();
            let s = BSplineSurface::new(degrees, knots, counts, points).expect("accepted");
            for cache in [
                SurfaceCache::for_trimmed_surface(&s),
                SurfaceCache::build(
                    &s,
                    SurfaceCacheSettings {
                        chord_height_tolerance: 0.1,
                        angle_tolerance_rad: 0.5,
                        min_side_length_ratio_uv: 0.001,
                    },
                ),
            ] {
                assert!(matches!(
                    cache,
                    Err(CacheRefusal::Failed(Fault::InvalidInput))
                ));
            }
        }
    }
}
