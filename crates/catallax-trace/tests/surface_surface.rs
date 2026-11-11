#![allow(missing_docs)]
//! Differential gate for SMLib's local surface-surface solvers
//! (`SmSurface::LocalSurfaceIntersect` and `LocalSurfaceSolve`) on B-spline
//! surfaces and planes, against saved outputs of the pinned C++ SMLib.

mod common;

use catallax_trace::extent::Extent2;
use catallax_trace::global_solver::{SolverOperation, SolverRefusal};
use catallax_trace::spline::BSplineSurface;
use catallax_trace::surface_surface::{SurfacePair, local_surface_intersect, local_surface_solve};
use catallax_trace::vector::Vector2;
use common::{Input, Output, run_gate, status_code};

fn domain(input: &mut Input) -> Extent2 {
    let (a, b, c, d) = (input.d(), input.d(), input.d(), input.d());
    Extent2::new(Vector2::new(a, b), Vector2::new(c, d))
}

/// `s` and a surface as `Input::surface` reads it, or `p` and a plane as
/// the plane gate reads it.
fn surface(input: &mut Input) -> BSplineSurface {
    if input.word() == "s" {
        return input.surface();
    }
    let kind = input.i();
    let (origin, x, y) = (input.v3(), input.v3(), input.v3());
    let scale = input.v2();
    let bounds = domain(input);
    if kind == 1 {
        BSplineSurface::canonical_plane(origin, x, y).expect("plane")
    } else {
        BSplineSurface::plane(origin, x, y, scale, bounds).expect("plane")
    }
}

fn answer(out: &mut Output, result: Result<Option<SurfacePair>, SolverRefusal>) {
    match result {
        Ok(found) => {
            out.i(common::SUCCESS);
            out.flag(found.is_some());
            if let Some(pair) = found {
                out.v2(pair.uv);
                out.v2(pair.other_uv);
                out.d(pair.value);
            }
        }
        Err(SolverRefusal::Failed(fault)) => out.i(status_code(fault)),
        Err(SolverRefusal::NotPorted(what)) => panic!("not ported: {what:?}"),
    }
}

fn run(op: &str, input: &mut Input, out: &mut Output) {
    let a = surface(input);
    let da = domain(input);
    let b = surface(input);
    let db = domain(input);
    let (guess, other_guess) = (input.v2(), input.v2());
    match op {
        "lsi" => {
            let tolerance = input.d();
            let has_normal = input.i() != 0;
            let normal = input.v3();
            let has_previous = input.i() != 0;
            let previous = input.v2();
            let boundary = input.i() != 0;
            answer(
                out,
                local_surface_intersect(
                    &a,
                    da,
                    &b,
                    db,
                    tolerance,
                    guess,
                    other_guess,
                    has_normal.then_some(normal),
                    has_previous.then_some(previous),
                    boundary,
                ),
            );
        }
        "lss" => {
            let operation = SolverOperation::from_code(input.i() as usize).expect("operation");
            answer(
                out,
                local_surface_solve(&a, da, &b, db, operation, guess, other_guess),
            );
        }
        other => panic!("unknown operation {other}"),
    }
}

#[test]
fn local_surface_solvers_match_smlib() {
    let report = run_gate(
        include_str!("fixtures/surface-surface/cases.txt"),
        |_| 16,
        run,
    );
    report.assert_passed(&[
        "tangent",
        "touching",
        "singular",
        "near-singular",
        "boundary",
    ]);
}

/// Regressions from the step 7 part 3 review.
mod review {
    use catallax_trace::extent::Extent2;
    use catallax_trace::global_solver::{SolverOperation, SolverRefusal};
    use catallax_trace::spline::{BSplineSurface, ControlPoint, NO_WEIGHT};
    use catallax_trace::surface_surface::{local_surface_intersect, local_surface_solve};
    use catallax_trace::vector::{Vector2, Vector3};

    fn bilinear(points: [(f64, f64, f64); 4], weight: Option<f64>) -> BSplineSurface {
        let knots = vec![0.0, 0.0, 1.0, 1.0];
        let points = points
            .iter()
            .map(|&(x, y, z)| ControlPoint {
                coordinates: match weight {
                    Some(w) => [x * w, y * w, z * w, w],
                    None => [x, y, z, NO_WEIGHT],
                },
            })
            .collect();
        BSplineSurface::new([1, 1], [knots.clone(), knots], [2, 2], points).expect("bilinear")
    }

    fn domain(u1: f64, v1: f64) -> Extent2 {
        Extent2::new(Vector2::new(0.0, 0.0), Vector2::new(u1, v1))
    }

    /// The side-curve search must not answer outside the requested
    /// domains. (SMLib answers `(0, 0.75)` on the first surface, outside its
    /// requested `[0, 1] × [0, 0.5]`.)
    #[test]
    fn side_curve_answers_stay_in_the_requested_domains() {
        let s = bilinear(
            [
                (0.0, 0.0, 0.0),
                (0.0, 1.0, 0.0),
                (1.0, 0.0, 0.0),
                (1.0, 1.0, 0.0),
            ],
            None,
        );
        let t = bilinear(
            [
                (0.0, 0.75, -0.5),
                (0.0, 0.75, 0.5),
                (1.0, 0.75, -0.5),
                (1.0, 0.75, 0.5),
            ],
            None,
        );
        let (ds, dt) = (domain(1.0, 0.5), domain(1.0, 1.0));
        let found = local_surface_intersect(
            &s,
            ds,
            &t,
            dt,
            1e-6,
            Vector2::new(0.0, 0.25),
            Vector2::new(0.0, 0.5),
            None,
            None,
            true,
        )
        .expect("search");
        if let Some(pair) = found {
            assert!(ds.contains_point(pair.uv, 0.0), "{pair:?}");
            assert!(dt.contains_point(pair.other_uv, 0.0), "{pair:?}");
        }
    }

    /// An evaluator failure at the guess is a typed failure, not an answer.
    /// (SMLib's solver reports the guess as found, and SMLib answers it.)
    #[test]
    fn normal_pair_evaluator_failure_is_refused() {
        let zero = bilinear([(0.0, 0.0, 0.0); 4], Some(1.0));
        let above = bilinear(
            [
                (0.0, 0.0, 1.0),
                (0.0, 1.0, 1.0),
                (1.0, 0.0, 1.0),
                (1.0, 1.0, 1.0),
            ],
            None,
        );
        let guess = Vector2::new(0.5, 0.5);
        let result = local_surface_solve(
            &zero,
            domain(1.0, 1.0),
            &above,
            domain(1.0, 1.0),
            SolverOperation::Normalize,
            guess,
            guess,
        );
        assert!(
            matches!(result, Err(SolverRefusal::Failed(_))),
            "{result:?}"
        );
        let _ = Vector3::default();
    }
}
