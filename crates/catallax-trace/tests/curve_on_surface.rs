#![allow(missing_docs)]
//! Differential gate for `SmCrvOnSurf` (the point, the surface parameters and
//! the surface normal along a parameter-space curve) and for the branch of
//! `SmGlobalSolver::BranchMayContainAnswers` that samples a curve on a surface
//! (`SM_SO_SURFACE_PROJECTED_INTERSECT`), against saved outputs of the pinned
//! C++ SMLib.

mod common;

use catallax_trace::curve_on_surface::CurveOnSurface;
use catallax_trace::curve_solve::PointSolveCurve;
use catallax_trace::extent::{Extent1, Extent3};
use catallax_trace::global_solver::{
    AuxData, GlobalSolver, SearchTrees, SolverOperation, SurfaceCurveSpans, TreeNode,
};
use catallax_trace::hermite::Hermite;
use catallax_trace::solution::NoObjects;
use catallax_trace::spline::{BSplineCurve, ControlPoint};
use catallax_trace::vector::Vector2;
use common::{Input, Output, run_gate};

enum UvCurve {
    Hermite(Hermite),
    Spline(BSplineCurve),
}

impl UvCurve {
    fn read(input: &mut Input) -> Self {
        if input.i() == 0 {
            let mut v = [0.0; 8];
            for value in &mut v {
                *value = input.d();
            }
            Self::Hermite(Hermite::planar(
                Vector2::new(v[0], v[1]),
                Vector2::new(v[2], v[3]),
                Vector2::new(v[4], v[5]),
                Vector2::new(v[6], v[7]),
            ))
        } else {
            let _dimension = input.i();
            let degree = input.i() as usize;
            let count = input.i() as usize;
            let knots = (0..count + degree + 1).map(|_| input.d()).collect();
            let points = (0..count)
                .map(|_| ControlPoint {
                    coordinates: [input.d(), input.d(), input.d(), input.d()],
                })
                .collect();
            Self::Spline(BSplineCurve::new(degree, knots, points).expect("valid test curve"))
        }
    }

    fn as_curve(&self) -> &dyn PointSolveCurve {
        match self {
            Self::Hermite(curve) => curve,
            Self::Spline(curve) => curve,
        }
    }
}

fn read_box(input: &mut Input) -> Extent3 {
    let (min, max) = (input.v3(), input.v3());
    Extent3 { min, max }
}

fn run(op: &str, input: &mut Input, out: &mut Output) {
    let surface = input.surface();
    let uv_curve = UvCurve::read(input);
    let curve = CurveOnSurface::new(uv_curve.as_curve(), &surface);
    match op {
        "cos-eval" => {
            let t = input.d();
            let left = input.i() != 0;
            if let Some(point) = out.status(curve.evaluate_point(t)) {
                out.v3(point);
            }
            if let Some(uv) = out.status(curve.evaluate_uv_point(t)) {
                out.v2(uv);
            }
            if let Some(([position, tangent], normal)) =
                out.status(curve.evaluate_surface_normal(t, left))
            {
                out.v2(position);
                out.v2(tangent);
                out.v3(normal);
            }
        }
        "cos-branch" => {
            let index = input.i();
            let first = read_box(input);
            let tolerance = input.d();
            // What the oracle read off its curve cache's node `index`.
            let second = read_box(input);
            let interval = Extent1::new(input.d(), input.d());
            let _ = index;
            let mut trees = SearchTrees::new();
            let node = |bbox| TreeNode {
                bbox,
                children: None,
                aux: AuxData::None,
            };
            let a = trees.add(node(first));
            let b = trees.add(node(second));
            let mut solver = GlobalSolver::new(&trees, &NoObjects);
            solver.operation = SolverOperation::SurfaceProjectedIntersect;
            solver.set_num_trees(2);
            solver.tolerance = tolerance;
            solver.surface_curve = Some(SurfaceCurveSpans {
                curve,
                spans: vec![(b, interval)],
            });
            let branch = [Some(a), Some(b), None];
            out.flag(
                solver
                    .base_branch_may_contain_answers(&branch)
                    .expect("a branch SMLib answers"),
            );
        }
        other => panic!("unknown operation {other}"),
    }
}

const FIXTURE: &str = include_str!("fixtures/curve-on-surface/cases.txt");

#[test]
fn curve_on_surface_matches_smlib() {
    let report = run_gate(FIXTURE, |_| 0, run);
    assert!(
        report.cases >= 1000,
        "the gate lost cases: {}",
        report.cases
    );
    report.assert_passed(&["regular", "tangent", "touching", "singular", "boundary"]);
}
