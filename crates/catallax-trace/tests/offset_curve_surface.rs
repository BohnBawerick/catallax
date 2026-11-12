#![allow(missing_docs)]
//! Differential gate for SMLib's local curve-surface solvers on offset
//! surfaces (`SmSurface::LocalCurveIntersect` and `LocalCurveSolve` with an
//! `SmIsoCurve` of another `SmOffsetSurface`, as `SmFilletSolver::StepSolve`
//! asks them for a clipped step) and its singularity tests (`IsSingularity`,
//! `GetSingularities`) on offset surfaces, against saved outputs of the pinned
//! C++ SMLib.

mod common;

use catallax_trace::curve_surface::{
    CurveSurfacePoint, Singular, get_singularities, is_singularity, local_curve_intersect,
    local_curve_solve,
};
use catallax_trace::extent::Extent2;
use catallax_trace::global_solver::{SolverOperation, SolverRefusal};
use catallax_trace::scalar::ZONE_TOL_3D;
use catallax_trace::spline::{BSplineSurface, ControlPoint};
use catallax_trace::surface_offset::OffsetSurface;
use catallax_trace::surface_solve::SurfaceIsoCurve;
use catallax_trace::vector::Vector2;
use common::{Input, Output, run_gate, status_code};

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

fn found(out: &mut Output, result: Result<Option<CurveSurfacePoint>, SolverRefusal>) {
    match result {
        Ok(point) => {
            out.i(common::SUCCESS);
            out.flag(point.is_some());
            if let Some(point) = point {
                out.d(point.t);
                out.v2(point.uv);
                out.d(point.distance);
            }
        }
        Err(SolverRefusal::Failed(fault)) => out.i(status_code(fault)),
        Err(SolverRefusal::NotPorted(operation)) => panic!("not ported: {operation:?}"),
    }
}

fn run(op: &str, input: &mut Input, out: &mut Output) {
    if op == "oc-sing" {
        let distance = input.d();
        let base = raw_surface(input);
        let offset = OffsetSurface::new(&base, distance);
        let uv = input.v2();
        let tolerance = input.d();
        let point_only = input.i() != 0;
        let direction = is_singularity(&offset, uv, tolerance, ZONE_TOL_3D, point_only)
            .expect("singularity test");
        out.flag(direction.is_some());
        out.i(match direction {
            Some(Singular::U) => 0,
            Some(Singular::V) => 1,
            Some(Singular::Both) => 3,
            None => 4,
        });
        let (flags, poles) = get_singularities(&offset, ZONE_TOL_3D);
        out.i(i64::from(flags));
        for pole in poles {
            match pole {
                Some(point) => out.v3(point),
                None => out.i(0),
            }
        }
        return;
    }
    let distance_a = input.d();
    let base_a = raw_surface(input);
    let distance_b = input.d();
    let base_b = raw_surface(input);
    let (offset_a, offset_b) = (
        OffsetSurface::new(&base_a, distance_a),
        OffsetSurface::new(&base_b, distance_b),
    );
    let axis = input.i() as usize;
    let iso = input.d();
    let clamp = input.i() != 0;
    let curve = SurfaceIsoCurve::new(&offset_a, axis, iso, true, clamp);
    let (low, high) = (input.v2(), input.v2());
    let domain = Extent2::new(low, high);
    let interval = curve.natural_interval();
    match op {
        "oc-intersect" => {
            let (tolerance, t, uv) = (input.d(), input.d(), input.v2());
            found(
                out,
                local_curve_intersect(
                    &offset_b,
                    domain,
                    &curve,
                    interval,
                    tolerance,
                    (t, uv),
                    ZONE_TOL_3D,
                ),
            );
        }
        "oc-solve" => {
            let operation =
                SolverOperation::from_code(input.i() as usize).expect("SMLib operation code");
            let (_tolerance, t, uv): (f64, f64, Vector2) = (input.d(), input.d(), input.v2());
            found(
                out,
                local_curve_solve(&offset_b, domain, &curve, interval, operation, (t, uv)),
            );
        }
        other => panic!("unknown operation {other}"),
    }
}

#[test]
fn offset_curve_surface_solvers_match_smlib() {
    let report = run_gate(
        include_str!("fixtures/offset-curve-surface/cases.txt"),
        |_| 0,
        run,
    );
    assert!(report.cases >= 450, "the gate lost cases: {}", report.cases);
    report.assert_passed(&[
        "regular",
        "tangent",
        "touching",
        "singular",
        "near-singular",
        "boundary",
    ]);
}
