//! Private-helper witnesses; not a claim that the public solver emits these ranges.
use super::*;
use crate::spline::ControlPoint;
use crate::surface_eval::evaluate_point;

fn cylinder(nonuniform: bool) -> crate::spline::BSplineSurface {
    let r = 0.5_f64.sqrt();
    let ring = [
        (1., 0., 1.),
        (1., 1., r),
        (0., 1., 1.),
        (-1., 1., r),
        (-1., 0., 1.),
        (-1., -1., r),
        (0., -1., 1.),
        (1., -1., r),
        (1., 0., 1.),
    ];
    let heights = if nonuniform {
        [0., 1., 1. + 1e-7]
    } else {
        [0., 0.5, 1.]
    };
    let points = ring
        .into_iter()
        .flat_map(|(x, y, w)| {
            heights.map(move |z| ControlPoint {
                coordinates: [x * w, y * w, z * w, w],
            })
        })
        .collect();
    crate::spline::BSplineSurface::new(
        [2, 1],
        [
            vec![0., 0., 0., 1., 1., 2., 2., 3., 3., 4., 4., 4.],
            vec![0., 0., 0.5, 1., 1.],
        ],
        [9, 3],
        points,
    )
    .unwrap()
}
fn helper(nonuniform: bool, pretrim: bool, reclip: bool, swap: bool) {
    let cylinder = cylinder(nonuniform);
    let plane = crate::spline::BSplineSurface::canonical_plane(
        Vector3::new(1., 0., 0.),
        Vector3::new(0., 1., 0.),
        Vector3::new(0., 0., 1.),
    )
    .unwrap();
    let domain = Extent2::new(
        Vector2::new(0., if pretrim { 0.125 } else { 0. }),
        Vector2::new(4., 1.),
    );
    let pdomain = Extent2::new(
        Vector2::new(-2., if reclip { 0.25 } else { -2. }),
        Vector2::new(2., 2.),
    );
    let owner = usize::from(swap);
    let intersector = if swap {
        Intersector::new(&plane, pdomain, &cylinder, domain)
    } else {
        Intersector::new(&cylinder, domain, &plane, pdomain)
    };
    let sides = match intersector.side_curves(owner) {
        Ok(sides) => sides,
        Err(SolverRefusal::NotPorted(_)) if nonuniform && pretrim => {
            println!("B02 P2 stage=source-pretrim result=typed-refusal");
            return;
        }
        other => panic!("side construction: {:?}", other.err()),
    };
    let mut solution = crate::solution::Solution::default();
    solution.start.parameters[0] = if reclip { 0. } else { 0.125 };
    solution.end.parameters[0] = 1.;
    let mut curves = Vec::new();
    let result = Tracer::default().side_range_curve(
        &intersector,
        &sides[0],
        (owner, 0),
        solution,
        &mut curves,
        &mut Vec::new(),
    );
    match result {
        Err(SolverRefusal::NotPorted(_)) if nonuniform => {
            assert!(curves.is_empty());
            println!(
                "B02 P2 stage={} swap={swap} result=typed-refusal",
                if reclip {
                    "opposite-drop-reclip"
                } else {
                    "initial-clip"
                }
            );
        }
        Ok(()) => {
            assert_eq!(curves.len(), 1);
            let curve = &curves[0];
            let interval = curve.curve.natural_interval();
            let mut max_gap = 0_f64;
            for i in 0..=32 {
                let t = interval.evaluate(f64::from(i) / 32.);
                let p = curve.curve.evaluate_point(t).unwrap();
                let uv = curve.uv_curves[owner]
                    .as_ref()
                    .unwrap()
                    .evaluate_point(t)
                    .unwrap();
                let q = evaluate_point(&cylinder, Vector2::new(uv.x, uv.y)).unwrap();
                max_gap = max_gap.max(p.distance(q));
            }
            let t = 0.56;
            let p = curve.curve.evaluate_point(t).unwrap();
            let uv = curve.uv_curves[owner]
                .as_ref()
                .unwrap()
                .evaluate_point(t)
                .unwrap();
            let q = evaluate_point(&cylinder, Vector2::new(uv.x, uv.y)).unwrap();
            println!(
                "B02 P2 nonuniform={nonuniform} pretrim={pretrim} reclip={reclip} swap={swap} t=.56 spatial_z={} lifted_z={} max_gap={max_gap}",
                p.z, q.z
            );
            assert!(max_gap <= 1e-10, "owning UV lost native parameter identity");
        }
        other => panic!("unexpected helper result: {other:?}"),
    }
}
#[test]
fn initial_clip_keeps_correspondence_or_refuses() {
    helper(true, false, false, false);
}
#[test]
fn source_pretrim_keeps_correspondence_or_refuses() {
    helper(true, true, false, false);
}
#[test]
fn opposite_side_reclip_keeps_correspondence_or_refuses() {
    for swap in [false, true] {
        helper(true, false, true, swap);
    }
}
#[test]
fn uniform_parameter_control_still_builds() {
    for swap in [false, true] {
        helper(false, false, false, swap);
        helper(false, true, false, swap);
        helper(false, false, true, swap);
    }
}
