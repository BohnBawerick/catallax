#![allow(missing_docs)]
//! Differential gate for SMLib's surface cache (`SmTrimSrfCache` built
//! through `CacheMakeOrValidate` without a face), `sm_ComputeNetConstants`,
//! `SmSurface::CalculateContinuities` and `SmBSplineSurface::
//! CalculateBoundingBox` against saved outputs of the pinned C++ SMLib.

mod common;

use catallax_trace::extent::{Extent2, Extent3};
use catallax_trace::geometry::Continuity;
use catallax_trace::pseudo_box::PseudoBox;
use catallax_trace::surface_cache::{
    BoundaryPlanes, CacheRefusal, SurfaceCache, SurfaceNodeData, boundary_planes, bounding_boxes,
    continuities, net_constants,
};
use common::{Input, Output, run_gate, status_code};

fn extent3(out: &mut Output, e: Extent3) {
    out.v3(e.min);
    out.v3(e.max);
}

fn extent2(out: &mut Output, e: Extent2) {
    out.v2(e.min);
    out.v2(e.max);
}

fn pseudo(out: &mut Output, b: &PseudoBox) {
    for basis in b.basis {
        out.v3(basis);
    }
    for interval in b.intervals {
        out.d(interval.min);
        out.d(interval.max);
    }
}

fn classes(out: &mut Output, c: &[Continuity]) {
    out.i(c.len() as i64);
    for class in c {
        out.i(*class as i64);
    }
}

fn dump(out: &mut Output, cache: &SurfaceCache, node: usize) {
    let tree = cache
        .trees
        .get(catallax_trace::solution::NodeId(node as u32))
        .expect("node");
    out.i(tree.aux as i64);
    extent3(out, tree.bbox);
    match cache.node_data[node] {
        SurfaceNodeData::Patch(index) => {
            let patch = &cache.patches[index];
            extent2(out, patch.uv_domain);
            pseudo(out, &patch.pseudo_box);
            for class in patch.u_continuities.iter().chain(&patch.v_continuities) {
                out.i(*class as i64);
            }
            out.i(patch.split_axis as i64);
            out.flag(patch.chord_height_satisfied);
            out.flag(patch.angle_satisfied);
        }
        SurfaceNodeData::Head(index) => {
            let head = &cache.heads[index];
            extent2(out, head.uv_domain);
            out.i(head.split_axis as i64);
        }
    }
    out.flag(tree.children.is_some());
    if let Some([a, b]) = tree.children {
        dump(out, cache, a.0 as usize);
        dump(out, cache, b.0 as usize);
    }
}

fn evaluate(op: &str, input: &mut Input, out: &mut Output) {
    let surface = input.surface();
    match op {
        "cache" => match SurfaceCache::for_trimmed_surface(&surface) {
            Ok(cache) => {
                out.i(common::SUCCESS);
                classes(out, &cache.u_continuities);
                classes(out, &cache.v_continuities);
                dump(out, &cache, cache.top.0 as usize);
            }
            Err(CacheRefusal::Failed(fault)) => out.i(status_code(fault)),
        },
        "net" => {
            let chords = net_constants(&surface, true, false, true);
            out.i(common::SUCCESS);
            out.d(chords.u_chord_height);
            out.d(chords.v_chord_height);
            out.d(chords.uv_chord_height);
            let angles = net_constants(&surface, false, true, false);
            out.i(common::SUCCESS);
            out.d(angles.u_angle_deg);
            out.d(angles.v_angle_deg);
        }
        "cont" => {
            for axis in 0..2 {
                if let Some((minimum, c)) = out.status(continuities(&surface, axis)) {
                    out.i(minimum as i64);
                    classes(out, &c);
                }
            }
        }
        "planes" => {
            out.i(common::SUCCESS);
            let planes =
                boundary_planes(&surface, surface.natural_domain()).unwrap_or(BoundaryPlanes {
                    planes: [None; 4],
                    degenerate: [false; 4],
                });
            for side in 0..4 {
                out.i(common::SUCCESS);
                out.flag(planes.planes[side].is_some());
                out.flag(planes.degenerate[side]);
                if let Some(plane) = planes.planes[side] {
                    out.v3(plane.base);
                    out.v3(plane.normal);
                    out.d(plane.tolerance);
                    for sample in plane.samples {
                        out.v3(sample);
                    }
                }
            }
        }
        "bbox" => {
            let (bbox, pseudo_box) = bounding_boxes(&surface, None);
            out.i(common::SUCCESS);
            extent3(out, bbox);
            pseudo(out, &pseudo_box);
        }
        other => panic!("unknown operation {other}"),
    }
}

#[test]
fn surface_cache_matches_smlib() {
    let report = run_gate(
        include_str!("fixtures/surface-cache/cases.txt"),
        |_| 16,
        evaluate,
    );
    report.assert_passed(&["singular", "near-singular", "boundary"]);
}
