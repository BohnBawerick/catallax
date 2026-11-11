#![allow(missing_docs)]
//! Differential gate for SMLib's curve cache (`SmCurveCache` as
//! `SmCurve::CacheMakeOrValidate` builds it) against saved outputs of the
//! pinned C++ SMLib.

mod common;

use catallax_trace::curve_cache::{CurveCache, NodeData, SpanData};
use catallax_trace::solution::NodeId;
use common::{Input, Output, run_gate};

/// SMLib's approximation tolerance for a curve without an edge: half the
/// default zone tolerance. The oracle prints it first, so the gate checks it.
const APPROX_TOLERANCE: f64 = 5.0e-6;

fn span(out: &mut Output, span: &SpanData) {
    out.i(span.left as i64);
    out.i(span.right as i64);
    out.i(span.continuities[0] as i64);
    out.i(span.continuities[1] as i64);
    out.d(span.max_length);
    out.d(span.max_chord_height_squared);
    out.d(span.max_turning_angle_deg);
    out.d(span.interval.min);
    out.d(span.interval.max);
}

fn node(out: &mut Output, cache: &CurveCache, id: NodeId) {
    let tree = cache.trees.get(id).expect("tree node");
    out.v3(tree.bbox.min);
    out.v3(tree.bbox.max);
    let (kind, index) = match cache.node_data[id.0 as usize] {
        NodeData::Leaf(i) => (1, i as i64),
        NodeData::Parent(i) => (2, i as i64),
        NodeData::None => (0, -1),
    };
    out.i(kind);
    out.i(index);
    out.flag(tree.children.is_some());
    if let Some([a, b]) = tree.children {
        node(out, cache, a);
        node(out, cache, b);
    }
}

fn evaluate(op: &str, input: &mut Input, out: &mut Output) {
    assert_eq!(op, "cache");
    let curve = input.curve();
    out.d(APPROX_TOLERANCE);
    let Ok(cache) = CurveCache::for_curve(&curve, APPROX_TOLERANCE) else {
        out.i(0);
        return;
    };
    out.i(1);
    out.i(cache.continuities.len() as i64);
    for class in &cache.continuities {
        out.i(*class as i64);
    }
    out.i(cache.spans.len() as i64);
    for (data, pseudo) in cache.spans.iter().zip(&cache.pseudo_boxes) {
        span(out, data);
        for basis in pseudo.basis {
            out.v3(basis);
        }
        for interval in pseudo.intervals {
            out.d(interval.min);
            out.d(interval.max);
        }
    }
    for parent in &cache.parents {
        span(out, parent);
    }
    node(out, &cache, cache.top);
}

#[test]
fn curve_cache_matches_smlib() {
    let report = run_gate(
        include_str!("fixtures/curve-cache/cases.txt"),
        |_| 16,
        evaluate,
    );
    report.assert_passed(&["tangent", "touching", "singular"]);
}
