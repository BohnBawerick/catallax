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

/// Regressions from the step 4 part 2 review.
mod review {
    use catallax_trace::Fault;
    use catallax_trace::curve::evaluate_geometric;
    use catallax_trace::curve_cache::CurveCache;
    use catallax_trace::spline::{BSplineCurve, ControlPoint, NO_WEIGHT};

    fn point(x: f64, y: f64, z: f64) -> ControlPoint {
        ControlPoint {
            coordinates: [x, y, z, NO_WEIGHT],
        }
    }

    /// The review's unclamped quadratic. SMLib's Bezier extraction assumes
    /// clamped knots and misplaces the spans of this curve (its span over
    /// `[2, 3]` gives `(1.5, 0.25, 0)` at 2.5 where the curve is at
    /// `(1.75, 0.25, 0)`), so the cache refuses it rather than build spans
    /// that leave the curve.
    #[test]
    fn cache_spans_follow_the_curve_or_refuse() {
        let curve = BSplineCurve::new(
            2,
            vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
            vec![
                point(0.0, 0.0, 0.0),
                point(2.0, 0.0, 0.0),
                point(2.0, 2.0, 0.0),
            ],
        )
        .expect("accepted by the constructor");
        match CurveCache::for_curve(&curve, 5.0e-6) {
            Err(fault) => assert_eq!(fault, Fault::InvalidInput),
            Ok(cache) => {
                for (span, bezier) in cache.spans.iter().zip(&cache.beziers) {
                    let t = span.interval.evaluate(0.5);
                    let on_curve = curve.evaluate_point(t).expect("curve");
                    let on_span = bezier.evaluate_point(t).expect("span");
                    assert!(
                        on_curve.distance(on_span) < 1e-12,
                        "span {:?} at {t}: {on_span:?} vs curve {on_curve:?}",
                        span.interval
                    );
                }
            }
        }
    }

    /// `EvaluateGeometric` builds at most the third Frenet vector; larger
    /// counts are refused instead of overrunning the four-slot frame.
    #[test]
    fn geometric_evaluation_refuses_counts_above_three() {
        let line = BSplineCurve::new(
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![point(0.0, 0.0, 0.0), point(1.0, 0.0, 0.0)],
        )
        .expect("line");
        for count in 0..=1 {
            assert!(
                evaluate_geometric(&line, 0.5, count, false).is_ok(),
                "{count}"
            );
        }
        // A line has no normal: counts 2 and 3 fail as SMLib does, without
        // a refusal for the count itself.
        for count in 2..=3 {
            assert_ne!(
                evaluate_geometric(&line, 0.5, count, false),
                Err(Fault::InvalidInput),
                "{count}"
            );
        }
        assert_eq!(
            evaluate_geometric(&line, 0.5, 4, false),
            Err(Fault::InvalidInput)
        );
    }
}
