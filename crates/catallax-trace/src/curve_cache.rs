// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: the cache owns its Bezier spans, span data and
// search-tree arena instead of SMLib's memory managers and cache manager;
// statuses become `Fault`; the approximation tolerance is the caller's
// argument instead of SMLib's global tolerance state.
// Translated from SMLib SmCurveCache.cpp (Tessellate, BuildTree,
// OutputSegment, sm_ComputeConstants, sm_ComputeParentData,
// sm_CreateBezierSegmentBlock), SmNurbsCrv.cpp (sm_CreateBezSegments,
// sm_CrvTrim) and SmCurve.cpp (CacheMakeOrValidate); sm_CrvTrim also adapts
// rust-brep crates/kernel/src/curve/sm_nurbs.rs, commit
// 4c83df59d777ccfcb14f441b1c8d984d1b5d93ab.
//! SMLib's curve cache: a curve cut into Bezier spans small enough for the
//! global search, with a binary tree of boxes over them (step 4 of the M25
//! port order).

use crate::Fault;
use crate::curve::{
    CurveDistanceExtrema, approximate_length, calculate_continuities, natural_knot_count, polygon,
};
use crate::extent::{Extent1, Extent3};
use crate::geometry::{Continuity, line_closest_point};
use crate::global_solver::{AuxData, SearchTrees, TreeNode};
use crate::local_solve::{EvalFunction, LocalSolve1d};
use crate::pseudo_box::PseudoBox;
use crate::scalar::{EFF_ZERO, EFF_ZERO_SQ, EFF_ZERO_SQRT, PI, acos, is_zero, sqrt};
use crate::solution::NodeId;
use crate::spline::{BSplineCurve, ControlPoint, combine};
use crate::vector::Vector3;

/// `SmBezierAux1d`: what a node of the curve tree knows about the spans
/// below it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpanData {
    /// Index of the leftmost leaf span below.
    pub left: usize,
    /// Index of the rightmost leaf span below.
    pub right: usize,
    /// Continuity with the left and right neighbouring spans.
    pub continuities: [Continuity; 2],
    /// Control-polygon base-line length.
    pub max_length: f64,
    /// Largest squared distance of a control point from the base line.
    pub max_chord_height_squared: f64,
    /// Accumulated control-polygon turning angle, in degrees.
    pub max_turning_angle_deg: f64,
    /// Parameter interval covered.
    pub interval: Extent1,
}

/// What a tree node's data refers to.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NodeData {
    /// A leaf: the index of its Bezier span.
    Leaf(usize),
    /// An internal node with propagated data: the index into the parents.
    Parent(usize),
    /// An internal node without data.
    None,
}

/// The tessellation controls of an `SmCurveCache`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CacheSettings {
    /// Keep every internal node (`bMakeFullTree`).
    pub make_full_tree: bool,
    /// Largest control-point distance from a span's base line; 0 ignores.
    pub chord_height_tolerance: f64,
    /// Largest accumulated turning angle of a span, in degrees; 0 ignores.
    pub angle_tolerance_deg: f64,
    /// Largest second box dimension of a span's polygon; 0 ignores.
    pub off_axis_tolerance: f64,
    /// Fewest spans; 0 ignores.
    pub min_segments: usize,
    /// Propagate span data up to internal nodes.
    pub compute_aux: bool,
    /// Largest span base-line length; 0 ignores.
    pub max_3d_distance: f64,
    /// Smallest span interval as a fraction of the whole; 0 ignores.
    pub min_param_ratio: f64,
}

/// `SmCurveCache` over a B-spline curve.
#[derive(Clone, Debug)]
pub struct CurveCache {
    /// The parameter interval covered.
    pub interval: Extent1,
    /// The continuity class at each unique knot, when the full tree is kept.
    pub continuities: Vec<Continuity>,
    /// Leaf span data, in parameter order.
    pub spans: Vec<SpanData>,
    /// The pseudo box of each leaf span.
    pub pseudo_boxes: Vec<PseudoBox>,
    /// The Bezier curve of each leaf span.
    pub beziers: Vec<BSplineCurve>,
    /// Internal node data, in creation order.
    pub parents: Vec<SpanData>,
    /// The tree: leaves first, then internal nodes as they are built.
    pub trees: SearchTrees,
    /// The data of each tree node, by node identity.
    pub node_data: Vec<NodeData>,
    /// The root of the tree.
    pub top: NodeId,
}

/// `N_BasisFindSpanAndMult` from the left: the span index and, when `u` is
/// a knot, its multiplicity.
fn span_and_multiplicity_left(curve: &BSplineCurve, u: f64) -> Result<(usize, usize), Fault> {
    let knots = curve.knots();
    let p = curve.degree();
    let span = curve.find_span(u, true)?;
    if u == knots[knots.len() - 1 - p] {
        return Ok((span, p + 1));
    }
    let mut k = span as isize;
    let mut multiplicity = 0;
    while k >= 0 && u == knots[k as usize] {
        multiplicity += 1;
        k -= 1;
    }
    Ok((span, multiplicity))
}

/// `sm_CrvTrim`: the curve restricted to `[ul, ur]` by knot insertion.
pub(crate) fn trim(curve: &BSplineCurve, ul: f64, ur: f64) -> Result<BSplineCurve, Fault> {
    if ur <= ul {
        return Err(Fault::InvalidInput);
    }
    let p = curve.degree() as isize;
    let knots = curve.knots();
    let m = knots.len() as isize - 1;
    let (spl, mll) = span_and_multiplicity_left(curve, ul)?;
    let (mut spr, mut mlr) = span_and_multiplicity_left(curve, ur)?;
    let (spl, mll) = (spl as isize, mll as isize);
    let is = spl - p;
    if ur == knots[(m - p) as usize] {
        spr = m as usize;
        mlr = p as usize + 1;
    }
    let (spr, mlr) = (spr as isize, mlr as isize);
    let ie = spr - mlr;
    let n = ie - is;
    if n < 0 {
        return Err(Fault::Degenerate);
    }
    let pw = curve.control_points();
    let mut qw: Vec<ControlPoint> = (is..=ie).map(|i| pw[i as usize]).collect();
    let u = |i: isize| knots[i as usize];
    let ll = spl - p;
    for i in 1..=p - mll {
        for j in 0..=p - i - mll {
            let left = u(ll + i + j);
            let alpha = (ul - left) / (u(spl + j + 1) - left);
            qw[j as usize] = combine(alpha, qw[j as usize + 1], 1.0 - alpha, qw[j as usize]);
        }
    }
    let lr = spr - p;
    let lk = n - p + mlr;
    for i in 1..=p - mlr {
        for j in (0..=p - i - mlr).rev() {
            let k = (lk + i + j) as usize;
            let mut left = u(lr + i + j);
            if left < ul {
                left = ul;
            }
            let alpha = (ur - left) / (u(spr + j + 1) - left);
            qw[k] = combine(alpha, qw[k], 1.0 - alpha, qw[k - 1]);
        }
    }
    let mut uq = vec![ul; p as usize + 1];
    uq.extend((spl + 1..=spr - mlr).map(u));
    uq.extend(std::iter::repeat_n(ur, p as usize + 1));
    BSplineCurve::new(p as usize, uq, qw)
}

/// Whether both ends of the knot vector repeat `degree + 1` times. SMLib's
/// Bezier extraction (`sm_CreateBezSegments`) and trim (`sm_CrvTrim`) assume
/// it, as SMLib's curves have it; the cache refuses other knot vectors the
/// constructor accepts.
pub(crate) fn is_clamped(curve: &BSplineCurve) -> bool {
    let (knots, p) = (curve.knots(), curve.degree());
    let n = knots.len();
    n >= 2 * (p + 1)
        && knots[..=p].iter().all(|&k| k == knots[0])
        && knots[n - 1 - p..].iter().all(|&k| k == knots[n - 1])
}

/// `N_BasisGetSpanCount`: the nonzero spans of the active knot range.
fn span_count(curve: &BSplineCurve) -> usize {
    let (knots, p) = (curve.knots(), curve.degree());
    (p..knots.len() - 1 - p)
        .filter(|&i| knots[i] != knots[i + 1])
        .count()
}

/// `sm_CreateBezSegments`: one Bezier curve per nonzero span, by knot
/// insertion.
fn bezier_segments(curve: &BSplineCurve) -> Result<Vec<BSplineCurve>, Fault> {
    let spans = span_count(curve);
    if spans < 1 {
        return Err(Fault::Degenerate);
    }
    let p = curve.degree();
    let up = curve.knots();
    let pw = curve.control_points();
    let m = up.len() - 1;
    let mut qws: Vec<Vec<ControlPoint>> = vec![vec![pw[0]; p + 1]; spans];
    let mut uqs: Vec<Vec<f64>> = vec![vec![0.0; 2 * p + 2]; spans];
    let (mut is, mut ie) = (p, p + 1);
    while ie < m && up[ie - 1] == up[ie] {
        is += 1;
        ie += 1;
    }
    for i in 0..=p {
        qws[0][i] = pw[i + is - p];
    }
    let mut iq = 0usize;
    let mut first = true;
    let mut alphas = vec![0.0; p];
    let mut betas = vec![0.0; p];
    while ie < m {
        if first {
            first = false;
        } else {
            iq += 1;
        }
        if iq >= spans {
            break;
        }
        let i0 = ie;
        while ie < m && up[ie] == up[ie + 1] {
            ie += 1;
        }
        let multiplicity = ie - i0 + 1;
        let r = p as isize - multiplicity as isize;
        if multiplicity < p {
            let numerator = up[ie] - up[is];
            for i in (multiplicity + 1..=p).rev() {
                alphas[i - multiplicity - 1] = numerator / (up[is + i] - up[is]);
                betas[i - multiplicity - 1] = 1.0 - alphas[i - multiplicity - 1];
            }
            for i in 1..=r as usize {
                let s = multiplicity + i;
                let save = r as usize - i;
                for j in (s..=p).rev() {
                    qws[iq][j] = combine(alphas[j - s], qws[iq][j], betas[j - s], qws[iq][j - 1]);
                }
                if ie < m && iq + 1 < spans {
                    qws[iq + 1][save] = qws[iq][p];
                }
            }
        }
        for i in 0..=p {
            uqs[iq][i] = up[is];
            uqs[iq][i + p + 1] = up[ie];
        }
        if ie < m && iq + 1 < spans {
            for i in r.max(0) as usize..=p {
                qws[iq + 1][i] = pw[ie - p + i];
            }
        }
        is = ie;
        ie += 1;
    }
    qws.into_iter()
        .zip(uqs)
        .map(|(points, knots)| BSplineCurve::new(p, knots, points))
        .collect()
}

/// The span measures of `sm_ComputeConstants`.
#[derive(Clone, Copy, Debug, Default)]
struct SpanMeasures {
    max_length: f64,
    max_chord_height_squared: f64,
    max_turning_angle_deg: f64,
    off_axis_distance: f64,
    curve_distance: f64,
    null_curve: bool,
}

/// `sm_ComputeConstants`: control-polygon measures of `bezier`, stopping at
/// the first exceeded limit, then the sampled distance from `original`.
fn compute_constants(
    bezier: &BSplineCurve,
    original: &BSplineCurve,
    interval: Extent1,
    settings: &CacheSettings,
    approx_tolerance: f64,
) -> Result<SpanMeasures, Fault> {
    let mut out = SpanMeasures::default();
    let chord_tolerance_squared = settings.chord_height_tolerance * settings.chord_height_tolerance;
    let angle_tolerance_rad = settings.angle_tolerance_deg * PI / 180.0;
    let poly = polygon(bezier, bezier.natural_interval())?;
    if poly.len() < 2 {
        return Err(Fault::Degenerate);
    }
    let size = poly.len();
    let line_point = poly[0];
    let mut line_vector = poly[size - 1] - line_point;
    out.max_length = line_vector.length();
    if settings.max_3d_distance != 0.0 && out.max_length > settings.max_3d_distance {
        return Ok(out);
    }
    let mut try_index = size - 2;
    while line_vector.length_squared() < EFF_ZERO * EFF_ZERO {
        if try_index == 0 {
            out.null_curve = true;
            return Ok(out);
        }
        line_vector = poly[try_index] - line_point;
        try_index -= 1;
    }
    let mut last_vector = Vector3::new(0.0, 0.0, 0.0);
    let mut total_angle = 0.0;
    let mut max_segment_chord = 0.0;
    let mut bbox = Extent3::default();
    bbox.add_point(poly[0]);
    for i in 1..size {
        let (current, next) = (poly[i - 1], poly[i]);
        bbox.add_point(next);
        if settings.off_axis_tolerance != 0.0 {
            let d = bbox.max - bbox.min;
            let over = [d.x, d.y, d.z]
                .into_iter()
                .filter(|&v| v > settings.off_axis_tolerance)
                .count();
            if over > 1 {
                out.off_axis_distance = if d.x >= d.y && d.x >= d.z {
                    crate::scalar::max(d.y, d.z)
                } else if d.y >= d.x && d.y >= d.z {
                    crate::scalar::max(d.x, d.z)
                } else {
                    crate::scalar::max(d.x, d.y)
                };
                return Ok(out);
            }
        }
        let parameter = line_closest_point(line_point, line_vector, next);
        let on_line = line_point + line_vector * parameter;
        let chord = (next - on_line).length_squared();
        if chord > max_segment_chord {
            max_segment_chord = chord;
        }
        if chord_tolerance_squared != 0.0 && max_segment_chord > chord_tolerance_squared {
            out.max_chord_height_squared = max_segment_chord;
            return Ok(out);
        }
        let mut step = next - current;
        let distance = step.length();
        let scale = 1.0 + next.max_dimension();
        if distance < EFF_ZERO * scale {
            continue;
        }
        step = step.divided(distance);
        if last_vector.length_squared() > EFF_ZERO_SQ {
            let cosine = step.dot(last_vector);
            total_angle += if cosine >= 1.0 {
                0.0
            } else if cosine <= -1.0 {
                PI
            } else {
                acos(cosine)
            };
        }
        if total_angle != 0.0 && total_angle > angle_tolerance_rad {
            out.max_turning_angle_deg = total_angle * 180.0 / PI;
            return Ok(out);
        }
        last_vector = step;
    }
    out.max_chord_height_squared = max_segment_chord;
    out.max_turning_angle_deg = total_angle * 180.0 / PI;

    let tolerance = 10.0 * EFF_ZERO;
    let mut solver = LocalSolve1d::with_interval(
        CurveDistanceExtrema {
            first: bezier,
            second: original,
            tolerance,
        },
        interval,
        false,
    );
    let steps = 20;
    let approx_squared = approx_tolerance * approx_tolerance;
    let (mut f, mut f1, mut found) = (0.0, 0.0, false);
    let mut previous_t = interval.min;
    let mut previous_f = 0.0;
    let _ = solver
        .evaluator
        .evaluate(previous_t, &mut previous_f, &mut f1, &mut found, true);
    let (mut point1, mut point2) = (Vector3::default(), Vector3::default());
    for ii in 1..=steps {
        let parameter = interval.evaluate(ii as f64 / f64::from(steps));
        let _ = solver
            .evaluator
            .evaluate(parameter, &mut f, &mut f1, &mut found, true);
        let compare = if found {
            true
        } else if f * previous_f < 0.0 {
            let guess = previous_t - previous_f * (parameter - previous_t) / (f - previous_f);
            solver.solve_it(guess, tolerance).is_ok()
        } else {
            false
        };
        if compare {
            // SMLib measures at the sample parameter, not at the root found.
            if let Ok(p) = bezier.evaluate_point(parameter) {
                point1 = p;
            }
            if let Ok(p) = original.evaluate_point(parameter) {
                point2 = p;
            }
            let distance_squared = point1.distance_squared(point2);
            if distance_squared > approx_squared {
                out.curve_distance = sqrt(distance_squared);
                break;
            }
        }
        previous_t = parameter;
        previous_f = f;
    }
    Ok(out)
}

/// Whether span measures pass the tessellation tests.
fn passes(
    measures: &SpanMeasures,
    length: f64,
    max_step: f64,
    settings: &CacheSettings,
    approx_tolerance: f64,
) -> bool {
    let chord_squared = settings.chord_height_tolerance * settings.chord_height_tolerance;
    !(measures.max_turning_angle_deg > 90.0
        || length > max_step
        || (settings.chord_height_tolerance != 0.0
            && measures.max_chord_height_squared > chord_squared)
        || (settings.angle_tolerance_deg != 0.0
            && measures.max_turning_angle_deg > settings.angle_tolerance_deg)
        || (settings.off_axis_tolerance != 0.0
            && measures.off_axis_distance > settings.off_axis_tolerance)
        || (settings.max_3d_distance != 0.0 && measures.max_length > settings.max_3d_distance)
        || measures.curve_distance > approx_tolerance)
}

impl CurveCache {
    /// `SmCurve::CacheMakeOrValidate`: the cache SMLib builds for a curve:
    /// the full tree with propagated data, a chord tolerance of the
    /// approximate length over `min(250, 30 * spans)`, a 15 degree angle
    /// tolerance and a minimum span ratio of 0.0001.
    pub fn for_curve(curve: &BSplineCurve, approx_tolerance: f64) -> Result<Self, Fault> {
        let interval = curve.natural_interval();
        let length = approximate_length(curve, interval, 5);
        let spans = natural_knot_count(curve) as f64 - 1.0;
        let factor = crate::scalar::min(250.0, 30.0 * spans);
        let settings = CacheSettings {
            make_full_tree: true,
            chord_height_tolerance: length / factor,
            angle_tolerance_deg: 15.0,
            off_axis_tolerance: 0.0,
            min_segments: 0,
            compute_aux: true,
            max_3d_distance: 0.0,
            min_param_ratio: 0.0001,
        };
        Self::tessellate(curve, interval, &settings, approx_tolerance)
    }

    /// `SmCurveCache::Tessellate` and `BuildTree`: Bezier spans of `curve`
    /// over `interval`, halved until each passes the tests of `settings`,
    /// and the tree over them.
    pub fn tessellate(
        curve: &BSplineCurve,
        interval: Extent1,
        settings: &CacheSettings,
        approx_tolerance: f64,
    ) -> Result<Self, Fault> {
        if is_zero(settings.chord_height_tolerance) && is_zero(settings.angle_tolerance_deg) {
            return Err(Fault::InvalidInput);
        }
        if !is_clamped(curve) {
            return Err(Fault::InvalidInput);
        }
        let continuities = if settings.make_full_tree {
            calculate_continuities(curve).1
        } else {
            Vec::new()
        };
        let continuity = |i: usize| {
            if settings.make_full_tree {
                continuities
                    .get(i)
                    .copied()
                    .unwrap_or(Continuity::Undefined)
            } else {
                Continuity::C1C2
            }
        };
        let beziers = bezier_segments(curve)?;
        let max_step = if settings.min_segments > 1 {
            (interval.max - interval.min) / settings.min_segments as f64
        } else {
            interval.max - interval.min
        };
        let mut spans: Vec<SpanData> = Vec::new();
        let mut leaves: Vec<BSplineCurve> = Vec::new();
        // `OutputSegment`: a leaf is its own leftmost and rightmost span.
        fn output(
            spans: &mut Vec<SpanData>,
            leaves: &mut Vec<BSplineCurve>,
            bezier: BSplineCurve,
            measures: &SpanMeasures,
            continuities: [Continuity; 2],
        ) {
            let index = spans.len();
            let knots = bezier.knots();
            spans.push(SpanData {
                left: index,
                right: index,
                continuities,
                max_length: measures.max_length,
                max_chord_height_squared: measures.max_chord_height_squared,
                max_turning_angle_deg: measures.max_turning_angle_deg,
                interval: Extent1::new(knots[0], knots[knots.len() - 1]),
            });
            leaves.push(bezier);
        }
        let count = beziers.len();
        for (i, bezier) in beziers.into_iter().enumerate() {
            let knots = bezier.knots();
            let mut span_interval = Extent1::new(knots[0], knots[knots.len() - 1]);
            if interval.are_disjoint(span_interval, 0.0) {
                continue;
            }
            let mut bezier = bezier;
            if !span_interval.is_contained_by(interval, EFF_ZERO) {
                span_interval = Extent1 {
                    min: crate::scalar::max(interval.min, span_interval.min),
                    max: crate::scalar::min(interval.max, span_interval.max),
                };
                let scaled = EFF_ZERO * (1.0 + span_interval.max.abs() + span_interval.min.abs());
                if span_interval.max - span_interval.min < scaled {
                    continue;
                }
                bezier = trim(&bezier, span_interval.min, span_interval.max)?;
            }
            let measures =
                compute_constants(&bezier, curve, span_interval, settings, approx_tolerance)?;
            if measures.null_curve && (!spans.is_empty() || i + 1 < count) {
                continue;
            }
            let mut pass = passes(
                &measures,
                span_interval.length(),
                max_step,
                settings,
                approx_tolerance,
            );
            if (settings.min_param_ratio != 0.0
                && span_interval.length() / 2.0 < interval.length() * settings.min_param_ratio)
                || measures.null_curve
            {
                pass = true;
            }
            if pass {
                output(
                    &mut spans,
                    &mut leaves,
                    bezier,
                    &measures,
                    [continuity(i), continuity(i + 1)],
                );
                continue;
            }
            let mut stack: Vec<BSplineCurve> = Vec::new();
            let [low, high] = bezier.split(span_interval.mid())?;
            stack.push(high);
            stack.push(low);
            let mut first_in_segment = true;
            while let Some(piece) = stack.pop() {
                let piece_interval = piece.natural_interval();
                let measures =
                    compute_constants(&piece, curve, piece_interval, settings, approx_tolerance)?;
                if measures.null_curve {
                    continue;
                }
                let mut pass = passes(
                    &measures,
                    piece_interval.length(),
                    max_step,
                    settings,
                    approx_tolerance,
                );
                if settings.min_param_ratio != 0.0
                    && piece_interval.length() / 2.0 < interval.length() * settings.min_param_ratio
                {
                    pass = true;
                }
                if piece_interval.length() < EFF_ZERO_SQRT {
                    pass = true;
                }
                if pass {
                    let left = if settings.make_full_tree && first_in_segment {
                        first_in_segment = false;
                        continuity(i)
                    } else {
                        Continuity::C1C2
                    };
                    let right = if settings.make_full_tree && stack.is_empty() {
                        continuity(i + 1)
                    } else {
                        Continuity::C1C2
                    };
                    output(&mut spans, &mut leaves, piece, &measures, [left, right]);
                    continue;
                }
                let [low, high] = piece.split(piece_interval.mid())?;
                stack.push(high);
                stack.push(low);
            }
        }
        Self::build_tree(curve, interval, settings, continuities, spans, leaves)
    }

    /// `BuildTree`: leaves with their boxes, then pairs joined level by
    /// level, an odd last node carried up.
    fn build_tree(
        curve: &BSplineCurve,
        interval: Extent1,
        settings: &CacheSettings,
        continuities: Vec<Continuity>,
        spans: Vec<SpanData>,
        beziers: Vec<BSplineCurve>,
    ) -> Result<Self, Fault> {
        if spans.is_empty() {
            return Err(Fault::Degenerate);
        }
        let chord_squared = settings.chord_height_tolerance * settings.chord_height_tolerance;
        let max_step = if settings.min_segments > 1 {
            (interval.max - interval.min) / settings.min_segments as f64
        } else {
            interval.max - interval.min
        };
        let mut trees = SearchTrees::new();
        let mut node_data = Vec::new();
        let mut pseudo_boxes = Vec::new();
        let mut current = Vec::new();
        for (i, bezier) in beziers.iter().enumerate() {
            let (bbox, pseudo) = bezier.bounding_boxes(bezier.natural_interval())?;
            pseudo_boxes.push(pseudo);
            current.push(trees.add(TreeNode {
                bbox,
                children: None,
                aux: AuxData::BezierCurve,
            }));
            node_data.push(NodeData::Leaf(i));
        }
        let mut parents: Vec<SpanData> = Vec::new();
        let span_of = |data: NodeData, parents: &[SpanData]| match data {
            NodeData::Leaf(i) => spans[i],
            NodeData::Parent(i) => parents[i],
            NodeData::None => spans[0],
        };
        while current.len() > 1 {
            let mut next = Vec::new();
            for &[a, b] in current.as_chunks::<2>().0 {
                let bbox = trees
                    .get(a)
                    .expect("node")
                    .bbox
                    .union(trees.get(b).expect("node").bbox);
                let mut node = TreeNode {
                    bbox,
                    children: Some([a, b]),
                    aux: AuxData::None,
                };
                let mut data = NodeData::None;
                if settings.compute_aux {
                    let left = span_of(node_data[a.0 as usize], &parents);
                    let right = span_of(node_data[b.0 as usize], &parents);
                    let parent = parent_data(curve, &spans, &left, &right);
                    node.aux = AuxData::CurveData;
                    data = NodeData::Parent(parents.len());
                    if !settings.make_full_tree {
                        let length = parent.interval.length();
                        let mut pass = !(parent.max_turning_angle_deg > 90.0
                            || (settings.chord_height_tolerance != 0.0
                                && parent.max_chord_height_squared > chord_squared)
                            || (settings.angle_tolerance_deg != 0.0
                                && parent.max_turning_angle_deg > settings.angle_tolerance_deg)
                            || (settings.max_3d_distance != 0.0
                                && parent.max_length > settings.max_3d_distance)
                            || length > max_step);
                        if settings.min_param_ratio != 0.0
                            && length < interval.length() * settings.min_param_ratio
                        {
                            pass = true;
                        }
                        if pass {
                            node.children = None;
                        }
                    }
                    parents.push(parent);
                }
                next.push(trees.add(node));
                node_data.push(data);
            }
            if current.len() % 2 == 1 {
                next.push(current[current.len() - 1]);
            }
            current = next;
        }
        Ok(Self {
            interval,
            continuities,
            spans,
            pseudo_boxes,
            beziers,
            parents,
            trees,
            node_data,
            top: current[0],
        })
    }
}

/// `sm_ComputeParentData` with the full tree kept: the children's data
/// joined, plus the turning angle across a join that is not tangent.
fn parent_data(
    curve: &BSplineCurve,
    spans: &[SpanData],
    left: &SpanData,
    right: &SpanData,
) -> SpanData {
    let mut parent = SpanData {
        left: left.left,
        right: right.right,
        continuities: [left.continuities[0], right.continuities[1]],
        max_length: left.max_length + right.max_length,
        max_chord_height_squared: 0.0,
        max_turning_angle_deg: left.max_turning_angle_deg + right.max_turning_angle_deg,
        interval: Extent1 {
            min: left.interval.min,
            max: right.interval.max,
        },
    };
    if left.continuities[1] < Continuity::G1 {
        let end = spans[left.right].interval.evaluate(1.0);
        let start = spans[right.left].interval.evaluate(0.0);
        if let (Ok(a), Ok(b)) = (
            curve.evaluate(end, 1, false, true),
            curve.evaluate(start, 1, true, true),
        ) && a[1].length_squared() > EFF_ZERO_SQ
            && b[1].length_squared() > EFF_ZERO_SQ
            && let Ok(angle) = a[1].angle_between(b[1])
        {
            parent.max_turning_angle_deg += angle * 180.0 / PI;
        }
    }
    parent
}
