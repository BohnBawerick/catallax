//! Display sampling of admitted native carriers. No exact Model is fabricated.
//!
//! The carriers are sampled in binary64 through their display variants: the
//! preview takes thousands of samples per face and claims no bound on any of
//! them, and an enclosing evaluation of each costs a body's whole allowance.
use std::{collections::BTreeMap, io::Write};

use catallax_import::{
    Curve,
    geometry_source::Topology,
    tolerant::{
        TolerantFace, TolerantUse,
        admission::{NativeTrim, TolerantImport},
        native::{DisplayCurve, DisplaySurface, NativeLocation, NativeSurface},
    },
};
use catallax_numbers::{Decider, Rational, Refusal, refund_on_refusal};
use spade::Triangulation;

use super::imported::{
    Cdt, FaceOutput, GRID_TRIM_REFUSAL, Node, SourceChart, SurfaceSamples, finish_face, grid_face,
    insert, push, seed_rectangles,
};
use crate::{PreviewEdge, PreviewMesh, PreviewOptions};

type P = [f64; 3];
type Uv = [f64; 2];
fn failure(detail: impl Into<String>) -> Refusal {
    Refusal::new("mesh.tolerant-preview-unresolved", detail)
}
fn source_error(error: catallax_import::tolerant::Error) -> Refusal {
    match error {
        catallax_import::tolerant::Error::Numerical(error) => error,
        other => failure(format!("{other:?}")),
    }
}
fn distance(a: P, b: P) -> f64 {
    (0..3).fold(0.0, |norm, axis| norm.hypot(a[axis] - b[axis]))
}
fn point(value: &[Rational; 3]) -> Result<P, Refusal> {
    let point = value.each_ref().map(Rational::to_f64_nearest);
    if point.iter().all(|coordinate| coordinate.is_finite()) {
        Ok(point)
    } else {
        Err(failure("vertex exceeds binary64 range"))
    }
}

/// For each parameter axis along which a native spline closes on itself
/// without declaring a period, the axis's lower limit and its width.
type Seams = [Option<[f64; 2]>; 2];

enum Carrier {
    /// The retained native surface, made for display, and its seams.
    Native(Box<DisplaySurface>, Seams),
    Inferred(SourceChart),
}
impl Carrier {
    /// A parameter past either end of a closed axis is read one turn back.
    /// Parameters inside the axis are left exactly as they are.
    fn wrap(&self, mut uv: Uv) -> Uv {
        if let Self::Native(_, seams) = self {
            for (value, seam) in uv.iter_mut().zip(seams) {
                if let Some([low, width]) = *seam
                    && !(low..=low + width).contains(value)
                {
                    *value = low + (*value - low).rem_euclid(width);
                }
            }
        }
        uv
    }
}
impl SurfaceSamples for Carrier {
    fn position(&self, uv: Uv, dec: &mut Decider) -> Result<P, Refusal> {
        match self {
            Self::Native(surface, _) => surface
                .sample(surface.location(self.wrap(uv)).map_err(source_error)?, dec)
                .map(|sample| sample.point_mm)
                .map_err(source_error),
            Self::Inferred(surface) => SurfaceSamples::position(surface, uv, dec),
        }
    }
    fn normal(&self, uv: Uv, dec: &mut Decider) -> Result<P, Refusal> {
        let Self::Native(surface, _) = self else {
            let Self::Inferred(surface) = self else {
                unreachable!()
            };
            return SurfaceSamples::normal(surface, uv, dec);
        };
        let location = surface.location(self.wrap(uv)).map_err(source_error)?;
        surface
            .sampled_normal(location, dec)
            .map_err(source_error)?
            .ok_or_else(|| failure("native sampled normal collapsed"))
    }
}

enum EdgeCarrier {
    /// The retained native curve, made for display, and its end parameters.
    Native(Box<DisplayCurve>, [f64; 2]),
    Line([P; 2]),
}
impl EdgeCarrier {
    fn sample(&self, fraction: f64, dec: &mut Decider) -> Result<P, Refusal> {
        match self {
            Self::Native(curve, range) => {
                let parameter = range[0] + fraction * (range[1] - range[0]);
                curve
                    .sample(curve.location(parameter).map_err(source_error)?, dec)
                    .map(|sample| sample.point_mm)
                    .map_err(source_error)
            }
            Self::Line(ends) => Ok(std::array::from_fn(|axis| {
                ends[0][axis] + fraction * (ends[1][axis] - ends[0][axis])
            })),
        }
    }
}

fn edge_samples(
    carrier: &EdgeCarrier,
    ends: [usize; 2],
    target: f64,
    mesh: &mut PreviewMesh<String, String>,
    options: &PreviewOptions,
    dec: &mut Decider,
) -> Result<Vec<(f64, usize)>, Refusal> {
    let mut output = vec![(0.0, ends[0])];
    let mut stack = vec![(0.0, 1.0, ends, 0)];
    while let Some((a, b, indices, depth)) = stack.pop() {
        dec.charge_work(1)?;
        let mid = (a + b) * 0.5;
        let actual = carrier.sample(mid, dec)?;
        let linear = std::array::from_fn(|axis| {
            (mesh.positions[indices[0]][axis] + mesh.positions[indices[1]][axis]) * 0.5
        });
        if distance(actual, linear) <= target * 0.5 {
            output.push((b, indices[1]));
        } else {
            if depth == 24 || mid <= a || mid >= b {
                return Err(failure("edge subdivision reached its sampling limit"));
            }
            let index = push(mesh, actual, options)?;
            stack.push((mid, b, [index, indices[1]], depth + 1));
            stack.push((a, mid, [indices[0], index], depth + 1));
        }
    }
    Ok(output)
}

fn grid_edge(
    carrier: &EdgeCarrier,
    ends: [usize; 2],
    divisions: u32,
    mesh: &mut PreviewMesh<String, String>,
    options: &PreviewOptions,
    dec: &mut Decider,
) -> Result<Vec<(f64, usize)>, Refusal> {
    let divisions = match carrier {
        EdgeCarrier::Line(_) => 1,
        EdgeCarrier::Native(..) if ends[0] == ends[1] => divisions.max(4),
        EdgeCarrier::Native(..) => divisions.max(2),
    };
    dec.charge_traversal(u64::from(divisions))?;
    let mut output = Vec::with_capacity(divisions as usize + 1);
    output.push((0.0, ends[0]));
    for i in 1..divisions {
        let fraction = f64::from(i) / f64::from(divisions);
        let index = push(mesh, carrier.sample(fraction, dec)?, options)?;
        output.push((fraction, index));
    }
    output.push((1.0, ends[1]));
    Ok(output)
}

fn trim_seed(trim: &NativeTrim, fraction: f64) -> NativeLocation {
    let next = trim
        .samples
        .partition_point(|sample| sample.fraction < fraction)
        .min(trim.samples.len() - 1)
        .max(1);
    let first = next - 1;
    let width = trim.samples[next].fraction - trim.samples[first].fraction;
    let weight = ((fraction - trim.samples[first].fraction) / width).clamp(0.0, 1.0);
    let a = trim.samples[first].surface;
    let b = trim.samples[next].surface;
    NativeLocation {
        chart: a.chart,
        parameters: std::array::from_fn(|axis| {
            a.parameters[axis] + weight * (b.parameters[axis] - a.parameters[axis])
        }),
    }
}

fn domain_loop(
    carrier: &Carrier,
    bounds: [[f64; 2]; 2],
    mesh: &mut PreviewMesh<String, String>,
    options: &PreviewOptions,
    dec: &mut Decider,
) -> Result<Vec<Node>, Refusal> {
    let corners = [
        [bounds[0][0], bounds[1][0]],
        [bounds[0][1], bounds[1][0]],
        [bounds[0][1], bounds[1][1]],
        [bounds[0][0], bounds[1][1]],
    ];
    let mut nodes = Vec::new();
    for side in 0..4 {
        let a = corners[side];
        let b = corners[(side + 1) % 4];
        let start = carrier.position(a, dec)?;
        let end = carrier.position(b, dec)?;
        let mut stack = vec![(a, b, start, end, 0)];
        while let Some((a, b, start, end, depth)) = stack.pop() {
            dec.charge_work(1)?;
            let mid = std::array::from_fn(|axis| (a[axis] + b[axis]) * 0.5);
            let actual = carrier.position(mid, dec)?;
            let linear = std::array::from_fn(|axis| (start[axis] + end[axis]) * 0.5);
            if distance(actual, linear) <= options.chord_tolerance.to_f64_nearest() * 0.5 {
                nodes.push(Node {
                    uv: a,
                    mesh: push(mesh, start, options)?,
                });
            } else {
                if depth == 24
                    || mid.map(f64::to_bits) == a.map(f64::to_bits)
                    || mid.map(f64::to_bits) == b.map(f64::to_bits)
                {
                    return Err(failure(
                        "natural seam subdivision reached its sampling limit",
                    ));
                }
                stack.push((mid, b, actual, end, depth + 1));
                stack.push((a, mid, start, actual, depth + 1));
            }
        }
    }
    Ok(nodes)
}

#[allow(clippy::too_many_arguments)]
fn cut_to_pole(
    body: &mut PreviewMesh<String, String>,
    rim: Vec<Node>,
    source: &NativeSurface,
    (periods, steps): ([Option<f64>; 2], u32),
    carrier: &Carrier,
    closure: (bool, Option<(usize, f64)>),
    options: &PreviewOptions,
    dec: &mut Decider,
) -> Result<Vec<Node>, Refusal> {
    let failure = || failure("single periodic rim has no finite collapsed pole");
    let (same_sense, vertex_loop) = closure;
    let [first, last] = [rim[0].uv, rim.last().unwrap().uv];
    let axis = (0..2)
        .find(|&i| periods[i].is_some_and(|p| (last[i] - first[i]).abs() > p * 0.5))
        .ok_or_else(failure)?;
    let period = periods[axis].unwrap();
    let normal = 1 - axis;
    let positive = last[axis] > first[axis];
    let side = usize::from(positive == (same_sense == (axis == 0)));
    let bound = (0..source.chart_count()).filter_map(|i| source.domain(i)[normal][side].as_ref());
    let bound = if side == 0 { bound.min() } else { bound.max() };
    let mut uv = first;
    let (pole_v, on_oriented_side) = if let Some(bound) = bound {
        (bound.to_f64_nearest(), true)
    } else {
        // STEP cones have an unbounded axial chart. Intersect two sampled
        // meridians, then check that the other meridians meet there too.
        let radial = |h, dec: &mut Decider| -> Result<P, Refusal> {
            let mut a = first;
            a[normal] = h;
            let mut b = a;
            b[axis] += period * 0.5;
            let (a, b) = (carrier.position(a, dec)?, carrier.position(b, dec)?);
            Ok(std::array::from_fn(|i| b[i] - a[i]))
        };
        let h = first[normal];
        let step = h.abs().max(1.0) * if side == 0 { -1.0 } else { 1.0 };
        let a = radial(h, dec)?;
        let b = radial(h + step, dec)?;
        let i = (0..3)
            .max_by(|&i, &j| (b[i] - a[i]).abs().total_cmp(&(b[j] - a[j]).abs()))
            .unwrap();
        let t = a[i] / (a[i] - b[i]);
        if !t.is_finite() || t == 0.0 || (t < 0.0 && vertex_loop.is_none()) {
            return Err(failure());
        }
        (h + t * step, t > 0.0)
    };
    uv[normal] = pole_v;
    let pole = carrier.position(uv, dec)?;
    let rounding = body.positions[rim[0].mesh]
        .iter()
        .chain(&pole)
        .fold(1.0_f64, |scale, x| scale.max(x.abs()))
        * 128.0
        * f64::EPSILON;
    for station in [0.25, 0.5, 0.75] {
        let mut sample = uv;
        sample[axis] += station * period;
        if distance(carrier.position(sample, dec)?, pole) > rounding {
            return Err(failure());
        }
    }
    let chord = options.chord_tolerance.to_f64_nearest();
    let matching = vertex_loop.filter(|&(i, tolerance)| {
        distance(pole, body.positions[i]) <= tolerance.max(rounding).max(chord)
    });
    let mesh = match matching {
        Some((index, _)) => index,
        None if on_oriented_side => push(body, pole, options)?,
        None => return Err(failure()),
    };
    close_rim(rim, axis, pole_v, mesh, carrier, steps, body, options, dec)
}

/// Carrier points strictly between two parameter points of one seam line,
/// at `steps` equal steps, as new mesh points.
fn seam_line(
    carrier: &Carrier,
    [from, to]: [Uv; 2],
    steps: u32,
    body: &mut PreviewMesh<String, String>,
    options: &PreviewOptions,
    dec: &mut Decider,
) -> Result<Vec<Node>, Refusal> {
    (1..steps)
        .map(|step| {
            let t = f64::from(step) / f64::from(steps);
            let uv = std::array::from_fn(|i| from[i] + (to[i] - from[i]) * t);
            Ok(Node {
                uv,
                mesh: push(body, carrier.position(uv, dec)?, options)?,
            })
        })
        .collect()
}

/// The same seam points one turn along `axis`.
fn turned(nodes: &[Node], axis: usize, turn: f64) -> impl DoubleEndedIterator<Item = Node> + '_ {
    nodes.iter().map(move |node| {
        let mut node = *node;
        node.uv[axis] += turn;
        node
    })
}

/// Closes a winding rim at a pole parameter of the other axis.
///
/// The region runs up the seam from the rim's last station to the pole,
/// back under the rim along the pole line with one node per rim station, all
/// at the pole's one mesh vertex, and down the seam to the rim's first
/// station. Both seam sides are one row of carrier points, and the pole line
/// is split column by column. A side or pole line left as one chord makes
/// the grid triangulation fan whole rows of stations to its two ends, and on
/// a curved surface such a fan stands up as a fin.
#[allow(clippy::too_many_arguments)]
fn close_rim(
    mut rim: Vec<Node>,
    axis: usize,
    pole: f64,
    mesh: usize,
    carrier: &Carrier,
    steps: u32,
    body: &mut PreviewMesh<String, String>,
    options: &PreviewOptions,
    dec: &mut Decider,
) -> Result<Vec<Node>, Refusal> {
    let [first, last] = [rim[0].uv, rim.last().unwrap().uv];
    let mut top = first;
    top[1 - axis] = pole;
    let side = seam_line(carrier, [first, top], steps, body, options, dec)?;
    let turn = last[axis] - first[axis];
    let under: Vec<Node> = rim
        .iter()
        .rev()
        .map(|node| {
            let mut uv = node.uv;
            uv[1 - axis] = pole;
            Node { uv, mesh }
        })
        .collect();
    rim.extend(turned(&side, axis, turn));
    rim.extend(under);
    rim.extend(side.into_iter().rev());
    Ok(rim)
}

/// Whether a refusal is the caller's limit or cancellation rather than this
/// face's geometry.
fn control(reason: &Refusal) -> bool {
    matches!(
        reason.code,
        "numbers.deadline"
            | "numbers.operation-budget"
            | "numbers.budget-exhausted"
            | "mesh.preview-resource-limit"
    ) || catallax_numbers::current_operation_work().exhausted
}

/// A position, or `None` where the face's own geometry has none.
fn optional<T>(result: Result<T, Refusal>) -> Result<Option<T>, Refusal> {
    match result {
        Ok(value) => Ok(Some(value)),
        Err(reason) if control(&reason) => Err(reason),
        Err(_) => Ok(None),
    }
}

/// The finite limits of a native surface's charts along one axis.
fn limits(source: &NativeSurface, axis: usize) -> Option<[f64; 2]> {
    let mut range: Option<[f64; 2]> = None;
    for chart in 0..source.chart_count() {
        let [Some(low), Some(high)] = &source.domain(chart)[axis] else {
            return None;
        };
        let [low, high] = [low, high].map(catallax_numbers::Rational::to_f64_nearest);
        range = Some(range.map_or([low, high], |[a, b]| [a.min(low), b.max(high)]));
    }
    range.filter(|[low, high]| low < high)
}

/// The axes along which a native carrier with no declared period closes on
/// itself, as a spline whose first and last control rows coincide does. Its
/// first and last parameter lines must be one curve: nine samples along
/// them agree within `target`. Trims on such a carrier cross the seam
/// between those lines, and are read on one turn of the closed axis.
fn seams(
    source: &NativeSurface,
    carrier: &Carrier,
    target: f64,
    dec: &mut Decider,
) -> Result<Seams, Refusal> {
    let mut seams = [None; 2];
    for (axis, seam) in seams.iter_mut().enumerate() {
        if source.periods()[axis].is_some() {
            continue;
        }
        let (Some([low, high]), Some([from, to])) =
            (limits(source, axis), limits(source, 1 - axis))
        else {
            continue;
        };
        let mut closed = true;
        for station in 0..=8 {
            let mut start = [0.0; 2];
            start[1 - axis] = from + (to - from) * f64::from(station) / 8.0;
            start[axis] = low;
            let mut end = start;
            end[axis] = high;
            let (Some(start), Some(end)) = (
                optional(carrier.position(start, dec))?,
                optional(carrier.position(end, dec))?,
            ) else {
                closed = false;
                break;
            };
            let gap = distance(start, end);
            if gap.is_nan() || gap > target {
                closed = false;
                break;
            }
        }
        if closed {
            *seam = Some([low, high - low]);
        }
    }
    Ok(seams)
}

/// Closes a lone winding rim at the vertex its face declares with a
/// VERTEX_LOOP, where sampling alone cannot place the pole.
///
/// The pole is a parameter on the rim's normal axis: a finite chart limit,
/// or where the rim's sampled meridians meet, found by secant steps on the
/// rim's widest diameter. The carrier at that parameter, on four meridians a
/// quarter turn apart, must reach the declared point within its vertex
/// tolerance or the chord target. That point then closes the face.
#[allow(clippy::too_many_arguments)]
fn close_at_vertex(
    body: &mut PreviewMesh<String, String>,
    rim: Vec<Node>,
    (periods, steps): ([Option<f64>; 2], u32),
    source: &NativeSurface,
    carrier: &Carrier,
    (vertex, tolerance): (usize, f64),
    options: &PreviewOptions,
    dec: &mut Decider,
) -> Result<Vec<Node>, Refusal> {
    let unreached = || failure("single periodic rim does not reach its declared vertex");
    let [first, last] = [rim[0].uv, rim.last().unwrap().uv];
    let axis = (0..2)
        .find(|&i| periods[i].is_some_and(|p| (last[i] - first[i]).abs() > p * 0.5))
        .ok_or_else(unreached)?;
    let period = periods[axis].unwrap();
    let normal = 1 - axis;
    let range = limits(source, normal);
    let mut candidates: Vec<f64> = range.into_iter().flatten().collect();
    let diameter = |h: f64, dec: &mut Decider| -> Result<Option<P>, Refusal> {
        let mut a = first;
        a[normal] = h;
        let mut b = a;
        b[axis] += period * 0.5;
        let (Some(a), Some(b)) = (
            optional(carrier.position(a, dec))?,
            optional(carrier.position(b, dec))?,
        ) else {
            return Ok(None);
        };
        Ok(Some(std::array::from_fn(|i| b[i] - a[i])))
    };
    let start = first[normal];
    if let Some(across) = diameter(start, dec)? {
        let i = (0..3)
            .max_by(|&i, &j| across[i].abs().total_cmp(&across[j].abs()))
            .unwrap();
        let mut a = (start, across[i]);
        let step = start.abs().max(1.0);
        if let Some(next) = diameter(start + step, dec)? {
            let mut b = (start + step, next[i]);
            for _ in 0..64 {
                dec.charge_work(1)?;
                let next = b.0 - b.1 * (b.0 - a.0) / (b.1 - a.1);
                if !next.is_finite() || next.to_bits() == b.0.to_bits() {
                    break;
                }
                let Some(value) = diameter(next, dec)? else {
                    break;
                };
                a = b;
                b = (next, value[i]);
            }
            if b.0.is_finite() {
                candidates.push(range.map_or(b.0, |[low, high]| b.0.clamp(low, high)));
            }
        }
    }
    let target = body.positions[vertex];
    let reach = tolerance.max(options.chord_tolerance.to_f64_nearest());
    let mut best: Option<(f64, f64)> = None;
    for pole in candidates {
        let mut farthest = 0.0f64;
        for station in 0..4 {
            let mut uv = first;
            uv[normal] = pole;
            uv[axis] += period * f64::from(station) * 0.25;
            let Some(point) = optional(carrier.position(uv, dec))? else {
                farthest = f64::INFINITY;
                break;
            };
            farthest = farthest.max(distance(point, target));
        }
        if farthest <= reach && best.is_none_or(|(_, old)| farthest < old) {
            best = Some((pole, farthest));
        }
    }
    let (pole, _) = best.ok_or_else(unreached)?;
    close_rim(rim, axis, pole, vertex, carrier, steps, body, options, dec)
}

/// Whether every edge loop of a face runs out along its edges and straight
/// back over the same parameter trace, as a face does whose two edges
/// admission welded into one. Such a face has no area at any sampling. Two
/// uses of one seam edge, a period apart, do not cancel.
fn encloses_nothing(face: &TolerantFace, trims: &FaceTrims<'_>) -> bool {
    let same_trace = |a: &NativeTrim, b: &NativeTrim| {
        a.samples.len() >= 2
            && b.samples.len() >= 2
            && (0..=8).all(|station| {
                let fraction = f64::from(station) / 8.0;
                let [p, q] = [trim_seed(a, fraction), trim_seed(b, 1.0 - fraction)];
                (0..2).all(|i| {
                    (p.parameters[i] - q.parameters[i]).abs()
                        <= 1e-9 * (1.0 + p.parameters[i].abs())
                })
            })
    };
    let mut edge_loops = 0;
    for (boundary_index, boundary) in face.loops.iter().enumerate() {
        if boundary.edges.is_empty() {
            continue;
        }
        edge_loops += 1;
        let Some(rows) = trims.get(&(face.name.as_str(), boundary_index)) else {
            return false;
        };
        if rows.len() != boundary.edges.len() {
            return false;
        }
        let cancels = |a: (&TolerantUse, &NativeTrim), b: (&TolerantUse, &NativeTrim)| {
            a.0.edge == b.0.edge && a.0.forward != b.0.forward && same_trace(a.1, b.1)
        };
        let mut open = Vec::new();
        for use_and_trim in boundary.edges.iter().zip(rows.iter().copied()) {
            if open.last().is_some_and(|&last| cancels(last, use_and_trim)) {
                open.pop();
            } else {
                open.push(use_and_trim);
            }
        }
        while open.len() >= 2 && cancels(open[0], *open.last().unwrap()) {
            open.pop();
            open.remove(0);
        }
        if !open.is_empty() {
            return false;
        }
    }
    edge_loops > 0
}

/// Source trims of each face boundary, in use order.
type FaceTrims<'a> = BTreeMap<(&'a str, usize), Vec<&'a NativeTrim>>;

/// Each drawn edge's samples: fraction along the edge and mesh point.
type EdgeRows<'a> = BTreeMap<&'a str, Vec<(f64, usize)>>;

/// A face left out of its body's display because it has no area. A face
/// drawn from its carrier's natural domain is never one.
fn sliver(import: &TolerantImport, body: &str, face: &TolerantFace, trims: &FaceTrims<'_>) -> bool {
    !import
        .evidence()
        .native_domains
        .iter()
        .any(|domain| domain.body == body && domain.face == face.name)
        && encloses_nothing(face, trims)
}

/// The trim chords, as the mesh indices at their ends, that cross or touch
/// another trim chord in the parameter plane.
fn crossing_chords(loops: &[Vec<Node>], dec: &mut Decider) -> Result<Vec<[usize; 2]>, Refusal> {
    let mut chords = Vec::new();
    for ring in loops {
        for i in 0..ring.len() {
            let [a, b] = [ring[i], ring[(i + 1) % ring.len()]];
            if a.mesh != b.mesh || a.uv.map(f64::to_bits) != b.uv.map(f64::to_bits) {
                chords.push([a, b]);
            }
        }
    }
    let span = |c: &[Node; 2]| [c[0].uv[0].min(c[1].uv[0]), c[0].uv[0].max(c[1].uv[0])];
    chords.sort_by(|a, b| span(a)[0].total_cmp(&span(b)[0]));
    let orient =
        |a: Uv, b: Uv, c: Uv| (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
    let within = |a: Uv, b: Uv, c: Uv| {
        (a[0].min(b[0])..=a[0].max(b[0])).contains(&c[0])
            && (a[1].min(b[1])..=a[1].max(b[1])).contains(&c[1])
    };
    let mut crossing = std::collections::BTreeSet::new();
    let mut tests = 0u64;
    for (i, p) in chords.iter().enumerate() {
        let reach = span(p)[1];
        for q in chords[i + 1..].iter().take_while(|q| span(q)[0] <= reach) {
            tests += 1;
            let shared = |n: &Node| {
                q.iter()
                    .any(|m| m.mesh == n.mesh && m.uv.map(f64::to_bits) == n.uv.map(f64::to_bits))
            };
            if p.iter().any(shared) {
                continue;
            }
            let [p0, p1, q0, q1] = [p[0].uv, p[1].uv, q[0].uv, q[1].uv];
            let d = [
                orient(q0, q1, p0),
                orient(q0, q1, p1),
                orient(p0, p1, q0),
                orient(p0, p1, q1),
            ];
            let proper = d[0] * d[1] < 0.0 && d[2] * d[3] < 0.0;
            let touching = (d[0] == 0.0 && within(q0, q1, p0))
                || (d[1] == 0.0 && within(q0, q1, p1))
                || (d[2] == 0.0 && within(p0, p1, q0))
                || (d[3] == 0.0 && within(p0, p1, q1));
            if proper || touching {
                for chord in [p, q] {
                    let [a, b] = [chord[0].mesh, chord[1].mesh];
                    crossing.insert([a.min(b), a.max(b)]);
                }
            }
        }
    }
    dec.charge_traversal(tests + chords.len() as u64)?;
    Ok(crossing.into_iter().collect())
}

/// Asks for the source curve's own point halfway along every chord of a
/// curved edge that crosses another trim chord, by edge and fraction. A
/// shared edge gains the same samples in every face that uses it.
fn request_splits(
    loops: &[Vec<Node>],
    rows: &EdgeRows<'_>,
    curved: &std::collections::BTreeSet<&str>,
    requests: &mut BTreeMap<String, Vec<f64>>,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let mut chords = BTreeMap::new();
    for (name, row) in rows.iter().filter(|(name, _)| curved.contains(*name)) {
        for pair in row.windows(2) {
            let [a, b] = [pair[0].1, pair[1].1];
            chords
                .entry([a.min(b), a.max(b)])
                .or_insert((*name, (pair[0].0 + pair[1].0) * 0.5));
        }
    }
    for chord in crossing_chords(loops, dec)? {
        if let Some((name, middle)) = chords.get(&chord) {
            requests
                .entry((*name).to_owned())
                .or_default()
                .push(*middle);
        }
    }
    Ok(())
}

// Cut two noncontractible rims along one common source-parameter seam. The
// additional cut nodes are display samples, not new pickable B-rep edges.
fn cut_band(
    mut rims: Vec<Vec<Node>>,
    carrier: &Carrier,
    mesh: &mut PreviewMesh<String, String>,
    options: &PreviewOptions,
    dec: &mut Decider,
    grid: Option<u32>,
) -> Result<Vec<Node>, Refusal> {
    if rims.len() != 2 {
        return Err(failure(format!(
            "periodic trim has {} winding rims; two must be paired",
            rims.len()
        )));
    }
    let delta = |rim: &[Node]| {
        [
            rim.last().unwrap().uv[0] - rim[0].uv[0],
            rim.last().unwrap().uv[1] - rim[0].uv[1],
        ]
    };
    let first = delta(&rims[0]);
    let axis = usize::from(first[1].abs() > first[0].abs());
    let period = first[axis].abs();
    if !period.is_finite() || period == 0.0 {
        return Err(failure("periodic rim has no finite winding"));
    }
    let cut = rims[0][0].uv[axis];
    for rim in &mut rims {
        let winding = delta(rim);
        if (winding[axis].abs() - period).abs() > period * 1e-6
            || winding[1 - axis].abs() > period * 1e-6
        {
            return Err(failure("paired rims have incompatible native winding"));
        }
        if winding[axis] < 0.0 {
            rim.reverse();
        }
        if grid.is_some() {
            rim.pop();
            let distance = |node: &Node| {
                let delta = (node.uv[axis] - cut) / period;
                (delta - delta.round()).abs()
            };
            let split = rim
                .iter()
                .enumerate()
                .min_by(|(_, a), (_, b)| distance(a).total_cmp(&distance(b)))
                .unwrap()
                .0;
            let shift = ((cut - rim[split].uv[axis]) / period).round() * period;
            for (i, node) in rim.iter_mut().enumerate() {
                node.uv[axis] += shift + if i < split { period } else { 0.0 };
            }
            rim.rotate_left(split);
            let mut last = rim[0];
            last.uv[axis] += period;
            rim.push(last);
            continue;
        }
        let shift = ((cut - rim[0].uv[axis]) / period).floor() * period;
        for node in rim.iter_mut() {
            node.uv[axis] += shift;
        }
        let target = if cut < rim[0].uv[axis] {
            cut + period
        } else {
            cut
        };
        let split = rim
            .windows(2)
            .position(|pair| pair[0].uv[axis] <= target && target <= pair[1].uv[axis])
            .ok_or_else(|| failure("native periodic rim is not monotone at the chosen seam"))?;
        let a = rim[split];
        let b = rim[split + 1];
        let fraction = if (b.uv[axis] - a.uv[axis]).abs() == 0.0 {
            0.0
        } else {
            (target - a.uv[axis]) / (b.uv[axis] - a.uv[axis])
        };
        let uv = std::array::from_fn(|i| a.uv[i] + fraction * (b.uv[i] - a.uv[i]));
        let cut_node = if fraction == 0.0 {
            a
        } else if fraction.to_bits() == 1.0_f64.to_bits() {
            b
        } else {
            Node {
                uv,
                mesh: push(mesh, carrier.position(uv, dec)?, options)?,
            }
        };
        let mut rotated = Vec::with_capacity(rim.len() + 1);
        rotated.push(cut_node);
        rotated.extend(rim[split + 1..rim.len() - 1].iter().copied());
        for node in &rim[..=split] {
            let mut node = *node;
            node.uv[axis] += period;
            rotated.push(node);
        }
        let mut end = cut_node;
        end.uv[axis] += period;
        rotated.push(end);
        rotated.dedup_by(|a, b| a.uv.iter().zip(b.uv).all(|(a, b)| (*a - b).abs() == 0.0));
        let shift = cut - rotated[0].uv[axis];
        if (shift / period - (shift / period).round()).abs() > 1e-6 {
            return Err(failure(
                "periodic seam alignment changed source parameterization",
            ));
        }
        for node in &mut rotated {
            node.uv[axis] += shift;
        }
        *rim = rotated;
    }
    let mut first = rims.remove(0);
    let second = rims.remove(0);
    let Some(steps) = grid else {
        first.extend(second.into_iter().rev());
        return Ok(first);
    };
    // The seam between the rims is one row of carrier points on both sides
    // of the region, as at a pole (see `close_rim`).
    let side = seam_line(
        carrier,
        [second[0].uv, first[0].uv],
        steps,
        mesh,
        options,
        dec,
    )?;
    first.extend(turned(&side, axis, period).rev());
    first.extend(second.into_iter().rev());
    first.extend(side);
    Ok(first)
}

struct Omission {
    kind: &'static str,
    name: String,
    reason: Refusal,
}

fn omit(
    omissions: &mut Option<&mut Vec<Omission>>,
    kind: &'static str,
    name: &str,
    reason: Refusal,
) -> Result<(), Refusal> {
    let Some(omissions) = omissions.as_mut() else {
        return Err(reason);
    };
    if control(&reason) {
        return Err(reason);
    }
    omissions.push(Omission {
        kind,
        name: name.into(),
        reason,
    });
    Ok(())
}

/// Preview one admitted native definition in its source-body coordinates.
/// Assembly placement remains in the source document, not in this body-local
/// preview. The returned triangles have no watertightness or deviation claim.
/// If a retained trim cannot be drawn, returns the usable faces and shared
/// edges and logs the omissions. Draw the edges even when `triangles` is
/// empty: some admitted bodies have a visible boundary but no drawable face.
/// This preview never certifies a body.
///
/// A face whose every boundary runs out along its edges and back over the
/// same trace has no area, so it is left out of the face set, and the log
/// says so.
pub fn preview_tolerant_body(
    import: &TolerantImport,
    index: usize,
    options: &PreviewOptions,
    dec: &mut Decider,
) -> Result<PreviewMesh<String, String>, Refusal> {
    super::validate(options)?;
    if index >= import.bodies().len() {
        return Err(failure("body index is out of range"));
    }
    let preview = preview_body(import, index, options, dec);
    if preview.is_ok() {
        let body = import.bodies()[index].definition();
        let (_, trims) = body_trims(import, body.name());
        let mut log = std::io::stderr().lock();
        for face in body.shells().iter().flat_map(|shell| &shell.faces) {
            if sliver(import, body.name(), face, &trims) {
                let _ = writeln!(
                    log,
                    "native preview dropped face {}: its boundary runs out and back along one trace, so it has no area to draw",
                    face.name,
                );
            }
        }
    }
    preview
}

fn preview_body(
    import: &TolerantImport,
    index: usize,
    options: &PreviewOptions,
    dec: &mut Decider,
) -> Result<PreviewMesh<String, String>, Refusal> {
    if let Ok(mesh) = refund_on_refusal(|| build(import, index, options, None, None, dec)) {
        return Ok(mesh);
    }
    let complete = refund_on_refusal(|| {
        super::coarsen(options, dec, |divisions, _, dec| {
            build(import, index, options, Some(divisions), None, dec)
        })
    });
    if complete.is_ok() {
        return complete;
    }
    let mut omissions = Vec::new();
    let mesh = super::coarsen(options, dec, |divisions, _, dec| {
        omissions.clear();
        build(
            import,
            index,
            options,
            Some(divisions),
            Some(&mut omissions),
            dec,
        )
    })?;
    if mesh.triangles.is_empty()
        && !mesh.edges.iter().any(|edge| {
            edge.vertices
                .windows(2)
                .any(|pair| distance(mesh.positions[pair[0]], mesh.positions[pair[1]]) > 0.0)
        })
    {
        return complete;
    }
    let mut log = std::io::stderr().lock();
    for omitted in omissions {
        let _ = writeln!(
            log,
            "native preview omitted {} {}: {}: {}",
            omitted.kind, omitted.name, omitted.reason.code, omitted.reason.detail,
        );
    }
    Ok(mesh)
}

/// A body's source trims: the first for each edge, and each face boundary's
/// in use order.
fn body_trims<'a>(
    import: &'a TolerantImport,
    body: &str,
) -> (BTreeMap<&'a str, &'a NativeTrim>, FaceTrims<'a>) {
    let mut first_trim = BTreeMap::new();
    let mut face_trims = FaceTrims::new();
    for trim in import
        .evidence()
        .native_trims
        .iter()
        .filter(|trim| trim.body == body)
    {
        first_trim.entry(trim.edge.as_str()).or_insert(trim);
        face_trims
            .entry((&trim.face, trim.boundary))
            .or_default()
            .push(trim);
    }
    for trims in face_trims.values_mut() {
        trims.sort_by_key(|trim| trim.use_index);
    }
    (first_trim, face_trims)
}

/// Rounds of chord splitting a last pass may spend on trims that cross.
const SPLIT_ROUNDS: usize = 8;

/// Builds the body's display. The last pass, which omits what it cannot
/// draw, first splits trim chords that cross: each round adds the source
/// curve's own midpoint to every crossing chord of a shared edge and builds
/// again, until nothing crosses or no chord can be split.
fn build(
    import: &TolerantImport,
    index: usize,
    options: &PreviewOptions,
    grid: Option<u32>,
    mut omissions: Option<&mut Vec<Omission>>,
    dec: &mut Decider,
) -> Result<PreviewMesh<String, String>, Refusal> {
    let mut refined = BTreeMap::<String, Vec<f64>>::new();
    for round in 0.. {
        if let Some(omissions) = omissions.as_deref_mut() {
            omissions.clear();
        }
        let mut requests = BTreeMap::new();
        let splitting = omissions.is_some() && round < SPLIT_ROUNDS;
        let mesh = build_once(
            import,
            index,
            options,
            grid,
            omissions.as_deref_mut(),
            (&refined, splitting.then_some(&mut requests)),
            dec,
        )?;
        let mut added = false;
        for (edge, fractions) in requests {
            let row = refined.entry(edge).or_default();
            for fraction in fractions {
                let at = row.partition_point(|f| *f < fraction);
                if row
                    .get(at)
                    .is_none_or(|f| f.to_bits() != fraction.to_bits())
                {
                    row.insert(at, fraction);
                    added = true;
                }
            }
        }
        if !added {
            return Ok(mesh);
        }
    }
    unreachable!("the rounds above end")
}

/// Extra edge samples, by edge, and where a pass may ask for more.
type Splits<'a> = (
    &'a BTreeMap<String, Vec<f64>>,
    Option<&'a mut BTreeMap<String, Vec<f64>>>,
);

fn build_once(
    import: &TolerantImport,
    index: usize,
    options: &PreviewOptions,
    grid: Option<u32>,
    mut omissions: Option<&mut Vec<Omission>>,
    (refined, mut requests): Splits<'_>,
    dec: &mut Decider,
) -> Result<PreviewMesh<String, String>, Refusal> {
    let body = import.bodies()[index].definition();
    let (first_trim, face_trims) = body_trims(import, body.name());
    let mut mesh = PreviewMesh::default();
    let vertices = body
        .vertices()
        .iter()
        .map(|vertex| push(&mut mesh, point(&vertex.source.point)?, options))
        .collect::<Result<Vec<_>, _>>()?;
    let mut edges = BTreeMap::new();
    for edge in body.edges() {
        let Some(trim) = first_trim.get(edge.name.as_str()) else {
            continue;
        };
        let position_start = mesh.positions.len();
        let sample = |mesh: &mut PreviewMesh<String, String>, dec: &mut Decider| {
            let ends = edge.vertices.map(|vertex| vertices[vertex]);
            let carrier = if edge.degenerate || matches!(edge.curve, Curve::Line(..)) {
                EdgeCarrier::Line(ends.map(|index| mesh.positions[index]))
            } else {
                let source_id = edge
                    .name
                    .strip_prefix('#')
                    .and_then(|name| name.parse().ok());
                let curve = match source_id
                    .map(|entity| import.parsed_source().topology(entity, dec))
                    .transpose()
                    .map_err(|error| failure(format!("{error:?}")))?
                {
                    Some(Topology::Edge {
                        curve: Some(curve), ..
                    }) => import.native_curve(curve),
                    _ => None,
                }
                .ok_or_else(|| failure("active edge has no retained native carrier"))?
                .for_display(dec)
                .map(Box::new)
                .map_err(source_error)?;
                let first = trim.samples.first().unwrap().curve.parameters[0];
                let last = trim.samples.last().unwrap().curve.parameters[0];
                EdgeCarrier::Native(
                    curve,
                    if trim.forward {
                        [first, last]
                    } else {
                        [last, first]
                    },
                )
            };
            let mut samples = if let Some(divisions) = grid {
                grid_edge(&carrier, ends, divisions, mesh, options, dec)
            } else {
                edge_samples(
                    &carrier,
                    ends,
                    options
                        .chord_tolerance
                        .to_f64_nearest()
                        .max(edge.tolerance_mm().to_f64_nearest()),
                    mesh,
                    options,
                    dec,
                )
            }?;
            // Points of the source curve itself where trim chords crossed.
            for &fraction in refined.get(edge.name.as_str()).into_iter().flatten() {
                let at = samples.partition_point(|(f, _)| *f < fraction);
                if samples
                    .get(at)
                    .is_none_or(|(f, _)| f.to_bits() != fraction.to_bits())
                {
                    let index = push(mesh, carrier.sample(fraction, dec)?, options)?;
                    samples.insert(at, (fraction, index));
                }
            }
            Ok(samples)
        };
        let samples = match sample(&mut mesh, dec) {
            Ok(samples) => samples,
            Err(reason) => {
                omit(&mut omissions, "edge", &edge.name, reason)?;
                mesh.positions.truncate(position_start);
                continue;
            }
        };
        mesh.edges.push(PreviewEdge {
            edge: edge.name.clone(),
            vertices: samples.iter().map(|(_, index)| *index).collect(),
        });
        edges.insert(edge.name.as_str(), samples);
    }
    // Only a curved edge's own points can move its chords; a straight
    // edge's chord is the edge. Only the last pass asks.
    let curved: std::collections::BTreeSet<&str> = body
        .edges()
        .iter()
        .filter(|edge| requests.is_some() && !edge.degenerate)
        .filter(|edge| !matches!(edge.curve, Curve::Line(..)))
        .map(|edge| edge.name.as_str())
        .collect();
    for face in body.shells().iter().flat_map(|shell| &shell.faces) {
        if sliver(import, body.name(), face, &face_trims) {
            continue;
        }
        let start = (mesh.positions.len(), mesh.triangles.len(), mesh.faces.len());
        let mut draw = |mesh: &mut PreviewMesh<String, String>, dec: &mut Decider| {
            let collapsed_vertices = face
                .loops
                .iter()
                .flat_map(|boundary| &boundary.edges)
                .map(|usage| &body.edges()[usage.edge])
                .filter(|edge| edge.degenerate)
                .flat_map(|edge| edge.vertices.map(|vertex| vertices[vertex]))
                .collect::<std::collections::BTreeSet<_>>();
            let source_id = face
                .name
                .strip_prefix('#')
                .and_then(|name| name.parse().ok());
            let surface = match source_id
                .map(|entity| import.parsed_source().topology(entity, dec))
                .transpose()
                .map_err(|error| failure(format!("{error:?}")))?
            {
                Some(Topology::Face {
                    surface: Some(surface),
                    ..
                }) => import.native_surface(surface),
                _ => None,
            };
            let mut carrier = match surface {
                Some(surface) => Carrier::Native(
                    Box::new(surface.for_display(dec).map_err(source_error)?),
                    [None; 2],
                ),
                None => Carrier::Inferred(SourceChart::new(&face.surface, dec)?),
            };
            let mut face_options = options.clone();
            if let Some(tolerance) = import.evidence().tolerances.iter().find(|entry| {
                entry.body == body.name() && entry.kind == "face" && entry.entity == face.name
            }) {
                face_options.chord_tolerance = face_options.chord_tolerance.max(
                    Rational::from_f64_exact(tolerance.tolerance_mm)
                        .ok_or_else(|| failure("nonfinite face tolerance"))?,
                );
            }
            let chord = face_options.chord_tolerance.to_f64_nearest();
            if let Some(source) = surface {
                let closed = seams(source, &carrier, chord * 0.125, dec)?;
                if let Carrier::Native(_, seams) = &mut carrier {
                    *seams = closed;
                }
            }
            // Declared periods, and one turn of each closed spline axis.
            let periods: [Option<f64>; 2] = match &carrier {
                Carrier::Native(surface, seams) => std::array::from_fn(|axis| {
                    surface.periods()[axis].or(seams[axis].map(|[_, width]| width))
                }),
                Carrier::Inferred(_) => [None; 2],
            };
            let mut loops = Vec::new();
            let mut winding = Vec::new();
            for (boundary_index, boundary) in face.loops.iter().enumerate() {
                if boundary.edges.is_empty() {
                    continue;
                }
                let trim_rows = face_trims
                    .get(&(face.name.as_str(), boundary_index))
                    .ok_or_else(|| failure("native trim evidence is missing"))?;
                let mut nodes: Vec<Node> = Vec::new();
                for offset in 0..boundary.edges.len() {
                    let usage_index = if boundary.orientation {
                        offset
                    } else {
                        boundary.edges.len() - 1 - offset
                    };
                    let usage = &boundary.edges[usage_index];
                    let trim = trim_rows[usage_index];
                    let edge = &body.edges()[usage.edge];
                    let row = edges.get(edge.name.as_str()).ok_or_else(|| {
                        failure(format!("shared edge {} has no preview samples", edge.name))
                    })?;
                    // What the body records for this edge and the vertices it
                    // was made to end at, on top of the face's own chord target.
                    let recorded = edge
                        .vertices
                        .iter()
                        .map(|vertex| body.vertices()[*vertex].tolerance_mm())
                        .fold(edge.tolerance_mm(), |most, each| most.max(each))
                        .to_f64_nearest();
                    let reach = face_options.chord_tolerance.to_f64_nearest() + recorded;
                    for offset in 0..row.len() {
                        let (fraction, index) = row[if usage.forward == boundary.orientation {
                            offset
                        } else {
                            row.len() - 1 - offset
                        }];
                        let collapsed = collapsed_vertices.contains(&index);
                        if !collapsed && nodes.last().is_some_and(|node| node.mesh == index) {
                            continue;
                        }
                        let seed = trim_seed(
                            trim,
                            if usage.forward {
                                fraction
                            } else {
                                1.0 - fraction
                            },
                        );
                        // A coarse pass reads parameters off the trim evidence,
                        // except on an edge split where trims crossed. Read off
                        // the evidence, a split sample would land on the chord
                        // between two evidence samples, so it is projected onto
                        // the carrier instead.
                        let mut uv = match &carrier {
                            Carrier::Native(surface, _)
                                if grid.is_none() || refined.contains_key(edge.name.as_str()) =>
                            {
                                let projection = surface
                                    .project(mesh.positions[index], Some(seed), dec)
                                    .map_err(source_error)?;
                                // The face's triangles are measured on its carrier,
                                // so a boundary sample's own miss is checked here,
                                // once, against what the body admits for it.
                                if projection.residual_mm > reach {
                                    return Err(failure(format!(
                                        "a sample of edge {} is {:e} mm from face {}, past the {recorded:e} mm recorded for the edge",
                                        edge.name, projection.residual_mm, face.name
                                    )));
                                }
                                projection.location.parameters
                            }
                            _ => seed.parameters,
                        };
                        // Across a seam, keep the trim on the turn of its
                        // previous sample.
                        if let Carrier::Native(_, seams) = &carrier
                            && let Some(previous) = nodes.last()
                        {
                            for (axis, seam) in seams.iter().enumerate() {
                                if let Some([_, width]) = *seam {
                                    uv[axis] +=
                                        ((previous.uv[axis] - uv[axis]) / width).round() * width;
                                }
                            }
                        }
                        if nodes.last().is_some_and(|node| {
                            node.mesh == index
                                && (node.uv[0] - uv[0]).hypot(node.uv[1] - uv[1]) <= 1e-8
                        }) {
                            continue;
                        }
                        nodes.push(Node { uv, mesh: index });
                    }
                }
                if nodes.len() < 4 || nodes[0].mesh != nodes.last().unwrap().mesh {
                    return Err(failure("native trim loop is not closed"));
                }
                let first = nodes[0];
                let last = *nodes.last().unwrap();
                let first_usage = &boundary.edges[if boundary.orientation {
                    0
                } else {
                    boundary.edges.len() - 1
                }];
                let vertex = body.edges()[first_usage.edge].vertices
                    [usize::from(first_usage.forward != boundary.orientation)];
                let tolerance = body.vertices()[vertex].tolerance_mm().to_f64_nearest();
                for endpoint in [first, last].into_iter().filter(|_| grid.is_none()) {
                    if distance(
                        carrier.position(endpoint.uv, dec)?,
                        mesh.positions[first.mesh],
                    ) > tolerance
                    {
                        return Err(failure(
                            "native UV closure exceeds its measured vertex tolerance",
                        ));
                    }
                }
                let periodic = match &carrier {
                    Carrier::Native(..) => periods.into_iter().enumerate().any(|(axis, period)| {
                        period.is_some_and(|period| {
                            let delta = last.uv[axis] - first.uv[axis];
                            let turns = (delta / period).round();
                            turns.abs() >= 1.0 && (delta - turns * period).abs() <= period * 1e-6
                        })
                    }),
                    Carrier::Inferred(surface) if grid.is_some() => surface
                        .periodic()
                        .into_iter()
                        .enumerate()
                        .any(|(axis, periodic)| {
                            periodic && (last.uv[axis] - first.uv[axis]).abs() > 0.5
                        }),
                    Carrier::Inferred(_) => {
                        (first.uv[0] - last.uv[0]).hypot(first.uv[1] - last.uv[1]) > 1e-8
                    }
                };
                if periodic {
                    winding.push(nodes);
                } else {
                    nodes.pop();
                    loops.push(nodes);
                }
            }
            let natural = import
                .evidence()
                .native_domains
                .iter()
                .find(|domain| domain.body == body.name() && domain.face == face.name);
            if let Some(domain) = natural {
                loops.push(domain_loop(
                    &carrier,
                    domain.bounds,
                    mesh,
                    &face_options,
                    dec,
                )?);
            }
            if grid.is_some()
                && winding.len() == 1
                && let Some(source) = surface
            {
                let mut poles = face.loops.iter().filter_map(|boundary| boundary.vertex);
                let explicit = poles.next().filter(|_| poles.next().is_none()).map(|i| {
                    (
                        vertices[i],
                        body.vertices()[i].tolerance_mm().to_f64_nearest(),
                    )
                });
                let rim = winding.pop().unwrap();
                let steps = grid.unwrap_or(1);
                let closed = cut_to_pole(
                    mesh,
                    rim.clone(),
                    source,
                    (periods, steps),
                    &carrier,
                    (face.same_sense, explicit),
                    &face_options,
                    dec,
                );
                loops.push(match (closed, explicit) {
                    (Ok(ring), _) => ring,
                    (Err(reason), Some(vertex)) if !control(&reason) => close_at_vertex(
                        mesh,
                        rim,
                        (periods, steps),
                        source,
                        &carrier,
                        vertex,
                        &face_options,
                        dec,
                    )?,
                    (Err(reason), _) => return Err(reason),
                });
            } else if !winding.is_empty() {
                loops.push(cut_band(winding, &carrier, mesh, &face_options, dec, grid)?);
            }
            let first = loops
                .first()
                .and_then(|nodes| nodes.first())
                .ok_or_else(|| failure("face has no sampled parameter region"))?;
            let mut bounds = [first.uv; 2];
            for node in loops.iter().flatten() {
                for axis in 0..2 {
                    bounds[0][axis] = bounds[0][axis].min(node.uv[axis]);
                    bounds[1][axis] = bounds[1][axis].max(node.uv[axis]);
                }
            }
            if let Some(divisions) = grid {
                let trims = requests.is_some().then(|| loops.clone());
                let drawn = grid_face(
                    &face.name,
                    face.same_sense,
                    &carrier,
                    bounds,
                    &mut loops,
                    mesh,
                    divisions,
                    &face_options,
                    dec,
                );
                if let (Err(reason), Some(requests), Some(trims)) =
                    (&drawn, requests.as_deref_mut(), trims)
                    && reason.code == GRID_TRIM_REFUSAL
                {
                    request_splits(&trims, &edges, &curved, requests, dec)?;
                }
                return drawn;
            }
            let mut cdt = Cdt::new();
            for boundary in &loops {
                for i in 0..boundary.len() {
                    dec.charge_work(1)?;
                    let a = insert(&mut cdt, boundary[i])?;
                    let b = insert(&mut cdt, boundary[(i + 1) % boundary.len()])?;
                    if a == b || !cdt.can_add_constraint(a, b) {
                        return Err(failure("native preview trim constraints cross or collapse"));
                    }
                    cdt.add_constraint(a, b);
                }
            }
            seed_rectangles(&carrier, bounds, &loops, &mut cdt, mesh, &face_options, dec)?;
            finish_face(
                FaceOutput {
                    name: &face.name,
                    same_sense: face.same_sense,
                    skip_collapsed: natural.is_some()
                        || face
                            .loops
                            .iter()
                            .flat_map(|boundary| &boundary.edges)
                            .any(|usage| body.edges()[usage.edge].degenerate),
                    // A vertex of a tolerant body is where its edges were made to
                    // meet, up to its tolerance away from each face's carrier.
                    // The projection above has bounded that miss. An inferred
                    // carrier takes its parameters from the trim evidence without
                    // a projection, and is measured as it was.
                    off_carrier_nodes: matches!(carrier, Carrier::Native(..)),
                },
                &carrier,
                cdt,
                &loops,
                mesh,
                &face_options,
                dec,
            )
        };
        if let Err(reason) = draw(&mut mesh, dec) {
            omit(&mut omissions, "face", &face.name, reason)?;
            mesh.positions.truncate(start.0);
            mesh.triangles.truncate(start.1);
            mesh.faces.truncate(start.2);
        }
    }
    dec.check_deadline()?;
    Ok(mesh)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two straight degree-one splines between the same two vertices, 1.5 um
    /// apart, which admission welds into one edge.
    const WELDED_DIGON: &str = "ISO-10303-21;
HEADER;
FILE_DESCRIPTION(('welded digon'),'2;1');
FILE_NAME('digon.step','',(''),(''),'','','');
FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));
ENDSEC;
DATA;
#1=CARTESIAN_POINT('',(1.,2.,0.));
#2=DIRECTION('',(0.,0.,1.));
#3=DIRECTION('',(1.,0.,0.));
#4=AXIS2_PLACEMENT_3D('',#1,#2,#3);
#5=PLANE('',#4);
#6=CARTESIAN_POINT('',(1.,2.,0.));
#7=CARTESIAN_POINT('',(1.2,2.01,0.));
#8=VERTEX_POINT('',#6);
#9=VERTEX_POINT('',#7);
#10=B_SPLINE_CURVE_WITH_KNOTS('',1,(#6,#7),.UNSPECIFIED.,.F.,.F.,(2,2),(0.,1.),.UNSPECIFIED.);
#11=CARTESIAN_POINT('',(1.2,2.0115,0.));
#12=CARTESIAN_POINT('',(1.,2.0015,0.));
#13=B_SPLINE_CURVE_WITH_KNOTS('',1,(#11,#12),.UNSPECIFIED.,.F.,.F.,(2,2),(0.,1.),.UNSPECIFIED.);
#14=EDGE_CURVE('',#8,#9,#10,.T.);
#15=EDGE_CURVE('',#9,#8,#13,.T.);
#16=ORIENTED_EDGE('',*,*,#14,.T.);
#17=ORIENTED_EDGE('',*,*,#15,.T.);
#18=EDGE_LOOP('',(#16,#17));
#19=FACE_BOUND('',#18,.T.);
#20=FACE_SURFACE('',(#19),#5,.T.);
#21=CLOSED_SHELL('',(#20));
#22=SHELL_BASED_SURFACE_MODEL('',(#21));
#23=( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.) );
#24=( GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNIT_ASSIGNED_CONTEXT((#23)) REPRESENTATION_CONTEXT('','') );
#25=MANIFOLD_SURFACE_SHAPE_REPRESENTATION('',(#22),#24);
ENDSEC;
END-ISO-10303-21;";

    #[test]
    fn a_welded_digon_is_left_out_rather_than_refused() {
        let source = WELDED_DIGON.as_bytes();
        let import = catallax_import::tolerant::admission::import_step(
            source,
            &catallax_import::sha256_hex(source),
            Default::default(),
            &mut Decider::default(),
        )
        .unwrap();
        let body = import.bodies()[0].definition();
        let (_, trims) = body_trims(&import, body.name());
        let [face] = body.shells()[0].faces.as_slice() else {
            panic!("one face");
        };
        let [boundary] = face.loops.as_slice() else {
            panic!("one boundary");
        };
        // Admission made the two source edges one, used out and back.
        assert_eq!(boundary.edges.len(), 2);
        assert_eq!(boundary.edges[0].edge, boundary.edges[1].edge);
        assert!(sliver(&import, body.name(), face, &trims));
        // The detailed first pass draws the body: its edge and no face.
        let mesh = build(
            &import,
            0,
            &PreviewOptions::default(),
            None,
            None,
            &mut Decider::default(),
        )
        .unwrap();
        assert!(mesh.faces.is_empty() && mesh.triangles.is_empty());
        assert_eq!(mesh.edges.len(), 1);
    }
}
