use std::collections::VecDeque;
use std::ops::Range;

use super::{DisplayFace, DisplayLoop};
use spade::{
    ConstrainedDelaunayTriangulation, HasPosition, Point2, PositionInTriangulation, Triangulation,
};

use TrimRefusal as FaceRefusal;

/// Why a sampled trim could not be triangulated.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum TrimRefusal {
    /// Nonfinite bounds, invalid period or zero grid detail.
    InvalidChart,
    /// The face has a trim that is not a closed polygon or paired band.
    UnsupportedTrim,
    /// Trim topology cannot be triangulated without changing a shared chord.
    InvalidTrim(&'static str),
    /// The requested triangulation exceeds its size limit.
    MeshSize,
}

/// An interior grid sample, rather than a sample from an input edge.
pub const NO_BOUNDARY: u32 = u32::MAX;

/// UV in the grid's normalized rectangle, with an optional shared edge sample.
#[derive(Clone, Copy)]
pub struct TrimVertex {
    /// Coordinates relative to [`Triangulated::bounds`].
    pub uv: [f64; 2],
    /// Edge sample index, or [`NO_BOUNDARY`].
    pub boundary: u32,
}

// In units of a chart span or period, not millimetres.
const UV_EPSILON: f64 = 1e-10;

type EdgeRange<'a> = &'a dyn Fn(usize) -> Option<Range<u32>>;
type SamePosition<'a> = &'a dyn Fn(u32, u32) -> bool;
type Cdt = ConstrainedDelaunayTriangulation<Site>;

struct Site {
    point: Point2<f64>,
    input: u32,
}

impl HasPosition for Site {
    type Scalar = f64;

    fn position(&self) -> Point2<f64> {
        self.point
    }
}

#[derive(Clone, Copy)]
struct TrimPoint {
    uv: [f64; 2],
    boundary: u32,
    seam_twin: Option<usize>,
}

struct Ring {
    points: Vec<TrimPoint>,
    winding: [i32; 2],
}

/// Trimmed grid indices, retaining incident edges and paired seam points.
pub struct Triangulated {
    /// Parameter or boundary sample for each vertex.
    pub vertices: Vec<TrimVertex>,
    /// Triangle vertex indices, three per triangle.
    pub indices: Vec<u32>,
    /// Actual parameter rectangle containing the trim.
    pub bounds: [[f64; 2]; 2],
}

/// Triangulate the trim and insert a regular grid without splitting its chords.
/// With `snap_joints`, the shared 3D edge indices are authoritative at joins:
/// independently admitted pcurves need not agree in UV to binary64 precision.
pub fn triangulate<F>(
    face: &DisplayFace<F>,
    edge_range: EdgeRange<'_>,
    same_position: SamePosition<'_>,
    detail: [u32; 2],
    max_vertices: usize,
    max_indices: usize,
    snap_joints: bool,
) -> Result<Triangulated, FaceRefusal> {
    let span = face.bounds.map(|[lower, upper]| upper - lower);
    if !span.iter().all(|span| span.is_finite() && *span > 0.0)
        || detail.contains(&0)
        || face
            .periods
            .iter()
            .flatten()
            .any(|p| !p.is_finite() || *p <= 0.0)
    {
        return Err(FaceRefusal::InvalidChart);
    }
    let tolerance = span.map(|s| s * UV_EPSILON);
    let mut rings = Vec::with_capacity(face.loops.len());
    let mut trim_count = 0;
    for trim in &face.loops {
        let ring = collect_ring(
            trim,
            edge_range,
            same_position,
            face.periods,
            tolerance,
            max_vertices - trim_count,
            snap_joints,
        )?;
        trim_count += ring.points.len();
        rings.push(ring);
    }
    if rings.is_empty() {
        return Err(FaceRefusal::UnsupportedTrim);
    }
    let mut wrapped = rings
        .iter()
        .enumerate()
        .filter_map(|(i, ring)| (ring.winding != [0, 0]).then_some(i));
    let mut seam_count = 0;
    if let Some(a) = wrapped.next() {
        let b = wrapped.next().ok_or(FaceRefusal::UnsupportedTrim)?;
        if wrapped.next().is_some() {
            return Err(FaceRefusal::UnsupportedTrim);
        }
        let b = rings.remove(b);
        let a = rings.remove(a);
        seam_count = detail.into_iter().max().unwrap() as usize - 1;
        if seam_count.saturating_mul(2) > max_vertices - trim_count {
            return Err(FaceRefusal::MeshSize);
        }
        trim_count += seam_count * 2;
        rings.push(join_band(a, b, face.periods, seam_count)?);
    }

    let grid_count = (detail[0] as usize - 1).saturating_mul(detail[1] as usize - 1);
    if trim_count.saturating_add(grid_count) > max_vertices {
        return Err(FaceRefusal::MeshSize);
    }
    let mut bounds = [[f64::INFINITY, f64::NEG_INFINITY]; 2];
    for point in rings.iter().flat_map(|ring| &ring.points) {
        for (axis, [lower, upper]) in bounds.iter_mut().enumerate() {
            *lower = lower.min(point.uv[axis]);
            *upper = upper.max(point.uv[axis]);
        }
    }
    let span = bounds.map(|[lower, upper]| upper - lower);
    if span.iter().any(|s| !s.is_finite() || *s <= 0.0) {
        return Err(FaceRefusal::InvalidTrim("trim has zero chart extent"));
    }
    let normalized =
        |uv: [f64; 2]| std::array::from_fn(|axis| (uv[axis] - bounds[axis][0]) / span[axis]);

    let mut cdt = Cdt::new();
    let mut inputs = Vec::with_capacity(trim_count + grid_count - seam_count);
    let mut seams = vec![None; seam_count];
    let mut constraints = Vec::with_capacity(trim_count);
    let mut handles = Vec::new();
    for ring in &rings {
        handles.clear();
        for point in &ring.points {
            let uv = normalized(point.uv);
            let mut push = || {
                let index = inputs.len() as u32;
                inputs.push(TrimVertex {
                    uv,
                    boundary: point.boundary,
                });
                index
            };
            let input = match point.seam_twin {
                Some(twin) => *seams[twin].get_or_insert_with(&mut push),
                None => push(),
            };
            let before = cdt.num_vertices();
            let handle = cdt
                .insert(Site {
                    point: Point2::new(uv[0], uv[1]),
                    input,
                })
                .map_err(|_| FaceRefusal::InvalidTrim("trim point is outside the CDT range"))?;
            if cdt.num_vertices() == before {
                return Err(FaceRefusal::InvalidTrim(
                    "distinct trim samples have the same UV",
                ));
            }
            handles.push(handle);
        }
        for i in 0..handles.len() {
            constraints.push((handles[i], handles[(i + 1) % handles.len()]));
        }
    }
    for &(a, b) in &constraints {
        if !cdt.can_add_constraint(a, b) || !cdt.add_constraint(a, b) {
            return Err(FaceRefusal::InvalidTrim("trim constraints cross"));
        }
        // A third boundary point on this segment would silently split the
        // shared chord differently on this face than on its neighbour.
        if !cdt
            .get_edge_from_neighbors(a, b)
            .is_some_and(|e| e.is_constraint_edge())
        {
            return Err(FaceRefusal::InvalidTrim(
                "another trim sample splits a shared chord",
            ));
        }
    }
    for v in 1..detail[1] {
        for u in 1..detail[0] {
            let uv = [
                f64::from(u) / f64::from(detail[0]),
                f64::from(v) / f64::from(detail[1]),
            ];
            let point = Point2::new(uv[0], uv[1]);
            match cdt.locate(point) {
                PositionInTriangulation::OnVertex(_)
                | PositionInTriangulation::OutsideOfConvexHull(_) => continue,
                PositionInTriangulation::OnEdge(edge)
                    if cdt.directed_edge(edge).is_constraint_edge() =>
                {
                    continue;
                }
                _ => {}
            }
            let input = inputs.len() as u32;
            cdt.insert(Site { point, input })
                .map_err(|_| FaceRefusal::InvalidTrim("grid point is outside the CDT range"))?;
            inputs.push(TrimVertex {
                uv,
                boundary: NO_BOUNDARY,
            });
        }
    }
    let labels = label_triangles(&cdt)?;
    let count = labels.iter().filter(|&&label| label == Some(true)).count();
    if count == 0 || count > max_indices / 3 {
        return Err(if count == 0 {
            FaceRefusal::InvalidTrim("trim encloses no triangles")
        } else {
            FaceRefusal::MeshSize
        });
    }
    let mut used = vec![u32::MAX; inputs.len()];
    let mut vertices = Vec::new();
    let mut indices = Vec::with_capacity(count * 3);
    for triangle in cdt
        .inner_faces()
        .filter(|face| labels[face.index()] == Some(true))
    {
        let mut ids = triangle
            .vertices()
            .map(|vertex| vertex.data().input as usize);
        if ids[0] == ids[1] || ids[1] == ids[2] || ids[0] == ids[2] {
            return Err(FaceRefusal::InvalidTrim("triangle repeats a trim vertex"));
        }
        if face.reversed {
            ids.swap(1, 2);
        }
        indices.extend(ids.map(|id| {
            if used[id] == u32::MAX {
                used[id] = vertices.len() as u32;
                vertices.push(inputs[id]);
            }
            used[id]
        }));
    }
    Ok(Triangulated {
        vertices,
        indices,
        bounds,
    })
}

fn collect_ring(
    trim: &DisplayLoop,
    edge_range: EdgeRange<'_>,
    same_position: SamePosition<'_>,
    periods: [Option<f64>; 2],
    tolerance: [f64; 2],
    limit: usize,
    snap_joints: bool,
) -> Result<Ring, FaceRefusal> {
    if trim.vertex.is_some() || trim.uses.is_empty() {
        return Err(FaceRefusal::UnsupportedTrim);
    }
    let mut points = Vec::<TrimPoint>::new();
    for edge in &trim.uses {
        let range = edge_range(edge.edge).ok_or(FaceRefusal::InvalidTrim("missing shared edge"))?;
        let get = |i: usize| {
            let sample = &edge.samples[i];
            let boundary = u32::try_from(sample.point)
                .ok()
                .and_then(|p| range.start.checked_add(p))
                .filter(|&p| p < range.end)
                .ok_or(FaceRefusal::InvalidTrim("missing shared sample"))?;
            if !sample.uv.iter().all(|uv| uv.is_finite()) {
                return Err(FaceRefusal::InvalidTrim("nonfinite trim parameter"));
            }
            Ok(TrimPoint {
                uv: sample.uv,
                boundary,
                seam_twin: None,
            })
        };
        if edge.samples.len() < 2 {
            return Err(FaceRefusal::InvalidTrim(
                "shared edge has fewer than two samples",
            ));
        }
        let first = get(0)?;
        let shift = if let Some(last) = points.last() {
            if !same_position(last.boundary, first.boundary) {
                return Err(FaceRefusal::InvalidTrim(
                    "shared edge endpoints are not sewn",
                ));
            }
            let (_, shift) = wraps(first.uv, last.uv, periods, tolerance, snap_joints)?;
            shift
        } else {
            points.push(first);
            [0.0; 2]
        };
        if points.len().saturating_add(edge.samples.len() - 1) > limit.saturating_add(1) {
            return Err(FaceRefusal::MeshSize);
        }
        for i in 1..edge.samples.len() {
            let mut point = get(i)?;
            for axis in 0..2 {
                point.uv[axis] += shift[axis];
            }
            points.push(point);
        }
    }
    let first = points[0];
    let last = *points.last().unwrap();
    if !same_position(first.boundary, last.boundary) {
        return Err(FaceRefusal::InvalidTrim(
            "trim does not close at a shared vertex",
        ));
    }
    let (winding, _) = wraps(first.uv, last.uv, periods, tolerance, snap_joints)?;
    if winding == [0, 0] {
        points.pop();
    }
    if points.len() < 3 {
        return Err(FaceRefusal::InvalidTrim("trim has fewer than three nodes"));
    }
    if points.len() > limit {
        return Err(FaceRefusal::MeshSize);
    }
    Ok(Ring { points, winding })
}

fn wraps(
    from: [f64; 2],
    to: [f64; 2],
    periods: [Option<f64>; 2],
    tolerance: [f64; 2],
    snap_joints: bool,
) -> Result<([i32; 2], [f64; 2]), FaceRefusal> {
    let mut turns = [0; 2];
    let mut shift = [0.0; 2];
    for axis in 0..2 {
        let d = to[axis] - from[axis];
        if !d.is_finite() {
            return Err(FaceRefusal::InvalidTrim("nonfinite trim join"));
        }
        if let Some(period) = periods[axis] {
            let n = (d / period).round();
            if n.abs() > f64::from(i32::MAX) {
                return Err(FaceRefusal::InvalidTrim("trim join exceeds period range"));
            }
            turns[axis] = n as i32;
            shift[axis] = n * period;
        }
        if !snap_joints && (d - shift[axis]).abs() > tolerance[axis] {
            return Err(FaceRefusal::InvalidTrim(
                "trim parameters do not meet at a joint",
            ));
        }
    }
    Ok((turns, shift))
}

// A periodic band has two rings of opposite winding. Cut it at paired
// endpoints and give each artificial seam point one vertex index on both sides.
fn join_band(
    a: Ring,
    mut b: Ring,
    periods: [Option<f64>; 2],
    seam_count: usize,
) -> Result<Ring, FaceRefusal> {
    let axis = usize::from(a.winding[0] == 0);
    if a.winding[axis].abs() != 1
        || a.winding[1 - axis] != 0
        || b.winding[axis] != -a.winding[axis]
        || b.winding[1 - axis] != 0
    {
        return Err(FaceRefusal::UnsupportedTrim);
    }
    let period = periods[axis].ok_or(FaceRefusal::InvalidChart)?;
    let a_last = *a.points.last().unwrap();
    b.points.pop();
    let distance = |p: &TrimPoint| {
        let d = (p.uv[axis] - a_last.uv[axis]) / period;
        (d - d.round()).abs()
    };
    let start = b
        .points
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| distance(a).total_cmp(&distance(b)))
        .unwrap()
        .0;
    let shift = ((a_last.uv[axis] - b.points[start].uv[axis]) / period).round() * period;
    let turn = f64::from(b.winding[axis]) * period;
    for (i, point) in b.points.iter_mut().enumerate() {
        point.uv[axis] += shift + if i < start { turn } else { 0.0 };
    }
    b.points.rotate_left(start);
    let mut b_last = b.points[0];
    b_last.uv[axis] += turn;
    let mut points = a.points;
    let a_first = points[0];
    let b_first = b.points[0];
    let interpolate = |p: TrimPoint, q: TrimPoint, twin: usize| {
        let t = (twin + 1) as f64 / (seam_count + 1) as f64;
        TrimPoint {
            uv: std::array::from_fn(|axis| p.uv[axis] + (q.uv[axis] - p.uv[axis]) * t),
            boundary: NO_BOUNDARY,
            seam_twin: Some(twin),
        }
    };
    points.extend((0..seam_count).map(|i| interpolate(a_last, b_first, i)));
    points.extend(b.points);
    points.push(b_last);
    points.extend(
        (0..seam_count)
            .rev()
            .map(|i| interpolate(a_first, b_last, i)),
    );
    Ok(Ring {
        points,
        winding: [0, 0],
    })
}

// Crossing a constraint toggles the inside label, including at holes. Visit
// every triangle once, rather than testing every centroid against every trim.
fn label_triangles(cdt: &Cdt) -> Result<Vec<Option<bool>>, FaceRefusal> {
    let mut labels = vec![None; cdt.num_all_faces()];
    let mut queue = VecDeque::new();
    for triangle in cdt.inner_faces() {
        for edge in triangle.adjacent_edges() {
            if edge.rev().face().as_inner().is_none() {
                mark(
                    &mut labels[triangle.index()],
                    edge.is_constraint_edge(),
                    || queue.push_back(triangle.fix()),
                )?;
            }
        }
    }
    while let Some(face) = queue.pop_front() {
        let inside = labels[face.index()].unwrap();
        for edge in cdt.face(face).adjacent_edges() {
            if let Some(neighbour) = edge.rev().face().as_inner() {
                mark(
                    &mut labels[neighbour.index()],
                    inside ^ edge.is_constraint_edge(),
                    || queue.push_back(neighbour.fix()),
                )?;
            }
        }
    }
    Ok(labels)
}

fn mark(label: &mut Option<bool>, inside: bool, queue: impl FnOnce()) -> Result<(), FaceRefusal> {
    match label {
        Some(old) if *old != inside => {
            Err(FaceRefusal::InvalidTrim("trim has inconsistent winding"))
        }
        Some(_) => Ok(()),
        None => {
            *label = Some(inside);
            queue();
            Ok(())
        }
    }
}
