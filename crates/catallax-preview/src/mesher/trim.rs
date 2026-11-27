use std::collections::VecDeque;
use std::ops::Range;

use catallax_mesh::float::{DisplayFace, DisplayLoop, DisplaySurface};
use spade::{
    ConstrainedDelaunayTriangulation, HasPosition, Point2, PositionInTriangulation, Triangulation,
};

use super::{FaceRefusal, InputVertex, NO_BOUNDARY, position32};

type EdgeRanges = [Result<Range<u32>, FaceRefusal>];
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

pub(super) struct Triangulated {
    pub vertices: Vec<InputVertex>,
    pub indices: Vec<u32>,
    pub bounds: [[f64; 2]; 2],
}

pub(super) fn triangulate<F>(
    face: &DisplayFace<F>,
    edges: &EdgeRanges,
    positions: &mut Vec<[f32; 4]>,
    max_positions: usize,
    detail: [u32; 2],
    max_vertices: usize,
    max_indices: usize,
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
    let mut rings = Vec::with_capacity(face.loops.len());
    let mut trim_count = 0;
    let mut poles = Vec::new();
    for trim in &face.loops {
        if let Some((point, uv)) = trim.vertex {
            if !trim.uses.is_empty() || uv.iter().any(|x| !x.is_finite()) {
                return Err(FaceRefusal::InvalidTrim);
            }
            poles.push((point, uv));
            continue;
        }
        let ring = collect_ring(
            trim,
            edges,
            positions,
            face.periods,
            max_vertices - trim_count,
        )?;
        trim_count += ring.points.len();
        rings.push(ring);
    }
    let mut seam_count = 0;
    let pole_start = positions.len() as u32;
    let wrapped: Vec<_> = rings
        .iter()
        .enumerate()
        .filter_map(|(i, ring)| (ring.winding != [0, 0]).then_some(i))
        .collect();
    if !poles.is_empty()
        || matches!(face.surface, DisplaySurface::Sphere { .. }) && wrapped.len() < 2
    {
        let needed = if rings.is_empty() { 2 } else { 1 };
        if (rings.is_empty() || wrapped.len() == 1)
            && poles.len() < needed
            && let DisplaySurface::Sphere { frame, radius } = face.surface
        {
            let direction = wrapped.first().map_or(-1, |&i| {
                rings[i].winding[0] * if face.reversed { -1 } else { 1 }
            });
            for sign in [direction, -direction]
                .into_iter()
                .take(needed - poles.len())
            {
                let sign = f64::from(sign);
                poles.push((
                    std::array::from_fn(|i| frame.origin[i] + sign * radius * frame.z[i]),
                    [face.bounds[0][0], sign * 0.25],
                ));
            }
        }
    }
    if !poles.is_empty() {
        let Some(period) = face.periods[0] else {
            return Err(FaceRefusal::UnsupportedTrim);
        };
        if face.periods[1].is_some()
            || !matches!((wrapped.as_slice(), poles.len()), ([_], 1) | ([], 2))
            || wrapped.is_empty() && !rings.is_empty()
        {
            return Err(FaceRefusal::UnsupportedTrim);
        }
        if poles.len() > max_positions - positions.len() {
            return Err(FaceRefusal::MeshSize);
        }
        seam_count = detail[1] as usize - 1;
        let mut pole_rings = poles.iter().enumerate().map(|(i, &(p, uv))| {
            let boundary = positions.len() as u32;
            positions.push(position32(p)?);
            Ok((uv, boundary, seam_count + i))
        });
        let a = if let Some(&i) = wrapped.first() {
            rings.remove(i)
        } else {
            let (uv, boundary, twin) = pole_rings.next().unwrap()?;
            let direction = if uv[1] < poles[1].1[1] { 1 } else { -1 };
            pole_ring(uv, boundary, twin, direction, period)
        };
        if a.winding[0].abs() != 1 || a.winding[1] != 0 {
            return Err(FaceRefusal::UnsupportedTrim);
        }
        let (mut uv, boundary, twin) = pole_rings.next().unwrap()?;
        uv[0] = a.points.last().unwrap().uv[0];
        let b = pole_ring(uv, boundary, twin, -a.winding[0], period);
        trim_count = trim_count.saturating_add(poles.len() * 2);
        rings.push(join_band(a, b, face.periods, seam_count)?);
    } else {
        match wrapped.as_slice() {
            [] if rings.is_empty() => return Err(FaceRefusal::UnsupportedTrim),
            [] => (),
            &[a, b] => {
                let b = rings.remove(b);
                let a = rings.remove(a);
                seam_count = detail.into_iter().max().unwrap() as usize - 1;
                rings.push(join_band(a, b, face.periods, seam_count)?);
            }
            _ => return Err(FaceRefusal::UnsupportedTrim),
        }
    }
    trim_count = trim_count.saturating_add(seam_count.saturating_mul(2));

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
        return Err(FaceRefusal::InvalidTrim);
    }
    let normalized =
        |uv: [f64; 2]| std::array::from_fn(|axis| (uv[axis] - bounds[axis][0]) / span[axis]);
    // The outer coedge loop gives the draw orientation. In charts whose
    // parameterization flips, a face's geometric sense alone does not.
    let outer_area = rings
        .iter()
        .map(|ring| {
            let uv = ring.points.iter().map(|p| normalized(p.uv));
            uv.clone()
                .zip(uv.cycle().skip(1))
                .map(|([u0, v0], [u1, v1])| u0 * v1 - u1 * v0)
                .sum::<f64>()
        })
        .max_by(|a, b| a.abs().total_cmp(&b.abs()))
        .ok_or(FaceRefusal::InvalidTrim)?;
    let reversed = outer_area < 0.0;

    let mut cdt = Cdt::new();
    let mut inputs = Vec::with_capacity(trim_count + grid_count - seam_count);
    let mut seams = vec![None; seam_count + poles.len()];
    let mut constraints = Vec::with_capacity(trim_count);
    let mut handles = Vec::new();
    for ring in &rings {
        handles.clear();
        for point in &ring.points {
            let uv = normalized(point.uv);
            let mut push = || {
                let index = inputs.len() as u32;
                inputs.push(InputVertex::new(uv, point.boundary));
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
                .map_err(|_| FaceRefusal::InvalidTrim)?;
            if cdt.num_vertices() == before {
                trim_diagnostic("duplicate projected vertex", [point.uv, uv]);
                return Err(FaceRefusal::InvalidTrim);
            }
            handles.push(handle);
        }
        for i in 0..handles.len() {
            constraints.push((handles[i], handles[(i + 1) % handles.len()]));
        }
    }
    for &(a, b) in &constraints {
        if !cdt.can_add_constraint(a, b) || !cdt.add_constraint(a, b) {
            trim_diagnostic(
                "crossing constraint",
                [cdt.vertex(a).position(), cdt.vertex(b).position()],
            );
            return Err(FaceRefusal::InvalidTrim);
        }
        // A third boundary point on this segment would silently split the
        // shared chord differently on this face than on its neighbour.
        if !cdt
            .get_edge_from_neighbors(a, b)
            .is_some_and(|e| e.is_constraint_edge())
        {
            trim_diagnostic(
                "split constraint",
                [cdt.vertex(a).position(), cdt.vertex(b).position()],
            );
            return Err(FaceRefusal::InvalidTrim);
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
                .map_err(|_| FaceRefusal::InvalidTrim)?;
            inputs.push(InputVertex::new(uv, NO_BOUNDARY));
        }
    }
    let labels = label_triangles(&cdt)?;
    let count = labels.iter().filter(|&&label| label == Some(true)).count();
    if count == 0 || count > max_indices / 3 {
        return Err(if count == 0 {
            FaceRefusal::InvalidTrim
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
        if (0..3).any(|i| ids[i] == ids[(i + 1) % 3]) {
            if (0..3).any(|i| {
                ids[i] == ids[(i + 1) % 3]
                    && (pole_start..pole_start + poles.len() as u32)
                        .contains(&inputs[ids[i]].boundary)
            }) {
                continue;
            }
            return Err(FaceRefusal::InvalidTrim);
        }
        if reversed {
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
    if indices.is_empty() {
        return Err(FaceRefusal::InvalidTrim);
    }
    Ok(Triangulated {
        vertices,
        indices,
        bounds,
    })
}

fn trim_diagnostic(message: &str, points: impl std::fmt::Debug) {
    if std::env::var_os("CATALLAX_TRIM_DIAGNOSTICS").is_some() {
        eprintln!("{message}: {points:?}");
    }
}

fn collect_ring(
    trim: &DisplayLoop,
    edges: &EdgeRanges,
    positions: &[[f32; 4]],
    periods: [Option<f64>; 2],
    limit: usize,
) -> Result<Ring, FaceRefusal> {
    if trim.vertex.is_some() || trim.uses.is_empty() {
        return Err(FaceRefusal::UnsupportedTrim);
    }
    let mut points = Vec::<TrimPoint>::new();
    for edge in &trim.uses {
        let range = edges
            .get(edge.edge)
            .ok_or(FaceRefusal::InvalidTrim)?
            .as_ref()
            .map_err(|e| *e)?;
        let get = |i: usize| {
            let sample = &edge.samples[i];
            let boundary = u32::try_from(sample.point)
                .ok()
                .and_then(|p| range.start.checked_add(p))
                .filter(|&p| p < range.end)
                .ok_or(FaceRefusal::InvalidTrim)?;
            if !sample.uv.iter().all(|uv| uv.is_finite()) {
                return Err(FaceRefusal::InvalidTrim);
            }
            Ok(TrimPoint {
                uv: sample.uv,
                boundary,
                seam_twin: None,
            })
        };
        if edge.samples.len() < 2 {
            return Err(FaceRefusal::InvalidTrim);
        }
        let first = get(0)?;
        let shift = if let Some(last) = points.last() {
            // The shared point is authoritative even when reprojection chose
            // different coordinates along a collapsed chart boundary.
            same_position(*last, first, positions)?;
            periodic_shift(first.uv, last.uv, periods)?.1
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
    same_position(first, last, positions)?;
    let (winding, shift) = periodic_shift(first.uv, last.uv, periods)?;
    if winding == [0, 0] {
        points.pop();
    } else {
        points.last_mut().unwrap().uv = std::array::from_fn(|i| first.uv[i] + shift[i]);
    }
    if points.len() < 3 {
        return Err(FaceRefusal::InvalidTrim);
    }
    if points.len() > limit {
        return Err(FaceRefusal::MeshSize);
    }
    Ok(Ring { points, winding })
}

fn same_position(a: TrimPoint, b: TrimPoint, positions: &[[f32; 4]]) -> Result<(), FaceRefusal> {
    if positions[a.boundary as usize].map(f32::to_bits)
        == positions[b.boundary as usize].map(f32::to_bits)
    {
        Ok(())
    } else {
        Err(FaceRefusal::InvalidTrim)
    }
}

fn periodic_shift(
    from: [f64; 2],
    to: [f64; 2],
    periods: [Option<f64>; 2],
) -> Result<([i32; 2], [f64; 2]), FaceRefusal> {
    let mut turns = [0; 2];
    let mut shift = [0.0; 2];
    for axis in 0..2 {
        let d = to[axis] - from[axis];
        if !d.is_finite() {
            return Err(FaceRefusal::InvalidTrim);
        }
        if let Some(period) = periods[axis] {
            let n = (d / period).round();
            if n.abs() > f64::from(i32::MAX) {
                return Err(FaceRefusal::InvalidTrim);
            }
            turns[axis] = n as i32;
            shift[axis] = n * period;
        }
    }
    Ok((turns, shift))
}

fn pole_ring(uv: [f64; 2], boundary: u32, twin: usize, direction: i32, period: f64) -> Ring {
    let first = TrimPoint {
        uv,
        boundary,
        seam_twin: Some(twin),
    };
    let mut last = first;
    last.uv[0] += f64::from(direction) * period;
    Ring {
        points: vec![first, last],
        winding: [direction, 0],
    }
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
        Some(old) if *old != inside => Err(FaceRefusal::InvalidTrim),
        Some(_) => Ok(()),
        None => {
            *label = Some(inside);
            queue();
            Ok(())
        }
    }
}
