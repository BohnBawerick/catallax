use std::collections::VecDeque;
use std::ops::Range;

use super::{DisplayEdgeUse, DisplayFace, DisplayFrame, DisplaySurface};
use spade::{
    ConstrainedDelaunayTriangulation, HasPosition, Point2, PositionInTriangulation, Triangulation,
};

/// The reason a sampled surface trim could not be sewn.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum TrimRefusal {
    /// Chart bounds, periods, or detail are invalid.
    InvalidChart,
    /// The parameter layout cannot be cut into a finite chart.
    UnsupportedTrim,
    /// The sampled trim is not a sewn polygon in the selected chart.
    InvalidTrim(&'static str),
    /// The triangulation exceeds its size limit.
    MeshSize,
    /// A generated pole position cannot be represented.
    InvalidPosition,
}

/// An interior grid vertex, without a shared boundary index.
pub const NO_BOUNDARY: u32 = u32::MAX;

/// One unproved sample in normalized grid coordinates.
#[derive(Clone, Copy)]
pub struct TrimVertex {
    /// Normalized coordinates within the returned bounds.
    pub uv: [f64; 2],
    /// Shared boundary index, or NO_BOUNDARY for an interior sample.
    pub boundary: u32,
}

/// Shared trim samples in the consumer's coordinate precision.
pub trait TrimBoundaries {
    /// Contiguous indices corresponding to one edge's samples.
    fn edge_range(&self, edge: usize) -> Option<Range<u32>>;
    /// Whether the non-interior indices designate the same position.
    fn same_position(&self, a: u32, b: u32) -> bool;
    /// Rounded position of one shared sample.
    fn position(&self, sample: u32) -> [f64; 3];
    /// Add a pole used by a singular chart boundary.
    fn add_pole(&mut self, position: [f64; 3]) -> Result<u32, TrimRefusal>;
}
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

/// Trimmed grid shared by CPU and GPU preview consumers.
pub struct Triangulated {
    /// Normalized UV or shared boundary reference for each vertex.
    pub vertices: Vec<TrimVertex>,
    /// Triangle vertex indices in groups of three.
    pub indices: Vec<u32>,
    /// Finite UV bounds of the triangulated region.
    pub bounds: [[f64; 2]; 2],
    /// Stereographic frame when a spherical trim has been recharted.
    pub sphere_chart: Option<DisplayFrame>,
}

/// Sew a trim and insert grid stations without splitting shared chords.
pub fn triangulate<F>(
    face: &DisplayFace<F>,
    boundary: &mut impl TrimBoundaries,
    detail: [u32; 2],
    max_vertices: usize,
    max_indices: usize,
) -> Result<Triangulated, TrimRefusal> {
    let span = face.bounds.map(|[lower, upper]| upper - lower);
    if !span.iter().all(|span| span.is_finite() && *span > 0.0)
        || detail.contains(&0)
        || face
            .periods
            .iter()
            .flatten()
            .any(|p| !p.is_finite() || *p <= 0.0)
    {
        return Err(TrimRefusal::InvalidChart);
    }
    let mut rings = Vec::with_capacity(face.loops.len());
    let mut trim_count = 0;
    let mut poles = Vec::new();
    for trim in &face.loops {
        if let Some((point, uv)) = trim.vertex {
            if !trim.uses.is_empty() || uv.iter().any(|x| !x.is_finite()) {
                return Err(TrimRefusal::InvalidTrim("invalid pole trim"));
            }
            poles.push((point, uv));
            continue;
        }
        if trim.uses.is_empty() {
            return Err(TrimRefusal::UnsupportedTrim);
        }
        let split = trim.uses.len() > 1 && trim.uses.iter().all(|edge| closed_use(edge, boundary));
        let width = if split { 1 } else { trim.uses.len() };
        for uses in trim.uses.chunks(width) {
            let ring = collect_ring(
                uses,
                boundary,
                face.periods,
                span,
                max_vertices - trim_count,
            )?;
            trim_count += ring.points.len();
            rings.push(ring);
        }
    }
    let sphere_chart = match face.surface {
        DisplaySurface::Sphere { frame, radius }
            if !rings.is_empty()
                && face.periods[0].is_some()
                && rings.iter().filter(|r| r.winding != [0, 0]).count() < 2 =>
        {
            spherical_chart(frame, radius, face.reversed, &mut rings, boundary)
        }
        _ => None,
    };
    if sphere_chart.is_some() {
        poles.clear();
    }
    let mut seam_count = 0;
    let wrapped: Vec<_> = rings
        .iter()
        .enumerate()
        .filter_map(|(i, ring)| (ring.winding != [0, 0]).then_some(i))
        .collect();
    if sphere_chart.is_none()
        && (!poles.is_empty()
            || matches!(face.surface, DisplaySurface::Sphere { .. }) && wrapped.len() < 2)
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
            return Err(TrimRefusal::UnsupportedTrim);
        };
        if face.periods[1].is_some()
            || !matches!((wrapped.as_slice(), poles.len()), ([_], 1) | ([], 2))
            || wrapped.is_empty() && !rings.is_empty()
        {
            return Err(TrimRefusal::UnsupportedTrim);
        }
        seam_count = detail[1] as usize - 1;
        let mut pole_rings = poles
            .iter()
            .enumerate()
            .map(|(i, &(p, uv))| Ok((uv, boundary.add_pole(p)?, seam_count + i)));
        let a = if let Some(&i) = wrapped.first() {
            rings.remove(i)
        } else {
            let (uv, boundary, twin) = pole_rings.next().unwrap()?;
            let direction = if uv[1] < poles[1].1[1] { 1 } else { -1 };
            pole_ring(uv, boundary, twin, direction, period)
        };
        if a.winding[0].abs() != 1 || a.winding[1] != 0 {
            return Err(TrimRefusal::UnsupportedTrim);
        }
        let (mut uv, boundary, twin) = pole_rings.next().unwrap()?;
        uv[0] = a.points.last().unwrap().uv[0];
        let b = pole_ring(uv, boundary, twin, -a.winding[0], period);
        trim_count = trim_count.saturating_add(poles.len() * 2);
        let band = join_band(a, b, &mut rings, face.periods, face.reversed, seam_count)?;
        rings.push(band);
    } else {
        match wrapped.as_slice() {
            [] if rings.is_empty() => return Err(TrimRefusal::UnsupportedTrim),
            [] => (),
            &[a, b] => {
                let b = rings.remove(b);
                let a = rings.remove(a);
                seam_count = detail.into_iter().max().unwrap() as usize - 1;
                let band = join_band(a, b, &mut rings, face.periods, face.reversed, seam_count)?;
                rings.push(band);
            }
            _ => return Err(TrimRefusal::UnsupportedTrim),
        }
    }
    trim_count = trim_count.saturating_add(seam_count.saturating_mul(2));

    let grid_count = (detail[0] as usize - 1).saturating_mul(detail[1] as usize - 1);
    if trim_count.saturating_add(grid_count) > max_vertices {
        return Err(TrimRefusal::MeshSize);
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
        trim_diagnostic("empty normalized chart", bounds);
        return Err(TrimRefusal::InvalidTrim("trim has zero chart extent"));
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
        .ok_or(TrimRefusal::InvalidTrim("invalid trim chart scale"))?;
    let reversed = outer_area < 0.0;

    let mut cdt = Cdt::new();
    let mut inputs: Vec<TrimVertex> = Vec::with_capacity(trim_count + grid_count - seam_count);
    let mut seams = vec![None; seam_count + poles.len()];
    let mut constraints = Vec::with_capacity(trim_count);
    let mut handles = Vec::new();
    for ring in &rings {
        handles.clear();
        for (i, point) in ring.points.iter().enumerate() {
            let uv = normalized(point.uv);
            let mut site = spade::mitigate_underflow(Point2::new(uv[0], uv[1]));
            if let PositionInTriangulation::OnVertex(existing) = cdt.locate(site) {
                let other = inputs[cdt.vertex(existing).data().input as usize].boundary;
                if same_boundary(point.boundary, other, boundary) {
                    if handles.last() != Some(&existing) {
                        handles.push(existing);
                    }
                    continue;
                }
                // Keep the 3D chord authoritative when two distinct samples
                // project to the same UV at a nearly collapsed chart end.
                let next = normalized(ring.points[(i + 1) % ring.points.len()].uv);
                let previous =
                    normalized(ring.points[(i + ring.points.len() - 1) % ring.points.len()].uv);
                let dx = next[0].midpoint(previous[0]) - site.x;
                let dy = next[1].midpoint(previous[1]) - site.y;
                let scale = dx.abs().max(dy.abs());
                if scale <= f64::MIN_POSITIVE {
                    trim_diagnostic(
                        "coincident projected trim neighbourhood",
                        (site, previous, next),
                    );
                    return Err(TrimRefusal::InvalidTrim(
                        "distinct trim samples have no separating UV",
                    ));
                }
                site = spade::mitigate_underflow(Point2::new(
                    site.x + 1e-11 * dx / scale,
                    site.y + 1e-11 * dy / scale,
                ));
            }
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
            let handle = cdt.insert(Site { point: site, input }).map_err(|error| {
                trim_diagnostic("invalid triangulation coordinate", (site, error));
                TrimRefusal::InvalidTrim("trim point is outside the CDT range")
            })?;
            if cdt.num_vertices() == before {
                trim_diagnostic("duplicate projected vertex", [point.uv, uv]);
                return Err(TrimRefusal::InvalidTrim(
                    "trim samples share a projected vertex",
                ));
            }
            handles.push(handle);
        }
        if handles.first() == handles.last() {
            handles.pop();
        }
        if handles.len() < 3 {
            return Err(TrimRefusal::InvalidTrim(
                "trim generated no chord constraints",
            ));
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
            return Err(TrimRefusal::InvalidTrim("trim constraints cross"));
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
            return Err(TrimRefusal::InvalidTrim(
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
                .map_err(|_| TrimRefusal::InvalidTrim("grid point is outside the CDT range"))?;
            inputs.push(TrimVertex {
                uv,
                boundary: NO_BOUNDARY,
            });
        }
    }
    let labels =
        label_triangles(&cdt).inspect_err(|_| trim_diagnostic("ambiguous trim parity", bounds))?;
    let count = labels.iter().filter(|&&label| label == Some(true)).count();
    if count == 0 || count > max_indices / 3 {
        trim_diagnostic("no interior triangles or index limit", (bounds, count));
        return Err(if count == 0 {
            TrimRefusal::InvalidTrim("trim encloses no triangles")
        } else {
            TrimRefusal::MeshSize
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
        if (0..3).any(|i| {
            let a = ids[i];
            let b = ids[(i + 1) % 3];
            a == b || same_boundary(inputs[a].boundary, inputs[b].boundary, boundary)
        }) {
            continue;
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
        trim_diagnostic("all triangles collapsed", bounds);
        return Err(TrimRefusal::InvalidTrim("all sampled triangles collapsed"));
    }
    Ok(Triangulated {
        vertices,
        indices,
        bounds,
        sphere_chart,
    })
}

fn trim_diagnostic(message: &str, points: impl std::fmt::Debug) {
    if std::env::var_os("CATALLAX_TRIM_DIAGNOSTICS").is_some() {
        eprintln!("{message}: {points:?}");
    }
}

fn closed_use(edge: &DisplayEdgeUse, boundary: &impl TrimBoundaries) -> bool {
    let Some(range) = boundary.edge_range(edge.edge) else {
        return false;
    };
    let sample_index = |i: usize| {
        u32::try_from(i)
            .ok()
            .and_then(|i| range.start.checked_add(i))
            .filter(|&i| i < range.end)
    };
    edge.samples
        .first()
        .zip(edge.samples.last())
        .is_some_and(|(a, b)| {
            sample_index(a.point)
                .zip(sample_index(b.point))
                .is_some_and(|(a, b)| same_boundary(a, b, boundary))
        })
}

fn collect_ring(
    uses: &[DisplayEdgeUse],
    boundary: &impl TrimBoundaries,
    periods: [Option<f64>; 2],
    span: [f64; 2],
    limit: usize,
) -> Result<Ring, TrimRefusal> {
    let mut points = Vec::<TrimPoint>::new();
    for edge in uses {
        let range = boundary
            .edge_range(edge.edge)
            .ok_or(TrimRefusal::InvalidTrim("missing shared edge"))?;
        let get = |i: usize| {
            let sample = &edge.samples[i];
            let boundary = u32::try_from(sample.point)
                .ok()
                .and_then(|p| range.start.checked_add(p))
                .filter(|&p| p < range.end)
                .ok_or(TrimRefusal::InvalidTrim("missing shared edge sample"))?;
            if !sample.uv.iter().all(|uv| uv.is_finite()) {
                return Err(TrimRefusal::InvalidTrim("nonfinite trim parameter"));
            }
            Ok(TrimPoint {
                uv: sample.uv,
                boundary,
                seam_twin: None,
            })
        };
        if edge.samples.len() < 2 {
            return Err(TrimRefusal::InvalidTrim(
                "shared edge has fewer than two samples",
            ));
        }
        let first = get(0)?;
        let shift = if let Some(last) = points.last() {
            // The shared point is authoritative even when reprojection chose
            // different coordinates along a collapsed chart boundary.
            same_position(*last, first, boundary)?;
            periodic_shift(first.uv, last.uv, periods)?.1
        } else {
            points.push(first);
            [0.0; 2]
        };
        if points.len().saturating_add(edge.samples.len() - 1) > limit.saturating_add(1) {
            return Err(TrimRefusal::MeshSize);
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
    same_position(first, last, boundary)?;
    let (winding, shift) = periodic_shift(first.uv, last.uv, periods)?;
    let lifted: [f64; 2] = std::array::from_fn(|i| first.uv[i] + shift[i]);
    let same_uv = (0..2).all(|i| (last.uv[i] - lifted[i]).abs() <= 1e-9 * span[i]);
    if same_uv {
        if winding == [0, 0] {
            points.pop();
        } else {
            points.last_mut().unwrap().uv = lifted;
        }
    } else if winding != [0, 0] {
        points.push(TrimPoint {
            uv: lifted,
            ..first
        });
    }
    if points.len() < 3 {
        return Err(TrimRefusal::InvalidTrim("trim has fewer than three nodes"));
    }
    if points.len() > limit {
        return Err(TrimRefusal::MeshSize);
    }
    Ok(Ring { points, winding })
}

fn same_position(
    a: TrimPoint,
    b: TrimPoint,
    boundary: &impl TrimBoundaries,
) -> Result<(), TrimRefusal> {
    if same_boundary(a.boundary, b.boundary, boundary) {
        Ok(())
    } else {
        trim_diagnostic("discontinuous shared trim point", (a.uv, b.uv));
        Err(TrimRefusal::InvalidTrim(
            "trim does not close at a shared vertex",
        ))
    }
}

fn same_boundary(a: u32, b: u32, boundary: &impl TrimBoundaries) -> bool {
    a != NO_BOUNDARY && b != NO_BOUNDARY && boundary.same_position(a, b)
}

fn periodic_shift(
    from: [f64; 2],
    to: [f64; 2],
    periods: [Option<f64>; 2],
) -> Result<([i32; 2], [f64; 2]), TrimRefusal> {
    let mut turns = [0; 2];
    let mut shift = [0.0; 2];
    for axis in 0..2 {
        let d = to[axis] - from[axis];
        if !d.is_finite() {
            return Err(TrimRefusal::InvalidTrim("nonfinite trim join"));
        }
        if let Some(period) = periods[axis] {
            let n = (d / period).round();
            if n.abs() > f64::from(i32::MAX) {
                return Err(TrimRefusal::InvalidTrim("trim period is outside range"));
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

// A geographic pole can be inside a face or on its boundary. Project from
// the antipode of the oriented patch instead, so the pole has no special status.
fn spherical_chart(
    frame: DisplayFrame,
    radius: f64,
    reversed: bool,
    rings: &mut [Ring],
    boundary: &impl TrimBoundaries,
) -> Option<DisplayFrame> {
    let dot = |a: [f64; 3], b: [f64; 3]| a.into_iter().zip(b).map(|(a, b)| a * b).sum::<f64>();
    let cross =
        |[x, y, z]: [f64; 3], [u, v, w]: [f64; 3]| [y * w - z * v, z * u - x * w, x * v - y * u];
    let unit = |a: [f64; 3]| {
        let length = dot(a, a).sqrt();
        (length.is_finite() && length > 0.0).then(|| a.map(|x| x / length))
    };
    let radial = |point: &TrimPoint| {
        let p = boundary.position(point.boundary);
        std::array::from_fn(|i| (p[i] - frame.origin[i]) / radius)
    };
    let mut area = [0.0; 3];
    for ring in &*rings {
        let mut previous = radial(ring.points.last()?);
        for point in &ring.points {
            let next = radial(point);
            let contribution = cross(previous, next);
            for i in 0..3 {
                area[i] += contribution[i];
            }
            previous = next;
        }
    }
    let z = unit(if reversed { area.map(|x| -x) } else { area })?;
    let reference = if z[0].abs() < 0.8 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    let y = unit(cross(z, reference))?;
    let x = cross(y, z);
    if rings.iter().flat_map(|ring| &ring.points).any(|p| {
        let r = radial(p);
        let length_sq = dot(r, r);
        !length_sq.is_finite() || length_sq <= f64::MIN_POSITIVE || dot(r, z) < -0.99
    }) {
        return None;
    }
    for ring in rings {
        if ring.winding != [0, 0] {
            ring.points.pop();
            ring.winding = [0, 0];
        }
        for p in &mut ring.points {
            let r = unit(radial(p)).unwrap();
            let scale = 1.0 / (1.0 + dot(r, z));
            p.uv = [dot(r, x) * scale, dot(r, y) * scale];
        }
    }
    Some(DisplayFrame {
        origin: frame.origin,
        x,
        y,
        z,
    })
}

// A periodic band has two rings of opposite winding. Cut it at paired
// endpoints and give each artificial seam point one vertex index on both sides.
fn join_band(
    mut a: Ring,
    mut b: Ring,
    holes: &mut [Ring],
    periods: [Option<f64>; 2],
    reversed: bool,
    seam_count: usize,
) -> Result<Ring, TrimRefusal> {
    let axis = usize::from(a.winding[0] == 0);
    if a.winding[axis].abs() != 1
        || a.winding[1 - axis] != 0
        || b.winding[axis] != -a.winding[axis]
        || b.winding[1 - axis] != 0
    {
        return Err(TrimRefusal::UnsupportedTrim);
    }
    let period = periods[axis].ok_or(TrimRefusal::InvalidChart)?;
    let across = 1 - axis;
    let direction = f64::from(a.winding[axis]);
    let start = seam_start(&a, holes, axis, period);
    rotate_ring(&mut a, start, axis, period, 0.0);
    let a_last = *a.points.last().unwrap();
    if b.points.len() == 2 && b.points[0].seam_twin.is_some() {
        b.points[0].uv[axis] = a_last.uv[axis];
    }
    // Score B's canonical points: the repeated endpoint is the first point.
    let distance = |p: &TrimPoint| {
        let d = (p.uv[axis] - a_last.uv[axis]) / period;
        (d - d.round()).abs()
    };
    let start = b.points[..b.points.len() - 1]
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| distance(a).total_cmp(&distance(b)))
        .unwrap()
        .0;
    let phase = ((a_last.uv[axis] - b.points[start].uv[axis]) / period).round() * period;
    rotate_ring(&mut b, start, axis, period, phase);
    if let Some(across_period) = periods[across] {
        // The opposite rail has two lifts on a doubly periodic carrier. Pick
        // the one on the interior side of the oriented first rail.
        let sign = direction * if (axis == 1) ^ reversed { -1.0 } else { 1.0 };
        let delta = b.points[0].uv[across] - a_last.uv[across];
        let offset = sign * (sign * delta).rem_euclid(across_period) - delta;
        for p in &mut b.points {
            p.uv[across] += offset;
        }
    }
    let b_last = b.points.pop().unwrap();
    let mut points = a.points;
    let a_first = points[0];
    let b_first = b.points[0];
    let middle = a_first.uv[axis].midpoint(a_last.uv[axis]);
    for hole in holes {
        let anchor = hole.points[0].uv[axis];
        let shift = ((middle - anchor) / period).round() * period;
        for p in &mut hole.points {
            p.uv[axis] += shift;
        }
    }
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

fn rotate_ring(ring: &mut Ring, start: usize, axis: usize, period: f64, shift: f64) {
    ring.points.pop();
    let turn = f64::from(ring.winding[axis]) * period;
    for (i, p) in ring.points.iter_mut().enumerate() {
        p.uv[axis] += shift + if i < start { turn } else { 0.0 };
    }
    ring.points.rotate_left(start);
    let mut last = ring.points[0];
    last.uv[axis] += turn;
    ring.points.push(last);
}

fn seam_start(ring: &Ring, holes: &[Ring], axis: usize, period: f64) -> usize {
    let intervals: Vec<_> = holes
        .iter()
        .map(|hole| {
            hole.points
                .iter()
                .fold([f64::INFINITY, f64::NEG_INFINITY], |[lo, hi], p| {
                    [lo.min(p.uv[axis]), hi.max(p.uv[axis])]
                })
        })
        .collect();
    let points = &ring.points[..ring.points.len() - 1];
    let count = points.len();
    let midpoint = count / 2;
    points
        .iter()
        .enumerate()
        .filter(|(_, p)| {
            intervals.iter().all(|&[lo, hi]| {
                let local = (p.uv[axis] - lo).rem_euclid(period);
                local > hi - lo + 0.02 * period && local < 0.98 * period
            })
        })
        .min_by_key(|(i, _)| i.abs_diff(midpoint))
        .map_or(0, |(i, _)| i)
}

// Crossing a constraint toggles the inside label, including at holes. Visit
// every triangle once, rather than testing every centroid against every trim.
pub(crate) fn label_triangles<V: HasPosition<Scalar = f64>>(
    cdt: &ConstrainedDelaunayTriangulation<V>,
) -> Result<Vec<Option<bool>>, TrimRefusal> {
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

fn mark(label: &mut Option<bool>, inside: bool, queue: impl FnOnce()) -> Result<(), TrimRefusal> {
    match label {
        Some(old) if *old != inside => {
            Err(TrimRefusal::InvalidTrim("trim has inconsistent winding"))
        }
        Some(_) => Ok(()),
        None => {
            *label = Some(inside);
            queue();
            Ok(())
        }
    }
}
