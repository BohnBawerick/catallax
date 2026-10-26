use std::collections::{BTreeMap, BTreeSet, VecDeque};

use catallax_numbers::{Answer, Answered, Decider, Interval, Rational, Refusal};
use catallax_shapes::Surface;
use catallax_topology::{BodyId, FaceId, Loop, Model, Sense};
use spade::{
    ConstrainedDelaunayTriangulation, HasPosition, HierarchyHintGenerator, Point2,
    PositionInTriangulation, Triangulation,
};

use crate::{
    FaceMesh, Mesh, MeshOptions, Triangle,
    chart::{self, Chart, Uv, align, chord_bound, key},
    edges::{self, Boundaries},
    failure,
};

#[derive(Clone, Copy)]
struct Node {
    uv: Point2<f64>,
    mesh: usize,
    error: usize,
}
impl HasPosition for Node {
    type Scalar = f64;
    fn position(&self) -> Point2<f64> {
        self.uv
    }
}
type Cdt = ConstrainedDelaunayTriangulation<Node, (), (), (), HierarchyHintGenerator<f64>>;

struct FaceBuilder<'a> {
    chart: &'a Chart,
    options: &'a MeshOptions,
    positions: &'a mut Vec<[f64; 3]>,
    rails: &'a edges::Rails,
    cdt: Cdt,
    canonical: BTreeMap<[u64; 2], (usize, usize)>,
    errors: Vec<Rational>,
    filtered_errors: Vec<Interval>,
    interior_left: BTreeSet<(usize, usize)>,
    prior_chords: &'a BTreeSet<(usize, usize)>,
    bounds: [[f64; 2]; 2],
    hessian: [Rational; 3],
    filtered_hessian: [Interval; 3],
    dec: &'a Decider,
}

impl FaceBuilder<'_> {
    fn node(&mut self, uv: Uv, boundary: Option<usize>) -> Result<Node, Refusal> {
        if self.positions.len() >= self.options.max_vertices {
            return Err(failure("face triangulation reached max_vertices"));
        }
        let canonical = self.chart.canonical(uv);
        let (mesh, error) = if let Some(&(mesh, error)) = self.canonical.get(&key(canonical)) {
            if boundary.is_some_and(|index| index != mesh) {
                return Err(failure(
                    "distinct boundary vertices collapse in the display chart",
                ));
            }
            (mesh, error)
        } else {
            let (p, residual) = self.chart.sample(
                canonical,
                boundary.map(|index| self.positions[index]),
                &self.options.chord_tolerance,
                self.dec,
            )?;
            let mesh = if let Some(index) = boundary {
                index
            } else {
                let index = self.positions.len();
                self.positions.push(p);
                index
            };
            let error = self.errors.len();
            self.filtered_errors
                .push(Interval::from_rational(&residual));
            self.errors.push(residual);
            self.canonical.insert(key(canonical), (mesh, error));
            (mesh, error)
        };
        Ok(Node {
            uv: Point2::new(uv[0], uv[1]),
            mesh,
            error,
        })
    }

    fn insert(
        &mut self,
        uv: Uv,
        boundary: Option<usize>,
    ) -> Result<spade::handles::FixedVertexHandle, Refusal> {
        if self.cdt.num_vertices() >= self.options.max_vertices {
            return Err(failure("face triangulation reached max_vertices"));
        }
        let node = self.node(uv, boundary)?;
        self.cdt
            .insert(node)
            .map_err(|e| failure(format!("chart vertex insertion: {e}")))
    }

    fn queue(
        &mut self,
        nodes: &mut Vec<Node>,
        indices: &mut BTreeMap<[u64; 2], usize>,
        uv: Uv,
        boundary: Option<usize>,
    ) -> Result<usize, Refusal> {
        if let Some(&index) = indices.get(&key(uv)) {
            if boundary.is_some_and(|mesh| mesh != nodes[index].mesh) {
                return Err(failure(
                    "distinct boundary vertices collapse in the display chart",
                ));
            }
            return Ok(index);
        }
        if nodes.len() >= self.options.max_vertices {
            return Err(failure("face triangulation reached max_vertices"));
        }
        let node = self.node(uv, boundary)?;
        let index = nodes.len();
        nodes.push(node);
        indices.insert(key(uv), index);
        Ok(index)
    }

    fn bulk_load(&mut self, nodes: Vec<Node>) -> Result<(), Refusal> {
        let constraints = self.interior_left.iter().map(|&(a, b)| [a, b]).collect();
        let mut conflict = false;
        // Input UVs are unique. Spade's stable bulk loader therefore preserves
        // the indices used by the oriented boundary constraints.
        let cdt = Cdt::try_bulk_load_cdt(nodes, constraints, |_| conflict = true)
            .map_err(|error| failure(format!("chart bulk insertion: {error}")))?;
        if conflict {
            return Err(failure("boundary polylines cross at this display accuracy"));
        }
        self.cdt = cdt;
        Ok(())
    }

    fn shifts(&self) -> impl Iterator<Item = Uv> + use<> {
        let u = if self.chart.periodic[0] {
            -1..=1
        } else {
            0..=0
        };
        let v = if self.chart.periodic[1] {
            -1..=1
        } else {
            0..=0
        };
        u.flat_map(move |i| v.clone().map(move |j| [f64::from(i), f64::from(j)]))
    }

    fn boundary(
        &mut self,
        segments: &[(Uv, Uv, usize, usize)],
        sense: Sense,
    ) -> Result<(), Refusal> {
        let mut nodes = Vec::new();
        let mut indices = BTreeMap::new();
        for &(mut a, mut b, first, last) in segments {
            for i in 0..2 {
                if self.chart.periodic[i] {
                    let shift = libm::floor(a[i]);
                    a[i] -= shift;
                    b[i] -= shift;
                }
            }
            for shift in self.shifts() {
                let a = self.queue(
                    &mut nodes,
                    &mut indices,
                    [a[0] + shift[0], a[1] + shift[1]],
                    Some(first),
                )?;
                let b = self.queue(
                    &mut nodes,
                    &mut indices,
                    [b[0] + shift[0], b[1] + shift[1]],
                    Some(last),
                )?;
                if a == b {
                    return Err(failure("boundary segment collapsed in parameter space"));
                }
                let pair = if sense == Sense::Same { (a, b) } else { (b, a) };
                self.interior_left.insert(pair);
            }
        }
        self.bulk_load(nodes)
    }

    fn free_point(&mut self, uv: Uv) -> Result<bool, Refusal> {
        let p = Point2::new(uv[0], uv[1]);
        match self.cdt.locate(p) {
            PositionInTriangulation::OnVertex(_) => return Ok(false),
            PositionInTriangulation::OnEdge(e)
                if self.cdt.is_constraint_edge(e.as_undirected()) =>
            {
                return Ok(false);
            }
            _ => {}
        }
        self.insert(uv, None)?;
        Ok(true)
    }

    fn grid(&mut self) -> Result<(), Refusal> {
        let widths = self.bounds.map(|b| b[1] - b[0]);
        let bounded_torus =
            matches!(self.chart.surface, Surface::Torus(_)) && !self.chart.periodic[1];
        let angular_step = self
            .options
            .angular_deviation
            .as_ref()
            .map(|limit| {
                catallax_numbers::Ball::pi_checked(chart::BITS, self.dec)
                    .map(|pi| limit / &pi.hi().mul_pow2(1))
            })
            .transpose()?;
        let rail_grid = matches!(self.chart.surface, Surface::Cone(_))
            || (bounded_torus && angular_step.is_some());
        let mut rail_u: Vec<_> = if rail_grid {
            self.cdt
                .vertices()
                .map(|vertex| vertex.data().uv.x)
                .collect()
        } else {
            Vec::new()
        };
        rail_u.sort_unstable_by(f64::total_cmp);
        // Keep distinct shared samples, but not nearly coincident interior
        // columns caused by rounding one longitude to adjacent chart quanta.
        rail_u.dedup_by(|a, b| (*a - *b).abs() <= 16.0 / 281_474_976_710_656.0);
        let rail_du = rail_u
            .windows(2)
            .map(|pair| (pair[1] - pair[0]).next_up())
            .fold(0.0_f64, f64::max);
        // Shared rails already resolve longitude more finely than a free grid.
        // Use their actual maximum gap, then refine only the other parameter.
        // The strict margin guarantees that the fixed longitude can meet both
        // controls. Sparse/partial trims retain the general two-dimensional grid.
        let align_rails = rail_du > 0.0
            && chord_bound(&[[0.0, 0.0], [rail_du, 0.0]], &self.hessian)
                <= self.options.chord_tolerance.mul_pow2(-2)
            && angular_step
                .as_ref()
                .is_none_or(|step| chart::q(rail_du) <= step.mul_pow2(-1));
        if !align_rails {
            rail_u.clear();
        }
        let mut counts = [1usize; 2];
        for (i, count) in counts.iter_mut().enumerate() {
            if self.chart.periodic[i] {
                *count = 8;
            }
        }
        // Reserve a quarter of the chord budget for endpoint replacement.
        // finish still certifies every facet at the full requested limit.
        let chord_limit = if align_rails {
            &self.options.chord_tolerance * &Rational::from_ratio(3, 4).unwrap()
        } else {
            self.options.chord_tolerance.mul_pow2(-1)
        };
        loop {
            let du = if align_rails {
                rail_du
            } else {
                widths[0] / counts[0] as f64
            };
            let dv = widths[1] / counts[1] as f64;
            let chord_reached = chord_bound(&[[0.0, 0.0], [du, dv]], &self.hessian) <= chord_limit;
            let (angular_reached, angular_dimension) = if let Some(step) = &angular_step {
                let u = match self.chart.surface {
                    Surface::Plane(_) => Rational::zero(),
                    Surface::EllipticCylinder(_) => &chart::q(du) * &self.chart.normal_scale,
                    _ => chart::q(du),
                };
                let v = if matches!(self.chart.surface, Surface::Sphere(_) | Surface::Torus(_)) {
                    chart::q(dv)
                } else {
                    Rational::zero()
                };
                (&u + &v <= *step, usize::from(v > u))
            } else {
                (true, 0)
            };
            if chord_reached && angular_reached {
                break;
            }
            let dimension = if chord_reached {
                // A thin surface can meet its chord bound before its normals
                // are resolved. Seed both angular directions rather than
                // repeatedly inserting centroids into a chord-only grid.
                angular_dimension
            } else {
                let u = &self.hessian[0] * &chart::q(du) * chart::q(du);
                let v = &self.hessian[2] * &chart::q(dv) * chart::q(dv);
                usize::from(v > u)
            };
            // A cone's mixed derivative also needs axial subdivisions.
            let dimension = if align_rails {
                1
            } else if !bounded_torus
                && !self.hessian[1].is_zero()
                && counts[dimension] > counts[1 - dimension] * 2
            {
                1 - dimension
            } else {
                dimension
            };
            counts[dimension] = counts[dimension]
                .checked_mul(2)
                .ok_or_else(|| failure("grid size overflow"))?;
            let columns = if align_rails { rail_u.len() } else { counts[0] };
            if columns.saturating_mul(counts[1]) > self.options.max_vertices {
                return Err(failure("tolerance exceeds the deterministic grid budget"));
            }
        }
        let ranges: [(i64, i64); 2] = std::array::from_fn(|i| {
            if self.chart.periodic[i] {
                (-(counts[i] as i64), 2 * counts[i] as i64)
            } else {
                (0, counts[i] as i64)
            }
        });
        let u_start = self.bounds[0][0];
        let parameters_u = (ranges[0].0..=ranges[0].1)
            .take(if align_rails { 0 } else { usize::MAX })
            .map(|i| u_start + widths[0] * i as f64 / counts[0] as f64)
            .chain(rail_u);
        let mut queued = None;
        for u in parameters_u {
            for j in ranges[1].0..=ranges[1].1 {
                // Align interior rows with shared rails. A staggered major-circle
                // chord can dominate a thin tube's minor-direction tangent, and
                // centroid refinement only makes that triangle's altitude smaller.
                // Boundary samples already own the end rows, including their
                // projection residuals; do not seed nearly coincident rail copies.
                let endpoint = j == 0 || j == counts[1] as i64;
                let apex = matches!(self.chart.surface, Surface::Cone(_))
                    && self.bounds[1][usize::from(j != 0)] == 0.0;
                if align_rails && endpoint && !apex {
                    continue;
                }
                let uv = [
                    u,
                    self.bounds[1][0] + widths[1] * j as f64 / counts[1] as f64,
                ];
                // Query the fixed boundary triangulation before bulk loading.
                // A new point on a rail would split only this face's constraint
                // and destroy the shared polyline on its adjacent face.
                match self.cdt.locate(Point2::new(uv[0], uv[1])) {
                    PositionInTriangulation::OnVertex(_) => continue,
                    PositionInTriangulation::OnEdge(edge)
                        if self.cdt.is_constraint_edge(edge.as_undirected()) =>
                    {
                        continue;
                    }
                    _ => {}
                }
                let (nodes, indices) = queued.get_or_insert_with(|| {
                    let nodes: Vec<_> = self.cdt.vertices().map(|vertex| *vertex.data()).collect();
                    let indices = nodes
                        .iter()
                        .enumerate()
                        .map(|(i, node)| (key([node.uv.x, node.uv.y]), i))
                        .collect();
                    (nodes, indices)
                });
                self.queue(nodes, indices, uv, None)?;
            }
        }
        if let Some((nodes, indices)) = queued {
            drop(indices);
            self.bulk_load(nodes)
        } else {
            Ok(())
        }
    }

    fn labels(&self) -> Result<Vec<bool>, Refusal> {
        // Classify the periodic quotient, not the finite copies' outer hull.
        // Incomplete far-away loop copies must not connect interior and exterior.
        let triangle_key = |points: [Point2<f64>; 3]| {
            // Preserve the relative lifts of all three vertices. Reducing each
            // vertex modulo one aliases triangles on opposite sides of a seam.
            let center = [
                (points[0].x + points[1].x + points[2].x) / 3.0,
                (points[0].y + points[1].y + points[2].y) / 3.0,
            ];
            let mut keys = points.map(|p| {
                let mut uv = [p.x, p.y];
                for i in 0..2 {
                    if self.chart.periodic[i] {
                        uv[i] -= libm::floor(center[i]);
                    }
                    if uv[i] == 0.0 {
                        uv[i] = 0.0;
                    }
                }
                key(uv)
            });
            keys.sort_unstable();
            keys
        };
        let mut central = BTreeMap::new();
        for f in self.cdt.inner_faces() {
            let points = f.positions();
            let center = [
                (points[0].x + points[1].x + points[2].x) / 3.0,
                (points[0].y + points[1].y + points[2].y) / 3.0,
            ];
            if self.central(center) {
                central.insert(triangle_key(points), f.fix());
            }
        }
        let mut labels = vec![None; self.cdt.num_all_faces()];
        for &start in central.values() {
            if labels[start.index()].is_some() {
                continue;
            }
            let mut group = vec![];
            let mut queue = VecDeque::from([start]);
            let mut seen = BTreeSet::from([start.index()]);
            let mut seed = None;
            while let Some(f) = queue.pop_front() {
                group.push(f.index());
                for e in self.cdt.face(f).adjacent_edges() {
                    if e.is_constraint_edge() {
                        let pair = (e.from().index(), e.to().index());
                        let left = self.interior_left.contains(&pair);
                        let right = self.interior_left.contains(&(pair.1, pair.0));
                        if left == right {
                            return Err(failure("a boundary constraint was split or overlapped"));
                        }
                        if seed.is_some_and(|s| s != left) {
                            return Err(failure(
                                "inconsistent oriented boundary regions at this accuracy",
                            ));
                        }
                        seed = Some(left);
                    } else if let Some(other) = e.rev().face().as_inner()
                        && let Some(&other) = central.get(&triangle_key(other.positions()))
                        && seen.insert(other.index())
                    {
                        queue.push_back(other);
                    }
                }
            }
            let inside = seed.unwrap_or(self.interior_left.is_empty());
            for index in group {
                labels[index] = Some(inside);
            }
        }
        Ok(labels.into_iter().map(|v| v.unwrap_or(false)).collect())
    }

    fn central(&self, uv: Uv) -> bool {
        (0..2).all(|i| !self.chart.periodic[i] || (0.0..1.0).contains(&uv[i]))
    }

    fn ear_point(&self, uvs: &[Uv; 3], center: Uv) -> Result<Uv, Refusal> {
        let mut point = self.chart.canonical(center);
        for (value, periodic) in point.iter_mut().zip(self.chart.periodic) {
            if periodic {
                *value = libm::round(*value * 281_474_976_710_656.0) / 281_474_976_710_656.0;
            }
        }
        let point = self.chart.canonical(point);
        let lifted = align(point, center, self.chart.periodic);
        // Strict orientation certifies an interior chart point, rather than a
        // rounded point on the shared polyline. Uncertain filters remain exact.
        let inside = (0..3).all(|i| chart::positive_turn(uvs[i], uvs[(i + 1) % 3], lifted));
        if !inside {
            return Err(failure(
                "shared-rail ear has no interior point at chart precision",
            ));
        }
        Ok(point)
    }

    fn finish(
        &mut self,
        face: FaceId,
        sense: Sense,
        dec: &Decider,
    ) -> Result<(Vec<Triangle>, Rational, Rational, Option<Rational>), Refusal> {
        let filtered_tolerance = Interval::from_rational(&self.options.chord_tolerance);
        for _ in 0..48 {
            dec.check_deadline()?;
            let labels = self.labels()?;
            let mut pending = BTreeMap::new();
            let mut triangles = Vec::new();
            let mut maximum = Rational::zero();
            let mut aspect = Rational::zero();
            for f in self.cdt.inner_faces() {
                if !labels[f.index()] {
                    continue;
                }
                let nodes = f.vertices().map(|v| *v.data());
                let uvs = nodes.map(|v| [v.uv.x, v.uv.y]);
                let center = std::array::from_fn(|i| (uvs[0][i] + uvs[1][i] + uvs[2][i]) / 3.0);
                if !self.central(center) {
                    continue;
                }
                let mut vertices = nodes.map(|n| n.mesh);
                if vertices[0] == vertices[1]
                    || vertices[1] == vertices[2]
                    || vertices[2] == vertices[0]
                {
                    continue;
                }
                if self.rails.contains_triangle(vertices) {
                    let point = self.ear_point(&uvs, center)?;
                    pending.insert(key(point), point);
                    continue;
                }
                // A face-local diagonal must not become a second shared edge.
                // Two incident faces can otherwise emit the same boundary-only
                // triangle in opposite directions, producing four edge uses.
                let mut split = false;
                for edge in f.adjacent_edges() {
                    let a = *edge.from().data();
                    let b = *edge.to().data();
                    if !edge.is_constraint_edge()
                        && self
                            .prior_chords
                            .contains(&(a.mesh.min(b.mesh), a.mesh.max(b.mesh)))
                    {
                        let midpoint = self
                            .chart
                            .canonical([(a.uv.x + b.uv.x) * 0.5, (a.uv.y + b.uv.y) * 0.5]);
                        pending.insert(key(midpoint), midpoint);
                        split = true;
                    }
                }
                if split {
                    continue;
                }
                let residual = nodes.iter().map(|n| &self.errors[n.error]).max().unwrap();
                let filtered_residual = nodes
                    .iter()
                    .map(|n| self.filtered_errors[n.error].hi())
                    .fold(0.0_f64, f64::max);
                let vmax = uvs[0][1].abs().max(uvs[1][1].abs()).max(uvs[2][1].abs());
                let filtered = if matches!(self.chart.surface, Surface::Cone(_)) {
                    let scale = Interval::point(vmax).unwrap()
                        / Interval::point(self.bounds[1][0].abs().max(self.bounds[1][1].abs()))
                            .unwrap();
                    let h = [
                        self.filtered_hessian[0] * scale,
                        self.filtered_hessian[1],
                        self.filtered_hessian[2],
                    ];
                    chart::filtered_chord_bound(&uvs, &h)
                        + Interval::point(filtered_residual).unwrap()
                } else {
                    chart::filtered_chord_bound(&uvs, &self.filtered_hessian)
                        + Interval::point(filtered_residual).unwrap()
                };
                let error = if filtered.is_finite() && filtered.hi() <= filtered_tolerance.lo() {
                    chart::q(filtered.hi())
                } else {
                    let hessian = self.chart.hessian(vmax, dec)?;
                    &chord_bound(&uvs, &hessian) + residual
                };
                if let Some(limit) = &self.options.max_edge_length {
                    let lengths = std::array::from_fn::<_, 3, _>(|i| {
                        crate::controls::length_squared(
                            self.positions[vertices[i]],
                            self.positions[vertices[(i + 1) % 3]],
                        )
                    });
                    let longest = (0..3).max_by_key(|&i| &lengths[i]).unwrap();
                    if lengths[longest] > limit * limit {
                        let next = (longest + 1) % 3;
                        let midpoint = self.chart.canonical(std::array::from_fn(|i| {
                            (uvs[longest][i] + uvs[next][i]) * 0.5
                        }));
                        pending.insert(key(midpoint), midpoint);
                        continue;
                    }
                }
                if error > self.options.chord_tolerance {
                    let mut center = self.chart.canonical(center);
                    for (value, periodic) in center.iter_mut().zip(self.chart.periodic) {
                        if periodic {
                            *value =
                                libm::round(*value * 281_474_976_710_656.0) / 281_474_976_710_656.0;
                        }
                    }
                    pending.insert(key(center), center);
                    continue;
                }
                let (ratio, angular_reached) = crate::controls::quality(
                    self.chart,
                    &uvs,
                    vertices.map(|i| self.positions[i]),
                    self.options.angular_deviation.as_ref(),
                    dec,
                )?;
                let ratio = ratio.ok_or_else(|| failure("display triangle has zero area"))?;
                if !angular_reached {
                    let mut center = self.chart.canonical(center);
                    for (value, periodic) in center.iter_mut().zip(self.chart.periodic) {
                        if periodic {
                            *value =
                                libm::round(*value * 281_474_976_710_656.0) / 281_474_976_710_656.0;
                        }
                    }
                    pending.insert(key(center), center);
                    continue;
                }
                aspect = aspect.max(ratio);
                maximum = maximum.max(error);
                if sense == Sense::Reversed {
                    vertices.swap(1, 2);
                }
                let min = (0..3).min_by_key(|&i| vertices[i]).unwrap();
                vertices.rotate_left(min);
                triangles.push(Triangle { vertices, face });
            }
            if pending.is_empty() {
                triangles.sort_unstable_by_key(|t| t.vertices);
                return Ok((
                    triangles,
                    maximum,
                    aspect,
                    self.options.angular_deviation.clone(),
                ));
            }
            let mut added = false;
            for center in pending.into_values() {
                for shift in self.shifts() {
                    added |= self.free_point([center[0] + shift[0], center[1] + shift[1]])?;
                }
            }
            if !added {
                return Err(failure(
                    "coordinate precision cannot meet the requested chord tolerance",
                ));
            }
        }
        Err(failure(
            "face refinement reached its deterministic subdivision limit",
        ))
    }
}

fn face_mesh(
    model: &Model,
    face_id: FaceId,
    chart: &mut Chart,
    boundaries: &mut Boundaries,
    options: &MeshOptions,
    prior_chords: &BTreeSet<(usize, usize)>,
    dec: &Decider,
) -> Result<(), Refusal> {
    let face = model.face(face_id);
    let mut segments = Vec::new();
    let mut projected = Vec::new();
    for &l in &face.loops {
        match model.loop_(l) {
            Loop::Vertex(v) => {
                projected.push(chart.project(boundaries.mesh.positions[boundaries.vertices[v]]))
            }
            Loop::Coedges(cs) => {
                for &c in cs {
                    let c = model.coedge(c);
                    let edge = &boundaries.mesh.edges[boundaries.edges[&c.edge]];
                    for pair in edge.vertices.windows(2) {
                        let pair = if c.sense == Sense::Same {
                            [pair[0], pair[1]]
                        } else {
                            [pair[1], pair[0]]
                        };
                        let a = chart.project(boundaries.mesh.positions[pair[0]]);
                        let b = align(
                            chart.project(boundaries.mesh.positions[pair[1]]),
                            a,
                            chart.periodic,
                        );
                        projected.extend([a, b]);
                        segments.push((a, b, pair[0], pair[1]));
                    }
                }
            }
        }
    }
    let bounds = std::array::from_fn(|i| {
        if chart.periodic[i] {
            [0.0, 1.0]
        } else if i == 1 && matches!(chart.surface, Surface::Sphere(_)) {
            [-0.25, 0.25]
        } else {
            [
                projected.iter().map(|p| p[i]).fold(f64::INFINITY, f64::min),
                projected
                    .iter()
                    .map(|p| p[i])
                    .fold(f64::NEG_INFINITY, f64::max),
            ]
        }
    });
    if bounds
        .iter()
        .any(|b| !b[0].is_finite() || !b[1].is_finite() || b[0] >= b[1])
    {
        return Err(failure(
            "face has no finite two-dimensional parameter region",
        ));
    }
    let hessian = chart.hessian(bounds[1][0].abs().max(bounds[1][1].abs()), dec)?;
    let mut builder = FaceBuilder {
        chart,
        options,
        positions: &mut boundaries.mesh.positions,
        rails: &boundaries.rails,
        cdt: Cdt::new(),
        canonical: BTreeMap::new(),
        errors: Vec::new(),
        filtered_errors: Vec::new(),
        interior_left: BTreeSet::new(),
        prior_chords,
        bounds,
        filtered_hessian: hessian.each_ref().map(Interval::from_rational),
        hessian: hessian.into_owned(),
        dec,
    };
    builder.boundary(&segments, face.sense)?;
    builder.grid()?;
    let (triangles, max_chord_error, max_aspect_ratio, max_angular_deviation) = builder
        .finish(face_id, face.sense, dec)
        .map_err(|e| Refusal::new(e.code, format!("face {}: {}", face_id.index(), e.detail)))?;
    let start = boundaries.mesh.triangles.len();
    boundaries.mesh.triangles.extend(triangles);
    if options
        .max_aspect_ratio
        .as_ref()
        .is_some_and(|limit| max_aspect_ratio > *limit)
    {
        boundaries
            .mesh
            .unavailable_controls
            .push(crate::UnavailableControl {
                face: face_id,
                control: crate::MeshControl::AspectRatio,
                reason: "the certified aspect-ratio bound exceeds the requested limit",
            });
    }
    boundaries.mesh.faces.push(FaceMesh {
        face: face_id,
        triangles: start..boundaries.mesh.triangles.len(),
        max_chord_error,
        max_aspect_ratio: Some(max_aspect_ratio),
        max_angular_deviation,
    });
    Ok(())
}

fn compact(mesh: &mut Mesh) {
    let mut used = vec![false; mesh.positions.len()];
    for t in &mesh.triangles {
        for &v in &t.vertices {
            used[v] = true;
        }
    }
    for e in &mesh.edges {
        for &v in &e.vertices {
            used[v] = true;
        }
    }
    let mut remap = vec![0; used.len()];
    let mut next = 0;
    mesh.positions.retain(|_| {
        let keep = used[next];
        next += 1;
        keep
    });
    let mut count = 0;
    for (i, used) in used.into_iter().enumerate() {
        if used {
            remap[i] = count;
            count += 1;
        }
    }
    for t in &mut mesh.triangles {
        for v in &mut t.vertices {
            *v = remap[*v];
        }
    }
    for e in &mut mesh.edges {
        for v in &mut e.vertices {
            *v = remap[*v];
        }
    }
}

pub(crate) fn check_closed(mesh: &Mesh) -> Result<(), Refusal> {
    let mut uses = BTreeMap::<(usize, usize), (usize, i32)>::new();
    for t in &mesh.triangles {
        for i in 0..3 {
            let (a, b) = (t.vertices[i], t.vertices[(i + 1) % 3]);
            let entry = uses.entry((a.min(b), a.max(b))).or_default();
            entry.0 += 1;
            entry.1 += if a < b { 1 } else { -1 };
        }
    }
    if uses.values().any(|&(n, winding)| n != 2 || winding != 0) {
        return Err(failure(
            "periodic or trimmed triangulation did not close consistently",
        ));
    }
    Ok(())
}

fn regular_chart(
    model: &Model,
    face: FaceId,
    options: &MeshOptions,
    dec: &Decider,
) -> Result<Chart, Refusal> {
    let surface = model.geom().surface(model.face(face).surface);
    let Surface::Sphere(sphere) = surface else {
        return Chart::new(surface, dec);
    };
    let equations: Vec<_> = model
        .face(face)
        .loops
        .iter()
        .flat_map(|l| match model.loop_(*l) {
            Loop::Coedges(cs) => cs
                .iter()
                .map(|c| model.edge(model.coedge(*c).edge).curve)
                .collect(),
            Loop::Vertex(_) => Vec::new(),
        })
        .map(|c| {
            model
                .geom()
                .curve(c)
                .equations(model.geom())
                .ok_or_else(|| {
                    failure("sphere pole separation requires a finite parametric boundary solver")
                })
        })
        .collect::<Result<Vec<_>, _>>()?;
    if equations.is_empty() {
        return Chart::new(surface, dec);
    }
    for choice in 0..options
        .max_vertices
        .min(equations.len().saturating_mul(64).saturating_add(1))
    {
        dec.check_deadline()?;
        let k = Rational::from_i64(choice as i64);
        let rotation = catallax_shapes::Rotation::from_quaternion(
            Rational::one(),
            k.clone(),
            &k * &k,
            k.pow(4).unwrap(),
        )?;
        let mut candidate = sphere.clone();
        candidate.axis = rotation.apply(&sphere.axis);
        candidate.ref_dir = rotation.apply(&sphere.ref_dir);
        let chart = Chart::new(&Surface::Sphere(candidate), dec)?;
        let poles = [
            chart.evaluate([0.0, -0.25], dec)?,
            chart.evaluate([0.0, 0.25], dec)?,
        ];
        let mut clear = true;
        'poles: for pole in &poles {
            for equations in &equations {
                let mut apart = false;
                for equation in equations {
                    if !equation.eval_ball(pole, dec)?.contains_zero() {
                        apart = true;
                        break;
                    }
                }
                if !apart {
                    clear = false;
                    break 'poles;
                }
            }
        }
        if clear {
            return Ok(chart);
        }
    }
    Err(failure(
        "no sphere chart poles separated from the boundary carriers within the vertex budget",
    ))
}

pub(crate) fn build(
    model: &Model,
    body: BodyId,
    options: &MeshOptions,
    dec: &mut Decider,
) -> Answered<Mesh> {
    let faces = model.faces_of_body(body);
    if faces.iter().any(|face| {
        matches!(
            model.geom().surface(model.face(*face).surface),
            Surface::BSpline(_)
        )
    }) {
        return crate::spline_prism::build(model, body, &faces, options, dec);
    }
    let mut charts = faces
        .iter()
        .map(|&f| Ok((f, regular_chart(model, f, options, dec)?)))
        .collect::<Result<BTreeMap<_, _>, Refusal>>()?;
    let mut boundaries = edges::sample(model, body, &charts, options, dec)?;
    let boundary_limit = boundaries.mesh.positions.len();
    let mut prior_chords = BTreeSet::new();
    for face in faces {
        let start = boundaries.mesh.triangles.len();
        if let Surface::Plane(plane) = &charts[&face].surface {
            crate::planar::build(model, face, plane, &mut boundaries, options, dec).map_err(
                |e| Refusal::new(e.code, format!("face {}: {}", face.index(), e.detail)),
            )?;
        } else {
            if let Some(center) =
                crate::bounded::minor_angle_center(model, face, &charts[&face], &boundaries, dec)?
            {
                charts.get_mut(&face).unwrap().restrict(1, center);
            }
            let strip = matches!(
                charts[&face].surface,
                Surface::Cylinder(_)
                    | Surface::EllipticCylinder(_)
                    | Surface::Cone(_)
                    | Surface::Torus(_)
            ) && crate::ruled::build(
                model,
                face,
                &charts[&face],
                &mut boundaries,
                options,
                dec,
            )?;
            if !strip {
                face_mesh(
                    model,
                    face,
                    charts.get_mut(&face).unwrap(),
                    &mut boundaries,
                    options,
                    &prior_chords,
                    dec,
                )?;
            }
        }
        let mut once = BTreeSet::new();
        for triangle in &boundaries.mesh.triangles[start..] {
            for i in 0..3 {
                let (a, b) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
                if a < boundary_limit && b < boundary_limit {
                    let pair = (a.min(b), a.max(b));
                    if !once.insert(pair) {
                        prior_chords.insert(pair);
                    }
                }
            }
        }
    }
    check_closed(&boundaries.mesh)?;
    compact(&mut boundaries.mesh);
    Ok(Answer {
        value: boundaries.mesh,
        marks: boundaries.marks,
    })
}
