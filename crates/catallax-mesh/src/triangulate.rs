use std::collections::{BTreeMap, BTreeSet, VecDeque};

use catallax_numbers::{Answer, Answered, Decider, Rational, Refusal};
use catallax_shapes::Surface;
use catallax_topology::{BodyId, FaceId, Loop, Model, Sense};
use spade::{
    ConstrainedDelaunayTriangulation, HasPosition, Point2, PositionInTriangulation, Triangulation,
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
type Cdt = ConstrainedDelaunayTriangulation<Node>;

struct FaceBuilder<'a> {
    chart: &'a Chart,
    options: &'a MeshOptions,
    positions: &'a mut Vec<[f64; 3]>,
    cdt: Cdt,
    canonical: BTreeMap<[u64; 2], (usize, usize)>,
    errors: Vec<Rational>,
    interior_left: BTreeSet<(usize, usize)>,
    bounds: [[f64; 2]; 2],
    hessian: [Rational; 3],
    dec: &'a Decider,
}

impl FaceBuilder<'_> {
    fn insert(
        &mut self,
        uv: Uv,
        boundary: Option<usize>,
    ) -> Result<spade::handles::FixedVertexHandle, Refusal> {
        if self.cdt.num_vertices() >= self.options.max_vertices
            || self.positions.len() >= self.options.max_vertices
        {
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
            let ball = self.chart.evaluate(canonical, self.dec)?;
            let mesh = if let Some(index) = boundary {
                index
            } else {
                let p = chart::position(&ball)?;
                let index = self.positions.len();
                self.positions.push(p);
                index
            };
            let error = self.errors.len();
            self.errors
                .push(chart::residual(&self.positions[mesh], &ball));
            self.canonical.insert(key(canonical), (mesh, error));
            (mesh, error)
        };
        self.cdt
            .insert(Node {
                uv: Point2::new(uv[0], uv[1]),
                mesh,
                error,
            })
            .map_err(|e| failure(format!("chart vertex insertion: {e}")))
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
        for &(mut a, mut b, first, last) in segments {
            for i in 0..2 {
                if self.chart.periodic[i] {
                    let shift = libm::floor(a[i]);
                    a[i] -= shift;
                    b[i] -= shift;
                }
            }
            for shift in self.shifts() {
                let a = self.insert([a[0] + shift[0], a[1] + shift[1]], Some(first))?;
                let b = self.insert([b[0] + shift[0], b[1] + shift[1]], Some(last))?;
                if a == b {
                    return Err(failure("boundary segment collapsed in parameter space"));
                }
                if !self.cdt.can_add_constraint(a, b) {
                    return Err(failure("boundary polylines cross at this display accuracy"));
                }
                self.cdt.add_constraint(a, b);
                let pair = if sense == Sense::Same {
                    (a.index(), b.index())
                } else {
                    (b.index(), a.index())
                };
                self.interior_left.insert(pair);
            }
        }
        Ok(())
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
        let mut counts = [1usize; 2];
        for (i, count) in counts.iter_mut().enumerate() {
            if self.chart.periodic[i] {
                *count = 8;
            }
        }
        loop {
            let du = widths[0] / counts[0] as f64;
            let dv = widths[1] / counts[1] as f64;
            if chord_bound(&[[0.0, 0.0], [du, dv]], &self.hessian)
                <= self.options.chord_tolerance.mul_pow2(-1)
            {
                break;
            }
            let u = &self.hessian[0] * &chart::q(du) * chart::q(du);
            let v = &self.hessian[2] * &chart::q(dv) * chart::q(dv);
            let dimension = if u >= v { 0 } else { 1 };
            // A cone's mixed derivative also needs axial subdivisions.
            let dimension =
                if !self.hessian[1].is_zero() && counts[dimension] > counts[1 - dimension] * 2 {
                    1 - dimension
                } else {
                    dimension
                };
            counts[dimension] = counts[dimension]
                .checked_mul(2)
                .ok_or_else(|| failure("grid size overflow"))?;
            if counts[0].saturating_mul(counts[1]) > self.options.max_vertices {
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
        for i in ranges[0].0..=ranges[0].1 {
            for j in ranges[1].0..=ranges[1].1 {
                let uv = std::array::from_fn(|k| {
                    self.bounds[k][0] + widths[k] * [i, j][k] as f64 / counts[k] as f64
                });
                self.free_point(uv)?;
            }
        }
        Ok(())
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

    fn finish(
        &mut self,
        face: FaceId,
        sense: Sense,
        dec: &Decider,
    ) -> Result<(Vec<Triangle>, Rational, Rational, Option<Rational>), Refusal> {
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
                let residual = nodes.iter().map(|n| &self.errors[n.error]).max().unwrap();
                let error = &chord_bound(&uvs, &self.hessian) + residual;
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
    chart: &Chart,
    boundaries: &mut Boundaries,
    options: &MeshOptions,
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
    let mut builder = FaceBuilder {
        chart,
        options,
        positions: &mut boundaries.mesh.positions,
        cdt: Cdt::new(),
        canonical: BTreeMap::new(),
        errors: Vec::new(),
        interior_left: BTreeSet::new(),
        bounds,
        hessian: chart.hessian(bounds[1][0].abs().max(bounds[1][1].abs()), dec)?,
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
        max_aspect_ratio,
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

fn check_closed(mesh: &Mesh) -> Result<(), Refusal> {
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
        .map(|c| model.geom().curve(c).equations(model.geom()))
        .collect();
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
    let charts = faces
        .iter()
        .map(|&f| Ok((f, regular_chart(model, f, options, dec)?)))
        .collect::<Result<BTreeMap<_, _>, Refusal>>()?;
    let mut boundaries = edges::sample(model, body, &charts, options, dec)?;
    for face in faces {
        if let Surface::Plane(plane) = &charts[&face].surface {
            crate::planar::build(model, face, plane, &mut boundaries, options, dec).map_err(
                |e| Refusal::new(e.code, format!("face {}: {}", face.index(), e.detail)),
            )?;
        } else {
            face_mesh(model, face, &charts[&face], &mut boundaries, options, dec)?;
        }
    }
    check_closed(&boundaries.mesh)?;
    compact(&mut boundaries.mesh);
    Ok(Answer {
        value: boundaries.mesh,
        marks: boundaries.marks,
    })
}
