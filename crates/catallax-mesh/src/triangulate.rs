use std::collections::{BTreeMap, BTreeSet, VecDeque};

use catallax_numbers::{
    Answer, Answered, Decider, Interval, Outcome, Rational, Real, RealPoly, Refusal, StrictSign,
};
use catallax_shapes::{Curve, Sphere, Surface, Vector3};
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
    unproved: bool,
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
        } else if self.unproved {
            let mesh = if let Some(index) = boundary {
                index
            } else {
                let point = self.chart.preview_position(canonical, self.dec)?;
                let index = self.positions.len();
                self.positions.push(point);
                index
            };
            self.canonical.insert(key(canonical), (mesh, 0));
            (mesh, 0)
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

    /// Grid points strictly inside a boundary segment's diametral disk. A seed
    /// there leaves a thin facet on a fixed chord, and centroid refinement only
    /// flattens that facet further from the surface normal. The floating test
    /// selects seeds and proves nothing: `finish` still certifies every facet.
    fn encroaching_seeds(
        &self,
        ranges: [(i64, i64); 2],
        counts: [usize; 2],
        parameter: impl Fn(usize, i64) -> f64,
    ) -> Result<BTreeSet<[u64; 2]>, Refusal> {
        let mut protected = BTreeSet::new();
        for edge in self.cdt.undirected_edges() {
            self.dec.check_deadline()?;
            if !edge.is_constraint_edge() {
                continue;
            }
            let [a, b] = edge.positions();
            let radius = ((a.x - b.x).powi(2) + (a.y - b.y).powi(2)).sqrt() * 0.5;
            let centers = [(a.x + b.x) * 0.5, (a.y + b.y) * 0.5];
            // One index of padding absorbs the rounding of this mapping.
            let spans: [(i64, i64); 2] = std::array::from_fn(|i| {
                let [start, end] = self.bounds[i];
                let scale = counts[i] as f64 / (end - start);
                let first = libm::floor((centers[i] - radius - start) * scale) as i64;
                let last = libm::ceil((centers[i] + radius - start) * scale) as i64;
                (
                    first.saturating_sub(1).max(ranges[i].0),
                    last.saturating_add(1).min(ranges[i].1),
                )
            });
            for i in spans[0].0..=spans[0].1 {
                let u = parameter(0, i);
                for j in spans[1].0..=spans[1].1 {
                    self.dec.check_deadline()?;
                    let v = parameter(1, j);
                    if (u - a.x) * (u - b.x) + (v - a.y) * (v - b.y) < 0.0 {
                        protected.insert(key([u, v]));
                    }
                }
            }
        }
        Ok(protected)
    }

    /// Whether a `du` by `dv` chart cell of a cylinder already meets the edge limit,
    /// or leaves any excess to a single split of its diagonal.
    ///
    /// A cylinder chart needs fine columns for its chord but no rows, and its turns
    /// and millimetres are not one metric. A long edge left to `finish` therefore
    /// splits off the columns, and Delaunay refinement then spreads toward square
    /// chart cells: near a million vertices on a strip 10 mm long. Seeding rows keeps
    /// the lattice instead. This estimate only sizes that seed grid; `finish` still
    /// measures every returned edge.
    fn cylinder_cell_fits(&self, du: f64, dv: f64) -> bool {
        let Some(limit) = self.edge_limit() else {
            return true;
        };
        let around = self.hessian[0].to_f64_nearest() * du / std::f64::consts::TAU;
        dv <= limit * 0.5 || libm::sqrt(around * around + dv * dv) <= limit
    }

    fn edge_limit(&self) -> Option<f64> {
        let limit = self.options.max_edge_length.as_ref()?;
        matches!(
            self.chart.surface,
            Surface::Cylinder(_) | Surface::EllipticCylinder(_)
        )
        .then(|| limit.to_f64_nearest())
    }

    /// The shared columns a cylinder's interior rows need, and their widest gap
    /// inside the trim. A rail sampled for a tighter neighbouring surface must not
    /// repeat its samples on every row, so one period keeps only the columns its
    /// own chord, angle and length ask for. Each kept column is a rail sample, and
    /// the periodic copies stay bit-identical shifts.
    fn row_columns(
        &self,
        columns: &[f64],
        widest_gap: f64,
        angular_step: Option<&Rational>,
    ) -> (Vec<f64>, f64) {
        let tolerance = self.options.chord_tolerance.to_f64_nearest();
        let curvature = self.hessian[0].to_f64_nearest();
        let chord = libm::sqrt(4.0 * tolerance / curvature);
        let angle = angular_step.map_or(f64::INFINITY, |step| {
            step.to_f64_nearest() * 0.5 / self.chart.normal_scale.to_f64_nearest()
        });
        let length = self.edge_limit().map_or(f64::INFINITY, |limit| {
            limit * 0.75 * std::f64::consts::TAU / curvature
        });
        let target = chord.min(angle).min(length);
        let period: Vec<_> = columns
            .iter()
            .copied()
            .filter(|u| (0.0..1.0).contains(u))
            .collect();
        let mut kept = Vec::<f64>::new();
        let mut du = 0.0_f64;
        for (i, &u) in period.iter().enumerate() {
            let next = period.get(i + 1).copied();
            // A gap wider than any boundary segment lies outside the trim.
            let after_empty = i > 0 && u - period[i - 1] > widest_gap;
            let before_empty = next.is_some_and(|next| next - u > widest_gap);
            let needed = next
                .zip(kept.last())
                .is_some_and(|(next, last)| next - last > target);
            if i == 0 || next.is_none() || after_empty || before_empty || needed {
                if let Some(last) = kept.last().filter(|_| !after_empty) {
                    du = du.max((u - last).next_up());
                }
                kept.push(u);
            }
        }
        let mut lifted: Vec<_> = self
            .shifts()
            .flat_map(|shift| kept.iter().map(move |u| u + shift[0]).collect::<Vec<_>>())
            .collect();
        lifted.sort_unstable_by(f64::total_cmp);
        (lifted, du)
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
        // No column gap inside a trim is wider than its widest boundary segment.
        let widest_segment = self.edge_limit().map(|_| {
            self.cdt
                .undirected_edges()
                .filter(|edge| edge.is_constraint_edge())
                .map(|edge| {
                    let [a, b] = edge.positions();
                    (a.x - b.x).abs()
                })
                .fold(0.0_f64, f64::max)
        });
        let edge_rows = widest_segment.is_some_and(|du| !self.cylinder_cell_fits(du, widths[1]));
        let rail_grid = matches!(self.chart.surface, Surface::Cone(_))
            || (bounded_torus && angular_step.is_some())
            || edge_rows;
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
        // A partial face's rows must not count the empty arc beyond its trim.
        let widest_gap = widest_segment
            .filter(|_| edge_rows)
            .unwrap_or(f64::INFINITY);
        let mut rail_du = rail_u
            .windows(2)
            .map(|pair| pair[1] - pair[0])
            .filter(|gap| *gap <= widest_gap)
            .map(f64::next_up)
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
        let edge_rows = edge_rows && align_rails;
        if edge_rows {
            (rail_u, rail_du) = self.row_columns(&rail_u, widest_gap, angular_step.as_ref());
        }
        let mut counts = [1usize; 2];
        for (i, count) in counts.iter_mut().enumerate() {
            if self.chart.periodic[i] {
                *count = 8;
            }
        }
        if self.unproved && !matches!(self.chart.surface, Surface::Plane(_)) {
            // Curved end trims need interior axial rows even on a linear ruling.
            counts[1] = counts[1].max(8);
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
            let edge_reached = !edge_rows || self.cylinder_cell_fits(du, dv);
            if chord_reached && angular_reached && edge_reached {
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
        // An untrimmed periodic grid has no constraints or excluded seeds.
        // Count its temporary copies before allocating them, not while queuing.
        if self.cdt.num_vertices() == 0
            && self.chart.corner_planes.is_none()
            && !self.chart.great_circle_sector
            && matches!(self.chart.surface, Surface::Sphere(_) | Surface::Torus(_))
            && ranges
                .iter()
                .map(|(lo, hi)| (hi - lo + 1) as usize)
                .fold(1usize, usize::saturating_mul)
                > self.options.max_vertices
        {
            return Err(failure("face triangulation reached max_vertices"));
        }
        let bounds = self.bounds;
        let parameter = move |dimension: usize, index: i64| {
            bounds[dimension][0] + widths[dimension] * index as f64 / counts[dimension] as f64
        };
        let protected = if self.chart.great_circle_sector {
            self.encroaching_seeds(ranges, counts, parameter)?
        } else {
            BTreeSet::new()
        };
        let parameters_u = (ranges[0].0..=ranges[0].1)
            .take(if align_rails { 0 } else { usize::MAX })
            .map(|i| parameter(0, i))
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
                let uv = [u, parameter(1, j)];
                if protected.contains(&key(uv)) || self.chart.outside_corner(uv, self.dec)? {
                    continue;
                }
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

    /// Samples facets for preview refinement. These checks are not bounds.
    fn finish_preview(&mut self, face: FaceId, sense: Sense) -> Result<Vec<Triangle>, Refusal> {
        let tolerance = self.options.chord_tolerance.to_f64_nearest();
        let angle = self
            .options
            .angular_deviation
            .as_ref()
            .map(|v| v.to_f64_nearest());
        for _ in 0..32 {
            self.dec.check_deadline()?;
            let labels = self.labels()?;
            let mut pending = BTreeMap::new();
            let mut triangles = Vec::new();
            for triangle in self.cdt.inner_faces() {
                self.dec.check_deadline()?;
                if !labels[triangle.index()] {
                    continue;
                }
                let uv = triangle.positions().map(|p| [p.x, p.y]);
                let center = std::array::from_fn(|i| (uv[0][i] + uv[1][i] + uv[2][i]) / 3.0);
                if !self.central(center) {
                    continue;
                }
                let mut vertices = triangle.vertices().map(|vertex| vertex.data().mesh);
                if vertices[0] == vertices[1]
                    || vertices[1] == vertices[2]
                    || vertices[0] == vertices[2]
                {
                    continue;
                }
                let p = vertices.map(|i| self.positions[i]);
                let ab: [f64; 3] = std::array::from_fn(|i| p[1][i] - p[0][i]);
                let ac: [f64; 3] = std::array::from_fn(|i| p[2][i] - p[0][i]);
                let cross: [f64; 3] = std::array::from_fn(|i| {
                    ab[(i + 1) % 3] * ac[(i + 2) % 3] - ab[(i + 2) % 3] * ac[(i + 1) % 3]
                });
                let norm = cross.iter().map(|x| x * x).sum::<f64>().sqrt();
                if norm == 0.0 {
                    continue;
                }
                let mut refine = false;
                for weights in [
                    [0.5, 0.5, 0.0],
                    [0.0, 0.5, 0.5],
                    [0.5, 0.0, 0.5],
                    [1.0 / 3.0; 3],
                ] {
                    let parameter =
                        std::array::from_fn(|i| (0..3).map(|j| weights[j] * uv[j][i]).sum::<f64>());
                    let actual = self
                        .chart
                        .preview_position(self.chart.canonical(parameter), self.dec)?;
                    let linear: [f64; 3] =
                        std::array::from_fn(|i| (0..3).map(|j| weights[j] * p[j][i]).sum::<f64>());
                    let distance = (0..3)
                        .map(|i| (actual[i] - linear[i]).powi(2))
                        .sum::<f64>()
                        .sqrt();
                    refine |= distance > tolerance;
                    if let Some(angle) = angle {
                        let (normal, _) = self
                            .chart
                            .filtered_normal(parameter)
                            .ok_or_else(|| failure("preview normal is not finite"))?;
                        let normal = normal.map(|v| v.mid());
                        let length = normal.iter().map(|x| x * x).sum::<f64>().sqrt();
                        let cosine =
                            (0..3).map(|i| cross[i] * normal[i]).sum::<f64>() / (norm * length);
                        refine |= !cosine.is_finite() || cosine.clamp(-1.0, 1.0).acos() > angle;
                    }
                }
                if refine {
                    let edge = triangle
                        .adjacent_edges()
                        .into_iter()
                        .filter(|edge| !edge.is_constraint_edge())
                        .map(|edge| {
                            let a = edge.from().position();
                            let b = edge.to().position();
                            let span = (a.x - b.x).abs()
                                + if self.chart.periodic[1] {
                                    (a.y - b.y).abs()
                                } else {
                                    0.0
                                };
                            (span, [(a.x + b.x) * 0.5, (a.y + b.y) * 0.5])
                        })
                        .max_by(|a, b| a.0.total_cmp(&b.0));
                    let point = edge.map_or(center, |(_, point)| point);
                    let point = self.chart.canonical(point);
                    pending.insert(key(point), point);
                    continue;
                }
                if sense == Sense::Reversed {
                    vertices.swap(1, 2);
                }
                let first = (0..3).min_by_key(|&i| vertices[i]).unwrap();
                vertices.rotate_left(first);
                triangles.push(Triangle { vertices, face });
            }
            if pending.is_empty() {
                triangles.sort_unstable_by_key(|t| t.vertices);
                return Ok(triangles);
            }
            let mut added = false;
            for point in pending.into_values() {
                for shift in self.shifts() {
                    added |= self.free_point([point[0] + shift[0], point[1] + shift[1]])?;
                }
            }
            if !added {
                return Err(failure(
                    "preview sampling cannot resolve its visual targets",
                ));
            }
        }
        Err(failure("preview sampling reached its refinement limit"))
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
        unproved: false,
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

/// Whether one loop of three finite arcs bounds this face, each on a great circle
/// of its sphere. This only chooses a meshing algorithm, so another trim stays
/// valid and an inconclusive comparison keeps the general route.
fn great_circle_sector(model: &Model, face: FaceId, sphere: &Sphere, dec: &mut Decider) -> bool {
    let [boundary] = model.face(face).loops.as_slice() else {
        return false;
    };
    let Loop::Coedges(coedges) = model.loop_(*boundary) else {
        return false;
    };
    if coedges.len() != 3 {
        return false;
    }
    let radius_sq = sphere.radius.square();
    let mut equal = |a: &Real, b: &Real| matches!(dec.sign(&(a - b)), Outcome::ProvenCoincident(_));
    coedges.iter().all(|&coedge| {
        let edge = model.edge(model.coedge(coedge).edge);
        let Curve::Circle(circle) = model.geom().curve(edge.curve) else {
            return false;
        };
        edge.start.is_some()
            && edge.end.is_some()
            && equal(&circle.center.x, &sphere.center.x)
            && equal(&circle.center.y, &sphere.center.y)
            && equal(&circle.center.z, &sphere.center.z)
            && equal(&circle.radius_sq, &radius_sq)
    })
}

/// The local frame and bounding hemispheres of one convex corner of a sphere.
struct ConvexCorner {
    axis: Vector3,
    ref_dir: Vector3,
    planes: [Vector3; 3],
}

/// Proves that a face passing [`great_circle_sector`] is the positive cone of
/// its three corner rays. This only chooses a meshing algorithm, so a negative or
/// inconclusive sign keeps the general route, as the same arcs around any other
/// region must.
///
/// Travel along an arc turns about `m`, the carrier normal when travel follows
/// the curve. Seen from outside the sphere, `m` points to the left of travel. For
/// an arc from ray `s` to ray `e`, both on the carrier, `m . (s x e) > 0` makes it
/// the minor arc and `m` a positive multiple of `s x e`. With each opposite ray on
/// the left, the three left hemispheres meet exactly in the rays' positive cone,
/// which the three minor arcs bound. Material lies left of the loop on a face with
/// the surface's sense, so the face is that cone and not its complement.
fn convex_corner(
    model: &Model,
    face: FaceId,
    sphere: &Sphere,
    dec: &mut Decider,
) -> Option<ConvexCorner> {
    let definition = model.face(face);
    let [boundary] = definition.loops.as_slice() else {
        return None;
    };
    let Loop::Coedges(coedges) = model.loop_(*boundary) else {
        return None;
    };
    let [first, second, third] = coedges.as_slice() else {
        return None;
    };
    if definition.sense != Sense::Same {
        return None;
    }
    let arc = |coedge| {
        let coedge = model.coedge(coedge);
        let edge = model.edge(coedge.edge);
        let Curve::Circle(circle) = model.geom().curve(edge.curve) else {
            return None;
        };
        let (from, to) = if coedge.sense == Sense::Same {
            (edge.start?, edge.end?)
        } else {
            (edge.end?, edge.start?)
        };
        let normal = if coedge.sense == edge.sense {
            circle.normal.clone()
        } else {
            circle.normal.scaled(&Real::from(-1))
        };
        Some((from, to, normal))
    };
    let arcs = [arc(*first)?, arc(*second)?, arc(*third)?];
    let starts = arcs.each_ref().map(|arc| arc.0);
    if starts[0] == starts[1]
        || starts[1] == starts[2]
        || starts[2] == starts[0]
        || (0..3).any(|i| arcs[i].1 != starts[(i + 1) % 3])
    {
        return None;
    }
    let rays = starts.map(|vertex| {
        model
            .geom()
            .position(model.vertex(vertex).point)
            .minus(&sphere.center)
    });
    let mut zero = |x: Real| matches!(dec.sign(&x), Outcome::ProvenCoincident(_));
    for i in 0..3 {
        let normal = &arcs[i].2;
        if !zero(normal.dot(&rays[i])) || !zero(normal.dot(&rays[(i + 1) % 3])) {
            return None;
        }
    }
    let mut positive = |x: Real| matches!(dec.sign(&x), Outcome::Proven(StrictSign::Positive));
    if !positive(rays[0].dot(&rays[1].cross(&rays[2]))) {
        return None;
    }
    for i in 0..3 {
        let [start, end, opposite] = std::array::from_fn(|k| &rays[(i + k) % 3]);
        let normal = &arcs[i].2;
        if !positive(normal.dot(&start.cross(end))) || !positive(normal.dot(opposite)) {
            return None;
        }
    }
    // Vertex handles fix one order. In any order, the mixed-sign coordinates of
    // `a - b + c` in the independent ray basis keep both poles out of the cone,
    // and none of them is zero, so neither pole lies on a boundary plane.
    let mut ordered: [_; 3] = std::array::from_fn(|i| (starts[i], &rays[i]));
    ordered.sort_unstable_by_key(|(vertex, _)| *vertex);
    let [(_, a), (_, b), (_, c)] = ordered;
    let axis = a.minus(b).plus(c);
    let ref_dir = axis.cross(a);
    let [(_, _, x), (_, _, y), (_, _, z)] = arcs;
    Some(ConvexCorner {
        axis,
        ref_dir,
        planes: [x, y, z],
    })
}

fn poles_clear(chart: &Chart, equations: &[Vec<RealPoly>], dec: &Decider) -> Result<bool, Refusal> {
    for pole in [[0.0, -0.25], [0.0, 0.25]] {
        let pole = chart.evaluate(pole, dec)?;
        for equations in equations {
            let mut apart = false;
            for equation in equations {
                if !equation.eval_ball(&pole, dec)?.contains_zero() {
                    apart = true;
                    break;
                }
            }
            if !apart {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

fn regular_chart(
    model: &Model,
    face: FaceId,
    options: &MeshOptions,
    dec: &mut Decider,
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
    let sector =
        options.angular_deviation.is_some() && great_circle_sector(model, face, sphere, dec);
    if sector && let Some(corner) = convex_corner(model, face, sphere, dec) {
        let mut local = sphere.clone();
        local.axis = corner.axis;
        local.ref_dir = corner.ref_dir;
        let mut chart = Chart::new(&Surface::Sphere(local), dec)?;
        if poles_clear(&chart, &equations, dec)? {
            let [x, y, z] = corner
                .planes
                .each_ref()
                .map(|normal| crate::evaluate::unit(normal, dec));
            chart.great_circle_sector = true;
            chart.corner_planes = Some([x?, y?, z?]);
            return Ok(chart);
        }
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
        let mut chart = Chart::new(&Surface::Sphere(candidate), dec)?;
        if poles_clear(&chart, &equations, dec)? {
            chart.great_circle_sector = sector;
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
    let junction = crate::junction::prepare(model, &faces, options, dec)?;
    if junction.is_none()
        && faces.iter().any(|face| {
            matches!(
                model.geom().surface(model.face(*face).surface),
                Surface::BSpline(_)
            )
        })
    {
        return crate::spline_prism::build(model, body, &faces, options, dec);
    }
    let fitted = |face: FaceId| junction.as_ref().is_some_and(|j| j.contains(face));
    let mut charts = faces
        .iter()
        .filter(|&&f| !fitted(f))
        .map(|&f| Ok((f, regular_chart(model, f, options, dec)?)))
        .collect::<Result<BTreeMap<_, _>, Refusal>>()?;
    let mut boundaries = edges::sample(model, body, &charts, junction.as_ref(), options, dec)?;
    let boundary_limit = boundaries.mesh.positions.len();
    let mut prior_chords = BTreeSet::new();
    for face in faces {
        let start = boundaries.mesh.triangles.len();
        if fitted(face) {
            junction
                .as_ref()
                .expect("fitted face")
                .triangulate(face, &mut boundaries, options)?;
        } else if let Surface::Plane(plane) = &charts[&face].surface {
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
                let positions = boundaries.mesh.positions.len();
                let result = face_mesh(
                    model,
                    face,
                    charts.get_mut(&face).unwrap(),
                    &mut boundaries,
                    options,
                    &prior_chords,
                    dec,
                );
                match result {
                    Err(error)
                        if error.code == "mesh.precision-or-budget"
                            && crate::periodic::eligible(model, face, &charts[&face], options) =>
                    {
                        boundaries.mesh.positions.truncate(positions);
                        crate::periodic::build(
                            model,
                            face,
                            &charts[&face],
                            &mut boundaries,
                            options,
                            dec,
                        )?;
                    }
                    other => other?,
                }
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

/// Shares the chart and boundary machinery, but never runs facet certification.
pub(crate) fn build_preview(
    model: &Model,
    body: BodyId,
    options: &crate::PreviewOptions,
    dec: &mut Decider,
) -> Result<crate::PreviewMesh, Refusal> {
    let visual = options;
    let options = MeshOptions {
        chord_tolerance: options.chord_tolerance.clone(),
        angular_deviation: options.angular_deviation.clone(),
        max_vertices: options.max_vertices,
        max_edge_length: None,
        max_aspect_ratio: None,
    };
    let faces = model.faces_of_body(body);
    let mut charts = faces
        .iter()
        .map(|&face| Ok((face, regular_chart(model, face, &options, dec)?)))
        .collect::<Result<BTreeMap<_, _>, Refusal>>()?;
    let mut boundaries = crate::preview::edges::build(model, body, &charts, visual, dec)?;
    for face in faces {
        crate::preview::strips::restrict(model, face, charts.get_mut(&face).unwrap(), &boundaries);
        if crate::preview::strips::build(model, face, &charts[&face], &mut boundaries, visual, dec)?
        {
            continue;
        }
        let chart = &charts[&face];
        let definition = model.face(face);
        let mut segments = Vec::new();
        let mut projected = Vec::new();
        for &boundary in &definition.loops {
            match model.loop_(boundary) {
                Loop::Vertex(vertex) => projected
                    .push(chart.project(boundaries.mesh.positions[boundaries.vertices[vertex]])),
                Loop::Coedges(coedges) => {
                    for &coedge in coedges {
                        let coedge = model.coedge(coedge);
                        let edge = &boundaries.mesh.edges[boundaries.edges[&coedge.edge]];
                        for pair in edge.vertices.windows(2) {
                            let [first, last] = if coedge.sense == Sense::Same {
                                [pair[0], pair[1]]
                            } else {
                                [pair[1], pair[0]]
                            };
                            if first == last && boundaries.merged_vertices {
                                continue;
                            }
                            let a = chart.project(boundaries.mesh.positions[first]);
                            let b = align(
                                chart.project(boundaries.mesh.positions[last]),
                                a,
                                chart.periodic,
                            );
                            segments.push((a, b, first, last));
                        }
                    }
                }
            }
        }
        if boundaries.merged_vertices {
            // Only simplified boundaries need cancellation. Preserve the usual
            // insertion order and avoid a second index for ordinary previews.
            let mut surviving = BTreeMap::new();
            for (a, b, first, last) in segments {
                // Include chart lifts: periodic arcs can share mesh endpoints
                // without tracing the same segment in the chart.
                if surviving.remove(&(last, first, key(b), key(a))).is_none() {
                    surviving.insert((first, last, key(a), key(b)), (a, b, first, last));
                }
            }
            segments = surviving.into_values().collect();
        }
        if boundaries.merged_vertices
            && segments.is_empty()
            && definition
                .loops
                .iter()
                .any(|&boundary| matches!(model.loop_(boundary), Loop::Coedges(_)))
        {
            let end = boundaries.mesh.triangles.len();
            boundaries.mesh.faces.push(crate::PreviewFace {
                face,
                triangles: end..end,
            });
            continue;
        }
        projected.extend(segments.iter().flat_map(|&(a, b, _, _)| [a, b]));
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
            return Err(failure("preview face has no finite parameter region"));
        }
        let hessian = chart.hessian(bounds[1][0].abs().max(bounds[1][1].abs()), dec)?;
        let prior_chords = BTreeSet::new();
        let rails = edges::Rails::default();
        let mut builder = FaceBuilder {
            chart,
            options: &options,
            unproved: true,
            positions: &mut boundaries.mesh.positions,
            rails: &rails,
            cdt: Cdt::new(),
            canonical: BTreeMap::new(),
            errors: Vec::new(),
            filtered_errors: Vec::new(),
            interior_left: BTreeSet::new(),
            prior_chords: &prior_chords,
            bounds,
            filtered_hessian: hessian.each_ref().map(Interval::from_rational),
            hessian: hessian.into_owned(),
            dec,
        };
        builder.boundary(&segments, definition.sense)?;
        builder.grid()?;
        let triangles = builder.finish_preview(face, definition.sense)?;
        let start = boundaries.mesh.triangles.len();
        boundaries.mesh.triangles.extend(triangles);
        boundaries.mesh.faces.push(crate::PreviewFace {
            face,
            triangles: start..boundaries.mesh.triangles.len(),
        });
    }
    Ok(boundaries.mesh)
}
