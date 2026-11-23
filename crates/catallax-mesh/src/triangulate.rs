use std::collections::{BTreeMap, BTreeSet, VecDeque};

use catallax_numbers::{
    Answer, Answered, Ball, Decider, Interval, Outcome, Rational, Real, RealPoly, Refusal,
    StrictSign,
};
use catallax_shapes::{
    Curve, Sphere, SplineSupport, Surface, Vector3, tessellation::checked_triangle_bound,
};
use catallax_topology::{BodyId, EdgeId, FaceId, Loop, Model, Sense};
use spade::{
    ConstrainedDelaunayTriangulation, FloatTriangulation, HasPosition, HierarchyHintGenerator,
    Point2, PositionInTriangulation, Triangulation,
};

use crate::{
    FaceMesh, Mesh, MeshOptions, Triangle,
    chart::{self, Chart, CollarBridge, NormalProof, SplineChart, Uv, align, chord_bound, key},
    edges::{self, Boundaries},
    failure,
};

#[derive(Clone, Copy)]
struct Node {
    uv: Point2<f64>,
    // Where the triangulation sees this node. See `FaceBuilder::stretch`.
    at: Point2<f64>,
    mesh: usize,
    error: usize,
}
impl HasPosition for Node {
    type Scalar = f64;
    fn position(&self) -> Point2<f64> {
        self.at
    }
}
type Cdt = ConstrainedDelaunayTriangulation<Node, (), (), (), HierarchyHintGenerator<f64>>;

type Segment = (Uv, Uv, usize, usize);

#[derive(Clone, Copy)]
struct BandCut {
    ends: [(Uv, usize); 2],
    sense: Sense,
}

impl BandCut {
    // The rails already have an authoritative set of samples. Rotate each
    // loop to an existing sample and lift its edges into one finite band.
    fn open(segments: &mut [Segment]) -> Result<Self, Refusal> {
        let choose = |row: f64, near: f64| {
            segments
                .iter()
                .filter(|(a, _, _, _)| a[1].to_bits() == row.to_bits())
                .map(|&(uv, _, index, _)| (align(uv, [near, row], [true, false]), index))
                .min_by(|(a, _), (b, _)| (a[0] - near).abs().total_cmp(&(b[0] - near).abs()))
                .ok_or_else(|| failure("periodic band is missing a sampled rail"))
        };
        let lower = choose(0.0, 0.0)?;
        let upper = choose(1.0, lower.0[0])?;
        let direction = segments
            .iter()
            .find(|(a, _, _, _)| a[1] == 0.0)
            .is_some_and(|(a, b, _, _)| a[0] < b[0]);
        let sense = if direction {
            Sense::Same
        } else {
            Sense::Reversed
        };
        let ends = [lower, upper];
        let mut winding = [0.0; 2];
        for (a, b, _, _) in segments {
            let row = usize::from(a[1].to_bits() == 1.0_f64.to_bits());
            if a[1].to_bits() != (row as f64).to_bits() || b[1].to_bits() != a[1].to_bits() {
                return Err(failure("periodic band rail leaves its chart row"));
            }
            let step = b[0] - a[0];
            if step == 0.0 || (step > 0.0) != (direction != (row == 1)) {
                return Err(failure(
                    "periodic band rails must wind in opposite directions",
                ));
            }
            winding[row] += step;
            let left = ends[row].0[0];
            let shift = libm::ceil(left - a[0].min(b[0]));
            a[0] += shift;
            b[0] += shift;
            if a[0].max(b[0]) > left + 1.0 {
                return Err(failure("periodic band seam crosses a rail segment"));
            }
        }
        if winding[0].abs().to_bits() != 1.0_f64.to_bits()
            || winding[1].abs().to_bits() != 1.0_f64.to_bits()
        {
            return Err(failure("periodic band rails must each wind once"));
        }
        Ok(Self { ends, sense })
    }

    fn lift(self, mut uv: Uv) -> Uv {
        let [lo, hi] = self.ends.map(|(end, _)| end[0]);
        let left = lo + (hi - lo) * uv[1];
        uv[0] += libm::ceil(left - uv[0]);
        uv
    }
}

struct FaceBuilder<'a> {
    chart: &'a Chart,
    options: &'a MeshOptions,
    band: Option<BandCut>,
    collapsed_station: Option<f64>,
    // A cylinder trimmed obliquely, which gets ruled columns and rows.
    ruled: bool,
    // Delaunay distances use the seed-cell metric, not mixed turns and
    // millimetres. Samples and all certificates retain their original UVs.
    stretch: f64,
    unproved: bool,
    positions: &'a mut Vec<[f64; 3]>,
    rails: &'a edges::Rails,
    cdt: Cdt,
    canonical: BTreeMap<[u64; 2], (usize, usize)>,
    errors: Vec<Rational>,
    filtered_errors: Vec<Interval>,
    interior_left: BTreeSet<(usize, usize)>,
    collar_edges: BTreeSet<[usize; 2]>,
    protected_trim: Vec<[Point2<f64>; 2]>,
    separation: f64,
    prior_chords: &'a BTreeSet<(usize, usize)>,
    bounds: [[f64; 2]; 2],
    hessian: [Rational; 3],
    filtered_hessian: [Interval; 3],
    dec: &'a mut Decider,
}

// A band refinement point may enter only if every fixed segment's disk
// excludes it. The interval comparison is allowed to decline a tangent point.
fn clears_constraints(at: Point2<f64>, segments: &[[Point2<f64>; 2]]) -> bool {
    let point = |x| Interval::point(x).unwrap();
    let x = point(at.x);
    let y = point(at.y);
    segments.iter().all(|[a, b]| {
        if (at.x <= a.x.min(b.x) || at.x >= a.x.max(b.x))
            && (at.y <= a.y.min(b.y) || at.y >= a.y.max(b.y))
        {
            return true;
        }
        let dot = (x - point(a.x)) * (x - point(b.x)) + (y - point(a.y)) * (y - point(b.y));
        dot.lo() >= 0.0
    })
}

impl FaceBuilder<'_> {
    fn protects_band(&self) -> bool {
        self.band.is_some()
            && matches!(
                self.chart.surface,
                Surface::BSpline(SplineSupport::RollingSaddle(_))
            )
    }

    fn spline_position(&mut self, uv: Uv, at: Point2<f64>) -> Result<[Interval; 3], Refusal> {
        let spline = self.chart.spline().expect("spline chart");
        if self.collapsed_station.is_some() {
            spline.evaluate_collapsed([at.x, at.y], self.dec)
        } else {
            spline.evaluate(uv, self.dec)
        }
    }

    fn node(&mut self, uv: Uv, boundary: Option<usize>) -> Result<Node, Refusal> {
        if self.positions.len() >= self.options.max_vertices {
            return Err(failure("face triangulation reached max_vertices"));
        }
        let canonical = self.chart.canonical(uv);
        let at = self.place(uv);
        let (mesh, error) = if let Some(&(mesh, error)) = self.canonical.get(&key(canonical)) {
            if boundary.is_some_and(|index| index != mesh) {
                return Err(failure(
                    "distinct boundary vertices collapse in the display chart",
                ));
            }
            (mesh, error)
        } else if self.unproved && self.chart.spline().is_none() {
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
            let replacement = boundary.map(|index| self.positions[index]);
            let (p, residual, filtered) = match self.chart.spline() {
                Some(spline) => {
                    let (p, bound) = if self.collapsed_station.is_some() {
                        let enclosed = self.spline_position(canonical, at)?;
                        SplineChart::rounded_vertex(enclosed, replacement)
                    } else {
                        spline.vertex(canonical, replacement, self.dec)?
                    };
                    if !bound.is_finite() {
                        return Err(failure("spline sample residual exceeds binary64"));
                    }
                    let residual = if self.unproved {
                        Rational::zero()
                    } else {
                        chart::q(bound.hi())
                    };
                    let filtered = if self.unproved || spline.patch_sampling() {
                        bound
                    } else {
                        Interval::from_rational(&residual)
                    };
                    (p, residual, filtered)
                }
                None => {
                    let (p, residual) = self.chart.sample(
                        canonical,
                        replacement,
                        &self.options.chord_tolerance,
                        self.dec,
                    )?;
                    let filtered = Interval::from_rational(&residual);
                    (p, residual, filtered)
                }
            };
            let mesh = if let Some(index) = boundary {
                index
            } else {
                let index = self.positions.len();
                self.positions.push(p);
                index
            };
            let error = self.errors.len();
            self.filtered_errors.push(filtered);
            self.errors.push(residual);
            self.canonical.insert(key(canonical), (mesh, error));
            (mesh, error)
        };
        Ok(Node {
            uv: Point2::new(uv[0], uv[1]),
            at,
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
        let at = self.place(uv);
        let placed_key = key([at.x, at.y]);
        if let Some(&index) = indices.get(&placed_key) {
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
        indices.insert(placed_key, index);
        Ok(index)
    }

    fn bulk_load(&mut self, nodes: Vec<Node>) -> Result<(), Refusal> {
        let constraints = self
            .interior_left
            .iter()
            .map(|&(a, b)| [a, b])
            .chain(self.collar_edges.iter().copied())
            .collect();
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

    fn tiled(&self) -> [bool; 2] {
        if self.band.is_some() {
            [false; 2]
        } else {
            self.chart.periodic
        }
    }

    fn refinement_uv(&self, uv: Uv, round: bool) -> Uv {
        let mut uv = if self.band.is_some() {
            uv
        } else {
            self.chart.canonical(uv)
        };
        if round || self.band.is_some() {
            for (value, periodic) in uv.iter_mut().zip(self.chart.periodic) {
                if periodic {
                    *value = chart::quantize_turn(*value);
                }
            }
        }
        uv
    }

    fn shifts(&self) -> impl Iterator<Item = Uv> + use<> {
        let u = if self.tiled()[0] { -1..=1 } else { 0..=0 };
        let v = if self.tiled()[1] { -1..=1 } else { 0..=0 };
        u.flat_map(move |i| v.clone().map(move |j| [f64::from(i), f64::from(j)]))
    }

    fn boundary(&mut self, segments: &[Segment], sense: Sense) -> Result<(), Refusal> {
        let mut nodes = Vec::new();
        let mut indices = BTreeMap::new();
        for &(mut a, mut b, first, last) in segments {
            for i in 0..2 {
                if self.tiled()[i] {
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
        let anchors = if self.protects_band() {
            segments
                .iter()
                .map(|&(a, b, _, _)| {
                    SplineChart::rail_anchor(a, b, self.stretch)
                        .ok_or_else(|| failure("shared edge has no planned collar altitude"))
                })
                .collect::<Result<Vec<_>, _>>()?
        } else {
            Vec::new()
        };
        let seam_collars = if let Some(band) = self.band {
            let caps = if self.protects_band() && self.options.angular_deviation.is_some() {
                Some(Self::seam_collar_points(band, segments, &anchors)?)
            } else {
                None
            };
            let seam = self.seam(band, &anchors, caps, &mut nodes, &mut indices)?;
            caps.map(|points| (points, seam))
        } else {
            None
        };
        if self.protects_band() {
            self.protected_trim = self
                .interior_left
                .iter()
                .map(|&(a, b)| [nodes[a].at, nodes[b].at])
                .collect();
            self.separation = self
                .protected_trim
                .iter()
                .map(|[a, b]| libm::hypot(b.x - a.x, b.y - a.y))
                .fold(f64::INFINITY, f64::min)
                / 64.0;
            if !self.separation.is_finite() || self.separation <= 0.0 {
                return Err(failure("periodic band has no positive trim separation"));
            }
            self.collar(segments, anchors, &mut nodes, &mut indices, seam_collars)?;
        }
        self.bulk_load(nodes)
    }

    fn seam_collar_points(
        band: BandCut,
        rails: &[Segment],
        anchors: &[Uv],
    ) -> Result<[Uv; 2], Refusal> {
        let [lo, hi] = band.ends.map(|(uv, _)| uv[0]);
        let at = |v| lo + (hi - lo) * v;
        let mut caps = [[0.0; 2]; 2];
        for (row, (start, _)) in band.ends.into_iter().enumerate() {
            let find_anchor = |end| {
                let end = key(end);
                rails
                    .iter()
                    .zip(anchors)
                    .find_map(|(&(a, b, _, _), &apex)| {
                        (key(a) == end || key(b) == end).then_some(apex)
                    })
            };
            let right = find_anchor(start)
                .ok_or_else(|| failure("seam collar has no first rail segment"))?;
            let mut left = find_anchor([start[0] + 1.0, start[1]])
                .ok_or_else(|| failure("seam collar has no last rail segment"))?;
            left[0] -= 1.0;
            let x = at(left[1]) - left[0];
            let y = right[0] - at(right[1]);
            if x <= 0.0 || y <= 0.0 {
                return Err(failure("seam collar apex is outside its copy of the band"));
            }
            let v = left[1] + (right[1] - left[1]) * (x / (x + y));
            caps[row] = [chart::quantize_turn(at(v)), v];
        }
        Ok(caps)
    }

    fn collar(
        &mut self,
        rails: &[Segment],
        anchors: Vec<Uv>,
        nodes: &mut Vec<Node>,
        indices: &mut BTreeMap<[u64; 2], usize>,
        seam_collars: Option<([Uv; 2], Vec<[usize; 2]>)>,
    ) -> Result<(), Refusal> {
        let band = self.band.expect("saddle band");
        let [lo, hi] = band.ends.map(|(uv, _)| uv[0]);
        #[derive(Clone, Copy)]
        struct Neighbor {
            apex: usize,
            endpoint: usize,
        }
        let mut neighbors = BTreeMap::<usize, [Option<Neighbor>; 2]>::new();
        for (&(a, b, _, _), inward) in rails.iter().zip(anchors) {
            let at = self.place(inward);
            let rho = self.separation;
            let x = (at.x - rho).next_down()..=(at.x + rho).next_up();
            let y = (at.y - rho).next_down()..=(at.y + rho).next_up();
            if !chart::positive_turn([hi, 1.0], [lo, 0.0], inward)
                || !chart::positive_turn([lo + 1.0, 0.0], [hi + 1.0, 1.0], inward)
                || !clears_constraints(at, &self.protected_trim)
                || !(a[0].min(b[0]) < inward[0] && inward[0] < a[0].max(b[0]))
                || nodes
                    .iter()
                    .any(|node| x.contains(&node.at.x) && y.contains(&node.at.y))
            {
                return Err(failure("periodic band collar encroaches on its trim"));
            }
            let ends = [a, b].map(|uv| {
                let at = self.place(uv);
                indices[&key([at.x, at.y])]
            });
            let middle = self.queue(nodes, indices, inward, None)?;
            for end in ends {
                self.collar_edges.insert([end.min(middle), end.max(middle)]);
            }
            if self.options.angular_deviation.is_some() {
                let [left, right] = if a[0] < b[0] {
                    ends
                } else {
                    [ends[1], ends[0]]
                };
                neighbors.entry(left).or_default()[1] = Some(Neighbor {
                    apex: middle,
                    endpoint: right,
                });
                neighbors.entry(right).or_default()[0] = Some(Neighbor {
                    apex: middle,
                    endpoint: left,
                });
            }
        }
        let Some(angle) = &self.options.angular_deviation else {
            return Ok(());
        };
        let cosine = crate::controls::cosine_enclosure(angle, self.dec)?;
        let spline = self.chart.spline().expect("saddle band chart");
        let uv = |i: usize| [nodes[i].uv.x, nodes[i].uv.y];
        if let Some((caps, seam)) = seam_collars {
            for (row, cap) in caps.into_iter().enumerate() {
                let tip = seam
                    .iter()
                    .position(|pair| (nodes[pair[0]].uv.y - cap[1]) == 0.0)
                    .ok_or_else(|| failure("seam collar has no planned tip"))?;
                let end = if row == 0 { 0 } else { seam.len() - 1 };
                for side in 0..2 {
                    let apex = neighbors[&seam[end][side]]
                        .into_iter()
                        .flatten()
                        .next()
                        .expect("one collar apex at each seam end")
                        .apex;
                    for pair in seam[end.min(tip)..=end.max(tip)].windows(2) {
                        let mut triangle = [pair[0][side], pair[1][side], apex];
                        let [a, b, c] = triangle.map(uv);
                        if !chart::positive_turn(a, b, c) {
                            triangle.swap(0, 1);
                        }
                        let vertices = triangle.map(|i| self.positions[nodes[i].mesh]);
                        let residual = triangle
                            .iter()
                            .map(|&i| self.filtered_errors[nodes[i].error].hi())
                            .fold(0.0_f64, f64::max);
                        if !spline.collar_facet_fits(
                            &triangle.map(uv),
                            vertices,
                            residual,
                            self.options,
                            cosine,
                            self.dec,
                        )? {
                            return Err(failure(format!(
                                "periodic seam collar exceeds its facet certificate at {:?}",
                                triangle.map(uv),
                            )));
                        }
                        for endpoint in [pair[0][side], pair[1][side]] {
                            self.collar_edges
                                .insert([apex.min(endpoint), apex.max(endpoint)]);
                        }
                    }
                }
            }
        }
        for (rail, [left, right]) in neighbors {
            let (Some(left), Some(right)) = (left, right) else {
                continue;
            };
            let apex_indices = [left.apex, right.apex];
            let triangle = [rail, left.apex, right.apex];
            let bridge = CollarBridge {
                rail: [left.endpoint, rail, right.endpoint].map(uv),
                apices: apex_indices.map(uv),
                positions: triangle.map(|i| self.positions[nodes[i].mesh]),
                residual: triangle
                    .into_iter()
                    .map(|i| self.filtered_errors[nodes[i].error].hi())
                    .fold(0.0_f64, f64::max),
            };
            if spline.collar_bridge_fits(bridge, self.options, cosine, self.dec)? {
                self.collar_edges
                    .insert([left.apex.min(right.apex), left.apex.max(right.apex)]);
            }
        }
        Ok(())
    }

    fn outside_collar(
        &self,
        face: spade::handles::FixedFaceHandle<spade::handles::InnerTag>,
    ) -> bool {
        self.cdt
            .face(face)
            .adjacent_edges()
            .iter()
            .any(|edge| !edge.is_constraint_edge())
    }

    fn separated(&self, at: Point2<f64>) -> bool {
        let rho = self.separation;
        let lower = Point2::new((at.x - rho).next_down(), (at.y - rho).next_down());
        let upper = Point2::new((at.x + rho).next_up(), (at.y + rho).next_up());
        self.cdt
            .get_vertices_in_rectangle(lower, upper)
            .next()
            .is_none()
    }

    fn admissible(&self, uv: Uv, labels: &[bool]) -> bool {
        let at = self.place(uv);
        let inside = match self.cdt.locate(at) {
            PositionInTriangulation::OnFace(f) => labels[f.index()] && self.outside_collar(f),
            PositionInTriangulation::OnEdge(e) => {
                !self.cdt.is_constraint_edge(e.as_undirected())
                    && labels[self.cdt.directed_edge(e).face().index()]
            }
            _ => false,
        };
        inside && clears_constraints(at, &self.protected_trim) && self.separated(at)
    }

    fn seam(
        &mut self,
        band: BandCut,
        anchors: &[Uv],
        caps: Option<[Uv; 2]>,
        nodes: &mut Vec<Node>,
        indices: &mut BTreeMap<[u64; 2], usize>,
    ) -> Result<Vec<[usize; 2]>, Refusal> {
        let spline = self
            .chart
            .spline()
            .expect("a finite band has a spline chart");
        let chord = Interval::from_rational(&self.options.chord_tolerance.mul_pow2(-3));
        let cosine = self
            .options
            .angular_deviation
            .as_ref()
            .map(|angle| {
                if self.protects_band() {
                    crate::controls::cosine_enclosure(angle, self.dec)
                } else {
                    crate::controls::cosine_limit(angle, self.dec)
                        .map(|q| Interval::from_rational(&q))
                }
            })
            .transpose()?;
        let [(lo, lo_index), (hi, hi_index)] = band.ends;
        let anchor_positions: Vec<_> = anchors.iter().map(|&uv| self.place(uv)).collect();
        let mut points = vec![lo];
        let mut pending = if let Some([a, b]) = caps {
            vec![(b, hi, 2), (a, b, 2), (lo, a, 2)]
        } else {
            vec![(lo, hi, 0)]
        };
        while let Some((a, b, depth)) = pending.pop() {
            let middle = [
                chart::quantize_turn(a[0].midpoint(b[0])),
                a[1].midpoint(b[1]),
            ];
            let within_length = if let Some(limit) = &self.options.max_edge_length {
                let [a, b] = [a, b].map(|uv| {
                    let uv = self.chart.canonical(uv);
                    self.canonical.get(&key(uv)).map_or_else(
                        || spline.evaluate(uv, self.dec).map(|p| p.map(|x| x.mid())),
                        |&(index, _)| Ok(self.positions[index]),
                    )
                });
                crate::controls::length_squared(a?, b?) <= limit * limit
            } else {
                true
            };
            let cuts = [
                [self.place(a), self.place(b)],
                [
                    self.place([a[0] + 1.0, a[1]]),
                    self.place([b[0] + 1.0, b[1]]),
                ],
            ];
            // Both seam copies use the stricter subdivision. Resolve their
            // encroachment before either one becomes a fixed constraint.
            let clear = anchor_positions
                .iter()
                .all(|&at| clears_constraints(at, &cuts));
            if clear
                && within_length
                && spline.seam_within([a, middle, b], chord, cosine, self.dec)?
            {
                points.push(b);
                continue;
            }
            if depth == 48
                || points.len().saturating_add(pending.len()) >= self.options.max_vertices
            {
                return Err(failure("periodic band seam exceeds the subdivision limit"));
            }
            pending.push((middle, b, depth + 1));
            pending.push((a, middle, depth + 1));
        }

        // Both copies use the same canonical UV and mesh vertex. In parameter
        // space they are independent boundary constraints with opposite winding.
        let last = points.len() - 1;
        let mut seam = Vec::with_capacity(points.len());
        let mut previous: Option<[usize; 2]> = None;
        for (i, uv) in points.into_iter().enumerate() {
            let boundary = if i == 0 {
                Some(lo_index)
            } else if i == last {
                Some(hi_index)
            } else {
                None
            };
            let left = self.queue(nodes, indices, uv, boundary)?;
            let right = self.queue(nodes, indices, [uv[0] + 1.0, uv[1]], boundary)?;
            if let Some([a, b]) = previous {
                self.interior_left.extend([(left, a), (b, right)]);
            }
            previous = Some([left, right]);
            seam.push([left, right]);
        }
        Ok(seam)
    }

    fn place(&self, uv: Uv) -> Point2<f64> {
        let point = if let Some(station) = self.collapsed_station {
            let y = uv[1] - station;
            Point2::new(uv[0] * y.abs(), if y == 0.0 { 0.0 } else { y })
        } else {
            Point2::new(uv[0], uv[1] * self.stretch)
        };
        if self.unproved {
            // Keep the sampled UV; only the CDT coordinates need Spade's range.
            spade::mitigate_underflow(point)
        } else {
            point
        }
    }

    fn unplace_collapsed(station: f64, [x, y]: Uv) -> Uv {
        if y == 0.0 {
            return [0.0, station];
        }
        [(x / y.abs()).clamp(0.0, 1.0), (station + y).clamp(0.0, 1.0)]
    }

    fn midpoint(station: Option<f64>, a: Node, b: Node) -> Uv {
        if let Some(station) = station {
            Self::unplace_collapsed(station, [a.at.x.midpoint(b.at.x), a.at.y.midpoint(b.at.y)])
        } else {
            [a.uv.x.midpoint(b.uv.x), a.uv.y.midpoint(b.uv.y)]
        }
    }

    fn free_point(&mut self, uv: Uv) -> Result<bool, Refusal> {
        match self.cdt.locate(self.place(uv)) {
            PositionInTriangulation::OnVertex(_) => return Ok(false),
            PositionInTriangulation::OnEdge(e)
                if self.cdt.is_constraint_edge(e.as_undirected()) =>
            {
                return Ok(false);
            }
            PositionInTriangulation::OutsideOfConvexHull(_) if self.collapsed_station.is_some() => {
                return Ok(false);
            }
            _ => {}
        }
        self.insert(uv, None)?;
        Ok(true)
    }

    /// The same rule for a lattice of columns and rows, sorted, which the
    /// triangulation sees through `stretch`.
    fn encroached_lattice(
        &self,
        columns: &[f64],
        rows: &[f64],
    ) -> Result<BTreeSet<[u64; 2]>, Refusal> {
        let mut protected = BTreeSet::new();
        for edge in self.cdt.undirected_edges() {
            self.dec.check_deadline()?;
            if !edge.is_constraint_edge() {
                continue;
            }
            let [a, b] = edge.positions();
            let radius = ((a.x - b.x).powi(2) + (a.y - b.y).powi(2)).sqrt() * 0.5;
            let [x, y] = [(a.x + b.x) * 0.5, (a.y + b.y) * 0.5];
            let first = columns.partition_point(|u| *u < x - radius);
            for &u in columns[first..].iter().take_while(|u| **u <= x + radius) {
                let low = rows.partition_point(|v| *v * self.stretch < y - radius);
                for &v in rows[low..]
                    .iter()
                    .take_while(|v| **v * self.stretch <= y + radius)
                {
                    let at = v * self.stretch;
                    if (u - a.x) * (u - b.x) + (at - a.y) * (at - b.y) < 0.0 {
                        protected.insert(key([u, v]));
                    }
                }
            }
        }
        Ok(protected)
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
                let stretch = if i == 0 { 1.0 } else { self.stretch };
                let first = libm::floor(((centers[i] - radius) / stretch - start) * scale) as i64;
                let last = libm::ceil(((centers[i] + radius) / stretch - start) * scale) as i64;
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
                    let at = v * self.stretch;
                    if (u - a.x) * (u - b.x) + (at - a.y) * (at - b.y) < 0.0 {
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
        let mut empty_after = Vec::<bool>::new();
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
                empty_after.push(before_empty);
            }
        }
        // With no edge limit a gap between kept columns can be wider than the
        // target, so free columns fill it evenly. The last gap wraps the period.
        if self.edge_limit().is_none() {
            let mut filled = Vec::with_capacity(kept.len());
            du = 0.0;
            for (i, &u) in kept.iter().enumerate() {
                let next = kept.get(i + 1).copied().unwrap_or(kept[0] + 1.0);
                let gap = next - u;
                filled.push(u);
                if empty_after[i] || (i + 1 == kept.len() && gap > widest_gap) {
                    continue;
                }
                let parts = libm::ceil(gap / target).max(1.0);
                for part in 1..parts as usize {
                    filled.push(u + gap * part as f64 / parts);
                }
                du = du.max((gap / parts).next_up());
            }
            kept = filled;
        }
        let mut lifted: Vec<_> = self
            .shifts()
            .flat_map(|shift| kept.iter().map(move |u| u + shift[0]).collect::<Vec<_>>())
            .collect();
        lifted.sort_unstable_by(f64::total_cmp);
        (lifted, du)
    }

    fn grid(&mut self) -> Result<(), Refusal> {
        if let Some(spline) = self.chart.spline() {
            let (chord, remaining) = if spline.patch_sampling() {
                (
                    self.options.chord_tolerance.clone(),
                    self.options.max_vertices,
                )
            } else {
                (
                    self.options.chord_tolerance.mul_pow2(-1),
                    self.options
                        .max_vertices
                        .saturating_sub(self.positions.len()),
                )
            };
            let mut grid = spline.seed_grid(&chord, remaining, self.dec)?;
            for uv in &mut grid {
                *uv = self.chart.canonical(*uv);
                if let Some(band) = self.band {
                    *uv = band.lift(*uv);
                }
            }
            // Use the fixed trim before inserting any seeds. Remove those
            // outside it and protect the shared segments' diametral disks,
            // as for analytic trim lattices.
            let labels = self.labels()?;
            let protect_band = self.protects_band();
            grid.retain(|&uv| {
                if protect_band {
                    self.admissible(uv, &labels)
                } else {
                    match self.cdt.locate(self.place(uv)) {
                        PositionInTriangulation::OnFace(f) => labels[f.index()],
                        PositionInTriangulation::OnEdge(e) => {
                            !self.cdt.is_constraint_edge(e.as_undirected())
                                && labels[self.cdt.directed_edge(e).face().index()]
                        }
                        _ => false,
                    }
                }
            });
            let at_x = |uv: Uv| {
                if self.collapsed_station.is_some() {
                    self.place(uv).x
                } else {
                    uv[0]
                }
            };
            grid.sort_by(|a, b| at_x(*a).total_cmp(&at_x(*b)).then(a[1].total_cmp(&b[1])));
            let mut protected = BTreeSet::new();
            for edge in self
                .cdt
                .undirected_edges()
                .filter(|edge| !protect_band && edge.is_constraint_edge())
            {
                self.dec.check_deadline()?;
                let [a, b] = edge.positions();
                let radius = libm::hypot(a.x - b.x, a.y - b.y) * 0.5;
                let middle = (a.x + b.x) * 0.5;
                let first = grid.partition_point(|uv| at_x(*uv) < middle - radius);
                for &uv in grid[first..]
                    .iter()
                    .take_while(|uv| at_x(**uv) <= middle + radius)
                {
                    let at = self.place(uv);
                    if (at.x - a.x) * (at.x - b.x) + (at.y - a.y) * (at.y - b.y) < 0.0 {
                        protected.insert(key(uv));
                    }
                }
            }
            for uv in grid {
                if protect_band && !self.separated(self.place(uv)) {
                    continue;
                }
                if !protected.contains(&key(uv)) {
                    for shift in self.shifts() {
                        self.free_point([uv[0] + shift[0], uv[1] + shift[1]])?;
                    }
                }
            }
            return Ok(());
        }
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
        // Columns follow a cylinder's straight rulings, so a facet between two of
        // them lies flat however long it is. A free grid staggers its points
        // against a trim that runs along the rows, and its slivers stand up off
        // the surface, where centroid refinement only makes them thinner.
        let ruled = self.ruled && self.cdt.num_vertices() > 0;
        // No column gap inside a trim is wider than its widest boundary segment.
        let widest_segment = (ruled || self.edge_limit().is_some()).then(|| {
            self.cdt
                .undirected_edges()
                .filter(|edge| edge.is_constraint_edge())
                .map(|edge| {
                    let [a, b] = edge.positions();
                    (a.x - b.x).abs()
                })
                .fold(0.0_f64, f64::max)
        });
        let edge_rows =
            ruled || widest_segment.is_some_and(|du| !self.cylinder_cell_fits(du, widths[1]));
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
            && (ruled
                || chord_bound(&[[0.0, 0.0], [rail_du, 0.0]], &self.hessian)
                    <= self.options.chord_tolerance.mul_pow2(-2)
                    && angular_step
                        .as_ref()
                        .is_none_or(|step| chart::q(rail_du) <= step.mul_pow2(-1)));
        if !align_rails {
            rail_u.clear();
        }
        let edge_rows = edge_rows && align_rails;
        if edge_rows {
            (rail_u, rail_du) = self.row_columns(&rail_u, widest_gap, angular_step.as_ref());
        }
        let mut counts = [1usize; 2];
        for (i, count) in counts.iter_mut().enumerate() {
            if self.tiled()[i] {
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
            // Rows as far apart as columns keep a lattice, whose every triangle
            // spans one column, and leave Delaunay nothing to stagger.
            let around = self.hessian[0].to_f64_nearest() * du / std::f64::consts::TAU;
            let edge_reached = if ruled {
                dv <= around || around <= 0.0
            } else {
                !edge_rows || self.cylinder_cell_fits(du, dv)
            };
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
        let du = if align_rails {
            rail_du
        } else {
            widths[0] / counts[0] as f64
        };
        let stretch = du / (widths[1] / counts[1] as f64);
        if stretch.is_finite() && stretch > 0.0 && stretch.to_bits() != self.stretch.to_bits() {
            self.stretch = stretch;
            let nodes = self
                .cdt
                .vertices()
                .map(|vertex| {
                    let mut node = *vertex.data();
                    node.at = self.place([node.uv.x, node.uv.y]);
                    node
                })
                .collect();
            self.bulk_load(nodes)?;
        }
        let ranges: [(i64, i64); 2] = std::array::from_fn(|i| {
            if self.tiled()[i] {
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
        } else if ruled && align_rails {
            let rows: Vec<_> = (ranges[1].0..=ranges[1].1)
                .map(|j| parameter(1, j))
                .collect();
            self.encroached_lattice(&rail_u, &rows)?
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
                match self.cdt.locate(self.place(uv)) {
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
                        .map(|(i, node)| (key([node.at.x, node.at.y]), i))
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
                    if self.tiled()[i] {
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
            let points = f.vertices().map(|vertex| vertex.data().uv);
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
                    let pair = (e.from().index(), e.to().index());
                    let collar = self
                        .collar_edges
                        .contains(&[pair.0.min(pair.1), pair.0.max(pair.1)]);
                    if e.is_constraint_edge() && !collar {
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
                        && let Some(&other) = central.get(&triangle_key(
                            other.vertices().map(|vertex| vertex.data().uv),
                        ))
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
        (0..2).all(|i| !self.tiled()[i] || (0.0..1.0).contains(&uv[i]))
    }

    fn ear_point(&self, uvs: &[Uv; 3], center: Uv) -> Result<Uv, Refusal> {
        let point = self.refinement_uv(center, true);
        let point = self.refinement_uv(point, false);
        let lifted = align(point, center, self.tiled());
        // Strict orientation certifies an interior chart point, rather than a
        // rounded point on the shared polyline. Uncertain filters remain exact.
        let placed = |uv| {
            if self.collapsed_station.is_some() {
                let p = self.place(uv);
                [p.x, p.y]
            } else {
                uv
            }
        };
        let inside = (0..3).all(|i| {
            chart::positive_turn(placed(uvs[i]), placed(uvs[(i + 1) % 3]), placed(lifted))
        });
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
                let uv = triangle.vertices().map(|v| [v.data().uv.x, v.data().uv.y]);
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
                            let a = edge.from().data().uv;
                            let b = edge.to().data().uv;
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
    ) -> Result<(Vec<Triangle>, Rational, Rational, Option<Rational>), Refusal> {
        let filtered_tolerance = Interval::from_rational(&self.options.chord_tolerance);
        let spline = self.chart.spline();
        let collapsed_station = self.collapsed_station.is_some();
        let protect_band = self.protects_band();
        let cosine = if spline.is_some() {
            self.options
                .angular_deviation
                .as_ref()
                .map(|angle| {
                    if spline.is_some_and(SplineChart::has_normal_cover) {
                        crate::controls::cosine_enclosure(angle, self.dec)
                    } else {
                        crate::controls::cosine_limit(angle, self.dec)
                            .map(|q| Interval::from_rational(&q))
                    }
                })
                .transpose()?
        } else {
            None
        };
        let mut certificates = BTreeMap::<[usize; 3], (Rational, Rational)>::new();
        let mut angular_failures = BTreeMap::new();
        for _ in 0..48 {
            self.dec.check_deadline()?;
            let labels = self.labels()?;
            let mut pending = BTreeMap::new();
            let mut angular_rows = BTreeSet::new();
            let mut triangles = Vec::new();
            let mut maximum = Rational::zero();
            let mut aspect = Rational::zero();
            for f in self.cdt.inner_faces() {
                if !labels[f.index()] {
                    continue;
                }
                let nodes = f.vertices().map(|v| *v.data());
                let uvs = nodes.map(|v| [v.uv.x, v.uv.y]);
                let at = nodes.map(|n| [n.at.x, n.at.y]);
                let average =
                    |uv: [Uv; 3]| std::array::from_fn(|i| (uv[0][i] + uv[1][i] + uv[2][i]) / 3.0);
                let center = if let Some(station) = self.collapsed_station {
                    Self::unplace_collapsed(station, average(at))
                } else {
                    average(uvs)
                };
                if !self.central(center) {
                    continue;
                }
                let mut vertices = nodes.map(|n| n.mesh);
                if vertices[0] == vertices[1]
                    || vertices[1] == vertices[2]
                    || vertices[2] == vertices[0]
                {
                    if protect_band {
                        return Err(failure("periodic band facet repeats a mesh vertex"));
                    }
                    continue;
                }
                if protect_band && !chart::positive_turn(uvs[0], uvs[1], uvs[2]) {
                    return Err(failure("periodic band has no positive parameter area"));
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
                        let midpoint = if collapsed_station {
                            Self::midpoint(self.collapsed_station, a, b)
                        } else {
                            [(a.uv.x + b.uv.x) * 0.5, (a.uv.y + b.uv.y) * 0.5]
                        };
                        let midpoint = self.refinement_uv(midpoint, false);
                        pending.insert(key(midpoint), midpoint);
                        split = true;
                    }
                }
                if split {
                    continue;
                }
                // A later CDT round does not move an unchanged spline facet.
                // Key by the cyclic vertex order to retain its proved normal.
                let first = (0..3).min_by_key(|&i| vertices[i]).unwrap();
                let mut certificate_key = vertices;
                certificate_key.rotate_left(first);
                if let Some((error, ratio)) = certificates.get(&certificate_key) {
                    maximum = maximum.max(error.clone());
                    aspect = aspect.max(ratio.clone());
                    vertices = certificate_key;
                    if sense == Sense::Reversed {
                        vertices.swap(1, 2);
                    }
                    triangles.push(Triangle { vertices, face });
                    continue;
                }
                let residual = nodes.iter().map(|n| &self.errors[n.error]).max().unwrap();
                let filtered_residual = nodes
                    .iter()
                    .map(|n| self.filtered_errors[n.error].hi())
                    .fold(0.0_f64, f64::max);
                let vmax = uvs[0][1].abs().max(uvs[1][1].abs()).max(uvs[2][1].abs());
                let mut normal_bounds = None;
                let filtered = if let Some(spline) = spline {
                    let replacement = Interval::point(filtered_residual).unwrap();
                    let chord = if collapsed_station {
                        let (chord, bounds) = spline.collapsed_facet(at, self.dec)?;
                        normal_bounds = Some(bounds);
                        chord
                    } else {
                        spline.chord_with_limit(
                            &uvs,
                            (filtered_tolerance - replacement).lo(),
                            self.dec,
                        )?
                    };
                    chord + replacement
                } else if matches!(self.chart.surface, Surface::Cone(_)) {
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
                        let midpoint = if collapsed_station {
                            Self::midpoint(self.collapsed_station, nodes[longest], nodes[next])
                        } else {
                            std::array::from_fn(|i| (uvs[longest][i] + uvs[next][i]) * 0.5)
                        };
                        let midpoint = self.refinement_uv(midpoint, false);
                        pending.insert(key(midpoint), midpoint);
                        continue;
                    }
                }
                let station = self.collapsed_station;
                let split_edge = self.band.is_some() || collapsed_station;
                let refine_at = |builder: &Self, prefer_transverse: bool| {
                    if !split_edge {
                        return center;
                    }
                    if prefer_transverse {
                        let split = f
                            .adjacent_edges()
                            .into_iter()
                            .filter(|edge| !edge.is_constraint_edge())
                            .filter_map(|edge| {
                                let a = *edge.from().data();
                                let b = *edge.to().data();
                                let span = (a.uv.y - b.uv.y).abs();
                                if span == 0.0 {
                                    return None;
                                }
                                // The midpoint can encroach while a point
                                // farther along this same edge clears the trim.
                                [0.5, 0.25, 0.75, 0.125, 0.875, 0.0625, 0.9375]
                                    .into_iter()
                                    .map(|weight| {
                                        let uv = [
                                            a.uv.x + weight * (b.uv.x - a.uv.x),
                                            a.uv.y + weight * (b.uv.y - a.uv.y),
                                        ];
                                        builder.refinement_uv(uv, true)
                                    })
                                    .find(|&point| {
                                        matches!(
                                            builder.cdt.locate(builder.place(point)),
                                            PositionInTriangulation::OnEdge(found)
                                                if found.as_undirected() == edge.fix().as_undirected()
                                        ) && builder.admissible(point, &labels)
                                    })
                                    .map(|point| (span, point))
                            })
                            .max_by(|a, b| a.0.total_cmp(&b.0));
                        if let Some((_, point)) = split {
                            // Inserting on the free edge removes it from the CDT.
                            return point;
                        }
                    }
                    // A centroid can leave the long edge of a thin facet
                    // intact. Split a free edge so its chord span shrinks.
                    let point = f
                        .adjacent_edges()
                        .into_iter()
                        .filter(|edge| !edge.is_constraint_edge())
                        .map(|edge| {
                            let a = *edge.from().data();
                            let b = *edge.to().data();
                            let length = libm::hypot(a.at.x - b.at.x, a.at.y - b.at.y);
                            (length, Self::midpoint(station, a, b))
                        })
                        .max_by(|a, b| a.0.total_cmp(&b.0))
                        .map_or(center, |(_, point)| point);
                    if !protect_band {
                        return point;
                    }
                    let at = f.circumcenter();
                    let circumcenter = [at.x, at.y / builder.stretch];
                    [circumcenter, point, center]
                        .map(|uv| builder.refinement_uv(uv, true))
                        .into_iter()
                        .find(|&uv| {
                            uv[0].is_finite()
                                && uv[1].is_finite()
                                && builder.admissible(uv, &labels)
                        })
                        .unwrap_or(center)
                };
                let mut error = if self.unproved {
                    (filtered.hi() <= filtered_tolerance.lo()).then(Rational::zero)
                } else {
                    let (interpolation, replacement) = if spline.is_some()
                        || filtered.is_finite() && filtered.hi() <= filtered_tolerance.lo()
                    {
                        if !filtered.is_finite() {
                            return Err(failure("spline interpolation bound exceeds binary64"));
                        }
                        (chart::q(filtered.hi()), Rational::zero())
                    } else {
                        let hessian = self.chart.hessian(vmax, self.dec)?;
                        (chord_bound(&uvs, &hessian), residual.clone())
                    };
                    checked_triangle_bound(
                        &interpolation,
                        &replacement,
                        &self.options.chord_tolerance,
                    )
                };
                let positions = vertices.map(|i| self.positions[i]);
                if error.is_none()
                    && protect_band
                    && cosine.is_some()
                    && let Some(spline) = spline
                {
                    let bound = spline.split_chord(&uvs, positions, filtered_residual, self.dec)?;
                    if bound.is_finite() {
                        error = checked_triangle_bound(
                            &chart::q(bound.hi()),
                            &Rational::zero(),
                            &self.options.chord_tolerance,
                        );
                    }
                }
                let Some(error) = error else {
                    let center = self.refinement_uv(refine_at(self, false), true);
                    pending.insert(key(center), center);
                    continue;
                };
                let angular_limit = if spline.is_some() {
                    None
                } else {
                    self.options.angular_deviation.as_ref()
                };
                let (ratio, mut angular_reached) = if self.unproved {
                    // The preview omits the aspect report. Its chart chord
                    // and rounded-facet normal still qualify the geometry.
                    (Rational::zero(), true)
                } else {
                    let (ratio, angular_reached) = crate::controls::quality(
                        self.chart,
                        &uvs,
                        positions,
                        angular_limit,
                        self.dec,
                    )?;
                    let ratio = ratio.ok_or_else(|| failure("display triangle has zero area"))?;
                    (ratio, angular_reached)
                };
                let mut angle_violation = false;
                if let Some(spline) = spline {
                    angular_reached = if let Some(bounds) = normal_bounds {
                        spline.facet_normals_within(bounds, positions, cosine, self.dec)?
                    } else {
                        let mut chart_key = uvs.map(key);
                        chart_key.rotate_left(first);
                        let key = (certificate_key, chart_key);
                        let proof = if let Some(&proof) = angular_failures.get(&key) {
                            proof
                        } else {
                            let proof = spline.normal_proof(&uvs, positions, cosine, self.dec)?;
                            if protect_band && cosine.is_some() && proof != NormalProof::Proven {
                                angular_failures.insert(key, proof);
                            }
                            proof
                        };
                        angle_violation = proof == NormalProof::Violation;
                        proof == NormalProof::Proven
                    };
                }
                if !angular_reached {
                    let center =
                        self.refinement_uv(refine_at(self, protect_band && angle_violation), true);
                    if protect_band
                        && angle_violation
                        && let Some(spline) = spline
                        && let PositionInTriangulation::OnEdge(handle) =
                            self.cdt.locate(self.place(center))
                    {
                        let edge = self.cdt.directed_edge(handle);
                        let a = edge.from().data().uv;
                        let b = edge.to().data().uv;
                        if (a.x - b.x) == 0.0 && (a.y - b.y) != 0.0 {
                            angular_rows.insert((
                                spline.longitudinal_span(a.x),
                                key([a.y.min(b.y), a.y.max(b.y)]),
                            ));
                        }
                    }
                    pending.insert(key(center), center);
                    continue;
                }
                if spline.is_some() {
                    certificates.insert(certificate_key, (error.clone(), ratio.clone()));
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
            if !angular_rows.is_empty()
                && let Some(spline) = spline
            {
                // A split can otherwise transfer a long transverse span to
                // the next column of an aligned grid. Split the matching
                // free edges in the strip and the matching column immediately
                // beyond each end, while their endpoints still exist.
                let active_ranges: BTreeSet<_> =
                    angular_rows.iter().map(|&(_, range)| range).collect();
                let mut columns = BTreeMap::<_, Vec<_>>::new();
                for edge in self.cdt.undirected_edges() {
                    if edge.is_constraint_edge() {
                        continue;
                    }
                    let [a, b] = edge.vertices().map(|v| v.data().uv);
                    let range = key([a.y.min(b.y), a.y.max(b.y)]);
                    if (a.x - b.x) == 0.0 && active_ranges.contains(&range) {
                        let active = angular_rows.contains(&(spline.longitudinal_span(a.x), range));
                        columns.entry(range).or_default().push((
                            [a.x, a.y.midpoint(b.y)],
                            edge.fix(),
                            active,
                        ));
                    }
                }
                for columns in columns.values_mut() {
                    columns.sort_unstable_by(|a, b| a.0[0].total_cmp(&b.0[0]));
                    for (i, &(point, target, active)) in columns.iter().enumerate() {
                        let neighbor = i.checked_sub(1).is_some_and(|j| columns[j].2)
                            || columns.get(i + 1).is_some_and(|column| column.2);
                        if !active && !neighbor {
                            continue;
                        }
                        if matches!(
                            self.cdt.locate(self.place(point)),
                            PositionInTriangulation::OnEdge(found) if found.as_undirected() == target
                        ) && self.admissible(point, &labels)
                        {
                            pending.insert(key(point), point);
                        }
                    }
                }
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
            if protect_band {
                pending.retain(|_, point| self.admissible(*point, &labels));
                if pending.is_empty() {
                    return Err(failure(
                        "periodic band has no refinement point clear of its fixed trim",
                    ));
                }
            }
            let mut added = false;
            for center in pending.into_values() {
                if protect_band && !self.separated(self.place(center)) {
                    continue;
                }
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

// A cylinder with no edge limit whose trim is an intersection or spline curve. A
// circle or straight trim runs along the chart, and its grid stays as it is.
fn ruled(model: &Model, face: FaceId, chart: &Chart, options: &MeshOptions) -> bool {
    matches!(
        chart.surface,
        Surface::Cylinder(_) | Surface::EllipticCylinder(_)
    ) && options.max_edge_length.is_none()
        && model.face(face).loops.iter().any(|boundary| {
            matches!(model.loop_(*boundary), Loop::Coedges(coedges) if coedges.iter().any(|coedge| {
                let edge = model.edge(model.coedge(*coedge).edge);
                matches!(
                    model.geom().curve(edge.curve),
                    Curve::Intersection(_) | Curve::BSpline(_)
                )
            }))
        })
}

fn face_mesh(
    model: &Model,
    face_id: FaceId,
    chart: &Chart,
    boundaries: &mut Boundaries,
    plan: Plan<'_>,
    prior_chords: &BTreeSet<(usize, usize)>,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let Plan { options, level } = plan;
    let face = model.face(face_id);
    let mut segments = Vec::new();
    let mut trim_error = Rational::zero();
    let mut projected = Vec::new();
    let mut parameters = BTreeMap::new();
    let goal = options.chord_tolerance.mul_pow2(-6).to_f64_nearest();
    let mut project = |index, seed, dec: &mut Decider| -> Result<Uv, Refusal> {
        let position = boundaries.mesh.positions[index];
        if let Some(spline) = chart.spline() {
            if let Some(&uv) = parameters.get(&index) {
                return Ok(uv);
            }
            let uv = spline.project(position, seed, goal, dec)?;
            parameters.insert(index, uv);
            Ok(uv)
        } else {
            Ok(chart.project(position))
        }
    };
    for &l in &face.loops {
        match model.loop_(l) {
            Loop::Vertex(v) => projected.push(project(boundaries.vertices[v], None, dec)?),
            Loop::Coedges(cs) => {
                let mut previous = None;
                for &c in cs {
                    let c = model.coedge(c);
                    let edge = &boundaries.mesh.edges[boundaries.edges[&c.edge]];
                    if chart.spline().is_some() && level != Level::Preview {
                        trim_error = trim_error.max(edge.max_chord_error.clone());
                    }
                    let trim = chart.trim(c.edge);
                    let mut parameter = |index, seed, dec: &mut Decider| {
                        let Some(trim) = trim else {
                            return project(index, seed, dec);
                        };
                        if Some(&index) == edge.vertices.first() {
                            return Ok(trim.ends[0]);
                        }
                        if Some(&index) == edge.vertices.last() {
                            return Ok(trim.ends[1]);
                        }
                        let mut uv = project(index, seed, dec)?;
                        uv[trim.fixed.0] = trim.fixed.1;
                        Ok(uv)
                    };
                    for pair in edge.vertices.windows(2) {
                        let pair = if c.sense == Sense::Same {
                            [pair[0], pair[1]]
                        } else {
                            [pair[1], pair[0]]
                        };
                        let a = parameter(pair[0], previous, dec)?;
                        let b = align(parameter(pair[1], Some(a), dec)?, a, chart.periodic);
                        previous = Some(b);
                        projected.extend([a, b]);
                        segments.push((a, b, pair[0], pair[1]));
                    }
                }
            }
        }
    }
    let band =
        if chart.spline().is_some() && chart.periodic == [true, false] && face.loops.len() == 2 {
            Some(BandCut::open(&mut segments)?)
        } else {
            None
        };
    let sense = band.map_or(face.sense, |band| band.sense);
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
    let hessian = if chart.spline().is_some() {
        std::array::from_fn(|_| Rational::zero())
    } else {
        chart
            .hessian(bounds[1][0].abs().max(bounds[1][1].abs()), dec)?
            .into_owned()
    };
    let is_ruled = ruled(model, face_id, chart, options);
    // Spade uses one metric per face. Scale its second coordinate by the
    // physical speeds while evaluating and certifying in the original chart.
    let protect_band = band.is_some()
        && matches!(
            chart.surface,
            Surface::BSpline(SplineSupport::RollingSaddle(_))
        );
    let stretch = chart.spline().map_or(1.0, |spline| {
        // Refinement can put several samples on the same UV line. A rounded
        // metric product can turn them into a facet with no parameter area.
        if spline.has_normal_cover() {
            spline.exact_metric_stretch()
        } else {
            spline.triangulation_stretch()
        }
    });
    let mut builder = FaceBuilder {
        chart,
        options,
        band,
        collapsed_station: chart.collapsed_station(),
        ruled: is_ruled,
        stretch,
        unproved: level == Level::Preview && chart.spline().is_some() && !protect_band,
        positions: &mut boundaries.mesh.positions,
        rails: &boundaries.rails,
        cdt: Cdt::new(),
        canonical: BTreeMap::new(),
        errors: Vec::new(),
        filtered_errors: Vec::new(),
        interior_left: BTreeSet::new(),
        collar_edges: BTreeSet::new(),
        protected_trim: Vec::new(),
        separation: 0.0,
        prior_chords,
        bounds,
        filtered_hessian: hessian.each_ref().map(Interval::from_rational),
        hessian,
        dec,
    };
    builder.boundary(&segments, sense)?;
    builder.grid()?;
    let (triangles, max_chord_error, max_aspect_ratio, max_angular_deviation) =
        builder.finish(face_id, sense).map_err(|mut e| {
            e.detail = format!("face {}: {}", face_id.index(), e.detail);
            e
        })?;
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
        max_chord_error: max_chord_error.max(trim_error),
        max_aspect_ratio: (level != Level::Preview).then_some(max_aspect_ratio),
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
    check_boundary(&mesh.triangles, &BTreeSet::new())
}

pub(crate) fn check_surface(model: &Model, body: BodyId, mesh: &Mesh) -> Result<(), Refusal> {
    check_surface_faces(model, model.faces_of_body(body), mesh)
}

pub(crate) fn check_surface_faces(
    model: &Model,
    faces: impl IntoIterator<Item = FaceId>,
    mesh: &Mesh,
) -> Result<(), Refusal> {
    let mut uses = BTreeMap::<EdgeId, usize>::new();
    for face in faces {
        for boundary in &model.face(face).loops {
            if let Loop::Coedges(coedges) = model.loop_(*boundary) {
                for coedge in coedges {
                    *uses.entry(model.coedge(*coedge).edge).or_default() += 1;
                }
            }
        }
    }
    let traces: BTreeMap<_, _> = mesh
        .edges
        .iter()
        .map(|edge| (edge.edge, &edge.vertices))
        .collect();
    let mut boundary = BTreeSet::new();
    for (edge, count) in uses {
        if count != 1 {
            continue;
        }
        let Some(vertices) = traces.get(&edge).filter(|vertices| vertices.len() >= 2) else {
            return Err(failure("a surface B-rep boundary has no mesh polyline"));
        };
        for ends in vertices.windows(2) {
            boundary.insert((ends[0].min(ends[1]), ends[0].max(ends[1])));
        }
    }
    check_boundary(&mesh.triangles, &boundary)
}

fn check_boundary(
    triangles: &[Triangle],
    boundary: &BTreeSet<(usize, usize)>,
) -> Result<(), Refusal> {
    let mut uses = BTreeMap::<(usize, usize), (usize, i32)>::new();
    for t in triangles {
        for i in 0..3 {
            let (a, b) = (t.vertices[i], t.vertices[(i + 1) % 3]);
            let entry = uses.entry((a.min(b), a.max(b))).or_default();
            entry.0 += 1;
            entry.1 += if a < b { 1 } else { -1 };
        }
    }
    if uses.iter().any(|(edge, &(n, winding))| {
        if boundary.contains(edge) {
            n != 1
        } else {
            n != 2 || winding != 0
        }
    }) || boundary.iter().any(|edge| !uses.contains_key(edge))
    {
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

fn poles_clear(
    chart: &Chart,
    equations: &[Vec<RealPoly>],
    hulls: &[[Ball; 3]],
    dec: &Decider,
) -> Result<bool, Refusal> {
    for pole in [[0.0, -0.25], [0.0, 0.25]] {
        let pole = chart.evaluate(pole, dec)?;
        if hulls.iter().any(|hull| {
            hull.iter()
                .zip(&pole)
                .all(|(h, p)| p.lo() <= h.hi() && p.hi() >= h.lo())
        }) {
            return Ok(false);
        }
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
    let curves = model
        .face(face)
        .loops
        .iter()
        .flat_map(|l| match model.loop_(*l) {
            Loop::Coedges(cs) => cs
                .iter()
                .map(|c| model.edge(model.coedge(*c).edge).curve)
                .collect(),
            Loop::Vertex(_) => Vec::new(),
        });
    let mut equations = Vec::new();
    let mut hulls = Vec::new();
    for c in curves {
        let curve = model.geom().curve(c);
        if let Some(carrier) = curve.equations(model.geom()) {
            equations.push(carrier);
        } else if let Curve::BSpline(spline) = curve {
            // Positive weights enclose the entire rail in its control hull.
            // The sphere chart only needs both pole points outside that hull.
            let mut hull: Option<[Ball; 3]> = None;
            for p in spline.controls() {
                let p = [
                    dec.enclose_at(&p.x, chart::BITS)?,
                    dec.enclose_at(&p.y, chart::BITS)?,
                    dec.enclose_at(&p.z, chart::BITS)?,
                ];
                hull = Some(hull.map_or_else(
                    || p.clone(),
                    |old| std::array::from_fn(|i| old[i].hull(&p[i])),
                ));
            }
            hulls.push(hull.expect("nonempty spline control hull"));
        } else {
            return Err(failure(
                "sphere pole separation requires a finite parametric boundary solver",
            ));
        }
    }
    if equations.is_empty() && hulls.is_empty() {
        return Chart::new(surface, dec);
    }
    let sector =
        options.angular_deviation.is_some() && great_circle_sector(model, face, sphere, dec);
    if sector && let Some(corner) = convex_corner(model, face, sphere, dec) {
        let mut local = sphere.clone();
        local.axis = corner.axis;
        local.ref_dir = corner.ref_dir;
        let mut chart = Chart::new(&Surface::Sphere(local), dec)?;
        if poles_clear(&chart, &equations, &hulls, dec)? {
            let [x, y, z] = corner
                .planes
                .each_ref()
                .map(|normal| crate::evaluate::unit(normal, dec));
            chart.great_circle_sector = true;
            chart.corner_planes = Some([x?, y?, z?]);
            return Ok(chart);
        }
    }
    for choice in 0..options.max_vertices.min(
        (equations.len() + hulls.len())
            .saturating_mul(64)
            .saturating_add(1),
    ) {
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
        if poles_clear(&chart, &equations, &hulls, dec)? {
            chart.great_circle_sector = sector;
            return Ok(chart);
        }
    }
    Err(failure(
        "no sphere chart poles separated from the boundary carriers within the vertex budget",
    ))
}

fn spline_chart(surface: &Surface) -> Option<&SplineSupport> {
    match surface {
        Surface::BSpline(
            support @ (SplineSupport::BezierExtrusion(_)
            | SplineSupport::EllipticOffset(_)
            | SplineSupport::LinearFillet(_)
            | SplineSupport::ObliqueDraft(_)
            | SplineSupport::RationalPatch(_)
            | SplineSupport::RollingPipe(_)
            | SplineSupport::RollingSphere(_)
            | SplineSupport::SplineFillet(_)),
        ) => Some(support),
        Surface::BSpline(support @ SplineSupport::RollingSaddle(saddle))
            if !matches!(
                saddle.definition().law().ring,
                catallax_shapes::SaddleRing::Crossing { .. }
            ) =>
        {
            Some(support)
        }
        Surface::BSpline(support @ SplineSupport::RationalAtlas(atlas))
            if atlas.charts().len() == 1 =>
        {
            Some(support)
        }
        _ => None,
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Level {
    Preview,
    Bounded,
    Proven,
}

#[derive(Clone, Copy)]
struct Plan<'a> {
    options: &'a MeshOptions,
    level: Level,
}

pub(crate) fn shared_spline_preview(
    model: &Model,
    body: BodyId,
    options: &crate::PreviewOptions,
    dec: &mut Decider,
) -> Option<Result<crate::PreviewMesh, Refusal>> {
    let source = model;
    let expansion = source
        .chart_expansion(body)
        .filter(|cover| cover.needs_cells());
    let prepared =
        expansion.and_then(|cover| crate::rational_atlas::expanded_surface(source, body, cover));
    let (model, body) = prepared.as_ref().map_or_else(
        || (source, expansion.map_or(body, |cover| cover.body())),
        |(model, body)| (model, *body),
    );
    let faces = model.faces_of_body(body);
    let mut found_spline = false;
    for face in faces {
        let surface = model.geom().surface(model.face(face).surface);
        if matches!(surface, Surface::BSpline(_)) {
            spline_chart(surface)?;
            found_spline = true;
        }
    }
    if !found_spline {
        return None;
    }
    let options = MeshOptions {
        chord_tolerance: options.chord_tolerance.clone(),
        angular_deviation: options.angular_deviation.clone(),
        max_vertices: options.max_vertices,
        max_edge_length: None,
        max_aspect_ratio: None,
    };
    Some(
        build_junction(
            model,
            body,
            Plan {
                options: &options,
                level: Level::Preview,
            },
            None,
            None,
            dec,
        )
        .and_then(|(mut mesh, _)| {
            if let Some(cover) = expansion {
                crate::rational_atlas::remap(source, cover, &mut mesh)?;
            }
            Ok(crate::PreviewMesh::from_mesh(mesh))
        }),
    )
}

pub(crate) fn build(
    model: &Model,
    body: BodyId,
    options: &MeshOptions,
    previous: Option<(&Model, &Mesh)>,
    dec: &mut Decider,
) -> Answered<Mesh> {
    let faces = model.faces_of_body(body);
    let use_spline_charts = faces.iter().all(|&face| {
        let surface = model.geom().surface(model.face(face).surface);
        !matches!(surface, Surface::BSpline(_)) || spline_chart(surface).is_some()
    });
    let junction = if use_spline_charts {
        None
    } else {
        crate::junction::prepare(model, &faces, options, dec)?
    };
    let previous = previous
        .filter(|_| {
            !faces.is_empty()
                && faces.iter().all(|&face| {
                    !matches!(
                        model.geom().surface(model.face(face).surface),
                        Surface::BSpline(
                            SplineSupport::EllipticOffset(_) | SplineSupport::ObliqueDraft(_)
                        )
                    )
                })
        })
        .map(|(model, mesh)| crate::reuse::Previous::new(model, mesh, dec))
        .transpose()?;
    let (value, marks) = build_junction(
        model,
        body,
        Plan {
            options,
            level: Level::Proven,
        },
        junction,
        previous.as_ref(),
        dec,
    )?;
    Ok(Answer { value, marks })
}

pub(crate) fn build_position(
    model: &Model,
    body: BodyId,
    options: &MeshOptions,
    dec: &mut Decider,
) -> Result<(Mesh, catallax_numbers::MarkLog), Refusal> {
    if let Some(expansion) = model
        .chart_expansion(body)
        .filter(|cover| cover.needs_cells())
    {
        let (mut mesh, marks) = build_position(model, expansion.body(), options, dec)?;
        crate::rational_atlas::remap(model, expansion, &mut mesh)?;
        return Ok((mesh, marks));
    }
    let faces = model.faces_of_body(body);
    let junction = crate::junction::prepare_position(model, &faces, options, dec)?;
    build_junction(
        model,
        body,
        Plan {
            options,
            level: Level::Bounded,
        },
        Some(junction),
        None,
        dec,
    )
}

fn build_junction(
    model: &Model,
    body: BodyId,
    plan: Plan<'_>,
    junction: Option<crate::junction::Junction>,
    previous: Option<&crate::reuse::Previous<'_>>,
    dec: &mut Decider,
) -> Result<(Mesh, catallax_numbers::MarkLog), Refusal> {
    let Plan { options, level } = plan;
    let faces = model.faces_of_body(body);
    let fitted = |face: FaceId| junction.as_ref().is_some_and(|j| j.contains(face));
    let mut charts = faces
        .iter()
        .filter(|&&f| !fitted(f))
        .map(|&f| {
            let chart = match spline_chart(model.geom().surface(model.face(f).surface)) {
                Some(support) => Chart::from_spline(model, f, support, dec)?,
                _ => regular_chart(model, f, options, dec)?,
            };
            Ok((f, chart))
        })
        .collect::<Result<BTreeMap<_, _>, Refusal>>()?;
    constrain_lens_chord(model, &mut charts, options, dec)?;
    let mut boundaries = edges::sample(model, body, &charts, junction.as_ref(), options, dec)?;
    let boundary_limit = boundaries.mesh.positions.len();
    let mut prior_chords = BTreeSet::new();
    for face in faces {
        let local_options = charts
            .get(&face)
            .and_then(|chart| chart.chord_goal.as_ref())
            .map(|chord| MeshOptions {
                chord_tolerance: chord.clone(),
                ..options.clone()
            });
        let options = local_options.as_ref().unwrap_or(options);
        let plan = Plan { options, level };
        let start = boundaries.mesh.triangles.len();
        let reused = previous
            .map(|previous| {
                previous.append(model, face, &mut boundaries, &prior_chords, options, dec)
            })
            .transpose()?
            .unwrap_or(false);
        if !reused {
            if fitted(face) {
                junction.as_ref().expect("fitted face").triangulate(
                    face,
                    &mut boundaries,
                    options,
                    level == Level::Bounded,
                )?;
            } else if let Surface::Plane(plane) = &charts[&face].surface {
                crate::planar::build(model, face, plane, &mut boundaries, options, dec).map_err(
                    |mut e| {
                        e.detail = format!("face {}: {}", face.index(), e.detail);
                        e
                    },
                )?;
            } else {
                if let Some(center) = crate::bounded::minor_angle_center(
                    model,
                    face,
                    &charts[&face],
                    &boundaries,
                    dec,
                )? {
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
                        &charts[&face],
                        &mut boundaries,
                        plan,
                        &prior_chords,
                        dec,
                    );
                    match result {
                        Err(error)
                            if error.code == "mesh.precision-or-budget"
                                && crate::periodic::eligible(
                                    model,
                                    face,
                                    &charts[&face],
                                    options,
                                ) =>
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
    if model.body(body).kind == catallax_topology::BodyKind::Surface {
        check_surface(model, body, &boundaries.mesh)?;
    } else {
        check_closed(&boundaries.mesh)?;
    }
    compact(&mut boundaries.mesh);
    Ok((boundaries.mesh, boundaries.marks))
}

fn constrain_lens_chord(
    model: &Model,
    charts: &mut BTreeMap<FaceId, Chart>,
    options: &MeshOptions,
    dec: &Decider,
) -> Result<(), Refusal> {
    let mut limits = BTreeMap::<FaceId, Rational>::new();
    for (&face, chart) in charts.iter() {
        let Surface::BSpline(SplineSupport::RollingSaddle(saddle)) = &chart.surface else {
            continue;
        };
        let law = saddle.definition().law();
        if law.ring != catallax_shapes::SaddleRing::Overlap || law.sides != [-1, -1] {
            continue;
        }
        // In an overlap inside both cylinders, the cylinder facets and the
        // rounded strip cut into the same small volume. Sample all three
        // surfaces at a chord tied to the roll radius, even in a coarse view.
        let goal = dec.enclose_at(&law.radius, chart::BITS)?.lo().mul_pow2(-3);
        if goal >= options.chord_tolerance || goal <= Rational::zero() {
            continue;
        }
        let rails = catallax_topology::saddle_rows(model, face)
            .ok_or_else(|| failure("overlap saddle has no checked rails"))?;
        let incident = rails.into_iter().flat_map(|(edge, _)| {
            model
                .coedges_of_edge(edge)
                .into_iter()
                .filter_map(|coedge| model.owner_of_coedge(coedge).map(|(_, face)| face))
        });
        for face in incident.chain(std::iter::once(face)) {
            let limit = limits.entry(face).or_insert_with(|| goal.clone());
            if goal < *limit {
                *limit = goal.clone();
            }
        }
    }
    for (face, limit) in limits {
        if let Some(chart) = charts.get_mut(&face) {
            chart.chord_goal = Some(limit);
        }
    }
    Ok(())
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
        let is_ruled = ruled(model, face, chart, &options);
        let mut builder = FaceBuilder {
            chart,
            options: &options,
            band: None,
            collapsed_station: None,
            ruled: is_ruled,
            stretch: 1.0,
            unproved: true,
            positions: &mut boundaries.mesh.positions,
            rails: &rails,
            cdt: Cdt::new(),
            canonical: BTreeMap::new(),
            errors: Vec::new(),
            filtered_errors: Vec::new(),
            interior_left: BTreeSet::new(),
            collar_edges: BTreeSet::new(),
            protected_trim: Vec::new(),
            separation: 0.0,
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
    if boundaries.merged_vertices && model.body(body).kind == catallax_topology::BodyKind::Solid {
        check_boundary(&boundaries.mesh.triangles, &BTreeSet::new())?;
    }
    Ok(boundaries.mesh)
}
