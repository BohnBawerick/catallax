use std::{
    collections::{BTreeMap, BTreeSet},
    sync::Arc,
};

use catallax_intersect::{TracePiece, trace};
use catallax_numbers::{Ball, Decider, Interval, MarkLog, Rational, Real, Refusal};
use catallax_shapes::{
    DisplayCut, FilletCutTrace, SplineSupport, Surface, tessellation::RationalSpan,
};
use catallax_topology::{BodyId, EdgeId, Loop, Model, Sense, VertexId};

use crate::{
    EdgeMesh, Mesh, MeshOptions,
    chart::{self, BITS, Chart, Trim, Uv, align},
    failure,
};

// Keep the shared trace and parameter, not a separate exact expression graph
// for every sample. Planar filters construct exact coordinates only on demand.
#[derive(Clone)]
pub(crate) enum Position {
    Exact([Real; 3]),
    Trace(Arc<TracePiece>, Rational),
}

impl Position {
    pub(crate) fn exact(&self, dec: &mut Decider) -> Result<[Real; 3], Refusal> {
        match self {
            Self::Exact(p) => Ok(p.clone()),
            Self::Trace(piece, t) => {
                let p = piece.evaluate_exact(t, dec)?;
                Ok([p.x, p.y, p.z])
            }
        }
    }

    pub(crate) fn enclose(&self, dec: &Decider) -> Result<[Ball; 3], Refusal> {
        match self {
            Self::Exact(p) => Ok([
                dec.enclose_at(&p[0], BITS)?,
                dec.enclose_at(&p[1], BITS)?,
                dec.enclose_at(&p[2], BITS)?,
            ]),
            Self::Trace(piece, t) => piece.evaluate(&Ball::from_rational(t, BITS), dec),
        }
    }
}

/// Membership in shared topological polylines, independent of chart copies.
#[derive(Default)]
pub(crate) struct Rails {
    entries: Vec<(usize, EdgeId)>,
    end: usize,
}

impl Rails {
    fn new(edges: &[EdgeMesh]) -> Self {
        let mut entries: Vec<_> = edges
            .iter()
            .flat_map(|edge| edge.vertices.iter().map(move |&vertex| (vertex, edge.edge)))
            .collect();
        entries.sort_unstable();
        entries.dedup();
        let end = entries.last().map_or(0, |&(vertex, _)| vertex + 1);
        Self { entries, end }
    }

    pub(crate) fn contains_triangle(&self, vertices: [usize; 3]) -> bool {
        if vertices.iter().any(|&vertex| vertex >= self.end) {
            return false;
        }
        let start = self
            .entries
            .partition_point(|&(vertex, _)| vertex < vertices[0]);
        self.entries[start..]
            .iter()
            .take_while(|&&(vertex, _)| vertex == vertices[0])
            .any(|&(_, edge)| {
                self.entries.binary_search(&(vertices[1], edge)).is_ok()
                    && self.entries.binary_search(&(vertices[2], edge)).is_ok()
            })
    }
}

pub(crate) struct Boundaries {
    pub(crate) mesh: Mesh,
    pub(crate) edges: BTreeMap<EdgeId, usize>,
    pub(crate) vertices: BTreeMap<VertexId, usize>,
    pub(crate) exact_positions: BTreeMap<usize, Position>,
    pub(crate) enclosures: BTreeMap<usize, [Ball; 3]>,
    pub(crate) rails: Rails,
    pub(crate) marks: MarkLog,
}

impl Boundaries {
    pub(crate) fn vertex(
        &mut self,
        model: &Model,
        id: VertexId,
        options: &MeshOptions,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        if let Some(&v) = self.vertices.get(&id) {
            return Ok(v);
        }
        reserve_position(self.mesh.positions.len(), options)?;
        let point = model.geom().position(model.vertex(id).point);
        let xyz = [&point.x, &point.y, &point.z]
            .map(|x| dec.enclose_at(x, BITS))
            .into_iter()
            .collect::<Result<Vec<_>, _>>()?;
        let p = chart::position(&xyz.try_into().unwrap())?;
        let index = self.mesh.positions.len();
        self.mesh.positions.push(p);
        self.vertices.insert(id, index);
        Ok(index)
    }
}

fn reserve_position(count: usize, options: &MeshOptions) -> Result<(), Refusal> {
    if count >= options.max_vertices {
        return Err(Refusal::new(
            "mesh.resource-limit",
            "boundary sampling reached max_vertices",
        ));
    }
    Ok(())
}

struct Cap<'a> {
    trace: &'a FilletCutTrace,
    chart: DisplayCut,
    forward: bool,
    coefficient: Rational,
}

#[derive(Clone)]
struct EdgePoint {
    xyz: [Ball; 3],
    cap_hint: Option<Uv>,
    edge_end: Option<usize>,
}

struct RailLeaf {
    piece: Arc<TracePiece>,
    range: [Rational; 2],
    indices: [usize; 2],
    endpoints: [EdgePoint; 2],
    depth: u32,
    bound: Rational,
}

struct Sampler<'a> {
    piece: &'a Arc<TracePiece>,
    cap: Option<Cap<'a>>,
    spline: Option<RationalSpan>,
    charts: &'a [IncidentChart<'a>],
    spline_cosine: Option<Interval>,
    options: &'a MeshOptions,
    dec: &'a mut Decider,
    positions: &'a mut Vec<[f64; 3]>,
    indices: &'a mut Vec<usize>,
    error: &'a mut Rational,
    exact_positions: Option<&'a mut BTreeMap<usize, Position>>,
    rail: Option<&'a mut Vec<RailLeaf>>,
    split_leaf: bool,
}

struct IncidentChart<'a> {
    chart: &'a Chart,
    fixed_v: Option<f64>,
    trim: Option<Trim>,
}

impl Sampler<'_> {
    fn acceleration(&self, interval: &Ball) -> Result<Rational, Refusal> {
        if let Some(cap) = &self.cap {
            return Ok(cap.coefficient.clone());
        }
        if let Some(spline) = &self.spline {
            // The same chord formula uses this scaled Bernstein coefficient.
            // Positive vertex weights reparameterize the straight chord.
            return Ok(chart::q(spline.chord_coefficients()[0].hi()).mul_pow2(-1));
        }
        if let Some(jets) = self
            .piece
            .taylor_interval(interval.to_interval(), 2, self.dec)?
        {
            let sum = jets[2]
                .iter()
                .fold(Interval::point(0.0).unwrap(), |sum, x| {
                    sum + Interval::point(x.lo().abs().max(x.hi().abs())).unwrap()
                });
            if sum.is_finite() {
                return Ok(chart::q(sum.hi()));
            }
        }
        Ok(self
            .piece
            .taylor(interval, 2, self.dec)?
            .pop()
            .unwrap()
            .iter()
            .fold(Rational::zero(), |s, x| s + x.lo().abs().max(x.hi().abs())))
    }

    fn insert_middle(&mut self, mut position: [f64; 3], t: &Rational) -> Result<usize, Refusal> {
        let index = self.positions.len();
        if let Some(points) = &mut self.exact_positions {
            let point = if let Some(cap) = &self.cap {
                let cut = cap.trace.cut();
                let mut coordinates = position.map(chart::q);
                // Spell the point exactly on the cap plane for its trim
                // predicates. The chord check includes this replacement.
                let axis = (0..3).max_by_key(|&i| cut.normal[i].abs()).unwrap();
                let numerator = (0..3)
                    .filter(|&i| i != axis)
                    .fold(cut.offset.clone(), |value, i| {
                        value - &cut.normal[i] * &coordinates[i]
                    });
                coordinates[axis] = &numerator / &cut.normal[axis];
                position[axis] = coordinates[axis].to_f64_nearest();
                if !position[axis].is_finite() {
                    return Err(failure("cap sample exceeds finite binary64 coordinates"));
                }
                Position::Exact(coordinates.map(Real::from_rational))
            } else {
                Position::Trace(self.piece.clone(), t.clone())
            };
            points.insert(index, point);
        }
        self.positions.push(position);
        Ok(index)
    }

    fn enclose_point(&mut self, t: &Rational) -> Result<EdgePoint, Refusal> {
        let goal = Interval::from_rational(&self.options.chord_tolerance.mul_pow2(-4)).lo();
        let parameter = Interval::from_rational(t);
        let mut cap_hint = None;
        let filtered = if let Some(cap) = &self.cap {
            let theta = if cap.forward {
                parameter
            } else {
                Interval::point(1.0).unwrap() - parameter
            };
            cap.chart.enclose(theta, goal, self.dec)?.map(|sample| {
                cap_hint = Some(sample.parameters.map(|value| value.mid()));
                sample.position
            })
        } else if let Some(spline) = &self.spline {
            Some(spline.evaluate_box([parameter, Interval::point(0.0).unwrap()], self.dec)?)
        } else {
            self.piece
                .taylor_interval(parameter, 0, self.dec)?
                .map(|jets| jets[0])
        };
        if let Some(values) = filtered {
            let width = values.iter().fold(Interval::point(0.0).unwrap(), |sum, x| {
                sum + Interval::point(x.width()).unwrap()
            });
            if width.is_finite() && width.hi() <= goal {
                return Ok(EdgePoint {
                    xyz: values.map(|x| {
                        Ball::from_rational(&chart::q(x.lo()), BITS)
                            .hull(&Ball::from_rational(&chart::q(x.hi()), BITS))
                    }),
                    cap_hint,
                    edge_end: None,
                });
            }
        }
        Ok(EdgePoint {
            xyz: self
                .piece
                .evaluate(&Ball::from_rational(t, BITS), self.dec)?,
            cap_hint,
            edge_end: None,
        })
    }

    fn segment(
        &mut self,
        lo: Rational,
        hi: Rational,
        indices: [usize; 2],
        depth: u32,
        inherited_acceleration: Option<&Rational>,
        endpoints: [&EdgePoint; 2],
    ) -> Result<(), Refusal> {
        self.dec.check_deadline()?;
        let split_leaf = std::mem::take(&mut self.split_leaf);
        let [first, last] = indices;
        let [a, b] = endpoints;
        let interval = Ball::from_rational(&lo, BITS).hull(&Ball::from_rational(&hi, BITS));
        let width = &hi - &lo;
        let residual = chart::residual(&self.positions[first], &a.xyz)
            .max(chart::residual(&self.positions[last], &b.xyz));
        let distance =
            |acceleration: &Rational| &(&(acceleration * &width) * &width).mul_pow2(-2) + &residual;
        let acceleration = if let Some(acceleration) = inherited_acceleration
            && distance(acceleration) <= self.options.chord_tolerance.mul_pow2(-3)
        {
            acceleration.clone()
        } else {
            self.acceleration(&interval)?
        };
        let bound = distance(&acceleration);
        let mid = (&lo + &hi).mul_pow2(-1);
        let middle_point = self.enclose_point(&mid)?;
        let middle = chart::position(&middle_point.xyz)?;
        let mut small = bound <= self.options.chord_tolerance.mul_pow2(-3);
        if let Some(limit) = &self.options.max_edge_length {
            small &= crate::controls::length_squared(self.positions[first], self.positions[last])
                <= limit * limit;
        }
        for incident in self.charts {
            let chart = incident.chart;
            if let (Some(spline), Surface::BSpline(support)) = (chart.spline(), &chart.surface) {
                if small {
                    let hints = if let Some(trim) = incident.trim {
                        [a, &middle_point, b].map(|point| point.edge_end.map(|i| trim.ends[i]))
                    } else if let (Some(cap), SplineSupport::SplineFillet(fillet)) =
                        (&self.cap, support)
                        && cap.trace.fillet().shares_definition_nodes(fillet)
                    {
                        [a.cap_hint, middle_point.cap_hint, b.cap_hint]
                    } else {
                        let carrier = Chart::carrier(support);
                        [
                            self.piece
                                .display_boundary_parameters(carrier, &lo, self.dec)?,
                            self.piece
                                .display_boundary_parameters(carrier, &mid, self.dec)?,
                            self.piece
                                .display_boundary_parameters(carrier, &hi, self.dec)?,
                        ]
                    };
                    small = spline.edge_within(
                        [self.positions[first], middle, self.positions[last]],
                        hints,
                        incident
                            .trim
                            .map(|trim| trim.fixed)
                            .or_else(|| incident.fixed_v.map(|v| (1, v))),
                        self.options,
                        self.spline_cosine,
                        self.dec,
                    )?;
                }
                continue;
            }
            let a = chart.project(self.positions[first]);
            let b = align(chart.project(self.positions[last]), a, chart.periodic);
            let m = align(chart.project(middle), a, chart.periodic);
            small &= (0..2).all(|i| !chart.periodic[i] || (b[i] - a[i]).abs() <= 0.125);
            if let Some(angle) = &self.options.angular_deviation {
                let mut span = (crate::chart::q(b[0]) - crate::chart::q(a[0])).abs();
                if matches!(
                    chart.surface,
                    catallax_shapes::Surface::Sphere(_) | catallax_shapes::Surface::Torus(_)
                ) {
                    span = span + (crate::chart::q(b[1]) - crate::chart::q(a[1])).abs();
                }
                if let catallax_shapes::Surface::EllipticCylinder(s) = &chart.surface {
                    let ratio = (&s.major_sq / &s.minor_sq).sqrt();
                    span = &span * &self.dec.enclose_at(&ratio, BITS)?.hi();
                }
                if !matches!(chart.surface, catallax_shapes::Surface::Plane(_)) {
                    let span = span * Rational::from_i64(16);
                    let span =
                        if matches!(chart.surface, catallax_shapes::Surface::EllipticCylinder(_)) {
                            &span * &chart.normal_scale
                        } else {
                            span
                        };
                    small &= span <= *angle;
                }
            }
            small &= chart.edge_chord_within(
                &[a, b, m],
                &self.options.chord_tolerance.mul_pow2(-3),
                self.dec,
            )?;
        }
        if small && !split_leaf {
            if let Some(rail) = &mut self.rail {
                rail.push(RailLeaf {
                    piece: Arc::clone(self.piece),
                    range: [lo, hi],
                    indices,
                    endpoints: [a.clone(), b.clone()],
                    depth,
                    bound: bound.clone(),
                });
            }
            *self.error = self.error.clone().max(bound);
            self.indices.push(last);
            return Ok(());
        }
        if depth >= 48 || self.positions.len() >= self.options.max_vertices {
            return Err(failure(
                "edge refinement exhausted the deterministic vertex/precision limit",
            ));
        }
        let index = self.insert_middle(middle, &mid)?;
        self.segment(
            lo,
            mid.clone(),
            [first, index],
            depth + 1,
            Some(&acceleration),
            [a, &middle_point],
        )?;
        self.segment(
            mid,
            hi,
            [index, last],
            depth + 1,
            Some(&acceleration),
            [&middle_point, b],
        )
    }
}

fn plan_collar_rail(
    mut leaves: Vec<RailLeaf>,
    charts: &[IncidentChart<'_>],
    edge_cosine: Option<Interval>,
    options: &MeshOptions,
    positions: &mut Vec<[f64; 3]>,
    dec: &mut Decider,
) -> Result<(Vec<usize>, Rational), Refusal> {
    let angle = options.angular_deviation.as_ref().expect("angular rail");
    let cosine = crate::controls::cosine_enclosure(angle, dec)?;
    let mut checked = BTreeSet::new();
    loop {
        let mut split = None;
        for i in 0..leaves.len() {
            dec.charge_work(1)?;
            let j = (i + leaves.len() - 1) % leaves.len();
            let [a, b] = leaves[j].indices;
            let [start, c] = leaves[i].indices;
            debug_assert_eq!(b, start);
            let key = [a, b, c];
            if checked.contains(&key) {
                continue;
            }
            let samples = key.map(|index| positions[index]);
            let mut fits = true;
            for incident in charts {
                dec.charge_work(1)?;
                if let Some(row) = incident.fixed_v
                    && matches!(
                        incident.chart.surface,
                        Surface::BSpline(SplineSupport::RollingSaddle(_))
                    )
                {
                    let chart = incident.chart.spline().expect("saddle band chart");
                    fits &=
                        chart.rail_pair_fits(samples, row, options, cosine, edge_cosine, dec)?;
                }
            }
            if !fits {
                // Splitting a leaf reruns the ordinary edge checks for both
                // incident charts. The rest of the shared polyline stays fixed.
                let left = crate::controls::length_squared(samples[0], samples[1]);
                let right = crate::controls::length_squared(samples[1], samples[2]);
                split = Some(if left >= right { j } else { i });
                break;
            }
            checked.insert(key);
        }
        let Some(i) = split else {
            dec.charge_work(leaves.len() as u64)?;
            let mut indices = vec![leaves[0].indices[0]];
            let mut error = Rational::zero();
            for leaf in leaves {
                indices.push(leaf.indices[1]);
                error = error.max(leaf.bound);
            }
            return Ok((indices, error));
        };
        let leaf = &leaves[i];
        let mut children = Vec::new();
        let mut unused_indices = Vec::new();
        let mut error = Rational::zero();
        let spline = leaf.piece.display_spline_span(dec)?;
        Sampler {
            piece: &leaf.piece,
            cap: None,
            spline,
            charts,
            spline_cosine: edge_cosine,
            options,
            dec,
            positions,
            indices: &mut unused_indices,
            error: &mut error,
            exact_positions: None,
            rail: Some(&mut children),
            split_leaf: true,
        }
        .segment(
            leaf.range[0].clone(),
            leaf.range[1].clone(),
            leaf.indices,
            leaf.depth,
            None,
            leaf.endpoints.each_ref(),
        )?;
        dec.charge_work((leaves.len() - i - 1) as u64)?;
        leaves.splice(i..=i, children);
    }
}

// One envelope row, kept on the envelope grid so both its faces share it.
fn seeded(
    model: &Model,
    id: EdgeId,
    seed: &crate::junction::Seed,
    options: &MeshOptions,
    dec: &mut Decider,
    result: &mut Boundaries,
) -> Result<(), Refusal> {
    let edge = model.edge(id);
    let (Some(start), Some(end)) = (edge.start, edge.end) else {
        return seeded_ring(id, seed, options, result);
    };
    let first = result.vertex(model, start, options, dec)?;
    let last = result.vertex(model, end, options, dec)?;
    let mut error = seed.error.clone();
    let mut vertices = Vec::with_capacity(seed.points.len());
    vertices.push(first);
    let inner = seed.points.len() - 2;
    for i in 0..inner {
        let source = if seed.forward { i + 1 } else { inner - i };
        let point = &seed.points[source];
        reserve_position(result.mesh.positions.len(), options)?;
        let position = chart::position(point)?;
        error = error.max(&seed.error + &chart::residual(&position, point));
        if let Some(exact) = seed.exact.get(source) {
            result
                .exact_positions
                .insert(result.mesh.positions.len(), Position::Exact(exact.clone()));
        }
        vertices.push(result.mesh.positions.len());
        result.mesh.positions.push(position);
    }
    vertices.push(last);
    result.edges.insert(id, result.mesh.edges.len());
    result.mesh.edges.push(EdgeMesh {
        edge: id,
        vertices,
        max_chord_error: error,
    });
    Ok(())
}

// One closed band rail. The polyline starts at the chart's seam and returns to
// it, walking the columns in the edge's direction.
fn seeded_ring(
    id: EdgeId,
    seed: &crate::junction::Seed,
    options: &MeshOptions,
    result: &mut Boundaries,
) -> Result<(), Refusal> {
    let count = seed.points.len() - 1;
    let mut error = seed.error.clone();
    let start = result.mesh.positions.len();
    for i in 0..count {
        reserve_position(result.mesh.positions.len(), options)?;
        let point = &seed.points[i];
        let position = chart::position(point)?;
        error = error.max(&seed.error + &chart::residual(&position, point));
        if let Some(exact) = seed.exact.get(i) {
            result
                .exact_positions
                .insert(result.mesh.positions.len(), Position::Exact(exact.clone()));
        }
        result.mesh.positions.push(position);
    }
    let mut vertices = Vec::with_capacity(count + 1);
    vertices.push(start);
    for i in 1..count {
        vertices.push(start + if seed.forward { i } else { count - i });
    }
    vertices.push(start);
    result.edges.insert(id, result.mesh.edges.len());
    result.mesh.edges.push(EdgeMesh {
        edge: id,
        vertices,
        max_chord_error: error,
    });
    Ok(())
}

pub(crate) fn sample(
    model: &Model,
    body: BodyId,
    charts: &BTreeMap<catallax_topology::FaceId, Chart>,
    junction: Option<&crate::junction::Junction>,
    options: &MeshOptions,
    dec: &mut Decider,
) -> Result<Boundaries, Refusal> {
    sample_selected(
        model,
        body,
        model.edges_of_body(body),
        charts,
        junction,
        options,
        dec,
    )
}

pub(crate) fn sample_selected(
    model: &Model,
    body: BodyId,
    selected: impl IntoIterator<Item = EdgeId>,
    charts: &BTreeMap<catallax_topology::FaceId, Chart>,
    junction: Option<&crate::junction::Junction>,
    options: &MeshOptions,
    dec: &mut Decider,
) -> Result<Boundaries, Refusal> {
    let mut result = Boundaries {
        mesh: Mesh::default(),
        edges: BTreeMap::new(),
        vertices: BTreeMap::new(),
        exact_positions: BTreeMap::new(),
        enclosures: BTreeMap::new(),
        rails: Rails::default(),
        marks: model.body(body).marks.clone(),
    };
    let mut adjacent = BTreeMap::<EdgeId, Vec<IncidentChart>>::new();
    for (&face, chart) in charts {
        let band_rails = match &chart.surface {
            Surface::BSpline(SplineSupport::RollingPipe(_)) => {
                let rows = match catallax_topology::pipe_rows(model, face, dec)? {
                    Some(rows) => Some(rows),
                    None => catallax_topology::partial_pipe_rows(model, face, dec)?,
                };
                rows.map(|rows| rows.map(|(edge, _)| edge))
            }
            Surface::BSpline(SplineSupport::EllipticOffset(_)) => {
                Some(catallax_topology::SplinePrism::check(model, body, dec)?.rings())
            }
            Surface::BSpline(SplineSupport::RollingSaddle(_) | SplineSupport::RollingSphere(_)) => {
                Some(
                    catallax_topology::saddle_rows(model, face)
                        .ok_or_else(|| failure("periodic rolling band has no checked rails"))?
                        .map(|(edge, _)| edge),
                )
            }
            Surface::BSpline(SplineSupport::ObliqueDraft(_)) => {
                let cell = catallax_topology::DraftedCylinder::check(model, body, dec)?;
                if face == cell.faces()[2] {
                    Some(cell.rings())
                } else {
                    cell.hole().map(|(_, _, rings)| rings)
                }
            }
            _ => None,
        };
        for &l in &model.face(face).loops {
            match model.loop_(l) {
                Loop::Coedges(cs) => {
                    for &c in cs {
                        let edge = model.coedge(c).edge;
                        let fixed_v = band_rails.and_then(|rows| {
                            rows.iter()
                                .position(|&rail| rail == edge)
                                .map(|row| row as f64)
                        });
                        let trim = chart.trim(edge).map(|mut trim| {
                            if model.edge(edge).sense == Sense::Reversed {
                                trim.ends.reverse();
                            }
                            trim
                        });
                        adjacent.entry(edge).or_default().push(IncidentChart {
                            chart,
                            fixed_v,
                            trim,
                        });
                    }
                }
                Loop::Vertex(v) => {
                    result.vertex(model, *v, options, dec)?;
                }
            }
        }
    }
    let spline_cosine = if charts.values().any(|chart| chart.spline().is_some()) {
        // Keep the pairwise spread within the angular tolerance by testing
        // a cone of half that angle around the middle sample's direction.
        options
            .angular_deviation
            .as_ref()
            .map(|angle| {
                crate::controls::cosine_limit(&angle.mul_pow2(-1), dec)
                    .map(|cosine| Interval::from_rational(&cosine))
            })
            .transpose()?
    } else {
        None
    };
    for edge_id in selected {
        let local_options = adjacent.get(&edge_id).and_then(|charts| {
            let chord = charts
                .iter()
                .filter_map(|incident| incident.chart.chord_goal.as_ref())
                .min()?;
            Some(MeshOptions {
                chord_tolerance: chord.clone(),
                ..options.clone()
            })
        });
        let options = local_options.as_ref().unwrap_or(options);
        let spline_cosine = if adjacent.get(&edge_id).is_some_and(|charts| {
            charts.iter().any(|incident| {
                matches!(
                    incident.chart.surface,
                    Surface::BSpline(SplineSupport::RollingSaddle(_))
                )
            })
        }) {
            options
                .angular_deviation
                .as_ref()
                .map(|angle| crate::controls::cosine_enclosure(&angle.mul_pow2(-1), dec))
                .transpose()?
        } else {
            spline_cosine
        };
        if let Some(seed) = junction.and_then(|junction| junction.seed(edge_id)) {
            seeded(model, edge_id, seed, options, dec, &mut result)?;
            continue;
        }
        let edge = model.edge(edge_id);
        let (from, to) = if edge.sense == Sense::Same {
            (edge.start, edge.end)
        } else {
            (edge.end, edge.start)
        };
        let traced = trace(
            model.geom(),
            edge.curve,
            from.map(|v| model.vertex(v).point),
            to.map(|v| model.vertex(v).point),
            dec,
        )?;
        result.marks.extend(&traced.marks);
        let pieces: Vec<_> = traced
            .value
            .pieces()
            .iter()
            .cloned()
            .map(Arc::new)
            .collect();
        let initial_values = pieces
            .first()
            .map(|piece| piece.evaluate(&Ball::from_rational(&Rational::zero(), BITS), dec))
            .transpose()?;
        let planar = adjacent[&edge_id]
            .iter()
            .any(|incident| matches!(incident.chart.surface, Surface::Plane(_)));
        let first = if let Some(v) = from {
            result.vertex(model, v, options, dec)?
        } else {
            reserve_position(result.mesh.positions.len(), options)?;
            let p = chart::position(initial_values.as_ref().unwrap())?;
            let index = result.mesh.positions.len();
            result.mesh.positions.push(p);
            if planar {
                result
                    .exact_positions
                    .insert(index, Position::Trace(pieces[0].clone(), Rational::zero()));
            }
            index
        };
        let end = if let Some(v) = to {
            result.vertex(model, v, options, dec)?
        } else {
            first
        };
        let mut indices = vec![first];
        let mut error = Rational::zero();
        let plan_collar = from == to
            && !planar
            && options.angular_deviation.is_some()
            && adjacent[&edge_id].iter().any(|incident| {
                incident.fixed_v.is_some()
                    && matches!(
                        incident.chart.surface,
                        Surface::BSpline(SplineSupport::RollingSaddle(_))
                    )
            });
        let mut rail = Vec::new();
        for (i, piece) in pieces.iter().enumerate() {
            let next_values;
            let first_values = if i == 0 {
                initial_values.as_ref().unwrap()
            } else {
                next_values = piece.evaluate(&Ball::from_rational(&Rational::zero(), BITS), dec)?;
                &next_values
            };
            let last_values = piece.evaluate(&Ball::from_rational(&Rational::one(), BITS), dec)?;
            let last = if i + 1 == pieces.len() {
                end
            } else {
                reserve_position(result.mesh.positions.len(), options)?;
                let p = chart::position(&last_values)?;
                let index = result.mesh.positions.len();
                result.mesh.positions.push(p);
                if planar {
                    result
                        .exact_positions
                        .insert(index, Position::Trace(piece.clone(), Rational::one()));
                }
                index
            };
            let start = *indices.last().unwrap();
            let cap = piece
                .fillet_cut()
                .map(|(trace, forward)| {
                    let chart = trace.display_chart(dec)?;
                    let coefficient = chart.acceleration_bound(dec)?.mul_pow2(-1);
                    Ok::<_, Refusal>(Cap {
                        trace,
                        chart,
                        forward,
                        coefficient,
                    })
                })
                .transpose()?;
            let endpoint = |i: usize, xyz, whole_edge: bool| EdgePoint {
                xyz,
                cap_hint: cap.as_ref().map(|cap| {
                    let side = if cap.forward { i } else { 1 - i };
                    [cap.trace.cut().ends[side].to_f64_nearest(), side as f64]
                }),
                edge_end: whole_edge.then_some(i),
            };
            let first_point = endpoint(0, first_values.clone(), i == 0);
            let last_point = endpoint(1, last_values, i + 1 == pieces.len());
            let spline = if adjacent[&edge_id]
                .iter()
                .any(|incident| incident.chart.spline().is_some())
            {
                piece.display_spline_span(dec)?
            } else {
                None
            };
            Sampler {
                piece,
                cap,
                spline,
                charts: &adjacent[&edge_id],
                spline_cosine,
                options,
                dec,
                positions: &mut result.mesh.positions,
                indices: &mut indices,
                error: &mut error,
                exact_positions: planar.then_some(&mut result.exact_positions),
                rail: plan_collar.then_some(&mut rail),
                split_leaf: false,
            }
            .segment(
                Rational::zero(),
                Rational::one(),
                [start, last],
                0,
                None,
                [&first_point, &last_point],
            )?;
        }
        if plan_collar {
            (indices, error) = plan_collar_rail(
                rail,
                &adjacent[&edge_id],
                spline_cosine,
                options,
                &mut result.mesh.positions,
                dec,
            )?;
        }
        // Chart seams can add a trace sample closer to an endpoint than the
        // periodic display grid can distinguish. Remove only interior mesh
        // samples, retaining the exact edge and both topological endpoints.
        let mut i = 1;
        while i < indices.len() {
            let (a, b) = (indices[i - 1], indices[i]);
            let collapses = adjacent[&edge_id].iter().any(|incident| {
                let chart = incident.chart;
                chart.spline().is_none()
                    && chart::key(chart.canonical(chart.project(result.mesh.positions[a])))
                        == chart::key(chart.canonical(chart.project(result.mesh.positions[b])))
            });
            if collapses && !(a == first && b == end) {
                let displacement = (0..3).fold(Rational::zero(), |sum, k| {
                    sum + (chart::q(result.mesh.positions[a][k])
                        - chart::q(result.mesh.positions[b][k]))
                    .abs()
                });
                let bound = &error + &displacement;
                if bound <= options.chord_tolerance.mul_pow2(-3) {
                    error = bound;
                    if b == end {
                        indices.remove(i - 1);
                        i = i.saturating_sub(1).max(1);
                    } else {
                        indices.remove(i);
                    }
                    continue;
                }
            }
            i += 1;
        }
        if edge.sense == Sense::Reversed {
            indices.reverse();
        }
        result.edges.insert(edge_id, result.mesh.edges.len());
        result.mesh.edges.push(EdgeMesh {
            edge: edge_id,
            vertices: indices,
            max_chord_error: error,
        });
    }
    result.rails = Rails::new(&result.mesh.edges);
    Ok(result)
}
