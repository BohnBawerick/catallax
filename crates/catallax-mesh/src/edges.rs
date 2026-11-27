use std::{
    collections::{BTreeMap, BTreeSet},
    rc::Rc,
    sync::Arc,
};

use catallax_intersect::{TracePiece, trace};
use catallax_numbers::{Ball, Decider, Interval, MarkLog, Rational, Real, Refusal};
use catallax_shapes::{
    Curve, DisplayCut, FilletCutTrace, SectionTracePiece, SplineBoundary, SplineSupport, Surface,
    tessellation::RationalSpan,
};
use catallax_topology::{BodyId, EdgeId, FaceId, Loop, Model, Sense, VertexId};

use crate::{
    EdgeMesh, Mesh, MeshOptions,
    chart::{self, BITS, Chart, RailSample, RailSpan, SplineChart, Trim, Uv, align},
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
    pub(crate) rails: Rails,
    pub(crate) marks: MarkLog,
}

impl Boundaries {
    /// Distinct shared edge samples of a face, in body-wide index order.
    pub(crate) fn face_edge_samples(&self, model: &Model, face: FaceId) -> BTreeSet<usize> {
        let mut samples = BTreeSet::new();
        for &boundary in &model.face(face).loops {
            if let Loop::Coedges(coedges) = model.loop_(boundary) {
                for &coedge in coedges {
                    if let Some(&edge) = self.edges.get(&model.coedge(coedge).edge) {
                        samples.extend(self.mesh.edges[edge].vertices.iter().copied());
                    }
                }
            }
        }
        samples
    }

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

enum EdgeEnclosure {
    Trace([Ball; 3]),
    Rail(RailSample),
}

struct EdgePoint {
    enclosure: EdgeEnclosure,
    cap_hint: Option<Uv>,
    edge_end: Option<usize>,
}

impl EdgePoint {
    fn position(&self) -> Result<[f64; 3], Refusal> {
        match &self.enclosure {
            EdgeEnclosure::Trace(xyz) => chart::position(xyz),
            EdgeEnclosure::Rail(sample) => {
                let position = sample.xyz.map(|x| x.mid() + 0.0);
                if position.iter().all(|x| x.is_finite()) {
                    Ok(position)
                } else {
                    Err(failure("rail exceeds finite binary64 coordinates"))
                }
            }
        }
    }

    fn residual(&self, position: [f64; 3], dec: &mut Decider) -> Result<Rational, Refusal> {
        match &self.enclosure {
            EdgeEnclosure::Trace(xyz) => Ok(chart::residual(&position, xyz)),
            EdgeEnclosure::Rail(sample) => {
                dec.charge_traversal(3)?;
                let error = position.into_iter().zip(sample.xyz).fold(
                    Interval::point(0.0).unwrap(),
                    |sum, (x, y)| {
                        let delta = Interval::point(x).unwrap() - y;
                        sum + Interval::new(0.0, delta.lo().abs().max(delta.hi().abs())).unwrap()
                    },
                );
                if !error.is_finite() {
                    return Err(failure("rail residual exceeds finite binary64 coordinates"));
                }
                Ok(chart::q(error.hi()))
            }
        }
    }
}

#[derive(Clone, Copy)]
enum EdgePiece<'a> {
    Trace(&'a Arc<TracePiece>),
    Rail {
        span: &'a RailSpan,
        chart: &'a SplineChart,
        row: usize,
    },
}

impl EdgePiece<'_> {
    fn trace(&self) -> Option<&Arc<TracePiece>> {
        match self {
            Self::Trace(piece) => Some(piece),
            Self::Rail { .. } => None,
        }
    }

    fn enclose_end(self, t: &Rational, dec: &mut Decider) -> Result<EdgeEnclosure, Refusal> {
        match self {
            Self::Trace(piece) => Ok(EdgeEnclosure::Trace(
                piece.evaluate(&Ball::from_rational(t, BITS), dec)?,
            )),
            Self::Rail { span, row, .. } => Ok(EdgeEnclosure::Rail(span.sample(
                row,
                Interval::from_rational(t),
                dec,
            )?)),
        }
    }

    fn admit_sample(
        self,
        point: &EdgePoint,
        position: [f64; 3],
        options: &MeshOptions,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        if let Self::Rail { chart, .. } = self
            && let EdgeEnclosure::Rail(sample) = &point.enclosure
        {
            let goal = Interval::from_rational(&options.chord_tolerance.mul_pow2(-6)).lo();
            chart.admit_rail_sample(sample, position, goal, dec)?;
        }
        Ok(())
    }
}

struct RailLeaf<'a> {
    piece: EdgePiece<'a>,
    spline: Option<Rc<RationalSpan>>,
    range: [Rational; 2],
    indices: [usize; 2],
    endpoints: [Rc<EdgePoint>; 2],
    depth: u32,
    bound: Rational,
}

struct PlanarLeaf {
    trace: Arc<TracePiece>,
    spline: Option<Rc<RationalSpan>>,
    range: [Rational; 2],
    indices: [usize; 2],
    endpoints: [Rc<EdgePoint>; 2],
    depth: u32,
}

struct Sampler<'a, 'span> {
    piece: EdgePiece<'span>,
    cap: Option<Cap<'a>>,
    spline: Option<Rc<RationalSpan>>,
    charts: &'a [IncidentChart<'a>],
    spline_cosine: Option<Interval>,
    options: &'a MeshOptions,
    dec: &'a mut Decider,
    positions: &'a mut Vec<[f64; 3]>,
    indices: &'a mut Vec<usize>,
    error: &'a mut Rational,
    exact_positions: Option<&'a mut BTreeMap<usize, Position>>,
    rail: Option<&'a mut Vec<RailLeaf<'span>>>,
    planar: Option<&'a mut BTreeMap<[usize; 2], PlanarLeaf>>,
    split_leaf: bool,
}

struct IncidentChart<'a> {
    chart: &'a Chart,
    fixed_v: Option<f64>,
    trim: Option<Trim>,
}

impl Sampler<'_, '_> {
    fn acceleration(&self, lo: &Rational, hi: &Rational) -> Result<Rational, Refusal> {
        if let Some(cap) = &self.cap {
            return Ok(cap.coefficient.clone());
        }
        let spline = match self.piece {
            EdgePiece::Rail { span, row, .. } => Some(span.curve(row)),
            EdgePiece::Trace(_) => self.spline.as_deref(),
        };
        if let Some(spline) = spline {
            // The same chord formula uses this scaled Bernstein coefficient.
            // Positive vertex weights reparameterize the straight chord.
            return Ok(chart::q(spline.chord_coefficients()[0].hi()).mul_pow2(-1));
        }
        let trace = self.piece.trace().expect("generic edge trace");
        let interval = Ball::from_rational(lo, BITS).hull(&Ball::from_rational(hi, BITS));
        if let Some(jets) = trace.taylor_interval(interval.to_interval(), 2, self.dec)? {
            let sum = jets[2]
                .iter()
                .fold(Interval::point(0.0).unwrap(), |sum, x| {
                    sum + Interval::point(x.lo().abs().max(x.hi().abs())).unwrap()
                });
            if sum.is_finite() {
                return Ok(chart::q(sum.hi()));
            }
        }
        Ok(trace
            .taylor(&interval, 2, self.dec)?
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
                Position::Trace(
                    self.piece.trace().expect("planar edge trace").clone(),
                    t.clone(),
                )
            };
            points.insert(index, point);
        }
        self.positions.push(position);
        Ok(index)
    }

    fn enclose_point(&mut self, t: &Rational) -> Result<EdgePoint, Refusal> {
        if let EdgePiece::Rail { .. } = self.piece {
            return Ok(EdgePoint {
                enclosure: self.piece.enclose_end(t, self.dec)?,
                cap_hint: None,
                edge_end: None,
            });
        }
        let trace = self.piece.trace().expect("generic edge trace");
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
            trace
                .taylor_interval(parameter, 0, self.dec)?
                .map(|jets| jets[0])
        };
        if let Some(values) = filtered {
            let width = values.iter().fold(Interval::point(0.0).unwrap(), |sum, x| {
                sum + Interval::point(x.width()).unwrap()
            });
            if width.is_finite() && width.hi() <= goal {
                return Ok(EdgePoint {
                    enclosure: EdgeEnclosure::Trace(values.map(|x| {
                        Ball::from_rational(&chart::q(x.lo()), BITS)
                            .hull(&Ball::from_rational(&chart::q(x.hi()), BITS))
                    })),
                    cap_hint,
                    edge_end: None,
                });
            }
        }
        Ok(EdgePoint {
            enclosure: EdgeEnclosure::Trace(
                trace.evaluate(&Ball::from_rational(t, BITS), self.dec)?,
            ),
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
        endpoints: [&Rc<EdgePoint>; 2],
    ) -> Result<(), Refusal> {
        self.dec.check_deadline()?;
        let split_leaf = std::mem::take(&mut self.split_leaf);
        let [first, last] = indices;
        let [a, b] = endpoints;
        let width = &hi - &lo;
        let residual = a
            .residual(self.positions[first], self.dec)?
            .max(b.residual(self.positions[last], self.dec)?);
        let distance =
            |acceleration: &Rational| &(&(acceleration * &width) * &width).mul_pow2(-2) + &residual;
        let acceleration = if let Some(acceleration) = inherited_acceleration
            && distance(acceleration) <= self.options.chord_tolerance.mul_pow2(-3)
        {
            Some(acceleration.clone())
        } else {
            match self.acceleration(&lo, &hi) {
                Ok(acceleration) => Some(acceleration),
                // Section isolation proves a regular branch. A wide power
                // enclosure can still include zero; check its subarcs.
                Err(error) if SectionTracePiece::can_subdivide_jet_refusal(&error) => None,
                Err(error) => return Err(error),
            }
        };
        let bound = acceleration.as_ref().map(distance);
        let mid = (&lo + &hi).mul_pow2(-1);
        let middle_point = Rc::new(self.enclose_point(&mid)?);
        let middle = middle_point.position()?;
        self.piece
            .admit_sample(&middle_point, middle, self.options, self.dec)?;
        let mut small = bound
            .as_ref()
            .is_some_and(|bound| *bound <= self.options.chord_tolerance.mul_pow2(-3));
        if let Some(limit) = &self.options.max_edge_length {
            small &= crate::controls::length_squared(self.positions[first], self.positions[last])
                <= limit * limit;
        }
        for incident in self.charts.iter().filter(|_| acceleration.is_some()) {
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
                    } else if let Some(piece) = self.piece.trace() {
                        let carrier = Chart::carrier(support);
                        [
                            piece.display_boundary_parameters(carrier, &lo, self.dec)?,
                            piece.display_boundary_parameters(carrier, &mid, self.dec)?,
                            piece.display_boundary_parameters(carrier, &hi, self.dec)?,
                        ]
                    } else {
                        [None; 3]
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
        if let Some(bound) = bound.filter(|_| small && !split_leaf) {
            if let Some(planar) = &mut self.planar {
                planar.insert(
                    [first.min(last), first.max(last)],
                    PlanarLeaf {
                        trace: self.piece.trace().expect("a planar trim trace").clone(),
                        spline: self.spline.clone(),
                        range: [lo.clone(), hi.clone()],
                        indices,
                        endpoints: [a.clone(), b.clone()],
                        depth,
                    },
                );
            }
            if let Some(rail) = &mut self.rail {
                rail.push(RailLeaf {
                    piece: self.piece,
                    spline: self.spline.clone(),
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
            acceleration.as_ref(),
            [a, &middle_point],
        )?;
        self.segment(
            mid,
            hi,
            [index, last],
            depth + 1,
            acceleration.as_ref(),
            [&middle_point, b],
        )
    }
}

fn plan_collar_rail(
    leaves: Vec<RailLeaf<'_>>,
    charts: &[IncidentChart<'_>],
    edge_cosine: Option<Interval>,
    options: &MeshOptions,
    positions: &mut Vec<[f64; 3]>,
    dec: &mut Decider,
) -> Result<(Vec<usize>, Rational), Refusal> {
    let angle = options.angular_deviation.as_ref().expect("angular rail");
    let cosine = crate::controls::cosine_enclosure(angle, dec)?;
    let mut pending = BTreeSet::new();
    let mut ordered = BTreeMap::new();
    for (ordinal, leaf) in leaves.into_iter().enumerate() {
        dec.charge_traversal(1)?;
        let key = (ordinal, leaf.range[0].clone());
        pending.insert(key.clone());
        ordered.insert(key, leaf);
    }
    while let Some(right_key) = pending.pop_first() {
        dec.charge_traversal(2)?;
        let (left_key, left) = ordered
            .range(..right_key.clone())
            .next_back()
            .or_else(|| ordered.last_key_value())
            .expect("nonempty rail");
        let right = &ordered[&right_key];
        let [a, b] = left.indices;
        let [start, c] = right.indices;
        debug_assert_eq!(b, start);
        let samples = [a, b, c].map(|index| positions[index]);
        let mut fits = true;
        for incident in charts {
            dec.charge_traversal(1)?;
            if let Some(row) = incident.fixed_v
                && matches!(
                    incident.chart.surface,
                    Surface::BSpline(SplineSupport::RollingSaddle(_))
                )
            {
                let chart = incident.chart.spline().expect("saddle band chart");
                fits &= chart.rail_pair_fits(samples, row, options, cosine, edge_cosine, dec)?;
            }
        }
        if fits {
            continue;
        }
        let left_length = crate::controls::length_squared(samples[0], samples[1]);
        let right_length = crate::controls::length_squared(samples[1], samples[2]);
        let split = if left_length >= right_length {
            left_key.clone()
        } else {
            right_key
        };
        // Only the new leaves and the following pair change. Process them
        // in the same cyclic order as the original scan.
        dec.charge_traversal(1)?;
        let successor = ordered
            .range((
                std::ops::Bound::Excluded(&split),
                std::ops::Bound::Unbounded,
            ))
            .next()
            .or_else(|| ordered.first_key_value())
            .expect("nonempty rail")
            .0;
        pending.insert(successor.clone());
        let leaf = ordered.remove(&split).expect("pending leaf");
        let mut children = Vec::new();
        let mut unused_indices = Vec::new();
        let mut error = Rational::zero();
        dec.charge_traversal(1)?;
        Sampler {
            piece: leaf.piece,
            cap: None,
            spline: leaf.spline,
            charts,
            spline_cosine: edge_cosine,
            options,
            dec,
            positions,
            indices: &mut unused_indices,
            error: &mut error,
            exact_positions: None,
            rail: Some(&mut children),
            planar: None,
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
        for child in children {
            dec.charge_traversal(1)?;
            let key = (split.0, child.range[0].clone());
            pending.insert(key.clone());
            ordered.insert(key, child);
        }
    }
    let mut indices = vec![ordered.first_key_value().expect("nonempty rail").1.indices[0]];
    let mut error = Rational::zero();
    for leaf in ordered.into_values() {
        dec.charge_traversal(1)?;
        indices.push(leaf.indices[1]);
        error = error.max(leaf.bound);
    }
    Ok((indices, error))
}

type CheckedRail<'a> = (&'a [RailSpan], &'a SplineChart, usize);

fn checked_rail<'a>(
    model: &Model,
    edge_id: EdgeId,
    incident: &[IncidentChart<'a>],
    dec: &mut Decider,
) -> Result<Option<CheckedRail<'a>>, Refusal> {
    let edge = model.edge(edge_id);
    // A topological seam vertex can lie inside a source knot span. Leave
    // partial arcs and rails on planar faces on their generic trace.
    if edge.start.is_some() || edge.end.is_some() {
        return Ok(None);
    }
    let Curve::BSpline(curve) = model.geom().curve(edge.curve) else {
        return Ok(None);
    };
    let mut chart = None;
    for use_ in incident {
        dec.charge_traversal(1)?;
        if matches!(use_.chart.surface, Surface::Plane(_)) {
            return Ok(None);
        }
        if use_.chart.spline().is_some() && chart.replace(use_).is_some() {
            return Ok(None);
        }
    }
    let Some(use_) = chart else { return Ok(None) };
    let Surface::BSpline(support @ SplineSupport::RollingSaddle(_)) = &use_.chart.surface else {
        return Ok(None);
    };
    let Some(row) = use_.fixed_v else {
        return Ok(None);
    };
    let row = usize::from(row > 0.0);
    let side = [SplineBoundary::VMin, SplineBoundary::VMax][row];
    let surface = Chart::carrier(support);
    if !curve.is_boundary_of(surface, side) {
        let boundary = catallax_shapes::BSplineCurve::boundary(surface, side);
        dec.charge_traversal((curve.knots().len() + curve.control_count()) as u64)?;
        if !curve.same_definition(&boundary) {
            return Ok(None);
        }
    }
    if !catallax_intersect::curve_is_closed(model.geom(), edge.curve, dec)? {
        return Ok(None);
    }
    let chart = use_.chart.spline().expect("checked saddle chart");
    Ok(chart
        .checked_rail_spans(curve, dec)?
        .map(|spans| (spans, chart, row)))
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
    sample_edges(
        model,
        body,
        model.edges_of_body(body),
        charts,
        junction,
        options,
        true,
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
    sample_edges(model, body, selected, charts, junction, options, false, dec)
}

#[allow(clippy::too_many_arguments)]
fn sample_edges(
    model: &Model,
    body: BodyId,
    selected: impl IntoIterator<Item = EdgeId>,
    charts: &BTreeMap<catallax_topology::FaceId, Chart>,
    junction: Option<&crate::junction::Junction>,
    options: &MeshOptions,
    check_crossings: bool,
    dec: &mut Decider,
) -> Result<Boundaries, Refusal> {
    let mut result = Boundaries {
        mesh: Mesh::default(),
        edges: BTreeMap::new(),
        vertices: BTreeMap::new(),
        exact_positions: BTreeMap::new(),
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
    let mut planar_leaves = BTreeMap::<EdgeId, BTreeMap<[usize; 2], PlanarLeaf>>::new();
    for edge_id in selected {
        let incident = &adjacent[&edge_id];
        let local_options = edge_options(incident, options);
        let options = local_options.as_ref().unwrap_or(options);
        let spline_cosine = edge_cosine(incident, options, spline_cosine, dec)?;
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
        let traced;
        let pieces: Vec<_> = if let Some((spans, chart, row)) =
            checked_rail(model, edge_id, &adjacent[&edge_id], dec)?
        {
            dec.charge_traversal(spans.len() as u64)?;
            spans
                .iter()
                .map(|span| EdgePiece::Rail { span, chart, row })
                .collect()
        } else {
            let answer = trace(
                model.geom(),
                edge.curve,
                from.map(|v| model.vertex(v).point),
                to.map(|v| model.vertex(v).point),
                dec,
            )?;
            result.marks.extend(&answer.marks);
            traced = answer
                .value
                .pieces()
                .iter()
                .cloned()
                .map(Arc::new)
                .collect::<Vec<_>>();
            traced.iter().map(EdgePiece::Trace).collect()
        };
        let zero = Rational::zero();
        let one = Rational::one();
        let mut initial = Some(EdgePoint {
            enclosure: pieces[0].enclose_end(&zero, dec)?,
            cap_hint: None,
            edge_end: Some(0),
        });
        let planar = adjacent[&edge_id]
            .iter()
            .any(|incident| matches!(incident.chart.surface, Surface::Plane(_)));
        let first = if let Some(v) = from {
            result.vertex(model, v, options, dec)?
        } else {
            reserve_position(result.mesh.positions.len(), options)?;
            let p = initial.as_ref().unwrap().position()?;
            let index = result.mesh.positions.len();
            result.mesh.positions.push(p);
            if planar {
                result.exact_positions.insert(
                    index,
                    Position::Trace(
                        pieces[0].trace().expect("planar edge trace").clone(),
                        zero.clone(),
                    ),
                );
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
        let needs_spline = adjacent[&edge_id]
            .iter()
            .any(|incident| incident.chart.spline().is_some());
        let mut leaves = (check_crossings
            && planar
            && !matches!(model.geom().curve(edge.curve), Curve::Line(_)))
        .then(|| planar_leaves.entry(edge_id).or_default());
        for (i, &piece) in pieces.iter().enumerate() {
            let chart_trace = piece.trace();
            let cap = chart_trace
                .map(|piece| cap(piece, dec))
                .transpose()?
                .flatten();
            let hint = |i: usize| {
                cap.as_ref().map(|cap| {
                    let side = if cap.forward { i } else { 1 - i };
                    [cap.trace.cut().ends[side].to_f64_nearest(), side as f64]
                })
            };
            let mut first_point = if let Some(point) = initial.take() {
                point
            } else {
                EdgePoint {
                    enclosure: piece.enclose_end(&zero, dec)?,
                    cap_hint: None,
                    edge_end: None,
                }
            };
            first_point.cap_hint = hint(0);
            let first_point = Rc::new(first_point);
            let last_point = Rc::new(EdgePoint {
                enclosure: piece.enclose_end(&one, dec)?,
                cap_hint: hint(1),
                edge_end: (i + 1 == pieces.len()).then_some(1),
            });
            let last = if i + 1 == pieces.len() {
                end
            } else {
                reserve_position(result.mesh.positions.len(), options)?;
                let p = last_point.position()?;
                let index = result.mesh.positions.len();
                result.mesh.positions.push(p);
                if planar {
                    result.exact_positions.insert(
                        index,
                        Position::Trace(
                            chart_trace.expect("planar edge trace").clone(),
                            one.clone(),
                        ),
                    );
                }
                index
            };
            let start = *indices.last().unwrap();
            piece.admit_sample(&first_point, result.mesh.positions[start], options, dec)?;
            piece.admit_sample(&last_point, result.mesh.positions[last], options, dec)?;
            let spline = if let Some(piece) = chart_trace.filter(|_| needs_spline) {
                piece.display_spline_span(dec)?.map(Rc::new)
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
                planar: leaves.as_deref_mut(),
                split_leaf: false,
            }
            .segment(
                zero.clone(),
                one.clone(),
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
    loop {
        let mut crossings = BTreeSet::new();
        if check_crossings {
            for (&face, chart) in charts {
                if let Surface::Plane(plane) = &chart.surface {
                    crossings.extend(crate::planar::crossing_edges(
                        model, face, plane, &result, dec,
                    )?);
                }
            }
        }
        if crossings.is_empty() {
            break;
        }
        let mut split = false;
        for (edge, pair) in crossings {
            let Some(leaves) = planar_leaves.get_mut(&edge) else {
                continue;
            };
            let Some(leaf) = leaves.remove(&pair) else {
                continue;
            };
            let boundary = &mut result.mesh.edges[result.edges[&edge]];
            let at = boundary
                .vertices
                .windows(2)
                .position(|ends| [ends[0].min(ends[1]), ends[0].max(ends[1])] == pair)
                .ok_or_else(|| failure("crossing planar chord is missing from its shared edge"))?;
            let incident = &adjacent[&edge];
            let local_options = edge_options(incident, options);
            let options = local_options.as_ref().unwrap_or(options);
            let spline_cosine = edge_cosine(incident, options, spline_cosine, dec)?;
            let mut indices = Vec::new();
            Sampler {
                piece: EdgePiece::Trace(&leaf.trace),
                cap: cap(&leaf.trace, dec)?,
                spline: leaf.spline,
                charts: incident,
                spline_cosine,
                options,
                dec,
                positions: &mut result.mesh.positions,
                indices: &mut indices,
                error: &mut boundary.max_chord_error,
                exact_positions: Some(&mut result.exact_positions),
                rail: None,
                planar: Some(leaves),
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
            indices.pop();
            if boundary.vertices[at] != leaf.indices[0] {
                indices.reverse();
            }
            boundary.vertices.splice(at + 1..at + 1, indices);
            split = true;
        }
        if !split {
            return Err(failure("crossing planar trim chords cannot be refined"));
        }
    }
    result.rails = Rails::new(&result.mesh.edges);
    Ok(result)
}

fn cap<'a>(piece: &'a TracePiece, dec: &mut Decider) -> Result<Option<Cap<'a>>, Refusal> {
    piece
        .fillet_cut()
        .map(|(trace, forward)| {
            let chart = trace.display_chart(dec)?;
            let coefficient = chart.acceleration_bound(dec)?.mul_pow2(-1);
            Ok(Cap {
                trace,
                chart,
                forward,
                coefficient,
            })
        })
        .transpose()
}

fn edge_options(incident: &[IncidentChart<'_>], options: &MeshOptions) -> Option<MeshOptions> {
    incident
        .iter()
        .filter_map(|chart| chart.chart.chord_goal.as_ref())
        .min()
        .map(|chord| MeshOptions {
            chord_tolerance: chord.clone(),
            ..options.clone()
        })
}

fn edge_cosine(
    incident: &[IncidentChart<'_>],
    options: &MeshOptions,
    cosine: Option<Interval>,
    dec: &mut Decider,
) -> Result<Option<Interval>, Refusal> {
    if incident.iter().any(|chart| {
        matches!(
            chart.chart.surface,
            Surface::BSpline(SplineSupport::RollingSaddle(_))
        )
    }) {
        options
            .angular_deviation
            .as_ref()
            .map(|angle| crate::controls::cosine_enclosure(&angle.mul_pow2(-1), dec))
            .transpose()
    } else {
        Ok(cosine)
    }
}
