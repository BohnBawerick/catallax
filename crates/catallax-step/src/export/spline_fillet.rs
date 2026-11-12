//! A spline fillet face between two patch faces, cut by end planes.
//!
//! The fillet's carrier is written with normalised knots, as a rational
//! patch is. Its rails are the carrier's own boundary rows between two of
//! its nodes, so a rail cannot drift from the face it bounds. On a trimmed
//! patch, a rail's chart curve is its exact cubic parameter curves. The
//! patch's shortened border rows are its own rows cut at the rail points.
//! Each end plane's cut is the carrier composed with a chart curve linear in
//! the second parameter, within the bound the shapes crate proves.

use super::*;
use crate::graph::Field;
use catallax_shapes::{
    BSplineSurface, CutCurve, RationalPatch, SplineBoundary, SplineFillet, SplineSupport,
};

/// How one edge of a fillet or of a trimmed patch is written.
#[derive(Clone)]
pub(super) enum Plan {
    /// A rail: row `side` of the fillet's carrier between `nodes`, shared
    /// with the trimmed patch face on `patch`.
    Rail {
        fillet: SurfaceId,
        side: usize,
        nodes: [usize; 2],
        patch: SurfaceId,
        /// Whether the edge's start is at the first node.
        forward: bool,
    },
    /// An end plane's cut of the fillet.
    Cut {
        fillet: SurfaceId,
        curve: Rc<CutCurve>,
        /// Whether the edge's start is at the first rail.
        forward: bool,
    },
    /// Part of a patch's border row, where parameter `fixed` is `at`; the
    /// other parameter runs from `ends[0]` at the edge's start to `ends[1]`.
    Segment {
        patch: SurfaceId,
        fixed: usize,
        at: Rational,
        ends: [Rational; 2],
    },
}

fn bad(detail: &str) -> Refusal {
    Refusal::new("step.spline-fillet-unproven", detail)
}

/// Written controls of a curve, with the rounding's bound on it.
struct WrittenControls {
    ids: Vec<Id>,
    points: Vec<Point>,
    weights: Vec<String>,
    error: Rational,
    bounds: Bounds,
}

/// A face's one-span patch, from a patch or a one-chart atlas.
pub(super) fn patch_of(model: &Model, surface: SurfaceId) -> Option<&RationalPatch> {
    match model.geom().surface(surface) {
        Surface::BSpline(SplineSupport::RationalPatch(patch)) => Some(patch),
        Surface::BSpline(SplineSupport::RationalAtlas(atlas)) if atlas.charts().len() == 1 => {
            Some(atlas.charts()[0].patch())
        }
        _ => None,
    }
}

fn fillet_of(model: &Model, surface: SurfaceId) -> Option<&SplineFillet> {
    match model.geom().surface(surface) {
        Surface::BSpline(SplineSupport::SplineFillet(fillet)) => Some(fillet),
        _ => None,
    }
}

/// A carrier parameter on the written, normalised knots.
fn normalised(carrier: &BSplineSurface, axis: usize, value: &Rational) -> Rational {
    let knots = carrier.knots()[axis];
    let origin = &knots[0];
    let width = knots.last().expect("clamped axis") - origin;
    (value - origin).checked_div(&width).expect("a finite axis")
}

fn same(a: &Point3, b: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
    match dec.sign(&a.minus(b).norm_squared()) {
        Outcome::Proven(_) => Ok(false),
        Outcome::ProvenCoincident(_) => Ok(true),
        Outcome::NotExact(error) => Err(Refusal::from(error)),
        Outcome::Refused(error) => Err(error),
    }
}

/// Whether `support` is `patch`'s own map on `patch`'s domain: both are one
/// rational span, so agreement at enough grid points is identity.
fn same_map(
    support: &RationalPatch,
    patch: &RationalPatch,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let [pu, pv] = support.carrier().degrees();
    let [qu, qv] = patch.carrier().degrees();
    let counts = [pu + qu + 1, pv + qv + 1];
    let domain = patch.domain();
    for i in 0..=counts[0] {
        for j in 0..=counts[1] {
            let at = [(0, i), (1, j)].map(|(axis, k)| {
                let fraction = Rational::from_ratio(k as i64, counts[axis] as i64).expect("count");
                Real::from_rational(
                    &domain[axis][0] + &(&fraction * &(&domain[axis][1] - &domain[axis][0])),
                )
            });
            if !same(
                &support.point_at(&at, dec)?,
                &patch.point_at(&at, dec)?,
                dec,
            )? {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

/// The exact coefficients of a plane's equation `n . x = c`.
fn plane_of(model: &Model, surface: SurfaceId) -> Option<([Rational; 3], Rational)> {
    let Surface::Plane(plane) = model.geom().surface(surface) else {
        return None;
    };
    let n = [
        plane.normal.x.as_rational()?,
        plane.normal.y.as_rational()?,
        plane.normal.z.as_rational()?,
    ];
    let o = [
        plane.origin.x.as_rational()?,
        plane.origin.y.as_rational()?,
        plane.origin.z.as_rational()?,
    ];
    let offset = n
        .iter()
        .zip(&o)
        .fold(Rational::zero(), |sum, (n, o)| sum + n * o);
    Some((n, offset))
}

fn boundary_side(fixed: usize, at_start: bool) -> SplineBoundary {
    match (fixed, at_start) {
        (0, true) => SplineBoundary::UMin,
        (0, false) => SplineBoundary::UMax,
        (_, true) => SplineBoundary::VMin,
        (_, false) => SplineBoundary::VMax,
    }
}

fn faces_using(model: &Model, body: BodyId) -> BTreeMap<EdgeId, Vec<FaceId>> {
    let mut uses: BTreeMap<EdgeId, Vec<FaceId>> = BTreeMap::new();
    for face in model.faces_of_body(body) {
        for lp in &model.face(face).loops {
            if let Loop::Coedges(coedges) = model.loop_(*lp) {
                for coedge in coedges {
                    uses.entry(model.coedge(*coedge).edge)
                        .or_default()
                        .push(face);
                }
            }
        }
    }
    uses
}

/// Bernstein coefficients of the same polynomial on `t = a + x (b - a)`.
fn restrict(coefficients: &[Rational], a: &Rational, b: &Rational) -> Vec<Rational> {
    // Blossom `a` taken `n - i` times and `b` taken `i` times.
    let n = coefficients.len() - 1;
    (0..=n)
        .map(|i| {
            let mut work = coefficients.to_vec();
            for step in 0..n {
                let t = if step < n - i { a } else { b };
                for j in 0..n - step {
                    work[j] = &work[j] + &(t * &(&work[j + 1] - &work[j]));
                }
            }
            work[0].clone()
        })
        .collect()
}

impl Exporter<'_> {
    /// Finds each fillet face, its rails and cuts, and its trimmed patch
    /// faces' shortened rows. Returns the trimmed faces' complete rows, to
    /// be registered after the complete patch faces.
    pub(super) fn plan_spline_fillets(
        &mut self,
        bodies: impl Iterator<Item = BodyId>,
    ) -> Result<Vec<rational_patch::Row>, Refusal> {
        let model = self.model;
        let mut full_rows = Vec::new();
        for body in bodies {
            let uses = faces_using(model, body);
            let fillets: Vec<FaceId> = model
                .faces_of_body(body)
                .into_iter()
                .filter(|face| fillet_of(model, model.face(*face).surface).is_some())
                .collect();
            for face in fillets {
                self.plan_fillet_face(face, &uses)?;
            }
            let trimmed: Vec<FaceId> = model
                .faces_of_body(body)
                .into_iter()
                .filter(|face| self.trimmed_faces.contains(face))
                .collect();
            for face in trimmed {
                full_rows.extend(self.plan_trimmed_face(face)?);
            }
        }
        Ok(full_rows)
    }

    fn plan_fillet_face(
        &mut self,
        face: FaceId,
        uses: &BTreeMap<EdgeId, Vec<FaceId>>,
    ) -> Result<(), Refusal> {
        let model = self.model;
        let surface = model.face(face).surface;
        let fillet = fillet_of(model, surface).expect("a fillet face");
        // The written chart's normal must agree with the surface's own,
        // which points to the ball.
        if !fillet.faces_ball(self.dec)? {
            return Err(bad(
                "the fillet carrier's chart normal points away from the ball",
            ));
        }
        let [lp] = model.face(face).loops.as_slice() else {
            return Err(bad("a fillet face needs one boundary loop"));
        };
        let Loop::Coedges(coedges) = model.loop_(*lp) else {
            return Err(bad("a fillet face needs edges"));
        };
        let nodes: Vec<Rational> = fillet.nodes().cloned().collect();
        let tolerance = self.options.curve_tolerance.mul_pow2(-1);
        for coedge in coedges.clone() {
            let id = model.coedge(coedge).edge;
            let edge = model.edge(id);
            let (Some(start), Some(end)) = (edge.start, edge.end) else {
                return Err(bad("a fillet edge needs two ends"));
            };
            let [start, end] = [start, end].map(|v| model.geom().position(model.vertex(v).point));
            let others: Vec<FaceId> = uses
                .get(&id)
                .into_iter()
                .flatten()
                .copied()
                .filter(|other| *other != face)
                .collect();
            let [other] = others[..] else {
                return Err(bad("a fillet edge needs one neighbour"));
            };
            let other_surface = model.face(other).surface;
            match model.geom().curve(edge.curve) {
                catallax_shapes::Curve::BSpline(curve) => {
                    let mut found = None;
                    for side in 0..2 {
                        let mut at = [None, None];
                        for (index, tau) in nodes.iter().enumerate() {
                            let point = fillet.rail_point(side, tau).expect("a node");
                            for (slot, end) in at.iter_mut().zip([&start, &end]) {
                                if slot.is_none() && same(&point, end, self.dec)? {
                                    *slot = Some(index);
                                }
                            }
                        }
                        if let [Some(a), Some(b)] = at {
                            found = Some((side, a, b));
                            break;
                        }
                    }
                    let (side, a, b) =
                        found.ok_or_else(|| bad("a fillet rail does not end at its nodes"))?;
                    let [lo, hi] = if a < b { [a, b] } else { [b, a] };
                    let rail = fillet.rail_curve(side, &nodes[lo], &nodes[hi], self.dec)?;
                    if !rail.same_definition(curve) {
                        return Err(bad("a fillet edge is not its carrier's row"));
                    }
                    let patch = patch_of(model, other_surface)
                        .ok_or_else(|| bad("a fillet rail's neighbour is not a patch face"))?;
                    if !same_map(&fillet.supports()[side], patch, self.dec)? {
                        return Err(bad("a fillet's support is not its neighbour's patch"));
                    }
                    self.trimmed_faces.insert(other);
                    self.fillet_edges.insert(
                        id,
                        Plan::Rail {
                            fillet: surface,
                            side,
                            nodes: [lo, hi],
                            patch: other_surface,
                            forward: a < b,
                        },
                    );
                }
                catallax_shapes::Curve::Intersection(_) => {
                    let (normal, offset) = plane_of(model, other_surface)
                        .ok_or_else(|| bad("a fillet cut's neighbour is not a rational plane"))?;
                    let cut = fillet
                        .plane_cut(normal, offset, self.dec)?
                        .ok_or_else(|| bad("an end plane does not cut its fillet"))?;
                    let curve = fillet.cut_curve(&cut, &tolerance, self.dec)?;
                    let first = &curve.controls[0].0;
                    let last = &curve.controls.last().expect("controls").0;
                    let forward = if same(first, &start, self.dec)? && same(last, &end, self.dec)? {
                        true
                    } else if same(first, &end, self.dec)? && same(last, &start, self.dec)? {
                        false
                    } else {
                        return Err(bad("a fillet cut does not join its edge's ends"));
                    };
                    self.fillet_edges.insert(
                        id,
                        Plan::Cut {
                            fillet: surface,
                            curve: Rc::new(curve),
                            forward,
                        },
                    );
                }
                _ => return Err(bad("a fillet edge is neither a rail nor a cut")),
            }
        }
        Ok(())
    }

    /// The rows of a trimmed patch face: each edge besides its rail is a
    /// complete border row, or part of one ending at a rail point.
    fn plan_trimmed_face(&mut self, face: FaceId) -> Result<Vec<rational_patch::Row>, Refusal> {
        let model = self.model;
        let surface = model.face(face).surface;
        let patch = patch_of(model, surface).expect("a trimmed patch face");
        let domain = patch.domain().clone();
        // Parameters of the corners and of the rail points.
        let mut known: Vec<(VertexId, [Rational; 2])> = Vec::new();
        let [lp] = model.face(face).loops.as_slice() else {
            return Err(bad("a trimmed patch face needs one boundary loop"));
        };
        let Loop::Coedges(coedges) = model.loop_(*lp) else {
            return Err(bad("a trimmed patch face needs edges"));
        };
        for coedge in coedges {
            let id = model.coedge(*coedge).edge;
            if let Some(Plan::Rail {
                fillet,
                side,
                nodes,
                ..
            }) = self.fillet_edges.get(&id)
            {
                let fillet = fillet_of(model, *fillet).expect("a fillet");
                let taus: Vec<Rational> = fillet.nodes().cloned().collect();
                let edge = model.edge(id);
                for node in nodes {
                    let parameters = fillet
                        .rail_parameters(*side, &taus[*node])
                        .expect("a rail node");
                    let point = fillet.rail_point(*side, &taus[*node]).expect("a node");
                    for vertex in [edge.start, edge.end].into_iter().flatten() {
                        let position = model.geom().position(model.vertex(vertex).point);
                        if same(&point, &position, self.dec)? {
                            known.push((vertex, parameters.clone()));
                        }
                    }
                }
            }
        }
        for u in &domain[0] {
            for v in &domain[1] {
                let corner = patch.point_at(
                    &[
                        Real::from_rational(u.clone()),
                        Real::from_rational(v.clone()),
                    ],
                    self.dec,
                )?;
                for coedge in coedges {
                    let edge = model.edge(model.coedge(*coedge).edge);
                    for vertex in [edge.start, edge.end].into_iter().flatten() {
                        let position = model.geom().position(model.vertex(vertex).point);
                        if !known.iter().any(|(k, _)| *k == vertex)
                            && same(&corner, &position, self.dec)?
                        {
                            known.push((vertex, [u.clone(), v.clone()]));
                        }
                    }
                }
            }
        }
        let parameters = |vertex: Option<VertexId>| {
            vertex.and_then(|v| known.iter().find(|(k, _)| *k == v).map(|(_, p)| p.clone()))
        };
        let mut full = Vec::new();
        for coedge in coedges {
            let use_ = model.coedge(*coedge);
            let id = use_.edge;
            if self.fillet_edges.contains_key(&id) {
                continue;
            }
            let edge = model.edge(id);
            let (Some(a), Some(b)) = (parameters(edge.start), parameters(edge.end)) else {
                return Err(bad(
                    "a trimmed patch edge does not end at corners or rail points",
                ));
            };
            let fixed = (0..2)
                .find(|&axis| a[axis] == b[axis] && domain[axis].contains(&a[axis]))
                .ok_or_else(|| bad("a trimmed patch edge is not on a border row"))?;
            let at = a[fixed].clone();
            let side = boundary_side(fixed, at == domain[fixed][0]);
            let forward = patch
                .boundary_orientation(side, model.geom().curve(edge.curve), self.dec)?
                .ok_or_else(|| bad("a trimmed patch edge is not its border row"))?;
            let free = 1 - fixed;
            let whole = [&a[free], &b[free]];
            let complete = (whole == [&domain[free][0], &domain[free][1]])
                || (whole == [&domain[free][1], &domain[free][0]]);
            if complete {
                let row_forward = whole[0] == &domain[free][0];
                if (edge.sense == Sense::Same) != (forward == row_forward) {
                    return Err(bad("a trimmed patch row reverses its parameter"));
                }
                full.push(rational_patch::Row::new(
                    id,
                    surface,
                    rational_patch::boundary_row(side),
                    row_forward,
                ));
            } else {
                self.fillet_edges.insert(
                    id,
                    Plan::Segment {
                        patch: surface,
                        fixed,
                        at,
                        ends: [a[free].clone(), b[free].clone()],
                    },
                );
            }
        }
        Ok(full)
    }

    fn knot_text(&self, value: &Rational) -> String {
        number::rational(value, self.digits()).text
    }

    /// A B-spline curve entity, rational when `weights` are given.
    fn spline_entity(
        &mut self,
        degree: usize,
        points: Vec<Id>,
        weights: Option<Vec<String>>,
        multiplicities: Vec<usize>,
        knots: Vec<String>,
    ) -> Id {
        let multiplicities = Field::List(
            multiplicities
                .into_iter()
                .map(|m| text(m.to_string()))
                .collect(),
        );
        let knots = Field::List(knots.into_iter().map(text).collect());
        match weights {
            None => self.graph.add(
                "B_SPLINE_CURVE_WITH_KNOTS",
                vec![
                    text("''"),
                    text(degree.to_string()),
                    refs(points),
                    text(".UNSPECIFIED."),
                    boolean(false),
                    boolean(false),
                    multiplicities,
                    knots,
                    text(".UNSPECIFIED."),
                ],
            ),
            Some(weights) => self.graph.add_complex(vec![
                ("BOUNDED_CURVE", vec![]),
                (
                    "B_SPLINE_CURVE",
                    vec![
                        text(degree.to_string()),
                        refs(points),
                        text(".UNSPECIFIED."),
                        boolean(false),
                        boolean(false),
                    ],
                ),
                (
                    "B_SPLINE_CURVE_WITH_KNOTS",
                    vec![multiplicities, knots, text(".UNSPECIFIED.")],
                ),
                ("CURVE", vec![]),
                ("GEOMETRIC_REPRESENTATION_ITEM", vec![]),
                (
                    "RATIONAL_B_SPLINE_CURVE",
                    vec![Field::List(weights.into_iter().map(text).collect())],
                ),
                ("REPRESENTATION_ITEM", vec![text("''")]),
            ]),
        }
    }

    /// A polynomial B-spline pcurve on a written surface. Returns it and the
    /// largest rounding of its controls, in chart units.
    fn spline_pcurve(
        &mut self,
        surface: Id,
        degree: usize,
        controls: &[[Rational; 2]],
        multiplicities: Vec<usize>,
        knots: Vec<String>,
    ) -> (Id, Rational) {
        let digits = self.digits();
        let mut rounding = Rational::zero();
        let points = controls
            .iter()
            .map(|[u, v]| {
                let [u, v] = [u, v].map(|value| number::rational(value, digits));
                rounding = rounding.clone().max(&u.error + &v.error);
                self.graph.add(
                    "CARTESIAN_POINT",
                    vec![text("''"), text(format!("({},{})", u.text, v.text))],
                )
            })
            .collect();
        let curve = self.spline_entity(degree, points, None, multiplicities, knots);
        let context = self.parameter_context();
        let representation = self.graph.add(
            "DEFINITIONAL_REPRESENTATION",
            vec![text("''"), refs([curve]), context.into()],
        );
        (
            self.graph.add(
                "PCURVE",
                vec![text("''"), surface.into(), representation.into()],
            ),
            rounding,
        )
    }

    /// Written points and weights of exact rational controls, with the
    /// rounding's bound on the curve: each point moves by its own rounding,
    /// and a relative weight change `e` moves a convex combination by at most
    /// `4 e` times the hull's size, as for written surfaces.
    fn written_controls(
        &mut self,
        controls: &[(Point3, Rational)],
    ) -> Result<WrittenControls, Refusal> {
        let digits = self.digits();
        let mut ids = Vec::with_capacity(controls.len());
        let mut points = Vec::with_capacity(controls.len());
        let mut weights = Vec::with_capacity(controls.len());
        let mut position = Rational::zero();
        let mut relative = Rational::zero();
        let mut bounds = None;
        for (point, weight) in controls {
            let xyz = [&point.x, &point.y, &point.z].map(|value| {
                value
                    .as_rational()
                    .map(|value| number::rational(&value, digits))
            });
            let [Some(x), Some(y), Some(z)] = xyz else {
                return Err(bad("a fillet curve control is not rational"));
            };
            let point = Point { xyz: [x, y, z] };
            let weight = number::rational(weight, digits);
            let lower = &weight.value - &weight.error;
            if lower <= Rational::zero() {
                return Err(number::budget("a curve weight lacks a separated decimal"));
            }
            relative = relative.max(&weight.error / &lower);
            position = position.max(point.error());
            include_bounds(&mut bounds, &point_bounds(&point));
            ids.push(self.point_entity(&point));
            points.push(point);
            weights.push(weight.text);
        }
        let bounds: Bounds = bounds.ok_or_else(unsupported)?;
        let diameter = (0..3).fold(Rational::zero(), |sum, i| {
            sum + (&bounds[i][1] - &bounds[i][0])
        });
        let error = &position + &(&relative.mul_pow2(2) * &diameter);
        Ok(WrittenControls {
            ids,
            points,
            weights,
            error,
            bounds,
        })
    }

    fn checked_ends(&self, ends: [&Point; 2], vertices: [&Vertex; 2]) -> Result<(), Refusal> {
        for (end, vertex) in ends.into_iter().zip(vertices) {
            let gap = (0..3).fold(Rational::zero(), |sum, axis| {
                sum + (&vertex.point.xyz[axis].value - &end.xyz[axis].value).abs()
            });
            if gap > self.options.curve_tolerance {
                return Err(number::budget(
                    "a rounded fillet curve end exceeds the file's sewing tolerance",
                ));
            }
        }
        Ok(())
    }

    pub(super) fn fillet_edge(
        &mut self,
        id: EdgeId,
        topology: &mut Topology,
    ) -> Result<Rc<Edge>, Refusal> {
        let plan = self.fillet_edges[&id].clone();
        let source = self.model.edge(id).clone();
        let (Some(a), Some(b)) = (source.start, source.end) else {
            return Err(unsupported());
        };
        let start = self.vertex(a, topology)?;
        let end = self.vertex(b, topology)?;
        let (curve, same, bounds) = match plan {
            Plan::Rail {
                fillet,
                side,
                nodes,
                patch,
                forward,
            } => self.rail_curve(fillet, side, nodes, patch, forward, [&start, &end])?,
            Plan::Cut {
                fillet,
                curve,
                forward,
            } => self.cut_curve(fillet, &curve, forward, [&start, &end])?,
            Plan::Segment {
                patch,
                fixed,
                at,
                ends,
            } => self.segment_curve(patch, fixed, &at, &ends, [&start, &end])?,
        };
        let edge = Rc::new(Edge {
            id: self.edge_entity(start.id, end.id, curve, same),
            bounds,
        });
        topology.edges.insert(id, edge.clone());
        Ok(edge)
    }

    fn rail_curve(
        &mut self,
        fillet_surface: SurfaceId,
        side: usize,
        nodes: [usize; 2],
        patch_surface: SurfaceId,
        forward: bool,
        vertices: [&Vertex; 2],
    ) -> Result<(Id, bool, Bounds), Refusal> {
        let model = self.model;
        let fillet = fillet_of(model, fillet_surface).expect("a fillet");
        let patch = patch_of(model, patch_surface).expect("a patch");
        let written = self.written_envelope(fillet_surface)?;
        let parameters = written.rational.as_ref().ok_or_else(unsupported)?;
        let carrier = fillet.carrier();
        let degree = carrier.degrees()[0];
        let taus: Vec<Rational> = fillet.nodes().cloned().collect();
        let knots: Vec<String> = (nodes[0]..=nodes[1])
            .map(|n| self.knot_text(&normalised(carrier, 0, &taus[n])))
            .collect();
        let mut multiplicities = vec![degree; knots.len()];
        multiplicities[0] = degree + 1;
        *multiplicities.last_mut().expect("two nodes") = degree + 1;
        let row = if side == 0 { 0 } else { parameters.rows - 1 };
        let indices: Vec<usize> = (nodes[0] * degree..=nodes[1] * degree)
            .map(|i| i * parameters.rows + row)
            .collect();
        let curve = self.spline_row(
            &written,
            &indices,
            degree,
            Field::List(multiplicities.iter().map(|m| text(m.to_string())).collect()),
            Field::List(knots.iter().cloned().map(text).collect()),
            false,
        );
        let first = &written.xyz[indices[0]];
        let last = &written.xyz[*indices.last().expect("a row")];
        let ends = if forward {
            [first, last]
        } else {
            [last, first]
        };
        self.checked_ends(ends, vertices)?;
        // On the fillet: its row, along `tau`.
        let theta = if side == 0 {
            Rational::zero()
        } else {
            Rational::one()
        };
        let along = [
            [normalised(carrier, 0, &taus[nodes[0]]), theta.clone()],
            [normalised(carrier, 0, &taus[nodes[1]]), theta],
        ];
        let (on_fillet, _) = self.spline_pcurve(
            written.id,
            1,
            &along,
            vec![2, 2],
            vec![knots[0].clone(), knots.last().expect("two nodes").clone()],
        );
        // On the patch: the rail's cubic parameter curves, span by span.
        let patch_written = self.written_envelope(patch_surface)?;
        let mut controls: Vec<[Rational; 2]> = Vec::new();
        for span in fillet.rail(side) {
            if span.tau[0] < taus[nodes[0]] || span.tau[1] > taus[nodes[1]] {
                continue;
            }
            let chart: Vec<[Rational; 2]> = (0..4)
                .map(|k| {
                    [0, 1].map(|axis| normalised(patch.carrier(), axis, &span.parameters[axis][k]))
                })
                .collect();
            if let Some(previous) = controls.last() {
                if previous != &chart[0] {
                    return Err(bad("consecutive rail spans do not share their node"));
                }
                controls.extend(chart.into_iter().skip(1));
            } else {
                controls.extend(chart);
            }
        }
        let mut cubic = vec![3; knots.len()];
        cubic[0] = 4;
        *cubic.last_mut().expect("two nodes") = 4;
        let (on_patch, rounding) =
            self.spline_pcurve(patch_written.id, 3, &controls, cubic, knots.clone());
        // A chart move of `rounding` moves the patch by at most its degree
        // sum times the written hull's size, per unit of normalised chart.
        let hull = patch_written.xyz.iter().fold(None, |mut bounds, point| {
            include_bounds(&mut bounds, &point_bounds(point));
            bounds
        });
        let hull: Bounds = hull.ok_or_else(unsupported)?;
        let size = (0..3).fold(Rational::zero(), |sum, i| sum + (&hull[i][1] - &hull[i][0]));
        let [pu, pv] = patch.carrier().degrees();
        let chart_error = &rounding * &(&size * &Rational::from_i64((pu + pv) as i64));
        self.deviation(&(&written.error + &patch_written.error) + &chart_error);
        let bounded = self.graph.add(
            "SURFACE_CURVE",
            vec![
                text("''"),
                curve.into(),
                refs([on_fillet, on_patch]),
                text(".CURVE_3D."),
            ],
        );
        let mut bounds = None;
        for index in &indices {
            include_bounds(&mut bounds, &point_bounds(&written.xyz[*index]));
        }
        Ok((bounded, forward, bounds.ok_or_else(unsupported)?))
    }

    fn cut_curve(
        &mut self,
        fillet_surface: SurfaceId,
        curve: &CutCurve,
        forward: bool,
        vertices: [&Vertex; 2],
    ) -> Result<(Id, bool, Bounds), Refusal> {
        let model = self.model;
        let fillet = fillet_of(model, fillet_surface).expect("a fillet");
        let written = self.written_envelope(fillet_surface)?;
        let WrittenControls {
            ids,
            points,
            weights,
            error: rounding,
            bounds,
        } = self.written_controls(&curve.controls)?;
        let ends = if forward {
            [&points[0], points.last().expect("controls")]
        } else {
            [points.last().expect("controls"), &points[0]]
        };
        self.checked_ends(ends, vertices)?;
        let knots: Vec<String> = curve.theta.iter().map(|t| self.knot_text(t)).collect();
        let mut multiplicities = vec![curve.degree; knots.len()];
        multiplicities[0] = curve.degree + 1;
        *multiplicities.last_mut().expect("two breakpoints") = curve.degree + 1;
        let entity = self.spline_entity(
            curve.degree,
            ids,
            Some(weights),
            multiplicities,
            knots.clone(),
        );
        // On the fillet: the chart curve, linear between breakpoints.
        let chart: Vec<[Rational; 2]> = curve
            .tau
            .iter()
            .zip(&curve.theta)
            .map(|(tau, theta)| [normalised(fillet.carrier(), 0, tau), theta.clone()])
            .collect();
        let mut linear = vec![1; knots.len()];
        linear[0] = 2;
        *linear.last_mut().expect("two breakpoints") = 2;
        let (on_fillet, _) = self.spline_pcurve(written.id, 1, &chart, linear, knots);
        self.report.approximated_curves += 1;
        self.deviation(&curve.upper + &rounding);
        let bounded = self.graph.add(
            "SURFACE_CURVE",
            vec![
                text("''"),
                entity.into(),
                refs([on_fillet]),
                text(".CURVE_3D."),
            ],
        );
        Ok((bounded, forward, bounds))
    }

    fn segment_curve(
        &mut self,
        patch_surface: SurfaceId,
        fixed: usize,
        at: &Rational,
        ends: &[Rational; 2],
        vertices: [&Vertex; 2],
    ) -> Result<(Id, bool, Bounds), Refusal> {
        let model = self.model;
        let patch = patch_of(model, patch_surface).expect("a patch");
        let free = 1 - fixed;
        let side = boundary_side(fixed, at == &patch.domain()[fixed][0]);
        let row = patch.boundary(side);
        let domain = row.domain().map(Rational::clone);
        let width = &domain[1] - &domain[0];
        let local = |value: &Rational| (value - &domain[0]).checked_div(&width).expect("width");
        // The row's exact homogeneous controls, restricted to the segment.
        let homogeneous: Vec<[Rational; 4]> = row
            .controls()
            .zip(row.weights())
            .map(|(point, weight)| {
                let w = weight.as_rational().ok_or_else(unsupported)?;
                let mut out: [Rational; 4] = Default::default();
                for (c, value) in [&point.x, &point.y, &point.z].into_iter().enumerate() {
                    out[c] = &value.as_rational().ok_or_else(unsupported)? * &w;
                }
                out[3] = w;
                Ok(out)
            })
            .collect::<Result<_, Refusal>>()?;
        let [a, b] = [local(&ends[0]), local(&ends[1])];
        let restricted: [Vec<Rational>; 4] = std::array::from_fn(|c| {
            let values: Vec<Rational> = homogeneous.iter().map(|h| h[c].clone()).collect();
            restrict(&values, &a, &b)
        });
        let controls: Vec<(Point3, Rational)> = (0..restricted[3].len())
            .map(|i| {
                let w = restricted[3][i].clone();
                let [x, y, z] =
                    [0, 1, 2].map(|c| restricted[c][i].checked_div(&w).expect("positive weight"));
                (Point3::new(x.into(), y.into(), z.into()), w)
            })
            .collect();
        let WrittenControls {
            ids,
            points,
            weights,
            error: rounding,
            bounds,
        } = self.written_controls(&controls)?;
        self.checked_ends([&points[0], points.last().expect("controls")], vertices)?;
        let degree = row.degree();
        let entity = self.spline_entity(
            degree,
            ids,
            Some(weights),
            vec![degree + 1, degree + 1],
            vec!["0.".into(), "1.".into()],
        );
        // On the patch: the row's parameter from one end to the other.
        let carrier = patch.carrier();
        let mut chart = [
            [Rational::zero(), Rational::zero()],
            [Rational::zero(), Rational::zero()],
        ];
        for (point, value) in chart.iter_mut().zip(ends) {
            point[fixed] = normalised(carrier, fixed, at);
            point[free] = normalised(carrier, free, value);
        }
        let written = self.written_envelope(patch_surface)?;
        let (on_patch, _) = self.spline_pcurve(
            written.id,
            1,
            &chart,
            vec![2, 2],
            vec!["0.".into(), "1.".into()],
        );
        self.deviation(&rounding + &written.error);
        let bounded = self.graph.add(
            "SURFACE_CURVE",
            vec![
                text("''"),
                entity.into(),
                refs([on_patch]),
                text(".CURVE_3D."),
            ],
        );
        Ok((bounded, true, bounds))
    }
}
