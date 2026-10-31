//! Export plan of a body with rolling-envelope faces. Its cylinder and cone
//! faces end in chains of trimmed edges instead of complete rings, so their
//! seams start at a chain vertex and need no boundary split. The envelope's own
//! boundary rows carry every shared edge, which keeps their parameters equal.
//!
//! A ring fillet that only touches its mouths has no envelope, yet its torus
//! and cones end in the same chains. Its plan reuses these chart cuts.

use std::collections::BTreeMap;
use std::rc::Rc;

use catallax_intersect::{CurveForm, curve_form};
use catallax_numbers::{Ball, Decider, Outcome, Rational, Real, Refusal, Sign, StrictSign};
use catallax_shapes::{
    BoundaryRow, Point3, RollingEnvelope, SplineSupport, Surface, SurfaceId, Vector3,
};
use catallax_topology::{
    BodyId, CoedgeId, EdgeId, FaceId, Loop, Model, Sense, VertexId, envelope_rows,
};

use crate::band::{self, Association, Chain, Parameter, Plan, Ring, SourceBand};

const BITS: u32 = 128;

/// The row's neighbouring analytic chart: its u is `u + u_sign * angle` and its
/// v is `v + v_sign * h`, where h is the rail height or the meridian's angle.
pub(crate) struct Neighbour {
    pub surface: SurfaceId,
    pub u: Ball,
    pub u_sign: i64,
    pub v: Ball,
    pub v_sign: i64,
    /// Radius that turns an angular error into a distance on this chart.
    pub radius: Real,
}

pub(crate) struct EnvelopeEdge {
    pub surface: SurfaceId,
    pub row: BoundaryRow,
    pub neighbour: Neighbour,
}

fn faces(model: &Model, bodies: impl Iterator<Item = BodyId>) -> Vec<FaceId> {
    bodies
        .flat_map(|body| model.body(body).lumps.iter())
        .flat_map(|lump| model.lump(*lump).shells.iter())
        .flat_map(|shell| model.shell(*shell).faces.iter().copied())
        .collect()
}

fn envelope(model: &Model, face: FaceId) -> Option<&RollingEnvelope> {
    match model.geom().surface(model.face(face).surface) {
        Surface::BSpline(SplineSupport::RollingEnvelope(envelope)) => Some(envelope),
        _ => None,
    }
}

/// Whether any face is a rolling envelope. No other scene reaches this plan,
/// so its encoding cannot change another file's bytes.
pub(crate) fn present(model: &Model, bodies: impl Iterator<Item = BodyId>) -> bool {
    faces(model, bodies)
        .into_iter()
        .any(|face| envelope(model, face).is_some())
}

/// Whether a torus ends in a chain of arcs instead of a complete parallel. Only
/// a ring fillet touching its mouths has one, so this route changes no other
/// file's bytes.
pub(crate) fn tangent_present(model: &Model, bodies: impl Iterator<Item = BodyId>) -> bool {
    faces(model, bodies).into_iter().any(|face| {
        matches!(
            model.geom().surface(model.face(face).surface),
            Surface::Torus(_)
        ) && ring_and_chain(model, face).is_some()
    })
}

fn unsupported(detail: &'static str) -> Refusal {
    Refusal::new("step.unsupported-geometry", detail)
}

fn sign(value: &Real, dec: &mut Decider) -> Result<i64, Refusal> {
    match dec.sign(value) {
        Outcome::Proven(StrictSign::Positive) => Ok(1),
        Outcome::Proven(StrictSign::Negative) => Ok(-1),
        Outcome::ProvenCoincident(_) => Err(unsupported("a junction frame is degenerate")),
        Outcome::NotExact(n) => Err(Refusal::new(n.code, n.detail)),
        Outcome::Refused(refusal) => Err(refusal),
    }
}

// A direction's written chart angle in [0, 2 pi), for chart frame x, y.
fn chart_angle(
    direction: &Vector3,
    x: &Vector3,
    y: &Vector3,
    dec: &mut Decider,
) -> Result<Ball, Refusal> {
    let across = dec.enclose_at(&direction.dot(y), BITS)?;
    let along = dec.enclose_at(&direction.dot(x), BITS)?;
    let pi = Ball::pi_checked(BITS, dec)?;
    // The enclosed arctangent has its branch on the negative x axis, which a
    // pattern's opposite corner can hit exactly.
    if along.hi() < Rational::zero() {
        let turned = Ball::atan2(&-&across, &-&along)
            .ok_or_else(|| unsupported("a junction direction has no chart angle"))?;
        return Ok(pi + turned);
    }
    let angle = Ball::atan2(&across, &along)
        .ok_or_else(|| unsupported("a junction direction has no chart angle"))?;
    // A zero angle stays zero; only a proved negative one moves up a turn.
    Ok(if angle.hi() < Rational::zero() {
        angle + pi * Ball::from_rational(&Rational::from_i64(2), BITS)
    } else {
        angle
    })
}

fn chart_axis(surface: &Surface) -> Option<Vector3> {
    match surface {
        Surface::Cylinder(s) => Some(crate::number::unit(&s.axis)),
        Surface::Cone(s) => Some(crate::number::unit(&s.axis)),
        Surface::Torus(s) => Some(crate::number::unit(&s.axis)),
        _ => None,
    }
}

// The row's neighbouring chart data. `reference` is the surface's written x.
fn neighbour(
    model: &Model,
    surface: SurfaceId,
    envelope: &RollingEnvelope,
    row: BoundaryRow,
    reference: &Vector3,
    dec: &mut Decider,
) -> Result<Neighbour, Refusal> {
    let carrier = model.geom().surface(surface);
    let axis = chart_axis(carrier)
        .ok_or_else(|| unsupported("an envelope edge needs a curved neighbour"))?;
    let x = reference.clone();
    let y = axis.cross(&x);
    let def = envelope.definition();
    let pi = Ball::pi_checked(BITS, dec)?;
    let turn_sign = sign(&def.across().dot(&axis.cross(def.toward())), dec)?;
    let lift = sign(&def.axis().dot(&axis), dec)?;
    let enclose = |value: &Real, dec: &mut Decider| dec.enclose_at(value, BITS);
    match (row, carrier) {
        (BoundaryRow::Rail, Surface::Cylinder(cylinder)) => Ok(Neighbour {
            surface,
            u: chart_angle(def.toward(), &x, &y, dec)?,
            u_sign: turn_sign,
            v: enclose(&def.origin().minus(&cylinder.origin).dot(&axis), dec)?,
            v_sign: lift,
            radius: def.body().clone(),
        }),
        (BoundaryRow::Mouth, Surface::Cone(cone)) => {
            // The near point is half a turn from the mouth's `toward` direction.
            let near = def.toward().scaled(&Real::from_i64(-1));
            let center = def.origin().plus(&def.toward().scaled(def.offset()));
            Ok(Neighbour {
                surface,
                u: chart_angle(&near, &x, &y, dec)?,
                u_sign: turn_sign,
                v: enclose(&center.minus(&cone.apex).dot(&axis), dec)?,
                v_sign: 0,
                radius: def.mouth().clone(),
            })
        }
        (BoundaryRow::Negative | BoundaryRow::Positive, Surface::Torus(_)) => {
            let side = if row == BoundaryRow::Negative {
                Sign::Negative
            } else {
                Sign::Positive
            };
            let [foot, _] = def.corners(side);
            let radial = foot.minus(def.origin());
            let radial = radial.minus(&def.axis().scaled(&radial.dot(def.axis())));
            Ok(Neighbour {
                surface,
                u: chart_angle(&radial, &x, &y, dec)?,
                u_sign: 0,
                v: pi,
                v_sign: lift,
                radius: def.radius().clone(),
            })
        }
        _ => Err(unsupported("an envelope row has the wrong neighbour")),
    }
}

// A face's complete ring and its closed multi-edge boundary.
fn ring_and_chain(model: &Model, face: FaceId) -> Option<(CoedgeId, &[CoedgeId])> {
    let [a, b] = model.face(face).loops.as_slice() else {
        return None;
    };
    let (Loop::Coedges(a), Loop::Coedges(b)) = (model.loop_(*a), model.loop_(*b)) else {
        return None;
    };
    match (a.as_slice(), b.as_slice()) {
        ([ring], chain) if chain.len() > 1 => Some((*ring, chain)),
        (chain, [ring]) if chain.len() > 1 => Some((*ring, chain)),
        _ => None,
    }
}

fn chain_uses(model: &Model, chain: &[CoedgeId]) -> Vec<(EdgeId, bool)> {
    chain
        .iter()
        .map(|c| {
            let c = model.coedge(*c);
            (c.edge, c.sense == Sense::Same)
        })
        .collect()
}

fn chain_start(model: &Model, use_: &(EdgeId, bool)) -> Option<VertexId> {
    let edge = model.edge(use_.0);
    if use_.1 { edge.start } else { edge.end }
}

// Turns a chain round to walk the other way, keeping the seam vertex first.
fn turn_round(model: &Model, uses: &mut [(EdgeId, bool)], vertex: VertexId) {
    uses.reverse();
    for use_ in uses.iter_mut() {
        use_.1 = !use_.1;
    }
    let first = uses
        .iter()
        .position(|use_| chain_start(model, use_) == Some(vertex))
        .expect("closed chain");
    uses.rotate_left(first);
}

// The closed multi-edge boundary of a junction-side cylinder or cone, with its
// complete opposite ring.
fn chain_band(
    model: &Model,
    face: FaceId,
    dec: &mut Decider,
) -> Result<Option<SourceBand>, Refusal> {
    let definition = model.face(face);
    let Some((ring, chain)) = ring_and_chain(model, face) else {
        return Ok(None);
    };
    let surface = model.geom().surface(definition.surface);
    let Some((origin, axis, _, slant)) = band::carrier_frame(surface) else {
        return Ok(None);
    };
    let Some(ring) = band::source_ring(model, face, model.coedge(ring).edge, dec)? else {
        return Ok(None);
    };
    let mut uses = chain_uses(model, chain);
    // Every chain carries circular arcs of the carrier's parallel, whose
    // normal fixes the walking direction about the axis.
    let mut circle = None;
    for &(edge, same) in &uses {
        if let CurveForm::Circle {
            center,
            normal,
            radius_sq,
            ..
        } = curve_form(model.geom(), model.edge(edge).curve, dec)?
        {
            let reversed = (model.edge(edge).sense == Sense::Same) != same;
            circle = Some((center, normal, radius_sq, reversed));
            break;
        }
    }
    let Some((center, normal, radius_sq, reversed)) = circle else {
        return Ok(None);
    };
    let along = match dec.sign(&normal.dot(&axis)) {
        Outcome::Proven(sign) => (sign == StrictSign::Positive) != reversed,
        Outcome::Refused(refusal) => return Err(refusal),
        _ => return Ok(None),
    };
    let Some(vertex) = chain_start(model, &uses[0]) else {
        return Ok(None);
    };
    let point = model.geom().position(model.vertex(vertex).point);
    let v = point.minus(origin).dot(&axis);
    let at_low = match dec.compare(&v, &ring.v) {
        Outcome::Proven(catallax_numbers::StrictOrder::Less) => true,
        Outcome::Proven(catallax_numbers::StrictOrder::Greater) => false,
        Outcome::Refused(refusal) => return Err(refusal),
        _ => return Ok(None),
    };
    if along != at_low {
        turn_round(model, &mut uses, vertex);
    }
    let radius = crate::number::sqrt_real(&radius_sq);
    let reference = point.minus(&center).scaled(&(&Real::one() / &radius));
    Ok(Some(SourceBand {
        face,
        surface: definition.surface,
        axis,
        reference,
        v_scale: slant,
        rings: vec![ring],
        apex: None,
        chain: Some((
            Chain {
                uses,
                vertex,
                at_low,
            },
            v,
            radius,
        )),
    }))
}

// A torus between its innermost parallel, a complete ring, and a quarter-turn
// parallel split into arcs. Two parallels bound two bands: the low boundary
// runs with the chart's u, read on the face's side, only on the quarter
// between them.
fn torus_chain_band(
    model: &Model,
    face: FaceId,
    dec: &mut Decider,
) -> Result<Option<SourceBand>, Refusal> {
    let definition = model.face(face);
    let Surface::Torus(torus) = model.geom().surface(definition.surface) else {
        return Ok(None);
    };
    let Some((ring, chain)) = ring_and_chain(model, face) else {
        return Ok(None);
    };
    let mut uses = chain_uses(model, chain);
    let axis = crate::number::unit(&torus.axis);
    let ring_edge = model.coedge(ring).edge;
    let first = uses[0].0;
    let (Some(ring), Some(parallel)) = (
        band::torus_ring(model, torus, &axis, ring_edge, dec)?,
        band::torus_ring(model, torus, &axis, first, dec)?,
    ) else {
        return Ok(None);
    };
    let curve = model.edge(first).curve;
    if ring.turns != Some(2)
        || parallel.turns == Some(2)
        || uses.iter().any(|use_| model.edge(use_.0).curve != curve)
    {
        return Ok(None);
    }
    let Some(vertex) = chain_start(model, &uses[0]) else {
        return Ok(None);
    };
    let at_low = parallel.turns < ring.turns;
    let along = model.edge(first).sense == Sense::Same;
    let with_chart = parallel.increasing == (along == uses[0].1);
    let low_with_chart = if at_low {
        with_chart
    } else {
        let along = model.edge(ring_edge).sense == Sense::Same;
        ring.increasing == (along == band::ring_use(model, face, ring_edge)?)
    };
    if low_with_chart != (definition.sense == Sense::Same) {
        return Ok(None);
    }
    if with_chart != at_low {
        turn_round(model, &mut uses, vertex);
    }
    let point = model.geom().position(model.vertex(vertex).point);
    let reference = point
        .minus(&parallel.center)
        .scaled(&(&Real::one() / &parallel.radius));
    Ok(Some(SourceBand {
        face,
        surface: definition.surface,
        axis,
        reference,
        v_scale: torus.minor.clone(),
        rings: vec![ring],
        apex: None,
        chain: Some((
            Chain {
                uses,
                vertex,
                at_low,
            },
            Real::zero(),
            parallel.radius,
        )),
    }))
}

fn add_arc(plan: &mut Plan, edge: EdgeId, radius: &Real, association: Association) {
    let ring = plan.arcs.entry(edge).or_insert_with(|| {
        Rc::new(Ring {
            origin: Point3::origin(),
            radius: radius.clone(),
            associations: Vec::with_capacity(2),
        })
    });
    Rc::get_mut(ring)
        .expect("arcs stay unshared while planning")
        .associations
        .push(association);
}

// Exact lines of a face's circular arcs. Each arc's circle takes the chart's
// reference, so its parameter is the chart angle.
fn chain_arcs(
    model: &Model,
    face: FaceId,
    reference: &Vector3,
    plan: &mut Plan,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    for lp in &model.face(face).loops {
        let Loop::Coedges(coedges) = model.loop_(*lp) else {
            continue;
        };
        if coedges.len() < 2 {
            continue;
        }
        for coedge in coedges {
            let id = model.coedge(*coedge).edge;
            if plan.envelope_edges.contains_key(&id) {
                continue;
            }
            let Some(ring) = band::source_ring(model, face, id, dec)? else {
                return Err(unsupported("a junction arc is not a parallel of its face"));
            };
            let surface = model.geom().surface(model.face(face).surface);
            let (.., v_scale) = band::carrier_frame(surface).expect("cylinder or cone");
            plan.curves.insert(
                ring.curve,
                CurveForm::Circle {
                    center: ring.center.clone(),
                    normal: ring.normal.clone(),
                    radius_sq: ring.radius_sq.clone(),
                    ref_dir: reference.clone(),
                },
            );
            add_arc(
                plan,
                id,
                &ring.radius,
                Association {
                    surface: model.face(face).surface,
                    v: Parameter::Exact(ring.v.clone()),
                    v_scale,
                    increasing: ring.increasing,
                },
            );
        }
    }
    Ok(())
}

fn torus_arcs(
    model: &Model,
    face: FaceId,
    reference: &Vector3,
    plan: &mut Plan,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let surface = model.face(face).surface;
    let Surface::Torus(torus) = model.geom().surface(surface) else {
        return Err(unsupported("a junction piece is not a torus"));
    };
    let axis = crate::number::unit(&torus.axis);
    let mut arcs = Vec::new();
    for lp in &model.face(face).loops {
        let Loop::Coedges(coedges) = model.loop_(*lp) else {
            return Err(unsupported("a torus piece has an edge boundary"));
        };
        if coedges.len() > 1 {
            arcs.extend(coedges.iter().copied());
        }
    }
    for coedge in arcs {
        let id = model.coedge(coedge).edge;
        if plan.envelope_edges.contains_key(&id) {
            continue;
        }
        let curve = model.edge(id).curve;
        let CurveForm::Circle {
            center,
            normal,
            radius_sq,
            ..
        } = curve_form(model.geom(), curve, dec)?
        else {
            return Err(unsupported("a torus piece arc is not a circle"));
        };
        let radius = crate::number::sqrt_real(&radius_sq);
        // The cylinder contact is the torus's innermost parallel; the plane
        // contact lies a tube radius along the axis from the core circle.
        let inner = match dec.compare(&radius, &torus.major) {
            Outcome::Proven(catallax_numbers::StrictOrder::Less) => true,
            Outcome::ProvenCoincident(_) => false,
            Outcome::Refused(refusal) => return Err(refusal),
            _ => return Err(unsupported("a torus piece arc is not a contact parallel")),
        };
        let turns = if inner {
            2
        } else if sign(&center.minus(&torus.center).dot(&axis), dec)? > 0 {
            1
        } else {
            3
        };
        let increasing = sign(&normal.dot(&axis), dec)? > 0;
        plan.curves.entry(curve).or_insert(CurveForm::Circle {
            center,
            normal,
            radius_sq,
            ref_dir: reference.clone(),
        });
        add_arc(
            plan,
            id,
            &radius,
            Association {
                surface,
                v: Parameter::QuarterTurns(turns),
                v_scale: torus.minor.clone(),
                increasing,
            },
        );
    }
    Ok(())
}

/// The plan of a ring fillet touching its mouths, with every other body in the
/// write. A scene without such a ring keeps the default plan. A face that has
/// no seam plan beside one refuses, since dropping the plan for the whole write
/// would lose the seams of every periodic face. That includes a later Boolean:
/// sewing merges the contact caps back, and no file can carry the pinched plane.
pub(crate) fn tangent_plan(
    model: &Model,
    bodies: impl Iterator<Item = BodyId> + Clone,
    dec: &mut Decider,
) -> Result<Plan, Refusal> {
    let mut sources = Vec::new();
    let mut tori = Vec::new();
    let mut chains = Vec::new();
    let mut unplanned = false;
    for face in faces(model, bodies) {
        let definition = model.face(face);
        let band = match model.geom().surface(definition.surface) {
            Surface::Plane(_) => continue,
            Surface::Sphere(_) | Surface::Torus(_) if definition.loops.is_empty() => continue,
            Surface::Torus(_) => match torus_chain_band(model, face, dec)? {
                Some(band) => {
                    tori.push(face);
                    Some(band)
                }
                None => band::torus_band(model, face, dec)?,
            },
            Surface::Cylinder(_) | Surface::Cone(_) => match band::source_band(model, face, dec)? {
                Some(band) => Some(band),
                None => {
                    chains.push(face);
                    chain_band(model, face, dec)?
                }
            },
            _ => None,
        };
        match band {
            Some(band) => sources.push(band),
            None => unplanned = true,
        }
    }
    if tori.is_empty() {
        return Ok(Plan::default());
    }
    if unplanned {
        return Err(unsupported(
            "a face written beside a tangent fillet ring has no seam plan",
        ));
    }
    let mut plan = band::assemble(model, sources);
    for face in chains {
        let reference = plan.references[&model.face(face).surface].clone();
        chain_arcs(model, face, &reference, &mut plan, dec)?;
    }
    for face in tori {
        let reference = plan.references[&model.face(face).surface].clone();
        torus_arcs(model, face, &reference, &mut plan, dec)?;
    }
    Ok(plan)
}

pub(crate) fn plan(
    model: &Model,
    bodies: impl Iterator<Item = BodyId> + Clone,
    dec: &mut Decider,
) -> Result<Plan, Refusal> {
    let faces = faces(model, bodies);
    let mut owners = BTreeMap::<EdgeId, Vec<FaceId>>::new();
    for &face in &faces {
        for lp in &model.face(face).loops {
            if let Loop::Coedges(coedges) = model.loop_(*lp) {
                for coedge in coedges {
                    owners
                        .entry(model.coedge(*coedge).edge)
                        .or_default()
                        .push(face);
                }
            }
        }
    }
    let mut sources = Vec::new();
    let mut tori = Vec::new();
    let mut envelopes = Vec::new();
    let mut chains = Vec::new();
    for &face in &faces {
        let definition = model.face(face);
        match model.geom().surface(definition.surface) {
            Surface::Plane(_) => {}
            Surface::Torus(_) => tori.push(face),
            Surface::Cylinder(_) | Surface::Cone(_) => {
                let band = match band::source_band(model, face, dec)? {
                    Some(band) => Some(band),
                    None => {
                        chains.push(face);
                        chain_band(model, face, dec)?
                    }
                };
                sources.push(band.ok_or_else(|| {
                    unsupported("a junction neighbour has no cut chart boundary")
                })?);
            }
            Surface::BSpline(SplineSupport::RollingEnvelope(envelope)) => {
                envelopes.push((face, envelope))
            }
            _ => return Err(unsupported("no STEP form beside a rolling envelope")),
        }
    }
    let mut plan = band::assemble(model, sources);
    // Every torus piece shares the chain cylinder's cut, so an arc's circle has
    // one parameter on both faces and each piece stays in one chart branch.
    let mut reference = None;
    for &face in &chains {
        let surface = model.face(face).surface;
        if matches!(model.geom().surface(surface), Surface::Cylinder(_)) {
            reference = plan.references.get(&surface).cloned();
        }
    }
    let reference = reference.ok_or_else(|| unsupported("a junction has no cylinder chain"))?;
    for &face in &tori {
        plan.references
            .insert(model.face(face).surface, reference.clone());
    }
    for (face, envelope) in envelopes {
        for (edge, row, _) in envelope_rows(model, face, dec)? {
            let other = owners[&edge]
                .iter()
                .copied()
                .find(|owner| *owner != face)
                .ok_or_else(|| unsupported("an envelope edge has no neighbouring face"))?;
            let surface = model.face(other).surface;
            let written = plan
                .references
                .get(&surface)
                .ok_or_else(|| unsupported("an envelope neighbour has no chart"))?
                .clone();
            let neighbour = neighbour(model, surface, envelope, row, &written, dec)?;
            plan.envelope_edges.insert(
                edge,
                Rc::new(EnvelopeEdge {
                    surface: model.face(face).surface,
                    row,
                    neighbour,
                }),
            );
        }
    }
    for &face in &chains {
        let surface = model.face(face).surface;
        let reference = plan.references[&surface].clone();
        chain_arcs(model, face, &reference, &mut plan, dec)?;
    }
    for &face in &tori {
        torus_arcs(model, face, &reference, &mut plan, dec)?;
    }
    // A meridian on the torus cut takes the end of the turn its piece lies on.
    let turn = Ball::pi_checked(BITS, dec)? * Ball::from_rational(&Rational::from_i64(2), BITS);
    for &face in &tori {
        let meridians: Vec<EdgeId> = owners
            .iter()
            .filter(|(edge, owner)| owner.contains(&face) && plan.envelope_edges.contains_key(edge))
            .map(|(edge, _)| *edge)
            .collect();
        let [a, c] = meridians.as_slice() else {
            return Err(unsupported("a torus piece has two meridians"));
        };
        let angles = [a, c].map(|edge| plan.envelope_edges[edge].neighbour.u.clone());
        for (edge, other) in [(a, &angles[1]), (c, &angles[0])] {
            let own = &plan.envelope_edges[edge].neighbour.u;
            if (other - own).lo() > Ball::pi_checked(BITS, dec)?.hi() {
                let entry = plan.envelope_edges.get_mut(edge).expect("meridian");
                let entry = Rc::get_mut(entry).expect("edges stay unshared while planning");
                entry.neighbour.u = &entry.neighbour.u + &turn;
            }
        }
    }
    Ok(plan)
}
