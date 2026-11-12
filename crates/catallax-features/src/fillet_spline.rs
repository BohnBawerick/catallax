//! A fillet on an edge between two spline faces (M11).
//!
//! The edge runs along a whole border row of each face's one-span patch, each
//! face is its whole patch with no other trims, and the edge is convex. Each
//! end vertex meets one more face, a rational plane that holds the side
//! faces' border edges there, as a loft's caps do. The ball rolls inside the
//! material.
//!
//! Both patches are extended past their faces along the edge, so the fitted
//! face reaches past both planes, with a node of `tau` where each rail meets
//! its face's border row. Each plane's proved cut joins its two rail points.
//! The side faces are trimmed along the rails, their border edges end at
//! the rail points, and each plane gains its cut where the edge's vertex
//! was. The blend face is named by the edge, so a radius edit keeps it.

use super::*;
use crate::{extrude::face_on, fillet::UNSUPPORTED_SPLINE};
use catallax_numbers::GEOMETRIC_GAP_UNPROVEN;
use catallax_shapes::{
    BSplineCurve, IntersectionCurve, Meet, MeetOf, Position, RailBorder, RationalPatch,
    SPLINE_FILLET_RAIL_DOMAIN, SplineFillet, SplineFilletDefinition, SurfaceId,
};
use catallax_topology::RationalPatchFace;

/// A face's one-span patch, from a patch or a one-chart atlas.
fn patch_of(model: &Model, face: FaceId) -> Option<RationalPatch> {
    match model.geom().surface(model.face(face).surface) {
        Surface::BSpline(SplineSupport::RationalPatch(patch)) => Some(patch.clone()),
        Surface::BSpline(SplineSupport::RationalAtlas(atlas)) if atlas.charts().len() == 1 => {
            Some(atlas.charts()[0].patch().clone())
        }
        _ => None,
    }
}

/// Whether the selection is one edge between two patch faces, which this
/// path rounds.
pub(super) fn applies(model: &Model, body: BodyId, selected: &[EdgeId]) -> bool {
    let [edge] = selected else {
        return false;
    };
    match uses(model, body).get(edge).map(Vec::as_slice) {
        Some([(a, _), (b, _)]) => patch_of(model, *a).is_some() && patch_of(model, *b).is_some(),
        _ => false,
    }
}

/// A side face along the edge, and the border row of its patch the edge
/// runs along.
struct Side {
    face: FaceId,
    /// The input edge's coedge sense in this face.
    coedge: Sense,
    patch: RationalPatch,
    /// The parameter held on the row, and its value there.
    fixed: usize,
    at: Rational,
    /// The other parameter's values at the edge's start and end.
    ends: [Rational; 2],
}

impl Side {
    fn along(&self) -> usize {
        1 - self.fixed
    }

    fn point(
        &self,
        fixed: &Rational,
        along: &Rational,
        dec: &mut Decider,
    ) -> Result<Point3, Refusal> {
        let mut at = [Real::zero(), Real::zero()];
        at[self.fixed] = Real::from_rational(fixed.clone());
        at[self.along()] = Real::from_rational(along.clone());
        self.patch.point_at(&at, dec)
    }

    /// The ball's side of the patch, against the face's outward normal.
    fn ball(&self, model: &Model) -> Sign {
        if model.face(self.face).sense == Sense::Same {
            Sign::Negative
        } else {
            Sign::Positive
        }
    }

    /// Binary64 directions at the edge's middle: along the edge as the
    /// patch's own parameter increases, into the face, and out of the
    /// material. They only choose which construction to try; every choice
    /// is then proved.
    fn directions(&self, model: &Model, dec: &mut Decider) -> Result<[[f64; 3]; 3], Refusal> {
        let mut at = [Real::zero(), Real::zero()];
        at[self.fixed] = Real::from_rational(self.at.clone());
        let middle = (&self.ends[0] + &self.ends[1]).mul_pow2(-1);
        at[self.along()] = Real::from_rational(middle);
        let partial = self
            .patch
            .derivatives_at(&at, dec)?
            .map(|d| [d.x.to_f64(), d.y.to_f64(), d.z.to_f64()]);
        let normal = cross(partial[0], partial[1]);
        let outward = if model.face(self.face).sense == Sense::Same {
            normal
        } else {
            normal.map(|x| -x)
        };
        let into = if self.at == self.patch.domain()[self.fixed][0] {
            partial[self.fixed]
        } else {
            partial[self.fixed].map(|x| -x)
        };
        Ok([partial[self.along()], into, outward])
    }
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Whether a face is its whole patch: one loop of its four complete border
/// rows, checked from the edges. Then a rail inside the patch's rectangle is
/// inside the face, whatever else the face's neighbours are.
fn whole(model: &Model, face: FaceId, dec: &mut Decider) -> Result<(), Refusal> {
    let not_whole = || {
        Refusal::new(
            UNSUPPORTED_SPLINE,
            "a spline fillet's side face must be its whole patch, with no other trims",
        )
    };
    match RationalPatchFace::of_chart_face(model, face, dec) {
        Ok(Some(_)) => Ok(()),
        Ok(None) => Err(not_whole()),
        Err(refusal) if refusal.code == GEOMETRIC_GAP_UNPROVEN => Err(not_whole()),
        Err(refusal) => Err(refusal),
    }
}

/// The side for a face whose patch has the edge as a whole border row.
fn side(
    model: &Model,
    face: FaceId,
    coedge: Sense,
    ends: [&Point3; 2],
    dec: &mut Decider,
) -> Result<Side, Refusal> {
    let patch =
        patch_of(model, face).ok_or_else(|| bad("a spline fillet needs one-span patch faces"))?;
    whole(model, face, dec)?;
    let domain = patch.domain().clone();
    for fixed in 0..2 {
        for at in &domain[fixed] {
            let mut side = Side {
                face,
                coedge,
                patch: patch.clone(),
                fixed,
                at: at.clone(),
                ends: domain[1 - fixed].clone(),
            };
            let corners = [
                side.point(at, &domain[1 - fixed][0], dec)?,
                side.point(at, &domain[1 - fixed][1], dec)?,
            ];
            if same(&corners[0], ends[0], dec)? && same(&corners[1], ends[1], dec)? {
                return Ok(side);
            }
            if same(&corners[1], ends[0], dec)? && same(&corners[0], ends[1], dec)? {
                side.ends.reverse();
                return Ok(side);
            }
        }
    }
    Err(bad(
        "a spline fillet's edge must be a whole border row of both patches",
    ))
}

/// The exact coefficients of a plane face's equation `n . x = c`.
fn plane_of(model: &Model, face: FaceId) -> Option<([Rational; 3], Rational)> {
    let Surface::Plane(plane) = model.geom().surface(model.face(face).surface) else {
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

/// The one face besides the sides that meets a vertex.
fn cap_at(model: &Model, uses: &Uses, vertex: VertexId, sides: [FaceId; 2]) -> Option<FaceId> {
    let mut faces = BTreeSet::new();
    for (edge, users) in uses {
        let e = model.edge(*edge);
        if e.start == Some(vertex) || e.end == Some(vertex) {
            faces.extend(users.iter().map(|(face, _)| *face));
        }
    }
    let mut others = faces.into_iter().filter(|face| !sides.contains(face));
    match (others.next(), others.next()) {
        (Some(face), None) => Some(face),
        _ => None,
    }
}

/// The one edge of `face` besides `edge` that ends at `vertex`.
fn neighbour(model: &Model, face: FaceId, edge: EdgeId, vertex: VertexId) -> Option<EdgeId> {
    let mut found = None;
    for lp in &model.face(face).loops {
        let Loop::Coedges(cs) = model.loop_(*lp) else {
            continue;
        };
        for c in cs {
            let other = model.coedge(*c).edge;
            let e = model.edge(other);
            if other != edge
                && (e.start == Some(vertex) || e.end == Some(vertex))
                && found.replace(other).is_some()
            {
                return None;
            }
        }
    }
    found
}

/// Whether a point lies on an edge's line strictly between its ends.
fn inside_line(
    model: &Model,
    edge: EdgeId,
    point: &Point3,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let e = model.edge(edge);
    let (Some(start), Some(end), Curve::Line(line)) = (e.start, e.end, model.geom().curve(e.curve))
    else {
        return Ok(false);
    };
    let [start, end] = [start, end].map(|v| model.geom().position(model.vertex(v).point));
    let off = point.minus(&line.origin).cross(&line.dir);
    let zero = |value: &Real, dec: &mut Decider| Ok::<_, Refusal>(sign(value, dec)? == Sign::Zero);
    if !(zero(&off.x, dec)? && zero(&off.y, dec)? && zero(&off.z, dec)?) {
        return Ok(false);
    }
    let between = point.minus(&start).dot(&end.minus(point));
    Ok(sign(&between, dec)? == Sign::Positive)
}

/// Rounds one convex edge between two spline faces.
pub(super) fn build(
    model: &mut Model,
    body: BodyId,
    edge: EdgeId,
    radius: &Real,
    step: &StepId,
    dec: &mut Decider,
) -> Result<(BodyId, Vec<FilletDefinition>), Refusal> {
    let radius = radius
        .as_rational()
        .ok_or_else(|| bad("a spline fillet needs a rational radius"))?;
    let uses = uses(model, body);
    let input = model.edge(edge).clone();
    let (Some(start), Some(end)) = (input.start, input.end) else {
        return Err(bad("a spline fillet needs an edge with two ends"));
    };
    let ends = [start, end];
    let points = ends.map(|v| model.geom().position(model.vertex(v).point));
    let Some(&[(fa, sa), (fb, sb)]) = uses.get(&edge).map(Vec::as_slice) else {
        return Err(bad("a spline fillet needs an edge between two faces"));
    };
    let mut first = side(model, fa, sa, [&points[0], &points[1]], dec)?;
    let mut second = side(model, fb, sb, [&points[0], &points[1]], dec)?;
    // The edge must be convex. Then order the supports so the carrier's own
    // normal, along `tau` and from the first rail to the second, points to
    // the ball; the fitted face proves it below.
    let [a, b] = [
        first.directions(model, dec)?,
        second.directions(model, dec)?,
    ];
    if dot(a[1], b[2]) >= 0.0 || dot(b[1], a[2]) >= 0.0 {
        return Err(Refusal::new(
            UNSUPPORTED_SPLINE,
            "only a convex edge between spline faces is rounded",
        ));
    }
    let across: [f64; 3] = std::array::from_fn(|i| b[1][i] - a[1][i]);
    let ball: [f64; 3] = std::array::from_fn(|i| a[1][i] + b[1][i]);
    if dot(cross(a[0], across), ball) < 0.0 {
        std::mem::swap(&mut first, &mut second);
    }
    // Each patch extended by a quarter of the edge's width each way along
    // it, and the face fitted an eighth of the width past each end.
    let extended = |side: &Side, dec: &mut Decider| {
        let mut domain = side.patch.domain().clone();
        let along = side.along();
        let quarter = (&domain[along][1] - &domain[along][0]).mul_pow2(-2);
        domain[along] = [&domain[along][0] - &quarter, &domain[along][1] + &quarter];
        side.patch.extended(domain, dec)
    };
    let supports = [extended(&first, dec)?, extended(&second, dec)?];
    let face = &first.patch.domain()[first.along()];
    let eighth = (&face[1] - &face[0]).mul_pow2(-3);
    let range = [&face[0] - &eighth, &face[1] + &eighth];
    let mut borders = Vec::with_capacity(4);
    for (support, side) in [&first, &second].into_iter().enumerate() {
        for value in &side.ends {
            borders.push(RailBorder {
                support,
                parameter: side.along(),
                value: value.clone(),
            });
        }
    }
    let definition = SplineFilletDefinition::new(
        supports,
        radius,
        [first.ball(model), second.ball(model)],
        first.along(),
        range,
    )?
    .with_borders(borders)?;
    let fillet = SplineFillet::fit(definition, dec)?;
    if !fillet.faces_ball(dec)? {
        return Err(Refusal::new(
            UNSUPPORTED_SPLINE,
            "the fitted face's own normal points away from the ball",
        ));
    }
    // The nodes where each rail meets its face's border rows, at the edge's
    // start and end. Between them each rail must lie on its face.
    let tau: [[Rational; 2]; 2] = {
        let b = fillet.borders();
        [[b[0].clone(), b[1].clone()], [b[2].clone(), b[3].clone()]]
    };
    let forward = tau[0][0] < tau[0][1];
    for (rail, side) in [&first, &second].into_iter().enumerate() {
        if (tau[rail][0] < tau[rail][1]) != forward {
            return Err(bad("the rails cross the face's ends in opposite orders"));
        }
        let [lo, hi] = if forward {
            [&tau[rail][0], &tau[rail][1]]
        } else {
            [&tau[rail][1], &tau[rail][0]]
        };
        if !fillet.rail_within(rail, [lo, hi], side.patch.domain(), dec)? {
            return Err(Refusal::new(
                SPLINE_FILLET_RAIL_DOMAIN,
                "a rail between its face's ends is not proved on the face",
            ));
        }
    }
    // Each end's plane must cut the face from one rail's border node to the
    // other's, and hold the side faces' border edges there.
    let mut caps = [None; 2];
    let mut borders_at = [[edge; 2]; 2];
    let mut rail_points: [[Option<Point3>; 2]; 2] = Default::default();
    for (end, vertex) in ends.iter().enumerate() {
        let cap = cap_at(model, &uses, *vertex, [first.face, second.face])
            .ok_or_else(|| bad("a spline fillet's edge must end at one more face"))?;
        let (normal, offset) = plane_of(model, cap)
            .ok_or_else(|| bad("a spline fillet ends only at rational planes"))?;
        let cut = fillet
            .plane_cut(normal, offset, dec)?
            .ok_or_else(|| bad("an end plane does not cut the fillet"))?;
        if cut.ends != [tau[0][end].clone(), tau[1][end].clone()] {
            return Err(bad(
                "an end plane cuts the fillet away from the rails' borders",
            ));
        }
        for (rail, side) in [&first, &second].into_iter().enumerate() {
            let border = neighbour(model, side.face, edge, *vertex)
                .ok_or_else(|| bad("a side face must meet the end plane in one edge"))?;
            let point = fillet
                .rail_point(rail, &tau[rail][end])
                .expect("a border is a node");
            if !inside_line(model, border, &point, dec)? {
                return Err(bad(
                    "a rail must end inside its face's straight border edge",
                ));
            }
            borders_at[end][rail] = border;
            rail_points[rail][end] = Some(point);
        }
        caps[end] = Some(cap);
    }
    let caps = caps.map(|cap| cap.expect("each end has a cap"));
    let rail_vertices =
        rail_points.map(|row| row.map(|point| vertex(model, point.expect("each rail point"))));

    let mut edit = Replacement::default();
    // The rails, each running from the edge's start to its end as the input
    // edge did, so the side faces keep their coedge senses.
    let mut rails = [edge; 2];
    for (rail, side) in [&first, &second].into_iter().enumerate() {
        let [lo, hi] = if forward {
            [&tau[rail][0], &tau[rail][1]]
        } else {
            [&tau[rail][1], &tau[rail][0]]
        };
        let curve: BSplineCurve = fillet.rail_curve(rail, lo, hi, dec)?;
        let curve = model.geom_mut().add_curve(Curve::BSpline(curve));
        rails[rail] = model.add_edge(Edge {
            curve,
            sense: if forward {
                Sense::Same
            } else {
                Sense::Reversed
            },
            start: Some(rail_vertices[rail][0]),
            end: Some(rail_vertices[rail][1]),
        });
        edit.edges.insert((side.face, edge), rails[rail]);
    }
    // Each plane's cut, from the first rail to the second, and the side
    // faces' border edges shortened to end at the rails.
    let surface = model
        .geom_mut()
        .add_surface(Surface::BSpline(SplineSupport::SplineFillet(fillet)));
    let mut cuts = [edge; 2];
    for (end, vertex) in ends.iter().enumerate() {
        let plane: SurfaceId = model.face(caps[end]).surface;
        let at = model
            .geom()
            .position(model.vertex(rail_vertices[0][end]).point);
        let witness = model.geom_mut().add_point(PointDef::Meet(Meet {
            of: MeetOf::Surfaces(vec![plane, surface]),
            position: Position::Closed(at),
        }));
        let curve = model
            .geom_mut()
            .add_curve(Curve::Intersection(IntersectionCurve {
                surfaces: [plane, surface],
                witness,
            }));
        cuts[end] = model.add_edge(Edge {
            curve,
            sense: Sense::Same,
            start: Some(rail_vertices[0][end]),
            end: Some(rail_vertices[1][end]),
        });
        edit.gaps.insert((caps[end], *vertex), Gap::Edge(cuts[end]));
        for (rail, side) in [&first, &second].into_iter().enumerate() {
            let old = borders_at[end][rail];
            let mut shortened = model.edge(old).clone();
            if shortened.start == Some(*vertex) {
                shortened.start = Some(rail_vertices[rail][end]);
            } else {
                shortened.end = Some(rail_vertices[rail][end]);
            }
            let shortened = model.add_edge(shortened);
            edit.edges.insert((side.face, old), shortened);
            edit.edges.insert((caps[end], old), shortened);
        }
    }
    // The blend's loop runs each rail against its side face's use, and each
    // cut between them.
    let first_use = (rails[0], reverse(first.coedge));
    let reaches = usize::from(first_use.1 == Sense::Same);
    let path = vec![
        first_use,
        (cuts[reaches], Sense::Same),
        (rails[1], reverse(second.coedge)),
        (cuts[1 - reaches], Sense::Reversed),
    ];
    for (i, (e, s)) in path.iter().enumerate() {
        let (next, t) = path[(i + 1) % path.len()];
        if edge_blend::ends(model, *e, *s).1 != edge_blend::ends(model, next, t).0 {
            return Err(bad("the fillet's boundary does not close"));
        }
    }
    let name = model
        .edge_name(body, edge)
        .ok_or_else(|| bad("a selected edge needs an unambiguous origin-pair name"))?;
    let origin = Origin {
        step: step.clone(),
        role: Role::Blend { edge: name },
    };
    // The surface's normal points to the ball, so the face's points out of
    // the material.
    let face = face_on(
        model,
        surface,
        Sense::Reversed,
        vec![path],
        None,
        vec![origin],
    );
    edit.faces.push((first.face, face));
    let output = edge_blend::apply(model, body, &edit)?;
    model.record_blend_base(output, body);
    Ok((output, Vec::new()))
}
