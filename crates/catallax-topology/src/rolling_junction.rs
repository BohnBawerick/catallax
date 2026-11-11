//! Boundary rows of a rolling-envelope face.

use catallax_numbers::{Ball, Decider, Outcome, Rational, Refusal, Sign, StrictSign};
use catallax_shapes::{
    BSplineCurve, BoundaryRow, Curve, Point3, RollingEnvelope, RollingPipe, SplineBoundary,
    SplineSupport, Surface,
};

use crate::{EdgeId, FaceId, Loop, Model, Sense};

fn unsupported(detail: &'static str) -> Refusal {
    Refusal::new("topology.rolling-envelope-boundary", detail)
}

fn same_point(a: &Point3, c: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
    match dec.sign(&a.minus(c).norm_squared()) {
        Outcome::ProvenCoincident(_) => Ok(true),
        Outcome::Proven(_) => Ok(false),
        Outcome::NotExact(n) => Err(Refusal::from(n)),
        Outcome::Refused(refusal) => Err(refusal),
    }
}

fn zero(value: &catallax_numbers::Real, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(match dec.sign(value) {
        Outcome::ProvenCoincident(_) => true,
        Outcome::Proven(_) => false,
        Outcome::NotExact(n) => return Err(Refusal::from(n)),
        Outcome::Refused(refusal) => return Err(refusal),
    })
}

// The edge is the carrier's whole row, not a piece or a look-alike. Each row
// that lies on a circle is the minor arc between its corners.
fn is_row(
    model: &Model,
    envelope: &RollingEnvelope,
    row: BoundaryRow,
    edge: EdgeId,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let edge = model.edge(edge);
    let curve = model.geom().curve(edge.curve);
    let Some((centre, radius_sq)) = envelope.boundary_circle(row) else {
        let rail = BSplineCurve::boundary(envelope.carrier(), SplineBoundary::VMin);
        return Ok(matches!(curve, Curve::BSpline(actual) if actual.same_definition(&rail)));
    };
    let Curve::Circle(circle) = curve else {
        return Ok(false);
    };
    if !same_point(&circle.center, &centre, dec)? || !zero(&(&circle.radius_sq - &radius_sq), dec)?
    {
        return Ok(false);
    }
    let (Some(start), Some(end)) = (edge.start, edge.end) else {
        return Ok(false);
    };
    let [first, last] = if edge.sense == Sense::Same {
        [start, end]
    } else {
        [end, start]
    }
    .map(|v| model.geom().position(model.vertex(v).point).minus(&centre));
    match dec.sign(&first.cross(&last).dot(&circle.normal)) {
        Outcome::Proven(catallax_numbers::StrictSign::Positive) => Ok(true),
        Outcome::Proven(_) | Outcome::ProvenCoincident(_) => Ok(false),
        Outcome::NotExact(n) => Err(Refusal::from(n)),
        Outcome::Refused(refusal) => Err(refusal),
    }
}

/// Whether the face retains its envelope's complete finite chart: one outer
/// loop of four coedges that walks the carrier's four whole boundary rows.
/// `false` claims nothing about a point; such a face owes a trimmed-domain
/// answer, as support incidence alone cannot make a point a member of it.
pub fn envelope_retains_chart(
    model: &Model,
    face: FaceId,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let Surface::BSpline(SplineSupport::RollingEnvelope(envelope)) =
        model.geom().surface(model.face(face).surface)
    else {
        return Ok(false);
    };
    let [boundary] = model.face(face).loops.as_slice() else {
        return Ok(false);
    };
    let Loop::Coedges(coedges) = model.loop_(*boundary) else {
        return Ok(false);
    };
    if coedges.len() != 4 {
        return Ok(false);
    }
    let [low, high] =
        [Sign::Negative, Sign::Positive].map(|sign| envelope.definition().corners(sign));
    let rows = [
        (BoundaryRow::Rail, [&low[0], &high[0]]),
        (BoundaryRow::Mouth, [&low[1], &high[1]]),
        (BoundaryRow::Negative, [&low[0], &low[1]]),
        (BoundaryRow::Positive, [&high[0], &high[1]]),
    ];
    let mut seen = [false; 4];
    for (index, coedge) in coedges.iter().enumerate() {
        let use_ = model.coedge(*coedge);
        let edge = model.edge(use_.edge);
        let (Some(start), Some(end)) = (edge.start, edge.end) else {
            return Ok(false);
        };
        let next = model.coedge(coedges[(index + 1) % 4]);
        let next_edge = model.edge(next.edge);
        let leaving = if next.sense == Sense::Same {
            next_edge.start
        } else {
            next_edge.end
        };
        if leaving
            != Some(if use_.sense == Sense::Same {
                end
            } else {
                start
            })
        {
            return Ok(false);
        }
        let points = [start, end].map(|v| model.geom().position(model.vertex(v).point));
        let mut row = None;
        for (at, (candidate, corners)) in rows.iter().enumerate() {
            if (same_point(&points[0], corners[0], dec)?
                && same_point(&points[1], corners[1], dec)?)
                || (same_point(&points[0], corners[1], dec)?
                    && same_point(&points[1], corners[0], dec)?)
            {
                row = Some((at, *candidate));
                break;
            }
        }
        let Some((at, row)) = row else {
            return Ok(false);
        };
        if std::mem::replace(&mut seen[at], true) || !is_row(model, envelope, row, use_.edge, dec)?
        {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Each edge of an envelope face with its chart row and its use by the face.
/// Meridians end on the exact corners of one transition. The rail is the only
/// spline edge and the mouth is the remaining arc.
pub fn envelope_rows(
    model: &Model,
    face: FaceId,
    dec: &mut Decider,
) -> Result<Vec<(EdgeId, BoundaryRow, Sense)>, Refusal> {
    let Surface::BSpline(SplineSupport::RollingEnvelope(envelope)) =
        model.geom().surface(model.face(face).surface)
    else {
        return Err(unsupported("the face is not a rolling envelope"));
    };
    let negative = envelope.definition().corners(Sign::Negative);
    let positive = envelope.definition().corners(Sign::Positive);
    let [boundary] = model.face(face).loops.as_slice() else {
        return Err(unsupported("an envelope face has one outer boundary"));
    };
    let Loop::Coedges(coedges) = model.loop_(*boundary) else {
        return Err(unsupported("an envelope face has one outer boundary"));
    };
    let mut rows = Vec::with_capacity(4);
    for coedge in coedges {
        let use_ = model.coedge(*coedge);
        let edge = model.edge(use_.edge);
        let (Some(start), Some(end)) = (edge.start, edge.end) else {
            return Err(unsupported("an envelope edge has two ends"));
        };
        let points = [start, end].map(|v| model.geom().position(model.vertex(v).point));
        let mut on = |corners: &[Point3; 2]| -> Result<bool, Refusal> {
            Ok((same_point(&points[0], &corners[0], dec)?
                && same_point(&points[1], &corners[1], dec)?)
                || (same_point(&points[0], &corners[1], dec)?
                    && same_point(&points[1], &corners[0], dec)?))
        };
        let row = if matches!(model.geom().curve(edge.curve), Curve::BSpline(_)) {
            BoundaryRow::Rail
        } else if on(&negative)? {
            BoundaryRow::Negative
        } else if on(&positive)? {
            BoundaryRow::Positive
        } else {
            BoundaryRow::Mouth
        };
        rows.push((use_.edge, row, use_.sense));
    }
    Ok(rows)
}

// A closed ring of one edge on the given rail: on its exact ellipse, or a
// cone's wall row itself.
fn pipe_ring(
    model: &Model,
    rails: &[(Curve, bool); 2],
    coedge: crate::CoedgeId,
    dec: &mut Decider,
) -> Result<Option<(EdgeId, bool, Sense)>, Refusal> {
    let use_ = model.coedge(coedge);
    let edge = model.edge(use_.edge);
    if edge.start != edge.end {
        return Ok(None);
    }
    let actual = match model.geom().curve(edge.curve) {
        Curve::Ellipse(actual) => actual,
        Curve::BSpline(actual) => {
            for (rail, is_plane) in rails {
                if let Curve::BSpline(rail) = rail
                    && actual.same_definition(rail)
                {
                    return Ok(Some((use_.edge, *is_plane, use_.sense)));
                }
            }
            return Ok(None);
        }
        _ => return Ok(None),
    };
    for (rail, is_plane) in rails {
        let Curve::Ellipse(rail) = rail else {
            continue;
        };
        if same_point(&actual.center, &rail.center, dec)?
            && zero(&(&actual.major_sq - &rail.major_sq), dec)?
            && zero(&(&actual.minor_sq - &rail.minor_sq), dec)?
            && zero(&actual.normal.cross(&rail.normal).norm_squared(), dec)?
        {
            return Ok(Some((use_.edge, *is_plane, use_.sense)));
        }
    }
    Ok(None)
}

/// A pipe face's two rings, plane rail first, when the face keeps the whole
/// band: two closed loops of one edge each, on the pipe's two exact rails. Each
/// sense is the loop's use of its edge.
pub fn pipe_rows(
    model: &Model,
    face: FaceId,
    dec: &mut Decider,
) -> Result<Option<[(EdgeId, Sense); 2]>, Refusal> {
    let Surface::BSpline(SplineSupport::RollingPipe(pipe)) =
        model.geom().surface(model.face(face).surface)
    else {
        return Ok(None);
    };
    let rails = rails(pipe, dec)?;
    let [first, second] = model.face(face).loops.as_slice() else {
        return Ok(None);
    };
    let mut found = [None; 2];
    for lp in [first, second] {
        let Loop::Coedges(coedges) = model.loop_(*lp) else {
            return Ok(None);
        };
        let [coedge] = coedges.as_slice() else {
            return Ok(None);
        };
        let Some((edge, is_plane, sense)) = pipe_ring(model, &rails, *coedge, dec)? else {
            return Ok(None);
        };
        found[usize::from(!is_plane)] = Some((edge, sense));
    }
    Ok(match found {
        [Some(plane), Some(wall)] => Some([plane, wall]),
        _ => None,
    })
}

/// The complete band's two rails, with each face's coedge sense: a saddle's
/// wide rail first, and a sphere band's sphere rail first. A regular ring has
/// two one-edge loops. A crossing band instead has one two-edge loop whose
/// rails share both pole vertices.
pub fn saddle_rows(model: &Model, face: FaceId) -> Option<[(EdgeId, Sense); 2]> {
    let (rails, crossing) = match model.geom().surface(model.face(face).surface) {
        Surface::BSpline(SplineSupport::RollingSaddle(saddle)) => (
            saddle.rails(),
            matches!(
                saddle.definition().law().ring,
                catallax_shapes::SaddleRing::Crossing { .. }
            ),
        ),
        Surface::BSpline(SplineSupport::RollingSphere(band)) => (band.rails(), false),
        _ => return None,
    };
    let coedges = if crossing {
        let [lp] = model.face(face).loops.as_slice() else {
            return None;
        };
        let Loop::Coedges(coedges) = model.loop_(*lp) else {
            return None;
        };
        let [first, second] = coedges.as_slice() else {
            return None;
        };
        let a = model.edge(model.coedge(*first).edge);
        let b = model.edge(model.coedge(*second).edge);
        if a.start.is_none()
            || a.end.is_none()
            || a.start == a.end
            || !((a.start == b.start && a.end == b.end) || (a.start == b.end && a.end == b.start))
        {
            return None;
        }
        [*first, *second]
    } else {
        let [first, second] = model.face(face).loops.as_slice() else {
            return None;
        };
        let single = |lp| {
            let Loop::Coedges(coedges) = model.loop_(lp) else {
                return None;
            };
            let [coedge] = coedges.as_slice() else {
                return None;
            };
            Some(*coedge)
        };
        [single(*first)?, single(*second)?]
    };
    let mut found = [None; 2];
    for coedge in coedges {
        let use_ = model.coedge(coedge);
        let edge = model.edge(use_.edge);
        if !crossing && edge.start != edge.end {
            return None;
        }
        let Curve::BSpline(actual) = model.geom().curve(edge.curve) else {
            return None;
        };
        let row = rails.iter().position(|rail| actual.same_definition(rail))?;
        found[row] = Some((use_.edge, use_.sense));
    }
    Some([found[0]?, found[1]?])
}

/// The cap of a sphere band's sphere rail that holds a point of the sphere
/// centred at `centre`: `Some(true)` on the side the rail turns
/// counterclockwise about, its band's turn axis, and `Some(false)` on the
/// other. A mirror reverses that axis against the cylinder's, so reading the
/// turn rather than the cylinder axis keeps a mirrored copy's caps. The rail winds once round the axis, so it cuts
/// the sphere into a cap on each side, one holding each pole. The meridian arc
/// from the point to a pole meets the rail only in spans whose control hulls
/// reach the arc's half-plane, so a point above or below every control of
/// those spans along the axis lies in that pole's cap. Positive weights keep
/// each span in its hull. `None` when `edge` is not a whole sphere rail or the
/// hulls do not separate the point from the rail.
pub fn band_cap(
    model: &Model,
    edge: EdgeId,
    centre: &Point3,
    point: &Point3,
    dec: &mut Decider,
) -> Result<Option<bool>, Refusal> {
    const BITS: u32 = 96;
    let source = model.edge(edge);
    if source.start != source.end {
        return Ok(None);
    }
    let Curve::BSpline(rail) = model.geom().curve(source.curve) else {
        return Ok(None);
    };
    let band = model
        .coedges_of_edge(edge)
        .into_iter()
        .filter_map(|coedge| model.owner_of_coedge(coedge))
        .find_map(
            |(_, owner)| match model.geom().surface(model.face(owner).surface) {
                Surface::BSpline(SplineSupport::RollingSphere(band))
                    if band.rails()[0].same_definition(rail) =>
                {
                    Some(band)
                }
                _ => None,
            },
        );
    let Some(band) = band else {
        return Ok(None);
    };
    let axis = &band.definition().turn_axis();
    let relative = point.minus(centre);
    let at = dec.enclose_at(&relative.dot(axis), BITS)?;
    // The arc's half-plane holds the axis and the point's offset from it.
    let ahead = relative.minus(&axis.scaled(&relative.dot(axis)));
    let normal = axis.cross(&ahead);
    let positive = |value: &Ball| value.lo() > Rational::zero();
    let negative = |value: &Ball| value.hi() < Rational::zero();
    let degree = rail.degree();
    let (mut above, mut below) = (true, true);
    for (index, _) in rail.spans() {
        let controls = index - degree..=index;
        let mut sides = Vec::with_capacity(degree + 1);
        let mut fronts = Vec::with_capacity(degree + 1);
        let mut heights = Vec::with_capacity(degree + 1);
        for i in controls {
            dec.check_deadline()?;
            let offset = rail.control(i).minus(centre);
            sides.push(dec.enclose_at(&offset.dot(&normal), BITS)?);
            fronts.push(dec.enclose_at(&offset.dot(&ahead), BITS)?);
            heights.push(dec.enclose_at(&offset.dot(axis), BITS)?);
        }
        let clear =
            sides.iter().all(positive) || sides.iter().all(negative) || fronts.iter().all(negative);
        if clear {
            continue;
        }
        for height in &heights {
            above &= at.lo() > height.hi();
            below &= at.hi() < height.lo();
        }
    }
    Ok(if above {
        Some(true)
    } else if below {
        Some(false)
    } else {
        None
    })
}

/// Whether a ring band's pocket adds material. A cylinder or sphere face keeps
/// its side of the material through any isometry, and the ball sits in the
/// material exactly at a convex edge, so a rail's neighbour settles it.
pub fn saddle_added(model: &Model, face: FaceId) -> Result<bool, Refusal> {
    let refused = || {
        Refusal::new(
            "topology.saddle-band",
            "a ring band's rail has no neighbour on its support",
        )
    };
    // The first rail's support: a saddle's wide cylinder, or a sphere band's
    // sphere. Both keep their material inside under `Sense::Same`.
    let (side, sphere) = match model.geom().surface(model.face(face).surface) {
        Surface::BSpline(SplineSupport::RollingSaddle(saddle)) => {
            (saddle.definition().law().sides[0], false)
        }
        Surface::BSpline(SplineSupport::RollingSphere(band)) => {
            (band.definition().sides()[0], true)
        }
        _ => return Err(refused()),
    };
    let [(rail, _), _] = saddle_rows(model, face).ok_or_else(refused)?;
    for coedge in model.coedges_of_edge(rail) {
        let Some((_, neighbour)) = model.owner_of_coedge(coedge) else {
            continue;
        };
        let boundary = model.face(neighbour);
        let support = match model.geom().surface(boundary.surface) {
            Surface::Cylinder(_) => !sphere,
            Surface::Sphere(_) => sphere,
            _ => false,
        };
        if neighbour != face && support {
            let material_inside = boundary.sense == Sense::Same;
            let ball_inside = side < 0;
            return Ok(material_inside != ball_inside);
        }
    }
    Err(refused())
}

/// Whether a band's outward side is along its contact rays, away from the ball
/// centre. The plane that shares its plane rail fixes this.
pub fn pipe_outward(
    model: &Model,
    face: FaceId,
    rail: EdgeId,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let Surface::BSpline(SplineSupport::RollingPipe(pipe)) =
        model.geom().surface(model.face(face).surface)
    else {
        return Err(no_plane());
    };
    let law = pipe.law();
    for coedge in model.coedges_of_edge(rail) {
        let Some((_, neighbour)) = model.owner_of_coedge(coedge) else {
            continue;
        };
        let boundary = model.face(neighbour);
        if neighbour == face {
            continue;
        }
        if let Surface::Plane(plane) = model.geom().surface(boundary.surface) {
            let ray = &law.normal.dot(&plane.normal) * &law.reach;
            let along = match dec.sign(&ray) {
                Outcome::Proven(StrictSign::Positive) => true,
                Outcome::Proven(StrictSign::Negative) => false,
                Outcome::Refused(refusal) => return Err(refusal),
                Outcome::NotExact(n) => return Err(no_plane().keeping(n)),
                Outcome::ProvenCoincident(_) => return Err(no_plane()),
            };
            return Ok(along == (boundary.sense == Sense::Same));
        }
    }
    Err(no_plane())
}

fn no_plane() -> Refusal {
    Refusal::new(
        "topology.pipe-band",
        "a pipe band's plane rail has no plane neighbour",
    )
}

fn rails(pipe: &RollingPipe, dec: &mut Decider) -> Result<[(Curve, bool); 2], Refusal> {
    Ok([
        (Curve::Ellipse(pipe.plane_rail(dec)?), true),
        (pipe.wall_curve(dec)?, false),
    ])
}
