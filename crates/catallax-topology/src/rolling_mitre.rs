//! The whole band a ball sweeps over a mitre, read from a face's edges.

use catallax_numbers::{Decider, Outcome, Rational, Real, Refusal};
use catallax_shapes::{BoundaryRow, Curve, Point3, RollingMitre, SplineSupport, Surface, Vector3};

use crate::{BodyId, EdgeId, FaceId, Loop, Model, Sense};

// An edge, its coedge sense in the face, and whether it starts at its row's
// first chart point.
type Row = (EdgeId, Sense, bool);

/// A face that keeps its [`RollingMitre`] carrier's whole chart: one loop of
/// the end section's arc and both contact rails, meeting at the band's two
/// feet and its apex.
#[derive(Clone, Copy, Debug)]
pub struct RollingMitreBand<'a> {
    /// The band's face.
    pub face: FaceId,
    /// Its exact support.
    pub support: &'a RollingMitre,
    // The contact rails, then the end section.
    rows: [Row; 3],
}

fn zero(value: Real, dec: &mut Decider) -> bool {
    matches!(dec.sign(&value), Outcome::ProvenCoincident(_))
}
// Equal definitions prove equal values without any refinement; only
// different constructions of one value need a proof.
fn equal(a: &Real, b: &Real, dec: &mut Decider) -> bool {
    a.same_definition(b) || zero(a - b, dec)
}
fn identical([a, b]: [[&Real; 3]; 2]) -> bool {
    (0..3).all(|i| a[i].same_definition(b[i]))
}
fn same(a: &Point3, b: &Point3, dec: &mut Decider) -> bool {
    identical([[&a.x, &a.y, &a.z], [&b.x, &b.y, &b.z]]) || zero(a.minus(b).norm_squared(), dec)
}
fn parallel(a: &Vector3, b: &Vector3, dec: &mut Decider) -> bool {
    identical([[&a.x, &a.y, &a.z], [&b.x, &b.y, &b.z]]) || zero(a.cross(b).norm_squared(), dec)
}

impl<'a> RollingMitreBand<'a> {
    /// Checks that the face is exactly the carrier's complete band. Other
    /// trims, a missing boundary or a moved end remain unproved.
    pub fn of_face(
        model: &'a Model,
        face: FaceId,
        dec: &mut Decider,
    ) -> Result<Option<Self>, Refusal> {
        let f = model.face(face);
        let Surface::BSpline(SplineSupport::RollingMitre(support)) =
            model.geom().surface(f.surface)
        else {
            return Ok(None);
        };
        let [lp] = f.loops[..] else {
            return Ok(None);
        };
        let Loop::Coedges(coedges) = model.loop_(lp) else {
            return Ok(None);
        };
        if coedges.len() != 3 {
            return Ok(None);
        }
        let section = support.section(dec)?;
        let rails = [support.rail(0, dec)?, support.rail(1, dec)?];
        let feet =
            [Rational::zero(), Rational::one()].map(|v| support.point(&Rational::zero(), &v));
        let apex = support.apex();
        let mut seen = [false; 3];
        let mut rows: [Option<Row>; 3] = [None; 3];
        for coedge in coedges {
            let use_ = model.coedge(*coedge);
            let edge = model.edge(use_.edge);
            let (Some(start), Some(end)) = (edge.start, edge.end) else {
                return Ok(None);
            };
            let ends = [start, end].map(|v| model.geom().position(model.vertex(v).point));
            let joins = |a: &Point3, b: &Point3, dec: &mut Decider| {
                (same(&ends[0], a, dec) && same(&ends[1], b, dec))
                    || (same(&ends[0], b, dec) && same(&ends[1], a, dec))
            };
            match model.geom().curve(edge.curve) {
                Curve::Circle(circle)
                    if !seen[2]
                        && same(&circle.center, &section.center, dec)
                        && equal(&circle.radius_sq, &section.radius_sq, dec)
                        && parallel(&circle.normal, &section.normal, dec)
                        && joins(&feet[0], &feet[1], dec) =>
                {
                    seen[2] = true;
                    rows[2] = Some((use_.edge, use_.sense, same(&ends[0], &feet[0], dec)));
                }
                Curve::Ellipse(ellipse) => {
                    let mut matched = false;
                    for (i, rail) in rails.iter().enumerate() {
                        if !seen[i]
                            && same(&ellipse.center, &rail.center, dec)
                            && parallel(&ellipse.normal, &rail.normal, dec)
                            && equal(&ellipse.major_sq, &rail.major_sq, dec)
                            && joins(&feet[i], &apex, dec)
                        {
                            seen[i] = true;
                            rows[i] = Some((use_.edge, use_.sense, same(&ends[0], &feet[i], dec)));
                            matched = true;
                            break;
                        }
                    }
                    if !matched {
                        return Ok(None);
                    }
                }
                _ => return Ok(None),
            }
        }
        let [Some(a), Some(b), Some(c)] = rows else {
            return Ok(None);
        };
        Ok((seen == [true; 3]).then_some(Self {
            face,
            support,
            rows: [a, b, c],
        }))
    }
    /// Whether a body holds a band over a mitre and all its spline faces are
    /// whole exact bands: such bands or rolling pipes. A band that is not
    /// its carrier's whole chart, or a body that fails its checks, refuses.
    pub fn in_body(model: &'a Model, body: BodyId, dec: &mut Decider) -> Result<bool, Refusal> {
        let mut found = false;
        let mut other_spline = false;
        for face in model.faces_of_body(body) {
            match model.geom().surface(model.face(face).surface) {
                Surface::BSpline(SplineSupport::RollingMitre(_)) => {
                    if Self::of_face(model, face, dec)?.is_none() {
                        return Err(Refusal::new(
                            catallax_numbers::GEOMETRIC_GAP_UNPROVEN,
                            "a rolling mitre band is not its carrier's whole chart",
                        ));
                    }
                    found = true;
                }
                Surface::BSpline(SplineSupport::RollingPipe(_)) => {}
                Surface::BSpline(_) => other_spline = true,
                _ => {}
            }
        }
        if found && !model.check(body).is_empty() {
            return Err(Refusal::new(
                catallax_numbers::GEOMETRIC_GAP_UNPROVEN,
                "a body with a rolling mitre band has an unproved boundary",
            ));
        }
        Ok(found && !other_spline)
    }
    /// The boundary rows in display order: the first rail along chart `u` at
    /// `v = 0`, the second at `v = 1`, and the end section along `v` at
    /// `u = 0`. The apex column has no edge. The sense is the face's coedge
    /// sense; `forward` starts at that row's first chart point.
    pub fn rows(&self) -> impl Iterator<Item = (EdgeId, BoundaryRow, Sense, bool)> + '_ {
        self.rows
            .iter()
            .zip([BoundaryRow::Rail, BoundaryRow::Mouth, BoundaryRow::Negative])
            .map(|((edge, sense, forward), row)| (*edge, row, *sense, *forward))
    }
}
