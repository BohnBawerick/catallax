//! Complete oriented domains of exact affine-radius contact-section bands.

use catallax_numbers::{Decider, GEOMETRIC_GAP_UNPROVEN, Outcome, Real, Refusal, Sign, StrictSign};
use catallax_shapes::{BoundaryRow, Curve, LinearFillet, Point3, SplineSupport, Surface};

use crate::{BodyId, Edge, EdgeId, FaceId, Loop, Model, Sense, VertexId};

fn bad(detail: &str) -> Refusal {
    Refusal::new(GEOMETRIC_GAP_UNPROVEN, detail)
}
fn sign(x: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    if let Some(q) = x.rational_ref() {
        return Ok(q.sign());
    }
    match dec.sign(x) {
        Outcome::Proven(StrictSign::Positive) => Ok(Sign::Positive),
        Outcome::Proven(StrictSign::Negative) => Ok(Sign::Negative),
        Outcome::ProvenCoincident(_) => Ok(Sign::Zero),
        Outcome::Refused(r) => Err(r),
        Outcome::NotExact(r) => Err(Refusal::from(r)),
    }
}
fn same(a: &Point3, b: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(sign(&a.minus(b).norm_squared(), dec)? == Sign::Zero)
}
fn ends(
    model: &Model,
    edge: &Edge,
    wanted: [&Point3; 2],
    dec: &mut Decider,
) -> Result<Option<([VertexId; 2], bool)>, Refusal> {
    let (Some(a), Some(b)) = (edge.start, edge.end) else {
        return Ok(None);
    };
    let p = model.geom().position(model.vertex(a).point);
    let q = model.geom().position(model.vertex(b).point);
    if same(&p, wanted[0], dec)? && same(&q, wanted[1], dec)? {
        return Ok(Some(([a, b], true)));
    }
    if same(&p, wanted[1], dec)? && same(&q, wanted[0], dec)? {
        return Ok(Some(([b, a], false)));
    }
    Ok(None)
}

type Row = (EdgeId, Sense, bool);

/// A whole exact variable band, including the quotient at a zero endpoint.
/// This checks physical edges and orientation, not only a carrier's hull.
#[derive(Clone, Copy, Debug)]
pub struct LinearFilletBand<'a> {
    face: FaceId,
    support: &'a LinearFillet,
    rails: [Row; 2],
    caps: [Option<Row>; 2],
}
impl<'a> LinearFilletBand<'a> {
    /// Checks that the face is exactly the carrier's complete contact band.
    /// Additional trims, missing boundaries, wrong senses and duplicated rows
    /// remain unproved. An exact zero end must share one physical apex vertex.
    pub fn of_face(
        model: &'a Model,
        face: FaceId,
        dec: &mut Decider,
    ) -> Result<Option<Self>, Refusal> {
        let f = model.face(face);
        let Surface::BSpline(SplineSupport::LinearFillet(support)) =
            model.geom().surface(f.surface)
        else {
            return Ok(None);
        };
        let [lp] = f.loops.as_slice() else {
            return Err(bad("a complete variable band has one boundary loop"));
        };
        let Loop::Coedges(cs) = model.loop_(*lp) else {
            return Err(bad("a variable band needs its physical boundary edges"));
        };
        let count = 4 - usize::from(support.zero_end(0)) - usize::from(support.zero_end(1));
        if cs.len() != count {
            return Err(bad(
                "variable band boundary is incomplete or additionally trimmed",
            ));
        }
        let mut rails: [Option<Row>; 2] = [None, None];
        let mut caps: [Option<Row>; 2] = [None, None];
        let mut vertices = [[None; 2]; 2];
        for c in cs {
            dec.check_deadline()?;
            let c = model.coedge(*c);
            let e = model.edge(c.edge);
            match model.geom().curve(e.curve) {
                Curve::Line(line) => {
                    let mut found = None;
                    for side in 0..2 {
                        let wanted = [support.foot(side, 0), support.foot(side, 1)];
                        if let Some((ends, forward)) = ends(model, e, wanted, dec)? {
                            if found.replace(side).is_some() {
                                return Err(bad("variable band rails are not distinct"));
                            }
                            let direction = wanted[1].minus(wanted[0]);
                            if sign(&line.dir.norm_squared(), dec)? != Sign::Positive
                                || sign(&line.dir.cross(&direction).norm_squared(), dec)?
                                    != Sign::Zero
                                || sign(
                                    &line
                                        .origin
                                        .minus(wanted[0])
                                        .cross(&direction)
                                        .norm_squared(),
                                    dec,
                                )? != Sign::Zero
                            {
                                return Err(bad(
                                    "variable contact rail does not lie on its exact boundary line",
                                ));
                            }
                            let positive = sign(&line.dir.dot(&direction), dec)? == Sign::Positive;
                            if (e.sense == Sense::Same) != (positive == forward) {
                                return Err(bad("variable rail has an inconsistent carrier sense"));
                            }
                            if rails[side].replace((c.edge, c.sense, forward)).is_some() {
                                return Err(bad("variable band repeats a contact rail"));
                            }
                            vertices[side] = ends.map(Some);
                        }
                    }
                    if found.is_none() {
                        return Err(bad("variable band contains an unproved line trim"));
                    }
                }
                Curve::Circle(circle) => {
                    let mut found = None;
                    for (end, cap) in caps.iter_mut().enumerate() {
                        let Some(section) = support.section(end) else {
                            continue;
                        };
                        let wanted = [support.foot(0, end), support.foot(1, end)];
                        if let Some((_, forward)) = ends(model, e, wanted, dec)? {
                            if found.replace(end).is_some() {
                                return Err(bad("variable end sections are not distinct"));
                            }
                            if !same(&circle.center, &section.center, dec)?
                                || sign(&(&circle.radius_sq - &section.radius_sq), dec)?
                                    != Sign::Zero
                                || sign(&circle.normal.norm_squared(), dec)? != Sign::Positive
                                || sign(&circle.normal.cross(&section.normal).norm_squared(), dec)?
                                    != Sign::Zero
                            {
                                return Err(bad(
                                    "variable band end is not its exact contact-section circle",
                                ));
                            }
                            let positive =
                                sign(&circle.normal.dot(&section.normal), dec)? == Sign::Positive;
                            if (e.sense == Sense::Same) != (positive == forward) {
                                return Err(bad(
                                    "variable end selects the long arc or reverses its carrier",
                                ));
                            }
                            if cap.replace((c.edge, c.sense, forward)).is_some() {
                                return Err(bad("variable band repeats an end section"));
                            }
                        }
                    }
                    if found.is_none() {
                        return Err(bad("variable band contains an unproved circular trim"));
                    }
                }
                _ => {
                    return Err(bad(
                        "variable band has a boundary outside its proved exact rows",
                    ));
                }
            }
        }
        let [Some(a), Some(b)] = rails else {
            return Err(bad("variable band is missing a contact rail"));
        };
        let rails = [a, b];
        for end in 0..2 {
            if support.zero_end(end) {
                if caps[end].is_some() || vertices[0][end] != vertices[1][end] {
                    return Err(bad(
                        "zero-radius boundary does not collapse to one sewn apex",
                    ));
                }
            } else {
                let Some((edge, _, forward)) = caps[end] else {
                    return Err(bad("positive-radius endpoint is missing its exact section"));
                };
                let e = model.edge(edge);
                let ends = if forward {
                    [e.start, e.end]
                } else {
                    [e.end, e.start]
                };
                if ends != [vertices[0][end], vertices[1][end]] {
                    return Err(bad("end section is not sewn to both contact rails"));
                }
            }
        }
        let ccw = f.sense == Sense::Same;
        for (row, wanted) in [(rails[0], false), (rails[1], true)] {
            if ((row.1 == Sense::Same) == row.2) != (wanted == ccw) {
                return Err(bad(
                    "variable band rail loop disagrees with its chart orientation",
                ));
            }
        }
        for (end, row) in caps.iter().enumerate() {
            if let Some(row) = row
                && ((row.1 == Sense::Same) == row.2) != ((end == 0) == ccw)
            {
                return Err(bad(
                    "variable band end loop disagrees with its chart orientation",
                ));
            }
        }
        Ok(Some(Self {
            face,
            support,
            rails,
            caps,
        }))
    }
    /// Checks every variable band and the complete structural boundary of a body.
    /// `false` means the body contains no bands from this exact family.
    pub fn in_body(model: &'a Model, body: BodyId, dec: &mut Decider) -> Result<bool, Refusal> {
        let mut found = false;
        let mut other_spline = false;
        for face in model.faces_of_body(body) {
            match model.geom().surface(model.face(face).surface) {
                Surface::BSpline(SplineSupport::LinearFillet(_)) => {
                    Self::of_face(model, face, dec)?;
                    found = true;
                }
                Surface::BSpline(_) => other_spline = true,
                _ => (),
            }
        }
        if found && (!model.check(body).is_empty() || other_spline) {
            return Err(bad(
                "exact variable-band body has an unproved complete boundary",
            ));
        }
        Ok(found)
    }
    /// The logical face whose entire parameter domain was checked.
    pub fn face(&self) -> FaceId {
        self.face
    }
    /// Retained exact geometric definition and rational carrier.
    pub fn support(&self) -> &'a LinearFillet {
        self.support
    }
    /// Shared physical boundary rows, in station-major display order. The sense
    /// is the face's coedge sense; `forward` starts at that row's first station.
    pub fn rows(&self) -> impl Iterator<Item = (EdgeId, BoundaryRow, Sense, bool)> + '_ {
        let rails = self
            .rails
            .iter()
            .copied()
            .zip([BoundaryRow::Rail, BoundaryRow::Mouth])
            .map(|((edge, sense, forward), row)| (edge, row, sense, forward));
        let caps = self
            .caps
            .iter()
            .copied()
            .zip([BoundaryRow::Negative, BoundaryRow::Positive])
            .filter_map(|(value, row)| {
                value.map(|(edge, sense, forward)| (edge, row, sense, forward))
            });
        rails.chain(caps)
    }
}
