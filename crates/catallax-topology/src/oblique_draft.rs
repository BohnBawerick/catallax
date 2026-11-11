//! Complete finite boundary checks for obliquely drafted cylindrical cells.

use crate::{BodyId, BodyKind, EdgeId, FaceId, Loop, Model, Sense};
use catallax_numbers::{Decider, GEOMETRIC_GAP_UNPROVEN, Outcome, Real, Refusal, Sign, StrictSign};
use catallax_shapes::{
    BSplineCurve, Curve, ObliqueDraft, Plane, Point3, SplineBoundary, SplineSupport, Surface,
};

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
fn sense(s: Sense) -> i64 {
    if s == Sense::Same { 1 } else { -1 }
}
fn ring(model: &Model, face: FaceId, index: usize) -> Result<(EdgeId, Sense), Refusal> {
    let lp = model
        .face(face)
        .loops
        .get(index)
        .ok_or_else(|| bad("missing ruled-cell trim"))?;
    let Loop::Coedges(cs) = model.loop_(*lp) else {
        return Err(bad("ruled trim is not a physical ring"));
    };
    let [c] = cs.as_slice() else {
        return Err(bad("ruled trim contains additional edges"));
    };
    let c = model.coedge(*c);
    let e = model.edge(c.edge);
    if e.start.is_some() || e.end.is_some() {
        return Err(bad("closed ruled ring has seam vertices"));
    }
    Ok((c.edge, c.sense))
}

/// A complete represented solid with one embedded ruled annulus and two disks.
/// The checked carrier proves star-shaped, convex horizontal sections; it does
/// not claim the entire solid is convex or replace a missing trim by its hull.
#[derive(Clone, Copy, Debug)]
pub struct DraftedCylinder<'a> {
    model: &'a Model,
    support: &'a ObliqueDraft,
    faces: [FaceId; 3],
    rings: [EdgeId; 2],
}
impl<'a> DraftedCylinder<'a> {
    /// Checks all shells, faces, exact shared rows and oriented cap domains.
    pub fn check(model: &'a Model, body: BodyId, dec: &mut Decider) -> Result<Self, Refusal> {
        dec.check_deadline()?;
        let entity = model.body(body);
        let [lump] = entity.lumps.as_slice() else {
            return Err(bad("ruled cell needs one material lump"));
        };
        let [shell] = model.lump(*lump).shells.as_slice() else {
            return Err(bad("ruled cell needs one complete shell"));
        };
        let fs = &model.shell(*shell).faces;
        if entity.kind != BodyKind::Solid
            || fs.len() != 3
            || !model.check(body).is_empty()
            || fs[0] == fs[1]
            || fs[0] == fs[2]
            || fs[1] == fs[2]
        {
            return Err(bad("ruled cell is not one valid three-face solid"));
        }
        let mut side = None;
        for &f in fs {
            match model.geom().surface(model.face(f).surface) {
                Surface::BSpline(SplineSupport::ObliqueDraft(s)) if side.is_none() => {
                    side = Some((f, s))
                }
                Surface::Plane(_) => (),
                _ => return Err(bad("ruled cell contains an extra or unsupported carrier")),
            }
        }
        let (side, support) = side.ok_or_else(|| bad("missing checked oblique draft support"))?;
        if model.face(side).sense != Sense::Same || model.face(side).loops.len() != 2 {
            return Err(bad("ruled side has incorrect sense or incomplete trims"));
        }
        let mut caps = [None, None];
        for &f in fs {
            if f == side {
                continue;
            }
            let face = model.face(f);
            let Surface::Plane(p) = model.geom().surface(face.surface) else {
                return Err(bad("ruled cap is not planar"));
            };
            Plane::new(p.origin.clone(), p.normal.clone(), p.ref_dir.clone(), dec)?;
            if sign(&p.normal.cross(support.axis()).norm_squared(), dec)? != Sign::Zero {
                return Err(bad("ruled cap is not perpendicular to pull"));
            }
            let h = p.origin.minus(support.origin()).dot(support.axis());
            let slot = if sign(&(&h - &support.heights()[0]), dec)? == Sign::Zero {
                0
            } else if sign(&(&h - &support.heights()[1]), dec)? == Sign::Zero {
                1
            } else {
                return Err(bad("ruled cap does not meet its exact boundary row"));
            };
            let normal = match sign(&p.normal.dot(support.axis()), dec)? {
                Sign::Positive => 1,
                Sign::Negative => -1,
                Sign::Zero => return Err(bad("zero cap normal")),
            };
            if caps[slot].replace(f).is_some()
                || face.loops.len() != 1
                || normal * sense(face.sense) != [-1, 1][slot]
            {
                return Err(bad("ruled cap is duplicated, incomplete or reversed"));
            }
        }
        let [Some(lower), Some(upper)] = caps else {
            return Err(bad("ruled cell is missing a cap"));
        };
        let faces = [lower, upper, side];
        let mut rings = Vec::with_capacity(2);
        for i in 0..2 {
            let (e, s) = ring(model, faces[i], 0)?;
            let edge = model.edge(e);
            let Curve::BSpline(curve) = model.geom().curve(edge.curve) else {
                return Err(bad("ruled cap ring is not the exact spline row"));
            };
            let expected = BSplineCurve::boundary(
                support.carrier(),
                [SplineBoundary::VMin, SplineBoundary::VMax][i],
            );
            if !curve.same_definition(&expected)
                || curve.is_closed() != Some(true)
                || sense(edge.sense) * sense(s) != [-1, 1][i]
            {
                return Err(bad(
                    "ruled cap ring differs from its complete oriented boundary row",
                ));
            }
            rings.push(e);
        }
        if rings[0] == rings[1] {
            return Err(bad("ruled cell uses one edge for different caps"));
        }
        let mut used = [false, false];
        for i in 0..2 {
            let (e, s) = ring(model, side, i)?;
            let j = rings
                .iter()
                .position(|x| *x == e)
                .ok_or_else(|| bad("ruled side is detached from its cap"))?;
            if std::mem::replace(&mut used[j], true)
                || sense(model.edge(e).sense) * sense(s) != [1, -1][j]
            {
                return Err(bad(
                    "ruled side does not use both rings opposite their cap uses",
                ));
            }
        }
        Ok(Self {
            model,
            support,
            faces,
            rings: [rings[0], rings[1]],
        })
    }
    /// Returns a checked cell only when this body's own faces contain its side.
    pub fn of_body(
        model: &'a Model,
        body: BodyId,
        dec: &mut Decider,
    ) -> Result<Option<Self>, Refusal> {
        if !model.faces_of_body(body).iter().any(|f| {
            matches!(
                model.geom().surface(model.face(*f).surface),
                Surface::BSpline(SplineSupport::ObliqueDraft(_))
            )
        }) {
            return Ok(None);
        }
        Self::check(model, body, dec).map(Some)
    }
    /// Checks the complete solid that owns the queried face, without borrowing
    /// a neighbouring solid's cap or inferring a domain from a bare support.
    pub fn of_face(
        model: &'a Model,
        face: FaceId,
        dec: &mut Decider,
    ) -> Result<Option<Self>, Refusal> {
        for body in model.bodies_with_face(face) {
            if let Some(cell) = Self::of_body(model, body, dec)? {
                return Ok(Some(cell));
            }
        }
        Ok(None)
    }
    /// Negative inside, zero on boundary, positive outside the represented cell.
    pub fn locate(&self, p: &Point3, dec: &mut Decider) -> Result<Sign, Refusal> {
        let band = self.band(p, dec)?;
        if band == Sign::Positive {
            return Ok(band);
        }
        Ok(match self.support.side_at(p, dec)? {
            Sign::Negative => band,
            s => s,
        })
    }
    fn band(&self, p: &Point3, dec: &mut Decider) -> Result<Sign, Refusal> {
        let z = self.support.coordinates(p)[2].clone();
        Ok(
            match (
                sign(&(&z - &self.support.heights()[0]), dec)?,
                sign(&(&z - &self.support.heights()[1]), dec)?,
            ) {
                (Sign::Negative, _) | (_, Sign::Positive) => Sign::Positive,
                (Sign::Zero, _) | (_, Sign::Zero) => Sign::Zero,
                _ => Sign::Negative,
            },
        )
    }
    /// Trim membership for a point already proved on the selected support.
    pub fn locate_on_face(
        &self,
        f: FaceId,
        p: &Point3,
        dec: &mut Decider,
    ) -> Result<Sign, Refusal> {
        if f == self.faces[2] {
            self.band(p, dec)
        } else {
            self.support.side_at(p, dec)
        }
    }
    /// Exact retained ruled support.
    pub fn support(&self) -> &'a ObliqueDraft {
        self.support
    }
    /// Lower cap, upper cap and ruled side, in material order.
    pub fn faces(&self) -> [FaceId; 3] {
        self.faces
    }
    /// Complete lower and upper physical rings.
    pub fn rings(&self) -> [EdgeId; 2] {
        self.rings
    }
    /// Immutable checked model.
    pub fn model(&self) -> &'a Model {
        self.model
    }
}
